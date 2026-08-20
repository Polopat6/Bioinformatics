"""
qc_certificate_manager.py

Builds a single, self-contained "QC Certificate" for a completed (or
errored) Advanced Mode / Monitor Mode pipeline run -- a machine-readable
JSON record plus a standalone HTML overview report -- summarizing every
QC/flag signal already produced across the run's stages into ONE place.

--- Why this exists ---
Today, QC information is real but scattered: FastQC/MultiQC flags,
fastp's poly-tail auto-fix events, per-sample mapping-rate quality tiers
(classify_mapping_rate), and Monitor's own folder-validation warnings
all exist independently, with no single "did this pass, and if not,
what exactly was flagged?" answer. For Auto Mode, that means a solo
researcher has to dig through several separate report files before
trusting their own run. For Monitor Mode, it's a bigger problem: a
counts matrix produced by a monitor is frequently handed off to SOMEONE
ELSE (a collaborator, a core lab, a PI) who relies on it being
QC'd and ready to go, but has no practical way to open FastQC/MultiQC
HTML reports themselves before starting their own downstream analysis.

This module produces exactly one machine-readable file
(qc_certificate.json) and one human-readable file (qc_report.html,
fully self-contained -- inline CSS, no external assets, no server
needed -- opens standalone in any browser) that travel WITH the
project's results package (see project_actions.py), so the QC status
is visible the moment someone opens the delivered data, not something
they have to come back to the portal to check.

--- Scope ---
This module only READS already-produced artifacts (FastQC/MultiQC
output, fastp JSON reports, the samplesheet, the counts matrix, the
saved reference choice) -- it does not run any new QC tool itself, and
has no Streamlit dependency, so it can be called identically from the
non-interactive orchestrator (advanced_mode_orchestrator.py, for both
Auto and Monitor-launched runs) and, if ever wanted, from an
interactive workspace too.

--- Overall status ---
A run's overall_status is "flagged" if ANY per-stage flag exists,
"passed" otherwise. This is deliberately a simple, conservative rollup
(any single flag -> not a clean pass) rather than a weighted/scored
system, since the audience for this (a downstream collaborator, or a
solo researcher deciding whether to trust their own run before DESeq2)
needs an unambiguous yes/no signal first -- the underlying detail
(exactly which stage, which sample, which check) is still fully
preserved in both the JSON and the HTML for anyone who wants to look
closer.

--- Reference provenance ---
The "reference" section of the certificate records which species/
assembly/source was used, and (best-effort, from whatever is currently
on disk) when the reference files were downloaded/built. This section
is intentionally best-effort for now: once the planned private,
permanent per-run reference copy (Monitor: monitor-owned; Auto:
per-project) lands, this section should be extended to read that
copy's own dedicated provenance manifest (exact source URL/assembly,
download timestamp, build parameters) instead of inferring timestamps
from shared-directory file mtimes as it does today.
"""

import json
import os
from datetime import datetime

import pandas as pd

import project_manager as pm
import reference_manager as rm
import quantification_manager as qm
import ingestion_manager as ingest
import fastqc_manager as fastqc
import fastp_manager as fastp


# ---------------------------------------------------------------------------
# Per-stage data collection -- each function reads whatever artifacts that
# stage already produced on disk and returns a plain, JSON-serializable
# dict with the same shape: {"status": "ok"|"flagged"|"unavailable",
# "flags": [str, ...], "detail": {...stage-specific...}}
# ---------------------------------------------------------------------------

def _collect_ingest_summary(project_name):
    samplesheet_path = pm.samplesheet_path(project_name)
    if not os.path.isfile(samplesheet_path):
        return {"status": "unavailable", "flags": [], "detail": {}}

    flags = []
    matched_samples = []
    try:
        samplesheet_df = pd.read_csv(samplesheet_path)
        matched_samples = sorted(samplesheet_df["sample"].astype(str).tolist()) if "sample" in samplesheet_df.columns else []
    except Exception:
        pass

    # Best-effort reconciliation against the original metadata, mirroring
    # ingestion_manager.build_match_table's own logic, so unmatched
    # samples (declared in metadata but with no FASTQ data, or vice
    # versa) are surfaced here too, not just at ingestion time.
    metadata_path = pm.metadata_path(project_name)
    if os.path.isfile(metadata_path):
        try:
            meta_df, read_error = ingest.read_metadata_file(metadata_path)
            if not read_error:
                fastq_names = ingest.list_existing_fastq(pm.fastq_dir(project_name))
                sample_pairs = ingest.validate_sample_pairs(fastq_names)
                match_table = ingest.build_match_table(sample_pairs, meta_df)
                unmatched = match_table[~match_table["Status"].astype(str).str.contains("✅")]
                for _, row in unmatched.iterrows():
                    flags.append(f"Sample '{row.get('Sample', '?')}' did not fully match: {row.get('Status', 'unknown status')}")
        except Exception:
            pass

    return {
        "status": "flagged" if flags else "ok",
        "flags": flags,
        "detail": {"n_matched_samples": len(matched_samples), "matched_samples": matched_samples},
    }


def _collect_pretrim_qc_summary(project_name):
    fastqc_dir = pm.fastqc_dir(project_name)
    if not os.path.isdir(fastqc_dir):
        return {"status": "unavailable", "flags": [], "detail": {}}

    flags = []
    per_file = []
    try:
        summary_df = fastqc.parse_fastqc_summaries(fastqc_dir)
        if summary_df is not None and not summary_df.empty:
            overview_df, details_by_file = fastqc.build_quality_flags(summary_df)
            for _, row in overview_df.iterrows():
                file_name = str(row.get("File", "?"))
                quality = str(row.get("Overall Quality", ""))
                per_file.append({"file": file_name, "overall_quality": quality})
                # Anything not a clean "good"/"pass" result is treated as
                # a flag -- checked defensively via the emoji/keyword
                # this app's own convention uses elsewhere (✅ = clean),
                # rather than assuming one fixed string.
                if "✅" not in quality:
                    explanations = details_by_file.get(file_name) or []
                    detail_text = "; ".join(explanations) if explanations else quality
                    flags.append(f"FastQC flagged '{file_name}': {detail_text}")
    except Exception as e:
        return {"status": "unavailable", "flags": [], "detail": {"error": str(e)}}

    return {
        "status": "flagged" if flags else "ok",
        "flags": flags,
        "detail": {"per_file": per_file},
    }


def _collect_trimming_summary(project_name):
    reports_dir = pm.fastp_reports_dir(project_name)
    if not os.path.isdir(reports_dir):
        return {"status": "unavailable", "flags": [], "detail": {}}

    flags = []
    per_sample = []
    try:
        summary_df = fastp.parse_fastp_reports(reports_dir)
        if summary_df is not None and not summary_df.empty:
            per_sample = summary_df.drop(columns=["_raw_adapter_cutting_json"], errors="ignore").to_dict(orient="records")
    except Exception:
        pass

    poly_tail_events = []
    try:
        poly_tail_issues = fastp.scan_all_samples_for_poly_tail_issues(reports_dir)
        for sample_name, issues in (poly_tail_issues or {}).items():
            for issue in issues:
                poly_tail_events.append({
                    "sample": sample_name, "read": issue.get("read"), "base": issue.get("base"),
                    "tail_pct_after": issue.get("tail_pct_after"),
                })
                # Note: this is informational, not necessarily a failure --
                # the orchestrator auto-fixes these when auto_fix_poly_tails
                # is enabled (the default). Still worth surfacing so a
                # downstream reviewer knows a residual-tail correction
                # happened, rather than only seeing "trimming: ok".
                flags.append(
                    f"Sample '{sample_name}': residual poly-{issue.get('base')} tail detected on {issue.get('read')} "
                    f"({issue.get('tail_pct_after')}% in the tail window) -- auto-corrected if poly-tail auto-fix was enabled."
                )
    except Exception:
        pass

    return {
        "status": "flagged" if flags else "ok",
        "flags": flags,
        "detail": {"per_sample": per_sample, "poly_tail_events": poly_tail_events},
    }


def _collect_reference_summary(project_name):
    species_key, is_custom = pm.get_reference_choice(project_name)
    if is_custom:
        reference_dir = pm.reference_dir(project_name)
        label, source, assembly = "Custom (user-uploaded) reference", "custom_upload", None
    elif species_key:
        reference_dir = pm.shared_reference_dir(species_key)
        entry = rm.REFERENCE_CATALOG.get(species_key, {})
        label = entry.get("label", species_key)
        source = entry.get("source")
        assembly = entry.get("assembly")
    else:
        return {"status": "unavailable", "flags": [], "detail": {}}

    # Best-effort provenance: earliest mtime among the reference's own
    # files as a stand-in for "when this reference was prepared" --
    # see this module's docstring's "Reference provenance" note for why
    # this is a placeholder pending the planned per-run private-copy
    # provenance manifest.
    prepared_at = None
    if os.path.isdir(reference_dir):
        try:
            mtimes = [
                os.path.getmtime(os.path.join(root, f))
                for root, _dirs, files in os.walk(reference_dir) for f in files
            ]
            if mtimes:
                prepared_at = datetime.fromtimestamp(min(mtimes)).isoformat(timespec="seconds")
        except OSError:
            pass

    return {
        "status": "ok",
        "flags": [],
        "detail": {
            "species_key": species_key, "label": label, "source": source,
            "assembly": assembly, "is_custom": is_custom,
            "reference_dir": reference_dir, "prepared_at": prepared_at,
        },
    }


def _collect_quantification_summary(project_name):
    method = pm.get_alignment_method(project_name)
    samplesheet_path = pm.samplesheet_path(project_name)
    if not method or not os.path.isfile(samplesheet_path):
        return {"status": "unavailable", "flags": [], "detail": {}}

    flags = []
    per_sample = []
    try:
        samplesheet_df = pd.read_csv(samplesheet_path)
        trimmed_dir = pm.trimmed_fastq_dir(project_name)
        manifest = qm.build_sample_manifest(samplesheet_df, trimmed_dir)

        if method == "salmon":
            quant_dir = pm.salmon_quant_dir(project_name)
            for entry in manifest:
                sample_out_dir = os.path.join(quant_dir, entry["sample"])
                rate = qm.parse_salmon_mapping_rate(sample_out_dir)
                classification = qm.classify_mapping_rate(rate)
                per_sample.append({"sample": entry["sample"], "mapping_rate_pct": rate, "quality_tier": classification["tier"]})
                if classification["tier"] in ("poor", "caution"):
                    flags.append(f"Sample '{entry['sample']}': {classification['message']}")
        else:
            align_dir = pm.star_align_dir(project_name)
            for entry in manifest:
                sample_out_dir = os.path.join(align_dir, entry["sample"])
                rate = qm.parse_star_mapping_rate(sample_out_dir, entry["sample"])
                classification = qm.classify_mapping_rate(rate)
                per_sample.append({"sample": entry["sample"], "mapping_rate_pct": rate, "quality_tier": classification["tier"]})
                if classification["tier"] in ("poor", "caution"):
                    flags.append(f"Sample '{entry['sample']}': {classification['message']}")
    except Exception as e:
        return {"status": "unavailable", "flags": [], "detail": {"error": str(e)}}

    return {
        "status": "flagged" if flags else "ok",
        "flags": flags,
        "detail": {"alignment_method": method, "per_sample": per_sample},
    }


def _collect_counts_matrix_summary(project_name):
    counts_path = pm.counts_matrix_path(project_name)
    if not os.path.isfile(counts_path):
        return {"status": "unavailable", "flags": [], "detail": {}}

    flags = []
    n_genes, n_samples, matrix_samples = 0, 0, []
    try:
        matrix_df = pd.read_csv(counts_path)
        n_genes = len(matrix_df)
        matrix_samples = [c for c in matrix_df.columns if c != "gene_id"]
        n_samples = len(matrix_samples)
    except Exception as e:
        return {"status": "unavailable", "flags": [], "detail": {"error": str(e)}}

    samplesheet_path = pm.samplesheet_path(project_name)
    if os.path.isfile(samplesheet_path):
        try:
            samplesheet_df = pd.read_csv(samplesheet_path)
            expected_samples = set(samplesheet_df["sample"].astype(str)) if "sample" in samplesheet_df.columns else set()
            missing = sorted(expected_samples - set(matrix_samples))
            for s in missing:
                flags.append(f"Sample '{s}' was matched during ingestion but has no quantification output in the final counts matrix.")
        except Exception:
            pass

    return {
        "status": "flagged" if flags else "ok",
        "flags": flags,
        "detail": {"n_genes": n_genes, "n_samples": n_samples, "samples": matrix_samples},
    }


_STAGE_COLLECTORS = [
    ("ingest", "1. FASTQ Ingestion & Sample Matching", _collect_ingest_summary),
    ("qc", "2. Pre-trim QC (FastQC + MultiQC)", _collect_pretrim_qc_summary),
    ("trimming", "3. Trimming (fastp)", _collect_trimming_summary),
    ("reference", "4. Reference Genome & Index", _collect_reference_summary),
    ("quantification", "5. Quantification (Salmon/STAR)", _collect_quantification_summary),
    ("counts_matrix", "6. Gene Counts Matrix", _collect_counts_matrix_summary),
]


# ---------------------------------------------------------------------------
# Certificate assembly
# ---------------------------------------------------------------------------

def build_qc_certificate(project_name, run_mode="auto", monitor_id=None):
    """
    Assemble the full QC certificate dict for project_name by reading
    every already-produced artifact across all six pipeline stages.

    run_mode: "auto" or "monitor" -- recorded in the certificate so a
        downstream reviewer knows which path produced this project,
        and (for "monitor") which monitor launched it.
    monitor_id: the launching monitor's ID, if run_mode == "monitor".

    Returns a plain, JSON-serializable dict -- safe to call
    json.dumps() on directly, and safe to call repeatedly/idempotently
    (e.g. once after a successful run, once after an errored run -- see
    save_qc_certificate below).
    """
    stages = {}
    all_flags = []
    for stage_key, stage_label, collector in _STAGE_COLLECTORS:
        result = collector(project_name)
        stages[stage_key] = {"label": stage_label, **result}
        for flag in result.get("flags", []):
            all_flags.append(f"[{stage_label}] {flag}")

    overall_status = "flagged" if all_flags else "passed"

    return {
        "project_name": project_name,
        "run_mode": run_mode,
        "monitor_id": monitor_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "overall_status": overall_status,
        "flags": all_flags,
        "stages": stages,
    }


# ---------------------------------------------------------------------------
# HTML report rendering -- fully self-contained (inline CSS, no external
# assets), so it opens standalone in any browser with no server needed,
# suitable for handing off alongside the project package zip.
# ---------------------------------------------------------------------------

_HTML_STYLE = """
<style>
  body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 2rem; background: #f6f7f9; color: #1a1a1a; }
  .container { max-width: 900px; margin: 0 auto; }
  h1 { font-size: 1.6rem; margin-bottom: 0.25rem; }
  .subtitle { color: #555; margin-bottom: 1.5rem; }
  .banner { padding: 1rem 1.25rem; border-radius: 8px; font-size: 1.15rem; font-weight: 600; margin-bottom: 1.5rem; }
  .banner.passed { background: #e6f4ea; color: #1e7a34; border: 1px solid #b7e0c3; }
  .banner.flagged { background: #fdecea; color: #a33; border: 1px solid #f3c2bd; }
  .stage { background: white; border-radius: 8px; border: 1px solid #e2e4e8; margin-bottom: 1rem; overflow: hidden; }
  .stage-header { padding: 0.85rem 1.1rem; font-weight: 600; display: flex; justify-content: space-between; align-items: center; }
  .stage-header.ok { background: #f0faf3; }
  .stage-header.flagged { background: #fff5f4; }
  .stage-header.unavailable { background: #f4f4f5; color: #888; }
  .badge { font-size: 0.85rem; padding: 0.15rem 0.6rem; border-radius: 999px; }
  .badge.ok { background: #d6f0dd; color: #1e7a34; }
  .badge.flagged { background: #fbdad6; color: #a33; }
  .badge.unavailable { background: #e5e5e7; color: #777; }
  .stage-body { padding: 0.9rem 1.1rem; font-size: 0.92rem; }
  .flag-item { padding: 0.4rem 0; border-bottom: 1px solid #f0f0f0; }
  .flag-item:last-child { border-bottom: none; }
  table { width: 100%; border-collapse: collapse; margin-top: 0.5rem; font-size: 0.88rem; }
  th, td { text-align: left; padding: 0.35rem 0.6rem; border-bottom: 1px solid #eee; }
  th { color: #666; font-weight: 600; }
  .meta { color: #777; font-size: 0.85rem; margin-bottom: 1.5rem; }
  .no-flags { color: #1e7a34; }
</style>
"""


def _render_stage_html(stage_key, stage):
    status = stage.get("status", "unavailable")
    flags = stage.get("flags", [])
    label = stage.get("label", stage_key)

    body_html = ""
    if status == "unavailable":
        body_html = "<p>No data available for this stage yet.</p>"
    elif flags:
        items = "".join(f'<div class="flag-item">⚠️ {_escape(f)}</div>' for f in flags)
        body_html = items
    else:
        body_html = '<p class="no-flags">✅ No issues detected.</p>'

    detail = stage.get("detail", {})
    if stage_key == "quantification" and detail.get("per_sample"):
        rows = "".join(
            f"<tr><td>{_escape(str(r.get('sample')))}</td>"
            f"<td>{_escape(str(r.get('mapping_rate_pct')))}%</td>"
            f"<td>{_escape(str(r.get('quality_tier')))}</td></tr>"
            for r in detail["per_sample"]
        )
        body_html += (
            "<table><tr><th>Sample</th><th>Mapping Rate</th><th>Quality</th></tr>"
            + rows + "</table>"
        )
    elif stage_key == "counts_matrix" and detail.get("n_genes"):
        body_html += f"<p>{detail['n_genes']:,} gene(s) × {detail['n_samples']} sample(s) in the final matrix.</p>"
    elif stage_key == "reference" and detail.get("label"):
        prepared_note = f" (prepared {detail['prepared_at']})" if detail.get("prepared_at") else ""
        body_html += f"<p>{_escape(detail['label'])}{prepared_note}</p>"

    return f"""
    <div class="stage">
      <div class="stage-header {status}">
        <span>{_escape(label)}</span>
        <span class="badge {status}">{status.upper()}</span>
      </div>
      <div class="stage-body">{body_html}</div>
    </div>
    """


def _escape(text):
    return (
        str(text)
        .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def render_html_report(certificate):
    """
    Render a QC certificate dict (from build_qc_certificate) into a
    single, self-contained HTML string -- no external CSS/JS/image
    assets, so the resulting file opens correctly on its own in any
    browser, including when extracted from a downloaded zip on a
    machine with no network access.
    """
    overall = certificate.get("overall_status", "passed")
    run_mode_label = "Monitor Mode" if certificate.get("run_mode") == "monitor" else "Auto Mode"
    monitor_note = f" (Monitor: {_escape(certificate['monitor_id'])})" if certificate.get("monitor_id") else ""

    banner_text = (
        "✅ PASSED -- no issues detected across any pipeline stage."
        if overall == "passed"
        else f"⚠️ FLAGGED -- {len(certificate.get('flags', []))} issue(s) detected. Review before relying on these results."
    )

    stages_html = "".join(
        _render_stage_html(stage_key, stage)
        for stage_key, stage in certificate.get("stages", {}).items()
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>QC Certificate -- {_escape(certificate.get('project_name', ''))}</title>
{_HTML_STYLE}
</head>
<body>
  <div class="container">
    <h1>QC Certificate</h1>
    <div class="subtitle">Project: <strong>{_escape(certificate.get('project_name', ''))}</strong> &middot; {run_mode_label}{monitor_note}</div>
    <div class="meta">Generated {_escape(certificate.get('generated_at', ''))}</div>
    <div class="banner {overall}">{banner_text}</div>
    {stages_html}
  </div>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Public entry point -- called by advanced_mode_orchestrator.py at the end
# of BOTH a successful and an errored run_pipeline(), so a certificate
# always exists reflecting however far the run actually got.
# ---------------------------------------------------------------------------

def save_qc_certificate(project_name, run_mode="auto", monitor_id=None):
    """
    Build and write both qc_certificate.json and qc_report.html into
    project_name's own project directory (see project_manager.py's
    qc_certificate_path()/qc_report_html_path()). Safe to call multiple
    times for the same project (e.g. once per stage-completion re-run,
    or once after a resumed run) -- always overwrites with the current,
    complete picture rather than appending.

    Returns (json_path, html_path).
    """
    certificate = build_qc_certificate(project_name, run_mode=run_mode, monitor_id=monitor_id)
    html_report = render_html_report(certificate)

    json_path = pm.qc_certificate_path(project_name)
    html_path = pm.qc_report_html_path(project_name)

    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as f:
        json.dump(certificate, f, indent=2)
    with open(html_path, "w") as f:
        f.write(html_report)

    return json_path, html_path
