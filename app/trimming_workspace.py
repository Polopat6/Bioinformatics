"""
trimming_workspace.py

Adapter Trimming & Post-Trimming Quality Control workspace.

This picks up exactly where bulk_rnaseq_workspace.py leaves off: it reuses
the *same active project* (matched samples + FastQC/MultiQC results
already on disk) rather than requiring anything to be re-uploaded.

Trimming is done with fastp, which conveniently generates its own
per-sample QC report (a .json + .html) as a side effect of trimming —
no separate FastQC run is needed post-trim. MultiQC has a built-in parser
for fastp's JSON output, so we get the same "one combined report" result
as Step 4 of the Bulk RNA-Seq workspace, just sourced from fastp instead
of FastQC.

Design goal: same as bulk_rnaseq_workspace.py — assume the user has little
to no bioinformatics background, explain each step in plain language.

This module is fully self-contained. All trimming/post-trim QC development
should happen here — editing this file has zero effect on
spatial_workspace.py or bulk_rnaseq_workspace.py.
"""

import json
import os
import shutil
import subprocess

import pandas as pd
import streamlit as st

import project_manager as pm

# Reuse the same workspace_key as bulk_rnaseq_workspace.py so the active
# project selection (st.session_state["bulk_rnaseq_project"]) is shared
# automatically between the two pages — no need to re-pick a project when
# jumping here via the "Proceed to Trimming" button.
WORKSPACE_KEY = "bulk_rnaseq"

# Known adapter sequences (or distinctive prefixes) for common library
# prep kits, used to translate fastp's raw detected adapter sequence into
# a recognizable kit name for users without a bioinformatics background.
# Matching is done on a prefix basis since fastp may report slightly
# different lengths depending on read length / overlap detection.
KNOWN_ADAPTERS = {
    "AGATCGGAAGAGC": "Illumina TruSeq / Universal Adapter",
    "CTGTCTCTTATACACATCT": "Nextera Transposase Adapter",
    "AAAAAAAAAAAAAAAAAAAAAAA": "Poly-A tail (not a true adapter — common in mRNA-seq)",
    "TGGAATTCTCGGGTGCCAAGG": "Illumina Small RNA 3' Adapter",
    "GATCGGAAGAGCACACGTCT": "Illumina TruSeq (alternate reported form)",
}

# fastp's JSON schema for the adapter_cutting section has varied slightly
# across versions, and paired-end overlap-based trimming (fastp's default
# PE method) does not always populate a specific adapter sequence field,
# since overlap trimming doesn't require knowing the adapter sequence at
# all. We try several known/possible key names defensively rather than
# assuming one fixed schema.
ADAPTER_SEQUENCE_KEYS_R1 = [
    "read1_adapter_sequence", "adapter_sequence", "adapter_sequence_r1",
]
ADAPTER_SEQUENCE_KEYS_R2 = [
    "read2_adapter_sequence", "adapter_sequence_r2",
]


def _identify_adapter(sequence):
    """
    Given a raw adapter sequence string reported by fastp, try to match
    it to a known library prep kit's adapter for a friendlier label.
    Falls back to a generic message if no match is found — fastp may
    still have correctly detected a real (just less common) adapter.
    """
    if not sequence:
        return None
    seq_upper = sequence.strip().upper()
    for known_seq, kit_name in KNOWN_ADAPTERS.items():
        if seq_upper.startswith(known_seq) or known_seq.startswith(seq_upper[:len(known_seq)]):
            return kit_name
    return "Unrecognized/custom adapter (fastp detected it automatically)"


def _extract_adapter_sequences(adapter_section):
    """
    Defensively pull R1/R2 adapter sequence strings out of fastp's
    'adapter_cutting' JSON section, trying several known key name
    variants since this has differed across fastp versions and
    detection modes (PE overlap-based trimming in particular may not
    populate any sequence field at all, even though trimming occurred).
    """
    seq_r1 = ""
    for key in ADAPTER_SEQUENCE_KEYS_R1:
        if adapter_section.get(key):
            seq_r1 = adapter_section[key]
            break

    seq_r2 = ""
    for key in ADAPTER_SEQUENCE_KEYS_R2:
        if adapter_section.get(key):
            seq_r2 = adapter_section[key]
            break

    return seq_r1, seq_r2


# ---------------------------------------------------------------------------
# Residual poly-G / poly-A tail detection
# ---------------------------------------------------------------------------
#
# fastp has two SEPARATE, independent settings for this, and neither is
# passed by _run_fastp_for_sample's default command above:
#   - "-g"/"--trim_poly_g": trims poly-G tails specifically. fastp only
#     auto-enables this on its own for reads whose header indicates a
#     NextSeq/NovaSeq-style 2-color chemistry instrument (which calls
#     "no signal" as a high-confidence G, producing artificial poly-G
#     tails) -- for any other platform (or a read header fastp doesn't
#     recognize), this never fires unless explicitly requested.
#   - "-x"/"--trim_poly_x": the general-purpose version, which also
#     catches poly-A tails (common in mRNA-seq libraries from poly-A
#     tail read-through). This is NEVER on by default -- it must be
#     explicitly requested every time.
#
# So a library with either issue can easily sail through Step 1 without
# fastp doing anything about it, then show up as an "Overrepresented
# sequences" warning in the combined MultiQC report -- often as a long
# run of G's or A's with no real adapter match. Rather than relying on
# that (comparatively coarse, sequence-frequency-based) signal, we
# detect this directly and more precisely from data fastp ALREADY wrote
# to its own JSON report as a side effect of Step 1 -- no second data
# source (e.g. a separate FastQC pass) is needed:
#
# fastp's JSON includes a "read1_after_filtering" (and, for paired-end
# data, "read2_after_filtering") section with a "content_curves" object
# -- per-cycle base-composition percentages (A/T/C/G/N) across the
# entire read length, AFTER trimming. A residual poly-G/poly-A tail
# shows up as an unmistakable spike in G% or A% concentrated in the
# LAST several cycles of the read, clearly departing from the
# ~25%-per-base baseline expected for random sequence elsewhere in the
# read. The matching "*_before_filtering" section is also read (when
# present) purely to show "this existed in the raw reads too" context
# in the flag message -- detection itself is based on the AFTER curve,
# since that's what tells us whether the issue actually survived
# trimming.
POLY_TAIL_BASES = ("G", "A")

# How many cycles at the very end of the read to treat as "the tail"
# when checking for a residual poly-G/poly-A run. 15 is a reasonable
# default for typical 75-150bp short-read data -- long enough to catch
# a real tail's full extent, short enough to not accidentally include
# a large fraction of a short read as "the tail".
POLY_TAIL_WINDOW = 15

# A tail is flagged only if BOTH conditions hold:
#   1. the tail's average composition for this base is at least this
#      fraction (0.35 = 35%) -- comfortably above the ~25% expected by
#      chance for any one of the 4 bases in random sequence, so a
#      modest, unremarkable fluctuation doesn't trigger a false alarm.
#   2. the tail's average is at least this many times higher than the
#      rest of the read's (non-tail) average for the same base -- so a
#      library that's simply naturally GC-rich (or AT-rich) throughout
#      isn't mistaken for a tail-specific artifact; the RATIO is what
#      actually identifies "something specific to the tail", not the
#      base's raw abundance.
POLY_TAIL_ABSOLUTE_FLOOR = 0.35
POLY_TAIL_RATIO_THRESHOLD = 2.0


def _normalize_fraction_curve(curve):
    """
    fastp's content_curves values are fractions (0.0-1.0) in the
    versions we've tested against, but this defensively rescales to
    fractions if a curve's values look like they're already expressed
    as percentages (0-100) instead -- so the threshold comparisons in
    _detect_residual_poly_tail_for_curve stay unit-consistent
    regardless of exactly which fastp version produced the JSON.
    """
    if not curve:
        return curve
    if max(curve) > 1.5:
        return [v / 100.0 for v in curve]
    return curve


def _get_content_curve(json_data, read_section_key, base):
    """
    Defensively extract a single base's per-cycle content curve (e.g.
    the G% value at every cycle/position along the read) from one of
    fastp's "read1_after_filtering" / "read1_before_filtering" /
    "read2_after_filtering" / "read2_before_filtering" JSON sections.

    Returns None if that section or curve isn't present at all --
    expected for "read2_*" on single-end data, and gracefully handled
    (rather than raising) for any fastp version whose JSON schema
    doesn't include content_curves for some reason.
    """
    section = json_data.get(read_section_key, {})
    if not isinstance(section, dict):
        return None
    curves = section.get("content_curves", {})
    curve = curves.get(base)
    if not curve:
        return None
    return curve


def _detect_residual_poly_tail_for_curve(curve_after, curve_before=None, tail_window=POLY_TAIL_WINDOW):
    """
    Check a single (read, base) content curve for a residual poly-tail
    signature: an abnormally high concentration of one base in the last
    `tail_window` cycles, relative to the rest of the read -- see the
    module-level POLY_TAIL_* constants above for exactly how "abnormal"
    is defined (an absolute floor AND a relative-to-body ratio, so
    neither a modest fluctuation nor a library that's simply naturally
    GC/AT-rich throughout gets flagged).

    curve_before, if given, is used only to report what this same
    tail's composition looked like in the RAW (pre-trim) reads, for
    context in the eventual user-facing message -- it does not affect
    the flagged/not-flagged decision itself, which is based purely on
    curve_after (the actual, current state of the trimmed reads).

    Returns None if curve_after is missing or too short to meaningfully
    split into "tail" vs. "body" regions (fewer than 5 cycles).
    Otherwise returns a dict:
        {
            "flagged": bool,
            "tail_pct_after": float,   # 0-100, this tail's avg composition, after trimming
            "body_pct_after": float,   # 0-100, the rest of the read's avg composition
            "ratio": float or None,    # tail_pct_after / body_pct_after
            "tail_pct_before": float or None,  # same tail region, BEFORE trimming (context only)
        }
    """
    if not curve_after or len(curve_after) < 5:
        return None

    curve_after = _normalize_fraction_curve(curve_after)
    window = min(tail_window, max(1, len(curve_after) // 3))
    tail_vals = curve_after[-window:]
    body_vals = curve_after[:-window] if len(curve_after) > window else curve_after

    tail_avg = sum(tail_vals) / len(tail_vals)
    body_avg = sum(body_vals) / len(body_vals) if body_vals else tail_avg
    ratio = (tail_avg / body_avg) if body_avg > 0 else float("inf")

    flagged = tail_avg >= POLY_TAIL_ABSOLUTE_FLOOR and ratio >= POLY_TAIL_RATIO_THRESHOLD

    result = {
        "flagged": flagged,
        "tail_pct_after": round(tail_avg * 100, 1),
        "body_pct_after": round(body_avg * 100, 1),
        "ratio": round(ratio, 2) if ratio != float("inf") else None,
        "tail_pct_before": None,
    }

    if curve_before:
        curve_before_norm = _normalize_fraction_curve(curve_before)
        tail_vals_before = curve_before_norm[-window:] if len(curve_before_norm) >= window else curve_before_norm
        if tail_vals_before:
            result["tail_pct_before"] = round(100 * sum(tail_vals_before) / len(tail_vals_before), 1)

    return result


def _detect_poly_tail_issues(json_data):
    """
    Check ONE sample's fastp JSON (already loaded from disk) for a
    residual poly-G and/or poly-A tail that survived Step 1's trimming
    -- across both R1 and R2 (if paired-end) and both bases fastp can
    specifically trim for this (see POLY_TAIL_BASES and the module
    docstring above for -g/--trim_poly_g vs. -x/--trim_poly_x).

    Returns a list of dicts (empty if nothing was flagged, or if this
    JSON doesn't have the needed content_curves data at all -- e.g. an
    older fastp version), one entry per (read, base) combination that
    was actually flagged:
        [{"read": "R1"|"R2", "base": "G"|"A", "flagged": True,
          "tail_pct_after": ..., "body_pct_after": ..., "ratio": ...,
          "tail_pct_before": ... or None}, ...]
    """
    issues = []
    for read_label, key_after, key_before in [
        ("R1", "read1_after_filtering", "read1_before_filtering"),
        ("R2", "read2_after_filtering", "read2_before_filtering"),
    ]:
        if key_after not in json_data:
            continue
        for base in POLY_TAIL_BASES:
            curve_after = _get_content_curve(json_data, key_after, base)
            curve_before = _get_content_curve(json_data, key_before, base)
            result = _detect_residual_poly_tail_for_curve(curve_after, curve_before)
            if result and result["flagged"]:
                issues.append({"read": read_label, "base": base, **result})
    return issues


def _scan_all_samples_for_poly_tail_issues(reports_dir):
    """
    Run _detect_poly_tail_issues across every *.fastp.json report
    currently on disk for this project, right after Step 1's trimming
    -- so the check reflects fastp's ACTUAL trimming result, not a
    stale/prior run.

    Returns a dict {sample_name: [issue_dict, ...]} -- only samples
    with at least one flagged (read, base) combination are included, so
    callers can simply check `if poly_tail_issues:` to decide whether to
    show anything at all.
    """
    issues_by_sample = {}
    if not os.path.isdir(reports_dir):
        return issues_by_sample

    for entry in os.listdir(reports_dir):
        if not entry.endswith(".fastp.json"):
            continue
        sample_name = entry[: -len(".fastp.json")]
        json_path = os.path.join(reports_dir, entry)
        try:
            with open(json_path) as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue

        sample_issues = _detect_poly_tail_issues(data)
        if sample_issues:
            issues_by_sample[sample_name] = sample_issues

    return issues_by_sample


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _tools_available():
    """Check whether fastp and multiqc are installed and on PATH."""
    return shutil.which("fastp") is not None, shutil.which("multiqc") is not None


def sample_already_trimmed(sample_row, trimmed_dir):
    """
    Check whether a sample's trimmed output file(s) already exist on
    disk, using the exact naming convention _run_fastp_for_sample writes
    (paired: {sample}_R1.trimmed.fastq.gz + {sample}_R2.trimmed.fastq.gz;
    single-end: {sample}.trimmed.fastq.gz).

    Used to let a re-run skip samples that already succeeded, rather
    than reprocessing (and re-writing to disk) every sample every time —
    important for larger datasets where disk space is a real constraint:
    if 14 of 16 samples already trimmed successfully and only 2 failed
    (e.g. due to running out of disk space mid-run), there's no reason
    to burn disk space and time re-writing the 14 that are already done
    just to retry the 2 that aren't.

    Returns True only if the expected output file(s) for this sample's
    layout (paired vs. single-end) are BOTH/ALL present — a partial
    write (e.g. only R1 exists because the process was killed mid-write)
    correctly returns False, so a genuinely incomplete/corrupted prior
    attempt is NOT mistaken for a successful one.
    """
    sample_name = sample_row["sample"]
    fastq_2 = sample_row.get("fastq_2", "")
    is_paired = isinstance(fastq_2, str) and fastq_2.strip() not in ("", "nan")

    if is_paired:
        out1 = os.path.join(trimmed_dir, f"{sample_name}_R1.trimmed.fastq.gz")
        out2 = os.path.join(trimmed_dir, f"{sample_name}_R2.trimmed.fastq.gz")
        return os.path.exists(out1) and os.path.exists(out2)
    else:
        out1 = os.path.join(trimmed_dir, f"{sample_name}.trimmed.fastq.gz")
        return os.path.exists(out1)


def _run_fastp_for_sample(sample_row, trimmed_dir, reports_dir, threads=3,
                           force_trim_poly_g=False, force_trim_poly_x=False):
    """
    Run fastp on a single sample's reads (paired-end or single-end),
    writing trimmed FASTQ output to trimmed_dir and a per-sample
    fastp.json/fastp.html report to reports_dir.

    `sample_row` is a row (as a dict-like) from the matched samplesheet,
    expected to have at least "sample", "fastq_1", and optionally
    "fastq_2" columns — these already contain full file paths as written
    by bulk_rnaseq_workspace.py's _write_matched_samplesheet().

    threads: fastp's own "-w"/"--thread" flag, controlling how many
    worker threads fastp uses *within processing a single sample*. Note
    this is a different axis of parallelism than running multiple
    samples at once — each call to this function still processes one
    sample fully before returning, but fastp can use multiple threads
    internally to speed up that one sample's trimming. fastp's own
    documentation notes diminishing returns (and an internal hard cap)
    beyond ~16 threads due to I/O bottlenecks rather than CPU
    availability, so this is not raised arbitrarily high by default.

    force_trim_poly_g, force_trim_poly_x: add fastp's "-g"/"--trim_poly_g"
    and/or "-x"/"--trim_poly_x" flags respectively (see the "Residual
    poly-G / poly-A tail detection" section above for what these do and
    why neither is passed by default). Both default to False, so a
    normal Step 1 run behaves exactly as before -- these are only set
    True by the Step 2 "auto re-trim flagged sample(s)" action, and only
    for the specific sample(s)/base(s) actually flagged, rather than
    unconditionally enabling extra trimming for every sample every time
    (which could clip real, biologically meaningful sequence for
    samples that don't have this issue).

    Always re-runs fastp on this sample's ORIGINAL raw FASTQ input(s) --
    never chains on top of already-trimmed output -- so a "re-trim with
    poly-tail removal" pass reflects a single, complete, reproducible
    fastp invocation (raw reads + all requested settings in one go),
    the same as a normal first-time trim, rather than compounding two
    separate trims' worth of decisions on the same reads.

    Returns (success: bool, log: str).
    """
    sample_name = sample_row["sample"]
    fastq_1 = sample_row["fastq_1"]
    fastq_2 = sample_row.get("fastq_2", "")
    is_paired = isinstance(fastq_2, str) and fastq_2.strip() not in ("", "nan")

    os.makedirs(trimmed_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    json_path = os.path.join(reports_dir, f"{sample_name}.fastp.json")
    html_path = os.path.join(reports_dir, f"{sample_name}.fastp.html")

    extra_flags = []
    if force_trim_poly_g:
        # --force_polyg_tail_trimming overrides fastp's own platform-
        # based auto-detection (which otherwise only trims poly-G on
        # reads whose header it recognizes as NextSeq/NovaSeq-style),
        # so this reliably applies poly-G trimming regardless of
        # whether fastp thinks this platform needs it.
        extra_flags.extend(["--trim_poly_g", "--force_polyg_tail_trimming"])
    if force_trim_poly_x:
        extra_flags.append("--trim_poly_x")

    if is_paired:
        out1 = os.path.join(trimmed_dir, f"{sample_name}_R1.trimmed.fastq.gz")
        out2 = os.path.join(trimmed_dir, f"{sample_name}_R2.trimmed.fastq.gz")
        cmd = [
            "fastp",
            "-i", fastq_1,
            "-I", fastq_2,
            "-o", out1,
            "-O", out2,
            "-j", json_path,
            "-h", html_path,
            "--detect_adapter_for_pe",
            "--json", json_path,
            "-w", str(threads),
        ] + extra_flags
    else:
        out1 = os.path.join(trimmed_dir, f"{sample_name}.trimmed.fastq.gz")
        cmd = [
            "fastp",
            "-i", fastq_1,
            "-o", out1,
            "-j", json_path,
            "-h", html_path,
            "-w", str(threads),
        ] + extra_flags

    try:
        result = subprocess.run(
            # Increased from the original 1 hour: with fastp threading
            # now enabled this should generally be faster, not slower,
            # but a larger safety ceiling avoids nuisance timeouts on
            # bigger real-world datasets (e.g. human RNA-seq samples)
            # that are meaningfully larger than our original small-scale
            # (E. coli) testing.
            cmd, capture_output=True, text=True, check=True, timeout=7200
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, f"fastp timed out after 2 hours on sample '{sample_name}'."


def _run_multiqc(reports_dir, multiqc_dir):
    """
    Run MultiQC on the fastp reports directory, combining every
    individual fastp JSON report into a single unified HTML report.
    MultiQC has a built-in module that automatically recognizes and
    parses fastp's JSON output format.

    Returns (success: bool, log: str).
    """
    os.makedirs(multiqc_dir, exist_ok=True)
    cmd = ["multiqc", reports_dir, "-o", multiqc_dir, "-f"]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=1800
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "MultiQC timed out after 30 minutes."


def _parse_fastp_reports(reports_dir):
    """
    Read every *.fastp.json file in reports_dir and extract a beginner-
    friendly summary: reads before/after, % passed filter, adapter
    trimming stats, and Q30 quality rate before vs. after.

    Returns a DataFrame with one row per sample.
    """
    rows = []
    if not os.path.isdir(reports_dir):
        return pd.DataFrame(rows)

    for entry in os.listdir(reports_dir):
        if not entry.endswith(".fastp.json"):
            continue
        sample_name = entry[: -len(".fastp.json")]
        json_path = os.path.join(reports_dir, entry)
        try:
            with open(json_path) as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue

        before = data.get("summary", {}).get("before_filtering", {})
        after = data.get("summary", {}).get("after_filtering", {})
        filtering = data.get("filtering_result", {})
        adapter = data.get("adapter_cutting", {})

        reads_before = before.get("total_reads", 0)
        reads_after = after.get("total_reads", 0)
        pct_reads_kept = round(100 * reads_after / reads_before, 1) if reads_before else None

        q30_before = round(before.get("q30_rate", 0) * 100, 1)
        q30_after = round(after.get("q30_rate", 0) * 100, 1)

        adapter_trimmed_reads = adapter.get("adapter_trimmed_reads", 0)

        # fastp reports the actual adapter sequence(s) it detected, when
        # available. See ADAPTER_SEQUENCE_KEYS_R1/R2 above for why this
        # is extracted defensively rather than assuming one fixed key.
        adapter_seq_r1, adapter_seq_r2 = _extract_adapter_sequences(adapter)

        adapter_kit_r1 = _identify_adapter(adapter_seq_r1) if adapter_seq_r1 else None
        adapter_kit_r2 = _identify_adapter(adapter_seq_r2) if adapter_seq_r2 else None

        if adapter_trimmed_reads == 0:
            adapter_display = "None detected"
        elif adapter_kit_r1 and adapter_kit_r2:
            adapter_display = adapter_kit_r1 if adapter_kit_r1 == adapter_kit_r2 else f"{adapter_kit_r1} (R1) / {adapter_kit_r2} (R2)"
        elif adapter_kit_r1:
            adapter_display = adapter_kit_r1
        else:
            # Adapter trimming happened (adapter_trimmed_reads > 0) but no
            # specific sequence was reported in the JSON — this is
            # expected for fastp's default PE overlap-based trimming
            # method, which trims based on read overlap position rather
            # than matching a known adapter sequence.
            adapter_display = "Adapter trimmed via overlap analysis (no specific sequence reported by fastp)"

        rows.append({
            "Sample": sample_name,
            "Reads Before": reads_before,
            "Reads After": reads_after,
            "% Reads Kept": pct_reads_kept,
            "Q30 Before": f"{q30_before}%",
            "Q30 After": f"{q30_after}%",
            "Reads w/ Adapter Trimmed": adapter_trimmed_reads,
            "Adapter Detected": adapter_display,
            "Adapter Sequence (R1)": adapter_seq_r1 or "—",
            "Adapter Sequence (R2)": adapter_seq_r2 or "—",
            "Low Quality Reads Removed": filtering.get("low_quality_reads", 0),
            "Too Short Reads Removed": filtering.get("too_short_reads", 0),
            "_raw_adapter_cutting_json": json.dumps(adapter),  # for diagnostics only
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("🧪 Adapter Trimming & Post-Trimming QC")
    st.markdown(
        "This workspace trims low-quality bases and leftover adapter "
        "sequences from your reads, then re-checks quality afterward to "
        "confirm the trimming actually helped. **No bioinformatics "
        "experience required** — follow the steps below in order."
    )
    st.markdown("---")

    # -----------------------------------------------------------------
    # Project selection — shared with the Bulk RNA-Seq workspace
    # -----------------------------------------------------------------
    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)

    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    st.markdown("---")

    # -----------------------------------------------------------------
    # Gate: require matched samples + QC to already be done
    # -----------------------------------------------------------------
    samplesheet_path = pm.samplesheet_path(project)
    qc_done = pm.has_completed_step(project, "qc_complete")

    if not os.path.exists(samplesheet_path):
        st.warning(
            "⚠️ This project doesn't have a matched sample list yet. "
            "Go to the **🧬 Bulk RNA-Seq Pipeline** page first, upload "
            "your FASTQ files and metadata, and save your matched sample "
            "list before trimming."
        )
        # Note: navigation uses the "nav_request" indirection (see
        # app.py) rather than setting st.session_state["assay_choice_radio"]
        # directly, since that widget's key is already instantiated
        # earlier in this same run.
        if st.button("⬅️ Go to Bulk RNA-Seq Pipeline", key="trim_gate1_back_btn"):
            st.session_state["nav_request"] = "🧬 Bulk RNA-Seq Pipeline"
            st.rerun()
        return

    if not qc_done:
        st.warning(
            "⚠️ Quality control hasn't been run for this project yet. "
            "We recommend reviewing raw read quality before trimming, so "
            "you know what (if anything) actually needs to be trimmed. "
            "Go to the **🧬 Bulk RNA-Seq Pipeline** page and run Step 4 "
            "(FastQC + MultiQC) first."
        )
        if st.button("⬅️ Go to Bulk RNA-Seq Pipeline", key="trim_gate2_back_btn"):
            st.session_state["nav_request"] = "🧬 Bulk RNA-Seq Pipeline"
            st.rerun()
        return

    samplesheet_df = pd.read_csv(samplesheet_path)
    st.success(f"✅ Using matched samples from project `{project}` ({len(samplesheet_df)} sample(s)).")
    st.dataframe(samplesheet_df, use_container_width=True, hide_index=True)

    trimmed_dir = pm.trimmed_fastq_dir(project)
    fastp_reports_dir = pm.fastp_reports_dir(project)
    posttrim_multiqc_dir = pm.posttrim_multiqc_dir(project)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 1: Trim adapters with fastp
    # -----------------------------------------------------------------
    st.header("Step 1: Trim Adapters & Low-Quality Bases")

    with st.expander("ℹ️ What is trimming, and why do it? (click to learn more)"):
        st.markdown(
            "**Trimming** removes parts of your reads that could "
            "interfere with downstream analysis:\n"
            "- Leftover **adapter sequences** (short synthetic sequences "
            "used during library prep that sometimes get sequenced by "
            "accident)\n"
            "- **Low-quality bases**, usually near the ends of reads\n\n"
            "Trimming isn't always necessary — if your earlier quality "
            "report came back clean, you may be able to skip straight to "
            "alignment. But if you saw adapter content warnings or poor "
            "end-of-read quality, trimming can meaningfully improve your "
            "results.\n\n"
            "This tool uses **fastp**, which automatically detects and "
            "removes adapters and trims low-quality bases using sensible "
            "defaults — no manual configuration needed. It also produces "
            "its own quality report for each sample as it trims, which "
            "we combine into a single report in Step 2."
        )

    fastp_ok, multiqc_ok = _tools_available()
    trimming_done = pm.has_completed_step(project, "trimming_complete")

    if not fastp_ok:
        st.error(
            "⚠️ fastp was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
    else:
        # Auto-detected default thread count, same pattern used for
        # FastQC in bulk_rnaseq_workspace.py — recommends a sensible
        # starting point based on this machine's actual CPU core count
        # rather than a fixed guess. fastp's own documentation notes
        # diminishing returns beyond ~16 threads (I/O-bound, not
        # CPU-bound past that point), so the slider is capped there
        # regardless of how many cores are detected.
        detected_cores, recommended_threads = pm.get_recommended_thread_count(max_default=16)
        fastp_threads = st.slider(
            "fastp threads (per sample):",
            min_value=1, max_value=16, value=min(recommended_threads, 16),
            help=(
                f"How many threads fastp uses while processing each "
                f"sample. Detected {detected_cores} CPU core(s) on this "
                f"machine, so {min(recommended_threads, 16)} is suggested. "
                "fastp itself reports diminishing returns beyond ~16 "
                "threads regardless of how many cores your machine has, "
                "since I/O (reading/writing files) becomes the limiting "
                "factor rather than CPU power at that point."
            ),
            key="fastp_threads_slider",
        )

        # --- Determine which samples still need trimming ---
        # Checked every render (not just on click) so the "already
        # trimmed" vs. "still needed" counts shown to the user are always
        # accurate and up to date, e.g. immediately reflecting a prior
        # partial run (some samples succeeded, some failed/were
        # interrupted by a disk-space error) without requiring a click
        # first just to see the current state.
        already_trimmed_rows = []
        needs_trimming_rows = []
        for _, row in samplesheet_df.iterrows():
            if sample_already_trimmed(row, trimmed_dir):
                already_trimmed_rows.append(row)
            else:
                needs_trimming_rows.append(row)

        if already_trimmed_rows and needs_trimming_rows:
            st.info(
                f"ℹ️ **{len(already_trimmed_rows)} of {len(samplesheet_df)} "
                f"sample(s)** already have trimmed output on disk and "
                f"will be **skipped** by default (no need to redo work "
                f"that already succeeded — this also avoids using disk "
                f"space re-writing files that are already there). "
                f"**{len(needs_trimming_rows)} sample(s)** still need "
                "trimming."
            )
        elif already_trimmed_rows and not needs_trimming_rows and trimming_done:
            st.success("✅ Trimming has already been run for this project — all samples have trimmed output on disk.")

        force_retrim_all = False
        if already_trimmed_rows:
            force_retrim_all = st.checkbox(
                "Force re-trim ALL samples (including ones already done)",
                value=False,
                help=(
                    "Check this only if you specifically want to redo "
                    "trimming for every sample — for example, after "
                    "changing the thread count above shouldn't require "
                    "this (thread count doesn't change the output), but "
                    "you might want this after deciding to use different "
                    "fastp settings in the future. Leaving this unchecked "
                    "(recommended) only processes samples that don't "
                    "already have trimmed output on disk."
                ),
                key="force_retrim_all_checkbox",
            )

        rows_to_process = list(samplesheet_df.iterrows()) if force_retrim_all else [
            (i, row) for i, row in samplesheet_df.iterrows() if row["sample"] in
            {r["sample"] for r in needs_trimming_rows}
        ]
        n_to_process = len(rows_to_process)

        if n_to_process == 0 and not force_retrim_all:
            trim_button_label = "✅ All Samples Already Trimmed"
            trim_button_disabled = True
        elif trimming_done or already_trimmed_rows:
            trim_button_label = f"✂️ Trim Remaining {n_to_process} Sample(s)" if not force_retrim_all else f"🔄 Force Re-trim All {n_to_process} Sample(s)"
            trim_button_disabled = False
        else:
            trim_button_label = "✂️ Trim All Matched Samples"
            trim_button_disabled = False

        trim_clicked = st.button(trim_button_label, disabled=trim_button_disabled)

        if trim_clicked:
            progress_bar = st.progress(0, text="Starting trimming...")
            all_success = True
            failed_samples = []

            for i, (_, row) in enumerate(rows_to_process):
                progress_bar.progress(
                    i / n_to_process,
                    text=f"Trimming sample {i + 1} of {n_to_process}: {row['sample']}...",
                )
                success, log = _run_fastp_for_sample(row, trimmed_dir, fastp_reports_dir, threads=fastp_threads)
                if not success:
                    all_success = False
                    failed_samples.append((row["sample"], log))

            progress_bar.progress(1.0, text="Trimming complete.")

            if failed_samples:
                st.error(f"⚠️ fastp failed on {len(failed_samples)} sample(s):")
                for sample_name, log in failed_samples:
                    with st.expander(f"Error details for {sample_name}"):
                        st.code(log)

            total_now_trimmed = len(already_trimmed_rows) + (n_to_process - len(failed_samples)) if not force_retrim_all else (n_to_process - len(failed_samples))

            if all_success:
                st.success(f"✅ Successfully trimmed {n_to_process} sample(s). ({total_now_trimmed} of {len(samplesheet_df)} total now complete.)")
                pm.mark_step_complete(project, "trimming_complete")
                trimming_done = True
            elif len(failed_samples) < n_to_process:
                st.warning(
                    f"Trimmed {n_to_process - len(failed_samples)} of "
                    f"{n_to_process} attempted sample(s) successfully "
                    f"({total_now_trimmed} of {len(samplesheet_df)} total "
                    "now complete). Review the errors above for the rest — "
                    "a common cause is running out of disk space mid-run; "
                    "free up space and click the button again to retry "
                    "just the remaining sample(s)."
                )
                pm.mark_step_complete(project, "trimming_complete")
                trimming_done = True

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 2: Post-Trimming Quality Control (via fastp's own reports)
    # -----------------------------------------------------------------
    st.header("Step 2: Post-Trimming Quality Control")

    with st.expander("ℹ️ What am I looking at here? (click to learn more)"):
        st.markdown(
            "Since fastp already generates a quality report for each "
            "sample while trimming, this step simply combines all of "
            "those individual reports into **one unified MultiQC report** "
            "— the same style of combined report you saw in Step 4 of "
            "the Bulk RNA-Seq workspace, just built from fastp's data "
            "instead of a separate FastQC run.\n\n"
            "The summary table below shows, for each sample: how many "
            "reads were kept after trimming, how many had adapters "
            "removed, and whether overall quality (Q30 rate) improved."
        )

    if not trimming_done:
        st.info("Complete Step 1 (trim your samples) to see post-trimming quality results here.")
        return

    if not multiqc_ok:
        st.error(
            "⚠️ MultiQC was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
        return

    multiqc_html_path = os.path.join(posttrim_multiqc_dir, "multiqc_report.html")
    qc_report_exists = os.path.exists(multiqc_html_path)

    qc_button_label = "🔄 Re-generate Combined Report" if qc_report_exists else "📊 Generate Combined Post-Trim Report"
    qc_clicked = st.button(qc_button_label)

    if qc_report_exists and not qc_clicked:
        st.success("✅ A combined post-trim quality report already exists for this project.")

    if qc_clicked:
        with st.spinner("Combining fastp reports with MultiQC..."):
            multiqc_success, multiqc_log = _run_multiqc(fastp_reports_dir, posttrim_multiqc_dir)

        if not multiqc_success:
            st.error("MultiQC failed to run. Details below:")
            st.code(multiqc_log)
            qc_report_exists = False
        else:
            st.success("✅ Combined post-trim quality report generated.")
            qc_report_exists = True

    if qc_report_exists:
        summary_df = _parse_fastp_reports(fastp_reports_dir)
        if not summary_df.empty:
            st.subheader("📋 Trimming Summary")

            # Show the main summary table without the raw adapter DNA
            # sequences / debug JSON — those are long, technical strings
            # that would clutter the table. The friendly "Adapter
            # Detected" column covers the everyday question of "what was
            # trimmed", and raw details are available below for anyone
            # who wants to verify or troubleshoot.
            display_cols = [
                c for c in summary_df.columns
                if c not in ("Adapter Sequence (R1)", "Adapter Sequence (R2)", "_raw_adapter_cutting_json")
            ]
            st.dataframe(summary_df[display_cols], use_container_width=True, hide_index=True)
            st.caption(
                "**% Reads Kept** close to 100% means little needed to be "
                "removed. A **Q30 After** rate higher than **Q30 Before** "
                "confirms trimming improved overall quality. **Adapter "
                "Detected** shows which library prep kit's adapter fastp "
                "found and removed, if it was able to identify one — for "
                "paired-end data, fastp often trims adapters using "
                "overlap analysis without needing to report a specific "
                "sequence, which is completely normal."
            )

            with st.expander("🔬 View raw adapter details per sample (for troubleshooting)"):
                st.markdown(
                    "These are the exact DNA sequences and raw JSON data "
                    "fastp reported for adapter detection. Most users "
                    "won't need this — it's here for verification or to "
                    "share with support if something looks wrong."
                )
                st.dataframe(
                    summary_df[["Sample", "Adapter Sequence (R1)", "Adapter Sequence (R2)"]],
                    use_container_width=True, hide_index=True,
                )
                for _, row in summary_df.iterrows():
                    st.markdown(f"**{row['Sample']}** — raw `adapter_cutting` JSON:")
                    st.code(row["_raw_adapter_cutting_json"], language="json")

            # --- Residual poly-G / poly-A tail check ---
            # See the "Residual poly-G / poly-A tail detection" section
            # near the top of this file for the full rationale: fastp
            # only auto-trims poly-G on NextSeq/NovaSeq-recognized reads,
            # and never auto-trims poly-A -- so this specific issue can
            # easily survive Step 1 untouched for other platforms/library
            # types. Checked here, using data fastp already wrote to its
            # own JSON report, rather than requiring a separate FastQC
            # pass on the trimmed reads.
            st.subheader("🧬 Poly-G / Poly-A Tail Check")
            with st.expander("ℹ️ What is this, and why does it need a separate check? (click to learn more)"):
                st.markdown(
                    "Some sequencing platforms (mainly Illumina's "
                    "2-color chemistry instruments, like NextSeq/"
                    "NovaSeq) can produce reads with a run of **G**'s at "
                    "the end when there's no real signal left to call. "
                    "Separately, mRNA-seq libraries can sometimes read "
                    "through into the transcript's actual **poly-A "
                    "tail**, producing a run of **A**'s instead. Both "
                    "are normal, well-understood artifacts -- but fastp "
                    "only automatically trims poly-G, and only when it "
                    "recognizes your read headers as coming from a "
                    "NextSeq/NovaSeq-style instrument. Poly-A is **never** "
                    "trimmed automatically. This means either issue can "
                    "slip through Step 1 untouched depending on your "
                    "platform and library type.\n\n"
                    "This check looks directly at fastp's own per-cycle "
                    "base-composition data (already generated during "
                    "Step 1 -- no extra tool needed) for an abnormal "
                    "concentration of G's or A's specifically in the "
                    "last several cycles of your trimmed reads, which is "
                    "the signature of a tail that survived trimming."
                )

            poly_tail_issues = _scan_all_samples_for_poly_tail_issues(fastp_reports_dir)

            if not poly_tail_issues:
                st.success("✅ No residual poly-G or poly-A tails detected in any sample's trimmed reads.")
            else:
                affected_samples = sorted(poly_tail_issues.keys())
                st.warning(
                    f"⚠️ **{len(affected_samples)} of {len(summary_df)} sample(s)** "
                    "still show a residual poly-G and/or poly-A tail after "
                    "trimming:"
                )
                for sample_name in affected_samples:
                    for issue in poly_tail_issues[sample_name]:
                        before_note = (
                            f" (was {issue['tail_pct_before']}% in the raw reads)"
                            if issue["tail_pct_before"] is not None else ""
                        )
                        st.caption(
                            f"- **{sample_name}** ({issue['read']}): "
                            f"**{issue['base']}** makes up **{issue['tail_pct_after']}%** "
                            f"of the last {POLY_TAIL_WINDOW} cycles, vs. "
                            f"{issue['body_pct_after']}% for the rest of the "
                            f"read{before_note}."
                        )

                st.markdown(
                    "This can be fixed by re-running fastp on just these "
                    "sample(s) with the matching poly-tail trimming option "
                    "explicitly enabled (`--trim_poly_g` for poly-G, "
                    "`--trim_poly_x` for poly-A) -- everything else about "
                    "the trim (adapter detection, quality filtering) stays "
                    "the same."
                )

                if st.button(
                    f"🔁 Auto Re-trim {len(affected_samples)} Flagged Sample(s)",
                    key="poly_tail_retrim_btn",
                ):
                    sample_row_lookup = {
                        row["sample"]: row for _, row in samplesheet_df.iterrows()
                    }
                    retrim_progress = st.progress(0, text="Starting poly-tail re-trim...")
                    retrim_failures = []

                    for i, sample_name in enumerate(affected_samples):
                        retrim_progress.progress(
                            i / len(affected_samples),
                            text=f"Re-trimming {sample_name} ({i + 1} of {len(affected_samples)})...",
                        )
                        bases_flagged = {issue["base"] for issue in poly_tail_issues[sample_name]}
                        sample_row = sample_row_lookup.get(sample_name)
                        if sample_row is None:
                            retrim_failures.append((sample_name, "Sample no longer found in the matched sample list."))
                            continue
                        success, log = _run_fastp_for_sample(
                            sample_row, trimmed_dir, fastp_reports_dir, threads=fastp_threads,
                            force_trim_poly_g=("G" in bases_flagged),
                            force_trim_poly_x=("A" in bases_flagged),
                        )
                        if not success:
                            retrim_failures.append((sample_name, log))

                    retrim_progress.progress(1.0, text="Poly-tail re-trim complete.")

                    if retrim_failures:
                        st.error(f"⚠️ Re-trimming failed for {len(retrim_failures)} sample(s):")
                        for sample_name, log in retrim_failures:
                            with st.expander(f"Error details for {sample_name}"):
                                st.code(log)
                    n_succeeded = len(affected_samples) - len(retrim_failures)
                    if n_succeeded:
                        st.success(
                            f"✅ Re-trimmed {n_succeeded} of {len(affected_samples)} "
                            "flagged sample(s) with poly-tail removal enabled."
                        )
                        # Automatically regenerate the combined MultiQC
                        # report too, so the summary table/report above
                        # reflects the fix immediately on the next
                        # rerun, rather than requiring a separate manual
                        # click of the "Re-generate Combined Report"
                        # button above just to see the result.
                        with st.spinner("Refreshing the combined report..."):
                            _run_multiqc(fastp_reports_dir, posttrim_multiqc_dir)
                        st.rerun()

        if os.path.exists(multiqc_html_path):
            st.subheader("📊 Full Combined Report")
            with open(multiqc_html_path, "rb") as f:
                report_bytes = f.read()
            st.download_button(
                "⬇️ Download Full MultiQC Report (.html)",
                data=report_bytes,
                file_name="multiqc_posttrim_report.html",
                mime="text/html",
            )
            st.caption(
                "Download the file above, then open it in your web "
                "browser to view the full interactive report with "
                "detailed charts for every sample."
            )

        st.markdown("---")
        st.success(
            f"🎉 Project `{project}` now has trimmed reads and post-"
            "trimming quality results saved. This project is ready "
            "for the next step: RNA alignment and gene counting."
        )

        if st.button("➡️ Proceed to RNA Alignment & Counts", type="primary", key="trim_proceed_align_btn"):
            # Same nav_request indirection used by the other "Proceed
            # to X" buttons in this app (see bulk_rnaseq_workspace.py
            # and app.py's module docstring for why a plain session
            # key is used here instead of directly setting
            # st.session_state["assay_choice_radio"]).
            st.session_state["nav_request"] = "🧮 RNA Alignment & Counts"
            st.rerun()
