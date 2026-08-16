"""
fastp_manager.py

Streamlit-independent fastp trimming execution and result-parsing logic
for the Bulk RNA-Seq pipeline's trimming + post-trim QC step.

Extracted out of trimming_workspace.py so this logic can be called both
from that interactive workspace AND from a non-interactive orchestrator
(Advanced Mode / Monitor Mode), which need to run the same per-sample
fastp -> MultiQC -> parse -> poly-tail-scan pipeline in the background
without any Streamlit UI calls in the way.

This module owns:
  - checking whether the fastp/multiqc CLI tools are installed
  - checking whether a sample's trimmed output already exists (resume
    support -- lets a re-run skip samples that already succeeded)
  - running fastp for a single sample (paired or single-end)
  - running MultiQC over fastp's JSON reports to build one combined
    post-trim report
  - parsing fastp's per-sample JSON into a beginner-friendly summary
    (reads before/after, adapter detected, Q30 before/after, etc.)
  - detecting residual poly-G/poly-A tails that survived trimming,
    across every sample's fastp JSON report

It does NOT know about Streamlit, project_manager step-tracking, or
anything UI-related -- callers (trimming_workspace.py today, an
Advanced/Monitor Mode orchestrator in the future) are responsible for
deciding when to call these functions, showing progress/errors, and
marking the project's "trimming_complete" step once every sample has
been trimmed successfully.
"""

import json
import os
import shutil
import subprocess

import pandas as pd

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


def identify_adapter(sequence):
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


def extract_adapter_sequences(adapter_section):
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
# passed by run_fastp_for_sample's default command below:
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
# So a library with either issue can easily sail through trimming
# without fastp doing anything about it, then show up as an
# "Overrepresented sequences" warning in the combined MultiQC report --
# often as a long run of G's or A's with no real adapter match. Rather
# than relying on that (comparatively coarse, sequence-frequency-based)
# signal, we detect this directly and more precisely from data fastp
# ALREADY wrote to its own JSON report as a side effect of trimming --
# no second data source (e.g. a separate FastQC pass) is needed:
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
# in the eventual flag message -- detection itself is based on the
# AFTER curve, since that's what tells us whether the issue actually
# survived trimming.
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


def normalize_fraction_curve(curve):
    """
    fastp's content_curves values are fractions (0.0-1.0) in the
    versions we've tested against, but this defensively rescales to
    fractions if a curve's values look like they're already expressed
    as percentages (0-100) instead -- so the threshold comparisons in
    detect_residual_poly_tail_for_curve stay unit-consistent regardless
    of exactly which fastp version produced the JSON.
    """
    if not curve:
        return curve
    if max(curve) > 1.5:
        return [v / 100.0 for v in curve]
    return curve


def get_content_curve(json_data, read_section_key, base):
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


def detect_residual_poly_tail_for_curve(curve_after, curve_before=None, tail_window=POLY_TAIL_WINDOW):
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

    curve_after = normalize_fraction_curve(curve_after)
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
        curve_before_norm = normalize_fraction_curve(curve_before)
        tail_vals_before = curve_before_norm[-window:] if len(curve_before_norm) >= window else curve_before_norm
        if tail_vals_before:
            result["tail_pct_before"] = round(100 * sum(tail_vals_before) / len(tail_vals_before), 1)

    return result


def detect_poly_tail_issues(json_data):
    """
    Check ONE sample's fastp JSON (already loaded from disk) for a
    residual poly-G and/or poly-A tail that survived trimming -- across
    both R1 and R2 (if paired-end) and both bases fastp can
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
            curve_after = get_content_curve(json_data, key_after, base)
            curve_before = get_content_curve(json_data, key_before, base)
            result = detect_residual_poly_tail_for_curve(curve_after, curve_before)
            if result and result["flagged"]:
                issues.append({"read": read_label, "base": base, **result})
    return issues


def scan_all_samples_for_poly_tail_issues(reports_dir):
    """
    Run detect_poly_tail_issues across every *.fastp.json report
    currently on disk for this project, right after trimming -- so the
    check reflects fastp's ACTUAL trimming result, not a stale/prior
    run.

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

        sample_issues = detect_poly_tail_issues(data)
        if sample_issues:
            issues_by_sample[sample_name] = sample_issues

    return issues_by_sample


# ---------------------------------------------------------------------------
# Trimming execution + report parsing
# ---------------------------------------------------------------------------

def tools_available():
    """Check whether fastp and multiqc are installed and on PATH."""
    return shutil.which("fastp") is not None, shutil.which("multiqc") is not None


def sample_already_trimmed(sample_row, trimmed_dir):
    """
    Check whether a sample's trimmed output file(s) already exist on
    disk, using the exact naming convention run_fastp_for_sample writes
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


def run_fastp_for_sample(sample_row, trimmed_dir, reports_dir, threads=3,
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
    normal trimming run behaves exactly as before -- these are only set
    True by the "auto re-trim flagged sample(s)" action, and only for
    the specific sample(s)/base(s) actually flagged, rather than
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
            # Increased from an original 1 hour: with fastp threading
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


def run_multiqc(reports_dir, multiqc_dir):
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


def parse_fastp_reports(reports_dir):
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
        adapter_seq_r1, adapter_seq_r2 = extract_adapter_sequences(adapter)

        adapter_kit_r1 = identify_adapter(adapter_seq_r1) if adapter_seq_r1 else None
        adapter_kit_r2 = identify_adapter(adapter_seq_r2) if adapter_seq_r2 else None

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
