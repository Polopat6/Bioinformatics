"""
fastqc_manager.py

Streamlit-independent FastQC/MultiQC execution and result-parsing logic
for the Bulk RNA-Seq pipeline's pre-trim QC step.

Extracted out of bulk_rnaseq_workspace.py so this logic can be called
both from that interactive workspace AND from a non-interactive
orchestrator (Advanced Mode / Monitor Mode), which need to run the same
FastQC -> MultiQC -> parse -> flag pipeline in the background without
any Streamlit UI calls in the way.

This module owns:
  - checking whether the fastqc/multiqc CLI tools are installed
  - running FastQC across a batch of FASTQ files
  - running MultiQC over FastQC's output to build one combined report
  - parsing FastQC's per-file summary.txt (bundled in each _fastqc.zip)
    into a tidy DataFrame
  - turning that raw PASS/WARN/FAIL table into a beginner-friendly
    per-sample overview with plain-language guidance

It does NOT know about Streamlit, project_manager step-tracking, or
anything UI-related -- callers (bulk_rnaseq_workspace.py today, an
Advanced/Monitor Mode orchestrator in the future) are responsible for
deciding when to call these functions, showing progress/errors, and
marking the project's "qc_complete" step once run_fastqc + run_multiqc
both succeed.
"""

import os
import shutil
import subprocess
import zipfile

import pandas as pd

# Plain-language explanations for each FastQC module, shown when a sample
# gets a WARN or FAIL for that module. Written for someone with little to
# no bioinformatics background -- explains what it means and what (if
# anything) to do about it.
FASTQC_MODULE_GUIDANCE = {
    "Per base sequence quality": (
        "Some positions in the reads have lower confidence scores. A dip "
        "at the very end of reads is common and usually fine. If quality "
        "is poor across most of the read, consider trimming low-quality "
        "ends before alignment (e.g. with fastp or Trimmomatic)."
    ),
    "Per tile sequence quality": (
        "Quality varied across the physical sequencer flowcell. This "
        "usually points to a technical issue during sequencing rather "
        "than something wrong with your sample."
    ),
    "Per sequence quality scores": (
        "A portion of reads have low overall quality. These reads may "
        "contribute noise to downstream analysis; a trimming/filtering "
        "step can remove them."
    ),
    "Per base sequence content": (
        "The proportion of A/T/G/C is uneven at certain read positions. "
        "This is expected and normal for RNA-seq (especially at the start "
        "of reads) and usually isn't a concern on its own."
    ),
    "Per sequence GC content": (
        "The GC content distribution looks different than expected. This "
        "can indicate contamination from another organism, or is simply "
        "normal variation depending on your species/library type."
    ),
    "Per base N content": (
        "Some positions have a high number of unresolved bases (N). A "
        "small amount is normal; a large amount can indicate a sequencing "
        "problem."
    ),
    "Sequence Length Distribution": (
        "Reads vary in length. This is expected if adapter trimming was "
        "already done upstream, and typically isn't a concern."
    ),
    "Sequence Duplication Levels": (
        "A high proportion of duplicate reads was detected. Some "
        "duplication is normal in RNA-seq (highly expressed genes "
        "naturally produce many identical reads), so this is less "
        "concerning than it would be for DNA sequencing."
    ),
    "Overrepresented sequences": (
        "Some exact sequences appear far more often than expected. This "
        "can be normal (highly expressed transcripts) or can indicate "
        "leftover adapter sequences or contamination — worth a quick look "
        "at the full report."
    ),
    "Adapter Content": (
        "Leftover sequencing adapter fragments were detected in the "
        "reads. It's recommended to trim adapters (e.g. with fastp or "
        "Trimmomatic) before proceeding to alignment, as adapters can "
        "interfere with accurate mapping."
    ),
}


def tools_available():
    """Check whether fastqc and multiqc are installed and on PATH."""
    return shutil.which("fastqc") is not None, shutil.which("multiqc") is not None


def run_fastqc(fastq_paths, fastqc_dir, threads=4):
    """
    Run FastQC on the given list of FASTQ file paths, writing output
    (an .html report and a _fastqc.zip per file) to the project's
    fastqc_dir.

    threads: FastQC's own built-in "-t" flag, which lets it process
    multiple files *simultaneously within this single subprocess call*
    (FastQC allocates one file per thread rather than splitting a single
    file's work across threads). Defaults to 4, a reasonable middle
    ground — high enough to meaningfully speed up multi-sample batches
    (e.g. 16 samples / 32 paired-end files, as with larger datasets like
    "airway"), without assuming the host machine has a large number of
    CPU cores available. This is a real speedup, not just a cosmetic
    flag: with only 1 effective thread (the previous implicit default),
    32 files are processed one at a time in sequence; with 4, up to 4
    files run concurrently, cutting wall-clock time roughly proportionally.

    Returns (success: bool, log: str).
    """
    os.makedirs(fastqc_dir, exist_ok=True)
    cmd = ["fastqc", "-o", fastqc_dir, "-t", str(threads)] + fastq_paths
    try:
        result = subprocess.run(
            # Increased from an original 30-minute limit that was set
            # while only testing against small bacterial (E. coli) FASTQ
            # files. Real-world datasets with many samples/large
            # human-genome-scale reads (e.g. 16 paired-end airway
            # samples = 32 files) can legitimately need well beyond 30
            # minutes even with threading, so this is bumped to 2 hours
            # as a safer ceiling that still catches genuinely stuck/hung
            # processes rather than merely slow ones.
            cmd, capture_output=True, text=True, check=True, timeout=7200
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "FastQC timed out after 2 hours."


def run_multiqc(fastqc_dir, multiqc_dir):
    """
    Run MultiQC on the FastQC output directory, combining every individual
    FastQC report into a single unified HTML report in the project's
    multiqc_dir.

    Returns (success: bool, log: str).
    """
    os.makedirs(multiqc_dir, exist_ok=True)
    cmd = ["multiqc", fastqc_dir, "-o", multiqc_dir, "-f"]
    try:
        result = subprocess.run(
            # Increased alongside FastQC's timeout above, for the same
            # reason: MultiQC scanning/combining many more (and larger)
            # FastQC reports than our original small-scale E. coli
            # testing can legitimately take longer than 30 minutes,
            # though MultiQC itself is generally much faster than FastQC
            # since it's just aggregating already-computed reports.
            cmd, capture_output=True, text=True, check=True, timeout=3600
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "MultiQC timed out after 1 hour."


def parse_fastqc_summaries(fastqc_dir):
    """
    Read the summary.txt file bundled inside each *_fastqc.zip output by
    FastQC. Each line is tab-separated: STATUS<TAB>MODULE_NAME<TAB>FILENAME
    where STATUS is one of PASS / WARN / FAIL.

    Returns a long-format DataFrame: filename, module, status.
    """
    rows = []
    if not os.path.isdir(fastqc_dir):
        return pd.DataFrame(columns=["filename", "module", "status"])

    for entry in os.listdir(fastqc_dir):
        if entry.endswith("_fastqc.zip"):
            zip_path = os.path.join(fastqc_dir, entry)
            try:
                with zipfile.ZipFile(zip_path) as zf:
                    summary_name = [n for n in zf.namelist() if n.endswith("summary.txt")]
                    if not summary_name:
                        continue
                    with zf.open(summary_name[0]) as f:
                        for line in f.read().decode("utf-8").splitlines():
                            parts = line.strip().split("\t")
                            if len(parts) == 3:
                                status, module, filename = parts
                                rows.append({
                                    "filename": filename,
                                    "module": module,
                                    "status": status,
                                })
            except zipfile.BadZipFile:
                continue

    return pd.DataFrame(rows)


def build_quality_flags(summary_df):
    """
    Convert the raw per-module PASS/WARN/FAIL table into a beginner-
    friendly per-sample overview: an overall status badge, plus plain-
    language explanations for anything that wasn't a clean PASS.
    """
    if summary_df.empty:
        return pd.DataFrame(columns=["File", "Overall Quality", "Details"]), {}

    overview_rows = []
    details_by_file = {}

    for filename, group in summary_df.groupby("filename"):
        has_fail = (group["status"] == "FAIL").any()
        has_warn = (group["status"] == "WARN").any()

        if has_fail:
            overall = "🔴 Needs attention"
        elif has_warn:
            overall = "🟡 Minor issues"
        else:
            overall = "🟢 Good quality"

        flagged = group[group["status"] != "PASS"]
        explanations = []
        for _, row in flagged.iterrows():
            icon = "🔴" if row["status"] == "FAIL" else "🟡"
            guidance = FASTQC_MODULE_GUIDANCE.get(
                row["module"],
                "No additional guidance available for this check."
            )
            explanations.append(f"{icon} **{row['module']}**: {guidance}")

        overview_rows.append({
            "File": filename,
            "Overall Quality": overall,
            "Details": f"{len(flagged)} item(s) flagged" if flagged.shape[0] > 0 else "All checks passed",
        })
        details_by_file[filename] = explanations

    return pd.DataFrame(overview_rows), details_by_file
