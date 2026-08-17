"""
single_cell/singlecell_ingestion_manager.py

FASTQ discovery/pairing and sample metadata handling for the Single-cell
RNA-Seq pipeline's Step 1 (Ingestion). Deliberately separate from the Bulk
RNA-Seq ingestion_manager.py -- single-cell FASTQ pairing is NOT symmetric
the way bulk paired-end reads are: R1 (cell barcode + UMI, short/fixed
length) and R2 (cDNA, the actual biological read) play fundamentally
different roles, and this module's find_r1_r2_pairs() reflects that by
returning them as distinct, clearly-labeled roles rather than a generic
"paired-end mate 1/mate 2" pair.

--- FASTQ naming convention ---
Follows the standard 10x/Illumina bcl2fastq convention (confirmed against
10x Genomics' own documentation), which every chemistry in
chemistry_manager.py's CHEMISTRY_CATALOG is generated under:
    [Sample Name]_S1_L00[Lane]_[R1|R2]_001.fastq.gz
    [Sample Name]_S1_[R1|R2]_001.fastq.gz   (no-lane variant)
Multiple lanes for the same sample (e.g. _L001_, _L002_) are grouped
together under that sample's name -- STARsolo natively accepts a
comma-separated list of per-lane FASTQs for the same sample, so lanes are
NOT concatenated by this code; they're passed through as-is and let
STARsolo/alevin-fry handle them directly.

--- Sample-level, not cell-level, metadata ---
metadata.csv here is exactly analogous to the Bulk RNA-Seq pipeline's
metadata: one row per SAMPLE (condition/treatment/donor), not per cell.
Per-cell metadata (cluster assignment, cell type, QC-pass/fail) is
generated automatically downstream, by the pipeline itself, once cells
have actually been identified -- it is never something a user uploads at
this stage.
"""
import os
import re

import pandas as pd

# Matches: <sample>_S<N>_L<lane>_<R1|R2>_001.fastq(.gz)  OR  <sample>_S<N>_<R1|R2>_001.fastq(.gz)
_FASTQ_PATTERN = re.compile(
    r"^(?P<sample>.+?)_S\d+(?:_L(?P<lane>\d+))?_(?P<read>R1|R2)_001\.f(ast)?q(\.gz)?$"
)


def find_fastq_files_in_directory(directory):
    "Return a sorted list of every *_R1_001/*_R2_001-style FASTQ filename found directly inside directory (non-recursive, matching find_fastq_filenames_in_directory's behavior in the bulk ingestion_manager.py)."
    if not os.path.isdir(directory):
        return []
    return sorted([
        f for f in os.listdir(directory)
        if _FASTQ_PATTERN.match(f)
    ])


def find_r1_r2_pairs(directory):
    """
    Scan directory for FASTQ files following the standard 10x/bcl2fastq
    naming convention and group them by sample name into R1/R2 pairs
    (grouping multiple lanes together under the same sample, each lane
    kept as a SEPARATE file path rather than concatenated).

    Returns a dict:
        {
            "<sample_name>": {
                "r1_files": [path, ...],  # sorted by lane
                "r2_files": [path, ...],  # sorted by lane, same order as r1_files
            },
            ...
        }

    Any file that doesn't match the expected naming convention at all is
    silently skipped here -- callers should separately surface
    find_unmatched_files() output to the user so an unexpectedly-named
    file isn't just invisibly dropped without explanation.
    """
    samples = {}
    for fname in find_fastq_files_in_directory(directory):
        m = _FASTQ_PATTERN.match(fname)
        sample = m.group("sample")
        lane = int(m.group("lane")) if m.group("lane") else 0
        read = m.group("read")
        entry = samples.setdefault(sample, {"r1_files": {}, "r2_files": {}})
        entry["r1_files" if read == "R1" else "r2_files"][lane] = os.path.join(directory, fname)

    result = {}
    for sample, entry in samples.items():
        lanes_sorted = sorted(set(entry["r1_files"].keys()) | set(entry["r2_files"].keys()))
        result[sample] = {
            "r1_files": [entry["r1_files"][l] for l in lanes_sorted if l in entry["r1_files"]],
            "r2_files": [entry["r2_files"][l] for l in lanes_sorted if l in entry["r2_files"]],
        }
    return result


def find_unmatched_files(directory):
    "Files directly inside directory that look like FASTQ (.fastq/.fq, possibly .gz) but don't match the expected 10x/bcl2fastq naming convention -- surfaced to the user rather than silently ignored, since an unexpectedly-named file is a common real-world mistake (e.g. a stray downloaded file, or a non-10x naming scheme)."
    if not os.path.isdir(directory):
        return []
    fastq_like = re.compile(r"\.f(ast)?q(\.gz)?$", re.IGNORECASE)
    return sorted([
        f for f in os.listdir(directory)
        if fastq_like.search(f) and not _FASTQ_PATTERN.match(f)
    ])


def validate_pairs(pairs):
    """
    Check every sample's R1/R2 lane lists for consistency (same lane
    count on both sides -- a sample with 2 R1 lanes but only 1 R2 lane
    indicates a missing/corrupted file, not a valid multi-lane sample).

    Returns a dict {sample_name: warning_message} for any sample with a
    mismatch -- empty dict if everything lines up. Does NOT raise; the
    caller (singlecell_workspace.py) is expected to surface these as
    warnings and let the user decide whether to proceed, exclude, or fix
    the underlying files.
    """
    warnings = {}
    for sample, entry in pairs.items():
        n_r1, n_r2 = len(entry["r1_files"]), len(entry["r2_files"])
        if n_r1 == 0 or n_r2 == 0:
            warnings[sample] = f"Missing {'R1' if n_r1 == 0 else 'R2'} file(s) entirely -- this sample cannot be processed."
        elif n_r1 != n_r2:
            warnings[sample] = f"Mismatched lane counts: {n_r1} R1 file(s) but {n_r2} R2 file(s) -- check for a missing or extra file."
    return warnings


def read_metadata_file(uploaded_file):
    """
    Parse an uploaded sample-level metadata file (.csv/.txt/.xlsx/.xls).
    Returns (dataframe, error_message) -- exactly mirroring the bulk
    ingestion_manager.py's read_metadata_file() return shape, so
    singlecell_workspace.py's error-handling code looks identical to the
    bulk workspace's.
    """
    name = getattr(uploaded_file, "name", "")
    try:
        if name.lower().endswith((".xlsx", ".xls")):
            df = pd.read_excel(uploaded_file)
        else:
            df = pd.read_csv(uploaded_file)
        return df, None
    except Exception as e:  # noqa: BLE001 -- surface any parse failure directly to the UI
        return None, f"Could not read metadata file: {e}"
