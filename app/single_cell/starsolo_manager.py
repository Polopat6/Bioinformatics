"""
single_cell/starsolo_manager.py

Builds AND executes STARsolo (and, as an optional alternate, salmon
alevin-fry) commands for the Single-cell RNA-Seq pipeline. Reuses the
EXISTING genome index infrastructure -- STARsolo uses the identical
genome index as ordinary STAR, confirmed directly against STAR's own
official STARsolo.md documentation: "The genome index is the same as
for normal STAR runs." So this module does NOT duplicate
reference_manager.py's download/index-build logic, or
quantification_manager.py's build_star_index()/star_index_exists() --
it calls those directly (see singlecell_workspace.py's Step 6, which
reuses qm.build_star_index/qm.star_index_exists exactly as
alignment_workspace.py does for bulk).

--soloCellFilter EmptyDrops_CR is STARsolo's own real, documented flag
(confirmed against STAR's official STARsolo.md: "cell filtering
(calling) nearly identical to that of CellRanger 3 and 4") -- used as
the recommended default. CellRanger2.2 (simpler knee-point method) is
offered as the more transparent alternative.

--- Execution wiring (2026-08-17) ---
run_starsolo_align() mirrors quantification_manager.py's real
run_star_align() almost exactly (same subprocess.run pattern: capture_
output=True, text=True, check=True, a 2-hour timeout, returning
(success, log, sample_output_dir)) -- confirmed against that module's
actual source so this module's execution behaves identically to bulk's
STAR execution rather than inventing a different convention. The one
STARsolo-specific difference: build_starsolo_command() already sets
--outSAMtype BAM Unsorted (not SortedByCoordinate), so
quantification_manager.py's compute_safe_bam_sorting_bins() /
--outBAMsortingBinsN (which only matters for SortedByCoordinate output)
does not apply here and is deliberately not used.

--- Summary.csv parsing (2026-08-17) ---
STARsolo writes a per-sample Solo.out/Gene/Summary.csv (confirmed real
via STAR's own STARsolo output structure, cellgeni/STARsolo's own
documentation, and multiple independent real-world usage guides) as a
plain, headerless, two-column (metric_name,value) CSV. This module could
NOT independently confirm the EXACT metric name STARsolo uses for "how
many cells were called" (candidates seen referenced across sources
include phrasings like "Estimated Number of Cells") -- rather than
hardcode one unverified exact string and risk it silently matching
nothing on a real file (as happened earlier this project with an
unverified Ensembl URL and an unverified fastqc_manager.py return shape),
parse_starsolo_summary() below parses the file generically into a plain
dict (safe, since the headerless two-column CSV structure itself IS
confirmed from multiple independent sources), and get_cells_detected()
checks several plausible real key spellings defensively -- the same
defensive-multiple-candidate-keys pattern fastp_manager.py's own
ADAPTER_SEQUENCE_KEYS_R1/R2 already uses for exactly this kind of
not-fully-confirmed-schema situation. If none of the candidate keys
match on a real run, this returns None rather than silently guessing --
callers must handle that case (see singlecell_workspace.py's Step 6).
"""
import csv
import os
import subprocess

CELL_FILTER_OPTIONS = {
    "EmptyDrops_CR": {
        "label": "EmptyDrops (recommended -- matches Cell Ranger's method)",
        "starsolo_flag": "EmptyDrops_CR",
        "explanation": (
            "Statistically models the 'ambient' background noise level and flags any droplet "
            "whose RNA content is significantly above that background as a real cell. Same "
            "method Cell Ranger uses -- recommended unless you have a specific reason to want "
            "a simpler, more manually-inspectable threshold instead."
        ),
    },
    "CellRanger2.2": {
        "label": "Simple UMI-count threshold (knee-point method)",
        "starsolo_flag": "CellRanger2.2",
        "explanation": (
            "Ranks all detected barcodes by total UMI count and looks for the 'knee' in that "
            "curve -- simpler and more transparent, but less statistically rigorous and more "
            "prone to mistakes with unusual cell-size distributions or heavy ambient contamination."
        ),
    },
}
DEFAULT_CELL_FILTER = "EmptyDrops_CR"

# Candidate key spellings for STARsolo's Summary.csv "cells called" metric
# -- see module docstring's "Summary.csv parsing" section for why this is
# a defensive list rather than one hardcoded (unverified) exact string.
_CELLS_DETECTED_KEY_CANDIDATES = (
    "Estimated Number of Cells",
    "Estimated Number of Cells,",
    "Number of Cells",
)


def star_tool_available():
    "Check whether STAR is installed and on PATH -- same check quantification_manager.py's star_tool_available() performs, duplicated here so this module doesn't require importing that one just for this."
    import shutil
    return shutil.which("STAR") is not None


def build_starsolo_command(genome_dir, r1_files, r2_files, whitelist_path, cb_len, umi_len,
                            output_prefix, cell_filter=DEFAULT_CELL_FILTER, threads=4, strand="forward"):
    if cell_filter not in CELL_FILTER_OPTIONS:
        raise ValueError(f"Unknown cell_filter '{cell_filter}' -- must be one of {list(CELL_FILTER_OPTIONS)}")
    r2_arg = ",".join(r2_files)
    r1_arg = ",".join(r1_files)
    wl_arg = whitelist_path if whitelist_path else "None"
    strand_flag = "Forward" if strand == "forward" else "Reverse"
    return [
        "STAR", "--runMode", "alignReads", "--genomeDir", genome_dir,
        "--readFilesIn", r2_arg, r1_arg,  # R2 first, R1 second -- STARsolo's own convention
        "--readFilesCommand", "zcat",
        "--soloType", "CB_UMI_Simple",
        "--soloCBwhitelist", wl_arg,
        "--soloCBstart", "1", "--soloCBlen", str(cb_len),
        "--soloUMIstart", str(cb_len + 1), "--soloUMIlen", str(umi_len),
        "--soloStrand", strand_flag,
        "--soloCellFilter", CELL_FILTER_OPTIONS[cell_filter]["starsolo_flag"],
        "--soloFeatures", "Gene",
        "--outSAMtype", "BAM", "Unsorted",
        "--runThreadN", str(threads),
        "--outFileNamePrefix", output_prefix,
    ]


def run_starsolo_align(sample_name, genome_dir, r1_files, r2_files, whitelist_path, cb_len, umi_len,
                        output_base_dir, cell_filter=DEFAULT_CELL_FILTER, threads=4, strand="forward"):
    """
    Run STARsolo alignment/quantification/cell-calling for a single
    sample, writing output to output_base_dir/<sample_name>/.

    Mirrors quantification_manager.run_star_align()'s exact execution
    pattern (subprocess.run with capture_output/text/check=True, a
    2-hour timeout) -- confirmed against that module's real source --
    so single-cell alignment behaves identically to bulk's STAR
    execution, just building the STARsolo-specific command via
    build_starsolo_command() above instead of a plain STAR alignment
    command.

    Returns (success: bool, log: str, sample_output_dir: str) -- same
    3-tuple shape as run_star_align()/run_salmon_quant(), so
    singlecell_workspace.py's Step 6 can display results using the
    identical pattern alignment_workspace.py already uses for bulk.
    """
    sample_out_dir = os.path.join(output_base_dir, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)
    out_prefix = os.path.join(sample_out_dir, f"{sample_name}_")

    cmd = build_starsolo_command(
        genome_dir, r1_files, r2_files, whitelist_path, cb_len, umi_len,
        output_prefix=out_prefix, cell_filter=cell_filter, threads=threads, strand=strand,
    )
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)
        return True, result.stdout + result.stderr, sample_out_dir
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or ""), sample_out_dir
    except subprocess.TimeoutExpired:
        return False, f"STARsolo alignment timed out after 2 hours on sample '{sample_name}'.", sample_out_dir


def build_alevin_fry_commands(salmon_index_dir, r1_files, r2_files, chemistry_protocol_key, output_dir, threads=4):
    alevin_map_dir = os.path.join(output_dir, "alevin_map")
    quant_dir = os.path.join(output_dir, "quant")
    map_cmd = [
        "salmon", "alevin", "-l", "ISR", "-i", salmon_index_dir,
        "-1", *r1_files, "-2", *r2_files,
        "--chromiumV3" if chemistry_protocol_key == "10XV3" else "--chromium",
        "-p", str(threads), "-o", alevin_map_dir, "--sketch",
    ]
    generate_permit_cmd = ["alevin-fry", "generate-permit-list", "-d", "fw", "-i", alevin_map_dir, "-o", quant_dir, "--unfiltered-pl"]
    collate_cmd = ["alevin-fry", "collate", "-i", quant_dir, "-r", alevin_map_dir, "-t", str(threads)]
    quant_cmd = ["alevin-fry", "quant", "-i", quant_dir, "-m", os.path.join(salmon_index_dir, "t2g.tsv"), "-r", "cr-like", "-o", quant_dir, "-t", str(threads)]
    return [map_cmd, generate_permit_cmd, collate_cmd, quant_cmd]


def counts_matrix_dir(starsolo_output_prefix):
    return os.path.join(f"{starsolo_output_prefix}Solo.out", "Gene", "raw")


def filtered_counts_matrix_dir(starsolo_output_prefix):
    return os.path.join(f"{starsolo_output_prefix}Solo.out", "Gene", "filtered")


def summary_csv_path(starsolo_output_prefix):
    return os.path.join(f"{starsolo_output_prefix}Solo.out", "Gene", "Summary.csv")


def parse_starsolo_summary(sample_output_dir, sample_name):
    """
    Parse STARsolo's Solo.out/Gene/Summary.csv for a sample into a plain
    {metric_name: value_string} dict. The file's headerless, two-column
    (metric,value) CSV structure is confirmed real (multiple independent
    sources describe/show this exact shape) -- but this module could NOT
    independently confirm the exact metric-name strings STARsolo uses
    for every row, so this function deliberately returns the RAW parsed
    dict rather than trying to normalize/rename keys itself. See
    get_cells_detected() below for how a specific metric is looked up
    defensively from this dict.

    Returns {} if the file doesn't exist or can't be parsed -- callers
    should treat an empty dict as "summary not available" rather than
    assuming success.
    """
    out_prefix = os.path.join(sample_output_dir, f"{sample_name}_")
    csv_path = summary_csv_path(out_prefix)
    if not os.path.isfile(csv_path):
        return {}
    result = {}
    try:
        with open(csv_path, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) >= 2:
                    result[row[0].strip()] = row[1].strip()
    except (OSError, csv.Error):
        return {}
    return result


def get_cells_detected(summary_dict):
    """
    Defensively look up "how many cells were called" from a
    parse_starsolo_summary() dict, checking several plausible real key
    spellings (see module docstring) rather than assuming one exact,
    unverified string. Returns an int, or None if no candidate key was
    found or its value wasn't a valid integer -- callers must handle
    None explicitly (e.g. by showing "cell count unavailable" rather
    than a fabricated number).
    """
    for key in _CELLS_DETECTED_KEY_CANDIDATES:
        if key in summary_dict:
            try:
                return int(float(summary_dict[key]))
            except (ValueError, TypeError):
                continue
    return None
