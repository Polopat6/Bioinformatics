"""
quantification_manager.py

Handles the actual index-building and quantification/alignment execution
for the RNA Alignment & Counts workspace (alignment_workspace.py), for
both Salmon and STAR.

Kept as its own module (separate from alignment_workspace.py) since this
involves subprocess execution, file layout detection, and log/output
parsing that is logically distinct from the Streamlit UI/workflow code —
same reasoning as reference_manager.py.
"""

import json
import os
import re
import shutil
import subprocess


# ---------------------------------------------------------------------------
# Paired-end vs. single-end detection (from trimmed FASTQ file layout)
# ---------------------------------------------------------------------------

def detect_trimmed_read_layout(trimmed_dir, sample_name):
    """
    Determine whether a sample's *trimmed* reads are paired-end or
    single-end by checking which files actually exist on disk in
    trimmed_dir, using the naming convention trimming_workspace.py
    writes:
        paired:     {sample}_R1.trimmed.fastq.gz + {sample}_R2.trimmed.fastq.gz
        single-end: {sample}.trimmed.fastq.gz

    We check the actual files on disk (not the original samplesheet's
    fastq_2 column) since that's the ground truth of what will actually
    be fed into Salmon/STAR — this also makes detection self-correcting
    if trimming happened to only produce one file for some reason.

    Returns a dict: {"read_type": "paired"|"single"|"missing",
                      "r1": path or None, "r2": path or None}
    """
    r1_path = os.path.join(trimmed_dir, f"{sample_name}_R1.trimmed.fastq.gz")
    r2_path = os.path.join(trimmed_dir, f"{sample_name}_R2.trimmed.fastq.gz")
    se_path = os.path.join(trimmed_dir, f"{sample_name}.trimmed.fastq.gz")

    if os.path.exists(r1_path) and os.path.exists(r2_path):
        return {"read_type": "paired", "r1": r1_path, "r2": r2_path}
    elif os.path.exists(se_path):
        return {"read_type": "single", "r1": se_path, "r2": None}
    elif os.path.exists(r1_path) and not os.path.exists(r2_path):
        return {"read_type": "missing_mate", "r1": r1_path, "r2": None}
    else:
        return {"read_type": "missing", "r1": None, "r2": None}


def build_sample_manifest(samplesheet_df, trimmed_dir):
    """
    Build a per-sample manifest describing the detected read layout for
    every sample in the matched samplesheet, based on actual trimmed
    files present on disk.

    Returns a list of dicts, one per sample:
        {"sample": ..., "read_type": ..., "r1": ..., "r2": ...}
    """
    manifest = []
    for _, row in samplesheet_df.iterrows():
        sample_name = row["sample"]
        layout = detect_trimmed_read_layout(trimmed_dir, sample_name)
        manifest.append({"sample": sample_name, **layout})
    return manifest


# ---------------------------------------------------------------------------
# Salmon: indexing
# ---------------------------------------------------------------------------

def salmon_tool_available():
    return shutil.which("salmon") is not None


def salmon_index_exists(index_dir):
    """
    Check whether a Salmon index already exists at index_dir. Salmon
    writes an "info.json" file into the index directory on successful
    completion, which we use as the signal that indexing finished
    (rather than just checking the directory exists, which could be
    true even if a previous indexing attempt failed partway through).
    """
    return os.path.exists(os.path.join(index_dir, "info.json"))


def build_salmon_index(transcriptome_fasta, index_dir, threads=4):
    """
    Build a Salmon index from a transcriptome FASTA.

    Returns (success: bool, log: str).
    """
    os.makedirs(index_dir, exist_ok=True)
    cmd = [
        "salmon", "index",
        "-t", transcriptome_fasta,
        "-i", index_dir,
        "-p", str(threads),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=3600)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "Salmon indexing timed out after 1 hour."


# ---------------------------------------------------------------------------
# Salmon: quantification
# ---------------------------------------------------------------------------

def run_salmon_quant(sample_entry, index_dir, output_base_dir, threads=4):
    """
    Run Salmon quantification for a single sample (paired or single-end),
    writing output to output_base_dir/<sample_name>/.

    Returns (success: bool, log: str, sample_output_dir: str).
    """
    sample_name = sample_entry["sample"]
    sample_out_dir = os.path.join(output_base_dir, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)

    base_cmd = ["salmon", "quant", "-i", index_dir, "-l", "A", "-p", str(threads), "-o", sample_out_dir]

    if sample_entry["read_type"] == "paired":
        cmd = base_cmd + ["-1", sample_entry["r1"], "-2", sample_entry["r2"]]
    elif sample_entry["read_type"] == "single":
        cmd = base_cmd + ["-r", sample_entry["r1"]]
    else:
        return False, f"Sample '{sample_name}' has an incomplete/missing trimmed read layout ({sample_entry['read_type']}); skipped.", sample_out_dir

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=3600)
        return True, result.stdout + result.stderr, sample_out_dir
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or ""), sample_out_dir
    except subprocess.TimeoutExpired:
        return False, f"Salmon quant timed out after 1 hour on sample '{sample_name}'.", sample_out_dir


def parse_salmon_mapping_rate(sample_output_dir):
    """
    Read Salmon's aux_info/meta_info.json to extract the overall mapping
    rate (% of reads that mapped to the transcriptome) for a sample.

    Returns a float percentage, or None if the file couldn't be read.
    """
    meta_path = os.path.join(sample_output_dir, "aux_info", "meta_info.json")
    if not os.path.exists(meta_path):
        return None
    try:
        with open(meta_path) as f:
            data = json.load(f)
        return round(data.get("percent_mapped", 0), 2)
    except (json.JSONDecodeError, OSError):
        return None


def read_salmon_quant_sf(sample_output_dir):
    """
    Read a sample's quant.sf output file (Salmon's per-transcript
    quantification table). Returns a pandas DataFrame, or None if the
    file doesn't exist.

    Columns: Name, Length, EffectiveLength, TPM, NumReads
    """
    import pandas as pd
    quant_path = os.path.join(sample_output_dir, "quant.sf")
    if not os.path.exists(quant_path):
        return None
    return pd.read_csv(quant_path, sep="\t")


# ---------------------------------------------------------------------------
# STAR: indexing
# ---------------------------------------------------------------------------

def star_tool_available():
    return shutil.which("STAR") is not None


def star_index_exists(index_dir):
    """
    Check whether a STAR genome index already exists. STAR writes a
    "SAindex" file (among others) into the index directory on successful
    completion.
    """
    return os.path.exists(os.path.join(index_dir, "SAindex"))


def build_star_index(genome_fasta, gtf_path, index_dir, threads=4, sjdb_overhang=100):
    """
    Build a STAR genome index from a genome FASTA + GTF annotation.

    sjdb_overhang should ideally be set to (read_length - 1); 100 is a
    reasonable default for typical 100-150bp reads and rarely matters
    much in practice.

    Returns (success: bool, log: str).
    """
    os.makedirs(index_dir, exist_ok=True)
    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--genomeDir", index_dir,
        "--genomeFastaFiles", genome_fasta,
        "--sjdbGTFfile", gtf_path,
        "--sjdbOverhang", str(sjdb_overhang),
        "--runThreadN", str(threads),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "STAR genome indexing timed out after 2 hours."


# ---------------------------------------------------------------------------
# STAR: alignment + built-in gene counting
# ---------------------------------------------------------------------------

def run_star_align(sample_entry, index_dir, output_base_dir, threads=4):
    """
    Run STAR alignment for a single sample (paired or single-end), using
    --quantMode GeneCounts to get gene-level counts directly from STAR
    without needing a separate counting tool (e.g. featureCounts).

    Output includes:
        <sample>_Aligned.sortedByCoord.out.bam
        <sample>_ReadsPerGene.out.tab   <- gene-level counts
        <sample>_Log.final.out         <- mapping statistics

    Returns (success: bool, log: str, sample_output_dir: str).
    """
    sample_name = sample_entry["sample"]
    sample_out_dir = os.path.join(output_base_dir, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)
    out_prefix = os.path.join(sample_out_dir, f"{sample_name}_")

    if sample_entry["read_type"] == "paired":
        read_files_arg = [sample_entry["r1"], sample_entry["r2"]]
    elif sample_entry["read_type"] == "single":
        read_files_arg = [sample_entry["r1"]]
    else:
        return False, f"Sample '{sample_name}' has an incomplete/missing trimmed read layout ({sample_entry['read_type']}); skipped.", sample_out_dir

    cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--genomeDir", index_dir,
        "--readFilesIn", *read_files_arg,
        "--readFilesCommand", "zcat",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--quantMode", "GeneCounts",
        "--outFileNamePrefix", out_prefix,
        "--runThreadN", str(threads),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)
        return True, result.stdout + result.stderr, sample_out_dir
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or ""), sample_out_dir
    except subprocess.TimeoutExpired:
        return False, f"STAR alignment timed out after 2 hours on sample '{sample_name}'.", sample_out_dir


def parse_star_mapping_rate(sample_output_dir, sample_name):
    """
    Parse STAR's Log.final.out to extract the uniquely-mapped read
    percentage for a sample.

    Returns a float percentage, or None if the file couldn't be read.
    """
    log_path = os.path.join(sample_output_dir, f"{sample_name}_Log.final.out")
    if not os.path.exists(log_path):
        return None
    try:
        with open(log_path) as f:
            content = f.read()
        match = re.search(r"Uniquely mapped reads % \|\s+([\d.]+)%", content)
        return float(match.group(1)) if match else None
    except OSError:
        return None


def read_star_gene_counts(sample_output_dir, sample_name, strandedness_column=1):
    """
    Read a sample's ReadsPerGene.out.tab (STAR's built-in gene counts
    output). This file has 4 columns per gene: gene_id, unstranded
    count, forward-stranded count, reverse-stranded count. The first 4
    rows are summary stats (N_unmapped, N_multimapping, etc.) rather
    than genes, and are skipped.

    strandedness_column: which count column to use —
        1 = unstranded (default, safest choice when library prep
            strandedness is unknown)
        2 = forward-stranded
        3 = reverse-stranded

    Returns a pandas DataFrame with columns ["gene_id", "count"], or
    None if the file doesn't exist.
    """
    import pandas as pd
    counts_path = os.path.join(sample_output_dir, f"{sample_name}_ReadsPerGene.out.tab")
    if not os.path.exists(counts_path):
        return None

    df = pd.read_csv(counts_path, sep="\t", header=None, skiprows=4,
                      names=["gene_id", "unstranded", "forward", "reverse"])
    col_map = {1: "unstranded", 2: "forward", 3: "reverse"}
    df = df[["gene_id", col_map[strandedness_column]]].rename(
        columns={col_map[strandedness_column]: "count"}
    )
    return df


# ---------------------------------------------------------------------------
# tximport: gene-level collapsing of Salmon's transcript-level output
# ---------------------------------------------------------------------------
#
# Salmon quantifies at the transcript level. For organisms with no
# introns (e.g. bacteria), each gene typically has exactly one
# transcript, so transcript-level counts are already gene-level. For
# eukaryotes with multiple transcripts per gene, transcript-level counts
# need to be collapsed to gene-level using a transcript-to-gene (tx2gene)
# mapping — this is exactly what tximport does, and is the standard,
# recommended approach for feeding Salmon output into DESeq2. See
# reference_manager.py's extract_tx2gene_from_ensembl_fasta /
# extract_tx2gene_from_gtf / build_identity_tx2gene + save_tx2gene_csv
# for how the tx2gene mapping itself is built for each reference source.

# This R script is written to a temp file and executed via Rscript. It
# takes three positional arguments: a tx2gene CSV path (columns TXNAME,
# GENEID), a sample manifest CSV path (columns sample, path — one row
# per sample pointing at that sample's quant.sf file), and an output
# path for the resulting gene-level counts matrix CSV.
#
# tximport's default countsFromAbundance="no" setting returns the raw
# summed read counts per gene (the appropriate input for DESeq2), rather
# than a length-bias-corrected variant — this matches what DESeq2's own
# documentation recommends when importing via tximport without also
# supplying an average transcript length offset matrix.
#
# dropInfReps=TRUE skips loading Salmon's inferential replicate data
# (bootstrap/Gibbs samples used for uncertainty-aware DE methods like
# swish). Without this flag, tximport requires the R "jsonlite" package
# to parse those replicate files even when a user only wants standard
# gene-level counts for DESeq2 — an easy-to-miss extra dependency that
# has nothing to do with the actual counts being requested here. Setting
# dropInfReps=TRUE avoids that dependency entirely for this standard-DE
# use case.
_TXIMPORT_R_SCRIPT = r'''
suppressMessages(library(tximport))

args <- commandArgs(trailingOnly = TRUE)
tx2gene_path <- args[1]
samples_manifest_path <- args[2]
output_path <- args[3]

tx2gene <- read.csv(tx2gene_path, stringsAsFactors = FALSE)
samples <- read.csv(samples_manifest_path, stringsAsFactors = FALSE)

files <- samples$path
names(files) <- samples$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop(paste("Missing quant.sf file(s):", paste(missing_files, collapse = ", ")))
}

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no", dropInfReps = TRUE)

counts_df <- as.data.frame(round(txi$counts))
counts_df$gene_id <- rownames(counts_df)
counts_df <- counts_df[, c("gene_id", samples$sample)]

write.csv(counts_df, output_path, row.names = FALSE)
cat("tximport gene-level collapsing completed successfully.\n")
cat(paste("Genes:", nrow(counts_df), "| Samples:", length(samples$sample)), "\n")
'''


def tximport_available():
    """Check whether Rscript (and by extension, R + tximport) is available."""
    return shutil.which("Rscript") is not None


def run_tximport_gene_collapse(sample_quant_paths, tx2gene_path, output_path, work_dir):
    """
    Run tximport (via an Rscript subprocess) to collapse Salmon's
    transcript-level quant.sf output into gene-level counts, using a
    pre-built tx2gene mapping (see reference_manager.py's
    extract_tx2gene_from_ensembl_fasta / extract_tx2gene_from_gtf /
    build_identity_tx2gene + save_tx2gene_csv).

    sample_quant_paths: dict {sample_name: path_to_quant.sf}
    tx2gene_path: path to a CSV with columns TXNAME, GENEID
    output_path: where to write the resulting gene-level counts matrix
    work_dir: scratch directory for the temporary R script and sample
        manifest CSV this function creates

    Returns (success: bool, message_or_log: str).
    """
    import pandas as pd

    if not tximport_available():
        return False, (
            "Rscript was not found on this system. R with the tximport "
            "package needs to be installed in your environment (it's "
            "included in the project's Dockerfile) before gene-level "
            "collapsing can run."
        )

    os.makedirs(work_dir, exist_ok=True)

    # Write the sample manifest tximport's R script expects.
    manifest_path = os.path.join(work_dir, "tximport_samples.csv")
    manifest_df = pd.DataFrame({
        "sample": list(sample_quant_paths.keys()),
        "path": list(sample_quant_paths.values()),
    })
    manifest_df.to_csv(manifest_path, index=False)

    # Write the R script itself to a temp file rather than passing it
    # inline via `Rscript -e`, since the script is long enough that
    # shell-escaping it inline would be fragile/error-prone.
    r_script_path = os.path.join(work_dir, "run_tximport.R")
    with open(r_script_path, "w") as f:
        f.write(_TXIMPORT_R_SCRIPT)

    cmd = ["Rscript", r_script_path, tx2gene_path, manifest_path, output_path]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=1800
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, f"tximport failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, "tximport timed out after 30 minutes."
