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
import math
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


def estimate_genome_length(genome_fasta):
    """
    Estimate a genome FASTA's total sequence length (bp), by summing
    the length of every non-header line -- used to auto-compute STAR's
    --genomeSAindexNbases for small genomes (see
    compute_genome_sa_index_nbases below for why this matters).

    This is a simple, dependency-free line-scan rather than a full
    FASTA parse (we only need a total length, not per-sequence
    lengths), so it stays fast even for a multi-GB genome FASTA.
    """
    total_length = 0
    with open(genome_fasta, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                continue
            total_length += len(line.strip())
    return total_length


def compute_genome_sa_index_nbases(genome_length):
    """
    Compute STAR's recommended --genomeSAindexNbases value for a genome
    of the given length, following STAR's own manual formula:

        min(14, log2(GenomeLength) / 2 - 1)

    STAR's default (14) is sized for large genomes (human/mouse-scale,
    several Gb). For smaller genomes -- a single chromosome, a custom
    non-model organism reference, a bacterial genome -- using the
    default value produces an oversized, invalid, or badly-degraded
    suffix array index; STAR's own log output explicitly recommends
    reducing this value for small genomes. Since custom/non-model
    species uploads are a core supported workflow in this app (see
    alignment_workspace.py's custom reference path), and those are
    disproportionately likely to be small genomes, this must be
    computed automatically rather than left at STAR's large-genome
    default in every case.

    A floor of 4 is applied since going much lower rarely makes sense
    and protects against a degenerate (near-zero or negative) result
    for extremely small inputs (e.g. a single short custom contig).
    """
    if genome_length <= 0:
        return 14
    recommended = min(14, int(math.log2(genome_length) / 2 - 1))
    return max(4, recommended)


def detect_fastq_read_length(fastq_path, is_gzipped=True):
    """
    Detect the actual read length from a FASTQ file's first read, for
    auto-computing STAR's --sjdbOverhang (which should be read_length -
    1 per STAR's manual, not left at a fixed guess -- a mismatched
    value measurably hurts splice-junction detection accuracy,
    especially for read lengths that differ from the common
    100-150bp range this app's previous fixed default of 100 assumed).

    is_gzipped: whether fastq_path is gzip-compressed (true for every
        trimmed FASTQ this app produces, per trimming_workspace.py's
        fastp-based trimming step, which always writes .fastq.gz
        output).

    Returns the read length as an int, or None if the file couldn't be
    read/parsed (e.g. empty file) -- callers should fall back to a
    reasonable default (100) in that case rather than crashing.
    """
    try:
        if is_gzipped:
            import gzip
            opener = gzip.open
        else:
            opener = open

        with opener(fastq_path, "rt", errors="ignore") as f:
            f.readline()  # FASTQ record 1: header line (starts with @)
            sequence_line = f.readline().strip()  # FASTQ record 2: the actual sequence
            if sequence_line:
                return len(sequence_line)
        return None
    except (OSError, EOFError):
        return None


# ---------------------------------------------------------------------------
# STAR: file-descriptor-limit-aware BAM sorting bin count
# ---------------------------------------------------------------------------
#
# STAR's coordinate-sorted BAM output (--outSAMtype BAM SortedByCoordinate)
# works by splitting alignments into --outBAMsortingBinsN genome bins
# (default: 50) and sorting each bin's worth of reads independently
# before merging -- this parallelizes the sort across threads. The
# practical consequence: STAR needs roughly one open file handle PER
# BIN, PER THREAD simultaneously during this step. With STAR's default
# 50 bins and enough threads (e.g. 24, as used in real testing on a
# 32-core HPC node), this can reach into the hundreds to low thousands
# of simultaneously open files -- comfortably exceeding the DEFAULT
# per-process open-file limit on many Linux systems/HPC login or
# interactive nodes (commonly 1024, via `ulimit -n`), causing STAR to
# fail with a "could not create output file ... BAMsort" error. This was
# hit directly during real full-genome alignment testing (24 threads,
# default 50 bins, default ulimit -n 1024 -> STAR crashed; confirmed
# fixed by either raising ulimit -n in the shell, OR -- what this
# function does instead -- lowering the bin count so the ACTUAL number
# of simultaneously open files this run will need comfortably stays
# under whatever this machine's limit happens to be, with no shell
# configuration required at all).

def get_open_file_limit():
    """
    Return this process's current soft limit on simultaneously open
    file descriptors (equivalent to the shell's `ulimit -n`), via the
    standard library's resource module (POSIX-only -- this app already
    assumes a Linux/Mac environment for all its other external tools,
    so this introduces no new platform constraint).

    Returns an int, or a conservative fallback of 1024 (a very common
    real-world default -- confirmed via direct testing to be the exact
    value that caused a real STAR failure) if the limit can't be
    determined for any reason (e.g. on a non-POSIX platform where the
    resource module itself isn't available).
    """
    try:
        import resource
        soft_limit, _hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
        return soft_limit
    except (ImportError, ValueError, OSError):
        return 1024


def compute_safe_bam_sorting_bins(threads, file_descriptor_limit=None, default_bins=50,
                                   safety_margin=0.5):
    """
    Compute a safe value for STAR's --outBAMsortingBinsN, given how many
    threads will be used and this machine's actual open-file-descriptor
    limit -- see the module-level comment above for the full rationale.

    threads: the --runThreadN value this alignment run will actually
        use -- more threads means more simultaneously open files for
        any given bin count, so the safe bin count must scale down as
        threads scales up.

    file_descriptor_limit: this machine's open-file limit -- if None
        (the default), detected automatically via get_open_file_limit().
        Exposed as a parameter mainly to make this function easily
        testable without needing to actually change the test process's
        real ulimit.

    default_bins: STAR's own default (50) -- returned unchanged
        whenever the detected file-descriptor limit is generous enough
        to comfortably support it at the given thread count, so runs on
        typically-configured machines see no behavior change at all
        from this safety logic.

    safety_margin: what fraction of the raw file-descriptor limit to
        actually budget for STAR's BAM sorting specifically (default
        50%) -- deliberately conservative, since STAR is not the only
        thing that may need open file descriptors during this run
        (Python's own subprocess machinery, any other files this
        process has open, etc.), and because the true simultaneous
        file count during BAM sorting isn't an exact, guaranteed
        (threads x bins) in every case -- this margin is a buffer
        against undercounting.

    Returns an int bin count: min(default_bins, a safe value computed
    from the file-descriptor budget divided by thread count), with a
    floor of 1 so this never returns a nonsensical zero-or-negative
    value even at extremely low limits/high thread counts.
    """
    if file_descriptor_limit is None:
        file_descriptor_limit = get_open_file_limit()

    safe_bins = int((file_descriptor_limit * safety_margin) / max(1, threads))
    return max(1, min(default_bins, safe_bins))


# ---------------------------------------------------------------------------
# STAR: ENCODE-recommended options bundle
# ---------------------------------------------------------------------------
#
# These are the standard, widely-used "ENCODE options" for bulk RNA-seq
# alignment with STAR -- a well-established bundle (not our own
# invention) that the STAR manual itself references as the settings
# used in the ENCODE project's own RNA-seq pipelines. They're offered
# here as ONE combined toggle rather than exposing each flag
# individually, since a non-expert user has no meaningful way to decide
# on any one of these in isolation -- they're designed and validated as
# a set.
ENCODE_OPTIONS_FLAGS = [
    "--outFilterType", "BySJout",
    "--outFilterMultimapNmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "1",
    "--outFilterMismatchNmax", "999",
    "--outFilterMismatchNoverLmax", "0.04",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--alignMatesGapMax", "1000000",
]


def build_star_index(genome_fasta, gtf_path, index_dir, threads=4, sjdb_overhang=100,
                      genome_sa_index_nbases=None):
    """
    Build a STAR genome index from a genome FASTA + GTF annotation.

    sjdb_overhang should ideally be set to (read_length - 1); 100 is a
    reasonable default for typical 100-150bp reads and rarely matters
    much in practice for genomes with reads in that common range. For
    read lengths that differ meaningfully from that (e.g. this app's
    older/shorter-read test data), the CALLER should compute this via
    detect_fastq_read_length() on an actual trimmed FASTQ first, rather
    than relying on this default -- see alignment_workspace.py's index
    build step.

    genome_sa_index_nbases: if None (the default), this is
    AUTO-COMPUTED from the genome FASTA's actual length via
    estimate_genome_length() + compute_genome_sa_index_nbases() --
    critical for small (e.g. single-chromosome, bacterial, or other
    custom non-model organism) genomes, where STAR's built-in default
    (effectively 14) produces an oversized/invalid index. Pass an
    explicit int only to override this auto-detection.

    Returns (success: bool, log: str).
    """
    os.makedirs(index_dir, exist_ok=True)

    if genome_sa_index_nbases is None:
        genome_length = estimate_genome_length(genome_fasta)
        genome_sa_index_nbases = compute_genome_sa_index_nbases(genome_length)

    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--genomeDir", index_dir,
        "--genomeFastaFiles", genome_fasta,
        "--sjdbGTFfile", gtf_path,
        "--sjdbOverhang", str(sjdb_overhang),
        "--runThreadN", str(threads),
        "--genomeSAindexNbases", str(genome_sa_index_nbases),
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

def run_star_align(sample_entry, index_dir, output_base_dir, threads=4,
                    use_two_pass=False, use_encode_options=False,
                    add_strand_field=False, limit_bam_sort_ram=None,
                    bam_sorting_bins=None):
    """
    Run STAR alignment for a single sample (paired or single-end), using
    --quantMode GeneCounts to get gene-level counts directly from STAR
    without needing a separate counting tool (e.g. featureCounts).

    use_two_pass: adds --twopassMode Basic. STAR performs a first-pass
        alignment to discover novel splice junctions, then re-aligns
        using those junctions as additional annotation -- meaningfully
        improves novel splice junction sensitivity, at the cost of
        roughly doubling alignment time per sample. Most valuable for
        less-annotated/non-model organisms where the supplied GTF is
        likely incomplete; less impactful for well-annotated genomes
        like human/mouse where most real junctions are already known.

    use_encode_options: adds the standard ENCODE-recommended options
        bundle (see ENCODE_OPTIONS_FLAGS above) -- a well-established,
        commonly-used preset for bulk RNA-seq, not exposed as
        individual flags since they're designed/validated as a set.

    add_strand_field: adds --outSAMstrandField intronMotif, which
        annotates each alignment's strand in the output BAM based on
        intron motif -- needed only if the resulting BAM will be fed
        into a downstream tool that expects this strand tag (e.g.
        Cufflinks-family tools); has no effect on the gene counts this
        app itself uses, so most users can safely leave this off.

    limit_bam_sort_ram: optional int (bytes) to explicitly cap the
        memory STAR's BAM-sorting step is allowed to use, via
        --limitBAMsortRAM. Only relevant if a run is failing with a
        BAM-sorting memory error on a very large genome/very deep
        sample -- left as None (STAR's own default behavior) otherwise.

    bam_sorting_bins: controls --outBAMsortingBinsN, the number of
        genome bins STAR splits alignments into for parallelized
        coordinate-sorting during BAM output. If None (the default),
        AUTO-COMPUTED via compute_safe_bam_sorting_bins(threads) based
        on this machine's actual open-file-descriptor limit -- see the
        module-level comment above compute_safe_bam_sorting_bins for
        the full rationale (confirmed to fix a real crash hit during
        real full-genome alignment testing at higher thread counts).
        Pass an explicit int to override this auto-detection.

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

    if bam_sorting_bins is None:
        bam_sorting_bins = compute_safe_bam_sorting_bins(threads)

    cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--genomeDir", index_dir,
        "--readFilesIn", *read_files_arg,
        "--readFilesCommand", "zcat",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outBAMsortingBinsN", str(bam_sorting_bins),
        "--quantMode", "GeneCounts",
        "--outFileNamePrefix", out_prefix,
        "--runThreadN", str(threads),
    ]

    if use_two_pass:
        cmd += ["--twopassMode", "Basic"]
    if use_encode_options:
        cmd += ENCODE_OPTIONS_FLAGS
    if add_strand_field:
        cmd += ["--outSAMstrandField", "intronMotif"]
    if limit_bam_sort_ram:
        cmd += ["--limitBAMsortRAM", str(limit_bam_sort_ram)]

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
# Mapping-rate quality classification (plain-language QC flag)
# ---------------------------------------------------------------------------
#
# Discovered via a real usability test (a non-bioinformatics user ran
# through the full pipeline against a genuinely mismatched reference --
# a different E. coli STRAIN than the one the experiment actually used):
# every sample was shown as a plain "✅ Success" in the results table
# regardless of its actual mapping rate, even though every single sample
# had a mapping rate under 50%. The existing caption text ("rates above
# ~70-80% are typically considered good...") was present on the SAME
# PAGE but easy to miss/ignore, especially for a user with no prior
# reference point for what a "good" rate even looks like. This function
# provides a hard, automatic 3-tier classification so a badly mismatched
# reference (wrong species, wrong strain, contamination, degraded RNA)
# is impossible to overlook, rather than relying on the user to
# separately notice and interpret a raw percentage number themselves.

def classify_mapping_rate(rate_pct):
    """
    Classify a mapping rate (Salmon's overall "Mapping Rate", or STAR's
    "Uniquely Mapped %") into a plain-language quality tier.

    Thresholds:
      - rate_pct >= 80: "good" -- a healthy mapping rate, no action
        needed. This matches the commonly cited "80%+ is good for a
        well-annotated reference" guidance already used elsewhere in
        this app's own help text.
      - 50 <= rate_pct < 80: "caution" -- lower than ideal but not
        necessarily wrong; many legitimate real-world datasets
        (degraded FFPE samples, less-complete non-model organism
        references, etc.) genuinely fall in this range. Flagged so the
        user is prompted to double-check their setup, without being
        told outright that something is broken.
      - rate_pct < 50: "poor" -- a strong, hard-to-ignore signal that
        something is very likely wrong. The single most common real-
        world cause (confirmed directly via usability testing) is a
        MISMATCHED REFERENCE -- not necessarily the wrong species, but
        potentially the wrong STRAIN of the same species (e.g. two
        bacterial strains can differ enough in their genome sequence
        that most reads fail to align to the wrong one's reference).
        Other plausible causes (sample contamination, badly degraded
        RNA, a corrupted/incomplete reference download) are listed too,
        but reference mismatch is called out first since it's both the
        most common cause in practice and the easiest for a user to
        actually go back and check/fix.

    rate_pct may be None (e.g. if the mapping-rate value itself
    couldn't be parsed from the tool's own output) -- handled
    gracefully as its own "unknown" tier rather than raising, since a
    missing rate is a different (and separately already-surfaced)
    problem from a genuinely low one.

    Returns a dict:
        {
            "tier": "good" | "caution" | "poor" | "unknown",
            "icon": a single emoji appropriate for inline display
                    (e.g. directly inside a results table cell),
            "short_label": a brief tier description suitable for a
                    table cell (e.g. "Good", "Check reference"),
            "message": a full plain-language explanation, suitable for
                    display via st.success/st.warning/st.error
                    immediately below/near the results table for any
                    sample that isn't a clean "good".
        }
    """
    if rate_pct is None:
        return {
            "tier": "unknown", "icon": "⚪", "short_label": "Unknown",
            "message": "Mapping rate could not be determined for this sample.",
        }

    if rate_pct >= 80:
        return {
            "tier": "good", "icon": "🟢", "short_label": "Good",
            "message": f"{rate_pct}% -- a healthy mapping rate.",
        }
    elif rate_pct >= 50:
        return {
            "tier": "caution", "icon": "🟡", "short_label": "Check reference",
            "message": (
                f"{rate_pct}% is lower than the ~80%+ typically expected "
                "for a well-matched reference. This CAN still be a "
                "legitimate result (e.g. a partially-annotated non-model "
                "organism, or somewhat degraded samples), but it's worth "
                "double-checking that you selected the correct species "
                "(and, for organisms with multiple common strains/"
                "sub-species, the correct STRAIN) before trusting these "
                "results."
            ),
        }
    else:
        return {
            "tier": "poor", "icon": "🔴", "short_label": "Likely wrong reference",
            "message": (
                f"{rate_pct}% is unusually low and is a strong signal "
                "that something doesn't match. The single most common "
                "cause is a MISMATCHED REFERENCE -- this doesn't always "
                "mean the wrong species outright; it's also very common "
                "to select the right species but the wrong STRAIN (e.g. "
                "two E. coli strains, or two mouse sub-strains, can "
                "differ enough in their genome sequence that most reads "
                "fail to align to a mismatched one). Other possible "
                "causes: sample contamination, badly degraded RNA, or a "
                "corrupted/incomplete reference download. Please "
                "double-check your reference/species selection in Step "
                "2 before proceeding to use these results."
            ),
        }
