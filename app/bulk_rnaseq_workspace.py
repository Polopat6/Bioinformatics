"""
bulk_rnaseq_workspace.py

Bulk RNA-Seq workspace: guided FASTQ ingestion + metadata matching.

Scope (intentionally limited for now): get raw sequencing files uploaded,
validated, and correctly matched to sample metadata. No alignment, QC
tools, or DESeq2 yet — that comes after this foundation is solid.

Design goal: every step assumes the user may have little to no
bioinformatics background. Each section explains *what* to do and *why*,
using plain language and expandable "help" boxes rather than jargon-only
labels.

This module is fully self-contained. All Bulk RNA-Seq development should
happen here — editing this file has zero effect on the Spatial
Transcriptomics workspace (spatial_workspace.py).
"""

import os
import re
import shutil
import subprocess
import zipfile

import pandas as pd
import streamlit as st

import project_manager as pm
import sra_manager as sra
import file_browser as fb

# Plain-language explanations for each FastQC module, shown when a sample
# gets a WARN or FAIL for that module. Written for someone with little to
# no bioinformatics background — explains what it means and what (if
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

# Extensions the server-side directory browser previews/matches against
# when looking for FASTQ files -- includes the leading "." since
# file_browser.py's extension filtering matches directly against
# os.path.splitext-style suffixes (unlike st.file_uploader's `type`
# argument, which expects no leading dot).
FASTQ_BROWSE_EXTENSIONS = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _build_example_metadata_xlsx():
    """
    Build an example metadata template as an in-memory .xlsx file, for
    users who prefer working in Excel over plain CSV.
    """
    import io

    example_df = pd.DataFrame({
        "sample": ["PatientA", "PatientB"],
        "condition": ["treated", "control"],
    })
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        example_df.to_excel(writer, index=False, sheet_name="Sheet1")
    return buffer.getvalue()


def _read_metadata_file(uploaded_file):
    """
    Read an uploaded metadata file into a DataFrame, supporting both CSV
    (.csv/.txt) and Excel (.xlsx/.xls) formats since users without a
    bioinformatics background often work in Excel rather than plain text.

    Returns (dataframe_or_none, error_message_or_none). If reading fails,
    a plain-language error message is returned instead of raising, so the
    caller can show it directly in the UI.
    """
    filename = uploaded_file.name.lower()
    try:
        if filename.endswith((".xlsx", ".xls")):
            # Reads the first sheet by default, which matches the guidance
            # given to users in the help box above.
            return pd.read_excel(uploaded_file, sheet_name=0), None
        else:
            return pd.read_csv(uploaded_file), None
    except Exception as e:
        return None, (
            "⚠️ We couldn't read this file. Please double check that it's "
            "a valid .csv or .xlsx file with your sample data on the first "
            f"sheet/tab. (Technical detail: {e})"
        )


def _validate_sample_pairs(filenames):
    """
    Group FASTQ filenames by sample name and detect R1/R2 mate pairs.
    Supports multiple common naming conventions:
      - Simple naming:      sample_R1.fastq.gz / sample_R2.fastq.gz
      - Illumina-style:     sample_S1_L001_R1_001.fastq.gz
      - SRA/EBI-style:      sample_1.fastq.gz / sample_2.fastq.gz

    The "_1"/"_2" convention is only matched at the very end of the
    filename (before the extension) to avoid false positives from sample
    names that happen to contain a number, e.g. "Donor10.fastq.gz" is
    correctly treated as single-end, not misread as a "_1" mate pair.

    Takes a plain list of filename strings (not Streamlit UploadedFile
    objects) so it can be used both for freshly uploaded files and for
    files already sitting on disk from a previous session (including
    files symlinked in via the server-side directory browser -- see
    _symlink_fastq_files_from_directory).

    Returns: { sample_name: {"R1": filename, "R2": filename} }
             or { sample_name: {"SE": filename} } for single-end samples.
    """
    sample_pairs = {}
    for name in filenames:
        base = name.replace(".fastq.gz", "").replace(".fastq", "")

        if "_R1" in base:
            sample_name = base.split("_R1")[0]
            sample_pairs.setdefault(sample_name, {})["R1"] = name
        elif "_R2" in base:
            sample_name = base.split("_R2")[0]
            sample_pairs.setdefault(sample_name, {})["R2"] = name
        elif re.search(r"_1$", base):
            sample_name = base[: -len("_1")]
            sample_pairs.setdefault(sample_name, {})["R1"] = name
        elif re.search(r"_2$", base):
            sample_name = base[: -len("_2")]
            sample_pairs.setdefault(sample_name, {})["R2"] = name
        else:
            sample_pairs.setdefault(base, {})["SE"] = name

    return sample_pairs


def _list_existing_fastq(fastq_dir):
    """
    List FASTQ filenames already saved to disk for this project from a
    previous session. Used so reopening a project shows previously
    uploaded files instead of appearing empty.

    Uses os.listdir + a plain name filter (not os.path.isfile), which
    means this INTENTIONALLY also picks up symlinks to real FASTQ files
    elsewhere on disk (e.g. ones created by the server-side "browse for
    a directory of FASTQ files" option below) exactly the same as a
    directly-uploaded, physically-copied file -- from this function's
    (and every downstream consumer's) point of view, a valid symlink and
    a real file are indistinguishable and equally usable.
    """
    if not os.path.isdir(fastq_dir):
        return []
    return sorted([
        name for name in os.listdir(fastq_dir)
        if name.endswith(".fastq") or name.endswith(".fastq.gz") or name.endswith(".gz")
    ])


def _save_uploaded_files(uploaded_files, fastq_dir):
    """
    Persist uploaded FASTQ files to disk inside the active project's
    FASTQ folder. Streamlit's file_uploader only keeps files in memory;
    downstream tools (FastQC, and eventually Nextflow) need real files
    on disk.
    """
    os.makedirs(fastq_dir, exist_ok=True)
    saved_paths = []
    for f in uploaded_files:
        dest_path = os.path.join(fastq_dir, f.name)
        with open(dest_path, "wb") as out_file:
            out_file.write(f.getbuffer())
        saved_paths.append(dest_path)
    return saved_paths


def _symlink_fastq_files_from_directory(source_dir, fastq_dir):
    """
    Scan source_dir (a server-side directory the user confirmed via
    file_browser.render_server_directory_browser -- e.g. a folder that
    already contains a full set of raw FASTQ files for many samples)
    for FASTQ files, and SYMLINK (never copy) each one directly into
    fastq_dir under its original filename.

    Symlinking (rather than copying) is the whole point of this feature:
    raw FASTQ files are often large and numerous (dozens of samples x
    paired-end files can easily total many tens or hundreds of GB), so
    physically duplicating them into fastq_dir would be slow and waste
    a large amount of disk space for no benefit when the files already
    exist in a perfectly usable location on the very same machine
    Streamlit is running on. A symlink is created essentially instantly
    regardless of the target file's size, and every downstream piece of
    code in this app (FastQC, trimming, alignment) reads through the
    symlink exactly as if it were a real file -- no other code needed
    to change to support this.

    Skips (rather than erroring on) any filename that would collide
    with a file/symlink already present in fastq_dir -- e.g. from an
    earlier upload or a previous browse-and-symlink action for the same
    project -- since silently overwriting a different, already-in-use
    file under the same name could otherwise corrupt an existing
    project's data unexpectedly.

    Returns (n_linked: int, n_skipped: int, skipped_names: list[str]).
    """
    os.makedirs(fastq_dir, exist_ok=True)
    _, matching_files = _find_fastq_filenames_in_directory(source_dir)

    n_linked = 0
    skipped_names = []

    for filename in matching_files:
        source_path = os.path.join(source_dir, filename)
        dest_path = os.path.join(fastq_dir, filename)

        if os.path.exists(dest_path) or os.path.islink(dest_path):
            skipped_names.append(filename)
            continue

        os.symlink(source_path, dest_path)
        n_linked += 1

    return n_linked, len(skipped_names), skipped_names


def _find_fastq_filenames_in_directory(source_dir):
    """
    List just the FASTQ-looking filenames directly inside source_dir
    (non-recursive -- files in subdirectories are NOT included, matching
    the same "immediate contents only" convention file_browser.py's own
    directory listing uses, so what the user previewed while browsing
    is exactly what gets matched here).

    Returns (all_entries: list[str], fastq_filenames: list[str]) -- the
    full raw listing (for diagnostic/empty-directory messaging) and the
    subset that looks like a FASTQ file by extension.
    """
    try:
        all_entries = sorted(os.listdir(source_dir))
    except OSError:
        return [], []

    fastq_filenames = [
        name for name in all_entries
        if not name.startswith(".") and os.path.isfile(os.path.join(source_dir, name))
        and any(name.lower().endswith(ext) for ext in FASTQ_BROWSE_EXTENSIONS)
    ]
    return all_entries, fastq_filenames


def _build_match_table(sample_pairs, meta_df):
    """
    Build a plain-language matching table showing, for every sample
    detected from FASTQ filenames AND every sample listed in the metadata
    file, whether they successfully matched up.

    This is the core beginner-friendly output: instead of a cryptic
    error, the user sees exactly which samples matched, which FASTQ
    samples had no metadata row, and which metadata rows had no FASTQ
    files.
    """
    fastq_samples = set(sample_pairs.keys())
    meta_samples = set(meta_df["sample"].astype(str)) if meta_df is not None and "sample" in meta_df.columns else set()

    all_samples = sorted(fastq_samples | meta_samples)
    rows = []
    for sample_name in all_samples:
        has_fastq = sample_name in fastq_samples
        has_meta = sample_name in meta_samples

        if has_fastq and has_meta:
            status = "✅ Matched"
        elif has_fastq and not has_meta:
            status = "⚠️ FASTQ uploaded, no metadata row found"
        else:
            status = "⚠️ Metadata row found, no FASTQ uploaded"

        reads = sample_pairs.get(sample_name, {})
        if "R1" in reads and "R2" in reads:
            read_type = "Paired-end (R1 + R2)"
        elif "SE" in reads:
            read_type = "Single-end"
        elif "R1" in reads:
            read_type = "⚠️ Missing R2 mate"
        elif "R2" in reads:
            read_type = "⚠️ Missing R1 mate"
        else:
            read_type = "—"

        rows.append({
            "Sample": sample_name,
            "Status": status,
            "Read Type": read_type,
        })

    return pd.DataFrame(rows)


def _write_matched_samplesheet(sample_pairs, meta_df, fastq_dir, samplesheet_path):
    """
    Write out only the samples that fully matched (have both FASTQ files
    on disk AND a metadata row) to a clean samplesheet CSV inside the
    active project. This is the file future steps (QC, trimming,
    alignment) will read.
    """
    fastq_samples = set(sample_pairs.keys())
    meta_samples = set(meta_df["sample"].astype(str)) if "sample" in meta_df.columns else set()
    matched_samples = fastq_samples & meta_samples

    rows = []
    for sample_name in sorted(matched_samples):
        reads = sample_pairs[sample_name]
        row = {
            "sample": sample_name,
            "fastq_1": os.path.join(fastq_dir, reads.get("R1", reads.get("SE", ""))),
            "fastq_2": os.path.join(fastq_dir, reads["R2"]) if "R2" in reads else "",
        }
        meta_row = meta_df[meta_df["sample"].astype(str) == sample_name]
        for col in meta_df.columns:
            if col != "sample":
                row[col] = meta_row.iloc[0][col]
        rows.append(row)

    samplesheet_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(samplesheet_path), exist_ok=True)
    samplesheet_df.to_csv(samplesheet_path, index=False)
    return samplesheet_df


def _tools_available():
    """Check whether fastqc and multiqc are installed and on PATH."""
    return shutil.which("fastqc") is not None, shutil.which("multiqc") is not None


def _run_fastqc(fastq_paths, fastqc_dir, threads=4):
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
            # Timeout increased from the original 30 minutes: that limit
            # was set while only testing against small bacterial (E.
            # coli) FASTQ files. Real-world datasets with many
            # samples/large human-genome-scale reads (e.g. 16 paired-end
            # airway samples = 32 files) can legitimately need well
            # beyond 30 minutes even with threading, so this is bumped
            # to 2 hours as a safer ceiling that still catches genuinely
            # stuck/hung processes rather than merely slow ones.
            cmd, capture_output=True, text=True, check=True, timeout=7200
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, "FastQC timed out after 2 hours."


def _run_multiqc(fastqc_dir, multiqc_dir):
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


def _parse_fastqc_summaries(fastqc_dir):
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


def _build_quality_flags(summary_df):
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


def _render_server_directory_fastq_section(fastq_dir):
    """
    Render the "point to a directory of FASTQ files already on this
    server" section of Step 1 -- an alternative to both manual upload
    and SRA download, for the case where this app is running on a
    remote host (HPC node, shared lab server) that already has a
    directory full of raw FASTQ files sitting on its own disk.

    Files are SYMLINKED (never copied) into fastq_dir -- see
    _symlink_fastq_files_from_directory's docstring for why this
    matters for potentially very large, numerous raw FASTQ files.
    """
    with st.expander("📂 Or point to a directory of FASTQ files already on this server"):
        st.markdown(
            "If this app is running on a shared server or HPC that "
            "already has your raw FASTQ files sitting on its own disk "
            "(rather than on your local computer), you can point "
            "directly at the folder containing them here -- this "
            "avoids needing to transfer potentially very large files "
            "through your browser at all. The files are linked into "
            "this project directly from their existing location "
            "(not copied), so this works instantly regardless of how "
            "large or numerous they are."
        )

        selected_dir = fb.render_server_directory_browser(
            key_prefix="fastq_dir_browse",
            preview_extensions=FASTQ_BROWSE_EXTENSIONS,
            label="Browse for the directory containing your FASTQ files:",
        )

        if selected_dir:
            all_entries, fastq_filenames = _find_fastq_filenames_in_directory(selected_dir)

            if not fastq_filenames:
                st.warning(
                    f"⚠️ No FASTQ files (`.fastq`, `.fastq.gz`, `.fq`, "
                    f"`.fq.gz`) were found directly inside this "
                    f"directory ({len(all_entries)} other item(s) "
                    "present). Double-check you've browsed to the right "
                    "folder -- note that files in SUBdirectories of the "
                    "selected folder are not automatically included; "
                    "browse directly to the folder containing the "
                    "actual FASTQ files."
                )
            else:
                st.success(f"✅ Found {len(fastq_filenames)} FASTQ file(s) in this directory.")
                st.dataframe(
                    pd.DataFrame({"File": fastq_filenames}),
                    use_container_width=True, hide_index=True,
                )

                if st.button("🔗 Link These Files Into This Project", key="fastq_dir_browse_link_btn"):
                    n_linked, n_skipped, skipped_names = _symlink_fastq_files_from_directory(
                        selected_dir, fastq_dir
                    )
                    if n_linked:
                        st.success(f"✅ Linked {n_linked} new file(s) into this project.")
                    if n_skipped:
                        st.info(
                            f"ℹ️ Skipped {n_skipped} file(s) that already exist "
                            f"in this project (same filename already present): "
                            f"{', '.join(skipped_names)}."
                        )
                    if n_linked:
                        st.rerun()


def _render_sra_lookup_section(fastq_dir):
    """
    Render the "fetch from SRA/NCBI" section of Step 1. Lets a user look
    up any NCBI accession (run, study, or BioProject), see each run's
    metadata — critically including its official Library Strategy — and
    download selected runs directly into the project's fastq_dir.

    The Library Strategy check exists specifically to prevent a common
    mistake: downloading a public dataset that looks like RNA-seq FASTQ
    files but is actually a different assay type (e.g. whole-genome
    sequencing), which silently produces a low, confusing mapping rate
    much later in the pipeline instead of a clear warning up front.

    Lookup results (including each sample's NCBI-provided attributes)
    are cached in st.session_state["_sra_lookup_rows"] so that Step 2's
    "auto-fill metadata from SRA" feature can reuse them without a
    second lookup.

    Downloads of multiple selected runs are done in parallel (via
    sra_manager.download_sra_runs_parallel) rather than one at a time,
    since prefetch/fasterq-dump are I/O-bound (waiting on network
    transfer) — running several concurrently meaningfully speeds up
    multi-sample downloads. See the note on Streamlit thread-safety
    below for why the progress bar is updated from the main thread only.
    """
    with st.expander("🔎 Or fetch a public dataset from SRA/NCBI"):
        st.markdown(
            "Instead of uploading your own files, you can search for and "
            "download publicly available sequencing data directly from "
            "NCBI's Sequence Read Archive (SRA), using any valid "
            "accession:\n"
            "- A single **run** (e.g. `SRR12345678`)\n"
            "- A **study** (e.g. `SRP123456`)\n"
            "- A whole **BioProject** (e.g. `PRJNA123456`)\n\n"
            "⚠️ **Important:** not everything on public repositories is "
            "RNA-seq data, even if it looks like ordinary FASTQ files. "
            "We check each run's official \"Library Strategy\" annotation "
            "from NCBI and clearly flag anything that isn't RNA-Seq (e.g. "
            "whole-genome sequencing, ChIP-seq), since using the wrong "
            "assay type will produce a low, confusing mapping rate much "
            "later in the pipeline instead of a clear warning now."
        )

        prefetch_ok, fasterq_ok = sra.tools_available()
        if not (prefetch_ok and fasterq_ok):
            missing = []
            if not prefetch_ok:
                missing.append("prefetch")
            if not fasterq_ok:
                missing.append("fasterq-dump")
            st.error(
                f"⚠️ {', '.join(missing)} not found on this system. The "
                "SRA Toolkit needs to be installed in your environment "
                "(it's included in the project's Dockerfile) before this "
                "feature can be used."
            )
            return

        accession_input = st.text_input(
            "Enter an NCBI/SRA accession:",
            placeholder="e.g. SRR12345678, SRP123456, or PRJNA123456",
            key="sra_accession_input",
        )

        if st.button("🔍 Look Up", key="sra_lookup_btn") and accession_input:
            with st.spinner(f"Looking up '{accession_input}' on NCBI..."):
                success, rows, message = sra.lookup_accession(accession_input)

            if not success:
                st.error(f"⚠️ {message}")
            else:
                st.session_state["_sra_lookup_rows"] = rows
                st.success(f"✅ {message}")

        lookup_rows = st.session_state.get("_sra_lookup_rows")
        if lookup_rows:
            display_rows = []
            for row in lookup_rows:
                strategy = row.get("LibraryStrategy", "unknown")
                is_rna = sra.is_rna_seq(row)
                display_rows.append({
                    "Run": row.get("Run", "—"),
                    "Organism": row.get("ScientificName", "—"),
                    "Library Strategy": strategy if is_rna else f"⚠️ {strategy} (NOT RNA-Seq)",
                    "Layout": row.get("LibraryLayout", "—"),
                    "Reads": row.get("spots", "—"),
                    "Size (MB)": row.get("size_MB", "—"),
                })
            results_df = pd.DataFrame(display_rows)
            st.dataframe(results_df, use_container_width=True, hide_index=True)

            non_rna_count = sum(1 for r in lookup_rows if not sra.is_rna_seq(r))
            if non_rna_count > 0:
                st.warning(
                    f"⚠️ {non_rna_count} of {len(lookup_rows)} run(s) above "
                    "are **not** annotated as RNA-Seq. Downloading and "
                    "running these through the RNA-seq pipeline will "
                    "produce a low, misleading mapping rate later — they "
                    "are shown for transparency but are not recommended "
                    "for this workflow unless you specifically intend to "
                    "use them."
                )

            run_options = [row["Run"] for row in lookup_rows]
            selected_runs = st.multiselect(
                "Select which run(s) to download:",
                options=run_options,
                key="sra_selected_runs",
            )

            # Auto-detected default (used below for both concurrency
            # controls), same pattern used for FastQC/fastp/Salmon/STAR
            # elsewhere in the app. detected_cores is the machine's REAL
            # core count; recommended_threads is a conservative SUGGESTED
            # starting value only -- every slider below uses
            # detected_cores as its max_value (so a user on a large
            # machine, e.g. an HPC node with dozens of cores, can
            # actually select a high value) and recommended_threads only
            # to pre-fill each slider's starting point. A previous
            # version of this file capped BOTH sliders' max_value at a
            # low fixed number regardless of the actual machine's core
            # count -- confirmed via real testing on a 32-core HPC node,
            # where these sliders were stuck well below the machine's
            # real capacity.
            detected_cores, recommended_threads = pm.get_recommended_thread_count()

            max_workers_help = (
                "How many samples to download at the same time. Since "
                "downloading is mostly spent waiting on network transfer "
                "rather than heavy local computation, downloading several "
                "samples at once is meaningfully faster than one at a "
                "time. A moderate default is suggested as a safe starting "
                "point -- increase it if your internet connection is "
                "fast and you're downloading many samples; decrease it "
                "if downloads start failing."
            )
            max_workers = st.slider(
                "Parallel downloads:",
                min_value=1, max_value=detected_cores, value=min(3, detected_cores),
                help=max_workers_help,
                key="sra_max_workers",
            )

            threads_per_run = st.slider(
                "Threads per download (fasterq-dump):",
                min_value=1, max_value=detected_cores, value=recommended_threads,
                help=(
                    f"How many threads fasterq-dump uses while extracting "
                    f"*each* sample. Detected {detected_cores} CPU core(s) "
                    f"on this machine, so {recommended_threads} is "
                    "suggested as a starting point. This multiplies with "
                    "'Parallel downloads' above — e.g. 3 parallel "
                    "downloads x 4 threads each means up to 12 threads "
                    "active at once, so consider lowering one slider if "
                    "you raise the other, especially on a smaller "
                    "machine. Both sliders can be raised up to this "
                    f"machine's full {detected_cores}-core capacity if "
                    "you have the resources to support it."
                ),
                key="sra_threads_per_run",
            )

            if selected_runs and st.button("⬇️ Download Selected Run(s)", key="sra_download_btn"):
                n_runs = len(selected_runs)
                progress_bar = st.progress(0, text=f"Starting {n_runs} download(s) ({max_workers} at a time)...")
                status_area = st.empty()
                completed_lines = []

                # IMPORTANT: sra.download_sra_runs_parallel() runs each
                # download in a background worker thread, but Streamlit's
                # UI elements (st.*) are NOT thread-safe to update from
                # those worker threads directly. To stay safe, the actual
                # network calls happen inside download_sra_runs_parallel,
                # but all UI updates here happen through the
                # on_run_complete callback, which download_sra_runs_parallel
                # guarantees is invoked from the main thread (via
                # as_completed(), which yields control back to the caller's
                # thread for each completed future) — so it's safe to call
                # progress_bar.progress(...) and status_area.write(...)
                # from inside this callback.
                completed_count = {"n": 0}

                def _on_complete(accession, success, message):
                    completed_count["n"] += 1
                    icon = "✅" if success else "⚠️"
                    completed_lines.append(f"{icon} **{accession}**: {message}")
                    status_area.markdown("\n\n".join(completed_lines))
                    progress_bar.progress(
                        completed_count["n"] / n_runs,
                        text=f"Completed {completed_count['n']} of {n_runs}...",
                    )

                results = sra.download_sra_runs_parallel(
                    selected_runs, fastq_dir, max_workers=max_workers,
                    on_run_complete=_on_complete, threads_per_run=threads_per_run,
                )

                progress_bar.progress(1.0, text="Done.")

                n_failed = sum(1 for _, success, _ in results if not success)
                if n_failed == 0:
                    st.success(f"✅ All {n_runs} download(s) completed successfully.")
                else:
                    st.warning(f"⚠️ {n_failed} of {n_runs} download(s) failed — see details above.")

                st.info("Refresh below to see the newly downloaded files reflected in your project.")
                st.rerun()


def _render_sra_metadata_autofill_section(project, metadata_saved_path):
    """
    Render the "auto-fill metadata from SRA" section of Step 2. Reuses
    whatever lookup results are cached from Step 1's SRA search
    (st.session_state["_sra_lookup_rows"]) and lets the user build a
    metadata table directly from each sample's NCBI-provided
    SAMPLE_ATTRIBUTES (e.g. strain, media, treatment, genotype) instead
    of typing it in by hand.

    Only shown if there are cached SRA lookup results available — if the
    user hasn't used the SRA lookup feature in Step 1, this section is
    skipped entirely (manual upload remains the only path).

    Returns True if metadata was just auto-filled and saved this run (so
    the caller can skip re-rendering the "no metadata yet" message), or
    False otherwise.
    """
    lookup_rows = st.session_state.get("_sra_lookup_rows")
    if not lookup_rows:
        return False

    with st.expander("📋 Or auto-fill metadata from your SRA lookup"):
        st.markdown(
            "Since you looked up sample(s) from SRA/NCBI in Step 1, we "
            "can try to build a metadata table automatically from each "
            "sample's own NCBI record — using whatever characteristics "
            "fields the original depositor provided (e.g. `strain`, "
            "`media`, `treatment`, `genotype`). **Not every dataset "
            "includes the same fields**, so please review the table "
            "below carefully and edit anything that looks wrong or "
            "incomplete before saving."
        )

        run_options = [row["Run"] for row in lookup_rows]
        selected_for_metadata = st.multiselect(
            "Which run(s) should be included in the auto-filled metadata?",
            options=run_options,
            default=run_options,
            key="sra_metadata_selected_runs",
        )

        if not selected_for_metadata:
            st.info("Select at least one run above to preview auto-filled metadata.")
            return False

        metadata_rows = sra.build_metadata_dataframe(lookup_rows, selected_runs=selected_for_metadata)
        preview_df = pd.DataFrame(metadata_rows)

        if preview_df.empty or len(preview_df.columns) <= 1:
            st.warning(
                "⚠️ NCBI didn't provide any additional characteristics "
                "fields for these sample(s) beyond the accession itself. "
                "You'll need to add a `condition` (or similar) column "
                "manually — download this as a starting point and edit "
                "it in Excel/Sheets, or use the manual upload option "
                "above instead."
            )

        st.markdown("**Preview (edit directly in the table if needed):**")
        edited_df = st.data_editor(
            preview_df,
            use_container_width=True,
            num_rows="dynamic",
            key="sra_metadata_editor",
        )

        if st.button("💾 Save This as My Metadata", key="save_sra_metadata_btn"):
            if "sample" not in edited_df.columns or edited_df["sample"].isna().any():
                st.error("⚠️ Every row needs a value in the `sample` column. Please fix this and try again.")
                return False

            os.makedirs(os.path.dirname(metadata_saved_path), exist_ok=True)
            edited_df.to_csv(metadata_saved_path, index=False)
            pm.save_sample_column(project, "sample")
            st.success(
                f"✅ Metadata saved for {len(edited_df)} sample(s) from "
                "your SRA lookup. Scroll down to continue to Step 3."
            )
            st.rerun()

    return False


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("🧬 Bulk RNA-Seq: Upload & Sample Matching")
    st.markdown(
        "This workspace helps you upload your raw sequencing files and "
        "connect them to your sample information — the first step before "
        "any analysis can run. **No bioinformatics experience required** — "
        "follow the steps below in order."
    )
    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 0: Project setup
    # -----------------------------------------------------------------
    project = pm.render_project_selector(workspace_key="bulk_rnaseq")

    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    fastq_dir = pm.fastq_dir(project)
    samplesheet_path = pm.samplesheet_path(project)
    fastqc_dir = pm.fastqc_dir(project)
    multiqc_dir = pm.multiqc_dir(project)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 1: Upload FASTQ files
    # -----------------------------------------------------------------
    st.header("Step 1: Upload Your Sequencing Files (FASTQ)")

    with st.expander("ℹ️ What is a FASTQ file? (click to learn more)"):
        st.markdown(
            "A **FASTQ file** is the raw output from a DNA/RNA sequencing "
            "machine. It contains the actual sequence of each read, along "
            "with a quality score for every letter. Files are usually "
            "compressed and end in `.fastq.gz` (or sometimes just "
            "`.fastq`).\n\n"
            "**Paired-end vs. single-end:** many experiments sequence each "
            "sample from both ends, producing two files per sample. This "
            "tool recognizes both common naming styles:\n"
            "- `_R1` / `_R2` (most common)\n"
            "- `_1` / `_2` (common for files downloaded from public "
            "repositories like SRA/ENA)\n\n"
            "If your files don't have either suffix, they're likely "
            "**single-end** (one file per sample), which is also "
            "supported.\n\n"
            "**Example filenames:**\n"
            "- `PatientA_R1.fastq.gz` + `PatientA_R2.fastq.gz` (paired-end)\n"
            "- `SRR12345678_1.fastq.gz` + `SRR12345678_2.fastq.gz` (paired-end)\n"
            "- `PatientB.fastq.gz` (single-end)"
        )

    # Show anything already uploaded to this project in a previous session,
    # so reopening a project doesn't look empty even though the files are
    # sitting on disk.
    existing_fastq_names = _list_existing_fastq(fastq_dir)
    if existing_fastq_names:
        st.markdown(f"**{len(existing_fastq_names)} file(s) already in this project:**")
        st.dataframe(
            pd.DataFrame({"File": existing_fastq_names}),
            use_container_width=True, hide_index=True,
        )

    uploaded_fastq = st.file_uploader(
        "Upload one or more FASTQ files (.fastq or .fastq.gz):",
        type=["fastq", "gz"],
        accept_multiple_files=True,
        help="You can select multiple files at once. Hold Ctrl (or Cmd on Mac) to select several files in the upload dialog."
    )

    if uploaded_fastq:
        saved_paths = _save_uploaded_files(uploaded_fastq, fastq_dir)
        st.success(f"✅ {len(saved_paths)} new file(s) uploaded successfully.")

    # --- Point to a directory of FASTQ files already on this server ---
    # (an alternative to browser upload -- see this function's docstring
    # for why this matters when the app runs on a remote host that
    # already has the data on its own disk.)
    _render_server_directory_fastq_section(fastq_dir)

    # --- Fetch from SRA/NCBI (alternative to manual upload) ---
    _render_sra_lookup_section(fastq_dir)

    # Build sample_pairs from the union of files already on disk (from
    # this session's upload above, plus anything from previous sessions,
    # plus anything just downloaded from SRA or symlinked in from a
    # server-side directory) rather than only from this session's upload
    # widget. This is what makes reopening a project actually reflect
    # prior progress.
    all_fastq_names = _list_existing_fastq(fastq_dir)

    sample_pairs = {}
    if all_fastq_names:
        sample_pairs = _validate_sample_pairs(all_fastq_names)

        st.markdown("**What we detected from your file names:**")
        detected_rows = []
        for sample_name, reads in sample_pairs.items():
            if "R1" in reads and "R2" in reads:
                detected_rows.append({"Sample": sample_name, "Type": "Paired-end", "Status": "✅ Complete"})
            elif "SE" in reads:
                detected_rows.append({"Sample": sample_name, "Type": "Single-end", "Status": "✅ Complete"})
            elif "R1" in reads:
                detected_rows.append({"Sample": sample_name, "Type": "Paired-end", "Status": "⚠️ Missing R2 file"})
            elif "R2" in reads:
                detected_rows.append({"Sample": sample_name, "Type": "Paired-end", "Status": "⚠️ Missing R1 file"})
        st.dataframe(pd.DataFrame(detected_rows), use_container_width=True, hide_index=True)

        incomplete = [s for s, r in sample_pairs.items() if not (("R1" in r and "R2" in r) or "SE" in r)]
        if incomplete:
            st.warning(
                f"⚠️ These samples are missing a mate pair file: **{', '.join(incomplete)}**. "
                "If you meant to upload paired-end data, make sure both the "
                "`_R1` and `_R2` files for each sample are uploaded together."
            )
    else:
        st.info("No files uploaded yet. Use the box above to get started.")

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 2: Upload metadata
    # -----------------------------------------------------------------
    st.header("Step 2: Upload Your Sample Information (Metadata)")

    with st.expander("ℹ️ What is a metadata file, and how should I format it? (click to learn more)"):
        st.markdown(
            "Your **metadata file** describes each sample — for example, "
            "which experimental group it belongs to (control vs. treated, "
            "healthy vs. disease, etc.). This is what allows the tool to "
            "later compare groups and find differences in gene expression.\n\n"
            "It should be a simple spreadsheet with **one row per sample**, "
            "saved as either a `.csv` file **or an Excel `.xlsx` file** — "
            "whichever is easier for you. **You don't need to rename any "
            "columns** — after you upload it, you'll be able to pick which "
            "column contains your sample names from a dropdown.\n\n"
            "The values in your sample name column must **match the "
            "sample names from your FASTQ file names** in Step 1 (without "
            "`_R1`/`_R2`, `_1`/`_2`, or file extensions).\n\n"
            "**Example:**\n\n"
            "| sample | condition |\n"
            "|---|---|\n"
            "| PatientA | treated |\n"
            "| PatientB | control |\n\n"
            "You can add extra columns beyond `condition` (e.g. `batch`, "
            "`sex`, `age`) — they'll be carried along automatically.\n\n"
            "**A couple of Excel tips:** make sure your data is on the "
            "**first sheet** of the workbook, and that the first row "
            "contains your column headers."
        )

        template_col1, template_col2 = st.columns(2)
        with template_col1:
            example_csv = "sample,condition\nPatientA,treated\nPatientB,control\n"
            st.download_button(
                "📄 Download template (.csv)",
                data=example_csv,
                file_name="example_metadata_template.csv",
                mime="text/csv",
            )
        with template_col2:
            example_xlsx_bytes = _build_example_metadata_xlsx()
            st.download_button(
                "📊 Download template (.xlsx)",
                data=example_xlsx_bytes,
                file_name="example_metadata_template.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )

    metadata_saved_path = pm.metadata_path(project)

    # --- Auto-fill metadata from SRA (alternative to manual upload) ---
    # Only rendered if Step 1's SRA lookup was actually used this session.
    _render_sra_metadata_autofill_section(project, metadata_saved_path)

    uploaded_meta = st.file_uploader(
        "Upload your sample metadata file (.csv or .xlsx):",
        type=["csv", "txt", "xlsx", "xls"],
    )

    raw_meta_df = None
    read_error = None

    if uploaded_meta:
        # A fresh upload this session always takes priority.
        raw_meta_df, read_error = _read_metadata_file(uploaded_meta)
    elif os.path.exists(metadata_saved_path):
        # No new upload this session, but this project already has a
        # metadata file saved from a previous session (either a manual
        # upload or an SRA auto-fill) — reload it so the project doesn't
        # appear empty.
        raw_meta_df = pd.read_csv(metadata_saved_path)
        st.info(f"📋 Using previously saved metadata for this project ({len(raw_meta_df)} row(s)). Upload a new file above to replace it.")
    elif os.path.exists(samplesheet_path):
        # Self-healing fallback for projects created before metadata.csv
        # persistence existed: the matched samplesheet already contains
        # the sample name + original metadata columns merged together, so
        # we can reconstruct a usable metadata view from it instead of
        # showing the project as if metadata was never uploaded.
        legacy_df = pd.read_csv(samplesheet_path)
        reconstructed_cols = [c for c in legacy_df.columns if c not in ("fastq_1", "fastq_2")]
        raw_meta_df = legacy_df[reconstructed_cols]
        # Persist it going forward so this fallback only needs to run once.
        os.makedirs(os.path.dirname(metadata_saved_path), exist_ok=True)
        raw_meta_df.to_csv(metadata_saved_path, index=False)
        pm.save_sample_column(project, "sample")
        st.info(
            f"📋 Recovered metadata for this project from your previously "
            f"saved sample list ({len(raw_meta_df)} row(s)). Upload a new "
            "file above if you'd like to replace it."
        )

    meta_df = None
    if read_error:
        st.error(read_error)
    elif raw_meta_df is not None:
        if raw_meta_df.empty or len(raw_meta_df.columns) == 0:
            st.error("⚠️ This file appears to be empty. Please check it and re-upload.")
        else:
            if uploaded_meta:
                st.success("✅ Metadata file uploaded! Here's a preview:")
            st.dataframe(raw_meta_df, use_container_width=True, hide_index=True)

            st.markdown("**Which column contains your sample names?**")
            # Prefer a previously saved selection for this project (from an
            # earlier session); otherwise guess based on common column names.
            saved_column = pm.get_sample_column(project)
            column_options = list(raw_meta_df.columns)
            likely_defaults = ["sample", "sample_id", "sampleid", "sample name", "samplename", "id"]

            default_index = 0
            if saved_column in column_options:
                default_index = column_options.index(saved_column)
            else:
                for i, col in enumerate(column_options):
                    if col.strip().lower() in likely_defaults:
                        default_index = i
                        break

            sample_col = st.selectbox(
                "Sample name column:",
                options=column_options,
                index=default_index,
                help="Pick whichever column lists your sample names — it doesn't need to be called 'sample'.",
            )

            # Normalize: create a working copy with the chosen column
            # renamed to "sample", so all downstream matching logic stays
            # simple regardless of what the user's original column was
            # called.
            if sample_col != "sample" and "sample" in raw_meta_df.columns:
                # Avoid an accidental column name collision.
                meta_df = raw_meta_df.rename(columns={"sample": "sample_original", sample_col: "sample"})
            else:
                meta_df = raw_meta_df.rename(columns={sample_col: "sample"})

            n_missing = meta_df["sample"].isna().sum()
            if n_missing > 0:
                st.warning(
                    f"⚠️ {n_missing} row(s) have a blank value in the "
                    f"selected column (`{sample_col}`). Those rows will be "
                    "ignored when matching to your FASTQ files."
                )
                meta_df = meta_df.dropna(subset=["sample"])

            meta_df["sample"] = meta_df["sample"].astype(str).str.strip()

            # Persist the raw (pre-rename) metadata + the chosen column so
            # this project remembers it next time it's reopened.
            os.makedirs(os.path.dirname(metadata_saved_path), exist_ok=True)
            raw_meta_df.to_csv(metadata_saved_path, index=False)
            pm.save_sample_column(project, sample_col)
    else:
        st.info("No metadata file uploaded yet. Use the box above, or download the template if you're not sure how to format one.")

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 3: Match FASTQ files to metadata
    # -----------------------------------------------------------------
    st.header("Step 3: Check That Everything Matches Up")

    if not sample_pairs or meta_df is None:
        st.info(
            "Complete Step 1 (upload FASTQ files) and Step 2 (upload "
            "metadata) to see your matching results here."
        )
        return

    match_table = _build_match_table(sample_pairs, meta_df)
    st.dataframe(match_table, use_container_width=True, hide_index=True)

    n_matched = (match_table["Status"] == "✅ Matched").sum()
    n_total = len(match_table)

    if n_matched == n_total and n_total > 0:
        st.success(
            f"🎉 All {n_total} sample(s) matched successfully! Every FASTQ "
            "sample has a corresponding metadata row, and vice versa."
        )

        samples_already_saved = os.path.exists(samplesheet_path)
        save_label = "🔄 Re-save matched sample list" if samples_already_saved else "💾 Save matched sample list"

        if samples_already_saved and not st.session_state.get("bulk_rnaseq_matched"):
            st.info("✅ A matched sample list has already been saved for this project.")
            st.dataframe(pd.read_csv(samplesheet_path), use_container_width=True, hide_index=True)

        if st.button(save_label):
            samplesheet_df = _write_matched_samplesheet(sample_pairs, meta_df, fastq_dir, samplesheet_path)
            st.session_state["bulk_rnaseq_matched"] = True
            pm.mark_step_complete(project, "samples_matched")
            st.success(f"Saved to `{samplesheet_path}`.")
            st.dataframe(samplesheet_df, use_container_width=True, hide_index=True)

        # -------------------------------------------------------------
        # STEP 4: Run Quality Control (FastQC + MultiQC)
        # -------------------------------------------------------------
        if st.session_state.get("bulk_rnaseq_matched") or os.path.exists(samplesheet_path):
            st.markdown("---")
            st.header("Step 4: Check Read Quality (FastQC + MultiQC)")

            with st.expander("ℹ️ What is quality control, and why does it matter? (click to learn more)"):
                st.markdown(
                    "Before analyzing your data, it's important to check "
                    "that the raw sequencing reads are actually good "
                    "quality. **FastQC** examines each file individually "
                    "and checks things like base quality scores, leftover "
                    "adapter sequences, and unusual duplication levels. "
                    "**MultiQC** then combines every individual FastQC "
                    "report into one single, easy-to-browse report so you "
                    "don't have to open dozens of separate files.\n\n"
                    "Not every warning is a problem — some are completely "
                    "normal for RNA-seq data. Below, each flagged item "
                    "comes with a plain-language explanation of what it "
                    "means and whether you need to do anything about it."
                )

            fastqc_ok, multiqc_ok = _tools_available()
            multiqc_html_path = os.path.join(multiqc_dir, "multiqc_report.html")
            qc_already_done = os.path.exists(multiqc_html_path)

            if not (fastqc_ok and multiqc_ok):
                missing = []
                if not fastqc_ok:
                    missing.append("FastQC")
                if not multiqc_ok:
                    missing.append("MultiQC")
                st.error(
                    f"⚠️ {', '.join(missing)} not found on this system. "
                    "These tools need to be installed in your environment "
                    "(they're included in the project's Dockerfile) before "
                    "this step can run."
                )
            else:
                button_label = "🔄 Re-run Quality Control" if qc_already_done else "🔬 Run Quality Control on Matched Samples"

                # Let the user control how many files FastQC processes at
                # once, same pattern as the SRA parallel-download slider —
                # higher values speed up larger batches (e.g. 16 samples /
                # 32 paired-end files) at the cost of using more CPU/RAM
                # simultaneously.
                #
                # IMPORTANT: detected_cores (the machine's REAL core
                # count) and recommended_threads (a conservative
                # SUGGESTED starting value, deliberately capped low by
                # project_manager.py's get_recommended_thread_count) serve
                # two different purposes here. This slider uses
                # detected_cores for max_value (so a user on a large
                # machine -- e.g. an HPC node with dozens of cores -- can
                # actually select a high thread count) and
                # recommended_threads only for the slider's initial
                # `value` (a sensible, un-intimidating starting point). A
                # previous version of this file incorrectly used
                # max(8, recommended_threads) as max_value, which --
                # since recommended_threads is itself always <= 8 --
                # silently capped this slider at 8 regardless of how many
                # cores were actually detected (confirmed via real
                # testing on a 32-core machine, where this slider was
                # stuck at a max of 8 despite 32 cores being available).
                detected_cores, recommended_threads = pm.get_recommended_thread_count()
                fastqc_threads = st.slider(
                    "FastQC parallel threads:",
                    min_value=1, max_value=detected_cores, value=recommended_threads,
                    help=(
                        f"How many files FastQC processes at the same time. "
                        f"Detected {detected_cores} CPU core(s) on this "
                        f"machine, so {recommended_threads} is suggested as "
                        "a starting point (leaving some headroom for the "
                        "rest of the system). Higher values finish faster "
                        "for larger sample batches but use more CPU and "
                        "memory at once -- you can raise this up to your "
                        "machine's full core count if you have the "
                        "resources to support it; lower it if your "
                        "computer seems to be struggling."
                    ),
                    key="fastqc_threads_slider",
                )

                run_clicked = st.button(button_label)

                if qc_already_done and not run_clicked:
                    st.success("✅ Quality control has already been run for this project.")

                if run_clicked:
                    fastq_paths = []
                    for sample_name in sorted(set(sample_pairs.keys()) & set(meta_df["sample"].astype(str))):
                        reads = sample_pairs[sample_name]
                        for key in ("R1", "R2", "SE"):
                            if key in reads:
                                fastq_paths.append(os.path.join(fastq_dir, reads[key]))

                    with st.spinner(f"Running FastQC on {len(fastq_paths)} file(s) ({fastqc_threads} at a time)..."):
                        fastqc_success, fastqc_log = _run_fastqc(fastq_paths, fastqc_dir, threads=fastqc_threads)

                    if not fastqc_success:
                        st.error("FastQC failed to run. Details below:")
                        st.code(fastqc_log)
                        st.info(
                            "💡 You can click the button above again to "
                            "retry — if it keeps timing out, try lowering "
                            "the thread count (sometimes fewer, more "
                            "reliable threads beat many contended ones on "
                            "memory-constrained machines)."
                        )
                        qc_already_done = False
                    else:
                        st.success("✅ FastQC completed.")
                        with st.spinner("Combining reports with MultiQC..."):
                            multiqc_success, multiqc_log = _run_multiqc(fastqc_dir, multiqc_dir)

                        if not multiqc_success:
                            st.error("MultiQC failed to run. Details below:")
                            st.code(multiqc_log)
                            qc_already_done = False
                        else:
                            st.success("✅ MultiQC report generated.")
                            pm.mark_step_complete(project, "qc_complete")
                            qc_already_done = True

                # Show results (either freshly generated above, or from a
                # previous session) whenever the MultiQC report exists on
                # disk — this is what makes reopening a project actually
                # display prior QC results instead of requiring a re-run.
                if qc_already_done:
                    summary_df = _parse_fastqc_summaries(fastqc_dir)
                    overview_df, details_by_file = _build_quality_flags(summary_df)

                    if not overview_df.empty:
                        st.subheader("📋 Quality Overview")
                        st.dataframe(overview_df, use_container_width=True, hide_index=True)

                        for filename, explanations in details_by_file.items():
                            if explanations:
                                with st.expander(f"Details for {filename}"):
                                    for line in explanations:
                                        st.markdown(f"- {line}")

                    # --- Offer the combined MultiQC report as a download ---
                    # Note: we intentionally do NOT embed the report via
                    # st.components.v1.html. MultiQC reports contain their
                    # own internal navigation/links, which can hijack the
                    # embedded iframe and redirect the whole app. A plain
                    # download button avoids that entirely.
                    if os.path.exists(multiqc_html_path):
                        st.subheader("📊 Full Combined Report")
                        with open(multiqc_html_path, "rb") as f:
                            report_bytes = f.read()
                        st.download_button(
                            "⬇️ Download Full MultiQC Report (.html)",
                            data=report_bytes,
                            file_name="multiqc_report.html",
                            mime="text/html",
                        )
                        st.caption(
                            "Download the file above, then open it in your "
                            "web browser to view the full interactive report "
                            "with detailed charts for every sample."
                        )

                    st.markdown("---")
                    st.success(
                        f"🎉 Project `{project}` now has matched samples and "
                        "quality control results saved. This project is "
                        "ready for the next step: adapter trimming and "
                        "post-trimming quality control."
                    )

                    if st.button("➡️ Proceed to Trimming & Post-Trimming QC", type="primary"):
                        # Jump straight to the Trimming workspace, carrying
                        # over the same active project automatically (the
                        # project selection lives in
                        # st.session_state["bulk_rnaseq_project"], which is
                        # shared across workspace modules).
                        #
                        # Note: we can't set st.session_state["assay_choice_radio"]
                        # directly here, since app.py's sidebar radio widget
                        # (that key) has already been instantiated earlier in
                        # this same run. Instead we set a plain "nav_request"
                        # flag, which app.py applies to assay_choice_radio on
                        # the next run, before the radio widget is created.
                        st.session_state["nav_request"] = "🧪 Trimming & Post-Trim QC"
                        st.rerun()
    else:
        st.warning(
            f"⚠️ {n_matched} of {n_total} sample(s) matched. Review the "
            "table above:\n\n"
            "- **\"FASTQ uploaded, no metadata row found\"** — add a row "
            "for this sample in your metadata file, using the exact same "
            "sample name.\n"
            "- **\"Metadata row found, no FASTQ uploaded\"** — either "
            "upload the missing FASTQ file(s), or remove this row from "
            "your metadata file if it's not part of this analysis.\n\n"
            "Once every sample shows ✅ Matched, you'll be able to save "
            "your matched sample list."
        )
