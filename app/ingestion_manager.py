"""
ingestion_manager.py

Streamlit-independent FASTQ ingestion + sample/metadata matching logic
for the Bulk RNA-Seq pipeline's first workflow step.

Extracted out of bulk_rnaseq_workspace.py so this logic can be called
both from that interactive workspace AND from a non-interactive
orchestrator (Advanced Mode / Monitor Mode), which need to run the same
"figure out what FASTQ files exist, detect R1/R2 mate pairs, match them
against the metadata sheet, and write out a clean samplesheet.csv"
pipeline in the background without any Streamlit UI calls in the way.

This module owns:
  - reading a metadata file (CSV or Excel) into a DataFrame
  - listing/detecting FASTQ files already on disk for a project
  - persisting uploaded FASTQ files to disk
  - symlinking FASTQ files in from a server-side source directory
    (never copying -- see symlink_fastq_files_from_directory)
  - grouping FASTQ filenames into samples and detecting R1/R2 mate pairs
    across several common naming conventions
  - removing a sample's FASTQ file(s) from a project
  - building the plain-language FASTQ<->metadata match table
  - writing the final matched samplesheet.csv that every downstream
    step (QC, trimming, alignment, quantification) reads from
  - merging/detecting blanks in NCBI-derived metadata attribute rows

It does NOT know about Streamlit, project_manager step-tracking, or
anything UI-related -- callers (bulk_rnaseq_workspace.py today, an
Advanced/Monitor Mode orchestrator in the future) are responsible for
deciding when to call these functions, showing progress/errors, and
marking the project's "samples_matched" step once a samplesheet has
been written successfully.
"""

import os
import re

import pandas as pd

# Extensions the server-side directory browser previews/matches against
# when looking for FASTQ files -- includes the leading "." since
# file_browser.py's extension filtering matches directly against
# os.path.splitext-style suffixes (unlike st.file_uploader's `type`
# argument, which expects no leading dot).
FASTQ_BROWSE_EXTENSIONS = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]


def read_metadata_file(file_or_path):
    """
    Read a metadata file into a DataFrame, supporting both CSV
    (.csv/.txt) and Excel (.xlsx/.xls) formats since users without a
    bioinformatics background often work in Excel rather than plain text.

    file_or_path: either a Streamlit UploadedFile object (from
        st.file_uploader) OR a plain string path to an existing file on
        disk (e.g. one selected via a server-side file browser, or an
        Advanced Mode project's pre-supplied metadata path) -- pandas'
        read_csv/read_excel both accept either a file-like object or a
        path string transparently, so no special-casing is needed here
        beyond checking the filename to decide CSV vs. Excel parsing.

    Returns (dataframe_or_none, error_message_or_none). If reading fails,
    a plain-language error message is returned instead of raising, so the
    caller can show it directly in the UI (or log it, for a
    non-interactive caller).
    """
    filename = (
        file_or_path.name if hasattr(file_or_path, "name") else os.path.basename(file_or_path)
    ).lower()
    try:
        if filename.endswith((".xlsx", ".xls")):
            # Reads the first sheet by default, which matches the guidance
            # given to users in the interactive workspace's help box.
            return pd.read_excel(file_or_path, sheet_name=0), None
        else:
            return pd.read_csv(file_or_path), None
    except Exception as e:
        return None, (
            "⚠️ We couldn't read this file. Please double check that it's "
            "a valid .csv or .xlsx file with your sample data on the first "
            f"sheet/tab. (Technical detail: {e})"
        )


def validate_sample_pairs(filenames):
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
    symlink_fastq_files_from_directory).

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


def list_existing_fastq(fastq_dir):
    """
    List FASTQ filenames already saved to disk for this project from a
    previous session. Used so reopening a project shows previously
    uploaded files instead of appearing empty.

    Uses os.listdir + a plain name filter (not os.path.isfile), which
    means this INTENTIONALLY also picks up symlinks to real FASTQ files
    elsewhere on disk (e.g. ones created by the server-side "browse for
    a directory of FASTQ files" option) exactly the same as a
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


def save_uploaded_files(uploaded_files, fastq_dir):
    """
    Persist uploaded FASTQ files to disk inside the active project's
    FASTQ folder. Streamlit's file_uploader only keeps files in memory;
    downstream tools (FastQC, and alignment) need real files on disk.

    uploaded_files: any iterable of file-like objects exposing `.name`
        and `.getbuffer()` -- satisfied by Streamlit's UploadedFile, but
        not Streamlit-specific itself.
    """
    os.makedirs(fastq_dir, exist_ok=True)
    saved_paths = []
    for f in uploaded_files:
        dest_path = os.path.join(fastq_dir, f.name)
        with open(dest_path, "wb") as out_file:
            out_file.write(f.getbuffer())
        saved_paths.append(dest_path)
    return saved_paths


def find_fastq_filenames_in_directory(source_dir):
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


def symlink_fastq_files_from_directory(source_dir, fastq_dir):
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
    exist in a perfectly usable location on the very same machine this
    app is running on. A symlink is created essentially instantly
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
    _, matching_files = find_fastq_filenames_in_directory(source_dir)

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


def remove_fastq_samples(fastq_dir, sample_pairs, sample_names_to_remove):
    """
    Remove every FASTQ file (R1, R2, and/or SE) belonging to the given
    sample name(s) from fastq_dir -- the mechanism behind the
    interactive workspace's "Remove a sample from this project" action,
    added specifically in response to a real usability-testing finding:
    a non-bioinformatics user had no easy way to remove a sample that
    had been downloaded by mistake (e.g. the wrong SRA accession, or a
    duplicate/unwanted run) short of manually locating and deleting
    files on disk themselves.

    Uses os.remove (not shutil.rmtree or similar) on each individual
    file/symlink path -- this is deliberately SAFE for the symlinked-
    FASTQ case (see symlink_fastq_files_from_directory above):
    os.remove on a symlink deletes only the symlink itself, never the
    file it points to, so removing a sample that was linked in from a
    server-side directory (rather than uploaded/copied) never touches
    or deletes the user's original data elsewhere on disk -- only this
    project's reference to it.

    sample_pairs: the dict returned by validate_sample_pairs (i.e.
        {sample_name: {"R1": filename, ...}}), used to look up exactly
        which filename(s) belong to each sample being removed.
    sample_names_to_remove: list/set of sample names (keys of
        sample_pairs) to remove.

    Returns (removed_filenames: list[str], errors: list[str]) -- the
    filenames that were actually deleted, and any per-file error
    messages encountered (e.g. a permissions issue) -- errors don't
    stop the rest of the removal from proceeding, so one problematic
    file doesn't block removing everything else that CAN be removed.
    """
    removed_filenames = []
    errors = []

    for sample_name in sample_names_to_remove:
        reads = sample_pairs.get(sample_name, {})
        for key in ("R1", "R2", "SE"):
            if key not in reads:
                continue
            filename = reads[key]
            path = os.path.join(fastq_dir, filename)
            try:
                if os.path.islink(path) or os.path.isfile(path):
                    os.remove(path)
                    removed_filenames.append(filename)
            except OSError as e:
                errors.append(f"{filename}: {e}")

    return removed_filenames, errors


def build_match_table(sample_pairs, meta_df):
    """
    Build a plain-language matching table showing, for every sample
    detected from FASTQ filenames AND every sample listed in the metadata
    file, whether they successfully matched up.

    This is the core beginner-friendly output: instead of a cryptic
    error, the user (or an Advanced Mode QC summary) sees exactly which
    samples matched, which FASTQ samples had no metadata row, and which
    metadata rows had no FASTQ files.
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


def write_matched_samplesheet(sample_pairs, meta_df, fastq_dir, samplesheet_path):
    """
    Write out only the samples that fully matched (have both FASTQ files
    on disk AND a metadata row) to a clean samplesheet CSV inside the
    active project. This is the file every downstream step (QC,
    trimming, alignment, quantification) reads from.
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


def detect_blank_metadata_samples(raw_meta_df, sample_col):
    """
    Return the sample name(s) whose metadata row has NOTHING filled in
    besides the sample name itself (every other column is blank/NaN) --
    used to pre-select sensible defaults for an NCBI metadata re-fetch
    action, since these are exactly the rows a user (or an automated
    ingestion step) would want to re-fetch (e.g. rows just added via a
    metadata-sync "Add a blank row for" action, or any other row that's
    never had its details filled in).

    A "blank" cell is either a true NaN or a string that's empty/only
    whitespace after stripping -- covers both a truly missing pandas
    value and a manually-added blank string cell.
    """
    other_cols = [c for c in raw_meta_df.columns if c != sample_col]
    if not other_cols:
        # No other columns exist at all -- every sample is trivially
        # "blank" in this case, though there'd be nothing meaningful
        # for a re-fetch to fill in either.
        return raw_meta_df[sample_col].astype(str).tolist()

    def _is_blank_row(row):
        return all(pd.isna(row[c]) or str(row[c]).strip() == "" for c in other_cols)

    blank_mask = raw_meta_df.apply(_is_blank_row, axis=1)
    return raw_meta_df.loc[blank_mask, sample_col].astype(str).tolist()


def merge_ncbi_metadata_rows(raw_meta_df, sample_col, metadata_rows, overwrite_existing=False):
    """
    Merge NCBI-derived attribute data (from sra_manager.build_metadata_dataframe,
    which returns a list of dicts each keyed by "sample" -- the Run
    accession) into raw_meta_df's EXISTING rows for matching samples,
    rather than adding new rows or replacing the whole table.

    This targets specific, already-existing row(s) and fills in just
    their missing/blank fields -- exactly what's needed when a sample
    was already added to the metadata sheet (e.g. via a sync feature)
    with blank values, and the caller wants to go back and populate
    those blanks from NCBI without disturbing any other already-filled-
    in samples or columns.

    Matching is done by comparing metadata_rows' "sample" value against
    raw_meta_df[sample_col] -- this works correctly because, for any
    sample originally downloaded via this app's SRA feature, the sample
    name IS the Run accession (fasterq-dump's own output naming
    convention, which this app's FASTQ-matching logic already relies on
    elsewhere), so a lookup keyed by Run accession lines up directly
    with the project's existing sample names with no extra translation
    needed.

    overwrite_existing: if False (the default, and the expected common
        case), only currently-BLANK cells are filled in -- any column
        that already has a real value for that sample is left
        completely untouched, even if NCBI's data for it happens to
        differ. If True, every matched column is overwritten with
        NCBI's value regardless of what was there before -- offered as
        an explicit opt-in for a genuine "start over from NCBI" case,
        never the default (since silently overwriting a user's own
        manually-entered/edited values would be a bad, low-trust
        surprise).

    A brand new column present in NCBI's data but not yet in
    raw_meta_df (e.g. "strain", "treatment") is automatically added
    (defaulting to blank for every OTHER, non-matched sample), so
    running this repeatedly for different samples over time correctly
    accumulates whatever attribute columns each individual NCBI lookup
    happens to provide.

    Returns (updated_df, filled_sample_names: list[str]) -- the merged
    DataFrame, and the list of sample names that actually received at
    least one filled-in value (for a clear success message -- a sample
    present in metadata_rows but with, say, every value already
    non-blank and overwrite_existing=False would correctly NOT appear
    in this list, since nothing was actually changed for it).
    """
    updated_df = raw_meta_df.copy()
    filled_samples = []

    for entry in metadata_rows:
        sample_name = entry.get("sample")
        mask = updated_df[sample_col].astype(str) == str(sample_name)
        if not mask.any():
            continue
        idx = updated_df.index[mask][0]

        any_filled = False
        for col, val in entry.items():
            if col == "sample":
                continue
            if col not in updated_df.columns:
                updated_df[col] = pd.Series([pd.NA] * len(updated_df), dtype=object)
            elif not pd.api.types.is_object_dtype(updated_df[col]):
                # A column that's entirely blank for every sample so far
                # (e.g. never filled in yet) gets read back from the CSV
                # as a NUMERIC dtype (float64) by pandas, since an
                # all-empty column has nothing to infer a string type
                # from. Newer/stricter pandas correctly REFUSES to
                # silently upcast that column when a real string value
                # (e.g. "27") is written into it -- raising
                # "TypeError: Invalid value '27' for dtype 'float64'" --
                # rather than the old, looser behavior of quietly
                # converting the whole column on the fly. This explicit
                # conversion must happen BEFORE the assignment below.
                updated_df[col] = updated_df[col].astype(object)

            current_val = updated_df.at[idx, col]
            is_blank = pd.isna(current_val) or str(current_val).strip() == ""
            if overwrite_existing or is_blank:
                updated_df.at[idx, col] = val
                any_filled = True

        if any_filled:
            filled_samples.append(sample_name)

    return updated_df, filled_samples
