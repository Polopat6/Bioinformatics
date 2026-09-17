"""
bulk_bam_manager.py

BAM -> FASTQ recovery/conversion for the Bulk RNA-Seq pipeline, using
`samtools fastq` -- a DIFFERENT tool, solving a DIFFERENT problem, from
single_cell/sc_sra_manager.py's own `bamtofastq`-based recovery path.
Do not confuse the two, and do not assume this module's design mirrors
that one beyond the shared high-level shape (BAM in, FASTQ out).

--- Why this is a genuinely different problem than single-cell's
    BAM recovery, not the same feature reused (2026-08-24) ---
single-cell's bamtofastq recovery exists because SOME 10x BAM tags
(CB/CR/UB/UR -- cell barcode and UMI) are NOT reconstructible from a
plain FASTQ extraction at all; the original BAM is the ONLY place that
information still exists once it's been stripped by SRA's own
normalization. Ordinary bulk RNA-seq reads have no such tags to begin
with -- a bulk BAM's reads are just ordinary aligned (or unaligned)
paired-end cDNA reads, and a completely standard tool
(`samtools fastq`, NOT 10x's own `bamtofastq`) can extract them without
any special barcode-preservation logic.

So this module's real motivating use case is NOT "recovering lost
data" the way single-cell's is -- it's "a user already has a BAM
(received from a collaborator, a public dataset distributed only as
aligned BAM, or their own separate prior pipeline run) and wants to
re-extract FASTQ so it can be run through THIS pipeline's own
ingestion/QC/trimming/alignment steps." This is a genuinely useful,
general capability, but it is not solving the same problem
single-cell's BAM recovery solves, and should not be described or
marketed as such.

--- The one real, non-obvious correctness risk: sort order (2026-08-24)
    ---
`samtools fastq` requires its input BAM to be grouped by READ NAME
(queryname-sorted, or otherwise name-grouped) to correctly reconstruct
R1/R2 mate pairs -- running it on a COORDINATE-sorted BAM (the far more
common on-disk sort order for any BAM that's been through a standard
alignment + indexing pipeline) does NOT produce an error; it silently
produces WRONGLY PAIRED (or entirely unpaired/scrambled) R1/R2 output,
a real, dangerous correctness bug that could go completely unnoticed
until a downstream alignment step produces bizarre-looking results.

This module therefore NEVER assumes an input BAM's sort order is
already correct. get_bam_sort_order() inspects the BAM's own header
(`@HD` line's `SO:` tag) via `samtools view -H` -- if that tag is
explicitly "queryname", the BAM is used as-is; for EVERY other case
(explicitly "coordinate", "unsorted", or the tag being absent/
unparseable entirely) this module defensively re-sorts by name first
(`samtools sort -n`) before ever calling `samtools fastq` -- erring
firmly on the side of a real, current-file's actual sort order,
matching this project's own established design philosophy (see
sc_sra_manager.py's own R1/R2/I1 role-classification heuristics, or
reference_manager.py's own stale-index diagnostics) of never silently
trusting an unconfirmed assumption when getting it wrong would corrupt
downstream results without any visible error.

--- gzip output (2026-08-24) ---
Modern `samtools fastq` (confirmed from samtools' own documented
behavior, available since well before any version this project would
reasonably be pinned to) automatically writes bgzf-compressed output
whenever an output filename passed to -1/-2/-0/-s ends in ".gz" --
no separate compression step is needed the way sc_sra_manager.py's own
_gzip_and_remove() was needed for fasterq-dump's plain output. As a
defensive fallback ONLY (in case an unusually old samtools build on
some system doesn't honor this), convert_bam_to_fastq() checks whether
the expected .gz files actually exist after conversion and, if only
plain, uncompressed files were produced instead, gzips them itself
using the exact same _gzip_and_remove() pattern already established in
sc_sra_manager.py -- so this module never silently leaves an
uncompressed file behind regardless of the installed samtools
version's exact behavior.

--- Output naming: matches ingestion_manager.py's OWN existing
    "_1"/"_2" convention (2026-08-24) ---
ingestion_manager.py's validate_sample_pairs() already recognizes
THREE naming conventions for pairing FASTQ files by sample -- this
module deliberately targets the "SRA/EBI-style" convention
(`<sample>_1.fastq.gz` / `<sample>_2.fastq.gz`), matching that
function's own `re.search(r"_1$", base)` check EXACTLY (checked only
after stripping the .fastq/.fastq.gz extension, so this is unambiguous
regardless of what else the sample name contains). This means a BAM
converted by this module drops directly into a project's existing
fastq_dir() and is picked up by the SAME, already-working
validate_sample_pairs()/build_match_table()/write_matched_samplesheet()
pipeline with ZERO special-casing needed anywhere else in the app --
exactly the same "produce output the existing pipeline already knows
how to consume" principle sc_sra_manager.py's own
finalize_role_assignment() follows for single-cell.

--- Explicitly OUT OF SCOPE for this pass (2026-08-24) ---
Automatically locating a free, original-format BAM for a BULK SRA
accession (the single-cell counterpart being
sc_sra_manager.find_original_format_bam_url()) is NOT built here yet --
doing so correctly would require calling into bulk's own real
sra_manager.py's actual `_esearch_sra()`/`_efetch_sra_full_xml()`
implementations, and this module intentionally does not guess at or
duplicate that module's internals from partial/unconfirmed knowledge.
This module currently only supports converting a BAM the user ALREADY
HAS locally (via upload or an existing path) -- see
save_uploaded_bam()/run_bulk_bam_conversion_from_uploaded_file() below.
Automatic SRA-side BAM discovery for bulk accessions is a reasonable
FUTURE addition once bulk's own real sra_manager.py internals are
directly confirmed, not attempted here from assumption.
"""
import gzip
import os
import shutil

# ---------------------------------------------------------------------------
# Availability + sort-order detection
# ---------------------------------------------------------------------------

def samtools_available():
    "Check whether samtools is available on PATH."
    return shutil.which("samtools") is not None


def get_bam_sort_order(bam_path, subprocess_runner=None):
    """
    Inspect a BAM's own header for its declared sort order, via
    `samtools view -H` -- looks specifically at the `@HD` line's `SO:`
    tag (the BAM spec's own official field for this).

    Returns one of: "queryname", "coordinate", "unknown" (the header
    exists but SO: is absent, or set to something else like
    "unsorted"), or None (samtools itself failed to run, e.g. the file
    isn't a valid BAM at all, or samtools isn't installed).

    Deliberately returns "unknown" (a real, distinct value) rather than
    quietly defaulting to either "queryname" or "coordinate" when the
    tag can't be confirmed -- see this module's own docstring, "The one
    real, non-obvious correctness risk", for why every caller of this
    function treats anything other than a CONFIRMED "queryname" as
    requiring a defensive re-sort.
    """
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    if not os.path.isfile(bam_path):
        return None

    try:
        result = runner(
            ["samtools", "view", "-H", bam_path],
            capture_output=True, text=True, timeout=120,
        )
    except Exception:
        return None

    if result.returncode != 0:
        return None

    for line in (result.stdout or "").splitlines():
        if line.startswith("@HD"):
            for field in line.split("\t"):
                if field.startswith("SO:"):
                    value = field[len("SO:"):].strip().lower()
                    if value in ("queryname", "coordinate"):
                        return value
                    return "unknown"
            return "unknown"  # @HD line present but no SO: field at all

    return "unknown"  # no @HD line at all in the header


def ensure_queryname_sorted(bam_path, work_dir, threads=4, subprocess_runner=None):
    """
    Guarantee the returned BAM path is name-sorted, re-sorting via
    `samtools sort -n` INTO work_dir (never modifying/overwriting the
    original input file in place) whenever get_bam_sort_order() does
    not return a CONFIRMED "queryname" result -- see this module's own
    docstring for why "confirmed queryname, or else re-sort" is the
    only safe default here, rather than trying to optimize away a
    re-sort based on an unconfirmed guess.

    Returns (final_bam_path: str, was_resorted: bool, message: str).
    was_resorted=True means the RETURNED path is a NEW, re-sorted copy
    under work_dir (the original bam_path is left completely
    untouched) -- callers should treat that returned path, not the
    original, as the correct input for convert_bam_to_fastq() below.

    Raises no exception on a samtools failure during the sort itself --
    returns (None, False, error_message) instead, matching this
    project's own established "return an inspectable failure result,
    don't let a subprocess failure propagate as an uncaught exception"
    convention used throughout sc_sra_manager.py.
    """
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    sort_order = get_bam_sort_order(bam_path, subprocess_runner=subprocess_runner)
    if sort_order == "queryname":
        return bam_path, False, "BAM is already confirmed name-sorted -- no re-sort needed."

    os.makedirs(work_dir, exist_ok=True)
    sorted_path = os.path.join(work_dir, os.path.basename(bam_path) + ".namesorted.bam")

    reason = {
        "coordinate": "the BAM's header explicitly declares coordinate sort order",
        "unknown": "the BAM's header does not explicitly confirm queryname sort order",
        None: "the BAM's sort order could not be determined at all (samtools could not read its header)",
    }.get(sort_order, "the BAM's sort order could not be confirmed")

    cmd = ["samtools", "sort", "-n", "-@", str(threads), "-o", sorted_path, bam_path]
    try:
        result = runner(cmd, capture_output=True, text=True, timeout=7200)
    except Exception as e:
        return None, False, f"Failed to re-sort BAM by name: {e}"

    if result.returncode != 0:
        return None, False, f"samtools sort -n failed: {result.stderr}"

    return sorted_path, True, (
        f"Re-sorted BAM by read name before conversion, because {reason} -- "
        f"running samtools fastq directly on a non-name-sorted BAM silently "
        f"produces incorrectly paired R1/R2 output, so this re-sort is always "
        f"performed defensively rather than risking that."
    )


# ---------------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------------

def _gzip_and_remove(src_path):
    "Compress src_path to src_path + '.gz', then delete the original -- same pattern already established in sc_sra_manager.py, used here only as a defensive fallback (see this module's own docstring, 'gzip output')."
    gz_path = src_path + ".gz"
    with open(src_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(src_path)
    return gz_path


def convert_bam_to_fastq(bam_path, dest_dir, sample_name, paired=True, threads=4,
                          subprocess_runner=None, timeout=7200):
    """
    Run `samtools fastq` against an ALREADY name-sorted BAM (callers
    should run ensure_queryname_sorted() first -- see
    run_bulk_bam_conversion() below for the orchestrated version that
    does this automatically) and produce output named to match
    ingestion_manager.py's own existing "_1"/"_2" convention.

    paired: if True (the default), extracts R1/R2 as
        "<sample_name>_1.fastq.gz" / "<sample_name>_2.fastq.gz".
        Unpaired/singleton reads (present in the BAM but lacking a
        mapped mate) are explicitly discarded (`-0 /dev/null -s
        /dev/null`) rather than silently left out of the produced file
        list with no explanation -- a bulk RNA-seq samplesheet expects
        a clean R1/R2 pair per sample, and singleton reads have no
        well-defined place in that convention; a user who specifically
        needs those reads preserved should be aware they are not
        included in this conversion path today.
        If False, extracts every read as a single unpaired file
        ("<sample_name>.fastq.gz").

    Returns (success: bool, message: str, produced_files: list[str]).
    """
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    if not os.path.isfile(bam_path):
        return False, f"BAM file not found: {bam_path}", []

    os.makedirs(dest_dir, exist_ok=True)

    if paired:
        r1_path = os.path.join(dest_dir, f"{sample_name}_1.fastq.gz")
        r2_path = os.path.join(dest_dir, f"{sample_name}_2.fastq.gz")
        cmd = [
            "samtools", "fastq", "-@", str(threads),
            "-1", r1_path, "-2", r2_path,
            "-0", os.devnull, "-s", os.devnull,
            "-n", bam_path,
        ]
        expected_files = [r1_path, r2_path]
    else:
        se_path = os.path.join(dest_dir, f"{sample_name}.fastq.gz")
        cmd = ["samtools", "fastq", "-@", str(threads), "-0", se_path, "-n", bam_path]
        expected_files = [se_path]

    try:
        result = runner(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess_module.TimeoutExpired:
        return False, f"samtools fastq timed out after {timeout} seconds.", []
    except Exception as e:
        return False, f"samtools fastq failed to run: {e}", []

    if result.returncode != 0:
        return False, f"samtools fastq failed: {result.stderr}", []

    # Defensive fallback (see this module's own docstring, "gzip
    # output") -- only triggers on an unusually old samtools build that
    # didn't honor the .gz suffix; the normal case is that these files
    # already exist exactly as requested.
    final_files = []
    for expected_path in expected_files:
        if os.path.isfile(expected_path):
            final_files.append(expected_path)
            continue
        plain_path = expected_path[:-len(".gz")]  # strip trailing ".gz"
        if os.path.isfile(plain_path):
            final_files.append(_gzip_and_remove(plain_path))

    if not final_files:
        return False, (
            f"samtools fastq completed without an error, but none of the expected output "
            f"file(s) were found ({', '.join(expected_files)})."
        ), []

    return True, f"Converted {bam_path} to {len(final_files)} FASTQ file(s) for sample '{sample_name}'.", final_files


def run_bulk_bam_conversion(bam_path, project_fastq_dir, sample_name, paired=True,
                             threads=4, progress_callback=None, subprocess_runner=None):
    """
    End-to-end orchestration: check/ensure name-sort order, convert to
    FASTQ, and clean up any intermediate re-sorted BAM copy -- the
    single function most callers (UI or otherwise) should actually use,
    rather than calling ensure_queryname_sorted()/convert_bam_to_fastq()
    separately.

    Writes output DIRECTLY into project_fastq_dir using
    ingestion_manager.py's own recognized "_1"/"_2" naming convention --
    unlike single-cell's own BAM recovery path (which writes into a
    scratch _sra_bam_work/ subdirectory and requires a separate
    classify-then-finalize step, since single-cell's ROLE ambiguity --
    which file is barcode vs. cDNA -- has no bulk equivalent at all;
    bulk's own R1/R2 roles are unambiguous the moment they're
    extracted), this function's output is immediately usable by the
    existing ingestion pipeline with NO further step required.

    Returns (success: bool, message: str, produced_files: list[str]).
    """
    if not samtools_available():
        return False, "`samtools` was not found on PATH -- install it first (conda-forge/bioconda package: samtools).", []

    work_dir = os.path.join(project_fastq_dir, "_bam_conversion_work", sample_name)
    os.makedirs(work_dir, exist_ok=True)

    if progress_callback:
        progress_callback(f"Checking sort order for {sample_name}...")

    sorted_bam_path, was_resorted, sort_message = ensure_queryname_sorted(
        bam_path, work_dir, threads=threads, subprocess_runner=subprocess_runner,
    )
    if sorted_bam_path is None:
        return False, sort_message, []

    if progress_callback:
        progress_callback(f"{sort_message} ({sample_name})")
        progress_callback(f"Converting {sample_name} to FASTQ (samtools fastq)...")

    success, message, produced_files = convert_bam_to_fastq(
        sorted_bam_path, project_fastq_dir, sample_name, paired=paired,
        threads=threads, subprocess_runner=subprocess_runner,
    )

    # Clean up the intermediate re-sorted BAM copy (never the ORIGINAL
    # bam_path, which this function never modifies or deletes) --
    # regardless of success/failure, so a failed conversion doesn't
    # leave a large, orphaned re-sorted BAM copy sitting in the
    # project's own directory tree indefinitely.
    if was_resorted and os.path.isfile(sorted_bam_path):
        os.remove(sorted_bam_path)

    return success, message, produced_files


# ---------------------------------------------------------------------------
# Upload-your-own-BAM support
# ---------------------------------------------------------------------------

def save_uploaded_bam(uploaded_file, dest_dir, sample_name):
    """
    Save a Streamlit-uploaded BAM file object to disk -- same pattern
    as single_cell/sc_sra_manager.py's own save_uploaded_bam(), kept as
    a separate copy in this module (rather than a shared import between
    the two pipelines) since bulk_bam_manager.py and sc_sra_manager.py
    are otherwise intentionally independent modules with no coupling
    between the bulk and single-cell pipelines' own dependency graphs.

    uploaded_file: a Streamlit UploadedFile object (or any file-like
        object exposing .name and .read()).
    dest_dir: directory to save the uploaded BAM into.
    sample_name: used to build a predictable destination filename
        (f"{sample_name}.uploaded.bam") so a later re-upload for the
        same sample cleanly overwrites rather than accumulating
        differently-named duplicate BAM files.

    Returns the saved file's full path.
    """
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, f"{sample_name}.uploaded.bam")
    with open(dest_path, "wb") as f:
        f.write(uploaded_file.read())
    return dest_path


def run_bulk_bam_conversion_from_uploaded_file(uploaded_bam_path, project_fastq_dir, sample_name,
                                                paired=True, threads=4, progress_callback=None,
                                                subprocess_runner=None):
    """
    Thin, explicitly-named wrapper around run_bulk_bam_conversion() for
    the "user already has this BAM on hand" case (uploaded via the UI,
    or otherwise already sitting on disk) -- kept as its own function
    (identical in behavior to calling run_bulk_bam_conversion()
    directly) purely so a caller reading this module's own public
    function list can immediately tell which entry point corresponds to
    "convert an uploaded/local BAM" without needing to infer it from a
    more generically-named function.
    """
    return run_bulk_bam_conversion(
        uploaded_bam_path, project_fastq_dir, sample_name, paired=paired,
        threads=threads, progress_callback=progress_callback, subprocess_runner=subprocess_runner,
    )
