"""
single_cell/sc_sra_manager.py

NCBI/SRA fetching for the Single-cell RNA-Seq pipeline's Step 1, adapted
from (and reusing) the Bulk RNA-Seq pipeline's sra_manager.py -- but with
single-cell-specific handling this module adds on top, because single-cell
SRA depositions are NOT as simple as bulk's "_1.fastq = R1, _2.fastq = R2"
convention.

--- Why this needs its own module rather than just reusing sra_manager
    directly ---
Accession VALIDATION (esearch/efetch lookup, RNA-Seq library-strategy
checking, sample-metadata extraction) is IDENTICAL to bulk and is reused
here directly via `import sra_manager as sra`.

What's different, and why this module exists:

1. FILE ROLE IS AMBIGUOUS AFTER DOWNLOAD -- classify_output_files() below
   samples read lengths to guess R1/R2/I1 roles, always surfaced for user
   confirmation, never silently trusted.

2. SOME DEPOSITIONS ARE BAM, NOT FASTQ, AND NEED A DIFFERENT TOOL
   ENTIRELY -- see "Original-format BAM recovery" section below.

3. DOWNLOAD EXECUTION IS IMPLEMENTED HERE, NOT REUSED FROM sra_manager.

--- Batch/parallel download support (2026-08-24) ---
download_and_classify_runs_parallel() is a direct port of bulk
sra_manager.py's own download_sra_runs_parallel() concurrency pattern.

--- Missing gzip compression fix (2026-08-24) ---
download_and_classify_run() now gzips its fasterq-dump output (a direct
port of bulk sra_manager.py's own _gzip_and_remove() helper).

--- Original-format BAM recovery: real, honest scope (2026-08-24) ---
Some single-cell SRA depositions have NO usable barcode/UMI information
in their standard FASTQ extraction at all. In these cases, the ONLY way
to recover a usable FASTQ is to obtain the run's ORIGINAL (not SRA-
normalized/"ETL") BAM file and convert it back to FASTQ using 10x's own
`bamtofastq` tool. find_original_format_bam_url() inspects the run's OWN
SRA metadata XML for a real <SRAFile supertype="Original" ...> entry
rather than guessing/constructing a URL, and reports honestly if none is
found.

--- v1-chemistry bamtofastq output resolution (2026-08-24, confirmed
    against REAL recovered data) ---
CONFIRMED via a real, successful recovery of Kang et al. 2018 (GSE96583):
for BAMs from Cell Ranger 1.0/1.1 (the --cr11 bamtofastq flag),
bamtofastq's own output is a FOUR-FILE-PER-LANE split: R1=cDNA(98bp),
R2=barcode(14bp), R3=UMI(10bp), I1=sample index(8bp). STARsolo's own
barcode model requires ONE single barcode+UMI read per lane.
resolve_bamtofastq_v1_split() detects this exact signature and
automatically concatenates each lane's R2+R3 into a single combined
read, producing a plain 2-file-per-lane result classify_output_files()
can classify normally.

This resolution runs automatically inside BOTH
run_full_bam_recovery_pipeline() (SRA-based recovery) and
run_bam_recovery_from_uploaded_file() (upload/browse-your-own-BAM path)
-- both converge on the exact same downstream functions, so a v1-era BAM
is handled identically regardless of how it was obtained.

--- bamtofastq "quiet success, zero files" diagnostic gap fix
    (2026-08-24, later same day) ---
A real reported bug: run_bamtofastq() previously only captured/returned
bamtofastq's own stdout/stderr when the tool exited with a NON-ZERO
return code -- if bamtofastq exited 0 (success) but produced ZERO FASTQ
files (a real, confirmed-possible outcome, e.g. when an explicit
chemistry flag mismatches the BAM's actual origin), callers had ZERO
visibility into what bamtofastq itself actually printed, since only a
generic "completed successfully" message was ever returned in the
success case.

Fixed by having run_bamtofastq() ALWAYS capture and return bamtofastq's
own stdout+stderr (not just on failure), and by having
run_full_bam_recovery_pipeline()/run_bam_recovery_from_uploaded_file()
explicitly check for this "success but zero files" case and surface
that captured tool output directly in the returned error message.

--- "Already downloaded" detection gap fix (2026-08-31) ---
A real reported bug: the covid_test single-cell project re-attempted
prefetch/fasterq-dump for accessions that had already been successfully
downloaded in a previous session. Root cause: this module previously had
NO detection logic at all for "has this accession already been
downloaded?" -- download_and_classify_run() unconditionally re-ran
prefetch + fasterq-dump for every accession, every time it was called.

This is made worse by finalize_role_assignment() renaming raw,
accession-named output files (e.g. "SRR13734390_1.fastq.gz") into the
project's standard sample-based 10x naming -- which, for MOST callers,
would no longer contain the original accession string anywhere. Once a
run is finalized, filenames alone can no longer answer "which accession
produced this file?" in general.

Fixed with THREE detection tiers, in priority order:
  1. Registry (_sra_download_registry.json, in the project's fastq dir --
     same pattern as bulk's own processed_registry.json) -- the durable,
     general-purpose source of truth going forward, updated by
     mark_accessions_finalized() (which callers must invoke once, right
     after a successful finalize_role_assignment() call).
  2. Raw staging directory fallback (_sra_work/<accession>/) -- for a run
     downloaded but not yet finalized/moved.
  3. Filename-pattern fallback (2026-08-31, added same day) -- SPECIFIC
     to singlecell_workspace.py's own _render_sra_source() UI, which
     calls finalize_role_assignment(..., sample_name) with sample_name
     SET TO THE ACCESSION ITSELF for the SRA-download path -- meaning a
     PREVIOUSLY finalized SRA run's files (from before this registry
     existed, so no registry entry was ever written for them) are still
     named like "SRR13734390_S1_R1_001.fastq.gz": the accession string is
     STILL the literal filename prefix, even after finalization, for
     this specific caller. This tier is intentionally an exact-prefix
     match (never a loose substring search) to avoid any false-positive
     risk against an unrelated sample. It does NOT apply to the "convert
     from a BAM file I already have" path, since that lets the user
     choose an arbitrary sample_name unrelated to any SRA accession.

A registry hit (tier 1) is NEVER trusted blindly -- at least one file it
references must still actually exist on disk, or the entry is treated as
stale and ignored, so a registry out of sync with reality (e.g. someone
manually deleted files) can't cause a silent false-positive.

download_and_classify_run() and download_and_classify_runs_parallel()
both now accept an optional force_redownload=False parameter; when
False (the default), an accession already found via ANY of the three
tiers above is skipped with an explicit "skipped" result rather than
re-downloaded.

backfill_registry_from_existing_files() is a new, ONE-TIME-USE utility
that proactively scans a project's fastq dir for tier-3
(filename-pattern) matches across a whole list of accessions and writes
registry entries for every match found, so that subsequent lookups can
hit the fast, general-purpose tier-1 registry check directly rather
than re-scanning the filesystem by pattern every single time. This is
the recommended fix for existing projects (like covid_test) that
already have finalized runs on disk from BEFORE this fix existed.

--- "prefetch quiet success, no .sra file" diagnostic fix (2026-09-01)
---
A real, confirmed bug (CONFIRMED root cause of a real download failure
on SRR13734386, a 72GB run): prefetch's own default maximum download
size is 20GB. When an accession exceeds this, prefetch does NOT error --
it exits 0 (success) while explicitly SKIPPING the download ("Download
of some files was skipped because they are too large"), which then
produced a confusing, unrelated-looking downstream fasterq-dump VFS/404
error ("Cannot resolve accession") with no indication the REAL cause was
a silent size-limit skip one step earlier.

Fixed two ways: (1) build_prefetch_command() now always passes an
explicit --max-size (100G by default, raised from prefetch's own 20GB
default, and overridable via the new max_size parameter threaded through
download_and_classify_run()/download_and_classify_runs_parallel()); (2)
download_and_classify_run() now verifies a real .sra file actually
exists immediately after prefetch returns (regardless of its exit code),
and if not, returns a clear diagnostic message -- specifically detecting
and calling out the size-limit-skip signature in prefetch's own captured
output when that's the cause, rather than silently falling through to
fasterq-dump's much more confusing, seemingly-unrelated error.

--- Proactive resume scan + stale-cache refresh caveat (2026-09-01)
---
scan_pending_unfinalized_runs() (below) is a proactive filesystem scan
-- completely independent of Streamlit session_state -- that finds every
accession with real, already-downloaded raw FASTQ output still sitting
in its own _sra_work/<accession>/ staging directory, not yet finalized.
The intended caller pattern (see singlecell_workspace.py's own
_render_sra_source()) is to call this ONCE per render and MERGE its
result into st.session_state["sc_sra_download_results"], so a
previously-downloaded-but-unconfirmed run reappears automatically in the
UI the moment a project is reopened, with no need to re-validate/
re-select anything.

*** REAL, CONFIRMED CALLER-SIDE BUG this surfaced (2026-09-01), fixed in
    singlecell_workspace.py itself, NOT in this module ***: the
UI-layer merge must NOT use a plain dict.setdefault()-style "only add if
the key doesn't already exist" merge. Confirmed via a real reported case
(SRR13734384, a 348GB run): this module's own scan can genuinely be
called and cached MULTIPLE TIMES across a single large accession's
download lifetime -- e.g. once when only its small 8bp I1 file had
finished downloading (correctly reporting "could not classify" for that
1-file snapshot at the time), and again later once its two ~348GB R1/R2
files had ALSO finished. Every scan_pending_unfinalized_runs() result is
always tagged skipped=True by construction (see below) -- so it is safe,
and NECESSARY, for the UI-layer merge to overwrite ANY existing cached
entry that is ALSO tagged skipped=True (since that entry can only have
come from this same resume-scan path previously, and therefore may be a
stale, incomplete snapshot from earlier in the same download), while
still NEVER overwriting a real, just-completed LIVE download result from
the current session (tagged skipped=False), which must always be
preserved untouched. A caller using plain setdefault() will appear to
work initially but will silently freeze on a run's FIRST scanned
snapshot forever, even after the real underlying download completes and
more files appear on disk -- this exact symptom (a run stuck showing
only 1 of its real 3 files, indefinitely, across page reloads) is what
originally surfaced this caveat.
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import sra_manager as sra  # reused directly -- see docstring
import chemistry_manager as chem

ROLE_R1 = "R1"   # cell barcode + UMI (short, fixed length)
ROLE_R2 = "R2"   # cDNA (longer, variable-ish)
ROLE_I1 = "I1"   # sample index read (very short, typically 8-10bp) -- not used downstream
ROLE_UNKNOWN = "UNKNOWN"

# ---------------------------------------------------------------------------
# Assay-type classification (bulk vs. single-cell vs. spatial) from SRA
# metadata -- NOT from LibraryStrategy
# ---------------------------------------------------------------------------
_POSSIBLE_TEXT_FIELDS = (
    "Title", "ExperimentTitle", "experiment_title", "SampleTitle", "sample_title",
    "StudyTitle", "study_title", "Description", "description", "SampleName",
    "sample_name", "LibraryName", "library_name", "SRAStudyTitle",
)

_SINGLE_CELL_KEYWORDS = (
    "single cell", "single-cell", "single nucleus", "single-nucleus", "scrna", "sc-rna",
    "snrna-seq", "sn-rna", "10x genomics", "10x chromium", "chromium", "droplet-based",
    "drop-seq", "dropseq", "indrops", "in-drop", "bd rhapsody", "rhapsody", "split-seq",
    "splitseq", "sci-rna-seq", "microwell-seq", "seq-well", "cel-seq", "celseq",
    "mars-seq", "marsseq",
)
_SPATIAL_KEYWORDS = (
    "visium", "slide-seq", "slideseq", "merfish", "xenium", "geomx", "cosmx",
    "spatial transcriptomics", "spatial rna", "spatial gene expression", "st-seq",
    "spatial-seq", "spatial barcoding",
)
_BULK_HINT_KEYWORDS = (
    "bulk rna", "bulk-seq", "bulk rna-seq", "total rna-seq", "poly-a selected",
    "ribo-depleted", "ribodepletion", "ribo-zero",
)

ASSAY_SINGLE_CELL = "single_cell"
ASSAY_SPATIAL = "spatial"
ASSAY_BULK = "bulk"
ASSAY_UNKNOWN = "insufficient_info"


def classify_assay_type_from_metadata(row):
    combined_text = " ".join(
        str(row.get(field, "")) for field in _POSSIBLE_TEXT_FIELDS if row.get(field)
    ).lower()
    if not combined_text.strip():
        return {
            "assay_type": ASSAY_UNKNOWN, "matched_keyword": None, "matched_field": None,
            "message": (
                "No title/description text field was available in this "
                "accession's metadata to check -- library strategy alone "
                "cannot distinguish bulk from single-cell/spatial data. "
                "Consider checking this accession's own SRA/GEO page "
                "manually before proceeding."
            ),
        }

    for keyword in _SPATIAL_KEYWORDS:
        if keyword in combined_text:
            return {
                "assay_type": ASSAY_SPATIAL, "matched_keyword": keyword, "matched_field": "title/description",
                "message": f"⚠️ This looks like **spatial transcriptomics** data (matched '{keyword}'), not droplet-based single-cell -- this pipeline does not support spatial data.",
            }
    for keyword in _SINGLE_CELL_KEYWORDS:
        if keyword in combined_text:
            return {
                "assay_type": ASSAY_SINGLE_CELL, "matched_keyword": keyword, "matched_field": "title/description",
                "message": f"✅ This looks like single-cell data (matched '{keyword}').",
            }
    for keyword in _BULK_HINT_KEYWORDS:
        if keyword in combined_text:
            return {
                "assay_type": ASSAY_BULK, "matched_keyword": keyword, "matched_field": "title/description",
                "message": f"⚠️ This looks like **bulk** RNA-seq data (matched '{keyword}'), not single-cell -- double-check before downloading here.",
            }
    return {
        "assay_type": ASSAY_UNKNOWN, "matched_keyword": None, "matched_field": None,
        "message": (
            "Title/description text was available but didn't match any "
            "known single-cell, spatial, or bulk keyword -- inconclusive."
        ),
    }


_I1_MAX_LEN = 12


def classify_output_files(fastq_paths, n_reads=2000):
    lengths = {}
    for path in fastq_paths:
        lengths[path] = chem.dominant_read_length(path, n_reads=n_reads)

    result = {p: {"role": ROLE_UNKNOWN, "length": lengths[p], "chemistry_candidates": []} for p in fastq_paths}

    valid_paths = [p for p in fastq_paths if lengths[p] is not None]
    if len(valid_paths) not in (2, 3):
        return result

    by_length = sorted(valid_paths, key=lambda p: lengths[p])

    if len(by_length) == 3:
        i1_path, r1_path, r2_path = by_length[0], by_length[1], by_length[2]
        if lengths[i1_path] > _I1_MAX_LEN:
            return result
        result[i1_path]["role"] = ROLE_I1
    else:
        r1_path, r2_path = by_length[0], by_length[1]

    r1_len = lengths[r1_path]
    r1_candidates = chem._R1_LEN_INDEX.get(r1_len, [])
    if not r1_candidates:
        return result

    result[r1_path]["role"] = ROLE_R1
    result[r1_path]["chemistry_candidates"] = r1_candidates
    result[r2_path]["role"] = ROLE_R2
    return result


def detect_likely_bam_derived_issue(classification):
    if not classification:
        return None
    roles = {info["role"] for info in classification.values()}
    if roles == {ROLE_UNKNOWN}:
        return (
            "⚠️ Could not confidently identify cell-barcode/UMI vs. cDNA reads "
            "from this run's downloaded file(s). This can happen when an "
            "accession was originally deposited as a 10x Genomics BAM file "
            "rather than plain FASTQ, or when only the cDNA read was ever "
            "deposited to SRA at all (a known real gap for some older, "
            "v1-chemistry-era 10x depositions) -- OR, importantly, when NOT "
            "ALL of this run's expected files have finished downloading yet "
            "(e.g. only a small I1 index file is present so far, out of an "
            "expected 2-3 file set) -- always double-check this run's real, "
            "current file listing on disk before assuming this is a genuine "
            "BAM-derived issue rather than an in-progress download. Check "
            "this accession's own SRA page, or try this module's "
            "find_original_format_bam_url() / run_full_bam_recovery_pipeline() "
            "to see whether the ORIGINAL Cell Ranger BAM (which retains "
            "barcode/UMI tags a plain FASTQ extraction does not) is directly "
            "recoverable for this run."
        )
    return None


DEFAULT_PREFETCH_MAX_SIZE = "100G"


def build_prefetch_command(accession, output_dir, max_size=DEFAULT_PREFETCH_MAX_SIZE):
    """
    --- 20GB default download-limit fix (2026-09-01) ---
    A real, confirmed bug: prefetch's own default maximum download size
    is 20GB. When an accession exceeds this, prefetch does NOT error --
    it exits 0 (success) while explicitly SKIPPING the download
    ("Download of some files was skipped because they are too large"),
    which then caused a confusing downstream fasterq-dump VFS/404 error
    with no indication the real cause was a silent size-limit skip one
    step earlier. Confirmed real-world trigger: SRR13734386 (72GB), part
    of a COVID single-cell PBMC batch where run sizes vary considerably
    (later confirmed up to ~350GB per file for some runs in this same
    batch).

    Fixed by always passing an explicit --max-size, defaulting to 100GB.
    Callers needing a different ceiling (e.g. a known-huge deposition, or
    a deliberately LOWER limit to avoid accidentally filling scratch
    disk) can override max_size directly.
    """
    return [
        "prefetch", accession, "--output-directory", output_dir,
        "--max-size", max_size,
    ]


def build_fasterq_dump_command(accession, sra_file_dir, output_dir, threads=4):
    return [
        "fasterq-dump", "--split-files", "--include-technical",
        "--threads", str(threads), "--outdir", output_dir,
        os.path.join(sra_file_dir, accession),
    ]


def _pigz_available():
    "Check whether pigz (parallel gzip) is installed -- a much faster, multi-threaded alternative to Python's own single-threaded gzip module, especially valuable on a many-core HPC node."
    return shutil.which("pigz") is not None


def _gzip_and_remove(src_path, threads=4):
    """
    Compress src_path to src_path + '.gz', then delete the original.
    Uses pigz (parallel gzip) when installed, falling back to Python's
    own gzip module at compresslevel=6 (matching standard command-line
    gzip's own default) otherwise.
    """
    gz_path = src_path + ".gz"

    if _pigz_available():
        with open(gz_path, "wb") as f_out:
            result = subprocess.run(
                ["pigz", "-p", str(threads), "-c", src_path],
                stdout=f_out, stderr=subprocess.PIPE,
            )
        if result.returncode == 0:
            os.remove(src_path)
            return gz_path
        if os.path.exists(gz_path):
            os.remove(gz_path)

    with open(src_path, "rb") as f_in, gzip.open(gz_path, "wb", compresslevel=6) as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(src_path)
    return gz_path


def _find_prefetched_sra_file(accession, work_dir):
    """
    Look for the .sra file prefetch is expected to have produced.
    prefetch's own convention (with --output-directory <work_dir>) is to
    write <work_dir>/<accession>/<accession>.sra -- but different
    sra-tools versions/configurations have occasionally been observed to
    place it directly at <work_dir>/<accession>.sra instead, so both
    locations are checked.
    """
    nested_path = os.path.join(work_dir, accession, f"{accession}.sra")
    if os.path.isfile(nested_path):
        return nested_path
    flat_path = os.path.join(work_dir, f"{accession}.sra")
    if os.path.isfile(flat_path):
        return flat_path
    return None


# ---------------------------------------------------------------------------
# "Already downloaded" detection (2026-08-31)
# ---------------------------------------------------------------------------
DOWNLOAD_REGISTRY_FILENAME = "_sra_download_registry.json"


def _registry_path(project_fastq_dir):
    return os.path.join(project_fastq_dir, DOWNLOAD_REGISTRY_FILENAME)


def _load_download_registry(project_fastq_dir):
    path = _registry_path(project_fastq_dir)
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        # A corrupt/unreadable registry must NEVER be silently trusted as
        # "everything is downloaded" or crash the caller -- treat it as
        # empty (forcing fresh, honest re-detection for every accession)
        # rather than guessing at its intended contents.
        return {}


def _save_download_registry(project_fastq_dir, registry):
    os.makedirs(project_fastq_dir, exist_ok=True)
    path = _registry_path(project_fastq_dir)
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump(registry, f, indent=2, sort_keys=True)
    os.replace(tmp_path, path)  # atomic write -- avoids ever leaving a torn/partial registry file behind


def _now_iso():
    return datetime.now(timezone.utc).isoformat()


def _find_raw_staged_files(accession, project_fastq_dir):
    """
    Tier 2: look for this accession's own already-downloaded raw output
    still sitting in its _sra_work/<accession>/ staging directory -- i.e.
    downloaded (and possibly classified) previously, but not yet
    finalized/moved into the project's main fastq dir.
    """
    work_dir = os.path.join(project_fastq_dir, "_sra_work", accession)
    if not os.path.isdir(work_dir):
        return []
    return sorted([
        os.path.join(work_dir, f) for f in os.listdir(work_dir)
        if f.startswith(accession) and (f.endswith(".fastq") or f.endswith(".fastq.gz"))
    ])


def _find_finalized_files_by_naming_convention(accession, project_fastq_dir):
    """
    Tier 3: filesystem-pattern fallback for accessions that were ALREADY
    downloaded, classified, AND finalized/moved by a PRIOR session --
    i.e. before this module's registry (and mark_accessions_finalized())
    existed, so no registry entry was ever written for them.

    This works specifically because singlecell_workspace.py's own
    _render_sra_source() calls finalize_role_assignment(..., sample_name)
    with sample_name SET TO THE ACCESSION ITSELF for the SRA-download
    path -- so finalize_role_assignment()'s own standard 10x-style
    naming convention ("{sample_name}_S1_{role}_001.ext" or
    "{sample_name}_S1_L{lane}_{role}_001.ext") means a previously
    finalized SRA run's files are named like
    "SRR13734390_S1_R1_001.fastq.gz" -- i.e. the accession string is
    STILL the literal filename prefix, even after finalization.

    Intentionally an EXACT prefix match (never a loose substring search)
    to avoid any false-positive risk against an unrelated sample. Does
    NOT apply to the "convert from a BAM file I already have" path,
    since that path lets the user choose an arbitrary sample_name
    unrelated to any SRA accession.
    """
    if not os.path.isdir(project_fastq_dir):
        return []
    prefix_flat = f"{accession}_S1_"
    prefix_laned = f"{accession}_S1_L"
    matches = []
    for f in os.listdir(project_fastq_dir):
        if not (f.endswith(".fastq") or f.endswith(".fastq.gz")):
            continue
        if f.startswith(prefix_flat) or f.startswith(prefix_laned):
            matches.append(os.path.join(project_fastq_dir, f))
    return sorted(matches)


def is_accession_already_downloaded(accession, project_fastq_dir):
    """
    Determine whether `accession` has already been successfully
    downloaded, checking THREE tiers in order (see module docstring for
    full rationale of each):

      1. registry (_sra_download_registry.json)
      2. raw staging directory (_sra_work/<accession>/)
      3. filename-pattern match in the project's main fastq dir
         ("<accession>_S1_..." -- catches runs finalized before this
         registry existed)

    A registry hit (tier 1) is NEVER trusted blindly -- at least one
    file it references must still actually exist on disk, or the entry
    is treated as stale and ignored.

    Returns a dict:
        {"already_downloaded": bool, "source": "registry" | "raw_staging" | "finalized_naming" | None,
         "existing_files": [...], "detail": "<human-readable explanation>"}
    """
    registry = _load_download_registry(project_fastq_dir)
    entry = registry.get(accession)
    if entry:
        candidate_files = entry.get("finalized_files") or entry.get("raw_files") or []
        still_present = [p for p in candidate_files if os.path.exists(p)]
        if still_present:
            return {
                "already_downloaded": True,
                "source": "registry",
                "existing_files": still_present,
                "detail": (
                    f"'{accession}' is recorded in the download registry as "
                    f"'{entry.get('status', 'downloaded')}' (sample "
                    f"'{entry.get('sample_name', '?')}'), and {len(still_present)} of its "
                    f"recorded file(s) are still present on disk."
                ),
            }
        # Registry says downloaded, but every file it references is now
        # missing -- fall through to the remaining tiers rather than
        # trusting a stale entry.

    raw_files = _find_raw_staged_files(accession, project_fastq_dir)
    if raw_files:
        return {
            "already_downloaded": True,
            "source": "raw_staging",
            "existing_files": raw_files,
            "detail": (
                f"'{accession}' already has {len(raw_files)} raw downloaded FASTQ file(s) "
                f"sitting in its own _sra_work/{accession}/ staging directory from a previous "
                f"run, but this accession is not (or is no longer) in the download registry -- "
                f"these files have not been finalized/moved yet."
            ),
        }

    finalized_files = _find_finalized_files_by_naming_convention(accession, project_fastq_dir)
    if finalized_files:
        return {
            "already_downloaded": True,
            "source": "finalized_naming",
            "existing_files": finalized_files,
            "detail": (
                f"'{accession}' already has {len(finalized_files)} finalized FASTQ file(s) "
                f"in this project's fastq directory matching this accession's standard naming "
                f"pattern ('{accession}_S1_...'), from a run finalized before this project's "
                f"download registry was introduced."
            ),
        }

    return {
        "already_downloaded": False,
        "source": None,
        "existing_files": [],
        "detail": f"No existing downloaded output found for '{accession}'.",
    }


def mark_accessions_finalized(project_fastq_dir, accession_to_finalized_files):
    """
    Record, in the download registry, that one or more accessions'
    output has been finalized (moved into the project's main fastq dir
    by finalize_role_assignment()).

    Callers should invoke this AFTER a successful finalize_role_assignment()
    call, passing a dict mapping each accession involved to the list of
    final destination file paths finalize_role_assignment() returned for
    it (e.g. {"SRR13734390": ["/path/SRR13734390_S1_R1_001.fastq.gz",
    "/path/SRR13734390_S1_R2_001.fastq.gz"]}).
    """
    registry = _load_download_registry(project_fastq_dir)
    for accession, finalized_files in accession_to_finalized_files.items():
        entry = registry.get(accession, {})
        entry["status"] = "finalized"
        entry["finalized_files"] = list(finalized_files)
        entry["finalized_timestamp"] = _now_iso()
        registry[accession] = entry
    _save_download_registry(project_fastq_dir, registry)


def backfill_registry_from_existing_files(project_fastq_dir, accessions):
    """
    ONE-TIME-USE utility for existing projects (like covid_test) that
    already have finalized SRA runs on disk from BEFORE this registry
    existed. Proactively scans for tier-3 (filename-pattern) matches
    across the given list of accessions and writes registry entries for
    every match found -- so future lookups hit the fast, general-purpose
    registry check directly rather than re-scanning the filesystem by
    pattern every time.

    accessions: an iterable of SRA run accessions to check (e.g. every
        accession you've ever attempted to download for this project --
        it's harmless to include ones that were never actually
        downloaded; they're simply skipped).

    Returns a dict summarizing what was found:
        {"backfilled": [accessions with files found and now registered],
         "not_found": [accessions with no matching files on disk]}

    Safe to run multiple times -- re-running simply re-confirms/refreshes
    already-backfilled entries rather than duplicating anything.
    """
    backfilled = []
    not_found = []
    registry = _load_download_registry(project_fastq_dir)

    for accession in accessions:
        # Don't bother re-scanning if a live registry entry already
        # exists and still points at real files.
        existing = is_accession_already_downloaded(accession, project_fastq_dir)
        if existing["already_downloaded"] and existing["source"] == "registry":
            backfilled.append(accession)
            continue

        finalized_files = _find_finalized_files_by_naming_convention(accession, project_fastq_dir)
        if finalized_files:
            entry = registry.get(accession, {})
            entry["status"] = "finalized"
            entry["finalized_files"] = finalized_files
            entry["finalized_timestamp"] = _now_iso()
            entry.setdefault("sample_name", accession)
            entry["backfilled"] = True
            registry[accession] = entry
            backfilled.append(accession)
        else:
            not_found.append(accession)

    _save_download_registry(project_fastq_dir, registry)
    return {"backfilled": backfilled, "not_found": not_found}


def scan_pending_unfinalized_runs(project_fastq_dir):
    """
    Proactive "resume" scan for previously-downloaded-but-unconfirmed
    runs (2026-09-01). Scans _sra_work/ directly for ANY accession
    subdirectory with real raw FASTQ output not yet finalized --
    completely independent of session_state or of which accessions the
    user has typed/selected in the current session.

    Intended caller pattern (see singlecell_workspace.py's own
    _render_sra_source()): call this ONCE per render and MERGE its
    result into st.session_state["sc_sra_download_results"], so every
    previously downloaded-but-unconfirmed run reappears automatically,
    ready to confirm, the moment the project is reopened.

    Returns a dict of {accession: result_dict}, in EXACTLY the same
    shape download_and_classify_run() itself returns for a successful
    (non-skipped) run, EXCEPT this function's own results are ALWAYS
    tagged "skipped": True (see the important caller-side caveat below).

    *** IMPORTANT CALLER-SIDE CAVEAT (2026-09-01, confirmed via a real
        reported bug on SRR13734384, a 348GB run) ***
    This function can be called and its result cached MULTIPLE TIMES
    across a single large accession's own download lifetime -- e.g.
    once when only a small 8bp I1 file had finished downloading so far
    (correctly reporting "could not classify roles" for that 1-file
    snapshot, since classify_output_files() requires 2-3 files to
    attempt classification at all), and again later once the run's
    much larger R1/R2 files have ALSO finished.

    Every result this function returns is tagged skipped=True by
    construction -- callers MUST use this to distinguish "a cached
    result that can safely be REFRESHED with a newer scan" (skipped=True
    -- it can only have come from this same resume-scan path, and may be
    a stale, incomplete snapshot from earlier in the same download) from
    "a real, just-completed LIVE download result from the current
    session" (skipped=False -- must be preserved untouched, never
    overwritten by a resume-scan result). A caller using a plain
    dict.setdefault()-style "only add if the key doesn't already exist"
    merge will appear to work initially but will silently FREEZE on a
    run's FIRST-ever scanned snapshot forever, even after the real
    underlying download completes and more files appear on disk -- this
    exact symptom (a run stuck showing only 1 of its real 2-3 files,
    indefinitely, across page reloads, requiring only a page refresh to
    ever notice something was wrong) is what originally surfaced this
    caveat. The correct merge is: overwrite any existing cached entry
    for accession IF that entry is missing OR itself tagged
    skipped=True; never overwrite an existing entry tagged skipped=False.
    """
    work_root = os.path.join(project_fastq_dir, "_sra_work")
    if not os.path.isdir(work_root):
        return {}

    pending = {}
    for accession in sorted(os.listdir(work_root)):
        work_dir = os.path.join(work_root, accession)
        if not os.path.isdir(work_dir):
            continue
        raw_files = _find_raw_staged_files(accession, project_fastq_dir)
        if not raw_files:
            # Either nothing downloaded yet for this accession (e.g. a
            # prefetch-only, no-FASTQ-output failure like the confirmed
            # 20GB-size-limit case), or it was already finalized/moved
            # out of staging -- either way, nothing pending to surface.
            continue
        classification = classify_output_files(raw_files)
        bam_warning = detect_likely_bam_derived_issue(classification)
        pending[accession] = {
            "success": True,
            "skipped": True,
            "message": (
                f"Found {len(raw_files)} previously downloaded file(s) for {accession} "
                f"already sitting in _sra_work/{accession}/ -- not yet confirmed/finalized."
            ),
            "classification": classification,
            "bam_warning": bam_warning,
        }
    return pending


_MATCH_BASIS_CONFIDENCE = {
    "GEOAccession": "high",
    "Experiment": "high",
    "SampleTitle": "low",
}


def _group_key_priority(row):
    """
    Return an ordered list of (match_basis, group_key) candidates for a
    single SRA metadata row (as returned by sra.lookup_multiple_
    accessions() / sra.lookup_accession()), from strongest to weakest
    same-biological-sample identity signal. Used by
    detect_multi_run_sample_groups() below.
    """
    candidates = []
    geo = (row.get("GEOAccession") or "").strip()
    if geo:
        candidates.append(("GEOAccession", geo))
    experiment = (row.get("Experiment") or "").strip()
    if experiment and experiment != "—":
        candidates.append(("Experiment", experiment))
    title = (row.get("SampleTitle") or "").strip()
    if title:
        candidates.append(("SampleTitle", title))
    return candidates


def detect_multi_run_sample_groups(rows, requested_accessions=None):
    """
    Given SRA metadata rows (the same shape returned by
    sra_manager.lookup_multiple_accessions() / lookup_accession(), once
    PART 1's "Experiment" field is added), detect which runs likely
    belong to the SAME biological sample and should therefore be merged
    as multiple lanes of one sample_name (via finalize_role_
    assignment()'s existing multi-file-per-role lane-numbering support)
    rather than being treated as separate samples.

    Grouping is attempted using the STRONGEST available identity signal
    for each row, in priority order:

      1. GEOAccession (GSM) -- the most reliable cross-experiment
         "same biological sample" signal for GEO-derived datasets.
         Multiple SRX experiments deposited for one GSM sample (e.g. a
         sample resequenced or split across separate SRA experiment
         submissions -- the confirmed real-world pattern in GSE166992)
         all still point back to the same SAMPLE/GSM element.
      2. Experiment (SRX) -- multiple runs registered under ONE
         experiment/library (a single EXPERIMENT_PACKAGE's own
         RUN_SET/RUN can list more than one run) are guaranteed to be
         the same sequenced library, even for non-GEO datasets with no
         GSM at all.
      3. SampleTitle -- a much weaker fallback (free-text, so only an
         EXACT match is used, and even then flagged "low" confidence)
         for datasets with neither a GEO accession nor a shared
         experiment accession.

    A row is claimed by the FIRST (highest-priority) tier it matches
    another row on; it is never double-counted across tiers. A row
    that matches nothing is left as its own singleton -- i.e. today's
    existing 1-run-1-sample behavior, unchanged.

    requested_accessions, if given (e.g. the exact list of Run
    accessions the user actually pasted/uploaded), is used only to flag
    -- via "includes_unrequested_runs" -- when a detected group contains
    a sibling run the user did NOT originally ask for. This is a real,
    expected possibility: NCBI's own efetch response for an
    EXPERIMENT_PACKAGE includes EVERY run in that experiment's RUN_SET
    regardless of which specific run accession was searched for, so a
    sibling run can legitimately surface here even if the user only
    pasted one of the two. This must always be surfaced to the user for
    confirmation, never silently auto-included.

    Returns a dict:
        {
          "groups": [
              {
                  "group_key": str,
                  "match_basis": "GEOAccession" | "Experiment" | "SampleTitle",
                  "confidence": "high" | "low",
                  "runs": [run_accession, ...],   # sorted, stable order
                  "suggested_sample_name": str,
                  "includes_unrequested_runs": bool,
              },
              ...   # only groups with 2+ runs are included
          ],
          "singletons": [run_accession, ...],   # no detected match
        }
    """
    requested_set = set(requested_accessions) if requested_accessions is not None else None

    by_basis = {"GEOAccession": {}, "Experiment": {}, "SampleTitle": {}}
    row_by_run = {}
    for row in rows:
        run = row.get("Run")
        if not run:
            continue
        row_by_run[run] = row
        for basis, key in _group_key_priority(row):
            by_basis[basis].setdefault(key, []).append(run)

    grouped_runs = set()
    groups = []

    for basis in ("GEOAccession", "Experiment", "SampleTitle"):
        for key, runs in by_basis[basis].items():
            # Only runs not already claimed by a HIGHER-priority tier
            # are eligible here, and only 2+ *remaining* runs form a
            # real multi-run group.
            remaining = [r for r in dict.fromkeys(runs) if r not in grouped_runs]
            if len(remaining) < 2:
                continue

            remaining_sorted = sorted(remaining)
            grouped_runs.update(remaining_sorted)

            geo_for_name = row_by_run[remaining_sorted[0]].get("GEOAccession") or ""
            suggested_name = geo_for_name.strip() or key or remaining_sorted[0]

            includes_unrequested = (
                requested_set is not None
                and any(r not in requested_set for r in remaining_sorted)
            )

            groups.append({
                "group_key": key,
                "match_basis": basis,
                "confidence": _MATCH_BASIS_CONFIDENCE[basis],
                "runs": remaining_sorted,
                "suggested_sample_name": suggested_name,
                "includes_unrequested_runs": includes_unrequested,
            })

    singletons = sorted(r for r in row_by_run if r not in grouped_runs)

    return {"groups": groups, "singletons": singletons}


def build_accession_to_sample_name(detected, confirmed_group_names=None,
                                    manual_singleton_names=None,
                                    manual_group_overrides=None):
    """
    Turn detect_multi_run_sample_groups()'s output into the final
    accession -> sample_name mapping to pass into
    download_and_classify_runs_parallel() (and, downstream,
    finalize_role_assignment()) -- AFTER the user has reviewed and
    confirmed (or overridden) the detected grouping in the UI.

    confirmed_group_names: optional {group_key: sample_name} overriding
        a detected group's own suggested_sample_name (e.g. the user
        edited the auto-suggested name in the UI).

    manual_singleton_names: optional {run_accession: sample_name} for
        runs left as singletons -- default is sample_name == run
        accession (today's existing 1:1 behavior); this allows override.

    manual_group_overrides: optional list of run-accession lists the
        user manually grouped together in the UI -- this is the
        "manual override allowing users to indicate that multiple
        FASTQ files represent different lanes of a single sample" from
        the original design notes, for cases auto-detection misses
        entirely (e.g. two runs sharing none of the three identity
        signals above). e.g. [["SRR1", "SRR2"]]. Each inner list's runs
        are merged under one sample_name, taken from
        confirmed_group_names keyed by a synthetic "manual:<i>" key, or
        defaulting to the first (sorted) run's own accession.

    Every run accession from `detected` (groups + singletons) MUST
    appear in the returned mapping exactly once, unless it was moved
    into a manual_group_overrides group instead. This function performs
    a completeness check and raises ValueError if any run ends up
    unmapped -- a silently-dropped accession would mean a real sample's
    data never gets downloaded/finalized at all, which must never
    happen quietly.
    """
    confirmed_group_names = confirmed_group_names or {}
    manual_singleton_names = manual_singleton_names or {}
    manual_group_overrides = manual_group_overrides or []

    mapping = {}
    manually_placed = set()

    for i, run_list in enumerate(manual_group_overrides):
        synthetic_key = f"manual:{i}"
        sample_name = confirmed_group_names.get(synthetic_key) or sorted(run_list)[0]
        for run in run_list:
            mapping[run] = sample_name
            manually_placed.add(run)

    for group in detected.get("groups", []):
        sample_name = confirmed_group_names.get(group["group_key"], group["suggested_sample_name"])
        for run in group["runs"]:
            if run in manually_placed:
                continue
            mapping[run] = sample_name

    for run in detected.get("singletons", []):
        if run in manually_placed:
            continue
        mapping[run] = manual_singleton_names.get(run, run)

    all_known_runs = set(mapping) | manually_placed
    expected_runs = set(detected.get("singletons", []))
    for group in detected.get("groups", []):
        expected_runs.update(group["runs"])
    missing = expected_runs - all_known_runs
    if missing:
        raise ValueError(
            f"The following run accession(s) were not assigned a sample_name and would be "
            f"silently skipped: {sorted(missing)}. This should never happen -- please report it."
        )

    return mapping


def download_and_classify_run(accession, project_fastq_dir, sample_name, threads=4,
                               subprocess_runner=None, force_redownload=False,
                               max_size=DEFAULT_PREFETCH_MAX_SIZE):
    """
    Unless force_redownload=True, this now checks
    is_accession_already_downloaded() FIRST and, if the accession is
    already present (via any of the three detection tiers), skips
    prefetch/fasterq-dump entirely and returns a "skipped" result
    instead of silently re-downloading data that's already there.
    """
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    if not force_redownload:
        existing = is_accession_already_downloaded(accession, project_fastq_dir)
        if existing["already_downloaded"]:
            classification = None
            bam_warning = None
            if existing["source"] in ("raw_staging", "finalized_naming"):
                classification = classify_output_files(existing["existing_files"])
                bam_warning = detect_likely_bam_derived_issue(classification)
            return {
                "success": True,
                "skipped": True,
                "message": (
                    f"Skipped {accession} -- already downloaded. {existing['detail']} "
                    f"Pass force_redownload=True to re-download anyway."
                ),
                "classification": classification,
                "bam_warning": bam_warning,
            }

    work_dir = os.path.join(project_fastq_dir, "_sra_work", accession)
    os.makedirs(work_dir, exist_ok=True)

    prefetch_cmd = build_prefetch_command(accession, work_dir, max_size=max_size)
    result = runner(prefetch_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {"success": False, "skipped": False, "message": f"prefetch failed for {accession}: {result.stderr}", "classification": None, "bam_warning": None}

    # --- "prefetch quiet success, no .sra file" diagnostic fix
    # (2026-09-01) -- see module docstring for full rationale.
    sra_file = _find_prefetched_sra_file(accession, work_dir)
    if sra_file is None:
        captured_output = (result.stdout or "") + (result.stderr or "")
        size_limit_hit = "skipped" in captured_output.lower() and "large" in captured_output.lower()
        if size_limit_hit:
            hint = (
                f"⚠️ This looks like prefetch's own max-download-size limit was hit -- its "
                f"captured output mentions the file being skipped for being too large. This "
                f"should not recur going forward since build_prefetch_command() now always "
                f"passes an explicit --max-size ({DEFAULT_PREFETCH_MAX_SIZE} by default); if "
                f"you're still seeing this, the accession may exceed even that raised limit -- "
                f"check its real size on NCBI and pass a higher max_size if needed."
            )
        else:
            hint = (
                f"This usually means one of: (1) a transient NCBI network issue or rate-limit "
                f"that didn't trigger a non-zero exit code -- retrying often resolves this; "
                f"(2) this accession requires dbGaP/controlled-access authorization (an .ngc "
                f"file) that wasn't provided; or (3) this specific run has been withdrawn or "
                f"moved on NCBI's end."
            )
        return {
            "success": False,
            "skipped": False,
            "message": (
                f"prefetch for {accession} reported success (exit code 0), but no .sra file "
                f"was found afterward at the expected location(s) under {work_dir}. {hint} "
                f"prefetch's own captured output:\n\n{captured_output or '(no output captured)'}"
            ),
            "classification": None, "bam_warning": None,
        }

    dump_cmd = build_fasterq_dump_command(accession, work_dir, work_dir, threads=threads)
    result2 = runner(dump_cmd, capture_output=True, text=True)
    if result2.returncode != 0:
        return {"success": False, "skipped": False, "message": f"fasterq-dump failed for {accession}: {result2.stderr}", "classification": None, "bam_warning": None}

    raw_produced_files = sorted([
        os.path.join(work_dir, f) for f in os.listdir(work_dir)
        if f.startswith(accession) and (f.endswith(".fastq") or f.endswith(".fastq.gz"))
    ])
    if not raw_produced_files:
        return {"success": False, "skipped": False, "message": f"fasterq-dump produced no output files for {accession}.", "classification": None, "bam_warning": None}

    produced_files = []
    for path in raw_produced_files:
        if path.endswith(".fastq"):
            produced_files.append(_gzip_and_remove(path))
        else:
            produced_files.append(path)
    produced_files = sorted(produced_files)

    classification = classify_output_files(produced_files)
    bam_warning = detect_likely_bam_derived_issue(classification)

    # Record this successful raw download in the registry immediately --
    # BEFORE finalization -- so a re-run before finalization still gets
    # detected via the registry (in addition to the raw-staging fallback
    # check above).
    registry = _load_download_registry(project_fastq_dir)
    registry[accession] = {
        "sample_name": sample_name,
        "status": "raw_downloaded",
        "raw_files": produced_files,
        "finalized_files": registry.get(accession, {}).get("finalized_files", []),
        "timestamp": _now_iso(),
    }
    _save_download_registry(project_fastq_dir, registry)

    return {
        "success": True,
        "skipped": False,
        "message": f"Downloaded, compressed, and classified {len(produced_files)} file(s) for {accession}.",
        "classification": classification,
        "bam_warning": bam_warning,
    }


def download_and_classify_runs_parallel(accession_to_sample_name, project_fastq_dir,
                                         max_workers=3, threads_per_run=4,
                                         on_run_complete=None, force_redownload=False,
                                         max_size=DEFAULT_PREFETCH_MAX_SIZE):
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_accession = {
            executor.submit(
                download_and_classify_run, accession, project_fastq_dir, sample_name, threads_per_run,
                None, force_redownload, max_size,
            ): accession
            for accession, sample_name in accession_to_sample_name.items()
        }

        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                result = future.result()
            except Exception as e:
                result = {
                    "success": False,
                    "skipped": False,
                    "message": f"Unexpected error downloading {accession}: {e}",
                    "classification": None,
                    "bam_warning": None,
                }

            results[accession] = result
            if on_run_complete:
                on_run_complete(accession, result)

    return results


def finalize_role_assignment(classification_overrides, project_fastq_dir, sample_name, lane_map=None):
    lane_map = lane_map or {}

    paths_by_role = {}
    for path, role in classification_overrides.items():
        if role not in (ROLE_R1, ROLE_R2):
            continue
        paths_by_role.setdefault(role, []).append(path)

    destinations = {}
    for role, paths in paths_by_role.items():
        if len(paths) == 1:
            path = paths[0]
            ext = ".fastq.gz" if path.endswith(".gz") else ".fastq"
            dest = os.path.join(project_fastq_dir, f"{sample_name}_S1_{role}_001{ext}")
            shutil.move(path, dest)
            destinations[role] = dest
        else:
            sorted_paths = sorted(paths, key=lambda p: (lane_map.get(p, 0), p))
            dest_list = []
            for i, path in enumerate(sorted_paths, start=1):
                lane_num = lane_map.get(path, i)
                ext = ".fastq.gz" if path.endswith(".gz") else ".fastq"
                dest = os.path.join(project_fastq_dir, f"{sample_name}_S1_L{lane_num:03d}_{role}_001{ext}")
                shutil.move(path, dest)
                dest_list.append(dest)
            destinations[role] = dest_list

    return destinations


# ---------------------------------------------------------------------------
# Original-format BAM recovery
# ---------------------------------------------------------------------------

def bamtofastq_available():
    return shutil.which("bamtofastq") is not None


BAMTOFASTQ_CHEMISTRY_FLAGS = {
    "auto": {
        "label": "Automatic (Cell Ranger 1.2+ / Long Ranger 2.1+ -- recommended default)",
        "flag": None,
        "explanation": (
            "No special flag needed -- the BAM's own header fields describe which "
            "chemistry/read-layout was used, and bamtofastq reads this automatically. "
            "Correct for the vast majority of real-world 10x BAMs. ⚠️ If this BAM "
            "predates Cell Ranger 1.2 (e.g. an older public dataset), 'auto' will NOT "
            "correctly detect its chemistry -- bamtofastq may report success while "
            "writing zero output files. If that happens, try one of the specific "
            "older-version options below instead."
        ),
    },
    "cr11": {
        "label": "Cell Ranger 1.0-1.1 (--cr11)",
        "flag": "--cr11",
        "explanation": (
            "Required specifically for BAMs produced by Cell Ranger 1.0 or 1.1 -- "
            "these predate the self-describing BAM header fields 'auto' relies on. "
            "Confirmed relevant for the ORIGINAL Kang et al. 2018 dataset (GSE96583). "
            "IMPORTANT: this chemistry's bamtofastq output is a 4-file-per-lane split "
            "(cDNA + barcode + UMI + sample index, in SEPARATE files) rather than a "
            "combined barcode+UMI read -- this module automatically detects and "
            "resolves that split (see resolve_bamtofastq_v1_split()), no manual action "
            "needed."
        ),
    },
    "gemcode": {
        "label": "Original GemCode data (Long Ranger 1.0-1.3) (--gemcode)",
        "flag": "--gemcode",
        "explanation": "Required for BAMs from the original GemCode instrument/chemistry, predating the Chromium platform entirely.",
    },
    "lr20": {
        "label": "Long Ranger 2.0 (--lr20)",
        "flag": "--lr20",
        "explanation": "Required specifically for BAMs produced by Long Ranger 2.0 (linked-read/genomic data, not gene expression).",
    },
}
DEFAULT_BAMTOFASTQ_CHEMISTRY_FLAG = "auto"


def find_original_format_bam_url(accession):
    uid_list, error = sra._esearch_sra(accession)
    if error:
        return None, f"Could not look up '{accession}': {error}"
    if not uid_list:
        return None, f"No SRA record found for '{accession}'."

    root, error = sra._efetch_sra_full_xml(uid_list)
    if error:
        return None, f"Could not fetch metadata for '{accession}': {error}"

    for run_el in root.iter("RUN"):
        run_accession = run_el.get("accession", "")
        if run_accession and run_accession != accession:
            continue
        for file_el in run_el.iter("SRAFile"):
            supertype = (file_el.get("supertype") or "").strip().lower()
            url = file_el.get("url")
            if supertype == "original" and url:
                filename = file_el.get("filename", os.path.basename(url))
                return url, f"Found an original-format file for '{accession}': {filename}"

    return None, (
        f"No directly-downloadable original-format file was found for '{accession}' in its "
        f"own SRA metadata. This most likely means NCBI has moved this run's original "
        f"submission to cloud-only storage (AWS/GCP), which requires your OWN cloud-billing "
        f"account to retrieve -- this pipeline does not automate that path. Check this "
        f"accession's own SRA 'Data access' page directly to confirm."
    )


def download_original_bam(accession, dest_dir, url=None, chunk_size=1024 * 1024, progress_callback=None):
    if url is None:
        url, message = find_original_format_bam_url(accession)
        if url is None:
            return False, None, message

    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, f"{accession}.original.bam")

    if progress_callback:
        progress_callback(f"Downloading original-format BAM for {accession}...")

    try:
        with urllib.request.urlopen(url, timeout=60) as response, open(dest_path, "wb") as out_file:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        if os.path.exists(dest_path):
            os.remove(dest_path)
        return False, None, f"Failed to download original-format BAM for {accession}: {e}"

    return True, dest_path, f"Downloaded original-format BAM for {accession} to {dest_path}."


def run_bamtofastq(bam_path, dest_dir, chemistry_key=DEFAULT_BAMTOFASTQ_CHEMISTRY_FLAG,
                    subprocess_runner=None, timeout=3600):
    """
    Run 10x Genomics' own `bamtofastq` tool. ALWAYS captures and returns
    bamtofastq's own stdout+stderr in the message, not just on a
    non-zero return code.
    """
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    if chemistry_key not in BAMTOFASTQ_CHEMISTRY_FLAGS:
        return False, f"Unknown bamtofastq chemistry option: {chemistry_key!r}"

    if not os.path.isfile(bam_path):
        return False, f"BAM file not found: {bam_path}"

    os.makedirs(os.path.dirname(dest_dir) or ".", exist_ok=True)

    cmd = ["bamtofastq"]
    flag = BAMTOFASTQ_CHEMISTRY_FLAGS[chemistry_key]["flag"]
    if flag:
        cmd.append(flag)
    cmd.extend([bam_path, dest_dir])

    try:
        result = runner(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess_module.TimeoutExpired:
        return False, f"bamtofastq timed out after {timeout} seconds."

    captured_output = (result.stdout or "") + (result.stderr or "")

    if result.returncode != 0:
        return False, f"bamtofastq failed (exit code {result.returncode}): {captured_output}"

    return True, f"bamtofastq exited successfully (exit code 0). Tool output:\n{captured_output}"


def find_fastq_files_under(root_dir):
    found = []
    for dirpath, _dirnames, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.endswith(".fastq") or fname.endswith(".fastq.gz"):
                found.append(os.path.join(dirpath, fname))
    return sorted(found)


def find_fastq_files_under_with_settle_retry(root_dir, max_attempts=6, delay_seconds=10,
                                              sleep_fn=None, progress_callback=None):
    """
    Retry-with-backoff wrapper around find_fastq_files_under() -- for
    network/cluster scratch filesystem write-visibility lag after
    bamtofastq exits.
    """
    import time as time_module
    sleep = sleep_fn or time_module.sleep

    found = find_fastq_files_under(root_dir)
    if found:
        return found

    for attempt in range(1, max_attempts):
        if progress_callback:
            progress_callback(
                f"No output files visible yet under {root_dir} -- this can be normal on "
                f"network/scratch storage for a large file even after the conversion tool "
                f"itself has exited. Waiting and re-checking (attempt {attempt} of "
                f"{max_attempts - 1})..."
            )
        sleep(delay_seconds)
        found = find_fastq_files_under(root_dir)
        if found:
            return found

    return found


# ---------------------------------------------------------------------------
# v1-chemistry bamtofastq output resolution (2026-08-24)
# ---------------------------------------------------------------------------

_BAMTOFASTQ_ROLE_PATTERN = re.compile(r"_(R1|R2|R3|I1)_")


def detect_bamtofastq_v1_split(fastq_paths):
    if not any(_BAMTOFASTQ_ROLE_PATTERN.search(os.path.basename(p)) and "_R3_" in os.path.basename(p) for p in fastq_paths):
        return {}

    groups = {}
    for path in fastq_paths:
        basename = os.path.basename(path)
        match = _BAMTOFASTQ_ROLE_PATTERN.search(basename)
        if not match:
            continue
        role = match.group(1)
        lane_key = _BAMTOFASTQ_ROLE_PATTERN.sub("_{role}_", basename, count=1)
        groups.setdefault(lane_key, {})[role] = path

    complete_groups = {
        lane_key: roles for lane_key, roles in groups.items()
        if "R1" in roles and "R2" in roles and "R3" in roles
    }
    return complete_groups


def _read_fastq_record(file_handle):
    header = file_handle.readline()
    if not header:
        return None
    seq = file_handle.readline()
    plus_line = file_handle.readline()
    qual = file_handle.readline()
    if not (seq and plus_line and qual):
        raise ValueError(f"Truncated FASTQ record encountered (incomplete 4-line group) after header: {header.strip()!r}")
    return header.rstrip("\n"), seq.rstrip("\n"), plus_line.rstrip("\n"), qual.rstrip("\n")


def _read_id_without_mate_suffix(header):
    read_id = header.split()[0] if header else ""
    if len(read_id) >= 2 and read_id[-2] == "/" and read_id[-1] in "123":
        read_id = read_id[:-2]
    return read_id


def concatenate_barcode_umi_reads(barcode_path, umi_path, dest_path, verify_read_identity=True):
    if not os.path.isfile(barcode_path):
        return False, f"Barcode file not found: {barcode_path}", 0
    if not os.path.isfile(umi_path):
        return False, f"UMI file not found: {umi_path}", 0

    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)

    barcode_opener = gzip.open if barcode_path.endswith(".gz") else open
    umi_opener = gzip.open if umi_path.endswith(".gz") else open
    dest_opener = gzip.open if dest_path.endswith(".gz") else open

    n_combined = 0
    try:
        with barcode_opener(barcode_path, "rt") as bc_f, \
             umi_opener(umi_path, "rt") as umi_f, \
             dest_opener(dest_path, "wt") as out_f:

            while True:
                bc_record = _read_fastq_record(bc_f)
                umi_record = _read_fastq_record(umi_f)

                if bc_record is None and umi_record is None:
                    break
                if bc_record is None or umi_record is None:
                    if os.path.exists(dest_path):
                        os.remove(dest_path)
                    return False, (
                        f"Barcode file and UMI file have a MISMATCHED read count "
                        f"(one file ended before the other, after {n_combined} matched "
                        f"read(s) already combined) -- these two files do not appear to "
                        f"correspond to the same original reads. Aborting rather than "
                        f"producing a partially-combined, unreliable output file."
                    ), n_combined

                bc_header, bc_seq, _bc_plus, bc_qual = bc_record
                umi_header, umi_seq, _umi_plus, umi_qual = umi_record

                if verify_read_identity:
                    bc_id = _read_id_without_mate_suffix(bc_header)
                    umi_id = _read_id_without_mate_suffix(umi_header)
                    if bc_id != umi_id:
                        if os.path.exists(dest_path):
                            os.remove(dest_path)
                        return False, (
                            f"Barcode read and UMI read at position {n_combined + 1} do NOT "
                            f"share the same read ID ('{bc_id}' vs. '{umi_id}') -- these two "
                            f"files do not appear to be correctly paired/ordered. Aborting "
                            f"rather than producing a potentially mis-paired, unreliable "
                            f"output file."
                        ), n_combined

                combined_seq = bc_seq + umi_seq
                combined_qual = bc_qual + umi_qual
                out_f.write(f"{bc_header}\n{combined_seq}\n+\n{combined_qual}\n")
                n_combined += 1

    except ValueError as e:
        if os.path.exists(dest_path):
            os.remove(dest_path)
        return False, f"Malformed FASTQ encountered while combining barcode+UMI reads: {e}", n_combined

    return True, f"Combined {n_combined} barcode+UMI read pair(s) into {dest_path}.", n_combined


def resolve_bamtofastq_v1_split(fastq_paths, work_dir, verify_read_identity=True):
    v1_groups = detect_bamtofastq_v1_split(fastq_paths)
    if not v1_groups:
        return fastq_paths, {}, None

    os.makedirs(work_dir, exist_ok=True)

    resolved_paths = []
    lane_map = {}
    failures = []

    for lane_index, (lane_key, roles) in enumerate(sorted(v1_groups.items()), start=1):
        barcode_path = roles["R2"]
        umi_path = roles["R3"]
        cdna_path = roles["R1"]

        combined_dest = os.path.join(work_dir, f"lane{lane_index:03d}_combined_barcode_umi.fastq.gz")
        success, message, _n_reads = concatenate_barcode_umi_reads(
            barcode_path, umi_path, combined_dest, verify_read_identity=verify_read_identity,
        )
        if not success:
            failures.append(f"Lane {lane_index} ({os.path.basename(barcode_path)} + {os.path.basename(umi_path)}): {message}")
            continue

        resolved_paths.append(combined_dest)
        resolved_paths.append(cdna_path)
        lane_map[combined_dest] = lane_index
        lane_map[cdna_path] = lane_index

    warning = None
    if failures:
        warning = (
            f"⚠️ {len(failures)} of {len(v1_groups)} lane(s) could not be resolved and were "
            f"excluded from the recovered result:\n" + "\n".join(f"  - {f}" for f in failures)
        )

    if not resolved_paths:
        return fastq_paths, {}, (warning or "⚠️ Could not resolve any lane of the detected v1-chemistry file split.")

    return resolved_paths, lane_map, warning


def classify_resolved_files(resolved_paths, lane_map):
    if not lane_map:
        return classify_output_files(resolved_paths)

    paths_by_lane = {}
    for path in resolved_paths:
        lane = lane_map.get(path, 1)
        paths_by_lane.setdefault(lane, []).append(path)

    merged = {}
    for lane, paths in paths_by_lane.items():
        merged.update(classify_output_files(paths))
    return merged


def run_full_bam_recovery_pipeline(accession, project_fastq_dir, sample_name,
                                    chemistry_key=DEFAULT_BAMTOFASTQ_CHEMISTRY_FLAG,
                                    progress_callback=None):
    work_dir = os.path.join(project_fastq_dir, "_sra_bam_work", accession)
    os.makedirs(work_dir, exist_ok=True)

    success, bam_path, message = download_original_bam(accession, work_dir, progress_callback=progress_callback)
    if not success:
        return {"success": False, "message": message, "classification": None, "bam_warning": None, "lane_map": {}}

    if progress_callback:
        progress_callback(f"Converting BAM to FASTQ for {accession} (bamtofastq)...")

    fastq_dest_dir = os.path.join(work_dir, "fastq_output")
    success, bamtofastq_message = run_bamtofastq(bam_path, fastq_dest_dir, chemistry_key=chemistry_key)
    if not success:
        return {"success": False, "message": bamtofastq_message, "classification": None, "bam_warning": None, "lane_map": {}}

    produced_files = find_fastq_files_under_with_settle_retry(fastq_dest_dir, progress_callback=progress_callback)
    if not produced_files:
        return {
            "success": False,
            "message": (
                f"bamtofastq reported success, but no FASTQ files became visible under "
                f"{fastq_dest_dir} even after waiting. This could mean the selected chemistry "
                f"option ('{chemistry_key}') does not match this BAM's true origin (if this BAM "
                f"predates Cell Ranger 1.2, try 'Cell Ranger 1.0-1.1 (--cr11)' instead of "
                f"'Automatic') -- OR it could mean this specific storage location is taking "
                f"unusually long to make large written files visible; if you know files should "
                f"be present, wait a few more minutes and check the destination directory "
                f"directly before retrying. bamtofastq's own captured output:\n\n{bamtofastq_message}"
            ),
            "classification": None, "bam_warning": None, "lane_map": {},
        }

    if progress_callback:
        progress_callback(f"Checking for a legacy multi-file barcode/UMI split ({accession})...")
    resolve_work_dir = os.path.join(work_dir, "resolved")
    resolved_paths, lane_map, resolve_warning = resolve_bamtofastq_v1_split(produced_files, resolve_work_dir)
    if resolve_warning and progress_callback:
        progress_callback(resolve_warning)

    classification = classify_resolved_files(resolved_paths, lane_map)
    n_lanes_note = f" across {len(set(lane_map.values()))} lane(s)" if lane_map else ""
    message = f"Recovered {len(resolved_paths)} FASTQ file(s) for {accession}{n_lanes_note} via original-format BAM + bamtofastq."
    if resolve_warning:
        message += f" ({resolve_warning})"

    return {
        "success": True,
        "message": message,
        "classification": classification,
        "bam_warning": None,
        "lane_map": lane_map,
    }


def run_full_bam_recovery_pipeline_parallel(accession_to_sample_name, project_fastq_dir,
                                             chemistry_key=DEFAULT_BAMTOFASTQ_CHEMISTRY_FLAG,
                                             max_workers=2, on_run_complete=None):
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_accession = {
            executor.submit(
                run_full_bam_recovery_pipeline, accession, project_fastq_dir, sample_name, chemistry_key,
            ): accession
            for accession, sample_name in accession_to_sample_name.items()
        }

        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                result = future.result()
            except Exception as e:
                result = {
                    "success": False,
                    "message": f"Unexpected error recovering {accession} from original-format BAM: {e}",
                    "classification": None, "bam_warning": None, "lane_map": {},
                }

            results[accession] = result
            if on_run_complete:
                on_run_complete(accession, result)

    return results


def save_uploaded_bam(uploaded_file, dest_dir, sample_name):
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, f"{sample_name}.uploaded.bam")
    with open(dest_path, "wb") as f:
        f.write(uploaded_file.read())
    return dest_path


def run_bam_recovery_from_uploaded_file(uploaded_bam_path, project_fastq_dir, sample_name,
                                         chemistry_key=DEFAULT_BAMTOFASTQ_CHEMISTRY_FLAG,
                                         progress_callback=None):
    if not os.path.isfile(uploaded_bam_path):
        return {
            "success": False,
            "message": f"Uploaded BAM file not found at {uploaded_bam_path}.",
            "classification": None, "bam_warning": None, "lane_map": {},
        }

    work_dir = os.path.join(project_fastq_dir, "_sra_bam_work", sample_name)
    os.makedirs(work_dir, exist_ok=True)

    if progress_callback:
        progress_callback(f"Converting uploaded BAM to FASTQ for {sample_name} (bamtofastq)...")

    fastq_dest_dir = os.path.join(work_dir, "fastq_output")
    success, bamtofastq_message = run_bamtofastq(uploaded_bam_path, fastq_dest_dir, chemistry_key=chemistry_key)
    if not success:
        return {"success": False, "message": bamtofastq_message, "classification": None, "bam_warning": None, "lane_map": {}}

    produced_files = find_fastq_files_under_with_settle_retry(fastq_dest_dir, progress_callback=progress_callback)
    if not produced_files:
        return {
            "success": False,
            "message": (
                f"bamtofastq reported success, but no FASTQ files became visible under "
                f"{fastq_dest_dir} even after waiting. This could mean the selected chemistry "
                f"option ('{chemistry_key}') does not match this BAM's true origin (if this BAM "
                f"predates Cell Ranger 1.2, try 'Cell Ranger 1.0-1.1 (--cr11)' instead of "
                f"'Automatic') -- OR it could mean this specific storage location is taking "
                f"unusually long to make large written files visible; if you know files should "
                f"be present, wait a few more minutes and check the destination directory "
                f"directly before retrying. bamtofastq's own captured output:\n\n{bamtofastq_message}"
            ),
            "classification": None, "bam_warning": None, "lane_map": {},
        }

    if progress_callback:
        progress_callback(f"Checking for a legacy multi-file barcode/UMI split ({sample_name})...")
    resolve_work_dir = os.path.join(work_dir, "resolved")
    resolved_paths, lane_map, resolve_warning = resolve_bamtofastq_v1_split(produced_files, resolve_work_dir)
    if resolve_warning and progress_callback:
        progress_callback(resolve_warning)

    classification = classify_resolved_files(resolved_paths, lane_map)
    n_lanes_note = f" across {len(set(lane_map.values()))} lane(s)" if lane_map else ""
    message = f"Recovered {len(resolved_paths)} FASTQ file(s) for {sample_name}{n_lanes_note} from the uploaded BAM."
    if resolve_warning:
        message += f" ({resolve_warning})"

    return {
        "success": True,
        "message": message,
        "classification": classification,
        "bam_warning": None,
        "lane_map": lane_map,
    }
