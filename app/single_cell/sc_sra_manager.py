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
chemistry flag mismatches the BAM's actual origin -- confirmed
motivating case: a Cell Ranger 1.1 BAM converted with the "auto"
chemistry option instead of the required "--cr11" flag, since CR 1.1
BAMs predate the self-describing header fields "auto" relies on),
callers had ZERO visibility into what bamtofastq itself actually
printed, since only a generic "completed successfully" message was ever
returned in the success case.

Fixed by having run_bamtofastq() ALWAYS capture and return bamtofastq's
own stdout+stderr (not just on failure), and by having
run_full_bam_recovery_pipeline()/run_bam_recovery_from_uploaded_file()
explicitly check for this "success but zero files" case and surface
that captured tool output directly in the returned error message --
rather than the previous generic "bamtofastq completed but no FASTQ
files were found" message, which gave no actionable information for
diagnosing WHY. This is purely additive -- the normal, expected-output
case is completely unaffected; only the previously-silent diagnostic
gap for this specific anomalous outcome is closed.
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
    roles = {info["role"] for info in classification.values()}
    if roles == {ROLE_UNKNOWN}:
        return (
            "⚠️ Could not confidently identify cell-barcode/UMI vs. cDNA reads "
            "from this run's downloaded file(s). This can happen when an "
            "accession was originally deposited as a 10x Genomics BAM file "
            "rather than plain FASTQ, or when only the cDNA read was ever "
            "deposited to SRA at all (a known real gap for some older, "
            "v1-chemistry-era 10x depositions). Check this accession's own "
            "SRA page, or try this module's find_original_format_bam_url() / "
            "run_full_bam_recovery_pipeline() to see whether the ORIGINAL "
            "Cell Ranger BAM (which retains barcode/UMI tags a plain FASTQ "
            "extraction does not) is directly recoverable for this run."
        )
    return None


def build_prefetch_command(accession, output_dir):
    return ["prefetch", accession, "--output-directory", output_dir]


def build_fasterq_dump_command(accession, sra_file_dir, output_dir, threads=4):
    return [
        "fasterq-dump", "--split-files", "--include-technical",
        "--threads", str(threads), "--outdir", output_dir,
        os.path.join(sra_file_dir, accession),
    ]


def _gzip_and_remove(src_path):
    gz_path = src_path + ".gz"
    with open(src_path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(src_path)
    return gz_path


def download_and_classify_run(accession, project_fastq_dir, sample_name, threads=4, subprocess_runner=None):
    import subprocess as subprocess_module
    runner = subprocess_runner or subprocess_module.run

    work_dir = os.path.join(project_fastq_dir, "_sra_work", accession)
    os.makedirs(work_dir, exist_ok=True)

    prefetch_cmd = build_prefetch_command(accession, work_dir)
    result = runner(prefetch_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {"success": False, "message": f"prefetch failed for {accession}: {result.stderr}", "classification": None, "bam_warning": None}

    dump_cmd = build_fasterq_dump_command(accession, work_dir, work_dir, threads=threads)
    result2 = runner(dump_cmd, capture_output=True, text=True)
    if result2.returncode != 0:
        return {"success": False, "message": f"fasterq-dump failed for {accession}: {result2.stderr}", "classification": None, "bam_warning": None}

    raw_produced_files = sorted([
        os.path.join(work_dir, f) for f in os.listdir(work_dir)
        if f.startswith(accession) and (f.endswith(".fastq") or f.endswith(".fastq.gz"))
    ])
    if not raw_produced_files:
        return {"success": False, "message": f"fasterq-dump produced no output files for {accession}.", "classification": None, "bam_warning": None}

    produced_files = []
    for path in raw_produced_files:
        if path.endswith(".fastq"):
            produced_files.append(_gzip_and_remove(path))
        else:
            produced_files.append(path)
    produced_files = sorted(produced_files)

    classification = classify_output_files(produced_files)
    bam_warning = detect_likely_bam_derived_issue(classification)
    return {
        "success": True,
        "message": f"Downloaded, compressed, and classified {len(produced_files)} file(s) for {accession}.",
        "classification": classification,
        "bam_warning": bam_warning,
    }


def download_and_classify_runs_parallel(accession_to_sample_name, project_fastq_dir,
                                         max_workers=3, threads_per_run=4,
                                         on_run_complete=None):
    results = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_accession = {
            executor.submit(
                download_and_classify_run, accession, project_fastq_dir, sample_name, threads_per_run,
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
    Run 10x Genomics' own `bamtofastq` tool.

    --- "Quiet success, zero files" diagnostic fix (2026-08-24) ---
    ALWAYS captures and returns bamtofastq's own stdout+stderr in the
    message, not just on a non-zero return code -- see this module's own
    docstring for the real motivating bug this fixes (a chemistry-flag
    mismatch, e.g. "auto" against a pre-1.2 Cell Ranger BAM, can cause
    bamtofastq to exit 0 while writing zero files, and the caller
    previously had no way to see WHY since only failure-path output was
    ever captured).

    Returns (success: bool, message: str) -- success reflects ONLY
    bamtofastq's own return code (0 = True); it does NOT check whether
    any files were actually produced -- that check happens one level up,
    in run_full_bam_recovery_pipeline()/run_bam_recovery_from_uploaded_file(),
    which have access to dest_dir's actual contents AND this function's
    own captured tool output to build an informative combined message if
    zero files turn up despite a successful (exit 0) run.
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

    # Always capture tool output now, regardless of return code -- see
    # this function's own docstring for why the SUCCESS case specifically
    # needed this too.
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
    Retry-with-backoff wrapper around find_fastq_files_under() -- added
    2026-08-24 after a REAL reported bug: a completed bamtofastq run
    (confirmed exit code 0, and a real, large ~24GB original BAM) was
    followed by a single, immediate find_fastq_files_under() check that
    found ZERO files -- while a directory listing checked manually
    minutes later showed all 16 expected output files present, with
    timestamps spread across roughly 5-10 MINUTES after the subprocess
    call returned.

    Root cause (network/cluster scratch filesystem write-visibility
    lag, not a chemistry-flag or tool-configuration issue): the
    project's own data directory is confirmed to live under
    /disk/bioscratch/... -- a network-mounted (NFS/Lustre-style)
    scratch filesystem, per the user's own real directory browser
    screenshot. On this class of storage, a subprocess reporting "I
    have exited" does NOT guarantee that a DIFFERENT process's (or
    even the SAME process's) very next directory listing will reflect
    every file that process just wrote -- there can be a real,
    sometimes multi-minute lag between a large file's local write
    completing and that write becoming visible via a fresh os.walk()/
    os.listdir() call, particularly for large output files still being
    flushed/synced to network storage. bamtofastq's own output for a
    large BAM is written as MULTIPLE SEQUENTIAL CHUNKED FILES over an
    extended period (confirmed directly: the user's own real recovered
    files show creation timestamps spread across 5-10 minutes for a
    single accession), making a bare, single, immediate check
    especially fragile for exactly this tool/workload.

    This function retries find_fastq_files_under() up to max_attempts
    times, sleeping delay_seconds between attempts, and returns as soon
    as ANY files are found (it does not try to guess when writing is
    "fully complete" beyond that -- see the module docstring note
    below on why this is a deliberately simple, bounded heuon-fix
    rather than a more complex "wait until file count stops changing"
    stability check).

    sleep_fn: injectable for testing (defaults to time.sleep) -- so a
        test can pass a no-op or instrumented sleep function without
        this function actually blocking for real wall-clock time.
    progress_callback: if given, called with a plain-language status
        message before each retry (after the FIRST, immediate check
        already came back empty) -- so a caller with a live progress
        area (e.g. this module's own Step 1 UI) can show the user
        that this specific, expected wait is happening, rather than
        the page appearing to hang with no explanation for what could
        be several minutes.

    Returns the same sorted list of file paths find_fastq_files_under()
    itself returns -- an empty list if genuinely no files are found
    after exhausting all attempts.
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
    """
    End-to-end orchestration of the full original-format BAM recovery
    path for ONE accession.

    --- "Quiet success, zero files" diagnostic fix (2026-08-24) ---
    When bamtofastq itself reports success (exit 0) but produces ZERO
    FASTQ files, this now includes bamtofastq's OWN captured tool output
    (from run_bamtofastq()'s own message) directly in the returned error
    message, along with an explicit hint to check whether the selected
    chemistry_key actually matches this BAM's true origin -- rather than
    the previous generic "no files found" message with no diagnostic
    detail at all.
    """
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

    # --- Network/scratch-filesystem write-visibility settle-retry
    # (2026-08-24) -- see find_fastq_files_under_with_settle_retry()'s
    # own docstring for the full rationale: a bare, single,
    # immediate check here was confirmed, via a real reported case, to
    # return empty even though bamtofastq's own (large, multi-chunk)
    # output files were still in the process of becoming visible on
    # network/cluster scratch storage -- NOT a chemistry-flag issue,
    # despite that having been an initially reasonable hypothesis
    # before this was confirmed with real directory-listing timestamps.
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
    """
    Recover usable FASTQ from a BAM file the USER ALREADY HAS on hand
    (uploaded, or resolved via the server-directory browser) -- the
    direct complement to run_full_bam_recovery_pipeline() above.

    --- "Quiet success, zero files" diagnostic fix (2026-08-24) ---
    Same fix as run_full_bam_recovery_pipeline()'s own identical
    section above -- see that function's docstring for the full
    rationale (this is the exact code path that surfaced the real,
    motivating bug: a Cell Ranger 1.1 BAM converted with "auto"
    chemistry, producing a quiet, uninformative "no files found"
    failure).
    """
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

    # --- Network/scratch-filesystem write-visibility settle-retry
    # (2026-08-24) -- see find_fastq_files_under_with_settle_retry()'s
    # own docstring for the full rationale, and
    # run_full_bam_recovery_pipeline()'s own identical block above for
    # the real, confirmed motivating case (a genuine ~24GB BAM whose
    # bamtofastq output only became fully visible on this project's
    # network/cluster scratch storage several minutes after the
    # subprocess itself had already exited).
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
