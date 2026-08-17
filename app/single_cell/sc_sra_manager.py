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
here directly via `import sra_manager as sra` -- see
_reused_from_sra_manager below for exactly which functions.

What's different, and why this module exists:

1. FILE ROLE IS AMBIGUOUS AFTER DOWNLOAD. In bulk RNA-seq, fasterq-dump's
   _1/_2 output reliably corresponds to R1/R2. In single-cell, this is
   NOT guaranteed -- a real user-reported case (ncbi/sra-tools GitHub
   issue #418) needed to manually inspect read lengths and cross-
   reference SRA's own metadata to figure out which of three output
   files (_1/_2/_3) was the sample index (I1), which was the cell
   barcode+UMI (R1), and which was the cDNA read (R2) -- there is no
   dependable ordering convention. This module's classify_output_files()
   automates that same length-based inspection (reusing
   chemistry_manager.dominant_read_length(), the same technique
   detect_chemistry() itself uses), but ALWAYS surfaces its guess for
   user confirmation via the workspace UI -- exactly like chemistry
   auto-detection itself -- rather than silently trusting it.

2. SOME DEPOSITIONS ARE BAM, NOT FASTQ, AND NEED A DIFFERENT TOOL
   ENTIRELY. 10x Genomics BAMs (Cell Ranger/Space Ranger/Long Ranger
   output) store barcode/UMI information in special BAM tags that
   ordinary fasterq-dump does not understand -- reconstructing correct
   FASTQs from these requires 10x's own `bamtofastq` tool, confirmed
   against that tool's own README, which explicitly warns: "If you are
   downloading BAM files from SRA, be sure to get the 'Original format'
   BAM file with Type = TenX... Only those BAM files have all the
   original tags present that are required." This module does NOT
   attempt to detect "was this originally a TenX BAM" from esearch/
   efetch metadata alone -- that information isn't reliably queryable
   that way -- so instead, detect_likely_bam_derived_issue() applies a
   POST-HOC heuristic after a fasterq-dump attempt: if the download
   produces an unexpected file count/structure (e.g. only 1 output file
   where 2-3 were expected, or files that don't classify into any
   sensible R1/R2/I1 pattern at all), this is flagged to the user as a
   likely sign of a non-standard/BAM-derived submission needing manual
   investigation on that accession's own SRA page -- NOT silently
   treated as a successful, usable download. Automating bamtofastq
   itself is intentionally out of scope for this pass (see
   environment.yml's own note on this) -- adding it is a reasonable
   future addition if this heuristic turns out to trigger often in
   practice.

3. DOWNLOAD EXECUTION IS IMPLEMENTED HERE, NOT REUSED FROM sra_manager.
   This module calls prefetch/fasterq-dump directly (see
   download_and_classify_run() below) rather than assuming a specific
   download-orchestration function exists in sra_manager.py -- avoiding
   guessing at an unconfirmed function signature in that module. The
   validation/lookup functions THIS module reuses from sra_manager are
   listed explicitly below; nothing else about sra_manager's internals
   is assumed.
"""
import os
import shutil
import subprocess

import sra_manager as sra  # reused directly -- see docstring
import chemistry_manager as chem

# Functions reused AS-IS from the Bulk RNA-Seq pipeline's sra_manager.py
# (confirmed interface, from advanced_mode_workspace.py's own usage):
#   sra.tools_available() -> (prefetch_ok: bool, fasterq_ok: bool)
#   sra.lookup_accession(term) -> (success, rows, message)
#   sra.lookup_multiple_accessions(accessions) -> (success, rows, message, not_found)
#   sra.is_rna_seq(row) -> bool
#   sra.build_metadata_dataframe(rows, selected_runs=None) -> dict
#   sra.has_any_descriptive_metadata(rows, selected_runs=None) -> bool
#   sra._split_accession_list_text(text) -> list[str]
#   sra.parse_accessions_from_file(uploaded_file) -> (accessions, error)

ROLE_R1 = "R1"   # cell barcode + UMI (short, fixed length)
ROLE_R2 = "R2"   # cDNA (longer, variable-ish)
ROLE_I1 = "I1"   # sample index read (very short, typically 8-10bp) -- not used downstream
ROLE_UNKNOWN = "UNKNOWN"

# ---------------------------------------------------------------------------
# Assay-type classification (bulk vs. single-cell vs. spatial) from SRA
# metadata -- NOT from LibraryStrategy
# ---------------------------------------------------------------------------
# CONFIRMED (2026-08-17): SRA's LibraryStrategy field is "RNA-Seq" (or
# occasionally "OTHER") for bulk, single-cell, AND spatial transcriptomics
# alike -- there is no distinct enum value that separates them in SRA's own
# schema. This isn't a gap in this pipeline's own is_rna_seq() check
# specifically; it's a documented, real limitation of the metadata schema
# itself (SenNet's own RNA-seq metadata guide states plainly: "The RNAseq
# assay itself is the same for bulk and single-cell templates"). Real tools
# that DO solve this (e.g. cellgeni/fetch10xmeta, Bioconductor's GEOquery)
# don't use LibraryStrategy for this either -- they parse free-text fields
# (experiment/study/sample TITLE or description) for platform-specific
# keywords ("10x", "Chromium", "Visium", etc.), since LibrarySource/
# LibrarySelection are typically identical (TRANSCRIPTOMIC / cDNA) across
# all three assay types too.
#
# classify_assay_type_from_metadata() below applies that same keyword-
# matching approach. IMPORTANT CAVEAT: this pipeline has only directly
# confirmed a SUBSET of the fields present in sra_manager.py's row dicts
# (Run, ScientificName, LibraryStrategy, size_MB, from confirmed real
# usage) -- the exact set of title/description fields available (if any)
# has not been independently verified against sra_manager.py's real esearch/
# efetch parsing. This function therefore checks EVERY plausible field name
# a title/description might be stored under, defensively, and returns
# "insufficient_info" honestly (rather than fabricating confidence) if none
# of those keys are present or none contain a recognizable keyword. If a
# real accession lookup shows title text is available under a field name
# not in _POSSIBLE_TEXT_FIELDS below, add it there.
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
    """
    Best-effort classification of an SRA row's likely assay type (bulk vs.
    single-cell vs. spatial) from free-text title/description fields --
    see this module's header comment for why LibraryStrategy itself cannot
    do this. Returns a dict:
        {
            "assay_type": ASSAY_SINGLE_CELL | ASSAY_SPATIAL | ASSAY_BULK | ASSAY_UNKNOWN,
            "matched_keyword": str or None,
            "matched_field": str or None,
            "message": str,
        }
    Checks SPATIAL keywords first, then SINGLE_CELL, then BULK-hint --
    spatial methods are checked first since some (e.g. Visium) are
    themselves droplet/bead-based and could otherwise loosely overlap with
    single-cell terminology in a naive check. Returns ASSAY_UNKNOWN
    (explicitly, not a guess) if no text field with a recognizable keyword
    is found at all -- this is expected to happen often, depending on what
    sra_manager.py's row dicts actually contain (see module-level caveat
    above), and callers must treat ASSAY_UNKNOWN as genuinely uninformative,
    not as evidence of anything.
    """
    combined_text = " ".join(
        str(row.get(field, "")) for field in _POSSIBLE_TEXT_FIELDS if row.get(field)
    ).lower()
    if not combined_text.strip():
        return {
            "assay_type": ASSAY_UNKNOWN, "matched_keyword": None, "matched_field": None,
            "message": (
                "No title/description text field was available in this "
                "accession's metadata to check -- library strategy alone "
                "cannot distinguish bulk from single-cell/spatial data "
                "(SRA's schema doesn't record this distinction directly). "
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
            "known single-cell, spatial, or bulk keyword -- inconclusive. "
            "Consider checking this accession's own SRA/GEO page manually "
            "before proceeding."
        ),
    }

# A sample-index read (I1) is always short and NOT one of the real cell-
# barcode+UMI lengths in chemistry_manager's own catalog -- used as a
# quick first-pass filter before applying the R1-length heuristic to
# whatever remains.
_I1_MAX_LEN = 12


def classify_output_files(fastq_paths, n_reads=2000):
    """
    Given a list of FASTQ file paths from a single SRA run (e.g.
    fasterq-dump's --split-files output: SRRxxxx_1.fastq, _2.fastq,
    possibly _3.fastq), classify each by its dominant read length into
    ROLE_R1 / ROLE_R2 / ROLE_I1 / ROLE_UNKNOWN.

    Returns a dict {path: {"role": ROLE_*, "length": int or None,
    "chemistry_candidates": [...]}} -- chemistry_candidates is populated
    (from chemistry_manager.CHEMISTRY_CATALOG's own R1-length index) only
    for files classified as ROLE_R1, exactly mirroring what
    detect_chemistry() itself would report for that length.

    This is a HEURISTIC, exactly like detect_chemistry() -- the caller
    (singlecell_workspace.py) MUST surface these classifications for
    user confirmation/override, never treat them as certain. With only
    2 files, the shorter is classified R1 and the longer R2. With 3
    files, the shortest is classified I1, and of the remaining two the
    shorter is R1 and longer is R2. Any other file count, or any
    resulting R1 candidate whose length doesn't match ANY known
    chemistry in CHEMISTRY_CATALOG, is left as ROLE_UNKNOWN for every
    file in that run -- signaling the caller to fall back to manual
    per-file assignment and consider calling
    detect_likely_bam_derived_issue() for an explicit warning.
    """
    lengths = {}
    for path in fastq_paths:
        lengths[path] = chem.dominant_read_length(path, n_reads=n_reads)

    result = {p: {"role": ROLE_UNKNOWN, "length": lengths[p], "chemistry_candidates": []} for p in fastq_paths}

    valid_paths = [p for p in fastq_paths if lengths[p] is not None]
    if len(valid_paths) not in (2, 3):
        return result  # can't confidently classify -- see docstring

    by_length = sorted(valid_paths, key=lambda p: lengths[p])

    if len(by_length) == 3:
        i1_path, r1_path, r2_path = by_length[0], by_length[1], by_length[2]
        if lengths[i1_path] > _I1_MAX_LEN:
            return result  # shortest file isn't index-read-short -- don't force a guess
        result[i1_path]["role"] = ROLE_I1
    else:
        r1_path, r2_path = by_length[0], by_length[1]

    r1_len = lengths[r1_path]
    r1_candidates = chem._R1_LEN_INDEX.get(r1_len, [])
    if not r1_candidates:
        return result  # R1 candidate's length matches no known chemistry at all -- don't force a guess

    result[r1_path]["role"] = ROLE_R1
    result[r1_path]["chemistry_candidates"] = r1_candidates
    result[r2_path]["role"] = ROLE_R2
    return result


def detect_likely_bam_derived_issue(classification):
    """
    Given classify_output_files()'s output for one run, return a
    plain-language warning string if the result looks like a likely
    sign of a non-standard or BAM-derived SRA submission (see module
    docstring point 2) -- or None if classification succeeded normally.
    This is a POST-HOC heuristic on the classification OUTCOME, not a
    lookup of the accession's actual original deposited format (which
    is not reliably queryable via esearch/efetch metadata alone).
    """
    roles = {info["role"] for info in classification.values()}
    if roles == {ROLE_UNKNOWN}:
        return (
            "⚠️ Could not confidently identify cell-barcode/UMI vs. cDNA reads "
            "from this run's downloaded file(s). This can happen when an "
            "accession was originally deposited as a 10x Genomics BAM file "
            "rather than plain FASTQ -- those BAMs store barcode/UMI "
            "information in special tags that ordinary SRA download tools "
            "don't reconstruct correctly. Check this accession's own SRA "
            "page for its 'Type' field (Type = TenX indicates this case) -- "
            "if so, 10x's own `bamtofastq` tool is required instead of this "
            "pipeline's standard NCBI/SRA fetch, which is not yet automated "
            "here. You can still assign roles manually below if you believe "
            "these files are usable despite this warning."
        )
    return None


def build_prefetch_command(accession, output_dir):
    return ["prefetch", accession, "--output-directory", output_dir]


def build_fasterq_dump_command(accession, sra_file_dir, output_dir, threads=4):
    "fasterq-dump with --split-files (NOT --split-3, which silently drops the sample-index read some depositions include as a 3rd file) and --include-technical (keeps that index read as its own output file rather than discarding it, so classify_output_files() actually gets a chance to see and classify it)."
    return [
        "fasterq-dump", "--split-files", "--include-technical",
        "--threads", str(threads), "--outdir", output_dir,
        os.path.join(sra_file_dir, accession),
    ]


def download_and_classify_run(accession, project_fastq_dir, sample_name, threads=4, subprocess_runner=None):
    """
    Download one SRA run accession (prefetch -> fasterq-dump) and
    classify its output files into R1/R2/(I1) roles. Returns a dict:
        {
            "success": bool,
            "message": str,
            "classification": {path: {"role", "length", "chemistry_candidates"}} or None,
            "bam_warning": str or None,   # see detect_likely_bam_derived_issue()
        }

    subprocess_runner is injectable (defaults to subprocess.run) purely
    for testing without invoking the real prefetch/fasterq-dump
    binaries -- see this module's test suite.

    Does NOT rename files into the standard *_R1_001.fastq.gz naming
    convention itself -- singlecell_workspace.py does that only AFTER
    the user has confirmed (or corrected) the role classification this
    function proposes, exactly mirroring how chemistry auto-detection
    is proposed-then-confirmed rather than applied silently.
    """
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

    produced_files = sorted([
        os.path.join(work_dir, f) for f in os.listdir(work_dir)
        if f.startswith(accession) and (f.endswith(".fastq") or f.endswith(".fastq.gz"))
    ])
    if not produced_files:
        return {"success": False, "message": f"fasterq-dump produced no output files for {accession}.", "classification": None, "bam_warning": None}

    classification = classify_output_files(produced_files)
    bam_warning = detect_likely_bam_derived_issue(classification)
    return {
        "success": True,
        "message": f"Downloaded and classified {len(produced_files)} file(s) for {accession}.",
        "classification": classification,
        "bam_warning": bam_warning,
    }


def finalize_role_assignment(classification_overrides, project_fastq_dir, sample_name):
    """
    Given a FINAL, user-confirmed {path: role} mapping (after the UI has
    shown download_and_classify_run()'s proposed classification and let
    the user correct any of it), rename/move each R1 and R2 file into
    this project's standard fastq_dir() using the normal
    <sample>_S1_R1_001.fastq.gz / <sample>_S1_R2_001.fastq.gz naming
    convention -- so singlecell_ingestion_manager.find_r1_r2_pairs()
    picks them up exactly like any locally-uploaded or browsed FASTQ
    pair, with no special-casing needed anywhere downstream of this
    point. I1 (sample index) and UNKNOWN files are left in place
    (not moved into fastq_dir()) since they aren't used by anything
    downstream.

    Returns the dict of final {role: destination_path} for R1/R2 only.
    """
    destinations = {}
    for path, role in classification_overrides.items():
        if role not in (ROLE_R1, ROLE_R2):
            continue
        ext = ".fastq.gz" if path.endswith(".gz") else ".fastq"
        dest = os.path.join(project_fastq_dir, f"{sample_name}_S1_{role}_001{ext}")
        shutil.move(path, dest)
        destinations[role] = dest
    return destinations
