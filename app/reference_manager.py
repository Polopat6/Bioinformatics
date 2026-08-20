"""
reference_manager.py

Handles fetching or accepting reference genome/transcriptome files needed
by the RNA Alignment & Counts workspace (alignment_workspace.py), for
both pre-loaded model organisms and custom (non-model) species.

Kept as its own module (rather than folded into alignment_workspace.py)
since reference handling involves network downloads, file decompression,
and format-specific logic (Ensembl vs. NCBI URL structures) that is
logically distinct from the Streamlit UI/workflow code.

--- Reference sources used ---

Human, Mouse, Yeast, Drosophila, C. elegans, Zebrafish: Ensembl's main
FTP site, using its "current_fasta" and "current_gtf" symlink
directories. These always point to the latest Ensembl release for that
species, so URLs here don't go stale as new Ensembl releases come out.
The exact filename within that directory still includes a release number
that changes over time (e.g. "Homo_sapiens.GRCh38.113.gtf.gz"), so we
fetch the directory listing and pattern-match the correct filename
rather than hardcoding it.

E. coli: Ensembl's bacterial genome FTP structure is organized
differently (collection-based, less predictable), so instead we use
NCBI RefSeq's assembly-specific FTP path for E. coli K-12 MG1655
(GCF_000005845.2 / ASM584v2). NCBI keeps every assembly version at a
permanently fixed path, so this URL will not go stale.

--- Mitochondrial genome verification (2026-08-17) ---
A real reported concern: a preset (auto-downloaded) reference's Step 5
UI only ever checks "does this directory already have files in it?" to
decide whether to show "✅ already prepared" -- it never actually
verifies WHAT is in those files. Combined with a separate, now-fixed
single-cell Cell-level QC bug where mitochondrial-gene detection relied
solely on gene-symbol matching, this made it genuinely impossible for a
user to tell whether a "0% mitochondrial" QC result meant "your cells
are very healthy" or "this specific downloaded reference copy is
missing the mitochondrial chromosome/contig entirely" (e.g. from an
older cached copy, a partial/interrupted historical download, or an
upstream Ensembl/NCBI change).

genome_fasta_contains_mito_contig() below provides a fast, second,
INDEPENDENT verification signal -- checking the actual downloaded
genome FASTA's own sequence headers (not just the GTF single-cell's
own get_mito_gene_ids_from_gtf() already checks) for a recognized
mitochondrial contig name. This only scans header lines (each genome
FASTA's ">chromosome_name ..." lines), never the sequence data itself,
so it stays fast even for a multi-gigabyte primary-assembly file.
Wired into singlecell_workspace.py's Step 5 preset-reference section
(see that module for details) to surface a concrete, actionable
"re-download this reference" option specifically when mitochondrial
content can't be verified, rather than only offering a generic
"force re-download" as a buried advanced/undifferentiated action.
"""

import fcntl
import gzip
import os
import re
import shutil
import ssl
import subprocess
import time
import urllib.request
import urllib.error
from datetime import datetime

# On some systems (notably macOS with Python installed from python.org),
# Python's SSL module doesn't automatically use the operating system's
# trusted certificate store, causing every HTTPS download to fail with
# "CERTIFICATE_VERIFY_FAILED: self signed certificate in certificate
# chain" even though the remote server (Ensembl/NCBI) has a perfectly
# valid certificate. We work around this by explicitly building an SSL
# context using the `certifi` package's bundled CA certificates, which
# is independent of whatever the local Python installation's OS-level
# cert configuration looks like.
try:
    import certifi
    _SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
except ImportError:
    # certifi isn't installed -- fall back to the default context. This
    # will hit the same CERTIFICATE_VERIFY_FAILED error on affected
    # systems, but we don't want a missing optional dependency to be a
    # hard crash at import time.
    _SSL_CONTEXT = ssl.create_default_context()

# ---------------------------------------------------------------------------
# Model organism catalog
# ---------------------------------------------------------------------------

REFERENCE_CATALOG = {
    "human": {
        "label": "Human (Homo sapiens, GRCh38)",
        "source": "ensembl",
        "species_dir": "homo_sapiens",
        "species_name": "Homo_sapiens",
        "assembly": "GRCh38",
        "genome_fasta_pattern": r"Homo_sapiens\.GRCh38\.dna\.primary_assembly\.fa\.gz",
        "cdna_fasta_pattern": r"Homo_sapiens\.GRCh38\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Homo_sapiens\.GRCh38\.\d+\.gtf\.gz",
    },
    "mouse": {
        "label": "Mouse (Mus musculus, GRCm39)",
        "source": "ensembl",
        "species_dir": "mus_musculus",
        "species_name": "Mus_musculus",
        "assembly": "GRCm39",
        "genome_fasta_pattern": r"Mus_musculus\.GRCm39\.dna\.primary_assembly\.fa\.gz",
        "cdna_fasta_pattern": r"Mus_musculus\.GRCm39\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Mus_musculus\.GRCm39\.\d+\.gtf\.gz",
    },
    "yeast": {
        "label": "Yeast (Saccharomyces cerevisiae, R64-1-1)",
        "source": "ensembl",
        "species_dir": "saccharomyces_cerevisiae",
        "species_name": "Saccharomyces_cerevisiae",
        "assembly": "R64-1-1",
        # Yeast's small genome is distributed as "toplevel" rather than
        # "primary_assembly" on Ensembl.
        "genome_fasta_pattern": r"Saccharomyces_cerevisiae\.R64-1-1\.dna\.toplevel\.fa\.gz",
        "cdna_fasta_pattern": r"Saccharomyces_cerevisiae\.R64-1-1\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Saccharomyces_cerevisiae\.R64-1-1\.\d+\.gtf\.gz",
    },
    "drosophila": {
        "label": "Fruit Fly (Drosophila melanogaster, BDGP6.54)",
        "source": "ensembl",
        "species_dir": "drosophila_melanogaster",
        "species_name": "Drosophila_melanogaster",
        "assembly": "BDGP6.54",
        "genome_fasta_pattern": r"Drosophila_melanogaster\.BDGP6\.54\.dna\.toplevel\.fa\.gz",
        "cdna_fasta_pattern": r"Drosophila_melanogaster\.BDGP6\.54\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Drosophila_melanogaster\.BDGP6\.54\.\d+\.gtf\.gz",
    },
    "celegans": {
        "label": "Roundworm (Caenorhabditis elegans, WBcel235)",
        "source": "ensembl",
        "species_dir": "caenorhabditis_elegans",
        "species_name": "Caenorhabditis_elegans",
        "assembly": "WBcel235",
        "genome_fasta_pattern": r"Caenorhabditis_elegans\.WBcel235\.dna\.toplevel\.fa\.gz",
        "cdna_fasta_pattern": r"Caenorhabditis_elegans\.WBcel235\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Caenorhabditis_elegans\.WBcel235\.\d+\.gtf\.gz",
    },
    "zebrafish": {
        "label": "Zebrafish (Danio rerio, GRCz11)",
        "source": "ensembl",
        "species_dir": "danio_rerio",
        "species_name": "Danio_rerio",
        "assembly": "GRCz11",
        "genome_fasta_pattern": r"Danio_rerio\.GRCz11\.dna\.primary_assembly\.fa\.gz",
        "cdna_fasta_pattern": r"Danio_rerio\.GRCz11\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Danio_rerio\.GRCz11\.\d+\.gtf\.gz",
    },
    "ecoli": {
        "label": "E. coli (K-12 MG1655)",
        "source": "ncbi",
        # NCBI RefSeq assembly GCF_000005845.2 (ASM584v2) -- a fixed,
        # permanent path that does not change over time.
        "ncbi_base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
        "genome_fasta_filename": "GCF_000005845.2_ASM584v2_genomic.fna.gz",
        "gtf_filename": "GCF_000005845.2_ASM584v2_genomic.gtf.gz",
        "cdna_fasta_filename": None,
        "no_introns": True,
    },
}

ENSEMBL_BASE_URL = "https://ftp.ensembl.org/pub"


# ---------------------------------------------------------------------------
# Mitochondrial genome verification (2026-08-17) -- see module docstring
# ---------------------------------------------------------------------------

# Kept identical to sc_cellqc_manager.py's own _MITO_SEQNAME_ALIASES set
# (that module's the authoritative GTF-gene-level check; this one is a
# lightweight, FASTA-header-only cross-check) -- both independently
# maintained since they check different file types, but should be kept
# in sync if either is ever extended with new naming conventions.
_MITO_SEQNAME_ALIASES = {"mt", "chrm", "chrmt", "m", "mitochondrion", "mtdna"}


def genome_fasta_contains_mito_contig(fasta_path, max_headers_to_scan=5000):
    """
    Check whether a genome FASTA file's own sequence headers include a
    recognized mitochondrial contig name (e.g. ">MT dna:chromosome...",
    ">chrM ...") -- WITHOUT reading or scanning any actual sequence data,
    so this stays fast even for a multi-gigabyte primary-assembly FASTA.

    This exists specifically to let the UI answer "does this ALREADY-
    DOWNLOADED reference actually include the mitochondrial genome?"
    directly and quickly, rather than the previous behavior of only
    checking "does this directory have files in it at all?" -- which
    says nothing about their actual content, and left a real, reported
    gap where a user had no way to distinguish a genuinely mito-free
    reference (which would silently produce misleading "0% mitochondrial"
    single-cell QC results) from a healthy one, without re-downloading
    or manually inspecting the file.

    fasta_path: path to a genome FASTA file (may be a large primary-
        assembly/toplevel file -- this function is header-scan-only and
        does not load sequence content into memory).
    max_headers_to_scan: safety cap on how many ">"-prefixed header
        lines to examine before giving up -- a real genome FASTA's
        mitochondrial contig (if present) is essentially always
        encountered within the first few dozen headers (contigs are
        conventionally ordered chromosome 1, 2, 3, ..., MT/X/Y at the
        end, or similar), so this cap is a defensive safety valve, not
        a normal-path constraint.

    Returns True if a mitochondrial-named header was found, False
    otherwise (including if the file doesn't exist or can't be read).
    """
    if not fasta_path or not os.path.isfile(fasta_path):
        return False

    headers_scanned = 0
    try:
        with open(fasta_path, "r", errors="replace") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                headers_scanned += 1
                if headers_scanned > max_headers_to_scan:
                    break
                # A FASTA header's sequence ID is the first whitespace-
                # delimited token after ">" (e.g. ">MT dna:chromosome
                # MT:GRCh38:MT:1:16569:1 REF" -> "MT").
                seq_id = line[1:].split()[0].strip().lower()
                if seq_id in _MITO_SEQNAME_ALIASES:
                    return True
    except OSError:
        return False

    return False


def verify_preset_reference_mito_content(genome_fasta_path, gtf_path):
    """
    Combined mitochondrial-content verification for a PRESET (auto-
    downloaded) reference's already-on-disk genome FASTA + GTF -- used
    by singlecell_workspace.py's Step 5 to answer "can this already-
    downloaded reference's mitochondrial content actually be confirmed?"
    immediately, without needing to re-download or run a full Cell-level
    QC pass first.

    Uses TWO independent checks (mirroring the same two-method approach
    already used for single-cell QC's own mito-gene detection):
      1. genome_fasta_contains_mito_contig() above -- does the genome
         FASTA's own sequence headers include a mitochondrial contig?
      2. A local import of sc_cellqc_manager.get_mito_gene_ids_from_gtf()
         -- does the GTF contain any gene-level features on a
         mitochondrial contig? (Imported locally, not at module level,
         to avoid a circular/unnecessary import for every other
         reference_manager.py caller that has nothing to do with
         single-cell QC.)

    Returns a dict:
        {
            "fasta_has_mito_contig": bool,
            "gtf_mito_gene_count": int,
            "verified": bool,  # True only if BOTH checks succeed
        }
    """
    fasta_ok = genome_fasta_contains_mito_contig(genome_fasta_path)

    gtf_gene_count = 0
    if gtf_path and os.path.isfile(gtf_path):
        try:
            # Local import: this module (reference_manager.py) is used
            # by the Bulk RNA-Seq pipeline too, which has no reason to
            # depend on single_cell/sc_cellqc_manager.py -- importing
            # only here, only when this specific verification function
            # is actually called, avoids adding an unconditional
            # cross-pipeline import at module load time.
            import sc_cellqc_manager as _cellqc
            gtf_gene_count = len(_cellqc.get_mito_gene_ids_from_gtf(gtf_path))
        except ImportError:
            gtf_gene_count = 0

    return {
        "fasta_has_mito_contig": fasta_ok,
        "gtf_mito_gene_count": gtf_gene_count,
        "verified": fasta_ok and gtf_gene_count > 0,
    }


# ---------------------------------------------------------------------------
# Concurrency-safe shared resource preparation
# ---------------------------------------------------------------------------
#
# Used for PRESET species' shared, project-independent reference files
# and Salmon/STAR indices (see project_manager.py's shared_reference_dir()
# and friends) -- NOT needed for custom uploads, since those are always
# scoped to a single project and therefore can never race with another
# project's concurrent access to the same files.
#
# The core problem this solves: if two different projects (each
# potentially running in its own Streamlit session/process) select the
# same preset species at close to the same time, and neither has
# downloaded/built it yet, naively checking "does this shared path
# exist?" and downloading/building if not would race -- both processes
# could see "not there yet" and both start writing into the same shared
# directory at the same time, corrupting each other's output, or one
# process could start using a reference that a second process is
# simultaneously overwriting mid-download/mid-build.
#
# The fix has two parts, used together:
#
#   1. File locking via fcntl.flock -- an OS-level, cross-PROCESS lock
#      (unlike a Python threading.Lock, which only protects against
#      concurrent threads within a single process, not concurrent
#      Streamlit sessions/processes, which is how real concurrent
#      users of this app are actually structured). A lock file exists
#      per shared resource; only one process can hold it at a time,
#      and every other process trying to acquire it simply waits until
#      the current holder releases it (or crashes -- flock's lock is
#      automatically released by the OS if the holding process dies
#      for any reason, so there's no separate "detect and clean up a
#      stale lock" logic needed, unlike a naive lock-file-existence
#      scheme).
#
#   2. Build into a private temporary directory, then perform a single
#      atomic os.rename() into the final shared location ONLY once the
#      build fully succeeds. This guarantees the final shared directory
#      either doesn't exist yet, or contains a fully-completed prior
#      build -- there is no window where a concurrent reader could see
#      a half-downloaded/half-indexed shared directory, and a crashed
#      or failed build never leaves partial output behind at the final
#      path for a later attempt to mistake for a successful one.
#
# Together, these mean: the FIRST project to request a given preset
# species does the real work (download or index build) while holding
# the lock; every OTHER project requesting the same species at the
# same time waits (with live progress feedback via
# wait_message_callback) until the first one finishes, then simply
# reuses the now-ready shared files -- no duplicate downloads, no
# duplicate index builds, no risk of two processes corrupting the same
# directory.

def _resource_is_ready(resource_dir):
    """
    A shared resource is considered "ready" if its directory exists and
    contains at least one file/subdirectory -- this is checked both
    before attempting to acquire a lock (fast path: nothing to do) and
    immediately after acquiring one (in case another process finished
    building it while this one was waiting for the lock).
    """
    return os.path.isdir(resource_dir) and len(os.listdir(resource_dir)) > 0

"""
PATCH for reference_manager.py -- shared concurrency-safety fix.

WHERE TO INSERT: directly after the existing `_resource_is_ready()`
function and BEFORE `ensure_shared_resource()` (both currently live in
the "Concurrency-safe shared resource preparation" section).

WHY: ensure_shared_resource() already protects against two BUILDERS
racing (via the .lock file + fcntl.flock), but it does nothing to stop
a READER from seeing a stale result. A force=True rebuild builds into a
private temp directory and only replaces the OLD resource via a single
os.rename() at the very end -- so the old (stale) resource stays fully
present on disk, and reads as "ready" to every existing check
(_resource_is_ready, qm.star_index_exists, qm.salmon_index_exists), for
the ENTIRE rebuild duration. This is exactly what let a user navigate
away mid-rebuild and see a false "already built" status while the real
build was still running (confirmed via `top` on the HPC host).

The lock file that ensure_shared_resource() already maintains is the
perfect signal to close this gap -- we just need to let READERS probe
it non-blockingly, without needing to win the lock themselves.

No other function in this file needs to change. Every CALLER (in both
alignment_workspace.py and singlecell_workspace.py) needs to switch
from checking readiness via _resource_is_ready()/qm.*_index_exists()
ALONE to using is_shared_resource_ready() / resource_build_in_progress()
below instead -- see the accompanying updated alignment_workspace.py
and singlecell_workspace.py for exactly where.
"""


def resource_build_in_progress(resource_dir):
    """
    Non-blocking probe: is ANY process (this session, another browser
    tab, another pipeline, another user) currently holding the build
    lock for this shared resource RIGHT NOW?

    Every "is this shared resource ready to use?" check across BOTH the
    Bulk RNA-Seq and Single-cell pipelines must go through this (via
    is_shared_resource_ready() below) instead of a bare
    os.path.isdir()/os.listdir() check or a single marker-file check
    (qm.star_index_exists's "SAindex" check, qm.salmon_index_exists's
    "info.json" check, etc.) -- none of those can distinguish a
    genuinely finished build from one that's actively being replaced in
    the background, since the OLD resource is left fully intact right
    up until ensure_shared_resource()'s final atomic os.rename().

    Returns True if a build is in progress right now (lock currently
    held by someone), False if the lock is free or has never been
    taken at all (e.g. this resource has never been built).
    """
    lock_path = resource_dir.rstrip(os.sep) + ".lock"
    if not os.path.isfile(lock_path):
        return False  # nobody has ever attempted to build this -- nothing in progress
    lock_file = open(lock_path, "a+")
    try:
        # Non-blocking probe only -- if we get the lock instantly, no
        # one else holds it (free). We immediately release it again;
        # this function only ever reports status, it never builds
        # anything or holds the lock itself.
        fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        fcntl.flock(lock_file, fcntl.LOCK_UN)
        return False
    except BlockingIOError:
        return True  # someone else currently holds the lock -> a build is running
    finally:
        lock_file.close()


def is_shared_resource_ready(resource_dir):
    """
    THE authoritative "is this shared resource actually usable right
    now?" check -- combines file-presence with the in-progress-build
    probe above. A populated directory during an in-progress
    force-rebuild is explicitly NOT ready, even though
    _resource_is_ready() alone (or any single marker-file check) would
    incorrectly say it is.

    Use this (or resource_build_in_progress() directly, when you need
    to distinguish "not ready because still building" from "not ready
    because it's never been built at all" for messaging purposes) at
    every call site that currently only checks
    _resource_is_ready()/qm.star_index_exists()/qm.salmon_index_exists().
    """
    return _resource_is_ready(resource_dir) and not resource_build_in_progress(resource_dir)


def ensure_shared_resource(resource_dir, build_fn, wait_message_callback=None,
                            poll_interval=5, force=False):
    """
    Ensure a shared, project-independent resource (a preset species'
    downloaded reference files, or its built Salmon/STAR index) exists
    at resource_dir, safely handling the case where multiple projects
    request the same resource at close to the same time -- see the
    module-level comment above for the full concurrency rationale.

    resource_dir: the FINAL destination directory for the resource
        (e.g. project_manager.shared_cdna_fasta_dir("human")). This
        function guarantees resource_dir only ever comes into
        existence via a single atomic rename from a completed temp
        build -- callers can safely assume that if resource_dir exists
        and is non-empty, it represents a fully-completed prior build,
        never a partial/in-progress one.

    build_fn(temp_dir) -> (success: bool, message: str): callback that
        performs the actual download/extraction/index-build, writing
        its COMPLETE output into temp_dir (never directly into
        resource_dir) -- this function handles moving the finished
        result into resource_dir atomically only after build_fn
        reports success. build_fn should behave exactly like it would
        if temp_dir were the real final destination (e.g. an existing
        reference_manager.py function like get_transcriptome_fasta_for_salmon
        or download_genome_and_gtf, called with temp_dir as its
        dest_dir argument, needs no changes to support this).

    wait_message_callback(elapsed_seconds): optional callback invoked
        periodically (about every poll_interval seconds) while this
        call is blocked waiting for a DIFFERENT process's build of the
        SAME resource to finish -- lets the UI show live "another
        project is currently preparing this reference, please wait..."
        feedback instead of an unexplained silent pause. Never called
        if this process acquires the lock immediately (the common,
        no-contention case) or if the resource was already ready
        before any locking was attempted.

    force: if True, skips the "already ready, reuse it" fast path
        entirely and always rebuilds -- used to support an explicit
        user-requested "re-download"/"re-index" action on what is now
        a SHARED resource. This is still fully lock-protected (so it
        can never race with another process concurrently trying to
        build or reuse the exact same resource), but note the one risk
        this does NOT solve: if some OTHER project's alignment/
        quantification run is actively reading the OLD version of
        this shared resource while a force rebuild replaces it, that
        other run could fail or behave unexpectedly once the files
        underneath it change -- this function has no way to know
        whether a shared resource is currently "in use" elsewhere, so
        a force rebuild should be used thoughtfully (e.g. not while
        you know other projects/users are actively running alignment
        against the same species).

    Returns (success: bool, message: str, built_by_this_call: bool).
    built_by_this_call is True only if THIS call actually performed the
    download/build (as opposed to finding the resource already ready,
    or finding it ready after waiting for another process to finish
    it) -- useful for the caller to decide whether to show a "built
    successfully" vs. a "was already available" message.
    """
    # Fast path: already fully built (the common case after the first
    # project sets up a given species) -- no lock needed at all. Never
    # taken when force=True, since the whole point of force is to
    # rebuild even though a ready copy already exists.
    # Fast path: already fully built (the common case after the first
    # project sets up a given species) -- no lock needed at all. Never
    # taken when force=True, since the whole point of force is to
    # rebuild even though a ready copy already exists. Uses
    # is_shared_resource_ready() (not the bare _resource_is_ready())
    # so this fast path correctly does NOT trigger while another
    # process is actively mid-rebuild -- otherwise a caller could race
    # past a real rebuild and silently reuse the stale, about-to-be-
    # replaced old copy instead of waiting for the lock below.
    if not force and is_shared_resource_ready(resource_dir):
        return True, "Already available (shared reference, reused from a previous build).", False
    lock_path = resource_dir.rstrip(os.sep) + ".lock"
    os.makedirs(os.path.dirname(lock_path) or ".", exist_ok=True)

    lock_file = open(lock_path, "a+")
    start_time = time.time()
    try:
        # Try to acquire the lock without blocking first, so a process
        # that gets it immediately (no contention) never even enters
        # the wait-and-poll loop or invokes wait_message_callback at
        # all.
        acquired = False
        try:
            fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
            acquired = True
        except BlockingIOError:
            acquired = False

        while not acquired:
            elapsed = time.time() - start_time
            if wait_message_callback:
                wait_message_callback(elapsed)
            time.sleep(poll_interval)
            try:
                fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                acquired = True
            except BlockingIOError:
                acquired = False

        # Lock acquired. Re-check: another process may have finished
        # building this exact resource while we were waiting for the
        # lock -- the single most common reason to be waiting at all.
        # Skipped entirely when force=True, for the same reason as the
        # pre-lock fast path above.
        if not force and _resource_is_ready(resource_dir):
            return True, "Already available (another project finished preparing it while this one was waiting).", False

        # We hold the lock and either the resource isn't ready yet, or
        # force=True was requested -- build it ourselves, into a
        # private temp directory scoped to this process's PID (so two
        # sequential failed attempts, or a process crash mid-build,
        # never collide with a later attempt's temp directory).
        temp_dir = resource_dir.rstrip(os.sep) + f".building.{os.getpid()}"
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)
        os.makedirs(temp_dir, exist_ok=True)

        try:
            success, message = build_fn(temp_dir)
        except Exception as e:
            success, message = False, f"Unexpected error while building shared resource: {e}"

        if success:
            # Atomic on the same filesystem -- resource_dir either
            # doesn't exist yet, or (the normal case for force=True,
            # or the rare case of a previous build somehow leaving
            # something there without a lock) gets fully replaced in
            # one step. There is no window where a concurrent reader
            # could observe a half-populated resource_dir.
            if os.path.isdir(resource_dir):
                shutil.rmtree(resource_dir, ignore_errors=True)
            os.rename(temp_dir, resource_dir)
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)

        return success, message, True
    finally:
        fcntl.flock(lock_file, fcntl.LOCK_UN)
        lock_file.close()


# ---------------------------------------------------------------------------
# Low-level download helpers
# ---------------------------------------------------------------------------

def _fetch_url_text(url, timeout=30):
    """Fetch a URL's content as text (used for directory listing pages)."""
    with urllib.request.urlopen(url, timeout=timeout, context=_SSL_CONTEXT) as response:
        return response.read().decode("utf-8", errors="ignore")


def _resolve_ensembl_filename(directory_url, filename_pattern):
    """
    Fetch an Ensembl FTP directory listing page and find the filename
    matching filename_pattern (a regex). This avoids hardcoding release
    numbers that change with every new Ensembl release.

    Returns the matched filename, or None if no match was found.
    """
    try:
        html = _fetch_url_text(directory_url)
    except (urllib.error.URLError, TimeoutError) as e:
        raise RuntimeError(f"Could not reach Ensembl FTP at {directory_url}: {e}")

    matches = re.findall(filename_pattern, html)
    if not matches:
        return None
    # Directory listings may reference the same filename multiple times
    # (e.g. in both an href and display text) -- dedupe while preserving
    # order, then take the first (they should all be identical anyway).
    return matches[0]


def _download_file(url, dest_path, progress_callback=None):
    """
    Download a file from url to dest_path, optionally reporting progress
    via progress_callback(bytes_downloaded, total_bytes).

    Returns (success: bool, error_message: str or None).
    """
    try:
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        with urllib.request.urlopen(url, timeout=60, context=_SSL_CONTEXT) as response:
            total_size = int(response.headers.get("Content-Length", 0))
            downloaded = 0
            chunk_size = 1024 * 1024  # 1MB chunks

            with open(dest_path, "wb") as out_file:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    out_file.write(chunk)
                    downloaded += len(chunk)
                    if progress_callback:
                        progress_callback(downloaded, total_size)
        return True, None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return False, str(e)


def _decompress_gz(gz_path, dest_path):
    """Decompress a .gz file to dest_path. Returns (success, error_message)."""
    try:
        with gzip.open(gz_path, "rb") as f_in:
            with open(dest_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True, None
    except (OSError, gzip.BadGzipFile) as e:
        return False, str(e)


# ---------------------------------------------------------------------------
# Public download functions
# ---------------------------------------------------------------------------

def download_cdna_fasta(species_key, dest_dir, progress_callback=None):
    """
    Download (or otherwise prepare) a transcriptome/cDNA FASTA for the
    given species key, for use with Salmon.

    Returns (success: bool, fasta_path_or_None, message: str).
    """
    entry = REFERENCE_CATALOG[species_key]
    dest_path = os.path.join(dest_dir, f"{species_key}.cdna.fa")
    gz_path = dest_path + ".gz"

    if entry["source"] == "ensembl":
        directory_url = f"{ENSEMBL_BASE_URL}/current_fasta/{entry['species_dir']}/cdna/"
        filename = _resolve_ensembl_filename(directory_url, entry["cdna_fasta_pattern"])
        if not filename:
            return False, None, (
                f"Could not find a matching cDNA FASTA file at {directory_url}. "
                "Ensembl may have changed their file naming -- please check "
                "manually or use the custom upload option instead."
            )
        url = directory_url + filename

    elif entry["source"] == "ncbi":
        if entry.get("cdna_fasta_filename"):
            url = f"{entry['ncbi_base_url']}/{entry['cdna_fasta_filename']}"
        else:
            return False, None, (
                "NCBI does not provide a standalone transcriptome FASTA for "
                "this organism. Use the 'genome + GTF' download instead, "
                "which will automatically extract transcript sequences."
            )
    else:
        return False, None, f"Unknown reference source '{entry['source']}'."

    success, error = _download_file(url, gz_path, progress_callback)
    if not success:
        return False, None, f"Download failed: {error}"

    success, error = _decompress_gz(gz_path, dest_path)
    if not success:
        return False, None, f"Decompression failed: {error}"

    os.remove(gz_path)
    return True, dest_path, f"Downloaded and extracted cDNA FASTA from {url}"


def get_transcriptome_fasta_for_salmon(species_key, dest_dir, progress_callback=None):
    """
    High-level helper that guarantees a usable transcriptome FASTA ends
    up on disk for Salmon, regardless of whether the species' catalog
    entry has a ready-made cDNA FASTA or not.

    Returns (success: bool, fasta_path_or_None, message: str).
    """
    entry = REFERENCE_CATALOG[species_key]

    if entry.get("cdna_fasta_filename") or entry.get("cdna_fasta_pattern"):
        return download_cdna_fasta(species_key, dest_dir, progress_callback)

    success, paths, message = download_genome_and_gtf(species_key, dest_dir, progress_callback)
    if not success:
        return False, None, f"Reference download failed: {message}"

    genome_fasta_path, gtf_path = paths
    extracted_fasta_path = os.path.join(dest_dir, f"{species_key}.cdna.fa")

    if entry.get("no_introns"):
        success, extraction_message = extract_gene_level_transcripts(
            genome_fasta_path, gtf_path, extracted_fasta_path
        )
    else:
        success, extraction_message = extract_transcripts_with_gffread(
            genome_fasta_path, gtf_path, extracted_fasta_path
        )

    if not success:
        return False, None, (
            f"Downloaded genome + annotation successfully, but transcript "
            f"extraction failed: {extraction_message}"
        )

    extraction_tool_note = (
        "direct gene-coordinate extraction (appropriate for organisms "
        "with no introns)" if entry.get("no_introns") else "gffread"
    )
    return True, extracted_fasta_path, (
        f"{message}. Transcript sequences were then extracted from the "
        f"genome + annotation using {extraction_tool_note}, since this "
        "species does not publish a standalone transcriptome FASTA."
    )


# ---------------------------------------------------------------------------
# Ensembl GTF -> GFF3 fallback (2026-08-17)
# ---------------------------------------------------------------------------

def _download_ensembl_annotation(entry, dest_dir, progress_callback=None):
    """
    Resolve and download this species' gene annotation from Ensembl,
    preferring GTF but automatically falling back to GFF3 (+ conversion
    to GTF via gffread) if the "current_gtf" directory 404s.

    Returns (success: bool, gtf_path_or_None, message: str).
    """
    gtf_dir_url = f"{ENSEMBL_BASE_URL}/current_gtf/{entry['species_dir']}/"
    gtf_filename = None
    try:
        gtf_filename = _resolve_ensembl_filename(gtf_dir_url, entry["gtf_pattern"])
    except RuntimeError:
        gtf_filename = None

    if gtf_filename:
        gtf_url = gtf_dir_url + gtf_filename
        gtf_gz = os.path.join(dest_dir, "annotation.gtf.gz")
        gtf_path = os.path.join(dest_dir, "annotation.gtf")
        success, error = _download_file(gtf_url, gtf_gz, progress_callback)
        if not success:
            return False, None, f"GTF download failed: {error}"
        success, error = _decompress_gz(gtf_gz, gtf_path)
        if not success:
            return False, None, f"GTF decompression failed: {error}"
        os.remove(gtf_gz)
        return True, gtf_path, f"Downloaded annotation (GTF) from {gtf_url}"

    # --- GFF3 fallback ---
    gff3_dir_url = f"{ENSEMBL_BASE_URL}/current_gff3/{entry['species_dir']}/"
    gff3_pattern = entry["gtf_pattern"].replace("gtf", "gff3")
    try:
        gff3_filename = _resolve_ensembl_filename(gff3_dir_url, gff3_pattern)
    except RuntimeError as e:
        return False, None, (
            f"Could not find a GTF annotation at {gtf_dir_url} (Ensembl "
            f"is in the middle of restructuring their FTP site as of "
            f"2026), and the GFF3 fallback also failed: {e}"
        )
    if not gff3_filename:
        return False, None, (
            f"Could not find a matching GTF OR GFF3 annotation file for "
            f"this species (checked {gtf_dir_url} and {gff3_dir_url}). "
            "Ensembl may have changed their file naming further -- please "
            "check manually or use the custom upload option instead."
        )

    if shutil.which("gffread") is None:
        return False, None, (
            f"Found a GFF3 annotation at {gff3_dir_url}{gff3_filename}, but "
            "no GTF was available at Ensembl's usual location, and gffread "
            "(needed to convert GFF3 -> GTF) is not installed on this "
            "system. gffread is part of this project's environment.yml."
        )

    gff3_url = gff3_dir_url + gff3_filename
    gff3_gz = os.path.join(dest_dir, "annotation.gff3.gz")
    gff3_path = os.path.join(dest_dir, "annotation.gff3")
    gtf_path = os.path.join(dest_dir, "annotation.gtf")

    success, error = _download_file(gff3_url, gff3_gz, progress_callback)
    if not success:
        return False, None, f"GFF3 download failed: {error}"
    success, error = _decompress_gz(gff3_gz, gff3_path)
    if not success:
        return False, None, f"GFF3 decompression failed: {error}"
    os.remove(gff3_gz)

    try:
        subprocess.run(
            ["gffread", gff3_path, "-T", "-o", gtf_path],
            capture_output=True, text=True, check=True, timeout=1800,
        )
    except subprocess.CalledProcessError as e:
        return False, None, f"gffread GFF3->GTF conversion failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, None, "gffread GFF3->GTF conversion timed out after 30 minutes."

    return True, gtf_path, (
        f"Ensembl's 'current_gtf' directory was unavailable for this "
        f"species (Ensembl is in the middle of restructuring their FTP "
        f"site as of 2026) -- downloaded GFF3 from {gff3_url} instead and "
        f"converted it to GTF using gffread."
    )


def download_genome_and_gtf(species_key, dest_dir, progress_callback=None):
    """
    Download the genome FASTA + GTF annotation for the given species key,
    for use with STAR (or for extracting transcripts via gffread if used
    as a Salmon fallback).

    Returns (success: bool, (genome_fasta_path, gtf_path) or None, message: str).
    """
    entry = REFERENCE_CATALOG[species_key]

    genome_gz = os.path.join(dest_dir, f"{species_key}.genome.fa.gz")
    genome_fa = os.path.join(dest_dir, f"{species_key}.genome.fa")
    gtf_path_final = os.path.join(dest_dir, f"{species_key}.annotation.gtf")

    if entry["source"] == "ensembl":
        genome_dir_url = f"{ENSEMBL_BASE_URL}/current_fasta/{entry['species_dir']}/dna/"
        genome_filename = _resolve_ensembl_filename(genome_dir_url, entry["genome_fasta_pattern"])
        if not genome_filename:
            return False, None, (
                f"Could not find a matching genome FASTA file at {genome_dir_url}. "
                "Ensembl may have changed their file naming -- please check "
                "manually or use the custom upload option instead."
            )
        genome_url = genome_dir_url + genome_filename

        success, error = _download_file(genome_url, genome_gz, progress_callback)
        if not success:
            return False, None, f"Genome download failed: {error}"
        success, error = _decompress_gz(genome_gz, genome_fa)
        if not success:
            return False, None, f"Genome decompression failed: {error}"
        os.remove(genome_gz)

        annotation_success, annotation_gtf_path, annotation_message = _download_ensembl_annotation(entry, dest_dir, progress_callback)
        if not annotation_success:
            return False, None, annotation_message
        if annotation_gtf_path != gtf_path_final:
            shutil.move(annotation_gtf_path, gtf_path_final)

        return True, (genome_fa, gtf_path_final), f"Downloaded genome from {genome_url}. {annotation_message}"

    elif entry["source"] == "ncbi":
        genome_url = f"{entry['ncbi_base_url']}/{entry['genome_fasta_filename']}"
        gtf_url = f"{entry['ncbi_base_url']}/{entry['gtf_filename']}"

        gtf_gz = os.path.join(dest_dir, f"{species_key}.annotation.gtf.gz")

        success, error = _download_file(genome_url, genome_gz, progress_callback)
        if not success:
            return False, None, f"Genome download failed: {error}"
        success, error = _decompress_gz(genome_gz, genome_fa)
        if not success:
            return False, None, f"Genome decompression failed: {error}"
        os.remove(genome_gz)

        success, error = _download_file(gtf_url, gtf_gz, progress_callback)
        if not success:
            return False, None, f"GTF download failed: {error}"
        success, error = _decompress_gz(gtf_gz, gtf_path_final)
        if not success:
            return False, None, f"GTF decompression failed: {error}"
        os.remove(gtf_gz)

        return True, (genome_fa, gtf_path_final), f"Downloaded genome from {genome_url} and annotation from {gtf_url}"

    else:
        return False, None, f"Unknown reference source '{entry['source']}'."


def extract_transcripts_with_gffread(genome_fasta_path, gtf_path, dest_fasta_path):
    """
    Extract transcript (cDNA) sequences from a genome FASTA + GTF/GFF
    annotation using gffread.

    Returns (success: bool, message: str).
    """
    if shutil.which("gffread") is None:
        return False, (
            "gffread is not installed on this system. It's required to "
            "extract transcript sequences from a genome + annotation file. "
            "(It's included in the project's Dockerfile.)"
        )

    cmd = ["gffread", "-w", dest_fasta_path, "-g", genome_fasta_path, gtf_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=1800)
        return True, "Transcript sequences extracted successfully."
    except subprocess.CalledProcessError as e:
        return False, f"gffread failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, "gffread timed out after 30 minutes."


def _parse_fasta_sequences(fasta_path):
    """
    Minimal, dependency-free FASTA parser. Returns a dict mapping
    sequence ID (the first whitespace-delimited token after '>') to the
    full concatenated sequence string.
    """
    sequences = {}
    current_id = None
    current_chunks = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_chunks)
                current_id = line[1:].split()[0]
                current_chunks = []
            else:
                current_chunks.append(line.strip())

    if current_id is not None:
        sequences[current_id] = "".join(current_chunks)

    return sequences


_COMPLEMENT_TABLE = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _reverse_complement(seq):
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def _parse_gtf_attribute(attribute_string, key):
    """Extract a single 'key "value"' pair from a GTF attribute column."""
    match = re.search(rf'{key} "([^"]+)"', attribute_string)
    return match.group(1) if match else None


def extract_gene_level_transcripts(genome_fasta_path, gtf_path, dest_fasta_path):
    """
    Extract one sequence per gene directly from a genome FASTA using
    gene-level coordinates from a GTF file, bypassing gffread entirely.

    Returns (success: bool, message: str).
    """
    try:
        genome_seqs = _parse_fasta_sequences(genome_fasta_path)
    except OSError as e:
        return False, f"Could not read genome FASTA: {e}"

    if not genome_seqs:
        return False, "Genome FASTA appears to be empty or unreadable."

    n_extracted = 0
    n_skipped = 0

    try:
        with open(gtf_path, "r") as gtf_file, open(dest_fasta_path, "w") as out_file:
            for line in gtf_file:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue

                seqid, _, _, start, end, _, strand, _, attributes = fields
                gene_id = _parse_gtf_attribute(attributes, "gene_id")
                if not gene_id:
                    n_skipped += 1
                    continue

                if seqid not in genome_seqs:
                    n_skipped += 1
                    continue

                start_idx, end_idx = int(start) - 1, int(end)
                gene_seq = genome_seqs[seqid][start_idx:end_idx]

                if strand == "-":
                    gene_seq = _reverse_complement(gene_seq)

                if not gene_seq:
                    n_skipped += 1
                    continue

                out_file.write(f">{gene_id}\n{gene_seq}\n")
                n_extracted += 1
    except OSError as e:
        return False, f"Could not read GTF or write output: {e}"

    if n_extracted == 0:
        return False, (
            "No gene sequences could be extracted -- check that the GTF "
            "file's chromosome/contig names match the genome FASTA's "
            "sequence IDs."
        )

    message = f"Extracted {n_extracted} gene sequences directly from genome coordinates."
    if n_skipped:
        message += f" ({n_skipped} GTF gene records were skipped due to missing IDs or coordinate mismatches.)"
    return True, message


def validate_fasta_file(filepath):
    """
    Basic sanity check that an uploaded/downloaded file actually looks
    like a FASTA file (starts with '>'), to catch cases where a user
    uploads the wrong file type entirely.
    """
    try:
        with open(filepath, "r", errors="ignore") as f:
            first_line = f.readline().strip()
        return first_line.startswith(">")
    except OSError:
        return False


def validate_annotation_file(filepath):
    """
    Basic sanity check that an uploaded annotation file looks like a
    GTF/GFF file (has tab-separated columns, and isn't just a FASTA).
    """
    try:
        with open(filepath, "r", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                return len(line.split("\t")) >= 8
        return False
    except OSError:
        return False


def count_transcripts_in_fasta(filepath):
    """Count how many sequences (transcripts/genes) a FASTA file contains."""
    try:
        count = 0
        with open(filepath, "r", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    count += 1
        return count
    except OSError:
        return 0


# ---------------------------------------------------------------------------
# tx2gene mapping extraction (for tximport gene-level collapsing)
# ---------------------------------------------------------------------------

def extract_tx2gene_from_ensembl_fasta(fasta_path):
    """
    Parse a transcript-to-gene mapping directly from an Ensembl-style
    cDNA FASTA's header lines.

    Returns a dict {transcript_id: gene_id}.
    """
    tx2gene = {}
    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            tx_id = header.split()[0]
            match = re.search(r"gene:(\S+)", header)
            if match:
                tx2gene[tx_id] = match.group(1)
    return tx2gene


def extract_gene_symbol_map_from_ensembl_fasta(fasta_path):
    """
    Parse a gene_id -> gene_symbol mapping directly from an Ensembl-style
    cDNA FASTA's header lines.

    Returns a dict {gene_id: gene_symbol}.
    """
    gene_symbol = {}
    gene_ids_seen = set()
    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            gene_match = re.search(r"gene:(\S+)", header)
            if not gene_match:
                continue
            gene_id = gene_match.group(1)
            gene_ids_seen.add(gene_id)
            symbol_match = re.search(r"gene_symbol:(\S+)", header)
            if symbol_match:
                gene_symbol[gene_id] = symbol_match.group(1)
    for gene_id in gene_ids_seen:
        gene_symbol.setdefault(gene_id, gene_id)
    return gene_symbol


def extract_gene_symbol_map_from_gtf(gtf_path):
    """
    Parse a gene_id -> gene_symbol mapping from a GTF file's feature
    lines, using the "gene_name" attribute.

    Returns a dict {gene_id: gene_symbol}.
    """
    gene_symbol = {}
    gene_ids_seen = set()
    with open(gtf_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attributes = fields[8]
            gene_id = _parse_gtf_attribute(attributes, "gene_id")
            if not gene_id:
                continue
            gene_ids_seen.add(gene_id)
            gene_name = _parse_gtf_attribute(attributes, "gene_name")
            if gene_name:
                gene_symbol[gene_id] = gene_name
    for gene_id in gene_ids_seen:
        gene_symbol.setdefault(gene_id, gene_id)
    return gene_symbol


def save_gene_symbol_map_csv(gene_symbol_map, dest_path):
    """
    Write a gene_id -> gene_symbol mapping dict to a CSV with columns
    "gene_id", "gene_name".
    """
    import csv
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_id", "gene_name"])
        for gene_id, gene_name in sorted(gene_symbol_map.items()):
            writer.writerow([gene_id, gene_name])
    return dest_path


def extract_tx2gene_from_gtf(gtf_path):
    """
    Parse a transcript-to-gene mapping from a GTF file's 'transcript'
    feature lines.

    Returns a dict {transcript_id: gene_id}.
    """
    tx2gene = {}
    with open(gtf_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            gene_match = re.search(r'gene_id "([^"]+)"', attrs)
            tx_match = re.search(r'transcript_id "([^"]+)"', attrs)
            if gene_match and tx_match:
                tx2gene[tx_match.group(1)] = gene_match.group(1)
    return tx2gene


def build_identity_tx2gene(fasta_path):
    """
    Build an identity tx2gene mapping (transcript_id == gene_id) from a
    FASTA's header lines.
    """
    tx2gene = {}
    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                gene_id = line[1:].strip().split()[0]
                tx2gene[gene_id] = gene_id
    return tx2gene


def save_tx2gene_csv(tx2gene, dest_path):
    """
    Write a tx2gene mapping dict to a CSV with columns TXNAME, GENEID.
    """
    import csv

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["TXNAME", "GENEID"])
        for tx_id, gene_id in tx2gene.items():
            writer.writerow([tx_id, gene_id])
    return dest_path
