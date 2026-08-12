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
    # certifi isn't installed — fall back to the default context. This
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
        # Confirmed present on Ensembl's main FTP site (not the separate
        # Ensembl Metazoa/Genomes portal), so the same current_fasta/
        # current_gtf symlink pattern used for human/mouse/yeast applies
        # unchanged. The exact "toplevel" vs. "primary_assembly" naming
        # below is inferred from Drosophila's compact, well-finished
        # genome structure (similar to yeast) rather than independently
        # confirmed against a live directory listing — if this doesn't
        # match on first real use, _resolve_ensembl_filename() will
        # return a clear "could not find a matching file" error rather
        # than silently downloading the wrong thing, so this is safe to
        # verify and correct on first actual use rather than a real risk.
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
        # Same note as Drosophila above re: pattern verification on
        # first real use — also confirmed on Ensembl's main FTP site,
        # same current_fasta/current_gtf symlink pattern applies.
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
        # Zebrafish has a chromosome-level assembly like human/mouse, so
        # it's distributed as "primary_assembly" (not "toplevel", which
        # is used for the smaller/less-finished genomes above).
        "genome_fasta_pattern": r"Danio_rerio\.GRCz11\.dna\.primary_assembly\.fa\.gz",
        "cdna_fasta_pattern": r"Danio_rerio\.GRCz11\.cdna\.all\.fa\.gz",
        "gtf_pattern": r"Danio_rerio\.GRCz11\.\d+\.gtf\.gz",
    },
    "ecoli": {
        "label": "E. coli (K-12 MG1655)",
        "source": "ncbi",
        # NCBI RefSeq assembly GCF_000005845.2 (ASM584v2) — a fixed,
        # permanent path that does not change over time.
        "ncbi_base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
        "genome_fasta_filename": "GCF_000005845.2_ASM584v2_genomic.fna.gz",
        "gtf_filename": "GCF_000005845.2_ASM584v2_genomic.gtf.gz",
        # NCBI does not provide a separate transcriptome-only FASTA for
        # this assembly by default, so for Salmon we extract transcript
        # sequences from the genome + GTF instead.
        "cdna_fasta_filename": None,
        # E. coli (and bacteria generally) have no introns, so NCBI's
        # GTF only contains "gene"/"CDS" features with no separate
        # transcript/mRNA parent line. gffread's parent-child ID
        # matching (built for eukaryotic gene models) fails on this
        # structure with "no valid ID found for GFF record". Since
        # gene == transcript for organisms with no splicing, we extract
        # sequences directly by gene coordinates instead of relying on
        # gffread at all.
        "no_introns": True,
    },
}

ENSEMBL_BASE_URL = "https://ftp.ensembl.org/pub"


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
    if not force and _resource_is_ready(resource_dir):
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
    # (e.g. in both an href and display text) — dedupe while preserving
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
                "Ensembl may have changed their file naming — please check "
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

    For species with a direct cDNA FASTA (e.g. human, mouse, yeast,
    Drosophila, C. elegans via Ensembl), this simply downloads it. For
    species without one (e.g. E. coli via NCBI, which does not publish a
    standalone transcriptome FASTA for this assembly), this automatically
    falls back to downloading the genome + GTF and extracting transcript
    sequences — so the UI never has to expose that distinction to the
    user or hand them a dead-end error message telling them to do
    something manually that the interface doesn't otherwise support for
    preset organisms.

    For the extraction step, organisms flagged "no_introns" in the
    catalog (bacteria/archaea) use direct gene-coordinate extraction
    instead of gffread, since gffread's eukaryotic gene->transcript->
    exon parent-child ID matching fails on NCBI's bacterial GTF format
    (which has no separate transcript/mRNA feature line to anchor to).

    Returns (success: bool, fasta_path_or_None, message: str).
    """
    entry = REFERENCE_CATALOG[species_key]

    # Direct path: a standalone cDNA FASTA is available for this species.
    if entry.get("cdna_fasta_filename") or entry.get("cdna_fasta_pattern"):
        return download_cdna_fasta(species_key, dest_dir, progress_callback)

    # Fallback path: no standalone transcriptome FASTA published for this
    # species/assembly — download genome + GTF and extract transcripts.
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


def download_genome_and_gtf(species_key, dest_dir, progress_callback=None):
    """
    Download the genome FASTA + GTF annotation for the given species key,
    for use with STAR (or for extracting transcripts via gffread if used
    as a Salmon fallback).

    Returns (success: bool, (genome_fasta_path, gtf_path) or None, message: str).
    """
    entry = REFERENCE_CATALOG[species_key]

    if entry["source"] == "ensembl":
        genome_dir_url = f"{ENSEMBL_BASE_URL}/current_fasta/{entry['species_dir']}/dna/"
        gtf_dir_url = f"{ENSEMBL_BASE_URL}/current_gtf/{entry['species_dir']}/"

        genome_filename = _resolve_ensembl_filename(genome_dir_url, entry["genome_fasta_pattern"])
        if not genome_filename:
            return False, None, (
                f"Could not find a matching genome FASTA file at {genome_dir_url}. "
                "Ensembl may have changed their file naming — please check "
                "manually or use the custom upload option instead."
            )
        gtf_filename = _resolve_ensembl_filename(gtf_dir_url, entry["gtf_pattern"])
        if not gtf_filename:
            return False, None, (
                f"Could not find a matching GTF file at {gtf_dir_url}. "
                "Ensembl may have changed their file naming — please check "
                "manually or use the custom upload option instead."
            )

        genome_url = genome_dir_url + genome_filename
        gtf_url = gtf_dir_url + gtf_filename

    elif entry["source"] == "ncbi":
        genome_url = f"{entry['ncbi_base_url']}/{entry['genome_fasta_filename']}"
        gtf_url = f"{entry['ncbi_base_url']}/{entry['gtf_filename']}"
    else:
        return False, None, f"Unknown reference source '{entry['source']}'."

    genome_gz = os.path.join(dest_dir, f"{species_key}.genome.fa.gz")
    genome_fa = os.path.join(dest_dir, f"{species_key}.genome.fa")
    gtf_gz = os.path.join(dest_dir, f"{species_key}.annotation.gtf.gz")
    gtf_path = os.path.join(dest_dir, f"{species_key}.annotation.gtf")

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
    success, error = _decompress_gz(gtf_gz, gtf_path)
    if not success:
        return False, None, f"GTF decompression failed: {error}"
    os.remove(gtf_gz)

    return True, (genome_fa, gtf_path), f"Downloaded genome from {genome_url} and annotation from {gtf_url}"


def extract_transcripts_with_gffread(genome_fasta_path, gtf_path, dest_fasta_path):
    """
    Extract transcript (cDNA) sequences from a genome FASTA + GTF/GFF
    annotation using gffread. This is used when only a genome + annotation
    are available (common for custom/non-model species) but a
    transcriptome FASTA is needed for Salmon.

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

    Used for direct gene-coordinate extraction on genomes small enough
    to comfortably fit in memory (e.g. bacterial genomes), avoiding a
    dependency on BioPython or samtools/pyfaidx for this simple case.
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

    This is the appropriate approach for organisms with no introns
    (bacteria, archaea) where "gene" and "transcript" are equivalent —
    it sidesteps gffread's eukaryotic gene->transcript->exon parent-
    child ID matching, which fails on NCBI's bacterial GTF format
    (features are "gene"/"CDS" only, with no separate transcript/mRNA
    line for gffread to anchor to).

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

                # GTF coordinates are 1-based, inclusive; Python slicing
                # is 0-based, exclusive on the end — convert accordingly.
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
            "No gene sequences could be extracted — check that the GTF "
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
#
# Salmon quantifies at the transcript level. To collapse transcript-level
# counts to gene-level (the standard, recommended input for DESeq2), a
# transcript-to-gene (tx2gene) mapping is needed. The correct source for
# this mapping differs depending on how the reference was obtained:
#
#   - Ensembl-downloaded cDNA (human/mouse/yeast/Drosophila/C. elegans
#     presets): the gene ID is embedded directly in each FASTA header,
#     so no separate GTF parsing is needed — see
#     extract_tx2gene_from_ensembl_fasta.
#
#   - No-intron organisms (e.g. E. coli): our own gene-level extraction
#     (extract_gene_level_transcripts, above) already writes one
#     sequence per gene using the gene ID as the header, so transcript
#     ID == gene ID by construction — see build_identity_tx2gene.
#
#   - Custom eukaryote uploads processed via gffread: gffread's output
#     FASTA headers don't reliably carry gene information, but the
#     original GTF's 'transcript' feature lines always do — see
#     extract_tx2gene_from_gtf.

def extract_tx2gene_from_ensembl_fasta(fasta_path):
    """
    Parse a transcript-to-gene mapping directly from an Ensembl-style
    cDNA FASTA's header lines. Ensembl's cDNA headers embed the parent
    gene ID directly, e.g.:

        >ENST00000456328.2 cdna chromosome:GRCh38:1:11869:14409:1
        gene:ENSG00000223972.5 gene_biotype:lncRNA ...

    This means no separate GTF download is needed just to build a
    tx2gene mapping for Ensembl-sourced references (human/mouse/yeast/
    Drosophila/C. elegans presets) — the mapping is already present in
    the FASTA we already downloaded for Salmon.

    Returns a dict {transcript_id: gene_id}. Returns an empty dict if no
    headers matched the expected pattern (e.g. a non-Ensembl FASTA was
    passed in) — callers should treat an empty result as "mapping not
    available" rather than assuming success.
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
    cDNA FASTA's header lines -- the same file already downloaded for
    Salmon and already parsed for tx2gene by
    extract_tx2gene_from_ensembl_fasta, so no extra download is needed
    just to get gene symbols for preset (human/mouse/yeast/Drosophila/
    C. elegans/zebrafish) references.

    Ensembl's cDNA headers embed both the parent gene ID and (when
    available) a human-readable symbol, e.g.:

        >ENST00000456328.2 cdna chromosome:GRCh38:1:11869:14409:1
        gene:ENSG00000223972.5 gene_biotype:lncRNA transcript_biotype:...
        gene_symbol:DDX11L1 description:...

    Not every gene has a populated "gene_symbol" field -- this is common
    for well-curated model organisms like human/mouse/zebrafish, but
    many genes even in those species (and more so in less-studied
    organisms) only have a systematic/locus identifier and no separate
    symbol. For any gene_id with no gene_symbol found on any of its
    transcripts, the gene_id itself is used as the "symbol" so every
    gene still gets a usable, human-facing label in the volcano plot
    and gene tables rather than being silently dropped from the
    mapping.

    Returns a dict {gene_id: gene_symbol}. Returns an empty dict if no
    headers matched the expected "gene:" pattern at all (e.g. a
    non-Ensembl FASTA was passed in) -- callers should treat an empty
    result as "mapping not available" rather than assuming success.
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
    lines, using the "gene_name" attribute -- the conventional GTF field
    for a gene's human-readable symbol (e.g. "TP53"), separate from its
    stable "gene_id" (e.g. "ENSG00000141510"). Used for custom
    eukaryotic reference uploads, where a GTF/GFF3 annotation is
    available but no Ensembl-style "gene_symbol:" FASTA header field
    exists (see extract_gene_symbol_map_from_ensembl_fasta for that
    preset-species path).

    Falls back to the gene_id itself for any gene whose GTF entry has
    no gene_name attribute (common for less-annotated or non-model
    organisms), so every gene still gets a usable label.

    Returns a dict {gene_id: gene_symbol}. Returns an empty dict if no
    'gene_id' attributes were found at all (e.g. a malformed or
    unexpected annotation format).
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
    "gene_id", "gene_name" -- matching the exact column names the
    Differential Expression workspace's gene-name-mapping loader already
    expects (see differential_expression_workspace.py's gene ID -> gene
    name CSV format), so this file can be read straight back in without
    any renaming step.
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
    feature lines (columns: transcript_id, gene_id attributes). This is
    the appropriate source for custom eukaryotic references where
    transcripts were extracted from a genome + GTF via gffread — gffread's
    output FASTA header doesn't reliably carry gene information, but the
    original GTF always does.

    Returns a dict {transcript_id: gene_id}. Returns an empty dict if no
    'transcript' feature lines with both attributes were found.
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
    FASTA's header lines. Appropriate for no-intron organisms (bacteria,
    archaea) where our own gene-level extraction (see
    extract_gene_level_transcripts) already writes one sequence per gene
    using the gene ID as the FASTA header — so "collapsing" is a no-op
    by construction, but we still return an explicit mapping for
    consistency with the other tx2gene sources.
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
    Write a tx2gene mapping dict to a CSV with columns TXNAME, GENEID —
    the column names tximport's tx2gene argument expects by convention.
    """
    import csv

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["TXNAME", "GENEID"])
        for tx_id, gene_id in tx2gene.items():
            writer.writerow([tx_id, gene_id])
    return dest_path
