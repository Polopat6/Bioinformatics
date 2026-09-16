"""
single_cell/whitelist_manager.py

Downloads and manages 10x Genomics (and other droplet-based) barcode
whitelist files -- the SHARED, admin-managed resource
chemistry_manager.py's own CHEMISTRY_CATALOG expects to find under
chemistry_manager.SHARED_WHITELISTS_ROOT ("data/shared_whitelists").

--- Why this module exists (2026-09-01) ---
A real, confirmed gap: chemistry_manager.py's own module docstring
states whitelist files are "a SHARED, admin-managed resource (like a
reference genome)" -- but unlike reference genomes (which have a real,
working ensure_shared_resource()-backed "Download & Prepare Reference"
button in Step 5) and unlike the eggNOG database (which has its own
admin-gated download section on the Setup & Deployment page), NOTHING
in this codebase actually downloads a whitelist file automatically.
Step 1's own UI, when a whitelist is missing, only shows: "Contact your
admin to install it" -- with no further guidance and no download
mechanism anywhere. A user hit this directly and in doing so
discovered a required whitelist (737K-august-2016.txt, needed for 10x
5' v1.1/v2) simply wasn't present, despite two OTHER whitelist files
already having been manually placed.

Also confirmed as a real, hard blocker for this project's own stated
goal of possibly bundling these files directly into the git repo:
GitHub hard-blocks any single file over 100MB without Git LFS, and
this project's own 3M-february-2018.txt is already 115.5MB uncompressed
-- committing it directly would break `git push` outright. Downloading
on demand (exactly mirroring the reference-genome and eggNOG-database
patterns already used elsewhere in this app) avoids this entirely and
keeps the git repo itself small, at the one-time cost of a network
fetch the first time each whitelist is actually needed.

--- Source: a single, purpose-built public mirror ---
noamteyssier/10x_whitelist_mirror (github.com) hosts all five files
covered by WHITELIST_CATALOG below, gzip-compressed, under consistent
filenames matching 10x's own official naming (as independently
confirmed via 10x Genomics' own public Knowledge Base article "What is
a barcode inclusion list (formerly barcode whitelist)?", which is the
authoritative source for which whitelist file belongs to which
chemistry -- see each catalog entry's own comment). This mirror is used
as the download source purely for convenience/uptime (a small,
dedicated, single-purpose repo is less likely to reorganize/rename its
files than piggybacking on a much larger, general-purpose bioinformatics
tool's own repo) -- the underlying barcode DATA itself originates from
10x Genomics/Cell Ranger regardless of which mirror serves it.

--- Downloaded-and-decompressed, not left gzipped (2026-09-01) ---
chemistry_manager.py's own whitelist_confirm_match() opens whitelist
files with a plain open() call, NOT gzip.open() -- so this module
always decompresses the downloaded .gz file to a plain-text file at
EXACTLY the filename chemistry_manager.CHEMISTRY_CATALOG's own
"whitelist_file" field expects (e.g. "737K-august-2016.txt", no .gz
suffix), matching what that existing code already assumes.

--- Atomic download (2026-09-01) ---
download_whitelist() downloads to a temporary path first and only
renames it into its final location after both the download AND the
decompression succeed -- so a network failure or truncated download
partway through can never leave a corrupt, partially-written whitelist
file sitting at the real expected path looking like it's ready for use.
"""
import gzip
import os
import shutil
import urllib.error
import urllib.request
import app_paths
# Must exactly match chemistry_manager.SHARED_WHITELISTS_ROOT -- kept as
# its own independent constant here (rather than importing
# chemistry_manager just for this one string) to avoid this module
# needing to import chemistry_manager at all; the two are kept in sync
# manually since both values are simple, rarely-changed path strings.
SHARED_WHITELISTS_ROOT = app_paths.data_path("shared_whitelists")

_MIRROR_BASE_URL = "https://raw.githubusercontent.com/noamteyssier/10x_whitelist_mirror/main"

# --- Catalog of every whitelist file chemistry_manager.py's own
# CHEMISTRY_CATALOG can reference, with a real, working download URL for
# each. Chemistry-key associations and source-file naming below are
# independently confirmed against 10x Genomics' own public Knowledge
# Base article "What is a barcode inclusion list (formerly barcode
# whitelist)?" (kb.10xgenomics.com) -- the authoritative source for
# exactly which whitelist file belongs to which chemistry generation.
WHITELIST_CATALOG = {
    "737K-april-2014_rc.txt": {
        "url": f"{_MIRROR_BASE_URL}/737K-april-2014_rc.txt.gz",
        "chemistry_keys": ["10x_3p_v1"],
        "approx_decompressed_mb": 11,
        "description": "10x Genomics 3' v1 (discontinued)",
    },
    "737K-august-2016.txt": {
        "url": f"{_MIRROR_BASE_URL}/737K-august-2016.txt.gz",
        "chemistry_keys": ["10x_3p_v2", "10x_5p_v2", "10x_5p_v3"],
        "approx_decompressed_mb": 12,
        "description": "10x Genomics 3' v2, and 5' v1.1/v2/v3 (same barcode set)",
    },
    "3M-february-2018.txt": {
        "url": f"{_MIRROR_BASE_URL}/3M-february-2018.txt.gz",
        "chemistry_keys": ["10x_3p_v3"],
        "approx_decompressed_mb": 116,
        "description": "10x Genomics 3' v3 / v3.1",
    },
    "3M-3pgex-may-2023.txt": {
        "url": f"{_MIRROR_BASE_URL}/3M-3pgex-may-2023.txt.gz",
        "chemistry_keys": ["10x_3p_v4"],
        "approx_decompressed_mb": 116,
        "description": "10x Genomics 3' v4 (newest 3' barcode set)",
    },
    "3M-5pgex-jan-2023.txt": {
        "url": f"{_MIRROR_BASE_URL}/3M-5pgex-jan-2023.txt.gz",
        "chemistry_keys": ["10x_5p_v4"],
        "approx_decompressed_mb": 58,
        "description": "10x Genomics 5' v4 (newest 5' barcode set)",
    },
}


def shared_whitelist_dir():
    return SHARED_WHITELISTS_ROOT


def whitelist_status(filename):
    """
    Check whether a single catalog whitelist file is already present.

    Returns a dict:
        {"present": bool, "path": str, "size_mb": float|None}
    """
    path = os.path.join(SHARED_WHITELISTS_ROOT, filename)
    if os.path.isfile(path):
        size_mb = round(os.path.getsize(path) / (1024 * 1024), 1)
        return {"present": True, "path": path, "size_mb": size_mb}
    return {"present": False, "path": path, "size_mb": None}


def get_all_whitelist_statuses():
    "Convenience wrapper: whitelist_status() for every entry in WHITELIST_CATALOG at once."
    return {filename: whitelist_status(filename) for filename in WHITELIST_CATALOG}


def download_whitelist(filename, dest_dir=None, chunk_size=1024 * 1024, progress_callback=None):
    """
    Download and decompress a single whitelist file from its catalog
    entry's URL, writing the final plain-text file to dest_dir (defaults
    to SHARED_WHITELISTS_ROOT) under EXACTLY `filename` -- matching what
    chemistry_manager.py's own CHEMISTRY_CATALOG "whitelist_file" field
    expects to find.

    Atomic: downloads + decompresses to a temporary path first, only
    renaming into the final location once both steps succeed -- a
    network failure or truncated download can never leave a corrupt
    file sitting at the real expected path.

    Returns (success: bool, message: str).
    """
    if filename not in WHITELIST_CATALOG:
        return False, f"'{filename}' is not a known whitelist file in WHITELIST_CATALOG."

    spec = WHITELIST_CATALOG[filename]
    dest_dir = dest_dir or SHARED_WHITELISTS_ROOT
    os.makedirs(dest_dir, exist_ok=True)

    final_path = os.path.join(dest_dir, filename)
    tmp_gz_path = os.path.join(dest_dir, f".{filename}.download.gz.tmp")
    tmp_txt_path = os.path.join(dest_dir, f".{filename}.decompress.tmp")

    if progress_callback:
        progress_callback(f"Downloading {filename}.gz ({spec['approx_decompressed_mb']}MB uncompressed)...")

    try:
        with urllib.request.urlopen(spec["url"], timeout=60) as response, open(tmp_gz_path, "wb") as out_file:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        if os.path.exists(tmp_gz_path):
            os.remove(tmp_gz_path)
        return False, f"Failed to download {filename}: {e}"

    if progress_callback:
        progress_callback(f"Decompressing {filename}...")

    try:
        with gzip.open(tmp_gz_path, "rb") as f_in, open(tmp_txt_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    except (OSError, gzip.BadGzipFile) as e:
        for p in (tmp_gz_path, tmp_txt_path):
            if os.path.exists(p):
                os.remove(p)
        return False, f"Failed to decompress {filename}: {e}"

    os.remove(tmp_gz_path)

    # --- Verify BEFORE promoting the temp file (2026-09-16) ---
    # The download/decompress above is atomic, but nothing verified the
    # CONTENT. A truncated-but-valid gzip stream, or a mirror that
    # reorganizes and serves an HTML 404 page that happens to gunzip,
    # would otherwise land at the real path and report success -- and
    # chemistry_manager.whitelist_confirm_match() opens these with a
    # plain open(), so a bad file degrades silently into "chemistry
    # could only be confirmed by read length" much later, far from the
    # real cause. approx_decompressed_mb was already in the catalog but
    # only used for a progress message; it is a real expectation, so
    # check against it.
    actual_mb = os.path.getsize(tmp_txt_path) / (1024 * 1024)
    expected_mb = spec["approx_decompressed_mb"]
    if not (expected_mb * 0.8 <= actual_mb <= expected_mb * 1.2):
        os.remove(tmp_txt_path)
        return False, (
            f"{filename} downloaded but looks wrong: {actual_mb:.1f}MB, "
            f"expected about {expected_mb}MB. The mirror may have changed "
            f"its files, or the download was truncated. Nothing was "
            f"installed -- the previous file (if any) is untouched."
        )

    # A barcode whitelist is one fixed-length ACGT barcode per line.
    # Anything else (HTML, an error page, a README) fails here.
    with open(tmp_txt_path) as f:
        first_line = f.readline().strip()
    if not first_line or set(first_line) - set("ACGTN"):
        os.remove(tmp_txt_path)
        return False, (
            f"{filename} does not look like a barcode list -- its first "
            f"line is {first_line[:40]!r}, which is not a plain ACGT "
            f"barcode. Nothing was installed."
        )

    os.replace(tmp_txt_path, final_path)  # atomic rename into its final, real location

    size_mb = round(os.path.getsize(final_path) / (1024 * 1024), 1)
    return True, f"Downloaded and installed {filename} ({size_mb}MB)."



