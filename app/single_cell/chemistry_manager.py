"""
single_cell/chemistry_manager.py

Defines the droplet-based single-cell chemistries this pipeline supports,
and auto-detects which one a given FASTQ pair was most likely generated
with -- ALWAYS surfaced to the user as a suggestion to confirm or override
via dropdown, never silently trusted.

Scope: droplet-based UMI methods only (10x Genomics all chemistry
generations, Drop-seq, inDrops, BD Rhapsody) -- plate-based (Smart-seq)
and combinatorial-barcoding (Parse/split-seq) methods are intentionally
not included.

Confirmed chemistry specs (barcode/UMI lengths, strand orientation) are
taken from STARsolo's own community-maintained chemistry reference table.

Whitelist files are a SHARED, admin-managed resource (like a reference
genome) -- placed once into SHARED_WHITELISTS_ROOT, reused by every
project using that chemistry.
"""
import gzip
import os
from collections import Counter

SHARED_WHITELISTS_ROOT = "data/shared_whitelists"

# --- Recommended R2 (cDNA) sequencing length, per chemistry (2026-08-17) ---
# CONFIRMED directly from 10x Genomics' own support KB article ("Can I
# sequence longer than 91-98bp on Read 2..."): "98 bases (5' v1) or 91
# bases (3' v3, v3.1, 5' v2 chemistry) are sufficient to align
# confidently to the transcriptome... shorter reads have increased
# chances of multi-mapping." This is the TARGET/recommended full
# sequencing length for confident unique mapping -- NOT the same thing
# as singlecell_trim_manager.py's length_required floor (the minimum
# length to keep a read AFTER trimming, below which it's discarded
# outright as unusable). The floor does not need to vary by chemistry
# (it's a "don't keep hopeless reads" safety net, not a mapping-
# confidence target) -- this field is shown to the user as CONTEXT
# alongside that floor slider, not used to change the floor's default.
# recommended_r2_length_confirmed=True only for the exact chemistries
# 10x's KB article explicitly names (v3/v3.1, 5' v2, 5' v1) -- every
# other chemistry's value here is a reasonable same-ballpark estimate
# (e.g. v4 shares v3's barcode/UMI structure per Broad Institute's WARP
# docs, so likely behaves similarly) but is explicitly NOT independently
# confirmed, and is labeled as such wherever displayed.
CHEMISTRY_CATALOG = {
    "10x_3p_v1": {
        "label": "10x Genomics 3' v1 (discontinued, R1=24bp)",
        "protocol_key": "10XV1", "cb_len": 14, "umi_len": 10, "strand": "forward",
        "r1_len": 24, "whitelist_file": "737K-april-2014_rc.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": False,
    },
    "10x_3p_v2": {
        "label": "10x Genomics 3' v2 (R1=26bp)",
        "protocol_key": "10XV2", "cb_len": 16, "umi_len": 10, "strand": "forward",
        "r1_len": 26, "whitelist_file": "737K-august-2016.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": False,
    },
    "10x_3p_v3": {
        "label": "10x Genomics 3' v3 / v3.1 (R1=28bp)",
        "protocol_key": "10XV3", "cb_len": 16, "umi_len": 12, "strand": "forward",
        "r1_len": 28, "whitelist_file": "3M-february-2018.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": True,
    },
    "10x_3p_v4": {
        "label": "10x Genomics 3' v4 (R1=28bp, newest barcode set)",
        "protocol_key": "10XV4", "cb_len": 16, "umi_len": 12, "strand": "forward",
        "r1_len": 28, "whitelist_file": "3M-3pgex-may-2023.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": False,
    },
    "10x_5p_v2": {
        "label": "10x Genomics 5' v1.1/v2 (R1=26bp)",
        "protocol_key": "10XV2_5P", "cb_len": 16, "umi_len": 10, "strand": "reverse",
        "r1_len": 26, "whitelist_file": "737K-august-2016.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": True,
    },
    "10x_5p_v3": {
        "label": "10x Genomics 5' v3 (R1=28bp)",
        "protocol_key": "10XV3_5P", "cb_len": 16, "umi_len": 12, "strand": "reverse",
        "r1_len": 28, "whitelist_file": "737K-august-2016.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": False,
    },
    "10x_5p_v4": {
        "label": "10x Genomics 5' v4 (R1=28bp, newest barcode set)",
        "protocol_key": "10XV4_5P", "cb_len": 16, "umi_len": 12, "strand": "reverse",
        "r1_len": 28, "whitelist_file": "3M-5pgex-jan-2023.txt", "whitelist_required": True,
        "recommended_r2_length": 91, "recommended_r2_length_confirmed": False,
    },
    "dropseq": {
        "label": "Drop-seq (R1=20bp, no fixed whitelist -- barcodes called computationally)",
        "protocol_key": "DROPSEQ", "cb_len": 12, "umi_len": 8, "strand": "forward",
        "r1_len": 20, "whitelist_file": None, "whitelist_required": False,
        "recommended_r2_length": None, "recommended_r2_length_confirmed": False,
    },
    "indrops": {
        "label": "inDrops (R1 length varies by version -- confirm manually)",
        "protocol_key": "INDROPS", "cb_len": None, "umi_len": 6, "strand": "forward",
        "r1_len": None, "whitelist_file": None, "whitelist_required": False,
        "recommended_r2_length": None, "recommended_r2_length_confirmed": False,
    },
    "bd_rhapsody": {
        "label": "BD Rhapsody (R1=61bp, split cell-label + UMI structure)",
        "protocol_key": "BDRHAPSODY", "cb_len": 27, "umi_len": 8, "strand": "forward",
        "r1_len": 61, "whitelist_file": None, "whitelist_required": False,
        "recommended_r2_length": None, "recommended_r2_length_confirmed": False,
    },
}

_R1_LEN_INDEX = {}
for _key, _spec in CHEMISTRY_CATALOG.items():
    if _spec["r1_len"] is not None:
        _R1_LEN_INDEX.setdefault(_spec["r1_len"], []).append(_key)


def shared_whitelist_path(chemistry_key):
    spec = CHEMISTRY_CATALOG.get(chemistry_key)
    if not spec or not spec.get("whitelist_file"):
        return None
    return os.path.join(SHARED_WHITELISTS_ROOT, spec["whitelist_file"])


def whitelist_available(chemistry_key):
    path = shared_whitelist_path(chemistry_key)
    return path is not None and os.path.isfile(path)


def _sample_r1_lengths(r1_fastq_path, n_reads=2000):
    opener = gzip.open if r1_fastq_path.endswith(".gz") else open
    lengths = Counter()
    with opener(r1_fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i >= n_reads * 4:
                break
            if i % 4 == 1:
                lengths[len(line.strip())] += 1
    return lengths


def dominant_read_length(fastq_path, n_reads=2000):
    "Return just the single most common read length observed in fastq_path (or None if unreadable) -- reused by sc_sra_manager.py to classify SRA-downloaded files by length."
    lengths = _sample_r1_lengths(fastq_path, n_reads=n_reads)
    if not lengths:
        return None
    return lengths.most_common(1)[0][0]


def detect_chemistry(r1_fastq_path, n_reads=2000):
    """
    Sample R1 and suggest the most likely chemistry(ies) by read length --
    a HEURISTIC, always surfaced for user confirmation, never silently
    trusted. Returns {"candidates": [...], "observed_r1_len": int|None,
    "confidence": "high"|"ambiguous"|"unknown", "message": str}.
    """
    if not os.path.isfile(r1_fastq_path):
        return {"candidates": [], "observed_r1_len": None, "confidence": "unknown", "message": f"R1 file not found: {r1_fastq_path}"}
    lengths = _sample_r1_lengths(r1_fastq_path, n_reads=n_reads)
    if not lengths:
        return {"candidates": [], "observed_r1_len": None, "confidence": "unknown", "message": "No reads could be sampled from this R1 file."}
    observed_len, _count = lengths.most_common(1)[0]
    candidates = _R1_LEN_INDEX.get(observed_len, [])
    if len(candidates) == 1:
        return {"candidates": candidates, "observed_r1_len": observed_len, "confidence": "high",
                "message": f"R1 reads are {observed_len}bp, matching {CHEMISTRY_CATALOG[candidates[0]]['label']} -- please confirm this looks right below before continuing."}
    if len(candidates) > 1:
        labels = ", ".join(CHEMISTRY_CATALOG[c]["label"] for c in candidates)
        return {"candidates": candidates, "observed_r1_len": observed_len, "confidence": "ambiguous",
                "message": f"R1 reads are {observed_len}bp, which matches more than one known chemistry ({labels}) -- these can't be told apart by read length alone. Please select the correct one below."}
    return {"candidates": [], "observed_r1_len": observed_len, "confidence": "unknown",
            "message": (f"R1 reads are {observed_len}bp, which doesn't match any chemistry this pipeline currently "
                        "supports (10x Genomics 3'/5' v1-v4, Drop-seq, BD Rhapsody). If this is a plate-based method "
                        "like Smart-seq or a combinatorial method like Parse/split-seq, it isn't supported by this "
                        "pipeline yet -- please select manually below if you believe this is a supported chemistry regardless.")}


def whitelist_confirm_match(r1_fastq_path, chemistry_key, n_reads=2000, min_match_fraction=0.25):
    """
    Stronger confirmation than length-only: if this chemistry's whitelist
    file is actually installed locally, sample R1 barcodes and check what
    fraction match. Returns (checked: bool, match_fraction: float|None,
    passed: bool|None) -- checked=False if no whitelist requirement or
    file not found; callers fall back to the length-only heuristic then.
    """
    spec = CHEMISTRY_CATALOG.get(chemistry_key)
    if not spec or not spec.get("whitelist_required"):
        return False, None, None
    wl_path = shared_whitelist_path(chemistry_key)
    if not wl_path or not os.path.isfile(wl_path):
        return False, None, None
    cb_len = spec["cb_len"]
    with open(wl_path) as f:
        whitelist = set(line.strip() for line in f if line.strip())
    opener = gzip.open if r1_fastq_path.endswith(".gz") else open
    total, matches = 0, 0
    with opener(r1_fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i >= n_reads * 4:
                break
            if i % 4 == 1:
                barcode = line.strip()[:cb_len]
                total += 1
                if barcode in whitelist:
                    matches += 1
    if total == 0:
        return True, 0.0, False
    fraction = matches / total
    return True, fraction, fraction >= min_match_fraction
