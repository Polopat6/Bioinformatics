"""
single_cell/sc_cellqc_manager.py

Backend for Phase 2 of the Single-cell RNA-Seq pipeline -- Cell-level QC
(doublet detection, ambient RNA correction, and per-cell filtering),
invoked from singlecell_workspace.py's render_cell_qc().

Mirrors deseq2_manager.py's exact design pattern: a single R script
template, driven by a JSON "job spec" file, executed via a single
Rscript subprocess call using the same
subprocess.run(capture_output=True, text=True, check=True, timeout=...)
pattern as run_deseq2_analysis().

--- Design decisions (locked in 2026-08-17) ---

Doublet detection: BOTH scDblFinder (default) and DoubletFinder
(selectable alternative) are implemented, using FLAG-not-auto-remove
semantics. DoubletFinder runs its own minimal, self-contained Seurat
preprocessing pipeline purely as its own internal input -- NOT the same
thing as Phase 3's future real clustering workspace.

Ambient RNA correction: DecontX (default, per-cell estimate, no
external clustering required) or SoupX (selectable alternative, global
contamination fraction, needs the raw/unfiltered matrix).

Per-cell filtering thresholds: adaptive MAD-based (3 median-absolute-
deviations) is the default for mitochondrial %, total counts, and
detected genes.

--- Mitochondrial gene detection: real bug found and fixed (2026-08-17) ---
Mitochondrial gene identification previously relied ENTIRELY on
matching a "^MT-" gene-SYMBOL prefix -- if that silently found zero
matches for any reason (no Symbol column populated, non-standard gene
naming, OR the reference genuinely lacking a mitochondrial contig), a
misleadingly "healthy-looking" 0% mitochondrial threshold resulted with
no warning at all.

Fixed with get_mito_gene_ids_from_gtf() -- parses the reference GTF
directly by seqname (chromosome), independent of gene-symbol
annotation. The R script now computes gene-symbol-based AND GTF-derived
mito sets INDEPENDENTLY, reports each count separately, unions them for
actual QC use, and emits an explicit warning (both to the run log and a
dedicated mito_gene_diagnostic.json file) if the union is still empty.

--- Custom-reference mitochondrial gene resolution (2026-08-17, later
    same day) ---
For a CUSTOM (user-uploaded) reference whose GTF has zero identifiable
mitochondrial genes via direct seqname lookup, two additional resolution
paths are supported, surfaced via singlecell_workspace.py's Step 5 UI:
get_mito_gene_symbols_from_gtf() extracts mitochondrial gene SYMBOLS
from any reference GTF (typically a preset species' own already-
downloaded GTF); match_custom_genes_by_mito_symbol() searches a CUSTOM
reference's own GTF for genes whose gene_name matches one of those
symbols; find_gene_ids_matching_terms() supports fully manual gene
ID/symbol entry by scanning a GTF's ENTIRE gene_id space (not just
genes on a recognized mitochondrial contig), for the case where a
user manually knows which genes are mitochondrial in a reference whose
contig-naming convention wasn't recognized at all.

--- QC results package export (2026-08-17, later same day) ---
build_qc_summary_markdown() produces a self-contained plain-text/
markdown summary for bundling into a downloadable .zip alongside the
raw CSV/JSON outputs -- see singlecell_workspace.py's
_build_qc_package_zip().

--- Stale STAR index diagnostic (2026-08-17, later same day) ---
A real reported issue: reference_manager.verify_preset_reference_mito_content()
confirms a preset reference's genome FASTA + GTF DO include mitochondrial
content, yet re-running Cell-level QC on an existing sample still reports
zero mitochondrial genes even after re-indexing and re-aligning.

Root cause (see singlecell_workspace.py's own module docstring for the
full story): a separate real bug meant "force re-index" wasn't actually
rebuilding anything -- but even setting that aside, Cell-level QC ONLY
re-reads the EXISTING STARsolo output already on disk from a prior Step
6 alignment run; it does NOT re-align. STARsolo's own features.tsv (the
actual gene list making up the count matrix's rows) is generated once,
AT INDEX-BUILD TIME, baked into the STAR genome index itself -- not
re-read from the reference GTF at alignment time. So if a sample's
STAR index was built (once) from an OLDER version of the reference,
that sample's count matrix may never have included mitochondrial genes
as rows AT ALL, regardless of what the CURRENT reference GTF says.

diagnose_starsolo_matrix_for_mito() below provides a FAST (Python-only,
no R subprocess), sample-specific, INDEPENDENT check: does THIS sample's
own already-aligned features.tsv actually contain any mitochondrial gene
rows at all? Run automatically by singlecell_workspace.py's
render_cell_qc() the moment a sample is selected, BEFORE the (much
slower) full R-based Cell-level QC pipeline runs, so a stale-index
situation is caught and explained immediately with a specific,
actionable message (rebuild the index + re-run alignment, NOT just
re-run Cell-level QC again) rather than requiring a user to manually
grep features.tsv to figure out what's wrong.

--- STARsolo MTX loading ---
STARsolo's Solo.out/Gene/filtered/ directory mirrors Cell Ranger's own
10x-format MTX output layout, loaded via the standard Bioconductor
DropletUtils::read10xCounts() function.
"""
import gzip
import json
import os
import re
import shutil
import subprocess

import pandas as pd

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

def cellqc_tools_available():
    "Check whether Rscript is available."
    return shutil.which("Rscript") is not None


# ---------------------------------------------------------------------------
# Mitochondrial gene detection -- GTF-derived (2026-08-17)
# ---------------------------------------------------------------------------

_MITO_SEQNAME_ALIASES = {
    "mt", "chrm", "chrmt", "m", "mitochondrion", "mtdna",
    # Added: additional real-world naming conventions seen across
    # NCBI/Ensembl/manually-assembled non-model-organism references.
    "mito", "chrmito", "mt-genome", "mitochondrial",
    "mitochondrial_genome", "mitogenome",
}


def get_mito_gene_ids_from_gtf(gtf_path, max_lines=None):
    """
    Parse a GTF/GFF3 annotation file directly and return the set of
    gene_id values for every gene-level feature located on a
    mitochondrial contig (seqname matching _MITO_SEQNAME_ALIASES).

    Independent of gene SYMBOL/name annotation -- works even for
    references with no gene_name attributes, or gene_name values that
    don't follow the "MT-" naming convention. If this ALSO returns
    empty, that's strong direct evidence the reference genome/GTF
    itself has no mitochondrial contig at all.

    Returns a sorted list of gene_id strings (possibly empty).
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return []

    gene_ids = set()
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')

    with open(gtf_path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            seqname = fields[0].strip().lower()
            if seqname not in _MITO_SEQNAME_ALIASES:
                continue
            match = gene_id_pattern.search(fields[8])
            if match:
                gene_ids.add(match.group(1).strip())

    return sorted(gene_ids)


def get_mito_gene_symbols_from_gtf(gtf_path, max_lines=None):
    """
    Like get_mito_gene_ids_from_gtf(), but extracts gene SYMBOLS
    (gene_name attribute) rather than gene IDs, for genes located on a
    mitochondrial contig.

    Used to build a portable, cross-reference "known mitochondrial gene
    symbols" set from a species whose reference DOES clearly mark its
    mitochondrial genes (typically a preset species' own already-
    downloaded GTF) -- see match_custom_genes_by_mito_symbol() below.

    Returns a sorted list of gene_name strings (possibly empty).
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return []

    symbols = set()
    gene_name_pattern = re.compile(r'gene_name[\s=]+"?([^";]+)"?')

    with open(gtf_path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            seqname = fields[0].strip().lower()
            if seqname not in _MITO_SEQNAME_ALIASES:
                continue
            match = gene_name_pattern.search(fields[8])
            if match:
                symbols.add(match.group(1).strip())

    return sorted(symbols)


def match_custom_genes_by_mito_symbol(custom_gtf_path, known_mito_symbols):
    """
    Search a CUSTOM reference's own GTF for genes whose gene_name
    attribute case-insensitively matches one of known_mito_symbols (a
    list of mitochondrial gene symbols, typically obtained from a
    PRESET species' own GTF via get_mito_gene_symbols_from_gtf()).

    This is a LOWER-CONFIDENCE fallback than direct seqname-based
    detection (get_mito_gene_ids_from_gtf) -- it assumes the custom
    reference's gene-naming convention happens to overlap with the
    chosen preset species' naming, which is not guaranteed.

    Returns a sorted list of gene_id strings from the CUSTOM GTF whose
    gene_name matched (case-insensitively) any symbol in
    known_mito_symbols.
    """
    if not custom_gtf_path or not os.path.isfile(custom_gtf_path) or not known_mito_symbols:
        return []

    symbol_set_lower = {s.strip().lower() for s in known_mito_symbols}
    matched_gene_ids = set()
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')
    gene_name_pattern = re.compile(r'gene_name[\s=]+"?([^";]+)"?')

    with open(custom_gtf_path, "r", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            name_match = gene_name_pattern.search(fields[8])
            if not name_match:
                continue
            if name_match.group(1).strip().lower() not in symbol_set_lower:
                continue
            id_match = gene_id_pattern.search(fields[8])
            if id_match:
                matched_gene_ids.add(id_match.group(1).strip())

    return sorted(matched_gene_ids)


def find_gene_ids_matching_terms(gtf_path, terms):
    """
    Scan a GTF's ENTIRE gene_id attribute space (no chromosome/contig
    restriction) for exact, case-insensitive matches against a list of
    user-supplied terms.

    This is deliberately NOT the same as get_mito_gene_ids_from_gtf()
    (which only looks at genes already on a recognized mitochondrial
    contig) -- it exists specifically for the manual-entry mitochondrial
    resolution path in singlecell_workspace.py, where a user is manually
    identifying which genes (by their real gene_id) are mitochondrial in
    a reference whose contig-naming convention wasn't recognized at all.

    Returns a sorted list of gene_id strings from gtf_path that matched
    (case-insensitively) any entry in terms.
    """
    if not gtf_path or not os.path.isfile(gtf_path) or not terms:
        return []

    terms_lower = {t.strip().lower() for t in terms}
    matched = set()
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')

    with open(gtf_path, "r", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            match = gene_id_pattern.search(fields[8])
            if match and match.group(1).strip().lower() in terms_lower:
                matched.add(match.group(1).strip())

    return sorted(matched)

def get_gtf_contig_summary(gtf_path, max_contigs=None):
    """
    Scan a GTF/GFF3 file and return every distinct contig (seqname,
    column 1) it references, along with how many DISTINCT gene_id
    values are found on each -- powers the mitochondrial-contig PICKER
    fallback in singlecell_workspace.py's Step 5, shown when automatic
    mito detection (alias matching) finds zero mitochondrial genes.

    Counts gene_id occurrences from ANY feature-type line (gene,
    transcript, exon, CDS, etc.), not just explicit "gene" rows --
    manually curated/hand-added annotations (e.g. a custom
    mitochondrial contig with only transcript/exon lines and no
    separately generated "gene" line) are common for non-model-organism
    references and must not be silently invisible to this function.

    This only ever counts genes per CONTIG -- it never searches by
    gene name/identity, so it stays entirely separate from any
    DEG/volcano-plot gene-search utility elsewhere in this codebase.

    max_contigs: optional cap on how many contigs to return (already
        sorted by gene_count descending, so this keeps the highest-
        gene-count -- i.e. most "chromosome-like" -- contigs first).
        Leave unset to return the full inventory; the UI can layer its
        own "show top N, expand for all" behavior on top of the full
        list if desired.

    Returns a list of dicts, sorted by gene_count descending:
        [{"contig": "chr1", "gene_count": 4213}, ...]
    Returns an empty list if the file doesn't exist or has no
    recognizable feature lines at all.
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return []
    contig_gene_ids = {}
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')
    with open(gtf_path, "r", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            seqname = fields[0].strip()
            if not seqname:
                continue
            bucket = contig_gene_ids.setdefault(seqname, set())
            match = gene_id_pattern.search(fields[8])
            if match:
                bucket.add(match.group(1).strip())
    summary = [{"contig": contig, "gene_count": len(ids)} for contig, ids in contig_gene_ids.items()]
    summary.sort(key=lambda d: d["gene_count"], reverse=True)
    if max_contigs:
        summary = summary[:max_contigs]
    return summary



# Keyword/gene-name hints used ONLY to pre-highlight likely mitochondrial
# contigs in the picker UI -- never used to silently auto-select a contig
# without the user confirming. These are well-known, highly-conserved
# mitochondrial gene names/abbreviations that show up across essentially
# all animal (and most eukaryotic) mitochondrial genomes regardless of the
# organism or annotation pipeline used.
_MITO_KEYWORD_SUBSTRING = re.compile(r"mitochond", re.IGNORECASE)
_MITO_ID_SUBSTRING = re.compile(r"mito", re.IGNORECASE)
_MITO_GENE_NAME_HINTS = {
    "cox1", "cox2", "cox3", "co1", "co2", "co3", "cytb", "cob",
    "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6",
    "atp6", "atp8", "rrns", "rrnl", "12s", "16s", "trnf", "trnv",
}


def suggest_mito_contigs_by_keyword(gtf_path, min_hits=2):
    """
    Best-effort SUGGESTION (never a confident auto-detection) of which
    contig(s) are likely mitochondrial, by scanning every feature
    line's full attribute string (any feature type -- gene, transcript,
    exon, CDS, etc.) for:
      - the substring "mitochond" anywhere in the attribute string
        (catches gene_biotype/description/product annotations like
        "mitochondrially encoded ..." even when gene_name itself
        doesn't say so), OR
      - a gene_name OR gene_id matching a well-known, highly-conserved
        mitochondrial gene abbreviation (cox1, nd1, atp6, cytb, 12S/16S
        rRNA, etc.) -- checking gene_id in addition to gene_name
        matters for manually curated references that may have no
        gene_name attribute at all, only gene_id values like "COX1", OR
      - the substring "mito" (a shorter, narrower check than
        "mitochond" above) appearing specifically WITHIN a gene_id or
        transcript_id value -- catches manually curated/custom
        references where a curator has deliberately embedded "mito"
        into gene/transcript identifiers themselves (e.g. "genemito3",
        "rnamito4"), a common real-world convention for hand-annotated
        non-model-organism mitochondrial contigs. Deliberately checked
        ONLY against gene_id/transcript_id (not the full attribute
        string or free-text descriptions), since the short "mito"
        substring has real false-positive risk in free text -- e.g.
        "mitogen-activated protein kinase", "mitosis", "mitotic" all
        start with "mito" and are unrelated nuclear genes/processes,
        but this is far less likely to occur incidentally within a
        short structured identifier field.
    then tallying which CONTIG those hits fall on.

    This exists purely to pre-highlight candidates in the contig-picker
    UI (e.g. a "⭐ likely mitochondrial" badge next to a suggested
    contig) -- the user always makes the final selection themselves;
    this function's output is never used to auto-apply a mito contig
    choice on its own.

    min_hits: minimum number of keyword/gene-name hits required on a
        contig before it's included as a suggestion, to guard against
        one coincidental substring match on an unrelated contig (e.g.
        a nuclear pseudogene with "mitochondrial" in its description).

    Returns a list of dicts, sorted by hit_count descending:
        [{"contig": "NC_012920.1", "hit_count": 37}, ...]
    Returns an empty list if the file doesn't exist or nothing matched.
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return []
    hit_counts = {}
    gene_name_pattern = re.compile(r'gene_name[\s=]+"?([^";]+)"?')
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')
    transcript_id_pattern = re.compile(r'transcript_id[\s=]+"?([^";]+)"?')
    with open(gtf_path, "r", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            seqname = fields[0].strip()
            attributes = fields[8]
            hit = bool(_MITO_KEYWORD_SUBSTRING.search(attributes))
            if not hit:
                name_match = gene_name_pattern.search(attributes)
                if name_match and name_match.group(1).strip().lower() in _MITO_GENE_NAME_HINTS:
                    hit = True
            if not hit:
                id_match = gene_id_pattern.search(attributes)
                if id_match:
                    gid = id_match.group(1).strip()
                    if gid.lower() in _MITO_GENE_NAME_HINTS or _MITO_ID_SUBSTRING.search(gid):
                        hit = True
            if not hit:
                tx_match = transcript_id_pattern.search(attributes)
                if tx_match and _MITO_ID_SUBSTRING.search(tx_match.group(1).strip()):
                    hit = True
            if hit:
                hit_counts[seqname] = hit_counts.get(seqname, 0) + 1
    suggestions = [
        {"contig": contig, "hit_count": n}
        for contig, n in hit_counts.items()
        if n >= min_hits
    ]
    suggestions.sort(key=lambda d: d["hit_count"], reverse=True)
    return suggestions


def get_mito_gene_ids_from_gtf_by_contigs(gtf_path, contigs):
    """
    Like get_mito_gene_ids_from_gtf(), but restricted to an EXPLICIT set
    of user-CONFIRMED contig names from the contig-picker UI, rather
    than the built-in _MITO_SEQNAME_ALIASES set.

    Matches gene_id on ANY feature-type line (gene, transcript, exon,
    CDS, etc.), not just explicit "gene" rows -- so this correctly
    resolves genes for manually curated references whose mitochondrial
    annotation has no separately generated "gene" line, only
    transcript/exon/CDS lines.

    Used once a user has selected one or more contigs from
    get_gtf_contig_summary()'s full inventory (optionally pre-
    highlighted via suggest_mito_contigs_by_keyword()) as the real
    mitochondrial contig(s) for a reference whose naming convention
    wasn't recognized by automatic alias-based detection.

    contigs: list of exact contig/seqname strings as they appear in the
        GTF's column 1. These should come directly from
        get_gtf_contig_summary()'s own output (i.e. selected from a
        rendered list, not free-typed), so no case-normalization is
        applied here -- exact match only.

    Returns a sorted list of gene_id strings (possibly empty, if the
    selected contig(s) have no recognizable gene_id values at all).
    """
    if not gtf_path or not os.path.isfile(gtf_path) or not contigs:
        return []
    contig_set = set(contigs)
    gene_ids = set()
    gene_id_pattern = re.compile(r'gene_id[\s=]+"?([^";]+)"?')
    with open(gtf_path, "r", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[0].strip() not in contig_set:
                continue
            match = gene_id_pattern.search(fields[8])
            if match:
                gene_ids.add(match.group(1).strip())
    return sorted(gene_ids)

# ---------------------------------------------------------------------------
# Stale STAR index diagnostic (2026-08-17) -- see module docstring
# ---------------------------------------------------------------------------

def diagnose_starsolo_matrix_for_mito(filtered_matrix_dir, mito_gene_ids=None):
    """
    Fast, R-independent check of whether STARsolo's own features.tsv (the
    actual gene set baked into THIS sample's count matrix at alignment
    time) includes any mitochondrial genes at all -- checked via BOTH:
      1. gene-SYMBOL '^MT-' prefix matching on features.tsv's column 2
         (works for references with real gene symbols, e.g. most
         preset model organisms downloaded via Ensembl's direct GTF), and
      2. direct gene-ID membership against mito_gene_ids (the SAME
         resolved mito gene ID list -- see resolve_mito_gene_ids() above
         -- used everywhere else in this pipeline) on features.tsv's
         column 1 (works for references with NO gene symbols at all,
         e.g. any reference that went through the GFF3/gffread fallback
         download path, which strips gene_name genome-wide -- see
         reference_manager.py's backfill_gene_names_from_gff3() for the
         full story on why that happens).

    Checking BOTH independently, then unioning the results, mirrors the
    exact same "two independent detection methods, unioned" pattern
    already used successfully by the main R script's own mito detection
    (symbol_mito_set + gtf_mito_set) -- see this module's docstring.

    Without checking by ID as well as by symbol, this function would
    ALWAYS report zero mitochondrial genes (a false "stale index"
    signature) for any reference with no gene symbols at all --
    completely independent of whether that sample's index/alignment is
    actually stale or perfectly fresh -- which is exactly the bug this
    fix resolves.

    mito_gene_ids: the reference's resolved mitochondrial gene ID list
        (e.g. from resolve_mito_gene_ids()) -- pass None/empty if this
        isn't known/available yet (the check then falls back to
        symbol-only detection, same as this function's original
        behavior).

    Returns a dict:
        {
            "total_genes": int,
            "mito_genes_in_matrix": int,        # union of both methods
            "symbol_match_count": int,
            "id_match_count": int,
            "checked_via_id": bool,              # was an ID list available to check?
            "sample_mito_gene_names": [str, ...],  # up to 5, for display
            "likely_stale_index": bool,
        }
    Returns None if no features file could be found/read.

    likely_stale_index is ONLY set True when we actually HAD a positive
    expectation of mitochondrial genes for this reference (i.e.
    mito_gene_ids resolved at least one real gene ID from the CURRENT
    reference's GTF) but found NONE of them -- by symbol OR by ID -- as
    rows in this sample's own matrix. If no mito_gene_ids are known at
    all, an all-zero result is treated as "we don't have enough
    information to tell," NOT as evidence of staleness -- that
    genuinely different situation (no mito genes resolvable for this
    reference at all) is already separately surfaced by the main
    pipeline's own "no_mito_genes_found" diagnostic.
    """
    candidate_names = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    found_path = None
    for fname in candidate_names:
        candidate = os.path.join(filtered_matrix_dir, fname)
        if os.path.isfile(candidate):
            found_path = candidate
            break
    if found_path is None:
        return None

    mito_gene_id_set = set(mito_gene_ids) if mito_gene_ids else set()

    opener = gzip.open if found_path.endswith(".gz") else open
    total = 0
    symbol_hits = []
    id_hits = []
    try:
        with opener(found_path, "rt", errors="replace") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                total += 1
                gene_id = fields[0]
                gene_name = fields[1]
                if gene_name.lower().startswith("mt-"):
                    symbol_hits.append(gene_name)
                if gene_id in mito_gene_id_set:
                    id_hits.append(gene_id)
    except OSError:
        return None

    union_count = len(set(symbol_hits) | set(id_hits))
    likely_stale = bool(mito_gene_id_set) and total > 0 and union_count == 0

    return {
        "total_genes": total,
        "mito_genes_in_matrix": union_count,
        "symbol_match_count": len(symbol_hits),
        "id_match_count": len(id_hits),
        "checked_via_id": bool(mito_gene_id_set),
        "sample_mito_gene_names": (symbol_hits[:5] if symbol_hits else id_hits[:5]),
        "likely_stale_index": likely_stale,
    }



# ---------------------------------------------------------------------------
# Doublet detection method options
# ---------------------------------------------------------------------------

DOUBLET_METHOD_OPTIONS = {
    "scDblFinder": {
        "label": "scDblFinder (recommended -- fast, automated)",
        "implemented": True,
        "explanation": (
            "Simulates artificial doublets by combining pairs of real cells, then trains a "
            "classifier to score every real cell by how similar its expression profile looks to "
            "those simulated doublets. Fast, fully automated, and generally performs well across "
            "a wide range of dataset types and sizes."
        ),
    },
    "doublet_finder": {
        "label": "DoubletFinder (R/Seurat) -- slower, more manual, best raw accuracy",
        "implemented": True,
        "explanation": (
            "A well-established alternative that can achieve the best raw detection accuracy in "
            "independent benchmarks, at the cost of needing its own internal PCA/clustering step "
            "and a 'pK' parameter sweep -- both run automatically here, but take meaningfully "
            "longer than scDblFinder. Most useful with real ground-truth doublets to validate "
            "against, matching a prior published/consortium analysis, or a second opinion on an "
            "unusually complex dataset."
        ),
    },
}
DEFAULT_DOUBLET_METHOD = "scDblFinder"

DOUBLET_SIMULATION_MODES = {
    "random": {
        "label": "Random (recommended default)",
        "explanation": (
            "Simulated doublets are created by randomly pairing cells across the whole dataset. "
            "Tends to generalize better on complex/heterogeneous datasets, no clustering dependency."
        ),
    },
    "cluster": {
        "label": "Cluster-based",
        "explanation": (
            "Simulated doublets are created by pairing cells preferentially from DIFFERENT "
            "clusters. Requires a clustering assignment to already exist."
        ),
    },
}
DEFAULT_SIMULATION_MODE = "random"


def compute_expected_doublet_rate(n_cells):
    "Auto-compute a starting expected doublet rate (~0.8% per 1,000 cells loaded), returned as a fraction."
    return round((n_cells / 1000.0) * 0.008, 4)


# ---------------------------------------------------------------------------
# Ambient RNA correction method options
# ---------------------------------------------------------------------------

AMBIENT_METHOD_OPTIONS = {
    "decontx": {
        "label": "DecontX (recommended default)",
        "explanation": (
            "Part of the R/Bioconductor 'celda' package. Estimates ambient RNA contamination as a "
            "per-CELL fraction, and computes its own internal clustering automatically."
        ),
    },
    "soupx": {
        "label": "SoupX",
        "explanation": (
            "A mature, widely-used alternative. Estimates the ambient 'soup' profile from empty "
            "droplets in the raw (unfiltered) matrix, then applies a single GLOBAL contamination "
            "fraction to every cell."
        ),
    },
}
DEFAULT_AMBIENT_METHOD = "decontx"

DEFAULT_MAD_NMADS = 3
MITO_SOURCE_USER_SELECTED_CONTIGS = "user_selected_contigs"
QC_METRIC_LABELS = {
    "sum": "Total UMI counts per cell",
    "detected": "Genes detected per cell",
    "subsets_Mito_percent": "Mitochondrial % per cell",
}

DEFAULT_DOUBLETFINDER_N_PCS = 30
DEFAULT_DOUBLETFINDER_CLUSTER_RESOLUTION = 0.8
DOUBLETFINDER_PN = 0.25


# ---------------------------------------------------------------------------
# R script (job-spec driven)
# ---------------------------------------------------------------------------
_CELLQC_R_SCRIPT = r'''
suppressMessages({
  library(jsonlite)
  library(DropletUtils)
  library(scuttle)
})
args <- commandArgs(trailingOnly = TRUE)
job_spec_path <- args[1]
job <- fromJSON(job_spec_path)

sce <- read10xCounts(job$filtered_matrix_dir, col.names = TRUE)
cat(paste("Loaded matrix:", ncol(sce), "cells x", nrow(sce), "genes"), "\n")

# --- Mitochondrial gene detection: TWO independent methods ---
symbol_mito_set <- character(0)
if ("Symbol" %in% colnames(rowData(sce))) {
  symbol_mito_set <- rownames(sce)[grepl("^mt-", rowData(sce)$Symbol, ignore.case = TRUE)]
}
gtf_mito_set <- character(0)
if (!is.null(job$mito_gene_ids) && length(job$mito_gene_ids) > 0) {
  gtf_mito_set <- intersect(rownames(sce), job$mito_gene_ids)
}
mito_set <- union(symbol_mito_set, gtf_mito_set)

cat(paste0(
  "Mitochondrial gene detection: ", length(symbol_mito_set), " found via gene-symbol '^MT-' matching, ",
  length(gtf_mito_set), " found via GTF/override gene ID lookup, ",
  length(mito_set), " total (union) used for QC."
), "\n")

mito_diagnostic <- list(
  symbol_match_count = length(symbol_mito_set),
  gtf_derived_count = length(gtf_mito_set),
  union_count = length(mito_set),
  mito_source = job$mito_source,
  sample_gene_ids = if (length(mito_set) > 0) head(mito_set, 5) else list()
)
if (length(mito_set) == 0) {
  cat("⚠️ WARNING: ZERO mitochondrial genes identified by EITHER detection method.\n")
  mito_diagnostic$warning <- "no_mito_genes_found"
} else {
  mito_diagnostic$warning <- NULL
}
write(toJSON(mito_diagnostic, auto_unbox = TRUE), file.path(job$output_dir, "mito_gene_diagnostic.json"))

qc_metrics <- perCellQCMetrics(sce, subsets = list(Mito = mito_set))

low_counts <- isOutlier(qc_metrics$sum, nmads = job$nmads, type = "lower", log = TRUE)
low_genes <- isOutlier(qc_metrics$detected, nmads = job$nmads, type = "lower", log = TRUE)
high_mito <- isOutlier(qc_metrics$subsets_Mito_percent, nmads = job$nmads, type = "higher")

qc_thresholds <- list(
  sum_lower = as.numeric(attr(low_counts, "thresholds")["lower"]),
  detected_lower = as.numeric(attr(low_genes, "thresholds")["lower"]),
  mito_upper = as.numeric(attr(high_mito, "thresholds")["higher"])
)
write(toJSON(qc_thresholds, auto_unbox = TRUE), file.path(job$output_dir, "qc_thresholds.json"))

qc_df <- as.data.frame(qc_metrics)
qc_df$barcode <- colnames(sce)
qc_df$low_counts_flag <- low_counts
qc_df$low_genes_flag <- low_genes
qc_df$high_mito_flag <- high_mito
qc_df$adaptive_qc_fail <- low_counts | low_genes | high_mito
cat(paste("Cells flagged by adaptive QC (any metric):", sum(qc_df$adaptive_qc_fail), "of", nrow(qc_df)), "\n")

# --- Doublet detection ---
if (job$doublet_method == "scDblFinder") {
  suppressMessages(library(scDblFinder))
  set.seed(job$random_seed)
  if (job$simulation_mode == "cluster") {
    sce <- scDblFinder(sce, dbr = job$expected_doublet_rate, clusters = TRUE)
  } else {
    sce <- scDblFinder(sce, dbr = job$expected_doublet_rate)
  }
  qc_df$doublet_score <- sce$scDblFinder.score
  qc_df$predicted_doublet <- sce$scDblFinder.class == "doublet"
  cat(paste("Doublets flagged by scDblFinder:", sum(qc_df$predicted_doublet), "of", nrow(qc_df)), "\n")
} else if (job$doublet_method == "doublet_finder") {
  suppressMessages({
    library(Seurat)
    library(DoubletFinder)
  })
  set.seed(job$random_seed)

  seu <- CreateSeuratObject(counts = counts(sce))
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)

  n_pcs <- min(job$doubletfinder_n_pcs, ncol(seu) - 1)
  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = job$doubletfinder_cluster_resolution, verbose = FALSE)

  sweep_res <- paramSweep(seu, PCs = 1:n_pcs, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  write.csv(bcmvn, file.path(job$output_dir, "doubletfinder_pk_sweep.csv"), row.names = FALSE)

  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat(paste("DoubletFinder pK sweep complete -- selected pK =", optimal_pK), "\n")

  # --- NEW HARDENING #1: is the selected pK sitting at the edge of the
  # tested sweep range? This can indicate the sweep didn't actually
  # bracket a real optimum (the true best pK may lie outside what was
  # tested), rather than genuinely landing on the best value.
  all_pK_values <- as.numeric(as.character(bcmvn$pK))
  pk_range_min <- min(all_pK_values)
  pk_range_max <- max(all_pK_values)
  pk_at_edge <- isTRUE(all.equal(optimal_pK, pk_range_min)) || isTRUE(all.equal(optimal_pK, pk_range_max))
  if (pk_at_edge) {
    cat(paste0(
      "\u26a0\ufe0f WARNING: Selected pK (", optimal_pK, ") sits at the ",
      if (isTRUE(all.equal(optimal_pK, pk_range_min))) "MINIMUM" else "MAXIMUM",
      " edge of the tested sweep range [", pk_range_min, ", ", pk_range_max, "] -- ",
      "this can mean the true optimal pK lies outside what was tested, rather than ",
      "genuinely being the best value found."
    ), "\n")
  }

  annotations <- seu@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  nExp_poi <- round(job$expected_doublet_rate * ncol(seu))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

  seu <- doubletFinder(seu, PCs = 1:n_pcs, pN = job$doubletfinder_pn, pK = optimal_pK,
                        nExp = nExp_poi_adj, sct = FALSE)

  pann_col <- grep("^pANN", colnames(seu@meta.data), value = TRUE)[1]
  class_col <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)[1]

  # --- NEW HARDENING #2: DoubletFinder names its own output columns
  # dynamically (embedding pN/pK/nExp values into the column name
  # itself, e.g. "pANN_0.25_0.09_44") -- if the installed DoubletFinder
  # version's naming convention doesn't match what this grep() expects,
  # pann_col/class_col silently become NA rather than erroring, and
  # every downstream cell would get silently NA doublet scores. This
  # HARD FAILS with a clear, actionable message instead of letting
  # that propagate silently.
  if (is.na(pann_col) || is.na(class_col)) {
    stop(paste0(
      "DoubletFinder did not produce the expected output columns in Seurat's ",
      "metadata (looked for columns starting with 'pANN' and 'DF.classifications'). ",
      "Actual columns present: ", paste(colnames(seu@meta.data), collapse = ", "), ". ",
      "This usually means the installed DoubletFinder version uses a different ",
      "column-naming convention than this pipeline expects -- check DoubletFinder's ",
      "installed version/changelog."
    ))
  }

  match_idx <- match(qc_df$barcode, colnames(seu))

  # --- NEW HARDENING #3: Seurat is known to sometimes alter cell
  # barcode names during CreateSeuratObject() (e.g. "-" -> "_"
  # substitution in some versions/configs) -- if that happens here,
  # match() silently returns NA for every affected barcode, and those
  # cells get silently NA doublet scores/classifications with no error
  # at all. This computes and reports exactly how many/what fraction of
  # cells were affected, so a partial silent failure is visible rather
  # than hidden.
  n_cells_total <- length(match_idx)
  n_cells_unmatched <- sum(is.na(match_idx))
  pct_cells_unmatched <- round(100 * n_cells_unmatched / n_cells_total, 1)
  if (n_cells_unmatched > 0) {
    cat(paste0(
      "\u26a0\ufe0f WARNING: ", n_cells_unmatched, " of ", n_cells_total, " cells (",
      pct_cells_unmatched, "%) could not be matched between the original barcode list ",
      "and Seurat's own cell names -- these cells will have NA doublet scores/",
      "classifications. This can happen if Seurat modified barcode names during ",
      "CreateSeuratObject() (e.g. a '-' -> '_' substitution)."
    ), "\n")
  }

  qc_df$doublet_score <- seu@meta.data[[pann_col]][match_idx]
  qc_df$predicted_doublet <- seu@meta.data[[class_col]][match_idx] == "Doublet"
  qc_df$doubletfinder_cluster <- as.character(seu@meta.data$seurat_clusters)[match_idx]
  qc_df$doubletfinder_pK_used <- optimal_pK

  # --- Diagnostic JSON, read back by sc_cellqc_manager.py's
  # diagnose_doubletfinder_result() and rendered in the UI ---
  doubletfinder_diagnostic <- list(
    selected_pK = optimal_pK,
    pK_range_min = pk_range_min,
    pK_range_max = pk_range_max,
    pK_at_edge = pk_at_edge,
    n_cells_total = n_cells_total,
    n_cells_unmatched = n_cells_unmatched,
    pct_cells_unmatched = pct_cells_unmatched,
    n_expected_doublets = nExp_poi_adj,
    homotypic_proportion = round(homotypic_prop, 4)
  )
  write(toJSON(doubletfinder_diagnostic, auto_unbox = TRUE), file.path(job$output_dir, "doubletfinder_diagnostic.json"))

  cat(paste("Doublets flagged by DoubletFinder:", sum(qc_df$predicted_doublet, na.rm = TRUE), "of", nrow(qc_df)), "\n")
} else {
  stop(paste("Doublet detection method not implemented:", job$doublet_method))
}


# --- Ambient RNA correction ---
original_counts_mat <- counts(sce)
if (job$ambient_method == "decontx") {
  suppressMessages(library(celda))
  decon_result <- decontX(sce)
  qc_df$ambient_contamination <- decon_result$decontX_contamination
  corrected_counts <- round(decontXcounts(decon_result))
  Matrix::writeMM(as(corrected_counts, "CsparseMatrix"), file.path(job$output_dir, "corrected_counts.mtx"))
  cat(paste("DecontX mean per-cell contamination fraction:", round(mean(decon_result$decontX_contamination), 4)), "\n")} else if (job$ambient_method == "soupx") {
  suppressMessages(library(SoupX))
  raw_sce <- read10xCounts(job$raw_matrix_dir, col.names = TRUE)
  sc <- SoupChannel(counts(raw_sce), counts(sce))
  quick_clusters <- tryCatch({
    suppressMessages(library(scran))
    clusters <- quickCluster(sce)
    setNames(as.character(clusters), colnames(sce))
  }, error = function(e) {
    cat(paste("Note: could not compute clustering for SoupX via scran::quickCluster() --",
              conditionMessage(e), "-- SoupX requires the `scran` package to be installed."), "\n")
    NULL
  })
  if (is.null(quick_clusters)) {
    stop("SoupX ambient RNA correction requires the R package `scran` (for quickCluster()) to provide clustering information -- please install `scran` and try again.")
  }
  sc <- setClusters(sc, quick_clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  corrected_counts <- adjustCounts(sc)
  qc_df$ambient_contamination <- sc$metaData$rho[match(qc_df$barcode, rownames(sc$metaData))]
  Matrix::writeMM(as(corrected_counts, "CsparseMatrix"), file.path(job$output_dir, "corrected_counts.mtx"))
  cat(paste("SoupX global estimated contamination fraction (rho):", round(sc$fit$rhoEst, 4)), "\n")
} else {
  stop(paste("Ambient RNA correction method not implemented:", job$ambient_method))
}

gene_removed <- Matrix::rowSums(original_counts_mat) - Matrix::rowSums(corrected_counts)
gene_symbols <- if ("Symbol" %in% colnames(rowData(sce))) rowData(sce)$Symbol else rownames(sce)
gene_removed_df <- data.frame(
  gene_id = rownames(sce),
  symbol = gene_symbols,
  counts_removed = as.numeric(gene_removed)
)
gene_removed_df <- gene_removed_df[order(-gene_removed_df$counts_removed), ]
top_genes_df <- head(gene_removed_df, 20)
write.csv(top_genes_df, file.path(job$output_dir, "top_ambient_genes.csv"), row.names = FALSE)

write.csv(qc_df, file.path(job$output_dir, "cell_qc_metrics.csv"), row.names = FALSE)
cat("Cell-level QC completed successfully.\n")
'''


def _run_r_script(script_text, script_path, r_args, timeout=1800):
    "Write script_text to script_path and run it via Rscript with r_args."
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write(script_text)
    cmd = ["Rscript", script_path] + list(r_args)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, f"Cell-level QC failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, f"Cell-level QC timed out after {timeout} seconds."

def resolve_mito_gene_ids(mito_gtf_path=None, mito_gene_ids_override=None):
    """
    Shared mitochondrial-gene-ID resolution logic -- used by BOTH
    run_cellqc_analysis() (to build the R job spec) and
    diagnose_starsolo_matrix_for_mito() (the fast pre-flight check run
    the moment a sample is selected), so both paths always agree on
    exactly which gene IDs count as mitochondrial for a given
    reference, with no risk of the two independently drifting out of
    sync over time.

    Truthiness check (not "is not None") on mito_gene_ids_override --
    an empty list means "no real override was ever resolved," and must
    fall through to GTF auto-detection exactly like a genuinely absent
    override (None) would (see the earlier empty-list-override
    incident this same logic already fixes elsewhere in this module).
    """
    if mito_gene_ids_override:
        return list(mito_gene_ids_override)
    return get_mito_gene_ids_from_gtf(mito_gtf_path) if mito_gtf_path else []

def run_cellqc_analysis(filtered_matrix_dir, output_dir, work_dir,
                         doublet_method=DEFAULT_DOUBLET_METHOD,
                         simulation_mode=DEFAULT_SIMULATION_MODE,
                         expected_doublet_rate=None,
                         ambient_method=DEFAULT_AMBIENT_METHOD,
                         raw_matrix_dir=None,
                         nmads=DEFAULT_MAD_NMADS,
                         mito_gtf_path=None,
                         mito_gene_ids_override=None,
                         mito_source="gtf_auto_detect",
                         random_seed=42,
                         doubletfinder_n_pcs=DEFAULT_DOUBLETFINDER_N_PCS,
                         doubletfinder_cluster_resolution=DEFAULT_DOUBLETFINDER_CLUSTER_RESOLUTION,
                         timeout=None):
    """
    Run Phase 2 cell-level QC.

    mito_gtf_path: path to the reference GTF -- parsed here (via
        get_mito_gene_ids_from_gtf) to derive mitochondrial gene IDs
        directly, unless mito_gene_ids_override is provided instead.
    mito_gene_ids_override: an explicit, already-resolved list of
        mitochondrial gene IDs to use directly -- takes priority over
        mito_gtf_path if both are given.
    mito_source: a short label recorded in the mito_gene_diagnostic.json
        output describing HOW mito_gene_ids_override (if used) was
        obtained.

    Returns (success: bool, log: str).
    """
    if not cellqc_tools_available():
        return False, (
            "Rscript was not found on this system. R with the DropletUtils, scuttle, "
            "scDblFinder (or Seurat + DoubletFinder), and celda (DecontX) [and/or SoupX] "
            "packages needs to be installed in your environment before this step can run."
        )
    if not os.path.isdir(filtered_matrix_dir):
        return False, f"STARsolo filtered matrix directory not found: {filtered_matrix_dir}"
    if ambient_method == "soupx" and not (raw_matrix_dir and os.path.isdir(raw_matrix_dir)):
        return False, "SoupX requires the raw (unfiltered) STARsolo matrix directory, which was not found."

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    if expected_doublet_rate is None:
        expected_doublet_rate = 0.008

    if timeout is None:
        timeout = 5400 if doublet_method == "doublet_finder" else 1800

    if mito_gene_ids_override:
        # Truthiness check (not "is not None") -- an empty list means
        # "no real override was ever resolved," and must fall through
        # to GTF auto-detection exactly like a genuinely absent
        # override (None) would. A present-but-empty list previously
        # short-circuited auto-detection entirely, silently disabling
        # mitochondrial gene detection even on a fully correct,
        # properly-annotated reference. See this patch's own module
        # docstring for the full incident this fixes.
        mito_gene_ids = list(mito_gene_ids_override)
    else:
        mito_gene_ids = get_mito_gene_ids_from_gtf(mito_gtf_path) if mito_gtf_path else []

    job_spec = {
        "filtered_matrix_dir": os.path.abspath(filtered_matrix_dir),
        "raw_matrix_dir": os.path.abspath(raw_matrix_dir) if raw_matrix_dir else "",
        "output_dir": os.path.abspath(output_dir),
        "doublet_method": doublet_method,
        "simulation_mode": simulation_mode,
        "expected_doublet_rate": expected_doublet_rate,
        "ambient_method": ambient_method,
        "nmads": nmads,
        "mito_gene_ids": mito_gene_ids,
        "mito_source": mito_source,
        "random_seed": random_seed,
        "doubletfinder_n_pcs": doubletfinder_n_pcs,
        "doubletfinder_cluster_resolution": doubletfinder_cluster_resolution,
        "doubletfinder_pn": DOUBLETFINDER_PN,
    }
    job_spec_path = os.path.join(work_dir, "cellqc_job_spec.json")
    with open(job_spec_path, "w") as f:
        json.dump(job_spec, f, indent=2)

    script_path = os.path.join(work_dir, "run_cellqc.R")
    return _run_r_script(_CELLQC_R_SCRIPT, script_path, [job_spec_path], timeout=timeout)


# ---------------------------------------------------------------------------
# Reading back results
# ---------------------------------------------------------------------------

def read_cell_qc_metrics(output_dir):
    path = os.path.join(output_dir, "cell_qc_metrics.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def read_qc_thresholds(output_dir):
    path = os.path.join(output_dir, "qc_thresholds.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def read_mito_gene_diagnostic(output_dir):
    path = os.path.join(output_dir, "mito_gene_diagnostic.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def corrected_counts_path(output_dir):
    return os.path.join(output_dir, "corrected_counts.mtx")


def read_top_ambient_genes(output_dir):
    path = os.path.join(output_dir, "top_ambient_genes.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def read_doubletfinder_pk_sweep(output_dir):
    path = os.path.join(output_dir, "doubletfinder_pk_sweep.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)

def read_doubletfinder_diagnostic(output_dir):
    "Read the doubletfinder_diagnostic.json written by a DoubletFinder run, or None if not present (e.g. scDblFinder was used instead, or this project predates this diagnostic)."
    import os
    import json
    path = os.path.join(output_dir, "doubletfinder_diagnostic.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def diagnose_doubletfinder_result(diagnostic):
    """
    Build plain-language, tiered warnings from a DoubletFinder
    diagnostic dict (see read_doubletfinder_diagnostic) -- surfaces the
    real correctness risks hardened into the R script:
      1. pK selected at the edge of its tested sweep range (may not be
         a genuine optimum).
      2. A meaningful fraction of cells couldn't be matched back to
         Seurat's own cell names (silently NA doublet scores for those
         cells).

    Returns a dict:
        {
            "flagged": bool,          # True if ANY issue was found
            "messages": [str, ...],   # one plain-language message per issue found
        }
    Returns {"flagged": False, "messages": []} if diagnostic is None
    (nothing to check) or no issues were found.
    """
    if diagnostic is None:
        return {"flagged": False, "messages": []}

    messages = []

    if diagnostic.get("pK_at_edge"):
        edge = "minimum" if diagnostic.get("selected_pK") == diagnostic.get("pK_range_min") else "maximum"
        messages.append(
            f"⚠️ The selected pK ({diagnostic.get('selected_pK')}) sits at the **{edge}** edge of "
            f"the tested range [{diagnostic.get('pK_range_min')}, {diagnostic.get('pK_range_max')}] -- "
            "this can mean the true optimal pK lies outside what was tested, rather than genuinely "
            "being the best value. Results may still be usable, but this parameter choice is less "
            "certain than a pK selected comfortably within the tested range."
        )

    pct_unmatched = diagnostic.get("pct_cells_unmatched", 0)
    if pct_unmatched and pct_unmatched > 0:
        n_unmatched = diagnostic.get("n_cells_unmatched", 0)
        n_total = diagnostic.get("n_cells_total", 0)
        severity = "🔴" if pct_unmatched >= 5 else "🟡"
        messages.append(
            f"{severity} **{n_unmatched:,} of {n_total:,} cell(s) ({pct_unmatched}%)** could not be "
            "matched between the original barcode list and Seurat's own cell names during "
            "DoubletFinder's internal processing -- these cells have **no doublet score or "
            "classification at all** (silently NA), rather than a real result. This most often "
            "happens when Seurat's own barcode-parsing altered cell names slightly (e.g. a "
            "\"-\" to \"_\" substitution). "
            + ("This is a large enough fraction to treat DoubletFinder's results for this sample with real caution." if pct_unmatched >= 5 else "This is a small fraction, but still worth being aware of.")
        )

    return {"flagged": len(messages) > 0, "messages": messages}



# ---------------------------------------------------------------------------
# QC summary / tiered flagging
# ---------------------------------------------------------------------------

def summarize_cellqc_results(qc_df):
    n_cells = len(qc_df)
    if n_cells == 0:
        return None

    pct_qc_fail = round(100 * qc_df["adaptive_qc_fail"].sum() / n_cells, 1)
    pct_doublets = round(100 * qc_df["predicted_doublet"].sum() / n_cells, 1)
    mean_contam = round(100 * qc_df["ambient_contamination"].mean(), 1) if "ambient_contamination" in qc_df.columns else None

    def _tier(pct, caution_at, poor_at):
        if pct >= poor_at:
            return "poor"
        if pct >= caution_at:
            return "caution"
        return "good"

    qc_tier = _tier(pct_qc_fail, caution_at=15.0, poor_at=30.0)
    doublet_tier = _tier(pct_doublets, caution_at=12.0, poor_at=20.0)
    ambient_tier = _tier(mean_contam, caution_at=20.0, poor_at=35.0) if mean_contam is not None else "good"

    messages = {
        "adaptive_qc": {
            "poor": f"🔴 {pct_qc_fail}% of cells fail adaptive QC (low counts/genes or high mito%) -- unusually high; double-check sample/library quality before proceeding.",
            "caution": f"🟡 {pct_qc_fail}% of cells fail adaptive QC -- worth a look, but not unusual for many samples.",
            "good": f"🟢 Only {pct_qc_fail}% of cells fail adaptive QC -- healthy range.",
        }[qc_tier],
        "doublet": {
            "poor": f"🔴 {pct_doublets}% of cells flagged as likely doublets -- notably high; check whether cell loading concentration was higher than intended.",
            "caution": f"🟡 {pct_doublets}% of cells flagged as likely doublets -- a bit elevated, but within a plausible range depending on loading density.",
            "good": f"🟢 {pct_doublets}% of cells flagged as likely doublets -- consistent with typical 10x loading rates.",
        }[doublet_tier],
        "ambient": (
            {
                "poor": f"🔴 Mean ambient RNA contamination is {mean_contam}% -- quite high; results may benefit from closer inspection of soup-derived genes.",
                "caution": f"🟡 Mean ambient RNA contamination is {mean_contam}% -- moderate; normal for some tissue types.",
                "good": f"🟢 Mean ambient RNA contamination is {mean_contam}% -- low, healthy range.",
            }[ambient_tier] if mean_contam is not None else "ℹ️ Ambient contamination fraction not available for this run."
        ),
    }

    return {
        "n_cells": n_cells,
        "pct_adaptive_qc_fail": pct_qc_fail,
        "adaptive_qc_tier": qc_tier,
        "pct_doublets": pct_doublets,
        "doublet_tier": doublet_tier,
        "mean_ambient_contamination_pct": mean_contam,
        "ambient_tier": ambient_tier,
        "messages": messages,
    }


def build_qc_summary_markdown(sample_name, summary, mito_diagnostic, thresholds, doublet_method, ambient_method):
    """
    Build a self-contained markdown summary of a Cell-level QC run --
    for bundling into the downloadable QC package .zip.
    """
    lines = [f"# Cell-level QC Summary — {sample_name}", ""]
    if summary:
        lines += [
            f"**Cells analyzed:** {summary['n_cells']:,}",
            "",
            "## QC Flags",
            f"- {summary['messages']['adaptive_qc']}",
            f"- {summary['messages']['doublet']} (method: {doublet_method})",
            f"- {summary['messages']['ambient']} (method: {ambient_method})",
            "",
        ]
    if thresholds:
        lines += [
            "## Adaptive Thresholds Used",
            f"- Total counts > {thresholds.get('sum_lower', 0):.0f}",
            f"- Genes detected > {thresholds.get('detected_lower', 0):.0f}",
            f"- Mitochondrial % < {thresholds.get('mito_upper', 0):.1f}%",
            "",
        ]
    if mito_diagnostic:
        lines += ["## Mitochondrial Gene Detection"]
        if mito_diagnostic.get("warning") == "no_mito_genes_found":
            lines += [
                "⚠️ **Zero mitochondrial genes were identified by either detection method.** "
                "A 0% mitochondrial result for every cell is not a real biological finding -- "
                "it strongly suggests the reference genome/GTF used for alignment lacks the "
                "mitochondrial chromosome, uses an unrecognized naming convention, or this "
                "sample's alignment is stale (built before a reference fix). Mitochondrial QC "
                "filtering is not meaningful until this is resolved.",
                "",
            ]
        else:
            lines += [
                f"- Gene-symbol matches: {mito_diagnostic.get('symbol_match_count', 0)}",
                f"- GTF/override-derived matches: {mito_diagnostic.get('gtf_derived_count', 0)}",
                f"- Total (union) used for QC: {mito_diagnostic.get('union_count', 0)}",
                f"- Resolution method: {mito_diagnostic.get('mito_source', 'unknown')}",
                "",
            ]
    return "\n".join(lines)
