"""
gene_id_mapper.py

Converts gene IDs between different identifier namespaces (Ensembl,
Entrez, RefSeq, FlyBase, WormBase, gene symbol, etc.) for the
Differential Expression workspace's "Gene ID -> Gene Name Mapping"
panel.

Why this exists (rather than only relying on the fast, no-R-required
parsing already done in reference_manager.py's
extract_gene_symbol_map_from_ensembl_fasta / extract_gene_symbol_map_from_gtf):
those functions can only produce a human-readable name when Ensembl (or
an uploaded GTF) already embeds one directly in the file, which isn't
universal across every species/release -- many genes (and entire
less-curated organisms) only carry a systematic/locus identifier there,
with no separate "gene_symbol"/"gene_name" field at all. This module
does a proper annotation database lookup instead, using Bioconductor's
clusterProfiler::bitr() ("Biological ID TRanslator") function via an
Rscript subprocess -- the same general approach deseq2_manager.py
already uses for running DESeq2 itself. This also generalizes to letting
the user pick ANY target namespace (e.g. Entrez ID instead of Symbol),
not just "fill in the gaps" for symbols specifically.

--- Requires ---
An R installation (like deseq2_manager.py's DESeq2 step) with the
`clusterProfiler` package AND the specific Bioconductor "OrgDb"
annotation package for the species being converted (e.g. `org.Hs.eg.db`
for human) -- see ORGDB_PACKAGES below for the exact package expected
per preset species. Errors from a missing R package are surfaced back
to the UI verbatim from R's own error message, since we can't detect a
missing R package without actually attempting to load it (same
limitation as deseq2_manager.py's Rscript/DESeq2 check).

--- A note on "SYMBOL" not being universal ---
Every OrgDb package supports different keytypes(), and a few of this
app's preset species deviate from the "obvious" choice:
  - org.Sc.sgd.db (yeast) has NO "SYMBOL" keytype at all -- Bioconductor
    ships it with "GENENAME" holding what every other OrgDb package
    calls SYMBOL (e.g. "GAL4"), and a separate "COMMON"/"ORF" pair for
    other name types. SYMBOL_KEYTYPE_OVERRIDES below routes yeast
    conversions to GENENAME automatically so "convert to a readable
    name" does the right thing without the user needing to know this
    quirk up front.
  - org.Dm.eg.db (fly) and org.Ce.eg.db (worm) are Entrez-based packages
    that ALSO support keying by their species' native community
    database ID (FLYBASE / WORMBASE respectively) -- which is what this
    app's Ensembl-sourced references for those two species actually use
    as "gene_id" natively (Ensembl imports FlyBase/WormBase gene IDs
    directly rather than minting its own ENS-prefixed ID for these
    species). detect_id_type() below recognizes both natively.
"""
import csv
import json
import os
import re
import shutil
import subprocess

import pandas as pd


# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

def bitr_tools_available():
    """
    Check whether Rscript is available. Whether clusterProfiler and the
    specific OrgDb package are installed can only be confirmed by
    actually attempting to load them in the R script itself (same
    pattern as deseq2_manager.py's deseq2_tools_available/DESeq2 check)
    -- a missing package produces a clear R-side error message that's
    passed back to the caller verbatim rather than us trying to
    pre-validate every possible package.
    """
    return shutil.which("Rscript") is not None


# ---------------------------------------------------------------------------
# Species -> Bioconductor OrgDb annotation package
# ---------------------------------------------------------------------------

# Maps this app's existing reference_manager.REFERENCE_CATALOG species
# keys 1:1 to their Bioconductor "OrgDb" annotation package, so a
# project's already-selected reference species can be reused directly
# here without asking the user to specify their organism a second time.
ORGDB_PACKAGES = {
    "human": "org.Hs.eg.db",
    "mouse": "org.Mm.eg.db",
    "yeast": "org.Sc.sgd.db",
    "drosophila": "org.Dm.eg.db",
    "celegans": "org.Ce.eg.db",
    "zebrafish": "org.Dr.eg.db",
    "ecoli": "org.EcK12.eg.db",
}

# Which OrgDb keytype should be used when the user just wants "a
# readable name"/SYMBOL, without needing to know each package's exact
# naming quirks. Every preset species uses the standard "SYMBOL"
# keytype EXCEPT yeast -- see the module docstring's note on
# org.Sc.sgd.db above for why GENENAME is the correct substitute there.
SYMBOL_KEYTYPE_OVERRIDES = {
    "yeast": "GENENAME",
}


def symbol_keytype_for_species(species_key):
    """The OrgDb keytype that represents a readable gene symbol/name for
    this species (almost always "SYMBOL"; see SYMBOL_KEYTYPE_OVERRIDES)."""
    return SYMBOL_KEYTYPE_OVERRIDES.get(species_key, "SYMBOL")


def orgdb_species_choices(reference_manager_module):
    """
    Build {species_key: display_label} for every preset species that
    has a known OrgDb package, for use in a manual species-picker
    dropdown (e.g. when a project used a custom reference with no
    automatic species association). Reuses the existing human-readable
    labels already defined in reference_manager.REFERENCE_CATALOG rather
    than duplicating species display names in this module.
    """
    labels = {}
    for species_key, orgdb_pkg in ORGDB_PACKAGES.items():
        catalog_entry = reference_manager_module.REFERENCE_CATALOG.get(species_key, {})
        base_label = catalog_entry.get("label", species_key)
        labels[species_key] = f"{base_label} — {orgdb_pkg}"
    return labels


# Common Bioconductor "keytypes" offered in the UI for the FROM/TO ID
# type pickers, covering both the standard vertebrate set (ENSEMBL,
# ENTREZID, REFSEQ, UNIPROT, SYMBOL, GENENAME) and the native
# community-database ID types this app's fly/worm/yeast presets
# actually use as their raw gene_id (FLYBASE, WORMBASE, ORF). Not every
# keytype is available for every OrgDb package -- bitr() itself
# validates this at run time and returns a clear R-side error if an
# unsupported keytype is requested, rather than us trying to maintain a
# fully exhaustive/exact per-package list here.
COMMON_KEY_TYPES = [
    "SYMBOL", "GENENAME", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT",
    "ENTREZID", "REFSEQ", "UNIPROT", "ACCNUM",
    "FLYBASE", "WORMBASE", "ORF",
]


# ---------------------------------------------------------------------------
# ID type auto-detection
# ---------------------------------------------------------------------------

# Regex patterns for recognizing common gene ID namespaces directly from
# the ID strings themselves, checked in this order (most specific
# first) against a SAMPLE of a project's gene IDs. Ensembl gene ID
# prefixes are species-specific for vertebrates (ENSG for human,
# ENSMUSG for mouse, ENSDARG for zebrafish, ...) but all share the same
# "ENS" + optional 3-4 letter species code + "G" (gene)/"T" (transcript)
# structure, so one pattern covers every Ensembl vertebrate preset
# species this app supports. Several non-vertebrate preset species
# (fly, worm, yeast) use their native community database's own ID
# scheme instead, since Ensembl imports these directly from
# FlyBase/WormBase/SGD rather than minting an "ENS"-prefixed ID -- those
# get their own dedicated patterns below.
_ID_PATTERNS = [
    ("FLYBASE",      re.compile(r"^FBgn\d+$")),                       # Drosophila (FlyBase gene ID)
    ("WORMBASE",      re.compile(r"^WBGene\d+$")),                     # C. elegans (WormBase gene ID)
    ("ORF",          re.compile(r"^Y[A-P][LR]\d{3}[WC](-[A-Z])?$")),   # Yeast systematic ORF name (SGD)
    ("ENSEMBLTRANS", re.compile(r"^ENS[A-Z]*T\d+")),                   # ENST/ENSMUST/ENSDART...
    ("ENSEMBLPROT",  re.compile(r"^ENS[A-Z]*P\d+")),                   # ENSP/ENSMUSP/ENSDARP...
    ("ENSEMBL",      re.compile(r"^ENS[A-Z]*G\d+")),                   # ENSG/ENSMUSG/ENSDARG... (gene)
    ("REFSEQ",       re.compile(r"^[NXY][MPRC]_\d+")),                 # NM_/NP_/NR_/XM_/XP_/XR_/NC_...
    ("UNIPROT",      re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")),
    ("ENTREZID",     re.compile(r"^\d+$")),                            # Plain numeric Entrez Gene ID
]


def detect_id_type(gene_ids, sample_size=200):
    """
    Guess the identifier namespace of a list/Series of gene IDs by
    pattern-matching a sample of them against _ID_PATTERNS above.

    Returns a dict:
        {
            "detected_type": str,      # one of _ID_PATTERNS' names, or "SYMBOL"
            "match_fraction": float,   # 0.0-1.0, how many of the sampled IDs matched
            "example_ids": [str, ...], # a few example IDs, for the user to sanity-check
        }

    "SYMBOL" is returned when no database-ID pattern matches a strong
    majority of the sample -- i.e. these IDs don't look like any known
    database accession scheme, so they're most likely already
    human-readable gene symbols (e.g. "TP53", "ACTB"), which have no
    single consistent regex of their own.
    """
    ids = [str(g) for g in gene_ids if pd.notna(g)]
    if not ids:
        return {"detected_type": "SYMBOL", "match_fraction": 0.0, "example_ids": []}

    sample = ids[:sample_size]
    best_type, best_count = "SYMBOL", 0
    for type_name, pattern in _ID_PATTERNS:
        count = sum(1 for gid in sample if pattern.match(gid))
        if count > best_count:
            best_type, best_count = type_name, count

    match_fraction = (best_count / len(sample)) if sample else 0.0
    # Require a reasonably strong majority match before trusting a
    # database-ID guess over the SYMBOL fallback -- a handful of
    # coincidentally-numeric or short-alphanumeric IDs shouldn't
    # override "these already look like readable symbols".
    if match_fraction < 0.5:
        best_type = "SYMBOL"

    return {
        "detected_type": best_type,
        "match_fraction": match_fraction,
        "example_ids": ids[:5],
    }


# ---------------------------------------------------------------------------
# R script for bitr()-based ID conversion
# ---------------------------------------------------------------------------

# This script is written to a temp file and executed via Rscript,
# mirroring deseq2_manager.py's _DESEQ2_R_SCRIPT pattern: a single
# positional argument (a path to a JSON "job spec" file) rather than
# plain CLI args, since the gene ID list can be very large.
_BITR_R_SCRIPT = r'''
suppressMessages({
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
job_spec_path <- args[1]
job <- fromJSON(job_spec_path)

orgdb_package <- job$orgdb_package
if (!requireNamespace(orgdb_package, quietly = TRUE)) {
  stop(paste0(
    "The Bioconductor annotation package '", orgdb_package, "' is not ",
    "installed. Install it with: BiocManager::install('", orgdb_package, "')"
  ))
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("The 'clusterProfiler' package is not installed. Install it with: BiocManager::install('clusterProfiler')")
}
suppressMessages(library(clusterProfiler))

gene_ids <- unique(as.character(job$gene_ids))

result <- tryCatch({
  clusterProfiler::bitr(
    gene_ids,
    fromType = job$from_type,
    toType = job$to_type,
    OrgDb = orgdb_package,
    drop = TRUE
  )
}, error = function(e) {
  stop(paste0("bitr() conversion failed: ", conditionMessage(e)))
})

# bitr() can return more than one output-type row per input ID (e.g.
# multiple UNIPROT accessions for one gene) -- keep only the first
# match per input ID so the final mapping stays 1:1, which is what the
# volcano plot / results table display expects.
result <- result[!duplicated(result[[job$from_type]]), ]

colnames(result) <- c("gene_id", "converted_id")
write.csv(result, job$output_path, row.names = FALSE)

cat(sprintf(
  "Converted %d of %d unique input ID(s) from %s to %s using %s.\n",
  nrow(result), length(gene_ids), job$from_type, job$to_type, orgdb_package
))
'''


def run_bitr_conversion(gene_ids, from_type, to_type, orgdb_package, work_dir):
    """
    Convert a list of gene IDs from one identifier namespace to another
    using clusterProfiler::bitr(), via an Rscript subprocess.

    gene_ids: list/Series of gene ID strings to convert (e.g. every
        gene_id in the project's counts matrix)
    from_type, to_type: Bioconductor keytype strings (see
        COMMON_KEY_TYPES), e.g. "ENSEMBL" -> "SYMBOL"
    orgdb_package: the Bioconductor OrgDb package name for this species
        (see ORGDB_PACKAGES), e.g. "org.Hs.eg.db"
    work_dir: scratch directory for the temp R script + job spec JSON
        + output CSV

    Returns a dict:
        {
            "success": bool,
            "message": str,           # human-readable status/error
            "mapping": dict,          # {gene_id: converted_id_or_original_id}
                                       # -- covers EVERY input gene_id,
                                       # falling back to the original ID
                                       # for anything bitr() couldn't
                                       # match, so nothing is ever
                                       # silently dropped
            "n_converted": int,       # how many WERE successfully converted
            "n_total": int,           # total unique input IDs
        }
    """
    unique_ids = sorted(set(str(g) for g in gene_ids if pd.notna(g)))
    n_total = len(unique_ids)

    if not bitr_tools_available():
        return {
            "success": False,
            "message": (
                "Rscript was not found on this system. R with the "
                "clusterProfiler package (and the relevant OrgDb "
                "annotation package) needs to be installed in your "
                "environment before ID conversion can run."
            ),
            "mapping": {gid: gid for gid in unique_ids},
            "n_converted": 0,
            "n_total": n_total,
        }

    if n_total == 0:
        return {
            "success": False,
            "message": "No gene IDs were provided to convert.",
            "mapping": {},
            "n_converted": 0,
            "n_total": 0,
        }

    os.makedirs(work_dir, exist_ok=True)
    output_csv_path = os.path.join(work_dir, "bitr_output.csv")
    job_spec = {
        "gene_ids": unique_ids,
        "from_type": from_type,
        "to_type": to_type,
        "orgdb_package": orgdb_package,
        "output_path": output_csv_path,
    }
    job_spec_path = os.path.join(work_dir, "bitr_job_spec.json")
    with open(job_spec_path, "w") as f:
        json.dump(job_spec, f, indent=2)

    r_script_path = os.path.join(work_dir, "run_bitr.R")
    with open(r_script_path, "w") as f:
        f.write(_BITR_R_SCRIPT)

    cmd = ["Rscript", r_script_path, job_spec_path]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=600
        )
        log = result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "success": False,
            "message": f"ID conversion failed: {(e.stdout or '') + (e.stderr or '')}",
            "mapping": {gid: gid for gid in unique_ids},
            "n_converted": 0,
            "n_total": n_total,
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "message": "ID conversion timed out after 10 minutes.",
            "mapping": {gid: gid for gid in unique_ids},
            "n_converted": 0,
            "n_total": n_total,
        }

    if not os.path.exists(output_csv_path):
        return {
            "success": False,
            "message": f"Conversion script ran but produced no output. Log: {log}",
            "mapping": {gid: gid for gid in unique_ids},
            "n_converted": 0,
            "n_total": n_total,
        }

    result_df = pd.read_csv(output_csv_path)
    converted = dict(zip(result_df["gene_id"].astype(str), result_df["converted_id"].astype(str)))

    # Build the FULL mapping covering every input ID, falling back to
    # the original ID for anything bitr() couldn't match -- so callers
    # never need to handle "missing" entries separately, consistent
    # with reference_manager.py's extract_gene_symbol_map_from_*
    # fallback convention.
    mapping = {gid: converted.get(gid, gid) for gid in unique_ids}
    n_converted = len(converted)

    pct = (n_converted / n_total * 100) if n_total else 0.0
    message = (
        f"Converted {n_converted:,} of {n_total:,} gene ID(s) ({pct:.1f}%) "
        f"from {from_type} to {to_type}. The remaining "
        f"{n_total - n_converted:,} gene(s) could not be matched and will "
        f"display their original {from_type} ID."
    )

    return {
        "success": True,
        "message": message,
        "mapping": mapping,
        "n_converted": n_converted,
        "n_total": n_total,
    }
