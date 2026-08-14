"""
ontology_manager.py

Backend for the Ontology Analysis workspace: runs Gene Ontology (GO),
KEGG, and Reactome enrichment via R's clusterProfiler (+ ReactomePA),
in three flavors:
  - ORA  (Over-Representation Analysis): enrichGO / enrichKEGG /
    enrichPathway -- takes a fixed list of "significant" genes plus a
    background "universe" and asks "are these genes enriched for any
    term more than chance?"
  - GSEA (Gene Set Enrichment Analysis): gseGO / gseKEGG / gsePathway --
    takes EVERY tested gene ranked by a score (log2FoldChange here) and
    asks "do this term's genes cluster toward one end of the ranking?"
  - compareCluster: runs ORA independently across multiple contrasts/
    gene lists at once and returns one combined table with a "Cluster"
    column, for side-by-side comparison.

Design decision -- ALL PLOTTING HAPPENS IN PYTHON (Plotly), not R: see
ontology_workspace.py's module docstring for the full rationale.

--- ORA: up- and down-regulated genes analyzed SEPARATELY (default) ---
Running ORA on a single gene list that MIXES up- and down-regulated
genes together is a real statistical/interpretive problem: a GO term
could be flagged as "enriched" purely because it happens to contain a
mix of genes going in both directions, without that term actually
representing a single, coherent, DIRECTIONAL biological signal (e.g. a
pathway that's genuinely activated, vs. one that's genuinely
suppressed, can get conflated into one ambiguous "changed somehow"
result). To address this, run_ora_analysis/build_ora_r_script default
to split_by_direction=True: two independent gene lists (up-regulated
only: padj < threshold AND log2FC >= lfc_threshold; down-regulated
only: padj < threshold AND log2FC <= -lfc_threshold) are each run
through enrichGO/enrichKEGG/enrichPathway SEPARATELY, producing two
distinct result files per database (suffixed "_up"/"_down" -- see
_result_file_path). If one direction has zero significant genes at the
chosen threshold, that direction's enrichment is skipped for that
database (with a clear note in the run log) WITHOUT aborting the other
direction's analysis or any other requested database -- each
direction/database combination is independently guarded in the
generated R script.

split_by_direction=False preserves the ORIGINAL combined behavior (one
mixed-direction gene list, one result file per database, no suffix) --
offered as an explicit, clearly-labeled non-default option in the UI
for anyone who specifically wants the old behavior or needs to compare
directly against a combined-gene-list result.

--- Gene set size bounds, GO simplification, permutations note ---
(Preserved from earlier design -- see min_gs_size/max_gs_size and
simplify_go parameters below for the same rationale as before: gene
sets that are too small produce unstable statistics; too large are
often generic/uninformative; GO's own redundant hierarchy benefits from
optional semantic-similarity-based simplification via GOSemSim, only
meaningful for a single GO sub-ontology at a time.)

Kept as its own module for the same reason as deseq2_manager.py: R
subprocess execution and file I/O that's logically distinct from the
Streamlit UI/workflow code.
"""

import itertools
import os
import shutil
import subprocess

import pandas as pd


# ---------------------------------------------------------------------------
# Organism resolution: species_key -> KEGG organism code / Reactome name
# ---------------------------------------------------------------------------

KEGG_ORGANISM_CODES = {
    "human": "hsa",
    "mouse": "mmu",
    "rat": "rno",
    "zebrafish": "dre",
    "fly": "dme",
    "drosophila": "dme",
    "worm": "cel",
    "celegans": "cel",
    "yeast": "sce",
    "ecoli": "eco",
}

REACTOME_ORGANISM_NAMES = {
    "human": "human",
    "mouse": "mouse",
    "rat": "rat",
    "zebrafish": "zebrafish",
    "fly": "fly",
    "drosophila": "fly",
    "worm": "celegans",
    "celegans": "celegans",
    "yeast": "yeast",
}

DEFAULT_MIN_GS_SIZE = 10
DEFAULT_MAX_GS_SIZE = 500
DEFAULT_SIMPLIFY_CUTOFF = 0.7

# The three directions a single database's ORA result can come in.
# "combined" is the legacy/non-recommended mode (one mixed gene list);
# "up"/"down" are the default, statistically-preferred split mode.
ORA_DIRECTIONS = ["up", "down", "combined"]


def get_kegg_organism_code(species_key):
    return KEGG_ORGANISM_CODES.get(species_key)


def get_reactome_organism_name(species_key):
    return REACTOME_ORGANISM_NAMES.get(species_key)


def clusterprofiler_tools_available():
    return shutil.which("Rscript") is not None


DATABASE_OPTIONS = ["GO", "KEGG", "Reactome"]
GO_ONTOLOGY_OPTIONS = {
    "Biological Process (BP)": "BP",
    "Molecular Function (MF)": "MF",
    "Cellular Component (CC)": "CC",
    "All three (BP + MF + CC)": "ALL",
}


def databases_available_for_species(species_key):
    return {
        "GO": True,
        "KEGG": get_kegg_organism_code(species_key) is not None,
        "Reactome": get_reactome_organism_name(species_key) is not None,
    }


def simplify_available_for_ontology(go_ontology):
    return go_ontology != "ALL"


# ---------------------------------------------------------------------------
# R script: per-database block builders (shared by ORA's combined/split modes)
# ---------------------------------------------------------------------------
#
# Rather than maintaining separate large R template strings for every
# combination of {database} x {combined/split} x {simplify on/off}, each
# database's R snippet is built dynamically in Python from small,
# composable pieces -- this keeps the actual R logic for "run enrichGO
# and save the result" defined in exactly ONE place per database,
# regardless of how many different (gene_var, out_name, guard) contexts
# it's used in.

def _build_simplify_snippet(result_var, go_ontology, simplify_go, simplify_cutoff):
    """
    Build the (possibly empty) R snippet that calls
    clusterProfiler::simplify() on result_var, IF simplify_go is True
    AND go_ontology is a single sub-ontology (not "ALL" -- see module
    docstring). Wrapped in tryCatch so a missing GOSemSim package
    degrades gracefully (analysis still completes using the
    un-simplified result) rather than hard-failing.
    """
    if not (simplify_go and simplify_available_for_ontology(go_ontology)):
        return ""
    return f'''
if (!is.null({result_var}) && nrow(as.data.frame({result_var})) > 0) {{
  tryCatch({{
    suppressMessages(library(GOSemSim))
    {result_var} <- clusterProfiler::simplify({result_var}, cutoff = {simplify_cutoff}, by = "p.adjust", select_fun = min)
    cat(paste("GO terms simplified (redundancy removed): now", nrow(as.data.frame({result_var})), "term(s)"), "\\n")
  }}, error = function(e) {{
    cat(paste("Note: could not simplify GO results (GOSemSim may not be installed) --", conditionMessage(e)), "\\n")
  }})
}}
'''


def _wrap_with_direction_guard(body, gene_var, out_name, database_label):
    """
    Wrap a per-database R block in a `if (length(gene_var) >= 1) {...}
    else {...}` guard -- used for split-direction ORA, where one
    direction (up or down) genuinely having zero significant genes at
    the chosen threshold is an EXPECTED, non-fatal outcome for that one
    database/direction combination specifically, and should not abort
    any other direction or database's analysis (unlike combined mode,
    where zero total significant genes is treated as a fatal,
    top-level stop() -- see the top-level guard in
    build_ora_r_script's script template).
    """
    return f'''
if (length({gene_var}) >= 1) {{
{body}
}} else {{
  cat("Skipping {database_label} ({out_name}): no significant genes in this direction at the chosen threshold.\\n")
}}
'''


def _build_go_ora_block(result_var, gene_var, out_name, orgdb_package, go_ontology,
                         min_gs_size, max_gs_size, simplify_go, simplify_cutoff, guarded):
    simplify_snippet = _build_simplify_snippet(result_var, go_ontology, simplify_go, simplify_cutoff)
    body = f'''
{result_var} <- enrichGO(
  gene = {gene_var}, universe = universe, OrgDb = {orgdb_package},
  keyType = "ENTREZID", ont = "{go_ontology}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size},
  readable = TRUE
)
{simplify_snippet}
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "GO") if guarded else body


def _build_kegg_ora_block(result_var, gene_var, out_name, kegg_organism,
                           min_gs_size, max_gs_size, guarded):
    body = f'''
{result_var} <- enrichKEGG(
  gene = {gene_var}, universe = universe, organism = "{kegg_organism}",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "KEGG") if guarded else body


def _build_reactome_ora_block(result_var, gene_var, out_name, reactome_organism,
                               min_gs_size, max_gs_size, guarded):
    body = f'''
{result_var} <- enrichPathway(
  gene = {gene_var}, universe = universe, organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, readable = TRUE
)
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "Reactome") if guarded else body


# ---------------------------------------------------------------------------
# R script: ORA (enrichGO / enrichKEGG / enrichPathway)
# ---------------------------------------------------------------------------

_R_ORA_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_dir <- args[2]

df <- read.csv(input_path, stringsAsFactors = FALSE)
df <- df[!is.na(df$gene_id) & !is.na(df$padj) & !is.na(df$log2FoldChange), ]

from_type <- "{from_type}"

if (from_type == "ENTREZID") {{
  df$ENTREZID <- df$gene_id
}} else {{
  suppressMessages(
    conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
  )
  df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
}}
df <- df[!duplicated(df$ENTREZID), ]

universe <- unique(df$ENTREZID)

{gene_list_block}

write_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched terms found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "term(s)"), "\n")
}}

{run_blocks}

cat("ORA analysis completed.\n")
'''

_GENE_LIST_BLOCK_COMBINED = r'''
sig_genes <- unique(df$ENTREZID[df$padj < {padj_threshold} & abs(df$log2FoldChange) >= {lfc_threshold}])
cat(paste("Significant genes (up + down combined):", length(sig_genes)), "\n")
if (length(sig_genes) < 1) {{
  stop("No significant genes at the chosen thresholds -- cannot run ORA. Try a less strict threshold.")
}}
'''

_GENE_LIST_BLOCK_SPLIT = r'''
sig_genes_up <- unique(df$ENTREZID[df$padj < {padj_threshold} & df$log2FoldChange >= {lfc_threshold}])
sig_genes_down <- unique(df$ENTREZID[df$padj < {padj_threshold} & df$log2FoldChange <= -{lfc_threshold}])
cat(paste("Up-regulated significant genes:", length(sig_genes_up)), "\n")
cat(paste("Down-regulated significant genes:", length(sig_genes_down)), "\n")
if (length(sig_genes_up) < 1 && length(sig_genes_down) < 1) {{
  stop("No significant genes (up OR down) at the chosen thresholds -- cannot run ORA. Try a less strict threshold.")
}}
'''


def build_ora_r_script(orgdb_package, from_type, padj_threshold, lfc_threshold,
                        run_go, go_ontology, run_kegg, kegg_organism,
                        run_reactome, reactome_organism,
                        min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                        simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                        split_by_direction=True):
    """
    Build the complete R script text for an ORA run across whichever
    of GO/KEGG/Reactome are requested.

    split_by_direction: if True (the DEFAULT, and the statistically
        recommended choice -- see module docstring), runs each enabled
        database's enrichment TWICE: once against only the
        up-regulated significant genes, once against only the
        down-regulated significant genes, writing two separate result
        files per database (e.g. "ora_GO_up.csv"/"ora_GO_down.csv").
        Each direction/database combination is independently guarded --
        a direction with zero significant genes for a given database is
        cleanly skipped (logged, not fatal) without affecting any other
        direction or database.

        If False, reproduces the ORIGINAL combined behavior: one
        mixed-direction significant-gene list, one result file per
        database with no suffix (e.g. "ora_GO.csv") -- offered as an
        explicit, non-default option for anyone who specifically wants
        this or needs to compare directly against a combined-list
        result.
    """
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    if split_by_direction:
        gene_list_block = _GENE_LIST_BLOCK_SPLIT.format(padj_threshold=padj_threshold, lfc_threshold=lfc_threshold)
        direction_specs = [("up", "sig_genes_up"), ("down", "sig_genes_down")]
    else:
        gene_list_block = _GENE_LIST_BLOCK_COMBINED.format(padj_threshold=padj_threshold, lfc_threshold=lfc_threshold)
        direction_specs = [(None, "sig_genes")]

    run_blocks_parts = []
    for direction, gene_var in direction_specs:
        suffix = f"_{direction}" if direction else ""
        guarded = split_by_direction

        if run_go:
            result_var = f"go_result{suffix}"
            out_name = f"ora_GO{suffix}"
            run_blocks_parts.append(_build_go_ora_block(
                result_var, gene_var, out_name, orgdb_package, go_ontology,
                min_gs_size, max_gs_size, simplify_go, simplify_cutoff, guarded,
            ))
        if run_kegg:
            result_var = f"kegg_result{suffix}"
            out_name = f"ora_KEGG{suffix}"
            run_blocks_parts.append(_build_kegg_ora_block(
                result_var, gene_var, out_name, kegg_organism, min_gs_size, max_gs_size, guarded,
            ))
        if run_reactome:
            result_var = f"reactome_result{suffix}"
            out_name = f"ora_Reactome{suffix}"
            run_blocks_parts.append(_build_reactome_ora_block(
                result_var, gene_var, out_name, reactome_organism, min_gs_size, max_gs_size, guarded,
            ))

    return _R_ORA_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, gene_list_block=gene_list_block,
        run_blocks="\n".join(run_blocks_parts),
    )


# ---------------------------------------------------------------------------
# R script: GSEA (gseGO / gseKEGG / gsePathway)
# ---------------------------------------------------------------------------
#
# GSEA does NOT need up/down direction splitting the way ORA does: its
# entire method is built around a single, full ranking of every tested
# gene (most up-regulated to most down-regulated), and a gene set's
# Normalized Enrichment Score (NES) sign already directly tells you
# whether that set skews up (positive NES) or down (negative NES) --
# there's no equivalent "mixed-direction gene list" ambiguity to
# correct for here, since NO gene list is ever constructed in the first
# place; the whole ranked list is used as-is.

_R_GSEA_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_dir <- args[2]

df <- read.csv(input_path, stringsAsFactors = FALSE)
df <- df[!is.na(df$gene_id) & !is.na(df$log2FoldChange), ]

from_type <- "{from_type}"

if (from_type == "ENTREZID") {{
  df$ENTREZID <- df$gene_id
}} else {{
  suppressMessages(
    conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
  )
  df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
}}
df <- df[!duplicated(df$ENTREZID), ]

gene_list <- df$log2FoldChange
names(gene_list) <- df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat(paste("Ranked gene list size:", length(gene_list)), "\n")

write_gsea_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched gene sets found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "gene set(s)"), "\n")

  tryCatch({{
    suppressMessages(library(enrichplot))
    running_score_rows <- list()
    for (gs_id in res_df$ID) {{
      gs_info <- enrichplot:::gsInfo(result_obj, geneSetID = gs_id)
      gs_info$gene_set_id <- gs_id
      running_score_rows[[gs_id]] <- gs_info
    }}
    if (length(running_score_rows) > 0) {{
      combined <- do.call(rbind, running_score_rows)
      write.csv(combined, file.path(output_dir, paste0(out_name, "_running_score.csv")), row.names = FALSE)
    }}
  }}, error = function(e) {{
    cat(paste("Note: could not export running-score data for", out_name, "--", conditionMessage(e)), "\n")
  }})
}}

{run_go_block}
{run_kegg_block}
{run_reactome_block}

cat("GSEA analysis completed.\n")
'''

_R_GSEA_GO_BLOCK = r'''
go_result <- gseGO(
  geneList = gene_list, OrgDb = {orgdb_package}, keyType = "ENTREZID",
  ont = "{go_ontology}", pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size},
  eps = 0
)
{simplify_block}
write_gsea_result(go_result, "gsea_GO")
'''

_R_GSEA_KEGG_BLOCK = r'''
kegg_result <- gseKEGG(
  geneList = gene_list, organism = "{kegg_organism}", keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, eps = 0
)
write_gsea_result(kegg_result, "gsea_KEGG")
'''

_R_GSEA_REACTOME_BLOCK = r'''
reactome_result <- gsePathway(
  geneList = gene_list, organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, eps = 0
)
write_gsea_result(reactome_result, "gsea_Reactome")
'''


def build_gsea_r_script(orgdb_package, from_type, run_go, go_ontology,
                         run_kegg, kegg_organism, run_reactome, reactome_organism,
                         min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                         simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF):
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    simplify_block = _build_simplify_snippet("go_result", go_ontology, simplify_go, simplify_cutoff)

    run_go_block = _R_GSEA_GO_BLOCK.format(
        orgdb_package=orgdb_package, go_ontology=go_ontology,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_block=simplify_block,
    ) if run_go else ""
    run_kegg_block = _R_GSEA_KEGG_BLOCK.format(
        kegg_organism=kegg_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_kegg else ""
    run_reactome_block = _R_GSEA_REACTOME_BLOCK.format(
        reactome_organism=reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_reactome else ""

    return _R_GSEA_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, run_go_block=run_go_block, run_kegg_block=run_kegg_block,
        run_reactome_block=run_reactome_block,
    )


# ---------------------------------------------------------------------------
# R script: compareCluster (ORA across multiple contrasts/gene lists at once)
# ---------------------------------------------------------------------------
#
# Note: compareCluster's result object is not compatible with
# clusterProfiler::simplify(), and direction-splitting is NOT offered
# here either (out of scope for this pass -- the user's request was
# specifically about ORA; compareCluster's cross-contrast comparison
# use case is different enough that mixing in a direction axis too
# would substantially complicate the comparison view without being
# explicitly requested).

_R_COMPARE_CLUSTER_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
output_dir <- args[2]

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

from_type <- "{from_type}"

gene_lists <- list()
universes <- c()

for (i in seq_len(nrow(manifest))) {{
  contrast_name <- manifest$contrast_name[i]
  input_path <- manifest$input_path[i]

  df <- read.csv(input_path, stringsAsFactors = FALSE)
  df <- df[!is.na(df$gene_id) & !is.na(df$padj) & !is.na(df$log2FoldChange), ]

  if (from_type == "ENTREZID") {{
    df$ENTREZID <- df$gene_id
  }} else {{
    suppressMessages(
      conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
    )
    df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
  }}
  df <- df[!duplicated(df$ENTREZID), ]

  sig_genes <- unique(df$ENTREZID[df$padj < {padj_threshold} & abs(df$log2FoldChange) >= {lfc_threshold}])
  gene_lists[[contrast_name]] <- sig_genes
  universes <- c(universes, df$ENTREZID)
}}

universe <- unique(universes)
cat(paste("Contrasts compared:", length(gene_lists)), "\n")
for (nm in names(gene_lists)) {{
  cat(paste(" -", nm, ":", length(gene_lists[[nm]]), "significant gene(s)"), "\n")
}}

write_cc_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched terms found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "row(s)"), "\n")
}}

{run_go_block}
{run_kegg_block}
{run_reactome_block}

cat("compareCluster analysis completed.\n")
'''

_R_CC_GO_BLOCK = r'''
go_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichGO", universe = universe,
  OrgDb = {orgdb_package}, keyType = "ENTREZID", ont = "{go_ontology}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(go_cc, "compareCluster_GO")
'''

_R_CC_KEGG_BLOCK = r'''
kegg_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichKEGG", universe = universe,
  organism = "{kegg_organism}", keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(kegg_cc, "compareCluster_KEGG")
'''

_R_CC_REACTOME_BLOCK = r'''
reactome_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichPathway", universe = universe,
  organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(reactome_cc, "compareCluster_Reactome")
'''


def build_compare_cluster_r_script(orgdb_package, from_type, padj_threshold, lfc_threshold,
                                    run_go, go_ontology, run_kegg, kegg_organism,
                                    run_reactome, reactome_organism,
                                    min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE):
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    run_go_block = _R_CC_GO_BLOCK.format(
        orgdb_package=orgdb_package, go_ontology=go_ontology,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_go else ""
    run_kegg_block = _R_CC_KEGG_BLOCK.format(
        kegg_organism=kegg_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_kegg else ""
    run_reactome_block = _R_CC_REACTOME_BLOCK.format(
        reactome_organism=reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_reactome else ""

    return _R_COMPARE_CLUSTER_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, padj_threshold=padj_threshold, lfc_threshold=lfc_threshold,
        run_go_block=run_go_block, run_kegg_block=run_kegg_block,
        run_reactome_block=run_reactome_block,
    )


# ---------------------------------------------------------------------------
# Subprocess execution
# ---------------------------------------------------------------------------

def _run_r_script(script_text, script_path, r_args, timeout=1800):
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write(script_text)

    cmd = ["Rscript", script_path] + list(r_args)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, f"Analysis timed out after {timeout // 60} minutes."


def run_ora_analysis(input_csv_path, output_dir, work_dir, orgdb_package, from_type,
                      padj_threshold, lfc_threshold, run_go, go_ontology,
                      run_kegg, kegg_organism, run_reactome, reactome_organism,
                      min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                      simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                      split_by_direction=True):
    "Run an ORA analysis (enrichGO/enrichKEGG/enrichPathway) for one contrast. Returns (success, log)."
    os.makedirs(output_dir, exist_ok=True)
    script_text = build_ora_r_script(
        orgdb_package, from_type, padj_threshold, lfc_threshold,
        run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_go=simplify_go, simplify_cutoff=simplify_cutoff,
        split_by_direction=split_by_direction,
    )
    script_path = os.path.join(work_dir, "run_ora.R")
    return _run_r_script(script_text, script_path, [input_csv_path, output_dir])


def run_gsea_analysis(input_csv_path, output_dir, work_dir, orgdb_package, from_type,
                       run_go, go_ontology, run_kegg, kegg_organism,
                       run_reactome, reactome_organism,
                       min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                       simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF):
    os.makedirs(output_dir, exist_ok=True)
    script_text = build_gsea_r_script(
        orgdb_package, from_type, run_go, go_ontology, run_kegg, kegg_organism,
        run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_go=simplify_go, simplify_cutoff=simplify_cutoff,
    )
    script_path = os.path.join(work_dir, "run_gsea.R")
    return _run_r_script(script_text, script_path, [input_csv_path, output_dir])


def run_compare_cluster_analysis(contrast_input_paths, output_dir, work_dir, orgdb_package,
                                  from_type, padj_threshold, lfc_threshold, run_go, go_ontology,
                                  run_kegg, kegg_organism, run_reactome, reactome_organism,
                                  min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    manifest_path = os.path.join(work_dir, "compare_cluster_manifest.csv")
    manifest_df = pd.DataFrame({
        "contrast_name": list(contrast_input_paths.keys()),
        "input_path": list(contrast_input_paths.values()),
    })
    manifest_df.to_csv(manifest_path, index=False)

    script_text = build_compare_cluster_r_script(
        orgdb_package, from_type, padj_threshold, lfc_threshold,
        run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    )
    script_path = os.path.join(work_dir, "run_compare_cluster.R")
    return _run_r_script(script_text, script_path, [manifest_path, output_dir])


# ---------------------------------------------------------------------------
# Reading results back into Python
# ---------------------------------------------------------------------------

def _result_file_path(output_dir, analysis_type, database, direction=None):
    """
    direction: None (legacy/combined-mode filename, no suffix -- e.g.
        "ora_GO.csv") or "up"/"down" (split-mode filenames -- e.g.
        "ora_GO_up.csv"/"ora_GO_down.csv"). Only meaningful for
        analysis_type == "ora"; GSEA and compareCluster never use a
        direction suffix.
    """
    suffix = f"_{direction}" if direction else ""
    return os.path.join(output_dir, f"{analysis_type}_{database}{suffix}.csv")


def enrichment_result_exists(output_dir, analysis_type, database, direction=None):
    return os.path.exists(_result_file_path(output_dir, analysis_type, database, direction=direction))


def available_ora_directions(output_dir, database):
    """
    For a given database, return the list of ORA direction(s) that
    actually have a saved result file in output_dir -- some subset of
    ["up", "down", "combined"], in that fixed order. Used by the UI to
    detect whether a given ORA run used split-by-direction mode (up/
    down present) or the legacy combined mode (only "combined" present)
    -- or, in principle, both if a project's history includes runs from
    both modes at different times (each mode's files coexist safely on
    disk since they use different filenames).
    """
    available = []
    for direction in ("up", "down"):
        if enrichment_result_exists(output_dir, "ora", database, direction=direction):
            available.append(direction)
    if enrichment_result_exists(output_dir, "ora", database, direction=None):
        available.append("combined")
    return available


def read_enrichment_result(output_dir, analysis_type, database, direction=None):
    path = _result_file_path(output_dir, analysis_type, database, direction=direction)
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)

    if "GeneRatio" in df.columns and not pd.api.types.is_numeric_dtype(df["GeneRatio"]):
        def _parse_ratio(s):
            try:
                num, denom = str(s).split("/")
                return float(num) / float(denom) if float(denom) != 0 else 0.0
            except (ValueError, ZeroDivisionError):
                return None
        df["GeneRatio_numeric"] = df["GeneRatio"].apply(_parse_ratio)

    return df


def read_gsea_running_score(output_dir, database, gene_set_id):
    path = os.path.join(output_dir, f"gsea_{database}_running_score.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    subset = df[df["gene_set_id"].astype(str) == str(gene_set_id)]
    return subset if not subset.empty else None


def list_available_databases(output_dir, analysis_type):
    "For GSEA/compareCluster (no direction concept). For ORA, prefer available_ora_directions per-database instead."
    return [db for db in DATABASE_OPTIONS if enrichment_result_exists(output_dir, analysis_type, db)]


def list_available_ora_databases(output_dir):
    "Return databases with AT LEAST ONE ora result (any direction: up, down, or combined)."
    return [db for db in DATABASE_OPTIONS if available_ora_directions(output_dir, db)]


# ---------------------------------------------------------------------------
# Combining results across GO/KEGG/Reactome for a single combined plot
# ---------------------------------------------------------------------------

def build_combined_results(output_dir, analysis_type, databases, n_top_per_db=10, direction=None):
    """
    direction: for ORA results specifically, which direction's file to
        read for each database ("up", "down", "combined", or None --
        None falls back to whatever enrichment_result_exists/
        read_enrichment_result's own default resolves to, i.e. the
        legacy no-suffix filename, which is correct for GSEA/
        compareCluster or legacy-combined-mode ORA).
    """
    read_direction = None if direction == "combined" else direction
    frames = []
    for db in databases:
        df = read_enrichment_result(output_dir, analysis_type, db, direction=read_direction)
        if df is None or df.empty:
            continue
        df = df.sort_values("p.adjust").head(n_top_per_db).copy()
        df["Database"] = db
        frames.append(df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def derive_category_column(df):
    """
    Add a "Category" column to a combined-results DataFrame (see
    build_combined_results), used to drive the combined view's
    background-shading-by-category feature (see
    ontology_workspace.py's _plot_combined_multi_database).

    Rather than lumping every GO row together under the generic label
    "GO", this uses clusterProfiler's own "ONTOLOGY" column (present
    on any GO result -- populated with "BP"/"MF"/"CC" per row,
    regardless of whether the original query used a single sub-
    ontology or ont="ALL") when available, so the combined view can
    show BP/MF/CC as genuinely distinct categories, exactly as the user
    would see them if they'd run each sub-ontology separately -- this
    is strictly MORE informative than a flat "GO" bucket, whether or
    not "All three" was actually selected. Non-GO rows (KEGG, Reactome)
    simply keep their existing "Database" value as their category,
    since KEGG/Reactome pathways have no equivalent sub-division.

    Returns a NEW DataFrame (does not modify df in place).
    """
    result = df.copy()
    if "ONTOLOGY" in result.columns:
        result["Category"] = result["ONTOLOGY"].fillna(result["Database"])
    else:
        result["Category"] = result["Database"]
    return result


# ---------------------------------------------------------------------------
# Term-similarity network (for the enrichment map / "emapplot"-style plot)
# ---------------------------------------------------------------------------

def _parse_gene_id_list(gene_id_str, sep="/"):
    if pd.isna(gene_id_str) or not str(gene_id_str).strip():
        return set()
    return set(str(gene_id_str).split(sep))


def build_term_similarity_network(result_df, n_top=30, similarity_threshold=0.2,
                                   gene_id_column="geneID"):
    if result_df is None or result_df.empty:
        return pd.DataFrame(), pd.DataFrame(columns=["source", "target", "weight"])

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)

    nodes_cols = ["ID", "Description", "p.adjust"]
    if count_col:
        nodes_cols.append(count_col)
    nodes_df = df[[c for c in nodes_cols if c in df.columns]].copy()
    if count_col and count_col != "Count":
        nodes_df = nodes_df.rename(columns={count_col: "Count"})

    gene_sets = {row["ID"]: _parse_gene_id_list(row.get(gene_id_column, "")) for _, row in df.iterrows()}

    edges = []
    for id_a, id_b in itertools.combinations(gene_sets.keys(), 2):
        set_a, set_b = gene_sets[id_a], gene_sets[id_b]
        if not set_a or not set_b:
            continue
        intersection = len(set_a & set_b)
        union = len(set_a | set_b)
        jaccard = intersection / union if union > 0 else 0.0
        if jaccard >= similarity_threshold:
            edges.append({"source": id_a, "target": id_b, "weight": jaccard})

    edges_df = pd.DataFrame(edges, columns=["source", "target", "weight"])
    return nodes_df, edges_df


# ---------------------------------------------------------------------------
# Gene-concept network (for the "cnetplot"-style plot)
# ---------------------------------------------------------------------------

def build_gene_concept_network(result_df, gene_fc_map=None, n_top=10, gene_id_column="geneID"):
    if result_df is None or result_df.empty:
        empty_edges = pd.DataFrame(columns=["source", "target"])
        return pd.DataFrame(), pd.DataFrame(columns=["gene_id", "log2FoldChange"]), empty_edges

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)

    term_cols = ["ID", "Description", "p.adjust"]
    if count_col:
        term_cols.append(count_col)
    term_nodes_df = df[[c for c in term_cols if c in df.columns]].copy()
    if count_col and count_col != "Count":
        term_nodes_df = term_nodes_df.rename(columns={count_col: "Count"})

    gene_fc_map = gene_fc_map or {}
    all_genes = set()
    edges = []
    for _, row in df.iterrows():
        genes = _parse_gene_id_list(row.get(gene_id_column, ""))
        all_genes.update(genes)
        for gene in genes:
            edges.append({"source": row["ID"], "target": gene})

    gene_nodes_df = pd.DataFrame({
        "gene_id": sorted(all_genes),
        "log2FoldChange": [gene_fc_map.get(g) for g in sorted(all_genes)],
    })
    edges_df = pd.DataFrame(edges, columns=["source", "target"])
    return term_nodes_df, gene_nodes_df, edges_df


# ---------------------------------------------------------------------------
# UpSet plot data (term/gene-set intersection sizes)
# ---------------------------------------------------------------------------

def build_upset_data(result_df, n_top=10, gene_id_column="geneID"):
    if result_df is None or result_df.empty:
        return pd.DataFrame(columns=["combination", "size", "term_ids"])

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    term_ids = df["ID"].tolist()
    term_labels = dict(zip(df["ID"], df.get("Description", df["ID"])))
    gene_sets = {row["ID"]: _parse_gene_id_list(row.get(gene_id_column, "")) for _, row in df.iterrows()}

    all_genes = set().union(*gene_sets.values()) if gene_sets else set()
    membership_counts = {}
    for gene in all_genes:
        membership = tuple(sorted(tid for tid in term_ids if gene in gene_sets[tid]))
        if membership:
            membership_counts[membership] = membership_counts.get(membership, 0) + 1

    rows = []
    for membership, size in membership_counts.items():
        labels = [term_labels.get(tid, tid) for tid in membership]
        rows.append({
            "combination": " + ".join(labels),
            "size": size,
            "term_ids": membership,
        })

    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["combination", "size", "term_ids"])
    return result.sort_values("size", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Ridge plot data (GSEA rank-metric distribution per gene set)
# ---------------------------------------------------------------------------

def build_ridge_plot_data(output_dir, database, result_df, n_top=10):
    if result_df is None or result_df.empty:
        return []

    top_df = result_df.sort_values("p.adjust").head(n_top)
    rows = []
    for _, row in top_df.iterrows():
        running = read_gsea_running_score(output_dir, database, row["ID"])
        if running is None:
            continue
        member_values = running.loc[running["position"] == 1, "geneList"].tolist()
        if len(member_values) < 2:
            continue
        rows.append({
            "id": row["ID"],
            "description": row.get("Description", row["ID"]),
            "NES": row.get("NES"),
            "p.adjust": row["p.adjust"],
            "values": member_values,
        })

    rows.sort(key=lambda r: (r["NES"] if r["NES"] is not None else 0))
    return rows


def summarize_enrichment_result(result_df, padj_threshold=0.05, qvalue_threshold=None):
    if result_df is None or result_df.empty:
        return "No terms were returned for this database."

    sig_mask = result_df["p.adjust"] < padj_threshold
    if qvalue_threshold is not None and "qvalue" in result_df.columns:
        sig_mask = sig_mask & (result_df["qvalue"] < qvalue_threshold)

    n_sig = sig_mask.sum()
    if n_sig == 0:
        qval_note = f" and qvalue < {qvalue_threshold}" if qvalue_threshold is not None else ""
        return (
            f"No terms reached significance at padj < {padj_threshold}{qval_note} "
            f"({len(result_df)} term(s) were tested in total). This can "
            "happen with a small significant-gene list, or if your genes "
            "genuinely don't share a strong common biological theme at "
            "this threshold -- consider trying GSEA instead, which uses "
            "every tested gene rather than requiring a hard significance "
            "cutoff first."
        )

    sig_df = result_df[sig_mask]
    top_row = sig_df.sort_values("p.adjust").iloc[0]
    return (
        f"**{n_sig:,}** term(s)/pathway(s) significant at padj < {padj_threshold}"
        + (f" and qvalue < {qvalue_threshold}" if qvalue_threshold is not None and "qvalue" in result_df.columns else "")
        + f" (out of {len(result_df):,} tested). Most significant: "
        f"**{top_row.get('Description', top_row.get('ID', 'unknown'))}** "
        f"(padj = {top_row['p.adjust']:.2e})."
    )
