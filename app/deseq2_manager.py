"""
deseq2_manager.py

Handles differential expression analysis using DESeq2 (via an Rscript
subprocess), for the Differential Expression workspace
(differential_expression_workspace.py).

Kept as its own module for the same reason as reference_manager.py and
quantification_manager.py: this involves subprocess execution, R script
generation, and statistical output parsing that is logically distinct
from the Streamlit UI/workflow code.

--- Design notes ---

Low-count filtering: done in R as part of the DESeq2 script (not in
Python), so the exact same filtering logic that determines what "kept"
vs. "removed" means is applied consistently at analysis time. A
Python-side *preview* function (preview_low_count_filter) is provided
separately so the UI can show the user how many genes a given threshold
would keep/remove *before* committing to a full DESeq2 run, without
duplicating the actual filtering logic used at run time.

Batch effects: rather than "correcting" the counts data directly (which
DESeq2's own documentation advises against for the differential
expression test itself), a batch column -- if provided -- is added as a
covariate in the design formula (e.g. ~ batch + condition). This is the
statistically recommended approach: it lets DESeq2 account for batch
variation while estimating condition effects, without discarding
information the way a naive pre-correction would. For visualization
purposes only (so batch effects are visible to the user, not to imply
they've been "fixed"), a separate batch-adjusted PCA view is offered
using limma::removeBatchEffect on the variance-stabilized counts -- this
adjusted data is NEVER used for the actual statistical test, only for
the optional visualization.

Multivariate design / multiple contrasts: the user can include more than
one condition column in the design formula (e.g. ~ batch + genotype +
treatment), and can request multiple pairwise contrasts from the primary
condition column (e.g. treated_A vs control, treated_B vs control) in a
single DESeq2 run -- the model is fit once, then each contrast is
extracted from the same fitted model, which is both faster and
statistically more consistent than re-fitting per comparison.

Ambiguous/missing metadata values: real-world metadata frequently has
blank cells or placeholder text ("none", "N/A", "-") in condition/batch
columns -- often meaning "no treatment" (a legitimate control category)
rather than truly unknown data. detect_ambiguous_values() and
apply_missing_value_resolution() let the UI surface these explicitly and
have the user decide, per value, whether it represents a real category
or should cause those samples to be excluded -- rather than silently
passing a real NaN into DESeq2 (which would error or misbehave) or
guessing incorrectly on the user's behalf.

--- QC/diagnostics additions ---

Beyond the core DESeq2 fit + results, the R script also exports several
model-fit and sample-level diagnostic files that the Differential
Expression workspace surfaces as additional QC visuals (see
"read_dispersion_estimates", "read_size_factors",
"read_sample_distance_matrix", and the "summarize_*"/"classify_*"/
"flag_*"/"detect_*" functions below):

  - Dispersion estimates (gene-wise, fitted trend, and final/shrunk) --
    DESeq2's own recommended diagnostic (plotDispEsts in R) for whether
    the negative-binomial model fit reasonably to this dataset.
  - Size factors -- DESeq2's per-sample normalization factors; a
    sample whose size factor is a strong outlier relative to the others
    often indicates a failed/degraded library.
  - A full sample-to-sample distance matrix (computed on the same VST
    data already used for PCA) -- the standard companion to PCA for
    sample-level QC, since PCA alone only shows the first few
    components and can miss an outlier/mislabeled sample that a full
    distance-based comparison would catch.
  - Per-contrast Cook's-distance outlier flagging and independent-
    filtering status -- DESeq2 silently sets "pvalue" to NA for genes
    whose result looks driven by a single outlier sample (Cook's
    distance cutoff), and separately sets "padj" (but not "pvalue") to
    NA for low-mean genes excluded by independent filtering to improve
    power. These are two distinct, legitimate DESeq2 behaviors that are
    otherwise easy to mistake for "missing data" -- the R script now
    tags each with an explicit boolean column (and, for Cook's
    outliers, which specific sample had the largest Cook's distance for
    that gene) so the workspace can explain -- rather than just
    silently show -- every NA row in the results table.

None of this changes DESeq2's actual statistical test in any way -- it
only extracts and exposes data that DESeq2 already computes internally
as part of a normal run, the same way the existing PCA/VST export does.
"""
import csv
import os
import re
import subprocess

import pandas as pd

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

def deseq2_tools_available():
    """
    Check whether Rscript is available. Whether the DESeq2/limma R
    packages themselves are installed can only be confirmed by actually
    attempting to load them (checked at the start of the R script itself,
    which reports a clear error back to Python if a package is missing --
    same pattern as the tximport dependency issues encountered earlier).
    """
    import shutil
    return shutil.which("Rscript") is not None


# ---------------------------------------------------------------------------
# Ambiguous / missing metadata value detection and resolution
# ---------------------------------------------------------------------------

# Common string values that represent "missing" but were entered as text
# rather than a true blank/NaN cell in a spreadsheet (e.g. someone typed
# "none", "N/A", or "-" instead of leaving the cell truly empty).
# Comparison is case-insensitive and whitespace-trimmed.
MISSING_LIKE_STRINGS = {"none", "na", "n/a", "null", "nan", "", "-", "unknown", "missing"}


def detect_ambiguous_values(series):
    """
    Detect values in a metadata column that are either:
      - true missing (NaN / None) -- an actually blank cell
      - "missing-like" text (case-insensitive match against common
        placeholder strings like "none", "N/A", "-", etc.) -- a cell
        that has *something* in it, but that something likely means
        "not applicable" or "no treatment" rather than being a real
        category name on its own merit.

    This distinction matters because a truly blank cell almost always
    means "we don't know / this wasn't recorded" (a candidate for
    exclusion), whereas a text value like "none" often means something
    meaningful and should become an explicit category (e.g. a control
    group with "no_antibiotic") rather than being dropped.

    Returns a dict:
        {
            "true_na_count": int,
            "missing_like_values": {original_value: count, ...},
        }
    """
    true_na_mask = series.isna()
    true_na_count = int(true_na_mask.sum())
    non_null = series.dropna().astype(str)
    normalized = non_null.str.strip().str.lower()
    missing_like_mask = normalized.isin(MISSING_LIKE_STRINGS)
    missing_like_values = {}
    for original_val in non_null[missing_like_mask].unique():
        count = int((non_null[missing_like_mask] == original_val).sum())
        missing_like_values[original_val] = count
    return {
        "true_na_count": true_na_count,
        "missing_like_values": missing_like_values,
    }


def validate_design_columns(meta_df, columns):
    """
    Check each candidate design/batch column for issues that would
    break or invalidate a DESeq2 model, so these can be caught
    immediately in the UI rather than after a multi-minute R run fails:
      - constant (only 1 unique value across all samples): DESeq2 will
        hard-error on this ("design contains one or more variables with
        all samples having the same value"), since a variable with zero
        variation cannot be used to estimate any effect. This commonly
        happens with metadata auto-filled from public repositories
        (e.g. SRA), where fields like "strain" or "organism" are the
        same for every sample in a single-species study and were never
        meant to be a design factor.
      - identifier-like (as many unique values as there are samples):
        DESeq2 would technically accept this without erroring, but a
        column like this (e.g. a free-text "sample_title" description)
        is essentially never a meaningful experimental factor -- it
        typically means every sample has a distinct value, which makes
        the model unidentifiable/meaningless for that term.

    Returns a dict {column: {"n_unique": int, "n_samples": int,
                              "is_constant": bool,
                              "is_identifier_like": bool}}.
    """
    n_samples = len(meta_df)
    results = {}
    for col in columns:
        if col not in meta_df.columns:
            continue
        n_unique = meta_df[col].astype(str).nunique(dropna=True)
        results[col] = {
            "n_unique": n_unique,
            "n_samples": n_samples,
            "is_constant": n_unique <= 1,
            "is_identifier_like": n_unique >= n_samples and n_samples > 1,
        }
    return results


def check_replication(meta_df, design_columns, batch_column=None, interaction_terms=None):
    """
    Check whether there is enough sample replication for DESeq2 to fit
    the requested design and estimate dispersion (within-group
    variability).

    This mirrors what DESeq2's design matrix actually requires, rather
    than naively demanding replication for every possible combination
    of every selected column:

      - Main effects (any design column, or the batch column, that
        isn't part of an interaction term): only need >= 2 samples at
        *each individual level* of that column on its own (e.g.
        "ampicillin" needs >= 2 samples) -- NOT crossed with every
        other selected column. Crossing everything unconditionally is
        a common false positive: a standard randomized-block design
        (e.g. 3 batches x 2 conditions, 1 sample per cell = 6 samples
        total) is a perfectly valid, identifiable additive model
        (~ batch + condition) even though every individual
        batch:condition cell has only 1 sample.
      - Interaction terms (e.g. "genotype:treatment", passed via
        interaction_terms): DO require >= 2 samples per *combination*
        of that specific pair of columns, since estimating an
        interaction coefficient needs replication within each cell of
        that particular cross. Only the columns actually named in an
        interaction term are cross-checked this way; every other
        design/batch column is still checked individually above.

    Without at least 2 replicates where they're actually needed,
    DESeq2 has no data left to distinguish real biological differences
    from noise, and will hard-error with "the design matrix has the
    same number of samples and coefficients to fit ...
    checkForExperimentalReplicates". This is a very common issue when
    metadata columns are auto-filled from a public repository (e.g.
    SRA) and happen to have a unique combination of factor levels for
    every single downloaded sample -- catching it here avoids a
    multi-minute R run that errors out at the dispersion estimation
    step.

    interaction_terms: optional list of strings in "colA:colB" form
        (see build_interaction_term_options), using the same column
        names already present in design_columns.

    Returns a dict:
        {
            "is_valid": bool,
            "group_counts": {group_label: count, ...},
            "under_replicated_groups": [group_label, ...],
        }
    """
    interaction_terms = interaction_terms or []

    all_main_columns = list(design_columns)
    if batch_column:
        all_main_columns.append(batch_column)

    if not all_main_columns:
        return {"is_valid": True, "group_counts": {}, "under_replicated_groups": []}

    group_counts = {}
    under_replicated = []

    # --- Main effects: checked individually, per column, per level ---
    for col in all_main_columns:
        if col not in meta_df.columns:
            continue
        level_counts = meta_df.groupby(col, dropna=False).size()
        for level, count in level_counts.items():
            label = f"{col}={level}"
            group_counts[label] = int(count)
            if count < 2:
                under_replicated.append(label)

    # --- Interaction terms: checked crossed, per pair of columns ---
    for term in interaction_terms:
        cols = [c for c in term.split(":") if c]
        if len(cols) < 2 or not all(c in meta_df.columns for c in cols):
            continue
        grouped = meta_df.groupby(cols, dropna=False).size()
        for group_key, count in grouped.items():
            if isinstance(group_key, tuple):
                inner_label = " / ".join(f"{col}={val}" for col, val in zip(cols, group_key))
            else:
                inner_label = f"{cols[0]}={group_key}"
            label = f"Interaction {term}: {inner_label}"
            group_counts[label] = int(count)
            if count < 2:
                under_replicated.append(label)

    return {
        "is_valid": len(under_replicated) == 0,
        "group_counts": group_counts,
        "under_replicated_groups": under_replicated,
    }


def has_any_ambiguous_values(meta_df, columns):
    """
    Check whether any of the given metadata columns contain ambiguous
    (true-NaN or missing-like text) values. Used to decide whether to
    show the resolution UI at all for a given design/batch column
    selection -- most datasets won't need this step, so it should only
    appear when actually relevant.
    """
    for col in columns:
        result = detect_ambiguous_values(meta_df[col])
        if result["true_na_count"] > 0 or result["missing_like_values"]:
            return True
    return False


def apply_missing_value_resolution(meta_df, column, resolution_map, exclude_samples):
    """
    Apply the user's decisions about ambiguous values in a single
    metadata column.

    resolution_map: dict mapping each ambiguous value to a replacement
        category label the user chose (e.g. {"none": "no_antibiotic"}).
        Use the special key "__NaN__" to specify the replacement for
        true blank/NaN cells specifically, since NaN can't be used as a
        normal dict key for string comparison.
    exclude_samples: list of sample names (matching the "sample" column)
        whose rows should be dropped entirely -- used when the user
        decides an ambiguous value truly represents missing/unusable
        data rather than a meaningful category.

    Returns a new DataFrame with resolutions applied. Any ambiguous
    value NOT mentioned in resolution_map and not covered by
    exclude_samples is left as-is (e.g. if the user only resolves some
    of several ambiguous values found).
    """
    result_df = meta_df.copy()
    if "__NaN__" in resolution_map:
        result_df[column] = result_df[column].fillna(resolution_map["__NaN__"])
    for ambiguous_val, replacement in resolution_map.items():
        if ambiguous_val == "__NaN__":
            continue
        mask = result_df[column].astype(str).str.strip().str.lower() == ambiguous_val.strip().lower()
        result_df.loc[mask, column] = replacement
    if exclude_samples:
        result_df = result_df[~result_df["sample"].isin(exclude_samples)]
    return result_df


# ---------------------------------------------------------------------------
# Python-side low-count filter preview (does not affect the actual R run;
# purely so the UI can show "N genes would be kept" before committing)
# ---------------------------------------------------------------------------

def preview_low_count_filter(counts_df, min_count, min_samples):
    """
    Given a counts matrix DataFrame (first column "gene_id", remaining
    columns = samples), compute how many genes would pass a filter of
    "at least min_count reads in at least min_samples samples" --
    matching the exact filter logic applied inside the DESeq2 R script
    (see _DESEQ2_R_SCRIPT's filtering step).

    Returns a dict: {"total_genes": int, "genes_kept": int,
                      "genes_removed": int, "pct_kept": float}
    """
    sample_cols = [c for c in counts_df.columns if c != "gene_id"]
    total_genes = len(counts_df)
    if total_genes == 0 or not sample_cols:
        return {"total_genes": 0, "genes_kept": 0, "genes_removed": 0, "pct_kept": 0.0}
    passes_threshold = (counts_df[sample_cols] >= min_count).sum(axis=1) >= min_samples
    genes_kept = int(passes_threshold.sum())
    genes_removed = total_genes - genes_kept
    pct_kept = round(100 * genes_kept / total_genes, 1) if total_genes else 0.0
    return {
        "total_genes": total_genes,
        "genes_kept": genes_kept,
        "genes_removed": genes_removed,
        "pct_kept": pct_kept,
    }


# ---------------------------------------------------------------------------
# R script for DESeq2 analysis
# ---------------------------------------------------------------------------
#
# This script is written to a temp file and executed via Rscript. It
# takes a single positional argument: a path to a JSON "job spec" file
# (written by run_deseq2_analysis below), since the number of
# parameters here (design columns, multiple contrasts, optional batch
# column, filtering thresholds) is too complex to pass cleanly as plain
# positional CLI arguments.
_DESEQ2_R_SCRIPT = r'''
suppressMessages({
  library(jsonlite)
  library(DESeq2)
})
args <- commandArgs(trailingOnly = TRUE)
job_spec_path <- args[1]
job <- fromJSON(job_spec_path)
counts_df <- read.csv(job$counts_path, row.names = 1, check.names = FALSE)
meta_df <- read.csv(job$metadata_path, row.names = 1, check.names = FALSE)
# Align sample order between counts columns and metadata rows explicitly
# -- DESeqDataSetFromMatrix requires this and does not reorder for you.
common_samples <- intersect(colnames(counts_df), rownames(meta_df))
if (length(common_samples) < 2) {
  stop("Fewer than 2 samples in common between counts matrix and metadata after matching sample names.")
}
counts_df <- counts_df[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]
# Ensure every design/batch column is treated as a categorical factor,
# not a numeric or character vector, since DESeq2's design formula
# requires factor columns.
for (col in job$design_columns) {
  meta_df[[col]] <- as.factor(meta_df[[col]])
}
if (!is.null(job$batch_column) && job$batch_column != "") {
  meta_df[[job$batch_column]] <- as.factor(meta_df[[job$batch_column]])
}
# Low-count filtering: keep genes with at least min_count reads in at
# least min_samples samples. This mirrors preview_low_count_filter() on
# the Python side exactly, so the "N genes kept" preview shown to the
# user before running matches what actually happens here.
keep <- rowSums(counts_df >= job$min_count) >= job$min_samples
counts_df <- counts_df[keep, , drop = FALSE]
cat(paste("Genes after low-count filtering:", nrow(counts_df), "of", length(keep)), "\n")
if (nrow(counts_df) == 0) {
  stop("No genes remain after low-count filtering -- try a lower threshold.")
}
# Build the full design formula. If a batch column is provided, it's
# included as a covariate (e.g. ~ batch + condition) rather than used to
# "pre-correct" the counts -- this lets DESeq2 properly account for batch
# variation while estimating the condition effect(s) of interest.
# job$interaction_terms is a list of strings already in R interaction
# syntax (e.g. "genotype:treatment"), built on the Python side.
formula_terms <- job$design_columns
if (!is.null(job$batch_column) && job$batch_column != "") {
  formula_terms <- c(job$batch_column, formula_terms)
}
if (!is.null(job$interaction_terms) && length(job$interaction_terms) > 0) {
  formula_terms <- c(formula_terms, job$interaction_terms)
}
design_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))
cat(paste("Full design formula:", deparse(design_formula)), "\n")
# --- Fit the model ---
# Wald test (job$test_type == "wald"): the standard fit, used for every
# pairwise contrast in job$contrasts (results(dds, contrast = ...)).
#
# LRT (job$test_type == "lrt"): DESeq2's equivalent of an ANOVA omnibus
# test. Rather than asking "how does level A differ from the
# reference?" (what Wald/contrasts test), LRT asks "does this factor --
# across ALL its levels -- or this interaction term, matter at all?" by
# comparing the full model against a "reduced" model missing the
# term(s) being tested. Each entry in job$lrt_tests specifies its own
# reduced formula, so a separate DESeq() fit is performed per LRT test
# (dispersion estimation is re-run each time, which is the standard,
# statistically safe way to support multiple different reduced models
# rather than trying to reuse one fit across differently-specified
# tests).
dds <- DESeqDataSetFromMatrix(countData = round(counts_df), colData = meta_df, design = design_formula)
dds <- DESeq(dds)

# --- Model-fit diagnostics: dispersion estimates + size factors ---
# Exported once per project (from the primary Wald-fitted dds above, or
# -- for an LRT-only run -- from the first LRT test's own fit below,
# since dispersion estimation happens as part of every DESeq() call
# regardless of test type). These let the UI show DESeq2's own
# recommended dispersion-fit diagnostic (plotDispEsts in R) and a size-
# factor QC table without requiring any extra R computation beyond what
# DESeq2 already does internally as part of a normal run.
disp_df <- data.frame(
  gene_id = rownames(dds),
  baseMean = mcols(dds)$baseMean,
  dispGeneEst = mcols(dds)$dispGeneEst,
  dispFit = mcols(dds)$dispFit,
  dispersion = mcols(dds)$dispersion
)
write.csv(disp_df, file.path(job$output_dir, "dispersion_estimates.csv"), row.names = FALSE)

size_factor_df <- data.frame(
  sample = names(sizeFactors(dds)),
  size_factor = as.numeric(sizeFactors(dds))
)
write.csv(size_factor_df, file.path(job$output_dir, "size_factors.csv"), row.names = FALSE)

# --- Variance-stabilized data + PCA, for visualization ---
vsd <- vst(dds, blind = TRUE)
pca_input <- t(assay(vsd))
pca_result <- prcomp(pca_input, scale. = FALSE)
pct_var <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
pca_df <- as.data.frame(pca_result$x[, 1:min(4, ncol(pca_result$x))])
pca_df$sample <- rownames(pca_df)
# Carry along every metadata column so the UI can color the PCA plot by
# any of them (condition, batch, etc.) without a separate R round-trip.
for (col in colnames(meta_df)) {
  pca_df[[col]] <- meta_df[common_samples, col]
}
write.csv(pca_df, file.path(job$output_dir, "pca_coordinates.csv"), row.names = FALSE)
writeLines(paste(pct_var, collapse = ","), file.path(job$output_dir, "pca_percent_variance.txt"))

# --- Sample-to-sample distance matrix (sample-level QC, companion to PCA) ---
# Computed on the SAME (non-batch-adjusted) VST data used for the "raw"
# PCA view above -- this is the standard complementary QC view
# recommended alongside PCA, since PCA's first 2-4 components can
# sometimes miss an outlier/mislabeled sample that a full pairwise
# distance comparison catches. Exported as a plain sample x sample
# matrix (Euclidean distance on VST values); the UI performs
# hierarchical clustering client-side for the heatmap's row/column
# ordering, so this file only needs to carry the raw distances.
sample_dist_matrix <- as.matrix(dist(t(assay(vsd))))
sample_dist_df <- as.data.frame(sample_dist_matrix)
sample_dist_df$sample <- rownames(sample_dist_df)
sample_dist_df <- sample_dist_df[, c("sample", colnames(sample_dist_matrix))]
write.csv(sample_dist_df, file.path(job$output_dir, "sample_distance_matrix.csv"), row.names = FALSE)

# --- Optional batch-adjusted PCA view (visualization only) ---
# This uses limma::removeBatchEffect purely to produce a second PCA view
# where batch-driven variation is visually suppressed, so the user can
# compare "before" (raw VST) vs. "after" (batch-adjusted) and judge how
# much of the variance in the first PCA was attributable to batch. This
# adjusted matrix is NOT used for the actual DESeq2 statistical test --
# the real test already accounts for batch correctly via the design
# formula above, which is the statistically appropriate approach.
if (!is.null(job$batch_column) && job$batch_column != "") {
  suppressMessages(library(limma))
  batch_adjusted <- removeBatchEffect(assay(vsd), batch = meta_df[[job$batch_column]])
  pca_adj_result <- prcomp(t(batch_adjusted), scale. = FALSE)
  pct_var_adj <- round(100 * (pca_adj_result$sdev^2 / sum(pca_adj_result$sdev^2)), 1)
  pca_adj_df <- as.data.frame(pca_adj_result$x[, 1:min(4, ncol(pca_adj_result$x))])
  pca_adj_df$sample <- rownames(pca_adj_df)
  for (col in colnames(meta_df)) {
    pca_adj_df[[col]] <- meta_df[common_samples, col]
  }
  write.csv(pca_adj_df, file.path(job$output_dir, "pca_coordinates_batch_adjusted.csv"), row.names = FALSE)
  writeLines(paste(pct_var_adj, collapse = ","), file.path(job$output_dir, "pca_percent_variance_batch_adjusted.txt"))
}
# --- Normalized counts (for user reference / downstream use) ---
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$gene_id <- rownames(norm_counts)
norm_counts <- norm_counts[, c("gene_id", common_samples)]
write.csv(norm_counts, file.path(job$output_dir, "normalized_counts.csv"), row.names = FALSE)

# --- Helper: which sample had the largest Cook's distance for each gene ---
# Used below to annotate outlier-flagged genes in the results tables
# with WHICH specific sample's count value looks like it's driving that
# gene's result -- much more actionable than just knowing a gene was
# flagged. assays(dds_obj)[["cooks"]] is a genes x samples matrix of
# per-gene, per-sample Cook's distances, populated by DESeq() as part
# of its normal fitting process (same data DESeq2 itself uses
# internally to decide whether to null out a gene's p-value).
compute_max_cooks_sample <- function(dds_obj, gene_ids) {
  cooks_mat <- tryCatch(assays(dds_obj)[["cooks"]], error = function(e) NULL)
  if (is.null(cooks_mat)) {
    return(rep(NA_character_, length(gene_ids)))
  }
  vapply(gene_ids, function(g) {
    if (!(g %in% rownames(cooks_mat))) return(NA_character_)
    row_vals <- cooks_mat[g, ]
    if (all(is.na(row_vals))) return(NA_character_)
    colnames(cooks_mat)[which.max(row_vals)]
  }, character(1))
}

# --- Wald test: per-contrast pairwise results ---
# job$contrasts is a list of {name, column, level1, level2} objects,
# where level1 is the numerator (e.g. "treated") and level2 is the
# denominator (e.g. "control") -- matching DESeq2's results(contrast=...)
# convention of c(column, numerator, denominator). Every entry here uses
# the single Wald-fitted dds above (no re-fitting needed per contrast).
if (job$test_type == "wald" && !is.null(job$contrasts) && length(job$contrasts$name) > 0) {
  # NOTE: job$contrasts is serialized from the Python side in
  # column-oriented form (contrasts_df.to_dict(orient="list")), i.e. a
  # JSON object like {"name": [...], "column": [...], ...} rather than
  # a JSON array of per-row objects. jsonlite's fromJSON() only
  # auto-promotes a JSON *array of row objects* into a data.frame -- a
  # single JSON object whose fields are equal-length arrays (this
  # case) is parsed as a plain named list instead, NOT a data.frame.
  # nrow() on a plain list returns NULL rather than a row count, which
  # previously caused "missing value where TRUE/FALSE needed" once
  # NULL was compared with > 0 inside this if() condition. Using
  # length(job$contrasts$name) instead (mirroring how the LRT branch
  # below already correctly uses length(job$lrt_tests$name)) works
  # correctly regardless of whether jsonlite happens to return a
  # data.frame or a plain list here.
  for (i in seq_along(job$contrasts$name)) {
    # Indexed via $column_name[i] rather than row-subsetting
    # (job$contrasts[i, ]) first, to avoid any ambiguity in how a
    # data.frame with potential list-columns behaves under single-row
    # subsetting -- explicit per-column indexing is unambiguous.
    contrast_name <- job$contrasts$name[i]
    contrast_column <- job$contrasts$column[i]
    contrast_level1 <- job$contrasts$level1[i]
    contrast_level2 <- job$contrasts$level2[i]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)

    res <- results(dds, contrast = c(contrast_column, contrast_level1, contrast_level2))
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)

    # --- Distinguish the two legitimate DESeq2 causes of an NA result ---
    # Cook's-distance outlier: DESeq2 sets "pvalue" itself to NA when a
    # gene's result looks driven by a single outlier sample (detected
    # via cooksCutoff, DESeq2's default behavior for results()).
    # Independent filtering: DESeq2 separately sets "padj" (but leaves
    # "pvalue" intact) to NA for low-mean genes it excludes from
    # multiple-testing correction to improve power -- a completely
    # different, unrelated reason for an NA. Tagging both explicitly
    # here (rather than leaving a bare NA for the user to puzzle over)
    # is what lets the workspace explain -- not just display -- every
    # NA row.
    res_df$cooks_outlier_flagged <- is.na(res_df$pvalue) & !is.na(res_df$baseMean) & res_df$baseMean > 0
    res_df$low_count_filtered <- !is.na(res_df$pvalue) & is.na(res_df$padj)
    res_df$max_cooks_sample <- compute_max_cooks_sample(dds, res_df$gene_id)

    # Record the independent filtering threshold DESeq2 chose for this
    # specific contrast (it's optimized per-results-call, so it can
    # differ slightly between contrasts even from the same fitted dds).
    filter_threshold <- metadata(res)$filterThreshold
    if (!is.null(filter_threshold) && length(filter_threshold) > 0) {
      writeLines(as.character(as.numeric(filter_threshold)), file.path(job$output_dir, paste0("filter_threshold_", safe_name, ".txt")))
    }

    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                          "cooks_outlier_flagged", "low_count_filtered", "max_cooks_sample")]
    res_df <- res_df[order(res_df$padj), ]
    write.csv(res_df, file.path(job$output_dir, paste0("results_", safe_name, ".csv")), row.names = FALSE)
    cat(paste("Contrast", contrast_name, "(Wald) - genes with padj < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE)), "\n")
  }
}
# --- LRT: omnibus tests, each with its own reduced model ---
# job$lrt_tests is a list of {name, reduced_terms} objects, where
# reduced_terms is the full design's term list minus whatever term(s)
# are being tested (e.g. testing "treatment" overall: reduced_terms
# would be the full formula's terms minus "treatment"). A separate
# DESeq() fit (with test="LRT") is performed per entry, since each may
# specify a different reduced model.
#
# Unlike Wald results, LRT's results() output does not have a natural
# "log2FoldChange" for factors with more than 2 levels (there's no
# single fold-change when testing 3+ groups at once) -- log2FoldChange
# in this case reflects an arbitrary reference comparison internal to
# the model and should be interpreted cautiously; the "stat" and "padj"
# columns are the primary output of interest for LRT (they reflect the
# omnibus test itself, analogous to an ANOVA F-test p-value).
if (job$test_type == "lrt" && !is.null(job$lrt_tests) && length(job$lrt_tests$name) > 0) {
  for (i in seq_along(job$lrt_tests$name)) {
    # reduced_terms is a nested (variable-length) array per test, which
    # jsonlite represents as a list-column -- accessed via [[i]] to
    # extract the i-th test's character vector of terms, rather than via
    # row-subsetting the whole data.frame first (which behaves
    # ambiguously when list-columns are present).
    test_name <- job$lrt_tests$name[i]
    reduced_terms <- job$lrt_tests$reduced_terms[[i]]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", test_name)
    if (length(reduced_terms) == 0) {
      reduced_formula <- ~ 1
    } else {
      reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
    }
    dds_lrt <- DESeqDataSetFromMatrix(countData = round(counts_df), colData = meta_df, design = design_formula)
    dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = reduced_formula)
    res <- results(dds_lrt)
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)

    # Same NA-cause tagging as the Wald branch above -- see the comment
    # there for why these two columns exist. Uses dds_lrt's own Cook's
    # matrix, since this is a separately-refit model (its own
    # dispersion/Cook's distances, not shared with the Wald-fitted dds).
    res_df$cooks_outlier_flagged <- is.na(res_df$pvalue) & !is.na(res_df$baseMean) & res_df$baseMean > 0
    res_df$low_count_filtered <- !is.na(res_df$pvalue) & is.na(res_df$padj)
    res_df$max_cooks_sample <- compute_max_cooks_sample(dds_lrt, res_df$gene_id)

    filter_threshold <- metadata(res)$filterThreshold
    if (!is.null(filter_threshold) && length(filter_threshold) > 0) {
      writeLines(as.character(as.numeric(filter_threshold)), file.path(job$output_dir, paste0("filter_threshold_", safe_name, ".txt")))
    }

    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                          "cooks_outlier_flagged", "low_count_filtered", "max_cooks_sample")]
    res_df <- res_df[order(res_df$padj), ]
    write.csv(res_df, file.path(job$output_dir, paste0("results_", safe_name, ".csv")), row.names = FALSE)
    cat(paste("LRT test", test_name, "(reduced model:", deparse(reduced_formula), ") - genes with padj < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE)), "\n")
  }
}
cat("DESeq2 analysis completed successfully.\n")
'''


def summarize_deseq2_config(config):
    """
    Build a short, human-readable one-line summary of a saved DESeq2
    configuration dict (the same shape passed to
    project_manager.save_deseq2_config), for display when a returning
    user reopens a project that already has results on disk -- e.g.
    "design = treatment + cell_line (batch: cell_line) | Wald: 3
    contrast(s) vs. Untreated" -- so they can tell at a glance whether
    the existing results match what they're looking for before
    deciding whether to expand Steps 1-4 and re-run with different
    settings.

    Returns "" if config is None/empty.
    """
    if not config:
        return ""

    design_columns = config.get("design_columns") or []
    batch_column = config.get("batch_column")
    interaction_terms = config.get("interaction_terms") or []
    test_type = config.get("test_type", "wald")

    formula_terms = list(design_columns) + list(interaction_terms)
    design_part = " + ".join(formula_terms) if formula_terms else "(none)"
    batch_part = f" (batch: {batch_column})" if batch_column else ""

    if test_type == "wald":
        contrasts = config.get("contrasts") or []
        n = len(contrasts)
        ref_levels = sorted({c.get("level2", "") for c in contrasts if c.get("level2")})
        ref_part = f" vs. {', '.join(ref_levels)}" if ref_levels else ""
        test_part = f"Wald: {n} contrast{'s' if n != 1 else ''}{ref_part}"
    else:
        lrt_tests = config.get("lrt_tests") or []
        n = len(lrt_tests)
        test_part = f"LRT: {n} test{'s' if n != 1 else ''}"

    return f"design = {design_part}{batch_part} | {test_part}"


def build_full_formula_terms(design_columns, batch_column, interaction_terms):
    """
    Build the ordered list of terms for the full design formula, e.g.
    ["batch", "genotype", "treatment", "genotype:treatment"].

    interaction_terms: list of strings already in "colA:colB" form (see
    build_interaction_term_options below for how these are generated
    from user-selected column pairs).
    """
    terms = []
    if batch_column:
        terms.append(batch_column)
    terms.extend(design_columns)
    terms.extend(interaction_terms or [])
    return terms


def build_interaction_term_options(design_columns):
    """
    Given the selected design columns, return every possible pairwise
    interaction term string (e.g. "genotype:treatment") the user could
    choose to add to the model. Only pairwise (2-way) interactions are
    offered -- higher-order (3-way+) interactions are rarely
    interpretable/well-powered in a typical RNA-seq sample size and are
    intentionally not exposed here to avoid encouraging an
    under-powered, hard-to-interpret model.
    """
    options = []
    for i in range(len(design_columns)):
        for j in range(i + 1, len(design_columns)):
            options.append(f"{design_columns[i]}:{design_columns[j]}")
    return options


def build_reduced_formula_terms(full_terms, terms_to_drop):
    """
    For an LRT test, the reduced model = full model minus the term(s)
    being tested. terms_to_drop must exactly match entries in
    full_terms (a design column name, or an interaction like
    "genotype:treatment").
    """
    return [t for t in full_terms if t not in terms_to_drop]


def run_deseq2_analysis(counts_matrix_path, metadata_path, design_columns,
                         batch_column, min_count, min_samples,
                         output_dir, work_dir, test_type="wald",
                         interaction_terms=None, contrasts=None, lrt_tests=None):
    """
    Run a full DESeq2 analysis: low-count filtering, model fitting
    (optionally with a batch covariate and/or interaction terms), VST +
    PCA computation, and results extraction using either the Wald test
    (pairwise contrasts) or the Likelihood Ratio Test (LRT, an ANOVA-
    style omnibus test) -- all in a single R script invocation.

    counts_matrix_path: path to the gene counts matrix CSV (columns:
        gene_id, sample1, sample2, ...)
    metadata_path: path to a metadata CSV with a "sample" column
        matching the counts matrix's sample columns, plus condition/
        batch columns
    design_columns: list of column names (from metadata) to include as
        the variable(s) of interest in the design formula, e.g.
        ["treatment"] or ["genotype", "treatment"] for a multivariate
        design
    batch_column: column name to include as a covariate for batch
        effects, or None/"" if not applicable
    min_count, min_samples: low-count filtering thresholds -- a gene is
        kept if it has at least min_count reads in at least min_samples
        samples
    output_dir: where DESeq2's output CSVs are written
    work_dir: scratch directory for the temp R script + job spec JSON
    test_type: "wald" (pairwise contrasts) or "lrt" (omnibus/ANOVA-style
        test comparing full vs. reduced models)
    interaction_terms: list of strings in "colA:colB" form (see
        build_interaction_term_options), added to the full design
        formula alongside design_columns and batch_column
    contrasts: (required if test_type == "wald") list of dicts, each
        {"name": str, "column": str, "level1": str, "level2": str} --
        level1 is the numerator (e.g. "treated"), level2 is the
        denominator (e.g. "control")
    lrt_tests: (required if test_type == "lrt") list of dicts, each
        {"name": str, "reduced_terms": [list of term strings]} -- the
        reduced model's terms for that specific omnibus test (typically
        the full model's terms minus the term being tested)

    Returns (success: bool, log: str).
    """
    import json

    if not deseq2_tools_available():
        return False, (
            "Rscript was not found on this system. R with the DESeq2 "
            "package needs to be installed in your environment (it's "
            "included in the project's Dockerfile) before this step can "
            "run."
        )

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    # Normalize metadata to have "sample" as the index column expected
    # by the R script's read.csv(..., row.names=1).
    meta_df = pd.read_csv(metadata_path)
    if "sample" not in meta_df.columns:
        return False, "Metadata file is missing a 'sample' column."

    meta_for_r_path = os.path.join(work_dir, "metadata_for_deseq2.csv")
    meta_df.to_csv(meta_for_r_path, index=False)
    # DESeq2's read.csv(row.names=1) usage in the R script expects the
    # first column to become rownames -- re-save the metadata with
    # "sample" explicitly as the first column and everything else after,
    # so R's row.names=1 picks it up correctly regardless of original
    # column ordering in the source file.
    cols = ["sample"] + [c for c in meta_df.columns if c != "sample"]
    meta_df[cols].to_csv(meta_for_r_path, index=False)

    job_spec = {
        "counts_path": os.path.abspath(counts_matrix_path),
        "metadata_path": os.path.abspath(meta_for_r_path),
        "design_columns": design_columns,
        "batch_column": batch_column or "",
        "interaction_terms": interaction_terms or [],
        "min_count": min_count,
        "min_samples": min_samples,
        "output_dir": os.path.abspath(output_dir),
        "test_type": test_type,
    }

    if test_type == "wald":
        contrasts_df = pd.DataFrame(contrasts or [])
        job_spec["contrasts"] = contrasts_df.to_dict(orient="list")
        job_spec["lrt_tests"] = {"name": [], "reduced_terms": []}
    else:
        job_spec["contrasts"] = {"name": [], "column": [], "level1": [], "level2": []}
        # lrt_tests has a variable-length nested list (reduced_terms per
        # test), so it's built as a plain list of dicts rather than
        # forced through pandas' to_dict(orient="list") -- json.dump
        # handles this shape natively, and jsonlite on the R side parses
        # a JSON array of objects with a ragged array field into a
        # data.frame with a list-column, which the R script's
        # job$lrt_tests$reduced_terms[[i]] indexing expects.
        job_spec["lrt_tests"] = lrt_tests or []

    job_spec_path = os.path.join(work_dir, "deseq2_job_spec.json")
    with open(job_spec_path, "w") as f:
        json.dump(job_spec, f, indent=2)

    r_script_path = os.path.join(work_dir, "run_deseq2.R")
    with open(r_script_path, "w") as f:
        f.write(_DESEQ2_R_SCRIPT)

    cmd = ["Rscript", r_script_path, job_spec_path]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, timeout=3600
        )
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, f"DESeq2 failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, "DESeq2 timed out after 1 hour."


# ---------------------------------------------------------------------------
# Reading back DESeq2's output
# ---------------------------------------------------------------------------

def list_contrast_results(output_dir):
    """List available contrast result CSVs in output_dir, returning contrast names."""
    if not os.path.isdir(output_dir):
        return []
    names = []
    for fname in os.listdir(output_dir):
        match = re.match(r"^results_(.+)\.csv$", fname)
        if match:
            names.append(match.group(1))
    return sorted(names)


def read_contrast_results(output_dir, contrast_safe_name):
    """Read a single contrast's results CSV into a DataFrame."""
    path = os.path.join(output_dir, f"results_{contrast_safe_name}.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def read_pca_coordinates(output_dir, batch_adjusted=False):
    """
    Read PCA coordinates (+ percent variance explained per PC) computed
    during the DESeq2 run.

    Returns (pca_df_or_None, pct_variance_list_or_None).
    """
    suffix = "_batch_adjusted" if batch_adjusted else ""
    pca_path = os.path.join(output_dir, f"pca_coordinates{suffix}.csv")
    var_path = os.path.join(output_dir, f"pca_percent_variance{suffix}.txt")
    if not os.path.exists(pca_path):
        return None, None
    pca_df = pd.read_csv(pca_path)
    pct_variance = None
    if os.path.exists(var_path):
        with open(var_path) as f:
            pct_variance = [float(v) for v in f.read().strip().split(",")]
    return pca_df, pct_variance


def read_normalized_counts(output_dir):
    path = os.path.join(output_dir, "normalized_counts.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def compute_batch_variance_metric(pca_df, batch_column):
    """
    For a single PCA coordinates DataFrame (either the raw or the
    batch-adjusted view — see read_pca_coordinates), estimate what
    share of each principal component's own variance is attributable
    to the batch column, as a simple QC diagnostic for batch
    correction (used to build the "before vs. after" comparison shown
    in the Differential Expression workspace).

    For each available PC (up to PC4), this computes a one-way-ANOVA-
    style effect size — the between-batch-group sum of squares divided
    by the total sum of squares for that PC, expressed as a percentage
    (0-100). This is 0% when every batch group has an identical mean
    position on that PC (no batch effect along that axis), and
    approaches 100% when a sample's position on that PC is almost
    fully determined by which batch it's in.

    This is a simple descriptive/QC diagnostic, not a formal
    significance test — it's meant to give an intuitive, comparable
    "how much did this shrink after adjustment" number, not a p-value.

    Returns a dict {"PC1": pct, "PC2": pct, ...}, or an empty dict if
    pca_df is None/missing the batch column.
    """
    if pca_df is None or not batch_column or batch_column not in pca_df.columns:
        return {}

    groups = pca_df[batch_column].astype(str)
    pc_cols = [c for c in pca_df.columns if re.fullmatch(r"PC\d+", c)]

    result = {}
    for pc_col in sorted(pc_cols, key=lambda c: int(c[2:]))[:4]:
        values = pca_df[pc_col].astype(float)
        grand_mean = values.mean()
        total_ss = ((values - grand_mean) ** 2).sum()
        if total_ss == 0:
            result[pc_col] = 0.0
            continue
        between_ss = 0.0
        for _, sub in values.groupby(groups):
            n = len(sub)
            group_mean = sub.mean()
            between_ss += n * (group_mean - grand_mean) ** 2
        result[pc_col] = round(float(100 * between_ss / total_ss), 1)

    return result


# ---------------------------------------------------------------------------
# QC diagnostics: dispersion, size factors, sample distances
# ---------------------------------------------------------------------------

def read_dispersion_estimates(output_dir):
    """
    Read the per-gene dispersion diagnostic table (gene-wise estimate,
    fitted trend, and final/shrunk value) written by the R script,
    for building the dispersion plot (DESeq2's own recommended
    model-fit diagnostic, equivalent to plotDispEsts in R).

    Returns a DataFrame, or None if this project's DESeq2 run predates
    this file being written (older projects can simply re-run DESeq2
    to generate it).
    """
    path = os.path.join(output_dir, "dispersion_estimates.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def assess_dispersion_fit(dispersion_df):
    """
    A simple, descriptive QC check on the dispersion plot: what
    fraction of genes have a gene-wise dispersion estimate far above
    the fitted trend line (in log-space), which -- in large numbers --
    can indicate the fitted trend didn't capture the data well (e.g.
    due to strong unmodeled heterogeneity, contamination, or a
    sample-swap). A modest fraction of genes scattering above the
    trend is completely normal and expected (this is exactly why
    DESeq2 shrinks toward the trend in the first place) -- this check
    is only meant to flag an unusually large fraction, not any
    deviation at all.

    Returns a dict:
        {
            "pct_far_above_trend": float,   # 0-100
            "flagged": bool,                # True if this looks unusual
            "message": str,                 # plain-language summary
        }
    or a dict with flagged=False and an explanatory message if the
    dispersion data isn't usable (e.g. all-NaN).
    """
    if dispersion_df is None or dispersion_df.empty:
        return {"pct_far_above_trend": 0.0, "flagged": False, "message": "No dispersion data available."}

    valid = dispersion_df.dropna(subset=["dispGeneEst", "dispFit"])
    valid = valid[(valid["dispGeneEst"] > 0) & (valid["dispFit"] > 0)]
    if valid.empty:
        return {"pct_far_above_trend": 0.0, "flagged": False, "message": "No valid dispersion estimates to assess."}

    import math
    log_ratio = valid["dispGeneEst"].apply(math.log2) - valid["dispFit"].apply(math.log2)
    # "Far above" the trend: more than 2 log2-units (4x) higher than
    # fitted -- a fairly generous margin, since scatter around the
    # trend is expected; this only catches genes clearly detached from
    # it.
    far_above = (log_ratio > 2).sum()
    pct = round(100 * far_above / len(valid), 1)

    flagged = pct > 15.0
    if flagged:
        message = (
            f"⚠️ {pct}% of genes have a dispersion estimate well above the "
            "fitted trend line. This can happen with strong unmodeled "
            "variability (e.g. an unaccounted-for batch effect, sample "
            "contamination, or a mixed-condition sample) -- check the PCA "
            "and sample distance heatmap above for any samples that don't "
            "cluster as expected."
        )
    else:
        message = (
            f"✅ {pct}% of genes fall well above the fitted trend -- within "
            "the normal range. The dispersion model looks like a "
            "reasonable fit for this dataset."
        )
    return {"pct_far_above_trend": pct, "flagged": flagged, "message": message}


def read_size_factors(output_dir):
    """
    Read the per-sample size factor table (DESeq2's median-of-ratios
    normalization factors) written by the R script.

    Returns a DataFrame with columns "sample", "size_factor", or None
    if this project's DESeq2 run predates this file being written.
    """
    path = os.path.join(output_dir, "size_factors.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def flag_size_factor_outliers(size_factor_df, threshold_multiplier=2.0):
    """
    Flag any sample whose size factor is unusually far from the median
    of all samples in this project -- often an early sign of a failed,
    degraded, or contamination-affected library (very different total
    sequencing depth/composition from its peers), worth double-checking
    before trusting this sample's results.

    threshold_multiplier: a sample is flagged if its size factor is
        more than this many times higher OR lower than the median
        (e.g. 2.0 = flagged if >2x or <0.5x the median). This is a
        simple, descriptive heuristic (not a formal statistical test)
        -- intentionally conservative so it only catches genuinely
        large discrepancies rather than routine sample-to-sample
        variation.

    Returns a dict:
        {
            "median_size_factor": float,
            "flagged_samples": [{"sample": str, "size_factor": float, "ratio_to_median": float}, ...],
        }
    """
    if size_factor_df is None or size_factor_df.empty:
        return {"median_size_factor": None, "flagged_samples": []}

    median_sf = float(size_factor_df["size_factor"].median())
    flagged = []
    if median_sf > 0:
        for _, row in size_factor_df.iterrows():
            ratio = row["size_factor"] / median_sf
            if ratio > threshold_multiplier or ratio < (1.0 / threshold_multiplier):
                flagged.append({
                    "sample": row["sample"],
                    "size_factor": round(float(row["size_factor"]), 3),
                    "ratio_to_median": round(float(ratio), 2),
                })
    return {"median_size_factor": round(median_sf, 3), "flagged_samples": flagged}


def read_sample_distance_matrix(output_dir):
    """
    Read the full sample-to-sample Euclidean distance matrix (computed
    on VST-transformed counts) written by the R script -- the standard
    QC companion to PCA for sample-level quality control.

    Returns a DataFrame indexed by sample name, with one column per
    sample (a square matrix), or None if this project's DESeq2 run
    predates this file being written.
    """
    path = os.path.join(output_dir, "sample_distance_matrix.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df = df.set_index("sample")
    return df


def detect_sample_clustering_mismatch(distance_df, meta_df, group_column):
    """
    For each sample, check whether its single NEAREST neighbor (by the
    sample-to-sample distance matrix) belongs to the same group_column
    value as itself. A sample whose closest match is in a DIFFERENT
    group is a concrete, specific reason to double-check that sample's
    labeling or quality -- much more actionable than just eyeballing a
    PCA plot or heatmap for anything that looks "off".

    distance_df: square DataFrame from read_sample_distance_matrix
        (index and columns both sample names).
    meta_df: the project's metadata DataFrame (must have a "sample"
        column and group_column).
    group_column: which metadata column defines the "expected" grouping
        to check against (e.g. the primary condition column).

    Returns a list of dicts, one per mismatched sample:
        [{"sample": str, "own_group": str, "nearest_neighbor": str,
          "neighbor_group": str}, ...]
    Empty list if every sample's nearest neighbor shares its own group,
    or if inputs are missing/invalid.
    """
    if distance_df is None or meta_df is None or group_column not in meta_df.columns:
        return []

    sample_to_group = dict(zip(meta_df["sample"].astype(str), meta_df[group_column].astype(str)))
    mismatches = []
    for sample in distance_df.index:
        if sample not in sample_to_group:
            continue
        own_group = sample_to_group[sample]
        # Nearest neighbor = smallest distance excluding the sample
        # itself (which is always 0 on the diagonal).
        distances_to_others = distance_df.loc[sample].drop(labels=[sample], errors="ignore")
        if distances_to_others.empty:
            continue
        nearest = distances_to_others.idxmin()
        neighbor_group = sample_to_group.get(nearest, "?")
        if neighbor_group != own_group:
            mismatches.append({
                "sample": sample,
                "own_group": own_group,
                "nearest_neighbor": nearest,
                "neighbor_group": neighbor_group,
            })
    return mismatches


# ---------------------------------------------------------------------------
# QC diagnostics: p-value histogram shape + MA plot bias
# ---------------------------------------------------------------------------

def classify_pvalue_histogram_shape(results_df, n_bins=20):
    """
    Classify the overall shape of a contrast's p-value distribution --
    one of the most commonly recommended, cheapest sanity checks after
    any differential expression run, since specific shapes map to
    specific, well-known problems:
      - "healthy": a spike near 0 with an approximately flat/uniform
        tail elsewhere -- the expected shape when there IS real
        differential signal, mixed with a background of non-DE genes.
      - "flat/uniform": no spike near 0 at all -- consistent with
        little to no real differential expression signal in this
        comparison (not necessarily a problem with the analysis
        itself, but worth knowing before over-interpreting a small
        number of "significant" hits, which may just be false
        positives from multiple testing).
      - "hump in the middle" or "anti-conservative" (more small
        p-values than expected AND a hump away from 0/1): often
        indicates an unmodeled covariate/batch effect distorting the
        test.
      - "spike near 1": often indicates overdispersion the model
        underestimated, or a design/contrast specification issue.

    This is a simple heuristic based on comparing the actual histogram
    bin counts against what a uniform distribution would predict --
    not a formal statistical test -- meant to give an immediate,
    plain-language read rather than requiring the user to interpret
    the shape themselves.

    Returns a dict:
        {
            "shape": str,           # one of the categories above
            "message": str,         # plain-language explanation
            "counts": list[int],    # histogram bin counts, for plotting
            "bin_edges": list[float],
        }
    """
    pvalues = results_df["pvalue"].dropna()
    if len(pvalues) < 10:
        return {
            "shape": "insufficient_data",
            "message": "Not enough tested genes to assess the p-value distribution shape.",
            "counts": [], "bin_edges": [],
        }

    counts, bin_edges = _histogram(pvalues.tolist(), n_bins, 0.0, 1.0)
    n_total = len(pvalues)
    expected_per_bin = n_total / n_bins

    first_bin_count = counts[0]
    last_bin_count = counts[-1]
    # "Tail" bins exclude the first and last bin, used to judge whether
    # the middle of the distribution is roughly flat (as expected for
    # non-DE genes under the null) or itself elevated/humped. Both the
    # AVERAGE across tail bins and the single MOST-elevated tail bin are
    # checked (not just the average) -- a hump concentrated in just 1-2
    # bins near the center (a classic batch-confounding signature) can
    # otherwise be diluted into an unremarkable-looking average once
    # spread across all ~18 tail bins, even though it's visually obvious
    # and diagnostically meaningful.
    tail_counts = counts[1:-1]
    tail_mean = (sum(tail_counts) / len(tail_counts)) if tail_counts else 0
    tail_max = max(tail_counts) if tail_counts else 0

    # Checked in this specific order of precedence: a spike at either
    # end (0 or 1) is checked FIRST and independently of the other, since
    # a real spike at one end can coexist with an otherwise flat middle
    # -- checking "is the middle flat" before "is there a spike at 1"
    # would otherwise misclassify a genuine spike-near-1 pattern as
    # merely "flat_uniform" (the middle bins alone do look flat; only
    # the last bin is elevated). The "hump" check comes last since it
    # specifically describes an elevated middle/background, which is a
    # different, unrelated pattern from an end-spike.
    if last_bin_count > expected_per_bin * 2 and first_bin_count <= expected_per_bin * 1.5:
        shape = "spike_near_one"
        message = (
            "⚠️ There's an unusual spike near p = 1. This pattern often "
            "shows up when the model has underestimated variability "
            "(overdispersion) for some genes, or can indicate an issue "
            "with how this contrast was specified -- worth double-"
            "checking the design/contrast setup above."
        )
    elif first_bin_count > expected_per_bin * 2 and tail_mean <= expected_per_bin * 1.3:
        shape = "healthy"
        message = (
            "✅ This looks like a healthy p-value distribution: a clear "
            "spike near 0 (real differential signal) with a roughly flat "
            "background elsewhere -- the expected pattern for a "
            "comparison with genuine differentially expressed genes."
        )
    elif tail_mean > expected_per_bin * 1.5 or tail_max > expected_per_bin * 4:
        shape = "hump_or_anticonservative"
        message = (
            "⚠️ The middle/background of this distribution looks elevated "
            "rather than flat, which often points to an unmodeled source "
            "of variation (like a batch effect not included in the design "
            "formula) distorting the test. Check the PCA and sample "
            "distance heatmap above for any structure not explained by "
            "your design columns."
        )
    elif first_bin_count <= expected_per_bin * 1.3:
        shape = "flat_uniform"
        message = (
            "ℹ️ This distribution looks roughly flat/uniform, with little "
            "to no spike near 0 -- consistent with **little to no real "
            "differential expression signal** in this comparison. Any "
            "individually \"significant\" genes here may be more likely to "
            "be false positives from multiple testing -- interpret with "
            "extra caution."
        )
    else:
        shape = "healthy"
        message = (
            "✅ This p-value distribution looks reasonable, with no major "
            "red flags."
        )

    return {
        "shape": shape,
        "message": message,
        "counts": counts,
        "bin_edges": bin_edges,
    }


def _histogram(values, n_bins, lo, hi):
    "Minimal dependency-free histogram, avoiding a numpy/scipy requirement for this one calculation."
    width = (hi - lo) / n_bins
    counts = [0] * n_bins
    for v in values:
        if v < lo or v > hi:
            continue
        idx = int((v - lo) / width)
        if idx >= n_bins:
            idx = n_bins - 1
        counts[idx] += 1
    bin_edges = [lo + i * width for i in range(n_bins + 1)]
    return counts, bin_edges


def classify_ma_bias(results_df, low_mean_percentile=25):
    """
    A simple, descriptive check for systematic bias in an MA plot: do
    genes with a LOW mean expression show a meaningfully different
    average fold-change direction/magnitude than genes with a HIGH mean
    expression? A healthy MA plot is roughly symmetric around
    log2FoldChange = 0 across the full range of expression -- a
    systematic drift specifically among low-count genes is a classic
    sign of a normalization problem or a strong compositional effect
    (e.g. a small number of very highly expressed genes dominating the
    size factor estimate).

    low_mean_percentile: genes below this percentile of baseMean are
        treated as "low-expression" for this comparison.

    Returns a dict:
        {
            "low_mean_avg_lfc": float,
            "high_mean_avg_lfc": float,
            "flagged": bool,
            "message": str,
        }
    """
    valid = results_df.dropna(subset=["baseMean", "log2FoldChange"])
    valid = valid[valid["baseMean"] > 0]
    if len(valid) < 20:
        return {
            "low_mean_avg_lfc": None, "high_mean_avg_lfc": None,
            "flagged": False,
            "message": "Not enough genes with valid data to assess MA-plot bias.",
        }

    threshold = valid["baseMean"].quantile(low_mean_percentile / 100.0)
    low_group = valid[valid["baseMean"] <= threshold]
    high_group = valid[valid["baseMean"] > threshold]

    low_avg = float(low_group["log2FoldChange"].mean())
    high_avg = float(high_group["log2FoldChange"].mean())

    # Flag if the low-expression group's average fold-change is both
    # meaningfully non-zero AND clearly different from the
    # high-expression group's average -- a modest asymmetry is normal,
    # this threshold is set to only catch a fairly pronounced drift.
    flagged = abs(low_avg) > 0.5 and abs(low_avg - high_avg) > 0.5

    if flagged:
        message = (
            f"⚠️ Low-expression genes show a systematically different "
            f"average fold-change ({low_avg:+.2f}) than high-expression "
            f"genes ({high_avg:+.2f}). This pattern can indicate a "
            "normalization issue or a strong compositional effect (e.g. a "
            "few very highly expressed genes dominating size factor "
            "estimation) -- check the size factor QC table above for any "
            "outlier samples."
        )
    else:
        message = (
            f"✅ Low-expression genes ({low_avg:+.2f} average log2FC) and "
            f"high-expression genes ({high_avg:+.2f} average log2FC) look "
            "reasonably symmetric around zero -- no signs of a "
            "normalization issue."
        )
    return {
        "low_mean_avg_lfc": round(low_avg, 3),
        "high_mean_avg_lfc": round(high_avg, 3),
        "flagged": flagged,
        "message": message,
    }


# ---------------------------------------------------------------------------
# QC diagnostics: independent filtering + Cook's outlier explanations
# ---------------------------------------------------------------------------

def read_filter_threshold(output_dir, contrast_safe_name):
    """
    Read the independent-filtering mean-count threshold DESeq2 chose
    for a specific contrast/LRT test (the baseMean cutoff below which
    genes are excluded from multiple-testing correction, optimizing the
    number of genes that pass at the target FDR).

    Returns a float, or None if unavailable (e.g. an older project run
    before this file was written, or independent filtering was not
    applicable for this contrast).
    """
    path = os.path.join(output_dir, f"filter_threshold_{contrast_safe_name}.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        content = f.read().strip()
    try:
        return float(content)
    except ValueError:
        return None


def summarize_na_genes(results_df, filter_threshold=None):
    """
    Summarize the two distinct, legitimate reasons a gene can show NA
    values in a DESeq2 results table, using the "cooks_outlier_flagged"
    and "low_count_filtered" columns the R script now tags on every
    results CSV (see _DESEQ2_R_SCRIPT) -- so the workspace can EXPLAIN
    every NA row's cause in plain language rather than leaving the user
    to wonder why a gene they expected to see is missing or blank.

    filter_threshold: optional, the baseMean cutoff from
        read_filter_threshold(), included in the message for context.

    Returns a dict:
        {
            "n_cooks_outliers": int,
            "n_low_count_filtered": int,
            "cooks_outlier_genes": list[str],   # up to 20, for display
            "cooks_outlier_sample_counts": {sample: count},  # which
                samples are most often the outlier driving a flagged gene
            "filter_threshold": float or None,
            "message": str,
        }
    Gracefully handles older results files that predate these columns
    (returns zeros with an explanatory message) rather than raising.
    """
    if "cooks_outlier_flagged" not in results_df.columns or "low_count_filtered" not in results_df.columns:
        return {
            "n_cooks_outliers": 0, "n_low_count_filtered": 0,
            "cooks_outlier_genes": [], "cooks_outlier_sample_counts": {},
            "filter_threshold": filter_threshold,
            "message": (
                "ℹ️ This project's results were generated before outlier/"
                "filtering diagnostics were added -- re-run DESeq2 to see "
                "this breakdown for future analyses."
            ),
        }

    cooks_mask = results_df["cooks_outlier_flagged"].fillna(False)
    filtered_mask = results_df["low_count_filtered"].fillna(False)

    n_cooks = int(cooks_mask.sum())
    n_filtered = int(filtered_mask.sum())

    cooks_genes = results_df.loc[cooks_mask, "gene_id"].astype(str).tolist()

    sample_counts = {}
    if "max_cooks_sample" in results_df.columns:
        sample_series = results_df.loc[cooks_mask, "max_cooks_sample"].dropna()
        sample_counts = sample_series.value_counts().to_dict()

    parts = []
    if n_cooks:
        top_sample_note = ""
        if sample_counts:
            top_sample, top_count = max(sample_counts.items(), key=lambda kv: kv[1])
            top_sample_note = f" -- most often (**{top_count}** gene(s)) traced back to sample **{top_sample}**"
        parts.append(
            f"**{n_cooks:,} gene(s)** have their p-value set to NA because "
            f"a single sample's count looks like it's driving the result "
            f"(Cook's distance outlier){top_sample_note}. This is DESeq2's "
            "own built-in protection against one unusual sample distorting "
            "a gene's result."
        )
    if n_filtered:
        threshold_note = f" (mean normalized count below ~{filter_threshold:.1f})" if filter_threshold else ""
        parts.append(
            f"**{n_filtered:,} gene(s)** have their adjusted p-value "
            f"(padj) set to NA{threshold_note} because DESeq2's independent "
            "filtering excluded them from multiple-testing correction -- "
            "genes with very low average expression have little power to "
            "ever reach significance, so excluding them improves detection "
            "power for the genes that remain. Their raw p-value is still "
            "shown."
        )
    if not parts:
        message = "✅ No genes were flagged as Cook's-distance outliers or excluded by independent filtering in this contrast."
    else:
        message = " ".join(parts)

    return {
        "n_cooks_outliers": n_cooks,
        "n_low_count_filtered": n_filtered,
        "cooks_outlier_genes": cooks_genes[:20],
        "cooks_outlier_sample_counts": sample_counts,
        "filter_threshold": filter_threshold,
        "message": message,
    }


# ---------------------------------------------------------------------------
# Venn diagram set computation (pure Python -- no R needed for set logic)
# ---------------------------------------------------------------------------

def get_significant_gene_sets(output_dir, contrast_names, padj_threshold=0.05,
                                lfc_threshold=None):
    """
    For each of the given contrast names, load its results CSV and
    return the set of significant gene IDs (padj < padj_threshold, and
    optionally |log2FoldChange| >= lfc_threshold), for building a Venn
    diagram of overlap across contrasts.

    Returns a dict {contrast_name: set_of_gene_ids}.
    """
    gene_sets = {}
    for name in contrast_names:
        df = read_contrast_results(output_dir, name)
        if df is None:
            gene_sets[name] = set()
            continue
        mask = df["padj"] < padj_threshold
        if lfc_threshold is not None:
            mask &= df["log2FoldChange"].abs() >= lfc_threshold
        gene_sets[name] = set(df.loc[mask, "gene_id"])
    return gene_sets


def compute_venn_regions(gene_sets):
    """
    Given a dict {name: set_of_genes} with 2 or 3 sets, compute every
    overlap region needed to draw a Venn diagram (e.g. for 2 sets: only
    A, only B, A and B; for 3 sets: all 7 regions).

    Returns a dict {region_label: count}, where region_label is e.g.
    "A", "B", "A & B" for 2 sets, or similarly for 3 sets. Also returns
    the region gene ID sets themselves under the same labels, as a
    second dict, so the caller can export the actual overlapping genes
    if needed (e.g. for further clusterProfiler analysis on just the
    shared/unique gene lists).

    Only supports 2 or 3 sets, matching typical Venn diagram use cases --
    for more than 3 contrasts, an UpSet-style plot would be more
    appropriate than a Venn diagram, which becomes visually unreadable
    beyond 3 sets.
    """
    names = list(gene_sets.keys())
    if len(names) not in (2, 3):
        raise ValueError("compute_venn_regions only supports 2 or 3 gene sets.")
    if len(names) == 2:
        a, b = names
        set_a, set_b = gene_sets[a], gene_sets[b]
        regions = {
            a: set_a - set_b,
            b: set_b - set_a,
            f"{a} & {b}": set_a & set_b,
        }
    else:
        a, b, c = names
        set_a, set_b, set_c = gene_sets[a], gene_sets[b], gene_sets[c]
        regions = {
            a: set_a - set_b - set_c,
            b: set_b - set_a - set_c,
            c: set_c - set_a - set_b,
            f"{a} & {b}": (set_a & set_b) - set_c,
            f"{a} & {c}": (set_a & set_c) - set_b,
            f"{b} & {c}": (set_b & set_c) - set_a,
            f"{a} & {b} & {c}": set_a & set_b & set_c,
        }
    counts = {label: len(genes) for label, genes in regions.items()}
    return counts, regions


# ---------------------------------------------------------------------------
# clusterProfiler export
# ---------------------------------------------------------------------------

def classify_regulation(results_df, padj_threshold=0.05, lfc_threshold=1.0):
    """
    Classify every gene in a contrast's results as "Up-regulated",
    "Down-regulated", "Not significant", or "NA (not tested)", using
    the same thresholds as the volcano plot -- factored out here as the
    single shared source of truth so the volcano plot's colors/legend
    and any downloadable/displayed table always agree with each other
    exactly, rather than risking two separate copies of this logic
    drifting apart over time.

    "NA (not tested)" is used for rows where padj or log2FoldChange is
    NaN -- this is a real, meaningful DESeq2 outcome (its independent
    filtering step assigns NA to genes it excludes from multiple-
    testing correction, typically very low-mean genes), not the same
    thing as "tested and found not significant", so it's kept as its
    own explicit category rather than being silently folded into
    "Not significant".

    Returns a pandas Series of string labels, aligned to results_df's
    index (so it can be assigned directly as a new column, e.g.
    results_df["Regulation"] = classify_regulation(results_df, ...)).
    """
    padj = results_df["padj"]
    lfc = results_df["log2FoldChange"]

    labels = pd.Series("Not significant", index=results_df.index, dtype="object")
    labels[padj.isna() | lfc.isna()] = "NA (not tested)"

    is_sig = padj < padj_threshold
    labels[is_sig & (lfc >= lfc_threshold)] = "Up-regulated"
    labels[is_sig & (lfc <= -lfc_threshold)] = "Down-regulated"
    # Re-apply the NA mask last, in case a NaN row happened to also
    # satisfy one of the boolean comparisons above (NaN comparisons
    # are always False in pandas, so this is only a defensive
    # safeguard, not expected to actually change anything).
    labels[padj.isna() | lfc.isna()] = "NA (not tested)"

    return labels


def build_pca_qc_export(variance_rows, eta_rows):
    """
    Combine the two small QC summary tables shown in the Differential
    Expression workspace's batch-correction comparison (per-PC percent
    variance explained before/after, and per-PC percent of variance
    attributable to batch before/after) into a single DataFrame
    suitable for a CSV download, merged on "Principal Component".

    variance_rows / eta_rows: lists of dicts, same shape as built in
    differential_expression_workspace.py's render() (each dict has a
    "Principal Component" key plus its own metric columns). Either can
    be an empty list if that particular table wasn't shown (e.g. no
    batch column selected) -- an outer merge is used so no PC present
    in only one of the two is dropped.

    Returns a DataFrame, or an empty DataFrame if both inputs are empty.
    """
    variance_df = pd.DataFrame(variance_rows) if variance_rows else pd.DataFrame(columns=["Principal Component"])
    eta_df = pd.DataFrame(eta_rows) if eta_rows else pd.DataFrame(columns=["Principal Component"])

    if variance_df.empty and eta_df.empty:
        return pd.DataFrame()
    if variance_df.empty:
        return eta_df
    if eta_df.empty:
        return variance_df

    return variance_df.merge(eta_df, on="Principal Component", how="outer")


def build_venn_export(output_dir, contrast_names, regions, display_label_map=None):
    """
    Build a full-detail CSV-ready export for the Venn diagram overlap
    view: one row per gene that appears in at least one of the
    significant gene sets shown on the diagram, with which overlap
    region it belongs to, plus that gene's real DESeq2 statistics
    (log2FoldChange, padj) from EVERY selected contrast -- not just a
    plain list of gene IDs -- so the export is genuinely useful for
    downstream analysis rather than just a name list.

    contrast_names: list of original contrast safe-names (as used by
        read_contrast_results/list_contrast_results) -- i.e. the actual
        file-backed names, regardless of any custom display label the
        user may have set in the UI.
    regions: dict {region_label: set(gene_id)}, as returned by
        compute_venn_regions -- region_label here is expected to use
        whatever *display* labels were already applied (e.g. a user's
        custom-renamed contrast labels), since that's what should
        appear in the "Venn Region" column for consistency with the
        diagram itself.
    display_label_map: optional dict {original_contrast_name:
        display_label} used only to name each contrast's
        log2FoldChange/padj columns using the same label shown in the
        diagram/legend (falls back to the original name for any
        contrast not present in the map).

    Returns a DataFrame with columns: gene_id, Venn Region, then
    "log2FoldChange (<label>)" and "padj (<label>)" for each contrast
    in contrast_names -- or an empty DataFrame if no genes appear in
    any region.
    """
    display_label_map = display_label_map or {}

    gene_to_region = {}
    for label, genes in regions.items():
        for gene_id in genes:
            gene_to_region[gene_id] = label

    if not gene_to_region:
        return pd.DataFrame()

    export_df = pd.DataFrame({"gene_id": sorted(gene_to_region.keys())})
    export_df["Venn Region"] = export_df["gene_id"].map(gene_to_region)

    for name in contrast_names:
        res = read_contrast_results(output_dir, name)
        if res is None:
            continue
        display_label = display_label_map.get(name, name)
        res_subset = res[["gene_id", "log2FoldChange", "padj"]].rename(columns={
            "log2FoldChange": f"log2FoldChange ({display_label})",
            "padj": f"padj ({display_label})",
        })
        export_df = export_df.merge(res_subset, on="gene_id", how="left")

    return export_df


def build_clusterprofiler_export(output_dir, contrast_safe_name):
    """
    Build a clean, ranked gene list from a contrast's results, in the
    format clusterProfiler's enrichGO/gseGO functions expect: a gene ID
    column and a ranking metric (here, log2FoldChange, with genes sorted
    by it -- a common ranking choice for GSEA-style analysis; padj is
    also included for filtering to significant genes only for standard
    over-representation analysis via enrichGO).

    Note: clusterProfiler's enrichGO/gseGO typically expect Entrez IDs
    or gene symbols, not Ensembl transcript/gene IDs directly, depending
    on the organism annotation package (org.Hs.eg.db, org.Mm.eg.db,
    etc.) being used -- converting between ID types (e.g. via
    clusterProfiler's own bitr() function or biomaRt) is left to the
    user's downstream R analysis, since the correct annotation package
    depends on their organism and isn't something this tool can
    determine automatically.

    Returns a DataFrame with columns: gene_id, log2FoldChange, padj --
    sorted by log2FoldChange descending, or None if the contrast wasn't
    found.
    """
    df = read_contrast_results(output_dir, contrast_safe_name)
    if df is None:
        return None
    export_df = df[["gene_id", "log2FoldChange", "padj"]].dropna(subset=["log2FoldChange"])
    export_df = export_df.sort_values("log2FoldChange", ascending=False)
    return export_df
