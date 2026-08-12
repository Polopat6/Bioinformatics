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
    res <- results(dds, contrast = c(contrast_column, contrast_level1, contrast_level2))
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    res_df <- res_df[order(res_df$padj), ]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)
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
    res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    res_df <- res_df[order(res_df$padj), ]
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", test_name)
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


def compute_batch_variance_metric(pca_df, batch_column, pcs=("PC1", "PC2")):
    """
    Compute a simple eta-squared-style QC metric per principal
    component: the fraction of that PC's total variance explained by
    batch group membership (between-batch sum of squares / total sum of
    squares), expressed as a percentage.

    This is NOT a formal statistical test (no p-value, no correction
    for other design factors) -- it's a lightweight, at-a-glance
    indicator meant purely to help the user compare, side by side, how
    much batch appears to be driving sample separation on each axis
    *before* vs. *after* the visualization-only batch-adjusted PCA
    view. The actual DESeq2 statistical test already accounts for
    batch correctly via the design formula regardless of what this
    metric shows.

    pca_df: a PCA coordinates DataFrame as returned by
        read_pca_coordinates() (must include the batch_column and each
        column named in pcs).
    batch_column: name of the metadata column identifying batch group
        membership.
    pcs: which PC columns to compute the metric for.

    Returns a dict {pc_name: eta_squared_percentage}, or an empty dict
    if batch_column isn't present in pca_df.
    """
    results = {}
    if pca_df is None or batch_column not in pca_df.columns:
        return results
    for pc in pcs:
        if pc not in pca_df.columns:
            continue
        grand_mean = pca_df[pc].mean()
        ss_total = ((pca_df[pc] - grand_mean) ** 2).sum()
        if ss_total == 0:
            results[pc] = 0.0
            continue
        ss_between = 0.0
        for _, group in pca_df.groupby(batch_column):
            group_mean = group[pc].mean()
            ss_between += len(group) * (group_mean - grand_mean) ** 2
        results[pc] = round(100 * ss_between / ss_total, 1)
    return results


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
