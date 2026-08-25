"""
single_cell/sc_downstream_manager.py

Backend for Phase 3 of the Single-cell RNA-Seq pipeline -- multi-sample
downstream analysis (normalization, HVG selection, PCA, batch
correction/integration, clustering, UMAP/t-SNE, cell-type annotation,
pseudobulk aggregation, and compositional analysis), invoked from
singlecell_workspace.py (or a dedicated sc_downstream_workspace.py --
see that module for the UI layer).

--- Architectural shift from Phase 1/2 (2026-08-20) ---
Phase 1 (ingestion through alignment) and Phase 2 (Cell-level QC) both
operate strictly PER-SAMPLE -- each sample gets its own STARsolo run,
its own Cell-level QC run, its own output directory. Phase 3 is
fundamentally different: normalization, batch correction, clustering,
and every downstream comparison in this module operate on a SINGLE,
COMBINED object spanning every sample in the project at once. This
module's very first step (load_and_combine_samples) is the bridge
between these two worlds -- it reads each sample's already-completed
Cell-level QC output (STARsolo's filtered matrix + Phase 2's per-cell
QC metrics/flags) and merges them into one AnnData object, with a new
"sample" column in .obs so every downstream step can still distinguish
which sample each cell originally came from.

--- Tool choice: Scanpy (Python), not R (2026-08-20) ---
Unlike Phase 2 (built on R/Bioconductor packages with no good Python
equivalent -- scDblFinder, DecontX, SoupX), Phase 3's operations
(normalization, HVG, PCA, clustering, UMAP) are Scanpy's own core
competency and match this project's stated preference for Python-
native integration with this Streamlit app. Everything in this module
runs in-process via scanpy/anndata -- no R subprocess, no job-spec
JSON file, no Rscript call.

--- Where real methodological choices exist vs. one-size-fits-all
    (2026-08-20) ---
Mirroring Phase 2's DecontX-vs-SoupX pattern: some Phase 3 steps have
genuinely contested alternatives in the literature (confirmed via
direct research, not assumed), and are exposed as real user choices
with plain-language tradeoff explanations, the same way Phase 2's
ambient-RNA-correction and doublet-detection methods are:

  - Normalization: shifted-log (default) vs. Pearson residuals.
    A 2023 benchmark (Ahlmann-Eltze & Huber) comparing 22 transformations
    found the simple shifted-log method OUTPERFORMS SCTransform (the
    gold-standard default in Seurat-based pipelines) on a range of
    tasks -- so shifted-log is the default here, not an afterthought.
    Scanpy's analytic Pearson residuals (Lause et al. 2021) is offered
    as the alternative -- a pure-Python method designed to replicate
    SCTransform's variance-stabilizing approach, since this project's
    environment.yml does not include R's actual SCTransform package
    (would require a heavy new R<->AnnData round-trip dependency chain
    for a method that benchmarks worse than the simpler default anyway).

  - Batch correction/integration: Harmony (default) vs. scVI (NOT
    available in this environment -- see environment.yml's own note:
    deliberately deferred due to its substantial PyTorch dependency
    chain, added only once actually implemented) vs. "skip" (for
    single-sample projects, where there is nothing to integrate).
    Independent benchmarks (scIB, Luecken et al.) show no universal
    winner between Harmony/scVI/Seurat-CCA -- performance depends on
    dataset heterogeneity, so this is a genuine user choice, not a
    settled question.

  - Clustering: Leiden (default, current best practice for graph
    connectivity) vs. Louvain (still offered -- a rigorous 2025
    benchmark spanning 34.9M parameter combinations found NO
    significant performance difference between the two algorithms).
    Both run through the same python-igraph backend already in this
    environment (scanpy's own `flavor="igraph"` louvain implementation
    -- no separate `louvain` package needed).

  - Embeddings: UMAP (default) + t-SNE, both offered in 2D and 3D
    (Plotly Scatter3d) per explicit user request.

--- QC filter toggle (2026-08-20) ---
load_and_combine_samples() accepts an explicit apply_qc_filters flag
(default True) rather than silently baking in "drop adaptive_qc_fail
and doublet-flagged cells" -- standard practice, but the user
explicitly requested this remain a toggle rather than a hardcoded
behavior, mirroring this project's consistent design philosophy of
exposing real analytical choices rather than making them silently.

--- External per-cell metadata extension point (2026-08-23) ---
load_sample_as_anndata() now ALSO checks for an optional
"external_cell_metadata.csv" file inside cellqc_output_dir (alongside
the standard cell_qc_metrics.csv Phase 2 itself generates) and, if
present, merges its columns into .obs the exact same way. This is a
narrow, additive extension point for samples whose donor/condition
identity was NOT determined by this pipeline's own per-sample workflow
at all -- e.g. a pooled multi-donor 10x library that was split into
per-donor "samples" using an EXTERNALLY-computed barcode-to-donor
assignment (see sc_genetic_demux_manager.py), such as the demuxlet-
based donor/condition labels published alongside the Kang et al. 2018
IFN-beta PBMC dataset (GSE96583). Unlike cell_qc_metrics.csv (always
regenerated fresh by every Phase 2 run, and therefore never a safe
place to stash externally-sourced columns), this second file is never
written or overwritten by anything in this pipeline itself -- it is
either absent (the common case: an ordinary sample with no external
donor/condition data) or present-and-untouched (placed there once, by
whatever process performed the external split), so merging it here is
purely additive and changes nothing about how any EXISTING sample
without such a file behaves. A genuine column-name collision against
an already-present cell_qc_metrics.csv column is treated as a real
data-preparation error worth surfacing, not silently overwritten either
way -- see this function's own inline comment at the merge site.

--- Gene-symbol resolution overlay (2026-08-25, FIXED same day after
    real HPC diagnosis) ---
A real, confirmed bug: adata.var["gene_symbol"] was previously set
DIRECTLY from STARsolo's own features.tsv column 2 (the gene name
STAR itself baked into the index at build time) with no further
resolution attempted anywhere in this module -- meaning
find_cluster_markers(), get_dotplot_data(), get_violin_plot_data(),
get_feature_plot_data(), and get_marker_heatmap_data() ALL inherited
whatever features.tsv happened to contain. For any sample whose STAR
index was built BEFORE reference_manager.py's own tiered gene-name
backfill fix (Name -> Ensembl description -> bare Ensembl ID, "gene:"
GFF3-ID-attribute prefix stripped -- see that module's own docstring,
"Tiered gene-name backfill for the GFF3 fallback path"), features.tsv
itself still contains the OLD, unhelpful raw "gene:ENSG..." form as
its gene-name column, and this module had no way to recover a better
name after the fact -- a real reported case confirmed this exact
symptom directly in Step 8's own cluster marker gene table.

_resolve_better_gene_symbols() below fixes this WITHOUT requiring a
full re-alignment (rebuilding the STAR index + re-running Step 6 +
Phase 2 Cell-level QC again for an already-completed sample) -- it
reuses reference_manager.py's own extract_gene_symbol_map_from_gtf()
(already implemented there, not reimplemented here -- avoiding yet
another independent copy of the same GTF-parsing logic, the exact
class of bug this project has been bitten by more than once already)
to look up a real gene name directly from the project's CURRENT
reference GTF, and overlays it OVER (not instead of) STARsolo's own
symbol -- but ONLY for a gene whose EXISTING symbol still looks
unresolved (identical to its own gene_id, or still carrying the raw
"gene:" prefix). A gene that already has a genuinely good symbol from
features.tsv is left completely untouched, and a gene with STILL no
better name available even in the current GTF (e.g. a genuinely
unannotated/novel locus) is correctly left as-is rather than having a
name fabricated for it.

--- REAL BUG FOUND AND FIXED, same day (2026-08-25), via direct HPC
    diagnosis on the user's own real reference/data ---
The FIRST version of _resolve_better_gene_symbols() looked up each
BARE gene_id directly against
reference_manager.extract_gene_symbol_map_from_gtf()'s own returned
dict -- but that dict's OWN keys are whatever literal string appears
in the GTF's own `gene_id "..."` attribute, which is NOT guaranteed to
be the bare, unprefixed form. Confirmed directly via `grep` on the
user's real, actual 1GB human.annotation.gtf (a GFF3-fallback-derived
reference -- see reference_manager.py's own module docstring for that
fallback path): every gene_id in this real file is written as
"gene:ENSG00000108691" -- i.e. the SAME GFF3 "gene:" ID-attribute
prefix already known to cause confusion elsewhere in this project (see
gff3_gene_name_resolver.py's own docstring) -- NOT the bare
"ENSG00000108691" form that STARsolo's features.tsv (and therefore
this pipeline's own adata.var_names/gene_id convention) actually uses.
Because of this, EVERY lookup in the original implementation silently
missed, and re-running "Combine Samples" with gtf_path supplied
produced NO visible correction at all -- confirmed directly by the
user via a real re-combine attempt that still showed unresolved
"gene:ENSG..." symbols in Step 8's marker table afterward.

Fixed by checking BOTH the bare gene_id AND its "gene:"-prefixed form
when looking up a replacement in gtf_symbol_map -- mirroring the SAME
defensive both-forms check already used once before, for the exact
same underlying reason, in reference_manager.py's own
backfill_gene_names_from_gff3_tiered() (see that function's own
docstring, "ROBUSTNESS" note, added after a similar direct-testing
discovery that gffread's own real GTF output convention for this
attribute could not be assumed either way without checking). This is
the second, independent place in this project where exactly this same
"gene:" prefix ambiguity had to be defensively handled -- worth keeping
in mind if a THIRD such lookup site is ever added later.

load_sample_as_anndata() and load_and_combine_samples() both accept an
optional gtf_path parameter (default None, a complete no-op if
omitted, preserving this module's own prior behavior exactly for any
caller that doesn't pass it) -- the workspace UI layer is responsible
for supplying the project's own confirmed reference GTF path (already
available via sc_project_manager.get_reference_choice(project)
["custom_gtf"]) when calling load_and_combine_samples() in Step 1.
"""
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import reference_manager as ref


# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------

def harmony_available():
    "Check whether harmonypy (Harmony batch correction) is installed."
    try:
        import harmonypy  # noqa: F401
        return True
    except ImportError:
        return False


def scvi_available():
    "Check whether scvi-tools is installed -- NOT included in this project's environment.yml by design (see module docstring); this always returns False today but is checked defensively rather than hardcoded, in case it's added later without this module needing an update."
    try:
        import scvi  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Gene-symbol resolution overlay (2026-08-25, fixed same day) -- see
# module docstring
# ---------------------------------------------------------------------------

def _gene_symbol_looks_unresolved(gene_id, symbol):
    """
    Returns True if `symbol` looks like it was NEVER actually resolved
    to a real gene name -- i.e. it's just the bare gene_id itself
    (STARsolo's features.tsv correctly falling back to gene_id when no
    real name was available at STAR-index build time), OR it still
    carries the confusing raw GFF3 "gene:" ID-attribute prefix (a
    reference indexed BEFORE reference_manager.py's tiered gene-name
    backfill fix, OR one whose gene_id itself uses this prefix
    convention -- see this module's own docstring). See this module's
    own docstring, "Gene-symbol resolution overlay", for the full
    rationale.
    """
    if symbol is None or (isinstance(symbol, float) and pd.isna(symbol)):
        return True
    symbol = str(symbol)
    if symbol == gene_id:
        return True
    if symbol.startswith("gene:"):
        return True
    return False


def _resolve_better_gene_symbols(var_df, gtf_path):
    """
    Overlay a corrected gene_symbol wherever the EXISTING value (as
    loaded directly from STARsolo's own features.tsv, in
    load_sample_as_anndata() below) looks unresolved (see
    _gene_symbol_looks_unresolved() above) -- reusing
    reference_manager.py's own extract_gene_symbol_map_from_gtf() to
    look up a real name directly from the project's CURRENT reference
    GTF, rather than reimplementing GTF-parsing logic here.

    var_df: a DataFrame indexed by gene_id, with an existing
        "gene_symbol" column.
    gtf_path: path to the reference GTF used for this project's
        alignment -- if None, or the file doesn't exist, this function
        is a safe, complete no-op and returns var_df UNCHANGED (the
        exact same object, not a copy) -- so a caller that never
        passes gtf_path at all sees IDENTICAL behavior to before this
        overlay existed.

    --- REAL BUG FIX (2026-08-25), confirmed via direct HPC diagnosis
        ---
    reference_manager.extract_gene_symbol_map_from_gtf()'s own
    returned dict is keyed by whatever literal string the GTF's own
    `gene_id "..."` attribute actually contains -- this is NOT
    guaranteed to be the bare, unprefixed gene ID. Confirmed directly
    (via `grep` on the user's real reference GTF) that a
    GFF3-fallback-derived reference can have EVERY gene_id written in
    the "gene:ENSG..." GFF3-ID-attribute-prefixed form. Since
    var_df's own index (this pipeline's gene_id convention, taken
    directly from STARsolo's features.tsv) is always the BARE form,
    looking up ONLY the bare form against gtf_symbol_map silently
    missed every single gene for this real reference -- the original
    version of this function had exactly this bug, and produced NO
    visible correction at all when tested against real data. Fixed by
    checking BOTH the bare gene_id and its "gene:"-prefixed form
    against gtf_symbol_map below.

    Returns a NEW DataFrame (var_df itself is never mutated in place)
    with "gene_symbol" values corrected wherever a real GTF-derived
    name was found -- a gene with STILL no better name available in
    the GTF either is left exactly as it was (still gene_id, or still
    "gene:..."), never fabricated. A gene whose EXISTING symbol
    already looks genuinely resolved is left completely untouched,
    even if the GTF's own value happens to differ.
    """
    if not gtf_path or not os.path.isfile(gtf_path):
        return var_df

    gtf_symbol_map = ref.extract_gene_symbol_map_from_gtf(gtf_path)
    if not gtf_symbol_map:
        return var_df

    result = var_df.copy()
    for gene_id in result.index:
        current_symbol = result.at[gene_id, "gene_symbol"]
        if _gene_symbol_looks_unresolved(gene_id, current_symbol):
            # --- THE FIX: check BOTH the bare gene_id AND its
            # "gene:"-prefixed form -- gtf_symbol_map's own keys may be
            # EITHER, depending on how the source GTF's gene_id
            # attribute was actually written (confirmed, via real HPC
            # diagnosis, to be the prefixed form for this user's real
            # reference). See this function's own docstring, "REAL BUG
            # FIX", for the full story.
            better_symbol = gtf_symbol_map.get(gene_id)
            if better_symbol is None:
                better_symbol = gtf_symbol_map.get(f"gene:{gene_id}")
            if better_symbol and not _gene_symbol_looks_unresolved(gene_id, better_symbol):
                result.at[gene_id, "gene_symbol"] = better_symbol

    return result


# ---------------------------------------------------------------------------
# 3.0: Multi-sample loading + combination
# ---------------------------------------------------------------------------

DEFAULT_MIN_CELLS_PER_GENE = 3

# Filename (inside cellqc_output_dir, alongside cell_qc_metrics.csv) that
# load_sample_as_anndata() checks for optional EXTERNALLY-sourced per-cell
# metadata (e.g. donor/condition labels from a genetic-demultiplexing
# split) -- see this module's own docstring, "External per-cell metadata
# extension point", for the full rationale.
EXTERNAL_CELL_METADATA_FILENAME = "external_cell_metadata.csv"


def load_sample_as_anndata(starsolo_output_prefix, cellqc_output_dir, sample_name,
                            apply_qc_filters=True, gtf_path=None):
    """
    Load ONE sample's STARsolo filtered matrix + Phase 2 Cell-level QC
    results into an AnnData object, with per-cell QC metrics/flags
    attached as .obs columns.

    starsolo_output_prefix: this sample's STARsolo output prefix (e.g.
        os.path.join(starsolo_dir, sample_name, f"{sample_name}_")) --
        used to locate Solo.out/Gene/filtered/ via the same convention
        already used by starsolo_manager.py's filtered_counts_matrix_dir().
    cellqc_output_dir: this sample's Cell-level QC output directory
        (contains cell_qc_metrics.csv from Phase 2's
        sc_cellqc_manager.run_cellqc_analysis(), and OPTIONALLY
        EXTERNAL_CELL_METADATA_FILENAME -- see this module's own
        docstring, "External per-cell metadata extension point").
    apply_qc_filters: if True (the default), drop cells flagged by
        Phase 2's adaptive QC (adaptive_qc_fail) or predicted as
        doublets (predicted_doublet) -- if False, load every cell in
        the filtered matrix regardless of QC flags, and those flags
        remain available as .obs columns for the user's own inspection/
        filtering later. Exposed as an explicit, user-facing TOGGLE
        (not a hardcoded default) per this module's own docstring.
    gtf_path: OPTIONAL path to this project's current reference GTF --
        if given, overlays a corrected gene_symbol for any gene whose
        STARsolo-provided symbol looks unresolved (see this module's
        own docstring, "Gene-symbol resolution overlay"). Defaults to
        None, a complete no-op preserving this function's prior exact
        behavior for any caller that doesn't pass it.

    Returns an AnnData object for this ONE sample, with:
        .obs["sample"] = sample_name (for every cell)
        .obs[qc columns from cell_qc_metrics.csv, if available]
        .obs[columns from external_cell_metadata.csv, if present --
             see this module's own docstring]
        .var_names = gene IDs (raw counts in .X)

    Raises FileNotFoundError if the STARsolo filtered matrix directory
    doesn't exist (caller should check this and surface a clear error
    rather than let this propagate as a generic traceback).
    """
    filtered_dir = os.path.join(starsolo_output_prefix + "Solo.out", "Gene", "filtered")
    if not os.path.isdir(filtered_dir):
        raise FileNotFoundError(
            f"STARsolo filtered matrix directory not found for sample '{sample_name}': {filtered_dir}"
        )

    adata = sc.read_mtx(os.path.join(filtered_dir, "matrix.mtx.gz") if os.path.isfile(os.path.join(filtered_dir, "matrix.mtx.gz")) else os.path.join(filtered_dir, "matrix.mtx")).T

    barcodes_path = os.path.join(filtered_dir, "barcodes.tsv.gz")
    if not os.path.isfile(barcodes_path):
        barcodes_path = os.path.join(filtered_dir, "barcodes.tsv")
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].tolist()

    features_path = os.path.join(filtered_dir, "features.tsv.gz")
    if not os.path.isfile(features_path):
        features_path = os.path.join(filtered_dir, "features.tsv")
    features_df = pd.read_csv(features_path, header=None, sep="\t")
    gene_ids = features_df[0].tolist()
    gene_symbols = features_df[1].tolist() if features_df.shape[1] > 1 else gene_ids

    adata.obs_names = barcodes
    adata.var_names = gene_ids
    adata.var["gene_symbol"] = gene_symbols
    adata.obs["sample"] = sample_name

    # --- Gene-symbol resolution overlay (2026-08-25) -- see this
    # module's own docstring for the full rationale. A complete no-op
    # if gtf_path is None (the default), preserving this function's
    # exact prior behavior for any caller that doesn't pass it.
    adata.var = _resolve_better_gene_symbols(adata.var, gtf_path)

    qc_metrics_path = os.path.join(cellqc_output_dir, "cell_qc_metrics.csv")
    if os.path.isfile(qc_metrics_path):
        qc_df = pd.read_csv(qc_metrics_path).set_index("barcode")
        # Reindex to this sample's actual barcodes (in the SAME order as
        # adata.obs_names) -- guards against any barcode present in the
        # QC file but not in this specific matrix load (shouldn't happen
        # in practice, since both come from the same filtered matrix, but
        # defensive alignment is cheap and avoids a silent length
        # mismatch if it ever does).
        qc_df = qc_df.reindex(adata.obs_names)
        for col in qc_df.columns:
            adata.obs[col] = qc_df[col].values

        if apply_qc_filters:
            keep_mask = pd.Series(True, index=adata.obs_names)
            if "adaptive_qc_fail" in adata.obs.columns:
                keep_mask &= ~adata.obs["adaptive_qc_fail"].fillna(False).astype(bool)
            if "predicted_doublet" in adata.obs.columns:
                keep_mask &= ~adata.obs["predicted_doublet"].fillna(False).astype(bool)
            adata = adata[keep_mask].copy()

    # --- External per-cell metadata (2026-08-23) -- see this module's own
    # docstring, "External per-cell metadata extension point", for the
    # full rationale. Merged AFTER Phase 2's own QC-filter subsetting
    # above, so this reindex is always against adata's CURRENT (possibly
    # already QC-filtered) obs_names -- exactly mirroring how the
    # cell_qc_metrics.csv merge above is itself keyed to adata's own
    # barcodes at the time it runs.
    external_metadata_path = os.path.join(cellqc_output_dir, EXTERNAL_CELL_METADATA_FILENAME)
    if os.path.isfile(external_metadata_path):
        ext_df = pd.read_csv(external_metadata_path)
        barcode_col = "barcode" if "barcode" in ext_df.columns else ext_df.columns[0]
        ext_df = ext_df.set_index(barcode_col)
        # Reindex to adata's CURRENT barcodes -- any barcode present in
        # this external file but not in adata (e.g. already dropped by
        # QC filtering above, or simply never assigned an external label
        # at all) is silently absent from the join, matching the same
        # reindex-based merge pattern used for cell_qc_metrics.csv above.
        # Any adata barcode with NO row in this external file gets NaN
        # for these columns -- expected and harmless for columns that
        # are purely optional/informational (e.g. a "cell" reference
        # column), but callers relying on a REQUIRED external column
        # (e.g. a demultiplexing "sample"/"stim" label the whole point
        # of loading this file was to obtain) should check for
        # unexpected NaNs themselves, since this function does not
        # assume any particular column is mandatory.
        ext_df = ext_df.reindex(adata.obs_names)
        for col in ext_df.columns:
            # Only add a column if it isn't ALREADY present from
            # cell_qc_metrics.csv above -- an external file should never
            # silently overwrite this pipeline's own QC-derived columns;
            # a genuine name collision is a caller/data-preparation
            # mistake worth surfacing, not silently resolving one way.
            if col not in adata.obs.columns:
                adata.obs[col] = ext_df[col].values

    adata.var_names_make_unique()
    return adata


def load_and_combine_samples(sample_specs, apply_qc_filters=True,
                              min_cells_per_gene=DEFAULT_MIN_CELLS_PER_GENE,
                              gtf_path=None):
    """
    Load and combine MULTIPLE samples into a single project-level
    AnnData object -- the bridge between Phase 1/2's per-sample world
    and Phase 3's multi-sample downstream analysis. See this module's
    own docstring, "Architectural shift from Phase 1/2", for the full
    rationale.

    sample_specs: list of dicts, one per sample to include:
        [{"sample_name": str, "starsolo_output_prefix": str,
          "cellqc_output_dir": str}, ...]
        Callers (the UI layer) are responsible for deciding WHICH
        samples to include here -- e.g. offering a multiselect of every
        sample that has completed BOTH Step 6 alignment AND Phase 2
        Cell-level QC, since a sample missing either has nothing valid
        to combine.

    apply_qc_filters: see load_sample_as_anndata()'s own docstring --
        applied identically and independently to EACH sample before
        combining (not applied once to the combined object), so a
        per-sample QC threshold computed by Phase 2 (which is itself
        sample-specific, per this project's adaptive MAD-based design)
        is honored correctly per sample.

    min_cells_per_gene: after combining, drop genes detected in fewer
        than this many cells ACROSS THE WHOLE COMBINED dataset (a
        standard, lightweight gene-level filter -- distinct from Phase
        2's per-CELL filtering, which already happened above). Genes
        passing this bar in even one sample but not reaching it in the
        combined multi-sample total are dropped here, since a gene
        detected in only 1-2 cells total contributes essentially no
        information to HVG selection, PCA, or clustering, and mostly
        just adds computational overhead and sparsity.

    gtf_path: OPTIONAL path to this project's current reference GTF --
        passed straight through to load_sample_as_anndata() for EACH
        sample (see this module's own docstring, "Gene-symbol
        resolution overlay"). Defaults to None, a complete no-op
        preserving this function's prior exact behavior for any caller
        that doesn't pass it.

    Returns (combined_adata: AnnData, per_sample_summary: list[dict]).
    per_sample_summary is a list of {"sample_name", "n_cells_loaded",
    "n_cells_after_qc_filter"} dicts -- for display in the UI so a user
    can see exactly how many cells each sample contributed and how many
    were dropped by the QC filter toggle, BEFORE any further analysis
    runs.

    Raises FileNotFoundError (propagated from load_sample_as_anndata)
    if any sample's STARsolo output can't be found -- callers should
    catch this per-sample if they want to skip a broken sample and
    continue with the rest, rather than failing the whole combination.
    """
    if not sample_specs:
        raise ValueError("At least one sample must be provided to combine.")

    adatas = []
    per_sample_summary = []
    for spec in sample_specs:
        sample_name = spec["sample_name"]
        # Load WITHOUT the QC filter first, to get a true "cells loaded"
        # count, then apply the filter separately so both counts are
        # available for the summary -- avoids loading the matrix twice.
        raw_adata = load_sample_as_anndata(
            spec["starsolo_output_prefix"], spec["cellqc_output_dir"], sample_name,
            apply_qc_filters=False, gtf_path=gtf_path,
        )
        n_loaded = raw_adata.n_obs

        if apply_qc_filters and "adaptive_qc_fail" in raw_adata.obs.columns:
            keep_mask = pd.Series(True, index=raw_adata.obs_names)
            if "adaptive_qc_fail" in raw_adata.obs.columns:
                keep_mask &= ~raw_adata.obs["adaptive_qc_fail"].fillna(False).astype(bool)
            if "predicted_doublet" in raw_adata.obs.columns:
                keep_mask &= ~raw_adata.obs["predicted_doublet"].fillna(False).astype(bool)
            filtered_adata = raw_adata[keep_mask].copy()
        else:
            filtered_adata = raw_adata

        per_sample_summary.append({
            "sample_name": sample_name,
            "n_cells_loaded": n_loaded,
            "n_cells_after_qc_filter": filtered_adata.n_obs,
        })
        adatas.append(filtered_adata)

    if len(adatas) == 1:
        combined = adatas[0]
    else:
        combined = ad.concat(
            adatas, join="outer", label="_concat_batch", index_unique="-",
            fill_value=0, merge="same",
        )
        # merge="same" (NOT anndata's own default of None) -- confirmed
        # via direct testing that ad.concat's DEFAULT merge behavior
        # silently DROPS every .var column (including gene_symbol, set
        # per-sample in load_sample_as_anndata above) unless explicitly
        # told to keep columns that are IDENTICAL across every input
        # object. Since every sample in a project shares the exact same
        # reference/GTF download (confirmed: features.tsv is produced
        # once per reference, not per sample), gene_symbol SHOULD be
        # identical for every gene_id across every sample -- "same" is
        # therefore the correct merge strategy here, not an approximation.
        # Without this fix, gene_symbol would silently vanish from the
        # combined object, breaking CellTypist (which requires gene
        # SYMBOLS, not IDs) and every biologist-facing plot that should
        # show a real gene name instead of a raw Ensembl ID.
        #
        # ad.concat's own "_concat_batch" column duplicates our own
        # "sample" column (already set per-sample above) -- drop it to
        # avoid two redundant, differently-named copies of the same
        # information sitting in .obs.
        if "_concat_batch" in combined.obs.columns:
            combined.obs.drop(columns=["_concat_batch"], inplace=True)

    if min_cells_per_gene and min_cells_per_gene > 0:
        sc.pp.filter_genes(combined, min_cells=min_cells_per_gene)

    combined.layers["counts"] = combined.X.copy()

    return combined, per_sample_summary


# ---------------------------------------------------------------------------
# 3.1: Normalization
# ---------------------------------------------------------------------------

NORMALIZATION_METHOD_OPTIONS = {
    "shifted_log": {
        "label": "Shifted log (recommended default)",
        "explanation": (
            "Scales each cell's counts by a size factor (correcting for sequencing depth "
            "differences between cells), then applies a log1p transform. A rigorous 2023 "
            "benchmark comparing 22 different normalization/transformation methods "
            "(Ahlmann-Eltze & Huber) found this simple approach actually OUTPERFORMS "
            "SCTransform (the more complex, widely-used default in Seurat-based pipelines) "
            "across a range of downstream tasks -- so this is the recommended default here, "
            "not a simplified fallback."
        ),
    },
    "pearson_residuals": {
        "label": "Analytic Pearson residuals",
        "explanation": (
            "A statistically-motivated alternative (Lause et al. 2021) that models each "
            "gene's counts with a negative binomial distribution and computes standardized "
            "residuals -- designed to replicate SCTransform's variance-stabilizing "
            "properties in pure Python, without needing R/SCTransform itself. Can better "
            "preserve biological variance for genes with very different expression "
            "magnitudes, at the cost of being less battle-tested than the shifted-log "
            "default above."
        ),
    },
}
DEFAULT_NORMALIZATION_METHOD = "shifted_log"
DEFAULT_TARGET_SUM = None  # None -> scanpy uses the median total count across cells


def normalize(adata, method=DEFAULT_NORMALIZATION_METHOD, target_sum=DEFAULT_TARGET_SUM,
              pearson_theta=100):
    """
    Apply the chosen normalization method to adata.X (expects raw
    counts -- callers should pass the freshly-combined adata from
    load_and_combine_samples(), before any other transformation).

    method: "shifted_log" or "pearson_residuals" -- see
        NORMALIZATION_METHOD_OPTIONS for the full explanation of each.
    target_sum: only used for "shifted_log" -- the total count each
        cell is normalized to before log1p. None (the default) uses
        scanpy's own default of the median total count across cells in
        THIS dataset.
    pearson_theta: only used for "pearson_residuals" -- the negative
        binomial overdispersion parameter (scanpy's own default of 100
        is used unless overridden).

    Modifies adata in place (matching scanpy's own convention) and
    additionally stores adata.layers["lognorm"] (or
    adata.layers["pearson_residuals"]) as an explicit, named copy of
    the normalized result -- so a later step (e.g. re-running HVG
    selection with a different method) can always fall back to
    adata.layers["counts"] (set by load_and_combine_samples) without
    needing to reload/re-combine from scratch.

    Returns adata (for convenient chaining), though it is also modified
    in place.
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if method == "shifted_log":
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X.copy()
    elif method == "pearson_residuals":
        # Pearson residuals expects RAW counts as input -- ensure we're
        # operating on the untouched counts layer, not any prior
        # transformation, regardless of what adata.X currently holds.
        adata.X = adata.layers["counts"].copy()
        sc.experimental.pp.normalize_pearson_residuals(adata, theta=pearson_theta)
        adata.layers["pearson_residuals"] = adata.X.copy()
    else:
        raise ValueError(f"Unknown normalization method: {method!r}")

    return adata


# ---------------------------------------------------------------------------
# 3.2: Highly variable gene (HVG) selection
# ---------------------------------------------------------------------------

DEFAULT_N_TOP_GENES = 2000


def select_highly_variable_genes(adata, n_top_genes=DEFAULT_N_TOP_GENES, batch_key=None):
    """
    Flag the most informative genes for downstream steps (PCA,
    clustering) -- expects adata.X to already be normalized (see
    normalize() above).

    batch_key: optional .obs column name (typically "sample") -- when
        given, HVG selection is computed PER BATCH and then combined
        (scanpy's own batch-aware ranking), so a gene that's highly
        variable in only one sample doesn't get selected purely due to
        a batch-specific effect. Recommended whenever a project
        combines 2+ samples; leave None for a single-sample project.

    Modifies adata in place (adds adata.var["highly_variable"] and
    related columns, following scanpy's own convention). Returns adata
    for convenient chaining.
    """
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    return adata


def get_hvg_plot_data(adata):
    """
    Extract the data needed to reproduce scanpy's own HVG diagnostic
    plot (mean expression vs. dispersion, with selected HVGs
    highlighted) as a plain DataFrame -- kept separate from any
    specific plotting library so the UI layer's Plotly-based
    customizable-plot pattern (matching every other plot in this app)
    can build the actual figure however it wants.

    Returns a DataFrame with columns: gene_id, mean, dispersion (or
    "variance" for the pearson_residuals path, if that column exists
    instead), highly_variable (bool).

    Returns None if HVG selection hasn't been run yet (adata.var has no
    "highly_variable" column).
    """
    if "highly_variable" not in adata.var.columns:
        return None

    df = pd.DataFrame({"gene_id": adata.var_names})
    for col in ("means", "dispersions", "dispersions_norm", "variances", "variances_norm"):
        if col in adata.var.columns:
            df[col] = adata.var[col].values
    df["highly_variable"] = adata.var["highly_variable"].values
    return df


# ---------------------------------------------------------------------------
# 3.3: PCA
# ---------------------------------------------------------------------------

DEFAULT_N_PCS = 50


def run_pca(adata, n_comps=DEFAULT_N_PCS, use_highly_variable=True, random_state=0):
    """
    Run PCA on adata.X (expects normalized data -- see normalize()
    above), restricted to highly-variable genes by default (standard
    practice -- see select_highly_variable_genes() above; must be run
    first if use_highly_variable=True).

    Modifies adata in place (adds adata.obsm["X_pca"],
    adata.varm["PCs"], adata.uns["pca"]["variance_ratio"], following
    scanpy's own convention). Returns adata for convenient chaining.
    """
    if use_highly_variable and "highly_variable" not in adata.var.columns:
        raise ValueError(
            "use_highly_variable=True requires HVG selection to be run first "
            "(select_highly_variable_genes()) -- adata.var has no 'highly_variable' column yet."
        )
    mask_var = "highly_variable" if use_highly_variable else None
    sc.pp.pca(adata, n_comps=n_comps, mask_var=mask_var, random_state=random_state)
    return adata


def get_pca_variance_ratio(adata):
    """
    Return the fraction of variance explained by each principal
    component, as a plain list of floats -- for building the "elbow
    plot" (variance explained per PC, helping a user choose how many
    PCs to carry forward into batch correction/clustering).

    Returns None if PCA hasn't been run yet.
    """
    if "pca" not in adata.uns or "variance_ratio" not in adata.uns["pca"]:
        return None
    return list(adata.uns["pca"]["variance_ratio"])


def get_pca_coordinates(adata, n_pcs_to_include=10, color_by_columns=None):
    """
    Extract PCA coordinates (PC1, PC2, PC3, ...) as a plain DataFrame,
    one row per cell, along with any requested .obs columns for
    coloring.

    color_by_columns: list of .obs column names to include (e.g.
        ["sample", "condition"]) -- if None, includes every non-PC
        column already in .obs.

    Returns None if PCA hasn't been run yet.
    """
    if "X_pca" not in adata.obsm:
        return None

    n_pcs_available = adata.obsm["X_pca"].shape[1]
    n_pcs_to_include = min(n_pcs_to_include, n_pcs_available)

    df = pd.DataFrame(
        adata.obsm["X_pca"][:, :n_pcs_to_include],
        columns=[f"PC{i+1}" for i in range(n_pcs_to_include)],
        index=adata.obs_names,
    )
    df["cell_barcode"] = adata.obs_names

    cols_to_add = color_by_columns if color_by_columns is not None else list(adata.obs.columns)
    for col in cols_to_add:
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values

    return df


# ---------------------------------------------------------------------------
# 3.4: Batch correction / integration
# ---------------------------------------------------------------------------

BATCH_CORRECTION_METHOD_OPTIONS = {
    "none": {
        "label": "Skip (single-sample project, or no correction desired)",
        "explanation": (
            "No batch correction is applied -- appropriate for a single-sample project "
            "(nothing to integrate), or if you specifically want to see how much of your "
            "samples' separation is attributable to real biological differences vs. batch, "
            "unadjusted."
        ),
    },
    "harmony": {
        "label": "Harmony (recommended default for multi-sample projects)",
        "explanation": (
            "Iteratively adjusts each cell's PCA coordinates to remove sample/batch-driven "
            "separation while preserving genuine biological variation. Fast, well-established, "
            "and a strong general-purpose default -- independent benchmarks (scIB/Luecken et "
            "al.) show no universal winner among Harmony/scVI/Seurat's own integration methods, "
            "so this is a reasonable starting choice, not a settled 'best' answer for every "
            "dataset."
        ),
    },
    # scVI intentionally NOT offered as a real option here -- NOT installed
    # in this project's environment (see environment.yml's own explicit
    # note: deferred due to its substantial PyTorch dependency chain,
    # to be added only once this specific path is actually implemented).
    # scvi_available() below exists so the UI layer can defensively check
    # and hide/disable this option entirely rather than offer a choice
    # that would immediately fail.
}
DEFAULT_BATCH_CORRECTION_METHOD = "harmony"


def run_batch_correction(adata, method=DEFAULT_BATCH_CORRECTION_METHOD, batch_key="sample",
                          basis="X_pca", adjusted_basis="X_pca_harmony", **harmony_kwargs):
    """
    Apply the chosen batch-correction/integration method to adata,
    working on the existing PCA embedding (run_pca() must be called
    first).

    method: "none" or "harmony" -- see BATCH_CORRECTION_METHOD_OPTIONS.
        "none" is a genuine no-op (returns adata completely unchanged) --
        exposed as an explicit choice rather than the caller needing to
        conditionally skip calling this function at all, so UI code can
        always call run_batch_correction() and let the "none" choice
        handle the no-op case uniformly.
    batch_key: the .obs column identifying which sample/batch each cell
        belongs to (typically "sample", the column set by
        load_and_combine_samples()).
    basis: the existing PCA embedding to correct (adata.obsm key) --
        defaults to "X_pca", matching run_pca()'s own default output key.
    adjusted_basis: where to store the batch-corrected embedding
        (adata.obsm key) -- defaults to "X_pca_harmony".

    IMPORTANT implementation note: this calls harmonypy.run_harmony()
    DIRECTLY, rather than going through scanpy.external.pp.harmony_integrate().
    Confirmed via direct testing that scanpy 1.12.3's own wrapper
    unconditionally transposes harmonypy's output (`harmony_out.Z_corr.T`),
    but harmonypy 2.0.0's actual Z_corr is ALREADY shaped (n_cells, n_pcs)
    -- the wrapper's transpose is apparently left over from an older
    harmonypy version with a different output orientation, and produces a
    hard AnnData shape-validation crash with the versions actually pinned
    in this project's environment.yml. Calling harmonypy directly avoids
    depending on that broken intermediate wrapper entirely.

    Modifies adata in place (adds adata.obsm[adjusted_basis]). Returns
    adata for convenient chaining.
    """
    if method == "none":
        return adata
    if method != "harmony":
        raise ValueError(f"Unknown batch correction method: {method!r}")

    if basis not in adata.obsm:
        raise ValueError(
            f"Batch correction requires PCA to be run first (run_pca()) -- "
            f"adata.obsm has no '{basis}' embedding yet."
        )
    if not harmony_available():
        raise RuntimeError(
            "Harmony batch correction requires the `harmonypy` package, which is not installed."
        )

    import harmonypy

    ho = harmonypy.run_harmony(
        adata.obsm[basis], adata.obs, [batch_key], verbose=False, **harmony_kwargs
    )
    z_corr = ho.Z_corr
    # Defensive shape check, regardless of which orientation THIS
    # installed harmonypy version happens to return -- rather than
    # assuming (n_cells, n_pcs) unconditionally (confirmed correct for
    # harmonypy 2.0.0 via direct testing, but pinned unpinned in this
    # project's environment.yml, so a future harmonypy version could in
    # principle differ again, the same way it apparently already has
    # once before -- see this function's own docstring).
    if z_corr.shape[0] != adata.n_obs and z_corr.shape[1] == adata.n_obs:
        z_corr = z_corr.T
    if z_corr.shape[0] != adata.n_obs:
        raise RuntimeError(
            f"Harmony's output shape {ho.Z_corr.shape} could not be reconciled with "
            f"this dataset's {adata.n_obs} cells -- this may indicate an incompatible "
            f"harmonypy version."
        )
    adata.obsm[adjusted_basis] = np.asarray(z_corr)

    return adata


# ---------------------------------------------------------------------------
# 3.5: Clustering
# ---------------------------------------------------------------------------

CLUSTERING_METHOD_OPTIONS = {
    "leiden": {
        "label": "Leiden (recommended default)",
        "explanation": (
            "Graph-based community detection, current best practice for single-cell "
            "clustering -- explicitly guarantees well-connected clusters (an improvement "
            "Leiden makes over Louvain's own known 'badly connected communities' failure "
            "mode). A rigorous 2025 benchmark spanning 34.9 million parameter combinations "
            "found no significant performance difference between Leiden and Louvain on "
            "average -- so both are offered here as genuine, real choices, not "
            "'recommended vs. deprecated.'"
        ),
    },
    "louvain": {
        "label": "Louvain",
        "explanation": (
            "The original graph-based community detection algorithm Leiden was later built "
            "to improve upon. Still a reasonable, well-established choice per the same 2025 "
            "benchmark referenced above. Runs here via the same igraph backend Leiden uses "
            "(no separate 'louvain' package needed) -- note this specific backend "
            "('igraph' flavor) does NOT use the resolution parameter the way Leiden does; "
            "see resolution's own parameter note below."
        ),
    },
}
DEFAULT_CLUSTERING_METHOD = "leiden"
DEFAULT_CLUSTERING_RESOLUTION = 1.0
DEFAULT_N_NEIGHBORS = 15


def compute_neighbors(adata, use_rep=None, n_neighbors=DEFAULT_N_NEIGHBORS, random_state=0):
    """
    Compute the nearest-neighbor graph clustering AND embedding
    (UMAP/t-SNE) both depend on -- must be run once, after PCA (and
    batch correction, if used), before either clustering or embedding.

    use_rep: which .obsm embedding to compute neighbors on -- pass
        "X_pca_harmony" if batch correction was applied, or "X_pca" if
        not (or if method="none" was used in run_batch_correction()).
        If None, scanpy's own default behavior applies.

    Modifies adata in place (adds adata.obsp["distances"],
    adata.obsp["connectivities"], adata.uns["neighbors"]). Returns
    adata for convenient chaining.
    """
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, random_state=random_state)
    return adata


def run_clustering(adata, method=DEFAULT_CLUSTERING_METHOD, resolution=DEFAULT_CLUSTERING_RESOLUTION,
                    key_added=None, random_state=0):
    """
    Cluster cells using the chosen graph-based method -- requires
    compute_neighbors() to have been run first.

    method: "leiden" or "louvain" -- see CLUSTERING_METHOD_OPTIONS.
    resolution: higher values produce MORE, SMALLER clusters. Confirmed
        via direct testing: this parameter has a REAL effect for Leiden,
        but NO EFFECT AT ALL for Louvain when run through the "igraph"
        flavor backend used here (scanpy itself prints a warning to this
        effect) -- this is called out explicitly in this function's
        return value (see below) so the UI layer can surface it rather
        than silently accept a resolution value that does nothing.
    key_added: the .obs column name to store cluster labels under --
        defaults to the method name itself ("leiden" or "louvain") if
        not given, matching scanpy's own default convention.

    Returns (adata, resolution_had_effect: bool) -- resolution_had_effect
    is False specifically for method="louvain" (since this runs via the
    "igraph" flavor, which ignores resolution entirely), True for
    method="leiden".
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "Clustering requires the neighbor graph to be computed first (compute_neighbors())."
        )

    resolution_had_effect = True
    if method == "leiden":
        key_added = key_added or "leiden"
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added, flavor="igraph",
                     n_iterations=2, random_state=random_state)
    elif method == "louvain":
        key_added = key_added or "louvain"
        sc.tl.louvain(adata, resolution=resolution, key_added=key_added, flavor="igraph",
                      random_state=random_state)
        resolution_had_effect = False
    else:
        raise ValueError(f"Unknown clustering method: {method!r}")

    return adata, resolution_had_effect


def get_cluster_summary(adata, cluster_key):
    """
    Return a per-cluster cell count summary, as a plain DataFrame --
    for the "cluster size" bar chart in the UI layer.

    Returns None if cluster_key isn't in adata.obs.
    """
    if cluster_key not in adata.obs.columns:
        return None
    counts = adata.obs[cluster_key].value_counts().sort_index()
    return pd.DataFrame({"cluster": counts.index.astype(str), "n_cells": counts.values})


# ---------------------------------------------------------------------------
# 3.6: Embeddings (UMAP + t-SNE, 2D and 3D)
# ---------------------------------------------------------------------------

EMBEDDING_METHOD_OPTIONS = {
    "umap": {
        "label": "UMAP (recommended default)",
        "explanation": (
            "The standard visualization embedding for single-cell data -- built directly on "
            "the same neighbor graph used for clustering, so clusters that look separated in "
            "UMAP space genuinely reflect the graph structure clustering itself used, unlike "
            "t-SNE (see its own note below)."
        ),
    },
    "tsne": {
        "label": "t-SNE",
        "explanation": (
            "An older, still-widely-used alternative embedding. Tends to more aggressively "
            "separate distinct clusters visually, but distances BETWEEN separated clusters "
            "(and relative cluster sizes) are not as reliably meaningful as UMAP's -- "
            "generally treat t-SNE as a visualization aid rather than a source of "
            "quantitative distance/size claims."
        ),
    },
}
DEFAULT_EMBEDDING_METHOD = "umap"


def run_embedding(adata, method=DEFAULT_EMBEDDING_METHOD, use_rep=None, n_components=2,
                   random_state=0):
    """
    Compute a 2D or 3D visualization embedding.

    method: "umap" or "tsne" -- see EMBEDDING_METHOD_OPTIONS.
    use_rep: which .obsm embedding to run this on -- pass
        "X_pca_harmony" if batch correction was applied, else "X_pca".
    n_components: 2 (standard 2D plot) or 3 (3D plot, per explicit user
        request -- rendered via Plotly's Scatter3d in the UI layer).

    Modifies adata in place (adds adata.obsm["X_umap"] or
    adata.obsm["X_tsne"], both shaped (n_cells, n_components)). Returns
    adata for convenient chaining.
    """
    if method == "umap":
        if "neighbors" not in adata.uns:
            raise ValueError(
                "UMAP requires the neighbor graph to be computed first (compute_neighbors()) "
                "-- UMAP is built directly on the same neighbor graph used for clustering."
            )
        sc.tl.umap(adata, n_components=n_components, random_state=random_state)
    elif method == "tsne":
        if use_rep is None:
            raise ValueError("t-SNE requires an explicit use_rep (e.g. 'X_pca' or 'X_pca_harmony').")
        sc.tl.tsne(adata, use_rep=use_rep, n_pcs=None, random_state=random_state)
        # scanpy's tl.tsne does not natively support n_components=3 --
        # confirmed via its own signature (always produces a 2D
        # embedding). For 3D t-SNE, call sklearn.manifold.TSNE directly.
        if n_components == 3:
            from sklearn.manifold import TSNE
            rep_matrix = adata.obsm[use_rep]
            tsne_3d = TSNE(n_components=3, random_state=random_state).fit_transform(rep_matrix)
            adata.obsm["X_tsne"] = tsne_3d
    else:
        raise ValueError(f"Unknown embedding method: {method!r}")

    return adata


def get_embedding_coordinates(adata, method=DEFAULT_EMBEDDING_METHOD, color_by_columns=None):
    """
    Extract embedding coordinates (2D or 3D) as a plain DataFrame, one
    row per cell.

    method: "umap" or "tsne" -- determines which .obsm key is read.
    color_by_columns: list of .obs column names to include -- if None,
        includes every column currently in .obs.

    Returns a DataFrame with columns "Dim1", "Dim2", and (if the stored
    embedding has 3 components) "Dim3", plus "cell_barcode" and every
    requested .obs column. Returns None if this embedding hasn't been
    computed yet.
    """
    key = "X_umap" if method == "umap" else "X_tsne"
    if key not in adata.obsm:
        return None

    coords = adata.obsm[key]
    n_dims = coords.shape[1]
    dim_cols = [f"Dim{i+1}" for i in range(n_dims)]

    df = pd.DataFrame(coords, columns=dim_cols, index=adata.obs_names)
    df["cell_barcode"] = adata.obs_names

    cols_to_add = color_by_columns if color_by_columns is not None else list(adata.obs.columns)
    for col in cols_to_add:
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values

    return df


# ---------------------------------------------------------------------------
# 3.7: Cell-type annotation
# ---------------------------------------------------------------------------
#
# Three genuinely different annotation strategies exist in the field
# (confirmed via direct research, not assumed): manual marker-gene
# scoring, reference-based mapping (SingleR/Azimuth), and classifier-
# based prediction (CellTypist/scANVI). Multiple current sources
# explicitly recommend combining 2-3 methods and trusting the
# consensus, not picking just one.
#
# This project's environment.yml includes CellTypist (classifier-based)
# but NOT SingleR/Azimuth (both R/Bioconductor-based reference-mapping
# tools) -- so only manual marker scoring (the required baseline, no
# extra dependency) and CellTypist (optional automated first pass) are
# built here. SingleR/Azimuth support could be added later as a
# reference-based third option if that package is added to
# environment.yml -- flagged here rather than silently omitted.

ANNOTATION_METHOD_OPTIONS = {
    "manual_markers": {
        "label": "Manual marker-gene scoring (required baseline)",
        "explanation": (
            "You provide a marker gene panel per expected cell type (e.g. CD3D/CD3E for "
            "T cells, MS4A1/CD79A for B cells) -- each cell gets a per-cell score for each "
            "panel, which you can then use to decide each cluster's identity. Transparent, "
            "requires no external model/download, and forces you to actually look at your "
            "data -- the field's own literature considers this necessary even alongside "
            "automated tools, not a step you can skip."
        ),
    },
    "celltypist": {
        "label": "CellTypist (optional automated first-pass)",
        "explanation": (
            "A classifier trained on large reference atlases, predicting a cell type label "
            "for every cell directly. Fast and convenient as a STARTING GUESS to review and "
            "correct against your own marker-based judgment -- not a substitute for it. "
            "Requires downloading a pretrained model (one-time, per model) and remapping "
            "your data to gene SYMBOLS (CellTypist's own requirement, different from this "
            "pipeline's default gene-ID-indexed data)."
        ),
    },
}
DEFAULT_ANNOTATION_METHOD = "manual_markers"


def score_marker_gene_sets(adata, marker_gene_sets, layer="lognorm"):
    """
    Compute a per-cell score for each of one or more marker gene panels
    (e.g. {"T_cell": ["CD3D", "CD3E"], "B_cell": ["MS4A1", "CD79A"]}) --
    the manual marker-scoring baseline (see ANNOTATION_METHOD_OPTIONS).

    marker_gene_sets: dict {cell_type_label: [gene_symbol_or_id, ...]}.
        Gene identifiers are matched against BOTH adata.var_names
        (gene IDs, this pipeline's default index) AND
        adata.var["gene_symbol"] (if present) -- so a user can supply
        either familiar gene SYMBOLS (e.g. "CD3D") or raw gene IDs.
    layer: which layer to score from -- defaults to "lognorm".

    Adds one new .obs column per cell type, named f"score_{cell_type_label}"
    (via scanpy's own sc.tl.score_genes()). Silently skips (with a note
    in the returned dict) any marker set where NONE of its genes were
    found in this dataset at all.

    Returns a dict {cell_type_label: {"n_genes_found": int,
    "n_genes_requested": int, "genes_found": [str, ...]}}.
    """
    if layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found -- run normalize() first (creates 'lognorm' or "
            f"'pearson_residuals', depending on the method used)."
        )

    symbol_to_id = {}
    if "gene_symbol" in adata.var.columns:
        symbol_to_id = dict(zip(adata.var["gene_symbol"], adata.var_names))

    results = {}
    original_X = adata.X
    try:
        adata.X = adata.layers[layer]
        for cell_type_label, requested_genes in marker_gene_sets.items():
            resolved_ids = []
            for gene in requested_genes:
                if gene in adata.var_names:
                    resolved_ids.append(gene)
                elif gene in symbol_to_id:
                    resolved_ids.append(symbol_to_id[gene])

            resolved_ids = sorted(set(resolved_ids))
            results[cell_type_label] = {
                "n_genes_found": len(resolved_ids),
                "n_genes_requested": len(requested_genes),
                "genes_found": resolved_ids,
            }
            if not resolved_ids:
                continue

            score_name = f"score_{cell_type_label}"
            sc.tl.score_genes(adata, gene_list=resolved_ids, score_name=score_name)
    finally:
        adata.X = original_X

    return results


def find_cluster_markers(adata, groupby, method="wilcoxon", n_genes=25, layer="lognorm"):
    """
    Identify marker genes that distinguish each cluster from all other
    cells -- standard cluster-annotation workflow.

    IMPORTANT distinction from Phase 3.8's pseudobulk DESeq2 bridge:
    this function answers "which genes distinguish THIS CLUSTER from
    other cells, within this same combined dataset" -- NOT "which genes
    differ between conditions/samples" (a cross-sample comparison, for
    which treating individual cells as independent replicates produces
    severely inflated false-discovery rates -- Squair et al. 2021).

    Returns a DataFrame with columns: cluster, gene_id, gene_symbol (if
    available), score, pvalue, pvalue_adj, log2FoldChange.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs -- run clustering first.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found -- run normalize() first.")

    original_X = adata.X
    try:
        adata.X = adata.layers[layer]
        sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)
    finally:
        adata.X = original_X

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    symbol_map = {}
    if "gene_symbol" in adata.var.columns:
        symbol_map = adata.var["gene_symbol"].to_dict()

    rows = []
    for group in groups:
        for i in range(len(result["names"][group])):
            gene_id = result["names"][group][i]
            rows.append({
                "cluster": group,
                "gene_id": gene_id,
                "gene_symbol": symbol_map.get(gene_id, gene_id),
                "score": result["scores"][group][i],
                "pvalue": result["pvals"][group][i],
                "pvalue_adj": result["pvals_adj"][group][i],
                "log2FoldChange": result["logfoldchanges"][group][i],
            })
    return pd.DataFrame(rows)


# --- CellTypist (optional automated first-pass annotation) ---

def celltypist_available():
    "Check whether the celltypist package is installed."
    try:
        import celltypist  # noqa: F401
        return True
    except ImportError:
        return False


def get_celltypist_available_models():
    """
    List CellTypist model names available for download.

    Returns a dict {model_filename: description_string}, or an empty
    dict if celltypist isn't installed or its model index can't be
    read for any reason.
    """
    if not celltypist_available():
        return {}
    import celltypist
    try:
        descriptions = celltypist.models.models_description()
        if hasattr(descriptions, "to_dict"):
            return dict(zip(descriptions.iloc[:, 0], descriptions.iloc[:, 1]))
        return dict(descriptions)
    except Exception:
        return {}


def celltypist_model_is_downloaded(model_name):
    "Check whether a specific CellTypist model has already been downloaded locally."
    if not celltypist_available():
        return False
    import celltypist
    try:
        celltypist.models.Model.load(model_name)
        return True
    except Exception:
        return False


def download_celltypist_model(model_name):
    """
    Download a specific CellTypist model (requires network access).

    Returns (success: bool, message: str).
    """
    if not celltypist_available():
        return False, "The `celltypist` package is not installed."
    import celltypist
    try:
        celltypist.models.download_models(model=[model_name])
        return True, f"Successfully downloaded CellTypist model '{model_name}'."
    except Exception as e:
        return False, f"Failed to download CellTypist model '{model_name}': {e}"


def run_celltypist_annotation(adata, model_name="Immune_All_Low.pkl", majority_voting=True,
                               cluster_key=None, layer="lognorm"):
    """
    Run CellTypist's automated cell-type classifier on this dataset --
    see ANNOTATION_METHOD_OPTIONS for why this is offered as an
    OPTIONAL first-pass guess, not a replacement for manual marker
    review.

    IMPORTANT: CellTypist requires gene SYMBOLS as var_names -- this
    pipeline indexes adata.var_names by gene ID by default, with
    symbols stored separately in adata.var["gene_symbol"]. This
    function creates a TEMPORARY COPY of adata with var_names remapped
    to symbols purely for CellTypist's own internal use -- the
    ORIGINAL adata passed in is never mutated by this remapping; only
    the resulting predicted-label columns are written back onto it.

    Adds two new .obs columns to the ORIGINAL adata (not a copy):
        "celltypist_predicted_label"
        "celltypist_majority_voting" (only if majority_voting=True)

    Returns (success: bool, message: str).
    """
    if not celltypist_available():
        return False, "The `celltypist` package is not installed."
    if layer not in adata.layers:
        return False, f"Layer '{layer}' not found -- run normalize() first."
    if majority_voting and (cluster_key is None or cluster_key not in adata.obs.columns):
        return False, (
            "majority_voting=True requires cluster_key to be set to an existing .obs "
            "column with this project's own cluster labels (e.g. 'leiden') -- run "
            "clustering first."
        )
    if not celltypist_model_is_downloaded(model_name):
        return False, (
            f"CellTypist model '{model_name}' has not been downloaded yet -- call "
            f"download_celltypist_model('{model_name}') first."
        )

    import celltypist

    if "gene_symbol" not in adata.var.columns:
        return False, (
            "This dataset has no 'gene_symbol' column in .var -- CellTypist requires "
            "gene symbols and cannot run on gene IDs alone."
        )

    symbols = adata.var["gene_symbol"].astype(str)
    valid_mask = symbols.notna() & (symbols != "") & (symbols != "nan")
    temp = adata[:, valid_mask.values].copy()
    temp.var_names = symbols[valid_mask].values
    temp.var_names_make_unique()
    temp.X = temp.layers[layer]

    try:
        prediction = celltypist.annotate(
            temp, model=model_name,
            majority_voting=majority_voting,
            over_clustering=(adata.obs[cluster_key].astype(str).values if majority_voting else None),
        )
    except Exception as e:
        return False, f"CellTypist annotation failed: {e}"

    result_df = prediction.predicted_labels
    adata.obs["celltypist_predicted_label"] = result_df["predicted_labels"].reindex(adata.obs_names).values
    if majority_voting and "majority_voting" in result_df.columns:
        adata.obs["celltypist_majority_voting"] = result_df["majority_voting"].reindex(adata.obs_names).values

    return True, f"CellTypist annotation complete using model '{model_name}'."


# --- Visualization data extraction (dot plot, violin, feature plot, heatmap) ---

def get_dotplot_data(adata, marker_genes, groupby, layer="lognorm"):
    """
    Compute the data needed for a standard single-cell dot plot.

    Returns a DataFrame with columns: gene_symbol, group, pct_expressing,
    mean_expression -- one row per (gene, group) combination.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found -- run normalize() first.")

    symbol_to_id = {}
    if "gene_symbol" in adata.var.columns:
        symbol_to_id = dict(zip(adata.var["gene_symbol"], adata.var_names))
    symbol_map = adata.var["gene_symbol"].to_dict() if "gene_symbol" in adata.var.columns else {}

    resolved_ids = []
    for gene in marker_genes:
        if gene in adata.var_names:
            resolved_ids.append(gene)
        elif gene in symbol_to_id:
            resolved_ids.append(symbol_to_id[gene])
    resolved_ids = list(dict.fromkeys(resolved_ids))

    if not resolved_ids:
        return pd.DataFrame(columns=["gene_symbol", "group", "pct_expressing", "mean_expression"])

    expr_matrix = adata[:, resolved_ids].layers[layer]
    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    rows = []
    groups = adata.obs[groupby].astype(str)
    for group_val in sorted(groups.unique()):
        mask = (groups == group_val).values
        sub = expr_matrix[mask, :]
        for j, gene_id in enumerate(resolved_ids):
            gene_col = sub[:, j]
            pct_expr = float((gene_col > 0).mean() * 100)
            mean_expr_among_expressing = float(gene_col[gene_col > 0].mean()) if (gene_col > 0).any() else 0.0
            rows.append({
                "gene_symbol": symbol_map.get(gene_id, gene_id),
                "group": group_val,
                "pct_expressing": pct_expr,
                "mean_expression": mean_expr_among_expressing,
            })
    return pd.DataFrame(rows)


def get_violin_plot_data(adata, gene, groupby, layer="lognorm"):
    """
    Extract per-cell expression values for ONE gene, alongside its
    group membership -- for a violin plot.

    Returns a DataFrame with columns: expression, group -- one row per
    cell. Returns None if the gene couldn't be resolved.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found -- run normalize() first.")

    gene_id = gene
    if gene not in adata.var_names and "gene_symbol" in adata.var.columns:
        symbol_to_id = dict(zip(adata.var["gene_symbol"], adata.var_names))
        gene_id = symbol_to_id.get(gene)

    if gene_id is None or gene_id not in adata.var_names:
        return None

    expr = adata[:, gene_id].layers[layer]
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = np.asarray(expr).flatten()

    return pd.DataFrame({
        "expression": expr,
        "group": adata.obs[groupby].astype(str).values,
    })


def get_feature_plot_data(adata, gene, embedding_method=DEFAULT_EMBEDDING_METHOD, layer="lognorm"):
    """
    Extract embedding coordinates (UMAP or t-SNE) PLUS one gene's
    per-cell expression -- for a "feature plot".

    Returns a DataFrame with columns: Dim1, Dim2, (Dim3 if 3D),
    expression, cell_barcode. Returns None if the gene couldn't be
    resolved, or if the requested embedding hasn't been computed yet.
    """
    embed_df = get_embedding_coordinates(adata, method=embedding_method, color_by_columns=[])
    if embed_df is None:
        return None

    gene_id = gene
    if gene not in adata.var_names and "gene_symbol" in adata.var.columns:
        symbol_to_id = dict(zip(adata.var["gene_symbol"], adata.var_names))
        gene_id = symbol_to_id.get(gene)

    if gene_id is None or gene_id not in adata.var_names:
        return None
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found -- run normalize() first.")

    expr = adata[:, gene_id].layers[layer]
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    embed_df["expression"] = np.asarray(expr).flatten()
    return embed_df


def get_marker_heatmap_data(adata, marker_genes, groupby, layer="lognorm", n_cells_per_group=50):
    """
    Build a gene x cell z-scored expression matrix for a set of marker
    genes, grouped and (optionally) downsampled by cluster/group.

    Returns (z_df, ordered_cell_groups). Returns (None, None) if no
    requested genes could be resolved.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found -- run normalize() first.")

    symbol_to_id = {}
    symbol_map = {}
    if "gene_symbol" in adata.var.columns:
        symbol_to_id = dict(zip(adata.var["gene_symbol"], adata.var_names))
        symbol_map = adata.var["gene_symbol"].to_dict()

    resolved_ids = []
    for gene in marker_genes:
        if gene in adata.var_names:
            resolved_ids.append(gene)
        elif gene in symbol_to_id:
            resolved_ids.append(symbol_to_id[gene])
    resolved_ids = list(dict.fromkeys(resolved_ids))

    if not resolved_ids:
        return None, None

    groups = adata.obs[groupby].astype(str)
    selected_cell_indices = []
    rng = np.random.default_rng(0)
    for group_val in sorted(groups.unique()):
        group_indices = np.where((groups == group_val).values)[0]
        if n_cells_per_group is not None and len(group_indices) > n_cells_per_group:
            group_indices = rng.choice(group_indices, size=n_cells_per_group, replace=False)
        selected_cell_indices.extend(sorted(group_indices))

    subset = adata[selected_cell_indices, resolved_ids]
    expr_matrix = subset.layers[layer]
    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    gene_labels = [symbol_map.get(g, g) for g in resolved_ids]
    cell_barcodes = subset.obs_names

    df = pd.DataFrame(expr_matrix.T, index=gene_labels, columns=cell_barcodes)
    row_mean = df.mean(axis=1)
    row_std = df.std(axis=1).replace(0, 1)
    z_df = df.sub(row_mean, axis=0).div(row_std, axis=0)

    ordered_groups = groups.loc[cell_barcodes]
    return z_df, ordered_groups

# ---------------------------------------------------------------------------
# 3.8: Pseudobulk aggregation
# ---------------------------------------------------------------------------

DEFAULT_MIN_CELLS_PER_PSEUDOBULK_GROUP = 10


def get_pseudobulk_group_sizes(adata, groupby_columns):
    """
    Preview how many cells would contribute to each pseudobulk group.

    Returns a DataFrame with one row per group: the groupby_columns'
    own values, plus "n_cells" -- sorted by n_cells ascending.

    Raises ValueError if any requested column isn't in adata.obs.
    """
    missing = [c for c in groupby_columns if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"Column(s) not found in adata.obs: {missing}")

    sizes = adata.obs.groupby(groupby_columns, dropna=False, observed=True).size()
    sizes_df = sizes.reset_index(name="n_cells")
    return sizes_df.sort_values("n_cells", ascending=True).reset_index(drop=True)


def _build_pseudobulk_sample_name(group_key, groupby_columns):
    """
    Build a single, filesystem/CSV-column-safe pseudobulk sample name
    from a groupby group key.
    """
    if not isinstance(group_key, tuple):
        group_key = (group_key,)
    parts = [str(v).replace(" ", "-") for v in group_key]
    return "_".join(parts)


def aggregate_pseudobulk(adata, groupby_columns, layer="counts",
                          min_cells=DEFAULT_MIN_CELLS_PER_PSEUDOBULK_GROUP):
    """
    Aggregate raw per-cell counts into pseudobulk "samples" by summing
    within each group defined by groupby_columns.

    Returns (pseudobulk_counts_df, pseudobulk_metadata_df, excluded_groups).

    Raises ValueError if "counts" (or the specified layer) is missing
    from adata.layers, or if groupby_columns references a column not in
    adata.obs.
    """
    if layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found in adata.layers -- pseudobulk aggregation requires raw "
            f"counts (set automatically by load_and_combine_samples() as the 'counts' layer)."
        )
    missing_cols = [c for c in groupby_columns if c not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Column(s) not found in adata.obs: {missing_cols}")

    counts_matrix = adata.layers[layer]
    if hasattr(counts_matrix, "toarray"):
        counts_matrix = counts_matrix.toarray()

    groups = adata.obs.groupby(groupby_columns, dropna=False, observed=True).indices

    sample_columns = {}
    metadata_rows = []
    excluded_groups = []

    for group_key, cell_indices in groups.items():
        n_cells = len(cell_indices)
        group_key_tuple = group_key if isinstance(group_key, tuple) else (group_key,)

        if n_cells < min_cells:
            excluded_entry = dict(zip(groupby_columns, group_key_tuple))
            excluded_entry["n_cells"] = n_cells
            excluded_groups.append(excluded_entry)
            continue

        pseudobulk_sample_name = _build_pseudobulk_sample_name(group_key_tuple, groupby_columns)
        summed_counts = counts_matrix[cell_indices, :].sum(axis=0)
        summed_counts = _flatten_to_1d(summed_counts)
        sample_columns[pseudobulk_sample_name] = summed_counts

        metadata_row = dict(zip(groupby_columns, group_key_tuple))
        metadata_row["sample"] = pseudobulk_sample_name
        metadata_row["n_cells_aggregated"] = n_cells
        metadata_rows.append(metadata_row)

    if not sample_columns:
        raise ValueError(
            f"No groups met the min_cells={min_cells} threshold -- every group had fewer "
            f"contributing cells than this. Lower min_cells, or check your groupby_columns "
            f"choice with get_pseudobulk_group_sizes() first."
        )

    pseudobulk_counts_df = pd.DataFrame(sample_columns, index=adata.var_names)
    pseudobulk_counts_df.index.name = "gene_id"
    pseudobulk_counts_df = pseudobulk_counts_df.reset_index()

    pseudobulk_metadata_df = pd.DataFrame(metadata_rows)
    other_cols = [c for c in groupby_columns if c != "sample"]
    ordered_cols = ["sample"] + other_cols + ["n_cells_aggregated"]
    pseudobulk_metadata_df = pseudobulk_metadata_df[ordered_cols]

    return pseudobulk_counts_df, pseudobulk_metadata_df, excluded_groups


def _flatten_to_1d(summed_counts):
    "Normalize a numpy matrix/array sum() result into a plain 1D array."
    import numpy as np
    arr = np.asarray(summed_counts)
    return arr.flatten()


def save_pseudobulk_for_deseq2(pseudobulk_counts_df, pseudobulk_metadata_df, dest_dir):
    """
    Write a pseudobulk counts/metadata pair to disk in the exact file
    shapes deseq2_manager.run_deseq2_analysis() expects.

    Returns (counts_path, metadata_path).
    """
    os.makedirs(dest_dir, exist_ok=True)
    counts_path = os.path.join(dest_dir, "pseudobulk_counts.csv")
    metadata_path = os.path.join(dest_dir, "pseudobulk_metadata.csv")
    pseudobulk_counts_df.to_csv(counts_path, index=False)
    pseudobulk_metadata_df.to_csv(metadata_path, index=False)
    return counts_path, metadata_path

# ---------------------------------------------------------------------------
# 3.9: Compositional analysis -- method options + tradeoffs
# ---------------------------------------------------------------------------

COMPOSITIONAL_METHOD_OPTIONS = {
    "propeller": {
        "label": "propeller (recommended default -- requires R/speckle)",
        "explanation": (
            "Transforms each sample's cell-type proportions (logit transform by default) "
            "and fits a proper limma model treating each SAMPLE as one observation -- the "
            "statistically correct unit of replication, avoiding the inflated false-positive "
            "rates that come from treating individual cells as independent (Phipson et al. "
            "2022, Bioinformatics). Uses empirical Bayes variance shrinkage (borrowing "
            "statistical strength across cell types), which is particularly valuable with "
            "few samples or few cell types. Requires the R `speckle` Bioconductor package."
        ),
    },
    "simple_test": {
        "label": "Simple proportion test (no R required, Python-only fallback)",
        "explanation": (
            "The same arcsine-square-root variance-stabilizing transform propeller itself "
            "uses, followed by a plain t-test (2 groups) or one-way ANOVA (3+ groups) per "
            "cell type, with Benjamini-Hochberg FDR correction across cell types. Correctly "
            "avoids pseudoreplication (still one observation per sample, not per cell), but "
            "lacks propeller's empirical Bayes variance-shrinkage step -- results can be "
            "noisier/less stable, especially with few samples or few cell types. "
            "**Also requires stricter replication than propeller**: EVERY group needs at "
            "least 2 samples -- propeller can still produce a result in that specific case "
            "by borrowing variance information across cell types, this fallback cannot. "
            "Use this when R/speckle isn't available in this environment; otherwise "
            "propeller is recommended."
        ),
    },
}
DEFAULT_COMPOSITIONAL_METHOD = "propeller"
DEFAULT_PROPELLER_TRANSFORM = "logit"
DEFAULT_MIN_SAMPLES_PER_GROUP_FOR_REPLICATION = 2


def compositional_tools_available():
    """
    Check whether Rscript is available.
    """
    import shutil
    return shutil.which("Rscript") is not None


def check_compositional_replication(sample_group_df, group_column, require_all_groups_replicated=False):
    """
    Check whether there is enough sample-level replication to run a
    compositional comparison.

    Returns a dict:
        {
            "is_valid": bool,
            "n_groups": int,
            "group_counts": {group_label: n_samples, ...},
            "under_replicated_groups": [group_label, ...],
            "message": str,
        }
    """
    if group_column not in sample_group_df.columns:
        return {
            "is_valid": False, "n_groups": 0, "group_counts": {},
            "under_replicated_groups": [],
            "message": f"Column '{group_column}' not found in the provided sample-level data.",
        }

    group_counts_series = sample_group_df[group_column].value_counts()
    group_counts = {str(k): int(v) for k, v in group_counts_series.items()}
    n_groups = len(group_counts)
    under_replicated = [g for g, n in group_counts.items() if n < DEFAULT_MIN_SAMPLES_PER_GROUP_FOR_REPLICATION]

    if n_groups < 2:
        return {
            "is_valid": False, "n_groups": n_groups, "group_counts": group_counts,
            "under_replicated_groups": under_replicated,
            "message": "At least 2 groups are needed for a compositional comparison -- only "
                       f"{n_groups} group(s) found.",
        }

    if require_all_groups_replicated:
        is_valid = len(under_replicated) == 0
        floor_description = "EVERY group needs"
    else:
        is_valid = any(n >= DEFAULT_MIN_SAMPLES_PER_GROUP_FOR_REPLICATION for n in group_counts.values())
        floor_description = "at least one group needs"

    if is_valid:
        message = (
            f"✅ {n_groups} groups found, with sample counts: "
            f"{', '.join(f'{g}={n}' for g, n in group_counts.items())}. "
            "Sufficient biological replication to estimate variability."
        )
    else:
        message = (
            f"⚠️ Insufficient replication for this method: "
            f"{', '.join(f'{g}={n}' for g, n in group_counts.items())}. "
            f"This method requires {floor_description} at least "
            f"{DEFAULT_MIN_SAMPLES_PER_GROUP_FOR_REPLICATION} samples "
            + (
                "(propeller's empirical Bayes variance-borrowing across cell types allows this "
                "more lenient floor)."
                if not require_all_groups_replicated else
                "(the simple, no-R-required fallback test has no variance-borrowing mechanism, so "
                "a group with only 1 sample makes the underlying t-test/ANOVA mathematically "
                "undefined -- propeller can handle this specific case, this fallback cannot)."
            )
        )

    return {
        "is_valid": is_valid, "n_groups": n_groups, "group_counts": group_counts,
        "under_replicated_groups": under_replicated, "message": message,
    }


def get_sample_level_group_mapping(adata, sample_key, group_column):
    """
    Build a one-row-per-SAMPLE DataFrame mapping each sample to its
    group/condition value, de-duplicated down from the cell-level
    adata.obs.

    Returns a DataFrame with columns [sample_key, group_column], one
    row per unique sample.

    Raises ValueError if any single sample has more than one distinct
    group_column value across its own cells.
    """
    if sample_key not in adata.obs.columns or group_column not in adata.obs.columns:
        raise ValueError(f"'{sample_key}' and/or '{group_column}' not found in adata.obs.")

    per_sample_groups = adata.obs.groupby(sample_key, observed=True)[group_column].nunique()
    inconsistent_samples = per_sample_groups[per_sample_groups > 1].index.tolist()
    if inconsistent_samples:
        raise ValueError(
            f"Sample(s) {inconsistent_samples} have more than one distinct value in "
            f"'{group_column}' across their own cells -- a sample must belong to exactly "
            f"one group/condition. Check this metadata column for errors."
        )

    mapping = adata.obs[[sample_key, group_column]].drop_duplicates().reset_index(drop=True)
    return mapping


def compute_cell_type_proportions(adata, sample_key, cluster_key):
    """
    Compute per-sample cell-type/cluster proportions -- one row per
    sample, one column per cluster, values summing to 1 across each
    row.

    Returns (proportions_df, counts_df).

    Raises ValueError if either column isn't present in adata.obs.
    """
    if sample_key not in adata.obs.columns or cluster_key not in adata.obs.columns:
        raise ValueError(f"'{sample_key}' and/or '{cluster_key}' not found in adata.obs.")

    counts_df = pd.crosstab(adata.obs[sample_key], adata.obs[cluster_key])
    proportions_df = counts_df.div(counts_df.sum(axis=1), axis=0)
    return proportions_df, counts_df


def get_composition_barplot_data(proportions_df):
    """
    Reshape a wide proportions_df into long format for a stacked
    composition bar plot.

    Returns a DataFrame with columns: sample, cluster, proportion.
    """
    long_df = proportions_df.reset_index().melt(
        id_vars=proportions_df.index.name or "index",
        var_name="cluster", value_name="proportion",
    )
    long_df = long_df.rename(columns={proportions_df.index.name or "index": "sample"})
    return long_df


# ---------------------------------------------------------------------------
# Method 1: propeller (R/speckle)
# ---------------------------------------------------------------------------

_PROPELLER_R_SCRIPT = r'''
suppressMessages({
  library(jsonlite)
  library(speckle)
  library(limma)
})
args <- commandArgs(trailingOnly = TRUE)
job_spec_path <- args[1]
job <- fromJSON(job_spec_path)

cell_df <- read.csv(job$cell_level_path, stringsAsFactors = FALSE)
sample_group_df <- read.csv(job$sample_group_path, stringsAsFactors = FALSE)

clusters <- as.factor(cell_df$cluster)
sample_ids <- as.factor(cell_df$sample)

prop.list <- getTransformedProps(clusters = clusters, sample = sample_ids, transform = job$transform)

sample_order <- colnames(prop.list$Proportions)
rownames(sample_group_df) <- sample_group_df$sample
sample_group_df <- sample_group_df[sample_order, , drop = FALSE]

group <- as.factor(sample_group_df$group)
group_levels <- levels(group)
n_groups <- length(group_levels)

design <- model.matrix(~0 + group)
colnames(design) <- group_levels

cat(paste("Groups:", paste(group_levels, collapse = ", ")), "\n")
cat(paste("Samples (in prop.list order):", paste(sample_order, collapse = ", ")), "\n")

if (n_groups == 2) {
  contrast_name <- paste0(group_levels[1], "-", group_levels[2])
  contrasts <- makeContrasts(contrasts = contrast_name, levels = design)
  result <- propeller.ttest(
    prop.list = prop.list, design = design, contrasts = contrasts,
    robust = job$robust, trend = job$trend, sort = TRUE
  )
  cat(paste("Ran propeller.ttest():", contrast_name), "\n")
} else {
  result <- propeller.anova(
    prop.list = prop.list, design = design, coef = seq_len(n_groups),
    robust = job$robust, trend = job$trend, sort = TRUE
  )
  cat(paste("Ran propeller.anova() across", n_groups, "groups"), "\n")
}

result_df <- as.data.frame(result)
result_df$cluster <- rownames(result_df)
result_df <- result_df[, c("cluster", setdiff(colnames(result_df), "cluster"))]
write.csv(result_df, file.path(job$output_dir, "propeller_results.csv"), row.names = FALSE)

props_out <- as.data.frame(t(prop.list$Proportions))
props_out$sample <- rownames(props_out)
props_out <- props_out[, c("sample", setdiff(colnames(props_out), "sample"))]
write.csv(props_out, file.path(job$output_dir, "propeller_proportions.csv"), row.names = FALSE)

cat("propeller analysis completed successfully.\n")
'''


def run_propeller_analysis(adata, sample_key, cluster_key, group_column,
                            output_dir, work_dir, transform=DEFAULT_PROPELLER_TRANSFORM,
                            robust=True, trend=False, timeout=1800):
    """
    Run propeller (via the R `speckle` package) on this combined
    AnnData -- the recommended default compositional-analysis method.

    Returns (success: bool, log: str).
    """
    import json
    import subprocess

    if not compositional_tools_available():
        return False, (
            "Rscript was not found on this system. R with the `speckle` package needs to be "
            "installed in your environment (bioconductor-speckle) before this step can run."
        )

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    sample_group_df = get_sample_level_group_mapping(adata, sample_key, group_column)
    sample_group_df = sample_group_df.rename(columns={sample_key: "sample", group_column: "group"})
    replication_check = check_compositional_replication(sample_group_df, "group")
    if not replication_check["is_valid"]:
        return False, replication_check["message"]

    cell_level_df = adata.obs[[sample_key, cluster_key]].rename(
        columns={sample_key: "sample", cluster_key: "cluster"}
    )
    cell_level_path = os.path.join(work_dir, "propeller_cell_level.csv")
    sample_group_path = os.path.join(work_dir, "propeller_sample_group.csv")
    cell_level_df.to_csv(cell_level_path, index=False)
    sample_group_df.to_csv(sample_group_path, index=False)

    job_spec = {
        "cell_level_path": os.path.abspath(cell_level_path),
        "sample_group_path": os.path.abspath(sample_group_path),
        "transform": transform,
        "robust": bool(robust),
        "trend": bool(trend),
        "output_dir": os.path.abspath(output_dir),
    }
    job_spec_path = os.path.join(work_dir, "propeller_job_spec.json")
    with open(job_spec_path, "w") as f:
        json.dump(job_spec, f, indent=2)

    r_script_path = os.path.join(work_dir, "run_propeller.R")
    with open(r_script_path, "w") as f:
        f.write(_PROPELLER_R_SCRIPT)

    cmd = ["Rscript", r_script_path, job_spec_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, f"propeller analysis failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, f"propeller analysis timed out after {timeout // 60} minutes."


def read_propeller_results(output_dir):
    "Read propeller's results CSV (one row per cluster/cell type) written by run_propeller_analysis(). Returns None if not found."
    path = os.path.join(output_dir, "propeller_results.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


def read_propeller_proportions(output_dir):
    "Read the per-sample proportions CSV propeller itself computed and exported. Returns None if not found."
    path = os.path.join(output_dir, "propeller_proportions.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Method 2: simple_test (Python-only fallback, no R required)
# ---------------------------------------------------------------------------

def run_simple_proportion_test(proportions_df, sample_group_df, group_column="group"):
    """
    Dependency-free Python fallback for compositional testing.

    Returns a DataFrame with columns: cluster, statistic, pvalue,
    padj, test_used ("t-test" for 2 groups, "ANOVA" for 3+) -- one row
    per cluster, sorted by padj ascending.

    Raises ValueError if replication is insufficient.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    sample_group_renamed = sample_group_df.rename(columns={group_column: "group"}) if group_column != "group" else sample_group_df.copy()
    replication_check = check_compositional_replication(sample_group_renamed, "group", require_all_groups_replicated=True)
    if not replication_check["is_valid"]:
        raise ValueError(replication_check["message"])

    sample_to_group = dict(zip(sample_group_renamed["sample"], sample_group_renamed["group"]))
    groups_per_sample = proportions_df.index.map(lambda s: sample_to_group.get(s))
    if groups_per_sample.isna().any():
        missing = proportions_df.index[groups_per_sample.isna()].tolist()
        raise ValueError(f"No group mapping found for sample(s): {missing}")

    transformed = np.arcsin(np.sqrt(proportions_df.clip(lower=0, upper=1)))

    unique_groups = sorted(set(groups_per_sample))
    n_groups = len(unique_groups)

    rows = []
    for cluster in transformed.columns:
        values_by_group = [
            transformed.loc[groups_per_sample == g, cluster].values
            for g in unique_groups
        ]
        if n_groups == 2:
            statistic, pvalue = stats.ttest_ind(values_by_group[0], values_by_group[1], equal_var=False)
            test_used = "t-test (Welch)"
        else:
            statistic, pvalue = stats.f_oneway(*values_by_group)
            test_used = "ANOVA"
        rows.append({"cluster": cluster, "statistic": float(statistic), "pvalue": float(pvalue), "test_used": test_used})

    result_df = pd.DataFrame(rows)
    valid_pvalues = result_df["pvalue"].notna()
    result_df["padj"] = np.nan
    if valid_pvalues.any():
        _, padj_values, _, _ = multipletests(result_df.loc[valid_pvalues, "pvalue"], method="fdr_bh")
        result_df.loc[valid_pvalues, "padj"] = padj_values

    return result_df.sort_values("padj", na_position="last").reset_index(drop=True)
