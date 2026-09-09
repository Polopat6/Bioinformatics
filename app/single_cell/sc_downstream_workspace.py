"""
single_cell/sc_downstream_workspace.py

Streamlit UI for Phase 3 of the Single-cell RNA-Seq pipeline --
multi-sample downstream analysis (normalization, HVG selection, PCA,
batch correction/integration, clustering, UMAP/t-SNE, cell-type
annotation, pseudobulk aggregation, and compositional analysis), built
on top of sc_downstream_manager.py's backend functions.

This file is being built incrementally, section by section:
    SECTION 1 (delivered): Step 1 (Combine Samples) -> Step 2
        (Normalization) -> Step 3 (HVG Selection) -> Step 4 (PCA),
        plus the shared caching/checkpoint/download infrastructure
        every later section also uses.
    SECTION 2 (delivered): Step 5 (Batch Correction) -> Step 6
        (Clustering) -> Step 7 (Embeddings).
    SECTION 3 (delivered): Step 8 (Cell-type Annotation).
    SECTION 4 (planned): Step 9 (Pseudobulk -> DESeq2 Bridge) -> Step
        10 (Compositional Analysis).

--- Why Phase 3 needs its own caching design, distinct from Phase 1/2
    (2026-08-21) ---
Phase 1/2's workspace pages track progress with a simple boolean per
step (scpm.mark_step_complete/has_completed_step), because each of
THEIR steps produces its own independent, separately-checkable output
file. Phase 3 is fundamentally different: scanpy/anndata's AnnData
object ACCUMULATES state across steps -- there is only ONE growing
object, not a series of independent files. Every Phase 3 step's actual
parameter dict (not just a boolean) is saved via
sc_project_manager.save_downstream_step_recipe(), and reuse of a
cached result ONLY happens when the user's CURRENT parameter
selections EXACTLY match what's saved (see
scpm.check_downstream_step_current()) -- a mismatch is always
surfaced explicitly to the user, never silently resolved either way.

--- Caching model: session_state (fast) + a single overwritten disk
    file (durable) + optional named checkpoints (deliberate) ---
_get_cached_adata()/_save_adata_state() implement a two-tier cache:
the LIVE AnnData object lives in st.session_state during an active
session, and is ALSO written to a single, default "current state"
.h5ad file (scpm.downstream_adata_path()) every time a step completes,
so a fresh session can resume exactly where a previous session left
off. _render_checkpoint_controls() offers an explicit, user-named,
independently-persisted alternative for preserving a SPECIFIC
parameter combination, opted into via its own UI action.

--- Widget defaults are read from the LAST SAVED recipe, not
    hardcoded (2026-08-21) ---
_recipe_default() is used for EVERY widget in this file that
represents a tracked step parameter, so reopening a project with no
new changes correctly shows "already up to date" instead of a false
mismatch caused purely by a widget resetting to its base default on a
cold reload.

--- Step 8 design notes (2026-08-23) ---
Cell-type annotation is NOT modeled as a single linear "run once, get
one cached result" step the way Steps 1-7 are -- the field's own
literature (and this project's own ANNOTATION_METHOD_OPTIONS
docstring in sc_downstream_manager.py) explicitly treats it as an
iterative, multi-method REVIEW process: explore cluster markers, score
manual marker panels, optionally cross-check against CellTypist, then
make a final human judgment call per cluster. Recipe-gating (the
"is_current / re-run" pattern from Steps 1-7) is therefore only
applied to the two genuinely re-run-able sub-steps (CellTypist, and
the final label assignment itself) -- cluster-marker exploration and
manual marker scoring are treated as repeatable, additive exploration
actions instead, matching how a real analyst actually uses them.

--- Steps 1-4 explanatory content added (2026-08-25) ---
Per direct user request, each of Steps 1-4 (and later, Steps 5-10)
includes plain-language "what this step does / why it matters / how to
interpret the output" guidance, aimed at someone who understands the
biology but may not have hands-on scanpy/anndata experience.
Explanatory text is placed directly beneath each step's own header
(the "why/what" framing, before any widgets), and, where a plot/result
already exists, an additional "How to read this" caption/expander is
added right next to that specific result.

--- PCA plot: independent per-axis PC selection + optional 3D Z-axis
    (2026-08-25) ---
_render_pca_step()'s own PCA scatter plot previously hardcoded PC1 (x)
vs. PC2 (y), 2D only. Per direct user request, this is now fully
configurable: three independent selectboxes let the user choose ANY
available PC for the X axis, Y axis, and (optionally, if a checkbox is
enabled) a Z axis. Enabling the Z-axis switches the plot from
plotly.express.scatter to plotly.express.scatter_3d.

--- Gene-symbol resolution wiring: Step 1 passes gtf_path through
    (2026-08-25) ---
sc_downstream_manager.py's load_and_combine_samples() accepts an
optional gtf_path parameter that overlays a corrected gene_symbol for
any gene whose STARsolo-provided symbol looks unresolved (see that
module's own docstring, "Gene-symbol resolution overlay"). This is a
complete no-op unless a caller actually supplies gtf_path -- so
_render_combine_step() below fetches this project's own confirmed
reference GTF path (via scpm.get_reference_choice(project)
["custom_gtf"]) and passes it straight through to
dsm.load_and_combine_samples(). IMPORTANT: this wiring must exist on
THIS side (the UI/workspace layer) for the manager-side fix to ever
actually take effect at all -- a manager-only fix with no
corresponding call-site wiring here would silently never activate,
which looks identical to "the fix doesn't work" from the outside.

--- Cluster/group label natural-sort ordering fix (2026-08-25) ---
A real, confirmed bug, reported directly via screenshots: numeric-
looking cluster labels ("0", "1", ..., "16") were displaying in a
confusing, non-numeric order in TWO places in this file specifically:
  1. Step 8a's "View top markers for cluster" dropdown used plain
     `sorted(markers_df["cluster"].unique(), key=str)` -- a
     LEXICOGRAPHIC string sort, which puts "10" immediately after "1"
     and before "2" (confirmed directly matching a real reported
     screenshot showing exactly this order).
  2. Step 7's UMAP/t-SNE embedding scatter plot's legend used
     `sorted(embed_df[color_by].astype(str).unique())` to build
     group_values for the color-picker UI, but NEVER passed an
     explicit category_orders argument to px.scatter/scatter_3d --
     meaning Plotly's own default legend ordering (order of first
     appearance in the underlying dataframe) took over instead,
     producing an even more scrambled, neither-numeric-nor-alphabetical
     legend order (confirmed directly matching a second real reported
     screenshot).
Both are fixed here by using sc_downstream_manager.py's new
natural_sort_unique() helper (see that module's own docstring) instead
of a bare sorted()/sorted(..., key=str) call, AND by explicitly passing
category_orders={color_by: group_values} into the Step 7 scatter calls
so Plotly is told the exact intended legend order directly, rather than
falling back to its own "first appearance" default.
"""
import os
import tempfile

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import sc_project_manager as scpm
import sc_downstream_manager as dsm
import singlecell_workspace as scw
import math
import sc_marker_panel_io as mpio
import sc_marker_confidence as mconf
import gene_id_mapper as gim

def _format_scientific_pvalue(value, sig_figs=3):
    """
    Format a p-value for display in true scientific notation, honestly
    distinguishing a real (if tiny) value from a genuine float64
    underflow (an actual 0.0 in memory that can no longer be recovered
    to its true magnitude) -- see this file's own module docstring,
    "Scientific-notation p-value display fix", for the full rationale.
    """
    if value is None:
        return ""
    try:
        if math.isnan(value):
            return ""
    except TypeError:
        return ""
    if value == 0.0:
        return "< 1e-308"
    return f"{value:.{sig_figs}e}"
WORKSPACE_KEY = "sc_downstream"
EXCLUDED_CELL_TYPES_UNS_KEY = "excluded_cell_types"

# ---------------------------------------------------------------------------
# AnnData session/disk cache helpers
# ---------------------------------------------------------------------------

def _session_key(project):
    return f"sc_downstream_adata_{project}"


def _get_cached_adata(project):
    """
    Return the currently active AnnData object for this project,
    preferring the in-memory session_state copy (fast, zero disk I/O)
    and falling back to the on-disk cache file
    (scpm.downstream_adata_path()) if this is a fresh session that
    hasn't loaded it into memory yet this run.

    Returns None if neither exists (a brand-new project with no Phase
    3 progress saved at all).
    """
    key = _session_key(project)
    if key in st.session_state:
        return st.session_state[key]
    if scpm.downstream_adata_cache_exists(project):
        import anndata as ad
        with st.spinner("Loading previously saved progress for this project..."):
            adata = ad.read_h5ad(scpm.downstream_adata_path(project))
        st.session_state[key] = adata
        return adata
    return None

def _strip_gene_colon_prefix(gene_id):
    """
    Strip a leading 'gene:' prefix (the GFF3 ID-attribute convention
    confirmed present in this project's real reference GTF) from a
    gene_id string, so ID-type auto-detection and bitr() lookup both
    see the real underlying accession (e.g. "ENSG00000310526") rather
    than a prefixed, unrecognizable string that always falls through to
    detect_id_type()'s own "SYMBOL" default.

    Returns the input unchanged if it doesn't start with "gene:" (e.g.
    a reference that never had this prefix issue at all, or a gene_id
    that's already bare).
    """
    if gene_id.startswith("gene:"):
        return gene_id[len("gene:"):]
    return gene_id


def _render_bitr_gene_name_fallback(project, adata):
    """
    Optional "fill in more gene names via Bioconductor bitr()" panel --
    single-cell's own analogue of the Bulk RNA-Seq Differential
    Expression workspace's "Gene ID -> Gene Name Mapping" panel
    (differential_expression_workspace._render_gene_id_mapping_panel),
    reusing gene_id_mapper.py directly (no changes to that module were
    needed or made -- it's already fully pipeline-agnostic; the
    "gene:"-prefix handling below lives entirely in THIS function,
    since it's specific to this project's own GFF3-fallback-derived
    references, not something gene_id_mapper.py itself should need to
    know about).

    Only offered for a PRESET reference (a real species_key is needed
    to look up the correct Bioconductor OrgDb annotation package) --
    for a custom/non-model-organism reference, there is no known OrgDb
    package to consult, so this section is skipped entirely rather than
    offering a lookup that would just fail.
    """
    reference_cfg = scpm.get_reference_choice(project) or {}
    species_key = reference_cfg.get("species_key")
    if reference_cfg.get("is_custom") or not species_key:
        return  # no known OrgDb package for a custom/non-model reference

    unresolved_ids = dsm.get_unresolved_gene_ids(adata.var)
    n_total = adata.n_vars

    with st.expander(
        f"🏷️ Fill in more gene names via Bioconductor bitr() "
        f"({len(unresolved_ids):,} of {n_total:,} genes still unresolved)",
        expanded=False,
    ):
        st.caption(
            "This project's reference GTF has already been used to resolve as "
            "many gene names as it directly contains (see Step 1's own combine "
            "summary) -- but not every gene_id in a GTF has a matching "
            "gene-level record with a real name (e.g. transcript/exon-only "
            "entries), so some genes remain stuck at their raw Ensembl ID. This "
            "runs a genuine Bioconductor annotation-database lookup "
            "(`clusterProfiler::bitr()`), completely independent of this "
            "project's own GTF, to try to resolve real names for whatever's "
            "still left."
        )

        if not unresolved_ids:
            st.success("✅ Every gene already has a resolved name -- nothing left to look up.")
            return

        if not gim.bitr_tools_available():
            orgdb_package = gim.ORGDB_PACKAGES.get(species_key, "the relevant OrgDb package")
            st.warning(
                f"⚠️ Rscript was not found on this system -- this needs R with "
                f"the `clusterProfiler` package and `{orgdb_package}` installed."
            )
            return

        # --- "gene:" prefix fix (2026-08-25) -- see this module's own
        # docstring, "The fix", for the full rationale. Detection AND
        # the eventual bitr() lookup both operate on the STRIPPED id;
        # results are mapped back onto the ORIGINAL (possibly prefixed)
        # gene_id further below.
        stripped_to_original = {}
        for gid in unresolved_ids:
            stripped_to_original.setdefault(_strip_gene_colon_prefix(gid), []).append(gid)
        stripped_ids = sorted(stripped_to_original.keys())

        n_stripped = sum(1 for gid in unresolved_ids if gid.startswith("gene:"))
        if n_stripped:
            st.caption(
                f"ℹ️ Detected a 'gene:' prefix on {n_stripped:,} of "
                f"{len(unresolved_ids):,} unresolved ID(s) (e.g. "
                f"`{unresolved_ids[0]}`) -- this is stripped automatically "
                "before ID-type detection and lookup, and the real result "
                "is applied back onto the original (prefixed) gene."
            )

        detection = gim.detect_id_type(stripped_ids)
        detected_type = detection["detected_type"]
        st.caption(
            f"Auto-detected ID type for your unresolved gene IDs: "
            f"**{detected_type}** ({detection['match_fraction'] * 100:.0f}% of a "
            f"sample matched -- e.g. `{'`, `'.join(detection['example_ids'][:3])}`)."
        )

        orgdb_package = gim.ORGDB_PACKAGES.get(species_key)
        default_to_type = gim.symbol_keytype_for_species(species_key)
        to_type_options = gim.COMMON_KEY_TYPES
        to_default_idx = to_type_options.index(default_to_type) if default_to_type in to_type_options else 0

        col_from, col_to = st.columns(2)
        with col_from:
            from_type_options = gim.COMMON_KEY_TYPES
            # Correctly pre-selects the AUTO-DETECTED type now that
            # detection itself runs on the stripped IDs -- previously
            # this always landed on SYMBOL's index regardless of what
            # the real underlying IDs were, since detection itself was
            # broken by the unstripped "gene:" prefix.
            from_default_idx = from_type_options.index(detected_type) if detected_type in from_type_options else 0
            picked_from_type = st.selectbox(
                "My unresolved gene IDs are currently in this format:",
                options=from_type_options, index=from_default_idx,
                key="sc_bitr_from_type_select",
            )
        with col_to:
            picked_to_type = st.selectbox(
                "Convert them to:", options=to_type_options, index=to_default_idx,
                key="sc_bitr_to_type_select",
            )

        if st.button(
            f"✨ Fill in {len(unresolved_ids):,} missing name(s)",
            key="sc_bitr_fill_missing_btn", type="primary",
        ):
            work_dir = scpm.downstream_bitr_work_dir(project)
            with st.spinner(f"Converting {len(stripped_ids):,} gene ID(s) via bitr()..."):
                result = gim.run_bitr_conversion(
                    stripped_ids, picked_from_type, picked_to_type, orgdb_package, work_dir,
                )
            if result["success"]:
                # Map bitr()'s own {stripped_id: converted_value}
                # mapping back onto every ORIGINAL (possibly "gene:"-
                # prefixed) gene_id it came from. CRITICAL: a gene
                # bitr() genuinely couldn't resolve (its own no-op
                # convention: converted_value == stripped_id) must map
                # back to the ORIGINAL prefixed id -- NOT the bare
                # stripped id -- so apply_bitr_symbol_mapping()'s own
                # "did this actually change?" check correctly treats it
                # as still fully unresolved, rather than silently
                # stripping the "gene:" prefix without ever finding a
                # real name and miscounting that as a resolution.
                original_keyed_mapping = {}
                for stripped_id, original_ids in stripped_to_original.items():
                    converted_value = result["mapping"].get(stripped_id, stripped_id)
                    for original_id in original_ids:
                        if converted_value == stripped_id:
                            original_keyed_mapping[original_id] = original_id
                        else:
                            original_keyed_mapping[original_id] = converted_value

                updated_var, n_newly_resolved = dsm.apply_bitr_symbol_mapping(adata.var, original_keyed_mapping)
                adata.var = updated_var
                _save_adata_state(project, adata)
                st.success(
                    f"✅ {result['message']} ({n_newly_resolved:,} gene(s) in "
                    "this project newly resolved.)"
                )
                st.rerun()
            else:
                st.error(f"⚠️ {result['message']}")


def _save_adata_state(project, adata):
    """
    Persist adata as this project's new "current state" -- updates
    BOTH the in-memory session_state copy AND overwrites the on-disk
    cache file (so a future session/restart can resume from here).
    """
    st.session_state[_session_key(project)] = adata
    path = scpm.downstream_adata_path(project)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    adata.write_h5ad(path)


def _recipe_default(project, step_name, param_name, fallback):
    """
    Read a previously-saved parameter value for step_name from this
    project's saved recipe, falling back to `fallback` if no recipe
    exists yet for this step, or that specific parameter isn't present
    in it.
    """
    saved = scpm.get_downstream_step_recipe(project, step_name)
    if saved and param_name in saved.get("params", {}):
        return saved["params"][param_name]
    return fallback


# ---------------------------------------------------------------------------
# Project selector
# ---------------------------------------------------------------------------

def _render_project_selector():
    st.header("Step 0: Choose a Single-cell Project")
    st.markdown(
        "Phase 3 (downstream analysis) builds on top of one or more samples "
        "that have ALREADY completed both Step 6 (alignment/cell-calling) "
        "AND Phase 2 (Cell-level QC) within the same project."
    )
    existing = scpm.list_projects()
    if not existing:
        st.info(
            "No single-cell projects found yet -- start one on the "
            "🧫 Single-cell RNA-Seq page first."
        )
        return None

    session_key = "sc_downstream_active_project"
    default_index = 0
    if session_key in st.session_state and st.session_state[session_key] in existing:
        default_index = existing.index(st.session_state[session_key])

    chosen = st.selectbox(
        "Project:", options=existing, index=default_index, key="sc_downstream_project_select",
    )
    st.session_state[session_key] = chosen
    st.success(f"**Active project: `{chosen}`**")
    return chosen


# ---------------------------------------------------------------------------
# Step 1: Combine Samples
# ---------------------------------------------------------------------------

def _render_combine_step(project):
    st.header("Step 1: Combine Samples")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Loads each selected sample's STARsolo filtered "
            "cell x gene matrix, attaches that sample's Phase 2 Cell-level QC "
            "flags (doublet calls, ambient-RNA contamination, adaptive-QC "
            "pass/fail) as per-cell metadata, and concatenates every sample "
            "into ONE combined dataset -- one row per cell, one column per "
            "gene, tagged with which original sample it came from.\n\n"
            "**Why it matters:** Every later Phase 3 step (normalization, "
            "clustering, batch correction, differential composition, etc.) "
            "operates on cells from ALL selected samples at once, so cell "
            "types/states can be compared directly across conditions/donors "
            "rather than analyzed in isolation per sample. This is also the "
            "step where a genuinely bad sample (e.g. one that failed QC "
            "outright) can simply be left unselected, without needing to "
            "redo any earlier pipeline steps.\n\n"
            "**How to interpret the result:** The reported cell x gene "
            "dimensions are your starting point for everything downstream -- "
            "a MUCH lower cell count than expected (e.g. far fewer cells than "
            "the sum of each sample's own Step 6 'Cells Detected') usually "
            "means the QC filters above are excluding more cells than "
            "intended, worth a quick sanity check before proceeding."
        )

    eligible_samples = scpm.get_samples_with_completed_cellqc(project)
    if not eligible_samples:
        st.warning(
            "⚠️ No samples in this project have BOTH completed alignment "
            "(Step 6) AND completed Cell-level QC (Phase 2) yet -- at "
            "least one sample must have both before Phase 3 can begin. "
            "Complete these on the 🧫 Single-cell RNA-Seq / 🔬 SC "
            "Cell-level QC pages first."
        )
        return None

    st.caption(
        f"Found {len(eligible_samples)} sample(s) with completed alignment + "
        f"Cell-level QC: {', '.join(eligible_samples)}"
    )

    default_samples = _recipe_default(project, "combine", "selected_samples", eligible_samples)
    default_samples = [s for s in default_samples if s in eligible_samples] or eligible_samples

    selected_samples = st.multiselect(
        "Which sample(s) should be combined for downstream analysis?",
        options=eligible_samples, default=default_samples,
        key="sc_downstream_combine_samples_select",
    )
    if not selected_samples:
        st.info("Select at least one sample above to continue.")
        return None
    if len(selected_samples) == 1:
        st.info(
            "ℹ️ Only one sample selected -- batch correction (a later step) "
            "will have nothing to correct for and should be set to 'skip'."
        )

    default_apply_qc = _recipe_default(project, "combine", "apply_qc_filters", True)
    apply_qc_filters = st.checkbox(
        "Apply Phase 2's QC filters when loading (drop adaptive-QC-failed "
        "and predicted-doublet cells)",
        value=default_apply_qc, key="sc_downstream_apply_qc_filters",
        help=(
            "If unchecked, every cell in each sample's filtered matrix is "
            "loaded regardless of Phase 2's QC flags -- those flags remain "
            "available as columns for your own inspection/filtering later, "
            "rather than being silently baked in."
        ),
    )

    default_min_cells = _recipe_default(
        project, "combine", "min_cells_per_gene", dsm.DEFAULT_MIN_CELLS_PER_GENE,
    )
    min_cells_per_gene = st.number_input(
        "Minimum cells a gene must be detected in (across the WHOLE combined "
        "dataset) to be kept:",
        min_value=0, max_value=100, value=int(default_min_cells),
        key="sc_downstream_min_cells_per_gene",
        help=(
            "A standard, lightweight gene-level filter applied AFTER "
            "combining -- distinct from Phase 2's per-cell filtering."
        ),
    )

    current_params = {
        "selected_samples": sorted(selected_samples),
        "apply_qc_filters": apply_qc_filters,
        "min_cells_per_gene": int(min_cells_per_gene),
    }
    is_current = scpm.check_downstream_step_current(project, "combine", current_params)
    cached_adata = _get_cached_adata(project)
    combine_done = cached_adata is not None and is_current

    if combine_done:
        st.success(
            f"✅ Already combined with these exact settings -- "
            f"{cached_adata.n_obs:,} cells x {cached_adata.n_vars:,} genes."
        )
    elif cached_adata is not None and not is_current:
        st.warning(
            "⚠️ Your current sample/filter selections differ from what's "
            "cached -- click below to re-combine with the new settings. "
            "(Re-combining will also invalidate every later Phase 3 step "
            "already computed for this project.)"
        )

    run_label = "🔄 Re-combine Samples" if cached_adata is not None else "🚀 Combine Samples"
    if st.button(run_label, key="sc_downstream_combine_btn", type="primary"):
        sample_specs = [
            scpm.get_sample_spec_for_downstream(project, s) for s in selected_samples
        ]
        # --- Gene-symbol resolution overlay (2026-08-25) -- see this
        # module's own docstring, "Gene-symbol resolution wiring", and
        # sc_downstream_manager.py's own module docstring for the full
        # rationale: this project's already-confirmed reference GTF
        # (Step 5) is passed through so any gene whose STARsolo-
        # provided symbol still looks unresolved (e.g. a raw
        # "gene:ENSG..." form) can be corrected here, at combine time --
        # WITHOUT requiring the sample to be re-aligned from scratch.
        reference_cfg = scpm.get_reference_choice(project) or {}
        gtf_path = reference_cfg.get("custom_gtf")
        with st.spinner(f"Loading and combining {len(selected_samples)} sample(s)..."):
            try:
                combined, per_sample_summary = dsm.load_and_combine_samples(
                    sample_specs, apply_qc_filters=apply_qc_filters,
                    min_cells_per_gene=int(min_cells_per_gene),
                    gtf_path=gtf_path,
                )
            except FileNotFoundError as e:
                st.error(f"⚠️ {e}")
                return None
        _save_adata_state(project, combined)
        scpm.save_downstream_step_recipe(project, "combine", current_params)
        st.session_state["sc_downstream_combine_summary"] = per_sample_summary
        st.success(
            f"✅ Combined {len(selected_samples)} sample(s) -- "
            f"{combined.n_obs:,} cells x {combined.n_vars:,} genes."
        )
        st.rerun()

    summary = st.session_state.get("sc_downstream_combine_summary")
    if summary:
        st.markdown("**Per-sample cell counts (this session's most recent combine):**")
        st.dataframe(pd.DataFrame(summary), use_container_width=True, hide_index=True)
        st.caption(
            "💡 **How to interpret this table:** compare each sample's cell "
            "count here against its own 'Cells Detected' from Step 6 -- a "
            "large drop is expected if QC filtering is enabled above (that's "
            "the point of Phase 2 QC), but a drop of more than ~30-50% for "
            "any one sample is worth a second look; it may indicate that "
            "sample had unusually poor quality relative to the others."
        )
        
    if combine_done:
        _render_bitr_gene_name_fallback(project, cached_adata)
    return cached_adata if combine_done else None
  
SAMPLE_METADATA_COLUMNS_UNS_KEY = "sample_metadata_columns"


def _render_sample_metadata_step(project, adata):
    st.header("Step 1b: Label Samples with Metadata")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Attaches sample-level metadata (e.g. experimental "
            "condition, donor, batch) onto every cell belonging to that sample -- "
            "so a column like 'condition' becomes available throughout the rest "
            "of this workspace, exactly the same way 'sample' itself already is.\n\n"
            "**Why it matters:** Step 9's pseudobulk export and Step 10's "
            "compositional analysis BOTH need a real experimental condition "
            "column to compare across (e.g. 'Healthy' vs. 'COVID-19') -- without "
            "this step, there is nothing valid to choose from those dropdowns "
            "yet, since Phase 2's own per-cell QC columns (e.g. 'sum', a total "
            "UMI count) are NOT sample-level conditions and cannot be used this "
            "way (they vary cell-to-cell, even within one sample).\n\n"
            "**Where this metadata can come from:** if you already filled in a "
            "metadata table during this project's Step 1 ingestion (FASTQ "
            "ingestion + chemistry + metadata, on the 🧫 Single-cell RNA-Seq "
            "page), you can import it directly below with one click. You can "
            "also edit sample metadata manually, or upload a fresh file -- "
            "whichever is easiest."
        )

    sample_key = "sample"
    if sample_key not in adata.obs.columns:
        st.warning("⚠️ No 'sample' column found in this dataset -- this step cannot run.")
        return

    unique_samples = sorted(adata.obs[sample_key].astype(str).unique())
    existing_columns = list(adata.uns.get(SAMPLE_METADATA_COLUMNS_UNS_KEY, []))

    if existing_columns:
        st.success(
            f"✅ This project currently has {len(existing_columns)} sample-level "
            f"metadata column(s) attached: **{', '.join(existing_columns)}**."
        )
        preview_df = adata.obs[[sample_key] + existing_columns].drop_duplicates(subset=[sample_key])
        preview_df = preview_df.sort_values(sample_key).reset_index(drop=True)
        with st.expander("👀 Current sample metadata", expanded=False):
            st.dataframe(preview_df, use_container_width=True, hide_index=True)
    else:
        st.info(
            "ℹ️ No sample-level metadata attached yet -- add at least a "
            "'condition' column below before Step 10 (Compositional Analysis) "
            "can run."
        )

    working_key = f"sc_downstream_sample_metadata_working_df_{project}"
    if working_key not in st.session_state:
        if existing_columns:
            base_df = adata.obs[[sample_key] + existing_columns].drop_duplicates(subset=[sample_key])
            base_df = base_df.sort_values(sample_key).reset_index(drop=True)
        else:
            base_df = pd.DataFrame({sample_key: unique_samples})
        st.session_state[working_key] = base_df

    # --- Option 1: import directly from this project's own Phase 1
    # ingestion metadata.csv, if it exists.
    phase1_metadata_path = scpm.metadata_path(project)
    if os.path.isfile(phase1_metadata_path):
        with st.expander("📥 Import from this project's own ingestion metadata (Step 1)", expanded=not bool(existing_columns)):
            phase1_df = pd.read_csv(phase1_metadata_path)
            st.caption(
                f"Found `{os.path.basename(phase1_metadata_path)}` from this project's "
                f"own Step 1 ingestion -- preview below."
            )
            st.dataframe(phase1_df, use_container_width=True, hide_index=True)
            if sample_key not in phase1_df.columns:
                st.warning(
                    f"⚠️ This file has no column named exactly '{sample_key}' -- cannot "
                    f"import automatically. Use manual entry or file upload below instead."
                )
            else:
                importable_columns = [c for c in phase1_df.columns if c != sample_key]
                cols_to_import = st.multiselect(
                    "Which column(s) should be imported as sample metadata?",
                    options=importable_columns, default=importable_columns,
                    key="sc_downstream_sample_meta_phase1_cols",
                )
                if cols_to_import and st.button(
                    "📥 Import Selected Column(s)", key="sc_downstream_sample_meta_import_phase1_btn",
                ):
                    import_df = phase1_df[[sample_key] + cols_to_import].drop_duplicates(subset=[sample_key])
                    current_working = st.session_state[working_key]
                    merged_working = current_working.merge(import_df, on=sample_key, how="outer", suffixes=("", "_imported"))
                    for col in cols_to_import:
                        imported_col = f"{col}_imported"
                        if imported_col in merged_working.columns:
                            merged_working[col] = merged_working[imported_col].combine_first(merged_working.get(col))
                            merged_working = merged_working.drop(columns=[imported_col])
                    st.session_state[working_key] = merged_working
                    st.success(f"✅ Imported {len(cols_to_import)} column(s) -- review below, then save.")
                    st.rerun()

    # --- Option 2: upload a fresh CSV/TXT/XLSX file (sample + arbitrary columns).
    with st.expander("📤 Upload a sample metadata file", expanded=False):
        st.caption(
            "Upload a spreadsheet with a column named exactly 'sample' (matching "
            "this project's own sample names) plus any other columns you want "
            "(e.g. 'condition', 'donor', 'batch')."
        )
        uploaded_file = st.file_uploader(
            "Upload a metadata file:", type=["csv", "txt", "xlsx", "xls"],
            key="sc_downstream_sample_meta_upload",
        )
        if uploaded_file is not None:
            uploaded_df, error = ing.read_metadata_file(uploaded_file)
            if error:
                st.error(f"⚠️ {error}")
            elif sample_key not in uploaded_df.columns:
                st.error(f"⚠️ Uploaded file must have a column named exactly '{sample_key}'.")
            else:
                st.dataframe(uploaded_df, use_container_width=True, hide_index=True)
                if st.button("📥 Merge Uploaded File", key="sc_downstream_sample_meta_upload_merge_btn"):
                    current_working = st.session_state[working_key]
                    merged_working = current_working.merge(uploaded_df, on=sample_key, how="outer", suffixes=("", "_uploaded"))
                    for col in uploaded_df.columns:
                        if col == sample_key:
                            continue
                        uploaded_col = f"{col}_uploaded"
                        if uploaded_col in merged_working.columns:
                            merged_working[col] = merged_working[uploaded_col].combine_first(merged_working.get(col))
                            merged_working = merged_working.drop(columns=[uploaded_col])
                    st.session_state[working_key] = merged_working
                    st.success("✅ Merged uploaded file -- review below, then save.")
                    st.rerun()

    # --- Option 3: direct manual editing.
    st.markdown("**✏️ Edit sample metadata directly:**")
    working_df = st.session_state[working_key]
    edited_df = st.data_editor(
        working_df, use_container_width=True, hide_index=True, num_rows="fixed",
        key="sc_downstream_sample_metadata_editor",
        column_config={sample_key: st.column_config.TextColumn(disabled=True)},
    )
    st.session_state[working_key] = edited_df

    with st.expander("➕ Add a new column (e.g. 'condition', 'donor', 'batch')"):
        new_col_name = st.text_input("New column name:", key="sc_downstream_sample_meta_new_col_name")
        if st.button("➕ Add Column", key="sc_downstream_sample_meta_add_col_btn"):
            clean_name = new_col_name.strip()
            if not clean_name:
                st.error("⚠️ Enter a column name first.")
            elif clean_name in edited_df.columns:
                st.error(f"⚠️ A column named `{clean_name}` already exists.")
            else:
                new_df = edited_df.copy()
                new_df[clean_name] = ""
                st.session_state[working_key] = new_df
                st.rerun()

    if st.button("💾 Save Sample Metadata", key="sc_downstream_save_sample_metadata_btn", type="primary"):
        final_df = st.session_state[working_key]
        metadata_columns = [c for c in final_df.columns if c != sample_key]
        if not metadata_columns:
            st.error("⚠️ Add at least one metadata column before saving.")
        else:
            adata, report = dsm.merge_sample_level_metadata(adata, final_df, sample_col=sample_key)
            all_columns = sorted(set(existing_columns) | set(metadata_columns))
            adata.uns[SAMPLE_METADATA_COLUMNS_UNS_KEY] = all_columns
            _save_adata_state(project, adata)
            st.success(f"✅ Saved sample metadata -- columns now available: {', '.join(all_columns)}.")
            if report["samples_missing_metadata"]:
                st.warning(
                    f"⚠️ {len(report['samples_missing_metadata'])} sample(s) in this dataset have NO "
                    f"metadata row and will show blank/NaN values: "
                    f"{', '.join(report['samples_missing_metadata'])}"
                )
            if report["samples_not_found_in_data"]:
                st.info(
                    f"ℹ️ {len(report['samples_not_found_in_data'])} metadata row(s) didn't match any "
                    f"sample actually in this dataset (ignored): "
                    f"{', '.join(report['samples_not_found_in_data'])}"
                )
            st.rerun()
  


# ---------------------------------------------------------------------------
# Step 2: Normalization
# ---------------------------------------------------------------------------

def _render_normalize_step(project, adata):
    st.header("Step 2: Normalization")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Corrects for the fact that different cells "
            "were sequenced to different depths (some cells simply got more "
            "UMI counts than others, for entirely technical reasons unrelated "
            "to biology) so that expression values become comparable "
            "CELL-TO-CELL. Two approaches are offered: a classic 'shifted "
            "log' transform (scale each cell to a common total, then "
            "log-transform), or 'Pearson residuals' (a model-based approach "
            "that doesn't require log-transforming and tends to preserve "
            "biological variance more faithfully for highly-expressed genes).\n\n"
            "**Why it matters:** Without this step, a cell that happened to "
            "be sequenced twice as deeply as another would LOOK like it's "
            "expressing every gene twice as much -- a purely technical "
            "artifact that would otherwise dominate every downstream "
            "analysis (HVG selection, PCA, clustering) and could easily be "
            "mistaken for a real biological difference between cell "
            "populations.\n\n"
            "**How to interpret the result:** There isn't a plot to check "
            "here -- normalization succeeding just means a new expression "
            "layer (`lognorm` or `pearson_residuals`) now exists on the "
            "dataset, ready to feed into HVG selection and PCA next. If "
            "you're unsure which method to pick, 'shifted log' is the more "
            "traditional, widely-used default and a safe first choice."
        )

    method_keys = list(dsm.NORMALIZATION_METHOD_OPTIONS.keys())
    default_method = _recipe_default(project, "normalize", "method", dsm.DEFAULT_NORMALIZATION_METHOD)
    method = st.radio(
        "Normalization method:", method_keys,
        format_func=lambda k: dsm.NORMALIZATION_METHOD_OPTIONS[k]["label"],
        index=method_keys.index(default_method) if default_method in method_keys else 0,
        key="sc_downstream_norm_method",
    )
    st.caption(dsm.NORMALIZATION_METHOD_OPTIONS[method]["explanation"])

    target_sum = None
    pearson_theta = 100
    if method == "shifted_log":
        default_target_sum = _recipe_default(project, "normalize", "target_sum", None)
        use_custom_target = st.checkbox(
            "Use a custom target sum (advanced) instead of the dataset's "
            "own median total count",
            value=(default_target_sum is not None), key="sc_downstream_norm_custom_target",
        )
        if use_custom_target:
            target_sum = st.number_input(
                "Target sum:", min_value=100.0, max_value=1_000_000.0,
                value=float(default_target_sum) if default_target_sum else 10000.0,
                step=100.0, key="sc_downstream_norm_target_sum",
            )
    else:
        default_theta = _recipe_default(project, "normalize", "pearson_theta", 100)
        pearson_theta = st.number_input(
            "Overdispersion parameter (theta):", min_value=1.0, max_value=1000.0,
            value=float(default_theta), step=10.0, key="sc_downstream_norm_theta",
        )

    current_params = {"method": method, "target_sum": target_sum, "pearson_theta": pearson_theta}
    is_current = scpm.check_downstream_step_current(project, "normalize", current_params)

    layer_name = "lognorm" if method == "shifted_log" else "pearson_residuals"
    already_has_layer = layer_name in adata.layers

    if is_current and already_has_layer:
        st.success(f"✅ Already normalized with these exact settings (layer: `{layer_name}`).")
    elif already_has_layer and not is_current:
        st.warning(
            "⚠️ Your current normalization settings differ from what's "
            "cached -- re-run below to apply them. (This will also "
            "invalidate every later Phase 3 step already computed.)"
        )

    run_label = "🔄 Re-run Normalization" if already_has_layer else "▶️ Run Normalization"
    if st.button(run_label, key="sc_downstream_normalize_btn", type="primary"):
        with st.spinner("Normalizing..."):
            dsm.normalize(adata, method=method, target_sum=target_sum, pearson_theta=pearson_theta)
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "normalize", current_params)
        st.success(f"✅ Normalization complete (layer: `{layer_name}`).")
        st.rerun()

    return is_current and already_has_layer


# ---------------------------------------------------------------------------
# Step 3: Highly Variable Gene (HVG) Selection
# ---------------------------------------------------------------------------

def _render_hvg_step(project, adata):
    st.header("Step 3: Highly Variable Gene (HVG) Selection")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Flags a subset of genes (typically 1,000-5,000, "
            "out of tens of thousands measured) whose expression varies the "
            "MOST across cells, above and beyond what's expected from "
            "technical/random noise alone.\n\n"
            "**Why it matters:** Most genes in a typical dataset are either "
            "expressed at a nearly constant, 'housekeeping' level across "
            "every cell, or are simply too lowly/sparsely detected to carry "
            "reliable signal -- including them in PCA/clustering mostly adds "
            "noise, dilutes real biological signal, and slows computation "
            "down substantially. Restricting to the most variable genes lets "
            "PCA focus on the genes actually likely to distinguish different "
            "cell types/states from each other.\n\n"
            "**How to interpret the plot below:** Each point is one gene, "
            "positioned by how highly expressed it is (x-axis) and how much "
            "MORE variable it is than expected for that expression level "
            "(y-axis) -- genes selected as 'highly variable' (red) should "
            "form a cloud clearly above the bulk of genes (blue) at every "
            "expression level, not just the highest-expression genes alone. "
            "If the selected genes look like they're only capturing the "
            "very highest-expressed genes with no separation at lower "
            "expression levels, that can indicate too few genes were "
            "requested, or this data has unusually flat variance overall."
        )

    default_n_top = _recipe_default(project, "hvg", "n_top_genes", dsm.DEFAULT_N_TOP_GENES)
    n_top_genes = st.slider(
        "Number of top variable genes to select:", min_value=200, max_value=10000,
        value=int(default_n_top), step=100, key="sc_downstream_hvg_n_top",
    )

    n_samples_in_data = adata.obs["sample"].nunique() if "sample" in adata.obs.columns else 1
    batch_key = None
    if n_samples_in_data > 1:
        default_batch_key = _recipe_default(project, "hvg", "batch_key", "sample")
        use_batch_key = st.checkbox(
            "Compute HVGs per-sample, then combine (recommended for 2+ samples)",
            value=(default_batch_key is not None), key="sc_downstream_hvg_use_batch",
            help=(
                "When checked, HVG selection is computed separately per "
                "sample and then combined (scanpy's own batch-aware "
                "ranking), so a gene that's highly variable in only one "
                "sample doesn't get selected purely due to a sample-"
                "specific effect."
            ),
        )
        batch_key = "sample" if use_batch_key else None

    current_params = {"n_top_genes": n_top_genes, "batch_key": batch_key}
    is_current = scpm.check_downstream_step_current(project, "hvg", current_params)
    already_has_hvg = "highly_variable" in adata.var.columns

    if is_current and already_has_hvg:
        st.success(
            f"✅ Already selected with these exact settings -- "
            f"{int(adata.var['highly_variable'].sum()):,} HVGs."
        )
    elif already_has_hvg and not is_current:
        st.warning(
            "⚠️ Your current HVG settings differ from what's cached -- "
            "re-run below to apply them. (This will also invalidate PCA "
            "and every later step already computed.)"
        )

    run_label = "🔄 Re-run HVG Selection" if already_has_hvg else "▶️ Select Highly Variable Genes"
    if st.button(run_label, key="sc_downstream_hvg_btn", type="primary"):
        if "lognorm" not in adata.layers and "pearson_residuals" not in adata.layers:
            st.error("⚠️ Run normalization (Step 2) first.")
            return is_current and already_has_hvg
        with st.spinner("Selecting highly variable genes..."):
            dsm.select_highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "hvg", current_params)
        st.success(f"✅ Selected {int(adata.var['highly_variable'].sum()):,} highly variable genes.")
        st.rerun()

    if already_has_hvg:
        with st.expander("📊 HVG diagnostic plot", expanded=True):
            hvg_df = dsm.get_hvg_plot_data(adata)
            if hvg_df is not None:
                x_col = "means" if "means" in hvg_df.columns else None
                y_col = (
                    "dispersions_norm" if "dispersions_norm" in hvg_df.columns
                    else ("variances_norm" if "variances_norm" in hvg_df.columns else None)
                )
                if x_col and y_col:
                    plot_df = hvg_df.copy()
                    plot_df["highly_variable"] = plot_df["highly_variable"].astype(str)
                    fig = px.scatter(
                        plot_df, x=x_col, y=y_col, color="highly_variable",
                        color_discrete_map={"True": "#d62728", "False": "#636EFA"},
                        labels={x_col: "Mean expression", y_col: "Normalized dispersion/variance",
                                "highly_variable": "Highly variable"},
                        opacity=0.5,
                    )
                    fig.update_layout(height=450, margin=dict(l=10, r=10, t=30, b=10))
                    style = scw._render_plot_style_controls(
                        "sc_downstream_hvg_plot", group_values=["True", "False"],
                        default_colors={"True": "#d62728", "False": "#636EFA"},
                    )
                    scw._apply_plot_style(
                        fig, style, default_title="HVG selection",
                        default_x_label="Mean expression", default_y_label="Normalized dispersion",
                    )
                    scw._render_plotly_chart(fig)
                    st.caption(
                        "💡 **How to read this:** red points (selected HVGs) should sit "
                        "ABOVE the blue cloud (non-selected genes) at every expression "
                        "level, not just the right-hand (highly expressed) side of the "
                        "plot -- that separation across the full x-axis range is the "
                        "actual sign that variable genes are being captured correctly."
                    )
                    scw._render_pdf_export(fig, "sc_downstream_hvg_plot", "hvg_diagnostic_plot")
                    scw._render_csv_download(
                        hvg_df, "hvg_diagnostic_data", "sc_downstream_hvg_plot",
                        expander_label="⬇️ Download HVG diagnostic data (.csv)",
                    )

    return is_current and already_has_hvg


# ---------------------------------------------------------------------------
# Step 4: PCA
# ---------------------------------------------------------------------------

def _render_pca_step(project, adata):
    st.header("Step 4: Principal Component Analysis (PCA)")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Compresses the (potentially thousands of) "
            "highly variable genes from Step 3 down into a much smaller "
            "number of 'principal components' (PCs) -- each PC is a "
            "computed combination of many genes, ordered so that PC1 "
            "captures the single largest axis of variation in the data, "
            "PC2 the next-largest INDEPENDENT axis, and so on.\n\n"
            "**Why it matters:** Nearly every later step in this pipeline "
            "(batch correction, the neighbor graph, clustering, UMAP/t-SNE) "
            "operates on this PCA embedding rather than on raw gene "
            "expression directly -- it's a dramatically more compact, less "
            "noisy representation of the same underlying biological "
            "signal, and is what makes clustering thousands of cells "
            "computationally practical at all.\n\n"
            "**How to interpret the plots below:** The elbow plot shows how "
            "much variance each successive PC explains -- look for the "
            "point where the curve visibly flattens ('the elbow'); PCs "
            "beyond that point are typically capturing mostly noise, and "
            "that elbow location is a reasonable guide for how many PCs to "
            "actually use in later steps (clustering/embedding). The "
            "scatter plot lets you visually check for **batch effects**: "
            "if cells from different samples separate into distinct, "
            "non-overlapping clouds along early PCs even before any real "
            "biological clustering, that's a concrete signal to use "
            "Step 5's batch correction before proceeding."
        )

    default_use_hvg = _recipe_default(project, "pca", "use_highly_variable", True)
    use_hvg = st.checkbox(
        "Restrict PCA to highly variable genes only (recommended)",
        value=default_use_hvg, key="sc_downstream_pca_use_hvg",
    )
    default_n_comps = _recipe_default(project, "pca", "n_comps", dsm.DEFAULT_N_PCS)
    n_comps = st.slider(
        "Number of principal components:", min_value=10, max_value=100,
        value=int(default_n_comps), step=5, key="sc_downstream_pca_n_comps",
    )

    current_params = {"n_comps": n_comps, "use_highly_variable": use_hvg}
    is_current = scpm.check_downstream_step_current(project, "pca", current_params)
    already_has_pca = "X_pca" in adata.obsm

    if is_current and already_has_pca:
        st.success(
            f"✅ Already computed with these exact settings -- "
            f"{adata.obsm['X_pca'].shape[1]} PCs."
        )
    elif already_has_pca and not is_current:
        st.warning(
            "⚠️ Your current PCA settings differ from what's cached -- "
            "re-run below. (This will also invalidate every later step "
            "already computed.)"
        )

    run_label = "🔄 Re-run PCA" if already_has_pca else "▶️ Run PCA"
    if st.button(run_label, key="sc_downstream_pca_btn", type="primary"):
        if use_hvg and "highly_variable" not in adata.var.columns:
            st.error(
                "⚠️ Run HVG selection (Step 3) first, or uncheck "
                "'restrict to highly variable genes'."
            )
            return is_current and already_has_pca
        with st.spinner("Running PCA..."):
            dsm.run_pca(adata, n_comps=n_comps, use_highly_variable=use_hvg)
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "pca", current_params)
        st.success("✅ PCA complete.")
        st.rerun()

    if already_has_pca:
        with st.expander("📊 Variance explained (elbow plot)", expanded=True):
            variance_ratio = dsm.get_pca_variance_ratio(adata)
            if variance_ratio:
                elbow_df = pd.DataFrame({
                    "PC": [f"PC{i+1}" for i in range(len(variance_ratio))],
                    "Variance explained": variance_ratio,
                })
                fig = go.Figure(go.Scatter(
                    x=elbow_df["PC"], y=elbow_df["Variance explained"],
                    mode="lines+markers", name="Variance explained",
                ))
                fig.update_layout(height=400, margin=dict(l=10, r=10, t=30, b=10))
                style = scw._render_plot_style_controls(
                    "sc_downstream_pca_elbow",
                    single_color_default={"Variance explained": "#636EFA"},
                    show_legend_controls=False,
                )
                scw._apply_plot_style(
                    fig, style, default_title="PCA variance explained",
                    default_x_label="Principal component",
                    default_y_label="Fraction of variance explained",
                )
                scw._render_plotly_chart(fig)
                st.caption(
                    "💡 **How to read this:** find the 'elbow' -- the point where the "
                    "curve stops dropping steeply and flattens out. PCs to the LEFT of "
                    "the elbow capture real, substantial biological variation; PCs to "
                    "the RIGHT are typically dominated by noise. This elbow location is "
                    "a reasonable estimate for how many PCs to actually carry into "
                    "Step 6 (Clustering) and Step 7 (Embedding), rather than blindly "
                    "using every PC computed here."
                )
                scw._render_pdf_export(fig, "sc_downstream_pca_elbow", "pca_elbow_plot")
                scw._render_csv_download(
                    elbow_df, "pca_variance_explained", "sc_downstream_pca_elbow",
                    expander_label="⬇️ Download variance-explained data (.csv)",
                )

        with st.expander("📊 PCA scatter (choose which PCs to view)", expanded=True):
            n_pcs_available = adata.obsm["X_pca"].shape[1]
            pc_options = [f"PC{i+1}" for i in range(n_pcs_available)]

            axis_col1, axis_col2, axis_col3, axis_col4 = st.columns([1, 1, 1, 1])
            with axis_col1:
                x_pc = st.selectbox(
                    "X axis:", options=pc_options, index=0, key="sc_downstream_pca_scatter_x_pc",
                )
            with axis_col2:
                y_pc = st.selectbox(
                    "Y axis:", options=pc_options,
                    index=1 if len(pc_options) > 1 else 0,
                    key="sc_downstream_pca_scatter_y_pc",
                )
            with axis_col3:
                use_z_axis = st.checkbox(
                    "Add Z axis (3D)", value=False, key="sc_downstream_pca_scatter_use_z",
                    help="Enable to view a THIRD principal component as a Z axis, rendered as an interactive 3D scatter plot instead of a flat 2D one.",
                )
            with axis_col4:
                z_pc = None
                if use_z_axis:
                    z_pc = st.selectbox(
                        "Z axis:", options=pc_options,
                        index=2 if len(pc_options) > 2 else 0,
                        key="sc_downstream_pca_scatter_z_pc",
                    )

            n_pcs_needed = max(
                int(x_pc.replace("PC", "")), int(y_pc.replace("PC", "")),
                int(z_pc.replace("PC", "")) if z_pc else 0,
            )
            pca_df = dsm.get_pca_coordinates(
                adata, n_pcs_to_include=min(n_pcs_needed, n_pcs_available), color_by_columns=["sample"],
            )
            if pca_df is not None and x_pc in pca_df.columns and y_pc in pca_df.columns:
                has_sample_col = "sample" in pca_df.columns
                sample_values = dsm.natural_sort_unique(pca_df["sample"].astype(str)) if has_sample_col else None

                if use_z_axis and z_pc and z_pc in pca_df.columns:
                    fig2 = px.scatter_3d(
                        pca_df, x=x_pc, y=y_pc, z=z_pc,
                        color="sample" if has_sample_col else None, opacity=0.6,
                        category_orders=({"sample": sample_values} if sample_values else None),
                    )
                    default_title = f"PCA ({x_pc} vs. {y_pc} vs. {z_pc})"
                else:
                    fig2 = px.scatter(
                        pca_df, x=x_pc, y=y_pc,
                        color="sample" if has_sample_col else None, opacity=0.6,
                        category_orders=({"sample": sample_values} if sample_values else None),
                    )
                    default_title = f"PCA ({x_pc} vs. {y_pc})"

                fig2.update_layout(height=550, margin=dict(l=10, r=10, t=30, b=10))
                style2 = scw._render_plot_style_controls(
                    "sc_downstream_pca_scatter", group_values=sample_values,
                    show_legend_controls=bool(sample_values),
                )
                scw._apply_plot_style(
                    fig2, style2, default_title=default_title,
                    default_x_label=x_pc, default_y_label=y_pc,
                )
                scw._render_plotly_chart(fig2)
                st.caption(
                    "💡 **How to read this:** if points from different samples (colors) "
                    "form separate, non-overlapping clouds here -- rather than mixing "
                    "together -- that's a real batch effect signal, and Step 5 (Batch "
                    "Correction) should be used before clustering. This isn't limited "
                    "to PC1/PC2 -- a batch effect (or a rare, distinct cell population) "
                    "can sometimes only become visible on a LATER pair of PCs, which is "
                    "why every PC is selectable above rather than only the first two."
                )
                scw._render_pdf_export(fig2, "sc_downstream_pca_scatter", "pca_scatter_plot")
                scw._render_csv_download(
                    pca_df, "pca_coordinates", "sc_downstream_pca_scatter",
                    expander_label="⬇️ Download PCA coordinates (.csv)",
                )

    return is_current and already_has_pca


# ---------------------------------------------------------------------------
# Shared: named checkpoints + .h5ad download
# ---------------------------------------------------------------------------

def _render_checkpoint_controls(project, adata):
    with st.expander("💾 Save / Load Named Checkpoints (optional)"):
        st.caption(
            "Unlike the automatic cache above (which is overwritten every "
            "time a step completes), a checkpoint preserves a SPECIFIC "
            "configuration under a name you choose -- useful for comparing "
            "e.g. 'HVG=2000, Leiden res=1.0' against a different parameter "
            "combination later, without losing either."
        )
        col1, col2 = st.columns(2)
        with col1:
            checkpoint_name = st.text_input("Checkpoint name:", key="sc_downstream_checkpoint_name")
            if st.button("💾 Save Current State as Checkpoint", key="sc_downstream_save_checkpoint_btn"):
                if not checkpoint_name.strip():
                    st.error("Enter a checkpoint name first.")
                else:
                    try:
                        path = scpm.downstream_checkpoint_path(project, checkpoint_name)
                        os.makedirs(os.path.dirname(path), exist_ok=True)
                        with st.spinner("Saving checkpoint..."):
                            adata.write_h5ad(path)
                        st.success(f"✅ Saved checkpoint '{checkpoint_name}'.")
                        st.rerun()
                    except ValueError as e:
                        st.error(f"⚠️ {e}")
        with col2:
            existing_checkpoints = scpm.list_downstream_checkpoints(project)
            if existing_checkpoints:
                chosen_checkpoint = st.selectbox(
                    "Load a checkpoint:", options=existing_checkpoints,
                    key="sc_downstream_load_checkpoint_select",
                )
                if st.button(
                    "📂 Load This Checkpoint (replaces current state)",
                    key="sc_downstream_load_checkpoint_btn",
                ):
                    import anndata as ad
                    path = scpm.downstream_checkpoint_path(project, chosen_checkpoint)
                    with st.spinner("Loading checkpoint..."):
                        loaded = ad.read_h5ad(path)
                    _save_adata_state(project, loaded)
                    st.warning(
                        "⚠️ Loading a checkpoint replaces the current state, but does NOT "
                        "automatically restore this project's saved step recipes to match it -- "
                        "the parameter widgets above may show a 'settings differ from cache' "
                        "warning until you re-select the same settings the checkpoint was built "
                        "with, or simply re-run each step to establish a fresh, matching recipe."
                    )
                    st.success(f"✅ Loaded checkpoint '{chosen_checkpoint}' as the current state.")
                    st.rerun()
                if st.button("🗑️ Delete This Checkpoint", key="sc_downstream_delete_checkpoint_btn"):
                    scpm.delete_downstream_checkpoint(project, chosen_checkpoint)
                    st.success(f"✅ Deleted checkpoint '{chosen_checkpoint}'.")
                    st.rerun()
            else:
                st.caption("No checkpoints saved yet for this project.")


def _render_download_section(project, adata):
    st.markdown("**⬇️ Download Current State**")
    st.caption(
        "Download this project's combined AnnData object in its CURRENT "
        "state (.h5ad -- the standard single-cell data interchange format) "
        "-- includes every layer, embedding, and annotation computed so "
        "far. Can be loaded directly into scanpy (`sc.read_h5ad(...)`) for "
        "further analysis outside this app."
    )
    if st.button("Prepare Download", key="sc_downstream_prepare_h5ad_btn"):
        with st.spinner("Preparing file..."):
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
                tmp_path = tmp.name
            adata.write_h5ad(tmp_path)
            with open(tmp_path, "rb") as f:
                st.session_state["sc_downstream_h5ad_bytes"] = f.read()
            os.remove(tmp_path)

    h5ad_bytes = st.session_state.get("sc_downstream_h5ad_bytes")
    if h5ad_bytes:
        st.download_button(
            "⬇️ Download AnnData (.h5ad)", data=h5ad_bytes,
            file_name=f"{project}_downstream.h5ad", mime="application/octet-stream",
            key="sc_downstream_h5ad_download_btn",
        )


# ---------------------------------------------------------------------------
# Step 5: Batch Correction
# ---------------------------------------------------------------------------

def _render_batch_correction_step(project, adata):
    st.header("Step 5: Batch Correction / Integration")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Adjusts each sample's PCA coordinates to remove "
            "systematic, sample-specific technical variation (e.g. differences "
            "in reagent lot, processing day, or sequencing run) WITHOUT erasing "
            "genuine biological differences between conditions/donors -- "
            "Harmony, the method offered here, works by iteratively nudging "
            "cells from different samples that already look similar to sit "
            "closer together in PCA space.\n\n"
            "**Why it matters:** If Step 4's PCA scatter plot showed samples "
            "separating into distinct clouds even for what should be similar "
            "cell types, clustering (Step 6) would otherwise group cells "
            "primarily by WHICH SAMPLE they came from, rather than by their "
            "actual cell type/state -- effectively discovering 'batches' "
            "instead of biology. This is only a real concern with 2+ samples; "
            "with a single sample there is nothing to correct for at all.\n\n"
            "**How to interpret the result:** There's no new plot here "
            "specifically, but it's worth revisiting Step 4's own PCA scatter "
            "(or Step 7's embedding scatter, once computed) afterward, colored "
            "by sample -- if correction worked, samples should now visibly mix "
            "together within each real cell population rather than forming "
            "separate clouds. If samples STILL separate cleanly after running "
            "Harmony, that can also be a sign of a genuine biological "
            "difference between conditions (e.g. a cell type present in one "
            "condition but absent in another) rather than a technical batch "
            "effect -- correction is not designed to erase real biology, only "
            "technical noise."
        )

    n_samples_in_data = adata.obs["sample"].nunique() if "sample" in adata.obs.columns else 1
    method_keys = list(dsm.BATCH_CORRECTION_METHOD_OPTIONS.keys())
    default_method = _recipe_default(
        project, "batch_correction", "method",
        dsm.DEFAULT_BATCH_CORRECTION_METHOD if n_samples_in_data > 1 else "none",
    )
    if default_method not in method_keys:
        default_method = "none"

    method = st.radio(
        "Batch correction method:", method_keys,
        format_func=lambda k: dsm.BATCH_CORRECTION_METHOD_OPTIONS[k]["label"],
        index=method_keys.index(default_method), key="sc_downstream_batch_method",
    )
    st.caption(dsm.BATCH_CORRECTION_METHOD_OPTIONS[method]["explanation"])

    if n_samples_in_data == 1 and method != "none":
        st.info(
            "ℹ️ Only one sample is present in this combined dataset -- there is "
            "nothing to correct for. Consider selecting 'Skip' unless you have "
            "a specific reason to run Harmony on a single sample."
        )

    if method == "harmony" and not dsm.harmony_available():
        st.error(
            "⚠️ The `harmonypy` package is not installed in this environment -- "
            "select 'Skip' to continue, or install harmonypy (see the "
            "⚙️ Setup & Deployment page)."
        )
        return False, None

    current_params = {"method": method, "batch_key": "sample" if method != "none" else None}
    is_current = scpm.check_downstream_step_current(project, "batch_correction", current_params)
    already_done = (method == "none") or ("X_pca_harmony" in adata.obsm)
    output_basis = "X_pca_harmony" if method == "harmony" else "X_pca"

    if method == "none":
        if not is_current:
            scpm.save_downstream_step_recipe(project, "batch_correction", current_params)
        st.success("✅ Batch correction skipped -- downstream steps will use the raw PCA embedding (`X_pca`).")
        return True, output_basis

    if is_current and already_done:
        st.success("✅ Harmony batch correction already applied with these exact settings.")
    elif already_done and not is_current:
        st.warning(
            "⚠️ Your current batch-correction settings differ from what's "
            "cached -- re-run below. (This will also invalidate every later "
            "step already computed.)"
        )

    run_label = "🔄 Re-run Harmony" if already_done and "X_pca_harmony" in adata.obsm else "▶️ Run Harmony"
    if st.button(run_label, key="sc_downstream_batch_btn", type="primary"):
        if "X_pca" not in adata.obsm:
            st.error("⚠️ Run PCA (Step 4) first.")
            return is_current and already_done, output_basis
        with st.spinner("Running Harmony batch correction..."):
            try:
                dsm.run_batch_correction(adata, method="harmony", batch_key="sample")
            except (RuntimeError, ValueError) as e:
                st.error(f"⚠️ {e}")
                return False, output_basis
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "batch_correction", current_params)
        st.success("✅ Harmony batch correction complete.")
        st.rerun()

    return (is_current and "X_pca_harmony" in adata.obsm), output_basis


# ---------------------------------------------------------------------------
# Step 6: Clustering (neighbors + leiden/louvain)
# ---------------------------------------------------------------------------

def _render_clustering_step(project, adata, use_rep):
    st.header("Step 6: Clustering")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** First builds a 'nearest-neighbor graph' -- for "
            "each cell, it finds the N most similar OTHER cells (in PCA space) "
            "and connects them. Then a community-detection algorithm (Leiden "
            "or Louvain) partitions that graph into discrete groups -- cells "
            "that are densely interconnected with each other, but only "
            "sparsely connected to cells outside the group, become one "
            "cluster.\n\n"
            "**Why it matters:** This is what actually turns a continuous "
            "cloud of individual cells into discrete, nameable groups -- the "
            "clusters produced here are what Step 7's embedding will visualize, "
            "and what Step 8's cell-type annotation assigns real biological "
            "labels to. Leiden is the more modern, generally recommended "
            "default; Louvain is offered as a well-established alternative, "
            "though note the caption below about resolution having no effect "
            "on it in this app's specific implementation.\n\n"
            "**How to interpret the result:** More clusters is not "
            "automatically 'better' -- resolution directly controls this "
            "trade-off (higher = more, smaller clusters; lower = fewer, "
            "larger ones), and there's no single universally correct value. "
            "A reasonable approach is to start around the default, look at "
            "the cluster size bar chart below (a cluster with only a handful "
            "of cells may not be biologically meaningful, or may be a doublet "
            "population Phase 2 QC didn't fully catch), and cross-check "
            "against Step 8's marker gene exploration -- if two clusters share "
            "essentially the same top marker genes, that's a sign resolution "
            "may be too high for this dataset."
        )

    default_n_neighbors = _recipe_default(project, "neighbors", "n_neighbors", dsm.DEFAULT_N_NEIGHBORS)
    n_neighbors = st.slider(
        "Number of neighbors:", min_value=5, max_value=100,
        value=int(default_n_neighbors), key="sc_downstream_n_neighbors",
        help="Used to build the nearest-neighbor graph that BOTH clustering and UMAP are computed from.",
    )

    method_keys = list(dsm.CLUSTERING_METHOD_OPTIONS.keys())
    default_method = _recipe_default(project, "clustering", "method", dsm.DEFAULT_CLUSTERING_METHOD)
    method = st.radio(
        "Clustering method:", method_keys,
        format_func=lambda k: dsm.CLUSTERING_METHOD_OPTIONS[k]["label"],
        index=method_keys.index(default_method) if default_method in method_keys else 0,
        key="sc_downstream_cluster_method",
    )
    st.caption(dsm.CLUSTERING_METHOD_OPTIONS[method]["explanation"])

    default_resolution = _recipe_default(project, "clustering", "resolution", dsm.DEFAULT_CLUSTERING_RESOLUTION)
    resolution = st.slider(
        "Resolution (higher = more, smaller clusters):", min_value=0.1, max_value=3.0,
        value=float(default_resolution), step=0.1, key="sc_downstream_cluster_resolution",
    )
    if method == "louvain":
        st.caption(
            "⚠️ Resolution has **no effect** on Louvain when run through this app's "
            "`igraph` backend (confirmed via direct testing) -- it is still shown "
            "above for consistency, but changing it will not change the result for "
            "this method."
        )

    cluster_key = method  # matches dsm.run_clustering's own key_added default

    neighbors_params = {"n_neighbors": n_neighbors, "use_rep": use_rep}
    clustering_params = {"method": method, "resolution": resolution}
    neighbors_current = scpm.check_downstream_step_current(project, "neighbors", neighbors_params)
    clustering_current = scpm.check_downstream_step_current(project, "clustering", clustering_params)
    already_has_neighbors = "neighbors" in adata.uns
    already_has_cluster = cluster_key in adata.obs.columns

    fully_current = neighbors_current and clustering_current and already_has_neighbors and already_has_cluster

    if fully_current:
        n_clusters_found = adata.obs[cluster_key].nunique()
        st.success(f"✅ Already clustered with these exact settings -- {n_clusters_found} cluster(s) found.")
    elif (already_has_neighbors or already_has_cluster) and not fully_current:
        st.warning(
            "⚠️ Your current neighbor/clustering settings differ from what's "
            "cached -- re-run below. (This will also invalidate embeddings and "
            "every later step already computed.)"
        )

    run_label = "🔄 Re-run Clustering" if already_has_cluster else "▶️ Compute Neighbors & Cluster"
    if st.button(run_label, key="sc_downstream_cluster_btn", type="primary"):
        if use_rep not in adata.obsm:
            st.error(f"⚠️ Run Step 4 (PCA)" + (" and Step 5 (Batch Correction)" if use_rep != "X_pca" else "") + " first.")
            return False
        with st.spinner("Computing neighbor graph..."):
            dsm.compute_neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
        scpm.save_downstream_step_recipe(project, "neighbors", neighbors_params)
        with st.spinner(f"Clustering ({method})..."):
            adata, resolution_had_effect = dsm.run_clustering(adata, method=method, resolution=resolution)
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "clustering", clustering_params)
        if not resolution_had_effect:
            st.info("ℹ️ Note: resolution had no effect on this Louvain run, as expected (see caption above).")
        st.success(f"✅ Clustering complete -- {adata.obs[cluster_key].nunique()} cluster(s) found.")
        st.rerun()

    if already_has_cluster and fully_current:
        with st.expander("📊 Cluster size summary", expanded=True):
            summary_df = dsm.get_cluster_summary(adata, cluster_key)
            if summary_df is not None:
                # Cluster labels are ALREADY in true numeric order here
                # (get_cluster_summary() itself now applies
                # dsm.natural_sort_unique() -- see that function's own
                # docstring) -- pass an explicit category_orders so the
                # bar chart's own x-axis respects this order too, rather
                # than Plotly silently re-deriving its own order.
                fig = go.Figure(go.Bar(x=summary_df["cluster"], y=summary_df["n_cells"], name="n_cells"))
                fig.update_layout(
                    height=400, margin=dict(l=10, r=10, t=30, b=10),
                    xaxis=dict(type="category", categoryorder="array", categoryarray=summary_df["cluster"].tolist()),
                )
                style = scw._render_plot_style_controls(
                    "sc_downstream_cluster_summary",
                    single_color_default={"n_cells": "#636EFA"}, show_legend_controls=False,
                )
                scw._apply_plot_style(
                    fig, style, default_title="Cluster sizes",
                    default_x_label="Cluster", default_y_label="Number of cells",
                )
                scw._render_plotly_chart(fig)
                st.caption(
                    "💡 **How to read this:** clusters should generally be a "
                    "reasonable fraction of your total cell count -- a cluster "
                    "with only a handful of cells relative to the others may "
                    "represent a rare-but-real cell population, OR a technical "
                    "artifact (e.g. a doublet population Phase 2 QC didn't "
                    "fully remove), and is worth checking against Step 8's "
                    "marker genes before treating it as a distinct cell type."
                )
                scw._render_pdf_export(fig, "sc_downstream_cluster_summary", "cluster_sizes_bar")
                scw._render_csv_download(
                    summary_df, "cluster_sizes", "sc_downstream_cluster_summary",
                    expander_label="⬇️ Download cluster size data (.csv)",
                )
        # --- NEW: cluster x sample breakdown (2026-09-09) -- see this
        # function's own docstring for the full rationale.
        _render_cluster_sample_breakdown(adata, cluster_key)

    return fully_current


DEFAULT_SAMPLE_DOMINANCE_THRESHOLD = 0.7


def _render_cluster_sample_breakdown(adata, cluster_key, sample_key="sample",
                                      dominance_threshold=DEFAULT_SAMPLE_DOMINANCE_THRESHOLD):
    """
    Cross-tabulates each cluster against the 'sample' column, to help
    distinguish a genuinely rare-but-real cell population from a
    cluster that's actually a batch/technical artifact -- i.e. a
    cluster composed almost entirely of cells from just ONE sample,
    rather than being spread across many samples (and, ideally, across
    multiple experimental conditions).

    dominance_threshold: if a single sample accounts for this fraction
        (or more) of a cluster's total cells, that cluster is flagged
        for review. 0.7 (70%) is a reasonable starting heuristic, not a
        hard statistical cutoff -- a project with very few samples
        overall, or a genuinely rare cell type expected to concentrate
        in only one condition, may reasonably trip this flag without
        it indicating a real problem. Always cross-check a flagged
        cluster's own marker genes (Step 8a) before concluding it's a
        technical artifact rather than real, condition-specific
        biology.

    Renders nothing (returns immediately) if sample_key isn't present
    in adata.obs at all (e.g. a single-sample project with no
    meaningful cross-sample comparison to make here).
    """
    if sample_key not in adata.obs.columns or cluster_key not in adata.obs.columns:
        return

    n_samples_in_data = adata.obs[sample_key].nunique()
    if n_samples_in_data <= 1:
        return  # nothing to cross-tab against with only one sample

    with st.expander("📊 Cluster x sample breakdown (checking for batch artifacts)", expanded=True):
        st.caption(
            "Each cluster's cells should generally be MIXED across multiple samples "
            "(and ideally across conditions) -- a cluster made up almost entirely of "
            "ONE sample's cells is a signal worth double-checking before treating it "
            "as a real, reproducible cell type: it could be a residual batch effect, "
            "or a sample-specific technical artifact (e.g. a leftover doublet "
            "population) rather than genuine biology."
        )

        counts_df = pd.crosstab(adata.obs[cluster_key].astype(str), adata.obs[sample_key].astype(str))
        ordered_clusters = dsm.natural_sort_unique(adata.obs[cluster_key].astype(str))
        counts_df = counts_df.reindex(ordered_clusters)
        proportions_df = counts_df.div(counts_df.sum(axis=1), axis=0)

        flagged = []
        for cluster in proportions_df.index:
            row = proportions_df.loc[cluster]
            max_sample = row.idxmax()
            max_prop = row.max()
            if max_prop >= dominance_threshold:
                flagged.append({
                    "cluster": cluster,
                    "n_cells": int(counts_df.loc[cluster].sum()),
                    "dominant_sample": max_sample,
                    "fraction_from_dominant_sample": round(float(max_prop), 3),
                })

        if flagged:
            flagged_clusters_str = ", ".join(f"{f['cluster']}" for f in flagged)
            st.warning(
                f"⚠️ {len(flagged)} cluster(s) are dominated (\u2265{dominance_threshold*100:.0f}%) by a "
                f"single sample: **{flagged_clusters_str}**. Review these against Step 8a's marker "
                f"genes before assigning them a confident cell-type label -- see the table below."
            )
            st.dataframe(pd.DataFrame(flagged), use_container_width=True, hide_index=True)
        else:
            st.success(
                f"✅ No cluster is dominated by a single sample (using a "
                f"{dominance_threshold*100:.0f}% threshold) -- every cluster's cells are "
                f"reasonably mixed across multiple samples."
            )

        # Stacked bar plot: proportion of each cluster contributed by
        # each sample -- lets a user visually scan every cluster at
        # once, not just the ones that cross the numeric threshold above.
        long_df = proportions_df.reset_index().melt(
            id_vars=cluster_key if cluster_key in proportions_df.reset_index().columns else "index",
            var_name=sample_key, value_name="proportion",
        )
        long_df = long_df.rename(columns={long_df.columns[0]: "cluster"})
        sample_values = dsm.natural_sort_unique(long_df[sample_key].astype(str))
        fig = px.bar(
            long_df, x="cluster", y="proportion", color=sample_key,
            category_orders={"cluster": ordered_clusters, sample_key: sample_values},
        )
        fig.update_layout(
            height=450, barmode="stack", margin=dict(l=10, r=10, t=30, b=10),
            xaxis=dict(type="category", categoryorder="array", categoryarray=ordered_clusters),
        )
        scw._render_plotly_chart(fig)
        st.caption(
            "💡 **How to read this:** each bar is one cluster, and each colored segment is "
            "the fraction of that cluster's cells contributed by one sample. A bar that's "
            "almost entirely ONE color is the same signal flagged numerically above -- "
            "worth a closer look before finalizing that cluster's cell-type label."
        )
        scw._render_pdf_export(fig, "sc_downstream_cluster_sample_breakdown", "cluster_sample_breakdown")
        scw._render_csv_download(
            counts_df.reset_index().rename(columns={"index": "cluster"}),
            "cluster_sample_counts", "sc_downstream_cluster_sample_breakdown",
            expander_label="⬇️ Download cluster x sample counts (.csv)",
        )


# ---------------------------------------------------------------------------
# Step 7: Embeddings (UMAP / t-SNE)
# ---------------------------------------------------------------------------

def _render_embedding_step(project, adata, use_rep, cluster_key):
    st.header("Step 7: Embedding (UMAP / t-SNE)")
    with st.expander("ℹ️ What is this step, and why does it matter?", expanded=False):
        st.markdown(
            "**What it does:** Further compresses the PCA embedding (already "
            "10s of dimensions) down to just 2 or 3 dimensions, specifically "
            "for VISUALIZATION -- UMAP and t-SNE are both designed to preserve "
            "local neighborhood structure (cells similar to each other in PCA "
            "space stay close together after this compression) at the cost of "
            "distorting large-scale distances.\n\n"
            "**Why it matters:** This is purely a visualization aid -- "
            "clustering (Step 6) is already complete and does NOT depend on "
            "this step at all, so it's safe to skip straight to Step 8 if you "
            "don't need a plot. Where it IS genuinely useful: visually "
            "confirming that your clusters look like sensible, separated "
            "groups, and later, for Step 8's 'feature plot' (visualizing one "
            "gene's expression across this same 2D/3D layout).\n\n"
            "**How to interpret the result:** Distances BETWEEN well-separated "
            "clusters on a UMAP/t-SNE plot are NOT meaningful -- two clusters "
            "drawn far apart are not necessarily more biologically different "
            "than two clusters drawn close together; only the LOCAL "
            "neighborhood structure (which cells are near which other cells) "
            "is preserved by design. Avoid reading anything into the overall "
            "shape, size, or relative positioning of clusters beyond that."
        )

    method_keys = list(dsm.EMBEDDING_METHOD_OPTIONS.keys())
    default_method = _recipe_default(project, "embedding", "method", dsm.DEFAULT_EMBEDDING_METHOD)
    method = st.radio(
        "Embedding method:", method_keys,
        format_func=lambda k: dsm.EMBEDDING_METHOD_OPTIONS[k]["label"],
        index=method_keys.index(default_method) if default_method in method_keys else 0,
        key="sc_downstream_embed_method",
    )
    st.caption(dsm.EMBEDDING_METHOD_OPTIONS[method]["explanation"])

    default_n_components = _recipe_default(project, "embedding", "n_components", 2)
    n_components = st.radio(
        "Dimensions:", [2, 3], index=[2, 3].index(default_n_components) if default_n_components in (2, 3) else 0,
        format_func=lambda n: f"{n}D", key="sc_downstream_embed_n_components", horizontal=True,
    )

    embed_key = "X_umap" if method == "umap" else "X_tsne"
    current_params = {"method": method, "use_rep": use_rep, "n_components": n_components}
    is_current = scpm.check_downstream_step_current(project, "embedding", current_params)
    already_done = embed_key in adata.obsm

    if is_current and already_done:
        st.success(f"✅ Already computed with these exact settings ({n_components}D {method.upper()}).")
    elif already_done and not is_current:
        st.warning("⚠️ Your current embedding settings differ from what's cached -- re-run below.")

    run_label = "🔄 Re-run Embedding" if already_done else "▶️ Compute Embedding"
    if st.button(run_label, key="sc_downstream_embed_btn", type="primary"):
        if method == "umap" and "neighbors" not in adata.uns:
            st.error("⚠️ Run Step 6 (Clustering, which computes the neighbor graph) first.")
            return False
        if method == "tsne" and use_rep not in adata.obsm:
            st.error(f"⚠️ Run Step 4/5 first -- '{use_rep}' embedding not found.")
            return False
        with st.spinner(f"Computing {n_components}D {method.upper()}..."):
            dsm.run_embedding(adata, method=method, use_rep=use_rep, n_components=n_components)
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "embedding", current_params)
        st.success("✅ Embedding complete.")
        st.rerun()

    if already_done and is_current:
        with st.expander("📊 Embedding scatter plot", expanded=True):
            color_options = ["sample"] + ([cluster_key] if cluster_key and cluster_key in adata.obs.columns else [])
            color_by = st.selectbox(
                "Color by:", options=color_options, key="sc_downstream_embed_color_by",
            )
            embed_df = dsm.get_embedding_coordinates(adata, method=method, color_by_columns=[color_by])
            if embed_df is not None and "Dim1" in embed_df.columns and "Dim2" in embed_df.columns:
                # --- Cluster/group label natural-sort ordering fix
                # (2026-08-25) -- see this module's own docstring for
                # the full rationale. Previously used
                # sorted(embed_df[color_by].astype(str).unique()) (a
                # plain lexicographic sort) to build group_values for
                # the color-picker widgets, but NEVER passed an explicit
                # category_orders to px.scatter/scatter_3d -- so Plotly's
                # own "first appearance in the data" default legend
                # order took over instead, producing the confirmed real
                # scrambled legend order. Both are fixed here: natural
                # sort for group_values, AND an explicit category_orders
                # argument telling Plotly the exact intended order.
                group_values = dsm.natural_sort_unique(embed_df[color_by].astype(str)) if color_by in embed_df.columns else None
                category_orders = {color_by: group_values} if group_values else None
                if n_components == 3 and "Dim3" in embed_df.columns:
                    fig = px.scatter_3d(
                        embed_df, x="Dim1", y="Dim2", z="Dim3", color=color_by, opacity=0.6,
                        category_orders=category_orders,
                    )
                else:
                    fig = px.scatter(
                        embed_df, x="Dim1", y="Dim2", color=color_by, opacity=0.6,
                        category_orders=category_orders,
                    )
                fig.update_layout(height=550, margin=dict(l=10, r=10, t=30, b=10))
                style = scw._render_plot_style_controls(
                    "sc_downstream_embed_scatter", group_values=group_values,
                    show_legend_controls=bool(group_values),
                )
                scw._apply_plot_style(
                    fig, style, default_title=f"{method.upper()} ({n_components}D)",
                    default_x_label="Dim1", default_y_label="Dim2",
                )
                scw._render_plotly_chart(fig)
                scw._render_pdf_export(fig, "sc_downstream_embed_scatter", f"{method}_embedding_plot")
                scw._render_csv_download(
                    embed_df, f"{method}_coordinates", "sc_downstream_embed_scatter",
                    expander_label="⬇️ Download embedding coordinates (.csv)",
                )

    return is_current and already_done


# ---------------------------------------------------------------------------
# Step 8: Cell-type Annotation
# ---------------------------------------------------------------------------

def _render_cluster_marker_explorer(project, adata, cluster_key):
    st.subheader("8a. Explore Cluster Marker Genes")
    st.caption(
        "Before assigning cell-type labels, it helps to see which genes most "
        "distinguish each cluster from the rest of the dataset. This uses "
        "cell-level statistics appropriately -- to compare CLUSTERS within "
        "this dataset, not to compare experimental CONDITIONS (that "
        "cross-sample comparison belongs in Step 9's pseudobulk bridge)."
    )
    with st.expander("ℹ️ How to interpret the marker gene table", expanded=False):
        st.markdown(
            "Each row is one gene that ranks highly for distinguishing ONE "
            "cluster from every other cluster combined. A HIGH `score` "
            "(the ranking statistic itself) combined with a LOW `pvalue_adj` "
            "and a positive `log2FoldChange` means that gene is both "
            "statistically robust AND meaningfully upregulated in that "
            "specific cluster -- these are your best candidates for a "
            "recognizable, textbook marker gene (e.g. CD3 for T cells, MS4A1 "
            "for B cells). A cluster whose top markers are all weak/ambiguous "
            "(low scores, or genes with no obvious biological role) can "
            "sometimes indicate that cluster doesn't represent a truly "
            "distinct cell population -- worth cross-checking against the "
            "manual marker scoring and CellTypist sections below before "
            "assigning it a confident final label."
        )

    method_options = {
        "wilcoxon": "Wilcoxon rank-sum (recommended)",
        "t-test": "t-test",
        "t-test_overestim_var": "t-test (overestimated variance)",
        "logreg": "Logistic regression",
    }
    default_n_genes = _recipe_default(project, "cluster_markers", "n_genes", 25)
    default_method = _recipe_default(project, "cluster_markers", "method", "wilcoxon")

    col1, col2 = st.columns(2)
    with col1:
        n_genes = st.slider(
            "Top genes per cluster:", min_value=5, max_value=100,
            value=int(default_n_genes), step=5, key="sc_downstream_marker_n_genes",
        )
    with col2:
        method = st.selectbox(
            "Ranking method:", options=list(method_options.keys()),
            format_func=lambda k: method_options[k],
            index=list(method_options.keys()).index(default_method) if default_method in method_options else 0,
            key="sc_downstream_marker_method",
        )

    if st.button("🔎 Find Cluster Marker Genes", key="sc_downstream_find_markers_btn", type="primary"):
        with st.spinner("Ranking marker genes per cluster..."):
            markers_df = dsm.find_cluster_markers(adata, groupby=cluster_key, method=method, n_genes=n_genes)
        st.session_state["sc_downstream_cluster_markers_df"] = markers_df
        st.success(f"✅ Found top {n_genes} marker genes for each of {adata.obs[cluster_key].nunique()} clusters.")

    markers_df = st.session_state.get("sc_downstream_cluster_markers_df")
    if markers_df is not None:
        with st.expander("📋 Marker gene table", expanded=True):
            # --- Cluster/group label natural-sort ordering fix
            # (2026-08-25) -- see this module's own docstring for the
            # full rationale. This dropdown previously used
            # sorted(markers_df["cluster"].unique(), key=str), a plain
            # LEXICOGRAPHIC string sort that put "10" right after "1"
            # and before "2" -- confirmed directly matching a real
            # reported screenshot showing exactly this order. Now uses
            # dsm.natural_sort_unique() for true numeric order instead.
            cluster_options = dsm.natural_sort_unique(markers_df["cluster"])
            chosen_cluster = st.selectbox(
                "View top markers for cluster:", options=["All"] + list(cluster_options),
                key="sc_downstream_marker_table_cluster",
            )
            display_df = markers_df if chosen_cluster == "All" else markers_df[markers_df["cluster"] == chosen_cluster]
                        # --- Scientific-notation p-value display fix (2026-08-25) ---
            # A real, confirmed issue: pvalue/pvalue_adj showed as "0"
            # for every row -- partly because Streamlit's own default
            # numeric display rounds away small floats, and partly
            # because scanpy's Wilcoxon test genuinely UNDERFLOWS to a
            # true 0.0 in float64 for large cell counts (no formatting
            # can recover a value that's already lost). Both cases are
            # now shown honestly: a real small p-value displays in true
            # scientific notation, and a genuinely underflowed value is
            # explicitly labeled "< 1e-308" (float64's own smallest
            # representable magnitude) rather than showing a misleading
            # bare "0" or a fabricated tiny number. The RAW, full-
            # precision numeric values are still what's written to the
            # CSV download below -- only this on-screen table is
            # reformatted for readability.
            display_df_formatted = display_df.copy()
            for pcol in ("pvalue", "pvalue_adj"):
                if pcol in display_df_formatted.columns:
                    display_df_formatted[pcol] = display_df_formatted[pcol].apply(_format_scientific_pvalue)
            st.dataframe(display_df_formatted, use_container_width=True, hide_index=True)
            st.caption(
                "💡 p-values are shown in scientific notation. **'< 1e-308'** means "
                "the true p-value was so small it exceeded float64's own precision "
                "floor (a real statistical result, not a display bug) -- this is "
                "common with cluster-marker tests on large cell counts; see Step 8a's "
                "own explanation above for why the p-value's exact magnitude here "
                "matters less than `score` and `log2FoldChange` for judging a marker's "
                "quality."
            )
            scw._render_csv_download(
                markers_df, "cluster_marker_genes", "sc_downstream_marker_table",
                expander_label="⬇️ Download full marker gene table (.csv)",
            )
    return markers_df


def _render_manual_marker_scoring(project, adata, cluster_key):
    st.subheader("8b. Manual Marker-Gene Scoring (required baseline)")
    st.caption(dsm.ANNOTATION_METHOD_OPTIONS["manual_markers"]["explanation"])
    with st.expander("ℹ️ What this does, and how to interpret the mean-score table", expanded=False):
        st.markdown(
            "**What it does:** For each cell-type panel you define below (a "
            "name plus a list of known marker genes for that cell type), every "
            "cell in the dataset gets a single combined 'score' summarizing how "
            "strongly it expresses that WHOLE panel together -- not just any "
            "one gene in isolation.\n\n"
            "**How to interpret the mean-score table:** Each row is a cluster, "
            "each column is one of your defined cell-type panels -- the value "
            "is that cluster's AVERAGE score for that panel. A cluster with a "
            "clearly higher score for one cell type than all the others is a "
            "strong candidate for that label. A cluster with similarly "
            "moderate scores across MULTIPLE cell types is more ambiguous -- "
            "either it's a genuinely mixed/transitional population, your "
            "marker panels overlap biologically (e.g. two very closely related "
            "immune subsets), or clustering resolution (Step 6) may need "
            "adjusting to better separate them."
        )

    if "sc_downstream_marker_sets" not in st.session_state:
        st.session_state["sc_downstream_marker_sets"] = dict(adata.uns.get("marker_gene_sets", {}))

    marker_sets = st.session_state["sc_downstream_marker_sets"]

    with st.expander("✏️ Define marker gene panels", expanded=not bool(marker_sets)):
        # --- Input-clearing fix (2026-08-25) -- see this file's own
        # module docstring. Confirmed via direct testing that
        # st.form(clear_on_submit=True) does NOT reliably clear these
        # widgets in this app's real Streamlit version (1.53.1) -- so a
        # more robust "key-cycling" approach is used instead: an
        # incrementing counter is baked into each widget's own key, so
        # after a successful add, the NEXT render uses fresh widget
        # keys Streamlit has never seen before, which always render
        # blank -- confirmed working via direct testing.
        if "sc_downstream_panel_form_counter" not in st.session_state:
            st.session_state["sc_downstream_panel_form_counter"] = 0
        _form_counter = st.session_state["sc_downstream_panel_form_counter"]

        col1, col2 = st.columns([1, 2])
        with col1:
            new_cell_type = st.text_input(
                "Cell type name:", key=f"sc_downstream_new_celltype_name_{_form_counter}",
            )
        with col2:
            new_genes = st.text_input(
                "Marker genes (comma-separated, symbols or IDs):",
                key=f"sc_downstream_new_celltype_genes_{_form_counter}",
            )
        if st.button("➕ Add / Update Panel", key=f"sc_downstream_add_marker_panel_btn_{_form_counter}"):
            if not new_cell_type.strip() or not new_genes.strip():
                st.error("Enter both a cell type name and at least one marker gene.")
            else:
                genes_list = [g.strip() for g in new_genes.split(",") if g.strip()]
                marker_sets[new_cell_type.strip()] = genes_list
                st.session_state["sc_downstream_marker_sets"] = marker_sets
                st.session_state["sc_downstream_panel_form_counter"] += 1
                st.rerun()

        # --- File upload fix (2026-08-25) -- lets a user provide an
        # entire marker panel table at once (.csv/.txt/.xlsx), instead
        # of typing each cell type's genes one at a time. Supports BOTH
        # a "wide" layout (one row per cell type, genes comma/semicolon/
        # pipe-separated within one cell) and a "long" layout (one gene
        # per row, the SAME cell_type repeated across multiple rows) --
        # see sc_marker_panel_io.py's own module docstring for the full
        # parsing rules. Uploaded panels are MERGED into (not replacing)
        # any panels already defined above -- a cell type present in
        # both is overwritten by the file's own version, matching this
        # section's own "Add/Update" semantics.
        st.markdown("**Or upload a marker panel file:**")
        uploaded_panel_file = st.file_uploader(
            "Upload a .csv, .txt, or .xlsx file:", type=["csv", "txt", "xlsx"],
            key=f"sc_downstream_panel_upload_{_form_counter}",
            help=(
                "Expected format: a column identifying the cell type (e.g. "
                "'cell_type') and a column of gene(s) (e.g. 'genes') -- either "
                "one row per cell type with multiple genes separated by commas/"
                "semicolons/pipes in one cell, OR one gene per row with the same "
                "cell type repeated across multiple rows. Both layouts are "
                "auto-detected."
            ),
        )
        if uploaded_panel_file is not None:
            if st.button("➕ Add Panels From File", key=f"sc_downstream_add_panels_from_file_btn_{_form_counter}"):
                try:
                    uploaded_panels = mpio.load_marker_panel_file(
                        uploaded_panel_file.name, uploaded_panel_file.getvalue(),
                    )
                except ValueError as e:
                    st.error(f"⚠️ {e}")
                else:
                    marker_sets.update(uploaded_panels)
                    st.session_state["sc_downstream_marker_sets"] = marker_sets
                    st.session_state["sc_downstream_panel_form_counter"] += 1
                    st.success(f"✅ Added/updated {len(uploaded_panels)} panel(s) from '{uploaded_panel_file.name}'.")
                    st.rerun()

        if marker_sets:
            st.markdown("**Current panels:**")
            for cell_type, genes in list(marker_sets.items()):
                row_cols = st.columns([3, 1])
                row_cols[0].markdown(f"- **{cell_type}**: {', '.join(genes)}")
                if row_cols[1].button("🗑️ Remove", key=f"sc_downstream_remove_panel_{cell_type}"):
                    del marker_sets[cell_type]
                    st.session_state["sc_downstream_marker_sets"] = marker_sets
                    st.rerun()

    if not marker_sets:
        st.info("Add at least one marker gene panel above to run scoring.")
        return marker_sets

    layer_options = [l for l in ("lognorm", "pearson_residuals") if l in adata.layers]
    layer = layer_options[0] if layer_options else "lognorm"
    if len(layer_options) > 1:
        layer = st.selectbox(
            "Expression layer to score from:", options=layer_options,
            key="sc_downstream_marker_score_layer",
        )

    if st.button("🧮 Score Marker Panels", key="sc_downstream_score_markers_btn", type="primary"):
        with st.spinner("Scoring marker gene panels..."):
            match_report = dsm.score_marker_gene_sets(adata, marker_sets, layer=layer)
        adata.uns["marker_gene_sets"] = marker_sets
        _save_adata_state(project, adata)
        st.session_state["sc_downstream_marker_match_report"] = match_report
        st.success("✅ Marker panels scored.")
        st.rerun()

    match_report = st.session_state.get("sc_downstream_marker_match_report")
    if match_report:
        with st.expander("📋 Gene match report", expanded=False):
            report_rows = [
                {
                    "cell_type": ct, "genes_requested": r["n_genes_requested"],
                    "genes_found": r["n_genes_found"], "matched_genes": ", ".join(r["genes_found"]),
                }
                for ct, r in match_report.items()
            ]
            st.dataframe(pd.DataFrame(report_rows), use_container_width=True, hide_index=True)
            for ct, r in match_report.items():
                if r["n_genes_found"] == 0:
                    st.warning(f"⚠️ None of the requested genes for '{ct}' were found in this dataset.")
                elif r["n_genes_found"] < r["n_genes_requested"]:
                    st.info(f"ℹ️ Only {r['n_genes_found']} of {r['n_genes_requested']} requested genes for '{ct}' were found.")

    score_cols = [f"score_{ct}" for ct in marker_sets if f"score_{ct}" in adata.obs.columns]
    if score_cols and cluster_key in adata.obs.columns:
        with st.expander("📊 Mean marker score per cluster", expanded=True):
            mean_scores = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
            mean_scores.columns = [c.replace("score_", "") for c in mean_scores.columns]
            # Re-order rows into true numeric cluster order, matching
            # this file's other cluster-ordering fixes.
            ordered_index = dsm.natural_sort_unique(mean_scores.index.astype(str))
            mean_scores = mean_scores.reindex(ordered_index)
            st.dataframe(mean_scores.style.background_gradient(cmap="RdBu_r", axis=None), use_container_width=True)
            scw._render_csv_download(
                mean_scores.reset_index(), "mean_marker_scores_per_cluster", "sc_downstream_marker_scores",
                expander_label="⬇️ Download mean scores per cluster (.csv)",
            )
    return marker_sets


def _render_celltypist_step(project, adata, cluster_key):
    st.subheader("8c. CellTypist Automated First-Pass (optional)")
    st.caption(dsm.ANNOTATION_METHOD_OPTIONS["celltypist"]["explanation"])
    with st.expander("ℹ️ How to interpret CellTypist's suggestions", expanded=False):
        st.markdown(
            "CellTypist is a pre-trained classifier -- it was trained on "
            "OTHER published datasets, not yours, so treat its predictions as "
            "a helpful second opinion to CROSS-CHECK against your own manual "
            "marker scoring above, not as ground truth to accept "
            "automatically. It works best for well-studied tissue types "
            "(e.g. immune cells in blood) that closely match its training "
            "data, and can be less reliable for unusual tissues, disease "
            "states, or non-human organisms.\n\n"
            "**Majority voting** (the checkbox above) is recommended because "
            "it re-assigns EVERY cell in one of your own Step 6 clusters to "
            "that cluster's single most common CellTypist prediction -- "
            "smoothing out noisy, cell-by-cell disagreement into one clean "
            "suggestion per cluster, which is exactly what the breakdown table "
            "below shows. If a cluster's breakdown is split roughly evenly "
            "across several different CellTypist labels with no clear "
            "majority, that's a sign this cluster may not map cleanly onto "
            "any single cell type in CellTypist's own reference -- lean more "
            "heavily on your manual marker panels for that cluster instead."
        )

    if not dsm.celltypist_available():
        st.warning(
            "⚠️ The `celltypist` package is not installed in this environment -- "
            "skip this step, or install celltypist (see ⚙️ Setup & Deployment)."
        )
        return None

    available_models = dsm.get_celltypist_available_models()
    if not available_models:
        st.warning("⚠️ Could not retrieve CellTypist's model catalog.")
        return None

    default_model = _recipe_default(project, "celltypist", "model_name", "Immune_All_Low.pkl")
    model_options = list(available_models.keys())
    model_name = st.selectbox(
        "CellTypist model:", options=model_options,
        format_func=lambda m: f"{m} -- {available_models[m]}",
        index=model_options.index(default_model) if default_model in model_options else 0,
        key="sc_downstream_celltypist_model",
    )

    is_downloaded = dsm.celltypist_model_is_downloaded(model_name)
    if not is_downloaded:
        st.info(f"ℹ️ Model '{model_name}' has not been downloaded yet.")
        if st.button("⬇️ Download Model", key="sc_downstream_download_celltypist_btn"):
            with st.spinner(f"Downloading '{model_name}'..."):
                success, message = dsm.download_celltypist_model(model_name)
            if success:
                st.success(f"✅ {message}")
                st.rerun()
            else:
                st.error(f"⚠️ {message}")
        return None
    st.success(f"✅ Model '{model_name}' is ready to use.")

    default_majority_voting = _recipe_default(project, "celltypist", "majority_voting", True)
    majority_voting = st.checkbox(
        "Refine predictions with majority voting on this project's existing clusters (recommended)",
        value=default_majority_voting, key="sc_downstream_celltypist_majority_voting",
        help=(
            f"Uses this project's own '{cluster_key}' clusters (Step 6) as the "
            "over-clustering input, so CellTypist's predictions can be compared "
            "cluster-for-cluster against your manual marker scoring above."
        ),
    )

    current_params = {
        "model_name": model_name, "majority_voting": majority_voting,
        "cluster_key": cluster_key if majority_voting else None,
    }
    is_current = adata.uns.get("celltypist_params") == current_params
    already_done = "celltypist_predicted_label" in adata.obs.columns

    if is_current and already_done:
        st.success("✅ CellTypist annotation already run with these exact settings.")
    elif already_done and not is_current:
        st.warning("⚠️ Your current CellTypist settings differ from what's cached -- re-run below.")

    run_label = "🔄 Re-run CellTypist" if already_done else "▶️ Run CellTypist Annotation"
    if st.button(run_label, key="sc_downstream_run_celltypist_btn", type="primary"):
        with st.spinner("Running CellTypist..."):
            success, message = dsm.run_celltypist_annotation(
                adata, model_name=model_name, majority_voting=majority_voting,
                cluster_key=cluster_key if majority_voting else None,
            )
        if success:
            adata.uns["celltypist_params"] = current_params
            _save_adata_state(project, adata)
            st.success(f"✅ {message}")
            st.rerun()
        else:
            st.error(f"⚠️ {message}")

    if already_done and is_current and cluster_key in adata.obs.columns:
        with st.expander("📊 CellTypist label breakdown per cluster", expanded=True):
            label_col = (
                "celltypist_majority_voting" if "celltypist_majority_voting" in adata.obs.columns
                else "celltypist_predicted_label"
            )
            breakdown = pd.crosstab(adata.obs[cluster_key], adata.obs[label_col])
            # Re-order rows into true numeric cluster order.
            ordered_index = dsm.natural_sort_unique(breakdown.index.astype(str))
            breakdown = breakdown.reindex(ordered_index)
            st.dataframe(breakdown, use_container_width=True)
            scw._render_csv_download(
                breakdown.reset_index(), "celltypist_breakdown_per_cluster", "sc_downstream_celltypist_breakdown",
                expander_label="⬇️ Download breakdown (.csv)",
            )

    return is_current and already_done


def _render_annotation_visualizations(adata, cluster_key, marker_sets):
    st.subheader("8d. Marker Visualization")
    with st.expander("ℹ️ What each plot type shows, and when to use it", expanded=False):
        st.markdown(
            "**Dot Plot:** one row per group (cluster/cell type), one column "
            "per gene -- dot SIZE shows what fraction of cells in that group "
            "express the gene at all, dot COLOR shows the average expression "
            "level among cells that do. The single best plot for a quick "
            "overview of many genes across many groups at once -- look for "
            "genes that are both large AND darkly colored in exactly one "
            "group, and faint/small everywhere else.\n\n"
            "**Violin Plot:** the full expression DISTRIBUTION of one gene, "
            "split out per group. Useful when a dot plot's single average "
            "value hides something important -- e.g. a gene that's bimodal "
            "(some cells high, some near-zero) within what looks like one "
            "cluster can be a clue that cluster is actually a mix of two "
            "cell states.\n\n"
            "**Feature Plot:** one gene's expression overlaid directly on the "
            "Step 7 UMAP/t-SNE layout -- lets you see WHERE on the embedding "
            "a gene is expressed, which is especially useful for confirming "
            "that a marker gene's expression pattern actually lines up with a "
            "specific, visually distinct region/cluster rather than being "
            "scattered diffusely everywhere.\n\n"
            "**Heatmap:** many genes x many individual cells at once, grouped "
            "and separated by cluster/cell type -- the most detailed view, "
            "useful for checking that an entire marker PANEL (not just one "
            "gene) shows a clean, consistent block pattern per group."
        )

    color_by_options = [cluster_key]
    if "celltypist_predicted_label" in adata.obs.columns:
        color_by_options.append("celltypist_predicted_label")
    if "celltypist_majority_voting" in adata.obs.columns:
        color_by_options.append("celltypist_majority_voting")

    groupby = st.selectbox("Group cells by:", options=color_by_options, key="sc_downstream_annot_viz_groupby")

    default_genes = []
    if marker_sets:
        for genes in marker_sets.values():
            default_genes.extend(genes)
    default_genes = list(dict.fromkeys(default_genes))[:20]

    genes_input = st.text_area(
        "Genes to visualize (comma-separated, symbols or IDs):",
        value=", ".join(default_genes), key="sc_downstream_viz_genes_input",
    )
    genes = [g.strip() for g in genes_input.split(",") if g.strip()]
    if not genes:
        st.info("Enter at least one gene above to generate plots.")
        return

    layer_options = [l for l in ("lognorm", "pearson_residuals") if l in adata.layers]
    layer = layer_options[0] if layer_options else "lognorm"

    tab_dot, tab_violin, tab_feature, tab_heatmap = st.tabs(
        ["Dot Plot", "Violin Plot", "Feature Plot", "Heatmap"]
    )

    with tab_dot:
        dot_df = dsm.get_dotplot_data(adata, genes, groupby, layer=layer)
        if dot_df is not None and not dot_df.empty:
            fig = px.scatter(
                dot_df, x="gene_symbol", y="group", size="pct_expressing", color="mean_expression",
                color_continuous_scale="Reds", size_max=20,
            )
            fig.update_layout(height=max(300, 40 * dot_df["group"].nunique()), margin=dict(l=10, r=10, t=30, b=10))
            scw._render_plotly_chart(fig)
            scw._render_pdf_export(fig, "sc_downstream_dotplot", "marker_dotplot")
            scw._render_csv_download(
                dot_df, "dotplot_data", "sc_downstream_dotplot",
                expander_label="⬇️ Download dot plot data (.csv)",
            )
        else:
            st.info("None of the requested genes were found in this dataset.")

    with tab_violin:
        gene_for_violin = st.selectbox("Gene:", options=genes, key="sc_downstream_violin_gene")
        violin_df = dsm.get_violin_plot_data(adata, gene_for_violin, groupby, layer=layer)
        if violin_df is not None:
            fig = px.violin(violin_df, x="group", y="expression", box=True, points=False)
            fig.update_layout(height=450, margin=dict(l=10, r=10, t=30, b=10))
            scw._render_plotly_chart(fig)
            scw._render_pdf_export(fig, "sc_downstream_violin", f"{gene_for_violin}_violin")
            scw._render_csv_download(
                violin_df, f"{gene_for_violin}_violin_data", "sc_downstream_violin",
                expander_label="⬇️ Download violin data (.csv)",
            )
        else:
            st.info(f"Gene '{gene_for_violin}' not found in this dataset.")

    with tab_feature:
        gene_for_feature = st.selectbox("Gene:", options=genes, key="sc_downstream_feature_gene")
        embed_method_options = [m for m in ("umap", "tsne") if f"X_{m}" in adata.obsm]
        if not embed_method_options:
            st.info("Compute an embedding (Step 7) first.")
        else:
            embed_method = st.selectbox(
                "Embedding:", options=embed_method_options, key="sc_downstream_feature_embed_method",
            )
            feature_df = dsm.get_feature_plot_data(adata, gene_for_feature, embedding_method=embed_method, layer=layer)
            if feature_df is not None:
                fig = px.scatter(
                    feature_df, x="Dim1", y="Dim2", color="expression",
                    color_continuous_scale="Viridis", opacity=0.7,
                )
                fig.update_layout(height=500, margin=dict(l=10, r=10, t=30, b=10))
                scw._render_plotly_chart(fig)
                scw._render_pdf_export(fig, "sc_downstream_feature_plot", f"{gene_for_feature}_feature_plot")
                scw._render_csv_download(
                    feature_df, f"{gene_for_feature}_feature_data", "sc_downstream_feature_plot",
                    expander_label="⬇️ Download feature plot data (.csv)",
                )
            else:
                st.info(f"Gene '{gene_for_feature}' not found in this dataset.")

    with tab_heatmap:
        n_cells_per_group = st.slider(
            "Max cells per group (for readability):", min_value=10, max_value=200,
            value=50, step=10, key="sc_downstream_heatmap_n_cells",
        )
        z_df, ordered_groups = dsm.get_marker_heatmap_data(
            adata, genes, groupby, layer=layer, n_cells_per_group=n_cells_per_group,
        )
        if z_df is not None:
            fig = go.Figure(go.Heatmap(
                z=z_df.values, x=list(range(z_df.shape[1])), y=z_df.index.tolist(),
                colorscale="RdBu_r", zmid=0,
            ))
            fig.update_layout(
                height=max(400, 20 * len(z_df.index)), margin=dict(l=10, r=10, t=30, b=10),
                xaxis=dict(showticklabels=False, title=f"Cells (grouped by {groupby})"),
            )
            scw._render_plotly_chart(fig)
            scw._render_pdf_export(fig, "sc_downstream_marker_heatmap", "marker_heatmap")
        else:
            st.info("None of the requested genes were found in this dataset.")


def _render_final_celltype_assignment(project, adata, cluster_key):
    st.subheader("8e. Assign Final Cell-Type Labels")
    st.caption(
        "Review the marker scores / CellTypist suggestions above, then assign a "
        "final cell-type label to each cluster. This creates `adata.obs['cell_type']`, "
        "the column Step 9 (Pseudobulk -> DESeq2) and Step 10 (Compositional Analysis) "
        "will both use to group cells by cell type."
    )
    with st.expander("ℹ️ Why this step is a human judgment call, not an automated one", expanded=False):
        st.markdown(
            "This is deliberately the ONE point in Step 8 where you make the "
            "final call yourself, rather than the app choosing for you -- "
            "marker gene ranking (8a), manual scoring (8b), and CellTypist "
            "(8c) are all just EVIDENCE to weigh, since no single automated "
            "signal is reliably correct across every cluster in every "
            "dataset. It's entirely normal, and often correct, to give two "
            "different clusters the SAME cell-type label (e.g. 'CD4+ T cell' "
            "for two clusters that differ mainly in activation state rather "
            "than true cell identity) -- Step 10's compositional analysis "
            "will then correctly treat them as one combined population when "
            "comparing proportions across conditions.\n\n"
            "**A practical tip:** the `celltypist_suggestion` column is only "
            "a starting point pre-filled into `final_label` -- always confirm "
            "it against your own manual marker scores (8b) before saving, "
            "especially for any cluster where CellTypist's own breakdown "
            "table showed no clear majority."
        )

    clusters = dsm.natural_sort_unique(adata.obs[cluster_key].astype(str))
    saved_recipe = scpm.get_downstream_step_recipe(project, "annotation")
    saved_labels = dict(saved_recipe["params"].get("labels", {})) if saved_recipe else {}

    # --- New "top manual marker" column (2026-08-25) -- see
    # sc_marker_confidence.py's own module docstring for the full
    # rationale: this pulls Step 8b's own mean-score-per-cluster table
    # directly into THIS table, so a user never has to scroll back to
    # 8b to remember which panel scored highest for a given cluster --
    # the exact gap that caused a real, confirmed mislabeling before.
    mean_scores_by_cluster = {}
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if score_cols and cluster_key in adata.obs.columns:
        mean_scores_df = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
        mean_scores_df.columns = [c.replace("score_", "") for c in mean_scores_df.columns]
        for cluster_label, row in mean_scores_df.iterrows():
            sorted_scores = sorted(row.items(), key=lambda kv: -kv[1])
            mean_scores_by_cluster[str(cluster_label)] = sorted_scores

    rows = []
    for cluster in clusters:
        n_cells = int((adata.obs[cluster_key].astype(str) == cluster).sum())
        suggestion = ""
        if "celltypist_majority_voting" in adata.obs.columns:
            mode_vals = adata.obs.loc[adata.obs[cluster_key].astype(str) == cluster, "celltypist_majority_voting"].mode()
            suggestion = mode_vals.iloc[0] if not mode_vals.empty else ""
        elif "celltypist_predicted_label" in adata.obs.columns:
            mode_vals = adata.obs.loc[adata.obs[cluster_key].astype(str) == cluster, "celltypist_predicted_label"].mode()
            suggestion = mode_vals.iloc[0] if not mode_vals.empty else ""

        top_marker_name, confidence, reason = mconf.classify_marker_confidence(
            mean_scores_by_cluster.get(str(cluster), [])
        )
        top_marker_display = (
            f"{mconf.CONFIDENCE_ICONS[confidence]}: {top_marker_name}" if top_marker_name
            else mconf.CONFIDENCE_ICONS["No signal"]
        )

        rows.append({
            "cluster": str(cluster), "n_cells": n_cells,
            "celltypist_suggestion": suggestion,
            "top_manual_marker": top_marker_display,
            "final_label": saved_labels.get(str(cluster), top_marker_name or suggestion or str(cluster)),
        })
    edit_df = pd.DataFrame(rows)

    st.caption(
        "💡 **'Top manual marker'** pulls Step 8b's own mean-score table directly "
        "into this row, so you don't need to scroll back to compare -- 🟢 **High** "
        "means a clear, trustworthy standout; 🟡 **Potential** means a real but "
        "modest or closely-contested lead (worth double-checking); 🟠 **Low** means "
        "only a weak hint; ⚪ **No signal** means none of your defined panels scored "
        "positively for this cluster at all."
    )
    edited_df = st.data_editor(
        edit_df, use_container_width=True, hide_index=True, key="sc_downstream_final_label_editor",
        column_config={
            "cluster": st.column_config.TextColumn(disabled=True),
            "n_cells": st.column_config.NumberColumn(disabled=True),
            "celltypist_suggestion": st.column_config.TextColumn(disabled=True),
            "top_manual_marker": st.column_config.TextColumn("Top manual marker (from 8b)", disabled=True),
            "final_label": st.column_config.TextColumn("Final cell-type label"),
        },
     )
    # --- NEW: exclude specific cell-type label(s) from Step 9/10
    # (2026-09-09) -- see this file's own "excluded cell types" patch
    # notes for the full rationale. A stale label (e.g. renamed since
    # the last save) is automatically dropped from the default
    # selection below, rather than causing an error.
    final_label_options = sorted(edited_df["final_label"].astype(str).unique())
    default_excluded = [
        c for c in adata.uns.get(EXCLUDED_CELL_TYPES_UNS_KEY, []) if c in final_label_options
    ]
    excluded_cell_types = st.multiselect(
        "🚫 Exclude cell-type label(s) from Step 9 (Pseudobulk) & Step 10 (Compositional Analysis):",
        options=final_label_options, default=default_excluded,
        key="sc_downstream_excluded_cell_types",
        help=(
            "Cells with an excluded label are NOT removed from the dataset -- they "
            "remain fully visible everywhere else (plots, marker exploration, "
            "downloads). They are only filtered out of Step 9 and Step 10's own "
            "analyses. Use this for QC-artifact 'cell types' (e.g. ambient RNA/RBC "
            "contamination) that aren't real biological populations and would only "
            "add noise -- or a misleading result -- to a differential-expression or "
            "compositional comparison."
        ),
    )    
    current_label_map = dict(zip(edited_df["cluster"].astype(str), edited_df["final_label"].astype(str)))
    current_params = {
        "labels": current_label_map, "cluster_key": cluster_key,
        "excluded_cell_types": sorted(excluded_cell_types),   # <-- NEW
    }
    is_current = scpm.check_downstream_step_current(project, "annotation", current_params)
    already_done = "cell_type" in adata.obs.columns

    if already_done and is_current:
        st.success(f"✅ Final labels saved -- {adata.obs['cell_type'].nunique()} distinct cell type(s).")
    elif already_done and not is_current:
        st.warning("⚠️ Your edited labels above differ from what's saved -- save again to apply them.")

    if st.button("✅ Save Final Cell-Type Labels", key="sc_downstream_save_final_labels_btn", type="primary"):
        adata.obs["cell_type"] = adata.obs[cluster_key].astype(str).map(current_label_map)
        adata.uns[EXCLUDED_CELL_TYPES_UNS_KEY] = sorted(excluded_cell_types)   # <-- NEW
        _save_adata_state(project, adata)
        scpm.save_downstream_step_recipe(project, "annotation", current_params)
        st.success("✅ Final cell-type labels saved to `adata.obs['cell_type']`.")
        st.rerun()
    return already_done and is_current


def _render_celltype_annotation_step(project, adata, cluster_key):
    st.header("Step 8: Cell-type Annotation")
    st.markdown(
        "Combines cluster marker-gene exploration, manual marker-gene scoring "
        "(required baseline), and an optional CellTypist automated first pass -- "
        "review all three, then assign a final label per cluster below."
    )

    _render_cluster_marker_explorer(project, adata, cluster_key)
    st.markdown("---")
    marker_sets = _render_manual_marker_scoring(project, adata, cluster_key)
    st.markdown("---")
    _render_celltypist_step(project, adata, cluster_key)
    st.markdown("---")
    _render_annotation_visualizations(adata, cluster_key, marker_sets)
    st.markdown("---")
    annotation_ready = _render_final_celltype_assignment(project, adata, cluster_key)

    return annotation_ready


# ---------------------------------------------------------------------------
# Step 9: Pseudobulk -> DESeq2 Bridge
# ---------------------------------------------------------------------------

_TECHNICAL_OBS_COLUMN_PREFIXES = ("score_",)
_TECHNICAL_OBS_COLUMNS = {
    "leiden", "louvain", "celltypist_predicted_label", "celltypist_majority_voting",
    "adaptive_qc_fail", "predicted_doublet",
}


def _get_groupable_obs_columns(adata):
    "Every .obs column that could plausibly serve as a pseudobulk groupby dimension or compositional condition column -- excludes obviously technical/derived columns, but does not exclude 'sample' or cluster/cell_type columns (those ARE meaningful groupby choices)."
    cols = []
    for c in adata.obs.columns:
        if c in _TECHNICAL_OBS_COLUMNS:
            continue
        if any(c.startswith(p) for p in _TECHNICAL_OBS_COLUMN_PREFIXES):
            continue
        cols.append(c)
    return cols


def _render_pseudobulk_step(project, adata, cluster_key):
    st.header("Step 9: Pseudobulk -> DESeq2 Bridge")
    st.markdown(
        "Aggregates raw per-cell counts into pseudobulk 'samples' by summing "
        "within each group you choose below, producing a counts/metadata pair "
        "ready to hand directly to the Bulk RNA-Seq pipeline's own, already-"
        "tested DESeq2 workflow. Individual cells are **never** treated as "
        "independent replicates here -- see `aggregate_pseudobulk()`'s own "
        "docstring (Squair et al. 2021) for why that would produce severely "
        "inflated false-discovery rates."
    )
    with st.expander("ℹ️ Why aggregate cells at all, instead of testing per-cell directly?", expanded=False):
        st.markdown(
            "**The core problem this step solves:** cells from the SAME "
            "sample are not statistically independent of each other -- two "
            "cells from the same donor/mouse are more similar to each other "
            "than to cells from a different donor, purely because they share "
            "that donor's own biology, technical processing, etc. Treating "
            "thousands of individual cells as if they were thousands of "
            "independent replicates (a common mistake) makes a differential "
            "expression test wildly overconfident -- it can report a gene as "
            "'highly significant' when the real underlying difference is "
            "driven by just one or two samples, not a true, repeatable "
            "biological effect.\n\n"
            "**What 'pseudobulk' fixes:** by SUMMING all cells within one "
            "group (e.g. all cells from one sample, or one sample x cell-type "
            "combination) into a single pseudobulk 'sample', the true, "
            "correct unit of replication becomes the ACTUAL biological sample "
            "again -- exactly matching how a real bulk RNA-seq experiment "
            "would be analyzed, using DESeq2's own well-established, properly "
            "calibrated statistics.\n\n"
            "**Choosing a groupby:** `['sample']` alone gives one pseudobulk "
            "sample per original sample -- the simplest, most standard "
            "choice. Adding a cluster/cell-type column tests each cell type "
            "SEPARATELY (e.g. 'is gene X different in T cells specifically "
            "between treated vs. untreated samples?'), but requires enough "
            "cells of that type in EVERY sample to be statistically "
            "meaningful -- check the group-size preview table below before "
            "aggregating; any group with very few cells is a weak, "
            "unreliable pseudobulk sample even if it clears the minimum-cell "
            "threshold."
        )
     # --- NEW: respect Step 8e's excluded cell-type list (2026-09-09)
     # -- see this file's own "excluded cell types" patch notes for the
     # full rationale. Rebinding `adata` here is scoped to THIS
     # function call only -- confirmed safe since this function never
     # calls _save_adata_state(), so the caller's own combined
     # AnnData object is never affected.
    excluded_cell_types = list(adata.uns.get(EXCLUDED_CELL_TYPES_UNS_KEY, []))
    if excluded_cell_types and "cell_type" in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[~adata.obs["cell_type"].isin(excluded_cell_types)].copy()
        st.info(
            f"ℹ️ Excluding {len(excluded_cell_types)} cell-type label(s) set aside in "
            f"Step 8e (**{', '.join(excluded_cell_types)}**) -- {n_before - adata.n_obs:,} "
            f"of {n_before:,} cells are excluded from THIS step only; the full dataset "
            f"itself is unaffected."
        )


    groupable_columns = _get_groupable_obs_columns(adata)
    default_groupby = _recipe_default(project, "pseudobulk", "groupby_columns", ["sample"])
    default_groupby = [c for c in default_groupby if c in groupable_columns] or ["sample"]

    groupby_columns = st.multiselect(
        "Group cells into pseudobulk samples by:", options=groupable_columns,
        default=default_groupby, key="sc_downstream_pb_groupby",
        help=(
            "Choose ['sample'] alone for one pseudobulk sample per original "
            "sample (the simplest, most common case). Add a cluster/cell-type "
            "column (e.g. '" + cluster_key + "' or 'cell_type') for a "
            "cell-type-specific comparison -- one pseudobulk sample per "
            "sample x cell-type combination."
        ),
    )
    if not groupby_columns:
        st.info("Select at least one column above to continue.")
        return False

    default_min_cells = _recipe_default(
        project, "pseudobulk", "min_cells", dsm.DEFAULT_MIN_CELLS_PER_PSEUDOBULK_GROUP,
    )
    min_cells = st.number_input(
        "Minimum cells required per group (groups below this are excluded):",
        min_value=1, max_value=1000, value=int(default_min_cells),
        key="sc_downstream_pb_min_cells",
    )

    with st.expander("👀 Preview group sizes before aggregating", expanded=True):
        try:
            sizes_df = dsm.get_pseudobulk_group_sizes(adata, groupby_columns)
        except ValueError as e:
            st.error(f"⚠️ {e}")
            return False
        sizes_df = sizes_df.copy()
        sizes_df["would_be_excluded"] = sizes_df["n_cells"] < min_cells
        st.dataframe(sizes_df, use_container_width=True, hide_index=True)
        n_excluded_preview = int(sizes_df["would_be_excluded"].sum())
        if n_excluded_preview:
            st.warning(
                f"⚠️ {n_excluded_preview} group(s) above have fewer than {min_cells} "
                f"cells and would be EXCLUDED from the aggregated output."
            )

    default_layer = _recipe_default(project, "pseudobulk", "layer", "counts")
    layer_options = [l for l in ("counts",) if l in adata.layers]
    if not layer_options:
        st.error("⚠️ No raw 'counts' layer found -- this should have been set automatically in Step 1 (Combine Samples).")
        return False
    layer = layer_options[0]

    current_params = {
        "groupby_columns": sorted(groupby_columns), "min_cells": int(min_cells), "layer": layer,
    }
    is_current = scpm.check_downstream_step_current(project, "pseudobulk", current_params)
    last_result = st.session_state.get("sc_downstream_pseudobulk_result")

    if is_current and last_result is not None:
        st.success(
            f"✅ Already aggregated with these exact settings -- "
            f"{last_result['metadata_df'].shape[0]} pseudobulk sample(s)."
        )
    elif last_result is not None and not is_current:
        st.warning("⚠️ Your current settings differ from the last aggregation run above -- re-run below to apply them.")

    run_label = "🔄 Re-run Aggregation" if last_result is not None else "▶️ Aggregate to Pseudobulk"
    if st.button(run_label, key="sc_downstream_pb_aggregate_btn", type="primary"):
        with st.spinner("Aggregating pseudobulk counts..."):
            try:
                counts_df, metadata_df, excluded_groups = dsm.aggregate_pseudobulk(
                    adata, groupby_columns, layer=layer, min_cells=int(min_cells),
                )
            except ValueError as e:
                st.error(f"⚠️ {e}")
                return False
        st.session_state["sc_downstream_pseudobulk_result"] = {
            "counts_df": counts_df, "metadata_df": metadata_df, "excluded_groups": excluded_groups,
        }
        scpm.save_downstream_step_recipe(project, "pseudobulk", current_params)
        st.success(f"✅ Aggregated {metadata_df.shape[0]} pseudobulk sample(s).")
        if excluded_groups:
            st.warning(f"⚠️ {len(excluded_groups)} group(s) were excluded for having fewer than {min_cells} cells.")
        st.rerun()

    last_result = st.session_state.get("sc_downstream_pseudobulk_result")
    if last_result is not None:
        with st.expander("📋 Pseudobulk metadata (one row per pseudobulk sample)", expanded=True):
            st.dataframe(last_result["metadata_df"], use_container_width=True, hide_index=True)
            st.caption(
                "💡 This metadata table is exactly what you'll upload as the "
                "sample sheet in the Bulk RNA-Seq pipeline's own DESeq2 "
                "workflow -- double-check that each row's condition/group "
                "columns look correct here BEFORE moving on, since any error "
                "will carry straight through into that downstream analysis."
            )
        with st.expander("📋 Pseudobulk counts (genes x pseudobulk samples)", expanded=False):
            st.dataframe(last_result["counts_df"].head(50), use_container_width=True, hide_index=True)
            st.caption(f"Showing first 50 of {last_result['counts_df'].shape[0]:,} genes.")
        if last_result["excluded_groups"]:
            with st.expander(f"⚠️ {len(last_result['excluded_groups'])} excluded group(s)", expanded=False):
                st.dataframe(pd.DataFrame(last_result["excluded_groups"]), use_container_width=True, hide_index=True)

        st.markdown("**💾 Save this export for use in the Bulk RNA-Seq DESeq2 workflow**")
        export_name = st.text_input(
            "Export name:", key="sc_downstream_pb_export_name",
            placeholder="e.g. by_sample_and_celltype",
        )
        if st.button("💾 Save Pseudobulk Export", key="sc_downstream_pb_save_export_btn"):
            if not export_name.strip():
                st.error("Enter an export name first.")
            else:
                try:
                    dest_dir = scpm.downstream_pseudobulk_export_dir(project, export_name)
                    counts_path, metadata_path = dsm.save_pseudobulk_for_deseq2(
                        last_result["counts_df"], last_result["metadata_df"], dest_dir,
                    )
                    st.success(f"✅ Saved export '{export_name}' -- ready for the Bulk RNA-Seq DESeq2 workflow.")
                    st.caption(f"`{counts_path}`  \n`{metadata_path}`")
                except ValueError as e:
                    st.error(f"⚠️ {e}")

        scw._render_csv_download(
            last_result["counts_df"], "pseudobulk_counts", "sc_downstream_pb_counts",
            expander_label="⬇️ Download pseudobulk counts (.csv)",
        )
        scw._render_csv_download(
            last_result["metadata_df"], "pseudobulk_metadata", "sc_downstream_pb_metadata",
            expander_label="⬇️ Download pseudobulk metadata (.csv)",
        )

    existing_exports = scpm.list_downstream_pseudobulk_exports(project)
    if existing_exports:
        with st.expander(f"📂 Existing saved exports for this project ({len(existing_exports)})", expanded=False):
            for export in existing_exports:
                row_cols = st.columns([3, 1])
                row_cols[0].markdown(f"- `{export}`")
                if row_cols[1].button("🗑️ Delete", key=f"sc_downstream_pb_delete_export_{export}"):
                    scpm.delete_downstream_pseudobulk_export(project, export)
                    st.rerun()

    return is_current and last_result is not None


# ---------------------------------------------------------------------------
# Step 10: Compositional Analysis
# ---------------------------------------------------------------------------

def _render_compositional_step(project, adata, cluster_key):
    st.header("Step 10: Compositional Analysis")
    st.markdown(
        "Tests whether cell-type/cluster PROPORTIONS differ significantly "
        "between experimental groups/conditions -- e.g. 'does treatment shift "
        "the fraction of T cells vs. B cells?' Both methods below correctly "
        "treat each SAMPLE (not each cell) as the unit of replication."
    )
    with st.expander("ℹ️ What is this step, and why is it different from Step 9?", expanded=False):
        st.markdown(
            "**What it does:** Step 9 asks 'within a given cell type, are "
            "individual GENES expressed differently between conditions?' "
            "This step asks a completely different question: 'does the "
            "overall MIX of cell types itself shift between conditions?' -- "
            "e.g. does treatment cause a larger fraction of the tissue to "
            "become a particular immune cell type, independent of whether "
            "any individual gene within that cell type changed at all.\n\n"
            "**Why it matters:** These two questions can have opposite or "
            "unrelated answers -- a treatment could shift cell-type "
            "PROPORTIONS dramatically while barely changing gene expression "
            "within each cell type, or vice versa. Both are usually worth "
            "checking; skipping this step would miss a real category of "
            "biological effect that pure differential-expression testing "
            "cannot detect at all.\n\n"
            "**How to interpret the stacked bar plot:** each bar is one "
            "sample, and each colored segment is that sample's fraction "
            "belonging to one cell type -- visually compare the segment "
            "sizes for the SAME cell type across samples/conditions before "
            "even running a formal test; an obvious, consistent shift you "
            "can already see by eye is a good sign the statistical test "
            "below is picking up something real, not just noise.\n\n"
            "**Choosing a method:** 'propeller' (R/speckle) is the more "
            "statistically rigorous option, explicitly designed for this "
            "exact compositional-proportions problem, and is recommended "
            "when R is available. The simple proportion test is a lighter-"
            "weight fallback with fewer assumptions but generally less "
            "statistical power -- useful mainly when R/speckle isn't "
            "installed."
        )

    grouping_options = [cluster_key]
    if "cell_type" in adata.obs.columns:
        grouping_options.insert(0, "cell_type")
    cell_type_key = st.selectbox(
        "Cell type / cluster column to compute proportions over:",
        options=grouping_options, key="sc_downstream_comp_celltype_key",
        help="'cell_type' (Step 8's final labels) is offered first if available, since it's typically more interpretable than raw cluster numbers.",
    )

    # --- FIXED (2026-09-09): use the new, stricter
    # dsm.get_sample_level_condition_columns() instead of the broad
    # _get_groupable_obs_columns() -- see that function's own docstring
    # for the full rationale. This ONLY offers a column here if it is
    # genuinely CONSTANT within every sample (a real condition/donor/
    # batch attribute), which is what get_sample_level_group_mapping()
    # itself already requires further below -- so a per-cell QC metric
    # like "sum" (which varies cell-to-cell within a sample) can never
    # be offered here again, regardless of its name. cluster_key and
    # "cell_type" are explicitly excluded too, since those describe the
    # CELL being tested, not the sample's own experimental condition.
    candidate_condition_columns = dsm.get_sample_level_condition_columns(
        adata, sample_key="sample", excluded_columns={cluster_key, "cell_type"},
    )
    if not candidate_condition_columns:
        st.warning(
            "⚠️ No sample-level condition/group column found yet in this dataset's "
            "metadata -- compositional analysis requires a real per-sample "
            "experimental condition column (e.g. 'condition': Healthy vs. "
            "COVID-19). Complete **Step 1b (Label Samples with Metadata)** "
            "above to add one before this step can run."
        )
        return False

    group_column = st.selectbox(
        "Experimental group/condition column to compare across:",
        options=candidate_condition_columns, key="sc_downstream_comp_group_column",
    )

    sample_key = "sample"
    try:
        proportions_df, counts_df = dsm.compute_cell_type_proportions(adata, sample_key, cell_type_key)
    except ValueError as e:
        st.error(f"⚠️ {e}")
        return False

    with st.expander("📊 Composition overview (stacked bar plot)", expanded=True):
        barplot_df = dsm.get_composition_barplot_data(proportions_df)
        fig = px.bar(
            barplot_df, x="sample", y="proportion", color="cluster",
            labels={"cluster": cell_type_key},
        )
        fig.update_layout(height=450, barmode="stack", margin=dict(l=10, r=10, t=30, b=10))
        scw._render_plotly_chart(fig)
        scw._render_pdf_export(fig, "sc_downstream_composition_barplot", "composition_barplot")
        scw._render_csv_download(
            barplot_df, "composition_proportions", "sc_downstream_composition_barplot",
            expander_label="⬇️ Download composition data (.csv)",
        )

    try:
        sample_group_df = dsm.get_sample_level_group_mapping(adata, sample_key, group_column)
    except ValueError as e:
        st.error(f"⚠️ {e}")
        return False

    method_keys = list(dsm.COMPOSITIONAL_METHOD_OPTIONS.keys())
    default_method = _recipe_default(project, "compositional", "method", dsm.DEFAULT_COMPOSITIONAL_METHOD)
    propeller_ok = dsm.compositional_tools_available()
    if not propeller_ok and default_method == "propeller":
        default_method = "simple_test"

    method = st.radio(
        "Compositional testing method:", method_keys,
        format_func=lambda k: dsm.COMPOSITIONAL_METHOD_OPTIONS[k]["label"] + ("" if k != "propeller" or propeller_ok else " -- ⚠️ Rscript not found"),
        index=method_keys.index(default_method) if default_method in method_keys else 0,
        key="sc_downstream_comp_method",
    )
    st.caption(dsm.COMPOSITIONAL_METHOD_OPTIONS[method]["explanation"])

    if method == "propeller" and not propeller_ok:
        st.error("⚠️ Rscript was not found -- select 'Simple proportion test' instead, or install R + the `speckle` package.")
        return False

    require_all = (method == "simple_test")
    renamed_sample_group = sample_group_df.rename(columns={group_column: "group"})
    replication_check = dsm.check_compositional_replication(
        renamed_sample_group, "group", require_all_groups_replicated=require_all,
    )
    if replication_check["is_valid"]:
        st.success(replication_check["message"])
    else:
        st.warning(replication_check["message"])

    transform = dsm.DEFAULT_PROPELLER_TRANSFORM
    robust, trend = True, False
    if method == "propeller":
        col1, col2, col3 = st.columns(3)
        with col1:
            transform = st.radio("Transform:", ["logit", "asin"], key="sc_downstream_comp_transform", horizontal=True)
        with col2:
            robust = st.checkbox("Robust variance estimation", value=True, key="sc_downstream_comp_robust")
        with col3:
            trend = st.checkbox("Fit mean-variance trend", value=False, key="sc_downstream_comp_trend")

    current_params = {
        "method": method, "cell_type_key": cell_type_key, "group_column": group_column,
        "transform": transform if method == "propeller" else None,
        "robust": robust if method == "propeller" else None,
        "trend": trend if method == "propeller" else None,
    }
    is_current = scpm.check_downstream_step_current(project, "compositional", current_params)
    last_results_df = st.session_state.get("sc_downstream_compositional_results")

    if is_current and last_results_df is not None:
        st.success("✅ Already analyzed with these exact settings.")
    elif last_results_df is not None and not is_current:
        st.warning("⚠️ Your current settings differ from the last run above -- re-run below to apply them.")

    run_label = "🔄 Re-run Compositional Analysis" if last_results_df is not None else "▶️ Run Compositional Analysis"
    if st.button(run_label, key="sc_downstream_comp_run_btn", type="primary"):
        if not replication_check["is_valid"]:
            st.error("⚠️ Insufficient replication for this method -- see the message above.")
            return False
        if method == "propeller":
            output_dir = scpm.downstream_compositional_output_dir(project)
            work_dir = scpm.downstream_compositional_work_dir(project)
            with st.spinner("Running propeller (R/speckle)..."):
                success, log = dsm.run_propeller_analysis(
                    adata, sample_key, cell_type_key, group_column,
                    output_dir=output_dir, work_dir=work_dir,
                    transform=transform, robust=robust, trend=trend,
                )
            if not success:
                st.error(f"⚠️ {log}")
                return False
            results_df = dsm.read_propeller_results(output_dir)
            st.session_state["sc_downstream_compositional_log"] = log
        else:
            with st.spinner("Running simple proportion test..."):
                try:
                    results_df = dsm.run_simple_proportion_test(proportions_df, sample_group_df, group_column=group_column)
                except ValueError as e:
                    st.error(f"⚠️ {e}")
                    return False
        st.session_state["sc_downstream_compositional_results"] = results_df
        scpm.save_downstream_step_recipe(project, "compositional", current_params)
        st.success("✅ Compositional analysis complete.")
        st.rerun()

    last_results_df = st.session_state.get("sc_downstream_compositional_results")
    if last_results_df is not None:
        with st.expander("📋 Results", expanded=True):
            st.dataframe(last_results_df, use_container_width=True, hide_index=True)
            st.caption(
                "💡 **How to read this:** look for a small adjusted p-value "
                "(e.g. below 0.05) combined with a meaningful effect size/"
                "direction for a given cell type -- that combination is what "
                "indicates a real, statistically supported shift in that cell "
                "type's proportion between your groups, not just a single "
                "small p-value in isolation. As with any compositional "
                "analysis, a result driven by only one or two samples should "
                "be treated cautiously even if it clears a significance "
                "threshold -- cross-check against the stacked bar plot above "
                "to confirm the pattern looks consistent across samples "
                "within each group."
            )
            scw._render_csv_download(
                last_results_df, "compositional_results", "sc_downstream_comp_results",
                expander_label="⬇️ Download results (.csv)",
            )
        if method == "propeller" and st.session_state.get("sc_downstream_compositional_log"):
            with st.expander("📜 R log", expanded=False):
                st.code(st.session_state["sc_downstream_compositional_log"])

    return is_current and last_results_df is not None


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("📊 SC Downstream Analysis")
    st.markdown(
        "Phase 3: multi-sample normalization, dimensionality reduction, "
        "batch correction, clustering, and downstream analysis -- built on "
        "top of one or more samples that have already completed Step 6 "
        "(alignment) and Phase 2 (Cell-level QC)."
    )
    st.markdown("---")

    project = _render_project_selector()
    if not project:
        return
    st.markdown("---")

    existing_adata = _get_cached_adata(project)
    if existing_adata is not None:
        st.info(
            f"ℹ️ Resuming previous progress for this project -- "
            f"{existing_adata.n_obs:,} cells x {existing_adata.n_vars:,} genes "
            f"currently cached."
        )

    adata = _render_combine_step(project)
    if adata is None:
        return

    # --- NEW: Step 1b -- Label Samples with Metadata (2026-09-09) --
    # see _render_sample_metadata_step()'s own docstring for the full
    # rationale. Placed here (immediately after Combine succeeds, before
    # Normalization) since every later step -- including Step 9's
    # pseudobulk export and Step 10's compositional analysis -- can
    # benefit from a real per-sample condition column being available
    # as early as possible.
    st.markdown("---")
    _render_sample_metadata_step(project, adata)

    st.markdown("---")
    normalize_ready = _render_normalize_step(project, adata)
    if not normalize_ready:
        st.info("Complete Step 2 (Normalization) above to continue.")
        st.markdown("---")
        _render_checkpoint_controls(project, adata)
        _render_download_section(project, adata)
        return

    st.markdown("---")
    hvg_ready = _render_hvg_step(project, adata)
    if not hvg_ready:
        st.info("Complete Step 3 (HVG Selection) above to continue.")
        st.markdown("---")
        _render_checkpoint_controls(project, adata)
        _render_download_section(project, adata)
        return

    st.markdown("---")
    pca_ready = _render_pca_step(project, adata)
    if not pca_ready:
        st.info("Complete Step 4 (PCA) above to continue.")
        st.markdown("---")
        _render_checkpoint_controls(project, adata)
        _render_download_section(project, adata)
        return

    st.markdown("---")
    batch_ready, use_rep = _render_batch_correction_step(project, adata)
    if not batch_ready:
        st.info("Complete Step 5 (Batch Correction) above to continue.")
        st.markdown("---")
        _render_checkpoint_controls(project, adata)
        _render_download_section(project, adata)
        return

    st.markdown("---")
    cluster_key = _recipe_default(project, "clustering", "method", dsm.DEFAULT_CLUSTERING_METHOD)
    clustering_ready = _render_clustering_step(project, adata, use_rep)
    if not clustering_ready:
        st.info("Complete Step 6 (Clustering) above to continue.")
        st.markdown("---")
        _render_checkpoint_controls(project, adata)
        _render_download_section(project, adata)
        return

    st.markdown("---")
    _render_embedding_step(project, adata, use_rep, cluster_key)

    st.markdown("---")
    annotation_ready = _render_celltype_annotation_step(project, adata, cluster_key)

    st.markdown("---")
    pseudobulk_ready = _render_pseudobulk_step(project, adata, cluster_key)

    st.markdown("---")
    _render_compositional_step(project, adata, cluster_key)

    st.markdown("---")
    _render_checkpoint_controls(project, adata)
    _render_download_section(project, adata)

    if annotation_ready and pseudobulk_ready:
        st.markdown("---")
        st.success("🎉 Steps 1-10 complete for this project's current configuration.")
