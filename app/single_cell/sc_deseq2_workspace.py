"""
single_cell/sc_deseq2_workspace.py

Runs the Bulk RNA-Seq pipeline's own DESeq2 infrastructure (deseq2_manager.py)
on a single-cell pseudobulk export, and reuses differential_expression_workspace.py's
(dew) plotting/styling/export code DIRECTLY -- every dew._plot_*,
dew._render_plot_style_controls, dew._apply_plot_style, dew._render_csv_download,
dew._render_pdf_export, dew._render_group_label_renaming_controls,
dew._render_heatmap_color_controls, the gene-labeling/search system, and even
dew's own underscore-prefixed "_load_*" result-reader helpers (e.g.
dew._load_pca_coordinates) all take a plain directory path or DataFrame as
their argument -- NONE of them are coupled to project_manager (pm) -- so this
module can call them verbatim rather than reimplementing any plotting logic.

This is the companion to sc_ontology_workspace.py (split into two separate
pages/radio-buttons per the user's request, mirroring how dew/ow are two
separate pages in the bulk pipeline). Both share project/export/cell-type
selection via st.session_state keys (see _get_or_render_selection below), so
switching between the two pages does not force re-picking a project/export.

v1 SCOPE NOTE: the gene-ID-mapping panel here supports bitr conversion +
manual CSV upload only -- NOT the bulk pipeline's "layer 1" auto-derived
mapping from alignment_workspace.py's FASTA/GTF parsing, since that
infrastructure hasn't been verified to exist on the single-cell side. This
covers the core need (readable gene names in plots) without guessing at
unseen sc-side reference infrastructure.
"""
import os

import pandas as pd
import streamlit as st

import sc_project_manager as scpm
import deseq2_manager as dm
import gene_id_mapper as gim
import reference_manager as rm
import differential_expression_workspace as dew

WORKSPACE_KEY = "sc_deseq2"


# ---------------------------------------------------------------------------
# Shared project/export/cell-type selection (also used by
# sc_ontology_workspace.py via the same session_state keys, so picking a
# project/export on one page carries over to the other).
# ---------------------------------------------------------------------------
_SESSION_PROJECT_KEY = "_sc_de_bridge_project"
_SESSION_EXPORT_KEY = "_sc_de_bridge_export"
_SESSION_CELL_TYPE_SUFFIX_KEY = "_sc_de_bridge_cell_type_suffix"

CELL_TYPE_COLUMN_CANDIDATES = ["cell_type", "celltype", "cluster", "cluster_id", "cell_type_label"]


def _detect_cell_type_column(meta_df):
    lower_to_real = {c.lower(): c for c in meta_df.columns}
    for candidate in CELL_TYPE_COLUMN_CANDIDATES:
        if candidate in lower_to_real:
            return lower_to_real[candidate]
    return None


def render_project_export_and_cell_type_selector():
    """
    Steps 0 and 0.5: pick a single-cell project + saved pseudobulk export,
    then (if the export includes a cell-type-like column) require picking
    ONE cell type before continuing -- see this module's docstring and the
    original prototype's extensive rationale for why this is REQUIRED, not
    optional (pseudoreplication risk from pooling multiple cell types per
    donor into one model; Crowell et al. 2020's muscat package;
    Squair et al. 2021).

    Returns (project, export_name, counts_df, meta_df, cell_type_suffix)
    -- counts_df/meta_df are already subset to the chosen cell type if
    applicable. Any of the first two may be None if selection isn't
    complete yet, in which case the caller should stop rendering further.
    """
    st.header("Step 0: Choose Your Project and Pseudobulk Export")
    st.markdown(
        "This page runs DESeq2 on a **pseudobulk export** already "
        "created in Step 9 of the single-cell downstream analysis "
        "workflow -- aggregating counts across cells (e.g. per sample, "
        "or per sample + cell type) into a bulk-like counts matrix, "
        "then reusing the exact same DESeq2 infrastructure the Bulk "
        "RNA-Seq pipeline uses."
    )
    existing_projects = scpm.list_projects()
    if not existing_projects:
        st.info("No single-cell projects found yet. Start one from the 🧫 Single-cell RNA-Seq page first.")
        return None, None, None, None, ""

    project = st.selectbox(
        "Single-cell project:", options=existing_projects, key=f"{WORKSPACE_KEY}_project_select",
    )
    export_names = scpm.list_downstream_pseudobulk_exports(project)
    if not export_names:
        st.warning(
            f"⚠️ Project `{project}` doesn't have any saved pseudobulk exports yet. "
            "Go to **📊 SC Analysis: Clustering & Cell Annotation** and complete Step 9 first."
        )
        if st.button("⬅️ Go to SC Analysis: Clustering & Cell Annotation", key=f"{WORKSPACE_KEY}_gate_back_btn"):
            import singlecell_workspace as scw
            st.session_state["nav_request"] = scw.SC_DOWNSTREAM_ANALYSIS_OPTION
            st.rerun()
        return project, None, None, None, ""

    export_name = st.selectbox(
        "Pseudobulk export:", options=export_names, key=f"{WORKSPACE_KEY}_export_select",
        help="Each export represents one saved 'grouped by' choice from Step 9 (e.g. grouped by sample, or by sample + cell type).",
    )
    export_dir = scpm.downstream_pseudobulk_export_dir(project, export_name)
    counts_path = os.path.join(export_dir, "pseudobulk_counts.csv")
    metadata_path = os.path.join(export_dir, "pseudobulk_metadata.csv")
    if not (os.path.isfile(counts_path) and os.path.isfile(metadata_path)):
        st.error(f"⚠️ Could not find pseudobulk_counts.csv / pseudobulk_metadata.csv for export `{export_name}`.")
        return project, None, None, None, ""

    counts_df = pd.read_csv(counts_path)
    meta_df = pd.read_csv(metadata_path)
    st.caption(f"Loaded {len(counts_df):,} genes across {len(meta_df)} pseudobulk sample(s).")

    cell_type_col = _detect_cell_type_column(meta_df)
    cell_type_suffix = ""
    if cell_type_col is not None:
        st.markdown("---")
        st.header("Step 0.5: Choose a Cell Type to Analyze")
        st.warning(
            f"⚠️ **This pseudobulk export includes a cell-type column "
            f"(`{cell_type_col}`), so each of your original samples "
            f"contributed MULTIPLE rows to this table -- one per cell "
            f"type.**\n\n"
            "**Why you must pick just one cell type before continuing:** "
            "if all cell types were combined into a single DESeq2 model "
            "together with your condition of interest, the multiple rows "
            "coming from the same original sample/donor would not be "
            "independent observations -- they'd share the same donor "
            "biology, which risks *pseudoreplication* (artificially "
            "inflated significance). This isn't a batch effect (a "
            "technical/nuisance factor); cell type is a real biological "
            "factor, but the correct, published way to handle it "
            "(Crowell et al. 2020, *Nature Communications* -- the "
            "`muscat` R package; Squair et al. 2021, *Nature "
            "Communications*, \"Confronting false discoveries in "
            "single-cell differential expression\") is to test your "
            "condition of interest **separately, within each cell "
            "type**, using your original samples/donors as the unit of "
            "replication.\n\nPick a cell type below. Each cell type's "
            "results are saved separately and automatically."
        )
        available_cell_types = sorted(meta_df[cell_type_col].astype(str).dropna().unique().tolist())
        if not available_cell_types:
            st.error(f"⚠️ Column `{cell_type_col}` has no usable values.")
            st.stop()
        chosen_cell_type = st.selectbox(
            "Cell type to analyze:", options=available_cell_types, key=f"{WORKSPACE_KEY}_cell_type_select",
        )
        keep_mask = meta_df[cell_type_col].astype(str) == chosen_cell_type
        meta_df = meta_df.loc[keep_mask].reset_index(drop=True)
        if len(meta_df) < 2:
            st.error(
                f"⚠️ Only {len(meta_df)} pseudobulk sample(s) remain after filtering to "
                f"`{chosen_cell_type}` -- DESeq2 needs at least 2 samples per group."
            )
            st.stop()
        sample_cols_to_keep = ["gene_id"] + meta_df["sample"].tolist()
        counts_df = counts_df[[c for c in sample_cols_to_keep if c in counts_df.columns]]
        meta_df = meta_df.drop(columns=[cell_type_col])
        cell_type_suffix = f"__{chosen_cell_type}"
        st.success(f"✅ Subset to **{chosen_cell_type}**: {len(meta_df)} pseudobulk sample(s) remaining.")

    # Persist to session_state so sc_ontology_workspace.py can pick up the
    # exact same selection without the user re-choosing anything there.
    st.session_state[_SESSION_PROJECT_KEY] = project
    st.session_state[_SESSION_EXPORT_KEY] = export_name
    st.session_state[_SESSION_CELL_TYPE_SUFFIX_KEY] = cell_type_suffix

    return project, export_name, counts_df, meta_df, cell_type_suffix


def get_shared_selection():
    """
    Used by sc_ontology_workspace.py to read back whatever project/export/
    cell-type was last selected on THIS page, without re-rendering the
    selector UI. Returns (project, export_name, cell_type_suffix), any of
    which may be None/"" if nothing has been selected here yet.
    """
    return (
        st.session_state.get(_SESSION_PROJECT_KEY),
        st.session_state.get(_SESSION_EXPORT_KEY),
        st.session_state.get(_SESSION_CELL_TYPE_SUFFIX_KEY, ""),
    )


# ---------------------------------------------------------------------------
# Steps 1-4: Analysis Settings (filtering, design, contrasts, run) --
# adapted from dew._render_analysis_settings, export-scoped rather than
# project-scoped. Kept as a separate function (not inlined in render())
# for the exact same reason dew's own version is separate -- see that
# function's docstring: early `return`s here must only exit THIS function,
# never skip rendering Step 5 below when the expander happens to be
# auto-collapsed (i.e. when this export already has results).
# ---------------------------------------------------------------------------
def _render_sc_analysis_settings(project, export_name, cell_type_suffix, counts_df, meta_df):
    deseq2_out_dir = scpm.sc_deseq2_output_dir(project, export_name + cell_type_suffix)
    deseq2_work_dir = scpm.sc_deseq2_work_dir(project, export_name + cell_type_suffix)
    os.makedirs(deseq2_work_dir, exist_ok=True)
    counts_matrix_path = os.path.join(deseq2_work_dir, "filtered_pseudobulk_counts.csv")
    counts_df.to_csv(counts_matrix_path, index=False)

    existing_results = dm.list_contrast_results(deseq2_out_dir)
    saved_config = scpm.get_sc_deseq2_config(project, export_name + cell_type_suffix) or {}

    with st.expander(
        "⚙️ Steps 1-4: Analysis Settings (filtering, design, contrasts, run)",
        expanded=not (existing_results and saved_config),
    ):
        st.header("Step 1: Remove Low-Count Genes")
        col1, col2 = st.columns(2)
        with col1:
            min_count = st.number_input(
                "Minimum read count:", min_value=0, value=10, step=1, key=f"{WORKSPACE_KEY}_min_count_input",
            )
        with col2:
            min_samples = st.number_input(
                "Minimum number of samples:", min_value=1, value=min(2, len(meta_df)), step=1,
                max_value=len(meta_df), key=f"{WORKSPACE_KEY}_min_samples_input",
            )
        filter_preview = dm.preview_low_count_filter(counts_df, min_count, min_samples)
        st.caption(
            f"With these settings: **{filter_preview['genes_kept']:,} of "
            f"{filter_preview['total_genes']:,} genes** would be kept ({filter_preview['pct_kept']}%)."
        )
        st.markdown("---")

        st.header("Step 2: Set Up Your Experimental Design")
        available_columns = [c for c in meta_df.columns if c not in ("sample", "n_cells_aggregated")]
        if not available_columns:
            st.error(
                "⚠️ No usable condition column found in this pseudobulk export's metadata. "
                "This can happen if Step 9 was grouped by 'sample' only and no per-sample "
                "condition metadata (e.g. via merge_sample_level_metadata()) was ever "
                "attached to this project. Re-export from Step 9 after confirming your "
                "project has real condition metadata attached."
            )
            return
        design_columns = st.multiselect(
            "Which column(s) represent your experimental condition(s) of interest?",
            options=available_columns, default=[], key=f"{WORKSPACE_KEY}_design_columns_select",
        )
        batch_column_options = ["(none)"] + [c for c in available_columns if c not in design_columns]
        batch_column_choice = st.selectbox(
            "Does this project have a batch effect column? (optional)",
            options=batch_column_options, key=f"{WORKSPACE_KEY}_batch_column_select",
        )
        batch_column = None if batch_column_choice == "(none)" else batch_column_choice
        if not design_columns:
            st.warning("⚠️ Select at least one design column to continue.")
            return
        columns_to_validate = list(design_columns) + ([batch_column] if batch_column else [])
        validation = dm.validate_design_columns(meta_df, columns_to_validate)
        constant_columns = [col for col, info in validation.items() if info["is_constant"]]
        if constant_columns:
            st.error(
                f"⚠️ **{', '.join(constant_columns)}** has the exact same value for every "
                "pseudobulk sample, so it can't be used as a design/batch factor."
            )
            return

        columns_to_check_for_missing = list(design_columns) + ([batch_column] if batch_column else [])
        meta_df, all_resolved = dew._render_missing_value_resolution(meta_df, columns_to_check_for_missing)
        if not all_resolved:
            return

        interaction_terms = []
        if len(design_columns) >= 2:
            interaction_options = dm.build_interaction_term_options(design_columns)
            interaction_terms = st.multiselect(
                "Include interaction term(s) in the model? (optional)",
                options=interaction_options, default=[], key=f"{WORKSPACE_KEY}_interaction_terms_select",
            )
        full_formula_terms = dm.build_full_formula_terms(design_columns, batch_column, interaction_terms)
        st.caption(f"**Full model formula:** ~ {' + '.join(full_formula_terms)}")

        replication = dm.check_replication(meta_df, design_columns, batch_column, interaction_terms)
        if not replication["is_valid"]:
            st.error(
                f"⚠️ The following group(s) only have 1 pseudobulk sample: "
                f"**{', '.join(replication['under_replicated_groups'])}**. DESeq2 needs at "
                "least 2 samples per group -- with pseudobulk data, this usually means "
                "grouping by fewer columns in Step 9, or combining categories."
            )
            return
        st.markdown("---")

        st.header("Step 3: Choose Your Analysis Type")
        test_type_choice = st.radio(
            "Which type of test do you want to run?",
            ["Wald test (pairwise comparisons)", "LRT (omnibus / ANOVA-style test)"],
            key=f"{WORKSPACE_KEY}_test_type_radio",
        )
        test_type = "wald" if test_type_choice.startswith("Wald") else "lrt"
        contrasts, lrt_tests = [], []
        if test_type == "wald":
            primary_column = st.selectbox(
                "Which design column should contrasts be based on?",
                options=design_columns, key=f"{WORKSPACE_KEY}_primary_contrast_column",
            )
            levels = sorted(meta_df[primary_column].astype(str).dropna().unique().tolist())
            if len(levels) < 2:
                st.error(f"⚠️ Column '{primary_column}' needs at least 2 distinct values.")
                return
            reference_level = st.selectbox(
                "Which level is your reference/control group?", options=levels,
                key=f"{WORKSPACE_KEY}_reference_level_select",
            )
            comparison_levels = st.multiselect(
                "Which level(s) should be compared against the reference?",
                options=[lv for lv in levels if lv != reference_level],
                default=[lv for lv in levels if lv != reference_level],
                key=f"{WORKSPACE_KEY}_comparison_levels_select",
            )
            contrasts = [
                {"name": f"{lv}_vs_{reference_level}", "column": primary_column, "level1": lv, "level2": reference_level}
                for lv in comparison_levels
            ]
            if not contrasts:
                st.warning("⚠️ Select at least one comparison level to continue.")
                return
        else:
            lrt_term_options = list(design_columns) + interaction_terms
            terms_to_test = st.multiselect(
                "Which term(s) should be tested?", options=lrt_term_options,
                default=lrt_term_options[:1] if lrt_term_options else [], key=f"{WORKSPACE_KEY}_lrt_terms_select",
            )
            lrt_tests = [
                {"name": f"{term.replace(':', '_x_')}_effect",
                 "reduced_terms": dm.build_reduced_formula_terms(full_formula_terms, [term])}
                for term in terms_to_test
            ]
            if not lrt_tests:
                st.warning("⚠️ Select at least one term to test to continue.")
                return
        st.markdown("---")

        st.header("Step 4: Run DESeq2")
        if not dm.deseq2_tools_available():
            st.error("⚠️ Rscript was not found on this system.")
            return
        resolved_metadata_path = os.path.join(deseq2_work_dir, "resolved_metadata.csv")
        meta_df.to_csv(resolved_metadata_path, index=False)
        results_already_exist_now = len(dm.list_contrast_results(deseq2_out_dir)) > 0
        run_label = "🔄 Re-run DESeq2" if results_already_exist_now else "🚀 Run DESeq2 Analysis"
        if st.button(run_label, key=f"{WORKSPACE_KEY}_run_deseq2_btn"):
            with st.spinner("Running DESeq2... this may take a few minutes."):
                success, log = dm.run_deseq2_analysis(
                    counts_matrix_path, resolved_metadata_path, design_columns, batch_column,
                    min_count, min_samples, deseq2_out_dir, deseq2_work_dir,
                    test_type=test_type, interaction_terms=interaction_terms,
                    contrasts=contrasts, lrt_tests=lrt_tests,
                )
            if not success:
                st.error("DESeq2 failed. Details below:")
                st.code(log)
                return
            st.success("✅ DESeq2 completed successfully.")
            st.code(log)
            scpm.save_sc_deseq2_config(project, export_name + cell_type_suffix, {
                "design_columns": design_columns, "batch_column": batch_column,
                "interaction_terms": interaction_terms, "test_type": test_type,
                "contrasts": contrasts, "lrt_tests": lrt_tests,
                "min_count": min_count, "min_samples": min_samples,
            })


# ---------------------------------------------------------------------------
# v1-scope gene-ID-mapping panel (bitr conversion + manual CSV upload only
# -- see module docstring for why the bulk "layer 1 auto-derived" mapping
# isn't ported here).
# ---------------------------------------------------------------------------
def _render_sc_gene_id_mapping_panel(project, export_name, cell_type_suffix, counts_df):
    full_export_name = export_name + cell_type_suffix

    # Layer 1: GTF-auto-derived map from this project's reference (see
    # singlecell_workspace.py's Step 5 reference-confirmation step) --
    # loaded first, exactly mirroring the bulk pipeline's own 3-layer
    # design (auto -> bitr -> manual upload).
    auto_gene_name_map = {}
    reference_map_path = scpm.sc_reference_gene_symbol_map_path(project)
    if os.path.exists(reference_map_path):
        try:
            ref_map_df = pd.read_csv(reference_map_path)
            auto_gene_name_map = dict(zip(ref_map_df["gene_id"].astype(str), ref_map_df["gene_name"].astype(str)))
        except Exception:
            auto_gene_name_map = {}

    # Layer 2 override: any prior bitr conversion saved specifically for
    # THIS pseudobulk export, layered on top of layer 1.
    auto_map_path = scpm.sc_gene_symbol_map_path(project, full_export_name)
    if os.path.exists(auto_map_path):
        try:
            auto_map_df = pd.read_csv(auto_map_path)
            auto_gene_name_map.update(dict(zip(auto_map_df["gene_id"].astype(str), auto_map_df["gene_name"].astype(str))))
        except Exception:
            pass

    # Layer 2 override: any prior bitr conversion saved specifically for
    # THIS pseudobulk export, layered on top of layer 1.
    auto_map_path = scpm.sc_gene_symbol_map_path(project, full_export_name)
    if os.path.exists(auto_map_path):
        try:
            auto_map_df = pd.read_csv(auto_map_path)
            auto_gene_name_map = dict(zip(auto_map_df["gene_id"].astype(str), auto_map_df["gene_name"].astype(str)))
        except Exception:
            auto_gene_name_map = {}

    gene_ids_all = sorted(set(counts_df["gene_id"].astype(str)))
    n_total = len(gene_ids_all)
    unresolved_ids = [gid for gid in gene_ids_all if auto_gene_name_map.get(gid, gid) == gid]
    n_resolved = n_total - len(unresolved_ids)

    with st.expander("🏷️ Gene ID → Gene Name Mapping", expanded=False):
        if auto_gene_name_map and n_resolved > 0:
            pct = (n_resolved / n_total * 100) if n_total else 0.0
            st.success(f"✅ **{n_resolved:,} of {n_total:,} gene IDs ({pct:.1f}%)** have a readable name.")
        else:
            st.info("ℹ️ No gene names resolved yet -- gene IDs are shown as-is for now.")

        st.markdown("##### 🧬 Convert Gene IDs (Bioconductor `bitr`)")
        detection_pool = unresolved_ids if unresolved_ids else gene_ids_all
        detection = gim.detect_id_type(detection_pool)
        detected_type = detection["detected_type"]

        species_choices = gim.orgdb_species_choices(rm)
        species_keys = list(species_choices.keys())
        previous_species = scpm.get_sc_ontology_species(project, full_export_name)
        default_species_index = species_keys.index(previous_species) if previous_species in species_keys else 0

        st.caption(
            f"Auto-detected ID type: **{detected_type}** ({detection['match_fraction'] * 100:.0f}% matched -- "
            f"e.g. `{'`, `'.join(detection['example_ids'][:3])}`)."
        )
        col_species, col_from, col_to = st.columns(3)
        with col_species:
            picked_species = st.selectbox(
                "Species (annotation database):", options=species_keys, index=default_species_index,
                format_func=lambda k: species_choices[k], key=f"{WORKSPACE_KEY}_gene_map_species_select",
            )
        with col_from:
            from_type_options = gim.COMMON_KEY_TYPES
            from_default_idx = from_type_options.index(detected_type) if detected_type in from_type_options else 0
            picked_from_type = st.selectbox(
                "My gene IDs are currently in this format:", options=from_type_options, index=from_default_idx,
                key=f"{WORKSPACE_KEY}_gene_map_from_type_select",
            )
        with col_to:
            to_type_options = gim.COMMON_KEY_TYPES
            default_to_type = gim.symbol_keytype_for_species(picked_species)
            to_default_idx = to_type_options.index(default_to_type) if default_to_type in to_type_options else 0
            picked_to_type = st.selectbox(
                "Convert them to:", options=to_type_options, index=to_default_idx,
                key=f"{WORKSPACE_KEY}_gene_map_to_type_select",
            )

        orgdb_package = gim.ORGDB_PACKAGES.get(picked_species)
        work_dir = scpm.sc_gene_id_mapping_work_dir(project, full_export_name)

        col_fill, col_override = st.columns(2)
        with col_fill:
            fill_disabled = not unresolved_ids
            if st.button(
                f"✨ Fill in {len(unresolved_ids):,} missing name(s) only",
                key=f"{WORKSPACE_KEY}_bitr_fill_missing_btn", disabled=fill_disabled,
            ):
                with st.spinner(f"Converting {len(unresolved_ids):,} gene ID(s) via bitr()..."):
                    result = gim.run_bitr_conversion(unresolved_ids, picked_from_type, picked_to_type, orgdb_package, work_dir)
                if result["success"]:
                    auto_gene_name_map.update(result["mapping"])
                    rm.save_gene_symbol_map_csv(auto_gene_name_map, auto_map_path)
                    scpm.save_sc_gene_id_mapping_meta(project, full_export_name, {
                        "source": "bitr_fill", "from_type": picked_from_type, "to_type": picked_to_type,
                        "orgdb_package": orgdb_package, "n_converted": result["n_converted"], "n_total": result["n_total"],
                    })
                    st.success(f"✅ {result['message']}")
                else:
                    st.error(f"⚠️ {result['message']}")
        with col_override:
            if st.button("🔁 Convert ALL gene IDs (override existing mapping)", key=f"{WORKSPACE_KEY}_bitr_override_all_btn"):
                with st.spinner(f"Converting {n_total:,} gene ID(s) via bitr()..."):
                    result = gim.run_bitr_conversion(gene_ids_all, picked_from_type, picked_to_type, orgdb_package, work_dir)
                if result["success"]:
                    auto_gene_name_map = result["mapping"]
                    rm.save_gene_symbol_map_csv(auto_gene_name_map, auto_map_path)
                    scpm.save_sc_gene_id_mapping_meta(project, full_export_name, {
                        "source": "bitr_full_override", "from_type": picked_from_type, "to_type": picked_to_type,
                        "orgdb_package": orgdb_package, "n_converted": result["n_converted"], "n_total": result["n_total"],
                    })
                    st.success(f"✅ {result['message']}")
                else:
                    st.error(f"⚠️ {result['message']}")

        if not gim.bitr_tools_available():
            st.caption("ℹ️ Rscript wasn't found -- the buttons above need R with `clusterProfiler` and the relevant Bioconductor annotation package installed.")

        st.markdown("---")
        st.caption(
            "Alternatively, upload your own 2-column CSV (columns named exactly 'gene_id' "
            "and 'gene_name') -- this stays local to your current session and takes "
            "priority over everything above."
        )
        gene_map_file = st.file_uploader(
            "Gene ID → Gene Name CSV (session-only override):", type=["csv"],
            key=f"{WORKSPACE_KEY}_gene_name_mapping_upload",
        )
        if gene_map_file is not None:
            try:
                map_df = pd.read_csv(gene_map_file)
                if "gene_id" in map_df.columns and "gene_name" in map_df.columns:
                    st.session_state[f"{WORKSPACE_KEY}_gene_name_map"] = dict(
                        zip(map_df["gene_id"].astype(str), map_df["gene_name"].astype(str))
                    )
                    st.success(f"✅ Loaded {len(st.session_state[f'{WORKSPACE_KEY}_gene_name_map']):,} mapping(s) for this session.")
                else:
                    st.error("⚠️ This CSV must have columns named exactly 'gene_id' and 'gene_name'.")
            except Exception as e:
                st.error(f"⚠️ Could not read this file: {e}")

    return st.session_state.get(f"{WORKSPACE_KEY}_gene_name_map") or auto_gene_name_map


# ---------------------------------------------------------------------------
# Step 5: Explore Results -- ported wholesale from dew.render()'s own Step
# 5, swapping every pm.deseq2_* call for the already-built scpm.sc_deseq2_*
# equivalent. All plotting/styling/export calls below are dew.* function
# calls, reused verbatim -- see this module's docstring.
# ---------------------------------------------------------------------------
def _render_step5_results(project, export_name, cell_type_suffix, meta_df, saved_config, counts_df):
    full_export_name = export_name + cell_type_suffix
    deseq2_out_dir = scpm.sc_deseq2_output_dir(project, full_export_name)
    batch_column = saved_config.get("batch_column")

    st.header("Step 5: Explore Your Results")

    available_contrasts = dm.list_contrast_results(deseq2_out_dir)
    if not available_contrasts:
        st.info("No results found yet -- run DESeq2 above.")
        return

    # --- PCA ---
    st.subheader("🔬 PCA: Sample Similarity")
    with st.expander("ℹ️ How to read this plot"):
        st.markdown(
            "Each point is one pseudobulk sample. Samples more similar in overall gene "
            "expression appear closer together. Clear clustering by condition is a good "
            "sign; clustering by batch instead suggests a batch effect (already accounted "
            "for statistically in your results via the design formula)."
        )
    pca_df, pct_var = dew._load_pca_coordinates(deseq2_out_dir)
    if pca_df is not None:
        color_options = [c for c in pca_df.columns if c not in ("PC1", "PC2", "PC3", "PC4", "sample")]
        if color_options:
            color_by = st.selectbox("Color points by:", options=color_options, key=f"{WORKSPACE_KEY}_pca_color_select")
            group_values = sorted(pca_df[color_by].astype(str).unique())
            pca_group_label_map = dew._render_group_label_renaming_controls(
                f"{WORKSPACE_KEY}_pca_{color_by}", group_values,
                label=f"Rename `{color_by}` group labels shown in the legend (optional)",
            )

            pca_adj_df, pct_var_adj = (None, None)
            if batch_column:
                pca_adj_df, pct_var_adj = dew._load_pca_coordinates(deseq2_out_dir, batch_adjusted=True)

            if pca_adj_df is not None:
                st.info(
                    "🔀 Two views: **Before** and **After** batch correction "
                    "(visualization only -- limma's removeBatchEffect). The actual "
                    "DESeq2 test already accounts for batch via the design formula "
                    "regardless of what's shown here."
                )
                st.markdown("#### 🔹 Before Batch Correction")
                show_labels_before = st.checkbox(
                    "🏷️ Show sample name labels on points", value=True, key=f"{WORKSPACE_KEY}_pca_before_show_labels",
                    help="Turn off to declutter a crowded plot -- sample names still show on hover regardless.",
                )
                show_ellipse_before, confidence_before = dew._render_confidence_ellipse_controls(f"{WORKSPACE_KEY}_pca_before")
                fig_before, ellipse_stats_before = dew._plot_pca(
                    pca_df, pct_var, color_by, show_confidence_ellipse=show_ellipse_before, confidence_level=confidence_before,
                    show_labels=show_labels_before,
                )               
                style_before = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_pca_before", group_values=group_values)
                dew._apply_plot_style(fig_before, style_before, default_title="PCA: Sample Similarity — Before Batch Correction")
                dew._apply_group_label_renaming(fig_before, pca_group_label_map)
                for es in ellipse_stats_before:
                    es["group"] = pca_group_label_map.get(es["group"], es["group"])
                dew._render_plotly_chart(fig_before)
                dew._render_ellipse_stats(ellipse_stats_before, f"{WORKSPACE_KEY}_pca_before")
                dew._render_pdf_export(fig_before, f"{WORKSPACE_KEY}_pca_before", "sc_pca_before_batch_correction")
                dew._render_csv_download(pca_df, "sc_pca_coordinates_before", f"{WORKSPACE_KEY}_pca_before", expander_label="⬇️ Download PCA Coordinates (Before).csv")

                st.markdown("#### 🔸 After Batch Correction")
                show_labels_after = st.checkbox(
                    "🏷️ Show sample name labels on points", value=True,
                    key=f"{WORKSPACE_KEY}_pca_after_show_labels",
                    help="Turn off to declutter a crowded plot -- sample names still show on hover regardless.",
                )
                show_ellipse_after, confidence_after = dew._render_confidence_ellipse_controls(f"{WORKSPACE_KEY}_pca_after")
                fig_after, ellipse_stats_after = dew._plot_pca(
                    pca_adj_df, pct_var_adj, color_by, show_confidence_ellipse=show_ellipse_after, confidence_level=confidence_after,
                    show_labels=show_labels_after,
                )                
                style_after = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_pca_after", group_values=group_values)
                dew._apply_plot_style(fig_after, style_after, default_title="PCA: Sample Similarity — After Batch Correction (Visualization Only)")
                dew._apply_group_label_renaming(fig_after, pca_group_label_map)
                for es in ellipse_stats_after:
                    es["group"] = pca_group_label_map.get(es["group"], es["group"])
                dew._render_plotly_chart(fig_after)
                dew._render_ellipse_stats(ellipse_stats_after, f"{WORKSPACE_KEY}_pca_after")
                dew._render_pdf_export(fig_after, f"{WORKSPACE_KEY}_pca_after", "sc_pca_after_batch_correction")
                dew._render_csv_download(pca_adj_df, "sc_pca_coordinates_after", f"{WORKSPACE_KEY}_pca_after", expander_label="⬇️ Download PCA Coordinates (After).csv")

                variance_rows = []
                n_pcs = min(len(pct_var or []), len(pct_var_adj or []), 4)
                for i in range(n_pcs):
                    variance_rows.append({
                        "Principal Component": f"PC{i + 1}", "% Variance (Before)": pct_var[i],
                        "% Variance (After)": pct_var_adj[i], "Change": round(pct_var_adj[i] - pct_var[i], 1),
                    })
                if variance_rows:
                    st.dataframe(pd.DataFrame(variance_rows), use_container_width=True, hide_index=True)
                batch_eta_before = dm.compute_batch_variance_metric(pca_df, batch_column)
                batch_eta_after = dm.compute_batch_variance_metric(pca_adj_df, batch_column)
                eta_rows = []
                if batch_eta_before and batch_eta_after:
                    eta_rows = [
                        {"Principal Component": pc, "% of Variance Explained by Batch (Before)": batch_eta_before.get(pc, 0.0),
                         "% of Variance Explained by Batch (After)": batch_eta_after.get(pc, 0.0)}
                        for pc in ("PC1", "PC2") if pc in batch_eta_before
                    ]
                    st.dataframe(pd.DataFrame(eta_rows), use_container_width=True, hide_index=True)
                pca_qc_export_df = dm.build_pca_qc_export(variance_rows, eta_rows)
                if not pca_qc_export_df.empty:
                    dew._render_csv_download(pca_qc_export_df, "sc_pca_qc_summary", f"{WORKSPACE_KEY}_pca_qc_summary", expander_label="⬇️ Download PCA QC Summary.csv")
            else:
                show_labels_single = st.checkbox(
                    "🏷️ Show sample name labels on points", value=True,
                    key=f"{WORKSPACE_KEY}_pca_single_show_labels",
                    help="Turn off to declutter a crowded plot -- sample names still show on hover regardless.",
                )
                show_ellipse_single, confidence_single = dew._render_confidence_ellipse_controls(f"{WORKSPACE_KEY}_pca_single")
                fig_single, ellipse_stats_single = dew._plot_pca(
                    pca_df, pct_var, color_by, show_confidence_ellipse=show_ellipse_single, confidence_level=confidence_single,
                    show_labels=show_labels_single,
                )                
                style_single = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_pca_single", group_values=group_values)
                dew._apply_plot_style(fig_single, style_single, default_title="PCA: Sample Similarity")
                dew._apply_group_label_renaming(fig_single, pca_group_label_map)
                for es in ellipse_stats_single:
                    es["group"] = pca_group_label_map.get(es["group"], es["group"])
                dew._render_plotly_chart(fig_single)
                dew._render_ellipse_stats(ellipse_stats_single, f"{WORKSPACE_KEY}_pca_single")
                dew._render_pdf_export(fig_single, f"{WORKSPACE_KEY}_pca_single", "sc_pca_sample_similarity")
                dew._render_csv_download(pca_df, "sc_pca_coordinates", f"{WORKSPACE_KEY}_pca_single", expander_label="⬇️ Download PCA Coordinates.csv")
    st.markdown("---")

    # --- Sample distance heatmap ---
    st.subheader("🌡️ Sample Distance Heatmap")
    distance_df = dew._load_sample_distance_matrix(deseq2_out_dir)
    if distance_df is not None:
        heatmap_group_columns = [c for c in meta_df.columns if c != "sample"]
        default_group_index = 0
        if saved_config.get("design_columns") and saved_config["design_columns"][0] in heatmap_group_columns:
            default_group_index = heatmap_group_columns.index(saved_config["design_columns"][0])
        heatmap_group_column, heatmap_show_group_annotation = None, True
        if heatmap_group_columns:
            col_group_select, col_hide_annotation = st.columns([2, 1.2])
            with col_group_select:
                heatmap_group_column = st.selectbox(
                    "Color/label samples by which metadata column?", options=heatmap_group_columns,
                    index=default_group_index, key=f"{WORKSPACE_KEY}_sample_dist_check_column",
                )
            with col_hide_annotation:
                heatmap_show_group_annotation = not st.checkbox(
                    "Hide grouping bar and legend", value=False, key=f"{WORKSPACE_KEY}_sample_dist_hide_group_annotation",
                )
        dist_colorscale, dist_reverse = dew._render_heatmap_color_controls(
            f"{WORKSPACE_KEY}_sample_dist", dew.SEQUENTIAL_COLORSCALE_OPTIONS, default="Blues", default_reverse=True,
        )
        dist_group_label_map = {}
        if heatmap_group_column and heatmap_show_group_annotation:
            dist_group_values = sorted(meta_df[heatmap_group_column].astype(str).unique())
            dist_group_label_map = dew._render_group_label_renaming_controls(
                f"{WORKSPACE_KEY}_sample_dist_{heatmap_group_column}", dist_group_values,
                label=f"Rename `{heatmap_group_column}` group labels shown on this heatmap (optional)",
            )
        dist_fig = dew._plot_sample_distance_heatmap(
            distance_df, meta_df=meta_df, group_column=heatmap_group_column,
            show_group_annotation=heatmap_show_group_annotation, group_label_map=dist_group_label_map,
            colorscale=dist_colorscale, reverse_colorscale=dist_reverse,
        )
        dist_style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_sample_dist", group_values=None, show_legend_controls=False)
        dew._apply_plot_style(dist_fig, dist_style, default_title="Sample-to-Sample Distance")
        dew._render_plotly_chart(dist_fig)
        dew._render_pdf_export(dist_fig, f"{WORKSPACE_KEY}_sample_dist", "sc_sample_distance_heatmap")

        if heatmap_group_column:
            mismatches = dm.detect_sample_clustering_mismatch(distance_df, meta_df, heatmap_group_column)
            if mismatches:
                mismatch_lines = "\n".join(
                    f"- **{m['sample']}** (`{heatmap_group_column}` = {m['own_group']}) is most similar to "
                    f"**{m['nearest_neighbor']}** (`{heatmap_group_column}` = {m['neighbor_group']})"
                    for m in mismatches
                )
                st.warning(f"⚠️ **{len(mismatches)} sample(s)** cluster with a different group:\n\n{mismatch_lines}")
            else:
                st.success(f"✅ Every sample's closest match shares its own `{heatmap_group_column}` group.")
        dew._render_csv_download(distance_df.reset_index(), "sc_sample_distance_matrix", f"{WORKSPACE_KEY}_sample_dist", expander_label="⬇️ Download Sample Distance Matrix.csv")
    else:
        st.info("Sample distance data isn't available yet -- re-run DESeq2 above.")
    st.markdown("---")

    # --- Model-fit diagnostics ---
    st.subheader("🔧 Model Fit Diagnostics")
    dispersion_df = dew._load_dispersion_estimates(deseq2_out_dir)
    if dispersion_df is not None:
        disp_fig = dew._plot_dispersion(dispersion_df)
        disp_style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_dispersion", group_values=None, show_legend_controls=True)
        dew._apply_plot_style(disp_fig, disp_style, default_title="Dispersion Estimates")
        dew._render_plotly_chart(disp_fig)
        dew._render_pdf_export(disp_fig, f"{WORKSPACE_KEY}_dispersion", "sc_dispersion_estimates")
        fit_assessment = dm.assess_dispersion_fit(dispersion_df)
        (st.warning if fit_assessment["flagged"] else st.success)(fit_assessment["message"])
        dew._render_csv_download(dispersion_df, "sc_dispersion_estimates", f"{WORKSPACE_KEY}_dispersion", expander_label="⬇️ Download Dispersion Estimates.csv")
    else:
        st.info("Dispersion data isn't available yet -- re-run DESeq2 above.")

    size_factor_df = dew._load_size_factors(deseq2_out_dir)
    if size_factor_df is not None:
        st.markdown("##### ⚖️ Size Factors (Normalization QC)")
        st.dataframe(size_factor_df, use_container_width=True, hide_index=True)
        sf_outliers = dm.flag_size_factor_outliers(size_factor_df)
        if sf_outliers["flagged_samples"]:
            outlier_lines = "\n".join(
                f"- **{o['sample']}**: size factor {o['size_factor']} ({o['ratio_to_median']}x the median)"
                for o in sf_outliers["flagged_samples"]
            )
            st.warning(f"⚠️ **{len(sf_outliers['flagged_samples'])} sample(s)** have a notable size factor:\n\n{outlier_lines}")
        else:
            st.success("✅ All samples' size factors are within a reasonable range.")
        dew._render_csv_download(size_factor_df, "sc_size_factors", f"{WORKSPACE_KEY}_size_factors", expander_label="⬇️ Download Size Factors.csv")
    st.markdown("---")

    # --- Per-contrast results ---
    st.subheader("📊 Results by Contrast")
    selected_contrast = st.selectbox("Select a contrast to view:", options=available_contrasts, key=f"{WORKSPACE_KEY}_view_contrast_select")
    results_df = dew._load_contrast_results(deseq2_out_dir, selected_contrast)
    if results_df is not None:
        gene_name_map = _render_sc_gene_id_mapping_panel(project, export_name, cell_type_suffix, counts_df)

        padj_cutoff = st.slider("Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key=f"{WORKSPACE_KEY}_padj_slider")
        lfc_cutoff = st.slider("Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key=f"{WORKSPACE_KEY}_lfc_slider")
        n_sig = ((results_df["padj"] < padj_cutoff) & (results_df["log2FoldChange"].abs() >= lfc_cutoff)).sum()
        st.caption(f"**{n_sig:,}** significant genes at padj < {padj_cutoff} and |log2FC| >= {lfc_cutoff}.")

        regulation_label_map = dew._render_group_label_renaming_controls(
            f"{WORKSPACE_KEY}_regulation_{selected_contrast}",
            ["Up-regulated", "Down-regulated", "Not significant"],
            label="Rename Up-regulated/Down-regulated/Not-significant labels (optional)",
        )

        volcano_fig, volcano_df = dew._plot_volcano(results_df, padj_cutoff, lfc_cutoff, gene_name_map=gene_name_map)
        volcano_style = dew._render_plot_style_controls(
            f"{WORKSPACE_KEY}_volcano_{selected_contrast}",
            group_values=["Up-regulated", "Down-regulated", "Not significant"],
            default_colors={"Up-regulated": "#d62728", "Down-regulated": "#1f77b4", "Not significant": "#7f7f7f"},
        )
        dew._apply_plot_style(volcano_fig, volcano_style, default_title=f"Volcano Plot: {selected_contrast}")
        dew._apply_group_label_renaming(volcano_fig, regulation_label_map)

        labeled_genes, sig_gene_ids = dew._render_gene_label_controls(selected_contrast, volcano_df, gene_name_map=gene_name_map)
        gene_annotations = dew._build_gene_label_annotations(volcano_df, labeled_genes, gene_name_map=gene_name_map)
        if gene_annotations:
            volcano_fig.update_layout(annotations=list(volcano_fig.layout.annotations) + gene_annotations)

        selection_event = st.plotly_chart(
            volcano_fig, use_container_width=True, config=dew._PLOTLY_CHART_CONFIG,
            on_select="rerun", selection_mode=["points"], key=f"{WORKSPACE_KEY}_volcano_chart_{selected_contrast}",
        )
        newly_added = False
        if selection_event and selection_event.selection and selection_event.selection.points:
            for pt in selection_event.selection.points:
                gid = pt.get("customdata")
                if isinstance(gid, (list, tuple)):
                    gid = gid[0] if gid else None
                if gid is None:
                    continue
                gid = str(gid)
                if gid not in labeled_genes:
                    labeled_genes[gid] = dict(dew._DEFAULT_GENE_LABEL_STYLE)
                    newly_added = True
        if newly_added:
            st.rerun()
        dew._render_pdf_export(volcano_fig, f"{WORKSPACE_KEY}_volcano_{selected_contrast}", f"sc_volcano_{selected_contrast}")
        st.markdown("---")

        st.markdown("##### 📉 MA Plot")
        ma_fig = dew._plot_ma(results_df, padj_cutoff, lfc_cutoff, gene_name_map=gene_name_map)
        ma_style = dew._render_plot_style_controls(
            f"{WORKSPACE_KEY}_ma_{selected_contrast}",
            group_values=["Up-regulated", "Down-regulated", "Not significant"],
            default_colors={"Up-regulated": "#d62728", "Down-regulated": "#1f77b4", "Not significant": "#7f7f7f"},
        )
        dew._apply_plot_style(ma_fig, ma_style, default_title=f"MA Plot: {selected_contrast}")
        dew._apply_group_label_renaming(ma_fig, regulation_label_map)
        dew._render_plotly_chart(ma_fig)
        dew._render_pdf_export(ma_fig, f"{WORKSPACE_KEY}_ma_{selected_contrast}", f"sc_ma_plot_{selected_contrast}")
        ma_bias = dm.classify_ma_bias(results_df)
        (st.warning if ma_bias["flagged"] else st.success)(ma_bias["message"])
        st.markdown("---")

        st.markdown("##### 📊 P-value Distribution")
        pvalue_shape = dm.classify_pvalue_histogram_shape(results_df)
        if pvalue_shape["counts"]:
            pval_fig = dew._plot_pvalue_histogram(pvalue_shape)
            pval_style = dew._render_plot_style_controls(
                f"{WORKSPACE_KEY}_pval_{selected_contrast}", group_values=["P-values"],
                default_colors={"P-values": "#636EFA"}, show_legend_controls=False,
            )
            dew._apply_plot_style(pval_fig, pval_style, default_title=f"P-value Distribution: {selected_contrast}")
            dew._render_plotly_chart(pval_fig)
            dew._render_pdf_export(pval_fig, f"{WORKSPACE_KEY}_pval_{selected_contrast}", f"sc_pvalue_histogram_{selected_contrast}")
            if pvalue_shape["shape"] == "healthy":
                st.success(pvalue_shape["message"])
            elif pvalue_shape["shape"] == "flat_uniform":
                st.info(pvalue_shape["message"])
            else:
                st.warning(pvalue_shape["message"])
        else:
            st.info(pvalue_shape["message"])
        st.markdown("---")

        st.markdown("##### 🧬 Top Genes")
        norm_counts_df = dew._load_normalized_counts(deseq2_out_dir)
        if norm_counts_df is not None:
            top_genes_group_columns = [c for c in meta_df.columns if c != "sample"]
            top_genes_group_column, top_genes_group_samples = None, True
            col_n, col_group, col_order = st.columns([1, 1.4, 1.4])
            with col_n:
                n_top_genes = st.number_input("Number of top genes:", min_value=5, max_value=100, value=25, step=5, key=f"{WORKSPACE_KEY}_top_genes_n_{selected_contrast}")
            if top_genes_group_columns:
                with col_group:
                    default_tg_idx = 0
                    if saved_config.get("design_columns") and saved_config["design_columns"][0] in top_genes_group_columns:
                        default_tg_idx = top_genes_group_columns.index(saved_config["design_columns"][0])
                    top_genes_group_column = st.selectbox(
                        "Group/label samples by:", options=top_genes_group_columns, index=default_tg_idx,
                        key=f"{WORKSPACE_KEY}_top_genes_group_column_{selected_contrast}",
                    )
                with col_order:
                    top_genes_group_samples = not st.checkbox(
                        "Ignore grouping (cluster samples purely by expression)", value=False,
                        key=f"{WORKSPACE_KEY}_top_genes_pure_cluster_{selected_contrast}",
                    )
            top_genes_colorscale, top_genes_reverse = dew._render_heatmap_color_controls(
                f"{WORKSPACE_KEY}_top_genes_{selected_contrast}", dew.DIVERGING_COLORSCALE_OPTIONS, default="RdBu", default_reverse=True,
            )
            top_genes_label_map = {}
            if top_genes_group_column and top_genes_group_samples:
                tg_group_values = sorted(meta_df[top_genes_group_column].astype(str).unique())
                top_genes_label_map = dew._render_group_label_renaming_controls(
                    f"{WORKSPACE_KEY}_top_genes_{selected_contrast}_{top_genes_group_column}", tg_group_values,
                    label=f"Rename `{top_genes_group_column}` group labels (optional)",
                )
            top_genes_fig, top_genes_z_df = dew._plot_top_genes_heatmap(
                norm_counts_df, results_df, gene_name_map=gene_name_map, n_top=n_top_genes,
                meta_df=meta_df, group_column=top_genes_group_column, group_samples=top_genes_group_samples,
                group_label_map=top_genes_label_map, colorscale=top_genes_colorscale, reverse_colorscale=top_genes_reverse,
            )
            top_genes_style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_top_genes_{selected_contrast}", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(top_genes_fig, top_genes_style, default_title=f"Top {n_top_genes} Genes: {selected_contrast}")
            dew._render_plotly_chart(top_genes_fig)
            dew._render_pdf_export(top_genes_fig, f"{WORKSPACE_KEY}_top_genes_{selected_contrast}", f"sc_top_genes_heatmap_{selected_contrast}")
            dew._render_csv_download(
                top_genes_z_df.reset_index().rename(columns={"index": "gene_id"}),
                f"sc_top_{n_top_genes}_genes_zscores_{selected_contrast}", f"{WORKSPACE_KEY}_top_genes_{selected_contrast}",
                expander_label=f"⬇️ Download Top {n_top_genes} Genes Z-Scores.csv",
            )
        else:
            st.info("Normalized counts aren't available yet -- re-run DESeq2 above.")
        st.markdown("---")

        annotated_results_df = results_df.copy()
        annotated_results_df["Regulation"] = dm.classify_regulation(annotated_results_df, padj_cutoff, lfc_cutoff)
        sci_notation_threshold = st.select_slider(
            "Switch very small p-values to scientific notation below:",
            options=[1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10], value=1e-4,
            format_func=lambda v: f"{v:.0e}", key=f"{WORKSPACE_KEY}_sci_notation_threshold_{selected_contrast}",
        )
        display_results_df = dew._format_pvalues_for_display(annotated_results_df, columns=("pvalue", "padj"), sci_notation_threshold=sci_notation_threshold)
        if regulation_label_map:
            display_results_df["Regulation"] = display_results_df["Regulation"].map(lambda r: regulation_label_map.get(r, r))
        st.dataframe(display_results_df.head(200), use_container_width=True, hide_index=True)
        dew._render_csv_download(annotated_results_df, f"sc_deseq2_results_{selected_contrast}", f"{WORKSPACE_KEY}_results_{selected_contrast}", expander_label=f"⬇️ Download Full Results ({selected_contrast}).csv")

        with st.expander("🧬 Export for clusterProfiler / Ontology Analysis"):
            export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrast)
            if export_df is not None:
                display_export_df = dew._format_pvalues_for_display(export_df, columns=("pvalue", "padj"), sci_notation_threshold=sci_notation_threshold)
                st.dataframe(display_export_df.head(20), use_container_width=True, hide_index=True)
                dew._render_csv_download(export_df, f"sc_clusterprofiler_ranked_{selected_contrast}", f"{WORKSPACE_KEY}_cp_{selected_contrast}", expander_label="⬇️ Download clusterProfiler-Ready Gene List.csv")
    st.markdown("---")

    # --- Venn diagram across contrasts ---
    st.subheader("🔵 Overlap Between Contrasts (Venn Diagram)")
    if len(available_contrasts) < 2:
        st.info("Run at least 2 contrasts to see an overlap diagram.")
    else:
        venn_contrasts = st.multiselect(
            "Select 2 or 3 contrasts to compare:", options=available_contrasts,
            default=available_contrasts[:min(2, len(available_contrasts))], key=f"{WORKSPACE_KEY}_venn_contrast_select",
        )
        if len(venn_contrasts) not in (2, 3):
            st.info("Please select exactly 2 or 3 contrasts.")
        else:
            venn_padj = st.slider("Significance threshold for overlap:", 0.01, 0.20, 0.05, step=0.01, key=f"{WORKSPACE_KEY}_venn_padj_slider")
            with st.expander("✏️ Rename contrast labels (optional)"):
                venn_labels = {name: st.text_input(f"Label for `{name}`:", value=name, key=f"{WORKSPACE_KEY}_venn_label_{name}") for name in venn_contrasts}

            gene_sets = dm.get_significant_gene_sets(deseq2_out_dir, venn_contrasts, padj_threshold=venn_padj)
            seen_labels, display_gene_sets, orig_to_display_label = set(), {}, {}
            for name, gene_set in gene_sets.items():
                custom_label = (venn_labels.get(name, "") or "").strip() or name
                if custom_label in seen_labels:
                    st.warning(f"⚠️ The label \"{custom_label}\" is used more than once -- keeping \"{name}\" instead.")
                    custom_label = name
                seen_labels.add(custom_label)
                display_gene_sets[custom_label] = gene_set
                orig_to_display_label[name] = custom_label

            display_names = list(display_gene_sets.keys())
            counts, regions = dm.compute_venn_regions(display_gene_sets)
            venn_fig = dew._plot_venn(counts, display_names)
            venn_style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_venn", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(venn_fig, venn_style, default_title="Overlap Between Contrasts")
            dew._render_plotly_chart(venn_fig)
            dew._render_pdf_export(venn_fig, f"{WORKSPACE_KEY}_venn", "sc_venn_diagram_overlap")

            with st.expander("View overlapping gene lists"):
                for label, gene_set in regions.items():
                    st.markdown(f"**{label}** ({len(gene_set)} genes)")
                    if gene_set:
                        st.text(", ".join(sorted(gene_set)[:50]) + (" ..." if len(gene_set) > 50 else ""))

            venn_export_df = dm.build_venn_export(deseq2_out_dir, venn_contrasts, regions, display_label_map=orig_to_display_label)
            if not venn_export_df.empty:
                dew._render_csv_download(venn_export_df, "sc_venn_overlap_full_results", f"{WORKSPACE_KEY}_venn_full", expander_label="⬇️ Download Full Venn Results.csv")

    st.markdown("---")
    st.success(f"🎉 Pseudobulk export `{export_name}{cell_type_suffix}` has completed DESeq2 analysis.")
    if st.button("➡️ Proceed to Ontology Analysis", type="primary", key=f"{WORKSPACE_KEY}_proceed_ontology_btn"):
        st.session_state["nav_request"] = "🌋 SC Ontology Analysis"
        st.rerun()


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------
def render():
    st.title("🌋 SC DESeq2")
    st.markdown(
        "Runs the Bulk RNA-Seq pipeline's own DESeq2 infrastructure on a "
        "single-cell pseudobulk export -- the same PCA, volcano/MA plots, "
        "QC diagnostics, and results tables you'd get in the Bulk pipeline."
    )
    st.markdown("---")

    project, export_name, counts_df, meta_df, cell_type_suffix = render_project_export_and_cell_type_selector()
    if not (project and export_name and counts_df is not None and meta_df is not None):
        return
    st.markdown("---")

    _render_sc_analysis_settings(project, export_name, cell_type_suffix, counts_df, meta_df)
    st.markdown("---")

    saved_config = scpm.get_sc_deseq2_config(project, export_name + cell_type_suffix) or {}
    _render_step5_results(project, export_name, cell_type_suffix, meta_df, saved_config, counts_df)
