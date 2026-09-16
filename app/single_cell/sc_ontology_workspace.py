"""
single_cell/sc_ontology_workspace.py

Runs the Bulk RNA-Seq pipeline's own Ontology Analysis infrastructure
(ontology_manager.py) on a single-cell pseudobulk export's DESeq2 results,
reusing ontology_workspace.py's (ow) plotting/styling code DIRECTLY -- every
ow._plot_*, ow._render_*_control(s), and the shared dew._render_plot_style_controls/
_apply_plot_style/_render_csv_download/_render_pdf_export/_render_heatmap_color_controls
calls ow itself already relies on, all take plain DataFrames/paths as arguments
and have NO dependency on project_manager (pm) -- so this module calls them
verbatim rather than reimplementing any plotting logic.

This is the companion to sc_deseq2_workspace.py. Both share project/export/
cell-type selection via sc_deseq2_workspace's session_state keys (see
get_shared_selection()) -- if the user hasn't visited the SC DESeq2 page
yet this session, this page falls back to rendering its own copy of the
Step 0/0.5 selector.

v1 SCOPE NOTE (matching sc_deseq2_ontology_workspace.py's original
docstring): ORA and GSEA only -- compareCluster (multi-contrast side-by-
side comparison) and eggNOG-mapper mode (for non-model organisms with no
curated annotation package) are NOT ported here. Both are real, separate
features `ow` supports; adding them later would follow the identical
reuse pattern already established in this file.
"""
import os

import pandas as pd
import streamlit as st

import sc_project_manager as scpm
import sc_deseq2_workspace as sdw
import deseq2_manager as dm
import gene_id_mapper as gim
import reference_manager as rm
import ontology_manager as om
import differential_expression_workspace as dew
import ontology_workspace as ow

WORKSPACE_KEY = "sc_ontology"


def _ensure_project_export_selection():
    """
    Reuses whatever project/export/cell-type was already picked on the SC
    DESeq2 page this session (see sc_deseq2_workspace.get_shared_selection).
    If nothing has been picked yet, falls back to rendering that same
    Step 0/0.5 selector directly on this page, so this page is still
    usable standalone.
    """
    project, export_name, cell_type_suffix = sdw.get_shared_selection()
    if project and export_name:
        export_dir = scpm.downstream_pseudobulk_export_dir(project, export_name)
        counts_path = os.path.join(export_dir, "pseudobulk_counts.csv")
        metadata_path = os.path.join(export_dir, "pseudobulk_metadata.csv")
        if os.path.isfile(counts_path) and os.path.isfile(metadata_path):
            st.success(f"✅ Using project `{project}`, export `{export_name}{cell_type_suffix}` (selected on the SC DESeq2 page).")
            if st.button("🔄 Change project/export/cell type", key=f"{WORKSPACE_KEY}_change_selection_btn"):
                for key in ("_sc_de_bridge_project", "_sc_de_bridge_export", "_sc_de_bridge_cell_type_suffix"):
                    st.session_state.pop(key, None)
                st.rerun()
            return project, export_name, cell_type_suffix

    st.info("ℹ️ No project/export selected yet this session -- choose one below.")
    project, export_name, counts_df, meta_df, cell_type_suffix = sdw.render_project_export_and_cell_type_selector()
    if not (project and export_name):
        return None, None, ""
    return project, export_name, cell_type_suffix


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------
def render():
    st.title("🌋 SC Ontology Analysis")
    st.markdown(
        "This workspace takes your SC DESeq2 results and asks: **what do "
        "these genes have in common, biologically?** It runs enrichment "
        "analysis against **Gene Ontology (GO)**, **KEGG pathways**, and "
        "**Reactome pathways** -- reusing the exact same infrastructure "
        "the Bulk RNA-Seq pipeline's Ontology Analysis page uses."
    )
    st.markdown("---")

    project, export_name, cell_type_suffix = _ensure_project_export_selection()
    if not project:
        return
    full_export_name = export_name + cell_type_suffix
    st.markdown("---")

    deseq2_out_dir = scpm.sc_deseq2_output_dir(project, full_export_name)
    available_contrasts = dm.list_contrast_results(deseq2_out_dir)
    if not available_contrasts:
        st.warning(
            "⚠️ This pseudobulk export doesn't have DESeq2 results yet. "
            "Complete that step on the SC DESeq2 page first."
        )
        if st.button("⬅️ Go to SC DESeq2", key=f"{WORKSPACE_KEY}_gate_back_btn"):
            st.session_state["nav_request"] = "🌋 SC DESeq2"
            st.rerun()
        return

    if not om.clusterprofiler_tools_available():
        st.error(
            "⚠️ Rscript was not found on this system. R with clusterProfiler "
            "(and ReactomePA, for Reactome pathways) needs to be installed "
            "before this step can run."
        )
        return

    # --- Species selection (v1: always manual confirmation -- no
    # sc-side equivalent of pm.get_reference_choice()'s preset/custom
    # auto-detection has been verified to exist, so we always ask,
    # reusing the sc_ontology species-override helpers already built
    # in sc_project_manager.py). ---
    species_choices = gim.orgdb_species_choices(rm)
    species_keys = list(species_choices.keys())
    previous_species = scpm.get_sc_ontology_species(project, full_export_name)
    default_index = species_keys.index(previous_species) if previous_species in species_keys else 0
    st.markdown("**Which organism is this data from?**")
    picked_species = st.selectbox(
        "Species (annotation database):", options=species_keys, index=default_index,
        format_func=lambda k: species_choices[k], key=f"{WORKSPACE_KEY}_species_select",
    )
    scpm.save_sc_ontology_species(project, full_export_name, picked_species)
    orgdb_package = gim.ORGDB_PACKAGES.get(picked_species)
    species_label = species_choices.get(picked_species, picked_species)
    kegg_organism = om.get_kegg_organism_code(picked_species)
    reactome_organism = om.get_reactome_organism_name(picked_species)
    st.session_state["_ontology_orgdb_package_for_gosemsim_check"] = orgdb_package

    db_availability = om.databases_available_for_species(picked_species)
    st.success(f"✅ Organism: **{species_label}** (annotation package: `{orgdb_package}`).")
    unavailable = [db for db, ok in db_availability.items() if not ok]
    if unavailable:
        st.info(f"ℹ️ {' and '.join(unavailable)} {'is' if len(unavailable) == 1 else 'are'} not available for {species_label} and will be grayed out below.")

    ontology_out_dir = scpm.sc_ontology_output_dir(project, full_export_name)
    ontology_work_dir = scpm.sc_ontology_work_dir(project, full_export_name)

    st.markdown("---")
    st.header("Step 1: Choose Your Analysis Approach")
    with st.expander("ℹ️ ORA vs. GSEA -- which should I use? (click to learn more)"):
        st.markdown(
            "**🎯 ORA (Over-Representation Analysis)** -- takes a FIXED "
            "list of your significant genes and asks whether they're "
            "enriched for any GO term/pathway more than chance.\n\n"
            "**📈 GSEA (Gene Set Enrichment Analysis)** -- uses EVERY "
            "tested gene, ranked from most up- to most down-regulated, "
            "and asks whether a term's genes cluster toward one end of "
            "that ranking."
        )
    analysis_approach = st.radio(
        "Which approach would you like to use?",
        ["ORA (single contrast)", "GSEA (single contrast)"], key=f"{WORKSPACE_KEY}_analysis_approach_radio",
    )

    st.markdown("---")
    st.header("Step 2: Choose Your Contrast and Databases")
    selected_contrast = st.selectbox("Which contrast?", options=available_contrasts, key=f"{WORKSPACE_KEY}_single_contrast_select")
    selected_contrasts = [selected_contrast]

    db_col1, db_col2, db_col3 = st.columns(3)
    run_go = db_col1.checkbox("GO (Gene Ontology)", value=True, key=f"{WORKSPACE_KEY}_run_go", disabled=not db_availability["GO"])
    run_kegg = db_col2.checkbox("KEGG pathways", value=db_availability["KEGG"], key=f"{WORKSPACE_KEY}_run_kegg", disabled=not db_availability["KEGG"])
    run_reactome = db_col3.checkbox("Reactome pathways", value=db_availability["Reactome"], key=f"{WORKSPACE_KEY}_run_reactome", disabled=not db_availability["Reactome"])

    go_ontology = "ALL"
    simplify_go, simplify_measure, simplify_cutoff = False, om.DEFAULT_SIMPLIFY_MEASURE, om.DEFAULT_SIMPLIFY_CUTOFF
    if run_go:
        go_ontology_label = st.selectbox(
            "GO sub-ontology:", options=list(om.GO_ONTOLOGY_OPTIONS.keys()), key=f"{WORKSPACE_KEY}_go_subontology_select",
        )
        go_ontology = om.GO_ONTOLOGY_OPTIONS[go_ontology_label]
        simplify_go, simplify_measure, simplify_cutoff = ow._render_simplify_control(f"{WORKSPACE_KEY}_step2", go_ontology)

    if not (run_go or run_kegg or run_reactome):
        st.warning("⚠️ Select at least one database to continue.")
        return

    min_gs_size, max_gs_size = ow._render_gene_set_size_controls(f"{WORKSPACE_KEY}_step2")

    needs_threshold = analysis_approach.startswith("ORA")
    padj_threshold, lfc_threshold, split_by_direction = 0.05, 1.0, True
    if needs_threshold:
        col_p, col_l = st.columns(2)
        with col_p:
            padj_threshold = st.slider("Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key=f"{WORKSPACE_KEY}_padj_slider")
        with col_l:
            lfc_threshold = st.slider("Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key=f"{WORKSPACE_KEY}_lfc_slider")
        split_by_direction = ow._render_direction_split_control(f"{WORKSPACE_KEY}_step2")

    first_export_path = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
    from_type = "SYMBOL"
    if first_export_path is not None and not first_export_path.empty:
        detection = gim.detect_id_type(first_export_path["gene_id"].astype(str).tolist())
        from_type = detection["detected_type"]
        st.caption(f"Auto-detected gene ID type: **{from_type}** ({detection['match_fraction'] * 100:.0f}% matched).")

    st.markdown("---")
    st.header("Step 3: Run the Analysis")

    analysis_type = "ora" if analysis_approach.startswith("ORA") else "gsea"
    contrast_key_for_state = selected_contrasts[0]
    output_dir = os.path.join(ontology_out_dir, analysis_type, contrast_key_for_state)
    work_dir = os.path.join(ontology_work_dir, analysis_type, contrast_key_for_state)

    results_already_exist = any(
        om.enrichment_result_exists(output_dir, analysis_type, db)
        or (analysis_type == "ora" and om.available_ora_directions(output_dir, db))
        for db in ["GO", "KEGG", "Reactome"]
    )
    run_label = "🔄 Re-run Analysis" if results_already_exist else "🚀 Run Enrichment Analysis"
    if results_already_exist:
        st.success("✅ Results already exist for this exact configuration.")

    if st.button(run_label, key=f"{WORKSPACE_KEY}_run_btn_{analysis_type}"):
        with st.spinner(f"Running {analysis_type.upper()}... this may take a few minutes."):
            export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
            os.makedirs(work_dir, exist_ok=True)
            input_path = os.path.join(work_dir, "input_gene_list.csv")
            export_df.to_csv(input_path, index=False)
            if analysis_type == "ora":
                success, log = om.run_ora_analysis(
                    input_path, output_dir, work_dir, orgdb_package, from_type,
                    padj_threshold, lfc_threshold, run_go, go_ontology,
                    run_kegg, kegg_organism, run_reactome, reactome_organism,
                    min_gs_size=min_gs_size, max_gs_size=max_gs_size,
                    simplify_go=simplify_go, simplify_measure=simplify_measure,
                    simplify_cutoff=simplify_cutoff, split_by_direction=split_by_direction,
                )
            else:
                success, log = om.run_gsea_analysis(
                    input_path, output_dir, work_dir, orgdb_package, from_type,
                    run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
                    min_gs_size=min_gs_size, max_gs_size=max_gs_size,
                    simplify_go=simplify_go, simplify_measure=simplify_measure,
                    simplify_cutoff=simplify_cutoff,
                )
        if not success:
            st.error("Analysis failed. Details below:")
            st.code(log)
            return
        st.success("✅ Analysis completed successfully.")
        if simplify_go and run_go:
            simplify_outcomes = om.parse_simplify_outcomes_from_log(log)
            if simplify_outcomes:
                om.save_simplify_status(output_dir, simplify_outcomes)
        with st.expander("View run log"):
            st.code(log)
        results_already_exist = True

    if not results_already_exist:
        return

    st.markdown("---")
    st.header("Step 4: Explore Your Results")

    gene_name_map = {}
    auto_map_path = scpm.sc_gene_symbol_map_path(project, full_export_name)
    if os.path.exists(auto_map_path):
        try:
            auto_map_df = pd.read_csv(auto_map_path)
            gene_name_map = dict(zip(auto_map_df["gene_id"].astype(str), auto_map_df["gene_name"].astype(str)))
        except Exception:
            gene_name_map = {}

    if analysis_type == "ora":
        available_dbs = om.list_available_ora_databases(output_dir)
        if not available_dbs:
            st.info("No results were found for any database (or zero terms were returned).")
            return

        gene_fc_map = ow._get_gene_fc_map_for_contrast(deseq2_out_dir, selected_contrasts[0], gene_name_map=gene_name_map)
        combined_view_direction_by_db = {}

        for db in available_dbs:
            directions = om.available_ora_directions(output_dir, db)
            if "up" in directions and "down" in directions:
                st.markdown(f"### {db}")
                st.caption("Analyzed separately for up- and down-regulated genes.")
                ow._render_single_database_results(
                    output_dir, "ora", db, f"{WORKSPACE_KEY}_ora_{db}_up", padj_threshold,
                    gene_fc_map=gene_fc_map, direction="up", section_label=f"{db} results — Up-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                st.markdown("---")
                ow._render_single_database_results(
                    output_dir, "ora", db, f"{WORKSPACE_KEY}_ora_{db}_down", padj_threshold,
                    gene_fc_map=gene_fc_map, direction="down", section_label=f"{db} results — Down-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                st.markdown("---")
                ow._render_ora_up_down_comparison(output_dir, db, f"{WORKSPACE_KEY}_ora_{db}")
                combined_view_direction_by_db[db] = "up"
            elif "combined" in directions:
                ow._render_single_database_results(
                    output_dir, "ora", db, f"{WORKSPACE_KEY}_ora_{db}", padj_threshold,
                    gene_fc_map=gene_fc_map, direction=None, orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                combined_view_direction_by_db[db] = "combined"
            else:
                only_direction = directions[0]
                ow._render_single_database_results(
                    output_dir, "ora", db, f"{WORKSPACE_KEY}_ora_{db}_{only_direction}", padj_threshold,
                    gene_fc_map=gene_fc_map, direction=only_direction,
                    section_label=f"{db} results — {'Up' if only_direction == 'up' else 'Down'}-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                combined_view_direction_by_db[db] = only_direction
            st.markdown("---")

        if len(available_dbs) >= 2:
            st.subheader("🔀 Combined View Across Databases")
            show_combined = st.checkbox("📊 Show combined dot plot across selected databases", key=f"{WORKSPACE_KEY}_ora_combined_checkbox")
            if show_combined:
                combine_dbs = st.multiselect("Which databases to combine?", options=available_dbs, default=available_dbs, key=f"{WORKSPACE_KEY}_ora_combine_db_select")
                n_top_combined = st.slider("Top terms per database:", min_value=3, max_value=20, value=10, step=1, key=f"{WORKSPACE_KEY}_ora_combined_n_top")
                combined_colorscale, combined_reverse = dew._render_heatmap_color_controls(
                    f"{WORKSPACE_KEY}_ora_combined", ow.SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
                )
                frames = []
                for db in combine_dbs:
                    direction = combined_view_direction_by_db.get(db)
                    read_direction = None if direction == "combined" else direction
                    df = om.read_enrichment_result(output_dir, "ora", db, direction=read_direction)
                    if df is not None and not df.empty:
                        df = df.sort_values("p.adjust").head(n_top_combined).copy()
                        df["Database"] = db
                        frames.append(df)
                combined_df = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
                combined_df_categorized = om.derive_category_column(combined_df)
                selected_categories, category_label_map = ow._render_category_filter_and_label_controls(combined_df_categorized, f"{WORKSPACE_KEY}_ora_combined")
                combined_df_filtered = combined_df_categorized[combined_df_categorized["Category"].isin(selected_categories)]
                combined_label_map = ow._render_term_label_editor(combined_df_filtered, f"{WORKSPACE_KEY}_ora_combined_labels")
                fig, category_colors = (None, {}) if combined_df_filtered.empty else ow._plot_combined_multi_database(
                    combined_df_filtered, label_map=combined_label_map, category_label_map=category_label_map,
                    colorscale=combined_colorscale, reverse_colorscale=combined_reverse,
                )
                if fig is None:
                    st.info("No data available to combine -- check at least one category above.")
                else:
                    style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_ora_combined", group_values=None, show_legend_controls=False)
                    dew._apply_plot_style(fig, style, default_title="Combined Enrichment -- GO / KEGG / Reactome")
                    dew._render_plotly_chart(fig)
                    dew._render_pdf_export(fig, f"{WORKSPACE_KEY}_ora_combined", "sc_combined_multi_database")
                    dew._render_csv_download(combined_df_filtered, "sc_combined_ora_results", f"{WORKSPACE_KEY}_ora_combined_table", expander_label="⬇️ Download Combined Results.csv")

    else:  # gsea
        available_dbs = om.list_available_databases(output_dir, analysis_type)
        if not available_dbs:
            st.info("No results were found for any database (or zero terms were returned).")
            return
        gene_fc_map = ow._get_gene_fc_map_for_contrast(deseq2_out_dir, selected_contrasts[0], gene_name_map=gene_name_map)

        for db in available_dbs:
            ow._render_single_database_results(
                output_dir, analysis_type, db, f"{WORKSPACE_KEY}_{analysis_type}_{db}", padj_threshold,
                gene_fc_map=gene_fc_map, orgdb_package=orgdb_package, go_ontology=go_ontology,
            )
            st.markdown("---")

        if len(available_dbs) >= 2:
            st.subheader("🔀 Combined View Across Databases")
            show_combined = st.checkbox("📊 Show combined dot plot across selected databases", key=f"{WORKSPACE_KEY}_{analysis_type}_combined_checkbox")
            if show_combined:
                combine_dbs = st.multiselect("Which databases to combine?", options=available_dbs, default=available_dbs, key=f"{WORKSPACE_KEY}_{analysis_type}_combine_db_select")
                n_top_combined = st.slider("Top terms per database:", min_value=3, max_value=20, value=10, step=1, key=f"{WORKSPACE_KEY}_{analysis_type}_combined_n_top")
                combined_colorscale, combined_reverse = dew._render_heatmap_color_controls(
                    f"{WORKSPACE_KEY}_{analysis_type}_combined", ow.SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
                )
                combined_df = om.build_combined_results(output_dir, analysis_type, combine_dbs, n_top_per_db=n_top_combined)
                combined_df_categorized = om.derive_category_column(combined_df)
                selected_categories, category_label_map = ow._render_category_filter_and_label_controls(combined_df_categorized, f"{WORKSPACE_KEY}_{analysis_type}_combined")
                combined_df_filtered = combined_df_categorized[combined_df_categorized["Category"].isin(selected_categories)]
                combined_label_map = ow._render_term_label_editor(combined_df_filtered, f"{WORKSPACE_KEY}_{analysis_type}_combined_labels")
                fig, category_colors = (None, {}) if combined_df_filtered.empty else ow._plot_combined_multi_database(
                    combined_df_filtered, label_map=combined_label_map, category_label_map=category_label_map,
                    colorscale=combined_colorscale, reverse_colorscale=combined_reverse,
                )
                if fig is None:
                    st.info("No data available to combine -- check at least one category above.")
                else:
                    style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_{analysis_type}_combined", group_values=None, show_legend_controls=False)
                    dew._apply_plot_style(fig, style, default_title="Combined Enrichment -- GO / KEGG / Reactome")
                    dew._render_plotly_chart(fig)
                    dew._render_pdf_export(fig, f"{WORKSPACE_KEY}_{analysis_type}_combined", "sc_combined_multi_database")
                    dew._render_csv_download(combined_df_filtered, f"sc_combined_{analysis_type}_results", f"{WORKSPACE_KEY}_{analysis_type}_combined_table", expander_label="⬇️ Download Combined Results.csv")

    st.markdown("---")
    st.success(f"🎉 Pseudobulk export `{export_name}{cell_type_suffix}` has ontology enrichment results ready to explore above.")
