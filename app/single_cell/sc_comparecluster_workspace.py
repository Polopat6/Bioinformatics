"""
single_cell/sc_comparecluster_workspace.py

Runs the Bulk RNA-Seq pipeline's own compareCluster infrastructure
(ontology_manager.run_compare_cluster_analysis) on MULTIPLE contrasts from
a single-cell pseudobulk export's DESeq2 results at once, for direct
side-by-side comparison -- reusing ontology_workspace.py's (ow) own
_plot_compare_cluster_dotplot and the same shared dew.*
styling/export helpers already used by sc_deseq2_workspace.py and
sc_ontology_workspace.py.

Kept as its own standalone script/page (rather than folded into
sc_ontology_workspace.py) per explicit request, mirroring how compareCluster
is functionally distinct from single-contrast ORA/GSEA in the bulk
pipeline (it compares MULTIPLE contrasts side-by-side, has no direction-
split option, and has no GO-simplification control -- matching the bulk
compareCluster branch's own scope exactly, see ontology_workspace.py's
render()).

Shares project/export/cell-type selection with sc_deseq2_workspace.py via
that module's session_state keys (get_shared_selection()), same as
sc_ontology_workspace.py -- so a project/export picked on either of the
other two SC pages carries over here automatically.
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

WORKSPACE_KEY = "sc_comparecluster"


def _ensure_project_export_selection():
    """
    Same fallback pattern as sc_ontology_workspace._ensure_project_export_selection:
    reuse whatever was already picked on the SC DESeq2 page this session,
    or render that same Step 0/0.5 selector directly if nothing has been
    picked yet, so this page works standalone too.
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


def render():
    st.title("🔀 SC compareCluster")
    st.markdown(
        "Runs GO/KEGG/Reactome enrichment **independently across multiple "
        "contrasts at once**, for direct side-by-side comparison -- "
        "reusing the Bulk RNA-Seq pipeline's own compareCluster "
        "infrastructure and plotting code."
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

    if len(available_contrasts) < 2:
        st.warning(
            "⚠️ compareCluster needs at least 2 contrasts to compare. "
            "Run DESeq2 with 2+ contrasts on the SC DESeq2 page first."
        )
        return

    if not om.clusterprofiler_tools_available():
        st.error(
            "⚠️ Rscript was not found on this system. R with clusterProfiler "
            "(and ReactomePA, for Reactome pathways) needs to be installed "
            "before this step can run."
        )
        return

    # --- Species selection (same v1 pattern as sc_ontology_workspace.py:
    # always ask directly, no sc-side preset/custom auto-detection has
    # been verified to exist yet). ---
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

    db_availability = om.databases_available_for_species(picked_species)
    st.success(f"✅ Organism: **{species_label}** (annotation package: `{orgdb_package}`).")
    unavailable = [db for db, ok in db_availability.items() if not ok]
    if unavailable:
        st.info(f"ℹ️ {' and '.join(unavailable)} {'is' if len(unavailable) == 1 else 'are'} not available for {species_label} and will be grayed out below.")

    ontology_out_dir = scpm.sc_ontology_output_dir(project, full_export_name)
    ontology_work_dir = scpm.sc_ontology_work_dir(project, full_export_name)

    st.markdown("---")
    st.header("Step 1: Choose Your Contrasts and Databases")

    selected_contrasts = st.multiselect(
        "Which contrasts should be compared?", options=available_contrasts,
        default=available_contrasts[:min(2, len(available_contrasts))], key=f"{WORKSPACE_KEY}_contrasts_select",
    )
    if len(selected_contrasts) < 2:
        st.warning("⚠️ Select at least 2 contrasts to run compareCluster.")
        return

    db_col1, db_col2, db_col3 = st.columns(3)
    run_go = db_col1.checkbox("GO (Gene Ontology)", value=True, key=f"{WORKSPACE_KEY}_run_go", disabled=not db_availability["GO"])
    run_kegg = db_col2.checkbox("KEGG pathways", value=db_availability["KEGG"], key=f"{WORKSPACE_KEY}_run_kegg", disabled=not db_availability["KEGG"])
    run_reactome = db_col3.checkbox("Reactome pathways", value=db_availability["Reactome"], key=f"{WORKSPACE_KEY}_run_reactome", disabled=not db_availability["Reactome"])

    go_ontology = "ALL"
    if run_go:
        go_ontology_label = st.selectbox(
            "GO sub-ontology:", options=list(om.GO_ONTOLOGY_OPTIONS.keys()), key=f"{WORKSPACE_KEY}_go_subontology_select",
        )
        go_ontology = om.GO_ONTOLOGY_OPTIONS[go_ontology_label]
        st.caption(
            "ℹ️ GO term simplification isn't offered for compareCluster "
            "(matches the Bulk pipeline's own behavior) -- results below "
            "show the full, un-simplified GO term set."
        )

    if not (run_go or run_kegg or run_reactome):
        st.warning("⚠️ Select at least one database to continue.")
        return

    min_gs_size, max_gs_size = ow._render_gene_set_size_controls(f"{WORKSPACE_KEY}_step1")

    with st.expander("ℹ️ Why do I need to set a threshold here? (click to learn more)"):
        st.markdown(
            "These thresholds determine which genes from each contrast's "
            "Differential Expression results count as significant for "
            "the comparison below."
        )
    col_p, col_l = st.columns(2)
    with col_p:
        padj_threshold = st.slider(
            "Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key=f"{WORKSPACE_KEY}_padj_slider",
        )
    with col_l:
        lfc_threshold = st.slider(
            "Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key=f"{WORKSPACE_KEY}_lfc_slider",
        )

    first_export_path = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
    from_type = "SYMBOL"
    if first_export_path is not None and not first_export_path.empty:
        detection = gim.detect_id_type(first_export_path["gene_id"].astype(str).tolist())
        from_type = detection["detected_type"]
        st.caption(
            f"Auto-detected gene ID type: **{from_type}** "
            f"({detection['match_fraction'] * 100:.0f}% matched -- e.g. "
            f"`{'`, `'.join(detection['example_ids'][:3])}`)."
        )

    st.markdown("---")
    st.header("Step 2: Run the Analysis")

    contrast_key_for_state = "_".join(sorted(selected_contrasts))
    output_dir = os.path.join(ontology_out_dir, "compareCluster", contrast_key_for_state)
    work_dir = os.path.join(ontology_work_dir, "compareCluster", contrast_key_for_state)

    results_already_exist = any(
        om.enrichment_result_exists(output_dir, "compareCluster", db) for db in ["GO", "KEGG", "Reactome"]
    )
    run_label = "🔄 Re-run Analysis" if results_already_exist else "🚀 Run compareCluster Analysis"
    if results_already_exist:
        st.success("✅ Results already exist for this exact configuration.")

    if st.button(run_label, key=f"{WORKSPACE_KEY}_run_btn"):
        with st.spinner("Running compareCluster... this may take a few minutes."):
            os.makedirs(work_dir, exist_ok=True)
            contrast_input_paths = {}
            for c in selected_contrasts:
                export_df = dm.build_clusterprofiler_export(deseq2_out_dir, c)
                input_path = os.path.join(work_dir, f"input_{c}.csv")
                export_df.to_csv(input_path, index=False)
                contrast_input_paths[c] = input_path
            success, log = om.run_compare_cluster_analysis(
                contrast_input_paths, output_dir, work_dir, orgdb_package, from_type,
                padj_threshold, lfc_threshold, run_go, go_ontology,
                run_kegg, kegg_organism, run_reactome, reactome_organism,
                min_gs_size=min_gs_size, max_gs_size=max_gs_size,
            )
        if not success:
            st.error("Analysis failed. Details below:")
            st.code(log)
            return
        st.success("✅ Analysis completed successfully.")
        with st.expander("View run log"):
            st.code(log)
        results_already_exist = True

    if not results_already_exist:
        return

    st.markdown("---")
    st.header("Step 3: Explore Your Results")

    available_dbs = [db for db in ["GO", "KEGG", "Reactome"] if om.enrichment_result_exists(output_dir, "compareCluster", db)]
    if not available_dbs:
        st.info("No compareCluster results were found (or zero terms were returned for any database).")
        return

    for db in available_dbs:
        cc_df = om.read_enrichment_result(output_dir, "compareCluster", db)
        st.markdown(f"#### {db} -- compareCluster results")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "One column per contrast you compared, one row per "
                "enriched term -- dot size shows gene count, dot color "
                "shows significance."
            )
        cc_colorscale, cc_reverse = dew._render_heatmap_color_controls(
            f"{WORKSPACE_KEY}_{db}", ow.SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
        )
        n_top_cc = st.slider(
            f"Top terms per contrast ({db}):", min_value=3, max_value=30, value=10, step=1,
            key=f"{WORKSPACE_KEY}_{db}_n_top",
            help="How many of the most significant terms are shown PER contrast column.",
        )
        cc_label_map = ow._render_term_label_editor(
            cc_df.sort_values("p.adjust").groupby("Cluster").head(n_top_cc).drop_duplicates("ID"),
            f"{WORKSPACE_KEY}_{db}_labels",
        )
        fig = ow._plot_compare_cluster_dotplot(
            cc_df, n_top_per_cluster=n_top_cc, label_map=cc_label_map,
            colorscale=cc_colorscale, reverse_colorscale=cc_reverse,
        )
        if fig is None:
            st.info("No data available to plot.")
        else:
            group_values = sorted(cc_df["Cluster"].astype(str).unique())
            cc_cluster_label_map = dew._render_group_label_renaming_controls(
                f"{WORKSPACE_KEY}_{db}", group_values, label=f"Rename contrast labels shown for {db} (optional)",
            )
            style = dew._render_plot_style_controls(f"{WORKSPACE_KEY}_{db}", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(fig, style, default_title=f"{db} -- compareCluster")
            if cc_cluster_label_map:
                fig.update_xaxes(
                    tickvals=group_values, ticktext=[cc_cluster_label_map.get(v, v) for v in group_values],
                )
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{WORKSPACE_KEY}_{db}", f"sc_{db}_compareCluster_dotplot")

        with st.expander(f"📄 Full {db} compareCluster results table"):
            st.dataframe(cc_df.sort_values("p.adjust"), use_container_width=True, hide_index=True)
            dew._render_csv_download(
                cc_df, f"sc_compareCluster_{db}_results", f"{WORKSPACE_KEY}_{db}_table",
                expander_label=f"⬇️ Download Full {db} compareCluster Results.csv",
            )
        st.markdown("---")

    st.success(f"🎉 Pseudobulk export `{export_name}{cell_type_suffix}` has compareCluster results ready to explore above.")
