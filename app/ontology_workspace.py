"""
ontology_workspace.py

Ontology Analysis workspace.

This picks up where differential_expression_workspace.py leaves off: it
reuses the same active project (per-contrast DESeq2 results already on
disk) rather than requiring anything to be re-uploaded. The clusterProfiler-
ready gene list exports built in the Differential Expression workspace
(see deseq2_manager.build_clusterprofiler_export) are the expected input
for whatever GO/KEGG enrichment logic lands here.

Design goal: same as the other workspaces -- assume the user has little
to no bioinformatics/statistics background, explain each step in plain
language, while still exposing real statistical choices rather than
hiding them.

This module is fully self-contained. All Ontology Analysis development
should happen here -- editing this file has zero effect on the other
workspace modules.

--- Status ---
Placeholder/stub only. The actual enrichment workflow (e.g. clusterProfiler's
enrichGO/enrichKEGG or gseGO/gseKEGG, organism annotation package selection,
term/pathway visualization) has not been built yet -- render() currently
only reuses the shared project selector and reports the workspace as not
yet available, so app.py has somewhere valid to route to via the
Differential Expression workspace's "Proceed to Ontology Analysis" button
without erroring.
"""

import streamlit as st

import project_manager as pm

# Reuse the same workspace_key as the other Bulk RNA-Seq pipeline stages
# so the active project selection is shared automatically across pages.
WORKSPACE_KEY = "bulk_rnaseq"


def render():
    st.title("🧬 Ontology Analysis")
    st.markdown(
        "This workspace will take your significant gene lists from "
        "Differential Expression and run gene ontology (GO) / pathway "
        "enrichment analysis, to help interpret what those genes have "
        "in common biologically. **No bioinformatics experience "
        "required** -- follow the steps below in order."
    )
    st.markdown("---")

    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)
    if not project:
        return

    st.info(
        "🚧 This workspace is under construction. GO/pathway enrichment "
        "(e.g. clusterProfiler's enrichGO/enrichKEGG or gseGO/gseKEGG) "
        "isn't wired up yet -- check back soon."
    )

    if not pm.has_completed_step(project, "deseq2_complete"):
        st.warning(
            "⚠️ This project doesn't have Differential Expression results "
            "yet. Complete that step first so there's a significant gene "
            "list to run enrichment analysis on."
        )
