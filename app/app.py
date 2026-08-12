"""
app.py

Thin router / entry point for the Multi-Omics Bioinformatics Portal.

This file should rarely need to change. All workspace-specific logic lives
in its own module:
    - spatial_workspace.py       -> 10x Genomics Visium Viewer
    - bulk_rnaseq_workspace.py   -> Bulk RNA-Seq Pipeline (upload, matching, QC)
    - trimming_workspace.py      -> Adapter Trimming & Post-Trim QC
    - alignment_workspace.py     -> RNA Alignment & Counts (Salmon quantification)
    - differential_expression_workspace.py -> Differential Expression (DESeq2)
    - ontology_workspace.py      -> Ontology Analysis (GO/pathway enrichment)

Keeping workspaces in separate files means work on one pipeline (e.g. Bulk
RNA-Seq) can never accidentally break another (e.g. Spatial Transcriptomics).

Note on cross-page navigation: Streamlit does not allow a widget's
session_state value (here, key="assay_choice_radio") to be reassigned
after that widget has already been instantiated in the current run. So
other workspace modules (e.g. the "Proceed to Trimming" button at the end
of bulk_rnaseq_workspace.py) do NOT set st.session_state["assay_choice_radio"]
directly. Instead, they set a separate, plain (non-widget) key,
st.session_state["nav_request"], to the target page name and call
st.rerun(). On the next run, this file checks "nav_request" and applies it
to "assay_choice_radio" BEFORE the radio widget is created below, which is
allowed. See the "Apply any pending navigation request" block below.
"""

import streamlit as st

import spatial_workspace
import bulk_rnaseq_workspace
import trimming_workspace
import alignment_workspace
import differential_expression_workspace
import ontology_workspace

# 1. Global Setup Layout
st.set_page_config(layout="wide", page_title="Multi-Omics Bioinformatics Portal")

# Apply any pending navigation request from another workspace module
# (e.g. a "Proceed to X" button) before the radio widget below is
# instantiated. This indirection is required because Streamlit forbids
# setting a widget-bound session_state key after that widget already
# exists for the current run.
if "nav_request" in st.session_state:
    st.session_state["assay_choice_radio"] = st.session_state.pop("nav_request")

# 2. Main Portal Routing Menu
st.sidebar.title("🧬 Multi-Omics Portal")
st.sidebar.markdown("---")
assay_choice = st.sidebar.radio(
    "Select Analysis Workspace:",
    [
        "📊 Portal Home",
        "🧠 10x Genomics Visium Viewer",
        "🧬 Bulk RNA-Seq Pipeline",
        "🧪 Trimming & Post-Trim QC",
        "🧮 RNA Alignment & Counts",
        "🌋 Differential Expression",
        "🧬 Ontology Analysis",
    ],
    key="assay_choice_radio",
)

# ==========================================
# 🏠 WORKSPACE 1: PORTAL HOME
# ==========================================
if assay_choice == "📊 Portal Home":
    st.title("🎛️ Clinical Multi-Omics Orchestration Hub")
    st.markdown(
        "Welcome to the production-grade multi-omics interface. Select an "
        "active pipeline workspace from the sidebar menu to process "
        "structural data metrics or view interactive expression atlases."
    )

    col1, col2 = st.columns(2)
    with col1:
        st.info("### 🧠 Spatial Transcriptomics\nProcess micro-array slides, format H&E histology matrices, and map localized cluster counts using DuckDB + Plotly overlays.")
    with col2:
        st.info("### 🧬 Bulk RNA-Seq\nUpload raw sequencing runs, evaluate sample quality structures, and execute differential gene expression layouts.")

# ==========================================
# 🧠 WORKSPACE 2: SPATIAL TRANSCRIPTOMICS
# ==========================================
elif assay_choice == "🧠 10x Genomics Visium Viewer":
    spatial_workspace.render()

# ==========================================
# 🧬 WORKSPACE 3: BULK RNA-SEQ WORKFLOW
# ==========================================
elif assay_choice == "🧬 Bulk RNA-Seq Pipeline":
    bulk_rnaseq_workspace.render()

# ==========================================
# 🧪 WORKSPACE 4: TRIMMING & POST-TRIM QC
# ==========================================
elif assay_choice == "🧪 Trimming & Post-Trim QC":
    trimming_workspace.render()

# ==========================================
# 🧮 WORKSPACE 5: RNA ALIGNMENT & COUNTS
# ==========================================
elif assay_choice == "🧮 RNA Alignment & Counts":
    alignment_workspace.render()

# =========================================
# 🌋 WORKSPACE 6: DIFFERENTIAL EXPRESSION
# =========================================
elif assay_choice == "🌋 Differential Expression":
    differential_expression_workspace.render()

# =========================================
# 🧬 WORKSPACE 7: ONTOLOGY ANALYSIS
# =========================================
elif assay_choice == "🧬 Ontology Analysis":
    ontology_workspace.render()

