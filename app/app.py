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
  - setup_workspace.py         -> Setup & Deployment (environment check, HPC SSH connections)
  - single_cell/singlecell_workspace.py -> Single-cell RNA-Seq (Phase 1: ingestion
    through STARsolo alignment/cell-calling, plus a Phase 2 Cell-level QC
    page; see that module's own docstring for scope -- droplet-based UMI
    methods only)

Keeping workspaces in separate files means work on one pipeline (e.g. Bulk
RNA-Seq) can never accidentally break another (e.g. Spatial Transcriptomics).

### Sidebar structure

Rather than one long flat list of radio buttons, workspace steps are grouped
into collapsible sidebar drawers ("expanders") -- one per pipeline -- via
render_pipeline_section() below. "Portal Home" and "Setup & Deployment" both
sit outside any pipeline group as standalone buttons, since neither is a
pipeline step -- Setup & Deployment is a cross-cutting utility page (see
setup_workspace.py's own module docstring), exactly like Portal Home.

Because a workspace's step now lives inside a _per-pipeline_ radio widget
instead of one single radio, routing can no longer be driven directly by a
single widget's session_state value. Instead, a plain (non-widget) key,
st.session_state["active_workspace"], is the single source of truth used
for routing below. It is updated in one of two ways:
  1. A grouped radio's on_change callback, fired the moment the user
     picks a step inside an expanded pipeline drawer.
  2. The pending-nav-request block below, used when another workspace
     module (e.g. the "Proceed to Trimming" button at the end of
     bulk_rnaseq_workspace.py, singlecell_workspace.py's own "Proceed to
     Phase 2: Cell-level QC" button, or project_actions.py's "Begin DE
     Analysis" hand-off from Auto/Monitor Mode) requests navigation by
     setting st.session_state["nav_request"] and calling st.rerun().

Note on cross-page navigation and expander state: Streamlit's st.expander
has no key parameter (even as of Streamlit 1.52), so it cannot track its
own open/closed state in session_state, and manually toggling one via the
UI is not visible to the script until some other rerun occurs -- at which
point whatever expanded= value the code passes will win. To make this
work _with_ that constraint rather than fight it, each pipeline's drawer is
simply expanded whenever it contains the currently active_workspace, and
collapsed otherwise. In practice this means the drawer for whichever
pipeline you're currently working in stays open (revealing its steps),
while other pipelines stay tucked away until you navigate into them.

--- single_cell/ subfolder import note (2026-08-17) ---
singlecell_workspace.py and its supporting modules (sc_project_manager.py,
chemistry_manager.py, singlecell_ingestion_manager.py,
singlecell_trim_manager.py, starsolo_manager.py, sc_cellqc_manager.py)
live in a single_cell/ subfolder rather than directly alongside this
file, for cleaner file organization as this pipeline grows. Since
Streamlit is always launched via `cd repo/app && streamlit run app.py`
(see DEPLOYMENT.md), app.py's own working directory is repo/app/ --
sys.path.insert below adds single_cell/ onto the import path so `import
singlecell_workspace` resolves normally, the same as any other
top-level module in this app, without needing package-relative imports
or an __init__.py-based package structure.

--- Single-cell Phase 2 placeholder route added (2026-08-17) ---
Added "🔬 SC Cell-level QC" as a fourth option in the "single_cell"
pipeline group below, plus a matching elif branch calling
singlecell_workspace.render_cell_qc(). This was added specifically so
Step 6's (STARsolo alignment) "➡️ Proceed to Phase 2: Cell-level QC"
button has a real, working navigation target -- without a matching
PIPELINE_GROUPS entry AND elif branch, that button's nav_request would
set active_workspace to a value matching no branch below, silently
rendering a blank page. render_cell_qc() itself now runs real Phase 2
logic (scDblFinder doublet detection, DecontX/SoupX ambient RNA
correction, adaptive per-cell filtering) via sc_cellqc_manager.py.
"""
import os
import sys
from functools import partial

import streamlit as st

# Make single_cell/'s modules importable as plain top-level modules (see
# this file's own docstring note above) -- inserted once, at import time,
# before singlecell_workspace itself (or anything it depends on) is
# imported below.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "single_cell"))

import spatial_workspace
import bulk_rnaseq_workspace
import advanced_mode_workspace
import monitor_mode_workspace
import trimming_workspace
import alignment_workspace
import differential_expression_workspace
import ontology_workspace
import setup_workspace
import singlecell_workspace

# 1. Global Setup Layout

st.set_page_config(layout="wide", page_title="Multi-Omics Bioinformatics Portal")
HOME_OPTION = "📊 Portal Home"
SETUP_OPTION = "⚙️ Setup & Deployment"

# Pipeline groups: each becomes one collapsible sidebar drawer titled title,
# containing a radio button list of that pipeline's sequential workspace
# steps (options, in pipeline order). Adding a future pipeline is just one
# more entry here -- no other code below needs to change.

# --- Drawer order (2026-08-17) ---
# Reordered to: Portal Home (standalone button, unaffected by this dict's
# order) -> Bulk RNA-Seq -> Single-cell RNA-Seq -> Spatial Transcriptomics
# -> Advanced Modes -> Setup & Deployment (standalone button). Python
# dicts preserve insertion order, and render_pipeline_section() below
# iterates PIPELINE_GROUPS.items() in that same order, so this dict's
# literal key order below IS the sidebar's actual top-to-bottom order.
PIPELINE_GROUPS = {
    "bulk_rnaseq": {
        "title": "🧬 Bulk RNA-Seq",
        "options": [
            "🧬 Bulk RNA-Seq Pipeline",
            "🧪 Trimming & Post-Trim QC",
            "🧮 RNA Alignment & Counts",
            "🌋 Differential Expression",
            "🧬 Ontology Analysis",
        ],
    },
    # Single-cell RNA-Seq: its own top-level drawer, separate from
    # "bulk_rnaseq" above, since single-cell projects live in their own
    # namespace (sc_project_manager.SC_PROJECTS_ROOT) and follow a
    # meaningfully different step sequence (chemistry detection, R1/R2-
    # aware trimming, STARsolo cell-calling).
    #
    # --- Multi-option fix (2026-08-17) ---
    # This drawer previously had only ONE option ("🧫 Single-cell
    # RNA-Seq" covering the whole pipeline as a single page), which was a
    # real navigation bug, not just a stylistic difference from Bulk
    # RNA-Seq: Streamlit only fires a radio widget's on_change callback
    # when its value actually CHANGES. With only one possible option,
    # that radio's value could never change, so on_change could never
    # fire, so clicking it could never navigate there at all. Splitting
    # into genuinely distinct step options -- mirroring Bulk RNA-Seq's
    # own Pipeline/Ingestion -> Trimming & Post-Trim QC -> Alignment &
    # Counts grouping exactly -- fixes that bug as a direct consequence
    # of giving this drawer the same step-by-step radio flow Bulk
    # RNA-Seq already has, not as a separate patch on top of it.
    #
    # --- Phase 2 route added (2026-08-17) ---
    # "🔬 SC Cell-level QC" appended as a fourth option -- gives Step 6's
    # "Proceed to Phase 2" button (singlecell_workspace.py) a real
    # navigation target; see this file's own module docstring section
    # "Single-cell Phase 2 placeholder route added" above for why both
    # this entry AND the matching elif branch below are required
    # together. Later phases (clustering/annotation, pseudobulk DE) are
    # planned as additional options appended to this same list.
    "single_cell": {
        "title": "🧫 Single-cell RNA-Seq",
        "options": [
            "🧫 Single-cell RNA-Seq",
            "🧪 SC Trimming & Post-Trim QC",
            "🧬 SC Alignment & Cell-Calling",
            "🔬 SC Cell-level QC",
        ],
    },
    "spatial": {
        "title": "🧠 Spatial Transcriptomics",
        "options": [
            "🧠 10x Genomics Visium Viewer",
        ],
    },
    # "Advanced Modes" is its own top-level drawer, separate from any
    # single pipeline's own step-by-step drawer (e.g. "bulk_rnaseq"
    # above) -- since the "Auto" workspace lets a user pick WHICH
    # pipeline to run non-interactively via a dropdown on its own first
    # screen (see advanced_mode_workspace.py's render()), rather than
    # being one more step nested inside any one pipeline's own drawer.
    # Monitor Mode (background folder-watcher automation) lives in this
    # same drawer alongside Auto.
    "advanced_modes": {
        "title": "🚀 Advanced Modes",
        "options": [
            "🤖 Auto",
            "📁 Monitor Mode",
        ],
    },
}


def _activate(radio_key: str):
    """
    on_change callback for a grouped radio widget: copies that widget's
    current value into the single active_workspace routing key used by
    the elif-chain below, so picking a step inside ANY pipeline's drawer
    updates routing the same way regardless of which drawer it came from.
    """
    st.session_state["active_workspace"] = st.session_state[radio_key]


def render_pipeline_section(name: str, options: list[str], key: str, icon: str = ""):
    """
    Render one collapsible sidebar drawer (expander) titled name,
    containing a radio button list options for that pipeline's steps.
    The drawer is expanded whenever it contains the currently active
    workspace (so the pipeline you're actively working in stays open),
    and collapsed otherwise -- see this module's docstring for why this
    "derive expanded= from current state" approach is used instead of
    trying to track the expander's own open/closed state directly.
    """
    radio_key = f"_workspace_radio_{key}"
    active = st.session_state.get("active_workspace")
    is_active_group = active in options
    with st.sidebar.expander(name, expanded=is_active_group):
        st.radio(
            "Go to:", options,
            index=options.index(active) if active in options else 0,
            key=radio_key, label_visibility="collapsed",
            on_change=_activate, args=(radio_key,),
        )


# ------------------------------------------------------------------
# Apply any pending navigation request from another workspace module
# (e.g. a "Proceed to Trimming" button) BEFORE any widgets below are
# instantiated. Streamlit forbids setting a widget-bound session_state
# key after that widget already exists for the current run, so this
# indirection (nav_request -> the target step's own group radio key)
# must happen first, exactly as it did with the single flat
# assay_choice_radio previously.
# ------------------------------------------------------------------

if "nav_request" in st.session_state:
    target = st.session_state.pop("nav_request")
    st.session_state["active_workspace"] = target
    for group_key, group in PIPELINE_GROUPS.items():
        if target in group["options"]:
            st.session_state[f"_workspace_radio_{group_key}"] = target
            break

if "active_workspace" not in st.session_state:
    st.session_state["active_workspace"] = HOME_OPTION

# 2. Main Portal Routing Menu

st.sidebar.title("🧬 Multi-Omics Portal")
st.sidebar.markdown("---")

if st.sidebar.button(HOME_OPTION, key="home_button", use_container_width=True):
    st.session_state["active_workspace"] = HOME_OPTION

for group_key, group in PIPELINE_GROUPS.items():
    render_pipeline_section(group["title"], group["options"], key=group_key)

st.sidebar.markdown("---")

# Standalone button, same pattern as HOME_OPTION above -- Setup &
# Deployment is a cross-cutting utility page (environment/dependency
# checker, HPC SSH connection management), not a pipeline step, so it
# does not belong inside any PIPELINE_GROUPS drawer.
if st.sidebar.button(SETUP_OPTION, key="setup_button", use_container_width=True):
    st.session_state["active_workspace"] = SETUP_OPTION

st.sidebar.markdown("---")

assay_choice = st.session_state["active_workspace"]

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
# 🤖 ADVANCED MODES: AUTO
# ==========================================

elif assay_choice == "🤖 Auto":
    advanced_mode_workspace.render()

# ==========================================
# 📁 ADVANCED MODES: MONITOR MODE
# ==========================================

elif assay_choice == "📁 Monitor Mode":
    monitor_mode_workspace.render()

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

# =========================================
# ⚙️ WORKSPACE 8: SETUP & DEPLOYMENT
# =========================================

elif assay_choice == SETUP_OPTION:
    setup_workspace.render()

# =========================================
# 🧫 WORKSPACE 9: SINGLE-CELL RNA-SEQ
# =========================================

elif assay_choice == "🧫 Single-cell RNA-Seq":
    singlecell_workspace.render_ingestion()

elif assay_choice == "🧪 SC Trimming & Post-Trim QC":
    singlecell_workspace.render_trimming()

elif assay_choice == "🧬 SC Alignment & Cell-Calling":
    singlecell_workspace.render_alignment()

# Phase 2 route (2026-08-17) -- see this file's own module docstring
# section "Single-cell Phase 2 placeholder route added" and
# singlecell_workspace.render_cell_qc()'s own docstring for context.
elif assay_choice == "🔬 SC Cell-level QC":
    singlecell_workspace.render_cell_qc()
