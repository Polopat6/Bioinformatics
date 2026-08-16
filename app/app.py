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

Sidebar structure
------------------
Rather than one long flat list of radio buttons, workspace steps are grouped
into collapsible sidebar drawers ("expanders") -- one per pipeline -- via
render_pipeline_section() below. "Portal Home" sits outside any pipeline
group as a standalone button, since it isn't a pipeline step.

Because a workspace's step now lives inside a *per-pipeline* radio widget
instead of one single radio, routing can no longer be driven directly by a
single widget's session_state value. Instead, a plain (non-widget) key,
st.session_state["active_workspace"], is the single source of truth used
for routing below. It is updated in one of two ways:
    1. A grouped radio's on_change callback, fired the moment the user
       picks a step inside an expanded pipeline drawer.
    2. The pending-nav-request block below, used when another workspace
       module (e.g. the "Proceed to Trimming" button at the end of
       bulk_rnaseq_workspace.py) requests navigation by setting
       st.session_state["nav_request"] and calling st.rerun().

Note on cross-page navigation and expander state: Streamlit's st.expander
has no `key` parameter (even as of Streamlit 1.52), so it cannot track its
own open/closed state in session_state, and manually toggling one via the
UI is not visible to the script until some other rerun occurs -- at which
point whatever `expanded=` value the code passes will win. To make this
work *with* that constraint rather than fight it, each pipeline's drawer is
simply expanded whenever it contains the currently active_workspace, and
collapsed otherwise. In practice this means the drawer for whichever
pipeline you're currently working in stays open (revealing its steps),
while other pipelines stay tucked away until you navigate into them.
"""

from functools import partial

import streamlit as st

import spatial_workspace
import bulk_rnaseq_workspace
import advanced_mode_workspace
import monitor_mode_workspace
import trimming_workspace
import alignment_workspace
import differential_expression_workspace
import ontology_workspace

# 1. Global Setup Layout
st.set_page_config(layout="wide", page_title="Multi-Omics Bioinformatics Portal")

HOME_OPTION = "📊 Portal Home"

# Pipeline groups: each becomes one collapsible sidebar drawer titled `title`,
# containing a radio button list of that pipeline's sequential workspace
# steps (`options`, in pipeline order). Adding a future pipeline (e.g.
# Single-cell RNA-seq) is just one more entry here -- no other code below
# needs to change.
PIPELINE_GROUPS = {
    "spatial": {
        "title": "🧠 Spatial Transcriptomics",
        "options": [
            "🧠 10x Genomics Visium Viewer",
        ],
    },
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
    # "Advanced Modes" is its own top-level drawer, separate from any
    # single pipeline's own step-by-step drawer (e.g. "bulk_rnaseq"
    # above) -- since the "Auto" workspace lets a user pick WHICH
    # pipeline to run non-interactively via a dropdown on its own first
    # screen (see advanced_mode_workspace.py's render()), rather than
    # being one more step nested inside any one pipeline's own drawer.
    # Monitor Mode (background folder-watcher automation) is intended
    # to live in this same drawer once its workspace UI is built -- see
    # advanced_mode/README.md's "What's still needed" section -- so
    # this group is deliberately structured to hold more than one
    # option going forward, even though "Auto" is the only one wired
    # up today.
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
    on_change callback for a grouped radio widget.

    Promotes the radio's newly clicked value to `active_workspace`, the
    single non-widget session_state key that drives routing below. Fired
    by Streamlit before the script reruns, so by the time the routing
    block executes, active_workspace already reflects this click.
    """
    st.session_state["active_workspace"] = st.session_state[radio_key]


def render_pipeline_section(name: str, options: list[str], key: str, icon: str = ""):
    """
    Render one collapsible sidebar drawer (expander) titled `name`,
    containing a radio button list `options` for that pipeline's steps.

    Parameters
    ----------
    name : str
        Display name of the pipeline, e.g. "Bulk RNA-Seq".
    options : list[str]
        Ordered list of workspace step labels for this pipeline, exactly
        matching the strings used in the routing if/elif block below.
    key : str
        Unique key prefix for this pipeline's radio widget
        (e.g. "bulk_rnaseq"). Must be unique across all pipeline groups.
    icon : str, optional
        Optional emoji/icon prefix for the drawer title.

    Notes
    -----
    The drawer auto-expands whenever `active_workspace` is one of its own
    `options` (i.e. whenever the user is currently working in this
    pipeline), and stays collapsed otherwise. See the module docstring for
    why this -- rather than trying to persist a manual user toggle -- is
    the appropriate approach given st.expander's lack of a `key` argument.
    """
    radio_key = f"{key}_workspace_radio"
    is_active_group = st.session_state.get("active_workspace") in options

    with st.sidebar.expander(f"{icon} {name}".strip(), expanded=is_active_group):
        st.radio(
            "Select step:",
            options,
            # index=None means no option is pre-checked until the user (or a
            # nav_request, via session_state) picks one. This matters even
            # for groups with a single step (e.g. Spatial currently has
            # only one): with a pre-selected default, Streamlit's on_change
            # only fires when the value actually *changes*, so clicking an
            # already-selected (and, for single-step groups, *always*
            # already-selected) option would silently do nothing. Starting
            # unselected guarantees every click is a real change and fires
            # _activate.
            index=None,
            key=radio_key,
            on_change=partial(_activate, radio_key),
            label_visibility="collapsed",
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
            st.session_state[f"{group_key}_workspace_radio"] = target
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
