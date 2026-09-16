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
  - single_cell/sc_downstream_workspace.py -> Single-cell RNA-Seq Phase 3
    (multi-sample downstream analysis: normalization through clustering,
    embeddings, cell-type annotation, pseudobulk -> DESeq2 bridge, and
    compositional analysis; see that module's own docstring for scope)
  - auth_manager.py -> Authentication (login/bootstrap) and admin-managed
    custom roles/permissions -- see that module's own docstring for the
    full design (built-in admin/tech roles, admin-definable custom roles,
    per-action permission catalog).

Keeping workspaces in separate files means work on one pipeline (e.g. Bulk
RNA-Seq) can never accidentally break another (e.g. Spatial Transcriptomics).

### Authentication gate (2026-08-24) ---

Every run of this script now starts by calling auth.render_login_gate()
BEFORE any sidebar/routing/workspace code executes at all -- if that
call returns False (no one is logged in yet: either a brand-new
deployment showing the one-time "create first admin" bootstrap form, or
an existing deployment showing an ordinary login form), this script
returns immediately without rendering the sidebar, any workspace, or
leaking ANY information about what pipelines/projects exist to an
unauthenticated visitor. Only once render_login_gate() returns True
(a real, verified session) does the rest of this file's existing
routing logic run, completely unchanged from before this gate was
added -- auth_manager.py's own role/permission checks are enforced
individually, deeper inside each gated action (e.g.
project_manager.py's delete flow, setup_workspace.py's eggNOG/install/
HPC sections), not by this top-level gate, which only answers "is
ANYONE allowed to see the app at all right now."

auth.render_user_badge() is called once, in the sidebar, immediately
after the gate passes -- shows the current username/role and a logout
button on every page for the rest of the session.

A new standalone sidebar entry, "👥 User & Role Management" (mirroring
"⚙️ Setup & Deployment"'s own existing pattern -- a cross-cutting
utility page outside any pipeline drawer), routes to
auth.render_user_management(). That page independently re-enforces its
own admin-only check the moment it renders (see auth_manager.py's own
render_user_management() docstring) -- this router file does not
attempt its own separate admin check before routing to it, matching
this app's own established "let the destination page re-verify its own
authorization, don't rely solely on a hidden sidebar entry" convention
already used for the Setup & Deployment page's own internal admin-gated
sections.

--- single_cell/ subfolder import note (2026-08-17) ---
singlecell_workspace.py and its supporting modules (sc_project_manager.py,
chemistry_manager.py, singlecell_ingestion_manager.py,
singlecell_trim_manager.py, starsolo_manager.py, sc_cellqc_manager.py,
sc_downstream_manager.py, sc_downstream_workspace.py) live in a
single_cell/ subfolder rather than directly alongside this file, for
cleaner file organization as this pipeline grows. Since Streamlit is
always launched via `cd repo/app && streamlit run app.py` (see
DEPLOYMENT.md), app.py's own working directory is repo/app/ --
sys.path.insert below adds single_cell/ onto the import path so `import
singlecell_workspace` (and `import sc_downstream_workspace`) resolve
normally, the same as any other top-level module in this app, without
needing package-relative imports or an __init__.py-based package
structure.

--- Phase 3 sidebar entry renamed (2026-08-25, shortened again same day) ---
The Single-cell Phase 3 downstream-analysis sidebar entry was renamed
from the original generic "📊 SC Downstream Analysis" to a more
informative name -- per direct user request, to make the sidebar itself
more informative about what this step actually does, rather than
requiring a user to already know that "downstream analysis" means
clustering and cell-type annotation. The FIRST more-informative name
tried ("📊 SC Analysis: Normalization → Clustering → Cell-Type
Annotation") was itself then shortened, again per direct user feedback
that it was too long, to the current "📊 SC Analysis: Clustering & Cell
Annotation". This same exact string is also used by
singlecell_workspace.py's own "Proceed to SC Analysis..." button -- both
files reference the identical string value via a shared named constant
(singlecell_workspace.SC_DOWNSTREAM_ANALYSIS_OPTION, imported and reused
directly below rather than re-typed as a second, independent string
literal here) specifically so these two files can never silently drift
out of sync with each other if this name is ever changed again in the
future. This is a pure rename -- the underlying route still points at
the exact same sc_downstream_workspace.render() call as before; no
change to that workspace's own scope or behavior.
"""
import os
import sys
from functools import partial

import streamlit as st

import auth_manager as auth


# 1. Global Setup Layout -- set_page_config must run before ANY other
# Streamlit call in the script, including the login gate below, per
# Streamlit's own requirement that it be the first st.* call made.
st.set_page_config(layout="wide", page_title="Pretty Awesome Transcriptomics Tool")

# ---------------------------------------------------------------------
# Authentication gate -- must run before any sidebar/routing/workspace
# code below. See this file's own module docstring, "Authentication
# gate", for the full rationale.
# ---------------------------------------------------------------------
if not auth.render_login_gate():
    st.stop()

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
import sc_downstream_workspace
import sc_deseq2_workspace
import sc_ontology_workspace
import sc_comparecluster_workspace

HOME_OPTION = "📊 Portal Home"
SETUP_OPTION = "⚙️ Setup & Deployment"
USER_MANAGEMENT_OPTION = "👥 User & Role Management"

# --- Phase 3 sidebar entry renamed (2026-08-25) -- see this file's own
# module docstring, "Phase 3 sidebar entry renamed", for the full
# rationale. Reused DIRECTLY from singlecell_workspace.py's own shared
# named constant (rather than re-typed as a second, independent string
# literal here) so this router's own PIPELINE_GROUPS/routing entry and
# that module's own "Proceed to SC Analysis..." button can never
# silently drift out of sync with each other if this name is ever
# changed again in the future.
SC_DOWNSTREAM_ANALYSIS_OPTION = singlecell_workspace.SC_DOWNSTREAM_ANALYSIS_OPTION

# Pipeline groups: each becomes one collapsible sidebar drawer titled title,
# containing a radio button list of that pipeline's sequential workspace
# steps (options, in pipeline order). Adding a future pipeline is just one
# more entry here -- no other code below needs to change.

# --- Drawer order (2026-08-17) ---
# Reordered to: Portal Home (standalone button, unaffected by this dict's
# order) -> Bulk RNA-Seq -> Single-cell RNA-Seq -> Spatial Transcriptomics
# -> Advanced Modes -> Setup & Deployment (standalone button) -> User &
# Role Management (standalone button). Python dicts preserve insertion
# order, and render_pipeline_section() below iterates
# PIPELINE_GROUPS.items() in that same order, so this dict's literal key
# order below IS the sidebar's actual top-to-bottom order.
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
    # --- Phase 3 route (2026-08-23, renamed 2026-08-25) ---
    # SC_DOWNSTREAM_ANALYSIS_OPTION routes to
    # sc_downstream_workspace.render(), covering ALL of Phase 3 (Steps
    # 1-10) as a single page -- see that module's own docstring for why
    # Phase 3 is not split across multiple pages the way Phase 1/2 are.
    # See this file's own module docstring, "Phase 3 sidebar entry
    # renamed", for why this is now a shared named constant rather than
    # a bare string literal.
    "single_cell": {
        "title": "🧫 Single-cell RNA-Seq",
        "options": [
            "🧫 Single-cell RNA-Seq", "🧪 SC Trimming & Post-Trim QC",
            "🧬 SC Alignment & Cell-Calling", "🔬 SC Cell-level QC",
            SC_DOWNSTREAM_ANALYSIS_OPTION,
            "🌋 SC DESeq2", 
            "🧬 SC Ontology Analysis","🔀 SC compareCluster",
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

st.sidebar.title("🧬 Pretty Awesome Transcriptomics Tool")
auth.render_user_badge()
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

# Standalone button, same pattern as SETUP_OPTION above -- User & Role
# Management is likewise a cross-cutting utility page, not a pipeline
# step. Shown to EVERY logged-in session regardless of role (not just
# admins) -- auth.render_user_management() itself immediately shows an
# access-denied message and renders nothing further for a non-admin
# session (see that function's own docstring) rather than this router
# file trying to hide the button for non-admins ahead of time. This
# matches setup_workspace.py's own established convention of gating
# each admin-only SECTION independently, rather than hiding an entire
# page's sidebar entry based on role.
if st.sidebar.button(USER_MANAGEMENT_OPTION, key="user_management_button", use_container_width=True):
    st.session_state["active_workspace"] = USER_MANAGEMENT_OPTION

st.sidebar.markdown("---")

assay_choice = st.session_state["active_workspace"]

# ==========================================
# 🏠 WORKSPACE 1: PORTAL HOME
# ==========================================

if assay_choice == HOME_OPTION:
    st.title("🧬 Pretty Awesome Transcriptomics Tool")
    st.caption("Yes, that spells P.A.T. -- a little nod to the person who built it.")
    st.markdown(
        "Welcome! This portal walks you through complete RNA-seq analysis "
        "workflows -- from raw reads to differential expression and "
        "biological interpretation -- without requiring a bioinformatics "
        "background. Pick a workflow below to get started, or continue an "
        "existing project from the sidebar."
    )
    st.markdown("---")

    col1, col2, col3 = st.columns(3)

    with col1:
        with st.container(border=True):
            st.markdown("### 🧬 Bulk RNA-Seq")
            st.markdown(
                "Analyze gene expression averaged across a whole tissue or "
                "sample -- the classic RNA-seq workflow. Upload FASTQ files, "
                "trim and QC your reads, align and quantify against a "
                "reference genome, run DESeq2 to find differentially "
                "expressed genes between conditions, and finish with "
                "GO/KEGG/Reactome enrichment to understand the biology "
                "behind your results."
            )
            st.caption(
                "Best for: comparing conditions (e.g. treated vs. control) "
                "at the tissue/sample level."
            )
            if st.button(
                "Start Bulk RNA-Seq →", key="home_card_bulk_btn",
                use_container_width=True, type="primary",
            ):
                st.session_state["nav_request"] = "🧬 Bulk RNA-Seq Pipeline"
                st.rerun()

    with col2:
        with st.container(border=True):
            st.markdown("### 🧫 Single-cell RNA-Seq")
            st.markdown(
                "Analyze gene expression at INDIVIDUAL CELL resolution "
                "instead of a tissue-wide average. Ingest droplet-based "
                "(10x-style) FASTQ files, detect chemistry automatically, "
                "align with STARsolo, run cell-level QC (doublets, ambient "
                "RNA), then cluster, annotate cell types, and explore "
                "composition differences across your samples."
            )
            st.caption(
                "Best for: discovering which cell types exist and how they "
                "individually respond to a condition."
            )
            if st.button(
                "Start Single-cell RNA-Seq →", key="home_card_sc_btn",
                use_container_width=True, type="primary",
            ):
                st.session_state["nav_request"] = "🧫 Single-cell RNA-Seq"
                st.rerun()

    with col3:
        with st.container(border=True):
            st.markdown("### ⚙️ Setup & Deployment")
            st.markdown(
                "Check which external tools (STAR, Salmon, R/DESeq2, "
                "clusterProfiler, etc.) are available in this environment, "
                "install anything missing, and manage HPC SSH connections "
                "for running heavier jobs on a remote cluster."
            )
            st.caption(
                "Best for: first-time setup, or troubleshooting a 'tool not "
                "found' error on another page."
            )
            if st.button(
                "Go to Setup & Deployment →", key="home_card_setup_btn",
                use_container_width=True,
            ):
                st.session_state["nav_request"] = SETUP_OPTION
                st.rerun()

    st.markdown("---")
    st.caption(
        "🧠 Also available: Spatial Transcriptomics (10x Visium) and "
        "Advanced Modes (Auto / Monitor Mode) -- see the sidebar drawers on "
        "the left."
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
# 👥 WORKSPACE: USER & ROLE MANAGEMENT
# =========================================

elif assay_choice == USER_MANAGEMENT_OPTION:
    auth.render_user_management()

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

elif assay_choice == "🌋 SC DESeq2":
    sc_deseq2_workspace.render()
elif assay_choice == "🧬 SC Ontology Analysis":
  sc_ontology_workspace.render()
elif assay_choice == "🔀 SC compareCluster":
  sc_comparecluster_workspace.render()  
    # Phase 3 route (2026-08-23, renamed 2026-08-25) -- see this file's own
# module docstring section "Phase 3 sidebar entry renamed" and
# sc_downstream_workspace.render()'s own docstring for context. Routes
# to the SAME single page for all of Steps 1-10 (see that module's own
# docstring for why Phase 3 is not split across multiple pages the way
# Phase 1/2 are) -- this is a pure rename of the sidebar entry string;
# the underlying route/behavior is completely unchanged.
elif assay_choice == SC_DOWNSTREAM_ANALYSIS_OPTION:
    sc_downstream_workspace.render()
