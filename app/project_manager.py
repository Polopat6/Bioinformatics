"""
project_manager.py

Lightweight "project" system shared across Bulk RNA-Seq workflow steps.

A project is just a named folder under PROJECTS_ROOT that keeps every
uploaded file and result together, so a user can:
  - stop after any step and come back to it later without re-uploading
  - move from ingestion -> QC -> trimming -> post-trim QC -> alignment ->
    counts using the same fixed set of files and sample matches

This module only manages the project's *identity* and folder structure.
It doesn't know anything about FASTQ files, QC, trimming, or alignment —
that logic stays in each workspace module. Those modules ask this one for
"where do my files live for project X" and "what's already been done".
"""

import json
import os
from datetime import datetime

import streamlit as st

PROJECTS_ROOT = "data/projects"


# ---------------------------------------------------------------------------
# Path helpers — every workspace step should get its file locations from
# these functions rather than hardcoding paths, so everything stays
# scoped correctly to the active project.
# ---------------------------------------------------------------------------

def list_projects():
    """Return a sorted list of existing project names."""
    if not os.path.isdir(PROJECTS_ROOT):
        return []
    return sorted([
        name for name in os.listdir(PROJECTS_ROOT)
        if os.path.isdir(os.path.join(PROJECTS_ROOT, name))
    ])


def project_dir(project_name):
    return os.path.join(PROJECTS_ROOT, project_name)


def fastq_dir(project_name):
    return os.path.join(project_dir(project_name), "fastq")


def samplesheet_path(project_name):
    return os.path.join(project_dir(project_name), "samplesheet.csv")


def metadata_path(project_name):
    """
    Path where the user's uploaded metadata is persisted, normalized to
    CSV regardless of whether the original upload was .csv or .xlsx.
    Storing a normalized copy means we don't need to remember the
    original file format to reload it later.
    """
    return os.path.join(project_dir(project_name), "metadata.csv")


def fastqc_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "fastqc")


def multiqc_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "multiqc")


def trimmed_fastq_dir(project_name):
    """Trimmed FASTQ output from fastp."""
    return os.path.join(project_dir(project_name), "trimmed")


def fastp_reports_dir(project_name):
    """
    Per-sample fastp QC reports (fastp.json + fastp.html), generated
    automatically by fastp as a side effect of trimming. The JSON files
    here are what MultiQC reads to build the combined post-trim report.
    """
    return os.path.join(project_dir(project_name), "qc", "fastp")


def posttrim_multiqc_dir(project_name):
    """Combined MultiQC report built from the fastp JSON reports."""
    return os.path.join(project_dir(project_name), "qc", "multiqc_posttrim")


def reference_dir(project_name):
    """Where the downloaded/uploaded transcriptome reference + Salmon index live."""
    return os.path.join(project_dir(project_name), "reference")


def salmon_index_dir(project_name):
    """The built Salmon index, generated once per reference and reused across runs."""
    return os.path.join(reference_dir(project_name), "salmon_index")


def salmon_quant_dir(project_name):
    """Per-sample Salmon quantification output (one subfolder per sample)."""
    return os.path.join(project_dir(project_name), "quant", "salmon")


def star_index_dir(project_name):
    """The built STAR genome index, generated once per reference and reused across runs."""
    return os.path.join(reference_dir(project_name), "star_index")


def star_align_dir(project_name):
    """Per-sample STAR alignment output (BAM files + ReadsPerGene.out.tab)."""
    return os.path.join(project_dir(project_name), "quant", "star")


def counts_matrix_path(project_name):
    """The final combined gene-level counts matrix (all samples, one row per gene)."""
    return os.path.join(project_dir(project_name), "quant", "gene_counts_matrix.csv")

# --- Additions needed in project_manager.py for the Differential Expression workspace ---
# Add these functions alongside the existing path helpers (e.g. near counts_matrix_path).
# --- Addition for project_manager.py ---
# Add this function alongside the other shared utility functions (not
# project-specific, so it doesn't need a project_name argument).



def get_recommended_thread_count(reserve_cores=1, max_default=8):
    """
    Detect how many CPU cores are available on this machine and return
    a sensible recommended default thread/worker count for CPU-bound
    tools (FastQC, fastp, Salmon, STAR, etc.), used to pre-fill each
    tool's thread-count slider rather than hardcoding a fixed default
    (e.g. 4) that may be far too high for a small machine or
    unnecessarily conservative for a large one.

    reserve_cores: how many cores to leave free for the OS, the
        Streamlit process itself, and general system responsiveness —
        default 1, so a tool never recommends using literally every
        available core.
    max_default: an upper cap on the *recommended* value, independent of
        how many cores are actually detected. This exists because tools
        like fastp are documented to see no real benefit (and some
        internally hard-cap) beyond ~16 threads due to I/O bottlenecks
        rather than CPU availability — recommending an extremely high
        default on a big server would not meaningfully speed things up
        and could cause resource contention if several tools/samples run
        concurrently. Users with genuinely large machines can still
        manually raise the slider past this recommended default if they
        want to experiment.

    Returns (total_cores_detected: int, recommended_default: int).
    total_cores_detected is also returned (not just the recommendation)
    so the UI can show the user *why* a particular default was chosen
    (e.g. "Detected 8 CPU cores on this machine").

    Falls back to 1 core detected / 1 recommended if os.cpu_count()
    can't determine the core count at all (returns None) — this can
    happen in some restricted/containerized environments — rather than
    crashing or recommending an arbitrary guessed value.
    """
    total_cores = os.cpu_count() or 1
    recommended = max(1, total_cores - reserve_cores)
    recommended = min(recommended, max_default)
    return total_cores, recommended




def deseq2_dir(project_name):
    """Where DESeq2 analysis inputs/outputs live for this project."""
    return os.path.join(project_dir(project_name), "deseq2")


def deseq2_output_dir(project_name):
    """DESeq2's actual result CSVs (per-contrast results, PCA, normalized counts)."""
    return os.path.join(deseq2_dir(project_name), "output")


def deseq2_work_dir(project_name):
    """Scratch directory for the temporary R script + job spec JSON."""
    return os.path.join(deseq2_dir(project_name), "work")


def save_deseq2_config(project_name, config):
    """
    Remember the user's DESeq2 configuration (design columns, batch
    column, contrasts, filtering thresholds) so re-opening the project
    doesn't require re-entering all of it. `config` should be a plain
    JSON-serializable dict.
    """
    info = load_info(project_name)
    info["deseq2_config"] = config
    save_info(project_name, info)


def get_deseq2_config(project_name):
    info = load_info(project_name)
    return info.get("deseq2_config")


def save_reference_choice(project_name, species_key, is_custom):
    """
    Remember which reference organism (or "custom") was selected for
    this project, so re-opening it doesn't require re-selecting.
    """
    info = load_info(project_name)
    info["reference_species"] = species_key
    info["reference_is_custom"] = is_custom
    save_info(project_name, info)


def get_reference_choice(project_name):
    info = load_info(project_name)
    return info.get("reference_species"), info.get("reference_is_custom", False)


def save_alignment_method(project_name, method):
    """
    Remember which alignment/quantification method ("salmon" or "star")
    the user chose for this project, so re-opening it doesn't require
    re-selecting — the reference setup and execution steps differ
    significantly between the two methods.
    """
    info = load_info(project_name)
    info["alignment_method"] = method
    save_info(project_name, info)


def get_alignment_method(project_name):
    info = load_info(project_name)
    return info.get("alignment_method")


def info_path(project_name):
    return os.path.join(project_dir(project_name), "project_info.json")


# ---------------------------------------------------------------------------
# Project lifecycle
# ---------------------------------------------------------------------------

def create_project(project_name):
    """
    Create a new project's folder structure. Returns True if created,
    False if a project with that name already exists.
    """
    d = project_dir(project_name)
    if os.path.exists(d):
        return False
    os.makedirs(fastq_dir(project_name), exist_ok=True)
    os.makedirs(fastqc_dir(project_name), exist_ok=True)
    os.makedirs(multiqc_dir(project_name), exist_ok=True)
    save_info(project_name, {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "steps_completed": [],
    })
    return True


def load_info(project_name):
    path = info_path(project_name)
    if not os.path.exists(path):
        return {"created_at": None, "steps_completed": []}
    with open(path) as f:
        return json.load(f)


def save_info(project_name, info):
    os.makedirs(project_dir(project_name), exist_ok=True)
    with open(info_path(project_name), "w") as f:
        json.dump(info, f, indent=2)


def mark_step_complete(project_name, step_name):
    """
    Record that a workflow step (e.g. 'samples_matched', 'qc_complete',
    'trimming_complete', 'reference_ready') has been completed for this
    project. Future steps/modules can check this to know what stage a
    project is at.
    """
    info = load_info(project_name)
    if step_name not in info.get("steps_completed", []):
        info.setdefault("steps_completed", []).append(step_name)
        info["last_updated"] = datetime.now().isoformat(timespec="seconds")
    save_info(project_name, info)


def has_completed_step(project_name, step_name):
    info = load_info(project_name)
    return step_name in info.get("steps_completed", [])


def save_sample_column(project_name, column_name):
    """
    Remember which column the user picked as the 'sample name' column in
    their metadata file, so re-opening the project doesn't require
    re-selecting it (though the user can still change it).
    """
    info = load_info(project_name)
    info["metadata_sample_column"] = column_name
    save_info(project_name, info)


def get_sample_column(project_name):
    info = load_info(project_name)
    return info.get("metadata_sample_column")


# ---------------------------------------------------------------------------
# Reusable UI: project creation / selection widget
# ---------------------------------------------------------------------------

def render_project_selector(workspace_key="bulk_rnaseq"):
    """
    Render a "Step 0" project creation/selection widget. Returns the
    selected project name (str), or None if no project is active yet.

    The selection is stored in
    st.session_state[f"{workspace_key}_project"] so it persists across
    reruns within the session, and is also recoverable across sessions
    since the project folder itself lives on disk.
    """
    st.header("Step 0: Set Up Your Project")
    st.markdown(
        "A **project** keeps all your files together — uploaded reads, "
        "matched samples, and quality control results — in one place. "
        "This lets you pick up exactly where you left off, and lets "
        "later steps (like adapter trimming) automatically find the "
        "right files without needing to re-upload anything."
    )

    session_key = f"{workspace_key}_project"
    existing = list_projects()

    mode = st.radio(
        "Would you like to start a new project or continue an existing one?",
        ["🆕 Start a new project", "📂 Continue an existing project"],
        horizontal=True,
        index=0 if not existing else 1,
        key=f"{workspace_key}_project_mode",
    )

    if mode == "🆕 Start a new project":
        new_name = st.text_input(
            "Name your project (e.g. 'LiverStudy2026' or 'TCellExperiment'):",
            help="Use letters, numbers, dashes, or underscores — spaces and special characters will be removed automatically.",
            key=f"{workspace_key}_new_project_name",
        )
        if new_name:
            safe_name = "".join(c for c in new_name.strip() if c.isalnum() or c in ("-", "_"))
            if safe_name != new_name.strip():
                st.caption(f"Note: your project will be saved as `{safe_name}` (spaces/special characters removed).")

            if st.button("✅ Create Project", key=f"{workspace_key}_create_project_btn"):
                if not safe_name:
                    st.error("Please enter a valid project name (letters, numbers, dashes, or underscores).")
                elif safe_name in existing:
                    st.error(
                        f"A project named `{safe_name}` already exists. "
                        "Choose 'Continue an existing project' above instead, "
                        "or pick a different name."
                    )
                else:
                    create_project(safe_name)
                    st.session_state[session_key] = safe_name
                    st.success(f"Project `{safe_name}` created!")
                    st.rerun()
    else:
        if not existing:
            st.info("No existing projects found yet. Choose 'Start a new project' above to create one.")
        else:
            chosen = st.selectbox(
                "Select a project:",
                options=existing,
                key=f"{workspace_key}_existing_project_select",
            )
            if st.button("📂 Open Project", key=f"{workspace_key}_open_project_btn"):
                st.session_state[session_key] = chosen
                st.rerun()

    selected_project = st.session_state.get(session_key)

    if selected_project:
        info = load_info(selected_project)
        st.success(f"**Active project: `{selected_project}`**")
        completed = info.get("steps_completed", [])
        if completed:
            st.caption(f"✔️ Steps completed so far: {', '.join(completed)}")
        if st.button("🔄 Switch to a different project", key=f"{workspace_key}_switch_project_btn"):
            del st.session_state[session_key]
            st.rerun()

    return selected_project
