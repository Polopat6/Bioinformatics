"""
project_actions.py

Shared "completed project" actions used by both Auto Mode's and Monitor
Mode's launched-projects lists (advanced_mode_workspace.py /
monitor_mode_workspace.py). Kept here -- rather than duplicated in each
workspace -- so QC/trimming/alignment/counts packaging and the DESeq2
hand-off behave identically regardless of which mode launched the run.

--- Schema update (2026-08-16) ---
This module was originally written without project_manager.py's exact
on-disk directory schema in hand, so _resolve_project_dir() and a
PACKAGE_SUBDIR_CANDIDATES guessing table checked several candidate
names/attributes and picked whichever existed. project_manager.py's
real schema has now been confirmed (project_dir(), fastqc_dir(),
multiqc_dir(), posttrim_multiqc_dir(), trimmed_fastq_dir(),
fastp_reports_dir(), salmon_quant_dir(), star_align_dir(),
counts_matrix_path(), get_alignment_method(), info_path()) -- every
path below now calls those functions directly instead of guessing, so
this module always matches whatever project_manager.py actually
creates on disk, with no silent "guessed wrong, packaged nothing"
failure mode.

Note: counts_matrix_path() returns a FILE (the combined gene counts
matrix CSV), not a directory -- project_package_available() and
build_project_package_zip() both handle it as a file, unlike every
other path here.

--- Package scope (2026-08-16) ---
By design, the downloadable package includes ONLY: FastQC + MultiQC
(pre- and post-trim) reports, fastp's per-sample QC reports, alignment
scores/info (log/summary files from the Salmon or STAR output
directory -- e.g. Log.final.out/SJ.out.tab/ReadsPerGene.out.tab for
STAR, or the meta_info.json/lib_format_counts.json/logs for Salmon),
and the final gene counts matrix. It deliberately EXCLUDES all raw and
trimmed FASTQ read files and any aligned-read files (BAM/SAM/CRAM) --
see _READ_FILE_EXTENSIONS / _is_read_file() below -- since those are
large sequencing data files already living in the project's own
folders on disk, not "results" a user would want bundled into a
lightweight results/QC download.

--- Navigation hand-off update (2026-08-16) ---
app.py's real router has now been confirmed: request_navigation_to_deseq2()
below sets the actual "nav_request" session_state key app.py's top-level
script checks (popping it, setting "active_workspace", and syncing the
right pipeline drawer's radio widget) -- replacing this function's earlier
placeholder keys ("active_workspace_step" / "sidebar_radio_bulk_rnaseq"),
which never matched anything in app.py.

--- Project pre-selection CONFIRMED (2026-08-16) ---
project_manager.render_project_selector()'s real implementation has now
been inspected. It stores the active project in a PLAIN (non-widget-bound)
session_state key -- st.session_state[f"{workspace_key}_project"] -- that
is only ever set inside that function via a button click (Create Project /
Open Project), or, as done here, externally before the function runs.
Since it is NOT a widget's own key (no st.selectbox/st.text_input/etc. is
declared with this exact key), Streamlit's "cannot set a widget-bound key
after that widget is instantiated" restriction does NOT apply -- setting it
ahead of time is safe and is exactly what render_project_selector() checks
first (`selected_project = st.session_state.get(session_key)`) before
falling back to its own "start new" / "continue existing" UI.
differential_expression_workspace.py calls
render_project_selector(workspace_key="bulk_rnaseq") (its own
WORKSPACE_KEY = "bulk_rnaseq"), so the real, confirmed key to set is
"bulk_rnaseq_project" -- replacing the earlier placeholder
"pending_preselect_project" key that render_project_selector() never
actually read.
"""
import io
import json
import os
import zipfile

import project_manager as pm

# Real, confirmed project_manager.py path helpers this module packages up.
# Each entry maps an arcname prefix (used inside the zip) -> the actual
# directory on disk. Built once per call (rather than a fixed module-level
# dict) since the alignment-stage entry depends on which method
# (Salmon vs. STAR) this specific project used.
#
# --- Scope update (2026-08-16): reads excluded from the package ---
# The package is meant to be a lightweight "results + QC" bundle a user can
# email/archive/hand off -- NOT a full copy of the project's sequencing
# reads (those already live in the project's own fastq/ and trimmed/
# folders on disk, and re-zipping full FASTQ data, raw or trimmed, would
# make this "download package" button produce a multi-GB file for a
# typical RNA-seq run). Two changes enforce this:
#   1. trimmed_fastq_dir() (the trimmed FASTQ reads themselves) has been
#      removed from _package_dir_map() entirely -- it was the one entry
#      here that held actual sequencing reads.
#   2. _READ_FILE_EXTENSIONS below is checked while walking every
#      remaining directory (including the alignment directory), so any
#      raw/aligned read data that lands alongside legitimate summary
#      files -- e.g. STAR's Aligned.sortedByCoord.out.bam sitting right
#      next to its own Log.final.out / SJ.out.tab / ReadsPerGene.out.tab
#      summary files in the same per-sample folder -- is skipped, while
#      those summary/log/score files ARE still included (this is what
#      "alignment scores and info" resolves to on disk: everything in
#      the alignment directory except the alignment file itself).
#      raw_fastq_dir() (pm.fastq_dir()) was never included in this map to
#      begin with, so no separate change was needed there.
_READ_FILE_EXTENSIONS = (
    ".fastq", ".fastq.gz", ".fq", ".fq.gz",  # raw/trimmed sequencing reads
    ".bam", ".sam", ".cram",                  # aligned reads
)


def _resolve_project_dir(project_name):
    "This project's root directory on disk -- see project_manager.py's own project_dir()."
    return pm.project_dir(project_name)


def _is_read_file(filename):
    "True if filename looks like a raw/trimmed/aligned sequencing read file that should be excluded from the download package."
    lower = filename.lower()
    return lower.endswith(_READ_FILE_EXTENSIONS)


def _package_dir_map(project_name):
    """
    Build this project's {arcname_prefix: source_dir} map for zipping,
    using project_manager.py's real path helpers directly. Only
    QC/report directories, the alignment directory (scores/info, not the
    aligned reads themselves -- see _is_read_file), and the final counts
    matrix (added separately in build_project_package_zip) are included.
    The alignment entry is chosen based on whichever method
    (get_alignment_method) this project actually used -- a project
    only ever has ONE of salmon_quant_dir()/star_align_dir() populated,
    never both.
    """
    alignment_method = pm.get_alignment_method(project_name)
    alignment_dir = (
        pm.star_align_dir(project_name) if alignment_method == "star"
        else pm.salmon_quant_dir(project_name)
    )
    return {
        "qc/fastqc": pm.fastqc_dir(project_name),
        "qc/multiqc": pm.multiqc_dir(project_name),
        "qc/multiqc_posttrim": pm.posttrim_multiqc_dir(project_name),
        "qc/fastp": pm.fastp_reports_dir(project_name),
        f"quant/{alignment_method or 'alignment'}": alignment_dir,
    }


def project_package_available(project_name):
    """
    Quick, cheap check for whether a project has at least a finished
    gene counts matrix on disk -- used to decide whether to show the
    download-package button at all (vs. a run that errored out before
    completion). counts_matrix_path() is a FILE, not a directory, so
    this checks os.path.isfile rather than isdir.
    """
    try:
        return os.path.isfile(pm.counts_matrix_path(project_name))
    except Exception:
        return False


def build_project_package_zip(project_name):
    """
    Bundle a completed project's FastQC, MultiQC (pre- and post-trim),
    trimming (fastp reports only, not the trimmed reads themselves),
    alignment scores/info (Salmon or STAR, whichever this project used
    -- log/summary files only, not the aligned reads themselves), and
    the final gene counts matrix -- plus a project_info.json summary of
    recorded step-completion status -- into a single in-memory zip
    file.

    Returns raw zip bytes, ready to hand straight to
    st.download_button(data=...).
    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        for arcname_prefix, src_dir in _package_dir_map(project_name).items():
            if not src_dir or not os.path.isdir(src_dir):
                continue
            for root, _dirs, files in os.walk(src_dir):
                for fname in files:
                    if _is_read_file(fname):
                        # Skip raw/trimmed/aligned sequencing reads (e.g. a
                        # BAM file sitting in the alignment directory
                        # alongside its own Log.final.out/SJ.out.tab/
                        # ReadsPerGene.out.tab summary files) -- only the
                        # scores/info/logs from that same directory are
                        # packaged, per this module's "QC + MultiQC +
                        # alignment scores/info + counts table only" scope.
                        continue
                    abs_path = os.path.join(root, fname)
                    rel_path = os.path.relpath(abs_path, src_dir)
                    zf.write(abs_path, arcname=os.path.join(arcname_prefix, rel_path))
        counts_path = pm.counts_matrix_path(project_name)
        if os.path.isfile(counts_path):
            zf.write(counts_path, arcname=os.path.join("quant", os.path.basename(counts_path)))
        info = pm.load_info(project_name)
        zf.writestr("project_info.json", json.dumps(info, indent=2))
    buffer.seek(0)
    return buffer.getvalue()


# ---------------------------------------------------------------------------
# "Begin DE Analysis" hand-off to the DESeq2 workspace
# ---------------------------------------------------------------------------
# CONFIRMED (2026-08-16) against app.py's real router: app.py's top-level
# script (not inside any function, so fully visible/verified) does:
#
#   if "nav_request" in st.session_state:
#       target = st.session_state.pop("nav_request")
#       st.session_state["active_workspace"] = target
#       for group_key, group in PIPELINE_GROUPS.items():
#           if target in group["options"]:
#               st.session_state[f"_workspace_radio_{group_key}"] = target
#               break
#
# So the ONLY thing a caller needs to do to navigate to a workspace is set
# st.session_state["nav_request"] to that workspace's EXACT option string
# from PIPELINE_GROUPS (app.py's own bulk_rnaseq_workspace.py "Proceed to
# Trimming" button uses this same mechanism) and then st.rerun(). The
# previous "NAV_TARGET_SESSION_KEY = active_workspace_step" /
# "sidebar_radio_bulk_rnaseq" keys in an earlier draft of this function were
# placeholder GUESSES that do not exist anywhere in app.py -- replaced here
# with the real "nav_request" key and the real, exact
# "🌋 Differential Expression" option string (NOT "DESeq2").
DESEQ2_NAV_TARGET = "🌋 Differential Expression"

# CONFIRMED (2026-08-16) against project_manager.render_project_selector()'s
# real implementation: differential_expression_workspace.py calls
# render_project_selector(workspace_key="bulk_rnaseq"), and that function
# stores/reads the active project via the plain (non-widget) session_state
# key f"{workspace_key}_project" -- i.e. "bulk_rnaseq_project" here. Setting
# this key before render_project_selector() runs makes
# `selected_project = st.session_state.get(session_key)` immediately truthy,
# so the DESeq2 page opens with this project already active (showing its
# "Active project" banner directly) instead of requiring the user to
# manually pick "Continue an existing project" and select it themselves.
NAV_PROJECT_SESSION_KEY = "bulk_rnaseq_project"


def request_navigation_to_deseq2(project_name):
    """
    Set the session_state flags needed to land the user on the
    Differential Expression (DESeq2) workspace with `project_name`
    already active: app.py's real "nav_request" router (see above) for
    the page navigation itself, and project_manager's real
    f"{workspace_key}_project" key (NAV_PROJECT_SESSION_KEY =
    "bulk_rnaseq_project") so render_project_selector() picks up this
    project immediately rather than showing its own "start new/continue
    existing" chooser. Caller should st.rerun() immediately after calling
    this.
    """
    import streamlit as st
    st.session_state[NAV_PROJECT_SESSION_KEY] = project_name
    st.session_state["nav_request"] = DESEQ2_NAV_TARGET
