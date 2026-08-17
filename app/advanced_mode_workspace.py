"""
advanced_mode_workspace.py

Streamlit wizard UI for "Auto" -- one of possibly several workspaces
living in the top-level "Advanced Modes" sidebar drawer (see app.py's
PIPELINE_GROUPS; Monitor Mode joins this same drawer via its own
monitor_mode_workspace.py). Auto's own FIRST screen is a pipeline
picker (see PIPELINE_AUTO_HANDLERS / render() below) -- today that
picker only offers Bulk RNA-Seq, but the dispatch structure exists so
that Single-cell RNA-seq / Spatial RNA-seq can register their own Auto
handler here later without restructuring this file, exactly mirroring
how app.py's PIPELINE_GROUPS itself is structured to grow one entry at
a time.

For whichever pipeline is selected, this defines that project's FASTQ
source, metadata, reference genome, and pipeline options in ONE guided
screen, then launches the whole thing (for Bulk RNA-Seq: FASTQ
ingestion -> QC -> trimming -> reference setup -> quantification ->
gene counts matrix) as a single resumable BACKGROUND process via
advanced_mode_orchestrator.py, instead of walking through each step
interactively. DESeq2 contrasts and Ontology Analysis remain separate,
interactive steps performed afterward, exactly as in the normal
step-by-step workflow -- this wizard/orchestrator stop at the finished
counts matrix.

--- Why NCBI/SRA validation happens HERE, synchronously, in the UI ---
sra_manager.lookup_accession() / lookup_multiple_accessions() are cheap
esearch/efetch metadata calls (seconds, not hours) -- there's no reason
to defer these into the background run. Doing them here, before the
project's config is even saved, means a typo'd or non-existent
accession is caught immediately, and the resulting SAMPLE_ATTRIBUTES
can auto-populate the metadata table for the user to review -- exactly
the same NCBI record data the interactive workspace's SRA lookup
feature surfaces today, just consolidated into this one upfront screen.
The actual FASTQ DOWNLOAD (prefetch + fasterq-dump), by contrast, is
genuinely long-running (minutes to multiple hours per run depending on
file size and NCBI server load) -- that work is deferred entirely to
the background pipeline process (advanced_mode_orchestrator._run_ingest,
stage 1 of the run), which is why a prominent time-expectation warning
is shown right before the launch button whenever the FASTQ source is
NCBI/SRA -- see _ncbi_download_warning() below.

--- Post-completion behavior (2026-08-16) ---
Previously, once a project's background run finished, this workspace
kept rendering the full Steps 1-4 configuration wizard right below the
"Pipeline complete" status panel -- as if inviting the user to
reconfigure/relaunch a project that was already done, with no
indication that anything had actually finished. Now, once
orch.get_status(project)["pipeline_status"] == "complete", the wizard
is replaced with a completed-project actions panel (via
project_actions.py, the same shared actions Monitor Mode's
launched-projects list uses): a "Prepare/Download Project Package"
button (QC/MultiQC/alignment scores/counts table -- no raw or aligned
reads, see project_actions.py) and a "Begin DE Analysis" button that
hands off straight to the Differential Expression workspace with this
project selected. Steps 1-4 are still reachable on demand via a
"Reconfigure & Re-launch This Project" expander, for the (presumably
rarer) case of wanting to re-run a project that's already complete.
"""
import os

import pandas as pd
import streamlit as st

import project_manager as pm
import sra_manager as sra
import reference_manager as ref
import file_browser as fb
import ingestion_manager as ingest
import advanced_mode_orchestrator as orch
import project_actions as pa

STAGE_LABELS = {
    "ingest": "1. FASTQ ingestion + sample matching",
    "qc": "2. Pre-trim QC (FastQC + MultiQC)",
    "trimming": "3. Trimming (fastp) + post-trim QC",
    "reference": "4. Reference genome + index setup",
    "quantification": "5. Quantification (Salmon/STAR)",
    "counts_matrix": "6. Gene counts matrix",
}
STAGE_ICONS = {
    "pending": "⏳", "running": "🔄", "complete": "✅", "skipped": "⏭️", "error": "❌",
}
# quantification_manager.run_tximport_gene_collapse() collapses Salmon's
# transcript-level quant.sf output to gene-level counts (via tximport's
# countsFromAbundance="lengthScaledTPM"), which matters for any organism
# with alternative splicing (multiple transcripts per gene) -- see that
# function's module-level comment for the full rationale. Enabled by
# default; the orchestrator's counts_matrix stage automatically falls
# back to the direct transcript-level merge if tximport/R aren't
# available on this system, or if no tx2gene mapping exists for the
# selected reference (e.g. STAR runs, which are already gene-level and
# don't need this at all).
_USE_TXIMPORT = True


def _init_state():
    for key, default in {
        "adv_sra_lookup_rows": None,
        "adv_sra_lookup_message": None,
    }.items():
        st.session_state.setdefault(key, default)


def _ncbi_download_warning(n_runs):
    return (
        f"⚠️ **Downloading from NCBI/SRA can take a while.** Fetching "
        f"{n_runs} run(s) could take anywhere from a few minutes to "
        f"several hours *per run*, depending on file size and current "
        f"NCBI server load -- large human RNA-seq runs in particular can "
        f"take an hour or more each. This download happens as **stage 1 "
        f"of the background run** you're about to launch: once started, "
        f"you're free to close this tab and check back later -- progress "
        f"is shown in the status panel below whenever you return."
    )


# ---------------------------------------------------------------------------
# Step 1: FASTQ source
# ---------------------------------------------------------------------------
def _render_fastq_source_step(project):
    st.subheader("Step 1: FASTQ Source")
    fastq_dir = pm.fastq_dir(project)
    existing = ingest.list_existing_fastq(fastq_dir)
    if existing:
        st.success(f"✅ {len(existing)} FASTQ file(s) already present in this project.")
    source = st.radio(
        "How will you provide FASTQ files for this project?",
        ["📤 Upload from my computer", "📂 Browse a directory on this server", "🔎 Fetch from NCBI/SRA"],
        key="adv_fastq_source_radio", horizontal=True,
    )
    if source.startswith("📤"):
        uploaded_files = st.file_uploader(
            "Upload your raw FASTQ file(s):", type=["fastq", "gz", "fq"],
            accept_multiple_files=True, key="adv_fastq_upload",
        )
        if uploaded_files and st.button("💾 Save Uploaded Files", key="adv_fastq_upload_save_btn"):
            saved = ingest.save_uploaded_files(uploaded_files, fastq_dir)
            st.success(f"✅ Saved {len(saved)} file(s) to this project.")
            st.rerun()
        # Files land directly in the project's own fastq_dir, so the
        # background run's "directory" ingest path is pointed at that
        # SAME directory -- ingestion_manager.symlink_fastq_files_from_directory
        # safely no-ops (skips) any file that's already there, so this
        # is just a harmless confirmation pass, not a real copy/link.
        return {"fastq_source": "directory", "fastq_source_dir": fastq_dir}
    elif source.startswith("📂"):
        selected_dir = fb.render_server_directory_browser(
            key_prefix="adv_fastq_dir_browse",
            preview_extensions=[".fastq", ".fastq.gz", ".fq", ".fq.gz"],
            label="Browse for the directory containing your FASTQ files:",
        )
        if selected_dir:
            _, fastq_filenames = ingest.find_fastq_filenames_in_directory(selected_dir)
            if fastq_filenames:
                st.success(f"✅ Found {len(fastq_filenames)} FASTQ file(s) in this directory.")
                return {"fastq_source": "directory", "fastq_source_dir": selected_dir}
            st.warning("⚠️ No FASTQ files found directly inside this directory.")
        return None
    else:
        return _render_ncbi_source(project)


def _render_ncbi_source(project):
    """
    Validate NCBI/SRA accession(s) SYNCHRONOUSLY (cheap esearch/efetch),
    right here in the UI -- see module docstring. The actual download
    is deferred to the background pipeline run (stage "ingest").
    """
    st.markdown(
        "Enter a **study/BioProject accession** (e.g. `SRP123456` or "
        "`PRJNA123456`), and/or paste **individual run accessions** "
        "(e.g. `SRR1234567`, separated by commas/spaces/newlines)."
    )
    prefetch_ok, fasterq_ok = sra.tools_available()
    if not (prefetch_ok and fasterq_ok):
        st.error("⚠️ The SRA Toolkit (`prefetch`/`fasterq-dump`) isn't available on this system yet.")
        return None
    single_term_input = st.text_input(
        "Study or BioProject accession:", placeholder="e.g. SRP123456 or PRJNA123456",
        key="adv_sra_single_term_input",
    )
    accession_list_text = st.text_area(
        "Or paste one or more individual run accessions:",
        placeholder="e.g.\nSRR1039508, SRR1039509, SRR1039512",
        key="adv_sra_accession_list_input", height=90,
    )
    uploaded_accession_file = st.file_uploader(
        "Or upload a file listing accessions (.txt, .csv, or .xlsx):",
        type=["txt", "csv", "xlsx", "xls"], key="adv_sra_accession_file_upload",
    )
    if st.button("🔍 Validate Accession(s) on NCBI", key="adv_sra_validate_btn"):
        combined_accessions = list(sra._split_accession_list_text(accession_list_text))
        file_read_error = None
        if uploaded_accession_file is not None:
            file_accessions, file_read_error = sra.parse_accessions_from_file(uploaded_accession_file)
            combined_accessions.extend(file_accessions)
        seen = set()
        combined_accessions = [a for a in combined_accessions if not (a in seen or seen.add(a))]
        if file_read_error:
            st.error(f"⚠️ {file_read_error}")
        elif combined_accessions:
            with st.spinner(f"Validating {len(combined_accessions)} accession(s) on NCBI..."):
                success, rows, message, not_found = sra.lookup_multiple_accessions(combined_accessions)
            if not success:
                st.error(f"⚠️ {message}")
                st.session_state["adv_sra_lookup_rows"] = None
            else:
                st.session_state["adv_sra_lookup_rows"] = rows
                st.success(f"✅ {message} — all valid, ready to download.")
                if not_found:
                    st.warning(f"⚠️ Not found on NCBI (check for typos): {', '.join(not_found)}")
        elif single_term_input:
            with st.spinner(f"Validating '{single_term_input}' on NCBI..."):
                success, rows, message = sra.lookup_accession(single_term_input)
            if not success:
                st.error(f"⚠️ {message}")
                st.session_state["adv_sra_lookup_rows"] = None
            else:
                st.session_state["adv_sra_lookup_rows"] = rows
                st.success(f"✅ {message} — valid, ready to download.")
        else:
            st.warning("⚠️ Please enter an accession, paste a list, or upload a file to validate.")
    lookup_rows = st.session_state.get("adv_sra_lookup_rows")
    if not lookup_rows:
        return None
    display_rows = [{
        "Run": row.get("Run", "—"),
        "Organism": row.get("ScientificName", "—"),
        "Library Strategy": row.get("LibraryStrategy", "unknown") if sra.is_rna_seq(row)
                            else f"⚠️ {row.get('LibraryStrategy', 'unknown')} (NOT RNA-Seq)",
        "Size (MB)": row.get("size_MB", "—"),
    } for row in lookup_rows]
    st.dataframe(pd.DataFrame(display_rows), use_container_width=True, hide_index=True)
    non_rna_count = sum(1 for r in lookup_rows if not sra.is_rna_seq(r))
    if non_rna_count:
        st.warning(f"⚠️ {non_rna_count} of {len(lookup_rows)} run(s) are not annotated as RNA-Seq.")
    run_options = [row["Run"] for row in lookup_rows]
    selected_runs = st.multiselect(
        "Select which validated run(s) to download:", options=run_options,
        default=run_options, key="adv_sra_selected_runs",
    )
    if not selected_runs:
        return None
    # Same "parallel downloads" control the interactive workspace's own
    # SRA lookup section offers (bulk_rnaseq_workspace.py's
    # _render_sra_lookup_section) -- since downloading is mostly spent
    # waiting on network transfer rather than heavy local computation,
    # fetching several runs at once is meaningfully faster than one at
    # a time. This value is threaded through the saved config
    # ("sra_max_workers") and read by the background pipeline's ingest
    # stage (advanced_mode_orchestrator._run_ingest) at download time --
    # it has no effect on anything in this synchronous validation UI
    # itself.
    detected_cores, _ = pm.get_recommended_thread_count()
    max_workers = st.slider(
        "Parallel downloads:", min_value=1, max_value=detected_cores,
        value=min(3, detected_cores), key="adv_sra_max_workers",
        help=(
            "How many runs to download at the same time during the "
            "background pipeline's FASTQ ingestion stage. A moderate "
            "default is suggested — increase it if your internet "
            "connection is fast and you're downloading many runs; "
            "decrease it if downloads start failing."
        ),
    )
    st.info(_ncbi_download_warning(len(selected_runs)))
    return {"fastq_source": "sra", "sra_accessions": selected_runs, "sra_max_workers": max_workers}


# ---------------------------------------------------------------------------
# Step 2: Metadata
# ---------------------------------------------------------------------------
def _render_metadata_step(project, fastq_config):
    st.subheader("Step 2: Sample Metadata")
    lookup_rows = st.session_state.get("adv_sra_lookup_rows")
    selected_runs = fastq_config.get("sra_accessions") if fastq_config else None
    if fastq_config and fastq_config.get("fastq_source") == "sra" and lookup_rows and selected_runs:
        st.markdown(
            "Auto-built from each validated sample's own NCBI record -- "
            "review/edit before continuing (works even before the "
            "download itself has run)."
        )
        default_df = pd.DataFrame(sra.build_metadata_dataframe(lookup_rows, selected_runs=selected_runs))
        if not sra.has_any_descriptive_metadata(lookup_rows, selected_runs=selected_runs):
            st.warning("⚠️ NCBI provided no characteristics beyond the run accession — add your own columns below.")
    else:
        input_method = st.radio(
            "How would you like to provide your metadata file?",
            ["📤 Upload from your computer", "📂 Browse files already on this server"],
            key="adv_metadata_input_method_radio", horizontal=True,
        )
        if input_method.startswith("📤"):
            metadata_file = st.file_uploader(
                "Upload your sample metadata file (.csv or .xlsx):",
                type=["csv", "txt", "xlsx", "xls"], key="adv_metadata_upload",
            )
        else:
            metadata_file = fb.render_server_file_browser(
                key_prefix="adv_metadata_browse", file_extensions=[".csv", ".txt", ".xlsx", ".xls"],
                label="Browse for your sample metadata file already on this server:",
            )
        default_df = None
        if metadata_file is not None:
            parsed_df, error = ingest.read_metadata_file(metadata_file)
            if error:
                st.error(error)
            else:
                default_df = parsed_df
    if default_df is None:
        st.info("Provide a metadata source above to continue.")
        return None
    if "sample" not in default_df.columns:
        st.error("⚠️ Your metadata must have a column named exactly `sample` matching FASTQ-derived sample names.")
    return st.data_editor(default_df, num_rows="dynamic", use_container_width=True, key="adv_metadata_editor")


# ---------------------------------------------------------------------------
# Step 3: Genome & pipeline options
# ---------------------------------------------------------------------------
def _render_genome_and_options_step(project):
    st.subheader("Step 3: Reference Genome & Pipeline Options")
    is_custom = st.checkbox("Use a custom (non-preset) reference genome", key="adv_ref_is_custom")
    reference_cfg = {"is_custom": is_custom}
    if not is_custom:
        species_options = list(ref.REFERENCE_CATALOG.keys())
        species_labels = {k: v["label"] for k, v in ref.REFERENCE_CATALOG.items()}
        species_choice = st.selectbox(
            "Reference organism:", options=species_options,
            format_func=lambda k: species_labels[k], key="adv_species_choice",
        )
        reference_cfg["species_key"] = species_choice
    else:
        st.caption("Provide a genome FASTA and GTF already on this server (used for STAR, and for Salmon unless a pre-extracted transcript FASTA is given).")
        genome_fasta = fb.render_server_file_browser(
            key_prefix="adv_ref_genome_browse", file_extensions=[".fa", ".fasta", ".fa.gz", ".fna"],
            label="Browse for the genome FASTA:",
        )
        gtf_path = fb.render_server_file_browser(
            key_prefix="adv_ref_gtf_browse", file_extensions=[".gtf", ".gtf.gz"],
            label="Browse for the GTF annotation:",
        )
        reference_cfg.update({
            "species_key": None, "custom_genome_fasta": genome_fasta,
            "custom_gtf": gtf_path, "custom_transcript_fasta": None,
        })
    alignment_method = st.radio(
        "Alignment/quantification method:", ["salmon", "star"],
        format_func=lambda m: "Salmon (pseudo-alignment, faster)" if m == "salmon" else "STAR (splice-aware alignment)",
        key="adv_alignment_method", horizontal=True,
    )
    # Same three options + same help text as alignment_workspace.py's
    # "⚙️ Advanced STAR options (optional)" expander -- previously
    # missing from this wizard entirely, so a STAR run launched via
    # Advanced Mode always silently used the defaults (all three off)
    # with no way to turn them on. Only shown when STAR is selected,
    # since these flags have no meaning for a Salmon run.
    star_options_cfg = {"use_two_pass": False, "use_encode_options": False, "add_strand_field": False}
    if alignment_method == "star":
        with st.expander("⚙️ Advanced STAR options (optional)"):
            star_options_cfg["use_two_pass"] = st.checkbox(
                "Two-pass mode (better novel splice junction detection)",
                value=False, key="adv_star_two_pass_checkbox",
                help=(
                    "STAR aligns your reads once to discover splice junctions, "
                    "then re-aligns using those junctions as extra annotation -- "
                    "meaningfully improves detection of splice junctions not "
                    "already in your GTF. Most useful for less-annotated/"
                    "non-model organisms; roughly doubles alignment time per "
                    "sample. Well-annotated genomes (human/mouse) benefit less "
                    "since most real junctions are usually already known."
                ),
            )
            star_options_cfg["use_encode_options"] = st.checkbox(
                "Use ENCODE-recommended filtering options",
                value=False, key="adv_star_encode_options_checkbox",
                help=(
                    "Applies the standard set of alignment filtering settings "
                    "used in the ENCODE project's own RNA-seq pipelines "
                    "(multimapping tolerance, splice junction overhang "
                    "minimums, mismatch tolerance, intron length bounds). A "
                    "well-established, commonly-used preset bundle rather "
                    "than a single tunable setting."
                ),
            )
            star_options_cfg["add_strand_field"] = st.checkbox(
                "Add strand field to BAM output",
                value=False, key="adv_star_strand_field_checkbox",
                help=(
                    "Annotates each alignment's strand in the output BAM file "
                    "based on intron motif. Only needed if you plan to feed "
                    "this project's BAM files into a downstream tool that "
                    "expects this strand tag (e.g. Cufflinks-family tools) -- "
                    "has no effect on the gene counts this app itself "
                    "produces, so most users can safely leave this unchecked."
                ),
            )
    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    threads = st.slider(
        "Threads for each tool (FastQC/fastp/Salmon/STAR):",
        min_value=1, max_value=detected_cores, value=recommended_threads, key="adv_pipeline_threads",
    )
    thread_cfg = {k: threads for k in ("fastqc", "fastp", "salmon_index", "salmon_quant", "star_index", "star_align")}
    # Placed right next to the thread-count picker since both affect
    # HOW the trimming (fastp) stage runs, rather than WHAT data goes
    # in -- controls whether the background pipeline's trimming stage
    # automatically detects AND fixes residual poly-G/poly-A tails
    # after the main fastp trim (see
    # advanced_mode_orchestrator._run_trimming) -- the same detect-
    # then-retrim logic the interactive Trimming workspace exposes via
    # a manual "Auto Re-trim Flagged Sample(s)" button, just fully
    # automatic here since there's no user present to click anything
    # during an unattended background run.
    auto_fix_poly_tails = st.toggle(
        "Automatically detect and fix residual poly-G/poly-A tails",
        value=True, key="adv_auto_fix_poly_tails",
    )
    if auto_fix_poly_tails:
        st.caption(
            "✅ **Recommended (on by default).** fastp does **not** trim poly-G "
            "or poly-A tails by default — poly-G trimming only self-enables "
            "for reads it recognizes as NextSeq/NovaSeq-style 2-color "
            "chemistry, and poly-A trimming is never automatic at all, even "
            "though poly-A tail read-through is common in ordinary mRNA-seq "
            "libraries. Left unchecked, these tails can inflate low-quality/"
            "low-complexity read fractions and show up as confusing "
            "\"overrepresented sequence\" warnings later. With this on, the "
            "background run automatically re-checks every sample's trimmed "
            "reads for this specific issue and re-trims only the affected "
            "sample(s) with the correct fix — unless you specifically need "
            "those tails intact for some downstream purpose (e.g. certain "
            "TSS/5′-end mapping protocols where the tail itself is "
            "informative), there's generally no reason to turn this off."
        )
    else:
        st.caption(
            "⚠️ Poly-G/poly-A tails will be left as fastp's default settings "
            "produce them, with no automatic detection or fix. Only turn "
            "this off if you have a specific reason to preserve these tails."
        )
    return {
        "reference": reference_cfg,
        "alignment_method": alignment_method,
        "threads": thread_cfg,
        "salmon_count_type": "NumReads",
        "use_tximport": _USE_TXIMPORT,
        "star_options": star_options_cfg,
        "auto_fix_poly_tails": auto_fix_poly_tails,
    }


# ---------------------------------------------------------------------------
# Status panel
# ---------------------------------------------------------------------------
def _render_status_panel(project):
    status = orch.get_status(project)
    if status.get("pipeline_status") == "not_started":
        return
    st.subheader("📡 Background Run Status")
    overall = status["pipeline_status"]
    if overall == "queued":
        # Written synchronously by launch_background_run() the instant
        # the background process is spawned -- see that function's
        # docstring. The subprocess itself hasn't necessarily started
        # doing real work yet (Python/pandas/etc. import overhead),
        # but this status is what makes the panel show up immediately
        # after pressing Launch, rather than needing a manual refresh.
        st.info("⏳ Queued — the background process has been launched and will begin shortly...")
    elif overall == "running":
        st.info(f"🔄 Running — currently on: {STAGE_LABELS.get(status.get('current_stage'), '—')}")
    elif overall == "complete":
        st.success("✅ Pipeline complete — gene counts matrix is ready.")
    elif overall == "error":
        st.error("❌ Pipeline stopped with an error — see details below.")
    for stage_key, _ in orch.PIPELINE_STAGES:
        entry = status["stages"][stage_key]
        icon = STAGE_ICONS.get(entry["status"], "•")
        st.markdown(f"{icon} **{STAGE_LABELS[stage_key]}** — {entry.get('message') or entry['status']}")
    if status.get("error"):
        with st.expander("Technical error details"):
            st.code(status["error"].get("traceback", ""))
    if orch.is_run_in_progress(project):
        if st.button("🔄 Refresh Status", key="adv_status_refresh_btn"):
            st.rerun()
    elif overall == "error":
        st.warning("You can fix the underlying issue and re-launch below — completed stages will be skipped (resumed), not redone.")


# ---------------------------------------------------------------------------
# Completed-project actions (shown in place of Steps 1-4 once done)
# ---------------------------------------------------------------------------
def _render_completed_project_actions(project):
    """
    Replaces the Steps 1-4 configuration wizard once this project's
    background run has actually finished (pipeline_status == "complete"),
    using the same shared project_actions.py actions Monitor Mode's
    launched-projects list offers: a packaged QC/MultiQC/alignment-
    scores/counts-table download (no raw or aligned reads -- see
    project_actions.py's module docstring), and a "Begin DE Analysis"
    hand-off straight to the Differential Expression workspace.

    Previously, the full Steps 1-4 wizard kept rendering unconditionally
    below the status panel even after "Pipeline complete" -- as if
    inviting reconfiguration/relaunch of an already-finished project,
    with no obvious "what do I do now that it's done" action. Steps 1-4
    are still reachable on demand via the "Reconfigure & Re-launch"
    expander below, for the less common case of wanting to re-run a
    completed project (e.g. with different options).
    """
    st.success(f"🎉 **{project}** is ready for downstream analysis.")
    if pa.project_package_available(project):
        col1, col2 = st.columns(2)
        with col1:
            zip_key = f"adv_pkg_zip_{project}"
            if st.button("📦 Prepare Download Package", key=f"{zip_key}_prep_btn"):
                st.session_state[zip_key] = pa.build_project_package_zip(project)
            if st.session_state.get(zip_key):
                st.download_button(
                    "⬇️ Download Project Package (.zip)",
                    data=st.session_state[zip_key],
                    file_name=f"{project}_package.zip",
                    mime="application/zip",
                    key=f"{zip_key}_dl_btn",
                    help="FastQC/MultiQC reports, alignment scores/info, and the gene counts matrix -- no raw or aligned read files.",
                )
        with col2:
            if st.button("🌋 Begin DE Analysis", key=f"adv_begin_de_{project}", type="primary"):
                pa.request_navigation_to_deseq2(project)
                st.rerun()
    else:
        st.caption("⏳ Gene counts matrix not found yet -- packaging/DE hand-off will appear once it's ready.")
    st.markdown("---")
    with st.expander("🔁 Reconfigure & Re-launch This Project"):
        st.caption(
            "This project already completed a run. Reconfiguring below "
            "and launching again will re-run the pipeline from Step 1 "
            "-- already-completed stages are resumed/skipped where "
            "possible (see advanced_mode_orchestrator's resumability), "
            "not blindly redone from scratch."
        )
        _render_configuration_wizard(project)


def _render_configuration_wizard(project):
    "Steps 1-4: FASTQ source -> metadata -> genome/options -> launch."
    fastq_config = _render_fastq_source_step(project)
    st.markdown("---")
    edited_metadata_df = _render_metadata_step(project, fastq_config or {})
    st.markdown("---")
    genome_and_options = _render_genome_and_options_step(project)
    st.markdown("---")
    st.subheader("Step 4: Launch")
    ready = fastq_config is not None and edited_metadata_df is not None and "sample" in edited_metadata_df.columns
    if not ready:
        st.info("Complete Steps 1–2 above (a valid FASTQ source and metadata with a `sample` column) to launch.")
        return
    if fastq_config.get("fastq_source") == "sra":
        st.warning(_ncbi_download_warning(len(fastq_config.get("sra_accessions", []))))
    if st.button("🚀 Save Configuration & Launch Background Run", key="adv_launch_btn", type="primary"):
        metadata_path = pm.metadata_path(project)
        edited_metadata_df.to_csv(metadata_path, index=False)
        config = {**fastq_config, "metadata_path": metadata_path, **genome_and_options}
        pid = orch.launch_background_run(project, config)
        st.success(
            f"✅ Background run launched (process ID {pid}). You can close "
            "this tab or navigate away — progress will be shown in the "
            "status panel above whenever you return to this page."
        )
        st.rerun()


# ---------------------------------------------------------------------------
# Bulk RNA-Seq's Auto pipeline (FASTQ -> ... -> gene counts matrix)
# ---------------------------------------------------------------------------
def _render_bulk_rnaseq_auto():
    """
    Everything that used to be this module's top-level render() before
    Auto became a pipeline-picker: the full Bulk RNA-Seq Auto wizard
    (FASTQ source -> metadata -> genome/options -> launch), unchanged
    in behavior -- only pulled into its own function so render() below
    can dispatch to it (or, in the future, to a Single-cell/Spatial
    equivalent) based on the pipeline dropdown.

    Once the project's run is complete, Steps 1-4 are replaced by
    _render_completed_project_actions() (download package + Begin DE
    Analysis) rather than continuing to show the configuration wizard
    -- see this module's docstring under "Post-completion behavior".
    """
    st.markdown(
        "Define your project's FASTQ source, metadata, reference genome, "
        "and pipeline options up front, then run FASTQ → trimming → "
        "alignment → gene counts matrix **non-interactively as a single "
        "background process**. DESeq2 contrasts and Ontology Analysis "
        "remain separate, interactive steps you run afterward, same as "
        "the step-by-step workflow."
    )
    _init_state()
    project = pm.render_project_selector(workspace_key="advanced_mode")
    if not project:
        return
    st.markdown("---")
    _render_status_panel(project)
    run_in_progress = orch.is_run_in_progress(project)
    if run_in_progress:
        st.markdown("---")
        st.info("A background run is currently in progress for this project. The configuration below is disabled until it finishes.")
        return
    status = orch.get_status(project)
    st.markdown("---")
    if status.get("pipeline_status") == "complete":
        _render_completed_project_actions(project)
        return
    _render_configuration_wizard(project)


# ---------------------------------------------------------------------------
# Pipeline dispatch -- Auto's own first screen
# ---------------------------------------------------------------------------
#
# Maps a pipeline's display label -> the function that renders its Auto
# wizard, or None for a pipeline that's on this app's roadmap (see this
# codebase's other planned-feature notes on Single-cell/Spatial RNA-seq)
# but doesn't have an Auto workflow built yet -- selecting one of those
# shows a plain "not built yet" message rather than omitting it from the
# dropdown entirely, so the picker's options consistently reflect every
# pipeline this portal supports (or will support), not just the ones
# with a finished Auto path today.
PIPELINE_AUTO_HANDLERS = {
    "🧬 Bulk RNA-Seq": _render_bulk_rnaseq_auto,
    "🧠 Single-cell RNA-Seq": None,
    "🧠 Spatial RNA-Seq": None,
}


def render():
    st.title("🤖 Auto")
    st.markdown(
        "Run any of this portal's pipelines **non-interactively as a "
        "single background process**, instead of walking through each "
        "step by hand -- pick which pipeline below, then define its "
        "FASTQ source, metadata, reference genome, and options up front "
        "on one guided screen."
    )
    pipeline_choice = st.selectbox(
        "Which pipeline would you like to run in Auto mode?",
        options=list(PIPELINE_AUTO_HANDLERS.keys()),
        key="adv_pipeline_choice",
    )
    handler = PIPELINE_AUTO_HANDLERS.get(pipeline_choice)
    st.markdown("---")
    if handler is None:
        st.info(
            f"⏳ An Auto workflow for **{pipeline_choice}** isn't built yet -- "
            "currently only Bulk RNA-Seq is supported here. Select "
            "\"🧬 Bulk RNA-Seq\" above, or use that pipeline's own "
            "step-by-step workspace in the sidebar in the meantime."
        )
        return
    handler()
