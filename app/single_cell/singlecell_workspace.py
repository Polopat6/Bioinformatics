"""
single_cell/singlecell_workspace.py

Streamlit UI for the Single-cell RNA-Seq pipeline -- Phase 1: Steps 1-6
(Ingestion + Chemistry -> Pre-trim QC -> Trimming -> Post-trim QC ->
Reference Setup -> STARsolo Alignment/Quantification/Cell-calling).

Scope: droplet-based UMI methods only -- 10x Genomics (all chemistry
generations), Drop-seq, inDrops, BD Rhapsody -- via STARsolo. Plate-based
(Smart-seq) and combinatorial-barcoding (Parse/split-seq) methods are NOT
supported; see chemistry_manager.py's module docstring for why.

--- Step 6 execution wiring (2026-08-17) ---
Previously, Step 6's "Run Alignment" button only showed an st.info()
placeholder describing what it WOULD do. This is now wired to real
execution, confirmed against the Bulk RNA-Seq pipeline's actual
quantification_manager.py and alignment_workspace.py source (pasted
directly rather than inferred -- both files' real interfaces were
independently confirmed before writing this):

- STARsolo uses the IDENTICAL genome index as ordinary STAR (confirmed
  directly from STAR's own official STARsolo.md: "The genome index is
  the same as for normal STAR runs") -- so this reuses
  quantification_manager.py's REAL qm.build_star_index()/
  qm.star_index_exists()/qm.detect_fastq_read_length() directly, exactly
  the same functions alignment_workspace.py's STAR path already calls,
  rather than reimplementing index-building logic in this pipeline.
- Preset species reuse the bulk pipeline's OWN shared, project-
  independent STAR index (project_manager.shared_star_index_dir(),
  confirmed real from alignment_workspace.py's exact same shared-vs-
  custom index selection logic) via reference_manager.ensure_shared_resource()
  -- a human genome's STAR index built by a bulk project is
  automatically reused here, and vice versa. Custom (non-preset)
  references get a per-project index instead (sc_project_manager.
  star_index_dir()), mirroring bulk's own project_manager.star_index_dir()
  convention for custom references.
- The actual alignment step calls starsolo_manager.run_starsolo_align(),
  which mirrors qm.run_star_align()'s exact confirmed subprocess pattern
  (capture_output/text/check=True, 2-hour timeout, returning
  (success, log, sample_output_dir)).
- Mapping-rate QC reuses quantification_manager.qm.parse_star_mapping_rate()
  and qm.classify_mapping_rate() UNCHANGED -- confirmed via direct testing
  (a real local subprocess run) that STARsolo writes Log.final.out in the
  exact same format as regular STAR, so bulk's existing parsing/
  classification logic works on single-cell's STARsolo output with zero
  modification. Cell-count reporting (STARsolo-specific, no bulk
  equivalent) uses starsolo_manager.parse_starsolo_summary()/
  get_cells_detected() instead.

--- Aligner-choice / genome-index ordering fix (2026-08-17, later same day) ---
The STAR genome-index build UI was moved into Step 5 earlier today (see
_render_genome_index_step's docstring) so it sits right next to genome
selection instead of being buried in Step 6. That move, on its own,
introduced a NEW ordering bug: Step 6's aligner choice (STARsolo vs.
alevin-fry) was still selected AFTER Step 5's index was already built --
so Step 5 would unconditionally build a STAR genomeGenerate index even
if the user was about to pick alevin-fry, which needs a completely
different artifact (a salmon index over a splici/spliceu reference, per
starsolo_manager.py's build_alevin_fry_commands()), not a STAR index.
Fixed by moving the aligner radio itself to the top of Step 5 (see
_render_aligner_choice()) and gating the genome-index build UI on that
choice -- STARsolo's STAR index only builds when STARsolo is selected;
alevin-fry shows its existing "not wired up yet" placeholder instead of
building the wrong index. Step 6 reads the confirmed choice back via
sc_project_manager.get_aligner_choice(project) rather than re-rendering
the radio a second time (which would raise a duplicate-widget-key error).

--- Aligner choice persistence (2026-08-17, same day) ---
The aligner choice is saved to the project's own project_info.json via
sc_project_manager.save_aligner_choice()/get_aligner_choice(), the same
"save unconditionally on every rerun" pattern already used for chemistry
(scpm.save_chemistry_choice, in _render_step1) -- so a previously
confirmed choice survives reopening the project after a restart instead
of always resetting to the STARsolo default, and Step 6 reads this same
persisted value directly rather than relying on st.session_state (which
would reset on a restart, and which _render_step5 rendering before
_render_step6 within the same run only accidentally makes work).

alevin-fry itself remains intentionally NOT fully wired up (no splici/
salmon-index build step, no real subprocess execution, no MTX/QC parsing
compatible with its output format) -- this was evaluated and deliberately
deferred as a larger follow-up piece of work, since STARsolo is already
the fully-functional recommended default and nothing downstream (Phase 2
cell-level QC) depends on alevin-fry existing.

Since single-cell projects use their own sc_project_manager.py rather
than the bulk pm module quantification_manager.py's functions were
originally written against, this workspace calls those qm functions
directly with single-cell-specific paths (sc_project_manager's
directory helpers) passed in as plain arguments -- qm's functions take
plain paths, not a project_manager module reference, so this works with
zero changes needed to quantification_manager.py itself.
"""
import os

import pandas as pd
import streamlit as st

import project_manager as pm  # reused directly for get_recommended_thread_count()
import reference_manager as ref
import file_browser as fb
import fastqc_manager as fastqc  # reused directly from the Bulk RNA-Seq pipeline
import fastp_manager as fastp_mgr  # reused directly for post-trim QC
import quantification_manager as qm  # reused directly for STAR index build + mapping-rate QC -- see module docstring

import sra_manager as sra  # reused directly for accession validation -- see sc_sra_manager.py

import sc_project_manager as scpm
import chemistry_manager as chem
import singlecell_ingestion_manager as ing
import singlecell_trim_manager as trim
import starsolo_manager as star
import sc_sra_manager as scsra

STAGE_LABELS = {
    "ingest": "1. FASTQ ingestion + chemistry + metadata",
    "pretrim_qc": "2. Pre-trim QC (FastQC + MultiQC)",
    "trimming": "3. Trimming (fastp, R2-only) + post-trim QC",
    "reference": "4. Reference genome + whitelist setup",
    "alignment": "5. STARsolo alignment + quantification + cell-calling",
}


def _render_project_selector():
    st.header("Step 0: Set Up Your Single-cell Project")
    st.markdown(
        "A single-cell **project** keeps your raw FASTQ files, chemistry "
        "choice, metadata, and every pipeline step's output together in "
        "one place -- separate from any Bulk RNA-Seq projects, even if "
        "they share the same name."
    )
    session_key = "sc_active_project"
    existing = scpm.list_projects()
    mode = st.radio(
        "Would you like to start a new project or continue an existing one?",
        ["🆕 Start a new project", "📂 Continue an existing project"],
        horizontal=True, index=0 if not existing else 1, key="sc_project_mode",
    )
    if mode == "🆕 Start a new project":
        new_name = st.text_input(
            "Name your project:", key="sc_new_project_name",
            help="Letters, numbers, dashes, or underscores -- spaces/special characters removed automatically.",
        )
        if new_name:
            safe_name = "".join(c for c in new_name.strip() if c.isalnum() or c in ("-", "_"))
            if st.button("✅ Create Project", key="sc_create_project_btn"):
                if not safe_name:
                    st.error("Please enter a valid project name.")
                elif safe_name in existing:
                    st.error(f"A project named `{safe_name}` already exists.")
                else:
                    scpm.create_project(safe_name)
                    st.session_state[session_key] = safe_name
                    st.rerun()
    else:
        if not existing:
            st.info("No existing single-cell projects found yet.")
        else:
            chosen = st.selectbox("Select a project:", options=existing, key="sc_existing_project_select")
            if st.button("📂 Open Project", key="sc_open_project_btn"):
                st.session_state[session_key] = chosen
                st.rerun()

    selected = st.session_state.get(session_key)
    if selected:
        info = scpm.load_info(selected)
        st.success(f"**Active single-cell project: `{selected}`**")
        completed = info.get("steps_completed", [])
        if completed:
            st.caption(f"✔️ Steps completed: {', '.join(completed)}")
        if st.button("🔄 Switch to a different project", key="sc_switch_project_btn"):
            del st.session_state[session_key]
            st.rerun()
    return selected


# ---------------------------------------------------------------------------
# Step 1 (NCBI/SRA source option)
# ---------------------------------------------------------------------------
def _render_sra_source(project, fastq_dir):
    st.markdown(
        "Enter a **study/BioProject accession** (e.g. `PRJNA474047`) "
        "and/or paste **individual run accessions** (e.g. `SRR1234567`)."
    )
    prefetch_ok, fasterq_ok = sra.tools_available()
    if not (prefetch_ok and fasterq_ok):
        st.error("⚠️ The SRA Toolkit (`prefetch`/`fasterq-dump`) isn't available on this system yet.")
        return

    st.info(
        "ℹ️ **Single-cell SRA data is messier than bulk.** Unlike bulk "
        "RNA-seq, the downloaded files' order does NOT reliably "
        "correspond to R1/R2 -- this pipeline classifies each file by "
        "read length after downloading and asks you to confirm the "
        "result below. Some accessions are deposited as 10x-specific "
        "BAM files rather than plain FASTQ, which this pipeline cannot "
        "yet reconstruct automatically -- you'll see a clear warning if "
        "that looks like what happened."
    )

    single_term = st.text_input("Study or BioProject accession:", placeholder="e.g. PRJNA474047", key="sc_sra_single_term")
    accession_list_text = st.text_area(
        "Or paste one or more individual run accessions:", placeholder="e.g.\nSRR1234567, SRR1234568",
        key="sc_sra_accession_list", height=90,
    )
    if st.button("🔍 Validate Accession(s) on NCBI", key="sc_sra_validate_btn"):
        combined = list(sra._split_accession_list_text(accession_list_text))
        if combined:
            with st.spinner(f"Validating {len(combined)} accession(s) on NCBI..."):
                success, rows, message, not_found = sra.lookup_multiple_accessions(combined)
            if success:
                st.session_state["sc_sra_lookup_rows"] = rows
                st.success(f"✅ {message}")
                if not_found:
                    st.warning(f"⚠️ Not found: {', '.join(not_found)}")
            else:
                st.error(f"⚠️ {message}")
                st.session_state["sc_sra_lookup_rows"] = None
        elif single_term:
            with st.spinner(f"Validating '{single_term}' on NCBI..."):
                success, rows, message = sra.lookup_accession(single_term)
            if success:
                st.session_state["sc_sra_lookup_rows"] = rows
                st.success(f"✅ {message}")
            else:
                st.error(f"⚠️ {message}")
                st.session_state["sc_sra_lookup_rows"] = None
        else:
            st.warning("⚠️ Please enter an accession or paste a list to validate.")

    lookup_rows = st.session_state.get("sc_sra_lookup_rows")
    if not lookup_rows:
        return

    assay_classifications = {row["Run"]: scsra.classify_assay_type_from_metadata(row) for row in lookup_rows}
    display_rows = [{
        "Run": row.get("Run", "—"), "Organism": row.get("ScientificName", "—"),
        "Library Strategy": row.get("LibraryStrategy", "unknown") if sra.is_rna_seq(row)
                            else f"⚠️ {row.get('LibraryStrategy', 'unknown')} (NOT RNA-Seq)",
        "Likely Assay Type (best-effort)": {
            scsra.ASSAY_SINGLE_CELL: "✅ Single-cell", scsra.ASSAY_SPATIAL: "⚠️ Spatial (unsupported here)",
            scsra.ASSAY_BULK: "⚠️ Looks like bulk", scsra.ASSAY_UNKNOWN: "❔ Unclear from metadata",
        }[assay_classifications[row["Run"]]["assay_type"]],
        "Size (MB)": row.get("size_MB", "—"),
    } for row in lookup_rows]
    st.dataframe(pd.DataFrame(display_rows), use_container_width=True, hide_index=True)
    st.caption(
        "ℹ️ \"Likely Assay Type\" is a best-effort guess from this accession's title/description "
        "text -- SRA's own metadata schema cannot distinguish bulk from single-cell/spatial data "
        "directly (Library Strategy is \"RNA-Seq\" for all three). Treat this as a helpful hint, "
        "not a certainty -- check the accession's own SRA/GEO page if you're unsure."
    )
    non_single_cell = [row["Run"] for row in lookup_rows if assay_classifications[row["Run"]]["assay_type"] in (scsra.ASSAY_SPATIAL, scsra.ASSAY_BULK)]
    if non_single_cell:
        st.warning(f"⚠️ {len(non_single_cell)} run(s) ({', '.join(non_single_cell)}) look like they might NOT be single-cell data based on their title/description -- you can still proceed if you believe this guess is wrong, but double-check first.")

    run_options = [row["Run"] for row in lookup_rows]
    selected_runs = st.multiselect("Select which validated run(s) to download:", options=run_options, default=run_options, key="sc_sra_selected_runs")
    if not selected_runs:
        return

    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    threads = st.slider("Threads for download:", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_sra_threads")

    if st.button("⬇️ Download & Classify Selected Run(s)", key="sc_sra_download_btn", type="primary"):
        results = {}
        progress = st.progress(0.0)
        for i, run in enumerate(selected_runs):
            with st.spinner(f"Downloading {run}..."):
                results[run] = scsra.download_and_classify_run(run, fastq_dir, run, threads=threads)
            progress.progress((i + 1) / len(selected_runs))
        st.session_state["sc_sra_download_results"] = results
        st.rerun()

    download_results = st.session_state.get("sc_sra_download_results")
    if not download_results:
        return

    st.markdown("---")
    st.markdown("**📋 Confirm File Roles Before Continuing**")
    st.caption("Each downloaded run's files were classified by read length. **Please confirm or correct these before continuing** -- an incorrect assignment here will break every downstream step.")
    any_finalized = False
    for run, result in download_results.items():
        with st.expander(f"Run: {run}", expanded=True):
            if not result["success"]:
                st.error(f"❌ {result['message']}")
                continue
            if result["bam_warning"]:
                st.warning(result["bam_warning"])
            overrides = {}
            for path, info in result["classification"].items():
                role_options = [scsra.ROLE_R1, scsra.ROLE_R2, scsra.ROLE_I1, scsra.ROLE_UNKNOWN]
                default_role = info["role"] if info["role"] in role_options else scsra.ROLE_UNKNOWN
                chosen_role = st.selectbox(
                    f"`{os.path.basename(path)}` ({info['length']}bp):", options=role_options,
                    index=role_options.index(default_role), key=f"sc_sra_role_{run}_{os.path.basename(path)}",
                )
                overrides[path] = chosen_role
            if st.button(f"✅ Confirm & Use These Files ({run})", key=f"sc_sra_confirm_btn_{run}"):
                destinations = scsra.finalize_role_assignment(overrides, fastq_dir, run)
                if scsra.ROLE_R1 in destinations and scsra.ROLE_R2 in destinations:
                    st.success(f"✅ {run}: R1/R2 files moved into project -- ready for chemistry detection below.")
                    any_finalized = True
                else:
                    st.error(f"⚠️ {run}: both an R1 and an R2 role must be assigned before this run can be used.")
    if any_finalized:
        st.rerun()


# ---------------------------------------------------------------------------
# Step 1: FASTQ ingestion + chemistry + metadata
# ---------------------------------------------------------------------------
def _render_step1(project):
    st.subheader("Step 1: FASTQ Files, Chemistry, and Sample Metadata")
    fastq_dir = scpm.fastq_dir(project)
    os.makedirs(fastq_dir, exist_ok=True)

    st.markdown("**Provide your FASTQ files:**")
    source = st.radio(
        "How will you provide FASTQ files?",
        ["📤 Upload from my computer", "📂 Browse a directory on this server", "🔎 Fetch from NCBI/SRA"],
        key="sc_fastq_source_radio", horizontal=True,
    )
    if source.startswith("📤"):
        uploaded = st.file_uploader(
            "Upload R1 + R2 FASTQ files (standard 10x naming: *_R1_001.fastq.gz / *_R2_001.fastq.gz):",
            type=["fastq", "gz", "fq"], accept_multiple_files=True, key="sc_fastq_upload",
        )
        if uploaded and st.button("💾 Save Uploaded Files", key="sc_fastq_upload_save_btn"):
            for f in uploaded:
                with open(os.path.join(fastq_dir, f.name), "wb") as out:
                    out.write(f.getbuffer())
            st.success(f"✅ Saved {len(uploaded)} file(s).")
            st.rerun()
        active_dir = fastq_dir
    elif source.startswith("📂"):
        active_dir = fb.render_server_directory_browser(
            key_prefix="sc_fastq_dir_browse", preview_extensions=[".fastq", ".fastq.gz", ".fq", ".fq.gz"],
            label="Browse for the directory containing your FASTQ files:",
        ) or fastq_dir
    else:
        _render_sra_source(project, fastq_dir)
        active_dir = fastq_dir

    pairs = ing.find_r1_r2_pairs(active_dir)
    unmatched = ing.find_unmatched_files(active_dir)
    if not pairs:
        st.info("No R1/R2 FASTQ pairs found yet in this location.")
        return None

    st.success(f"✅ Found {len(pairs)} sample(s): {', '.join(pairs.keys())}")
    warnings = ing.validate_pairs(pairs)
    for sample, msg in warnings.items():
        st.warning(f"⚠️ **{sample}**: {msg}")
    if unmatched:
        with st.expander(f"⚠️ {len(unmatched)} file(s) found but not recognized as standard R1/R2 FASTQ pairs"):
            for f in unmatched:
                st.caption(f)

    st.markdown("---")
    st.markdown("**🧪 Chemistry Detection**")
    st.caption(
        "This pipeline supports droplet-based UMI methods (10x Genomics all chemistries, "
        "Drop-seq, inDrops, BD Rhapsody) via STARsolo. Plate-based methods like Smart-seq (no "
        "cell barcodes at all) or combinatorial methods like Parse/split-seq are not supported."
    )
    first_sample = next(iter(pairs))
    first_r1 = pairs[first_sample]["r1_files"][0]
    detection = chem.detect_chemistry(first_r1)
    (st.success if detection["confidence"] == "high" else st.warning if detection["confidence"] == "ambiguous" else st.error)(detection["message"])

    all_chem_options = list(chem.CHEMISTRY_CATALOG.keys())
    default_idx = 0
    if detection["candidates"]:
        default_idx = all_chem_options.index(detection["candidates"][0])
    chosen_chemistry = st.selectbox(
        "Chemistry (confirm the auto-detected suggestion above, or pick manually):",
        options=all_chem_options, format_func=lambda k: chem.CHEMISTRY_CATALOG[k]["label"],
        index=default_idx, key="sc_chemistry_select",
    )
    was_auto = chosen_chemistry in detection["candidates"]
    if chem.CHEMISTRY_CATALOG[chosen_chemistry].get("whitelist_required"):
        if chem.whitelist_available(chosen_chemistry):
            checked, fraction, passed = chem.whitelist_confirm_match(first_r1, chosen_chemistry)
            if checked:
                icon = "✅" if passed else "⚠️"
                st.caption(f"{icon} Whitelist barcode match: {fraction*100:.1f}% of sampled R1 reads matched this chemistry's known barcode list.")
        else:
            st.caption("ℹ️ This chemistry's barcode whitelist file isn't installed on this server yet -- detection above is based on read length only. Contact your admin to install it for a stronger confirmation (see DEPLOYMENT.md).")
    scpm.save_chemistry_choice(project, chosen_chemistry, was_auto_detected=was_auto, user_confirmed=True)

    st.markdown("---")
    st.markdown("**📋 Sample Metadata** (one row per sample -- condition/treatment/donor, NOT per-cell data)")
    metadata_path = scpm.metadata_path(project)
    sra_lookup_rows = st.session_state.get("sc_sra_lookup_rows")
    sra_selected_runs = st.session_state.get("sc_sra_selected_runs")
    if os.path.exists(metadata_path):
        default_df = pd.read_csv(metadata_path)
    elif sra_lookup_rows and sra_selected_runs:
        default_df = pd.DataFrame(sra.build_metadata_dataframe(sra_lookup_rows, selected_runs=sra_selected_runs))
        if not sra.has_any_descriptive_metadata(sra_lookup_rows, selected_runs=sra_selected_runs):
            st.caption("ℹ️ NCBI provided no characteristics beyond the run accession for these runs -- add your own columns below.")
    else:
        default_df = pd.DataFrame({"sample": list(pairs.keys())})
    edited_df = st.data_editor(default_df, num_rows="dynamic", use_container_width=True, key="sc_metadata_editor")
    if "sample" not in edited_df.columns:
        st.error("⚠️ Metadata must have a column named exactly `sample` matching your FASTQ sample names.")
        return None
    if st.button("💾 Save Metadata", key="sc_save_metadata_btn"):
        edited_df.to_csv(metadata_path, index=False)
        st.success("✅ Metadata saved.")

    if os.path.exists(metadata_path):
        scpm.mark_step_complete(project, "ingest")
        scpm.save_fastq_source_dir(project, active_dir)
        return pairs
    return None


# ---------------------------------------------------------------------------
# Step 2: Pre-trim QC
# ---------------------------------------------------------------------------
def _render_pretrim_qc_step(project, pairs):
    st.subheader("Step 2: Pre-trim QC (FastQC + MultiQC)")
    st.caption(
        "Run on R2 (the biological cDNA read). R1 (cell barcode + UMI) "
        "can look 'abnormal' in these reports -- e.g. failing per-base "
        "sequence content checks -- this is EXPECTED and not a problem, "
        "since R1 is a structured barcode, not random biological sequence."
    )

    fastqc_ok, multiqc_ok = fastqc.tools_available()
    if not (fastqc_ok and multiqc_ok):
        st.error("⚠️ FastQC/MultiQC aren't available on this system yet -- see DEPLOYMENT.md / the ⚙️ Setup & Deployment page.")
        return

    r2_files = [path for entry in pairs.values() for path in entry["r2_files"]]
    fastqc_dir = scpm.fastqc_dir(project)
    multiqc_dir = scpm.multiqc_dir(project)

    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    threads = st.slider("Threads for FastQC/MultiQC (pre-trim):", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_pretrim_qc_threads")

    if st.button("▶️ Run Pre-trim QC", key="sc_run_pretrim_qc_btn", type="primary"):
        os.makedirs(fastqc_dir, exist_ok=True)
        os.makedirs(multiqc_dir, exist_ok=True)
        with st.spinner(f"Running FastQC on {len(r2_files)} R2 file(s)..."):
            fastqc_success, fastqc_log = fastqc.run_fastqc(r2_files, fastqc_dir, threads=threads)
        if not fastqc_success:
            st.error("❌ FastQC failed for pre-trim QC.")
            with st.expander("📜 FastQC log"):
                st.code(fastqc_log or "(no output captured)", language="text")
        else:
            with st.spinner("Running MultiQC..."):
                multiqc_success, multiqc_log = fastqc.run_multiqc(fastqc_dir, multiqc_dir)
            if not multiqc_success:
                st.error("❌ MultiQC failed.")
                with st.expander("📜 MultiQC log"):
                    st.code(multiqc_log or "(no output captured)", language="text")
            else:
                st.session_state["sc_pretrim_qc_ran"] = True
                scpm.mark_step_complete(project, "pretrim_qc")
                st.success(f"✅ Pre-trim QC complete -- {len(r2_files)} R2 file(s) checked.")

    if st.session_state.get("sc_pretrim_qc_ran") or scpm.has_completed_step(project, "pretrim_qc"):
        summary_df = fastqc.parse_fastqc_summaries(fastqc_dir)
        if summary_df is not None and not summary_df.empty:
            overview_df, details_by_file = fastqc.build_quality_flags(summary_df)
            st.markdown("**Pre-trim QC Overview:**")
            st.dataframe(overview_df, use_container_width=True, hide_index=True)
            for filename, explanations in details_by_file.items():
                if explanations:
                    with st.expander(f"ℹ️ Details for `{filename}`"):
                        for explanation in explanations:
                            st.markdown(explanation)
        multiqc_report = os.path.join(multiqc_dir, "multiqc_report.html")
        if os.path.isfile(multiqc_report):
            st.caption(f"📄 Full combined report: `{multiqc_report}`")

        st.markdown("---")
        action_col1, action_col2 = st.columns(2)
        with action_col1:
            if os.path.isfile(multiqc_report):
                with open(multiqc_report, "rb") as f:
                    st.download_button(
                        "⬇️ Download MultiQC Report (.html)",
                        data=f.read(), file_name=f"{project}_pretrim_multiqc_report.html",
                        mime="text/html", key="sc_pretrim_qc_multiqc_dl_btn",
                    )
        with action_col2:
            if st.button("➡️ Proceed to Trimming", key="sc_proceed_to_trimming_btn", type="primary"):
                st.session_state["nav_request"] = "🧪 SC Trimming & Post-Trim QC"
                st.rerun()


# ---------------------------------------------------------------------------
# Step 4: Post-trim QC
# ---------------------------------------------------------------------------
def _render_posttrim_qc_step(project, pairs):
    st.subheader("Step 4: Post-trim QC (fastp reports + MultiQC)")
    st.caption(
        "Reuses fastp's own per-sample JSON report -- generated automatically as a "
        "side effect of Step 3's trimming -- rather than re-running FastQC on the "
        "trimmed reads. MultiQC has a built-in parser for fastp's JSON output, so "
        "this gives the same 'one combined report' result without redundant "
        "recomputation, exactly matching the Bulk RNA-Seq pipeline's own design."
    )

    _fastp_ok, multiqc_ok = fastp_mgr.tools_available()
    if not multiqc_ok:
        st.error("⚠️ MultiQC isn't available on this system yet -- see DEPLOYMENT.md / the ⚙️ Setup & Deployment page.")
        return

    fastp_reports_dir = scpm.fastp_reports_dir(project)
    posttrim_multiqc_dir = scpm.posttrim_multiqc_dir(project)

    if not os.path.isdir(fastp_reports_dir) or not any(f.endswith(".fastp.json") for f in os.listdir(fastp_reports_dir)):
        st.info("Run Step 3 (Trimming) first -- no fastp reports found yet.")
        return

    if st.button("▶️ Run Post-trim QC", key="sc_run_posttrim_qc_btn", type="primary"):
        os.makedirs(posttrim_multiqc_dir, exist_ok=True)
        with st.spinner("Running MultiQC on fastp's reports..."):
            multiqc_success, multiqc_log = fastp_mgr.run_multiqc(fastp_reports_dir, posttrim_multiqc_dir)
        if not multiqc_success:
            st.error("❌ MultiQC failed.")
            with st.expander("📜 MultiQC log"):
                st.code(multiqc_log or "(no output captured)", language="text")
        else:
            st.session_state["sc_posttrim_qc_ran"] = True
            scpm.mark_step_complete(project, "posttrim_qc")
            st.success("✅ Post-trim QC complete.")

    if st.session_state.get("sc_posttrim_qc_ran") or scpm.has_completed_step(project, "posttrim_qc"):
        summary_df = fastp_mgr.parse_fastp_reports(fastp_reports_dir)
        if summary_df is not None and not summary_df.empty:
            display_df = summary_df.drop(columns=["_raw_adapter_cutting_json"], errors="ignore")
            st.markdown("**Post-trim Summary (per sample, from fastp's own report):**")
            st.dataframe(display_df, use_container_width=True, hide_index=True)

        poly_tail_issues = fastp_mgr.scan_all_samples_for_poly_tail_issues(fastp_reports_dir)
        if poly_tail_issues:
            st.warning(f"⚠️ {len(poly_tail_issues)} sample(s) show a residual poly-G/poly-A tail signature even after trimming:")
            for sample_name, issues in poly_tail_issues.items():
                with st.expander(f"⚠️ `{sample_name}` -- {len(issues)} issue(s) flagged"):
                    for issue in issues:
                        before_note = f" (was {issue['tail_pct_before']}% before trimming)" if issue.get("tail_pct_before") is not None else ""
                        st.markdown(
                            f"- **{issue['read']}**, poly-**{issue['base']}**: tail region averages "
                            f"{issue['tail_pct_after']}% {issue['base']} after trimming (vs. "
                            f"{issue['body_pct_after']}% in the rest of the read, ~{issue['ratio']}x higher)"
                            f"{before_note}."
                        )
                    if sample_name in pairs:
                        if st.button(f"🔁 Re-trim `{sample_name}` with poly-tail removal enabled", key=f"sc_retrim_polytail_{sample_name}"):
                            r1, r2 = pairs[sample_name]["r1_files"][0], pairs[sample_name]["r2_files"][0]
                            trimmed_dir = scpm.trimmed_fastq_dir(project)
                            paths = trim.fastp_output_paths(trimmed_dir, fastp_reports_dir, sample_name)
                            quality = st.session_state.get("sc_trim_quality", trim.DEFAULT_QUALIFIED_QUALITY_PHRED)
                            length_required = st.session_state.get("sc_trim_length", trim.DEFAULT_LENGTH_REQUIRED)
                            re_trim_threads = st.session_state.get("sc_trim_threads", 4)
                            with st.spinner(f"Re-trimming {sample_name} with poly-tail removal enabled..."):
                                success, message = trim.run_trim_sample(
                                    r1, r2, paths, qualified_quality_phred=quality,
                                    length_required=length_required, threads=re_trim_threads,
                                    auto_fix_poly_tails=True,
                                )
                            (st.success if success else st.error)(f"{sample_name}: {message}")
                            if success:
                                st.session_state["sc_posttrim_qc_ran"] = False
                                st.rerun()

        posttrim_multiqc_report = os.path.join(posttrim_multiqc_dir, "multiqc_report.html")

        st.markdown("---")
        action_col1, action_col2 = st.columns(2)
        with action_col1:
            if os.path.isfile(posttrim_multiqc_report):
                with open(posttrim_multiqc_report, "rb") as f:
                    st.download_button(
                        "⬇️ Download MultiQC Report (.html)",
                        data=f.read(), file_name=f"{project}_posttrim_multiqc_report.html",
                        mime="text/html", key="sc_posttrim_qc_multiqc_dl_btn",
                    )
        with action_col2:
            if st.button("➡️ Proceed to Alignment & Cell-Calling", key="sc_proceed_to_alignment_btn", type="primary"):
                st.session_state["nav_request"] = "🧬 SC Alignment & Cell-Calling"
                st.rerun()


# ---------------------------------------------------------------------------
# Step 3: Trimming
# ---------------------------------------------------------------------------
def _render_step3(project, pairs):
    st.subheader("Step 3: Trimming (fastp -- R2 only)")
    st.markdown(
        "**R1 (cell barcode + UMI) is never trimmed** -- only R2 (the "
        "cDNA read) gets adapter/quality trimming. This is essential: "
        "any length change to R1 would corrupt barcode/UMI extraction "
        "for every downstream step."
    )
    st.caption(
        "ℹ️ Because R1 is protected by never being passed to fastp at all, R2 is trimmed alone in "
        "single-end mode. This means fastp's adapter detection uses its single-end method rather "
        "than its more robust paired-end overlap-based detection -- a real, unavoidable trade-off "
        "for keeping R1 completely intact."
    )

    quality = st.slider(
        "Minimum quality score (Phred) for R2:", min_value=5, max_value=30,
        value=trim.DEFAULT_QUALIFIED_QUALITY_PHRED, key="sc_trim_quality",
        help=("Bases below this quality are considered low-quality. Single-cell data commonly "
              "uses a gentler threshold (Q15) than bulk RNA-seq's usual Q20/Q30."),
    )

    chemistry_info = scpm.get_chemistry_choice(project)
    chemistry_key = chemistry_info["chemistry_key"] if chemistry_info else None
    chemistry_spec = chem.CHEMISTRY_CATALOG.get(chemistry_key) if chemistry_key else None

    with st.expander("ℹ️ What does \"minimum R2 length\" mean, and does chemistry matter?", expanded=False):
        st.markdown(
            "**This slider is a FLOOR, not a target.** Any R2 read shorter than this "
            "value *after* quality/adapter trimming is discarded outright.\n\n"
            "**Does chemistry affect this floor?** Not directly -- there's no documented "
            "chemistry-specific value for the floor itself. **What chemistry DOES affect "
            "is the recommended target sequencing length**, confirmed from 10x's own "
            "documentation: **91bp** for 3' v3/v3.1 and 5' v2, **98bp** for 5' v1."
        )
        if chemistry_spec and chemistry_spec.get("recommended_r2_length"):
            confirmed = chemistry_spec.get("recommended_r2_length_confirmed")
            icon = "✅" if confirmed else "ℹ️"
            confidence_note = "confirmed directly from 10x's own documentation for this exact chemistry" if confirmed else "a reasonable estimate, not independently confirmed for this one"
            st.markdown(f"{icon} **Your selected chemistry** (`{chemistry_spec['label']}`) has a recommended target R2 length of **~{chemistry_spec['recommended_r2_length']}bp** ({confidence_note}).")
        elif chemistry_key:
            st.caption("ℹ️ No specific target R2 length is documented for your selected chemistry.")

    length_required = st.slider(
        "Minimum R2 read length after trimming (floor -- see explanation above):",
        min_value=10, max_value=50, value=trim.DEFAULT_LENGTH_REQUIRED, key="sc_trim_length",
    )

    st.markdown("---")
    auto_fix_poly_tails = st.toggle(
        "Automatically detect and trim poly-A/poly-G tails on R2",
        value=trim.DEFAULT_AUTO_FIX_POLY_TAILS, key="sc_trim_auto_fix_poly_tails",
    )
    if auto_fix_poly_tails:
        st.caption(
            "✅ **Recommended (on by default).** This matters *more* for single-cell 3' data "
            "specifically: a peer-reviewed 2022 study (Svoboda et al., *NAR Genomics and "
            "Bioinformatics*) confirms internal oligo(dT) priming causes systematic poly-A "
            "read-through contamination in single-cell RNA-seq, explicitly naming 10x."
        )
    else:
        st.caption("⚠️ Poly-A/poly-G tails will be left as fastp's default settings produce them.")

    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    threads = st.slider("Threads for fastp:", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_trim_threads")
    if st.button("▶️ Run Trimming", key="sc_run_trim_btn"):
        trimmed_dir = scpm.trimmed_fastq_dir(project)
        fastp_dir = scpm.fastp_reports_dir(project)
        os.makedirs(trimmed_dir, exist_ok=True)
        os.makedirs(fastp_dir, exist_ok=True)
        progress = st.progress(0.0)
        for i, (sample, entry) in enumerate(pairs.items()):
            r1, r2 = entry["r1_files"][0], entry["r2_files"][0]
            paths = trim.fastp_output_paths(trimmed_dir, fastp_dir, sample)
            success, message = trim.run_trim_sample(
                r1, r2, paths, qualified_quality_phred=quality, length_required=length_required,
                threads=threads, auto_fix_poly_tails=auto_fix_poly_tails,
            )
            (st.success if success else st.error)(f"{sample}: {message}")
            progress.progress((i + 1) / len(pairs))
        scpm.mark_step_complete(project, "trimming")


# ---------------------------------------------------------------------------
# Step 5: Reference + whitelist setup
# ---------------------------------------------------------------------------
_ANNOTATION_GFF_EXTENSIONS = (".gff", ".gff3", ".gff.gz", ".gff3.gz")
_ANNOTATION_ALL_EXTENSIONS = (".gtf", ".gtf.gz") + _ANNOTATION_GFF_EXTENSIONS
_FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")

# Aligner choice keys -- shared between Step 5 (where the choice is made
# and used to gate genome-index building) and Step 6 (where the choice
# determines which quantification path runs). See the module docstring's
# "Aligner-choice / genome-index ordering fix" section for why this
# choice now lives in Step 5 rather than Step 6.
ALIGNER_STARSOLO = "starsolo"
ALIGNER_ALEVIN_FRY = "alevin_fry"


def _detect_annotation_format(path_or_filename):
    'Return "gff3" if path_or_filename ends in a GFF/GFF3 extension, else "gtf".'
    lower = (path_or_filename or "").lower()
    return "gff3" if lower.endswith(_ANNOTATION_GFF_EXTENSIONS) else "gtf"


def _render_gff3_warning_if_needed(reference_cfg):
    if reference_cfg.get("annotation_format") != "gff3":
        return
    st.warning(
        "⚠️ **GFF/GFF3 detected**, not GTF. STAR needs an extra flag "
        "(`--sjdbGTFtagExonParentTranscript Parent`) to build a genome index "
        "correctly from GFF3. The safer path is converting to GTF first:\n\n"
        "```\ngffread your_file.gff3 -T -o your_file.gtf\n```\n\n"
        "(`gffread` is already available in this portal's environment)."
    )


def _render_reference_file_validation(genome_fasta, gtf_path):
    if genome_fasta and os.path.isfile(genome_fasta):
        if ref.validate_fasta_file(genome_fasta):
            st.caption(f"✅ `{os.path.basename(genome_fasta)}` looks like a valid FASTA file.")
        else:
            st.error(f"⚠️ `{os.path.basename(genome_fasta)}` does not look like a valid FASTA file.")
    if gtf_path and os.path.isfile(gtf_path):
        if ref.validate_annotation_file(gtf_path):
            st.caption(f"✅ `{os.path.basename(gtf_path)}` looks like a valid GTF/GFF annotation file.")
        else:
            st.error(f"⚠️ `{os.path.basename(gtf_path)}` does not look like a valid GTF/GFF annotation file.")


_ALIGNER_RADIO_OPTIONS = ["STARsolo (recommended)", "alevin-fry (faster, alignment-free)"]
_ALIGNER_LABEL_TO_KEY = {
    _ALIGNER_RADIO_OPTIONS[0]: ALIGNER_STARSOLO,
    _ALIGNER_RADIO_OPTIONS[1]: ALIGNER_ALEVIN_FRY,
}


def _render_aligner_choice(project):
    """
    Alignment/quantification method choice (STARsolo vs. alevin-fry) --
    rendered at the TOP of Step 5, immediately before genome/reference
    selection (moved here 2026-08-17, same day as the genome-index build
    UI's own move into Step 5).

    WHY THIS MOVED: the genome-index build UI was moved into Step 5
    earlier today so it sits next to genome selection instead of being
    buried in Step 6 -- but that move alone introduced a NEW ordering
    bug, since Step 6's aligner radio was still selected AFTER Step 5's
    index was already built. STARsolo and alevin-fry need fundamentally
    different index artifacts (a STAR genomeGenerate index vs. a salmon
    index built from a splici/spliceu reference, per
    starsolo_manager.py's build_alevin_fry_commands()) -- so building
    "the" genome index before knowing which aligner was chosen risked
    silently building the wrong one. Moving the choice here, ABOVE the
    index-build UI, and gating that UI on this choice, fixes this.

    PERSISTENCE (added same day): the choice is now saved to this
    project's own project_info.json via sc_project_manager.py's
    save_aligner_choice(), the same "save unconditionally on every
    rerun" pattern _render_step1 already uses for chemistry
    (scpm.save_chemistry_choice) -- so a previously-made choice survives
    reopening the project after a restart, rather than always resetting
    to the STARsolo default. The radio's own default index is likewise
    read back from scpm.get_aligner_choice() rather than hardcoded, so a
    returning project shows its last confirmed choice pre-selected.

    Returns one of ALIGNER_STARSOLO / ALIGNER_ALEVIN_FRY. Step 6 reads
    this same value back out via scpm.get_aligner_choice(project) --
    the persisted source of truth -- rather than out of
    st.session_state, so Step 6 is correct even if it's ever reached
    without Step 5 re-rendering first in the same run (e.g. a future
    direct-navigation entry point), not just within the current
    Streamlit session.

    Also warns (mirroring the Bulk RNA-Seq pipeline's own Salmon<->STAR
    switch warning in alignment_workspace.py) if this choice differs
    from the previously persisted one -- switching aligners means Step 5
    (reference/index setup) and Step 6 (alignment/quantification/cell-
    calling) need to be redone for the newly selected method; previous
    results from the other method are left on disk but won't be reused.
    """
    st.markdown("**🧬 Alignment/Quantification Method**")

    saved_aligner = scpm.get_aligner_choice(project)
    default_index = 1 if saved_aligner == ALIGNER_ALEVIN_FRY else 0

    aligner_choice = st.radio(
        "Alignment/quantification method:",
        _ALIGNER_RADIO_OPTIONS,
        index=default_index,
        key="sc_aligner_choice", horizontal=True,
        help=("STARsolo (alignment-based) is the recommended default -- 2025/2026 benchmarking found it "
              "preserves more accurate cell identities and marker expression than alignment-free "
              "alternatives like alevin-fry. Chosen here, before reference/index setup below, since "
              "the two methods need different index types."),
    )
    aligner = _ALIGNER_LABEL_TO_KEY[aligner_choice]

    if saved_aligner and saved_aligner != aligner:
        saved_label = "STARsolo" if saved_aligner == ALIGNER_STARSOLO else "alevin-fry"
        st.warning(
            f"⚠️ This project previously used **{saved_label}**. "
            "Switching methods means Step 5 (reference/index setup) and Step 6 "
            "(alignment/quantification/cell-calling) will need to be run again for "
            "the new method — previous results from the other method are kept "
            "on disk but won't be used going forward."
        )

    scpm.save_aligner_choice(project, aligner)

    if aligner == ALIGNER_ALEVIN_FRY:
        st.info(
            "ℹ️ alevin-fry execution isn't wired up yet in this pipeline (no splici-reference/"
            "salmon-index build step, no alignment execution). The genome/index setup below is "
            "STARsolo-specific and will be skipped -- switch back to STARsolo to continue "
            "through Step 5/6, or come back once alevin-fry support is added."
        )
    return aligner


def _render_step5(project, pairs):
    st.subheader("Step 5: Reference Genome & Whitelist Setup")

    aligner = _render_aligner_choice(project)
    st.markdown("---")

    st.caption(
        "STARsolo uses the exact same genome index as ordinary STAR -- "
        "provide a genome FASTA + gene annotation (GTF or GFF/GFF3) from "
        "a preset (auto-downloaded), a location already on this server/"
        "HPC, or by uploading directly from your computer."
    )

    existing_reference = scpm.get_reference_choice(project)
    if existing_reference:
        if existing_reference.get("is_custom"):
            fmt = existing_reference.get("annotation_format", "gtf").upper()
            st.caption(f"ℹ️ Currently confirmed: custom reference (`{os.path.basename(existing_reference.get('custom_genome_fasta') or '—')}` + {fmt} annotation).")
        else:
            species_key = existing_reference.get("species_key")
            label = ref.REFERENCE_CATALOG.get(species_key, {}).get("label", species_key)
            st.caption(f"ℹ️ Currently confirmed: preset reference -- {label}.")

    source = st.radio(
        "Where should the reference genome come from?",
        ["📚 Use a preset genome (auto-downloaded)", "📂 Browse a directory on this server (HPC)", "📤 Upload from my computer"],
        key="sc_ref_source_radio", horizontal=True,
    )
    reference_cfg = {}

    if source.startswith("📚"):
        species_options = list(ref.REFERENCE_CATALOG.keys())
        species_labels = {k: v["label"] for k, v in ref.REFERENCE_CATALOG.items()}
        species_choice = st.selectbox("Reference organism:", options=species_options, format_func=lambda k: species_labels[k], key="sc_species_choice")
        reference_cfg.update({"is_custom": False, "species_key": species_choice})

        target_dir = pm.shared_genome_dir(species_choice)
        already_ready = os.path.isdir(target_dir) and len(os.listdir(target_dir)) > 0
        genome_fasta_expected = os.path.join(target_dir, f"{species_choice}.genome.fa")
        gtf_expected = os.path.join(target_dir, f"{species_choice}.annotation.gtf")

        if already_ready and os.path.isfile(genome_fasta_expected) and os.path.isfile(gtf_expected):
            st.success(f"✅ Reference files already prepared (shared, reused across projects): `{os.path.basename(genome_fasta_expected)}` + `{os.path.basename(gtf_expected)}`")
            reference_cfg.update({"custom_genome_fasta": genome_fasta_expected, "custom_gtf": gtf_expected, "annotation_format": "gtf"})
        else:
            st.info("ℹ️ This reference hasn't been downloaded on this server yet. Downloading can take several minutes -- this only needs to happen ONCE per species.")
            if st.button("⬇️ Download & Prepare Reference", key="sc_ref_download_btn", type="primary"):
                def _build_fn(temp_dir, _species=species_choice):
                    success, _paths, message = ref.download_genome_and_gtf(_species, temp_dir)
                    return success, message

                wait_placeholder = st.empty()

                def _wait_message(elapsed_seconds, _ph=wait_placeholder):
                    _ph.info(f"⏳ Another project is currently preparing this same reference -- waiting ({elapsed_seconds:.0f}s so far)...")

                with st.spinner(f"Preparing reference for {species_labels[species_choice]} (this may take several minutes)..."):
                    success, message, built_by_this_call = ref.ensure_shared_resource(
                        target_dir, build_fn=_build_fn, wait_message_callback=_wait_message,
                    )
                wait_placeholder.empty()
                if success:
                    verb = "Downloaded" if built_by_this_call else "Confirmed"
                    st.success(f"✅ {verb}: {message}")
                    st.rerun()
                else:
                    st.error(f"❌ Reference preparation failed: {message}")
            reference_cfg.update({"custom_genome_fasta": None, "custom_gtf": None})

    elif source.startswith("📂"):
        genome_fasta = fb.render_server_file_browser(
            key_prefix="sc_ref_genome_browse_hpc", file_extensions=list(_FASTA_EXTENSIONS),
            label="Browse for the genome FASTA already on this server/HPC:",
        )
        gtf_path = fb.render_server_file_browser(
            key_prefix="sc_ref_gtf_browse_hpc", file_extensions=list(_ANNOTATION_ALL_EXTENSIONS),
            label="Browse for the gene annotation (GTF or GFF/GFF3) already on this server/HPC:",
        )
        reference_cfg.update({"is_custom": True, "species_key": None, "custom_genome_fasta": genome_fasta, "custom_gtf": gtf_path})
        if gtf_path:
            reference_cfg["annotation_format"] = _detect_annotation_format(gtf_path)
        _render_reference_file_validation(genome_fasta, gtf_path)

    else:  # "📤 Upload from my computer"
        custom_dir = scpm.custom_reference_dir(project)
        os.makedirs(custom_dir, exist_ok=True)
        uploaded_fasta = st.file_uploader("Upload genome FASTA:", type=["fa", "fasta", "fna", "gz"], key="sc_ref_fasta_upload")
        uploaded_annotation = st.file_uploader("Upload gene annotation (GTF or GFF/GFF3):", type=["gtf", "gff", "gff3", "gz"], key="sc_ref_annotation_upload")
        if (uploaded_fasta or uploaded_annotation) and st.button("💾 Save Uploaded Reference File(s)", key="sc_ref_upload_save_btn"):
            if uploaded_fasta:
                with open(os.path.join(custom_dir, uploaded_fasta.name), "wb") as out:
                    out.write(uploaded_fasta.getbuffer())
            if uploaded_annotation:
                with open(os.path.join(custom_dir, uploaded_annotation.name), "wb") as out:
                    out.write(uploaded_annotation.getbuffer())
            st.success("✅ Saved uploaded reference file(s).")
            st.rerun()

        existing_files = os.listdir(custom_dir) if os.path.isdir(custom_dir) else []
        existing_fasta = next((f for f in existing_files if f.lower().endswith(_FASTA_EXTENSIONS)), None)
        existing_annotation = next((f for f in existing_files if f.lower().endswith(_ANNOTATION_ALL_EXTENSIONS)), None)
        if existing_fasta:
            st.success(f"✅ Genome FASTA on file: `{existing_fasta}`")
        if existing_annotation:
            st.success(f"✅ Gene annotation on file: `{existing_annotation}`")
        genome_fasta = os.path.join(custom_dir, existing_fasta) if existing_fasta else None
        gtf_path = os.path.join(custom_dir, existing_annotation) if existing_annotation else None
        reference_cfg.update({"is_custom": True, "species_key": None, "custom_genome_fasta": genome_fasta, "custom_gtf": gtf_path})
        if gtf_path:
            reference_cfg["annotation_format"] = _detect_annotation_format(gtf_path)
        _render_reference_file_validation(genome_fasta, gtf_path)

    _render_gff3_warning_if_needed(reference_cfg)

    # --- Genome index status/build -- STARsolo-specific, gated on aligner choice ---
    # Moved into Step 5 on 2026-08-17 so it sits right after genome
    # selection instead of being buried in Step 6. Later the same day,
    # gated on `aligner == ALIGNER_STARSOLO` -- see _render_aligner_choice's
    # docstring and the module docstring's "Aligner-choice / genome-index
    # ordering fix" section for why: STARsolo and alevin-fry need
    # different index artifacts, so this STAR-specific build UI must not
    # run unconditionally before the aligner is even chosen.
    st.markdown("---")
    if aligner == ALIGNER_STARSOLO:
        _render_genome_index_step(project, reference_cfg, pairs)
    else:
        st.markdown("**🔧 alevin-fry Reference & Index**")
        st.info(
            "ℹ️ alevin-fry needs a different reference artifact than STARsolo -- a salmon index "
            "built from a splici/spliceu transcriptome (not a STAR genome index) -- and that "
            "build step isn't implemented in this pipeline yet. Switch back to **STARsolo** "
            "above to continue setting up a reference and index."
        )

    chemistry_info = scpm.get_chemistry_choice(project)
    if chemistry_info:
        chemistry_key = chemistry_info["chemistry_key"]
        spec = chem.CHEMISTRY_CATALOG[chemistry_key]
        st.markdown(f"**Whitelist for `{spec['label']}`:**")
        if spec.get("whitelist_required"):
            if chem.whitelist_available(chemistry_key):
                st.success(f"✅ Whitelist file found: `{spec['whitelist_file']}`")
            else:
                st.error(f"❌ Whitelist file `{spec['whitelist_file']}` not found on this server. Contact your admin to install it, or alignment cannot proceed for this chemistry.")
        else:
            st.info(f"ℹ️ {spec['label']} does not use a fixed whitelist -- barcodes are called computationally.")
    else:
        st.warning("⚠️ Complete Step 1 (chemistry selection) before this step.")

    if st.button("✅ Confirm Reference & Whitelist Setup", key="sc_confirm_reference_btn"):
        missing = []
        if reference_cfg.get("is_custom"):
            if not reference_cfg.get("custom_genome_fasta"):
                missing.append("genome FASTA")
            if not reference_cfg.get("custom_gtf"):
                missing.append("gene annotation (GTF/GFF)")
        if missing:
            st.error(f"⚠️ Missing: {', '.join(missing)} -- provide these above before confirming.")
        else:
            scpm.save_reference_choice(project, reference_cfg)
            scpm.mark_step_complete(project, "reference")
            st.success("✅ Reference setup confirmed.")


# ---------------------------------------------------------------------------
# Step 6: STARsolo alignment -- REAL EXECUTION (2026-08-17)
# ---------------------------------------------------------------------------
def _resolve_star_index_dir(project, reference_cfg):
    """
    Resolve (WITHOUT building) the STAR genome index location this
    project's alignment will use -- reusing quantification_manager.py's
    REAL qm.star_index_exists() directly (confirmed against that
    module's actual source), since STARsolo's genome index is identical
    to ordinary STAR's (confirmed from STAR's own official STARsolo.md).
    Preset species reuse the bulk pipeline's OWN shared index location
    (project_manager.shared_star_index_dir()) -- mirroring
    alignment_workspace.py's exact same shared-vs-custom index selection
    logic (confirmed against that real module's source) -- so a human
    STAR index already built by a bulk project is reused here
    automatically, and vice versa. Custom references get a per-project
    index instead (sc_project_manager.star_index_dir()).

    Returns (index_dir: str, ready: bool, index_is_shared: bool).
    """
    is_custom = reference_cfg.get("is_custom", False)
    species_key = reference_cfg.get("species_key")
    index_is_shared = not is_custom and bool(species_key)
    index_dir = pm.shared_star_index_dir(species_key) if index_is_shared else scpm.star_index_dir(project)
    ready = qm.star_index_exists(index_dir)
    return index_dir, ready, index_is_shared


def _render_genome_index_step(project, reference_cfg, pairs):
    """
    Genome index status + build UI -- shown in Step 5 immediately after
    genome selection (moved here 2026-08-17; previously lived buried in
    Step 6, well below the aligner/cell-calling choices, disconnected
    from the genome selection it actually depends on).

    STARsolo-specific: only called from _render_step5 when the aligner
    choice (rendered above, at the top of Step 5 -- see
    _render_aligner_choice) is ALIGNER_STARSOLO. See the module
    docstring's "Aligner-choice / genome-index ordering fix" section for
    why this gating was added.

    Requires reference_cfg to already have a real genome_fasta + gtf_path
    on disk (from Step 5's genome selection above this call) -- if
    either is missing, shows nothing at all rather than a confusing
    "index status" for a reference that isn't even selected yet.
    """
    genome_fasta = reference_cfg.get("custom_genome_fasta")
    gtf_path = reference_cfg.get("custom_gtf")
    if not genome_fasta or not gtf_path or not os.path.isfile(genome_fasta) or not os.path.isfile(gtf_path):
        return

    st.markdown("**🔧 STAR/STARsolo Genome Index**")
    st.caption(
        "STARsolo uses the identical genome index as ordinary STAR -- this is a one-time build "
        "per reference (shared across every project using this species, if using a preset)."
    )

    index_dir, index_ready, index_is_shared = _resolve_star_index_dir(project, reference_cfg)
    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    index_threads = st.slider("Threads for indexing:", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_star_index_threads")

    # --- Detect R2 read length for STAR's --sjdbOverhang -- confirmed via
    # quantification_manager.qm.detect_fastq_read_length(), same function
    # alignment_workspace.py already uses for bulk. Detected from R2 (the
    # biological cDNA read), NOT R1 (fixed-length barcode+UMI structure,
    # which would give a meaningless "read length" for this purpose).
    detected_read_length = None
    for entry in pairs.values():
        if entry["r2_files"]:
            detected_read_length = qm.detect_fastq_read_length(entry["r2_files"][0], is_gzipped=entry["r2_files"][0].endswith(".gz"))
            if detected_read_length:
                break
    sjdb_overhang = (detected_read_length - 1) if detected_read_length else 100
    if detected_read_length:
        st.caption(f"📏 Detected R2 read length: **{detected_read_length}bp** -- using `--sjdbOverhang {sjdb_overhang}` for the genome index.")
    else:
        st.caption("📏 Could not detect R2 read length automatically -- falling back to the default `--sjdbOverhang 100`.")

    if index_ready:
        st.success("✅ STAR genome index already built " + ("(shared across every project using this species)." if index_is_shared else "for this project's reference."))
    else:
        st.info("ℹ️ **No genome index found yet -- build it below before alignment can run.** This is a one-time step per reference and can take a while for larger genomes.")

    build_clicked = False
    if not index_ready:
        build_clicked = st.button("🔧 Build Genome Index", key="sc_build_star_index_btn", type="primary")
    else:
        with st.expander("🔁 Force a fresh re-index (advanced)"):
            st.caption("⚠️ Only do this if you have a specific reason to believe the current index is corrupted.")
            if st.button("🔄 Re-build Index", key="sc_build_star_index_force_btn"):
                build_clicked = True

    if build_clicked:
        def _build_index_impl(temp_dir_or_final):
            return qm.build_star_index(genome_fasta, gtf_path, temp_dir_or_final, threads=index_threads, sjdb_overhang=sjdb_overhang)

        if index_is_shared:
            status_placeholder = st.empty()

            def _wait_cb(elapsed):
                status_placeholder.info(f"⏳ Another project is currently building this index. Waiting... ({int(elapsed)}s elapsed)")

            with st.spinner("Building genome index... this may take a while for larger genomes."):
                success, message, built = ref.ensure_shared_resource(index_dir, _build_index_impl, wait_message_callback=_wait_cb)
            status_placeholder.empty()
        else:
            with st.spinner("Building genome index... this may take a while for larger genomes."):
                success, message = qm.build_star_index(genome_fasta, gtf_path, index_dir, threads=index_threads, sjdb_overhang=sjdb_overhang)

        if success:
            st.success("✅ Genome index built successfully.")
        else:
            st.error("Genome indexing failed. Details below:")
            st.code(message)


def _render_mapping_rate_qc(results_df, rate_column):
    """
    Render the 3-tier plain-language mapping-rate QC summary (poor/
    caution/good) for a results_df that has a numeric rate_column (e.g.
    "Uniquely Mapped (%)") -- factored out of _render_step6's run-block
    2026-08-17 so the exact same tiering logic can be reused both right
    after a fresh run AND when re-displaying previously PERSISTED
    results on reopening a project (see _render_step6's restructuring
    for why persisted results need to be shown independent of a fresh
    run). Mirrors the Bulk RNA-Seq pipeline's own
    alignment_workspace._render_mapping_rate_qc pattern exactly.
    """
    poor = [(row["Sample"], qm.classify_mapping_rate(row[rate_column])) for _, row in results_df.iterrows() if qm.classify_mapping_rate(row[rate_column])["tier"] == "poor"]
    caution = [(row["Sample"], qm.classify_mapping_rate(row[rate_column])) for _, row in results_df.iterrows() if qm.classify_mapping_rate(row[rate_column])["tier"] == "caution"]
    if poor:
        lines = "\n".join(f"- **{s}**: {c['message']}" for s, c in poor)
        st.error(f"🔴 **{len(poor)} sample(s) have a very low mapping rate (below 50%)** -- a strong signal something doesn't match (wrong reference organism/strain most commonly):\n\n{lines}")
    if caution:
        lines = "\n".join(f"- **{s}**: {c['message']}" for s, c in caution)
        st.warning(f"🟡 **{len(caution)} sample(s) have a lower-than-ideal mapping rate (50-79%)**:\n\n{lines}")
    if not poor and not caution and len(results_df) > 0:
        st.success("🟢 All samples have a healthy mapping rate (80%+).")


def _render_starsolo_run_controls(project, pairs, cell_filter):
    """
    STAR/STARsolo tool + reference + index availability checks, and the
    actual "Run/Re-run Alignment" button + execution -- factored out of
    _render_step6 2026-08-17 as part of a bug fix (see _render_step6's
    docstring/comments): previously these live-availability checks ran
    BEFORE checking whether alignment had already completed, so if
    STAR/the genome index/reference files were ever momentarily
    unreachable (e.g. right after an environment/session restart, before
    a scratch mount or conda env was back in place), Step 6 would bail
    out early and never even show that alignment had already succeeded
    -- even though that completion flag was safely persisted the whole
    time. Now this function is ONLY invoked when the user is actually
    about to run (or re-run) alignment, not merely to view prior status.

    On a successful run, persists the per-sample results to disk via
    scpm.save_alignment_results() (new 2026-08-17) in addition to the
    existing boolean scpm.mark_step_complete() flag -- previously only
    the boolean was saved, so the actual metrics table (Cells Detected,
    Uniquely Mapped %, Quality) had no way to reappear after a reload
    without re-running alignment from scratch.
    """
    if not star.star_tool_available():
        st.error("⚠️ STAR was not found on this system. See DEPLOYMENT.md / the ⚙️ Setup & Deployment page.")
        return

    reference_cfg = scpm.get_reference_choice(project)
    if not reference_cfg:
        st.warning("⚠️ Complete Step 5 (reference setup) before this step.")
        return
    chemistry_info = scpm.get_chemistry_choice(project)
    if not chemistry_info:
        st.warning("⚠️ Complete Step 1 (chemistry selection) before this step.")
        return
    chemistry_key = chemistry_info["chemistry_key"]
    chemistry_spec = chem.CHEMISTRY_CATALOG[chemistry_key]

    genome_fasta = reference_cfg.get("custom_genome_fasta")
    gtf_path = reference_cfg.get("custom_gtf")
    if not genome_fasta or not gtf_path or not os.path.isfile(genome_fasta) or not os.path.isfile(gtf_path):
        st.error("⚠️ Could not find a genome + annotation for this project. Please complete Step 5 first.")
        return

    whitelist_path = chem.shared_whitelist_path(chemistry_key)
    if chemistry_spec.get("whitelist_required") and not (whitelist_path and os.path.isfile(whitelist_path)):
        st.error(f"❌ Whitelist file for `{chemistry_spec['label']}` isn't installed on this server yet. Revisit Step 5 for details.")
        return

    detected_cores, recommended_threads = pm.get_recommended_thread_count()

    # --- Genome index readiness check ONLY -- the actual build UI (with
    # its own thread slider) now lives in Step 5, immediately after
    # genome selection and gated on the STARsolo aligner choice (see
    # _render_genome_index_step, moved there 2026-08-17 so it's not
    # buried below the aligner/cell-calling choices, disconnected from
    # the genome selection it depends on).
    index_dir, index_ready, index_is_shared = _resolve_star_index_dir(project, reference_cfg)
    if not index_ready:
        st.warning("⚠️ No genome index has been built yet for this reference. Revisit **Step 5** (right after genome selection) to build it before alignment can run.")
        return

    st.markdown("---")

    # --- Alignment execution ---
    st.markdown("**🚀 Run STARsolo Alignment**")
    align_dir = scpm.starsolo_output_dir(project)
    align_done = scpm.has_completed_step(project, "alignment")

    align_threads = st.slider("Threads per sample (alignment):", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_star_align_threads")

    align_label = "🔄 Re-run Alignment" if align_done else "🚀 Align & Call Cells for All Samples"

    if st.button(align_label, key="sc_run_alignment_btn", type="primary"):
        progress_bar = st.progress(0, text="Starting alignment...")
        results = []
        sample_names = list(pairs.keys())

        for i, sample_name in enumerate(sample_names):
            entry = pairs[sample_name]
            progress_bar.progress(i / len(sample_names), text=f"Aligning {i + 1}/{len(sample_names)}: {sample_name}...")
            success, log, sample_out_dir = star.run_starsolo_align(
                sample_name, index_dir, entry["r1_files"], entry["r2_files"], whitelist_path,
                cb_len=chemistry_spec["cb_len"], umi_len=chemistry_spec["umi_len"],
                output_base_dir=align_dir, cell_filter=cell_filter, threads=align_threads,
                strand=chemistry_spec["strand"],
            )
            # Reuses the REAL bulk qm.parse_star_mapping_rate() UNCHANGED --
            # confirmed via direct testing that STARsolo writes Log.final.out
            # in the exact same format regular STAR does.
            mapping_rate = qm.parse_star_mapping_rate(sample_out_dir, sample_name) if success else None
            quality = qm.classify_mapping_rate(mapping_rate) if success else None
            summary_dict = star.parse_starsolo_summary(sample_out_dir, sample_name) if success else {}
            cells_detected = star.get_cells_detected(summary_dict) if summary_dict else None
            results.append({
                "Sample": sample_name,
                "Status": "✅ Success" if success else "❌ Failed",
                "Cells Detected": cells_detected if cells_detected is not None else "—",
                "Uniquely Mapped %": f"{mapping_rate}%" if mapping_rate is not None else "—",
                "Uniquely Mapped (%)": mapping_rate,
                "Quality": f"{quality['icon']} {quality['short_label']}" if quality else "—",
            })
            if not success:
                with st.expander(f"Error details for {sample_name}"):
                    st.code(log)

        progress_bar.progress(1.0, text="Alignment complete.")
        results_df = pd.DataFrame(results)
        st.dataframe(results_df.drop(columns=["Uniquely Mapped (%)"]), use_container_width=True, hide_index=True)

        if (results_df["Status"] == "✅ Success").any():
            scpm.mark_step_complete(project, "alignment")
            # Persist the actual per-sample metrics, not just the boolean
            # completion flag -- see this function's docstring for why.
            scpm.save_alignment_results(project, results)
            st.success("✅ STARsolo alignment & cell-calling complete for at least one sample.")

        _render_mapping_rate_qc(results_df, rate_column="Uniquely Mapped (%)")


def _render_step6(project, pairs):
    st.subheader("Step 6: Alignment, Quantification & Cell-calling (STARsolo)")

    # Aligner choice is now made at the top of Step 5 (see
    # _render_aligner_choice) rather than rendered here -- re-rendering
    # the same radio here (with the same key) would raise a duplicate-
    # widget-key error, so we just read the confirmed choice back
    # instead. Read from the PERSISTED project setting
    # (scpm.get_aligner_choice()) rather than st.session_state, so this
    # is correct even across a restart or a fresh session where Step 5
    # hasn't rendered yet in the current run -- not just within a single
    # continuous Streamlit session. Defaults to STARsolo if never set
    # (e.g. a project created before this setting existed).
    aligner = scpm.get_aligner_choice(project) or ALIGNER_STARSOLO
    aligner_label = "STARsolo (recommended)" if aligner == ALIGNER_STARSOLO else "alevin-fry (faster, alignment-free)"
    st.caption(f"Using alignment/quantification method chosen in Step 5: **{aligner_label}**.")

    if aligner == ALIGNER_ALEVIN_FRY:
        st.info("ℹ️ alevin-fry execution isn't wired up yet in this pipeline -- STARsolo is fully functional below. Switch back to STARsolo in Step 5 to continue.")
        return

    cell_filter = st.radio(
        "Cell-calling method:", list(star.CELL_FILTER_OPTIONS.keys()),
        format_func=lambda k: star.CELL_FILTER_OPTIONS[k]["label"],
        index=list(star.CELL_FILTER_OPTIONS.keys()).index(star.DEFAULT_CELL_FILTER), key="sc_cell_filter_choice",
    )
    st.caption(star.CELL_FILTER_OPTIONS[cell_filter]["explanation"])

    st.markdown("---")

    # --- Completion status + persisted results -- checked and shown
    # FIRST, before any live STAR/reference/index availability checks.
    # THIS IS THE FIX for a real reported bug (2026-08-17): previously,
    # _render_step6 ran star.star_tool_available()/genome+GTF file
    # existence/whitelist/genome-index-readiness checks BEFORE ever
    # looking at scpm.has_completed_step(project, "alignment") -- so if
    # ANY of those live checks failed (e.g. right after an environment
    # reload, before STAR was back on PATH or a scratch-mounted
    # reference/index was remounted), the function returned early with a
    # warning/error and never reached the "already completed" check at
    # all, even though that completion flag -- and now, the actual
    # per-sample results too -- were safely persisted in
    # project_info.json/alignment_results.json the whole time. Viewing
    # PAST results no longer requires STAR, the index, or reference
    # files to currently be reachable; those are only required to
    # actually (re-)run alignment, which now lives in the "Re-run"
    # expander below via _render_starsolo_run_controls().
    align_done = scpm.has_completed_step(project, "alignment")
    persisted_results = scpm.get_alignment_results(project) if align_done else None

    if align_done:
        st.success("✅ Alignment & cell-calling has already been run for this project.")
        if persisted_results:
            results_df = pd.DataFrame(persisted_results)
            st.dataframe(results_df.drop(columns=["Uniquely Mapped (%)"], errors="ignore"), use_container_width=True, hide_index=True)
            if "Uniquely Mapped (%)" in results_df.columns:
                _render_mapping_rate_qc(results_df, rate_column="Uniquely Mapped (%)")
        else:
            st.caption(
                "ℹ️ Detailed per-sample metrics from that run weren't saved for this project "
                "(it may have completed before this was tracked) -- re-run below to regenerate them."
            )

        st.markdown("---")
        # nav_request target "🔬 SC Cell-level QC" -- confirmed 2026-08-17
        # against the real app.py: this requires a matching entry in
        # PIPELINE_GROUPS["single_cell"]["options"] (so the nav_request ->
        # radio-key lookup in app.py actually finds it) AND a matching
        # elif branch in app.py's routing chain calling render_cell_qc()
        # below -- both added alongside this button. Without both, this
        # button would set active_workspace to a value matching no elif
        # branch, silently rendering a blank page.
        if st.button("➡️ Proceed to Phase 2: Cell-level QC", key="sc_proceed_to_cellqc_btn", type="primary"):
            st.session_state["nav_request"] = "🔬 SC Cell-level QC"
            st.rerun()

        st.markdown("---")
        with st.expander("🔁 Re-run alignment (advanced)"):
            st.caption(
                "Only do this if you have a specific reason -- e.g. you changed chemistry, "
                "reference, or STAR options since the last run. This will overwrite the results "
                "shown above once it completes."
            )
            _render_starsolo_run_controls(project, pairs, cell_filter)
    else:
        _render_starsolo_run_controls(project, pairs, cell_filter)


def _require_project_with_source_dir():
    project = _render_project_selector()
    if not project:
        return None, None
    st.markdown("---")
    source_dir = scpm.get_fastq_source_dir(project)
    if not source_dir:
        st.warning("⚠️ Complete **Step 1** on the **🧫 Single-cell RNA-Seq** page first (FASTQ ingestion + chemistry + metadata) before continuing here.")
        return project, None
    pairs = ing.find_r1_r2_pairs(source_dir)
    if not pairs:
        st.error(f"❌ No R1/R2 FASTQ pairs found anymore in this project's previously-configured source directory (`{source_dir}`) -- it may have been moved, renamed, or emptied since Step 1 ran. Please revisit the **🧫 Single-cell RNA-Seq** page to reselect your FASTQ source.")
        return project, None
    for sample, msg in ing.validate_pairs(pairs).items():
        st.warning(f"⚠️ **{sample}**: {msg}")
    return project, pairs


def render_ingestion():
    st.title("🧫 Single-cell RNA-Seq")
    st.markdown(
        "Process droplet-based single-cell RNA-seq data: FASTQ ingestion "
        "+ chemistry detection, then pre-trim QC. **Supports 10x Genomics "
        "(all chemistries), Drop-seq, inDrops, and BD Rhapsody** -- "
        "plate-based (Smart-seq) and combinatorial-barcoding (Parse/"
        "split-seq) methods are not yet supported."
    )
    st.markdown("---")
    project = _render_project_selector()
    if not project:
        return
    st.markdown("---")
    pairs = _render_step1(project)
    st.markdown("---")
    if not pairs:
        st.info("Complete Step 1 above to continue to pre-trim QC.")
        return
    _render_pretrim_qc_step(project, pairs)


def render_trimming():
    st.title("🧪 Single-cell Trimming & Post-Trim QC")
    st.markdown("Trim adapter/quality issues from R2 only (R1's cell barcode + UMI is never modified), then re-run QC on the trimmed reads.")
    st.markdown("---")
    project, pairs = _require_project_with_source_dir()
    if not pairs:
        return
    st.markdown("---")
    _render_step3(project, pairs)
    st.markdown("---")
    _render_posttrim_qc_step(project, pairs)


def render_alignment():
    st.title("🧬 Single-cell Alignment & Cell-Calling")
    st.markdown("Set up the reference genome + barcode whitelist, then run STARsolo (or, optionally, alevin-fry) to align, quantify, and call real cells from ambient/empty droplets.")
    st.markdown("---")
    project, pairs = _require_project_with_source_dir()
    if not pairs:
        return
    st.markdown("---")
    _render_step5(project, pairs)
    st.markdown("---")
    _render_step6(project, pairs)


def render_cell_qc():
    """
    Phase 2 entry point -- Cell-level QC (scDblFinder doublet detection,
    SoupX/DecontX ambient RNA correction, interactive per-cell
    filtering). Added 2026-08-17 as a real (not dead-end) navigation
    target for Step 6's "➡️ Proceed to Phase 2: Cell-level QC" button --
    intentionally a placeholder for now, pending the actual Phase 2
    manager module(s) and UI (not yet built). Requires alignment
    (Step 6) to have completed for the active project, mirroring every
    other render_*() entry point's own step-gating pattern in this file.
    """
    st.title("🔬 Single-cell Cell-level QC")
    st.markdown(
        "Doublet detection (scDblFinder), ambient RNA correction "
        "(SoupX/DecontX), and interactive per-cell filtering -- applied "
        "to STARsolo's filtered cell x gene matrix from Step 6."
    )
    st.markdown("---")
    project, pairs = _require_project_with_source_dir()
    if not pairs:
        return
    if not scpm.has_completed_step(project, "alignment"):
        st.warning("⚠️ Complete **Step 6** (🧬 SC Alignment & Cell-Calling) first -- this page needs a completed STARsolo run's filtered cell x gene matrix to work from.")
        return
    st.info(
        "🚧 **Coming soon.** Phase 2 (cell-level QC) hasn't been built yet -- this page is a "
        "placeholder so navigation from Step 6 works correctly in the meantime. Your completed "
        "alignment results are safe and waiting here for when this is implemented."
    )
