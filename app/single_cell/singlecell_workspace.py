"""
single_cell/singlecell_workspace.py

Streamlit UI for the Single-cell RNA-Seq pipeline -- Phase 1: Steps 1-6
(Ingestion + Chemistry -> Pre-trim QC -> Trimming -> Post-trim QC ->
Reference Setup -> STARsolo Alignment/Quantification/Cell-calling), plus
Phase 2: Cell-level QC (doublet detection, ambient RNA correction,
per-cell filtering with standard visualizations and a downloadable
results package).

--- Preset-reference mitochondrial verification + redownload (2026-08-17) ---
A real reported gap: Step 5's preset ("auto") reference section only
ever checked "does this directory already have files?" to decide
whether to show "✅ already prepared" -- never what's actually IN those
files. Combined with the mito-detection bug fixed the same day (see
sc_cellqc_manager.py's own module docstring), this made it genuinely
impossible to tell whether an already-downloaded preset reference
actually includes the mitochondrial genome, without either
re-downloading blindly or manually inspecting files on disk.

Fixed by calling reference_manager.verify_preset_reference_mito_content()
immediately after an "already prepared" preset reference is found,
surfacing the result prominently (not buried in an unrelated "advanced"
expander), and offering a DEDICATED "🔄 Re-download to include
mitochondrial genome" button specifically when verification fails --
distinct from (and more discoverable than) the pre-existing generic
"force re-download" advanced option, which existed but gave no
indication of WHY someone might want to use it.

--- Custom-reference mitochondrial gene resolution UX (2026-08-17,
    later same day) ---
For a CUSTOM (user-uploaded) reference, Step 5 now runs the same direct
GTF-based mito auto-detection immediately after GTF confirmation. If
that finds zero genes, three explicit resolution paths are offered
(rather than silently reporting a misleading 0%):
  1. "Not applicable for this organism" -- explicit opt-out.
  2. "Try matching against a preset organism's known mitochondrial gene
     symbols" -- uses sc_cellqc_manager.get_mito_gene_symbols_from_gtf()
     on a chosen PRESET species' own already-downloaded GTF, then
     match_custom_genes_by_mito_symbol() to search the custom GTF for
     matching gene_name values. Only offered for preset species whose
     reference has actually already been downloaded on this system
     (checked via reference_manager). Always shown with an explicit
     lower-confidence disclaimer, since this assumes shared gene-naming
     convention between the custom reference and the chosen preset --
     not a verified identity.
  3. "Specify manually" -- free-text gene ID/symbol entry, matched
     directly against the custom GTF.
Whichever path is used, the resolved gene ID list (and a "mito_source"
label describing which path was used) is persisted into this project's
saved reference_choice dict, so Phase 2's Cell-level QC page doesn't
need to re-resolve it on every run.

--- Downloadable QC results package (2026-08-17, later same day) ---
render_cell_qc() now offers a "📦 Download QC Package (.zip)" button
after a completed run, bundling the per-cell metrics CSV, QC thresholds,
mitochondrial-gene diagnostic, top-ambient-genes table, DoubletFinder's
pK sweep (if used), a self-contained markdown summary (via
sc_cellqc_manager.build_qc_summary_markdown()), and PNG renders of
every QC visualization (via Plotly's kaleido-based image export,
already a project dependency) into a single zip file -- so a completed
run's full results can be shared or archived without this app open.

--- Genome-index force-rebuild bug fix + stale-index diagnostic
    (2026-08-17, later same day) ---
A real reported bug, discovered after a user re-indexed a preset human
reference (following the mitochondrial-verification fix above), saw the
"reindex" complete in ~15 seconds (versus the genuine 10-15 minutes a
real STAR human-genome index build takes), re-ran Step 6 alignment, and
STILL saw zero mitochondrial genes in Cell-level QC afterward.

Root cause: _render_genome_index_step's "🔄 Re-build Index" (force
re-index) button set the SAME build_clicked flag used for "build an
index that doesn't exist yet", and BOTH cases fell through to the
identical ref.ensure_shared_resource(index_dir, _build_index_impl,
wait_message_callback=_wait_cb) call -- with NO force=True ever passed.
ensure_shared_resource()'s own documented behavior is to return
immediately with "already available" WITHOUT ever invoking build_fn if
force is not explicitly set AND the resource directory already has
files in it -- which is exactly the case for a force-reindex click,
since the whole point of that button is that an index already exists.
So "Re-build Index" silently did nothing beyond a fast existence check,
every time it was clicked, for every preset/shared reference in this
app -- explaining both the suspiciously-fast "reindex" and why
mitochondrial genes never appeared afterward (the STALE index, built
before the reference's mitochondrial content was confirmed/fixed, was
never actually replaced, so STARsolo's own features.tsv -- generated
once, AT INDEX-BUILD TIME, and never revisited at alignment time --
still had no mitochondrial gene rows in it, regardless of what the
CURRENT reference GTF says).

Fixed by tracking build_clicked and force_rebuild as TWO SEPARATE
flags and passing force=force_rebuild explicitly to
ref.ensure_shared_resource() -- this now exactly mirrors the Bulk
RNA-Seq pipeline's own real, original alignment_workspace.py pattern
(confirmed against that module's actual source), which already tracked
build_clicked/force_index separately and passed force=force_index
correctly. The single-cell version had collapsed this into one flag
during an earlier reconstruction, introducing this regression.

Also added: sc_cellqc_manager.diagnose_starsolo_matrix_for_mito(), a
fast, R-independent, PRE-FLIGHT check run automatically the moment a
sample is selected in render_cell_qc() -- BEFORE the (much slower)
full R-based Cell-level QC pipeline runs. It reads STARsolo's own
features.tsv directly to check whether THIS SAMPLE's already-aligned
count matrix has any mitochondrial gene rows at all, independent of
what the reference GTF currently says. This distinguishes, immediately
and without needing a full QC re-run:
  1. The matrix already has mitochondrial genes -> proceed normally.
  2. The matrix has genes overall but ZERO are mitochondrial -> a
     STALE INDEX/alignment (the exact bug above) -- surfaced with a
     specific, actionable message pointing at Step 5's real rebuild
     button and Step 6's re-run-alignment action, NOT a suggestion to
     just re-run Cell-level QC again (which cannot fix this, since QC
     only re-reads existing STARsolo output and never re-aligns).

--- Stale-resource-during-rebuild race fix (2026-08-18) ---
A real reported bug: a user force-rebuilding a shared STAR genome index
navigated between radio-button pages while the rebuild was running (on
an HPC host, confirmed still running via `top`) and, on navigating
back, saw the index reported as "✅ already built" even though the real
rebuild was still in progress. Root cause: reference_manager.py's
ensure_shared_resource() builds a fresh copy into a private temp
directory and only replaces the OLD index via a single atomic
os.rename() at the very END of a successful build -- so the old
(stale) index stays fully present on disk, and any existence-only
check (qm.star_index_exists) reports "ready" for the ENTIRE rebuild
duration, with no way to distinguish a finished build from one that's
actively being replaced.

Fixed by using reference_manager.resource_build_in_progress() (a
non-blocking probe on the resource's existing .lock file) as part of
the authoritative readiness check in _resolve_star_index_dir(), and by
showing an explicit "🔄 build in progress, please wait" status in
_render_genome_index_step() that blocks further build actions while
true. Also, CUSTOM (non-shared, per-project) index builds now route
through ref.ensure_shared_resource() too -- previously that branch
called qm.build_star_index() DIRECTLY against index_dir with NO
temp-dir staging, atomic replacement, or lock-file protection at all,
so a custom reference's rebuild had the identical stale-appears-ready
symptom with even less of a safety net than the shared path. Both
branches now share the same temp-directory + atomic-rename + lock-file
infrastructure -- there is no longer a separate "shared" vs. "custom"
code path for the actual build call, only for which directory is
targeted. Because _render_starsolo_run_controls() in Step 6 also calls
_resolve_star_index_dir() for its own "is the index ready?" gate, this
fix automatically blocks alignment too while a rebuild is in progress,
with no separate Step 6 change needed.

See sc_project_manager.py and sc_cellqc_manager.py's own module
docstrings for full detail on all other Phase 1/Phase 2 design
decisions and fixes made earlier the same day.
"""
import io
import os
import zipfile

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import project_manager as pm
import reference_manager as ref
import file_browser as fb
import fastqc_manager as fastqc
import fastp_manager as fastp_mgr
import quantification_manager as qm

import sra_manager as sra

import sc_project_manager as scpm
import chemistry_manager as chem
import singlecell_ingestion_manager as ing
import singlecell_trim_manager as trim
import starsolo_manager as star
import sc_sra_manager as scsra
import sc_cellqc_manager as cellqc

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


def _render_preset_mito_verification(species_choice, species_labels, genome_fasta_expected, gtf_expected):
    """
    Verify + surface whether a PRESET reference's already-downloaded
    files actually include the mitochondrial genome -- see this module's
    own docstring, "Preset-reference mitochondrial verification +
    redownload", for the full rationale.

    Rendered immediately after Step 5's "✅ Reference files already
    prepared" message for a preset species -- NOT buried in the
    unrelated generic "force re-download" advanced expander below it,
    since this is a concrete, actionable, and fairly likely-to-matter
    result (as opposed to the generic re-download option, which exists
    for a vaguer "I suspect corruption" scenario).
    """
    verification = ref.verify_preset_reference_mito_content(genome_fasta_expected, gtf_expected)

    if verification["verified"]:
        st.success(
            f"✅ Mitochondrial genome verified in this reference: the genome FASTA includes a "
            f"mitochondrial contig, and {verification['gtf_mito_gene_count']} mitochondrial "
            f"gene(s) were found in the annotation."
        )
        return

    st.error(
        "🔴 **Could not verify mitochondrial genome content in this already-downloaded "
        "reference.** "
        f"Genome FASTA mitochondrial contig found: {'✅ Yes' if verification['fasta_has_mito_contig'] else '❌ No'}. "
        f"Mitochondrial genes found in annotation: {verification['gtf_mito_gene_count']}.\n\n"
        "If this reference was downloaded before this check existed, or from a partial/"
        "interrupted earlier download, it may be missing the mitochondrial chromosome entirely -- "
        "this would make single-cell mitochondrial QC metrics meaningless (a misleading 0% for "
        "every cell), and could also affect alignment accuracy for any real mitochondrial-"
        "origin reads in both this and the Bulk RNA-Seq pipeline, since they share this same "
        "downloaded reference."
    )
    if st.button(
        f"🔄 Re-download {species_labels[species_choice]} Reference (to include mitochondrial genome)",
        key="sc_ref_mito_redownload_btn", type="primary",
    ):
        target_dir = pm.shared_genome_dir(species_choice)

        def _build_fn(temp_dir, _species=species_choice):
            success, _paths, message = ref.download_genome_and_gtf(_species, temp_dir)
            return success, message

        wait_placeholder = st.empty()

        def _wait_message(elapsed_seconds, _ph=wait_placeholder):
            _ph.info(f"⏳ Another project is currently preparing this same reference -- waiting ({elapsed_seconds:.0f}s so far)...")

        with st.spinner(f"Re-downloading reference for {species_labels[species_choice]}..."):
            success, message, _built = ref.ensure_shared_resource(
                target_dir, build_fn=_build_fn, wait_message_callback=_wait_message, force=True,
            )
        wait_placeholder.empty()
        if success:
            st.success(f"✅ {message}")
            st.warning(
                "⚠️ **Important:** re-downloading the reference does NOT rebuild the STAR genome "
                "index or re-run alignment automatically -- both use a separately-cached, one-time "
                "build. Scroll down to **rebuild the genome index** (using the real force-rebuild fix "
                "below), then revisit **Step 6** to re-run alignment for every sample, before this "
                "fix will actually show up in Cell-level QC results."
            )
            st.rerun()
        else:
            st.error(f"❌ Re-download failed: {message}")


_MITO_RESOLUTION_OPTIONS = [
    "Not applicable for this organism",
    "Try matching against a preset organism's known mitochondrial gene symbols",
    "Specify manually",
]


def _render_custom_mito_resolution(project, custom_gtf_path):
    """
    For a CUSTOM reference whose GTF has zero mitochondrial genes found
    via direct seqname lookup, offer three explicit resolution paths --
    see this module's own docstring, "Custom-reference mitochondrial
    gene resolution UX", for the full rationale.

    Persists the resolved gene ID list (plus a "mito_source" label) into
    this project's saved reference_choice dict via
    scpm.save_reference_choice(), so Phase 2 doesn't need to re-resolve
    it on every Cell-level QC run.
    """
    st.markdown("**🧬 Mitochondrial Gene Resolution**")
    st.warning(
        "⚠️ No mitochondrial genes could be identified directly from this reference's GTF "
        "(by chromosome/contig lookup). A mitochondrial percentage of 0% for every cell in "
        "Phase 2's Cell-level QC would not be a meaningful biological result if left "
        "unresolved -- choose how to proceed:"
    )
    resolution_choice = st.radio(
        "How would you like to handle mitochondrial gene identification for this reference?",
        _MITO_RESOLUTION_OPTIONS, key="sc_custom_mito_resolution_choice",
    )

    resolved_gene_ids = []
    mito_source = "not_applicable"

    if resolution_choice == _MITO_RESOLUTION_OPTIONS[0]:
        st.caption("ℹ️ Mitochondrial QC will be skipped for this project -- mitochondrial % will always read 0%, and adaptive mito-based filtering will have no effect.")

    elif resolution_choice == _MITO_RESOLUTION_OPTIONS[1]:
        st.caption(
            "⚠️ **Lower confidence than a direct match.** This assumes your custom reference's "
            "gene-naming convention happens to overlap with the chosen preset species' naming -- "
            "not a verified identity. Different annotation sources, assembly versions, or a "
            "sufficiently divergent organism could produce few or no matches, or in rare cases "
            "an incorrect match. **Interpret mitochondrial QC metrics from this project with "
            "extra caution if you use this option.**"
        )
        species_options = list(ref.REFERENCE_CATALOG.keys())
        species_labels = {k: v["label"] for k, v in ref.REFERENCE_CATALOG.items()}
        # Only offer species whose reference has actually already been
        # downloaded on this system -- matching against a species that
        # hasn't been prepared yet would have nothing real to compare
        # against.
        available_species = [
            k for k in species_options
            if os.path.isdir(pm.shared_genome_dir(k)) and len(os.listdir(pm.shared_genome_dir(k))) > 0
        ]
        if not available_species:
            st.info("ℹ️ No preset species reference has been downloaded on this server yet -- prepare one via a Bulk RNA-Seq or single-cell project first to unlock this option.")
        else:
            preset_species_choice = st.selectbox(
                "Match against which preset species' known mitochondrial gene symbols?",
                options=available_species, format_func=lambda k: species_labels[k], key="sc_custom_mito_preset_species",
            )
            preset_gtf_path = os.path.join(pm.shared_genome_dir(preset_species_choice), f"{preset_species_choice}.annotation.gtf")
            preset_symbols = cellqc.get_mito_gene_symbols_from_gtf(preset_gtf_path)
            st.caption(f"Found {len(preset_symbols)} known mitochondrial gene symbol(s) for {species_labels[preset_species_choice]}.")
            if st.button("🔍 Search Custom Reference for Matching Gene Symbols", key="sc_custom_mito_symbol_match_btn"):
                matched = cellqc.match_custom_genes_by_mito_symbol(custom_gtf_path, preset_symbols)
                st.session_state["sc_custom_mito_symbol_matches"] = matched
            matched = st.session_state.get("sc_custom_mito_symbol_matches")
            if matched is not None:
                if matched:
                    st.success(f"✅ Found {len(matched)} matching gene(s) in your custom reference.")
                    resolved_gene_ids = matched
                    mito_source = f"custom_symbol_match:{preset_species_choice}"
                else:
                    st.warning("⚠️ No matching gene symbols found -- this reference's naming convention doesn't appear to overlap with the chosen preset species.")

    else:  # "Specify manually"
        st.caption(
            "Enter gene IDs or gene symbols (one per line, or comma-separated) known to be "
            "mitochondrial for your organism -- matched against your custom GTF's gene_id and "
            "gene_name fields (case-insensitive)."
        )
        manual_text = st.text_area("Gene IDs or symbols:", key="sc_custom_mito_manual_entry", height=100)
        if manual_text and st.button("🔍 Match Manual Entry Against Reference", key="sc_custom_mito_manual_match_btn"):
            manual_terms = [t.strip() for t in manual_text.replace(",", "\n").splitlines() if t.strip()]
            # Match manual entries against BOTH gene symbols (gene_name,
            # via match_custom_genes_by_mito_symbol) and, separately,
            # gene IDs across the ENTIRE GTF (via
            # find_gene_ids_matching_terms -- NOT the contig-restricted
            # get_mito_gene_ids_from_gtf, which would trivially return
            # empty here since that's exactly why this manual-entry path
            # was reached in the first place). A user may reasonably
            # supply either symbols or IDs.
            symbol_matched = cellqc.match_custom_genes_by_mito_symbol(custom_gtf_path, manual_terms)
            direct_id_matches = cellqc.find_gene_ids_matching_terms(custom_gtf_path, manual_terms)
            combined = sorted(set(symbol_matched) | set(direct_id_matches))
            st.session_state["sc_custom_mito_manual_matches"] = combined
        manual_matches = st.session_state.get("sc_custom_mito_manual_matches")
        if manual_matches is not None:
            if manual_matches:
                st.success(f"✅ Found {len(manual_matches)} matching gene(s).")
                resolved_gene_ids = manual_matches
                mito_source = "custom_manual_entry"
            else:
                st.warning("⚠️ No matching genes found for the entered terms.")

    if st.button("💾 Save Mitochondrial Gene Resolution", key="sc_custom_mito_save_btn"):
        existing_cfg = scpm.get_reference_choice(project) or {}
        existing_cfg["mito_gene_ids"] = resolved_gene_ids
        existing_cfg["mito_source"] = mito_source
        scpm.save_reference_choice(project, existing_cfg)
        st.success(f"✅ Saved ({len(resolved_gene_ids)} mitochondrial gene(s), source: {mito_source}).")


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

            # --- Mitochondrial genome verification (2026-08-17) ---
            _render_preset_mito_verification(species_choice, species_labels, genome_fasta_expected, gtf_expected)
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

        with st.expander("🔁 Force a fresh re-download (advanced -- generic, e.g. suspected corruption)"):
            st.caption(
                "⚠️ This reference is shared across every project using this species. "
                "Re-downloading replaces it for everyone. Use the dedicated mitochondrial-"
                "genome re-download button above instead if that's specifically your concern."
            )
            if st.button("🔄 Re-download Reference for All Projects", key="sc_ref_download_force_btn"):
                def _build_fn(temp_dir, _species=species_choice):
                    success, _paths, message = ref.download_genome_and_gtf(_species, temp_dir)
                    return success, message
                with st.spinner("Re-downloading..."):
                    success, message, _built = ref.ensure_shared_resource(target_dir, build_fn=_build_fn, force=True)
                if success:
                    st.success(f"✅ {message}")
                    st.rerun()
                else:
                    st.error(f"❌ {message}")

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

        if gtf_path and reference_cfg.get("annotation_format") == "gtf":
            mito_ids = cellqc.get_mito_gene_ids_from_gtf(gtf_path)
            if mito_ids:
                st.success(f"✅ {len(mito_ids)} mitochondrial gene(s) auto-detected directly from this GTF.")
                reference_cfg["mito_gene_ids"] = mito_ids
                reference_cfg["mito_source"] = "gtf_auto_detect"
            else:
                _render_custom_mito_resolution(project, gtf_path)

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

        if gtf_path and reference_cfg.get("annotation_format") == "gtf":
            mito_ids = cellqc.get_mito_gene_ids_from_gtf(gtf_path)
            if mito_ids:
                st.success(f"✅ {len(mito_ids)} mitochondrial gene(s) auto-detected directly from this GTF.")
                reference_cfg["mito_gene_ids"] = mito_ids
                reference_cfg["mito_source"] = "gtf_auto_detect"
            else:
                _render_custom_mito_resolution(project, gtf_path)

    _render_gff3_warning_if_needed(reference_cfg)

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
            # Preserve any already-resolved mito_gene_ids/mito_source
            # (e.g. set via _render_custom_mito_resolution's own save
            # button) rather than overwriting them with this call, which
            # doesn't itself carry that resolution state.
            prior_cfg = scpm.get_reference_choice(project) or {}
            for key in ("mito_gene_ids", "mito_source"):
                if key in prior_cfg and key not in reference_cfg:
                    reference_cfg[key] = prior_cfg[key]
            scpm.save_reference_choice(project, reference_cfg)
            scpm.mark_step_complete(project, "reference")
            st.success("✅ Reference setup confirmed.")


# ---------------------------------------------------------------------------
# Step 6: STARsolo alignment -- REAL EXECUTION
# ---------------------------------------------------------------------------

def _resolve_star_index_dir(project, reference_cfg):
    is_custom = reference_cfg.get("is_custom", False)
    species_key = reference_cfg.get("species_key")
    index_is_shared = not is_custom and bool(species_key)
    index_dir = pm.shared_star_index_dir(species_key) if index_is_shared else scpm.star_index_dir(project)
    # files_present alone (qm.star_index_exists -- checks for the
    # "SAindex" marker file) cannot tell a genuinely finished index
    # apart from an OLD index that's actively being replaced by a
    # force-rebuild in progress right now -- ensure_shared_resource()
    # builds into a private temp dir and only swaps it in via a single
    # atomic os.rename() at the very END of a successful build, so the
    # old index stays fully present (and looks "ready") for the entire
    # rebuild duration. See reference_manager.resource_build_in_progress()
    # for the full rationale (this was a real reported bug: navigating
    # away mid-rebuild and back showed a false "already built" status
    # while the real rebuild was still running, confirmed via `top`).
    #
    # This check now applies uniformly to BOTH shared (preset) and
    # custom (per-project) indexes, since _render_genome_index_step
    # below now routes BOTH build paths through
    # ref.ensure_shared_resource() -- previously only the shared path
    # had a lock file to probe at all.
    files_present = qm.star_index_exists(index_dir)
    ready = files_present and not ref.resource_build_in_progress(index_dir)
    return index_dir, ready, index_is_shared


def _render_genome_index_step(project, reference_cfg, pairs):
    """
    Genome index status + build UI -- shown in Step 5 immediately after
    genome selection, STARsolo-specific (gated on aligner choice).

    --- Force-rebuild bug fix (2026-08-17) ---
    build_clicked (no index yet) and force_rebuild (an index already
    exists, user explicitly wants it replaced) are tracked as TWO
    SEPARATE flags -- see this module's earlier docstring for the
    original flag-collapse regression this fixed.

    --- Stale-resource-during-rebuild race fix (2026-08-18) ---
    A real reported bug: a user force-rebuilding a shared STAR genome
    index navigated between radio-button pages while the rebuild was
    running (confirmed still running via `top` on the HPC host) and,
    on navigating back, saw the index reported as "✅ already built"
    even though the real rebuild was still in progress. Root cause:
    reference_manager.ensure_shared_resource() builds into a private
    temp directory and only replaces the OLD index via a single atomic
    os.rename() at the very END of a successful build -- so the old
    (stale) index stays fully present on disk, and any existence-only
    check (qm.star_index_exists) reports "ready" for the ENTIRE rebuild
    duration, with no way to distinguish a finished build from one
    that's actively being replaced.

    Fixed by using ref.resource_build_in_progress() (a non-blocking
    probe on the resource's existing .lock file) as part of the
    authoritative readiness check in _resolve_star_index_dir() above,
    and by showing an explicit "🔄 build in progress, please wait"
    status here that blocks further build actions while true.

    Also: CUSTOM (non-shared, per-project) index builds now route
    through ref.ensure_shared_resource() too, exactly like the shared
    path -- previously the custom branch called qm.build_star_index()
    DIRECTLY against index_dir with NO temp-dir staging, atomic
    replacement, or lock-file protection at all, so a custom
    reference's rebuild had the identical stale-appears-ready symptom
    with even less of a safety net than the shared path had. Both
    branches now share the same temp-directory + atomic-rename +
    lock-file infrastructure -- there is no longer a separate "shared"
    vs. "custom" code path for the actual build call, only for which
    directory is targeted.
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
    build_in_progress = ref.resource_build_in_progress(index_dir)

    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    index_threads = st.slider("Threads for indexing:", min_value=1, max_value=detected_cores, value=recommended_threads, key="sc_star_index_threads")

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

    if build_in_progress:
        st.warning(
            "🔄 **A genome index build/rebuild is currently in progress** (started by "
            "this or another session/pipeline) -- please wait. You cannot proceed to "
            "alignment until it finishes; re-check this page periodically rather than "
            "assuming a quick 'complete' means it actually is. (A genuine rebuild for a "
            "full genome legitimately takes several minutes or more -- if this status "
            "disappears within seconds, something may not have started correctly.)"
        )
    elif index_ready:
        st.success("✅ STAR genome index already built " + ("(shared across every project using this species)." if index_is_shared else "for this project's reference."))
    else:
        st.info("ℹ️ **No genome index found yet -- build it below before alignment can run.** This is a one-time step per reference and can take a while for larger genomes.")

    # --- THE FIX: two separate flags, not one -- AND both builds now
    # route through ensure_shared_resource() (see below) ---
    build_clicked = False
    force_rebuild = False
    if build_in_progress:
        pass  # no build controls shown while a build is already running
    elif not index_ready:
        build_clicked = st.button("🔧 Build Genome Index", key="sc_build_star_index_btn", type="primary")
    else:
        with st.expander("🔁 Force a fresh re-index (advanced)"):
            st.caption(
                "⚠️ Only do this if you have a specific reason to believe the current index is "
                "corrupted or out of date (e.g. the reference's genome/annotation was "
                "re-downloaded since this index was built -- such as after using the "
                "mitochondrial-genome re-download button above)."
                + (" This index is shared across every project using this species -- rebuilding "
                   "replaces it for everyone." if index_is_shared else "")
            )
            if st.button("🔄 Re-build Index", key="sc_build_star_index_force_btn"):
                force_rebuild = True
            st.caption(
                "ℹ️ A genuine STAR index rebuild for a full genome (e.g. human/mouse) legitimately "
                "takes several minutes -- if this completes in just a few seconds, something did "
                "not rebuild correctly; please report this."
            )

    if build_clicked or force_rebuild:
        def _build_index_impl(temp_dir):
            return qm.build_star_index(genome_fasta, gtf_path, temp_dir, threads=index_threads, sjdb_overhang=sjdb_overhang)

        # Both shared (preset) AND custom (per-project) index builds
        # now go through ensure_shared_resource() -- previously only
        # the shared branch got temp-dir + atomic-rename + lock
        # protection here; a custom reference's rebuild wrote STAR's
        # output directly into index_dir with no protection at all.
        # See this function's docstring for the full rationale.
        status_placeholder = st.empty()

        def _wait_cb(elapsed):
            status_placeholder.info(f"⏳ Another project/session is currently building this index. Waiting... ({int(elapsed)}s elapsed)")

        with st.spinner("Building genome index... this may take a while for larger genomes."):
            success, message, built = ref.ensure_shared_resource(
                index_dir, _build_index_impl, wait_message_callback=_wait_cb, force=force_rebuild,
            )
        status_placeholder.empty()

        if success:
            st.success("✅ Genome index built successfully.")
            if force_rebuild:
                st.warning(
                    "⚠️ **Important:** rebuilding the index does NOT automatically re-run Step 6 "
                    "alignment for any sample. Revisit Step 6 and re-run alignment for every "
                    "sample before this fix will actually show up in downstream results (e.g. "
                    "Cell-level QC's mitochondrial gene detection)."
                )
        else:
            st.error("Genome indexing failed. Details below:")
            st.code(message)


def _render_mapping_rate_qc(results_df, rate_column):
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

    index_dir, index_ready, index_is_shared = _resolve_star_index_dir(project, reference_cfg)
    if not index_ready:
        st.warning("⚠️ No genome index has been built yet for this reference. Revisit **Step 5** (right after genome selection) to build it before alignment can run.")
        return

    st.markdown("---")

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
            scpm.save_alignment_results(project, results)
            st.success("✅ STARsolo alignment & cell-calling complete for at least one sample.")

        _render_mapping_rate_qc(results_df, rate_column="Uniquely Mapped (%)")


def _render_step6(project, pairs):
    st.subheader("Step 6: Alignment, Quantification & Cell-calling (STARsolo)")

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


# ---------------------------------------------------------------------------
# Phase 2: Cell-level QC visualizations
# ---------------------------------------------------------------------------
def _render_stale_index_diagnostic(diagnostic, sample_name):
    """
    Render sc_cellqc_manager.diagnose_starsolo_matrix_for_mito()'s result
    -- a fast, pre-flight check of THIS sample's already-aligned matrix,
    run BEFORE the (much slower) full Cell-level QC pipeline, so a
    stale-index situation is caught and explained immediately rather than
    only discovered after re-running the whole R pipeline again. See this
    module's own docstring, "Genome-index force-rebuild bug fix +
    stale-index diagnostic", for the full rationale.
    """
    if diagnostic is None:
        return  # no features.tsv found -- _render_step6/filtered_dir check already covers this
    if not diagnostic["likely_stale_index"]:
        return  # matrix has mito genes -- nothing to warn about here

    st.error(
        f"🔴 **This sample's (`{sample_name}`) already-aligned matrix has "
        f"{diagnostic['total_genes']:,} genes total, but ZERO of them are mitochondrial.** "
        "This is the specific signature of a **stale STAR index/alignment** -- Cell-level QC "
        "only re-reads this existing matrix; it does NOT re-align, so running (or re-running) "
        "QC on this sample **cannot fix this** no matter what mitochondrial gene resolution is "
        "configured in Step 5.\n\n"
        "**To actually fix this:**\n"
        "1. Confirm the reference itself has mitochondrial content (Step 5's verification, above "
        "the reference selection).\n"
        "2. If needed, use Step 5's **'🔁 Force a fresh re-index'** to genuinely rebuild the STAR "
        "genome index (a real rebuild for a full genome takes several minutes -- if it finishes "
        "in seconds, it did not actually rebuild).\n"
        "3. Revisit **Step 6** and **re-run alignment** for this sample -- this regenerates the "
        "matrix with the corrected gene set.\n"
        "4. **Then** re-run Cell-level QC here."
    )


def _render_mito_diagnostic(diagnostic):
    if diagnostic is None:
        st.caption(
            "ℹ️ This project's cell-level QC run predates the mitochondrial-gene detection "
            "diagnostic -- re-run to see this breakdown."
        )
        return

    symbol_count = diagnostic.get("symbol_match_count", 0)
    gtf_count = diagnostic.get("gtf_derived_count", 0)
    union_count = diagnostic.get("union_count", 0)
    mito_source = diagnostic.get("mito_source", "gtf_auto_detect")

    if diagnostic.get("warning") == "no_mito_genes_found" or union_count == 0:
        st.error(
            "🔴 **Zero mitochondrial genes were identified in this reference, by either detection "
            "method.** A mitochondrial percentage of 0% for every cell is almost never a real "
            "biological result -- it's the expected symptom of this. This most likely means the "
            "**reference genome/GTF used for this project's alignment does not include the "
            "mitochondrial chromosome/contig at all**, its chromosome-naming convention wasn't "
            "recognized, OR this sample's STAR index/alignment is STALE (built before the "
            "reference's mitochondrial content was confirmed/fixed) -- see the diagnostic above, "
            "if shown, for that specific check. **Mitochondrial QC filtering is not meaningful "
            "until this is resolved.**\n\n"
            f"- Gene-symbol ('^MT-' prefix) matches: {symbol_count}\n"
            f"- GTF/override-derived matches: {gtf_count}\n"
        )
    else:
        sample_ids = diagnostic.get("sample_gene_ids") or []
        sample_note = f" (e.g. `{'`, `'.join(sample_ids[:3])}`)" if sample_ids else ""
        source_note = {
            "gtf_auto_detect": "direct GTF chromosome lookup",
        }.get(mito_source, mito_source.replace("_", " "))
        st.success(
            f"🟢 **{union_count} mitochondrial gene(s)** identified in this reference{sample_note} "
            f"-- {symbol_count} via gene-symbol matching, {gtf_count} via {source_note}."
        )


def _render_cellqc_visualizations(qc_df, output_dir, doublet_method_used):
    st.markdown("### 📊 QC Visualizations")

    figures = {}

    st.markdown("**Per-cell QC metric distributions**")
    violin_cols = st.columns(3)
    metric_specs = [
        ("sum", "Total UMI counts", violin_cols[0]),
        ("detected", "Genes detected", violin_cols[1]),
        ("subsets_Mito_percent", "Mitochondrial %", violin_cols[2]),
    ]
    for col_name, label, col in metric_specs:
        if col_name not in qc_df.columns:
            continue
        fig = go.Figure()
        fig.add_trace(go.Violin(y=qc_df[col_name], box_visible=True, points="outliers", name=label, meanline_visible=True))
        fig.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10), title=label, showlegend=False)
        with col:
            st.plotly_chart(fig, use_container_width=True)
        figures[f"qc_violin_{col_name}"] = fig

    st.caption(
        "Standard first-look QC plots: too-low counts/genes suggest empty droplets or "
        "low-quality cells; unusually high counts/genes can indicate doublets; elevated "
        "mitochondrial % suggests stressed or dying cells."
    )

    if all(c in qc_df.columns for c in ("sum", "detected", "subsets_Mito_percent")):
        st.markdown("**Total counts vs. genes detected** (colored by mitochondrial %)")
        scatter_fig = px.scatter(
            qc_df, x="sum", y="detected", color="subsets_Mito_percent",
            color_continuous_scale="Inferno_r",
            labels={"sum": "Total UMI counts", "detected": "Genes detected", "subsets_Mito_percent": "Mito %"},
            opacity=0.6,
        )
        scatter_fig.update_layout(height=420, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(scatter_fig, use_container_width=True)
        figures["qc_scatter_counts_vs_genes"] = scatter_fig
        st.caption(
            "Cells that are both low-count AND high-mito (bright points in the lower-left) are "
            "the strongest candidates for exclusion -- low complexity alone can still be "
            "biologically real for some cell types."
        )

    if "doublet_score" in qc_df.columns:
        st.markdown(f"**Doublet score distribution** ({doublet_method_used})")
        hist_fig = px.histogram(
            qc_df, x="doublet_score", color="predicted_doublet" if "predicted_doublet" in qc_df.columns else None,
            nbins=50, labels={"doublet_score": "Doublet score", "predicted_doublet": "Predicted doublet"},
            color_discrete_map={True: "#d62728", False: "#1f77b4"},
        )
        hist_fig.update_layout(height=380, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(hist_fig, use_container_width=True)
        figures["doublet_score_histogram"] = hist_fig
        st.caption(
            "A healthy doublet score distribution is typically **bimodal** -- most cells scored "
            "near 0 (singlets), with a distinct, separated group scored near 1 (real doublets)."
        )

        if "sum" in qc_df.columns and "predicted_doublet" in qc_df.columns:
            st.markdown("**Total UMI counts, split by doublet call**")
            split_fig = go.Figure()
            for is_doublet, label, color in [(False, "Singlet", "#1f77b4"), (True, "Doublet", "#d62728")]:
                subset = qc_df[qc_df["predicted_doublet"] == is_doublet]
                if not subset.empty:
                    split_fig.add_trace(go.Violin(y=subset["sum"], name=label, box_visible=True, meanline_visible=True, line_color=color))
            split_fig.update_layout(height=380, margin=dict(l=10, r=10, t=30, b=10), yaxis_title="Total UMI counts")
            st.plotly_chart(split_fig, use_container_width=True)
            figures["doublet_split_violin"] = split_fig
            st.caption(
                "Predicted doublets should generally show HIGHER total UMI counts than singlets."
            )

    if "ambient_contamination" in qc_df.columns and qc_df["ambient_contamination"].notna().any():
        st.markdown("**Ambient RNA contamination fraction distribution**")
        contam_fig = px.histogram(qc_df, x="ambient_contamination", nbins=40, labels={"ambient_contamination": "Contamination fraction"})
        contam_fig.update_layout(height=350, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(contam_fig, use_container_width=True)
        figures["ambient_contamination_histogram"] = contam_fig
        st.caption(
            "DecontX reports a genuine per-cell distribution here; SoupX applies one global "
            "estimate to every cell, so this histogram will show a single spike instead."
        )

    top_ambient_df = cellqc.read_top_ambient_genes(output_dir)
    if top_ambient_df is not None and not top_ambient_df.empty:
        st.markdown("**Top genes contributing to ambient RNA signal**")
        bar_fig = px.bar(
            top_ambient_df.head(15).sort_values("counts_removed"),
            x="counts_removed", y="symbol", orientation="h",
            labels={"counts_removed": "Total UMI counts removed", "symbol": "Gene"},
        )
        bar_fig.update_layout(height=420, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(bar_fig, use_container_width=True)
        figures["top_ambient_genes_bar"] = bar_fig
        st.caption(
            "These are typically a small number of highly-expressed marker genes from the most "
            "abundant cell type(s) in the sample, leaking into every droplet as background 'soup'."
        )

    return figures


def _build_qc_package_zip(sample_name, qc_df, output_dir, mito_diagnostic, thresholds,
                           doublet_method, ambient_method, figures):
    """
    Bundle a completed Cell-level QC run's full results into a single
    in-memory zip file -- see this module's own docstring, "Downloadable
    QC results package", for the full rationale.

    Includes: the raw per-cell metrics CSV, QC thresholds JSON, mito
    diagnostic JSON, top-ambient-genes CSV, DoubletFinder's pK sweep CSV
    (if this run used it), a self-contained markdown summary, and a PNG
    render of every QC visualization figure passed in (via Plotly's
    kaleido-based image export). Image export failures (e.g. kaleido
    unavailable for some reason) are caught individually and noted in
    the summary rather than failing the whole package.

    Returns the zip file's raw bytes, ready for st.download_button.
    """
    buffer = io.BytesIO()
    summary = cellqc.summarize_cellqc_results(qc_df)

    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("cell_qc_metrics.csv", qc_df.to_csv(index=False))

        if thresholds:
            import json as _json
            zf.writestr("qc_thresholds.json", _json.dumps(thresholds, indent=2))
        if mito_diagnostic:
            import json as _json
            zf.writestr("mito_gene_diagnostic.json", _json.dumps(mito_diagnostic, indent=2))

        top_ambient_df = cellqc.read_top_ambient_genes(output_dir)
        if top_ambient_df is not None:
            zf.writestr("top_ambient_genes.csv", top_ambient_df.to_csv(index=False))

        pk_sweep_df = cellqc.read_doubletfinder_pk_sweep(output_dir)
        if pk_sweep_df is not None:
            zf.writestr("doubletfinder_pk_sweep.csv", pk_sweep_df.to_csv(index=False))

        summary_md = cellqc.build_qc_summary_markdown(
            sample_name, summary, mito_diagnostic, thresholds, doublet_method, ambient_method,
        )
        zf.writestr("SUMMARY.md", summary_md)

        image_errors = []
        for fig_name, fig in (figures or {}).items():
            try:
                png_bytes = fig.to_image(format="png", width=1000, height=600, scale=2)
                zf.writestr(f"plots/{fig_name}.png", png_bytes)
            except Exception as e:
                image_errors.append(f"{fig_name}: {e}")
        if image_errors:
            zf.writestr(
                "plots/_export_errors.txt",
                "Some plots could not be exported as images (kaleido issue?):\n" + "\n".join(image_errors),
            )

    buffer.seek(0)
    return buffer.getvalue()


def render_cell_qc():
    """
    Phase 2 entry point -- Cell-level QC (scDblFinder/DoubletFinder
    doublet detection, DecontX/SoupX ambient RNA correction, adaptive
    per-cell filtering with standard visualizations and a
    downloadable results package).
    """
    st.title("🔬 Single-cell Cell-level QC")
    st.markdown(
        "Doublet detection, ambient RNA correction, and adaptive per-cell "
        "filtering -- applied to STARsolo's filtered cell x gene matrix "
        "from Step 6."
    )
    st.markdown("---")
    project, pairs = _require_project_with_source_dir()
    if not pairs:
        return
    if not scpm.has_completed_step(project, "alignment"):
        st.warning("⚠️ Complete **Step 6** (🧬 SC Alignment & Cell-Calling) first -- this page needs a completed STARsolo run's filtered cell x gene matrix to work from.")
        return

    if not cellqc.cellqc_tools_available():
        st.error(
            "⚠️ Rscript was not found on this system. R with the DropletUtils, scuttle, "
            "scDblFinder, and celda (DecontX) [and, if selected, SoupX] packages needs to be "
            "installed in your environment before this step can run."
        )
        return

    align_dir = scpm.starsolo_output_dir(project)
    sample_names = list(pairs.keys())
    if not sample_names:
        st.warning("⚠️ No samples found for this project.")
        return
    sample_name = st.selectbox("Sample to run cell-level QC on:", options=sample_names, key="sc_cellqc_sample_select")

    sample_out_dir = os.path.join(align_dir, sample_name)
    output_prefix = os.path.join(sample_out_dir, f"{sample_name}_")
    filtered_dir = star.filtered_counts_matrix_dir(output_prefix)
    raw_dir = star.counts_matrix_dir(output_prefix)

    if not os.path.isdir(filtered_dir):
        st.error(f"⚠️ Could not find STARsolo's filtered matrix directory for `{sample_name}` at `{filtered_dir}`. Re-run Step 6 for this sample if this is unexpected.")
        return

    # --- Stale-index pre-flight diagnostic (2026-08-17) ---
    # Fast, R-independent check run the MOMENT a sample is selected --
    # before the (much slower) full Cell-level QC R pipeline below --
    # see this module's own docstring, "Genome-index force-rebuild bug
    # fix + stale-index diagnostic", for the full rationale.
    stale_diagnostic = cellqc.diagnose_starsolo_matrix_for_mito(filtered_dir)
    _render_stale_index_diagnostic(stale_diagnostic, sample_name)

    reference_cfg = scpm.get_reference_choice(project) or {}
    mito_gtf_path = reference_cfg.get("custom_gtf")
    # Prefer an already-resolved override list (e.g. from custom-
    # reference symbol matching or manual entry in Step 5) over a fresh
    # GTF auto-parse -- see run_cellqc_analysis()'s own docstring.
    mito_gene_ids_override = reference_cfg.get("mito_gene_ids")
    mito_source = reference_cfg.get("mito_source", "gtf_auto_detect")

    st.markdown("---")
    st.markdown("**🧬 Doublet Detection**")
    doublet_method_keys = list(cellqc.DOUBLET_METHOD_OPTIONS.keys())
    doublet_method = st.radio(
        "Doublet detection method:", doublet_method_keys,
        format_func=lambda k: cellqc.DOUBLET_METHOD_OPTIONS[k]["label"],
        index=doublet_method_keys.index(cellqc.DEFAULT_DOUBLET_METHOD), key="sc_doublet_method_choice",
    )
    st.caption(cellqc.DOUBLET_METHOD_OPTIONS[doublet_method]["explanation"])
    if not cellqc.DOUBLET_METHOD_OPTIONS[doublet_method]["implemented"]:
        st.warning("⚠️ This method isn't implemented yet -- select scDblFinder to continue.")
        return

    simulation_mode = cellqc.DEFAULT_SIMULATION_MODE
    doubletfinder_n_pcs = cellqc.DEFAULT_DOUBLETFINDER_N_PCS
    doubletfinder_cluster_resolution = cellqc.DEFAULT_DOUBLETFINDER_CLUSTER_RESOLUTION

    if doublet_method == "scDblFinder":
        sim_mode_keys = list(cellqc.DOUBLET_SIMULATION_MODES.keys())
        simulation_mode = st.radio(
            "Doublet simulation mode:", sim_mode_keys,
            format_func=lambda k: cellqc.DOUBLET_SIMULATION_MODES[k]["label"],
            index=sim_mode_keys.index(cellqc.DEFAULT_SIMULATION_MODE), key="sc_doublet_sim_mode_choice",
        )
        st.caption(cellqc.DOUBLET_SIMULATION_MODES[simulation_mode]["explanation"])
    else:
        st.info(
            "⏱️ **DoubletFinder runs noticeably slower than scDblFinder** -- it computes its own "
            "internal PCA/clustering, then sweeps many candidate 'pK' parameter values to find the "
            "best one automatically. Expect this to take several minutes."
        )
        with st.expander("⚙️ DoubletFinder internal preprocessing settings (advanced)"):
            st.caption(
                "These control DoubletFinder's own required internal PCA/clustering step -- "
                "**not** the pipeline's future Phase 3 clustering feature."
            )
            doubletfinder_n_pcs = st.slider(
                "Number of principal components:", min_value=5, max_value=50,
                value=cellqc.DEFAULT_DOUBLETFINDER_N_PCS, key="sc_doubletfinder_n_pcs",
            )
            doubletfinder_cluster_resolution = st.slider(
                "Clustering resolution:", min_value=0.1, max_value=2.0,
                value=cellqc.DEFAULT_DOUBLETFINDER_CLUSTER_RESOLUTION, step=0.1, key="sc_doubletfinder_cluster_res",
            )

    persisted_results = scpm.get_alignment_results(project) or []
    cells_detected = next((r.get("Cells Detected") for r in persisted_results if r.get("Sample") == sample_name), None)
    default_dbr = cellqc.compute_expected_doublet_rate(cells_detected) if isinstance(cells_detected, (int, float)) else 0.008
    expected_doublet_rate = st.slider(
        "Expected doublet rate:", min_value=0.0, max_value=0.30, value=float(default_dbr), step=0.001,
        format="%.3f", key="sc_expected_doublet_rate",
        help=(
            f"Auto-computed as ~0.8% per 1,000 cells loaded" +
            (f" ({cells_detected:,} cells detected in Step 6 -> {default_dbr:.3f})" if cells_detected else "") +
            "."
        ),
    )

    st.markdown("---")
    st.markdown("**🧪 Ambient RNA Correction**")
    ambient_method_keys = list(cellqc.AMBIENT_METHOD_OPTIONS.keys())
    ambient_method = st.radio(
        "Ambient RNA correction method:", ambient_method_keys,
        format_func=lambda k: cellqc.AMBIENT_METHOD_OPTIONS[k]["label"],
        index=ambient_method_keys.index(cellqc.DEFAULT_AMBIENT_METHOD), key="sc_ambient_method_choice",
    )
    st.caption(cellqc.AMBIENT_METHOD_OPTIONS[ambient_method]["explanation"])
    if ambient_method == "soupx" and not os.path.isdir(raw_dir):
        st.error(f"⚠️ SoupX requires STARsolo's raw (unfiltered) matrix, which wasn't found at `{raw_dir}`.")
        return

    st.markdown("---")
    st.markdown("**🎚️ Per-cell Filtering Thresholds**")
    st.caption(
        "Adaptive thresholds (3 median-absolute-deviations from the median) are computed "
        "automatically per-sample for total counts, genes detected, and mitochondrial %."
    )
    nmads = st.slider(
        "MAD sensitivity (lower = stricter):", min_value=1.0, max_value=5.0, value=float(cellqc.DEFAULT_MAD_NMADS), step=0.5,
        key="sc_cellqc_nmads",
    )
    if reference_cfg.get("is_custom") and mito_gene_ids_override is not None:
        st.caption(f"ℹ️ Using {len(mito_gene_ids_override)} previously-resolved mitochondrial gene(s) for this custom reference (source: {mito_source}).")

    output_dir = os.path.join(sample_out_dir, "cellqc")
    work_dir = os.path.join(sample_out_dir, "cellqc_work")

    spinner_text = (
        "Running DoubletFinder (PCA + clustering + pK sweep) and ambient RNA correction... "
        "this can take several minutes." if doublet_method == "doublet_finder" else
        "Running doublet detection + ambient RNA correction... this may take a few minutes."
    )
    if st.button("▶️ Run Cell-level QC", key="sc_run_cellqc_btn", type="primary"):
        with st.spinner(spinner_text):
            success, log = cellqc.run_cellqc_analysis(
                filtered_matrix_dir=filtered_dir, output_dir=output_dir, work_dir=work_dir,
                doublet_method=doublet_method, simulation_mode=simulation_mode,
                expected_doublet_rate=expected_doublet_rate, ambient_method=ambient_method,
                raw_matrix_dir=raw_dir if ambient_method == "soupx" else None, nmads=nmads,
                mito_gtf_path=mito_gtf_path, mito_gene_ids_override=mito_gene_ids_override,
                mito_source=mito_source,
                doubletfinder_n_pcs=doubletfinder_n_pcs,
                doubletfinder_cluster_resolution=doubletfinder_cluster_resolution,
            )
        if not success:
            st.error("Cell-level QC failed. Details below:")
            st.code(log)
        else:
            st.success("✅ Cell-level QC completed.")
            with st.expander("📜 Run log"):
                st.code(log)

    qc_df = cellqc.read_cell_qc_metrics(output_dir)
    if qc_df is not None:
        st.markdown("---")
        st.markdown("**📊 Results**")

        mito_diagnostic = cellqc.read_mito_gene_diagnostic(output_dir)
        _render_mito_diagnostic(mito_diagnostic)

        summary = cellqc.summarize_cellqc_results(qc_df)
        if summary:
            st.markdown(summary["messages"]["adaptive_qc"])
            st.markdown(summary["messages"]["doublet"])
            st.markdown(summary["messages"]["ambient"])
        thresholds = cellqc.read_qc_thresholds(output_dir)
        if thresholds:
            st.caption(
                f"Adaptive thresholds used: total counts > {thresholds.get('sum_lower', 0):.0f}, "
                f"genes detected > {thresholds.get('detected_lower', 0):.0f}, "
                f"mitochondrial % < {thresholds.get('mito_upper', 0):.1f}%"
            )

        pk_sweep_df = cellqc.read_doubletfinder_pk_sweep(output_dir)
        if pk_sweep_df is not None and not pk_sweep_df.empty:
            with st.expander("🔍 DoubletFinder pK parameter sweep details"):
                st.caption(
                    "Each row is a candidate 'pK' value DoubletFinder tested; the one with the "
                    "highest **BCmetric** was automatically selected and used for the final "
                    "classification above."
                )
                if "doubletfinder_pK_used" in qc_df.columns and not qc_df["doubletfinder_pK_used"].empty:
                    st.caption(f"✅ Selected pK: **{qc_df['doubletfinder_pK_used'].iloc[0]}**")
                chart_df = pk_sweep_df.copy()
                if "pK" in chart_df.columns:
                    chart_df["pK"] = chart_df["pK"].astype(str)
                    chart_df = chart_df.set_index("pK")
                if "BCmetric" in chart_df.columns:
                    st.bar_chart(chart_df["BCmetric"])
                st.dataframe(pk_sweep_df, use_container_width=True, hide_index=True)

        st.markdown("---")
        figures = _render_cellqc_visualizations(qc_df, output_dir, doublet_method)

        st.markdown("---")
        st.markdown("**📦 Download Full QC Package**")
        st.caption(
            "Bundles the per-cell metrics table, QC thresholds, mitochondrial-gene diagnostic, "
            "top-ambient-genes table, DoubletFinder's pK sweep (if used), a plain-text summary, "
            "and PNG copies of every plot above into a single .zip file."
        )
        try:
            zip_bytes = _build_qc_package_zip(
                sample_name, qc_df, output_dir, mito_diagnostic, thresholds,
                doublet_method, ambient_method, figures,
            )
            st.download_button(
                "📦 Download QC Package (.zip)", data=zip_bytes,
                file_name=f"{sample_name}_cellqc_package.zip", mime="application/zip",
                key="sc_cellqc_download_package_btn",
            )
        except Exception as e:
            st.error(f"⚠️ Could not build the QC package: {e}")

        st.markdown("---")
        with st.expander("📋 Full per-cell QC table"):
            st.dataframe(qc_df, use_container_width=True, hide_index=True)


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
