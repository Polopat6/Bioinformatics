"""
trimming_workspace.py

Adapter Trimming & Post-Trimming Quality Control workspace.

This picks up exactly where bulk_rnaseq_workspace.py leaves off: it reuses
the *same active project* (matched samples + FastQC/MultiQC results
already on disk) rather than requiring anything to be re-uploaded.

Trimming is done with fastp, which conveniently generates its own
per-sample QC report (a .json + .html) as a side effect of trimming —
no separate FastQC run is needed post-trim. MultiQC has a built-in parser
for fastp's JSON output, so we get the same "one combined report" result
as Step 4 of the Bulk RNA-Seq workspace, just sourced from fastp instead
of FastQC.

Design goal: same as bulk_rnaseq_workspace.py — assume the user has little
to no bioinformatics background, explain each step in plain language.

This module is fully self-contained. All trimming/post-trim QC development
should happen here — editing this file has zero effect on
spatial_workspace.py or bulk_rnaseq_workspace.py.
"""

import os

import pandas as pd
import streamlit as st

import project_manager as pm
import fastp_manager as fastp

# Reuse the same workspace_key as bulk_rnaseq_workspace.py so the active
# project selection (st.session_state["bulk_rnaseq_project"]) is shared
# automatically between the two pages — no need to re-pick a project when
# jumping here via the "Proceed to Trimming" button.
WORKSPACE_KEY = "bulk_rnaseq"

# Known adapter sequences (or distinctive prefixes) for common library
# prep kits, used to translate fastp's raw detected adapter sequence into
# a recognizable kit name for users without a bioinformatics background.
# Matching is done on a prefix basis since fastp may report slightly
# different lengths depending on read length / overlap detection.


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("🧪 Adapter Trimming & Post-Trimming QC")
    st.markdown(
        "This workspace trims low-quality bases and leftover adapter "
        "sequences from your reads, then re-checks quality afterward to "
        "confirm the trimming actually helped. **No bioinformatics "
        "experience required** — follow the steps below in order."
    )
    st.markdown("---")

    # -----------------------------------------------------------------
    # Project selection — shared with the Bulk RNA-Seq workspace
    # -----------------------------------------------------------------
    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)

    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    st.markdown("---")

    # -----------------------------------------------------------------
    # Gate: require matched samples + QC to already be done
    # -----------------------------------------------------------------
    samplesheet_path = pm.samplesheet_path(project)
    qc_done = pm.has_completed_step(project, "qc_complete")

    if not os.path.exists(samplesheet_path):
        st.warning(
            "⚠️ This project doesn't have a matched sample list yet. "
            "Go to the **🧬 Bulk RNA-Seq Pipeline** page first, upload "
            "your FASTQ files and metadata, and save your matched sample "
            "list before trimming."
        )
        # Note: navigation uses the "nav_request" indirection (see
        # app.py) rather than setting st.session_state["assay_choice_radio"]
        # directly, since that widget's key is already instantiated
        # earlier in this same run.
        if st.button("⬅️ Go to Bulk RNA-Seq Pipeline", key="trim_gate1_back_btn"):
            st.session_state["nav_request"] = "🧬 Bulk RNA-Seq Pipeline"
            st.rerun()
        return

    if not qc_done:
        st.warning(
            "⚠️ Quality control hasn't been run for this project yet. "
            "We recommend reviewing raw read quality before trimming, so "
            "you know what (if anything) actually needs to be trimmed. "
            "Go to the **🧬 Bulk RNA-Seq Pipeline** page and run Step 4 "
            "(FastQC + MultiQC) first."
        )
        if st.button("⬅️ Go to Bulk RNA-Seq Pipeline", key="trim_gate2_back_btn"):
            st.session_state["nav_request"] = "🧬 Bulk RNA-Seq Pipeline"
            st.rerun()
        return

    samplesheet_df = pd.read_csv(samplesheet_path)
    st.success(f"✅ Using matched samples from project `{project}` ({len(samplesheet_df)} sample(s)).")
    st.dataframe(samplesheet_df, use_container_width=True, hide_index=True)

    trimmed_dir = pm.trimmed_fastq_dir(project)
    fastp_reports_dir = pm.fastp_reports_dir(project)
    posttrim_multiqc_dir = pm.posttrim_multiqc_dir(project)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 1: Trim adapters with fastp
    # -----------------------------------------------------------------
    st.header("Step 1: Trim Adapters & Low-Quality Bases")

    with st.expander("ℹ️ What is trimming, and why do it? (click to learn more)"):
        st.markdown(
            "**Trimming** removes parts of your reads that could "
            "interfere with downstream analysis:\n"
            "- Leftover **adapter sequences** (short synthetic sequences "
            "used during library prep that sometimes get sequenced by "
            "accident)\n"
            "- **Low-quality bases**, usually near the ends of reads\n\n"
            "Trimming isn't always necessary — if your earlier quality "
            "report came back clean, you may be able to skip straight to "
            "alignment. But if you saw adapter content warnings or poor "
            "end-of-read quality, trimming can meaningfully improve your "
            "results.\n\n"
            "This tool uses **fastp**, which automatically detects and "
            "removes adapters and trims low-quality bases using sensible "
            "defaults — no manual configuration needed. It also produces "
            "its own quality report for each sample as it trims, which "
            "we combine into a single report in Step 2."
        )

    fastp_ok, multiqc_ok = fastp.tools_available()
    trimming_done = pm.has_completed_step(project, "trimming_complete")

    if not fastp_ok:
        st.error(
            "⚠️ fastp was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
    else:
        # Auto-detected default thread count, same pattern used for
        # FastQC in bulk_rnaseq_workspace.py — recommends a sensible
        # starting point based on this machine's actual CPU core count
        # rather than a fixed guess. fastp's own documentation notes
        # diminishing returns beyond ~16 threads (I/O-bound, not
        # CPU-bound past that point), so the slider is capped there
        # regardless of how many cores are detected.
        detected_cores, recommended_threads = pm.get_recommended_thread_count(max_default=16)
        fastp_threads = st.slider(
            "fastp threads (per sample):",
            min_value=1, max_value=16, value=min(recommended_threads, 16),
            help=(
                f"How many threads fastp uses while processing each "
                f"sample. Detected {detected_cores} CPU core(s) on this "
                f"machine, so {min(recommended_threads, 16)} is suggested. "
                "fastp itself reports diminishing returns beyond ~16 "
                "threads regardless of how many cores your machine has, "
                "since I/O (reading/writing files) becomes the limiting "
                "factor rather than CPU power at that point."
            ),
            key="fastp_threads_slider",
        )

        # --- Determine which samples still need trimming ---
        # Checked every render (not just on click) so the "already
        # trimmed" vs. "still needed" counts shown to the user are always
        # accurate and up to date, e.g. immediately reflecting a prior
        # partial run (some samples succeeded, some failed/were
        # interrupted by a disk-space error) without requiring a click
        # first just to see the current state.
        already_trimmed_rows = []
        needs_trimming_rows = []
        for _, row in samplesheet_df.iterrows():
            if fastp.sample_already_trimmed(row, trimmed_dir):
                already_trimmed_rows.append(row)
            else:
                needs_trimming_rows.append(row)

        if already_trimmed_rows and needs_trimming_rows:
            st.info(
                f"ℹ️ **{len(already_trimmed_rows)} of {len(samplesheet_df)} "
                f"sample(s)** already have trimmed output on disk and "
                f"will be **skipped** by default (no need to redo work "
                f"that already succeeded — this also avoids using disk "
                f"space re-writing files that are already there). "
                f"**{len(needs_trimming_rows)} sample(s)** still need "
                "trimming."
            )
        elif already_trimmed_rows and not needs_trimming_rows and trimming_done:
            st.success("✅ Trimming has already been run for this project — all samples have trimmed output on disk.")

        force_retrim_all = False
        if already_trimmed_rows:
            force_retrim_all = st.checkbox(
                "Force re-trim ALL samples (including ones already done)",
                value=False,
                help=(
                    "Check this only if you specifically want to redo "
                    "trimming for every sample — for example, after "
                    "changing the thread count above shouldn't require "
                    "this (thread count doesn't change the output), but "
                    "you might want this after deciding to use different "
                    "fastp settings in the future. Leaving this unchecked "
                    "(recommended) only processes samples that don't "
                    "already have trimmed output on disk."
                ),
                key="force_retrim_all_checkbox",
            )

        rows_to_process = list(samplesheet_df.iterrows()) if force_retrim_all else [
            (i, row) for i, row in samplesheet_df.iterrows() if row["sample"] in
            {r["sample"] for r in needs_trimming_rows}
        ]
        n_to_process = len(rows_to_process)

        if n_to_process == 0 and not force_retrim_all:
            trim_button_label = "✅ All Samples Already Trimmed"
            trim_button_disabled = True
        elif trimming_done or already_trimmed_rows:
            trim_button_label = f"✂️ Trim Remaining {n_to_process} Sample(s)" if not force_retrim_all else f"🔄 Force Re-trim All {n_to_process} Sample(s)"
            trim_button_disabled = False
        else:
            trim_button_label = "✂️ Trim All Matched Samples"
            trim_button_disabled = False

        trim_clicked = st.button(trim_button_label, disabled=trim_button_disabled)

        if trim_clicked:
            progress_bar = st.progress(0, text="Starting trimming...")
            all_success = True
            failed_samples = []

            for i, (_, row) in enumerate(rows_to_process):
                progress_bar.progress(
                    i / n_to_process,
                    text=f"Trimming sample {i + 1} of {n_to_process}: {row['sample']}...",
                )
                success, log = fastp.run_fastp_for_sample(row, trimmed_dir, fastp_reports_dir, threads=fastp_threads)
                if not success:
                    all_success = False
                    failed_samples.append((row["sample"], log))

            progress_bar.progress(1.0, text="Trimming complete.")

            if failed_samples:
                st.error(f"⚠️ fastp failed on {len(failed_samples)} sample(s):")
                for sample_name, log in failed_samples:
                    with st.expander(f"Error details for {sample_name}"):
                        st.code(log)

            total_now_trimmed = len(already_trimmed_rows) + (n_to_process - len(failed_samples)) if not force_retrim_all else (n_to_process - len(failed_samples))

            if all_success:
                st.success(f"✅ Successfully trimmed {n_to_process} sample(s). ({total_now_trimmed} of {len(samplesheet_df)} total now complete.)")
                pm.mark_step_complete(project, "trimming_complete")
                trimming_done = True
            elif len(failed_samples) < n_to_process:
                st.warning(
                    f"Trimmed {n_to_process - len(failed_samples)} of "
                    f"{n_to_process} attempted sample(s) successfully "
                    f"({total_now_trimmed} of {len(samplesheet_df)} total "
                    "now complete). Review the errors above for the rest — "
                    "a common cause is running out of disk space mid-run; "
                    "free up space and click the button again to retry "
                    "just the remaining sample(s)."
                )
                pm.mark_step_complete(project, "trimming_complete")
                trimming_done = True

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 2: Post-Trimming Quality Control (via fastp's own reports)
    # -----------------------------------------------------------------
    st.header("Step 2: Post-Trimming Quality Control")

    with st.expander("ℹ️ What am I looking at here? (click to learn more)"):
        st.markdown(
            "Since fastp already generates a quality report for each "
            "sample while trimming, this step simply combines all of "
            "those individual reports into **one unified MultiQC report** "
            "— the same style of combined report you saw in Step 4 of "
            "the Bulk RNA-Seq workspace, just built from fastp's data "
            "instead of a separate FastQC run.\n\n"
            "The summary table below shows, for each sample: how many "
            "reads were kept after trimming, how many had adapters "
            "removed, and whether overall quality (Q30 rate) improved."
        )

    if not trimming_done:
        st.info("Complete Step 1 (trim your samples) to see post-trimming quality results here.")
        return

    if not multiqc_ok:
        st.error(
            "⚠️ MultiQC was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
        return

    multiqc_html_path = os.path.join(posttrim_multiqc_dir, "multiqc_report.html")
    qc_report_exists = os.path.exists(multiqc_html_path)

    qc_button_label = "🔄 Re-generate Combined Report" if qc_report_exists else "📊 Generate Combined Post-Trim Report"
    qc_clicked = st.button(qc_button_label)

    if qc_report_exists and not qc_clicked:
        st.success("✅ A combined post-trim quality report already exists for this project.")

    if qc_clicked:
        with st.spinner("Combining fastp reports with MultiQC..."):
            multiqc_success, multiqc_log = fastp.run_multiqc(fastp_reports_dir, posttrim_multiqc_dir)

        if not multiqc_success:
            st.error("MultiQC failed to run. Details below:")
            st.code(multiqc_log)
            qc_report_exists = False
        else:
            st.success("✅ Combined post-trim quality report generated.")
            qc_report_exists = True

    if qc_report_exists:
        summary_df = fastp.parse_fastp_reports(fastp_reports_dir)
        if not summary_df.empty:
            st.subheader("📋 Trimming Summary")

            # Show the main summary table without the raw adapter DNA
            # sequences / debug JSON — those are long, technical strings
            # that would clutter the table. The friendly "Adapter
            # Detected" column covers the everyday question of "what was
            # trimmed", and raw details are available below for anyone
            # who wants to verify or troubleshoot.
            display_cols = [
                c for c in summary_df.columns
                if c not in ("Adapter Sequence (R1)", "Adapter Sequence (R2)", "_raw_adapter_cutting_json")
            ]
            st.dataframe(summary_df[display_cols], use_container_width=True, hide_index=True)
            st.caption(
                "**% Reads Kept** close to 100% means little needed to be "
                "removed. A **Q30 After** rate higher than **Q30 Before** "
                "confirms trimming improved overall quality. **Adapter "
                "Detected** shows which library prep kit's adapter fastp "
                "found and removed, if it was able to identify one — for "
                "paired-end data, fastp often trims adapters using "
                "overlap analysis without needing to report a specific "
                "sequence, which is completely normal."
            )

            with st.expander("🔬 View raw adapter details per sample (for troubleshooting)"):
                st.markdown(
                    "These are the exact DNA sequences and raw JSON data "
                    "fastp reported for adapter detection. Most users "
                    "won't need this — it's here for verification or to "
                    "share with support if something looks wrong."
                )
                st.dataframe(
                    summary_df[["Sample", "Adapter Sequence (R1)", "Adapter Sequence (R2)"]],
                    use_container_width=True, hide_index=True,
                )
                for _, row in summary_df.iterrows():
                    st.markdown(f"**{row['Sample']}** — raw `adapter_cutting` JSON:")
                    st.code(row["_raw_adapter_cutting_json"], language="json")

            # --- Residual poly-G / poly-A tail check ---
            # See the "Residual poly-G / poly-A tail detection" section
            # near the top of this file for the full rationale: fastp
            # only auto-trims poly-G on NextSeq/NovaSeq-recognized reads,
            # and never auto-trims poly-A -- so this specific issue can
            # easily survive Step 1 untouched for other platforms/library
            # types. Checked here, using data fastp already wrote to its
            # own JSON report, rather than requiring a separate FastQC
            # pass on the trimmed reads.
            st.subheader("🧬 Poly-G / Poly-A Tail Check")
            with st.expander("ℹ️ What is this, and why does it need a separate check? (click to learn more)"):
                st.markdown(
                    "Some sequencing platforms (mainly Illumina's "
                    "2-color chemistry instruments, like NextSeq/"
                    "NovaSeq) can produce reads with a run of **G**'s at "
                    "the end when there's no real signal left to call. "
                    "Separately, mRNA-seq libraries can sometimes read "
                    "through into the transcript's actual **poly-A "
                    "tail**, producing a run of **A**'s instead. Both "
                    "are normal, well-understood artifacts -- but fastp "
                    "only automatically trims poly-G, and only when it "
                    "recognizes your read headers as coming from a "
                    "NextSeq/NovaSeq-style instrument. Poly-A is **never** "
                    "trimmed automatically. This means either issue can "
                    "slip through Step 1 untouched depending on your "
                    "platform and library type.\n\n"
                    "This check looks directly at fastp's own per-cycle "
                    "base-composition data (already generated during "
                    "Step 1 -- no extra tool needed) for an abnormal "
                    "concentration of G's or A's specifically in the "
                    "last several cycles of your trimmed reads, which is "
                    "the signature of a tail that survived trimming."
                )

            poly_tail_issues = fastp.scan_all_samples_for_poly_tail_issues(fastp_reports_dir)

            if not poly_tail_issues:
                st.success("✅ No residual poly-G or poly-A tails detected in any sample's trimmed reads.")
            else:
                affected_samples = sorted(poly_tail_issues.keys())
                st.warning(
                    f"⚠️ **{len(affected_samples)} of {len(summary_df)} sample(s)** "
                    "still show a residual poly-G and/or poly-A tail after "
                    "trimming:"
                )
                for sample_name in affected_samples:
                    for issue in poly_tail_issues[sample_name]:
                        before_note = (
                            f" (was {issue['tail_pct_before']}% in the raw reads)"
                            if issue["tail_pct_before"] is not None else ""
                        )
                        st.caption(
                            f"- **{sample_name}** ({issue['read']}): "
                            f"**{issue['base']}** makes up **{issue['tail_pct_after']}%** "
                            f"of the last {fastp.POLY_TAIL_WINDOW} cycles, vs. "
                            f"{issue['body_pct_after']}% for the rest of the "
                            f"read{before_note}."
                        )

                st.markdown(
                    "This can be fixed by re-running fastp on just these "
                    "sample(s) with the matching poly-tail trimming option "
                    "explicitly enabled (`--trim_poly_g` for poly-G, "
                    "`--trim_poly_x` for poly-A) -- everything else about "
                    "the trim (adapter detection, quality filtering) stays "
                    "the same."
                )

                if st.button(
                    f"🔁 Auto Re-trim {len(affected_samples)} Flagged Sample(s)",
                    key="poly_tail_retrim_btn",
                ):
                    sample_row_lookup = {
                        row["sample"]: row for _, row in samplesheet_df.iterrows()
                    }
                    retrim_progress = st.progress(0, text="Starting poly-tail re-trim...")
                    retrim_failures = []

                    for i, sample_name in enumerate(affected_samples):
                        retrim_progress.progress(
                            i / len(affected_samples),
                            text=f"Re-trimming {sample_name} ({i + 1} of {len(affected_samples)})...",
                        )
                        bases_flagged = {issue["base"] for issue in poly_tail_issues[sample_name]}
                        sample_row = sample_row_lookup.get(sample_name)
                        if sample_row is None:
                            retrim_failures.append((sample_name, "Sample no longer found in the matched sample list."))
                            continue
                        success, log = fastp.run_fastp_for_sample(
                            sample_row, trimmed_dir, fastp_reports_dir, threads=fastp_threads,
                            force_trim_poly_g=("G" in bases_flagged),
                            force_trim_poly_x=("A" in bases_flagged),
                        )
                        if not success:
                            retrim_failures.append((sample_name, log))

                    retrim_progress.progress(1.0, text="Poly-tail re-trim complete.")

                    if retrim_failures:
                        st.error(f"⚠️ Re-trimming failed for {len(retrim_failures)} sample(s):")
                        for sample_name, log in retrim_failures:
                            with st.expander(f"Error details for {sample_name}"):
                                st.code(log)
                    n_succeeded = len(affected_samples) - len(retrim_failures)
                    if n_succeeded:
                        st.success(
                            f"✅ Re-trimmed {n_succeeded} of {len(affected_samples)} "
                            "flagged sample(s) with poly-tail removal enabled."
                        )
                        # Automatically regenerate the combined MultiQC
                        # report too, so the summary table/report above
                        # reflects the fix immediately on the next
                        # rerun, rather than requiring a separate manual
                        # click of the "Re-generate Combined Report"
                        # button above just to see the result.
                        with st.spinner("Refreshing the combined report..."):
                            fastp.run_multiqc(fastp_reports_dir, posttrim_multiqc_dir)
                        st.rerun()

        if os.path.exists(multiqc_html_path):
            st.subheader("📊 Full Combined Report")
            with open(multiqc_html_path, "rb") as f:
                report_bytes = f.read()
            st.download_button(
                "⬇️ Download Full MultiQC Report (.html)",
                data=report_bytes,
                file_name="multiqc_posttrim_report.html",
                mime="text/html",
            )
            st.caption(
                "Download the file above, then open it in your web "
                "browser to view the full interactive report with "
                "detailed charts for every sample."
            )

        st.markdown("---")
        st.success(
            f"🎉 Project `{project}` now has trimmed reads and post-"
            "trimming quality results saved. This project is ready "
            "for the next step: RNA alignment and gene counting."
        )

        if st.button("➡️ Proceed to RNA Alignment & Counts", type="primary", key="trim_proceed_align_btn"):
            # Same nav_request indirection used by the other "Proceed
            # to X" buttons in this app (see bulk_rnaseq_workspace.py
            # and app.py's module docstring for why a plain session
            # key is used here instead of directly setting
            # st.session_state["assay_choice_radio"]).
            st.session_state["nav_request"] = "🧮 RNA Alignment & Counts"
            st.rerun()
