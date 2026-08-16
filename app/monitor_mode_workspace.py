"""
monitor_mode_workspace.py

Streamlit UI for Monitor Mode -- the second workspace living in the
top-level "Advanced Modes" sidebar drawer alongside "Auto" (see
app.py's PIPELINE_GROUPS). Where Auto launches ONE non-interactive
Bulk RNA-Seq run on demand, Monitor Mode configures one or MORE
long-lived background watchers, each automatically launching a run for
every new sample folder dropped into its own chosen directory.

--- Multiple, independent monitors ---
Each monitor is configured and run almost like its own mini "project":
its own watch directory, its own sample-ID column convention, its own
pipeline preset (genome/alignment method/thread counts/STAR options/
poly-tail handling), and its own notification settings -- started,
stopped, and reconfigured completely independently of every other
monitor (see monitor_manager.py's module docstring for the full
rationale). This page renders every existing monitor as its own
expander (status, controls, activity feed, launched-projects table),
followed by a form to configure and start a brand-new one.

Every project a monitor launches is a completely normal project
afterward -- browsable in Trimming/Alignment/DESeq2/Ontology like any
other -- since Monitor Mode launches through the exact same
advanced_mode_orchestrator engine Auto uses.
"""
import streamlit as st

import project_manager as pm
import reference_manager as ref
import file_browser as fb
import monitor_manager as mm
import notification_manager as notif
import advanced_mode_orchestrator as orch

STATUS_ICONS = {"watching": "👀", "rejected": "⚠️", "launched": "🔄", "complete": "✅", "pipeline_error": "❌"}
ACTIVITY_ICONS = {
    "first_seen": "👁️", "still_copying": "⏳", "rejected": "⚠️", "warning": "⚠️",
    "launched": "🚀", "pipeline_complete": "✅", "pipeline_failed": "❌",
}


# ---------------------------------------------------------------------------
# Shared sub-forms (parameterized by key_prefix so the SAME widgets can be
# rendered independently for multiple monitors without session_state key
# collisions -- see app.py's own note on why grouped-widget state needs
# distinct keys across otherwise-identical screens)
# ---------------------------------------------------------------------------
def _render_preset_options_step(key_prefix, existing=None):
    """
    The reference genome + alignment/quantification options for a
    monitor's pipeline preset -- configured ONCE per monitor and
    applied to every sample folder that monitor finds (individual
    folders don't get their own genome/options screen, since there's
    no user present to fill one out when a folder is dropped
    automatically).
    Returns a dict matching the shape advanced_mode_orchestrator.py's
    config expects for "reference"/"alignment_method"/"threads"/
    "star_options", plus fixed "salmon_count_type"/"use_tximport"
    defaults and the poly-tail toggle.
    """
    existing = existing or {}
    st.subheader("Pipeline Preset")
    is_custom = st.checkbox(
        "Use a custom (non-preset) reference genome",
        value=existing.get("reference", {}).get("is_custom", False), key=f"{key_prefix}_ref_is_custom",
    )
    reference_cfg = {"is_custom": is_custom}
    if not is_custom:
        species_options = list(ref.REFERENCE_CATALOG.keys())
        species_labels = {k: v["label"] for k, v in ref.REFERENCE_CATALOG.items()}
        default_species = existing.get("reference", {}).get("species_key") or species_options[0]
        species_choice = st.selectbox(
            "Reference organism:", options=species_options,
            index=species_options.index(default_species) if default_species in species_options else 0,
            format_func=lambda k: species_labels[k], key=f"{key_prefix}_species_choice",
        )
        reference_cfg["species_key"] = species_choice
    else:
        st.caption("Provide a genome FASTA and GTF already on this server (used for STAR, and for Salmon unless a pre-extracted transcript FASTA is given).")
        genome_fasta = fb.render_server_file_browser(
            key_prefix=f"{key_prefix}_ref_genome_browse", file_extensions=[".fa", ".fasta", ".fa.gz", ".fna"],
            label="Browse for the genome FASTA:",
        )
        gtf_path = fb.render_server_file_browser(
            key_prefix=f"{key_prefix}_ref_gtf_browse", file_extensions=[".gtf", ".gtf.gz"],
            label="Browse for the GTF annotation:",
        )
        reference_cfg.update({
            "species_key": None, "custom_genome_fasta": genome_fasta,
            "custom_gtf": gtf_path, "custom_transcript_fasta": None,
        })
    alignment_options = ["salmon", "star"]
    default_alignment = existing.get("alignment_method", "salmon")
    alignment_method = st.radio(
        "Alignment/quantification method:", alignment_options,
        index=alignment_options.index(default_alignment) if default_alignment in alignment_options else 0,
        format_func=lambda m: "Salmon (pseudo-alignment, faster)" if m == "salmon" else "STAR (splice-aware alignment)",
        key=f"{key_prefix}_alignment_method", horizontal=True,
    )
    star_options_cfg = {"use_two_pass": False, "use_encode_options": False, "add_strand_field": False}
    if alignment_method == "star":
        with st.expander("⚙️ Advanced STAR options (optional)"):
            star_options_cfg["use_two_pass"] = st.checkbox(
                "Two-pass mode (better novel splice junction detection)",
                value=existing.get("star_options", {}).get("use_two_pass", False), key=f"{key_prefix}_star_two_pass_checkbox",
                help=(
                    "STAR aligns your reads once to discover splice junctions, "
                    "then re-aligns using those junctions as extra annotation -- "
                    "roughly doubles alignment time per sample; most useful for "
                    "less-annotated/non-model organisms."
                ),
            )
            star_options_cfg["use_encode_options"] = st.checkbox(
                "Use ENCODE-recommended filtering options",
                value=existing.get("star_options", {}).get("use_encode_options", False), key=f"{key_prefix}_star_encode_options_checkbox",
                help="Applies the standard set of alignment filtering settings used in the ENCODE project's own RNA-seq pipelines.",
            )
            star_options_cfg["add_strand_field"] = st.checkbox(
                "Add strand field to BAM output",
                value=existing.get("star_options", {}).get("add_strand_field", False), key=f"{key_prefix}_star_strand_field_checkbox",
                help="Only needed for downstream tools that expect this strand tag (e.g. Cufflinks-family tools).",
            )
    detected_cores, recommended_threads = pm.get_recommended_thread_count()
    threads = st.slider(
        "Threads for each tool (FastQC/fastp/Salmon/STAR):",
        min_value=1, max_value=detected_cores,
        value=existing.get("threads", {}).get("fastqc", recommended_threads), key=f"{key_prefix}_pipeline_threads",
    )
    thread_cfg = {k: threads for k in ("fastqc", "fastp", "salmon_index", "salmon_quant", "star_index", "star_align")}
    auto_fix_poly_tails = st.toggle(
        "Automatically detect and fix residual poly-G/poly-A tails",
        value=existing.get("auto_fix_poly_tails", True), key=f"{key_prefix}_auto_fix_poly_tails",
    )
    if auto_fix_poly_tails:
        st.caption("✅ **Recommended (on by default).** fastp does not trim poly-G/poly-A tails by default -- leftover tails are automatically detected and re-trimmed after the main fastp pass.")
    else:
        st.caption("⚠️ Poly-G/poly-A tails will be left as fastp's default settings produce them, with no automatic detection or fix.")
    return {
        "reference": reference_cfg,
        "alignment_method": alignment_method,
        "threads": thread_cfg,
        "salmon_count_type": "NumReads",
        "use_tximport": True,
        "star_options": star_options_cfg,
        "auto_fix_poly_tails": auto_fix_poly_tails,
    }


def _render_notification_settings(key_prefix, existing=None):
    """
    Render this monitor's "notify me by email/text" settings -- a fast
    ✅/⚠️ SMTP-readiness check up front (so a user immediately knows
    whether this deployment even supports it), then optional email +
    phone/carrier inputs, three "notify me when..." toggles, and a
    "send test notification" button so delivery can be confirmed
    before relying on it for an unattended run.
    Returns a dict matching monitor_manager.py's "notifications" config
    shape.
    """
    existing = existing or {}
    st.subheader("🔔 Notifications")
    configured, status_message = notif.get_smtp_config_status()
    (st.caption if configured else st.warning)(status_message)
    enabled = st.checkbox(
        "Enable email/text notifications for this monitor",
        value=existing.get("enabled", False), key=f"{key_prefix}_notif_enabled",
    )
    result = {
        "enabled": enabled, "email_address": "", "sms_phone_number": "", "sms_carrier": "",
        "notify_on_warning": True, "notify_on_error": True, "notify_on_success": False,
    }
    if not enabled:
        return result
    col1, col2 = st.columns(2)
    with col1:
        result["email_address"] = st.text_input(
            "Email address (optional):", value=existing.get("email_address", ""), key=f"{key_prefix}_notif_email",
        )
    with col2:
        phone = st.text_input(
            "Phone number for text alerts (optional):", value=existing.get("sms_phone_number", ""),
            key=f"{key_prefix}_notif_phone", placeholder="e.g. 503-555-0100",
        )
        carrier_options = [""] + list(notif.CARRIER_SMS_GATEWAYS.keys())
        existing_carrier = existing.get("sms_carrier", "")
        carrier = st.selectbox(
            "Carrier (for text alerts):", options=carrier_options,
            index=carrier_options.index(existing_carrier) if existing_carrier in carrier_options else 0,
            key=f"{key_prefix}_notif_carrier",
        )
        result["sms_phone_number"] = phone
        result["sms_carrier"] = carrier
        if phone and not carrier:
            st.caption("⚠️ Select a carrier to enable text alerts -- a phone number alone can't be routed to a text.")
    st.caption(
        "ℹ️ Text alerts are sent via your carrier's free email-to-SMS "
        "gateway -- convenient and zero-cost, but not as reliable/"
        "guaranteed as a dedicated SMS API. Some carriers may delay, "
        "truncate, or occasionally drop these; email is the more "
        "dependable channel if that matters for your use case."
    )
    st.markdown("**Notify me when:**")
    c1, c2, c3 = st.columns(3)
    with c1:
        result["notify_on_warning"] = st.checkbox(
            "⚠️ A warning occurs", value=existing.get("notify_on_warning", True), key=f"{key_prefix}_notif_warn",
            help="E.g. a dropped folder had more FASTQ samples than its metadata declared.",
        )
    with c2:
        result["notify_on_error"] = st.checkbox(
            "❌ An error occurs", value=existing.get("notify_on_error", True), key=f"{key_prefix}_notif_error",
            help="A folder was rejected outright, or a launched pipeline run failed partway through.",
        )
    with c3:
        result["notify_on_success"] = st.checkbox(
            "✅ A run completes", value=existing.get("notify_on_success", False), key=f"{key_prefix}_notif_success",
            help="A launched pipeline run finished successfully (gene counts matrix ready).",
        )
    sms_address = notif.build_sms_gateway_address(phone, carrier) if phone and carrier else None
    if (result["email_address"] or sms_address) and st.button("📨 Send test notification", key=f"{key_prefix}_notif_test_btn"):
        outcomes = notif.send_test_notification(email_address=result["email_address"] or None, sms_address=sms_address)
        for channel, outcome in outcomes.items():
            if outcome is None:
                continue
            success, message = outcome
            (st.success if success else st.error)(f"{channel}: {message}")
    return result


# ---------------------------------------------------------------------------
# Per-monitor panel (existing monitor)
# ---------------------------------------------------------------------------
def _render_monitor_activity_feed(monitor_id):
    st.markdown("**📋 Activity Feed**")
    entries = mm.read_activity_log(monitor_id, n_most_recent=50)
    if not entries:
        st.caption("No activity yet.")
        return
    for entry in entries:
        icon = ACTIVITY_ICONS.get(entry.get("event"), "•")
        st.markdown(f"{icon} **{entry.get('timestamp', '—')}** — `{entry.get('folder', '—')}`: {entry.get('event', '—')} — {entry.get('detail', '')}")


def _render_monitor_launched_projects(monitor_id):
    st.markdown("**🚀 Launched Projects**")
    registry = mm.read_registry_summary(monitor_id)
    launched = {k: v for k, v in registry.items() if v.get("status") in ("launched", "complete", "pipeline_error")}
    if not launched:
        st.caption("No projects have been launched yet.")
        return
    for folder_name, info in launched.items():
        project_name = info.get("project_name", folder_name)
        pm_status = pm.load_info(project_name)
        completed_steps = pm_status.get("steps_completed", [])
        icon = STATUS_ICONS.get(info.get("status"), "•")
        warning_note = f" ⚠️ {info['warning']}" if info.get("warning") else ""
        st.markdown(
            f"{icon} **{folder_name}** → project `{project_name}` "
            f"(launched {info.get('launched_at', '—')}) — "
            f"{len(completed_steps)} step(s) completed.{warning_note}"
        )


def _render_existing_monitor(monitor_id):
    config = mm.load_monitor_config(monitor_id)
    running = mm.is_monitor_running(monitor_id)
    status_label = "🟢 Running" if running else "⚪ Stopped"
    with st.expander(f"**{monitor_id}** — {status_label} — watching `{config['watch_dir']}`", expanded=False):
        col1, col2, col3 = st.columns([2, 1, 1])
        with col1:
            st.caption(
                f"Sample-ID column: `{config.get('sample_id_column', mm.DEFAULT_SAMPLE_ID_COLUMN)}` · "
                f"Poll every {config.get('poll_interval_seconds', 30)}s"
            )
        with col2:
            if running:
                if st.button("⏹️ Stop", key=f"stop_{monitor_id}"):
                    mm.stop_monitor(monitor_id)
                    st.success(f"Stopped monitor '{monitor_id}'. Any pipeline run(s) it already launched keep running independently.")
                    st.rerun()
            else:
                if st.button("▶️ Start", key=f"start_{monitor_id}"):
                    pid = mm.launch_monitor_daemon(monitor_id)
                    st.success(f"Started monitor '{monitor_id}' (process ID {pid}).")
                    st.rerun()
        with col3:
            if not running:
                if st.button("🗑️ Delete", key=f"delete_{monitor_id}"):
                    st.session_state[f"_confirm_delete_{monitor_id}"] = True
        if st.session_state.get(f"_confirm_delete_{monitor_id}"):
            st.warning(f"⚠️ Permanently delete monitor '{monitor_id}'? This does NOT delete any project(s) it already launched.")
            cc1, cc2 = st.columns(2)
            with cc1:
                if st.button("Yes, delete this monitor", key=f"confirm_delete_{monitor_id}"):
                    mm.delete_monitor(monitor_id)
                    st.session_state.pop(f"_confirm_delete_{monitor_id}", None)
                    st.success(f"Monitor '{monitor_id}' deleted.")
                    st.rerun()
            with cc2:
                if st.button("Cancel", key=f"cancel_delete_{monitor_id}"):
                    st.session_state.pop(f"_confirm_delete_{monitor_id}", None)
                    st.rerun()
        st.markdown("---")
        _render_monitor_activity_feed(monitor_id)
        st.markdown("---")
        _render_monitor_launched_projects(monitor_id)


# ---------------------------------------------------------------------------
# New monitor setup form
# ---------------------------------------------------------------------------
def _render_new_monitor_form():
    st.subheader("➕ Configure a New Monitor")
    st.caption(
        "Each monitor watches its own directory with its own settings -- "
        "almost like its own mini-project -- and can be started, stopped, "
        "or reconfigured completely independently of any other monitor."
    )
    monitor_name_raw = st.text_input("Monitor name:", key="new_mon_name", placeholder="e.g. human-rnaseq-intake")
    monitor_id = mm.sanitize_monitor_id(monitor_name_raw)
    if monitor_name_raw and not monitor_id:
        st.error("Monitor name must contain at least one letter, number, dash, or underscore.")
    elif monitor_id and mm.monitor_exists(monitor_id):
        st.error(f"A monitor named '{monitor_id}' already exists -- choose a different name.")
    watch_dir = fb.render_server_directory_browser(
        key_prefix="new_mon_watch_dir_browse",
        label="Browse for the directory to watch for new sample folders:",
    )
    sample_id_column = st.text_input(
        "Metadata column containing sample IDs:", value=mm.DEFAULT_SAMPLE_ID_COLUMN, key="new_mon_sample_id_column",
        help=(
            "Each dropped folder must contain exactly one .csv/.txt/.xlsx/"
            ".xls file with 'metadata' in its filename. This is the exact "
            "column name in that file holding each sample's ID -- defaults "
            f"to '{mm.DEFAULT_SAMPLE_ID_COLUMN}' if your metadata already "
            "uses that convention."
        ),
    )
    poll_interval = st.slider(
        "How often to check the watched directory (seconds):",
        min_value=10, max_value=300, value=30, key="new_mon_poll_interval",
        help=(
            "A new folder must be seen as unchanged (same total file size "
            "and count) across two consecutive checks before this monitor "
            "will launch a run against it -- guards against triggering on "
            "a folder that's still being copied/synced into place."
        ),
    )
    st.info(
        "**Before a run is launched, each folder is also checked that:** "
        "① it contains a metadata file as described above, and ② the "
        "number of samples detected from FASTQ filenames matches the "
        "number of samples listed in that metadata column. Fewer FASTQ "
        "samples than the metadata declares blocks the folder entirely "
        "(with a notification, if enabled); MORE FASTQ samples than "
        "declared only triggers a warning notification, and the run "
        "still proceeds using the samples that do match."
    )
    st.markdown("---")
    preset = _render_preset_options_step("new_mon")
    st.markdown("---")
    notifications = _render_notification_settings("new_mon")
    st.markdown("---")
    ready = bool(monitor_id) and not mm.monitor_exists(monitor_id) and bool(watch_dir)
    if not ready:
        st.info("Provide a valid, unused monitor name and a watch directory above to continue.")
        return
    if st.button("▶️ Save & Start Monitor", key="new_mon_save_start_btn", type="primary"):
        config = {
            "monitor_id": monitor_id,
            "watch_dir": watch_dir,
            "poll_interval_seconds": poll_interval,
            "sample_id_column": sample_id_column.strip() or mm.DEFAULT_SAMPLE_ID_COLUMN,
            "pipeline_config_template": preset,
            "notifications": notifications,
        }
        mm.create_monitor(monitor_id, config)
        pid = mm.launch_monitor_daemon(monitor_id)
        st.success(f"✅ Monitor '{monitor_id}' started (process ID {pid}), watching `{watch_dir}`. You can close this tab -- it keeps running in the background.")
        st.rerun()


def render():
    st.title("📁 Monitor Mode")
    st.markdown(
        "Watch one or more directories for newly-dropped sample folders "
        "(FASTQ files + a metadata sheet) and **automatically launch a "
        "Bulk RNA-Seq Auto run for each one** -- FASTQ → trimming → "
        "alignment → gene counts matrix. Each watched directory is "
        "configured as its own independent monitor, with its own "
        "pipeline preset and notification settings. DESeq2 contrasts "
        "and Ontology Analysis remain separate, interactive steps you "
        "run afterward for each resulting project, same as Auto."
    )
    st.markdown("---")
    existing_monitors = mm.list_monitors()
    if existing_monitors:
        st.subheader(f"Configured Monitors ({len(existing_monitors)})")
        for monitor_id in existing_monitors:
            _render_existing_monitor(monitor_id)
        st.markdown("---")
    else:
        st.caption("No monitors configured yet -- set one up below.")
    _render_new_monitor_form()
