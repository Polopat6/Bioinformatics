"""
monitor_manager.py

"Monitor Mode": watch one or more directories for newly-dropped sample
folders (each expected to contain FASTQ files + a metadata sheet), and
automatically launch the same Bulk RNA-Seq Auto pipeline against each
one -- FASTQ ingestion -> QC -> trimming -> reference setup ->
quantification -> gene counts matrix -- exactly the same scope and
engine as the "Auto" workspace's Bulk RNA-Seq handler.

--- Multiple, independent monitors (each configured like its own project) ---
This module supports MANY separate monitors at once, each watching its
own directory with its own settings -- a monitor is identified by a
short, user-chosen monitor_id (sanitized the same way
project_manager.py sanitizes project names) and gets its own
subdirectory under MONITOR_ROOT:

    data/monitor/<monitor_id>/
        config.json             (this monitor's settings -- see
                                  MonitorConfig shape below)
        processed_registry.json (this monitor's "already handled" folders)
        activity.jsonl          (this monitor's human-readable activity feed)
        monitor.pid             (this monitor's daemon process ID)
        monitor.log             (this monitor's daemon's stdout/stderr)

Each monitor runs as its OWN detached background daemon process (see
launch_monitor_daemon), so any one monitor can be started, stopped, or
reconfigured completely independently of every other monitor -- e.g. a
"human RNA-seq intake" monitor watching /scratch/incoming_human and a
"mouse RNA-seq intake" monitor watching /scratch/incoming_mouse, each
with their own reference/alignment preset, sample-ID column
convention, and notification settings, running side by side.

--- Reuses the SAME execution engine as Auto, deliberately ---
Every folder any monitor launches is handed off to
project_manager.create_project() + advanced_mode_orchestrator.
launch_background_run() -- the exact same functions Auto's Bulk
RNA-Seq wizard calls. This keeps a SINGLE source of truth for "what
does a non-interactive Bulk RNA-Seq run actually do": a bug fix or new
feature in advanced_mode_orchestrator.py automatically applies to
every monitor (and to Auto) without needing to be implemented twice.
Launched projects are prefixed with their monitor_id
(see _build_project_name) so two different monitors can never
accidentally collide on the same project name even if they happen to
see identically-named dropped folders.

--- Folder "stability" check ---
A folder can appear in the watched directory before every file has
finished being copied/synced into it (e.g. an rsync/SCP transfer still
in progress). Triggering a run against a still-copying folder would
either crash immediately or -- worse -- silently run against
truncated/corrupt FASTQ data. This is guarded against by requiring a
folder's (total size, file count) fingerprint to be UNCHANGED across
two consecutive poll cycles before it's ever considered a candidate
for validation/launch -- see _process_one_folder below.

--- Folder validation: metadata detection + sample count reconciliation ---
On top of the stability check above, a candidate folder must pass two
further checks before a run is launched (see check_dropped_folder):

  1. Metadata file detection: the folder must contain EXACTLY ONE file
     whose extension is .csv/.txt/.xlsx/.xls AND whose filename
     contains "metadata" (case-insensitive) anywhere in it -- e.g.
     "sample_metadata.csv", "Metadata_v2.xlsx". No other naming
     convention is accepted, and zero or multiple matches means the
     folder is rejected as ambiguous.

  2. Sample-ID column + FASTQ/metadata count reconciliation: each
     monitor is configured with a "sample_id_column" (default
     "sample_id"; the user may point this at any exact column name
     their metadata actually uses). The number of DISTINCT samples
     detected from FASTQ filenames (via
     ingestion_manager.validate_sample_pairs, which already correctly
     collapses R1/R2 mate pairs into one sample) is compared against
     the number of distinct values in that metadata column:
       - FEWER FASTQ samples than metadata rows -> REJECTED. At least
         one sample the metadata declares has no FASTQ data at all,
         which almost always indicates an incomplete/still-wrong drop
         rather than something safe to silently proceed with.
       - MORE FASTQ samples than metadata rows -> a WARNING is logged
         and (if configured) a notification is sent, but the run is
         still LAUNCHED -- the extra, undeclared FASTQ sample(s) are
         simply not part of this analysis (the orchestrator's own
         ingest stage still does exact name-based matching and only
         carries forward samples present in both).
       - Equal counts -> proceeds normally, no warning.

--- Persisted "already processed" registry, per monitor ---
Avoids re-triggering a run against a folder a monitor already launched
(or rejected) in an earlier poll cycle -- or in a PREVIOUS run of that
monitor's daemon, e.g. after a restart -- persisted to disk per
monitor, since the whole point of a long-running watcher is to survive
being restarted without "forgetting" what it already handled.

--- Pipeline completion tracking + notifications ---
Because advanced_mode_orchestrator.launch_background_run() hands a
run off to its own fully detached OS process, this module's watch loop
also periodically polls advanced_mode_orchestrator.get_status() for
every folder it has already launched but not yet seen finish, so it
can send a success/failure notification the moment that background run
actually completes or errors -- see _poll_launched_projects. This is
in addition to (not instead of) the folder-level warning/error
notifications sent at launch time.
"""
import json
import os
import subprocess
import sys
import time
from datetime import datetime

import project_manager as pm
import ingestion_manager as ingest
import advanced_mode_orchestrator as orch
import notification_manager as notif

MONITOR_ROOT = "data/monitor"
FASTQ_EXTENSIONS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
METADATA_EXTENSIONS = (".csv", ".txt", ".xlsx", ".xls")
DEFAULT_SAMPLE_ID_COLUMN = "sample_id"

# Registry statuses that mean "never touch this folder again" (skip it
# outright in the folder-scanning loop). "watching" is deliberately NOT
# included here -- a folder in that state is still being checked for
# fingerprint stability every poll cycle. "launched" IS included (the
# folder itself is done being processed), but its underlying pipeline
# run may still be polled for completion separately -- see
# _poll_launched_projects, which iterates the whole registry
# independent of this skip-list.
_TERMINAL_FOLDER_STATUSES = ("launched", "rejected", "complete", "pipeline_error")

# ---------------------------------------------------------------------------
# Path helpers -- one subdirectory (and one full set of state files) per
# monitor, mirroring project_manager.py's one-function-per-path convention.
# ---------------------------------------------------------------------------
def sanitize_monitor_id(raw_name):
    """Same sanitization rule as project_manager's project-name field -- letters, numbers, dashes, underscores only."""
    return "".join(c for c in str(raw_name).strip() if c.isalnum() or c in ("-", "_"))


def monitor_dir(monitor_id):
    return os.path.join(MONITOR_ROOT, monitor_id)


def monitor_config_path(monitor_id):
    return os.path.join(monitor_dir(monitor_id), "config.json")


def monitor_registry_path(monitor_id):
    """Persisted {folder_name: {...}} of every folder this monitor has already handled -- survives daemon restarts."""
    return os.path.join(monitor_dir(monitor_id), "processed_registry.json")


def monitor_pid_path(monitor_id):
    return os.path.join(monitor_dir(monitor_id), "monitor.pid")


def monitor_log_path(monitor_id):
    return os.path.join(monitor_dir(monitor_id), "monitor.log")


def monitor_activity_log_path(monitor_id):
    """Human-readable, append-only log of every folder this monitor has seen and what it did -- shown in its Activity Feed."""
    return os.path.join(monitor_dir(monitor_id), "activity.jsonl")


# ---------------------------------------------------------------------------
# Monitor config (plain dict, JSON-serializable -- needs to survive being
# written to disk and read back by a completely separate OS process).
#
# {
#   "monitor_id": str,
#   "watch_dir": str,
#   "poll_interval_seconds": int,          # default 30
#   "sample_id_column": str,               # default "sample_id" -- the exact
#                                          # metadata column holding sample IDs
#   "pipeline_config_template": dict,      # everything advanced_mode_orchestrator's
#                                          # config needs EXCEPT fastq_source/
#                                          # fastq_source_dir/metadata_path,
#                                          # which get filled in per-folder.
#   "notifications": {
#       "enabled": bool,
#       "email_address": str,              # "" if unused
#       "sms_phone_number": str,           # "" if unused
#       "sms_carrier": str,                # key into notification_manager.CARRIER_SMS_GATEWAYS, "" if unused
#       "notify_on_warning": bool,         # default True
#       "notify_on_error": bool,           # default True
#       "notify_on_success": bool,         # default False
#   },
# }
# ---------------------------------------------------------------------------
def list_monitors():
    """Return a sorted list of every existing monitor_id (any subdirectory of MONITOR_ROOT containing a config.json)."""
    if not os.path.isdir(MONITOR_ROOT):
        return []
    return sorted([
        name for name in os.listdir(MONITOR_ROOT)
        if os.path.isfile(os.path.join(MONITOR_ROOT, name, "config.json"))
    ])


def monitor_exists(monitor_id):
    return os.path.isfile(monitor_config_path(monitor_id))


def create_monitor(monitor_id, config):
    """
    Create a brand-new monitor's on-disk state and save its config.
    Returns True if created, False if a monitor with this ID already
    exists (mirrors project_manager.create_project's return contract).
    """
    if monitor_exists(monitor_id):
        return False
    os.makedirs(monitor_dir(monitor_id), exist_ok=True)
    save_monitor_config(monitor_id, config)
    return True


def save_monitor_config(monitor_id, config):
    os.makedirs(monitor_dir(monitor_id), exist_ok=True)
    with open(monitor_config_path(monitor_id), "w") as f:
        json.dump(config, f, indent=2)


def load_monitor_config(monitor_id):
    path = monitor_config_path(monitor_id)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def delete_monitor(monitor_id):
    """
    Permanently delete a monitor's entire on-disk state (config,
    registry, activity log, PID/log files) -- does NOT touch any
    project(s) it already launched, since those live under
    project_manager's own PROJECTS_ROOT, completely independent of
    MONITOR_ROOT. Caller is responsible for confirming the monitor is
    stopped first (see is_monitor_running) -- deleting a running
    monitor's config out from under it is not guarded against here.
    """
    import shutil
    d = monitor_dir(monitor_id)
    if os.path.isdir(d):
        shutil.rmtree(d)


# ---------------------------------------------------------------------------
# Notification config helper
# ---------------------------------------------------------------------------
def _effective_notification_config(monitor_config):
    """
    Translate a monitor's raw "notifications" config block (phone
    number + carrier label, as entered in the UI) into the plain
    {"enabled", "email_address", "sms_address", ...} shape
    notification_manager.py's dispatch functions expect -- building the
    actual carrier-gateway email address once here rather than
    scattering that translation across every call site.
    Returns None if notifications are disabled for this monitor.
    """
    raw = monitor_config.get("notifications") or {}
    if not raw.get("enabled"):
        return None
    sms_address = None
    if raw.get("sms_phone_number") and raw.get("sms_carrier"):
        sms_address = notif.build_sms_gateway_address(raw["sms_phone_number"], raw["sms_carrier"])
    return {
        "enabled": True,
        "email_address": raw.get("email_address") or None,
        "sms_address": sms_address,
        "notify_on_warning": raw.get("notify_on_warning", True),
        "notify_on_error": raw.get("notify_on_error", True),
        "notify_on_success": raw.get("notify_on_success", False),
    }


# ---------------------------------------------------------------------------
# Registry (persisted "already processed" set), per monitor
# ---------------------------------------------------------------------------
def _load_registry(monitor_id):
    path = monitor_registry_path(monitor_id)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def _save_registry(monitor_id, registry):
    path = monitor_registry_path(monitor_id)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(registry, f, indent=2)
    os.replace(tmp, path)


def _log_activity(monitor_id, entry):
    entry["timestamp"] = datetime.now().isoformat(timespec="seconds")
    path = monitor_activity_log_path(monitor_id)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def read_activity_log(monitor_id, n_most_recent=50):
    """Read back a monitor's most recent N activity entries, for its live feed -- most recent first."""
    path = monitor_activity_log_path(monitor_id)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        lines = f.readlines()
    entries = [json.loads(line) for line in lines[-n_most_recent:]]
    return list(reversed(entries))


def read_registry_summary(monitor_id):
    """Public read-only accessor for the Monitor Mode UI's per-monitor project list -- returns the plain registry dict."""
    return _load_registry(monitor_id)


# ---------------------------------------------------------------------------
# Folder detection, stability, and validation
# ---------------------------------------------------------------------------
def _folder_fingerprint(folder_path):
    """
    A cheap (total size, file count) fingerprint used to detect whether
    a folder's contents are still changing (e.g. an in-progress copy).
    Deliberately lightweight (not a full hash), since this runs on
    every poll cycle for every candidate folder across every monitor.
    """
    total_size = 0
    file_count = 0
    for root, _dirs, files in os.walk(folder_path):
        for fname in files:
            try:
                total_size += os.path.getsize(os.path.join(root, fname))
                file_count += 1
            except OSError:
                continue  # file vanished mid-walk (e.g. a temp file) -- ignore
    return (total_size, file_count)


def _find_metadata_file(folder_path):
    """
    Find this folder's metadata file: exactly one file whose extension
    is .csv/.txt/.xlsx/.xls AND whose filename contains "metadata"
    (case-insensitive) anywhere in it -- e.g. "sample_metadata.csv",
    "Metadata.xlsx", "RNAseq_METADATA_v2.txt" all match.
    Returns the file's path, or None if zero or more than one such
    file is found (an ambiguous folder is rejected, not guessed at).
    """
    candidates = []
    for fname in os.listdir(folder_path):
        full_path = os.path.join(folder_path, fname)
        if not os.path.isfile(full_path):
            continue
        name_no_ext, ext = os.path.splitext(fname)
        if ext.lower() not in METADATA_EXTENSIONS:
            continue
        if "metadata" in fname.lower():
            candidates.append(full_path)
    if len(candidates) == 1:
        return candidates[0]
    return None


def check_dropped_folder(folder_path, sample_id_column):
    """
    Validate that a dropped folder actually looks like a complete,
    usable sample set BEFORE launching a pipeline run against it -- see
    this module's docstring's "Folder validation" section for the full
    rules. Any blocking issue here means "reject this folder and log
    why", never "guess and proceed"; a non-blocking issue (more FASTQ
    samples than the metadata declares) is surfaced as a warning
    instead, and the caller proceeds with the launch.

    Returns a dict:
        {
            "ready": bool,
            "reason": str,               # human-readable status/rejection reason
            "fastq_sample_count": int,
            "metadata_sample_count": int or None,
            "metadata_path": str or None,
            "warning": str or None,      # non-blocking issue, if any
        }
    """
    _, fastq_files = ingest.find_fastq_filenames_in_directory(folder_path)
    if not fastq_files:
        return {
            "ready": False, "reason": "No FASTQ files found in this folder.",
            "fastq_sample_count": 0, "metadata_sample_count": None,
            "metadata_path": None, "warning": None,
        }
    metadata_path = _find_metadata_file(folder_path)
    if not metadata_path:
        return {
            "ready": False,
            "reason": (
                "No metadata file found -- this folder must contain exactly "
                "ONE file with a .csv/.txt/.xlsx/.xls extension whose "
                "filename includes 'metadata' (e.g. 'sample_metadata.csv')."
            ),
            "fastq_sample_count": len(fastq_files), "metadata_sample_count": None,
            "metadata_path": None, "warning": None,
        }
    meta_df, read_error = ingest.read_metadata_file(metadata_path)
    if read_error:
        return {
            "ready": False, "reason": f"Metadata file could not be read: {read_error}",
            "fastq_sample_count": len(fastq_files), "metadata_sample_count": None,
            "metadata_path": metadata_path, "warning": None,
        }
    if sample_id_column not in meta_df.columns:
        return {
            "ready": False,
            "reason": (
                f"Metadata file has no '{sample_id_column}' column "
                f"(configured sample-ID column for this monitor). Found "
                f"columns: {list(meta_df.columns)}."
            ),
            "fastq_sample_count": len(fastq_files), "metadata_sample_count": None,
            "metadata_path": metadata_path, "warning": None,
        }
    sample_pairs = ingest.validate_sample_pairs(fastq_files)
    n_fastq_samples = len(sample_pairs)
    n_metadata_samples = meta_df[sample_id_column].dropna().astype(str).nunique()
    if n_fastq_samples < n_metadata_samples:
        return {
            "ready": False,
            "reason": (
                f"FASTQ/metadata sample count mismatch: only {n_fastq_samples} "
                f"sample(s) detected from FASTQ filenames, but the metadata's "
                f"'{sample_id_column}' column lists {n_metadata_samples} "
                "sample(s). Every sample declared in the metadata needs "
                "matching FASTQ file(s) before this folder can be processed."
            ),
            "fastq_sample_count": n_fastq_samples, "metadata_sample_count": n_metadata_samples,
            "metadata_path": metadata_path, "warning": None,
        }
    warning = None
    if n_fastq_samples > n_metadata_samples:
        extra = n_fastq_samples - n_metadata_samples
        warning = (
            f"{n_fastq_samples} FASTQ sample(s) were detected, but the "
            f"metadata's '{sample_id_column}' column only lists "
            f"{n_metadata_samples} sample(s) -- {extra} extra FASTQ "
            "sample(s) are not declared in the metadata and will be "
            "excluded from the analysis. Proceeding anyway."
        )
    return {
        "ready": True, "reason": "OK",
        "fastq_sample_count": n_fastq_samples, "metadata_sample_count": n_metadata_samples,
        "metadata_path": metadata_path, "warning": warning,
    }


def _build_project_name(monitor_id, folder_name):
    """
    Build this folder's project name, prefixed with its monitor_id so
    two different monitors can never collide on the same project name
    even if they happen to see identically-named dropped folders (this
    also makes it immediately obvious, from the project list alone,
    which monitor launched which project -- "almost like its own
    project" per monitor, as intended).
    """
    raw = f"{monitor_id}_{folder_name}"
    return "".join(c for c in raw.strip() if c.isalnum() or c in ("-", "_"))


def _build_project_config_for_folder(monitor_config, folder_path, project_name):
    """
    Build the advanced_mode_orchestrator config dict for one dropped
    folder, starting from the monitor's preset pipeline_config_template
    and filling in this specific folder's FASTQ source + metadata.
    fastq_source_dir is pointed directly at the dropped folder itself
    (NOT copied into the project first) -- advanced_mode_orchestrator's
    own ingest stage already symlinks from whatever directory it's
    given into the project's own fastq_dir.
    """
    config = dict(monitor_config.get("pipeline_config_template", {}))
    config["fastq_source"] = "directory"
    config["fastq_source_dir"] = folder_path
    config["metadata_path"] = pm.metadata_path(project_name)
    return config


def _process_one_folder(monitor_id, monitor_config, folder_name, registry):
    """
    Handle exactly one candidate folder: check fingerprint stability,
    validate it (metadata detection + sample-ID column + FASTQ/
    metadata count reconciliation), and if it passes, create a project
    and launch it via advanced_mode_orchestrator -- exactly the same
    launch path Auto's own launch button uses. Every outcome (still
    copying, rejected, warned-but-launched, launched) is written to
    both the persisted registry and the human-readable activity log,
    and warnings/errors additionally trigger this monitor's configured
    notification channel(s).
    """
    notification_config = _effective_notification_config(monitor_config)
    folder_path = os.path.join(monitor_config["watch_dir"], folder_name)
    fingerprint_now = _folder_fingerprint(folder_path)
    last_seen = registry.get(folder_name, {}).get("last_fingerprint")
    if last_seen is None:
        # First time seeing this folder at all -- record its
        # fingerprint and wait for the NEXT poll cycle to check
        # stability, rather than risk validating/launching against a
        # still-copying folder.
        registry[folder_name] = {"status": "watching", "last_fingerprint": list(fingerprint_now)}
        _log_activity(monitor_id, {"folder": folder_name, "event": "first_seen", "detail": str(fingerprint_now)})
        return
    if tuple(last_seen) != fingerprint_now:
        registry[folder_name]["last_fingerprint"] = list(fingerprint_now)
        _log_activity(monitor_id, {"folder": folder_name, "event": "still_copying", "detail": str(fingerprint_now)})
        return
    # Fingerprint stable across two consecutive polls -- proceed to validate.
    sample_id_column = monitor_config.get("sample_id_column") or DEFAULT_SAMPLE_ID_COLUMN
    check = check_dropped_folder(folder_path, sample_id_column)
    if not check["ready"]:
        registry[folder_name] = {"status": "rejected", "reason": check["reason"]}
        _log_activity(monitor_id, {"folder": folder_name, "event": "rejected", "detail": check["reason"]})
        notif.notify_error(notification_config, monitor_id, folder_name, check["reason"])
        return
    if check["warning"]:
        _log_activity(monitor_id, {"folder": folder_name, "event": "warning", "detail": check["warning"]})
        notif.notify_warning(notification_config, monitor_id, folder_name, check["warning"])
    project_name = _build_project_name(monitor_id, folder_name)
    if not project_name:
        reason = f"Folder name '{folder_name}' produced an empty/invalid project name after sanitization."
        registry[folder_name] = {"status": "rejected", "reason": reason}
        _log_activity(monitor_id, {"folder": folder_name, "event": "rejected", "detail": reason})
        notif.notify_error(notification_config, monitor_id, folder_name, reason)
        return
    if project_name in pm.list_projects():
        # A project with this name already exists -- most likely this
        # folder was already processed in a previous daemon lifetime
        # and the registry itself is what's stale. Rather than
        # silently launching a duplicate/conflicting run against an
        # existing project, skip and log clearly so a user can
        # investigate.
        reason = f"A project named '{project_name}' already exists -- skipping to avoid conflicting with it."
        registry[folder_name] = {"status": "rejected", "reason": reason}
        _log_activity(monitor_id, {"folder": folder_name, "event": "rejected", "detail": reason})
        notif.notify_error(notification_config, monitor_id, folder_name, reason)
        return
    pm.create_project(project_name)
    # Copy the found metadata file into the project's own metadata_path
    # (rather than pointing advanced_mode_orchestrator directly at the
    # file still sitting in the watched folder), giving the project its
    # own hermetic copy -- consistent with how Auto's own wizard saves
    # the user's reviewed/edited metadata, and safe against the watched
    # folder later being modified or cleaned up externally.
    #
    # The rest of the pipeline (ingestion_manager /
    # advanced_mode_orchestrator) expects the sample-ID column to be
    # named exactly "sample" -- rather than threading a configurable
    # column name through every downstream module, the monitor's
    # chosen sample_id_column is renamed to "sample" HERE, once, in
    # this hermetic copy, so every existing downstream step keeps
    # working completely unmodified.
    meta_df, _ = ingest.read_metadata_file(check["metadata_path"])
    if sample_id_column != "sample":
        meta_df = meta_df.rename(columns={sample_id_column: "sample"})
    meta_df.to_csv(pm.metadata_path(project_name), index=False)
    config = _build_project_config_for_folder(monitor_config, folder_path, project_name)
    pid = orch.launch_background_run(project_name, config)
    registry[folder_name] = {
        "status": "launched", "project_name": project_name,
        "pid": pid, "launched_at": datetime.now().isoformat(timespec="seconds"),
        "notified_complete": False,
        "warning": check["warning"],
    }
    _log_activity(monitor_id, {
        "folder": folder_name, "event": "launched",
        "detail": f"project={project_name}, pid={pid}, {check['fastq_sample_count']} FASTQ sample(s)"
                  + (f" (warning: {check['warning']})" if check["warning"] else ""),
    })


def _poll_launched_projects(monitor_id, monitor_config, registry):
    """
    For every folder this monitor has already launched but whose
    background pipeline run hasn't yet been confirmed finished, check
    advanced_mode_orchestrator's on-disk status for that project and,
    the moment it completes or errors, send the configured success/
    failure notification and record the outcome -- so a user gets
    notified about the ACTUAL multi-hour pipeline result, not just
    "a folder was launched".
    """
    notification_config = _effective_notification_config(monitor_config)
    for folder_name, entry in registry.items():
        if entry.get("status") != "launched" or entry.get("notified_complete"):
            continue
        project_name = entry.get("project_name")
        if not project_name:
            continue
        status = orch.get_status(project_name)
        pipeline_status = status.get("pipeline_status")
        if pipeline_status == "complete":
            notif.notify_pipeline_success(notification_config, monitor_id, project_name)
            entry["status"] = "complete"
            entry["notified_complete"] = True
            _log_activity(monitor_id, {"folder": folder_name, "event": "pipeline_complete", "detail": f"project={project_name}"})
        elif pipeline_status == "error":
            error_detail = status.get("error") or {}
            stage = error_detail.get("stage", "unknown")
            message = error_detail.get("message", "Unknown error.")
            notif.notify_pipeline_error(notification_config, monitor_id, project_name, stage, message)
            entry["status"] = "pipeline_error"
            entry["notified_complete"] = True
            _log_activity(monitor_id, {"folder": folder_name, "event": "pipeline_failed", "detail": f"project={project_name}, stage={stage}: {message}"})
        # else: still "queued"/"running" -- check again next poll cycle.


# ---------------------------------------------------------------------------
# Main watch loop
# ---------------------------------------------------------------------------
def run_monitor_loop(monitor_config):
    """
    The actual infinite loop for ONE monitor -- meant to be the entire
    body of that monitor's daemon process (see monitor_daemon_runner.py
    / launch_monitor_daemon below). Never returns under normal
    operation; stopped externally via SIGTERM (see stop_monitor).
    """
    monitor_id = monitor_config["monitor_id"]
    while True:
        registry = _load_registry(monitor_id)
        watch_dir = monitor_config["watch_dir"]
        if os.path.isdir(watch_dir):
            for entry in sorted(os.listdir(watch_dir)):
                entry_path = os.path.join(watch_dir, entry)
                if not os.path.isdir(entry_path):
                    continue
                if registry.get(entry, {}).get("status") in _TERMINAL_FOLDER_STATUSES:
                    continue  # already handled, terminal state -- never re-process
                _process_one_folder(monitor_id, monitor_config, entry, registry)
        _poll_launched_projects(monitor_id, monitor_config, registry)
        _save_registry(monitor_id, registry)
        time.sleep(monitor_config.get("poll_interval_seconds", 30))


# ---------------------------------------------------------------------------
# Background-launch helpers for the daemon itself (one independent
# detached process PER MONITOR -- separate from any individual pipeline
# run's own process; see advanced_mode_orchestrator.launch_background_run
# for the analogous per-run version and its docstring's rationale, which
# applies equally here).
# ---------------------------------------------------------------------------
def launch_monitor_daemon(monitor_id, python_executable=None):
    """
    Launch monitor_daemon_runner.py as its own, fully independent OS
    process for this specific monitor_id, via subprocess.Popen -- kept
    running even if the Streamlit server restarts or every browser tab
    is closed, exactly like an individual pipeline run's own background
    process. Assumes this monitor's config has already been saved (via
    create_monitor/save_monitor_config).
    Returns the launched process's PID (int).
    """
    python_executable = python_executable or sys.executable
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, "monitor_daemon_runner.py")
    log_path = monitor_log_path(monitor_id)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    log_file = open(log_path, "a")
    log_file.write(f"\n\n===== Launching Monitor '{monitor_id}' daemon at {datetime.now().isoformat(timespec='seconds')} =====\n")
    log_file.flush()
    popen_kwargs = {}
    if os.name == "posix":
        popen_kwargs["start_new_session"] = True
    proc = subprocess.Popen(
        [python_executable, script_path, monitor_config_path(monitor_id)],
        stdout=log_file, stderr=subprocess.STDOUT,
        **popen_kwargs,
    )
    with open(monitor_pid_path(monitor_id), "w") as f:
        f.write(str(proc.pid))
    return proc.pid


def get_monitor_pid(monitor_id):
    path = monitor_pid_path(monitor_id)
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return int(f.read().strip())
    except (ValueError, OSError):
        return None


def is_monitor_running(monitor_id):
    """Zombie-safe liveness check for this monitor's daemon process -- reuses advanced_mode_orchestrator's own zombie-aware is_process_alive rather than duplicating that logic."""
    return orch.is_process_alive(get_monitor_pid(monitor_id))


def stop_monitor(monitor_id):
    """
    Gracefully stop this monitor's watcher (SIGTERM, not SIGKILL -- the
    loop has no in-progress state that needs flushing). Returns True if
    a running daemon was found and signaled, False if nothing was
    running.
    IMPORTANT: this stops the WATCHER itself, not any pipeline run(s)
    it has already launched -- those are independent detached
    processes and continue running to completion even after the
    watcher stops.
    """
    if not is_monitor_running(monitor_id):
        return False
    pid = get_monitor_pid(monitor_id)
    try:
        os.kill(pid, 15)  # SIGTERM
        return True
    except OSError:
        return False


def any_monitor_running():
    """Convenience check used by the UI to decide whether ANY monitor is currently active, across all configured monitors."""
    return any(is_monitor_running(mid) for mid in list_monitors())
