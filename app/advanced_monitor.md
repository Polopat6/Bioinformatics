# Advanced Mode + Monitor Mode — Core Infrastructure

This describes the ACTUAL, currently-implemented design (this file
previously described an earlier, more elaborate plan --
`pipeline_config.py`, `pipeline_orchestrator.py`, `email_notifier.py`,
`background_launcher.py`, and `run_pipeline_headless.py` -- that was
never actually built; delete any references to those files/classes.
The real implementation below is simpler, dict-based throughout (no
`PipelineConfig`/`MonitorConfig` classes), and is what's live in the
app today.

## Files (dependency order)

- **project_manager.py** / **advanced_mode_orchestrator.py** — existing,
  unchanged engine: `create_project()`, step-tracking, and the
  resumable, detached-background-process pipeline runner (FASTQ →
  QC → trimming → reference → quantification → counts matrix, then
  STOPS -- DESeq2/Ontology remain separate interactive steps). Both
  Auto and every Monitor use this exact engine via
  `launch_background_run()` / `get_status()`.
- **notification_manager.py** — plain-`smtplib`-based email
  notifications, plus "text message" delivery via free carrier
  email-to-SMS gateways (`CARRIER_SMS_GATEWAYS`; no paid SMS API
  involved, so no new dependency). SMTP server/credentials are
  deployment-level environment variables (`PIPELINE_SMTP_HOST`,
  `_PORT`, `_USER`, `_PASSWORD`, `_FROM_ADDRESS`, `_USE_TLS`) — NOT
  stored per-monitor. `get_smtp_config_status()` gives a fast ✅/⚠️
  readiness check for the setup UI; `send_test_notification()` backs
  a "send test notification" button. Every public function returns
  `(success, message)` and never raises — a broken/misconfigured relay
  must never crash a pipeline run or the monitor daemon. Distinct
  `notify_warning` / `notify_error` / `notify_pipeline_error` /
  `notify_pipeline_success` functions cover, respectively: a
  non-blocking folder issue (extra undeclared FASTQ samples), a
  blocking folder rejection, a background pipeline run failing
  mid-run, and a background pipeline run finishing successfully.
- **monitor_manager.py** — Monitor Mode's directory watcher, now
  supporting **multiple independent monitors** at once. Each monitor
  is identified by a short `monitor_id` and gets its own state
  directory (`data/monitor/<monitor_id>/`: `config.json`,
  `processed_registry.json`, `activity.jsonl`, `monitor.pid`,
  `monitor.log`) and its own detached daemon process — start/stop/
  reconfigure any one monitor without touching any other. Each
  monitor's config carries its own `watch_dir`,
  `poll_interval_seconds`, `sample_id_column`,
  `pipeline_config_template`, and `notifications` block.
  Per candidate folder, in order:
  1. **Stability check** — folder's (size, file count) fingerprint
     must be unchanged across two consecutive polls (guards against a
     still-copying/still-syncing folder).
  2. **Metadata detection** — the folder must contain *exactly one*
     file with a `.csv`/`.txt`/`.xlsx`/`.xls` extension whose filename
     contains `"metadata"` (case-insensitive) anywhere in it. Zero or
     multiple matches → rejected as ambiguous.
  3. **Sample-ID column + count reconciliation** — the monitor's
     configured `sample_id_column` (default `"sample_id"`) must exist
     in that metadata file. The number of distinct FASTQ-derived
     samples (via `ingestion_manager.validate_sample_pairs`, which
     already collapses R1/R2 pairs) is compared against the number of
     distinct values in that column:
     - **fewer** FASTQ samples than metadata rows → **rejected**
       (blocking) — some declared sample has no FASTQ data at all.
     - **more** FASTQ samples than metadata rows → **warning only**,
       run still launches (extra undeclared samples are simply
       excluded downstream by the orchestrator's own exact-name
       matching).
     - equal → proceeds silently.
  Launched projects are named `<monitor_id>_<folder_name>` so two
  monitors can never collide on the same project name. After launch,
  the watch loop also polls `advanced_mode_orchestrator.get_status()`
  for every not-yet-finished launched project every cycle, firing a
  success/failure notification the moment that background run actually
  completes or errors (`_poll_launched_projects`).
- **monitor_daemon_runner.py** — tiny, Streamlit-free CLI entry point
  for ONE monitor's daemon: `python monitor_daemon_runner.py
  <that_monitor's_config.json_path>`. Launched by
  `monitor_manager.launch_monitor_daemon(monitor_id)`.
- **monitor_mode_workspace.py** — Streamlit UI. Lists every configured
  monitor as its own expander (status, Start/Stop, Delete, activity
  feed, launched-projects table), followed by a "Configure a New
  Monitor" form (watch directory, sample-ID column name, poll
  interval, pipeline preset, and notification setup — email address
  and/or phone+carrier, three `notify_on_*` toggles, and a "send test
  notification" button using `notification_manager`).

## Known stale/orphaned files — recommended cleanup

- **`run_monitor_headless.py`** references `from monitor_manager import
  MonitorConfig, run_monitor_loop` — `MonitorConfig` does not exist
  anywhere in the real (dict-based) implementation, and this file is
  never imported by `app.py` or `monitor_manager.py`. It's a leftover
  from the earlier, unbuilt class-based design described above.
  **Recommend deleting this file** — `monitor_daemon_runner.py` is the
  real, working daemon entry point, and now takes a per-monitor config
  path (see above).

## Required deployment setup (for notifications)

Set these environment variables wherever the app/monitor daemons run
(Docker env, HPC job script, local `.env`) — without them, email/text
notifications silently no-op (logged, never a hard error):

```
PIPELINE_SMTP_HOST=smtp.yourorg.example.com
PIPELINE_SMTP_PORT=587                  # optional, defaults to 587
PIPELINE_SMTP_USER=your_smtp_username    # optional, if auth required
PIPELINE_SMTP_PASSWORD=your_smtp_password
PIPELINE_SMTP_FROM_ADDRESS=portal-noreply@yourorg.example.com
PIPELINE_SMTP_USE_TLS=true               # optional, defaults to true
```

"Text message" notifications reuse this same SMTP relay via free
carrier email-to-SMS gateways (see `notification_manager.
CARRIER_SMS_GATEWAYS`) — no separate SMS account/API key needed, but
also not as guaranteed-delivery as a real SMS API; surfaced as a
caveat directly in the setup UI.

## What's still worth doing next

- Consider extending per-project notification hooks (success/failure
  emails) to **Auto** runs too, not just Monitor-launched ones — Auto
  currently has no notification config at all since it's normally
  watched interactively, but a long unattended Auto run (e.g. a large
  NCBI/SRA fetch) could benefit from the same failure alert.
- `sample_id_column` mismatches are currently caught only at the
  per-folder validation step; consider surfacing the SAME check as a
  live preview in the "Configure a New Monitor" form (e.g. "browse a
  sample folder now to test your settings") before a monitor is ever
  started, so a misconfigured column name is caught immediately rather
  than on the first real dropped folder.
- Carrier email-to-SMS gateway domains occasionally change or get
  discontinued without notice — worth a periodic sanity check against
  `notification_manager.CARRIER_SMS_GATEWAYS`.
