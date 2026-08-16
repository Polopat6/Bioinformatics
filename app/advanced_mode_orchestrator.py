"""
advanced_mode_orchestrator.py

Non-interactive pipeline runner for Advanced Mode / Monitor Mode.

This is the piece that ties together every Streamlit-independent
manager module (ingestion_manager, fastqc_manager, fastp_manager,
reference_manager, quantification_manager, counts_matrix_manager) into
a single, resumable, background-runnable pipeline:

    FASTQ ingestion -> pre-trim QC -> trimming -> reference setup ->
    quantification -> gene counts matrix

...and STOPS THERE. DESeq2 contrasts and Ontology Analysis are
explicitly OUT OF SCOPE for this orchestrator -- they remain
interactive steps a user performs afterward in the normal step-by-step
workspaces, once this pipeline has produced a finished counts matrix.

--------------------------------------------------------------------
WHY THIS MODULE EXISTS / HOW IT'S MEANT TO BE RUN
--------------------------------------------------------------------
A full run of this pipeline (FASTQ -> counts matrix) can take many
hours for real datasets (reference download, index building, and
alignment/quantification dominate). It is NOT meant to run inside a
Streamlit request/response cycle -- Streamlit's own process model isn't
built for multi-hour background work, and a user closing their browser
tab shouldn't kill the run.

Instead, this module is designed to be launched as its OWN OS process,
independent of the Streamlit app's process, e.g.:

    subprocess.Popen(
        ["python", "advanced_mode_orchestrator.py", project_name, config_json_path],
        stdout=open(log_path, "a"), stderr=subprocess.STDOUT,
    )

The Streamlit app then POLLS the pipeline's on-disk JSON status file
(see status_path() below) to render live progress, rather than holding
a connection open to the running process itself. This also gives free
resumability: if the orchestrator process is killed (server restart,
crash, etc.) partway through, simply launching it again for the same
project safely picks up from the last completed step (see "RESUMABILITY"
below) instead of starting over.

--------------------------------------------------------------------
RESUMABILITY
--------------------------------------------------------------------
Before running each stage, the orchestrator checks
project_manager.has_completed_step() for that stage's step name and
skips it entirely if already marked complete -- the exact same
step-tracking mechanism the interactive workspaces already use, so an
Advanced Mode run and a manually-run interactive session are always
mutually consistent (a project doesn't care WHICH path completed a
step, only that it's done). Within a stage, per-sample work (trimming,
quantification) additionally uses each manager module's own
already-done check (e.g. fastp_manager.sample_already_trimmed) so a
partially-completed stage resumes at the sample level too, not just
the pipeline-stage level.

--------------------------------------------------------------------
STATUS FILE
--------------------------------------------------------------------
Written to <project_dir>/advanced_mode_status.json throughout the run.
Shape:
    {
      "pipeline_status": "running" | "complete" | "error",
      "current_stage": "trimming" | ... | None,
      "started_at": "...", "updated_at": "...",
      "error": null or {"stage": ..., "message": ...},
      "stages": {
        "ingest": {"status": "pending"|"running"|"complete"|"skipped"|"error",
                   "message": "...", "started_at": ..., "finished_at": ...},
        ... one entry per PIPELINE_STAGES below ...
      }
    }
"""

import json
import os
import subprocess
import sys
import traceback
from datetime import datetime

import pandas as pd

import project_manager as pm
import reference_manager as rm
import quantification_manager as qm
import counts_matrix_manager as cmm
import ingestion_manager as ingest
import fastqc_manager as fastqc
import fastp_manager as fastp
import sra_manager as sra

# Ordered list of (stage_key, project_manager step_name) -- this order
# IS the pipeline's execution order, and matches the step-tracking
# vocabulary alignment_workspace.py / bulk_rnaseq_workspace.py /
# trimming_workspace.py already use, so Advanced Mode and manual/
# interactive progress are always readable from the exact same
# project_info.json.
PIPELINE_STAGES = [
    ("ingest", "samples_matched"),
    ("qc", "qc_complete"),
    ("trimming", "trimming_complete"),
    ("reference", "reference_ready"),
    ("quantification", "quantification_complete"),
    ("counts_matrix", "counts_matrix_complete"),
]


# ---------------------------------------------------------------------------
# Status file helpers
# ---------------------------------------------------------------------------

def status_path(project_name):
    return os.path.join(pm.project_dir(project_name), "advanced_mode_status.json")


def _now():
    return datetime.now().isoformat(timespec="seconds")


def _load_status(project_name):
    path = status_path(project_name)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {
        "pipeline_status": "not_started",
        "current_stage": None,
        "started_at": None,
        "updated_at": None,
        "error": None,
        "stages": {
            key: {"status": "pending", "message": None, "started_at": None, "finished_at": None}
            for key, _ in PIPELINE_STAGES
        },
    }


def _save_status(project_name, status):
    status["updated_at"] = _now()
    os.makedirs(pm.project_dir(project_name), exist_ok=True)
    with open(status_path(project_name), "w") as f:
        json.dump(status, f, indent=2)


def _set_stage_status(project_name, status, stage_key, stage_status, message=None):
    """Update one stage's entry in-place and persist the whole status file."""
    entry = status["stages"][stage_key]
    entry["status"] = stage_status
    if message is not None:
        entry["message"] = message
    if stage_status == "running" and entry["started_at"] is None:
        entry["started_at"] = _now()
    if stage_status in ("complete", "skipped", "error"):
        entry["finished_at"] = _now()
    status["current_stage"] = stage_key if stage_status == "running" else status["current_stage"]
    _save_status(project_name, status)


def get_status(project_name):
    """Public read-only accessor for a Streamlit polling UI."""
    return _load_status(project_name)


# ---------------------------------------------------------------------------
# Background-launch helpers -- for a Streamlit UI to kick off a run as an
# independent OS process (see module docstring's "HOW THIS IS MEANT TO BE
# RUN" section) and later check whether it's still going, without ever
# holding a direct handle to the running process across Streamlit reruns
# (which would be lost the moment the script reruns anyway).
# ---------------------------------------------------------------------------

def config_path(project_name):
    return os.path.join(pm.project_dir(project_name), "advanced_mode_config.json")


def log_path(project_name):
    return os.path.join(pm.project_dir(project_name), "advanced_mode_log.txt")


def launch_info_path(project_name):
    return os.path.join(pm.project_dir(project_name), "advanced_mode_launch_info.json")


def save_config(project_name, config):
    os.makedirs(pm.project_dir(project_name), exist_ok=True)
    with open(config_path(project_name), "w") as f:
        json.dump(config, f, indent=2)


def load_config(project_name):
    path = config_path(project_name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _save_launch_info(project_name, pid):
    info = {"pid": pid, "launched_at": _now()}
    with open(launch_info_path(project_name), "w") as f:
        json.dump(info, f, indent=2)


def get_launch_info(project_name):
    path = launch_info_path(project_name)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _is_zombie(pid):
    """
    Check whether pid is currently a ZOMBIE process (already exited,
    but not yet reaped by a parent calling wait()/poll()) -- this
    matters because os.kill(pid, 0) below returns successfully for a
    zombie too, not just a genuinely still-running process. Since every
    Streamlit rerun is a fresh top-level script execution with no
    memory of the original subprocess.Popen object (and Monitor Mode's
    watch loop launches many runs in sequence over a potentially very
    long lifetime, making this edge case meaningfully more likely to
    actually occur than a single manual Advanced Mode launch), nothing
    in this process naturally calls .wait()/.poll() on these detached
    children to reap them -- so, left unhandled, a run that finished
    hours ago could be reported as "still running" indefinitely.

    Tries /proc first (Linux -- covers Docker and HPC, this app's two
    primary deployment targets), falling back to `ps` (e.g. macOS for
    local development) if /proc isn't available.
    """
    proc_status_path = f"/proc/{pid}/status"
    if os.path.exists(proc_status_path):
        try:
            with open(proc_status_path) as f:
                for line in f:
                    if line.startswith("State:"):
                        return "zombie" in line.lower() or line.strip().split()[1] == "Z"
        except OSError:
            return False
        return False

    try:
        result = subprocess.run(
            ["ps", "-o", "stat=", "-p", str(pid)],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip().startswith("Z")
    except (OSError, subprocess.TimeoutExpired):
        return False


def is_process_alive(pid):
    """
    Best-effort check for whether a PID is still a running process on
    this machine, using only the standard library (no psutil
    dependency) -- os.kill(pid, 0) sends no actual signal, it just
    probes whether the OS will let us signal that PID, which fails
    with ProcessLookupError/OSError if it doesn't exist. This can
    still be fooled by PID reuse (the original process exited and a
    new, unrelated process later got assigned the same PID), which is
    why is_run_in_progress() below layers the pipeline_status file on
    top of this rather than trusting this check alone.

    Also treats a ZOMBIE process (see _is_zombie above) as NOT alive --
    from this function's caller's perspective ("has this pipeline run
    actually finished?"), a zombie has already finished its real work;
    it just hasn't been formally reaped by the OS yet, which is a
    bookkeeping detail the caller shouldn't need to care about.
    """
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    if _is_zombie(pid):
        return False
    return True


def is_run_in_progress(project_name):
    """
    True if the status file says the pipeline is "queued" (launched but
    the background process hasn't yet reached its first status write --
    see the module docstring's note in launch_background_run about why
    "queued" exists) OR "running", AND the PID recorded at launch time
    still appears to be alive -- guards against a stale "running"/
    "queued" status left behind by a process that was killed/crashed
    without getting a chance to write an "error" status (e.g. a server
    restart or OOM-kill), which would otherwise make is_run_in_progress
    wrongly report True forever and block the user from ever
    relaunching.
    """
    status = _load_status(project_name)
    if status["pipeline_status"] not in ("queued", "running"):
        return False
    launch_info = get_launch_info(project_name)
    if not launch_info:
        return False
    return is_process_alive(launch_info.get("pid"))


def launch_background_run(project_name, config, python_executable=None):
    """
    Save `config` to disk and launch this module's CLI entry point
    (see the "if __name__ == '__main__'" block below) as its own,
    fully independent OS process via subprocess.Popen -- NOT a
    Streamlit-managed subprocess, so it keeps running even if the
    Streamlit server restarts or every browser tab is closed.

    start_new_session=True (POSIX) detaches the child into its own
    process group/session, so it is not killed by signals sent to (or
    terminal hangup of) the parent Streamlit process's process group.

    stdout/stderr are redirected to log_path(project_name) (appended,
    not truncated, so re-launches after a resume don't destroy earlier
    log context) rather than inherited, since a detached background
    process has no terminal to write to anyway.

    IMPORTANT -- "queued" status written BEFORE Popen: a freshly
    spawned Python subprocess needs a moment to actually start
    executing (interpreter startup + importing pandas/reference_manager/
    quantification_manager/etc. inside run_pipeline can easily take a
    second or more) before it reaches run_pipeline()'s own first
    _save_status() call. If a caller immediately reruns/re-renders a
    Streamlit page right after calling this function (the normal
    pattern -- see advanced_mode_workspace.py's launch button), a
    naive implementation would read back a stale, pre-existing status
    file (or none at all, i.e. "not_started") for that brief window,
    incorrectly showing the full configuration form again instead of
    the status panel -- confirmed directly via real usability testing,
    where the status screen didn't appear until a manual page refresh.
    Writing an explicit "queued" status HERE, synchronously, before
    Popen returns, closes that window entirely: by the time this
    function returns (and thus by the time any immediately-following
    st.rerun() re-renders), the status file already reflects a
    real, in-progress run.

    Returns the launched process's PID (int).
    """
    save_config(project_name, config)

    status = _load_status(project_name)
    status["pipeline_status"] = "queued"
    status["current_stage"] = None
    status["started_at"] = None
    status["error"] = None
    for stage_key, _ in PIPELINE_STAGES:
        status["stages"][stage_key] = {"status": "pending", "message": None, "started_at": None, "finished_at": None}
    _save_status(project_name, status)

    python_executable = python_executable or sys.executable
    this_script = os.path.abspath(__file__)

    log_file = open(log_path(project_name), "a")
    log_file.write(f"\n\n===== Launching Advanced Mode run at {_now()} =====\n")
    log_file.flush()

    popen_kwargs = {}
    if os.name == "posix":
        popen_kwargs["start_new_session"] = True

    proc = subprocess.Popen(
        [python_executable, this_script, project_name, config_path(project_name)],
        stdout=log_file, stderr=subprocess.STDOUT,
        **popen_kwargs,
    )
    _save_launch_info(project_name, proc.pid)
    return proc.pid


# ---------------------------------------------------------------------------
# Config shape (a plain dict, saved to disk as JSON by whatever UI/CLI
# creates the Advanced Mode project, then loaded here) -- see this
# module's docstring for the intended subprocess.Popen invocation.
#
# {
#   "fastq_source": "directory" | "sra",
#   "fastq_source_dir": str,       # required if fastq_source == "directory"
#   "sra_accessions": [str, ...],  # required if fastq_source == "sra"
#   "sra_max_workers": int,        # optional, default 3 -- how many SRA
#                                  # runs to download in parallel (only
#                                  # used if fastq_source == "sra")
#   "metadata_path": str,          # path to a .csv/.xlsx metadata file
#   "reference": {
#       "is_custom": bool,
#       "species_key": str or None,          # preset catalog key, e.g. "human"
#       "custom_genome_fasta": str or None,  # required for custom + STAR
#       "custom_gtf": str or None,           # required for custom + STAR,
#                                             # and for custom + Salmon
#                                             # unless custom_transcript_fasta given
#       "custom_transcript_fasta": str or None,  # optional pre-extracted
#                                                 # transcript FASTA for custom + Salmon
#       "no_introns": bool,   # custom organism has no introns (identity tx2gene)
#   },
#   "alignment_method": "salmon" | "star",
#   "threads": {  # all optional, sensible defaults applied if omitted
#       "fastqc": 4, "fastp": 3,
#       "salmon_index": 4, "salmon_quant": 4,
#       "star_index": 4, "star_align": 4,
#   },
#   "auto_fix_poly_tails": bool,   # default True -- automatically detect
#                                  # AND re-trim residual poly-G/poly-A
#                                  # tails after the main trimming step
#                                  # (see _run_trimming below)
#   "salmon_count_type": "NumReads" | "TPM",   # default "NumReads"
#   "use_tximport": bool,                       # default True
#   "star_options": {
#       "use_two_pass": bool, "use_encode_options": bool,
#       "add_strand_field": bool,
#   },
# }
# ---------------------------------------------------------------------------

DEFAULT_THREADS = {
    "fastqc": 4, "fastp": 3,
    "salmon_index": 4, "salmon_quant": 4,
    "star_index": 4, "star_align": 4,
}


def _threads(config, key):
    return config.get("threads", {}).get(key, DEFAULT_THREADS[key])


# ---------------------------------------------------------------------------
# Stage 1: FASTQ ingestion + sample/metadata matching
# ---------------------------------------------------------------------------

def _run_ingest(project_name, config):
    fastq_dir = pm.fastq_dir(project_name)
    os.makedirs(fastq_dir, exist_ok=True)

    source = config["fastq_source"]
    if source == "directory":
        n_linked, n_skipped, skipped = ingest.symlink_fastq_files_from_directory(
            config["fastq_source_dir"], fastq_dir
        )
        if n_linked == 0 and n_skipped == 0:
            raise RuntimeError(
                f"No FASTQ files found in source directory '{config['fastq_source_dir']}'."
            )
    elif source == "sra":
        accessions = config["sra_accessions"]
        # Returns a list of (accession, success, message) tuples -- see
        # sra_manager.download_sra_runs_parallel's docstring.
        # sra_max_workers is user-configurable (see the "Parallel
        # downloads" slider in advanced_mode_workspace.py's NCBI/SRA
        # source step) since the right value depends heavily on the
        # user's own network connection and how many runs are being
        # fetched -- defaults to 3 (a moderate, safe starting point,
        # matching the interactive workspace's own default) if omitted.
        results = sra.download_sra_runs_parallel(
            accessions, fastq_dir,
            max_workers=config.get("sra_max_workers", 3),
            threads_per_run=_threads(config, "fastp"),
        )
        failed = [(acc, msg) for acc, success, msg in results if not success]
        if len(failed) == len(accessions):
            raise RuntimeError(f"All {len(accessions)} SRA download(s) failed: {failed}")
    else:
        raise ValueError(f"Unknown fastq_source: {source!r}")

    # Metadata
    raw_meta_df, read_error = ingest.read_metadata_file(config["metadata_path"])
    if read_error:
        raise RuntimeError(f"Failed to read metadata file: {read_error}")
    if "sample" not in raw_meta_df.columns:
        raise RuntimeError(
            "Metadata file must have a 'sample' column matching FASTQ-derived sample names."
        )

    # Match + write samplesheet
    all_fastq_names = ingest.list_existing_fastq(fastq_dir)
    sample_pairs = ingest.validate_sample_pairs(all_fastq_names)
    match_table = ingest.build_match_table(sample_pairs, raw_meta_df)

    n_matched = (match_table["Status"] == "✅ Matched").sum()
    if n_matched == 0:
        raise RuntimeError(
            "No samples matched between FASTQ files and metadata rows. "
            f"Match table:\n{match_table.to_string(index=False)}"
        )

    samplesheet_path = pm.samplesheet_path(project_name)
    samplesheet_df = ingest.write_matched_samplesheet(
        sample_pairs, raw_meta_df, fastq_dir, samplesheet_path
    )

    pm.mark_step_complete(project_name, "samples_matched")
    return f"{n_matched} sample(s) matched and written to samplesheet.csv"


# ---------------------------------------------------------------------------
# Stage 2: Pre-trim QC (FastQC + MultiQC)
# ---------------------------------------------------------------------------

def _run_qc(project_name, config):
    fastqc_ok, multiqc_ok = fastqc.tools_available()
    if not (fastqc_ok and multiqc_ok):
        raise RuntimeError("FastQC and/or MultiQC not found on PATH.")

    fastq_dir = pm.fastq_dir(project_name)
    fastqc_dir = pm.fastqc_dir(project_name)
    multiqc_dir = pm.multiqc_dir(project_name)

    samplesheet_df = pd.read_csv(pm.samplesheet_path(project_name))
    fastq_paths = []
    for _, row in samplesheet_df.iterrows():
        for col in ("fastq_1", "fastq_2"):
            val = row.get(col, "")
            if isinstance(val, str) and val.strip():
                fastq_paths.append(val)

    success, log = fastqc.run_fastqc(fastq_paths, fastqc_dir, threads=_threads(config, "fastqc"))
    if not success:
        raise RuntimeError(f"FastQC failed:\n{log}")

    success, log = fastqc.run_multiqc(fastqc_dir, multiqc_dir)
    if not success:
        raise RuntimeError(f"MultiQC failed:\n{log}")

    pm.mark_step_complete(project_name, "qc_complete")
    return f"FastQC + MultiQC completed for {len(fastq_paths)} file(s)."


# ---------------------------------------------------------------------------
# Stage 3: Trimming (fastp) + post-trim QC
# ---------------------------------------------------------------------------

def _run_trimming(project_name, config):
    fastp_ok, multiqc_ok = fastp.tools_available()
    if not (fastp_ok and multiqc_ok):
        raise RuntimeError("fastp and/or MultiQC not found on PATH.")

    samplesheet_df = pd.read_csv(pm.samplesheet_path(project_name))
    trimmed_dir = pm.trimmed_fastq_dir(project_name)
    reports_dir = pm.fastp_reports_dir(project_name)
    multiqc_dir = pm.posttrim_multiqc_dir(project_name)

    n_trimmed = 0
    n_skipped = 0
    failures = []
    for _, row in samplesheet_df.iterrows():
        if fastp.sample_already_trimmed(row, trimmed_dir):
            n_skipped += 1
            continue
        success, log = fastp.run_fastp_for_sample(
            row, trimmed_dir, reports_dir, threads=_threads(config, "fastp")
        )
        if success:
            n_trimmed += 1
        else:
            failures.append((row["sample"], log))

    if n_trimmed == 0 and n_skipped == 0:
        raise RuntimeError(f"Trimming failed for every sample: {failures}")

    # --- Automatic residual poly-G/poly-A tail detection + re-trim ---
    # Mirrors the interactive Trimming workspace's "Auto Re-trim Flagged
    # Sample(s)" button (trimming_workspace.py) -- but since there's no
    # user present in a background Advanced Mode run to click that
    # button, this does the same detect-then-fix work AUTOMATICALLY,
    # gated by config["auto_fix_poly_tails"] (default True -- see
    # advanced_mode_workspace.py's toggle and its help text for why
    # leaving this on is recommended for virtually all use cases).
    #
    # fastp.scan_all_samples_for_poly_tail_issues() re-analyzes the
    # fastp JSON reports ALREADY written by the trimming loop above
    # (no extra tool/pass needed) for an abnormal G/A concentration in
    # the last several cycles of the TRIMMED reads -- the signature of
    # a poly-G (2-color chemistry artifact) or poly-A (mRNA poly-A
    # tail read-through) tail that fastp's default settings didn't
    # catch (see fastp_manager.py's module-level comment for exactly
    # why neither is trimmed by default). Any sample flagged here is
    # re-trimmed from its ORIGINAL raw FASTQ (never chained on top of
    # the already-trimmed output -- see run_fastp_for_sample's
    # docstring) with ONLY the specific base(s)/flag(s) that sample was
    # actually flagged for, so samples without this issue are never
    # unnecessarily re-processed or over-trimmed.
    n_poly_fixed = 0
    poly_fix_failures = []
    if config.get("auto_fix_poly_tails", True):
        poly_tail_issues = fastp.scan_all_samples_for_poly_tail_issues(reports_dir)
        if poly_tail_issues:
            sample_row_lookup = {row["sample"]: row for _, row in samplesheet_df.iterrows()}
            for sample_name, issues in poly_tail_issues.items():
                sample_row = sample_row_lookup.get(sample_name)
                if sample_row is None:
                    poly_fix_failures.append((sample_name, "Sample no longer found in the matched sample list."))
                    continue
                bases_flagged = {issue["base"] for issue in issues}
                success, log = fastp.run_fastp_for_sample(
                    sample_row, trimmed_dir, reports_dir, threads=_threads(config, "fastp"),
                    force_trim_poly_g=("G" in bases_flagged), force_trim_poly_x=("A" in bases_flagged),
                )
                if success:
                    n_poly_fixed += 1
                else:
                    poly_fix_failures.append((sample_name, log))

    # MultiQC is (re-)generated AFTER any poly-tail re-trims above, so
    # the combined report reflects each sample's FINAL trimmed state --
    # not the pre-fix intermediate result for any sample that got
    # re-trimmed.
    success, log = fastp.run_multiqc(reports_dir, multiqc_dir)
    if not success:
        raise RuntimeError(f"Post-trim MultiQC failed:\n{log}")

    pm.mark_step_complete(project_name, "trimming_complete")
    msg = f"{n_trimmed} sample(s) trimmed, {n_skipped} already done (resumed)."
    if failures:
        msg += f" {len(failures)} sample(s) FAILED: {[f[0] for f in failures]}"
    if config.get("auto_fix_poly_tails", True):
        if n_poly_fixed:
            msg += f" Auto-fixed residual poly-G/poly-A tail(s) in {n_poly_fixed} sample(s)."
        if poly_fix_failures:
            msg += f" Poly-tail re-trim FAILED for {len(poly_fix_failures)} sample(s): {[f[0] for f in poly_fix_failures]}"
    return msg


# ---------------------------------------------------------------------------
# Stage 4: Reference setup (download/build index for preset OR custom)
# ---------------------------------------------------------------------------

def _run_reference(project_name, config):
    ref_cfg = config["reference"]
    method = config["alignment_method"]
    is_custom = ref_cfg.get("is_custom", False)
    species_key = ref_cfg.get("species_key")

    pm.save_alignment_method(project_name, method)
    pm.save_reference_choice(project_name, species_key, is_custom=is_custom)

    if not is_custom and species_key:
        reference_dir = pm.shared_reference_dir(species_key)
    else:
        reference_dir = pm.reference_dir(project_name)
    os.makedirs(reference_dir, exist_ok=True)

    if method == "salmon":
        if not is_custom and species_key:
            cdna_dir = pm.shared_cdna_fasta_dir(species_key)

            def _build_cdna(temp_dir):
                success, fasta_path, message = rm.get_transcriptome_fasta_for_salmon(species_key, temp_dir)
                if success and not rm.validate_fasta_file(fasta_path):
                    return False, "Downloaded file failed FASTA validation."
                return success, message

            success, message, built = rm.ensure_shared_resource(cdna_dir, _build_cdna)
            if not success:
                raise RuntimeError(f"Reference download failed: {message}")
            transcriptome_fasta = os.path.join(cdna_dir, f"{species_key}.cdna.fa")
            index_dir = pm.shared_salmon_index_dir(species_key)
        else:
            transcriptome_fasta = ref_cfg.get("custom_transcript_fasta")
            if not transcriptome_fasta:
                genome_fasta = ref_cfg["custom_genome_fasta"]
                gtf_path = ref_cfg["custom_gtf"]
                transcriptome_fasta = os.path.join(reference_dir, "custom_extracted_transcripts.fa")
                ok, extraction_message = rm.extract_transcripts_with_gffread(genome_fasta, gtf_path, transcriptome_fasta)
                if not ok or not os.path.exists(transcriptome_fasta):
                    raise RuntimeError(f"gffread transcript extraction failed for custom reference: {extraction_message}")
            index_dir = pm.salmon_index_dir(project_name)

        if not qm.salmon_index_exists(index_dir):
            def _build_index(temp_dir_or_final):
                return qm.build_salmon_index(transcriptome_fasta, temp_dir_or_final, threads=_threads(config, "salmon_index"))

            if not is_custom and species_key:
                success, message, built = rm.ensure_shared_resource(index_dir, _build_index)
                if not success:
                    raise RuntimeError(f"Salmon indexing failed: {message}")
            else:
                success, log = qm.build_salmon_index(transcriptome_fasta, index_dir, threads=_threads(config, "salmon_index"))
                if not success:
                    raise RuntimeError(f"Salmon indexing failed:\n{log}")

    else:  # star
        if not is_custom and species_key:
            genome_dir = pm.shared_genome_dir(species_key)

            def _build_genome(temp_dir):
                # download_genome_and_gtf returns a 3-tuple
                # (success, (genome_path, gtf_path)_or_None, message),
                # but ensure_shared_resource's build_fn contract
                # requires exactly (success, message) -- unwrap here.
                success, _paths, message = rm.download_genome_and_gtf(species_key, temp_dir)
                return success, message

            success, message, built = rm.ensure_shared_resource(genome_dir, _build_genome)
            if not success:
                raise RuntimeError(f"Genome download failed: {message}")
            genome_fasta = os.path.join(genome_dir, f"{species_key}.genome.fa")
            gtf_path = os.path.join(genome_dir, f"{species_key}.annotation.gtf")
            index_dir = pm.shared_star_index_dir(species_key)
        else:
            genome_fasta = ref_cfg["custom_genome_fasta"]
            gtf_path = ref_cfg["custom_gtf"]
            index_dir = pm.star_index_dir(project_name)

        if not qm.star_index_exists(index_dir):
            # Read length affects sjdbOverhang; use a representative
            # trimmed FASTQ if trimming has already run, otherwise fall
            # back to STAR/quantification_manager's own default.
            trimmed_dir = pm.trimmed_fastq_dir(project_name)
            sjdb_overhang = 100
            if os.path.isdir(trimmed_dir):
                sample_files = os.listdir(trimmed_dir)
                if sample_files:
                    read_len = qm.detect_fastq_read_length(os.path.join(trimmed_dir, sample_files[0]))
                    if read_len:
                        sjdb_overhang = max(read_len - 1, 1)

            def _build_index(temp_dir_or_final):
                return qm.build_star_index(
                    genome_fasta, gtf_path, temp_dir_or_final,
                    threads=_threads(config, "star_index"), sjdb_overhang=sjdb_overhang,
                )

            if not is_custom and species_key:
                success, message, built = rm.ensure_shared_resource(index_dir, _build_index)
                if not success:
                    raise RuntimeError(f"STAR indexing failed: {message}")
            else:
                success, log = qm.build_star_index(
                    genome_fasta, gtf_path, index_dir,
                    threads=_threads(config, "star_index"), sjdb_overhang=sjdb_overhang,
                )
                if not success:
                    raise RuntimeError(f"STAR indexing failed:\n{log}")

    pm.mark_step_complete(project_name, "reference_ready")
    return f"Reference + {method} index ready at {reference_dir}"


# ---------------------------------------------------------------------------
# Stage 5: Quantification (Salmon or STAR, per sample)
# ---------------------------------------------------------------------------

def _run_quantification(project_name, config):
    method = config["alignment_method"]
    species_key, is_custom = pm.get_reference_choice(project_name)
    samplesheet_df = pd.read_csv(pm.samplesheet_path(project_name))
    trimmed_dir = pm.trimmed_fastq_dir(project_name)
    manifest = qm.build_sample_manifest(samplesheet_df, trimmed_dir)

    n_success = 0
    failures = []

    if method == "salmon":
        index_dir = pm.shared_salmon_index_dir(species_key) if (not is_custom and species_key) else pm.salmon_index_dir(project_name)
        quant_dir = pm.salmon_quant_dir(project_name)
        for entry in manifest:
            sample_out = os.path.join(quant_dir, entry["sample"], "quant.sf")
            if os.path.exists(sample_out):
                n_success += 1
                continue
            success, log, _ = qm.run_salmon_quant(entry, index_dir, quant_dir, threads=_threads(config, "salmon_quant"))
            if success:
                n_success += 1
            else:
                failures.append((entry["sample"], log))
    else:
        index_dir = pm.shared_star_index_dir(species_key) if (not is_custom and species_key) else pm.star_index_dir(project_name)
        align_dir = pm.star_align_dir(project_name)
        star_opts = config.get("star_options", {})
        for entry in manifest:
            sample_out = os.path.join(align_dir, entry["sample"], f"{entry['sample']}_ReadsPerGene.out.tab")
            if os.path.exists(sample_out):
                n_success += 1
                continue
            success, log, _ = qm.run_star_align(
                entry, index_dir, align_dir, threads=_threads(config, "star_align"),
                use_two_pass=star_opts.get("use_two_pass", False),
                use_encode_options=star_opts.get("use_encode_options", False),
                add_strand_field=star_opts.get("add_strand_field", False),
            )
            if success:
                n_success += 1
            else:
                failures.append((entry["sample"], log))

    if n_success == 0:
        raise RuntimeError(f"Quantification failed for every sample: {failures}")

    pm.mark_step_complete(project_name, "quantification_complete")
    msg = f"{n_success}/{len(manifest)} sample(s) quantified successfully."
    if failures:
        msg += f" FAILED: {[f[0] for f in failures]}"
    return msg


# ---------------------------------------------------------------------------
# Stage 6: Combine into a gene counts matrix (final Advanced Mode output)
# ---------------------------------------------------------------------------

def _run_counts_matrix(project_name, config):
    method = config["alignment_method"]
    species_key, is_custom = pm.get_reference_choice(project_name)
    samplesheet_df = pd.read_csv(pm.samplesheet_path(project_name))
    sample_names = sorted(samplesheet_df["sample"].astype(str).tolist())
    reference_dir = pm.shared_reference_dir(species_key) if (not is_custom and species_key) else pm.reference_dir(project_name)
    counts_path = pm.counts_matrix_path(project_name)

    if method == "salmon":
        count_column = "NumReads" if config.get("salmon_count_type", "NumReads") == "NumReads" else "TPM"
        use_tximport = config.get("use_tximport", True) and count_column == "NumReads"
        quant_dir = pm.salmon_quant_dir(project_name)

        matrix_df = None
        missing = []
        if use_tximport:
            tx2gene, tx2gene_source = cmm.get_tx2gene_mapping(project_name, reference_dir)
            if tx2gene:
                sample_quant_paths = {}
                for s in sample_names:
                    p = os.path.join(quant_dir, s, "quant.sf")
                    if os.path.exists(p):
                        sample_quant_paths[s] = p
                    else:
                        missing.append(s)
                if sample_quant_paths:
                    tx2gene_path = rm.save_tx2gene_csv(tx2gene, os.path.join(reference_dir, "tx2gene", "tx2gene.csv"))
                    work_dir = os.path.join(reference_dir, "tximport_work")
                    success, log = qm.run_tximport_gene_collapse(sample_quant_paths, tx2gene_path, counts_path, work_dir)
                    if success:
                        matrix_df = pd.read_csv(counts_path)
                    # else: fall through to direct merge below

        if matrix_df is None:
            matrix_df, extra_missing = cmm.merge_salmon_counts(sample_names, quant_dir, count_column=count_column)
            missing.extend(extra_missing)
    else:
        align_dir = pm.star_align_dir(project_name)
        matrix_df, missing = cmm.merge_star_counts(sample_names, align_dir)

    if matrix_df is None or matrix_df.empty:
        raise RuntimeError(f"No quantification output found for any sample. Missing: {missing}")

    os.makedirs(os.path.dirname(counts_path), exist_ok=True)
    matrix_df.to_csv(counts_path, index=False)
    pm.mark_step_complete(project_name, "counts_matrix_complete")

    # Best-effort gene symbol mapping, same as the interactive workspace.
    gene_symbol_map, gene_symbol_source = cmm.get_gene_symbol_mapping(project_name, reference_dir, method)
    if gene_symbol_map:
        rm.save_gene_symbol_map_csv(gene_symbol_map, pm.gene_symbol_map_path(project_name))
        pm.save_gene_id_mapping_meta(project_name, {"source": "auto_parse", "detail": gene_symbol_source})

    n_genes = len(matrix_df)
    n_samples = len(matrix_df.columns) - 1
    msg = f"Counts matrix built: {n_genes:,} gene(s) × {n_samples} sample(s)."
    if missing:
        msg += f" Skipped (no quant output): {missing}"
    return msg


STAGE_FUNCS = {
    "ingest": _run_ingest,
    "qc": _run_qc,
    "trimming": _run_trimming,
    "reference": _run_reference,
    "quantification": _run_quantification,
    "counts_matrix": _run_counts_matrix,
}


# ---------------------------------------------------------------------------
# Main pipeline runner
# ---------------------------------------------------------------------------

def run_pipeline(project_name, config):
    """
    Run every pipeline stage in order for project_name, resuming from
    wherever project_manager's step-tracking says this project last
    left off, and STOPPING after counts_matrix_complete (never runs
    DESeq2 or Ontology Analysis -- see module docstring).

    config: a dict matching the shape documented above this module's
        "Config shape" comment block.

    Returns the final status dict (same shape written to
    status_path(project_name)). Raises nothing -- all stage errors are
    caught, recorded in the status file's "error" field, and cause this
    function to return with pipeline_status == "error" rather than
    propagating, so a caller running this in a background process can
    simply check the returned/persisted status rather than needing a
    try/except around this call.
    """
    if not pm.list_projects() or project_name not in pm.list_projects():
        pm.create_project(project_name)

    status = _load_status(project_name)
    status["pipeline_status"] = "running"
    if status["started_at"] is None:
        status["started_at"] = _now()
    status["error"] = None
    _save_status(project_name, status)

    for stage_key, step_name in PIPELINE_STAGES:
        if pm.has_completed_step(project_name, step_name):
            _set_stage_status(project_name, status, stage_key, "skipped", "Already completed previously.")
            continue

        _set_stage_status(project_name, status, stage_key, "running", "In progress...")
        try:
            message = STAGE_FUNCS[stage_key](project_name, config)
            _set_stage_status(project_name, status, stage_key, "complete", message)
        except Exception as e:
            error_detail = {"stage": stage_key, "message": str(e), "traceback": traceback.format_exc()}
            status["error"] = error_detail
            status["pipeline_status"] = "error"
            _set_stage_status(project_name, status, stage_key, "error", str(e))
            return status

    status["pipeline_status"] = "complete"
    status["current_stage"] = None
    _save_status(project_name, status)
    return status


# ---------------------------------------------------------------------------
# CLI entry point -- for launching via subprocess.Popen as its own
# independent OS process (see module docstring).
#
#   python advanced_mode_orchestrator.py <project_name> <config_json_path>
#
# Exit code 0 on pipeline_status == "complete", 1 otherwise (error, or
# unexpected top-level exception) -- lets a launcher/wrapper script
# check success via the process return code in addition to polling the
# status file.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python advanced_mode_orchestrator.py <project_name> <config_json_path>")
        sys.exit(2)

    _project_name, _config_path = sys.argv[1], sys.argv[2]
    with open(_config_path) as _f:
        _config = json.load(_f)

    _final_status = run_pipeline(_project_name, _config)
    sys.exit(0 if _final_status["pipeline_status"] == "complete" else 1)
