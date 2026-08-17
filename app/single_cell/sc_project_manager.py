"""
single_cell/sc_project_manager.py

Lightweight "project" system for the Single-cell RNA-Seq pipeline --
deliberately mirrors project_manager.py's exact interface (list_projects,
project_dir, create_project, load_info/save_info, mark_step_complete,
has_completed_step) so anyone already familiar with the Bulk RNA-Seq
project system recognizes this immediately, but points at its OWN root
(SC_PROJECTS_ROOT) so single-cell and bulk projects can never collide by
name or get mixed up in the same listing.

get_recommended_thread_count() is intentionally NOT duplicated here --
single_cell_workspace.py imports and reuses project_manager's own copy
directly, since detecting available CPU cores has nothing pipeline-
specific about it.

--- Path layout for a single-cell project ---
data/singlecell_projects/<project>/
    fastq/                      raw FASTQ files (R1 = barcode+UMI, R2 = cDNA)
    metadata.csv                sample-level metadata (condition/treatment/donor)
    qc/fastqc/                  pre-trim FastQC (R2-focused)
    qc/multiqc/                 pre-trim MultiQC
    trimmed/                    post-trim FASTQ (R1 passed through untouched,
                                 R2 actually trimmed -- see singlecell_trim_manager.py)
    qc/fastp/                   fastp's own per-sample JSON/HTML reports
    qc/multiqc_posttrim/        post-trim MultiQC
    custom_reference/           uploaded custom genome FASTA + annotation
    star_index/                 per-project STAR/STARsolo index -- ONLY used
                                 for a CUSTOM (non-preset) reference; preset
                                 species reuse the SAME shared index the bulk
                                 pipeline builds (pm.shared_star_index_dir())
    starsolo/                   STARsolo output (per-sample Solo.out/ dirs)
    alevin_fry/                 alevin-fry output, when that alternate
                                 aligner is selected instead of STARsolo
    project_info.json           steps_completed, chemistry choice, aligner
                                 choice, reference choice, etc.
"""
import json
import os
from datetime import datetime

SC_PROJECTS_ROOT = "data/singlecell_projects"


def list_projects():
    "Return a sorted list of existing single-cell project names."
    if not os.path.isdir(SC_PROJECTS_ROOT):
        return []
    return sorted([
        name for name in os.listdir(SC_PROJECTS_ROOT)
        if os.path.isdir(os.path.join(SC_PROJECTS_ROOT, name))
    ])


def project_dir(project_name):
    return os.path.join(SC_PROJECTS_ROOT, project_name)


def fastq_dir(project_name):
    return os.path.join(project_dir(project_name), "fastq")


def metadata_path(project_name):
    "Sample-level metadata (condition/treatment/donor) -- NOT per-cell metadata, which is generated downstream by the pipeline itself."
    return os.path.join(project_dir(project_name), "metadata.csv")


def fastqc_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "fastqc")


def multiqc_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "multiqc")


def trimmed_fastq_dir(project_name):
    "Post-trim FASTQ output -- R1 copied through untouched, R2 actually trimmed. See singlecell_trim_manager.py's module docstring for why R1 must never be modified."
    return os.path.join(project_dir(project_name), "trimmed")


def fastp_reports_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "fastp")


def posttrim_multiqc_dir(project_name):
    return os.path.join(project_dir(project_name), "qc", "multiqc_posttrim")


def posttrim_fastqc_dir(project_name):
    "Post-trim FastQC output (separate from the pre-trim fastqc_dir() above), used for the Step 4 post-trim QC re-check."
    return os.path.join(project_dir(project_name), "qc", "fastqc_posttrim")


def custom_reference_dir(project_name):
    """
    Where a CUSTOM (non-preset, user-uploaded) reference's genome FASTA
    + gene annotation (GTF/GFF/GFF3) are saved when provided via Step
    5's "Upload from my computer" source option -- project-scoped, since
    two different single-cell projects' custom uploads have no
    guarantee of being the same organism/assembly (mirrors the Bulk
    RNA-Seq pipeline's own reference_dir() -- deliberately NOT shared
    across projects the way preset/catalog species are).
    """
    return os.path.join(project_dir(project_name), "custom_reference")


def star_index_dir(project_name):
    """
    Per-project STAR/STARsolo genome index location -- used ONLY for a
    CUSTOM (non-preset) reference. Preset species reuse the bulk
    pipeline's OWN shared, project-independent index location instead
    (project_manager.shared_star_index_dir(species_key)) -- see
    singlecell_workspace.py's Step 6, which mirrors
    alignment_workspace.py's exact same shared-vs-per-project index
    selection logic (confirmed against that real module's source).
    Mirrors project_manager.py's own star_index_dir(project) naming
    convention for a custom-reference bulk project.
    """
    return os.path.join(project_dir(project_name), "star_index")


def starsolo_output_dir(project_name):
    "Root directory for STARsolo's per-sample Solo.out/ output (cell x gene matrices, cell-calling summary, etc.)."
    return os.path.join(project_dir(project_name), "starsolo")


def alevin_fry_output_dir(project_name):
    "Root directory for alevin-fry's per-sample output, when that alternate aligner is selected instead of STARsolo."
    return os.path.join(project_dir(project_name), "alevin_fry")


def info_path(project_name):
    return os.path.join(project_dir(project_name), "project_info.json")


def create_project(project_name):
    "Create a new single-cell project's folder structure. Returns True if created, False if a project with that name already exists."
    d = project_dir(project_name)
    if os.path.exists(d):
        return False
    os.makedirs(fastq_dir(project_name), exist_ok=True)
    os.makedirs(fastqc_dir(project_name), exist_ok=True)
    os.makedirs(multiqc_dir(project_name), exist_ok=True)
    save_info(project_name, {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "steps_completed": [],
        "pipeline_type": "single_cell",
    })
    return True


def load_info(project_name):
    path = info_path(project_name)
    if not os.path.exists(path):
        return {"created_at": None, "steps_completed": [], "pipeline_type": "single_cell"}
    with open(path) as f:
        return json.load(f)


def save_info(project_name, info):
    os.makedirs(project_dir(project_name), exist_ok=True)
    with open(info_path(project_name), "w") as f:
        json.dump(info, f, indent=2)


def mark_step_complete(project_name, step_name):
    info = load_info(project_name)
    if step_name not in info.get("steps_completed", []):
        info.setdefault("steps_completed", []).append(step_name)
        info["last_updated"] = datetime.now().isoformat(timespec="seconds")
    save_info(project_name, info)


def has_completed_step(project_name, step_name):
    info = load_info(project_name)
    return step_name in info.get("steps_completed", [])


def save_chemistry_choice(project_name, chemistry_key, was_auto_detected, user_confirmed):
    info = load_info(project_name)
    info["chemistry"] = {
        "chemistry_key": chemistry_key,
        "was_auto_detected": was_auto_detected,
        "user_confirmed": user_confirmed,
    }
    save_info(project_name, info)


def get_chemistry_choice(project_name):
    info = load_info(project_name)
    return info.get("chemistry")


def save_sample_column(project_name, column_name):
    info = load_info(project_name)
    info["metadata_sample_column"] = column_name
    save_info(project_name, info)


def get_sample_column(project_name):
    info = load_info(project_name)
    return info.get("metadata_sample_column")


def save_fastq_source_dir(project_name, directory):
    """
    Remember which directory this project's raw FASTQ files actually live
    in -- needed because Step 1 (ingestion) allows EITHER uploading files
    directly into this project's own fastq_dir(), OR browsing an
    arbitrary external directory elsewhere on the server. Later steps
    (trimming, alignment) are separate sidebar pages/reruns and need to
    re-discover the same R1/R2 pairs Step 1 found, without re-running
    Step 1's own upload/browse UI -- persisting the resolved directory
    here is what makes that possible.
    """
    info = load_info(project_name)
    info["fastq_source_dir"] = directory
    save_info(project_name, info)


def get_fastq_source_dir(project_name):
    info = load_info(project_name)
    return info.get("fastq_source_dir")


def save_aligner_choice(project_name, aligner):
    """
    Persist Step 5's confirmed alignment/quantification method choice
    (ALIGNER_STARSOLO / ALIGNER_ALEVIN_FRY, from singlecell_workspace.py)
    -- added 2026-08-17 alongside the fix that moved this choice to the
    top of Step 5 (see singlecell_workspace.py's module docstring,
    "Aligner-choice / genome-index ordering fix").

    Without this, reopening a project after a restart would always fall
    back to whatever the radio's hardcoded default is (STARsolo), silently
    discarding a previously-confirmed alevin-fry selection and forcing the
    user to re-pick it every time even though nothing else about the
    project changed.

    Saved as its own top-level "aligner" key (NOT nested under
    "reference", even though the choice affects reference/index setup
    downstream) since it's conceptually a distinct pipeline setting --
    mirrors save_chemistry_choice's own separate top-level "chemistry"
    key above. Also saved unconditionally on every rerun the moment the
    radio renders (see singlecell_workspace.py's _render_aligner_choice),
    the same way save_chemistry_choice is called every rerun in Step 1
    rather than gated behind a separate "confirm" button -- so the very
    latest selection is always what's persisted, with no extra click
    required.
    """
    info = load_info(project_name)
    info["aligner"] = aligner
    save_info(project_name, info)


def get_aligner_choice(project_name):
    "Return the previously confirmed aligner choice ('starsolo'/'alevin_fry'), or None if never set (e.g. a project created before this setting existed, or a brand-new project)."
    info = load_info(project_name)
    return info.get("aligner")


def alignment_results_path(project_name):
    "Path to the persisted per-sample STARsolo alignment/cell-calling results (Cells Detected, Uniquely Mapped %, Quality, etc.)."
    return os.path.join(project_dir(project_name), "alignment_results.json")


def save_alignment_results(project_name, results):
    """
    Persist Step 6's per-sample STARsolo results table (list of dicts --
    Sample/Status/Cells Detected/Uniquely Mapped %/Uniquely Mapped (%)/
    Quality) to disk, separate from the boolean "alignment" entry in
    steps_completed.

    Added 2026-08-17 to fix a real reported bug: previously, ONLY the
    boolean completion flag was persisted -- the actual metrics table
    shown right after a run existed solely in that run's local variable
    and vanished on any later rerun (e.g. reopening the project, or an
    environment/session restart). Combined with singlecell_workspace.py's
    _render_step6 now checking `has_completed_step(project, "alignment")`
    BEFORE (rather than after) any live STAR/reference/index availability
    checks, this makes a previously-completed run's results visible again
    immediately on reopening a project, without needing STAR, the genome
    index, or reference files to still be reachable just to VIEW them
    (those are now only required to actually re-run alignment).
    """
    info = load_info(project_name)
    with open(alignment_results_path(project_name), "w") as f:
        json.dump(results, f, indent=2)
    info["alignment_results_saved"] = True
    save_info(project_name, info)


def get_alignment_results(project_name):
    "Return the previously persisted per-sample alignment results (list of dicts), or None if never saved (e.g. a project whose alignment completed before this feature existed)."
    path = alignment_results_path(project_name)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def save_reference_choice(project_name, reference_cfg):
    """
    Persist Step 5's confirmed reference configuration -- is_custom,
    species_key (preset case) or custom_genome_fasta/custom_gtf paths +
    detected annotation_format ("gtf"/"gff3", custom cases). Without
    this, re-opening a project or a later Step 6 rerun would have no
    way to know which reference was actually confirmed.
    """
    info = load_info(project_name)
    info["reference"] = reference_cfg
    save_info(project_name, info)


def get_reference_choice(project_name):
    info = load_info(project_name)
    return info.get("reference")
