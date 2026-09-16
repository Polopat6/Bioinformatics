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
    downstream/                 Phase 3 downstream-analysis state (see its
                                 own path helpers below):
        current_state.h5ad          single overwritten "current state" cache
        downstream_recipe.json      per-step parameter fingerprints
        checkpoints/                 user-named, independently-persisted .h5ad snapshots
        pseudobulk/<export_name>/    Step 9 pseudobulk counts+metadata CSVs, one
                                     subdirectory per user-named export (since a
                                     user may want several groupby_columns choices
                                     side by side, e.g. "by_sample" vs.
                                     "by_sample_and_celltype", without one
                                     silently overwriting the other)
        compositional/output/       Step 10 propeller/simple-test results CSVs
        compositional/work/         Step 10 R job-spec/script staging files
    project_info.json           steps_completed, chemistry choice, aligner
                                 choice, reference choice, etc.
"""
import json
import os
from datetime import datetime
import app_paths
import atomic_io

SC_PROJECTS_ROOT = app_paths.data_path("singlecell_projects")


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

def sc_reference_gene_symbol_map_path(project_name):
    """
    Where this project's GTF-auto-derived gene_id -> gene_name map lives
    (layer 1 of the SC gene-ID-mapping panel's 3-layer lookup, mirroring
    the bulk pipeline's own layer 1). Project-scoped (not export-scoped):
    this depends only on the project's confirmed reference GTF, which is
    shared across every pseudobulk export in this project.
    """
    return os.path.join(custom_reference_dir(project_name), "gene_symbol_map.csv")

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
    return atomic_io.read_json(
        info_path(project_name),
        default={"created_at": None, "steps_completed": [],
                 "pipeline_type": "single_cell"},
        on_corrupt="raise",
    )


def save_info(project_name, info):
    atomic_io.atomic_write_json(info_path(project_name), info)


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
    atomic_io.atomic_write_json(alignment_results_path(project_name), results)
    info["alignment_results_saved"] = True
    save_info(project_name, info)


def get_alignment_results(project_name):
    "Return the previously persisted per-sample alignment results (list of dicts), or None if never saved (e.g. a project whose alignment completed before this feature existed)."
    return atomic_io.read_json(
        alignment_results_path(project_name), default=None, on_corrupt="default",
    )


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

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def downstream_dir(project_name):
    "Root directory for all Phase 3 (single-cell downstream analysis) state for this project."
    return os.path.join(project_dir(project_name), "downstream")


def downstream_adata_path(project_name):
    """
    The single, default "current state" AnnData cache file -- OVERWRITTEN
    by each completed step (see this module addition's own top-of-file
    docstring, "Caching model", for the full rationale on why this is
    one file, not one-per-step).
    """
    return os.path.join(downstream_dir(project_name), "current_state.h5ad")


def downstream_checkpoint_dir(project_name):
    "Directory holding user-named, explicitly-saved checkpoint .h5ad files (see save_downstream_checkpoint_path below) -- separate from the single default cache file above."
    return os.path.join(downstream_dir(project_name), "checkpoints")


def _sanitize_name(name, what="name"):
    """
    Shared sanitization helper -- letters, numbers, dashes, and
    underscores only (matches this project's existing convention for
    project/sample/checkpoint names). Raises ValueError if nothing
    valid remains after sanitizing, so an empty/whitespace-only input
    never silently collapses into a nonsensical bare filename.
    """
    safe_name = "".join(c for c in (name or "").strip() if c.isalnum() or c in ("-", "_"))
    if not safe_name:
        raise ValueError(f"{what} must contain at least one letter, number, dash, or underscore.")
    return safe_name


def downstream_checkpoint_path(project_name, checkpoint_name):
    """
    Build the file path for a user-named checkpoint. checkpoint_name is
    sanitized the same way project/sample names already are elsewhere
    in this app (letters, numbers, dashes, underscores only).

    Raises ValueError if checkpoint_name has no valid characters at all
    once sanitized -- an empty/whitespace-only name would otherwise
    silently collapse to a nonsensical bare ".h5ad" filename.
    """
    safe_name = _sanitize_name(checkpoint_name, what="Checkpoint name")
    return os.path.join(downstream_checkpoint_dir(project_name), f"{safe_name}.h5ad")


def list_downstream_checkpoints(project_name):
    "Return a sorted list of existing checkpoint names (without the .h5ad extension) for this project."
    d = downstream_checkpoint_dir(project_name)
    if not os.path.isdir(d):
        return []
    return sorted(f[:-5] for f in os.listdir(d) if f.endswith(".h5ad"))


def delete_downstream_checkpoint(project_name, checkpoint_name):
    "Delete a named checkpoint. Returns True if it existed and was deleted, False if there was nothing to delete."
    path = downstream_checkpoint_path(project_name, checkpoint_name)
    if not os.path.isfile(path):
        return False
    os.remove(path)
    return True


def downstream_adata_cache_exists(project_name):
    "Quick check: does this project have ANY Phase 3 progress saved at all (a current-state cache file on disk)? Used by the workspace to decide whether to show a 'resume where you left off' message."
    return os.path.isfile(downstream_adata_path(project_name))


# --- Step 9: Pseudobulk export directories ---
#
# A user may reasonably want SEVERAL pseudobulk exports side by side for
# the same project (e.g. one grouped by ["sample"] alone, another by
# ["sample", "cell_type"] for a cell-type-specific comparison) --
# mirrors downstream_checkpoint_path()'s own "user-named, independently
# persisted" pattern above, rather than a single overwritten location
# the way the main AnnData cache works (that single-overwrite design is
# only correct for the growing AnnData object itself, where every prior
# step's state is already redundantly contained in the latest state --
# it does NOT apply here, since two different pseudobulk exports are
# genuinely independent artifacts a user may want to compare later, not
# a strict superset of one another).

def downstream_pseudobulk_root_dir(project_name):
    "Root directory holding every named pseudobulk export for this project."
    return os.path.join(downstream_dir(project_name), "pseudobulk")


def downstream_pseudobulk_export_dir(project_name, export_name):
    """
    Build the directory path for one named pseudobulk export -- passed
    directly as dsm.save_pseudobulk_for_deseq2()'s own dest_dir argument
    (which itself writes "pseudobulk_counts.csv" and
    "pseudobulk_metadata.csv" inside it).

    export_name is sanitized the same way checkpoint names are (letters,
    numbers, dashes, underscores only) -- raises ValueError if nothing
    valid remains.
    """
    safe_name = _sanitize_name(export_name, what="Pseudobulk export name")
    return os.path.join(downstream_pseudobulk_root_dir(project_name), safe_name)


def list_downstream_pseudobulk_exports(project_name):
    "Return a sorted list of existing pseudobulk export names for this project."
    d = downstream_pseudobulk_root_dir(project_name)
    if not os.path.isdir(d):
        return []
    return sorted([
        name for name in os.listdir(d)
        if os.path.isdir(os.path.join(d, name))
    ])


def delete_downstream_pseudobulk_export(project_name, export_name):
    "Delete a named pseudobulk export directory (counts + metadata CSVs). Returns True if it existed and was deleted, False otherwise."
    import shutil
    path = downstream_pseudobulk_export_dir(project_name, export_name)
    if not os.path.isdir(path):
        return False
    shutil.rmtree(path)
    return True

# ---------------------------------------------------------------------------
# Step 9 -> DE bridge: DESeq2 + Ontology Analysis on a saved pseudobulk export
# ---------------------------------------------------------------------------
#
# Mirrors project_manager.py's own deseq2_dir()/deseq2_output_dir()/
# deseq2_work_dir()/ontology_dir()/ontology_output_dir()/ontology_work_dir()
# naming convention exactly, but scoped one level deeper -- under a
# SPECIFIC pseudobulk export (downstream_pseudobulk_export_dir(project,
# export_name)) rather than directly under the project root. This is
# necessary, not just cosmetic: a single-cell project can have MULTIPLE
# independent pseudobulk exports (e.g. grouped by "sample" vs. by
# "sample"+"cell_type"), and each is really its own independent
# "bulk-like" counts matrix that needs its own independent DESeq2/
# Ontology results -- exactly the same way two different bulk projects
# would never share one set of DESeq2 results.


def sc_reference_gene_symbol_map_path(project_name):
    """
    Where this project's GTF-auto-derived gene_id -> gene_name map lives
    (layer 1 of the SC gene-ID-mapping panel's 3-layer lookup, mirroring
    the bulk pipeline's own layer 1 -- see sc_deseq2_workspace.py's
    module docstring). Project-scoped (not export-scoped): this depends
    only on the project's confirmed reference GTF, which is shared
    across every pseudobulk export in this project.
    """
    return os.path.join(custom_reference_dir(project_name), "gene_symbol_map.csv")

def sc_deseq2_dir(project_name, export_name):
    "Where DESeq2 analysis inputs/outputs live for ONE pseudobulk export of this single-cell project."
    return os.path.join(downstream_pseudobulk_export_dir(project_name, export_name), "deseq2")

def sc_deseq2_output_dir(project_name, export_name):
    "DESeq2's actual result CSVs for this pseudobulk export -- same file shapes deseq2_manager.py writes for a bulk project."
    return os.path.join(sc_deseq2_dir(project_name, export_name), "output")

def sc_deseq2_work_dir(project_name, export_name):
    "Scratch directory for this pseudobulk export's temporary DESeq2 R script + job spec JSON."
    return os.path.join(sc_deseq2_dir(project_name, export_name), "work")

def sc_deseq2_config_path(project_name, export_name):
    """
    Where this pseudobulk export's confirmed DESeq2 configuration is
    saved -- mirrors project_manager.save_deseq2_config()'s own info,
    but stored as its own small JSON file INSIDE this export's deseq2
    directory rather than nested inside project_info.json, since a
    single-cell project's DIFFERENT pseudobulk exports each need their
    own completely independent DESeq2 configuration.
    """
    return os.path.join(sc_deseq2_dir(project_name, export_name), "deseq2_config.json")

def save_sc_deseq2_config(project_name, export_name, config):
    "Persist this pseudobulk export's confirmed DESeq2 configuration -- see sc_deseq2_config_path()'s own docstring."
    atomic_io.atomic_write_json(
        sc_deseq2_config_path(project_name, export_name), config,
    )

def get_sc_deseq2_config(project_name, export_name):
    "Return this pseudobulk export's previously saved DESeq2 configuration, or None if DESeq2 has never been run for it."
    return atomic_io.read_json(
        sc_deseq2_config_path(project_name, export_name),
        default=None, on_corrupt="default",
    )

def sc_ontology_dir(project_name, export_name):
    "Where Ontology Analysis (GO/KEGG/Reactome enrichment) inputs/outputs live for this pseudobulk export."
    return os.path.join(downstream_pseudobulk_export_dir(project_name, export_name), "ontology")

def sc_ontology_output_dir(project_name, export_name):
    "Root directory for this pseudobulk export's Ontology Analysis result files."
    return os.path.join(sc_ontology_dir(project_name, export_name), "output")

def sc_ontology_work_dir(project_name, export_name):
    "Scratch directory for this pseudobulk export's temporary Ontology Analysis R scripts + input gene list CSVs."
    return os.path.join(sc_ontology_dir(project_name, export_name), "work")

def save_sc_ontology_species(project_name, export_name, species_key):
    """
    Persist a manually-confirmed organism choice for this pseudobulk
    export's Ontology Analysis. Deliberately NOT auto-derived from this
    project's alignment reference_choice -- that dict's exact shape
    differs between preset and custom references, and guessing at its
    keys here risks silently resolving the wrong species. Asking the
    user to confirm once (mirroring ontology_workspace.py's own
    _render_species_override_picker pattern for an unrecognized bulk
    reference) is the safer, always-correct approach.
    """
    path = os.path.join(sc_ontology_dir(project_name, export_name), "species_choice.json")
    atomic_io.atomic_write_json(path, {"species_key": species_key})

def get_sc_ontology_species(project_name, export_name):
    path = os.path.join(sc_ontology_dir(project_name, export_name), "species_choice.json")
    saved = atomic_io.read_json(path, default=None, on_corrupt="default")
    return saved.get("species_key") if saved else None
# --- Step 10: Compositional analysis output/work directories ---
#
# Unlike pseudobulk exports above, compositional analysis results are
# NOT modeled as multiple independently-named exports -- a single
# "current result" location is sufficient (mirrors the main AnnData
# cache's own single-overwrite pattern), since re-running with different
# parameters (transform, cluster_key, group_column) is expected to
# simply replace the prior comparison rather than be kept alongside it.
# If side-by-side comparison of multiple compositional runs becomes a
# real need later, this can be extended to a named-export pattern
# identical to the pseudobulk helpers above without disrupting existing
# projects (a fresh sub-directory naming scheme, not a breaking change
# to this single default location).

def downstream_compositional_dir(project_name):
    "Root directory for Step 10 (Compositional Analysis) state for this project."
    return os.path.join(downstream_dir(project_name), "compositional")


def downstream_compositional_output_dir(project_name):
    "Where propeller/simple-test result CSVs are written -- passed directly as dsm.run_propeller_analysis()'s own output_dir argument."
    return os.path.join(downstream_compositional_dir(project_name), "output")


def downstream_compositional_work_dir(project_name):
    "Where the propeller R job-spec JSON and staged R script are written -- passed directly as dsm.run_propeller_analysis()'s own work_dir argument. Kept separate from output_dir so a user Browse of results isn't cluttered with intermediate job-staging files."
    return os.path.join(downstream_compositional_dir(project_name), "work")
# --- bitr() gene-name-conversion work directory ---
#
# Scratch location for the bitr()-based gene-name fallback's own R
# job-spec JSON + staged R script + output CSV (see
# gene_id_mapper.run_bitr_conversion()'s own dest_dir/work_dir
# parameter) -- mirrors downstream_compositional_work_dir()'s own
# "kept separate from any results directory, since these are pure
# intermediate staging files" rationale. No corresponding "output_dir"
# is needed here (unlike Step 10's own output/work split) since a
# bitr() conversion's real, meaningful result is applied directly onto
# adata.var["gene_symbol"] in memory (then persisted via the project's
# own single current_state.h5ad cache) rather than written out as a
# separate standalone results file a user would browse to later.
def downstream_bitr_work_dir(project_name):
    "Where the bitr() gene-name-conversion R job-spec JSON and staged R script are written -- passed directly as gene_id_mapper.run_bitr_conversion()'s own work_dir argument."
    return os.path.join(downstream_dir(project_name), "bitr_work")

# ---------------------------------------------------------------------------
# Per-step parameter "recipe" (fingerprint) system
# ---------------------------------------------------------------------------

# Defines both the SET of tracked Phase 3 steps and, critically, their
# PIPELINE ORDER -- this order is what drives cascading invalidation in
# save_downstream_step_recipe() below. "pseudobulk" and "compositional"
# are placed last since they are terminal, side-analysis consumers of
# whatever clustering/annotation state already exists -- nothing in
# this pipeline is computed ON TOP OF their own output, so they don't
# need to trigger any further cascading invalidation themselves (though
# they still get their own recipe-tracking like every other step, so a
# user changing THEIR OWN parameters -- e.g. a different groupby column
# for pseudobulk -- is still correctly detected).
DOWNSTREAM_STEP_ORDER = [
    "combine",
    "normalize",
    "hvg",
    "pca",
    "batch_correction",
    "neighbors",
    "clustering",
    "embedding",
    "annotation",
    "pseudobulk",
    "compositional",
]


def downstream_recipe_path(project_name):
    return os.path.join(downstream_dir(project_name), "downstream_recipe.json")


def load_downstream_recipe(project_name):
    """
    Load this project's full Phase 3 recipe dict -- {step_name: {"params":
    {...}, "computed_at": "..."}}. Returns an empty dict if no recipe has
    been saved yet at all (e.g. a brand-new project, or one that hasn't
    started Phase 3).
    """
    return atomic_io.read_json(
        downstream_recipe_path(project_name), default={}, on_corrupt="default",
    )


def _save_downstream_recipe_file(project_name, recipe):
    atomic_io.atomic_write_json(downstream_recipe_path(project_name), recipe)


def save_downstream_step_recipe(project_name, step_name, params):
    """
    Save step_name's actual parameter dict (params), stamped with the
    current time, AND cascade-invalidate (delete the saved recipe for)
    every step that comes AFTER step_name in DOWNSTREAM_STEP_ORDER --
    see this module addition's own top-of-file docstring, "Cascading
    invalidation", for the full rationale on why this must happen
    automatically and unconditionally, not as an opt-in behavior.

    params: a plain, JSON-serializable dict of EVERY parameter that
        materially affects this step's computed result (e.g. for "hvg":
        {"n_top_genes": 2000, "batch_key": "sample"}). Callers (the
        workspace layer) are responsible for including every parameter
        that matters here -- an omitted parameter would let a real
        change go undetected by check_downstream_step_current() below,
        silently serving a stale cached result. This is a genuine,
        ongoing engineering discipline requirement (see this module
        addition's own top-of-file docstring) each NEW parameter added
        to a step in the future must also be added here.

    Raises ValueError if step_name isn't a recognized member of
    DOWNSTREAM_STEP_ORDER (a real usage error from the calling UI code,
    not something to silently ignore or guess at).
    """
    from datetime import datetime

    if step_name not in DOWNSTREAM_STEP_ORDER:
        raise ValueError(f"Unknown downstream step: {step_name!r} -- must be one of {DOWNSTREAM_STEP_ORDER}")

    recipe = load_downstream_recipe(project_name)
    recipe[step_name] = {
        "params": params,
        "computed_at": datetime.now().isoformat(timespec="seconds"),
    }

    step_index = DOWNSTREAM_STEP_ORDER.index(step_name)
    for later_step in DOWNSTREAM_STEP_ORDER[step_index + 1:]:
        recipe.pop(later_step, None)

    _save_downstream_recipe_file(project_name, recipe)


def get_downstream_step_recipe(project_name, step_name):
    "Return step_name's saved {'params': {...}, 'computed_at': str} dict, or None if never saved (or invalidated by an earlier step being re-run since)."
    recipe = load_downstream_recipe(project_name)
    return recipe.get(step_name)


def check_downstream_step_current(project_name, step_name, current_params):
    """
    The authoritative "is it safe to reuse the cached result for this
    step, or must it be recomputed?" check.

    Returns True ONLY if step_name has a saved recipe AND its saved
    params dict is EXACTLY EQUAL (plain Python dict equality) to
    current_params -- i.e. every tracked parameter matches precisely,
    with no partial-match/fuzzy-match logic at all. Returns False
    otherwise (no recipe saved yet, OR any parameter differs, OR this
    step was cascade-invalidated by an earlier step being re-run).

    This function makes NO decision about what to DO with a False
    result -- per this module addition's own top-of-file docstring,
    the workspace layer must surface a mismatch explicitly to the user
    rather than silently recomputing OR silently reusing a stale cache.
    """
    saved = get_downstream_step_recipe(project_name, step_name)
    if saved is None:
        return False
    return saved.get("params") == current_params


def clear_downstream_step_and_after(project_name, step_name):
    """
    Explicit, MANUAL cascade-invalidation of step_name and every step
    after it in DOWNSTREAM_STEP_ORDER -- exposed for a "force re-run
    from this step onward" UI action, distinct from the AUTOMATIC
    cascade that save_downstream_step_recipe() already performs
    whenever a step is genuinely recomputed with new/different params.
    Useful e.g. if a user wants to force a fresh run with the EXACT
    SAME parameters as before (which check_downstream_step_current()
    would otherwise correctly report as "still current" and thus skip).

    Raises ValueError if step_name isn't a recognized step.
    """
    if step_name not in DOWNSTREAM_STEP_ORDER:
        raise ValueError(f"Unknown downstream step: {step_name!r} -- must be one of {DOWNSTREAM_STEP_ORDER}")
    recipe = load_downstream_recipe(project_name)
    step_index = DOWNSTREAM_STEP_ORDER.index(step_name)
    for s in DOWNSTREAM_STEP_ORDER[step_index:]:
        recipe.pop(s, None)
    _save_downstream_recipe_file(project_name, recipe)


# ---------------------------------------------------------------------------
# Sample eligibility for Phase 3 combination
# ---------------------------------------------------------------------------

def get_samples_with_completed_cellqc(project_name):
    """
    Return a list of sample names that have BOTH a completed STARsolo
    filtered matrix AND a completed Phase 2 Cell-level QC run (i.e.
    cell_qc_metrics.csv exists) -- the set of samples eligible to be
    combined into Phase 3's multi-sample AnnData object.

    A sample must have BOTH: STARsolo alignment produces the raw
    material Cell-level QC reads FROM, but it's Cell-level QC's own
    output (cell_qc_metrics.csv, containing doublet/ambient-RNA/
    adaptive-QC flags) that Phase 3's own
    sc_downstream_manager.load_sample_as_anndata() actually needs to
    attach as .obs columns -- a sample with completed alignment but no
    completed Cell-level QC run has no QC flags available to load, so
    it is correctly excluded here even though its raw matrix technically
    already exists on disk.

    --- Real bug fixed (2026-08-25): duplicated STARsolo output-path
        logic silently diverged from the authoritative source ---
    This function previously reconstructed the STARsolo filtered-matrix
    directory path via its OWN hand-written string concatenation
    (os.path.join(output_prefix + "Solo.out", "Gene", "filtered")) --
    a SEPARATE, independent reimplementation of the exact same path
    convention that singlecell_workspace.py's own render_cell_qc() and
    Step 6 alignment code already compute via
    starsolo_manager.filtered_counts_matrix_dir(output_prefix).

    This is precisely the class of bug this project has been bitten by
    before (see sc_cellqc_manager.py's own module docstring for an
    earlier, unrelated incident where duplicated GTF-parsing logic
    silently drifted apart between two modules) -- a real reported case
    confirmed this exact symptom: a project where alignment (Step 6)
    and Cell-level QC (Phase 2) had BOTH genuinely completed
    successfully for a sample (confirmed via real, correct mitochondrial-
    gene detection and QC metrics displayed on-screen) still had this
    function report ZERO eligible samples for Phase 3, because its own
    hardcoded path reconstruction here did not actually match whatever
    the real starsolo_manager.filtered_counts_matrix_dir() function
    itself computes.

    Fixed by calling starsolo_manager.filtered_counts_matrix_dir()
    DIRECTLY (the exact same function singlecell_workspace.py's own
    Step 6/Cell-level QC code already calls) instead of maintaining a
    second, independent copy of this path logic here -- this makes it
    structurally impossible for this eligibility check and the actual
    STARsolo output location to silently diverge again in the future,
    regardless of how starsolo_manager.py's own path convention is
    implemented or changes later.

    Imported locally (inside this function), not at module level, to
    avoid introducing a module-level import dependency into every
    OTHER function in this file that has nothing to do with STARsolo
    paths specifically -- mirrors this project's own established
    "import locally where a dependency is genuinely narrow/optional"
    convention (e.g. reference_manager.py's own local import of
    sc_cellqc_manager inside verify_preset_reference_mito_content()).
    """
    import starsolo_manager as star

    align_dir = starsolo_output_dir(project_name)
    if not os.path.isdir(align_dir):
        return []
    eligible = []
    for sample_name in sorted(os.listdir(align_dir)):
        sample_dir = os.path.join(align_dir, sample_name)
        if not os.path.isdir(sample_dir):
            continue
        cellqc_metrics_path = os.path.join(sample_dir, "cellqc", "cell_qc_metrics.csv")
        output_prefix = os.path.join(sample_dir, f"{sample_name}_")
        filtered_dir = star.filtered_counts_matrix_dir(output_prefix)
        if os.path.isfile(cellqc_metrics_path) and os.path.isdir(filtered_dir):
            eligible.append(sample_name)
    return eligible


def get_sample_spec_for_downstream(project_name, sample_name):
    """
    Build ONE entry of the sample_specs list
    sc_downstream_manager.load_and_combine_samples() expects, using this
    project's own known STARsolo/Cell-level-QC directory conventions --
    so the workspace layer never needs to hand-construct these paths
    itself.

    Returns {"sample_name": str, "starsolo_output_prefix": str,
    "cellqc_output_dir": str}.
    """
    sample_dir = os.path.join(starsolo_output_dir(project_name), sample_name)
    return {
        "sample_name": sample_name,
        "starsolo_output_prefix": os.path.join(sample_dir, f"{sample_name}_"),
        "cellqc_output_dir": os.path.join(sample_dir, "cellqc"),
    }
# ---------------------------------------------------------------------------
# Append to sc_project_manager.py -- gene-ID-mapping tracking for the SC
# DESeq2/Ontology bridge (v1 scope: bitr conversion + manual CSV upload
# only -- see sc_deseq2_workspace.py's module docstring for why the bulk
# pipeline's "layer 1 auto-derived-from-alignment-reference" gene map
# isn't ported here yet).
# ---------------------------------------------------------------------------

def sc_gene_symbol_map_path(project_name, export_name):
    "Where this pseudobulk export's gene_id -> gene_name mapping (from bitr conversion) is saved."
    return os.path.join(sc_deseq2_dir(project_name, export_name), "gene_symbol_map.csv")


def get_sc_gene_id_mapping_meta(project_name, export_name):
    "Return this export's saved gene-ID-mapping metadata (which conversion was run, when), or None."
    path = os.path.join(sc_deseq2_dir(project_name, export_name), "gene_id_mapping_meta.json")
    return atomic_io.read_json(path, default=None, on_corrupt="default")


def save_sc_gene_id_mapping_meta(project_name, export_name, meta):
    "Persist this export's gene-ID-mapping metadata."
    path = os.path.join(sc_deseq2_dir(project_name, export_name), "gene_id_mapping_meta.json")
    atomic_io.atomic_write_json(path, meta)


def sc_gene_id_mapping_work_dir(project_name, export_name):
    "Scratch directory for this export's gene-ID-mapping (bitr) intermediate files."
    return os.path.join(sc_deseq2_work_dir(project_name, export_name), "gene_id_mapping")
