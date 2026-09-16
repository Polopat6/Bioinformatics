"""
marker_library_manager.py

Custom tissue/cell-type marker set library for Phase 3.7 single-cell
annotation -- for users who already have curated marker gene knowledge
for their own tissue or organism (from literature, prior experience, or
a lab's own established panels), as a complement to (NOT a replacement
for) the eggNOG-mapper ortholog-bridging path in eggnog_manager.py.

--- Why this is a separate, simpler feature from eggNOG ---

eggNOG-mapper (eggnog_manager.py) solves "I only know HUMAN/MOUSE
marker genes, and need to find their non-model-organism equivalent" --
it requires a real, per-reference one-time setup (protein extraction +
a ~49GB database + an emapper.py run) before it produces anything
usable.

This module solves a DIFFERENT, much simpler problem: "I ALREADY KNOW
the correct marker genes for my own organism/tissue (e.g. from a prior
publication specifically on Austrofundulus limnaeus, or hard-won lab
knowledge), and want to save/reuse/import that list" -- no eggNOG
database, no ortholog bridging, no per-reference setup at all. The two
features are complementary, not alternatives: a user might build an
initial marker set FROM an eggNOG-derived symbol map, then save/refine
it here as their own curated library entry for reuse across future
projects with the same organism/tissue.

Every marker set saved here is consumed by
sc_downstream_manager.score_marker_gene_sets() using EXACTLY the same
{cell_type_label: [gene_list]} dict shape that function already
accepts and already resolves flexibly (by gene ID OR gene symbol,
against whatever a given project's own reference happens to use) --
this module adds NO new resolution logic of its own, purely
persistence, organization, and CSV import/export around marker sets a
user already has or is actively building.

--- Storage design ---

Marker sets are stored as a SHARED library (data/marker_libraries/),
independent of any single project -- mirroring the reasoning already
applied to the eggNOG database and reference genomes: a marker panel
for "T cell" or "killifish liver cell types" is inherently reusable
knowledge, not something that should be silently duplicated or
re-typed for every new project that happens to study the same
tissue/organism. Each library is a single named JSON file (e.g.
"killifish_liver_v1.json"), so a user can maintain several distinct,
independently-named panels (e.g. one per tissue, one per organism,
one from a specific publication) without them colliding or being
forced into one giant undifferentiated list.

A project can still layer its OWN one-off additions/edits on top of a
loaded shared library (see get_project_marker_sets()/
save_project_marker_overrides() below) without modifying the shared
library file itself -- the same "shared resource + per-project
override" pattern already used for e.g. mitochondrial gene resolution
in singlecell_workspace.py.
"""
import csv
import json
import os
import app_paths
from datetime import datetime


MARKER_LIBRARY_ROOT = app_paths.data_path("marker_libraries")


# ---------------------------------------------------------------------------
# Shared marker library CRUD
# ---------------------------------------------------------------------------

def list_marker_libraries():
    "Return a sorted list of existing shared marker library names (without the .json extension)."
    if not os.path.isdir(MARKER_LIBRARY_ROOT):
        return []
    return sorted([
        f[:-5] for f in os.listdir(MARKER_LIBRARY_ROOT)
        if f.endswith(".json") and os.path.isfile(os.path.join(MARKER_LIBRARY_ROOT, f))
    ])


def _library_path(library_name):
    return os.path.join(MARKER_LIBRARY_ROOT, f"{library_name}.json")


def load_marker_library(library_name):
    """
    Load a shared marker library by name.

    Returns a dict:
        {
            "marker_sets": {cell_type_label: [gene, ...], ...},
            "metadata": {"description": str, "organism": str,
                         "created_at": str, "last_updated": str},
        }
    Returns None if this library doesn't exist.
    """
    path = _library_path(library_name)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def save_marker_library(library_name, marker_sets, description="", organism=""):
    """
    Save (create or overwrite) a shared marker library.

    library_name: a filesystem-safe name (e.g. "killifish_liver_v1") --
        NOT sanitized here; the caller (UI layer) is expected to
        sanitize user-provided names the same way
        sc_project_manager.create_project() already does for project
        names, before calling this function.
    marker_sets: dict {cell_type_label: [gene_id_or_symbol, ...]}.
    description, organism: free-text metadata shown in the library
        picker UI to help a user identify the right library later
        (e.g. "Killifish liver cell types, from Smith et al. 2025" /
        "Austrofundulus limnaeus").

    Preserves the original "created_at" timestamp if this library
    already exists (only "last_updated" changes on a re-save) -- so
    re-saving an edited library doesn't lose its original creation
    date.
    """
    existing = load_marker_library(library_name)
    created_at = existing["metadata"]["created_at"] if existing else datetime.now().isoformat(timespec="seconds")

    data = {
        "marker_sets": marker_sets,
        "metadata": {
            "description": description,
            "organism": organism,
            "created_at": created_at,
            "last_updated": datetime.now().isoformat(timespec="seconds"),
        },
    }
    os.makedirs(MARKER_LIBRARY_ROOT, exist_ok=True)
    with open(_library_path(library_name), "w") as f:
        json.dump(data, f, indent=2)
    return data


def delete_marker_library(library_name):
    "Delete a shared marker library. Returns True if it existed and was deleted, False if it didn't exist."
    path = _library_path(library_name)
    if not os.path.isfile(path):
        return False
    os.remove(path)
    return True


# ---------------------------------------------------------------------------
# CSV import / export
# ---------------------------------------------------------------------------

def import_marker_sets_from_csv(csv_path):
    """
    Parse a user-uploaded CSV into a marker_sets dict, ready to pass
    directly to save_marker_library() or
    sc_downstream_manager.score_marker_gene_sets().

    Expected CSV format: two columns, header row required, named
    EXACTLY "cell_type" and "gene" (case-insensitive) -- ONE ROW PER
    (cell_type, gene) PAIR (long format), e.g.:

        cell_type,gene
        T_cell,CD3D
        T_cell,CD3E
        B_cell,MS4A1
        B_cell,CD79A

    This long format (rather than one-row-per-cell-type with a
    delimited gene list packed into a single cell) is deliberately
    chosen to be the simplest possible spreadsheet shape for a
    non-technical user to build directly in Excel/Google Sheets --
    every row is a single, unambiguous fact ("this gene is a marker for
    this cell type"), with no need to know or agree on an internal
    delimiter convention (comma vs. semicolon vs. pipe) for packing
    multiple genes into one cell.

    Returns (marker_sets: dict, warnings: list[str]) -- marker_sets is
    {cell_type_label: [gene, ...]} with genes de-duplicated per cell
    type but ORDER-PRESERVED (first occurrence order); warnings lists
    any skipped rows (e.g. missing cell_type or gene value) so the
    caller can surface exactly what was silently dropped, rather than
    the row count just not matching what the user expects with no
    explanation.

    Raises ValueError if the CSV doesn't have the required column
    names at all (a clear, immediate failure rather than silently
    returning an empty/wrong result from a misnamed-column file).
    """
    marker_sets = {}
    warnings = []

    with open(csv_path, newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("This CSV file appears to be empty.")

        normalized_fields = {name.strip().lower(): name for name in reader.fieldnames}
        if "cell_type" not in normalized_fields or "gene" not in normalized_fields:
            raise ValueError(
                "This CSV must have columns named exactly 'cell_type' and 'gene' "
                f"(found columns: {reader.fieldnames})."
            )
        cell_type_col = normalized_fields["cell_type"]
        gene_col = normalized_fields["gene"]

        for i, row in enumerate(reader, start=2):  # start=2: row 1 is the header
            cell_type = (row.get(cell_type_col) or "").strip()
            gene = (row.get(gene_col) or "").strip()
            if not cell_type or not gene:
                warnings.append(f"Row {i}: skipped (missing cell_type or gene value).")
                continue
            marker_sets.setdefault(cell_type, [])
            if gene not in marker_sets[cell_type]:
                marker_sets[cell_type].append(gene)

    return marker_sets, warnings


def export_marker_sets_to_csv(marker_sets, dest_path):
    """
    Write a marker_sets dict back out to the same long-format CSV
    convention import_marker_sets_from_csv() reads -- so a library
    built/edited in the UI can be downloaded, shared with a colleague,
    or archived outside this app, and later re-imported unchanged.

    Returns dest_path.
    """
    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
    with open(dest_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_type", "gene"])
        for cell_type, genes in marker_sets.items():
            for gene in genes:
                writer.writerow([cell_type, gene])
    return dest_path


# ---------------------------------------------------------------------------
# Per-project overrides on top of a shared library
# ---------------------------------------------------------------------------
#
# Mirrors the exact "shared resource + per-project override" pattern
# already used for mitochondrial gene resolution in
# singlecell_workspace.py's reference_choice dict -- a project can add,
# remove, or edit individual cell-type marker sets on top of a loaded
# shared library, WITHOUT modifying that shared library file itself
# (which other projects may also be using unmodified).

def get_project_marker_sets(project_info, base_library_name=None):
    """
    Resolve the EFFECTIVE marker_sets dict for one project: starts from
    a shared library (if base_library_name is given and exists), then
    applies any project-specific overrides already saved onto
    project_info (see save_project_marker_overrides() below) on top --
    an override entry for a given cell_type_label REPLACES that
    cell_type's gene list entirely (not merged/unioned with the base
    library's own list for that same label), so a project can cleanly
    correct/replace a single cell type's marker list without needing
    to also carry along every other untouched cell type's overrides.

    project_info: this project's own loaded info dict (e.g. from
        sc_project_manager.load_info(project)) -- reads its
        "marker_set_overrides" key if present (added by
        save_project_marker_overrides() below); does not require any
        particular caller to have set this key already (defaults to no
        overrides at all if absent).
    base_library_name: name of a shared library to start from, or None
        to start from an empty dict (a project with entirely its own,
        library-independent marker sets).

    Returns the effective marker_sets dict, ready to pass directly to
    sc_downstream_manager.score_marker_gene_sets().
    """
    effective = {}
    if base_library_name:
        library = load_marker_library(base_library_name)
        if library:
            effective.update(library["marker_sets"])

    overrides = (project_info or {}).get("marker_set_overrides", {})
    effective.update(overrides)

    return effective


def save_project_marker_overrides(project_info, marker_set_overrides, base_library_name=None):
    """
    Persist project-specific marker set overrides (and optionally which
    shared library they're layered on top of) onto a project's own
    info dict -- mirrors sc_project_manager.py's own
    save_reference_choice()/save_chemistry_choice() convention of
    storing a labeled sub-dict under a clear top-level key.

    project_info: this project's own loaded info dict -- MODIFIED IN
        PLACE and also returned, so the caller can immediately pass it
        to sc_project_manager.save_info(project, project_info) to
        actually persist it to disk (this function does not itself
        write to disk -- it has no dependency on sc_project_manager.py,
        keeping this module usable independently of the single-cell
        project system specifically).
    marker_set_overrides: dict {cell_type_label: [gene, ...]} -- the
        project's own additions/edits, layered on top of
        base_library_name (if any) by get_project_marker_sets() above.
    base_library_name: name of the shared library this project's
        overrides are based on, or None if this project uses no shared
        library at all (entirely its own marker sets).

    Returns the modified project_info dict.
    """
    project_info["marker_set_overrides"] = marker_set_overrides
    project_info["marker_set_base_library"] = base_library_name
    return project_info
