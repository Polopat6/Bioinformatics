"""
single_cell/sc_genetic_demux_manager.py

Support for CONSUMING an externally/precomputed genetic-demultiplexing
result (e.g. demuxlet/popscle output) to split ONE pooled, multi-donor
STARsolo sample into N donor-labeled "samples" that behave, to every
other part of this pipeline, exactly like ordinary per-donor STARsolo
output.

--- Why this module exists (2026-08-23) ---
Some public datasets (the motivating case: Kang et al. 2018, GSE96583 --
8 lupus donors' PBMCs pooled per condition and sequenced as just 2 10x
libraries -- "batch 2 control" / "batch 2 stim") pool MULTIPLE donors
into a SINGLE 10x library, and rely on genetic multiplexing (demuxlet,
using each donor's own SNP genotypes) to recover per-cell donor identity
AFTER sequencing. Pulling such an accession's raw FASTQ from NCBI/SRA
and running it through this pipeline's normal Phase 1 (ingestion ->
STARsolo) produces exactly ONE aligned sample per accession -- e.g. 2
samples total for Kang (one per condition) -- with NO donor-level
replication at all, which defeats the actual purpose of Steps 9/10
(pseudobulk DESeq2, compositional analysis), both of which require
real per-sample (here: per-donor) replication to run meaningfully.

This module does NOT implement genetic demultiplexing itself (running
demuxlet/popscle requires each donor's own reference genotype VCF,
which is privacy-sensitive patient data and is NOT part of any public
GEO deposit for datasets like this -- see this project's own recorded
design note on this). Instead, it CONSUMES an already-computed,
publicly available barcode-to-donor assignment table (for Kang: GEO's
own GSE96583_batch2.total.tsne.df.tsv.gz supplementary file, or the
identical data via muscData's Kang18_8vs8() accessor) and uses it to
split one pooled STARsolo sample's filtered matrix into N real,
independent, donor-labeled samples -- each written out as an ordinary
10x-format MTX triplet, indistinguishable from genuine per-donor
STARsolo output to every downstream step (Phase 2 Cell-level QC, Phase
3 combine/annotation/pseudobulk/compositional).

This same mechanism is intentionally dataset-agnostic: ANY externally-
sourced barcode-to-donor/condition assignment table (a future
genotype-free souporcell/vireo run, a CITE-seq hashtag assignment file,
etc.) can be consumed the same way, via the same functions -- this
module has no Kang-specific logic baked into its core split/match
functions; only the convenience wrapper at the bottom (documented as
such) is shaped around Kang's own specific column conventions.

--- Where donor/condition labels actually end up (2026-08-23) ---
Each split-off donor sample gets a small "cellqc/external_cell_metadata.csv"
file written into what will become that sample's own Cell-level QC
output directory -- consumed by sc_downstream_manager.py's
load_sample_as_anndata(), which merges it into .obs alongside (never
overwriting) Phase 2's own real cell_qc_metrics.csv, once a user runs
REAL Cell-level QC on that split sample through the existing, unmodified
Phase 2 UI. See that function's own docstring, "External per-cell
metadata extension point", for the full merge rationale -- nothing in
THIS module writes to cell_qc_metrics.csv itself, or bypasses Phase 2
in any way; a split sample must still go through genuine Cell-level QC
like any other sample before it's eligible for Phase 3 combination (see
sc_project_manager.get_samples_with_completed_cellqc()).

--- SoupX ambient-RNA correction is NOT meaningfully available for a
    split sample (2026-08-23) ---
A genetic demultiplexing assignment (like demuxlet's) is only computed
for CALLED cells that passed cell-calling in the first place -- empty
droplets in STARsolo's own raw/unfiltered matrix have no donor identity
at all (there was never enough signal to call one). This means a split
sample has no valid per-donor "raw" matrix to hand SoupX (which
specifically needs the raw/unfiltered matrix's empty-droplet population
to estimate the ambient RNA profile) -- DecontX (Phase 2's default,
which only needs the filtered matrix) remains fully valid and is the
recommended choice for any split sample. This module does not attempt
to fabricate or split a "raw" matrix for this reason; it only ever
produces a filtered/-shaped output.

--- Barcode suffix handling (2026-08-23) ---
STARsolo's own barcodes.tsv output and a third-party barcode table (like
GEO's) may or may not share the same "-1" GEM-well suffix convention --
this is NOT guaranteed to match by construction (confirmed: Cell
Ranger's own convention always appends "-N"; STARsolo's default output
format has varied across versions/settings). match_barcodes() below
tries an exact match first, and if that yields a poor match rate,
automatically retries after normalizing both sides' barcodes (stripping
any trailing "-<digits>" suffix) -- but ALWAYS reports which strategy
was used and the resulting match rate, rather than silently guessing.
A caller seeing an unexpectedly low match rate even after normalization
should treat that as a real red flag (e.g. wrong external table matched
to the wrong pooled sample), not proceed blindly.
"""
import gzip
import os
import re

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

# Must exactly match sc_downstream_manager.EXTERNAL_CELL_METADATA_FILENAME.
# Duplicated here as a plain string constant (rather than importing
# sc_downstream_manager at module level) to keep this module's only real
# dependency on that one narrow, well-documented filename convention --
# not a full cross-module import coupling -- matching this project's own
# established "dependency injection over hidden coupling" convention
# (e.g. eggnog_manager.download_eggnog_database's own
# ensure_shared_resource_fn parameter). If this filename is ever changed
# in sc_downstream_manager.py, it must be changed here too.
EXTERNAL_METADATA_FILENAME = "external_cell_metadata.csv"


# ---------------------------------------------------------------------------
# Low-level: reading/writing a 10x-format MTX triplet directly
# ---------------------------------------------------------------------------

def read_starsolo_filtered_matrix(filtered_dir):
    """
    Load a STARsolo (or Cell Ranger-compatible) filtered/ directory's
    raw MTX triplet directly -- deliberately NOT going through
    scanpy/anndata here, since this module only needs to SUBSET and
    RE-WRITE the raw matrix, not analyze it.

    Returns (matrix, barcodes, features_df):
        matrix: scipy sparse matrix, shape (n_genes, n_cells) -- the
            SAME raw orientation STARsolo's own matrix.mtx uses (genes
            as rows, cells as columns) -- NOT transposed, so this can
            be written back out (see write_starsolo_style_filtered_matrix
            below) in the exact same convention real STARsolo output
            uses, and so sc_downstream_manager.load_sample_as_anndata's
            own `sc.read_mtx(...).T` continues to work identically on
            split output as on genuine STARsolo output.
        barcodes: list[str], length n_cells, in the SAME column order
            as matrix.
        features_df: DataFrame with columns "gene_id", "gene_symbol"
            (and "feature_type" if a 3rd column was present in the
            original features.tsv), in the SAME row order as matrix.

    Raises FileNotFoundError if any of the three required files
    (barcodes/features/matrix, gzipped or not) can't be found.
    """
    def _find(*candidates):
        for c in candidates:
            path = os.path.join(filtered_dir, c)
            if os.path.isfile(path):
                return path
        return None

    barcodes_path = _find("barcodes.tsv.gz", "barcodes.tsv")
    features_path = _find("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
    matrix_path = _find("matrix.mtx.gz", "matrix.mtx")

    missing = [
        name for name, path in
        (("barcodes", barcodes_path), ("features", features_path), ("matrix", matrix_path))
        if path is None
    ]
    if missing:
        raise FileNotFoundError(
            f"Could not find required 10x-format file(s) ({', '.join(missing)}) in {filtered_dir}"
        )

    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].tolist()

    features_raw = pd.read_csv(features_path, header=None, sep="\t")
    features_df = pd.DataFrame({"gene_id": features_raw[0]})
    features_df["gene_symbol"] = features_raw[1] if features_raw.shape[1] > 1 else features_raw[0]
    if features_raw.shape[1] > 2:
        features_df["feature_type"] = features_raw[2]

    matrix = sio.mmread(matrix_path).tocsc()

    if matrix.shape[1] != len(barcodes) or matrix.shape[0] != len(features_df):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match barcodes ({len(barcodes)}) / "
            f"features ({len(features_df)}) counts in {filtered_dir} -- files may be "
            f"corrupted or mismatched."
        )

    return matrix, barcodes, features_df


def write_starsolo_style_filtered_matrix(dest_dir, matrix, barcodes, features_df):
    """
    Write a (genes x cells) sparse matrix + barcode list + feature table
    out as a standard 10x-format MTX triplet
    (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz) at dest_dir --
    indistinguishable, to any downstream reader (including this
    pipeline's own sc_downstream_manager.load_sample_as_anndata()),
    from genuine STARsolo output.

    matrix: scipy sparse matrix, shape (n_genes, n_cells) -- same raw
        orientation as read_starsolo_filtered_matrix()'s own return
        value; NOT transposed.
    barcodes: list[str], length n_cells, matching matrix's column order.
    features_df: DataFrame with "gene_id" and "gene_symbol" columns
        (and optionally "feature_type"), matching matrix's row order.

    Overwrites any existing files at dest_dir with the same names.
    """
    os.makedirs(dest_dir, exist_ok=True)

    barcodes_path = os.path.join(dest_dir, "barcodes.tsv.gz")
    with gzip.open(barcodes_path, "wt") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    features_path = os.path.join(dest_dir, "features.tsv.gz")
    with gzip.open(features_path, "wt") as f:
        for _, row in features_df.iterrows():
            cols = [str(row["gene_id"]), str(row["gene_symbol"])]
            if "feature_type" in features_df.columns:
                cols.append(str(row["feature_type"]))
            else:
                cols.append("Gene Expression")
            f.write("\t".join(cols) + "\n")

    matrix_path_gz = os.path.join(dest_dir, "matrix.mtx.gz")
    # scipy.io.mmwrite doesn't gzip directly -- write to an uncompressed
    # temp path in the same directory, then gzip it in place, so a
    # crash mid-write never leaves a half-written matrix.mtx.gz sitting
    # at the final filename.
    tmp_path = os.path.join(dest_dir, "_matrix_tmp.mtx")
    sio.mmwrite(tmp_path, sparse.coo_matrix(matrix))
    with open(tmp_path, "rb") as f_in, gzip.open(matrix_path_gz, "wb") as f_out:
        f_out.writelines(f_in)
    os.remove(tmp_path)


# ---------------------------------------------------------------------------
# Loading an external barcode-level assignment table
# ---------------------------------------------------------------------------

def load_external_barcode_table(table_path, sep=None):
    """
    Load an external barcode-level metadata table (e.g. a downloaded
    GEO supplementary file, or any similarly-shaped barcode-to-donor/
    condition assignment table) -- deliberately flexible about exact
    formatting, since this file's real-world shape (delimiter, column
    naming/capitalization, presence of an index column) is often not
    100% documented ahead of time and should never be silently assumed.

    sep: explicit delimiter to use ("\\t" or ","). If None (the
        default), auto-detects by trying tab first, falling back to
        comma if that produces only a single column (a common sign the
        wrong delimiter was tried).

    Returns the raw DataFrame, completely unmodified otherwise (no
    column renaming/filtering happens here -- see
    standardize_barcode_table() below for that) -- so a caller can
    always inspect df.columns themselves first if the table's real
    shape is still uncertain.

    Raises FileNotFoundError if table_path doesn't exist.
    """
    if not os.path.isfile(table_path):
        raise FileNotFoundError(f"External barcode table not found: {table_path}")

    if sep is not None:
        return pd.read_csv(table_path, sep=sep)

    df = pd.read_csv(table_path, sep="\t")
    if df.shape[1] <= 1:
        df = pd.read_csv(table_path, sep=",")
    return df


def standardize_barcode_table(df, barcode_col=None, donor_col=None, condition_col=None,
                               extra_cols=None):
    """
    Normalize an arbitrary external barcode-level table into a
    predictable shape: columns "barcode", "donor", (optionally)
    "condition", plus any requested extra_cols carried through
    unchanged -- so every other function in this module can rely on
    consistent column names regardless of the source file's own real
    header naming.

    barcode_col: name of the column (or the DataFrame's own index, if
        None and no column named "barcode"/"Unnamed: 0" is found)
        holding each row's cell barcode. If None, this function tries,
        in order: a column literally named "barcode", then "Unnamed: 0"
        (pandas' own default name for an unlabeled first CSV column,
        common when a table was exported with its barcode as a row
        index), then finally falls back to the DataFrame's own index
        if it looks barcode-like (see _looks_like_barcode_series below).
    donor_col: name of the column holding donor/individual identity
        (e.g. "ind" for the Kang et al. 2018 table). REQUIRED -- there
        is no universal default column name for this across possible
        source tables, so this must always be specified explicitly by
        the caller rather than guessed.
    condition_col: optional name of the column holding experimental
        condition (e.g. "stim"). May be None if the source table has no
        condition column at all (e.g. a dataset with only one
        condition, where condition is instead implied by which pooled
        accession/sample this table applies to).
    extra_cols: optional list of additional column names to carry
        through unchanged (e.g. ["cell", "cluster"] for Kang -- the
        original paper's own cell-type/cluster annotation, kept as an
        optional REFERENCE-ONLY cross-check column, never used
        functionally by this pipeline itself).

    Returns a DataFrame with columns ["barcode", "donor"] plus
    "condition" (if condition_col was given) plus every requested
    extra_cols -- one row per barcode in the ORIGINAL table (not yet
    matched/subset against any particular pooled sample's own
    barcodes -- see match_barcodes() for that).

    Raises ValueError if donor_col isn't found in df, or if barcode_col
    couldn't be resolved by any of the above strategies.
    """
    if donor_col not in df.columns:
        raise ValueError(f"donor_col '{donor_col}' not found in table columns: {list(df.columns)}")
    if condition_col is not None and condition_col not in df.columns:
        raise ValueError(f"condition_col '{condition_col}' not found in table columns: {list(df.columns)}")

    resolved_barcode_col = barcode_col
    barcode_series = None
    if resolved_barcode_col is not None:
        if resolved_barcode_col not in df.columns:
            raise ValueError(f"barcode_col '{resolved_barcode_col}' not found in table columns: {list(df.columns)}")
        barcode_series = df[resolved_barcode_col]
    else:
        for candidate in ("barcode", "Barcode", "Unnamed: 0"):
            if candidate in df.columns:
                barcode_series = df[candidate]
                break
        if barcode_series is None and _looks_like_barcode_series(df.index):
            barcode_series = df.index.to_series().reset_index(drop=True)
        if barcode_series is None:
            raise ValueError(
                "Could not resolve a barcode column automatically -- no column named "
                "'barcode'/'Barcode'/'Unnamed: 0' was found, and the table's own index "
                "doesn't look like cell barcodes either. Pass barcode_col explicitly."
            )

    result = pd.DataFrame({
        "barcode": barcode_series.astype(str).values,
        "donor": df[donor_col].astype(str).values,
    })
    if condition_col is not None:
        result["condition"] = df[condition_col].astype(str).values
    if extra_cols:
        for col in extra_cols:
            if col in df.columns:
                result[col] = df[col].values
    return result


_BARCODE_LIKE_PATTERN = re.compile(r"^[ACGTN]{8,20}(-\d+)?$", re.IGNORECASE)


def _looks_like_barcode_series(index_like, sample_size=20):
    "Heuristic: does this Index/Series look like it's made of DNA cell barcodes (ACGT, 8-20bp, optional -N suffix)? Used only to decide whether a table's own row index is usable as its barcode column when no explicit barcode column can be found -- never used to VALIDATE real barcode content elsewhere."
    values = list(index_like)[:sample_size]
    if not values:
        return False
    matches = sum(1 for v in values if _BARCODE_LIKE_PATTERN.match(str(v)))
    return matches / len(values) >= 0.8


# ---------------------------------------------------------------------------
# Barcode matching (pooled STARsolo barcodes <-> external table barcodes)
# ---------------------------------------------------------------------------

_TRAILING_SUFFIX_PATTERN = re.compile(r"-\d+$")


def _strip_suffix(barcode):
    return _TRAILING_SUFFIX_PATTERN.sub("", barcode)


def match_barcodes(pooled_barcodes, standardized_table_df, min_acceptable_match_rate=0.5):
    """
    Match a pooled sample's own STARsolo barcodes against an external
    table's barcodes, trying an exact match first and automatically
    falling back to suffix-normalized matching if the exact match rate
    is poor -- see this module's own docstring, "Barcode suffix
    handling", for the full rationale on why this can't be assumed to
    match by construction.

    pooled_barcodes: list[str], from read_starsolo_filtered_matrix().
    standardized_table_df: output of standardize_barcode_table().
    min_acceptable_match_rate: if the EXACT match rate is below this
        threshold, automatically retry with suffix-normalized barcodes
        on both sides before giving up.

    Returns a dict:
        {
            "strategy_used": "exact" | "suffix_normalized",
            "match_rate": float,             # fraction of pooled_barcodes matched
            "n_matched": int,
            "n_unmatched": int,
            "matched_donor_map": {pooled_barcode: donor_row_dict, ...},
                # donor_row_dict is the FULL matched row from
                # standardized_table_df (as a dict) -- includes "donor",
                # "condition" (if present), and any extra_cols.
            "unmatched_barcodes": [str, ...],  # up to 20, for display
        }

    Does NOT raise on a poor match rate -- always returns its result
    and lets the caller (typically preview_donor_split(), and
    ultimately the UI layer) decide whether to proceed, since a
    genuinely low match rate might still be intentional/expected in
    some edge case, and this module's own job is to report clearly,
    not to silently block.
    """
    table_by_exact = {row["barcode"]: row for row in standardized_table_df.to_dict("records")}

    def _try_match(barcode_transform):
        matched = {}
        lookup = {barcode_transform(k): v for k, v in table_by_exact.items()}
        for bc in pooled_barcodes:
            key = barcode_transform(bc)
            if key in lookup:
                matched[bc] = lookup[key]
        return matched

    exact_matched = _try_match(lambda b: b)
    exact_rate = len(exact_matched) / len(pooled_barcodes) if pooled_barcodes else 0.0

    if exact_rate >= min_acceptable_match_rate:
        matched = exact_matched
        strategy = "exact"
        match_rate = exact_rate
    else:
        suffix_matched = _try_match(_strip_suffix)
        suffix_rate = len(suffix_matched) / len(pooled_barcodes) if pooled_barcodes else 0.0
        if suffix_rate > exact_rate:
            matched = suffix_matched
            strategy = "suffix_normalized"
            match_rate = suffix_rate
        else:
            matched = exact_matched
            strategy = "exact"
            match_rate = exact_rate

    unmatched = [bc for bc in pooled_barcodes if bc not in matched]

    return {
        "strategy_used": strategy,
        "match_rate": round(match_rate, 4),
        "n_matched": len(matched),
        "n_unmatched": len(unmatched),
        "matched_donor_map": matched,
        "unmatched_barcodes": unmatched[:20],
    }


# ---------------------------------------------------------------------------
# Preview (dry-run) + actual split
# ---------------------------------------------------------------------------

def preview_donor_split(filtered_dir, standardized_table_df, singlets_only=True,
                         singlet_flag_col=None, singlet_flag_value="singlet",
                         min_cells_per_donor=1):
    """
    Dry-run preview of how a pooled sample would split by donor --
    reads the pooled matrix's barcodes ONLY (does not load the full
    matrix into memory or write anything), matches them against
    standardized_table_df, and reports per-donor cell counts, unmatched
    barcodes, and (if requested) how many cells would be excluded as
    non-singlets -- all BEFORE committing to the actual split, mirroring
    this project's own established "preview before run" pattern (e.g.
    sc_downstream_manager.get_pseudobulk_group_sizes()).

    singlets_only: if True (the default) and singlet_flag_col is
        provided, cells whose singlet_flag_col value != singlet_flag_value
        are excluded from the per-donor counts below and reported
        separately -- e.g. dropping demuxlet-flagged "doublet"/
        "ambiguous" cells before ever creating a donor split for them
        (a GENETIC multiplet call, entirely independent from and
        complementary to Phase 2's own EXPRESSION-based doublet
        detection that will still run separately on each resulting
        donor sample).
    singlet_flag_col: which column in standardized_table_df holds the
        singlet/doublet/ambiguous flag (e.g. "multiplets" for Kang) --
        required if singlets_only=True.

    Returns a dict:
        {
            "match_info": {...},               # see match_barcodes()
            "per_donor_counts": {donor: n_cells, ...},
            "n_excluded_non_singlet": int,      # 0 if singlets_only=False
            "n_cells_total_in_matrix": int,
        }
    """
    _, pooled_barcodes, _ = read_starsolo_filtered_matrix(filtered_dir)
    match_info = match_barcodes(pooled_barcodes, standardized_table_df)

    per_donor_counts = {}
    n_excluded_non_singlet = 0
    for bc, row in match_info["matched_donor_map"].items():
        if singlets_only and singlet_flag_col:
            flag_value = row.get(singlet_flag_col)
            if flag_value != singlet_flag_value:
                n_excluded_non_singlet += 1
                continue
        donor = row["donor"]
        per_donor_counts[donor] = per_donor_counts.get(donor, 0) + 1

    per_donor_counts = {
        d: n for d, n in per_donor_counts.items() if n >= min_cells_per_donor
    }

    return {
        "match_info": match_info,
        "per_donor_counts": per_donor_counts,
        "n_excluded_non_singlet": n_excluded_non_singlet,
        "n_cells_total_in_matrix": len(pooled_barcodes),
    }


def split_pooled_sample_by_donor(filtered_dir, standardized_table_df, dest_root_dir,
                                  sample_name_fn, singlets_only=True,
                                  singlet_flag_col=None, singlet_flag_value="singlet",
                                  min_cells_per_donor=1):
    """
    Actually perform the donor split -- reads the full pooled matrix
    ONCE, subsets it per donor, and writes each donor's own
    STARsolo-style filtered/ directory plus its
    "cellqc/{EXTERNAL_METADATA_FILENAME}" barcode-level metadata file.

    filtered_dir: the pooled sample's own STARsolo filtered/ directory.
    standardized_table_df: output of standardize_barcode_table().
    dest_root_dir: root directory under which each donor's own
        "{sample_name}/{sample_name}_Solo.out/Gene/filtered/" directory
        tree is created -- pass sc_project_manager.starsolo_output_dir(
        project_name) directly so each split donor "sample" lands
        exactly where sc_project_manager.get_samples_with_completed_cellqc()
        and get_sample_spec_for_downstream() already expect to find any
        ordinary sample.
    sample_name_fn: callable(donor: str) -> str, building each output
        sample's name (e.g. lambda donor: f"{donor}_ctrl" for a control-
        condition pooled accession). Kept as an injected callable
        (rather than a fixed naming template) so the caller controls
        exactly how condition/prefix information is folded into the
        final sample name, without this function needing to know
        anything about condition naming conventions itself.
    singlets_only / singlet_flag_col / singlet_flag_value /
        min_cells_per_donor: see preview_donor_split()'s own docstring
        -- applied identically here.

    Returns a dict:
        {
            "per_donor_summary": [
                {"donor": str, "sample_name": str, "n_cells": int,
                 "output_dir": str}, ...
            ],
            "match_info": {...},               # see match_barcodes()
            "n_excluded_non_singlet": int,
        }

    Does NOT touch any project-level metadata.csv or project_info.json
    -- registering each new sample's donor/condition into the project's
    own sample-level metadata is handled separately by
    append_sample_metadata_rows() below, kept as an explicit, separate
    step so this function's own responsibility stays narrowly scoped to
    "read one pooled matrix, write N donor matrices."
    """
    matrix, pooled_barcodes, features_df = read_starsolo_filtered_matrix(filtered_dir)
    match_info = match_barcodes(pooled_barcodes, standardized_table_df)

    barcode_to_col_index = {bc: i for i, bc in enumerate(pooled_barcodes)}

    donor_to_barcodes = {}
    donor_to_rows = {}
    n_excluded_non_singlet = 0
    for bc, row in match_info["matched_donor_map"].items():
        if singlets_only and singlet_flag_col:
            flag_value = row.get(singlet_flag_col)
            if flag_value != singlet_flag_value:
                n_excluded_non_singlet += 1
                continue
        donor = row["donor"]
        donor_to_barcodes.setdefault(donor, []).append(bc)
        # IMPORTANT: "barcode": bc must come AFTER **row in this dict
        # literal, not before -- row (a matched record from
        # standardized_table_df) already carries its OWN "barcode" key,
        # populated from the EXTERNAL table's own barcode value (which,
        # per this module's own "Barcode suffix handling" docstring
        # section, may differ from the pooled matrix's real barcode by
        # exactly the kind of "-1" suffix mismatch match_barcodes() is
        # designed to tolerate). In a Python dict literal, later keys
        # win when merging -- {"barcode": bc, **row} would therefore let
        # row's own (possibly suffix-mismatched) "barcode" value SILENTLY
        # overwrite the correct, real pooled-matrix barcode (bc),
        # producing an external_cell_metadata.csv whose "barcode" column
        # does not actually match its own sibling matrix's real
        # barcodes.tsv at all -- confirmed via direct testing to be a
        # real bug, not just a hypothetical concern. {**row, "barcode": bc}
        # (row first, explicit barcode key last) is the correct order.
        donor_to_rows.setdefault(donor, []).append({**row, "barcode": bc})

    per_donor_summary = []
    for donor, donor_barcodes in donor_to_barcodes.items():
        if len(donor_barcodes) < min_cells_per_donor:
            continue

        sample_name = sample_name_fn(donor)
        col_indices = [barcode_to_col_index[bc] for bc in donor_barcodes]
        donor_matrix = matrix[:, col_indices]

        output_dir = os.path.join(
            dest_root_dir, sample_name, f"{sample_name}_Solo.out", "Gene", "filtered",
        )
        write_starsolo_style_filtered_matrix(output_dir, donor_matrix, donor_barcodes, features_df)

        cellqc_dir = os.path.join(dest_root_dir, sample_name, "cellqc")
        os.makedirs(cellqc_dir, exist_ok=True)
        external_metadata_df = pd.DataFrame(donor_to_rows[donor])
        external_metadata_df.to_csv(
            os.path.join(cellqc_dir, EXTERNAL_METADATA_FILENAME), index=False,
        )

        per_donor_summary.append({
            "donor": donor,
            "sample_name": sample_name,
            "n_cells": len(donor_barcodes),
            "output_dir": os.path.dirname(os.path.dirname(os.path.dirname(output_dir))),
        })

    return {
        "per_donor_summary": per_donor_summary,
        "match_info": match_info,
        "n_excluded_non_singlet": n_excluded_non_singlet,
    }


# ---------------------------------------------------------------------------
# Sample-level metadata.csv registration
# ---------------------------------------------------------------------------

def append_sample_metadata_rows(metadata_csv_path, rows, sample_col="sample"):
    """
    Append (or update, if a sample_name already has a row) new rows to
    a project's sample-level metadata.csv -- matching
    singlecell_ingestion_manager.py's own "one row per SAMPLE, not per
    cell" metadata.csv convention exactly, so every split donor sample
    is registered as a completely ordinary sample from the rest of the
    pipeline's perspective.

    rows: list of dicts, each MUST include sample_col (the sample name
        key, matching each entry's own directory name under
        starsolo_output_dir()) -- e.g.
        [{"sample": "donor101_ctrl", "donor": "donor101", "condition": "ctrl"}, ...]

    If metadata_csv_path doesn't exist yet, it's created fresh with
    exactly the columns present in rows. If it already exists, new
    columns present in rows but not in the existing file are added
    (filled with NaN for pre-existing rows); rows whose sample_col value
    already exists in the file are UPDATED in place (not duplicated) --
    since re-running a split (e.g. after fixing a bad barcode-table
    path) should cleanly replace that sample's own metadata row, not
    accumulate duplicate rows for the same sample name.

    Returns the final combined DataFrame that was written.
    """
    new_df = pd.DataFrame(rows)
    if sample_col not in new_df.columns:
        raise ValueError(f"Every row must include a '{sample_col}' key.")

    if os.path.isfile(metadata_csv_path):
        existing_df = pd.read_csv(metadata_csv_path)
        if sample_col not in existing_df.columns:
            raise ValueError(
                f"Existing metadata.csv at {metadata_csv_path} has no '{sample_col}' column -- "
                f"cannot safely merge new rows into it."
            )
        existing_df = existing_df[~existing_df[sample_col].isin(new_df[sample_col])]
        combined_df = pd.concat([existing_df, new_df], ignore_index=True, sort=False)
    else:
        combined_df = new_df

    os.makedirs(os.path.dirname(metadata_csv_path) or ".", exist_ok=True)
    combined_df.to_csv(metadata_csv_path, index=False)
    return combined_df


# ---------------------------------------------------------------------------
# High-level orchestration (project-aware convenience wrapper)
# ---------------------------------------------------------------------------

def split_pooled_project_sample(project_name, pooled_sample_name, external_table_df,
                                 donor_col, condition_value, condition_col=None,
                                 barcode_col=None, extra_cols=None,
                                 singlets_only=True, singlet_flag_col=None,
                                 singlet_flag_value="singlet", min_cells_per_donor=1,
                                 sample_name_template="{donor}_{condition}"):
    """
    Project-aware convenience wrapper around split_pooled_sample_by_donor()
    + append_sample_metadata_rows() -- locates the pooled sample's own
    filtered/ directory via sc_project_manager's existing path
    conventions, performs the split, and registers each resulting donor
    sample into the project's own metadata.csv, all in one call.

    Imports sc_project_manager LOCALLY (function-level, not at module
    top) rather than as a hard top-of-file dependency -- this module's
    own core split/match/write logic above has NO dependency on
    sc_project_manager at all and is independently usable/testable;
    only this one convenience wrapper needs to know about this
    project's specific directory-layout conventions.

    project_name: the single-cell project this pooled sample belongs to.
    pooled_sample_name: the EXISTING sample name (already aligned via
        Step 6 STARsolo) representing the pooled, multi-donor library
        to split -- e.g. "batch2_control" or "batch2_stim".
    external_table_df: RAW external barcode table (e.g. from
        load_external_barcode_table()) -- standardized internally via
        standardize_barcode_table() using the donor_col/condition_col/
        barcode_col/extra_cols arguments below.
    donor_col: column in external_table_df holding donor identity
        (e.g. "ind" for Kang).
    condition_value: the FIXED condition label for every cell in THIS
        pooled sample (e.g. "ctrl" or "stim") -- used to build each
        split sample's own name (see sample_name_template) and written
        into each donor's own metadata row, even if condition_col is
        also present in external_table_df and could in principle be
        read per-cell. This is deliberately the caller-supplied,
        AUTHORITATIVE condition value for this whole pooled sample
        (matching the real-world Kang design, where an entire pooled
        10x library is ALWAYS one single condition -- see this
        module's own docstring) -- condition_col (if given) is
        additionally read and preserved per-cell purely as a
        cross-check/reference column, not as the source of truth.
    condition_col: optional column in external_table_df ALSO holding a
        per-cell condition label, carried through as a reference-only
        column (see condition_value above for why it isn't the
        authoritative source).
    barcode_col / extra_cols: see standardize_barcode_table().
    singlets_only / singlet_flag_col / singlet_flag_value /
        min_cells_per_donor: see preview_donor_split()/
        split_pooled_sample_by_donor()'s own docstrings.
    sample_name_template: a Python format string with "{donor}" and
        "{condition}" placeholders, used to build each split sample's
        final name (e.g. "{donor}_{condition}" -> "donor101_ctrl").

    Returns the same dict shape as split_pooled_sample_by_donor(),
    augmented with "metadata_df" (the full, updated project metadata.csv
    contents after registration).

    Raises FileNotFoundError if pooled_sample_name has no completed
    STARsolo alignment yet (no filtered/ directory found).
    """
    import sc_project_manager as scpm

    standardized_df = standardize_barcode_table(
        external_table_df, barcode_col=barcode_col, donor_col=donor_col,
        condition_col=condition_col, extra_cols=extra_cols,
    )

    pooled_spec = scpm.get_sample_spec_for_downstream(project_name, pooled_sample_name)
    filtered_dir = os.path.join(
        pooled_spec["starsolo_output_prefix"] + "Solo.out", "Gene", "filtered",
    )
    if not os.path.isdir(filtered_dir):
        raise FileNotFoundError(
            f"No completed STARsolo alignment found for pooled sample '{pooled_sample_name}' "
            f"in project '{project_name}' -- expected: {filtered_dir}"
        )

    dest_root_dir = scpm.starsolo_output_dir(project_name)

    def _sample_name_fn(donor):
        safe_donor = re.sub(r"[^A-Za-z0-9_-]", "", str(donor))
        return sample_name_template.format(donor=safe_donor, condition=condition_value)

    split_result = split_pooled_sample_by_donor(
        filtered_dir, standardized_df, dest_root_dir, _sample_name_fn,
        singlets_only=singlets_only, singlet_flag_col=singlet_flag_col,
        singlet_flag_value=singlet_flag_value, min_cells_per_donor=min_cells_per_donor,
    )

    metadata_rows = [
        {
            "sample": entry["sample_name"],
            "donor": entry["donor"],
            "condition": condition_value,
            "source_pooled_sample": pooled_sample_name,
        }
        for entry in split_result["per_donor_summary"]
    ]
    metadata_df = append_sample_metadata_rows(
        scpm.metadata_path(project_name), metadata_rows, sample_col="sample",
    )
    split_result["metadata_df"] = metadata_df
    return split_result


# ---------------------------------------------------------------------------
# Kang et al. 2018 (GSE96583) -- dataset-specific convenience notes
# ---------------------------------------------------------------------------
#
# NOT a separate code path -- this dataset uses every function above
# exactly as documented, with these specific column-name arguments:
#
#   external_table_df = load_external_barcode_table(
#       "GSE96583_batch2.total.tsne.df.tsv.gz")  # or the muscData-derived equivalent
#   split_pooled_project_sample(
#       project_name, "batch2_control", external_table_df,
#       donor_col="ind", condition_value="ctrl", condition_col="stim",
#       extra_cols=["cell", "cluster"],
#       singlets_only=True, singlet_flag_col="multiplets", singlet_flag_value="singlet",
#       sample_name_template="{donor}_{condition}",
#   )
#   # ... and again for the "batch2_stim" pooled sample with condition_value="stim"
#
# The exact column names/capitalization above (particularly whether the
# raw GEO file literally uses "ind"/"stim"/"cell"/"cluster"/"multiplets"
# verbatim, vs. some other casing/naming in the RAW GEO file specifically
# vs. the muscData-curated copy) should be CONFIRMED by inspecting the
# actual downloaded file's own header row before running this for real --
# see load_external_barcode_table()'s own docstring for why this module
# does not hard-assume any of Kang's own column names internally.
