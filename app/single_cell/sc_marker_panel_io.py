"""
single_cell/sc_marker_panel_io.py

Standalone file-parsing helper for Step 8b's "upload a marker panel
file" feature (2026-08-25, per direct user request) -- lets a user
provide an entire marker gene panel table at once (.csv, .txt, or
.xlsx), instead of typing each cell type's genes one at a time via the
manual "Cell type name" / "Marker genes" text inputs.

--- Supported layouts ---
Both common layouts are auto-detected, so a user doesn't need to know
in advance which one their file uses:

  1. "Wide" layout -- one row per cell type, with multiple genes
     packed into a single cell, separated by comma, semicolon, or pipe:
         cell_type,genes
         CD4 T cell,"IL7R, CD3D, CD3E, CD3G, CD4"
         B cell,"MS4A1, CD79A, CD79B, CD19"

  2. "Long" layout -- one gene per row, with the SAME cell_type value
     repeated across multiple rows (a common export shape from a
     spreadsheet or database query):
         cell_type,gene
         B cell,MS4A1
         B cell,CD79A
         B cell,CD79B
         B cell,CD19

Both are handled by the SAME parsing logic: each row's gene-column
value is split on comma/semicolon/pipe (a genuinely single-gene value
simply has nothing to split on), and results are grouped/de-duplicated
by cell_type while preserving first-seen gene order.

--- Column name detection ---
Recognizes common header name variants for both the cell-type column
("cell_type", "celltype", "cell type", "type", "label", "name") and the
gene column ("genes", "gene", "markers", "marker_genes", "marker
genes", "marker gene") -- case-insensitively. If neither is recognized
but the file has EXACTLY 2 columns, falls back to a positional
interpretation (first column = cell type, second column = genes)
rather than failing outright, since a simple two-column file without a
"proper" header is a common, reasonable case to still support.

--- File format support ---
.csv and .txt are both read via pandas with sep=None, engine="python"
(auto-detects comma, tab, semicolon, etc.) -- .xlsx via pandas'
read_excel (requires the openpyxl package, already used elsewhere in
this project's spreadsheet-handling code).
"""
import io
import re

import pandas as pd

_CELL_TYPE_COLUMN_ALIASES = {"cell_type", "celltype", "cell type", "type", "label", "name"}
_GENE_COLUMN_ALIASES = {"genes", "gene", "markers", "marker_genes", "marker genes", "marker gene"}
_GENE_SPLIT_PATTERN = re.compile(r"[;,|]+")


def _find_column(df_columns, aliases):
    "Find the first column in df_columns whose (lowercased, stripped) name matches one of `aliases`, or None if none match."
    normalized = {c: c.strip().lower() for c in df_columns}
    for original, norm in normalized.items():
        if norm in aliases:
            return original
    return None


def parse_marker_panel_dataframe(df):
    """
    Parse an already-loaded DataFrame into a {cell_type: [gene, gene,
    ...]} dict -- see this module's own docstring for the full
    supported-layout and column-detection rules.

    Raises ValueError with a clear, actionable message if the required
    columns can't be identified, or if no valid rows are found at all.
    """
    if df.shape[1] < 2:
        raise ValueError(
            "Expected at least 2 columns (a cell type column and a gene column) -- "
            f"found only {df.shape[1]}."
        )

    cell_type_col = _find_column(df.columns, _CELL_TYPE_COLUMN_ALIASES)
    gene_col = _find_column(df.columns, _GENE_COLUMN_ALIASES)

    if cell_type_col is None or gene_col is None:
        # Fall back to positional (first two columns) -- a common case
        # for a simple two-column file without a recognized header name.
        if df.shape[1] == 2:
            cell_type_col, gene_col = df.columns[0], df.columns[1]
        else:
            raise ValueError(
                "Could not identify a cell-type column and a gene column. Expected "
                "column headers like 'cell_type' and 'genes' (case-insensitive), or "
                "exactly 2 columns with cell type in the first and genes in the second."
            )

    panels = {}
    for _, row in df.iterrows():
        cell_type = str(row[cell_type_col]).strip()
        gene_value = row[gene_col]
        if not cell_type or cell_type.lower() == "nan":
            continue
        if pd.isna(gene_value):
            continue
        gene_value = str(gene_value).strip()
        if not gene_value:
            continue

        # Handles BOTH wide format (multiple genes in one cell, split
        # by comma/semicolon/pipe) AND long format (a single gene per
        # row, multiple rows sharing the same cell_type) uniformly -- a
        # single-gene value simply has nothing to split on.
        genes_in_this_row = [g.strip() for g in _GENE_SPLIT_PATTERN.split(gene_value) if g.strip()]

        existing = panels.setdefault(cell_type, [])
        for gene in genes_in_this_row:
            if gene not in existing:
                existing.append(gene)

    if not panels:
        raise ValueError("No valid (cell type, gene) rows were found in this file.")

    return panels


def load_marker_panel_file(filename, file_bytes):
    """
    Load a marker-panel file (.csv, .txt, or .xlsx) from raw bytes
    (e.g. Streamlit's own UploadedFile.getvalue()) and parse it into a
    {cell_type: [gene, ...]} dict.

    filename: used ONLY to determine the file extension/format -- not
        read from disk (file_bytes is the actual content).
    file_bytes: the raw file content as bytes.

    Raises ValueError (with a clear, user-facing message) for any
    parsing failure, including an unsupported extension or a missing
    Excel-reading engine.
    """
    ext = filename.lower().rsplit(".", 1)[-1] if "." in filename else ""

    if ext in ("csv", "txt"):
        try:
            df = pd.read_csv(io.BytesIO(file_bytes), sep=None, engine="python")
        except Exception as e:
            raise ValueError(f"Could not parse this file as CSV/TSV: {e}")
    elif ext in ("xlsx", "xls"):
        try:
            df = pd.read_excel(io.BytesIO(file_bytes))
        except ImportError:
            raise ValueError(
                "Reading .xlsx files requires the `openpyxl` package, which is not "
                "installed in this environment."
            )
        except Exception as e:
            raise ValueError(f"Could not parse this Excel file: {e}")
    else:
        raise ValueError(f"Unsupported file type '.{ext}' -- please upload a .csv, .txt, or .xlsx file.")

    return parse_marker_panel_dataframe(df)
