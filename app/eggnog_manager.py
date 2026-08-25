"""
eggnog_manager.py

Shared infrastructure for eggNOG-mapper -- orthology-based functional
annotation for ANY organism, used as a bridge for two distinct
downstream needs across this app's pipelines:

  1. Single-cell Phase 3.7 cell-type annotation: for a NON-MODEL
     organism (e.g. Austrofundulus limnaeus/killifish), CellTypist and
     SingleR/celldex are both fundamentally human/mouse-centric (see
     the design discussion this module grew out of) and provide no
     real annotation path. eggNOG-mapper's own output includes a
     "Preferred_name" column -- a human-readable gene symbol
     transferred via orthology (often literally matching well-known
     human/mouse gene nomenclature, e.g. "cd3d") -- letting a familiar
     marker panel (CD3D, MS4A1, ...) resolve directly against a
     non-model organism's own gene IDs, without needing a species-
     specific reference atlas at all.

  2. Bulk RNA-Seq Ontology Analysis: any organism lacking a curated
     Bioconductor org.*.eg.db package (i.e. anything outside this
     project's REFERENCE_CATALOG preset species) currently has NO path
     to GO/KEGG/Reactome-style enrichment at all. eggNOG-mapper's own
     annotation output includes GO term and KEGG orthology (KO)
     assignments per gene, which can be fed into clusterProfiler's
     generic enricher() function (via a TERM2GENE table) as a direct
     substitute for enrichGO()/enrichKEGG(), which both require a
     curated org.db/KEGG organism code that simply doesn't exist for
     most non-model organisms.

This is deliberately ONE shared piece of infrastructure consumed
differently by each pipeline, rather than two separate builds -- both
downstream needs are fed by the exact same one-time, per-reference
eggNOG-mapper run and its resulting .emapper.annotations file.

--- Design decisions (2026-08-20) ---

Local database install, NOT the public eggNOG-mapper web service --
matches this project's existing 2026-08-21 decision (recorded
separately) and the same reasoning already applied elsewhere in this
app: the public web service has no documented automated submission API
and reported reliability issues, so a local database is used instead.

The eggNOG core database is genuinely large (~45GB core annotation/taxonomy
database + ~4GB DIAMOND search database = ~49GB total, confirmed
directly from eggNOG-mapper's own documentation) -- treated as a
SHARED, project-independent, one-time resource using the EXACT SAME
ensure_shared_resource() locking/atomic-build infrastructure already
used by reference_manager.py for reference genomes and indices, NOT
duplicated per project/user/organism. This module intentionally does
NOT implement its own separate locking scheme -- it imports and reuses
reference_manager.ensure_shared_resource() directly.

Database setup is an explicit, deliberate, ADMIN-style action (a large,
slow, disk-heavy one-time operation) -- NOT triggered implicitly by any
per-project analysis step. Since this project has no user role/profile
system yet (a known, explicitly deferred piece of future work), "admin-
only" is enforced here only as a UI-level convention (a clearly-labeled,
separately-gated setup action with an explicit disk-space warning and
confirmation step) rather than a real technical permission check --
this module's own functions do not themselves enforce any authorization,
they just make the action's cost/scope impossible to trigger by
accident.

check_disk_space() is called BEFORE attempting the download and
returns a clear, actionable warning (not a vague failure partway
through a 45GB download) specifically flagging this as most likely to
matter for a LOCAL/laptop install -- an HPC deployment with terabytes
of scratch space (this project's own primary testing environment) is
very unlikely to be constrained by this, but a local install easily
could be.

--- Reference GTF -> protein FASTA -> eggNOG-mapper pipeline ---

eggNOG-mapper requires PROTEIN (or CDS) sequences as input, not a
genome + GTF directly. extract_protein_fasta_for_eggnog() uses gffread
(already an existing dependency in this project's environment.yml,
used elsewhere for GFF3->GTF conversion and transcript extraction) via
its "-y" flag, confirmed directly from gffread's own documentation to
produce protein-translated FASTA output from a genome + GTF/GFF.

run_eggnog_mapper() is a one-time, PER-REFERENCE operation (like the
genome download/index build it sits alongside) -- its result (the
parsed .emapper.annotations file) is what actually gets reused by both
downstream consumers, not the raw protein FASTA or intermediate search
files.
"""
import csv
import json
import os
import re
import shutil
import subprocess


# ---------------------------------------------------------------------------
# Availability + disk space checks
# ---------------------------------------------------------------------------

# Confirmed directly from eggNOG-mapper's own documentation: ~45GB for
# the core eggnog.db + eggnog.taxa.db, plus ~4GB for the DIAMOND search
# database (the default/recommended search mode) -- ~49GB total for a
# standard install. A generous safety margin is added on top (rather
# than warning right at the exact boundary) since decompression and
# intermediate download artifacts transiently use additional space
# during the download itself, not just the final installed size.
EGGNOG_CORE_DB_SIZE_GB = 45
EGGNOG_DIAMOND_DB_SIZE_GB = 4
EGGNOG_TOTAL_REQUIRED_GB = EGGNOG_CORE_DB_SIZE_GB + EGGNOG_DIAMOND_DB_SIZE_GB
EGGNOG_RECOMMENDED_FREE_GB = 65  # total + safety margin for in-progress downloads


def eggnog_mapper_available():
    "Check whether emapper.py (eggNOG-mapper's main CLI) is available on PATH."
    return shutil.which("emapper.py") is not None


def download_eggnog_data_script_available():
    "Check whether download_eggnog_data.py (eggNOG-mapper's own database-download script) is available on PATH."
    return shutil.which("download_eggnog_data.py") is not None


def check_disk_space(target_dir):
    """
    Check available disk space at (or near) target_dir BEFORE
    attempting the eggNOG database download -- surfaced as an explicit,
    actionable warning rather than letting a ~49GB download fail
    partway through with a generic "no space left on device" error.

    target_dir: the directory the database will be installed into (or
        its intended parent, if it doesn't exist yet -- this function
        walks up to the nearest existing ancestor directory to check,
        so a not-yet-created target directory doesn't cause an error
        here).

    Returns a dict:
        {
            "free_gb": float,           # currently available, in GB
            "required_gb": int,         # EGGNOG_TOTAL_REQUIRED_GB
            "recommended_free_gb": int, # EGGNOG_RECOMMENDED_FREE_GB
            "sufficient": bool,         # free_gb >= required_gb
            "comfortable": bool,        # free_gb >= recommended_free_gb
            "message": str,             # plain-language summary, tiered
                                         # by severity (see below)
        }
    """
    check_path = target_dir
    while check_path and not os.path.isdir(check_path):
        parent = os.path.dirname(check_path.rstrip(os.sep))
        if parent == check_path:
            check_path = "/"
            break
        check_path = parent

    total, used, free = shutil.disk_usage(check_path)
    free_gb = free / (1024 ** 3)

    sufficient = free_gb >= EGGNOG_TOTAL_REQUIRED_GB
    comfortable = free_gb >= EGGNOG_RECOMMENDED_FREE_GB

    if not sufficient:
        message = (
            f"🔴 Only {free_gb:.1f} GB free at {check_path} -- the eggNOG database "
            f"needs approximately {EGGNOG_TOTAL_REQUIRED_GB} GB ({EGGNOG_CORE_DB_SIZE_GB} GB core "
            f"annotation database + {EGGNOG_DIAMOND_DB_SIZE_GB} GB DIAMOND search database). "
            "This download will very likely fail partway through or fill this disk -- "
            "this is especially worth double-checking if you're running this on a local "
            "machine/laptop rather than an HPC system with dedicated scratch space."
        )
    elif not comfortable:
        message = (
            f"🟡 {free_gb:.1f} GB free at {check_path} -- technically enough for the "
            f"~{EGGNOG_TOTAL_REQUIRED_GB} GB database itself, but with little margin for "
            "temporary download/decompression overhead. Recommended: at least "
            f"{EGGNOG_RECOMMENDED_FREE_GB} GB free before proceeding, particularly on a "
            "local machine/laptop where this space can't easily be expanded."
        )
    else:
        message = (
            f"🟢 {free_gb:.1f} GB free at {check_path} -- comfortably sufficient for the "
            f"~{EGGNOG_TOTAL_REQUIRED_GB} GB eggNOG database."
        )

    return {
        "free_gb": round(free_gb, 1),
        "required_gb": EGGNOG_TOTAL_REQUIRED_GB,
        "recommended_free_gb": EGGNOG_RECOMMENDED_FREE_GB,
        "sufficient": sufficient,
        "comfortable": comfortable,
        "message": message,
    }


def eggnog_database_is_installed(data_dir):
    """
    Check whether the core eggNOG database files already exist at
    data_dir -- specifically eggnog.db (the large core annotation
    database) and eggnog_proteins.dmnd (the DIAMOND search database,
    the default/recommended search mode). Does NOT check optional
    databases (PFAM, MMseqs2, HMMER-per-taxon) -- those are separate,
    optional add-ons this module does not manage.

    Returns True only if BOTH required files are present -- a partial
    download (e.g. interrupted mid-way) correctly reports as NOT
    installed, rather than a misleading partial-success reading.
    """
    core_db = os.path.join(data_dir, "eggnog.db")
    diamond_db = os.path.join(data_dir, "eggnog_proteins.dmnd")
    return os.path.isfile(core_db) and os.path.isfile(diamond_db)


# ---------------------------------------------------------------------------
# Database setup (shared, admin-style action)
# ---------------------------------------------------------------------------

def download_eggnog_database(data_dir, ensure_shared_resource_fn, wait_message_callback=None,
                              timeout=None):
    """
    Download the core eggNOG database (eggnog.db + eggnog.taxa.db +
    the DIAMOND search database) to data_dir, using the SAME shared-
    resource locking/atomic-build infrastructure already used by
    reference_manager.py for genome references and indices -- see this
    module's own docstring for the full rationale on why this is
    treated as one shared, project-independent resource rather than
    duplicated per project/organism.

    ensure_shared_resource_fn: pass reference_manager.ensure_shared_resource
        directly (dependency-injected rather than imported at module
        level here, to keep this module's only real coupling to
        reference_manager.py explicit and visible at the call site,
        rather than a hidden top-of-file import -- this mirrors how
        sc_cellqc_manager.py's own diagnose_starsolo_matrix_for_mito()
        takes mito_gene_ids as an explicit parameter rather than
        reaching for global state).

    data_dir: the FINAL shared destination directory for the database
        (e.g. a dedicated data/shared_resources/eggnog/ path, analogous
        to reference_manager.py's shared_star_index_dir() convention).

    IMPORTANT -- caller responsibility: check_disk_space(data_dir)
    should ALWAYS be called and its result shown to the user BEFORE
    calling this function -- this function does NOT itself block a
    download on insufficient disk space (consistent with this
    project's general "check-and-warn, don't silently prevent" pattern
    elsewhere), but the caller (UI layer) is expected to require an
    explicit user confirmation if check_disk_space() reports "not
    sufficient".

    Returns (success: bool, message: str, built_by_this_call: bool) --
    same three-tuple shape as reference_manager.ensure_shared_resource()
    itself, since this function is a thin wrapper around it.
    """
    if not download_eggnog_data_script_available():
        return False, (
            "download_eggnog_data.py was not found on this system -- the `eggnog-mapper` "
            "package needs to be installed in your environment before the database can be "
            "downloaded."
        ), False

    def _build_fn(temp_dir):
        cmd = ["download_eggnog_data.py", "-y", "--data_dir", temp_dir]
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True,
                timeout=timeout or 21600,  # 6 hours -- a genuinely large, slow download
            )
            return True, result.stdout + result.stderr
        except subprocess.CalledProcessError as e:
            return False, f"eggNOG database download failed: {(e.stdout or '') + (e.stderr or '')}"
        except subprocess.TimeoutExpired:
            return False, "eggNOG database download timed out."

    return ensure_shared_resource_fn(data_dir, build_fn=_build_fn, wait_message_callback=wait_message_callback)


# ---------------------------------------------------------------------------
# Protein extraction (genome + GTF -> protein FASTA, via gffread)
# ---------------------------------------------------------------------------

def extract_protein_fasta_for_eggnog(genome_fasta_path, gtf_path, dest_fasta_path):
    """
    Extract protein sequences (translated CDS) from a genome FASTA +
    GTF/GFF annotation using gffread's "-y" flag -- confirmed directly
    from gffread's own documentation to produce protein-translated
    FASTA output, one sequence per transcript. This is the required
    input FORMAT for eggNOG-mapper (protein or CDS sequences,
    "--itype proteins" or "--itype CDS") -- eggNOG-mapper cannot be run
    directly against a genome + GTF the way STAR/Salmon can.

    Already-installed dependency: gffread is already listed in this
    project's environment.yml (used elsewhere for GFF3->GTF conversion
    and transcript extraction) -- no new external tool dependency is
    introduced by this function.

    Returns (success: bool, message: str).
    """
    if shutil.which("gffread") is None:
        return False, (
            "gffread is not installed on this system. It's required to extract protein "
            "sequences for eggNOG-mapper. (It's included in this project's environment.yml.)"
        )

    cmd = ["gffread", gtf_path, "-g", genome_fasta_path, "-y", dest_fasta_path]
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=3600)
        return True, "Protein sequences extracted successfully for eggNOG-mapper."
    except subprocess.CalledProcessError as e:
        return False, f"gffread protein extraction failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, "gffread protein extraction timed out after 1 hour."


# ---------------------------------------------------------------------------
# Running eggNOG-mapper itself
# ---------------------------------------------------------------------------

def run_eggnog_mapper(protein_fasta_path, output_dir, output_prefix, data_dir,
                       cpu=4, timeout=None):
    """
    Run eggNOG-mapper (emapper.py) on a protein FASTA -- a one-time,
    PER-REFERENCE operation (like a genome index build), NOT re-run
    per-project or per-sample. The resulting .emapper.annotations file
    is what both downstream consumers (single-cell marker resolution,
    bulk GO/KEGG enrichment) actually read -- see
    parse_eggnog_annotations() below.

    protein_fasta_path: output of extract_protein_fasta_for_eggnog().
    output_dir: directory to write emapper's output files into.
    output_prefix: filename prefix (emapper.py's own "-o" argument) --
        the resulting file will be
        {output_dir}/{output_prefix}.emapper.annotations.
    data_dir: the eggNOG database directory (must already be installed
        -- see download_eggnog_database()/eggnog_database_is_installed()
        above).
    cpu: number of CPU threads to use (emapper.py's own "--cpu" flag).

    Returns (success: bool, annotations_path_or_None: str, message: str).
    """
    if not eggnog_mapper_available():
        return False, None, (
            "emapper.py was not found on this system. The `eggnog-mapper` package needs to "
            "be installed in your environment before this step can run."
        )
    if not eggnog_database_is_installed(data_dir):
        return False, None, (
            f"The eggNOG database has not been installed yet at {data_dir} -- run the "
            "one-time database setup first (see download_eggnog_database())."
        )
    if not os.path.isfile(protein_fasta_path):
        return False, None, f"Protein FASTA not found: {protein_fasta_path}"

    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        "emapper.py",
        "-i", protein_fasta_path,
        "--itype", "proteins",
        "-o", output_prefix,
        "--output_dir", output_dir,
        "--data_dir", data_dir,
        "--cpu", str(cpu),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True,
            timeout=timeout or 21600,  # 6 hours -- a real, potentially large genome-wide run
        )
        log = result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, None, f"eggNOG-mapper failed: {(e.stdout or '') + (e.stderr or '')}"
    except subprocess.TimeoutExpired:
        return False, None, "eggNOG-mapper timed out."

    annotations_path = os.path.join(output_dir, f"{output_prefix}.emapper.annotations")
    if not os.path.isfile(annotations_path):
        return False, None, (
            "eggNOG-mapper completed without an error, but the expected output file "
            f"({annotations_path}) was not found -- this may indicate a version mismatch "
            "or an unexpected output naming convention."
        )

    return True, annotations_path, f"eggNOG-mapper annotation complete.\n{log}"


# ---------------------------------------------------------------------------
# Parsing eggNOG-mapper's output
# ---------------------------------------------------------------------------

# Confirmed directly from eggNOG-mapper's own documented output column
# schema (the .emapper.annotations TSV). File structure: comment lines
# prefixed with "##" (metadata), a header line prefixed with "#", data
# rows, and (unless --no_file_comments is used) footer comment lines --
# all of which must be correctly skipped/parsed around.
_EGGNOG_ANNOTATION_COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


def parse_eggnog_annotations(annotations_path):
    """
    Parse an eggNOG-mapper .emapper.annotations TSV file into a plain
    DataFrame -- handles the "##"-prefixed metadata/footer comment
    lines and the "#"-prefixed header line correctly (confirmed via
    eggNOG-mapper's own documented output file structure), rather than
    naively passing the file straight to pandas.read_csv (which would
    otherwise choke on or misparse the comment lines).

    Missing-value placeholders ("-", eggNOG-mapper's own convention for
    "no annotation available" in a given column) are normalized to
    real empty strings/NaN, so downstream code can use simple
    truthiness/notna() checks rather than needing to know eggNOG's own
    specific sentinel value.

    Returns a DataFrame with (at least) the columns in
    _EGGNOG_ANNOTATION_COLUMNS. Returns None if the file doesn't exist.
    """
    import pandas as pd

    if not os.path.isfile(annotations_path):
        return None

    header_line_idx = None
    with open(annotations_path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("#query") or (line.startswith("#") and "query" in line.split("\t")[0]):
                header_line_idx = i
                break
            if line.startswith("##"):
                continue
            if not line.startswith("#"):
                break

    if header_line_idx is None:
        # Fallback: no explicit header row with a "#query" prefix was
        # found (older/newer eggNOG-mapper versions have varied this
        # slightly) -- read all non-"##"-prefixed lines and apply this
        # module's own known column list directly, on the assumption
        # that data rows still appear in the standard column order.
        rows = []
        with open(annotations_path, "r", errors="replace") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                rows.append(line.rstrip("\n").split("\t"))
        df = pd.DataFrame(rows)
        n_cols = min(len(_EGGNOG_ANNOTATION_COLUMNS), df.shape[1])
        df = df.iloc[:, :n_cols]
        df.columns = _EGGNOG_ANNOTATION_COLUMNS[:n_cols]
    else:
        df = pd.read_csv(
            annotations_path, sep="\t", skiprows=header_line_idx,
            comment=None, header=0,
        )
        df.columns = [c.lstrip("#") for c in df.columns]
        # Drop any trailing "##"-prefixed footer summary lines that
        # pandas may have read in as literal data rows (they appear
        # AFTER the real data, so they show up as extra rows whose
        # first column starts with "##").
        first_col = df.columns[0]
        df = df[~df[first_col].astype(str).str.startswith("##")].reset_index(drop=True)

    # Normalize eggNOG's own "-" missing-value sentinel to real NaN
    # across every column, so downstream code can use ordinary
    # pandas truthiness/notna() checks.
    df = df.replace("-", pd.NA)

    return df


# ---------------------------------------------------------------------------
# Consumer 1: gene_symbol_map for single-cell marker resolution
# ---------------------------------------------------------------------------

def build_gene_symbol_map_from_eggnog(annotations_df, tx2gene_map=None):
    """
    Build a gene_id -> gene_symbol mapping from eggNOG-mapper's own
    "Preferred_name" column -- the KEY BRIDGE this module provides for
    non-model-organism single-cell marker resolution (see this
    module's own docstring, point 1).

    annotations_df: output of parse_eggnog_annotations().
    tx2gene_map: optional dict {transcript_id: gene_id} (e.g. from
        reference_manager.extract_tx2gene_from_gtf()) -- eggNOG-mapper's
        "query" column identifies each INPUT sequence, which (since
        extract_protein_fasta_for_eggnog() extracts one protein per
        TRANSCRIPT via gffread) is typically a transcript ID, not a
        gene ID. If tx2gene_map is provided, query IDs are mapped back
        to their gene ID via this table; if a given query ID isn't
        found in tx2gene_map (or if tx2gene_map is None entirely), the
        query ID itself is used as the key directly -- this keeps the
        function usable even without a tx2gene mapping on hand, at the
        cost of keying by transcript ID instead of gene ID in that
        case (the caller should be aware this could then produce
        multiple rows per gene, one per transcript, unless collapsed
        beforehand).

    Only rows with a REAL (non-missing) Preferred_name are included --
    a query with no orthology-derived name at all (a common, expected
    outcome for many genes, not an error) is simply not present in the
    output rather than being included with an empty/placeholder value.
    When multiple transcripts for the same gene resolve to DIFFERENT
    Preferred_name values (can happen with alternate/ambiguous
    orthology calls), the FIRST one encountered is kept -- callers
    wanting a specific tie-breaking rule beyond "first" should
    pre-filter annotations_df themselves before calling this function.

    Returns a dict {gene_id: gene_symbol} -- directly compatible with
    reference_manager.save_gene_symbol_map_csv()'s own expected input
    shape, so a caller can do:
        mapping = build_gene_symbol_map_from_eggnog(df, tx2gene)
        rm.save_gene_symbol_map_csv(mapping, pm.gene_symbol_map_path(project))
    to persist this exactly like any other gene symbol map already
    used elsewhere in this app.
    """
    if annotations_df is None or "Preferred_name" not in annotations_df.columns:
        return {}

    mapping = {}
    for _, row in annotations_df.iterrows():
        preferred_name = row.get("Preferred_name")
        if preferred_name is None or (hasattr(preferred_name, "__len__") and len(str(preferred_name)) == 0):
            continue
        try:
            import pandas as pd
            if pd.isna(preferred_name):
                continue
        except Exception:
            pass

        query_id = row.get("query")
        if query_id is None:
            continue

        gene_id = query_id
        if tx2gene_map and query_id in tx2gene_map:
            gene_id = tx2gene_map[query_id]

        if gene_id not in mapping:
            mapping[gene_id] = str(preferred_name)

    return mapping


# ---------------------------------------------------------------------------
# Consumer 2: GO / KEGG TERM2GENE tables for clusterProfiler::enricher()
# ---------------------------------------------------------------------------

def build_go_term2gene_from_eggnog(annotations_df, tx2gene_map=None):
    """
    Build a long-format (TERM2GENE-shaped) DataFrame of GO term ->
    gene_id assignments from eggNOG-mapper's "GOs" column -- for
    clusterProfiler::enricher()'s generic, org.db-independent
    enrichment path (see this module's own docstring, point 2). This
    is the direct substitute for clusterProfiler::enrichGO(), which
    REQUIRES a curated Bioconductor org.*.eg.db package that simply
    doesn't exist for most non-model organisms.

    The "GOs" column is a comma-separated list of GO term IDs per gene
    (eggNOG-mapper's own format) -- this function explodes that into
    ONE ROW PER (GO_term, gene) PAIR, the long format
    clusterProfiler::enricher()'s own TERM2GENE argument expects.

    annotations_df: output of parse_eggnog_annotations().
    tx2gene_map: see build_gene_symbol_map_from_eggnog()'s own
        docstring -- same optional transcript-to-gene collapsing logic
        applies here.

    Returns a DataFrame with columns "term" (GO ID, e.g. "GO:0005515")
    and "gene" (gene_id) -- ready to be written to a CSV/passed
    directly as clusterProfiler::enricher()'s TERM2GENE argument.
    Returns an empty DataFrame (same columns, zero rows) if no GO
    annotations are present at all.
    """
    import pandas as pd

    if annotations_df is None or "GOs" not in annotations_df.columns:
        return pd.DataFrame(columns=["term", "gene"])

    rows = []
    for _, row in annotations_df.iterrows():
        go_terms = row.get("GOs")
        if go_terms is None:
            continue
        try:
            if pd.isna(go_terms):
                continue
        except Exception:
            pass

        query_id = row.get("query")
        if query_id is None:
            continue
        gene_id = query_id
        if tx2gene_map and query_id in tx2gene_map:
            gene_id = tx2gene_map[query_id]

        for term in str(go_terms).split(","):
            term = term.strip()
            if term:
                rows.append({"term": term, "gene": gene_id})

    return pd.DataFrame(rows, columns=["term", "gene"])


def build_kegg_term2gene_from_eggnog(annotations_df, tx2gene_map=None):
    """
    Build a long-format (TERM2GENE-shaped) DataFrame of KEGG Orthology
    (KO) term -> gene_id assignments from eggNOG-mapper's "KEGG_ko"
    column -- the direct substitute for
    clusterProfiler::enrichKEGG()'s own organism-code-based lookup,
    which (like enrichGO()) has no path for most non-model organisms.

    Same shape/logic as build_go_term2gene_from_eggnog() above, applied
    to the "KEGG_ko" column instead of "GOs" -- eggNOG-mapper's KEGG_ko
    values look like "ko:K04525" (with a "ko:" prefix); this prefix is
    stripped so the resulting term column contains bare KO identifiers
    (e.g. "K04525"), matching the convention clusterProfiler's own KEGG-
    related functions/reference tables typically expect.

    Returns a DataFrame with columns "term" (bare KO ID) and "gene"
    (gene_id). Returns an empty DataFrame (same columns, zero rows) if
    no KEGG_ko annotations are present at all.
    """
    import pandas as pd

    if annotations_df is None or "KEGG_ko" not in annotations_df.columns:
        return pd.DataFrame(columns=["term", "gene"])

    rows = []
    for _, row in annotations_df.iterrows():
        ko_terms = row.get("KEGG_ko")
        if ko_terms is None:
            continue
        try:
            if pd.isna(ko_terms):
                continue
        except Exception:
            pass

        query_id = row.get("query")
        if query_id is None:
            continue
        gene_id = query_id
        if tx2gene_map and query_id in tx2gene_map:
            gene_id = tx2gene_map[query_id]

        for term in str(ko_terms).split(","):
            term = term.strip()
            if term.startswith("ko:"):
                term = term[len("ko:"):]
            if term:
                rows.append({"term": term, "gene": gene_id})

    return pd.DataFrame(rows, columns=["term", "gene"])


def save_term2gene_csv(term2gene_df, dest_path):
    """
    Write a TERM2GENE-shaped DataFrame (columns "term", "gene") to a
    CSV -- mirrors reference_manager.save_gene_symbol_map_csv()'s own
    simple "write this mapping to a persistent file" convention, so
    this can be saved alongside a reference's other derived files and
    reloaded on a later Ontology Analysis run without needing to
    re-parse the raw .emapper.annotations file every time.
    """
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    term2gene_df.to_csv(dest_path, index=False)
    return dest_path
# ---------------------------------------------------------------------------
# Source-organism attribution (from seed_ortholog's taxid prefix)
# ---------------------------------------------------------------------------

# A modest, commonly-relevant set of NCBI taxonomy IDs -- covers the
# organisms most likely to appear as a "source" in eggNOG's own broad
# orthologous-group database, and every organism in this app's OWN
# REFERENCE_CATALOG/org.*.eg.db preset list (so a preset-species source
# is always resolved to a real, familiar name, not a bare numeral).
# This is NOT an exhaustive list of eggNOG's full ~12,500+ reference
# taxa -- an unrecognized taxid is displayed AS ITS RAW NUMBER rather
# than guessed at (see resolve_taxid_to_organism_name() below), so a
# missing entry here degrades gracefully rather than showing something
# actively wrong.
COMMON_TAXID_TO_ORGANISM = {
    "9606": "Homo sapiens (human)",
    "10090": "Mus musculus (mouse)",
    "10116": "Rattus norvegicus (rat)",
    "7955": "Danio rerio (zebrafish)",
    "7227": "Drosophila melanogaster (fly)",
    "6239": "Caenorhabditis elegans (roundworm)",
    "4932": "Saccharomyces cerevisiae (yeast)",
    "559292": "Saccharomyces cerevisiae S288C (yeast)",
    "511145": "Escherichia coli str. K-12 (E. coli)",
    "83333": "Escherichia coli K-12 (E. coli)",
    "3702": "Arabidopsis thaliana (thale cress)",
    "9031": "Gallus gallus (chicken)",
    "9598": "Pan troglodytes (chimpanzee)",
    "9544": "Macaca mulatta (rhesus macaque)",
    "8355": "Xenopus laevis (African clawed frog)",
    "8364": "Xenopus tropicalis (western clawed frog)",
    "9615": "Canis lupus familiaris (dog)",
    "9823": "Sus scrofa (pig)",
    "9913": "Bos taurus (cow)",
    "7719": "Ciona intestinalis (sea squirt)",
    "6669": "Daphnia pulex (water flea)",
    "6945": "Ixodes scapularis (tick)",
    "7460": "Apis mellifera (honeybee)",
    "7159": "Aedes aegypti (mosquito)",
}


def parse_seed_ortholog_taxid(seed_ortholog_value):
    """
    Extract the NCBI taxonomy ID prefix from an eggNOG-mapper
    "seed_ortholog" value (format: "{taxid}.{protein_id}", e.g.
    "9606.ENSP00000000233" -> "9606") -- confirmed directly from
    eggNOG-mapper's own documented output format.

    Returns the taxid as a string, or None if seed_ortholog_value
    doesn't match the expected "{digits}.{anything}" pattern at all
    (e.g. missing/malformed value) -- callers should treat None as
    "source organism unknown", not an error.
    """
    if not seed_ortholog_value or not isinstance(seed_ortholog_value, str):
        return None
    match = re.match(r"^(\d+)\.", seed_ortholog_value)
    return match.group(1) if match else None


def resolve_taxid_to_organism_name(taxid):
    """
    Resolve a taxid string to a readable organism name via
    COMMON_TAXID_TO_ORGANISM -- falls back to a plainly-labeled "Taxon
    <id>" string (rather than guessing, or silently omitting it) for
    any taxid not in this app's own curated common list, so a real but
    less-common source organism is still visibly reported, just
    without a resolved common name.
    """
    if not taxid:
        return "Unknown"
    return COMMON_TAXID_TO_ORGANISM.get(str(taxid), f"Taxon {taxid} (not in this app's common-name list)")


def build_gene_symbol_map_with_source_from_eggnog(annotations_df, tx2gene_map=None):
    """
    Like build_gene_symbol_map_from_eggnog() (see that function's own
    docstring for the core logic this mirrors), but ALSO resolves and
    returns each gene's SOURCE ORGANISM (parsed from seed_ortholog's
    taxid prefix) -- answering "which organism did this gene's name
    actually come from?", which matters most under the default "auto"
    tax_scope, where the source can genuinely vary gene-by-gene.

    Returns a dict {gene_id: {"gene_symbol": str, "source_organism":
    str, "source_taxid": str_or_None}} -- a richer shape than
    build_gene_symbol_map_from_eggnog()'s plain {gene_id: symbol} dict,
    intended for a UI DISPLAY table (e.g. "gene X = CD3D, sourced from
    Homo sapiens") rather than direct use as a
    save_gene_symbol_map_csv()-compatible mapping -- callers wanting
    THAT simpler shape should still use
    build_gene_symbol_map_from_eggnog() instead (or derive it from this
    function's output by taking only the "gene_symbol" sub-value).
    """
    import pandas as pd

    if annotations_df is None or "Preferred_name" not in annotations_df.columns:
        return {}

    result = {}
    for _, row in annotations_df.iterrows():
        preferred_name = row.get("Preferred_name")
        try:
            if pd.isna(preferred_name):
                continue
        except Exception:
            if preferred_name is None:
                continue
        if not str(preferred_name).strip():
            continue

        query_id = row.get("query")
        if query_id is None:
            continue
        gene_id = query_id
        if tx2gene_map and query_id in tx2gene_map:
            gene_id = tx2gene_map[query_id]

        if gene_id in result:
            continue  # first-encountered wins, matching build_gene_symbol_map_from_eggnog's own tie-breaking

        taxid = parse_seed_ortholog_taxid(row.get("seed_ortholog"))
        result[gene_id] = {
            "gene_symbol": str(preferred_name),
            "source_organism": resolve_taxid_to_organism_name(taxid),
            "source_taxid": taxid,
        }

    return result


# ---------------------------------------------------------------------------
# --target_taxa scope options (Automatic / Closely-related / Human)
# ---------------------------------------------------------------------------

TARGET_TAXA_SCOPE_OPTIONS = {
    "auto": {
        "label": "Automatic (recommended default)",
        "target_taxa_arg": None,  # no --target_taxa restriction at all
        "explanation": (
            "eggNOG-mapper pulls annotation evidence from whichever curated ortholog is "
            "evolutionarily closest for EACH gene individually -- the most scientifically "
            "defensible choice per-gene, but the SOURCE organism behind any given gene's name "
            "will vary unpredictably from gene to gene (this app reports each gene's actual "
            "source organism so you can see this directly, rather than it being hidden). "
            "Generally gives the BEST overall coverage (fewest genes left completely "
            "unannotated), since it never refuses evidence from a real, valid ortholog just "
            "because it's not from your preferred organism."
        ),
    },
    "closely_related": {
        "label": "A specific, closely-related organism",
        "target_taxa_arg": "custom",  # user supplies their own taxid
        "explanation": (
            "Restricts annotation evidence to ONLY the organism you specify (e.g. zebrafish "
            "for a fish). More evolutionarily appropriate gene-name/ortholog transfer per "
            "gene, since divergence time is smaller -- but ⚠️ IMPORTANT: this does NOT, by "
            "itself, unlock standard GO/KEGG enrichment (enrichGO()/enrichKEGG()) unless that "
            "organism ALREADY has its own curated Bioconductor org.*.eg.db package installed "
            "in this app (check this app's own preset organism list first) -- otherwise, use "
            "the eggNOG-derived GO/KEGG TERM2GENE enrichment path instead, which works "
            "regardless. Also note: gene-name/functional-database coverage for most "
            "non-mammalian, non-human organisms is comparatively sparse, so restricting here "
            "will likely leave MORE genes completely unannotated than 'Automatic' or 'Human' "
            "would."
        ),
    },
    "human": {
        "label": "Human (Homo sapiens)",
        "target_taxa_arg": "9606",
        "explanation": (
            "Restricts annotation evidence to ONLY genes with a human ortholog. Real upside: "
            "human GO/KEGG/pathway curation is by far the deepest and most complete of any "
            "organism, so downstream Ontology Analysis will have the richest, most detailed "
            "term set to test against. Real downside: any gene whose orthologous group has "
            "NO human member at all will receive NO annotation whatsoever, rather than "
            "falling back to the next-closest available evidence -- this is a genuine "
            "coverage-vs-depth tradeoff (fewer genes annotated overall, but richer functional "
            "detail for the ones that are), not a strictly-better choice than Automatic."
        ),
    },
}
DEFAULT_TARGET_TAXA_SCOPE = "auto"


def build_emapper_target_taxa_args(scope_key, custom_taxid=None):
    """
    Build the actual emapper.py CLI argument list for the chosen
    TARGET_TAXA_SCOPE_OPTIONS entry -- kept as a single, explicit
    function (rather than inlined in run_eggnog_mapper()) so the exact
    command-line flag/value is fully visible and testable independent
    of the subprocess-running logic itself.

    scope_key: one of TARGET_TAXA_SCOPE_OPTIONS's own keys ("auto",
        "closely_related", "human").
    custom_taxid: REQUIRED when scope_key == "closely_related" -- a
        plain NCBI taxonomy ID string (e.g. "7955" for zebrafish) the
        user has supplied themselves; this app does not maintain its
        own "closely related organism" picker/catalog, since the
        correct choice is entirely dependent on the user's own
        specific non-model organism and cannot be meaningfully
        defaulted or guessed at.

    Returns a list of CLI argument strings to append to the emapper.py
    command (e.g. ["--target_taxa", "9606"]), or an empty list for
    scope_key == "auto" (no restriction at all).

    Raises ValueError if scope_key == "closely_related" but
    custom_taxid was not provided, or if scope_key isn't recognized at
    all -- both are real usage errors from the calling UI code, not
    something to silently default around.
    """
    if scope_key not in TARGET_TAXA_SCOPE_OPTIONS:
        raise ValueError(f"Unknown target_taxa scope: {scope_key!r}")

    if scope_key == "auto":
        return []
    if scope_key == "human":
        return ["--target_taxa", TARGET_TAXA_SCOPE_OPTIONS["human"]["target_taxa_arg"]]
    if scope_key == "closely_related":
        if not custom_taxid:
            raise ValueError(
                "scope_key='closely_related' requires custom_taxid to be provided "
                "(the user's own chosen organism's NCBI taxonomy ID)."
            )
        return ["--target_taxa", str(custom_taxid)]

    raise ValueError(f"Unknown target_taxa scope: {scope_key!r}")


# ---------------------------------------------------------------------------
# Annotation coverage summary
# ---------------------------------------------------------------------------

def build_annotation_coverage_summary(annotations_df, target_taxid=None):
    """
    Build a plain-language coverage report for one eggNOG-mapper run --
    "how many genes/transcripts did we actually get useful annotation
    for?" -- surfaced explicitly rather than leaving a user to
    discover coverage gaps only indirectly (e.g. by noticing their
    gene-symbol map or GO/KEGG term counts seem smaller than expected).

    annotations_df: output of parse_eggnog_annotations().
    target_taxid: if a specific --target_taxa restriction was used for
        this run (see build_emapper_target_taxa_args()), pass that
        same taxid here to additionally report how many of the
        annotated genes' evidence actually came from THAT organism
        specifically vs. genes that received SOME annotation but not
        from the requested organism (relevant mainly for a
        "closely_related" scope choice where the search wasn't
        strictly limited to one single taxon at the eggNOG-mapper
        level, or for auditing how much a "human" restriction actually
        cost in coverage vs. what 'auto' would have found). Pass None
        (the default) to skip this specific breakdown.

    Returns a dict:
        {
            "total_queries": int,
            "n_with_any_annotation": int,
            "n_with_gene_name": int,      # non-missing Preferred_name
            "n_with_go_terms": int,        # non-missing GOs
            "n_with_kegg_terms": int,      # non-missing KEGG_ko
            "n_from_target_taxid": int_or_None,  # only if target_taxid given
            "message": str,                # plain-language summary
        }
    Returns None if annotations_df is None or empty.
    """
    import pandas as pd

    if annotations_df is None or len(annotations_df) == 0:
        return None

    def _count_non_missing(col):
        if col not in annotations_df.columns:
            return 0
        return int(annotations_df[col].notna().sum())

    total = len(annotations_df)
    n_named = _count_non_missing("Preferred_name")
    n_go = _count_non_missing("GOs")
    n_kegg = _count_non_missing("KEGG_ko")

    # "Any annotation at all" -- a query counts if it has EITHER a name
    # OR GO terms OR KEGG terms (some queries may have one but not the
    # others, e.g. a KEGG assignment with no clean Preferred_name).
    any_mask = pd.Series(False, index=annotations_df.index)
    for col in ("Preferred_name", "GOs", "KEGG_ko"):
        if col in annotations_df.columns:
            any_mask = any_mask | annotations_df[col].notna()
    n_any = int(any_mask.sum())

    n_from_target = None
    if target_taxid and "seed_ortholog" in annotations_df.columns:
        n_from_target = int(
            annotations_df["seed_ortholog"]
            .apply(lambda v: parse_seed_ortholog_taxid(v) == str(target_taxid))
            .sum()
        )

    pct_any = round(100 * n_any / total, 1) if total else 0.0
    pct_named = round(100 * n_named / total, 1) if total else 0.0

    message_lines = [
        f"**{n_any:,} of {total:,} genes/transcripts ({pct_any}%)** received SOME eggNOG annotation.",
        f"- **{n_named:,} ({pct_named}%)** received a usable gene name (Preferred_name) -- these feed this reference's gene-symbol map.",
        f"- **{n_go:,}** received GO term annotation -- these feed GO enrichment.",
        f"- **{n_kegg:,}** received KEGG Orthology (KO) annotation -- these feed KEGG enrichment.",
    ]
    if n_from_target is not None:
        pct_target = round(100 * n_from_target / total, 1) if total else 0.0
        message_lines.append(
            f"- Of the annotated genes, **{n_from_target:,} ({pct_target}%)** specifically had "
            f"evidence sourced from the requested organism (taxid {target_taxid}) -- the "
            f"remainder either had no ortholog there at all, or (if a non-'human'/'auto' "
            f"scope was used) fell back to a different available organism's evidence."
        )

    return {
        "total_queries": total,
        "n_with_any_annotation": n_any,
        "n_with_gene_name": n_named,
        "n_with_go_terms": n_go,
        "n_with_kegg_terms": n_kegg,
        "n_from_target_taxid": n_from_target,
        "message": "\n".join(message_lines),
    }
