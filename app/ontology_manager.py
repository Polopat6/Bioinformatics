"""
ontology_manager.py

Backend for the Ontology Analysis workspace: runs Gene Ontology (GO),
KEGG, and Reactome enrichment via R's clusterProfiler (+ ReactomePA),
in three flavors:
  - ORA  (Over-Representation Analysis): enrichGO / enrichKEGG /
    enrichPathway -- takes a fixed list of "significant" genes plus a
    background "universe" and asks "are these genes enriched for any
    term more than chance?"
  - GSEA (Gene Set Enrichment Analysis): gseGO / gseKEGG / gsePathway --
    takes EVERY tested gene ranked by a score (log2FoldChange here) and
    asks "do this term's genes cluster toward one end of the ranking?"
  - compareCluster: runs ORA independently across multiple contrasts/
    gene lists at once and returns one combined table with a "Cluster"
    column, for side-by-side comparison.

Design decision -- ALL PLOTTING HAPPENS IN PYTHON (Plotly), not R: see
ontology_workspace.py's module docstring for the full rationale.
"""

import itertools
import json
import os
import re
import shutil
import subprocess

import pandas as pd


KEGG_ORGANISM_CODES = {
    "human": "hsa",
    "mouse": "mmu",
    "rat": "rno",
    "zebrafish": "dre",
    "fly": "dme",
    "drosophila": "dme",
    "worm": "cel",
    "celegans": "cel",
    "yeast": "sce",
    "ecoli": "eco",
}

REACTOME_ORGANISM_NAMES = {
    "human": "human",
    "mouse": "mouse",
    "rat": "rat",
    "zebrafish": "zebrafish",
    "fly": "fly",
    "drosophila": "fly",
    "worm": "celegans",
    "celegans": "celegans",
    "yeast": "yeast",
}

DEFAULT_MIN_GS_SIZE = 10
DEFAULT_MAX_GS_SIZE = 500
DEFAULT_SIMPLIFY_CUTOFF = 0.7

ORA_DIRECTIONS = ["up", "down", "combined"]

# ---------------------------------------------------------------------------
# Semantic-/term-similarity measures
#
# The SAME five GOSemSim measures (plus, for the network/emap step only, a
# sixth "JC" gene-overlap option) show up at TWO separate points in this
# pipeline:
#   1. clusterProfiler::simplify()   -- collapses redundant GO terms in a
#      single database's result table (_build_simplify_snippet below).
#   2. The term-similarity network feeding the "enrichment map" plot
#      (build_term_similarity_network below) -- this workspace's existing
#      implementation already computes one of these ("JC", gene-set
#      Jaccard overlap) natively in Python; the other five require a live
#      GOSemSim call (compute_go_semantic_similarity_matrix) since they
#      depend on GO's own graph structure/annotation statistics, which
#      only clusterProfiler/GOSemSim (R) has access to.
#
# Both call sites accept the SAME "Resnik"/"Lin"/"Rel"/"Jiang"/"Wang"
# measure keys so a user only has to learn one small vocabulary, even
# though which function actually consumes that measure differs.
# ---------------------------------------------------------------------------

# clusterProfiler's OWN default for simplify() has drifted across package
# versions -- "Rel" in older releases, "Wang" in current ones. Rather than
# relying on whatever the installed version happens to default to (which
# would silently change this app's behavior on a clusterProfiler upgrade),
# this app pins its own explicit default here.
DEFAULT_SIMPLIFY_MEASURE = "Wang"
DEFAULT_SIMILARITY_METHOD = "JC"

# Measures valid for BOTH simplify() and the network-similarity step.
SEMANTIC_SIMILARITY_MEASURES = ["Resnik", "Lin", "Rel", "Jiang", "Wang"]

# Measures that require GOSemSim to compute Information Content (IC) over
# the GO annotation corpus first (godata(..., computeIC = TRUE)). "Wang" is
# the only GOSemSim measure that does NOT need this -- it only uses GO
# graph topology -- so it can be built faster (computeIC = FALSE).
IC_BASED_SIMILARITY_MEASURES = ["Resnik", "Lin", "Rel", "Jiang"]

# Dropdown options for the simplify() measure picker (ontology_workspace.py).
SIMPLIFY_MEASURE_OPTIONS = {
    "Wang (graph-based -- fastest, current clusterProfiler default)": "Wang",
    "Rel / Relevance (Schlicker -- former clusterProfiler default)": "Rel",
    "Resnik": "Resnik",
    "Lin": "Lin",
    "Jiang & Conrath": "Jiang",
}

# Dropdown options for the network-similarity (enrichment map) method
# picker -- includes "JC", which is NOT valid for simplify().
SIMILARITY_METHOD_OPTIONS = {
    "Gene overlap / Jaccard (default -- fast, works for GO/KEGG/Reactome)": "JC",
    "Wang (GO semantic similarity, graph-based)": "Wang",
    "Rel / Relevance (GO semantic similarity)": "Rel",
    "Resnik (GO semantic similarity)": "Resnik",
    "Lin (GO semantic similarity)": "Lin",
    "Jiang & Conrath (GO semantic similarity)": "Jiang",
}

# Plain-language explanations shown next to each measure's picker in the
# UI (ontology_workspace.py) -- formula strings are plain-text/LaTeX-ish
# (renderable via st.latex if desired) rather than markdown, since they're
# meant to be read as math.
#
# "practical_effect": a SECOND, non-mathematical explanation focused on
# what a user should expect to SEE change in their results/plot if they
# switch TO this measure from another one -- distinct from "summary"
# (which explains the math) and "details" (which explains the theory) --
# added specifically so a non-statistician user has a concrete, actionable
# sense of "what happens if I pick this" rather than only a formula.
SIMILARITY_MEASURE_INFO = {
    "Resnik": {
        "category": "Information-content (IC) based",
        "requires_compute_ic": True,
        "formula": r"sim_{Resnik}(t_1, t_2) = IC(MICA)",
        "summary": "Similarity = information content of the terms' most informative common ancestor (MICA).",
        "details": (
            "The oldest, simplest information-content method. Looks only at how "
            "specific (rare) the SHARED ancestor term is -- it ignores how "
            "specific t1 and t2 themselves are, which is its well-known "
            "weakness (motivated the Lin/Rel corrections below)."
        ),
        "practical_effect": (
            "Tends to be the most GENEROUS measure -- it can call two terms "
            "\"similar\" even if they're fairly different in scope, as long as "
            "they share a reasonably specific ancestor. Expect MORE terms "
            "merged together by simplify(), or MORE connections drawn on the "
            "enrichment map, compared to Lin/Rel/Jiang."
        ),
        "citation": "Resnik, P. (1999). Semantic Similarity in a Taxonomy. J. Artif. Intell. Res., 11, 95-130.",
    },
    "Lin": {
        "category": "Information-content (IC) based",
        "requires_compute_ic": True,
        "formula": r"sim_{Lin}(t_1, t_2) = \dfrac{2 \cdot IC(MICA)}{IC(t_1) + IC(t_2)}",
        "summary": "Resnik's MICA score, normalized by how specific t1 and t2 individually are.",
        "details": (
            "Fixes Resnik's main blind spot by normalizing against the sum of "
            "the two terms' own information content, penalizing pairs where "
            "t1/t2 are very generic even if they share a specific ancestor."
        ),
        "practical_effect": (
            "STRICTER than Resnik -- two broad/generic terms will no longer "
            "look similar just because they share an ancestor. Expect FEWER "
            "terms merged by simplify(), and a SPARSER enrichment map, than "
            "with Resnik."
        ),
        "citation": "Lin, D. (1998). An Information-Theoretic Definition of Similarity. ICML.",
    },
    "Rel": {
        "category": "Information-content (IC) based",
        "requires_compute_ic": True,
        "formula": r"sim_{Rel}(t_1, t_2) = \dfrac{2 \cdot IC(MICA)\,(1 - p(MICA))}{IC(t_1) + IC(t_2)}",
        "summary": "Lin's formula plus a (1 - p(MICA)) term that further discounts overly generic ancestors.",
        "details": (
            "Proposed by Schlicker et al. (2006) as a hybrid of Resnik and Lin. "
            "This was clusterProfiler::simplify()'s default in older package "
            "versions (current versions default to Wang instead)."
        ),
        "practical_effect": (
            "Similar in spirit to Lin, but pushes similarity scores down "
            "further whenever the shared ancestor is a very COMMON GO term "
            "(e.g. \"biological_process\"). A safe middle-ground choice if "
            "you want IC-based results but distrust an overly generic shared "
            "ancestor driving the outcome."
        ),
        "citation": "Schlicker, A. et al. (2006). BMC Bioinformatics, 7, 302.",
    },
    "Jiang": {
        "category": "Information-content (IC) based",
        "requires_compute_ic": True,
        "formula": r"sim_{Jiang}(t_1, t_2) = 1 - \min(1,\; IC(t_1) + IC(t_2) - 2\cdot IC(MICA))",
        "summary": "A distance-based reformulation: penalizes the 'extra' information each term carries beyond its shared ancestor.",
        "details": (
            "Instead of a ratio (like Lin/Rel), sums each term's 'private' "
            "information content beyond the shared MICA and converts that "
            "distance back into a bounded similarity."
        ),
        "practical_effect": (
            "Behaves similarly to Lin/Rel in most cases, but reacts "
            "differently when t1 and t2 have very UNEQUAL specificity (e.g. "
            "one broad term, one narrow term) -- can produce a noticeably "
            "different set of \"redundant\" terms in that specific situation. "
            "Worth trying if Lin/Rel results look off for your data."
        ),
        "citation": "Jiang, J.J. & Conrath, D.W. (1997). ROCLING X.",
    },
    "Wang": {
        "category": "Graph-structure based",
        "requires_compute_ic": False,
        "formula": r"S_{GO}(A,B) = \dfrac{\sum_{t \in T_A \cap T_B} (S_A(t) + S_B(t))}{SV(A) + SV(B)}",
        "summary": "No information content needed -- uses only the GO graph topology and edge-type weights.",
        "details": (
            "The only measure here that skips GO annotation-frequency "
            "statistics entirely (no computeIC step needed) -- each term's "
            "'semantic value' is built by propagating weighted contributions "
            "from its ancestors in the GO graph. Fastest to prepare, and the "
            "current clusterProfiler::simplify() default."
        ),
        "practical_effect": (
            "Similarity is based purely on shared POSITION in the GO "
            "hierarchy, not on how common/rare terms are in your organism's "
            "annotations -- so results won't shift if your organism's "
            "annotation database is unusually sparse or dense. Generally "
            "gives a balanced, moderate amount of term-merging/connections -- "
            "a reasonable default if you're not sure which to pick."
        ),
        "citation": "Wang, J.Z. et al. (2007). Bioinformatics, 23(10), 1274-1281.",
    },
    "JC": {
        "category": "Gene-overlap based",
        "requires_compute_ic": False,
        "formula": r"J(t_1, t_2) = \dfrac{|Genes(t_1) \cap Genes(t_2)|}{|Genes(t_1) \cup Genes(t_2)|}",
        "summary": "Not a GOSemSim measure at all -- pure gene-set overlap (Jaccard). Works for GO, KEGG, and Reactome alike.",
        "details": (
            "Ignores GO graph structure entirely -- just measures how much two "
            "enriched terms' underlying gene sets overlap. This is what this "
            "workspace's enrichment map has always used, and needs no R/"
            "GOSemSim call at all, unlike the other five measures here."
        ),
        "practical_effect": (
            "Two terms connect ONLY if they were triggered by many of the "
            "SAME genes in YOUR specific result -- completely independent of "
            "GO's hierarchy/wording. This is the only measure that also "
            "works for KEGG/Reactome. It can miss connections between terms "
            "that are conceptually/hierarchically related but happen to be "
            "driven by different genes in your data."
        ),
        "citation": "Jaccard, P. (1912). New Phytologist, 11(2), 37-50.",
    },
}


def get_kegg_organism_code(species_key):
    return KEGG_ORGANISM_CODES.get(species_key)


def get_reactome_organism_name(species_key):
    return REACTOME_ORGANISM_NAMES.get(species_key)


def clusterprofiler_tools_available():
    return shutil.which("Rscript") is not None


DATABASE_OPTIONS = ["GO", "KEGG", "Reactome"]
GO_ONTOLOGY_OPTIONS = {
    "Biological Process (BP)": "BP",
    "Molecular Function (MF)": "MF",
    "Cellular Component (CC)": "CC",
    "All three (BP + MF + CC)": "ALL",
}


def databases_available_for_species(species_key):
    return {
        "GO": True,
        "KEGG": get_kegg_organism_code(species_key) is not None,
        "Reactome": get_reactome_organism_name(species_key) is not None,
    }


def simplify_available_for_ontology(go_ontology):
    return go_ontology != "ALL"


_gosemsim_availability_cache = {}


def check_gosemsim_available(force_recheck=False):
    if not force_recheck and "result" in _gosemsim_availability_cache:
        return _gosemsim_availability_cache["result"]

    if not clusterprofiler_tools_available():
        result = (False, "Rscript was not found on this system, so GOSemSim's availability could not be checked.")
        _gosemsim_availability_cache["result"] = result
        return result

    check_script = (
        'suppressMessages(ok <- requireNamespace("GOSemSim", quietly = TRUE)); '
        'if (ok) { '
        'cat("AVAILABLE version=", as.character(utils::packageVersion("GOSemSim")), "\\n", sep="") '
        '} else { '
        'cat("NOT_AVAILABLE\\n") '
        '}'
    )
    try:
        proc = subprocess.run(
            ["Rscript", "-e", check_script],
            capture_output=True, text=True, timeout=60,
        )
        output = (proc.stdout or "") + (proc.stderr or "")
        if "AVAILABLE version=" in output:
            version_match = re.search(r"AVAILABLE version=([\w.\-]+)", output)
            version = version_match.group(1) if version_match else "unknown version"
            outcome = (True, f"GOSemSim is installed (version {version}).")
        elif "NOT_AVAILABLE" in output:
            outcome = (False, "GOSemSim does not appear to be installed in this R environment.")
        else:
            outcome = (False, f"Could not determine GOSemSim availability from Rscript's output: {output.strip()[:300]!r}")
    except subprocess.TimeoutExpired:
        outcome = (False, "Timed out while checking GOSemSim availability.")
    except OSError as e:
        outcome = (False, f"Error while checking GOSemSim availability: {e}")

    _gosemsim_availability_cache["result"] = outcome
    return outcome


def verify_gosemsim_functional(orgdb_package, go_ontology="BP", timeout=180):
    if not clusterprofiler_tools_available():
        return False, "Rscript was not found on this system, so this check could not be run."

    script = f'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
tryCatch({{
  suppressMessages(library(GOSemSim))
  toy_universe <- AnnotationDbi::keys({orgdb_package}, keytype = "ENTREZID")
  toy_genes <- head(toy_universe, 50)
  toy_result <- suppressMessages(enrichGO(
    gene = toy_genes, universe = toy_universe, OrgDb = {orgdb_package},
    keyType = "ENTREZID", ont = "{go_ontology}",
    pvalueCutoff = 1, qvalueCutoff = 1
  ))
  if (!is.null(toy_result) && nrow(as.data.frame(toy_result)) > 1) {{
    invisible(clusterProfiler::simplify(toy_result, cutoff = 0.7, by = "p.adjust", select_fun = min))
    cat("FUNCTIONAL\\n")
  }} else {{
    cat("INCONCLUSIVE: the toy gene list did not produce enough enriched terms to meaningfully test simplify() in this environment.\\n")
  }}
}}, error = function(e) {{
  cat(paste0("NOT_FUNCTIONAL: ", conditionMessage(e)), "\\n")
}})
'''
    try:
        proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, timeout=timeout)
        output = (proc.stdout or "") + (proc.stderr or "")
        if "FUNCTIONAL" in output and "NOT_FUNCTIONAL" not in output:
            return True, "GOSemSim successfully simplified a live test result for this organism/ontology."
        if "INCONCLUSIVE" in output:
            detail_match = re.search(r"INCONCLUSIVE: (.*)", output)
            detail = detail_match.group(1) if detail_match else output.strip()[:500]
            return False, f"Inconclusive: {detail}"
        detail_match = re.search(r"NOT_FUNCTIONAL: (.*)", output)
        detail = detail_match.group(1) if detail_match else output.strip()[:500]
        return False, f"GOSemSim failed a live functional test: {detail}"
    except subprocess.TimeoutExpired:
        return False, (
            f"Functional check timed out after {timeout}s -- this can happen on "
            "the very first check while GOSemSim builds/caches its semantic "
            "similarity data for this organism. Trying again may be faster."
        )
    except OSError as e:
        return False, f"Error while running functional check: {e}"


def _direction_from_result_var(result_var):
    if result_var is None:
        return None
    if result_var.endswith("_up"):
        return "up"
    if result_var.endswith("_down"):
        return "down"
    return None


def parse_simplify_outcomes_from_log(log_text):
    outcomes = []
    for line in log_text.splitlines():
        simplified_match = re.search(r"GO terms simplified for (go_result\w*): (\d+) -> (\d+) term", line)
        if simplified_match:
            result_var, n_before, n_after = simplified_match.groups()
            outcomes.append({
                "outcome": "simplified", "detail": line.strip(),
                "direction": _direction_from_result_var(result_var),
                "n_before": int(n_before), "n_after": int(n_after),
            })
            continue
        if "Note: could not simplify GO results for " in line:
            match = re.search(r"for (go_result\w*) \(GOSemSim", line)
            result_var = match.group(1) if match else None
            outcomes.append({
                "outcome": "fallback", "detail": line.strip(),
                "direction": _direction_from_result_var(result_var),
                "n_before": None, "n_after": None,
            })
    return outcomes


def save_simplify_status(output_dir, outcomes):
    path = os.path.join(output_dir, "simplify_status.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(path, "w") as f:
        json.dump(outcomes, f, indent=2)


def load_simplify_status(output_dir):
    path = os.path.join(output_dir, "simplify_status.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _build_semdata_setup_snippet(orgdb_package, go_ontology, measure, semdata_var):
    """
    Build the R snippet that prepares the ONE shared GOSemSim semData
    object a script run's simplify() call(s) will use, assigned to
    semdata_var (e.g. ".simplify_semdata") -- built ONCE per script run
    (not once per direction block), since up/down-split ORA runs
    otherwise would have redundantly rebuilt this twice.

    "Wang" does not need Information Content (IC) statistics at all, so
    it's built with computeIC = FALSE (materially faster); the other
    four measures (Resnik/Lin/Rel/Jiang) require computeIC = TRUE.

    Wrapped in tryCatch so a missing/broken GOSemSim install degrades
    to semdata_var being NULL, rather than aborting the whole script --
    _build_simplify_snippet's own tryCatch below then reports this
    (via its existing "Note: could not simplify..." fallback message)
    without affecting any other database/direction block.
    """
    needs_ic = "TRUE" if measure in IC_BASED_SIMILARITY_MEASURES else "FALSE"
    return f'''
{semdata_var} <- tryCatch({{
  suppressMessages(library(GOSemSim))
  GOSemSim::godata("{orgdb_package}", ont = "{go_ontology}", computeIC = {needs_ic})
}}, error = function(e) {{
  cat(paste("Note: could not prepare GOSemSim data for simplify() (measure={measure}) --", conditionMessage(e)), "\\n")
  NULL
}})
'''


def _build_simplify_snippet(result_var, go_ontology, simplify_go, simplify_cutoff,
                             measure=DEFAULT_SIMPLIFY_MEASURE, semdata_var=".simplify_semdata"):
    """
    measure: one of SEMANTIC_SIMILARITY_MEASURES ("Resnik", "Lin", "Rel",
        "Jiang", "Wang") -- passed straight through to
        clusterProfiler::simplify()'s own measure= argument. Explicitly
        pinned rather than left to simplify()'s own (version-dependent)
        default -- see DEFAULT_SIMPLIFY_MEASURE's comment above.
    semdata_var: name of the R variable holding the shared semData
        object built once per script run by _build_semdata_setup_snippet
        -- reused across every simplify() call in this run rather than
        rebuilt per direction/database block.
    """
    if not (simplify_go and simplify_available_for_ontology(go_ontology)):
        return ""
    return f'''
if (!is.null({result_var}) && nrow(as.data.frame({result_var})) > 0) {{
  tryCatch({{
    if (is.null({semdata_var})) stop("GOSemSim semantic data is not available for measure '{measure}' (see earlier note).")
    n_before_simplify <- nrow(as.data.frame({result_var}))
    {result_var} <- clusterProfiler::simplify({result_var}, cutoff = {simplify_cutoff}, by = "p.adjust", select_fun = min, measure = "{measure}", semData = {semdata_var})
    n_after_simplify <- nrow(as.data.frame({result_var}))
    cat(paste("GO terms simplified for {result_var}:", n_before_simplify, "->", n_after_simplify, "term(s)"), "\\n")
  }}, error = function(e) {{
    cat(paste("Note: could not simplify GO results for {result_var} (GOSemSim may not be installed) --", conditionMessage(e)), "\\n")
  }})
}}
'''


def _wrap_with_direction_guard(body, gene_var, out_name, database_label):
    return f'''
if (length({gene_var}) >= 1) {{
{body}
}} else {{
  cat("Skipping {database_label} ({out_name}): no significant genes in this direction at the chosen threshold.\\n")
}}
'''


def _build_go_ora_block(result_var, gene_var, out_name, orgdb_package, go_ontology,
                         min_gs_size, max_gs_size, simplify_go, simplify_cutoff, guarded,
                         simplify_measure=DEFAULT_SIMPLIFY_MEASURE, semdata_var=".simplify_semdata"):
    simplify_snippet = _build_simplify_snippet(
        result_var, go_ontology, simplify_go, simplify_cutoff,
        measure=simplify_measure, semdata_var=semdata_var,
    )
    body = f'''
{result_var} <- enrichGO(
  gene = {gene_var}, universe = universe, OrgDb = {orgdb_package},
  keyType = "ENTREZID", ont = "{go_ontology}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size},
  readable = TRUE
)
{simplify_snippet}
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "GO") if guarded else body


def _build_kegg_ora_block(result_var, gene_var, out_name, kegg_organism,
                           min_gs_size, max_gs_size, guarded):
    body = f'''
{result_var} <- enrichKEGG(
  gene = {gene_var}, universe = universe, organism = "{kegg_organism}",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "KEGG") if guarded else body


def _build_reactome_ora_block(result_var, gene_var, out_name, reactome_organism,
                               min_gs_size, max_gs_size, guarded):
    body = f'''
{result_var} <- enrichPathway(
  gene = {gene_var}, universe = universe, organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, readable = TRUE
)
write_result({result_var}, "{out_name}")
'''
    return _wrap_with_direction_guard(body, gene_var, out_name, "Reactome") if guarded else body


_R_ORA_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_dir <- args[2]

df <- read.csv(input_path, stringsAsFactors = FALSE)
df <- df[!is.na(df$gene_id) & !is.na(df$padj) & !is.na(df$log2FoldChange), ]

from_type <- "{from_type}"

if (from_type == "ENTREZID") {{
  df$ENTREZID <- df$gene_id
}} else {{
  suppressMessages(
    conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
  )
  df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
}}
df <- df[!duplicated(df$ENTREZID), ]

universe <- unique(df$ENTREZID)

{gene_list_block}

{semdata_setup_block}

write_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched terms found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "term(s)"), "\n")
}}

{run_blocks}

cat("ORA analysis completed.\n")
'''

_GENE_LIST_BLOCK_COMBINED = r'''
sig_genes <- unique(df$ENTREZID[df$padj < {padj_threshold} & abs(df$log2FoldChange) >= {lfc_threshold}])
cat(paste("Significant genes (up + down combined):", length(sig_genes)), "\n")
if (length(sig_genes) < 1) {{
  stop("No significant genes at the chosen thresholds -- cannot run ORA. Try a less strict threshold.")
}}
'''

_GENE_LIST_BLOCK_SPLIT = r'''
sig_genes_up <- unique(df$ENTREZID[df$padj < {padj_threshold} & df$log2FoldChange >= {lfc_threshold}])
sig_genes_down <- unique(df$ENTREZID[df$padj < {padj_threshold} & df$log2FoldChange <= -{lfc_threshold}])
cat(paste("Up-regulated significant genes:", length(sig_genes_up)), "\n")
cat(paste("Down-regulated significant genes:", length(sig_genes_down)), "\n")
if (length(sig_genes_up) < 1 && length(sig_genes_down) < 1) {{
  stop("No significant genes (up OR down) at the chosen thresholds -- cannot run ORA. Try a less strict threshold.")
}}
'''


def build_ora_r_script(orgdb_package, from_type, padj_threshold, lfc_threshold,
                        run_go, go_ontology, run_kegg, kegg_organism,
                        run_reactome, reactome_organism,
                        min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                        simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                        simplify_measure=DEFAULT_SIMPLIFY_MEASURE,
                        split_by_direction=True):
    """
    simplify_measure: one of SEMANTIC_SIMILARITY_MEASURES ("Resnik",
        "Lin", "Rel", "Jiang", "Wang") -- see _build_simplify_snippet's
        docstring. Only relevant when simplify_go=True and run_go=True
        and go_ontology is a single sub-ontology.
    """
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    if split_by_direction:
        gene_list_block = _GENE_LIST_BLOCK_SPLIT.format(padj_threshold=padj_threshold, lfc_threshold=lfc_threshold)
        direction_specs = [("up", "sig_genes_up"), ("down", "sig_genes_down")]
    else:
        gene_list_block = _GENE_LIST_BLOCK_COMBINED.format(padj_threshold=padj_threshold, lfc_threshold=lfc_threshold)
        direction_specs = [(None, "sig_genes")]

    # Build the shared semData object ONCE for the whole script run
    # (reused by every direction block's simplify() call below) rather
    # than once per direction -- see _build_semdata_setup_snippet's
    # docstring. Only needed at all when GO + simplify are both active.
    semdata_var = ".simplify_semdata"
    needs_simplify_setup = run_go and simplify_go and simplify_available_for_ontology(go_ontology)
    semdata_setup_block = (
        _build_semdata_setup_snippet(orgdb_package, go_ontology, simplify_measure, semdata_var)
        if needs_simplify_setup else ""
    )

    run_blocks_parts = []
    for direction, gene_var in direction_specs:
        suffix = f"_{direction}" if direction else ""
        guarded = split_by_direction

        if run_go:
            result_var = f"go_result{suffix}"
            out_name = f"ora_GO{suffix}"
            run_blocks_parts.append(_build_go_ora_block(
                result_var, gene_var, out_name, orgdb_package, go_ontology,
                min_gs_size, max_gs_size, simplify_go, simplify_cutoff, guarded,
                simplify_measure=simplify_measure, semdata_var=semdata_var,
            ))
        if run_kegg:
            result_var = f"kegg_result{suffix}"
            out_name = f"ora_KEGG{suffix}"
            run_blocks_parts.append(_build_kegg_ora_block(
                result_var, gene_var, out_name, kegg_organism, min_gs_size, max_gs_size, guarded,
            ))
        if run_reactome:
            result_var = f"reactome_result{suffix}"
            out_name = f"ora_Reactome{suffix}"
            run_blocks_parts.append(_build_reactome_ora_block(
                result_var, gene_var, out_name, reactome_organism, min_gs_size, max_gs_size, guarded,
            ))

    return _R_ORA_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, gene_list_block=gene_list_block,
        semdata_setup_block=semdata_setup_block,
        run_blocks="\n".join(run_blocks_parts),
    )


_R_GSEA_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_dir <- args[2]

df <- read.csv(input_path, stringsAsFactors = FALSE)
df <- df[!is.na(df$gene_id) & !is.na(df$log2FoldChange), ]

from_type <- "{from_type}"

if (from_type == "ENTREZID") {{
  df$ENTREZID <- df$gene_id
}} else {{
  suppressMessages(
    conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
  )
  df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
}}
df <- df[!duplicated(df$ENTREZID), ]

gene_list <- df$log2FoldChange
names(gene_list) <- df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat(paste("Ranked gene list size:", length(gene_list)), "\n")

{semdata_setup_block}

write_gsea_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched gene sets found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "gene set(s)"), "\n")

  tryCatch({{
    suppressMessages(library(enrichplot))
    running_score_rows <- list()
    for (gs_id in res_df$ID) {{
      gs_info <- enrichplot:::gsInfo(result_obj, geneSetID = gs_id)
      gs_info$gene_set_id <- gs_id
      running_score_rows[[gs_id]] <- gs_info
    }}
    if (length(running_score_rows) > 0) {{
      combined <- do.call(rbind, running_score_rows)
      write.csv(combined, file.path(output_dir, paste0(out_name, "_running_score.csv")), row.names = FALSE)
    }}
  }}, error = function(e) {{
    cat(paste("Note: could not export running-score data for", out_name, "--", conditionMessage(e)), "\n")
  }})
}}

{run_go_block}
{run_kegg_block}
{run_reactome_block}

cat("GSEA analysis completed.\n")
'''

_R_GSEA_GO_BLOCK = r'''
go_result <- gseGO(
  geneList = gene_list, OrgDb = {orgdb_package}, keyType = "ENTREZID",
  ont = "{go_ontology}", pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size},
  eps = 0
)
{simplify_block}
write_gsea_result(go_result, "gsea_GO")
'''

_R_GSEA_KEGG_BLOCK = r'''
kegg_result <- gseKEGG(
  geneList = gene_list, organism = "{kegg_organism}", keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, eps = 0
)
write_gsea_result(kegg_result, "gsea_KEGG")
'''

_R_GSEA_REACTOME_BLOCK = r'''
reactome_result <- gsePathway(
  geneList = gene_list, organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}, eps = 0
)
write_gsea_result(reactome_result, "gsea_Reactome")
'''


def build_gsea_r_script(orgdb_package, from_type, run_go, go_ontology,
                         run_kegg, kegg_organism, run_reactome, reactome_organism,
                         min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                         simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                         simplify_measure=DEFAULT_SIMPLIFY_MEASURE):
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    semdata_var = ".simplify_semdata"
    needs_simplify_setup = run_go and simplify_go and simplify_available_for_ontology(go_ontology)
    semdata_setup_block = (
        _build_semdata_setup_snippet(orgdb_package, go_ontology, simplify_measure, semdata_var)
        if needs_simplify_setup else ""
    )

    simplify_block = _build_simplify_snippet(
        "go_result", go_ontology, simplify_go, simplify_cutoff,
        measure=simplify_measure, semdata_var=semdata_var,
    )

    run_go_block = _R_GSEA_GO_BLOCK.format(
        orgdb_package=orgdb_package, go_ontology=go_ontology,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_block=simplify_block,
    ) if run_go else ""
    run_kegg_block = _R_GSEA_KEGG_BLOCK.format(
        kegg_organism=kegg_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_kegg else ""
    run_reactome_block = _R_GSEA_REACTOME_BLOCK.format(
        reactome_organism=reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_reactome else ""

    return _R_GSEA_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, semdata_setup_block=semdata_setup_block,
        run_go_block=run_go_block, run_kegg_block=run_kegg_block,
        run_reactome_block=run_reactome_block,
    )


_R_COMPARE_CLUSTER_SCRIPT_TEMPLATE = r'''
suppressMessages(library(clusterProfiler))
suppressMessages(library({orgdb_package}))
{reactome_library_line}

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
output_dir <- args[2]

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

from_type <- "{from_type}"

gene_lists <- list()
universes <- c()

for (i in seq_len(nrow(manifest))) {{
  contrast_name <- manifest$contrast_name[i]
  input_path <- manifest$input_path[i]

  df <- read.csv(input_path, stringsAsFactors = FALSE)
  df <- df[!is.na(df$gene_id) & !is.na(df$padj) & !is.na(df$log2FoldChange), ]

  if (from_type == "ENTREZID") {{
    df$ENTREZID <- df$gene_id
  }} else {{
    suppressMessages(
      conv <- bitr(df$gene_id, fromType = from_type, toType = "ENTREZID", OrgDb = {orgdb_package})
    )
    df <- merge(df, conv, by.x = "gene_id", by.y = from_type)
  }}
  df <- df[!duplicated(df$ENTREZID), ]

  sig_genes <- unique(df$ENTREZID[df$padj < {padj_threshold} & abs(df$log2FoldChange) >= {lfc_threshold}])
  gene_lists[[contrast_name]] <- sig_genes
  universes <- c(universes, df$ENTREZID)
}}

universe <- unique(universes)
cat(paste("Contrasts compared:", length(gene_lists)), "\n")
for (nm in names(gene_lists)) {{
  cat(paste(" -", nm, ":", length(gene_lists[[nm]]), "significant gene(s)"), "\n")
}}

write_cc_result <- function(result_obj, out_name) {{
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {{
    cat(paste("No enriched terms found for:", out_name), "\n")
    return(invisible(NULL))
  }}
  res_df <- as.data.frame(result_obj)
  write.csv(res_df, file.path(output_dir, paste0(out_name, ".csv")), row.names = FALSE)
  cat(paste("Saved:", out_name, "-", nrow(res_df), "row(s)"), "\n")
}}

{run_go_block}
{run_kegg_block}
{run_reactome_block}

cat("compareCluster analysis completed.\n")
'''

_R_CC_GO_BLOCK = r'''
go_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichGO", universe = universe,
  OrgDb = {orgdb_package}, keyType = "ENTREZID", ont = "{go_ontology}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(go_cc, "compareCluster_GO")
'''

_R_CC_KEGG_BLOCK = r'''
kegg_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichKEGG", universe = universe,
  organism = "{kegg_organism}", keyType = "ncbi-geneid",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(kegg_cc, "compareCluster_KEGG")
'''

_R_CC_REACTOME_BLOCK = r'''
reactome_cc <- compareCluster(
  geneClusters = gene_lists, fun = "enrichPathway", universe = universe,
  organism = "{reactome_organism}",
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = {min_gs_size}, maxGSSize = {max_gs_size}
)
write_cc_result(reactome_cc, "compareCluster_Reactome")
'''


def build_compare_cluster_r_script(orgdb_package, from_type, padj_threshold, lfc_threshold,
                                    run_go, go_ontology, run_kegg, kegg_organism,
                                    run_reactome, reactome_organism,
                                    min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE):
    reactome_library_line = "suppressMessages(library(ReactomePA))" if run_reactome else ""

    run_go_block = _R_CC_GO_BLOCK.format(
        orgdb_package=orgdb_package, go_ontology=go_ontology,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_go else ""
    run_kegg_block = _R_CC_KEGG_BLOCK.format(
        kegg_organism=kegg_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_kegg else ""
    run_reactome_block = _R_CC_REACTOME_BLOCK.format(
        reactome_organism=reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    ) if run_reactome else ""

    return _R_COMPARE_CLUSTER_SCRIPT_TEMPLATE.format(
        orgdb_package=orgdb_package, reactome_library_line=reactome_library_line,
        from_type=from_type, padj_threshold=padj_threshold, lfc_threshold=lfc_threshold,
        run_go_block=run_go_block, run_kegg_block=run_kegg_block,
        run_reactome_block=run_reactome_block,
    )


def _run_r_script(script_text, script_path, r_args, timeout=1800):
    os.makedirs(os.path.dirname(script_path), exist_ok=True)
    with open(script_path, "w") as f:
        f.write(script_text)

    cmd = ["Rscript", script_path] + list(r_args)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
        return True, result.stdout + result.stderr
    except subprocess.CalledProcessError as e:
        return False, (e.stdout or "") + (e.stderr or "")
    except subprocess.TimeoutExpired:
        return False, f"Analysis timed out after {timeout // 60} minutes."


def run_ora_analysis(input_csv_path, output_dir, work_dir, orgdb_package, from_type,
                      padj_threshold, lfc_threshold, run_go, go_ontology,
                      run_kegg, kegg_organism, run_reactome, reactome_organism,
                      min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                      simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                      simplify_measure=DEFAULT_SIMPLIFY_MEASURE,
                      split_by_direction=True):
    os.makedirs(output_dir, exist_ok=True)
    script_text = build_ora_r_script(
        orgdb_package, from_type, padj_threshold, lfc_threshold,
        run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_go=simplify_go, simplify_cutoff=simplify_cutoff,
        simplify_measure=simplify_measure,
        split_by_direction=split_by_direction,
    )
    script_path = os.path.join(work_dir, "run_ora.R")
    return _run_r_script(script_text, script_path, [input_csv_path, output_dir])


def run_gsea_analysis(input_csv_path, output_dir, work_dir, orgdb_package, from_type,
                       run_go, go_ontology, run_kegg, kegg_organism,
                       run_reactome, reactome_organism,
                       min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE,
                       simplify_go=False, simplify_cutoff=DEFAULT_SIMPLIFY_CUTOFF,
                       simplify_measure=DEFAULT_SIMPLIFY_MEASURE):
    os.makedirs(output_dir, exist_ok=True)
    script_text = build_gsea_r_script(
        orgdb_package, from_type, run_go, go_ontology, run_kegg, kegg_organism,
        run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
        simplify_go=simplify_go, simplify_cutoff=simplify_cutoff,
        simplify_measure=simplify_measure,
    )
    script_path = os.path.join(work_dir, "run_gsea.R")
    return _run_r_script(script_text, script_path, [input_csv_path, output_dir])


def run_compare_cluster_analysis(contrast_input_paths, output_dir, work_dir, orgdb_package,
                                  from_type, padj_threshold, lfc_threshold, run_go, go_ontology,
                                  run_kegg, kegg_organism, run_reactome, reactome_organism,
                                  min_gs_size=DEFAULT_MIN_GS_SIZE, max_gs_size=DEFAULT_MAX_GS_SIZE):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(work_dir, exist_ok=True)

    manifest_path = os.path.join(work_dir, "compare_cluster_manifest.csv")
    manifest_df = pd.DataFrame({
        "contrast_name": list(contrast_input_paths.keys()),
        "input_path": list(contrast_input_paths.values()),
    })
    manifest_df.to_csv(manifest_path, index=False)

    script_text = build_compare_cluster_r_script(
        orgdb_package, from_type, padj_threshold, lfc_threshold,
        run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
        min_gs_size=min_gs_size, max_gs_size=max_gs_size,
    )
    script_path = os.path.join(work_dir, "run_compare_cluster.R")
    return _run_r_script(script_text, script_path, [manifest_path, output_dir])


def _result_file_path(output_dir, analysis_type, database, direction=None):
    suffix = f"_{direction}" if direction else ""
    return os.path.join(output_dir, f"{analysis_type}_{database}{suffix}.csv")


def enrichment_result_exists(output_dir, analysis_type, database, direction=None):
    return os.path.exists(_result_file_path(output_dir, analysis_type, database, direction=direction))


def available_ora_directions(output_dir, database):
    available = []
    for direction in ("up", "down"):
        if enrichment_result_exists(output_dir, "ora", database, direction=direction):
            available.append(direction)
    if enrichment_result_exists(output_dir, "ora", database, direction=None):
        available.append("combined")
    return available


def read_enrichment_result(output_dir, analysis_type, database, direction=None):
    path = _result_file_path(output_dir, analysis_type, database, direction=direction)
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)

    if "GeneRatio" in df.columns and not pd.api.types.is_numeric_dtype(df["GeneRatio"]):
        def _parse_ratio(s):
            try:
                num, denom = str(s).split("/")
                return float(num) / float(denom) if float(denom) != 0 else 0.0
            except (ValueError, ZeroDivisionError):
                return None
        df["GeneRatio_numeric"] = df["GeneRatio"].apply(_parse_ratio)

    return df


def read_gsea_running_score(output_dir, database, gene_set_id):
    path = os.path.join(output_dir, f"gsea_{database}_running_score.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    subset = df[df["gene_set_id"].astype(str) == str(gene_set_id)]
    return subset if not subset.empty else None


def list_available_databases(output_dir, analysis_type):
    return [db for db in DATABASE_OPTIONS if enrichment_result_exists(output_dir, analysis_type, db)]


def list_available_ora_databases(output_dir):
    return [db for db in DATABASE_OPTIONS if available_ora_directions(output_dir, db)]


def build_combined_results(output_dir, analysis_type, databases, n_top_per_db=10, direction=None):
    read_direction = None if direction == "combined" else direction
    frames = []
    for db in databases:
        df = read_enrichment_result(output_dir, analysis_type, db, direction=read_direction)
        if df is None or df.empty:
            continue
        df = df.sort_values("p.adjust").head(n_top_per_db).copy()
        df["Database"] = db
        frames.append(df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def derive_category_column(df):
    result = df.copy()
    if "ONTOLOGY" in result.columns:
        result["Category"] = result["ONTOLOGY"].fillna(result["Database"])
    else:
        result["Category"] = result["Database"]
    return result


def _parse_gene_id_list(gene_id_str, sep="/"):
    if pd.isna(gene_id_str) or not str(gene_id_str).strip():
        return set()
    return set(str(gene_id_str).split(sep))


# ---------------------------------------------------------------------------
# GOSemSim-based term similarity (for the enrichment map, as an alternative
# to the gene-overlap/Jaccard similarity computed natively in Python below)
# ---------------------------------------------------------------------------
#
# In-memory cache (lives for the Streamlit process's lifetime, same pattern
# as _gosemsim_availability_cache above) keyed on the exact term list +
# organism/ontology/measure combination -- avoids re-spawning an R
# subprocess on every UI rerun (e.g. a user only tweaking a plot's color
# scale, not its term list or similarity method) for what would otherwise
# be an identical, expensive recomputation. A failed computation (None) is
# NOT cached, so a transient error (or installing GOSemSim mid-session)
# doesn't get "stuck" without a retry.
_semantic_similarity_cache = {}

_R_SEMANTIC_SIMILARITY_SCRIPT_TEMPLATE = r'''
suppressMessages(library(GOSemSim))

args <- commandArgs(trailingOnly = TRUE)
term_ids_path <- args[1]
output_path <- args[2]

term_ids <- readLines(term_ids_path)

semData <- {godata_call}

sim_matrix <- GOSemSim::mgoSim(term_ids, term_ids, semData = semData, measure = "{measure}", combine = NULL)

write.csv(sim_matrix, output_path, row.names = TRUE)
cat("Semantic similarity matrix computed.\n")
'''


def compute_go_semantic_similarity_matrix(term_ids, orgdb_package, go_ontology,
                                           measure="Wang", work_dir=None, timeout=300,
                                           use_cache=True):
    """
    Compute a full pairwise GOSemSim similarity matrix for term_ids (GO
    IDs, e.g. "GO:0006955") via a live Rscript call to
    GOSemSim::mgoSim(), for use as an ALTERNATIVE to this workspace's
    default gene-overlap/Jaccard term similarity (see
    build_term_similarity_network below).

    measure: one of SEMANTIC_SIMILARITY_MEASURES ("Resnik", "Lin",
        "Rel", "Jiang", "Wang") -- NOT "JC", which is the gene-overlap
        measure computed natively in Python and has no GOSemSim
        equivalent to call out to here.

    Only meaningful for GO term IDs specifically -- unlike gene-overlap
    similarity (which works for any database's terms, since it only
    needs each term's member gene list), these measures rely on GO's
    own graph structure/annotation statistics, so KEGG/Reactome term
    IDs cannot be passed here.

    Returns a pandas DataFrame (term_ids x term_ids) of pairwise
    similarity scores, or None if the computation failed for any
    reason (Rscript missing, GOSemSim not installed, timeout, etc.) --
    callers should treat None as "fall back to gene-overlap (JC)
    similarity instead" rather than a hard error, consistent with this
    module's general graceful-degradation philosophy for GOSemSim-
    dependent features (see check_gosemsim_available's docstring).
    """
    if not term_ids or not clusterprofiler_tools_available():
        return None

    cache_key = (tuple(sorted(term_ids)), orgdb_package, go_ontology, measure)
    if use_cache and cache_key in _semantic_similarity_cache:
        return _semantic_similarity_cache[cache_key]

    work_dir = work_dir or os.path.join(os.getcwd(), ".ontology_semsim_tmp")
    os.makedirs(work_dir, exist_ok=True)
    term_ids_path = os.path.join(work_dir, "term_ids.txt")
    output_path = os.path.join(work_dir, "sim_matrix.csv")
    with open(term_ids_path, "w") as f:
        f.write("\n".join(term_ids))

    needs_ic = "TRUE" if measure in IC_BASED_SIMILARITY_MEASURES else "FALSE"
    godata_call = f'GOSemSim::godata("{orgdb_package}", ont = "{go_ontology}", computeIC = {needs_ic})'
    script_text = _R_SEMANTIC_SIMILARITY_SCRIPT_TEMPLATE.format(measure=measure, godata_call=godata_call)
    script_path = os.path.join(work_dir, "compute_semantic_similarity.R")

    success, _log = _run_r_script(script_text, script_path, [term_ids_path, output_path], timeout=timeout)
    if not success or not os.path.exists(output_path):
        return None

    try:
        sim_df = pd.read_csv(output_path, index_col=0)
        sim_df.index = sim_df.index.astype(str)
        sim_df.columns = sim_df.columns.astype(str)
    except Exception:
        return None

    if use_cache:
        _semantic_similarity_cache[cache_key] = sim_df
    return sim_df


def build_term_similarity_network(result_df, n_top=30, similarity_threshold=0.2,
                                   gene_id_column="geneID", similarity_method="JC",
                                   orgdb_package=None, go_ontology=None):
    """
    similarity_method: "JC" (DEFAULT -- gene-overlap/Jaccard, computed
        directly in Python from each term's member gene list; this is
        this workspace's original, unchanged behavior, and the only
        option that works for GO, KEGG, AND Reactome results alike) or
        one of SEMANTIC_SIMILARITY_MEASURES ("Resnik", "Lin", "Rel",
        "Jiang", "Wang") -- these use GO's own graph structure/
        annotation statistics INSTEAD of gene overlap, via a live
        GOSemSim call (compute_go_semantic_similarity_matrix), and
        therefore ONLY apply to GO results. Callers MUST pass
        orgdb_package and go_ontology when using one of these five.

    If a GOSemSim-based computation is requested but fails for any
    reason (not installed, R error, timeout), this function gracefully
    falls back to the same gene-overlap/Jaccard similarity "JC" would
    have produced, rather than returning an empty network -- consistent
    with this module's GOSemSim graceful-degradation philosophy
    elsewhere (see check_gosemsim_available's docstring).
    """
    if result_df is None or result_df.empty:
        return pd.DataFrame(), pd.DataFrame(columns=["source", "target", "weight"])

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)

    nodes_cols = ["ID", "Description", "p.adjust"]
    if count_col:
        nodes_cols.append(count_col)
    nodes_df = df[[c for c in nodes_cols if c in df.columns]].copy()
    if count_col and count_col != "Count":
        nodes_df = nodes_df.rename(columns={count_col: "Count"})

    term_ids = df["ID"].tolist()

    if similarity_method != "JC" and similarity_method in SEMANTIC_SIMILARITY_MEASURES:
        if not (orgdb_package and go_ontology):
            raise ValueError(
                f"similarity_method='{similarity_method}' requires orgdb_package and "
                "go_ontology (GOSemSim measures only apply to GO results)."
            )
        sim_matrix = compute_go_semantic_similarity_matrix(term_ids, orgdb_package, go_ontology, measure=similarity_method)
        if sim_matrix is not None:
            edges = []
            for id_a, id_b in itertools.combinations(term_ids, 2):
                try:
                    score = float(sim_matrix.loc[id_a, id_b])
                except (KeyError, ValueError, TypeError):
                    continue
                if pd.notna(score) and score >= similarity_threshold:
                    edges.append({"source": id_a, "target": id_b, "weight": score})
            edges_df = pd.DataFrame(edges, columns=["source", "target", "weight"])
            return nodes_df, edges_df
        # else: fall through to gene-overlap (JC) below -- GOSemSim
        # computation failed/unavailable, degrade gracefully rather
        # than returning an empty network.

    gene_sets = {row["ID"]: _parse_gene_id_list(row.get(gene_id_column, "")) for _, row in df.iterrows()}

    edges = []
    for id_a, id_b in itertools.combinations(gene_sets.keys(), 2):
        set_a, set_b = gene_sets[id_a], gene_sets[id_b]
        if not set_a or not set_b:
            continue
        intersection = len(set_a & set_b)
        union = len(set_a | set_b)
        jaccard = intersection / union if union > 0 else 0.0
        if jaccard >= similarity_threshold:
            edges.append({"source": id_a, "target": id_b, "weight": jaccard})

    edges_df = pd.DataFrame(edges, columns=["source", "target", "weight"])
    return nodes_df, edges_df


def build_gene_concept_network(result_df, gene_fc_map=None, n_top=10, gene_id_column="geneID"):
    if result_df is None or result_df.empty:
        empty_edges = pd.DataFrame(columns=["source", "target"])
        return pd.DataFrame(), pd.DataFrame(columns=["gene_id", "log2FoldChange"]), empty_edges

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)

    term_cols = ["ID", "Description", "p.adjust"]
    if count_col:
        term_cols.append(count_col)
    term_nodes_df = df[[c for c in term_cols if c in df.columns]].copy()
    if count_col and count_col != "Count":
        term_nodes_df = term_nodes_df.rename(columns={count_col: "Count"})

    gene_fc_map = gene_fc_map or {}
    all_genes = set()
    edges = []
    for _, row in df.iterrows():
        genes = _parse_gene_id_list(row.get(gene_id_column, ""))
        all_genes.update(genes)
        for gene in genes:
            edges.append({"source": row["ID"], "target": gene})

    gene_nodes_df = pd.DataFrame({
        "gene_id": sorted(all_genes),
        "log2FoldChange": [gene_fc_map.get(g) for g in sorted(all_genes)],
    })
    edges_df = pd.DataFrame(edges, columns=["source", "target"])
    return term_nodes_df, gene_nodes_df, edges_df


def build_upset_data(result_df, n_top=10, gene_id_column="geneID"):
    if result_df is None or result_df.empty:
        return pd.DataFrame(columns=["combination", "size", "term_ids"])

    df = result_df.sort_values("p.adjust").head(n_top).reset_index(drop=True)
    term_ids = df["ID"].tolist()
    term_labels = dict(zip(df["ID"], df.get("Description", df["ID"])))
    gene_sets = {row["ID"]: _parse_gene_id_list(row.get(gene_id_column, "")) for _, row in df.iterrows()}

    all_genes = set().union(*gene_sets.values()) if gene_sets else set()
    membership_counts = {}
    for gene in all_genes:
        membership = tuple(sorted(tid for tid in term_ids if gene in gene_sets[tid]))
        if membership:
            membership_counts[membership] = membership_counts.get(membership, 0) + 1

    rows = []
    for membership, size in membership_counts.items():
        labels = [term_labels.get(tid, tid) for tid in membership]
        rows.append({
            "combination": " + ".join(labels),
            "size": size,
            "term_ids": membership,
        })

    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=["combination", "size", "term_ids"])
    return result.sort_values("size", ascending=False).reset_index(drop=True)


def build_ridge_plot_data(output_dir, database, result_df, n_top=10):
    if result_df is None or result_df.empty:
        return []

    top_df = result_df.sort_values("p.adjust").head(n_top)
    rows = []
    for _, row in top_df.iterrows():
        running = read_gsea_running_score(output_dir, database, row["ID"])
        if running is None:
            continue
        member_values = running.loc[running["position"] == 1, "geneList"].tolist()
        if len(member_values) < 2:
            continue
        rows.append({
            "id": row["ID"],
            "description": row.get("Description", row["ID"]),
            "NES": row.get("NES"),
            "p.adjust": row["p.adjust"],
            "values": member_values,
        })

    rows.sort(key=lambda r: (r["NES"] if r["NES"] is not None else 0))
    return rows


def summarize_enrichment_result(result_df, padj_threshold=0.05, qvalue_threshold=None):
    if result_df is None or result_df.empty:
        return "No terms were returned for this database."

    sig_mask = result_df["p.adjust"] < padj_threshold
    if qvalue_threshold is not None and "qvalue" in result_df.columns:
        sig_mask = sig_mask & (result_df["qvalue"] < qvalue_threshold)

    n_sig = sig_mask.sum()
    if n_sig == 0:
        qval_note = f" and qvalue < {qvalue_threshold}" if qvalue_threshold is not None else ""
        return (
            f"No terms reached significance at padj < {padj_threshold}{qval_note} "
            f"({len(result_df)} term(s) were tested in total). This can "
            "happen with a small significant-gene list, or if your genes "
            "genuinely don't share a strong common biological theme at "
            "this threshold -- consider trying GSEA instead, which uses "
            "every tested gene rather than requiring a hard significance "
            "cutoff first."
        )

    sig_df = result_df[sig_mask]
    top_row = sig_df.sort_values("p.adjust").iloc[0]
    return (
        f"**{n_sig:,}** term(s)/pathway(s) significant at padj < {padj_threshold}"
        + (f" and qvalue < {qvalue_threshold}" if qvalue_threshold is not None and "qvalue" in result_df.columns else "")
        + f" (out of {len(result_df):,} tested). Most significant: "
        f"**{top_row.get('Description', top_row.get('ID', 'unknown'))}** "
        f"(padj = {top_row['p.adjust']:.2e})."
    )
