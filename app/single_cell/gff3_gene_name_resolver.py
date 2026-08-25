"""
single_cell/gff3_gene_name_resolver.py

--- IMPORTANT: this file must live in the same directory as
    singlecell_workspace.py (i.e. repo/app/single_cell/), NOT alongside
    reference_manager.py in repo/app/ -- app.py's own sys.path.insert()
    setup only adds the single_cell/ subfolder onto the import path, so
    this module (imported as a plain top-level name, "import
    gff3_gene_name_resolver", by BOTH singlecell_workspace.py directly
    and reference_manager.py's own backfill_gene_names_from_gff3_tiered()
    via a deferred, function-body-level import) can only resolve
    correctly if it's placed there. ---

Standalone, self-contained tiered gene-name resolution for Ensembl-style
GFF3 files -- built to address a real, confirmed gap in this project's
existing GFF3-fallback reference pipeline (see reference_manager.py's
own backfill_gene_names_from_gff3(), which restores gene_name from
GFF3's "Name" attribute after gffread's default GTF conversion silently
drops it genome-wide).

--- The real gap this closes (2026-08-25) ---
backfill_gene_names_from_gff3() only recovers a gene_name where the
source GFF3 actually HAS a "Name" attribute for that gene. Ensembl's own
GFF3 format is confirmed (directly, via a real published example) to
ALSO carry a separate "description" attribute on the same gene line --
e.g.:
    ID=gene:ENSTGUG00000013637;Name=DCBLD2;biotype=protein_coding;
    description=discoidin, CUB and LCCL domain containing 2
    [Source:NCBI gene;Acc:100218192];gene_id=ENSTGUG00000013637
For genes where "Name" is genuinely absent (a real, known gap for some
novel/predicted loci, certain non-coding RNAs, and pseudogenes in
Ensembl's own annotation), "description" is very often STILL present
and contains real, useful information that was previously discarded
entirely -- the affected gene fell straight through to displaying its
bare Ensembl ID (or, prior to a separate earlier fix, a confusing
"gene:ENSG..." hybrid form) with no attempt to use this second,
independently-populated attribute at all.

--- Resolution tiers (confirmed with the user, 2026-08-25) ---
For each gene, in order:
  1. "Name" attribute (the real, curated gene symbol) -- HIGHEST
     confidence, used as-is.
  2. "description" attribute, if "Name" is absent -- a real, sourced
     piece of information (not a fabrication), but NOT a clean gene
     symbol -- typically a full free-text description (e.g.
     "discoidin, CUB and LCCL domain containing 2"), with any trailing
     "[Source:...]" annotation-provenance suffix stripped off, and any
     percent-encoded characters (e.g. "%2C" for a literal comma, "%3B"
     for a literal semicolon -- confirmed directly from a real Ensembl
     GFF3 example) properly decoded.
  3. The bare Ensembl gene ID itself (e.g. "ENSG00000251562", with any
     GFF3 "gene:" ID-attribute-style prefix stripped) -- the final,
     honest fallback when NEITHER of the above is available.

--- Explicitly OUT of scope here (deferred, 2026-08-25) ---
Attempting a THIRD real-symbol-recovery path via clusterProfiler::bitr()
(Ensembl ID -> gene symbol conversion, for model organisms with an
available org.*.eg.db annotation package) was explicitly discussed and
DEFERRED -- this is considered squarely an Ontology Analysis-time
concern (this project's own bulk RNA-seq Ontology Analysis workflow
already uses bitr() for exactly this kind of ID conversion), not
something this single-cell reference-setup/display fix should also take
on. If a gene falls through to tier 3 (bare ID) here, that is
considered an acceptable, honest end state for THIS context.

--- Summary/flagging (2026-08-25) ---
summarize_gene_name_resolution() counts how many genes in a GFF3 landed
in each tier, for surfacing as a plain-language summary at Step 5
(reference setup) immediately after a GFF3-fallback reference is
confirmed -- e.g. "38,214 of 41,063 genes (93.1%) have a curated gene
symbol; 2,451 (6.0%) will display using their Ensembl description;
398 (0.9%) have neither and will display their bare Ensembl ID."
"""
import os
import re
import urllib.parse

_ID_ATTR_PATTERN = re.compile(r'(?:^|;)ID=([^;]+)')
_GENE_ID_ATTR_PATTERN = re.compile(r'(?:^|;)gene_id=([^;]+)')
_NAME_ATTR_PATTERN = re.compile(r'(?:^|;)Name=([^;]+)')
_DESCRIPTION_ATTR_PATTERN = re.compile(r'(?:^|;)description=([^;]+)')

RESOLUTION_TIER_NAME = "name"
RESOLUTION_TIER_DESCRIPTION = "description"
RESOLUTION_TIER_ID_ONLY = "id_only"


def _strip_gene_prefix(gene_id):
    "Strip a GFF3 'gene:' ID-attribute-style prefix, if present -- e.g. 'gene:ENSG00000251562' -> 'ENSG00000251562'."
    return gene_id[len("gene:"):] if gene_id.startswith("gene:") else gene_id


def _strip_source_suffix(description):
    """
    Strip Ensembl's own trailing "[Source:...]" annotation-provenance
    suffix from a description value -- e.g.
    "discoidin, CUB and LCCL domain containing 2 [Source:NCBI gene;Acc:100218192]"
    -> "discoidin, CUB and LCCL domain containing 2"
    """
    return re.sub(r"\s*\[Source:[^\]]*\]\s*$", "", description).strip()


def resolve_gene_names_from_gff3(gff3_path, max_lines=None):
    """
    Parse an Ensembl-style GFF3 file directly and, for every gene-level
    feature found, resolve a display name using the three-tier
    priority described in this module's own docstring: Name ->
    description (Source-suffix-stripped, percent-decoded) -> bare
    Ensembl gene ID (gene:-prefix-stripped).

    Only "gene" feature-type lines (column 3) are considered.

    Returns a dict: {gene_id: {"name": str, "tier": str}}
      - gene_id: the gene's OWN ID, with any "gene:" prefix stripped.
      - "name": the resolved display name for that gene.
      - "tier": one of RESOLUTION_TIER_NAME, RESOLUTION_TIER_DESCRIPTION,
        RESOLUTION_TIER_ID_ONLY.
    """
    if not gff3_path or not os.path.isfile(gff3_path):
        return {}

    resolved = {}

    with open(gff3_path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            feature_type = fields[2].strip().lower()
            if feature_type != "gene":
                continue

            attributes = fields[8]

            gene_id_match = _GENE_ID_ATTR_PATTERN.search(attributes)
            if gene_id_match:
                gene_id = gene_id_match.group(1).strip()
            else:
                id_match = _ID_ATTR_PATTERN.search(attributes)
                if not id_match:
                    continue
                gene_id = _strip_gene_prefix(id_match.group(1).strip())

            name_match = _NAME_ATTR_PATTERN.search(attributes)
            if name_match and name_match.group(1).strip():
                resolved[gene_id] = {"name": name_match.group(1).strip(), "tier": RESOLUTION_TIER_NAME}
                continue

            description_match = _DESCRIPTION_ATTR_PATTERN.search(attributes)
            if description_match and description_match.group(1).strip():
                raw_description = urllib.parse.unquote(description_match.group(1).strip())
                cleaned_description = _strip_source_suffix(raw_description)
                if cleaned_description:
                    resolved[gene_id] = {"name": cleaned_description, "tier": RESOLUTION_TIER_DESCRIPTION}
                    continue

            resolved[gene_id] = {"name": gene_id, "tier": RESOLUTION_TIER_ID_ONLY}

    return resolved


def summarize_gene_name_resolution(gff3_path, max_lines=None):
    """
    Build a plain-language-ready summary of how many genes in gff3_path
    resolved to each tier.

    Returns a dict:
        {
            "total_genes": int,
            "name_count": int,
            "description_count": int,
            "id_only_count": int,
            "name_pct": float,
            "description_pct": float,
            "id_only_pct": float,
        }
    Returns None if gff3_path doesn't exist/has no genes at all.
    """
    resolved = resolve_gene_names_from_gff3(gff3_path, max_lines=max_lines)
    total = len(resolved)
    if total == 0:
        return None

    name_count = sum(1 for v in resolved.values() if v["tier"] == RESOLUTION_TIER_NAME)
    description_count = sum(1 for v in resolved.values() if v["tier"] == RESOLUTION_TIER_DESCRIPTION)
    id_only_count = sum(1 for v in resolved.values() if v["tier"] == RESOLUTION_TIER_ID_ONLY)

    return {
        "total_genes": total,
        "name_count": name_count,
        "description_count": description_count,
        "id_only_count": id_only_count,
        "name_pct": round(100 * name_count / total, 1),
        "description_pct": round(100 * description_count / total, 1),
        "id_only_pct": round(100 * id_only_count / total, 1),
    }


def build_gene_name_resolution_summary_message(summary):
    """
    Build a single plain-language markdown string from
    summarize_gene_name_resolution()'s own output -- for direct display
    in Step 5's UI. Returns None if summary is None (nothing to show).
    """
    if summary is None:
        return None

    total = summary["total_genes"]
    lines = [
        f"**Gene name resolution for this reference** ({total:,} gene(s) total):",
        f"- ✅ {summary['name_count']:,} ({summary['name_pct']}%) have a curated gene symbol.",
    ]
    if summary["description_count"] > 0:
        lines.append(
            f"- ℹ️ {summary['description_count']:,} ({summary['description_pct']}%) have no curated symbol, "
            "but will display using their Ensembl-provided description instead."
        )
    if summary["id_only_count"] > 0:
        lines.append(
            f"- ⚠️ {summary['id_only_count']:,} ({summary['id_only_pct']}%) have neither a symbol nor a "
            "description available, and will display their bare Ensembl gene ID."
        )
    return "\n".join(lines)
