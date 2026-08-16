"""
counts_matrix_manager.py

Streamlit-independent gene counts matrix assembly logic for the Bulk
RNA-Seq pipeline's final "combine into a gene counts matrix" step.

Extracted out of alignment_workspace.py so this logic can be called
both from that interactive workspace AND from a non-interactive
orchestrator (Advanced Mode / Monitor Mode), which need to merge
per-sample Salmon/STAR quantification output into one matrix, and
auto-detect a transcript-to-gene / gene-symbol mapping for the active
reference, without any Streamlit UI calls in the way.

This module owns:
  - merging per-sample Salmon quant.sf files into one wide counts
    matrix (transcript-level, or gene-level once tximport has already
    collapsed them -- tximport itself still runs via
    quantification_manager.run_tximport_gene_collapse)
  - merging per-sample STAR ReadsPerGene.out.tab files into one wide
    gene-level counts matrix
  - auto-detecting a transcript-to-gene (tx2gene) mapping for the
    active reference, so Salmon's transcript-level output can be
    collapsed to gene-level
  - auto-detecting a gene_id -> gene_symbol mapping for the active
    reference, so downstream Differential Expression can show readable
    gene symbols without a manual upload

It does NOT know about Streamlit, project_manager step-tracking, or
anything UI-related -- callers (alignment_workspace.py today, an
Advanced/Monitor Mode orchestrator in the future) are responsible for
deciding when to call these functions, showing progress/errors, saving
the resulting matrix/mapping to disk, and marking the project's
"counts_matrix_complete" step once a matrix has been written
successfully.
"""

import os

import pandas as pd

import project_manager as pm
import reference_manager as rm


def merge_salmon_counts(sample_names, salmon_quant_dir, count_column="NumReads"):
    """
    Merge per-sample Salmon quant.sf files into one wide gene/transcript
    counts matrix (rows = gene/transcript IDs, columns = samples).

    count_column: "NumReads" (raw estimated read counts, appropriate for
    DESeq2/differential expression) or "TPM" (normalized, better for
    direct cross-sample expression-level comparisons/visualization but
    not appropriate as DESeq2 input, which expects raw counts).

    Note on gene vs. transcript level: Salmon quantifies at the
    transcript level. For organisms with no introns/alternative splicing
    (e.g. bacteria like E. coli), each gene typically has exactly one
    transcript, so transcript-level and gene-level counts are
    effectively the same thing, and this direct merge is appropriate.
    For eukaryotes with multiple transcripts per gene, transcript-level
    counts should normally be collapsed to gene-level first (typically
    via tximport in R, using a transcript-to-gene mapping) before
    differential expression — this direct merge does NOT do that
    collapsing, so it will produce one row per transcript rather than
    per gene for such organisms. See get_tx2gene_mapping below and
    quantification_manager.run_tximport_gene_collapse for the
    gene-collapsing path.

    Returns (merged_df_or_None, missing_samples_list).
    """
    merged = None
    missing_samples = []

    for sample_name in sample_names:
        quant_path = os.path.join(salmon_quant_dir, sample_name, "quant.sf")
        if not os.path.exists(quant_path):
            missing_samples.append(sample_name)
            continue

        df = pd.read_csv(quant_path, sep="\t")
        df = df[["Name", count_column]].rename(columns={"Name": "gene_id", count_column: sample_name})

        merged = df if merged is None else merged.merge(df, on="gene_id", how="outer")

    if merged is not None and count_column == "NumReads":
        # Salmon's NumReads is fractional (due to resolving multi-mapping
        # reads probabilistically across transcripts) — round to whole
        # read counts, which is what DESeq2 and most downstream tools
        # expect.
        for col in merged.columns:
            if col != "gene_id":
                merged[col] = merged[col].fillna(0).round(0).astype(int)

    return merged, missing_samples


def merge_star_counts(sample_names, star_align_dir, strandedness_column=1):
    """
    Merge per-sample STAR ReadsPerGene.out.tab files into one wide gene
    counts matrix. This is already gene-level (STAR's --quantMode
    GeneCounts assigns reads directly to genes using the GTF, unlike
    Salmon's transcript-level output), so no separate gene-collapsing
    step is needed here.

    strandedness_column: which of STAR's three count columns to use —
        1 = unstranded (default; safest choice when the library prep's
            strandedness protocol is unknown)
        2 = forward-stranded
        3 = reverse-stranded
    STAR's ReadsPerGene.out.tab format has 4 columns (gene_id,
    unstranded, forward, reverse) and its first 4 rows are summary
    statistics (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
    rather than genes — these are skipped.

    Returns (merged_df_or_None, missing_samples_list).
    """
    col_map = {1: "unstranded", 2: "forward", 3: "reverse"}
    use_col = col_map.get(strandedness_column, "unstranded")

    merged = None
    missing_samples = []

    for sample_name in sample_names:
        path = os.path.join(star_align_dir, sample_name, f"{sample_name}_ReadsPerGene.out.tab")
        if not os.path.exists(path):
            missing_samples.append(sample_name)
            continue

        df = pd.read_csv(
            path, sep="\t", header=None, skiprows=4,
            names=["gene_id", "unstranded", "forward", "reverse"],
        )
        df = df[["gene_id", use_col]].rename(columns={use_col: sample_name})

        merged = df if merged is None else merged.merge(df, on="gene_id", how="outer")

    return merged, missing_samples


def get_tx2gene_mapping(project, reference_dir):
    """
    Attempt to build a tx2gene mapping for this project's reference,
    auto-detecting the correct source based on how the reference was
    obtained:
      - Preset Ensembl species (human/mouse/yeast): parsed directly from
        the cDNA FASTA's headers (no GTF needed).
      - Preset no-intron species (E. coli): identity mapping, since our
        own gene-level extraction already writes one sequence per gene.
      - Custom uploads: tries the GTF first (correct for eukaryotes
        processed via gffread); falls back to identity mapping if the
        user flagged their organism as having no introns.

    project: the project name, used only to look up which
        species/custom choice was saved for it (see
        project_manager.get_reference_choice) -- this function itself
        does not know about project_manager's directory layout beyond
        that.

    reference_dir: the EFFECTIVE reference directory already resolved
        by the caller -- reference_manager.shared_reference_dir(species_key)
        for a preset organism, or project_manager.reference_dir(project)
        for a custom upload.

    Returns (tx2gene_dict_or_None, source_description_str). An empty/None
    mapping means gene-level collapsing isn't available for this
    reference and the caller should fall back to the direct transcript-
    level merge.
    """
    species_key, is_custom = pm.get_reference_choice(project)

    if not is_custom and species_key:
        entry = rm.REFERENCE_CATALOG.get(species_key, {})
        # Preset species' downloaded files live one level deeper than
        # reference_dir itself, in a "cdna" subdirectory -- this is
        # because reference_dir here is the SHARED, project-independent
        # location (see project_manager.py's shared_reference_dir()),
        # and reference_manager.py's ensure_shared_resource() manages
        # that "cdna" subdirectory as one atomically-built-and-renamed
        # unit (so a concurrent reader never sees a half-downloaded
        # file). Custom (per-project) references have no such
        # subdirectory -- their files sit flat directly in
        # reference_dir, matching the original convention exactly.
        cdna_path = os.path.join(reference_dir, "cdna", f"{species_key}.cdna.fa")

        if entry.get("no_introns"):
            if os.path.exists(cdna_path):
                return rm.build_identity_tx2gene(cdna_path), "identity mapping (no-intron organism)"
        elif os.path.exists(cdna_path):
            tx2gene = rm.extract_tx2gene_from_ensembl_fasta(cdna_path)
            if tx2gene:
                return tx2gene, "parsed from Ensembl cDNA FASTA headers"

        return None, None

    # Custom species path
    custom_gtf_path = os.path.join(reference_dir, "custom_input.gtf")
    extracted_fasta_path = os.path.join(reference_dir, "custom_extracted_transcripts.fa")
    custom_fasta_path = os.path.join(reference_dir, "custom_input.fa")

    if os.path.exists(custom_gtf_path):
        tx2gene = rm.extract_tx2gene_from_gtf(custom_gtf_path)
        if tx2gene:
            return tx2gene, "parsed from your uploaded GTF/GFF3 annotation"

    # If flagged as no-intron and we have an extracted (identity) FASTA,
    # fall back to identity mapping.
    if os.path.exists(extracted_fasta_path):
        return rm.build_identity_tx2gene(extracted_fasta_path), "identity mapping (no-intron organism)"
    if os.path.exists(custom_fasta_path):
        return rm.build_identity_tx2gene(custom_fasta_path), "identity mapping (assuming transcript FASTA already provided)"

    return None, None


def get_gene_symbol_mapping(project, reference_dir, method):
    """
    Attempt to build a gene_id -> gene_symbol mapping for this project's
    reference, auto-detecting the correct source based on how the
    reference was obtained -- mirrors get_tx2gene_mapping's source
    detection above, but for human-readable gene symbols rather than
    transcript-to-gene collapsing:
      - Preset Ensembl species (human/mouse/yeast/Drosophila/
        C. elegans/zebrafish): parsed directly from the cDNA FASTA's
        "gene_symbol:" header field when that FASTA is already on disk
        (from the Salmon path); falls back to the downloaded genome GTF's
        "gene_name" attribute if only the STAR path was used and no cDNA
        FASTA was ever downloaded for this project.
      - Preset no-intron species (E. coli): identity mapping (gene_id
        used as its own symbol), since these organisms are identified by
        locus tag rather than a separate common name.
      - Custom uploads: parsed from the uploaded GTF/GFF3's "gene_name"
        attribute; falls back to identity mapping if no GTF is present
        or no gene_name attributes were found.

    method: "salmon" or "star" -- currently unused directly by this
        function's own logic (both paths probe the same candidate file
        locations regardless), but kept as an explicit parameter since
        callers already have it on hand and a future reference source
        may need to branch on it.

    Returns (gene_symbol_dict_or_None, source_description_str). A None
    mapping means no symbol source was available for this reference --
    the caller should skip saving a mapping file in that case, and the
    Differential Expression workspace will fall back to showing raw
    gene IDs.
    """
    species_key, is_custom = pm.get_reference_choice(project)

    if not is_custom and species_key:
        entry = rm.REFERENCE_CATALOG.get(species_key, {})
        # See get_tx2gene_mapping's matching comment above for why
        # preset species' files live in "cdna"/"genome" subdirectories
        # of the shared reference_dir, rather than flat inside it.
        cdna_path = os.path.join(reference_dir, "cdna", f"{species_key}.cdna.fa")
        gtf_path = os.path.join(reference_dir, "genome", f"{species_key}.annotation.gtf")

        if entry.get("no_introns"):
            if os.path.exists(cdna_path):
                return rm.build_identity_tx2gene(cdna_path), "identity mapping (no-intron organism)"
        elif os.path.exists(cdna_path):
            gene_symbol = rm.extract_gene_symbol_map_from_ensembl_fasta(cdna_path)
            if gene_symbol:
                return gene_symbol, "parsed from Ensembl cDNA FASTA headers"
        if os.path.exists(gtf_path):
            gene_symbol = rm.extract_gene_symbol_map_from_gtf(gtf_path)
            if gene_symbol:
                return gene_symbol, "parsed from downloaded genome annotation (GTF)"

        return None, None

    # Custom species path
    custom_gtf_path = os.path.join(reference_dir, "custom_input.gtf")
    extracted_fasta_path = os.path.join(reference_dir, "custom_extracted_transcripts.fa")
    custom_fasta_path = os.path.join(reference_dir, "custom_input.fa")

    if os.path.exists(custom_gtf_path):
        gene_symbol = rm.extract_gene_symbol_map_from_gtf(custom_gtf_path)
        if gene_symbol:
            return gene_symbol, "parsed from your uploaded GTF/GFF3 annotation"

    if os.path.exists(extracted_fasta_path):
        return rm.build_identity_tx2gene(extracted_fasta_path), "identity mapping (no gene_name found)"
    if os.path.exists(custom_fasta_path):
        return rm.build_identity_tx2gene(custom_fasta_path), "identity mapping (no gene_name found)"

    return None, None
