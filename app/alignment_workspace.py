"""
alignment_workspace.py

RNA Alignment & Counts workspace.

This picks up where trimming_workspace.py leaves off: it reuses the same
active project (trimmed reads already on disk) rather than requiring
anything to be re-uploaded.

Two alignment/quantification methods are offered, since they serve
genuinely different needs:

  - Salmon (pseudo-alignment against a transcriptome):
      Fast, low memory, no novel isoform detection. Needs only a
      transcriptome FASTA. Best for standard differential expression
      when you don't need novel splice junction/isoform discovery.
      Output plugs into tximport + DESeq2.

  - STAR (splice-aware genome alignment):
      Slower, memory-hungry (needs a full genome index), but detects
      novel splice junctions/isoforms and produces BAM files usable by
      other tools (e.g. IGV). Needs a genome FASTA + GTF annotation.
      Gene-level counts come directly from STAR's built-in
      --quantMode GeneCounts output (ReadsPerGene.out.tab per sample),
      avoiding the need for a separate counting tool like featureCounts.

Design goal: same as the other workspaces — assume the user has little
to no bioinformatics background, explain each step (and this choice) in
plain language.

This module is fully self-contained. All alignment/counting development
should happen here — editing this file has zero effect on
spatial_workspace.py, bulk_rnaseq_workspace.py, or trimming_workspace.py.
"""

import os

import pandas as pd
import streamlit as st

import project_manager as pm
import reference_manager as rm
import quantification_manager as qm
import file_browser as fb


# ---------------------------------------------------------------------------
# Step 4 helpers: merging per-sample quantification output into one matrix
# ---------------------------------------------------------------------------

def _merge_salmon_counts(sample_names, salmon_quant_dir, count_column="NumReads"):
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
    per gene for such organisms. See the "how to interpret this" help
    text shown in the UI for this caveat.

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


def _merge_star_counts(sample_names, star_align_dir, strandedness_column=1):
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


def _get_tx2gene_mapping(project, reference_dir):
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


def _get_gene_symbol_mapping(project, reference_dir, method):
    """
    Attempt to build a gene_id -> gene_symbol mapping for this project's
    reference, auto-detecting the correct source based on how the
    reference was obtained -- mirrors _get_tx2gene_mapping's source
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

    Returns (gene_symbol_dict_or_None, source_description_str). A None
    mapping means no symbol source was available for this reference --
    the caller should skip saving a mapping file in that case, and the
    Differential Expression workspace will fall back to showing raw
    gene IDs.
    """
    species_key, is_custom = pm.get_reference_choice(project)

    if not is_custom and species_key:
        entry = rm.REFERENCE_CATALOG.get(species_key, {})
        # See _get_tx2gene_mapping's matching comment above for why
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


def _render_counts_matrix_step(project, method, samplesheet_df, reference_dir):
    """
    Step 4 UI: merge per-sample quantification output (Salmon quant.sf
    or STAR ReadsPerGene.out.tab) into a single combined gene counts
    matrix, ready for differential expression.

    reference_dir: the EFFECTIVE reference directory already resolved
        by render() -- pm.shared_reference_dir(species_key) for a
        preset organism, or pm.reference_dir(project) for a custom
        upload. Passed in explicitly (rather than recomputed here)
        so this function correctly uses the SHARED location for
        preset species' tx2gene mapping/gene symbol mapping, instead
        of always looking in this project's own private reference
        folder (which would be empty/wrong for a preset species after
        the shared-reference change).
    """
    st.header("Step 4: Combine Into a Gene Counts Matrix")

    with st.expander("ℹ️ What is a gene counts matrix, and why do I need one? (click to learn more)"):
        st.markdown(
            "Right now you have **one separate result per sample** — "
            "each telling you how much each gene was expressed in that "
            "one sample. A **gene counts matrix** combines all of your "
            "samples into a single table: one row per gene, one column "
            "per sample. This is the standard input format for "
            "differential expression tools like DESeq2, which compare "
            "expression levels *across* samples to find genes that "
            "differ between your experimental conditions."
        )
        if method == "salmon":
            st.markdown(
                "\n**A note on gene vs. transcript level:** Salmon "
                "quantifies at the transcript level. For organisms with "
                "no introns (e.g. bacteria like E. coli), each gene "
                "usually has exactly one transcript, so this matrix is "
                "effectively already gene-level. For eukaryotes with "
                "multiple transcripts per gene, transcript counts "
                "technically need an extra collapsing step (normally "
                "done with `tximport` in R) to become proper gene-level "
                "counts — this tool can do that collapsing automatically "
                "below if a transcript-to-gene mapping can be determined "
                "for your reference."
            )

    quantification_done = pm.has_completed_step(project, "quantification_complete")
    if not quantification_done:
        st.info("Complete Step 3 (run quantification) to unlock this step.")
        return

    sample_names = sorted(samplesheet_df["sample"].astype(str).tolist())

    counts_matrix_path = pm.counts_matrix_path(project)
    matrix_already_built = os.path.exists(counts_matrix_path)

    count_type_help = (
        "Raw estimated read counts (recommended for differential "
        "expression tools like DESeq2) vs. TPM (a normalized value "
        "that's easier to compare directly across samples/genes for "
        "visualization, but is NOT appropriate as input to DESeq2, "
        "which expects raw counts)."
    )

    use_tximport = False
    tx2gene = None
    tx2gene_source = None
    # reference_dir is now a parameter passed in by render() (already
    # correctly resolved to the shared or custom-per-project location)
    # rather than recomputed here -- see this function's docstring.

    if method == "salmon":
        count_type = st.radio(
            "Which value should the matrix contain?",
            ["Raw read counts (for DESeq2)", "TPM (normalized, for visualization only)"],
            help=count_type_help,
            key="counts_matrix_type_radio",
        )
        count_column = "NumReads" if count_type.startswith("Raw") else "TPM"

        # Gene-level collapsing via tximport is only offered for raw
        # counts (the DESeq2-relevant path) — TPM is typically viewed at
        # whatever level the user finds useful, and tximport's own TPM
        # summarization has different semantics we don't replicate here.
        if count_column == "NumReads":
            tx2gene, tx2gene_source = _get_tx2gene_mapping(project, reference_dir)

            if tx2gene:
                collapse_help = (
                    "Recommended for eukaryotes with multiple transcripts "
                    "per gene. Uses tximport (R) to properly sum "
                    "transcript-level counts up to gene-level using a "
                    "transcript-to-gene mapping. For organisms with no "
                    "introns (like E. coli), this has no effect since "
                    "each gene already has exactly one transcript."
                )
                use_tximport = st.checkbox(
                    "Collapse to gene-level using tximport (recommended)",
                    value=True,
                    help=collapse_help,
                    key="use_tximport_checkbox",
                )
                if use_tximport:
                    st.caption(f"tx2gene mapping source: {tx2gene_source}")
            else:
                st.info(
                    "ℹ️ A transcript-to-gene mapping couldn't be determined "
                    "automatically for this reference, so the matrix below "
                    "will be at the transcript level rather than "
                    "gene-level. This is fine for no-intron organisms, but "
                    "for eukaryotes this means some genes may appear as "
                    "multiple rows (one per transcript)."
                )
    else:
        count_column = None  # STAR's ReadsPerGene.out.tab is always raw counts

    build_label = "🔄 Re-build Counts Matrix" if matrix_already_built else "📊 Build Counts Matrix"
    if matrix_already_built and not st.session_state.get("_counts_matrix_clicked"):
        st.success("✅ A counts matrix has already been built for this project.")

    if st.button(build_label, key="build_counts_matrix_btn"):
        st.session_state["_counts_matrix_clicked"] = True
        missing_samples = []
        matrix_df = None

        if method == "salmon" and use_tximport and tx2gene:
            quant_dir = pm.salmon_quant_dir(project)
            sample_quant_paths = {}
            for sample_name in sample_names:
                quant_path = os.path.join(quant_dir, sample_name, "quant.sf")
                if os.path.exists(quant_path):
                    sample_quant_paths[sample_name] = quant_path
                else:
                    missing_samples.append(sample_name)

            if sample_quant_paths:
                tx2gene_dir = os.path.join(reference_dir, "tx2gene")
                tx2gene_path = rm.save_tx2gene_csv(tx2gene, os.path.join(tx2gene_dir, "tx2gene.csv"))
                work_dir = os.path.join(reference_dir, "tximport_work")

                with st.spinner("Running tximport to collapse transcripts to gene-level..."):
                    success, log = qm.run_tximport_gene_collapse(
                        sample_quant_paths, tx2gene_path, counts_matrix_path, work_dir
                    )

                if not success:
                    st.error("tximport failed. Details below:")
                    st.code(log)
                    st.info(
                        "💡 Falling back to a direct transcript-level merge "
                        "instead. You can uncheck the tximport option above "
                        "and rebuild if you'd prefer that outcome directly."
                    )
                    matrix_df, extra_missing = _merge_salmon_counts(sample_names, quant_dir, count_column=count_column)
                    missing_samples.extend(extra_missing)
                else:
                    st.success(f"✅ {log.strip().splitlines()[-1] if log.strip() else 'tximport completed.'}")
                    matrix_df = pd.read_csv(counts_matrix_path)

        elif method == "salmon":
            quant_dir = pm.salmon_quant_dir(project)
            matrix_df, missing_samples = _merge_salmon_counts(sample_names, quant_dir, count_column=count_column)
        else:
            align_dir = pm.star_align_dir(project)
            matrix_df, missing_samples = _merge_star_counts(sample_names, align_dir)

        if missing_samples:
            st.warning(
                f"⚠️ These samples were skipped because their "
                f"quantification output wasn't found: **{', '.join(missing_samples)}**. "
                "Re-run Step 3 for these samples if this is unexpected."
            )

        if matrix_df is None or matrix_df.empty:
            st.error(
                "⚠️ No quantification output could be found for any "
                "sample. Please make sure Step 3 completed successfully "
                "before building the counts matrix."
            )
        else:
            os.makedirs(os.path.dirname(counts_matrix_path), exist_ok=True)
            matrix_df.to_csv(counts_matrix_path, index=False)
            pm.mark_step_complete(project, "counts_matrix_complete")

            # Auto-build a gene_id -> gene_symbol mapping for this
            # project's reference/species, so the Differential Expression
            # workspace's volcano plot and results tables show readable
            # gene symbols automatically -- no manual upload required.
            # This is best-effort: if no symbol source is available
            # (e.g. a custom reference with no gene_name in its GTF), no
            # mapping file is written and DE falls back to raw gene IDs.
            gene_symbol_map, gene_symbol_source = _get_gene_symbol_mapping(project, reference_dir, method)
            if gene_symbol_map:
                rm.save_gene_symbol_map_csv(gene_symbol_map, pm.gene_symbol_map_path(project))
                pm.save_gene_id_mapping_meta(project, {
                    "source": "auto_parse",
                    "detail": gene_symbol_source,
                })
                st.caption(f"🏷️ Gene symbol mapping ready ({gene_symbol_source}) -- will be used automatically in Differential Expression. For genes that fall back to their raw reference ID, the Differential Expression workspace can convert them via a proper annotation database lookup.")

            n_genes = len(matrix_df)
            n_samples = len(matrix_df.columns) - 1  # minus the gene_id column
            st.success(f"✅ Counts matrix built: {n_genes:,} gene(s) × {n_samples} sample(s).")
            st.dataframe(matrix_df, use_container_width=True, hide_index=True)

            csv_bytes = matrix_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download Counts Matrix (.csv)",
                data=csv_bytes,
                file_name="gene_counts_matrix.csv",
                mime="text/csv",
            )

    elif matrix_already_built:
        # Show the previously built matrix even if the button wasn't
        # just clicked this run, so reopening the project displays it.
        existing_df = pd.read_csv(counts_matrix_path)
        n_genes = len(existing_df)
        n_samples = len(existing_df.columns) - 1
        st.caption(f"Current matrix: {n_genes:,} gene(s) × {n_samples} sample(s).")
        st.dataframe(existing_df, use_container_width=True, hide_index=True)

        csv_bytes = existing_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download Counts Matrix (.csv)",
            data=csv_bytes,
            file_name="gene_counts_matrix.csv",
            mime="text/csv",
            key="download_existing_matrix_btn",
        )

        # Backfill the gene symbol mapping for projects whose counts
        # matrix was built before this feature existed, so reopening an
        # older project doesn't require rebuilding the matrix just to
        # get gene symbols in Differential Expression.
        if not os.path.exists(pm.gene_symbol_map_path(project)):
            gene_symbol_map, gene_symbol_source = _get_gene_symbol_mapping(project, reference_dir, method)
            if gene_symbol_map:
                rm.save_gene_symbol_map_csv(gene_symbol_map, pm.gene_symbol_map_path(project))
                pm.save_gene_id_mapping_meta(project, {
                    "source": "auto_parse",
                    "detail": gene_symbol_source,
                })
                st.caption(f"🏷️ Gene symbol mapping ready ({gene_symbol_source}) -- will be used automatically in Differential Expression. For genes that fall back to their raw reference ID, the Differential Expression workspace can convert them via a proper annotation database lookup.")

    if os.path.exists(counts_matrix_path):
        st.markdown("---")
        st.success(
            f"🎉 Project `{project}` now has a combined gene counts "
            "matrix, ready for differential expression analysis."
        )

        if st.button("➡️ Proceed to Differential Expression", type="primary", key="align_proceed_de_btn"):
            # Same nav_request indirection used by the other "Proceed to X"
            # buttons in this app (see bulk_rnaseq_workspace.py and app.py's
            # module docstring for why a plain session key is used here
            # instead of directly setting st.session_state["assay_choice_radio"]).
            st.session_state["nav_request"] = "🌋 Differential Expression"
            st.rerun()

# Reuse the same workspace_key as the other Bulk RNA-Seq pipeline stages
# so the active project selection is shared automatically across pages.
WORKSPACE_KEY = "bulk_rnaseq"

# Accepted FASTA file extensions for custom reference uploads. ".fna"
# ("FASTA Nucleic Acid") is included because it's the standard extension
# NCBI uses for genome/RefSeq downloads — a very common source for
# non-model organism references.
FASTA_UPLOAD_TYPES = ["fa", "fasta", "fna", "fa.gz", "fasta.gz", "fna.gz"]
ANNOTATION_UPLOAD_TYPES = ["gtf", "gff", "gff3", "gtf.gz", "gff.gz", "gff3.gz"]

# Extension lists used by the server-side file browser (see
# file_browser.py) to filter its file listing to plausible candidates
# for each file type -- includes the leading "." (unlike
# FASTA_UPLOAD_TYPES/ANNOTATION_UPLOAD_TYPES above, which st.file_uploader
# expects WITHOUT a leading dot), since file_browser.render_server_file_browser
# matches these directly against os.path.splitext-style suffixes.
FASTA_BROWSE_EXTENSIONS = [f".{ext}" for ext in FASTA_UPLOAD_TYPES]
ANNOTATION_BROWSE_EXTENSIONS = [f".{ext}" for ext in ANNOTATION_UPLOAD_TYPES]


def _render_reference_file_input(purpose_label, key_prefix, file_extensions, dest_path,
                                  help_text=None):
    """
    Render a single reference-file input point that offers the user a
    choice between two ways to provide a file, and returns the
    resulting on-disk path once the file is available -- or None if
    nothing has been provided yet this run.

    This exists specifically to solve a real friction point discovered
    during actual deployment testing: when this app runs on a remote
    host (an HPC node, a shared lab server) that ALREADY has the needed
    reference file sitting on its own disk, forcing the user through
    st.file_uploader means needlessly transferring a potentially
    multi-GB file from that SAME machine, through a browser, back to
    itself -- pointlessly slow (or outright impractical for a
    multi-gigabyte genome FASTA) compared to simply pointing the app at
    the file's existing path.

    purpose_label: what this file is for, used in the radio button's
        section label (e.g. "genome FASTA", "GTF/GFF3 annotation").
    key_prefix: unique Streamlit widget/session_state key scoping for
        this specific input point (e.g. "custom_fasta", "custom_gtf")
        -- required so two calls to this function on the same page (one
        for the FASTA, one for the GTF) don't collide.
    file_extensions: list of extensions WITHOUT a leading dot (Streamlit
        file_uploader's own convention), e.g. FASTA_UPLOAD_TYPES -- used
        both for the upload widget's `type` argument and (converted to
        a leading-dot form) to filter the server-browse file listing.
    dest_path: where an UPLOADED file should be saved to on disk (server-
        browsed files are used directly from wherever they already are
        -- see file_browser.py's module docstring -- and are never
        copied here).
    help_text: optional help string shown under the section label.

    Returns the resulting file's path (str) if one is available after
    this render (either just-uploaded-and-saved, or an existing
    server-browsed selection), or None otherwise.
    """
    if help_text:
        st.caption(help_text)

    input_method = st.radio(
        f"How would you like to provide the {purpose_label}?",
        ["📤 Upload from your computer", "📂 Browse files already on this server"],
        key=f"{key_prefix}_input_method_radio",
        horizontal=True,
    )

    if input_method.startswith("📤"):
        uploaded_file = st.file_uploader(
            f"Upload {purpose_label}",
            type=file_extensions,
            key=f"{key_prefix}_upload",
        )
        if uploaded_file is not None:
            with open(dest_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            return dest_path
        # An upload from an EARLIER run may already be saved on disk at
        # dest_path even though nothing was just re-uploaded THIS run
        # (Streamlit's file_uploader doesn't persist a "previously
        # uploaded" state across reruns the way a plain value would) --
        # so fall back to whatever's already there, if anything.
        return dest_path if os.path.exists(dest_path) else None

    else:
        browse_extensions = [f".{ext}" for ext in file_extensions]
        selected_path = fb.render_server_file_browser(
            key_prefix=f"{key_prefix}_browse",
            file_extensions=browse_extensions,
            label=f"Browse for the {purpose_label} already on this server:",
        )
        # Server-browsed files are used directly from their existing
        # location -- never copied into dest_path -- since the entire
        # point of this option is to avoid needlessly duplicating a
        # potentially multi-GB file that's already sitting somewhere
        # perfectly usable.
        return selected_path


def render():
    st.title("🧮 RNA Alignment & Counts")
    st.markdown(
        "This workspace turns your trimmed reads into a gene-level "
        "counts table — how much each gene was expressed in each "
        "sample — ready for differential expression analysis. **No "
        "bioinformatics experience required** — follow the steps below "
        "in order."
    )
    st.markdown("---")

    # -----------------------------------------------------------------
    # Project selection — shared with the other Bulk RNA-Seq stages
    # -----------------------------------------------------------------
    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)

    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    st.markdown("---")

    # -----------------------------------------------------------------
    # Gate: require trimming to already be done
    # -----------------------------------------------------------------
    trimming_done = pm.has_completed_step(project, "trimming_complete")

    if not trimming_done:
        st.warning(
            "⚠️ This project doesn't have trimmed reads yet. Go to the "
            "**🧪 Trimming & Post-Trim QC** page first and complete "
            "Step 1 (trim your samples) before proceeding to alignment."
        )
        if st.button("⬅️ Go to Trimming & Post-Trim QC", key="align_gate_back_btn"):
            st.session_state["nav_request"] = "🧪 Trimming & Post-Trim QC"
            st.rerun()
        return

    st.success(f"✅ Using trimmed reads from project `{project}`.")

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 1: Choose an alignment/quantification method
    # -----------------------------------------------------------------
    st.header("Step 1: Choose an Alignment Method")

    st.markdown(
        "There are two common ways to measure gene expression from your "
        "reads. Both give you a valid gene counts table for differential "
        "expression — the right choice depends on what else you need "
        "from your data."
    )

    with st.expander("📊 Compare Salmon vs. STAR (click to see the full comparison)"):
        comparison_df = pd.DataFrame([
            {
                "Aspect": "What it needs",
                "Salmon": "Transcriptome FASTA only",
                "STAR": "Genome FASTA + gene annotation (GTF)",
            },
            {
                "Aspect": "Speed",
                "Salmon": "Fast (minutes per sample)",
                "STAR": "Slower (can be an hour+ per sample)",
            },
            {
                "Aspect": "Memory needed",
                "Salmon": "Low (a few GB)",
                "STAR": "High (~16–30GB depending on species)",
            },
            {
                "Aspect": "Detects novel splice junctions / isoforms",
                "Salmon": "❌ No — limited to known transcripts",
                "STAR": "✅ Yes",
            },
            {
                "Aspect": "Produces BAM files (for viewing in IGV, etc.)",
                "Salmon": "❌ No",
                "STAR": "✅ Yes",
            },
            {
                "Aspect": "Best for",
                "Salmon": "Standard differential expression, fast turnaround",
                "STAR": "Discovery-focused studies, novel isoform/fusion detection, publication pipelines needing BAM files",
            },
        ])
        st.dataframe(comparison_df, use_container_width=True, hide_index=True)
        st.caption(
            "**Bottom line:** if you're mainly asking \"which genes are "
            "up or down between my conditions\", Salmon is faster and "
            "just as statistically valid for that question. Choose STAR "
            "if you specifically need novel isoform/splice detection or "
            "BAM files for other tools."
        )

    saved_method = pm.get_alignment_method(project)
    method_options = ["Salmon (fast, transcriptome-based)", "STAR (splice-aware, genome-based)"]
    default_index = 0
    if saved_method == "star":
        default_index = 1

    method_choice = st.radio(
        "Which method would you like to use for this project?",
        method_options,
        index=default_index,
        key="alignment_method_radio",
    )
    method = "star" if method_choice.startswith("STAR") else "salmon"

    if saved_method and saved_method != method:
        st.warning(
            f"⚠️ This project previously used **{saved_method.upper()}**. "
            "Switching methods means Step 2 (reference setup) and Step 3 "
            "(quantification) will need to be run again for the new "
            "method — previous results from the other method are kept "
            "on disk but won't be used going forward."
        )

    pm.save_alignment_method(project, method)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 2: Reference setup (differs by method)
    # -----------------------------------------------------------------
    st.header("Step 2: Set Up Your Reference")

    if method == "salmon":
        with st.expander("ℹ️ What is a transcriptome reference? (click to learn more)"):
            st.markdown(
                "Salmon needs to know the sequence of every known "
                "transcript in your organism — this is called a "
                "**transcriptome reference** (sometimes labeled 'cDNA' "
                "on download sites like GENCODE or Ensembl). Salmon "
                "compares your trimmed reads against this reference to "
                "estimate how much each transcript/gene was expressed.\n\n"
                "This file is specific to a species (e.g. human vs. "
                "mouse) — using the wrong species' reference will "
                "produce meaningless results."
            )
    else:
        with st.expander("ℹ️ What does STAR need, and why is it more involved? (click to learn more)"):
            st.markdown(
                "STAR aligns your reads directly to the full **genome**, "
                "not just known transcripts — this is what lets it "
                "detect novel splice junctions. For this it needs two "
                "files:\n"
                "- A **genome FASTA** (the full DNA sequence of every "
                "chromosome)\n"
                "- A **GTF annotation file** (which tells STAR where "
                "genes and exons are located on the genome)\n\n"
                "STAR then builds a **genome index** from these files — "
                "a one-time setup step per species that can take a "
                "significant amount of memory (often 16–30GB depending "
                "on the organism) and take a while to complete, but only "
                "needs to be done once per reference."
            )

    reference_dir = pm.reference_dir(project)
    os.makedirs(reference_dir, exist_ok=True)

    saved_species, saved_is_custom = pm.get_reference_choice(project)

    ref_source_options = ["🧬 Use a pre-loaded model organism", "📁 Upload my own reference (custom species)"]
    default_source_index = 1 if saved_is_custom else 0
    ref_source = st.radio(
        "How would you like to provide your reference?",
        ref_source_options,
        index=default_source_index,
        key="ref_source_radio",
    )
    is_custom = ref_source.startswith("📁")

    # -------------------------------------------------------------
    # PRESET MODEL ORGANISM PATH
    # -------------------------------------------------------------
    if not is_custom:
        species_labels = {key: entry["label"] for key, entry in rm.REFERENCE_CATALOG.items()}
        species_keys = list(species_labels.keys())
        default_species_index = 0
        if saved_species in species_keys and not saved_is_custom:
            default_species_index = species_keys.index(saved_species)

        species_choice_label = st.selectbox(
            "Select your organism:",
            options=[species_labels[k] for k in species_keys],
            index=default_species_index,
            key="species_select",
        )
        species_key = species_keys[[species_labels[k] for k in species_keys].index(species_choice_label)]

        pm.save_reference_choice(project, species_key, is_custom=False)

        # Preset species reuse a SHARED, project-independent reference
        # location (see project_manager.py's shared_reference_dir() and
        # reference_manager.py's ensure_shared_resource()) -- the first
        # project to request a given species downloads/builds it once;
        # every other project (including ones in other Streamlit
        # sessions running at the same time) safely reuses the same
        # files rather than each downloading/building their own
        # separate, redundant copy. Custom uploads (the `else` branch
        # further below) are NOT shared -- see reference_dir there.
        reference_dir = pm.shared_reference_dir(species_key)
        os.makedirs(reference_dir, exist_ok=True)

        if method == "salmon":
            shared_cdna_dir = pm.shared_cdna_fasta_dir(species_key)
            already_downloaded = rm._resource_is_ready(shared_cdna_dir)

            if already_downloaded:
                st.success(
                    f"✅ Reference already available for **{species_labels[species_key]}** "
                    "(shared across every project using this species)."
                )
            else:
                st.info(
                    f"ℹ️ No project has downloaded a reference for "
                    f"**{species_labels[species_key]}** yet. The first "
                    "download will be saved in a shared location and "
                    "reused automatically by this and any other project "
                    "using this species -- no need to re-download it "
                    "again later."
                )

            dl_label = "⬇️ Download Reference" if not already_downloaded else None
            force_dl = False
            if already_downloaded:
                with st.expander("🔁 Force a fresh re-download (advanced)"):
                    st.caption(
                        "⚠️ This reference is shared across every project "
                        "using this species. Re-downloading replaces it "
                        "for everyone -- only do this if you have a "
                        "specific reason to believe the current copy is "
                        "corrupted or out of date, and ideally not while "
                        "another project may be actively using it."
                    )
                    force_dl = st.button("🔄 Re-download Reference for All Projects", key="download_cdna_force_btn")

            if (dl_label and st.button(dl_label, key="download_cdna_btn")) or force_dl:
                status_placeholder = st.empty()
                progress_bar = st.progress(0, text="Starting download...")

                def _update_progress(downloaded, total):
                    if total > 0:
                        pct = min(downloaded / total, 1.0)
                        progress_bar.progress(pct, text=f"Downloading... {downloaded // (1024*1024)}MB / {total // (1024*1024)}MB")

                def _wait_callback(elapsed):
                    status_placeholder.info(
                        f"⏳ Another project is currently preparing this "
                        f"reference. Waiting for it to finish... "
                        f"({int(elapsed)}s elapsed)"
                    )

                def _build_cdna_impl(temp_dir):
                    success, fasta_path, message = rm.get_transcriptome_fasta_for_salmon(species_key, temp_dir, _update_progress)
                    if success and not rm.validate_fasta_file(fasta_path):
                        return False, "The downloaded file doesn't look like a valid FASTA file."
                    return success, message

                success, message, built = rm.ensure_shared_resource(
                    shared_cdna_dir, _build_cdna_impl, wait_message_callback=_wait_callback, force=force_dl,
                )
                progress_bar.progress(1.0, text="Done.")
                status_placeholder.empty()

                if not success:
                    st.error(f"⚠️ {message}")
                else:
                    st.success(f"✅ {message}")
                    pm.mark_step_complete(project, "reference_ready")

        else:  # STAR
            shared_genome_dir = pm.shared_genome_dir(species_key)
            already_downloaded = rm._resource_is_ready(shared_genome_dir)

            if already_downloaded:
                st.success(
                    f"✅ Reference already available for **{species_labels[species_key]}** "
                    "(shared across every project using this species)."
                )
            else:
                st.info(
                    f"ℹ️ No project has downloaded a reference for "
                    f"**{species_labels[species_key]}** yet. The first "
                    "download will be saved in a shared location and "
                    "reused automatically by this and any other project "
                    "using this species -- no need to re-download it "
                    "again later."
                )

            st.caption(
                "⏳ This download includes the full genome sequence and may "
                "take a while depending on species size and your internet "
                "connection."
            )

            dl_label = "⬇️ Download Genome + Annotation" if not already_downloaded else None
            force_dl = False
            if already_downloaded:
                with st.expander("🔁 Force a fresh re-download (advanced)"):
                    st.caption(
                        "⚠️ This reference is shared across every project "
                        "using this species. Re-downloading replaces it "
                        "for everyone -- only do this if you have a "
                        "specific reason to believe the current copy is "
                        "corrupted or out of date, and ideally not while "
                        "another project may be actively using it."
                    )
                    force_dl = st.button("🔄 Re-download Reference for All Projects", key="download_genome_gtf_force_btn")

            if (dl_label and st.button(dl_label, key="download_genome_gtf_btn")) or force_dl:
                status_placeholder = st.empty()
                progress_bar = st.progress(0, text="Starting download...")

                def _update_progress(downloaded, total):
                    if total > 0:
                        pct = min(downloaded / total, 1.0)
                        progress_bar.progress(pct, text=f"Downloading... {downloaded // (1024*1024)}MB / {total // (1024*1024)}MB")

                def _wait_callback(elapsed):
                    status_placeholder.info(
                        f"⏳ Another project is currently preparing this "
                        f"reference. Waiting for it to finish... "
                        f"({int(elapsed)}s elapsed)"
                    )

                def _build_genome_impl(temp_dir):
                    success, paths, message = rm.download_genome_and_gtf(species_key, temp_dir, _update_progress)
                    if not success:
                        return False, message
                    genome_fa, gtf_file = paths
                    if not rm.validate_fasta_file(genome_fa):
                        return False, "The downloaded genome file doesn't look like a valid FASTA file."
                    if not rm.validate_annotation_file(gtf_file):
                        return False, "The downloaded annotation file doesn't look like a valid GTF file."
                    return True, message

                success, message, built = rm.ensure_shared_resource(
                    shared_genome_dir, _build_genome_impl, wait_message_callback=_wait_callback, force=force_dl,
                )
                progress_bar.progress(1.0, text="Done.")
                status_placeholder.empty()

                if not success:
                    st.error(f"⚠️ {message}")
                else:
                    st.success(f"✅ {message}")
                    pm.mark_step_complete(project, "reference_ready")

    # -------------------------------------------------------------
    # CUSTOM (NON-MODEL) SPECIES PATH
    # -------------------------------------------------------------
    else:
        pm.save_reference_choice(project, "custom", is_custom=True)

        # Custom uploads are always scoped to THIS project's own
        # private reference_dir -- never shared -- since two different
        # projects' custom uploads have no guarantee of actually being
        # the same organism/assembly even if they happen to look
        # similar, so sharing them would risk silently mixing up
        # unrelated references across projects. Note this applies
        # regardless of whether the file arrived via upload or via the
        # server-side browser below -- ONLY uploaded files get COPIED
        # into this directory; server-browsed files are referenced
        # directly from wherever they already live (see
        # _render_reference_file_input's docstring) and are never
        # duplicated into reference_dir at all.
        reference_dir = pm.reference_dir(project)
        os.makedirs(reference_dir, exist_ok=True)

        st.markdown(
            "Provide your reference files below. You'll need a "
            "**genome or transcriptome FASTA** file, plus a **GTF or "
            "GFF3 annotation** file. For each one, you can either "
            "upload it from your own computer, or -- if this app is "
            "running on a shared server/HPC that already has the file "
            "sitting on its disk -- browse for it directly there "
            "instead, without needing to transfer it through your "
            "browser at all (much faster for large genome files)."
        )

        with st.expander("ℹ️ Not sure which files you need? (click to learn more)"):
            st.markdown(
                "- **For Salmon:** ideally a transcriptome/cDNA FASTA "
                "(one sequence per transcript). If you only have a "
                "genome FASTA + annotation, that's fine too — we'll "
                "automatically extract transcript sequences for you "
                "using a tool called `gffread` (for organisms with "
                "introns) or direct gene-coordinate extraction (for "
                "organisms without introns, like bacteria).\n"
                "- **For STAR:** a genome FASTA (full chromosome "
                "sequences) + a GTF/GFF3 annotation file are both "
                "required.\n\n"
                "FASTA files usually have extensions like `.fa`, "
                "`.fasta`, or `.fna` (`.fna` is common for genomes "
                "downloaded from NCBI). Annotation files usually end in "
                "`.gtf`, `.gff`, `.gff3`, or their gzipped versions."
            )

        custom_fasta_dest = os.path.join(reference_dir, "custom_input.fa")
        custom_gtf_dest = os.path.join(reference_dir, "custom_input.gtf")

        fasta_label = "genome FASTA" if method == "star" else "FASTA file (transcriptome preferred; genome FASTA also accepted if paired with an annotation below)"
        st.markdown(f"**{'Genome FASTA' if method == 'star' else 'FASTA file'}** ({'required for STAR' if method == 'star' else 'transcriptome preferred'}):")
        fasta_path = _render_reference_file_input(
            purpose_label=fasta_label,
            key_prefix="custom_fasta",
            file_extensions=FASTA_UPLOAD_TYPES,
            dest_path=custom_fasta_dest,
        )

        st.markdown("---")
        st.markdown("**Annotation file** (GTF or GFF3):")
        annotation_help = (
            "Required for STAR. For Salmon, only required if you provided a "
            "genome FASTA rather than a transcriptome FASTA."
        )
        gtf_path_input = _render_reference_file_input(
            purpose_label="GTF/GFF3 annotation file",
            key_prefix="custom_gtf",
            file_extensions=ANNOTATION_UPLOAD_TYPES,
            dest_path=custom_gtf_dest,
            help_text=annotation_help,
        )

        # Let the user tell us directly whether their organism has
        # introns, rather than guessing — this determines whether custom
        # Salmon extraction (when only a genome + annotation is given)
        # uses gffread or direct gene-coordinate extraction.
        is_no_intron_organism = False
        if method == "salmon":
            is_no_intron_organism = st.checkbox(
                "This organism has no introns (e.g. bacteria, archaea, viruses)",
                value=False,
                key="no_intron_checkbox",
                help=(
                    "Check this if your species doesn't splice its genes "
                    "(most prokaryotes). This changes how transcript "
                    "sequences are extracted if you upload a genome + "
                    "annotation instead of a ready-made transcriptome — "
                    "gffread often fails on no-intron GTF/GFF files "
                    "because they lack a separate transcript/mRNA feature "
                    "line for it to anchor to."
                ),
            )

        st.markdown("---")
        if fasta_path and st.button("💾 Confirm Reference Files", key="save_custom_ref_btn"):
            # fasta_path/gtf_path_input may point to a file we just
            # copied into reference_dir (upload path) OR to a file
            # that's already sitting somewhere ELSE entirely on this
            # server (browse path) -- either way, by this point they're
            # simply "the path to use", so the validation/completion
            # logic below doesn't need to distinguish between the two
            # any further. We DO need to remember which path was
            # actually chosen, though, since every downstream step
            # (Salmon/STAR index building, tx2gene/gene-symbol mapping)
            # currently assumes the fixed conventional filenames
            # custom_input.fa/custom_input.gtf inside reference_dir --
            # so a server-browsed file's path is persisted explicitly
            # rather than silently relying on those fixed names.
            if not rm.validate_fasta_file(fasta_path):
                st.error("⚠️ That FASTA file doesn't look valid. Please check it and try again.")
            else:
                st.success(f"✅ FASTA file confirmed ({os.path.basename(fasta_path)}).")
                st.session_state["_custom_fasta_actual_path"] = fasta_path

                gtf_saved = False
                if gtf_path_input:
                    if not rm.validate_annotation_file(gtf_path_input):
                        st.error("⚠️ That annotation file doesn't look like a valid GTF/GFF file.")
                    else:
                        st.success(f"✅ Annotation file confirmed ({os.path.basename(gtf_path_input)}).")
                        st.session_state["_custom_gtf_actual_path"] = gtf_path_input
                        gtf_saved = True

                if method == "star" and not gtf_saved:
                    st.error("⚠️ STAR requires both a genome FASTA and a GTF/GFF3 annotation file. Please provide both.")
                else:
                    pm.mark_step_complete(project, "reference_ready")
                    st.session_state["_custom_ref_saved"] = True
                    st.session_state["_custom_no_intron"] = is_no_intron_organism

        # Resolve the ACTUAL fasta/gtf paths to use for extraction/
        # indexing below -- prefer whatever was just confirmed this run
        # (st.session_state["_custom_fasta_actual_path"], set above),
        # falling back to the conventional custom_input.fa/gtf location
        # for backward compatibility with projects set up before this
        # server-browse option existed (which always used that fixed
        # path via a plain upload).
        resolved_custom_fasta = st.session_state.get("_custom_fasta_actual_path", custom_fasta_dest)
        resolved_custom_gtf = st.session_state.get("_custom_gtf_actual_path", custom_gtf_dest)

        # If this looks like a genome FASTA (not already a transcriptome)
        # and we're using Salmon, offer to extract transcripts.
        if method == "salmon" and st.session_state.get("_custom_ref_saved") and os.path.exists(resolved_custom_gtf):
            st.markdown("---")
            st.info(
                "Since you provided a genome + annotation, Salmon needs "
                "transcript-level sequences extracted from them first."
            )
            if st.button("🧬 Extract Transcript Sequences", key="extract_btn"):
                extracted_path = os.path.join(reference_dir, "custom_extracted_transcripts.fa")
                use_no_intron_method = st.session_state.get("_custom_no_intron", False)

                with st.spinner("Extracting transcript sequences..."):
                    if use_no_intron_method:
                        success, message = rm.extract_gene_level_transcripts(
                            resolved_custom_fasta, resolved_custom_gtf, extracted_path
                        )
                    else:
                        success, message = rm.extract_transcripts_with_gffread(
                            resolved_custom_fasta, resolved_custom_gtf, extracted_path
                        )

                if success:
                    st.success(f"✅ {message}")
                else:
                    st.error(f"⚠️ {message}")
                    if not use_no_intron_method:
                        st.info(
                            "💡 If your organism has no introns (bacteria, "
                            "archaea, viruses), check the box above and "
                            "try extraction again — gffread often fails on "
                            "no-intron annotation files because they lack "
                            "a separate transcript/mRNA feature line."
                        )

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 3: Run Quantification
    # -----------------------------------------------------------------
    st.header("Step 3: Run Quantification")

    reference_ready = pm.has_completed_step(project, "reference_ready")
    if not reference_ready:
        st.info("Complete Step 2 (set up your reference) to unlock this step.")
        return

    samplesheet_path = pm.samplesheet_path(project)
    if not os.path.exists(samplesheet_path):
        st.warning("⚠️ No matched sample list found for this project. Please check Step 3 of the Bulk RNA-Seq Pipeline page.")
        return

    samplesheet_df = pd.read_csv(samplesheet_path)
    trimmed_dir = pm.trimmed_fastq_dir(project)

    # --- Auto-detect paired-end vs. single-end, and show the user ---
    st.subheader("📋 Detected Read Layout")
    st.markdown(
        "Before running anything, here's what we detected for each "
        "sample based on your trimmed files. **Please double-check this "
        "matches what you expect** — if a sample shows up as \"missing\" "
        "or the wrong type, something may have gone wrong in trimming."
    )

    manifest = qm.build_sample_manifest(samplesheet_df, trimmed_dir)
    manifest_display = []
    read_type_labels = {
        "paired": "✅ Paired-end (R1 + R2)",
        "single": "✅ Single-end",
        "missing_mate": "⚠️ R1 found, but R2 missing",
        "missing": "❌ No trimmed files found",
    }
    for entry in manifest:
        manifest_display.append({
            "Sample": entry["sample"],
            "Detected Layout": read_type_labels.get(entry["read_type"], entry["read_type"]),
        })
    st.dataframe(pd.DataFrame(manifest_display), use_container_width=True, hide_index=True)

    problem_samples = [e["sample"] for e in manifest if e["read_type"] in ("missing", "missing_mate")]
    if problem_samples:
        st.error(
            f"⚠️ These samples have incomplete trimmed data and will be "
            f"skipped during quantification: **{', '.join(problem_samples)}**. "
            "Revisit the Trimming & Post-Trim QC page if this is unexpected."
        )

    usable_manifest = [e for e in manifest if e["read_type"] in ("paired", "single")]
    if not usable_manifest:
        st.warning("No samples with usable trimmed reads were found — nothing to quantify yet.")
        return

    st.markdown("---")

    if method == "salmon":
        _render_salmon_quantification(project, reference_dir, trimmed_dir, usable_manifest)
    else:
        _render_star_quantification(project, reference_dir, trimmed_dir, usable_manifest)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 4: Combine Into a Gene Counts Matrix
    # -----------------------------------------------------------------
    _render_counts_matrix_step(project, method, samplesheet_df, reference_dir)


def _render_salmon_quantification(project, reference_dir, trimmed_dir, manifest):
    """Step 3 UI for the Salmon path: index check/build, then quantify."""
    if not qm.salmon_tool_available():
        st.error(
            "⚠️ Salmon was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
        return

    cdna_candidates = [
        os.path.join(reference_dir, "custom_extracted_transcripts.fa"),
    ]
    species_key, is_custom = pm.get_reference_choice(project)
    if not is_custom and species_key:
        # Preset species' shared cDNA FASTA lives in a "cdna"
        # subdirectory of reference_dir -- see _get_tx2gene_mapping's
        # comment above for why.
        cdna_candidates.insert(0, os.path.join(reference_dir, "cdna", f"{species_key}.cdna.fa"))
    else:
        # A custom reference's FASTA may be either the conventional
        # copied-in path (custom_input.fa, for an uploaded file) or an
        # arbitrary server-side path (for a browsed file, remembered in
        # session_state by _render_reference_file_input's caller) --
        # prefer whichever was actually confirmed, falling back to the
        # conventional path for backward compatibility with projects
        # set up before the server-browse option existed.
        cdna_candidates.insert(0, st.session_state.get(
            "_custom_fasta_actual_path", os.path.join(reference_dir, "custom_input.fa")
        ))

    transcriptome_fasta = next((p for p in cdna_candidates if os.path.exists(p)), None)
    if not transcriptome_fasta:
        st.error("⚠️ Could not find a transcriptome FASTA for this project. Please complete Step 2 first.")
        return

    n_transcripts = rm.count_transcripts_in_fasta(transcriptome_fasta)
    st.caption(f"Using transcriptome reference: `{os.path.basename(transcriptome_fasta)}` ({n_transcripts:,} sequences)")

    # Auto-detected thread default, same pattern as FastQC/fastp — used
    # for both the index-build step and the quantification step below,
    # since both accept Salmon's "-p" flag via quantification_manager.py.
    #
    # IMPORTANT: detected_cores (the machine's REAL core count) and
    # recommended_threads (a conservative SUGGESTED starting value,
    # deliberately capped low by get_recommended_thread_count) serve two
    # different purposes here. Every slider below uses detected_cores
    # for max_value (so a user on a large machine -- e.g. an HPC node --
    # can actually select a high thread count) and recommended_threads
    # only for the slider's initial `value` (a sensible, un-intimidating
    # starting point). A previous version of this file incorrectly used
    # max(8, recommended_threads) as max_value, which -- since
    # recommended_threads is itself always <= 8 -- silently capped every
    # slider at 8 regardless of how many cores were actually detected
    # (confirmed via real testing on a 32-core machine, where every
    # slider was stuck at a max of 8 despite 32 cores being available).
    detected_cores, recommended_threads = pm.get_recommended_thread_count()

    # --- Index step ---
    st.subheader("🔧 Salmon Index")
    with st.expander("ℹ️ When would I need to re-index? (click to learn more)"):
        st.markdown(
            "The Salmon index only needs to be built **once** per "
            "reference and is reused for every sample after that. You "
            "would want to **re-index** if:\n"
            "- You changed your reference (uploaded a different FASTA, "
            "or switched species/organism)\n"
            "- You re-extracted transcripts (e.g. re-ran gffread with "
            "different settings)\n"
            "- The index seems corrupted or quantification is failing "
            "with unexplained errors\n\n"
            "You do **not** need to re-index just because you added new "
            "samples using the *same* reference — the existing index "
            "will work fine for those."
        )

    # Preset species reuse a SHARED, project-independent Salmon index
    # (built once per species, reused by every project) -- custom
    # uploads keep their existing per-project index, since there's no
    # guarantee two different projects' custom references are actually
    # the same thing even if named similarly. See reference_manager.py's
    # ensure_shared_resource() for how concurrent index-build requests
    # across different projects/sessions are made safe.
    index_is_shared = not is_custom and bool(species_key)
    index_dir = pm.shared_salmon_index_dir(species_key) if index_is_shared else pm.salmon_index_dir(project)
    # qm.salmon_index_exists just checks whether a valid-looking Salmon
    # index is present at a given path -- this works identically
    # whether that path happens to be the shared or per-project
    # location, so the same check is used for both.
    index_ready = qm.salmon_index_exists(index_dir)

    index_threads = st.slider(
        "Threads for indexing:",
        min_value=1, max_value=detected_cores, value=recommended_threads,
        help=(
            f"Detected {detected_cores} CPU core(s) on this machine, so "
            f"{recommended_threads} is suggested as a starting point. "
            "Indexing is a one-time step per reference, so it's usually "
            "fine to use most (or all) of your available cores here -- "
            "you can raise this slider up to your machine's full core "
            "count if you'd like to speed this up further."
        ),
        key="salmon_index_threads_slider",
    )

    if index_ready:
        st.success(
            "✅ Salmon index already built "
            + ("(shared across every project using this species)." if index_is_shared else "for this project's reference.")
        )

    build_clicked = False
    force_index = False
    if not index_ready:
        build_clicked = st.button("🔧 Build Salmon Index", key="build_salmon_index_btn")
    elif index_is_shared:
        with st.expander("🔁 Force a fresh re-index (advanced)"):
            st.caption(
                "⚠️ This index is shared across every project using this "
                "species. Rebuilding replaces it for everyone -- only do "
                "this if you have a specific reason to believe the "
                "current index is corrupted, and ideally not while "
                "another project may be actively quantifying against it."
            )
            force_index = st.button("🔄 Re-build Index for All Projects", key="build_salmon_index_force_btn")
    else:
        build_clicked = st.button("🔄 Re-build Index", key="build_salmon_index_btn")

    if build_clicked or force_index:
        def _build_salmon_index_impl(temp_dir):
            return qm.build_salmon_index(transcriptome_fasta, temp_dir, threads=index_threads)

        if index_is_shared:
            status_placeholder = st.empty()

            def _wait_callback(elapsed):
                status_placeholder.info(
                    f"⏳ Another project is currently building this index. "
                    f"Waiting for it to finish... ({int(elapsed)}s elapsed)"
                )

            with st.spinner("Building Salmon index... this may take a few minutes."):
                success, log, built = rm.ensure_shared_resource(
                    index_dir, _build_salmon_index_impl, wait_message_callback=_wait_callback, force=force_index,
                )
            status_placeholder.empty()
        else:
            with st.spinner("Building Salmon index... this may take a few minutes."):
                success, log = qm.build_salmon_index(transcriptome_fasta, index_dir, threads=index_threads)

        if success:
            st.success("✅ Salmon index built successfully.")
            index_ready = True
        else:
            st.error("Salmon indexing failed. Details below:")
            st.code(log)
            return

    if not index_ready:
        st.info("Build the Salmon index above before running quantification.")
        return

    st.markdown("---")

    # --- Quantification step ---
    st.subheader("🚀 Run Quantification")
    quant_dir = pm.salmon_quant_dir(project)
    quant_done = pm.has_completed_step(project, "quantification_complete")

    quant_threads = st.slider(
        "Threads per sample (quantification):",
        min_value=1, max_value=detected_cores, value=recommended_threads,
        help=(
            f"Detected {detected_cores} CPU core(s) on this machine, so "
            f"{recommended_threads} is suggested as a starting point. "
            "Since samples are quantified one at a time (not "
            "simultaneously), this controls how many threads Salmon "
            "uses *within* each sample's quantification run -- you can "
            "raise this up to your machine's full core count."
        ),
        key="salmon_quant_threads_slider",
    )

    quant_label = "🔄 Re-run Quantification" if quant_done else "🚀 Quantify All Samples"
    if quant_done and not st.session_state.get("_salmon_quant_clicked"):
        st.success("✅ Quantification has already been run for this project.")

    if st.button(quant_label, key="run_salmon_quant_btn"):
        st.session_state["_salmon_quant_clicked"] = True
        progress_bar = st.progress(0, text="Starting quantification...")
        results = []
        n_samples = len(manifest)

        for i, entry in enumerate(manifest):
            progress_bar.progress(i / n_samples, text=f"Quantifying {i + 1}/{n_samples}: {entry['sample']}...")
            success, log, sample_out_dir = qm.run_salmon_quant(entry, index_dir, quant_dir, threads=quant_threads)
            mapping_rate = qm.parse_salmon_mapping_rate(sample_out_dir) if success else None
            results.append({
                "Sample": entry["sample"],
                "Status": "✅ Success" if success else "❌ Failed",
                "Mapping Rate": f"{mapping_rate}%" if mapping_rate is not None else "—",
            })
            if not success:
                with st.expander(f"Error details for {entry['sample']}"):
                    st.code(log)

        progress_bar.progress(1.0, text="Quantification complete.")
        results_df = pd.DataFrame(results)
        st.dataframe(results_df, use_container_width=True, hide_index=True)

        if (results_df["Status"] == "✅ Success").any():
            pm.mark_step_complete(project, "quantification_complete")
            st.success("✅ Quantification complete for at least one sample.")
            st.caption(
                "**Mapping Rate** shows the % of reads that matched a "
                "known transcript. Rates above ~70-80% are typically "
                "considered good for well-annotated organisms; lower "
                "rates can indicate contamination, wrong species "
                "reference, or degraded samples."
            )


def _render_star_quantification(project, reference_dir, trimmed_dir, manifest):
    """Step 3 UI for the STAR path: genome index check/build, then align + count."""
    if not qm.star_tool_available():
        st.error(
            "⚠️ STAR was not found on this system. It needs to be "
            "installed in your environment (included in the project's "
            "Dockerfile) before this step can run."
        )
        return

    species_key, is_custom = pm.get_reference_choice(project)
    if not is_custom and species_key:
        # Preset species' shared genome + GTF live in a "genome"
        # subdirectory of reference_dir -- see _get_tx2gene_mapping's
        # comment above for why.
        genome_fasta = os.path.join(reference_dir, "genome", f"{species_key}.genome.fa")
        gtf_path = os.path.join(reference_dir, "genome", f"{species_key}.annotation.gtf")
    else:
        # A custom reference's FASTA/GTF may be either the conventional
        # copied-in path (custom_input.fa/gtf, for an uploaded file) or
        # an arbitrary server-side path (for browsed files) -- prefer
        # whichever was actually confirmed in Step 2, falling back to
        # the conventional path for backward compatibility with
        # projects set up before the server-browse option existed.
        genome_fasta = st.session_state.get(
            "_custom_fasta_actual_path", os.path.join(reference_dir, "custom_input.fa")
        )
        gtf_path = st.session_state.get(
            "_custom_gtf_actual_path", os.path.join(reference_dir, "custom_input.gtf")
        )

    if not (os.path.exists(genome_fasta) and os.path.exists(gtf_path)):
        st.error("⚠️ Could not find a genome + annotation for this project. Please complete Step 2 first.")
        return

    st.caption(f"Using genome reference: `{os.path.basename(genome_fasta)}` + `{os.path.basename(gtf_path)}`")

    # Auto-detected thread default, same pattern as Salmon above — STAR
    # is generally more CPU- and memory-hungry than Salmon, so having an
    # accurate, machine-specific default matters even more here. See
    # the matching comment in _render_salmon_quantification above for
    # why detected_cores (not recommended_threads) is used as every
    # slider's max_value below.
    detected_cores, recommended_threads = pm.get_recommended_thread_count()

    # --- Auto-detect read length from an actual trimmed sample, for
    # STAR's --sjdbOverhang (ideally read_length - 1 per STAR's manual).
    # Previously this app always used a fixed default of 100 regardless
    # of the actual data -- fine for typical 100-150bp reads, but wrong
    # for anything else (e.g. our own airway test data's 63bp reads,
    # where the correct value is 62, not 99). Detected here (once, from
    # the first usable sample) rather than left as a fixed guess, so
    # the index build below is correctly sized for THIS project's
    # actual read length automatically.
    detected_read_length = None
    for entry in manifest:
        detected_read_length = qm.detect_fastq_read_length(entry["r1"], is_gzipped=True)
        if detected_read_length:
            break

    if detected_read_length:
        sjdb_overhang = detected_read_length - 1
        st.caption(
            f"📏 Detected read length: **{detected_read_length}bp** "
            f"(from `{os.path.basename(manifest[0]['r1'])}`) -- using "
            f"`--sjdbOverhang {sjdb_overhang}` for the genome index below."
        )
    else:
        sjdb_overhang = 100
        st.caption(
            "📏 Could not automatically detect read length from your "
            "trimmed reads -- falling back to the default "
            "`--sjdbOverhang 100` (appropriate for typical 100-150bp "
            "reads; if your reads are a very different length, this "
            "may reduce splice-junction detection accuracy slightly)."
        )

    # --- Index step ---
    st.subheader("🔧 STAR Genome Index")
    with st.expander("ℹ️ When would I need to re-index? (click to learn more)"):
        st.markdown(
            "The STAR genome index only needs to be built **once** per "
            "reference and is reused for every sample after that. You "
            "would want to **re-index** if:\n"
            "- You changed your reference (uploaded a different genome "
            "FASTA/GTF, or switched species/organism)\n"
            "- Your read length changed substantially from previous "
            "samples (the index is optimized using a value called "
            "`sjdbOverhang`, ideally read length minus 1 — a big change "
            "in read length can benefit from a fresh index; this app "
            "now detects your read length automatically each time, so "
            "you'll see a note above if it changed)\n"
            "- The index seems corrupted or alignment is failing with "
            "unexplained errors\n\n"
            "You do **not** need to re-index just because you added new "
            "samples using the *same* reference and similar read "
            "lengths — the existing index will work fine for those.\n\n"
            "⚠️ Building a STAR index can require significant memory "
            "(16–30GB depending on species) and take a while — this is "
            "normal and only happens once per reference. Note that this "
            "memory requirement is driven by genome size, not thread "
            "count -- raising the thread count below speeds up the "
            "CPU-bound portion of indexing, but won't reduce (or "
            "increase) how much RAM is needed.\n\n"
            "ℹ️ This app also automatically sizes STAR's "
            "`--genomeSAindexNbases` parameter based on your genome's "
            "actual size -- important for smaller genomes (a single "
            "chromosome, a bacterial genome, or another small custom "
            "reference), where STAR's large-genome default would "
            "otherwise produce an oversized or invalid index."
        )

    # Preset species reuse a SHARED, project-independent STAR index --
    # same rationale as the Salmon index above.
    index_is_shared = not is_custom and bool(species_key)
    index_dir = pm.shared_star_index_dir(species_key) if index_is_shared else pm.star_index_dir(project)
    index_ready = qm.star_index_exists(index_dir)

    index_threads = st.slider(
        "Threads for indexing:",
        min_value=1, max_value=detected_cores, value=recommended_threads,
        help=(
            f"Detected {detected_cores} CPU core(s) on this machine, so "
            f"{recommended_threads} is suggested as a starting point. "
            "STAR's indexing step needs a fairly fixed amount of memory "
            "regardless of thread count (16-30GB depending on species) "
            "-- but more threads DOES meaningfully speed up the "
            "CPU-bound portion of indexing (STAR's suffix array "
            "construction parallelizes well), so raising this toward "
            "your machine's full core count is worthwhile if you have "
            "the RAM to support it."
        ),
        key="star_index_threads_slider",
    )

    if index_ready:
        st.success(
            "✅ STAR genome index already built "
            + ("(shared across every project using this species)." if index_is_shared else "for this project's reference.")
        )

    build_clicked = False
    force_index = False
    if not index_ready:
        build_clicked = st.button("🔧 Build STAR Index", key="build_star_index_btn")
    elif index_is_shared:
        with st.expander("🔁 Force a fresh re-index (advanced)"):
            st.caption(
                "⚠️ This index is shared across every project using this "
                "species. Rebuilding replaces it for everyone -- only do "
                "this if you have a specific reason to believe the "
                "current index is corrupted, and ideally not while "
                "another project may be actively aligning against it."
            )
            force_index = st.button("🔄 Re-build Index for All Projects", key="build_star_index_force_btn")
    else:
        build_clicked = st.button("🔄 Re-build Index", key="build_star_index_btn")

    if build_clicked or force_index:
        def _build_star_index_impl(temp_dir):
            # genome_sa_index_nbases is left as None (the default) so
            # build_star_index auto-computes it from the genome's
            # actual length -- see quantification_manager.py's
            # compute_genome_sa_index_nbases() docstring for why this
            # matters, especially for smaller genomes.
            return qm.build_star_index(
                genome_fasta, gtf_path, temp_dir, threads=index_threads,
                sjdb_overhang=sjdb_overhang,
            )

        if index_is_shared:
            status_placeholder = st.empty()

            def _wait_callback(elapsed):
                status_placeholder.info(
                    f"⏳ Another project is currently building this index. "
                    f"Waiting for it to finish... ({int(elapsed)}s elapsed)"
                )

            with st.spinner("Building STAR genome index... this can take a while, especially for larger genomes."):
                success, log, built = rm.ensure_shared_resource(
                    index_dir, _build_star_index_impl, wait_message_callback=_wait_callback, force=force_index,
                )
            status_placeholder.empty()
        else:
            with st.spinner("Building STAR genome index... this can take a while, especially for larger genomes."):
                success, log = qm.build_star_index(
                    genome_fasta, gtf_path, index_dir, threads=index_threads,
                    sjdb_overhang=sjdb_overhang,
                )

        if success:
            st.success("✅ STAR genome index built successfully.")
            index_ready = True
        else:
            st.error("STAR indexing failed. Details below:")
            st.code(log)
            return

    if not index_ready:
        st.info("Build the STAR genome index above before running alignment.")
        return

    st.markdown("---")

    # --- Alignment + counting step ---
    st.subheader("🚀 Run Alignment & Gene Counting")
    align_dir = pm.star_align_dir(project)
    quant_done = pm.has_completed_step(project, "quantification_complete")

    align_threads = st.slider(
        "Threads per sample (alignment):",
        min_value=1, max_value=detected_cores, value=recommended_threads,
        help=(
            f"Detected {detected_cores} CPU core(s) on this machine, so "
            f"{recommended_threads} is suggested as a starting point. "
            "Since samples are aligned one at a time (not "
            "simultaneously), this controls how many threads STAR uses "
            "*within* each sample's alignment run -- you can raise this "
            "up to your machine's full core count."
        ),
        key="star_align_threads_slider",
    )

    # --- Advanced STAR alignment options ---
    # Exposed as a small set of well-established, commonly-used toggles
    # rather than individual low-level flags -- a non-expert user has no
    # meaningful basis to decide on any single underlying flag in
    # isolation, so each option here bundles a validated, documented
    # STAR feature/preset instead. See quantification_manager.py's
    # run_star_align docstring for exactly what each one does.
    with st.expander("⚙️ Advanced STAR options (optional)"):
        use_two_pass = st.checkbox(
            "Two-pass mode (better novel splice junction detection)",
            value=False,
            key="star_two_pass_checkbox",
            help=(
                "STAR aligns your reads once to discover splice "
                "junctions, then re-aligns using those junctions as "
                "extra annotation -- meaningfully improves detection of "
                "splice junctions not already in your GTF. Most useful "
                "for less-annotated/non-model organisms; roughly "
                "doubles alignment time per sample. Well-annotated "
                "genomes (human/mouse) benefit less since most real "
                "junctions are usually already known."
            ),
        )
        use_encode_options = st.checkbox(
            "Use ENCODE-recommended filtering options",
            value=False,
            key="star_encode_options_checkbox",
            help=(
                "Applies the standard set of alignment filtering "
                "settings used in the ENCODE project's own RNA-seq "
                "pipelines (multimapping tolerance, splice junction "
                "overhang minimums, mismatch tolerance, intron length "
                "bounds). A well-established, commonly-used preset "
                "bundle rather than a single tunable setting."
            ),
        )
        add_strand_field = st.checkbox(
            "Add strand field to BAM output",
            value=False,
            key="star_strand_field_checkbox",
            help=(
                "Annotates each alignment's strand in the output BAM "
                "file based on intron motif. Only needed if you plan to "
                "feed this project's BAM files into a downstream tool "
                "that expects this strand tag (e.g. Cufflinks-family "
                "tools) -- has no effect on the gene counts this app "
                "itself produces, so most users can safely leave this "
                "unchecked."
            ),
        )

    quant_label = "🔄 Re-run Alignment" if quant_done else "🚀 Align & Count All Samples"
    if quant_done and not st.session_state.get("_star_align_clicked"):
        st.success("✅ Alignment & counting has already been run for this project.")

    if st.button(quant_label, key="run_star_align_btn"):
        st.session_state["_star_align_clicked"] = True
        progress_bar = st.progress(0, text="Starting alignment...")
        results = []
        n_samples = len(manifest)

        for i, entry in enumerate(manifest):
            progress_bar.progress(i / n_samples, text=f"Aligning {i + 1}/{n_samples}: {entry['sample']}...")
            success, log, sample_out_dir = qm.run_star_align(
                entry, index_dir, align_dir, threads=align_threads,
                use_two_pass=use_two_pass, use_encode_options=use_encode_options,
                add_strand_field=add_strand_field,
            )
            mapping_rate = qm.parse_star_mapping_rate(sample_out_dir, entry["sample"]) if success else None
            results.append({
                "Sample": entry["sample"],
                "Status": "✅ Success" if success else "❌ Failed",
                "Uniquely Mapped %": f"{mapping_rate}%" if mapping_rate is not None else "—",
            })
            if not success:
                with st.expander(f"Error details for {entry['sample']}"):
                    st.code(log)

        progress_bar.progress(1.0, text="Alignment complete.")
        results_df = pd.DataFrame(results)
        st.dataframe(results_df, use_container_width=True, hide_index=True)

        if (results_df["Status"] == "✅ Success").any():
            pm.mark_step_complete(project, "quantification_complete")
            st.success("✅ Alignment & gene counting complete for at least one sample.")
            st.caption(
                "**Uniquely Mapped %** shows the % of reads that aligned "
                "to exactly one location in the genome. Rates above "
                "~70-80% are typically considered good; lower rates can "
                "indicate contamination, wrong species reference, or "
                "degraded samples."
            )
