"""
differential_expression_workspace.py

Differential Expression workspace: runs DESeq2 on a project's gene
counts matrix (built in the RNA Alignment & Counts workspace), with
support for:
  - User-controlled low-count gene filtering, with a live preview
  - Detection and resolution of ambiguous/missing metadata values in
    design/batch columns (blank cells, "none", "N/A", etc.) before they
    can silently break or skew the analysis
  - Multivariate experimental designs (multiple condition columns)
  - Batch effects handled as a design covariate (statistically correct
    approach), with optional batch-adjusted PCA for visualization
  - Multiple pairwise contrasts from a single fitted model
  - PCA, volcano plots, and Venn diagrams (2-3 contrasts) of DE gene
    overlap
  - Clean, ranked gene list export formatted for clusterProfiler

This picks up where alignment_workspace.py's Step 4 (counts matrix)
leaves off, reusing the same active project.

Design goal: same as the other workspaces — assume the user has little
to no bioinformatics/statistics background, explain each step in plain
language, while still exposing real statistical choices (filtering
thresholds, batch handling, contrasts) rather than hiding them.

This module is fully self-contained. All differential expression
development should happen here — editing this file has zero effect on
the other workspace modules.
"""

import math
import os

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import project_manager as pm
import deseq2_manager as dm

WORKSPACE_KEY = "bulk_rnaseq"


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_pca(pca_df, pct_variance, color_by):
    """Build a PC1 vs PC2 scatter plot, colored by the given metadata column."""
    pc1_label = f"PC1 ({pct_variance[0]}% variance)" if pct_variance else "PC1"
    pc2_label = f"PC2 ({pct_variance[1]}% variance)" if pct_variance else "PC2"

    fig = go.Figure()
    for group_value in sorted(pca_df[color_by].astype(str).unique()):
        subset = pca_df[pca_df[color_by].astype(str) == group_value]
        fig.add_trace(go.Scatter(
            x=subset["PC1"], y=subset["PC2"],
            mode="markers+text",
            text=subset["sample"],
            textposition="top center",
            name=group_value,
            marker=dict(size=12),
        ))
    fig.update_layout(
        xaxis_title=pc1_label,
        yaxis_title=pc2_label,
        height=500,
        legend_title=color_by,
    )
    return fig


def _plot_volcano(results_df, padj_threshold=0.05, lfc_threshold=1.0):
    """
    Build a volcano plot (log2FoldChange vs. -log10(padj)), coloring
    points by significance (padj < threshold AND |log2FC| >= threshold).
    """
    df = results_df.dropna(subset=["padj", "log2FoldChange"]).copy()
    # -log10(padj): clip padj away from exactly 0 first (log10(0) is
    # undefined) using a tiny floor value, so a gene with a vanishingly
    # small p-value still gets a large-but-finite y-position rather than
    # causing a math error.
    df["neg_log10_padj"] = -df["padj"].clip(lower=1e-300).apply(math.log10)
    df["significant"] = (df["padj"] < padj_threshold) & (df["log2FoldChange"].abs() >= lfc_threshold)

    fig = go.Figure()
    for is_sig, label, color in [(True, "Significant", "#d62728"), (False, "Not significant", "#7f7f7f")]:
        subset = df[df["significant"] == is_sig]
        fig.add_trace(go.Scatter(
            x=subset["log2FoldChange"], y=subset["neg_log10_padj"],
            mode="markers",
            name=label,
            marker=dict(size=6, color=color, opacity=0.6),
            text=subset["gene_id"],
            hovertemplate="%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}<extra></extra>",
        ))

    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")
    # Horizontal significance threshold line, drawn at -log10(padj_threshold)
    # — e.g. padj_threshold=0.05 -> line at y=1.301, matching the same
    # -log10 transform used for the plotted points above.
    padj_threshold_line_y = -math.log10(padj_threshold)
    fig.add_hline(y=padj_threshold_line_y, line_dash="dash", line_color="gray")

    fig.update_layout(
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(adjusted p-value)",
        height=550,
    )
    return fig


# Standard 2-circle and 3-circle Venn layout coordinates (unit circles at
# fixed classic positions), used to draw a lightweight Venn diagram with
# Plotly shapes rather than depending on an external Venn-plotting
# library (e.g. matplotlib-venn), which may not be installed.
_VENN_2_CIRCLES = [
    {"x0": -0.6, "y0": -0.5, "x1": 0.9, "y1": 1.0, "label_x": -0.35, "label_y": 0.25},
    {"x0": -0.3, "y0": -0.5, "x1": 1.2, "y1": 1.0, "label_x": 0.65, "label_y": 0.25},
]
_VENN_3_CIRCLES = [
    {"x0": -0.9, "y0": -0.6, "x1": 0.7, "y1": 1.0, "label_x": -0.55, "label_y": 0.55},
    {"x0": -0.3, "y0": -0.6, "x1": 1.3, "y1": 1.0, "label_x": 0.85, "label_y": 0.55},
    {"x0": -0.6, "y0": -1.1, "x1": 1.0, "y1": 0.5, "label_x": 0.15, "label_y": -0.65},
]


def _plot_venn(counts, names):
    """
    Draw a simple 2- or 3-circle Venn diagram using Plotly shapes, with
    region counts annotated at approximate positions. This intentionally
    keeps to a lightweight, dependency-free approach (no matplotlib-venn
    required) — good enough for a quick visual overlap check, though not
    as precisely proportional as a dedicated Venn library.
    """
    n = len(names)
    circles = _VENN_2_CIRCLES if n == 2 else _VENN_3_CIRCLES
    colors = ["rgba(99,110,250,0.4)", "rgba(239,85,59,0.4)", "rgba(0,204,150,0.4)"]

    fig = go.Figure()
    for i, name in enumerate(names):
        c = circles[i]
        fig.add_shape(
            type="circle", xref="x", yref="y",
            x0=c["x0"], y0=c["y0"], x1=c["x1"], y1=c["y1"],
            fillcolor=colors[i], line=dict(color=colors[i].replace("0.4", "1.0")),
        )
        fig.add_annotation(x=c["label_x"], y=c["label_y"] + 0.35, text=f"<b>{name}</b>", showarrow=False)

    # Annotate each region's count at a reasonable position. Exact
    # centroid placement for arbitrary overlaps is nontrivial with plain
    # shapes, so positions below are hand-tuned for the classic 2/3-circle
    # layout rather than computed generically.
    if n == 2:
        a, b = names
        positions = {
            a: (-0.35, 0.25), b: (0.65, 0.25), f"{a} & {b}": (0.15, 0.25),
        }
    else:
        a, b, c_ = names
        positions = {
            a: (-0.55, 0.15), b: (0.85, 0.15), c_: (0.15, -0.65),
            f"{a} & {b}": (0.15, 0.35),
            f"{a} & {c_}": (-0.15, -0.25),
            f"{b} & {c_}": (0.45, -0.25),
            f"{a} & {b} & {c_}": (0.15, -0.05),
        }

    for label, (x, y) in positions.items():
        fig.add_annotation(x=x, y=y, text=str(counts.get(label, 0)), showarrow=False, font=dict(size=16))

    fig.update_xaxes(visible=False, range=[-1.3, 1.6])
    fig.update_yaxes(visible=False, range=[-1.4, 1.3], scaleanchor="x", scaleratio=1)
    fig.update_layout(height=500, showlegend=False, plot_bgcolor="white")
    return fig


# ---------------------------------------------------------------------------
# Ambiguous / missing value resolution UI
# ---------------------------------------------------------------------------

def _render_missing_value_resolution(meta_df, columns_to_check):
    """
    For each column in columns_to_check, detect ambiguous values (true
    blanks/NaN, or missing-like text like "none"/"N/A") and — if any are
    found — show a resolution UI letting the user assign each one a real
    category label, or mark affected samples for exclusion.

    Returns (resolved_meta_df, all_resolved: bool). If any ambiguous
    value in any checked column hasn't been explicitly resolved yet,
    all_resolved is False and the caller should block proceeding to the
    next step — this prevents a real NaN or an unresolved "none" from
    silently reaching DESeq2, where it would either error out or (worse)
    be misinterpreted.
    """
    any_ambiguous = dm.has_any_ambiguous_values(meta_df, columns_to_check)
    if not any_ambiguous:
        return meta_df, True

    st.warning(
        "⚠️ We found some blank or ambiguous values in your selected "
        "column(s) below. Please tell us what each one means before "
        "continuing — for example, a blank cell or \"none\" in an "
        "antibiotic column often means \"no antibiotic was used\" (a "
        "real control group), not missing data."
    )

    resolved_df = meta_df.copy()
    all_resolved = True

    for column in columns_to_check:
        detection = dm.detect_ambiguous_values(meta_df[column])
        if detection["true_na_count"] == 0 and not detection["missing_like_values"]:
            continue

        st.markdown(f"**Column: `{column}`**")
        resolution_map = {}
        exclude_samples = []

        if detection["true_na_count"] > 0:
            affected_samples = meta_df.loc[meta_df[column].isna(), "sample"].tolist()
            st.caption(
                f"🔲 **{detection['true_na_count']} blank cell(s)** found "
                f"— affected sample(s): {', '.join(affected_samples)}"
            )
            blank_choice = st.radio(
                f"What does a blank value in `{column}` mean?",
                ["This represents a real category (e.g. no treatment/control)", "These samples should be excluded from analysis"],
                key=f"blank_choice_{column}",
                horizontal=False,
            )
            if blank_choice.startswith("This represents"):
                label = st.text_input(
                    f"What should this category be called?",
                    value=f"no_{column}",
                    key=f"blank_label_{column}",
                )
                if label.strip():
                    resolution_map["__NaN__"] = label.strip()
                else:
                    all_resolved = False
            else:
                exclude_samples.extend(affected_samples)

        for ambiguous_val, count in detection["missing_like_values"].items():
            affected_samples = meta_df.loc[
                meta_df[column].astype(str).str.strip().str.lower() == ambiguous_val.strip().lower(),
                "sample"
            ].tolist()
            st.caption(
                f"🔲 **\"{ambiguous_val}\"** found in {count} sample(s) "
                f"— affected sample(s): {', '.join(affected_samples)}"
            )
            value_choice = st.radio(
                f"What does \"{ambiguous_val}\" mean in `{column}`?",
                [f"This represents a real category (e.g. no treatment/control)", "These samples should be excluded from analysis"],
                key=f"value_choice_{column}_{ambiguous_val}",
                horizontal=False,
            )
            if value_choice.startswith("This represents"):
                label = st.text_input(
                    f"What should \"{ambiguous_val}\" be relabeled as?",
                    value=f"no_{column}",
                    key=f"value_label_{column}_{ambiguous_val}",
                )
                if label.strip():
                    resolution_map[ambiguous_val] = label.strip()
                else:
                    all_resolved = False
            else:
                exclude_samples.extend(affected_samples)

        resolved_df = dm.apply_missing_value_resolution(resolved_df, column, resolution_map, exclude_samples)
        st.markdown("---")

    if all_resolved:
        st.success("✅ All ambiguous values have been resolved.")
        with st.expander("Preview resolved metadata"):
            st.dataframe(resolved_df, use_container_width=True, hide_index=True)

    return resolved_df, all_resolved


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("📈 Differential Expression (DESeq2)")
    st.markdown(
        "This workspace compares gene expression between your "
        "experimental conditions to find genes that are significantly "
        "up- or down-regulated. **No statistics background required** — "
        "follow the steps below in order, and each choice is explained "
        "as you go."
    )
    st.markdown("---")

    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)
    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    st.markdown("---")

    counts_matrix_done = pm.has_completed_step(project, "counts_matrix_complete")
    if not counts_matrix_done:
        st.warning(
            "⚠️ This project doesn't have a gene counts matrix yet. Go "
            "to the **🧮 RNA Alignment & Counts** page and complete "
            "Step 4 (build the counts matrix) first."
        )
        if st.button("⬅️ Go to RNA Alignment & Counts", key="de_gate_back_btn"):
            st.session_state["nav_request"] = "🧮 RNA Alignment & Counts"
            st.rerun()
        return

    counts_matrix_path = pm.counts_matrix_path(project)
    metadata_path = pm.metadata_path(project)

    if not os.path.exists(counts_matrix_path) or not os.path.exists(metadata_path):
        st.error("⚠️ Could not find the counts matrix or metadata file for this project.")
        return

    counts_df = pd.read_csv(counts_matrix_path)
    meta_df = pd.read_csv(metadata_path)

    st.success(f"✅ Using counts matrix ({len(counts_df):,} genes) and metadata ({len(meta_df)} samples) from project `{project}`.")

    deseq2_out_dir = pm.deseq2_output_dir(project)
    deseq2_work_dir = pm.deseq2_work_dir(project)

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 1: Low-count gene filtering
    # -----------------------------------------------------------------
    st.header("Step 1: Remove Low-Count Genes")

    with st.expander("ℹ️ Why remove low-count genes? (click to learn more)"):
        st.markdown(
            "Genes with very few reads across your samples don't carry "
            "enough information to reliably detect real differences — "
            "including them mostly adds statistical noise and can make "
            "DESeq2's multiple-testing correction less powerful for the "
            "genes that *do* have enough data. A common approach is to "
            "keep only genes with a minimum number of reads in at least "
            "a handful of samples."
        )

    col1, col2 = st.columns(2)
    with col1:
        min_count = st.number_input(
            "Minimum read count:", min_value=0, value=10, step=1,
            help="A gene must have at least this many reads in a sample for that sample to 'count' toward the threshold below.",
            key="min_count_input",
        )
    with col2:
        min_samples = st.number_input(
            "Minimum number of samples:", min_value=1, value=min(2, len(meta_df)), step=1,
            max_value=len(meta_df),
            help="The gene must meet the read count threshold in at least this many samples to be kept.",
            key="min_samples_input",
        )

    filter_preview = dm.preview_low_count_filter(counts_df, min_count, min_samples)
    st.caption(
        f"With these settings: **{filter_preview['genes_kept']:,} of "
        f"{filter_preview['total_genes']:,} genes** would be kept "
        f"({filter_preview['pct_kept']}%), removing "
        f"{filter_preview['genes_removed']:,} low-count genes."
    )

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 2: Experimental design
    # -----------------------------------------------------------------
    st.header("Step 2: Set Up Your Experimental Design")

    with st.expander("ℹ️ What is an experimental design, and what about batch effects? (click to learn more)"):
        st.markdown(
            "Your **design** tells DESeq2 which column(s) in your "
            "metadata represent the experimental conditions you want to "
            "compare (e.g. `treatment`). If you have more than one "
            "relevant factor (e.g. `genotype` AND `treatment`), you can "
            "select multiple columns — DESeq2 will account for all of "
            "them together.\n\n"
            "**Batch effects** happen when samples processed at "
            "different times, by different people, or on different "
            "equipment show systematic differences unrelated to your "
            "actual biological question. If you have a column "
            "identifying batches (e.g. `batch`, `sequencing_run`, "
            "`prep_date`), select it below — DESeq2 will account for "
            "batch variation *while* estimating your condition effects, "
            "which is the statistically correct way to handle this "
            "(rather than 'correcting' the data directly, which can "
            "discard real information).\n\n"
            "**Important — what this selection does and doesn't do:** "
            "your project's full metadata file (every column, every "
            "sample) stays exactly as it is — this step does **not** "
            "remove or overwrite anything. You're only telling DESeq2 "
            "*for this specific analysis run* which column(s) to treat "
            "as the experimental condition(s). If you later want to "
            "test a **different** condition (e.g. you tested "
            "`antibiotic` this time, but later want to test `media` "
            "instead), or test conditions in a **different combination** "
            "(e.g. `antibiotic` alone, then `antibiotic` + `media` "
            "together), you'll come back to this exact step and simply "
            "select different column(s) — each combination you want to "
            "test is its own pass through Steps 2-4, and each will save "
            "its own separate results without affecting the others."
        )

    available_columns = [c for c in meta_df.columns if c != "sample"]

    st.caption(
        "💡 Tip: if your metadata came from an auto-fill (e.g. SRA "
        "lookup), it often includes extra descriptive columns (like a "
        "sample title, strain, or organism name) that **aren't** "
        "meant to be experimental conditions — only select the "
        "column(s) that actually vary in the way you want to compare "
        "(e.g. `treatment`, `antibiotic`, `genotype`)."
    )
    design_columns = st.multiselect(
        "Which column(s) represent your experimental condition(s) of interest?",
        options=available_columns,
        default=[],
        help="Select more than one for a multivariate design (e.g. genotype + treatment). Leave unselected if unsure, then check the validation below before proceeding.",
        key="design_columns_select",
    )

    batch_column_options = ["(none)"] + [c for c in available_columns if c not in design_columns]
    batch_column_choice = st.selectbox(
        "Does this project have a batch effect column? (optional)",
        options=batch_column_options,
        key="batch_column_select",
    )
    batch_column = None if batch_column_choice == "(none)" else batch_column_choice

    # --- Live preview of the selected columns ---
    # Shown right after the pickers so the user can visually confirm
    # what's actually in each selected column (values, spelling,
    # capitalization, blanks) before committing — rather than having to
    # scroll back up to the raw metadata or guess from memory.
    preview_columns = ["sample"] + design_columns
    if batch_column:
        preview_columns.append(batch_column)

    if design_columns or batch_column:
        st.caption("Preview of your selected column(s) (first 5 rows):")
        st.dataframe(meta_df[preview_columns].head(5), use_container_width=True, hide_index=True)

    if not design_columns:
        st.warning("⚠️ Select at least one design column to continue.")
        return

    # --- Validate column choices before allowing progression ---
    # Catches two common, easy-to-miss problems immediately (rather than
    # after a multi-minute DESeq2 run fails):
    #   1. A column with the same value for every sample (e.g. "strain"
    #      in a single-species study) — DESeq2 hard-errors on this.
    #   2. A column that's essentially a per-sample identifier (as many
    #      unique values as samples, e.g. a free-text sample
    #      description) — technically won't crash, but is never a
    #      meaningful design factor and usually indicates the wrong
    #      column was selected.
    columns_to_validate = list(design_columns)
    if batch_column:
        columns_to_validate.append(batch_column)

    validation = dm.validate_design_columns(meta_df, columns_to_validate)
    constant_columns = [col for col, info in validation.items() if info["is_constant"]]
    identifier_columns = [col for col, info in validation.items() if info["is_identifier_like"]]

    if constant_columns:
        st.error(
            f"⚠️ **{', '.join(constant_columns)}** has the exact same "
            "value for every sample in this project, so it can't be "
            "used as a design/batch factor — DESeq2 requires some "
            "variation across samples to estimate an effect. Please "
            "deselect it above."
        )
        return

    if identifier_columns:
        st.warning(
            f"⚠️ **{', '.join(identifier_columns)}** has a different "
            "value for almost every sample (like a free-text "
            "description or ID), which usually means it's not a real "
            "experimental condition — including it will likely make "
            "your model meaningless. If this wasn't intentional, "
            "deselect it above."
        )

    # --- Check for sufficient replication before proceeding ---
    # DESeq2 needs at least 2 samples per group (combination of design
    # + batch column levels) to estimate dispersion (within-group
    # variability) at all — with only 1 sample per group, there's no
    # data left to distinguish real biological signal from noise, and
    # DESeq2 hard-errors during dispersion estimation
    # ("checkForExperimentalReplicates"). This is a very common issue
    # when metadata comes from a small, hand-picked subset of samples
    # (e.g. a handful of runs pulled from a larger public dataset) where
    # each one happens to represent a different, unreplicated condition.
    replication = dm.check_replication(meta_df, design_columns, batch_column)
    if not replication["is_valid"]:
        st.error(
            f"⚠️ The following group(s) only have 1 sample: "
            f"**{', '.join(replication['under_replicated_groups'])}**. "
            "DESeq2 needs **at least 2 samples per group** to estimate "
            "natural variability and distinguish real differences from "
            "noise — with only 1 sample, this isn't statistically "
            "possible.\n\n"
            "To fix this, you'll need to either add more samples for "
            "the under-replicated group(s) (e.g. download additional "
            "runs from the same study/condition), or choose a different "
            "condition/combination of conditions where every group has "
            "at least 2 samples."
        )
        with st.expander("View sample counts per group"):
            counts_df = pd.DataFrame([
                {"Group": label, "Samples": count}
                for label, count in replication["group_counts"].items()
            ])
            st.dataframe(counts_df, use_container_width=True, hide_index=True)
        return

    # --- Detect and resolve ambiguous/missing values before proceeding ---
    # Checked across every column actually used in the model (design
    # columns + batch, if selected), since an unresolved blank/"none"
    # value in any of these would either crash DESeq2 or silently create
    # a spurious factor level.
    columns_to_check = list(design_columns)
    if batch_column:
        columns_to_check.append(batch_column)

    meta_df, values_resolved = _render_missing_value_resolution(meta_df, columns_to_check)
    if not values_resolved:
        st.warning("⚠️ Please resolve all ambiguous values above before continuing.")
        return

    # --- Interaction terms (optional, only relevant with 2+ design columns) ---
    interaction_terms = []
    if len(design_columns) >= 2:
        with st.expander("ℹ️ What is an interaction term? (click to learn more)"):
            st.markdown(
                "An **interaction term** tests whether the effect of "
                "one variable *depends on* another. For example, "
                "instead of just asking \"does treatment matter?\" and "
                "\"does genotype matter?\" separately, an interaction "
                "between `genotype` and `treatment` asks: \"is the "
                "*effect* of treatment different depending on genotype?"
                "\" (e.g. treatment works in wild-type mice but not in "
                "knockout mice).\n\n"
                "Only add an interaction if you have a specific "
                "biological reason to expect one — adding unnecessary "
                "interaction terms reduces statistical power for your "
                "main comparisons, especially with smaller sample sizes."
            )
        interaction_options = dm.build_interaction_term_options(design_columns)
        interaction_terms = st.multiselect(
            "Include interaction term(s) in the model? (optional)",
            options=interaction_options,
            default=[],
            key="interaction_terms_select",
        )

    full_formula_terms = dm.build_full_formula_terms(design_columns, batch_column, interaction_terms)
    st.caption(f"**Full model formula:** ~ {' + '.join(full_formula_terms)}")

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 3: Choose Your Analysis Type
    # -----------------------------------------------------------------
    st.header("Step 3: Choose Your Analysis Type")

    with st.expander("ℹ️ Wald test vs. LRT (Likelihood Ratio Test) — which should I use? (click to learn more)"):
        st.markdown(
            "DESeq2 offers two different statistical approaches, and "
            "which one you want depends on your question:\n\n"
            "**🎯 Wald test (pairwise comparisons)** — tests one "
            "specific comparison at a time, e.g. \"is drugA different "
            "from control?\". This is what most people want when they "
            "have a clear reference group and want to know how each "
            "other group differs from it. You can run several pairwise "
            "comparisons (contrasts) from one model fit.\n\n"
            "**📐 LRT (Likelihood Ratio Test)** — DESeq2's equivalent of "
            "an **ANOVA**. Instead of comparing two specific groups, it "
            "asks a broader question: \"does this factor matter *at "
            "all*, across all of its levels?\" or \"does this "
            "interaction matter at all?\". It works by comparing a "
            "**full model** (with the term you're testing) against a "
            "**reduced model** (without it) and testing whether removing "
            "that term makes the model fit significantly worse. This is "
            "especially useful when a factor has 3+ levels and you want "
            "one overall p-value per gene (\"does this gene respond "
            "differently to *any* of the treatments?\") rather than "
            "several separate pairwise p-values, or when you want to "
            "test whether an **interaction term** matters at all.\n\n"
            "**Rule of thumb:** if you're asking \"how does A differ "
            "from B\", use Wald. If you're asking \"does this variable "
            "matter overall\" or \"does this interaction matter\", use "
            "LRT."
        )

    test_type_choice = st.radio(
        "Which type of test do you want to run?",
        ["Wald test (pairwise comparisons)", "LRT (omnibus / ANOVA-style test)"],
        key="test_type_radio",
    )
    test_type = "wald" if test_type_choice.startswith("Wald") else "lrt"

    contrasts = []
    lrt_tests = []

    if test_type == "wald":
        st.subheader("Set Up Your Pairwise Comparisons")

        primary_column = st.selectbox(
            "Which design column should contrasts be based on?",
            options=design_columns,
            key="primary_contrast_column",
        )

        levels = sorted(meta_df[primary_column].astype(str).dropna().unique().tolist())
        if len(levels) < 2:
            st.error(f"⚠️ Column '{primary_column}' needs at least 2 distinct values to define a contrast.")
            return

        reference_level = st.selectbox(
            "Which level is your reference/control group?",
            options=levels,
            key="reference_level_select",
        )
        comparison_levels = st.multiselect(
            "Which level(s) should be compared against the reference?",
            options=[lv for lv in levels if lv != reference_level],
            default=[lv for lv in levels if lv != reference_level],
            key="comparison_levels_select",
        )

        contrasts = [
            {
                "name": f"{lv}_vs_{reference_level}",
                "column": primary_column,
                "level1": lv,
                "level2": reference_level,
            }
            for lv in comparison_levels
        ]

        if contrasts:
            st.markdown("**Contrasts that will be run:**")
            st.dataframe(pd.DataFrame(contrasts)[["name", "level1", "level2"]].rename(
                columns={"name": "Contrast", "level1": "Comparison", "level2": "Reference"}
            ), use_container_width=True, hide_index=True)
        else:
            st.warning("⚠️ Select at least one comparison level to continue.")
            return

    else:  # LRT
        st.subheader("Set Up Your Omnibus Test(s)")
        st.markdown(
            "Select which term(s) in your model you want to test "
            "overall. Each one you select becomes its own LRT test, "
            "comparing the full model against a reduced model with that "
            "term removed."
        )

        lrt_term_options = list(design_columns) + interaction_terms
        terms_to_test = st.multiselect(
            "Which term(s) should be tested?",
            options=lrt_term_options,
            default=lrt_term_options[:1] if lrt_term_options else [],
            help="Select a design column to test its overall effect (like an ANOVA F-test), or an interaction term to test whether it matters at all.",
            key="lrt_terms_select",
        )

        lrt_tests = [
            {
                "name": f"{term.replace(':', '_x_')}_effect",
                "reduced_terms": dm.build_reduced_formula_terms(full_formula_terms, [term]),
            }
            for term in terms_to_test
        ]

        if lrt_tests:
            st.markdown("**LRT test(s) that will be run:**")
            preview_rows = [
                {"Test": t["name"], "Full Model": " + ".join(full_formula_terms), "Reduced Model": " + ".join(t["reduced_terms"]) if t["reduced_terms"] else "(intercept only)"}
                for t in lrt_tests
            ]
            st.dataframe(pd.DataFrame(preview_rows), use_container_width=True, hide_index=True)
            st.caption(
                "⚠️ Note: for a factor with more than 2 levels, LRT's "
                "`log2FoldChange` column in the results reflects only "
                "one internal reference comparison and should be "
                "interpreted cautiously — the **p-value/adjusted "
                "p-value** is the meaningful output of an LRT test (it "
                "tells you whether the term matters overall), similar "
                "to an ANOVA F-test p-value."
            )
        else:
            st.warning("⚠️ Select at least one term to test to continue.")
            return

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 4: Run DESeq2
    # -----------------------------------------------------------------
    st.header("Step 4: Run DESeq2")

    if not dm.deseq2_tools_available():
        st.error(
            "⚠️ Rscript was not found on this system. R with DESeq2 "
            "needs to be installed in your environment (included in the "
            "project's Dockerfile) before this step can run."
        )
        return

    # Since the resolved metadata (with ambiguous values fixed and/or
    # samples excluded) needs to be what DESeq2 actually reads, we write
    # it out to the same metadata_path the R script will load — this
    # keeps run_deseq2_analysis's interface simple (it just takes a path)
    # while ensuring resolved data is what's actually used.
    resolved_metadata_path = os.path.join(deseq2_work_dir, "resolved_metadata.csv")
    os.makedirs(deseq2_work_dir, exist_ok=True)
    meta_df.to_csv(resolved_metadata_path, index=False)

    results_already_exist = len(dm.list_contrast_results(deseq2_out_dir)) > 0
    run_label = "🔄 Re-run DESeq2" if results_already_exist else "🚀 Run DESeq2 Analysis"

    if results_already_exist and not st.session_state.get("_deseq2_run_clicked"):
        st.success("✅ DESeq2 has already been run for this project.")

    if st.button(run_label, key="run_deseq2_btn"):
        st.session_state["_deseq2_run_clicked"] = True

        with st.spinner("Running DESeq2... this may take a few minutes."):
            success, log = dm.run_deseq2_analysis(
                counts_matrix_path, resolved_metadata_path, design_columns,
                batch_column, min_count, min_samples, deseq2_out_dir, deseq2_work_dir,
                test_type=test_type, interaction_terms=interaction_terms,
                contrasts=contrasts, lrt_tests=lrt_tests,
            )

        if not success:
            st.error("DESeq2 failed. Details below:")
            st.code(log)
            return
        else:
            st.success("✅ DESeq2 completed successfully.")
            st.code(log)
            pm.mark_step_complete(project, "deseq2_complete")
            pm.save_deseq2_config(project, {
                "design_columns": design_columns,
                "batch_column": batch_column,
                "interaction_terms": interaction_terms,
                "test_type": test_type,
                "contrasts": contrasts,
                "lrt_tests": lrt_tests,
                "min_count": min_count,
                "min_samples": min_samples,
            })
            results_already_exist = True

    if not results_already_exist:
        return

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 5: Explore Results
    # -----------------------------------------------------------------
    st.header("Step 5: Explore Your Results")

    available_contrasts = dm.list_contrast_results(deseq2_out_dir)
    if not available_contrasts:
        st.info("No results found yet — run DESeq2 above.")
        return

    # --- PCA ---
    st.subheader("🔬 PCA: Sample Similarity")
    with st.expander("ℹ️ How to read this plot"):
        st.markdown(
            "Each point is one sample. Samples that are more similar in "
            "overall gene expression appear closer together. If your "
            "samples cluster clearly by condition, that's a good sign — "
            "if they instead cluster by batch, that suggests a strong "
            "batch effect (which DESeq2 has already accounted for "
            "statistically in your results, via the design formula)."
        )

    pca_df, pct_var = dm.read_pca_coordinates(deseq2_out_dir)
    if pca_df is not None:
        color_options = [c for c in pca_df.columns if c not in ("PC1", "PC2", "PC3", "PC4", "sample")]
        color_by = st.selectbox("Color points by:", options=color_options, key="pca_color_select")
        st.plotly_chart(_plot_pca(pca_df, pct_var, color_by), use_container_width=True)

        if batch_column:
            pca_adj_df, pct_var_adj = dm.read_pca_coordinates(deseq2_out_dir, batch_adjusted=True)
            if pca_adj_df is not None:
                st.caption(
                    "**Batch-adjusted view** (for visualization only — "
                    "the actual statistical test above already accounts "
                    "for batch via the design formula, this view just "
                    "helps you *see* how much batch was contributing to "
                    "sample separation):"
                )
                st.plotly_chart(_plot_pca(pca_adj_df, pct_var_adj, color_by), use_container_width=True)

    st.markdown("---")

    # --- Per-contrast results + volcano ---
    st.subheader("📊 Results by Contrast")
    selected_contrast = st.selectbox("Select a contrast to view:", options=available_contrasts, key="view_contrast_select")

    results_df = dm.read_contrast_results(deseq2_out_dir, selected_contrast)
    if results_df is not None:
        padj_cutoff = st.slider("Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key="padj_slider")
        lfc_cutoff = st.slider("Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key="lfc_slider")

        n_sig = ((results_df["padj"] < padj_cutoff) & (results_df["log2FoldChange"].abs() >= lfc_cutoff)).sum()
        st.caption(f"**{n_sig:,}** significant genes at padj < {padj_cutoff} and |log2FC| ≥ {lfc_cutoff}.")

        st.plotly_chart(_plot_volcano(results_df, padj_cutoff, lfc_cutoff), use_container_width=True)

        st.dataframe(results_df.head(200), use_container_width=True, hide_index=True)
        st.caption("Showing top 200 rows (sorted by adjusted p-value). Download the full table below.")

        csv_bytes = results_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            f"⬇️ Download Full Results ({selected_contrast}).csv",
            data=csv_bytes,
            file_name=f"deseq2_results_{selected_contrast}.csv",
            mime="text/csv",
            key="download_results_btn",
        )

        # --- clusterProfiler export ---
        with st.expander("🧬 Export for clusterProfiler (gene ontology enrichment)"):
            st.markdown(
                "This exports a clean, ranked gene list (gene ID, "
                "log2 fold change, adjusted p-value) suitable for use "
                "with R's `clusterProfiler` package "
                "(`enrichGO`/`gseGO`).\n\n"
                "⚠️ **Note:** `clusterProfiler` typically expects Entrez "
                "IDs or gene symbols rather than Ensembl transcript/gene "
                "IDs directly, depending on which organism annotation "
                "package you use (e.g. `org.Hs.eg.db` for human, "
                "`org.Mm.eg.db` for mouse). You may need to convert gene "
                "IDs first using clusterProfiler's own `bitr()` function "
                "or `biomaRt` — this depends on your organism and isn't "
                "something we can determine automatically."
            )
            export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrast)
            if export_df is not None:
                st.dataframe(export_df.head(20), use_container_width=True, hide_index=True)
                export_csv = export_df.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "⬇️ Download clusterProfiler-Ready Gene List (.csv)",
                    data=export_csv,
                    file_name=f"clusterprofiler_ranked_{selected_contrast}.csv",
                    mime="text/csv",
                    key="download_clusterprofiler_btn",
                )

    st.markdown("---")

    # --- Venn diagram across contrasts ---
    st.subheader("🔵 Overlap Between Contrasts (Venn Diagram)")
    if len(available_contrasts) < 2:
        st.info("Run at least 2 contrasts to see an overlap diagram.")
    else:
        venn_contrasts = st.multiselect(
            "Select 2 or 3 contrasts to compare:",
            options=available_contrasts,
            default=available_contrasts[:min(2, len(available_contrasts))],
            key="venn_contrast_select",
        )
        if len(venn_contrasts) not in (2, 3):
            st.info("Please select exactly 2 or 3 contrasts (Venn diagrams beyond 3 sets become unreadable).")
        else:
            venn_padj = st.slider("Significance threshold for overlap:", 0.01, 0.20, 0.05, step=0.01, key="venn_padj_slider")
            gene_sets = dm.get_significant_gene_sets(deseq2_out_dir, venn_contrasts, padj_threshold=venn_padj)
            counts, regions = dm.compute_venn_regions(gene_sets)
            st.plotly_chart(_plot_venn(counts, venn_contrasts), use_container_width=True)

            with st.expander("View overlapping gene lists"):
                for label, gene_set in regions.items():
                    st.markdown(f"**{label}** ({len(gene_set)} genes)")
                    if gene_set:
                        st.text(", ".join(sorted(gene_set)[:50]) + (" ..." if len(gene_set) > 50 else ""))

    st.markdown("---")
    st.success(
        f"🎉 Project `{project}` has completed differential expression "
        "analysis. Results, PCA, volcano plots, and Venn diagrams are "
        "available above, along with clusterProfiler-ready exports for "
        "downstream gene ontology enrichment analysis."
    )
