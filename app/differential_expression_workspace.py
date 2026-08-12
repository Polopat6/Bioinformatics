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
  - Per-plot customization (title, axis labels, fonts, colors, themes,
    draggable legend) and PDF export with size selection

This picks up where alignment_workspace.py's Step 4 (counts matrix)
leaves off, reusing the same active project.

Design goal: same as the other workspaces -- assume the user has little
to no bioinformatics/statistics background, explain each step in plain
language, while still exposing real statistical choices (filtering
thresholds, batch handling, contrasts) rather than hiding them.

This module is fully self-contained. All differential expression
development should happen here -- editing this file has zero effect on
the other workspace modules.
"""
import math
import os
import re

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import project_manager as pm
import deseq2_manager as dm
import reference_manager as rm
import gene_id_mapper as gim

WORKSPACE_KEY = "bulk_rnaseq"

# ---------------------------------------------------------------------------
# Plot customization + export helpers
# ---------------------------------------------------------------------------
#
# These are intentionally decoupled from the figure-building functions
# below (_plot_pca, _plot_volcano, _plot_venn): each "_plot_*" function
# only builds the raw data traces, and _apply_plot_style() is layered
# on top afterward to handle cosmetics (title, axis labels, font,
# theme, legend placement, per-group colors). This keeps the data
# logic and the styling logic independent, so styling controls can be
# reused across every plot in this workspace without duplicating code.

# Curated subset of Plotly's built-in named templates -- all valid
# values for `fig.update_layout(template=...)` without needing any
# extra dependency.
PLOTLY_THEME_OPTIONS = {
    "Default": "plotly",
    "Plotly White": "plotly_white",
    "Plotly Dark": "plotly_dark",
    "Seaborn": "seaborn",
    "Simple White": "simple_white",
    "ggplot2-style": "ggplot2",
    "Presentation": "presentation",
}

# Web-safe font families that render consistently across browsers
# without requiring any font files to be installed locally.
FONT_FAMILY_OPTIONS = [
    "Arial", "Helvetica", "Times New Roman", "Courier New", "Georgia", "Verdana",
]

# Matches Plotly's own default "plotly" template colorway exactly, so
# that if a user never touches the color pickers, the plot looks
# identical to before this customization panel existed.
_DEFAULT_COLOR_PALETTE = [
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
]

# Common output sizes for PDF export, in pixels. "Print quality"
# presets are sized at 300 DPI for their nominal paper dimensions.
PDF_SIZE_PRESETS = {
    "Small (600 x 450 px)": (600, 450),
    "Medium (900 x 675 px)": (900, 675),
    "Large (1200 x 900 px)": (1200, 900),
    "Presentation Widescreen (1920 x 1080 px)": (1920, 1080),
    "US Letter Portrait, print quality (2550 x 3300 px)": (2550, 3300),
    "US Letter Landscape, print quality (3300 x 2550 px)": (3300, 2550),
    "A4 Portrait, print quality (2480 x 3508 px)": (2480, 3508),
    "Custom size": None,
}

# Plotly.js "edits" config that enables click-and-drag repositioning
# directly on the rendered chart in the browser, in addition to the
# preset position dropdowns offered in the style panel:
#   - legendPosition: drag the legend anywhere on the chart.
#   - annotationPosition: drag any annotation -- this is what makes the
#     plot title draggable (see _apply_plot_style, which renders the
#     title as an annotation rather than Plotly's native layout.title,
#     since Plotly.js has no built-in way to drag-reposition a real
#     chart title) AND, as a side effect, makes every other annotation
#     in a figure draggable too -- including the Venn diagram's circle
#     labels and region-count labels, letting the user drag those
#     directly on the chart.
#   - titleText: lets the user double-click the title annotation/axis
#     titles to edit their text directly on the chart, as a bonus on
#     top of the text inputs already offered in the style panel.
#   - annotationTail: per Plotly.js's own "edits" reference, this "has
#     only an effect for annotations with arrows" and "enables changing
#     the length and direction of the arrow" -- i.e. dragging moves
#     just the label/text end (ax, ay) while the arrowhead (x, y,
#     pinned to the gene's actual data point) stays put. This is what
#     makes each volcano gene label's position adjustable (see
#     _build_gene_label_annotations) while its arrow keeps pointing at
#     the correct dot. (annotationPosition also has an effect on
#     arrowed annotations -- moving the whole arrow+label together,
#     anchor included -- but annotationTail is the one that gives
#     independent label repositioning, which is what's wanted here.)
_PLOTLY_CHART_CONFIG = {
    "displaylogo": False,
    "edits": {
        "legendPosition": True,
        "annotationPosition": True,
        "annotationTail": True,
        "titleText": True,
    },
}


def _pdf_export_available():
    "Check whether the `kaleido` package (Plotly's static image export backend) is installed."
    try:
        import kaleido  # noqa: F401
        return True
    except ImportError:
        return False


def _render_plotly_chart(fig):
    "Render a Plotly figure with click-and-drag legend repositioning enabled."
    st.plotly_chart(fig, use_container_width=True, config=_PLOTLY_CHART_CONFIG)


def _render_plot_style_controls(key_prefix, group_values=None, default_colors=None,
                                 show_legend_controls=True):
    """
    Render a reusable "customize this plot's appearance" panel: a
    custom title, axis label overrides, font family/size, a theme
    picker, a starting legend-position preset, and (if group_values is
    given) one color picker per group/category so the user can recolor
    specific series to their preference.

    group_values: list of distinct group/category labels this plot has
        one trace per (e.g. unique metadata levels for a PCA plot, or
        ["Significant", "Not significant"] for a volcano plot). Pass
        None/empty to skip the color-picker section entirely (e.g. for
        plots without a meaningful per-group color scheme, like the
        Venn diagram).
    default_colors: optional dict {group_label: hex_color} giving the
        *existing* default color for each group, so that if the user
        never touches a color picker, the plot's appearance doesn't
        change from before this panel existed. Falls back to Plotly's
        own default qualitative color sequence for any group not
        listed here.
    show_legend_controls: set False for plots that don't display a
        legend at all (e.g. the Venn diagram, which uses text
        annotations instead).

    Returns a style options dict consumed by _apply_plot_style().
    """
    default_colors = default_colors or {}

    with st.expander("🎨 Customize this plot's appearance"):
        col1, col2 = st.columns(2)
        with col1:
            custom_title = st.text_input(
                "Custom title (optional):", value="", key=f"{key_prefix}_title",
                help="Leave blank to use the default title.",
            )
            font_family = st.selectbox(
                "Font style:", options=FONT_FAMILY_OPTIONS, key=f"{key_prefix}_font_family",
            )
            font_size = st.slider(
                "Font size:", min_value=8, max_value=28, value=13, key=f"{key_prefix}_font_size",
            )
        with col2:
            x_axis_label = st.text_input(
                "Custom X-axis label (optional):", value="", key=f"{key_prefix}_xlabel",
                help="Leave blank to use the default axis label.",
            )
            y_axis_label = st.text_input(
                "Custom Y-axis label (optional):", value="", key=f"{key_prefix}_ylabel",
                help="Leave blank to use the default axis label.",
            )
            theme_choice = st.selectbox(
                "Theme:", options=list(PLOTLY_THEME_OPTIONS.keys()), key=f"{key_prefix}_theme",
            )

        legend_position = "Right (default)"
        if show_legend_controls:
            st.markdown("**Legend placement**")
            legend_position = st.selectbox(
                "Starting legend position:",
                options=["Right (default)", "Top", "Bottom", "Left"],
                key=f"{key_prefix}_legend_pos",
            )
            st.caption(
                "💡 Tip: you can also click and drag the legend directly "
                "on the chart itself to fine-tune its position."
            )

        color_map = {}
        if group_values:
            st.markdown("**Colors**")
            n_cols = min(4, max(1, len(group_values)))
            color_cols = st.columns(n_cols)
            for i, group_val in enumerate(group_values):
                fallback_color = _DEFAULT_COLOR_PALETTE[i % len(_DEFAULT_COLOR_PALETTE)]
                default_color = default_colors.get(str(group_val), fallback_color)
                with color_cols[i % n_cols]:
                    color_map[str(group_val)] = st.color_picker(
                        str(group_val), value=default_color, key=f"{key_prefix}_color_{group_val}",
                    )

    return {
        "title": custom_title,
        "font_family": font_family,
        "font_size": font_size,
        "x_axis_label": x_axis_label,
        "y_axis_label": y_axis_label,
        "theme": PLOTLY_THEME_OPTIONS[theme_choice],
        "legend_position": legend_position,
        "color_map": color_map,
    }


def _apply_plot_style(fig, style, default_title="", default_x_label=None,
                       default_y_label=None):
    """
    Apply a style dict (from _render_plot_style_controls) to an
    existing Plotly figure: theme/template, title (falling back to
    default_title if the user left it blank), axis label overrides,
    font family/size, a starting legend placement, and per-trace
    recoloring (matched by trace name against style["color_map"]).

    The title is deliberately NOT set via Plotly's native
    fig.update_layout(title=...) -- Plotly.js has no built-in way to
    drag-reposition a native chart title. Instead, the title is added
    as a plain annotation (showarrow=False, positioned just above the
    plot area in paper coordinates). Combined with the
    edits.annotationPosition config flag set on every chart in this
    workspace (see _PLOTLY_CHART_CONFIG), this makes the title fully
    click-and-drag repositionable directly on the rendered chart, the
    same way the legend already is.

    No title_suffix parameter -- any "Before/After Batch Correction"
    -style default should be baked directly into default_title by the
    caller (e.g. default_title="PCA: Sample Similarity — Before Batch
    Correction"), so the user can freely retype/rename it via the
    style panel's title field with no forced, unremovable suffix.
    """
    title_text = (style.get("title") or "").strip() or default_title

    x_label = (style.get("x_axis_label") or "").strip() or default_x_label
    y_label = (style.get("y_axis_label") or "").strip() or default_y_label

    legend_layout = {}
    pos = style.get("legend_position", "Right (default)")
    if pos == "Top":
        legend_layout = dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5)
    elif pos == "Bottom":
        legend_layout = dict(orientation="h", yanchor="top", y=-0.2, xanchor="center", x=0.5)
    elif pos == "Left":
        legend_layout = dict(yanchor="middle", y=0.5, xanchor="right", x=-0.25)

    font_family = style.get("font_family", "Arial")
    font_size = style.get("font_size", 13)

    layout_kwargs = dict(
        template=style.get("theme", "plotly"),
        font=dict(family=font_family, size=font_size),
        # Native title cleared -- replaced by the draggable annotation
        # added below, so nothing renders twice.
        title="",
    )
    if title_text:
        # Reserve a bit of extra top margin so the title annotation
        # (placed above the plot area, at paper y=1.10) isn't clipped
        # -- otherwise a plot with tight default margins could cut it
        # off, especially at larger font sizes.
        layout_kwargs["margin"] = dict(t=max(70, int(font_size * 3.2)))
    if legend_layout:
        layout_kwargs["legend"] = legend_layout
    fig.update_layout(**layout_kwargs)

    if title_text:
        fig.add_annotation(
            text=title_text,
            xref="paper", yref="paper",
            x=0.5, y=1.10,
            xanchor="center", yanchor="bottom",
            showarrow=False,
            font=dict(family=font_family, size=font_size + 5),
            name="plot_title",
        )

    if x_label:
        fig.update_xaxes(title_text=x_label)
    if y_label:
        fig.update_yaxes(title_text=y_label)

    color_map = style.get("color_map") or {}
    if color_map:
        for trace in fig.data:
            trace_name = getattr(trace, "name", None)
            if trace_name in color_map:
                if hasattr(trace, "marker") and trace.marker is not None:
                    trace.marker.color = color_map[trace_name]
                if hasattr(trace, "line") and trace.line is not None and trace.line.color is not None:
                    trace.line.color = color_map[trace_name]

    return fig


def _sanitize_filename(custom_name, fallback, extension):
    """
    Sanitize a user-provided file name for use in a download_button:
    strip path separators (so a name like "../../etc/passwd" can't
    escape the intended download location), strip whitespace, and
    strip a trailing extension the user may have typed themselves so
    it isn't doubled up when the real extension is appended by the
    caller. Falls back to `fallback` if the result is empty.

    extension: the file extension WITHOUT a leading dot (e.g. "pdf",
    "csv") -- used only to strip a redundant user-typed extension;
    the caller is still responsible for appending "f.{extension}"
    themselves when building the final file_name.

    Returns the sanitized name, WITHOUT the extension.
    """
    safe = (custom_name or fallback).strip()
    safe = re.sub(r"[\\/]+", "_", safe)
    safe = re.sub(rf"\.{re.escape(extension)}$", "", safe, flags=re.IGNORECASE).strip()
    if not safe:
        safe = fallback
    return safe


def _render_csv_download(df, default_filename, key_prefix, expander_label,
                          help_text=None, button_label="⬇️ Download CSV"):
    """
    Render a "download this table as CSV" panel: a file name text
    input (pre-filled with a sensible default, but fully editable) and
    a download button using the (sanitized) chosen name -- mirroring
    the same expander + editable-filename pattern already used by
    _render_pdf_export, so every downloadable file in this workspace
    (PDF or CSV) offers the same rename-before-download experience
    rather than being locked to a fixed name.

    df: the DataFrame to export. Converted to CSV bytes fresh on every
        render -- cheap enough that no "Generate" button/step is needed
        here (unlike PDF export, which has real rendering cost via
        kaleido).
    expander_label: the visible label for the wrapping expander (e.g.
        the same text the download button used to show before this
        helper existed), so the download stays easy to find even
        though it's now tucked in an expander alongside the rename
        field.
    """
    with st.expander(expander_label):
        custom_name = st.text_input(
            "File name:", value=default_filename, key=f"{key_prefix}_csv_filename_input",
            help="The .csv extension is added automatically -- no need to type it.",
        )
        safe_name = _sanitize_filename(custom_name, default_filename, "csv")
        csv_bytes = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            button_label,
            data=csv_bytes,
            file_name=f"{safe_name}.csv",
            mime="text/csv",
            key=f"{key_prefix}_csv_download_btn",
            help=help_text,
        )


def _render_pdf_export(fig, key_prefix, filename_base):
    """
    Render a "save this plot as a PDF" panel for a given Plotly figure:
    a size preset selector (or custom width/height in pixels), a
    generate button, and (once generated) a download button for the
    resulting PDF bytes.

    Requires the `kaleido` package (Plotly's static image export
    backend) to be installed in the environment -- if it isn't, a
    clear message is shown instead of a broken button, matching the
    pattern used elsewhere in this app for optional external tools
    (e.g. Rscript, Salmon, STAR availability checks).
    """
    with st.expander("📄 Save this plot as a PDF"):
        if not _pdf_export_available():
            st.warning(
                "⚠️ PDF export requires the `kaleido` package, which "
                "isn't installed in this environment. Add `kaleido` to "
                "your project's Python requirements (and rebuild the "
                "Docker image) to enable this feature."
            )
            return

        size_choice = st.selectbox(
            "PDF size:", options=list(PDF_SIZE_PRESETS.keys()), key=f"{key_prefix}_pdf_size",
        )
        if size_choice == "Custom size":
            c1, c2 = st.columns(2)
            with c1:
                width_px = st.number_input(
                    "Width (pixels):", min_value=200, max_value=6000, value=1000, step=50,
                    key=f"{key_prefix}_pdf_width",
                )
            with c2:
                height_px = st.number_input(
                    "Height (pixels):", min_value=200, max_value=6000, value=700, step=50,
                    key=f"{key_prefix}_pdf_height",
                )
        else:
            width_px, height_px = PDF_SIZE_PRESETS[size_choice]
            st.caption(f"Output size: {width_px} x {height_px} pixels.")

        # User-editable file name, pre-filled with the sensible default
        # (filename_base) but fully overridable -- e.g. so a "Dexamethasone_vs_Untreated"
        # volcano plot could be saved as "figure_3_dex_volcano" instead.
        custom_filename = st.text_input(
            "File name:", value=filename_base, key=f"{key_prefix}_pdf_filename_input",
            help="The .pdf extension is added automatically -- no need to type it.",
        )

        if st.button("Generate PDF", key=f"{key_prefix}_pdf_generate_btn"):
            try:
                pdf_bytes = fig.to_image(format="pdf", width=width_px, height=height_px)
                st.session_state[f"{key_prefix}_pdf_bytes"] = pdf_bytes
            except Exception as e:
                st.session_state.pop(f"{key_prefix}_pdf_bytes", None)
                st.error(f"⚠️ Could not generate the PDF. Details: {e}")

        pdf_bytes = st.session_state.get(f"{key_prefix}_pdf_bytes")
        if pdf_bytes:
            safe_filename = _sanitize_filename(custom_filename, filename_base, "pdf")

            st.download_button(
                "⬇️ Download PDF",
                data=pdf_bytes,
                file_name=f"{safe_filename}.pdf",
                mime="application/pdf",
                key=f"{key_prefix}_pdf_download_btn",
            )


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_pca(pca_df, pct_variance, color_by):
    "Build a PC1 vs PC2 scatter plot, colored by the given metadata column."
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


def _plot_volcano(results_df, padj_threshold=0.05, lfc_threshold=1.0, gene_name_map=None):
    """
    Build a volcano plot (log2FoldChange vs. -log10(padj)), coloring
    points into three categories so up- and down-regulated genes are
    visually distinguishable at a glance rather than lumped into one
    "Significant" color:
      - "Up-regulated": significant AND log2FoldChange >= lfc_threshold
      - "Down-regulated": significant AND log2FoldChange <= -lfc_threshold
      - "Not significant": everything else
    ("Significant" here means padj < padj_threshold AND
    |log2FoldChange| >= lfc_threshold, same definition as before.)

    gene_name_map: optional dict {gene_id: human_readable_name} (see
        the "Gene ID -> Gene Name Mapping" uploader in render()). When
        given, hover text shows the readable name alongside the raw
        gene ID; when omitted, the raw gene ID is shown as-is.

    Returns (fig, df) -- the figure AND the prepared DataFrame (with
    the "neg_log10_padj" and "direction" columns added, NaN rows
    already dropped) so the caller can reuse it directly for building
    gene-label annotations (_build_gene_label_annotations) without
    recomputing the same classification a second time.

    Each trace also carries customdata=gene_id (as a plain, one-column
    array) so that Streamlit's st.plotly_chart(..., on_select="rerun")
    click-selection events can identify exactly which gene a clicked
    point corresponds to, regardless of which of the three
    up/down/not-significant traces it belongs to.
    """
    df = results_df.dropna(subset=["padj", "log2FoldChange"]).copy()
    # -log10(padj): clip padj away from exactly 0 first (log10(0) is
    # undefined) using a tiny floor value, so a gene with a vanishingly
    # small p-value still gets a large-but-finite y-position rather than
    # causing a math error.
    df["neg_log10_padj"] = -df["padj"].clip(lower=1e-300).apply(math.log10)

    # Shared classification logic (also used to annotate the
    # downloadable/displayed full results table below), so the
    # volcano's colors/legend and any exported "Regulation" column
    # always agree exactly. NaN rows are already dropped above (a gene
    # with no padj/log2FC has no valid plot position anyway), so
    # classify_regulation's "NA (not tested)" category never actually
    # appears here.
    df["direction"] = dm.classify_regulation(df, padj_threshold, lfc_threshold)

    gene_name_map = gene_name_map or {}
    df["display_name"] = df["gene_id"].astype(str).map(lambda g: gene_name_map.get(g, g))

    fig = go.Figure()
    for label, color in [
        ("Up-regulated", "#d62728"),
        ("Down-regulated", "#1f77b4"),
        ("Not significant", "#7f7f7f"),
    ]:
        subset = df[df["direction"] == label]
        fig.add_trace(go.Scatter(
            x=subset["log2FoldChange"], y=subset["neg_log10_padj"],
            mode="markers",
            name=label,
            marker=dict(size=6, color=color, opacity=0.6),
            text=subset["display_name"],
            customdata=subset["gene_id"],
            hovertemplate="%{text}<br>log2FC: %{x:.2f}<br>-log10(padj): %{y:.2f}<extra></extra>",
        ))

    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")
    # Horizontal significance threshold line, drawn at -log10(padj_threshold)
    # -- e.g. padj_threshold=0.05 -> line at y=1.301, matching the same
    # -log10 transform used for the plotted points above.
    padj_threshold_line_y = -math.log10(padj_threshold)
    fig.add_hline(y=padj_threshold_line_y, line_dash="dash", line_color="gray")

    fig.update_layout(
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(adjusted p-value)",
        height=550,
    )
    return fig, df


# Default style applied to a gene label the first time it's added
# (via "Label All" or by clicking its point) -- the user can then
# override any of these per-gene in the "Manage Gene Labels" panel.
_DEFAULT_GENE_LABEL_STYLE = {
    "font_family": "Arial",
    "font_size": 11,
    "bold": False,
    "color": "#000000",
}


def _build_gene_label_annotations(volcano_df, labeled_genes, gene_name_map=None):
    """
    Build one Plotly annotation per labeled gene: an arrow pointing
    from a small offset label position to that gene's actual data
    point (log2FoldChange, -log10(padj)) on the volcano plot, styled
    using that specific gene's own font/size/bold/color settings.

    volcano_df: the prepared DataFrame returned by _plot_volcano (must
        have "gene_id", "log2FoldChange", "neg_log10_padj" columns).
    labeled_genes: dict {gene_id: style_dict}, where style_dict has
        keys "font_family", "font_size", "bold", "color" (see
        _DEFAULT_GENE_LABEL_STYLE) -- this is the same dict rendered/
        edited by the "Manage Gene Labels" panel in render().
    gene_name_map: optional dict {gene_id: human_readable_name}; falls
        back to the raw gene_id if not provided or if a given gene_id
        isn't present in the mapping.

    Bold is applied via a literal "<b>...</b>" wrapper around the
    label text, since Plotly annotation "font" objects support family/
    size/color but have no separate boolean "weight"/"bold" property --
    wrapping text in Plotly's small supported pseudo-HTML subset is the
    documented way to achieve bold text in an annotation.

    Each annotation's tail (ax, ay -- the label's own position) is
    draggable independently of its arrowhead (which stays pinned to
    the gene's actual data point) via the "annotationTail" chart edit
    flag set in _PLOTLY_CHART_CONFIG.

    Returns a list of annotation dicts (NOT yet applied to any
    figure -- the caller merges this with any other annotations
    already on the figure, e.g. the plot title added by
    _apply_plot_style, via fig.update_layout(annotations=...)).
    """
    gene_name_map = gene_name_map or {}
    lookup = volcano_df.set_index(volcano_df["gene_id"].astype(str))

    annotations = []
    for gene_id, style in labeled_genes.items():
        gene_id_str = str(gene_id)
        if gene_id_str not in lookup.index:
            continue
        row = lookup.loc[gene_id_str]
        # A gene_id could theoretically appear more than once in the
        # results table (shouldn't happen with well-formed DESeq2
        # output, but .loc on a duplicate index returns a DataFrame
        # instead of a Series) -- defensively take the first row only.
        if hasattr(row, "iloc") and row.ndim == 2:
            row = row.iloc[0]

        display_name = gene_name_map.get(gene_id_str, gene_id_str)
        text = f"<b>{display_name}</b>" if style.get("bold") else display_name

        annotations.append(dict(
            x=row["log2FoldChange"], y=row["neg_log10_padj"],
            xref="x", yref="y",
            text=text,
            showarrow=True,
            arrowhead=2, arrowsize=1, arrowwidth=1,
            arrowcolor=style.get("color", "#000000"),
            ax=30, ay=-30,
            font=dict(
                family=style.get("font_family", "Arial"),
                size=style.get("font_size", 11),
                color=style.get("color", "#000000"),
            ),
            name=f"gene_label_{gene_id_str}",
        ))
    return annotations


def _render_gene_label_controls(contrast_key, volcano_df, gene_name_map=None):
    """
    Render the full gene-labeling control set for a volcano plot:
      - a "Label all significant genes" checkbox (adds every
        up/down-regulated gene to the label set; unchecking removes
        only the ones that were added this way, leaving any
        individually clicked genes alone)
      - an expander listing every currently labeled gene, each with
        its own font family, font size, bold checkbox, color picker,
        and a remove button

    Click-to-label itself is handled by the caller right after
    rendering the chart (since the click/selection event is only
    available as st.plotly_chart's return value, which requires the
    figure -- built using THIS function's returned labels dict -- to
    already exist first). See render()'s "Results by Contrast" section
    for exactly how the two halves (this function, and the post-render
    click handling) fit together.

    Returns (labeled_genes: dict, sig_gene_ids: set) -- the current
    {gene_id: style} dict (mutated in place as the user adjusts
    controls) and the set of gene_ids currently classified as
    Up-/Down-regulated (needed by the caller to implement "Label All").
    """
    labels_key = f"_volcano_gene_labels_{contrast_key}"
    auto_key = f"_volcano_gene_labels_auto_{contrast_key}"
    labeled_genes = st.session_state.setdefault(labels_key, {})
    auto_added = st.session_state.setdefault(auto_key, set())

    sig_gene_ids = set(
        volcano_df.loc[volcano_df["direction"].isin(["Up-regulated", "Down-regulated"]), "gene_id"].astype(str)
    )

    label_all = st.checkbox(
        "🏷️ Label all significant genes (up + down regulated)",
        key=f"label_all_{contrast_key}",
        help="Adds a name label to every up- or down-regulated point on the plot below. You can still remove or restyle individual labels afterward.",
    )

    if label_all:
        for gid in sig_gene_ids:
            if gid not in labeled_genes:
                labeled_genes[gid] = dict(_DEFAULT_GENE_LABEL_STYLE)
                auto_added.add(gid)
    else:
        # Only remove genes that were added BY "label all" (tracked via
        # auto_added) -- any gene the user individually clicked on stays
        # labeled even if this checkbox is unchecked afterward.
        for gid in list(auto_added):
            labeled_genes.pop(gid, None)
        auto_added.clear()

    if labeled_genes:
        gene_name_map = gene_name_map or {}
        with st.expander(f"🎨 Manage Gene Labels ({len(labeled_genes)} labeled)"):
            st.caption(
                "Adjust each labeled gene's text individually below, or "
                "remove it. Drag any label directly on the chart above "
                "to reposition it -- its arrow will keep pointing at "
                "the correct data point."
            )
            for gid in sorted(labeled_genes.keys(), key=lambda g: gene_name_map.get(g, g)):
                style = labeled_genes[gid]
                display_name = gene_name_map.get(gid, gid)
                cols = st.columns([2.2, 1.6, 1.2, 0.9, 1, 0.6])
                with cols[0]:
                    st.markdown(f"**{display_name}**")
                with cols[1]:
                    current_family = style.get("font_family", "Arial")
                    family_index = FONT_FAMILY_OPTIONS.index(current_family) if current_family in FONT_FAMILY_OPTIONS else 0
                    style["font_family"] = st.selectbox(
                        "Font", options=FONT_FAMILY_OPTIONS, index=family_index,
                        key=f"gene_label_font_{contrast_key}_{gid}", label_visibility="collapsed",
                    )
                with cols[2]:
                    style["font_size"] = st.slider(
                        "Size", min_value=6, max_value=28, value=style.get("font_size", 11),
                        key=f"gene_label_size_{contrast_key}_{gid}", label_visibility="collapsed",
                    )
                with cols[3]:
                    style["bold"] = st.checkbox(
                        "Bold", value=style.get("bold", False),
                        key=f"gene_label_bold_{contrast_key}_{gid}",
                    )
                with cols[4]:
                    style["color"] = st.color_picker(
                        "Color", value=style.get("color", "#000000"),
                        key=f"gene_label_color_{contrast_key}_{gid}", label_visibility="collapsed",
                    )
                with cols[5]:
                    if st.button("✖", key=f"gene_label_remove_{contrast_key}_{gid}", help=f"Remove label for {display_name}"):
                        labeled_genes.pop(gid, None)
                        auto_added.discard(gid)
                        st.rerun()

    return labeled_genes, sig_gene_ids


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
    required) -- good enough for a quick visual overlap check, though not
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

def _render_gene_id_mapping_panel(project, counts_df):
    """
    Renders the "Gene ID -> Gene Name Mapping" expander and returns the
    resulting gene_name_map dict to use for the volcano plot's hover
    text/gene labels for the rest of render().

    Three layers, each overriding the one before it:
      1. The fast, no-R-required mapping auto-built in
         alignment_workspace.py (Ensembl FASTA/GTF parsing) -- this is
         what's on disk at pm.gene_symbol_map_path(project) by default.
         Genes it couldn't resolve fall back to their raw reference ID
         (see reference_manager.py's extract_gene_symbol_map_from_*).
      2. An optional Bioconductor clusterProfiler::bitr() conversion
         (see gene_id_mapper.py), run on-demand from this panel, either
         to "fill in" just the still-unresolved genes from layer 1, or
         to fully re-convert every gene ID to a different target
         namespace entirely (e.g. Entrez ID instead of Symbol). Both
         modes persist their result back to the same
         gene_symbol_map.csv on disk, so they survive reopening the
         project later.
      3. A manually uploaded CSV, which only lives in this session
         (st.session_state) and always wins over layers 1-2 -- useful
         for a mapping sourced from somewhere entirely outside this
         app's reference/annotation pipeline.
    """
    auto_map_path = pm.gene_symbol_map_path(project)
    auto_gene_name_map = {}
    if os.path.exists(auto_map_path):
        try:
            auto_map_df = pd.read_csv(auto_map_path)
            auto_gene_name_map = dict(
                zip(auto_map_df["gene_id"].astype(str), auto_map_df["gene_name"].astype(str))
            )
        except Exception:
            auto_gene_name_map = {}

    mapping_meta = pm.get_gene_id_mapping_meta(project) or {}
    gene_ids_all = sorted(set(counts_df["gene_id"].astype(str)))
    n_total = len(gene_ids_all)
    # A gene counts as "still unresolved" if the mapping has no entry
    # for it, or its entry is just the gene_id echoed back at itself
    # (the fallback convention used by every mapping source in this
    # app -- see reference_manager.py and gene_id_mapper.py).
    unresolved_ids = [gid for gid in gene_ids_all if auto_gene_name_map.get(gid, gid) == gid]
    n_resolved = n_total - len(unresolved_ids)

    with st.expander("🏷️ Gene ID → Gene Name Mapping", expanded=False):
        if auto_gene_name_map and n_resolved > 0:
            pct = (n_resolved / n_total * 100) if n_total else 0.0
            source_detail = mapping_meta.get("detail") or mapping_meta.get("source", "")
            st.success(
                f"✅ **{n_resolved:,} of {n_total:,} gene IDs ({pct:.1f}%)** "
                f"have a readable name{f' ({source_detail})' if source_detail else ''}. "
                f"**{len(unresolved_ids):,}** gene(s) are still shown by "
                "their raw reference ID below."
            )
        else:
            st.info(
                "ℹ️ No gene names could be resolved automatically from "
                "this project's reference yet -- gene IDs are shown as-is "
                "for now."
            )

        # --- Bioconductor bitr() conversion ---
        st.markdown("##### 🧬 Convert Gene IDs (Bioconductor `bitr`)")

        detection_pool = unresolved_ids if unresolved_ids else gene_ids_all
        detection = gim.detect_id_type(detection_pool)
        detected_type = detection["detected_type"]

        species_key, is_custom = pm.get_reference_choice(project)
        species_choices = gim.orgdb_species_choices(rm)
        species_keys = list(species_choices.keys())
        default_species_index = (
            species_keys.index(species_key) if species_key in species_keys else 0
        )

        st.caption(
            f"Auto-detected ID type for your {'unresolved ' if unresolved_ids else ''}"
            f"gene IDs: **{detected_type}** ({detection['match_fraction'] * 100:.0f}% "
            f"of a sample matched -- e.g. `{'`, `'.join(detection['example_ids'][:3])}`). "
            "Adjust the options below if this doesn't look right, or if "
            "you'd like a different output format entirely (e.g. Entrez "
            "ID instead of Symbol)."
        )

        col_species, col_from, col_to = st.columns(3)
        with col_species:
            picked_species = st.selectbox(
                "Species (annotation database):",
                options=species_keys,
                index=default_species_index,
                format_func=lambda k: species_choices[k],
                key="gene_id_map_species_select",
            )
        with col_from:
            from_type_options = gim.COMMON_KEY_TYPES
            from_default_idx = from_type_options.index(detected_type) if detected_type in from_type_options else 0
            picked_from_type = st.selectbox(
                "My gene IDs are currently in this format:",
                options=from_type_options,
                index=from_default_idx,
                key="gene_id_map_from_type_select",
            )
        with col_to:
            to_type_options = gim.COMMON_KEY_TYPES
            default_to_type = gim.symbol_keytype_for_species(picked_species)
            to_default_idx = to_type_options.index(default_to_type) if default_to_type in to_type_options else 0
            picked_to_type = st.selectbox(
                "Convert them to:",
                options=to_type_options,
                index=to_default_idx,
                key="gene_id_map_to_type_select",
            )

        orgdb_package = gim.ORGDB_PACKAGES.get(picked_species)
        work_dir = pm.gene_id_mapping_work_dir(project)

        col_fill, col_override = st.columns(2)
        with col_fill:
            fill_disabled = not unresolved_ids
            if st.button(
                f"✨ Fill in {len(unresolved_ids):,} missing name(s) only",
                key="bitr_fill_missing_btn",
                disabled=fill_disabled,
                help="Recommended: only converts the genes that don't already have a readable name, leaving everything else untouched." if not fill_disabled else "Every gene already has a readable name -- nothing to fill in.",
            ):
                with st.spinner(f"Converting {len(unresolved_ids):,} gene ID(s) via bitr()..."):
                    result = gim.run_bitr_conversion(
                        unresolved_ids, picked_from_type, picked_to_type, orgdb_package, work_dir,
                    )
                if result["success"]:
                    # Only overwrite entries that were previously
                    # unresolved -- genes that already had a name from
                    # the fast auto-parse are left exactly as they were.
                    auto_gene_name_map.update(result["mapping"])
                    rm.save_gene_symbol_map_csv(auto_gene_name_map, auto_map_path)
                    pm.save_gene_id_mapping_meta(project, {
                        "source": "bitr_fill",
                        "detail": f"bitr: {picked_from_type} → {picked_to_type} ({species_choices[picked_species]})",
                        "from_type": picked_from_type,
                        "to_type": picked_to_type,
                        "orgdb_package": orgdb_package,
                        "n_converted": result["n_converted"],
                        "n_total": result["n_total"],
                    })
                    st.success(f"✅ {result['message']}")
                else:
                    st.error(f"⚠️ {result['message']}")

        with col_override:
            if st.button(
                "🔁 Convert ALL gene IDs (override existing mapping)",
                key="bitr_override_all_btn",
                help="Re-converts every gene ID from scratch to your chosen output format, replacing the current mapping entirely -- use this if you want a different ID type than what's currently shown (e.g. Entrez instead of Symbol).",
            ):
                with st.spinner(f"Converting {n_total:,} gene ID(s) via bitr()..."):
                    result = gim.run_bitr_conversion(
                        gene_ids_all, picked_from_type, picked_to_type, orgdb_package, work_dir,
                    )
                if result["success"]:
                    auto_gene_name_map = result["mapping"]
                    rm.save_gene_symbol_map_csv(auto_gene_name_map, auto_map_path)
                    pm.save_gene_id_mapping_meta(project, {
                        "source": "bitr_full_override",
                        "detail": f"bitr: {picked_from_type} → {picked_to_type} ({species_choices[picked_species]})",
                        "from_type": picked_from_type,
                        "to_type": picked_to_type,
                        "orgdb_package": orgdb_package,
                        "n_converted": result["n_converted"],
                        "n_total": result["n_total"],
                    })
                    st.success(f"✅ {result['message']}")
                else:
                    st.error(f"⚠️ {result['message']}")

        if not gim.bitr_tools_available():
            st.caption(
                "ℹ️ Rscript wasn't found on this system -- the buttons "
                "above need R with the `clusterProfiler` package and the "
                "relevant Bioconductor annotation package (e.g. "
                f"`{orgdb_package}` for {species_choices.get(picked_species, picked_species)}) "
                "installed."
            )

        st.markdown("---")
        st.caption(
            "Alternatively, upload your own 2-column CSV (columns named "
            "exactly 'gene_id' and 'gene_name') -- this stays local to "
            "your current session and takes priority over everything "
            "above."
        )
        gene_map_file = st.file_uploader(
            "Gene ID → Gene Name CSV (session-only override):",
            type=["csv"], key="gene_name_mapping_upload",
        )
        if gene_map_file is not None:
            try:
                map_df = pd.read_csv(gene_map_file)
                if "gene_id" in map_df.columns and "gene_name" in map_df.columns:
                    st.session_state["_gene_name_map"] = dict(
                        zip(map_df["gene_id"].astype(str), map_df["gene_name"].astype(str))
                    )
                    st.success(f"✅ Loaded {len(st.session_state['_gene_name_map']):,} gene ID → name mapping(s) (overriding everything above, for this session only).")
                else:
                    st.error("⚠️ This CSV must have columns named exactly 'gene_id' and 'gene_name'.")
            except Exception as e:
                st.error(f"⚠️ Could not read this file: {e}")

    # A manually uploaded mapping (this session) takes precedence over
    # the auto/bitr-derived one; otherwise use whatever's currently in
    # auto_gene_name_map (freshly updated above if a conversion just
    # ran), and finally an empty dict (raw gene IDs shown as-is).
    return st.session_state.get("_gene_name_map") or auto_gene_name_map


def _render_missing_value_resolution(meta_df, columns_to_check):
    """
    For each column in columns_to_check, detect ambiguous values (true
    blanks/NaN, or missing-like text like "none"/"N/A") and -- if any are
    found -- show a resolution UI letting the user assign each one a real
    category label, or mark affected samples for exclusion.

    Returns (resolved_meta_df, all_resolved: bool). If any ambiguous
    value in any checked column hasn't been explicitly resolved yet,
    all_resolved is False and the caller should block proceeding to the
    next step -- this prevents a real NaN or an unresolved "none" from
    silently reaching DESeq2, where it would either error out or (worse)
    be misinterpreted.
    """
    any_ambiguous = dm.has_any_ambiguous_values(meta_df, columns_to_check)
    if not any_ambiguous:
        return meta_df, True

    st.warning(
        "⚠️ We found some blank or ambiguous values in your selected "
        "column(s) below. Please tell us what each one means before "
        "continuing -- for example, a blank cell or \"none\" in an "
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
                f"-- affected sample(s): {', '.join(affected_samples)}"
            )
            blank_choice = st.radio(
                f"What does a blank value in `{column}` mean?",
                ["This represents a real category (e.g. no treatment/control)", "These samples should be excluded from analysis"],
                key=f"blank_choice_{column}",
                horizontal=False,
            )
            if blank_choice.startswith("This represents"):
                label = st.text_input(
                    "What should this category be called?",
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
                f"-- affected sample(s): {', '.join(affected_samples)}"
            )
            value_choice = st.radio(
                f"What does \"{ambiguous_val}\" mean in `{column}`?",
                ["This represents a real category (e.g. no treatment/control)", "These samples should be excluded from analysis"],
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

def _render_analysis_settings(project, deseq2_out_dir, deseq2_work_dir,
                               counts_matrix_path, counts_df, meta_df,
                               existing_results, saved_config):
    """
    Renders Steps 1-4 (low-count filtering, experimental design,
    analysis type, and running DESeq2) inside the collapsible
    "Analysis Settings" expander.

    IMPORTANT: this is a separate function (rather than inlined
    directly in render()) specifically so that its many internal
    early `return` statements (guarding against unselected design
    columns, unresolved ambiguous metadata values, under-replicated
    groups, Rscript not being available, etc.) only exit THIS
    function -- not the caller. Previously, when this logic lived
    directly in render() and one of these checks failed while the
    expander was auto-collapsed (which happens exactly when a
    project already has existing results -- see
    expanded=not(existing_results and saved_config) below), the
    resulting warning/error was hidden inside the collapsed
    expander AND the early `return` prevented "Step 5: Explore
    Your Results" from ever rendering below it -- so a returning
    user would see the "this project already has DESeq2 results"
    success message immediately above, then nothing, with no
    visible explanation. Extracting this into its own function
    means Step 5 in render() always gets a chance to render
    (gated only by whether result files actually exist on disk),
    regardless of whether re-validating Steps 1-4 on this
    particular rerun happens to succeed.
    """
    with st.expander(
        "⚙️ Steps 1-4: Analysis Settings (filtering, design, contrasts, run)",
        expanded=not (existing_results and saved_config),
    ):
        # -----------------------------------------------------------------
        # STEP 1: Low-count gene filtering
        # -----------------------------------------------------------------
        st.header("Step 1: Remove Low-Count Genes")

        with st.expander("ℹ️ Why remove low-count genes? (click to learn more)"):
            st.markdown(
                "Genes with very few reads across your samples don't carry "
                "enough information to reliably detect real differences -- "
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
                "select multiple columns -- DESeq2 will account for all of "
                "them together.\n\n"
                "**Batch effects** happen when samples processed at "
                "different times, by different people, or on different "
                "equipment show systematic differences unrelated to your "
                "actual biological question. If you have a column "
                "identifying batches (e.g. `batch`, `sequencing_run`, "
                "`prep_date`), select it below -- DESeq2 will account for "
                "batch variation *while* estimating your condition effects, "
                "which is the statistically correct way to handle this "
                "(rather than 'correcting' the data directly, which can "
                "discard real information).\n\n"
                "**Important -- what this selection does and doesn't do:** "
                "your project's full metadata file (every column, every "
                "sample) stays exactly as it is -- this step does **not** "
                "remove or overwrite anything. You're only telling DESeq2 "
                "*for this specific analysis run* which column(s) to treat "
                "as the experimental condition(s). If you later want to "
                "test a **different** condition (e.g. you tested "
                "`antibiotic` this time, but later want to test `media` "
                "instead), or test conditions in a **different combination** "
                "(e.g. `antibiotic` alone, then `antibiotic` + `media` "
                "together), you'll come back to this exact step and simply "
                "select different column(s) -- each combination you want to "
                "test is its own pass through Steps 2-4, and each will save "
                "its own separate results without affecting the others."
            )

        available_columns = [c for c in meta_df.columns if c != "sample"]

        st.caption(
            "💡 Tip: if your metadata came from an auto-fill (e.g. SRA "
            "lookup), it often includes extra descriptive columns (like a "
            "sample title, strain, or organism name) that **aren't** "
            "meant to be experimental conditions -- only select the "
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
        # capitalization, blanks) before committing -- rather than having to
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
        #      in a single-species study) -- DESeq2 hard-errors on this.
        #   2. A column that's essentially a per-sample identifier (as many
        #      unique values as samples, e.g. a free-text sample
        #      description) -- technically won't crash, but is never a
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
                "used as a design/batch factor -- DESeq2 requires some "
                "variation across samples to estimate an effect. Please "
                "deselect it above."
            )
            return

        if identifier_columns:
            st.warning(
                f"⚠️ **{', '.join(identifier_columns)}** has a different "
                "value for almost every sample (like a free-text "
                "description or ID), which usually means it's not a real "
                "experimental condition -- including it will likely make "
                "your model meaningless. If this wasn't intentional, "
                "deselect it above."
            )

        # --- Interaction terms (optional, only relevant with 2+ design columns) ---
        # Chosen here, before the replication check below, so that check
        # knows exactly which column pairs (if any) actually need
        # cross-replication -- an additive model needs far less replication
        # than one with an explicit interaction term.
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
                    "biological reason to expect one -- adding unnecessary "
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

        # --- Check for sufficient replication before proceeding ---
        # DESeq2 needs at least 2 samples per group to estimate dispersion
        # (within-group variability) at all -- with only 1 sample per group,
        # there's no data left to distinguish real biological signal from
        # noise, and DESeq2 hard-errors during dispersion estimation
        # ("checkForExperimentalReplicates"). This is a very common issue
        # when metadata comes from a small, hand-picked subset of samples
        # (e.g. a handful of runs pulled from a larger public dataset) where
        # each one happens to represent a different, unreplicated condition.
        #
        # Note this checks each main-effect column (design/batch) for >= 2
        # samples per level *individually*, and only cross-checks specific
        # column pairs when an interaction term between them was actually
        # selected above -- NOT every combination of every selected column.
        # A standard additive model (e.g. ~ batch + condition) does not
        # require 2+ samples in every batch:condition cell, so crossing
        # everything unconditionally would incorrectly reject valid,
        # standard experimental designs (e.g. a randomized-block layout).
        replication = dm.check_replication(meta_df, design_columns, batch_column, interaction_terms)
        if not replication["is_valid"]:
            st.error(
                f"⚠️ The following group(s) only have 1 sample: "
                f"**{', '.join(replication['under_replicated_groups'])}**. "
                "DESeq2 needs **at least 2 samples per group** to estimate "
                "natural variability and distinguish real differences from "
                "noise -- with only 1 sample, this isn't statistically "
                "possible.\n\n"
                "To fix this, you'll need to either add more samples for "
                "the under-replicated group(s) (e.g. download additional "
                "runs from the same study/condition), choose a different "
                "condition/combination of conditions where every group has "
                "at least 2 samples, or remove the interaction term causing "
                "the under-replicated combination."
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

        st.markdown("---")

        # -----------------------------------------------------------------
        # STEP 3: Choose Your Analysis Type
        # -----------------------------------------------------------------
        st.header("Step 3: Choose Your Analysis Type")

        with st.expander("ℹ️ Wald test vs. LRT (Likelihood Ratio Test) -- which should I use? (click to learn more)"):
            st.markdown(
                "DESeq2 offers two different statistical approaches, and "
                "which one you want depends on your question:\n\n"
                "**🎯 Wald test (pairwise comparisons)** -- tests one "
                "specific comparison at a time, e.g. \"is drugA different "
                "from control?\". This is what most people want when they "
                "have a clear reference group and want to know how each "
                "other group differs from it. You can run several pairwise "
                "comparisons (contrasts) from one model fit.\n\n"
                "**📐 LRT (Likelihood Ratio Test)** -- DESeq2's equivalent of "
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
                    "interpreted cautiously -- the **p-value/adjusted "
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
        # it out to the same metadata_path the R script will load -- this
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


def render():
    st.title("📈 Differential Expression (DESeq2)")
    st.markdown(
        "This workspace compares gene expression between your "
        "experimental conditions to find genes that are significantly "
        "up- or down-regulated. **No statistics background required** -- "
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
    # --- Load previously saved DESeq2 configuration, if any ---
    # Lets a returning user land back on their existing results without
    # re-selecting the exact same design/contrasts in Steps 1-3. Values
    # are pushed into each widget's session_state *once* per project per
    # session (guarded by config_applied_key below) -- Streamlit uses
    # session_state as a widget's initial value whenever the widget's
    # key already has one, which is the standard, documented way to
    # pre-fill a widget programmatically. After this first application,
    # the user's own edits take over normally on every later rerun, the
    # same as any other widget.
    saved_config = pm.get_deseq2_config(project)
    existing_results = dm.list_contrast_results(deseq2_out_dir)
    config_applied_key = f"_deseq2_config_applied_{project}"

    if saved_config and not st.session_state.get(config_applied_key):
        st.session_state["min_count_input"] = saved_config.get("min_count", 10)
        st.session_state["min_samples_input"] = saved_config.get("min_samples", min(2, len(meta_df)))
        st.session_state["design_columns_select"] = saved_config.get("design_columns", [])
        st.session_state["batch_column_select"] = saved_config.get("batch_column") or "(none)"
        st.session_state["interaction_terms_select"] = saved_config.get("interaction_terms", [])
        st.session_state["test_type_radio"] = (
            "Wald test (pairwise comparisons)"
            if saved_config.get("test_type", "wald") == "wald"
            else "LRT (omnibus / ANOVA-style test)"
        )

        saved_contrasts = saved_config.get("contrasts") or []
        if saved_contrasts:
            st.session_state["primary_contrast_column"] = saved_contrasts[0].get("column")
            st.session_state["reference_level_select"] = saved_contrasts[0].get("level2")
            st.session_state["comparison_levels_select"] = [c.get("level1") for c in saved_contrasts]

        saved_lrt = saved_config.get("lrt_tests") or []
        if saved_lrt:
            # Recover the original tested term(s) per LRT test by diffing
            # the full formula's terms against that test's reduced_terms
            # -- more robust than trying to reverse the "_x_"-joined test
            # name back into a column name, which would be lossy if a
            # column name itself happened to contain "_x_".
            full_terms_saved = dm.build_full_formula_terms(
                saved_config.get("design_columns", []),
                saved_config.get("batch_column"),
                saved_config.get("interaction_terms", []),
            )
            recovered_terms = []
            for lrt_test in saved_lrt:
                reduced = lrt_test.get("reduced_terms", [])
                dropped = [t for t in full_terms_saved if t not in reduced]
                if dropped:
                    recovered_terms.append(dropped[0])
            if recovered_terms:
                st.session_state["lrt_terms_select"] = recovered_terms

        st.session_state[config_applied_key] = True

    if existing_results and saved_config:
        st.success(
            f"✅ This project already has DESeq2 results — "
            f"{dm.summarize_deseq2_config(saved_config)}. Jump straight to "
            "**Step 5** below to explore them, or expand **Steps 1-4** "
            "below if you'd like to change the analysis settings and "
            "re-run."
        )

    st.markdown("---")

    _render_analysis_settings(
        project, deseq2_out_dir, deseq2_work_dir, counts_matrix_path,
        counts_df, meta_df, existing_results, saved_config,
    )

    # Re-read the saved config after _render_analysis_settings runs, since
    # it may have just been created/updated by a fresh DESeq2 run above --
    # this is how Step 5 below gets its own copy of batch_column (needed
    # for the before/after batch-correction PCA views) without that
    # variable having to leak out of _render_analysis_settings' local
    # scope. Falls back to an empty dict for a first-time run where saving
    # happens inside _render_analysis_settings but nothing was on disk
    # before this call.
    saved_config = pm.get_deseq2_config(project) or {}
    batch_column = saved_config.get("batch_column")

    st.markdown("---")

    # -----------------------------------------------------------------
    # STEP 5: Explore Results
    # -----------------------------------------------------------------
    st.header("Step 5: Explore Your Results")

    available_contrasts = dm.list_contrast_results(deseq2_out_dir)
    if not available_contrasts:
        st.info("No results found yet -- run DESeq2 above.")
        return

    # --- PCA ---
    st.subheader("🔬 PCA: Sample Similarity")
    with st.expander("ℹ️ How to read this plot"):
        st.markdown(
            "Each point is one sample. Samples that are more similar in "
            "overall gene expression appear closer together. If your "
            "samples cluster clearly by condition, that's a good sign -- "
            "if they instead cluster by batch, that suggests a strong "
            "batch effect (which DESeq2 has already accounted for "
            "statistically in your results, via the design formula)."
        )

    pca_df, pct_var = dm.read_pca_coordinates(deseq2_out_dir)
    if pca_df is not None:
        color_options = [c for c in pca_df.columns if c not in ("PC1", "PC2", "PC3", "PC4", "sample")]
        color_by = st.selectbox("Color points by:", options=color_options, key="pca_color_select")
        group_values = sorted(pca_df[color_by].astype(str).unique())

        pca_adj_df, pct_var_adj = (None, None)
        if batch_column:
            pca_adj_df, pct_var_adj = dm.read_pca_coordinates(deseq2_out_dir, batch_adjusted=True)

        if pca_adj_df is not None:
            st.info(
                "🔀 Two views are shown below: **Before** batch correction "
                "(the raw variance-stabilized data) and **After** batch "
                "correction (for visualization only -- limma's "
                "`removeBatchEffect` applied purely so you can *see* how "
                "much of the sample separation was attributable to "
                "batch). The actual DESeq2 statistical test already "
                "accounts for batch correctly via the design formula, "
                "regardless of what's shown here.\n\n"
                "💡 Each plot below has its own independent style panel "
                "-- including its title -- so you can fully rename "
                "either one (e.g. drop \"(Visualization Only)\") if you "
                "prefer."
            )
            col_before, col_after = st.columns(2)
            with col_before:
                st.markdown("#### 🔹 Before Batch Correction")
                fig_before = _plot_pca(pca_df, pct_var, color_by)
                pca_before_style = _render_plot_style_controls("pca_before", group_values=group_values)
                _apply_plot_style(
                    fig_before, pca_before_style,
                    default_title="PCA: Sample Similarity — Before Batch Correction",
                )
                _render_plotly_chart(fig_before)
                _render_pdf_export(fig_before, "pca_before", "pca_before_batch_correction")
                _render_csv_download(
                    pca_df, "pca_coordinates_before_batch_correction", "pca_before",
                    expander_label="⬇️ Download PCA Coordinates (Before).csv",
                    help_text="Every sample's PC1-PC4 coordinates plus every metadata column available for coloring. Percent variance explained and batch-effect size are in the QC summary download below.",
                )
            with col_after:
                st.markdown("#### 🔸 After Batch Correction")
                fig_after = _plot_pca(pca_adj_df, pct_var_adj, color_by)
                pca_after_style = _render_plot_style_controls("pca_after", group_values=group_values)
                _apply_plot_style(
                    fig_after, pca_after_style,
                    default_title="PCA: Sample Similarity — After Batch Correction (Visualization Only)",
                )
                _render_plotly_chart(fig_after)
                _render_pdf_export(fig_after, "pca_after", "pca_after_batch_correction")
                _render_csv_download(
                    pca_adj_df, "pca_coordinates_after_batch_correction", "pca_after",
                    expander_label="⬇️ Download PCA Coordinates (After).csv",
                    help_text="Every sample's batch-adjusted PC1-PC4 coordinates (visualization only) plus every metadata column available for coloring.",
                )

            # --- Additional QC insights from the batch correction ---
            st.markdown("##### 📐 Batch-Correction QC Insights")

            variance_rows = []
            n_pcs = min(len(pct_var or []), len(pct_var_adj or []), 4)
            for i in range(n_pcs):
                variance_rows.append({
                    "Principal Component": f"PC{i + 1}",
                    "% Variance (Before)": pct_var[i],
                    "% Variance (After)": pct_var_adj[i],
                    "Change": round(pct_var_adj[i] - pct_var[i], 1),
                })
            if variance_rows:
                st.caption(
                    "Variance explained per principal component, before "
                    "vs. after the batch-adjustment view:"
                )
                st.dataframe(pd.DataFrame(variance_rows), use_container_width=True, hide_index=True)

            batch_eta_before = dm.compute_batch_variance_metric(pca_df, batch_column)
            batch_eta_after = dm.compute_batch_variance_metric(pca_adj_df, batch_column)
            eta_rows = []
            if batch_eta_before and batch_eta_after:
                eta_rows = [
                    {
                        "Principal Component": pc,
                        "% of Variance Explained by Batch (Before)": batch_eta_before.get(pc, 0.0),
                        "% of Variance Explained by Batch (After)": batch_eta_after.get(pc, 0.0),
                    }
                    for pc in ("PC1", "PC2")
                    if pc in batch_eta_before
                ]
                st.caption(
                    "Estimated share of each PC's variance attributable "
                    f"to your batch column (`{batch_column}`) -- a simple "
                    "between-group/total variance ratio, not a formal "
                    "statistical test, but a useful before/after sanity "
                    "check on how much batch-driven separation the "
                    "visualization adjustment removed:"
                )
                st.dataframe(pd.DataFrame(eta_rows), use_container_width=True, hide_index=True)

            pca_qc_export_df = dm.build_pca_qc_export(variance_rows, eta_rows)
            if not pca_qc_export_df.empty:
                _render_csv_download(
                    pca_qc_export_df, "pca_qc_summary_batch_correction", "pca_qc_summary",
                    expander_label="⬇️ Download PCA QC Summary (variance + batch effect size).csv",
                )
        else:
            st.markdown("#### PCA: Sample Similarity")
            fig_single = _plot_pca(pca_df, pct_var, color_by)
            pca_single_style = _render_plot_style_controls("pca_single", group_values=group_values)
            _apply_plot_style(fig_single, pca_single_style, default_title="PCA: Sample Similarity")
            _render_plotly_chart(fig_single)
            _render_pdf_export(fig_single, "pca_single", "pca_sample_similarity")

            single_qc_rows = [
                {"Principal Component": f"PC{i + 1}", "% Variance": pct_var[i]}
                for i in range(min(len(pct_var or []), 4))
            ]
            col_coords, col_variance = st.columns(2)
            with col_coords:
                _render_csv_download(
                    pca_df, "pca_coordinates", "pca_single",
                    expander_label="⬇️ Download PCA Coordinates.csv",
                    help_text="Every sample's PC1-PC4 coordinates plus every metadata column available for coloring.",
                )
            with col_variance:
                if single_qc_rows:
                    _render_csv_download(
                        pd.DataFrame(single_qc_rows), "pca_percent_variance", "pca_single_variance",
                        expander_label="⬇️ Download % Variance Explained.csv",
                    )

    st.markdown("---")

    # --- Per-contrast results + volcano ---
    st.subheader("📊 Results by Contrast")

    gene_name_map = _render_gene_id_mapping_panel(project, counts_df)

    selected_contrast = st.selectbox("Select a contrast to view:", options=available_contrasts, key="view_contrast_select")

    results_df = dm.read_contrast_results(deseq2_out_dir, selected_contrast)
    if results_df is not None:
        padj_cutoff = st.slider("Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key="padj_slider")
        lfc_cutoff = st.slider("Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key="lfc_slider")

        n_sig = ((results_df["padj"] < padj_cutoff) & (results_df["log2FoldChange"].abs() >= lfc_cutoff)).sum()
        st.caption(f"**{n_sig:,}** significant genes at padj < {padj_cutoff} and |log2FC| >= {lfc_cutoff}.")

        volcano_fig, volcano_df = _plot_volcano(results_df, padj_cutoff, lfc_cutoff, gene_name_map=gene_name_map)
        volcano_style = _render_plot_style_controls(
            f"volcano_{selected_contrast}",
            group_values=["Up-regulated", "Down-regulated", "Not significant"],
            default_colors={
                "Up-regulated": "#d62728",
                "Down-regulated": "#1f77b4",
                "Not significant": "#7f7f7f",
            },
        )
        _apply_plot_style(
            volcano_fig, volcano_style,
            default_title=f"Volcano Plot: {selected_contrast}",
        )

        # --- Gene labeling: "Label All" + per-gene style controls ---
        # Built (and its resulting annotations attached to volcano_fig)
        # BEFORE the chart is rendered below, so labels added on a
        # previous run/rerun show up immediately. Newly CLICKED genes
        # are handled just after rendering (see the on_select handling
        # below it), since a click's result is only available as
        # st.plotly_chart's return value -- that discovery triggers one
        # extra st.rerun() so the clicked gene's label appears using
        # the same code path as everything else, keeping this logic in
        # one place rather than duplicated.
        labeled_genes, sig_gene_ids = _render_gene_label_controls(
            selected_contrast, volcano_df, gene_name_map=gene_name_map,
        )
        gene_annotations = _build_gene_label_annotations(volcano_df, labeled_genes, gene_name_map=gene_name_map)
        if gene_annotations:
            volcano_fig.update_layout(annotations=list(volcano_fig.layout.annotations) + gene_annotations)

        st.caption(
            "💡 Tip: click any point on the chart below to add a label "
            "for that gene. Drag any label to reposition it -- its "
            "arrow will keep pointing at the correct dot."
        )
        selection_event = st.plotly_chart(
            volcano_fig, use_container_width=True, config=_PLOTLY_CHART_CONFIG,
            on_select="rerun", selection_mode=["points"],
            key=f"volcano_chart_{selected_contrast}",
        )

        # Process any newly clicked point(s): each trace's customdata
        # was set to gene_id in _plot_volcano, so a click's gene
        # identity survives regardless of which of the three
        # up/down/not-significant traces the point belongs to.
        newly_added = False
        if selection_event and selection_event.selection and selection_event.selection.points:
            for pt in selection_event.selection.points:
                gid = pt.get("customdata")
                if isinstance(gid, (list, tuple)):
                    gid = gid[0] if gid else None
                if gid is None:
                    continue
                gid = str(gid)
                if gid not in labeled_genes:
                    labeled_genes[gid] = dict(_DEFAULT_GENE_LABEL_STYLE)
                    newly_added = True
        if newly_added:
            st.rerun()

        _render_pdf_export(volcano_fig, f"volcano_{selected_contrast}", f"volcano_{selected_contrast}")

        # Annotate a copy of the full (unfiltered, every gene) results
        # table with a "Regulation" column at the current threshold
        # settings, using the exact same classification the volcano
        # plot uses above -- so both the on-screen preview and the
        # downloadable CSV always show the complete DESeq2 output for
        # every gene (not just the significant/up/down subset), with
        # the up/down/not-significant/not-tested status included as an
        # extra, clearly-labeled column rather than filtering rows out.
        annotated_results_df = results_df.copy()
        annotated_results_df["Regulation"] = dm.classify_regulation(annotated_results_df, padj_cutoff, lfc_cutoff)

        st.dataframe(annotated_results_df.head(200), use_container_width=True, hide_index=True)
        st.caption(
            "Showing top 200 rows (sorted by adjusted p-value). The "
            "**Regulation** column reflects the significance threshold "
            "and fold-change cutoff sliders above. Download the full, "
            "unfiltered table (all genes, every DESeq2 statistic, plus "
            "this Regulation column) below."
        )

        _render_csv_download(
            annotated_results_df, f"deseq2_results_{selected_contrast}", f"results_{selected_contrast}",
            expander_label=f"⬇️ Download Full Results ({selected_contrast}).csv",
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
                "or `biomaRt` -- this depends on your organism and isn't "
                "something we can determine automatically."
            )
            export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrast)
            if export_df is not None:
                st.dataframe(export_df.head(20), use_container_width=True, hide_index=True)
                _render_csv_download(
                    export_df, f"clusterprofiler_ranked_{selected_contrast}", f"clusterprofiler_{selected_contrast}",
                    expander_label="⬇️ Download clusterProfiler-Ready Gene List (.csv)",
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

            with st.expander("✏️ Rename contrast labels (optional)"):
                st.caption(
                    "Customize how each contrast is labeled on the "
                    "diagram below -- e.g. shorten "
                    "\"Albuterol_Dexamethasone_vs_Untreated\" to \"Alb + "
                    "Dex\". This only changes the display label; the "
                    "underlying contrast and its results are unaffected. "
                    "You can also drag each label directly on the chart "
                    "once it's drawn."
                )
                venn_labels = {}
                for name in venn_contrasts:
                    venn_labels[name] = st.text_input(
                        f"Label for `{name}`:", value=name, key=f"venn_label_{name}",
                    )

            gene_sets = dm.get_significant_gene_sets(deseq2_out_dir, venn_contrasts, padj_threshold=venn_padj)

            # Remap to the user's custom display labels (falling back to
            # the original contrast name if left blank) *before*
            # computing regions, so every downstream label -- the Venn
            # circles/annotations and the "View overlapping gene lists"
            # expander below -- shows the custom name consistently.
            # Guards against two contrasts accidentally being renamed to
            # the *same* label, which would otherwise silently merge two
            # different gene sets together under one key.
            seen_labels = set()
            display_gene_sets = {}
            orig_to_display_label = {}
            for name, gene_set in gene_sets.items():
                custom_label = (venn_labels.get(name, "") or "").strip() or name
                if custom_label in seen_labels:
                    st.warning(
                        f"⚠️ The label \"{custom_label}\" is used more "
                        f"than once -- keeping the original name "
                        f"\"{name}\" for that contrast instead, so the "
                        "two aren't combined."
                    )
                    custom_label = name
                seen_labels.add(custom_label)
                display_gene_sets[custom_label] = gene_set
                orig_to_display_label[name] = custom_label

            display_names = list(display_gene_sets.keys())
            counts, regions = dm.compute_venn_regions(display_gene_sets)

            venn_fig = _plot_venn(counts, display_names)
            venn_style = _render_plot_style_controls("venn", group_values=None, show_legend_controls=False)
            _apply_plot_style(venn_fig, venn_style, default_title="Overlap Between Contrasts")
            st.caption(
                "💡 Tip: click and drag any label directly on the chart "
                "below (circle names, region counts, or the title) to "
                "reposition it."
            )
            _render_plotly_chart(venn_fig)
            _render_pdf_export(venn_fig, "venn", "venn_diagram_overlap")

            with st.expander("View overlapping gene lists"):
                for label, gene_set in regions.items():
                    st.markdown(f"**{label}** ({len(gene_set)} genes)")
                    if gene_set:
                        st.text(", ".join(sorted(gene_set)[:50]) + (" ..." if len(gene_set) > 50 else ""))

            # Full-detail CSV export: every gene shown on the diagram,
            # which Venn region it's in, AND its actual DESeq2 stats
            # (log2FoldChange, padj) from every selected contrast --
            # not just the truncated, name-only list shown above.
            venn_export_df = dm.build_venn_export(
                deseq2_out_dir, venn_contrasts, regions,
                display_label_map=orig_to_display_label,
            )
            if not venn_export_df.empty:
                _render_csv_download(
                    venn_export_df, "venn_overlap_full_results", "venn_full",
                    expander_label="⬇️ Download Full Venn Results (all genes + stats per contrast).csv",
                    help_text="One row per gene shown on the diagram, with its Venn region and log2FoldChange/padj from every selected contrast -- not just gene names.",
                )

    st.markdown("---")
    st.success(
        f"🎉 Project `{project}` has completed differential expression "
        "analysis. Results, PCA, volcano plots, and Venn diagrams are "
        "available above, along with clusterProfiler-ready exports for "
        "downstream gene ontology enrichment analysis."
    )

    if st.button("➡️ Proceed to Ontology Analysis", type="primary", key="de_proceed_ontology_btn"):
        # Same nav_request indirection used by the other "Proceed to X"
        # buttons in this app (see bulk_rnaseq_workspace.py and app.py's
        # module docstring for why a plain session key is used here
        # instead of directly setting st.session_state["assay_choice_radio"]).
        st.session_state["nav_request"] = "🧬 Ontology Analysis"
        st.rerun()
