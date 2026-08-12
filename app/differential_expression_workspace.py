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
# Cached disk reads
# ---------------------------------------------------------------------------
#
# Streamlit reruns this entire script top-to-bottom on EVERY widget
# interaction (moving a slider, checking a box, typing a character) --
# this is fundamental to how Streamlit works, not something specific to
# this file. Without caching, that means every one of this workspace's
# CSV reads (counts matrix, metadata, PCA coordinates, sample distance
# matrix, dispersion estimates, size factors, per-contrast results,
# normalized counts) gets re-read from disk and re-parsed by pandas on
# every single interaction, even ones that only affect a completely
# unrelated widget elsewhere on the page (e.g. moving the volcano
# plot's font-size slider still re-reads the sample distance matrix).
# For a small project this is unnoticeable; for a larger one (many
# thousands of genes and/or many samples) these repeated reads/parses
# are the most likely source of any perceptible lag -- see
# https://docs.streamlit.io/develop/concepts/architecture/caching.
#
# Each wrapper below is keyed on the underlying file's modification
# time (mtime) IN ADDITION to its normal arguments, so Streamlit's
# cache is used on every interaction EXCEPT when the file has actually
# changed on disk (e.g. after re-running DESeq2 overwrites it) --
# there's no need to manually clear the cache or worry about serving
# stale results after a re-run. If a file doesn't exist yet, mtime is
# None, which is still a perfectly valid (if unchanging) cache key.
def _mtime_or_none(path):
    "Get a file's modification time for use as an automatic cache-invalidation key, or None if it doesn't exist."
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


#
# IMPORTANT: none of the cache-key parameters below (mtime, coords_mtime,
# var_mtime) may start with a leading underscore. Streamlit's
# @st.cache_data specifically treats a leading-underscore parameter name
# as "do NOT hash this / do NOT use it as part of the cache key" (that
# convention exists so unhashable objects, like a DB connection, can
# still be passed through a cached function safely). Since the entire
# point of passing an mtime here is to have it INCLUDED in the cache
# key so a changed file correctly invalidates the cache, these
# parameters must be named WITHOUT a leading underscore.
@st.cache_data(show_spinner=False)
def _cached_read_csv(path, mtime):
    "Cached pd.read_csv, invalidated automatically whenever the file's mtime changes."
    return pd.read_csv(path)


@st.cache_data(show_spinner=False)
def _cached_read_pca_coordinates(deseq2_out_dir, batch_adjusted, coords_mtime, var_mtime):
    return dm.read_pca_coordinates(deseq2_out_dir, batch_adjusted=batch_adjusted)


@st.cache_data(show_spinner=False)
def _cached_read_sample_distance_matrix(deseq2_out_dir, mtime):
    return dm.read_sample_distance_matrix(deseq2_out_dir)


@st.cache_data(show_spinner=False)
def _cached_read_dispersion_estimates(deseq2_out_dir, mtime):
    return dm.read_dispersion_estimates(deseq2_out_dir)


@st.cache_data(show_spinner=False)
def _cached_read_size_factors(deseq2_out_dir, mtime):
    return dm.read_size_factors(deseq2_out_dir)


@st.cache_data(show_spinner=False)
def _cached_read_contrast_results(deseq2_out_dir, contrast_name, mtime):
    return dm.read_contrast_results(deseq2_out_dir, contrast_name)


@st.cache_data(show_spinner=False)
def _cached_read_normalized_counts(deseq2_out_dir, mtime):
    return dm.read_normalized_counts(deseq2_out_dir)


def _load_counts_df(counts_matrix_path):
    "Cached read of the project's gene counts matrix CSV."
    return _cached_read_csv(counts_matrix_path, _mtime_or_none(counts_matrix_path))


def _load_meta_df(metadata_path):
    "Cached read of the project's metadata CSV."
    return _cached_read_csv(metadata_path, _mtime_or_none(metadata_path))


def _load_pca_coordinates(deseq2_out_dir, batch_adjusted=False):
    "Cached wrapper for dm.read_pca_coordinates, auto-invalidated when the underlying CSV/txt files change."
    suffix = "_batch_adjusted" if batch_adjusted else ""
    coords_path = os.path.join(deseq2_out_dir, f"pca_coordinates{suffix}.csv")
    var_path = os.path.join(deseq2_out_dir, f"pca_percent_variance{suffix}.txt")
    return _cached_read_pca_coordinates(
        deseq2_out_dir, batch_adjusted, _mtime_or_none(coords_path), _mtime_or_none(var_path),
    )


def _load_sample_distance_matrix(deseq2_out_dir):
    "Cached wrapper for dm.read_sample_distance_matrix, auto-invalidated when the underlying CSV changes."
    path = os.path.join(deseq2_out_dir, "sample_distance_matrix.csv")
    return _cached_read_sample_distance_matrix(deseq2_out_dir, _mtime_or_none(path))


def _load_dispersion_estimates(deseq2_out_dir):
    "Cached wrapper for dm.read_dispersion_estimates, auto-invalidated when the underlying CSV changes."
    path = os.path.join(deseq2_out_dir, "dispersion_estimates.csv")
    return _cached_read_dispersion_estimates(deseq2_out_dir, _mtime_or_none(path))


def _load_size_factors(deseq2_out_dir):
    "Cached wrapper for dm.read_size_factors, auto-invalidated when the underlying CSV changes."
    path = os.path.join(deseq2_out_dir, "size_factors.csv")
    return _cached_read_size_factors(deseq2_out_dir, _mtime_or_none(path))


def _load_contrast_results(deseq2_out_dir, contrast_name):
    "Cached wrapper for dm.read_contrast_results, auto-invalidated when that contrast's CSV changes."
    path = os.path.join(deseq2_out_dir, f"results_{contrast_name}.csv")
    return _cached_read_contrast_results(deseq2_out_dir, contrast_name, _mtime_or_none(path))


def _load_normalized_counts(deseq2_out_dir):
    "Cached wrapper for dm.read_normalized_counts, auto-invalidated when the underlying CSV changes."
    path = os.path.join(deseq2_out_dir, "normalized_counts.csv")
    return _cached_read_normalized_counts(deseq2_out_dir, _mtime_or_none(path))

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
                options=["Right (default)", "Top", "Bottom", "Left", "Bottom-right (inside plot)"],
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
    elif pos == "Bottom-right (inside plot)":
        # Anchored INSIDE the plot area's bottom-right corner, rather
        # than off to the right of the whole figure (Plotly's default)
        # -- useful for plots that already have a colorbar occupying
        # the space just outside the right edge (e.g. the heatmaps in
        # this workspace), where the default right-side legend would
        # otherwise render directly beneath/overlapping the colorbar.
        legend_layout = dict(
            yanchor="bottom", y=0.01, xanchor="right", x=0.99,
            bgcolor="rgba(255,255,255,0.7)", bordercolor="rgba(0,0,0,0.2)", borderwidth=1,
        )

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


_CONFIDENCE_LEVEL_OPTIONS = {"90%": 0.90, "95%": 0.95, "99%": 0.99}


def _render_confidence_ellipse_controls(key_prefix):
    """
    Render a "Show confidence ellipses around groups" checkbox and (if
    checked) a confidence-level selector, for a PCA plot.

    Returns (show_ellipse: bool, confidence_level: float) -- pass these
    directly as _plot_pca's show_confidence_ellipse/confidence_level
    arguments.
    """
    col_toggle, col_level = st.columns([2, 1])
    with col_toggle:
        show_ellipse = st.checkbox(
            "🔵 Show confidence ellipses around groups",
            value=False,
            key=f"{key_prefix}_show_ellipse",
            help="Draws a dashed ellipse around each group's samples, sized to the selected confidence level (assumes each group's samples are roughly bivariate-normal on PC1/PC2). Groups with fewer than 3 samples are skipped.",
        )
    confidence_level = 0.95
    if show_ellipse:
        with col_level:
            level_choice = st.selectbox(
                "Confidence level:",
                options=list(_CONFIDENCE_LEVEL_OPTIONS.keys()),
                index=1,
                key=f"{key_prefix}_ellipse_confidence",
                label_visibility="collapsed",
            )
            confidence_level = _CONFIDENCE_LEVEL_OPTIONS[level_choice]
    return show_ellipse, confidence_level


def _render_ellipse_stats(ellipse_stats, key_prefix):
    """
    Display a small stats table for the confidence ellipses drawn on a
    PCA plot (see _compute_confidence_ellipse) -- one row per group,
    with its sample count, confidence level, semi-major/semi-minor axis
    lengths, rotation angle, and enclosed area -- plus a CSV download of
    the same table. Renders nothing if ellipse_stats is empty (e.g. the
    user hasn't enabled ellipses, or every group had too few samples).
    """
    if not ellipse_stats:
        return
    rows = [
        {
            "Group": s["group"],
            "N Samples": s["n_samples"],
            "Confidence Level": f"{int(s['confidence'] * 100)}%",
            "Semi-Major Axis": round(s["semi_major_axis"], 3),
            "Semi-Minor Axis": round(s["semi_minor_axis"], 3),
            "Angle (degrees)": round(s["angle_degrees"], 1),
            "Area": round(s["area"], 3),
        }
        for s in ellipse_stats
    ]
    stats_df = pd.DataFrame(rows)
    st.caption(
        "📐 Confidence ellipse statistics per group (PC1/PC2 units -- "
        "semi-major/semi-minor axis lengths, rotation angle, and "
        "enclosed area of each dashed ellipse above):"
    )
    st.dataframe(stats_df, use_container_width=True, hide_index=True)
    _render_csv_download(
        stats_df, f"{key_prefix}_confidence_ellipse_stats", f"{key_prefix}_ellipse_stats",
        expander_label="⬇️ Download Confidence Ellipse Statistics (.csv)",
    )


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

# Chi-square critical values for 2 degrees of freedom, at the confidence
# levels offered in the PCA plot's "Show confidence ellipses" control --
# this is the standard way to size a 2D confidence ellipse from a
# covariance matrix (the ellipse boundary is the set of points at
# Mahalanobis distance sqrt(chi2_critical_value) from the group's
# centroid). Hardcoded here (rather than computed via scipy.stats.chi2)
# so this feature works even in an environment without scipy installed,
# matching this workspace's existing "graceful without scipy" approach
# for hierarchical clustering.
_CHI2_CRITICAL_VALUES_2DOF = {
    0.90: 4.605,
    0.95: 5.991,
    0.99: 9.210,
}


def _compute_confidence_ellipse(x_vals, y_vals, confidence=0.95, n_points=100):
    """
    Compute the boundary points of a 2D confidence ellipse for a set of
    (x, y) points, assuming an underlying bivariate normal distribution
    -- the standard "confidence ellipse" construction used across
    statistical plotting tools (e.g. matplotlib/seaborn's
    confidence_ellipse recipes, R's stat_ellipse): the ellipse is
    centered at the group's centroid, oriented along the eigenvectors
    of its covariance matrix, with semi-axis lengths scaled by the
    chi-square critical value for 2 degrees of freedom at the requested
    confidence level (see _CHI2_CRITICAL_VALUES_2DOF).

    Requires at least 3 points to compute a meaningful (non-degenerate)
    covariance matrix -- returns None for smaller groups rather than
    drawing a misleading ellipse from too little data.

    Returns a dict:
        {
            "ellipse_x": [...], "ellipse_y": [...],  # boundary points, for plotting
            "center_x": float, "center_y": float,      # centroid (also the ellipse's own center)
            "semi_major_axis": float, "semi_minor_axis": float,
            "angle_degrees": float,   # rotation of the semi-major axis from the x-axis
            "area": float,
            "n_samples": int,
            "confidence": float,
        }
    or None if there aren't enough points.
    """
    import numpy as np

    x_arr = np.asarray(x_vals, dtype=float)
    y_arr = np.asarray(y_vals, dtype=float)
    n = len(x_arr)
    if n < 3:
        return None

    center_x, center_y = x_arr.mean(), y_arr.mean()
    cov = np.cov(x_arr, y_arr)

    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # np.linalg.eigh returns eigenvalues in ASCENDING order -- reverse so
    # index 0 is the largest (major axis), matching the semi_major/
    # semi_minor naming below.
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Guard against a degenerate (zero-variance) direction -- e.g. every
    # sample in this group happens to share an identical PC value on one
    # axis -- which would otherwise produce a zero-width ellipse or a
    # math domain error taking its square root.
    eigenvalues = np.clip(eigenvalues, a_min=1e-10, a_max=None)

    chi2_val = _CHI2_CRITICAL_VALUES_2DOF.get(confidence, 5.991)
    semi_major = float(np.sqrt(eigenvalues[0] * chi2_val))
    semi_minor = float(np.sqrt(eigenvalues[1] * chi2_val))
    angle_rad = float(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))

    theta = np.linspace(0, 2 * np.pi, n_points)
    ellipse_local_x = semi_major * np.cos(theta)
    ellipse_local_y = semi_minor * np.sin(theta)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    ellipse_x = center_x + ellipse_local_x * cos_a - ellipse_local_y * sin_a
    ellipse_y = center_y + ellipse_local_x * sin_a + ellipse_local_y * cos_a

    return {
        "ellipse_x": ellipse_x.tolist(),
        "ellipse_y": ellipse_y.tolist(),
        "center_x": center_x,
        "center_y": center_y,
        "semi_major_axis": semi_major,
        "semi_minor_axis": semi_minor,
        "angle_degrees": float(np.degrees(angle_rad)),
        "area": float(math.pi * semi_major * semi_minor),
        "n_samples": n,
        "confidence": confidence,
    }


def _plot_pca(pca_df, pct_variance, color_by, show_confidence_ellipse=False, confidence_level=0.95):
    """
    Build a PC1 vs PC2 scatter plot, colored by the given metadata
    column.

    show_confidence_ellipse, confidence_level: if show_confidence_ellipse
        is True, draws a dashed confidence ellipse (see
        _compute_confidence_ellipse) around each group's samples,
        colored to match that group's markers (using this workspace's
        default color palette, matched by index in the same order the
        marker traces are built -- if the user later recolors a group
        via the style panel, _apply_plot_style's color_map recoloring
        picks up the ellipse trace too, since it shares its group's
        exact trace "name"). Groups with fewer than 3 samples are
        skipped (too little data for a meaningful ellipse) rather than
        causing an error.

    Returns (fig, ellipse_stats) -- the figure, and a list of per-group
    ellipse statistics dicts (see _compute_confidence_ellipse) for
    display below the plot, or an empty list if
    show_confidence_ellipse is False or no group had enough samples.
    """
    pc1_label = f"PC1 ({pct_variance[0]}% variance)" if pct_variance else "PC1"
    pc2_label = f"PC2 ({pct_variance[1]}% variance)" if pct_variance else "PC2"

    fig = go.Figure()
    ellipse_stats = []
    group_labels = sorted(pca_df[color_by].astype(str).unique())
    for i, group_value in enumerate(group_labels):
        subset = pca_df[pca_df[color_by].astype(str) == group_value]
        fig.add_trace(go.Scatter(
            x=subset["PC1"], y=subset["PC2"],
            mode="markers+text",
            text=subset["sample"],
            textposition="top center",
            name=group_value,
            marker=dict(size=12),
        ))

        if show_confidence_ellipse:
            ellipse = _compute_confidence_ellipse(subset["PC1"], subset["PC2"], confidence=confidence_level)
            if ellipse is not None:
                ellipse_color = _DEFAULT_COLOR_PALETTE[i % len(_DEFAULT_COLOR_PALETTE)]
                fig.add_trace(go.Scatter(
                    x=ellipse["ellipse_x"], y=ellipse["ellipse_y"],
                    mode="lines",
                    line=dict(color=ellipse_color, dash="dot", width=2),
                    name=group_value,
                    showlegend=False,
                    hoverinfo="skip",
                ))
                ellipse_stats.append(dict(ellipse, group=group_value))

    fig.update_layout(
        xaxis_title=pc1_label,
        yaxis_title=pc2_label,
        height=500,
        legend_title=color_by,
    )
    return fig, ellipse_stats


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
      - a search box to find and label ANY gene by its ID or
        human-readable symbol, regardless of where it falls on the
        plot or whether it's significant (see _render_gene_search)
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

    _render_gene_search(contrast_key, volcano_df, labeled_genes, gene_name_map=gene_name_map)

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


def _render_gene_search(contrast_key, volcano_df, labeled_genes, gene_name_map=None):
    """
    Render a "search for a gene" box that lets the user find and label
    ANY gene on the volcano plot by typing either its reference ID
    (e.g. "ENSG00000141510", version suffix optional) or its
    human-readable symbol (e.g. "TP53", if a gene name mapping is
    active -- see _render_gene_id_mapping_panel), regardless of where
    that gene falls on the plot (deep in "Not significant", visually
    overlapping other points, off in a crowded cluster, etc.) or
    whether it happens to be significant.

    Matching is substring-based and case-insensitive, checked against
    both the raw gene_id AND its mapped display name (if any), so a
    partial symbol (e.g. "tp53") or a partial/unversioned ID still
    finds the right gene. If more than one gene matches (e.g. an
    ambiguous partial symbol shared by a gene family), every match is
    shown in a follow-up picker so the user can disambiguate rather
    than a guess being made on their behalf.

    Adds a match directly into `labeled_genes` (the same dict used by
    click-to-label and "Label All" -- see _render_gene_label_controls),
    so a gene found this way gets identical treatment: it shows up in
    the "Manage Gene Labels" panel, can be restyled or removed there,
    and its label can be dragged on the chart like any other.
    """
    gene_name_map = gene_name_map or {}

    st.markdown("**🔎 Find a specific gene**")
    search_col, button_col = st.columns([3, 1])
    with search_col:
        search_term = st.text_input(
            "Search by gene ID or symbol:",
            key=f"gene_search_input_{contrast_key}",
            placeholder="e.g. TP53 or ENSG00000141510",
            label_visibility="collapsed",
            help="Finds and labels a gene anywhere on the plot, even if it's not significant or is hard to click precisely.",
        )
    with button_col:
        search_clicked = st.button("Find", key=f"gene_search_btn_{contrast_key}", use_container_width=True)

    if not search_term or not search_term.strip():
        return

    term = search_term.strip().lower()
    all_gene_ids = volcano_df["gene_id"].astype(str)

    # Match against BOTH the raw gene_id and its mapped display name
    # (if a gene name mapping is active), so the user can search by
    # whichever one they happen to know -- e.g. typing "tp53" finds a
    # gene whose gene_id is "ENSG00000141510" via its mapped symbol,
    # and typing an (optionally versioned) Ensembl ID finds it directly
    # even if a mapping never resolved a name for it.
    def _matches(gid):
        gid_str = str(gid)
        if term in gid_str.lower():
            return True
        display_name = gene_name_map.get(gid_str, "")
        return bool(display_name) and term in display_name.lower()

    matching_ids = [gid for gid in all_gene_ids if _matches(gid)]
    # Preserve first-seen order but drop duplicates (a gene_id should
    # already be unique in volcano_df, but this is a cheap safeguard).
    seen = set()
    matching_ids = [gid for gid in matching_ids if not (gid in seen or seen.add(gid))]

    if not matching_ids:
        st.warning(f"⚠️ No gene found matching \"{search_term}\". Check the spelling, or try just part of the ID/symbol.")
        return

    if len(matching_ids) == 1:
        gid_to_add = matching_ids[0]
        display_name = gene_name_map.get(gid_to_add, gid_to_add)
        already_labeled = gid_to_add in labeled_genes
        st.caption(f"✅ Found **{display_name}** (`{gid_to_add}`)" + (" -- already labeled." if already_labeled else "."))
        if search_clicked and not already_labeled:
            labeled_genes[gid_to_add] = dict(_DEFAULT_GENE_LABEL_STYLE)
            st.success(f"🏷️ Labeled **{display_name}** on the plot below.")
            st.rerun()
        elif search_clicked and already_labeled:
            st.info(f"**{display_name}** is already labeled -- adjust it in \"Manage Gene Labels\" below.")
    else:
        st.caption(f"Found **{len(matching_ids)}** genes matching \"{search_term}\" -- pick one:")
        options_map = {
            f"{gene_name_map.get(gid, gid)} ({gid})": gid for gid in matching_ids
        }
        picked_label = st.selectbox(
            "Select the gene you meant:",
            options=list(options_map.keys()),
            key=f"gene_search_disambiguate_{contrast_key}",
            label_visibility="collapsed",
        )
        gid_to_add = options_map[picked_label]
        already_labeled = gid_to_add in labeled_genes
        if st.button(
            "🏷️ Label this gene" if not already_labeled else "Already labeled",
            key=f"gene_search_confirm_{contrast_key}",
            disabled=already_labeled,
        ):
            labeled_genes[gid_to_add] = dict(_DEFAULT_GENE_LABEL_STYLE)
            st.success(f"🏷️ Labeled **{gene_name_map.get(gid_to_add, gid_to_add)}** on the plot below.")
            st.rerun()


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
# QC diagnostic plots
# ---------------------------------------------------------------------------
#
# These plots surface data DESeq2 already computes internally as part of
# a normal run (dispersion estimates, size factors, sample distances,
# p-values, mean expression) -- see deseq2_manager.py's module docstring
# ("QC/diagnostics additions") for what the R script now exports and why.
# Each one pairs with a plain-language classification/flag from
# deseq2_manager.py (assess_dispersion_fit, classify_pvalue_histogram_shape,
# classify_ma_bias, flag_size_factor_outliers, detect_sample_clustering_mismatch)
# so the user gets an informed read on the plot, not just the raw chart.

def _cluster_order_from_distance_matrix(distance_df):
    """
    Compute a hierarchical-clustering leaf order for a square sample
    distance matrix, for use as the row/column order of the sample
    distance heatmap -- grouping similar samples adjacently makes any
    clustering pattern (or lack thereof) far easier to read at a glance
    than an arbitrary/alphabetical ordering.

    Uses scipy's hierarchical clustering (average linkage) if
    available; gracefully falls back to the distance matrix's existing
    (alphabetical) sample order if scipy isn't installed, so this
    optional enhancement never blocks the heatmap itself from
    rendering.

    Returns a list of sample names in clustered order.
    """
    samples = list(distance_df.index)
    if len(samples) < 3:
        return samples
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
        # squareform requires a symmetric matrix with an exact-zero
        # diagonal -- guard against tiny floating-point asymmetry from
        # R's dist()/CSV round-trip by symmetrizing explicitly first.
        matrix = distance_df.values
        matrix = (matrix + matrix.T) / 2
        import numpy as np
        np.fill_diagonal(matrix, 0.0)
        condensed = squareform(matrix, checks=False)
        Z = linkage(condensed, method="average")
        order = leaves_list(Z)
        return [samples[i] for i in order]
    except Exception:
        return samples


def _plot_sample_distance_heatmap(distance_df, meta_df=None, group_column=None, show_group_annotation=True):
    """
    Build a sample-to-sample distance heatmap (Euclidean distance on
    VST-transformed counts), with rows/columns ordered by hierarchical
    clustering -- the standard companion QC view to PCA, since PCA's
    first 2-4 components can sometimes miss an outlier or mislabeled
    sample that a full pairwise distance comparison catches. Darker
    (lower-distance) cells indicate more similar samples.

    meta_df, group_column: if given (a metadata DataFrame with a
        "sample" column, and the name of the column to group by), each
        axis tick label is augmented with that sample's group value
        (e.g. "sample1 (control)") and a slim color-coded annotation
        strip is drawn directly above the heatmap, one colored cell per
        sample, so which samples belong to which group is visible at a
        glance from color alone -- without the user needing to already
        recognize every sample name. If either argument is omitted,
        falls back to plain sample-name-only labels (the original
        behavior), so this enhancement is fully optional.
    show_group_annotation: set False to suppress the color-coded strip,
        legend, and label augmentation entirely -- even if meta_df/
        group_column are provided -- falling back to plain sample-name
        labels, same as if those arguments were omitted. This is
        separate from whether group_column itself is set, so a caller
        (see render()'s "Hide grouping bar and legend" checkbox) can
        keep group_column selected for OTHER purposes (e.g. the
        clustering-mismatch smart flag below the heatmap) while still
        hiding the visual annotation on the plot itself.
    """
    order = _cluster_order_from_distance_matrix(distance_df)
    ordered_df = distance_df.loc[order, order]

    sample_to_group = {}
    if show_group_annotation and meta_df is not None and group_column and group_column in meta_df.columns:
        sample_to_group = dict(
            zip(meta_df["sample"].astype(str), meta_df[group_column].astype(str))
        )

    # Build display labels that include each sample's group (when
    # available), e.g. "sample1 (control)" -- shown on both axes and in
    # the hover text, so groupings are readable even without the color
    # strip (colorblind-safe fallback, and useful in the PDF export
    # where hover text isn't available).
    if sample_to_group:
        display_labels = [
            f"{s} ({sample_to_group.get(s, '?')})" for s in ordered_df.index
        ]
    else:
        display_labels = list(ordered_df.index)

    fig = go.Figure(data=go.Heatmap(
        z=ordered_df.values,
        x=display_labels,
        y=display_labels,
        colorscale="Blues_r",
        colorbar=dict(title="Distance", x=1.02),
        hovertemplate="%{y} vs %{x}<br>Distance: %{z:.1f}<extra></extra>",
    ))

    if sample_to_group:
        # --- Color-coded group annotation strip ---
        # A slim second heatmap trace, one column wide, positioned via
        # its own x-axis (xaxis2) directly to the right of the main
        # heatmap -- one colored cell per sample (in the same clustered
        # order), so which samples share a group is visible purely from
        # color, at a glance, without reading every label. This mirrors
        # the standard "annotation sidebar" convention used by
        # heatmap tools like pheatmap/ComplexHeatmap for exactly this
        # purpose.
        groups_in_order = [sample_to_group.get(s, "?") for s in ordered_df.index]
        unique_groups = sorted(set(groups_in_order))
        group_to_code = {g: i for i, g in enumerate(unique_groups)}
        z_strip = [[group_to_code[g]] for g in groups_in_order]

        n_groups = max(len(unique_groups), 1)
        strip_colors = [
            _DEFAULT_COLOR_PALETTE[i % len(_DEFAULT_COLOR_PALETTE)] for i in range(n_groups)
        ]
        # Build a discrete colorscale (Plotly heatmaps only support
        # continuous colorscales natively) by giving each integer group
        # code a hard-edged color band rather than a smooth gradient.
        if n_groups == 1:
            colorscale = [[0, strip_colors[0]], [1, strip_colors[0]]]
        else:
            colorscale = []
            for i, color in enumerate(strip_colors):
                colorscale.append([i / n_groups, color])
                colorscale.append([(i + 1) / n_groups, color])

        fig.add_trace(go.Heatmap(
            z=z_strip,
            x=["Group"],
            y=display_labels,
            xaxis="x2",
            colorscale=colorscale,
            zmin=0, zmax=n_groups,
            showscale=False,
            text=[[g] for g in groups_in_order],
            hovertemplate="%{y}<br>Group: %{text}<extra></extra>",
        ))

        fig.update_layout(
            xaxis=dict(domain=[0, 0.94], tickangle=-45),
            xaxis2=dict(domain=[0.96, 1.0], side="top"),
        )

        # Manual legend (since the annotation strip's own colorbar is
        # hidden above via showscale=False -- a discrete legend reads
        # far more clearly here than a continuous color axis for a
        # handful of category labels) -- one invisible-marker scatter
        # trace per group, purely so its name/color shows up in the
        # figure's standard legend. Positioned INSIDE the plot area's
        # bottom-right corner by default (rather than Plotly's default
        # right-of-figure placement), since the main heatmap's colorbar
        # already occupies that space just outside the right edge --
        # stacking a right-side legend directly under/beside it there
        # looked crowded. The user's own style panel choice (see
        # _apply_plot_style's "legend_position" handling) still
        # overrides this if they pick something else.
        fig.update_layout(legend=dict(
            yanchor="bottom", y=0.01, xanchor="right", x=0.99,
            bgcolor="rgba(255,255,255,0.7)", bordercolor="rgba(0,0,0,0.2)", borderwidth=1,
        ))
        for group_name, color in zip(unique_groups, strip_colors):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=10, color=color),
                name=str(group_name),
                showlegend=True,
            ))
    else:
        fig.update_layout(xaxis=dict(tickangle=-45))

    fig.update_layout(height=550)
    return fig


def _plot_dispersion(dispersion_df):
    """
    Build DESeq2's standard dispersion diagnostic plot (equivalent to
    plotDispEsts in R): each gene's raw gene-wise dispersion estimate
    (light gray points) plotted against its mean expression, the fitted
    mean-dispersion trend line DESeq2 estimated across all genes (solid
    line), and the final shrunken dispersion value actually used in the
    statistical test (colored points) -- so the user can see both how
    much shrinkage was applied and whether the overall trend looks like
    a reasonable fit for this dataset.

    Both axes are log-scaled, matching the standard convention for this
    plot (dispersion and mean expression both span several orders of
    magnitude).
    """
    df = dispersion_df.dropna(subset=["baseMean"]).copy()
    df = df[df["baseMean"] > 0]

    fig = go.Figure()

    gene_est = df.dropna(subset=["dispGeneEst"])
    gene_est = gene_est[gene_est["dispGeneEst"] > 0]
    fig.add_trace(go.Scatter(
        x=gene_est["baseMean"], y=gene_est["dispGeneEst"],
        mode="markers", name="Gene-wise estimate",
        marker=dict(size=4, color="#B0B0B0", opacity=0.5),
        text=gene_est["gene_id"],
        hovertemplate="%{text}<br>Mean expression: %{x:.1f}<br>Gene-wise dispersion: %{y:.4f}<extra></extra>",
    ))

    fitted = df.dropna(subset=["dispFit"]).sort_values("baseMean")
    fitted = fitted[fitted["dispFit"] > 0]
    fig.add_trace(go.Scatter(
        x=fitted["baseMean"], y=fitted["dispFit"],
        mode="lines", name="Fitted trend",
        line=dict(color="#d62728", width=2),
    ))

    final = df.dropna(subset=["dispersion"])
    final = final[final["dispersion"] > 0]
    fig.add_trace(go.Scatter(
        x=final["baseMean"], y=final["dispersion"],
        mode="markers", name="Final (shrunk) estimate",
        marker=dict(size=4, color="#1f77b4", opacity=0.6),
        text=final["gene_id"],
        hovertemplate="%{text}<br>Mean expression: %{x:.1f}<br>Final dispersion: %{y:.4f}<extra></extra>",
    ))

    fig.update_xaxes(title_text="Mean of normalized counts", type="log")
    fig.update_yaxes(title_text="Dispersion", type="log")
    fig.update_layout(height=500)
    return fig


def _plot_pvalue_histogram(pvalue_shape_result, bar_color="#636EFA"):
    """
    Build a bar-chart histogram of a contrast's p-value distribution,
    using the same bin counts/edges already computed by
    deseq2_manager.classify_pvalue_histogram_shape() -- so the plotted
    chart and its accompanying plain-language classification/flag are
    guaranteed to describe the exact same data.

    bar_color: the bars' fill color. The trace is given the name
        "P-values" specifically so _apply_plot_style's color_map
        (matched by trace name -- see _render_plot_style_controls'
        per-group color pickers) can recolor it like any other plot in
        this workspace, even though a histogram only has a single
        "group" rather than multiple colored series.
    """
    counts = pvalue_shape_result["counts"]
    bin_edges = pvalue_shape_result["bin_edges"]
    bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(counts))]
    bin_width = bin_edges[1] - bin_edges[0] if len(bin_edges) > 1 else 0.05

    fig = go.Figure(data=go.Bar(
        x=bin_centers, y=counts, width=bin_width * 0.95,
        name="P-values",
        marker=dict(color=bar_color),
        hovertemplate="p-value: %{x:.2f}<br>Gene count: %{y}<extra></extra>",
    ))
    fig.update_xaxes(title_text="p-value", range=[0, 1])
    fig.update_yaxes(title_text="Number of genes")
    fig.update_layout(height=400, showlegend=False)
    return fig


def _plot_ma(results_df, padj_threshold=0.05, lfc_threshold=1.0, gene_name_map=None):
    """
    Build an MA plot (mean expression vs. log2 fold change), using the
    exact same up/down/not-significant classification and color scheme
    as the volcano plot (via dm.classify_regulation), so the two plots'
    legends/colors always agree with each other. Where the volcano plot
    answers "which genes changed and by how much", the MA plot is the
    correct diagnostic for whether normalization worked: a healthy MA
    plot is roughly symmetric around log2FoldChange = 0 across the full
    range of expression -- systematic drift specifically among
    low-expression genes usually indicates a normalization or
    compositional issue that the volcano plot alone can't reveal (see
    deseq2_manager.classify_ma_bias for the accompanying plain-language
    check).

    Mean expression (baseMean) is shown on a log-scaled x-axis, the
    standard convention for this plot.

    Returns the figure only (no separate "df" return value needed here,
    since this reuses results_df directly rather than building its own
    filtered/derived copy the way _plot_volcano does for gene labeling).
    """
    df = results_df.dropna(subset=["baseMean", "log2FoldChange"]).copy()
    df = df[df["baseMean"] > 0]
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
            x=subset["baseMean"], y=subset["log2FoldChange"],
            mode="markers",
            name=label,
            marker=dict(size=5, color=color, opacity=0.5),
            text=subset["display_name"],
            hovertemplate="%{text}<br>Mean expression: %{x:.1f}<br>log2FC: %{y:.2f}<extra></extra>",
        ))

    fig.add_hline(y=0, line_dash="dash", line_color="gray")

    fig.update_xaxes(title_text="Mean of normalized counts", type="log")
    fig.update_yaxes(title_text="log2 Fold Change")
    fig.update_layout(height=500)
    return fig


def _cluster_order_generic(matrix_df, axis=0):
    """
    Compute a hierarchical-clustering leaf order for either the rows
    (axis=0) or columns (axis=1) of an arbitrary numeric DataFrame,
    using average-linkage on Euclidean distance -- the same clustering
    approach as _cluster_order_from_distance_matrix, generalized to a
    non-square matrix (e.g. genes x samples) rather than requiring a
    pre-computed square distance matrix.

    Gracefully falls back to the existing (unclustered) order if scipy
    isn't installed or clustering otherwise fails, so this is always
    safe to call.

    Returns a list of labels (row or column names) in clustered order.
    """
    labels = list(matrix_df.index) if axis == 0 else list(matrix_df.columns)
    if len(labels) < 3:
        return labels
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        data = matrix_df.values if axis == 0 else matrix_df.values.T
        Z = linkage(data, method="average", metric="euclidean")
        order = leaves_list(Z)
        return [labels[i] for i in order]
    except Exception:
        return labels


def _plot_top_genes_heatmap(norm_counts_df, results_df, gene_name_map=None, n_top=25,
                             meta_df=None, group_column=None, group_samples=True):
    """
    Build a z-scored expression heatmap of the top N most significant
    genes (by adjusted p-value) from the current contrast, across every
    sample in the project -- the third of the standard "big three"
    RNA-seq QC/results plots (alongside the volcano and MA plots).
    Where the volcano/MA plots show WHICH genes are significant and by
    how much, this shows what those genes' ACTUAL expression pattern
    looks like across every sample -- a genuine sanity check that the
    statistical hits correspond to a believable biological pattern
    (e.g. do they separate samples the way you'd expect?), and a chance
    to spot subgroups or residual structure the model didn't fully
    absorb.

    Values are log2(x+1) transformed, then each gene (row) is z-scored
    -- mean 0, unit variance -- ACROSS SAMPLES, so genes with very
    different absolute expression levels sit on a comparable color
    scale (standard practice for this type of heatmap; otherwise a
    handful of very highly expressed genes would visually dominate the
    whole plot regardless of how interesting their pattern actually
    is). Rows (genes) are always reordered by hierarchical clustering,
    so genes with similar patterns end up adjacent to each other.

    norm_counts_df: the project's full normalized counts table (see
        deseq2_manager.read_normalized_counts) -- columns: gene_id,
        then one column per sample.
    results_df: the current contrast's results table (needs gene_id
        and padj columns) -- used only to rank/select which genes are
        "top".
    n_top: how many top genes (by smallest padj) to include.

    meta_df, group_column: if given (a metadata DataFrame with a
        "sample" column, and the name of the column to group by), each
        sample's group is appended to its column label (e.g. "sample1
        (control)"), and a color-coded annotation strip is drawn just
        above the sample columns -- the same convention already used
        by _plot_sample_distance_heatmap -- so it's visually obvious
        which samples belong to which group without needing to
        recognize every sample name.
    group_samples: if True (the default) AND a group_column is given,
        samples are physically ordered by group FIRST (grouped
        together, in the group's first-seen order in meta_df), with
        hierarchical clustering only used to order samples WITHIN each
        group -- directly answering "do my top genes actually separate
        my experimental groups?" by literally placing each group's
        samples together rather than letting pure expression
        similarity potentially interleave them. Set False to instead
        order every sample by pure expression-similarity clustering
        across the whole dataset, ignoring group membership entirely
        (useful if you specifically want to check whether clustering
        recovers your groups on its own, without them being forced
        together).

    Returns (fig, z_df) -- the figure AND the z-scored matrix itself
    (genes x samples, in the final displayed order, gene_id as the
    index) so the caller can offer it as a downloadable CSV.
    """
    import numpy as np

    ranked = results_df.dropna(subset=["padj"]).sort_values("padj")
    top_gene_ids = ranked["gene_id"].astype(str).head(n_top).tolist()

    norm_counts_indexed = norm_counts_df.set_index(norm_counts_df["gene_id"].astype(str))
    sample_cols = [c for c in norm_counts_df.columns if c != "gene_id"]

    available_ids = [g for g in top_gene_ids if g in norm_counts_indexed.index]
    subset = norm_counts_indexed.loc[available_ids, sample_cols].astype(float)

    log_subset = np.log2(subset + 1)
    row_mean = log_subset.mean(axis=1)
    row_std = log_subset.std(axis=1).replace(0, 1)  # avoid divide-by-zero for a rare constant-expression gene
    z_df = log_subset.sub(row_mean, axis=0).div(row_std, axis=0)

    gene_order = _cluster_order_generic(z_df, axis=0)

    # sample_to_group is only built when group_samples is True -- when the
    # user checks "Ignore grouping" (group_samples=False), this stays empty
    # on purpose, so every "if sample_to_group:" block below (display
    # labels, the color-coded annotation strip, and its legend) correctly
    # skips rendering entirely. Showing a group-colored strip/legend next
    # to samples that were explicitly ordered WITHOUT regard to group
    # membership would be misleading (the colors would appear scattered
    # rather than blocked, since that's the whole point of this mode).
    sample_to_group = {}
    if group_samples and meta_df is not None and group_column and group_column in meta_df.columns:
        sample_to_group = dict(
            zip(meta_df["sample"].astype(str), meta_df[group_column].astype(str))
        )

    if sample_to_group and group_samples:
        # --- Group samples together, cluster only WITHIN each group ---
        # Groups are ordered by first appearance in meta_df (usually
        # matches the order the user entered/uploaded their metadata,
        # e.g. control before treated) rather than alphabetically,
        # which can otherwise produce a counter-intuitive order (e.g.
        # "treated" before "control" alphabetically).
        seen_groups = []
        for s in z_df.columns:
            g = sample_to_group.get(s)
            if g is not None and g not in seen_groups:
                seen_groups.append(g)

        sample_order = []
        for g in seen_groups:
            members = [s for s in z_df.columns if sample_to_group.get(s) == g]
            if len(members) >= 3:
                sub_order = _cluster_order_generic(z_df[members], axis=1)
            else:
                sub_order = members
            sample_order.extend(sub_order)
        # Any sample with no group mapping (shouldn't normally happen)
        # is appended at the end rather than silently dropped.
        ungrouped = [s for s in z_df.columns if s not in sample_order]
        sample_order.extend(ungrouped)
    else:
        # No group column available, or the user opted for pure
        # expression-based clustering across the whole dataset.
        sample_order = _cluster_order_generic(z_df, axis=1)

    z_df = z_df.loc[gene_order, sample_order]

    gene_name_map = gene_name_map or {}
    display_gene_labels = [gene_name_map.get(g, g) for g in z_df.index]

    if sample_to_group:
        display_sample_labels = [
            f"{s} ({sample_to_group.get(s, '?')})" for s in z_df.columns
        ]
    else:
        display_sample_labels = list(z_df.columns)

    fig = go.Figure(data=go.Heatmap(
        z=z_df.values,
        x=display_sample_labels,
        y=display_gene_labels,
        colorscale="RdBu_r",
        zmid=0,
        colorbar=dict(title="Z-score", x=1.08),
        hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>",
    ))

    n_genes = len(z_df.index)
    fig_height = max(400, 22 * n_genes)

    if sample_to_group:
        # --- Color-coded group annotation strip, above the columns ---
        # A slim second heatmap trace spanning one row, positioned via
        # its own y-axis (yaxis2) directly above the main heatmap --
        # one colored cell per sample (in the same group-then-cluster
        # order as the main heatmap), so which samples share a group is
        # visible purely from color, at a glance. Mirrors the same
        # "annotation sidebar" convention already used by
        # _plot_sample_distance_heatmap, just rotated 90 degrees since
        # here the samples are COLUMNS rather than rows.
        groups_in_order = [sample_to_group.get(s, "?") for s in z_df.columns]
        unique_groups = seen_groups if (sample_to_group and group_samples) else sorted(set(groups_in_order))
        group_to_code = {g: i for i, g in enumerate(unique_groups)}
        z_strip = [[group_to_code.get(g, 0) for g in groups_in_order]]

        n_groups = max(len(unique_groups), 1)
        strip_colors = [
            _DEFAULT_COLOR_PALETTE[i % len(_DEFAULT_COLOR_PALETTE)] for i in range(n_groups)
        ]
        if n_groups == 1:
            colorscale = [[0, strip_colors[0]], [1, strip_colors[0]]]
        else:
            colorscale = []
            for i, color in enumerate(strip_colors):
                colorscale.append([i / n_groups, color])
                colorscale.append([(i + 1) / n_groups, color])

        fig.add_trace(go.Heatmap(
            z=z_strip,
            x=display_sample_labels,
            y=["Group"],
            yaxis="y2",
            colorscale=colorscale,
            zmin=0, zmax=n_groups,
            showscale=False,
            text=[groups_in_order],
            hovertemplate="%{x}<br>Group: %{text}<extra></extra>",
        ))

        # Reserve a slim band at the top of the figure for the strip
        # (yaxis2, above) and give the main heatmap the rest (yaxis,
        # below) -- domains are fractions of the TOTAL figure height,
        # so a fixed ~30px-equivalent band works reasonably regardless
        # of how many genes are shown.
        strip_fraction = min(0.12, 30 / fig_height)
        fig.update_layout(
            yaxis=dict(domain=[0, 1 - strip_fraction - 0.02]),
            yaxis2=dict(domain=[1 - strip_fraction, 1], showticklabels=True),
            xaxis=dict(tickangle=-45),
        )

        # Manual legend (discrete group labels read more clearly here
        # than the strip's own hidden color axis) -- one invisible-
        # marker scatter trace per group, purely so its name/color
        # shows up in the figure's standard legend. Positioned INSIDE
        # the main heatmap's bottom-right corner by default (rather
        # than Plotly's default right-of-figure placement, which sits
        # directly under/beside the colorbar occupying that space
        # already) -- the user's own style panel choice still overrides
        # this if they pick something else (see _apply_plot_style).
        fig.update_layout(legend=dict(
            yanchor="bottom", y=0.01, xanchor="right", x=0.99,
            bgcolor="rgba(255,255,255,0.7)", bordercolor="rgba(0,0,0,0.2)", borderwidth=1,
        ))
        for group_name, color in zip(unique_groups, strip_colors):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=10, color=color),
                name=str(group_name),
                showlegend=True,
            ))
    else:
        fig.update_layout(xaxis=dict(tickangle=-45))

    fig.update_layout(height=fig_height)
    return fig, z_df


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

    counts_df = _load_counts_df(counts_matrix_path)
    meta_df = _load_meta_df(metadata_path)

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

    pca_df, pct_var = _load_pca_coordinates(deseq2_out_dir)
    if pca_df is not None:
        color_options = [c for c in pca_df.columns if c not in ("PC1", "PC2", "PC3", "PC4", "sample")]
        color_by = st.selectbox("Color points by:", options=color_options, key="pca_color_select")
        group_values = sorted(pca_df[color_by].astype(str).unique())

        pca_adj_df, pct_var_adj = (None, None)
        if batch_column:
            pca_adj_df, pct_var_adj = _load_pca_coordinates(deseq2_out_dir, batch_adjusted=True)

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
            # Stacked vertically (rather than side-by-side columns) so
            # each PCA plot has full width to breathe -- two narrow
            # half-width plots side by side made both harder to read,
            # especially once group labels/legends/confidence ellipses
            # are added.
            st.markdown("#### 🔹 Before Batch Correction")
            show_ellipse_before, confidence_before = _render_confidence_ellipse_controls("pca_before")
            fig_before, ellipse_stats_before = _plot_pca(
                pca_df, pct_var, color_by,
                show_confidence_ellipse=show_ellipse_before, confidence_level=confidence_before,
            )
            pca_before_style = _render_plot_style_controls("pca_before", group_values=group_values)
            _apply_plot_style(
                fig_before, pca_before_style,
                default_title="PCA: Sample Similarity — Before Batch Correction",
            )
            _render_plotly_chart(fig_before)
            _render_ellipse_stats(ellipse_stats_before, "pca_before")
            _render_pdf_export(fig_before, "pca_before", "pca_before_batch_correction")
            _render_csv_download(
                pca_df, "pca_coordinates_before_batch_correction", "pca_before",
                expander_label="⬇️ Download PCA Coordinates (Before).csv",
                help_text="Every sample's PC1-PC4 coordinates plus every metadata column available for coloring. Percent variance explained and batch-effect size are in the QC summary download below.",
            )

            st.markdown("---")

            st.markdown("#### 🔸 After Batch Correction")
            show_ellipse_after, confidence_after = _render_confidence_ellipse_controls("pca_after")
            fig_after, ellipse_stats_after = _plot_pca(
                pca_adj_df, pct_var_adj, color_by,
                show_confidence_ellipse=show_ellipse_after, confidence_level=confidence_after,
            )
            pca_after_style = _render_plot_style_controls("pca_after", group_values=group_values)
            _apply_plot_style(
                fig_after, pca_after_style,
                default_title="PCA: Sample Similarity — After Batch Correction (Visualization Only)",
            )
            _render_plotly_chart(fig_after)
            _render_ellipse_stats(ellipse_stats_after, "pca_after")
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
            show_ellipse_single, confidence_single = _render_confidence_ellipse_controls("pca_single")
            fig_single, ellipse_stats_single = _plot_pca(
                pca_df, pct_var, color_by,
                show_confidence_ellipse=show_ellipse_single, confidence_level=confidence_single,
            )
            pca_single_style = _render_plot_style_controls("pca_single", group_values=group_values)
            _apply_plot_style(fig_single, pca_single_style, default_title="PCA: Sample Similarity")
            _render_plotly_chart(fig_single)
            _render_ellipse_stats(ellipse_stats_single, "pca_single")
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

    # --- Sample-to-sample distance heatmap (companion QC view to PCA) ---
    st.subheader("🌡️ Sample Distance Heatmap")
    with st.expander("ℹ️ How to read this plot"):
        st.markdown(
            "This is the standard companion view to the PCA plot above: "
            "every pair of samples is compared directly (using the same "
            "variance-stabilized data as PCA), rather than only looking "
            "at the first few principal components. Samples are grouped "
            "by similarity (hierarchical clustering) so any expected -- "
            "or unexpected -- grouping is easy to spot at a glance. "
            "Darker cells mean more similar samples. Each sample's group "
            "(picked below) is also shown as a colored strip alongside "
            "the heatmap and appended to its label, so groupings are "
            "readable even without recognizing every sample name."
        )

    distance_df = _load_sample_distance_matrix(deseq2_out_dir)
    if distance_df is not None:
        # --- Pick which metadata column to color/label samples by ---
        # Shown BEFORE building the heatmap figure (rather than after,
        # like the clustering-mismatch check below) since this same
        # selection now drives the heatmap's own color strip and axis
        # labels, not just the mismatch check.
        heatmap_group_columns = [c for c in meta_df.columns if c != "sample"]
        default_group_index = 0
        if saved_config.get("design_columns") and saved_config["design_columns"][0] in heatmap_group_columns:
            default_group_index = heatmap_group_columns.index(saved_config["design_columns"][0])
        heatmap_group_column = None
        heatmap_show_group_annotation = True
        if heatmap_group_columns:
            col_group_select, col_hide_annotation = st.columns([2, 1.2])
            with col_group_select:
                heatmap_group_column = st.selectbox(
                    "Color/label samples by which metadata column?",
                    options=heatmap_group_columns,
                    index=default_group_index,
                    key="sample_dist_check_column",
                    help="Colors a strip alongside the heatmap and appends this column's value to each sample's label, so groupings are visible at a glance -- also used below to flag any sample whose closest match is in a different group.",
                )
            with col_hide_annotation:
                heatmap_show_group_annotation = not st.checkbox(
                    "Hide grouping bar and legend",
                    value=False,
                    key="sample_dist_hide_group_annotation",
                    help="Removes the color-coded strip and legend from the plot below, showing plain sample names instead -- the metadata column selected above is still used for the clustering-mismatch check further down.",
                )

        dist_fig = _plot_sample_distance_heatmap(
            distance_df, meta_df=meta_df, group_column=heatmap_group_column,
            show_group_annotation=heatmap_show_group_annotation,
        )
        dist_style = _render_plot_style_controls("sample_dist", group_values=None, show_legend_controls=False)
        _apply_plot_style(dist_fig, dist_style, default_title="Sample-to-Sample Distance")
        _render_plotly_chart(dist_fig)
        _render_pdf_export(dist_fig, "sample_dist", "sample_distance_heatmap")

        # --- Smart flag: does any sample's nearest neighbor belong to a different group? ---
        if heatmap_group_column:
            mismatches = dm.detect_sample_clustering_mismatch(distance_df, meta_df, heatmap_group_column)
            if mismatches:
                mismatch_lines = "\n".join(
                    f"- **{m['sample']}** (`{heatmap_group_column}` = {m['own_group']}) is most similar to "
                    f"**{m['nearest_neighbor']}** (`{heatmap_group_column}` = {m['neighbor_group']})"
                    for m in mismatches
                )
                st.warning(
                    f"⚠️ **{len(mismatches)} sample(s)** cluster more closely with a "
                    f"*different* `{heatmap_group_column}` group than their own:\n\n{mismatch_lines}\n\n"
                    "This doesn't necessarily mean something is wrong -- but it's "
                    "worth double-checking these samples' labeling, or considering "
                    "whether they should be excluded, especially if this doesn't "
                    "match your experimental expectations."
                )
            else:
                st.success(f"✅ Every sample's closest match shares its own `{heatmap_group_column}` group -- no clustering surprises.")

        _render_csv_download(
            distance_df.reset_index(), "sample_distance_matrix", "sample_dist",
            expander_label="⬇️ Download Sample Distance Matrix (.csv)",
            help_text="The full sample x sample Euclidean distance matrix (computed on variance-stabilized counts).",
        )
    else:
        st.info("Sample distance data isn't available for this project yet -- re-run DESeq2 above to generate it.")

    st.markdown("---")

    # --- Model-fit diagnostics: dispersion plot + size factor QC ---
    st.subheader("🔧 Model Fit Diagnostics")

    dispersion_df = _load_dispersion_estimates(deseq2_out_dir)
    if dispersion_df is not None:
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "DESeq2 estimates how much each gene's counts vary "
                "across your samples beyond simple random (Poisson) "
                "noise -- this is its **dispersion**. Because per-gene "
                "estimates (gray points) are unreliable with typical "
                "sample sizes, DESeq2 shrinks them toward a fitted trend "
                "(red line) shared across all genes, producing the final "
                "values (blue points) actually used in the statistical "
                "test. Most genes scattering loosely around the trend is "
                "completely normal -- this plot (and the flag below) is "
                "only meant to catch an unusually poor overall fit."
            )
        disp_fig = _plot_dispersion(dispersion_df)
        disp_style = _render_plot_style_controls("dispersion", group_values=None, show_legend_controls=True)
        _apply_plot_style(disp_fig, disp_style, default_title="Dispersion Estimates")
        _render_plotly_chart(disp_fig)
        _render_pdf_export(disp_fig, "dispersion", "dispersion_estimates")

        fit_assessment = dm.assess_dispersion_fit(dispersion_df)
        if fit_assessment["flagged"]:
            st.warning(fit_assessment["message"])
        else:
            st.success(fit_assessment["message"])

        _render_csv_download(
            dispersion_df, "dispersion_estimates", "dispersion",
            expander_label="⬇️ Download Dispersion Estimates (.csv)",
            help_text="Per-gene mean expression, gene-wise dispersion estimate, fitted trend value, and final (shrunk) dispersion used in the statistical test.",
        )
    else:
        st.info("Dispersion data isn't available for this project yet -- re-run DESeq2 above to generate it.")

    size_factor_df = _load_size_factors(deseq2_out_dir)
    if size_factor_df is not None:
        st.markdown("##### ⚖️ Size Factors (Normalization QC)")
        st.caption(
            "DESeq2's per-sample normalization factors -- these correct "
            "for differences in sequencing depth/composition between "
            "samples before any comparison is made. A size factor far "
            "from the others can indicate a failed, degraded, or "
            "unusually composed library."
        )
        st.dataframe(size_factor_df, use_container_width=True, hide_index=True)

        sf_outliers = dm.flag_size_factor_outliers(size_factor_df)
        if sf_outliers["flagged_samples"]:
            outlier_lines = "\n".join(
                f"- **{o['sample']}**: size factor {o['size_factor']} "
                f"({o['ratio_to_median']}x the median of {sf_outliers['median_size_factor']})"
                for o in sf_outliers["flagged_samples"]
            )
            st.warning(
                f"⚠️ **{len(sf_outliers['flagged_samples'])} sample(s)** have "
                f"a size factor notably different from the rest:\n\n{outlier_lines}\n\n"
                "Consider double-checking these samples' sequencing depth "
                "and library quality."
            )
        else:
            st.success("✅ All samples' size factors are within a reasonable range of each other.")

        _render_csv_download(
            size_factor_df, "size_factors", "size_factors",
            expander_label="⬇️ Download Size Factors (.csv)",
        )
    else:
        st.info("Size factor data isn't available for this project yet -- re-run DESeq2 above to generate it.")

    st.markdown("---")

    # --- Per-contrast results + volcano ---
    st.subheader("📊 Results by Contrast")

    selected_contrast = st.selectbox("Select a contrast to view:", options=available_contrasts, key="view_contrast_select")

    results_df = _load_contrast_results(deseq2_out_dir, selected_contrast)
    if results_df is not None:
        # --- Gene ID -> Gene Name Mapping + significance/fold-change
        # thresholds, placed directly together right above the volcano
        # plot -- both of these most directly shape what the volcano plot
        # actually shows (readable gene names in hover/labels, and the
        # dashed threshold lines + up/down/not-significant coloring), so
        # they live right next to it rather than further up the page.
        gene_name_map = _render_gene_id_mapping_panel(project, counts_df)

        st.caption(
            "The two thresholds below set the volcano plot's dashed "
            "threshold lines and coloring (and are also reused further "
            "down for the MA plot's coloring and the results table's "
            "Regulation column)."
        )
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

        # --- Gene labeling: "Label All" + search + per-gene style controls ---
        # Built (and its resulting annotations attached to volcano_fig)
        # BEFORE the chart is rendered below, so labels added on a
        # previous run/rerun (including via the search box inside
        # _render_gene_label_controls) show up immediately. Newly
        # CLICKED genes are handled just after rendering (see the
        # on_select handling below it), since a click's result is only
        # available as st.plotly_chart's return value -- that discovery
        # triggers one extra st.rerun() so the clicked gene's label
        # appears using the same code path as everything else, keeping
        # this logic in one place rather than duplicated.
        labeled_genes, sig_gene_ids = _render_gene_label_controls(
            selected_contrast, volcano_df, gene_name_map=gene_name_map,
        )
        gene_annotations = _build_gene_label_annotations(volcano_df, labeled_genes, gene_name_map=gene_name_map)
        if gene_annotations:
            volcano_fig.update_layout(annotations=list(volcano_fig.layout.annotations) + gene_annotations)

        st.caption(
            "💡 Tip: click any point on the chart below to add a label "
            "for that gene, or use the search box above to find and "
            "label a gene by ID/symbol regardless of where it falls on "
            "the plot. Drag any label to reposition it -- its arrow "
            "will keep pointing at the correct dot."
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

        st.markdown("---")

        # --- MA plot: the correct diagnostic for normalization issues ---
        # Reuses the exact same padj_cutoff/lfc_cutoff (and gene_name_map)
        # already set above via the sliders next to the volcano plot, so
        # this plot's coloring/hover names always match your current
        # selection with no separate controls needed here.
        st.markdown("##### 📉 MA Plot")
        with st.expander("ℹ️ How to read this plot, and why it's different from the volcano plot"):
            st.markdown(
                "The volcano plot above answers "
                "**\"which genes changed, and by how much?\"**. This plot "
                "answers a different question: **\"did normalization work "
                "correctly across the full range of expression levels?\"** "
                "A healthy MA plot is roughly symmetric around "
                "log2FoldChange = 0 at every expression level (x-axis). "
                "Genes systematically drifting away from zero "
                "specifically at LOW expression levels usually signals a "
                "normalization or compositional issue that the volcano "
                "plot alone can't reveal."
            )
        ma_fig = _plot_ma(results_df, padj_cutoff, lfc_cutoff, gene_name_map=gene_name_map)
        ma_style = _render_plot_style_controls(
            f"ma_{selected_contrast}",
            group_values=["Up-regulated", "Down-regulated", "Not significant"],
            default_colors={
                "Up-regulated": "#d62728",
                "Down-regulated": "#1f77b4",
                "Not significant": "#7f7f7f",
            },
        )
        _apply_plot_style(ma_fig, ma_style, default_title=f"MA Plot: {selected_contrast}")
        _render_plotly_chart(ma_fig)
        _render_pdf_export(ma_fig, f"ma_{selected_contrast}", f"ma_plot_{selected_contrast}")

        ma_bias = dm.classify_ma_bias(results_df)
        if ma_bias["flagged"]:
            st.warning(ma_bias["message"])
        else:
            st.success(ma_bias["message"])

        st.markdown("---")

        # --- P-value histogram: a cheap, high-value sanity check ---
        # Independent of any significance/fold-change threshold -- this
        # looks at the raw p-value distribution itself, so it doesn't
        # depend on the sliders above at all.
        st.markdown("##### 📊 P-value Distribution")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "One of the simplest, most informative sanity checks "
                "after any differential expression test. A **healthy** "
                "result usually shows a spike of small p-values near 0 "
                "(genuine differential signal) sitting on top of a "
                "roughly flat background elsewhere (non-differential "
                "genes, which are uniformly distributed by definition). "
                "Other shapes -- a flat distribution with no spike, an "
                "elevated hump in the middle, or a spike near 1 -- each "
                "point to a specific, well-known issue, explained below."
            )
        pvalue_shape = dm.classify_pvalue_histogram_shape(results_df)
        if pvalue_shape["counts"]:
            pval_fig = _plot_pvalue_histogram(pvalue_shape)
            pval_style = _render_plot_style_controls(
                f"pval_{selected_contrast}",
                group_values=["P-values"],
                default_colors={"P-values": "#636EFA"},
                show_legend_controls=False,
            )
            _apply_plot_style(pval_fig, pval_style, default_title=f"P-value Distribution: {selected_contrast}")
            _render_plotly_chart(pval_fig)
            _render_pdf_export(pval_fig, f"pval_{selected_contrast}", f"pvalue_histogram_{selected_contrast}")

            if pvalue_shape["shape"] in ("healthy",):
                st.success(pvalue_shape["message"])
            elif pvalue_shape["shape"] == "flat_uniform":
                st.info(pvalue_shape["message"])
            else:
                st.warning(pvalue_shape["message"])
        else:
            st.info(pvalue_shape["message"])

        st.markdown("---")

        # --- Top Genes heatmap: what do the actual hits look like? ---
        # The third of the standard "big three" RNA-seq plots (alongside
        # the volcano and MA plots above) -- see _plot_top_genes_heatmap's
        # docstring for the full rationale. Uses the project's full
        # normalized counts table (already exported by the DESeq2 R
        # script for other purposes -- see deseq2_manager.read_normalized_counts),
        # so no extra computation is needed beyond what DESeq2 already
        # produces during a normal run.
        st.markdown("##### 🧬 Top Genes")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "The volcano and MA plots above tell you WHICH genes are "
                "significant and by how much -- this plot shows what "
                "those top genes' ACTUAL expression pattern looks like "
                "across every sample. Each row is one gene, each column "
                "is one sample, and color shows that gene's expression "
                "relative to its own average (a **z-score**: red = "
                "higher than usual for that gene, blue = lower), so "
                "genes with very different absolute expression levels "
                "are still shown on a comparable scale. Genes are "
                "reordered by similarity (hierarchical clustering); by "
                "default, samples are grouped together by the metadata "
                "column you pick below (with a colored strip and label "
                "showing each sample's group), so you can directly see "
                "whether these top genes actually separate your "
                "experimental groups the way you'd expect -- not just "
                "cluster by generic expression similarity.\n\n"
                "This is a genuine sanity check, not just a "
                "restatement of the volcano plot: if a sample's "
                "expression pattern doesn't match its group, or you "
                "see an unexpected subgroup, that's worth investigating "
                "further."
            )

        norm_counts_df = _load_normalized_counts(deseq2_out_dir)
        if norm_counts_df is not None:
            top_genes_group_columns = [c for c in meta_df.columns if c != "sample"]
            top_genes_group_column = None
            top_genes_group_samples = True
            col_n, col_group, col_order = st.columns([1, 1.4, 1.4])
            with col_n:
                n_top_genes = st.number_input(
                    "Number of top genes:", min_value=5, max_value=100, value=25, step=5,
                    key=f"top_genes_n_{selected_contrast}",
                    help="Genes are ranked by smallest adjusted p-value (padj) in this contrast.",
                )
            if top_genes_group_columns:
                with col_group:
                    default_tg_group_index = 0
                    if saved_config.get("design_columns") and saved_config["design_columns"][0] in top_genes_group_columns:
                        default_tg_group_index = top_genes_group_columns.index(saved_config["design_columns"][0])
                    top_genes_group_column = st.selectbox(
                        "Group/label samples by:",
                        options=top_genes_group_columns,
                        index=default_tg_group_index,
                        key=f"top_genes_group_column_{selected_contrast}",
                        help="Colors a strip above the sample columns and appends this column's value to each sample's label.",
                    )
                with col_order:
                    top_genes_group_samples = not st.checkbox(
                        "Ignore grouping (cluster samples purely by expression)",
                        value=False,
                        key=f"top_genes_pure_cluster_{selected_contrast}",
                        help="By default, samples are physically grouped together by the column above, clustering only within each group. Check this to instead order every sample purely by expression similarity across the whole dataset, ignoring group membership.",
                    )

            top_genes_fig, top_genes_z_df = _plot_top_genes_heatmap(
                norm_counts_df, results_df, gene_name_map=gene_name_map, n_top=n_top_genes,
                meta_df=meta_df, group_column=top_genes_group_column, group_samples=top_genes_group_samples,
            )
            top_genes_style = _render_plot_style_controls(
                f"top_genes_{selected_contrast}", group_values=None, show_legend_controls=False,
            )
            _apply_plot_style(top_genes_fig, top_genes_style, default_title=f"Top {n_top_genes} Genes: {selected_contrast}")
            _render_plotly_chart(top_genes_fig)
            _render_pdf_export(top_genes_fig, f"top_genes_{selected_contrast}", f"top_genes_heatmap_{selected_contrast}")

            _render_csv_download(
                top_genes_z_df.reset_index().rename(columns={"index": "gene_id"}),
                f"top_{n_top_genes}_genes_zscores_{selected_contrast}", f"top_genes_{selected_contrast}",
                expander_label=f"⬇️ Download Top {n_top_genes} Genes Z-Scores (.csv)",
                help_text="Log2-transformed, per-gene z-scored expression values for the top genes shown above, across every sample.",
            )
        else:
            st.info("Normalized counts aren't available for this project yet -- re-run DESeq2 above to generate them.")

        st.markdown("---")

        # --- Explain NA genes: Cook's-distance outliers + independent filtering ---
        # These are two distinct, legitimate DESeq2 behaviors that silently
        # produce NA values in the results table right below -- surfaced
        # here explicitly (see deseq2_manager.summarize_na_genes) so the
        # user understands WHY a gene is missing/blank rather than just
        # seeing an unexplained NA. Placed here (using gene_name_map,
        # already computed above for the volcano plot) rather than earlier
        # on the page, since it explains the table that follows directly
        # below.
        na_summary = dm.summarize_na_genes(
            results_df, filter_threshold=dm.read_filter_threshold(deseq2_out_dir, selected_contrast),
        )
        if na_summary["n_cooks_outliers"] or na_summary["n_low_count_filtered"]:
            with st.expander(
                f"ℹ️ Why do some genes show NA? ({na_summary['n_cooks_outliers'] + na_summary['n_low_count_filtered']} gene(s) affected)"
            ):
                st.markdown(na_summary["message"])
                if na_summary["cooks_outlier_genes"]:
                    display_genes = [gene_name_map.get(g, g) for g in na_summary["cooks_outlier_genes"]]
                    st.caption(
                        "Genes flagged as Cook's-distance outliers "
                        f"(showing up to 20): {', '.join(display_genes)}"
                    )

        st.markdown("---")

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
