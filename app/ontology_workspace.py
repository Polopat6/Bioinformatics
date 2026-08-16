"""
ontology_workspace.py

Ontology Analysis workspace: runs GO / KEGG / Reactome enrichment via
clusterProfiler (+ ReactomePA) on a project's DESeq2 results, in three
flavors: ORA, GSEA, and compareCluster (multi-contrast comparison).

This picks up where differential_expression_workspace.py leaves off,
reusing the same active project and the clusterProfiler-ready gene list
CSVs already exported there.

This module is fully self-contained beyond the DE-workspace styling
import -- all Ontology Analysis development should happen here.
"""

import os

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import project_manager as pm
import deseq2_manager as dm
import reference_manager as rm
import gene_id_mapper as gim
import ontology_manager as om
import differential_expression_workspace as dew

WORKSPACE_KEY = "bulk_rnaseq"

SEQUENTIAL_COLORSCALE_OPTIONS = dew.SEQUENTIAL_COLORSCALE_OPTIONS
DIVERGING_COLORSCALE_OPTIONS = dew.DIVERGING_COLORSCALE_OPTIONS

_CATEGORY_BAND_PALETTE = dew._DEFAULT_COLOR_PALETTE


def _hex_to_rgba(hex_color, alpha=0.15):
    "Convert a '#RRGGBB' hex color string to an 'rgba(r,g,b,alpha)' string, for use as a low-opacity background fill."
    hex_color = hex_color.lstrip("#")
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# ---------------------------------------------------------------------------
# Network layout helper (for the enrichment map + gene-concept network plots)
# ---------------------------------------------------------------------------

def _compute_network_layout(node_ids, edges_df, seed=42):
    node_ids = list(node_ids)
    try:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(node_ids)
        if edges_df is not None and not edges_df.empty:
            for _, row in edges_df.iterrows():
                weight = row["weight"] if "weight" in edges_df.columns else 1.0
                G.add_edge(row["source"], row["target"], weight=weight)
        pos = nx.spring_layout(G, seed=seed, weight="weight" if "weight" in (edges_df.columns if edges_df is not None else []) else None)
        return {node: (float(coord[0]), float(coord[1])) for node, coord in pos.items()}
    except ImportError:
        import math
        import random

        n = len(node_ids)
        if n == 0:
            return {}
        rng = random.Random(seed)
        start_angle = rng.uniform(0, 2 * math.pi)
        return {
            node: (
                math.cos(start_angle + 2 * math.pi * i / n),
                math.sin(start_angle + 2 * math.pi * i / n),
            )
            for i, node in enumerate(node_ids)
        }


def _render_layout_seed_control(key_prefix):
    seed_key = f"{key_prefix}_layout_seed"
    if seed_key not in st.session_state:
        st.session_state[seed_key] = 42
    if st.button("🔀 Try a different layout", key=f"{key_prefix}_regenerate_layout_btn",
                 help="Re-arranges the network using a different starting layout -- useful if the current arrangement looks too cluttered or overlapping."):
        st.session_state[seed_key] += 1
    return st.session_state[seed_key]


def _build_node_label_annotations(node_ids, positions, display_labels, colors=None, font_size=11,
                                   ax=25, ay=-20):
    colors = colors or ["black"] * len(node_ids)
    annotations = []
    for i, nid in enumerate(node_ids):
        if nid not in positions:
            continue
        x, y = positions[nid]
        annotations.append(dict(
            x=x, y=y, xref="x", yref="y",
            text=str(display_labels[i]),
            showarrow=True, arrowhead=1, arrowsize=1, arrowwidth=1,
            arrowcolor="rgba(120,120,120,0.6)",
            ax=ax, ay=ay,
            font=dict(size=font_size, color=colors[i]),
            name=f"node_label_{nid}",
        ))
    return annotations


# ---------------------------------------------------------------------------
# Term label editor (manual shortening/renaming of GO/KEGG/Reactome term names)
# ---------------------------------------------------------------------------

def _render_term_label_editor(display_df, key_prefix, id_col="ID", desc_col="Description",
                               label="✏️ Shorten/rename term labels shown on this plot (optional)"):
    if display_df is None or display_df.empty:
        return {}

    with st.expander(label):
        st.caption(
            "Edit the \"Custom Label\" column below to shorten or "
            "rename any term as it appears on the plot -- this only "
            "changes what's DISPLAYED, your actual results are "
            "unaffected."
        )
        editor_df = pd.DataFrame({
            "ID": display_df[id_col].astype(str),
            "Original Description": display_df[desc_col].astype(str) if desc_col in display_df.columns else display_df[id_col].astype(str),
        })
        editor_df["Custom Label"] = editor_df["Original Description"]

        edited = st.data_editor(
            editor_df,
            use_container_width=True,
            hide_index=True,
            disabled=["ID", "Original Description"],
            key=f"{key_prefix}_term_label_editor",
        )

        label_map = {}
        for _, row in edited.iterrows():
            custom = str(row["Custom Label"]).strip()
            original = str(row["Original Description"]).strip()
            if custom and custom != original:
                label_map[row["ID"]] = custom
        return label_map


def _render_gene_label_editor(gene_ids, key_prefix):
    if not gene_ids:
        return {}

    with st.expander("✏️ Shorten/rename gene labels shown on this plot (optional)"):
        st.caption(
            "Edit the \"Custom Label\" column below to rename any "
            "gene's displayed label (e.g. if you'd prefer a different "
            "symbol/alias) -- this only changes what's DISPLAYED."
        )
        editor_df = pd.DataFrame({"Gene ID": sorted(gene_ids)})
        editor_df["Custom Label"] = editor_df["Gene ID"]

        edited = st.data_editor(
            editor_df, use_container_width=True, hide_index=True,
            disabled=["Gene ID"], key=f"{key_prefix}_gene_label_editor",
        )

        label_map = {}
        for _, row in edited.iterrows():
            custom = str(row["Custom Label"]).strip()
            original = str(row["Gene ID"]).strip()
            if custom and custom != original:
                label_map[row["Gene ID"]] = custom
        return label_map


def _apply_term_labels(df, label_map, id_col="ID", desc_col="Description"):
    "Return a copy of df with desc_col replaced by any custom label(s) from label_map. Non-destructive."
    if not label_map:
        return df
    result = df.copy()
    result[desc_col] = [
        label_map.get(row_id, original_desc)
        for row_id, original_desc in zip(result[id_col], result[desc_col])
    ]
    return result


# ---------------------------------------------------------------------------
# Sort-by and x-axis-value controls (dot plot / bar plot)
# ---------------------------------------------------------------------------

def _get_sort_options(is_gsea, has_generatio):
    options = {"Significance (p.adjust, most significant first)": ("p.adjust", True)}
    if is_gsea:
        options["Gene set size"] = ("setSize", False)
        options["Enrichment Score (NES, most positive/up first)"] = ("NES", False)
        options["Enrichment Score (NES, most negative/down first)"] = ("NES", True)
    else:
        options["Gene count"] = ("Count", False)
        if has_generatio:
            options["Gene ratio"] = ("GeneRatio_numeric", False)
    return options


def _get_x_axis_options(is_gsea, has_generatio):
    if is_gsea:
        options = {"Normalized Enrichment Score (NES)": "NES"}
    else:
        options = {}
        if has_generatio:
            options["Gene Ratio"] = "GeneRatio_numeric"
        options["Gene Count"] = "Count"
    options["-log10(adjusted p-value)"] = "neg_log10_padj"
    return options


def _prepare_ranked_plot_df(result_df, sort_column, ascending, n_top):
    import math

    df = result_df.copy()
    if "neg_log10_padj" not in df.columns and "p.adjust" in df.columns:
        df["neg_log10_padj"] = -df["p.adjust"].clip(lower=1e-300).apply(math.log10)

    df = df.dropna(subset=[sort_column]).sort_values(sort_column, ascending=ascending).head(n_top)
    return df.iloc[::-1]


def _render_ranked_plot_controls(result_df, is_gsea, key_prefix):
    has_generatio = "GeneRatio_numeric" in result_df.columns

    n_top = st.slider(
        "Number of top terms to show:", min_value=5, max_value=50, value=20, step=5,
        key=f"{key_prefix}_n_top",
        help="How many terms are plotted, selected using the 'rank/select top terms by' choice below. Increasing this shows more (and less significant) terms; decreasing it focuses the plot on only the strongest hits.",
    )

    sort_options = _get_sort_options(is_gsea, has_generatio)
    sort_label = st.selectbox(
        "Rank/select top terms by:", options=list(sort_options.keys()),
        key=f"{key_prefix}_sort_by",
        help="Determines BOTH which terms count as the 'top N' above, AND their display order.",
    )
    sort_column, sort_ascending = sort_options[sort_label]
    if sort_column not in result_df.columns:
        st.info(f"ℹ️ '{sort_label}' isn't available for this result (missing column) -- falling back to significance.")
        sort_column, sort_ascending = "p.adjust", True

    prepared_df = _prepare_ranked_plot_df(result_df, sort_column, sort_ascending, n_top)

    x_axis_options = _get_x_axis_options(is_gsea, has_generatio)
    x_label = st.selectbox(
        "X-axis value:", options=list(x_axis_options.keys()),
        key=f"{key_prefix}_x_axis",
        help="What each term's horizontal position represents -- purely a display choice, does not change which terms are shown or their significance.",
    )
    x_column = x_axis_options[x_label]

    return prepared_df, x_column, x_label


# ---------------------------------------------------------------------------
# Plot builders
# ---------------------------------------------------------------------------

def _plot_dotplot(df, x_column, x_label, colorscale="Viridis", reverse_colorscale=True):
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)
    sizes = df[count_col] if count_col else pd.Series([10] * len(df), index=df.index)
    if count_col and sizes.max() > sizes.min():
        scaled_sizes = 8 + 22 * (sizes - sizes.min()) / (sizes.max() - sizes.min())
    else:
        scaled_sizes = pd.Series([16] * len(df), index=df.index)

    labels = df["Description"].astype(str) if "Description" in df.columns else df["ID"].astype(str)

    fig = go.Figure(data=go.Scatter(
        x=df[x_column], y=labels,
        mode="markers",
        marker=dict(
            size=scaled_sizes, sizemode="diameter",
            color=df["p.adjust"], colorscale=colorscale, reversescale=reverse_colorscale,
            colorbar=dict(title="Adjusted<br>p-value"),
            line=dict(width=1, color="DarkSlateGrey"),
        ),
        text=[
            f"{row.get('Description', row.get('ID', ''))}<br>"
            f"padj: {row['p.adjust']:.2e}<br>"
            f"{count_col or 'Count'}: {row.get(count_col, '—') if count_col else '—'}"
            for _, row in df.iterrows()
        ],
        hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title="",
        height=max(400, 28 * len(df)),
        margin=dict(l=10),
    )
    return fig


def _plot_barplot(df, x_column, x_label, colorscale="Blues"):
    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)
    labels = df["Description"].astype(str) if "Description" in df.columns else df["ID"].astype(str)

    fig = go.Figure(data=go.Bar(
        x=df[x_column], y=labels,
        orientation="h",
        marker=dict(
            color=df[count_col] if count_col else "#636EFA",
            colorscale=colorscale if count_col else None,
            colorbar=dict(title=count_col or "Count") if count_col else None,
        ),
        text=[f"padj: {p:.2e}" for p in df["p.adjust"]],
        hovertemplate="%{y}<br>%{text}<extra></extra>",
    ))
    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title="",
        height=max(400, 28 * len(df)),
        margin=dict(l=10),
    )
    return fig


def _plot_enrichment_map(result_df, n_top=30, similarity_threshold=0.2, label_map=None,
                          colorscale="Viridis", reverse_colorscale=True, layout_seed=42,
                          similarity_method="JC", orgdb_package=None, go_ontology=None):
    nodes_df, edges_df = om.build_term_similarity_network(
        result_df, n_top=n_top, similarity_threshold=similarity_threshold,
        similarity_method=similarity_method, orgdb_package=orgdb_package, go_ontology=go_ontology,
    )
    if nodes_df.empty or len(nodes_df) < 2:
        return None, []

    nodes_df = _apply_term_labels(nodes_df, label_map)

    positions = _compute_network_layout(nodes_df["ID"].tolist(), edges_df, seed=layout_seed)

    fig = go.Figure()

    if not edges_df.empty:
        edge_x, edge_y = [], []
        for _, row in edges_df.iterrows():
            x0, y0 = positions[row["source"]]
            x1, y1 = positions[row["target"]]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=1, color="rgba(150,150,150,0.5)"),
            hoverinfo="skip", showlegend=False,
        ))

    count_col = "Count" if "Count" in nodes_df.columns else None
    sizes = nodes_df[count_col] if count_col else pd.Series([15] * len(nodes_df))
    if count_col and sizes.max() > sizes.min():
        scaled_sizes = 15 + 35 * (sizes - sizes.min()) / (sizes.max() - sizes.min())
    else:
        scaled_sizes = pd.Series([25] * len(nodes_df))

    node_x = [positions[nid][0] for nid in nodes_df["ID"]]
    node_y = [positions[nid][1] for nid in nodes_df["ID"]]

    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode="markers",
        marker=dict(
            size=scaled_sizes, sizemode="diameter",
            color=nodes_df["p.adjust"], colorscale=colorscale, reversescale=reverse_colorscale,
            colorbar=dict(title="Adjusted<br>p-value"),
            line=dict(width=1, color="DarkSlateGrey"),
        ),
        hovertext=[
            f"{row['Description']}<br>padj: {row['p.adjust']:.2e}"
            + (f"<br>{count_col}: {row[count_col]}" if count_col else "")
            for _, row in nodes_df.iterrows()
        ],
        hovertemplate="%{hovertext}<extra></extra>",
        showlegend=False,
        name="Terms",
    ))

    label_annotations = _build_node_label_annotations(
        nodes_df["ID"].tolist(), positions,
        nodes_df["Description"].astype(str).str.slice(0, 30).tolist(),
    )
    fig.update_layout(annotations=label_annotations)

    fig.update_xaxes(visible=False, showgrid=False, zeroline=False)
    fig.update_yaxes(visible=False, showgrid=False, zeroline=False)
    fig.update_layout(height=650, plot_bgcolor="white")
    return fig, nodes_df["ID"].tolist()


def _plot_gene_concept_network(result_df, gene_fc_map=None, n_top=10, term_label_map=None,
                                gene_label_map=None, colorscale="RdBu", reverse_colorscale=True,
                                layout_seed=42):
    term_nodes_df, gene_nodes_df, edges_df = om.build_gene_concept_network(
        result_df, gene_fc_map=gene_fc_map, n_top=n_top,
    )
    if term_nodes_df.empty:
        return None, [], []

    term_nodes_df = _apply_term_labels(term_nodes_df, term_label_map)

    all_node_ids = term_nodes_df["ID"].tolist() + gene_nodes_df["gene_id"].tolist()
    positions = _compute_network_layout(all_node_ids, edges_df, seed=layout_seed)

    fig = go.Figure()

    if not edges_df.empty:
        edge_x, edge_y = [], []
        for _, row in edges_df.iterrows():
            x0, y0 = positions[row["source"]]
            x1, y1 = positions[row["target"]]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
        fig.add_trace(go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=1, color="rgba(150,150,150,0.4)"),
            hoverinfo="skip", showlegend=False,
        ))

    term_x = [positions[tid][0] for tid in term_nodes_df["ID"]]
    term_y = [positions[tid][1] for tid in term_nodes_df["ID"]]
    fig.add_trace(go.Scatter(
        x=term_x, y=term_y, mode="markers",
        marker=dict(size=32, color="#FDBF6F", line=dict(width=1.5, color="DarkSlateGrey")),
        hovertext=[f"{row['Description']}<br>padj: {row['p.adjust']:.2e}" for _, row in term_nodes_df.iterrows()],
        hovertemplate="%{hovertext}<extra></extra>",
        name="Enriched terms", showlegend=True,
    ))

    gene_x = [positions[gid][0] for gid in gene_nodes_df["gene_id"]]
    gene_y = [positions[gid][1] for gid in gene_nodes_df["gene_id"]]
    has_fc = gene_nodes_df["log2FoldChange"].notna().any()
    fig.add_trace(go.Scatter(
        x=gene_x, y=gene_y, mode="markers",
        marker=dict(
            size=14,
            color=gene_nodes_df["log2FoldChange"] if has_fc else "#AEC6CF",
            colorscale=colorscale if has_fc else None, reversescale=reverse_colorscale if has_fc else False,
            cmid=0 if has_fc else None,
            colorbar=dict(title="log2FC", x=1.1) if has_fc else None,
            line=dict(width=1, color="DarkSlateGrey"),
        ),
        hovertext=[
            f"{row['gene_id']}" + (f"<br>log2FC: {row['log2FoldChange']:.2f}" if pd.notna(row["log2FoldChange"]) else "")
            for _, row in gene_nodes_df.iterrows()
        ],
        hovertemplate="%{hovertext}<extra></extra>",
        name="Genes", showlegend=True,
    ))

    gene_label_map = gene_label_map or {}
    gene_display_labels = [gene_label_map.get(g, g) for g in gene_nodes_df["gene_id"]]

    term_label_annotations = _build_node_label_annotations(
        term_nodes_df["ID"].tolist(), positions,
        term_nodes_df["Description"].astype(str).str.slice(0, 28).tolist(),
        ax=25, ay=-20,
    )
    gene_label_annotations = _build_node_label_annotations(
        gene_nodes_df["gene_id"].tolist(), positions, gene_display_labels,
        font_size=9, ax=-15, ay=15,
    )
    fig.update_layout(annotations=term_label_annotations + gene_label_annotations)

    fig.update_xaxes(visible=False, showgrid=False, zeroline=False)
    fig.update_yaxes(visible=False, showgrid=False, zeroline=False)
    fig.update_layout(height=650, plot_bgcolor="white")
    return fig, term_nodes_df["ID"].tolist(), gene_nodes_df["gene_id"].tolist()


def _plot_gsea_running_score(running_score_df, term_label):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=running_score_df["x"], y=running_score_df["runningScore"],
        mode="lines", name="Running enrichment score",
        line=dict(color="#2CA02C", width=2),
        hovertemplate="Rank: %{x}<br>Running score: %{y:.3f}<extra></extra>",
    ))
    fig.add_hline(y=0, line_dash="dash", line_color="gray")

    member_ranks = running_score_df.loc[running_score_df["position"] == 1, "x"]
    y_min = running_score_df["runningScore"].min()
    tick_y0 = y_min - 0.05
    tick_y1 = y_min - 0.02
    for rank in member_ranks:
        fig.add_shape(
            type="line", x0=rank, x1=rank, y0=tick_y0, y1=tick_y1,
            line=dict(color="black", width=1),
        )

    fig.update_layout(
        xaxis_title="Rank in ordered gene list",
        yaxis_title="Running enrichment score",
        height=450,
    )
    return fig


def _plot_ridge(ridge_data, title_metric_label="log2 Fold Change"):
    import numpy as np

    if not ridge_data:
        return None

    fig = go.Figure()
    n = len(ridge_data)
    offset_step = 1.2

    x_min = min(min(d["values"]) for d in ridge_data)
    x_max = max(max(d["values"]) for d in ridge_data)
    x_grid = np.linspace(x_min - 0.5, x_max + 0.5, 200)

    for i, entry in enumerate(ridge_data):
        values = np.array(entry["values"])
        offset = i * offset_step

        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(values)
            density = kde(x_grid)
        except Exception:
            counts, edges = np.histogram(values, bins=15, range=(x_min - 0.5, x_max + 0.5))
            bin_centers = (edges[:-1] + edges[1:]) / 2
            density = np.interp(x_grid, bin_centers, counts, left=0, right=0)

        if density.max() > 0:
            density = density / density.max() * offset_step * 0.85

        color = "#d62728" if (entry.get("NES") or 0) > 0 else "#1f77b4"

        fig.add_trace(go.Scatter(
            x=x_grid, y=density + offset,
            mode="lines", fill="tonexty" if i == 0 else "tozeroy",
            line=dict(color=color, width=1.5),
            fillcolor=color.replace(")", ", 0.4)").replace("#d62728", "rgba(214,39,40,0.4)").replace("#1f77b4", "rgba(31,119,180,0.4)"),
            name=str(entry["description"])[:40],
            hovertemplate=f"{entry['description']}<br>NES: {entry.get('NES', 'N/A')}<extra></extra>",
        ))
        fig.add_hline(y=offset, line_width=0.5, line_color="lightgray")

    fig.update_layout(
        xaxis_title=title_metric_label,
        yaxis=dict(
            tickmode="array",
            tickvals=[i * offset_step for i in range(n)],
            ticktext=[str(d["description"])[:35] for d in ridge_data],
        ),
        height=max(400, 60 * n),
        showlegend=False,
        margin=dict(l=200),
    )
    return fig


def _plot_upset(upset_df, n_show=15, bar_color="#4C72B0"):
    if upset_df is None or upset_df.empty:
        return None

    df = upset_df.head(n_show).iloc[::-1]

    fig = go.Figure(data=go.Bar(
        x=df["size"], y=df["combination"],
        orientation="h",
        marker=dict(color=bar_color),
        hovertemplate="%{y}<br>Genes: %{x}<extra></extra>",
    ))
    fig.update_layout(
        xaxis_title="Number of genes",
        yaxis_title="",
        height=max(350, 30 * len(df)),
        margin=dict(l=10),
    )
    return fig


def _plot_compare_cluster_dotplot(cc_df, n_top_per_cluster=10, label_map=None,
                                   colorscale="Viridis", reverse_colorscale=True):
    if cc_df is None or cc_df.empty:
        return None

    cc_df = _apply_term_labels(cc_df, label_map)

    top_ids_per_cluster = (
        cc_df.sort_values("p.adjust").groupby("Cluster").head(n_top_per_cluster)["ID"].unique()
    )
    df = cc_df[cc_df["ID"].isin(top_ids_per_cluster)].copy()

    id_to_desc = dict(zip(df["ID"], df.get("Description", df["ID"])))
    term_order = df.groupby("ID")["p.adjust"].min().sort_values().index.tolist()
    term_order = term_order[::-1]

    count_col = "Count" if "Count" in df.columns else None
    sizes = df[count_col] if count_col else pd.Series([15] * len(df))
    if count_col and sizes.max() > sizes.min():
        scaled_sizes = 8 + 22 * (sizes - sizes.min()) / (sizes.max() - sizes.min())
    else:
        scaled_sizes = pd.Series([16] * len(df), index=df.index)

    fig = go.Figure(data=go.Scatter(
        x=df["Cluster"], y=df["ID"].map(id_to_desc),
        mode="markers",
        marker=dict(
            size=scaled_sizes, sizemode="diameter",
            color=df["p.adjust"], colorscale=colorscale, reversescale=reverse_colorscale,
            colorbar=dict(title="Adjusted<br>p-value"),
            line=dict(width=1, color="DarkSlateGrey"),
        ),
        text=[
            f"{row.get('Description', row['ID'])}<br>Cluster: {row['Cluster']}<br>padj: {row['p.adjust']:.2e}"
            for _, row in df.iterrows()
        ],
        hovertemplate="%{text}<extra></extra>",
    ))
    fig.update_layout(
        xaxis_title="Contrast",
        yaxis=dict(categoryorder="array", categoryarray=[id_to_desc[t] for t in term_order]),
        yaxis_title="",
        height=max(400, 28 * len(term_order)),
        margin=dict(l=10),
    )
    return fig


def _render_category_filter_and_label_controls(combined_df, key_prefix):
    if combined_df is None or combined_df.empty or "Category" not in combined_df.columns:
        return set(), {}

    available_categories = sorted(combined_df["Category"].unique())

    st.markdown("**Which categories should be included?**")
    st.caption(
        "Unchecking a category removes its terms from the plot "
        "entirely (not just its color) -- the plot adjusts "
        "automatically to only the categories you keep checked."
    )
    selected_categories = set()
    cols = st.columns(min(len(available_categories), 5))
    for i, cat in enumerate(available_categories):
        with cols[i % len(cols)]:
            if st.checkbox(cat, value=True, key=f"{key_prefix}_category_checkbox_{cat}"):
                selected_categories.add(cat)

    category_label_map = {}
    with st.expander("✏️ Rename category labels shown in the legend (optional)"):
        st.caption(
            "Edit the \"Custom Label\" column below to rename how each "
            "category appears in the plot's legend (e.g. \"BP\" -> "
            "\"Biological Process\") -- this only changes what's "
            "DISPLAYED, it does not affect which terms are included."
        )
        editor_df = pd.DataFrame({"Category": available_categories})
        editor_df["Custom Label"] = editor_df["Category"]

        edited = st.data_editor(
            editor_df, use_container_width=True, hide_index=True,
            disabled=["Category"], key=f"{key_prefix}_category_label_editor",
        )
        for _, row in edited.iterrows():
            custom = str(row["Custom Label"]).strip()
            original = str(row["Category"]).strip()
            if custom and custom != original:
                category_label_map[row["Category"]] = custom

    return selected_categories, category_label_map


def _plot_combined_multi_database(combined_df, label_map=None, category_label_map=None,
                                   colorscale="Viridis", reverse_colorscale=True):
    import math

    if combined_df is None or combined_df.empty:
        return None, {}

    combined_df = _apply_term_labels(combined_df, label_map)
    category_label_map = category_label_map or {}

    df = combined_df.dropna(subset=["p.adjust"]).sort_values("p.adjust").copy()
    df["neg_log10_padj"] = -df["p.adjust"].clip(lower=1e-300).apply(math.log10)
    df = df.iloc[::-1].reset_index(drop=True)

    count_col = "Count" if "Count" in df.columns else ("setSize" if "setSize" in df.columns else None)
    sizes = df[count_col] if count_col else pd.Series([15] * len(df))
    if count_col and sizes.max() > sizes.min():
        scaled_sizes = 8 + 22 * (sizes - sizes.min()) / (sizes.max() - sizes.min())
    else:
        scaled_sizes = pd.Series([16] * len(df), index=df.index)

    labels = df["Description"].astype(str) if "Description" in df.columns else df["ID"].astype(str)

    unique_categories = sorted(df["Category"].unique())
    category_colors = {
        cat: _hex_to_rgba(_CATEGORY_BAND_PALETTE[i % len(_CATEGORY_BAND_PALETTE)], alpha=0.18)
        for i, cat in enumerate(unique_categories)
    }

    def _display_category(cat):
        return category_label_map.get(cat, cat)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["neg_log10_padj"], y=labels,
        mode="markers",
        marker=dict(
            size=scaled_sizes, sizemode="diameter",
            color=df["p.adjust"], colorscale=colorscale, reversescale=reverse_colorscale,
            colorbar=dict(title="Adjusted<br>p-value"),
            line=dict(width=1, color="DarkSlateGrey"),
        ),
        text=[
            f"{row['Database']} / {_display_category(row['Category'])}<br>padj: {row['p.adjust']:.2e}"
            for _, row in df.iterrows()
        ],
        hovertemplate="%{y}<br>%{text}<extra></extra>",
        showlegend=False,
    ))

    shapes = [
        dict(
            type="rect", xref="paper", yref="y",
            x0=0, x1=1, y0=i - 0.5, y1=i + 0.5,
            fillcolor=category_colors[row["Category"]],
            line=dict(width=0), layer="below",
        )
        for i, (_, row) in enumerate(df.iterrows())
    ]
    fig.update_layout(shapes=shapes)

    displayed_category_colors = {}
    for cat in unique_categories:
        display_name = _display_category(cat)
        displayed_category_colors[display_name] = category_colors[cat]
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=14, color=category_colors[cat], symbol="square"),
            name=display_name,
            showlegend=True,
            hoverinfo="skip",
        ))

    fig.update_layout(
        xaxis_title="-log10(adjusted p-value)",
        yaxis_title="",
        height=max(400, 28 * len(df)),
        margin=dict(l=10),
        legend=dict(title="Category"),
    )
    return fig, displayed_category_colors


# ---------------------------------------------------------------------------
# Small shared UI helpers
# ---------------------------------------------------------------------------

def _species_and_orgdb_for_project(project):
    override_key = pm.get_ontology_species_override(project)
    if override_key:
        orgdb_package = gim.ORGDB_PACKAGES.get(override_key)
        species_label = rm.REFERENCE_CATALOG.get(override_key, {}).get("label", override_key)
        return override_key, orgdb_package, species_label

    species_key, is_custom = pm.get_reference_choice(project)
    if is_custom or not species_key:
        return None, None, None
    orgdb_package = gim.ORGDB_PACKAGES.get(species_key)
    species_label = rm.REFERENCE_CATALOG.get(species_key, {}).get("label", species_key)
    return species_key, orgdb_package, species_label


def _render_species_override_picker(project):
    st.warning(
        "⚠️ This project's reference was provided as a custom upload, "
        "so we can't automatically determine which Bioconductor "
        "annotation package to use for GO/KEGG/Reactome enrichment. "
        "If your data is actually one of this app's supported preset "
        "organisms, select it below to continue -- ontology enrichment "
        "only needs to know the correct SPECIES, not the actual "
        "reference files themselves."
    )

    species_choices = gim.orgdb_species_choices(rm)
    species_keys = list(species_choices.keys())

    previous_override = pm.get_ontology_species_override(project)
    default_index = species_keys.index(previous_override) if previous_override in species_keys else 0

    picked_key = st.selectbox(
        "Which organism is this data actually from?",
        options=species_keys,
        index=default_index,
        format_func=lambda k: species_choices[k],
        key="ontology_species_override_select",
    )

    if st.button("✅ Use This Organism for Ontology Analysis", key="ontology_species_override_confirm_btn"):
        pm.save_ontology_species_override(project, picked_key)
        st.rerun()

    if previous_override in species_keys:
        orgdb_package = gim.ORGDB_PACKAGES.get(previous_override)
        species_label = species_choices.get(previous_override, previous_override)
        return previous_override, orgdb_package, species_label

    return None, None, None


def _render_plot_checkbox_gate(available_plots, key_prefix):
    st.markdown("**Which plot(s) would you like to see?**")
    selected = set()
    n_cols = 2
    cols = st.columns(n_cols)
    for i, (plot_key, label, help_text) in enumerate(available_plots):
        with cols[i % n_cols]:
            if st.checkbox(label, key=f"{key_prefix}_show_{plot_key}", help=help_text):
                selected.add(plot_key)
    return selected


def _get_gene_fc_map_for_contrast(deseq2_out_dir, contrast_name, gene_name_map=None):
    results_df = dm.read_contrast_results(deseq2_out_dir, contrast_name)
    if results_df is None:
        return {}
    gene_name_map = gene_name_map or {}
    fc_map = {}
    for _, row in results_df.dropna(subset=["log2FoldChange"]).iterrows():
        gid = str(row["gene_id"])
        display_name = gene_name_map.get(gid, gid)
        fc_map[display_name] = row["log2FoldChange"]
        fc_map[gid] = row["log2FoldChange"]
    return fc_map


def _render_gene_set_size_controls(key_prefix):
    with st.expander("ℹ️ What are 'gene set size' limits, and why set them? (click to learn more)"):
        st.markdown(
            "Every GO term/KEGG pathway/Reactome pathway is really just "
            "a predefined LIST of genes believed to share some "
            "biological role. These two settings control which of "
            "those lists are even considered for testing:\n\n"
            "- **Minimum gene set size** -- excludes very SMALL gene "
            "sets (e.g. just 1-2 genes), which produce a statistically "
            "unstable, unreliable result.\n"
            "- **Maximum gene set size** -- excludes very LARGE, "
            "overly generic gene sets (e.g. thousands of genes), which "
            "are frequently uninformative and can dominate results "
            "without saying anything specific about your experiment.\n\n"
            "The defaults below (10-500 genes) match clusterProfiler's "
            "own standard defaults."
        )
    col1, col2 = st.columns(2)
    with col1:
        min_gs_size = st.number_input(
            "Minimum gene set size:", min_value=1, max_value=1000,
            value=om.DEFAULT_MIN_GS_SIZE, step=1,
            key=f"{key_prefix}_min_gs_size",
            help="Raising this excludes small, statistically unstable gene sets from being tested at all -- fewer total terms tested, but each result is more reliable.",
        )
    with col2:
        max_gs_size = st.number_input(
            "Maximum gene set size:", min_value=1, max_value=10000,
            value=om.DEFAULT_MAX_GS_SIZE, step=10,
            key=f"{key_prefix}_max_gs_size",
            help="Lowering this excludes very large, generic gene sets from being tested -- helps prevent a handful of huge, uninformative categories from dominating your results.",
        )
    if min_gs_size >= max_gs_size:
        st.warning("⚠️ Minimum gene set size must be smaller than the maximum. Using the defaults (10-500) instead.")
        return om.DEFAULT_MIN_GS_SIZE, om.DEFAULT_MAX_GS_SIZE
    return int(min_gs_size), int(max_gs_size)


def _render_direction_split_control(key_prefix):
    with st.expander("ℹ️ Why analyze up- and down-regulated genes separately? (click to learn more)"):
        st.markdown(
            "By default, this workspace runs ORA **twice** for each "
            "database you select: once using only your UP-regulated "
            "significant genes, and once using only your DOWN-"
            "regulated significant genes -- rather than lumping both "
            "directions into a single mixed gene list.\n\n"
            "**Why this matters:** if up- and down-regulated genes are "
            "combined into one list, a GO term/pathway can get flagged "
            "as \"enriched\" purely because it happens to contain a mix "
            "of genes going in both directions -- without that term "
            "actually representing one coherent, DIRECTIONAL biological "
            "signal. For example, a pathway that's genuinely "
            "**activated** and a completely different pathway that's "
            "genuinely **suppressed** could get blurred together into "
            "one ambiguous \"changed somehow\" result. Analyzing each "
            "direction separately gives you a cleaner, more "
            "interpretable answer: \"here's what's enriched among genes "
            "going UP\" and, independently, \"here's what's enriched "
            "among genes going DOWN.\"\n\n"
            "You can still directly compare the two directions "
            "side-by-side once both are run (see Step 4 below). "
            "Unchecking this reverts to the older, combined-list "
            "behavior, offered here for anyone who specifically wants "
            "it or needs to compare against a combined-list result."
        )
    return st.checkbox(
        "🔀 Analyze up- and down-regulated genes separately (recommended)",
        value=True,
        key=f"{key_prefix}_split_direction",
    )


def _render_qvalue_threshold_control(key_prefix):
    use_qvalue = st.checkbox(
        "Also require a q-value below a threshold (optional, stricter)",
        value=False,
        key=f"{key_prefix}_use_qvalue_filter",
        help=(
            "q-value is a DIFFERENT multiple-testing correction than "
            "the adjusted p-value (padj) used above. Checking this "
            "adds q-value as an ADDITIONAL requirement on top of padj."
        ),
    )
    if not use_qvalue:
        return None
    return st.slider(
        "Q-value significance threshold:", 0.01, 0.20, 0.05, step=0.01,
        key=f"{key_prefix}_qvalue_threshold",
        help="Lowering this makes the significance summary stricter (fewer terms count as significant); raising it is more lenient (more terms count). This only affects the significance SUMMARY caption -- it does not remove any rows from the underlying result table or plots.",
    )


def _render_simplify_control(key_prefix, go_ontology):
    """
    Render the "simplify GO results" checkbox, ALONGSIDE a proactive
    GOSemSim availability check, a similarity-measure picker (WHICH
    GOSemSim measure decides which terms count as redundant), and a
    similarity-cutoff slider (HOW similar two terms must be before
    they're collapsed together).

    Returns (simplify_go: bool, simplify_measure: str, simplify_cutoff:
    float) -- simplify_measure/simplify_cutoff default to
    om.DEFAULT_SIMPLIFY_MEASURE / om.DEFAULT_SIMPLIFY_CUTOFF and are
    only meaningful when simplify_go is True, but are always returned
    so callers can unconditionally unpack the tuple.
    """
    available_ontology = om.simplify_available_for_ontology(go_ontology)
    with st.expander("ℹ️ What does 'simplify GO results' do? (click to learn more)"):
        st.markdown(
            "GO terms are organized in a hierarchy, so a real result "
            "often contains many terms that are highly REDUNDANT with "
            "each other -- e.g. \"immune response\" and \"regulation of "
            "immune response\" both showing up separately, when really "
            "they're describing the same underlying signal. Checking "
            "this option collapses semantically similar terms down to "
            "just the single most significant representative from each "
            "group, using:\n\n"
            "- **which similarity measure** decides whether two terms "
            "count as \"similar enough\" (picker below), and\n"
            "- **how similar is similar enough** (cutoff slider below).\n\n"
            "**Only available for a single GO sub-ontology at a time "
            "(BP, MF, or CC)** -- not available when \"All three (BP + "
            "MF + CC)\" is selected, since there's no single similarity "
            "measure that meaningfully compares terms across different "
            "sub-ontologies."
        )

    gosemsim_ok, gosemsim_detail = om.check_gosemsim_available()

    col_checkbox, col_status = st.columns([3, 2])
    with col_checkbox:
        simplify_go = st.checkbox(
            "🧹 Simplify GO results (remove redundant/highly similar terms)",
            value=False,
            disabled=not available_ontology,
            key=f"{key_prefix}_simplify_go",
            help="Requires the GOSemSim package -- if missing, the analysis still completes successfully using un-simplified results." if available_ontology else "Only available when a single GO sub-ontology (BP, MF, or CC) is selected.",
        )
    with col_status:
        if available_ontology:
            if gosemsim_ok:
                st.caption("✅ GOSemSim detected")
            else:
                st.caption("⚠️ GOSemSim not detected")

    simplify_measure = om.DEFAULT_SIMPLIFY_MEASURE
    simplify_cutoff = om.DEFAULT_SIMPLIFY_CUTOFF
    if available_ontology:
        measure_display_names = list(om.SIMPLIFY_MEASURE_OPTIONS.keys())
        default_measure_display = next(
            (k for k, v in om.SIMPLIFY_MEASURE_OPTIONS.items() if v == om.DEFAULT_SIMPLIFY_MEASURE),
            measure_display_names[0],
        )
        chosen_measure_display = st.selectbox(
            "Similarity measure used to detect redundant GO terms:",
            options=measure_display_names,
            index=measure_display_names.index(default_measure_display),
            key=f"{key_prefix}_simplify_measure",
            help="Only affects the result if the checkbox above is checked. Different measures can merge a different number/set of terms -- see the explanation panel below for what to expect from each.",
        )
        simplify_measure = om.SIMPLIFY_MEASURE_OPTIONS[chosen_measure_display]

        measure_info = om.SIMILARITY_MEASURE_INFO[simplify_measure]
        with st.expander(f"ℹ️ About \"{chosen_measure_display}\" ({measure_info['category']})"):
            st.latex(measure_info["formula"])
            st.markdown(f"**What it computes:** {measure_info['summary']}")
            st.markdown(f"**Background:** {measure_info['details']}")
            st.markdown(f"**What to expect in YOUR results:** {measure_info['practical_effect']}")
            st.caption(
                "📊 Requires computing Information Content over the GO corpus first (done automatically)"
                if measure_info["requires_compute_ic"]
                else "⚡ No Information Content computation needed -- uses only GO graph structure"
            )
            st.caption(f"Reference: {measure_info['citation']}")

        simplify_cutoff = st.slider(
            "Similarity cutoff (how similar is 'redundant'?):",
            min_value=0.1, max_value=1.0, value=om.DEFAULT_SIMPLIFY_CUTOFF, step=0.05,
            key=f"{key_prefix}_simplify_cutoff",
            help=(
                "Two GO terms are collapsed into one (keeping only the more "
                "significant of the pair) when their similarity score, using "
                "the measure chosen above, is at or above this value.\n\n"
                "• LOWER cutoff (e.g. 0.4) = STRICTER redundancy check -- more "
                "terms get merged/removed, leaving a shorter, more curated list.\n\n"
                "• HIGHER cutoff (e.g. 0.9) = LOOSER redundancy check -- only "
                "near-identical terms get merged, leaving more terms in the "
                "final result.\n\n"
                "clusterProfiler's own standard default is 0.7."
            ),
        )

    if available_ontology and not gosemsim_ok:
        st.warning(
            f"⚠️ {gosemsim_detail} Checking the box above will NOT "
            "actually simplify your results -- the analysis will "
            "still complete successfully, just using the un-simplified "
            "GO terms (this is expected, graceful behavior, not an "
            "error). Install GOSemSim in this environment to enable "
            "simplification, then re-check below."
        )

    if available_ontology:
        with st.expander("🔬 Run a deeper GOSemSim check (optional, slower)"):
            st.caption(
                "The status above only confirms whether the GOSemSim "
                "package is *installed*. This runs a genuinely deeper "
                "test -- actually simplifying a small real result for "
                "your organism -- to catch a different, rarer failure "
                "mode: GOSemSim installed but unable to build its "
                "semantic-similarity data for this specific organism/"
                "ontology. This can take anywhere from a few seconds to "
                "over a minute (longer the first time, since GOSemSim "
                "caches this data after building it once)."
            )
            col_recheck, col_functional = st.columns(2)
            with col_recheck:
                if st.button("🔄 Re-check availability", key=f"{key_prefix}_recheck_gosemsim_btn"):
                    om.check_gosemsim_available(force_recheck=True)
                    st.rerun()
            with col_functional:
                run_functional_check = st.button(
                    "🧪 Run live functional test", key=f"{key_prefix}_run_functional_check_btn",
                )
            if run_functional_check:
                orgdb_package = st.session_state.get("_ontology_orgdb_package_for_gosemsim_check")
                test_ontology = go_ontology if go_ontology != "ALL" else "BP"
                if not orgdb_package:
                    st.error("⚠️ Could not determine which organism to test against -- select an organism above first.")
                else:
                    with st.spinner(f"Running a live GOSemSim test for {orgdb_package} ({test_ontology})..."):
                        functional, detail = om.verify_gosemsim_functional(orgdb_package, test_ontology)
                    if functional:
                        st.success(f"✅ {detail}")
                    else:
                        st.error(f"⚠️ {detail}")

    return simplify_go, simplify_measure, simplify_cutoff


def _render_similarity_method_control(key_prefix, database):
    """
    Render the "how should term similarity be measured?" control for the
    enrichment map plot -- an INDEPENDENT setting from the "simplify GO
    results" measure above (a user might simplify with Wang but still
    view the map with the original gene-overlap similarity, or vice
    versa).

    For GO results, offers the full choice between gene-overlap (JC --
    this workspace's original, unchanged default) and GOSemSim's five
    semantic-similarity measures (Resnik/Lin/Rel/Jiang/Wang). For KEGG/
    Reactome results, GOSemSim's measures don't apply (no GO graph
    structure to use), so no dropdown is rendered at all -- gene-overlap
    similarity is used silently, exactly as this workspace has always
    done for those databases.

    Returns the chosen measure key (e.g. "JC", "Wang", "Rel", ...).
    """
    if database != "GO":
        return "JC"

    display_names = list(om.SIMILARITY_METHOD_OPTIONS.keys())
    default_display = next(
        (k for k, v in om.SIMILARITY_METHOD_OPTIONS.items() if v == om.DEFAULT_SIMILARITY_METHOD),
        display_names[0],
    )
    chosen_display = st.selectbox(
        "How should term similarity be measured?",
        options=display_names,
        index=display_names.index(default_display),
        key=f"{key_prefix}_similarity_method",
        help="Determines which terms are drawn close together and connected in the network below. Changing this can noticeably change the plot's layout and which connections appear -- see the explanation panel for what to expect.",
    )
    method = om.SIMILARITY_METHOD_OPTIONS[chosen_display]

    info = om.SIMILARITY_MEASURE_INFO[method]
    with st.expander(f"ℹ️ About \"{chosen_display}\" ({info['category']})"):
        st.latex(info["formula"])
        st.markdown(f"**What it computes:** {info['summary']}")
        st.markdown(f"**Background:** {info['details']}")
        st.markdown(f"**What to expect in YOUR plot:** {info['practical_effect']}")
        st.caption(
            "📊 Requires building GO semantic-similarity data (GOSemSim)"
            if method != "JC"
            else "⚡ No GOSemSim/R call needed -- computed directly from each term's gene list"
        )
        st.caption(f"Reference: {info['citation']}")

    if method != "JC":
        st.caption(
            "If GOSemSim isn't installed/available (or the computation "
            "fails for any reason), this plot will automatically fall "
            "back to gene-overlap (Jaccard) similarity instead."
        )

    return method


def _render_simplify_status(output_dir, direction=None):
    outcomes = om.load_simplify_status(output_dir)
    if not outcomes:
        return

    if direction is not None:
        outcomes = [o for o in outcomes if o.get("direction") == direction]
        if not outcomes:
            return

    n_simplified = sum(1 for o in outcomes if o["outcome"] == "simplified")
    n_fallback = sum(1 for o in outcomes if o["outcome"] == "fallback")

    if n_fallback == 0:
        st.success(
            f"✅ GO term simplification was applied successfully "
            f"({n_simplified} of {len(outcomes)} GO result(s) simplified)."
        )
    elif n_simplified == 0:
        st.warning(
            f"⚠️ GO term simplification did **not** run for this result "
            f"({n_fallback} of {len(outcomes)} GO result(s)) -- GOSemSim "
            "may not be installed in this environment. The GO results "
            "shown below are the full, UN-simplified set (this is "
            "expected, graceful behavior -- the analysis itself "
            "completed successfully either way)."
        )
        with st.expander("Show technical detail"):
            for o in outcomes:
                st.code(o["detail"])
    else:
        st.warning(
            f"⚠️ GO term simplification partially applied: "
            f"{n_simplified} of {len(outcomes)} GO result(s) were "
            f"simplified, but {n_fallback} fell back to un-simplified "
            "results (see technical detail below)."
        )
        with st.expander("Show technical detail"):
            for o in outcomes:
                st.code(f"[{o['outcome']}] {o['detail']}")


# ---------------------------------------------------------------------------
# Result rendering: a single database's ORA or GSEA result
# ---------------------------------------------------------------------------

def _render_single_database_results(output_dir, analysis_type, database, key_prefix,
                                     padj_threshold, gene_fc_map=None, direction=None,
                                     section_label=None, orgdb_package=None, go_ontology=None):
    """
    orgdb_package, go_ontology: forwarded from render()'s Step 2 state --
        used only by the enrichment-map section's (new) similarity-method
        picker, when database == "GO" and a GOSemSim measure (rather than
        the default gene-overlap "JC") is selected.
    """
    read_direction = direction if analysis_type == "ora" else None
    result_df = om.read_enrichment_result(output_dir, analysis_type, database, direction=read_direction)
    if result_df is None or result_df.empty:
        st.info(f"No {database} results were found (or zero terms were returned).")
        return

    st.markdown(f"#### {section_label or f'{database} results'}")

    if database == "GO":
        _render_simplify_status(output_dir, direction=direction)

    qvalue_threshold = _render_qvalue_threshold_control(f"{key_prefix}_qval")
    st.caption(om.summarize_enrichment_result(result_df, padj_threshold=padj_threshold, qvalue_threshold=qvalue_threshold))

    is_gsea = analysis_type == "gsea"

    available_plots = [
        ("dotplot", "🔵 Dot plot", "The classic clusterProfiler summary view: one row per term, dot size = gene count, dot color = significance. Good first stop for any result."),
        ("barplot", "📊 Bar plot", "Bars ranked by your chosen sort/x-axis value -- a simpler alternative to the dot plot when gene-count nuance matters less."),
        ("emap", "🕸️ Enrichment map (term network)", "Connects enriched terms that share many of the same genes, so clusters of related/redundant terms are visible at a glance instead of scrolling a long list of similarly-worded terms."),
        ("cnet", "🔗 Gene-concept network", "Directly connects enriched terms to the individual genes driving them -- shows which specific genes matter, and whether any gene is shared across multiple terms."),
    ]
    if is_gsea:
        available_plots.append(
            ("running_score", "📈 GSEA running-score plot", "The classic GSEA enrichment curve for one specific gene set at a time -- shows exactly where that gene set's members fall along your full ranked gene list."),
        )
        available_plots.append(
            ("ridge", "🏔️ Ridge plot", "Shows the distribution of fold-change values among each top gene set's own member genes, stacked as ridges -- lets you compare which gene sets skew strongly up/down vs. which are more centered."),
        )
    available_plots.append(
        ("upset", "📶 UpSet plot (term overlap)", "Shows exactly how many genes belong to each specific COMBINATION of enriched terms -- the standard alternative to a Venn diagram once you have more than 2-3 terms to compare."),
    )

    selected_plots = _render_plot_checkbox_gate(available_plots, key_prefix)

    if not selected_plots:
        st.caption("Check one or more boxes above to display the corresponding plot(s).")
        return

    if "dotplot" in selected_plots:
        st.markdown("##### Dot plot")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Each row is one enriched term, ordered by whatever "
                "you choose to rank/sort by below. **Dot size** "
                "reflects how many genes support that term. **Dot "
                "color** shows the adjusted p-value -- pick a "
                "different color scale below if you prefer."
            )
        dot_prepared_df, dot_x_column, dot_x_label = _render_ranked_plot_controls(
            result_df, is_gsea, f"{key_prefix}_dotplot",
        )
        dot_label_map = _render_term_label_editor(dot_prepared_df, f"{key_prefix}_dotplot_labels")
        dot_labeled_df = _apply_term_labels(dot_prepared_df, dot_label_map)
        dot_colorscale, dot_reverse = dew._render_heatmap_color_controls(
            f"{key_prefix}_dotplot", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
        )
        fig = _plot_dotplot(dot_labeled_df, dot_x_column, dot_x_label, colorscale=dot_colorscale, reverse_colorscale=dot_reverse)
        style = dew._render_plot_style_controls(f"{key_prefix}_dotplot", group_values=None, show_legend_controls=False)
        dew._apply_plot_style(fig, style, default_title=f"{database} {analysis_type.upper()} -- Dot Plot")
        dew._render_plotly_chart(fig)
        dew._render_pdf_export(fig, f"{key_prefix}_dotplot", f"{database}_{analysis_type}_dotplot")

    if "barplot" in selected_plots:
        st.markdown("##### Bar plot")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Bars show whatever value you choose below for each "
                "term, ordered by your chosen ranking. Bar color "
                "reflects how many genes support that term."
            )
        bar_prepared_df, bar_x_column, bar_x_label = _render_ranked_plot_controls(
            result_df, is_gsea, f"{key_prefix}_barplot",
        )
        bar_label_map = _render_term_label_editor(bar_prepared_df, f"{key_prefix}_barplot_labels")
        bar_labeled_df = _apply_term_labels(bar_prepared_df, bar_label_map)
        bar_colorscale, _ = dew._render_heatmap_color_controls(
            f"{key_prefix}_barplot", SEQUENTIAL_COLORSCALE_OPTIONS, default="Blues", default_reverse=False,
        )
        fig = _plot_barplot(bar_labeled_df, bar_x_column, bar_x_label, colorscale=bar_colorscale)
        style = dew._render_plot_style_controls(f"{key_prefix}_barplot", group_values=None, show_legend_controls=False)
        dew._apply_plot_style(fig, style, default_title=f"{database} {analysis_type.upper()} -- Bar Plot")
        dew._render_plotly_chart(fig)
        dew._render_pdf_export(fig, f"{key_prefix}_barplot", f"{database}_{analysis_type}_barplot")

    if "emap" in selected_plots:
        st.markdown("##### Enrichment map")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Each node is one enriched term. Terms that are similar "
                "to each other (however similarity is defined below) are "
                "drawn close together and connected by a line. Node size "
                "= gene count; node color = significance.\n\n"
                "💡 **Labels are individually draggable** -- click and "
                "drag any term's label to reposition it if labels "
                "overlap. The node itself stays fixed at its computed "
                "position -- only the label moves, with its arrow "
                "always pointing back at the correct node. Try the "
                "\"different layout\" button below for an alternative "
                "automatic arrangement."
            )
        n_top_emap = st.slider(
            "Number of top terms to include:", min_value=5, max_value=50, value=30, step=5,
            key=f"{key_prefix}_emap_n_top",
            help="How many of the most significant terms are included as nodes. More terms = a bigger, potentially busier network; fewer terms = a smaller, more focused one.",
        )
        emap_colorscale, emap_reverse = dew._render_heatmap_color_controls(
            f"{key_prefix}_emap", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
        )
        emap_similarity_method = _render_similarity_method_control(f"{key_prefix}_emap", database)
        emap_threshold_label = (
            "Minimum semantic similarity to draw a connection:"
            if emap_similarity_method != "JC"
            else "Minimum gene overlap (Jaccard similarity) to draw a connection:"
        )
        similarity_threshold = st.slider(
            emap_threshold_label,
            0.0, 1.0, 0.2, step=0.05, key=f"{key_prefix}_emap_threshold",
            help=(
                "Two terms are only connected by a line if their similarity "
                "score (using the method chosen above) is at or above this "
                "value.\n\n"
                "• LOWER threshold (e.g. 0.05) = MORE connections drawn -- "
                "denser network, easier to spot broad clusters, but can look "
                "cluttered.\n\n"
                "• HIGHER threshold (e.g. 0.5) = FEWER connections drawn -- "
                "sparser network showing only the strongest relationships; "
                "some terms may end up with no connections at all (isolated "
                "nodes) if nothing else is similar enough."
            ),
        )
        emap_seed = _render_layout_seed_control(f"{key_prefix}_emap")
        fig, emap_node_ids = _plot_enrichment_map(
            result_df, n_top=n_top_emap, similarity_threshold=similarity_threshold,
            colorscale=emap_colorscale, reverse_colorscale=emap_reverse, layout_seed=emap_seed,
            similarity_method=emap_similarity_method, orgdb_package=orgdb_package, go_ontology=go_ontology,
        )
        if fig is None:
            st.info("Not enough overlapping terms to draw a meaningful network at this threshold -- try lowering the similarity slider above.")
        else:
            emap_label_map = _render_term_label_editor(
                result_df[result_df["ID"].isin(emap_node_ids)], f"{key_prefix}_emap_labels",
            )
            if emap_label_map:
                fig, _ = _plot_enrichment_map(
                    result_df, n_top=n_top_emap, similarity_threshold=similarity_threshold,
                    label_map=emap_label_map, colorscale=emap_colorscale,
                    reverse_colorscale=emap_reverse, layout_seed=emap_seed,
                    similarity_method=emap_similarity_method, orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
            style = dew._render_plot_style_controls(f"{key_prefix}_emap", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(fig, style, default_title=f"{database} {analysis_type.upper()} -- Enrichment Map")
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{key_prefix}_emap", f"{database}_{analysis_type}_enrichment_map")

    if "cnet" in selected_plots:
        st.markdown("##### Gene-concept network")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Larger orange nodes are enriched terms; smaller nodes "
                "are the individual genes that belong to them, colored "
                "by their own fold change when available.\n\n"
                "💡 **Both term labels AND gene labels are individually "
                "draggable** -- drag any label to declutter overlapping "
                "text; its arrow will keep pointing at the correct "
                "node. Try the \"different layout\" button below for an "
                "alternative automatic arrangement."
            )
        cnet_colorscale, cnet_reverse = dew._render_heatmap_color_controls(
            f"{key_prefix}_cnet", DIVERGING_COLORSCALE_OPTIONS, default="RdBu", default_reverse=True,
        )
        n_cnet_terms = st.slider(
            "Number of top terms to include:", min_value=3, max_value=20, value=8, step=1,
            key=f"{key_prefix}_cnet_n_top",
            help="How many of the most significant terms are included. More terms pull in more genes too (since every gene belonging to an included term is shown), which can make the network larger and busier.",
        )
        cnet_seed = _render_layout_seed_control(f"{key_prefix}_cnet")
        fig, cnet_term_ids, cnet_gene_ids = _plot_gene_concept_network(
            result_df, gene_fc_map=gene_fc_map, n_top=n_cnet_terms,
            colorscale=cnet_colorscale, reverse_colorscale=cnet_reverse, layout_seed=cnet_seed,
        )
        if fig is None:
            st.info("No data available to build this network.")
        else:
            cnet_term_label_map = _render_term_label_editor(
                result_df[result_df["ID"].isin(cnet_term_ids)], f"{key_prefix}_cnet_labels",
                label="✏️ Shorten/rename TERM labels shown on this plot (optional)",
            )
            cnet_gene_label_map = _render_gene_label_editor(cnet_gene_ids, f"{key_prefix}_cnet")
            if cnet_term_label_map or cnet_gene_label_map:
                fig, _, _ = _plot_gene_concept_network(
                    result_df, gene_fc_map=gene_fc_map, n_top=n_cnet_terms,
                    term_label_map=cnet_term_label_map, gene_label_map=cnet_gene_label_map,
                    colorscale=cnet_colorscale, reverse_colorscale=cnet_reverse, layout_seed=cnet_seed,
                )
            style = dew._render_plot_style_controls(f"{key_prefix}_cnet", group_values=None, show_legend_controls=True)
            dew._apply_plot_style(fig, style, default_title=f"{database} {analysis_type.upper()} -- Gene-Concept Network")
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{key_prefix}_cnet", f"{database}_{analysis_type}_gene_concept_network")

    if is_gsea and "running_score" in selected_plots:
        st.markdown("##### GSEA running-score plot")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "This is the classic GSEA enrichment curve for ONE gene "
                "set at a time. The green curve tracks a running "
                "'enrichment score' -- it goes UP whenever it passes "
                "one of this gene set's member genes. The black tick "
                "marks along the bottom show exactly where each member "
                "gene falls."
            )
        term_options = result_df.sort_values("p.adjust")["ID"].tolist()
        term_labels = dict(zip(result_df["ID"], result_df.get("Description", result_df["ID"])))
        picked_term = st.selectbox(
            "Which gene set?", options=term_options,
            format_func=lambda tid: f"{term_labels.get(tid, tid)} ({tid})",
            key=f"{key_prefix}_running_score_term",
        )
        running_df = om.read_gsea_running_score(output_dir, database, picked_term)
        if running_df is None:
            st.info("Running-score data isn't available for this gene set.")
        else:
            fig = _plot_gsea_running_score(running_df, term_labels.get(picked_term, picked_term))
            style = dew._render_plot_style_controls(f"{key_prefix}_running_score", group_values=None, show_legend_controls=True)
            dew._apply_plot_style(fig, style, default_title=f"GSEA Running Score: {term_labels.get(picked_term, picked_term)}")
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{key_prefix}_running_score", f"{database}_gsea_running_score_{picked_term}")

    if is_gsea and "ridge" in selected_plots:
        st.markdown("##### Ridge plot")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Each ridge shows the distribution of fold-change "
                "values among ONE gene set's own member genes, stacked "
                "top to bottom by Normalized Enrichment Score (NES)."
            )
        n_ridge = st.slider(
            "Number of gene sets to show:", min_value=3, max_value=20, value=10, step=1,
            key=f"{key_prefix}_ridge_n_top",
            help="How many of the most significant gene sets get their own ridge. More gene sets = a taller plot with more ridges to compare.",
        )
        ridge_data = om.build_ridge_plot_data(output_dir, database, result_df, n_top=n_ridge)
        if not ridge_data:
            st.info("Not enough data available to build a ridge plot for this result.")
        else:
            fig = _plot_ridge(ridge_data)
            style = dew._render_plot_style_controls(f"{key_prefix}_ridge", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(fig, style, default_title=f"{database} GSEA -- Ridge Plot")
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{key_prefix}_ridge", f"{database}_gsea_ridge_plot")

    if "upset" in selected_plots:
        st.markdown("##### UpSet plot")
        with st.expander("ℹ️ How to read this plot"):
            st.markdown(
                "Each bar shows how many genes belong to EXACTLY that "
                "specific combination of terms -- not just any overlap, "
                "but genes unique to precisely that intersection."
            )
        upset_color = st.color_picker("Bar color:", value="#4C72B0", key=f"{key_prefix}_upset_color")
        n_upset_terms = st.slider(
            "Number of top terms to include:", min_value=3, max_value=15, value=8, step=1,
            key=f"{key_prefix}_upset_n_top",
            help="How many of the most significant terms are compared for overlap. More terms means more possible combinations shown, which can make the plot longer.",
        )
        upset_data = om.build_upset_data(result_df, n_top=n_upset_terms)
        fig = _plot_upset(upset_data, bar_color=upset_color)
        if fig is None:
            st.info("Not enough overlapping terms to build an UpSet plot.")
        else:
            style = dew._render_plot_style_controls(f"{key_prefix}_upset", group_values=None, show_legend_controls=False)
            dew._apply_plot_style(fig, style, default_title=f"{database} {analysis_type.upper()} -- Term Overlap (UpSet)")
            dew._render_plotly_chart(fig)
            dew._render_pdf_export(fig, f"{key_prefix}_upset", f"{database}_{analysis_type}_upset")

    with st.expander(f"📄 Full {database} results table"):
        display_df = result_df.sort_values("p.adjust")
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        dew._render_csv_download(
            result_df, f"{analysis_type}_{database}_results{'_' + direction if direction else ''}", f"{key_prefix}_table",
            expander_label=f"⬇️ Download Full {database} Results (.csv)",
        )


def _render_ora_up_down_comparison(output_dir, database, key_prefix):
    up_df = om.read_enrichment_result(output_dir, "ora", database, direction="up")
    down_df = om.read_enrichment_result(output_dir, "ora", database, direction="down")
    if (up_df is None or up_df.empty) or (down_df is None or down_df.empty):
        return

    with st.expander(f"📊 Compare {database} Up vs Down side-by-side"):
        st.caption(
            "Uses ONE shared set of settings for both sides, so the "
            "two dot plots are directly comparable to each other."
        )
        n_top = st.slider(
            "Number of top terms to show (each side):", min_value=5, max_value=30, value=15, step=5,
            key=f"{key_prefix}_compare_n_top",
        )
        colorscale, reverse = dew._render_heatmap_color_controls(
            f"{key_prefix}_compare", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
        )

        col_up, col_down = st.columns(2)
        with col_up:
            st.markdown("**⬆️ Up-regulated**")
            up_prepared = _prepare_ranked_plot_df(up_df, "p.adjust", True, n_top)
            fig_up = _plot_dotplot(up_prepared, "neg_log10_padj", "-log10(padj)", colorscale=colorscale, reverse_colorscale=reverse)
            fig_up.update_layout(height=max(350, 25 * len(up_prepared)))
            dew._render_plotly_chart(fig_up)
        with col_down:
            st.markdown("**⬇️ Down-regulated**")
            down_prepared = _prepare_ranked_plot_df(down_df, "p.adjust", True, n_top)
            fig_down = _plot_dotplot(down_prepared, "neg_log10_padj", "-log10(padj)", colorscale=colorscale, reverse_colorscale=reverse)
            fig_down.update_layout(height=max(350, 25 * len(down_prepared)))
            dew._render_plotly_chart(fig_down)


# ---------------------------------------------------------------------------
# Main render function
# ---------------------------------------------------------------------------

def render():
    st.title("🧬 Ontology Analysis")
    st.markdown(
        "This workspace takes your Differential Expression results and "
        "asks: **what do these genes have in common, biologically?** "
        "It runs enrichment analysis against three widely-used "
        "databases -- **Gene Ontology (GO)**, **KEGG pathways**, and "
        "**Reactome pathways** -- using two complementary statistical "
        "approaches (explained below). **No statistics background "
        "required** -- each choice and every plot is explained in "
        "plain language as you go."
    )
    st.markdown("---")

    project = pm.render_project_selector(workspace_key=WORKSPACE_KEY)
    if not project:
        st.info("⬆️ Create or select a project above to get started.")
        return

    st.markdown("---")

    if not pm.has_completed_step(project, "deseq2_complete"):
        st.warning(
            "⚠️ This project doesn't have Differential Expression "
            "results yet. Complete that step first so there's a "
            "significant gene list to run enrichment analysis on."
        )
        if st.button("⬅️ Go to Differential Expression", key="ontology_gate_back_btn"):
            st.session_state["nav_request"] = "📈 Differential Expression"
            st.rerun()
        return

    if not om.clusterprofiler_tools_available():
        st.error(
            "⚠️ Rscript was not found on this system. R with "
            "clusterProfiler (and ReactomePA, for Reactome pathways) "
            "needs to be installed in your environment before this "
            "step can run."
        )
        return

    deseq2_out_dir = pm.deseq2_output_dir(project)
    available_contrasts = dm.list_contrast_results(deseq2_out_dir)
    if not available_contrasts:
        st.warning("⚠️ No Differential Expression contrast results were found for this project.")
        return

    species_key, orgdb_package, species_label = _species_and_orgdb_for_project(project)
    if not orgdb_package:
        species_key, orgdb_package, species_label = _render_species_override_picker(project)
        if not orgdb_package:
            return

    st.session_state["_ontology_orgdb_package_for_gosemsim_check"] = orgdb_package

    db_availability = om.databases_available_for_species(species_key)
    kegg_organism = om.get_kegg_organism_code(species_key)
    reactome_organism = om.get_reactome_organism_name(species_key)

    st.success(f"✅ Organism: **{species_label}** (annotation package: `{orgdb_package}`).")
    unavailable = [db for db, ok in db_availability.items() if not ok]
    if unavailable:
        st.info(
            f"ℹ️ {' and '.join(unavailable)} {'is' if len(unavailable) == 1 else 'are'} "
            f"not available for {species_label} and will be grayed out below."
        )

    ontology_out_dir = pm.ontology_output_dir(project)
    ontology_work_dir = pm.ontology_work_dir(project)

    st.markdown("---")

    st.header("Step 1: Choose Your Analysis Approach")

    with st.expander("ℹ️ ORA vs. GSEA vs. compareCluster -- which should I use? (click to learn more)"):
        st.markdown(
            "**🎯 ORA (Over-Representation Analysis)** -- takes a FIXED "
            "list of your significant genes and asks: *\"are these "
            "genes enriched for any GO term/pathway more than you'd "
            "expect by chance?\"*\n\n"
            "**📈 GSEA (Gene Set Enrichment Analysis)** -- uses EVERY "
            "gene that was tested, ranked from most up- to most "
            "down-regulated, and asks whether a term's genes cluster "
            "toward one end of that ranking.\n\n"
            "**🔀 compareCluster** -- runs ORA independently across "
            "**multiple contrasts at once**, for direct side-by-side "
            "comparison."
        )

    analysis_approach = st.radio(
        "Which approach would you like to use?",
        ["ORA (single contrast)", "GSEA (single contrast)", "compareCluster (multiple contrasts)"],
        key="ontology_analysis_approach_radio",
    )

    st.markdown("---")

    st.header("Step 2: Choose Your Contrast(s) and Databases")

    if analysis_approach.startswith("compareCluster"):
        selected_contrasts = st.multiselect(
            "Which contrasts should be compared?",
            options=available_contrasts,
            default=available_contrasts[:min(2, len(available_contrasts))],
            key="ontology_cc_contrasts_select",
        )
        if len(selected_contrasts) < 2:
            st.warning("⚠️ Select at least 2 contrasts to run compareCluster.")
            return
    else:
        selected_contrast = st.selectbox(
            "Which contrast?", options=available_contrasts, key="ontology_single_contrast_select",
        )
        selected_contrasts = [selected_contrast]

    db_col1, db_col2, db_col3 = st.columns(3)
    run_go = db_col1.checkbox("GO (Gene Ontology)", value=True, key="ontology_run_go", disabled=not db_availability["GO"])
    run_kegg = db_col2.checkbox("KEGG pathways", value=db_availability["KEGG"], key="ontology_run_kegg", disabled=not db_availability["KEGG"])
    run_reactome = db_col3.checkbox("Reactome pathways", value=db_availability["Reactome"], key="ontology_run_reactome", disabled=not db_availability["Reactome"])

    go_ontology = "BP"
    simplify_go = False
    simplify_measure = om.DEFAULT_SIMPLIFY_MEASURE
    simplify_cutoff = om.DEFAULT_SIMPLIFY_CUTOFF
    if run_go:
        with st.expander("ℹ️ What's the difference between Biological Process, Molecular Function, and Cellular Component? (click for a detailed explanation)"):
            st.markdown(
                "Gene Ontology organizes its terms into three "
                "completely separate, complementary sub-ontologies:\n\n"
                "**🔬 Biological Process (BP)** -- describes a larger "
                "biological GOAL or PROGRAM a gene contributes to. "
                "Examples: *\"inflammatory response\"*, *\"DNA "
                "replication\"*, *\"apoptotic process\"*. **This is the "
                "most commonly used sub-ontology for RNA-seq "
                "enrichment**, since it most directly answers \"what "
                "biological programs changed in my experiment?\"\n\n"
                "**⚙️ Molecular Function (MF)** -- describes the "
                "specific BIOCHEMICAL ACTIVITY a gene's protein product "
                "performs. Examples: *\"ATP binding\"*, *\"kinase "
                "activity\"*, *\"DNA-binding transcription factor "
                "activity\"*.\n\n"
                "**📍 Cellular Component (CC)** -- describes WHERE in "
                "the cell a gene's protein product is located. "
                "Examples: *\"mitochondrion\"*, *\"nucleus\"*, *\"plasma "
                "membrane\"*.\n\n"
                "**\"All three (BP + MF + CC)\"** combines all three -- "
                "the most complete picture, but GO term "
                "\"simplification\" isn't available when combining all "
                "three, since it requires a similarity measure specific "
                "to one sub-ontology."
            )
        go_ontology_label = st.selectbox(
            "GO sub-ontology:", options=list(om.GO_ONTOLOGY_OPTIONS.keys()), key="ontology_go_subontology_select",
        )
        go_ontology = om.GO_ONTOLOGY_OPTIONS[go_ontology_label]

        if not analysis_approach.startswith("compareCluster"):
            simplify_go, simplify_measure, simplify_cutoff = _render_simplify_control("ontology_step2", go_ontology)

    if not (run_go or run_kegg or run_reactome):
        st.warning("⚠️ Select at least one database to continue.")
        return

    min_gs_size, max_gs_size = _render_gene_set_size_controls("ontology_step2")

    needs_threshold = not analysis_approach.startswith("GSEA")
    padj_threshold, lfc_threshold = 0.05, 1.0
    split_by_direction = True
    if needs_threshold:
        with st.expander("ℹ️ Why do I need to set a threshold here, and what does it mean? (click to learn more)"):
            st.markdown(
                "These thresholds determine exactly which genes from "
                "your Differential Expression results count as "
                "significant for ORA testing."
            )
        col_p, col_l = st.columns(2)
        with col_p:
            padj_threshold = st.slider(
                "Significance threshold (adjusted p-value):", 0.01, 0.20, 0.05, step=0.01, key="ontology_padj_slider",
                help="Genes must have an adjusted p-value BELOW this to count as significant. Lowering it (e.g. 0.01) is stricter -- fewer, more confident genes go into the enrichment test. Raising it (e.g. 0.10) is more lenient -- more genes are included, which can surface more (but potentially noisier) enriched terms.",
            )
        with col_l:
            lfc_threshold = st.slider(
                "Minimum |log2 fold change|:", 0.0, 4.0, 1.0, step=0.1, key="ontology_lfc_slider",
                help="Genes must change by at least this much (in either direction) to count as significant. Raising it keeps only genes with a LARGER effect size, giving a smaller, more strongly-changed gene list; lowering it (down to 0) includes genes with any detectable change, however small.",
            )

        if analysis_approach.startswith("ORA"):
            split_by_direction = _render_direction_split_control("ontology_step2")

    first_export_path = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
    from_type = "SYMBOL"
    if first_export_path is not None and not first_export_path.empty:
        detection = gim.detect_id_type(first_export_path["gene_id"].astype(str).tolist())
        from_type = detection["detected_type"]
        st.caption(
            f"Auto-detected gene ID type: **{from_type}** "
            f"({detection['match_fraction'] * 100:.0f}% matched -- e.g. "
            f"`{'`, `'.join(detection['example_ids'][:3])}`)."
        )

    st.markdown("---")

    st.header("Step 3: Run the Analysis")

    if analysis_approach.startswith("compareCluster"):
        analysis_type = "compareCluster"
        contrast_key_for_state = "_".join(sorted(selected_contrasts))
    elif analysis_approach.startswith("ORA"):
        analysis_type = "ora"
        contrast_key_for_state = selected_contrasts[0]
    else:
        analysis_type = "gsea"
        contrast_key_for_state = selected_contrasts[0]

    output_dir = os.path.join(ontology_out_dir, analysis_type, contrast_key_for_state)
    work_dir = os.path.join(ontology_work_dir, analysis_type, contrast_key_for_state)

    results_already_exist = any(
        om.enrichment_result_exists(output_dir, analysis_type, db)
        or (analysis_type == "ora" and om.available_ora_directions(output_dir, db))
        for db in ["GO", "KEGG", "Reactome"]
    )
    run_label = "🔄 Re-run Analysis" if results_already_exist else "🚀 Run Enrichment Analysis"

    if results_already_exist and not st.session_state.get(f"_ontology_run_clicked_{contrast_key_for_state}"):
        st.success("✅ Results already exist for this exact configuration.")

    if st.button(run_label, key=f"ontology_run_btn_{analysis_type}"):
        st.session_state[f"_ontology_run_clicked_{contrast_key_for_state}"] = True

        with st.spinner(f"Running {analysis_type.upper()}... this may take a few minutes."):
            if analysis_type == "ora":
                export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
                input_path = os.path.join(work_dir, "input_gene_list.csv")
                os.makedirs(work_dir, exist_ok=True)
                export_df.to_csv(input_path, index=False)
                success, log = om.run_ora_analysis(
                    input_path, output_dir, work_dir, orgdb_package, from_type,
                    padj_threshold, lfc_threshold, run_go, go_ontology,
                    run_kegg, kegg_organism, run_reactome, reactome_organism,
                    min_gs_size=min_gs_size, max_gs_size=max_gs_size,
                    simplify_go=simplify_go, simplify_measure=simplify_measure,
                    simplify_cutoff=simplify_cutoff,
                    split_by_direction=split_by_direction,
                )
            elif analysis_type == "gsea":
                export_df = dm.build_clusterprofiler_export(deseq2_out_dir, selected_contrasts[0])
                input_path = os.path.join(work_dir, "input_gene_list.csv")
                os.makedirs(work_dir, exist_ok=True)
                export_df.to_csv(input_path, index=False)
                success, log = om.run_gsea_analysis(
                    input_path, output_dir, work_dir, orgdb_package, from_type,
                    run_go, go_ontology, run_kegg, kegg_organism, run_reactome, reactome_organism,
                    min_gs_size=min_gs_size, max_gs_size=max_gs_size,
                    simplify_go=simplify_go, simplify_measure=simplify_measure,
                    simplify_cutoff=simplify_cutoff,
                )
            else:  # compareCluster
                os.makedirs(work_dir, exist_ok=True)
                contrast_input_paths = {}
                for c in selected_contrasts:
                    export_df = dm.build_clusterprofiler_export(deseq2_out_dir, c)
                    input_path = os.path.join(work_dir, f"input_{c}.csv")
                    export_df.to_csv(input_path, index=False)
                    contrast_input_paths[c] = input_path
                success, log = om.run_compare_cluster_analysis(
                    contrast_input_paths, output_dir, work_dir, orgdb_package, from_type,
                    padj_threshold, lfc_threshold, run_go, go_ontology,
                    run_kegg, kegg_organism, run_reactome, reactome_organism,
                    min_gs_size=min_gs_size, max_gs_size=max_gs_size,
                )

        if not success:
            st.error("Analysis failed. Details below:")
            st.code(log)
            return
        else:
            st.success("✅ Analysis completed successfully.")

            if simplify_go and run_go:
                simplify_outcomes = om.parse_simplify_outcomes_from_log(log)
                if simplify_outcomes:
                    om.save_simplify_status(output_dir, simplify_outcomes)

            with st.expander("View run log"):
                st.code(log)
            results_already_exist = True

    if not results_already_exist:
        return

    st.markdown("---")

    st.header("Step 4: Explore Your Results")

    gene_name_map = {}
    auto_map_path = pm.gene_symbol_map_path(project)
    if os.path.exists(auto_map_path):
        try:
            auto_map_df = pd.read_csv(auto_map_path)
            gene_name_map = dict(zip(auto_map_df["gene_id"].astype(str), auto_map_df["gene_name"].astype(str)))
        except Exception:
            gene_name_map = {}

    if analysis_type == "compareCluster":
        available_dbs = [db for db in ["GO", "KEGG", "Reactome"] if om.enrichment_result_exists(output_dir, "compareCluster", db)]
        if not available_dbs:
            st.info("No compareCluster results were found (or zero terms were returned for any database).")
            return

        for db in available_dbs:
            cc_df = om.read_enrichment_result(output_dir, "compareCluster", db)
            st.markdown(f"#### {db} -- compareCluster results")
            with st.expander("ℹ️ How to read this plot"):
                st.markdown(
                    "One column per contrast you compared, one row per "
                    "enriched term -- dot size shows gene count, dot "
                    "color shows significance."
                )
            cc_colorscale, cc_reverse = dew._render_heatmap_color_controls(
                f"cc_{db}", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
            )
            n_top_cc = st.slider(
                f"Top terms per contrast ({db}):", min_value=3, max_value=30, value=10, step=1,
                key=f"cc_{db}_n_top",
                help="How many of the most significant terms are shown PER contrast column. More terms per contrast = a taller plot showing more of each contrast's results.",
            )
            cc_label_map = _render_term_label_editor(
                cc_df.sort_values("p.adjust").groupby("Cluster").head(n_top_cc).drop_duplicates("ID"),
                f"cc_{db}_labels",
            )
            fig = _plot_compare_cluster_dotplot(
                cc_df, n_top_per_cluster=n_top_cc, label_map=cc_label_map,
                colorscale=cc_colorscale, reverse_colorscale=cc_reverse,
            )
            if fig is None:
                st.info("No data available to plot.")
            else:
                group_values = sorted(cc_df["Cluster"].astype(str).unique())
                cc_cluster_label_map = dew._render_group_label_renaming_controls(
                    f"cc_{db}", group_values, label=f"Rename contrast labels shown for {db} (optional)",
                )
                style = dew._render_plot_style_controls(f"cc_{db}", group_values=None, show_legend_controls=False)
                dew._apply_plot_style(fig, style, default_title=f"{db} -- compareCluster")
                if cc_cluster_label_map:
                    fig.update_xaxes(
                        tickvals=group_values,
                        ticktext=[cc_cluster_label_map.get(v, v) for v in group_values],
                    )
                dew._render_plotly_chart(fig)
                dew._render_pdf_export(fig, f"cc_{db}", f"{db}_compareCluster_dotplot")

            with st.expander(f"📄 Full {db} compareCluster results table"):
                st.dataframe(cc_df.sort_values("p.adjust"), use_container_width=True, hide_index=True)
                dew._render_csv_download(
                    cc_df, f"compareCluster_{db}_results", f"cc_{db}_table",
                    expander_label=f"⬇️ Download Full {db} compareCluster Results (.csv)",
                )
            st.markdown("---")

    elif analysis_type == "ora":
        available_dbs = om.list_available_ora_databases(output_dir)
        if not available_dbs:
            st.info("No results were found for any database (or zero terms were returned).")
            return

        gene_fc_map = _get_gene_fc_map_for_contrast(deseq2_out_dir, selected_contrasts[0], gene_name_map=gene_name_map)

        combined_view_direction_by_db = {}

        for db in available_dbs:
            directions = om.available_ora_directions(output_dir, db)

            if "up" in directions and "down" in directions:
                st.markdown(f"### {db}")
                st.caption(
                    "Analyzed separately for up- and down-regulated "
                    "genes (see Step 2's explanation) -- each has its "
                    "own full, independently customizable results "
                    "section below."
                )
                _render_single_database_results(
                    output_dir, "ora", db, f"ora_{db}_up", padj_threshold,
                    gene_fc_map=gene_fc_map, direction="up",
                    section_label=f"{db} results — Up-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                st.markdown("---")
                _render_single_database_results(
                    output_dir, "ora", db, f"ora_{db}_down", padj_threshold,
                    gene_fc_map=gene_fc_map, direction="down",
                    section_label=f"{db} results — Down-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                st.markdown("---")
                _render_ora_up_down_comparison(output_dir, db, f"ora_{db}")
                combined_view_direction_by_db[db] = "up"
            elif "combined" in directions:
                _render_single_database_results(
                    output_dir, "ora", db, f"ora_{db}", padj_threshold,
                    gene_fc_map=gene_fc_map, direction=None,
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                combined_view_direction_by_db[db] = "combined"
            else:
                only_direction = directions[0]
                _render_single_database_results(
                    output_dir, "ora", db, f"ora_{db}_{only_direction}", padj_threshold,
                    gene_fc_map=gene_fc_map, direction=only_direction,
                    section_label=f"{db} results — {'Up' if only_direction == 'up' else 'Down'}-regulated genes",
                    orgdb_package=orgdb_package, go_ontology=go_ontology,
                )
                combined_view_direction_by_db[db] = only_direction

            st.markdown("---")

        if len(available_dbs) >= 2:
            st.subheader("🔀 Combined View Across Databases")
            with st.expander("ℹ️ Why combine GO, KEGG, and Reactome results?"):
                st.markdown(
                    "GO, KEGG, and Reactome each describe biology from "
                    "a different angle. Seeing your top hits from all "
                    "three together, on one plot -- with each row's "
                    "BACKGROUND shaded by its category (GO's BP/MF/CC, "
                    "or KEGG/Reactome, explained in the plot's own "
                    "legend) -- can reveal a more complete picture, "
                    "while dot color is free to show significance "
                    "directly. You can check/uncheck which categories "
                    "to include, and rename each category's legend "
                    "text, below."
                )
            show_combined = st.checkbox(
                "📊 Show combined dot plot across selected databases", key="ora_combined_checkbox",
            )
            if show_combined:
                combine_dbs = st.multiselect(
                    "Which databases to combine?", options=available_dbs, default=available_dbs,
                    key="ora_combine_db_select",
                )
                n_top_combined = st.slider(
                    "Top terms per database:", min_value=3, max_value=20, value=10, step=1,
                    key="ora_combined_n_top",
                    help="How many of the most significant terms are pulled in FROM EACH database before combining. More terms per database = a longer combined plot.",
                )
                combined_colorscale, combined_reverse = dew._render_heatmap_color_controls(
                    "ora_combined", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
                )
                frames = []
                for db in combine_dbs:
                    direction = combined_view_direction_by_db.get(db)
                    read_direction = None if direction == "combined" else direction
                    df = om.read_enrichment_result(output_dir, "ora", db, direction=read_direction)
                    if df is not None and not df.empty:
                        df = df.sort_values("p.adjust").head(n_top_combined).copy()
                        df["Database"] = db
                        frames.append(df)
                combined_df = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()

                combined_df_categorized = om.derive_category_column(combined_df)
                selected_categories, category_label_map = _render_category_filter_and_label_controls(
                    combined_df_categorized, "ora_combined",
                )
                combined_df_filtered = combined_df_categorized[
                    combined_df_categorized["Category"].isin(selected_categories)
                ]

                combined_label_map = _render_term_label_editor(combined_df_filtered, "ora_combined_labels")
                if combined_df_filtered.empty:
                    fig, category_colors = None, {}
                else:
                    fig, category_colors = _plot_combined_multi_database(
                        combined_df_filtered, label_map=combined_label_map,
                        category_label_map=category_label_map,
                        colorscale=combined_colorscale, reverse_colorscale=combined_reverse,
                    )
                if fig is None:
                    st.info("No data available to combine -- check at least one category above.")
                else:
                    style = dew._render_plot_style_controls(
                        "ora_combined", group_values=None, show_legend_controls=False,
                    )
                    dew._apply_plot_style(fig, style, default_title="Combined Enrichment -- GO / KEGG / Reactome")
                    dew._render_plotly_chart(fig)
                    dew._render_pdf_export(fig, "ora_combined", "combined_multi_database")
                    dew._render_csv_download(
                        combined_df_filtered, "combined_ora_results", "ora_combined_table",
                        expander_label="⬇️ Download Combined Results (.csv)",
                    )

    else:  # gsea
        available_dbs = om.list_available_databases(output_dir, analysis_type)
        if not available_dbs:
            st.info("No results were found for any database (or zero terms were returned).")
            return

        gene_fc_map = _get_gene_fc_map_for_contrast(deseq2_out_dir, selected_contrasts[0], gene_name_map=gene_name_map)

        for db in available_dbs:
            _render_single_database_results(
                output_dir, analysis_type, db, f"{analysis_type}_{db}",
                padj_threshold, gene_fc_map=gene_fc_map,
                orgdb_package=orgdb_package, go_ontology=go_ontology,
            )
            st.markdown("---")

        if len(available_dbs) >= 2:
            st.subheader("🔀 Combined View Across Databases")
            with st.expander("ℹ️ Why combine GO, KEGG, and Reactome results?"):
                st.markdown(
                    "GO, KEGG, and Reactome each describe biology from "
                    "a different angle. Seeing your top hits from all "
                    "three together, on one plot -- with each row's "
                    "BACKGROUND shaded by its category, explained in "
                    "the plot's own legend -- can reveal a more "
                    "complete picture. You can check/uncheck which "
                    "categories to include, and rename each category's "
                    "legend text, below."
                )
            show_combined = st.checkbox(
                "📊 Show combined dot plot across selected databases", key=f"{analysis_type}_combined_checkbox",
            )
            if show_combined:
                combine_dbs = st.multiselect(
                    "Which databases to combine?", options=available_dbs, default=available_dbs,
                    key=f"{analysis_type}_combine_db_select",
                )
                n_top_combined = st.slider(
                    "Top terms per database:", min_value=3, max_value=20, value=10, step=1,
                    key=f"{analysis_type}_combined_n_top",
                    help="How many of the most significant terms are pulled in FROM EACH database before combining. More terms per database = a longer combined plot.",
                )
                combined_colorscale, combined_reverse = dew._render_heatmap_color_controls(
                    f"{analysis_type}_combined", SEQUENTIAL_COLORSCALE_OPTIONS, default="Viridis", default_reverse=True,
                )
                combined_df = om.build_combined_results(output_dir, analysis_type, combine_dbs, n_top_per_db=n_top_combined)

                combined_df_categorized = om.derive_category_column(combined_df)
                selected_categories, category_label_map = _render_category_filter_and_label_controls(
                    combined_df_categorized, f"{analysis_type}_combined",
                )
                combined_df_filtered = combined_df_categorized[
                    combined_df_categorized["Category"].isin(selected_categories)
                ]

                combined_label_map = _render_term_label_editor(combined_df_filtered, f"{analysis_type}_combined_labels")
                if combined_df_filtered.empty:
                    fig, category_colors = None, {}
                else:
                    fig, category_colors = _plot_combined_multi_database(
                        combined_df_filtered, label_map=combined_label_map,
                        category_label_map=category_label_map,
                        colorscale=combined_colorscale, reverse_colorscale=combined_reverse,
                    )
                if fig is None:
                    st.info("No data available to combine -- check at least one category above.")
                else:
                    style = dew._render_plot_style_controls(
                        f"{analysis_type}_combined", group_values=None, show_legend_controls=False,
                    )
                    dew._apply_plot_style(fig, style, default_title="Combined Enrichment -- GO / KEGG / Reactome")
                    dew._render_plotly_chart(fig)
                    dew._render_pdf_export(fig, f"{analysis_type}_combined", "combined_multi_database")
                    dew._render_csv_download(
                        combined_df_filtered, f"combined_{analysis_type}_results", f"{analysis_type}_combined_table",
                        expander_label="⬇️ Download Combined Results (.csv)",
                    )

    st.markdown("---")
    st.success(
        f"🎉 Project `{project}` has ontology enrichment results ready "
        "to explore above."
    )
