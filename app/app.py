import streamlit as st
import duckdb
import plotly.express as px
import os
import json
import pandas as pd
from PIL import Image

st.set_page_config(layout="wide", page_title="Spatial Transcriptomics Atlas")
st.title("🧬 Spatial Transcriptomics Tissue Explorer")
st.markdown("---")

db_path = "spatial_atlas.db"
image_path = "detected_tissue_image.jpg"  # High-Res Aligned 10x Space Ranger asset
json_path = "scalefactors_json.json"

if not os.path.exists(db_path):
    st.error(f"⚠️ Database file `{db_path}` not found! Please execute your Nextflow pipeline first.")
else:
    con = duckdb.connect(database=db_path, read_only=True)
    
    st.sidebar.header("📊 Visualization Settings")
    view_mode = st.sidebar.selectbox(
        "Select Field View Mode:",
        ["Expression and Tissue", "Expression Only", "Tissue Only"]
    )
    
    genes_df = con.execute("SELECT DISTINCT gene FROM gene_expression ORDER BY gene").df()
    gene_list = genes_df['gene'].tolist()
    selected_gene = st.sidebar.selectbox("Select a marker gene to visualize:", gene_list)

    cluster_df = con.execute("SELECT spot_id, imagecol, imagerow, seurat_clusters FROM spatial_metadata").df()

    expr_query = f"""
        SELECT m.spot_id, m.imagecol, m.imagerow, e.expression
        FROM spatial_metadata m
        JOIN gene_expression e ON m.spot_id = e.spot_id
        WHERE e.gene = '{selected_gene}'
    """
    expr_df = con.execute(expr_query).df()
    con.close()

    # 📊 SUMMARY METRICS ROWS
    total_spots = len(cluster_df)
    total_clusters = cluster_df['seurat_clusters'].nunique()
    avg_expression = round(expr_df['expression'].mean(), 2)
    max_expression = round(expr_df['expression'].max(), 2)

    m_col1, m_col2, m_col3, m_col4 = st.columns(4)
    with m_col1:
        st.metric(label="Total Tissue Spots Indexed", value=f"{total_spots:,}")
    with m_col2:
        st.metric(label="Identified Tissue Clusters", value=total_clusters)
    with m_col3:
        st.metric(label=f"Average {selected_gene} Intensity", value=avg_expression)
    with m_col4:
        st.metric(label=f"Peak {selected_gene} Value", value=max_expression)
        
    st.markdown("---")

    # 🖼️ 📊 OFFICIAL 10X ALIGNED TISSUE ALIGNMENT MATRIX ENGINE
    show_background = view_mode in ["Expression and Tissue", "Tissue Only"] and os.path.exists(image_path)
    
    if os.path.exists(image_path):
        raw_img = Image.open(image_path)
        img_width, img_height = raw_img.size
    else:
        img_width, img_height = 2000, 2000
        raw_img = None

    if show_background and os.path.exists(json_path):
        with open(json_path, 'r') as f:
            scale_factors = json.load(f)
        
        # Extract the high-res scaling factor
        scaling_multiplier = scale_factors.get('tissue_hires_scalef', 1.0)
        
        # Translate coordinates straight onto the high-resolution pixel footprint
        cluster_df['aligned_x'] = cluster_df['imagecol'] * scaling_multiplier
        cluster_df['aligned_y'] = cluster_df['imagerow'] * scaling_multiplier
        expr_df['aligned_x'] = expr_df['imagecol'] * scaling_multiplier
        expr_df['aligned_y'] = expr_df['imagerow'] * scaling_multiplier
        
        plot_x_col, plot_y_col = 'aligned_x', 'aligned_y'
        
        # Force graph coordinate constraints to lock strictly onto image borders
        x_min_bound, x_max_bound = 0, img_width
        y_min_bound, y_max_bound = 0, img_height
    else:
        plot_x_col, plot_y_col = 'imagecol', 'imagerow'
        cluster_df['aligned_x'] = cluster_df['imagecol']
        cluster_df['aligned_y'] = cluster_df['imagerow']
        expr_df['aligned_x'] = expr_df['imagecol']
        expr_df['aligned_y'] = expr_df['imagerow']
        x_min_bound, x_max_bound = int(cluster_df['imagecol'].min() * 0.95), int(cluster_df['imagecol'].max() * 1.05)
        y_min_bound, y_max_bound = int(cluster_df['imagerow'].min() * 0.95), int(cluster_df['imagerow'].max() * 1.05)

    # 📈 DUAL GRAPHS COLUMNS
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📍 Tissue Clusters Overlay")
        
        # In Tissue Only mode, we hide the cluster scatter dots by forcing opacity to 0
        fig_clusters = px.scatter(
            cluster_df, x=plot_x_col, y=plot_y_col, color="seurat_clusters",
            labels={plot_x_col: "X Axis (Image Pixels)", plot_y_col: "Y Axis (Image Pixels)", "seurat_clusters": "Cluster"},
            hover_data=["spot_id"],
            color_discrete_sequence=px.colors.qualitative.Safe,
            opacity=0.0 if view_mode == "Tissue Only" else 0.65
        )
        
        fig_clusters.update_traces(marker=dict(size=8 if show_background else 6))
        
        if show_background:
            fig_clusters.update_layout(
                images=[dict(
                    source=raw_img, xref="x", yref="y", x=0, y=0,
                    sizex=img_width, sizey=img_height,
                    sizing="stretch", xanchor="left", yanchor="top",
                    opacity=0.8, layer="below"
                )],
                xaxis=dict(range=[x_min_bound, x_max_bound], showgrid=False, zeroline=False, visible=False),
                yaxis=dict(range=[y_max_bound, y_min_bound], showgrid=False, zeroline=False, visible=False),
                margin=dict(l=5, r=5, t=5, b=5),
                legend=dict(orientation="h", yanchor="top", y=-0.02, xanchor="center", x=0.5)
            )
            # Hide the legend box entirely if the dots are invisible
            if view_mode == "Tissue Only":
                fig_clusters.update_layout(showlegend=False)
        else:
            fig_clusters.update_yaxes(autorange="reversed")
        st.plotly_chart(fig_clusters, use_container_width=True)

    with col2:
        st.subheader(f"✨ Expression Heatmap Overlay: {selected_gene}")
        
        if view_mode == "Tissue Only":
            st.info("💡 Field View Mode set to 'Tissue Only'. Gene expression scatter points are hidden.")
            # 💡 FIX SECURED: Added proper bracket arrays [0] to create a valid single-row dataframe safely
            dummy_df = pd.DataFrame({plot_x_col: [0], plot_y_col: [0], 'expression': [0]})
            fig_expr = px.scatter(dummy_df, x=plot_x_col, y=plot_y_col, opacity=0.0)
            fig_expr.update_coloraxes(showscale=False)  # Hide the viridis sidebar gradient
        else:
            fig_expr = px.scatter(
                expr_df, x=plot_x_col, y=plot_y_col, color="expression",
                labels={plot_x_col: "X Axis (Image Pixels)", plot_y_col: "Y Axis (Image Pixels)", "expression": "Level"},
                hover_data=["spot_id"],
                color_continuous_scale="Viridis",
                opacity=0.75 if show_background else 1.0
            )
            fig_expr.update_traces(marker=dict(size=8 if show_background else 6))
            
        if show_background:
            fig_expr.update_layout(
                images=[dict(
                    source=raw_img, xref="x", yref="y", x=0, y=0,
                    sizex=img_width, sizey=img_height,
                    sizing="stretch", xanchor="left", yanchor="top",
                    opacity=0.8, layer="below"
                )],
                xaxis=dict(range=[x_min_bound, x_max_bound], showgrid=False, zeroline=False, visible=False),
                yaxis=dict(range=[y_max_bound, y_min_bound], showgrid=False, zeroline=False, visible=False),
                margin=dict(l=5, r=5, t=5, b=5)
            )
        else:
            fig_expr.update_yaxes(autorange="reversed")
        st.plotly_chart(fig_expr, use_container_width=True)
