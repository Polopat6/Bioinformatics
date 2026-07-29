import streamlit as st
import duckdb
import plotly.express as px
import os

# Set a wide layout for our dashboards side-by-side
st.set_page_config(layout="wide", page_title="Spatial Transcriptomics Atlas")

st.title("🧬 Spatial Transcriptomics Tissue Explorer")
st.markdown("---")

db_path = "spatial_atlas.db"

# Check if the database file exists yet
if not os.path.exists(db_path):
    st.error(f"⚠️ Database file `{db_path}` not found! Please run your Nextflow pipeline or R script first to generate it.")
else:
    # 1. Connect to the shared built-in DuckDB file (Read-Only)
    con = duckdb.connect(database=db_path, read_only=True)

    # 2. Sidebar configuration and data fetching
    st.sidebar.header("📊 Visualization Settings")
    
    # Get a list of unique genes stored in the database for the dropdown
    genes_df = con.execute("SELECT DISTINCT gene FROM gene_expression ORDER BY gene").df()
    gene_list = genes_df['gene'].tolist()
    
    selected_gene = st.sidebar.selectbox("Select a marker gene to visualize:", gene_list)

    # 3. Execute two independent, optimized SQL queries
    # Query A: Get all spatial spots and their clusters (unfiltered by gene)
    cluster_df = con.execute("""
        SELECT spot_id, imagecol, imagerow, seurat_clusters 
        FROM spatial_metadata
    """).df()

    # Query B: Get expression coordinates ONLY for the chosen gene
    expr_query = f"""
        SELECT m.spot_id, m.imagecol, m.imagerow, e.expression
        FROM spatial_metadata m
        JOIN gene_expression e ON m.spot_id = e.spot_id
        WHERE e.gene = '{selected_gene}'
    """
    expr_df = con.execute(expr_query).df()
    con.close()

    # 4. Display layouts using columns
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("📍 Tissue Clusters")
        st.caption("Complete structural anatomical zones calculated by your pipeline.")
        
        fig_clusters = px.scatter(
            cluster_df, x="imagecol", y="imagerow", color="seurat_clusters",
            labels={"imagecol": "X Pixel", "imagerow": "Y Pixel", "seurat_clusters": "Cluster"},
            hover_data=["spot_id"],
            color_discrete_sequence=px.colors.qualitative.Safe
        )
        fig_clusters.update_yaxes(autorange="reversed")
        st.plotly_chart(fig_clusters, use_container_width=True)

    with col2:
        st.subheader(f"✨ Expression: {selected_gene}")
        st.caption(f"Continuous expression heatmap layout across the tissue slide for {selected_gene}.")
        
        fig_expr = px.scatter(
            expr_df, x="imagecol", y="imagerow", color="expression",
            labels={"imagecol": "X Pixel", "imagerow": "Y Pixel", "expression": "Level"},
            hover_data=["spot_id"],
            color_continuous_scale="Viridis"
        )
        fig_expr.update_yaxes(autorange="reversed")
        st.plotly_chart(fig_expr, use_container_width=True)
