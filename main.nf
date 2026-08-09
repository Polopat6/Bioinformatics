nextflow.enable.dsl=2

process RUN_SPATIAL_ANALYTICS {
    publishDir "${projectDir}", mode: 'copy'

    input:
    path 'data'

    output:
    path "spatial_atlas.db"
    path "detected_tissue_image.jpg"  // 💡 Add this output rule
    path "scalefactors_json.json"

    script:
    """
    process_spatial.R
    """
}

workflow {
    data_ch = Channel.fromPath("${projectDir}/data", type: 'dir')
    RUN_SPATIAL_ANALYTICS(data_ch)
}
