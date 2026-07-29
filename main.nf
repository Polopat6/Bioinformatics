nextflow.enable.dsl=2

process RUN_SPATIAL_ANALYTICS {
    // Copy the final database out of the working folder into your main directory
    publishDir "${projectDir}", mode: 'copy'

    output:
    path "spatial_atlas.db"

    script:
    """
    process_spatial.R
    """
}

workflow {
    RUN_SPATIAL_ANALYTICS()
}
