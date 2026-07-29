nextflow.enable.dsl=2

process RUN_SPATIAL_ANALYTICS {
    publishDir "${projectDir}", mode: 'copy'

    input:
    path data_folder  // Tells Nextflow to safely deliver the data directory into the sandbox

    output:
    path "spatial_atlas.db"

    script:
    """
    process_spatial.R
    """
}

workflow {
    // Create a path channel pointing directly to your local downloaded folder
    data_ch = Channel.fromPath("${projectDir}/data", type: 'dir')
    
    RUN_SPATIAL_ANALYTICS(data_ch)
}
