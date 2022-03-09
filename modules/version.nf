process getVersions {

    label 'version'

    publishDir "${params.outdir}/", mode: 'copy'
    publishDir "${params.outdir}/store/", mode: 'copy'

    output:
    path("*.txt")

    script:
    """
    echo Version: ${workflow.manifest.version} > pipeline.info.txt
    echo Primers: ${params.primers} >> pipeline.info.txt
    echo Adapters: ${params.adapters} >> pipeline.info.txt
    echo Shiver initialization directory: ${params.initDir} >> pipeline.info.txt
    echo Pipeline config: ${params.config} >> pipeline.info.txt
    echo Host genome: ${params.hostGenome}/${params.hostGenomeBase} >> pipeline.info.txt
    echo Coverage threshold: ${params.coverage_threshold_eti} >> pipeline.info.txt
    """
}
