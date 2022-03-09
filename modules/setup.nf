process referenceDownload {

    label 'referenceDownload'

    publishDir "${params.hostGenome}/", mode: 'copy'

    output:
    path("${params.hostGenomeBase}*")

    script:
    """
    wget ${params.hostURL}
    bwa index -p ${params.hostGenomeBase} *
    """
}
