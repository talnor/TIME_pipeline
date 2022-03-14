process buildDatabase {

    label 'buildDatabase'

    publishDir "${params.outdir}/database/", mode: 'copy'
    publishDir "${params.hostGenome}", mode: 'copy'

    input:
    path(fasta)

    output:
    path("${params.hostGenomeBase}*")

    script:
    """
    bwa index -p ${params.hostGenomeBase} ${fasta}
    """
}
