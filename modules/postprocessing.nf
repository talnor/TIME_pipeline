process infectionEstimation {

    label 'infectionEstimation'

    publishDir "${params.outdir}/infectionEstimation/", mode: 'copy'
    publishDir "${params.outdir}/store/", pattern: "eti*.csv", mode: 'copy'
    publishDir "${params.outdir}/store/eti_calculations/", pattern: "*pairwise_distance*.csv", mode: 'copy'

    input:
    file(frequencies)
    val(samples)

    output:
    path("*.csv")

    script:
    """
    calculate_eti.py . \
    ${params.coverage_threshold_eti} \
    ${params.ticket} \
    eti_calculations_${params.ticket}_${params.coverage_threshold_eti}X_\$(date +'%Y%m%d-%H%M%S').csv \
    ${samples}
    """
}

process plotCoverage {

    label 'plotCoverage'

    publishDir "${params.outdir}/coverage/", mode: 'copy'
    publishDir "${params.outdir}/store/${sample}/", pattern: "*.png", mode: 'copy'

    input:
    tuple val(sample), path(coverageFile)

    output:
    path("*.png")

    script:
    """
    plot_coverage.py ${coverageFile} ${sample}
    """
}