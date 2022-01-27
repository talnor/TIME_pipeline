process infectionEstimation {

    label 'infectionEstimation'

    publishDir "${params.outdir}/infectionEstimation", mode: 'copy'

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

//#python calculate_eti.py . ${params.min_cov_eti} \
//#> eti_calculations_${params.ticket}_${params.min_cov_eti}X_\$(date +'%Y%m%d-%H%M%S').csv

/*
process plotCoverage {

    label 'plotCoverage'

    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:

    output:
    *.csv

    script:
    """
    python bin/plot_coverage.py $folder coverage.csv ${sample}
    """
}*/