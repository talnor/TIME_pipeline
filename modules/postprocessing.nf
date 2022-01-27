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
    echo ${frequencies}
    echo ${samples}
    touch asd.csv
    calculate_eti.py . ${params.min_cov_eti} \
    > eti_calculations_${params.ticket}_${params.min_cov_eti}X_\$(date +'%Y%m%d-%H%M%S').csv
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