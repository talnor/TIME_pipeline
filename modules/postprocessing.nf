process infectionEstimation {

    label 'infectionEstimation'

    publishDir "${params.outdir}/infectionEstimation", mode: 'copy'

    input:

    output:
    *.csv

    script:
    """
    python calculate_eti.py ${resultsDir} ${min_cov_eti} \
    > eti_calculations_${ticket}_${min_cov_eti}X_${DATE}.csv
    """
}

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
}