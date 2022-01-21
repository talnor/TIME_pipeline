process shiver {

    label 'shiver'

    publishDir "${params.outdir}/shiver", mode: 'copy'

    input:
    #1 = the contigs.fasta path
    #2 = a nice name as base for output
    #3 = forward reads
    #4 = reverse reads
    #5 = an output folder
    #6 = the shiver initDir
    #7 = the config file

    output:


    script:
    """
    # Align contigs with Shiver
    shiver_align_contigs.sh \
    ${initDir} \
    ${shiverConfig} \
    ${contigs} \
    ${sample}

    # Map reads with Shiver
    shiver_map_reads.sh \
    ${initDir} \
    ${shiverConfig} \
    ${contigs} \
    ${sample} \
    ${sample}.blast \
    ${sample}_cut_wRefs.fasta \
    ${forward} \
    ${reverse}
    """
}
