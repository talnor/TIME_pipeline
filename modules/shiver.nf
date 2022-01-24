process shiver {

    label 'shiver'

    publishDir "${params.outdir}/shiver", mode: 'copy'

    input:
    //assembly.out.contigs.combine(ch_shiverInitDir, ch_shiverConfigFile)
    //  .join(hostDepletion.out.filtered)).view()

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

    # Calculate depth of coverage
    samtools depth -a -d 1000000 {sample}_remap.bam > ${sample}_remap_depth.csv
    """
}
