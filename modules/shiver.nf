process shiver {

    maxForks 1

    label 'shiver'

    publishDir "${params.outdir}/${sample}/shiver/", mode: 'copy'

    input:
    tuple val(sample), path(contigs), path(shiverInitDir), path(shiverConfigFile), path(reads)

    output:
    path("${sample}*"), emit: output
    path("${sample}_remap_BaseFreqs_WithHXB2.csv"), emit: basefreq

    script:
    """
    # Align contigs with Shiver
    shiver_align_contigs.sh \
    ${shiverInitDir} \
    ${shiverConfigFile} \
    ${contigs} \
    ${sample}

    # Map reads with Shiver
    shiver_map_reads.sh \
    ${shiverInitDir} \
    ${shiverConfigFile} \
    ${contigs} \
    ${sample} \
    ${sample}.blast \
    ${sample}_cut_wRefs.fasta \
    ${reads[0]} \
    ${reads[1]}

    # Calculate depth of coverage
    samtools depth -a -d 1000000 ${sample}_remap.bam > ${sample}_remap_depth.csv
    """
}

process ShiverInitialisation {

    label 'init'

    publishDir "${params.outdir}/initdir", mode: 'copy'

    input:


    output:
    path(${shiverInitDir})

    script:
    """
    shiver_init.sh \
    ${shiverInitDir} \
    ${shiverConfigFile} \
    ${references} \
    ${adapters} \
    ${primers}

    echo $date > ${shiverInitDir}/setup.log
    echo shiver_init.sh \
    ${shiverInitDir} \
    ${shiverConfigFile} \
    ${references} \
    ${adapters} \
    ${primers} >> ${shiverInitDir}/setup.log
    """
}
