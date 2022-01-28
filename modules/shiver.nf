process shiver {

    maxForks 1

    label 'shiver'

    publishDir "${params.outdir}/${sample}/shiver/", mode: 'copy'
    publishDir "${params.outdir}/store/${sample}/", pattern: "*.{csv,txt,fasta,fai}", mode: 'copy'

    input:
    tuple val(sample), path(contigs), path(shiverInitDir), path(shiverConfigFile), path(reads)

    output:
    path("${sample}*"), emit: output
    path("${sample}_remap_BaseFreqs_WithHXB2.csv"), emit: basefreq
    tuple val(sample), path("${sample}_remap_depth.csv"), emit: depth

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

    publishDir "${params.outdir}/initdir/", mode: 'copy'

    input:
    path(primers)
    path(adapters)
    path(config)
    path(references)

    output:
    path("shiverInitDir")

    script:
    """
    shiver_init.sh \
    ./shiverInitDir \
    ${config} \
    ${references} \
    ${adapters} \
    ${primers}

    echo shiver_init.sh \
    ./shiverInitDir \
    ${config} \
    ${references} \
    ${adapters} \
    ${primers} > ./shiverInitDir/setup.log
    cat .command.out .command.err > ./shiverInitDir/setup_build.log
    """
}
