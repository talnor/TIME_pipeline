process trimming {

    label 'trimming'

    publishDir "${params.outdir}/trimming/${sample}", mode: 'copy'

    input:
    set val(sample), file(forward), file(reverse), file(adapters)

    output:
    tuple sample, path("${sample}_trimmed_1.fastq.gz"), path("${sample}_trimmed_2.fastq.gz"), emit: trim
    tuple sample, path("${sample}_trimmed_unpaired_1.fastq.gz"), path("${sample}_trimmed_unpaired_2.fastq.gz"), emit: unpaired
    path "trimmomatic.log.txt", emit: log

    script:
    """
    trimmomatic PE -threads ${params.cpus} -trimlog trimmomatic.log.txt \
    ${forward} ${reverse} \
    ${sample}_trimmed_1.fastq.gz ${sample}_trimmed_unpaired1.fastq.gz \
    ${sample}_trimmed_2.fastq.gz ${sample}_trimmed_unpaired_2.fastq.gz \
    ILLUMINACLIP:${adapters}:2:10:7:1:false \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50
    """
}

/*
process hostDepletion {

    label 'hostDepletion'

    publishDir "${params.outdir}/hostdepletion", mode: 'copy'

    input:
    reads
    genome/transcriptome

    output:
    original bam
    mapped.bam
    filtered fastq.gz

    script:
    """
    # Map with bwa mem
    bwa mem -t ${params.cpus} ${database} \
    ${forward} ${reverse} > ${sample}.sam

    # Filter out all mapped reads
    samtools view -bh ${sample}.sam > ${sample}.bam
    samtools view -bh -F 4 ${sample}.bam > ${sample}_mapped.bam
    samtools sort ${sample}_mapped.bam -o ${sample}_mapped_sorted.bam

    # Filter out all unmapped read pairs
    samtools view -bh -f 12 ${sample}.bam > ${sample}_unmapped.bam #unmapped=both reads in a pair are unmapped
    samtools sort -n ${sample}_unmapped.bam -o ${sample}_unmapped_sorted.bam
    bedtools bamtofastq -i ${sample}_unmapped_sorted.bam -fq ${sample}_filtered_1.fastq -fq2 ${sample}_filtered_1.fastq
    gzip ${sample}_filtered_1.fastq
    gzip ${sample}_filtered_2.fastq
    """
}

process hostStats {

    label 'hostStats'

    publishDir "${params.outdir}/hostStats", mode: 'copy'

    input:
    mapped bam
    original bam

    output:
    *.txt
    *.pdf

    script:
    """
    # Remove duplicates
    picard MarkDuplicates \
    I=${sample}_mapped_sorted.bam \
    O=${sample}_mapped_sorted_dedup.bam \
    M=${sample}_mapped_sorted_dedup_metrics.txt \
    REMOVE_DUPLICATES=true

    # Remove supplementary alignments
    samtools view -b -F 2048 ${sample}_mapped_sorted_dedup.bam > ${sample}_mapped_sorted_dedup_nosuppl.bam

    # Count mapped reads (excluding supplementary alignments)
    samtools view -c -F 2048 ${sample}_mapped.bam > ${sample}_mapped_reads.txt
    samtools view -c -F 2048 ${sample}_mapped_sorted_dedup.bam > ${sample}_mapped_dedup_reads.txt
    samtools view -c -F 2048 ${sample}.sam > ${sample}_total_reads.txt

    # Calculate insert sizes
    picard CollectInsertSizeMetrics \
    I=${sample}_mapped_sorted_dedup_nosuppl.bam \
    O=${sample}_mapped_sorted_dedup_nosuppl_insertsizes.txt \
    H=${sample}_mapped_sorted_dedup_nosuppl_hist.pdf
    grep -A2  "## METRICS CLASS" ${sample}_mapped_sorted_dedup_nosuppl_insertsizes.txt \
    > ${sample}_mapped_sorted_dedup_nosuppl_insertsizes_short.txt
    """
}

process assembly {

    label 'assembly'

    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    filtered fastq.gz

    output:
    contigs.fasta

    script:
    """
    spades.py -t ${params.cpus} \
    -1 ${forward}  -2 ${reverse} \
    -o .
    """
}*/