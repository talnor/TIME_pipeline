process trimming {

    label 'trimming'

    publishDir "${params.outdir}/${sample}/trimming/", mode: 'copy'

    input:
    tuple val(sample), path(reads), path(adapters), path(primers)

    output:
    tuple val(sample), path("${sample}_trimmed_{1,2}.fastq.gz"), emit: trim
    tuple val(sample), path("${sample}_trimmed_unpaired_{1,2}.fastq.gz"), emit: unpaired

    script:
    """
    cat ${adapters} ${primers} > adapters_primers.fasta
    trimmomatic PE -threads ${task.cpus} \
    ${reads[0]} ${reads[1]} \
    ${sample}_trimmed_1.fastq.gz ${sample}_trimmed_unpaired_1.fastq.gz \
    ${sample}_trimmed_2.fastq.gz ${sample}_trimmed_unpaired_2.fastq.gz \
    ILLUMINACLIP:adapters_primers.fasta:2:10:7:1:false \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:4:20 \
    MINLEN:50
    """
}

process hostDepletion {

    label 'hostDepletion'

    publishDir "${params.outdir}/${sample}/hostdepletion/", mode: 'copy'

    input:
    tuple val(sample), path(reads), path(database)

    output:
    tuple val(sample), path("${sample}_mapped_sorted.bam"), path("${sample}.bam"), emit: bam
    tuple val(sample), path("${sample}_filtered_{1,2}.fastq.gz"), emit: filtered

    script:
    """
    # Map with bwa mem
    bwa mem -t ${task.cpus} ${database}/${params.hostGenomeBase} \
    ${reads[0]} ${reads[1]} > ${sample}.sam
    samtools view -bh ${sample}.sam > ${sample}.bam

    # Filter out all mapped reads
    samtools view -bh -F 4 ${sample}.bam > ${sample}_mapped.bam
    samtools sort ${sample}_mapped.bam -o ${sample}_mapped_sorted.bam

    # Filter out all unmapped read pairs (both reads in a pair are unmapped)
    samtools view -bh -f 12 ${sample}.bam > ${sample}_unmapped.bam
    samtools sort -n ${sample}_unmapped.bam -o ${sample}_unmapped_sorted.bam
    bedtools bamtofastq -i ${sample}_unmapped_sorted.bam -fq ${sample}_filtered_1.fastq -fq2 ${sample}_filtered_2.fastq
    gzip ${sample}_filtered_*.fastq
    """
}

process hostStats {

    label 'hostStats'

    publishDir "${params.outdir}/${sample}/hostStats/", mode: 'copy'
    publishDir "${params.outdir}/store/${sample}/", pattern: "*.{txt,pdf}", mode: 'copy'

    input:
    tuple val(sample), path("mapped.bam"), path(bam)

    output:
    tuple path("*.txt"), path("*.pdf")

    script:
    """
    # Remove duplicates
    picard MarkDuplicates \
    I=mapped.bam \
    O=${sample}_mapped_sorted_dedup.bam \
    M=${sample}_mapped_sorted_dedup_metrics.txt \
    REMOVE_DUPLICATES=true

    # Remove supplementary alignments
    samtools view -b -F 2048 ${sample}_mapped_sorted_dedup.bam > ${sample}_mapped_sorted_dedup_nosuppl.bam

    # Count mapped reads (excluding supplementary alignments)
    samtools view -c -F 2048 ${bam} > ${sample}_total_reads.txt
    samtools view -c -F 2048 mapped.bam > ${sample}_mapped_reads.txt
    samtools view -c -F 2048 ${sample}_mapped_sorted_dedup.bam > ${sample}_mapped_dedup_reads.txt

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

    errorStrategy { task.exitStatus == 1 ? 'ignore' : 'retry' }

    publishDir "${params.outdir}/${sample}/assembly/", mode: 'copy'
    publishDir "${params.outdir}/store/${sample}/", pattern: "*.fasta", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_contigs.fasta"), emit: contigs

    script:
    """
    spades.py -t ${task.cpus} \
    -1 ${reads[0]} -2 ${reads[1]} \
    -o .
    mv contigs.fasta ${sample}_contigs.fasta
    """
}
