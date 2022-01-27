// enable dsl2
nextflow.enable.dsl=2

// import modules
include {trimming} from '../modules/preprocessing.nf'
include {hostDepletion} from '../modules/preprocessing.nf'
include {hostStats} from '../modules/preprocessing.nf'
include {assembly} from '../modules/preprocessing.nf'
include {shiver} from '../modules/shiver.nf'
include {infectionEstimation} from '../modules/postprocessing.nf'

workflow timeAnalysis {
    take:
      ch_fastq
      ch_primersFile
      ch_adaptersFile
      ch_shiverConfigFile
      ch_shiverInitDir
      ch_hostGenome

    main:
    //ch_fastq.combine(ch_adaptersFile).view()
    trimming(ch_fastq.combine(ch_adaptersFile))
    hostDepletion(trimming.out.trim.combine(ch_hostGenome))
    hostStats(hostDepletion.out.bam)
    assembly(hostDepletion.out.filtered)
    shiver(assembly.out.contigs.combine(ch_shiverInitDir).combine(ch_shiverConfigFile)
        .join(hostDepletion.out.filtered))
    shiver.out.basefreq.collect().view()
    ch_fastq.map{ it[0] }.toSortedList().view()
    infectionEstimation(shiver.out.basefreq.collect(), ch_fastq.map{ it[0] }.collect())
}