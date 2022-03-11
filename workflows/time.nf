// enable dsl2
nextflow.enable.dsl=2

// import modules
include {getInfo} from '../modules/version.nf'
include {trimming} from '../modules/preprocessing.nf'
include {hostDepletion} from '../modules/preprocessing.nf'
include {hostStats} from '../modules/preprocessing.nf'
include {assembly} from '../modules/preprocessing.nf'
include {shiver} from '../modules/shiver.nf'
include {infectionEstimation} from '../modules/postprocessing.nf'
include {plotCoverage} from '../modules/postprocessing.nf'

workflow timeAnalysis {
    take:
      ch_fastq
      ch_primersFile
      ch_adaptersFile
      ch_shiverConfigFile
      ch_shiverInitDir
      ch_hostGenome

    main:
        getInfo()
        trimming(ch_fastq.combine(ch_adaptersFile).combine(ch_primersFile))
        hostDepletion(trimming.out.trim.combine(ch_hostGenome))
        hostStats(hostDepletion.out.bam)
        assembly(hostDepletion.out.filtered)
        shiver(assembly.out.contigs.combine(ch_shiverInitDir).combine(ch_shiverConfigFile)
            .join(hostDepletion.out.filtered))
        infectionEstimation(shiver.out.basefreq.collect(), ch_fastq.map{ it[0] }.collect())
        plotCoverage(shiver.out.depth)
}