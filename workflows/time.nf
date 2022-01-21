// import modules
include {trimming} from '../modules/preprocessing.nf'
//include {hostDepletion} from '../modules/preprocessing.nf'
//include {hostStats} from '../modules/preprocessing.nf'
//include {assembly} from '../modules/preprocessing.nf'
//include {shiver} from '../modules/shiver.nf'
//include {infectionEstimation} from '../modules/postprocessing.nf'

workflow timeAnalysis {
    take:
      ch_fastq
      ch_primersFile
      ch_adaptersFile
      ch_shiverConfigFile
      ch_shiverInitDir

    main:
    trimming(ch_filePairs.combine(ch_adaptersFile))
    //hostDepletion()
    //hostStats()
    //assembly()
    //shiver()
    //infectionEstimation()
}