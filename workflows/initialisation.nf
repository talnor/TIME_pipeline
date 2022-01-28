// enable dsl2
nextflow.enable.dsl=2

// import modules
include {ShiverInitialisation} from '../modules/shiver.nf'

workflow initialisation {
    take:
      ch_primersFile
      ch_adaptersFile
      ch_shiverConfigFile
      ch_referenceFile

    main:
        ShiverInitialisation(ch_primersFile, ch_adaptersFile, ch_shiverConfigFile, ch_referenceFile)

}
