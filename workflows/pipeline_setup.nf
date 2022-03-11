// enable dsl2
nextflow.enable.dsl=2

// import modules
include {buildDatabase} from '../modules/setup.nf'

workflow setup {
    take:
      ch_hostFasta

    main:
        buildDatabase(ch_hostFasta)

}
