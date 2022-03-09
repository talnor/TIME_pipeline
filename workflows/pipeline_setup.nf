// enable dsl2
nextflow.enable.dsl=2

// import modules
include {referenceDownload} from '../modules/setup.nf'

workflow setup {

    main:
        referenceDownload()

}
