

include {
    RUN_PEPSICKLE as RUN_PEPSICKLE_INPUT
    RUN_PEPSICKLE as RUN_PEPSICKLE_REF
    GENERATE_PEPTIDES
} from './modules/processes/wrapper_processes'

include { 
    TAP_WORKFLOW
} from './modules/workflows/tap_workflow'

workflow {
    
    main: 

        def pepsickle_inFasta_ch = channel.fromPath(params.inputFasta)
                                .splitFasta(by: params.batchSize, file: true)
        
        def refFasta_ch = channel.fromPath(params.refProteome)
                                .splitFasta(by: params.batchSize, file: true)

        RUN_PEPSICKLE_INPUT(pepsickle_inFasta_ch)
        RUN_PEPSICKLE_REF(refFasta_ch)

        def pepsickle_inFasta_out_ch  = RUN_PEPSICKLE_INPUT.out.collectFile(name: 'pepsickleTest-nc-id-udp.txt')
        def pepsickle_refFasta_out_ch = RUN_PEPSICKLE_REF.out.collectFile(name: 'pepsickleTest-nc-id-ref.txt')


        // def tap_in_ch = out_ch.splitText(by: params.batchSize*3, file: true)
        //TAP_WORKFLOW(params.tapInput)

    publish:
        pepsickle_udp_sublist = pepsickle_inFasta_out_ch
        pepsickle_ref_sublist = pepsickle_refFasta_out_ch

        //tap = TAP_WORKFLOW.out.collectFile(name: params.tapOutName)
}

output {
    pepsickle_udp_sublist {
        path 'pepsickle'
        mode 'copy'
    }
    pepsickle_ref_sublist {
        path 'pepsickle'
        mode 'copy'
    }
    // excludedGenes {
    //     path 'excludedGenes'
    //     mode 'copy'
    // }

    // intermediates {
    //     path 'intermediates'
    //     mode 'copy'
    // }
    // tap {
    //     path 'tap'
    //     mode 'copy'
    // }
}