

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

        def pepsickle_inFasta_out_ch  = RUN_PEPSICKLE_INPUT.out.collectFile(name: 'pepsickleInFasta.txt')
        def pepsickle_refFasta_out_ch = RUN_PEPSICKLE_REF.out.collectFile(name: 'pepsickleRefFasta.txt')

        // Process convert protein fasta dbs to peptide dbs
        // TODO if statement for the cleave parameter
        // if params.cleave != None:
        //GENERATE_PEPTIDES(inFasta_ch.peptides, params.refProteome, pepsickle_out_ch, params.pepsickleCleavageThreshold, params.nmers, params.noncanonical)
        // TODO - issue here if the pepsickle_out_ch is empty if we don't want to use pepsickle cleavage??
        // Perhaps change in the python file, that cleavage must be used, but can set threshold to 1, and then a simple if to skip the cleavage filtering in the python file
        // Still think about what this means for the RUN_PEPSICKLE process

        // TODO (BUG) splitting inFasta into diff processes, means duplication of peptides in the results
        // peptideTranscriptMap is only built locally for each process, so many peptides coming from multiple transcripts will be duplicated
        // GENERATE_PEPTIDES(params.inputFasta.combine(pepsickle_inFasta_out_ch).combine(pepsickle_refFasta_out_ch), 
        //                 params.refProteome, params.refGff, params.genesToRemoveFromRef, 
        //                 params.pepsickleCleavageThreshold, params.nmers, params.noncanonical)
        
        def genPep_input_ch = channel.fromPath(params.inputFasta) 

        GENERATE_PEPTIDES(genPep_input_ch.combine(pepsickle_inFasta_out_ch).combine(pepsickle_refFasta_out_ch), 
                        params.refProteome, params.refGff, params.genesToRemoveFromRef, 
                        params.pepsickleCleavageThreshold, params.nmers, params.noncanonical)

        def out_ch = GENERATE_PEPTIDES.out.txt

        def tap_in_ch = out_ch.splitText(by: params.batchSize*3, file: true)
        TAP_WORKFLOW(tap_in_ch)

    publish:
        allNmers = out_ch.collectFile(name: 'udp-large-c-04-allNmers.txt')
        //excludedGenes = GENERATE_PEPTIDES.out.log.collectFile(name: 'excluded_genes-symbols.log')
        //intermediates = GENERATE_PEPTIDES.out.txt
        tap = TAP_WORKFLOW.out.collectFile(name: 'udp-large-c-04-allNmers-TapOut.txt')
}

output {
    allNmers {
        path 'allNmers'
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
    tap {
        path 'tap'
        mode 'copy'
    }
}