// TODO
// Header - author etc
include {
    RUN_PEPSICKLE
    GENERATE_PEPTIDES
} from './modules/processes/wrapper_processes'

include { 
    TAP_WORKFLOW
} from './modules/workflows/tap_workflow'

workflow {

    main:
    
    // TODO if statement for if we're cleaving peptides
    RUN_PEPSICKLE(params.udp_out)

    // Process convert protein fasta dbs to peptide dbs
    // TODO if statement for the cleave parameter
    // if params.cleave != None:
    GENERATE_PEPTIDES(params.udp_out, params.refProteome, RUN_PEPSICKLE.out, 
                        params.pepsickleCleavageThreshold, params.nmers, params.noncanonical)

    // TODO
    // Expression Workflow
    // Argument for pre-filtered or not
    // Add check for abundance.tsv file if pre-filtered is indicated
    // Run Kallisto if not pre-filtered

    // TODO
    // Binding Affinity Workflow (TAP) 
    // - explore adding more tools (CLTAP) or develop own
    //TAP_WORKFLOW(GENERATE_PEPTIDES.out.txt)

    // TODO
    // Binding Affinity Workflow (MHC)

    // TODO
    // Binding Stability Workflow (MHC)

    // TODO
    // Immunogenicity Prediction Workflow (pMHC-TCR binding affinity)

    // TODO
    // Peptide Prioritisation Workflow

    publish:
        pepsickle = RUN_PEPSICKLE.out
        allNmersTxt = GENERATE_PEPTIDES.out.txt
        allNmersCsv = GENERATE_PEPTIDES.out.Csv
        //tap = TAP_WORKFLOW.out
}

output {
    pepsickle {
        path 'pepsickle'
        mode 'copy'
    }
    allNmersTxt {
        path 'allNmers'
        mode 'copy'
    }
    allNmersCsv {
        path 'allNmers'
        mode 'copy'
    }
    // tap {
    //     path 'tap'
    //     mode 'copy'
    // }
}

// TODO
// Sub-Workflows
// potential to place these in different files - look at docs and convention