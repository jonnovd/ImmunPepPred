// TODO
// Header - author etc
include {
    GENERATE_PEPTIDES
} from './modules/processes/wrapper_processes'

include { 
    TAP_WORKFLOW
} from './modules/workflows/tap_workflow'

workflow {

    main:
    // TODO
    // Process convert protein fasta dbs to peptide dbs
    GENERATE_PEPTIDES(params.udp_out)

    // TODO
    // Expression Workflow
    // Argument for pre-filtered or not
    // Add check for abundance.tsv file if pre-filtered is indicated
    // Run Kallisto if not pre-filtered

    // TODO
    // Binding Affinity Workflow (TAP) 
    // - explore adding more tools (CLTAP) or develop own
    TAP_WORKFLOW(GENERATE_PEPTIDES.out.txt)

    // TODO
    // Hla Prediction Workflow

    // TODO
    // Binding Affinity Workflow (MHC)

    // TODO
    // Binding Stability Workflow (MHC)

    // TODO
    // Immunogenicity Prediction Workflow (pMHC-TCR binding affinity)

    // TODO
    // Peptide Prioritisation Workflow

    publish:
        allNmers = GENERATE_PEPTIDES.out.csv
        tap = TAP_WORKFLOW.out
}

output {
    allNmers {
        path 'allNmers'
        mode 'copy'
    }
    tap {
        path 'tap'
        mode 'copy'
    }
}

// TODO
// Sub-Workflows
// potential to place these in different files - look at docs and convention