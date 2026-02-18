// TODO
// Header - author etc

include { 
    TAP_WORKFLOW
} from './modules/workflows/tap_workflow'

// TODO
// Main Workflow

workflow {

    main:
    // TODO
    // Process convert protein fasta dbs to peptide dbs

    // TODO
    // Expression Workflow
    // Argument for pre-filtered or not
    // Add check for abundance.tsv file if pre-filtered is indicated
    // Run Kallisto if not pre-filtered

    // TODO
    // Hla Prediction Workflow

    // TODO
    // Binding Affinity Workflow (TAP) - explore adding more tools or develop own
    TAP_WORKFLOW()

    // TODO
    // Binding Affinity Workflow (MHC)

    // TODO
    // Binding Stability Workflow (MHC)

    // TODO
    // Immunogenicity Prediction Workflow (pMHC-TCR binding affinity)

    // TODO
    // Peptide Prioritisation Workflow

    publish:
    tap = TAP_WORKFLOW.out
}

output {
    tap {
        path 'tap'
        mode 'copy'
    }
}

// TODO
// Sub-Workflows
// potential to place these in different files - look at docs and convention