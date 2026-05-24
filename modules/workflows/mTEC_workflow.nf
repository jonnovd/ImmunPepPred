// Implements Aho-Corasick python implementation for efficient searching of many peptides against a db
// For 20 Mil peptides against mTEC translated read files with 97 Mil x 26 AA peps: Runtime ~ 15 mins
// The process shouldn't need more than 8Gb RAM but worth checking

nextflow.enable.dsl=2

// IMPORT MODULES AND PROCESSES
include {
    BUILD_AUTOMATON
    SEARCH_MTECS
    MERGE_RESULTS
} from '../processes/mTEC_processes.nf'

// MAIN WORKFLOW
workflow MTEC_WORKFLOW {
    take:
        peptide_file_ch

    main:
        automaton_ch = BUILD_AUTOMATON(peptide_file_ch)

        mtec_files_ch   = Channel.fromPath(params.mtec_files)
        //                            .splitText(by: params.mtec_batchSize)
        // Convert to Channel containing: file, batch_number
        // BUG: If you decide to split the mtec read files, we may get a bug here,
        // depending on how you name the split files
        mtec_files_ch = mtec_files_ch.map { file -> [file, ("${file.baseName}" =~ /_(\d+)$/)[0][1]]}
                                    //.view() // To print to terminal for debug

        search_ch       = automaton_ch.combine(peptide_file_ch).combine(mtec_files_ch)
        counts_file_ch  = SEARCH_MTECS(search_ch)

        merged_results_ch = MERGE_RESULTS(counts_file_ch.collect(), params.mtec_expression_threshold)

    emit:
        mtec_raw_counts = merged_results_ch.mtec_raw_counts
        mtec_expression_classification = merged_results_ch.mtec_expression_classification

}