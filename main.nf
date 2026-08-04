include {
    GET_ALL_ALLELES
    RUN_PEPSICKLE as RUN_PEPSICKLE_INPUT
    RUN_PEPSICKLE as RUN_PEPSICKLE_REF
    GENERATE_PEPTIDES
    GET_INPUT_PEPTIDES
    RUN_PATHOGENICITY
    RUN_SELF_SIMILARITY
    RUN_REPITOPE
    PROCESS_HLA_WORKFLOW_OUTPUT
    BUILD_FEATURE_TABLE
} from './modules/processes/wrapper_processes'

include {
    RUN_DEEPIMMUNO
    RUN_PRIME
} from './modules/processes/immunogenicitytools.nf'

include { 
    TAP_WORKFLOW
} from './modules/workflows/tap_workflow'

// Only using the PREPPEPTIDES process that was created in this nf file
// TODO move this process into wrapper_processes
include {
    PREPPEPTIDES
} from './modules/processes/hlapredictiontools'

include { 
    HLA_WORKFLOW
} from './modules/workflows/hla-pred_workflow'

include {
    MTEC_WORKFLOW
} from './modules/workflows/mTEC_workflow'

include {
    IMMUNOGENICITY_WORKFLOW
} from './modules/workflows/immunogenicity_workflow'

workflow {
    
    main: 

        GET_ALL_ALLELES(params.common_hla_alleles, params.input_hla_alleles)
        
        if (params.use_generatePeptides) {
            def pepsickle_input_ch = channel.fromPath(params.inputFasta)
                                            .splitFasta(by: params.batchSize, file: true)
            RUN_PEPSICKLE_INPUT(pepsickle_input_ch)
            // def refFasta_ch = channel.fromPath(params.refProteome)
            //                         .splitFasta(by: params.batchSize, file: true)
            // Design choice to remove any peptides matching human proteome despite cleavage score
            // as we are therefore certain it is not a self peptide 
            // (we could lose a few hits here, but hopefully not important ones)
            // RUN_PEPSICKLE_REF(refFasta_ch)
            //def pepsickle_refFasta_out_ch = //RUN_PEPSICKLE_REF.out.collectFile(name: 'pepsickleRefFasta.txt')

            def genPep_input_ch = channel.fromPath(params.inputFasta)
            // This process is not intended to be parallelised:
            // peptideTranscriptMap is built locally for each process, 
            // so peptides coming from multiple transcripts would be duplicated in the output
            GENERATE_PEPTIDES(genPep_input_ch.combine(RUN_PEPSICKLE_INPUT.out.collectFile()), //.combine(pepsickle_refFasta_out_ch), 
                            params.refProteome, params.refGff, params.genesToRemoveFromRef, 
                            params.pepsickleCleavageThreshold, params.peptide_lengths, params.filter)
        }

        def generatedPeps_ch = params.use_generatePeptides
            ? GENERATE_PEPTIDES.out.txt
            : Channel.empty()

        def customPeps_ch = params.use_customPeptides
            ? Channel.fromPath(params.customPeptides)
            : Channel.empty()

        if (params.use_customPeptides) {
            // Large batch_size ensures only 1 peptide file is output
            customPeps_ch = PREPPEPTIDES(customPeps_ch, params.min_peptide_length, params.max_peptide_length, 900000000).out
        }

        // Merge different input peptide channels
        all_peptides_ch = generatedPeps_ch.mix(customPeps_ch).collect()
        GET_INPUT_PEPTIDES(all_peptides_ch)

        outFiles = GET_INPUT_PEPTIDES.out.collectFile(name: 'in-peptides_all.txt')
        if (params.use_generatePeptides) {
            outFiles = outFiles.mix(GENERATE_PEPTIDES.out.txt.collectFile(name: 'in-peptides_generated-peptides.txt'))
            // Uncomment to store Pepsickle file in output
            //outFiles = outFiles.mix(RUN_PEPSICKLE_INPUT.out.collectFile(name: 'proteasomal-cleavage_pepsickle-out.txt'))
        }

        // HLA-pred
        // Handles parallelisation inside the workflow
        if (params.use_hlapred) {
            hlaPred_output_ch = [ out: Channel.empty() ] as Object
            def hlaPredPeptideFile_ch   = GET_INPUT_PEPTIDES.out
            hlaPred_output_ch           = HLA_WORKFLOW(hlaPredPeptideFile_ch, GET_ALL_ALLELES.out, GET_ALL_ALLELES.out)out.collect()

            // To deal with calculating summary columns and concatenating peptides over split allele channels
            PROCESS_HLA_WORKFLOW_OUTPUT(hlaPred_output_ch, GET_ALL_ALLELES.out, params.common_hla_alleles)

            outFiles = outFiles.mix(PROCESS_HLA_WORKFLOW_OUTPUT.out)
        }

        if (params.use_tap) {
            def tap_input_ch = GET_INPUT_PEPTIDES.out.splitText(by: params.batchSize*3, file: true)
            TAP_WORKFLOW(tap_input_ch)

            //outFiles = outFiles.mix(TAP_WORKFLOW.out.collectFile(name: 'tap_binding_file.csv', storeDir: "${params.outputDir}/tap"))
            // TODO: Collect files below, otherwise we're printing several separate tap files?
            outFiles = outFiles.mix(TAP_WORKFLOW.out.collectFile(name: "tap_deeptap_prediction.csv", keepHeader: true, skip: 1))
        }

        // Process for mTEC expression
        if (params.use_mtec) {
            MTEC_WORKFLOW(GET_INPUT_PEPTIDES.out)

            outFiles = outFiles.mix(MTEC_WORKFLOW.out)//.csv.collectFile(name: 'mtec_expression_file.csv'))
            //outFiles = outFiles.mix(MTEC_WORKFLOW.out.tsv.collectFile(name: 'mtec_expression_file.tsv'))
        }

        if (params.use_pathogenicity) { // TODO: Why did I put params.hla-pred here as well?
            // Uses 107 Gb of RAM when running against 24k pathogen peptides
            RUN_PATHOGENICITY(GET_INPUT_PEPTIDES.out.splitText(by: 400000, file: true), params.iedb_peptides) // params.batchSize*5

            outFiles = outFiles.mix(RUN_PATHOGENICITY.out.collectFile(name: 'pathogenicity_pyHex_out.csv'))
        }

        if (params.use_selfsimilarity) { 
            // 1.3M self peptide sequences
            RUN_SELF_SIMILARITY(GET_INPUT_PEPTIDES.out.splitText(by: 20000, file: true), params.benign_self_peptides) // params.batchSize*5

            outFiles = outFiles.mix(RUN_SELF_SIMILARITY.out.collectFile(name: "selfsimilarity_pyHex_out.csv"))
        }

        if (params.use_deepimmuno || params.use_prime) {
            IMMUNOGENICITY_WORKFLOW(GET_INPUT_PEPTIDES.out.splitText(by: 5000000, file: true), GET_ALL_ALLELES.out)
            outFiles = outFiles.mix(IMMUNOGENICITY_WORKFLOW.out)
        }

        // TODO - Create Workflow to include PREP_REPITOPE
        // if (params.use_repitope) {
        //     RUN_REPITOPE(GET_INPUT_PEPTIDES.out, "${params.peptide_lengths}")
        //     outFiles = outFiles.mix(RUN_REPITOPE.out)
        // }

        if (params.create_feature_table && params.use_deepimmuno && params.use_prime && params.use_pathogenicity && params.use_mtec && params.use_hlapred) {
            immuno_split = IMMUNOGENICITY_WORKFLOW.out.branch {
                prime:      it.name.contains('PRIME')
                deepimmuno: it.name.contains('deepimmuno')
            }
            
            BUILD_FEATURE_TABLE(PROCESS_HLA_WORKFLOW_OUTPUT.out, 
                                MTEC_WORKFLOW.out.mtec_expression_classification, 
                                RUN_PATHOGENICITY.out.collectFile(name: 'pathogenicity_pyHex_out.csv'), 
                                immuno_split.deepimmuno,
                                immuno_split.prime)//,
                                // params.common_hla_alleles,
                                // GET_ALL_ALLELES.out)
            
            outFiles = outFiles.mix(BUILD_FEATURE_TABLE.out)
        }

        // TODO Prediction model
        // if (params.use_ml_model) {
        //     PRIORITISE_IMMUNONGENIC_PEPTIDES(GET_INPUT_PEPTIDES.out.collectFile(name: 'in-peptides_all.txt'),
        //                                     BUILD_FEATURE_TABLE.out,
        //                                     params.prioritisation_model)

        //     outFiles = outFiles.mix(PRIORITISE_IMMUNONGENIC_PEPTIDES.out)
        // }

    publish:
        outFiles = outFiles
        // hex = RUN_HEX.out.collectFile(name: 'hexOut.csv')
        // allNmers = GET_INPUT_PEPTIDES.out.collectFile(name: 'anyXnmers.txt') // Just a way to name the file, not particularly correct
        //excludedGenes = GENERATE_PEPTIDES.out.log.collectFile(name: 'excluded_genes-symbols.log')
        //intermediates = GENERATE_PEPTIDES.out.txt
        // tap = TAP_WORKFLOW.out.collectFile(name: 'test.txt')
}

output {
    outFiles {
        // path { file -> "${file}/"}
        path { file -> 
            def dirName = file.baseName.contains('_') ? file.baseName.split('_')[0] : file.baseName
            "${params.outputDir}/${dirName}/"
        }
        mode 'copy'
        
    } 
}