
/* 
----------------------------------------------------------------------------------------
 James Wright -  May 2025 - V0.0.1 
 HLA Peptide Binding Workflow - MHCFlow
----------------------------------------------------------------------------------------
*/ 

nextflow.enable.dsl=2

// IMPORT MODULES AND PROCESSES
include {
    PREPPEPTIDES
    PREPMHCFLURRY
    MHCFLURRY
    NETMHC
    NETMHCSTABPAN
    MIXMHCPRED
    MHCNUGGETS
    MERGERESULTS
} from '../processes/hlapredictiontools'

// MAIN WORKFLOW
workflow HLA_WORKFLOW {
    
    // TODO: Add MHC class II processing
    // TODO: Add GIBBS CLUSTERING Process
    
    //Workflow Header Print Summary of Settings
    log.info """\
         .                                                  .

         ===================================================
         I C R - MHC-FLOW - HLA BINDING PREDICTION  
         ===================================================

         alleles file       : ${params.hla_alleles}
         min peptide length : ${params.min_peptide_length}
         max peptide length : ${params.max_peptide_length}
         batch size         : ${params.hla_batch_size}

         ===================================================

         .                                                  .
        """
        .stripIndent()

    take:
        ch_rawPeptides

    main:
        //Create channel for the initial peptide list file
        // ch_rawPeptideFile = Channel
        //     .fromPath(params.peptides)

        //Create channel for the HLA allele list file
        ch_hlaAlleleFile = Channel
            .fromPath(params.hla_alleles)

        //Create channel of HLA alleles in ch_hlaAlleleFile
        ch_hlaAlleles = Channel
            .fromPath(params.hla_alleles)
            .splitText()
            .map { it.trim() }
            .toList()
            //.collate(2)
            .map { it -> it.join(',') }
            .view()

        //Process to filter unsuitable peptides and batch peptides into sets of x
        ch_peptideFile = PREPPEPTIDES(ch_rawPeptides, params.min_peptide_length, params.max_peptide_length, params.hla_batch_size)

        // Convert ch_peptideFile channel to a tuple with batch number and file path
        ch_peptide = ch_peptideFile.out
            .flatten()
            .map{ it -> [ ("${it.baseName}" =~ /_(\d+)$/)[0][1], it ] }
            .view()

        // old code placement
        //Process to prepare the input file for MHCFlurry
        //ch_mhcflurry_input = PREPMHCFLURRY(ch_peptide, ch_hlaAlleleFile)

        //Process to run MHCFlurry
        ch_mhcflurry_output = [ out: Channel.empty() ] as Object
        if (params.use_mhcflurry) {
            //Process to prepare the input file for MHCFlurry
            ch_mhcflurry_input = PREPMHCFLURRY(ch_peptide, ch_hlaAlleleFile)
            ch_mhcflurry_output = MHCFLURRY(ch_mhcflurry_input.out)
        }

        //Process to run NetMHC
        ch_mhcnuggets_output = [ out: Channel.empty() ] as Object
        if (params.use_mhcnuggets){
            ch_mhcnuggets_output = MHCNUGGETS(ch_peptide, ch_hlaAlleleFile)
        }

        //Default peptide length 8-11 will need to add param to increase this for longer peptides
        ch_netmhcpan_output = [ out: Channel.empty() ] as Object
        if (params.use_netmhcpan){
            ch_netmhcpan_output = NETMHC(ch_peptide.combine(ch_hlaAlleles))
        }

        //Maximum peptide length is 14!
        ch_mixmhcpred_output = [ out: Channel.empty() ] as Object
        if (params.use_mixmhcpred){
            ch_mixmhcpred_output = MIXMHCPRED(ch_peptide.combine(ch_hlaAlleles))
        }

        ch_combined_results = ch_mhcflurry_output.out
                .mix(ch_netmhcpan_output.out)
                .mix(ch_mixmhcpred_output.out) 
                .mix(ch_mhcnuggets_output.out)
                .groupTuple(by: 0)

        // Merge and reformat results from all tools, find number of weak and strong binding alleles for each peptide
        ch_merged_results = MERGERESULTS(
            ch_combined_results
        )

        ch_netmhcstabpan_output = [ out: Channel.empty() ] as Object
        if (params.use_mhcstabpan){
            ch_netmhcstabpan_output = NETMHCSTABPAN(ch_peptide.combine(ch_hlaAlleles))
        }
        
    emit:
        // Merge batch results into single table
        out = ch_merged_results.out
                                .collectFile(keepHeader: true, skip: 1, name: 'hla_final_prediction_results.tsv')
}