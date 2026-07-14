nextflow.enable.dsl=2

// IMPORT MODULES AND PROCESSES
include {
    PREP_DEEPIMMUNO_INPUT
    RUN_DEEPIMMUNO
    PROCESS_DEEPIMMUNO_OUTPUT
    RUN_PRIME
} from '../processes/immunogenicitytools.nf'

// MAIN WORKFLOW
workflow IMMUNOGENICITY_WORKFLOW {

    take:
        peptide_file_ch

    main:

        outFiles = Channel.empty()

        // ONLY RUNS ON 9 and 10mers
        if (params.use_deepimmuno) {

            PREP_DEEPIMMUNO_INPUT(peptide_file_ch.splitText(by: 5000000, file: true), params.hla_alleles)

            RUN_DEEPIMMUNO(PREP_DEEPIMMUNO_INPUT.out, file("${params.deepimmuno_dir}/data"), file("${params.deepimmuno_dir}/models"))
            
            PROCESS_DEEPIMMUNO_OUTPUT(RUN_DEEPIMMUNO.out.collectFile())
            
            outFiles = outFiles.mix(PROCESS_DEEPIMMUNO_OUTPUT.out)
        }

        if (params.use_prime) {
            // prime_allele_ch = Channel.fromPath(params.hla_alleles)
            //                         .splitText()
            //                         .collect { allele ->
            //                             "$allele, "
            //                         }
            RUN_PRIME(peptide_file_ch.splitText(by: 5000000, file: true), params.hla_alleles)
            outFiles = outFiles.mix(RUN_PRIME.out.collectFile(name: "immunogenicity_PRIME_results.txt", skip: 12, keepHeader: true))
        }
    emit:
        outFiles

}