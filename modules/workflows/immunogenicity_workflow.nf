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
        all_hla_alleles_ch

    main:

        outFiles = Channel.empty()

        // ONLY RUNS ON 9 and 10mers
        if (params.use_deepimmuno) {

            PREP_DEEPIMMUNO_INPUT(
                peptide_file_ch.splitText(by: 5000000, file: true),
                all_hla_alleles_ch,
                params.min_peptide_length,
                params.max_peptide_length
            )

            deepimmuno_inputs = PREP_DEEPIMMUNO_INPUT.out.mer9
                .mix(PREP_DEEPIMMUNO_INPUT.out.mer10, PREP_DEEPIMMUNO_INPUT.out.mer11)

            RUN_DEEPIMMUNO(
                deepimmuno_inputs,
                file("${params.deepimmuno_dir}/data"),
                file("${params.deepimmuno_dir}/models")
            )
            
            outFiles = outFiles.mix(RUN_DEEPIMMUNO.out.deepimmuno_out
                                                .collectFile(name: "immunogenicity_deepimmuno_all.csv", keepHeader: true, skip: 1))
        }

        if (params.use_prime) {
            RUN_PRIME(peptide_file_ch.splitText(by: 5000000, file: true), all_hla_alleles_ch)
            outFiles = outFiles.mix(RUN_PRIME.out.collectFile(name: "immunogenicity_PRIME_results.txt", skip: 12, keepHeader: true))
        }
    emit:
        outFiles

}