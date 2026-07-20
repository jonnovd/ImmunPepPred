
// process PREP_DEEPIMMUNO_INPUT {
//     conda "${params.env_deepimmuno}"

//     input:
//         path peptide_file
//         path allele_file
//         val min_peptide_length
//         val max_peptide_length

//     output:
//         path "deepimmuno_in_9mers.csv", optional: true
//         path "deepimmuno_in_10mers.csv", optional: true
//         path "deepimmuno_in_11mers.csv", optional: true

//     script:
//     """
//         python ${params.prep_deepimmuno} -p $peptide_file -a $allele_file --output deepimmuno_in_9_10mers.csv --output11 deepimmuno_in_11mers.csv
//     """
// }

process PREP_DEEPIMMUNO_INPUT {
    conda "${params.env_deepimmuno}"

    input:
        path peptide_file
        path allele_file
        val min_peptide_length
        val max_peptide_length

    output:
        path "deepimmuno_in_9mers.csv",  optional: true, emit: mer9
        path "deepimmuno_in_10mers.csv", optional: true, emit: mer10
        path "deepimmuno_in_11mers.csv", optional: true, emit: mer11

    script:
        def minLen = min_peptide_length.toInteger()
        def maxLen = max_peptide_length.toInteger()

        def output9  = (minLen <= 9  && maxLen >= 9)  ? "--output9 deepimmuno_in_9mers.csv"   : ""
        def output10 = (minLen <= 10 && maxLen >= 10) ? "--output10 deepimmuno_in_10mers.csv" : ""
        def output11 = (minLen <= 11 && maxLen >= 11) ? "--output11 deepimmuno_in_11mers.csv" : ""
    """
    python ${params.prep_deepimmuno} -p $peptide_file -a $allele_file --minlength $min_peptide_length \
    --maxlength $max_peptide_length $output9 $output10 $output11
    """
}

// process RUN_DEEPIMMUNO {
//     conda "${params.env_deepimmuno}"

//     input:
//         path input_peptide_allele_csv
//         path deepimmuno_data    // Just to stage the data file in the work dir
//         path deepimmuno_models  // same
//     output:
//         path "immunogenicity_deepimmuno-out.txt", emit: deepimmuno_out
//     script:
//     """
//         python ${params.deepimmuno_dir}/deepimmuno-cnn.py --mode "multiple" --intdir $input_peptide_allele_csv --outdir .
//         # mv "deepimmuno-cnn-result.txt" immunogenicity_deepimmuno-out.txt
//     """
// }

process RUN_DEEPIMMUNO {
    conda "${params.env_deepimmuno}"

    input:
        path input_peptide_allele_csv
        path deepimmuno_data    // Just to stage the data file in the work dir
        path deepimmuno_models  // same

    output:
        path "immunogenicity_deepimmuno_*.csv", emit: deepimmuno_out

    script:
        def lengroupMatch = (input_peptide_allele_csv.baseName =~ /(\d+)mers/)
        if (!lengroupMatch) {
            error "Could not determine peptide length group from input filename: ${input_peptide_allele_csv.name}"
        }
        def lengroup = lengroupMatch[0][1]
    """
    python ${params.deepimmuno_dir}/deepimmuno-cnn.py --mode "multiple" --intdir $input_peptide_allele_csv --outdir .
    python ${params.process_deepimmuno_out} -f deepimmuno-cnn-result.txt --lengroup $lengroup --output immunogenicity_deepimmuno_${lengroup}mers.csv
    """
}

process PROCESS_DEEPIMMUNO_OUTPUT {
    conda "${params.env_deepimmuno}"

    input:
        path deepimmuno_csv

    output:
        path "immunogenicity_deepimmuno.csv"

    script:
    """
        python ${params.process_deepimmuno_out} -f $deepimmuno_csv --output immunogenicity_deepimmuno.csv
    """
}

process RUN_PRIME {
    conda "${params.mixmhcpred_env}"

    input:
        path peptideFile
        path allelesFile
    output:
        path "immunogenicity_prime_out.txt"
    script:
        """
        alleles=\$(cut -d'-' -f2 ${allelesFile} | tr -d ':' | tr '\n' ',' | sed 's/,\$//')
        bash ${params.prime} -i $peptideFile -o "immunogenicity_prime_out.txt" -a \$alleles -mix "${params.mixmhcpred_binary}"
        """
}