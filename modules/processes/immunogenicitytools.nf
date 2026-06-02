
process RUN_DEEPIMMUNO {
    conda "${params.env_deepimmuno}"

    input:
        path input_peptide_allele_csv
        path deepimmuno_data    // Just to stage the data file in the work dir
        path deepimmuno_models  // same
    output:
        path "immuno_deepimmuno-out.txt", emit: deepimmuno_out
    script:
    """
        python ${params.deepimmuno_dir}/deepimmuno-cnn.py --mode "multiple" --intdir $input_peptide_allele_csv --outdir .
        mv "deepimmuno-cnn-result.txt" immuno_deepimmuno-out.txt
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