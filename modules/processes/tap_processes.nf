

process PREP_DEEPTAP_INPUT {
    conda "${params.env_deeptap}"

    input:
        path peptides_csv
    output: 
        path "${peptides_csv.baseName}_DeepTapIn.csv"
    script:
    """
        { echo "peptide"; cat $peptides_csv; } > "${peptides_csv.baseName}_DeepTapIn.csv"
    """

}

process RUN_DEEPTAP {
    conda "${params.env_deeptap}"

    input:
        path peptides_csv
    output:
        path "${peptides_csv.baseName}_DeepTAP_cla_predresult.csv", emit: deeptap_out
    script:
    """
        python ${params.deeptap} -t cla -f $peptides_csv
    """
}