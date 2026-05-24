
process PREP_DEEPTAP_INPUT {
    input:
        path peptides_csv
    output: 
        path "${peptides_csv.baseName}_DeepTapIn.csv", emit: deeptap_in
    script:
    """
        { echo "peptide"; cat $peptides_csv; } > "${peptides_csv.baseName}_DeepTapIn.csv"
    """

}

process RUN_DEEPTAP {
    conda "${params.env_deeptap}"

    input:
        path peptides_csv
        //path 'deeptap*.csv'
    output:
        path "tap_DeepTAP_cla_predresult_rank.csv", emit: deeptap_out
        //path "deeptap*_DeepTAP_cla_predresult_rank.csv", emit: deeptap_out
    script:
    """
        # python ${params.deeptap} -t cla -f '$peptides_csv' -o "."
        python ${params.rundeeptap} --deeptap ${params.deeptap} -t cla -f '$peptides_csv' -o "."
        mv "${peptides_csv.baseName}_DeepTAP_cla_predresult_rank.csv" "tap_DeepTAP_cla_predresult_rank.csv"
    """
}