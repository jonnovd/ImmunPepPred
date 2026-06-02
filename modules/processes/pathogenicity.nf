
process RUN_HEX {
    conda "${params.env_pathogenicity}"

    input:
        path peptides
        path iedb_peptides
    output:
        path "pathogenicity_pyHex_out.csv", emit: pyhex_out
        //path "deeptap*_DeepTAP_cla_predresult_rank.csv", emit: deeptap_out
    script:
    """
        python ${params.pyHex} --peptides $peptides --reference $iedb_peptides --output pathogenicity_pyHex_out.csv --workers 4
    """
}