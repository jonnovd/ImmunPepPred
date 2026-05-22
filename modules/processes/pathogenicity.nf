
process RUN_HEX {
    conda "${params.env_pathogenicity}"

    input:
        path peptides
        path iedb_peptides
    output:
        path "pyHex_out.csv", emit: pyhex_out
        //path "deeptap*_DeepTAP_cla_predresult_rank.csv", emit: deeptap_out
    script:
    """
        python ${params.pyHex} --peptides $peptides --iedb $iedb_peptides --output pyHex_out.csv --workers 4
    """
}