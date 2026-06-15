process BUILD_AUTOMATON {
    conda "${params.env_mtec}"
    input:
        path peptide_file

    output:
        path "mtec_automaton.pickle", emit: automaton

    """
        python ${params.build_automaton} -p $peptide_file -o "mtec_automaton.pickle"
    """

}

process SEARCH_MTECS {
    conda "${params.env_mtec}"
    input:
        tuple path(automaton), path(peptide_file), path(mtec_file)//, val(batch_number)//tuple path(mtec_file), val(batch_number)

    output:
        path "mtec_peptide_counts_${mtec_file.baseName}.tsv", emit: pepCounts

    """
        python ${params.search_mtecs} -a $automaton -p $peptide_file -s $mtec_file -o "mtec_peptide_counts_${mtec_file.baseName}.tsv"
    """

}

process MERGE_RESULTS {
    conda "${params.env_mtec}"
    input:
        path peptide_counts_files
        val threshold
    output:
        path "mtec_counts_merged.tsv", emit: mtec_raw_counts
        path "mtec_expression_classification_merged.csv", emit: mtec_expression_classification

    """
        python ${params.merge_mtec_files} -i $peptide_counts_files -o "mtec_counts_merged.tsv" -c "mtec_expression_classification_merged.csv" -t $threshold
    """
}