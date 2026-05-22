process BUILD_AUTOMATON {
    conda "${params.env_mtec}"
    input:
        path peptide_file

    output:
        path "automaton.pickle"

    """
        python ${params.build_automaton} -p $peptide_file -o "automaton.pickle"
    """

}

process SEARCH_MTECS {
    conda "${params.env_mtec}"
    input:
        tuple path(automaton), path(peptide_file), path(mtec_file), val(batch_number)//tuple path(mtec_file), val(batch_number)

    output:
        path "peptide_counts_${batch_number}.tsv"

    """
        python ${params.search_mtecs} -a $automaton -p $peptide_file -s $mtec_file -o "peptide_counts_${batch_number}.tsv"
    """

}

process MERGE_RESULTS {
    conda "${params.env_mtec}"
    input:
        path peptide_counts_files
        val threshold
    output:
        path "mtec_counts_merged.tsv", emit: tsv
        path "peptide_expression_in_mtecs_merged.tsv", emit: csv

    """
        python ${params.merge_mtec_files} -i $peptide_counts_files -o "mtec_counts_merged.tsv" -c "peptide_expression_in_mtecs_merged.tsv" -t $threshold
    """
}