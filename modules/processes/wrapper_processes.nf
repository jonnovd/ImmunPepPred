process GENERATE_PEPTIDES {
    // TODO Do I need an environment for this?
    input:
        path fasta
    output:
        path "${fasta.baseName}_allNmers.txt", emit: txt

    script:
    """
        python3 ${params.generatePeptides} -i ${fasta} -p ${fasta.baseName}_allNmers.txt
    """
}