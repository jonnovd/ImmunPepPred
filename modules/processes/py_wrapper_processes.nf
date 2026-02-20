process GENERATE_PEPTIDES {
    // TODO Do I need an environment for this?
    input:
        path udp_out
    output:
        path 'allNmersFile.txt'

    script:
    """
        python3 ${params.generatePeptides} -f "${udp_out}" -o "allNmersFile.txt"
    """
}