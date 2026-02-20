process GENERATE_PEPTIDES {
    // TODO Do I need an environment for this?
    input:
        path udp_out
    output:
        path 'allNmersFile.txt', emit: txt
        path 'allNmersWithTranscriptIDFile.txt', emit: csv

    script:
    """
        python3 ${params.generatePeptides} -i "${udp_out}" -p "allNmersFile.txt" -c 'allNmersWithTranscriptIDFile.txt'
    """
}