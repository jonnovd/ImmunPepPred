process RUN_PEPSICKLE {
    conda "${params.env_pepsickle}"
    
    input:
        path fasta
    output:
        path "cleavedPeps.txt"
    script:
        """
        pepsickle -f $fasta -o "cleavedPeps.txt"
        """
}

// Required input: fasta, cleavageThreshold, nmers, filter
// Optional inputs can be left as NO_FILE to skip those functions in python script
process GENERATE_PEPTIDES {
    // TODO Do I need an environment for this? - probably for the python command
    //conda "${params.env_pepsickle}"
    input:
        tuple path (fasta), path (cleavagePredictions)//, path (refCleavagePredictions)
        path refProteome
        path refGff
        path genesToRemoveFromRef
        val cleavageThreshold
        val nmers
        val filter

    output:
        path "${fasta.baseName}_${filter}_nmers_${nmers}.txt", emit: txt
        path "${fasta.baseName}_${filter}_nmers_${nmers}.csv", emit: csv
        path "${fasta.baseName}_excluded_genes.log", emit: log, optional: true

    // TODO Check argument changes
    script:
    def ref = filter == 'noRefPeps' ? "-r $refProteome" : "" //--reference-cleavage-prediction $refCleavagePredictions" : ""
    def cleave = cleavageThreshold != 'null' ? "-c $cleavagePredictions -t $cleavageThreshold" : ""
    // TODO Decide if this functionality is necessary
    def gff = refGff.name != 'NO_reffGFF' ? "-g $refGff -e $genesToRemoveFromRef" : ""
    def nmers_string = nmers.replace('-', ' ')
    
    """
    python ${params.generatePeptides} \
        -i ${fasta} \
        ${ref} \
        ${cleave} \
        ${gff} \
        -o ${fasta.baseName}_${filter}_nmers_${nmers}.txt \
        -O ${fasta.baseName}_${filter}_nmers_${nmers}.csv \
        -l ${nmers_string}
    """
}

process GET_INPUT_PEPTIDES {
    input:
        path peptideFiles  // receives a list: [file1] or [file1, file2]

    output:
        path 'peptides.txt'

    script:
    """
    for f in ${peptideFiles}; do
        sed -e '\$a\\' "\$f" >> peptides.txt
    done
    """
}

process RUN_PATHOGENICITY {
    conda "${params.env_pathogenicity}"

    input:
        path peptides
        path iedb_peptides
    output:
        path "pathogenicity_pyHex_out.csv", emit: pyhex_out
    script:
    """
        python ${params.pyHexPathogenicity} --peptides $peptides --reference $iedb_peptides --magic-number ${params.pyHex_weight} --output pathogenicity_pyHex_out.csv --workers 4
    """
}

process RUN_SELF_SIMILARITY {
    conda "${params.env_pathogenicity}"

    input:
        path peptides
        path benign_self_peptides
    output:
        path "selfsimilarity_pyHex_out.csv", emit: pyhex_out
    script:
    """
        python ${params.pyHexSelfSimilarity} --peptides $peptides --reference $benign_self_peptides --magic-number ${params.pyHex_weight} --output selfsimilarity_pyHex_out.csv --workers 4
    """
}

process RUN_REPITOPE {
    container "${params.repitope_container}"
    input:
        path peptide_file
        val pepLens

    output:
        path "repitope_out.csv"

    script:
    def pept_len_range = pepLens.tokenize('-').with { it.first() + ':' + it.last() }
    """
        python ${params.repitope} \
            --input ${peptide_file} \
            --home /data/repitope \
            --output "repitope_out.csv" \
            --pept_len_range ${pept_len_range}
    """
}