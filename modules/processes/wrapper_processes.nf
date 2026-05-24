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