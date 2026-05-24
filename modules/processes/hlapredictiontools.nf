
//TODO: Add description of the processes
//TODO: set cluster process names for each process

process PREPPEPTIDES {

    /* 
     * This process takes a list of peptides and prepares them for input.
     * It filters the peptides to only include those of length 8-11 and removes any
     * selenocystine (U) or non-standard amino acid containing peptides from the peptide input file.
     * The process also batches the peptides into sets of x.
     *
     * The output is files with the filtered peptides.
     * 
     * Input: 
     * - peptides: A file containing a list of peptides.
     * - min_peptide_length: Minimum length of peptides to include.
     * - max_peptide_length: Maximum length of peptides to include.
     * - batch_size: Number of peptides to include in each batch.
     *
     * Output:
     * - Files containing the filtered peptides in batches of x

    */

    input:
    path peptides
    val min_peptide_length
    val max_peptide_length
    val batch_size

    output:
    path "peptides*.txt", emit: out
    path "*.log", emit: hidden // only for logging into the log directory

    //publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
    """
    python ${params.hlaScripts_path}/process_peptides.py ${peptides} ${min_peptide_length} ${max_peptide_length} ${batch_size} > peptide_processing.log
    """

}

process PREPMHCFLURRY {

    input:
    tuple val(batch_number), path(peptides)
    each path(hla_alleles)

    output:
    tuple val(batch_number), path("mhcflurry_input_${batch_number}.csv"), emit: out

    script:
    """
    python ${params.hlaScripts_path}/generate-mhcflurry_input.py ${peptides} ${hla_alleles} mhcflurry_input_${batch_number}.csv 
    """
}

process MHCFLURRY {

    conda "${params.mhcflurry_env}"

    input:
    tuple val(batch_number), path(input)

    output:
    tuple val(batch_number), path ("mhcflurry_results_${batch_number}.csv"), emit: out
    path "*.log", emit: hidden // only for logging into the log directory

    //publishDir "${params.out_dir}", mode: 'copy', pattern: 'mhcflurry_results*.csv'
    //publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    // Error: mhcflurry-predict not found
    // Old command: mhcflurry-predict ${input} --out mhcflurry_results_${batch_number}.csv > mhcflurry.log
    script:
    """
    ${params.mhcflurry_env}/bin/mhcflurry-predict ${input} --out mhcflurry_results_${batch_number}.csv > mhcflurry.log
    """
}

process NETMHC {
    input:
    tuple val(batch_number), path(peptides), val(hla)

    output:
    tuple val(batch_number), path ("netmhc_results_${batch_number}.txt"), emit: out

    //publishDir "${params.out_dir}", mode: 'copy', pattern: 'netmhc_results*.txt'

    script:
    """
    ${params.netmhc_path} -p ${peptides} -a ${hla} > netmhc_results_${batch_number}.txt
    """
}

// TODO: Can only take input of one peptide length at a time
// TODO: All predictions for non-9mers are approximations
// idk if it's worth using this tool
process NETMHCSTABPAN {
    input:
    tuple val(batch_number), path(peptides), val(hla)

    output:
    tuple val(batch_number), path ("netmhcstabpan_results_${batch_number}.txt"), emit: out

    //publishDir "$projectDir/results/netmhcstabpan", mode: 'copy', pattern: 'netmhcstabpan_results*.txt'

    script:
    def pepLengths_string = params.peptide_lengths.replace('-', ',')
    """
    ${params.netmhcstabpan_path} -p ${peptides} -a ${hla} -l ${pepLengths_string} > netmhcstabpan_results_${batch_number}.txt
    """
}

process MIXMHCPRED {

    conda "${params.mixmhcpred_env}"

    input:
    tuple val(batch_number), path(peptides), val(hla)

    output:
    tuple val(batch_number), path ("mixmhcpred_results_${batch_number}.txt"), emit: out
    path "*.log", emit: hidden // only for logging into the log directory

    //publishDir "${params.out_dir}", mode: 'copy', pattern: 'mixmhcpred_results*.txt'
    //publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
    """
    ${params.hlaScripts_path}/mixmhcpred/MixMHCpred -i ${peptides} -a ${hla} -o mixmhcpred_results_${batch_number}.txt > mixmhcpred_${batch_number}.log
    """
}

process MHCNUGGETS {

    conda "${params.mhcnuggets_env}"

    input:
    tuple val(batch_number), path(peptides)
    each path(hla_alleles)

    output:
    tuple val(batch_number), path ("mhcnuggets_results_${batch_number}.txt"), emit: out
    path "*.log", emit: hidden // only for logging into the log directory

    //publishDir "${params.out_dir}", mode: 'copy', pattern: 'mhcnuggets_results*.txt'
    //publishDir "${params.log_dir}", mode: 'copy', pattern: '*.log'

    script:
    """
    python ${params.hlaScripts_path}/run_mhc-nuggets.py --peptides_file $peptides --hla_file $hla_alleles --output mhcnuggets_results_${batch_number}.txt > mhcnuggets_${batch_number}.log
    """
}

process MERGERESULTS {

    conda "${params.mixmhcpred_env}"

    input:
    tuple val(batch_number), path(files)

    output:
    path "*_merged_results.txt", emit: out

    publishDir "${params.hla_out_dir}", mode: 'copy', pattern: '*_merged_results.txt'

    script:
    """
    python ${params.hlaScripts_path}/merge_hla_prediction_results_improved.py --input $files --output ${batch_number}_merged_results.txt
    """
}