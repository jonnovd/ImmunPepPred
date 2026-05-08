import re
import pandas as pd
import argparse

def load_mhcflurry_results(file_path, merged_df):
    """
    Load MHCFlurry prediction results from a file and merge them into the main DataFrame.
    Assumes the file contains columns: 
    peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_processing_score,mhcflurry_presentation_score,mhcflurry_presentation_percentile.
    """
    print(f"Loading MHCFlurry results from {file_path}")
    try:
        df = pd.read_csv(file_path)
        df = df[['peptide', 'allele', 'mhcflurry_presentation_percentile']]
        df.rename(columns={'mhcflurry_presentation_percentile': 'MHCFlurry_%Rank'}, inplace=True)
        merged_df = pd.merge(merged_df, df, on=['peptide', 'allele'], how='outer')

    except Exception as e:
        print(f"Error loading MHCFlurry results from {file_path}: {e}")
    return merged_df

def load_netmhcpan_results(file_path, merged_df):
    """
    Load NetMHCpan prediction results from a file and merge them into the main DataFrame.
    Assumes the file contains columns: peptide,allele,netmhcpan_affinity,netmhcpan_affinity_percentile.
    """
    print(f"Loading NetMHCpan results from {file_path}")
    try:
        rows = []
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("   1 HLA-"):
                    cols = re.split(r'\s+', line.strip())

                    rows.append({
                        'peptide': cols[2],
                        'allele': cols[1].replace("*", ""),
                        'NetMHCpan_%Rank': cols[12]
                    })
        if rows:
            df = pd.DataFrame(rows)
            return merged_df.merge(df, on=['peptide', 'allele'], how='outer')
            # return pd.concat([merged_df, new_df], ignore_index=True)

    except Exception as e:
        print(f"Error loading NETMHC results from {file_path}: {e}")
    return merged_df

def reformat_allele(allele):
    """
    Reformat an allele from 'A0301' to 'HLA-A03:01'.
    """
    if len(allele) == 5:  # Ensure the input is in the expected format
        return f"HLA-{allele[:1]}{allele[1:3]}:{allele[3:]}"
    else:
        raise ValueError(f"Invalid allele format: {allele}")

def load_mixmhcpred_results(file_path, merged_df):
    """
    Load MixMHCpred prediction results from a file and merge them into the main DataFrame.
    Assumes the file contains columns: 
    Peptide, Score_bestAllele, BestAllele, %Rank_bestAllele, Score_A0301, %Rank_A0301, Score_A0101, %Rank_A0101, etc.
    """
    print(f"Loading MixMHCpred results from {file_path}")
    try:
        df = pd.read_csv(file_path, sep="\t", comment='#')
        # Extract the peptide column
        peptide_col = df['Peptide']
        
        # Initialize a list to store temporary DataFrames
        temp_dfs = []
        
        # Iterate over all columns to extract %Rank_ columns and their corresponding alleles
        for col in df.columns:
            if re.match(r'%Rank_[A-Z]\d+', col):
                hla = col.split('_')[1]  # Extract allele name from column header
                rank_col = df[col]
                allele = reformat_allele(hla)  # Reformat the allele
                
                # Create a temporary DataFrame for the current allele
                temp_df = pd.DataFrame({
                    'peptide': peptide_col,
                    'allele': allele,
                    'MixMHCpred_%Rank': rank_col
                })
                temp_dfs.append(temp_df)
        
        # Concatenate all temporary DataFrames into a single DataFrame
        final_temp_df = pd.concat(temp_dfs, ignore_index=True)
        
        # Merge the final temporary DataFrame with the main merged_df
        merged_df = pd.merge(merged_df, final_temp_df, on=['peptide', 'allele'], how='outer')
    except Exception as e:
        print(f"Error loading MixMHCpred results from {file_path}: {e}")
    return merged_df

def load_mhcnuggets_results(file_path, merged_df):
    """
    Load MHCnuggets prediction results from a file and merge them into the main DataFrame.
    Assumes the file contains columns: Peptide,HLA,%Rank.
    """
    print(f"Loading MHCnuggets results from {file_path}")
    try:
        df = pd.read_csv(file_path, sep="\t")
        df.rename(columns={'%Rank': 'MHCnuggets_%Rank'}, inplace=True)
        df.rename(columns={'HLA': 'allele'}, inplace=True)
        df.rename(columns={'Peptide': 'peptide'}, inplace=True)
        
        # Multiply the %Rank value by 100
        df['MHCnuggets_%Rank'] = df['MHCnuggets_%Rank'] * 100
        
        merged_df = pd.merge(merged_df, df, on=['peptide', 'allele'], how='outer')
    except Exception as e:
        print(f"Error loading MHCnuggets results from {file_path}: {e}")
    return merged_df

def main():
    parser = argparse.ArgumentParser(description="Merge HLA binding prediction results from multiple tools.")
    parser.add_argument('--input', nargs='+', help="Path to results files as list")
    parser.add_argument('--mhcflurry',  help="Path to mhcflurry results file")
    parser.add_argument('--netmhcpan',  help="Path to netmhcpan results file")
    parser.add_argument('--mixmhcpred', help="Path to mixmhcpred results file")
    parser.add_argument('--mhcnuggets', help="Path to mhcnuggets results file")
    parser.add_argument('--output', default="merged_hla_prediction_results.tsv", help="Output file path")

    args = parser.parse_args()

    merged_df = pd.DataFrame(columns=['peptide', 'allele'])

    used_cols = []

    if args.input:
        for f in args.input:
            if f.startswith("mhcnuggets"):
                merged_df = load_mhcnuggets_results(f, merged_df)
                used_cols.append('MHCnuggets_%Rank')
            elif f.startswith("netmhc"):
                merged_df = load_netmhcpan_results(f, merged_df)
                used_cols.append('NetMHCpan_%Rank')
            elif f.startswith("mixmhcpred"):
                merged_df = load_mixmhcpred_results(f, merged_df)
                used_cols.append('MixMHCpred_%Rank')
            elif f.startswith("mhcflurry"):
                merged_df = load_mhcflurry_results(f, merged_df)
                used_cols.append('MHCFlurry_%Rank')
            # Add additional conditions for other tools as needed
    else:
        if args.mhcflurry:
            merged_df = load_mhcflurry_results(args.mhcflurry, merged_df)
            used_cols.append('MHCFlurry_%Rank')

        if args.netmhcpan:
            merged_df = load_netmhcpan_results(args.netmhcpan, merged_df)
            used_cols.append('NetMHCpan_%Rank')

        if args.mixmhcpred:
            merged_df = load_mixmhcpred_results(args.mixmhcpred, merged_df)
            used_cols.append('MixMHCpred_%Rank')

        if args.mhcnuggets:
            merged_df = load_mhcnuggets_results(args.mhcnuggets, merged_df)
            used_cols.append('MHCnuggets_%Rank')

    
    # Calculate the number of alleles with %Rank below thresholds
    # Ensure all %Rank columns are converted to floats and format values to 3 decimal places
    for col in used_cols:
        if col in merged_df.columns:
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').round(3)

    # Pivot the DataFrame to have one row per peptide and columns for %Rank values from each tool for each allele
    pivoted_df = merged_df.pivot_table(index='peptide', columns='allele', values=used_cols, aggfunc='first')

    # Flatten the multi-level columns
    pivoted_df.columns = [f"{col[1]}_{col[0]}" for col in pivoted_df.columns]

    # Calculate the minimum %Rank across tools
    merged_df['%Rank'] = merged_df[used_cols].min(axis=1)
    weak_binders = merged_df[merged_df['%Rank'] < 2].groupby('peptide').size()
    strong_binders = merged_df[merged_df['%Rank'] < 0.5].groupby('peptide').size()

    # Add columns for weak and strong binders
    pivoted_df['Weak_Binders_Count'] = pivoted_df.index.map(weak_binders).fillna(0).astype(int)
    pivoted_df['Strong_Binders_Count'] = pivoted_df.index.map(strong_binders).fillna(0).astype(int)

    # Determine the best binding allele (lowest %Rank) for each peptide
    best_binding_allele = merged_df.loc[merged_df.groupby('peptide')['%Rank'].idxmin()]
    pivoted_df['Best_Binding_Allele'] = pivoted_df.index.map(best_binding_allele.set_index('peptide')['allele'])

    # Determine the minimum rank for each tool for each peptide - For ML Algorithm
    tools = [col.split('_')[0] for col in used_cols]
    for tool in tools:
        tool_cols = [col for col in pivoted_df.columns if col.endswith(f'_{tool}_%Rank')]
        if tool_cols:
            pivoted_df[f'{tool}_bestRank'] = pivoted_df[tool_cols].min(axis=1)
            pivoted_df[f'{tool}_bestRankAllele'] = pivoted_df[tool_cols].idxmin(axis=1).str.split('_').str[0]

    # Save the pivoted results to a file
    pivoted_df.to_csv(args.output, sep="\t")
    print(f"Merged results saved to {args.output}")

if __name__ == "__main__":
    main()