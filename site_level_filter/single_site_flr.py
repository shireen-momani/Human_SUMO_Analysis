import os
import pandas as pd
pd.set_option('expand_frame_repr', False)

"""
Designed for experiments with a single SUMOylation footprint.

Workflow:
    - Filter rows by False Localization Rate (FLR) threshold
    - Remove contaminants and decoy peptides
    - Collapse peptidoforms to site-level

"""
def site_based_one_input(input_file, flr_cutoff, identified_isoform, output_file):

    df_mzid_flr = pd.read_csv(input_file, sep=',')
    # Drop unnecessary columns and rename columns
    df_mzid_flr.drop(columns=['Protein', 'Protein position', 'PTM'], inplace=True)
    df_mzid_flr.rename(columns={'All_Proteins': 'Protein', 'All_PTM_protein_positions': 'Protein position', 'All_PTMs': 'PTM'}, inplace=True)
    
    # Split the 'Protein' column by ':' into a list
    df_mzid_flr['Protein'] = df_mzid_flr['Protein'].str.split(':')
    df_mzid_flr['Protein position'] = df_mzid_flr['Protein position'].str.split(':')
    
    # Explode the 'Protein' and 'Protein position' columns and create new rows
    df_protein_exploded = pd.DataFrame({
        'Protein': [item for sublist in df_mzid_flr['Protein'] for item in sublist],
        'Protein position': [item for sublist in df_mzid_flr['Protein position'] for item in sublist],
        'Peptide_mod': df_mzid_flr['Peptide_mod'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'Peptide_mod_pos': df_mzid_flr['Peptide_mod_pos'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'Binomial_final_score': df_mzid_flr['Binomial_final_score'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'Peptide': df_mzid_flr['Peptide'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'PTM positions': df_mzid_flr['PTM positions'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'PTM': df_mzid_flr['PTM'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True),
        'pA_q_value_BA': df_mzid_flr['pA_q_value_BA'].repeat(df_mzid_flr['Protein'].apply(len)).reset_index(drop=True)
    })
    
    df_protein_exploded['Protein position'] = df_protein_exploded['Protein position'].str.split(';')
    df_protein_exploded['PTM'] = df_protein_exploded['PTM'].str.split(';')
    
    # Explode the 'PTM' and 'Protein' columns and create new rows
    df_ptm_exploded = pd.DataFrame({
        'Protein': df_protein_exploded['Protein'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'Protein position': [item for sublist in df_protein_exploded['Protein position'] for item in sublist],
        'Peptide': df_protein_exploded['Peptide'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'Peptide_mod': df_protein_exploded['Peptide_mod'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'Peptide_mod_pos': df_protein_exploded['Peptide_mod_pos'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'PTM positions': df_protein_exploded['PTM positions'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'PTM': [item for sublist in df_protein_exploded['PTM'] for item in sublist],
        'Binomial_final_score': df_protein_exploded['Binomial_final_score'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'pA_q_value_BA': df_protein_exploded['pA_q_value_BA'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True)
    })
    
    df_sumo_only = df_ptm_exploded[df_ptm_exploded['PTM'] == 'Sumo']
    df = df_sumo_only[~df_sumo_only['Protein'].str.contains('DECOY|CONTAM')].copy()
    # Filter data with FLR <= flr_cutoff
    df = df[df['pA_q_value_BA'].astype(float) <= flr_cutoff]
    
    
    # Sort the data by Site Q-Value in ascending order
    df.sort_values(by='pA_q_value_BA', ascending=True, inplace=True)
    
    # Remove duplicate sites, keeping the entry with the lowest FLR
    df.drop_duplicates(subset=['Protein', 'Protein position'], keep='first', inplace=True)
    
    # Extract modified K or A from 'Peptide' column based on positions from 'PTM positions' column
    extracted_characters = df.apply(lambda row: row['Peptide'][int(row['PTM positions']) - 1], axis=1)
    
    # Create 'PTM_residue' column 
    df['PTM_residue'] = extracted_characters.apply(lambda char: 1 if char == 'A' else 0)
    df['PTM_residue'] = df['PTM_residue'].replace({1: 'A', 0: 'K'})
    
    df['Protein position'] = df['Protein position'].astype(int)
    
    df['Identified SUMO isoform'] = identified_isoform

    if not os.path.exists('outputs'):
        os.makedirs('outputs')

    df.rename(columns={'Binomial_final_score': 'BinomialScore','pA_q_value_BA': 'FLR'}, inplace=True)
    
    # Select and reorder columns
    df = df[[
        'Protein',
        'Protein position',
        'Peptide',
        'PTM_residue',
        'Identified SUMO isoform',
        'Peptide_mod_pos',
        'BinomialScore',
        'FLR'
    ]]
    
    # Write to a csv file
    df.to_csv(output_file, sep=',', index=False, header=True)
    
    # Remove Decoy PTM
    df_sites_without_decoy = df[df['PTM_residue'] != 'A']
    PTM_decoy = len(df) - len(df_sites_without_decoy)
   
    # Count of SUMO sites (without decoy PTM) and SUMOylated proteins
    df_unique_proteins = df_sites_without_decoy.drop_duplicates(subset=['Protein'], keep='first')

    print(f'\n=== Analysis Summary ===\n'
      f'FLR ≤ {flr_cutoff}\n'
      f'  - SUMO sites (decoys removed): {len(df_sites_without_decoy):,}\n'
      f'  - SUMOylated proteins: {len(df_unique_proteins):,}\n'
      f'  - Decoy PTMs removed: {PTM_decoy:,}')

    
# Example usage
site_based_one_input(
    input_file='example_input/binomial_peptidoform_collapsed_FLR.csv',
    flr_cutoff=0.05, 
    identified_isoform='SUMO3', 
    output_file='outputs/site_based_flr.csv'
)
