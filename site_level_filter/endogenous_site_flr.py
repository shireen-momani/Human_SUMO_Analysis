import os
import pandas as pd
from Bio import SeqIO

pd.set_option('expand_frame_repr', False)

# To include additional column if SUMO sites is on peptide C-terminal lysines
def determine_peptide_c_terminal(file_path):
    df = pd.read_csv(file_path, sep=',')
    records = df.to_dict('records')
    c_terminal = []
    for i in range(len(records)):
        if len(records[i]['Peptide']) - records[i]['PTM positions'] == 0:
            c_terminal.append(True)
        else:
            c_terminal.append(False)
    df['peptide C-terminal'] = c_terminal
    return df
"""
Designed for endogenous datasets generated with Lys-C/Asp-N digestion strategy.

Workflow:
    - Filter rows by False Localization Rate (FLR) threshold
    - Remove contaminants and decoy peptides
    - Remove Asp-N derived SUMO sites on peptide C-terminal lysines unless the next residue is aspartic acid (D) or glutamic acid (E).
    - Collapse peptidoforms to site-level

Sample input files are provided for testing:
   - Peptidoform-site mzidFLR output file
   - Human proteome FASTA
"""
def site_based_Endogenous(input_file, proteome, flr_cutoff, identified_isoform, output_file):
    df = determine_peptide_c_terminal(input_file)
    df.drop(columns=['Protein', 'Protein position', 'PTM'], inplace=True)
    df.rename(columns={'All_Proteins': 'Protein', 'All_PTM_protein_positions': 'Protein position', 'All_PTMs': 'PTM'}, inplace=True)
    # Split the 'Protein' column by ':' into a list
    df['Protein'] = df['Protein'].str.split(':') 
    df['Protein position'] = df['Protein position'].str.split(':')
    # Explode the 'Protein' and 'Protein position' columns and create new rows
    df_protein_exploded = pd.DataFrame({
        'Protein': [item for sublist in df['Protein'] for item in sublist],
        'Protein position': [item for sublist in df['Protein position'] for item in sublist],
        'Peptide_mod_pos': df['Peptide_mod_pos'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'Peptide_mod': df['Peptide_mod'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'Binomial_final_score': df['Binomial_final_score'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'Peptide': df['Peptide'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'PTM positions': df['PTM positions'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'PTM': df['PTM'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'peptide C-terminal': df['peptide C-terminal'].repeat(df['Protein'].apply(len)).reset_index(drop=True),
        'pA_q_value_BA': df['pA_q_value_BA'].repeat(df['Protein'].apply(len)).reset_index(drop=True)
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
        'peptide C-terminal': df_protein_exploded['peptide C-terminal'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True),
        'pA_q_value_BA': df_protein_exploded['pA_q_value_BA'].repeat(df_protein_exploded['PTM'].apply(len)).reset_index(drop=True)
    })
    df_sumo_only = df_ptm_exploded[df_ptm_exploded['PTM'] == 'Sumo']
    df = df_sumo_only[~df_sumo_only['Protein'].str.contains('DECOY|CONTAM')].copy()
    # group by 'Protein and 'Protein position' columns and get 'Protein position' column as lists in rows
    df_Peptide_mod = df.groupby(['Protein', 'Protein position'])['Peptide_mod'].apply(list).reset_index(name='Peptide_mod')
    # group by 'Protein columns and get Protein position' and 'Peptide_mod columns as lists in rows
    df_Peptide_mod = df_Peptide_mod.groupby('Protein').apply(lambda x: x[['Protein position','Peptide_mod']].values.tolist()).reset_index(name='col')
    # A smaller fasta file containing SUMOylated proteins
    fastaFile = proteome
    df_txt = df.drop_duplicates(subset=['Protein'], keep='last')#Remove duplicate id
    id_file = df_txt['Protein']#take only protein id column
    wanted = id_file.tolist()#put it in a list
    outputfile = 'UP000005640_9606_less_seqs.fasta'
    records = (r for r in SeqIO.parse(fastaFile, format= 'fasta') if r.id in wanted)
    count = SeqIO.write(records, outputfile, 'fasta')
    print('Saved %i records from %s to %s' % (count, fastaFile, outputfile))
    if count < len(wanted):
        print('Warning %i IDs not found in %s' % (len(wanted)-count, fastaFile))
    fastaFile2 = 'UP000005640_9606_less_seqs.fasta'
    dict_df_Peptide_mod = df_Peptide_mod.to_dict('records')
    records = list(SeqIO.parse(fastaFile2, format='fasta'))
    data = []
    for i in range(len(dict_df_Peptide_mod)):
        for e in dict_df_Peptide_mod[i]['col']:
            for j in e[1]:
                acc = dict_df_Peptide_mod[i]['Protein']
                pos = int(e[0])
                record = next((r for r in records if r.id == acc), None)
                if record is not None:
                    if pos < len(record.seq):
                        data.append({
                            'Protein': acc,
                            'Protein position': pos,
                            'Peptide_mod': j,
                            'next_residue': record.seq[pos]
                        })
                    elif int(pos) == len(record.seq):
                        data.append({
                            'Protein': acc,
                            'Protein position': pos,
                            'Peptide_mod': j,
                            'next_residue': '-'
                        })
    df_next_residue = pd.DataFrame(data)
    df['Protein position'] = df['Protein position'].astype(int)
    # Merge the two dataframes
    df = pd.merge(df, df_next_residue, on=['Protein', 'Protein position','Peptide_mod'], how='outer')
    df = df[df['pA_q_value_BA'].astype(float) <= flr_cutoff]

    # Sort the data by Site Q-Value in ascending order
    df.sort_values(by='pA_q_value_BA', ascending=True, inplace=True)
    # Keep first occurrence of each group    
    df.drop_duplicates(subset=['Protein', 'Protein position'], keep='first', inplace=True)
    
    # Remove rows where the value in column 'peptide C-terminal' is TRUE and the value in column 'next_residue' is not equal to 'E or D'
    df_filtered = df.drop(df[(df['peptide C-terminal'] == True) & ~((df['next_residue'] == 'E') | (df['next_residue'] == 'D'))].index)
    count = len(df) - len(df_filtered)
    # Extract modified K or A from 'Peptide' column based on positions from 'PTM positions' column
    extracted_characters = df_filtered.apply(lambda row: row['Peptide'][row['PTM positions'] - 1], axis=1)
    df_filtered['PTM_residue'] = extracted_characters.apply(lambda char: 1 if char == 'A' else 0)
    df_filtered['PTM_residue'] = df_filtered['PTM_residue'].replace({1: 'A', 0: 'K'})

    df_filtered['Identified SUMO isoform'] = identified_isoform

    if not os.path.exists('outputs'):
        os.makedirs('outputs')
    df_filtered.rename(columns={'Binomial_final_score': 'BinomialScore','pA_q_value_BA': 'FLR'}, inplace=True)
    
    # Select and reorder columns
    df_final = df_filtered[[
        'Protein',
        'Protein position',
        'Peptide',
        'PTM_residue',
        'Identified SUMO isoform',
        'Peptide_mod_pos',
        'BinomialScore',
        'FLR'
    ]]
    df_final.to_csv(output_file, sep=',', index=False, header=True)
    # Remove Decoy PTM
    df_final_count = df_final[df_final['PTM_residue'] != "A"]
    PTM_decoy = len(df_final) - len(df_final_count)

    # Print count of SUMO sites (without decoy PTM) and SUMOylated proteins
    df_final_prot = df_final_count.drop_duplicates(subset=['Protein'], keep='first')
    print(f'FLR <= {flr_cutoff}')
    print(f'  - SUMO sites (decoys removed): {len(df_final_count):,}')
    print(f'  - SUMOylated proteins: {len(df_final_prot):,}')
    print(f'  - Decoy PTMs removed: {PTM_decoy:,}')
    print(f'  - Asp-N filter: {count:,} sites removed (peptide C-terminal lysines not followed by D or E)')



# Example usage
site_based_Endogenous(
    input_file='example_input/binomial_peptidoform_collapsed_FLR_endogenous.csv',
    proteome='example_input/UP000005640_9606.fasta', 
    flr_cutoff=0.05, 
    identified_isoform='SUMO2/3', 
    output_file='outputs/endogenous_site_based_flr.csv'
)

# Remove temporary files
os.remove('UP000005640_9606_less_seqs.fasta')