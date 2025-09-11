import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def process_sumo_data(datasets):
    
    def read_and_process_csv(dataset):
        file = f'example_input/{dataset}_site_based_flr.csv'
        df = pd.read_csv(file, sep=',')
        df = df[['Protein', 'Protein position', 'PTM_residue', 'Peptide_mod_pos', 'FLR', 'BinomialScore', 'Identified SUMO isoform']]
        df.columns = [
            'Protein',
            'Protein position',
            'PTM_residue',
            f'{dataset}_Peptide_mod_pos',
            f'{dataset}_FLR',
            f'{dataset}_BinomialScore',
            f'{dataset}_Identified_SUMO_isoform'
        ]
        return df

    # Read and merge all datasets
    df_list = [read_and_process_csv(dataset) for dataset in datasets]
    merged_df = df_list[0]
    for i in range(1, len(df_list)):
        merged_df = pd.merge(merged_df, df_list[i], on=['Protein', 'Protein position', 'PTM_residue'], how='outer')

    # For each site, count datasets with FLR below 1%
    condition = [f'{dataset}_FLR' for dataset in datasets]
    merged_df['Count of datasets below 1% FLR'] = (merged_df[condition] < 0.01).sum(axis=1)

    # Divide into Gold, Silver, Bronze sets
    """
    The Gold set includes sites identified in two or more datasets with an FLR below 1%. 
    The Silver set includes sites found in one dataset with an FLR below 1%. 
    The Bronze set includes sites that do not meet the criteria for inclusion in the Gold or Silver sets but still exhibit an FLR below 5%
    """
    conditions = [
    merged_df['Count of datasets below 1% FLR'] >= 2,
    merged_df['Count of datasets below 1% FLR'] == 1,
    merged_df['Count of datasets below 1% FLR'] == 0
    ]

    choices = ['Gold', 'Silver', 'Bronze']

    merged_df['PTM_FLR_category'] = np.select(conditions, choices, default='NA')

    if not os.path.exists('output'):
        os.makedirs('output')  
    
    # Make 'PTM_FLR_category' the fourth column
    cols = merged_df.columns.tolist()
    
    # Move 'PTM_FLR_category' to the fourth position
    cols.insert(3, cols.pop(cols.index('PTM_FLR_category')))
    merged_df = merged_df[cols]

    # Remove column
    merged_df = merged_df.drop(columns=['Count of datasets below 1% FLR'])

    merged_df.to_csv('output/human_sumo_gold_silver_bronze.csv', sep=',', index=False, header=True)
    

# Run function
datasets = ['dataset_1','dataset_2','dataset_3'] # We analysed thirteen dataset but for demo purpose, we only include three datasets here
process_sumo_data(datasets)
