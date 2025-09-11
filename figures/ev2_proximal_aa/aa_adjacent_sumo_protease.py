import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def annotate_bars_stacked(ax, bars1, bars2, p_values):
    for bar1, bar2, p_val in zip(bars1, bars2, p_values):
        height1 = bar1.get_height()
        height2 = bar2.get_height()
            
        max_height = max(height1, height2) + 0.1 

        x_position = (bar1.get_x() + bar2.get_x()) / 2 + bar1.get_width() / 2  # Center of the two bars

        # Add asterisk
        if p_val < 0.0001:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.25, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.5, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif p_val < 0.001:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.25, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif p_val < 0.05:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif p_val >= 0.05:
            ax.text(x_position, max_height, '', ha='center', va='bottom', fontsize=10)

def plot_aa_frequencies(input_files, output_file):
    positions = ['-2', '-1', '+1', '+2']
    
    titles = ['Trypsin SUMOylated', 'Trypsin unSUMOylated', 'Lys-C/Asp-N SUMOylated', 'Lys-C/Asp-N unSUMOylated', 'Lys-C or Lys-C/Glu-C SUMOylated', 'Lys-C or Lys-C/Glu-C unSUMOylated']
    custom_colors = [
        (0.12, 0.47, 0.71),  # dark blue for Trypsin SUMOylated
        (0.68, 0.85, 0.90),  # Light blue for Trypsin unSUMOylated
        (0.17, 0.63, 0.17),  # dark green for Lys-C/Asp-N SUMOylated
        (0.70, 0.87, 0.54),  # Light green for Lys-C/Asp-N unSUMOylated
        (0.84, 0.15, 0.16),  # dark red for Lys-C or Lys-C/Glu-C SUMOylated
        (1.0, 0.71, 0.76),   # Light red for Lys-C or Lys-C/Glu-C unSUMOylated
    ]

    y_labels = ['', '', '', '']
    
 
    amino_acids = None

    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(7.09, 7.09))  # 180mm width, 180mm height
    
    # Set font
    plt.rcParams.update({'font.size': 8, 'font.family': 'Arial'})
    axs = axs.ravel()

    for i, file_input in enumerate(input_files):
        df = pd.read_csv(file_input)
        
        norm_cols = ['Trypsin SUMOylated norm.', 'Trypsin not SUMOylated norm.',
                     'Asp-N SUMOylated norm.', 'Asp-N not SUMOylated norm.',
                     'Glu-C SUMOylated norm.', 'Glu-C not SUMOylated norm.']
        df[norm_cols] = df[norm_cols].fillna(0)
        
        pval_cols = ['Asp-N SUMOylated vs Non SUMOylated adj_pval',
                     'Glu-C SUMOylated vs Non SUMOylated adj_pval',
                     'Trypsin SUMOylated vs Non SUMOylated adj_pval']
        df[pval_cols] = df[pval_cols].fillna(1.0)
        
        # Get amino acids present in the current DataFrame
        current_amino_acids = list(df['Amino Acid'].unique())
        
        if amino_acids is None:
            amino_acids = current_amino_acids

        x = np.arange(len(current_amino_acids)) * 0.85  
        bar_width = 0.11  

        # Plot the bars
        bars1 = axs[i].bar(x - 1.5 * bar_width, df['Trypsin SUMOylated norm.'], width=bar_width, color=custom_colors[0], label=titles[0], edgecolor='black')
        bars2 = axs[i].bar(x - 0.5 * bar_width, df['Trypsin not SUMOylated norm.'], width=bar_width, color=custom_colors[1], label=titles[1], edgecolor='black')
        bars3 = axs[i].bar(x + 0.5 * bar_width, df['Asp-N SUMOylated norm.'], width=bar_width, color=custom_colors[2], label=titles[2], edgecolor='black')
        bars4 = axs[i].bar(x + 1.5 * bar_width, df['Asp-N not SUMOylated norm.'], width=bar_width, color=custom_colors[3], label=titles[3], edgecolor='black')
        bars5 = axs[i].bar(x + 2.5 * bar_width, df['Glu-C SUMOylated norm.'], width=bar_width, color=custom_colors[4], label=titles[4], edgecolor='black')
        bars6 = axs[i].bar(x + 3.5 * bar_width, df['Glu-C not SUMOylated norm.'], width=bar_width, color=custom_colors[5], label=titles[5], edgecolor='black')

        # Get the p-values and annotate bars with significance
        p_values1 = df['Trypsin SUMOylated vs Non SUMOylated adj_pval']  # First p-value column for bars1 and bars2
        p_values2 = df['Asp-N SUMOylated vs Non SUMOylated adj_pval']  # Second p-value column for bars3 and bars4
        p_values3 = df['Glu-C SUMOylated vs Non SUMOylated adj_pval']  # Third p-value column for bars5 and bars6

        # Annotate the first two bars with the first column of p-values
        annotate_bars_stacked(axs[i], bars1, bars2, p_values1)

        # Annotate the second two bars with the second column of p-values
        annotate_bars_stacked(axs[i], bars3, bars4, p_values2)

        # Annotate the third two bars with the third column of p-values
        annotate_bars_stacked(axs[i], bars5, bars6, p_values3)

        # Add horizontal line for normalization
        axs[i].axhline(y=1, linestyle='--', color='darkgray', linewidth=1)

        # Formatting
        axs[i].set_xticks(x)
        axs[i].set_xticklabels(current_amino_acids, rotation=0, fontsize=9)
        axs[i].set_yticks(np.arange(0, 6, 1))
        axs[i].set_xlabel('Lys {}'.format(positions[i]), fontweight='bold', fontsize=10)
        axs[i].set_ylabel('Normalized counts',  fontsize=10)
        axs[i].set_xlim([-0.5, x[-1] + 0.5]) 
        axs[i].text(-0.05, 1.05, y_labels[i], transform=axs[i].transAxes,
                    fontsize=18, fontweight='bold', va='center', ha='center')
        axs[i].set_ylim([0, 5.5])
        # Remove top and right spines
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
    # Add legend
    handles = [Patch(facecolor=custom_colors[i], edgecolor='black', label=titles[i]) for i in range(len(titles))]
    
    axs[0].legend(handles=handles, fontsize=11, ncol=2, bbox_to_anchor=(0.5, 1.15), loc='center', frameon=False) 

    if not os.path.exists('output'):
        os.makedirs('output')  
    
    plt.tight_layout()
  
    fig.savefig(output_file, dpi=600)
    pdf_file = output_file.replace('.png', '.pdf')
    fig.savefig(pdf_file, dpi=600)
    plt.close(fig) 
    


input_files = [
    'example_input/-2_sites_relative_to_Lys.csv',
    'example_input/-1_sites_relative_to_Lys.csv',
    'example_input/+1_sites_relative_to_Lys.csv',
    'example_input/+2_sites_relative_to_Lys.csv'
]


# Output file 
output_file = 'output/counts_of_proximal_aa_protease.png'

# Generate the plot
plot_aa_frequencies(input_files, output_file)