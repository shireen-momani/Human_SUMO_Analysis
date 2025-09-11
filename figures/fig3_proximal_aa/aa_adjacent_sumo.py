import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def annotate_bars_stacked(ax, bars1, bars4, p_values):
    for bar1, bar4, p_val in zip(bars1, bars4, p_values):
        height1 = bar1.get_height()
        height4 = bar4.get_height()
        max_height = max(height1, height4) + 0.1 

        x_position = bar1.get_x() + bar1.get_width() / 2 

        # Add the significance asterisk vertically
        if p_val < 0.0001:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.25, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.5, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif p_val < 0.001:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
            ax.text(x_position, max_height + 0.25, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif p_val < 0.05:
            ax.text(x_position, max_height, '*', ha='center', va='bottom', fontsize=10, fontweight='bold')

def plot_aa_frequencies(file_inputs, output_file):
    positions = ['-2', '-1', '+1', '+2']
    custom_colors = ['#FFD700', '#C0C0C0', '#CD7F32', '#FF0000']
    titles = ['Gold', 'Gold-Silver', 'Gold-Silver-Bronze', 'Non-SUMOylated']

    amino_acids = None

    # Set figure size to (180mm x 180mm)
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(7.09, 7.09), edgecolor='black')  # 180mm width, 180mm height
    axs = axs.ravel()

    plt.rcParams.update({'font.size': 8, 'font.family': 'Arial'})

    for i, file_input in enumerate(file_inputs):
        df = pd.read_csv(file_input)
        if amino_acids is None:
            amino_acids = list(df['Amino Acid'].unique())

        x = np.arange(len(amino_acids))  
        bar_width = 0.2 

        # Plot the bars
        bars1 = axs[i].bar(x - 1.5 * bar_width, df['Gold norm.'], width=bar_width, color=custom_colors[0], label=titles[0], edgecolor='black')
        bars2 = axs[i].bar(x - 0.5 * bar_width, df['Gold-Silver norm.'], width=bar_width, color=custom_colors[1], label=titles[1], edgecolor='black')
        bars3 = axs[i].bar(x + 0.5 * bar_width, df['Gold-Silver-Bronze norm.'], width=bar_width, color=custom_colors[2], label=titles[2], edgecolor='black')
        bars4 = axs[i].bar(x + 1.5 * bar_width, df['non-SUMOylated norm.'], width=bar_width, color=custom_colors[3], label=titles[3], edgecolor='black')

        annotate_bars_stacked(axs[i], bars1, bars4, df['Gold vs non-SUMOylated adj_pval'])

        axs[i].axhline(y=1, linestyle='--', color='darkgray', linewidth=1) # Add horizontal line for normalization
        axs[i].set_xticks(x)
        axs[i].set_xticklabels(amino_acids, rotation=0, fontsize=9)
        axs[i].set_yticks(np.arange(0, 8, 1))
        axs[i].tick_params(axis='y', labelsize=9)
        axs[i].set_xlabel('Lys {}'.format(positions[i]), fontweight='bold', fontsize=10)
        axs[i].set_ylabel('Normalized counts', fontsize=10) # Add y-axis label to each subplot
        axs[i].set_ylim([0, 4])
        axs[i].set_xlim([-1, len(amino_acids) - 0.5])
        axs[i].spines['top'].set_visible(False)  # Remove top and right spines
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['left'].set_linewidth(1.5) # Set left and bottom spine line width
        axs[i].spines['bottom'].set_linewidth(1.5)

    # Add legend
    handles = [Patch(facecolor=custom_colors[i], edgecolor='black', label=titles[i]) for i in range(len(titles))]
    axs[0].legend(handles=handles, fontsize=10.5, ncol=4, bbox_to_anchor=(0.5, 1.15), loc='center', frameon=False)

    if not os.path.exists('output'):
        os.makedirs('output')  

    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) 
    fig.savefig('{}.png'.format(output_file), 
                dpi=600, 
                bbox_inches='tight') 
    plt.savefig("{}.pdf".format(output_file),
              format='pdf',
              dpi=600,
              bbox_inches='tight')



file_paths_updated = [
    'example_input/-2_sites_relative_to_lys.csv',
    'example_input/-1_sites_relative_to_lys.csv',
    'example_input/+1_sites_relative_to_lys.csv',
    'example_input/+2_sites_relative_to_lys.csv'
]


output_file = 'output/Counts_of_proximal_AA'

# Generate the plot
plot_aa_frequencies(file_paths_updated, output_file)