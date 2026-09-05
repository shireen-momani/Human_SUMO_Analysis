import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import seaborn as sns
import numpy as np
import warnings
from matplotlib.patches import Rectangle

pd.set_option("expand_frame_repr", False)
import matplotlib
matplotlib.rcParams['ps.fonttype'] = 3
matplotlib.rcParams['pdf.fonttype'] = 3
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica']

# Create one figure with 2 columns (each subdivided into 2 rows).
fig = plt.figure(figsize=(180/25.4, 6.0), dpi=600)

# Outer grid with more spacing
outer_gs = gridspec.GridSpec(1, 2, width_ratios=[48, 52], wspace=0.45, figure=fig)

# Create subgrids for left and right columns
left_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_gs[0], hspace=0.75)
right_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_gs[1], hspace=0.75)

# Create subplots
ax1 = fig.add_subplot(left_gs[0])  # top-left (A)
ax3 = fig.add_subplot(left_gs[1])  # bottom-left (C)
ax2 = fig.add_subplot(right_gs[0])  # top-right (B)
ax4 = fig.add_subplot(right_gs[1])  # bottom-right (D)


# Label subplots as A, B, C, D
ax1.text(-0.18, 1.15, 'A', transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left', clip_on=False)
ax2.text(-0.18, 1.15, 'B', transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left', clip_on=False)
ax3.text(-0.18, 1.15, 'C', transform=ax3.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left', clip_on=False)
ax4.text(-0.18, 1.15, 'D', transform=ax4.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left', clip_on=False)

# Set the theme for the plots
sns.set_style("white")
sns.set_context("paper")


# Subplot A
groups = ['Gold', 'Silver', 'Bronze']
count_k = [8639, 14050, 13032]
count_a = [46, 243, 2080] 

color_k = '#4A90E2'  #  Blue
color_a = '#E67E22'  #  Orange

ax1.bar(groups, count_k, color=color_k, label="K", edgecolor="none", width=0.7)
ax1.bar(groups, count_a, bottom=count_k, color=color_a, label="A", edgecolor="none", width=0.7)

# Add text labels 
# For K  
for i, (group, k_val) in enumerate(zip(groups, count_k)):
    ax1.text(i, k_val / 2, str(k_val), ha='center', va='center', 
             color='white', size=8, fontweight='normal')

# For A 
for i, (group, k_val, a_val) in enumerate(zip(groups, count_k, count_a)):
    ax1.text(i, k_val + a_val + 300, str(a_val), ha='center', va='bottom', 
             color='black', size=8, fontweight='normal')

ax1.set_ylabel("Number of Sites", fontsize=10, fontweight='bold')
ax1.set_xticks(np.arange(len(groups)))
ax1.set_xticklabels(groups, fontsize=9)
ax1.tick_params(axis='x', labelsize=9)
ax1.tick_params(axis='y', labelsize=9)
ax1.set_yticks([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000])
ax1.set_yticklabels(['0', '2k', '4k', '6k', '8k', '10k', '12k', '14k', '16k'])
ax1.legend(frameon=False, fontsize=9, loc='upper left', bbox_to_anchor=(0, 1.0))
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_linewidth(1.2)
ax1.spines['bottom'].set_linewidth(1.2)
ax1.grid(False)


# Subplot B 
categories = ['Gold', 'Gold-Silver', 'Gold-Silver-Bronze', 'Non-SUMOylated']
proportions = [
   
    (4601/8639) * 100,
    (10318/22695) * 100,
    (15540/35795) * 100,
    (19273/65010) * 100
]

ax2.bar(categories, proportions, color="#34495E", edgecolor="none", width=0.5)

# Add asterisks
for i in range(3):  # Only the first three bars (Gold, Gold-Silver, Gold-Silver-Bronze)
    ax2.text(i, proportions[i] + 1.5, '*\n*\n*', ha='center', va='bottom', fontsize=8, fontweight='bold', linespacing=0.5)

proteome_percent = (210746/653505) * 100
ax2.axhline(y=proteome_percent, linestyle='--', color='#F39C12', linewidth=1.0)

ax2.set_ylabel('Disorder > 0.5 (%)', fontsize=10, fontweight='bold')
ax2.set_xticks(np.arange(len(categories)))
ax2.set_xticklabels(categories, fontsize=9, rotation=12, ha='right')
ax2.tick_params(axis='x', labelsize=9)
ax2.tick_params(axis='y', labelsize=9)
ax2.set_ylim(0, 60)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_linewidth(1.2)
ax2.spines['bottom'].set_linewidth(1.2)
ax2.grid(False)


# Subplot C
sets = ["Gold", "Gold-Silver", "Gold-Silver-Bronze", "Non-SUMOylated"]
structures = ["Total","α-helix", "β-strand", "Coiled coil", "Turn"]

df_results = pd.read_csv('example_input/secondary_structure_proportions.csv')
pivot = df_results.pivot(index='Filename', columns='Secondary Structure', values='Percentage')

filenames = [
    "Gold",
    "Gold_Silver",
    "Gold_Silver_Bronze",
    "Non-SUMOylated_Lysines_on_peptides_Below_1%_FDR"
]

values = {
    "α-helix": pivot.loc[filenames, "Helix"].tolist(),
    "β-strand": pivot.loc[filenames, "Beta strand"].tolist(),
    "Coiled coil": pivot.loc[filenames, "Coiled coil"].tolist(),
    "Turn": pivot.loc[filenames, "Turn"].tolist(),
    "Total": pivot.loc[filenames, ["Helix", "Beta strand", "Coiled coil", "Turn"]].sum(axis=1).tolist()
}

color_list = [
    '#2C3E50',   # Total (dark charcoal)
    '#3498DB',   # α-helix (light blue)
    '#16A085',   # β-strand (teal green)
    '#9B59B6',   # Coiled coil (purple)
    '#E74C3C'    # Turn (red)
]
line_order = ["Total", "α-helix", "β-strand", "Coiled coil", "Turn"] 

for struct in line_order:
    ydata = values[struct]
    idx = structures.index(struct)
    ax3.plot(sets, ydata, marker='o', linestyle='-', 
             label=struct, linewidth=1.5, markersize=4, 
             color=color_list[idx], markeredgecolor='black', markeredgewidth=0.5)

# Baseline lines with matching colors 
ax3.axhline(y=19.02, linestyle='--', color='#2C3E50', linewidth=0.8)
ax3.axhline(y=8.81, linestyle='--', color='#3498DB', linewidth=0.8)
ax3.axhline(y=4.1, linestyle='--', color='#16A085', linewidth=0.8)
ax3.axhline(y=5.26, linestyle='--', color='#9B59B6', linewidth=0.8)
ax3.axhline(y=0.85, linestyle='--', color='#E74C3C', linewidth=0.8)

# Add asterisks to the "Total" line 
total_y = values["Total"]
for i in range(3): 
    if i == 0:  # Gold
        y_pos = total_y[i] + 3.0
    elif i == 1:  # Gold-Silver
        y_pos = total_y[i] + 2.0
    else:  # Gold-Silver-Bronze
        y_pos = total_y[i] + 1.5 
    ax3.text(i, y_pos, '*\n*\n*', ha='center', va='bottom', fontsize=8, fontweight='bold', color='black', linespacing=0.5)
ax3.set_ylabel("Secondary Structure (%)", fontsize=10, fontweight='bold')
ax3.set_xticks(np.arange(len(sets)))
ax3.set_xticklabels(sets, fontsize=9, rotation=12, ha='right')
ax3.tick_params(axis='x', labelsize=9)
ax3.tick_params(axis='y', labelsize=9)
ax3.legend(
    fontsize=8,
    loc='upper center',
    bbox_to_anchor=(0.55, 1.12),
    ncol=3,
    frameon=False,
    columnspacing=1.2
)
ax3.set_ylim(0, 30)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['left'].set_linewidth(1.2)
ax3.spines['bottom'].set_linewidth(1.2)
ax3.grid(False)


# Subplot D
positions = ["-2", "-1", "K", "+1", "+2"]
gold = np.array([0.066097, 0.064709, 0.057175, 0.064083, 0.062500]) * 100
gold_silver = np.array([0.063289, 0.060772, 0.053508, 0.061352, 0.062009]) * 100
gold_silver_bronze = np.array([0.062137, 0.058021, 0.049807, 0.060374, 0.060288]) * 100
unSUMO = np.array([0.056680, 0.057020, 0.046897, 0.059652, 0.060670]) * 100

bar_width = 0.022
x_adjusted = np.array([-0.19, -0.095, 0.0, 0.095, 0.19])
offsets = np.linspace(-1.5 * bar_width, 1.5 * bar_width, 4)

gold_color = "#F1C40F"       # Gold
silver_color = "#95A5A6"     # Silver
bronze_color = "#E67E22"     # Dark orange
nonsumo_color = "#E74C3C"    # Bright red

ax4.bar(x_adjusted + offsets[0], gold, bar_width, label="Gold", color=gold_color, edgecolor="black", linewidth=0.3)
ax4.bar(x_adjusted + offsets[1], gold_silver, bar_width, label="Gold-Silver", color=silver_color, edgecolor="black", linewidth=0.3)
ax4.bar(x_adjusted + offsets[2], gold_silver_bronze, bar_width, label="Gold-Silver-Bronze", color=bronze_color, edgecolor="black", linewidth=0.3)
ax4.bar(x_adjusted + offsets[3], unSUMO, bar_width, label="Non-SUMOylated", color=nonsumo_color, edgecolor="black", linewidth=0.3)

significant_positions = {
    ("Gold", "-2", "*", offsets[0]),
    ("Gold", "-1", "*", offsets[0]),
    ("Gold", "K", "*", offsets[0]),
    ("Gold-Silver", "-2", "*", offsets[1]),
    ("Gold-Silver", "K", "*", offsets[1]),
    ("Gold-Silver-Bronze", "-2", "*", offsets[2])
}

for category, pos, sig, offset in significant_positions:
    pos_index = positions.index(pos)
    x_pos = x_adjusted[pos_index] + offset
    y_val = {
        "Gold": gold,
        "Gold-Silver": gold_silver,
        "Gold-Silver-Bronze": gold_silver_bronze,
        "UnSUMOylated": unSUMO
    }[category][pos_index]

    ax4.text(x_pos, y_val + 0.2, sig, ha="center", va="bottom", fontsize=8, fontweight="bold")

ax4.set_xticks(x_adjusted)
ax4.set_xticklabels(positions, fontsize=9)
ax4.set_ylabel("Disease variants (%)", fontsize=10, fontweight='bold')
ax4.set_xlabel("Position relative to lysine (residues)", fontsize=10, fontweight='bold')
ax4.tick_params(axis='y', labelsize=9)
ax4.tick_params(axis='x', labelsize=9)
ax4.legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=2, frameon=False, columnspacing=1.0)
ax4.set_ylim(0, 10)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['left'].set_linewidth(1.2)
ax4.spines['bottom'].set_linewidth(1.2)
ax4.grid(False)


if not os.path.exists('output'):
        os.makedirs('output')
        
plt.savefig('output/fig1_GBS_site_counts_structural_context_combined.pdf', format='pdf', dpi=600, bbox_inches='tight', facecolor='white')
plt.savefig('output/fig1_GBS_site_counts_structural_context_combined.png', format='png', dpi=600, bbox_inches='tight', facecolor='white')
# plt.show()