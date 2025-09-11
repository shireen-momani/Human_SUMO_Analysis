import matplotlib.pyplot as plt
from collections import Counter
import os
import pandas as pd

non_sumo_gold_relative_distribution =  [(0, 91.12266585699105), (1, 7.53377865492637), 
                                        (2, 1.1310156368604827), (3, 0.18217701533323213), 
                                        (4, 0.026567481402763018), (5, 0.003795354486109003)]
sumo_relative_distribution =  [(0, 83.93332561639079), (1, 13.68213913647413), 
                               (2, 2.0951499016089827), (3, 0.26623451788401437), 
                               (4, 0.023150827642088204)]
non_sumo_relative_distribution = [(0, 93.49330872173512), (1, 5.608367943393324), 
                                  (2, 0.7768035686817413), (3, 0.1045992924165513), 
                                  (4, 0.015382248884786957), (5, 0.0015382248884786955)]

# Font parameters
import matplotlib
matplotlib.rcParams['ps.fonttype'] = 3
matplotlib.rcParams['pdf.fonttype'] = 3
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4

# Create figure and axes
fig, ax = plt.subplots(figsize=(3.5, 2.8), dpi=600)
marker_size = 5
line_width = 1.2
alpha_value = 0.9


# SUMOylated sites in Gold set
sumo_x = [item[0] for item in sumo_relative_distribution if item[0] >= 1]
sumo_y = [item[1] for item in sumo_relative_distribution if item[0] >= 1]
ax.plot(
    sumo_x,
    sumo_y,
    marker='D',                    
    markersize=marker_size,
    markeredgecolor='white',
    markeredgewidth=0.3,
    linestyle='-',
    linewidth=line_width,
    color='#0066CC',
    alpha=alpha_value,
    label='SUMOylated (Gold)',
    zorder=3
)


# Non-SUMOylated sites in Gold set
non_sumo_gold_x = [item[0] for item in non_sumo_gold_relative_distribution if item[0] >= 1]
non_sumo_gold_y = [item[1] for item in non_sumo_gold_relative_distribution if item[0] >= 1]
ax.plot(
    non_sumo_gold_x,
    non_sumo_gold_y,
    marker='s',                     
    markersize=marker_size-0.5,
    markeredgecolor='white',
    markeredgewidth=0.3,
    linestyle='-',
    linewidth=line_width,
    color='#009E73',
    alpha=alpha_value,
    label='Non-SUMOylated (Gold)',
    zorder=2
)

# Non-SUMOylated sites at any FLR
non_sumo_x = [item[0] for item in non_sumo_relative_distribution if item[0] >= 1]
non_sumo_y = [item[1] for item in non_sumo_relative_distribution if item[0] >= 1]
ax.plot(
    non_sumo_x,
    non_sumo_y,
    marker='o',                     
    markersize=marker_size-0.5,
    markeredgecolor='white',
    markeredgewidth=0.3,
    linestyle='-',
    linewidth=line_width,
    color="#E69200" ,
    alpha=alpha_value,
    label='Non-SUMOylated (All)',
    zorder=1
)


# x-axis ticks to show only values ≥ 1
max_x = max(max(non_sumo_x) if non_sumo_x else 0, 
            max(sumo_x) if sumo_x else 0, 
            max(non_sumo_gold_x) if non_sumo_gold_x else 0)
ax.set_xticks(range(1, max_x + 1))
ax.set_xlim(0.5, max_x + 0.5)

ax.set_xlabel("Number of phosphosites within ±10 residues", fontsize=10, fontweight='bold')
ax.set_ylabel("Relative frequency\nof lysine residues (%)", fontsize=10, fontweight='bold')
max_y = max(max(sumo_y) if sumo_y else 0,
            max(non_sumo_y) if non_sumo_y else 0,
            max(non_sumo_gold_y) if non_sumo_gold_y else 0)
ax.set_ylim(0, max_y * 1.1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_position(('outward', 5))
ax.spines["bottom"].set_position(('outward', 5))
ax.grid(True, axis='y', linestyle='-', alpha=0.15, linewidth=0.5)
ax.set_axisbelow(True)
ax.tick_params(axis='x', labelsize=9)
ax.tick_params(axis='y', labelsize=9)

legend = ax.legend(
    loc='upper right',
    frameon=False,
    fontsize=9,
    handlelength=1.5,
    handletextpad=0.5,
    borderpad=0.3,
    columnspacing=0.5,
    labelspacing=0.3
)

fig.tight_layout(pad=0.5)

if not os.path.exists('output'):
        os.makedirs('output')

plt.savefig("output/phosphosites_flanking_lysine_distribution.png",
            dpi=600, 
            bbox_inches="tight",
            transparent=False,
            facecolor='white')
plt.savefig("output/phosphosites_flanking_lysine_distribution.pdf",
            dpi=600,
            bbox_inches="tight",
            transparent=False,
            facecolor='white')
