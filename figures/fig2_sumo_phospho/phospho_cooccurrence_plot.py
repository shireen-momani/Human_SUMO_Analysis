import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # use a non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Load pre-calculated data
sumo_relative_distribution = [(0, 83.93332561639079), (1, 13.68213913647413), 
                              (2, 2.0951499016089827), (3, 0.26623451788401437), 
                              (4, 0.023150827642088204)]
non_sumo_relative_distribution = [(0, 93.49330872173512), (1, 5.608367943393324), 
                                  (2, 0.7768035686817413), (3, 0.1045992924165513), 
                                  (4, 0.015382248884786957), (5, 0.0015382248884786955)]
positions = [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
sumo_proportions = [4.318322023442319, 4.441702652683529, 
                    5.367057371992598, 4.009870450339297, 
                    4.256631708821715, 6.909315237507712, 
                    4.009870450339297, 5.120296113510179, 
                    5.737199259716225, 0.9253547193090685, 
                    5.798889574336829, 3.51634793337446, 
                    5.798889574336829, 4.133251079580505, 
                    6.539173349784084, 5.922270203578038, 
                    5.1819864281307835, 6.47748303516348, 
                    6.230721776681062, 5.305367057371992]
non_sumo_proportions = [5.790010193679918, 5.99388379204893, 
                        6.992864424057084, 6.055045871559633, 
                        6.564729867482161, 5.953109072375128, 
                        4.770642201834862, 5.4841997961264015, 
                        4.016309887869521, 0.9785932721712538, 
                        3.9551478083588174, 2.4872579001019366, 
                        4.994903160040775, 4.994903160040775, 
                        4.729867482161061, 5.300713557594292, 
                        5.524974515800204, 5.178389398572885, 
                        4.954128440366973, 5.28032619775739]
proteome_proportions = [2.5702231597646277, 2.8644387698456755, 
                        3.3362939935605644, 3.4750749417120015, 
                        3.663817031197957, 4.2522482513600535, 
                        4.424336627067836, 4.879538137004552, 
                        4.785167092261575, 5.123792605751083, 
                        5.834351060286444, 5.739980015543466, 
                        5.401354502053958, 5.5234817364272235, 
                        5.245919840124348, 5.623404019096259, 
                        5.50682802264905, 5.467969357166648, 
                        5.667813922504719, 5.41800821583213]
reject_list = [False, False, False,  True,  True, False, False, 
               False,  True, False,  True, False, False, False,  
               True, False, False, False, False, False]
# Plotting style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.dpi'] = 300

# Color scheme
sumo_color = '#1f77b4'  # Blue
non_sumo_color = '#E69200'  # Orange 
baseline_color = '#6b7280'  # Grey

# Create figure with (180mm width)
fig = plt.figure(figsize=(7.09, 3.54))
gs = gridspec.GridSpec(1, 2, width_ratios=[30, 70], wspace=0.3)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

marker_size = 5
line_width = 1.5

sumo_x = [item[0] for item in sumo_relative_distribution if item[0] >= 1]
sumo_y = [item[1] for item in sumo_relative_distribution if item[0] >= 1]


ax1.text(-0.30, 1.1, 'A', transform=ax1.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')
ax2.text(-0.11, 1.1, 'B', transform=ax2.transAxes,
         fontsize=14, fontweight='bold', va='top', ha='left')

ax1.plot(sumo_x, sumo_y, marker='D', linestyle='-', linewidth=line_width,
         color=sumo_color, label='SUMOylated (Gold set)', markersize=marker_size,
         markeredgecolor='white', markeredgewidth=0.8,
         zorder=3)

non_sumo_x = [item[0] for item in non_sumo_relative_distribution if item[0] >= 1]
non_sumo_y = [item[1] for item in non_sumo_relative_distribution if item[0] >= 1]

ax1.plot(non_sumo_x, non_sumo_y, marker='o', linestyle='-', linewidth=line_width,
         color=non_sumo_color, label='Non-SUMOylated', markersize=marker_size,
         markeredgecolor='white', markeredgewidth=0.8,
         zorder=2)

ax1.set_xlabel('Number of phosphosites', fontsize=10, fontweight='bold')
ax1.set_ylabel('Relative frequency of lysine residues (%)', fontsize=10, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_linewidth(1.2)
ax1.spines['left'].set_linewidth(1.2)

# Add grid
ax1.grid(True, axis='y', alpha=0.15, linestyle='-', linewidth=0.5)
ax1.set_axisbelow(True)

ax1.set_ylim(-0.5, 14.5)
ax1.set_xlim(0.7, 5.3)
ax1.set_xticks(range(1, 6))

ax1.tick_params(which='major', length=3, width=0.8, labelsize=8)
ax1.tick_params(which='minor', length=2, width=0.6)

# Positional enrichment bar plot
x_bar_positions = []
for i, pos in enumerate(positions):
    if pos < 0:
        x_bar_positions.append(i)
    else:
        x_bar_positions.append(i + 1) 

x = np.array(x_bar_positions)
width = 0.4

# SUMO bars
ax2.bar(
    x - width/2,
    sumo_proportions,
    width=width,
    color=sumo_color,
    edgecolor='white',
    linewidth=0.5,
    label='SUMOylated (Gold set)',
)

# Non-SUMO bars
ax2.bar(
    x + width/2,
    non_sumo_proportions,
    width=width,
    color=non_sumo_color,
    edgecolor='white',
    linewidth=0.5,
    label='Non-SUMOylated',
)

# Proteome proportions baseline
ax2.fill_between(
    x,
    proteome_proportions,
    color='#808080',
    alpha=0.25,
    zorder=0
)
# Create modified x-axis with space for K at position 0
x_ticks = []
x_labels = []
for i, pos in enumerate(positions):
    if pos < 0:
        x_ticks.append(i)
        x_labels.append(str(pos))
    else:
        x_ticks.append(i + 1) 
        x_labels.append(f'+{pos}')

# Add K position between -1 and +1
k_position = 9.5 + 0.5 
x_ticks.insert(10, k_position)
x_labels.insert(10, 'K')

ax2.set_xticks(x_ticks)
ax2.set_xticklabels(x_labels)

for i, label in enumerate(ax2.get_xticklabels()):
    if label.get_text() == 'K':
        label.set_fontweight('bold')
ax2.set_xlabel('Position relative to lysine (residues)', fontsize=10, fontweight='bold')
ax2.set_ylabel('Relative frequency of phosphosites (%)', fontsize=10, fontweight='bold')
ax2.set_ylim(0, max(sumo_proportions + non_sumo_proportions + proteome_proportions)+1)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_linewidth(1.2)
ax2.spines['left'].set_linewidth(1.2)
ax2.grid(True, axis='y', alpha=0.15, linestyle='-', linewidth=0.5)
ax2.set_axisbelow(True)

legend2 = ax2.legend(loc='upper right', bbox_to_anchor=(0.98, 1.1), frameon=False, fontsize=9,
                    handlelength=1.5, handletextpad=0.6, ncol=1)

# Add major and minor ticks
ax2.tick_params(which='major', length=3, width=0.8, labelsize=8)
ax2.tick_params(which='minor', length=2, width=0.6)

significant_positions = []
for i, (pos, reject) in enumerate(zip(positions, reject_list)):
    if reject:
        significant_positions.append(i)

for idx in significant_positions:
    max_height = max(sumo_proportions[idx], non_sumo_proportions[idx])
    x_pos = x_bar_positions[idx]
    
    ax2.text(
        x_pos,
        max_height + 0.08,
        '*',
        ha='center',
        va='bottom',
        color='black',
        fontsize=12,
        fontweight='bold'
    )


if not os.path.exists('output'):
        os.makedirs('output')

plt.savefig('output/phospho_sumo_co_occurrence_subplot.png', dpi=600, bbox_inches='tight', 
            facecolor='white', edgecolor='none', format='png')
plt.savefig('output/phospho_sumo_co_occurrence_subplot.pdf', bbox_inches='tight', 
            facecolor='white', edgecolor='none', format='pdf')
plt.close('all')
