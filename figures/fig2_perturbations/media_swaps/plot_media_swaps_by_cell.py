from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn

# This script will plot Scarlet-LAMP1 photon count and lifetime data
# from the media swap data. Processing of raw images was
# conducted by media_swaps_by_cell.py, which lives with the raw data.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[2]
data_path_1 = manuscript_path/'source_data'/'media_swaps'
data_path = data_path_1/'media_swaps_cell_means.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
# fix the things that ruin seaborn
mpl.rcParams['lines.markersize'] = 0
mpl.rcParams['lines.marker'] = '.'

# load in the results
results = pd.read_csv(data_path)
DMEM = results.loc[results['imaging_media'] != 'HPLM']

## basic diagnostic plot - by condition
fig1 = plt.figure(figsize=(3, 3), dpi=300)
h = [Size.Fixed(1.0), Size.Fixed(1.5)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
plot_order = ['gibco_gmax_fresh', 'gibco_fresh_gln', 'corning_gln',
              'gibco_gmax', 'gibco_base']
label_list = ['GMax', 'Fresh Gln', 'Orig. Gln',
              'Spent GMax', 'DMEM Base']
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
sns.stripplot(data=DMEM, x='imaging_media', y='mean_tau_ns',
              ax=fig1.axes[0], order=plot_order, size=2.5)
pp = sns.pointplot(data=DMEM, x='imaging_media', y='mean_tau_ns',
              ax=fig1.axes[0], order=plot_order,
              ci="sd", estimator=np.median, scale=0.5,
              join=False, errwidth=1, capsize=0.1, color='k')
plt.setp(pp.lines, zorder=100)
plt.setp(pp.collections, zorder=100, label="")

# some formatting
axs1.set_ylim(2.0, 3.6)
axs1.set_yticks(np.linspace(2.0, 3.6, 5))
axs1.set_xticklabels(label_list, rotation=45, ha='right')
for i, ax in enumerate(fig1.axes):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Lifetime (ns)')

# some statistics, more for show than anything as data are visually interpretable
data = [DMEM.loc[ids, 'mean_tau_ns'].values for
        ids in DMEM.groupby('imaging_media').groups.values()]

KW = kruskal(*data)
posthocs = posthoc_dunn(DMEM, val_col='mean_tau_ns', group_col='imaging_media',
                        p_adjust='holm')
print(KW)
posthocs.to_csv('media_swaps_dunn_posthocs.csv')

# add significance metrics to the plot (note that the values will update
# as the values above change, but the selection of which lines to draw
# was done manually via inspection of the posthocs table)
plt.hlines([3.2, 3.4], [0, 0], [2, 3], colors='k')
plt.text(1, 3.21,
         'p=%s' % format(posthocs.at[plot_order[0], plot_order[2]], '.1e'),
         ha='center', va='bottom', fontsize=8)
plt.text(1.5, 3.41,
         'p=%s' % format(posthocs.at[plot_order[0], plot_order[3]], '.1e'),
         ha='center', va='bottom', fontsize=8)

fig1.savefig('media_swaps_by_cell.pdf', bbox_inches='tight',
             transparent=True)

# also print out the median lifetime in each group
print(DMEM.groupby('imaging_media')['mean_tau_ns'].median())

plt.show()



