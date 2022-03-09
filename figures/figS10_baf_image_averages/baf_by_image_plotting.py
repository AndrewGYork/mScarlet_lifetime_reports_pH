from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from matplotlib import ticker
import pandas as pd
import seaborn as sns

# This script plots the "by cell"/"by image" data from the bafilomycin
# time series and endpoint analysis as a complement to the by-lysosome
# analysis present elsewhere in the manuscript.

current_dir = Path.cwd()
data_path = current_dir.parents[1] / 'source_data'
still_path = data_path / 'bafilomycin_endpoint' / 'baf_endpoint_by_cell.csv'
ts_path = data_path / 'bafilomycin_time_series' / 'baf_time_series_by_image.csv'
still_df = pd.read_csv(still_path)
ts_df = pd.read_csv(ts_path)
ts_df['time'] = (ts_df['frame_ID']-1)*5 - 8

# basic plot setup
plt.style.use(current_dir.parents[1] / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


# FIGURE 1: TIME SERIES - ALL CELLS EXPLICITLY SHOWN
fig1 = plt.figure(figsize=(4,4), dpi=300)
# generate fixed size axes, 2 inches square
h = [Size.Fixed(1.0), Size.Fixed(2)]
v = [Size.Fixed(0.7), Size.Fixed(2)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

# identify the various recordings & the condition they are from
data = [ts_df.loc[ids, 'mean_tau_ns'].values for
        ids in ts_df.groupby(['date','position']).groups.values()]
times = [ts_df.loc[ids, 'time'].values for
         ids in ts_df.groupby(['date','position']).groups.values()]
drug_list = [ts_df.loc[ids, 'condition'].values[0] for
             ids in ts_df.groupby(['date','position']).groups.values()]
time_scale = range(-8, 52, 5)

# make the plots
labelled_D = False
labelled_B = False
for i, (time_scale, trace) in enumerate(zip(times, data)):
    if drug_list[i] == 'DMSO':
        plot_color = cycle[0]
        if labelled_D:
            lab = ''
        else:
            lab = 'DMSO'
            labelled_D = True
    else:
        plot_color = cycle[1]
        if labelled_B:
            lab = ''
        else:
            lab = 'BafA'
            labelled_B = True
        
    axs1.plot(time_scale, trace, markersize=2, linewidth=1, marker='o',
              color=plot_color, label=lab)
    
time_label = range(-10, 60, 10)
for i, ax in enumerate(fig1.axes):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Lifetime (ns)')
    ax.set_xlabel('Time (min)')
    ax.set_xticks(time_label)
    ax.set_ylim(1.8, 3.8)
    ax.set_yticks(np.arange(1.8, 3.80001, 0.4))
plt.legend()
    
fig1.savefig('u2os_time_series_by_cell.pdf', bbox_inches='tight',
             transparent=True)

# FIGURE 2: TIME SERIES - AVERAGED ACROSS ALL CELLS
fig2 = plt.figure(figsize=(4,4), dpi=300)
# generate fixed size axes, 2 inches square
h = [Size.Fixed(1.0), Size.Fixed(2)]
v = [Size.Fixed(0.7), Size.Fixed(2)]
divider = Divider(fig2, (0, 0, 1, 1), h, v, aspect=False)
axs2 = fig2.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
sns.lineplot(data=ts_df, x='time', y='mean_tau_ns', hue='condition',
             estimator=np.median, ci='sd')
   
for i, ax in enumerate(fig2.axes):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Lifetime (ns)')
    ax.set_xlabel('Time (min)')
    ax.set_ylim(1.8, 3.8)
    ax.set_xlim(-10, 50)
    ax.xaxis.set_major_locator(ticker.LinearLocator(7))
    ax.set_yticks(np.arange(1.8, 3.80001, 0.4))
plt.legend()
    
fig2.savefig('u2os_time_series_all_cell_avg.pdf', bbox_inches='tight',
             transparent=True)

# FIGURE 3: ENDPOINT IMAGES: DMSO VS BAF
fig3 = plt.figure(figsize=(4,4), dpi=300)
# generate fixed size axes, 2 inches square
h = [Size.Fixed(1.0), Size.Fixed(2)]
v = [Size.Fixed(0.7), Size.Fixed(2)]
divider = Divider(fig3, (0, 0, 1, 1), h, v, aspect=False)
axs3 = fig3.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

sns.stripplot(data=still_df, x='drug', y='mean_tau_ns', size=3)
    
for i, ax in enumerate(fig3.axes):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Lifetime (ns)')
    ax.set_ylim(1.8, 3.8)
    ax.set_yticks(np.arange(1.8, 3.80001, 0.4))
    ax.set_xlabel('Condition')
    ax.set_xticklabels(['DMSO', '100 nM BafA'])
axs3.set_title('5 hr Incubation')
    
fig3.savefig('u2os_baf_endpoint_by_cell.pdf', bbox_inches='tight',
             transparent=True)

# count the total number of cells in each condition
init = ts_df.loc[ts_df['frame_ID'] == 1]
print('TS counts:', init.groupby('condition').size())
print()
print('Still counts:', still_df.groupby('drug').size())
print()
print('Mean lifetimes:', still_df.groupby('drug')['mean_tau_ns'].mean())

plt.show()

