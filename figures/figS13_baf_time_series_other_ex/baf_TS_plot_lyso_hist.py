from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd

# This script plots histograms of pHlys by lysosome from two particular
# recordings (DMSO 10/29 01, pos 5; Baf 10/29 01 pos 9) for use in main
# text figure 3.

current_dir = Path.cwd()
man_dir = current_dir.parents[1]
data_path = man_dir / 'source_data' / 'bafilomycin_time_series' / 'baf_time_series_individ_lyso_results.csv'
results = pd.read_csv(data_path)

# basic plot setup
plt.style.use(man_dir / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# locate the particular recording and position
df1029 = results.loc[results['date'] == 20211029]
df1215 = results.loc[results['date'] == 20211215]
dmso1 = df1029.loc[df1029['position'] == 2]
dmso2 = df1215.loc[df1215['position'] == 3]
baf1 = df1029.loc[df1029['position'] == 8]
baf2 = df1215.loc[df1215['position'] == 7]
# check that my histogram bins are OK
assert np.max(dmso1['mean_tau_ns'].values) < 5.5
assert np.max(baf1['mean_tau_ns'].values) < 5.5
assert np.max(dmso2['mean_tau_ns'].values) < 5.5
assert np.max(baf2['mean_tau_ns'].values) < 5.5
assert np.min(dmso1['mean_tau_ns'].values) > 1
assert np.min(baf1['mean_tau_ns'].values) > 1
assert np.min(dmso2['mean_tau_ns'].values) > 1
assert np.min(baf2['mean_tau_ns'].values) > 1
bin_list = np.linspace(1, 5.5, num=45)
for x in range(12):
    # figure setup
    fig1 = plt.figure(figsize=(3,3), dpi=300)
    h = [Size.Fixed(1.0), Size.Fixed(1)]
    v = [Size.Fixed(0.7), Size.Fixed(1)]
    divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
    axs1 = fig1.add_axes(divider.get_position(),
                         axes_locator=divider.new_locator(nx=1, ny=1))
    # find and plot the correct data
    this_frame_dmso1 = dmso1.loc[dmso1['frame_ID'] == x+1]
    axs1.hist(this_frame_dmso1['mean_tau_ns'].values, bins=bin_list, alpha=0.7,
              label='DMSO_1', color='#56B4E9')
    this_frame_dmso2 = dmso2.loc[dmso2['frame_ID'] == x+1]
    axs1.hist(this_frame_dmso2['mean_tau_ns'].values, bins=bin_list, alpha=0.7,
              label='DMSO_2', color='#009E73')
    this_frame_baf1 = baf1.loc[baf1['frame_ID'] == x+1]
    axs1.hist(this_frame_baf1['mean_tau_ns'].values, bins=bin_list, alpha=0.7,
              label='BafA_1', histtype=u'step', linewidth=0.75, color='#D55E00')
    this_frame_baf2 = baf2.loc[baf2['frame_ID'] == x+1]
    axs1.hist(this_frame_baf2['mean_tau_ns'].values, bins=bin_list, alpha=0.7,
              label='BafA_2', histtype=u'step', linewidth=0.75, color='#E69F00')
    # formatting
    axs1.set_ylabel('# Lysosomes')
    axs1.set_xlabel('Lifetime (ns)')
    axs1.set_xlim(1, 5.5)
    axs1.set_ylim(0, 100)
    axs1.set_title('%d min' % (x * 5 - 8)) # time 0 = when baf was added
    axs1.legend()
    out_path = current_dir / ('combined_hist_t%d.pdf'%x)
    fig1.savefig(out_path, bbox_inches='tight',
                 transparent=True)

plt.show()

