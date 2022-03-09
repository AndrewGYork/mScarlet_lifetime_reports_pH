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
temp = results.loc[results['date'] == 20211029]
temp1 = temp.loc[temp['recording'] == 1]
baf = temp1.loc[temp1['position'] == 9]
dmso = temp1.loc[temp1['position'] == 5]
assert np.max(dmso['mean_tau_ns'].values) < 5
assert np.max(baf['mean_tau_ns'].values) < 5
assert np.min(dmso['mean_tau_ns'].values) > 1
assert np.min(baf['mean_tau_ns'].values) > 1
bin_list = np.linspace(1, 5, num=40)
for x in range(12):
    # figure setup
    fig1 = plt.figure(figsize=(3,3), dpi=300)
    h = [Size.Fixed(1.0), Size.Fixed(1)]
    v = [Size.Fixed(0.7), Size.Fixed(1)]
    divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
    axs1 = fig1.add_axes(divider.get_position(),
                         axes_locator=divider.new_locator(nx=1, ny=1))
    # find and plot the correct data
    this_frame_dmso = dmso.loc[dmso['frame_ID'] == x+1]
    axs1.hist(this_frame_dmso['mean_tau_ns'].values, bins=bin_list, alpha=0.5,
              label='DMSO')
    this_frame_baf = baf.loc[baf['frame_ID'] == x+1]
    axs1.hist(this_frame_baf['mean_tau_ns'].values, bins=bin_list, alpha=0.5,
              label='BafA')
    # formatting
    axs1.set_ylabel('# Lysosomes')
    axs1.set_xlabel('Lifetime (ns)')
    axs1.set_xlim(1, 5)
    axs1.set_ylim(0, 60)
    axs1.set_title('%d min' % (x * 5 - 8)) # time 0 = when baf was added
    axs1.legend()
    out_path = current_dir / 'one_fov' / ('baf_dmso_hist_oneFOV_t%d.pdf'%x)
    fig1.savefig(out_path, bbox_inches='tight',
                 transparent=True)

plt.show()

