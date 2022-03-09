from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd

# This script plots histograms of pHlys by lysosome from the six sample
# images in the extra images of bafilomycin endpoint supplementary
# figure.

current_dir = Path.cwd()
man_dir = current_dir.parents[1]
data_path = man_dir / 'source_data' / 'bafilomycin_endpoint' / 'baf_endpoint_by_lysosome.csv'
results = pd.read_csv(data_path)

# basic plot setup
plt.style.use(man_dir / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# locate the particular recording and position
df1029 = results.loc[results['date'] == 20211029]
df1117 = results.loc[results['date'] == 20211117]
dmso1 = df1029.loc[(df1029['sample_ID'] == 5) & (df1029['image_ID'] == 2)]
dmso2 = df1029.loc[(df1029['sample_ID'] == 3) & (df1029['image_ID'] == 4)]
dmso3 = df1117.loc[(df1117['sample_ID'] == 13) & (df1117['image_ID'] == 4)]
baf1 = df1029.loc[(df1029['sample_ID'] == 6) & (df1029['image_ID'] == 1)]
baf2 = df1117.loc[(df1117['sample_ID'] == 14) & (df1117['image_ID'] == 3)]
baf3 = df1117.loc[(df1117['sample_ID'] == 16) & (df1117['image_ID'] == 1)]

# check that my histogram bins are OK
for data in [dmso1, dmso2, dmso3, baf1, baf2, baf3]:
    assert np.max(data['mean_tau_ns'].values) < 5.5
    assert np.min(data['mean_tau_ns'].values) > 1

bin_list = np.linspace(1, 5.5, num=45)
# figure setup
fig1 = plt.figure(figsize=(6,3), dpi=300)
h = [Size.Fixed(1.0), Size.Fixed(3)]
v = [Size.Fixed(0.7), Size.Fixed(1)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

# find and plot the correct data
color_list = ['#56B4E9','#009E73','#8E72B2','#D55E00','#E69F00', '#F0E442']
label_list = ['DMSO 1', 'DMSO 2', 'DMSO 3', 'BafA 1', 'BafA 2', 'BafA 3']
for i, (data, lab, col) in enumerate(zip([dmso1, dmso2, dmso3, baf1, baf2, baf3],
                                         label_list, color_list)):
    axs1.hist(data['mean_tau_ns'].values, bins=bin_list, alpha=0.7,
              label=lab, histtype=u'step', color=col, linewidth=0.75,
              zorder=100-i) # this just makes the lighter blue show up on top

# formatting
axs1.set_ylabel('# Lysosomes')
axs1.set_xlabel('Lifetime (ns)')
axs1.set_xlim(1, 5.5)
axs1.set_ylim(0, 60)
axs1.legend()
out_path = current_dir / ('combined_hist.pdf')
fig1.savefig(out_path, bbox_inches='tight', transparent=True)

plt.show()

