from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd
import seaborn as sns
from scipy.stats import linregress

# This script plots data calculated for each lysosome in the baf endpoint
# dataset. See script "baf_endpoint_by_lysosome_grinstein.py"
# for more information on segmentation. This script is separate to avoid
# doing the expensive computation (segmentation by lysosome) multiple
# times during visualization.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[2]
data_path_1 = manuscript_path/'source_data'/'bafilomycin_endpoint'
data_path = data_path_1/'baf_endpoint_by_lysosome.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
# fix the things that ruin violin plots
mpl.rcParams['lines.markersize'] = 0
mpl.rcParams['lines.marker'] = '.'

#read in the data
results = pd.read_csv(data_path)

# FIGURE 1 - BAF VS DMSO
baf = results.loc[results['drug'] == '100nM_bafA1']
dmso = results.loc[results['drug'] == '0.1%_DMSO']
fig1 = plt.figure(figsize=(6,3), dpi=300)
# generate 2:1 aspect ratio
h = [Size.Fixed(1.0), Size.Fixed(3)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
assert min(results['mean_tau_ns'].values) > 1.5
assert max(results['mean_tau_ns'].values) < 5
tau_edges = np.linspace(1.5, 5, 35)
plt.hist(dmso['mean_tau_ns'].values, bins=tau_edges, label='DMSO', alpha=0.5)
plt.hist(baf['mean_tau_ns'].values, bins=tau_edges, label='100nM BafA', alpha=0.5)
axs1.set_xlabel('Lifetime (ns)')
axs1.set_ylabel('Lysosome Count')
plt.legend()

fig1.savefig('baf_dmso_lyso_hist.pdf', bbox_inches='tight',
             transparent=True)

# print the SD of the distribution under each condition
print('standard deviation by condition')
print(results.groupby('drug')['mean_tau_ns'].std())

plt.show()

