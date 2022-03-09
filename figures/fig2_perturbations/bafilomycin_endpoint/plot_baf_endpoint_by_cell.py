from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd
import seaborn as sns

# This script will plot Scarlet-LAMP1 photon count and lifetime data
# from the endpoint bafilomycin data. Processing of raw images was
# conducted by baf_endpoint_by_cell.py, which lives with the raw data.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[2]
data_path_1 = manuscript_path/'source_data'/'bafilomycin_endpoint'
data_path = data_path_1/'baf_endpoint_by_cell.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# load in the results
results = pd.read_csv(data_path)

## basic diagnostic plot - by condition
fig1 = plt.figure(figsize=(3, 3), dpi=300)
h = [Size.Fixed(1.0), Size.Fixed(1.5)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
sns.stripplot(data=results, x='drug', y='mean_tau_ns',
              ax=fig1.axes[0], order=['0.1%_DMSO', '100nM_bafA1'],
              size=2.5)
axs1.set_ylim(1.9, 3.5)
axs1.set_yticks(np.linspace(1.9, 3.5, 5))
for i, ax in enumerate(fig1.axes):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylabel('Lifetime (ns)')

fig1.savefig('baf_vs_DMSO_by_cell.png', bbox_inches='tight')


plt.show()



