from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import seaborn as sns

# This script plots data calculated for each lysosome in the baf endpoint
# dataset. See script "baf_endpoint_by_lysosome_grinstein.py"
# for more information on segmentation. This script is separate to avoid
# doing the expensive computation (segmentation by lysosome) multiple
# times during visualization.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[1]
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

# FIGURE 1 - ANALYSIS BY ZONES A LA GRINSTEIN
fig1 = plt.figure(figsize=(4,4), dpi=300)
# generate 1:1 aspect ratio, 2 inch axes
h = [Size.Fixed(1.0), Size.Fixed(2)]
v = [Size.Fixed(0.7), Size.Fixed(2)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
results = results.astype({'zone': 'int32'}) # bit of a hack to get display right
sns.violinplot(data=results, x='zone', y='mean_tau_ns', hue='drug',
               inner='quartile')
plt.legend(frameon=True)
axs1.set_xlabel('Cellular Region')
axs1.set_ylabel('Lifetime (ns)')

fig1.savefig('baf_dmso_by_zone.pdf', bbox_inches='tight', transparent=True)

# calculate summary statistics by zone, incl. sample size
gb_zone = results.groupby(['drug','zone'])
print(gb_zone['mean_tau_ns'].count())
print(gb_zone['mean_tau_ns'].mean())


plt.show()

