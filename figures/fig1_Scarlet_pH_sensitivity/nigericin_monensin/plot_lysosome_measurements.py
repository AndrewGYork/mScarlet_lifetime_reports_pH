from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import seaborn as sns
import tifffile as tf
import pandas as pd

# Run this file after nigericin_monensin_combo_by_lysosome.py, which
# will calculate the mean and standard deviation # of photons and mean
# arrival time per lysosomal region. Note: it will do some filtering for
# size/# of photons - be sure you are accounting for this in thinking
# about the lysosomal pH distribution.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[2]
data_path = manuscript_path / 'source_data' / 'nigericin_monensin_U2OS' / 'nigericin_monensin_lyso_2pxmin.csv'

#read in the data
results = pd.read_csv(data_path)

# import style sheet and overwrite the parameters that will mess with violinplot
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
mpl.rcParams['lines.markersize'] = 0
mpl.rcParams['lines.marker'] = '.'

# generate a figure
fig1 = plt.figure(figsize=(3,3), dpi=300)
# generate fixed size axes, 1.5 inches square
h = [Size.Fixed(1.0), Size.Fixed(1.5)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

sns.violinplot(data=results, x='buffer_pH', y='mean_tau_ns', inner='quartile',
               color='#248cd6', linewidth=0.75)

# labels and ranges
axs1.set_ylim(1,5)
axs1.set_xticks(np.arange(0, 8, 2))
axs1.set_xticklabels([4, 5, 6, 7])
axs1.set_ylabel('Lifetime (ns)')
axs1.set_xlabel('Buffer pH')
axs1.set_title('Individual Lysosomes')

fig1.savefig('lysosome_calibration_violin.pdf', bbox_inches='tight',
             transparent=True)

# print out the number of lysosomes per measurement for captions
print('Number of lysosomes per pH')
print(results.groupby('buffer_pH').size())

# print the median by pH for inclusion in main text
print('median by pH')
print(results.groupby('buffer_pH')['mean_tau_ns'].median())

# print the standard deviation by pH for inclusion in main text
print('standard deviation by pH')
print(results.groupby('buffer_pH')['mean_tau_ns'].std())

plt.show()
