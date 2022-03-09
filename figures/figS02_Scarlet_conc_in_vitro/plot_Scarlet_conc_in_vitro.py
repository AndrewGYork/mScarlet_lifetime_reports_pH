from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import linregress

# This script plots the in vitro data taken for purified Scarlet-His.
# Data were analyzed with SPT (mean arrival time). Median and stdev of
# the mean arrival time for each image was calculated by
# generate_in_vitro_summary_table.py.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[1]
data_path_1 = manuscript_path / 'source_data' / 'in_vitro_characterization' 
data_path = data_path_1 / 'combined_in_vitro_mean_arrival_time.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# load in the results
results = pd.read_csv(data_path)

# now let's generate the graph
fig1 = plt.figure(figsize=(5, 2.1))
h = [Size.Fixed(0.5), Size.Fixed(1.5)]
v = [Size.Fixed(0.45), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

conc = results.loc[results['series'] == 'Scarlet_conc']
conc_median = conc.groupby(['Scarlet_mgmL'])['median_tau_ns'].mean()

# convert from mg/mL to mM
mg_mL = conc_median.index.values
molar = mg_mL / 27952 # theoretical mass = 27952 Da

fig1.axes[0].plot(molar, conc_median.values)

# axis formatting
fig1.axes[0].set_ylim(3, 4)
fig1.axes[0].set_ylabel('Lifetime (ns)')
fig1.axes[0].set_xlabel('mScarlet Concentration (M)')
fig1.axes[0].set_xscale('log')
fig1.axes[0].set_xticks([1e-7, 1e-6, 1e-5, 1e-4])

fig1.savefig('scarlet_conc.png')

plt.show()
