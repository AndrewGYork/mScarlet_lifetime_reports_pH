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
# generate_in_vitro_summary_table.py in the source_data folder.

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

# Generate the figure
fig2, axs2 = plt.subplots(2, 2, sharey=True, figsize=(4, 4))
plt.subplots_adjust(hspace=0.6, wspace=0.2)

# subplot a - ionic strength
IS_order = [0.05, 0.5, 1, 1.5, 2, 5]
ionic = results.loc[results['series'] == 'PBS_ionic_strength']
ionic_higher = ionic.loc[ionic['PBS_strength'] > 0.3]
fig2.axes[0].plot(ionic_higher['PBS_strength'], ionic_higher['median_tau_ns'])
fig2.axes[0].set_title('Ionic Strength')
fig2.axes[0].set_xlabel('PBS Dilution')
fig2.axes[0].set_ylabel('Lifetime (ns)')
fig2.axes[0].set_ylim(3.0, 4.0)

# subplot b - KCl
kcl = results.loc[results['series'] == 'KCl']
fig2.axes[1].plot(kcl['KCl_mM'], kcl['median_tau_ns'])
fig2.axes[1].set_title('KCl')
fig2.axes[1].set_xlabel('KCl (mM)')
fig2.axes[1].set_ylabel('Lifetime (ns)')

# subplot c - MgCl2
mgcl2 = results.loc[results['series'] == 'MgCl2']
fig2.axes[2].plot(mgcl2['MgCl2_mM'], mgcl2['median_tau_ns'])
fig2.axes[2].set_title('MgCl$_{2}$')
fig2.axes[2].set_xlabel('MgCl$_{2}$ (mM)')
fig2.axes[2].set_ylabel('Lifetime (ns)')

# subplot d - CaCl2
cacl2 = results.loc[results['series'] == 'CaCl2']
fig2.axes[3].plot(cacl2['CaCl2_mM'], cacl2['median_tau_ns'])
fig2.axes[3].set_title('CaCl$_{2}$')
fig2.axes[3].set_xlabel('CaCl$_{2}$ (mM)')
fig2.axes[3].set_ylabel('Lifetime (ns)')

# add panel labels
plt.figtext(0.01, 0.9, 'A', weight='bold', fontsize=11)
plt.figtext(0.5, 0.9, 'B', weight='bold', fontsize=11)
plt.figtext(0.01, 0.425, 'C', weight='bold', fontsize=11)
plt.figtext(0.5, 0.425, 'D', weight='bold', fontsize=11)

fig2.savefig('ionic_dependencies.png', bbox_inches='tight')

plt.show()
