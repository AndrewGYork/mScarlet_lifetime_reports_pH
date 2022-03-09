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

#####################################################
## ANALYSIS OF VISCOSITY AND REFRACTIVE INDEX IN GLYCEROL
glycerol = results.loc[results['series'] == 'glycerol']
glyc_summary = glycerol.groupby(['%_glycerol'])['median_tau_ns'].mean()
# literature viscosity: http://www.met.reading.ac.uk/~sws04cdw/viscosity_calc.html
# literature RI: "physical properties of glycerine and its solutions"
lit_viscosity = np.array([0.935, 1.2021, 1.890,
                          5.357]) # in cP: 0, 10, 25, 50% glyc w/w at 23C
lit_RI = np.array([1.33303, 1.34481, 1.36404,
                   1.39809]) # 0, 10, 25, 50% glyc w/w at 20C

# set up the figures for plotting
fig3, axs3 = plt.subplots(1, 3, figsize=(7, 2))
plt.subplots_adjust(wspace=0.5)

# panel a: lifetime vs. % glycerol
fig3.axes[0].plot(glyc_summary.index.values, glyc_summary.values)
fig3.axes[0].set_xlabel('% glycerol (w/v)')
fig3.axes[0].set_ylabel('Lifetime (ns)')
fig3.axes[0].set_ylim(3, 4)
fig3.axes[0].set_title('% Glycerol')

# panel b: lifetime vs. viscosity for corresponding % glycerol
fig3.axes[1].plot(lit_viscosity, glyc_summary.values)
fig3.axes[1].set_xlabel('Dynamic viscosity (cP)')
fig3.axes[1].set_ylabel('Lifetime (ns)')
fig3.axes[1].set_ylim(3, 4)
fig3.axes[1].set_title('Viscosity')

# panel c: strickler-berg predicted, 1/lifetime vs RI^2
fig3.axes[2].plot(lit_RI**2, 1/glyc_summary.values, linewidth=0)
fig3.axes[2].set_xlabel('n$^2$')
fig3.axes[2].set_ylabel('1 / Lifetime (ns$^{-1}$)')
fig3.axes[2].set_xlim(1.7, 2.3)
fig3.axes[2].set_yticks(np.linspace(0.27, 0.30, 4))
fig3.axes[2].set_ylim(0.27, 0.30)
fig3.axes[2].set_title('Refractive Index')

# add a line of best fit to panel c
slope, inter, rval, *other = linregress(
    lit_RI**2, 1/glyc_summary.values)
plt.text(1.95, 0.285, "y=%0.3fx+%0.3f" % (slope, inter))
plt.text(1.95, 0.280, "r$^2$=%0.3f" % rval **2)
x_range = [lit_RI[0]**2, lit_RI[-1]**2]
y_range = [slope * x + inter for x in x_range]
fig3.axes[2].plot(x_range, y_range, color=cycle[0], label='')

# add panel labels
plt.figtext(0.05, 0.9, 'A', weight='bold', fontsize=11)
plt.figtext(0.34, 0.9, 'B', weight='bold', fontsize=11)
plt.figtext(0.62, 0.9, 'C', weight='bold', fontsize=11)

fig3.savefig('glycerol.png', bbox_inches='tight')

plt.show()
