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

def logistic4(pH, min_val, hill, pKa, max_val):
    return ((min_val - max_val) / (1.0 + ((pH/pKa)**hill))) + max_val

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[2]
data_path_1 = manuscript_path / 'source_data' / 'in_vitro_characterization' 
data_path = data_path_1 / 'combined_in_vitro_mean_arrival_time.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# load in the results
results = pd.read_csv(data_path)

#########################################################################
###FIGURE 5. In vitro in PBS at 35C (this alone for main text
fig5 = plt.figure(figsize=(3, 3))
# generate fixed size axes, 1.5 inches
h = [Size.Fixed(1.0), Size.Fixed(1.5)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig5, (0, 0, 1, 1), h, v, aspect=False)
axs5 = fig5.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))
# now let's find the data we want to plot
pH_temp = results.loc[results['series'] == 'pH_temp']
summary = pH_temp.groupby(['buffer', 'temp_C', 'pH'])['median_tau_ns'].mean()
PBS_35 = summary.at['PBS', 35]
init_params = [1.7, 15, 5, 3.5]
# generate the fit again for plotting
popt, pcov = curve_fit(logistic4,
                       xdata = PBS_35.index.tolist(),
                       ydata = PBS_35.tolist(),
                       p0 = init_params,
                       maxfev=10000)
min_val, hill, pKa, max_val = popt

# now plot the means and the fit curves
axs5.plot(PBS_35.index.tolist(), PBS_35.tolist(), linewidth=0)
pH_plotting = np.linspace(4, 7.5, num=500)
axs5.plot(pH_plotting, logistic4(pH_plotting, min_val, hill, pKa, max_val),
          label='', marker=None, markersize=0, color=cycle[0])

# some plot formatting
axs5.set_title('In Vitro')
axs5.set_ylim(1.8, 3.8)
axs5.set_yticks(np.linspace(1.8, 3.8, 6))
for ax in fig5.axes:
    ax.set_xlabel('pH')
    ax.set_ylabel('Lifetime (ns)')

# save the results
fig5.savefig('PBS_35C_for_main_text.pdf', bbox_inches='tight',
             transparent=True)

plt.show()
