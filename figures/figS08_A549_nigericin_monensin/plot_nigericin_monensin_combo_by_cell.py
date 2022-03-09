from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from scipy.optimize import curve_fit
import pandas as pd

# This script will plot Scarlet-LAMP1 photon count and lifetime data
# from 2021-10-27 and 2021-12-03 in A549.
#
# This script reads in the output of
# 'nigericin_monensin_combo_by_cell.py', which processes the lifetime
# image data for the nigericin and monensin calibration. This script
# then generates a pH titration curve and fits it to a 4- parameter
# logistic function.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[1]
data_path = manuscript_path/'source_data'/'nigericin_monensin_A549'/'A549_nigericin_monensin_cell_means.csv'

# basic plot setup
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# load in the results
results = pd.read_csv(data_path)

# Calculate some rudimentary stats
gb_pH = results.groupby(['buffer_pH'])
tau_means = gb_pH['mean_tau_ns'].mean().tolist()
tau_stds = gb_pH['mean_tau_ns'].std().tolist()

fig1 = plt.figure(figsize=(3,3), dpi=300)
# generate fixed size axes, 1.5 inches
h = [Size.Fixed(1.0), Size.Fixed(1.5)]
v = [Size.Fixed(0.7), Size.Fixed(1.5)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

# fit to a 4 parameter logistic function
def logistic4(pH, min_val, hill, pKa, max_val):
    return ((min_val - max_val) / (1.0 + ((pH/pKa)**hill))) + max_val
# set up a dataframe to store the fit outputs
fit_output = pd.DataFrame(columns={'condition','model', 'temperature', 'min_val',
                                   'min_val_SD', 'max_val', 'max_val_SD', 'hill',
                                   'hill_SD', 'pKa', 'pKa_SD'})
# perform the fitting
pH_range = gb_pH['mean_tau_ns'].mean().index.tolist()
init_params = [1.7, 15, 5, 3.5]
popt, pcov = curve_fit(logistic4,
                       xdata = pH_range,
                       ydata = tau_means,
                       p0 = init_params,
                       maxfev=10000)
fit_output.at[0, 'condition'] = 'A549'
fit_output.at[0, 'temperature'] = 35
fit_output.at[0, 'model'] = 'mean_arrival'
min_val, hill, pKa, max_val = popt
fit_output.at[0, 'min_val'] = min_val
fit_output.at[0, 'max_val'] = max_val
fit_output.at[0, 'hill'] = hill
fit_output.at[0, 'pKa'] = pKa
perr = np.sqrt(np.diag(pcov))
fit_output.at[0, 'min_val_SD'] = perr[0]
fit_output.at[0, 'max_val_SD'] = perr[3]
fit_output.at[0, 'hill_SD'] = perr[1]
fit_output.at[0, 'pKa_SD'] = perr[2]

# now plot the means and the fit curves
pH_plotting = np.linspace(4, 7.5, num=500)
axs1.plot(pH_plotting, logistic4(pH_plotting, min_val, hill, pKa, max_val),
          label='', marker=None, markersize=0, color=cycle[0])

# medians +/- stdev
axs1.errorbar(pH_range, tau_means, tau_stds, linewidth=0, elinewidth=1,
              markersize=4, marker='.', capthick=1, color=cycle[0])
axs1.spines['right'].set_visible(False)
axs1.spines['top'].set_visible(False)
axs1.set_ylabel('Lifetime (ns)')
axs1.set_xlabel('Buffer pH')
axs1.set_ylim(1.6, 3.6)
axs1.set_yticks(np.linspace(1.6, 3.6, 6))
axs1.set_title('A549 Cells')

fig1.savefig('A549_nigericin_monensin_means_whole_cell.pdf',
             bbox_inches='tight', transparent=True)
fit_output.to_csv('A549_4PL_fits.csv')

plt.show()



