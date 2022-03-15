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

def logistic4(pH, min_val, hill, pKa, max_val):
    return ((min_val - max_val) / (1.0 + ((pH/pKa)**hill))) + max_val

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

# set up a dataframe to store the fit outputs
fit_output = pd.DataFrame(columns={'condition','model', 'temperature', 'min_val',
                                   'min_val_SD', 'max_val', 'max_val_SD', 'hill',
                                   'hill_SD', 'pKa', 'pKa_SD'})

##################################################################
## ANALYSIS OF pH AND TEMP DEPENDENCE in PBS, carmody's
# with mean arrival time models
# Take the mean of both replicates for each condition
fig1, axs1 = plt.subplots(1, 2, sharey=True, figsize=(6, 3))
pH_temp = results.loc[results['series'] == 'pH_temp']
summary = pH_temp.groupby(['buffer', 'temp_C', 'pH'])['median_tau_ns'].mean()
buffer_list = ['PBS', 'PBS', 'Carmody_IB', 'Carmody_IB']
temp_list = [23, 35, 23, 35]
num_fits = 0
init_params = [1.7, 15, 5, 3.5]
for i, (buff, temp) in enumerate(zip(buffer_list, temp_list)):
    pH_range = summary.at[buff, temp].index.tolist()
    popt, pcov = curve_fit(logistic4,
                           xdata = pH_range,
                           ydata = summary.at[buff, temp].tolist(),
                           p0 = init_params,
                           maxfev=10000)
    fit_output.at[num_fits, 'condition'] = buff
    fit_output.at[num_fits, 'temperature'] = temp
    fit_output.at[num_fits, 'model'] = 'mean_arrival'
    min_val, hill, pKa, max_val = popt
    fit_output.at[num_fits, 'min_val'] = min_val
    fit_output.at[num_fits, 'max_val'] = max_val
    fit_output.at[num_fits, 'hill'] = hill
    fit_output.at[num_fits, 'pKa'] = pKa
    perr = np.sqrt(np.diag(pcov))
    fit_output.at[num_fits, 'min_val_SD'] = perr[0]
    fit_output.at[num_fits, 'max_val_SD'] = perr[3]
    fit_output.at[num_fits, 'hill_SD'] = perr[1]
    fit_output.at[num_fits, 'pKa_SD'] = perr[2]
    num_fits += 1

    # now plot the means and the fit curves
    if buff == 'PBS':
        this_ax = fig1.axes[0]
    elif buff == 'Carmody_IB':
        this_ax = fig1.axes[1]
    if temp == 23:
        this_color = cycle[0]
    elif temp == 35:
        this_color = cycle[1]
    this_ax.plot(pH_range,
                 summary.at[buff, temp].tolist(),
                 label=u"%02.f\N{DEGREE SIGN}C" % temp,
                 linewidth=0)
    pH_plotting = np.linspace(4, 7.5, num=500)
    this_ax.plot(pH_plotting, logistic4(pH_plotting, min_val, hill, pKa, max_val),
                 label='', marker=None, markersize=0,
                 color=this_color)

# some plot formatting
plt.legend()
fig1.axes[0].set_title('PBS')
fig1.axes[1].set_title('Carmody Imaging Buffer')
fig1.axes[0].set_ylim(1.8, 3.8)
fig1.axes[0].set_yticks(np.linspace(1.8, 3.8, 6))
for ax in fig1.axes:
    ax.set_xlabel('pH')
    ax.set_ylabel('Lifetime (ns)')

# add panel labels
plt.figtext(0.05, 0.9, 'A', weight='bold', fontsize=11)
plt.figtext(0.5, 0.9, 'B', weight='bold', fontsize=11)

# save the results
fit_output.to_csv('in_vitro_4PL_fits.csv')
fig1.savefig('PBS_Carm_temperature.png', bbox_inches='tight')

plt.show()
