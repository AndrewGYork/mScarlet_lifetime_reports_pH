import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from pathlib import Path

# Run this script after extract_decays.py
# It generates a figure with normalized lifetime decays at 8 pHs in vitro.

# read in the decays parsed by the other script
# axis order is pH, time
raw_data = np.loadtxt('2021-11-18_PBS_35C_decays.csv', dtype='float64')
pH_range = [4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5]

# basic plot setup
manuscript_path = Path.cwd().parents[1]
plt.style.use(manuscript_path / 'figures' / 'default.mplstyle')
colormap = plt.get_cmap('plasma')
# get rid of some settings we don't want
mpl.rcParams['lines.markersize'] = 0
mpl.rcParams['lines.marker'] = '.'
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['lines.linewidth'] = 0.75

# SEMILOG Y AXIS
fig1 = plt.figure(figsize=(2,2), dpi=300)
# generate fixed size axes, 1 inch square
h = [Size.Fixed(0.5), Size.Fixed(1.2)]
v = [Size.Fixed(0.5), Size.Fixed(1.2)]
divider = Divider(fig1, (0, 0, 1, 1), h, v, aspect=False)
axs1 = fig1.add_axes(divider.get_position(),
                     axes_locator=divider.new_locator(nx=1, ny=1))

time = np.linspace(0, 25 - 25/400, 400)
# actually plot the data
for x, pH in enumerate(pH_range):
    if x == 0: # pH 4 trace
        this_max = np.max(raw_data[x, :])
        plt.semilogy(time, raw_data[x, :] / this_max,
                     color=colormap(x/8),
                     label=pH,
                     zorder=10 + (-1*x))
    elif x == 6: # pH 7 trace
        this_max = np.max(raw_data[x, :])
        plt.semilogy(time, raw_data[x, :] / this_max,
                     color=colormap(5/8),
                     label=pH,
                     zorder=10 + (-1*x))
# plotting cleanup
axs1.spines['top'].set_visible(False)
axs1.spines['right'].set_visible(False)
axs1.legend(fontsize=8, ncol=2, frameon=False)
axs1.set_xlabel('Time (ns)')
axs1.set_ylabel('Norm. Photon Count')

fig1.savefig('2021-11-18_PBS_35C_pH7p5_pH4.pdf', bbox_inches='tight',
             transparent=True)

### SAME THING BUT LINEAR Y AXIS
##fig2, axs2 = plt.subplots(1, 1, figsize=(4,4), dpi=300)
### actually plot the data
##for x, pH in enumerate(pH_range):
##    this_max = np.max(raw_data[x, :])
##    plt.plot(time, raw_data[x, :] / this_max, color=colormap(x/8), label=pH,
##             zorder=10 + (-1*x))
### plotting cleanup
##axs2.spines['top'].set_visible(False)
##axs2.spines['right'].set_visible(False)
##axs2.legend(fontsize=10, ncol=2, frameon=False)
##axs2.set_xlabel('Time (ns)')
##axs2.set_ylabel('Norm. Photon Count')
##
##fig2.savefig('2021-11-18_PBS_35C_decays_linear.png', bbox_inches='tight')

plt.show()
