import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Run this script after extract_decays.py
# It generates a figure with normalized lifetime decays at 8 pHs in vitro.

# read in the decays parsed by the other script
# axis order is pH, time
raw_data = np.loadtxt('2021-11-18_PBS_35C_decays.csv', dtype='float64')
pH_range = [4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5]

# SEMILOG Y AXIS
fig1, axs1 = plt.subplots(1, 1, figsize=(6,4), dpi=300)
colormap = plt.get_cmap('plasma')
time = np.linspace(0, 25 - 25/400, 400)
# actually plot the data
for x, pH in enumerate(pH_range):
    this_max = np.max(raw_data[x, :])
    plt.semilogy(time, raw_data[x, :] / this_max,
                 color=colormap(x/8),
                 label=pH,
                 zorder=10 + (-1*x))
# plotting cleanup
axs1.spines['top'].set_visible(False)
axs1.spines['right'].set_visible(False)
axs1.legend(fontsize=10, ncol=2, frameon=False)
axs1.set_xlabel('Time (ns)')
axs1.set_ylabel('Norm. Photon Count')

fig1.savefig('2021-11-18_PBS_35C_decays_semilog.png', bbox_inches='tight')

# SAME THING BUT LINEAR Y AXIS
fig2, axs2 = plt.subplots(1, 1, figsize=(6,4), dpi=300)
# actually plot the data
for x, pH in enumerate(pH_range):
    this_max = np.max(raw_data[x, :])
    plt.plot(time, raw_data[x, :] / this_max, color=colormap(x/8), label=pH,
             zorder=10 + (-1*x))
# plotting cleanup
axs2.spines['top'].set_visible(False)
axs2.spines['right'].set_visible(False)
axs2.legend(fontsize=10, ncol=2, frameon=False)
axs2.set_xlabel('Time (ns)')
axs2.set_ylabel('Norm. Photon Count')

fig2.savefig('2021-11-18_PBS_35C_decays_linear.png', bbox_inches='tight')

plt.show()
