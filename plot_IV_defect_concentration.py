import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

# the inset represents the zoom-up of the small portion in the parent axes.
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import numpy as np


lst_dirs = ['device_666/IV_pristine.txt', 'device_666/IV_mo_both_int.txt', 'device_888/IV_with_optim/IV_pristine.txt', 'device_888/IV_with_optim/IV_mo_both_int.txt', 'device_1110/IV_pristine.txt', 'device_1110/IV_mo_both_int.txt']
lst_labels = ['width 6L - pristine', 'width 6L - V_Mo', 'width 8L - pristine', 'width 8L - V_Mo', 'width 10L - pristine', 'width 10L - V_Mo']
lst_colors = ['black', 'black', 'blue', 'blue', 'red', 'red']

fig, ax = plt.subplots(figsize=[12,18])#$, squeeze=True)

# ******************************************************
# ******************************************************
# parameters for the plot
ax.set_xlim(0.0, 1.5)
ax.set_xlabel('Voltage', fontsize=20)
ax.set_ylim(0.0, 4.0*(1e-0)) # *(1e-(10-9)) since Current is in nanoAmper
ax.set_ylabel('Current (nA)', fontsize=20) #($10^{-10}$)

# ******************************************************
# ******************************************************
# create the box for the inset plots
ax_inset = zoomed_inset_axes(ax, 8, bbox_to_anchor=(0.1, 0.30), bbox_transform=ax.transAxes, loc='lower left')

# parameters for the inset box and plot
x1, x2, y1, y2 = 1.35, 1.5, -0.001, 0.1
ax_inset.set_xlim(x1, x2)
ax_inset.set_ylim(y1, y2)
ax_inset.xaxis.set_major_locator(ticker.MultipleLocator(0.3))
ax_inset.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
ax_inset.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax_inset.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
ax_inset.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax_inset.yaxis.offsetText.set_fontsize(20)
ax_inset.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)
ax_inset.yaxis.set_major_locator(ticker.MultipleLocator(0.5))#*(1e-9)))
#ax_inset.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax_inset.tick_params(axis='y', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)

plt.xticks(visible=False)
plt.yticks(visible=False)

# read in and plot the IV inputs
#scale_factor = 1e10
#print(scale_factor)
# ******************************************************
# ******************************************************
idx = 0
for dr, lab, cl in zip(lst_dirs, lst_labels, lst_colors):
    # read in the IV inputs
    fh = np.loadtxt(dr, unpack=True)
    print(dr, fh[0], fh[1]*(1e+9))
    # plot the IV
    if idx % 2 == 0:
        ax.plot(fh[0], fh[1]*(1e+9), 'o--', linewidth=4, markersize=12, color=cl, label=lab)
        # plot the inset image
        ax_inset.plot(fh[0], fh[1]*(1e+9), 'o--', linewidth=4, markersize=12, color=cl, label=lab)
    
    else:
        ax.plot(fh[0], fh[1]*(1e+9), 'o-', linewidth=4, markersize=12, color=cl, label=lab)
        # plot the inset image
        ax_inset.plot(fh[0], fh[1]*(1e+9), 'o-', linewidth=4, markersize=12, color=cl, label=lab)
        
    idx += 1
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax, ax_inset, loc1=1, loc2=3, fc="none", ec="0.5")

# Set labels, tick labels, parameters for all plots.
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.3))
ax.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', left=True, right=True)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', left=True, right=True)
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.yaxis.offsetText.set_fontsize(20)
ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))#*(1e-9)))
#ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5*(1e-9)))
#ax.tick_params(axis='y', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)

#plt.title('different concentration of $V_{Mo}$', fontsize=23)
#handles, labels = ax.get_legend_handles_labels()
ax.legend(fontsize=23, loc='upper left')
plt.savefig('IV_defect_concentration.png')

plt.show()
