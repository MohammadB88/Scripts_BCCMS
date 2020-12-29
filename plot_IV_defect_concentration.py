import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

# the inset represents the zoom-up of the small portion in the parent axes.
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import numpy as np

# read in the IV inputs
# pristine device
f_pris_w6 = np.loadtxt('device_666/IV_pristine.txt', unpack=True)
f_pris_w8 = np.loadtxt('device_888/IV_with_optim/IV_pristine.txt', unpack=True)
f_pris_w10 = np.loadtxt('device_1110/IV_pristine.txt', unpack=True)
# device with V_Mo
f_w6 = np.loadtxt('device_666/IV_mo_both_int.txt', unpack=True)
f_w8 = np.loadtxt('device_888/IV_with_optim/IV_mo_both_int.txt', unpack=True)
f_w10 = np.loadtxt('device_1110/IV_mo_both_int.txt', unpack=True)

# ******************************************************
# ******************************************************
# parameters for the plot
fig = plt.figure(figsize=[20, 12])
plt.xlim(0.0, 1.5)
plt.xlabel('Voltage', fontsize=20)
plt.ylim(0.0, 4.0*(1e-1)) # *(1e-(10-9)) since Current is in nanoAmper
plt.ylabel('Current (nA)', fontsize=20) #($10^{-10}$)

# Set labels, tick labels, parameters for all plots.
plt.axes().xaxis.set_major_locator(ticker.MultipleLocator(0.1))
plt.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', left=True, right=True)
plt.axes().xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
plt.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', left=True, right=True)
plt.axes().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.axes().yaxis.offsetText.set_fontsize(20)
plt.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', left=True, right=True)
plt.axes().yaxis.set_major_locator(ticker.MultipleLocator(0.5*(1e-9)))
#plt.axes().yaxis.set_minor_locator(ticker.MultipleLocator(0.5*(1e-9)))
#plt.tick_params(axis='y', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', left=True, right=True)

# ******************************************************
# ******************************************************
# create the box for the inset plots
fig_inset = zoomed_inset_axes(fig, 1.5, bbox_to_anchor=(0.03, 0.35), loc='lower left') #, bbox_transform=fig.transAxes

# parameters for the inset box and plot
x1, x2, y1, y2 = 0.8, 1.5, -0.001, 0.5
fig_inset.set_xlim(x1, x2)
fig_inset.set_ylim(y1, y2)
fig_inset.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
fig_inset.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
fig_inset.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
fig_inset.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
fig_inset.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
fig_inset.yaxis.offsetText.set_fontsize(20)
fig_inset.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', left=True, right=True)
fig_inset.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#fig_inset.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#fig_inset.tick_params(axis='y', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', left=True, right=True)


# plot the IV curves
#scale_factor = 1e10
#print(scale_factor)
# ******************************************************
# ******************************************************
# Pristine device 
plt.plot(f_pris_w6[0], f_pris_w6[1]*(1e+9), 'o--', linewidth=3, markersize=12, color='black', label='device width 6 L')
plt.plot(f_pris_w8[0], f_pris_w8[1]*(1e+9), 'o--', linewidth=3, markersize=12, color='blue', label='device width 8 L')
plt.plot(f_pris_w10[0], f_pris_w10[1]*(1e+9), 'o--', linewidth=3, markersize=12, color='red', label='device width 10 L')
# device with V_Mo
plt.plot(f_w6[0], f_w6[1]*(1e+9), 'o-', linewidth=3, markersize=12, color='black', label='device width 6 L')
plt.plot(f_w8[0], f_w8[1]*(1e+9), 'o-', linewidth=3, markersize=12, color='blue', label='device width 8 L')
plt.plot(f_w10[0], f_w10[1]*(1e+9), 'o-', linewidth=3, markersize=12, color='red', label='device width 10 L')

# ******************************************************
# ******************************************************
# plot the inset image
fig_inset.plot(f_pris_w6[0], f_pris_w6[1]*(1e+9), marker='o', linewidth=4, markersize=12, label='{}'.format(label))
# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(fig, fig_inset, loc1=1, loc2=3, fc="none", ec="0.5")

plt.title('different concentration of $V_{Mo}$', fontsize=23)
plt.legend(fontsize=23, loc='upper left')
plt.savefig('IV_defect_concentration.png')

plt.show()
