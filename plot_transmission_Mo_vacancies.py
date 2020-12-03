import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

# the inset represents the zoom-up of the small portion in the parent axes.
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#import glob # matches filenames and directories knowing some part of the name
import numpy as np

direc = ['TS_0.00/', 'TS_1.40/']

#layer = [2]

lst_vacancies = ['pristine', 'int_LR_v_mo', 'int_LR_v_mo3s', 'int_LR_v_mo6s']
lst_labels = ['Pristine', '$V_{Mo}$', '$V_{Mo3s}$', '$V_{Mo6s}$']

# parameters for the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), squeeze=True)

# set the title for left and right plots
ax1.set_title("Bias 0.00 eV", fontsize=20, loc='center', pad=-25.0)
ax2.set_title("Bias 1.40 eV", fontsize=20, loc='center', pad=-25.0)

# create the box for the inset plots
ax1_insetl = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(0.33, 0.55), bbox_transform=ax1.transAxes, loc='lower left')
ax1_insetr = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(0.70, 0.55), bbox_transform=ax1.transAxes, loc='lower left')

ax2_insetl = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(1.30, 0.55), bbox_transform=ax1.transAxes, loc='lower left')
ax2_insetr = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(1.70, 0.55), bbox_transform=ax1.transAxes, loc='lower left')

# Set labels, tick labels, parameters for all plots.
for ax in (ax1,ax2):
    ax.set_xlim(-2.0,2.0)
    ax.set_xlabel('Energy (eV)', fontsize=20)
    ax.set_ylim(0, 3.3)
    ax.set_ylabel('Transmission', fontsize=20) #($10^{-10}$)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
    
# **********************  AX1 **************************
# ******************************************************
# parameters for the inset box and plot at the LEFT side
x1, x2, y1, y2 = -0.8, -0.3, -0.005, 0.5
ax1_insetl.set_xlim(x1, x2)
ax1_insetl.set_ylim(y1, y2)
ax1_insetl.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1_insetl.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax1_insetl.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1_insetl.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax1_insetl.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1_insetl.yaxis.offsetText.set_fontsize(20)
ax1_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
ax1_insetl.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax1_insetl.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax1_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

# parameters for the inset box and plot at the RIGHT side
x1, x2, y1, y2 = 1.0, 1.5, -0.005, 0.5
ax1_insetr.set_xlim(x1, x2)
ax1_insetr.set_ylim(y1, y2)
ax1_insetr.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1_insetr.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax1_insetr.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1_insetr.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax1_insetr.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1_insetr.yaxis.offsetText.set_fontsize(20)
ax1_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
ax1_insetr.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax1_insetr.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax1_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

plt.xticks(visible=False)
plt.yticks(visible=False)

# **********************  AX2 **************************
# ******************************************************
# parameters for the inset box and plot at the LEFT side
x1, x2, y1, y2 = -1.3, -0.6, -0.005, 0.5
ax2_insetl.set_xlim(x1, x2)
ax2_insetl.set_ylim(y1, y2)
ax2_insetl.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2_insetl.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax2_insetl.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2_insetl.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax2_insetl.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2_insetl.yaxis.offsetText.set_fontsize(20)
ax2_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
ax2_insetl.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax2_insetl.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax2_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

# parameters for the inset box and plot at the RIGHT side
x1, x2, y1, y2 = 1.4, 2.0, -0.005, 0.5
ax2_insetr.set_xlim(x1, x2)
ax2_insetr.set_ylim(y1, y2)
ax2_insetr.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2_insetr.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax2_insetr.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2_insetr.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax2_insetr.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2_insetr.yaxis.offsetText.set_fontsize(20)
ax2_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
ax2_insetr.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax2_insetr.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax2_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
plt.xticks(visible=False)
plt.yticks(visible=False)

# *************************************
# read in and plot the pdos for TS_0.00
for vacancy,label in zip(lst_vacancies,lst_labels):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0])))
    #f_l = np.loadtxt(glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]))[0], unpack=True)
    #print('PRINT', '{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]))
    f_l = np.loadtxt('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]), unpack=True)
    ax1.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # plot the inset image
    ax1_insetl.plot(f_l[0], f_l[1], label='{}'.format(label))
    ax1_insetr.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax1, ax1_insetl, loc1=2, loc2=4, fc="none", ec="0.5")
    mark_inset(ax1, ax1_insetr, loc1=2, loc2=4, fc="none", ec="0.5")

# read in and plot the pdos for TS_1.40
for vacancy,label in zip(lst_vacancies,lst_labels):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1])))
    #f_l = np.loadtxt(glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]))[0], unpack=True)
    #print('PRINT', '{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]))
    f_l = np.loadtxt('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]), unpack=True)
    ax2.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # plot the inset image
    ax2_insetl.plot(f_l[0], f_l[1], label='{}'.format(label))
    ax2_insetr.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax2, ax2_insetl, loc1=2, loc2=4, fc="none", ec="0.5")
    mark_inset(ax2, ax2_insetr, loc1=2, loc2=4, fc="none", ec="0.5")

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in (ax1,ax2):
    ax.label_outer()


handles, labels = ax1.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.00), ncol=4, fancybox=True, shadow=True, fontsize=20)

#fig.text(0.5, 0.92, 'Projected DOS onto the layers at the left interface', fontsize=20, horizontalalignment='center', verticalalignment='top')

plt.subplots_adjust(right=1.00, wspace=0.00, hspace=0.0, bottom=0.0, top=0.005)
plt.tight_layout()

plt.savefig('transmission_Mo_vacancies.png')

plt.show()
