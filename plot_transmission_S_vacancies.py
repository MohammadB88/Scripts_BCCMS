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
lst_vacancies = ['pristine', 'int_LR_v_s', 'int_LR_v_2stop', 'int_LR_v_2spar']
lst_labels = ['Pristine', '$V_{S}$', '$V_{2S-top}$', '$V_{2S-par}$']


# parameters for the plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 18), sharex=True, squeeze=True)

## add a big axis, hide frame
#fig.add_subplot(111, frameon=False)
## hide tick and tick label of the big axis
#plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#plt.xlabel('Energy (eV)', fontsize=20, labelpad=12.0)
#plt.ylabel('Projected Local DOS (Left Interface)', fontsize=20, labelpad=15.0)

# set the title for left and right plots
ax1.set_title("V = 0.00", fontsize=20, loc='center', pad=-26.0)
ax2.set_title("V = 1.40", fontsize=20, loc='center', pad=-26.0)

# create the box for the inset plots
ax1_insetl = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(0.4, 0.45), bbox_transform=ax1.transAxes, loc='lower left')
#ax1_insetr = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(0.70, 0.55), bbox_transform=ax1.transAxes, loc='lower left')

ax2_insetl = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(0.4, -0.85), bbox_transform=ax1.transAxes, loc='lower left')
#ax2_insetr = zoomed_inset_axes(ax1, 2, bbox_to_anchor=(1.70, 0.55), bbox_transform=ax1.transAxes, loc='lower left')

# Set labels, tick labels, parameters for all plots.
for ax in (ax1,ax2):
    ax.set_xlim(-2.0,2.0)
    ax.set_xlabel('Energy (eV)', fontsize=20)
    ax.set_ylim(0, 3.5)
    ax.set_ylabel('Transmission', fontsize=20) #($10^{-10}$)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)
    #ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)

# **********************  AX1 **************************
# ******************************************************
# parameters for the inset box and plot at the LEFT side
x1, x2, y1, y2 = -1.0, -0.5, -0.005, 0.7
ax1_insetl.set_xlim(x1, x2)
ax1_insetl.set_ylim(y1, y2)
ax1_insetl.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1_insetl.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
ax1_insetl.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1_insetl.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
ax1_insetl.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1_insetl.yaxis.offsetText.set_fontsize(20)
ax1_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)
ax1_insetl.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax1_insetl.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax1_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

## parameters for the inset box and plot at the RIGHT side
#x1, x2, y1, y2 = 1.0, 1.4, -0.005, 0.5
#ax1_insetr.set_xlim(x1, x2)
#ax1_insetr.set_ylim(y1, y2)
#ax1_insetr.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax1_insetr.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
#ax1_insetr.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
#ax1_insetr.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
#ax1_insetr.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#ax1_insetr.yaxis.offsetText.set_fontsize(20)
#ax1_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
#ax1_insetr.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
##ax1_insetr.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
##ax1_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

plt.xticks(visible=False)
plt.yticks(visible=False)

# **********************  AX2 **************************
# ******************************************************
# parameters for the inset box and plot at the LEFT side
x1, x2, y1, y2 = -1.5, -0.8, -0.005, 0.7
ax2_insetl.set_xlim(x1, x2)
ax2_insetl.set_ylim(y1, y2)
ax2_insetl.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2_insetl.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
ax2_insetl.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2_insetl.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
ax2_insetl.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2_insetl.yaxis.offsetText.set_fontsize(20)
ax2_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', bottom=True, top=True)
ax2_insetl.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax2_insetl.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
#ax2_insetl.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

## parameters for the inset box and plot at the RIGHT side
#x1, x2, y1, y2 = 1.4, 1.9, -0.005, 0.5
#ax2_insetr.set_xlim(x1, x2)
#ax2_insetr.set_ylim(y1, y2)
#ax2_insetr.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
#ax2_insetr.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
#ax2_insetr.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
#ax2_insetr.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
#ax2_insetr.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#ax2_insetr.yaxis.offsetText.set_fontsize(20)
#ax2_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
#ax2_insetr.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
##ax2_insetr.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
##ax2_insetr.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

plt.xticks(visible=False)
plt.yticks(visible=False)

# read in and plot the pdos for TS_0.00
for vacancy,label in zip(lst_vacancies, lst_labels):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0])))
    #f_l = np.loadtxt(glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]))[0], unpack=True)
    #print('PRINT', '{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]))
    f_l = np.loadtxt('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[0]), unpack=True)
    ax1.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # plot the inset image
    ax1_insetl.plot(f_l[0], f_l[1], label='{}'.format(label))
    #ax1_insetr.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax1, ax1_insetl, loc1=2, loc2=4, fc="none", ec="0.5")
    #mark_inset(ax1, ax1_insetr, loc1=2, loc2=4, fc="none", ec="0.5")


# read in and plot the pdos for TS_1.40
for vacancy,label in zip(lst_vacancies, lst_labels):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1])))
    #f_l = np.loadtxt(glob.glob('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]))[0], unpack=True)
    #print('PRINT', '{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]))
    f_l = np.loadtxt('{}/{}MoS2_dev.TBT.AVTRANS_Left-Right'.format(vacancy,direc[1]), unpack=True)
    ax2.plot(f_l[0], f_l[1], label='{}'.format(label))

    # plot the inset image
    ax2_insetl.plot(f_l[0], f_l[1], label='{}'.format(label))
    #ax2_insetr.plot(f_l[0], f_l[1], label='{}'.format(label))
    
    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax2, ax2_insetl, loc1=2, loc2=4, fc="none", ec="0.5")
    #mark_inset(ax2, ax2_insetr, loc1=2, loc2=4, fc="none", ec="0.5")
    
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in (ax1,ax2):
    ax.label_outer()


handles, labels = ax1.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.00), ncol=4, fancybox=True, shadow=True, fontsize=20)

#fig.text(0.5, 0.92, 'Projected DOS onto the layers at the left interface', fontsize=20, horizontalalignment='center', verticalalignment='top')

plt.tight_layout(rect=(-0.01,0,1.00,0.94))
#plt.subplots_adjust(left=0.01, right=0.02, wspace=0.00, hspace=0.0, bottom=0.0,top=0.35)

plt.savefig('transmission_S_vacancies.png')

plt.show()
