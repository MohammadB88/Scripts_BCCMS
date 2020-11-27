import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

import glob # matches filenames and directories knowing some part of the name
import numpy as np

direc = ['TS_0.00/', 'TS_1.40/']

layer = [2]
lst_vacancies = ['pristine', 'int_LR_v_s_with_mo', 'int_LR_v_2stop_with_mo', 'int_LR_v_mo_with_s', 'int_LR_v_mo_with_2stop', 'int_LR_v_mo_split']
lst_labels = ['Pristine', '$Mo_{S}$', '$Mo_{2S-top}$', '$S_{Mo}$', '$2S-top_{Mo}$', '$Mo_{split}$']

# parameters for the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), squeeze=True)

# Set labels, tick labels, parameters for all plots.
for ax in (ax1,ax2):
    ax.set_xlim(-2.0,2.0)
    ax.set_xlabel('Energy (eV)', fontsize=20)
    ax.set_ylim(0, 100)
    ax.set_ylabel('Projected DOS', fontsize=20) #($10^{-10}$)

ax1.set_title("Bias 0.00 eV", fontsize=20, loc='center', pad=-25.0)
ax2.set_title("Bias 1.40 eV", fontsize=20, loc='center', pad=-25.0)


ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.yaxis.offsetText.set_fontsize(20)
ax1.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax1.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)

# read in and plot the pdos for TS_0.00
for vacancy,label in zip(lst_vacancies,lst_labels):
    #print(vacancy, label)
    print('PRINT',glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[0],layer[0])))
    f_l = np.loadtxt(glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[0],layer[0]))[0], unpack=True)
    ax1.plot(f_l[0], f_l[1], label='{}'.format(label))

#ax2.set_xlim(-2.0,2.0)
#ax2.set_xlabel('Energy (eV)', fontsize=20)
#ax2.set_ylim(0, 100)
#ax2.set_ylabel('Projected DOS', fontsize=20) #($10^{-10}$)
#ax2.set_title("Bias 1.40 eV", fontsize=20, loc='center', pad=-25.0)

ax2.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax2.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20)
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14)
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.yaxis.offsetText.set_fontsize(20)
ax2.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
#ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax2.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)


# read in and plot the pdos for TS_1.40
for vacancy,label in zip(lst_vacancies,lst_labels):
    #print(vacancy, label)
    print('PRINT',glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[1],layer[0])))
    f_l = np.loadtxt(glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[1],layer[0]))[0], unpack=True)
    ax2.plot(f_l[0], f_l[1], label='{}'.format(label))

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in (ax1,ax2):
    ax.label_outer()


handles, labels = ax1.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.00), ncol=6, fancybox=True, shadow=True, fontsize=20)

fig.text(0.5, 0.92, 'Projected DOS onto the layers at the left interface', fontsize=20, horizontalalignment='center', verticalalignment='top')

plt.subplots_adjust(right=1.00, wspace=0.00, hspace=0.0, bottom=0.0, top=0.005)
plt.tight_layout()


plt.savefig('PDOS_antisites.png')

plt.show()
