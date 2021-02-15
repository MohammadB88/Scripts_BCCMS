import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter

import glob # matches filenames and directories knowing some part of the name
import numpy as np

direc = ['TS_0.00/', 'TS_1.40/']

layer = [2]
lst_vacancies = ['pristine', 'int_LR_v_mo', 'int_LR_v_mo3s', 'int_LR_v_mo6s']
lst_labels = ['Perfect Int', '$V_{Mo}$', '$V_{Mo3s}$', '$V_{Mo6s}$']
fermi_shifts = [00.0000000000, -0.009075999999999862, 0.0009380000000001054, -0.024042999999999815]

# parameters for the plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True, squeeze=True)

# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel('$E$ - $E_{F}$ (eV)', fontsize=20, labelpad=12.0)
#plt.ylabel('Projected Local DOS (Left Interface)', fontsize=20, labelpad=15.0)
plt.ylabel('LDOS (arb. units)', fontsize=20, labelpad=15.0)

# set the title for left and right plots
ax1.set_title("V = 0.00", fontsize=19, loc='center', pad=-26.0, fontweight="bold")
ax2.set_title("V = 1.40", fontsize=19, loc='center', pad=-26.0, fontweight="bold")

# Set labels, tick labels, parameters for all plots.
for ax in (ax1,ax2):
    ax.set_xlim(-1.0,1.5)
    #ax.set_xlabel('Energy (eV)', fontsize=20)
    ax.set_ylim(0, 40)
    #ax.set_ylabel('Projected Local DOS (Left Interface)', fontsize=20) #($10^{-10}$)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=20, direction='in', bottom=True, top=True)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=14, direction='in', bottom=True, top=True)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.offsetText.set_fontsize(20)
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22, direction='in', left=True, right=True)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    #ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=22)
    
    # make spines (Frames) thicker
    ax.spines["top"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["right"].set_linewidth(2)

# read in and plot the pdos for TS_0.00
for vacancy, label, f_shift in zip(lst_vacancies, lst_labels, fermi_shifts):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[0],layer[0])))
    f_l = np.loadtxt(glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[0],layer[0]))[0], unpack=True)
    #print(f_shift)
    #ax1.plot(f_l[0] + f_shift, f_l[1], label='{}'.format(label))
    ax1.plot(f_l[0], f_l[1], label='{}'.format(label))

# read in and plot the pdos for TS_1.40
for vacancy, label, f_shift in zip(lst_vacancies, lst_labels, fermi_shifts):
    #print(vacancy, label)
    #print('PRINT',glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[1],layer[0])))
    f_l = np.loadtxt(glob.glob('{}/{}dos_Left_int_{}L_*.dat'.format(vacancy,direc[1],layer[0]))[0], unpack=True)
    #ax2.plot(f_l[0] + f_shift, f_l[1], label='{}'.format(label))
    ax2.plot(f_l[0], f_l[1], label='{}'.format(label))

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in (ax1,ax2):
    ax.label_outer()


handles, labels = ax1.get_legend_handles_labels()

fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=4, fancybox=True, shadow=True, fontsize=20)

#fig.text(0.5, 0.92, 'Projected Local DOS onto the layers at the left interface', fontsize=20, horizontalalignment='center', verticalalignment='top')

plt.tight_layout(rect=(-0.04,-0.05,1.015,0.925))
#plt.subplots_adjust(left=0.01, right=0.02, wspace=0.00, hspace=0.0, bottom=0.0,top=0.35)


plt.savefig('PLDOS_Mo_vacancies.png')
#plt.savefig('PLDOS_Mo_vacancies_shift.png')

#plt.show()
