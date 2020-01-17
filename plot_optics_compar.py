# PLotting the figures of the optical spectra of thedefective monolayers

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

fig, ax = plt.subplots(ncols=3, figsize=(20, 10), squeeze=True)
fig.suptitle('Left to Right: the effect of strain in directions X, Y, Shear T1', fontsize=20, fontweight='bold')

#ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       #title='About as simple as it gets, folks')
#ax.grid()


def peak_position(a):
    peaks, properties = find_peaks(a, height=0.1)
    return peaks

def plt_numpy(sub_idx, fn1, fn2, fn3):
    ### Reading in three filenames and plot them on a single figure
    ###
    
    # plotting the spectra
    x1, y1 = np.loadtxt(fn1, skiprows=8, unpack=True)
    ax[sub_idx].plot(x1, y1, label='0.0%', color='black')
    x2, y2 = np.loadtxt(fn2, skiprows=8, unpack=True)
    ax[sub_idx].plot(x2, y2+0.6, label='+2.0%', color='red')
    x3, y3 = np.loadtxt(fn3, skiprows=8, unpack=True)
    ax[sub_idx].plot(x3, y3+1.0, label='-2.0%', color='blue')
    
    # highligh the peak positions and plot the vertical lines at the peack's energies
    print(peak_position(y1[:850]))
    #print(peak_position(y2[:850]))
    #print(peak_position(y3[:850]))
    
    for idx in peak_position(y1[:850]):
        ax[sub_idx].axvline(x=x1[idx], ymin=0, ymax=2, color='black', linestyle='--', lw=2)
    #for idx in peak_position(y2[:850]):
        #ax[sub_idx].axvline(x=x2[idx], ymin=0, ymax=2, color='red', linestyle='--', lw=2)
    #for idx in peak_position(y3[:850]):
        #ax[sub_idx].axvline(x=x3[idx], ymin=0, ymax=2, color='blue', linestyle='--', lw=2)
        
    # the appearance of the plots
    ax[sub_idx].set_xlabel('Energy (eV)', fontsize = 20)
    ax[sub_idx].set_xlim(0,2.1)
    ax[sub_idx].set_ylim(0,2.1)
    ax[sub_idx].tick_params(axis='x', which='major', labelsize=22)
    ax[0].set_ylabel('Intensity', fontsize = 20)
    ax[0].tick_params(axis='y', which='major', labelsize=22)
    ax[1].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    ax[2].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    
    return ax[sub_idx]

axes_0 = plt_numpy(0, 'strain0.0/WSe2-z-unpolarized.EPSIMG', 'D1/strain2.0/WSe2-z-unpolarized.EPSIMG', 'D1/strain-2.0/WSe2-z-unpolarized.EPSIMG')
#print(axes_0)
axes_1 = plt_numpy(1, 'strain0.0/WSe2-z-unpolarized.EPSIMG', 'D2/strain2.0/WSe2-z-unpolarized.EPSIMG', 'D2/strain-2.0/WSe2-z-unpolarized.EPSIMG')

axes_2 = plt_numpy(2, 'strain0.0/WSe2-z-unpolarized.EPSIMG', 'D4/strain2.0/WSe2-z-unpolarized.EPSIMG', 'D4/strain-2.0/WSe2-z-unpolarized.EPSIMG')
#print(handles,handles)

#handles, labels = axes_0.get_legend_handles_labels()

#handles_labels = []
#handles_labels.append(axes_0.get_legend_handles_labels())
#handles_labels.append(axes_1.get_legend_handles_labels())
#handles_labels.append(axes_2.get_legend_handles_labels())

#handles, labels = [sum(lol, []) for lol in zip(*handles_labels)]
#print(handles,labels)

#fig.legend(handles,labels, loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
axes_0.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=3, fancybox=True, shadow=True, fontsize = 19)
axes_1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=3, fancybox=True, shadow=True, fontsize = 19)
axes_2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=3, fancybox=True, shadow=True, fontsize = 19)
plt.subplots_adjust(bottom=0.0, top=1.3, right=0.98, wspace=0.00, hspace=0.0)
plt.tight_layout()
plt.savefig('./compar_strains.png')
plt.show()
