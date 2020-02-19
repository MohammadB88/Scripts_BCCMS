# Plotting the figures of the optical spectra of thedefective monolayers

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

fig, ax = plt.subplots(ncols=3, figsize=(20, 10), squeeze=True)
#fig.suptitle('Left to Right: the effect of strain in directions X, Y, Shear T1', fontsize=20, fontweight='bold')

#ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       #title='About as simple as it gets, folks')
#ax.grid()


def peak_position(a):
    peaks, properties = find_peaks(a, height=0.01)
    return peaks

def plt_numpy_only0(sub_idx, fn1, fn2, fn3, fn4):#, fn5, fn6):
    ### Reading in three filenames and plot them on a single figure
    ###
    
    # plotting the spectra
    x1, y1 = np.loadtxt(fn1, skiprows=8, unpack=True)
    ax[sub_idx].plot(x1, y1, label='0.0%', color='black')
    x2, y2 = np.loadtxt(fn2, skiprows=8, unpack=True)
    ax[sub_idx].plot(x2, y2+0.6, label='+2.0%', color='red')
    x3, y3 = np.loadtxt(fn3, skiprows=8, unpack=True)
    ax[sub_idx].plot(x3, y3+1.0, label='-2.0%', color='blue')
    #x4, y4 = np.loadtxt(fn4, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x4, y4, color='black', ls='--')
    #x5, y5 = np.loadtxt(fn5, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x5, y5+0.6, color='red', ls='--')
    #x6, y6 = np.loadtxt(fn6, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x6, y6+1.0, color='blue', ls='--')
    
    # highligh the peak positions and plot the vertical lines at the peack's energies
    output = open('peak_position_1mov.txt', 'w') 
    print('\n Here are the indices and energies for peaks at {} strain.\n'.format(sub_idx))
    print(peak_position(y1[:850]))
    print([x1[idx] for idx in peak_position(y1[:850]) if x1[idx] < 1.0])
    output.write('strain 0.0% \n ')
    output.write('{}\n'.format([x1[idx] for idx in peak_position(y1[:850]) if x1[idx] < 1.0]))
    output.write('-----------------------------------\n')
    print(peak_position(y2[:850]))
    print([x2[idx] for idx in peak_position(y2[:850]) if x1[idx] < 1.25])
    output.write('strain +2.0% \n ')
    output.write('{}\n'.format([x2[idx] for idx in peak_position(y2[:850]) if x1[idx] < 1.25]))
    output.write('-----------------------------------\n')
    print(peak_position(y3[:850]))
    print([x3[idx] for idx in peak_position(y3[:850]) if x1[idx] < 1.25])
    output.write('strain -2.0% \n ')
    output.write('{}\n'.format([x3[idx] for idx in peak_position(y3[:850]) if x1[idx] < 1.25]))
    output.write('-----------------------------------\n')
    
    #for idx in peak_position(y1[:850]):
        #if x1[idx] < 1.0:
            #ax[sub_idx].axvline(x=x1[idx], ymin=0, ymax=2, color='black', linestyle='--', lw=2)
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
    ax[0].tick_params(axis='y', which='major', width=1.00, length=5.0, labelsize=22)
    ax[0].tick_params(axis='y', which='minor', width=0.75, length=2.5, labelsize=10)
    
    # add the labels to the plots
    ax[0].text(0.1, 1.95, 'X-direction', fontsize=22, style='oblique')
    ax[1].text(0.1, 1.95, 'Y-direction', fontsize=22, style='oblique')
    ax[2].text(0.1, 1.95, 'Shear T1', fontsize=22, style='oblique')
    
    # arguments for the inset plot
    axins = ax[0].inset_axes([0.15, 0.7, 0.3, 0.2])
    x4, y4 = np.loadtxt(fn4, skiprows=8, unpack=True)
    axins.plot(x4, y4, color='black', ls='--')
    axins.set_xlim(0,0.7)
    axins.set_ylim(0,0.3)
    axins.xaxis.set_major_locator(ticker.MultipleLocator(0.3))
    axins.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    axins.set_xlabel('Energy (eV)', fontsize = 14)
    axins.set_ylabel('Intensity', fontsize = 14)
    
    ax[sub_idx].yaxis.set_major_locator(ticker.MultipleLocator(0.4))
    ax[sub_idx].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    #ax[0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x}"))
    ax[1].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    #ax[1].tick_params(axis='y', which='minor', labelsize=16, labelleft=False)
    ax[2].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    #ax[2].tick_params(axis='y', which='minor', labelsize=16, labelleft=False)
    
    return ax[sub_idx]

def plt_numpy(sub_idx, fn1, fn2, fn3):#, fn4, fn5, fn6):
    ### Reading in three filenames and plot them on a single figure
    ###
    
    # plotting the spectra
    x1, y1 = np.loadtxt(fn1, skiprows=8, unpack=True)
    ax[sub_idx].plot(x1, y1, label='0.0%', color='black')
    x2, y2 = np.loadtxt(fn2, skiprows=8, unpack=True)
    ax[sub_idx].plot(x2, y2+0.6, label='+2.0%', color='red')
    x3, y3 = np.loadtxt(fn3, skiprows=8, unpack=True)
    ax[sub_idx].plot(x3, y3+1.0, label='-2.0%', color='blue')
    #x4, y4 = np.loadtxt(fn4, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x4, y4, color='black', ls='--')
    #x5, y5 = np.loadtxt(fn5, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x5, y5+0.6, color='red', ls='--')
    #x6, y6 = np.loadtxt(fn6, skiprows=8, unpack=True)
    #ax[sub_idx].plot(x6, y6+1.0, color='blue', ls='--')
    
    # highligh the peak positions and plot the vertical lines at the peack's energies
    print('\n Here are the indices and energies for peaks at {} strain.\n'.format(sub_idx))
    print(peak_position(y1[:850]))
    print([x1[idx] for idx in peak_position(y1[:850]) if x1[idx] < 1.0])
    print(peak_position(y2[:850]))
    print([x2[idx] for idx in peak_position(y2[:850]) if x1[idx] < 1.25])
    print(peak_position(y3[:850]))
    print([x3[idx] for idx in peak_position(y3[:850]) if x1[idx] < 1.25])
    
    #for idx in peak_position(y1[:850]):
        #if x1[idx] < 1.0:
            #ax[sub_idx].axvline(x=x1[idx], ymin=0, ymax=2, color='black', linestyle='--', lw=2)
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
    ax[0].tick_params(axis='y', which='major', width=1.00, length=5.0, labelsize=22)
    ax[0].tick_params(axis='y', which='minor', width=0.75, length=2.5, labelsize=10)
    ax[sub_idx].yaxis.set_major_locator(ticker.MultipleLocator(0.4))
    ax[sub_idx].yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    #ax[0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x}"))
    ax[1].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    #ax[1].tick_params(axis='y', which='minor', labelsize=16, labelleft=False)
    ax[2].tick_params(axis='y', which='major', labelsize=22, labelleft=False)
    #ax[2].tick_params(axis='y', which='minor', labelsize=16, labelleft=False)
    
    return ax[sub_idx]

axes_0 = plt_numpy_only0(0, '1mov-D1/z-unpolarized/MoS2.EPSIMG-strain+00', '1mov-D1/z-unpolarized/MoS2.EPSIMG-strain+20', '1mov-D1/z-unpolarized/MoS2.EPSIMG-strain-20', '1mov-D1/z-polarized/MoS2.EPSIMG-strain+00')#, '1mov-D1/z-polarized/MoS2.EPSIMG-strain+20', '1mov-D1/z-polarized/MoS2.EPSIMG-strain-20')
#print(axes_0)
axes_1 = plt_numpy(1, '1mov-D2/z-unpolarized/MoS2.EPSIMG-strain+00', '1mov-D2/z-unpolarized/MoS2.EPSIMG-strain+20', '1mov-D2/z-unpolarized/MoS2.EPSIMG-strain-20')#, '1mov-D2/z-polarized/MoS2.EPSIMG-strain+00', '1mov-D2/z-polarized/MoS2.EPSIMG-strain+20', '1mov-D2/z-polarized/MoS2.EPSIMG-strain-20')

axes_2 = plt_numpy(2, '1mov-D4/z-unpolarized/MoS2.EPSIMG-strain+00', '1mov-D4/z-unpolarized/MoS2.EPSIMG-strain+20', '1mov-D4/z-unpolarized/MoS2.EPSIMG-strain-20')#, '1mov-D4/z-polarized/MoS2.EPSIMG-strain+00', '1mov-D4/z-polarized/MoS2.EPSIMG-strain+20', '1mov-D4/z-polarized/MoS2.EPSIMG-strain-20')
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
plt.subplots_adjust(right=0.98, wspace=0.00, hspace=0.0) #bottom=0.0, top=1.3, 
plt.tight_layout()
plt.savefig('./compar_strains_1mov.png')
plt.show()
