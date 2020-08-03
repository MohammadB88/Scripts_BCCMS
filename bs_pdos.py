# Plotting the band structures the MoS2 monolayers along with the corresponding orbitals (fatbands)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import linecache


num_kpoins = 360

# number of considered bands for ML: 15, for bulk: 21
num_bands = 21 

# fermi level for ML: -4.3487, for bulk: -3.2951453224165num_bands0
fermi_level = -3.2951453224165150

# ylimit for ML: y_min, y_max = -4, 4 ||| for bulk: y_min, y_max = -4, 4
y_min, y_max = -4, 4

# read bs and PDOS from Mo orbitals
bs_Mo = './Mo_4d.dat'
bs_data_Mo = np.loadtxt(bs_Mo, skiprows=3, unpack=True)
print(len(bs_data_Mo), len(bs_data_Mo[0]))
print(type(bs_data_Mo[2]))
print(bs_data_Mo[2].size)
#print(bs_data_Mo[0,0:num_kpoins])
#print(bs_data_Mo[1,0:num_kpoins])

bs_dx2y2 = './Mo_4dx2-y2.dat'
bs_data_dx2y2 = np.loadtxt(bs_dx2y2, skiprows=3, unpack=True)
bs_dxy = './Mo_4dxy.dat'
bs_data_dxy = np.loadtxt(bs_dxy, skiprows=3, unpack=True)
bs_dxz = './Mo_4dxz.dat'
bs_data_dxz = np.loadtxt(bs_dxz, skiprows=3, unpack=True)
bs_dyz = './Mo_4dyz.dat'
bs_data_dyz = np.loadtxt(bs_dyz, skiprows=3, unpack=True)
bs_dz2 = './Mo_4dz2.dat'
bs_data_dz2 = np.loadtxt(bs_dz2, skiprows=3, unpack=True)

# read bs and PDOS from S orbitals
bs_S = './S_3p.dat'
bs_data_S = np.loadtxt(bs_S, skiprows=3, unpack=True)
bs_px = './S_3px.dat'
bs_data_px = np.loadtxt(bs_px, skiprows=3, unpack=True)
bs_py = './S_3py.dat'
bs_data_py = np.loadtxt(bs_py, skiprows=3, unpack=True)
bs_pz = './S_3pz.dat'
bs_data_pz = np.loadtxt(bs_pz, skiprows=3, unpack=True)


fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(20, 10), squeeze=True)

# Set the properties for the left plot
ax1.set_xlim(0,1.652) #1.666106
ax1.set_xticks([0.0, 0.701543, 1.052315, 1.652]) #1.666106
ax1.set_xticklabels((r'$\Gamma$', 'K', 'M', r'$\Gamma$'), size=38, style='oblique')
ax1.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=34)
ax1.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=22)
ax1.set_ylim(y_min, y_max)
ax1.set_ylabel('Energy (eV)', fontsize = 38)
ax1.tick_params(axis='y', which='major', width=2.00, length=5.0, labelsize=34)
ax1.tick_params(axis='y', which='minor', width=1.00, length=3.5, labelsize=22)
#ax1.set_aspect(0.1)

# Set the properties for the right plot
ax2.set_xlim(0,1.652)#1.666106
ax2.set_xticks([0.0, 0.701543, 1.052315, 1.652]) #1.666106
ax2.set_xticklabels((r'$\Gamma$', 'K', 'M', r'$\Gamma$'), size=38, style='oblique')
ax2.tick_params(axis='x', which='major', width=2.00, length=5.0, labelsize=34)
ax2.tick_params(axis='x', which='minor', width=1.00, length=3.5, labelsize=22)
ax2.set_ylim(y_min, y_max)
#ax2.set_ylabel('Energy (eV)', fontsize = 38)
ax2.tick_params(axis='y', which='major', width=2.00, length=5.0, labelleft=False, labelsize=34)
ax2.tick_params(axis='y', which='minor', width=2.00, length=5.0, labelleft=False, labelsize=22)
#ax2.set_aspect(0.1)

#inp_orbs = [(bs_data_Mo, 'black'), (bs_data_dx2y2, 'blue'), bs_data_dxy, bs_data_dxz, bs_data_dyz, bs_data_dz2, bs_data_S, bs_data_px, bs_data_py, bs_data_pz]    
inp_orbs = [bs_data_dx2y2, bs_data_dxy, bs_data_dxz, bs_data_dyz, bs_data_dz2,
            bs_data_px, bs_data_py, bs_data_pz]

# ************************************************************************
# BS as base and three sets of plots to include Mo orbitals in the left plot 
# *********************  dxy+dx2y2   *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = inp_orbs[0][0,0+step:num_kpoins+step]
    y = inp_orbs[0][1,0+step:num_kpoins+step] + (-1) * fermi_level
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    #print(points.shape)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments.shape)
    lwidths = inp_orbs[0][2,0+step:num_kpoins+step] + inp_orbs[1][2,0+step:num_kpoins+step] 
    #lc = LineCollection(segments, linewidth=25*lwidths, color='red', label=r'$d_{xy} + d_{x^2 + y^2}$')
    if i == 4:
        lc = LineCollection(segments, linewidth=25*lwidths, color='red', label=r'$d_{xy} + d_{x^2 + y^2}$')
    else:
        lc = LineCollection(segments, linewidth=25*lwidths, color='red')
    ax1.add_collection(lc)
    ax1.autoscale_view(True,True,True)
# *********************  dxz+dyz   *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = inp_orbs[2][0,0+step:num_kpoins+step]
    y = inp_orbs[2][1,0+step:num_kpoins+step] + -fermi_level 
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[2][2,0+step:num_kpoins+step] + inp_orbs[3][2,0+step:num_kpoins+step] 
    #lc = LineCollection(segments, linewidth=25*lwidths, color='green', label=r'$d_{xz} + d_{yz}$')
    if i == 6:
        lc = LineCollection(segments, linewidth=25*lwidths, color='green', label=r'$d_{xz} + d_{yz}$')
    else:
        lc = LineCollection(segments, linewidth=25*lwidths, color='green')
    ax1.add_collection(lc)
    ax1.autoscale_view(True,True,True)
# *********************   dz2    *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = inp_orbs[4][0,0+step:num_kpoins+step]
    y = inp_orbs[4][1,0+step:num_kpoins+step] + -fermi_level
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[4][2,0+step:num_kpoins+step]
    #lc = LineCollection(segments, linewidth=25*lwidths, color='blue', label=r'$d_{z^2}$')
    if i == 8:
        lc = LineCollection(segments, linewidth=25*lwidths, color='blue', label=r'$d_{z^2}$')
    else:
        lc = LineCollection(segments, linewidth=25*lwidths, color='blue')
    ax1.add_collection(lc)
    ax1.autoscale_view(True,True,True)
# ********************* plot bs as base for the plot *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = bs_data_Mo[0,0+step:num_kpoins+step]
    y = bs_data_Mo[1,0+step:num_kpoins+step] + -fermi_level
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    #lwidths = bs_data_Mo[2,0+step:num_kpoins+step]
    lc = LineCollection(segments, linewidth=2, color='black')
    ax1.add_collection(lc)
    ax1.autoscale_view(True,True,True)
    
# ************************************************************************
# BS as base and two sets of plots to include S orbitals in the right plot 
# *********************   px+py   *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = inp_orbs[5][0,0+step:num_kpoins+step]
    y = inp_orbs[5][1,0+step:num_kpoins+step] + -fermi_level 
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[5][2,0+step:num_kpoins+step] + inp_orbs[6][2,0+step:num_kpoins+step] 
    if i == 11:
        lc = LineCollection(segments, linewidth=25*lwidths, color='indianred', label=r'$p_x + p_y$')
    else:
        lc = LineCollection(segments, linewidth=25*lwidths, color='indianred')
    ax2.add_collection(lc)
    ax2.autoscale_view(True,True,True)
# *********************   pz   *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = inp_orbs[7][0,0+step:num_kpoins+step]
    y = inp_orbs[7][1,0+step:num_kpoins+step] + -fermi_level 
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[7][2,0+step:num_kpoins+step]
    if i == 2:
        lc = LineCollection(segments, linewidth=25*lwidths, color='steelblue', label=r'$p_z$')
    else:
        lc = LineCollection(segments, linewidth=25*lwidths, color='steelblue')
    ax2.add_collection(lc)
    ax2.autoscale_view(True,True,True)
# ********************* plot bs as base for the plot *********************
for i in range(num_bands):
    step = i*num_kpoins
    x = bs_data_S[0,0+step:num_kpoins+step]
    y = bs_data_S[1,0+step:num_kpoins+step] + -fermi_level
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    #lwidths = bs_data_S[2,0+step:num_kpoins+step]
    lc = LineCollection(segments, linewidth=3, color='black')
    ax2.add_collection(lc)
    ax2.autoscale_view(True,True,True)

#(lc1, lc2, lc3), (r'$d_{xy} + d_{x^2 + y^2}$', r'$d_{xz} + d_{yz}$', r'$d_{z^2}$')
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=3, fancybox=True, shadow=True, fontsize = 32)
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=2, fancybox=True, shadow=True, fontsize = 32)
plt.subplots_adjust(right=1.00, left=0.0, wspace=0.00, hspace=0.0) #bottom=0.0, top=1.3, 
plt.tight_layout()
plt.savefig('./bs_pdos_bulk.png')
plt.show()

