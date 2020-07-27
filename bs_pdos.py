# Plotting the band structures the MoS2 monolayers along with the corresponding orbitals (fatbands)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import linecache


#fig, ax = plt.subplots(ncols=3, figsize=(20, 10), squeeze=True)

# read bs and PDOS from Mo orbitals
bs_Mo = './Mo_4d.dat'
bs_data_Mo = np.loadtxt(bs_Mo, skiprows=3, unpack=True)
print(len(bs_data_Mo), len(bs_data_Mo[0]))
print(type(bs_data_Mo[2]))
print(bs_data_Mo[2].size)
#print(bs_data_Mo[0,0:90])
#print(bs_data_Mo[1,0:90])

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


fig, a = plt.subplots()
a.set_xlim(0,1.65)
a.set_ylim(-8,-0.2)

#inp_orbs = [(bs_data_Mo, 'black'), (bs_data_dx2y2, 'blue'), bs_data_dxy, bs_data_dxz, bs_data_dyz, bs_data_dz2, bs_data_S, bs_data_px, bs_data_py, bs_data_pz]    
inp_orbs = [bs_data_dx2y2, bs_data_dxy, bs_data_dxz, bs_data_dyz, bs_data_dz2,
            bs_data_px, bs_data_py, bs_data_pz]

# ************************************************************************
# BS as base and three sets of plots to include Mo orbitals in the left plot 
# plot bs as base for the plot
for i in range(15):
    step = i*90
    x = bs_data_Mo[0,0+step:90+step]
    y = bs_data_Mo[1,0+step:90+step]
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    #lwidths = bs_data_Mo[2,0+step:90+step]
    lc = LineCollection(segments, linewidth=2, color='black')
    #fig, a = plt.subplots()
    a.add_collection(lc)
    a.autoscale_view(True,True,True)

# *********************  dxy+dx2y2   *********************
for i in range(15):
    step = i*90
    x = inp_orbs[0][0,0+step:90+step]
    y = inp_orbs[0][1,0+step:90+step] 
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[0][2,0+step:90+step] + inp_orbs[1][2,0+step:90+step] 
    lc = LineCollection(segments, linewidth=20*lwidths, color='red')
    a.add_collection(lc)
    a.autoscale_view(True,True,True)
# *********************  dxz+dyz   *********************
for i in range(15):
    step = i*90
    x = inp_orbs[2][0,0+step:90+step]
    y = inp_orbs[2][1,0+step:90+step] 
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[2][2,0+step:90+step] + inp_orbs[3][2,0+step:90+step] 
    lc = LineCollection(segments, linewidth=20*lwidths, color='green')
    a.add_collection(lc)
    a.autoscale_view(True,True,True)
# *********************   dz2    *********************
for i in range(15):
    step = i*90
    x = inp_orbs[4][0,0+step:90+step]
    y = inp_orbs[4][1,0+step:90+step]
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[4][2,0+step:90+step]
    lc = LineCollection(segments, linewidth=20*lwidths, color='blue')
    a.add_collection(lc)
    a.autoscale_view(True,True,True)

# ************************************************************************
# BS as base and two sets of plots to include S orbitals in the right plot 
# plot bs as base for the plot
for i in range(15):
    step = i*90
    x = bs_data_S[0,0+step:90+step]
    y = bs_data_S[1,0+step:90+step]
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    #lwidths = bs_data_S[2,0+step:90+step]
    lc = LineCollection(segments, linewidth=2, color='black')
    #fig, a = plt.subplots()
    a.add_collection(lc)
    a.autoscale_view(True,True,True)

# *********************   px+py   *********************
for i in range(15):
    step = i*90
    x = inp_orbs[0][0,0+step:90+step]
    y = inp_orbs[0][1,0+step:90+step] 
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[0][2,0+step:90+step] + inp_orbs[1][2,0+step:90+step] 
    lc = LineCollection(segments, linewidth=20*lwidths, color='red')
    a.add_collection(lc)
    a.autoscale_view(True,True,True)
# *********************   pz   *********************
for i in range(15):
    step = i*90
    x = inp_orbs[2][0,0+step:90+step]
    y = inp_orbs[2][1,0+step:90+step] 
    #print(x)
    #print(y)
    lwidths=1+x[:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    #print(segments, segments.shape)
    lwidths = inp_orbs[2][2,0+step:90+step] + inp_orbs[3][2,0+step:90+step] 
    lc = LineCollection(segments, linewidth=20*lwidths, color='green')
    a.add_collection(lc)
    a.autoscale_view(True,True,True)



plt.show()
plt.savefig('./fatbands_test.png')
inp = input('press 1: ')
if inp == '1':
    print('finished')
    
#plt.plot(symmetryline, band_engs, color='black')
#plt.fill_between(0.048132, pdos_engs+2, pdos_engs-2, color='red', alpha=0.5)
#plt.savefig('./compar_strains_1mov.png')
#fig.show()

