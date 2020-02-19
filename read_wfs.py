# Read the wfs-coefficient file to obtain the major orbitals of atoms surrounding the defect
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import linecache

PATH = '/mnt/local/mb1988_data/mos2/plot-wavefunction/88-supercell/1mov/D1/strain2.0/wavefunction-coefficient/wavefunction-coefficient.txt'
#PATH = '/mnt/local/mb1988_data/mos2/plot-wavefunction/88-supercell/1mov/strain0.0/wavefunction-coefficient/wavefunction-coefficient.txt'

#line = linecache.getline(PATH,12+3377).rstrip().split()
#linecache.clearcache()
#print(line)

# Store the band number of the DLs and band edges
list_bands = []
for i in range(571,571+7):
    list_bands.append(i)

# Read the idx of the above bands in the wfs file
idx_bands = []
with open(PATH, 'r') as fh:
    for i in range(1,70000):
        line = fh.readline().rstrip().split()
        if 'Wavefunction' in line:
            if int(line[2]) >= list_bands[0] and int(line[2]) <= list_bands[-1]:
                #print(i, line[2])
                idx_bands.append(i)
print(idx_bands)

# List of atomic numbers accoring to the geometry
atom_idx = [37,38,45,46,53,54,91,99,100,155,163,164]

# Read the orbitals of atoms around the vacancy involved the above DLs
output = open('required_wfs_strain+20.txt', 'w')
with open(PATH, 'r') as fh:
    line = fh.readlines()
    print(len(idx_bands))
    for i in range(len(idx_bands)):
        output.write('-----------------------------------------------------------------------\n')
        output.write(line[idx_bands[i]-1])
        output.write(line[idx_bands[i]])
        output.write('-----------------------------------------------------------------------\n')
        for j in range(idx_bands[i]+10,idx_bands[i]+3000):
            words = line[j].rstrip().split()
            if '-----------------------------------------------------------------------' in words:
                break
            #print(words)
            #any(words[0] == '{}'.format(num) for num in atom_idx) and 
            if abs(float(words[5])) >= 0.075:
                output.write('{}\n'.format([words[i] for i in [0,1,4,5]]))
                #print(words)q

output.close()
