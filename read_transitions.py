# Read the bands file to obtain the transistions (energy differences) at kpoints \Gamma and K
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import linecache


PATH = '/mnt/local/mb1988_data/mos2/strain-on-mos2/supercell88/1mov/D1/strain0.0/elec-prop/MoS2.bands'

# Read the Fermi energy
F_eng = linecache.getline(PATH,1).split()[0]
linecache.clearcache()
print(F_eng)

# Let's first print the 30 energies around the Fermi level at two kpoints
print('Energies at Gamma point: \n')
gamma_eng_empty = []
gamma_eng_full = []
for l in range(62,63):
    line = linecache.getline(PATH,l).rstrip().split()
    linecache.clearcache()
    for eng in line[:7]:
        if eng < F_eng:
            gamma_eng_full.append(float(eng))
        else:
            gamma_eng_empty.append(float(eng))
    print(line)
    
#print(gamma_eng_full)
#print(gamma_eng_empty)

k_eng_empty = []
k_eng_full = []
print('Energies at K point: \n')
for l in range(377*26+62,377*26+63):
    line = linecache.getline(PATH,l).rstrip().split()
    linecache.clearcache()
    for eng in line[:7]:
        if eng < F_eng:
            k_eng_full.append(float(eng))
        else:
            k_eng_empty.append(float(eng))
    print(line)

#print(k_eng_full)
#print(k_eng_empty)

print('Transistions at Gamma point: \n')
for idx_empty in range(len(k_eng_empty)):
    for idx_full in range(len(k_eng_full)):
        gamma_delta_eng = abs(gamma_eng_empty[idx_empty]-gamma_eng_full[idx_full])
        if gamma_delta_eng > 0.70 and gamma_delta_eng < 1.00:
            print(gamma_eng_empty[idx_empty], gamma_eng_full[idx_full], gamma_delta_eng)
            
print('Transistions at K point: \n')
for idx_empty in range(len(k_eng_empty)):
    for idx_full in range(len(k_eng_full)):
        k_delta_eng = abs(k_eng_empty[idx_empty]-k_eng_full[idx_full])
        if k_delta_eng > 0.70 and k_delta_eng < 1.00:
            print(k_eng_empty[idx_empty], k_eng_full[idx_full], k_delta_eng)
