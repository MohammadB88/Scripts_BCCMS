import os
import sys
import io
import numpy as np
import linecache
#from decimal import Decimal
import math
from pathlib import Path
import sisl

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#import pylab as plt

from scipy.interpolate import griddata

#LIST_DIR = ['arm0_small_interface', 'arm1', 'arm2', 'arm3', 'arm4', 'arm5']
#system_dims = np.empty((70,4))

#ref_geom_name = input('Please enter the name of the file to read the reference data: \n')

#input_geom_name = input('Please enter the name of the file to read the input data: \n')

#output_file_name = input('Please enter the name of the output file (in ".txt" format): \n')

os.mkdir('disp_files')

ref_geom_name = sisl.io.siesta.xvSileSiesta('MoS2.XV')
input_geom_name = sisl.io.siesta.xvSileSiesta('SZP/MoS2.XV')

fh_ref_geom = sisl.io.siesta.fdfSileSiesta('disp_files/ref_geom.fdf', 'w')
fh_ref_geom.write_geometry(ref_geom_name.read_geometry())
fh_input_geom = sisl.io.siesta.fdfSileSiesta('disp_files/inp_geom.fdf', 'w')
fh_input_geom.write_geometry(input_geom_name.read_geometry())

ref_geom = 'disp_files/ref_geom.fdf'
inp_geom = 'disp_files/inp_geom.fdf'

output_file_name = 'disp_files/test_disp.txt'

# Read in the number of atoms.
with open('disp_files/ref_geom.fdf', 'r') as fh_NoAtoms:
    for line in fh_NoAtoms:
        if 'NumberOfAtoms' in line:
            # line_NoAtoms is the line where the Number Of Atoms is written.
            line_NoAtoms = (line.rstrip()).split()
            NoAtoms = int(line_NoAtoms[1])

# read and store the index of the first element in files of the initial and final geometry
my_file = Path("{}".format(output_file_name))
if my_file.exists():
    print(my_file.exists())
    exit()
## Open the reference and the input files.
#fh_ref_geom = open(, 'r')
#fh_input_geom = open(, 'r')

# Read in the index number where the reference geometry starts in the input file.
idx_ref_geom = 0
with open('{}'.format(ref_geom)) as fh_ref:
    for line in fh_ref:
        idx_ref_geom = idx_ref_geom + 1
        if '%block AtomicCoordinatesAndAtomicSpecies' in line:
            break
idx_ref_geom = idx_ref_geom + 1  # Starting index of the reference geometry is stored.

# Read in the index number where the input geometry starts in the input file.
idx_inp_geom = 0
with open('{}'.format(inp_geom)) as fh_input:
    for line in fh_input:
        idx_inp_geom = idx_inp_geom + 1
        if '%block AtomicCoordinatesAndAtomicSpecies' in line:
            break
idx_inp_geom = idx_inp_geom + 1  # Starting index of the input geometry is stored.

print('There are *** {} *** number of atoms in the system.'.format(NoAtoms))
print('Initial geoemtry starts from: {}'.format(idx_ref_geom))
print('Final geoemtry starts from: {}'.format(idx_inp_geom))
    
#fhand_initial_geom.close()
#fhand_final_geom.close()

# I will create an array to store reference atomic positions and the displacement.
# I use reference atoms becasue they act as reference for different systems with different displacements.
system_dims = np.empty((NoAtoms,4))
print(type(system_dims))
for atom_idx in range(0, NoAtoms):
    atom_ref = ((linecache.getline('{}'.format(ref_geom), idx_ref_geom+atom_idx)).rstrip()).split()
    linecache.clearcache()
    #print(atom_ref)
    system_dims[atom_idx, 0] = float(atom_ref[0]) # atomic coordinate in X_direction
    system_dims[atom_idx, 1] = float(atom_ref[1]) # atomic coordinate in Y_direction
    system_dims[atom_idx, 2] = float(atom_ref[2]) # atomic coordinate in Z_direction
    atom_inp = ((linecache.getline('{}'.format(inp_geom), idx_inp_geom+atom_idx)).rstrip()).split()
    linecache.clearcache()
    #print(atom_inp)
    #sys_dis[atom,0] = float(atom_inp[0]) # atomic coordinate in X_direction
    #sys_dis[atom,1] = float(atom_inp[1]) # atomic coordinate in Y_direction
    #sys_dis[atom,2] = float(atom_inp[2]) # atomic coordinate in Z_direction
    # calculate the displacement between the optimized and initial atomic positions
    system_dims[atom_idx,3] = math.sqrt((float(atom_inp[0])-float(atom_ref[0]))**2 +(float(atom_inp[1])-float(atom_ref[1]))**2 +(float(atom_inp[2])-float(atom_ref[2]))**2)
    #print(system_dims[atom_idx,3])
np.savetxt('{}'.format(output_file_name), system_dims);
#print('The atomic displacements for the {} input are calculated and stored in {}.'.format(input_geom_name, output_file_name))

# Start plotting the displacement map.
os.chdir('disp_files')

x, y, z, disp = np.loadtxt('test_disp.txt', unpack=True)

XZ = np.vstack((x,z)).T
 
print(min(x),max(x))
print(min(z),max(z))
#print(disp)


#Lx = 5.50099000
#Lz = 98.45600000

xgrid, zgrid = np.mgrid[0.43:5.51, 0.74:98.46]


grid_z0 = griddata(XZ, disp, (xgrid, zgrid), method='nearest', rescale=True)
grid_z1 = griddata(XZ, disp, (xgrid, zgrid), method='linear', rescale=True)
#grid_z2 = griddata(XZ, disp, (xgrid, zgrid), method='cubic', rescale=True)

#print(np.shape(grid_z0), grid_z0)
#print(np.shape(grid_z1), grid_z1)
#print(np.shape(grid_z2), grid_z2)


plt.figure(figsize=(2,6))

#plt.subplot(111)
plt.imshow(grid_z0.T, extent=(0,1,0,1), aspect='auto', origin='lower')
plt.title('Nearest')
plt.colorbar()
plt.show()
#plt.subplot(112)
plt.figure(figsize=(2,6))
plt.imshow(grid_z1.T, extent=(0,1,0,1), aspect='auto', origin='lower')
plt.title('Linear')
plt.colorbar()
plt.show()
#plt.subplot(113)
#plt.imshow(grid_z2.T, extent=(0,1,0,1), aspect='auto', origin='lower')
#plt.title('Cubic')
##plt.gcf().set_size_inches(6, 6)
#plt.colorbar()
#plt.show()

os.chdir('../')
