# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:54:29 2019

@author: mb1988
"""

# Implementation of photon-electron interaction up to the first perturbation
# In the later versions, I will put this script in a module/function

import sisl
import math
import numpy as np
from scipy import constants

# Assign the constant variables
c = constants.value('speed of light in vacuum')
N = 10 # number of photons (it will be asked from the user in the later versions)
e_charge = constants.value('elementary charge') # elementary charge
pi_value = constants.pi # value of pi constant
mu_r = 2.0*(10**(-5)) # relative magnetic susceptibilit (PRL110,247201 (2013))
epsilon_r = 15.2 # relative electric constant (experimental from PRB91, 125304 (2015) and PRB3, 4286 (1971))
epsilon = 1 # dielectric constant (for simplicity let's assume it is in vacuum)

distance_phsource = 5*(10**(-9)) # distance from light source to the monolayer surface 
# reading the coordinates of the atoms at the device's boundary
geom_fh = sisl.io.siesta.xvSileSiesta('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.XV')
geom = geom_fh.read_geometry(species_Z=True)
x_beg = geom[183,0]
z_beg = geom[180,2]
x_end = geom[182,0]
z_end = geom[377,2]
print('These are the device boundaries: \n x_beg:{}, x_end:{}, z_beg:{}, z_end:{}'.format(x_beg,z_beg,x_end,z_end)) # these coordinates are in Angstrom
area_phproj = ((x_end-x_beg)*(z_end-z_beg))*(10**(-20)) # area on which the light is shined.
V = distance_phsource * area_phproj # the volume in which light propagates from source to the monolayer surface under the projected light.

h_bar = constants.value('Planck constant in eV s')/(2*pi_value) # the Planck constant devided by 2*pi
E_gap = 1.7 # energy gap calculated for 2H-phase of MoS2 monolayer (from Transmission or DOS outputs)
omega = E_gap / h_bar# the frequency of the light (for simplicity let's consider it equal to E_gap = = h_bar * omega)

print(h_bar)
#print(constants.value('Planck constant over 2 pi in eV s'))

# read in and assign the device Hamiltonian using sisl (let's just try it out)
H_0 = sisl.Hamiltonian.read('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')

# it enforces the energy boundary of the perturbation
#if energy_level = (-1)* h_bar * omega: which energy_level ????????????????
#    delta_energy = 1
delta_energy = 1 # for simplicity assume to be 1 meaning only energy_level equals (-1)* h_bar * omega

# the photon flux defined as the number of photonsper unit time per unit area
I_omega = (N * c)/(V * math.sqrt(mu_r * epsilon_r))

print(I_omega)

orbital_coordinate = np.zeros((9*198,2))
new_i = 0
for ai in range(0,198): # ai atom index
   for oi in range(0,9): # oi orbital index
       orbital_coordinate[new_i,0] = H_0.geometry.a2o(0, True)[oi]
       orbital_coordinate[new_i,1] = geom[180+ai,2]
       new_i = new_i + 1 # I need this so that these loops iterates and write all the atoms and their coordinate correctly.
#
#print(orbital_coordinate)
#print(geom[180+197,2])

atom_distance = np.zeros((198,9*198))
for ai in range(0,198): 
   for oi in range(0,9*198):
       atom_distance[ai,oi] = geom[180+ai,2] - orbital_coordinate[oi,1] # z_m - z_l 

# it enforces to consider only Hamiltonian elements which are representing "occupied Orbital m and unoccupied Orbital l"
AdaggerA = 1 # to just see if the code works


# read in and assign the device density matrix using sisl (let's just try it out)
#DM = sisl.physics.DensityMatrix.read('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')
fdf = sisl.get_sile('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')
DM = fdf.read_density_matrix(order=['TSDE'])
print(DM)
print(DM[0,1])
print(np.shape(DM))
M = DM.copy() # M stores the mulliken charges
M._csr._D[:,:] += DM._csr._D[:,-1].reshape(-1,1)
print(M._csr._D[:,:])
print(np.shape(M._csr._D[:,:]))
print((DM._csr._D[:,-1].reshape(-1,1))[10,0])
#M._csr._D[:,:] = np.multiply(M._csr._D[:,:], DM._csr._D[:,-1].reshape(-1,1))

print(np.shape(M))
for i in range(1):
    for j in range(50):
        if M[i,j].any() != 0:
            print(M[i,j])
print(M)

#print(H_0[3,4])
#print(H_0.geometry.a2o(0, True))
#print(H_0.geometry.a2o(197, True))
#print(9*197)
#print(H_0.geometry.a2o(2250, True))

# just to see if everything is working
z_m = 10
z_l = 5

print(H_0)
print(H_0[1,1,1])
print(np.shape(H_0))
# perturbation Hamiltonian
H_perturbation = np.zeros(np.shape(H_0))
H_perturbation = ((2 * pi_value * e_charge)/(h_bar)) * (z_m - z_l) * (((h_bar)/(2 * omega * epsilon * V))**(1/2)) * (N * delta_energy) * (H_0) * (AdaggerA)

#print(H_perturbation)
