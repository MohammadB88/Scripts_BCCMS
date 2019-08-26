# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:54:29 2019

@author: mb1988
"""

# Implementation of photon-electron interaction up to the first perturbation
# In the later versions, I will put this script in a module/function

import sisl
import math
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
area_phproj = 8*(10**(-18)) # area on which the light is shined.
V = distance_phsource * area_phproj # the volume in which light propagates from source to the monolayer surface under the projected light.
h_bar = constants.value('Planck constant in eV s')/(2*pi_value) # the Planck constant devided by 2*pi
E_gap = 1.7 # energy gap calculated for 2H-phase of MoS2 monolayer (from Transmission or DOS outputs)
omega = E_gap / h_bar# the frequency of the light (for simplicity let's consider it equal to E_gap = = h_bar * omega)

print(h_bar)
#print(constants.value('Planck constant over 2 pi in eV s'))

# it enforces the energy boundary of the perturbation
#if energy_level = (-1)* h_bar * omega: which energy_level ????????????????
#    delta_energy = 1
delta_energy = 1 # for simplicity assume to be 1 meaning only energy_level equals (-1)* h_bar * omega

# the photon flux defined as the number of photonsper unit time per unit area
I_omega = (N * c)/(V * math.sqrt(mu_r * epsilon_r))

print(I_omega)

# it enforces to consider only Hamiltonian elements which are representing "occupied Orbital m and unoccupied Orbital l"
#if m == occupied and l == unoccupied:
#    AdaggerA = 1

AdaggerA = 1 # to just see if the code works

# read in and assign the device Hamiltonian using sisl (let's just try it out)
H_0 = sisl.Hamiltonian.read('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')

print([200 -- 300])
print(H_0[3,4])
print(H_0.geometry.a2o(0, True))
print(H_0.geometry.a2o(227, True))
print(9*227)
print(H_0.geometry.a2o(2250, True))

# just to see if everything is working
z_m = 10
z_l = 5

# perturbation Hamiltonian
H_perturbation = ((2 * pi_value * e_charge)/(h_bar)) * (z_m - z_l) * (((h_bar)/(2 * omega * epsilon * V))**(1/2)) * (N * delta_energy) * (H_0[3,4]) * (AdaggerA)

print(H_perturbation)
