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
from scipy import constants, sparse, io
import time


time_begin = time.time()
# assign the constant variables
print('Start assigning the constant variables. \n')
c = constants.value('speed of light in vacuum')
#N = 10 #*(10**(6)) # number of photons (it will be asked from the user in the later versions)
e_charge = constants.value('elementary charge') # elementary charge
pi_value = constants.pi # value of pi constant
mu = constants.mu_0 # the magnetic constant // magnetic susceptibility === m kg s-2 A-2
mu_r = 1.0 # 2.0*(10**(-6)) # relative magnetic susceptibility (It's DIMENSIONLESS!!!) (PRL110,247201 (2013))
epsilon = constants.epsilon_0*((55.26349406*(10**(6)))/constants.epsilon_0) # dielectric constant === A2·s4·kg−1·m−3 => e2.GeV-1.fm-1
epsilon_r = 15.0 # relative electric constant (It's DIMENSIONLESS!!!) outofplane = 6.4 inplane = 15.1  (experimental from PRB91, 125304 (2015) and PRB3, 4286 (1971)) (calculation from Laturia, Akash npj 2D Materials and Applications 2018, DOI:10.1038/s41699-018-0050-x)

print(mu, epsilon)
# reading the coordinates of the atoms at the device's boundary
geom_fh = sisl.io.siesta.xvSileSiesta('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.XV')
geom = geom_fh.read_geometry()
x_beg = geom[183,0]
z_beg = geom[180,2]
x_end = geom[182,0]
z_end = geom[377,2]
#print(geom)
print('These are the device boundaries: \n x_beg:{}, x_end:{}, z_beg:{}, z_end:{}'.format(x_beg,z_beg,x_end,z_end)) # these coordinates are in Angstrom
area_phproj = ((x_end-x_beg)*(z_end-z_beg))*(10**(0)) # area on which the light is shined.

h = constants.value('Planck constant in eV s') # the Planck constant in eV
h_bar = h/(2*pi_value) # the Planck constant devided by 2*pi in eV
print('h_bar is equal to {} eV. '.format(h_bar))

# assign the energy gap of MoS2 monolayer
E_gap = 1.73 # energy gap calculated for 2H-phase of MoS2 monolayer (from DFT bandgap calculation)

# read in the photon energy and power density from user
power_density = float(input('Please enter the power density in (W/m^2): '))
photon_energy = float(input('Please enter the photon energy in (eV): '))

print('These are photon properties: \n')
if photon_energy >= E_gap:
    print('photon energy {} is greater than or equal to the bandgap {}.\n'.format(photon_energy,E_gap))
else:
    print('photon energy {} is smaller than the bandgap {}.\n'.format(photon_energy,E_gap))

omega = photon_energy / h_bar # the frequency of the projected light
lambda_value = (2 * pi_value * c)/omega # the wavelength of the projected light
print('photon frequency (times 2pi), omega = {} s^-1, and photon wavelength, lambda = {} m. \n'.format(omega, lambda_value))

# the photon flux is defined as the number of photons per unit time per unit area
I_omega = power_density / (e_charge * photon_energy)

distance_phsource = 1000*(10**(-3)) # distance from light source to the monolayer surface 
V = distance_phsource * area_phproj # * (400*400*(10**(-12))*3.14) # the volume in which light propagates from source to the monolayer surface under the projected light.
#N = (1.0/c) * I_omega * V * math.sqrt(mu_r * epsilon_r) # number of photons

# bose-einstein distribution
boltzman_constant = constants.value('Boltzmann constant in eV/K')
photon_temperature = 3000 # Prof. Jacky's suggestion
N = 1./(math.exp(photon_energy/(boltzman_constant * photon_temperature))-1)
print('Photon flux is equal to {} N_ph/m^2.s, giving {} number of photons hiting the area of {} Ang^2 over the scattering region. \n'.format(I_omega, N, area_phproj))

# the photon intensity is defined as !!!!!??????
#lambda_value = 740*(10**(-9)) # 500*(10**(-9)) # visible light wavelength (380 - 740 nm) TOO SMALL!!??!?
#Intensity = I_omega * (constants.value('Planck constant') * c)/(lambda_value) # h in Joul => Intensity (light power) in W/m^2
#print(constants.value('Planck constant'), c, constants.value('Planck constant') * c)
#print('Photon intensity is equal to {} (N_ph/m^2.s)*(J/s). '.format(Intensity))

# DENSITY MATRIX
#
# read in and assign the device density matrix using sisl (let's just try it out)
print('\nStart reading in the density matrix. \n')
fdf = sisl.get_sile('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')
DM = fdf.read_density_matrix(order=['TSDE'])
#print(DM)
#print(DM[0,1])
print('DM shape (from .TSDE file): {}'.format(np.shape(DM)))
print('\nStart calculating Mulliken populations from the input DM.')
M = DM.copy() # M stores the mulliken charges
print('The sum of absolute value of overlap terms is equal to {}. '.format(np.abs(M._csr._D[:, -1]).sum()))
print('***** If the above number in complex Hamiltonians as SIESTA\'s is ZERO, something is wrong!!!! *****\n') # Later on I will put a try: ... except ... to control this term.
print(np.abs(M._csr._D[:, -1]).sum())
M._csr._D[:,:] *= DM._csr._D[:,-1].reshape(-1,1)

#for i in range(8*9,8*9+10):
    #for j in range(8*9,8*9+10):
        #print(i, j, M[i,j])
            
#print(M)
tmp_nonzero = 0 
tmp_zero = 0 
threshold = 0.1 # amount of Mulliken charge to consider an orbital occupied.
a_m = [] # occupied states = 1, unoccupied states = 0
adagger_l = [] # occupied states = 0, unoccupied states = 1
for i in range(180*9,180*9+198*9):
    if M[i,i][0][0] >= threshold:
        tmp_nonzero += 1
        #print(M[i,i][0][0])
        a_m.append(1)
        adagger_l.append(0) 
    elif M[i,i][0][0] < threshold:
        a_m.append(0)
        adagger_l.append(1)
        tmp_zero += 1

print(np.count_nonzero(a_m), len(a_m), a_m[0:10])
for k in range(1620,1630):
    print(M[k,k][0][0])
print(np.count_nonzero(adagger_l), len(adagger_l), adagger_l[:10])
#nonzero_mulliken = sparse.csr_matrix.count_nonzero(M)
print('There are {} and {} number of terms with M[i,j][0][0] >= {} and M[i,j][0][0] < {}, respectively, in the first image. \n'.format(tmp_nonzero, tmp_zero, threshold, threshold)) 
#print('\nThe first 10 indecis with nonzero Mulliken charges are {}. '.format((np.where(a_m == 1))))

# it enforces to consider only Hamiltonian elements which are representing "occupied Orbital m and unoccupied Orbital l"
#AdaggerA = 1 # to just see if the code works
beginforloop = time.time()
print('a_m:',np.shape(np.array(a_m).reshape(1,len(a_m))),'adagger_l:',np.shape(np.array(adagger_l).reshape(len(adagger_l),1)))
adaggera_value = np.array(adagger_l).reshape(len(adagger_l),1).dot(np.array(a_m).reshape(1,len(a_m)))
#print(np.shape(adaggera_value))
AdaggerA = np.zeros((198*9,198*9,1))
for l in range(180*9,180*9+198*9):
    for m in range(180*9,180*9+198*9):
        if adaggera_value[l-180*9,m-180*9] == 1: # M[m,l][0][0] < threshold
            AdaggerA[l-180*9,m-180*9,0] = 1  # there is possibility "one" for transition from orbital "m" to orbital "l"
        else:
            AdaggerA[l-180*9,m-180*9,0] = 0  # there is possibility "zero" for transition from orbital "m" to orbital "l"
endforloop = time.time()
print('TIME for the for loop: {}'.format(endforloop-beginforloop))

print('This is a 10x10 sample of the matrix AdaggerA to make sure it is working as expected.')
print(np.array(AdaggerA[0:10,0:10,0].reshape(10,10)))

#AdaggerA_nonzero_elm = np.count_nonzero(np.count_nonzero(AdaggerA, axis=2))
AdaggerA_nonzero_elm = np.count_nonzero(np.count_nonzero(AdaggerA, axis=2))
print('There are {} number of nonzero elements in matrix AdaggerA. '.format(AdaggerA_nonzero_elm))

# HAMILTONIAN  
#
# read in and assign the device Hamiltonian using sisl (let's just try it out)
print('\nStart reading in the unperturbed Hamiltonian (H_0). \n')
H_0 = sisl.Hamiltonian.read('/mnt/local/mb1988_data/mos2/devices/1t-2h-interface/armchair/arm2/device_width6/devgeom_constraint_SZP/scat_pristine/MoS2.fdf')

print(H_0)
print('The Hamiltonian term on the first orbital of the first atom in the device region is equal to {}.'.format(H_0[180*9+0,180*9+0,0]))
print(np.shape(H_0))

print('Selecting part of the unperturbed Hamiltonian that belongs to the device region (H_0_lm).')
#try:
    #H_0_lm = np.loadtxt('H_0_lm.dat')
#except:
    #np.save('H_0_lm.dat.txt', H_0_lm)

H_0_lm = np.zeros((198*9,198*9,1))
for l in range(0,198*9):
    for m in range(0,198*9):
        H_0_lm[l,m] = H_0[180*9+l,180*9+m,0]

#print(H_0_lm)
print('It tests if the script correctly read the H_0 in the device region: {}.'.format(H_0_lm[0,0,0]))
print('H_0_lm shape: {}'.format(np.shape(H_0_lm)))
H_0_lm_nonzero_elm = np.count_nonzero(np.count_nonzero(H_0_lm, axis=2))
print('There are {} nonzero elements in matrix H_0_lm. '.format(H_0_lm_nonzero_elm))

# it enforces the energy boundary of the perturbation
#if energy_level = (-1)* h_bar * omega: which energy_level ????????????????
#    delta_energy = 1
print('\n***** I have to also think about how to include the effect of bosonic operators. *****\n')
delta_energy = 1 # for simplicity assume to be 1 meaning only energy_level equals (-1)* h_bar * omega

# ATOM_DISTANCE (z_m - z_l)
#
orbital_coordinate = np.zeros((9*198,2))
new_i = 0
for ai in range(0,198): # ai atom index
   for oi in range(0,9): # oi orbital index
       orbital_coordinate[new_i,0] = H_0.geometry.a2o(0, True)[oi]
       orbital_coordinate[new_i,1] = geom[180+ai,2]
       new_i = new_i + 1 # I need this so that these loops iterates and write all the atoms and their coordinate correctly.
       
#print('Here is a list of orbital indices and their corresponding coordinate: .\n')
#print(orbital_coordinate)
#print(geom[180+197,2])

atom_distance = np.zeros((198*9,198*9,1))
for aoi1 in range(0,198*9): # aoi1 and aoi2 are atom_orbital indices = atom*orbital
   for aoi2 in range(0,198*9):
       atom_distance[aoi1,aoi2] = (orbital_coordinate[aoi1,1] - orbital_coordinate[aoi2,1])*(10**(0)) # z_m - z_l 

print('atom_distance shape: {} '.format(np.shape(atom_distance)))
atom_distance_nonzero_elm = np.count_nonzero(np.count_nonzero(atom_distance, axis=2))
print('There are {} nonzero elements in matrix atom_distance. '.format(atom_distance_nonzero_elm))

# perturbation Hamiltonian
H_perturbation = np.zeros(((180*9+198*9+180*9),(180*9+198*9+180*9)*9,1), dtype=np.complex64)
print(np.shape(H_perturbation))
H_perturbation_device = np.zeros((198*9,198*9,1), dtype=np.complex64)
print('H_perturbation_device shape before assignment: {}'.format(np.shape(H_perturbation_device)))
#print(H_perturbation_device)
#H_perturbation = ((2 * pi_value * e_charge = 1)/(h_bar)) * (z_m - z_l) * (((h_bar)/(2 * omega * epsilon * V))**(1/2)) * (N * delta_energy) * (H_0) * (AdaggerA)
for l in range(0,198*9):
    for m in range(0,198*9):
        H_perturbation_device[l,m,0] = (((1j)* 1)/(h_bar)) * (((h_bar * math.sqrt(mu_r * epsilon_r) * I_omega)/(2 * N * omega * epsilon * c))**(1/2)) * ((2 * pi_value * N)**(1/2) * delta_energy) * atom_distance[l,m,0] * H_0_lm[l,m,0] * AdaggerA[l,m,0] # atom_distance = (z_m - z_l)
        #print(atom_distance[l,m,0])
        #print(H_0_lm[l,m,0])
        #print(AdaggerA[l,m,0])
        #print(H_perturbation_device[l,m,0])

print('\n***** I may have to devide this perturbation by "2" because <l|H_pertb|m> = <m|H_pertb|l>. *****')

#print(H_perturbation_device)
print('H_perturbation_device shape after assignment: {}'.format(np.shape(H_perturbation_device)))

H_pertb_device_nonzero_elm = np.count_nonzero(np.count_nonzero(H_perturbation_device, axis=2))          
print('There are {} number of non-zero elements in the H_perturbation_device. '.format(H_pertb_device_nonzero_elm))

# This is only for the first image. 
#H_perturbation[180*9:180*9+198*9,180*9:180*9+198*9,0] = H_perturbation_device[:,:,0] #(-1j)*
# I sould generalize it to include the effect of perturbation on other images.
for n_img in range(0,9):
    print((180*9+558*9*n_img)-(180*9+198*9+558*9*n_img))
    H_perturbation[180*9:180*9+198*9,180*9+558*9*n_img:180*9+198*9+558*9*n_img,0] = H_perturbation_device[:,:,0] #(-1j)*

#print(180*9, 180*9, H_perturbation[180*9,180*9,:])

H_pertb_nonzero_elm = np.count_nonzero(np.count_nonzero(H_perturbation, axis=2))          
print('There are {} number of non-zero elements in the H_perturbation. '.format(H_pertb_nonzero_elm))

print(np.shape(H_perturbation))


print(H_perturbation[1800,16200,0])
#print(H_perturbation[:,:,1])


# HAMILTONIAN ASSIGNMENT
#
# write the H_perturbation into a Sparse Matrix and then to a .nc file using sisl
beg_final_assign = time.time()

test_csr = sparse.csr_matrix(H_perturbation[:,:,0])
print(sparse.csr_matrix.count_nonzero(test_csr))
print(test_csr)
# A for loop to assign a range of energies to multiple deltaH files.
#for energies in range(1.6547,1.7,1):#np.linspace(-0.1,2.1,11):
    #print('*********************'+str(energies)+'*********************')
    #enrg = int(10*energies)
    #dH = sisl.get_sile('deltaH_{}.dH.nc'.format(enrg), 'w')
    #final_H_pertb = sisl.Hamiltonian(fdf.read_geometry(), dtype = np.complex128)
    ##print(np.shape(H_perturbation[:,0,0].reshape((5022,1))))
    ##final_H_pertb = sparse.csr_matrix((H_perturbation[:,:,0].reshape(5022*5022*9,1), (H_perturbation[:,0,0].reshape((5022,1)),H_perturbation[0,:,0].reshape((5022*9,1)))), shape=(5022,5022*9))
    #print(np.shape(final_H_pertb))
    #final_H_pertb = final_H_pertb.fromsp(fdf.read_geometry(),test_csr)
    #print(final_H_pertb)
    #dH.write_delta(final_H_pertb, E = energies)#, K=[0,0,0])


dH = sisl.get_sile('deltaH_multiE.dH.nc', 'w')
final_H_pertb = sisl.Hamiltonian(fdf.read_geometry(), dtype = np.complex128)
#print(np.shape(H_perturbation[:,0,0].reshape((5022,1))))
#final_H_pertb = sparse.csr_matrix((H_perturbation[:,:,0].reshape(5022*5022*9,1), (H_perturbation[:,0,0].reshape((5022,1)),H_perturbation[0,:,0].reshape((5022*9,1)))), shape=(5022,5022*9))
print(np.shape(final_H_pertb))
final_H_pertb = final_H_pertb.fromsp(fdf.read_geometry(),test_csr)
print(final_H_pertb)
#for energies in np.linspace(-0.1,2.1,121):
    ##print('*********************'+str(energies)+'*********************')
    #for kpoints in np.linspace(0,0.5,5):
        ##print('*********************'+str(kpoints)+'*********************')
        #dH.write_delta(final_H_pertb, E = energies, k = [kpoints,0,0])
for energies in np.linspace(-2.1,2.1,4001):
    #print('*********************'+str(energies)+'*********************')
    dH.write_delta(final_H_pertb, E = energies)

#for kpoints in np.linspace(0,0.5,21):
    ##print('*********************'+str(kpoints)+'*********************')
    #dH.write_delta(final_H_pertb, k = [kpoints,0,0])

#dH.write_delta(final_H_pertb, k = [0.14285714,0,0])

#klist = []
#for kpoints0 in np.linspace(0,10,1):
    #klist = [kpoints0,0,0]
    #print('*********************'+str(kpoints0)+'*********************')
    #dH.write_delta(final_H_pertb, K = klist)#, K=[0,0,0])
#print(klist)
end_final_assign = time.time()

print('Total time for the final assignment of 11 the perturbation Hamiltonian is {}. '.format(end_final_assign-beg_final_assign))


time_end = time.time()

print('It takes {} for this script to calcualte, assign, and write the H_perturbation to an external .nc file.'.format(time_end - time_begin))
