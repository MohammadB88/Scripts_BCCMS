# Scripts_BCCMS
These scripts are written and used during my time as a PhD student at the BCCMS for data processing.
I will explain the functionality of each script along with their names.

# displacement.py
It calculates the displacement of atoms in a structure against a chosen reference.
In this script, one needs to import "sisl" module along with python's built-in modules.

# pc_implementation.py
With this script, I will be able to include the effect of photon-electron interactions into my transport calculations. Thus, making possible a perturbative approach in studying photodetector devices based on 2D materials.

# calculate_IV.py
This script calculates the I(V) curves for non-perturbed device obtained from two methods, calculating all non-equilibrium Hamiltonians and interpolating the Hamiltonian for all the biases. It also calculates the photocurrent coming from the e-ph interaction. At the end, it plots all of these currents in a single image. 

# collect_fatbands.py
I need this script to produce the dat files containing fatbands and select an energy and its fatband of a specific K-point.
