# Written by Brian Day
# Last updated: 19 April 2019

# --------------------------------------------------
# ----- Import Python Packages ---------------------
# --------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import ase as ase
from ase import Atoms, io, spacegroup, build, visualize
import pandas as pd
from math import sqrt
from random import seed, random

# --------------------------------------------------
# ----- Manual inputs ------------------------------
# --------------------------------------------------
# 'Experimental' (DFT) data
energy_file = 'placeholder.csv'
dft_energy = pd.read_csv(energy_file)

# Pure molecule / surface file
mol_file = 'molecules/sotolon.xyz'
mol_atoms = ase.io.read(mol_file)
mol_name = mol_file.split('/')[1].split('.')[0]
surf_file ='surfaces/graphite.cif'
surf_name = surf_file.split('/')[1].split('.')[0]
xyz_file_dir = mol_name+'_'+surf_name

# Configuration Files
num_configs = 50
config_files_path = xyz_file_dir

# Define the bodies here
# Sotolon 2-body, a
b1 = np.array([1,2,3,6,9,17])-1
b2 = np.array([4,5,7,8,10,11,12,13,14,15,16])-1
all_bodies = np.array([b1,b2])

# --------------------------------------------------
# ----- User-defined Functions ---------------------
# --------------------------------------------------
# Create atoms object for a set of coordinates with associated atom types
def create_ase_atoms(coords, types):
    atoms_object = Atoms(positions=coords.values, symbols=types.values)
    return atoms_object

# Calculate the remiander
def mod(a,b):
    remainder = a%b
    return remainder

# Obtain coordinates of molecule / surface within a configuration
def get_coordinates(mol_file, surf_file, mixed_file):
    mol_atoms_file = ase.io.read(mol_file)
    surf_atoms_file = ase.io.read(surf_file)
    mixed_atoms_file = ase.io.read(mixed_file)
    mol_atoms = mixed_atoms_file[range(len(mol_atoms_file))]
    surf_atoms = mixed_atoms_file[range(len(mol_atoms_file),len(mixed_atoms_file))]
    return mol_atoms, surf_atoms

# Calculate the center of positions (non-weighted center of mass)
def get_center_of_positions(self, scaled=False):
    num_atoms = len(self)
    all_xyz = self.get_positions()
    avg_xyz = np.sum(all_xyz, axis=0)/num_atoms
    return avg_xyz

# Duplicate surface atoms for calculating atoms within the r_cutoff
# Only needed for periodic calculations
def duplicate_surface(surf_atoms):
    # surf_atoms should be an ase atoms object
    return surf_mult

# Create distinct bodies from a set of atoms
def create_bodies(mol_atoms,body_list):
    # mol_atoms should be an ase objects
    # body list should be an array of arrays listing each atom number in each body
    body_cop = np.zeros((len(body_list),3))
    body_com = np.zeros((len(body_list),3))
    for i in range(len(body_list)):
        body = Atoms()
        for j in range(len(body_list[i])):
            body = body.extend(mol_atoms[body_list[i][j]])
        body_cop[i,:] = get_center_of_positions(body)
        body_com[i,:] = body.get_center_of_mass()
    return body_cop, body_com

# Get the distance between two objects
def get_distance_between_coords(p1, p2):
    squared_dist = 0
    if len(p1) != len(p2):
        print('Invalid input: Postion vectors of different dimesnions')
    else:
        for i in range(len(p1)):
            squared_dist += (p1[i] - p2[i])**2
        dist = sqrt(squared_dist)
    return dist

# --------------------------------------------------
# ----- Prepare data / configuration files ---------
# --------------------------------------------------

for i in range(num_configs)

# --------------------------------------------------
# ----- Test data w/ noise -------------------------
# --------------------------------------------------
# # Test Data (i.e. Perfect Fit) - 1 Body
# num_configs = 101
# num_surf_atoms = 1
# num_bodies = 1
# num_params = 2*num_bodies
#
# sigma = 3.5
# epsilon = 0.066
# A_LJ = 4 * epsilon * (sigma**12)
# B_LJ = 4 * epsilon * (sigma**6)
# AB_LJ = np.array([A_LJ,B_LJ]).reshape(2,1)
#
# rlim = ([3,8])
# rvect = np.ones((num_configs, num_surf_atoms, num_bodies))
# rvect[:,0,0] = np.linspace(rlim[0], rlim[1], num=num_configs, endpoint=True)
# alpha = rvect**(-12)
# beta = rvect**(-6)
#
# E_DFT = np.ones((num_configs))
# E_DFT = alpha[:,0,0]*A_LJ - beta[:,0,0]*B_LJ
# len(E_DFT)
# abs(min(E_DFT))
# noise = np.random.normal(0,abs(min(E_DFT)),len(E_DFT))
# E_DFT_wnoise = E_DFT + noise;
#
# # Plot the Test Data
# plt.plot(rvect[:,0,0],E_DFT[:],'g-o',rvect[:,0,0],E_DFT_wnoise[:],'r-o')
# plt.xlabel('Distance [Angstrom]')
# plt.ylabel('Energy [kcal/mol]')
# plt.title('Lennard-Jones Test Data')
# plt.show()
#
# E_DFT = E_DFT_wnoise

# --------------------------------------------------
# ----- Parameter Fitting --------------------------
# --------------------------------------------------
# num_configs is a fixed value determined prior to any fitting.
# num_params is a fixed value based on the number of bodies.
# num_surf_atoms changes depending on the configuration, and hence so do alpha
# and beta. Add a filtering step which determines the number of surface atoms
# within the cutoff, as well as their position and/or distance between each body.
# Store the surf atoms critical to each config/body as an array of arrays.
A = np.zeros((num_params,num_params))
B = np.zeros((num_params,1))

# Calculate the A and B Matricies for the system of equations to be solved.
# Form of equation: Ax = B --> Solution: x = inv(A)*B
for n in range(num_configs):
    for m in range(num_surf_atoms[n,b]):
    # If switching to an array of arrays to take advantage of r_cutoff for speed
    # purposes, can change to in range (len(specific row/column/etc)).
    # Should add a preprocessing step for the data to handle this.

        for p in range(num_params):
            if mod(p,2) == 0: # A type parameter
                b = int(p/2)
                B_nmb = 2*E_DFT[n]*alpha[n,m,b]
                B[p] += B_nmb

                for pp in range(num_params):
                    # print('n: ', n, '\tm: ', m, '\tp: ', p, '\tpp: ', pp)
                    if pp == p: # A type param with itself
                        bp = int(pp/2)
                        A_nmb = 2*alpha[n,m,b]*alpha[n,m,bp]
                        A[p,pp] += A_nmb
                    else:
                        if mod(pp,2) == 0: # A type param with another A type param
                            bp = int(pp/2)
                            A_nmb = 2*alpha[n,m,b]*alpha[n,m,bp]
                            A[p,pp] += A_nmb
                        elif mod(pp,2) == 1: # A type param with B type param
                            bp = int((pp-1)/2)
                            A_nmb = -2*alpha[n,m,b]*beta[n,m,bp]
                            A[p,pp] += A_nmb

            elif mod(p,2) == 1: # B type parameter
                b = int((p-1)/2)
                B_nmb = -2*E_DFT[n] * beta[n,m,b]
                B[p] += B_nmb

                for pp in range(num_params):
                    # print('n: ', n, '\tm: ', m, '\tp: ', p, '\tpp: ', pp)
                    if pp == p: # B type param with itself
                        bp = int((pp-1)/2)
                        A_nmb = 2*beta[n,m,b]*beta[n,m,bp]
                        A[p,pp] += A_nmb
                    else:
                        if mod(pp,2) == 0: # B type param with A type param
                            bp = int(pp/2)
                            A_nmb = -2*beta[n,m,b]*alpha[n,m,bp]
                            A[p,pp] += A_nmb
                        elif mod(pp,2) == 1: # B type param with another B type param
                            bp = int((pp-1)/2)
                            A_nmb = 2*beta[n,m,b]*beta[n,m,bp]
                            A[p,pp] += A_nmb

# Print the Resulting Matricies
print('\n----- Built Matricies -----')
print('A matrix =\n', A)
print('B matrix =\n', B)
print('---------------------------')

# Calculate the inverse of the coefficient matrix
# Calculate the solution
A_inv = np.linalg.inv(A)
AB_fit = np.matmul(A_inv,B)
print('\n\n----- Solutions -----')
print('AB_LJ = \n', AB_LJ)
print('AB_fit = \n', AB_fit)

# --------------------------------------------------
# ----- Extract Solution ---------------------------
# --------------------------------------------------
# First calculate the mixed sigma(s) and epsilon(s) from the fitted values of
# A and B. Then use mixing rules to extract body specific sigma(s) and epsilon(s)
# after parameter fitting. Multiple mixing rules exist which can vary the
# sigma and epsilon parameters calculated. Note that the mixing rules do NOT
# impact the fitted parameters, A and B.

# Calculate mixed sigma(s) / epsilon(s)
for i in range(num_bodies):
    sigma_fit = (AB_fit[i]/AB_fit[i+1])**(1/6)
    epsilon_fit = (AB_fit[i+1])**2/(4*AB_fit[i])
    print('sigma_fit = ',sigma_fit)
    print('epsilon_fit =',epsilon_fit)

print('---------------------\n')

# --------------------------------------------------
# ----- Plot the results ---------------------------
# --------------------------------------------------
E_DFT = alpha[:,0,0]*A_LJ - beta[:,0,0]*B_LJ
E_fit = alpha[:,0,0]*AB_fit[0] - beta[:,0,0]*AB_fit[1]
plt.plot(rvect[:,0,0],E_DFT[:],'g-o',rvect[:,0,0],E_DFT_wnoise[:],'r-o',rvect[:,0,0],E_fit,'b-o')
plt.xlabel('Distance [Angstrom]')
plt.ylabel('Energy [kcal/mol]')
plt.title('Lennard-Jones Test Data')
plt.gca().legend(('Test Data','Test Data w/ Noise','Fit to Noisy Data'))
plt.show()
