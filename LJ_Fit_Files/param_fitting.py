# Written by Brian Day

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
# This section calls on real data for DFT calculations.
# --------------------------------------------------
# Pure molecule / surface file
mol_file = '../Config_files/_molecules/argon.xyz'
mol_name = mol_file.split('/')[1].split('.')[0]
surf_file ='surfaces/graphite.cif'
surf_name = surf_file.split('/')[1].split('.')[0]
xyz_file_dir = mol_name+'_'+surf_name

# Define the bodies here
# N.B. -1 shifts the indexing to start at 0
b1 = np.array([1])-1
all_bodies = np.array([b1])
num_bodies = len(all_bodies)
num_params = 2*num_bodies

# Configuration Files
num_configs = 31
config_files_path = "../Config_files/argon_graphite"

# 'Experimental' (DFT) data
energy_file = 'energy_argon_graphite.txt'
dft_energy = pd.read_csv(energy_file,header=None)
E_DFT_orig = np.array(dft_energy)
E_DFT = E_DFT_orig

# --------------------------------------------------
# ----- Test data w/ noise -------------------------
# Generated by providing a sigma and epsilon value and
# calulating the energy and adding noise. Skip sections
# to run correctly.---------------------------------
# --------------------------------------------------
# Test Data (i.e. Perfect Fit) - 1 Body
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
# k = 0
# E_DFT = np.zeros((num_configs,num_bodies))
# E_DFT_wnoise = np.zeros((num_configs,num_bodies))
# for i in range(num_configs):
#     for j in range(num_bodies):
#             E_DFT[i][j] = alpha[i,j,k]*A_LJ - beta[i,j,k]*B_LJ
#
# noise = np.random.normal(0,float(abs(min(E_DFT))),len(E_DFT))
#
# for i in range(len(E_DFT)):
#     E_DFT_wnoise[i] = E_DFT[i] + 15*noise[i]
#
# # Plot the Test Data
# plt.plot(rvect[:,0],E_DFT,'go',rvect[:,0],E_DFT_wnoise,'rx')
# plt.xlabel('Distance [Angstrom]')
# plt.ylabel('Energy [kcal/mol]')
# plt.title('Lennard-Jones Test Data')
# plt.show()
#
# E_DFT = E_DFT_wnoise

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
def get_coordinates(mol_file, mixed_file):
    mol_atoms_file = ase.io.read(mol_file)
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
    # mol_atoms must be an ase object
    # body list must be an array of arrays listing all atom numbers in each body
    # Even for one body, must be array of arrays to work with this function
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

# Create an array of only the surface atoms which are within the cutoff distance
# from the body(s) of interest for a single configuration.
def get_filtered_list(body_pos,surf_atoms,r_cutoff):
    # Prepare empty lists for each output
    surf_atoms_list_per_body = []
    surf_atoms_pos_per_body = []
    distances_per_body = []
    for i in range(len(body_pos)):
        p_body = body_pos[i]
        # Prepare empty lists to be filled during inner loop, and appended to
        # main output. Cleared in between each body loop.
        surf_atoms_list_filtered = []
        surf_atoms_pos_filtered = []
        distance = []
        for j in range(len(surf_atoms)):
            p_surf_atom = surf_atoms.get_positions()[j]
            r_dist = get_distance_between_coords(p_body,p_surf_atom)
            if r_dist <= r_cutoff:
                surf_atoms_list_filtered.append(j)
                surf_atoms_pos_filtered.append(p_surf_atom)
                distance.append(r_dist)
        surf_atoms_list_per_body.append(surf_atoms_list_filtered)
        surf_atoms_pos_per_body.append(surf_atoms_pos_filtered)
        distances_per_body.append(distance)
    return surf_atoms_list_per_body, surf_atoms_pos_per_body, distances_per_body

# Raise a list of floats to a power
def list_to_power(list_input,power):
    import copy
    list_output = copy.deepcopy(list_input)
    for i in range(len(list_output[:])):
        for j in range(len(list_output[i][:])):
            for k in range(len(list_output[i][j][:])):
                list_output[i][j][k] = list_output[i][j][k]**power
    return list_output

# --------------------------------------------------
# ----- Prepare data / configuration files ---------
# --------------------------------------------------
r_cutoff = 10
surf_atoms_index_allconfigs = []
surf_atoms_pos_allconfigs = []
distances_allconfigs = []
for i in range(num_configs):
    # Clear lists in between each configuration
    surf_atoms_index = []
    surf_atoms_pos = []
    distances = []
    # Define file to be worked on
    mixed_file_name = "argon_graphite"+'_'+str(i+1)+'.xyz'
    mixed_file = config_files_path+'/'+mixed_file_name
    # Break apart file and find relevant atoms
    mol_atoms,surf_atoms = get_coordinates(mol_file,mixed_file)
    mol_bodies_cop, mol_bodies_com = create_bodies(mol_atoms,all_bodies)
    surf_atoms_index, surf_atoms_pos, distances = get_filtered_list(mol_bodies_cop,surf_atoms,r_cutoff)
    # Append to lists
    surf_atoms_index_allconfigs.append(surf_atoms_index)
    surf_atoms_pos_allconfigs.append(surf_atoms_pos)
    distances_allconfigs.append(distances)

rvect = np.zeros(len(distances_allconfigs))
for i in range(len(distances_allconfigs)):
    rvect[i] = np.min(distances_allconfigs[i][:])

alpha = list_to_power(distances_allconfigs, -12)
beta = list_to_power(distances_allconfigs, -6)

# --------------------------------------------------
# ----- Test data w/ noise -------------------------
# Generated by calculating the LJ energy from real config
# files with artifical sigma/epsilon values, plus added
# noise. Full code should be used for fitting. -----
# --------------------------------------------------
# sigma = 3.5
# epsilon = 0.000066
# A_LJ = 4 * epsilon * (sigma**12)
# B_LJ = 4 * epsilon * (sigma**6)
# AB_LJ = np.array([A_LJ,B_LJ]).reshape(2,1)
#
# t1 = np.zeros((num_configs,1))
# t2 = np.zeros((num_configs,1))
# E_DFT_orig = np.zeros((num_configs,1))
# for n in range(num_configs): # num_configs is a predetermined, fixed value.
#     t1[n] = np.sum(alpha[n][:][:])*AB_LJ[0]
#     t2[n] = np.sum(beta[n][:][:])*AB_LJ[1]
#     E_DFT_orig[n] = t1[n]-t2[n]
#
# noise = np.random.normal(0,float(abs(min(E_DFT))),len(E_DFT))
# E_DFT_wnoise = np.zeros((num_configs,1))
# for i in range(len(E_DFT)):
#     E_DFT_wnoise[i] = E_DFT_orig[i] + 10*noise[i]
#
# plt.plot(rvect,E_DFT_orig,'go',rvect,E_DFT_wnoise,'rx')
# E_DFT = E_DFT_wnoise

# --------------------------------------------------
# ----- Parameter Fitting --------------------------
# --------------------------------------------------
A = np.zeros((num_params,num_params))
B = np.zeros((num_params,1))

# Calculate the A and B Matricies for the system of equations to be solved.
# Form of equation: Ax = B --> Solution: x = inv(A)*B

# Notes:
# num_configs is a predetermined, fixed value.
# num_params is a predetermined, fixed value.
# len(alpha[n][b]) is equivalent to num_surf_atoms, but changes with specific
# body and configuration. Note that alpha and beta have exact same shape, so
# alpha will be the default for generic indexing purposes.

for n in range(num_configs):
    if E_DFT[n] >= 0.101:
        continue

    # Primary Body Iteration
    for p in range(num_params):
        b = int(p/2)
        for m in range(len(alpha[n][b])):


            if mod(p,2) == 0: # A type parameter
                bp = int(p/2)
                B_nbm = 2*E_DFT[n]*alpha[n][b][m]
                B[p] += B_nbm

                # Secondary Body Iteration (A type params)
                for pp in range(2):
                    bp = int(pp/2)
                    for mp in range(len(alpha[n][bp])):
                        if pp == p: # A type param with itself
                            A_nbm = 2*alpha[n][b][m]*alpha[n][bp][mp]
                            A[p,pp] += A_nbm
                        else:
                            if mod(pp,2) == 0: # A type param with another A type param
                                A_nbm = 2*alpha[n][b][m]*alpha[n][bp][mp]
                                A[p,pp] += A_nbm
                            elif mod(pp,2) == 1: # A type param with B type param
                                bp = int((pp-1)/2)
                                A_nbm = -2*alpha[n][b][m]*beta[n][bp][mp]
                                A[p,pp] += A_nbm

            elif mod(p,2) == 1: # B type parameter
                b = int((p-1)/2)
                B_nbm = -2*E_DFT[n] * beta[n][bp][m]
                B[p] += B_nbm

                # Secondary Body Iteration (B type params)
                for pp in range(num_params):
                    bp = int((pp-1)/2)
                    for mp in range(len(alpha[n][bp])):
                        if pp == p: # B type param with itself
                            A_nbm = 2*beta[n][b][m]*beta[n][bp][mp]
                            A[p,pp] += A_nbm
                        else:
                            if mod(pp,2) == 0: # B type param with A type param
                                bp = int(pp/2)
                                A_nbm = -2*beta[n][b][m]*alpha[n][bp][mp]
                                A[p,pp] += A_nbm
                            elif mod(pp,2) == 1: # B type param with another B type param
                                bp = int((pp-1)/2)
                                A_nbm = 2*beta[n][b][m]*beta[n][bp][mp]
                                A[p,pp] += A_nbm

# --------------------------------------------------
# ----- Extract Solution ---------------------------
# --------------------------------------------------
# First calculate the mixed sigma(s) and epsilon(s) from the fitted values of
# A and B. Then use mixing rules to extract body specific sigma(s) and epsilon(s)
# after parameter fitting. Multiple mixing rules exist which can vary the
# sigma and epsilon parameters calculated. Note that the mixing rules do NOT
# impact the fitted parameters, A and B.

# Print the Resulting Matricies
print('\n----- Built Matricies -----')
print('A matrix =\n', A)
print('B matrix =\n', B)
print('---------------------------\n')

# Calculate the inverse of the coefficient matrix
# Calculate the solution
A_inv = np.linalg.inv(A)
AB_fit = np.matmul(A_inv,B)
print('\n----- Solutions -----')
# print('AB_LJ = \n', AB_LJ)
print('AB_fit = \n', AB_fit)

# Calculate mixed sigma(s) / epsilon(s)
for i in range(num_bodies):
    sigma_fit = (AB_fit[i]/AB_fit[i+1])**(1/6)
    epsilon_fit = (AB_fit[i+1])**2/(4*AB_fit[i])
    print('sigma_fit = ',sigma_fit)
    print('epsilon_fit =',epsilon_fit)

print('---------------------\n')

# --------------------------------------------------
# ----- Calculate the Lennard-Jones Energy ---------
# --------------------------------------------------
t1 = np.zeros((num_configs,1))
t2 = np.zeros((num_configs,1))
E_LJ = np.zeros((num_configs,1))
for n in range(num_configs): # num_configs is a predetermined, fixed value.
    for b in range(num_bodies): # num_params is a predetermined, fixed value.
        t1[n] = np.sum(alpha[n][b]*AB_fit[0])
        t2[n] = np.sum(beta[n][b]*AB_fit[1])
        E_LJ[n] = t1[n]-t2[n]

residuals = E_DFT-E_LJ
SSE = np.sum(residuals**2)
print('Sum of Squared Error = ',SSE)
print('---------------------\n')

# --------------------------------------------------
# ----- Plot the results ---------------------------
# --------------------------------------------------
# E_DFT = alpha[:,0,0]*A_LJ - beta[:,0,0]*B_LJ
# E_fit = alpha[:,0,0]*AB_fit[0] - beta[:,0,0]*AB_fit[1]
# plt.plot(rvect[:,0,0],E_DFT,'go',rvect[:,0,0],E_DFT_wnoise,'rx',rvect[:,0,0],E_LJ,'bx') # Use with arbitrary surf atoms test data
# plt.plot(rvect,E_DFT_orig,'go',rvect,E_DFT_wnoise,'rx',rvect,E_LJ,'bx') # Use with real config files (test data)
plt.plot(rvect,E_DFT_orig,'go',rvect,E_LJ,'bx') # Use with real config files (real data)
plt.xlabel('Distance [Angstrom]')
plt.ylabel('Energy [eV]')
plt.title('Lennard-Jones Test Data')
plt.gca().legend(("Experimental Argon Data","Lennard-Jones Energy"))
plt.show()
