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
# Pure molecule / surface file
mol_file = 'molecules/sotolon.xyz'
mol_name = mol_file.split('/')[1].split('.')[0]
surf_file ='surfaces/graphite.cif'
surf_name = surf_file.split('/')[1].split('.')[0]
xyz_file_dir = mol_name+'_'+surf_name

# Define the bodies here
# N.B. -1 shifts the indexing to start at 0
# Sotolon 2-body, a
b1 = np.array([1,2,3,6,9,17])-1
b2 = np.array([4,5,7,8,10,11,12,13,14,15,16])-1
all_bodies = np.array([b1,b2])
num_bodies = len(all_bodies)
num_params = 2*num_bodies

# Configuration Files
num_configs = 50
config_files_path = xyz_file_dir

# 'Experimental' (DFT) data
energy_file = 'placeholder.csv'
dft_energy = pd.read_csv(energy_file)

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
r_cutoff = 8
surf_atoms_index_allconfigs = []
surf_atoms_pos_allconfigs = []
distances_allconfigs = []
for i in range(num_configs):
    # Clear lists in between each configuration
    surf_atoms_index = []
    surf_atoms_pos = []
    distances = []
    # Define file to be worked on
    mixed_file_name = config_files_path+'_'+str(i+1)+'.xyz'
    mixed_file = config_files_path+'/'+mixed_file_name
    # Break apart file and find relevant atoms
    mol_atoms,surf_atoms = get_coordinates(mol_file,mixed_file)
    mol_bodies_cop, mol_bodies_com = create_bodies(mol_atoms_test,all_bodies)
    surf_atoms_index, surf_atoms_pos, distances = get_filtered_list(mol_bodies_cop,surf_atoms_test,r_cutoff)
    # Append to lists
    surf_atoms_index_allconfigs.append(surf_atoms_index)
    surf_atoms_pos_allconfigs.append(surf_atoms_pos)
    distances_allconfigs.append(distances)

alpha = list_to_power(distances_allconfigs, -12)
beta = list_to_power(distances_allconfigs, -6)

# --------------------------------------------------
# ----- Parameter Fitting --------------------------
# --------------------------------------------------
A = np.zeros((num_params,num_params))
B = np.zeros((num_params,1))

# Calculate the A and B Matricies for the system of equations to be solved.
# Form of equation: Ax = B --> Solution: x = inv(A)*B
for n in range(num_configs): # num_configs is a predetermined, fixed value.
    for p in range(num_params): # num_params is a predetermined, fixed value.
        b = int(p/2)
        for m in range(len(alpha[n][b])):
        # len(alpha[n][b]) is equivalent to num_surf_atoms, but changes with
        # specific body and configuration. Note that alpha and beta have exact
        # same shape, so alpha will be the default for generic indexing purposes.
            if mod(p,2) == 0: # A type parameter
            # Cycles through all surface atoms for that body / config.
            # No indexing issues.
                bp = int(p/2)
                B_nbm = 2*E_DFT[n]*alpha[n][b][m]
                B[p] += B_nbm

                for pp in range(2):
                    if pp == p: # A type param with itself
                    # Since it's a parameter with itself, len(r[n][b]) will
                    # always = len(r[n][bp]). No indexing issues.
                        bp = int(pp/2)
                        A_nbm = 2*alpha[n][b][m]*alpha[n][bp][m]
                        A[p,pp] += A_nbm
                    else:
                        if mod(pp,2) == 0: # A type param with another A type param
                        # Possible indexing issues from body to body due to
                        # different number of surface atoms within cutoff.
                            bp = int(pp/2)
                            if m >= len(alpha[n][bp]):
                                pass
                            else:
                                A_nbm = 2*alpha[n][b][m]*alpha[n][bp][m]
                                A[p,pp] += A_nbm
                        elif mod(pp,2) == 1: # A type param with B type param
                        # Possible indexing issues from body to body due to
                        # different number of surface atoms within cutoff.
                            bp = int((pp-1)/2)
                            if m >= len(beta[n][bp]):
                                pass
                            else:
                                A_nbm = -2*alpha[n][b][m]*beta[n][bp][m]
                                A[p,pp] += A_nbm

            # Same notes regarding indexing above apply here.
            elif mod(p,2) == 1: # B type parameter
                b = int((p-1)/2)
                B_nbm = -2*E_DFT[n] * beta[n,m,b]
                B[p] += B_nbm

                for pp in range(num_params):
                    if pp == p: # B type param with itself
                        bp = int((pp-1)/2)
                        A_nbm = 2*beta[n][b][m]*beta[n][bp][m]
                        A[p,pp] += A_nbm
                    else:
                        if mod(pp,2) == 0: # B type param with A type param
                            bp = int(pp/2)
                            if m >= len(alpha[n][bp]):
                                pass
                            else:
                                A_nbm = -2*beta[n][b][m]*alpha[n][bp][m]
                                A[p,pp] += A_nbm
                        elif mod(pp,2) == 1: # B type param with another B type param
                            bp = int((pp-1)/2)
                            if m >= len(beta[n][bp]):
                                pass
                            else:
                                A_nbm = 2*beta[n][b][m]*beta[n][bp][m]
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
print('---------------------------')

# Calculate the inverse of the coefficient matrix
# Calculate the solution
A_inv = np.linalg.inv(A)
AB_fit = np.matmul(A_inv,B)
print('\n\n----- Solutions -----')
print('AB_LJ = \n', AB_LJ)
print('AB_fit = \n', AB_fit)

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
