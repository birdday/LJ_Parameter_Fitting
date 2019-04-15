# --------------------------------------------------
# ----- Import Python Packages ---------------------
# --------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

from pandas import read_csv
from math import sqrt
# from ase import Atoms

# # --------------------------------------------------
# # ----- User-defined Functions ---------------------
# # --------------------------------------------------
# # Get the distance between two objects
# def distance_between_coords(p1, p2):
#     dist = sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
#     return dist
#
# # Get the coordinates for the moleucle and surface
# def get_coordinates(dir_lists, v, b, f):
#     # Create DataFrame of coordinates for each configuration
#     for j, dir in enumerate(tqdm(dir_list, leave=False)):
#         df = read_csv(dir + '.xyz', header=2, names=['type', 'x', 'y', 'z'], sep='\s+')
#         all_xyz = df[['x', 'y', 'z']]
#         all_types = df['type']
#         surf_z = coords['z'].iloc[0]; n = len(coords['z'])
#
#     # Separate molecule and surface by checking z coordinate
#     for i in range(n):
#         if df['z'].iloc[i] != surf_z:
#             molecule_xyz = coords.iloc[i:,:]
#             molecule_types = types.iloc[i:]
#             surface_xyz = coords.iloc[:i,:]
#             surface_types = types.iloc[:i]
#
#     return all_xyz, all_types, molecule_xyz, molecule_types, surface_xyz, surface_types
#
# def create_ase_atoms(mol_coords, mol_types, surf_coords, surf_types):
#     mol_atoms = Atoms(positions=mol_coords.values, symbols=mol_types.values)
#     surf_atoms = Atoms(positions=surf_coords.values, symbols=surf_types.values)
#     return mol_atoms, surf_atoms
#
# # Calculate the 'internal' energy of the surface
# def get_surface_energy(f, surf_list):
#     E_surf = 0
#     for i, p1 in enumerate(surf_list):
#         for j, p2 in enumerate(surf_list[(i+1):]):
#             dist = distance_between_coords(p1, p2)
#             if dist > cutoff:
#                 continue;
#             e = f.subs([(r, dist), (sig, surface_sig), (eps, surface_eps)])
#             E_surf = E_surf + e
#     return E_surf
#
# # Calculate the 'internal' energy of the molecule
# def get_molecule_energy(f, mol_list, surf_list):
#     E_mol = 1;
#     return E_mol

# Calculate the energy between the molecule bodies and surface atoms
f = 1; cutoff = 6;
def get__mixed_energy(f, mol_list, surf_list):
    E_mixed = 0
    for i, p1 in enumerate(mol_list):
        for j, p2 in enumerate(surf_list):
            dist = distance_between_coords(p1, p2)
            if dist > cutoff:
                continue;
            e = f.subs(r, dist)
            E_mixed = E_mixed + e
    return E_mixed

# Calculate the remiander
def mod(a,b):
    remainder = a%b
    return remainder

# --------------------------------------------------
# ----- 'Experimental' (DFT) Data ------------------
# --------------------------------------------------
# Test Data (i.e. Perfect Fit) - 1 Body
num_configs = 101
num_surfatoms = 1
num_bodies = 1
num_params = 2*num_bodies

sigma = 3.5
epsilon = 0.066
A_LJ = 4 * epsilon * (sigma**12)
B_LJ = 4 * epsilon * (sigma**6)
AB_LJ = np.array([A_LJ,B_LJ]).reshape(2,1)

rlim = ([3,8])
rvect = np.ones((num_configs, num_surfatoms, num_bodies))
rvect[:,0,0] = np.linspace(rlim[0], rlim[1], num=num_configs, endpoint=True)
alpha = rvect**(-12)
beta = rvect**(-6)

E_DFT = np.ones((num_configs))
E_DFT = alpha[:,0,0]*A_LJ - beta[:,0,0]*B_LJ

# Plot the Test Data
plt.plot(rvect[:,0,0],E_DFT[:],'g-o')
plt.xlabel('Distance [Angstrom]')
plt.ylabel('Energy [kcal/mol]')
plt.title('Lennard-Jones Test Data')
plt.show()

# --------------------------------------------------
# ----- Parameter Fitting --------------------------
# --------------------------------------------------
# Currently (14 April 2019, Brian Day) attempting to set up a system of linear
# equations which can be solved thorugh Gaussian elimination. The minimuzation
# function is the Sum of Squared Errors (SSE).
A = np.zeros((num_params,num_params))
B = np.zeros((num_params,1))

# Calculate the A and B Matricies for the system of equations to be solved.
# Form of equation: Ax = B --> Solution: x = inv(A)*B
for n in range(num_configs):
    for m in range(num_surfatoms):
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
