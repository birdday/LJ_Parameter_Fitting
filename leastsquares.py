# --------------------------------------------------
# ----- Import Python Packages ---------------------
# --------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

from pandas import read_csv
from math import sqrt
from ase import Atoms

# ----------------------------------------------------------------------------------------------------
# --------------------------------------------------
# ----- Function Parameters ------------------------
# --------------------------------------------------
# Molecule Representation
bodies = 1;
# Initial guess
params = [(1e-5, 5.3)]

# --------------------------------------------------
# ----- Experimental / Literature Values -----------
# --------------------------------------------------
graphene_vac = -305.49689755541885; graphene_pervac = -7.39805875741705;
sotolon_vac = -86.96520401259681; sotolon_pervac = -86.96759558341823;
carbon_sig = 3.55; carbon_eps = 0.066 / 627.509 #kcal/mol -> Hartree
surface_sig = carbon_sig; surface_eps = carbon_eps;

# --------------------------------------------------
# ----- 'Experimental' (DFT) Data ------------------
# --------------------------------------------------
# 'Experimental' Data
files = open('filenames.txt')
dirs = [line.strip('\n') for line in files.readlines()]
files.close()

# Read energies and radii
energy_file = open('energies.txt'); v = []; rad = []; lines = energy_file.readlines()
r_list = [float(l.strip('\n').split('\t')[1]) for l in lines]
e_list = [(float(l.strip('\n').split('\t')[2]) - sotolon_vac - graphene_vac) / 27.2114 for l in lines]
energy_file.close();

# Filter out high-energy (non-LJ) data points
e_min = min(e_list); directories = []; print(e_min)
for i, vv in enumerate(e_list):
    if vv < -0.007:
        v.append(vv); directories.append(dirs[i]); rad.append(r_list[i])

# Plot potential energy vs. distance
# N.B. Treats molecule of interest as a single body
# plt.plot(rad, v, 'go'); plt.show()
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------
# ----- User-defined Functions ---------------------
# --------------------------------------------------

# Get the distance between two objects
def distance_between_coords(p1, p2):
    dist = sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
    return dist

# Get the coordinates for the moleucle and surface
def get_coordinates(dir_lists, v, b, f):
    # Create DataFrame of coordinates for each configuration
    for j, dir in enumerate(tqdm(dir_list, leave=False)):
        df = read_csv(dir + '.xyz', header=2, names=['type', 'x', 'y', 'z'], sep='\s+')
        all_xyz = df[['x', 'y', 'z']]
        all_types = df['type']
        surf_z = coords['z'].iloc[0]; n = len(coords['z'])

    # Separate molecule and surface by checking z coordinate
    for i in range(n):
        if df['z'].iloc[i] != surf_z:
            molecule_xyz = coords.iloc[i:,:]
            molecule_types = types.iloc[i:]
            surface_xyz = coords.iloc[:i,:]
            surface_types = types.iloc[:i]

    return all_xyz, all_types, molecule_xyz, molecule_types, surface_xyz, surface_types

def create_ase_atoms(mol_coords, mol_types, surf_coords, surf_types):
    mol_atoms = Atoms(positions=mol_coords.values, symbols=mol_types.values)
    surf_atoms = Atoms(positions=surf_coords.values, symbols=surf_types.values)
    return mol_atoms, surf_atoms

# Calculate the 'internal' energy of the surface
def get_surface_energy(f, surf_list):
    E_surf = 0
    for i, p1 in enumerate(surf_list):
        for j, p2 in enumerate(surf_list[(i+1):]):
            dist = distance_between_coords(p1, p2)
            if dist > cutoff:
                continue;
            e = f.subs([(r, dist), (sig, surface_sig), (eps, surface_eps)])
            E_surf = E_surf + e
    return E_surf

# Calculate the 'internal' energy of the molecule
def get_molecule_energy(f, mol_list, surf_list):
    return E_molecule

# Calculate the energy between the molecule bodies and surface atoms
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

# --------------------------------------------------
# ----- Parameter Fitting --------------------------
# --------------------------------------------------
# Currently (11 April 2019, Brian Day) using a quasi-Newton method to optimize
# the multidimensional problem presented here.
# Other options include: (1) conjugate-gradient method, (2) ???
num_params = 2*bodies
A =[]; B =[];

# Calculate the first derivates of the MSE with respect to each variable
# Using the simplified form of the 12-6 Lennard-Jones equation, this yields:
for i in range(num_datapoints):
    for j in range(num_bodies):
        dMSE_Aij = ( A[] / (r[ij]**18) ) + ( B[] / r[ij] ** 12)
        dMSE_Bij = ( A[] / (r[ij]**18) ) + ( B[] / r[ij] ** 12)
        dMSE_Aij_total += dMSE_Aij
        dMSE_Bij_total += dMSE_Bij
    dMSE_Ai = E_DFT[i] + dMSE_Aij_total
    dMSE_Bi = E_DFT[i] + dMSE_Bij_total
dMSE_A_total = E_DFT[i] - dMSE_Ai
dMSE_B_total = E_DFT[i] - dMSE_Bi

# --------------------------------------------------
# ----- Extract Solution ---------------------------
# --------------------------------------------------
# Use mixing rules to extract body specific sigma and epsilon after parameter
# fitting. Multiple mixing rules exist which can vary the parameters calculated.
# Note that the mixing rules do NOT impact the fitted parameters.
