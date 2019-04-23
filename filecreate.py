# Orinally written by Spencer Smith
# Modifications by Brian Day
# Last updated: 19 April 2019

# --------------------------------------------------
# ----- Import Python Packages ---------------------
# --------------------------------------------------
import ase as ase
from ase import io, spacegroup, build, visualize
import numpy as np
import math as math
from math import cos, sin, pow
from os import mkdir, chdir
from random import uniform
from shutil import copy

# --------------------------------------------------
# ----- Manual inputs ------------------------------
# --------------------------------------------------
# Read in molecule xyz file and surface cif file
mol_file = 'molecules/p_sarin.mol'
mol_orig = ase.io.read(mol_file)
mol_name = mol_file.split('/')[1].split('.')[0]

surf_file ='surfaces/graphite.cif'
surf_cell = ase.io.read(surf_file)
surf_name = surf_file.split('/')[1].split('.')[0]
xyz_file_dir = mol_name+'_'+surf_name
xyz_file_dir_images = mol_name+'_'+surf_name+'_images'

surf_cell_spacegroup = 'P 1'
buffer = 8
N_configs = 50
zrange = [1,8]

# --------------------------------------------------
# ----- User defined functions ---------------------
# --------------------------------------------------
# Calculate the center of positions (non-weighted center of mass)
def get_center_of_positions(self, scaled=False):
    num_atoms = len(self)
    all_xyz = self.get_positions()
    avg_xyz = np.sum(all_xyz, axis=0)/num_atoms
    return avg_xyz

# Calculate the remiander
def mod(a,b):
    remainder = a%b
    return remainder

# --------------------------------------------------
# ----- Prepare surface slab -----------------------
# --------------------------------------------------
# Get surface crystal properties
surf_cell_params = surf_cell.get_cell_lengths_and_angles()
surf_cell_params = np.round(1E8*np.array(surf_cell_params))*1E-8
[a, b, c, alpha, beta, gamma] = surf_cell_params
[alpha, beta, gamma] = np.deg2rad([alpha, beta, gamma])
surf_cell_vects = np.array([a,b,c])
surf_cell_angles = np.array([alpha, beta, gamma])
surf_cell_atomtypes = surf_cell.get_chemical_symbols()
surf_cell_frac_coords = surf_cell.get_scaled_positions()

# Define python ase equivalent of cif file
surf_crystal = ase.spacegroup.crystal(
surf_cell_atomtypes, surf_cell_frac_coords,
spacegroup=surf_cell_spacegroup, cellpar=surf_cell_params)

# Define an orthoganalization matrix - Is this general for use with cut?
cy = (cos(alpha)-cos(gamma)*cos(beta))/sin(gamma)
M_orth = np.array([[1,0,0],[-cos(gamma),1,0],[cos(beta),cy,1]])

# Calculate box distance for a given molecule size / buffer
# Currently a rough approximate of the size of the molecule, fix later
dists = np.max(mol_orig.positions.max(0) - mol_orig.positions.min(0))
box_size = dists + 2*buffer;

# Calculate number of unit cells
# Must be even number so file is periodic
num_cells = np.ceil(box_size/surf_cell_vects)
num_cells
for i in range(len(num_cells)):
    if mod(num_cells[i],2) == 1:
        num_cells[i] = num_cells[i]+1

# Build orthogonalized surface slab
# N.B. - If seemingly stuck in loop while using cut, set tolerance parameter to 0.
# Should not be a problem for most common crystals.
ao = np.int(num_cells[0])*M_orth[0,:]
bo = np.int(num_cells[1])*M_orth[1,:]
co = np.int(num_cells[2])*M_orth[2,:]
surf_slab_orig = ase.build.cut(surf_crystal, a=ao, b=bo, origo=(0,0,0), nlayers=1, tolerance=0.01)

# --------------------------------------------------
# ----- Align molecule and surface -----------------
# --------------------------------------------------
# Place center of surface and molecule at origin
surf_shift = get_center_of_positions(surf_slab_orig)
mol_shift = get_center_of_positions(mol_orig)
mol_orig.translate(surf_shift-mol_shift)

# --------------------------------------------------
# ----- Create all configurations ------------------
# --------------------------------------------------
mkdir(xyz_file_dir)
mkdir(xyz_file_dir_images)
for j in range(1,N_configs+1):
    # Reset molecule position / Reset surface file
    mol = mol_orig.copy()
    surf_slab = surf_slab_orig.copy()

    # Generate random rotation angles / distances
    # N.B. If trying to create additional orientations, shift the range at the
    # top of the for loop, otherwise, with random seeding, you will simply
    # generate the exact same set of configurations.
    rand_seed = np.random.seed(j)
    rand_dist = np.random.rand(1)*(zrange[1]-zrange[0])+zrange[0]
    rand_angles = np.random.rand(3)*360
    phi = rand_angles[0]
    theta = rand_angles[1]
    psi = rand_angles[2]

    # Rotate / Translate the molecule
    # Translate after rotating to eliminate clipping of molecule/surface.
    # Since center of molecule is at origin, minimum z will be negative.
    mol.euler_rotate(phi=phi, theta=theta, psi=psi, center='COP')
    pos_after_rot = mol.get_positions()
    extra_dist = np.min(pos_after_rot[:,2])
    mol.translate([0,0,rand_dist-extra_dist])

    # Join molecule and surface, Write the xyz file
    filename = (xyz_file_dir+'_'+str(j))
    output_file = mol.extend(surf_slab)
    chdir(xyz_file_dir)
    ase.io.write(filename+'.xyz',output_file)
    chdir('../'+xyz_file_dir_images)
    ase.io.write(filename+'.png',output_file, rotation='-90x')
    chdir('..')
