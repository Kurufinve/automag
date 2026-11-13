"""
automag.3_monte_carlo.2_run_vampire.py
======================================

Script which runs Vampire for Monte Carlo simulation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

# default values for some input variables
calculator = 'vasp'
struct_suffix = ''
script_header = '#!/bin/bash'

import os,sys

cwd = os.getcwd()

from common.input_loader import apply_stage_input

apply_stage_input(globals())

    
import numpy as np

from pymatgen.core.structure import Structure

from common.results import find_trials_dir, resolve_configuration, final_magmoms
from common.magnetic import resolve_magnetic_elements

# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
input_structure = Structure.from_file(path_to_poscar)

formula = input_structure.formula.replace(' ','')
path_to_automag = os.environ.get('AUTOMAG_PATH')
# The collinear results may have been produced with a different calculator
# than the one used here (mixed-calculator pipeline).
coll_calculator = globals().get('coll_calculator') or calculator
configuration = resolve_configuration(
    configuration, path_to_automag, formula, struct_suffix, coll_calculator)
path_to_trials = str(find_trials_dir(
    path_to_automag, formula, struct_suffix, coll_calculator))
print(f'Path to trials: {path_to_trials}')


setting = 1
magmom = None
while os.path.isfile(f'{path_to_trials}/configurations{setting:03d}.txt'):
    with open(f'{path_to_trials}/configurations{setting:03d}.txt', 'rt') as f:
        for line in f:
            values = line.split()
            if values[0] == configuration:
                magmom = [int(item) for item in values[1:]]
                break
        if magmom is not None:
            break
    setting += 1

if magmom is None:
    raise IOError(f'configuration {configuration} not found.')

# get materials from magmom
materials = [0 if item >= 0 else 1 for item in magmom]

# create a pymatgen Structure object
structure = Structure.from_file(f'{path_to_trials}/setting{setting:03d}.vasp')

# find out which atoms are magnetic: 'magnetic_atoms' is optional and
# otherwise derived from spin_values / pymatgen's DEFAULT_MAGMOMS.
magnetic_atoms = resolve_magnetic_elements(
    structure, globals().get('spin_values'), globals().get('magnetic_atoms'))
for element in structure.composition.elements:
    element.is_magnetic = element.name in magnetic_atoms
print(f'Magnetic elements: {magnetic_atoms}')

# The converged moments of the calculation are reported for the full cell,
# while 'structure' keeps only the magnetic atoms - remember the mask.
magnetic_mask = np.array([site.specie.name in magnetic_atoms
                          for site in structure])
non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)
print('Magnetic part of the structure')
print(structure)
center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)

# get thresholds for comparing atomic distances
# 'distances_between_neighbors' and 'coupling_constants' are appended to
# input.py by 1_coupling_constants.py when the Heisenberg model could be
# fitted; without them the VAMPIRE unit cell cannot be written.
missing = [name for name in ('distances_between_neighbors', 'coupling_constants')
           if name not in globals()]
if missing:
    raise SystemExit(
        'ERROR: ' + ' and '.join(missing) + ' not found in input.py.\n'
        '       1_coupling_constants.py did not fit the Heisenberg model, so\n'
        '       there are no coupling constants to write into vamp.ucf.\n'
        '       See its output above (usually "SYSTEM OF N EQUATION(S) IN M\n'
        '       UNKNOWNS"): raise supercell_size in the 2_coll stage, lower\n'
        '       cutoff_radius, or drop the monte_carlo stage.')
thresholds = [(a + b) / 2 for a, b in zip(distances_between_neighbors[:-1], distances_between_neighbors[1:])]
thresholds = [0.0] + thresholds + [100.0]

try:
    os.mkdir('vampire')
except:
    pass
# print unit cell size for VAMPIRE input
with open('vampire/vamp.ucf', 'w') as f:
    f.write('# Unit cell size:\n')
    f.write(f'{structure.lattice.a:19.16f}  {structure.lattice.b:19.16f}  {structure.lattice.c:19.16f}\n')

    f.write('# Unit cell vectors:\n')
    for i, (vector, modulus) in enumerate(zip(structure.lattice.matrix, structure.lattice.abc)):
        v = vector / modulus
        f.write(f'{v[0]: 10.8e}   {v[1]: 10.8e}   {v[2]: 10.8e}\n')

    # print fractional coordinates for VAMPIRE input
    f.write('# Atoms num_atoms num_materials; id cx cy cz mat cat hcat\n')
    f.write(f'{len(structure):d} {len(np.unique(materials))}\n')
    for i, (coord, material) in enumerate(zip(structure.frac_coords, materials)):
        coord += 0.     # gets rid of the minus zero values
        for j, item in enumerate(coord):
            if item < 0:
                coord[j] += 1
            elif item >= 1:
                coord[j] -= 1
        f.write(f'{i:2d}   {coord[0]:18.16f}  {coord[1]:18.16f}  {coord[2]:18.16f}  {material} 0 0\n')

    f.write('# Interactions n exctype; id i j dx dy dz Jij\n')
    f.write(f'{len(center_indices)} isotropic\n')
    for i, (atom1, atom2, offset, distance) in enumerate(zip(center_indices, point_indices, offset_vectors, distances)):
        for low_lim, high_lim, coupling_constant in zip(thresholds[:-1], thresholds[1:], coupling_constants):
            if low_lim < distance < high_lim:
                f.write(f'{i:3d}   {atom1:2d}  {atom2:2d}  '
                        f'{int(offset[0]):2d} {int(offset[1]):2d} {int(offset[2]):2d}   {coupling_constant: 6.4e}\n')

# ---------------------------------------------------------------------------
# vamp.mat carries the *converged* magnetic moment of every sublattice.  The
# bundled template has 4.2 muB hard-coded, which is wrong for any real
# calculation: take the moments of the chosen configuration from the 2_coll
# results (the same source that 4_mae uses).
# ---------------------------------------------------------------------------
moments = final_magmoms(path_to_automag, formula, configuration,
                        struct_suffix, coll_calculator)
if moments is not None:
    moments = np.asarray(moments, dtype=float)
    if len(moments) == len(magnetic_mask):
        moments = moments[magnetic_mask]
    else:
        print(f'WARNING: the 2_coll results list {len(moments)} moments while '
              f'the cell has {len(magnetic_mask)} sites; using spin_values')
        moments = None
if moments is None:
    # Fall back to the requested spin values (no calculation to read).
    moments = np.array([float((globals().get('spin_values') or {}).get(
        site.specie.symbol, [0])[0]) for site in structure])
    print('WARNING: no converged moments found, vamp.mat uses the requested '
          'spin values')

material_ids = np.array(materials)
number_of_materials = int(len(np.unique(material_ids)))
with open('vampire/vamp.mat', 'w') as f:
    f.write('#---------------------------------------------------\n')
    f.write('# Number of Materials\n')
    f.write('#---------------------------------------------------\n')
    f.write(f'material:num-materials={number_of_materials}\n')
    for index in range(number_of_materials):
        atoms = np.where(material_ids == index)[0]
        element = (structure[int(atoms[0])].specie.symbol if len(atoms)
                   else sorted(magnetic_atoms)[0])
        if len(atoms):
            moment = float(np.mean(np.abs(moments[atoms])))
        else:
            moment = 0.0
        suffix = 'up' if index % 2 == 0 else 'dn'
        direction = '0,0,1' if index % 2 == 0 else '0,0,-1'
        f.write('#---------------------------------------------------\n')
        f.write(f'# Material {index + 1} {element}\n')
        f.write('#---------------------------------------------------\n')
        f.write(f'material[{index + 1}]:material-name={element}-{suffix}\n')
        f.write(f'material[{index + 1}]:atomic-spin-moment={moment:.4f} !muB\n')
        f.write(f'material[{index + 1}]:material-element={element}\n')
        f.write(f'material[{index + 1}]:uniaxial-anisotropy-constant=0.0\n')
        f.write(f'material[{index + 1}]:initial-spin-direction={direction}\n')
summary = ', '.join(
    f'{np.abs(moments[material_ids == index]).mean():.3f}'
    for index in range(number_of_materials))
print(f'vamp.mat written: {number_of_materials} material(s), '
      f'moments {summary} muB')
