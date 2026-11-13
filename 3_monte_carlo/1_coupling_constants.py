"""
automag.3_monte_carlo.1_coupling_constants.py
=============================================

Script which computes the coupling constants between magnetic atoms.

.. first codeauthor:: Michele Galasso <m.galasso@yandex.com>
.. second codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
calculator = 'vasp'
struct_suffix = ''

import os,sys

cwd = os.getcwd()
struct_suffix = ''

from common.input_loader import apply_stage_input

apply_stage_input(globals())


import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from pymatgen.core.structure import Structure

from common.results import coll_bundle
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
path_to_coll = str(coll_bundle(path_to_automag, formula, struct_suffix,
                               coll_calculator))
print(f'Path to the results from collinear calculations: {path_to_coll}')



# initialize variables
structure = None
states = None
energies = None

# read input files from previous step
# 2_plot_results.py names the artefacts after the collinear calculator.
for item in os.listdir(path_to_coll):
    rel_path = os.path.join(path_to_coll, item)
    if os.path.isfile(rel_path):
        if item.startswith(f'{formula}{struct_suffix}_{coll_calculator}_setting') and item.endswith('.vasp'):
            structure = Structure.from_file(rel_path)
        if item.startswith(f'{formula}{struct_suffix}_{coll_calculator}_states') and item.endswith('.txt'):
            with open(rel_path, 'rt') as f:
                states = json.load(f)
        if item.startswith(f'{formula}{struct_suffix}_{coll_calculator}_energies') and item.endswith('.txt'):
            with open(rel_path, 'rt') as f:
                energies = json.load(f)

if structure is None:
    raise IOError(f'No setting file found in {path_to_coll} folder.')
if states is None:
    raise IOError(f'No states file found in {path_to_coll} folder.')
if energies is None:
    raise IOError(f'No energies file found in {path_to_coll} folder.')

# find out which atoms are magnetic: 'magnetic_atoms' is optional and
# otherwise derived from spin_values / pymatgen's DEFAULT_MAGMOMS.
magnetic_atoms = resolve_magnetic_elements(
    structure, globals().get('spin_values'), globals().get('magnetic_atoms'))
for element in structure.composition.elements:
    element.is_magnetic = element.name in magnetic_atoms
print(f'Magnetic elements: {magnetic_atoms}')

# from eV/atom to total energy of the unit cell
energies = structure.num_sites * np.array(energies)



# The number of unknowns grows with cutoff_radius (one coupling per neighbour
# shell plus the reference energy), while the number of independent equations
# is fixed by the magnetic configurations available.  When the system is
# underdetermined the cutoff is shrunk by CUTOFF_STEP and refitted, never below
# the nearest-neighbour shell.
CUTOFF_STEP = 0.1


def fit_model(cutoff):
    """Least-squares fit of the Heisenberg model for one cutoff radius."""
    center_indices, point_indices, offset_vectors, distances = \
        structure.get_neighbor_list(cutoff)
    shells, counts = np.unique(np.around(distances, 3), return_counts=True)

    def system(configurations):
        matrix = []
        for item in configurations:
            equation = [1]
            for distance in shells:
                count = 0
                for atom1, atom2, d in zip(center_indices, point_indices, distances):
                    if np.isclose(d, distance, atol=0.02):
                        count += item[atom1] * item[atom2]
                equation.append(-count // 2)
            matrix.append(equation)
        return np.array(matrix)

    A = system(configurations_fit)
    values = np.linalg.lstsq(A, energies_fit, rcond=None)
    return {
        'cutoff': float(cutoff),
        'shells': shells,
        'counts': counts,
        'rank': int(np.linalg.matrix_rank(A)),
        'unknowns': int(len(shells) + 1),
        'values': values,
        'system': system,
        'A': A,
    }


non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
structure.remove_species(non_magnetic_atoms)

# create fit and control group
fit_group_size = 1 - control_group_size
index = int(round(len(states) * fit_group_size))
configurations_fit = np.array(states[:index])
configurations_control = np.array(states[index:])
energies_fit = energies[:index]
energies_control = energies[index:]

cutoff_radius_requested = float(cutoff_radius)
nearest_distance = float(np.min(
    structure.get_neighbor_list(cutoff_radius_requested)[3]))

model = fit_model(cutoff_radius_requested)
cutoff = cutoff_radius_requested
while ((model['rank'] < model['unknowns'] or len(model['shells']) < 1)
       and cutoff - CUTOFF_STEP >= nearest_distance):
    print(f'WARNING: system of {model["rank"]} independent equation(s) in '
          f'{model["unknowns"]} unknowns at cutoff_radius = {cutoff:.2f} A.')
    cutoff = round(cutoff - CUTOFF_STEP, 10)
    print(f'         reducing cutoff_radius to {cutoff:.2f} A and refitting')
    model = fit_model(cutoff)

unique_distances = model['shells']
counts = model['counts']
cutoff_radius_used = model['cutoff']
if abs(cutoff_radius_used - cutoff_radius_requested) > 1e-9:
    print(f'WARNING: cutoff_radius was reduced from {cutoff_radius_requested:.2f} '
          f'to {cutoff_radius_used:.2f} A so that {len(unique_distances)} '
          f'coupling(s) can be fitted with {model["rank"]} independent '
          f'equation(s).')

values = model['values']
B = model['system'](configurations_control)
predictions = B @ values[0]
PCC = np.corrcoef(predictions, energies_control)

plt.rcParams.update({'font.size': 13})
plt.locator_params(axis='x', nbins=5)
plt.xlabel(f'Heisenberg model energy (eV)')
plt.ylabel(f'DFT energy (eV)')

# print results
if model['rank'] == model['unknowns'] and len(unique_distances) >= 1:
    coupling_constants = values[0][1:] * 1.60218e-19
    print(f'distances between neighbors: {unique_distances.tolist()}')
    print(f'counts: {counts.tolist()}')
    print(f'coupling constants: {np.array2string(coupling_constants, precision=8, separator=", ")}')

    print(f'append_coupling_constants = {append_coupling_constants}')
    if append_coupling_constants:
        print('Appending coupling constants into input.py')
        with open('input.py', 'a') as f:
            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ndistances_between_neighbors = {unique_distances.tolist()}\n')

            f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
            f.write(f'\ncoupling_constants = {np.array2string(coupling_constants, precision=8, separator=", ")}\n')

            if abs(cutoff_radius_used - cutoff_radius_requested) > 1e-9:
                f.write('\n# LINE ADDED BY THE SCRIPT 1_coupling_constants.py')
                f.write('\n# cutoff_radius was reduced to fit the Heisenberg model')
                f.write(f'\ncutoff_radius = {cutoff_radius_used}\n')

    plt.scatter(predictions, energies_control, label=f'PCC: {PCC[0, 1]:.2f}')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    plt.legend()
    # plt.show()
    plt.savefig('model.png', bbox_inches='tight')
    print(f'PCC: {PCC[0, 1]:.2f}')
else:
    # Still underdetermined at the smallest sensible cutoff: nothing can be
    # passed to the Monte Carlo backend.
    print(f'ERROR: SYSTEM OF {model["rank"]} INDEPENDENT EQUATION(S) '
          f'IN {model["unknowns"]} UNKNOWNS!')
    print(f'       cutoff_radius was already reduced from '
          f'{cutoff_radius_requested:.2f} to {cutoff_radius_used:.2f} A '
          f'(nearest neighbour is at {nearest_distance:.2f} A).')
    print()
    print('The Heisenberg model cannot be determined from these magnetic')
    print('configurations, so there is nothing to pass to the Monte Carlo')
    print('backend.  Ways out:')
    print('  * generate more configurations: raise "supercell_size" in the')
    print('    2_coll stage (e.g. 2 instead of 1), which multiplies the number')
    print('    of magnetic atoms and orderings;')
    print('  * lower "cutoff_radius" in the monte_carlo stage explicitly;')
    print('  * or run the pipeline without the monte_carlo stage.')
    raise SystemExit(1)
