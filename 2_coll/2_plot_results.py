"""
automag.2_coll.2_plot_results
=============================

Script which plots results of magnetic relaxations.

.. first codeauthor:: Michele Galasso <m.galasso@yandex.com>
.. second codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
use_fireworks = True
calculator = 'vasp'
create_cc_input = True # if create input.py for coupling constants calculation in 3_monte_carlo folder
create_mae_input = True # if create input.py for MAE calculation in 4_mae folder

import os,sys

cwd = os.getcwd()

try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    try:
        exec(f"from {input_file.split('.')[0]} import *")
    except:
        from input import *
except IndexError:
    print(f'Using the input.py file from folder: {cwd}')
    from input import *


import json
import shutil
import numpy as np
import matplotlib.pyplot as plt

from copy import copy
from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def better_sort(state):
    if state[0] == 'nm':
        return 'aa'
    elif state[0][:2] == 'fm':
        return 'aaa' + f'{int(state[0][2:]):3d}'
    else:
        return state[0][:3] + f'{int(state[0][3:]):3d}'


# take care of the case when lower_cutoff has not been specified
if 'lower_cutoff' not in globals():
    lower_cutoff = 0

# setting default tolerance in determining if magmoms can be considered the same
if 'tolerance' not in globals():
    tolerance = 0.08

# create an ase atoms object
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
atoms = read(path_to_poscar)

# get the multiplicities of each Wyckoff position
structure = Structure.from_file(path_to_poscar)
analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = analyzer.get_symmetrized_structure()
multiplicities = np.array([len(item) for item in symmetrized_structure.equivalent_indices])

# path to results file with data
calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
# output_file = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal', empirical=True)}_singlepoint.txt")
output_file = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}_singlepoint_{calculator}.txt")

print(f'Reading energies and mamgoms from file {output_file}')

trials_path = f"trials"

formula_calculator_template = f"{structure.formula.replace(' ','')}_{calculator}"


# exit if no trials folder
if not os.path.isdir(trials_path):
    print (f"No trials folder found. Seeking for trials_{structure.formula.replace(' ','')} folder")
    trials_path = f"trials_{structure.formula.replace(' ','')}"

    if not os.path.isdir(trials_path):
        print (f"No trials_{structure.formula.replace(' ','')} folder found. Seeking for {formula_calculator_template}/trials_{structure.formula.replace(' ','')} folder")
        trials_path = f"{formula_calculator_template}/trials_{structure.formula.replace(' ','')}"

        if not os.path.isdir(trials_path):
            raise IOError(f"No {structure.formula.replace(' ','')}/trials_{structure.formula.replace(' ','')} folder found.")

# read lines
lines, maginfos, presents = [], [], []
with open(output_file, 'rt') as f:
    for line in f:
        if line[0] != ' ':
            presents.append(line.split()[0])
            lines.append(line)
        else:
            maginfos.append(line)

data = {}
setting = 1
not_found = []
while os.path.isfile(f"{trials_path}/configurations{setting:03d}.txt"):
    with open(f"{trials_path}/configurations{setting:03d}.txt", 'rt') as f:
        for line in f:
            values = line.split()
            init_state = values[0]
            if init_state in presents:
                dct = {
                    'setting': setting,
                    'init_spins': [int(item) for item in values[1:]],
                }
                data[init_state] = dct
            else:
                not_found.append(init_state)
    setting += 1

# keep the maximum value of setting
max_setting = copy(setting)

# extract the results
red = []
not_converged = []
all_final_magmoms = []
for line, maginfo in zip(lines, maginfos):
    values = line.split()
    init_state = values[0]

    if values[-4].split('=')[1] == 'NONCONVERGED':
        not_converged.append(init_state)
        del data[init_state]
    else:
        initial, final = maginfo.split('final_magmoms=')
        initial = initial[initial.index('[') + 1:initial.index(']')]
        final = final[final.index('[') + 1:final.index(']')]
        initial = np.array(initial.split(), dtype=float)
        final = np.array(final.split(), dtype=float)

        all_final_magmoms.extend(final.tolist())

        # exclude low-spin configurations
        flag = False
        mask = np.nonzero(initial)
        if np.all(np.abs(final[mask]) > lower_cutoff) or init_state == 'nm':
            flag = True

        magnification = len(initial) // sum(multiplicities)
        current_multiplicities = magnification * multiplicities

        start = 0
        kept = True
        for multiplicity in current_multiplicities:
            w_initial = initial[start:start + multiplicity]
            w_final = final[start:start + multiplicity]
            start += multiplicity

            if not np.any(w_initial):
                if np.any(np.around(w_final)):
                    kept = False
            else:
                prod = np.multiply(np.sign(w_initial), w_final)
                # if prod.max() - prod.min() > 0.08:
                if prod.max() - prod.min() > tolerance:
                    kept = False

        if kept and flag:
            data[init_state]['kept_magmoms'] = True
        else:
            data[init_state]['kept_magmoms'] = False
            red.append(init_state)

        data[init_state]['energy'] = float(values[-1].split('=')[1]) / len(initial)     # energy per atom

if len(red) != 0:
    print(f"The following {len(red)} configuration(s) did not keep the original magmoms and will be marked in red "
          f"on the graph: {', '.join(red)}")
if len(not_converged) != 0:
    print(f"The energy calculation of the following {len(not_converged)} configuration(s) did not converge and will "
          f"not be shown on the graph: {', '.join(not_converged)}")
if len(not_found) != 0:
    print(f"The following {len(not_found)} configuration(s) reported an error during energy calculation and will "
          f"not be shown on the graph: {', '.join(not_found)}")

plt.hist(np.abs(all_final_magmoms), bins=40)
plt.savefig('spin_distribution.png')

setting = 1
final_states = []
final_setting = setting
final_energies = []
final_min = np.inf
while setting < max_setting:
    current_states = []
    current_energies = []
    for init_state, value in data.items():
        if init_state != 'nm' and value['setting'] == setting and value['kept_magmoms']:
            current_states.append(np.sign(value['init_spins']).tolist())
            current_energies.append(value['energy'])
    if min(current_energies) < final_min:
        final_setting = setting
        final_states = current_states
        final_energies = current_energies
        final_min = min(current_energies)
    setting += 1

# filter out configurations containing NM states
tc_states = []
tc_energies = []
for state, energy in zip(final_states, final_energies):
    if 0 not in state:
        tc_states.append(state)
        tc_energies.append(energy)

# write states to file
with open(f"{formula_calculator_template}_states{final_setting:03d}.txt", 'wt') as f:
    json.dump(tc_states, f)

# write energies to file
with open(f"{formula_calculator_template}_energies{final_setting:03d}.txt", 'wt') as f:
    json.dump(tc_energies, f)

# copy setting file with geometry
shutil.copy(f"{trials_path}/setting{final_setting:03d}.vasp", f"./{formula_calculator_template}_setting{final_setting:03d}.vasp")

# extract values for plot
bar_labels = []
energies = []
kept_magmoms = []
for key, value in sorted(data.items(), key=better_sort):
    bar_labels.append(key)
    energies.append(value['energy'])
    kept_magmoms.append(value['kept_magmoms'])

energies = np.array(energies)
kept_magmoms = np.array(kept_magmoms)

# energies from eV/atom to meV/atom
energies *= 1000

# normalize energies for plot
energies -= min(energies)
toplim = max(energies) * 1.15
bottomlim = -0.1 * max(energies)

energies += 0.1 * max(energies)
x_axis = np.arange(1, len(energies) + 1)

# split into chunks if too many configurations
n_chunks = len(energies) // 14 + 1
x_axis_split = np.array_split(x_axis, n_chunks)
energies_split = np.array_split(energies, n_chunks)
kept_magmoms_split = np.array_split(kept_magmoms, n_chunks)
bar_labels_split = np.array_split(bar_labels, n_chunks)

for i, (X, Y, kept_magmoms_chunk, bar_labels_chunk) in \
        enumerate(zip(x_axis_split, energies_split, kept_magmoms_split, bar_labels_split)):
    # increase matplotlib pyplot font size
    plt.rcParams.update({'font.size': 20})

    # set figure size
    plt.figure(figsize=(16, 9))

    # plot results
    plt.bar(X[~kept_magmoms_chunk], Y[~kept_magmoms_chunk], bottom=bottomlim, color='r')
    plt.bar(X[kept_magmoms_chunk], Y[kept_magmoms_chunk], bottom=bottomlim, color='b')

    # label bars
    ax = plt.gca()
    rects = ax.patches
    rects = sorted(rects, key=lambda x: x.get_x())
    for bar_label, rect in zip(bar_labels_chunk, rects):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2, height + 0.85 * bottomlim, bar_label,
                fontsize='medium', ha='center', va='bottom', rotation='vertical')

    # label axes
    plt.xlabel('configurations')
    plt.ylabel(f'relative energy (meV/atom)')
    plt.xticks(X)
    plt.ylim(top=toplim)

    # save or show bar chart
    plt.savefig(f"{formula_calculator_template}_stability{i + 1:02d}.png", bbox_inches='tight')
    plt.close()

print(f'The most stable configuration is {bar_labels[np.argmin(energies)]}.')
with open(f'{formula_calculator_template}_most_stable_configuration.txt','w') as f:
    f.write(bar_labels[np.argmin(energies)])

with open(f'{formula_calculator_template}_configurations_stability_ascending.txt','w') as f:
    for i in np.argsort(energies):
        f.write(f'{bar_labels[i]} {energies[i]}\n')



if np.logical_not(kept_magmoms[np.argmin(energies)]):
    print('WARNING: values of initial and final magnetic moments of the most stable configuration '
          'significantly differ.')


# moving files and figures into a separate folder
os.system(f"mkdir {formula_calculator_template}")
os.system(f"cp -f input.py {formula_calculator_template}/input_{calculator}.py")
os.system(f"mv -f {formula_calculator_template}_stability*.png {formula_calculator_template}")
try:
    os.system(f"mv -f trials_{structure.formula.replace(' ','')} {formula_calculator_template}")
except:
    pass
os.system(f"mv -f {formula_calculator_template}_setting*.vasp {formula_calculator_template}")
os.system(f"mv -f {formula_calculator_template}_energies*.txt {formula_calculator_template}")
os.system(f"mv -f {formula_calculator_template}_states*.txt {formula_calculator_template}")
os.system(f"mv -f {formula_calculator_template}_most_stable_configuration.txt {formula_calculator_template}")
os.system(f"mv -f {formula_calculator_template}_configurations_stability_ascending.txt {formula_calculator_template}")

# acrhiving the folder
os.system(f"tar -zcvf {formula_calculator_template}.tar.gz {formula_calculator_template}")
