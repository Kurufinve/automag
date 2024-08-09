import re
import sys
import os
import numpy as np
from ase.io import read
from ase.calculators.vasp import Vasp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.io.vasp import Incar, Outcar

# default parameters for reading vasp output
read_energy_convergence = False
read_enthalpy = False
read_magmoms = False

try:
    mode = str(sys.argv[1])
except:
    mode = 'magmoms' # default mode is read magmoms

if mode == 'magmoms':
    read_enthalpy = True
    read_magmoms = True

filename = 'singlepoint.txt'

atoms_final = read('OUTCAR')

formula = atoms_final.get_chemical_formula(mode='metal')



# print(np.array2string(magmoms, 10000))

cwd = os.getcwd()
state = cwd.split('/')[-2]
# calculator = cwd.split('/')[-3]
calculatro = 'vasp'
compound = cwd.split('/')[-4]
mode = cwd.split('/')[-1]
try:
    calcid = re.search(r'\d+', state).group()
except:
    calcid = 0

if mode == 'singlepoint':
    calc_singlepoint = Vasp(directory=cwd)
    is_singlepoint_converged = calc_singlepoint.read_convergence()
    if is_singlepoint_converged: singlepoint_converged_string = 'convergedconverged'
    elif not is_singlepoint_converged: singlepoint_converged_string = 'NONCONVERGEDNONCONVERGED'
    output_line = f'{state}       {calcid}  singlepoint={singlepoint_converged_string}  {formula}   '
elif mode == 'recalc':
    calc_singlepoint = Vasp(directory=cwd.replace('recalc','singlepoint'))
    is_singlepoint_converged = calc_singlepoint.read_convergence()
    calc_recalc = Vasp(directory=cwd)
    is_recalc_converged = calc_recalc.read_convergence()
    # print(f'is_singlepoint_converged: {is_singlepoint_converged}')
    # print(f'is_recalc_converged: {is_recalc_converged}')
    if is_singlepoint_converged: singlepoint_converged_string = 'convergedconverged'
    elif not is_singlepoint_converged: singlepoint_converged_string = 'NONCONVERGEDNONCONVERGED'
    if is_recalc_converged: recalc_converged_string = 'convergedconverged'
    elif not is_recalc_converged: recalc_converged_string = 'NONCONVERGEDNONCONVERGED'
    output_line = f'{state}       {calcid}  singlepoint={singlepoint_converged_string}  recalc={recalc_converged_string}   {formula}    '


structure = Structure.from_file('CONTCAR')
analyzer = SpacegroupAnalyzer(structure)
output_line += '{:10s}  '.format(analyzer.get_space_group_symbol())


if read_energy_convergence:
    errors = []
    with open('OUTCAR', 'r') as f:
        for line in f:
            if 'kinetic energy error' in line:
                errors.append(float(line.split()[5]))

    _, indices, counts = np.unique(atoms_final.numbers, return_index=True, return_counts=True)
    num_ions = counts[np.argsort(indices)]
    correction = sum(np.multiply(errors, num_ions))
else:
    correction = 0

if read_enthalpy:
    # get enthalpy from OUTCAR
    with open('OUTCAR', 'r') as f:
        for line in f:
            if 'enthalpy' in line:
                enthalpy_line = line
    try:
        output_line += 'enthalpy={}\n'.format(float(enthalpy_line.split()[4]) + correction)
    except:
        output_line += 'energy={}\n'.format(atoms_final.get_potential_energy() + correction)
else:
    output_line += 'energy={}\n'.format(atoms_final.get_potential_energy() + correction)

if read_magmoms:
    with open('../singlepoint/initial_magmoms.txt','r') as f:
        lines = f.readlines()
        print(lines[0])
        magmoms = np.array([float(i) for i in lines[0].split()])
    output_line += '          initial_magmoms={}  '.format(np.array2string(magmoms, 10000))

    try:
        magmoms_final = atoms_final.get_magnetic_moments()
    except:
        magmoms_final = np.zeros(len(magmoms))

    output_line += 'final_magmoms={}\n'.format(np.array2string(magmoms_final, 10000))

    filename = f"{atoms_final.get_chemical_formula(mode='metal')}_singlepoint_{calculator}.txt"

with open(os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold', filename), 'a') as f:
    f.write(output_line)