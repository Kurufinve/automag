import re
import sys
import os
import numpy as np
from ase.io import read
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

# print(np.array2string(magmoms, 10000))


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

    output_line += 'enthalpy={}\n'.format(float(enthalpy_line.split()[4]) + correction)
else:
    output_line += 'energy={}\n'.format(atoms_final.get_potential_energy() + correction)

if read_magmoms:
    with open('../singlepoint/initial_magmoms.txt','r') as f:
        mamgoms = np.array(f.readlines().split(''))
        # magmoms = np.array(self['initial_magmoms'])
    output_line += '          initial_magmoms={}  '.format(np.array2string(magmoms, 10000))

    try:
        magmoms_final = atoms_final.get_magnetic_moments()
    except:
        magmoms_final = np.zeros(len(magmoms))

    output_line += 'final_magmoms={}\n'.format(np.array2string(magmoms_final, 10000))

    filename = '_'.join(atoms_final.get_chemical_formula(mode='metal', empirical=True),'singlepoint.txt')

with open(os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold', filename), 'a') as f:
    f.write(output_line)



















with open(input_file, 'r') as file:
    file_content = file.read()

# Define the pattern to search for
pattern = r'MAGMOM = ([\d\.\*\s-]+)'
match = re.search(pattern, file_content)

if match:
    magmom_str = match.group(1)
    print("Found MAGMOM string:", magmom_str)

    # Replace the numbers after "=" with numbers from a NumPy array
    # Construct the new MAGMOM string
    new_magmom_str = f'{np.array2string(magmoms, 10000)}'.replace('[','').replace(']', '')
    new_content = file_content.replace(magmom_str, new_magmom_str)

    # Write the updated content back to the file
    with open(input_file, 'w') as file:
        file.write(new_content)
        print("Replacement completed and saved to the file.")
else:
    print("MAGMOM string not found in the file.")