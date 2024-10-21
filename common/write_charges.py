import re
import sys
import os
import shutil
import numpy as np
from ase.io import read
from ase.calculators.vasp import Vasp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Incar, Outcar

# small helping function
def is_numeric(value):
    try:
        float(value)  # Try to convert to float
        return True
    except ValueError:
        return False

# default parameters for reading vasp output
read_energy_convergence = False
read_enthalpy = False
read_magmoms = False

try:
    dummy_position = int(sys.argv[1])
except:
    dummy_position = 0 # default calculator is vasp
    
print(f'dummy_position is {dummy_position}')


try:
    calculator = str(sys.argv[2])
except:
    calculator = 'vasp' # default calculator is vasp


# filename = 'singlepoint.txt'

if calculator == 'vasp':
    incar_bare = Incar.from_file('singlepoint/INCAR')
    outcar_bare = Outcar('singlepoint/OUTCAR')
    bare_magmoms = [item['tot'] for item, ref in zip(outcar_bare.magnetization, incar_bare['MAGMOM']) if ref != 0]
    atoms_final = read(f'singlepoint/OUTCAR')

# atoms_final_nsc_list = []
# atoms_final_sc_list  = []

listdir = os.listdir('.')

pert_folders = [float(item) for item in listdir if is_numeric(item)]

# Sort the numeric folders in ascending order
pert_folders.sort()

# If you want to keep the original string representation, convert back to strings
sorted_pert_folders = [str(num) for num in pert_folders]

filename = f"{atoms_final.get_chemical_formula(mode='metal')}_charges_{calculator}.txt"

for pert_folder in sorted_pert_folders:
    if calculator == 'vasp':
        # atoms_final_nsc = read(f'{pert_folder}/nsc/OUTCAR')
        # atoms_final_sc = read(f'{pert_folder}/sc/OUTCAR')
        # structure = Structure.from_file('CONTCAR')

        pert_value = float(pert_folder)
        # get dummy index
        with open(f'{pert_folder}/nsc/POSCAR', 'rt') as f:
            poscar_lines = f.readlines()

        # elements list from POSCAR file where elements are grouped 
        elem_list = poscar_lines[5].split()

        # formula = atoms_final.get_chemical_formula(mode='metal')

        write_output = True


        try:
            num_ions = [int(item) for item in poscar_lines[5].split()]
        except ValueError:
            num_ions = [int(item) for item in poscar_lines[6].split()]

        # Create the resulting full elements list (where elements are not grouped) using a list comprehension
        elem_list_full = [elem for elem, count in zip(elem_list, num_ions) for _ in range(count)]

        dummy_index = dummy_position
        dummy_atom = elem_list_full[dummy_index]

        for elem, amount in zip(elem_list, num_ions):
            if elem == dummy_atom:
                if amount == 1:
                    break
                else:
                    raise ValueError('More than one dummy atom in the structure')
            else:
                # dummy_index += amount
                pass

        # get charges
        charges = []

        for step in ['nsc', 'sc']:
            outcar = Outcar(f'{pert_folder}/{step}/OUTCAR')

            calc = Vasp(directory=f'{pert_folder}/{step}')
            is_converged = calc.read_convergence()

            if 'f' in outcar.charge[dummy_index]:
                charges.append(outcar.charge[dummy_index]['f'])
            else:
                charges.append(outcar.charge[dummy_index]['d'])

            final_magmoms = [item['tot'] for item, ref in zip(outcar.magnetization, incar_bare['MAGMOM']) if ref != 0]
            for bare_magmom, final_magmom in zip(bare_magmoms, final_magmoms):
                if bare_magmom != 0:
                    # if the magnetic moment changes too much, do not write charges in output
                    if final_magmom / bare_magmom < 0.5 or final_magmom / bare_magmom > 2.0 or is_converged == False:
                        write_output = False

        if write_output:
            with open(os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold', filename), 'a') as f:
                f.write(f"{pert_value:5.2f}  {charges[0]}  {charges[1]}\n")
