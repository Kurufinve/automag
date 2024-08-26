import re
import numpy as np
from ase.io import read
# from pymatgen.io.vasp import Poscar
# from pymatgen.io.vasp import Incar, Outcar

folder = 'singlepoint'
input_file = 'INCAR'


try:
    calculator = str(sys.argv[1])
except:
    calculator = 'vasp' # default calculator is vasp


if calculator == 'vasp':
    atoms = read(f'../{folder}/OUTCAR')


# magmoms_vasp = atoms_vasp.get_magnetic_moments()

try:
    magmoms = atoms.get_magnetic_moments()
except:
    magmoms = np.zeros(len(atoms.get_number_of_atoms))

print(np.array2string(magmoms, 10000))

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