"""
automag.1_lin_response.2_plot_results
=====================================

Script which plots results of linear response U calculation.

.. codeauthor:: Michele Galasso <m.galasso@yandex.com>
"""

# default values for some input variables
use_fireworks = False
calculator = 'vasp'

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


import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

from scipy import stats

# getting formula
# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
# create ase.Atoms object
atoms = read(path_to_poscar)

filename_template = f"{atoms.get_chemical_formula(mode='metal')}_Ucalc_{calculator}"


# increase matplotlib pyplot font size
plt.rcParams.update({'font.size': 20})

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
data = np.loadtxt(os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}_charges_{calculator}.txt"))

perturbations = data[:, 0]
nscf = data[:, 1]
scf = data[:, 2]

plt.figure(figsize=(16, 9))

slope_nscf, intercept_nscf, r_value_nscf, p_value_nscf, std_err_nscf = stats.linregress(perturbations, nscf)
slope_scf, intercept_scf, r_value_scf, p_value_scf, std_err_scf = stats.linregress(perturbations, scf)

print(f'U = {(1/slope_scf) - (1/slope_nscf):4.2f}')
with open(f"{filename_template}.txt",'w') as f:
    f.write(f'U = {(1/slope_scf) - (1/slope_nscf):4.2f}')

plt.plot(perturbations, nscf, 'ro', label=f'NSCF (slope {slope_nscf:.4f})')
plt.plot(perturbations, scf, 'bo', label=f'SCF (slope {slope_scf:.4f})')

nscf_line = [value * slope_nscf + intercept_nscf for value in perturbations]
scf_line = [value * slope_scf + intercept_scf for value in perturbations]
plt.plot(perturbations, nscf_line, 'r-')
plt.plot(perturbations, scf_line, 'b-')

plt.xlabel('α (eV)')
plt.ylabel('d-electrons on first Fe site')
plt.legend()
# plt.show()
plt.savefig(f"{filename_template}.png", bbox_inches='tight')
# f"{structure.formula.replace(' ','')}_{calculator}

try:
    os.mkdir(f"{filename_template}")
except:
    pass

os.system(f'mv {filename_template}.* {filename_template}')