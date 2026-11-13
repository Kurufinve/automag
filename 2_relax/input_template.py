# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# suffix for the structure (arbitrary string distinguishing structures with
# the same composition)
struct_suffix = ''

# configuration to relax; if omitted, the most stable configuration found by
# 2_coll is used
# configuration = 'afm1'

# spin values used if the converged moments from 2_coll are not available.
# OPTIONAL: when empty, pymatgen's DEFAULT_MAGMOMS are used.
spin_values = {
    'Fe': [5],
}

# specify the calculator: 'vasp', 'espresso', 'ase:<dotted.path>', ...
calculator = 'vasp'

# base calculation parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 830,
    'ediff': 1e-6,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 200,
    'kpts': 20,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1, -1],
    'ldauu': [5.17, 0, 0],
    'ldauj': [0, 0, 0],
    'ldauprint': 2,
    'lreal': 'Auto',
}

# relaxation-specific overrides
relax_params = {
    'ibrion': 2,
    'isif': 3,
    'nsw': 100,
    'ediffg': -0.01,
}

# copy the relaxed geometry into <workdir>/geometries as <name>_relaxed.vasp
update_structure = False

# specify if use fireworks database for managing calculations or not
use_fireworks = False

# if we do not intend to use fireworks we need to specify additional variables
jobheader = """#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=1-00:00:00
#SBATCH --partition=lenovo
#SBATCH --output=vasp-%j.out
#SBATCH --error=vasp-%j.error
"""

calculator_command = "export UCX_TLS=ud,sm,self; module load vasp/6.4.3; mpirun vasp_std"

environment_activate = "source /home/dpoletaev/miniconda3/bin/activate automag"

environment_deactivate = "conda deactivate"

parallel_over_configurations = True
