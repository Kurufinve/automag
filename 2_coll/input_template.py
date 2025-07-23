# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Zr2AlFe3.vasp'

# maximum supercell size for generating distinct magnetic configurations
supercell_size = 1

# choose the absolute values given to up and down spins
spin_values = {
    'Fe': [5],
}

# specify the calculator: 'vasp, 'qe', or 'fplo'

calculator = 'vasp'

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 2,
    'encut': 650,
    'ediff': 1e-6,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 200,
    'kpts': 50,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    # 'ldaul': [2, -1, -1],
    'ldaul': [-1, -1, 2],
    # 'ldauu': [5.17, 0, 0],
    'ldauu': [0, 0, 4.16],
    'ldauj': [0, 0, 0],
    'ldauprint': 2,
    'lreal': 'Auto',
}

# specify a cutoff for picking only high-spin configurations from output
# lower_cutoff = 1.7

# specify if use fireworks database for managing calculations or not
use_fireworks = False

# if we do not itended to use fireworks we need to specify additional variables
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

exec_command = "sbatch"
