# parameters for MAE calculation
# number of grid points for phi angle rotation (0-360 degrees)
Nph = 20

# number of grid points for theta angle rotation (0-180 degrees)
Nth = 10

# nomber of grid points for MAE curve (0-360 degrees)
N_MAE = 20

""" 
poscar_file and struct_suffix are needed for
finding the appropriate structure_folder in ../2_coll   
The initial structure is readed from file  
../2_coll/<structure_folder>/trials/setting{setting:03d}.vasp
"""
# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# suffix for the structure (arbitrary string for distinguishing structures with the same composition)
# struct_suffix = '_mp-1221736'
struct_suffix = ''

# choose the configuration to use for MAE calculation
configuration = 'afm1'

# try to find primitive cell from the given (super)cell taking into account the magnetic moments  
standardize_cell = True

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
#SBATCH --reservation=dpoletaev_28
""" 

calculator_command = "export UCX_TLS=ud,sm,self; module load vasp/6.4.3; mpirun vasp_std"

environment_activate = "source /home/dpoletaev/miniconda3/bin/activate automag"

environment_deactivate = "conda deactivate"




