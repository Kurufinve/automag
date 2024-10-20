# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# define the atom to introduce as dummy atom
dummy_atom = 'Zn'

# define the position of the dummy atom (from 0 to N_ATOMS-1)
dummy_position = 0

# specify the calculator: 'vasp, 'qe', or 'fplo'

calculator = 'vasp'

# define the perturbations in eV to apply to the dummy atom
perturbations = [-0.08, -0.05, -0.02, 0.02, 0.05, 0.08]

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 820,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'kpts': 20,
    'lmaxmix': 4,
    'nelm': 200,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Fe', 'Co']

# choose the magnetic configuration to use for U calculation (default FM-HS)
# configuration = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

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
