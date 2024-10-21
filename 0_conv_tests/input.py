# choose the desired mode: 'encut' or 'kgrid'
mode = 'kgrid'
#mode = 'encut'

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Zr2AlFe3.vasp'

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'istart': 0,
    'prec': 'Accurate',
    'ncore': 4,
    'ediff': 1e-6,
    'ismear': 1,
    'lcharg': False,
    'lwave': False,
    #'sigma': 0.2,
    #'kpts': 50,
    'encut': 850,
}

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Fe', 'Co']

# choose the magnetic configuration to use for convergence tests (default FM-HS)
# configuration = 6 * [4.0] + 6 * [-4.0] + 18 * [0.0]

# choose the trial values for ENCUT (default from 500 to 1000 eV at steps of 10 eV)
encut_values = range(500, 1050, 50)

# choose the trial values for SIGMA (default from 0.05 to 0.2 eV at steps of 0.05 eV)
sigma_values = [item / 100 for item in range(5, 25, 5)]

# choose the trial values for R_k (default from 20 to 100 Ang^-1 at steps of 10 Ang^-1)
kpts_values = range(20, 80, 10)

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