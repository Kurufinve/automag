# choose the configuration to use for Monte Carlo simulation
# 'lowest' = the most stable configuration found by 2_coll
configuration = 'lowest'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 4.1

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = False

# Monte Carlo backend run by 5_run_mc.py: 'vampire' or 'espins'
mc_backend = 'vampire'
vampire_run_command = 'vampire'

# VAMPIRE input/material files (defaults are the bundled files in
# 3_monte_carlo/vampire_input)
# vampire_input_file = '3_monte_carlo/vampire_input/input'
# vampire_material_file = '3_monte_carlo/vampire_input/vamp.mat'
# vampire_output_file = 'output'

# ESpinS backend:
# espins_run_command = 'module load intel2021/mkl/latest; /path/to/mc.x'
# espins_output_file = 'result.dat'
# espins_output_converter = 'python convert_espins.py'

# optional timeout for the external Monte Carlo run (seconds)
# mc_timeout = 7 * 24 * 3600

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# name of the calculator with which singlepoint energies were computed
calculator = 'vasp'

struct_suffix = '_mp-1221736'