# choose the configuration to use for Monte Carlo simulation
configuration = 'afm1'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 4.1

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = False

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_conventional.vasp'

# name of the calculator with which singlepoint energies were computed
calculator = 'vasp'

struct_suffix = '_mp-1221736'