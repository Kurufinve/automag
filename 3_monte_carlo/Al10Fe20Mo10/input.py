# choose the configuration to use for Monte Carlo simulation
configuration = 'afm448'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 12.1

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = False

# choose the atomic types to be considered magnetic (default transition metals)
magnetic_atoms = ['Fe']

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Al10Fe20Mo10.vasp'

# name of the calculator with which singlepoint energies were computed
calculator = 'vasp'

struct_suffix = ''