# choose the configuration to use for Monte Carlo simulation
configuration = 'afm448'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 10.1

# choose the size of the control group
control_group_size = 0.2

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = True

# choose the atomic types to be considered magnetic (default transition metals)
magnetic_atoms = ['Fe']

# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Al10Fe20Mo10.vasp'

# name of the calculator with which singlepoint energies were computed
calculator = 'vasp'

struct_suffix = ''
# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
distances_between_neighbors = [2.93, 4.144, 5.076, 5.861, 6.552, 7.178, 8.288, 8.791, 9.267, 9.719]

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
coupling_constants = [ 4.89176601e-22, -2.22882271e-21, -1.52146090e-21,  2.57554824e-21,
  6.36806047e-22, -5.76565977e-22,  1.14336801e-21,  4.61271149e-23,
  8.22520718e-22, -4.80218439e-22]
