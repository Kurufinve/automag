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

# command for running espins
espins_run_command = 'module load intel2021/mkl/latest; /home/dpoletaev/soft/ESpinS/mc.x'

# name of the calculator with which singlepoint energies were computed
calculator = 'vasp'
# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
distances_between_neighbors = [2.898, 2.969, 3.362, 3.703, 3.981]

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
coupling_constants = [-3.40156121e-22, -3.72027152e-22, -5.33069829e-21, -3.56077373e-21,
 -1.12825138e-21]
