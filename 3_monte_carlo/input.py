# choose the configuration to use for Monte Carlo simulation
configuration = 'afm1'

# choose the cutoff radius in Angstrom for neighbor search
cutoff_radius = 3.8

# choose the size of the control group
control_group_size = 0.4

# append coupling constant before launching 2_write_vampire_ucf.py
append_coupling_constants = True

# choose the atomic types to be considered magnetic (default transition metals)
# magnetic_atoms = ['Mn']

# run command for ESpinS
espins_run_command = "module load compilers/intel; mpirun -np 1 /trinity/home/d.poletaev/soft/ESpinS/mc.x"

# run command for VAMPIRE
vampire_run_command = "module load apps/vampire/5.0.0; vampire-serial"

# # LINE ADDED BY THE SCRIPT 1_coupling_constants.py
# distances_between_neighbors = [2.898, 2.969, 3.362, 3.703, 3.981]

# # LINE ADDED BY THE SCRIPT 1_coupling_constants.py
# coupling_constants = [-1.57576668e-21, -1.47983623e-21, -9.05409734e-21, -5.70586508e-21,
#  -1.73453735e-21]

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
distances_between_neighbors = [2.898, 2.969, 3.362, 3.703]

# LINE ADDED BY THE SCRIPT 1_coupling_constants.py
coupling_constants = [-1.95765117e-21, -1.68784500e-21, -9.10171057e-21, -5.66138072e-21]
