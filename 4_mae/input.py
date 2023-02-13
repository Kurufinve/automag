# name of the poscar file to use in the automag/geometries folder
poscar_file = 'Fe2O3-alpha_primitive.vasp'

# maximum supercell size for generating distinct magnetic configurations
supercell_size = 2

# choose the absolute values given to up and down spins
spin_values = {
    'Fe': [2],
}

# define the VASP parameters
params = {
    'xc': 'PBE',
    'setups': 'recommended',
    'prec': 'Accurate',
    'ncore': 4,
    'encut': 850,
    'ediff': 1e-6,
    'ismear': 1,
    'sigma': 0.1,
    'nelm': 200,
    'kpts': 40,
    'lmaxmix': 4,
    'lcharg': False,
    'lwave': False,
    'isym': 0,
    'ldau': True,
    'ldautype': 2,
    'ldaul': [2, -1],
    'ldauu': [2.6, 0],
    'ldauj': [0, 0],
    'ldauprint': 2,
}

# specify a cutoff for picking only high-spin configurations from output
# lower_cutoff = 1.7
