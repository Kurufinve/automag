"""
automag.4_mae.2_plot_results
=============================

Script which plots results of finding Phi and Theta angles of z axis corresponind to min and max energy.

.. first codeauthor:: Michele Galasso <m.galasso@yandex.com>
.. second codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
use_fireworks = True
calculator = 'vasp'
create_mae_curve_input = True # if create input.py for MAE calculation in 4_mae folder

import os,sys

cwd = os.getcwd()

from common.input_loader import apply_stage_input

apply_stage_input(globals())


import json
import glob
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from copy import copy
from ase.io import read
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.vasp.outputs import Eigenval
from common.results import find_trials_dir, resolve_configuration
from common.sweep import convergence_report, read_manifest, read_parameters


Bh = 9.274009994*1e-24
mu0 = 4*np.pi*1e-7
eV = 1.602176634*1e-19
Ang = 1e-10
rad2deg = 57.2958

def better_sort(state):
    if state[0] == 'nm':
        return 'aa'
    elif state[0][:2] == 'fm':
        return 'aaa' + f'{int(state[0][2:]):3d}'
    else:
        return state[0][:3] + f'{int(state[0][3:]):3d}'

def readPOSCAR(path_to_poscar):
    struct = Structure.from_file(path_to_poscar)
    myV = struct.volume
    myAtoms_symbols = struct.symbol_set
    mycomp = struct.composition
    mycomp_atoms_num = mycomp.num_atoms
    myesp_num = []
    for atoms in myAtoms_symbols:
        myesp_num.append(mycomp.get_atomic_fraction(atoms)*mycomp_atoms_num)
    return myV, myAtoms_symbols, myesp_num, mycomp_atoms_num

def Zero_to_one(i):
    if abs(i) < 1e-10:
        return 1.0
    else:
        return i

def NormedCross(myvec1, myvec2):
    return(np.cross(myvec1, myvec2)/Zero_to_one(np.linalg.norm(np.cross(myvec1, myvec2))))

def fmt_kpts(kpts):
    """Compact label for a k-point specification (scalar or explicit mesh)."""
    if isinstance(kpts, (list, tuple)):
        return 'x'.join(str(k) for k in kpts)
    return str(kpts)

def load_cell_info():
    """Read the sidecar JSON written by 1_submit.py (may be absent)."""
    for setting_number in range(1, 10):
        name = f'setting{setting_number:03d}_{configuration}_processed.json'
        if os.path.isfile(name):
            with open(name, 'rt') as f:
                return json.load(f)
    return None

def get_calc_cell(mae_dir, cell_info, path_to_poscar):
    """Return (volume, n_atoms, formula) of the cell that was calculated.

    All MJ/m^3 and per-atom conversions must use this cell. Using the input
    geometry instead (which may be a different supercell) silently scales the
    MAE by the volume/atom ratio of the two cells.
    """
    candidates = [
        os.path.join(mae_dir, 'z', 'singlepoint', 'POSCAR'),
        os.path.join(mae_dir, 'z', 'POSCAR'),
    ]
    if cell_info and cell_info.get('structure_file'):
        candidates.append(cell_info['structure_file'])
    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            struct = Structure.from_file(candidate)
            return (struct.volume, len(struct),
                    struct.composition.reduced_formula.replace(' ', ''))
    struct = Structure.from_file(path_to_poscar)
    print('WARNING: the calculated cell could not be found; falling back to the '
          f'input geometry ({path_to_poscar}) for the MJ/m^3 conversion.')
    return (struct.volume, len(struct),
            struct.composition.reduced_formula.replace(' ', ''))

# conversion factor: (eV/A^3) -> (MJ/m^3)
E2MJ = eV/(Ang**3)*(10.0**-6.0)

# take care of the case when lower_cutoff has not been specified
if 'lower_cutoff' not in globals():
    lower_cutoff = 0

# setting default tolerance in determining if magmoms can be considered the same
if 'tolerance' not in globals():
    tolerance = 0.08

# full path to poscar file
rel_path_to_poscar = '/geometries/' + poscar_file
path_to_automag = os.environ.get('AUTOMAG_PATH')
path_to_poscar = path_to_automag + rel_path_to_poscar
input_structure = Structure.from_file(path_to_poscar)

formula = input_structure.formula.replace(' ','')
path_to_automag = os.environ.get('AUTOMAG_PATH')
# The collinear results may have been produced with a different calculator
# than the one used for the MAE calculations (mixed-calculator pipeline).
coll_calculator = globals().get('coll_calculator') or calculator
configuration = resolve_configuration(
    configuration, path_to_automag, formula, struct_suffix, coll_calculator)
path_to_trials = str(find_trials_dir(
    path_to_automag, formula, struct_suffix, coll_calculator))
print(f'Path to trials: {path_to_trials}')

cell_info = load_cell_info()


# finding configuration and settings

setting = 1
magmom = None
while os.path.isfile(f'{path_to_trials}/configurations{setting:03d}.txt'):
    with open(f'{path_to_trials}/configurations{setting:03d}.txt', 'rt') as f:
        for line in f:
            values = line.split()
            if values[0] == configuration:
                magmom = [int(item) for item in values[1:]]
                break
        if magmom is not None:
            break
    setting += 1

if cell_info is not None and cell_info.get('setting'):
    setting = int(cell_info['setting'])

if magmom is None and cell_info is None:
    raise IOError(f'configuration {configuration} not found and no processed-cell sidecar file is present.')
if magmom is None:
    print(f'No 2_coll configuration data for {configuration}: using the cell from the sidecar file.')

# # # get materials from magmom
# # materials = [0 if item >= 0 else 1 for item in magmom]

# # create a pymatgen Structure object
# pmg_structure = Structure.from_file(f'{path_to_trials}/setting{setting:03d}.vasp')
# ase_structure = pmg_structure.to_ase_atoms()

# # assigning magnetic moments to the atoms

# calc_results = f"{ase_structure.get_chemical_formula(mode='metal')}{struct_suffix}_singlepoint_{calculator}.txt"

# path_to_calc_results = path_to_automag + '/CalcFold/' + calc_results
# full_magmom = None
# print(f'Reading full magmoms for configuration {configuration} from file {path_to_calc_results}')
# with open(path_to_calc_results, 'rt') as f:
#     lines = f.readlines()
#     i = 0
#     for line in lines:
#         values = line.split()
#         if values[0] == configuration:
#             print(f'Found configuration: {values[0]}')
#             magmoms_str = lines[i+1]
#             break
#         i += 1
#     full_magmom_str = magmoms_str.split('final_magmoms=')[1].strip().strip('[]')
#     number_magmoms = full_magmom_str.split()
#     full_magmom = [float(num) for num in number_magmoms]

# if full_magmom is None:
#     raise IOError(f'Magmoms for configuration {configuration} not found.')


# print(f'Magnetic moments of all atoms for configuration {configuration}:')
# print(full_magmom)
# pmg_structure.add_site_property("magmom", full_magmom)

# standardized_structure = pmg_structure.get_primitive_structure(tolerance=0.2, use_site_props=True)



standardized_structure_name = f'setting{setting:03d}_{configuration}_standardized.vasp'

if cell_info and cell_info.get('structure_file'):
    processed_structure_name = cell_info['structure_file']
else:
    processed_structure_name = standardized_structure_name

atoms = read(processed_structure_name)

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
compound_dir = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}{struct_suffix}")
state_dir = os.path.join(compound_dir, f'{calculator}/{configuration}')

def process_mae_dir(mae_dir, params, tag=''):
    """Analyse one MAE parameter set and write its table.

    ``mae_dir`` is the ``mae_*`` directory of the set, ``params`` its
    calculation parameters and ``tag`` the suffix that distinguishes the
    sets of a convergence study (empty for a single set, so the historical
    file names are kept).
    """
    U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
    J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
    encut = params['encut']
    kpts = params['kpts']
    # mae_dir = os.path.join(state_dir, 'mae')

    # the cell that VASP actually calculated (used for every normalization below)
    volume, n_atoms, calc_formula = get_calc_cell(mae_dir, cell_info, path_to_poscar)
    print(f'Calculated cell: {calc_formula}, {n_atoms} atoms, V = {volume:.4f} A^3')


    phi = np.linspace(0.0001, 2 * np.pi-0.0001, Nph+1)
    theta = np.linspace(0+0.0001, np.pi-0.00001, Nth+1)


    ### plotting the results of calculations

    try: 
        os.mkdir(f'outputs_{configuration}')
    except:
        pass

    # U = params['ldauu'][params['ldaul'].index(2)]
    # U = params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]
    # J = 0

    fig = plt.figure()
    fig.suptitle(r'Enthalpy as a function of $\theta$ for different $\phi$ for U='+str(U)+', J='+str(J))
    plt.rcParams['axes.grid'] = True
    plt.xlabel(r'$\theta$', size=15)
    plt.ylabel('Enthapy $(MJ/m^3)$', size=15)
    all_phi = []
    all_theta = []
    all_E = []
    colors = [cm.rainbow(i) for i in np.linspace(0, 1, Nph+1)]
    ref_oszicar = f"{mae_dir}/z/singlepoint/OSZICAR"
    if not os.path.isfile(ref_oszicar):
        raise IOError(f'{ref_oszicar} not found: run 1_submit.py and the VASP jobs first.')
    E_ref = float((Oszicar(ref_oszicar).all_energies[-1])[-2])
    for phin, c in zip(phi, colors):
        TOTEN = []
        theta_read = []
        for thetan in theta:
            folder = mae_dir+'/PhTh_'+str(np.round((phin/np.pi)*180,2))+'_'+str(np.round((thetan/np.pi)*180,2))+'/singlepoint'
            print('Folder is:'+folder)
            if os.path.isdir(folder) == False:
                print(folder, ' doesn\'t exist')
                continue
            elif os.path.isdir(folder) == True:
                print("folder ", folder," was found!")
                with open(folder+"/OUTCAR") as file:
                    FinCond = 0
                    for line in file:
                        if "General timing and accounting informations for this job:" in line:
                            FinCond = 1
                            #print('energies: ', Oszicar("./"+folder+"/OSZICAR").all_energies)
                            TOTEN.append(float((Oszicar(folder+"/OSZICAR").all_energies[-1])[-2]))
                            theta_read.append(thetan)
                            all_theta.append(thetan)
                            all_phi.append(phin)
                            all_E.append(float((Oszicar(folder+"/OSZICAR").all_energies[-1])[-2]))
                            print('Energy for Phi='+str(np.round(phin,2))+' and theta='+str(np.round(thetan,2))+'  '+str(TOTEN[-1]))
                    if FinCond == 0:
                        print('job for Phi='+str(np.round(phin,2))+' and theta='+str(np.round(thetan,2))+' is not finished and will be skipped!')
        #print('TOTEN: ', TOTEN)
        # plt.plot(np.array(theta_read)*180.0/np.pi, (np.array(TOTEN)-E_ref)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), '-o',color=c, label=r'$\phi$='+str(np.round(phin,3)))
        plt.plot(np.array(theta_read)*180.0/np.pi, (np.array(TOTEN)-E_ref)/volume*E2MJ/n_atoms, '-o',color=c, label=r'$\phi$='+str(np.round(phin*180.0/np.pi))+r'$^{\circ}$')
        np.save(f"./outputs_{configuration}/TOTEN_"+str(np.round(phin,3))+f"_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}{tag}.npy", (np.array(TOTEN)-E_ref)/volume*E2MJ/n_atoms)
        np.save(f"./outputs_{configuration}/theta_"+str(np.round(phin,3))+f"_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}{tag}.npy", np.array(theta_read)*180.0/np.pi)
        #plt.plot(theta, TOTEN)
    plt.legend(prop={'size': 8})
    plt.savefig(f"./outputs_{configuration}/EvsTheta_curves_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}{tag}.png")

    if not all_E:
        raise IOError(f'No completed PhTh calculations found under {mae_dir}. '
                      'Run 1_submit.py and the VASP jobs first.')

    table_name = f'MAE_theta_phi_{configuration}{tag}_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}.txt'
    outputs = open(table_name,'w')

    MAE = (max(all_E)-min(all_E))
    K1 = (max(all_E)-min(all_E))

    for e1, p1, t1 in zip(all_E, all_phi, all_theta):
        if e1==min(all_E):
            phi_min = p1; theta_min = t1
        if e1 == max(all_E):
            phi_max = p1; theta_max = t1
    outputs.write('phi_min: '+str(np.round(phi_min*180/np.pi,2))+',   theta_min: '+str(np.round(theta_min*180/np.pi,2))+'\n')
    outputs.write('phi_max: '+str(np.round(phi_max*180/np.pi,2))+',   theta_max: '+str(np.round(theta_max*180/np.pi,2))+'\n')
    print('phi_min: ', np.round(phi_min*180/np.pi,2), 'theta_min: ', np.round(theta_min*180/np.pi,2))
    print('phi_max: ', np.round(phi_max*180/np.pi,2), 'theta_max: ', np.round(theta_max*180/np.pi,2))
    min_vec = [np.sin(theta_min) * np.cos(phi_min), np.sin(theta_min) * np.sin(phi_min), np.cos(theta_min)]
    max_vec = [np.sin(theta_max) * np.cos(phi_max), np.sin(theta_max) * np.sin(phi_max), np.cos(theta_max)]
    print('max_vec: ', np.round(max_vec,2))
    print('min_vec: ', np.round(min_vec,2))
    outputs.write('max_vec: '+str([float(v) for v in np.round(max_vec,2)])+'\n')
    outputs.write('min_vec: '+str([float(v) for v in np.round(min_vec,2)])+'\n')

    MAE_x = np.array(min_vec[:])
    # MAE_z is the normal of the plane spanned by the easy and hard axes.
    # If they are (numerically) parallel - which happens when both sit on the
    # theta~0 pole - fall back to a deterministic perpendicular axis.  The old
    # random fallback made the MAE curve irreproducible.
    cross_vec = np.array(np.cross(min_vec, max_vec), dtype=float)
    if np.linalg.norm(cross_vec) > 0.0001:
        MAE_z = cross_vec/np.linalg.norm(cross_vec)
    else:
        reference = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(reference, MAE_x)) > 0.9:
            reference = np.array([0.0, 1.0, 0.0])
        MAE_z = np.cross(MAE_x, reference)
        MAE_z = MAE_z/np.linalg.norm(MAE_z)
        print('WARNING: easy and hard axis are (anti)parallel; using a deterministic '
              f'perpendicular axis {np.round(MAE_z, 3)} for the MAE curve.')
    MAE_y = np.array(np.cross(MAE_z, MAE_x))
    MAE_y = MAE_y/np.linalg.norm(MAE_y)
    outputs.write('MAE_x: '+str([float(v) for v in np.round(MAE_x,3)])+', norm='+str(np.round(np.linalg.norm(MAE_x),3))+', easy axis\n')
    outputs.write('MAE_z: '+str([float(v) for v in np.round(MAE_z,3)])+', norm='+str(np.round(np.linalg.norm(MAE_z),3))+'\n')
    outputs.write('MAE_y: '+str([float(v) for v in np.round(MAE_y,3)])+', norm='+str(np.round(np.linalg.norm(MAE_y),3))+', hard axis\n')
    outputs.write(f'MAE: {MAE*1000:.6f} meV, {MAE/volume*E2MJ:.6f} MJ/m^3 (per cell), '
                  f'{MAE/volume*E2MJ/n_atoms:.6f} MJ/m^3 (per atom)\n')
    outputs.write(f'calculated cell: {calc_formula}, {n_atoms} atoms, V = {volume:.4f} A^3\n')
    print('MAE_x: ', np.round(MAE_x,3), ', norm', np.round(np.linalg.norm(MAE_x),3), '(easy axis)')
    print('MAE_z: ', np.round(MAE_z,3), ', norm', np.round(np.linalg.norm(MAE_z),3))
    print('MAE_y: ', np.round(MAE_y,3), ', norm', np.round(np.linalg.norm(MAE_y),3), '(hard axis)')
    print(f'MAE = {MAE*1000:.6f} meV = {MAE/volume*E2MJ:.6f} MJ/m^3 (per cell) = '
          f'{MAE/volume*E2MJ/n_atoms:.6f} MJ/m^3 (per atom)')


    outputs.close()

    with open(f'MAE_theta_phi_{configuration}{tag}.json', 'wt') as handle:
        json.dump({
            'tag': tag,
            'table': table_name,
            'mae_dir': mae_dir,
            'params': params,
            'configuration': configuration,
            'MAE_eV': MAE,
            'MAE_MJ_m3_cell': MAE / volume * E2MJ,
            'MAE_MJ_m3_atom': MAE / volume * E2MJ / n_atoms,
            'MAE_x': [float(v) for v in MAE_x],
            'MAE_y': [float(v) for v in MAE_y],
            'MAE_z': [float(v) for v in MAE_z],
            'phi_min': float(phi_min), 'theta_min': float(theta_min),
            'phi_max': float(phi_max), 'theta_max': float(theta_max),
            'volume': float(volume), 'n_atoms': int(n_atoms),
        }, handle, indent=2, default=float)
    return MAE, MAE / volume * E2MJ

# ---------------------------------------------------------------------------
# Analyse every parameter set of the MAE stage.
#
# One set -> the historical file names (no suffix).  Several sets (a
# convergence study, see common/sweep.py) -> one table per set plus a summary
# that compares them.
# ---------------------------------------------------------------------------
def discover_sets():
    directories = read_manifest(os.getcwd())
    if directories:
        return directories
    return [path for path in sorted(glob.glob(os.path.join(state_dir, 'mae_*')))
            if os.path.isdir(path)]


def set_metadata(directory):
    set_params, tag, label = read_parameters(directory)
    return set_params or params, tag, label


directories = discover_sets()
if not directories:
    raise IOError(f'No MAE calculation directories found under {state_dir}: '
                  'run 1_submit.py and the calculations first.')

results = []
for directory in directories:
    set_params, tag, label = set_metadata(directory)
    mae_e, mae_per_cell = process_mae_dir(directory, set_params, tag)
    results.append(((tag or 'single', label), mae_per_cell, {}))

if len(results) > 1:
    report = convergence_report(results, quantity='MAE (MJ/m^3, per cell)')
    print(report)
    with open(f'MAE_convergence_{configuration}.txt', 'wt') as handle:
        handle.write(report + '\n')
