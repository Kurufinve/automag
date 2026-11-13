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

import os,sys
import re
from scipy.optimize import curve_fit
import subprocess

cwd = os.getcwd()

from common.input_loader import apply_stage_input

apply_stage_input(globals())


import glob
import json
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
from common.sweep import (convergence_report, read_manifest,
                          read_parameters)


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
    """Return (volume, n_atoms, formula, Structure) of the calculated cell.

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
                    struct.composition.reduced_formula.replace(' ', ''), struct)
    struct = Structure.from_file(path_to_poscar)
    print('WARNING: the calculated cell could not be found; falling back to the '
          f'input geometry ({path_to_poscar}) for the MJ/m^3 conversion.')
    return (struct.volume, len(struct),
            struct.composition.reduced_formula.replace(' ', ''), struct)

# conversion factor: (eV/A^3) -> (MJ/m^3)
E2MJ = eV/(Ang**3)*(10.0**-6.0)

def magnetization(oszicar_path):
    """Total |M| in Bohr magnetons from the last ``mag=`` line of an OSZICAR.

    Handles both collinear (one component) and non-collinear (three
    components) files.
    """
    mags = []
    with open(oszicar_path) as f:
        for line in f:
            if 'mag=' in line:
                tail = line.partition('mag=')[2]
                values = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", tail)]
                if values:
                    mags = values
    if not mags:
        raise IOError(f'No "mag=" line found in {oszicar_path}')
    components = mags[:3] if len(mags) >= 3 else mags
    return float(np.sqrt(sum(c * c for c in components)))

def FunMAE(x, a, b, c):
    return a*((np.sin(x))**2)+b*((np.sin(x))**4)+c

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


standardized_structure_name = f'setting{setting:03d}_{configuration}_standardized.vasp'

if cell_info and cell_info.get('structure_file'):
    processed_structure_name = cell_info['structure_file']
else:
    processed_structure_name = standardized_structure_name

atoms = read(processed_structure_name)

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
compound_dir = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}{struct_suffix}")
state_dir = os.path.join(compound_dir, f'{calculator}/{configuration}')

def process_set(mae_dir, params, tag=''):
    """Analyse one MAE curve parameter set.

    ``mae_dir`` is the ``mae_*`` directory, ``params`` the calculation
    settings of the set and ``tag`` the suffix distinguishing the sets of a
    convergence study (empty for a single set).
    """
    U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
    J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
    encut = params['encut']
    kpts = params['kpts']
    # mae_dir = os.path.join(state_dir, 'mae')
    K = params['kpts']

    # the cell that VASP actually calculated (used for every normalization below)
    volume, n_atoms, calc_formula, calc_structure = get_calc_cell(mae_dir, cell_info, path_to_poscar)
    print(f'Calculated cell: {calc_formula}, {n_atoms} atoms, V = {volume:.4f} A^3')

    phi = np.linspace(0.0001, 2 * np.pi-0.0001, Nph+1)
    theta = np.linspace(0+0.0001, np.pi-0.00001, Nth+1)


    ### plotting the results of calculations

    try:
        os.mkdir(f'outputs_{configuration}')
    except:
        pass

    fig = plt.figure()
    fig.suptitle('MAE curve (Enthalpy) for U='+str(U)+', J='+str(J)+'. Atoms number = '+str(n_atoms)+', Volume = '+str(np.round(volume,1))+' $A^3$ ', fontsize=10)
    plt.rcParams['axes.grid'] = True
    plt.xlabel('$\\alpha$', size=15)
    plt.ylabel('Enthalpy $(MJ/m^3)$', size=15)
    MAE_vs_K = []

    E_alpha = []
    alpha_read = []
    MAE_E = []
    ref_oszicar = f"{mae_dir}/z/singlepoint/OSZICAR"
    if not os.path.isfile(ref_oszicar):
        raise IOError(f'{ref_oszicar} not found: run 1_submit.py and the VASP jobs first.')
    E_ref_MAE = float((Oszicar(ref_oszicar).all_energies[-1])[-2])
    folder_MAE_curve = f"{mae_dir}/z/singlepoint"
    # E_ref_MAE = float((Oszicar(mae_dir+"/z/singlepoint/OSZICAR").all_energies[-1])[-2])
    for alpha in  np.linspace(0, 2 * np.pi, N_MAE+1):
        folder = mae_dir+'/K_'+fmt_kpts(K)+'_RtMAE_'+str(np.round((alpha/np.pi)*180,2))+'/singlepoint'
        print('Folder is:'+folder)
        if os.path.isdir(folder) == False:
            print(folder, ' doesn\'t exist')
            continue
        elif os.path.isdir(folder) == True:
            print(folder," was found!")
            with open(folder+"/OUTCAR") as file:
                FinCond = 0
                for line in file:
                    if "General timing and accounting informations for this job:" in line:
                        FinCond = 1
                        #print('energies: ', Oszicar(folder+"/OSZICAR").all_energies)
                        E_alpha.append(float((Oszicar(folder+"/OSZICAR").all_energies[-1])[-2]))
                        alpha_read.append(alpha)
                        print('Energy for alpha='+str(np.round((alpha/np.pi)*180,2))+'  '+str(E_alpha[-1]))
                if FinCond == 0:
                    print('job for alpha='+str(np.round((alpha/np.pi)*180,2))+' did not finish successfully!')
        #print('TOTEN: ', TOTEN)

    if not E_alpha:
        raise IOError(f'No completed MAE curve calculations found under {mae_dir}. '
                      'Run 3_submit.py and the VASP jobs first.')

    outputs = open(f'MAE_curve_{configuration}{tag}_U{U:.1f}_J{J:.1f}_K{fmt_kpts(kpts)}_EN{encut}.txt','w')

    # per cell volume (MJ/m^3) and the MAE.py per-atom convention dE/(V*N)
    mae_curve_per_cell = (np.array(E_alpha)-E_ref_MAE)/volume*E2MJ
    mae_curve_per_atom = mae_curve_per_cell/n_atoms

    plt.plot(np.array(alpha_read)*180.0/np.pi, mae_curve_per_atom, '-o', label = "K = "+fmt_kpts(K))
    plt.legend()
    # if K == min(kpoints):
    np.save(f"./outputs_{configuration}/K_{fmt_kpts(K)}_MAE_curves_alpha_U"+str(U)+"_J"+str(J)+tag+".npy", np.array(alpha_read)*180.0/np.pi)
    np.save(f"./outputs_{configuration}/K_{fmt_kpts(K)}_MAE_curves_mae_U"+str(U)+"_J"+str(J)+tag+".npy", mae_curve_per_cell)
    np.save(f"./outputs_{configuration}/K_"+fmt_kpts(K)+"_raw_MAE_E_U"+str(U)+"_J"+str(J)+tag+".npy", np.array(E_alpha))
    np.save(f"./outputs_{configuration}/K_"+fmt_kpts(K)+"_MAE_E_U"+str(U)+"_J"+str(J)+tag+".npy", mae_curve_per_atom)
    np.save(f"./outputs_{configuration}/K_"+fmt_kpts(K)+"_MAE_alpha_U"+str(U)+"_J"+str(J)+tag+".npy", np.array(alpha_read)*180.0/np.pi)
    plt.savefig(f"./outputs_{configuration}/MAE_curves_U"+str(U)+"_J"+str(J)+tag+".png")
    MAE = (max(E_alpha)-min(E_alpha))
    K1 = (max(E_alpha)-min(E_alpha))
    MAE_per_cell = MAE/volume*E2MJ
    MAE_per_atom = MAE_per_cell/n_atoms
    MAE_vs_K.append(MAE_per_cell)

    M0_bohr = magnetization(folder_MAE_curve+'/OSZICAR')
    M0_A_per_m = M0_bohr*(Bh/(Ang**3))/volume      # A/m (per cell volume)
    M0_MA_per_m = M0_A_per_m*1e-6                  # MA/m
    BH_max = (1e-3)*0.25*mu0*M0_A_per_m**2         # kJ/m^3
    mu0_Ha = 2*(K1*eV)/(Bh*M0_bohr)                # T
    hardness = np.sqrt((K1/volume*(eV/(Ang**3)))/(mu0*M0_A_per_m**2))

    mycomp = calc_structure.composition
    symbols = calc_structure.symbol_set
    numbers = [mycomp.get_atomic_fraction(a)*mycomp.num_atoms for a in symbols]

    outputs.write('\n---------------------------------------------------------------------------------------\nPOSCAR information (calculated cell): \n')
    outputs.write('Volume = '+str(volume)+' \n')
    outputs.write('composition = '+str(symbols)+' \n')
    outputs.write('number of each element = '+str(numbers)+' \n')
    outputs.write('total number of atoms = '+str(n_atoms)+' \n')
    outputs.write('\n---------------------------------------------------------------------------------------\nVASP unit outputs: \n')
    outputs.write('MAE = '+str(MAE)+' eV\n')
    outputs.write('K1 = '+str(K1)+' eV\n')
    outputs.write('M0 = '+str(M0_bohr)+' Bohr Magneton\n')
    outputs.write('\n---------------------------------------------------------------------------------------\nSI unit outputs: \n')
    outputs.write("MAE = "+str(MAE_per_cell)+" (MJ/m^3, per cell)\n")
    outputs.write("MAE = "+str(MAE_per_atom)+" (MJ/m^3, per atom = dE/(V*N))\n")
    param, param_cov = curve_fit(FunMAE, np.array(alpha_read), (np.array(E_alpha)-min(np.array(E_alpha)))/volume*E2MJ)
    outputs.write("K1  = "+str(param[0])+" (MJ/m^3)\n")
    outputs.write("K2  = "+str(param[1])+" (MJ/m^3)\n")
    outputs.write("M0  = "+str(M0_MA_per_m)+" MA/m\n")
    outputs.write("(BH)max = "+str(BH_max)+" KJ/m^3\n")
    outputs.write("mu0xHa = "+str(mu0_Ha)+" T\n")
    outputs.write("Hardenss= "+str(hardness)+" \n")

    print("MAE = ", MAE_per_cell, " (MJ/m^3, per cell)")
    print("MAE = ", MAE_per_atom, " (MJ/m^3, per atom = dE/(V*N))")
    print("K1  = ", param[0], " (MJ/m^3)")
    print("K2  = ", param[1], " (MJ/m^3)")
    print("M0  = ", M0_MA_per_m, " MA/m")
    print("(BH)max = ", BH_max, " KJ/m^3")
    print("mu0xHa = ", mu0_Ha, " T")
    print("Hardenss= ", hardness)
    # subprocess.call('$('+'mkdir '+folder_MAE_curve+'_RtMAE_IsDone)', shell=True)

    outputs.write("MAE vs K = "+str(MAE_vs_K)+" (MJ/m^3, per cell)\n")
    # outputs.write("KPOINTS = "+str(np.array(nkpt_read)/readPOSCAR(path_to_poscar)[3])+" (number of Kpoints)/atom\n")
    outputs.close()

    results = {
        'tag': tag,
        'mae_dir': mae_dir,
        'params': params,
        'configuration': configuration,
        'MAE_eV': MAE,
        'MAE_MJ_m3_cell': MAE_per_cell,
        'MAE_MJ_m3_atom': MAE_per_atom,
        'K1_MJ_m3': float(param[0]),
        'K2_MJ_m3': float(param[1]),
        'M0_MA_per_m': float(M0_MA_per_m),
        'BH_max_kJ_m3': float(BH_max),
        'mu0_Ha_T': float(mu0_Ha),
        'hardness': float(hardness),
        'M0_bohr': float(M0_bohr),
        'volume': float(volume),
        'n_atoms': int(n_atoms),
        'alpha_deg': [float(v) for v in alpha_read],
        'E_alpha': [float(v) for v in E_alpha],
    }
    with open(f'MAE_curve_{configuration}{tag}.json', 'wt') as handle:
        json.dump(results, handle, indent=2, default=float)
    return results

# ---------------------------------------------------------------------------
# Analyse every parameter set of the MAE curve stage.
# ---------------------------------------------------------------------------
def discover_sets():
    """``(directory, parameters, tag, label)`` of every MAE parameter set."""
    directories = read_manifest(os.getcwd())
    if directories:
        entries = []
        for directory in directories:
            set_params, tag, label = read_parameters(directory)
            entries.append((directory, set_params or params, tag, label))
        return entries

    entries = []
    for summary in sorted(glob.glob(f'MAE_theta_phi_{configuration}*.json')):
        with open(summary, 'rt') as handle:
            meta = json.load(handle)
        directory = meta.get('mae_dir')
        if not directory or not os.path.isdir(directory):
            index = next(i for i, x in enumerate(params['ldaul']) if x > 0)
            U = np.round(float(params['ldauu'][index]), 1)
            J = np.round(float(params['ldauj'][index]), 1)
            directory = os.path.join(
                state_dir, f'mae_U{U:.1f}_J{J:.1f}'
                           f'_K{fmt_kpts(params["kpts"])}_EN{params["encut"]}')
        entries.append((directory, meta.get('params') or params,
                        meta.get('tag', ''),
                        meta.get('label', meta.get('tag') or 'single')))
    return entries


summaries = discover_sets()
if not summaries:
    raise IOError('No MAE theta/phi tables found: run 2_plot_results.py first.')

collected = []
for directory, set_params, tag, label in summaries:
    collected.append((tag, label, process_set(directory, set_params, tag)))

if len(collected) > 1:
    entries = [((tag or 'single', label), value['MAE_MJ_m3_cell'], value)
               for tag, label, value in collected]
    report = convergence_report(entries, quantity='MAE (MJ/m^3, per cell)')
    entries_k1 = [((tag or 'single', label), value['K1_MJ_m3'], value)
                  for tag, label, value in collected]
    report += '\n\n' + convergence_report(entries_k1,
                                           quantity='K1 (MJ/m^3)')
    print(report)
    with open(f'MAE_curve_convergence_{configuration}.txt', 'wt') as handle:
        handle.write(report + '\n')
