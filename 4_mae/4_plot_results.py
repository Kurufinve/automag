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

try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    try:
        exec(f"from {input_file.split('.')[0]} import *")
    except:
        from input import *
except IndexError:
    print(f'Using the input.py file from folder: {cwd}')
    from input import *


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

def magnetization(fileOutCar):
    Mag = 'Nan'
    with open(fileOutCar) as file:#open("./z/OUTCAR") as file:
        for line in file:
            if "General timing and accounting informations for this job:" in line:  
                with open(f"{mae_dir}/z/singlepoint/OSZICAR") as file:
                    for line in file:
                        if 'mag=' in line:
                            head, sep, tail = line.partition('mag=')
                            Mag = [float(j) for j in re.findall(r"[-+]?\d*\.?\d+|\d+", tail)]
                M0 = np.sqrt((Mag[0]**2)+(Mag[1]**2)+(Mag[2]**2))
                #print('total magnetization:', np.sqrt((Mag[0]**2)+(Mag[1]**2)+(Mag[2]**2)))
    return M0

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
# path to the folder with results from collinear calculations
path_to_coll = path_to_automag + '/2_coll/' + f'{formula}{struct_suffix}/{formula}_{calculator}/'
if not os.path.isdir(path_to_coll):
    path_to_coll = path_to_automag + '/2_coll/' + f'{formula}_{calculator}/'
    if not os.path.isdir(path_to_coll):
        path_to_coll = path_to_automag + '/2_coll/'

print(f'Path to the results from collinear calculations: {path_to_coll}')

path_to_trials = path_to_coll + f'trials_{formula}/'


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

if magmom is None:
    raise IOError(f'configuration {configuration} not found.')


standardized_structure_name = f'setting{setting:03d}_{configuration}_standardized.vasp'

atoms = read(standardized_structure_name)

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
compound_dir = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}{struct_suffix}")
state_dir = os.path.join(compound_dir, f'{calculator}/{configuration}')

U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
encut = params['encut']
kpts = params['kpts']
mae_dir = os.path.join(state_dir, f'mae_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}')
# mae_dir = os.path.join(state_dir, 'mae')
K = params['kpts']

phi = np.linspace(0.0001, 2 * np.pi-0.0001, Nph+1)
theta = np.linspace(0+0.0001, np.pi-0.00001, Nth+1)


### plotting the results of calculations

try: 
    os.mkdir(f'outputs_{configuration}')
except:
    pass
outputs = open(f'MAE_curve_{configuration}_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.txt','w')


fig = plt.figure()
fig.suptitle('MAE curve (Enthalpy) for U='+str(U)+', J='+str(J)+'. Atoms number = '+str(readPOSCAR(path_to_poscar)[3])+', Volume = '+str(np.round(readPOSCAR(path_to_poscar)[0],1))+' $A^3$ ', fontsize=10)
plt.rcParams['axes.grid'] = True
plt.xlabel('$\\alpha$', size=15)
plt.ylabel('Enthalpy $(MJ/m^3)$', size=15)
MAE_vs_K = []

E_alpha = []
alpha_read = []
MAE_E = []
E_ref_MAE = float((Oszicar(f"{mae_dir}/z/singlepoint/OSZICAR").all_energies[-1])[-2])
folder_MAE_curve = f"{mae_dir}/z/singlepoint"
# E_ref_MAE = float((Oszicar(mae_dir+"/z/singlepoint/OSZICAR").all_energies[-1])[-2])
for alpha in  np.linspace(0, 2 * np.pi, N_MAE+1):
    folder = mae_dir+'/K_'+str(K)+'_RtMAE_'+str(np.round((alpha/np.pi)*180,2))+'/singlepoint'
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
plt.plot(np.array(alpha_read)*180.0/np.pi, (np.array(E_alpha)-E_ref_MAE)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), '-o', label = "K = "+str(K))
plt.legend()
# if K == min(kpoints):
np.save(f"./outputs_{configuration}/K_{K}_MAE_curves_alpha_U"+str(U)+"_J"+str(J)+".npy", np.array(alpha_read)*180.0/np.pi)
np.save(f"./outputs_{configuration}/K_{K}_MAE_curves_mae_U"+str(U)+"_J"+str(J)+".npy", (np.array(E_alpha)-E_ref_MAE)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0))
np.save(f"./outputs_{configuration}/K_"+str(K)+"_raw_MAE_E_U"+str(U)+"_J"+str(J)+".npy", np.array(E_alpha))
np.save(f"./outputs_{configuration}/K_"+str(K)+"_MAE_E_U"+str(U)+"_J"+str(J)+".npy", (np.array(E_alpha)-E_ref_MAE)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]))
np.save(f"./outputs_{configuration}/K_"+str(K)+"_MAE_alpha_U"+str(U)+"_J"+str(J)+".npy", np.array(alpha_read)*180.0/np.pi)
plt.savefig(f"./outputs_{configuration}/MAE_curves_U"+str(U)+"_J"+str(J)+".png")
MAE = (max(E_alpha)-min(E_alpha))
K1 = (max(E_alpha)-min(E_alpha))
MAE_vs_K.append(MAE*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0))

outputs.write('\n---------------------------------------------------------------------------------------\nPOSCAR information: \n')
outputs.write('Volume = '+str(readPOSCAR(path_to_poscar)[0])+' \n')
outputs.write('composition = '+str(readPOSCAR(path_to_poscar)[1])+' \n')
outputs.write('number of each element = '+str(readPOSCAR(path_to_poscar)[2])+' \n')
outputs.write('total number of atoms = '+str(readPOSCAR(path_to_poscar)[3])+' \n')
outputs.write('\n---------------------------------------------------------------------------------------\nVASP unit outputs: \n')
outputs.write('MAE = '+str(MAE)+' eV\n')
outputs.write('K1 = '+str(K1)+' eV\n')
outputs.write('M0 = '+str(magnetization(folder_MAE_curve+'/OUTCAR'))+' Bohr Magneton\n')
outputs.write('\n---------------------------------------------------------------------------------------\nSI unit outputs (MAE, K1 and M0 are calculated per atom): \n')
outputs.write("MAE = "+str(MAE*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0))+" (MJ/m^3)\n")
param, param_cov = curve_fit(FunMAE, np.array(alpha_read), (np.array(E_alpha)-min(np.array(E_alpha)))*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0))
outputs.write("K1  = "+str(param[0])+" (MJ/m^3)\n")
outputs.write("K2  = "+str(param[1])+" (MJ/m^3)\n")
outputs.write("M0  = "+str(magnetization(folder_MAE_curve+'/OUTCAR')*(1e-6)*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0])*(1))+" MA/m\n")
outputs.write("(BH)max = "+str((1e-3)*0.25*mu0*((magnetization(folder_MAE_curve+'/OUTCAR')*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0]))**2))+" KJ/m^3\n")
outputs.write("mu0xHa = "+str(2*(K1*eV)/(Bh*magnetization(folder_MAE_curve+'/OUTCAR')))+" T\n")
outputs.write("Hardenss= "+ str(np.sqrt((K1*((1.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0))/(mu0*((magnetization(folder_MAE_curve+'/OUTCAR')*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0]))**2))))+" \n")

#print(max((np.array(E_alpha)-E_ref_MAE)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]))-min((np.array(E_alpha)-E_ref_MAE)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3])))
print("MAE = ", MAE*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), " (MJ/m^3)")
print("K1  = ", K1*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), " (MJ/m^3)")
print("M0  = ", magnetization(folder_MAE_curve+'/OUTCAR')*(1e-6)*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0])*(1/readPOSCAR(path_to_poscar)[3]), " MA/m")
print("(BH)max = ", (1e-3)*0.25*mu0*((magnetization(folder_MAE_curve+'/OUTCAR')*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0]))**2), " KJ/m^3")
print("mu0xHa = ", 2*(K1*eV)/(Bh*magnetization(folder_MAE_curve+'/OUTCAR')), " T")
print("Hardenss= ", np.sqrt((K1*((1.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]))/(mu0*((magnetization(folder_MAE_curve+'/OUTCAR')*(Bh/(Ang**3))*(1/readPOSCAR(path_to_poscar)[0]))**2))))
# subprocess.call('$('+'mkdir '+folder_MAE_curve+'_RtMAE_IsDone)', shell=True)

outputs.write("MAE vs K = "+str(MAE_vs_K)+" (MJ/m^3)\n")
# outputs.write("KPOINTS = "+str(np.array(nkpt_read)/readPOSCAR(path_to_poscar)[3])+" (number of Kpoints)/atom\n")
outputs.close()
# fig = plt.figure()
# plt.rcParams['axes.grid'] = True
# fig.suptitle('MAE as a function of number of K-points for U='+str(U)+', J='+str(J))
# plt.xlabel(r'(number of Kpoints)', size=15)
# plt.ylabel(r'MAE $(MJ/m^3)$', size=15)
# np.save(f"./outputs_{configuration}/MAE_vs_Kpoints_curve_U"+str(U)+"_J"+str(J)+"_x.npy", np.array(nkpt_read))
# np.save(f"./outputs_{configuration}/MAE_vs_Kpoints_curve_U"+str(U)+"_J"+str(J)+"_y.npy", MAE_vs_K)
# plt.plot(np.array(nkpt_read), MAE_vs_K, '-o',label = "U= "+str(U)+", J="+str(J))
# plt.savefig(f"./outputs_{configuration}/MAE_vs_Kpoints_curve_U"+str(U)+"_J"+str(J)+".png")

