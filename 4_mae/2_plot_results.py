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

atoms = read(standardized_structure_name)

calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
compound_dir = os.path.join(calcfold, f"{atoms.get_chemical_formula(mode='metal')}{struct_suffix}")
state_dir = os.path.join(compound_dir, f'{calculator}/{configuration}')

mae_dir = os.path.join(state_dir, 'mae')


phi = np.linspace(0.0001, 2 * np.pi-0.0001, Nph+1)
theta = np.linspace(0+0.0001, np.pi-0.00001, Nth+1)


### plotting the results of calculations

try: 
    os.mkdir('outputs')
except:
    pass

# U = params['ldauu'][params['ldaul'].index(2)]
U = params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]
J = 0

fig = plt.figure()
fig.suptitle(r'Enthalpy as a function of $\theta$ for different $\phi$ for U='+str(U)+', J='+str(J))
plt.rcParams['axes.grid'] = True
plt.xlabel(r'$\theta$', size=15)
plt.ylabel('Enthapy $(MJ/m^3)$', size=15)
all_phi = []
all_theta = []
all_E = []
colors = [cm.rainbow(i) for i in np.linspace(0, 1, Nph+1)]
E_ref = float((Oszicar(f"{mae_dir}/z/singlepoint/OSZICAR").all_energies[-1])[-2])
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
                    all_E.append(np.log(-2)) #to ignore in plotting
                    print('job for Phi='+str(np.round(phin,2))+' and theta='+str(np.round(thetan,2))+' is not finished!')
    #print('TOTEN: ', TOTEN)
    # plt.plot(np.array(theta_read)*180.0/np.pi, (np.array(TOTEN)-E_ref)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), '-o',color=c, label=r'$\phi$='+str(np.round(phin,3)))
    plt.plot(np.array(theta_read)*180.0/np.pi, (np.array(TOTEN)-E_ref)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]), '-o',color=c, label=r'$\phi$='+str(np.round(phin*rad2deg))+r'$^{\circ}$')
    np.save("./outputs/TOTEN_"+str(np.round(phin,3))+".npy", (np.array(TOTEN)-E_ref)*((10.0**-6.0)/readPOSCAR(path_to_poscar)[0])*(eV/(Ang**3))*(1.0/readPOSCAR(path_to_poscar)[3]))
    np.save("./outputs/theta_"+str(np.round(phin,3))+".npy", np.array(theta_read)*180.0/np.pi)
    #plt.plot(theta, TOTEN)
plt.legend(prop={'size': 8})
plt.savefig("./outputs/EvsTheta_curves_U"+str(U)+"_J"+str(J)+".png")



outputs = open('MAE_theta_phi.txt','w')

MAE = (max(all_E)-min(all_E))
K1 = (max(all_E)-min(all_E))

for e1, p1,t1 in zip(all_E, all_phi, all_theta):
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
outputs.write('max_vec: '+str(list(np.round(max_vec,2)))+'\n')
outputs.write('min_vec: '+str(list(np.round(min_vec,2)))+'\n')

MAE_x = np.array(min_vec[:])
if max(abs(np.array(np.cross(min_vec, max_vec))))>0.0001:
    MAE_z = np.array(NormedCross(min_vec, max_vec))
else:
    if os.path.isfile('./outputs/MAE_Ran_vec.npy') == False:
        Ran_vec = np.random.rand(3)
        Ran_vec = Ran_vec/np.linalg.norm(Ran_vec)
        MAE_z = np.array(np.cross(min_vec, Ran_vec))/np.linalg.norm(np.cross(min_vec, Ran_vec))
        np.save('./outputs/MAE_Ran_vec.npy', MAE_z)
    if os.path.isfile('./outputs/MAE_Ran_vec.npy') == True:
        MAE_z = np.load('./outputs/MAE_Ran_vec.npy')
MAE_y = np.array(np.cross(MAE_z, MAE_x))/np.linalg.norm(np.cross(MAE_z, MAE_x))
outputs.write('MAE_x: '+str(list(np.round(MAE_x,3)))+', norm='+str(np.round(np.linalg.norm(MAE_x)))+'\n')
outputs.write('MAE_z: '+str(list(np.round(MAE_z,3)))+', norm='+str(np.round(np.linalg.norm(MAE_y)))+'\n')
outputs.write('MAE_y: '+str(list(np.round(MAE_y,3)))+', norm='+str(np.round(np.linalg.norm(MAE_z)))+'\n')
print('MAE_x: ', np.round(MAE_x,3), ', norm', np.round(np.linalg.norm(MAE_x)))
print('MAE_z: ', np.round(MAE_z,3), ', norm', np.round(np.linalg.norm(MAE_y)))
print('MAE_y: ', np.round(MAE_y,3), ', norm', np.round(np.linalg.norm(MAE_z)))


outputs.close()




























# # extract the results
# red = []
# not_converged = []
# all_final_magmoms = []
# for line, maginfo in zip(lines, maginfos):
#     values = line.split()
#     init_state = values[0]

#     if values[-4].split('=')[1] == 'NONCONVERGED':
#         not_converged.append(init_state)
#         del data[init_state]
#     else:
#         initial, final = maginfo.split('final_magmoms=')
#         initial = initial[initial.index('[') + 1:initial.index(']')]
#         final = final[final.index('[') + 1:final.index(']')]
#         initial = np.array(initial.split(), dtype=float)
#         final = np.array(final.split(), dtype=float)

#         all_final_magmoms.extend(final.tolist())

#         # exclude low-spin configurations
#         flag = False
#         mask = np.nonzero(initial)
#         if np.all(np.abs(final[mask]) > lower_cutoff) or init_state == 'nm':
#             flag = True

#         magnification = len(initial) // sum(multiplicities)
#         current_multiplicities = magnification * multiplicities

#         start = 0
#         kept = True
#         for multiplicity in current_multiplicities:
#             w_initial = initial[start:start + multiplicity]
#             w_final = final[start:start + multiplicity]
#             start += multiplicity

#             if not np.any(w_initial):
#                 if np.any(np.around(w_final)):
#                     kept = False
#             else:
#                 prod = np.multiply(np.sign(w_initial), w_final)
#                 # if prod.max() - prod.min() > 0.08:
#                 if prod.max() - prod.min() > tolerance:
#                     kept = False

#         if kept and flag:
#             data[init_state]['kept_magmoms'] = True
#         else:
#             data[init_state]['kept_magmoms'] = False
#             red.append(init_state)

#         data[init_state]['energy'] = float(values[-1].split('=')[1]) / len(initial)     # energy per atom

# if len(red) != 0:
#     print(f"The following {len(red)} configuration(s) did not keep the original magmoms and will be marked in red "
#           f"on the graph: {', '.join(red)}")
# if len(not_converged) != 0:
#     print(f"The energy calculation of the following {len(not_converged)} configuration(s) did not converge and will "
#           f"not be shown on the graph: {', '.join(not_converged)}")
# if len(not_found) != 0:
#     print(f"The following {len(not_found)} configuration(s) reported an error during energy calculation and will "
#           f"not be shown on the graph: {', '.join(not_found)}")

# # plt.hist(np.abs(all_final_magmoms), bins=40)
# # plt.savefig('spin_distribution.png')

# setting = 1
# final_states = []
# final_setting = setting
# final_energies = []
# final_min = np.inf
# while setting < max_setting:
#     current_states = []
#     current_energies = []
#     for init_state, value in data.items():
#         if init_state != 'nm' and value['setting'] == setting and value['kept_magmoms']:
#             current_states.append(np.sign(value['init_spins']).tolist())
#             current_energies.append(value['energy'])
#     if min(current_energies) < final_min:
#         final_setting = setting
#         final_states = current_states
#         final_energies = current_energies
#         final_min = min(current_energies)
#     setting += 1

# # filter out configurations containing NM states
# tc_states = []
# tc_energies = []
# for state, energy in zip(final_states, final_energies):
#     if 0 not in state:
#         tc_states.append(state)
#         tc_energies.append(energy)

# # write states to file
# with open(f"{structure.formula.replace(' ','')}_{calculator}_states{final_setting:03d}.txt", 'wt') as f:
#     json.dump(tc_states, f)

# # write energies to file
# with open(f"{structure.formula.replace(' ','')}_{calculator}_energies{final_setting:03d}.txt", 'wt') as f:
#     json.dump(tc_energies, f)

# # copy setting file with geometry
# shutil.copy(f"{trials_path}/setting{final_setting:03d}.vasp", f"./{structure.formula.replace(' ','')}_{calculator}_setting{final_setting:03d}.vasp")

# # extract values for plot
# bar_labels = []
# energies = []
# kept_magmoms = []
# for key, value in sorted(data.items(), key=better_sort):
#     bar_labels.append(key)
#     energies.append(value['energy'])
#     kept_magmoms.append(value['kept_magmoms'])

# energies = np.array(energies)
# kept_magmoms = np.array(kept_magmoms)

# # energies from eV/atom to meV/atom
# energies *= 1000

# # normalize energies for plot
# energies -= min(energies)
# toplim = max(energies) * 1.15
# bottomlim = -0.1 * max(energies)

# energies += 0.1 * max(energies)
# x_axis = np.arange(1, len(energies) + 1)

# # split into chunks if too many configurations
# n_chunks = len(energies) // 14 + 1
# x_axis_split = np.array_split(x_axis, n_chunks)
# energies_split = np.array_split(energies, n_chunks)
# kept_magmoms_split = np.array_split(kept_magmoms, n_chunks)
# bar_labels_split = np.array_split(bar_labels, n_chunks)

# for i, (X, Y, kept_magmoms_chunk, bar_labels_chunk) in \
#         enumerate(zip(x_axis_split, energies_split, kept_magmoms_split, bar_labels_split)):
#     # increase matplotlib pyplot font size
#     plt.rcParams.update({'font.size': 20})

#     # set figure size
#     plt.figure(figsize=(16, 9))

#     # plot results
#     plt.bar(X[~kept_magmoms_chunk], Y[~kept_magmoms_chunk], bottom=bottomlim, color='r')
#     plt.bar(X[kept_magmoms_chunk], Y[kept_magmoms_chunk], bottom=bottomlim, color='b')

#     # label bars
#     ax = plt.gca()
#     rects = ax.patches
#     rects = sorted(rects, key=lambda x: x.get_x())
#     for bar_label, rect in zip(bar_labels_chunk, rects):
#         height = rect.get_height()
#         ax.text(rect.get_x() + rect.get_width() / 2, height + 0.85 * bottomlim, bar_label,
#                 fontsize='medium', ha='center', va='bottom', rotation='vertical')

#     # label axes
#     plt.xlabel('configurations')
#     plt.ylabel(f'relative energy (meV/atom)')
#     plt.xticks(X)
#     plt.ylim(top=toplim)

#     # save or show bar chart
#     plt.savefig(f"{structure.formula.replace(' ','')}_{calculator}_stability{i + 1:02d}.png", bbox_inches='tight')
#     plt.close()

# print(f'The most stable configuration is {bar_labels[np.argmin(energies)]}.')
# if np.logical_not(kept_magmoms[np.argmin(energies)]):
#     print('WARNING: values of initial and final magnetic moments of the most stable configuration '
#           'significantly differ.')

# # moving files and figures into a separate folder
# os.system(f"mkdir {structure.formula.replace(' ','')}_{calculator}")
# os.system(f"cp -f input.py {structure.formula.replace(' ','')}_{calculator}/input_{calculator}.py")
# os.system(f"mv -f {structure.formula.replace(' ','')}_{calculator}_stability*.png {structure.formula.replace(' ','')}_{calculator}")
# try:
#     os.system(f"mv -f trials_{structure.formula.replace(' ','')} {structure.formula.replace(' ','')}_{calculator}")
# except:
#     pass
# os.system(f"mv -f {structure.formula.replace(' ','')}_{calculator}_setting*.vasp {structure.formula.replace(' ','')}_{calculator}")
# os.system(f"mv -f {structure.formula.replace(' ','')}_{calculator}_energies*.txt {structure.formula.replace(' ','')}_{calculator}")
# os.system(f"mv -f {structure.formula.replace(' ','')}_{calculator}_states*.txt {structure.formula.replace(' ','')}_{calculator}")

# # acrhiving the folder
# os.system(f"tar -zcvf {structure.formula.replace(' ','')}_{calculator}.tar.gz {structure.formula.replace(' ','')}_{calculator}")
