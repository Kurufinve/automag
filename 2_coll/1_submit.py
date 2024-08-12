"""
automag.2_coll.1_submit
=======================

Script which runs enumlib and submits calculations.

.. first codeauthor:: Michele Galasso <m.galasso@yandex.com>
.. second codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

# default values for some input variables
use_fireworks = True
calculator = 'vasp'

cwd = os.getcwd()

try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    from input_file import *
except IndexError:
    print(f'Using the default input.py file from folder: {cwd}')
    from input import *

import os
import subprocess
import numpy as np
import shutil

from itertools import product
from pymatgen.io.vasp import Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from datetime import datetime

if use_fireworks:
    from common.SubmitFirework import SubmitFirework
else:
    from common.SubmitManual import SubmitManual


def launch_enumlib(count, split):
    os.mkdir(f'enumlib{count}')
    os.chdir(f'enumlib{count}')

    with open('struct_enum.in', 'w') as f:
        f.write('generated by Automag\n')
        f.write('bulk\n')

        for lat_vector in symmetrized_structure.lattice.matrix:
            for component in lat_vector:
                f.write(f'{component:14.10f}        ')
            f.write('\n')

        case = len(split) + sum(split)
        f.write(f'  {case} -nary case\n')
        f.write(f'    {symmetrized_structure.num_sites} # Number of points in the multilattice\n')

        offset = 0
        for i, (s, wyckoff) in enumerate(zip(split, symmetrized_structure.equivalent_sites)):
            if s:
                offset += 1

            for atom in wyckoff:
                for component in atom.coords:
                    f.write(f'{component:14.10f}        ')
                if s:
                    f.write(f'{i + offset - 1}/{i + offset}\n')
                else:
                    f.write(f'{i + offset}\n')

        f.write(f'    1 {supercell_size}   # Starting and ending cell sizes for search\n')
        f.write('0.10000000E-06 # Epsilon (finite precision parameter)\n')
        f.write('full list of labelings\n')
        f.write('# Concentration restrictions\n')

        for s, wyckoff in zip(split, symmetrized_structure.equivalent_sites):
            if s:
                for _ in range(2):
                    f.write(f'{len(wyckoff):4d}')
                    f.write(f'{len(wyckoff):4d}')
                    f.write(f'{symmetrized_structure.num_sites * 2:4d}\n')
            else:
                f.write(f'{len(wyckoff) * 2:4d}')
                f.write(f'{len(wyckoff) * 2:4d}')
                f.write(f'{symmetrized_structure.num_sites * 2:4d}\n')

    # process = subprocess.Popen('/home/michele/softs/enumlib/src/enum.x')
    process = subprocess.Popen('enum.x')
    try:
        process.wait(timeout=60)
    except subprocess.TimeoutExpired:
        process.kill()

    # solve a bug of enumlib which produces an extra new line
    with open('struct_enum.out', 'rt') as file:
        lines = file.readlines()

    line_number = 0
    os.remove('struct_enum.out')
    with open('struct_enum.out', 'wt') as file:
        while line_number < len(lines):
            line = lines[line_number]
            if '(Non)Equivalency list' not in line:
                file.write(line)
                line_number += 1
            else:
                file.write(line.strip('\n'))

                offset = 0
                while not lines[line_number + offset + 1].startswith('start'):
                    file.write(' ' + lines[line_number + offset + 1].strip('\n'))
                    offset += 1

                file.write('\n')
                line_number += offset + 1

    # os.system('/home/michele/softs/enumlib/aux_src/makeStr.py 1 500')
    os.system('makeStr.py 1 500')

    for j in range(501):
        if os.path.isfile(f'vasp.{j + 1}'):
            conf_poscar = Poscar.from_file(f'vasp.{j + 1}')
        else:
            break

        if conf_poscar.structure.lattice in lattices:
            index = lattices.index(conf_poscar.structure.lattice)
            current_coords = conf_poscar.structure.frac_coords.tolist()
            reference_coords = coordinates[index].tolist()
            mapping = [np.nonzero([np.allclose(cur, ref) for cur in current_coords])[0][0] for ref in reference_coords]
        else:
            lattices.append(conf_poscar.structure.lattice)
            coordinates.append(conf_poscar.structure.frac_coords)
            configurations.append([])
            index = len(lattices) - 1
            mapping = list(range(len(conf_poscar.structure.frac_coords)))

        counter = 0
        split_groups = []
        site_magmoms = []
        for to_split, states in zip(split, wyckoff_magmoms):
            if to_split:
                split_groups.append([counter, counter + 1])
                for _ in range(2):
                    site_magmoms.append([1, -1])
                    counter += 1
            else:
                site_magmoms.append(states)
                counter += 1

        # adapt equivalent_multipliers in case of supercells
        coefficient = len(conf_poscar.structure) // len(structure)

        current_multipliers = []
        for multiplier in equivalent_multipliers:
            current_multipliers.append(np.repeat(multiplier, coefficient))

        # use a flag to check that all sites which have been split are AFM
        for conf in product(*site_magmoms):
            flag = True
            conf_array = np.array(conf)
            for group in split_groups:
                if sum(conf_array[group]) != 0:
                    flag = False
            if flag:
                configuration = np.repeat(conf, conf_poscar.natoms)
                transformed_configuration = configuration[mapping]

                for mult in current_multipliers:
                    candidate_conf = np.multiply(transformed_configuration, mult).tolist()

                    # if all spins are zero do not include the NM configuration
                    if len(np.nonzero(candidate_conf)[0]) > 0:
                        first_nonzero_index = np.nonzero(candidate_conf)[0][0]

                        # if the first non-zero spin is negative flip all spins
                        if candidate_conf[first_nonzero_index] < 0:
                            candidate_conf = [-item for item in candidate_conf]

                        # add to list of configurations for the current settings
                        if candidate_conf not in configurations[index]:
                            configurations[index].append(candidate_conf)

    os.chdir('..')


# full path to poscar file
# path_to_poscar = '../geometries/' + poscar_file
path_to_poscar = os.path.join(os.environ.get('AUTOMAG_PATH'),f'/geometries/{poscar_file}')

# create Structure and SymmetrizedStructure objects
structure = Structure.from_file(path_to_poscar)
analyzer = SpacegroupAnalyzer(structure)
symmetrized_structure = analyzer.get_symmetrized_structure()

# find out which atoms are magnetic
for element in structure.composition.elements:
    if element.name in spin_values:
        element.is_magnetic = True
    else:
        element.is_magnetic = False

# print(dir(structure))
# print(structure.formula.replace(' ',''))

if os.path.exists(f"trials_{structure.formula.replace(' ','')}"):
    # print('Cannot create a folder named trials: an object with the same name already exists.')
    # exit()
    # time = datetime.now().strftime('%Y-%m-%d%H_%M_%S')
    # print(f'A folder named trials already exists. Moving it to the folder: trials_{structure.formula}_{time}')
    # os.system(f"mv trials trials_{structure.formula}_{time}")
    print(f'A folder named trials already exists. moving it to another folder!')
    os.system(f"rm -rf trials_{structure.formula.replace(' ','')}_old")
    try:
        shutil.move(f"trials_{structure.formula.replace(' ','')}", f"trials_{structure.formula.replace(' ','')}_old")
    except:
        os.system(f"rm -rf trials_{structure.formula.replace(' ','')}_old")
        shutil.move(f"trials_{structure.formula.replace(' ','')}", f"trials_{structure.formula.replace(' ','')}_old")
    # os.system(f"mv trials_{structure.formula.replace(' ','')} trials_{structure.formula.replace(' ','')}_old")
    os.system(f"rm -rf trials_{structure.formula.replace(' ','')}")



os.mkdir(f"trials_{structure.formula.replace(' ','')}")
os.chdir(f"trials_{structure.formula.replace(' ','')}")

# geometrical settings and respective lists of magnetic configurations
lattices = [structure.lattice]
coordinates = [structure.frac_coords]
configurations = [[]]

# get the multiplicities of each Wyckoff position
multiplicities = [len(item) for item in symmetrized_structure.equivalent_indices]

wyckoff_magmoms = []
equivalent_multipliers = []
for multiplicity, wyckoff in zip(multiplicities, symmetrized_structure.equivalent_sites):
    if wyckoff[0].specie.is_magnetic:
        wyckoff_magmoms.append([1, 0, -1])

        if len(spin_values[wyckoff[0].specie.name]) == 2:
            val1 = spin_values[wyckoff[0].specie.name][0]
            val2 = spin_values[wyckoff[0].specie.name][1]

            if len(equivalent_multipliers) == 0:
                equivalent_multipliers.append(np.repeat(val1, multiplicity))
                equivalent_multipliers.append(np.repeat(val2, multiplicity))
            else:
                new_equivalent_multipliers = []
                for multiplier in equivalent_multipliers:
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val1, multiplicity)))
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val2, multiplicity)))

                equivalent_multipliers = new_equivalent_multipliers

        elif len(spin_values[wyckoff[0].specie.name]) == 1:
            val1 = spin_values[wyckoff[0].specie.name][0]

            if len(equivalent_multipliers) == 0:
                equivalent_multipliers.append(np.repeat(val1, multiplicity))
            else:
                new_equivalent_multipliers = []
                for multiplier in equivalent_multipliers:
                    new_equivalent_multipliers.append(np.append(multiplier, np.repeat(val1, multiplicity)))

                equivalent_multipliers = new_equivalent_multipliers

        else:
            raise ValueError('Max 2 spin values for each magnetic species.')

    else:
        wyckoff_magmoms.append([0])

        if len(equivalent_multipliers) == 0:
            equivalent_multipliers.append(np.repeat(0, multiplicity))
        else:
            new_equivalent_multipliers = []
            for multiplier in equivalent_multipliers:
                new_equivalent_multipliers.append(np.append(multiplier, np.repeat(0, multiplicity)))

            equivalent_multipliers = new_equivalent_multipliers

# get all possible configurations without any splitting
for conf in product(*wyckoff_magmoms):
    configuration = np.repeat(conf, multiplicities)

    for mult in equivalent_multipliers:
        candidate_conf = np.multiply(configuration, mult).tolist()

        # if the configuration is not NM, get the index of the first non-zero element
        if len(np.nonzero(candidate_conf)[0]) > 0:
            first_nonzero_index = np.nonzero(candidate_conf)[0][0]

            # if the first non-zero spin is negative flip all spins
            if candidate_conf[first_nonzero_index] < 0:
                candidate_conf = [-item for item in candidate_conf]

        # add to list of configurations for the current settings
        if candidate_conf not in configurations[0]:
            configurations[0].append(candidate_conf)

# split all possible combinations of Wyckoff positions
splits = []
possibilities = [[0, 1] if len(item) > 1 else [0] for item in wyckoff_magmoms]
for split in product(*possibilities):
    if sum(split) != 0:
        splits.append(split)

for i, split in enumerate(splits):
    launch_enumlib(i + 1, split)

# merge the first and second settings if they are equivalent
if len(lattices) > 1:
    transformation_matrix = np.dot(lattices[1].matrix, np.linalg.inv(lattices[0].matrix))
    inv_transformation_matrix = np.linalg.inv(transformation_matrix)
    determinant = int(round(np.linalg.det(transformation_matrix), 6))

    if determinant == 1:
        del_flag = True
        origin_shift = np.around(coordinates[1][0] - np.dot(coordinates[0][0], inv_transformation_matrix), decimals=6)
        for coord1, coord2 in zip(coordinates[0], coordinates[1]):
            if not np.allclose((np.dot(coord1, inv_transformation_matrix) + origin_shift) % 1, coord2):
                del_flag = False

        if del_flag:
            del lattices[0]
            del coordinates[0]
            configurations[0].extend(configurations[1])
            del configurations[1]

# write output and submit calculations
fm_count = 1
afm_count = 1
fim_count = 1
original_ch_symbols = [atom.name for atom in structure.species]
for i, (lattice, frac_coords, confs) in enumerate(zip(lattices, coordinates, configurations)):
    magnification = len(frac_coords) // len(structure.frac_coords)
    ch_symbols = np.repeat(original_ch_symbols, magnification)
    setting = Structure(lattice, ch_symbols, frac_coords)
    setting.to(fmt='poscar', filename=f'setting{i + 1:03d}.vasp')
    mask = [item.is_magnetic for item in setting.species]

    for conf in confs:
        conf_array = np.array(conf)
        with open(f'configurations{i + 1:03d}.txt', 'a') as f:
            if np.sum(np.abs(conf)) == 0:
                state = 'nm'
            elif min(conf) >= 0:
                state = 'fm' + str(fm_count)
                fm_count += 1
            elif np.sum(conf) == 0:
                state = 'afm' + str(afm_count)
                afm_count += 1
            else:
                state = 'fim' + str(fim_count)
                fim_count += 1

            f.write(f'{state:>6s}  ')
            f.write(' '.join(f'{e:2d}' for e in conf_array[mask]))
            f.write('\n')

        if use_fireworks:
            run = SubmitFirework(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
                                 name=state)
            run.submit()

        else:
            run = SubmitManual(f'setting{i + 1:03d}.vasp', mode='singlepoint', fix_params=params, magmoms=conf,
                                 name=state, calculator = calculator, jobheader = jobheader, calculator_command = calculator_command,
                                 environment_activate = environment_activate, environment_deactivate = environment_deactivate)
            run.submit()
