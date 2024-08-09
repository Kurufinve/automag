"""
automag.common.SubmitManual
=============================

Class which takes care of writing input files for ab initio calculations.

.. codeauthor:: Daniil Poletaev <poletaev.dan@gmail.com>
"""

import numpy as np

from copy import copy
from ase.io import read
from typing import Union
from itertools import product
import os
from ase.calculators.vasp import Vasp


# from common.utilities import atoms_to_encode #, VaspCalculationTask, WriteOutputTask, WriteChargesTask

class SubmitManual(object):
    def __init__(self, poscar_file: str, mode: str, fix_params: dict, magmoms: list,
                 encut_values: Union[list, range] = None, sigma_values: Union[list, range] = None,
                 kpts_values: Union[list, range] = None, pert_values: Union[list, range] = None,
                 name: str = None, dummy_atom: str = None, dummy_position: int = None, 
                 calculator: str = 'vasp', jobheader: str = '#!/bin/bash', calculator_command: str = None, 
                 environment_activate: str = None, environment_deactivate: str = None):
        if mode == 'encut':
            assert encut_values is not None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = product(encut_values)
        elif mode == 'kgrid':
            assert encut_values is None
            assert sigma_values is not None
            assert kpts_values is not None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = product(sigma_values, kpts_values)
        elif mode == 'perturbations':
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is not None
            assert dummy_atom is not None
            assert dummy_position is not None
            self.var_params = []
        elif mode == 'singlepoint':
            assert encut_values is None
            assert sigma_values is None
            assert kpts_values is None
            assert pert_values is None
            assert dummy_atom is None
            assert dummy_position is None
            self.var_params = []
        else:
            raise ValueError(f'Value of mode = {mode} not understood.')

        if mode == 'encut':
            self.energy_convergence = True
        else:
            self.energy_convergence = False

        self.magmoms = np.array(magmoms)
        self.mode = mode
        self.fix_params = fix_params
        self.poscar_file = poscar_file
        self.name = name
        self.pert_values = pert_values
        self.dummy_atom = dummy_atom
        self.dummy_position = dummy_position
        self.calculator = calculator
        self.jobheader = jobheader
        self.calculator_command = calculator_command
        self.environment_activate = environment_activate
        self.environment_deactivate = environment_deactivate

    def submit(self):
        params = copy(self.fix_params)
        if self.var_params:
            for values in self.var_params:
                if len(values) == 1:
                    if self.mode == 'encut':
                        params['encut'] = values[0]
                    # else:
                    #     params['ldauu'] = [values[0], 0, 0]
                    #     params['ldauj'] = [values[0], 0, 0]
                    name = self.mode + str(values[0])
                elif len(values) == 2:
                    params['sigma'] = values[0]
                    params['kpts'] = values[1]
                    name = self.mode + str(values[0]) + '-' + str(values[1])
                else:
                    raise ValueError('Convergence tests for three parameters simultaneously are not supported.')

                self.write_input(params, name)
        else:
            if self.name is not None:
                self.write_input(params, f'{self.name}')
            else:
                self.write_input(params, self.mode)

    def write_input(self, params, name):


        # create an atoms object as a class atribute to make it accesible across the class methods and encode it
        self.atoms = read(self.poscar_file)
        if self.mode == 'perturbations':
            ch_symbols = self.atoms.get_chemical_symbols()
            atom_ucalc = ch_symbols[self.dummy_position]
            ch_symbols[self.dummy_position] = self.dummy_atom
            self.atoms.set_chemical_symbols(ch_symbols)
            # encode = atoms_to_encode(self.atoms)
        else:
            # encode = atoms_to_encode(self.atoms)
            pass

        # create working directory for manual submission
        calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
        compound_dir = os.path.join(calcfold, f"{self.atoms.get_chemical_formula(mode='metal')}")
        state_dir = os.path.join(compound_dir, f'{self.calculator}/{name}')
        try:
            os.mkdir(state_dir)
        except:
            pass

        if self.mode == 'singlepoint'
            if name == 'nm':

                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(state_dir,'singlepoint',params)



            else:

                if self.calculator == "vasp":
                    # single-point run + recalc run
                    self.write_vasp_input_files(state_dir,'singlepoint+recalc',params)
                    # # recalc run with magmoms from previous run
                    # self.write_vasp_input_files(state_dir,'recalc',params)

            # appending a line to the submit script 
            script_name = f'{compound_dir}/{self.calculator}/run_all.sh'
            if not os.path.exists(script_name):
                with open(script_name,'w') as f:
                    f.write('#!/bin/sh\n')
                os.system(f'chmod +x {script_name}')
            submit_command = f'cd {state_dir}; sbatch jobscript\n'
            write_submit_script(script_name,submit_command)


        elif self.mode == 'perturbations':
            next_id = 2
            nsc_fireworks = []
            sc_fireworks = []
            out_fireworks = []
            for perturbation in self.pert_values:
                nsc_firetask = VaspCalculationTask(
                    calc_params=params,
                    encode=encode,
                    magmoms=self.magmoms,
                    pert_step='NSC',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                    atom_ucalc=atom_ucalc,
                )

                nsc_firework = Firework(
                    [nsc_firetask],
                    name='nsc',
                    spec={'_pass_job_info': True},
                    fw_id=next_id,
                )

                nsc_fireworks.append(nsc_firework)
                next_id += 1

                sc_firetask = VaspCalculationTask(
                    calc_params=params,
                    encode=encode,
                    magmoms=self.magmoms,
                    pert_step='SC',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                    atom_ucalc=atom_ucalc,
                )

                sc_firework = Firework(
                    [sc_firetask],
                    name='sc',
                    spec={'_pass_job_info': True},
                    fw_id=next_id,
                )

                sc_fireworks.append(sc_firework)
                next_id += 1

                out_firetask = WriteChargesTask(
                    filename='charges.txt',
                    pert_value=perturbation,
                    dummy_atom=self.dummy_atom,
                )
                out_firework = Firework(
                    [out_firetask],
                    name='write_charges',
                    spec={'_queueadapter': {'ntasks': 1, 'walltime': '00:30:00'}},
                    fw_id=next_id,
                )

                out_fireworks.append(out_firework)
                next_id += 1

            fireworks.append(nsc_fireworks)
            fireworks.append(sc_fireworks)
            fireworks.append(out_fireworks)
        else:

    def set_magmoms(self,params):

        if self.magmoms.any() == 'previous':
            # assert 'encode' not in self
            magmoms = self.atoms.get_magnetic_moments().round()
        else:
            magmoms = np.array(self.magmoms)

        if self.magmoms.any():
            # add the necessary VASP parameters
            params['lorbit'] = 11
            params['ispin'] = 2
            self.atoms.set_initial_magnetic_moments(magmoms)

        return params


    def set_pert_step(self):
        atom_ucalc = Element(self['atom_ucalc'])
        elem_list_sorted, indices = np.unique(self.atoms.get_chemical_symbols(), return_index=True)
        elem_list = elem_list_sorted[np.argsort(indices)]

        self['calc_params']['ldaul'] = []
        self['calc_params']['ldauu'] = []
        self['calc_params']['ldauj'] = []
        for atom in elem_list:
            if atom == self['dummy_atom']:
                if atom_ucalc.is_transition_metal:
                    self['calc_params']['ldaul'].append(2)
                elif atom_ucalc.is_lanthanoid or atom_ucalc.is_actinoid:
                    self['calc_params']['ldaul'].append(3)
                else:
                    raise ValueError(f"Cannot calculate U for the element {self['atom_ucalc']}")

                self['calc_params']['ldauu'].append(self['pert_value'])
                self['calc_params']['ldauj'].append(self['pert_value'])
            else:
                self['calc_params']['ldaul'].append(-1)
                self['calc_params']['ldauu'].append(0)
                self['calc_params']['ldauj'].append(0)        


    def write_vasp_input_files(self,state_dir,innermode='singlepoint',params = None):

        # if magmoms in input
        if 'magmoms' in self.__dict__.keys():
            # if self.magmoms != None:
            # self.set_magmoms()
            params = self.set_magmoms(params)

        # if perturbation calculation
        if 'pert_step' in self.__dict__.keys():
            # if self.pert_step != None:
            self.set_pert_step()


        write_magmoms_script_name = 'get_magmoms_vasp.py'
        write_magmoms_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{write_magmoms_script_name}')

        write_vasp_output_script_name = 'write_output_vasp.py'
        write_vasp_output_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{write_vasp_output_script_name}')

        if innermode == 'singlepoint' or innermode == 'recalc':
            workdir =  os.path.join(state_dir, mode)

            try:
                os.mkdir(workdir)
            except:
                pass

            os.system(f'cp {write_magmoms_script_path} {workdir}')
            os.system(f'cp {write_vasp_output_script_path} {workdir}')

            # creating ase.calculator object and writing input   
            calc = Vasp(atoms=self.atoms,
                        directory=workdir,
                        **params)

            calc.write_input(self.atoms)

            # create input script
            if innermode == 'singlepoint':
                self.write_initial_magmoms(f'{workdir}/initial_magmoms.txt')
                with open(f'{state_dir}/jobscript','w') as f:
                    f.write(self.jobheader)
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                    f.write('\n')
                    f.write('\n')
                    f.write(f'cd {workdir}\n')
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write(f'cp {write_vasp_output_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    f.write(f'python {write_vasp_output_script_name}\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)

            elif innermode == 'recalc':
                # write_magmoms_script_name = 'get_magmoms_vasp.py'
                # write_magmoms_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{write_magmoms_script_name}')
                with open(f'{state_dir}/jobscript','w') as f:
                    f.write(self.jobheader)
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                    f.write('\n')
                    f.write('\n')
                    f.write(f'cd {workdir}\n')
                    f.write('\n')
                    f.write(f'cp {write_magmoms_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    f.write(f'python {write_magmoms_script_name}\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write('\n')
                    f.write(f'cp {write_vasp_output_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    f.write(f'python {write_vasp_output_script_name}\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)


        elif innermode == 'singlepoint+recalc':

            workdir_singlepoint =  os.path.join(state_dir, 'singlepoint')
            workdir_recalc =  os.path.join(state_dir, 'recalc')

            try:
                os.mkdir(workdir_singlepoint)
                os.mkdir(workdir_recalc)
            except:
                pass


            # creating ase.calculator object for singlepoint run and writing input   
            calc_singlepoint = Vasp(atoms=self.atoms,
                        directory=workdir_singlepoint,
                        **params)

            calc_singlepoint.write_input(self.atoms)

            # creating ase.calculator object for recalc run and writing input   
            calc_recalc = Vasp(atoms=self.atoms,
                        directory=workdir_recalc,
                        **params)

            calc_recalc.write_input(self.atoms)



            os.system(f'cp {write_magmoms_script_path} {workdir_recalc}')
            os.system(f'cp {write_vasp_output_script_path} {workdir_recalc}')


            self.write_initial_magmoms(f'{workdir_singlepoint}/initial_magmoms.txt')
            with open(f'{state_dir}/jobscript','w') as f:
                f.write(self.jobheader)
                f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                f.write('\n')
                f.write('\n')
                f.write(f'cd {workdir_singlepoint}\n')
                f.write('\n')
                f.write(self.calculator_command)
                f.write('\n')
                f.write(f'cd {workdir_recalc}\n')
                f.write('\n')
                f.write(f'cp {write_magmoms_script_path} .\n')
                f.write(self.environment_activate)
                f.write('\n')
                f.write(f'python {write_magmoms_script_name}\n')
                f.write('\n')
                f.write(self.environment_deactivate)
                f.write('\n')
                f.write(self.calculator_command)
                f.write('\n')
                f.write(f'cp {write_vasp_output_script_path} .\n')
                f.write(self.environment_activate)
                f.write('\n')
                f.write(f'python {write_vasp_output_script_name}\n')
                f.write('\n')
                f.write(self.environment_deactivate)
        return

    def write_initial_magmoms(self,path_to_file):

        with open(path_to_file,'w') as f:
            for magmom in self.magmoms:
                f.write(f'{magmom} ')

        return

""" Helping functions """

def write_submit_script(script_name,submit_command):

    with open(script_name,'a') as f:
        f.write(submit_command)