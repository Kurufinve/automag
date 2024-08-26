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
from pymatgen.core.periodic_table import Element


# from common.utilities import atoms_to_encode #, VaspCalculationTask, WriteOutputTask, WriteChargesTask

class SubmitManual(object):
    def __init__(self, poscar_file: str, mode: str, fix_params: dict, magmoms: list,
                 encut_values: Union[list, range] = None, sigma_values: Union[list, range] = None,
                 kpts_values: Union[list, range] = None, pert_values: Union[list, range] = None,
                 name: str = None, dummy_atom: str = None, dummy_position: int = None, 
                 calculator: str = 'vasp', jobheader: str = '#!/bin/bash', calculator_command: str = None, 
                 environment_activate: str = None, environment_deactivate: str = None, exec_command: str = 'exec'):
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
        self.exec_command = exec_command
        self.calculator_command = calculator_command
        self.environment_activate = environment_activate
        self.environment_deactivate = environment_deactivate

        self.write_charges_script_name = 'write_charges.py'
        self.write_charges_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_charges_script_name}')

        self.write_magmoms_script_name = 'write_magmoms_input.py'
        self.write_magmoms_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_magmoms_script_name}')

        self.write_output_script_name = 'write_output.py'
        self.write_output_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_output_script_name}')        

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

        # create working directory for manual submission
        calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
        compound_dir = os.path.join(calcfold, f"{self.atoms.get_chemical_formula(mode='metal')}")
        state_dir = os.path.join(compound_dir, f'{self.calculator}/{name}')
        try:
            os.mkdir(state_dir)
        except:
            pass

        if self.mode == 'singlepoint':
            if name == 'nm':

                # write 
                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(state_dir,'singlepoint',params)

                # write jobscript
                with open(f'{state_dir}/jobscript','w') as f:
                    f.write(self.jobheader)
                    if '#SBATCH' in self.jobheader:
                        print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}_{self.name} calculation")
                        f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                    f.write('\n')
                    f.write('\n')
                    f.write(f'cd {state_dir}/singlepoint \n')
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write(f'rm WAVECAR CHG\n') # remove large files
                    f.write(f'cp {self.write_output_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    if self.calculator == "vasp":
                        f.write(f'python {self.write_output_script_name} vasp\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)
            else:

                if self.calculator == "vasp":
                    # single-point run + recalc run
                    self.write_vasp_input_files(state_dir,'singlepoint+recalc',params)

                # write jobscript
                with open(f'{state_dir}/jobscript','w') as f:
                    f.write(self.jobheader)
                    if '#SBATCH' in self.jobheader:
                        print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}_{self.name} calculation")
                        f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                    f.write('\n')
                    f.write('\n')
                    f.write(f'cd {state_dir}/singlepoint\n')
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write(f'rm WAVECAR CHG\n') # remove large files
                    f.write(f'cd {state_dir}/recalc\n')
                    f.write('\n')
                    f.write(f'cp {self.write_magmoms_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    if self.calculator == "vasp":
                        f.write(f'python {self.write_magmoms_script_name} vasp\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write(f'rm WAVECAR CHG\n') # remove large files
                    f.write(f'cp {self.write_output_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    if self.calculator == "vasp":
                        f.write(f'python {self.write_output_script_name} vasp\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)

            # appending a line to the submit script 
            script_name = f'{compound_dir}/{self.calculator}/run_all.sh'
            if not os.path.exists(script_name):
                with open(script_name,'w') as f:
                    f.write('#!/bin/sh\n')
                os.system(f'chmod +x {script_name}')
            submit_command = f'cd {state_dir}; {self.exec_command} jobscript\n'
            write_submit_script(script_name,submit_command)


        elif self.mode == 'perturbations':

            self.ch_symbols = self.atoms.get_chemical_symbols()
            self.atom_ucalc = self.ch_symbols[self.dummy_position]
            self.ch_symbols[self.dummy_position] = self.dummy_atom
            # state_dir_pert = os.path.join(state_dir, 'perturbations')
            state_dir_pert = state_dir
            self.atoms.set_chemical_symbols(self.ch_symbols)

            print(f' atom_ucalc = {self.atom_ucalc}')
            print(f' dummy_atom = {self.dummy_atom}')

            # resolving the path to VASP pseudos
            if self.calculator == "vasp":
                if 'xc' in params.keys():
                    if params['xc'] == 'PBE':
                        pp_path = os.path.join(os.environ.get('VASP_PP_PATH'), 'potpaw_PBE')
                    elif params['xc'] == 'LDA':
                        pp_path = os.path.join(os.environ.get('VASP_PP_PATH'), 'potpaw_LDA')
                else:
                    pp_path = os.path.join(os.environ.get('VASP_PP_PATH'), 'potpaw_PBE')

                # checking if the '_sv', '_pv', '_GW' or other version of PAW potential is used
                if 'setups' in params.keys():
                    if isinstance(params['setups'], dict):
                        if self.atom_ucalc in params['setups'].keys():
                            potpaw_name = self.atom_ucalc + params['setups'][self.atom_ucalc]
                        else:
                            potpaw_name = self.atom_ucalc
                    else:
                        potpaw_name = self.atom_ucalc
                else:
                    potpaw_name = self.atom_ucalc

                atom_ucalc_pp_path = os.path.join(pp_path, f'{potpaw_name}/POTCAR')
                dummy_atom_pp_path = os.path.join(pp_path, f'{self.dummy_atom}/POTCAR')

                dummy_atom_pp_path_temp = os.path.join(pp_path, f'{self.dummy_atom}/_POTCAR')

                print(f'atom_ucalc_pp_path = {atom_ucalc_pp_path}')
                print(f'dummy_atom_pp_path = {dummy_atom_pp_path}')

                # temporary renaming the dummy atom POTCAR
                os.system(f'mv {dummy_atom_pp_path} {dummy_atom_pp_path_temp}')
                
                # copying the ucalc atom POTCAR to the dummy atom folder 
                os.system(f'cp {atom_ucalc_pp_path} {dummy_atom_pp_path}')




            # write files for common singlepoint calculation for WAVECAR and CHGCAR
            self.write_vasp_input_files(state_dir_pert,'singlepoint',params)

            # start writing common jobscript for singlepoint calculation and calculations for all perturbations
            with open(f'{state_dir}/jobscript','w') as f:
                f.write(self.jobheader)
                if '#SBATCH' in self.jobheader:
                    print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}_{self.name} calculations")
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                f.write('\n')
                f.write('\n')
                f.write(f'cd {state_dir_pert}/singlepoint\n')
                f.write('\n')
                f.write(self.calculator_command)
                f.write('\n')

            # write files for perturbation calculations at different pert values
            for perturbation in self.pert_values:
                print(f'perturbation: {perturbation}')
                pert_dir = os.path.join(state_dir_pert, f'{perturbation}')

                self.write_vasp_input_files(pert_dir,'nsc',params,pert_value=perturbation)
                self.write_vasp_input_files(pert_dir,'sc',params,pert_value=perturbation)

                # appending the part to the jobscript
                with open(f'{state_dir}/jobscript','a') as f:
                    f.write('\n')
                    f.write(f'cd {pert_dir}/nsc\n')
                    f.write(f'cp ../../singlepoint/WAVECAR . \n')
                    f.write(f'cp ../../singlepoint/CHGCAR . \n')
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')
                    f.write(f'cd {pert_dir}/sc\n')
                    f.write(f'cp ../../singlepoint/WAVECAR . \n')
                    f.write('\n')
                    f.write(self.calculator_command)
                    f.write('\n')

            # os.system(f'cp {self.write_charges_script_path} {state_dir_pert}')
            with open(f'{state_dir}/jobscript','a') as f:
                f.write(f'cd {state_dir}\n')
                f.write(f'cp {self.write_charges_script_path} .\n')
                f.write(self.environment_activate)
                f.write('\n')
                if self.calculator == "vasp":
                    f.write(f'python {self.write_charges_script_name} vasp\n')
                f.write('\n')
                f.write(self.environment_deactivate)

            # restoring the name of the dummy atom POTCAR
            os.system(f'mv {dummy_atom_pp_path_temp} {dummy_atom_pp_path}')

        else:
            pass

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


    def write_vasp_input_files(self,state_dir,innermode='singlepoint',params = None, pert_value = None):

        # if magmoms in input
        if 'magmoms' in self.__dict__.keys():
            # if self.magmoms != None:
            # self.set_magmoms()
            params = self.set_magmoms(params)

        # # if perturbation calculation
        # if 'pert_step' in self.__dict__.keys():
        #     # if self.pert_step != None:
        #     self.set_pert_step()



        if innermode == 'singlepoint' or innermode == 'recalc':
            workdir =  os.path.join(state_dir, innermode)

            try:
                os.mkdir(workdir)
            except:
                pass

            # os.system(f'cp {write_magmoms_script_path} {workdir}')
            # os.system(f'cp {write_vasp_output_script_path} {workdir}')

            # creating ase.calculator object and writing input   
            calc = Vasp(atoms=self.atoms,
                        directory=workdir,
                        **params)

            calc.write_input(self.atoms)

            # create input script
            if innermode == 'singlepoint':

                if self.mode == 'singlepoint':
                    self.write_initial_magmoms(f'{workdir}/initial_magmoms.txt')


            elif innermode == 'recalc':
                # # write_magmoms_script_name = 'get_magmoms_vasp.py'
                # # write_magmoms_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{write_magmoms_script_name}')
                # with open(f'{state_dir}/jobscript','w') as f:
                #     f.write(self.jobheader)
                #     f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}_{self.name}")
                #     f.write('\n')
                #     f.write('\n')
                #     f.write(f'cd {workdir}\n')
                #     f.write('\n')
                #     f.write(f'cp {self.write_magmoms_script_path} .\n')
                #     f.write(self.environment_activate)
                #     f.write('\n')
                #     f.write(f'python {self.write_magmoms_script_name}\n')
                #     f.write('\n')
                #     f.write(self.environment_deactivate)
                #     f.write('\n')
                #     f.write(self.calculator_command)
                #     f.write('\n')
                #     f.write('\n')
                #     # write output information only in the case of 2_coll run
                #     if self.mode == 'singlepoint':
                #         f.write(f'cp {self.write_vasp_output_script_path} .\n')
                #         f.write(self.environment_activate)
                #         f.write('\n')
                #         f.write(f'python {self.write_vasp_output_script_name}\n')
                #         f.write('\n')
                #         f.write(self.environment_deactivate)
                pass


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

            self.write_initial_magmoms(f'{workdir_singlepoint}/initial_magmoms.txt')


        elif innermode == 'nsc' or innermode == 'sc':

            workdir =  os.path.join(state_dir, innermode)

            try:
                os.mkdir(workdir)
            except:
                pass

            atom_ucalc = Element(self.atom_ucalc)
            elem_list_sorted, indices = np.unique(self.atoms.get_chemical_symbols(), return_index=True)
            elem_list = elem_list_sorted[np.argsort(indices)]

            params['ldaul'] = []
            params['ldauu'] = []
            params['ldauj'] = []
            for atom in elem_list:
                if atom == self.dummy_atom:
                    if atom_ucalc.is_transition_metal:
                        params['ldaul'].append(2)
                    elif atom_ucalc.is_lanthanoid or atom_ucalc.is_actinoid:
                        params['ldaul'].append(3)
                    else:
                        raise ValueError(f"Cannot calculate U for the element {self.atom_ucalc}")

                    params['ldauu'].append(pert_value)
                    params['ldauj'].append(pert_value)
                else:
                    params['ldaul'].append(-1)
                    params['ldauu'].append(0)
                    params['ldauj'].append(0)

            params['ldau'] = True
            params['ldautype'] = 3

            if innermode == 'nsc':
                params['icharg'] = 11
            elif  innermode == 'sc':
                params['icharg'] = 0

            # creating ase.calculator object and writing input   
            calc = Vasp(atoms=self.atoms,
                        directory=workdir,
                        **params)

            calc.write_input(self.atoms)

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