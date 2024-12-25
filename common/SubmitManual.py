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

from pymatgen.io.vasp.outputs import Oszicar


# from common.utilities import atoms_to_encode #, VaspCalculationTask, WriteOutputTask, WriteChargesTask

class SubmitManual(object):
    def __init__(self, poscar_file: str, mode: str, fix_params: dict, magmoms: list, struct_suffix = '', 
                 encut_values: Union[list, range] = None, sigma_values: Union[list, range] = None,
                 kpts_values: Union[list, range] = None, pert_values: Union[list, range] = None,
                 Nth = None, Nph = None, N_MAE = None, MAE_x = None, MAE_y = None, MAE_z = None, use_symmetries = False,
                 name: str = None, dummy_atom: str = None, dummy_position: int = None, 
                 calculator: str = 'vasp', jobheader: str = '#!/bin/bash', calculator_command: str = None, 
                 environment_activate: str = None, environment_deactivate: str = None, 
                 parallel_over_configurations = True):
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
        elif mode == 'mae_theta_phi':
            assert Nth is not None
            assert Nph is not None
            self.var_params = []
        elif mode == 'mae_curve':
            assert N_MAE is not None
            assert MAE_x is not None
            assert MAE_y is not None
            assert MAE_z is not None
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
        self.struct_suffix = struct_suffix
        self.pert_values = pert_values
        self.dummy_atom = dummy_atom
        self.dummy_position = dummy_position
        self.Nth = Nth
        self.Nph = Nph
        self.N_MAE = N_MAE # number of points on MAE curve (from lowest energy SAXIS to highest energy SAXIS)
        self.MAE_x = MAE_x 
        self.MAE_z = MAE_z 
        self.MAE_y = MAE_y 
        self.use_symmetries = use_symmetries # if using crystal symmetries for MAE calculations
        self.calculator = calculator
        self.jobheader = jobheader
        self.calculator_command = calculator_command
        self.environment_activate = environment_activate
        self.environment_deactivate = environment_deactivate
        self.parallel_over_configurations = parallel_over_configurations

        self.write_charges_script_name = 'write_charges.py'
        self.write_charges_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_charges_script_name}')

        self.write_magmoms_script_name = 'write_magmoms_input.py'
        self.write_magmoms_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_magmoms_script_name}')

        self.write_output_script_name = 'write_output.py'
        self.write_output_script_path = os.path.join(os.environ.get('AUTOMAG_PATH'), f'common/{self.write_output_script_name}')        

    def submit(self):
        params = copy(self.fix_params)
        if self.var_params:

            self.atoms = read(self.poscar_file)
            calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
            compound_dir = os.path.join(calcfold, f"{self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}")
            
            # creating compound_dir if it does not exist
            try:
                os.mkdir(compound_dir)
            except:
                pass

            calculator_dir = os.path.join(compound_dir, f'{self.calculator}')

            # creating calculator_dir if it does not exist
            try:
                os.mkdir(calculator_dir)
            except:
                pass


            convergence_dir = os.path.join(calculator_dir, 'convergence')

            # creating convergence_dir if it does not exist
            try:
                os.mkdir(convergence_dir)
            except:
                pass

            if self.mode == 'encut':
                self.conv_value_dir = os.path.join(convergence_dir, 'encut')
            elif self.mode == 'kgrid':
                self.conv_value_dir = os.path.join(convergence_dir, 'kgrid')
            try:
                os.mkdir(self.conv_value_dir)
            except:
                pass

            # write jobscript
            with open(f'{self.conv_value_dir}/jobscript','w') as f:
                f.write(self.jobheader)
                if '#SBATCH' in self.jobheader:
                    print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.mode} convergence calculation")
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.mode}")
                    f.write("\n")

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

    def write_input(self, params, name, rewrite=False):

        ### TODO: make a check if calculations folders are exist and calculations are finished successfully 

        # create an atoms object as a class atribute to make it accesible across the class methods and encode it

        # create working directory for manual submission
        if self.var_params:
            state_dir = os.path.join(self.conv_value_dir, f'{name}')
        else:
            self.atoms = read(self.poscar_file)
            calcfold = os.path.join(os.environ.get('AUTOMAG_PATH'), 'CalcFold')
            compound_dir = os.path.join(calcfold, f"{self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}")
            state_dir = os.path.join(compound_dir, f'{self.calculator}/{name}')

        try:
            os.mkdir(state_dir)
        except:
            pass

        def create_mae_dir(state_dir):
            # directory with all MAE calculations
            U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
            J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]),1)
            encut = params['encut']
            kpts = params['kpts']
            mae_dir = os.path.join(state_dir, f'mae_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}')
            try:
                os.mkdir(mae_dir)
            except:
                pass


            return mae_dir

        def check_z_dir():
            chgcar_found = False
            # checking the WAVECAR and CHGCAR in z_dir
            if os.path.isfile(f'{z_dir}/singlepoint/CHGCAR')==True and os.path.isfile(f'{z_dir}/singlepoint/WAVECAR')==True:
                if os.path.getsize(f'{z_dir}/singlepoint/CHGCAR')>0 and os.path.getsize(f'{z_dir}/singlepoint/WAVECAR')>0: 
                    print(f"Collinear CHGCAR and WAVECAR are found in {z_dir}/singlepoint directory")
                    chgcar_found = True
            # checking the WAVECAR and CHGCAR in state_dir/recalc
            elif os.path.isfile(f'{state_dir}/recalc/CHGCAR')==True and os.path.isfile(f'{state_dir}/recalc/WAVECAR')==True:
                if os.path.getsize(f'{state_dir}/recalc/CHGCAR')>0 and os.path.getsize(f'{state_dir}/recalc/WAVECAR')>0: 
                    print(f"Collinear CHGCAR and WAVECAR are found in {state_dir}/recalc directory")
                    os.system(f'cp -r {state_dir}/recalc {z_dir}/singlepoint')
                    chgcar_found = True
            # checking the WAVECAR and CHGCAR in state_dir/singelpoint
            elif os.path.isfile(f'{state_dir}/singlepoint/CHGCAR')==True and os.path.isfile(f'{state_dir}/singlepoint/WAVECAR')==True:
                if os.path.getsize(f'{state_dir}/singlepoint/CHGCAR')>0 and os.path.getsize(f'{state_dir}/singlepoint/WAVECAR')>0: 
                    print(f"Collinear CHGCAR and WAVECAR are found in {state_dir}/singlepoint directory")
                    os.system(f'cp -r {state_dir}/singlepoint {z_dir}/singlepoint')
                    chgcar_found = True

            else:
                print(f"Collinear CHGCAR and/or WAVECAR are not found in {state_dir}/recalc nor in {state_dir}/singlepoint directory")
                chgcar_found = False
                try:
                    os.mkdir(z_dir)
                except:
                    pass


            if chgcar_found == False:
                params['lnoncollinear'] = True
                params['lsorbit'] = True
                params['icharg'] = 2
                params['istart'] = 0
                params['lcharg'] = True
                params['lwave'] = True
                # writing input files
                self.write_vasp_input_files(z_dir,'singlepoint',params)
                # writing the jobscript one single point calculation without recalc
                self.write_jobscript(jobdir=z_dir,
                                     jobscript_name='jobscript_mae_theta_phi',
                                     mode='singlepoint',
                                     write_output=False,
                                     keep_large_files=True)

            return chgcar_found


        if self.mode == 'singlepoint':
            if name == 'nm':

                # write 
                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(state_dir,'singlepoint',params)

                # write jobscript
                self.write_jobscript(jobdir=state_dir,mode='singlepoint')

            else:

                if self.calculator == "vasp":
                    # single-point run + recalc run
                    self.write_vasp_input_files(state_dir,'singlepoint+recalc',params)

                # write jobscript
                self.write_jobscript(jobdir=state_dir,mode='singlepoint+recalc')


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

            # write the jobscript
            self.write_jobscript(jobdir=state_dir_pert,
                                 params=params,mode=self.mode)

            # restoring the name of the dummy atom POTCAR
            os.system(f'mv {dummy_atom_pp_path_temp} {dummy_atom_pp_path}')


        elif self.mode == 'encut' or self.mode == 'kgrid':

                # write input files
                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(state_dir,'singlepoint',params)

                # write the jobscript
                self.write_jobscript(jobdir=state_dir,mode=self.mode)


        elif self.mode == 'mae_theta_phi':

            ### Finding SAXIS that correspond to the lowest and highest energy

            mae_dir = create_mae_dir(state_dir)

            # directory with collinear calculation with SAXIS || z
            z_dir = os.path.join(mae_dir, 'z')
            try:
                os.mkdir(z_dir)
            except:
                pass


            if self.calculator == "vasp":
                chgcar_found = check_z_dir()




            phi = np.linspace(0.0001, 2 * np.pi-0.0001, self.Nph+1)
            theta = np.linspace(0+0.0001, np.pi-0.00001, self.Nth+1)
            THETA, PHI = np.meshgrid(theta, phi)

            # ncl_magmom = np.array([(0,0,m) for m in list(self.magmoms)])
            # self.magmoms = ncl_magmom

            params['voskown'] = 1
            params['lnoncollinear'] = True
            params['lsorbit'] = True
            params['icharg'] = 11
            params['istart'] = 1
            params['gga_compat'] = False
            params['lcharg'] = False
            params['lwave'] = False
            # params['magmom'] = ncl_magmom
            # self.set_magmoms(params)
            # self.atoms.set_initial_magnetic_moments(ncl_magmom)


            for phin, thetan in zip(PHI.flatten(),THETA.flatten()):
                # creating folder for calculation at fixed Phi and Theta
                phin_thetan_dir_name = 'PhTh_'+str(np.round((phin/np.pi)*180,2))+'_'+str(np.round((thetan/np.pi)*180,2))
                phin_thetan_dir = os.path.join(mae_dir, phin_thetan_dir_name)
                try:
                    os.mkdir(phin_thetan_dir)
                except:
                    pass

                X,Y,Z = RotatingAngles(thetan, phin)
                params['saxis'] = [np.round(X,3),np.round(Y,3),np.round(Z,3)]

                # write input files
                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(phin_thetan_dir,'singlepoint',params)

                # jobdir = phin_thetan_dir

                calcname = f'_MAE_{phin_thetan_dir_name}'

                self.write_jobscript(jobdir=phin_thetan_dir,
                                     jobdir_prev=f'{z_dir}/singlepoint',
                                     jobscript_name='jobscript_mae_theta_phi',
                                     calcname=f'_MAE_{phin_thetan_dir_name}',
                                     mode='nsc',write_output=False,keep_large_files=False)


        elif self.mode == 'mae_curve':

            mae_dir = create_mae_dir(state_dir)

            # directory with collinear calculation with SAXIS || z
            z_dir = os.path.join(mae_dir, 'z')
            try:
                os.mkdir(z_dir)
            except:
                pass

            N_MAE = self.N_MAE
            MAE_x = self.MAE_x
            MAE_z = self.MAE_z
            MAE_y = self.MAE_y

            params['voskown'] = 1
            params['lnoncollinear'] = True
            params['lsorbit'] = True
            params['icharg'] = 11
            params['istart'] = 1
            params['gga_compat'] = False
            params['lcharg'] = False
            params['lwave'] = False

            if self.calculator == "vasp":
                chgcar_found = check_z_dir()

                # E_ref_MAE = float((Oszicar(f"{z_dir}/singlepoint/OSZICAR").all_energies[-1])[-2])

            for alpha in  np.linspace(0, 2*np.pi, N_MAE+1):
                alpha_dir_name = 'K_'+str(params['kpts'])+'_RtMAE_'+str(np.round((alpha/np.pi)*180,2))
                alpha_dir = os.path.join(mae_dir, alpha_dir_name)
                try:
                    os.mkdir(alpha_dir)
                except:
                    pass

                SAX = np.sin(alpha)*MAE_y+np.cos(alpha)*MAE_x
                # searchExp = 'SAXIS'
                print('alpha: ', np.round(alpha*180/np.pi,3), ', SAXIS: ', np.round(SAX,3))
                params['saxis'] = np.round(SAX,3).tolist()

                # write input files
                if self.calculator == "vasp":
                    # only single-point run
                    self.write_vasp_input_files(alpha_dir,'singlepoint',params)

                # jobdir = phin_thetan_dir

                calcname = f'_MAE_{alpha_dir_name}'

                self.write_jobscript(jobdir=alpha_dir,
                                     jobdir_prev=f'{z_dir}/singlepoint',
                                     jobscript_name='jobscript_mae_curve',
                                     calcname=f'_MAE_curve_{alpha_dir_name}',
                                     mode='nsc',write_output=False,keep_large_files=False)

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


        if innermode == 'singlepoint' or innermode == 'recalc':
            workdir =  os.path.join(state_dir, innermode)

            try:
                os.mkdir(workdir)
            except:
                pass

            # creating ase.calculator object and writing input   
            calc = Vasp(atoms=self.atoms,
                        directory=workdir,
                        **params)

            calc.write_input(self.atoms)

            # create input script
            if innermode == 'singlepoint':

                if self.mode == 'singlepoint':
                    self.write_initial_magmoms(f'{workdir}/initial_magmoms.txt')

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


    def write_jobscript(self,jobdir,jobdir_prev=None,jobscript_name='jobscript',calcname='',params=None,mode='singlepoint',write_output=True,keep_large_files=False):

        """

        jobdir : directory for job
        jobdir_prev : directory with previous job 
        calcaname : name of calculation
        
        """

        # taking a higher level path to the jobdir
        job_top_dir = os.path.dirname(jobdir)
        # job_top_dir = os.path.dirname(os.path.dirname(jobdir))


        # writing a jobscript for singlepoint run (used in 2_coll stage)
        if mode == 'singlepoint' or mode == 'singlepoint+recalc':

            if self.parallel_over_configurations:
                # writing separate jobscripts for every configuration
                f = open(f'{jobdir}/{jobscript_name}','w')
                f.write(self.jobheader)
                if '#SBATCH' in self.jobheader:
                    print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.name} calculation")
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.name}")   



            else:
                # writing one jobscript for all configurations
                if os.path.exists(f'{job_top_dir}/{jobscript_name}'):
                    # open existing jobscript file for appending lines
                    f = open(f'{job_top_dir}/{jobscript_name}','a')
                else:
                    # create new jobscript file and writing header
                    f = open(f'{job_top_dir}/{jobscript_name}','w')
                    f.write(self.jobheader)
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{mode}")  

            f.write('\n')
            f.write('\n')
            f.write(f'cd {jobdir}/singlepoint \n')
            f.write('\n')
            f.write(self.calculator_command)
            f.write('\n')
            if not keep_large_files:
                f.write(f'rm WAVECAR CHG\n') # remove large files

            if mode == 'singlepoint+recalc':
                # do recalculation with final magnetic moments 
                f.write(f'cd {jobdir}/recalc\n')
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
            if not keep_large_files:
                f.write(f'rm WAVECAR CHG\n') # remove large files

            if write_output:
                f.write(f'cp {self.write_output_script_path} .\n')
                f.write(self.environment_activate)
                f.write('\n')
                if self.calculator == "vasp":
                    f.write(f'python {self.write_output_script_name} magmoms vasp None {self.struct_suffix}\n')
                f.write('\n')
                f.write(self.environment_deactivate)
                f.write('\n')

            f.write('\n')
            f.write('\n')
            f.close()

        elif mode == 'perturbations':

            # start writing common jobscript for singlepoint calculation and calculations for all perturbations
            with open(f'{job_top_dir}/{jobscript_name}','w') as f:
                f.write(self.jobheader)
                if '#SBATCH' in self.jobheader:
                    print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{mode}_{self.name} calculations")
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{mode}_{self.name}")
                f.write('\n')
                f.write('\n')
                f.write(f'cd {jobdir}/singlepoint\n')
                f.write('\n')
                f.write(self.calculator_command)
                f.write('\n')

            # write files for perturbation calculations at different pert values
            for perturbation in self.pert_values:
                print(f'perturbation: {perturbation}')
                pert_dir = os.path.join(jobdir, f'{perturbation}')

                self.write_vasp_input_files(pert_dir,'nsc',params,pert_value=perturbation)
                self.write_vasp_input_files(pert_dir,'sc',params,pert_value=perturbation)

                # appending the part to the jobscript
                with open(f'{job_top_dir}/{jobscript_name}','a') as f:
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

            if write_output:
            # os.system(f'cp {self.write_charges_script_path} {state_dir_pert}')
                with open(f'{job_top_dir}/{jobscript_name}','a') as f:
                    f.write(f'cd {job_top_dir}\n')
                    f.write(f'cp {self.write_charges_script_path} .\n')
                    f.write(self.environment_activate)
                    f.write('\n')
                    # if self.calculator == "vasp":
                    f.write(f'python {self.write_charges_script_name} {self.dummy_position} {self.calculator}\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)            
                    f.write('\n')

        elif mode == 'encut' or mode == 'kgrid':
            # add the lines to the jobscript
            with open(f'{self.conv_value_dir}/{jobscript_name}','a') as f:
                f.write(f'cd {jobdir}/singlepoint\n')
                f.write(self.calculator_command)
                f.write('\n')
                f.write(f'rm WAVECAR CHG\n') # remove large files
                f.write(f'cp {self.write_output_script_path} .\n')
                f.write(self.environment_activate)
                f.write('\n')
                if write_output:
                    if self.calculator == "vasp":
                        f.write(f'python {self.write_output_script_name} energy vasp {self.mode} {self.struct_suffix}\n')
                    f.write('\n')
                    f.write(self.environment_deactivate)
                f.write('\n')            

        elif mode == 'nsc':

            # write jobscript
            if self.parallel_over_configurations:
                # writing separate jobscripts for every configuration
                f = open(f'{jobdir}/{jobscript_name}','w')
                f.write(self.jobheader)
                if '#SBATCH' in self.jobheader:
                    print(f"Using SLURM queue system for {self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.name}{calcname} MAE calculation")
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_{self.name}{calcname}")                    
            else:
                # writing one jobscript for all configurations
                if os.path.exists(f'{job_top_dir}/{jobscript_name}'):
                    # open existing jobscript file for appending lines
                    f = open(f'{job_top_dir}/{jobscript_name}','a')
                else:
                    # create new jobscript file and writing header
                    f = open(f'{job_top_dir}/{jobscript_name}','w')
                    f.write(self.jobheader)
                    f.write(f"#SBATCH --job-name={self.atoms.get_chemical_formula(mode='metal')}{self.struct_suffix}_MAE")  

            f.write('\n')
            f.write('\n')
            f.write(f'cd {jobdir}/singlepoint \n')
            f.write('\n')
            f.write(f'cp {jobdir_prev}/CHGCAR ./CHGCAR \n')
            f.write(f'ln -s {jobdir_prev}/WAVECAR ./WAVECAR \n')
            f.write('\n')
            f.write(self.calculator_command)
            f.write('\n')
            if not keep_large_files:
                f.write(f'rm WAVECAR CHG\n') # remove large files
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.close()



        # write helping script for submitting the calculations through the job system
        if self.parallel_over_configurations:
            # appending a line to the submit script 
            script_name = f'{job_top_dir}/run_all_{mode}.sh'
            if not os.path.exists(script_name):
                with open(script_name,'w') as f1:
                    f1.write('#!/bin/sh\n')
                os.system(f'chmod +x {script_name}')
            if '#SBATCH' in self.jobheader:
                # using slurm scheduler
                exec_command = 'sbatch' 
            elif '#PBS' in self.jobheader:
                # using pbs scheduler
                exec_command = 'qsub'
            else:
                # using sh command
                exec_command = 'sh'

            submit_command = f'cd {jobdir}; {exec_command} jobscript\n'
            write_submit_script(script_name,submit_command)

        return

""" Helping functions """

def write_submit_script(script_name,submit_command):

    with open(script_name,'a') as f:
        f.write(submit_command)


""" MAE functions """

def RotatingAngles(mytheta, myphi):
    x= np.sin(mytheta)*np.cos(myphi)
    y= np.sin(mytheta)*np.sin(myphi)
    z= np.cos(mytheta)
    return x,y,z