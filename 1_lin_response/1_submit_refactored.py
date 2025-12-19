"""
Refactored linear response U calculation submission script.

This demonstrates the SOLID-compliant architecture for perturbation calculations.
The old code is preserved in 1_submit.py for reference.

Follows:
- Single Responsibility Principle: Clear separation of concerns
- Dependency Inversion Principle: Depends on abstractions via factory
- Open/Closed Principle: Can switch workflow systems without code changes
"""

import os
import sys
from ase.io import read
from pymatgen.core.structure import Structure

# Import refactored components
from core.factories.workflow_factory import WorkflowSubmitterFactory
from core.services.calculation_service import CalculationService
from core.domain.calculation import CalculationParameters

# Default values
use_fireworks = False
calculator = 'vasp'
jobheader = """#!/bin/bash"""
calculator_command = "mpirun vasp_std"
environment_activate = "source .venv/bin/activate"
environment_deactivate = "deactivate"
struct_suffix = ''

# Get current directory
cwd = os.getcwd()

# Load input configuration
try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    exec(f"from {input_file.split('.')[0]} import *")
except IndexError:
    print(f'Using the input.py file from folder: {cwd}')
    from input import *


def main():
    """Main execution function following SOLID principles."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    calcfold_path = os.path.join(path_to_automag, 'CalcFold')
    
    # Load structure
    atoms = read(path_to_poscar)
    
    # Substitute dummy atom
    ch_symbols = atoms.get_chemical_symbols()
    atom_ucalc = ch_symbols[dummy_position]
    ch_symbols[dummy_position] = dummy_atom
    atoms.set_chemical_symbols(ch_symbols)
    
    print(f'Dummy atom: {dummy_atom}')
    print(f'Original atom at position {dummy_position}: {atom_ucalc}')
    
    # Create magnetic configuration if not specified
    if 'configuration' not in globals():
        configuration = []
        structure = Structure.from_file(path_to_poscar)
        for atom in structure.species:
            if 'magnetic_atoms' not in globals():
                if atom.is_transition_metal:
                    configuration.append(4.0)
                else:
                    configuration.append(0.0)
            else:
                if atom in magnetic_atoms:
                    configuration.append(4.0)
                else:
                    configuration.append(0.0)
    
    # Create calculation parameters (validated, immutable)
    base_params = CalculationParameters(**params)
    base_params.validate()
    
    # Swap POTCAR files (dummy atom uses original atom's pseudopotential)
    _swap_potcar_files(dummy_atom, atom_ucalc, params.get('xc', 'PBE'), 
                       params.get('setups', {}))
    
    try:
        # Create workflow submitter using factory (Dependency Injection)
        if use_fireworks:
            launchpad_file = os.path.join(
                os.path.expanduser('~'),
                '.fireworks/my_launchpad.yaml'
            )
            submitter = WorkflowSubmitterFactory.create_fireworks_submitter(launchpad_file)
        else:
            # For linear response, use input_cell by default
            cell_type = 'input_cell'
            
            submitter = WorkflowSubmitterFactory.create_manual_submitter(
                calcfold_path=calcfold_path,
                jobheader=jobheader,
                calculator_command=calculator_command,
                environment_activate=environment_activate,
                environment_deactivate=environment_deactivate,
                cell_type=cell_type
            )
        
        # Create service (Single Responsibility)
        service = CalculationService(submitter)
        
        # Submit perturbation workflow
        print(f"Submitting perturbation workflow with {len(perturbations)} perturbation values...")
        
        workflow_id = service.submit_perturbation_workflow(
            atoms=atoms,
            base_params=base_params,
            magmoms=configuration,
            perturbations=perturbations,
            dummy_atom=dummy_atom,
            dummy_position=dummy_position,
            atom_ucalc=atom_ucalc,
            workflow_name=f'perturbations_{dummy_atom}'
        )
        
        print(f"Submitted workflow: {workflow_id}")
        print(f"Perturbation values: {perturbations}")
        
    finally:
        # Restore original POTCAR file
        _restore_potcar_files(dummy_atom, params.get('xc', 'PBE'))


def _swap_potcar_files(dummy_atom: str, atom_ucalc: str, xc: str, setups: dict) -> None:
    """
    Temporarily swap POTCAR file for dummy atom with original atom's POTCAR.
    
    This is necessary for linear response calculations where we want to treat
    one atom independently while maintaining its electronic properties.
    """
    vasp_pp_path = os.environ.get('VASP_PP_PATH')
    
    # Determine pseudopotential directory
    if xc == 'PBE':
        pp_path = os.path.join(vasp_pp_path, 'potpaw_PBE')
    elif xc == 'LDA':
        pp_path = os.path.join(vasp_pp_path, 'potpaw_LDA')
    else:
        pp_path = os.path.join(vasp_pp_path, 'potpaw_PBE')
    
    # Check for special POTCAR versions (e.g., '_sv', '_pv')
    if isinstance(setups, dict) and atom_ucalc in setups:
        potpaw_name = atom_ucalc + setups[atom_ucalc]
    else:
        potpaw_name = atom_ucalc
    
    atom_ucalc_pp = os.path.join(pp_path, potpaw_name, 'POTCAR')
    dummy_atom_pp = os.path.join(pp_path, dummy_atom, 'POTCAR')
    dummy_atom_pp_backup = os.path.join(pp_path, dummy_atom, '_POTCAR')
    
    print(f'Swapping POTCAR: {dummy_atom} <- {potpaw_name}')
    
    # Backup dummy atom POTCAR
    if os.path.exists(dummy_atom_pp):
        os.system(f'mv {dummy_atom_pp} {dummy_atom_pp_backup}')
    
    # Copy original atom POTCAR to dummy atom location
    os.system(f'cp {atom_ucalc_pp} {dummy_atom_pp}')


def _restore_potcar_files(dummy_atom: str, xc: str) -> None:
    """Restore original POTCAR file for dummy atom."""
    vasp_pp_path = os.environ.get('VASP_PP_PATH')
    
    if xc == 'PBE':
        pp_path = os.path.join(vasp_pp_path, 'potpaw_PBE')
    elif xc == 'LDA':
        pp_path = os.path.join(vasp_pp_path, 'potpaw_LDA')
    else:
        pp_path = os.path.join(vasp_pp_path, 'potpaw_PBE')
    
    dummy_atom_pp = os.path.join(pp_path, dummy_atom, 'POTCAR')
    dummy_atom_pp_backup = os.path.join(pp_path, dummy_atom, '_POTCAR')
    
    # Restore original POTCAR
    if os.path.exists(dummy_atom_pp_backup):
        os.system(f'mv {dummy_atom_pp_backup} {dummy_atom_pp}')
        print(f'Restored original POTCAR for {dummy_atom}')


if __name__ == '__main__':
    main()
