"""
Example refactored convergence test submission script.

This demonstrates how to use the new SOLID-compliant architecture.
The old code is preserved in 1_submit.py.old for reference.
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
    
    # Create calculation parameters (immutable, validated)
    base_params = CalculationParameters(**params)
    base_params.validate()
    
    # Create workflow submitter using factory (Dependency Injection)
    if use_fireworks:
        launchpad_file = os.path.join(
            os.path.expanduser('~'),
            '.fireworks/my_launchpad.yaml'
        )
        submitter = WorkflowSubmitterFactory.create_fireworks_submitter(launchpad_file)
    else:
        # For convergence tests, use input_cell by default
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
    
    # Submit convergence tests based on mode
    if mode == 'encut':
        # Set default values if not specified
        if 'encut_values' not in globals():
            encut_values = range(500, 1010, 10)
        
        # Submit calculations
        print(f"Submitting {len(encut_values)} ENCUT convergence tests...")
        job_ids = service.submit_convergence_test(
            atoms=atoms,
            base_params=base_params,
            magmoms=configuration,
            test_parameter='encut',
            test_values=list(encut_values),
            mode_name='encut'
        )
        print(f"Submitted {len(job_ids)} calculations")
        
    elif mode == 'kgrid':
        # Set default values if not specified
        if 'sigma_values' not in globals():
            sigma_values = [item / 100 for item in range(5, 25, 5)]
        if 'kpts_values' not in globals():
            kpts_values = range(20, 110, 10)
        
        # Submit calculations for each combination
        print(f"Submitting k-grid convergence tests...")
        total_jobs = 0
        
        for sigma in sigma_values:
            for kpts in kpts_values:
                # Create parameters with both sigma and kpts
                test_params = CalculationParameters(**params)
                test_params_dict = test_params.to_dict()
                test_params_dict['sigma'] = sigma
                test_params_dict['kpts'] = kpts
                
                test_params_obj = CalculationParameters(**test_params_dict)
                
                # Submit single calculation
                from core.interfaces.workflow import CalculationConfig
                config = CalculationConfig(
                    calc_params=test_params_obj.to_dict(),
                    magmoms=configuration,
                    name=f"kgrid{sigma}-{kpts}"
                )
                job_id = submitter.submit_single_calculation(atoms, config)
                total_jobs += 1
        
        print(f"Submitted {total_jobs} calculations")


if __name__ == '__main__':
    main()
