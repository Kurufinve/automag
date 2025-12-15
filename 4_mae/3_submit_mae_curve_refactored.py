"""
Refactored MAE curve submission script.

Submits calculations along rotation path from easy to hard axis.

Run this after 2_analyze_results_refactored.py to get easy/hard axes.
"""

import os
import sys
import numpy as np
from ase.io import read

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
N_MAE = 20  # Number of points for MAE curve

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
    """Main execution function."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    calcfold_path = os.path.join(path_to_automag, 'CalcFold')
    
    # Check that easy/hard axes are defined
    if 'easy_axis' not in globals() or 'hard_axis' not in globals():
        print("ERROR: easy_axis and hard_axis must be defined in input.py")
        print("\nRun 2_analyze_results_refactored.py first to find these axes.")
        print("Then add them to input.py, for example:")
        print("  easy_axis = np.array([0.0, 0.0, 1.0])")
        print("  hard_axis = np.array([1.0, 0.0, 0.0])")
        return
    
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE CALCULATION - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    print(f"Easy axis: {easy_axis}")
    print(f"Hard axis: {hard_axis}")
    print(f"Number of points: {N_MAE + 1}")
    
    # Load processed structure (from 1_submit)
    try:
        # Try to find the processed structure file (new naming)
        import glob
        processed_files = glob.glob(f'*{configuration}_processed.vasp')
        standardized_files = glob.glob(f'*{configuration}_standardized.vasp')  # Legacy naming
        
        if processed_files:
            structure_path = processed_files[0]
        elif standardized_files:
            structure_path = standardized_files[0]
        else:
            print("ERROR: Structure file not found!")
            print(f"Expected: *{configuration}_processed.vasp or *{configuration}_standardized.vasp")
            print("Run 1_submit_refactored.py first.")
            return
    except Exception as e:
        print(f"ERROR loading structure: {e}")
        return
    
    atoms = read(structure_path)
    print(f"Structure loaded from: {structure_path}")
    
    # Load NCL magmoms from config file
    ncl_magmoms = None
    config_file = f'{configuration}_mae_config.txt'
    
    if os.path.exists(config_file):
        print(f"Reading configuration from: {config_file}")
        try:
            with open(config_file, 'r') as f:
                for line in f:
                    if 'NCL magmoms:' in line:
                        # Parse the magmoms - format is a list of tuples
                        magmoms_str = line.split('NCL magmoms:')[1].strip()
                        # Use eval to parse the list of tuples
                        # Format: [(0, 0, m1), (0, 0, m2), ...]
                        ncl_magmoms = eval(magmoms_str)
                        print(f"Loaded {len(ncl_magmoms)} NCL magnetic moments from config")
                        break
        except Exception as e:
            print(f"Warning: Could not parse NCL magmoms from config: {e}")
            print(f"Will try to use ncl_magmoms from input.py if available")
    else:
        print(f"Warning: Config file {config_file} not found")
        print(f"Will try to use ncl_magmoms from input.py if available")
    
    # Fall back to global scope if not loaded from config
    if ncl_magmoms is None:
        if 'ncl_magmoms' in globals():
            ncl_magmoms = globals()['ncl_magmoms']
            print(f"Using ncl_magmoms from input.py: {len(ncl_magmoms)} values")
        else:
            print("ERROR: ncl_magmoms not found!")
            print("\nncl_magmoms should be either:")
            print("  1. In the config file (created by 1_submit_refactored.py)")
            print(f"  2. Defined in input.py as: ncl_magmoms = [(0, 0, m1), (0, 0, m2), ...]")
            print("\nPlease run 1_submit_refactored.py first, or define ncl_magmoms manually.")
            return
    
    # Create base parameters for non-collinear calculations
    base_params_dict = params.copy()
    base_params_dict.update({
        'voskown': 1,
        'lnoncollinear': True,
        'lsorbit': True,
        'gga_compat': False,
    })
    
    base_params = CalculationParameters(**base_params_dict)
    base_params.validate()
    
    # Create workflow submitter
    if use_fireworks:
        launchpad_file = os.path.join(os.path.expanduser('~'), '.fireworks/my_launchpad.yaml')
        submitter = WorkflowSubmitterFactory.create_fireworks_submitter(launchpad_file)
    else:
        submitter = WorkflowSubmitterFactory.create_manual_submitter(
            calcfold_path=calcfold_path,
            jobheader=jobheader,
            calculator_command=calculator_command,
            environment_activate=environment_activate,
            environment_deactivate=environment_deactivate
        )
    
    # Create service
    service = CalculationService(submitter)
    
    # Submit MAE curve calculations
    print(f"\nSubmitting MAE curve calculations...")
    
    job_ids = service.submit_mae_curve(
        atoms=atoms,
        base_params=base_params,
        magmoms=ncl_magmoms,
        n_points=N_MAE,
        easy_axis=easy_axis,
        hard_axis=hard_axis,
        workflow_name=f'mae_curve_{configuration}'
    )
    
    print(f"\n✓ Submitted {len(job_ids)} MAE curve calculations")
    print(f"\nAfter calculations complete, run 4_plot_mae_curve_refactored.py")


if __name__ == '__main__':
    main()
