"""
Refactored MAE (Magnetocrystalline Anisotropy Energy) submission script.

This demonstrates the SOLID-compliant architecture for MAE calculations.
The old code is preserved in 1_submit.py and MAE.py for reference.

Follows:
- Single Responsibility Principle: Separate classes for MAE analysis
- Dependency Inversion Principle: Depends on abstractions
- Open/Closed Principle: Easy to extend with new MAE analysis methods
"""

import os
import sys
import numpy as np
from pathlib import Path

from ase.io import read
from pymatgen.core.structure import Structure

# Import refactored components
from core.factories.workflow_factory import WorkflowSubmitterFactory
from core.services.calculation_service import CalculationService
from core.services.mae_service import (
    MAEDirectionGenerator,
    MAEAnalyzer,
    MAEResultsLoader,
    MAEPlotter
)
from core.domain.calculation import CalculationParameters

# Default values
use_fireworks = False
calculator = 'vasp'
jobheader = """#!/bin/bash"""
calculator_command = "mpirun vasp_std"
environment_activate = "source .venv/bin/activate"
environment_deactivate = "deactivate"
parallel_over_configurations = True
struct_suffix = ''

# MAE calculation parameters
Nph = 20  # Number of phi points
Nth = 10  # Number of theta points
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


def load_configuration_from_collinear(path_to_automag: str,
                                     formula: str,
                                     configuration_name: str,
                                     struct_suffix: str = '',
                                     calculator_type: str = 'vasp') -> tuple:
    """
    Load magnetic configuration from collinear calculations.
    
    Returns:
        Tuple of (structure_path, magmoms, setting_number)
    """
    # Find collinear results
    path_to_coll = Path(path_to_automag) / '2_coll' / f'{formula}{struct_suffix}' / f'{formula}_{calculator_type}'
    if not path_to_coll.exists():
        path_to_coll = Path(path_to_automag) / '2_coll' / f'{formula}_{calculator_type}'
        if not path_to_coll.exists():
            path_to_coll = Path(path_to_automag) / '2_coll'
    
    print(f'Loading configuration from: {path_to_coll}')
    
    path_to_trials = path_to_coll / f'trials_{formula}'
    
    # Find configuration in trials
    setting = 1
    magmom = None
    
    while (path_to_trials / f'configurations{setting:03d}.txt').exists():
        with open(path_to_trials / f'configurations{setting:03d}.txt', 'r') as f:
            for line in f:
                values = line.split()
                if values[0] == configuration_name:
                    magmom = [int(item) for item in values[1:]]
                    break
            if magmom is not None:
                break
        setting += 1
    
    if magmom is None:
        raise IOError(f'Configuration {configuration_name} not found in {path_to_trials}')
    
    print(f'Found configuration {configuration_name} in setting{setting:03d}')
    print(f'Initial magmoms: {magmom}')
    
    # Load structure
    structure_path = path_to_trials / f'setting{setting:03d}.vasp'
    
    # Load final magnetic moments from calculation results
    calc_results_file = Path(path_to_automag) / 'CalcFold' / f"{formula}{struct_suffix}_singlepoint_{calculator_type}.txt"
    
    full_magmom = None
    print(f'Reading final magmoms from: {calc_results_file}')
    
    with open(calc_results_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            values = line.split()
            if values[0] == configuration_name:
                print(f'Found configuration: {values[0]}')
                magmoms_str = lines[i + 1]
                full_magmom_str = magmoms_str.split('final_magmoms=')[1].strip().strip('[]')
                number_magmoms = full_magmom_str.split()
                full_magmom = [float(num) for num in number_magmoms]
                break
    
    if full_magmom is None:
        raise IOError(f'Final magmoms for {configuration_name} not found')
    
    print(f'Final magnetic moments: {full_magmom}')
    
    return str(structure_path), full_magmom, setting


def standardize_structure(structure_path: str, 
                         magmoms: list,
                         output_name: str) -> tuple:
    """
    Standardize structure to primitive cell.
    
    Returns:
        Tuple of (standardized_structure_path, ncl_magmoms)
    """
    # Load pymatgen structure
    pmg_structure = Structure.from_file(structure_path)
    
    # Add magnetic moments as site property
    pmg_structure.add_site_property("magmom", magmoms)
    
    # Get primitive cell
    standardized_structure = pmg_structure.get_primitive_structure(
        tolerance=0.2,
        use_site_props=True
    )
    
    # Save standardized structure
    standardized_structure.to(filename=output_name, fmt='POSCAR')
    
    # Get standardized magnetic moments (collinear)
    standardized_magmom = list(standardized_structure.site_properties['magmom'])
    
    # Convert to non-collinear format (0, 0, m)
    ncl_magmom = [(0, 0, m) for m in standardized_magmom]
    
    print(f'Standardized structure saved to: {output_name}')
    print(f'Standardized magmoms (NCL): {ncl_magmom}')
    
    return output_name, ncl_magmom


def main():
    """Main execution function following SOLID principles."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    calcfold_path = os.path.join(path_to_automag, 'CalcFold')
    
    # Load input structure to get formula
    input_structure = Structure.from_file(path_to_poscar)
    formula = input_structure.formula.replace(' ', '')
    
    print(f"\n{'=' * 70}")
    print(f"MAE CALCULATION FOR {formula} - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    
    # Load configuration from collinear calculations
    structure_path, final_magmoms, setting = load_configuration_from_collinear(
        path_to_automag, formula, configuration, struct_suffix, calculator
    )
    
    # Standardize structure to primitive cell
    standardized_file = f'setting{setting:03d}_{configuration}_standardized.vasp'
    standardized_path, ncl_magmoms = standardize_structure(
        structure_path, final_magmoms, standardized_file
    )
    
    # Load standardized structure
    atoms = read(standardized_path)
    
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
        print("WARNING: FireWorks mode - MAE calculations with FireWorks not fully tested!")
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
    
    # Submit MAE theta-phi grid calculations
    print(f"\n{'=' * 70}")
    print(f"SUBMITTING MAE THETA-PHI GRID CALCULATIONS")
    print(f"Grid size: {Nth} theta points × {Nph} phi points = {(Nth+1)*(Nph+1)} calculations")
    print(f"{'=' * 70}\n")
    
    job_ids = service.submit_mae_theta_phi_grid(
        atoms=atoms,
        base_params=base_params,
        magmoms=ncl_magmoms,
        n_theta=Nth,
        n_phi=Nph,
        reference_dir='z',
        workflow_name=f'mae_grid_{configuration}'
    )
    
    print(f"✓ Submitted {len(job_ids)} MAE grid calculations")
    
    # Note about MAE curve calculations
    print(f"\n{'=' * 70}")
    print("MAE CURVE CALCULATION (Post-processing)")
    print(f"{'=' * 70}")
    print("After grid calculations complete:")
    print("1. Run 2_analyze_results.py to find easy/hard axes")
    print("2. Run 3_submit_mae_curve.py to calculate MAE along rotation path")
    print(f"{'=' * 70}\n")
    
    # Save configuration info
    with open(f'{configuration}_mae_config.txt', 'w') as f:
        f.write(f"Configuration: {configuration}\n")
        f.write(f"Structure: {standardized_file}\n")
        f.write(f"Formula: {formula}\n")
        f.write(f"Grid size: {Nth}×{Nph}\n")
        f.write(f"NCL magmoms: {ncl_magmoms}\n")
        f.write(f"Number of jobs: {len(job_ids)}\n")
    
    print(f"Configuration saved to: {configuration}_mae_config.txt")


if __name__ == '__main__':
    main()
