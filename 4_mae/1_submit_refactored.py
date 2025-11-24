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


def generate_submission_script(mae_dir: Path,
                               calc_dirs: list,
                               directions: list,
                               parallel: bool,
                               jobheader: str,
                               calculator_command: str,
                               environment_activate: str,
                               environment_deactivate: str):
    """
    Generate SLURM submission helper script.
    
    Args:
        mae_dir: Base MAE directory
        calc_dirs: List of calculation directories
        directions: List of MAE directions
        parallel: If True, create parallel submission script; if False, sequential
        jobheader: SLURM job header template
        calculator_command: Command to run VASP
        environment_activate: Command to activate environment
        environment_deactivate: Command to deactivate environment
    """
    
    if parallel:
        # Create parallel submission script (one click to submit all jobs)
        script_path = mae_dir / 'submit_mae_grid.sh'
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# MAE Grid Parallel Submission Script\n")
            f.write("# This script submits all MAE grid calculations in parallel\n\n")
            
            f.write(f"echo 'Submitting {len(calc_dirs)} MAE grid calculations in parallel...'\n\n")
            
            for calc_dir, direction in zip(calc_dirs, directions):
                job_name = f"mae_{direction.name}"
                rel_dir = calc_dir.relative_to(mae_dir)
                
                f.write(f"# Submit {direction.name}\n")
                f.write(f"cd {calc_dir}\n")
                f.write(f"cat > job.sh << 'EOF'\n")
                f.write(jobheader + "\n")
                f.write(f"#SBATCH -J {job_name}\n")
                f.write(f"#SBATCH -o {calc_dir}/slurm-%j.out\n")
                f.write(f"#SBATCH -e {calc_dir}/slurm-%j.err\n\n")
                f.write(f"{environment_activate}\n")
                f.write(f"{calculator_command}\n")
                f.write(f"{environment_deactivate}\n")
                f.write("EOF\n")
                f.write("sbatch job.sh\n")
                f.write(f"cd {mae_dir}\n\n")
            
            f.write("echo 'All jobs submitted!'\n")
            f.write("echo 'Monitor with: squeue -u $USER'\n")
        
        script_path.chmod(0o755)
        print(f"Created parallel submission script: {script_path}")
    
    else:
        # Create sequential submission script
        script_path = mae_dir / 'submit_mae_grid_sequential.sh'
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# MAE Grid Sequential Submission Script\n")
            f.write("# This script runs MAE grid calculations sequentially in a single SLURM job\n\n")
            
            f.write(jobheader + "\n")
            f.write("#SBATCH -J mae_grid_sequential\n")
            f.write(f"#SBATCH -o {mae_dir}/mae_sequential-%j.out\n")
            f.write(f"#SBATCH -e {mae_dir}/mae_sequential-%j.err\n\n")
            
            f.write(f"{environment_activate}\n\n")
            
            for calc_dir, direction in zip(calc_dirs, directions):
                f.write(f"echo 'Running calculation for {direction.name}...'\n")
                f.write(f"cd {calc_dir}\n")
                f.write(f"{calculator_command}\n")
                f.write("\n")
                f.write(f"# Check if calculation completed successfully\n")
                f.write(f"if [ $? -ne 0 ]; then\n")
                f.write(f"    echo 'ERROR: Calculation failed for {direction.name}'\n")
                f.write(f"    exit 1\n")
                f.write(f"fi\n")
                f.write(f"echo 'Completed {direction.name}'\n\n")
            
            f.write(f"{environment_deactivate}\n")
            f.write("echo 'All MAE grid calculations completed!'\n")
        
        script_path.chmod(0o755)
        print(f"Created sequential submission script: {script_path}")


def main():
    """Main execution function following SOLID principles."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    
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
    
    # Get U and J values for directory naming
    ldauu_val = params.get('ldauu', [0.0])
    ldauj_val = params.get('ldauj', [0.0])
    U = ldauu_val[next(i for i, x in enumerate(params.get('ldaul', [])) if x > 0)] if params.get('ldaul') else 0.0
    J = ldauj_val[next(i for i, x in enumerate(params.get('ldaul', [])) if x > 0)] if params.get('ldaul') else 0.0
    
    # Handle kpts and encut as lists or single values
    kpts_list = params['kpts'] if isinstance(params['kpts'], list) else [params['kpts']]
    encut_list = params['encut'] if isinstance(params['encut'], list) else [params['encut']]
    
    # Iterate over kpts and encut combinations
    for kpts_val in kpts_list:
        for encut_val in encut_list:
            # Create MAE calculation directory
            calcfold_path = Path(path_to_automag) / 'CalcFold'
            mae_base_dir = calcfold_path / f"{formula}{struct_suffix}" / calculator / configuration
            mae_dir = mae_base_dir / f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}"
            mae_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"\nCreating MAE directory: {mae_dir}")
            
            # Create base parameters for non-collinear calculations
            base_params_dict = params.copy()
            base_params_dict.update({
                'kpts': kpts_val,
                'encut': encut_val,
                'voskown': 1,
                'lnoncollinear': True,
                'lsorbit': True,
                'gga_compat': False,
            })
            
            base_params = CalculationParameters(**base_params_dict)
            base_params.validate()
            
            # Generate MAE directions
            from core.services.mae_service import MAEDirectionGenerator
            directions = MAEDirectionGenerator.generate_theta_phi_grid(Nth, Nph)
            
            # Create calculation directories and input files
            calc_dirs = []
            for direction in directions:
                calc_dir = mae_dir / direction.name
                calc_dir.mkdir(parents=True, exist_ok=True)
                calc_dirs.append(calc_dir)
                
                # Copy structure
                import shutil
                shutil.copy(standardized_path, calc_dir / 'POSCAR')
                
                # Create ASE atoms and write VASP inputs
                from ase.io import read as ase_read
                atoms_calc = ase_read(standardized_path)
                
                # Set non-collinear magnetic moments
                # For ASE with non-collinear, we need to flatten the tuples to a 1D array
                # ncl_magmoms is like [(0, 0, 5.0), (0, 0, 5.0), ...]
                # We need [0, 0, 5.0, 0, 0, 5.0, ...]
                magmom_flat = []
                for mx, my, mz in ncl_magmoms:
                    magmom_flat.extend([mx, my, mz])
                atoms_calc.set_initial_magnetic_moments(magmom_flat)
                
                # Prepare VASP parameters
                params_with_saxis = base_params_dict.copy()
                params_with_saxis['saxis'] = list(direction.saxis)
                
                # For grid calculations after reference, read WAVECAR and CHGCAR
                if direction != directions[0]:  # Not the first (reference) calculation
                    params_with_saxis['icharg'] = 11
                    params_with_saxis['istart'] = 1
                    params_with_saxis['lcharg'] = False
                    params_with_saxis['lwave'] = False
                
                # Write VASP input files using ASE
                from ase.calculators.vasp import Vasp
                calc = Vasp(**params_with_saxis)
                calc.write_input(atoms_calc, directory=str(calc_dir))
            
            # Generate helper submission script
            generate_submission_script(
                mae_dir, calc_dirs, directions, parallel_over_configurations,
                jobheader, calculator_command, environment_activate, environment_deactivate
            )
            
            print(f"✓ Created {len(calc_dirs)} calculation directories in {mae_dir}")
            print(f"✓ Generated submission script: {mae_dir / 'submit_mae_grid.sh'}")
    
    # Save configuration info
    with open(f'{configuration}_mae_config.txt', 'w') as f:
        f.write(f"Configuration: {configuration}\n")
        f.write(f"Structure: {standardized_file}\n")
        f.write(f"Formula: {formula}\n")
        f.write(f"Grid size: {Nth}×{Nph}\n")
        f.write(f"NCL magmoms: {ncl_magmoms}\n")
        f.write(f"Kpts values: {kpts_list}\n")
        f.write(f"Encut values: {encut_list}\n")
        f.write(f"U = {U:.1f}, J = {J:.1f}\n")
    
    print(f"\n✓ Configuration saved to: {configuration}_mae_config.txt")
    print(f"\n{'=' * 70}")
    print("To submit calculations, run the generated script:")
    print(f"  bash <mae_directory>/submit_mae_grid.sh")
    print(f"{'=' * 70}\n")


if __name__ == '__main__':
    main()
