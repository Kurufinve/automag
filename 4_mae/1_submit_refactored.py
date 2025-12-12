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
                         output_name: str,
                         to_primitive: bool = False) -> tuple:
    """
    Standardize structure (optionally to primitive cell).
    
    Args:
        structure_path: Path to input structure
        magmoms: Magnetic moments
        output_name: Output file name
        to_primitive: If True, convert to primitive cell (default: False)
    
    Returns:
        Tuple of (standardized_structure_path, ncl_magmoms)
    """
    # Load pymatgen structure
    pmg_structure = Structure.from_file(structure_path)
    
    # Add magnetic moments as site property
    pmg_structure.add_site_property("magmom", magmoms)
    
    if to_primitive:
        # Get primitive cell
        standardized_structure = pmg_structure.get_primitive_structure(
            tolerance=0.2,
            use_site_props=True
        )
        print(f'Converted to primitive cell: {len(pmg_structure)} -> {len(standardized_structure)} atoms')
    else:
        # Keep original cell
        standardized_structure = pmg_structure
        print(f'Using original cell: {len(standardized_structure)} atoms')
    
    # Save standardized structure
    standardized_structure.to(filename=output_name, fmt='POSCAR')
    
    # Get standardized magnetic moments (collinear)
    standardized_magmom = list(standardized_structure.site_properties['magmom'])
    
    # Convert to non-collinear format (0, 0, m)
    ncl_magmom = [(0, 0, m) for m in standardized_magmom]
    
    print(f'Standardized structure saved to: {output_name}')
    print(f'Standardized magmoms (NCL): {len(ncl_magmom)} atoms')
    
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
            
            f.write("# IMPORTANT: Run reference calculation first!\n")
            f.write("echo 'Step 1: Submitting reference calculation (z folder)...'\n\n")
            
            # Submit reference calculation first
            z_dir = mae_dir / 'z'
            f.write(f"cd {z_dir}\n")
            f.write(f"cat > job.sh << 'EOF'\n")
            f.write(jobheader + "\n")
            f.write(f"#SBATCH -J mae_ref_z\n")
            f.write(f"#SBATCH -o {z_dir}/slurm-%j.out\n")
            f.write(f"#SBATCH -e {z_dir}/slurm-%j.err\n\n")
            f.write(f"{environment_activate}\n")
            f.write(f"{calculator_command}\n")
            f.write(f"{environment_deactivate}\n")
            f.write("EOF\n")
            f.write("REF_JOB=$(sbatch job.sh | awk '{print $4}')\n")
            f.write("echo \"Reference job submitted: $REF_JOB\"\n")
            f.write(f"cd {mae_dir}\n\n")
            
            f.write(f"echo 'Step 2: Creating symlinks to WAVECAR/CHGCAR from z folder...'\n\n")
            
            # Create symlinks in all grid directories
            for calc_dir in calc_dirs:
                f.write(f"cd {calc_dir}\n")
                f.write(f"ln -sf ../z/WAVECAR ./WAVECAR\n")
                f.write(f"ln -sf ../z/CHGCAR ./CHGCAR\n")
                f.write(f"cd {mae_dir}\n")
            
            f.write(f"\necho 'Step 3: Submitting {len(calc_dirs)} MAE grid calculations...'\n")
            f.write(f"echo 'Grid calculations will wait for reference job to complete'\n\n")
            
            for calc_dir, direction in zip(calc_dirs, directions):
                job_name = f"mae_{direction.name}"
                
                f.write(f"# Submit {direction.name}\n")
                f.write(f"cd {calc_dir}\n")
                f.write(f"cat > job.sh << 'EOF'\n")
                f.write(jobheader + "\n")
                f.write(f"#SBATCH -J {job_name}\n")
                f.write(f"#SBATCH -o {calc_dir}/slurm-%j.out\n")
                f.write(f"#SBATCH -e {calc_dir}/slurm-%j.err\n")
                f.write(f"#SBATCH --dependency=afterok:$REF_JOB\n\n")  # Wait for reference
                f.write(f"{environment_activate}\n")
                f.write(f"{calculator_command}\n")
                f.write(f"{environment_deactivate}\n")
                f.write("EOF\n")
                f.write("sbatch job.sh\n")
                f.write(f"cd {mae_dir}\n\n")
            
            f.write("echo 'All jobs submitted!'\n")
            f.write("echo 'Reference job: $REF_JOB'\n")
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
            
            # First run reference calculation
            z_dir = mae_dir / 'z'
            f.write(f"echo 'Step 1: Running reference calculation (z folder)...'\n")
            f.write(f"cd {z_dir}\n")
            f.write(f"{calculator_command}\n")
            f.write("\n")
            f.write(f"if [ $? -ne 0 ]; then\n")
            f.write(f"    echo 'ERROR: Reference calculation failed in z folder'\n")
            f.write(f"    exit 1\n")
            f.write(f"fi\n")
            f.write(f"echo 'Reference calculation completed'\n\n")
            
            # Create symlinks
            f.write(f"echo 'Step 2: Creating symlinks to WAVECAR/CHGCAR...'\n")
            for calc_dir in calc_dirs:
                f.write(f"cd {calc_dir}\n")
                f.write(f"ln -sf ../z/WAVECAR ./WAVECAR\n")
                f.write(f"ln -sf ../z/CHGCAR ./CHGCAR\n")
            f.write("echo 'Symlinks created'\n\n")
            
            # Run grid calculations
            f.write(f"echo 'Step 3: Running {len(calc_dirs)} grid calculations...'\n\n")
            
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
    
    # Standardize structure (use original cell to preserve atom count)
    standardized_file = f'setting{setting:03d}_{configuration}_standardized.vasp'
    to_primitive = False  # Set to True if you want primitive cell
    standardized_path, ncl_magmoms = standardize_structure(
        structure_path, final_magmoms, standardized_file, to_primitive
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
            
            # First, create reference calculation (z folder) - self-consistent
            print(f"\nCreating reference self-consistent calculation in 'z' folder...")
            z_dir = mae_dir / 'z'
            z_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy structure to z folder
            import shutil
            shutil.copy(standardized_path, z_dir / 'POSCAR')
            
            # Create ASE atoms for reference calculation
            from ase.io import read as ase_read
            atoms_ref = ase_read(standardized_path)
            
            # Reference calculation parameters (self-consistent, no SAXIS rotation yet)
            params_ref = base_params_dict.copy()
            params_ref['saxis'] = [0, 0, 1]  # Reference direction (z-axis)
            params_ref.pop('magmom', None)
            
            # Keep default ICHARG for self-consistent calculation
            params_ref.pop('icharg', None)
            params_ref.pop('istart', None)
            
            # Ensure CHGCAR and WAVECAR are saved
            params_ref['lcharg'] = True
            params_ref['lwave'] = True
            
            # Write reference calculation input files
            from ase.calculators.vasp import Vasp
            calc_ref = Vasp(directory=str(z_dir), **params_ref)
            calc_ref.write_input(atoms_ref)
            
            # Add non-collinear MAGMOM to reference INCAR
            incar_ref_path = z_dir / 'INCAR'
            with open(incar_ref_path, 'r') as f:
                incar_ref_lines = f.readlines()
            
            magmom_added = False
            for i, line in enumerate(incar_ref_lines):
                if line.strip().startswith('MAGMOM'):
                    magmom_str = '  '.join([f'{mx} {my} {mz}' for mx, my, mz in ncl_magmoms])
                    incar_ref_lines[i] = f'MAGMOM = {magmom_str}\n'
                    magmom_added = True
                    break
            
            if not magmom_added:
                magmom_str = '  '.join([f'{mx} {my} {mz}' for mx, my, mz in ncl_magmoms])
                incar_ref_lines.append(f'MAGMOM = {magmom_str}\n')
            
            with open(incar_ref_path, 'w') as f:
                f.writelines(incar_ref_lines)
            
            print(f"✓ Created reference calculation in {z_dir}")
            
            # Now create grid calculations (non-self-consistent, reading from z)
            print(f"\nCreating MAE grid calculations (non-self-consistent)...")
            
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
                
                # Prepare VASP parameters (without magmom for now)
                params_with_saxis = base_params_dict.copy()
                params_with_saxis['saxis'] = list(direction.saxis)
                
                # Remove magmom if present (we'll add it manually to INCAR)
                params_with_saxis.pop('magmom', None)
                
                # All grid calculations are non-self-consistent, reading from z folder
                params_with_saxis['icharg'] = 11  # Read CHGCAR for charge density
                params_with_saxis['istart'] = 1   # Read WAVECAR for wavefunctions
                params_with_saxis['lcharg'] = False  # Don't write CHGCAR
                params_with_saxis['lwave'] = False   # Don't write WAVECAR
                
                # Write VASP input files using ASE (without NCL magmom)
                from ase.calculators.vasp import Vasp
                calc = Vasp(directory=str(calc_dir), **params_with_saxis)
                calc.write_input(atoms_calc)
                
                # Create symlinks to WAVECAR and CHGCAR from z folder
                # These will be created after z calculation completes
                # For now, just note in a README
                readme_path = calc_dir / 'README.txt'
                with open(readme_path, 'w') as f:
                    f.write(f"This calculation reads WAVECAR and CHGCAR from ../z/\n")
                    f.write(f"Run the reference calculation in ../z/ first!\n")
                    f.write(f"\nAfter z/ completes, create symlinks:\n")
                    f.write(f"  ln -s ../z/WAVECAR ./WAVECAR\n")
                    f.write(f"  ln -s ../z/CHGCAR ./CHGCAR\n")
                
                # Now manually add non-collinear MAGMOM to INCAR
                incar_path = calc_dir / 'INCAR'
                with open(incar_path, 'r') as f:
                    incar_lines = f.readlines()
                
                # Find and replace/add MAGMOM line
                magmom_added = False
                for i, line in enumerate(incar_lines):
                    if line.strip().startswith('MAGMOM'):
                        # Replace existing MAGMOM with NCL format
                        magmom_str = '  '.join([f'{mx} {my} {mz}' for mx, my, mz in ncl_magmoms])
                        incar_lines[i] = f'MAGMOM = {magmom_str}\n'
                        magmom_added = True
                        break
                
                if not magmom_added:
                    # Add MAGMOM if not present
                    magmom_str = '  '.join([f'{mx} {my} {mz}' for mx, my, mz in ncl_magmoms])
                    incar_lines.append(f'MAGMOM = {magmom_str}\n')
                
                # Write back the modified INCAR
                with open(incar_path, 'w') as f:
                    f.writelines(incar_lines)
            
            # Generate helper submission script
            generate_submission_script(
                mae_dir, calc_dirs, directions, parallel_over_configurations,
                jobheader, calculator_command, environment_activate, environment_deactivate
            )
            
            print(f"✓ Created reference calculation in z/ folder")
            print(f"✓ Created {len(calc_dirs)} grid calculation directories")
            print(f"✓ Generated submission script with dependencies")
    
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
