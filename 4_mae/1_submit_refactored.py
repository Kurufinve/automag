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

# Shared helpers (single source of truth for configuration lookup and cell
# processing, also used by 1_submit.py)
from mae_utils import (
    load_collinear_configuration,
    process_structure as _process_structure_shared,
)

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

# Structure processing options
standardize_cell = True  # Apply standardization using SpacegroupAnalyzer
use_primitive_cell = True  # Convert to primitive cell (only if standardize_cell=True)

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
    """Load a magnetic configuration from the collinear results.

    Thin wrapper around the shared ``mae_utils.load_collinear_configuration``
    so that every MAE script resolves configurations the same way.
    """
    return load_collinear_configuration(
        path_to_automag=path_to_automag,
        formula=formula,
        configuration=configuration_name,
        struct_suffix=struct_suffix,
        calculator=calculator_type,
    )


def process_structure(structure_path: str,
                     magmoms: list,
                     output_name: str,
                     symmetrize: bool = True,
                     to_primitive: bool = True,
                     symprec: float = 0.1) -> tuple:
    """Standardize (optionally primitive) and emit NCL magmoms.

    Thin wrapper around the shared ``mae_utils.process_structure``.
    """
    return _process_structure_shared(
        structure_path=structure_path,
        magmoms=magmoms,
        output_name=output_name,
        symmetrize=symmetrize,
        to_primitive=to_primitive,
        symprec=symprec,
    )


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
    
    # Load input structure to get original formula
    input_structure = Structure.from_file(path_to_poscar)
    original_formula = input_structure.formula.replace(' ', '')
    
    print(f"\n{'=' * 70}")
    print(f"MAE CALCULATION FOR {original_formula} - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    
    # Load configuration from collinear calculations
    structure_path, final_magmoms, setting = load_configuration_from_collinear(
        path_to_automag, original_formula, configuration, struct_suffix, calculator
    )
    
    # Process structure with standardization options
    processed_file = f'setting{setting:03d}_{configuration}_processed.vasp'
    processed_path, ncl_magmoms, processed_structure = process_structure(
        structure_path, 
        final_magmoms, 
        processed_file,
        symmetrize=standardize_cell,
        to_primitive=use_primitive_cell
    )
    
    # Get formula from processed structure (may differ from original if primitive)
    processed_formula = processed_structure.formula.replace(' ', '')
    
    print(f"\nUsing formula from processed structure: {processed_formula}")
    if processed_formula != original_formula:
        print(f"  (Original formula was: {original_formula})")
    
    # Load processed structure for ASE
    atoms = read(processed_path)
    
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
            # Determine cell type for folder naming
            if not standardize_cell:
                cell_type_folder = 'input_cell'
            elif use_primitive_cell:
                cell_type_folder = 'primitive_cell'
            else:
                cell_type_folder = 'conventional_cell'
            
    # Get atom count from processed structure
            n_atoms = len(processed_structure)
            
            # Get reduced formula for new directory pattern
            reduced_formula = processed_structure.composition.reduced_formula.replace(' ', '')
            
            # Create MAE calculation directory using new pattern
            calcfold_path = Path(path_to_automag) / 'CalcFold'
            mae_base_dir = calcfold_path / f"{reduced_formula}{struct_suffix}" / f"{processed_formula}_{cell_type_folder}" / calculator / configuration
            mae_dir = mae_base_dir / f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type_folder}_{n_atoms}atoms"
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
            shutil.copy(processed_path, z_dir / 'POSCAR')
            
            # Create ASE atoms for reference calculation
            from ase.io import read as ase_read
            atoms_ref = ase_read(processed_path)
            
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
                shutil.copy(processed_path, calc_dir / 'POSCAR')
                
                # Create ASE atoms and write VASP inputs
                from ase.io import read as ase_read
                atoms_calc = ase_read(processed_path)
                
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
    
    # Save configuration info with dynamic filename
    # Determine cell type for config filename
    if not standardize_cell:
        cell_type_folder = 'input_cell'
    elif use_primitive_cell:
        cell_type_folder = 'primitive_cell'
    else:
        cell_type_folder = 'conventional_cell'
    
    # Get first kpts and encut values for config filename (configs are per-parameter set)
    kpts_val_config = kpts_list[0]
    encut_val_config = encut_list[0]
    n_atoms_config = len(processed_structure)
    
    # Determine cell type for filename generation
    if not standardize_cell:
        cell_type = 'input_cell'
    elif use_primitive_cell:
        cell_type = 'primitive_cell'
    else:
        cell_type = 'conventional_cell'

    filename_suffix = f"{configuration}_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms"

    config_filename = f'mae_config_{filename_suffix}.txt'
    
    with open(config_filename, 'w') as f:
        f.write(f"Configuration: {configuration}\n")
        f.write(f"Original formula: {original_formula}\n")
        f.write(f"Processed formula: {processed_formula}\n")
        # f.write(f"Structure file: {processed_file}\n")
        f.write(f"Structure file: {processed_path}\n")
        f.write(f"Standardization: {standardize_cell}\n")
        f.write(f"Primitive cell: {use_primitive_cell if standardize_cell else 'N/A'}\n")
        

        f.write(f"Cell type: {cell_type}\n")
        
        f.write(f"Grid size: {Nth}×{Nph}\n")
        f.write(f"NCL magmoms: {ncl_magmoms}\n")
        f.write(f"Kpts values: {kpts_list}\n")
        f.write(f"Encut values: {encut_list}\n")
        f.write(f"U = {U:.1f}, J = {J:.1f}\n")
    
    print(f"\n✓ Configuration saved to: {config_filename}")
    print(f"\n{'=' * 70}")
    print("To submit calculations, run the generated script:")
    print(f"  bash <mae_directory>/submit_mae_grid.sh")
    print(f"{'=' * 70}\n")


if __name__ == '__main__':
    main()
