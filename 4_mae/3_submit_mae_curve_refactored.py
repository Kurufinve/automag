"""
Refactored MAE curve submission script.

Submits calculations along rotation path from easy to hard axis.
Uses the same workflow approach as 1_submit_refactored.py with:
- Proper directory structure and naming
- Reference + grid calculation workflow
- VASP input file generation
- Parallel/sequential submission scripts
- Symlink management for WAVECAR/CHGCAR

Run this after 2_analyze_results_refactored.py to get easy/hard axes.
"""

import os
import sys
import numpy as np
from pathlib import Path
from ase.io import read

from pymatgen.core.structure import Structure

# Import refactored components
from core.services.mae_service import MAEDirectionGenerator
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


def generate_mae_curve_submission_script(mae_curve_dir: Path,
                                          calc_dirs: list,
                                          direction_names: list,
                                          parallel: bool,
                                          jobheader: str,
                                          calculator_command: str,
                                          environment_activate: str,
                                          environment_deactivate: str):
    """
    Generate SLURM submission helper script for MAE curve calculations.
    
    Args:
        mae_curve_dir: Base MAE curve directory
        calc_dirs: List of calculation directories
        direction_names: List of direction names
        parallel: If True, create parallel submission script; if False, sequential
        jobheader: SLURM job header template
        calculator_command: Command to run VASP
        environment_activate: Command to activate environment
        environment_deactivate: Command to deactivate environment
    """
    
    if parallel:
        # Create parallel submission script (one click to submit all jobs)
        script_path = mae_curve_dir / 'submit_mae_curve.sh'
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# MAE Curve Parallel Submission Script\n")
            f.write("# This script submits all MAE curve calculations in parallel\n\n")
            
            f.write(f"echo 'Submitting {len(calc_dirs)} MAE curve calculations...'\n\n")
            
            for calc_dir, dir_name in zip(calc_dirs, direction_names):
                job_name = f"mae_curve_{dir_name}"
                
                f.write(f"# Submit {dir_name}\n")
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
                f.write(f"cd {mae_curve_dir}\n\n")
            
            f.write("echo 'All jobs submitted!'\n")
            f.write("echo 'Monitor with: squeue -u $USER'\n")
        
        script_path.chmod(0o755)
        print(f"Created parallel submission script: {script_path}")
    
    else:
        # Create sequential submission script
        script_path = mae_curve_dir / 'submit_mae_curve_sequential.sh'
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# MAE Curve Sequential Submission Script\n")
            f.write("# This script runs MAE curve calculations sequentially in a single SLURM job\n\n")
            
            f.write(jobheader + "\n")
            f.write("#SBATCH -J mae_curve_sequential\n")
            f.write(f"#SBATCH -o {mae_curve_dir}/mae_curve_sequential-%j.out\n")
            f.write(f"#SBATCH -e {mae_curve_dir}/mae_curve_sequential-%j.err\n\n")
            
            f.write(f"{environment_activate}\n\n")
            
            # Run calculations sequentially
            f.write(f"echo 'Running {len(calc_dirs)} MAE curve calculations sequentially...'\n\n")
            
            for calc_dir, dir_name in zip(calc_dirs, direction_names):
                f.write(f"echo 'Running calculation for {dir_name}...'\n")
                f.write(f"cd {calc_dir}\n")
                f.write(f"{calculator_command}\n")
                f.write("\n")
                f.write(f"# Check if calculation completed successfully\n")
                f.write(f"if [ $? -ne 0 ]; then\n")
                f.write(f"    echo 'ERROR: Calculation failed for {dir_name}'\n")
                f.write(f"    exit 1\n")
                f.write(f"fi\n")
                f.write(f"echo 'Completed {dir_name}'\n\n")
            
            f.write(f"{environment_deactivate}\n")
            f.write("echo 'All MAE curve calculations completed!'\n")
        
        script_path.chmod(0o755)
        print(f"Created sequential submission script: {script_path}")


def main():
    """Main execution function."""
    from ase.calculators.vasp import Vasp
    import glob
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    if not path_to_automag:
        print("ERROR: AUTOMAG_PATH environment variable not set!")
        return
    
    calcfold_path = Path(path_to_automag) / 'CalcFold'
    calcfold_path.mkdir(parents=True, exist_ok=True)
    
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
    structure_path = None
    try:
        # Try to find the processed structure file (new naming)
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
    pmg_structure = Structure.from_file(structure_path)
    print(f"Structure loaded from: {structure_path}")
    print(f"  Atoms: {len(atoms)}")
    print(f"  Formula: {pmg_structure.composition.reduced_formula}")
    
    # Load calculation parameters from config file
    ncl_magmoms = None
    cell_type = 'unknown_cell'
    ldauu_val = [0.0]
    ldauj_val = [0.0]
    kpts_val = None
    encut_val = None
    
    config_file = f'{configuration}_mae_config.txt'
    
    if os.path.exists(config_file):
        print(f"\nReading configuration from: {config_file}")
        try:
            with open(config_file, 'r') as f:
                for line in f:
                    if 'NCL magmoms:' in line:
                        # Parse the magmoms - format is a list of tuples
                        magmoms_str = line.split('NCL magmoms:')[1].strip()
                        # Use eval to parse the list of tuples
                        # Format: [(0, 0, m1), (0, 0, m2), ...]
                        ncl_magmoms = eval(magmoms_str)
                        print(f"  → Loaded {len(ncl_magmoms)} NCL magnetic moments")
                    elif 'Cell type:' in line:
                        cell_type = line.split('Cell type:')[1].strip()
                        print(f"  → Cell type: {cell_type}")
                    elif 'LDAUU:' in line:
                        ldauu_str = line.split('LDAUU:')[1].strip()
                        ldauu_val = eval(ldauu_str)
                    elif 'LDAUJ:' in line:
                        ldauj_str = line.split('LDAUJ:')[1].strip()
                        ldauj_val = eval(ldauj_str)
                    elif 'K-points:' in line:
                        kpts_val = int(line.split('K-points:')[1].strip())
                    elif 'ENCUT:' in line:
                        encut_val = int(line.split('ENCUT:')[1].strip())
        except Exception as e:
            print(f"Warning: Could not parse some parameters from config: {e}")
            print(f"Will try to use values from input.py if available")
    else:
        print(f"Warning: Config file {config_file} not found")
        print(f"Will try to use values from input.py if available")
    
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
    
    # Get U and J values
    ldaul_val = params.get('ldaul', [])
    if kpts_val is None:
        kpts_val = params['kpts'] if not isinstance(params['kpts'], list) else params['kpts'][0]
    if encut_val is None:
        encut_val = params['encut'] if not isinstance(params['encut'], list) else params['encut'][0]
    
    # Extract U and J for magnetic atoms
    U = ldauu_val[next((i for i, x in enumerate(ldaul_val) if x > 0), 0)] if ldaul_val else 0.0
    J = ldauj_val[next((i for i, x in enumerate(ldaul_val) if x > 0), 0)] if ldaul_val else 0.0
    
    print(f"\nCalculation parameters:")
    print(f"  U = {U}, J = {J}")
    print(f"  K-points = {kpts_val}")
    print(f"  ENCUT = {encut_val} eV")
    print(f"  Cell type = {cell_type}")
    
    # Generate MAE curve directions
    print(f"\nGenerating MAE curve directions...")
    directions = MAEDirectionGenerator.generate_mae_curve(
        n_points=N_MAE,
        easy_axis=easy_axis,
        hard_axis=hard_axis
    )
    print(f"  → Generated {len(directions)} directions")
    
    # Create directory structure
    formula = pmg_structure.composition.reduced_formula
    mae_curve_dirname = f'mae_curve_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}'
    
    # Build full path: CalcFold/{formula}{struct_suffix}/{calculator}/{configuration}/{mae_curve_dir}
    base_calc_dir = calcfold_path / f"{formula}{struct_suffix}" / calculator / configuration
    mae_curve_dir = base_calc_dir / mae_curve_dirname
    
    print(f"\nCreating calculation directories...")
    print(f"  Base path: {mae_curve_dir}")
    
    mae_curve_dir.mkdir(parents=True, exist_ok=True)
    
    # Create calculation directories for each point
    calc_dirs = []
    direction_names = []
    
    for direction in directions:
        calc_dir = mae_curve_dir / direction.name
        calc_dir.mkdir(parents=True, exist_ok=True)
        calc_dirs.append(calc_dir)
        direction_names.append(direction.name)
    
    print(f"  → Created {len(calc_dirs)} calculation directories")
    
    # Prepare base VASP parameters
    vasp_params = params.copy()
    vasp_params.update({
        'voskown': 1,
        'lnoncollinear': True,
        'lsorbit': True,
        'gga_compat': False,
        'icharg': 11,  # Read charge density from CHGCAR
        'istart': 1,   # Read wavefunction from WAVECAR
        'lcharg': False,  # Don't write CHGCAR for curve points
        'lwave': False,   # Don't write WAVECAR for curve points
    })
    
    # Write VASP input files for each direction
    print(f"\nWriting VASP input files...")
    
    for i, (direction, calc_dir) in enumerate(zip(directions, calc_dirs)):
        # Set magnetic moments with rotated direction
        rotated_magmoms = []
        for mx, my, mz in ncl_magmoms:
            # Rotate magnetic moment to point along direction
            m_magnitude = np.sqrt(mx**2 + my**2 + mz**2)
            # Convert saxis tuple to numpy array for element-wise multiplication
            saxis_array = np.array(direction.saxis)
            rotated_moment = saxis_array * m_magnitude
            rotated_magmoms.append(tuple(rotated_moment))
        
        # Update VASP parameters for this direction
        direction_params = vasp_params.copy()
        direction_params['saxis'] = list(direction.saxis)
        
        # Create VASP calculator
        calc = Vasp(
            directory=str(calc_dir),
            **direction_params
        )
        
        # Set atoms and magmoms
        atoms_copy = atoms.copy()
        atoms_copy.set_initial_magnetic_moments(
            [m[2] for m in ncl_magmoms]  # Use z-component as collinear magmom
        )
        atoms_copy.calc = calc
        
        # Write input files
        calc.initialize(atoms_copy)
        calc.write_input(atoms_copy)
        
        # Manually set NCL MAGMOM in INCAR
        incar_path = calc_dir / 'INCAR'
        with open(incar_path, 'r') as f:
            incar_lines = f.readlines()
        
        # Replace MAGMOM line with NCL format
        with open(incar_path, 'w') as f:
            for line in incar_lines:
                if line.strip().startswith('MAGMOM'):
                    # Write NCL MAGMOM format
                    magmom_str = ' '.join([f"{mx} {my} {mz}" for mx, my, mz in rotated_magmoms])
                    f.write(f"MAGMOM = {magmom_str}\n")
                else:
                    f.write(line)
        
        if (i + 1) % 5 == 0 or (i + 1) == len(directions):
            print(f"  → Written {i + 1}/{len(directions)} input files")
    
    print(f"  ✓ All input files created")
    
    # Create symlinks to reference WAVECAR and CHGCAR from MAE grid calculation
    # Look for the reference 'z' directory from the MAE grid calculation
    mae_grid_dir = base_calc_dir.parent / f'mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}'
    z_ref_dir = mae_grid_dir / 'z'
    
    print(f"\nCreating symlinks to reference files...")
    if z_ref_dir.exists():
        wavecar_src = z_ref_dir / 'WAVECAR'
        chgcar_src = z_ref_dir / 'CHGCAR'
        
        if wavecar_src.exists() and chgcar_src.exists():
            symlink_count = 0
            for calc_dir in calc_dirs:
                wavecar_dst = calc_dir / 'WAVECAR'
                chgcar_dst = calc_dir / 'CHGCAR'
                
                # Remove existing symlinks/files if they exist
                if wavecar_dst.exists() or wavecar_dst.is_symlink():
                    wavecar_dst.unlink()
                if chgcar_dst.exists() or chgcar_dst.is_symlink():
                    chgcar_dst.unlink()
                
                # Create symlinks
                try:
                    wavecar_dst.symlink_to(wavecar_src)
                    chgcar_dst.symlink_to(chgcar_src)
                    symlink_count += 1
                except Exception as e:
                    print(f"  Warning: Could not create symlink in {calc_dir.name}: {e}")
            
            print(f"  → Created symlinks in {symlink_count} directories")
            print(f"  → Reference: {z_ref_dir}")
        else:
            print(f"  Warning: WAVECAR/CHGCAR not found in {z_ref_dir}")
            print(f"  → Run MAE grid calculation first (1_submit_refactored.py)")
    else:
        print(f"  Warning: Reference directory not found: {z_ref_dir}")
        print(f"  → Run MAE grid calculation first (1_submit_refactored.py)")
    
    # Generate submission scripts
    print(f"\nGenerating submission scripts...")
    
    # Generate both parallel and sequential scripts
    generate_mae_curve_submission_script(
        mae_curve_dir=mae_curve_dir,
        calc_dirs=calc_dirs,
        direction_names=direction_names,
        parallel=True,
        jobheader=jobheader,
        calculator_command=calculator_command,
        environment_activate=environment_activate,
        environment_deactivate=environment_deactivate
    )
    
    generate_mae_curve_submission_script(
        mae_curve_dir=mae_curve_dir,
        calc_dirs=calc_dirs,
        direction_names=direction_names,
        parallel=False,
        jobheader=jobheader,
        calculator_command=calculator_command,
        environment_activate=environment_activate,
        environment_deactivate=environment_deactivate
    )
    
    # Save configuration file
    config_output_path = mae_curve_dir / f'{configuration}_mae_curve_config.txt'
    print(f"\nSaving configuration...")
    with open(config_output_path, 'w') as f:
        f.write(f"MAE Curve Configuration\n")
        f.write(f"={'=' * 50}\n")
        f.write(f"Configuration: {configuration}\n")
        f.write(f"Formula: {formula}\n")
        f.write(f"Number of atoms: {len(atoms)}\n")
        f.write(f"Easy axis: {easy_axis}\n")
        f.write(f"Hard axis: {hard_axis}\n")
        f.write(f"Number of curve points: {len(directions)}\n")
        f.write(f"NCL magmoms: {ncl_magmoms}\n")
        f.write(f"LDAUU: {ldauu_val}\n")
        f.write(f"LDAUJ: {ldauj_val}\n")
        f.write(f"K-points: {kpts_val}\n")
        f.write(f"ENCUT: {encut_val}\n")
        f.write(f"Cell type: {cell_type}\n")
        f.write(f"Structure file: {structure_path}\n")
        f.write(f"Reference directory: {z_ref_dir if z_ref_dir.exists() else 'Not found'}\n")
    
    print(f"  → Saved to: {config_output_path}")
    
    # Print summary
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE SETUP COMPLETE")
    print(f"{'=' * 70}")
    print(f"\nCalculation directory: {mae_curve_dir}")
    print(f"Number of calculations: {len(calc_dirs)}")
    print(f"\nTo submit calculations:")
    print(f"  cd {mae_curve_dir}")
    print(f"  ./submit_mae_curve.sh  # Parallel submission")
    print(f"  # OR")
    print(f"  sbatch submit_mae_curve_sequential.sh  # Sequential in one job")
    print(f"\nAfter calculations complete:")
    print(f"  cd {Path.cwd()}")
    print(f"  python 4_plot_mae_curve_refactored.py")
    print()


if __name__ == '__main__':
    main()
