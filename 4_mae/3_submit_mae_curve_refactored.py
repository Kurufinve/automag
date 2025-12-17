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
get_full_mae_curve = True  # If True, generate full 360° circle; if False, only easy-to-hard segment

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
    
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE CALCULATION - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    
    # Load processed structure file (matching approach from 2_analyze_results_refactored.py)
    # Determine expected cell type from input parameters to find the correct processed file
    import glob
    
    # Access standardization parameters from global namespace (matching other MAE scripts)
    standardize_cell_check = globals().get('standardize_cell', True)
    use_primitive_cell_check = globals().get('use_primitive_cell', True)
    
    # Determine expected cell type for file search
    if not standardize_cell_check:
        expected_cell_type = 'original'
    elif use_primitive_cell_check:
        expected_cell_type = 'primitive'
    else:
        expected_cell_type = 'conventional'
    
    # Search for processed file with the specific cell type
    # Pattern: setting*_{configuration}_{cell_type}.vasp
    processed_files = glob.glob(f'setting*_{configuration}_{expected_cell_type}.vasp')
    
    if processed_files:
        structure_path = processed_files[0]
        print(f"Found processed structure: {structure_path}")
        print(f"  Cell type: {expected_cell_type}")
    else:
        # Fallback: try legacy naming patterns
        print(f"Warning: No processed structure file found matching pattern: setting*_{configuration}_{expected_cell_type}.vasp")
        print(f"Trying legacy naming patterns...")
        
        # Try old "processed" naming (without cell type suffix)
        legacy_processed = glob.glob(f'*{configuration}_processed.vasp')
        legacy_standardized = glob.glob(f'*{configuration}_standardized.vasp')
        
        if legacy_processed:
            structure_path = legacy_processed[0]
            print(f"Using legacy processed file: {structure_path}")
        elif legacy_standardized:
            structure_path = legacy_standardized[0]
            print(f"Using legacy standardized file: {structure_path}")
        else:
            print("ERROR: Structure file not found!")
            print(f"Expected pattern: setting*_{configuration}_{expected_cell_type}.vasp")
            print(f"  Or legacy: *{configuration}_processed.vasp")
            print(f"  Or legacy: *{configuration}_standardized.vasp")
            print("\nPlease run 1_submit_refactored.py first to generate the structure file.")
            return
    
    # Load structure using both ASE and pymatgen
    try:
        atoms = read(structure_path)
        pmg_structure = Structure.from_file(structure_path)
        print(f"Structure loaded successfully")
        print(f"  Atoms: {len(atoms)}")
        print(f"  Formula: {pmg_structure.formula.replace(' ', '')}")
        print(f"  Reduced formula: {pmg_structure.composition.reduced_formula}")
    except Exception as e:
        print(f"ERROR loading structure from {structure_path}: {e}")
        return
    
    # Load calculation parameters from config file
    ncl_magmoms = None
    cell_type = None  # Will be determined from config or input params
    processed_formula = None  # Will be read from config or determined from structure
    ldauu_val = [0.0]
    ldauj_val = [0.0]
    kpts_val = None
    encut_val = None
    easy_axis = None  # Will be read from config file
    hard_axis = None  # Will be read from config file
    
    # Get U and J values for config filename (matching 1_submit_refactored.py exactly)
    ldauu_val = params.get('ldauu', [0.0])
    ldauj_val = params.get('ldauj', [0.0])
    ldaul_val = params.get('ldaul', [])
    
    # Extract U and J for config filename
    U = ldauu_val[next(i for i, x in enumerate(ldaul_val) if x > 0)] if ldaul_val else 0.0
    J = ldauj_val[next(i for i, x in enumerate(ldaul_val) if x > 0)] if ldaul_val else 0.0
    
    # Get kpts and encut for config filename
    kpts_val = params['kpts'] if not isinstance(params['kpts'], list) else params['kpts'][0]
    encut_val = params['encut'] if not isinstance(params['encut'], list) else params['encut'][0]
    
    # Determine cell_type for config filename
    standardize_cell = globals().get('standardize_cell', True)
    use_primitive_cell = globals().get('use_primitive_cell', True)
    
    if not standardize_cell:
        cell_type = 'input_cell'
    elif use_primitive_cell:
        cell_type = 'primitive_cell'
    else:
        cell_type = 'conventional_cell'
    
    # Get atom count for config filename
    n_atoms = len(atoms)
    
    # Construct dynamic config filename
    config_file = f'{configuration}_mae_config_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms.txt'
    
    # Try to find config file if exact match not found
    if not os.path.exists(config_file):
        # Try to find any matching config file with wildcard pattern
        import glob
        pattern = f'{configuration}_mae_config_*.txt'
        matching_configs = glob.glob(pattern)
        if matching_configs:
            config_file = matching_configs[0]
            print(f"Using config file: {config_file}")
        else:
            print(f"Warning: No config file found matching pattern: {pattern}")
    
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
                    elif 'Processed formula:' in line:
                        processed_formula = line.split('Processed formula:')[1].strip()
                        print(f"  → Processed formula: {processed_formula}")
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
                    elif line.startswith('Easy axis:'):
                        # Parse easy axis: format is "Easy axis: [x, y, z]"
                        axis_str = line.split('Easy axis:')[1].strip()
                        # Remove brackets and parse
                        axis_str = axis_str.strip('[]')
                        easy_axis = np.array([float(x.strip()) for x in axis_str.split(',')])
                        print(f"  → Loaded easy axis: {easy_axis}")
                    elif line.startswith('Hard axis:'):
                        # Parse hard axis: format is "Hard axis: [x, y, z]"
                        axis_str = line.split('Hard axis:')[1].strip()
                        # Remove brackets and parse
                        axis_str = axis_str.strip('[]')
                        hard_axis = np.array([float(x.strip()) for x in axis_str.split(',')])
                        print(f"  → Loaded hard axis: {hard_axis}")
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
    
    # Check easy_axis and hard_axis - try config first, then input.py
    if easy_axis is None:
        if 'easy_axis' in globals():
            easy_axis = globals()['easy_axis']
            print(f"Using easy_axis from input.py: {easy_axis}")
        else:
            print("ERROR: easy_axis not found!")
            print("\neasy_axis should be either:")
            print("  1. In the config file (created by 2_analyze_results_refactored.py)")
            print("  2. Defined in input.py as: easy_axis = np.array([x, y, z])")
            print("\nPlease run 2_analyze_results_refactored.py first to determine easy/hard axes.")
            return
    
    if hard_axis is None:
        if 'hard_axis' in globals():
            hard_axis = globals()['hard_axis']
            print(f"Using hard_axis from input.py: {hard_axis}")
        else:
            print("ERROR: hard_axis not found!")
            print("\nhard_axis should be either:")
            print("  1. In the config file (created by 2_analyze_results_refactored.py)")
            print("  2. Defined in input.py as: hard_axis = np.array([x, y, z])")
            print("\nPlease run 2_analyze_results_refactored.py first to determine easy/hard axes.")
            return
    
    # Determine processed formula if not loaded from config
    # Use the actual formula with atom counts, not the reduced formula
    if processed_formula is None:
        # Get formula from structure file - this includes actual atom counts
        processed_formula = pmg_structure.formula.replace(' ', '')
        print(f"  → Processed formula determined from structure: {processed_formula}")
    
    print(f"\nCalculation parameters:")
    print(f"  Formula = {processed_formula}")
    print(f"  U = {U}, J = {J}")
    print(f"  K-points = {kpts_val}")
    print(f"  ENCUT = {encut_val} eV")
    print(f"  Cell type = {cell_type}")
    print(f"  Number of atoms = {n_atoms}")
    
    print(f"\nMAE axes:")
    print(f"  Easy axis: {easy_axis}")
    print(f"  Hard axis: {hard_axis}")
    
    # Generate MAE curve directions
    # Check if user wants full circle or just easy-to-hard segment
    full_curve = globals().get('get_full_mae_curve', True)
    
    print(f"\nGenerating MAE curve directions...")
    if full_curve:
        print(f"  → Mode: Full 360° circular curve")
        directions = MAEDirectionGenerator.generate_mae_curve_full(
            n_points=N_MAE,
            easy_axis=easy_axis,
            hard_axis=hard_axis
        )
    else:
        print(f"  → Mode: Segment from easy to hard axis")
        directions = MAEDirectionGenerator.generate_mae_curve(
            n_points=N_MAE,
            easy_axis=easy_axis,
            hard_axis=hard_axis
        )
    print(f"  → Generated {len(directions)} directions")
    print(f"  → Angle range: {directions[0].angle:.1f}° to {directions[-1].angle:.1f}°")
    
    # Create MAE curve calculation directory using EXACT same pattern as 1_submit_refactored.py
    # Pattern: CalcFold/{formula}{struct_suffix}/{calculator}/{configuration}/mae_curve_U{U}_J{J}_K{kpts}_EN{encut}_{cell_type}_{n_atoms}atoms
    calcfold_path = Path(path_to_automag) / 'CalcFold'
    mae_base_dir = calcfold_path / f"{processed_formula}{struct_suffix}" / calculator / configuration
    mae_curve_dir = mae_base_dir / f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms"
    mae_curve_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nCreating MAE curve directory: {mae_curve_dir}")
    
    # Create calculation directories for each point
    calc_dirs = []
    direction_names = []
    
    for direction in directions:
        # Name each folder as RtMAE_{angle} where angle is in degrees
        dir_name = f"RtMAE_{direction.angle:.1f}"
        calc_dir = mae_curve_dir / dir_name
        calc_dir.mkdir(parents=True, exist_ok=True)
        calc_dirs.append(calc_dir)
        direction_names.append(dir_name)
    
    print(f"  → Created {len(calc_dirs)} RtMAE calculation directories inside mae_curve_dir")
    
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
        # IMPORTANT: Do NOT rotate magnetic moments!
        # Only SAXIS (spin quantization axis) rotates
        # Magnetic moments remain constant - only their orientation relative to SAXIS changes
        
        # Update VASP parameters for this direction
        direction_params = vasp_params.copy()
        direction_params['saxis'] = list(direction.saxis)  # Rotate SAXIS only
        
        # Create VASP calculator
        calc = Vasp(
            directory=str(calc_dir),
            **direction_params
        )
        
        # Set atoms and magmoms (use original NCL magmoms without rotation)
        atoms_copy = atoms.copy()
        atoms_copy.set_initial_magnetic_moments(
            [m[2] for m in ncl_magmoms]  # Use z-component as collinear magmom
        )
        atoms_copy.calc = calc
        
        # Write input files
        calc.initialize(atoms_copy)
        calc.write_input(atoms_copy)
        
        # Manually set NCL MAGMOM in INCAR (use original magmoms, not rotated)
        incar_path = calc_dir / 'INCAR'
        with open(incar_path, 'r') as f:
            incar_lines = f.readlines()
        
        # Replace MAGMOM line with NCL format
        with open(incar_path, 'w') as f:
            for line in incar_lines:
                if line.strip().startswith('MAGMOM'):
                    # Write NCL MAGMOM format with ORIGINAL (unrotated) moments
                    magmom_str = ' '.join([f"{mx} {my} {mz}" for mx, my, mz in ncl_magmoms])
                    f.write(f"MAGMOM = {magmom_str}\n")
                else:
                    f.write(line)
        
        if (i + 1) % 5 == 0 or (i + 1) == len(directions):
            print(f"  → Written {i + 1}/{len(directions)} input files")
    
    print(f"  ✓ All input files created")
    
    # Create symlinks to reference WAVECAR and CHGCAR from MAE grid calculation's z/ folder
    # Look for the reference 'z' directory from the MAE grid calculation (sibling directory)
    # Pattern: CalcFold/{formula}{struct_suffix}/{calculator}/{configuration}/mae_U{U}_J{J}_K{kpts}_EN{encut}_{cell_type}_{n_atoms}atoms/z/
    mae_grid_dir = mae_base_dir / f'mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms'
    z_ref_dir = mae_grid_dir / 'z'
    
    print(f"\nCreating symlinks to reference files...")
    print(f"  Reference directory: {z_ref_dir}")
    
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
                
                # Create symlinks (relative paths for portability)
                try:
                    # Use relative path from calc_dir to z_ref_dir
                    rel_wavecar = os.path.relpath(wavecar_src, calc_dir)
                    rel_chgcar = os.path.relpath(chgcar_src, calc_dir)
                    
                    os.symlink(rel_wavecar, wavecar_dst)
                    os.symlink(rel_chgcar, chgcar_dst)
                    symlink_count += 1
                except Exception as e:
                    print(f"  Warning: Could not create symlink in {calc_dir.name}: {e}")
            
            print(f"  → Created symlinks in {symlink_count}/{len(calc_dirs)} directories")
            print(f"  ✓ All RtMAE directories ready for non-self-consistent calculations")
        else:
            print(f"  Warning: WAVECAR/CHGCAR not found in {z_ref_dir}")
            print(f"  → Run MAE grid calculation first (1_submit_refactored.py)")
            print(f"  → Make sure the z/ calculation has completed successfully")
    else:
        print(f"  Warning: Reference directory not found: {z_ref_dir}")
        print(f"  → Run MAE grid calculation first (1_submit_refactored.py)")
    
    # Generate submission scripts in the mae_curve_dir (matching 1_submit_refactored.py pattern)
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
    
    # Save configuration file in the mae_curve_dir (matching 1_submit_refactored.py pattern)
    config_output_path = mae_curve_dir / f'{configuration}_mae_curve_config.txt'
    print(f"\nSaving configuration...")
    with open(config_output_path, 'w') as f:
        f.write(f"MAE Curve Configuration\n")
        f.write(f"={'=' * 50}\n")
        f.write(f"Configuration: {configuration}\n")
        f.write(f"Formula: {processed_formula}\n")
        f.write(f"Number of atoms: {len(atoms)}\n")
        f.write(f"Easy axis: {easy_axis}\n")
        f.write(f"Hard axis: {hard_axis}\n")
        f.write(f"Full MAE curve: {full_curve}\n")
        f.write(f"Number of curve points: {len(directions)}\n")
        f.write(f"Angle range: {directions[0].angle:.1f} to {directions[-1].angle:.1f} degrees\n")
        f.write(f"NCL magmoms: {ncl_magmoms}\n")
        f.write(f"LDAUU: {ldauu_val}\n")
        f.write(f"LDAUJ: {ldauj_val}\n")
        f.write(f"K-points: {kpts_val}\n")
        f.write(f"ENCUT: {encut_val}\n")
        f.write(f"Cell type: {cell_type}\n")
        f.write(f"Structure file: {structure_path}\n")
        f.write(f"Reference directory: {z_ref_dir if z_ref_dir.exists() else 'Not found'}\n")
        f.write(f"Directory pattern: mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms/RtMAE_{{angle}}\n")
    
    print(f"  → Saved to: {config_output_path}")
    
    # Print summary
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE SETUP COMPLETE")
    print(f"{'=' * 70}")
    print(f"\nMAE curve directory: {mae_curve_dir}")
    print(f"Curve mode: {'Full 360° circle' if full_curve else 'Easy-to-hard segment'}")
    print(f"Number of RtMAE calculations: {len(calc_dirs)}")
    print(f"Angle range: {directions[0].angle:.1f}° to {directions[-1].angle:.1f}°")
    print(f"\nDirectory structure (matching 1_submit_refactored.py):")
    print(f"  {mae_base_dir}/")
    print(f"    ├── mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms/  (MAE grid reference)")
    print(f"    │   ├── z/  (reference self-consistent calculation)")
    print(f"    │   ├── PhTh_*/  (grid calculations)")
    print(f"    │   └── ...")
    print(f"    └── mae_curve_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}_{cell_type}_{n_atoms}atoms/")
    print(f"        ├── RtMAE_{directions[0].angle:.1f}/  (starting point)")
    print(f"        ├── RtMAE_{{angle}}/  (...)")
    print(f"        └── RtMAE_{directions[-1].angle:.1f}/  (ending point)")
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
