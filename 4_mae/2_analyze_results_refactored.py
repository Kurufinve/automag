"""
Refactored MAE results analysis script.

Analyzes theta-phi grid results to find easy/hard axes and calculate MAE.

Follows SOLID principles with dedicated analyzer and plotter classes.
"""

import os
import sys
import numpy as np
from pathlib import Path

from pymatgen.core.structure import Structure

# Import refactored components
from core.services.mae_service import (
    MAEDirectionGenerator,
    MAEAnalyzer,
    MAEResultsLoader,
    MAEPlotter
)

# Default values
calculator = 'vasp'
struct_suffix = ''
Nph = 20
Nth = 10

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
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    
    # Load structure
    structure = Structure.from_file(path_to_poscar)
    formula = structure.formula.replace(' ', '')
    volume = structure.volume
    
    print(f"\n{'=' * 70}")
    print(f"MAE ANALYSIS FOR {formula} - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    
    # Generate directions that were calculated
    directions = MAEDirectionGenerator.generate_theta_phi_grid(Nth, Nph)
    print(f"Analyzing {len(directions)} calculations...")
    
    # Determine the base path for MAE calculations
    # The script can be run from:
    # 1. The MAE calculation directory itself (mae_U*_J*_K*_EN*/)
    # 2. The 4_mae directory (need to find the MAE directory)
    
    base_path = Path.cwd()
    
    # Check if we're in the MAE calculation directory (should have 'z' subfolder)
    if not (base_path / 'z').exists():
        # We're probably in 4_mae directory, need to find the MAE directory
        # Try to construct the path from input parameters
        try:
            # Get U and J values from params if available
            ldauu_val = params.get('ldauu', [0.0])
            ldauj_val = params.get('ldauj', [0.0])
            ldaul_val = params.get('ldaul', [])
            U = ldauu_val[next((i for i, x in enumerate(ldaul_val) if x > 0), 0)] if ldaul_val else 0.0
            J = ldauj_val[next((i for i, x in enumerate(ldaul_val) if x > 0), 0)] if ldaul_val else 0.0
            
            kpts_val = params['kpts'] if not isinstance(params['kpts'], list) else params['kpts'][0]
            encut_val = params['encut'] if not isinstance(params['encut'], list) else params['encut'][0]
            
            # Construct expected path
            calcfold_path = Path(path_to_automag) / 'CalcFold'
            mae_base_dir = calcfold_path / f"{formula}{struct_suffix}" / calculator / configuration
            mae_dir = mae_base_dir / f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}"
            
            if mae_dir.exists() and (mae_dir / 'z').exists():
                base_path = mae_dir
                print(f"Found MAE directory: {base_path}")
            else:
                print(f"ERROR: Could not find MAE calculation directory!")
                print(f"Expected path: {mae_dir}")
                print(f"\nPlease run this script from the MAE calculation directory:")
                print(f"  cd {mae_dir}")
                print(f"  python {Path(__file__).parent}/2_analyze_results_refactored.py")
                return
        except Exception as e:
            print(f"ERROR: Could not determine MAE directory: {e}")
            print(f"\nPlease run this script from within the MAE calculation directory.")
            print(f"The directory should contain 'z/' subfolder and 'PhTh_*' subfolders.")
            return
    else:
        print(f"Running from MAE directory: {base_path}")
    
    # Load results using dedicated loader (SRP)
    loader = MAEResultsLoader(base_path)
    
    # Load reference energy
    ref_energy = loader.load_reference_energy('z')
    if ref_energy is None:
        print("WARNING: Could not load reference energy from z/ directory")
        print("Check that the reference calculation completed successfully.")
        ref_energy = 0.0
    
    print(f"Reference energy: {ref_energy:.6f} eV")
    
    # Load grid results
    print("\nLoading grid calculation results...")
    results = loader.load_theta_phi_results(directions)
    
    print(f"Loaded {len(results)} results out of {len(directions)} calculations")
    
    if len(results) == 0:
        print("\nERROR: No results found!")
        print("\nPossible reasons:")
        print("  1. Calculations haven't completed yet")
        print("  2. OSZICAR files don't exist in PhTh_* directories")
        print("  3. Running from wrong directory")
        print(f"\nExpected structure:")
        print(f"  {base_path}/z/OSZICAR (reference)")
        print(f"  {base_path}/PhTh_0.01_0.01/OSZICAR")
        print(f"  {base_path}/PhTh_0.01_18.0/OSZICAR")
        print(f"  ... etc.")
        print(f"\nChecking a few directories:")
        for i, direction in enumerate(directions[:3]):  # Check first 3
            calc_dir = base_path / direction.name
            oszicar_path = calc_dir / 'OSZICAR'
            print(f"  {direction.name}/OSZICAR exists: {oszicar_path.exists()}")
        return
    
    # Extract energies in order
    energies = []
    valid_directions = []
    
    for direction in directions:
        if direction.name in results:
            energies.append(results[direction.name])
            valid_directions.append(direction)
    
    energies = np.array(energies)
    
    # Analyze using dedicated analyzer (SRP)
    analyzer = MAEAnalyzer()
    
    # Find easy and hard axes
    easy_dir, hard_dir = analyzer.find_easy_hard_axes(valid_directions, energies)
    
    # Calculate MAE
    mae_ev = analyzer.calculate_mae(np.min(energies), np.max(energies))
    mae_mj_m3 = analyzer.calculate_mae_per_volume(mae_ev, volume)
    
    # Display results
    print(f"\n{'=' * 70}")
    print("RESULTS")
    print(f"{'=' * 70}")
    print(f"Easy axis direction:")
    print(f"  Theta: {easy_dir.theta * 180 / np.pi:.2f}°")
    print(f"  Phi:   {easy_dir.phi * 180 / np.pi:.2f}°")
    print(f"  SAXIS: ({easy_dir.saxis[0]:.3f}, {easy_dir.saxis[1]:.3f}, {easy_dir.saxis[2]:.3f})")
    print(f"  Energy: {np.min(energies):.6f} eV")
    print()
    print(f"Hard axis direction:")
    print(f"  Theta: {hard_dir.theta * 180 / np.pi:.2f}°")
    print(f"  Phi:   {hard_dir.phi * 180 / np.pi:.2f}°")
    print(f"  SAXIS: ({hard_dir.saxis[0]:.3f}, {hard_dir.saxis[1]:.3f}, {hard_dir.saxis[2]:.3f})")
    print(f"  Energy: {np.max(energies):.6f} eV")
    print()
    print(f"MAE:")
    print(f"  {mae_ev:.6f} eV")
    print(f"  {mae_ev * 1000:.3f} meV")
    print(f"  {mae_mj_m3:.3f} MJ/m³")
    print(f"{'=' * 70}\n")
    
    # Reshape for plotting if we have complete grid
    expected_points = (Nth + 1) * (Nph + 1)
    if len(valid_directions) == expected_points:
        # Create grid for plotting
        theta_grid = np.zeros((Nph + 1, Nth + 1))
        phi_grid = np.zeros((Nph + 1, Nth + 1))
        energy_grid = np.zeros((Nph + 1, Nth + 1))
        
        idx = 0
        for i in range(Nph + 1):
            for j in range(Nth + 1):
                if idx < len(valid_directions):
                    theta_grid[i, j] = valid_directions[idx].theta
                    phi_grid[i, j] = valid_directions[idx].phi
                    energy_grid[i, j] = energies[idx] - np.min(energies)
                    idx += 1
        
        # Plot using dedicated plotter (SRP)
        plotter = MAEPlotter()
        
        # 1. 3D surface plot
        plotter.plot_theta_phi_surface(
            theta_grid, phi_grid, energy_grid,
            output_file=f'mae_surface_{configuration}.png'
        )
        
        # 2. 2D E vs Theta plot at different Phi values (in MJ/m³)
        plotter.plot_energy_vs_theta_at_phi(
            theta_grid, phi_grid, energy_grid,
            volume=volume,
            output_file=f'mae_theta_phi_{configuration}.png'
        )
    else:
        print(f"WARNING: Incomplete grid ({len(valid_directions)}/{expected_points})")
        print("Skipping plots. Complete all calculations first.")
    
    # Save results to file
    output_file = f'mae_results_{configuration}.txt'
    with open(output_file, 'w') as f:
        f.write(f"MAE Analysis Results for {formula} - {configuration}\n")
        f.write(f"{'=' * 70}\n\n")
        f.write(f"Easy axis:\n")
        f.write(f"  Theta: {easy_dir.theta * 180 / np.pi:.2f}°\n")
        f.write(f"  Phi:   {easy_dir.phi * 180 / np.pi:.2f}°\n")
        f.write(f"  SAXIS: ({easy_dir.saxis[0]:.3f}, {easy_dir.saxis[1]:.3f}, {easy_dir.saxis[2]:.3f})\n")
        f.write(f"  Energy: {np.min(energies):.6f} eV\n\n")
        f.write(f"Hard axis:\n")
        f.write(f"  Theta: {hard_dir.theta * 180 / np.pi:.2f}°\n")
        f.write(f"  Phi:   {hard_dir.phi * 180 / np.pi:.2f}°\n")
        f.write(f"  SAXIS: ({hard_dir.saxis[0]:.3f}, {hard_dir.saxis[1]:.3f}, {hard_dir.saxis[2]:.3f})\n")
        f.write(f"  Energy: {np.max(energies):.6f} eV\n\n")
        f.write(f"MAE:\n")
        f.write(f"  {mae_ev:.6f} eV\n")
        f.write(f"  {mae_ev * 1000:.3f} meV\n")
        f.write(f"  {mae_mj_m3:.3f} MJ/m³\n\n")
        f.write(f"For MAE curve calculation, use:\n")
        f.write(f"  easy_axis = np.array([{easy_dir.saxis[0]:.3f}, {easy_dir.saxis[1]:.3f}, {easy_dir.saxis[2]:.3f}])\n")
        f.write(f"  hard_axis = np.array([{hard_dir.saxis[0]:.3f}, {hard_dir.saxis[1]:.3f}, {hard_dir.saxis[2]:.3f}])\n")
    
    print(f"Results saved to: {output_file}")
    
    if len(valid_directions) == expected_points:
        print(f"\nPlots generated:")
        print(f"  1. 3D surface plot: mae_surface_{configuration}.png")
        print(f"  2. 2D E(θ) at different φ: mae_theta_phi_{configuration}.png")
    
    print(f"\nNext step: Use the easy/hard axes above for MAE curve calculation")
    print(f"Add to input.py and run 3_submit_mae_curve_refactored.py")


if __name__ == '__main__':
    main()
