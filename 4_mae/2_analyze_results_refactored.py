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
    
    # Load results using dedicated loader (SRP)
    base_path = Path.cwd()
    loader = MAEResultsLoader(base_path)
    
    # Load reference energy
    ref_energy = loader.load_reference_energy('z')
    if ref_energy is None:
        print("WARNING: Could not load reference energy from z/ directory")
        ref_energy = 0.0
    
    print(f"Reference energy: {ref_energy:.6f} eV")
    
    # Load grid results
    results = loader.load_theta_phi_results(directions)
    
    print(f"Loaded {len(results)} results out of {len(directions)} calculations")
    
    if len(results) == 0:
        print("ERROR: No results found!")
        print("Make sure calculations have completed and OSZICAR files exist.")
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
        plotter.plot_theta_phi_surface(
            theta_grid, phi_grid, energy_grid,
            output_file=f'mae_surface_{configuration}.png'
        )
    else:
        print(f"WARNING: Incomplete grid ({len(valid_directions)}/{expected_points})")
        print("Skipping surface plot. Complete all calculations first.")
    
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
    print(f"\nNext step: Use the easy/hard axes above for MAE curve calculation")
    print(f"Add to input.py and run 3_submit_mae_curve_refactored.py")


if __name__ == '__main__':
    main()
