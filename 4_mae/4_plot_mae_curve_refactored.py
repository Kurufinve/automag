"""
Refactored MAE curve plotting and analysis script.

Plots MAE curve results and calculates magnetic properties.

Run this after 3_submit_mae_curve_refactored.py completes.

Follows:
- Single Responsibility Principle: Separate classes for different tasks
- Dependency Inversion Principle: Depends on abstractions
- Open/Closed Principle: Extensible for new analysis methods
"""

import os
import sys
import re
import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Oszicar

# Import refactored components
from core.services.mae_service import (
    MAEDirection,
    MAEAnalyzer,
    MAEResultsLoader,
    MAEPlotter
)

# Physical constants
Bh = 9.274009994e-24  # Bohr magneton in J/T
mu0 = 4 * np.pi * 1e-7  # Vacuum permeability in H/m
eV = 1.602176634e-19  # eV to J
Ang = 1e-10  # Angstrom to m

# Default values
use_fireworks = False
calculator = 'vasp'
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


def read_structure_properties(path_to_poscar: str) -> dict:
    """
    Read structure properties from POSCAR.
    
    Args:
        path_to_poscar: Path to POSCAR file
        
    Returns:
        Dictionary with volume, composition, etc.
    """
    struct = Structure.from_file(path_to_poscar)
    
    return {
        'volume': struct.volume,
        'composition': struct.composition,
        'num_atoms': struct.composition.num_atoms,
        'formula': struct.formula.replace(' ', '')
    }


def extract_magnetization(outcar_path: str, oszicar_path: str) -> float:
    """
    Extract total magnetization from OUTCAR/OSZICAR.
    
    Args:
        outcar_path: Path to OUTCAR
        oszicar_path: Path to OSZICAR
        
    Returns:
        Total magnetization in Bohr magnetons
    """
    # Check if calculation completed
    with open(outcar_path) as f:
        for line in f:
            if "General timing and accounting informations for this job:" in line:
                # Read magnetization from OSZICAR
                with open(oszicar_path) as oszi:
                    for line in oszi:
                        if 'mag=' in line:
                            mag_str = line.partition('mag=')[2]
                            mag_components = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+|\d+", mag_str)]
                            return np.sqrt(sum(m**2 for m in mag_components))
    
    raise ValueError(f"Could not extract magnetization from {outcar_path}")


def mae_function(theta, k1, k2, c):
    """MAE fitting function: E(θ) = K₁sin²θ + K₂sin⁴θ + c"""
    return k1 * (np.sin(theta)**2) + k2 * (np.sin(theta)**4) + c


def calculate_magnetic_properties(mae_ev: float,
                                  k1_mj_m3: float,
                                  volume: float,
                                  magnetization_bohr: float) -> dict:
    """
    Calculate derived magnetic properties.
    
    Args:
        mae_ev: MAE in eV
        k1_mj_m3: K1 anisotropy constant in MJ/m³
        volume: Unit cell volume in ų
        magnetization_bohr: Magnetization in Bohr magnetons
        
    Returns:
        Dictionary of calculated properties
    """
    # Convert MAE to SI units (J/m³)
    k1_si = k1_mj_m3 * 1e6  # MJ/m³ to J/m³
    
    # Magnetization in A/m
    m0_am = magnetization_bohr * (Bh / (Ang**3)) * (1 / volume)
    
    # Maximum energy product (BH)max in kJ/m³
    bh_max = 0.25 * mu0 * (m0_am**2) * 1e-3
    
    # Anisotropy field in Tesla
    mu0_ha = 2 * k1_si / (Bh * magnetization_bohr)
    
    # Hardness parameter (dimensionless)
    hardness = np.sqrt(k1_si / (mu0 * m0_am**2))
    
    return {
        'M0_MA_per_m': m0_am * 1e-6,  # MA/m
        'BH_max_kJ_m3': bh_max,
        'mu0_Ha_T': mu0_ha,
        'hardness': hardness
    }


def main():
    """Main execution function."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    calcfold_path = os.path.join(path_to_automag, 'CalcFold')
    
    # Load original structure
    original_structure = Structure.from_file(path_to_poscar)
    original_formula = original_structure.formula.replace(' ', '')
    volume = original_structure.volume
    
    # Determine expected cell type from input parameters to find correct processed file
    standardize_cell = globals().get('standardize_cell', True)
    use_primitive_cell = globals().get('use_primitive_cell', True)
    
    if not standardize_cell:
        expected_cell_type = 'original'
        cell_type = 'input_cell'
    elif use_primitive_cell:
        expected_cell_type = 'primitive'
        cell_type = 'primitive_cell'
    else:
        expected_cell_type = 'conventional'
        cell_type = 'conventional_cell'
    
    # Try to load processed structure to get actual formula and atom count
    import glob
    processed_structure = None
    processed_files = glob.glob(f'setting*_{configuration}_{expected_cell_type}.vasp')
    
    if processed_files:
        try:
            processed_structure = Structure.from_file(processed_files[0])
            formula = processed_structure.formula.replace(' ', '')
            num_atoms = len(processed_structure)
            print(f"Using processed structure: {processed_files[0]}")
            print(f"Processed formula: {formula}")
            print(f"Cell type: {expected_cell_type}")
            if formula != original_formula:
                print(f"Original formula: {original_formula}")
        except Exception as e:
            print(f"Warning: Could not load processed structure: {e}")
            formula = original_formula
            num_atoms = len(original_structure)
            processed_structure = None
    else:
        print(f"Warning: No processed structure file found matching pattern: setting*_{configuration}_{expected_cell_type}.vasp")
        print(f"Using original structure")
        formula = original_formula
        num_atoms = len(original_structure)
    
    # Get U, J values
    U = np.round(float(params['ldauu'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]), 1)
    J = np.round(float(params['ldauj'][next(i for i, x in enumerate(params['ldaul']) if x > 0)]), 1)
    encut = params['encut']
    kpts = params['kpts']
    
    # Generate comprehensive filename suffix (matching other MAE scripts)
    filename_suffix = f"{configuration}_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}_{cell_type}_{num_atoms}atoms"
    
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE ANALYSIS - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    print(f"Formula: {formula}")
    print(f"Volume: {volume:.2f} ų")
    print(f"Number of atoms: {num_atoms}")
    print(f"U = {U}, J = {J}, ENCUT = {encut}, KPTS = {kpts}")
    
    # Determine results path
    compound_dir = Path(calcfold_path) / f"{formula}{struct_suffix}"
    state_dir = compound_dir / calculator / configuration
    mae_dir = state_dir / f'mae_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}'
    
    if not mae_dir.exists():
        print(f"ERROR: MAE directory not found: {mae_dir}")
        print("Run 3_submit_mae_curve_refactored.py first.")
        return
    
    print(f"MAE directory: {mae_dir}")
    
    # Load reference energy
    reference_path = mae_dir / 'z' / 'singlepoint'
    ref_oszicar = reference_path / 'OSZICAR'
    
    if not ref_oszicar.exists():
        print(f"ERROR: Reference OSZICAR not found: {ref_oszicar}")
        return
    
    oszicar_ref = Oszicar(str(ref_oszicar))
    e_ref = float(oszicar_ref.all_energies[-1][-2])
    
    print(f"Reference energy: {e_ref:.6f} eV")
    
    # Load MAE curve results
    energies = []
    angles = []
    
    for alpha in np.linspace(0, 2 * np.pi, N_MAE + 1):
        folder_name = f'K_{kpts}_RtMAE_{np.round((alpha / np.pi) * 180, 2)}'
        folder = mae_dir / folder_name / 'singlepoint'
        
        if not folder.exists():
            print(f"Warning: Folder not found: {folder}")
            continue
        
        oszicar_path = folder / 'OSZICAR'
        outcar_path = folder / 'OUTCAR'
        
        if not oszicar_path.exists():
            print(f"Warning: OSZICAR not found in {folder}")
            continue
        
        # Check if calculation completed
        completed = False
        if outcar_path.exists():
            with open(outcar_path) as f:
                for line in f:
                    if "General timing and accounting informations for this job:" in line:
                        completed = True
                        break
        
        if completed:
            oszicar = Oszicar(str(oszicar_path))
            energy = float(oszicar.all_energies[-1][-2])
            energies.append(energy)
            angles.append(alpha)
            print(f"α = {np.round((alpha / np.pi) * 180, 2):6.2f}°  E = {energy:.6f} eV")
        else:
            print(f"Warning: Calculation not completed in {folder}")
    
    if len(energies) == 0:
        print("\nERROR: No completed MAE curve calculations found!")
        return
    
    energies = np.array(energies)
    angles = np.array(angles)
    
    # Calculate MAE
    mae_ev = np.max(energies) - np.min(energies)
    
    # Convert to MJ/m³
    mae_mj_m3_total = mae_ev * (eV / (volume * Ang**3)) * 1e-6
    mae_mj_m3_per_atom = mae_mj_m3_total / num_atoms
    
    print(f"\n{'=' * 70}")
    print("MAE RESULTS")
    print(f"{'=' * 70}")
    print(f"MAE (total):     {mae_mj_m3_total:.6f} MJ/m³")
    print(f"MAE (per atom):  {mae_mj_m3_per_atom:.6f} MJ/m³")
    print(f"MAE (eV):        {mae_ev:.6f} eV")
    
    # Fit MAE curve
    try:
        energies_relative = (energies - np.min(energies)) * (eV / (volume * Ang**3)) * 1e-6
        popt, pcov = curve_fit(mae_function, angles, energies_relative)
        k1, k2, c = popt
        
        print(f"\nFitted parameters:")
        print(f"K1 = {k1:.6f} MJ/m³")
        print(f"K2 = {k2:.6f} MJ/m³")
        
        # Calculate fitted curve
        angles_fit = np.linspace(0, 2 * np.pi, 100)
        energies_fit = mae_function(angles_fit, k1, k2, c)
        
    except Exception as e:
        print(f"\nWarning: Could not fit MAE curve: {e}")
        k1 = mae_mj_m3_total
        k2 = 0.0
        energies_fit = None
        angles_fit = None
    
    # Extract magnetization
    try:
        ref_outcar = reference_path / 'OUTCAR'
        ref_oszicar = reference_path / 'OSZICAR'
        magnetization = extract_magnetization(str(ref_outcar), str(ref_oszicar))
        
        print(f"\nMagnetization:")
        print(f"M0 = {magnetization:.4f} μB (Bohr magnetons)")
        
        # Calculate derived properties
        mag_props = calculate_magnetic_properties(
            mae_ev, k1, volume, magnetization
        )
        
        print(f"\nDerived magnetic properties:")
        print(f"M0 = {mag_props['M0_MA_per_m']:.4f} MA/m")
        print(f"(BH)max = {mag_props['BH_max_kJ_m3']:.4f} kJ/m³")
        print(f"μ₀Hₐ = {mag_props['mu0_Ha_T']:.4f} T")
        print(f"Hardness = {mag_props['hardness']:.4f}")
        
    except Exception as e:
        print(f"\nWarning: Could not calculate magnetic properties: {e}")
        magnetization = None
        mag_props = None
    
    # Create output directory
    output_dir = Path(f'outputs_{configuration}')
    output_dir.mkdir(exist_ok=True)
    
    # Save numerical results with comprehensive naming
    np.save(output_dir / f'mae_curve_alpha_{filename_suffix}.npy', angles * 180 / np.pi)
    np.save(output_dir / f'mae_curve_raw_energy_{filename_suffix}.npy', energies)
    np.save(output_dir / f'mae_curve_energy_per_atom_{filename_suffix}.npy', 
            (energies - e_ref) * (eV / (volume * Ang**3)) * 1e-6 / num_atoms)
    
    # Plot MAE curve
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Convert to degrees and MJ/m³ per atom
    angles_deg = angles * 180 / np.pi
    energies_plot = (energies - e_ref) * (eV / (volume * Ang**3)) * 1e-6 / num_atoms
    
    ax.plot(angles_deg, energies_plot, 'o-', linewidth=2, markersize=8, 
            label=f'DFT (K={kpts})')
    
    if energies_fit is not None:
        angles_fit_deg = angles_fit * 180 / np.pi
        ax.plot(angles_fit_deg, energies_fit / num_atoms, '--', linewidth=2,
                label=f'Fit: K₁={k1:.3f}, K₂={k2:.3f} MJ/m³')
    
    ax.set_xlabel('Rotation angle α (degrees)', fontsize=14)
    ax.set_ylabel('Energy (MJ/m³ per atom)', fontsize=14)
    ax.set_title(f'MAE Curve - {configuration} (U={U}, J={J}, K={kpts}, {cell_type})', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'mae_curve_{filename_suffix}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Plot saved to: {output_dir / f'mae_curve_{filename_suffix}.png'}")
    
    # Write text report
    report_file = output_dir / f'mae_curve_report_{filename_suffix}.txt'
    
    with open(report_file, 'w') as f:
        f.write(f"{'=' * 70}\n")
        f.write(f"MAE Curve Analysis Report\n")
        f.write(f"{'=' * 70}\n\n")
        
        f.write("Structure Information:\n")
        f.write(f"  Formula: {formula}\n")
        f.write(f"  Configuration: {configuration}\n")
        f.write(f"  Volume: {volume:.4f} ų\n")
        f.write(f"  Number of atoms: {num_atoms}\n\n")
        
        f.write("Calculation Parameters:\n")
        f.write(f"  U = {U} eV\n")
        f.write(f"  J = {J} eV\n")
        f.write(f"  ENCUT = {encut} eV\n")
        f.write(f"  KPTS = {kpts}\n")
        f.write(f"  Cell type = {cell_type}\n\n")
        
        f.write("MAE Results:\n")
        f.write(f"  MAE = {mae_ev:.6f} eV\n")
        f.write(f"  MAE = {mae_mj_m3_total:.6f} MJ/m³ (total)\n")
        f.write(f"  MAE = {mae_mj_m3_per_atom:.6f} MJ/m³ (per atom)\n\n")
        
        f.write("Fitted Anisotropy Constants:\n")
        f.write(f"  K1 = {k1:.6f} MJ/m³\n")
        f.write(f"  K2 = {k2:.6f} MJ/m³\n\n")
        
        if mag_props is not None:
            f.write("Magnetic Properties:\n")
            f.write(f"  M0 = {magnetization:.4f} μB (Bohr magnetons)\n")
            f.write(f"  M0 = {mag_props['M0_MA_per_m']:.4f} MA/m\n")
            f.write(f"  (BH)max = {mag_props['BH_max_kJ_m3']:.4f} kJ/m³\n")
            f.write(f"  μ₀Hₐ = {mag_props['mu0_Ha_T']:.4f} T\n")
            f.write(f"  Hardness = {mag_props['hardness']:.4f}\n\n")
        
        f.write(f"{'=' * 70}\n")
    
    print(f"✓ Report saved to: {report_file}")
    
    print(f"\n{'=' * 70}")
    print("ANALYSIS COMPLETE")
    print(f"{'=' * 70}\n")


if __name__ == '__main__':
    main()
