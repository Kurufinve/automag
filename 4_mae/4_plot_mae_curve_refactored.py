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

# Import MAE utilities
# Import MAE utilities
from mae_utils import (
    load_processed_structure,
    extract_hubbard_uj_from_params,
    extract_convergence_params_from_params,
    construct_mae_directory_path,
    construct_config_filename
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


def calculate_magnetic_properties(k1_mj_m3: float,
                                  volume: float,
                                  magnetization_bohr: float,
                                  mae_ev: float = None) -> dict:
    """
    Calculate derived magnetic properties.
    
    Args:
        k1_mj_m3: K1 anisotropy constant in MJ/m³ (per cell)
        volume: Unit cell volume in ų
        magnetization_bohr: Magnetization in Bohr magnetons
        mae_ev: Total MAE in eV (per cell). MAE.py derives the anisotropy field
            and hardness from the *total* measured MAE, not from the fitted K1,
            so pass it here to reproduce the reference values.
        
    Returns:
        Dictionary of calculated properties
    """
    # Debug output
    print(f"\n[DEBUG] calculate_magnetic_properties inputs:")
    print(f"  K1: {k1_mj_m3:.6e} MJ/m³")
    print(f"  Volume: {volume:.6e} ų")
    print(f"  Magnetization: {magnetization_bohr:.6e} μB")
    
    # Convert K1 to SI units (J/m³)
    k1_si = k1_mj_m3 * 1e6  # MJ/m³ to J/m³
    print(f"  K1 (SI): {k1_si:.6e} J/m³")
    
    # Magnetization in A/m
    # Convert volume from ų to m³
    volume_m3 = volume * (Ang**3)  # ų to m³
    print(f"  Volume (m³): {volume_m3:.6e} m³")
    
    # Convert magnetization from Bohr magnetons to A/m
    # Bohr magneton: μB = 9.274 × 10⁻²⁴ J/T = 9.274 × 10⁻²⁴ A·m²
    # M (A/m) = (Total magnetic moment in μB × μB) / Volume in m³
    m0_am = magnetization_bohr * Bh / volume_m3  # A/m
    print(f"  M0: {m0_am:.6e} A/m")
    
    # Maximum energy product (BH)max in kJ/m³
    bh_max = 0.25 * mu0 * (m0_am**2) * 1e-3
    
    # Anisotropy field, matching MAE.py:
    #   mu0xHa = 2*(K1*eV)/(Bh*M0)
    # where MAE.py's K1 is the total MAE in eV and M0 is in Bohr magnetons,
    # i.e. it is the ratio 2*K1/M0 expressed in Tesla. The fitted K1 from the
    # sin²/sin⁴ model is used only when an explicit MAE (eV) is not supplied.
    if mae_ev is not None:
        k1_for_field_ev = mae_ev
    else:
        k1_for_field_ev = k1_si * volume_m3 / eV
    mu0_ha = 2 * (k1_for_field_ev * eV) / (Bh * magnetization_bohr)  # T
    print(f"  μ₀Ha: {mu0_ha:.6e} T")
    
    # Hardness parameter (dimensionless), also based on the total MAE (eV).
    # MAE.py: sqrt((K1_eV*eV/(V*Ang**3))/(mu0*M0_A_per_m**2))
    hardness = np.sqrt(((k1_for_field_ev * eV) / volume_m3) / (mu0 * m0_am**2))
    
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
    
    # Load processed structure using utility function
    structure_path, processed_structure, expected_cell_type, cell_type = load_processed_structure(
        configuration=configuration,
        standardize_cell=globals().get('standardize_cell', True),
        use_primitive_cell=globals().get('use_primitive_cell', True),
        fallback_to_legacy=False,
        verbose=True
    )
    
    # Use the volume from the processed structure if available (the one actually used in calculations)
    # Otherwise fall back to original structure
    if processed_structure is not None:
        formula = processed_structure.formula.replace(' ', '')
        reduced_formula = processed_structure.composition.reduced_formula.replace(' ', '')
        num_atoms = len(processed_structure)
        volume = processed_structure.volume  # Use volume of structure actually calculated
        if formula != original_formula:
            print(f"Original formula: {original_formula}")
            print(f"Using processed structure with volume: {volume:.2f} ų")
    else:
        print(f"Using original structure")
        formula = original_formula
        reduced_formula = original_structure.composition.reduced_formula.replace(' ', '')
        num_atoms = len(original_structure)
        volume = original_structure.volume
    
    # Extract parameters using utility functions
    U, J = extract_hubbard_uj_from_params(params)
    kpts, encut = extract_convergence_params_from_params(params, use_first=True)
    
    # Generate comprehensive filename suffix (matching other MAE scripts)
    filename_suffix = f"{configuration}_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}_{cell_type}_{num_atoms}atoms"
    
    print(f"\n{'=' * 70}")
    print(f"MAE CURVE ANALYSIS - Configuration: {configuration}")
    print(f"{'=' * 70}\n")
    print(f"Formula: {formula}")
    print(f"Volume: {volume:.2f} ų")
    print(f"Number of atoms: {num_atoms}")
    print(f"U = {U}, J = {J}, ENCUT = {encut}, KPTS = {kpts}")
    
    # Construct MAE directory path using utility function
    mae_dir = construct_mae_directory_path(
        path_to_automag=path_to_automag,
        formula=formula,
        reduced_formula=reduced_formula,
        calculator=calculator,
        configuration=configuration,
        U=U, J=J,
        kpts_val=kpts,
        encut_val=encut,
        cell_type=cell_type,
        n_atoms=num_atoms,
        struct_suffix=struct_suffix
    )
    
    if not mae_dir.exists():
        print(f"ERROR: MAE directory not found: {mae_dir}")
        print("Run 1_submit_refactored.py first.")
        return
    
    print(f"MAE directory: {mae_dir}")
    
    # Reference energies.
    #
    # MAE.py deliberately uses *two* references:
    #   - `z`            : SCF reference for the theta-phi grid (K = K_select)
    #   - `z_mae_curve`  : SCF reference for the MAE curve  (K = smallest K)
    # The curve reference can differ from the grid one because it is computed at
    # a different KSPACING, so substituting `z` for both shifts the curve
    # vertically. Fall back to `z` only if `z_mae_curve` was never produced.
    def _read_ref_energy(folder: Path):
        for candidate in (folder / 'OSZICAR', folder / 'singlepoint' / 'OSZICAR'):
            if candidate.exists():
                return float(Oszicar(str(candidate)).all_energies[-1][-2]), candidate
        return None, None

    e_ref, ref_oszicar = _read_ref_energy(mae_dir / 'z')
    if e_ref is None:
        print(f"ERROR: Reference OSZICAR not found under: {mae_dir / 'z'}")
        return
    reference_path = ref_oszicar.parent
    print(f"Grid reference energy (z/):           {e_ref:.6f} eV")

    e_ref_curve, ref_curve_oszicar = _read_ref_energy(mae_dir / 'z_mae_curve')
    if e_ref_curve is None:
        print("WARNING: No z_mae_curve/ reference found; using the grid z/ reference.")
        e_ref_curve, ref_curve_oszicar = e_ref, ref_oszicar
    else:
        print(f"MAE curve reference energy (z_mae_curve/): {e_ref_curve:.6f} eV")
    curve_reference_path = ref_curve_oszicar.parent
    
    # Load MAE curve results from existing RtMAE folders.
    # Script 3 tags each folder with its curve mode (RtMAE_full_{angle} or
    # RtMAE_seg_{angle}) because the easy->hard segment and the full 360° circle
    # share the same angle values. Only one mode may be plotted at a time; if
    # several are present we pick the one with the most points and say so.
    import glob
    import re
    
    print(f"\nSearching for RtMAE folders in: {mae_dir}")
    rtmae_folders = glob.glob(str(mae_dir / 'RtMAE_*'))
    
    if not rtmae_folders:
        print(f"\nERROR: No RtMAE folders found in {mae_dir}")
        print("Run 3_submit_mae_curve_refactored.py first to create MAE curve calculations.")
        return
    
    # Extract mode and angle from folder name: RtMAE_{mode}_{angle}, with the
    # untagged legacy form RtMAE_{angle} still accepted.
    all_folders = []
    for folder_path in rtmae_folders:
        folder_name = Path(folder_path).name
        match = re.match(r'RtMAE_(full|seg)_([0-9.]+)$', folder_name)
        if match:
            mode, angle_deg = match.group(1), float(match.group(2))
        else:
            match = re.match(r'RtMAE_([0-9.]+)$', folder_name)
            if not match:
                continue
            mode, angle_deg = 'legacy', float(match.group(1))
        all_folders.append((mode, angle_deg, Path(folder_path)))
    
    modes = {}
    for mode, _, _ in all_folders:
        modes[mode] = modes.get(mode, 0) + 1
    if len(modes) > 1:
        print(f"WARNING: RtMAE folders from multiple curve modes present: {modes}")
    selected_mode = max(modes, key=modes.get) if modes else None
    print(f"Using curve mode: {selected_mode} ({modes.get(selected_mode, 0)} folders)")

    # (angle, folder_path) for the selected mode, sorted by angle
    folder_data = [(angle_deg, folder_path)
                   for mode, angle_deg, folder_path in all_folders
                   if mode == selected_mode]
    folder_data.sort(key=lambda x: x[0])
    
    print(f"Found {len(folder_data)} RtMAE folders:")
    for angle_deg, folder in folder_data:
        print(f"  - {folder.name} ({angle_deg}°)")
    
    # Load energies from sorted folders
    energies = []
    angles = []
    
    for angle_deg, folder in folder_data:
        # Convert angle to radians
        angle_rad = angle_deg * np.pi / 180.0
        
        # Look for OSZICAR in folder or folder/singlepoint
        oszicar_path = folder / 'OSZICAR'
        outcar_path = folder / 'OUTCAR'
        
        if not oszicar_path.exists():
            # Try singlepoint subdirectory
            oszicar_path = folder / 'singlepoint' / 'OSZICAR'
            outcar_path = folder / 'singlepoint' / 'OUTCAR'
        
        if not oszicar_path.exists():
            print(f"Warning: OSZICAR not found in {folder.name}")
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
            angles.append(angle_rad)
            print(f"α = {angle_deg:6.1f}°  E = {energy:.6f} eV")
        else:
            print(f"Warning: Calculation not completed in {folder.name}")
    
    if len(energies) == 0:
        print("\nERROR: No completed MAE curve calculations found!")
        return
    
    energies = np.array(energies)
    angles = np.array(angles)
    
    # Calculate MAE
    mae_ev = np.max(energies) - np.min(energies)
    
    # Convert to MJ/m³
    volume_conversion = (eV / (volume * Ang**3)) * 1e-6  # eV -> MJ/m³ for the cell
    mae_mj_m3_total = mae_ev * volume_conversion
    mae_mj_m3_per_atom = mae_mj_m3_total / num_atoms
    
    print(f"\n{'=' * 70}")
    print("MAE RESULTS")
    print(f"{'=' * 70}")
    print(f"MAE (total):     {mae_mj_m3_total:.6f} MJ/m³")
    print(f"MAE (per atom):  {mae_mj_m3_per_atom:.6f} MJ/m³")
    print(f"MAE (eV):        {mae_ev:.6f} eV")
    
    # Energies used for the plot are relative to the MAE-curve reference
    # (z_mae_curve in MAE.py), keeping a single baseline so the fit curve is
    # not vertically offset from the data points.
    energies_rel_mj_m3 = (energies - e_ref_curve) * volume_conversion
    
    # Fit MAE curve.
    # MAE.py fits the sin²/sin⁴ model to the curve shifted so that its minimum
    # is zero: (E_alpha - min(E_alpha)). Using the reference-energy baseline
    # instead changes the fitted K1/K2 because the constant `c` then absorbs
    # the offset.
    energies_fit_baseline = (energies - np.min(energies)) * volume_conversion
    try:
        popt, pcov = curve_fit(mae_function, angles, energies_fit_baseline)
        k1, k2, c = popt
        
        print(f"\nFitted parameters (baseline: min(E_alpha) = 0):")
        print(f"K1 = {k1:.6f} MJ/m³ (per cell volume)")
        print(f"K2 = {k2:.6f} MJ/m³ (per cell volume)")
        
        # If K1 from fit is too small or negative, use MAE-based estimate instead
        if k1 <= 0 or abs(k1) < 1e-6:
            print(f"Warning: Fitted K1 is too small or negative ({k1:.6e}). Using MAE-based estimate.")
            k1 = mae_mj_m3_total  # Use total MAE as K1 estimate
            print(f"Using K1 (from MAE) = {k1:.6f} MJ/m³")
        
        # Fitted curve, shifted back onto the plotting baseline so that it
        # overlays the data points.
        angles_fit = np.linspace(0, 2 * np.pi, 100)
        offset_mj_m3 = (np.min(energies) - e_ref_curve) * volume_conversion
        energies_fit = mae_function(angles_fit, k1, k2, c) + offset_mj_m3
        
    except Exception as e:
        print(f"\nWarning: Could not fit MAE curve: {e}")
        k1 = mae_mj_m3_total
        k2 = 0.0
        energies_fit = None
        angles_fit = None
        print(f"Using K1 (from MAE) = {k1:.6f} MJ/m³")
    
    # Extract magnetization from the MAE-curve reference, as MAE.py does.
    try:
        ref_outcar = curve_reference_path / 'OUTCAR'
        ref_oszicar = curve_reference_path / 'OSZICAR'
        magnetization = extract_magnetization(str(ref_outcar), str(ref_oszicar))
        
        print(f"\nMagnetization:")
        print(f"M0 = {magnetization:.4f} μB (Bohr magnetons)")
        
        # Calculate derived properties. The anisotropy field and hardness use
        # the *total* MAE in eV, matching MAE.py.
        mag_props = calculate_magnetic_properties(
            k1, volume, magnetization, mae_ev=mae_ev
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
    # Generate config filename using utility function
    config_filename = construct_config_filename(
        configuration=configuration,
        U=U, J=J,
        kpts_val=kpts,
        encut_val=encut,
        cell_type=cell_type,
        n_atoms=num_atoms
    )

    # Extract filename suffix for output files (remove '_mae_config' prefix and '.txt' extension)
    filename_suffix = config_filename.replace(
        '_mae_config', '').replace('.txt', '')
    # filename_suffix = f"{configuration}_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}_{cell_type}_{num_atoms}atoms"
    output_dir = Path(f'outputs_{filename_suffix}')
    output_dir.mkdir(exist_ok=True)
    
    # Save numerical results with comprehensive naming
    np.save(output_dir / f'mae_curve_alpha_{filename_suffix}.npy', angles * 180 / np.pi)
    np.save(output_dir / f'mae_curve_raw_energy_{filename_suffix}.npy', energies)
    np.save(output_dir / f'mae_curve_energy_{filename_suffix}.npy', energies_rel_mj_m3)
    np.save(output_dir / f'mae_curve_energy_per_atom_{filename_suffix}.npy',
            energies_rel_mj_m3 / num_atoms)
    
    # Plot MAE curve
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Convert to degrees and MJ/m³ per atom
    angles_deg = angles * 180 / np.pi
    energies_plot = energies_rel_mj_m3 / num_atoms
    
    ax.plot(angles_deg, energies_plot, 'o-', linewidth=2, markersize=8, 
            label=f'DFT (K={kpts})')
    
    if energies_fit is not None:
        angles_fit_deg = angles_fit * 180 / np.pi
        # Same baseline (relative to reference energy) and same per-atom
        # normalization as the data points above.
        ax.plot(angles_fit_deg, energies_fit / num_atoms, '--', linewidth=2,
                label=f'Fit: K₁={k1:.3f}, K₂={k2:.3f} MJ/m³ (per cell)')
    
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
