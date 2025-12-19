# Script Porting and Debugging Analysis

## Overview

This document analyzes the differences between the original self-sufficient MAE.py script and the modular ported implementation consisting of 1_submit.py, 2_plot_results.py, 3_submit.py, 4_plot_results.py, and the common SubmitManual.py module. The goal is to identify why the results are not identical and provide solutions to ensure consistent behavior.

## Architecture Comparison

### Original MAE.py Architecture
```
graph TB
    A[MAE.py - Single Script] --> B[Parameter Definition]
    B --> C[VASP Input Generation]
    C --> D[Job Submission & Monitoring]
    D --> E[K-point Convergence Tests]
    E --> F[Collinear Reference Calculation]
    F --> G[Non-collinear Theta-Phi Grid]
    G --> H[MAE Curve Calculation]
    H --> I[Energy Fitting & Analysis]
    I --> J[Results Output]
```

### Ported Modular Architecture
```
graph TB
    A[input.py Configuration] --> B[1_submit.py]
    B --> C[VASP Theta-Phi Grid Calc]
    C --> D[2_plot_results.py]
    D --> E[Axis Identification]
    E --> F[3_submit.py]
    F --> G[VASP MAE Curve Calc]
    G --> H[4_plot_results.py]
    H --> I[Final Analysis]
    
    J[SubmitManual.py] --> B
    J --> F
```

## Critical Differences Identified

### 1. K-point Convergence Testing

**MAE.py Implementation:**
- Performs automatic k-point convergence testing with predefined values: `[0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2]`
- Automatically selects optimal k-point based on energy convergence criteria
- Tests both collinear (std) and non-collinear (ncl) calculations for each k-point

**Ported Implementation:**
- Uses fixed k-point value from `params['kpts']` in input configuration
- No automatic k-point optimization
- **Impact:** Different k-point selection can lead to different energy surfaces and MAE values

### 2. Energy Reference Selection

**MAE.py Implementation:**
```python
# Selects k-point based on convergence criteria
K_select_list = [y for x,y in zip(abs(E0_SI), Kp_read) if x<0.12 and y <0.11 and y>0.05]
K_select = max(K_select_list)
E_ref = float((Oszicar(f"./K_{K_select}_EN_{Encut_select}/OSZICAR").all_energies[-1])[-2])
```

**Ported Implementation:**
```python
# Uses z-direction reference directly
E_ref = float((Oszicar(f"{mae_dir}/z/singlepoint/OSZICAR").all_energies[-1])[-2])
```

**Impact:** Different reference energies affect all relative energy calculations.

### 3. Magnetic Moment Handling

**MAE.py Implementation:**
- Uses element-based magnetic moment assignment with predefined values
- `NAGMOM_ncl()` function creates non-collinear moments: `'0 0 5'` for magnetic atoms, `'0 0 2'` for others

**Ported Implementation:**
- Uses actual magnetic moments from previous collinear calculations
- Converts to non-collinear format: `[(0,0,m) for m in standardized_magmom]`

**Impact:** Different initial magnetic moments can affect convergence and final energies.

### 4. Directory Structure and File Organization

**MAE.py Implementation:**
- Creates directories with specific naming: `K_{K}_EN_{En}`, `PhTh_{phi}_{theta}`
- Manages all file operations internally

**Ported Implementation:**
- Uses CalcFold hierarchy: `{compound_dir}/{calculator}/{configuration}/mae_U{U}_J{J}_K{kpts}_EN{encut}`
- Relies on SubmitManual.py for directory creation and file management

### 5. Calculation Sequencing and Dependencies

**MAE.py Implementation:**
- Sequential execution with built-in job monitoring
- Automatic resubmission on failure
- Integrated convergence checking

**Ported Implementation:**
- Manual execution of separate scripts
- No automatic job monitoring between stages
- Requires user intervention between stages

## Specific Technical Issues

### Issue 1: SAXIS Vector Calculation

**MAE.py Implementation:**
```python
def RotatingAngles(mytheta, myphi):
    x = np.sin(mytheta)*np.cos(myphi)
    y = np.sin(mytheta)*np.sin(myphi)
    z = np.cos(mytheta)
    return x,y,z

# Applied directly to SAXIS
X,Y,Z = RotatingAngles(thetan, phin)
params['saxis'] = [np.round(X,3),np.round(Y,3),np.round(Z,3)]
```

**SubmitManual.py Implementation:**
```python
def RotatingAngles(mytheta, myphi):
    x = np.sin(mytheta)*np.cos(myphi)
    y = np.sin(mytheta)*np.sin(myphi)
    z = np.cos(mytheta)
    return x,y,z

# Applied in mae_theta_phi mode
X,Y,Z = RotatingAngles(thetan, phin)
params['saxis'] = [np.round(X,3),np.round(Y,3),np.round(Z,3)]
```

**Verification:** Both implementations use identical SAXIS calculations.

### Issue 2: Energy Unit Conversion

**MAE.py Implementation:**
```python
# Converts to MJ/m³
energy_SI = (energy - E_ref) * ((10.0**-6.0)/volume) * (eV/(Ang**3)) * (1.0/n_atoms)
```

**Ported Implementation:**
```python
# Same conversion in plot scripts
energy_SI = (energy - E_ref) * ((10.0**-6.0)/readPOSCAR()[0]) * (eV/(Ang**3)) * (1.0/readPOSCAR()[3])
```

**Verification:** Energy conversion formulas are identical.

### Issue 3: MAE Coordinate System Definition

**MAE.py Implementation:**
```python
MAE_x = np.array(min_vec[:])
if max(abs(np.array(np.cross(min_vec, max_vec))))>0.0001:
    MAE_z = np.array(NormedCross(min_vec, max_vec))
else:
    # Uses random vector approach
    Ran_vec = np.random.rand(3)
    MAE_z = np.array(np.cross(min_vec, Ran_vec))/np.linalg.norm(np.cross(min_vec, Ran_vec))
MAE_y = np.array(np.cross(MAE_z, MAE_x))/np.linalg.norm(np.cross(MAE_z, MAE_x))
```

**2_plot_results.py Implementation:**
```python
MAE_x = np.array(min_vec[:])
if max(abs(np.array(np.cross(min_vec, max_vec))))>0.0001:
    MAE_z = np.array(NormedCross(min_vec, max_vec))
else:
    # Uses saved random vector
    if os.path.isfile(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy'):
        MAE_z = np.load(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy')
MAE_y = np.array(np.cross(MAE_z, MAE_x))/np.linalg.norm(np.cross(MAE_z, MAE_x))
```

**Issue:** Random vector generation may differ between runs, affecting coordinate system consistency.

## Root Cause Analysis

### Primary Issues

1. **Structure Standardization Differences**
   - MAE.py: Uses original structure
   - Ported: Uses primitive cell from pymatgen
   - **Solution:** Ensure consistent structure handling

2. **Reference Energy Inconsistency**
   - Different calculation paths lead to different reference energies
   - **Solution:** Use same reference calculation approach

3. **Magnetic Moment Source**
   - MAE.py: Element-based defaults
   - Ported: Previous calculation results
   - **Solution:** Ensure consistent magnetic moment values

### Secondary Issues

1. **Random Vector Seeding**
   - Non-deterministic coordinate system generation
   - **Solution:** Use fixed seed or consistent vector selection

2. **File Path Differences**
   - Different directory structures affect file access
   - **Solution:** Normalize path handling

3. **Calculation Parameter Inheritance**
   - Some VASP parameters may differ between implementations
   - **Solution:** Verify all VASP parameters match exactly

## Recommended Solutions

### 1. Consistent Structure Handling

Ensure both implementations use the same structure:

```python
# Option 1: Disable standardization in ported version
standardize_cell = False

# Option 2: Apply same standardization in MAE.py
if standardize_cell:
    structure = structure.get_primitive_structure(tolerance=0.2, use_site_props=True)
```

### 2. Reference Energy Synchronization

Use identical reference calculation approach:

```python
# Ensure same WAVECAR/CHGCAR source
# Use same k-point for reference calculation
# Apply same calculation parameters
```

### 3. Deterministic Coordinate System

Fix random vector generation:

```python
# Use fixed seed
np.random.seed(42)

# Or use deterministic fallback vector
if max(abs(np.array(np.cross(min_vec, max_vec)))) <= 0.0001:
    # Use predetermined orthogonal vector instead of random
    fallback_vec = np.array([1, 0, 0]) if abs(min_vec[0]) < 0.9 else np.array([0, 1, 0])
    MAE_z = np.array(np.cross(min_vec, fallback_vec))
    MAE_z = MAE_z / np.linalg.norm(MAE_z)
```

### 4. Parameter Verification Matrix

Create verification system to ensure parameter consistency:

| Parameter | MAE.py Value | Ported Value | Status |
|-----------|--------------|--------------|---------|
| SAXIS | Calculated | Calculated | ✅ |
| MAGMOM | Element-based | History-based | ⚠️ |
| LDAU* | Fixed | Configurable | ⚠️ |

## Testing and Validation Strategy

### 1. Parameter Alignment Test

Run both implementations with identical:
- Structure file
- K-point value
- Magnetic moments
- All VASP parameters

### 2. Intermediate Results Comparison

Compare at each stage:
- Energy values for each (θ, φ) point
- MAE coordinate system vectors
- Reference energy values
- Final anisotropy constants

### 3. Numerical Precision Analysis

Check for:
- Floating-point precision differences
- Array ordering differences
- Unit conversion consistency

### 4. File I/O Verification

Ensure:
- Same POSCAR reading
- Identical VASP input generation
- Consistent output parsing

## Implementation Priority

1. **High Priority:** Reference energy alignment and magnetic moment source synchronization
2. **Medium Priority:** Consistent structure handling and deterministic coordinate system
3. **Low Priority:** Directory structure normalization

This systematic approach should identify and resolve the discrepancies between the original MAE.py script and the ported modular implementation.

## Practical Implementation Solutions

### Solution 1: Fix Deterministic Coordinate System

**Issue:** Random vector generation in MAE coordinate system causes inconsistent results.

**File to modify:** `4_mae/2_plot_results.py`

```python
# Replace the random vector section (around line 180-190)
# Original problematic code:
# if max(abs(np.array(np.cross(min_vec, max_vec))))>0.0001:
#     MAE_z = np.array(NormedCross(min_vec, max_vec))
# else:
#     if os.path.isfile(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy') == False:
#         Ran_vec = np.random.rand(3)
#         Ran_vec = Ran_vec/np.linalg.norm(Ran_vec)
#         MAE_z = np.array(np.cross(min_vec, Ran_vec))/np.linalg.norm(np.cross(min_vec, Ran_vec))
#         np.save(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy', MAE_z)
#     if os.path.isfile(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy') == True:
#         MAE_z = np.load(f'./outputs_{configuration}/MAE_Ran_vec_U{U:.1f}_J{J:.1f}_K{kpts}_EN{encut}.npy')

# Fixed deterministic code:
MAE_x = np.array(min_vec[:])
if max(abs(np.array(np.cross(min_vec, max_vec)))) > 0.0001:
    MAE_z = np.array(NormedCross(min_vec, max_vec))
else:
    # Use deterministic fallback vector (same as MAE.py)
    np.random.seed(42)  # Fixed seed for reproducibility
    Ran_vec = np.random.rand(3)
    Ran_vec = Ran_vec / np.linalg.norm(Ran_vec)
    MAE_z = np.array(np.cross(min_vec, Ran_vec)) / np.linalg.norm(np.cross(min_vec, Ran_vec))
MAE_y = np.array(np.cross(MAE_z, MAE_x)) / np.linalg.norm(np.cross(MAE_z, MAE_x))
```

### Solution 2: Ensure Consistent Magnetic Moment Handling

**Issue:** Different magnetic moment sources between MAE.py and ported version.

**File to modify:** `4_mae/1_submit.py` and `4_mae/3_submit.py`

```python
# Add option to use element-based magnetic moments like MAE.py
# Add this function to both files:

def get_element_based_magmoms(atoms, magnetic_elements={'Fe': 5.0}):
    """
    Generate magnetic moments based on element types (MAE.py style)
    """
    magmoms = []
    for symbol in atoms.get_chemical_symbols():
        if symbol in magnetic_elements:
            magmoms.append(magnetic_elements[symbol])
        else:
            magmoms.append(2.0)  # Default for non-magnetic elements
    return [(0, 0, m) for m in magmoms]

# Add configuration option in input files:
# use_element_based_magmoms = True  # Set to True to match MAE.py behavior
# magnetic_elements = {'Fe': 5.0}  # Define magnetic moments per element

# Then modify the magnetic moment assignment:
if 'use_element_based_magmoms' in globals() and use_element_based_magmoms:
    standardized_magmom = get_element_based_magmoms(ase_structure, magnetic_elements)
else:
    # Use existing collinear calculation results
    standardized_magmom = [(0,0,m) for m in standardized_magmom]
```

### Solution 3: Consistent Structure Handling

**Issue:** Different structure standardization between implementations.

**File to modify:** `4_mae/1_submit.py` and `4_mae/3_submit.py`

```python
# Add option to disable structure standardization
# Add this configuration in input.py:
# use_primitive_structure = False  # Set to False to match MAE.py behavior

# Modify structure processing:
if 'use_primitive_structure' in globals() and use_primitive_structure:
    standardized_structure = pmg_structure.get_primitive_structure(tolerance=0.2, use_site_props=True)
else:
    # Use original structure (MAE.py behavior)
    standardized_structure = pmg_structure

# Update structure output filename:
if 'use_primitive_structure' in globals() and use_primitive_structure:
    standardized_structure_name = f'setting{setting:03d}_{configuration}_standardized.vasp'
else:
    standardized_structure_name = f'setting{setting:03d}_{configuration}_original.vasp'
```

### Solution 4: Reference Energy Synchronization

**Issue:** Different reference energy calculations.

**File to modify:** `common/SubmitManual.py`

```python
# In the mae_theta_phi and mae_curve modes, ensure reference calculation matches MAE.py
# Modify the check_z_dir() function:

def check_z_dir():
    chgcar_found = False
    # First check if we should use the same reference as MAE.py
    # Look for K-point optimized calculation if available
    if 'use_kpoint_reference' in globals() and use_kpoint_reference:
        # Try to find k-point optimized reference calculation
        for k_val in [0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2]:
            ref_dir = f'../K_{k_val}_EN_{params["encut"]}'
            if os.path.isdir(ref_dir) and os.path.isfile(f'{ref_dir}/CHGCAR'):
                print(f"Using k-point reference from {ref_dir}")
                os.system(f'cp {ref_dir}/CHGCAR {z_dir}/singlepoint/')
                os.system(f'cp {ref_dir}/WAVECAR {z_dir}/singlepoint/')
                chgcar_found = True
                break
    
    # Original logic as fallback
    if not chgcar_found:
        # ... existing check_z_dir logic ...
        pass
    
    return chgcar_found
```

### Solution 5: VASP Parameter Consistency

**Issue:** Ensure all VASP parameters match between implementations.

**File to modify:** `4_mae/input_template.py`

```python
# Add MAE.py compatible parameter set
mae_py_compatible_params = {
    'voskown': 1,
    'lnoncollinear': True,
    'lsorbit': True,
    'gga_compat': False,
    'lcharg': False,
    'lwave': False,
    'icharg': 11,
    'istart': 1,
    # Add other MAE.py specific parameters
}

# Option to use MAE.py parameter set
use_mae_py_params = True  # Set to True for exact MAE.py compatibility

if use_mae_py_params:
    params.update(mae_py_compatible_params)
```

### Solution 6: Enhanced Debugging and Validation

**Create new file:** `4_mae/validate_consistency.py`

```python
#!/usr/bin/env python3
"""
Validation script to compare MAE.py and ported implementation results
"""

import numpy as np
import os
from pymatgen.io.vasp.outputs import Oszicar

def compare_energy_surfaces(mae_py_dir, ported_dir, tolerance=1e-6):
    """
    Compare energy surfaces from both implementations
    """
    differences = []
    
    # Compare theta-phi grid energies
    for phi_deg in np.linspace(0, 360, 21):
        for theta_deg in np.linspace(0, 180, 11):
            # MAE.py format
            mae_py_folder = f"{mae_py_dir}/PhTh_{phi_deg:.2f}_{theta_deg:.2f}"
            # Ported format  
            ported_folder = f"{ported_dir}/PhTh_{phi_deg:.2f}_{theta_deg:.2f}/singlepoint"
            
            if os.path.exists(f"{mae_py_folder}/OSZICAR") and os.path.exists(f"{ported_folder}/OSZICAR"):
                mae_py_energy = float(Oszicar(f"{mae_py_folder}/OSZICAR").all_energies[-1][-2])
                ported_energy = float(Oszicar(f"{ported_folder}/OSZICAR").all_energies[-1][-2])
                
                diff = abs(mae_py_energy - ported_energy)
                if diff > tolerance:
                    differences.append((phi_deg, theta_deg, diff, mae_py_energy, ported_energy))
    
    return differences

def validate_mae_vectors(mae_py_output, ported_output):
    """
    Compare MAE coordinate system vectors
    """
    # Parse vectors from output files
    # Implementation depends on output format
    pass

if __name__ == "__main__":
    # Usage: python validate_consistency.py
    print("Validating MAE implementation consistency...")
    
    # Add validation logic here
    differences = compare_energy_surfaces("./mae_py_results", "./ported_results")
    
    if differences:
        print(f"Found {len(differences)} energy differences above tolerance")
        for diff in differences[:10]:  # Show first 10
            print(f"Phi={diff[0]:.1f}, Theta={diff[1]:.1f}: diff={diff[2]:.2e} eV")
    else:
        print("Energy surfaces match within tolerance")
```

These practical solutions address the key issues identified in the analysis and provide specific code modifications to ensure consistent results between the original MAE.py script and the ported modular implementation.