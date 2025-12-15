# MAE Analysis Script Fix

## Problem

The `2_analyze_results_refactored.py` script was failing with:
```
WARNING: Could not load reference energy from z/ directory
Reference energy: 0.000000 eV
Loaded 0 results out of 231 calculations
ERROR: No results found!
```

## Root Cause

The `MAEResultsLoader` class was looking for OSZICAR files in the wrong locations:

1. **Grid calculations**: Looking for `PhTh_*/singlepoint/OSZICAR` instead of `PhTh_*/OSZICAR`
2. **Reference calculation**: The path was correct (`z/OSZICAR`), but error handling was insufficient

The mismatch occurred because:
- The old MAE workflow (from `MAE.py`) wrote results to `PhTh_*/singlepoint/`
- The refactored workflow (`1_submit_refactored.py`) writes directly to `PhTh_*/` and `z/`

## Solution

### 1. Fixed `MAEResultsLoader.load_theta_phi_results()` 

**File**: `core/services/mae_service.py`

Changed to check **multiple possible locations** for OSZICAR files:
```python
# Try multiple possible locations for OSZICAR
possible_paths = [
    folder / 'OSZICAR',  # Direct in calculation folder (refactored workflow)
    folder / 'singlepoint' / 'OSZICAR',  # Old workflow location
]

for oszicar_path in possible_paths:
    if oszicar_path.exists():
        try:
            oszicar = Oszicar(str(oszicar_path))
            energy = float(oszicar.all_energies[-1][-2])
            results[direction.name] = energy
            break  # Found it, stop searching
        except Exception as e:
            continue  # Try next path
```

This makes the loader **backward-compatible** with both old and new workflows.

### 2. Improved `MAEResultsLoader.load_reference_energy()`

**File**: `core/services/mae_service.py`

Added better error messages:
```python
if ref_path.exists():
    try:
        oszicar = Oszicar(str(ref_path))
        return float(oszicar.all_energies[-1][-2])
    except Exception as e:
        print(f"Error reading reference OSZICAR: {e}")
        pass
else:
    print(f"Reference OSZICAR not found at: {ref_path}")
```

### 3. Enhanced `2_analyze_results_refactored.py`

**File**: `4_mae/2_analyze_results_refactored.py`

#### Automatic MAE Directory Detection

The script now:
1. Checks if it's run from within the MAE directory (looks for `z/` subfolder)
2. If not, attempts to construct the path from input parameters
3. Provides helpful error messages with expected paths

```python
# Check if we're in the MAE calculation directory (should have 'z' subfolder)
if not (base_path / 'z').exists():
    # Construct expected path from input.py parameters
    mae_dir = calcfold_path / f"{formula}{struct_suffix}" / calculator / configuration / 
              f"mae_U{U:.1f}_J{J:.1f}_K{kpts_val}_EN{encut_val}"
    
    if mae_dir.exists() and (mae_dir / 'z').exists():
        base_path = mae_dir
        print(f"Found MAE directory: {base_path}")
    else:
        # Provide helpful error message
        print(f"ERROR: Could not find MAE calculation directory!")
        print(f"Expected path: {mae_dir}")
```

#### Better Error Diagnostics

When no results are found, the script now:
- Lists possible reasons
- Shows expected directory structure
- Checks and reports on first few directories

```python
if len(results) == 0:
    print("\nERROR: No results found!")
    print("\nPossible reasons:")
    print("  1. Calculations haven't completed yet")
    print("  2. OSZICAR files don't exist in PhTh_* directories")
    print("  3. Running from wrong directory")
    print(f"\nExpected structure:")
    print(f"  {base_path}/z/OSZICAR (reference)")
    print(f"  {base_path}/PhTh_0.01_0.01/OSZICAR")
    # ... check actual directories
```

## How to Use

### Option 1: Run from 4_mae directory (Recommended)

```bash
cd /path/to/automag-1/4_mae
python 2_analyze_results_refactored.py
```

The script will automatically find the MAE directory from `input.py` parameters.

### Option 2: Run from MAE directory

```bash
cd /path/to/CalcFold/Fe12O18/vasp/fm1/mae_U0.0_J0.0_K20_EN830
python /path/to/automag-1/4_mae/2_analyze_results_refactored.py
```

## Expected Directory Structure

After calculations complete, the structure should be:

```
CalcFold/Fe12O18/vasp/fm1/mae_U0.0_J0.0_K20_EN830/
├── z/
│   ├── POSCAR
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   ├── OSZICAR          ← Reference energy read from here
│   ├── OUTCAR
│   ├── WAVECAR
│   └── CHGCAR
├── PhTh_0.01_0.01/
│   ├── POSCAR
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   ├── OSZICAR          ← Grid energies read from here
│   ├── OUTCAR
│   ├── WAVECAR -> ../z/WAVECAR
│   └── CHGCAR -> ../z/CHGCAR
├── PhTh_0.01_18.0/
│   ├── ...
│   └── OSZICAR          ← Grid energies read from here
├── PhTh_0.01_36.0/
│   └── ...
...
```

## Testing

After applying these fixes, the script should:
1. Successfully load reference energy from `z/OSZICAR`
2. Load all 231 grid calculation results from `PhTh_*/OSZICAR`
3. Calculate MAE and identify easy/hard axes
4. Generate surface plot if all calculations completed
5. Save results to `mae_results_{configuration}.txt`

## Verification

Run the script and check for output like:
```
Found MAE directory: /path/to/CalcFold/.../mae_U0.0_J0.0_K20_EN830
Reference energy: -354.123456 eV
Loading grid calculation results...
Loaded 231 results out of 231 calculations
```

If you see `Loaded 0 results`, check:
1. VASP calculations actually completed (check OSZICAR files exist)
2. You're in the correct directory
3. The error diagnostic messages for specific paths
