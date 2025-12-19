# CalcFold Directory Naming Convention Update

## Overview

The folder naming convention within the CalcFold directory structure has been updated to provide better organization and clarity about the cell types used in calculations.

## Old Pattern
```
CalcFold/{formula}{struct_suffix}/{calculator}/{configuration}/
```

Example:
```
CalcFold/Fe12O18/vasp/fm1/
```

## New Pattern
```
CalcFold/{reduced_formula}{struct_suffix}/{formula}_{cell_type}/{calculator}/{configuration}/
```

Example:
```
CalcFold/Fe2O3/Fe12O18_primitive_cell/vasp/fm1/
```

## Key Components

### reduced_formula
- Obtained from `structure.composition.reduced_formula.replace(' ', '')`
- Represents the simplest formula unit (e.g., Fe₂O₃ instead of Fe₁₂O₁₈)
- Groups related calculations of different supercell sizes together

### formula
- Full chemical formula using ASE's metal mode: `atoms.get_chemical_formula(mode='metal')`
- Represents the actual cell size used in the calculation

### cell_type
One of four values:
- `input_cell` - Structure used as-is from input file (no standardization)
- `primitive_cell` - Standardized to primitive cell
- `conventional_cell` - Standardized to conventional cell
- `supercell` - Supercell expansion (implementation for usage to be developed later)

## Modified Files

### Core Components
1. **common/SubmitManual.py**
   - Added `cell_type` parameter to `__init__` (default: 'input_cell')
   - Updated directory creation to use new pattern
   - Imports pymatgen.Structure for reduced_formula calculation

2. **infrastructure/workflow/manual_adapter.py**
   - Added `cell_type` parameter to `__init__` (default: 'input_cell')
   - Updated `submit_single_calculation()` to use new pattern
   - Updated `submit_workflow()` to use new pattern
   - Imports pymatgen.Structure for reduced_formula calculation

3. **core/factories/workflow_factory.py**
   - Added `cell_type` parameter to `create_manual_submitter()` (default: 'input_cell')
   - Passes cell_type to ManualJobSubmitter

### MAE Module (4_mae/)
1. **mae_utils.py**
   - Updated `construct_mae_directory_path()` to accept `reduced_formula` parameter
   - Updated path construction to use new pattern
   - Updated documentation and examples

2. **1_submit_refactored.py**
   - Calculates reduced_formula from processed structure
   - Uses new directory pattern for MAE calculations

3. **2_plot_results.py**
   - Loads pymatgen Structure to get reduced_formula
   - Detects cell_type from filename (primitive/conventional/input_cell)
   - Uses new directory pattern for reading results

4. **4_plot_results.py**
   - Same updates as 2_plot_results.py

### Collinear Module (2_coll/)
1. **1_submit_refactored.py**
   - Determines cell_type based on standardize_cell and use_primitive_cell flags
   - Passes cell_type to workflow factory

### Convergence Tests (0_conv_tests/)
1. **1_submit_refactored.py**
   - Uses 'input_cell' as default cell_type
   - Passes cell_type to workflow factory

### Linear Response (1_lin_response/)
1. **1_submit_refactored.py**
   - Uses 'input_cell' as default cell_type
   - Passes cell_type to workflow factory

### Monte Carlo (3_monte_carlo/)
- No changes needed (doesn't use CalcFold structure)

### Common Utilities
- **write_output.py, write_charges.py, etc.**
  - No changes needed (write flat files directly to CalcFold, not using directory structure)

## Migration Guide

### For New Calculations
Simply use the updated scripts - the new naming convention will be applied automatically based on your `standardize_cell` and `use_primitive_cell` settings.

### For Existing Calculations
The old directory structure will continue to work. However, new calculations will use the new pattern. If you want to reorganize existing calculations:

1. Identify the cell type used:
   - If `standardize_cell=False` → use `input_cell`
   - If `standardize_cell=True` and `use_primitive_cell=True` → use `primitive_cell`
   - If `standardize_cell=True` and `use_primitive_cell=False` → use `conventional_cell`

2. Calculate reduced_formula:
   ```python
   from pymatgen.core.structure import Structure
   structure = Structure.from_file('your_structure.vasp')
   reduced_formula = structure.composition.reduced_formula.replace(' ', '')
   ```

3. Rename directories:
   ```bash
   # Old: CalcFold/Fe12O18/vasp/fm1/
   # New: CalcFold/Fe2O3/Fe12O18_primitive_cell/vasp/fm1/
   mv CalcFold/Fe12O18 CalcFold/Fe2O3/Fe12O18_primitive_cell
   ```

## Benefits

1. **Clear Cell Type Identification**: The cell_type is now explicit in the directory name
2. **Better Organization**: Related calculations (different supercell sizes) are grouped under reduced_formula
3. **Consistency**: Uniform naming across all calculation types
4. **Future-Proof**: Supercell option included for future development
5. **Backward Compatible**: Old scripts can be updated gradually; flat output files remain unchanged

## Cell Type Determination Logic

### For refactored scripts:
```python
if not standardize_cell:
    cell_type = 'input_cell'
elif use_primitive_cell:
    cell_type = 'primitive_cell'
else:
    cell_type = 'conventional_cell'
```

### For plotting/analysis scripts:
```python
# Detect from filename
cell_type = 'input_cell'  # Default
if 'primitive' in filename:
    cell_type = 'primitive_cell'
elif 'conventional' in filename:
    cell_type = 'conventional_cell'
```

## Example Directory Structures

### Convergence Tests
```
CalcFold/
└── Fe2O3/
    └── Fe12O18_input_cell/
        └── vasp/
            └── convergence/
                ├── encut/
                └── kgrid/
```

### Linear Response
```
CalcFold/
└── Fe2O3/
    └── Fe12O18_input_cell/
        └── vasp/
            └── perturbations_Mn/
                ├── singlepoint/
                ├── 0.1/
                └── 0.2/
```

### Collinear Magnetic Search
```
CalcFold/
└── Fe2O3/
    └── Fe12O18_primitive_cell/
        └── vasp/
            ├── fm1/
            ├── afm1/
            └── afm2/
```

### MAE Calculations
```
CalcFold/
└── Fe2O3/
    └── Fe12O18_primitive_cell/
        └── vasp/
            └── fm1/
                └── mae_U5.2_J0.0_K20_EN830_primitive_cell_10atoms/
                    ├── z/
                    ├── PhTh_0.0_0.0/
                    ├── PhTh_0.0_18.0/
                    └── ...
```

## Notes

- The `struct_suffix` parameter remains available for additional naming flexibility
- Output files (e.g., `Fe12O18_singlepoint_vasp.txt`) in the root CalcFold directory maintain their flat structure
- The supercell option is included in the implementation but its specific usage will be developed later
