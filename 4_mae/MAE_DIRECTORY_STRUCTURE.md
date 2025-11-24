# MAE Directory Structure and Submission Scripts

## Overview

The refactored MAE submission script ([`1_submit_refactored.py`](1_submit_refactored.py)) now creates a proper directory structure and generates helper SLURM submission scripts for easy job management.

## Directory Structure

The script creates directories following this pattern:

```
CalcFold/
└── {formula}{struct_suffix}/
    └── {calculator}/
        └── {configuration}/
            └── mae_U{U}_J{J}_K{kpts}_EN{encut}/
                ├── PhTh_{phi}_{theta}/        # Grid calculations
                │   ├── POSCAR
                │   ├── INCAR
                │   ├── KPOINTS
                │   └── POTCAR (symlink)
                ├── PhTh_...
                ├── ...
                └── submit_mae_grid.sh         # Generated submission script
```

### Example:

For a Fe₁₂O₁₈ system with configuration `fm1`, U=4.0, J=0.0, kpts=0.08, ENCUT=600:

```
CalcFold/
└── Fe12O18/
    └── vasp/
        └── fm1/
            └── mae_U4.0_J0.0_K0.08_EN600/
                ├── PhTh_0.0_0.0/
                ├── PhTh_0.0_18.0/
                ├── PhTh_0.0_36.0/
                ├── ...
                ├── PhTh_342.0_162.0/
                └── submit_mae_grid.sh
```

## Multiple K-points and ENCUT Values

If `kpts` and/or `encut` are provided as lists in `input.py`, the script creates separate directories for each combination:

```python
# input.py
params = {
    'kpts': [0.08, 0.1, 0.12],
    'encut': [600, 700],
    # ... other parameters
}
```

This creates 6 directories:
```
mae_U4.0_J0.0_K0.08_EN600/
mae_U4.0_J0.0_K0.08_EN700/
mae_U4.0_J0.0_K0.1_EN600/
mae_U4.0_J0.0_K0.1_EN700/
mae_U4.0_J0.0_K0.12_EN600/
mae_U4.0_J0.0_K0.12_EN700/
```

Each directory contains its own grid calculations and submission script.

## Generated Submission Scripts

### Parallel Mode (`parallel_over_configurations = True`)

Creates `submit_mae_grid.sh` that submits all grid calculations as separate SLURM jobs:

```bash
#!/bin/bash
# MAE Grid Parallel Submission Script
# This script submits all MAE grid calculations in parallel

echo 'Submitting 231 MAE grid calculations in parallel...'

# Submit PhTh_0.0_0.0
cd /path/to/CalcFold/.../mae_U4.0_J0.0_K0.08_EN600/PhTh_0.0_0.0
cat > job.sh << 'EOF'
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J mae_PhTh_0.0_0.0
#SBATCH -o /path/to/.../PhTh_0.0_0.0/slurm-%j.out
#SBATCH -e /path/to/.../PhTh_0.0_0.0/slurm-%j.err

source .venv/bin/activate
mpirun vasp_std
deactivate
EOF
sbatch job.sh
cd /path/to/mae_U4.0_J0.0_K0.08_EN600

# Submit PhTh_0.0_18.0
...

echo 'All jobs submitted!'
echo 'Monitor with: squeue -u $USER'
```

**Usage:**
```bash
cd CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600
bash submit_mae_grid.sh
```

### Sequential Mode (`parallel_over_configurations = False`)

Creates `submit_mae_grid_sequential.sh` that runs all calculations sequentially in a single SLURM job:

```bash
#!/bin/bash
# MAE Grid Sequential Submission Script
# This script runs MAE grid calculations sequentially in a single SLURM job

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 168:00:00
#SBATCH -J mae_grid_sequential
#SBATCH -o /path/to/mae_U4.0_J0.0_K0.08_EN600/mae_sequential-%j.out
#SBATCH -e /path/to/mae_U4.0_J0.0_K0.08_EN600/mae_sequential-%j.err

source .venv/bin/activate

echo 'Running calculation for PhTh_0.0_0.0...'
cd /path/to/.../PhTh_0.0_0.0
mpirun vasp_std

# Check if calculation completed successfully
if [ $? -ne 0 ]; then
    echo 'ERROR: Calculation failed for PhTh_0.0_0.0'
    exit 1
fi
echo 'Completed PhTh_0.0_0.0'

echo 'Running calculation for PhTh_0.0_18.0...'
cd /path/to/.../PhTh_0.0_18.0
mpirun vasp_std

...

deactivate
echo 'All MAE grid calculations completed!'
```

**Usage:**
```bash
cd CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600
sbatch submit_mae_grid_sequential.sh
```

## Input Parameters

### Required in `input.py`:

```python
# Structure and configuration
poscar_file = 'Fe12O18.vasp'
configuration = 'fm1'  # From collinear calculations
struct_suffix = ''  # Optional suffix

# MAE grid resolution
Nph = 20  # Number of phi points (0 to 2π)
Nth = 10  # Number of theta points (0 to π)

# VASP parameters
params = {
    'kpts': 0.08,  # Can be a list: [0.08, 0.1, 0.12]
    'encut': 600,  # Can be a list: [600, 700]
    'xc': 'PBE',
    'prec': 'Accurate',
    'sigma': 0.1,
    'ismear': 0,
    'ldaul': [2, -1],  # Fe: d-electrons, O: no +U
    'ldauu': [4.0, 0.0],
    'ldauj': [0.0, 0.0],
    'ldau': True,
    'ldautype': 2,
    'lorbit': 11,
    'voskown': 1,
    'lnoncollinear': True,  # Required for MAE
    'lsorbit': True,        # Required for MAE
    'gga_compat': False,
    # ... other VASP parameters
}

# Job submission settings
calculator = 'vasp'
parallel_over_configurations = True  # or False for sequential

jobheader = """#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00"""

calculator_command = "mpirun vasp_std"
environment_activate = "source .venv/bin/activate"
environment_deactivate = "deactivate"
```

## Workflow Steps

### 1. Submit Grid Calculations

```bash
cd 4_mae
python 1_submit_refactored.py
```

**What it does:**
- Loads configuration from collinear calculations (`2_coll` results)
- Standardizes structure to primitive cell
- Creates directory structure in `CalcFold/`
- For each kpts/encut combination:
  - Creates `mae_U{U}_J{J}_K{kpts}_EN{encut}/` directory
  - Creates grid of calculation directories (PhTh_*)
  - Writes VASP input files (INCAR with SAXIS, POSCAR, KPOINTS)
  - Generates submission script

**Output:**
```
CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600/
  ├── PhTh_0.0_0.0/ ... PhTh_342.0_162.0/  (231 directories)
  └── submit_mae_grid.sh
  
setting001_fm1_standardized.vasp
fm1_mae_config.txt
```

### 2. Submit Jobs

**Parallel mode:**
```bash
cd CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600
bash submit_mae_grid.sh
```

**Sequential mode:**
```bash
cd CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600
sbatch submit_mae_grid_sequential.sh
```

### 3. Monitor Jobs

**Parallel:**
```bash
squeue -u $USER
watch squeue -u $USER
```

**Sequential:**
```bash
tail -f mae_sequential-<jobid>.out
```

### 4. Analyze Results

After calculations complete:

```bash
cd 4_mae
python 2_analyze_results_refactored.py
```

## Advantages of This Structure

### 1. **Organized Directory Hierarchy**
- Clear path from formula → calculator → configuration → MAE parameters
- Easy to find and compare different parameter sets
- Follows VASP calculation best practices

### 2. **Multiple Parameter Sets**
- Test convergence with different k-points and ENCUT
- All results organized in parallel directories
- Easy comparison of MAE vs. computational cost

### 3. **Flexible Job Submission**
- **Parallel mode**: Fast if queue allows many jobs
- **Sequential mode**: Single long job if queue limits apply
- One-click submission via generated scripts

### 4. **Error Handling**
- Sequential mode stops on first error
- Each calculation has individual SLURM output files
- Easy to identify and rerun failed calculations

### 5. **Reproducibility**
- Configuration saved in `fm1_mae_config.txt`
- All VASP inputs preserved
- Directory names encode all parameters

## Configuration File Example

The script creates `{configuration}_mae_config.txt`:

```
Configuration: fm1
Structure: setting001_fm1_standardized.vasp
Formula: Fe12O18
Grid size: 10×20
NCL magmoms: [(0, 0, 5.0), (0, 0, 5.0), ..., (0, 0, -2.0)]
Kpts values: [0.08, 0.1, 0.12]
Encut values: [600, 700]
U = 4.0, J = 0.0
```

## Troubleshooting

### Issue: Directory already exists

**Solution:** The script uses `mkdir(parents=True, exist_ok=True)`, so existing directories are not overwritten. Delete manually if you want to regenerate:

```bash
rm -rf CalcFold/Fe12O18/vasp/fm1/mae_U4.0_J0.0_K0.08_EN600
```

### Issue: POTCAR not found

**Solution:** Ensure `VASP_PP_PATH` environment variable is set:

```bash
export VASP_PP_PATH=/path/to/vasp/pseudopotentials
```

### Issue: Too many jobs in parallel mode

**Solution:** Use sequential mode or submit in batches:

```python
parallel_over_configurations = False
```

Or modify the generated script to submit in chunks.

### Issue: Sequential job times out

**Solution:** Increase wall time in jobheader:

```python
jobheader = """#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 168:00:00"""  # 7 days
```

## Best Practices

1. **Start with coarse grid**: Test with `Nph=10, Nth=5` first
2. **Test one k-point**: Use single value initially, then expand
3. **Check convergence**: Compare MAE for different k-points/ENCUT
4. **Monitor disk space**: Each PhTh_* directory can be several MB
5. **Use sequential mode**: For long calculations or queue limits
6. **Save important files**: Keep OSZICAR, OUTCAR, final structures

## Performance Considerations

### Grid Size vs. Computational Cost

| Grid | Calculations | Parallel Jobs | Sequential Time (est.) |
|------|--------------|---------------|------------------------|
| 5×10 | 66 | 66 jobs | ~66 × 2h = 5.5 days |
| 10×20 | 231 | 231 jobs | ~231 × 2h = 19 days |
| 20×40 | 861 | 861 jobs | ~861 × 2h = 72 days |

**Recommendation:** Start with 5×10 or 10×20 for production calculations.

### Disk Space

- Each PhTh_* directory: ~5-50 MB (depends on system size)
- Total for 10×20 grid: ~1-10 GB
- With multiple k-points/ENCUT: multiply accordingly

---

**The directory structure and submission scripts make MAE calculations organized, reproducible, and easy to manage!**
