# Refactored Submission Scripts - Usage Guide

This document explains how to use the refactored submission scripts that follow SOLID principles.

## 📁 Available Refactored Scripts

| Original Script | Refactored Script | Purpose |
|----------------|-------------------|---------|
| `0_conv_tests/1_submit.py` | [`0_conv_tests/1_submit_refactored.py`](../0_conv_tests/1_submit_refactored.py) | Convergence tests (ENCUT, k-grid) |
| `1_lin_response/1_submit.py` | [`1_lin_response/1_submit_refactored.py`](../1_lin_response/1_submit_refactored.py) | Linear response U calculation |
| `2_coll/1_submit.py` | [`2_coll/1_submit_refactored.py`](../2_coll/1_submit_refactored.py) | Collinear magnetic search |
| `3_monte_carlo/1_coupling_constants.py` | [`3_monte_carlo/1_submit_refactored.py`](../3_monte_carlo/1_submit_refactored.py) | Coupling constants calculation |

## 🎯 Key Improvements

### SOLID Principles Applied

1. **Single Responsibility Principle (SRP)**
   - Each class has one clear purpose
   - `CalculationService` - orchestrates workflows
   - `CouplingConstantsCalculator` - computes coupling constants only
   - `ResultsLoader` - loads results only
   - `ResultsPlotter` - plots results only

2. **Open/Closed Principle (OCP)**
   - Easy to switch between FireWorks and Manual submission
   - Can add new calculators without modifying existing code
   - Factory pattern allows extension

3. **Liskov Substitution Principle (LSP)**
   - `FireworksWorkflowSubmitter` and `ManualJobSubmitter` are interchangeable
   - Both implement `WorkflowSubmitter` interface

4. **Interface Segregation Principle (ISP)**
   - Focused interfaces: `StructureSerializer`, `ResultParser`, `CalculationInputWriter`
   - Clients only depend on what they need

5. **Dependency Inversion Principle (DIP)**
   - Depend on abstractions via `WorkflowSubmitter` interface
   - Dependencies injected through factories
   - No hard-coded workflow system

## 📖 Usage Examples

### 1. Convergence Tests

```bash
cd 0_conv_tests
python 1_submit_refactored.py
```

**Key Features:**
- Automatically validates parameters
- Easy to switch between FireWorks and manual submission
- Clear separation between configuration and execution

**Configuration (input.py):**
```python
poscar_file = 'Fe12O18.vasp'
mode = 'encut'  # or 'kgrid'

params = {
    'sigma': 0.1,
    'kpts': 40,
    # ... other VASP parameters
}

# Optional: specify test values
encut_values = range(500, 1010, 10)
```

### 2. Linear Response U Calculation

```bash
cd 1_lin_response
python 1_submit_refactored.py
```

**Key Features:**
- Automatic POTCAR file swapping and restoration
- Validates dummy atom configuration
- Handles perturbation workflow automatically

**Configuration (input.py):**
```python
poscar_file = 'Fe12O18.vasp'
dummy_atom = 'Zn'
dummy_position = 0
perturbations = [0.1, 0.2, 0.3, 0.4, 0.5]

params = {
    'encut': 500,
    'sigma': 0.1,
    # ... other parameters
}
```

### 3. Collinear Magnetic Search

```bash
cd 2_coll
python 1_submit_refactored.py
```

**Key Features:**
- Clean separation of enumlib logic from submission
- Reusable `MagneticConfigurationGenerator` class
- Automatic classification (FM/AFM/FiM)

**Configuration (input.py):**
```python
poscar_file = 'Fe12O18.vasp'
supercell_size = 2

spin_values = {
    'Fe': [5.0],  # Single high-spin value
    # 'Fe': [5.0, 1.0],  # High-spin and low-spin
}

params = {
    'encut': 500,
    # ... other parameters
}
```

### 4. Monte Carlo Coupling Constants

```bash
cd 3_monte_carlo
python 1_submit_refactored.py
```

**Key Features:**
- Dedicated classes for calculation, loading, and plotting
- Clean data flow and error handling
- Automatic result validation

**Configuration (input.py):**
```python
poscar_file = 'Fe12O18.vasp'
configuration = 'fm1'  # From collinear search
cutoff_radius = 6.0
control_group_size = 0.2
append_coupling_constants = True
```

## 🔄 Switching Between FireWorks and Manual Submission

All refactored scripts support both workflow systems through a simple configuration change:

```python
# Option 1: Use FireWorks
use_fireworks = True

# Option 2: Use Manual Submission (SLURM/PBS)
use_fireworks = False
calculator_command = "module load vasp/6.4.3 && mpirun vasp_std"
jobheader = """#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 24:00:00"""
```

The **same code works with both!** This demonstrates the Liskov Substitution Principle.

## 🧪 Testing

Each refactored component can be unit tested independently:

```bash
# Run all tests
python -m pytest tests/test_refactored_components.py -v

# Run specific test
python -m pytest tests/test_refactored_components.py::TestCalculationService -v
```

## 📊 Comparison: Old vs New

### Before (Old Code)

```python
# 1_lin_response/1_submit.py
if use_fireworks:
    run = SubmitFirework(path_to_poscar, mode='perturbations', 
                         fix_params=params, pert_values=perturbations,
                         magmoms=configuration, dummy_atom=dummy_atom, 
                         dummy_position=dummy_position)
    run.submit()
else:
    run = SubmitManual(path_to_poscar, mode='perturbations', 
                       fix_params=params, pert_values=perturbations,
                       magmoms=configuration, name='perturbations', 
                       dummy_atom=dummy_atom, dummy_position=dummy_position,
                       calculator=calculator, jobheader=jobheader, 
                       calculator_command=calculator_command,
                       environment_activate=environment_activate, 
                       environment_deactivate=environment_deactivate)
    run.submit()
```

**Issues:**
- Violates DRY (Don't Repeat Yourself)
- Hard to test
- Tight coupling to specific implementations
- Mixed responsibilities

### After (Refactored Code)

```python
# 1_lin_response/1_submit_refactored.py
# Create submitter using factory
submitter = WorkflowSubmitterFactory.create_from_config(
    use_fireworks=use_fireworks,
    launchpad_file=launchpad_file if use_fireworks else None,
    calcfold_path=calcfold_path,
    calculator_command=calculator_command
)

# Create service
service = CalculationService(submitter)

# Submit workflow
workflow_id = service.submit_perturbation_workflow(
    atoms=atoms,
    base_params=base_params,
    magmoms=configuration,
    perturbations=perturbations,
    dummy_atom=dummy_atom,
    dummy_position=dummy_position,
    atom_ucalc=atom_ucalc,
    workflow_name=f'perturbations_{dummy_atom}'
)
```

**Benefits:**
- ✅ DRY - no code duplication
- ✅ Testable - can mock submitter
- ✅ Flexible - easy to add new workflow systems
- ✅ Clear - single responsibility per class

## 🛠️ Extending the Architecture

### Adding a New Calculator (e.g., Quantum ESPRESSO)

```python
# infrastructure/calculators/qe_adapter.py
from core.interfaces.workflow import CalculationInputWriter

class QEInputWriter(CalculationInputWriter):
    def write_input_files(self, atoms, params, workdir):
        # Write QE input files
        pass
    
    def set_magnetic_moments(self, atoms, magmoms):
        # Set starting magnetization
        pass
```

Then use it:

```python
from infrastructure.calculators.qe_adapter import QEInputWriter

input_writer = QEInputWriter()
submitter = ManualJobSubmitter(input_writer, ...)
```

### Adding a New Workflow System (e.g., Parsl)

```python
# infrastructure/workflow/parsl_adapter.py
from core.interfaces.workflow import WorkflowSubmitter

class ParslWorkflowSubmitter(WorkflowSubmitter):
    def submit_single_calculation(self, atoms, config):
        # Use Parsl to submit
        pass
```

## 📝 Migration Checklist

When migrating from old scripts to refactored ones:

- [ ] Review your `input.py` configuration
- [ ] Test with a small calculation first
- [ ] Verify output file locations are correct
- [ ] Check job submission to your queue system
- [ ] Compare results with old implementation
- [ ] Update any downstream analysis scripts

## ⚠️ Important Notes

1. **Backward Compatibility**: Old scripts still work! Both versions coexist.

2. **POTCAR Management**: Linear response script automatically handles POTCAR swapping.

3. **Directory Structure**: Refactored scripts create the same directory structure as originals.

4. **Result Files**: Output format matches original scripts for compatibility.

## 🐛 Troubleshooting

### Issue: "Module not found" error

**Solution:** Ensure you're running from the automag-1 root directory:
```bash
cd /path/to/automag-1
python 1_lin_response/1_submit_refactored.py
```

### Issue: FireWorks LaunchPad not found

**Solution:** Check your LaunchPad configuration:
```python
# Verify path in script or set explicitly
launchpad_file = '/path/to/my_launchpad.yaml'
```

### Issue: VASP parameters not recognized

**Solution:** Use `CalculationParameters` to validate:
```python
from core.domain.calculation import CalculationParameters

params = CalculationParameters(**your_params)
params.validate()  # Will raise error if invalid
```

## 📚 Further Reading

- [REFACTORING_GUIDE.md](../REFACTORING_GUIDE.md) - Complete refactoring overview
- [tests/test_refactored_components.py](../tests/test_refactored_components.py) - Unit test examples
- Architecture diagrams and design patterns

## 🤝 Contributing

When adding new features:

1. Follow SOLID principles
2. Add unit tests
3. Update this documentation
4. Maintain backward compatibility

## 📞 Support

For questions or issues with the refactored code:
1. Check this guide first
2. Review test examples
3. Compare with original implementation
4. Submit an issue with details

---

**Note:** The refactored architecture is designed to be maintainable and extensible. Take time to understand the design patterns used, as they will make future modifications much easier!
