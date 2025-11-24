# SOLID Refactoring Implementation Guide

## What Has Been Refactored

This refactoring follows SOLID principles to improve code quality, maintainability, and testability.

### New Architecture

```
automag-1/
├── core/                           # Core domain and interfaces
│   ├── domain/
│   │   └── calculation.py          # Domain models (SRP)
│   ├── interfaces/
│   │   ├── workflow.py             # Abstract interfaces (OCP, DIP)
│   │   └── serialization.py        # Serialization interface (ISP)
│   ├── services/
│   │   └── calculation_service.py  # Business logic (SRP, DIP)
│   └── factories/
│       └── workflow_factory.py     # Object creation (OCP)
│
├── infrastructure/                 # Framework-specific implementations
│   ├── workflow/
│   │   ├── fireworks_adapter.py    # FireWorks implementation (LSP)
│   │   └── manual_adapter.py       # Manual submission (LSP)
│   ├── calculators/
│   │   └── vasp_adapter.py         # VASP-specific code (SRP, ISP)
│   └── serialization/
│       └── atoms_serializer.py     # JSON serialization (SRP)
│
└── tests/
    └── test_refactored_components.py # Unit tests
```

## SOLID Principles Applied

### 1. **Single Responsibility Principle (SRP)**

**Before:** `SubmitFirework.py` (259 lines) handled:
- Workflow submission
- Parameter validation
- File I/O
- Configuration management

**After:** Split into focused classes:
- `CalculationParameters` - manages parameters only
- `CalculationService` - orchestrates workflows only
- `FireworksWorkflowSubmitter` - handles FireWorks submission only
- `VaspInputWriter` - writes VASP inputs only

### 2. **Open/Closed Principle (OCP)**

**Before:** Adding new workflow systems or calculators required modifying existing code.

**After:** 
- `WorkflowSubmitter` interface allows new implementations without modifying existing code
- Factory pattern enables adding new submitter types
- Strategy pattern for different calculation modes

### 3. **Liskov Substitution Principle (LSP)**

**Before:** `SubmitFirework` and `SubmitManual` couldn't be used interchangeably.

**After:** Both implement `WorkflowSubmitter` interface and can be substituted:

```python
# Works with either submitter
def run_calculation(submitter: WorkflowSubmitter, atoms, config):
    return submitter.submit_single_calculation(atoms, config)

# Can use FireWorks
submitter = WorkflowSubmitterFactory.create_fireworks_submitter(launchpad)
run_calculation(submitter, atoms, config)

# Or Manual - same code works!
submitter = WorkflowSubmitterFactory.create_manual_submitter(calcfold)
run_calculation(submitter, atoms, config)
```

### 4. **Interface Segregation Principle (ISP)**

**Before:** `utilities.py` forced clients to import unrelated functionality.

**After:** Focused interfaces:
- `StructureSerializer` - only serialization
- `CalculationInputWriter` - only input writing
- `ResultParser` - only result parsing

Clients only depend on what they need.

### 5. **Dependency Inversion Principle (DIP)**

**Before:** Direct dependencies on FireWorks, hard-coded paths, global configuration.

**After:** 
- Depend on `WorkflowSubmitter` abstraction, not concrete FireWorks
- Dependencies injected through constructors
- Configuration externalized

## How to Use the Refactored Code

### Example 1: Convergence Test (New Way)

```python
from core.factories.workflow_factory import WorkflowSubmitterFactory
from core.services.calculation_service import CalculationService
from core.domain.calculation import CalculationParameters
from ase.io import read

# Load structure
atoms = read('structure.vasp')

# Create parameters (validated)
params = CalculationParameters(sigma=0.1, kpts=40)

# Create submitter using factory
submitter = WorkflowSubmitterFactory.create_from_config(
    use_fireworks=False,
    calcfold_path='CalcFold',
    calculator_command='mpirun vasp_std'
)

# Create service
service = CalculationService(submitter)

# Submit convergence test
job_ids = service.submit_convergence_test(
    atoms=atoms,
    base_params=params,
    magmoms=[5.0, 5.0],
    test_parameter='encut',
    test_values=[500, 600, 700],
    mode_name='encut'
)
```

### Example 2: Switch Between FireWorks and Manual

```python
# Same code, different submitter!

# Option 1: FireWorks
submitter = WorkflowSubmitterFactory.create_fireworks_submitter(
    launchpad_file='~/.fireworks/my_launchpad.yaml'
)

# Option 2: Manual (SLURM)
submitter = WorkflowSubmitterFactory.create_manual_submitter(
    calcfold_path='CalcFold',
    jobheader='#!/bin/bash\n#SBATCH -N 1',
    calculator_command='module load vasp && mpirun vasp_std'
)

# Service code is identical!
service = CalculationService(submitter)
```

## Migration Path

### Phase 1: Use New Code Alongside Old (Current)

The refactored code exists alongside the old code:
- Old scripts: `0_conv_tests/1_submit.py`
- New example: `0_conv_tests/1_submit_refactored.py`

Both work, allowing gradual migration.

### Phase 2: Update Scripts Gradually

Update one workflow at a time:
1. ✅ Convergence tests (example provided)
2. ⏳ Linear response
3. ⏳ Collinear search
4. ⏳ Monte Carlo
5. ⏳ MAE calculations

### Phase 3: Deprecate Old Code

Once all workflows use new code, remove old classes.

## Benefits Achieved

### ✅ Testability
- Small, focused classes easy to unit test
- Dependencies can be mocked
- See `tests/test_refactored_components.py` for examples

### ✅ Maintainability
- Clear separation of concerns
- Changes localized to specific modules
- Easier to understand and modify

### ✅ Extensibility
- Add new calculators (Quantum ESPRESSO, etc.) without modifying existing code
- Support new workflow systems (Parsl, Prefect) easily
- Add calculation types with minimal changes

### ✅ Flexibility
- Easy to switch between FireWorks and manual submission
- Configuration changes don't require code changes
- Support different HPC environments

## Running Tests

```bash
# Install pytest if needed
pip install pytest

# Run tests
cd automag-1
python -m pytest tests/test_refactored_components.py -v
```

## Next Steps

1. **Try the refactored code:**
   ```bash
   python 0_conv_tests/1_submit_refactored.py
   ```

2. **Review the interfaces:**
   - `core/interfaces/workflow.py`
   - `core/interfaces/serialization.py`

3. **Adapt for your needs:**
   - Modify `core/factories/workflow_factory.py` for your setup
   - Customize `infrastructure/workflow/manual_adapter.py` for your queue system

4. **Gradually migrate:**
   - Update one script at a time
   - Test thoroughly before deprecating old code

## Questions?

The refactored code follows industry best practices and design patterns. Each component has clear documentation explaining its purpose and how it follows SOLID principles.
