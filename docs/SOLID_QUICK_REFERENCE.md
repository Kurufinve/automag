# SOLID Refactoring - Quick Reference Card

## 🎯 Quick Start

### Run Refactored Scripts

```bash
# Convergence tests
cd 0_conv_tests && python 1_submit_refactored.py

# Linear response U
cd 1_lin_response && python 1_submit_refactored.py

# Collinear search
cd 2_coll && python 1_submit_refactored.py

# Monte Carlo coupling
cd 3_monte_carlo && python 1_submit_refactored.py
```

### Switch Workflow Systems

```python
# In any refactored script:

# FireWorks
use_fireworks = True

# Manual (SLURM/PBS)
use_fireworks = False
```

## 📦 Architecture Layers

```
core/                      # Business logic (framework-independent)
├── domain/                # Data models
├── interfaces/            # Abstract interfaces
├── services/              # Business operations
└── factories/             # Object creation

infrastructure/            # Framework implementations
├── workflow/              # FireWorks, Manual adapters
├── calculators/           # VASP, QE adapters
└── serialization/         # JSON, etc.
```

## 🔑 Key Classes

| Class | Purpose | Principle |
|-------|---------|-----------|
| `CalculationParameters` | Immutable parameters with validation | SRP |
| `WorkflowSubmitter` | Abstract workflow interface | OCP, DIP |
| `CalculationService` | Orchestrates calculations | SRP, DIP |
| `WorkflowSubmitterFactory` | Creates submitters | OCP, DIP |
| `VaspInputWriter` | Writes VASP inputs only | SRP, ISP |
| `VaspResultParser` | Parses VASP outputs only | SRP, ISP |

## 🧩 Common Patterns

### Create a Calculation

```python
from core.domain.calculation import CalculationParameters
from core.factories.workflow_factory import WorkflowSubmitterFactory
from core.services.calculation_service import CalculationService

# 1. Create parameters
params = CalculationParameters(encut=500, sigma=0.1, kpts=40)
params.validate()

# 2. Create submitter (abstraction!)
submitter = WorkflowSubmitterFactory.create_from_config(
    use_fireworks=False,
    calcfold_path='CalcFold',
    calculator_command='mpirun vasp_std'
)

# 3. Create service
service = CalculationService(submitter)

# 4. Submit
job_id = service.submit_singlepoint_calculation(
    atoms, params, magmoms=[5.0], name='test'
)
```

### Parse Results

```python
from infrastructure.calculators.vasp_adapter import VaspResultParser

parser = VaspResultParser()

# Get what you need
energy = parser.parse_energy('path/to/workdir')
magmoms = parser.parse_magnetic_moments('path/to/workdir')
converged = parser.check_convergence('path/to/workdir')
```

### Serialize Structure

```python
from infrastructure.serialization.atoms_serializer import JSONAtomsSerializer

serializer = JSONAtomsSerializer()

# Serialize
json_str = serializer.serialize(atoms)

# Deserialize
atoms_restored = serializer.deserialize(json_str)
```

## 🧪 Testing

```python
# Unit test with mocks
from unittest.mock import Mock

mock_submitter = Mock()
mock_submitter.submit_single_calculation.return_value = "job_123"

service = CalculationService(mock_submitter)
job_id = service.submit_singlepoint_calculation(...)

assert mock_submitter.submit_single_calculation.called
```

## ✅ SOLID Checklist

When writing new code:

- [ ] **SRP**: Does this class have one clear responsibility?
- [ ] **OCP**: Can I extend without modifying existing code?
- [ ] **LSP**: Are implementations truly substitutable?
- [ ] **ISP**: Are interfaces focused (not bloated)?
- [ ] **DIP**: Do I depend on abstractions, not concretions?

## 📝 Code Templates

### New Calculator Adapter

```python
from core.interfaces.workflow import CalculationInputWriter

class MyCalculatorInputWriter(CalculationInputWriter):
    def write_input_files(self, atoms, params, workdir):
        # Write inputs
        pass
    
    def set_magnetic_moments(self, atoms, magmoms):
        # Set magmoms
        return atoms
```

### New Workflow Submitter

```python
from core.interfaces.workflow import WorkflowSubmitter

class MyWorkflowSubmitter(WorkflowSubmitter):
    def submit_single_calculation(self, atoms, config):
        # Submit
        return job_id
    
    def submit_workflow(self, atoms, configs, name):
        # Submit workflow
        return workflow_id
    
    def check_status(self, job_id):
        # Check status
        return status
```

### New Service Method

```python
class CalculationService:
    def submit_my_calculation(self, atoms, params, **kwargs):
        # Validate
        params.validate()
        
        # Create config
        config = CalculationConfig(...)
        
        # Submit via abstraction
        return self._submitter.submit_single_calculation(atoms, config)
```

## 🚨 Common Pitfalls

❌ **Don't:**
```python
# Hard-coded dependency
from infrastructure.workflow.fireworks_adapter import FireworksWorkflowSubmitter
submitter = FireworksWorkflowSubmitter(...)  # Violates DIP
```

✅ **Do:**
```python
# Use factory (dependency injection)
from core.factories.workflow_factory import WorkflowSubmitterFactory
submitter = WorkflowSubmitterFactory.create_from_config(...)  # Follows DIP
```

❌ **Don't:**
```python
# Mixed responsibilities
class MyClass:
    def do_everything(self):
        self.load_data()
        self.process_data()
        self.save_data()
        self.plot_data()
        # Violates SRP
```

✅ **Do:**
```python
# Focused classes
loader = DataLoader()
processor = DataProcessor()
saver = DataSaver()
plotter = DataPlotter()

data = loader.load()
processed = processor.process(data)
saver.save(processed)
plotter.plot(processed)
```

## 📊 Before/After Comparison

### Before (Violates SOLID)
```python
if use_fireworks:
    run = SubmitFirework(...)  # Different class
    run.submit()
else:
    run = SubmitManual(...)    # Different class
    run.submit()
# Violates LSP, OCP, DIP
```

### After (Follows SOLID)
```python
submitter = factory.create(use_fireworks)  # Common interface
service = CalculationService(submitter)     # Dependency injection
job_id = service.submit(...)                # Polymorphism
# Follows all SOLID principles!
```

## 🎓 Learning Resources

- **SOLID Principles**: See main REFACTORING_GUIDE.md
- **Design Patterns**: Factory, Strategy, Adapter patterns used
- **Unit Testing**: tests/test_refactored_components.py
- **Examples**: All *_refactored.py scripts

## 💡 Pro Tips

1. **Use factories** for object creation → enables DIP
2. **Inject dependencies** via constructors → testable code
3. **Program to interfaces** not implementations → flexibility
4. **Keep classes small** (~100 lines) → maintainability
5. **Test in isolation** with mocks → fast tests

## 🔗 File Reference

| What | Where |
|------|-------|
| Interfaces | `core/interfaces/` |
| Domain models | `core/domain/` |
| Services | `core/services/` |
| Factories | `core/factories/` |
| VASP adapter | `infrastructure/calculators/vasp_adapter.py` |
| FireWorks adapter | `infrastructure/workflow/fireworks_adapter.py` |
| Manual adapter | `infrastructure/workflow/manual_adapter.py` |
| Tests | `tests/test_refactored_components.py` |
| Examples | `*_refactored.py` scripts |

---

**Remember:** The goal isn't perfect code, but code that's easy to change, test, and understand!
