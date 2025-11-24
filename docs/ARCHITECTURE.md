# Automag-1 Refactored Architecture

## Overview

The refactored architecture follows Clean Architecture and SOLID principles, separating concerns into distinct layers.

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                        PRESENTATION LAYER                        │
│  (Submission Scripts: *_refactored.py)                          │
│                                                                  │
│  • 1_submit_refactored.py (convergence, linear response, etc.) │
│  • User input via input.py files                               │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                      APPLICATION LAYER                           │
│  (Use Cases & Orchestration)                                    │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  CalculationService                               │          │
│  │  • submit_convergence_test()                     │          │
│  │  • submit_singlepoint_calculation()              │          │
│  │  • submit_perturbation_workflow()                │          │
│  └──────────────────────────────────────────────────┘          │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Factories                                        │          │
│  │  • WorkflowSubmitterFactory                      │          │
│  │    - create_fireworks_submitter()                │          │
│  │    - create_manual_submitter()                   │          │
│  └──────────────────────────────────────────────────┘          │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                        DOMAIN LAYER                              │
│  (Business Logic - Framework Independent)                       │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Domain Models                                    │          │
│  │  • CalculationParameters                         │          │
│  │  • MagneticConfiguration                         │          │
│  │  • CalculationMode (enum)                        │          │
│  └──────────────────────────────────────────────────┘          │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Interfaces (Abstractions)                       │          │
│  │  • WorkflowSubmitter (interface)                 │          │
│  │  • CalculationInputWriter (interface)            │          │
│  │  • ResultParser (interface)                      │          │
│  │  • StructureSerializer (interface)               │          │
│  └──────────────────────────────────────────────────┘          │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                    INFRASTRUCTURE LAYER                          │
│  (Framework-Specific Implementations)                           │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Workflow Adapters                                │          │
│  │  • FireworksWorkflowSubmitter                    │          │
│  │  • ManualJobSubmitter                            │          │
│  └──────────────────────────────────────────────────┘          │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Calculator Adapters                              │          │
│  │  • VaspInputWriter                               │          │
│  │  • VaspResultParser                              │          │
│  └──────────────────────────────────────────────────┘          │
│                                                                  │
│  ┌──────────────────────────────────────────────────┐          │
│  │  Serialization                                    │          │
│  │  • JSONAtomsSerializer                           │          │
│  └──────────────────────────────────────────────────┘          │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│                      EXTERNAL SYSTEMS                            │
│                                                                  │
│  • FireWorks + MongoDB        • VASP Calculator                 │
│  • SLURM/PBS Queue System     • File System                     │
└─────────────────────────────────────────────────────────────────┘
```

## Dependency Flow

```
High Level (Abstract) ───────────────► Low Level (Concrete)
     ▲                                       │
     │                                       │
     │  Dependency Inversion                 │
     │  (Depend on Abstractions)            │
     │                                       │
     └───────────────────────────────────────┘
            (Injection via Factories)
```

## Component Relationships

```
┌─────────────────────────────────────────────────────────────────┐
│                   Client Code (Scripts)                          │
└───────────────────────────┬─────────────────────────────────────┘
                            │
                            │ uses
                            ▼
        ┌──────────────────────────────────────┐
        │   WorkflowSubmitterFactory           │
        │   (creates appropriate submitter)    │
        └─────────────┬────────────────────────┘
                      │
                      │ returns
                      ▼
        ┌──────────────────────────────────────┐
        │   WorkflowSubmitter (interface)       │
        │   • submit_single_calculation()      │
        │   • submit_workflow()                │
        │   • check_status()                   │
        └─────────────┬────────────────────────┘
                      │
                      │ implemented by
         ┌────────────┴────────────┐
         ▼                         ▼
┌────────────────────┐   ┌────────────────────┐
│ FireworksWorkflow  │   │  ManualJob         │
│ Submitter          │   │  Submitter         │
│                    │   │                    │
│ (uses FireWorks)   │   │ (uses SLURM/PBS)   │
└────────────────────┘   └────────────────────┘
```

## Data Flow Example: Convergence Test

```
1. User runs script
   └─► 1_submit_refactored.py

2. Load configuration
   └─► input.py (poscar_file, params, test_values)

3. Create parameters
   └─► CalculationParameters(**params)
       └─► validate()  ✓

4. Create workflow submitter via factory
   └─► WorkflowSubmitterFactory.create_from_config()
       └─► Returns: FireworksWorkflowSubmitter OR ManualJobSubmitter

5. Create service with dependency injection
   └─► CalculationService(submitter)

6. Submit convergence test
   └─► service.submit_convergence_test(...)
       └─► For each test value:
           └─► Create CalculationConfig
           └─► submitter.submit_single_calculation()
               └─► (FireWorks OR Manual - same interface!)

7. Job execution
   └─► VASP runs
   └─► Results collected

8. Analysis
   └─► 2_plot_results.py
```

## Key Design Patterns Used

### 1. Factory Pattern
```python
# WorkflowSubmitterFactory
submitter = WorkflowSubmitterFactory.create_from_config(
    use_fireworks=False,
    calcfold_path='CalcFold'
)
# Returns concrete implementation based on config
```

### 2. Strategy Pattern
```python
# Different workflow strategies
class WorkflowSubmitter(ABC):  # Strategy interface
    @abstractmethod
    def submit_single_calculation(...): pass

class FireworksWorkflowSubmitter(WorkflowSubmitter):  # Concrete strategy
    def submit_single_calculation(...): ...

class ManualJobSubmitter(WorkflowSubmitter):  # Concrete strategy
    def submit_single_calculation(...): ...
```

### 3. Adapter Pattern
```python
# VaspInputWriter adapts ASE Vasp calculator to our interface
class VaspInputWriter(CalculationInputWriter):
    def write_input_files(self, atoms, params, workdir):
        calc = Vasp(atoms=atoms, directory=workdir, **params)
        calc.write_input(atoms)
```

### 4. Dependency Injection
```python
# Service receives dependencies via constructor
class CalculationService:
    def __init__(self, submitter: WorkflowSubmitter):
        self._submitter = submitter  # Injected dependency
```

## SOLID Principles in Action

### Single Responsibility Principle
```
Before:
  SubmitFirework (259 lines)
    ├─ Workflow submission
    ├─ Parameter validation
    ├─ File I/O
    └─ Configuration

After:
  CalculationParameters    → Parameter validation only
  CalculationService       → Workflow orchestration only
  FireworksWorkflowSubmitter → FireWorks submission only
  VaspInputWriter          → VASP input generation only
```

### Open/Closed Principle
```
Want to add Quantum ESPRESSO support?

❌ Before: Modify SubmitFirework.py and SubmitManual.py

✓ After: 
  1. Create QEInputWriter(CalculationInputWriter)
  2. Create QEResultParser(ResultParser)
  3. Update factory
  4. Done! No changes to existing code
```

### Liskov Substitution Principle
```python
# Any WorkflowSubmitter can be substituted
def run_calc(submitter: WorkflowSubmitter):
    return submitter.submit_single_calculation(...)

# Works with either implementation
run_calc(FireworksWorkflowSubmitter(...))  ✓
run_calc(ManualJobSubmitter(...))          ✓
```

### Interface Segregation Principle
```
Before:
  utilities.py → Everything in one file

After:
  StructureSerializer     → Serialization only
  CalculationInputWriter  → Input writing only
  ResultParser            → Result parsing only
  
Clients depend only on what they need!
```

### Dependency Inversion Principle
```
High-level modules (CalculationService) depend on abstractions
(WorkflowSubmitter interface), not on low-level modules
(FireworksWorkflowSubmitter).

Dependencies are injected via constructors, not hard-coded.
```

## File Structure

```
automag-1/
├── core/                       # Domain & Application Layer
│   ├── domain/
│   │   └── calculation.py      # Business entities
│   ├── interfaces/
│   │   ├── workflow.py         # Abstract interfaces
│   │   └── serialization.py
│   ├── services/
│   │   └── calculation_service.py  # Business logic
│   └── factories/
│       └── workflow_factory.py # Object creation
│
├── infrastructure/             # Infrastructure Layer
│   ├── workflow/
│   │   ├── fireworks_adapter.py   # FireWorks impl
│   │   └── manual_adapter.py      # Manual impl
│   ├── calculators/
│   │   └── vasp_adapter.py        # VASP impl
│   └── serialization/
│       └── atoms_serializer.py    # JSON impl
│
├── tests/                      # Test Layer
│   └── test_refactored_components.py
│
├── *_submit_refactored.py      # Presentation Layer
│   (in each workflow directory)
│
└── docs/                       # Documentation
    ├── ARCHITECTURE.md (this file)
    ├── REFACTORED_SCRIPTS_GUIDE.md
    └── SOLID_QUICK_REFERENCE.md
```

## Testing Strategy

```
Unit Tests (tests/test_refactored_components.py)
├── Domain Models
│   ├── test_parameters_validation
│   └── test_magnetic_configuration
├── Serialization
│   └── test_serialize_deserialize_roundtrip
├── Service Layer
│   ├── test_submit_convergence_test
│   └── test_submit_singlepoint_with_recalc
└── Interface Compliance
    └── test_submitters_have_same_interface

All tests use mocks for dependencies (Dependency Injection benefit!)
```

## Extension Points

### Adding New Features

1. **New Calculator Support**
   - Implement `CalculationInputWriter`
   - Implement `ResultParser`
   - Register in factory

2. **New Workflow System**
   - Implement `WorkflowSubmitter`
   - Register in factory

3. **New Calculation Type**
   - Add method to `CalculationService`
   - Use existing infrastructure

4. **New Serialization Format**
   - Implement `StructureSerializer`
   - Swap in where needed

## Benefits Summary

| Aspect | Before | After |
|--------|--------|-------|
| **Testability** | Hard (tight coupling) | Easy (DI, mocks) |
| **Flexibility** | Rigid (hard-coded) | Flexible (pluggable) |
| **Maintainability** | Difficult (big classes) | Easy (small, focused) |
| **Extensibility** | Requires modifications | Add new implementations |
| **Understanding** | Complex (mixed concerns) | Clear (separation) |

## Conclusion

The refactored architecture provides:
- ✅ Clear separation of concerns
- ✅ Testable components
- ✅ Pluggable infrastructure
- ✅ Framework independence
- ✅ Future-proof design

**The architecture is production-ready and demonstrates industry best practices for maintainable scientific software.**
