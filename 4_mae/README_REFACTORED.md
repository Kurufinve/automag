# MAE (Magnetocrystalline Anisotropy Energy) - Refactored Workflow

## 📋 Overview

The MAE workflow has been completely refactored following SOLID principles. The original monolithic [`MAE.py`](MAE.py) (1086 lines) has been split into:

1. **Core Services** (`core/services/mae_service.py`) - Reusable MAE components
2. **Four Workflow Scripts** - Focused, single-purpose scripts

This refactoring makes the code:
- ✅ **Maintainable** - Each component has one clear responsibility
- ✅ **Testable** - Classes can be tested in isolation
- ✅ **Extensible** - Easy to add new MAE analysis methods
- ✅ **Reusable** - Components work across different projects

## 🎯 Refactored Scripts

| Script | Purpose | Lines | Status |
|--------|---------|-------|--------|
| [`1_submit_refactored.py`](1_submit_refactored.py) | Submit theta-phi grid calculations | 268 | ✅ Complete |
| [`2_analyze_results_refactored.py`](2_analyze_results_refactored.py) | Analyze grid, find easy/hard axes | 185 | ✅ Complete |
| [`3_submit_mae_curve_refactored.py`](3_submit_mae_curve_refactored.py) | Submit MAE curve calculations | 151 | ✅ Complete |
| [`4_plot_mae_curve_refactored.py`](4_plot_mae_curve_refactored.py) | Plot curve, calculate properties | 390 | ✅ Complete |

**Total refactored code:** ~994 lines (vs. 1086 lines original)
**Key improvement:** Code is now organized into focused, testable modules

## 🏗️ Architecture

### Core Services (`core/services/mae_service.py`)

Five focused classes following Single Responsibility Principle:

```python
@dataclass
class MAEDirection:
    """Represents a magnetization direction (θ, φ, SAXIS vector)"""
    theta: float
    phi: float
    saxis: Tuple[float, float, float]

class MAEDirectionGenerator:
    """Generates theta-phi grids and rotation curves"""
    @staticmethod
    def generate_theta_phi_grid(n_theta, n_phi) -> List[MAEDirection]
    @staticmethod
    def generate_mae_curve(n_points, easy_axis, hard_axis) -> List[MAEDirection]

class MAEAnalyzer:
    """Analyzes MAE results"""
    @staticmethod
    def find_easy_hard_axes(directions, energies) -> Tuple[MAEDirection, MAEDirection]
    @staticmethod
    def calculate_mae(min_energy, max_energy) -> float
    @staticmethod
    def calculate_mae_per_volume(mae, volume) -> float
    @staticmethod
    def fit_mae_curve(angles, energies) -> Tuple[K1, K2, c]

class MAEResultsLoader:
    """Loads OSZICAR/OUTCAR results"""
    def load_theta_phi_results(directions) -> Dict[str, float]
    def load_reference_energy(reference_dir) -> float

class MAEPlotter:
    """Creates visualizations"""
    @staticmethod
    def plot_theta_phi_surface(theta_grid, phi_grid, energy_grid)
    @staticmethod
    def plot_mae_curve(angles, energies, fitted_energies)
```

### Calculation Service Integration

The `CalculationService` class orchestrates MAE workflows:

```python
class CalculationService:
    def submit_mae_theta_phi_grid(
        atoms, base_params, magmoms, n_theta, n_phi
    ) -> List[str]:
        """Submit grid of non-collinear calculations"""
        directions = MAEDirectionGenerator.generate_theta_phi_grid(n_theta, n_phi)
        # Submit each direction with SAXIS set appropriately
        
    def submit_mae_curve(
        atoms, base_params, magmoms, n_points, easy_axis, hard_axis
    ) -> List[str]:
        """Submit rotation curve from easy to hard axis"""
        directions = MAEDirectionGenerator.generate_mae_curve(n_points, easy_axis, hard_axis)
        # Submit each direction
```

## 📖 Usage Guide

### Step 1: Submit Theta-Phi Grid

Explores all magnetization directions on a spherical grid.

```bash
cd 4_mae
python 1_submit_refactored.py
```

**What it does:**
1. Loads magnetic configuration from collinear calculations (`2_coll` results)
2. Standardizes structure to primitive cell
3. Converts collinear magnetic moments to non-collinear format
4. Generates (Nθ+1) × (Nφ+1) magnetization directions
5. Submits non-collinear DFT calculations with SAXIS for each direction

**Input (`input.py`):**
```python
poscar_file = 'Fe12O18.vasp'
configuration = 'fm1'  # From collinear search

# Grid resolution
Nph = 20  # Phi points: 0 to 2π
Nth = 10  # Theta points: 0 to π

params = {
    'encut': 600,
    'kpts': 0.08,
    'voskown': 1,
    'lnoncollinear': True,
    'lsorbit': True,
    'gga_compat': False,
    # ... other non-collinear parameters
}
```

**Output:**
- Structure: `setting{N}_{configuration}_standardized.vasp`
- Config: `{configuration}_mae_config.txt`
- Jobs: (21 × 11) = 231 calculations submitted

### Step 2: Analyze Grid Results

Finds easy axis (minimum energy) and hard axis (maximum energy).

```bash
python 2_analyze_results_refactored.py
```

**What it does:**
1. Loads all grid calculation results using `MAEResultsLoader`
2. Finds easy and hard magnetization axes using `MAEAnalyzer`
3. Calculates MAE = E(hard) - E(easy)
4. Converts to MJ/m³ units
5. Generates 3D surface plots using `MAEPlotter`

**Output:**
```
Easy axis: [0.000, 0.000, 1.000] (θ=0.0°, φ=0.0°)
Hard axis: [1.000, 0.000, 0.000] (θ=90.0°, φ=0.0°)
MAE = 2.456 MJ/m³ (0.0123 eV)
```

Files:
- `mae_analysis_{configuration}.txt` - Text report
- `mae_surface_{configuration}.png` - 3D surface plot

**Add to `input.py`:**
```python
# From analysis output
easy_axis = np.array([0.0, 0.0, 1.0])
hard_axis = np.array([1.0, 0.0, 0.0])
```

### Step 3: Submit MAE Curve

Calculates high-resolution energy curve along rotation path.

```bash
python 3_submit_mae_curve_refactored.py
```

**What it does:**
1. Uses easy/hard axes from Step 2
2. Generates rotation path using `MAEDirectionGenerator.generate_mae_curve()`
3. Submits (N_MAE+1) non-collinear calculations along the path

**Input (`input.py`):**
```python
easy_axis = np.array([0.0, 0.0, 1.0])  # From Step 2
hard_axis = np.array([1.0, 0.0, 0.0])  # From Step 2

N_MAE = 20  # Number of points along curve
```

**Output:**
- Jobs: 21 calculations along rotation path

### Step 4: Plot and Analyze MAE Curve

Plots MAE curve and calculates magnetic properties.

```bash
python 4_plot_mae_curve_refactored.py
```

**What it does:**
1. Loads curve results
2. Fits to E(α) = K₁sin²α + K₂sin⁴α + c
3. Extracts magnetization from reference calculation
4. Calculates derived properties:
   - Magnetization M₀
   - Maximum energy product (BH)ₘₐₓ
   - Anisotropy field μ₀Hₐ
   - Magnetic hardness parameter
5. Generates publication-quality plots

**Output:**
```
MAE RESULTS
======================================
MAE (total):     2.456789 MJ/m³
MAE (per atom):  0.123456 MJ/m³
MAE (eV):        0.012345 eV

Fitted parameters:
K1 = 2.345678 MJ/m³
K2 = 0.111111 MJ/m³

Magnetic Properties:
M0 = 45.6789 μB (Bohr magnetons)
M0 = 1.234 MA/m
(BH)max = 123.45 kJ/m³
μ₀Hₐ = 2.34 T
Hardness = 0.567
```

Files:
- `MAE_curve_U{U}_J{J}.png` - Plot with fitted curve
- `MAE_curve_{configuration}_U{U}_J{J}.txt` - Comprehensive report
- `.npy` files - Numerical data for further analysis

## 🎨 SOLID Principles in Action

### Single Responsibility Principle (SRP)

Each class has **one reason to change**:

| Class | Responsibility | Would Change If... |
|-------|---------------|-------------------|
| `MAEDirectionGenerator` | Generate directions | Direction generation algorithm changes |
| `MAEAnalyzer` | Analyze results | Analysis method changes |
| `MAEResultsLoader` | Load files | File format changes |
| `MAEPlotter` | Create plots | Visualization requirements change |

### Open/Closed Principle (OCP)

Easy to **extend without modifying existing code**:

```python
# Add new MAE analysis method WITHOUT modifying MAEAnalyzer
class SymmetryBasedMAEAnalyzer(MAEAnalyzer):
    @staticmethod
    def calculate_mae_from_symmetry(crystal_system, anisotropy_constants):
        # New method using symmetry
        pass
```

### Dependency Inversion Principle (DIP)

High-level workflow depends on **abstractions**, not concrete implementations:

```python
# CalculationService depends on WorkflowSubmitter interface
service = CalculationService(submitter)  # Works with ANY submitter

# Can switch between FireWorks and Manual seamlessly
submitter = FireworksWorkflowSubmitter(...)  # Or ManualJobSubmitter(...)
```

### Interface Segregation Principle (ISP)

Each interface is **focused and minimal**:

```python
class MAEResultsLoader:
    # Only loading methods - doesn't force analysis methods
    def load_theta_phi_results(...)
    def load_reference_energy(...)

class MAEAnalyzer:
    # Only analysis methods - doesn't force loading methods
    def find_easy_hard_axes(...)
    def calculate_mae(...)
```

## 📊 Comparison: Old vs. New

### Before: Monolithic Script

`MAE.py` (1086 lines) contained everything:
- ❌ Direction generation mixed with job submission
- ❌ Analysis mixed with plotting
- ❌ Hard to test individual components
- ❌ Difficult to reuse parts in other projects
- ❌ Complex interdependencies

### After: Modular Architecture

Split into focused components:
- ✅ Clear separation of concerns
- ✅ Each class easily testable
- ✅ Components reusable across projects
- ✅ Easy to extend with new features
- ✅ Simple, linear workflows

**Example - Testing MAE Analysis:**

```python
# Before: Had to run entire 1000+ line script
# After: Test just the analyzer
def test_mae_analyzer():
    directions = [MAEDirection(0, 0, (0,0,1)), ...]
    energies = [-100.5, -100.3, ...]
    
    easy, hard = MAEAnalyzer.find_easy_hard_axes(directions, energies)
    mae = MAEAnalyzer.calculate_mae(min(energies), max(energies))
    
    assert mae > 0
```

## 🧪 Testing

Each MAE component can be tested independently:

```python
# Test direction generation
def test_theta_phi_grid():
    dirs = MAEDirectionGenerator.generate_theta_phi_grid(10, 20)
    assert len(dirs) == 11 * 21  # (n+1) points
    assert all(np.isclose(np.linalg.norm(d.saxis), 1.0) for d in dirs)

# Test MAE calculation
def test_mae_calculation():
    mae = MAEAnalyzer.calculate_mae(-100.5, -100.3)
    assert mae == 0.2

# Test curve fitting
def test_mae_curve_fitting():
    angles = np.linspace(0, 2*np.pi, 21)
    energies = 2.0 * np.sin(angles)**2  # K1 = 2.0
    k1, k2, c = MAEAnalyzer.fit_mae_curve(angles, energies)
    assert np.isclose(k1, 2.0, atol=0.01)
```

## 🔧 Customization

### Change Grid Resolution

```python
# Fine grid (more accurate, slower)
Nph = 40
Nth = 20

# Coarse grid (faster, less accurate)
Nph = 10
Nth = 5
```

### Change Curve Resolution

```python
# High resolution curve
N_MAE = 40

# Low resolution curve
N_MAE = 10
```

### Use Different K-points for Curve

```python
# In Step 3, modify base_params
params['kpts'] = 0.05  # Finer k-point mesh for curve
```

## 🐛 Troubleshooting

### Issue: Easy/hard axes are identical

**Cause:** MAE is very small or grid is too coarse

**Solution:**
1. Increase grid resolution: `Nph = 40, Nth = 20`
2. Check if spin-orbit coupling is enabled: `lsorbit = True`
3. Verify non-collinear calculation: `lnoncollinear = True`

### Issue: Calculations not completing

**Cause:** Non-collinear calculations require more time/memory

**Solution:**
1. Increase wall time in job header
2. Use more nodes/cores
3. Reduce k-point density for initial grid

### Issue: Fitted curve doesn't match data

**Cause:** MAE may have higher-order anisotropy terms

**Solution:**
1. Increase curve resolution: `N_MAE = 40`
2. Check if easy/hard axes are correct
3. Consider symmetry-allowed terms

## 📚 Physics Background

### MAE Fundamentals

The magnetocrystalline anisotropy energy describes how the total energy varies with magnetization direction:

```
E(θ) = K₁ sin²θ + K₂ sin⁴θ + ...
```

Where:
- **K₁, K₂**: Anisotropy constants (MJ/m³)
- **θ**: Angle from easy axis
- **MAE**: Energy difference between hard and easy axes

### Non-Collinear DFT

VASP parameters for MAE calculations:

```python
lnoncollinear = True   # Enable non-collinear magnetism
lsorbit = True         # Enable spin-orbit coupling
saxis = [x, y, z]      # Magnetization direction
voskown = 1            # Vosko-Wilk-Nusair interpolation
gga_compat = False     # Don't use old PAW datasets
```

### Magnetic Properties

Derived quantities:

1. **Magnetization**: M₀ = √(Mₓ² + Mᵧ² + Mᵤ²)
2. **Energy Product**: (BH)ₘₐₓ = ¼ μ₀M₀²
3. **Anisotropy Field**: μ₀Hₐ = 2K₁/(M₀)
4. **Hardness**: √(K₁/(μ₀M₀²))

## 🌟 Best Practices

1. **Always run collinear calculations first** (Step 0: `2_coll`)
2. **Use converged k-points and ENCUT** (from `0_conv_tests`)
3. **Start with coarse grid** (Nph=10, Nth=5) to test
4. **Verify spin-orbit coupling** is enabled
5. **Check magnetic moments** are reasonable
6. **Plot surface** before running curve
7. **Save all output files** for reproducibility

## 🎯 Benefits Summary

| Aspect | Improvement |
|--------|------------|
| **Code Organization** | 5 focused classes vs. 1 monolithic script |
| **Testability** | Each component testable independently |
| **Maintainability** | Clear responsibilities, easy to modify |
| **Extensibility** | Add new analysis methods without touching core |
| **Reusability** | Components work across different projects |
| **Documentation** | Self-documenting class names and methods |
| **Error Handling** | Validation at each step |

## 📖 Further Reading

- [Core MAE Service](../core/services/mae_service.py) - Implementation details
- [Calculation Service](../core/services/calculation_service.py) - Workflow orchestration
- [REFACTORED_SCRIPTS_GUIDE.md](../docs/REFACTORED_SCRIPTS_GUIDE.md) - Complete guide
- [SOLID_QUICK_REFERENCE.md](../docs/SOLID_QUICK_REFERENCE.md) - Design principles

## 🤝 Contributing

To add new MAE analysis methods:

1. Add method to appropriate class in `mae_service.py`
2. Follow Single Responsibility Principle
3. Add unit tests
4. Update this documentation

Example:

```python
class MAEAnalyzer:
    @staticmethod
    def calculate_mae_temperature_dependence(
        k1: float, k2: float, temperature: float
    ) -> float:
        """Calculate MAE at given temperature."""
        # Implementation
        pass
```

---

**The refactored MAE workflow is production-ready and fully compatible with the existing Automag-1 infrastructure!** 🚀
