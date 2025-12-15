# MAE 2D Energy vs Theta Plot Feature

## Overview

A new 2D plotting feature has been added to the MAE analysis workflow that displays **Energy vs Theta** dependencies at different Phi values, with energy automatically converted to **MJ/m³** units.

## What Was Added

### 1. New Method in `MAEPlotter` Class

**File**: `core/services/mae_service.py`

**Method**: `plot_energy_vs_theta_at_phi()`

```python
@staticmethod
def plot_energy_vs_theta_at_phi(theta_grid: np.ndarray,
                                 phi_grid: np.ndarray,
                                 energy_grid: np.ndarray,
                                 volume: float,
                                 output_file: str = 'mae_energy_vs_theta.png'):
    """
    Plot 2D E vs Theta curves for different Phi values.
    
    Args:
        theta_grid: Theta values (Nph+1, Nth+1) in radians
        phi_grid: Phi values (Nph+1, Nth+1) in radians
        energy_grid: Energy values in eV (Nph+1, Nth+1)
        volume: Unit cell volume in Ų
        output_file: Output file name
    """
```

**Features**:
- Automatically converts energies from eV to MJ/m³
- Plots E vs θ curves for each φ value
- Uses different colors for each φ (viridis colormap)
- Displays energies relative to minimum (MAE convention)
- Handles legends intelligently (multi-column for many curves)
- Adds grid and proper axis labels

### 2. Integration into Analysis Script

**File**: `4_mae/2_analyze_results_refactored.py`

The script now generates **two plots** when analysis completes:

1. **3D Surface Plot** (existing): `mae_surface_{configuration}.png`
   - Shows energy landscape over entire θ-φ space
   - Energy in eV

2. **2D Line Plot** (new): `mae_theta_phi_{configuration}.png`
   - Shows E(θ) curves at different φ values
   - Energy in MJ/m³
   - Each φ value shown in different color

## Energy Unit Conversion

The conversion from eV to MJ/m³ is performed automatically:

```python
# Physical constants
eV = 1.602176634e-19  # Joules per eV
Ang = 1e-10  # meters per Angstrom

# Conversion formula
energy_rel_ev = energy_grid - np.min(energy_grid)  # Relative to minimum
energy_mj_m3 = (energy_rel_ev * eV / (volume * Ang**3)) * 1e-6
```

Where:
- `energy_grid`: Raw VASP energies in eV
- `volume`: Unit cell volume in Ų (from structure)
- Result: Energy density in MJ/m³

## Usage

Simply run the analysis script as before:

```bash
cd /path/to/automag-1/4_mae
python 2_analyze_results_refactored.py
```

**Output**:
```
MAE ANALYSIS FOR Fe12O18 - Configuration: fm1
======================================================================

Analyzing 231 calculations...
Found MAE directory: /path/to/CalcFold/.../mae_U0.0_J0.0_K20_EN830
Reference energy: -354.123456 eV
Loading grid calculation results...
Loaded 231 results out of 231 calculations

======================================================================
RESULTS
======================================================================
Easy axis direction:
  Theta: 0.00°
  Phi:   0.00°
  SAXIS: (0.000, 0.000, 1.000)
  Energy: -354.123456 eV

Hard axis direction:
  Theta: 90.00°
  Phi:   0.00°
  SAXIS: (1.000, 0.000, 0.000)
  Energy: -354.120123 eV

MAE:
  0.003333 eV
  3.333 meV
  125.456 MJ/m³
======================================================================

Surface plot saved to mae_surface_fm1.png
2D Energy vs Theta plot saved to mae_theta_phi_fm1.png
Results saved to: mae_results_fm1.txt

Plots generated:
  1. 3D surface plot: mae_surface_fm1.png
  2. 2D E(θ) at different φ: mae_theta_phi_fm1.png

Next step: Use the easy/hard axes above for MAE curve calculation
Add to input.py and run 3_submit_mae_curve_refactored.py
```

## Plot Interpretation

### 2D E vs θ Plot

The new plot shows:

- **X-axis**: Theta angle (0° to 180°)
- **Y-axis**: Energy in MJ/m³ (relative to minimum)
- **Multiple curves**: One for each φ value
- **Color coding**: From blue (φ=0°) to yellow (φ=360°) using viridis colormap

**Physical Interpretation**:

1. **Lowest curve**: Indicates the easy axis direction
2. **Highest curve**: Indicates the hard axis direction  
3. **Curve separation**: Shows anisotropy strength
4. **Curve shape**: Reveals symmetry of the magnetic system

For example:
- **Cubic symmetry**: All curves should overlap (isotropic)
- **Uniaxial symmetry**: Curves should show rotational symmetry around one axis
- **Biaxial/lower symmetry**: More complex curve patterns

## Example Output

For a typical MAE calculation with Nph=20, Nth=10 (231 points):

**mae_theta_phi_fm1.png** will show:
- 21 colored curves (one for each φ from 0° to 360°)
- θ from 0° to 180° on x-axis
- Energy in MJ/m³ on y-axis
- Legend showing φ values
- Grid for easy reading

## Advantages Over 3D Surface Plot

1. **Clearer trends**: Easier to see E(θ) dependence at specific φ
2. **Quantitative reading**: MJ/m³ units directly comparable to literature
3. **Better for publication**: Standard format in magnetism papers
4. **Easier analysis**: Can identify easy/hard axes by inspecting curves
5. **Symmetry analysis**: Curve patterns reveal magnetic symmetry

## Technical Details

### Grid Structure

The data is organized as:
```
energy_grid[i, j] where:
  i = phi index (0 to Nph)
  j = theta index (0 to Nth)
```

Each row `i` represents one φ value with varying θ.

### Color Mapping

Uses matplotlib's viridis colormap:
- Perceptually uniform
- Colorblind-friendly
- Print-friendly (converts well to grayscale)

### Legend Handling

- **≤10 curves**: Single column, font size 9
- **>10 curves**: Two columns, font size 8

This ensures readability even with many φ values.

## Customization

To customize the plot, edit the `plot_energy_vs_theta_at_phi()` method in `core/services/mae_service.py`:

```python
# Change colormap
colors = cm.plasma(np.linspace(0, 1, n_phi))  # Instead of viridis

# Change line style
ax.plot(theta_deg, energy_slice, '-', ...)  # No markers

# Change figure size
fig, ax = plt.subplots(figsize=(12, 8))  # Larger plot

# Change energy units (keep in eV instead of MJ/m³)
energy_to_plot = energy_rel_ev  # Don't convert
ax.set_ylabel('Energy (eV)', fontsize=14)
```

## Integration with Existing Workflow

This feature is **fully integrated** with the existing MAE analysis workflow:

1. **No breaking changes**: Existing functionality unchanged
2. **Automatic execution**: Plot generated if grid is complete
3. **Same input requirements**: Uses existing input.py parameters
4. **Consistent naming**: Follows existing file naming conventions
5. **SOLID principles**: Added as new method in `MAEPlotter` class

## Future Enhancements

Potential improvements:

1. **Interactive plots**: Use plotly for zoom/pan capabilities
2. **Export data**: Save curve data to CSV for external plotting
3. **Fit functions**: Overlay analytical fits (sin²θ, sin⁴θ)
4. **Multiple configurations**: Compare different magnetic states
5. **Polar plots**: Display in polar coordinates for symmetry visualization

## Troubleshooting

**Q: Plot not generated?**
- Check that all calculations completed (need complete θ-φ grid)
- Verify expected_points matches number of valid_directions

**Q: All curves overlap?**
- This is physical! Indicates isotropic or very small MAE
- Check MAE value - might be below numerical precision

**Q: Legend too crowded?**
- Reduce Nph in input.py (e.g., use Nph=10 instead of 20)
- Or customize legend in mae_service.py

**Q: Energy scale too small to see?**
- MAE might be very small (< 0.1 MJ/m³)
- This is valid for materials with weak anisotropy
- Consider logarithmic scale for very small values
