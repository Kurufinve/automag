"""
MAE (Magnetocrystalline Anisotropy Energy) calculation services.

Refactored from 4_mae/MAE.py following SOLID principles.

Follows:
- Single Responsibility Principle: Separate classes for different MAE tasks
- Dependency Inversion Principle: Depends on abstractions
- Open/Closed Principle: Extensible for new MAE analysis methods
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from dataclasses import dataclass

from pymatgen.core.structure import Structure
from ase import Atoms


@dataclass
class MAEDirection:
    """Represents a magnetization direction for MAE calculation."""
    theta: float  # Polar angle in radians
    phi: float    # Azimuthal angle in radians
    saxis: Tuple[float, float, float]  # Unit vector
    angle: float = 0.0  # Rotation angle from easy axis (for MAE curve)
    
    @classmethod
    def from_angles(cls, theta: float, phi: float) -> 'MAEDirection':
        """Create MAE direction from spherical coordinates."""
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        return cls(theta=theta, phi=phi, saxis=(x, y, z))
    
    @property
    def name(self) -> str:
        """Get name for this direction."""
        theta_deg = np.round((self.theta / np.pi) * 180, 2)
        phi_deg = np.round((self.phi / np.pi) * 180, 2)
        return f'PhTh_{phi_deg}_{theta_deg}'


class MAEDirectionGenerator:
    """
    Generator for MAE calculation directions.
    
    Follows Single Responsibility Principle - only generates directions.
    """
    
    @staticmethod
    def generate_theta_phi_grid(n_theta: int, n_phi: int) -> List[MAEDirection]:
        """
        Generate grid of directions for theta-phi exploration.
        
        Args:
            n_theta: Number of theta points (0 to pi)
            n_phi: Number of phi points (0 to 2*pi)
            
        Returns:
            List of MAE directions
        """
        phi_values = np.linspace(0.0001, 2 * np.pi - 0.0001, n_phi + 1)
        theta_values = np.linspace(0.0001, np.pi - 0.0001, n_theta + 1)
        
        directions = []
        for phi in phi_values:
            for theta in theta_values:
                directions.append(MAEDirection.from_angles(theta, phi))
        
        return directions
    
    @staticmethod
    def generate_mae_curve(n_points: int, 
                          easy_axis: np.ndarray,
                          hard_axis: np.ndarray) -> List[MAEDirection]:
        """
        Generate points along MAE curve from easy to hard axis.
        
        Args:
            n_points: Number of points
            easy_axis: Easy magnetization axis (unit vector)
            hard_axis: Hard magnetization axis (unit vector)
            
        Returns:
            List of MAE directions with rotation angles
        """
        directions = []
        
        # Ensure normalized
        easy_axis = easy_axis / np.linalg.norm(easy_axis)
        hard_axis = hard_axis / np.linalg.norm(hard_axis)
        
        # Generate rotation axis perpendicular to both
        rotation_axis = np.cross(easy_axis, hard_axis)
        if np.linalg.norm(rotation_axis) < 1e-10:
            # Axes are parallel, use arbitrary perpendicular axis
            if abs(easy_axis[2]) < 0.9:
                rotation_axis = np.cross(easy_axis, np.array([0, 0, 1]))
            else:
                rotation_axis = np.cross(easy_axis, np.array([1, 0, 0]))
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        
        # Calculate angle between easy and hard axes
        angle_max = np.arccos(np.clip(np.dot(easy_axis, hard_axis), -1, 1))
        
        for i, alpha in enumerate(np.linspace(0, angle_max, n_points + 1)):
            # Rotate from easy axis using Rodrigues' rotation formula
            saxis = (easy_axis * np.cos(alpha) + 
                    rotation_axis * np.sin(alpha) * np.dot(rotation_axis, easy_axis) +
                    np.cross(rotation_axis, easy_axis) * np.sin(alpha))
            saxis = saxis / np.linalg.norm(saxis)
            
            # Convert to spherical coordinates
            theta = np.arccos(np.clip(saxis[2], -1, 1))
            phi = np.arctan2(saxis[1], saxis[0])
            if phi < 0:
                phi += 2 * np.pi
            
            # Store angle in degrees for directory naming
            angle_deg = alpha * 180.0 / np.pi
            
            directions.append(MAEDirection(
                theta=theta, 
                phi=phi, 
                saxis=tuple(saxis),
                angle=angle_deg
            ))
        
        return directions


class MAEAnalyzer:
    """
    Analyzer for MAE calculation results.
    
    Follows Single Responsibility Principle - only analyzes MAE data.
    """
    
    @staticmethod
    def find_easy_hard_axes(directions: List[MAEDirection],
                           energies: List[float]) -> Tuple[MAEDirection, MAEDirection]:
        """
        Find easy (minimum energy) and hard (maximum energy) axes.
        
        Args:
            directions: List of MAE directions
            energies: Corresponding energies
            
        Returns:
            Tuple of (easy_direction, hard_direction)
        """
        min_idx = np.argmin(energies)
        max_idx = np.argmax(energies)
        
        return directions[min_idx], directions[max_idx]
    
    @staticmethod
    def calculate_mae(min_energy: float, max_energy: float) -> float:
        """
        Calculate MAE as energy difference.
        
        Args:
            min_energy: Minimum energy
            max_energy: Maximum energy
            
        Returns:
            MAE in eV
        """
        return max_energy - min_energy
    
    @staticmethod
    def calculate_mae_per_volume(mae: float, volume: float) -> float:
        """
        Calculate MAE per unit volume in MJ/m³.
        
        Args:
            mae: MAE in eV
            volume: Unit cell volume in Ų
            
        Returns:
            MAE in MJ/m³
        """
        eV = 1.602176634e-19  # J
        Ang = 1e-10  # m
        
        # Convert eV/cell to J/m³
        mae_si = mae * eV / (volume * Ang**3)
        
        # Convert to MJ/m³
        return mae_si * 1e-6
    
    @staticmethod
    def fit_mae_curve(angles: np.ndarray, 
                     energies: np.ndarray) -> Tuple[float, float, float]:
        """
        Fit MAE curve to sin²θ and sin⁴θ terms.
        
        E(θ) = K₁sin²θ + K₂sin⁴θ + c
        
        Args:
            angles: Rotation angles in radians
            energies: Corresponding energies
            
        Returns:
            Tuple of (K1, K2, constant)
        """
        from scipy.optimize import curve_fit
        
        def mae_function(theta, k1, k2, c):
            return k1 * (np.sin(theta)**2) + k2 * (np.sin(theta)**4) + c
        
        try:
            popt, _ = curve_fit(mae_function, angles, energies)
            return tuple(popt)
        except:
            # If fitting fails, return zeros
            return (0.0, 0.0, np.mean(energies))


class MAEResultsLoader:
    """
    Loader for MAE calculation results.
    
    Follows Single Responsibility Principle - only loads results.
    """
    
    def __init__(self, base_path: Path):
        """
        Initialize results loader.
        
        Args:
            base_path: Base path for MAE calculations
        """
        self.base_path = Path(base_path)
    
    def load_theta_phi_results(self, 
                               directions: List[MAEDirection]) -> Dict[str, float]:
        """
        Load energies for theta-phi grid calculations.
        
        Args:
            directions: List of MAE directions
            
        Returns:
            Dictionary mapping direction name to energy
        """
        from pymatgen.io.vasp.outputs import Oszicar
        
        results = {}
        
        for direction in directions:
            folder = self.base_path / direction.name
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
        
        return results
    
    def load_reference_energy(self, reference_dir: str = 'z') -> Optional[float]:
        """
        Load reference energy from self-consistent reference calculation.
        
        Args:
            reference_dir: Directory name for reference calculation
            
        Returns:
            Reference energy or None if not found
        """
        from pymatgen.io.vasp.outputs import Oszicar
        
        # Reference calculation OSZICAR is directly in the reference folder
        ref_path = self.base_path / reference_dir / 'OSZICAR'
        
        if ref_path.exists():
            try:
                oszicar = Oszicar(str(ref_path))
                return float(oszicar.all_energies[-1][-2])
            except Exception as e:
                print(f"Error reading reference OSZICAR: {e}")
                pass
        else:
            print(f"Reference OSZICAR not found at: {ref_path}")
        
        return None


class MAEPlotter:
    """
    Plotter for MAE results.
    
    Follows Single Responsibility Principle - only handles plotting.
    """
    
    @staticmethod
    def plot_theta_phi_surface(theta_grid: np.ndarray,
                               phi_grid: np.ndarray,
                               energy_grid: np.ndarray,
                               output_file: str = 'mae_surface.png'):
        """
        Plot 3D surface of MAE vs theta and phi.
        
        Args:
            theta_grid: Theta values
            phi_grid: Phi values
            energy_grid: Energy values
            output_file: Output file name
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Convert to degrees for plotting
        theta_deg = theta_grid * 180 / np.pi
        phi_deg = phi_grid * 180 / np.pi
        
        surf = ax.plot_surface(phi_deg, theta_deg, energy_grid,
                              cmap='viridis', alpha=0.8)
        
        ax.set_xlabel('Phi (degrees)', fontsize=12)
        ax.set_ylabel('Theta (degrees)', fontsize=12)
        ax.set_zlabel('Energy (eV)', fontsize=12)
        ax.set_title('MAE Energy Surface', fontsize=14)
        
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Surface plot saved to {output_file}")
    
    @staticmethod
    def plot_mae_curve(angles: np.ndarray,
                      energies: np.ndarray,
                      fitted_energies: Optional[np.ndarray] = None,
                      output_file: str = 'mae_curve.png'):
        """
        Plot MAE curve along rotation path.
        
        Args:
            angles: Rotation angles in radians
            energies: Energies at each angle
            fitted_energies: Fitted energies (optional)
            output_file: Output file name
        """
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        # Convert to degrees
        angles_deg = angles * 180 / np.pi
        
        # Relative to minimum
        energies_rel = (energies - np.min(energies)) * 1000  # meV
        
        plt.plot(angles_deg, energies_rel, 'o-', label='DFT', linewidth=2, markersize=8)
        
        if fitted_energies is not None:
            fitted_rel = (fitted_energies - np.min(energies)) * 1000
            plt.plot(angles_deg, fitted_rel, '--', label='Fitted', linewidth=2)
        
        plt.xlabel('Rotation Angle (degrees)', fontsize=14)
        plt.ylabel('Relative Energy (meV)', fontsize=14)
        plt.title('MAE Curve', fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"MAE curve saved to {output_file}")
    
    @staticmethod
    def plot_energy_vs_theta_at_phi(theta_grid: np.ndarray,
                                     phi_grid: np.ndarray,
                                     energy_grid: np.ndarray,
                                     volume: float,
                                     output_file: str = 'mae_energy_vs_theta.png'):
        """
        Plot 2D E vs Theta curves for different Phi values.
        
        Args:
            theta_grid: Theta values (Nph+1, Nth+1)
            phi_grid: Phi values (Nph+1, Nth+1)
            energy_grid: Energy values in eV (Nph+1, Nth+1)
            volume: Unit cell volume in Ų
            output_file: Output file name
        """
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        
        # Convert energy from eV to MJ/m³
        eV = 1.602176634e-19  # J
        Ang = 1e-10  # m
        
        # Subtract minimum energy to get relative energies
        energy_rel_ev = energy_grid - np.min(energy_grid)
        
        # Convert to MJ/m³
        energy_mj_m3 = (energy_rel_ev * eV / (volume * Ang**3)) * 1e-6
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 7))
        
        # Get number of phi values
        n_phi = phi_grid.shape[0]
        
        # Create colormap
        colors = cm.viridis(np.linspace(0, 1, n_phi))
        
        # Plot each phi slice
        for i in range(n_phi):
            theta_deg = theta_grid[i, :] * 180 / np.pi
            phi_deg = phi_grid[i, 0] * 180 / np.pi  # Phi is constant for each row
            energy_slice = energy_mj_m3[i, :]
            
            # Skip if all zeros (no data)
            if not np.allclose(energy_slice, 0, atol=1e-10):
                ax.plot(theta_deg, energy_slice, '-o', 
                       color=colors[i], 
                       label=f'φ = {phi_deg:.1f}°',
                       linewidth=1.5, 
                       markersize=4,
                       alpha=0.8)
        
        ax.set_xlabel('Theta (degrees)', fontsize=14)
        ax.set_ylabel('Energy (MJ/m³)', fontsize=14)
        ax.set_title('MAE: Energy vs Theta at Different Phi Values', fontsize=16)
        ax.grid(True, alpha=0.3)
        
        # Add legend with smaller font and multiple columns if needed
        if n_phi <= 10:
            ax.legend(fontsize=9, loc='best')
        else:
            # For many phi values, use smaller font and multiple columns
            ax.legend(fontsize=8, ncol=2, loc='best')
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"2D Energy vs Theta plot saved to {output_file}")
