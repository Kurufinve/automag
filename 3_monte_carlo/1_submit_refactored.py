"""
Refactored coupling constants calculation script for Monte Carlo simulation.

This demonstrates the SOLID-compliant architecture for analyzing collinear results
and computing Heisenberg model coupling constants.

The old code is preserved in 1_coupling_constants.py for reference.

Follows:
- Single Responsibility Principle: Separate classes for different tasks
- Dependency Inversion Principle: Depend on result parser interface
- Open/Closed Principle: Easy to extend with new analysis methods
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from pymatgen.core.structure import Structure

# Import refactored components
from infrastructure.calculators.vasp_adapter import VaspResultParser

# Default values
calculator = 'vasp'
struct_suffix = ''

# Get current directory
cwd = os.getcwd()

# Load input configuration
try:
    input_file = sys.argv[1]
    print(f'Using the {input_file} file as input')
    exec(f'from {input_file.split(".")[0]} import *')
except IndexError:
    print(f'Using the input.py file from folder: {cwd}')
    from input import *


class CouplingConstantsCalculator:
    """
    Calculator for Heisenberg model coupling constants.
    
    Follows Single Responsibility Principle - only handles coupling constant calculation.
    """
    
    def __init__(self, 
                 structure: Structure,
                 states: List[List[int]],
                 energies: np.ndarray,
                 cutoff_radius: float):
        """
        Initialize calculator.
        
        Args:
            structure: Magnetic structure (non-magnetic atoms removed)
            states: List of magnetic configurations
            energies: DFT energies (total, not per atom)
            cutoff_radius: Maximum distance for neighbor interactions
        """
        self.structure = structure
        self.states = states
        self.energies = energies
        self.cutoff_radius = cutoff_radius
        
        # Get neighbor list
        self.center_indices, self.point_indices, self.offset_vectors, self.distances = \
            structure.get_neighbor_list(cutoff_radius)
        
        # Get unique distances
        self.unique_distances, self.counts = np.unique(
            np.around(self.distances, 3), return_counts=True
        )
    
    def calculate(self, control_group_size: float = 0.2) -> Dict:
        """
        Calculate coupling constants using linear fit.
        
        Args:
            control_group_size: Fraction of data to use for validation (0-1)
            
        Returns:
            Dictionary with results including coupling constants and PCC
        """
        # Split into fit and control groups
        fit_size = 1 - control_group_size
        split_index = int(round(len(self.states) * fit_size))
        
        states_fit = np.array(self.states[:split_index])
        states_control = np.array(self.states[split_index:])
        energies_fit = self.energies[:split_index]
        energies_control = self.energies[split_index:]
        
        # Build system matrix
        A = self._build_system_matrix(states_fit)
        
        # Solve for coupling constants
        solution = np.linalg.lstsq(A, energies_fit, rcond=None)
        
        # Calculate predictions for control group
        B = self._build_system_matrix(states_control)
        predictions = B @ solution[0]
        
        # Calculate Pearson correlation coefficient
        pcc = np.corrcoef(predictions, energies_control)[0, 1]
        
        # Check if system is well-determined
        is_valid = np.linalg.matrix_rank(A) == len(self.unique_distances) + 1
        
        if is_valid:
            # Convert to SI units (J)
            coupling_constants = solution[0][1:] * 1.60218e-19
        else:
            coupling_constants = None
        
        return {
            'is_valid': is_valid,
            'coupling_constants': coupling_constants,
            'distances': self.unique_distances,
            'counts': self.counts,
            'pcc': pcc,
            'predictions': predictions,
            'energies_control': energies_control,
            'matrix_rank': np.linalg.matrix_rank(A),
            'num_unknowns': len(self.unique_distances) + 1
        }
    
    def _build_system_matrix(self, configurations: np.ndarray) -> np.ndarray:
        """Build system matrix for linear fit."""
        matrix = []
        
        for config in configurations:
            equation = [1]  # Constant term
            
            # For each unique distance, count interactions
            for distance in self.unique_distances:
                count = 0
                for atom1, atom2, d in zip(
                    self.center_indices, self.point_indices, self.distances
                ):
                    if np.isclose(d, distance, atol=0.02):
                        count += config[atom1] * config[atom2]
                
                equation.append(-count // 2)
            
            matrix.append(equation)
        
        return np.array(matrix)


class ResultsLoader:
    """
    Loader for collinear calculation results.
    
    Follows Single Responsibility Principle - only loads results.
    """
    
    def __init__(self, 
                 results_path: Path,
                 formula: str,
                 calculator: str = 'vasp'):
        """
        Initialize results loader.
        
        Args:
            results_path: Path to results directory
            formula: Chemical formula
            calculator: Calculator name
        """
        self.results_path = results_path
        self.formula = formula
        self.calculator = calculator
    
    def load(self) -> Tuple[Optional[Structure], Optional[List], Optional[List]]:
        """
        Load structure, states, and energies from previous calculations.
        
        Returns:
            Tuple of (structure, states, energies) or None if not found
        """
        structure = None
        states = None
        energies = None
        
        # Search for result files
        for item in os.listdir(self.results_path):
            file_path = self.results_path / item
            
            if not file_path.is_file():
                continue
            
            # Structure file
            if (item.startswith(f'{self.formula}_{self.calculator}_setting') and 
                item.endswith('.vasp')):
                structure = Structure.from_file(str(file_path))
            
            # States file
            if (item.startswith(f'{self.formula}_{self.calculator}_states') and 
                item.endswith('.txt')):
                with open(file_path, 'r') as f:
                    states = json.load(f)
            
            # Energies file
            if (item.startswith(f'{self.formula}_{self.calculator}_energies') and 
                item.endswith('.txt')):
                with open(file_path, 'r') as f:
                    energies = json.load(f)
        
        return structure, states, energies


class ResultsPlotter:
    """
    Plotter for coupling constants results.
    
    Follows Single Responsibility Principle - only handles plotting.
    """
    
    @staticmethod
    def plot_model_fit(predictions: np.ndarray,
                      actual: np.ndarray,
                      pcc: float,
                      output_file: str = 'model.png'):
        """
        Plot Heisenberg model predictions vs DFT energies.
        
        Args:
            predictions: Model predictions
            actual: Actual DFT energies
            pcc: Pearson correlation coefficient
            output_file: Output file name
        """
        plt.rcParams.update({'font.size': 13})
        plt.figure(figsize=(8, 6))
        plt.locator_params(axis='x', nbins=5)
        plt.xlabel('Heisenberg model energy (eV)', size=15)
        plt.ylabel('DFT energy (eV)', size=15)
        plt.grid(True, alpha=0.3)
        
        plt.scatter(predictions, actual, label=f'PCC: {pcc:.2f}', alpha=0.6)
        
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        
        plt.legend(prop={'size': 15})
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plot saved to {output_file}")


def main():
    """Main execution function following SOLID principles."""
    
    # Get paths
    path_to_automag = os.environ.get('AUTOMAG_PATH')
    path_to_poscar = os.path.join(path_to_automag, 'geometries', poscar_file)
    input_structure = Structure.from_file(path_to_poscar)
    
    formula = input_structure.formula.replace(' ', '')
    
    # Find results directory
    results_path = _find_results_directory(path_to_automag, formula, struct_suffix, calculator)
    print(f'Results path: {results_path}')
    
    # Load results using dedicated loader (SRP)
    loader = ResultsLoader(results_path, formula, calculator)
    structure, states, energies = loader.load()
    
    # Validate loaded data
    if structure is None:
        raise IOError(f'No setting file found in {results_path} folder.')
    if states is None:
        raise IOError(f'No states file found in {results_path} folder.')
    if energies is None:
        raise IOError(f'No energies file found in {results_path} folder.')
    
    print(f"Loaded {len(states)} configurations")
    
    # Mark magnetic atoms
    for element in structure.composition.elements:
        if 'magnetic_atoms' not in globals():
            element.is_magnetic = element.is_transition_metal
        else:
            element.is_magnetic = (element.name in magnetic_atoms)
    
    # Remove non-magnetic atoms from structure
    non_magnetic_atoms = [
        element.symbol 
        for element in structure.composition.elements 
        if not element.is_magnetic
    ]
    structure.remove_species(non_magnetic_atoms)
    
    # Convert energies from eV/atom to total energy
    energies_total = structure.num_sites * np.array(energies)
    
    # Calculate coupling constants (SRP - dedicated calculator class)
    print(f"\nCalculating coupling constants with cutoff={cutoff_radius} Å...")
    calculator = CouplingConstantsCalculator(
        structure=structure,
        states=states,
        energies=energies_total,
        cutoff_radius=cutoff_radius
    )
    
    results = calculator.calculate(control_group_size=control_group_size)
    
    # Display results
    if results['is_valid']:
        print(f"\n{'=' * 70}")
        print("COUPLING CONSTANTS CALCULATION RESULTS")
        print(f"{'=' * 70}")
        print(f"Distances between neighbors: {results['distances'].tolist()}")
        print(f"Neighbor counts: {results['counts'].tolist()}")
        print(f"Coupling constants (J): {np.array2string(results['coupling_constants'], precision=8, separator=', ')}")
        print(f"Pearson Correlation Coefficient: {results['pcc']:.4f}")
        print(f"{'=' * 70}\n")
        
        # Optionally append to input file
        if append_coupling_constants:
            print('Appending coupling constants to input.py...')
            with open('input.py', 'a') as f:
                f.write('\n# LINE ADDED BY 1_submit_refactored.py')
                f.write(f'\ndistances_between_neighbors = {results["distances"].tolist()}\n')
                f.write('\n# LINE ADDED BY 1_submit_refactored.py')
                f.write(f'\ncoupling_constants = {np.array2string(results["coupling_constants"], precision=8, separator=", ")}\n')
        
        # Plot results (SRP - dedicated plotter class)
        plotter = ResultsPlotter()
        plotter.plot_model_fit(
            predictions=results['predictions'],
            actual=results['energies_control'],
            pcc=results['pcc']
        )
        
    else:
        print(f"\nERROR: SYSTEM OF {results['matrix_rank']} INDEPENDENT EQUATION(S) "
              f"IN {results['num_unknowns']} UNKNOWNS!")
        print("Try adjusting cutoff_radius or using more configurations.")


def _find_results_directory(automag_path: str, 
                            formula: str,
                            suffix: str,
                            calc: str) -> Path:
    """Find the directory containing collinear calculation results."""
    base_path = Path(automag_path) / '2_coll'
    
    # Try different possible paths
    possible_paths = [
        base_path / f'{formula}{suffix}' / f'{formula}_{calc}',
        base_path / f'{formula}_{calc}',
        base_path
    ]
    
    for path in possible_paths:
        if path.exists() and path.is_dir():
            return path
    
    # Default to base path
    return base_path


if __name__ == '__main__':
    main()
