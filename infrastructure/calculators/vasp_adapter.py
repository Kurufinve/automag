"""
VASP calculator adapter.

Refactored from utilities.py and SubmitManual.py - follows SRP and ISP.
"""

import os
import numpy as np
from pathlib import Path
from typing import Dict, List
from ase import Atoms
from ase.calculators.vasp import Vasp

from core.interfaces.workflow import CalculationInputWriter, ResultParser


class VaspInputWriter(CalculationInputWriter):
    """
    VASP-specific input file writer.
    
    Follows:
    - Single Responsibility Principle: Only writes VASP inputs
    - Interface Segregation Principle: Focused interface
    """
    
    def write_input_files(self, atoms: Atoms, params: Dict, workdir: str) -> None:
        """
        Write VASP input files (INCAR, POSCAR, KPOINTS, POTCAR).
        
        Args:
            atoms: ASE Atoms object
            params: VASP calculation parameters
            workdir: Working directory
        """
        os.makedirs(workdir, exist_ok=True)
        
        # Convert lists to numpy arrays
        params_copy = params.copy()
        for key, value in params_copy.items():
            if isinstance(value, list):
                params_copy[key] = np.asarray(value)
        
        # Create VASP calculator and write input
        calc = Vasp(atoms=atoms, directory=workdir, **params_copy)
        calc.write_input(atoms)
    
    def set_magnetic_moments(self, atoms: Atoms, magmoms: List[float]) -> Atoms:
        """
        Set initial magnetic moments on atoms.
        
        Args:
            atoms: ASE Atoms object
            magmoms: List of magnetic moments
            
        Returns:
            Atoms object with magnetic moments set
        """
        atoms_copy = atoms.copy()
        atoms_copy.set_initial_magnetic_moments(np.array(magmoms))
        return atoms_copy


class VaspResultParser(ResultParser):
    """
    VASP-specific result parser.
    
    Follows:
    - Single Responsibility Principle: Only parses VASP outputs
    - Interface Segregation Principle: Clients only use what they need
    """
    
    def parse_energy(self, workdir: str) -> float:
        """
        Parse energy from VASP OUTCAR.
        
        Args:
            workdir: Working directory containing OUTCAR
            
        Returns:
            Final energy in eV
        """
        from ase.io import read
        atoms = read(os.path.join(workdir, 'OUTCAR'))
        return atoms.get_potential_energy()
    
    def parse_magnetic_moments(self, workdir: str) -> List[float]:
        """
        Parse final magnetic moments from VASP output.
        
        Args:
            workdir: Working directory containing OUTCAR
            
        Returns:
            List of final magnetic moments
        """
        from ase.io import read
        atoms = read(os.path.join(workdir, 'OUTCAR'))
        
        try:
            magmoms = atoms.get_magnetic_moments()
            return magmoms.tolist()
        except:
            # Return zeros if magnetic moments not available
            return [0.0] * len(atoms)
    
    def check_convergence(self, workdir: str) -> bool:
        """
        Check if VASP calculation converged.
        
        Args:
            workdir: Working directory
            
        Returns:
            True if converged
        """
        calc = Vasp(directory=workdir)
        return calc.read_convergence()
    
    def parse_enthalpy(self, workdir: str) -> float:
        """
        Parse enthalpy from VASP OUTCAR.
        
        Args:
            workdir: Working directory
            
        Returns:
            Enthalpy in eV
        """
        outcar_path = os.path.join(workdir, 'OUTCAR')
        enthalpy = None
        
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'enthalpy' in line:
                    enthalpy_line = line
        
        if enthalpy_line:
            enthalpy = float(enthalpy_line.split()[4])
        else:
            # Fallback to energy
            enthalpy = self.parse_energy(workdir)
        
        return enthalpy
    
    def parse_charges(self, workdir: str, atom_index: int) -> Dict[str, float]:
        """
        Parse charges for specific atom (for linear response).
        
        Args:
            workdir: Working directory
            atom_index: Index of atom to get charges for
            
        Returns:
            Dictionary with charge information
        """
        from pymatgen.io.vasp.outputs import Outcar
        
        outcar = Outcar(os.path.join(workdir, 'OUTCAR'))
        
        # Get d or f charges depending on orbital
        if 'f' in outcar.charge[atom_index]:
            charge = outcar.charge[atom_index]['f']
        else:
            charge = outcar.charge[atom_index]['d']
        
        return {'charge': charge}
