"""
Abstract interfaces for workflow submission.

Follows:
- Open/Closed Principle: Extensible without modification
- Dependency Inversion Principle: Depend on abstractions
- Liskov Substitution Principle: All implementations are substitutable
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional
from dataclasses import dataclass
from ase import Atoms


@dataclass
class CalculationConfig:
    """Configuration for a single calculation (immutable)."""
    calc_params: Dict
    magmoms: List[float]
    name: str
    encode: Optional[str] = None
    pert_step: Optional[str] = None
    pert_value: Optional[float] = None
    dummy_atom: Optional[str] = None
    atom_ucalc: Optional[str] = None


class WorkflowSubmitter(ABC):
    """
    Abstract interface for workflow submission.
    
    Implementations can use FireWorks, manual job submission, or any other system.
    This follows the Open/Closed Principle - new submission methods can be added
    without modifying existing code.
    """
    
    @abstractmethod
    def submit_single_calculation(self, atoms: Atoms, config: CalculationConfig) -> str:
        """
        Submit a single calculation.
        
        Args:
            atoms: ASE Atoms object representing the structure
            config: Configuration for the calculation
            
        Returns:
            Job identifier
        """
        pass
    
    @abstractmethod
    def submit_workflow(self, atoms: Atoms, configs: List[CalculationConfig], 
                       workflow_name: str) -> str:
        """
        Submit a workflow with multiple calculations.
        
        Args:
            atoms: ASE Atoms object representing the structure
            configs: List of calculation configurations
            workflow_name: Name for the workflow
            
        Returns:
            Workflow identifier
        """
        pass
    
    @abstractmethod
    def check_status(self, job_id: str) -> str:
        """
        Check the status of a submitted job/workflow.
        
        Args:
            job_id: Job or workflow identifier
            
        Returns:
            Status string (e.g., 'RUNNING', 'COMPLETED', 'FAILED')
        """
        pass


class CalculationInputWriter(ABC):
    """
    Abstract interface for writing calculation input files.
    
    This follows Interface Segregation Principle - focused on single responsibility.
    """
    
    @abstractmethod
    def write_input_files(self, atoms: Atoms, params: Dict, workdir: str) -> None:
        """
        Write input files for a calculation.
        
        Args:
            atoms: ASE Atoms object
            params: Calculation parameters
            workdir: Working directory for the calculation
        """
        pass
    
    @abstractmethod
    def set_magnetic_moments(self, atoms: Atoms, magmoms: List[float]) -> Atoms:
        """
        Set magnetic moments on atoms object.
        
        Args:
            atoms: ASE Atoms object
            magmoms: Magnetic moments to set
            
        Returns:
            Atoms object with magnetic moments set
        """
        pass


class ResultParser(ABC):
    """
    Abstract interface for parsing calculation results.
    
    Follows Interface Segregation Principle - clients only depend on what they need.
    """
    
    @abstractmethod
    def parse_energy(self, workdir: str) -> float:
        """Parse energy from calculation output."""
        pass
    
    @abstractmethod
    def parse_magnetic_moments(self, workdir: str) -> List[float]:
        """Parse final magnetic moments from calculation output."""
        pass
    
    @abstractmethod
    def check_convergence(self, workdir: str) -> bool:
        """Check if calculation converged."""
        pass
