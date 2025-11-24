"""
Service layer for orchestrating calculations.

Follows:
- Single Responsibility Principle: Each service has one purpose
- Dependency Inversion Principle: Depends on abstractions
"""

from typing import List, Dict, Optional, Tuple
from ase import Atoms
import numpy as np

from core.interfaces.workflow import WorkflowSubmitter, CalculationConfig
from core.domain.calculation import CalculationParameters, CalculationMode


class CalculationService:
    """
    Service for orchestrating calculation submissions.
    
    This service coordinates between different components following
    the Dependency Inversion Principle.
    """
    
    def __init__(self, submitter: WorkflowSubmitter):
        """
        Initialize calculation service.
        
        Args:
            submitter: Workflow submitter implementation (FireWorks or Manual)
        """
        self._submitter = submitter
    
    def submit_convergence_test(self,
                                atoms: Atoms,
                                base_params: CalculationParameters,
                                magmoms: List[float],
                                test_parameter: str,
                                test_values: List[float],
                                mode_name: str) -> List[str]:
        """
        Submit convergence test calculations.
        
        Args:
            atoms: Structure to test
            base_params: Base calculation parameters
            magmoms: Magnetic moments
            test_parameter: Parameter to test ('encut', 'sigma', 'kpts')
            test_values: Values to test
            mode_name: Mode name for output
            
        Returns:
            List of job IDs
        """
        job_ids = []
        
        for value in test_values:
            # Create parameters for this test value
            params_dict = base_params.to_dict()
            params_dict[test_parameter] = value
            
            # Create configuration
            config = CalculationConfig(
                calc_params=params_dict,
                magmoms=magmoms,
                name=f"{test_parameter}{value}"
            )
            
            # Submit
            job_id = self._submitter.submit_single_calculation(atoms, config)
            job_ids.append(job_id)
        
        return job_ids
    
    def submit_singlepoint_calculation(self,
                                      atoms: Atoms,
                                      params: CalculationParameters,
                                      magmoms: List[float],
                                      name: str,
                                      with_recalc: bool = True) -> str:
        """
        Submit single-point calculation (optionally with recalc).
        
        Args:
            atoms: Structure
            params: Calculation parameters
            magmoms: Magnetic moments
            name: Calculation name
            with_recalc: Whether to include recalculation step
            
        Returns:
            Job ID
        """
        if not with_recalc or name == 'nm':
            # Single calculation only
            config = CalculationConfig(
                calc_params=params.to_dict(),
                magmoms=magmoms,
                name=name
            )
            return self._submitter.submit_single_calculation(atoms, config)
        else:
            # Workflow with singlepoint + recalc
            configs = [
                CalculationConfig(
                    calc_params=params.to_dict(),
                    magmoms=magmoms,
                    name='singlepoint'
                ),
                CalculationConfig(
                    calc_params=params.to_dict(),
                    magmoms='previous',  # Use previous magmoms
                    name='recalc'
                )
            ]
            return self._submitter.submit_workflow(atoms, configs, name)
    
    def submit_perturbation_workflow(self,
                                    atoms: Atoms,
                                    base_params: CalculationParameters,
                                    magmoms: List[float],
                                    perturbations: List[float],
                                    dummy_atom: str,
                                    dummy_position: int,
                                    atom_ucalc: str,
                                    workflow_name: str) -> str:
        """
        Submit perturbation workflow for linear response U calculation.
        
        Args:
            atoms: Structure (with dummy atom substituted)
            base_params: Base calculation parameters
            magmoms: Magnetic moments
            perturbations: Perturbation values to apply
            dummy_atom: Dummy atom symbol
            dummy_position: Position of dummy atom
            atom_ucalc: Original atom type
            workflow_name: Name for workflow
            
        Returns:
            Workflow ID
        """
        # First: bare calculation
        configs = [
            CalculationConfig(
                calc_params=base_params.to_dict(),
                magmoms=magmoms,
                name='singlepoint'
            )
        ]
        
        # Add NSC and SC calculations for each perturbation
        for pert_value in perturbations:
            # Non-self-consistent
            configs.append(CalculationConfig(
                calc_params=base_params.to_dict(),
                magmoms=magmoms,
                name=f'nsc_{pert_value}',
                pert_step='NSC',
                pert_value=pert_value,
                dummy_atom=dummy_atom,
                atom_ucalc=atom_ucalc
            ))
            
            # Self-consistent
            configs.append(CalculationConfig(
                calc_params=base_params.to_dict(),
                magmoms=magmoms,
                name=f'sc_{pert_value}',
                pert_step='SC',
                pert_value=pert_value,
                dummy_atom=dummy_atom,
                atom_ucalc=atom_ucalc
            ))
        
        return self._submitter.submit_workflow(atoms, configs, workflow_name)
    
    def submit_mae_theta_phi_grid(self,
                                   atoms: Atoms,
                                   base_params: CalculationParameters,
                                   magmoms: List[Tuple[float, float, float]],
                                   n_theta: int,
                                   n_phi: int,
                                   reference_dir: str = 'z',
                                   workflow_name: str = 'mae_grid') -> List[str]:
        """
        Submit MAE calculations on theta-phi grid.
        
        Args:
            atoms: Structure
            base_params: Base calculation parameters (non-collinear)
            magmoms: Non-collinear magnetic moments (x,y,z) tuples
            n_theta: Number of theta points
            n_phi: Number of phi points
            reference_dir: Directory for reference collinear calculation
            workflow_name: Name for workflow
            
        Returns:
            List of job IDs
        """
        from core.services.mae_service import MAEDirectionGenerator
        
        # Generate grid directions
        directions = MAEDirectionGenerator.generate_theta_phi_grid(n_theta, n_phi)
        
        job_ids = []
        
        # Submit calculation for each direction
        for direction in directions:
            # Update SAXIS in parameters
            params_dict = base_params.to_dict()
            params_dict['saxis'] = list(direction.saxis)
            params_dict['icharg'] = 11  # Read charge density
            params_dict['istart'] = 1   # Read WAVECAR
            params_dict['lcharg'] = False
            params_dict['lwave'] = False
            
            config = CalculationConfig(
                calc_params=params_dict,
                magmoms=list(magmoms),
                name=direction.name
            )
            
            job_id = self._submitter.submit_single_calculation(atoms, config)
            job_ids.append(job_id)
        
        return job_ids
    
    def submit_mae_curve(self,
                        atoms: Atoms,
                        base_params: CalculationParameters,
                        magmoms: List[Tuple[float, float, float]],
                        n_points: int,
                        easy_axis: np.ndarray,
                        hard_axis: np.ndarray,
                        workflow_name: str = 'mae_curve') -> List[str]:
        """
        Submit MAE curve calculations along rotation path.
        
        Args:
            atoms: Structure
            base_params: Base calculation parameters (non-collinear)
            magmoms: Non-collinear magnetic moments
            n_points: Number of points along curve
            easy_axis: Easy magnetization axis
            hard_axis: Hard magnetization axis
            workflow_name: Name for workflow
            
        Returns:
            List of job IDs
        """
        from core.services.mae_service import MAEDirectionGenerator
        
        # Generate curve directions
        directions = MAEDirectionGenerator.generate_mae_curve(
            n_points, easy_axis, hard_axis
        )
        
        job_ids = []
        
        # Submit calculation for each direction
        for direction in directions:
            params_dict = base_params.to_dict()
            params_dict['saxis'] = list(direction.saxis)
            params_dict['icharg'] = 11
            params_dict['istart'] = 1
            params_dict['lcharg'] = False
            params_dict['lwave'] = False
            
            config = CalculationConfig(
                calc_params=params_dict,
                magmoms=list(magmoms),
                name=f'mae_curve_{direction.name}'
            )
            
            job_id = self._submitter.submit_single_calculation(atoms, config)
            job_ids.append(job_id)
        
        return job_ids
