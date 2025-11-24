"""
Factory for creating workflow submitters.

Follows:
- Open/Closed Principle: Can add new submitter types without modifying factory
- Dependency Inversion Principle: Returns abstractions
"""

from typing import Optional
from pathlib import Path

from core.interfaces.workflow import WorkflowSubmitter
from infrastructure.workflow.fireworks_adapter import FireworksWorkflowSubmitter
from infrastructure.workflow.manual_adapter import ManualJobSubmitter
from infrastructure.serialization.atoms_serializer import JSONAtomsSerializer
from infrastructure.calculators.vasp_adapter import VaspInputWriter


class WorkflowSubmitterFactory:
    """
    Factory for creating workflow submitters.
    
    This follows the Factory pattern and Open/Closed Principle.
    New submitter types can be added without modifying existing code.
    """
    
    @staticmethod
    def create_fireworks_submitter(launchpad_file: str) -> WorkflowSubmitter:
        """
        Create FireWorks-based submitter.
        
        Args:
            launchpad_file: Path to LaunchPad YAML file
            
        Returns:
            WorkflowSubmitter implementation
        """
        serializer = JSONAtomsSerializer()
        return FireworksWorkflowSubmitter(launchpad_file, serializer)
    
    @staticmethod
    def create_manual_submitter(
        calcfold_path: str,
        jobheader: str = "#!/bin/bash",
        calculator_command: str = "mpirun vasp_std",
        environment_activate: str = "",
        environment_deactivate: str = "",
        queue_system: str = "slurm"
    ) -> WorkflowSubmitter:
        """
        Create manual job submission submitter.
        
        Args:
            calcfold_path: Base path for calculations
            jobheader: Job script header
            calculator_command: Command to run calculator
            environment_activate: Environment activation command
            environment_deactivate: Environment deactivation command
            queue_system: Queue system type
            
        Returns:
            WorkflowSubmitter implementation
        """
        input_writer = VaspInputWriter()
        
        return ManualJobSubmitter(
            input_writer=input_writer,
            calcfold_path=calcfold_path,
            jobheader=jobheader,
            calculator_command=calculator_command,
            environment_activate=environment_activate,
            environment_deactivate=environment_deactivate,
            queue_system=queue_system
        )
    
    @staticmethod
    def create_from_config(
        use_fireworks: bool,
        launchpad_file: Optional[str] = None,
        calcfold_path: Optional[str] = None,
        **kwargs
    ) -> WorkflowSubmitter:
        """
        Create submitter based on configuration.
        
        Args:
            use_fireworks: Whether to use FireWorks
            launchpad_file: LaunchPad file for FireWorks
            calcfold_path: Calculation folder for manual submission
            **kwargs: Additional parameters for manual submission
            
        Returns:
            WorkflowSubmitter implementation
        """
        if use_fireworks:
            if launchpad_file is None:
                raise ValueError("launchpad_file required for FireWorks")
            return WorkflowSubmitterFactory.create_fireworks_submitter(launchpad_file)
        else:
            if calcfold_path is None:
                import os
                calcfold_path = os.path.join(
                    os.environ.get('AUTOMAG_PATH', '.'),
                    'CalcFold'
                )
            return WorkflowSubmitterFactory.create_manual_submitter(
                calcfold_path, **kwargs
            )
