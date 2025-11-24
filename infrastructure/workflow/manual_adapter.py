"""
Manual job submission adapter (SLURM/PBS).

Refactored from common/SubmitManual.py to follow SOLID principles.

Follows:
- Single Responsibility Principle: Only handles manual job submission
- Open/Closed Principle: Implements WorkflowSubmitter interface
- Liskov Substitution Principle: Can substitute for any WorkflowSubmitter
- Dependency Inversion Principle: Depends on abstractions
"""

import os
from typing import List, Dict, Optional
from pathlib import Path
from ase import Atoms
from ase.calculators.vasp import Vasp

from core.interfaces.workflow import (
    WorkflowSubmitter,
    CalculationConfig,
    CalculationInputWriter
)
from core.domain.calculation import CalculationMode


class ManualJobSubmitter(WorkflowSubmitter):
    """
    Manual job submission implementation (SLURM/PBS).
    
    This adapter allows the same workflow interface to work with
    manual job submission systems.
    """
    
    def __init__(self,
                 input_writer: CalculationInputWriter,
                 calcfold_path: str,
                 jobheader: str = "#!/bin/bash",
                 calculator_command: str = "mpirun vasp_std",
                 environment_activate: str = "",
                 environment_deactivate: str = "",
                 queue_system: str = "slurm"):
        """
        Initialize manual job submitter.
        
        Args:
            input_writer: Implementation for writing calculation inputs
            calcfold_path: Base path for calculations
            jobheader: Job script header
            calculator_command: Command to run calculator
            environment_activate: Command to activate environment
            environment_deactivate: Command to deactivate environment
            queue_system: Queue system (slurm, pbs, etc.)
        """
        self._input_writer = input_writer
        self._calcfold_path = Path(calcfold_path)
        self._jobheader = jobheader
        self._calculator_command = calculator_command
        self._env_activate = environment_activate
        self._env_deactivate = environment_deactivate
        self._queue_system = queue_system
    
    def submit_single_calculation(self, atoms: Atoms, config: CalculationConfig) -> str:
        """
        Submit a single calculation manually.
        
        Args:
            atoms: ASE Atoms object
            config: Calculation configuration
            
        Returns:
            Job directory path
        """
        # Create work directory
        formula = atoms.get_chemical_formula(mode='metal')
        workdir = self._calcfold_path / formula / config.name
        workdir.mkdir(parents=True, exist_ok=True)
        
        # Set magnetic moments if provided
        if config.magmoms:
            atoms = self._input_writer.set_magnetic_moments(atoms, config.magmoms)
        
        # Write input files
        self._input_writer.write_input_files(
            atoms, config.calc_params, str(workdir / 'singlepoint')
        )
        
        # Write job script
        self._write_jobscript(workdir, config.name, mode='singlepoint')
        
        return str(workdir)
    
    def submit_workflow(self, atoms: Atoms, configs: List[CalculationConfig],
                       workflow_name: str) -> str:
        """
        Submit a multi-step workflow manually.
        
        Args:
            atoms: ASE Atoms object
            configs: List of calculation configurations
            workflow_name: Name for the workflow
            
        Returns:
            Job directory path
        """
        # Create work directory
        formula = atoms.get_chemical_formula(mode='metal')
        workdir = self._calcfold_path / formula / workflow_name
        workdir.mkdir(parents=True, exist_ok=True)
        
        # Write inputs for each step
        for config in configs:
            if config.magmoms:
                atoms_with_mag = self._input_writer.set_magnetic_moments(
                    atoms, config.magmoms
                )
            else:
                atoms_with_mag = atoms
            
            step_dir = workdir / config.name
            self._input_writer.write_input_files(
                atoms_with_mag, config.calc_params, str(step_dir)
            )
        
        # Write master job script
        self._write_workflow_jobscript(workdir, configs, workflow_name)
        
        return str(workdir)
    
    def check_status(self, job_id: str) -> str:
        """
        Check status of a manual job.
        
        Args:
            job_id: Job directory path
            
        Returns:
            Status string
        """
        job_path = Path(job_id)
        
        # Check if output files exist
        if (job_path / 'OUTCAR').exists():
            # Check if completed
            with open(job_path / 'OUTCAR', 'r') as f:
                content = f.read()
                if 'General timing and accounting informations' in content:
                    return 'COMPLETED'
                else:
                    return 'RUNNING'
        
        return 'PENDING'
    
    def _write_jobscript(self, workdir: Path, job_name: str, mode: str) -> None:
        """Write job submission script."""
        jobscript_path = workdir / 'jobscript'
        
        with open(jobscript_path, 'w') as f:
            f.write(self._jobheader + '\n')
            
            if '#SBATCH' in self._jobheader:
                f.write(f'#SBATCH --job-name={job_name}\n')
            
            f.write('\n')
            f.write(f'cd {workdir / mode}\n')
            
            if self._env_activate:
                f.write(f'{self._env_activate}\n')
            
            f.write(f'{self._calculator_command}\n')
            
            if self._env_deactivate:
                f.write(f'{self._env_deactivate}\n')
    
    def _write_workflow_jobscript(self, workdir: Path, 
                                   configs: List[CalculationConfig],
                                   workflow_name: str) -> None:
        """Write job script for multi-step workflow."""
        jobscript_path = workdir / 'jobscript'
        
        with open(jobscript_path, 'w') as f:
            f.write(self._jobheader + '\n')
            
            if '#SBATCH' in self._jobheader:
                f.write(f'#SBATCH --job-name={workflow_name}\n')
            
            f.write('\n')
            
            for config in configs:
                f.write(f'cd {workdir / config.name}\n')
                
                if self._env_activate:
                    f.write(f'{self._env_activate}\n')
                
                f.write(f'{self._calculator_command}\n')
                
                if self._env_deactivate:
                    f.write(f'{self._env_deactivate}\n')
                
                f.write('\n')
