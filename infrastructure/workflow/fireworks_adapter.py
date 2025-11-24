"""
FireWorks-based workflow submission adapter.

Refactored from common/SubmitFirework.py to follow SOLID principles.

Follows:
- Single Responsibility Principle: Only handles FireWorks submission
- Open/Closed Principle: Implements WorkflowSubmitter interface
- Liskov Substitution Principle: Can substitute for any WorkflowSubmitter
- Dependency Inversion Principle: Depends on abstractions
"""

from typing import List, Dict
from ase import Atoms
from fireworks import LaunchPad, Firework, Workflow

from core.interfaces.workflow import (
    WorkflowSubmitter, 
    CalculationConfig,
    CalculationInputWriter
)
from core.interfaces.serialization import StructureSerializer
from common.utilities import VaspCalculationTask, WriteOutputTask, WriteChargesTask


class FireworksWorkflowSubmitter(WorkflowSubmitter):
    """
    FireWorks-based implementation of WorkflowSubmitter.
    
    This class follows the Adapter pattern to integrate FireWorks
    while conforming to the WorkflowSubmitter interface.
    """
    
    def __init__(self, 
                 launchpad_file: str,
                 serializer: StructureSerializer):
        """
        Initialize FireWorks submitter.
        
        Args:
            launchpad_file: Path to LaunchPad YAML file
            serializer: Structure serializer implementation
        """
        self._launchpad = LaunchPad.from_file(launchpad_file)
        self._serializer = serializer
    
    def submit_single_calculation(self, atoms: Atoms, config: CalculationConfig) -> str:
        """
        Submit a single calculation via FireWorks.
        
        Args:
            atoms: ASE Atoms object
            config: Calculation configuration
            
        Returns:
            Workflow name/ID
        """
        encode = self._serializer.serialize(atoms)
        
        # Create single firework
        task = VaspCalculationTask(
            calc_params=config.calc_params,
            encode=encode,
            magmoms=config.magmoms,
            pert_step=config.pert_step,
            pert_value=config.pert_value,
            dummy_atom=config.dummy_atom,
            atom_ucalc=config.atom_ucalc
        )
        
        firework = Firework([task], name=config.name, spec={'_pass_job_info': True})
        workflow = Workflow([firework], name=config.name)
        
        self._launchpad.add_wf(workflow)
        return workflow.name
    
    def submit_workflow(self, atoms: Atoms, configs: List[CalculationConfig], 
                       workflow_name: str) -> str:
        """
        Submit a multi-step workflow via FireWorks.
        
        Args:
            atoms: ASE Atoms object
            configs: List of calculation configurations
            workflow_name: Name for the workflow
            
        Returns:
            Workflow name/ID
        """
        encode = self._serializer.serialize(atoms)
        fireworks = []
        
        for i, config in enumerate(configs):
            if config.encode is None:
                config_dict = {
                    'calc_params': config.calc_params,
                    'encode': encode,
                    'magmoms': config.magmoms
                }
            else:
                config_dict = {
                    'calc_params': config.calc_params,
                    'magmoms': config.magmoms
                }
            
            if config.pert_step:
                config_dict.update({
                    'pert_step': config.pert_step,
                    'pert_value': config.pert_value,
                    'dummy_atom': config.dummy_atom,
                    'atom_ucalc': config.atom_ucalc
                })
            
            task = VaspCalculationTask(**config_dict)
            fw = Firework([task], name=config.name, 
                         spec={'_pass_job_info': True}, fw_id=i)
            fireworks.append(fw)
        
        # Create workflow with dependencies (each FW depends on previous)
        links_dict = {}
        for i in range(len(fireworks) - 1):
            links_dict[fireworks[i].fw_id] = [fireworks[i + 1].fw_id]
        
        workflow = Workflow(fireworks, name=workflow_name, links_dict=links_dict)
        self._launchpad.add_wf(workflow)
        
        return workflow.name
    
    def check_status(self, job_id: str) -> str:
        """
        Check status of a FireWorks workflow.
        
        Args:
            job_id: Workflow name
            
        Returns:
            Status string
        """
        # Query LaunchPad for workflow status
        wf_ids = self._launchpad.get_wf_ids({'name': job_id})
        if not wf_ids:
            return 'NOT_FOUND'
        
        wf = self._launchpad.get_wf_by_fw_id(wf_ids[0])
        return wf.state
