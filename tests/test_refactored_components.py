"""
Unit tests for refactored components.

These tests demonstrate the improved testability from SOLID refactoring.
"""

import pytest
from unittest.mock import Mock, MagicMock
import numpy as np
from ase import Atoms

from core.domain.calculation import CalculationParameters, MagneticConfiguration
from core.services.calculation_service import CalculationService
from core.interfaces.workflow import CalculationConfig
from infrastructure.serialization.atoms_serializer import JSONAtomsSerializer


class TestCalculationParameters:
    """Test domain models (Single Responsibility Principle)."""
    
    def test_parameters_validation_positive_encut(self):
        """Test ENCUT validation."""
        params = CalculationParameters(encut=500)
        assert params.validate() is True
    
    def test_parameters_validation_negative_encut_raises(self):
        """Test ENCUT validation with negative value."""
        params = CalculationParameters(encut=-100)
        with pytest.raises(ValueError, match="ENCUT must be positive"):
            params.validate()
    
    def test_parameters_to_dict_excludes_none(self):
        """Test conversion to dict excludes None values."""
        params = CalculationParameters(encut=500, sigma=0.1, kpts=None)
        result = params.to_dict()
        assert 'encut' in result
        assert 'sigma' in result
        assert 'kpts' not in result


class TestMagneticConfiguration:
    """Test magnetic configuration models."""
    
    def test_is_ferromagnetic(self):
        """Test ferromagnetic detection."""
        config = MagneticConfiguration(magmoms=[5.0, 5.0, 5.0])
        assert config.is_ferromagnetic() is True
    
    def test_is_antiferromagnetic(self):
        """Test antiferromagnetic detection."""
        config = MagneticConfiguration(magmoms=[5.0, -5.0, 5.0, -5.0])
        assert config.is_antiferromagnetic() is True
    
    def test_is_nonmagnetic(self):
        """Test non-magnetic detection."""
        config = MagneticConfiguration(magmoms=[0.0, 0.0, 0.0])
        assert config.is_nonmagnetic() is True


class TestJSONAtomsSerializer:
    """Test serialization (Interface Segregation Principle)."""
    
    def test_serialize_deserialize_roundtrip(self):
        """Test that serialization/deserialization preserves data."""
        # Create simple atoms object
        atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        
        serializer = JSONAtomsSerializer()
        
        # Serialize
        serialized = serializer.serialize(atoms)
        assert isinstance(serialized, str)
        
        # Deserialize
        atoms_restored = serializer.deserialize(serialized)
        
        # Check equality
        assert len(atoms) == len(atoms_restored)
        assert (atoms.numbers == atoms_restored.numbers).all()


class TestCalculationService:
    """Test service layer (Dependency Inversion Principle)."""
    
    def test_submit_convergence_test_calls_submitter(self):
        """Test that service delegates to submitter correctly."""
        # Create mock submitter
        mock_submitter = Mock()
        mock_submitter.submit_single_calculation.return_value = "job_123"
        
        # Create service with injected dependency
        service = CalculationService(mock_submitter)
        
        # Create test data
        atoms = Atoms('Fe', positions=[[0, 0, 0]])
        params = CalculationParameters(sigma=0.1, kpts=40)
        magmoms = [5.0]
        test_values = [500, 600, 700]
        
        # Execute
        job_ids = service.submit_convergence_test(
            atoms, params, magmoms, 'encut', test_values, 'encut'
        )
        
        # Verify
        assert len(job_ids) == 3
        assert mock_submitter.submit_single_calculation.call_count == 3
    
    def test_submit_singlepoint_with_recalc_creates_workflow(self):
        """Test that recalc option creates workflow."""
        mock_submitter = Mock()
        mock_submitter.submit_workflow.return_value = "workflow_123"
        
        service = CalculationService(mock_submitter)
        
        atoms = Atoms('Fe', positions=[[0, 0, 0]])
        params = CalculationParameters(encut=500)
        magmoms = [5.0]
        
        # Execute with recalc
        job_id = service.submit_singlepoint_calculation(
            atoms, params, magmoms, 'test', with_recalc=True
        )
        
        # Verify workflow was submitted
        assert mock_submitter.submit_workflow.called
        assert job_id == "workflow_123"


class TestLiskovSubstitution:
    """
    Test Liskov Substitution Principle.
    
    Both FireWorks and Manual submitters should be interchangeable.
    """
    
    def test_submitters_have_same_interface(self):
        """Test that both submitters implement the same methods."""
        from infrastructure.workflow.fireworks_adapter import FireworksWorkflowSubmitter
        from infrastructure.workflow.manual_adapter import ManualJobSubmitter
        
        # Check they have the required methods
        required_methods = [
            'submit_single_calculation',
            'submit_workflow',
            'check_status'
        ]
        
        for method in required_methods:
            assert hasattr(FireworksWorkflowSubmitter, method)
            assert hasattr(ManualJobSubmitter, method)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
