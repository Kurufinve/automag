"""
Abstract interface for structure serialization.

Follows Interface Segregation Principle - focused interface.
"""

from abc import ABC, abstractmethod
from ase import Atoms


class StructureSerializer(ABC):
    """Interface for serializing/deserializing atomic structures."""
    
    @abstractmethod
    def serialize(self, atoms: Atoms) -> str:
        """
        Serialize Atoms object to string.
        
        Args:
            atoms: ASE Atoms object
            
        Returns:
            Serialized string representation
        """
        pass
    
    @abstractmethod
    def deserialize(self, data: str) -> Atoms:
        """
        Deserialize string to Atoms object.
        
        Args:
            data: Serialized string
            
        Returns:
            ASE Atoms object
        """
        pass
