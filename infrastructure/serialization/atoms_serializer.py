"""
JSON-based Atoms serialization.

Refactored from common/utilities.py - follows Single Responsibility Principle.
"""

import json
from ase import Atoms
from core.interfaces.serialization import StructureSerializer


class JSONAtomsSerializer(StructureSerializer):
    """
    JSON-based implementation of structure serialization.
    
    Follows:
    - Single Responsibility Principle: Only handles serialization
    - Open/Closed Principle: Can be extended with different formats
    """
    
    def serialize(self, atoms: Atoms) -> str:
        """
        Serialize Atoms object to JSON string.
        
        Args:
            atoms: ASE Atoms object
            
        Returns:
            JSON string representation
        """
        data = {
            'cell': atoms.get_cell().tolist(),
            'scaled_positions': atoms.get_scaled_positions().tolist(),
            'numbers': atoms.get_atomic_numbers().tolist(),
            'pbc': atoms.get_pbc().tolist(),
        }
        return json.dumps(data)
    
    def deserialize(self, data: str) -> Atoms:
        """
        Deserialize JSON string to Atoms object.
        
        Args:
            data: JSON string
            
        Returns:
            ASE Atoms object
        """
        try:
            parsed_data = json.loads(data, encoding='utf-8')
        except TypeError:
            parsed_data = json.loads(data)
        
        atoms = Atoms(
            cell=parsed_data['cell'],
            scaled_positions=parsed_data['scaled_positions'],
            numbers=parsed_data['numbers'],
            pbc=parsed_data['pbc'],
        )
        
        return atoms
