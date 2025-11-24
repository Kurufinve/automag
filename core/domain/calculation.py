"""
Domain models for calculations.

Follows Single Responsibility Principle - each class has one clear purpose.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Union
from enum import Enum


class CalculationMode(Enum):
    """Enumeration of calculation modes."""
    ENCUT = "encut"
    KGRID = "kgrid"
    PERTURBATIONS = "perturbations"
    SINGLEPOINT = "singlepoint"
    MAE_THETA_PHI = "mae_theta_phi"
    MAE_CURVE = "mae_curve"


@dataclass
class CalculationParameters:
    """
    Immutable calculation parameters.
    
    Follows Single Responsibility Principle - only manages parameters.
    """
    encut: Optional[float] = None
    sigma: Optional[float] = None
    ismear: Optional[int] = None
    kpts: Optional[Union[int, float]] = None
    xc: Optional[str] = None
    setups: Optional[Dict[str, str]] = None
    prec: Optional[str] = None
    algo: Optional[str] = None
    ediff: Optional[float] = None
    ediffg: Optional[float] = None
    nelm: Optional[int] = None
    ibrion: Optional[int] = None
    isif: Optional[int] = None
    nsw: Optional[int] = None
    lreal: Optional[Union[bool, str]] = None
    ncore: Optional[int] = None
    npar: Optional[int] = None
    ldauu: Optional[List[float]] = None
    ldauj: Optional[List[float]] = None
    ldaul: Optional[List[int]] = None
    ldau: Optional[bool] = None
    ldautype: Optional[int] = None
    lorbit: Optional[int] = None
    ispin: Optional[int] = None
    lnoncollinear: Optional[bool] = None
    lsorbit: Optional[bool] = None
    icharg: Optional[int] = None
    istart: Optional[int] = None
    lcharg: Optional[bool] = None
    lwave: Optional[bool] = None
    saxis: Optional[List[float]] = None
    voskown: Optional[int] = None
    gga_compat: Optional[bool] = None
    lmaxmix: Optional[int] = None
    amix: Optional[float] = None
    bmix: Optional[float] = None
    amix_mag: Optional[float] = None
    bmix_mag: Optional[float] = None
    lasph: Optional[bool] = None
    
    def validate(self) -> bool:
        """
        Validate parameter constraints.
        
        Returns:
            True if valid
            
        Raises:
            ValueError: If validation fails
        """
        if self.encut is not None and self.encut < 0:
            raise ValueError("ENCUT must be positive")
        if self.sigma is not None and self.sigma < 0:
            raise ValueError("SIGMA must be positive")
        if self.kpts is not None and self.kpts <= 0:
            raise ValueError("KPTS must be positive")
        return True
    
    def to_dict(self) -> Dict:
        """Convert to dictionary, excluding None values."""
        return {k: v for k, v in self.__dict__.items() if v is not None}


@dataclass
class MagneticConfiguration:
    """
    Magnetic configuration for a structure.
    
    Follows Single Responsibility Principle.
    """
    magmoms: List[float]
    magnetic_atoms: Optional[List[str]] = None
    
    def is_ferromagnetic(self) -> bool:
        """Check if configuration is ferromagnetic."""
        return all(m >= 0 for m in self.magmoms) or all(m <= 0 for m in self.magmoms)
    
    def is_antiferromagnetic(self) -> bool:
        """Check if configuration is antiferromagnetic."""
        return sum(self.magmoms) == 0 and not all(m == 0 for m in self.magmoms)
    
    def is_nonmagnetic(self) -> bool:
        """Check if configuration is non-magnetic."""
        return all(m == 0 for m in self.magmoms)
