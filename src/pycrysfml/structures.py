"""
Crystallographic data structures for PyCrysFML08

This module defines Python dataclasses representing crystallographic entities:
- Crystal unit cells
- Atoms and their positions
- Space groups
- Magnetic structures

These structures are designed to be compatible with CrysFML08's Fortran types
while providing a Pythonic interface.
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Optional
import math
import numpy as np


@dataclass
class CrystalCell:
    """
    Crystal unit cell parameters.

    Attributes:
        a, b, c: Unit cell lengths in Angstroms
        alpha, beta, gamma: Unit cell angles in degrees

    Example:
        >>> cell = CrystalCell(a=4.76, b=4.76, c=13.02, alpha=90, beta=90, gamma=120)
        >>> print(f"Volume: {cell.volume:.2f} Å³")
        Volume: 256.62 Å³
    """
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0

    @property
    def lengths(self) -> Tuple[float, float, float]:
        """Return cell lengths as a tuple (a, b, c)."""
        return (self.a, self.b, self.c)

    @property
    def angles(self) -> Tuple[float, float, float]:
        """Return cell angles as a tuple (alpha, beta, gamma)."""
        return (self.alpha, self.beta, self.gamma)

    @property
    def volume(self) -> float:
        """
        Calculate unit cell volume using the metric tensor.

        Returns:
            Volume in cubic Angstroms
        """
        a, b, c = self.a, self.b, self.c
        alpha_rad = math.radians(self.alpha)
        beta_rad = math.radians(self.beta)
        gamma_rad = math.radians(self.gamma)

        cos_a = math.cos(alpha_rad)
        cos_b = math.cos(beta_rad)
        cos_g = math.cos(gamma_rad)

        vol = a * b * c * math.sqrt(
            1 - cos_a**2 - cos_b**2 - cos_g**2 + 2*cos_a*cos_b*cos_g
        )
        return vol

    @property
    def metric_tensor(self) -> np.ndarray:
        """
        Calculate the metric tensor G for the unit cell.

        The metric tensor is used for calculating distances and angles
        in the crystal coordinate system.

        Returns:
            3x3 numpy array representing the metric tensor
        """
        a, b, c = self.a, self.b, self.c
        alpha_rad = math.radians(self.alpha)
        beta_rad = math.radians(self.beta)
        gamma_rad = math.radians(self.gamma)

        cos_a = math.cos(alpha_rad)
        cos_b = math.cos(beta_rad)
        cos_g = math.cos(gamma_rad)

        G = np.array([
            [a*a,         a*b*cos_g,   a*c*cos_b],
            [a*b*cos_g,   b*b,         b*c*cos_a],
            [a*c*cos_b,   b*c*cos_a,   c*c]
        ])
        return G

    def d_spacing(self, h: int, k: int, l: int) -> float:
        """
        Calculate d-spacing for given Miller indices.

        Args:
            h, k, l: Miller indices

        Returns:
            d-spacing in Angstroms
        """
        # For general triclinic case, use reciprocal metric tensor
        G = self.metric_tensor
        G_inv = np.linalg.inv(G)
        hkl = np.array([h, k, l])
        d_inv_sq = hkl @ G_inv @ hkl
        return 1.0 / math.sqrt(d_inv_sq)

    def crystal_system(self) -> str:
        """
        Determine the crystal system from cell parameters.

        Returns:
            One of: 'cubic', 'tetragonal', 'orthorhombic', 'hexagonal',
                   'trigonal', 'monoclinic', 'triclinic'
        """
        tol = 0.01  # Tolerance for angle comparisons

        a_eq_b = abs(self.a - self.b) < tol * self.a
        b_eq_c = abs(self.b - self.c) < tol * self.b
        a_eq_c = abs(self.a - self.c) < tol * self.a

        alpha_90 = abs(self.alpha - 90) < tol
        beta_90 = abs(self.beta - 90) < tol
        gamma_90 = abs(self.gamma - 90) < tol
        gamma_120 = abs(self.gamma - 120) < tol

        if a_eq_b and b_eq_c and alpha_90 and beta_90 and gamma_90:
            return "cubic"
        elif a_eq_b and alpha_90 and beta_90 and gamma_90:
            return "tetragonal"
        elif alpha_90 and beta_90 and gamma_90:
            return "orthorhombic"
        elif a_eq_b and alpha_90 and beta_90 and gamma_120:
            return "hexagonal"
        elif a_eq_b and a_eq_c and abs(self.alpha - self.beta) < tol and abs(self.beta - self.gamma) < tol:
            return "trigonal"
        elif (alpha_90 and gamma_90) or (alpha_90 and beta_90) or (beta_90 and gamma_90):
            return "monoclinic"
        else:
            return "triclinic"


@dataclass
class Atom:
    """
    Atom in a crystal structure.

    Attributes:
        label: Unique identifier for the atom (e.g., "Fe1", "O2")
        symbol: Element symbol (e.g., "Fe", "O")
        x, y, z: Fractional coordinates
        occupancy: Site occupancy (0 to 1)
        b_iso: Isotropic displacement parameter (B-factor) in Å²
        multiplicity: Site multiplicity from space group

    Example:
        >>> fe = Atom(label="Fe1", symbol="Fe", x=0.0, y=0.0, z=0.0)
        >>> print(fe)
        Atom(Fe1: Fe at (0.000, 0.000, 0.000))
    """
    label: str
    symbol: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    b_iso: float = 0.0
    multiplicity: int = 1

    @property
    def position(self) -> Tuple[float, float, float]:
        """Return fractional coordinates as a tuple."""
        return (self.x, self.y, self.z)

    def __str__(self) -> str:
        return f"Atom({self.label}: {self.symbol} at ({self.x:.3f}, {self.y:.3f}, {self.z:.3f}))"


@dataclass
class MagneticAtom(Atom):
    """
    Magnetic atom with moment information.

    Extends Atom with magnetic moment components and propagation vector
    information for describing magnetic structures.

    Attributes:
        All attributes from Atom, plus:
        moment: Magnetic moment components [mx, my, mz] in Bohr magnetons
        form_factor: Magnetic form factor name (e.g., "JFE3" for Fe3+)
        phase: Phase factor for incommensurate structures

    Example:
        >>> fe_mag = MagneticAtom(
        ...     label="Fe1", symbol="Fe", x=0, y=0, z=0,
        ...     moment=[0, 0, 4.5], form_factor="JFE3"
        ... )
    """
    moment: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    form_factor: str = ""
    phase: float = 0.0

    @property
    def moment_magnitude(self) -> float:
        """Calculate the magnitude of the magnetic moment."""
        return math.sqrt(sum(m**2 for m in self.moment))


@dataclass
class SpaceGroup:
    """
    Space group information.

    Attributes:
        number: International Tables number (1-230)
        symbol: Hermann-Mauguin symbol (e.g., "Fm-3m")
        hall_symbol: Hall symbol (optional)
        crystal_system: Crystal system name

    Example:
        >>> sg = SpaceGroup(number=225, symbol="Fm-3m", crystal_system="cubic")
    """
    number: int
    symbol: str
    hall_symbol: str = ""
    crystal_system: str = ""

    def __str__(self) -> str:
        return f"SpaceGroup({self.number}: {self.symbol})"


@dataclass
class PropagationVector:
    """
    Magnetic propagation vector k.

    Describes the periodicity of a magnetic structure relative to
    the nuclear unit cell.

    Attributes:
        kx, ky, kz: Components in reciprocal lattice units

    Example:
        >>> k = PropagationVector(0.5, 0, 0.5)  # AFM doubling along a and c
    """
    kx: float = 0.0
    ky: float = 0.0
    kz: float = 0.0

    @property
    def components(self) -> Tuple[float, float, float]:
        """Return k-vector components as a tuple."""
        return (self.kx, self.ky, self.kz)

    @property
    def is_commensurate(self) -> bool:
        """Check if the propagation vector is commensurate."""
        tol = 1e-6
        for k in [self.kx, self.ky, self.kz]:
            # Check if k is a simple fraction (n/m where m <= 12)
            for m in range(1, 13):
                if abs(k * m - round(k * m)) < tol:
                    break
            else:
                return False
        return True


@dataclass
class MagneticStructure:
    """
    Complete magnetic structure description.

    Contains all information needed to describe a magnetic structure:
    the parent crystal structure, propagation vector(s), and magnetic atoms.

    Attributes:
        cell: Crystal unit cell
        spacegroup: Parent space group
        magnetic_spacegroup: Magnetic space group symbol (optional)
        kvectors: List of propagation vectors
        atoms: List of non-magnetic atoms
        magnetic_atoms: List of magnetic atoms with moments

    Example:
        >>> mag = MagneticStructure(
        ...     cell=CrystalCell(3.75, 5.73, 11.27, 90, 90, 90),
        ...     spacegroup=SpaceGroup(71, "Immm"),
        ...     kvectors=[PropagationVector(0.5, 0, 0.5)],
        ...     magnetic_atoms=[ho_atom, ni_atom]
        ... )
    """
    cell: CrystalCell
    spacegroup: SpaceGroup
    magnetic_spacegroup: str = ""
    kvectors: List[PropagationVector] = field(default_factory=list)
    atoms: List[Atom] = field(default_factory=list)
    magnetic_atoms: List[MagneticAtom] = field(default_factory=list)

    @property
    def is_commensurate(self) -> bool:
        """Check if all propagation vectors are commensurate."""
        return all(k.is_commensurate for k in self.kvectors)

    @property
    def total_moment(self) -> float:
        """Calculate total magnetic moment per formula unit."""
        return sum(atom.moment_magnitude for atom in self.magnetic_atoms)
