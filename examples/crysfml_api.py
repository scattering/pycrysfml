"""
CrysFML08 Python API via gfort2py

This module provides a Python interface to the CrysFML08 Fortran library
using gfort2py bindings. It wraps the low-level Fortran functions in
more Pythonic interfaces.

Example usage:
    from crysfml_api import CrysFML

    cfml = CrysFML()

    # Get atomic mass
    mass = cfml.get_atomic_mass("Fe")

    # Calculate formula mass
    formula_mass = cfml.formula_mass({"Al": 2, "O": 3})
"""

import os
import sys
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import re

# Find the build directory (parent of examples)
BUILD_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BUILD_DIR)

try:
    import gfort2py as gf
except ImportError:
    print("ERROR: gfort2py not installed. Run: pip install gfort2py")
    sys.exit(1)

# Library path
LIB_PATH = os.path.join(BUILD_DIR, "libcrysfml_test.dylib")
if not os.path.exists(LIB_PATH):
    LIB_PATH = os.path.join(BUILD_DIR, "libcrysfml_test.so")


@dataclass
class CrystalCell:
    """Crystal unit cell parameters"""
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float

    @property
    def lengths(self) -> Tuple[float, float, float]:
        return (self.a, self.b, self.c)

    @property
    def angles(self) -> Tuple[float, float, float]:
        return (self.alpha, self.beta, self.gamma)

    @property
    def volume(self) -> float:
        """Calculate unit cell volume using the metric tensor"""
        import math
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


@dataclass
class Atom:
    """Atom in a crystal structure"""
    label: str
    symbol: str
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    b_iso: float = 0.0
    multiplicity: int = 1


@dataclass
class SpaceGroup:
    """Space group information"""
    number: int
    symbol: str
    crystal_system: str = ""


class CrysFML:
    """
    Python interface to CrysFML08 library

    Provides access to:
    - Scattering tables (atomic masses, scattering lengths, etc.)
    - Mathematical functions
    - Crystallographic calculations
    """

    def __init__(self, lib_path: str = None):
        """Initialize CrysFML interface"""
        if lib_path is None:
            lib_path = LIB_PATH

        if not os.path.exists(lib_path):
            raise FileNotFoundError(f"Library not found: {lib_path}")

        self._lib_path = lib_path
        self._modules = {}

    def _get_module(self, name: str):
        """Load a module lazily"""
        if name not in self._modules:
            mod_path = os.path.join(BUILD_DIR, f"{name}.mod")
            if not os.path.exists(mod_path):
                raise FileNotFoundError(f"Module file not found: {mod_path}")
            self._modules[name] = gf.fFort(self._lib_path, mod_path)
        return self._modules[name]

    # ==================== Scattering Tables ====================

    def get_atomic_mass(self, symbol: str) -> float:
        """
        Get atomic mass for an element

        Args:
            symbol: Element symbol (e.g., "Fe", "O", "Al")

        Returns:
            Atomic mass in atomic mass units (amu)
        """
        lib = self._get_module("cfml_scattering_tables")
        # Pad symbol to 2 characters as Fortran expects
        symbol = symbol.strip().ljust(2)
        result = lib.get_atomic_mass(symbol)
        return result.result

    def get_covalent_radius(self, symbol: str) -> float:
        """
        Get covalent radius for an element

        Args:
            symbol: Element symbol

        Returns:
            Covalent radius in Angstroms
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_covalent_radius(symbol)
        return result.result

    def get_fermi_length(self, symbol: str) -> float:
        """
        Get neutron Fermi scattering length for an element

        Args:
            symbol: Element symbol

        Returns:
            Fermi scattering length in units of 10^-12 cm
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_fermi_length(symbol)
        return result.result

    def get_atomic_volume(self, symbol: str) -> float:
        """
        Get atomic volume for an element

        Args:
            symbol: Element symbol

        Returns:
            Atomic volume
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_atomic_vol(symbol)
        return result.result

    def formula_mass(self, composition: Dict[str, int]) -> float:
        """
        Calculate formula mass for a chemical composition

        Args:
            composition: Dictionary of element symbols to counts
                        e.g., {"Al": 2, "O": 3} for Al2O3

        Returns:
            Formula mass in amu
        """
        total = 0.0
        for symbol, count in composition.items():
            mass = self.get_atomic_mass(symbol)
            total += mass * count
        return total

    # ==================== Math Functions ====================

    def factorial(self, n: int) -> int:
        """
        Calculate factorial of n

        Args:
            n: Non-negative integer

        Returns:
            n!
        """
        lib = self._get_module("cfml_maths")
        result = lib.factorial_i(n)
        return result.result

    def gcd(self, a: int, b: int) -> int:
        """
        Calculate greatest common divisor of a and b

        Args:
            a, b: Integers

        Returns:
            GCD(a, b)
        """
        lib = self._get_module("cfml_maths")
        result = lib.gcd(a, b)
        return result.result

    def lcm(self, a: int, b: int) -> int:
        """
        Calculate least common multiple of a and b

        Args:
            a, b: Positive integers

        Returns:
            LCM(a, b)
        """
        lib = self._get_module("cfml_maths")
        result = lib.lcm(a, b)
        return result.result

    # ==================== CIF Parsing ====================

    def read_cif(self, filename: str) -> Tuple[CrystalCell, SpaceGroup, List[Atom]]:
        """
        Parse a CIF file and extract crystal structure information

        Note: This is a Python-based parser, not using CrysFML's Fortran parser
        due to gfort2py limitations with complex derived types.

        Args:
            filename: Path to CIF file

        Returns:
            Tuple of (CrystalCell, SpaceGroup, list of Atoms)
        """
        with open(filename, 'r') as f:
            content = f.read()

        # Parse cell parameters
        a = self._parse_cif_value(content, '_cell_length_a')
        b = self._parse_cif_value(content, '_cell_length_b')
        c = self._parse_cif_value(content, '_cell_length_c')
        alpha = self._parse_cif_value(content, '_cell_angle_alpha')
        beta = self._parse_cif_value(content, '_cell_angle_beta')
        gamma = self._parse_cif_value(content, '_cell_angle_gamma')

        cell = CrystalCell(a, b, c, alpha, beta, gamma)

        # Parse space group
        sg_number = int(self._parse_cif_value(content, '_symmetry_Int_Tables_number'))
        sg_symbol = self._parse_cif_string(content, '_symmetry_space_group_name_H-M')

        spacegroup = SpaceGroup(sg_number, sg_symbol)

        # Parse atoms
        atoms = self._parse_cif_atoms(content)

        return cell, spacegroup, atoms

    def _parse_cif_value(self, content: str, tag: str) -> float:
        """Parse a numeric value from CIF, handling uncertainties in parentheses"""
        pattern = rf'{tag}\s+([\d.]+)(?:\(\d+\))?'
        match = re.search(pattern, content)
        if match:
            return float(match.group(1))
        raise ValueError(f"CIF tag not found: {tag}")

    def _parse_cif_string(self, content: str, tag: str) -> str:
        """Parse a string value from CIF"""
        pattern = rf"{tag}\s+['\"]?([^'\"\n]+)['\"]?"
        match = re.search(pattern, content)
        if match:
            return match.group(1).strip()
        return ""

    def _parse_cif_atoms(self, content: str) -> List[Atom]:
        """Parse atom site loop from CIF"""
        atoms = []

        # Find the atom_site loop
        atom_loop_match = re.search(
            r'loop_\s*\n((?:_atom_site_\w+\s*\n)+)((?:.*\n)*?)(?=loop_|#|$)',
            content
        )

        if not atom_loop_match:
            return atoms

        # Parse column headers
        headers_text = atom_loop_match.group(1)
        headers = re.findall(r'_atom_site_(\w+)', headers_text)

        # Parse data rows
        data_text = atom_loop_match.group(2)

        # Find column indices
        label_idx = headers.index('label') if 'label' in headers else None
        symbol_idx = headers.index('type_symbol') if 'type_symbol' in headers else None
        x_idx = headers.index('fract_x') if 'fract_x' in headers else None
        y_idx = headers.index('fract_y') if 'fract_y' in headers else None
        z_idx = headers.index('fract_z') if 'fract_z' in headers else None
        occ_idx = headers.index('occupancy') if 'occupancy' in headers else None
        mult_idx = headers.index('symmetry_multiplicity') if 'symmetry_multiplicity' in headers else None
        b_idx = headers.index('B_iso_or_equiv') if 'B_iso_or_equiv' in headers else None

        # Parse each atom line
        for line in data_text.strip().split('\n'):
            line = line.strip()
            if not line or line.startswith('_') or line.startswith('loop_') or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < len(headers):
                continue

            try:
                label = parts[label_idx] if label_idx is not None else "?"
                symbol = parts[symbol_idx] if symbol_idx is not None else label[:2].rstrip('0123456789+-')
                # Clean up symbol (remove charge like 3+ or 2-)
                symbol = re.sub(r'[0-9+-]+', '', symbol)

                x = self._parse_cif_coord(parts[x_idx]) if x_idx is not None else 0.0
                y = self._parse_cif_coord(parts[y_idx]) if y_idx is not None else 0.0
                z = self._parse_cif_coord(parts[z_idx]) if z_idx is not None else 0.0
                occ = float(parts[occ_idx].rstrip('.')) if occ_idx is not None else 1.0
                mult = int(parts[mult_idx]) if mult_idx is not None else 1
                b_iso = float(parts[b_idx]) if b_idx is not None and parts[b_idx] != '.' else 0.0

                atom = Atom(
                    label=label,
                    symbol=symbol,
                    x=x, y=y, z=z,
                    occupancy=occ,
                    multiplicity=mult,
                    b_iso=b_iso
                )
                atoms.append(atom)
            except (ValueError, IndexError):
                continue

        return atoms

    def _parse_cif_coord(self, value: str) -> float:
        """Parse a coordinate value, handling uncertainties"""
        # Remove uncertainty in parentheses
        value = re.sub(r'\(\d+\)', '', value)
        return float(value)


def read_diffraction_data(filename: str, skip_header: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read diffraction pattern data (2theta, intensity)

    Args:
        filename: Path to data file
        skip_header: Number of header lines to skip

    Returns:
        Tuple of (2theta array, intensity array)
    """
    data = np.loadtxt(filename, skiprows=skip_header)
    if data.ndim == 1:
        # Single column - assume regular 2theta spacing
        return None, data
    return data[:, 0], data[:, 1]


def read_background(filename: str, skip_header: int = 5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read background spline data

    Args:
        filename: Path to background file
        skip_header: Number of header lines to skip

    Returns:
        Tuple of (x positions, background values)
    """
    data = np.loadtxt(filename, skiprows=skip_header)
    return data[:, 0], data[:, 1]
