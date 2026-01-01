"""
CrysFML08 API - Python interface to the Fortran library

This module provides the main CrysFML class that interfaces with the
compiled CrysFML08 Fortran library using gfort2py.

The CrysFML class provides access to:
- Scattering tables (atomic masses, scattering lengths, form factors)
- Mathematical functions (factorial, gcd, lcm)
- Crystallographic calculations

Note: Some CrysFML08 functionality involving complex Fortran derived types
(CLASS polymorphism) is not accessible through gfort2py. For these features,
pure Python implementations are provided in other modules.
"""

import os
import sys
from typing import Dict, Optional
from pathlib import Path

# Try to import gfort2py
try:
    import gfort2py as gf
    HAS_GFORT2PY = True
except ImportError:
    HAS_GFORT2PY = False
    gf = None


def find_library() -> Optional[Path]:
    """
    Find the CrysFML shared library.

    Searches in the following locations:
    1. PYCRYSFML_LIB environment variable
    2. Package directory
    3. Current working directory

    Returns:
        Path to library file, or None if not found
    """
    # Check environment variable
    env_path = os.environ.get('PYCRYSFML_LIB')
    if env_path and os.path.exists(env_path):
        return Path(env_path)

    # Check package directory
    pkg_dir = Path(__file__).parent
    for name in ['libcrysfml.dylib', 'libcrysfml.so', 'libcrysfml_test.dylib', 'libcrysfml_test.so']:
        lib_path = pkg_dir / name
        if lib_path.exists():
            return lib_path

    # Check parent directories (for development)
    for parent in [pkg_dir.parent, pkg_dir.parent.parent]:
        for name in ['libcrysfml.dylib', 'libcrysfml.so', 'libcrysfml_test.dylib', 'libcrysfml_test.so']:
            lib_path = parent / name
            if lib_path.exists():
                return lib_path

    # Check current directory
    cwd = Path.cwd()
    for name in ['libcrysfml.dylib', 'libcrysfml.so', 'libcrysfml_test.dylib', 'libcrysfml_test.so']:
        lib_path = cwd / name
        if lib_path.exists():
            return lib_path

    return None


def find_modules() -> Optional[Path]:
    """
    Find the directory containing Fortran .mod files.

    Returns:
        Path to directory containing .mod files, or None if not found
    """
    # Check environment variable
    env_path = os.environ.get('PYCRYSFML_MODS')
    if env_path and os.path.exists(env_path):
        return Path(env_path)

    # Check package directory and parents
    pkg_dir = Path(__file__).parent
    for check_dir in [pkg_dir, pkg_dir.parent, pkg_dir.parent.parent]:
        if (check_dir / 'cfml_maths.mod').exists():
            return check_dir

    # Check current directory
    if (Path.cwd() / 'cfml_maths.mod').exists():
        return Path.cwd()

    return None


class CrysFML:
    """
    Python interface to CrysFML08 library.

    Provides access to crystallographic functions from the CrysFML08
    Fortran library. Uses gfort2py for the Fortran-Python interface.

    Attributes:
        lib_path: Path to the shared library
        mod_dir: Path to the directory containing .mod files

    Example:
        >>> from pycrysfml import CrysFML
        >>> cfml = CrysFML()
        >>> mass = cfml.get_atomic_mass("Fe")
        >>> print(f"Iron: {mass:.3f} amu")
        Iron: 55.847 amu

    Note:
        The library and module files must be built before using this class.
        See the documentation for build instructions.
    """

    def __init__(self, lib_path: str = None, mod_dir: str = None):
        """
        Initialize CrysFML interface.

        Args:
            lib_path: Path to shared library (optional, auto-detected if not provided)
            mod_dir: Path to directory containing .mod files (optional)

        Raises:
            ImportError: If gfort2py is not installed
            FileNotFoundError: If library or module files are not found
        """
        if not HAS_GFORT2PY:
            raise ImportError(
                "gfort2py is required but not installed. "
                "Install with: pip install gfort2py"
            )

        # Find library
        if lib_path:
            self._lib_path = Path(lib_path)
        else:
            self._lib_path = find_library()

        if self._lib_path is None or not self._lib_path.exists():
            raise FileNotFoundError(
                "CrysFML shared library not found. "
                "Set PYCRYSFML_LIB environment variable or provide lib_path."
            )

        # Find module directory
        if mod_dir:
            self._mod_dir = Path(mod_dir)
        else:
            self._mod_dir = find_modules()
            if self._mod_dir is None:
                # Try same directory as library
                self._mod_dir = self._lib_path.parent

        # Cache for loaded modules
        self._modules: Dict[str, any] = {}

    def _get_module(self, name: str):
        """
        Load a Fortran module lazily.

        Args:
            name: Module name (without .mod extension)

        Returns:
            gfort2py module interface
        """
        if name not in self._modules:
            mod_path = self._mod_dir / f"{name}.mod"
            if not mod_path.exists():
                raise FileNotFoundError(f"Module file not found: {mod_path}")
            self._modules[name] = gf.fFort(str(self._lib_path), str(mod_path))
        return self._modules[name]

    # ==================== Scattering Tables ====================

    def get_atomic_mass(self, symbol: str) -> float:
        """
        Get atomic mass for an element.

        Args:
            symbol: Element symbol (e.g., "Fe", "O", "Al")

        Returns:
            Atomic mass in atomic mass units (amu)

        Example:
            >>> cfml.get_atomic_mass("Fe")
            55.847
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_atomic_mass(symbol)
        return result.result

    def get_covalent_radius(self, symbol: str) -> float:
        """
        Get covalent radius for an element.

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
        Get neutron Fermi scattering length for an element.

        The Fermi scattering length (also called coherent scattering length)
        is the bound coherent scattering length for thermal neutrons.

        Args:
            symbol: Element symbol

        Returns:
            Fermi scattering length in units of 10^-12 cm (fm * 10)

        Note:
            Some elements (H, Li, Ti, V, Mn) have negative scattering lengths.
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_fermi_length(symbol)
        return result.result

    def get_atomic_volume(self, symbol: str) -> float:
        """
        Get atomic volume for an element.

        Args:
            symbol: Element symbol

        Returns:
            Atomic volume in cubic Angstroms
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_atomic_vol(symbol)
        return result.result

    def get_absorption_xs(self, symbol: str) -> float:
        """
        Get neutron absorption cross section for an element.

        Args:
            symbol: Element symbol

        Returns:
            Absorption cross section in barns
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_abs_xs(symbol)
        return result.result

    def get_incoherent_xs(self, symbol: str) -> float:
        """
        Get neutron incoherent scattering cross section for an element.

        Args:
            symbol: Element symbol

        Returns:
            Incoherent cross section in barns
        """
        lib = self._get_module("cfml_scattering_tables")
        symbol = symbol.strip().ljust(2)
        result = lib.get_inc_xs(symbol)
        return result.result

    def formula_mass(self, composition: Dict[str, int]) -> float:
        """
        Calculate formula mass for a chemical composition.

        Args:
            composition: Dictionary of element symbols to counts
                        e.g., {"Al": 2, "O": 3} for Al2O3

        Returns:
            Formula mass in amu

        Example:
            >>> cfml.formula_mass({"Fe": 2, "O": 3})
            159.687
        """
        total = 0.0
        for symbol, count in composition.items():
            mass = self.get_atomic_mass(symbol)
            total += mass * count
        return total

    # ==================== Math Functions ====================

    def factorial(self, n: int) -> int:
        """
        Calculate factorial of n.

        Uses CrysFML's optimized Fortran implementation.

        Args:
            n: Non-negative integer

        Returns:
            n! (n factorial)

        Example:
            >>> cfml.factorial(5)
            120
        """
        lib = self._get_module("cfml_maths")
        result = lib.factorial_i(n)
        return result.result

    def gcd(self, a: int, b: int) -> int:
        """
        Calculate greatest common divisor of a and b.

        Args:
            a, b: Integers

        Returns:
            GCD(a, b)

        Example:
            >>> cfml.gcd(48, 18)
            6
        """
        lib = self._get_module("cfml_maths")
        result = lib.gcd(a, b)
        return result.result

    def lcm(self, a: int, b: int) -> int:
        """
        Calculate least common multiple of a and b.

        Args:
            a, b: Positive integers

        Returns:
            LCM(a, b)

        Example:
            >>> cfml.lcm(12, 18)
            36
        """
        lib = self._get_module("cfml_maths")
        result = lib.lcm(a, b)
        return result.result

    # ==================== Module Information ====================

    def list_modules(self) -> Dict[str, int]:
        """
        List available CrysFML modules and their symbol counts.

        Returns:
            Dictionary mapping module names to number of available symbols

        Example:
            >>> for mod, count in cfml.list_modules().items():
            ...     print(f"{mod}: {count} symbols")
        """
        modules = {}
        module_names = [
            "cfml_globaldeps",
            "cfml_maths",
            "cfml_strings",
            "cfml_rational",
            "cfml_metrics",
            "cfml_gspacegroups",
            "cfml_symmetry_tables",
            "cfml_scattering_tables",
        ]

        for name in module_names:
            mod_path = self._mod_dir / f"{name}.mod"
            if mod_path.exists():
                try:
                    lib = self._get_module(name)
                    symbols = [s for s in dir(lib) if not s.startswith('_')]
                    modules[name] = len(symbols)
                except Exception:
                    modules[name] = -1  # Error loading

        return modules

    def get_version(self) -> str:
        """
        Get the version of the loaded CrysFML library.

        Returns:
            Version string or "unknown" if not available
        """
        try:
            lib = self._get_module("cfml_globaldeps")
            # Try to access version info if available
            if hasattr(lib, 'cfml_version'):
                return str(lib.cfml_version)
        except Exception:
            pass
        return "unknown"
