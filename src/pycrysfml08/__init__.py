"""
PyCrysFML08 - Python bindings for CrysFML2008

A Python interface to the CrysFML2008 Fortran crystallographic library,
providing access to crystallographic calculations, space group operations,
and scattering data.

This package uses gfort2py to interface with a shared library compiled
from the CrysFML2008 Fortran source code.

CrysFML2008 is developed by:
    Juan Rodriguez-Carvajal (Institut Laue-Langevin)
    Javier Gonzalez-Platas (Universidad de La Laguna)

For more information about CrysFML, see:
    https://www.ill.eu/sites/fullprof/php/programs24b7.html

Example usage:
    >>> from pycrysfml08 import CrysFML
    >>> cfml = CrysFML()
    >>> mass = cfml.get_atomic_mass("Fe")
    >>> print(f"Iron atomic mass: {mass:.3f} amu")
    Iron atomic mass: 55.847 amu

License:
    This Python wrapper is provided under the same license terms as CrysFML.
    CrysFML is free software for academic use.
"""

__version__ = "0.1.0"
__author__ = "William Ratcliff"
__credits__ = [
    "Juan Rodriguez-Carvajal (CrysFML)",
    "Javier Gonzalez-Platas (CrysFML)",
    "Institut Laue-Langevin (ILL)",
]

from .api import CrysFML
from .structures import CrystalCell, Atom, SpaceGroup, MagneticAtom, MagneticStructure
from .io import read_cif, read_cfl, read_diffraction_data, read_background

__all__ = [
    "CrysFML",
    "CrystalCell",
    "Atom",
    "SpaceGroup",
    "MagneticAtom",
    "MagneticStructure",
    "read_cif",
    "read_cfl",
    "read_diffraction_data",
    "read_background",
]
