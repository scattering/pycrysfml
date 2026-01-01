"""
File I/O routines for PyCrysFML08

This module provides functions for reading crystallographic data files:
- CIF files (Crystallographic Information Framework)
- CFL files (CrysFML format with magnetic structure support)
- Diffraction data files
- Background files

The parsers are implemented in pure Python to avoid gfort2py limitations
with complex Fortran derived types.
"""

import re
import numpy as np
from typing import Tuple, List, Optional, Dict, Any
from pathlib import Path

from .structures import (
    CrystalCell, Atom, SpaceGroup, MagneticAtom,
    MagneticStructure, PropagationVector
)


def read_cif(filename: str) -> Tuple[CrystalCell, SpaceGroup, List[Atom]]:
    """
    Parse a CIF file and extract crystal structure information.

    Reads standard CIF format files and extracts:
    - Unit cell parameters
    - Space group information
    - Atom positions

    Args:
        filename: Path to CIF file

    Returns:
        Tuple of (CrystalCell, SpaceGroup, list of Atoms)

    Example:
        >>> cell, sg, atoms = read_cif("structure.cif")
        >>> print(f"Space group: {sg.symbol}")
        >>> print(f"Number of atoms: {len(atoms)}")
    """
    with open(filename, 'r') as f:
        content = f.read()

    # Parse cell parameters
    a = _parse_cif_value(content, '_cell_length_a')
    b = _parse_cif_value(content, '_cell_length_b')
    c = _parse_cif_value(content, '_cell_length_c')
    alpha = _parse_cif_value(content, '_cell_angle_alpha')
    beta = _parse_cif_value(content, '_cell_angle_beta')
    gamma = _parse_cif_value(content, '_cell_angle_gamma')

    cell = CrystalCell(a, b, c, alpha, beta, gamma)

    # Parse space group
    try:
        sg_number = int(_parse_cif_value(content, '_symmetry_Int_Tables_number'))
    except (ValueError, KeyError):
        sg_number = 0

    sg_symbol = _parse_cif_string(content, '_symmetry_space_group_name_H-M')
    if not sg_symbol:
        sg_symbol = _parse_cif_string(content, '_space_group_name_H-M_alt')

    spacegroup = SpaceGroup(sg_number, sg_symbol, crystal_system=cell.crystal_system())

    # Parse atoms
    atoms = _parse_cif_atoms(content)

    return cell, spacegroup, atoms


def read_cfl(filename: str) -> MagneticStructure:
    """
    Parse a CFL file (CrysFML format) with magnetic structure.

    CFL files are used by FullProf and CrysFML to describe crystal and
    magnetic structures. This parser extracts:
    - Unit cell parameters
    - Space group
    - Nuclear atoms
    - Magnetic atoms with moments and propagation vectors

    Args:
        filename: Path to CFL file

    Returns:
        MagneticStructure object containing all structure information

    Example:
        >>> mag_struct = read_cfl("magnetic.cfl")
        >>> for atom in mag_struct.magnetic_atoms:
        ...     print(f"{atom.label}: moment = {atom.moment_magnitude:.2f} μB")
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Initialize structure components
    cell = None
    spacegroup = None
    atoms = []
    magnetic_atoms = []
    kvectors = []
    magnetic_spacegroup = ""

    in_mag_structure = False
    current_matom = None

    for i, line in enumerate(lines):
        line = line.strip()

        # Skip comments and empty lines
        if not line or line.startswith('!'):
            continue

        # Parse TITLE (ignored)
        if line.upper().startswith('TITLE'):
            continue

        # Parse space group
        if line.upper().startswith('SPGR'):
            parts = line.split(None, 1)
            if len(parts) > 1:
                spacegroup = SpaceGroup(0, parts[1].strip())
            continue

        # Parse cell parameters
        if line.upper().startswith('CELL'):
            parts = line.split()
            if len(parts) >= 7:
                cell = CrystalCell(
                    a=float(parts[1]),
                    b=float(parts[2]),
                    c=float(parts[3]),
                    alpha=float(parts[4]),
                    beta=float(parts[5]),
                    gamma=float(parts[6])
                )
            continue

        # Parse nuclear atoms
        if line.upper().startswith('ATOM') and not in_mag_structure:
            parts = line.split()
            if len(parts) >= 7:
                atom = Atom(
                    label=parts[1],
                    symbol=parts[2],
                    x=float(parts[3]),
                    y=float(parts[4]),
                    z=float(parts[5]),
                    b_iso=float(parts[6]) if len(parts) > 6 else 0.0,
                    occupancy=float(parts[7]) if len(parts) > 7 else 1.0
                )
                atoms.append(atom)
            continue

        # Enter magnetic structure block
        if line.upper().startswith('MAG_STRUCTURE'):
            in_mag_structure = True
            continue

        # Exit magnetic structure block
        if line.upper().startswith('END_MAG_STRUCTURE'):
            if current_matom:
                magnetic_atoms.append(current_matom)
            in_mag_structure = False
            continue

        # Parse magnetic structure content
        if in_mag_structure:
            upper_line = line.upper()

            # Magnetic lattice
            if upper_line.startswith('LATTICE'):
                parts = line.split(None, 1)
                if len(parts) > 1:
                    magnetic_spacegroup = parts[1].strip()
                continue

            # Propagation vector
            if upper_line.startswith('KVECT'):
                parts = line.split()
                if len(parts) >= 4:
                    kvectors.append(PropagationVector(
                        kx=float(parts[1]),
                        ky=float(parts[2]),
                        kz=float(parts[3])
                    ))
                continue

            # Magnetic atom
            if upper_line.startswith('MATOM'):
                # Save previous magnetic atom if exists
                if current_matom:
                    magnetic_atoms.append(current_matom)

                parts = line.split()
                if len(parts) >= 7:
                    current_matom = MagneticAtom(
                        label=parts[1],
                        symbol=parts[1].rstrip('0123456789'),
                        form_factor=parts[2],
                        x=float(parts[3]),
                        y=float(parts[4]),
                        z=float(parts[5]),
                        b_iso=float(parts[6]) if len(parts) > 6 else 0.0,
                        occupancy=float(parts[7]) if len(parts) > 7 else 1.0
                    )
                continue

            # Basis vectors (BASR, BASI) for representation analysis
            if upper_line.startswith('BASR') or upper_line.startswith('BASI'):
                # These define the basis vectors for irreducible representations
                # For now, we just note they exist
                continue

            # Fourier coefficients (SKP lines)
            if upper_line.startswith('SKP') and current_matom:
                parts = line.split()
                # SKP format: SKP irep ik Rx Ry Rz Ix Iy Iz phase
                if len(parts) >= 9:
                    # Real and imaginary parts of Fourier coefficients
                    rx, ry, rz = float(parts[3]), float(parts[4]), float(parts[5])
                    ix, iy, iz = float(parts[6]), float(parts[7]), float(parts[8])
                    # For commensurate structures, moment is just the real part
                    current_matom.moment = [rx, ry, rz]
                    if len(parts) > 9:
                        current_matom.phase = float(parts[9])
                continue

            # BFCOEF (basis function coefficients)
            if upper_line.startswith('BFCOEF') and current_matom:
                parts = line.split()
                if len(parts) >= 5:
                    # bfcoef irep ibas c1 c2 phase
                    c1, c2 = float(parts[3]), float(parts[4])
                    current_matom.moment = [c1, c2, 0.0]  # Simplified
                continue

    # Create default cell if not found
    if cell is None:
        cell = CrystalCell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0)

    # Create default space group if not found
    if spacegroup is None:
        spacegroup = SpaceGroup(1, "P1")

    return MagneticStructure(
        cell=cell,
        spacegroup=spacegroup,
        magnetic_spacegroup=magnetic_spacegroup,
        kvectors=kvectors,
        atoms=atoms,
        magnetic_atoms=magnetic_atoms
    )


def read_diffraction_data(
    filename: str,
    skip_header: int = 1,
    columns: Tuple[int, int] = (0, 1)
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read diffraction pattern data from a file.

    Supports common formats with 2theta and intensity columns.

    Args:
        filename: Path to data file
        skip_header: Number of header lines to skip
        columns: Tuple of column indices for (2theta, intensity)

    Returns:
        Tuple of (2theta array, intensity array)

    Example:
        >>> tt, intensity = read_diffraction_data("pattern.dat")
        >>> print(f"2θ range: {tt.min():.1f} - {tt.max():.1f}°")
    """
    try:
        data = np.loadtxt(filename, skiprows=skip_header)
        if data.ndim == 1:
            # Single column - intensity only
            return None, data
        return data[:, columns[0]], data[:, columns[1]]
    except Exception as e:
        raise IOError(f"Error reading diffraction data from {filename}: {e}")


def read_background(
    filename: str,
    skip_header: int = 5
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read background spline data from a file.

    Background files typically contain (x, y) pairs defining
    spline points for background subtraction.

    Args:
        filename: Path to background file
        skip_header: Number of header lines to skip

    Returns:
        Tuple of (x positions, background values)

    Example:
        >>> bg_x, bg_y = read_background("background.bgr")
    """
    data = np.loadtxt(filename, skiprows=skip_header)
    return data[:, 0], data[:, 1]


def read_ill_data(
    filename: str,
    instrument: str = "",
    background_file: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Read diffraction data in ILL instrument format.

    Supports data from ILL neutron instruments like D1B, D2B, D20.

    Args:
        filename: Path to data file
        instrument: Instrument name (e.g., "D1B", "D2B")
        background_file: Optional path to background file

    Returns:
        Tuple of (2theta, intensity, error) arrays
    """
    # Simple implementation - read as standard format
    data = np.loadtxt(filename, skiprows=1)

    if data.ndim == 1:
        tt = None
        intensity = data
        error = None
    elif data.shape[1] >= 3:
        tt = data[:, 0]
        intensity = data[:, 1]
        error = data[:, 2]
    else:
        tt = data[:, 0]
        intensity = data[:, 1]
        error = None

    return tt, intensity, error


# ============== Private helper functions ==============

def _parse_cif_value(content: str, tag: str) -> float:
    """Parse a numeric value from CIF, handling uncertainties in parentheses."""
    pattern = rf'{tag}\s+([\d.]+)(?:\(\d+\))?'
    match = re.search(pattern, content)
    if match:
        return float(match.group(1))
    raise ValueError(f"CIF tag not found: {tag}")


def _parse_cif_string(content: str, tag: str) -> str:
    """Parse a string value from CIF."""
    pattern = rf"{tag}\s+['\"]?([^'\"\n]+)['\"]?"
    match = re.search(pattern, content)
    if match:
        return match.group(1).strip()
    return ""


def _parse_cif_atoms(content: str) -> List[Atom]:
    """Parse atom site loop from CIF."""
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
    def get_idx(name):
        return headers.index(name) if name in headers else None

    label_idx = get_idx('label')
    symbol_idx = get_idx('type_symbol')
    x_idx = get_idx('fract_x')
    y_idx = get_idx('fract_y')
    z_idx = get_idx('fract_z')
    occ_idx = get_idx('occupancy')
    mult_idx = get_idx('symmetry_multiplicity')
    b_idx = get_idx('B_iso_or_equiv')

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

            x = _parse_cif_coord(parts[x_idx]) if x_idx is not None else 0.0
            y = _parse_cif_coord(parts[y_idx]) if y_idx is not None else 0.0
            z = _parse_cif_coord(parts[z_idx]) if z_idx is not None else 0.0
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


def _parse_cif_coord(value: str) -> float:
    """Parse a coordinate value, handling uncertainties."""
    value = re.sub(r'\(\d+\)', '', value)
    return float(value)
