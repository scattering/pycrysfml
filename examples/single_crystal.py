#!/usr/bin/env python3
"""
Single Crystal Diffraction Example

This example demonstrates single crystal diffraction concepts
and calculations using PyCrysFML:

1. Structure factor considerations
2. Extinction corrections
3. hkl reflection generation
4. Intensity calculations

Based on single crystal data from the original pycrysfml examples.
"""

import os
import sys
import math
import numpy as np
from pathlib import Path
from typing import List, Tuple

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pycrysfml import CrysFML, read_cfl
from pycrysfml.structures import CrystalCell, Atom

# Data directory
DATA_DIR = Path(__file__).parent / "data"


def generate_hkl_list(h_max: int, k_max: int, l_max: int,
                      cell: CrystalCell,
                      d_min: float = 0.5) -> List[Tuple[int, int, int, float]]:
    """
    Generate list of (h, k, l) reflections up to given limits.

    Args:
        h_max, k_max, l_max: Maximum Miller indices
        cell: Crystal unit cell
        d_min: Minimum d-spacing to include

    Returns:
        List of (h, k, l, d) tuples sorted by d-spacing
    """
    reflections = []

    for h in range(-h_max, h_max + 1):
        for k in range(-k_max, k_max + 1):
            for l in range(-l_max, l_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue

                d = cell.d_spacing(h, k, l)
                if d >= d_min:
                    reflections.append((h, k, l, d))

    # Sort by d-spacing (descending)
    reflections.sort(key=lambda x: -x[3])

    return reflections


def sin_theta_over_lambda(d: float) -> float:
    """Calculate sin(θ)/λ from d-spacing."""
    return 1.0 / (2.0 * d)


def lorentz_factor(two_theta: float) -> float:
    """
    Calculate Lorentz factor for single crystal.

    L = 1 / sin(2θ)

    This accounts for the different times reflections spend
    in the diffraction condition.
    """
    two_theta_rad = math.radians(two_theta)
    return 1.0 / math.sin(two_theta_rad)


def polarization_factor(two_theta: float, polarization: float = 0.5) -> float:
    """
    Calculate polarization factor.

    For unpolarized beam (polarization = 0.5):
    P = (1 + cos²(2θ)) / 2

    Args:
        two_theta: Scattering angle in degrees
        polarization: Degree of polarization (0 to 1)

    Returns:
        Polarization factor
    """
    two_theta_rad = math.radians(two_theta)
    cos2 = math.cos(two_theta_rad) ** 2
    return (1.0 + cos2) / 2.0


def extinction_correction(intensity: float, extinction_param: float,
                         wavelength: float, volume: float) -> float:
    """
    Apply secondary extinction correction (Becker-Coppens Type I).

    The extinction reduces the observed intensity of strong reflections
    due to multiple scattering within the crystal.

    Args:
        intensity: Calculated intensity
        extinction_param: Extinction parameter (0 = no extinction)
        wavelength: Wavelength in Angstroms
        volume: Unit cell volume in Å³

    Returns:
        Extinction-corrected intensity
    """
    if extinction_param <= 0:
        return intensity

    # Simplified Type I isotropic extinction
    x = extinction_param * intensity * wavelength ** 3 / volume
    y = 1.0 / math.sqrt(1.0 + x)
    return intensity * y


def calculate_structure_factor_example(cell: CrystalCell, atoms: List[Atom],
                                       h: int, k: int, l: int,
                                       cfml: CrysFML) -> complex:
    """
    Calculate structure factor F(hkl) for a reflection.

    F(hkl) = Σ fⱼ · exp(2πi(hxⱼ + kyⱼ + lzⱼ)) · exp(-Bⱼ·sin²θ/λ²)

    This is a simplified calculation for demonstration.
    A full implementation would use CrysFML's structure factor routines.
    """
    d = cell.d_spacing(h, k, l)
    stol = sin_theta_over_lambda(d)
    stol_sq = stol ** 2

    F = complex(0, 0)

    for atom in atoms:
        # Get approximate scattering factor (proportional to Z)
        # In a full implementation, this would use X-ray form factor tables
        mass = cfml.get_atomic_mass(atom.symbol)
        Z_approx = int(round(mass / 2))  # Very rough approximation

        # Phase factor
        phase = 2 * math.pi * (h * atom.x + k * atom.y + l * atom.z)

        # Temperature factor
        temp_factor = math.exp(-atom.b_iso * stol_sq)

        # Add contribution
        F += Z_approx * atom.occupancy * temp_factor * complex(
            math.cos(phase), math.sin(phase)
        )

    return F


def main():
    """Main example routine."""

    print("=" * 70)
    print("Single Crystal Diffraction Analysis with PyCrysFML")
    print("=" * 70)

    # Initialize CrysFML
    cfml = CrysFML()

    # Example: Create a simple structure (NaCl type)
    print("\n--- Example: NaCl Structure ---")

    cell = CrystalCell(a=5.6402, b=5.6402, c=5.6402, alpha=90, beta=90, gamma=90)
    print(f"Unit cell: a = {cell.a:.4f} Å (cubic)")
    print(f"Volume: {cell.volume:.2f} Å³")
    print(f"Crystal system: {cell.crystal_system()}")

    # Define atoms (simplified NaCl)
    atoms = [
        Atom("Na1", "Na", 0.0, 0.0, 0.0, occupancy=1.0, b_iso=1.0),
        Atom("Cl1", "Cl", 0.5, 0.5, 0.5, occupancy=1.0, b_iso=1.0),
    ]

    print(f"\nAtoms:")
    for atom in atoms:
        mass = cfml.get_atomic_mass(atom.symbol)
        print(f"  {atom.label}: {atom.symbol} at ({atom.x:.2f}, {atom.y:.2f}, {atom.z:.2f}), "
              f"mass = {mass:.2f} amu")

    # Generate reflections
    print("\n--- Reflection List ---")
    reflections = generate_hkl_list(5, 5, 5, cell, d_min=1.0)

    # Filter for allowed reflections (FCC: h,k,l all even or all odd)
    allowed = []
    for h, k, l, d in reflections:
        parity = (h % 2, k % 2, l % 2)
        if parity == (0, 0, 0) or parity == (1, 1, 1):
            allowed.append((h, k, l, d))

    print(f"\nAllowed reflections (FCC): {len(allowed)}")
    print(f"\n{'(h k l)':<12} {'d (Å)':<10} {'sin(θ)/λ':<12} {'|F|²':>10}")
    print("-" * 50)

    for h, k, l, d in allowed[:15]:  # First 15 reflections
        stol = sin_theta_over_lambda(d)

        # Calculate structure factor
        F = calculate_structure_factor_example(cell, atoms, h, k, l, cfml)
        F_sq = abs(F) ** 2

        print(f"({h:2d} {k:2d} {l:2d})   {d:8.4f}   {stol:10.4f}   {F_sq:10.1f}")

    # Extinction effects demonstration
    print("\n--- Extinction Effects ---")
    print("Strong reflections are reduced by extinction (multiple scattering)")
    print("\nExtinction correction for (2 2 0):")

    h, k, l = 2, 2, 0
    d = cell.d_spacing(h, k, l)
    F = calculate_structure_factor_example(cell, atoms, h, k, l, cfml)
    I_calc = abs(F) ** 2

    print(f"Calculated |F|² = {I_calc:.1f}")
    for ext_param in [0.0, 0.001, 0.01, 0.1]:
        I_corr = extinction_correction(I_calc, ext_param, wavelength=0.71, volume=cell.volume)
        ratio = I_corr / I_calc if I_calc > 0 else 1.0
        print(f"  ext = {ext_param:.3f}: I_obs/I_calc = {ratio:.3f}")

    # Lorentz-polarization corrections
    print("\n--- Lorentz-Polarization Factors ---")
    print("These geometric corrections depend on scattering angle.")
    print(f"\n{'2θ (deg)':<10} {'Lorentz':<12} {'Polarization':<15} {'LP':<12}")
    print("-" * 55)

    for two_theta in [10, 30, 60, 90, 120, 150]:
        L = lorentz_factor(two_theta)
        P = polarization_factor(two_theta)
        LP = L * P
        print(f"{two_theta:<10} {L:<12.4f} {P:<15.4f} {LP:<12.4f}")

    # Read single crystal data if available
    sxtal_file = DATA_DIR / "sxtal" / "35Lcombined.cfl"
    if sxtal_file.exists():
        print(f"\n--- Single Crystal Data File ---")
        print(f"Reading: {sxtal_file}")

        mag_struct = read_cfl(str(sxtal_file))
        print(f"Cell: {mag_struct.cell.a:.3f} x {mag_struct.cell.b:.3f} x {mag_struct.cell.c:.3f} Å")
        print(f"Space group: {mag_struct.spacegroup.symbol}")
        print(f"Magnetic atoms: {len(mag_struct.magnetic_atoms)}")

        if mag_struct.kvectors:
            print("Propagation vectors:")
            for k in mag_struct.kvectors:
                print(f"  k = ({k.kx:.4f}, {k.ky:.4f}, {k.kz:.4f})")

    # Summary
    print("\n" + "=" * 70)
    print("Summary: Single Crystal vs Powder Diffraction")
    print("=" * 70)
    print("""
Single crystal diffraction provides:
- Complete 3D reciprocal space mapping
- Separation of overlapping reflections
- Precise intensity measurements for structure determination
- Access to magnetic satellite reflections

Key corrections needed:
1. Lorentz factor - accounts for reflection geometry
2. Polarization factor - for X-rays, depends on polarization state
3. Absorption correction - for sample absorption
4. Extinction correction - for strong reflections (Type I and II)
5. Thermal diffuse scattering - background from phonons

For magnetic structures:
- Magnetic form factor f(Q) must be included
- Polarized neutrons can separate nuclear and magnetic contributions
- Spherical neutron polarimetry provides moment direction
""")

    print("=" * 70)
    print("Example completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
