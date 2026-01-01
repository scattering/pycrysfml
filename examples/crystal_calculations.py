#!/usr/bin/env python3
"""
Crystal Calculations Example

This example demonstrates crystallographic calculations using CrysFML08 data:
- d-spacing calculations
- Bragg's law
- Structure factor considerations
- Unit cell geometry

These calculations combine Python-based crystallography with CrysFML's
accurate element data from the scattering tables.
"""

import os
import sys
import math
import numpy as np
from typing import Tuple, List

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crysfml_api import CrysFML, CrystalCell


def d_spacing_cubic(a: float, h: int, k: int, l: int) -> float:
    """Calculate d-spacing for cubic crystal"""
    return a / math.sqrt(h**2 + k**2 + l**2)


def d_spacing_hexagonal(a: float, c: float, h: int, k: int, l: int) -> float:
    """Calculate d-spacing for hexagonal crystal"""
    term1 = (4/3) * (h**2 + h*k + k**2) / a**2
    term2 = l**2 / c**2
    return 1.0 / math.sqrt(term1 + term2)


def d_spacing_orthorhombic(a: float, b: float, c: float, h: int, k: int, l: int) -> float:
    """Calculate d-spacing for orthorhombic crystal"""
    return 1.0 / math.sqrt((h/a)**2 + (k/b)**2 + (l/c)**2)


def d_spacing_general(cell: CrystalCell, h: int, k: int, l: int) -> float:
    """
    Calculate d-spacing for general triclinic crystal

    Uses the reciprocal metric tensor approach
    """
    a, b, c = cell.a, cell.b, cell.c
    alpha = math.radians(cell.alpha)
    beta = math.radians(cell.beta)
    gamma = math.radians(cell.gamma)

    # Reciprocal cell parameters
    cos_a, cos_b, cos_g = math.cos(alpha), math.cos(beta), math.cos(gamma)
    sin_a, sin_b, sin_g = math.sin(alpha), math.sin(beta), math.sin(gamma)

    V = cell.volume

    # Reciprocal cell
    a_star = b * c * sin_a / V
    b_star = a * c * sin_b / V
    c_star = a * b * sin_g / V

    cos_a_star = (cos_b * cos_g - cos_a) / (sin_b * sin_g)
    cos_b_star = (cos_a * cos_g - cos_b) / (sin_a * sin_g)
    cos_g_star = (cos_a * cos_b - cos_g) / (sin_a * sin_b)

    # d-spacing from reciprocal metric
    d_inv_sq = (h**2 * a_star**2 + k**2 * b_star**2 + l**2 * c_star**2
                + 2*h*k*a_star*b_star*cos_g_star
                + 2*k*l*b_star*c_star*cos_a_star
                + 2*h*l*a_star*c_star*cos_b_star)

    return 1.0 / math.sqrt(d_inv_sq)


def bragg_angle(d: float, wavelength: float) -> float:
    """
    Calculate Bragg angle (2theta) in degrees

    Args:
        d: d-spacing in Angstroms
        wavelength: X-ray/neutron wavelength in Angstroms

    Returns:
        2theta in degrees
    """
    sin_theta = wavelength / (2 * d)
    if abs(sin_theta) > 1:
        return float('nan')  # Reflection not accessible
    return 2 * math.degrees(math.asin(sin_theta))


def q_from_d(d: float) -> float:
    """Calculate Q = 2*pi/d from d-spacing"""
    return 2 * math.pi / d


def generate_hkl_list(h_max: int, k_max: int, l_max: int) -> List[Tuple[int, int, int]]:
    """Generate list of (h,k,l) indices up to given maxima"""
    hkl_list = []
    for h in range(-h_max, h_max + 1):
        for k in range(-k_max, k_max + 1):
            for l in range(-l_max, l_max + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                hkl_list.append((h, k, l))
    return hkl_list


def main():
    """Main example routine"""

    print("=" * 70)
    print("Crystallographic Calculations with CrysFML08")
    print("=" * 70)

    cfml = CrysFML()

    # Example 1: Al2O3 (Corundum) - Hexagonal/Rhombohedral
    print("\n--- Example 1: Al2O3 (Corundum, R-3c) ---")

    # From the CIF file
    al2o3_cell = CrystalCell(
        a=4.7698, b=4.7698, c=13.0243,
        alpha=90.0, beta=90.0, gamma=120.0
    )

    print(f"Cell: a={al2o3_cell.a:.4f}, c={al2o3_cell.c:.4f} A")
    print(f"Volume: {al2o3_cell.volume:.2f} A^3")

    # Cu K-alpha wavelength
    wavelength = 1.5406  # Cu K-alpha

    print(f"\nReflections for Cu K-alpha (lambda = {wavelength} A):")
    print(f"{'(h k l)':<12} {'d (A)':<10} {'2theta (deg)':<12} {'Q (A^-1)':<10}")
    print("-" * 50)

    # Some important reflections for corundum
    reflections = [(0, 1, 2), (1, 0, 4), (1, 1, 0), (1, 1, 3), (0, 2, 4), (1, 1, 6)]
    for hkl in reflections:
        d = d_spacing_hexagonal(al2o3_cell.a, al2o3_cell.c, *hkl)
        two_theta = bragg_angle(d, wavelength)
        q = q_from_d(d)
        print(f"({hkl[0]:2d} {hkl[1]:2d} {hkl[2]:2d})   {d:8.4f}   {two_theta:10.3f}   {q:8.4f}")

    # Example 2: NaCl (Rock Salt) - Cubic
    print("\n--- Example 2: NaCl (Rock Salt, Fm-3m) ---")

    a_nacl = 5.6402  # Angstroms

    print(f"Cell: a = {a_nacl:.4f} A (cubic)")
    print(f"Volume: {a_nacl**3:.2f} A^3")

    # Calculate density
    formula_mass = cfml.formula_mass({"Na": 1, "Cl": 1})
    Z = 4  # Formula units per cell
    avogadro = 6.02214076e23
    volume_cm3 = (a_nacl * 1e-8)**3
    density = (Z * formula_mass) / (avogadro * volume_cm3)
    print(f"Formula mass: {formula_mass:.4f} amu")
    print(f"Calculated density: {density:.3f} g/cm^3 (literature: 2.17 g/cm^3)")

    print(f"\nReflections for Cu K-alpha:")
    print(f"{'(h k l)':<12} {'d (A)':<10} {'2theta (deg)':<12}")
    print("-" * 40)

    # NaCl reflections (FCC, F-centered: h,k,l all even or all odd)
    nacl_reflections = [(1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1), (2, 2, 2), (4, 0, 0)]
    for hkl in nacl_reflections:
        d = d_spacing_cubic(a_nacl, *hkl)
        two_theta = bragg_angle(d, wavelength)
        if not math.isnan(two_theta):
            print(f"({hkl[0]:2d} {hkl[1]:2d} {hkl[2]:2d})   {d:8.4f}   {two_theta:10.3f}")

    # Example 3: Scattering contrast
    print("\n--- Example 3: Neutron Scattering Contrast ---")
    print("Comparing X-ray and neutron sensitivity:")

    elements = ["H", "D", "C", "N", "O", "Fe", "Pb"]
    print(f"{'Element':<10} {'Mass (amu)':<12} {'b_neutron':<12} {'Contrast':<15}")
    print("-" * 55)

    # Note: D (deuterium) isn't in the standard tables, showing H twice
    for elem in ["H", "C", "N", "O", "Fe", "Pb"]:
        mass = cfml.get_atomic_mass(elem)
        fermi = cfml.get_fermi_length(elem)
        # X-ray scattering roughly proportional to Z (number of electrons)
        Z_approx = round(mass / 2) if mass < 40 else round(mass / 2.5)
        contrast = "Light element" if mass < 20 else ("Medium" if mass < 60 else "Heavy")
        print(f"{elem:<10} {mass:<12.2f} {fermi:<12.4f} {contrast:<15}")

    print("\nNote: Neutrons can 'see' light elements like H clearly,")
    print("while X-rays are dominated by heavy elements like Pb.")

    # Example 4: Reciprocal lattice
    print("\n--- Example 4: Reciprocal Lattice Parameters ---")

    cells = [
        ("Cubic (NaCl)", CrystalCell(5.64, 5.64, 5.64, 90, 90, 90)),
        ("Hexagonal (Al2O3)", CrystalCell(4.77, 4.77, 13.02, 90, 90, 120)),
        ("Orthorhombic", CrystalCell(4.0, 5.0, 6.0, 90, 90, 90)),
    ]

    for name, cell in cells:
        a, b, c = cell.a, cell.b, cell.c
        V = cell.volume

        # For orthogonal cells, a* = 2*pi/a
        # This is simplified for orthogonal cases
        if cell.alpha == 90 and cell.beta == 90 and cell.gamma == 90:
            a_star = 2 * math.pi / a
            b_star = 2 * math.pi / b
            c_star = 2 * math.pi / c
            print(f"\n{name}:")
            print(f"  Direct: a={a:.3f}, b={b:.3f}, c={c:.3f} A")
            print(f"  Reciprocal: a*={a_star:.4f}, b*={b_star:.4f}, c*={c_star:.4f} A^-1")
        elif cell.gamma == 120:
            # Hexagonal
            a_star = 4 * math.pi / (math.sqrt(3) * a)
            c_star = 2 * math.pi / c
            print(f"\n{name}:")
            print(f"  Direct: a={a:.3f}, c={c:.3f} A, gamma=120 deg")
            print(f"  Reciprocal: a*={a_star:.4f}, c*={c_star:.4f} A^-1")

    print("\n" + "=" * 70)
    print("Calculations completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
