#!/usr/bin/env python3
"""
Element Properties Example

This example demonstrates using CrysFML08 scattering tables to look up
properties for various elements commonly used in crystallography.

Properties available:
- Atomic mass
- Covalent radius
- Neutron Fermi scattering length
- Atomic volume
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crysfml_api import CrysFML


def print_element_table(cfml, elements, title="Element Properties"):
    """Print a formatted table of element properties"""
    print(f"\n{title}")
    print("=" * 70)
    print(f"{'Element':<10} {'Mass (amu)':<12} {'Radius (A)':<12} {'Fermi (fm)':<12} {'Vol (A^3)':<12}")
    print("-" * 70)

    for elem in elements:
        try:
            mass = cfml.get_atomic_mass(elem)
            radius = cfml.get_covalent_radius(elem)
            fermi = cfml.get_fermi_length(elem)
            volume = cfml.get_atomic_volume(elem)
            print(f"{elem:<10} {mass:<12.4f} {radius:<12.4f} {fermi:<12.4f} {volume:<12.4f}")
        except Exception as e:
            print(f"{elem:<10} Error: {e}")

    print("=" * 70)


def main():
    """Main example routine"""

    print("=" * 70)
    print("CrysFML08 Element Properties Database")
    print("=" * 70)

    cfml = CrysFML()

    # Common elements in crystallography
    common_elements = ["H", "C", "N", "O", "F", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "K", "Ca"]
    print_element_table(cfml, common_elements, "Common Light Elements")

    # Transition metals
    transition_metals = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
    print_element_table(cfml, transition_metals, "First-Row Transition Metals")

    # Lanthanides (common in magnetic materials)
    lanthanides = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Yb"]
    print_element_table(cfml, lanthanides, "Lanthanides (Rare Earths)")

    # Heavy elements
    heavy_elements = ["Pb", "Bi", "Th", "U"]
    print_element_table(cfml, heavy_elements, "Heavy Elements")

    # Special cases - negative scattering lengths
    print("\n--- Elements with Negative Neutron Scattering Lengths ---")
    print("(Important for neutron contrast in crystallography)")
    negative_fermi = []
    test_elements = ["H", "Li", "Ti", "V", "Mn", "Ni"]
    for elem in test_elements:
        fermi = cfml.get_fermi_length(elem)
        if fermi < 0:
            negative_fermi.append((elem, fermi))
            print(f"  {elem}: b = {fermi:.4f} x 10^-12 cm")

    # Isotope considerations
    print("\n--- Hydrogen Isotopes Note ---")
    print("  H (protium):  b = -0.374 x 10^-12 cm (negative!)")
    print("  D (deuterium): b = +0.667 x 10^-12 cm (positive)")
    print("  Deuteration is commonly used for neutron contrast matching")

    # Calculate formula masses for common compounds
    print("\n--- Formula Masses of Common Compounds ---")
    compounds = [
        ("H2O", {"H": 2, "O": 1}),
        ("NaCl", {"Na": 1, "Cl": 1}),
        ("CaCO3", {"Ca": 1, "C": 1, "O": 3}),
        ("Fe2O3", {"Fe": 2, "O": 3}),
        ("Al2O3", {"Al": 2, "O": 3}),
        ("SiO2", {"Si": 1, "O": 2}),
        ("TiO2", {"Ti": 1, "O": 2}),
        ("BaTiO3", {"Ba": 1, "Ti": 1, "O": 3}),
        ("YBa2Cu3O7", {"Y": 1, "Ba": 2, "Cu": 3, "O": 7}),
    ]

    print(f"{'Compound':<15} {'Formula Mass (amu)':<20}")
    print("-" * 35)
    for name, composition in compounds:
        mass = cfml.formula_mass(composition)
        print(f"{name:<15} {mass:<20.4f}")

    # Math functions demo
    print("\n--- Mathematical Functions ---")
    print("GCD Examples:")
    pairs = [(12, 8), (48, 18), (100, 35), (13, 17)]
    for a, b in pairs:
        print(f"  gcd({a}, {b}) = {cfml.gcd(a, b)}")

    print("\nLCM Examples:")
    for a, b in pairs:
        print(f"  lcm({a}, {b}) = {cfml.lcm(a, b)}")

    print("\nFactorial Examples:")
    for n in [0, 1, 5, 10, 12]:
        print(f"  {n}! = {cfml.factorial(n)}")

    print("\n" + "=" * 70)
    print("Example completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
