#!/usr/bin/env python3
"""
Al2O3 (Corundum) Example

This example demonstrates using CrysFML08 Python bindings to:
1. Read a CIF file and extract crystal structure information
2. Look up scattering data for the elements
3. Calculate basic crystallographic properties

Based on the original pycrysfml Al2O3 example, adapted for gfort2py bindings.

Reference:
    Kondo, S.; Tateishi, K.; Ishizawa, N. (2008)
    Japanese Journal of Applied Physics 47, 616-619
    "Structural evolution of corundum at high temperatures"
"""

import os
import sys
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crysfml_api import CrysFML, read_diffraction_data, read_background

# Data files
DATAPATH = os.path.dirname(os.path.abspath(__file__))
CIF_FILE = os.path.join(DATAPATH, "Al2O3.cif")
DATA_FILE = os.path.join(DATAPATH, "Al2O3.dat")
BACKGROUND_FILE = os.path.join(DATAPATH, "Al2O3 Background.BGR")


def main():
    """Main example routine"""

    print("=" * 60)
    print("Al2O3 (Corundum) Crystal Structure Analysis")
    print("=" * 60)

    # Initialize CrysFML interface
    cfml = CrysFML()

    # Read CIF file
    print("\nReading CIF file:", CIF_FILE)
    cell, spacegroup, atoms = cfml.read_cif(CIF_FILE)

    # Print cell parameters
    print("\n--- Unit Cell Parameters ---")
    print(f"  a = {cell.a:.4f} A")
    print(f"  b = {cell.b:.4f} A")
    print(f"  c = {cell.c:.4f} A")
    print(f"  alpha = {cell.alpha:.2f} deg")
    print(f"  beta  = {cell.beta:.2f} deg")
    print(f"  gamma = {cell.gamma:.2f} deg")
    print(f"  Volume = {cell.volume:.2f} A^3")

    # Print space group
    print("\n--- Space Group ---")
    print(f"  Number: {spacegroup.number}")
    print(f"  Symbol: {spacegroup.symbol}")

    # Print atoms
    print("\n--- Atoms ---")
    print(f"  {'Label':<8} {'Symbol':<6} {'x':>8} {'y':>8} {'z':>8} {'Occ':>6} {'Mult':>5}")
    print("  " + "-" * 52)
    for atom in atoms:
        print(f"  {atom.label:<8} {atom.symbol:<6} {atom.x:8.4f} {atom.y:8.4f} {atom.z:8.4f} {atom.occupancy:6.2f} {atom.multiplicity:5d}")

    # Look up element properties using CrysFML scattering tables
    print("\n--- Element Properties (from CrysFML scattering tables) ---")
    elements = set(atom.symbol for atom in atoms)
    for elem in sorted(elements):
        mass = cfml.get_atomic_mass(elem)
        radius = cfml.get_covalent_radius(elem)
        fermi = cfml.get_fermi_length(elem)
        print(f"  {elem}:")
        print(f"    Atomic mass:      {mass:.4f} amu")
        print(f"    Covalent radius:  {radius:.4f} A")
        print(f"    Fermi length:     {fermi:.4f} x 10^-12 cm")

    # Calculate formula mass
    print("\n--- Formula Mass ---")
    composition = {"Al": 2, "O": 3}
    formula_mass = cfml.formula_mass(composition)
    print(f"  Al2O3: {formula_mass:.4f} amu")

    # Calculate density
    Z = 6  # Formula units per unit cell (from CIF)
    avogadro = 6.02214076e23
    # Volume in cm^3 (1 A = 1e-8 cm)
    volume_cm3 = cell.volume * 1e-24
    density = (Z * formula_mass) / (avogadro * volume_cm3)
    print(f"  Calculated density: {density:.3f} g/cm^3")
    print(f"  (Literature value: ~3.98 g/cm^3)")

    # Read and display diffraction data info
    print("\n--- Diffraction Data ---")
    if os.path.exists(DATA_FILE):
        tt, intensity = read_diffraction_data(DATA_FILE)
        if tt is not None:
            print(f"  Data file: {DATA_FILE}")
            print(f"  2theta range: {tt.min():.2f} - {tt.max():.2f} degrees")
            print(f"  Number of points: {len(tt)}")
            print(f"  Max intensity: {intensity.max():.1f}")
        else:
            print(f"  Intensity data loaded ({len(intensity)} points)")

    # Read background
    if os.path.exists(BACKGROUND_FILE):
        bg_x, bg_y = read_background(BACKGROUND_FILE)
        print(f"\n  Background file: {BACKGROUND_FILE}")
        print(f"  Background spline points: {len(bg_x)}")

    # Demonstrate math functions
    print("\n--- CrysFML Math Functions Demo ---")
    print(f"  factorial(5) = {cfml.factorial(5)}")
    print(f"  factorial(10) = {cfml.factorial(10)}")
    print(f"  gcd(48, 18) = {cfml.gcd(48, 18)}")
    print(f"  lcm(12, 18) = {cfml.lcm(12, 18)}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
