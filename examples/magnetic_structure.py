#!/usr/bin/env python3
"""
Magnetic Structure Example

This example demonstrates reading and analyzing magnetic structures
using PyCrysFML. It shows how to:

1. Read CFL files containing magnetic structure information
2. Parse propagation vectors and magnetic atoms
3. Calculate magnetic moments and analyze magnetic symmetry
4. Compare commensurate vs incommensurate structures

Examples based on:
- HoBaNiO (HOBK) - Antiferromagnetic with k=(1/2, 0, 1/2)
- DyMn6Ge6 - Complex magnetic structure with multiple k-vectors

Reference data from the original pycrysfml examples.
"""

import os
import sys
import math
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pycrysfml import CrysFML, read_cfl
from pycrysfml.structures import PropagationVector

# Data directory
DATA_DIR = Path(__file__).parent / "data"


def analyze_magnetic_structure(cfl_file: str, name: str):
    """Analyze a magnetic structure from a CFL file."""

    print(f"\n{'='*60}")
    print(f"Magnetic Structure: {name}")
    print(f"{'='*60}")

    # Read the CFL file
    mag_struct = read_cfl(cfl_file)

    # Print cell parameters
    cell = mag_struct.cell
    print(f"\n--- Unit Cell ---")
    print(f"  a = {cell.a:.4f} Å")
    print(f"  b = {cell.b:.4f} Å")
    print(f"  c = {cell.c:.4f} Å")
    print(f"  α = {cell.alpha:.2f}°, β = {cell.beta:.2f}°, γ = {cell.gamma:.2f}°")
    print(f"  Volume = {cell.volume:.2f} Å³")
    print(f"  Crystal system: {cell.crystal_system()}")

    # Print space group
    print(f"\n--- Space Group ---")
    print(f"  Symbol: {mag_struct.spacegroup.symbol}")
    if mag_struct.magnetic_spacegroup:
        print(f"  Magnetic space group: {mag_struct.magnetic_spacegroup}")

    # Print propagation vectors
    print(f"\n--- Propagation Vector(s) ---")
    if mag_struct.kvectors:
        for i, k in enumerate(mag_struct.kvectors, 1):
            print(f"  k{i} = ({k.kx:.4f}, {k.ky:.4f}, {k.kz:.4f})")
            if k.is_commensurate:
                print(f"       [Commensurate]")
            else:
                print(f"       [Incommensurate]")
    else:
        print("  k = (0, 0, 0) [Ferromagnetic or not specified]")

    # Print nuclear atoms
    print(f"\n--- Nuclear Atoms ({len(mag_struct.atoms)}) ---")
    if mag_struct.atoms:
        print(f"  {'Label':<8} {'Symbol':<6} {'x':>8} {'y':>8} {'z':>8}")
        print(f"  {'-'*42}")
        for atom in mag_struct.atoms:
            print(f"  {atom.label:<8} {atom.symbol:<6} {atom.x:8.4f} {atom.y:8.4f} {atom.z:8.4f}")

    # Print magnetic atoms
    print(f"\n--- Magnetic Atoms ({len(mag_struct.magnetic_atoms)}) ---")
    if mag_struct.magnetic_atoms:
        print(f"  {'Label':<8} {'Form Factor':<10} {'|M| (μB)':>10} {'Mx':>8} {'My':>8} {'Mz':>8}")
        print(f"  {'-'*58}")
        for atom in mag_struct.magnetic_atoms:
            m = atom.moment
            mag = atom.moment_magnitude
            print(f"  {atom.label:<8} {atom.form_factor:<10} {mag:10.3f} {m[0]:8.3f} {m[1]:8.3f} {m[2]:8.3f}")

        # Calculate total moment
        total_moment = sum(a.moment_magnitude for a in mag_struct.magnetic_atoms)
        print(f"\n  Total magnetic moment: {total_moment:.3f} μB/formula unit")

    # Magnetic structure classification
    print(f"\n--- Structure Classification ---")
    if mag_struct.is_commensurate:
        print("  Type: Commensurate magnetic structure")

        # Check for special cases
        if not mag_struct.kvectors or all(k.kx == 0 and k.ky == 0 and k.kz == 0 for k in mag_struct.kvectors):
            print("  Order: Ferromagnetic (k = 0)")
        elif any(k.kx == 0.5 or k.ky == 0.5 or k.kz == 0.5 for k in mag_struct.kvectors):
            print("  Order: Antiferromagnetic (unit cell doubling)")
    else:
        print("  Type: Incommensurate magnetic structure")
        print("  (Magnetic unit cell is not a simple multiple of nuclear cell)")


def compare_neutron_sensitivity(cfml: CrysFML, elements: list):
    """Compare neutron scattering sensitivity for elements."""

    print(f"\n{'='*60}")
    print("Neutron Scattering Properties for Magnetic Elements")
    print(f"{'='*60}")

    print(f"\n{'Element':<10} {'Mass (amu)':<12} {'b (fm)':<10} {'σ_abs (b)':<12}")
    print("-" * 50)

    for elem in elements:
        try:
            mass = cfml.get_atomic_mass(elem)
            fermi = cfml.get_fermi_length(elem)
            try:
                abs_xs = cfml.get_absorption_xs(elem)
            except:
                abs_xs = float('nan')

            print(f"{elem:<10} {mass:<12.3f} {fermi:<10.4f} {abs_xs:<12.3f}")
        except Exception as e:
            print(f"{elem:<10} Error: {e}")


def main():
    """Main example routine."""

    print("=" * 60)
    print("Magnetic Structure Analysis with PyCrysFML")
    print("=" * 60)

    # Initialize CrysFML
    cfml = CrysFML()

    # Example 1: HOBK - HoBaNiO antiferromagnet
    hobk_file = DATA_DIR / "hobk" / "hobk1.cfl"
    if hobk_file.exists():
        analyze_magnetic_structure(str(hobk_file), "HoBaNiO (HOBK)")
    else:
        print(f"\nNote: HOBK data file not found at {hobk_file}")
        print("Creating example magnetic structure...")

        # Create example structure programmatically
        from pycrysfml.structures import (
            CrystalCell, SpaceGroup, MagneticAtom, MagneticStructure
        )

        cell = CrystalCell(a=3.754, b=5.730, c=11.269, alpha=90, beta=90, gamma=90)
        sg = SpaceGroup(71, "Immm")

        ho = MagneticAtom(
            label="Ho", symbol="Ho", form_factor="JHO3",
            x=0.5, y=0.0, z=0.2026,
            moment=[0.127, 8.993, 0.0]
        )
        ni = MagneticAtom(
            label="Ni", symbol="Ni", form_factor="MNI2",
            x=0.0, y=0.0, z=0.0,
            moment=[0.584, -1.285, 0.0]
        )

        mag_struct = MagneticStructure(
            cell=cell,
            spacegroup=sg,
            magnetic_spacegroup="I -1",
            kvectors=[PropagationVector(0.5, 0.0, 0.5)],
            magnetic_atoms=[ho, ni]
        )

        print(f"\n--- Example: HOBK-type Structure ---")
        print(f"Cell: {cell.a:.3f} x {cell.b:.3f} x {cell.c:.3f} Å")
        print(f"Space group: {sg.symbol}")
        print(f"k-vector: (0.5, 0, 0.5) - antiferromagnetic")
        print(f"Ho moment: {ho.moment_magnitude:.2f} μB")
        print(f"Ni moment: {ni.moment_magnitude:.2f} μB")

    # Example 2: Dy compound
    dy_file = DATA_DIR / "dy" / "dy.cfl"
    if dy_file.exists():
        analyze_magnetic_structure(str(dy_file), "DyMn6Ge6")

    # Compare neutron sensitivity for magnetic elements
    magnetic_elements = ["Ho", "Dy", "Gd", "Fe", "Mn", "Ni", "Co"]
    compare_neutron_sensitivity(cfml, magnetic_elements)

    # Discussion of magnetic neutron scattering
    print(f"\n{'='*60}")
    print("Notes on Magnetic Neutron Scattering")
    print(f"{'='*60}")
    print("""
Neutron scattering is uniquely suited for studying magnetic structures because:

1. Neutrons have a magnetic moment that couples to unpaired electrons
2. Magnetic scattering cross section is comparable to nuclear scattering
3. Both ordered and fluctuating magnetic moments can be studied

The magnetic form factor f(Q) describes how the magnetic scattering
intensity decreases with Q = 4π sin(θ)/λ due to the spatial extent
of the unpaired electron cloud.

For rare earth elements (Ho, Dy, Gd), the 4f electrons are well-localized,
giving a slowly-decaying form factor.

For transition metals (Fe, Mn, Ni), the 3d electrons are more extended,
resulting in faster form factor decay.
""")

    print("=" * 60)
    print("Example completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
