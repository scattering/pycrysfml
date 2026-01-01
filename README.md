# PyCrysFML08

Python bindings for the CrysFML2008 crystallographic Fortran library.

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

## Overview

PyCrysFML08 provides Python access to the powerful CrysFML2008 Fortran library for crystallographic computing. This package enables:

- **Scattering data lookups**: Atomic masses, neutron scattering lengths, X-ray form factors
- **Crystal structure handling**: CIF and CFL file parsing, unit cell calculations
- **Magnetic structures**: Propagation vectors, magnetic moments, magnetic symmetry
- **Mathematical utilities**: Crystallographic calculations, space group operations

## Credits and Acknowledgments

### CrysFML2008

This package provides Python bindings to **CrysFML2008**, a comprehensive Fortran library for crystallographic computing developed at the Institut Laue-Langevin (ILL):

**Authors:**
- **Juan Rodriguez-Carvajal** - Institut Laue-Langevin, Grenoble, France
- **Javier Gonzalez-Platas** - Universidad de La Laguna, Tenerife, Spain

**CrysFML Homepage:** https://www.ill.eu/sites/fullprof/

**CrysFML Repository:** https://code.ill.fr/scientific-software/CrysFML2008

CrysFML is distributed under the LGPL license and is free for academic use.

### Original PyCrysFML

This project builds upon the original **pycrysfml** Python bindings:
- https://github.com/scattering/pycrysfml

### ILL PyCrysFML

The ILL has also developed Python bindings as part of their scientific software suite:
- https://code.ill.fr/scientific-software/pycrysfml

## Installation

### From PyPI (when available)

```bash
pip install pycrysfml08
```

### From Source

#### Prerequisites

1. **gfortran 14** (gfortran 15 module format not yet supported by gfort2py)
   ```bash
   # macOS
   brew install gcc@14

   # Ubuntu/Debian
   sudo apt install gfortran-14
   ```

2. **Python 3.9+** with matching architecture (ARM64 or x86_64)

#### Build Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/scattering/pycrysfml08.git
   cd pycrysfml08
   ```

2. Clone CrysFML2008 source:
   ```bash
   git clone https://code.ill.fr/scientific-software/CrysFML2008.git
   ```

3. Build the Fortran library (see [BUILDING.md](BUILDING.md) for detailed instructions):
   ```bash
   ./scripts/build_library.sh
   ```

4. Install the Python package:
   ```bash
   pip install -e .
   ```

## Quick Start

```python
from pycrysfml08 import CrysFML, read_cif, read_cfl

# Initialize the library
cfml = CrysFML()

# Look up element properties
mass = cfml.get_atomic_mass("Fe")
print(f"Iron atomic mass: {mass:.3f} amu")

b = cfml.get_fermi_length("Fe")
print(f"Iron neutron scattering length: {b:.4f} fm")

# Calculate formula mass
formula_mass = cfml.formula_mass({"Fe": 2, "O": 3})
print(f"Fe2O3 formula mass: {formula_mass:.3f} amu")

# Read a crystal structure
cell, spacegroup, atoms = read_cif("structure.cif")
print(f"Space group: {spacegroup.symbol}")
print(f"Cell volume: {cell.volume:.2f} Å³")

# Read a magnetic structure
mag_struct = read_cfl("magnetic.cfl")
for atom in mag_struct.magnetic_atoms:
    print(f"{atom.label}: moment = {atom.moment_magnitude:.2f} μB")
```

## Examples

See the `examples/` directory for complete examples:

- `al2o3_example.py` - Crystal structure analysis (Al₂O₃ corundum)
- `element_properties.py` - Element database lookups
- `crystal_calculations.py` - d-spacing, Bragg angles, density calculations
- `magnetic_structure.py` - Magnetic structure analysis (HOBK, Dy compounds)
- `single_crystal.py` - Single crystal diffraction

Run an example:
```bash
python examples/al2o3_example.py
```

## API Reference

### CrysFML Class

The main interface to the Fortran library:

```python
from pycrysfml08 import CrysFML

cfml = CrysFML()

# Scattering tables
cfml.get_atomic_mass(symbol)      # Atomic mass in amu
cfml.get_covalent_radius(symbol)  # Covalent radius in Å
cfml.get_fermi_length(symbol)     # Neutron scattering length
cfml.get_atomic_volume(symbol)    # Atomic volume in Å³
cfml.get_absorption_xs(symbol)    # Neutron absorption cross section
cfml.get_incoherent_xs(symbol)    # Incoherent scattering cross section
cfml.formula_mass(composition)    # Calculate formula mass

# Math functions
cfml.factorial(n)    # n!
cfml.gcd(a, b)       # Greatest common divisor
cfml.lcm(a, b)       # Least common multiple
```

### Data Structures

```python
from pycrysfml08 import CrystalCell, Atom, SpaceGroup, MagneticAtom

# Crystal cell
cell = CrystalCell(a=5.0, b=5.0, c=5.0, alpha=90, beta=90, gamma=90)
print(cell.volume)          # Cell volume
print(cell.crystal_system)  # 'cubic', 'hexagonal', etc.
print(cell.d_spacing(1,1,1)) # d-spacing for (111)

# Atoms
atom = Atom(label="Fe1", symbol="Fe", x=0.0, y=0.0, z=0.0)

# Magnetic atoms
mag_atom = MagneticAtom(
    label="Fe1", symbol="Fe", x=0, y=0, z=0,
    moment=[0, 0, 4.5], form_factor="JFE3"
)
```

### File I/O

```python
from pycrysfml08 import read_cif, read_cfl, read_diffraction_data

# CIF files
cell, spacegroup, atoms = read_cif("structure.cif")

# CFL files (CrysFML format with magnetic structures)
mag_struct = read_cfl("magnetic.cfl")

# Diffraction data
two_theta, intensity = read_diffraction_data("pattern.dat")
```

## Current Limitations

Due to gfort2py's handling of Fortran derived types and CLASS polymorphism, some CrysFML08 features have limited Python access:

| Feature | Status |
|---------|--------|
| Scattering tables | ✅ Full support |
| Math functions | ✅ Full support |
| CIF/CFL parsing | ✅ Python implementation |
| Space group lookup | ⚠️ Types accessible, some methods limited |
| Crystal cell calculations | ⚠️ Types accessible, some methods limited |
| Structure factors | ❌ Not yet implemented |
| Powder pattern simulation | ❌ Not yet implemented |

Pure Python implementations are provided where Fortran access is limited.

## Requirements

- Python ≥ 3.9
- NumPy ≥ 1.20
- gfort2py ≥ 2.6
- gfortran-14 (for building from source)

## License

This Python package is distributed under the **LGPL-3.0** license, consistent with CrysFML.

CrysFML2008 is free software for academic use. Commercial users should contact the authors.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Citation

If you use this software in your research, please cite:

**CrysFML:**
> Rodriguez-Carvajal, J. & Gonzalez-Platas, J. (2023). CrysFML2008: A Fortran Library for Crystallographic Computing. Institut Laue-Langevin.

**FullProf Suite:**
> Rodriguez-Carvajal, J. (1993). Recent advances in magnetic structure determination by neutron powder diffraction. Physica B, 192, 55-69.

## Support

- **Issues:** https://github.com/scattering/pycrysfml08/issues
- **CrysFML mailing list:** fullprof-crystallography@ill.fr

## Related Projects

- [CrysFML2008](https://code.ill.fr/scientific-software/CrysFML2008) - The Fortran library
- [FullProf Suite](https://www.ill.eu/sites/fullprof/) - Crystallographic software suite
- [pycrysfml (original)](https://github.com/scattering/pycrysfml) - Original SWIG-based bindings
- [gfort2py](https://github.com/rjfarmer/gfort2py) - Fortran-Python interface
