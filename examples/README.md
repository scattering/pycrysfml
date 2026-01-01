# CrysFML08 Python Examples

This folder contains examples demonstrating how to use the CrysFML08 Fortran library
from Python via the gfort2py bindings.

## Prerequisites

1. Build the CrysFML08 shared library (see `../BUILDING_CRYSFML_PYTHON.md`)
2. Activate the virtual environment:
   ```bash
   source ../venv/bin/activate
   ```

## Examples

### 1. `al2o3_example.py` - Crystal Structure Analysis
Demonstrates:
- Reading a CIF file
- Extracting cell parameters, space group, and atom positions
- Looking up element properties from CrysFML scattering tables
- Calculating formula mass and density
- Reading diffraction data files

```bash
python al2o3_example.py
```

### 2. `element_properties.py` - Element Database
Demonstrates:
- Looking up atomic masses for all elements
- Covalent radii
- Neutron Fermi scattering lengths
- Atomic volumes
- Calculating formula masses for compounds
- Using math functions (gcd, lcm, factorial)

```bash
python element_properties.py
```

## API Module

### `crysfml_api.py`
A Python wrapper providing:

```python
from crysfml_api import CrysFML

cfml = CrysFML()

# Scattering tables
mass = cfml.get_atomic_mass("Fe")       # 55.847 amu
radius = cfml.get_covalent_radius("Fe") # 1.17 Angstrom
fermi = cfml.get_fermi_length("Fe")     # 0.945 x 10^-12 cm

# Formula mass
mass = cfml.formula_mass({"Fe": 2, "O": 3})  # Fe2O3

# Math functions
cfml.factorial(5)  # 120
cfml.gcd(48, 18)   # 6
cfml.lcm(12, 18)   # 36

# CIF parsing
cell, spacegroup, atoms = cfml.read_cif("structure.cif")
```

## Data Files

- `Al2O3.cif` - Corundum crystal structure (R-3c, #167)
- `Al2O3.dat` - Diffraction pattern data
- `Al2O3 Background.BGR` - Background spline points
- `pbso4.cfl` - Lead sulfate structure (Pnma, #62)

## Current Limitations

Due to gfort2py's handling of Fortran derived types and CLASS polymorphism,
some CrysFML08 features are not directly accessible:

### What Works
- Scattering tables (atomic mass, radii, scattering lengths)
- Mathematical functions (factorial, gcd, lcm)
- Module loading and symbol access

### Limited Support
- Space group operations (types load but methods may fail)
- Crystal cell calculations (types load but some methods fail)
- Complex derived type operations

### Workarounds
- CIF parsing is done in Python rather than using CrysFML's parser
- Crystal calculations are done in Python using the cell parameters

## Porting from pycrysfml

The original pycrysfml used SWIG-generated wrappers that had full access to
Fortran derived types. If you're porting code from pycrysfml:

| pycrysfml | gfort2py version |
|-----------|------------------|
| `H.readInfo(cif_file)` | `cfml.read_cif(cif_file)` |
| `funcs.get_atomic_mass(sym)` | `cfml.get_atomic_mass(sym)` |
| `cell.volume` | `cell.volume` (property) |
| `SpaceGroup(name)` | Not directly supported |
| `diffPattern(...)` | Not yet implemented |

## Contributing

To add more examples:
1. Add data files to this folder
2. Create a new Python script following the existing patterns
3. Use the `crysfml_api.CrysFML` class for CrysFML functions
4. Add documentation to this README
