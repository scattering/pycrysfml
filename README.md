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
pip install pycrysfml
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
   git clone https://github.com/scattering/pycrysfml.git
   cd pycrysfml
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
from pycrysfml import CrysFML, read_cif, read_cfl

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
from pycrysfml import CrysFML

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
from pycrysfml import CrystalCell, Atom, SpaceGroup, MagneticAtom

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
from pycrysfml import read_cif, read_cfl, read_diffraction_data

# CIF files
cell, spacegroup, atoms = read_cif("structure.cif")

# CFL files (CrysFML format with magnetic structures)
mag_struct = read_cfl("magnetic.cfl")

# Diffraction data
two_theta, intensity = read_diffraction_data("pattern.dat")
```

## Current Limitations and gfort2py Technical Details

### Overview

This package uses [gfort2py](https://github.com/rjfarmer/gfort2py) to interface Python with the compiled Fortran library. gfort2py reads gfortran `.mod` files and the shared library to expose Fortran routines to Python. While this approach works well for many use cases, certain Fortran 2008 features present challenges.

### Feature Support Matrix

| Feature | Status | Notes |
|---------|--------|-------|
| Scattering tables | ✅ Full | `get_atomic_mass`, `get_fermi_length`, etc. |
| Math functions | ✅ Full | `factorial`, `gcd`, `lcm`, array operations |
| String utilities | ✅ Full | Character manipulation functions |
| CIF/CFL parsing | ✅ Python | Pure Python implementation |
| Crystal cell | ⚠️ Partial | Python dataclass with methods; Fortran type read-only |
| Space groups | ⚠️ Partial | Table lookups work; type-bound methods limited |
| Atom lists | ⚠️ Partial | Can read types; allocatable arrays need care |
| Structure factors | ❌ Planned | Requires derived type I/O |
| Powder simulation | ❌ Planned | Complex type hierarchies |
| Magnetic structure factors | ❌ Planned | CLASS polymorphism issues |

### What Works Well

**Simple functions with primitive types:**
```python
# These work perfectly - scalar in, scalar out
mass = cfml.get_atomic_mass("Fe")        # Returns float
n = cfml.factorial(5)                     # Returns int
gcd = cfml.gcd(12, 8)                     # Returns int
```

**Array operations:**
```python
# NumPy arrays map to Fortran arrays
import numpy as np
matrix = np.array([[1, 2], [3, 4]], dtype=np.float64)
det = cfml.determinant(matrix)
```

**Accessing module variables:**
```python
# Global constants and tables are accessible
lib.num_chem_info.get()  # Number of elements in database
```

### What Has Limitations

**1. Fortran Derived Types (TYPE)**

CrysFML08 uses complex derived types like `Cell_Type`, `SpGr_Type`, and `Atom_Type`. gfort2py can access these but with restrictions:

```fortran
! Fortran definition
type :: Cell_Type
    real(kind=cp) :: cell(6)      ! a, b, c, alpha, beta, gamma
    real(kind=cp) :: rcell(6)     ! Reciprocal cell
    real(kind=cp) :: Vol          ! Volume
    ! ... many more components
end type
```

```python
# Python access - reading works
cell = lib.some_cell_variable
a = cell.cell[0]  # Works
vol = cell.vol    # Works

# But creating new instances is complex
# We provide Python dataclasses instead:
from pycrysfml import CrystalCell
cell = CrystalCell(a=5.0, b=5.0, c=5.0, alpha=90, beta=90, gamma=90)
```

**2. Type-Bound Procedures (Methods)**

Fortran 2003+ allows methods on types. These are challenging to call via gfort2py:

```fortran
! Fortran - type-bound procedure
type :: Cell_Type
contains
    procedure :: get_volume => calc_volume
end type
```

```python
# Direct method calls may not work
# cell.get_volume()  # May fail

# Workaround: call module procedure directly or use Python implementation
volume = cell.a * cell.b * cell.c * volume_factor(cell.alpha, cell.beta, cell.gamma)
```

**3. CLASS Polymorphism**

CrysFML08 uses CLASS for polymorphic types (especially in space group handling):

```fortran
! Fortran polymorphic type
class(SpG_Type), allocatable :: SpaceGroup
```

gfort2py has difficulty with:
- Allocating polymorphic variables
- Determining runtime type
- Calling type-bound procedures on CLASS variables

**4. Allocatable Components**

Derived types with allocatable array components require special handling:

```fortran
type :: Atom_List_Type
    integer :: natoms
    type(Atom_Type), dimension(:), allocatable :: atom  ! Allocatable array
end type
```

```python
# Reading may work if allocated in Fortran first
# But allocating from Python is problematic
```

**5. Intent(OUT) Derived Types**

Functions that return derived types or have INTENT(OUT) derived type arguments often fail:

```fortran
subroutine Set_SpaceGroup(symbol, SpG)
    character(len=*), intent(in) :: symbol
    type(SpG_Type), intent(out) :: SpG  ! Problematic
end subroutine
```

### Our Workarounds

**1. Python Dataclasses**

We provide pure Python equivalents for common types:

```python
@dataclass
class CrystalCell:
    a: float
    b: float
    c: float
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0

    @property
    def volume(self) -> float:
        # Pure Python calculation
        ...

    def d_spacing(self, h: int, k: int, l: int) -> float:
        # Pure Python calculation
        ...
```

**2. Python File Parsers**

CIF and CFL parsing is implemented in pure Python, avoiding Fortran I/O issues:

```python
def read_cif(filename: str) -> Tuple[CrystalCell, SpaceGroup, List[Atom]]:
    # Pure Python parsing
    ...
```

**3. Wrapper Functions**

Where possible, we wrap Fortran calls with Python error handling:

```python
def get_atomic_mass(self, symbol: str) -> float:
    """Get atomic mass with proper error handling."""
    self._ensure_tables_loaded()
    result = self._scattering_lib.get_atomic_mass(symbol.strip())
    return result.result
```

### Future Improvements

1. **Structure Factor Calculations**: Implement Python wrappers that call Fortran for heavy computation
2. **Powder Pattern Simulation**: Build Python interface to core Fortran routines
3. **Space Group Operations**: Expose more symmetry operations via direct function calls
4. **Contribute to gfort2py**: Help improve derived type and CLASS support upstream

### Alternative Approaches

If you need full CrysFML08 functionality, consider:

1. **ILL's pycrysfml**: Uses f2py with manual interface files
   - https://code.ill.fr/scientific-software/pycrysfml

2. **Direct Fortran**: Write a small Fortran program that calls CrysFML08 and outputs results

3. **C Interoperability**: CrysFML08 has some C-compatible interfaces via `ISO_C_BINDING`

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

- **Issues:** https://github.com/scattering/pycrysfml/issues
- **CrysFML mailing list:** fullprof-crystallography@ill.fr

## Related Projects

- [CrysFML2008](https://code.ill.fr/scientific-software/CrysFML2008) - The Fortran library
- [FullProf Suite](https://www.ill.eu/sites/fullprof/) - Crystallographic software suite
- [pycrysfml (original)](https://github.com/scattering/pycrysfml) - Original SWIG-based bindings
- [gfort2py](https://github.com/rjfarmer/gfort2py) - Fortran-Python interface
