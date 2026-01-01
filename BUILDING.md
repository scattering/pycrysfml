# Building PyCrysFML08

This document provides detailed instructions for building the CrysFML08
Fortran library and the PyCrysFML08 Python bindings.

## Prerequisites

### Compiler Requirements

**gfortran 14.x is required** (not gfortran 15).

The gfort2py library currently supports gfortran module version 15, which is
produced by gfortran-14. gfortran-15 produces module version 16, which is
not yet supported.

#### macOS (Apple Silicon or Intel)

```bash
# Install gfortran-14 via Homebrew
brew install gcc@14

# Verify installation
/opt/homebrew/bin/gfortran-14 --version
# Should show: GNU Fortran (Homebrew GCC 14.x.x) 14.x.x
```

#### Ubuntu/Debian Linux

```bash
# Add GCC 14 repository if needed
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update

# Install gfortran-14
sudo apt install gfortran-14

# Verify
gfortran-14 --version
```

### Python Requirements

- **Python 3.9 or later**
- **Architecture must match**: ARM64 Python for ARM64 library, x86_64 for x86_64

**Important for macOS Apple Silicon users:**
Anaconda/Miniconda typically installs x86_64 Python via Rosetta. You need
native ARM64 Python to use ARM64-compiled libraries.

```bash
# Install native ARM64 Python via Homebrew
brew install python@3.12

# Create virtual environment with ARM64 Python
/opt/homebrew/bin/python3 -m venv venv

# Activate
source venv/bin/activate

# Install dependencies
pip install numpy gfort2py
```

## Obtaining CrysFML2008 Source

Clone the CrysFML2008 repository from ILL GitLab:

```bash
git clone https://code.ill.fr/scientific-software/CrysFML2008.git
```

## Build Process

### Step 1: Set Up Build Environment

```bash
# Set compiler path
export GFORTRAN="/opt/homebrew/bin/gfortran-14"  # macOS
# or
export GFORTRAN="/usr/bin/gfortran-14"  # Linux

# Set paths
export SRC="/path/to/CrysFML2008/Src"
export BUILDDIR="/path/to/build"

mkdir -p "$BUILDDIR"
cd "$BUILDDIR"
```

### Step 2: Copy Source Files

```bash
# Copy main module files
cp "$SRC"/CFML_GlobalDeps.f90 .
cp "$SRC"/CFML_Maths.f90 .
cp "$SRC"/CFML_Strings.f90 .
cp "$SRC"/CFML_Rational.f90 .
cp "$SRC"/CFML_Metrics.f90 .
cp "$SRC"/CFML_gSpaceGroups.f90 .
cp "$SRC"/CFML_Symmetry_Tables.f90 .
cp "$SRC"/CFML_Scattering_Tables.f90 .

# Copy submodule directories
cp -r "$SRC"/CFML_Maths .
cp -r "$SRC"/CFML_Strings .
cp -r "$SRC"/CFML_Rational .
cp -r "$SRC"/CFML_Metrics .
cp -r "$SRC"/CFML_gSpaceGroups .
cp -r "$SRC"/CFML_Scattering_Tables .
```

### Step 3: Compile Modules

Modules must be compiled in dependency order:

```bash
FFLAGS="-c -fPIC -ffree-line-length-none -O2"

# 1. Global dependencies (no dependencies)
$GFORTRAN $FFLAGS CFML_GlobalDeps.f90 -o CFML_GlobalDeps.o

# 2. Maths (depends on GlobalDeps)
$GFORTRAN $FFLAGS CFML_Maths.f90 -o CFML_Maths.o
for f in CFML_Maths/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done

# 3. Strings (depends on GlobalDeps, Maths)
$GFORTRAN $FFLAGS CFML_Strings.f90 -o CFML_Strings.o
for f in CFML_Strings/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done

# 4. Rational (depends on GlobalDeps, Maths, Strings)
$GFORTRAN $FFLAGS CFML_Rational.f90 -o CFML_Rational.o
for f in CFML_Rational/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done

# 5. Metrics (depends on GlobalDeps, Maths, Strings, Rational)
$GFORTRAN $FFLAGS CFML_Metrics.f90 -o CFML_Metrics.o
for f in CFML_Metrics/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done

# 6. Symmetry Tables (depends on GlobalDeps)
$GFORTRAN $FFLAGS CFML_Symmetry_Tables.f90 -o CFML_Symmetry_Tables.o

# 7. Scattering Tables (depends on GlobalDeps, Strings)
$GFORTRAN $FFLAGS CFML_Scattering_Tables.f90 -o CFML_Scattering_Tables.o
for f in CFML_Scattering_Tables/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done

# 8. gSpaceGroups (depends on most modules)
$GFORTRAN $FFLAGS CFML_gSpaceGroups.f90 -o CFML_gSpaceGroups.o
for f in CFML_gSpaceGroups/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done
```

### Step 4: Create Shared Library

```bash
# Collect all object files
OBJECTS=$(find . -name "*.o" | tr '\n' ' ')

# Link into shared library
# macOS:
$GFORTRAN -shared -o libcrysfml.dylib $OBJECTS

# Linux:
$GFORTRAN -shared -o libcrysfml.so $OBJECTS

# Verify
ls -la libcrysfml.*
file libcrysfml.*
```

### Step 5: Verify Module Files

```bash
# Check that .mod files were created
ls -la *.mod

# Expected files:
# cfml_globaldeps.mod
# cfml_maths.mod
# cfml_strings.mod
# cfml_rational.mod
# cfml_metrics.mod
# cfml_gspacegroups.mod
# cfml_symmetry_tables.mod
# cfml_scattering_tables.mod
```

## Testing the Build

### Quick Test

```python
import gfort2py as gf

lib = gf.fFort("libcrysfml.dylib", "cfml_maths.mod")
result = lib.factorial_i(5)
print(f"5! = {result.result}")  # Should print: 5! = 120
```

### Full Test Suite

```bash
# From the package directory
python -m pytest tests/
```

## Build Script

A convenience build script is provided:

```bash
#!/bin/bash
# scripts/build_library.sh

set -e

# Configuration
GFORTRAN="${GFORTRAN:-gfortran-14}"
SRC="${CRYSFML_SRC:-../CrysFML2008/Src}"
BUILDDIR="${BUILDDIR:-.}"

FFLAGS="-c -fPIC -ffree-line-length-none -O2 -J$BUILDDIR"

echo "Building CrysFML08 library..."
echo "Compiler: $GFORTRAN"
echo "Source: $SRC"
echo "Build dir: $BUILDDIR"

# [Build commands as above...]

echo "Build complete!"
echo "Library: $BUILDDIR/libcrysfml.dylib (or .so)"
```

## Troubleshooting

### "Unsupported module version"

```
gfort2py only supports module versions ['0', '15'], got 16
```

**Solution:** Use gfortran-14 instead of gfortran-15.

### "wrong architecture"

```
OSError: dlopen(...): mach-o file, but is an incompatible architecture
```

**Solution:** Ensure Python architecture matches library architecture:

```bash
# Check library
file libcrysfml.dylib  # Should show arm64 or x86_64

# Check Python
python -c "import platform; print(platform.machine())"

# If mismatched, create a venv with matching Python
```

### Undefined symbols

```
Undefined symbols for architecture arm64: "_cfml_strings_..."
```

**Solution:** Ensure all submodule .f90 files are compiled. Check that files
in CFML_Strings/*.f90, CFML_Maths/*.f90, etc. are included.

### "Module file not found"

**Solution:** Ensure .mod files are in the same directory as the library,
or set the PYCRYSFML08_MODS environment variable.

## Alternative: LFortran (Experimental)

LFortran is a modern Fortran compiler with potential for WebAssembly output.
However, it currently lacks support for some Fortran 2008 features used by
CrysFML:

- Internal file I/O (`write(unit=char_var, ...)`)
- Derived type array constructors
- Protected attribute

See `BUILDING_CRYSFML_PYTHON.md` for details on using LFortran.

## Environment Variables

| Variable | Description |
|----------|-------------|
| `PYCRYSFML08_LIB` | Path to libcrysfml.dylib/.so |
| `PYCRYSFML08_MODS` | Path to directory containing .mod files |
| `GFORTRAN` | Path to gfortran compiler |
| `CRYSFML_SRC` | Path to CrysFML2008/Src directory |
