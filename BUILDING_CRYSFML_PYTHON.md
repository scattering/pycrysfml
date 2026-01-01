# Building CrysFML08 Python Bindings

This document describes two approaches for creating Python bindings to the CrysFML08 Fortran library:
1. **LFortran** - Modern LLVM-based Fortran compiler (experimental, partially working)
2. **gfort2py** - Python bindings via gfortran (recommended, fully working)

---

## Approach 1: LFortran (Experimental)

### Current Status: Partially Working

LFortran is a modern Fortran compiler built on LLVM with potential for WebAssembly compilation and better tooling integration. However, it currently lacks support for several Fortran features used by CrysFML08.

### Installation

```bash
# Install via Homebrew
brew tap fortran-lang/fortran
brew install lfortran

# Verify installation
lfortran --version  # Should show 0.48.0 or later
```

### Known Issues and Workarounds

#### Issue 1: z3 Library Version Mismatch (macOS)
LFortran may be linked against an older z3 version than what's installed.

**Fix:**
```bash
# Create symlink if needed (adjust versions as necessary)
ln -sf /opt/homebrew/opt/z3/lib/libz3.4.15.dylib /opt/homebrew/opt/z3/lib/libz3.4.14.dylib
```

#### Issue 2: `protected` Attribute Not Supported
LFortran doesn't support the `protected` attribute on module variables.

**Location:** `CFML_Maths.f90` (1 occurrence)

**Fix:**
```bash
sed 's/, protected//' CFML_Maths.f90 > CFML_Maths_patched.f90
```

#### Issue 3: Internal File I/O Not Supported (BLOCKING)
LFortran doesn't support writing to internal character variables:
```fortran
write(unit=char_var, fmt="(i40)") integer_value  ! Not supported
```

**Location:** `CFML_Rational.f90`

**Workaround:** Create a pure Fortran helper module:

```fortran
! CFML_IntStr.f90 - Replacement for internal file I/O
Module CFML_IntStr
    use CFML_GlobalDeps, only: LI
    implicit none
    private
    public :: int_to_str, rational_to_str

Contains

    Pure Function int_to_str(n) Result(s)
        integer(kind=LI), intent(in) :: n
        character(len=40) :: s
        integer(kind=LI) :: absn, digit
        integer :: i, ndigits
        logical :: negative
        character(len=40) :: temp

        s = ' '
        if (n == 0_LI) then
            s = '0'
            return
        end if

        negative = (n < 0_LI)
        absn = abs(n)
        ndigits = 0
        temp = ' '

        do while (absn > 0_LI)
            ndigits = ndigits + 1
            digit = mod(absn, 10_LI)
            temp(41-ndigits:41-ndigits) = char(ichar('0') + int(digit))
            absn = absn / 10_LI
        end do

        if (negative) then
            ndigits = ndigits + 1
            temp(41-ndigits:41-ndigits) = '-'
        end if

        s = adjustl(temp(41-ndigits:40))
    End Function int_to_str

    Pure Function rational_to_str(num, denom) Result(s)
        integer(kind=LI), intent(in) :: num, denom
        character(len=80) :: s

        if (denom == 1_LI) then
            s = int_to_str(num)
        else
            s = trim(adjustl(int_to_str(num))) // '/' // trim(adjustl(int_to_str(denom)))
        end if
    End Function rational_to_str

End Module CFML_IntStr
```

Then patch CFML_Rational.f90:
```bash
sed -e 's/Use CFML_Strings,    only : Pack_String/Use CFML_Strings,    only : Pack_String\n    Use CFML_IntStr,     only : rational_to_str/' \
    -e 's/write(unit=line,fmt="(i40)") sr%numerator/line = rational_to_str(sr%numerator, 1_LI)/' \
    -e 's/write(unit=line,fmt="(i40,a,i40)") sr%numerator,"\/",sr%denominator/line = rational_to_str(sr%numerator, sr%denominator)/' \
    CFML_Rational_orig.f90 > CFML_Rational.f90
```

#### Issue 4: Derived Type Array Constructors Not Supported (BLOCKING, NO WORKAROUND)
LFortran doesn't support array constructors for derived types:
```fortran
type(Xray_Form_Type), dimension(n) :: data
data = [Xray_Form_Type("H", ...), Xray_Form_Type("He", ...), ...]  ! Not supported
```

**Location:** `CFML_Scattering_Tables.f90`

**Status:** No workaround available. This blocks compilation of the scattering tables module, which contains essential crystallographic data.

### Modules That Work with LFortran

| Module | Status | Notes |
|--------|--------|-------|
| CFML_GlobalDeps | Working | No modifications needed |
| CFML_Maths | Working | Remove `protected` attribute |
| CFML_Strings | Working | No modifications needed |
| CFML_Metrics | Working | No modifications needed |
| CFML_Random | Working | No modifications needed |
| CFML_Rational | Working | Requires CFML_IntStr workaround |
| CFML_Scattering_Tables | **BLOCKED** | Derived type array constructors |
| CFML_gSpaceGroups | **BLOCKED** | Depends on CFML_Scattering_Tables |

### What LFortran Needs to Support CrysFML08

1. **Internal file I/O** - `write(unit=char_var, ...)` syntax
2. **Derived type array constructors** - `[Type(...), Type(...), ...]`
3. **Protected attribute** - Minor, easily patched

### Future Potential

If LFortran adds support for the missing features:
- Native Python bindings could be generated
- WebAssembly compilation would enable browser-based crystallography tools
- Better debugging and error messages than gfortran

---

## Approach 2: gfort2py (Recommended)

### Current Status: Fully Working

gfort2py provides Python bindings to Fortran code by parsing gfortran's `.mod` files and using ctypes to call the compiled shared library.

### Requirements

- **gfortran 14.x** (gfortran 15 creates module version 16, not yet supported by gfort2py)
- **Python 3.x** (must match architecture of compiled library)
- **gfort2py 2.6.2+**
- **numpy**

### Installation

#### Step 1: Install gfortran-14

```bash
# macOS with Homebrew
brew install gcc@14

# Verify version
/opt/homebrew/bin/gfortran-14 --version
# Should show: GNU Fortran (Homebrew GCC 14.x.x) 14.x.x
```

#### Step 2: Create Python Virtual Environment

**Important:** On Apple Silicon Macs, you must use ARM Python, not x86_64 (Anaconda/Miniconda).

```bash
# Use Homebrew's Python (ARM native on Apple Silicon)
/opt/homebrew/bin/python3 -m venv venv

# Activate
source venv/bin/activate

# Install dependencies
pip install gfort2py numpy
```

### Build Process

#### Step 1: Set Up Build Directory

```bash
GFORTRAN="/opt/homebrew/bin/gfortran-14"
BUILDDIR="/path/to/build"
SRC="/path/to/CrysFML2008/Src"

mkdir -p "$BUILDDIR"
cd "$BUILDDIR"
```

#### Step 2: Copy and Prepare Source Files

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

#### Step 3: Compile Modules (Order Matters!)

```bash
FFLAGS="-c -fPIC -ffree-line-length-none -O2"

# 1. Base dependencies (no deps)
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

# 8. gSpaceGroups (depends on most modules above)
$GFORTRAN $FFLAGS CFML_gSpaceGroups.f90 -o CFML_gSpaceGroups.o
for f in CFML_gSpaceGroups/*.f90; do
    $GFORTRAN $FFLAGS "$f" -o "${f%.f90}.o"
done
```

#### Step 4: Create Shared Library

```bash
# Collect all object files
OBJECTS=$(find . -name "*.o" | tr '\n' ' ')

# Link into shared library
$GFORTRAN -shared -o libcrysfml.dylib $OBJECTS

# Verify
ls -la libcrysfml.dylib
file libcrysfml.dylib  # Should show: Mach-O 64-bit dynamically linked shared library arm64
```

### Using gfort2py

```python
#!/usr/bin/env python3
import gfort2py as gf
import os

# Load library and module
lib_path = "libcrysfml.dylib"
mod_path = "cfml_maths.mod"

lib = gf.fFort(lib_path, mod_path)

# List available symbols
symbols = [a for a in dir(lib) if not a.startswith('_')]
print(f"Available symbols: {len(symbols)}")

# Call a function
result = lib.factorial_i(5)
print(f"factorial_i(5) = {result.result}")  # Output: 120
```

### Working Modules and Symbol Counts

| Module | Symbols | Description |
|--------|---------|-------------|
| cfml_globaldeps | 29 | Global dependencies, precision kinds |
| cfml_maths | 57 | Mathematical functions |
| cfml_strings | 49 | String manipulation |
| cfml_rational | 16 | Rational number arithmetic |
| cfml_metrics | 44 | Crystal cell metrics |
| cfml_gspacegroups | 73 | Space group operations |
| cfml_symmetry_tables | 69 | Symmetry lookup tables |
| cfml_scattering_tables | 50 | X-ray/neutron scattering factors |
| **Total** | **387** | |

### Limitations of gfort2py

1. **gfortran 15 not supported** - Module version 16 is not yet recognized by gfort2py 2.6.2
2. **Architecture must match** - Library and Python must both be ARM64 or both x86_64
3. **Complex derived types** - Some complex Fortran types may require numpy arrays
4. **Class polymorphism** - Fortran CLASS types may have limited support
5. **Allocatable arrays** - Require proper memory management from Python side
6. **String handling** - Fixed-length Fortran strings need proper padding

### Troubleshooting

#### "Unsupported module version"
```
gfort2py only supports module versions ['0', '15'], got 16
```
**Solution:** Use gfortran-14 instead of gfortran-15

#### "wrong architecture" or "mach-o but wrong architecture"
```
OSError: dlopen(...): mach-o file, but is an incompatible architecture
```
**Solution:** Ensure Python and library have matching architecture:
```bash
file libcrysfml.dylib  # Check library arch
python -c "import platform; print(platform.machine())"  # Check Python arch
```

#### Missing symbols at link time
```
Undefined symbols for architecture arm64: "_cfml_strings_..."
```
**Solution:** Ensure all submodule .f90 files are compiled (check CFML_Strings/*.f90, etc.)

---

## Comparison Summary

| Feature | LFortran | gfort2py |
|---------|----------|----------|
| **Status** | Experimental | Production-ready |
| **Modules working** | 6 of 10+ | All tested |
| **Source modification** | Required | None |
| **Python bindings** | Future potential | Working now |
| **WebAssembly** | Potential | No |
| **Compiler required** | LFortran 0.48+ | gfortran-14 |

**Recommendation:** Use gfort2py with gfortran-14 for immediate work. Monitor LFortran development for future WebAssembly capabilities.
