#!/usr/bin/env python3
"""Test gfort2py with CrysFML08 modules"""

import os
import sys

# Add path to the build directory
build_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(build_dir)

try:
    import gfort2py as gf
    print(f"gfort2py version: {gf.__version__}")
except ImportError as e:
    print(f"Failed to import gfort2py: {e}")
    sys.exit(1)

# Try to load the library
lib_path = os.path.join(build_dir, "libcrysfml_test.dylib")
mod_path = os.path.join(build_dir, "cfml_gspacegroups.mod")

print(f"\nLibrary: {lib_path}")
print(f"Module:  {mod_path}")
print(f"Library exists: {os.path.exists(lib_path)}")
print(f"Module exists:  {os.path.exists(mod_path)}")

try:
    print("\nLoading library with gfort2py...")
    lib = gf.fFort(lib_path, mod_path)
    print("SUCCESS: Library loaded!")

    # List available functions/variables
    print("\nAvailable symbols:")
    attrs = [a for a in dir(lib) if not a.startswith('_')]
    for attr in attrs[:20]:  # First 20
        print(f"  - {attr}")
    if len(attrs) > 20:
        print(f"  ... and {len(attrs) - 20} more")

except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

# Also try loading a simpler module
print("\n--- Testing simpler module (cfml_maths) ---")
mod_path2 = os.path.join(build_dir, "cfml_maths.mod")
try:
    lib2 = gf.fFort(lib_path, mod_path2)
    print("SUCCESS: cfml_maths loaded!")
    attrs = [a for a in dir(lib2) if not a.startswith('_')]
    print(f"Found {len(attrs)} symbols")
    for attr in attrs[:10]:
        print(f"  - {attr}")
except Exception as e:
    print(f"ERROR: {e}")
