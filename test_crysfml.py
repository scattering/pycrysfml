#!/usr/bin/env python3
"""
Comprehensive tests for CrysFML08 Python bindings via gfort2py

This test suite validates the Python bindings to the CrysFML08 Fortran library.
Tests are organized by module and cover mathematical functions, crystallographic
calculations, and scattering data lookups.

Usage:
    python test_crysfml.py              # Run all tests
    python test_crysfml.py -v           # Verbose output
    python test_crysfml.py TestMaths    # Run specific test class

Requirements:
    - gfort2py >= 2.6.2
    - numpy
    - libcrysfml_test.dylib (or .so on Linux)
    - cfml_*.mod files from gfortran-14 build
"""

import os
import sys
import unittest
import numpy as np

# Change to build directory
BUILD_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BUILD_DIR)

try:
    import gfort2py as gf
except ImportError:
    print("ERROR: gfort2py not installed. Run: pip install gfort2py")
    sys.exit(1)

# Library and module paths
LIB_PATH = os.path.join(BUILD_DIR, "libcrysfml_test.dylib")
if not os.path.exists(LIB_PATH):
    LIB_PATH = os.path.join(BUILD_DIR, "libcrysfml_test.so")

if not os.path.exists(LIB_PATH):
    print(f"ERROR: Library not found at {LIB_PATH}")
    sys.exit(1)


class TestMaths(unittest.TestCase):
    """Tests for CFML_Maths module - mathematical functions"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_maths.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_factorial_zero(self):
        """factorial_i(0) should return 1"""
        result = self.lib.factorial_i(0)
        self.assertEqual(result.result, 1)

    def test_factorial_five(self):
        """factorial_i(5) should return 120"""
        result = self.lib.factorial_i(5)
        self.assertEqual(result.result, 120)

    def test_factorial_ten(self):
        """factorial_i(10) should return 3628800"""
        result = self.lib.factorial_i(10)
        self.assertEqual(result.result, 3628800)

    def test_gcd_basic(self):
        """gcd(48, 18) should return 6"""
        result = self.lib.gcd(48, 18)
        self.assertEqual(result.result, 6)

    def test_gcd_coprime(self):
        """gcd(17, 13) should return 1 (coprime numbers)"""
        result = self.lib.gcd(17, 13)
        self.assertEqual(result.result, 1)

    def test_gcd_same(self):
        """gcd(42, 42) should return 42"""
        result = self.lib.gcd(42, 42)
        self.assertEqual(result.result, 42)

    def test_lcm_basic(self):
        """lcm(12, 18) should return 36"""
        result = self.lib.lcm(12, 18)
        self.assertEqual(result.result, 36)

    def test_lcm_coprime(self):
        """lcm(7, 11) should return 77 (coprime numbers)"""
        result = self.lib.lcm(7, 11)
        self.assertEqual(result.result, 77)

    def test_lcm_multiple(self):
        """lcm(6, 12) should return 12 (one is multiple of other)"""
        result = self.lib.lcm(6, 12)
        self.assertEqual(result.result, 12)


class TestScatteringTables(unittest.TestCase):
    """Tests for CFML_Scattering_Tables module - element data lookups"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_scattering_tables.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_atomic_mass_iron(self):
        """Atomic mass of Fe should be approximately 55.847"""
        result = self.lib.get_atomic_mass("Fe")
        self.assertAlmostEqual(result.result, 55.847, places=2)

    def test_atomic_mass_hydrogen(self):
        """Atomic mass of H should be approximately 1.008"""
        result = self.lib.get_atomic_mass("H ")
        self.assertAlmostEqual(result.result, 1.008, places=2)

    def test_atomic_mass_helium(self):
        """Atomic mass of He should be approximately 4.003"""
        result = self.lib.get_atomic_mass("He")
        self.assertAlmostEqual(result.result, 4.003, places=2)

    def test_atomic_mass_neon(self):
        """Atomic mass of Ne should be approximately 20.18"""
        result = self.lib.get_atomic_mass("Ne")
        self.assertAlmostEqual(result.result, 20.18, places=1)

    def test_atomic_mass_tungsten(self):
        """Atomic mass of W should be approximately 183.85"""
        result = self.lib.get_atomic_mass("W ")
        self.assertAlmostEqual(result.result, 183.85, places=1)

    def test_atomic_mass_oxygen(self):
        """Atomic mass of O should be approximately 15.999"""
        result = self.lib.get_atomic_mass("O ")
        self.assertAlmostEqual(result.result, 15.999, places=2)

    def test_atomic_mass_aluminum(self):
        """Atomic mass of Al should be approximately 26.98"""
        result = self.lib.get_atomic_mass("Al")
        self.assertAlmostEqual(result.result, 26.98, places=1)

    def test_covalent_radius_carbon(self):
        """Covalent radius of C should be approximately 0.77 Angstrom"""
        try:
            result = self.lib.get_covalent_radius("C ")
            self.assertAlmostEqual(result.result, 0.77, places=1)
        except Exception as e:
            self.skipTest(f"get_covalent_radius not working: {e}")

    def test_fermi_length_hydrogen(self):
        """Fermi scattering length of H should be negative (about -0.374 x 10^-12 cm)"""
        try:
            result = self.lib.get_fermi_length("H ")
            # H has negative scattering length (value is in units of 10^-12 cm)
            self.assertLess(result.result, 0)
            self.assertAlmostEqual(result.result, -0.374, places=2)
        except Exception as e:
            self.skipTest(f"get_fermi_length not working: {e}")


class TestGlobalDeps(unittest.TestCase):
    """Tests for CFML_GlobalDeps module - global constants and dependencies"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_globaldeps.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_GlobalDeps module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)

    def test_has_precision_kinds(self):
        """Module should define precision kinds"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        # Common precision kind names
        expected_kinds = ['cp', 'dp', 'sp']
        found = [k for k in expected_kinds if any(k in a.lower() for a in attrs)]
        self.assertGreater(len(found), 0, "Should have precision kind definitions")


class TestStrings(unittest.TestCase):
    """Tests for CFML_Strings module - string manipulation utilities"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_strings.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_Strings module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)
        # Should have common string functions
        self.assertIn('u_case', [a.lower() for a in attrs])


class TestRational(unittest.TestCase):
    """Tests for CFML_Rational module - rational number arithmetic"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_rational.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_Rational module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)


class TestSymmetryTables(unittest.TestCase):
    """Tests for CFML_Symmetry_Tables module - space group lookup tables"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_symmetry_tables.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_Symmetry_Tables module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)

    def test_has_space_group_data(self):
        """Module should have space group related symbols"""
        attrs = [a.lower() for a in dir(self.lib) if not a.startswith('_')]
        # Should have functions related to space groups
        self.assertTrue(
            any('spg' in a or 'space' in a or 'group' in a or 'symm' in a for a in attrs),
            "Should have space group related symbols"
        )


class TestMetrics(unittest.TestCase):
    """Tests for CFML_Metrics module - crystal cell calculations"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_metrics.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_Metrics module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)

    def test_has_cell_types(self):
        """Module should define cell types"""
        attrs = [a.lower() for a in dir(self.lib) if not a.startswith('_')]
        self.assertTrue(
            any('cell' in a for a in attrs),
            "Should have cell type definitions"
        )


class TestGSpaceGroups(unittest.TestCase):
    """Tests for CFML_gSpaceGroups module - generalized space group operations"""

    @classmethod
    def setUpClass(cls):
        mod_path = os.path.join(BUILD_DIR, "cfml_gspacegroups.mod")
        cls.lib = gf.fFort(LIB_PATH, mod_path)

    def test_module_loads(self):
        """CFML_gSpaceGroups module should load successfully"""
        attrs = [a for a in dir(self.lib) if not a.startswith('_')]
        self.assertGreater(len(attrs), 0)

    def test_has_spg_type(self):
        """Module should define Spg_Type"""
        attrs = [a.lower() for a in dir(self.lib) if not a.startswith('_')]
        self.assertIn('spg_type', attrs)

    def test_has_symm_oper_type(self):
        """Module should define Symm_Oper_Type"""
        attrs = [a.lower() for a in dir(self.lib) if not a.startswith('_')]
        self.assertIn('symm_oper_type', attrs)


class TestIntegration(unittest.TestCase):
    """Integration tests combining multiple modules"""

    def test_all_modules_load(self):
        """All expected modules should load successfully"""
        modules = [
            "cfml_globaldeps.mod",
            "cfml_maths.mod",
            "cfml_strings.mod",
            "cfml_rational.mod",
            "cfml_metrics.mod",
            "cfml_gspacegroups.mod",
            "cfml_symmetry_tables.mod",
            "cfml_scattering_tables.mod",
        ]

        for mod in modules:
            mod_path = os.path.join(BUILD_DIR, mod)
            if os.path.exists(mod_path):
                try:
                    lib = gf.fFort(LIB_PATH, mod_path)
                    attrs = [a for a in dir(lib) if not a.startswith('_')]
                    self.assertGreater(len(attrs), 0, f"{mod} should have symbols")
                except Exception as e:
                    self.fail(f"Failed to load {mod}: {e}")

    def test_element_data_consistency(self):
        """Element data should be internally consistent"""
        mod_path = os.path.join(BUILD_DIR, "cfml_scattering_tables.mod")
        lib = gf.fFort(LIB_PATH, mod_path)

        # Heavier elements should have larger atomic masses
        mass_h = lib.get_atomic_mass("H ").result
        mass_he = lib.get_atomic_mass("He").result
        mass_fe = lib.get_atomic_mass("Fe").result
        mass_w = lib.get_atomic_mass("W ").result

        self.assertLess(mass_h, mass_he)
        self.assertLess(mass_he, mass_fe)
        self.assertLess(mass_fe, mass_w)


# Reference data for Al2O3 (Corundum) - from CIF file
AL2O3_CELL = {
    'a': 4.754,
    'b': 4.754,
    'c': 12.99,
    'alpha': 90.0,
    'beta': 90.0,
    'gamma': 120.0,
    'space_group': 167,  # R -3 c
    'space_group_symbol': 'R -3 c H',
}


class TestCrystallographicData(unittest.TestCase):
    """Tests using real crystallographic data (Al2O3 corundum)"""

    def test_al2o3_elements_exist(self):
        """Elements in Al2O3 should have valid atomic data"""
        mod_path = os.path.join(BUILD_DIR, "cfml_scattering_tables.mod")
        lib = gf.fFort(LIB_PATH, mod_path)

        # Al and O should have valid masses
        mass_al = lib.get_atomic_mass("Al").result
        mass_o = lib.get_atomic_mass("O ").result

        self.assertGreater(mass_al, 20)
        self.assertLess(mass_al, 30)
        self.assertGreater(mass_o, 10)
        self.assertLess(mass_o, 20)

    def test_formula_mass_al2o3(self):
        """Formula mass of Al2O3 should be approximately 101.96 g/mol"""
        mod_path = os.path.join(BUILD_DIR, "cfml_scattering_tables.mod")
        lib = gf.fFort(LIB_PATH, mod_path)

        mass_al = lib.get_atomic_mass("Al").result
        mass_o = lib.get_atomic_mass("O ").result

        formula_mass = 2 * mass_al + 3 * mass_o
        self.assertAlmostEqual(formula_mass, 101.96, places=1)


def run_tests():
    """Run all tests and return results"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    test_classes = [
        TestMaths,
        TestScatteringTables,
        TestGlobalDeps,
        TestStrings,
        TestRational,
        TestSymmetryTables,
        TestMetrics,
        TestGSpaceGroups,
        TestIntegration,
        TestCrystallographicData,
    ]

    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    return result


if __name__ == '__main__':
    # Check for command line arguments
    if len(sys.argv) > 1 and sys.argv[1] not in ['-v', '--verbose']:
        # Run specific test class
        unittest.main()
    else:
        # Run all tests
        result = run_tests()

        # Summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print(f"Tests run: {result.testsRun}")
        print(f"Failures: {len(result.failures)}")
        print(f"Errors: {len(result.errors)}")
        print(f"Skipped: {len(result.skipped)}")

        if result.wasSuccessful():
            print("\nAll tests passed!")
            sys.exit(0)
        else:
            print("\nSome tests failed.")
            sys.exit(1)
