import unittest
from pycrysfml import *
# TODO: populate tests
class TestCrysFMLIO(unittest.TestCase):
    def test_cfl_reader(self):
        pass
    def test_pcr_reader(self):
        pass
    def test_cif_reader(self):
        pass
    def test_mcif_reader(self):
        pass
class TestSpaceGroup(unittest.TestCase):
    pass
class TestAtom(unittest.TestCase):
    pass
class TestMagAtom(unittest.TestCase):
    pass
class TestAtomList(unittest.TestCase):
    pass
class TestCrystallCell(unittest.TestCase):
    pass
class TestMagAtomList(unittest.TestCase):
    pass
class TestMagSymmetry(unittest.TestCase):
    pass
class TestSatelliteGen(unittest.TestCase):
    pass
class TestReflection(unittest.TestCase):
    pass
class TestMagReflection(unittest.TestCase):
    pass
class TestReflectionList(unittest.TestCase):
    pass
class TestMagReflectionList(unittest.TestCase):
    pass
class TestMiscCalculations(unittest.TestCase):
    pass
class TestStrFactors(unittest.TestCase):
    def test_test(self):
        "test of green"
        self.assertEqual(type("hello"), str, msg="hello")
    def test_mag_orthohombic(self):
        "Test magnetic structure factor calculations for orthorhombic crystall cell"
        pass
    def test_mag_hexagonal(self):
        "Test magnetic structure factor calculations for hexagonal crystall cell"
        pass
    