import unittest
import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from pycrysfml import *
from fswig_hklgen import readInfo, SpaceGroup, AtomList, CrystalCell, Atom
# TODO: populate tests
class TestCrysFMLIO(unittest.TestCase):
    """Test of CrysFML I/O functions, crystal cell, atom, spacegroup, magnetic symmetry, and atomList types"""
    def test_cfl_reader(self):
        """Test of CrysFML cfl reader and Magnetic data types"""
        DATAPATH = os.path.dirname(os.path.abspath(__file__))
        cfl = os.path.join(DATAPATH,r"hobk.cfl")
    def test_pcr_reader(self):
        """Test of CrysFML PCR reader"""
        DATAPATH = os.path.dirname(os.path.abspath(__file__))
        pcr = os.path.join(DATAPATH,r"pbso4.pcr")
        spaceGroup, cell, atomList = readInfo(pcr)
        self.assertEqual(np.allclose(np.array([8.478, 5.3967, 6.958]), np.array(cell.length()), rtol=0.0001), True,
                         msg="Crystal Cell populated incorrectly (a, b, c)")
        self.assertEqual(np.allclose(np.array([90.0, 90.0, 90.0]), np.array(cell.angle())), True,
                         msg="Crystal Cell populated incorrectly (alpha, beta, gamma)")
        self.assertEqual(spaceGroup.xtalSystem.lower(), "orthorhombic", msg="SpaceGroup error: xtalSystem")
        self.assertEqual(spaceGroup.number, 62, msg="SpaceGroup error: number")
        self.assertEqual(len(atomList), 5, msg="AtomList length = "+str(len(atomList))+" Should be 5")
        atom = atomList[0]
        self.assertEqual(atom.label().lower(), "pb", msg="Wrong atom in atomList")
        self.assertEqual(np.allclose(np.array([0.188, 0.25, 0.167]), np.array(atom.coords()), rtol=0.01), True,
                         msg="Atom coordinates wrong")
        self.assertEqual(atom.occupancy(), 0.5, msg="Atom Occupancy wrong")
    def test_cif_reader(self):
        """ Test of CrysFML cif reader"""
        DATAPATH = os.path.dirname(os.path.abspath(__file__))
        cif = os.path.join(DATAPATH,r"al2o3.cif")
        spaceGroup, cell, atomList = readInfo(cif)
        self.assertEqual(np.allclose(np.array([4.7698, 4.769, 13.024]), np.array(cell.length()), rtol=0.001), True,
                         msg="Crystal Cell populated incorrectly (a, b, c)")
        self.assertEqual(np.allclose(np.array([90.0, 90.0, 120.0]), np.array(cell.angle())), True,
                         msg="Crystal Cell populated incorrectly (alpha, beta, gamma)")
        self.assertEqual(spaceGroup.xtalSystem.lower(), "trigonal", msg="SpaceGroup error: xtalSystem")
        self.assertEqual(spaceGroup.number, 167, msg="SpaceGroup error: number")
        self.assertEqual(len(atomList), 2, msg="AtomList length = "+str(len(atomList))+" Should be 2")
        atom = atomList[0]
        self.assertEqual(atom.label().lower(), "al1", msg="Wrong atom in atomList")
        self.assertEqual(np.allclose(np.array([0.0, 0.00, 0.352]), np.array(atom.coords()), rtol=0.01), True,
                         msg="Atom coordinates wrong")
        self.assertAlmostEqual(atom.occupancy(), 0.3333333, msg="Atom Occupancy wrong")        
    def test_mcif_reader(self):
        """ Test of CrysFML mcif reader and Magnetic Data types """
        DATAPATH = os.path.dirname(os.path.abspath(__file__))
        mcif = os.path.join(DATAPATH,r"lfo.mcif")
        pass
class TestSatelliteGen(unittest.TestCase):
    pass
class TestMiscCalculations(unittest.TestCase):
    pass
class TestStrFactors(unittest.TestCase):
    """ Test of nuclear and magnetic structure factor calculations for a variety of crystal cells """
    def test_mag_orthohombic(self):
        """Test magnetic structure factor calculations for orthorhombic crystall cell"""
        pass
    def test_mag_hexagonal(self):
        """Test magnetic structure factor calculations for hexagonal crystall cell"""
        pass
    