# hklGen.py
# Generates predicted diffraction peak positions for a given unit cell and space
#   group, using the Fortran CFML library.
# Also uses the program "bumps" to perform fitting of calculated diffraction
#   patterns to observed data.
# Created 6/4/2013
# Last edited 6/24/2013

import sys
from math import floor, sqrt, log, tan, radians
from string import rstrip
from ctypes import cdll, Structure, c_int, c_float, c_char, c_bool, c_char_p, \
                   c_void_p, c_ulong, addressof, sizeof, POINTER

import numpy as np
import pylab

#np.seterr(all='raise')
import lattice_calculator_procedural2 as latcalc

# This should be the location of the CFML library
lib = cdll["./libcrysfml.so"]

# DVDim: contains array dimension information. Attributes:
#   stride_mult - stride multiplier for the dimension
#   lower_bound - first index for the dimension
#   upper_bound - last index for the dimension
class DVDim(Structure):
    _fields_ = [("stride_mult", c_ulong), ("lower_bound", c_ulong),
                ("upper_bound", c_ulong)]

# DV: dope vector for gfortran that passes array information. Attributes:
#   base_addr   - base address for the array
#   base        - base offset
#   dtype       - contains the element size, type (3 bits), and rank (3 bits)
#   dim         - DVDim object for the vector
class DV(Structure):
    _fields_ = [("base_addr",c_void_p), ("base", c_void_p), ("dtype", c_ulong),
                ("dim", DVDim*7)]
  
    def __len__(self):
        return self.dim[0].upper_bound - self.dim[0].lower_bound + 1
    def elemsize(self):
        return self.dv.dtype >> 6
    def rank(self):
        return self.dv.dtype&7
    def data(self, dataType):
        return (dataType*len(self)).from_address(self.base_addr)

# dv_dtype: calculates the "dtype" attribute for a dope vector
def dv_dtype(size,type,rank): return size*64+type*8+rank

# build_struct_dv: constructs a dope vector for an array of derived types
def build_struct_dv(array):
    dv = DV()
    dv.base_addr = addressof(array)
    dv.base = c_void_p()
    dv.dtype = dv_dtype(sizeof(array[0]), 5, 1) # 5 = derived type
    dv.dim[0].stride_mult = 1
    dv.dim[0].lower_bound = 1
    dv.dim[0].upper_bound = len(array)
    return dv

# deconstruct_dv: converts a dope vector into an array of any type
def deconstruct_dv(dv, dataType):
#    dims = dv.dtype & 7
    elem_size = dv.dtype >> 6
    size = dv.dim[0].upper_bound - dv.dim[0].lower_bound + 1
    array = [None for i in xrange(size)]
    for i in xrange(size):
        offset = i * elem_size * dv.dim[0].stride_mult
        array[i] = dataType.from_address(dv.base_addr + offset)
    return array

# SymmetryOp Attributes:
#   rot     - rotational part of symmetry operator (3 by 3 matrix)
#   trans   - translational part of symmetry operator (vector)
class SymmetryOp(Structure):
    _fields_ = [("rot", c_int*9), ("trans", c_float*3)]

# WyckoffPos attributes:
#   multip      - multiplicity
#   site        - site symmetry
#   numElements - number of elements in orbit
#   origin      - origin
#   orbit       - strings containing orbit information
class WyckoffPos(Structure):
    _fields_ = [("multip", c_int), ("site", c_char*6),
                ("numElements", c_int), ("origin", c_char*40),
                ("orbit", c_char*40*48)]

# Wyckoff attributes:
#   numOrbits   - number of orbits
#   orbits      - list of Wyckoff position objects
class Wyckoff(Structure):
    _fields_ = [("numOrbits", c_int), ("orbits", WyckoffPos*26)]

# SpaceGroup attributes:
#   number          - space group number
#   symbol          - Hermann-Mauguin symbol
#   hallSymbol      - Hall symbol
#   xtalSystem      - Crystal system
#   laue            - Laue class
#   pointGroup      - corresponding point group
#   info            - additional information
#   setting         - space group setting information (IT, KO, ML, ZA, Table,
#                     Standard, or UnConventional)
#   hex             - true if space group is hexagonal
#   lattice         - lattice type
#   latticeSymbol   - lattice type symbol
#   latticeNum      - number of lattice points in a cell
#   latticeTrans    - lattice translations
#   bravais         - Bravais symbol and translations
#   centerInfo      - information about symmetry center
#   centerType      - 0 = centric (-1 not at origin), 1 = acentric,
#                     2 = centric (-1 at origin)
#   centerCoords    - fractional coordinates of inversion center
#   numOps          - number of symmetry operators in the reduced set
#   multip          - multiplicity of the general position
#   numGens         - minimum number of operators to generate the group
#   symmetryOps     - list of symmetry operators
#   symmetryOpsSymb - string form of symmetry operator objects
#   wyckoff         - object containing Wyckoff information
#   asymmetricUnit  - direct space parameters for the asymmetric unit
class SpaceGroup(Structure):
    # TODO: Create nice-looking constructors for these objects
    _fields_ = [("number", c_int), ("symbol", c_char*20),
                ("hallSymbol", c_char*16), ("xtalSystem", c_char*12),
                ("laue", c_char*5), ("pointGroup", c_char*5),
                ("info", c_char*5), ("setting", c_char*80),
                ("hex", c_int), ("lattice", c_char),
                ("latticeSymbol", c_char*2), ("latticeNum", c_int),
                ("latticeTrans", c_float*3*12), ("bravais", c_char*51),
                ("centerInfo", c_char*80), ("centerType", c_int),
                ("centerCoords", c_float*3), ("numOps", c_int),
                ("multip", c_int), ("numGens", c_int),
                ("symmetryOps", SymmetryOp*192),
                ("symmetryOpsSymb", c_char*40*192),
                ("wyckoff", Wyckoff), ("asymmetricUnit", c_float*6)]

# CrystalCell attributes:
#   length, angle           - arrays of unit cell parameters
#   lengthSD, angleSD       - standard deviations of parameters
#   rLength, rAngle         - arrays of reciprocal cell parameters
#   GD, GR                  - direct and reciprocal space metric tensors
#   xtalToOrth, orthToXtal  - matrices to convert between orthonormal and
#                             crystallographic bases
#   BLB, invBLB             - Busing-Levy B-matrix and its inverse
#   volume, rVolume         - direct and reciprocal cell volumes
#   cartType                - Cartesian reference frame type (cartType = 'A'
#                             designates x || a)
class CrystalCell(Structure):
    _fields_ = [("length", c_float*3), ("angle", c_float*3),
                ("lengthSD", c_float*3), ("angleSD", c_float*3),
                ("rLength", c_float*3), ("rAngle", c_float*3),
                ("GD", c_float*9), ("GR", c_float*9),
                ("xtalToOrth", c_float*9), ("orthToXtal", c_float*9),
                ("BLB", c_float*9), ("invBLB", c_float*9),
                ("volume", c_float), ("rVolume", c_float),
                ("cartType", c_char)]

# Atom attributes:
#   label       - label for the atom
#   element     - chemical symbol of the element 
#   sFactorSymb - chemical symbol used in the structure factor
#   active      - used for program control
#   atomicNum   - atomic number
#   multip      - site multiplicity
#   coords      - fractional coordinates
#   occupancy   - site occupancy
#   BIso        - isotropic temperature (Debye-Waller) factor
#   uType       - type of anisotropic thermal factor: "u_ij", "b_ij", "beta",
#                     or none
#   thType      - thermal factor type: "isotr", "aniso", or "other"
#   U           - matrix entries U11, U22, U33, U12, U13, and U23
#   UEquiv      - equivalent U
#   charge      - self-explanatory
#   moment      - magnetic moment
#   index       - index "for different purposes", whatever mysterious purposes
#                     those might be
#   numVars     - a variable completely lacking in documentation
#   freeVars    - free variables used "for different purposes" again
#   atomInfo    - string containing miscellaneous information
class Atom(Structure):
    _fields_ = [("label", c_char*10), ("element", c_char*2),
                ("sFactorSymb", c_char*4), ("active", c_int),
                ("atomicNum", c_int), ("multip", c_int), ("coords", c_float*3),
                ("coordsSD", c_float*3), ("multCoords", c_float*3),
                ("lsqCoords", c_int*3), ("occupancy", c_float),
                ("occupancySD", c_float), ("multOccupancy", c_float),
                ("lsqOccupancy", c_int), ("BIso", c_float), ("BIsoSD", c_float),
                ("multBIso", c_float), ("lsqBIso", c_float),
                ("uType", c_char*4), ("thType", c_char*5), ("U", c_float*6),
                ("USD", c_float*6), ("UEquiv", c_float), ("multU", c_float*6),
                ("lsqU", c_int*6), ("charge", c_float), ("moment", c_float),
                ("index", c_int*5), ("numVars", c_int),
                ("freeVars", c_float*10), ("atomInfo", c_char*40)]

# AtomList attributes:
#   numAtoms    - the number of atoms
#   atoms       - a list of Atom objects
class AtomList(Structure):
    _fields_ = [("numAtoms", c_int), ("atoms", DV)]

# Reflection attributes:
#   hkl         - list containing hkl indices for the reflection
#   multip      - multiplicity
#   FObs        - observed structure factor
#   FCalc       - calculated structure factor
#   FSD         - standard deviation of structure factor
#   s           - s = sin(theta)/lambda = 1/(2d) [No 4*pi factor!]
#   weight      - weight of reflection
#   phase       - phase angle in degrees
#   realPart    - real part of the structure factor
#   imPart      - imaginary part of the structure factor
#   aa          - currently unused
#   bb          - currently unused
class Reflection(Structure):
    _fields_ = [("hkl", c_int*3), ("multip", c_int), ("FObs", c_float), 
                ("FCalc", c_float), ("FSD", c_float), ("s", c_float),
                ("weight", c_float), ("phase", c_float), ("realPart", c_float),
                ("imPart", c_float), ("aa", c_float), ("bb", c_float)]

# ReflectionList attributes
#   numReflections  - the number of reflections
#   reflections     - list of Reflection objects
class ReflectionList(Structure):
    _fields_ = [("numReflections", c_int), ("reflections", DV)]

#    def __init__(self):
#        self.reflections.dtype = 0
#        self.reflections.base_addr = 0
#        self.reflections.base = 0

#    def __init__(self, reflections):
#        self.numReflections = c_int(len(reflections))
#        Reflection68 = Reflection*68
#        self.reflections = byref(build_struct_dv(Reflection68(*reflections)))
#        self.reflections = Reflection68(*reflections)

# Gaussian: represents a Gaussian function that can be evaluated at any
#   2*theta value. u, v, and w are fitting parameters.
class Gaussian(object):
    scaleFactor = 1    
    
    def __init__(self, center, u, v, w, I, hkl=[0,0,0]):
        self.C0 = 4*log(2)
        self.center = center    # 2*theta position
        self.u = u
        self.v = v
        self.w = w
        self.I = I
        try:
            self.H = sqrt(u*(tan(radians(center/2))**2)
                          + v*tan(radians(center/2)) + w)
            self.scale = self.I * sqrt(self.C0/np.pi)/self.H * Gaussian.scaleFactor
        except ValueError:
            self.H = 0
            self.scale = 0
        self.hkl = hkl

    # __call__: returns the value of the Gaussian at some 2*theta positions
    def __call__(self, x):
        return self.scale * np.exp(-self.C0*(x-self.center)**2/self.H**2)

    def add(self, v, x):
        # only add to nearby 2*theta positions
        idx = (x>self.center-self.H*3) & (x<self.center+self.H*3)
        v[idx] += self.__call__(x[idx])

# LinSpline: represents a linear spline function to be used for the background
class LinSpline(object):
    def __init__(self, arg1, arg2=None):
        if type(arg1) == type(np.array([])):
            # read in x and y coordinates from lists
            self.x = np.copy(arg1)
            self.y = np.copy(arg2)
        elif type(arg1) == str:
            # read in x and y coordinates from a file
            self.x, self.y = np.loadtxt(arg1, dtype=float, skiprows=5, unpack=True)

    # __call__: returns the interpolated y value at some x position
    def __call__(self, x):
        # locate the two points to interpolate between
        return np.interp(x, self.x, self.y)

    def __repr__(self):
        return "LinSpline(" + str(self.x) + ", " + str(self.y) + ")"

# twoTheta: converts a sin(theta)/lambda position to a 2*theta position
def twoTheta(s, wavelength):
    if (s*wavelength >= 1): return 180.0
    if (s*wavelength <= 0): return 0.0
    return 2*np.degrees(np.arcsin(s*wavelength))

# getS: converts a 2*theta position to a sin(theta)/lambda position
def getS(tt, wavelength):
    return np.sin(np.radians(tt/2))/wavelength

# hklString: converts an hkl list into a string
def hklString(hkl):
    return "%d %d %d" % tuple(hkl)

# setSpaceGroup: constructs a SpaceGroup object from a provided symbol/number
def setSpaceGroup(name, spaceGroup):
    #print >>sys.stderr, name, spaceGroup
    fn = lib.__cfml_crystallographic_symmetry_MOD_set_spacegroup
    fn.argtypes = [c_char_p, POINTER(SpaceGroup), c_void_p, POINTER(c_int),
                   c_char_p, c_char_p, c_int, c_int, c_int]
    fn.restype = None
    fn(name, spaceGroup, None, None, None, None, len(name), 0, 0)
    return

# setCrystalCell: constructs a CrystalCell object from provided parameters
def setCrystalCell(length, angle, cell):
    fn = lib.__cfml_crystal_metrics_MOD_set_crystal_cell
    float3 = c_float*3
    fn.argtypes = [POINTER(float3), POINTER(float3), POINTER(CrystalCell),
                   POINTER(c_char), POINTER(float3)]
    fn.restype = None
    fn(float3(*length), float3(*angle), cell, None, None)
    return

# getMaxNumRef: returns the maximum number of reflections for a given cell
def getMaxNumRef(sMax, volume, sMin=0.0, multip=2):
    fn = lib.__cfml_reflections_utilities_MOD_get_maxnumref
    fn.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_float),
                   POINTER(c_int)]
    fn.restype = c_int
    numref = fn(c_float(sMax), c_float(volume), c_float(sMin), c_int(multip))
    return numref

# hklUni: constructs a list of unique reflections in a specified range
#   If code == "r", then d-spacings are used
def hklUni(cell, spaceGroup, sMin, sMax, code, numRef, maxRef, getList=False):
    if (not getList):
        fn = lib.__cfml_reflections_utilities_MOD_hkl_uni_reflection
        c_ReflectionArray = Reflection*max(maxRef,1)
        reflections = c_ReflectionArray()
        fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup), POINTER(c_bool),
                       POINTER(c_float), POINTER(c_float), c_char_p,
                       POINTER(c_int), POINTER(DV), POINTER(c_bool)]
        fn.restype = None
        fn(cell, spaceGroup, c_bool(True), c_float(sMin), c_float(sMax),
           code, numRef, build_struct_dv(reflections), c_bool(False))
        return reflections
    else:
        fn = lib.__cfml_reflections_utilities_MOD_hkl_uni_refllist
        fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup),
                       POINTER(c_bool), POINTER(c_float), POINTER(c_float),
                       c_char_p, POINTER(c_int), POINTER(ReflectionList),
                       POINTER(c_bool)]
        fn.restype = None
        refList = ReflectionList()
        fn(cell, spaceGroup, c_bool(True), c_float(sMin), c_float(sMax), code,
           numRef, refList, c_bool(False))
        return refList

# hklGen: generates a list of reflections
def hklGen(spaceGroup, cell, wavelength, sMin, sMax, getList=False):
    # Calculate the reflection positions
    maxReflections = getMaxNumRef(sMax+0.2, cell.volume,
                                  multip=spaceGroup.multip)
    # Create a reference that will be modified by hklUni()
    reflectionCount = c_int(maxReflections)
    reflections = hklUni(cell, spaceGroup, sMin, sMax, 's', reflectionCount,
                         maxReflections, getList)
    if (type(reflections) == type([])):
        reflections = reflections[:reflectionCount.value]    
    return reflections

# setAtoms: creates an atom list object and places atoms in specified positions
def setAtoms(elements, atomicNums=None, coords=None, occupancies=None, BIsos=None):
    init = lib.__cfml_atom_typedef_MOD_allocate_atom_list
    init.argtypes = [POINTER(c_int), POINTER(AtomList), POINTER(c_int)]
    init.restype = None
    atomList = AtomList()
    numAtoms = len(elements)
    init(c_int(numAtoms), atomList, None)

    atomList.numAtoms = numAtoms
    atoms = atomList.atoms.data(Atom)
    c_array3 = c_float*3
    if isinstance(elements[0], Atom):
        atoms = np.copy(elements)
    else:
        for i in xrange(numAtoms):
            atoms[i].element = elements[i]
            atoms[i].atomicNum = c_int(atomicNums[i])
            atoms[i].coords = c_array3(*coords[i])
            atoms[i].occupancy = c_float(occupancies[i])
            atoms[i].BIso = c_float(BIsos[i])    
    return atomList

# calcStructFact: calculates the structure factor squared for a list of planes
#   using provided atomic positions
def calcStructFact(refList, atomList, spaceGroup, wavelength):
    init = lib.__cfml_structure_factors_MOD_init_calc_strfactors
    init.argtypes = [POINTER(ReflectionList), POINTER(AtomList),
                     POINTER(SpaceGroup), c_char_p, POINTER(c_float),
                     POINTER(c_int), c_int]
#    init.argtypes = [POINTER(AtomList), c_char_p, POINTER(c_float),
#                     POINTER(c_int), c_int]  
    init.restype = None
    init(refList, atomList, spaceGroup, 'NUC', c_float(wavelength), None, 3)
    
    calc = lib.__cfml_structure_factors_MOD_calc_strfactor
    calc.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_float),
                     POINTER(AtomList), POINTER(SpaceGroup), POINTER(c_float),
                     POINTER(DV), POINTER(c_float*2), c_int, c_int]
    calc.restype = None
    structFacts = [c_float() for i in xrange(refList.numReflections)]
    reflections = refList.reflections.data(Reflection)
    for i, reflection in enumerate(reflections):
        # calculates the square of the structure factor
        calc('P', 'NUC', c_int(i+1), c_float(reflection.s**2), atomList,
             spaceGroup, structFacts[i], None, None, 1, 3)
    # convert from c_float to float
    structFacts = [sf.value for sf in structFacts]
    return structFacts

# calcIntensity: calculates the intensity for a given set of reflections,
#   based on the structure factor
def calcIntensity(refList, atomList, spaceGroup, wavelength):
    reflections = refList.reflections.data(Reflection)
    sfs = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength))
    multips = np.array([ref.multip for ref in reflections])
    
    tt = np.radians(np.array([twoTheta(ref.s, wavelength) for ref in reflections]))
    lorentz = (1+np.cos(tt)**2) / (np.sin(tt)*np.sin(tt/2))
    return sfs * multips * lorentz

# inputInfo: requests space group and cell information, wavelength, and range
#   of interest
def inputInfo():
    # Input the space group name and create it
    groupName = raw_input("Enter space group "
                          "(HM symbol, Hall symbol, or number): ")
    spaceGroup = SpaceGroup()
    setSpaceGroup(groupName, spaceGroup)

    # Remove excess spaces from Fortran
    spaceGroup.xtalSystem = rstrip(spaceGroup.xtalSystem)
    length = [0, 0, 0]
    angle = [0, 0, 0]

    # Ask for parameters according to the crystal system and create the cell
    if (spaceGroup.xtalSystem.lower() == "triclinic"):
        length[0], length[1], length[2], angle[0], angle[1], angle[2] = \
            input("Enter cell parameters (a, b, c, alpha, beta, gamma): ")
    elif (spaceGroup.xtalSystem.lower() == "monoclinic"):
        angle[0] = angle[2] = 90
        length[0], length[1], length[2], angle[1] = \
            input("Enter cell parameters (a, b, c, beta): ")
    elif (spaceGroup.xtalSystem.lower() == "orthorhombic"):
        angle[0] = angle[1] = angle[2] = 90
        length[0], length[1], length[2] = \
            input("Enter cell parameters (a, b, c): ")
    elif (spaceGroup.xtalSystem.lower() == "tetragonal"):
        angle[0] = angle[1] = angle[2] = 90
        length[0], length[2] = input("Enter cell parameters (a, c): ")
        length[1] = length[0]
    elif (spaceGroup.xtalSystem.lower() in ["rhombohedral", "hexagonal"]):
        angle[0] = angle[1] = 90
        angle[2] = 120
        length[0], length[2] = input("Enter cell parameters (a, c): ")
        length[1] = length[0]
    elif (spaceGroup.xtalSystem.lower() == "cubic"):
        angle[0] = angle[1] = angle[2] = 90
        length[0] = input("Enter cell parameter (a): ")
        length[1] = length[2] = length[0]
    cell = CrystalCell()
    setCrystalCell(length, angle, cell)

    # Input the wavelength and range for hkl calculation (and adjust the range
    #   if necessary)
    wavelength = input("Enter the wavelength: ")
    sMin, sMax = input("Enter the sin(theta)/lambda interval: ")
    adjusted = False
    if (sMin < 0):
        sMin = 0
        adjusted = True
    if (sMax > 1.0/wavelength):
        sMax = 1.0/wavelength
        adjusted = True
    if (adjusted):
        print "sin(theta)/lambda interval adjusted to [%f, %f]" % (sMin, sMax)
    return (spaceGroup, cell, wavelength, sMin, sMax)

# printReflections: outputs spacegroup information and a list of reflections
def printReflections(reflections, spaceGroup, wavelength, sMin, sMax):
    reflectionCount = len(reflections)
    print
    print "Space group information"
    print "---------------------------------------"
    print "              Number: ", spaceGroup.number
    print "          H-M Symbol: ", spaceGroup.symbol
    print "         Hall Symbol: ", spaceGroup.hallSymbol
    print "      Crystal System: ", spaceGroup.xtalSystem
    print "          Laue Class: ", spaceGroup.laue
    print "         Point Group: ", spaceGroup.pointGroup
    print "General Multiplicity: ", spaceGroup.multip
    print "---------------------------------------"
    print
    print "%d reflections found for %f < sin(theta)/lambda < %f" \
            % (reflectionCount, sMin, sMax)
    print "                        (%f < 2*theta < %f)" % \
            (twoTheta(sMin, wavelength), twoTheta(sMax, wavelength))
    print
    print " h k l  mult  sin(theta)/lambda  2*theta"
    for refl in reflections:
        print "(%d %d %d)  %d        %f       %f" % \
                (refl.hkl[0], refl.hkl[1], refl.hkl[2], refl.multip, refl.s,
                 twoTheta(refl.s, wavelength))
    print
    return

# makeGaussians() creates a series of Gaussians to represent the powder
#   diffractionn pattern
def makeGaussians(reflections, coeffs, I, scale, wavelength):
    Gaussian.scaleFactor = scale
    gaussians = [Gaussian(twoTheta(rk.s, wavelength),
                          coeffs[0], coeffs[1], coeffs[2], Ik, rk.hkl)
                 for rk,Ik in zip(reflections,I)]
    return gaussians

# getIntensity: calculates the intensity at a given 2*theta position, or for an
#   array of 2*theta positions
def getIntensity(gaussians, background, tt):
    #return background(tt) + sum(g(tt) for g in gaussians)
    v = background(tt)
    for g in gaussians:
        g.add(v,tt)
    return v

# plotPattern: given a series of Gaussians and a background, plots the predicted
#   intensity at every 2*theta position in a specified range, as well as the
#   observed intensity everywhere on a given list of points
def plotPattern(gaussians, background, ttObs, observed, ttMin, ttMax, ttStep,
                exclusions=None, labels=None):
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    ttCalc = np.linspace(ttMin, ttMax, numPoints)
    if(exclusions != None): ttCalc = removeRange(ttCalc, exclusions)
    intensity = np.array(getIntensity(gaussians, background, ttCalc))
    pylab.plot(ttCalc, intensity, '-', label="Caclulated")
    pylab.plot(ttObs, observed, '-', label="Observed")
    pylab.xlabel(r"$2 \theta$")
    pylab.ylabel("Intensity")
    pylab.legend()
    if (labels == "hkl"):
        for g in gaussians:
            if (g.center <= ttMax):
                pylab.text(g.center, np.interp(g.center, ttCalc, intensity),
                           hklString(g.hkl),
                           ha="center", va="bottom", rotation="vertical")
    return

# removeRange: takes in an array of 2*theta intervals and removes them from
#   consideration for data analysis, with an optional argument for removing the
#   corresponding intensities as well
def removeRange(tt, remove, intensity=None):
    if (remove == None):
        if (intensity != None): return (tt, intensity)
        else: return tt
    if (type(remove[0]) not in [list, type(np.array([]))]):
        # single interval
        keepEntries = (tt < remove[0]) | (tt > remove[1])
        tt = tt[keepEntries]
        if (intensity != None):
            intensity = intensity[keepEntries]
            return (tt, intensity)
        else: return tt
    else:
        # array of intervals
        if (intensity != None):
            for interval in remove:
                tt, intensity = removeRange(tt, interval, intensity)
            return (tt, intensity)
        else:
            for interval in remove:
                tt  = removeRange(tt, interval)
            return tt

# readFile: acquires cell, space group, and atomic information from a .cif,
#   .cfl, .pcr, or .shx file
def readFile(filename):
    fn = lib.__cfml_io_formats_MOD_readn_set_xtal_structure_split
    fn.argtypes = [c_char_p, POINTER(CrystalCell), POINTER(SpaceGroup),
                   POINTER(AtomList), c_char_p, POINTER(c_int),
                   c_void_p, c_void_p, c_char_p, c_int, c_int, c_int]
    fn.restype = None
    cell = CrystalCell()
    spaceGroup = SpaceGroup()
    atomList = AtomList()
    fn(filename, cell, spaceGroup, atomList, filename[-3:], None, None, None,
       None, len(filename), 3, 0)
    print filename, filename[-3:]
    return (spaceGroup, cell, atomList)

#def chiSquare(observed, gaussians, background, tt):
#    expected = np.array(getIntensity(gaussians, background, tt))
#    return sum((observed-expected)**2/expected)

# Model: represents an object that can be used with bumps for optimization
#   purposes
class Model(object):

    def __init__(self, tt, observed, background, u, v, w,
                 wavelength, spaceGroupName, cell, atoms, exclusions=None):
        if (isinstance(spaceGroupName, SpaceGroup)):
            self.spaceGroup = spaceGroupName
        else:
            self.spaceGroup = SpaceGroup()
            setSpaceGroup(spaceGroupName, self.spaceGroup)
        self.tt = tt
        self.observed = observed
        self.background = background
        self.u = Parameter(u, name='u')
        self.v = Parameter(v, name='v')
        self.w = Parameter(w, name='w')
        self.scale = Parameter(1, name='scale')
        self.wavelength = wavelength
        self.cell = cell
        self.exclusions = exclusions
        self.ttMin = min(self.tt)
        self.ttMax = max(self.tt)
        self.sMin = getS(self.ttMin, self.wavelength)
        self.sMax = getS(self.ttMax, self.wavelength)

        self._set_reflections()
        self.atomListModel = AtomListModel(atoms)
        self.update()

    def _set_reflections(self):
        maxCell = CrystalCell()
        maxLattice = self.cell.getMaxLattice()
        setCrystalCell(maxLattice[:3], maxLattice[3:], maxCell)
        self.refList = hklGen(self.spaceGroup, maxCell,
                              self.wavelength, self.sMin, self.sMax, True)
        self.maxReflections = self.refList.reflections.data(Reflection)
        self.reflections = np.copy(self.maxReflections)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["refList"]
        return state
    def __setstate__(self, state):
        self.__dict__ = state
        self._set_reflections()

    def parameters(self):
        return {'u': self.u,
                'v': self.v,
                'w': self.w,
                'scale': self.scale,
                'cell': self.cell.parameters(),
                'atoms': self.atomListModel.parameters()
                }

    def numpoints(self):
        return len(self.observed)

    def theory(self):
        return getIntensity(self.gaussians, self.background, self.tt)

    def residuals(self):
        return (self.theory() - self.observed)/np.sqrt(self.observed)

    def nllf(self):
        return np.sum(self.residuals()**2)

    def plot(self, view="linear"):
        plotPattern(self.gaussians, self.background, self.tt, self.observed,
                    self.ttMin, self.ttMax, 0.01, self.exclusions)

#    def _cache_cell_pars(self):
#        self._cell_pars = dict((k,v.value) for k,v in self.cell.items())
    def update(self):
        self.cell.update()
        self.atomListModel.update()
        lattice = latcalc.Lattice(*self.cell.getLattice())
        h, k, l = np.array(zip(*[reflection.hkl for reflection in self.maxReflections]))
        ttPos = lattice.calc_twotheta(self.wavelength, h, k, l)
        # move nonexistent peaks out of the way to 2*theta = -20
        ttPos[np.isnan(ttPos)] = -20
        for i in xrange(len(self.reflections)):
            self.reflections[i].s = getS(ttPos[i], self.wavelength)
        self.intensities = calcIntensity(self.refList, self.atomListModel.atomList, 
                                         self.spaceGroup, self.wavelength)
        self.gaussians = makeGaussians(self.reflections,
                                       [self.u.value, self.v.value, self.w.value],
                                       self.intensities, self.scale.value,
                                       self.wavelength)

# Triclinic/Monoclinic/Orthorhombic/Tetragonal/Hexagonal/CubicCell: classes
#   that contain lattice information with refinable parameters to interface
#   with bumps
class TriclinicCell(object):
    def __init__(self, a, b, c, alpha, beta, gamma):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.b = Parameter(b, name='b')
        self.c = Parameter(c, name='c')
        self.alpha = Parameter(alpha, name='alpha')
        self.beta = Parameter(beta, name='beta')
        self.gamma = Parameter(gamma, name='gamma')
        self.update()
    def parameters(self):
        return {'a': self.a, 'b': self.b, 'c': self.c,
                'alpha': self.alpha, 'beta': self.beta, 'gamma': self.gamma}
    def update(self):
        a = self.a.value
        b = self.b.value
        c = self.c.value
        alpha = self.alpha.value
        beta = self.beta.value
        gamma = self.gamma.value
        setCrystalCell([a,b,c], [alpha, beta, gamma], self.cell)
    def getLattice(self):
        return [self.a.value, self.b.value, self.c.value,
                self.alpha, self.beta, self.gamma]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.b.bounds.limits[1],
                self.c.bounds.limits[1], self.alpha.bounds.limits[1],
                self.beta.bounds.limits[1], self.gamma.bounds.limits[1]]

class MonoclinicCell(object):
    def __init__(self, a, b, c, beta):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.b = Parameter(b, name='b')
        self.c = Parameter(c, name='c')
        self.beta = Parameter(beta, name='beta')
        self.update()
    def parameters(self):
        return {'a': self.a, 'b': self.b, 'c': self.c,
                'beta': self.beta}
    def update(self):
        a = self.a.value
        b = self.b.value
        c = self.c.value
        beta = self.beta.value
        setCrystalCell([a,b,c], [90, beta, 90], self.cell)
    def getLattice(self):
        return [self.a.value, self.b.value, self.c.value,
                90, self.beta.value, 90]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.b.bounds.limits[1],
                self.c.bounds.limits[1], 90, self.beta.bounds.limits[1], 90]

class OrthorhombicCell(object):
    def __init__(self, a, b, c):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.b = Parameter(b, name='b')
        self.c = Parameter(c, name='c')
        self.update()
    def parameters(self):
        return {'a': self.a, 'b': self.b, 'c': self.c}
    def update(self):
        a = self.a.value
        b = self.b.value
        c = self.c.value
        setCrystalCell([a,b,c], [90,90,90], self.cell)
    def getLattice(self):
        return [self.a.value, self.b.value, self.c.value, 90, 90, 90]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.b.bounds.limits[1],
                self.c.bounds.limits[1], 90, 90, 90]

class TetragonalCell(object):
    def __init__(self, a, c):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.c = Parameter(c, name='c')
        self.update()
    def parameters(self):
        return {'a': self.a, 'c': self.c}
    def update(self):
        a = self.a.value
        c = self.c.value
        setCrystalCell([a,a,c], [90,90,90], self.cell)
    def getLattice(self):
        return [self.a.value, self.a.value, self.c.value, 90, 90, 90]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.a.bounds.limits[1],
                self.c.bounds.limits[1], 90, 90, 90]

class HexagonalCell(object):
    def __init__(self, a, c):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.c = Parameter(c, name='c')
        self.update()
    def parameters(self):
        return {'a': self.a, 'c': self.c}
    def update(self):
        a = self.a.value
        c = self.c.value
#        print >>sys.stderr,"updating hex with",a,c
        setCrystalCell([a,a,c], [90,90,120], self.cell)
    def getLattice(self):
        return [self.a.value, self.a.value, self.c.value, 90, 90, 120]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.a.bounds.limits[1],
                self.c.bounds.limits[1], 90, 90, 120]

class CubicCell(object):
    def __init__(self, a):
        self.cell = CrystalCell()
        self.a = Parameter(a, name='a')
        self.update()
    def parameters(self):
        return {'a': self.a}
    def update(self):
        a = self.a.value
        setCrystalCell([a,a,a], [90,90,90], self.cell)
    def getLattice(self):
        return [self.a.value, self.a.value, self.a.value, 90, 90, 90]
    def getMaxLattice(self):
        return [self.a.bounds.limits[1], self.a.bounds.limits[1],
                self.a.bounds.limits[1], 90, 90, 90]

# makeCell: creates a bumps-compatible cell object from a CrystalCell object
def makeCell(cell, xtalSystem):
    (a, b, c) = cell.length
    (alpha, beta, gamma) = cell.angle
    if (xtalSystem.lower() == "triclinic"):
        newCell = TriclinicCell(a, b, c, alpha, beta, gamma)
    elif (xtalSystem.lower() == "monoclinic"):
        newCell = MonoclinicCell(a, b, c, beta)
    elif (xtalSystem.lower() == "orthorhombic"):
        newCell = OrthorhombicCell(a, b, c)
    elif (xtalSystem.lower() == "tetragonal"):
        newCell = TetragonalCell(a,c)
    elif (xtalSystem.lower() in ["rhombohedral", "hexagonal"]):
        newCell = HexagonalCell(a,c)
    elif (xtalSystem.lower() == "cubic"):
        newCell = CubicCell(a)
    return newCell

class AtomListModel(object):    
    # TODO: make occupancy a refinable parameter

    def __init__(self, atoms):
        self._make_atom_list(atoms)
        self.B = [None] * len(self.atoms)
        for i, atom in enumerate(self.atoms):
            self.B[i] = Parameter(atom.BIso, name=rstrip(atom.label)+"_B")
#            self.occ[i] = Parameter(atom.occupancy*atom)

    def _make_atom_list(self, atoms):
        if (isinstance(atoms, AtomList)):
            self.atomList = atoms
        elif (isinstance(atoms[0], Atom)):
            self.atomList = setAtoms(atoms)
        else:
            self.atomList = setAtoms(*zip(*atoms))
        self.atoms = self.atomList.atoms.data(Atom)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["atomList"]
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        self._make_atom_list(self.atoms)

    def parameters(self):
        return dict(zip([param.name for param in self.B],
                        [param for param in self.B]))
    def update(self):
        for i, atom in enumerate(self.atoms):
            atom.BIso = self.B[i].value

# createData: creates and saves a diffraction pattern with the given parameters
def createData(spaceGroup, cell, wavelength, ttMin, ttMax, ttStep, coeffs, I, backg, fileName):
    sMin = getS(ttMin, wavelength)
    sMax = getS(ttMax, wavelength)
    reflections = hklGen(spaceGroup, cell, wavelength, sMin, sMax)
    g = makeGaussians(reflections, coeffs, I, wavelength)
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    tt = np.linspace(ttMin, ttMax, numPoints)
    intensity = np.array(getIntensity(g, backg, tt))
    np.savetxt(fileName, (tt, intensity), delimiter=" ")
    return

# testInfo: loads up space group and cell information for testing purposes
def testInfo():
    # Information for Al2O3
    length = [4.7698, 4.7698, 13.0243]
    angle = [90,90,120]
    cell = CrystalCell()
    setCrystalCell(length, angle, cell)
    spaceGroup = SpaceGroup()
    setSpaceGroup("167", spaceGroup)
    wavelength = 1.5403
    sMin, sMax = getS(3, wavelength), getS(167.8, wavelength)
    atoms = [["Al", 13, (0,0,.35231), 1.0/3, .00523],
             ["O", 8, (.3061,0,.25), 1.0/2, .00585]]
    return (spaceGroup, cell, wavelength, sMin, sMax, atoms)

# testStructFact: tests the structure factor calculations
def testStructFact():
    print "starting test"
    (spaceGroup, cell, wavelength, sMin, sMax) = testInfo()
    refList = hklGen(spaceGroup, cell, wavelength, sMin, sMax, True)
    print str(refList.numReflections) + " reflections generated"
    atomList = setAtoms(["Al", "O"], [13, 8], [(0,0,.35231), (.3061,0,.25)],
                     [12.0/36, 18.0/36], [.00523, .00585])
    print "atom list created"
    sf = calcStructFact(refList, atomList, spaceGroup, wavelength)
    print sf
    reflections = refList.reflections.data(Reflection)
    for i in xrange(len(sf)):
        print hklString(refList.reflections.data(Reflection)[i].hkl[:]), sf[i],\
                sqrt(sf[i])
    out = np.array([[ref.hkl[0] for ref in reflections],
                    [ref.hkl[1] for ref in reflections],
                    [ref.hkl[2] for ref in reflections],
                    [sf[i] for i in xrange(len(sf))]])
    out = out.transpose()
    np.savetxt("structure factors.dat", out)

def fit():
    np.seterr(divide="ignore", invalid="ignore")
    backgFile = "/home/djq/Pycrysfml/LuFeO3 Background.BGR"
    observedFile = "/home/djq/Pycrysfml/lufep001.dat"
    infoFile = "/home/djq/Pycrysfml/LuFeMnO3.cif"
    (spaceGroup, crystalCell, atoms) = readFile(infoFile)
    spaceGroup.xtalSystem = rstrip(spaceGroup.xtalSystem)
    wavelength = 1.5406
#    spaceGroup = "185"
#    atoms = [["Al", 13, (0,0,.35231), 1.0/3, .00523],
#             ["O", 8, (.3061,0,.25), 1.0/2, .00585]]
#    atoms = [["Lu", 71, (0,0,.27207), 2.0/12, .0041],
#             ["Lu", 71, (.3333,.6667,.2332), 4.0/12, .0053],
#             ["Fe", 26, (.3332,0,0), 6.0/12, .0018],
#             ["O", 8, (.303,0,.1542), 6.0/12, .012],
#             ["O", 8, (.649,0,.332), 6.0/12, .012],
#             ["O", 8, (0,0,.472), 2.0/12, .012],
#             ["O", 8, (.3333,.6667,.017), 4.0/12, .012]]
    backg = LinSpline(backgFile)
#    backg = LinSpline(np.array([0,1]),np.array([0,0]))
    tt = np.linspace(3, 167.75, 3296)
#    exclusions = np.array([[60,66], [72.75, 75.75]])
    exclusions = np.array([[0,20],[37.3,39.1],[43.85,45.9],[64.25,66.3],
                           [76.15,79.8],[81.7,83.1],[89.68,99.9],[109.95,111.25],
                           [115.25,118.45],[133.95,141.25],[156.7,180]])
    data = np.loadtxt(observedFile, dtype=float, skiprows=3)
#    tt = data[0]
    observed = data.flatten()
    observed = observed[:2*len(tt):2]
    tt, observed = removeRange(tt, exclusions, observed)
#    observed = data.flatten()[:len(tt)]
#    cell = HexagonalCell(4.7698, 13.0243)
#    cell = HexagonalCell(5.965, 11.702)
    cell = makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
    cell.c.pm(0.5)
    m = Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
              atoms, exclusions)
    m.u.range(0,1)
    m.v.range(-1,0)
    m.w.range(0,10)
    m.scale.range(0,10)
    for i in xrange(len(m.atomListModel.B)):
        m.atomListModel.B[i].range(0, .1)
    M = FitProblem(m)
    #M.model_update()
    return M

def main():
#    (spaceGroup, cell, wavelength, sMin, sMax) = inputInfo()
#    (spaceGroup, cell, wavelength, sMin, sMax, atoms) = testInfo()
    (spaceGroup, cell, atomList) = readFile("LuFeMnO3.cif")
    atoms = deconstruct_dv(atomList.atoms, Atom)
    for atom in atoms: print atom.label, atom.occupancy, atom.multip
    wavelength = 1.5403
    sMin, sMax = getS(3, wavelength), getS(167.8, wavelength)
    backgFile = "LuFeO3 Background.BGR"
    observedFile = "Al2O3.dat"
#    sg = SpaceGroup()
#    cell = CrystalCell()
#    setSpaceGroup("81", sg)
#    setCrystalCell([3,3,5], [90,90,90], cell)
#    peakCount = 24
#    backg = LinSpline(np.array([0, 180]), np.array([100, 100]))
#    createData(sg, cell, 2, 10, 140, 0.05, [.01,-.01,1],
#               [2000-75*i for i in xrange(peakCount)], backg, "testData.txt")

    refList = hklGen(spaceGroup, cell, wavelength, sMin, sMax, True)
    reflections = list(refList.reflections.data(Reflection))
    printReflections(reflections, spaceGroup, wavelength, sMin, sMax)
#    atomList = setAtoms(*zip(*atoms))
    intensities = calcIntensity(refList, atomList, spaceGroup, wavelength)
    g = makeGaussians(reflections,[.347, -.278, .166], intensities, 1, wavelength)
    backg = LinSpline(backgFile)
#    backg = LinSpline(np.array([0,180]), np.array([100,100]))
#    exclusions = np.array([[60,66],[72.75,75.75]])
    exclusions = None
    data = np.loadtxt(observedFile, dtype=float, skiprows=1)
#    tt = data[0]
    tt = np.linspace(3, 167.75, 3296)
    observed = data.flatten()[:len(tt)]
    tt, observed = removeRange(tt, exclusions, observed)
    plotPattern(g, backg, tt, observed, twoTheta(sMin, wavelength),
                twoTheta(sMax, wavelength), .01, exclusions, labels="hkl")
    pylab.show()
    return

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    from bumps.names import Parameter, FitProblem
    problem = fit()


'''
Input data (LuFeO3):
185
5.97,11.7
2.4437
0,.314

Input data (test):
81
3,5
2
0.04,0.47

Input data (Al2O3):
167
4.7698, 13.0243
1.5403
3, 167.8
'''