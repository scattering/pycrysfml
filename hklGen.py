# hklGen.py
# Generates predicted diffraction peak positions for a given unit cell and space
#   group, using the Fortran CFML library.
# Also uses the program "bumps" to perform fitting of calculated diffraction
#   patterns to observed data.
# Created 6/4/2013
# Last edited 7/25/2013

import sys
import os
from math import floor, sqrt, log, tan, radians
from string import rstrip, ljust, rjust, center
from copy import deepcopy
from collections import OrderedDict
from ctypes import cdll, Structure, c_int, c_float, c_char, c_bool, c_char_p, \
                   c_void_p, c_ulong, addressof, sizeof, POINTER

import numpy as np
import pylab

try:
    from bumps.names import Parameter, FitProblem
except(ImportError):
    pass

# This should be the location of the CFML library
LIBFILE = "libcrysfml.so"
LIBPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), LIBFILE)
lib = cdll[LIBPATH]

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

# deconstruct_dv: converts a dope vector into an array of any type (currently
#   only implemented for 1-D arrays)
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
    _fields_ = [("rot", c_int*3*3), ("trans", c_float*3)]

# MagSymmetryOp attributes:
#   rot     - roatational part of symmetry operator
#   phase   - phase as a fraction of 2*pi
class MagSymmetryOp(Structure):
    _fields_ = [("rot", c_int*3*3), ("phase", c_float)]

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

    def __init__(self, groupName=None):
        if (groupName != None):
            fn = getattr(lib, "__cfml_crystallographic_symmetry_MOD_set_spacegroup")
            fn.argtypes = [c_char_p, POINTER(SpaceGroup), c_void_p, POINTER(c_int),
                           c_char_p, c_char_p, c_int, c_int, c_int]
            fn.restype = None
            fn(groupName, self, None, None, None, None, len(groupName), 0, 0)

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

    def __init__(self, length=None, angle=None):
        if (length != None):
            self.setCell(length, angle)
    
    def setCell(self, length, angle):
        fn = getattr(lib, "__cfml_crystal_metrics_MOD_set_crystal_cell")
        float3 = c_float*3
        fn.argtypes = [POINTER(float3), POINTER(float3), POINTER(CrystalCell),
                       POINTER(c_char), POINTER(float3)]
        fn.restype = None
        fn(float3(*length), float3(*angle), self, None, None)        

# MagSymmetry attributes:
# [corresponds to CFML MagSymm_k_type]
#   name            - name describing magnetic symmetry
#   SkType          - "Spherical_Frame" designates input Fourier coefficients
#                     should be in spherical coordinates
#   lattice         - lattice type
#   irrRepsNum      - number of irreducible representations (max 4)
#   magSymOpsNum    - number of magnetic symmetry operators per crystallographic
#                     symmetry operator (max 8)
#   centerType      - 0 = centric (-1 not at origin), 1 = acentric,
#                     2 = centric (-1 at origin)
#   magCenterType   - 1 = acentric magnetic symmetry, 2 = centric
#   numk            - number of independent propagation vectors
#   k               - propagation vectors
#   numCentVec      - number of centering lattice vectors
#   centVec         - centering lattice vectors
#   numSymOps       - number of crystallographic symmetry operators
#   multip          - multiplicity of the space group
#   numBasisFunc    - number of basis functions per representation (can be
#                     given a negative value to indicate a complex basis)
#   coeffType       - 0 = real basis function coefficients, 1 = pure imaginary
#   symOpsSymb      - symbolic form of symmetry operators
#   symOps          - crystallographic symmetry operators
#   magSymOpsSymb   - symbolic form of magnetic operators
#   magSymOps       - magnetic symmetry operators
#   basis           - coeffs of basis functions of irreducible representations
class MagSymmetry(Structure):
    _fields_ = [("name", c_char*31), ("SkType", c_char*15), ("lattice", c_char),
                ("numIrreps", c_int), ("magSymOpsNum", c_int),
                ("centerType", c_int), ("magCenterType", c_int),
                ("numk", c_int), ("k", c_float*3*12), ("numCentVec", c_int),
                ("centVec", c_float*3*4), ("numSymOps", c_int),
                ("multip", c_int), ("numBasisFunc", c_int*4),
                ("coeffType", c_int*12*4), ("symOpsSymb", c_char*40*48),
                ("symOps", SymmetryOp*48), ("magSymOpsSymb", c_char*40*48*8),
                ("magSymOps", MagSymmetryOp*48*8), ("basis", c_float*2*3*12*48*4)]

    def setBasis(self, irrRepNum, symOpNum, vectorNum, v):
        c_array2 = c_float*2
        self.basis[irrRepNum][symOpNum][vectorNum] = \
            (c_array2*3)(c_array2(v[0].real, v[0].imag),
                         c_array2(v[1].real, v[1].imag),
                         c_array2(v[2].real, v[2].imag))
                                                                
# Atom attributes:
#   label       - label for the atom
#   element     - chemical symbol of the element 
#   strFactSymb - chemical symbol used in the structure factor
#   active      - used for program control
#   atomicNum   - atomic number
#   multip      - site multiplicity
#   coords      - fractional coordinates
#   occupancy   - site occupancy
#   BIso        - isotropic temperature (Debye-Waller) factor
#   uType       - type of anisotropic thermal factor: "u_ij", "b_ij", "beta",
#                 or none
#   thType      - thermal factor type: "isotr", "aniso", or "other"
#   U           - matrix entries U11, U22, U33, U12, U13, and U23
#   UEquiv      - equivalent U
#   charge      - self-explanatory
#   moment      - magnetic moment
#   index       - index "for different purposes", whatever mysterious purposes
#                 those might be
#   numVars     - a variable completely lacking in documentation
#   freeVars    - free variables used "for different purposes" again
#   atomInfo    - string containing miscellaneous information
#   *** mult, lsq prefixes indicate multipliers and positions in the list of
#       LSQ parameters, respectively; SD indicates standard deviation
class Atom(Structure):
    _fields_ = [("label", c_char*10), ("element", c_char*2),
                ("strFactSymb", c_char*4), ("active", c_int),
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

    def __init__(self, *args):
        # construct an atom from a list of attributes
        if (len(args) == 6):
            init = getattr(lib, "__cfml_atom_typedef_MOD_init_atom_type")
            init.argtypes = [POINTER(Atom)]
            init.restype = None
            init(self)
            
            self.label = ljust(args[0], 10)
            self.element = ljust(args[1], 2)
            self.strFactSymb = ljust(self.element, 4)
            self.coords = (c_float*3)(*args[2])
            self.multip = args[3]
            self.occupancy = c_float(args[4])
            self.BIso = c_float(args[5])
#        # copy over attributes from a magnetic atom
#        for field in Atom._fields_:
#            setattr(self, field[0], getattr(magAtom, field[0]))

    def sameSite(self, other):
        # TODO: make this work for equivalent sites, not just identical ones
        # returns true if two atoms occupy the same position
        # Warning: they must be specified with identical starting coordinates
        eps = 0.001
        return all([approxEq(self.coords[i], other.coords[i], eps)
                    for i in xrange(3)])

# MagAtom attributes: same as atom, plus:
#   numkVectors     - number of propagation vectors (excluding -k)
#   irrepNum        - index of the irreducible representation to be used
#   SkReal          - real part of Fourier coefficient
#   SkRealSphere    - real part of Fourier coefficient (spherical coordinates)
#   SkIm            - imaginary part of Fourier coefficient
#   SkimSphere      - imaginary part of Fourier coefficient (spherical coords)
#   phase           - magnetic phase (fraction of 2*pi)
#   basis           - coefficients of the basis functions
class MagAtom(Structure):
    _fields_ = [("label", c_char*10), ("element", c_char*2),
                ("strFactSymb", c_char*4), ("active", c_int),
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
                ("freeVars", c_float*10), ("atomInfo", c_char*40),
                ("numkVectors", c_int), ("irrepNum", c_int*12),
                ("SkReal", c_float*3*12), ("SkRealSphere", c_float*3*12),
                ("multSkReal", c_float*3*12), ("lsqSkReal", c_int*3*12),
                ("SkIm", c_float*3*12), ("SkImSphere", c_float*3*12),
                ("multSkIm", c_float*3*12), ("lsqSkIm", c_int*3*12),
                ("phase", c_float*12), ("multPhase", c_float*12),
                ("lsqPhase", c_int*12), ("basis", c_float*12*12),
                ("multBasis", c_float*12*12), ("lsqBasis", c_int*12*12)]

    def sameSite(self, other):
        # returns true if two atoms occupy the same position
        # Warning: they must be specified with identical starting coordinates
        eps = 0.001
        return all([approxEq(self.coords[i], other.coords[i], eps)
                    for i in xrange(3)])

# AtomList attributes:
#   numAtoms    - the number of atoms
#   atoms       - a list of Atom objects
#   magnetic    - True if this is a list of MagAtoms
class AtomList(Structure):
    _fields_ = [("numAtoms", c_int), ("atoms", DV)]

    def __init__(self, atoms=None, magnetic=False):
        self.magnetic = magnetic
        if (atoms != None):
            init = getattr(lib, "__cfml_atom_typedef_MOD_allocate_atom_list")
            init.argtypes = [POINTER(c_int), POINTER(AtomList), POINTER(c_int)]
            init.restype = None
            numAtoms = len(atoms)
            init(c_int(numAtoms), self, None)
            self.numAtoms = c_int(numAtoms)

        # copy information from provided atom list
        for i, atom in enumerate(self):
            for field in atom._fields_:
                setattr(atom, field[0], getattr(atoms[i], field[0]))
#            atoms = self.atoms.data(Atom)
#            c_array3 = c_float*3
#            if not isinstance(elements[0], Atom):
#                # build atoms from scratch
#                for i in xrange(numAtoms):
#                    atoms[i].element = elements[i]
#                    atoms[i].atomicNum = c_int(atomicNums[i])
#                    atoms[i].coords = c_array3(*coords[i])
#                    atoms[i].occupancy = c_float(occupancies[i])
#                    atoms[i].BIso = c_float(BIsos[i])  
#                    atoms[i].label = elements[i]

    def __len__(self):
        return int(self.numAtoms)
    
    def __getitem__(self, index):
        if (index < 0): index += len(self)
        if (self.magnetic): dtype = MagAtom
        else: dtype = Atom
        return self.atoms.data(dtype)[index]

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

# MagReflection attributes:
# [corresponds to CFML MagH_type]
#   equalMinus      - True if k is equivalent to -k
#   multip          - multiplicity
#   knum            - index of the propagation vector (k)
#   signk           - equal to +1 for -k and -1 for +k, because somebody
#                     thought that labeling system was logical
#   s               - sin(theta)/lambda
#   magIntVecSq     - norm squared of the magnetic interaction vector
#   hkl             - reciprocal scattering vector +/- k
#   magStrFact      - magnetic structure factor
#   magIntVec       - magnetic interaction vector
#   magIntVecCart   - magnetic interaction vector (Cartesian coordinates)
class MagReflection(Structure):
    _fields_ = [("equalMinus", c_int),  ("multip", c_int), ("knum", c_int),
                ("signk", c_float), ("s", c_float), ("magIntVecSq", c_float),
                ("hkl", c_float*3), ("magStrFact", c_float*2*3),
                ("magIntVec", c_float*2*3), ("magIntVecCart", c_float*2*3)]

    # can initialize this from a regular (non-magnetic) reflection
    def __init__(self, reflection=None):
        if (reflection != None):
            self.hkl = reflection.hkl
            self.multip = reflection.multip
            self.s = reflection.s

# ReflectionList attributes
#   numReflections  - the number of reflections
#   reflections     - list of Reflection objects
#   magnetic        - True if this is a list of MagReflections
class ReflectionList(Structure):
    _fields_ = [("numReflections", c_int), ("reflections", DV)]

    def __init__(self, magnetic=False):
        self.magnetic = magnetic
    
    def __len__(self):
        return int(self.numReflections)
    
    def __getitem__(self, index):
        if (index < 0): index += len(self)
        if (self.magnetic): dtype = MagReflection
        else: dtype = Reflection
        return self.reflections.data(dtype)[index]

# FileList: represents a Fortran file object
class FileList(Structure):
    _fields_ = [("numLines", c_int), ("lines", DV)]
    
    def __init__(self, filename):
        makeList = getattr(lib, "__cfml_io_formats_MOD_file_to_filelist")
        makeList.argtypes = [c_char_p, POINTER(FileList), c_int]
        makeList.restype = None
        makeList(filename, self, len(filename))

# Gaussian: represents a Gaussian function that can be evaluated at any
#   2*theta value. u, v, and w are fitting parameters.
class Gaussian(object):
    # TODO: add option for psuedo-Voigt peak shape
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
    def __init__(self, arg1=None, arg2=None):
        if (arg1 == None):
            # create default uniform 0 background
            self.x = [0,1]
            self.y = [0,0]
        elif isSequence(arg1):
            # read in x and y coordinates from lists
            self.x = np.copy(arg1)
            self.y = np.copy(arg2)
        elif (type(arg1) == str):
            # read in x and y coordinates from a file
            self.x, self.y = np.loadtxt(arg1, dtype=float, skiprows=5, unpack=True)
        else:
            # create a uniform background with the numeric value of arg1
            self.x = [0,1]
            self.y = [arg1, arg1]

    # __call__: returns the interpolated y value at some x position
    def __call__(self, x):
        # locate the two points to interpolate between
        return np.interp(x, self.x, self.y)

    def __repr__(self):
        return "LinSpline(" + str(self.x) + ", " + str(self.y) + ")"

# approxEq: returns True if two floats are equal to within some tolerance
def approxEq(num1, num2, eps):
    return abs(num1-num2) <= eps

# isSequence: returns True if a value is a list, tuple, numpy array, etc. and 
#   False if it is a string or a non-sequence
def isSequence(x):
    return ((not hasattr(x, "strip") and hasattr(x, "__getslice__")) or
            hasattr(x, "__iter__"))

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
    try:
        for x in hkl:
            assert(x == int(x))
        return "%d %d %d" % tuple(hkl)
    except(AssertionError):
        return "%.1f %.1f %.1f" % tuple(hkl)

# calcS: calculates the sin(theta)/lambda value for a given list of planes
def calcS(cell, hkl):
    if isSequence(hkl[0]):
        # list of hkl positions
        return [calcS(cell, h) for h in hkl]
    else:
        # single hkl
        fn = lib.__cfml_reflections_utilities_MOD_hs_r
        float3 = c_float*3
        fn.argtypes = [POINTER(float3), POINTER(CrystalCell)]
        fn.restype = c_float
        return float(fn(float3(*hkl), cell))

# applySymOp: applies a symmetry operator to a given vector and normalizes the
#   resulting vector to stay within one unit cell
def applySymOp(v, symOp):
    rotMat = np.mat(symOp.rot)
    transMat = np.mat(symOp.trans)
    vMat = np.mat(v).T
    newV = rotMat * vMat + transMat
    return np.array(newV.T[0])%1

# dotProduct: returns the dot product of two complex vectors (stored as 3x2
#   Numpy arrays) or a list of two complex vectors
def dotProduct(v1, v2):
    if (not isSequence(v1[0][0]) or len(v1[0][0]) == 1):
        u1 = np.matrix([complex(v1[0][0], -v1[0][1]),
                        complex(v1[1][0], -v1[1][1]),
                        complex(v1[2][0], -v1[2][1])])
        u2 = np.matrix([[complex(v2[0][0], v2[0][1])],
                        [complex(v2[1][0], v2[1][1])],
                        [complex(v2[2][0], v2[2][1])]])
        dot = np.dot(u1, u2)
        return np.array(dot)[0][0]
    else:
        return np.array([dotProduct(u1,u2) for u1, u2 in zip(v1, v2)])

# inputInfo: requests space group and cell information, wavelength, and range
#   of interest
def inputInfo():
    # Input the space group name and create it
    groupName = raw_input("Enter space group "
                          "(HM symbol, Hall symbol, or number): ")
    spaceGroup = SpaceGroup(groupName)

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
    cell = CrystalCell(length, angle)

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

# readInfo: acquires cell, space group, and atomic information from a .cif,
#   .cfl, .pcr, or .shx file
def readInfo(filename):
    fn = lib.__cfml_io_formats_MOD_readn_set_xtal_structure_split
    fn.argtypes = [c_char_p, POINTER(CrystalCell), POINTER(SpaceGroup),
                   POINTER(AtomList), c_char_p, POINTER(c_int),
                   c_void_p, c_void_p, c_char_p, c_int, c_int, c_int]
    fn.restype = None
    cell = CrystalCell()
    spaceGroup = SpaceGroup()
    atomList = AtomList()
    ext = filename[-3:]
    fn(filename, cell, spaceGroup, atomList, ext, None, None, None,
       None, len(filename), len(ext), 0)
    return (spaceGroup, cell, atomList)

# readMagInfo: acquires cell, space group, atomic, and magnetic information
#   from a .cfl file
def readMagInfo(filename):
    fn = lib.__cfml_magnetic_symmetry_MOD_readn_set_magnetic_structure
    fn.argtypes = [POINTER(FileList), POINTER(c_int), POINTER(c_int),
                   POINTER(MagSymmetry), POINTER(AtomList), c_void_p,
                   c_void_p, POINTER(CrystalCell)]
    fn.restype = None

    info = readInfo(filename)
    spaceGroup = info[0]
    cell = info[1]
    symmetry = MagSymmetry()
    atomList = AtomList(magnetic=True)
    fileList = FileList(filename)

    print "reading information from " + filename
    fn(fileList, c_int(0), c_int(0), symmetry, atomList, None, None, cell)
    return (spaceGroup, cell, atomList, symmetry)

# readData: reads in a data file for the observed intensities and 2*theta
#   values
#   xy: 2*theta and intensity values in two-column format
#   y: data file only contains intensities for a linear 2*theta range
def readData(filename, kind="xy", skiplines=0, skipcols=0, colstep=1,
             start=None, stop=None, step=None, exclusions=None):
    if (kind == "xy"):
        tt, observed = np.loadtxt(filename, dtype=float, usecols=(0,1),
                                  skiprows=skiplines, unpack=True)
    elif (kind == "y"):
        tt = np.linspace(start, stop, round((stop-start)/float(step) + 1))
        data = np.loadtxt(filename, dtype=float, skiprows=skiplines)
        observed = data.flatten()
        observed = observed[skipcols:skipcols+colstep*len(tt):colstep]
    tt, observed = removeRange(tt, exclusions, observed)
    return (tt, observed)

# getMaxNumRef: returns the maximum number of reflections for a given cell
def getMaxNumRef(sMax, volume, sMin=0.0, multip=2):
    fn = lib.__cfml_reflections_utilities_MOD_get_maxnumref
    fn.argtypes = [POINTER(c_float), POINTER(c_float), POINTER(c_float),
                   POINTER(c_int)]
    fn.restype = c_int
    numref = fn(c_float(sMax), c_float(volume), c_float(sMin), c_int(multip))
    return numref

# hklGen: generates a list of reflections in a specified range
#   If getList is true, returns a ReflectionList object
def hklGen(spaceGroup, cell, sMin, sMax, getList=False, xtal=False):
    # Calculate the reflection positions
    maxReflections = getMaxNumRef(sMax+0.2, cell.volume, multip=spaceGroup.multip)
    # Create a reference that will be modified by calling Fortran
    reflectionCount = c_int(maxReflections)
    if (not getList):
        c_ReflectionArray = Reflection*max(maxReflections,1)
        reflections = c_ReflectionArray()        
        if xtal:
            # single crystal reflections (also used for magnetic structures)
            # This option is currently non-functioning (use genList=True instead)
            fn = lib.__cfml_reflections_utilities_MOD_hkl_gen_sxtal_reflection
            fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup),
                           POINTER(c_float), POINTER(c_float), POINTER(c_int),
                           POINTER(DV), POINTER(c_int*3), POINTER(c_int*3*2)]
            fn.restype = None
            fn(cell, spaceGroup, c_float(sMin), c_float(sMax), reflectionCount,
               build_struct_dv(reflections), None, None)
        else:
            # powder reflections
            fn = lib.__cfml_reflections_utilities_MOD_hkl_uni_reflection

            fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup), POINTER(c_bool),
                           POINTER(c_float), POINTER(c_float), c_char_p,
                           POINTER(c_int), POINTER(DV), POINTER(c_bool)]
            fn.restype = None
            fn(cell, spaceGroup, c_bool(True), c_float(sMin), c_float(sMax), 's',
               reflectionCount, build_struct_dv(reflections), c_bool(False))
    else:
        if xtal:
            fn = lib.__cfml_reflections_utilities_MOD_hkl_gen_sxtal_list
            fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup),
                           POINTER(c_float), POINTER(c_float), POINTER(c_int),
                           POINTER(ReflectionList), POINTER(c_int*3),
                           POINTER(c_int*3*2)]
            fn.restype = None
            reflections = ReflectionList()
            fn(cell, spaceGroup, c_float(sMin), c_float(sMax), reflectionCount,
               reflections, None, None)
        else:
            fn = lib.__cfml_reflections_utilities_MOD_hkl_uni_refllist
            fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup),
                           POINTER(c_bool), POINTER(c_float), POINTER(c_float),
                           c_char_p, POINTER(c_int), POINTER(ReflectionList),
                           POINTER(c_bool)]
            fn.restype = None
            reflections = ReflectionList()
            fn(cell, spaceGroup, c_bool(True), c_float(sMin), c_float(sMax), 's',
               reflectionCount, reflections, c_bool(False))
    
    if (not isinstance(reflections, ReflectionList)):
        reflections = reflections[:reflectionCount.value]    
    return reflections

# satelliteGen: generates a list of magnetic satellite reflections below a
#   maximum sin(theta)/lambda value
def satelliteGen(cell, symmetry, sMax):
    fn = lib.__cfml_magnetic_structure_factors_MOD_gen_satellites
    fn.argtypes = [POINTER(CrystalCell), POINTER(MagSymmetry), POINTER(c_float),
                   POINTER(ReflectionList), POINTER(c_bool), POINTER(c_bool),
                   POINTER(ReflectionList)]
    fn.restype = None
    refList = ReflectionList(True)
    fn(cell, symmetry, c_float(sMax), refList, c_bool(True), c_bool(True), None)
    return refList

# calcStructFact: calculates the structure factor squared for a list of planes
#   using provided atomic positions
def calcStructFact(refList, atomList, spaceGroup, wavelength):
    init = lib.__cfml_structure_factors_MOD_init_calc_strfactors
    init.argtypes = [POINTER(ReflectionList), POINTER(AtomList),
                     POINTER(SpaceGroup), c_char_p, POINTER(c_float),
                     POINTER(c_int), c_int]
    init.restype = None
    init(refList, atomList, spaceGroup, 'NUC', c_float(wavelength), None, 3)
    
    calc = lib.__cfml_structure_factors_MOD_calc_strfactor
    calc.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_float),
                     POINTER(AtomList), POINTER(SpaceGroup), POINTER(c_float),
                     POINTER(DV), POINTER(c_float*2), c_int, c_int]
    calc.restype = None
    structFacts = [c_float() for i in xrange(refList.numReflections)]
    reflections = refList[:]
    for i, reflection in enumerate(reflections):
        # calculates the square of the structure factor
        calc('P', 'NUC', c_int(i+1), c_float(reflection.s**2), atomList,
             spaceGroup, structFacts[i], None, None, 1, 3)
    # convert from c_float to float
    structFacts = [sf.value for sf in structFacts]
    return structFacts

# calcMagStructFact: calculates the magnetic structure factors around a list
#   of lattice reflections
def calcMagStructFact(refList, atomList, symmetry, cell):
    init = lib.__cfml_magnetic_structure_factors_MOD_init_mag_structure_factors
    init.argtypes = [POINTER(ReflectionList), POINTER(AtomList),
                     POINTER(MagSymmetry), POINTER(c_int)]
    init.restype = None
    init(refList, atomList, symmetry, None)

    calc = lib.__cfml_magnetic_structure_factors_MOD_mag_structure_factors
    calc.argtypes = [POINTER(AtomList), POINTER(MagSymmetry),
                     POINTER(ReflectionList)]
    calc.restype = None
    calc(atomList, symmetry, refList)

    # calculate the "magnetic interaction vector" (the square of which is
    #   proportional to the intensity)    
    calcMiv = lib.__cfml_magnetic_structure_factors_MOD_calc_mag_interaction_vector
    calcMiv.argtypes = [POINTER(ReflectionList), POINTER(CrystalCell)]
    calcMiv.restype = None
    calcMiv(refList, cell)
#    strFacts = np.array([ref.magStrFact for ref in refList])
    mivs = np.array([ref.magIntVec for ref in refList])
    return mivs

# Ye Olde Function:
#    fn = lib.__cfml_magnetic_structure_factors_MOD_calc_magnetic_strf_miv
#    fn.argtypes = [POINTER(CrystalCell), POINTER(MagSymmetry),
#                   POINTER(AtomList), POINTER(MagReflection)]
#    fn.restype = None
#
#    float3 = c_float*3
#    equiv = lib.__cfml_propagation_vectors_MOD_k_equiv_minus_k
#    equiv.argtypes = [POINTER(float3), c_char_p, c_int]
#    equiv.restype = c_bool
#    equalMinus = [equiv(float3(*symmetry.k[i]), symmetry.lattice,
#                        len(symmetry.lattice)) for i in xrange(symmetry.numk)]
#
#    magRefs = [MagReflection() for i in xrange(len(reflections)*2*symmetry.numk)]
#    for i in xrange(len(reflections)):
#        # loop through all k and -k vectors
#        for ii in xrange(-symmetry.numk, symmetry.numk+1):
#            if (ii == 0): continue
#            index = ii + symmetry.numk
#            if (ii > 0): index -= 1
#            ref = magRefs[i*2*symmetry.numk+index]
#            ref.signk = -copysign(1,ii)
#            ref.knum = abs(ii)
#            newhkl = np.array(reflections[i].hkl) - \
#                     ref.signk * np.array(symmetry.k[ref.knum-1])
#            ref.hkl = float3(*newhkl)
#            ref.s = calcS(cell, ref.hkl)
#            ref.multip = reflections[i].multip
#            ref.equalMinus = equalMinus[ref.knum-1]
#            fn(cell, symmetry, atomList, ref)
#    return magRefs

# calcIntensity: calculates the intensity for a given set of reflections,
#   based on the structure factor
def calcIntensity(refList, atomList, spaceGroup, wavelength, cell=None,
                  magnetic=False):
    # TODO: be smarter about determining whether the structure is magnetic
    # TODO: make sure magnetic phase factor is properly being taken into account
    if (refList.magnetic):
        sfs = calcMagStructFact(refList, atomList, spaceGroup, cell)
        sfs2 = np.array([np.sum(np.array(sf)**2) for sf in sfs])
    else:
        sfs2 = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength))
    multips = np.array([ref.multip for ref in refList])
    
    tt = np.radians(np.array([twoTheta(ref.s, wavelength) for ref in refList]))
    lorentz = (1+np.cos(tt)**2) / (np.sin(tt)*np.sin(tt/2))
    return sfs2 * multips * lorentz

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

# removeRange: takes in an array of 2*theta intervals and removes them from
#   consideration for data analysis, with an optional argument for removing the
#   corresponding intensities as well
def removeRange(tt, remove, intensity=None):
    if (remove == None):
        if (intensity != None): return (tt, intensity)
        else: return tt
    if (not isSequence(remove[0]) or len(remove[0]) == 1):
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

# diffPattern: generates a neutron diffraction pattern from a file containing
#   crystallographic information or from the same information generated
#   elsewhere
def diffPattern(infoFile=None, backgroundFile=None, wavelength=1.5403,
                ttMin=0, ttMax=180, ttStep=0.05, exclusions=None,
                spaceGroup=None, cell=None, atomList=None,
                symmetry=None, basisSymmetry=None, magAtomList=None,
                magnetic=False, info=False, plot=False, saveFile=None):
    background = LinSpline(backgroundFile)
    sMin, sMax = getS(ttMin, wavelength), getS(ttMax, wavelength)
    if magnetic:
        if (infoFile != None):
            info = readMagInfo(infoFile)
            if (spaceGroup == None): spaceGroup = info[0]
            if (cell == None): cell = info[1]
            if (magAtomList == None): magAtomList = info[2]
            if (symmetry == None): symmetry = info[3]
        if (basisSymmetry == None): basisSymmetry = symmetry
        # magnetic peaks
        magRefList = satelliteGen(cell, symmetry, sMax)
        print "length of reflection list " + str(len(magRefList))
        magIntensities = calcIntensity(magRefList, magAtomList, basisSymmetry,
                                       wavelength, cell, True)
        # add in structural peaks
        if (atomList == None): atomList = readInfo(infoFile)[2]
        refList = hklGen(spaceGroup, cell, sMin, sMax, True)
        intensities = calcIntensity(refList, atomList, spaceGroup, wavelength)
        reflections = magRefList[:] + refList[:]
        intensities = np.append(magIntensities, intensities)
    else:
        if (infoFile != None):
            info = readInfo(infoFile)
            if (spaceGroup == None): spaceGroup = info[0]
            if (cell == None): cell = info[1]
            if (atomList == None): atomList = info[2]
        refList = hklGen(spaceGroup, cell, sMin, sMax, True)
        reflections = refList[:]
        intensities = calcIntensity(refList, atomList, spaceGroup, wavelength)
    gaussians = makeGaussians(reflections,[0,0,1], intensities, 1, wavelength)
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    tt = np.linspace(ttMin, ttMax, numPoints)
    intensity = getIntensity(gaussians, background, tt)

    if info:
        if magnetic:
            printInfo(cell, spaceGroup, (atomList, magAtomList), (refList, magRefList),
                      wavelength, basisSymmetry)
        else:
            printInfo(cell, spaceGroup, atomList, refList, wavelength)
    if plot:
        plotPattern(gaussians, background, None, None, ttMin, ttMax, ttStep,
                    exclusions, "hkl")
        pylab.show()
    if saveFile:
        np.savetxt(saveFile, (tt, intensity), delimiter=" ")
    return (tt, intensity)

# printInfo: prints out information about the provided space group and atoms,
#   as well as the generated reflections
def printInfo(cell, spaceGroup, atomLists, refLists, wavelength, symmetry=None):
    if (isinstance(refLists, ReflectionList)):
        atomLists = (atomLists,)
        refLists = (refLists,)
    
    divider = "-" * 40
    print "Cell information (%s cell)" % rstrip(spaceGroup.xtalSystem)
    print divider
    print " a = %.3f   alpha = %.3f" % (cell.length[0], cell.angle[0])
    print " b = %.3f   beta  = %.3f" % (cell.length[1], cell.angle[1])
    print " c = %.3f   gamma = %.3f" % (cell.length[2], cell.angle[2])
    print divider
    print
    print "Space group information"
    print divider
    print "               Number: ", spaceGroup.number
    print "           H-M Symbol: ", spaceGroup.symbol
    print "          Hall Symbol: ", spaceGroup.hallSymbol
    print "       Crystal System: ", spaceGroup.xtalSystem
    print "           Laue Class: ", spaceGroup.laue
    print "          Point Group: ", spaceGroup.pointGroup
    print " General Multiplicity: ", spaceGroup.multip
    print divider
    print
    print "Atom information (%d atoms)" % len(atomLists[0])
    print divider
    atomList = atomLists[0]
    magnetic = atomList.magnetic
    label = [rstrip(atom.label) for atom in atomList]
    x, y, z = tuple(["%.3f" % atom.coords[i] for atom in atomList]
                    for i in xrange(3))
    multip = [str(atom.multip) for atom in atomList]
    occupancy = ["%.3f" % (atom.occupancy*spaceGroup.multip/atom.multip)
                 for atom in atomList]
    # Figure out what the width of each column should be
    width = OrderedDict([('label', max(len(max(label, key=len)), 5)),
                         ('x', len(max(x, key=len))),
                         ('y', len(max(y, key=len))),
                         ('z', len(max(z, key=len))),
                         ('mult', max(len(max(multip, key=len)), 4)),
                         ('occ', max(len(max(occupancy, key=len)), 3)),
                         ])
    print "%s   %s %s %s   %s  %s" % tuple([center(key, v) for key, v 
                                             in width.iteritems()])
    for i in xrange(len(atomList)):
        print "%s  (%s %s %s)  %s  %s" % (center(label[i], width["label"]),
                                          rjust(x[i], width["x"]),
                                          rjust(y[i], width["y"]),
                                          rjust(z[i], width["z"]),
                                          center(multip[i], width["mult"]),
                                          rjust(occupancy[i], width["occ"]))
    print divider
    print
    print "Reflection information (%d reflections)" % \
          sum([len(refList) for refList in refLists])
    print divider
    for atomList, refList in zip(atomLists, refLists):
        magnetic = refList.magnetic
        if magnetic: symmObject = symmetry
        else: symmObject = spaceGroup
        h, k, l = tuple([str(ref.hkl[i]) for ref in refList] for i in xrange(3))
        multip = [str(ref.multip) for ref in refList]
        tt = ["%.3f" % twoTheta(ref.s, wavelength) for ref in refList]
        intensity = ["%.3f" % I for I in calcIntensity(refList, atomList, symmObject,
                                                       wavelength, cell, magnetic)]
        # Figure out what the width of each column should be
        width = OrderedDict([('h', len(max(h, key=len))),
                             ('k', len(max(k, key=len))),
                             ('l', len(max(l, key=len))),
                             ('mult', max(len(max(multip, key=len)), 4)),
                             ('2*theta', max(len(max(tt, key=len)), 7)),
                             ('intensity', max(len(max(intensity, key=len)), 9))
                            ])
        print "  %s %s %s   %s  %s  %s" % tuple([center(key, v) for key, v 
                                                 in width.iteritems()])
        for i in xrange(len(refList)):
            print " (%s %s %s)  %s  %s  %s" % (rjust(h[i], width["h"]),
                                               rjust(k[i], width["k"]),
                                               rjust(l[i], width["l"]),
                                               center(multip[i], width["mult"]),
                                               rjust(tt[i], width["2*theta"]),
                                               rjust(intensity[i], width["intensity"]))
        print
    print divider
    print

# plotPattern: given a series of Gaussians and a background, plots the predicted
#   intensity at every 2*theta position in a specified range, as well as the
#   observed intensity everywhere on a given list of points
def plotPattern(gaussians, background, ttObs, observed, ttMin, ttMax, ttStep,
                exclusions=None, labels=None):
    # TODO: add option for residual plot
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    ttCalc = np.linspace(ttMin, ttMax, numPoints)
    if(exclusions != None): ttCalc = removeRange(ttCalc, exclusions)
    intensity = np.array(getIntensity(gaussians, background, ttCalc))
    pylab.plot(ttCalc, intensity, '-', label="Caclulated")
    if (observed != None): pylab.plot(ttObs, observed, '-', label="Observed")
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
        self.cell.setCell([a,b,c], [alpha, beta, gamma])
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
        self.cell.setCell([a,b,c], [90, beta, 90])
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
        self.cell.setCell([a,b,c], [90,90,90])
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
        self.cell.setCell([a,a,c], [90,90,90])
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
        self.cell.setCell([a,a,c], [90,90,120])
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
        self.cell.setCell([a,a,a], [90,90,90])
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
    elif (xtalSystem.lower() in ["rhombohedral", "hexagonal", "trigonal"]):
        newCell = HexagonalCell(a,c)
    elif (xtalSystem.lower() == "cubic"):
        newCell = CubicCell(a)
    return newCell

# Model: represents an object that can be used with bumps for optimization
#   purposes
class Model(object):

    def __init__(self, tt, observed, background, u, v, w,
                 wavelength, spaceGroupName, cell, atoms, exclusions=None,
                 magnetic=False, symmetry=None, newSymmetry=None):
        if (isinstance(spaceGroupName, SpaceGroup)):
            self.spaceGroup = spaceGroupName
        else:
            self.spaceGroup = SpaceGroup(spaceGroupName)
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
        self.magnetic = magnetic
        if magnetic:
            self.symmetry = symmetry
            # used for basis vector structure factor generation
            self.newSymmetry = newSymmetry
            self.atomListModel = AtomListModel(atoms, self.spaceGroup.multip,
                                               True, self.newSymmetry)            
        else:
            self.atomListModel = AtomListModel(atoms, self.spaceGroup.multip, False)
        self._set_reflections()
        self.update()

    def _set_reflections(self):
        maxLattice = self.cell.getMaxLattice()
        maxCell = CrystalCell(maxLattice[:3], maxLattice[3:])
        self.refList = hklGen(self.spaceGroup, maxCell,
                              self.sMin, self.sMax, True)
        self.reflections = self.refList[:]
        if self.magnetic:
            self.magRefList = satelliteGen(self.cell.cell, self.symmetry,
                                           self.sMax)
            self.magReflections = self.magRefList[:]

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["refList"]
        del state["magRefList"]
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
                    self.ttMin, self.ttMax, 0.01, self.exclusions, labels=None)

#    def _cache_cell_pars(self):
#        self._cell_pars = dict((k,v.value) for k,v in self.cell.items())
    def update(self):
        self.cell.update()
        self.atomListModel.update()
        
        hkls = [reflection.hkl for reflection in self.reflections]
        sList = calcS(self.cell.cell, hkls)
        ttPos = np.array([twoTheta(s, self.wavelength) for s in sList])
        # move nonexistent peaks (placed at 180) out of the way to 2*theta = -20
        ttPos[np.abs(ttPos - 180*np.ones_like(ttPos)) < 0.0001] = -20
        for i in xrange(len(self.reflections)):
            self.reflections[i].s = getS(ttPos[i], self.wavelength)
        self.intensities = calcIntensity(self.refList, self.atomListModel.atomList, 
                                         self.spaceGroup, self.wavelength)
        self.gaussians = makeGaussians(self.reflections,
                                       [self.u.value, self.v.value, self.w.value],
                                       self.intensities, self.scale.value,
                                       self.wavelength)
        if self.magnetic:
            # update magnetic reflections and add their peaks to the list of
            #   Gaussians
            hkls = [reflection.hkl for reflection in self.magReflections]
            sList = calcS(self.cell.cell, hkls)
            ttPos = np.array([twoTheta(s, self.wavelength) for s in sList])
            # move nonexistent peaks (placed at 180) out of the way to 2*theta = -20
            ttPos[np.abs(ttPos - 180*np.ones_like(ttPos)) < 0.0001] = -20
            for i in xrange(len(self.magReflections)):
                self.magReflections[i].s = getS(ttPos[i], self.wavelength)
            self.magIntensities = calcIntensity(self.magRefList,
                                                self.atomListModel.magAtomList, 
                                                self.newSymmetry, self.wavelength,
                                                self.cell.cell, True)
            self.gaussians.extend(makeGaussians(self.magReflections,
                                           [self.u.value, self.v.value, self.w.value],
                                           self.magIntensities, self.scale.value,
                                           self.wavelength))

class AtomListModel(object):
    # TODO: make occupancy constraints automatic

    def __init__(self, atoms, sgmultip, magnetic=False, symmetry=None):
        self.sgmultip = sgmultip
        self.magnetic = magnetic
        self.symmetry = symmetry
        self._rebuild_object(atoms)

    def _rebuild_object(self, atoms):
        if (not self.magnetic):
            # one list of atoms
            if (isinstance(atoms, AtomList)):
                self.atomList = atoms
                self.atoms = atoms[:]
            else:
                self.atomList = AtomList(atoms)
                self.atoms = atoms
        else:
            # a tuple containing a list of atoms and a list of magnetic atoms
            if (isinstance(atoms[0], AtomList)):
                self.atomList = atoms[0]
                self.magAtomList = atoms[1]
                self.atoms = self.atomList[:]
                self.magAtoms = self.magAtomList[:]
            else:
                self.atomList = AtomList(atoms[0])
                self.magAtomList = AtomList(atoms[1], magnetic=True)
                self.atoms = atoms[0]
                self.magAtoms = atoms[1]

#        modelAtoms = deepcopy(self.atoms)
#        for magAtom in self.magAtoms:
#            print >>sys.stderr, magAtom.label
#            for atom in modelAtoms:
#                if (magAtom == atom):
#                    print >>sys.stderr, magAtom.label + " = " + atom.label
        self.atomModels = [AtomModel(atom, self.sgmultip) for atom in self.atoms]
        if self.magnetic:
            # correct atom models to include magnetic atoms
            for magAtom in self.magAtoms:
                for model in self.atomModels:
                    if (magAtom.label.rstrip() == model.atom.label.rstrip() and \
                        magAtom.sameSite(model.atom)):
                        model.addMagAtom(magAtom, self.symmetry)
        self.modelsDict = dict([(am.atom.label, am) for am in self.atomModels])
#        print >>sys.stderr, "atom models: ", [(am.atom.label, am.magnetic)
#                                              for am in self.atomModels]

    def __getstate__(self):
        state = self.atoms, self.sgmultip, self.magnetic
        return state

    def __setstate__(self, state):
        self.atoms, self.sgmultip, self.magnetic = state
        self._rebuild_object(self.atoms)
        
    def parameters(self):
        params = dict(zip([atom.label for atom in self.atoms],
                          [am.parameters() for am in self.atomModels]))
#        if self.magnetic:
#            params.update(zip([atom.label for atom in self.magAtoms],
#                              [am.parameters() for am in self.magAtomModels]))
        return params

    def update(self):
#        print >>sys.stderr, len(self.parameters())
#        print >>sys.stderr, "label: |" + self.atoms[0].label + "|"        
#         if len(self.parameters()) == 1:
#            print >>sys.stderr, self.parameters()
        for atomModel in self.atomModels:
            atomModel.update()
#        if self.magnetic:
#            for atomModel in self.magAtomModels:
#                atomModel.update()

    def __len__(self):
        return len(self.modelsDict)
    def __getitem__(self, key):
        return self.modelsDict[key]

class AtomModel(object):
    # TODO: implement anisotropic thermal factors

    def __init__(self, atom, sgmultip):
#        if isSequence(atoms):
            # changes properties of both a regular atom and a magnetic atom
            #   (they should correspond to the same atomic position)
#            self.atom = atoms[0]
#            self.magAtom = atoms[1]
#            self.magAtom.label = rstrip(self.magAtom.label)
#            self.magnetic = True
#        else:
        self.atom = atom
        self.magnetic = False
        self.sgmultip = sgmultip
        self.atom.label = rstrip(self.atom.label)
        self.B = Parameter(self.atom.BIso, name=self.atom.label + " B")
        occ = self.atom.occupancy / self.atom.multip * self.sgmultip
        self.occ = Parameter(occ, name=self.atom.label + " occ")
        self.x = Parameter(self.atom.coords[0], name=self.atom.label + " x")
        self.y = Parameter(self.atom.coords[1], name=self.atom.label + " y")
        self.z = Parameter(self.atom.coords[2], name=self.atom.label + " z")
    
    def addMagAtom(self, magAtom, symmetry):
        # add a secondary magnetic atom object to the model
        self.magAtom = magAtom
        self.magAtom.label = rstrip(self.magAtom.label)
        self.symmetry = symmetry
        self.magnetic = True
        self.numVectors = self.symmetry.numBasisFunc[self.magAtom.irrepNum[0]-1]
        self.coeffs = [None] * self.numVectors
        for i in xrange(self.numVectors):
            self.coeffs[i] = Parameter(self.magAtom.basis[0][i],
                                       name=self.atom.label + " C" + str(i))
#        self.phase = Parameter(0, name=self.atom.label + " phase")
    
    def parameters(self):
        params = {self.B.name: self.B, self.occ.name: self.occ,
                  self.x.name: self.x, self.y.name: self.y, self.z.name: self.z}
        if self.magnetic:
            params.update([(coeff.name, coeff) for coeff in self.coeffs])
#            params.update([(self.phase.name, self.phase)])
        return params

    def update(self):
        self.atom.BIso = self.B.value
        occ = self.occ.value * self.atom.multip / self.sgmultip
        self.atom.occupancy = occ
        self.atom.coords = (c_float*3)(self.x, self.y, self.z)
        
        if self.magnetic:
            self.magAtom.BIso = self.B.value            
            self.magAtom.occupancy = occ        
            self.magAtom.coords = (c_float*3)(self.x, self.y, self.z)
            for i in xrange(self.numVectors):
                self.magAtom.basis[0][i] = self.coeffs[i].value
#                print >>sys.stderr, self.magAtom.label, self.magAtom.basis[0][i]
#            self.magAtom.phase[0] = self.phase.value

# testStructFact: tests the structure factor calculations
def testStructFact():
    print "starting test"
    spaceGroup, cell, atomList = readInfo("Data/Al2O3.cif")
    wavelength = 1.5403
    sMin, sMax = getS(3, wavelength), getS(167.8, wavelength)
    refList = hklGen(spaceGroup, cell, sMin, sMax, True)
    print str(refList.numReflections) + " reflections generated"
    sf = calcStructFact(refList, atomList, spaceGroup, wavelength)
    print sf
    for i in xrange(len(sf)):
        print hklString(refList[i].hkl[:]), sf[i], sqrt(sf[i])
    out = np.array([[ref.hkl[0] for ref in refList],
                    [ref.hkl[1] for ref in refList],
                    [ref.hkl[2] for ref in refList],
                    [sf[i] for i in xrange(len(sf))]])
    out = out.transpose()
    np.savetxt("structure factors.dat", out)
    print "test complete"

# testMagStructFact: tests the magnetic structure factor calculations
def testMagStructFact():
    print "starting test"
    diffPattern(infoFile="Data/hobanio.cfl", backgroundFile=50, wavelength=2.524,
                ttMin=0, ttMax=90, magnetic=True, info=True, plot=True)
    print "test complete"

# testMagBasis: tests the calculation of a magnetic diffraction pattern using
#   the basis of an irreducible representation
def testMagBasis():
    # Applicable variables (Fortran name) [Array dims]:
    #   MagAtom.numkVectors (nvk) - number of k vectors (max 12)
    #   MagAtom.irrepNum (imat) [12] - index of the irrep for each k vector
    #   MagAtom.basis (cbas) [12,12] - coefficients of basis vectors
    #   MagSymmetry.numIrreps (nirreps) - number of irreducible reps (max 4)
    #   MagSymmetry.numSymOps (numops) - number of symmetry operators (max 48)
    #   MagSymmetry.numBasisFunc (nbas) [4] - number of basis vectors per irrep
    #   MagSymmetry.basis (basf) [4,48,12,3,2] - basis vectors

    print "starting test"
    backgFile = 50
    observedFile = "Data/hobanio.dat"
    infoFile = "Data/Magnetic_Tests/2atom.cfl"
    spaceGroup, crystalCell, magAtomList, symmetry = readMagInfo(infoFile)
    atomList = readInfo(infoFile)[2]

    spaceGroup.xtalSystem = spaceGroup.xtalSystem.rstrip()
    wavelength = 1.5403
    ttMin = 10
    ttMax = 100
    exclusions = None

    basisSymmetry = deepcopy(symmetry)
    if (basisSymmetry.centerType == 2):
        # change to acentric; correct the number of symmetry operators
        basisSymmetry.centerType = 1
        basisSymmetry.numSymOps *= 2
    basisSymmetry.numIrreps = 1
    basisSymmetry.numBasisFunc[0] = 1
    basisSymmetry.numBasisFunc[1] = 1
    c_array2 = c_float*2
    c_array23 = c_float*2*3
    zero = c_array2()
#    print >>sys.stderr, "Number of symmetry operators: ", symmetry.numSymOps
#    print >>sys.stderr, "Symmetry operators: ", np.array([np.array(so.rot, dtype=float)
#                                                 for so in symmetry.symOps])
#    print >>sys.stderr, "Number of symmetry operators (space group): ", \
#                        spaceGroup.numOps
#    print >>sys.stderr, "Symmetry operators (space group): ", \
#                        np.array([np.array(so.rot, dtype=float) for so in spaceGroup.symmetryOps])
    for i in xrange(symmetry.numSymOps):
        # Basis 1 (from SARAh)
        basisSymmetry.basis[0][i][0] = c_array23(c_array2(1.0, 0.0), zero, zero)
#        basisSymmetry.basis[0][i][1] = c_array23(zero, c_array2(2.0, 0.0), zero)
#        basisSymmetry.basis[0][i][2] = c_array23(zero, zero, c_array2(2.0, 0.0))
        # Basis 2 (from SARAh)

    for magAtom in magAtomList:
        magAtom.numkVectors = 1
        if (magAtom.label.rstrip().lower() == "fe"):
            magAtom.irrepNum[0] = 1
            magAtom.basis[0][0] = 1
            magAtom.basis[0][1] = 0
            magAtom.basis[0][2] = 0
        else:
            magAtom.irrepNum[0] = 1
            magAtom.basis[0][0] = 1
#            magAtom.basis[0][1] = 0
#            magAtom.basis[0][2] = 0
#        print >>sys.stderr, magAtom.label.lower(), \
#                np.array(basisSymmetry.basis[magAtom.irrepNum[0]-1][0][1])

    diffPattern(infoFile=infoFile, backgroundFile=0, wavelength=1.5403,
                ttMin=ttMin, ttMax=ttMax, ttStep=0.05, exclusions=None,
                basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                magnetic=True, info=True, plot=False)
    return

    tt, observed = readData(observedFile, kind="y", skiplines=5, skipcols=1,
                            colstep=2, start=ttMin, stop=ttMax, step=0.2)
    backg = LinSpline(backgFile)
    cell = makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
    cell.b.pm(0.5)
    cell.c.pm(0.5)
    cell.alpha.pm(0)
    cell.beta.pm(0)
    cell.gamma.pm(0)
    m = Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
              (atomList, magAtomList), exclusions, magnetic=True,
              symmetry=symmetry, newSymmetry=basisSymmetry)
    m.u.range(0,2)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.scale.range(0,100)
    for atomModel in m.atomListModel.atomModels:
        atomModel.B.range(0, 10)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
        if (atomModel.magnetic):
            # fit coefficients of basis vectors
            for coeff in atomModel.coeffs:
                coeff.range(-10,10)
    M = FitProblem(m)
    M.model_update()
    return M
#    wavelength = 2.524
#    ttMin = 0
#    ttMax = 90
#    diffPattern(infoFile="Data/hobanio.cfl", backgroundFile=50, wavelength=wavelength,
#                symmetry=symmetry, magAtomList=magAtomList, ttMin=ttMin, ttMax=ttMax,
#                magnetic=True, info=True, plot=True)
#    symmetry.centerType = 2
#    print "test complete"

def main():
#    testMagStructFact()
#    testMagBasis()
    pass

if __name__ == "__main__":
    # program run normally
    main()
