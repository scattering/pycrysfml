# fswig_hklgen.py - Summer 2014
#   Generates predicted diffraction peak positions for a given unit cell and space
#   group, using the Fortran CFML library.
#   Also uses the program "bumps" to perform fitting of calculated diffraction
#   patterns to observed data.
#   
#   FortWrap/Swig Version based on original from 2013
#

from pycrysfml import *
import os
import numpy as np
import pylab
from math import floor, sqrt, log, tan, radians
from string import rstrip, ljust, rjust, center
from collections import OrderedDict
funcs = FortFuncs()

# class definitions:

# SymmetryOp Attributes:
#   rot     - rotational part of symmetry operator (3 by 3 matrix)
#   trans   - translational part of symmetry operator (vector)
class SymmetryOp(sym_oper_type):
    def __init__(self):
        sym_oper_type.__init__(self)
    
# MagSymmetryOp attributes:
#   rot     - roatational part of symmetry operator
#   phase   - phase as a fraction of 2*pi
class MagSymmetryOp(msym_oper_type):
    def __init__(self):
        msym_oper_type.__init__(self)

# WyckoffPos attributes:
#   multip      - multiplicity
#   site        - site symmetry
#   numElements - number of elements in orbit
#   origin      - origin
#   orbit       - strings containing orbit information
class WyckoffPos(wyck_pos_type):
    def __init__(self):
        wyck_pos_type.__init__(self)
        
# Wyckoff attributes:
#   numOrbits   - number of orbits
#   orbits      - list of Wyckoff position objects
class Wyckoff(wyckoff_type):
    def __init__(self):
        wyckoff_type.__init__(self)
    
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
class SpaceGroup(space_group_type):
    def __init__(self, groupName=None):
        space_group_type.__init__(self)
        if (groupName != None):
            groupName = str(groupName)
            funcs.set_spacegroup(groupName, self, None, None, None, None)

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
class CrystalCell(crystal_cell_type):
    def __init__(self, length=None, angle=None):
        crystal_cell_type.__init__(self)
        if (length != None):
            self.setCell(length, angle)
    def setCell(self, length, angle):
        funcs.set_crystal_cell(length, angle, self, None, None)     

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
class MagSymmetry(magsymm_k_type):
    def __init__(self):
        magsymm_k_type.__init__(self)    
    def setBasis(self, irrRepNum, symOpNum, vectorNum, v):
        # TODO: fix this method, sets basf array
        #c_array2 = c_float*2
        #self.basis[irrRepNum][symOpNum][vectorNum] = \
        #    (c_array2*3)(c_array2(v[0].real, v[0].imag),
        #                 c_array2(v[1].real, v[1].imag),
        #                 c_array2(v[2].real, v[2].imag))
        pass
                                                                
# Atom attributes:
#   lab --> label       - label for the atom
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
class Atom(atom_type):
    def __init__(self, *args):
        # construct an atom from a list of attributes
        atom_type.__init__(self)
        if (len(args) == 6):
            funcs.init_atom_type(self)
            self.set_atom_lab(ljust(args[0],20)) # set atom label
            self.label = ljust(args[0],20) # old field preserved for lack of fortstring get methods TODO: fix this
            self.set_atom_chemsymb(ljust(args[1], 2)) # set element
            self.element = ljust(args[1], 2) # see above
            self.set_atom_sfacsymb(ljust(self.element, 4))
            self.set_atom_x(args[2])
            self.coords = args[2] # preserved for lack of get method TODO: add fortwrap support for array return types or use fortran work-around
            self.set_atom_mult(args[3])
            self.set_atom_occ(float(args[4]))
            self.set_atom_biso(float(args[5]))
    def multip(self):
        return self.get_atom_mult()
    def occupancy(self):
        return self.get_atom_occ()
    def BIso(self):
        return self.get_atom_biso()
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
class MagAtom(matom_type):
    def __init__(self):
        matom_type.__init__(self)
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
class AtomList(atom_list_type, matom_list_type):
    def __init__(self, atoms=None, magnetic=False):
        self.magnetic = magnetic
        self.index = -1
        if magnetic:
            matom_list_type.__init__(self)
        else:
            atom_list_type.__init__(self)
        if (atoms != None):
            self.numAtoms = len(atoms)
            funcs.allocate_atom_list(numAtoms, self, None)
            self.numAtoms = numAtoms
            self.set_atom_list_natoms(numAtoms)
            
            # TODO: replace following with pycrysfml equivalent
            ## copy information from provided atom list
            #for i, atom in enumerate(self):
            #    for field in atom._fields_:
            #        setattr(atom, field[0], getattr(atoms[i], field[0]))
    def __len__(self):
        if self.magnetic:
            return self.get_matom_list_natoms()
        else:
            return self.get_atom_list_natoms()
    def __iter__(self):
        return self
    def next(self):
        self.index += 1
        if self.index == len(self):
            self.index = -1
            raise StopIteration
        return self[self.index]
    def __getitem__(self, index):
        if (index < 0): index += len(self)
        ind = intp()
        ind.assign(index)
        if self.magnetic:
            result = MagAtom()
            self.get_matom_list_element(result, ind)
            return result
        else:
            result = Atom()
            self.get_atom_list_element(result, ind)
            return result
    def __setitem__(self, index, value):
        ind = intp()
        ind.assign(index)
        if self.magnetic:
            self.set_matom_list_element(value, ind)
        else:
            self.set_atom_list_element(value, ind)

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
class Reflection(reflection_type):
    def __init__(self):
        reflection_type.__init__(self)

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
class MagReflection(magh_type):
    # can initialize this from a regular (non-magnetic) reflection
    def __init__(self, reflection=None):
        magh_type.__init__(self)
        if (reflection != None):
            self.hkl = reflection.hkl # preserved for lack of working getter
            self.set_magh_h(reflection.hkl)
            self.set_magh_mult(reflection.multip())
            self.set_magh_s(reflection.s())
    def multip(self):
        return self.get_magh_mult()
    
    def s(self):
        return self.get_magh_s()
    
# ReflectionList attributes
#   numReflections  - the number of reflections
#   reflections     - list of Reflection objects
#   magnetic        - True if this is a list of MagReflections
class ReflectionList(reflection_list_type, magh_list_type):
    # TODO: fix this... add get/set methods and rewrite python list operator methods
    #     : also add get/ set to crystal_cell_type
    def __init__(self, magnetic=False):
        self.index = -1
        if not magnetic:
            reflection_list_type.__init__(self)
        else:
            magh_list_type.__init__(self)
        self.magnetic = magnetic
    def __len__(self):
        #    return int(self.numReflections)
        if self.magnetic:
            return self.get_magh_list_nref()
        else:
            return self.get_reflection_list_nref()
    def __iter__(self):
        return self
    def next(self):
        self.index += 1
        if self.index == len(self):
            self.index = -1
            raise StopIteration
        return self[self.index]    
    def __getitem__(self, index):
        if (index < 0): index += len(self)
        ind = intp()
        ind.assign(index)
        if self.magnetic:
            result = MagReflection()
            self.get_magh_list_element(result, ind)
            return result
        else:
            result = Reflection()
            self.get_reflection_list_element(result, ind)
            return result
    def __setitem__(self, index, value):
        ind = intp()
        ind.assign(index)
        if self.magnetic:
            self.set_magh_list_element(value, ind)
        else:
            self.set_reflection_list_element(value, ind)    

# FileList: represents a Fortran file object
class FileList(file_list_type):    
    def __init__(self, filename):
        file_list_type.__init__(self)
        funcs.file_to_filelist(filename, self)

# function defs:

# readInfo: acquires cell, space group, and atomic information from a .cif,
#   .cfl, .pcr, or .shx file
def readInfo(filename):
    # read the file
    cell = CrystalCell()
    spaceGroup = SpaceGroup()
    atomList = AtomList()
    ext = filename.split(".")[len(filename.split("."))-1]
    funcs.readxtal_structure_file(filename, cell, spaceGroup, atomList, ext, None, None, None)
    return (spaceGroup, cell, atomList)

# Gaussian: represents a Gaussian function that can be evaluated at any
#   2*theta value. u, v, and w are fitting parameters.
class Gaussian(object):
    # TODO: add option for psuedo-Voigt and other peak shapes
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

# getMaxNumRef: returns the maximum number of reflections for a given cell
def getMaxNumRef(sMax, volume, sMin=0.0, multip=2):
    return funcs.get_maxnumref(sMax, volume, sMin, multip)

# hklGen: generates a list of reflections in a specified range
#   If getList is true, returns a ReflectionList object
def hklGen(spaceGroup, cell, sMin, sMax, getList=False, xtal=False):
    # Calculate the reflection positions
    maxReflections = getMaxNumRef(sMax+0.2, cell.volume, multip=spaceGroup.multip)
    # Create a reference that will be modified by calling Fortran
    reflectionCount = maxReflections
    if (not getList):
        # TODO: fix this
        #c_ReflectionArray = Reflection*max(maxReflections,1)
        #reflections = c_ReflectionArray()        
        #if xtal:
        #    # single crystal reflections (also used for magnetic structures)
        #    # This option is currently non-functioning (use genList=True instead)
        #    fn = lib.__cfml_reflections_utilities_MOD_hkl_gen_sxtal_reflection
        #    fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup),
        #                   POINTER(c_float), POINTER(c_float), POINTER(c_int),
        #                   POINTER(DV), POINTER(c_int*3), POINTER(c_int*3*2)]
        #    fn.restype = None
        #    fn(cell, spaceGroup, c_float(sMin), c_float(sMax), reflectionCount,
        #       build_struct_dv(reflections), None, None)
        #else:
        #    # powder reflections
        #    fn = lib.__cfml_reflections_utilities_MOD_hkl_uni_reflection
        #
        #    fn.argtypes = [POINTER(CrystalCell), POINTER(SpaceGroup), POINTER(c_bool),
        #                   POINTER(c_float), POINTER(c_float), c_char_p,
        #                   POINTER(c_int), POINTER(DV), POINTER(c_bool)]
        #    fn.restype = None
        #    fn(cell, spaceGroup, c_bool(True), c_float(sMin), c_float(sMax), 's',
        #       reflectionCount, build_struct_dv(reflections), c_bool(False))
        pass
    else:
        if xtal:
            reflections = ReflectionList()
            funcs.hklgen_sxtal_list(cell, spaceGroup, sMin, sMax, reflectionCount, reflections)
        else:
            reflections = ReflectionList()
            funcs.hkluni_refllist(cell, spaceGroup, True, sMin, sMax, 's', reflectionCount, reflections, False)
    if (not isinstance(reflections, ReflectionList)):
        reflections = reflections[:reflectionCount]    
    return reflections

# calcStructFact: calculates the structure factor squared for a list of planes
#   using provided atomic positions
def calcStructFact(refList, atomList, spaceGroup, wavelength):
    funcs.init_calc_strfactors(refList, atomList, spaceGroup, 'NUC', wavelength, None, 3)
    structFacts = [float() for i in xrange(refList.get_reflection_list_nref())]
    reflections = refList[:]
    for i, reflection in enumerate(reflections):
        # calculates the square of the structure factor
        funcs.calc_strfactor('P', 'NUC', i+1, float(reflection.s**2), atomList, spaceGroup, structFacts[i], None, None)
    return structFacts

# calcMagStructFact: calculates the magnetic structure factors around a list
#   of lattice reflections
# TODO: convert to pycrysfml
def calcMagStructFact(refList, atomList, symmetry, cell):
    funcs.init_mag_structure_factors(refList, atomList, symmetry, None)
    funcs.mag_structure_factors(atomList, symmetry, refList)
    # calculate the "magnetic interaction vector" (the square of which is
    #   proportional to the intensity)    
    funcs.calc_mag_interaction_vector(refList, cell)
    mivs = np.array([ref.magIntVec for ref in refList])
    return mivs


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
#    lorentz = (1+np.cos(tt)**2) / (np.sin(tt)*np.sin(tt/2))
    lorentz = (np.sin(tt)*np.sin(tt/2)) ** -1
    return sfs2 * multips * lorentz

# makeGaussians() creates a series of Gaussians to represent the powder
#   diffraction pattern
def makeGaussians(reflections, coeffs, I, scale, wavelength):
    Gaussian.scaleFactor = scale
    gaussians = [Gaussian(twoTheta(rk.s(), wavelength),
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
                uvw=[0,0,1], scale=1,
                magnetic=False, info=False, plot=False, saveFile=None,
                observedData=(None,None), labels=None):
    background = LinSpline(backgroundFile)
    sMin, sMax = getS(ttMin, wavelength), getS(ttMax, wavelength)
    if magnetic:
        #if (infoFile != None):
        #    info = readMagInfo(infoFile)
        #    if (spaceGroup == None): spaceGroup = info[0]
        #    if (cell == None): cell = info[1]
        #    if (magAtomList == None): magAtomList = info[2]
        #    if (symmetry == None): symmetry = info[3]
        #if (basisSymmetry == None): basisSymmetry = symmetry
        ## magnetic peaks
        #magRefList = satelliteGen(cell, symmetry, sMax)
        #print "length of reflection list " + str(len(magRefList))
        #magIntensities = calcIntensity(magRefList, magAtomList, basisSymmetry,
        #                               wavelength, cell, True)
        ## add in structural peaks
        #if (atomList == None): atomList = readInfo(infoFile)[2]
        #refList = hklGen(spaceGroup, cell, sMin, sMax, True)
        #intensities = calcIntensity(refList, atomList, spaceGroup, wavelength)
        #reflections = magRefList[:] + refList[:]
        #intensities = np.append(magIntensities, intensities)
        pass
    else:
        if (infoFile != None):
            info = readInfo(infoFile)
            if (spaceGroup == None): spaceGroup = info[0]
            if (cell == None): cell = info[1]
            if (atomList == None): atomList = info[2]
        refList = hklGen(spaceGroup, cell, sMin, sMax, True)
        reflections = refList[:]
        intensities = calcIntensity(refList, atomList, spaceGroup, wavelength)
    gaussians = makeGaussians(reflections, uvw, intensities, scale, wavelength)
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
        plotPattern(gaussians, background, observedData[0], observedData[1],
                    ttMin, ttMax, ttStep, exclusions, labels=labels)
        pylab.show()
    if saveFile:
        np.savetxt(saveFile, (tt, intensity), delimiter=" ")
    return (tt, intensity)

# TODO: fix print info
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
                exclusions=None, labels=None, residuals=False):
    # TODO: finish residual plot
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    ttCalc = np.linspace(ttMin, ttMax, numPoints)
    if(exclusions != None): ttCalc = removeRange(ttCalc, exclusions)
    intensity = np.array(getIntensity(gaussians, background, ttCalc))
    if (observed != None):
        if exclusions:
            ttObs, observed = removeRange(ttObs, exclusions, observed)
        pylab.plot(ttObs, observed, '-g', label="Observed",lw=1)
    pylab.plot(ttCalc, intensity, '-b', label="Calculated", lw=1)
#    pylab.fill_between(ttObs, observed, intensity, color="lightblue")
    pylab.xlabel(r"$2 \theta$")
    pylab.ylabel("Intensity")
    pylab.legend()
    if labels:
        for g in gaussians:
            if (g.center <= ttMax):
                pylab.text(g.center, np.interp(g.center, ttCalc, intensity),
                           hklString(g.hkl),
                           ha="center", va="bottom", rotation="vertical")
    if (residuals and observed):
        intensityCalc = np.array(getIntensity(gaussians, background, ttObs))
        resid = observed - intensityCalc
        pylab.plot(ttObs, resid, label="Residuals")
    return


if __name__ == '__main__':
    DATAPATH = os.path.dirname(os.path.abspath(__file__))
    infoFile = os.path.join(DATAPATH,"Al2O3.cif")
    (sG, cC, atoms) = readInfo(infoFile)
    print sG.get_space_group_numspg()
    #sG.get_space_group_spg_symb(s)
    for i in range(len(atoms)):
        print atoms[i].get_atom_occ()
    for atom in atoms:
        print atom.get_atom_occ()