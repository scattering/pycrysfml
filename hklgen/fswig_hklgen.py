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
try:
    import matplotlib.pylab as pylab
except:
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
    @property
    def number(self):
        return self.get_space_group_numspg()    
    @property
    def multip(self):
        return self.get_space_group_multip()
    @multip.setter
    def multip(self, value):
        self.set_space_group_multip(value)
    @property
    def xtalSystem(self):
        return getSpaceGroup_crystalsys(self)
    @xtalSystem.setter
    def xtalSystem(self, value):
        self.set_space_group_crystalsys(value)

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
    def length(self):
        LVec = FloatVector([0 for i in range(3)])
        self.get_crystal_cell_cell(LVec)
        return LVec
    def angle(self):
        AVec = FloatVector([0 for i in range(3)])
        self.get_crystal_cell_ang(AVec)
        return AVec
    @property
    def volume(self):
        return self.get_crystal_cell_cellvol()
    @volume.setter
    def volume(self, value):
        self.set_crystal_cell_cellvol(value)
    def setCell(self, length, angle):
        funcs.set_crystal_cell(FloatVector(length), FloatVector(angle), self, None, None)     

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
    def centerType(self):
        return self.get_magsymm_k_centred()
    def setCenterType(self, value):
        self.set_magsymm_k_centred(value)
    def numIrreps(self):
        return self.get_magsymm_k_nirreps()
    def setNumIrreps(self, value):
        self.set_magsymm_k_nirreps(value)
    def numBasisFunc(self):
        result = IntVector([0 for i in range(4)])
        self.get_magsymm_k_nbas(result)
        return list(result)
    def setNumBasisFunc(self, value):
        self.set_magsymm_k_nbas(IntVector(value))
    def setNumBasisFunc_ind(self, ind, value):
        result = self.numBasisFunc()
        result[ind] = value
        self.setNumBasisFunc(result)
    def getBasis(self, irrRepNum, symOpNum, vectorNum):
        result = FloatVector([0 for i in range(6)])
        self.get_basis_element(irrRepNum, symOpNum, vectorNum, result)
        #self.get_basis_element(vectorNum, symOpNum, irrRepNum, result)
        return list(result)
    def setBasis(self, irrRepNum, symOpNum, vectorNum, v):
        # TODO: fix this method, sets basf array
        #c_array2 = c_float*2
        #self.basis[irrRepNum][symOpNum][vectorNum] = \
        #    (c_array2*3)(c_array2(v[0].real, v[0].imag),
        #                 c_array2(v[1].real, v[1].imag),
        #                 c_array2(v[2].real, v[2].imag))
        ri_comps = []
        for num in v:
            ri_comps.append(num.real)
            ri_comps.append(num.imag)
        self.set_basis_element(irrRepNum, symOpNum, vectorNum, FloatVector(ri_comps))
        #self.set_basis_element(vectorNum, symOpNum, irrRepNum, FloatVector(ri_comps))
    def kvec(self, index):
        k = FloatVector([0,0,0])
        self.get_kvector(k, index)
        return np.array(k)
    @property
    def lattice(self):
        return getMagsymmK_latt(self)
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
            self.set_atom_chemsymb(ljust(args[1], 2)) # set element
            self.set_atom_sfacsymb(ljust(self.element, 4))
            self.set_atom_x(FloatVector(args[2]))
            self.set_atom_mult(args[3])
            self.set_atom_occ(float(args[4]))
            self.set_atom_biso(float(args[5]))
    def coords(self):
        CVec = FloatVector([0 for i in range(3)])
        self.get_atom_x(CVec)
        return list(CVec)
    def setCoords(self, value):
        self.set_atom_x(FloatVector(value))
    def multip(self):
        return self.get_atom_mult()
    def occupancy(self):
        return self.get_atom_occ()
    def setOccupancy(self, value):
        self.set_atom_occ(value)
    def BIso(self):
        return self.get_atom_biso()
    def setBIso(self, value):
        self.set_atom_biso(value)
    def label(self):
        return getAtom_lab(self)
    def setLabel(self, label):
        self.set_atom_lab(label)
    def sameSite(self, other):
        # TODO: make this work for equivalent sites, not just identical ones
        # returns true if two atoms occupy the same position
        # Warning: they must be specified with identical starting coordinates
        eps = 0.001
        return all([approxEq(self.coords()[i], other.coords()[i], eps)
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
    def coords(self):
            CVec = FloatVector([0 for i in range(3)])
            self.get_matom_x(CVec)
            return list(CVec)    
    def sameSite(self, other):
        # returns true if two atoms occupy the same position
        # Warning: they must be specified with identical starting coordinates
        eps = 0.001
        return all([approxEq(self.coords()[i], other.coords()[i], eps)
                    for i in xrange(3)])
    def setCoords(self, value):
        self.set_matom_x(FloatVector(value))    
    def multip(self):
        return self.get_matom_mult()
    def occupancy(self):
        return self.get_matom_occ()
    def setOccupancy(self, value):
        self.set_matom_occ(value)    
    def BIso(self):
        return self.get_matom_biso()
    def setBIso(self, value):
        self.set_matom_biso(value)
    def irrepNum(self):
        result = IntVector([0 for i in range(12)])
        self.get_matom_imat(result)
        return list(result)
    def basis(self):
        result = [[self.get_matom_basis_element(i, j) for j in range(12)] for i in range(12)]
        return list(result)
    def setBasis(self, value):
        for i in range(12):
            for j in range(12):
                self.set_matom_basis_element(i, j, value[i][j])
    def setBasis_ind(self, i, j, k):
        b = self.basis()
        b[j][i] = k
        self.setBasis(b)
    def setIrrepNum(self, value):
        self.set_matom_imat(IntVector(value))
    def setIrrepNum_ind(self, i, k):
        value = self.irrepNum()
        value[i] = k
        self.setIrrepNum(IntVector(value))        
    def numkVectors(self):
        return self.get_matom_nvk()
    def setNumkVectors(self, value):
        self.set_matom_nvk(value)
    def label(self):
        return getMatom_lab(self)
    def setLabel(self, label):
        self.set_matom_lab(label)
    @property
    def phase(self):
        ph = FloatVector([0 for i in range(12)])
        self.get_matom_mphas(ph)
        return list(ph)
    @property
    def phase_multiplier(self):
        phm = FloatVector([0 for i in range(12)])
        self.get_matom_mmphas(phm)
        return list(phm)
    def set_phase(self, value):
        #ph = FloatVector([value, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        ph = FloatVector([value for i in range(12)])
        self.set_matom_mphas(ph)
    def debug(self):
        funcs.printbasis(self)
# AtomList attributes:
#   numAtoms    - the number of atoms
#   atoms       - a list of Atom objects
#   magnetic    - True if this is a list of MagAtoms
# Requires manual resetting of elements when dealing with arrays for some reason. e.g.
#atm = atoms[0]
#atm.setCoords([0,0,0.15])
#atoms[0] = atm
# instead of:
#atoms[0].setCoords([0,0,0.15])
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
            funcs.allocate_atom_list(len(atoms), self, None)
            #self.numAtoms = numAtoms
            self.set_atom_list_natoms(len(atoms))
            
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
        if isinstance(index, int):
            if (index < 0): index += len(self)
            if self.magnetic:
                result = MagAtom()
                self.get_matom_list_element(result, index)
                return result
            else:
                result = Atom()
                self.get_atom_list_element(result, index)
                return result
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self))    # index is a slice
            L = []
            for i in range(start, stop, step):
                L.append(self.__getitem__(i))
            return L
        else:
            raise TypeError("index must be int or slice")        
    def __setitem__(self, index, value):
        if self.magnetic:
            self.set_matom_list_element(value, index)
            #print 'setting index', index, 'to', value
        else:
            self.set_atom_list_element(value, index)

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
    @property
    def hkl(self):
        hklVec = IntVector([0 for i in range(3)])
        self.get_reflection_h(hklVec)
        return list(hklVec)
    @hkl.setter
    def hkl(self, value):
        self.set_reflection_h(FloatVector(value))
    @property
    def s(self):
        return self.get_reflection_s()
    @s.setter
    def s(self, value):
        self.set_reflection_s(value)
    @property
    def multip(self):
        return self.get_reflection_mult()
    @multip.setter
    def multip(self, value):
        self.set_reflection_mult(value)
    

# MagReflection attributes:
# [corresponds to CFML MagH_type]
#   keq             - True if k is equivalent to -k
#   multip          - multiplicity
#   knum            - index of the propagation vector (k)
#   signk           - equal to +1 for -k and -1 for +k, because somebody
#                     thought that labeling system was logical
#   s               - sin(theta)/lambda
#   sqMiV     - norm squared of the magnetic interaction vector
#   hkl             - reciprocal scattering vector +/- k
#   magStrFact      - magnetic structure factor
#   magIntVec       - magnetic interaction vector
#   magIntVecCart   - magnetic interaction vector (Cartesian coordinates)
class MagReflection(magh_type):
    # can initialize this from a regular (non-magnetic) reflection
    def __init__(self, reflection=None):
        magh_type.__init__(self)
        if (reflection != None):
            self.hkl = reflection.hkl
            self.multip = reflection.multip
            self.set_magh_s(reflection.get_reflection_s())
    def __eq__(self, other):
        return self.s == other.s and self.sqMiV == other.sqMiV
    @property
    def multip(self):
        return self.get_magh_mult()
    @multip.setter
    def multip(self, value):
        self.set_magh_mult(value)
    @property
    def keq(self):
        return self.get_magh_keqv_minus()
    def magStrFact(self):
        rVec = FloatVector([0 for i in range(6)])
        self.get_msf(rVec)
        result = []
        for i in range(0, len(rVec), 2):
            result.append(complex(rVec[i],rVec[i+1]))
        return result
    def magIntVec(self):
        rVec = FloatVector([0 for i in range(6)])
        self.get_miv(rVec)
        result = []
        for i in range(0, len(rVec), 2):
            result.append(complex(rVec[i],rVec[i+1]))
        return result
    def setMagIntVec(self, value):
        rVec = []
        for num in value:
            rVec.append(num.real, num.imag)
        self.setMagIntVec(FloatVector(rVec))
    @property
    def s(self):
        return self.get_magh_s()
    @s.setter
    def s(self, value):
        self.set_magh_s(value)
    @property
    def hkl(self):
        hklVec = FloatVector([0,0,0])
        self.get_magh_h(hklVec)
        return list(hklVec)
    @hkl.setter
    def hkl(self, value):
        self.set_magh_h(FloatVector(value))
    @property
    def sqMiV(self):
        return self.get_magh_sqmiv()
    
# ReflectionList attributes
#   numReflections  - the number of reflections
#   reflections     - list of Reflection objects
#   magnetic        - True if this is a list of MagReflections
class ReflectionList(reflection_list_type, magh_list_type):
    def __init__(self, cast_list=None, magnetic=None):
        self.index = -1
        if cast_list != None:
            # copy constructor
            if isinstance(cast_list[0], magh_type):
                self.magnetic = True
                magh_list_type.__init__(self)
                self.set_magh_list_nref(len(cast_list))
                funcs.alloc_mhlist_array(self)
            elif isinstance(cast_list[0], reflection_type):
                self.magnetic = False
                reflection_list_type.__init__(self)
                self.set_reflection_list_nref(len(cast_list))
                funcs.alloc_refllist_array(self)
            for i in range(len(cast_list)):
                ref = cast_list[i]
                self[i]=ref
        else:
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
        if isinstance(index, int):
            if (index < 0): index += len(self)
            if self.magnetic:
                result = MagReflection()
                self.get_magh_list_element(result, index)
                return result
            else:
                result = Reflection()
                self.get_reflection_list_element(result, index)
                return result
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(self))    # index is a slice
            L = []
            for i in range(start, stop, step):
                L.append(self.__getitem__(i))
            return L
            
        else:
            raise TypeError("index must be int or slice")
    def __setitem__(self, index, value):
        ind = intp()
        ind.assign(index)
        if self.magnetic:
            self.set_magh_list_element(value, index)
        else:
            self.set_reflection_list_element(value, ind)
    @property
    def nref(self):
        if self.magnetic:
            return self.get_magh_list_nref()
        else:
            return self.get_reflection_list_nref()
    @nref.setter
    def nref(self, value):
        if self.magnetic:
            self.set_magh_list_nref(value)
        else:
            self.set_reflection_list_nref(value)

# FileList: represents a Fortran file object
class FileList(file_list_type):    
    def __init__(self, filename):
        file_list_type.__init__(self)
        funcs.file_to_filelist(filename, self)

# function defs:
# disposable pointer conversions for parameters
def int_to_p(i):
    point = intp()
    point.assign(i)
    return point
def float_to_p(f):
    point = floatp()
    point.assign(f)
    return point
# readInfo: acquires cell, space group, and atomic information from a .cif,
#   .cfl, .pcr, or .shx file
def readInfo(filename):
    # read the file
    cell = CrystalCell()
    spaceGroup = SpaceGroup()
    atomList = AtomList()
    ext = filename.split(".")[-1]
    funcs.readxtal_structure_file(filename, cell, spaceGroup, atomList, ext, None, None, None)
    return (spaceGroup, cell, atomList)

# Peak: represents a Peak shape function that can be evaluated at any
#   2*theta value. u, v, and w are fitting parameters.
class Peak(object):
    scaleFactor = 1    
    def __init__(self, center, u, v, w, I, hkl=[None, None, None], shape="Gaussian", eta=0):
        self.center = center    # 2*theta position
        self.u = u
        self.v = v
        self.w = w
        self.I = I
        self.shape = shape
        self.eta = eta
        try:
            self.H = sqrt(u*(tan(radians(center/2))**2)
                          + v*tan(radians(center/2)) + w)
            self.scale = self.I * Peak.scaleFactor
        except ValueError:
            self.H = 0
            self.scale = 0
        self.hkl = hkl

    # __call__: returns the value of the Peak at some 2*theta positions
    def __call__(self, x):
        if self.shape.lower() == 'gaussian':
            return [self.scale * funcs.calcgaussian(value-self.center, self.H) for value in x]
        elif self.shape.lower() == 'pseudovoigt':
            return [self.scale * funcs.calcpseudovoigt(value-self.center, self.H, float(self.eta)) for value in x]
        elif self.shape.lower() == 'lorentzian':
            return [self.scale * funcs.calclorentzian(value-self.center, self.H) for value in x]
        else:
            print "Unsupported Peak-shape: "+shape+"\n"
            quit()
    def add(self, v, x):
        # only add to nearby 2*theta positions
        idx = (x>self.center-self.H*3) & (x<self.center+self.H*3)
        try:
            v[idx] += self.__call__(np.array(x)[idx])
        except:
            pass
# LinSpline: represents a linear spline function to be used for the background
class LinSpline(object):
    def __init__(self, arg1=None, arg2=None, fn="LinSpline"):
        self.fn = fn
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
        if self.fn.lower() == 'linspline':
            # locate the two points to interpolate between
            return np.interp(x, self.x, self.y)
        elif self.fn.lower() == 'polynomial':
            pass
    #def __add__(self, other):
        #for i in range(len(self.y)):
            #self.y[i] += other
    def __repr__(self):
        return "LinSpline(" + str(self.x) + ", " + str(self.y) + ")"
# readMagInfo: acquires cell, space group, atomic, and magnetic information
#   from a .cfl file
def readMagInfo(filename):
    info = readInfo(filename)
    spaceGroup = info[0]
    cell = info[1]
    symmetry = MagSymmetry()
    atomList = AtomList(magnetic=True)
    fileList = FileList(filename)
    zerop = intp()
    zp = intp()
    zerop.assign(0)
    zp.assign(0)
    funcs.read_mag_cfl_file(fileList, zerop, zp, symmetry, atomList, None, None, cell)
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
def readIllData(filename, instrument, bacfile):
    data = diffraction_pattern_type()
    funcs.read_ill_data(filename, data, instrument)
    tt = FloatVector([0]*data.get_diffraction_pattern_npts())
    observed = FloatVector([0]*data.get_diffraction_pattern_npts())
    sigma = FloatVector([0]*data.get_diffraction_pattern_npts())
    #bkg = FloatVector([0]*data.get_diffraction_pattern_npts())
    #print "points",data.get_diffraction_pattern_npts()
    #return None,None
    data.get_diffraction_pattern_x(tt)
    data.get_diffraction_pattern_y(observed)
    data.get_diffraction_pattern_sigma(sigma)
    #error from ILL instruments is given as sigma squared according to crysfml documentation
    sigma = np.sqrt(np.array(sigma))
    #print data.get_diffraction_pattern_step()
    funcs.read_background_file(bacfile, "POL", data)
    #data.get_diffraction_pattern_bgr(bkg)
    #print list(bkg)
    #print min(observed)
    return list(tt), list(observed), list(sigma)
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
        return funcs.hkls_r(FloatVector(list(hkl)), cell)

# applySymOp: applies a symmetry operator to a given vector and normalizes the
#   resulting vector to stay within one unit cell
def applySymOp(v, symOp):
    rotMat = np.mat(symOp.get_sym_oper_rot())
    vMat = np.mat(v).T
    newV = np.array((rotMat * vMat).T) + np.array(symOp.get_sym_oper_tr())
    return np.array(newV)%1

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

# getMaxNumRef: returns the maximum number of reflections for a given cell
def getMaxNumRef(sMax, volume, sMin=0.0, multip=2):
    sMin_p = floatp()
    sMin_p.assign(sMin)
    multip_p = intp()
    multip_p.assign(multip)
    return funcs.get_maxnumref(sMax, volume, sMin_p, multip_p)

# hklGen: generates a list of reflections in a specified range
#   If getList is true, returns a ReflectionList object
def hklGen(spaceGroup, cell, sMin, sMax, getList=True, xtal=False):
    # Calculate the reflection positions
    maxReflections = getMaxNumRef(sMax+0.2, cell.volume, multip=spaceGroup.multip)
    # Create a reference that will be modified by calling Fortran
    reflectionCount = maxReflections
    if (not getList):
        # Should just pass back a ReflectionList() type, no reason to have python list of individual reflections
        # if else preserved for compatibility with ctypes model files which require the getList flag
        # TODO: Remove this and fix model files
        pass
    else:
        if xtal:
            reflections = ReflectionList()
            funcs.hklgen_sxtal_list(cell, spaceGroup, np.asscalar(sMin), np.asscalar(sMax), int_to_p(reflectionCount), reflections)
        else:
            reflections = ReflectionList()
            funcs.hkluni_refllist(cell, spaceGroup, True, np.asscalar(sMin), np.asscalar(sMax), 's', reflectionCount, reflections)
    if (not isinstance(reflections, ReflectionList)):
        reflections = reflections[:reflectionCount]    
    return reflections
# satelliteGen: generates a list of magnetic satellite reflections below a
#   maximum sin(theta)/lambda value
def satelliteGen(cell, symmetry, sMax, hkls=None):
    refList = ReflectionList(magnetic=True)
    if hkls != None:
        funcs.gen_satellites(cell, symmetry, sMax, refList, int_to_p(1), None, hkls)
    else:
        funcs.gen_satellites(cell, symmetry, sMax, refList, int_to_p(1))
    return refList
# satelliteGen_python: python implementation used for debugging
def satelliteGen_python(cell, sMax, hkls, kvec=[0.5,0,0.5]):
    kvec = [0,0.,0.]
    kvec2 = [0.,0,0.1651010]
    refList = ReflectionList(magnetic=True)#[0 for i in range(len(hkls)*2+2)]#ReflectionList(True)
    hkls = []
    for x in range(-7,7):
        for y in range(-7,7):
            for z in range(-7,7):
                #if x+y+z in [2*n for n in range(-7,7)]:
                refl = Reflection()
                refl.set_reflection_h(IntVector([x,y,z]))
                hkls.append(refl)
    print(len(hkls)), "python hkls"
    refList.set_magh_list_nref(len(hkls)*2+2)
    funcs.alloc_mhlist_array(refList)
    i = 0
    for reflection in hkls:
        hkl = reflection.hkl
        hkp = list(np.add(hkl, kvec))
        hkm = list(np.add(hkl,np.negative(kvec)))
        sp = calcS(cell, hkp)
        sm = calcS(cell, hkm)
        hkp2 = list(np.add(hkl, kvec2))
        hkm2 = list(np.add(hkl, np.negative(kvec2)))
        sp2 = calcS(cell, hkp2)
        sm2 = calcS(cell, hkm2)
       # if sum(hkl) in [2*n for n in range(-3,3)]:
        if sp <= sMax:
            mref = MagReflection(reflection)
            mref.set_magh_h(FloatVector(hkp))
            mref.set_magh_s(sp)
            mref.set_magh_num_k(1)
            mref.set_magh_signp(-1.0)
            mref.set_magh_keqv_minus(False)
            mref.set_magh_mult(2)
            refList[i] = mref
            i += 1
        if sm <= sMax:
            mref = MagReflection(reflection)
            mref.set_magh_h(FloatVector(hkm))
            mref.set_magh_s(sm)
            mref.set_magh_num_k(1)
            mref.set_magh_signp(1.0)
            mref.set_magh_keqv_minus(False)
            mref.set_magh_mult(2)
            refList[i] = mref
            i += 1            
        if sp2 <= sMax:
            mref = MagReflection(reflection)
            mref.set_magh_h(FloatVector(hkp2))
            mref.set_magh_s(sp2)
            mref.set_magh_num_k(2)
            mref.set_magh_signp(-1.0)
            mref.set_magh_keqv_minus(False)
            mref.set_magh_mult(1)
            refList[i] = mref
            i += 1
        if sm2 <= sMax:
            mref = MagReflection(reflection)
            mref.set_magh_h(FloatVector(hkm2))
            mref.set_magh_s(sm2)
            mref.set_magh_num_k(2)
            mref.set_magh_signp(1.0)
            mref.set_magh_keqv_minus(False)
            mref.set_magh_mult(1)
            refList[i] = mref
            i += 1        
                #if hkl[2] == 0 and hkl[0] <= 1:
                    ## this part causes double counting
                    #ks = np.add([-hkl[0],-hkl[1],0.0],kvec)
                    #mref = MagReflection()
                    #mref.set_magh_h(FloatVector(ks))
                    #mref.set_magh_s(calcS(cell, ks))
                    #mref.set_magh_num_k(1)
                    #mref.set_magh_signp(-1.0)
                    #mref.set_magh_keqv_minus(True)
                    #mref.set_magh_mult(2)
                    #refList[i] = mref
                    #i += 1                    
            #if sm <= sMax:
                ## minus k reflections
                #mref = MagReflection()
                #mref.set_magh_h(FloatVector(hkm))
                #mref.set_magh_s(sm)
                #mref.set_magh_num_k(1)
                #mref.set_magh_signp(1.0)
                #mref.set_magh_keqv_minus(True)
                #mref.set_magh_mult(2)
                #refList[i] = mref
                #i += 1
        #elif hkl[1] in [2*n for n in range(-3,3)] and hkl[2] in [2*n for n in range(-3,3)]:
            ## special condition
            #if sp <= sMax:
                #mref = MagReflection()
                #mref.set_magh_h(FloatVector(hkp))
                #mref.set_magh_s(sp)
                #mref.set_magh_num_k(1)
                #mref.set_magh_signp(-1.0)
                #mref.set_magh_keqv_minus(True)
                #mref.set_magh_mult(2)
                #refList[i] = mref
                #i += 1
                #print hkl
    refList.set_magh_list_nref(i)
    return refList
def AbsorptionCorrection(tt, muR):
    Sth2 = np.sin((tt/2.0)*np.pi/180.0)**2
    T0 = 16.0/(3.*np.pi)
    T1 = (25.99978-0.01911*Sth2**0.25)*np.exp(-0.024551*Sth2)+ \
            0.109561*np.sqrt(Sth2)-26.04556
    T2 = -0.02489-0.39499*Sth2+1.219077*Sth2**1.5- \
            1.31268*Sth2**2+0.871081*Sth2**2.5-0.2327*Sth2**3
    T3 = 0.003045+0.018167*Sth2-0.03305*Sth2**2
    Trns = -T0*muR-T1*muR**2-T2*muR**3-T3*muR**4
    return np.exp(Trns)
# calcStructFact: calculates the structure factor squared for a list of planes
#   using provided atomic positions
def calcStructFact(refList, atomList, spaceGroup, wavelength, xtal=False):
    wavelength_p = floatp()
    wavelength_p.assign(wavelength)
    funcs.init_calc_strfactors(refList, atomList, spaceGroup, 'NUC', wavelength_p)
    structFacts = [float() for i in xrange(refList.get_reflection_list_nref())]
    reflections = refList[:]
    if xtal: code = 'S' 
    else: code = 'P'
    for i, reflection in enumerate(reflections):
        # calculates the square of the structure factor
        sFpoint = floatp()
        sFpoint.assign(structFacts[i])
        funcs.calc_strfactor(code, 'NUC', i+1, float(reflection.get_reflection_s()**2), atomList, spaceGroup, sFpoint)
        structFacts[i] = sFpoint.value()
    return structFacts

## calcMagStructFact: calculates the magnetic structure factors around a list
##   of lattice reflections
#def calcMagStructFact(refList, atomList, symmetry, cell): 
    ##funcs.print_strfac_stuff(cell, atomList, symmetry, refList)
    #funcs.init_mag_structure_factors(refList, atomList, symmetry)
    #funcs.mag_structure_factors(cell, atomList, symmetry, refList)
    ## calculate the "magnetic interaction vector" (the square of which is
    ##   proportional to the intensity)    
    #funcs.calc_mag_interaction_vector(refList, cell)
    ##for atom in atomList:
        ##print atom.basis()
    #mivs = np.array([ref.magIntVec() for ref in refList])
    #return mivs
# calcMagStructFact: calculates the magnetic structure factors around a list
#   of lattice reflections
def calcMagStructFact(refList, atomList, symmetry, cell): 
    # calculate the "magnetic interaction vectors" (the square of which is
    #   proportional to the intensity)    
    sqMivs = []
    for ref in refList:
        funcs.calc_magnetic_strf_miv(cell, symmetry, atomList, ref)
        sqMivs.append(ref.sqMiV)
    return np.array(sqMivs)
# calcIntensity: calculates the intensity for a given set of reflections,
#   based on the structure factor
def calcIntensity(refList, atomList, spaceGroup, wavelength, cell=None,
                  magnetic=False, xtal=False, extinctions=None, scale=None, muR=None):
    # TODO: make sure magnetic phase factor is properly being taken into account
    if (refList.magnetic):
        sfs2 = calcMagStructFact(refList, atomList, spaceGroup, cell)#np.array([np.sum(np.array(sf)*np.conj(np.array(sf))) for sf in sfs])
        multips = np.array([ref.get_magh_mult() for ref in refList])
        tt = np.radians(np.array([twoTheta(ref.get_magh_s(), wavelength) for ref in refList]))
        #sfs2 *= multips
    else:
        sfs2 = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength, xtal=xtal))
        multips = np.array([ref.get_reflection_mult() for ref in refList])
        tt = np.radians(np.array([twoTheta(ref.get_reflection_s(), wavelength) for ref in refList]))
        sfs2 *= multips
#    lorentz = (1+np.cos(tt)**2) / (np.sin(tt)*np.sin(tt/2))
    lorentz = (np.sin(tt)*np.sin(tt/2)) ** -1
    if muR != None: return sfs2 *lorentz * AbsorptionCorrection(tt, muR)
    return sfs2 * lorentz

# makePeaks() creates a series of Peaks to represent the powder
#   diffraction pattern
def makePeaks(reflections, coeffs, I, scale, wavelength, shape="Gaussian", eta=None, base=0):
    Peak.scaleFactor = scale
    peaks = [Peak(twoTheta(rk.s, wavelength),
                          coeffs[0], coeffs[1], coeffs[2], Ik, rk.hkl, shape, eta)
                 for rk,Ik in zip(reflections,I)]
    return peaks

# getIntensity: calculates the intensity at a given 2*theta position, or for an
#   array of 2*theta positions
def getIntensity(peaks, background, tt, base=0):
    #return background(tt) + sum(g(tt) for g in peaks)
    v = background(tt)
    for g in peaks:
        g.add(v,tt)
    return np.array(v)+base


# removeRange: takes in an array of 2*theta intervals and removes them from
#   consideration for data analysis, with an optional argument for removing the
#   corresponding intensities as well
def removeRange(tt, remove, intensity=None):
    if (remove == None or len(remove) < 1):
        if (intensity != None): return (tt, intensity)
        else: return tt
    if (not isSequence(remove[0]) or len(remove[0]) == 1):
        # single interval
        keepEntries = (tt < remove[0]) | (tt > remove[1])
        tt = tt[keepEntries]
        if (intensity != None):
            intensity = np.array(intensity)[keepEntries]
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
                observedData=(None,None), labels=None, base=0, residuals=False, error=None, muR=None):
    background = LinSpline(backgroundFile)
    sMin, sMax = getS(ttMin, wavelength), getS(180.0, wavelength)
    if magnetic:
        if (infoFile != None):
            infofile = readMagInfo(infoFile)
            if (spaceGroup == None): spaceGroup = infofile[0]
            if (cell == None): cell = infofile[1]
            if (magAtomList == None): magAtomList = infofile[2]
            if (symmetry == None): symmetry = infofile[3]
        if (basisSymmetry == None): basisSymmetry = symmetry
        ## magnetic peaks
        # convert magnetic symmetry to space group
        latt = getMagsymmK_latt(basisSymmetry)
        if basisSymmetry.get_magsymm_k_mcentred() == 1: 
            latt+= " -1" 
        else:
            latt += " 1"
        spg = SpaceGroup()
        funcs.set_spacegroup(latt, spg)
        #funcs.write_spacegroup(spg)
        # use this space group to generate magnetic hkls (refList2)
        #for i in range(len(magAtomList)):
            #matom = magAtomList[i]
            #matom.set_matom_mmphas(FloatVector([1.0 for j in range(12)]))
            #magAtomList[i] = matom
        refList = hklGen(spaceGroup, cell, sMin, sMax, True, xtal=False)
        refList2 = hklGen(spg, cell, sMin, np.sin(179.5/2)/wavelength, True, xtal=True)
        magRefList = satelliteGen(cell, symmetry, sMax, hkls=refList2)#satelliteGen_python(cell, sMax, None)#
        #newList = []
        #for ref in magRefList:
            #if ref not in newList:
                #newList.append(ref)
        #magRefList = ReflectionList(newList)
        print "length of reflection list " + str(len(magRefList))
        magIntensities = calcIntensity(magRefList, magAtomList, basisSymmetry,
                                       wavelength, cell, True, muR=muR)
        # add in structural peaks
        if (atomList == None): atomList = readInfo(infoFile)[2]
        #refList = hklGen(spaceGroup, cell, sMin, sMax, True, xtal=xtal)
        intensities = calcIntensity(refList, atomList, spaceGroup, wavelength, muR=muR)
        reflections = magRefList[:] + refList[:]
        intensities = np.append(magIntensities, intensities)
    else:
        if (infoFile != None):
            infofile = readInfo(infoFile)
            if (spaceGroup == None): spaceGroup = infofile[0]
            if (cell == None): cell = infofile[1]
            if (atomList == None): atomList = infofile[2]         
        refList = hklGen(spaceGroup, cell, sMin, sMax, True, xtal=False)
        reflections = refList[:]
        intensities = calcIntensity(refList, atomList, spaceGroup, wavelength, muR=muR)
    peaks = makePeaks(reflections, uvw, intensities, scale, wavelength, base=base)
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    tt = np.linspace(ttMin, ttMax, numPoints)
    intensity = getIntensity(peaks, background, tt, base=base)

    if info:
        if magnetic:
            printInfo(cell, spaceGroup, (atomList, magAtomList), (refList, magRefList),
                      wavelength, basisSymmetry, muR=muR)
        else:
            printInfo(cell, spaceGroup, atomList, refList, wavelength, muR=muR)
    if plot:
        plotPattern(peaks, background, observedData[0], observedData[1],
                    ttMin, ttMax, ttStep, exclusions, labels=labels, base=base, residuals=residuals, error=error)
        pylab.show()
    if saveFile:
        np.savetxt(saveFile, (tt, intensity), delimiter=" ")
    return (tt, intensity)


# printInfo: prints out information about the provided space group and atoms,
#   as well as the generated reflections
def printInfo(cell, spaceGroup, atomLists, refLists, wavelength, symmetry=None, muR=None):
    print "Wavelength:", wavelength
    if (isinstance(refLists, ReflectionList)):
        atomLists = (atomLists,)
        refLists = (refLists,)
    
    divider = "-" * 40
    print "Cell information (%s cell)" % rstrip(getSpaceGroup_crystalsys(spaceGroup))
    print divider
    print " a = %.3f   alpha = %.3f" % (cell.length()[0], cell.angle()[0])
    print " b = %.3f   beta  = %.3f" % (cell.length()[1], cell.angle()[1])
    print " c = %.3f   gamma = %.3f" % (cell.length()[2], cell.angle()[2])
    print divider
    print
    print "Space group information"
    print divider
    print "               Number: ", spaceGroup.get_space_group_numspg()
    print "           H-M Symbol: ", getSpaceGroup_spg_symb(spaceGroup)
    print "          Hall Symbol: ", getSpaceGroup_hall(spaceGroup)
    print "       Crystal System: ", getSpaceGroup_crystalsys(spaceGroup)
    print "           Laue Class: ", getSpaceGroup_laue(spaceGroup)
    print "          Point Group: ", getSpaceGroup_pg(spaceGroup)
    print " General Multiplicity: ", spaceGroup.get_space_group_multip()
    print divider
    print
    print "Atom information (%d atoms)" % len(atomLists[0])
    print divider
    atomList = atomLists[0]
    magnetic = atomList.magnetic
    label = [rstrip(getAtom_lab(atom)) for atom in atomList]
    x, y, z = tuple(["%.3f" % atom.coords()[i] for atom in atomList]
                    for i in xrange(3))
    multip = [str(atom.get_atom_mult()) for atom in atomList]
    occupancy = ["%.3f" % (atom.get_atom_occ()*spaceGroup.get_space_group_multip()/atom.get_atom_mult())
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
        intensity = ["%.3f" % I for I in calcIntensity(refList, atomList, symmObject, wavelength, cell, magnetic, muR=muR)]
        #dtype = [('tt', float),('h', 'S10'), ('k', 'S10'), ('l','S10'), ('intensity', 'S10')]
        #array1 = np.array([(tt[i], str(float(h[i])+0.5),k[i],str(float(l[i])+0.5),intensity[i]) for i in range(len(tt))], dtype=dtype)
        #array2 = np.sort(array1, order='tt')
        #print array2
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

# plotPattern: given a series of Peaks and a background, plots the predicted
#   intensity at every 2*theta position in a specified range, as well as the
#   observed intensity everywhere on a given list of points
#   used for powder patterns
def plotPattern(peaks, background, ttObs, observed, ttMin, ttMax, ttStep,
                exclusions=None, labels=None, residuals=False, base=0, error=None):
    # TODO: scale residual plot
    numPoints = int(floor((ttMax-ttMin)/ttStep)) + 1
    ttCalc = np.linspace(ttMin, ttMax, numPoints)
    if(exclusions != None): ttCalc = removeRange(ttCalc, exclusions)
    intensity = np.array(getIntensity(peaks, background, ttCalc, base=base))
    pylab.subplot(211)
    if (observed != None):
        if exclusions:
            ttObs, observed = removeRange(ttObs, exclusions, observed)
        pylab.plot(ttObs, observed, '-go', linestyle="None", label="Observed",lw=1)
    pylab.plot(ttCalc, np.array(intensity), '-b', label="Calculated", lw=1)
    intensityCalc = np.array(getIntensity(peaks, background, ttObs, base=base))
    pylab.errorbar(ttObs, np.array(observed), yerr=error, fmt=None, ecolor='g')
#    pylab.fill_between(ttObs, observed, intensity, color="lightblue")
    pylab.xlabel(r"$2 \theta$")
    pylab.ylabel("Intensity")
    pylab.legend()
    if labels:
        for g in peaks:
            if (g.center <= ttMax):
                pylab.text(g.center, np.interp(g.center, ttCalc, intensity),
                           hklString(g.hkl),
                           ha="center", va="bottom", rotation="vertical")
    if (residuals):
        resid = observed - intensityCalc
        pylab.subplot(212)
        pylab.plot(ttObs, resid, label="Residuals")
        pylab.yticks(range(0,int(max(intensity)), int(max(intensity))/10))
    return
if __name__ == '__main__':
    DATAPATH = os.path.dirname(os.path.abspath(__file__))
    dataFile = os.path.join(DATAPATH,r"HOBK/hobk.dat")
    bkgFile = os.path.join(DATAPATH,r"HOBK/hobk_bas.bac")
    tt, observed = readIllData(dataFile, "D1B", bkgFile)
    print tt
    print observed
    #pylab.plot(tt,observed,marker='s',linestyle='None')
    #pylab.show()