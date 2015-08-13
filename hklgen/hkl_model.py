from pycrysfml import *
from fswig_hklgen import *
from string import rstrip, ljust, rjust, center
import sys
try:
    from bumps.names import Parameter, FitProblem
except(ImportError):
    pass
# TODO: Fix imports

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
    (a, b, c) = cell.length()
    (alpha, beta, gamma) = cell.angle()
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
                 magnetic=False, symmetry=None, newSymmetry=None, base=None, scale=1, eta=0,zero=None, error=None, hkls=None, muR=None):
        if (isinstance(spaceGroupName, SpaceGroup)):
            self.spaceGroup = spaceGroupName
        else:
            self.spaceGroup = SpaceGroup(spaceGroupName)
        self.tt = np.array(tt)
        self.observed = observed
        self.background = background
        self.u = Parameter(u, name='u')
        self.v = Parameter(v, name='v')
        self.w = Parameter(w, name='w')
        self.scale = Parameter(scale, name='scale')
        self.eta = Parameter(eta, name='eta')
        self.error = error
        self.muR = muR
        # hkls for single crystal only
        if hkls != None:
            self.refList = hkls
        if base != None:
            self.base = Parameter(base, name='base')
            self.has_base = True
        else:
            self.base=0
            self.has_base = False
        if zero != None:
            self.zero = Parameter(zero, name='zero')
            self.has_zero = True
        else:
            self.zero = 0
            self.has_zero = False
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
            self.atomListModel = AtomListModel(atoms, self.spaceGroup.get_space_group_multip(),
                                               True, self.newSymmetry)            
        else:
            self.atomListModel = AtomListModel(atoms, self.spaceGroup.get_space_group_multip(), False)
        self._set_reflections()           
        self.update()
    def _set_reflections(self):
        maxLattice = self.cell.getMaxLattice()
        maxCell = CrystalCell(maxLattice[:3], maxLattice[3:])
        # powder calculations
        self.refList = hklGen(self.spaceGroup, self.cell.cell, self.sMin, self.sMax, True, xtal=False)
        self.reflections = self.refList[:]
        if self.magnetic:
            # convert magnetic symmetry to space group
            latt = getMagsymmK_latt(self.symmetry)
            if self.symmetry.get_magsymm_k_mcentred() == 1: 
                latt+= " -1" 
            else:
                latt += " 1"
            spg = SpaceGroup()
            funcs.set_spacegroup(latt, spg)
            # use this space group to generate magnetic hkls      
            hkls = hklGen(spg, self.cell.cell, self.sMin, np.sin(179.5/2)/self.wavelength, True, xtal=True)
            self.magRefList = satelliteGen(self.cell.cell, self.symmetry, self.sMax, hkls=hkls)#satelliteGen_python(self.cell.cell, self.sMax, hkls)
            newList = []
            for ref in self.magRefList:
                if ref not in newList:
                    newList.append(ref)
            self.magRefList = ReflectionList(newList)                
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
        if self.has_base and self.has_zero:
            return {'u': self.u,
                    'v': self.v,
                    'w': self.w,
                    'scale': self.scale,
                    'eta': self.eta,
                    'base': self.base,
                    'zero' : self.zero,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
        elif self.has_base:
            return {'u': self.u,
                    'v': self.v,
                    'w': self.w,
                    'scale': self.scale,
                    'eta': self.eta,
                    'base': self.base,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
        elif self.has_zero:
            return {'u': self.u,
                    'v': self.v,
                    'w': self.w,
                    'scale': self.scale,
                    'eta': self.eta,
                    'zero' : self.zero,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
        else:
            return {'u': self.u,
                    'v': self.v,
                    'w': self.w,
                    'scale': self.scale,
                    'eta': self.eta,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
        
    def numpoints(self):
        return len(self.observed)

    def theory(self):
        # add check in case zero not defined
        if self.has_base and self.has_zero:
            return getIntensity(self.peaks, self.background, removeRange(self.tt, self.exclusions)-self.zero.value, base=self.base.value)
        elif self.has_base:
            return getIntensity(self.peaks, self.background, removeRange(self.tt, self.exclusions)-self.zero, base=self.base.value)
        elif self.has_zero:
            return getIntensity(self.peaks, self.background, removeRannge(self.tt, self.exclusions)-self.zero.value, base=self.base)
        else:
            return getIntensity(self.peaks, self.background, removeRange(self.tt, self.exclusions)-self.zero, base=self.base)

    def residuals(self):
        return (self.theory() - self.observed)/(np.sqrt(self.observed)+1)

    def nllf(self):
        return np.sum(self.residuals()**2)

    def plot(self, view="linear"):
        import pylab
        if self.has_base and self.has_zero:
            base, zero = self.base.value, self.zero.value
        elif self.has_base:
            base, zero = self.base.value, self.zero
        elif self.has_zero:
            base, zero = self.base, self.zero.value  
        else:
            base, zero = self.base, self.zero
        plotPattern(self.peaks, self.background, self.tt-zero, self.observed,
                            self.ttMin, self.ttMax, 0.01, self.exclusions, labels=None, base = base, residuals=True, error=self.error)          

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
            self.reflections[i].set_reflection_s(getS(ttPos[i], self.wavelength))
        self.intensities = calcIntensity(self.refList, self.atomListModel.atomList, self.spaceGroup, self.wavelength, muR=self.muR)
        self.peaks = makePeaks(self.reflections,
                                       [self.u.value, self.v.value, self.w.value],
                                       self.intensities, self.scale.value,
                                       self.wavelength, shape="pseudovoigt", eta=self.eta)
        if self.magnetic:
            # update magnetic reflections and add their peaks to the list of
            #   Peaks         
            hkls = [reflection.hkl for reflection in self.magReflections]
            sList = calcS(self.cell.cell, hkls)
            ttPos = np.array([twoTheta(s, self.wavelength) for s in sList])
            # move nonexistent peaks (placed at 180) out of the way to 2*theta = -20
            ttPos[np.abs(ttPos - 180*np.ones_like(ttPos)) < 0.0001] = -20
            for i in xrange(len(self.magReflections)):
                self.magReflections[i].set_magh_s(getS(ttPos[i], self.wavelength))            
            #printInfo(self.cell.cell, self.spaceGroup, [self.atomListModel.atomList, self.atomListModel.magAtomList], [self.refList,self.magRefList], self.wavelength, symmetry=self.newSymmetry)
            self.magIntensities = calcIntensity(self.magRefList,
                                                self.atomListModel.magAtomList, 
                                                self.newSymmetry, self.wavelength,
                                                self.cell.cell, True, muR=self.muR)
            #print self.magIntensities
            self.peaks.extend(makePeaks(self.magReflections,
                                           [self.u.value, self.v.value, self.w.value],
                                           self.magIntensities, self.scale.value,
                                           self.wavelength, shape="pseudovoigt", eta=self.eta))

class AtomListModel(object):
    # TODO: make occupancy constraints automatic
    
    def __init__(self, atoms, sgmultip, magnetic=False, symmetry=None):
        self.sgmultip = sgmultip
        self.magnetic = magnetic
        self.symmetry = symmetry
        self._rebuild_object(atoms)
        ## special parameters for changing basis
        #self.angle = [None] * symmetry.numSymOps
        #for i in xrange(symmetry.numSymOps):
            #self.angle[i] = Parameter(0, name="Ang"+str(i))
        #self.magnitude = Parameter(1, name="Mag")

    def _rebuild_object(self, atoms):
        if (not self.magnetic):
            # one list of atoms
            if (isinstance(atoms[0], AtomList)):
                self.atomList = atoms[0]
                self.atoms = atoms[0][:]
            else:
                self.atomList = AtomList(atoms[0])
                self.atoms = atoms
            self.atomModels = [AtomModel(atom, self.sgmultip) for atom in atoms[0]]
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
            self.atomModels = [AtomModel(atom, self.sgmultip) for atom in self.atoms]
        if self.magnetic:
            # correct atom models to include magnetic atoms
            for magAtom in self.magAtoms:
                for model in self.atomModels:
                    if (magAtom.label().rstrip() == model.atom.label().rstrip() and \
                        magAtom.sameSite(model.atom)):
                        model.addMagAtom(magAtom, self.symmetry)
                        print magAtom.label()
            index = 0
            for magAtom in self.magAtoms:
                print index, magAtom.label()
                self.magAtomList[index] = magAtom
                self.magAtoms[index] = magAtom
                index += 1              
        self.modelsDict = dict([(am.atom.label(), am) for am in self.atomModels])
#        print >>sys.stderr, "atom models: ", [(am.atom.label, am.magnetic)
#                                              for am in self.atomModels]

    def __getstate__(self):
        state = self.atoms, self.sgmultip, self.magnetic
        return state

    def __setstate__(self, state):
        self.atoms, self.sgmultip, self.magnetic = state
        self._rebuild_object(self.atoms)
        
    def parameters(self):
        params = dict(zip([atom.label() for atom in self.atoms],
                          [am.parameters() for am in self.atomModels]))
        # special parameters for changing basis
        #params.update(zip([angle.name for angle in self.angle],
                          #[angle for angle in self.angle]))
        #params.update({self.magnitude.name: self.magnitude})
        return params

    def update(self):
#        print >>sys.stderr, len(self.parameters())
#        print >>sys.stderr, "label: |" + self.atoms[0].label + "|"        
#         if len(self.parameters()) == 1:
#            print >>sys.stderr, self.parameters()
        index = 0
        I2 = 0
        for atomModel in self.atomModels:
            atomModel.update()
            self.atomList[I2] = atomModel.atom
            I2 += 1
            if hasattr(atomModel, 'magAtoms') and atomModel.magAtoms != None:
                for matom in atomModel.magAtoms:
                    self.magAtomList[index] = matom
                    self.magAtoms[index] = matom
                    #print atomModel.magAtom.basis()
                    index += 1
        # update basis vectors instead of coefficients (only needed in special
        #   circumstances)
        #if self.magnetic:
            #for i in xrange(self.symmetry.numSymOps):
                ## correct for hexagonal coordinates
                #xcomp = self.magnitude.value * np.cos(radians(self.angle[i].value))
                #ycomp = self.magnitude.value * np.sin(radians(self.angle[i].value))
                #bcomp = ycomp / np.sin(np.radians(120))
                #acomp = xcomp + bcomp / 2
                #self.symmetry.setBasis(0, i, 0, [acomp,0,0])
                #self.symmetry.setBasis(0, i, 1, [0,bcomp,0])

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
        self.magAtoms = None
        self.matom_index = -1
        self.coeffs = None
        self.phases = None
        self.sgmultip = sgmultip
        self.atom.set_atom_lab(rstrip(self.atom.label()))
        self.B = Parameter(self.atom.BIso(), name=self.atom.label() + " B")
        occ = self.atom.occupancy() / self.atom.multip() * self.sgmultip
        self.occ = Parameter(occ, name=getAtom_lab(self.atom) + " occ")
        self.x = Parameter(self.atom.coords()[0], name=self.atom.label() + " x")
        self.y = Parameter(self.atom.coords()[1], name=self.atom.label() + " y")
        self.z = Parameter(self.atom.coords()[2], name=self.atom.label() + " z")
    
    def addMagAtom(self, magAtom, symmetry):
        # add a secondary magnetic atom object to the model
        if self.magAtoms != None:
            self.magAtoms.append(magAtom)
        else:
            self.magAtoms = [magAtom]
        self.matom_index += 1
        matom = self.magAtoms[self.matom_index]
        matom.setLabel(rstrip(matom.label()))
        self.symmetry = symmetry
        self.magnetic = True
        self.numVectors = self.symmetry.numBasisFunc()[self.magAtoms[self.matom_index].irrepNum()[0]-1]
        if self.matom_index != 0:
            param_label = " kvect: " + str(self.matom_index+1)
        else:
            param_label = ""        
        if self.coeffs == None:
            self.coeffs = [None] * self.numVectors
            j = 0
        else:
            self.coeffs.extend([None] * self.numVectors)
            j = len(self.coeffs)
        for i in xrange(self.numVectors):
            self.coeffs[i+j] = Parameter(matom.basis()[0][i], name=self.atom.label() + " C" + str(i) + param_label)
        if self.phases == None:
            self.phases = [Parameter(matom.phase[0], name=self.atom.label() + " phase"+param_label)]
        else:
            self.phases.append(Parameter(matom.phase[0], name=self.atom.label() + " phase"+param_label))
        self.magAtoms[self.matom_index] = matom
    
    def parameters(self):
        params = {self.B.name: self.B, self.occ.name: self.occ,
                  self.x.name: self.x, self.y.name: self.y, self.z.name: self.z}
        if self.magnetic:
            params.update([(coeff.name, coeff) for coeff in self.coeffs])
            for phase in self.phases:
                params.update([(phase.name, phase)])
        return params

    def update(self):
        self.atom.setBIso(self.B.value)
        occ = self.occ.value * self.atom.multip() / self.sgmultip
        self.atom.setOccupancy(occ)
        self.atom.setCoords([float(self.x), float(self.y), float(self.z)])
        
        if self.magnetic:
            for i in range(len(self.magAtoms)):
                matom = self.magAtoms[i]
                matom.setBIso(self.B.value)
                matom.setCoords([float(self.x), float(self.y), float(self.z)])
                for j in xrange(self.numVectors):
                    #self.magAtom.basis[0][i] = self.coeffs[i].value
                    matom.setBasis_ind(0,j, self.coeffs[j+i].value)
                matom.set_phase(self.phases[i].value)
                self.magAtoms[i] = matom
            # following line breaks magnetism if magnetic occupancies are different from nuclear occupancies
            # It may however, be needed in order to fit the occupancy
           # self.magAtom.setOccupancy(occ)