from pycrysfml import *
from fswig_hklgen import *
from hkl_model import TriclinicCell, MonoclinicCell, OrthorhombicCell, TetragonalCell, HexagonalCell, CubicCell, makeCell, AtomListModel, AtomModel
from string import rstrip, ljust, rjust, center
import sys
try:
    from bumps.names import Parameter, FitProblem
except(ImportError):
    pass
# Class Objects
class sXtalPeak(object):
    def __init__(self, sfs2, svalue, hkl, error=None):
        self.sfs2 = sfs2
        self.svalue = svalue
        self.hkl = hkl
        self.error = error
    def __eq__(self, other):
        #if self.svalue == other.svalue: print self.svalue, " != ", other.svalue
        #return approxEq(self.svalue, other.svalue, 0.00001)
        return self.hkl == other.hkl
    def __str__(self):
        return "Single Crystal Peak at: "+ str(self.svalue) + " with intensity: " + str(self.sfs2)

# functions
def readIntFile(filename, skiplines=3, exclusions=None, kind="dat", cell=None):
    # TODO: implement exclusions
    if kind == "dat":
        HKLs = np.loadtxt(filename, dtype=int, usecols=(0,1,2), skiprows=skiplines, comments='!')
        data = np.loadtxt(filename, dtype=float, usecols=(3,4,6,7,8,9), skiprows=skiplines, comments='!')
        wavelength = np.loadtxt(filename, dtype=float, skiprows=skiplines-1, usecols=[0], comments='!')[0]
        refList = ReflectionList()
        refList.set_reflection_list_nref(len(HKLs[:,0]))
        funcs.alloc_refllist_array(refList)
        tt = data[:,2]
        for i in range(len(HKLs[:,0])):
            reflection = Reflection()
            hkl = IntVector(HKLs[i,:])
            reflection.set_reflection_h(hkl)
            reflection.set_reflection_s(getS(tt[i], wavelength))
            reflection.set_reflection_mult(1)
            refList[i] = reflection
        # return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
        return wavelength, refList, data[:,0], data[:,1], tt, data[:,3:]    
    else:
        HKLs = np.loadtxt(filename, dtype=int, usecols=(0,1,2), skiprows=skiplines, comments='!')
        data = np.loadtxt(filename, dtype=float, usecols=(3,4), skiprows=skiplines, comments='!')
        wavelength = np.loadtxt(filename, dtype=float, skiprows=skiplines-1, usecols=[0], comments='!')[0]
        refList = ReflectionList()
        refList.set_reflection_list_nref(len(HKLs[:,0]))
        funcs.alloc_refllist_array(refList)
        for i in range(len(HKLs[:,0])):
            reflection = Reflection()
            hkl = IntVector(HKLs[i,:])
            reflection.set_reflection_h(hkl)
            reflection.set_reflection_s(calcS(cell, hkl))
            reflection.set_reflection_mult(1)
            refList[i] = reflection
        # return wavelength, refList, sfs2, error
        return wavelength, refList, data[:,0], data[:,1]      
def readMagIntFile(filename, cell=None):
    f = [line.strip().split() for line in open(filename)]
    hkls = []
    sfs2 = []
    error = []
    kvec = []
    reflist = []
    for line in f:
        if len(line) == 7:
            hkls.append([float(line[0]), float(line[1]), float(line[2])])
            sfs2.append(float(line[4]))
            error.append(float(line[5]))
        elif len(line) == 4:
            kvec = np.array([float(line[1]), float(line[2]), float(line[3])])
    for hkl in hkls:
        # add check for k == -k
        hkl = hkl+kvec
        reflection = MagReflection()
        reflection.set_magh_h(FloatVector(hkl))
        reflection.set_magh_s(calcS(cell, FloatVector(hkl)))
        reflection.set_magh_mult(1)
        reflist.append(reflection)
    return ReflectionList(reflist), sfs2, error
# Create a list of single crystal peak objects from a list 
# of structure factors squared and sin(theta)/lambda values
def makeXtalPeaks(sfs2, svalues, refList, peaks=None, error=None):
    if peaks == None:
        peaks = []
    for i in range(len(svalues)):
        if error != None:
            p = sXtalPeak(sfs2[i], svalues[i], refList[i].hkl, error[i])
        else:
            p = sXtalPeak(sfs2[i], svalues[i], refList[i].hkl)
        if p not in peaks:
            peaks.append(p)
        else:
            #print "Peak at: ", svalues[i], "adding: ", peaks[peaks.index(p)].sfs2, " to: ", sfs2[i], " = ", peaks[peaks.index(p)].sfs2+sfs2[i]
            peaks[peaks.index(p)].sfs2 += sfs2[i]
            if peaks[peaks.index(p)].error != None:
                peaks[peaks.index(p)].error = np.sqrt(peaks[peaks.index(p)].error**2+error[i]**2)
            #peaks.append(p)
            pass
    return peaks
# Check if an intensity's sin(theta)/lambda value is approximately in a list of st/l values
def checkInt(value, sCalc):
    index = 0
    for s in sCalc:
        if approxEq(value, s, 0.001):
            return True, index
        index += 1
    return False, False
def getXtalIntensity(peaks, sList=None, background=None, exclusions=None, base=0, scale=1):
    if background == None:
        background = np.zeros(len(sList))
    if sList == None:
        return (np.array([peak.sfs2 for peak in peaks])*scale)
    else:
        intensities = []
        icalc = [peak.sfs2 for peak in peaks]
        scalc = [peak.svalue for peak in peaks]
        #for s in sList:
            #if checkInt(s, scalc)[0]:
                #intensities.append(icalc[checkInt(s, scalc)[1]])
            #else:
                #intensities.append(-10000.0)
        #for peak in peaks:
            #if peak.svalue in sList:
                #intensities.append(peak.sfs2)
        return np.array(icalc)*scale, np.array(scalc)
# calcIntensity: calculates the intensity for a given set of reflections,
#   based on the structure factor
def calcXtalIntensity(refList, atomList, spaceGroup, wavelength, cell=None,
                  magnetic=False, extinctions=None, scale=None):
    # TODO: make sure magnetic phase factor is properly being taken into account
    if (refList.magnetic):
        sfs2 = calcMagStructFact(refList, atomList, spaceGroup, cell)
        tt = np.radians(np.array([twoTheta(ref.get_magh_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_magh_s() for ref in refList])
    else:
        sfs2 = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength, xtal=True))
        tt = np.radians(np.array([twoTheta(ref.get_reflection_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_reflection_s() for ref in refList])
    if extinctions != None:
        extfs = []
        for i in range(len(sfs2)):
            #newsfs.append(extinctionFactor(sfs2[i], wavelength, tt[i], scale, extinctions))
            while len(extinctions) < 6:
                extinctions.append(0)
            extf = floatp()
            funcs.shelx_extinction(3, 1, wavelength, svalues[i]**2, FloatVector(refList[i].hkl), sfs2[i], FloatVector(extinctions), extf)
            extfs.append(extf.value())
        sfs2 *= np.array(extfs)
    return sfs2, svalues

# diffPatternXtal: generates a neutron diffraction pattern from a file containing
#   crystallographic information or from the same information generated
#   elsewhere Use this version for single crystal data
def diffPatternXtal(infoFile=None, backgroundFile=None, wavelength=1.5403,
                    tt=None, exclusions=None,
                    spaceGroup=None, cell=None, atomList=None,
                    symmetry=None, basisSymmetry=None, magAtomList=None, scale=1,
                    magnetic=False, nuclear=True, info=False, plot=False, saveFile=None,
                    obsIntensity=None, labels=None, base=0, residuals=False, error=None, refList=None, extinctions=None):
    background = None
    sMin, sMax = getS(min(tt), wavelength), getS(max(tt), wavelength)
    reflections = None
    peaks = None
    if magnetic:
        if (infoFile != None):
            infofile = readMagInfo(infoFile)
            if (spaceGroup == None): spaceGroup = infofile[0]
            if (cell == None): cell = infofile[1]
            if (magAtomList == None): magAtomList = infofile[2]
            if (symmetry == None): symmetry = infofile[3]
        if (basisSymmetry == None): basisSymmetry = symmetry
        magRefList = satelliteGen(cell, symmetry, np.sin(179.5/2)/wavelength, hkls=refList)
        print "length of reflection list " + str(len(magRefList))
        sfs2, svalues = calcXtalIntensity(magRefList, magAtomList, basisSymmetry,
                                      wavelength, cell, True, extinctions=extinctions, scale=scale)
        magpeaks = makeXtalPeaks(sfs2, svalues, magRefList)
        # add in structural peaks
        if (atomList == None): atomList = readInfo(infoFile)[2]
        if nuclear:
            sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, extinctions=extinctions, scale=scale)     
            peaks = makeXtalPeaks(sfs2, svalues, refList, peaks=magpeaks)
            reflections = magRefList[:] + refList[:]
            intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], background=background, exclusions=exclusions, base=base, scale=scale)
        else:
            peaks = magpeaks
            reflections = magRefList[:]
            intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], exclusions=exclusions, scale=scale)
    else:
        if (infoFile != None):
            infofile = readInfo(infoFile)
            if (spaceGroup == None): spaceGroup = infofile[0]
            if (cell == None): cell = infofile[1]
            if (atomList == None): atomList = infofile[2]         
        reflections = refList[:]
        sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, extinctions=extinctions, scale=scale)
        peaks = makeXtalPeaks(sfs2, svalues, refList, error=error)
        intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], background=background, exclusions=exclusions, base=base, scale=scale)
    if info:
        if magnetic:
            printInfo(cell, spaceGroup, (atomList, magAtomList), (refList, magRefList),
                      wavelength, basisSymmetry)
        else:
            printInfo(cell, spaceGroup, atomList, refList, wavelength)
    if plot:
        sObs = np.array([getS(value, wavelength) for value in tt])
        plotXtalPattern(peaks, sObs, obsIntensity, background=background, error=error, base=base, residuals=residuals,labels=labels, scale=scale, refList=refList)
        pylab.show()
    if saveFile:
        np.savetxt(saveFile, (tt, intensity), delimiter=" ")
    return
def plotXtalPattern(peaks, sList, obsIntensity, background=None, 
                    exclusions=None, labels=None, residuals=False, base=0, scale=1, error=None, refList=None):
    # plot single crystal intesities vs q as single points
    obspeaks = makeXtalPeaks(obsIntensity, sList, refList, error=error)
    sList = np.array([peak.svalue for peak in obspeaks])
    obsIntensity = np.array([peak.sfs2 for peak in obspeaks])
    calcIntensity, scalc = getXtalIntensity(peaks, sList=sList, exclusions=exclusions, scale=scale)
    pylab.subplot(211)
    if (obsIntensity != None):
        pylab.plot(sList*(4*np.pi), obsIntensity, '-go', linestyle="None", label="Observed",lw=1)
    pylab.plot(scalc*(4*np.pi), calcIntensity, '-bo', linestyle="None", label="Calculated", lw=1)
    error = [peak.error for peak in obspeaks]
    pylab.errorbar(sList*(4*np.pi), obsIntensity, yerr=error, fmt=None, ecolor='g')
    pylab.xlabel("Q")
    pylab.ylabel("Intensity")
    pylab.legend()
    if (residuals):
        calcPeaks = makeXtalPeaks(calcIntensity, scalc, peaks)
        resid = []
        resid_stl = []
        for peak in obspeaks:
            resid_stl.append(peak.svalue)
            if peak in calcPeaks:
                resid.append(peak.sfs2-calcPeaks[calcPeaks.index(peak)].sfs2)
            else:
                resid.append(peak.sfs2)
        #resid = obsIntensity - calcIntensity
        pylab.subplot(212)
        pylab.plot(np.array(resid_stl)*(4*np.pi), np.array(resid), '-bo', linestyle="None", label="Residuals")
    return

# printInfo: prints out information about the provided space group and atoms,
#   as well as the generated reflections
def printXtalInfo(cell, spaceGroup, atomLists, refLists, wavelength, symmetry=None):
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
        intensity = ["%.3f" % I for I in calcXtalIntensity(refList, atomList, symmObject, wavelength, cell, magnetic)[0]]
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
    
# bumps model
# Model: represents an object that can be used with bumps for optimization
#   purposes single crystal version
class Model(object):

    def __init__(self, tt, observed, background,
                 wavelength, spaceGroupName, cell, atoms, exclusions=None,
                 magnetic=False, symmetry=None, newSymmetry=None, scale=1, zero=None, error=None, hkls=None, extinction=[0]):
        if (isinstance(spaceGroupName, SpaceGroup)):
            self.spaceGroup = spaceGroupName
        else:
            self.spaceGroup = SpaceGroup(spaceGroupName)
        self.tt = np.array(tt)
        self.obspeaks = makeXtalPeaks(observed, [getS(ttval, wavelength) for ttval in self.tt], refList=hkls, error=error)
        self.sList = np.array([peak.svalue for peak in self.obspeaks])
        self.observed = np.array([peak.sfs2 for peak in self.obspeaks])
        self.background = background
        self.scale = Parameter(scale, name='scale')
        self.extinctions = [Parameter(extinction[i], name="Extinction"+str(i)) for i in range(len(extinction))]
        self.error = error
        self.refList = hkls
        self.base=0
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
        if self.magnetic:
            self.magRefList = satelliteGen(self.cell.cell, self.symmetry, np.sin(179.5/2)/self.wavelength, hkls=self.refList)
            self.magReflections = self.magRefList[:]
        self.reflections = self.refList
            
    def __getstate__(self):
        state = self.__dict__.copy()
        del state["refList"]
        del state["magRefList"]
        return state
    
    def __setstate__(self, state):
        self.__dict__ = state
        self._set_reflections()

    def parameters(self):
        if self.has_zero:
            params = {
                    'scale': self.scale,
                    'zero' : self.zero,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
            for p in self.extinctions:
                params[p.name] = p
            return params            
        else:
            params = { 
                    'scale': self.scale,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
            for p in self.extinctions:
                params[p.name] = p
            return params
        
    def numpoints(self):
        return len(self.observed)

    def theory(self):
        return getXtalIntensity(self.peaks, background=self.background, sList=self.sList, scale=self.scale.value)[0]

    def residuals(self):
        calcPeaks = makeXtalPeaks(self.theory(), self.sList, self.peaks)
        resid = []
        resid_stl = []
        for peak in self.obspeaks:
            resid_stl.append(peak.svalue)
            if peak in calcPeaks:
                resid.append(peak.sfs2-calcPeaks[calcPeaks.index(peak)].sfs2)
            else:
                resid.append(peak.sfs2)
        return np.array(resid)/(np.sqrt(self.observed)+1)

    def nllf(self):
        return np.sum(self.residuals()**2)

    def plot(self, view="linear"):
        import pylab
        if self.has_zero:
            zero = self.zero.value  
        else:
            zero = self.zero
        plotXtalPattern(self.peaks, self.sList, self.observed, 
                       background=self.background, 
                       exclusions=self.exclusions, 
                       residuals=True, refList=self.reflections,
                       error=self.error, scale=self.scale.value)
    def update(self):  
        self.cell.update()
        self.atomListModel.update()
        if self.magnetic:
            # update magnetic reflections and add their peaks to the list of
            #   Peaks         
            hkls = [reflection.hkl for reflection in self.magReflections]
            sList = calcS(self.cell.cell, hkls)
            for i in xrange(len(self.magReflections)):
                self.magReflections[i].set_magh_s(sList[i])
            sfs2, svalues = calcXtalIntensity(self.magRefList, self.atomListModel.magAtomList, self.symmetry, self.wavelength, magnetic=True, cell=self.cell.cell, extinctions=[ext.value for ext in self.extinctions], scale=self.scale.value)
            self.magIntensities = sfs2
            #print self.magIntensities
            self.peaks = makeXtalPeaks(sfs2, svalues, self.magRefList)
            #self.sList = np.array([peak.svalue for peak in self.peaks])
            #self.peaks.extend(makeXtalPeaks(sfs2, svalues))        
        hkls = [reflection.hkl for reflection in self.reflections]
        sList = calcS(self.cell.cell, hkls)
        for i in xrange(len(self.reflections)):
            self.reflections[i].set_reflection_s(sList[i])
        sfs2, svalues = calcXtalIntensity(self.refList, self.atomListModel.atomList, self.spaceGroup, self.wavelength, extinctions=[ext.value for ext in self.extinctions], scale=self.scale.value)
        self.intensities = sfs2
        if not self.magnetic: self.peaks = None
        self.peaks = makeXtalPeaks(sfs2, svalues, self.refList, peaks=self.peaks)
        #self.sList = svalues