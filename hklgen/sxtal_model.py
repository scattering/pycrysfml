from fswig_hklgen import *
# Class Objects
class sXtalPeak(object):
    def __init__(self, sfs2, svalue):
        self.sfs2 = sfs2
        self.svalue = svalue
    def __eq__(self, other):
        return self.svalue == other.svalue
# functions
def readIntFile(filename, skiplines=3, exclusions=None):
    # TODO: implement exclusions
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
# Create a list of single crystal peak objects from a list 
# of structure factors squared and sin(theta)/lambda values
def makeXtalPeaks(sfs2, svalues):
    peaks = []
    for i in range(len(svalues)):
        p = sXtalPeak(sfs2[i], svalues[i])
        if p not in peaks:
            peaks.append(p)
        else:
            peaks[peaks.index(p)].sfs2 += sfs2[i]
    return peaks
# extinctionFactor: applies extinction coeffiecients to the structure factor squared
#   currently implemented for only one extinction coefficient
def extinctionFactor(sfs2, wavelength, tt, scale, coeffs):
    return sfs2 * (scale*(1+(0.001*coeffs[0]*sfs2*wavelength**3)/np.sin(tt))**(1/4))**2
# calcIntensity: calculates the intensity for a given set of reflections,
#   based on the structure factor
def calcXtalIntensity(refList, atomList, spaceGroup, wavelength, cell=None,
                  magnetic=False, xtal=False, extinctions=None, scale=None):
    # TODO: make sure magnetic phase factor is properly being taken into account
    if (refList.magnetic):
        sfs = calcMagStructFact(refList, atomList, spaceGroup, cell)
        sfs2 = np.array([np.sum(np.array(sf)*np.conj(np.array(sf))) for sf in sfs])
        multips = np.array([ref.get_magh_mult() for ref in refList])
        tt = np.radians(np.array([twoTheta(ref.get_magh_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_magh_s() for ref in refList])
    else:
        sfs2 = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength, xtal=xtal))
        multips = np.array([ref.get_reflection_mult() for ref in refList])
        tt = np.radians(np.array([twoTheta(ref.get_reflection_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_reflection_s() for ref in refList])
        if not xtal: sfs2 *= multips
#    lorentz = (1+np.cos(tt)**2) / (np.sin(tt)*np.sin(tt/2))
    lorentz = (np.sin(tt)*np.sin(tt/2)) ** -1
    if not xtal: return sfs2 * lorentz
    if extinctions != None:
        newsfs = []
        for i in range(len(sfs2)):
            newsfs.append(extinctionFactor(sfs2[i], wavelength, tt[i], scale, extinctions))
        sfs2 = np.array(newsfs)
    return sfs2, svalues

# diffPatternXtal: generates a neutron diffraction pattern from a file containing
#   crystallographic information or from the same information generated
#   elsewhere Use this version for single crystal data
def diffPatternXtal(infoFile=None, backgroundFile=None, wavelength=1.5403,
                    tt=None, exclusions=None,
                    spaceGroup=None, cell=None, atomList=None,
                    symmetry=None, basisSymmetry=None, magAtomList=None, scale=1,
                    magnetic=False, info=False, plot=False, saveFile=None,
                    obsIntensity=None, labels=None, base=0, residuals=False, error=None, refList=None, extinctions=None):
    background = LinSpline(backgroundFile)
    sMin, sMax = getS(min(tt), wavelength), getS(max(tt), wavelength)
    reflections = None
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
                                      wavelength, cell, True, xtal=True, extinctions=extinctions, scale=scale)
        magIntensities = makeXtalPeaks(sfs2, svalues)
        # add in structural peaks
        if (atomList == None): atomList = readInfo(infoFile)[2]
        sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, xtal=True, extinctions=extinctions, scale=scale)     
        intensities = makeXtalPeaks(sfs2, svalues)
        reflections = magRefList[:] + refList[:]
        sfs2 = np.array([peak.sfs2 for peak in magIntensities])
        sfs2 = np.append(sfs2, np.array([peak.sfs2 for peak in intensities]))
        svals = np.array([peak.svalue for peak in magIntensities])
        svals = np.append(svals, np.array([peak.svalue for peak in intensities]))
        intensities = makeXtalPeaks(sfs2, svals)
    else:
        if (infoFile != None):
            infofile = readInfo(infoFile)
            if (spaceGroup == None): spaceGroup = infofile[0]
            if (cell == None): cell = infofile[1]
            if (atomList == None): atomList = infofile[2]         
        reflections = refList[:]
        sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, xtal=True, extinctions=extinctions, scale=scale)
        intensities = makeXtalPeaks(sfs2, svalues)
    if info:
        if magnetic:
            printInfo(cell, spaceGroup, (atomList, magAtomList), (refList, magRefList),
                      wavelength, basisSymmetry)
        else:
            printInfo(cell, spaceGroup, atomList, refList, wavelength)
    if plot:
        sObs = np.array([getS(value, wavelength) for value in tt])
        sCalc = np.array([peak.svalue for peak in intensities])
        Icalc = np.array([peak.sfs2 for peak in intensities])
        plotXtalPattern(sObs, sCalc, obsIntensity, Icalc*scale, background=background, error=error, base=base, residuals=residuals,labels=labels)
        pylab.show()
    if saveFile:
        np.savetxt(saveFile, (tt, intensity), delimiter=" ")
    return
def plotXtalPattern(sObs, sCalc, obsIntensity, calcIntensity, background=None, 
                    exclusions=None, labels=None, residuals=False, base=0, error=None):
    # plot single crystal intesities vs q as single points
    peaks = makeXtalPeaks(obsIntensity, sObs)
    obsIntensity = np.array([peak.sfs2 for peak in peaks])
    sObs = np.array([peak.svalue for peak in peaks])
    pylab.subplot(211)
    if (obsIntensity != None):
        pylab.plot(sObs*(4*np.pi), obsIntensity, '-go', linestyle="None", label="Observed",lw=1)
    pylab.plot(sCalc*(4*np.pi), calcIntensity, '-bo', linestyle="None", label="Calculated", lw=1)
    pylab.errorbar(sObs*(4*np.pi), obsIntensity, yerr=[0.1373,1.3232,0.2366,3.7006,0.4503,7.4395,0.4176,1.8760], fmt=None, ecolor='g')
    pylab.xlabel("Q")
    pylab.ylabel("Intensity")
    pylab.legend()
    if (residuals):
        resid = obsIntensity - calcIntensity
        pylab.subplot(212)
        pylab.plot(sObs*(4*np.pi), resid, '-bo', linestyle="None", label="Residuals")
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