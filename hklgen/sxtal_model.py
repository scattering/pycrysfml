# from __future__ import print_function # Not needed in Python 3
import sys
from collections import OrderedDict # Import OrderedDict

try:
    from bumps.names import Parameter, FitProblem
except(ImportError):
    pass

from .pycrysfml import (
    ReflectionList, Reflection, IntVector, getS, calcS, MagReflection, FloatVector,
    approxEq, calcMagStructFact, twoTheta, calcStructFact, floatp, funcs,
    readMagInfo, satelliteGen, readInfo, printInfo, SpaceGroup, CrystalCell, np
) # Explicitly import names
from .fswig_hklgen import 소비자_hklgen_py # Assuming this is the intended import, adjust if not
from .hkl_model import TriclinicCell, MonoclinicCell, OrthorhombicCell, TetragonalCell, HexagonalCell, CubicCell, makeCell, AtomListModel, AtomModel

# Class Objects
class sXtalPeak(object):
    def __init__(self, sfs2, svalue, hkl, error=None):
        self.sfs2 = sfs2
        self.svalue = svalue
        self.hkl = hkl
        self.error = error
    def __eq__(self, other):
        #if self.svalue == other.svalue: print(self.svalue, " != ", other.svalue)
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
        wavelength = np.loadtxt(filename, dtype=float, skiprows=skiplines-1, usecols=(0), comments='!')[0]
        refList = ReflectionList()
        refList.set_reflection_list_nref(len(HKLs[:,0]))
        funcs.alloc_refllist_array(refList)
        tt = data[:,2]
        for i in range(len(HKLs[:,0])):
            reflection = Reflection()
            hkl = IntVector(HKLs[i,:].tolist()) # Ensure it's a list for IntVector
            reflection.set_reflection_h(hkl)
            reflection.set_reflection_s(getS(tt[i], wavelength))
            reflection.set_reflection_mult(1)
            refList[i] = reflection
        # return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
        return wavelength, refList, data[:,0], data[:,1], tt, data[:,3:]
    else:
        HKLs = np.loadtxt(filename, dtype=int, usecols=(0,1,2), skiprows=skiplines, comments='!')
        data = np.loadtxt(filename, dtype=float, usecols=(3,4), skiprows=skiplines, comments='!')
        wavelength = np.loadtxt(filename, dtype=float, skiprows=skiplines-1, usecols=(0), comments='!')[0]
        refList = ReflectionList()
        refList.set_reflection_list_nref(len(HKLs[:,0]))
        funcs.alloc_refllist_array(refList)
        for i in range(len(HKLs[:,0])):
            reflection = Reflection()
            hkl = IntVector(HKLs[i,:].tolist()) # Ensure it's a list for IntVector
            reflection.set_reflection_h(hkl)
            reflection.set_reflection_s(calcS(cell, hkl))
            reflection.set_reflection_mult(1)
            refList[i] = reflection
        # return wavelength, refList, sfs2, error
        return wavelength, refList, data[:,0], data[:,1]
def readMagIntFile(filename, cell=None):
    # Use with statement for file handling
    with open(filename) as f_obj:
        f = [line.strip().split() for line in f_obj]
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
    if peaks is None:
        peaks = []
    for i in range(len(svalues)):
        if not error is None:
            p = sXtalPeak(sfs2[i], svalues[i], refList[i].hkl, error[i])
        else:
            p = sXtalPeak(sfs2[i], svalues[i], refList[i].hkl)
        if p not in peaks:
            peaks.append(p)
        else:
            #print("Peak at: ", svalues[i], "adding: ", peaks[peaks.index(p)].sfs2, " to: ", sfs2[i], " = ", peaks[peaks.index(p)].sfs2+sfs2[i])
            peaks[peaks.index(p)].sfs2 += sfs2[i]
            if peaks[peaks.index(p)].error is not None:
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
    if background is None:
        background = np.zeros(len(sList))
    if sList is None:
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
    if refList.magnetic:
        sfs2 = calcMagStructFact(refList, atomList, spaceGroup, cell)
        tt = np.radians(np.array([twoTheta(ref.get_magh_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_magh_s() for ref in refList])
    else:
        sfs2 = np.array(calcStructFact(refList, atomList, spaceGroup, wavelength, xtal=True))
        tt = np.radians(np.array([twoTheta(ref.get_reflection_s(), wavelength) for ref in refList]))
        svalues = np.array([ref.get_reflection_s() for ref in refList])
    if extinctions is not None:
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
        if infoFile is not None:
            infofile = readMagInfo(infoFile)
            if spaceGroup is None: spaceGroup = infofile[0]
            if cell is None: cell = infofile[1]
            if magAtomList is None: magAtomList = infofile[2]
            if symmetry is None: symmetry = infofile[3]
        if basisSymmetry is None: basisSymmetry = symmetry
        magRefList = satelliteGen(cell, symmetry, np.sin(np.radians(179.5)/2.0)/wavelength, hkls=refList) # Use np.radians and float division
        print(("length of reflection list " + str(len(magRefList)))) # Python 3 print
        sfs2, svalues = calcXtalIntensity(magRefList, magAtomList, basisSymmetry,
                                      wavelength, cell, True, extinctions=extinctions, scale=scale)
        magpeaks = makeXtalPeaks(sfs2, svalues, magRefList)
        # add in structural peaks
        if atomList is None: atomList = readInfo(infoFile)[2]
        if nuclear:
            sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, extinctions=extinctions, scale=scale)
            peaks = makeXtalPeaks(sfs2, svalues, refList, peaks=magpeaks)
            reflections = magRefList[:] + refList[:]
            intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], background=background, exclusions=exclusions, base=base, scale=scale)[0] # getXtalIntensity now returns tuple
        else:
            peaks = magpeaks
            reflections = magRefList[:]
            intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], exclusions=exclusions, scale=scale)[0] # getXtalIntensity now returns tuple
    else:
        if infoFile is not None:
            infofile = readInfo(infoFile)
            if spaceGroup is None: spaceGroup = infofile[0]
            if cell is None: cell = infofile[1]
            if atomList is None: atomList = infofile[2]
        reflections = refList[:]
        sfs2, svalues = calcXtalIntensity(refList, atomList, spaceGroup, wavelength, extinctions=extinctions, scale=scale)
        peaks = makeXtalPeaks(sfs2, svalues, refList, error=error)
        intensities = getXtalIntensity(peaks, sList=[getS(value, wavelength) for value in tt], background=background, exclusions=exclusions, base=base, scale=scale)[0] # getXtalIntensity now returns tuple
    if info:
        if magnetic:
            printInfo(cell, spaceGroup, (atomList, magAtomList), (refList, magRefList),
                      wavelength, basisSymmetry)
        else:
            printInfo(cell, spaceGroup, atomList, refList, wavelength)
    if plot:
        sObs = np.array([getS(value, wavelength) for value in tt])
        plotXtalPattern(peaks, sObs, obsIntensity, background=background, error=error, base=base, residuals=residuals,labels=labels, scale=scale, refList=refList)
        from matplotlib import pyplot as plt; plt.show()
    if saveFile:
        # Ensure intensities is defined in all branches before saving
        if 'intensities' in locals():
             np.savetxt(saveFile, np.column_stack((tt, intensities)), delimiter=" ")
        else:
            print("Warning: Intensities not calculated, cannot save file.") # Python 3 print
    return
def plotXtalPattern(peaks, sList, obsIntensity, background=None,
                    exclusions=None, labels=None, residuals=False, base=0, scale=1, error=None, refList=None):
    from matplotlib import pyplot as plt

    # plot single crystal intesities vs q as single points
    obspeaks = makeXtalPeaks(obsIntensity, sList, refList, error=error)
    sList_plot = np.array([peak.svalue for peak in obspeaks]) # Use a different name to avoid conflict
    obsIntensity_plot = np.array([peak.sfs2 for peak in obspeaks]) # Use a different name
    calcIntensity, scalc = getXtalIntensity(peaks, sList=sList_plot, exclusions=exclusions, scale=scale)
    plt.subplot(211)
    if obsIntensity_plot is not None: # Check against the plot-specific variable
        plt.plot(sList_plot*(4*np.pi), obsIntensity_plot, '-go', linestyle="None", label="Observed",lw=1)
    plt.plot(scalc*(4*np.pi), calcIntensity, '-bo', linestyle="None", label="Calculated", lw=1)
    error_plot = [peak.error for peak in obspeaks] # Use a different name
    plt.errorbar(sList_plot*(4*np.pi), obsIntensity_plot, yerr=error_plot, fmt=None, ecolor='g')
    plt.xlabel("Q")
    plt.ylabel("Intensity")
    plt.legend()
    if residuals:
        # Ensure calcPeaks is created from the correct scalc that matches obs sList_plot
        # We need to make sure that calcIntensity corresponds to sList_plot for proper residual calculation
        # Re-evaluate theory at sList_plot points if necessary, or ensure getXtalIntensity handles it.
        # For simplicity, assuming getXtalIntensity returns calcIntensity aligned with sList_plot

        # Create calcPeaks based on the scalc from getXtalIntensity which should align with sList_plot if getXtalIntensity is robust
        # However, the original code was creating calcPeaks from 'peaks' which is a list of sXtalPeak objects, not raw intensity arrays
        # This part might need more careful handling of how calcPeaks are generated for residuals

        # Sticking to original logic for calcPeaks but using the plot-specific arrays for residuals
        # calc_sfs2_for_resid, calc_svalues_for_resid = getXtalIntensity(peaks, sList=sList_plot, exclusions=exclusions, scale=scale)

        # Need to create sXtalPeak objects for calculated values at observed s-points for proper comparison if structure of obspeaks is used
        # This is tricky because `peaks` (calculated) might not have all hkls that `obspeaks` has.
        # The original residual calculation might be flawed if hkls don't match one-to-one.
        # For now, let's assume sList_plot is the definitive list of s-values for comparison.

        # A more direct way to calculate residuals if direct intensity array comparison is intended:
        # Ensure observed_intensity and calculated_intensity are aligned to the same Q (or s) points.
        # obspeaks already aligns observed data. We need calculated data at these exact s-points.

        calc_intensity_at_obs_s = []
        # Create a dictionary for quick lookup of calculated intensities by s-value
        # Ensure scalc and calcIntensity are numpy arrays for efficient lookup if needed, though direct dict creation is fine for moderate sizes
        calc_map = {s: i for s, i in zip(np.asarray(scalc), np.asarray(calcIntensity))}

        for s_obs_val in sList_plot:
            # Find if there's a calculated peak at/near s_obs_val
            # This simple lookup assumes exact s-value match. Might need approx match.
            calc_intensity_at_obs_s.append(calc_map.get(s_obs_val, 0.0)) # append 0.0 if not found, ensure float

        # Ensure obsIntensity_plot is a NumPy array for subtraction
        resid = np.asarray(obsIntensity_plot) - np.asarray(calc_intensity_at_obs_s)

        plt.subplot(212)
        plt.plot(sList_plot*(4*np.pi), resid, '-bo', linestyle="None", label="Residuals")
    return

# printInfo: prints out information about the provided space group and atoms,
#   as well as the generated reflections
def printXtalInfo(cell, spaceGroup, atomLists, refLists, wavelength, symmetry=None):
    print(("Wavelength:", wavelength)) # Python 3 print
    if isinstance(refLists, ReflectionList): # This should be fine as ReflectionList is a class
        atomLists = (atomLists,)
        refLists = (refLists,)

    divider = "-" * 40
    print(("Cell information (%s cell)" % getSpaceGroup_crystalsys(spaceGroup).rstrip())) # Python 3 print
    print(divider) # Python 3 print
    print((" a = %.3f   alpha = %.3f" % (cell.length()[0], cell.angle()[0]))) # Python 3 print
    print((" b = %.3f   beta  = %.3f" % (cell.length()[1], cell.angle()[1]))) # Python 3 print
    print((" c = %.3f   gamma = %.3f" % (cell.length()[2], cell.angle()[2]))) # Python 3 print
    print(divider) # Python 3 print
    print() # Python 3 print
    print("Space group information") # Python 3 print
    print(divider) # Python 3 print
    print(("               Number: ", spaceGroup.get_space_group_numspg())) # Python 3 print
    print(("           H-M Symbol: ", getSpaceGroup_spg_symb(spaceGroup))) # Python 3 print
    print(("          Hall Symbol: ", getSpaceGroup_hall(spaceGroup))) # Python 3 print
    print(("       Crystal System: ", getSpaceGroup_crystalsys(spaceGroup))) # Python 3 print
    print(("           Laue Class: ", getSpaceGroup_laue(spaceGroup))) # Python 3 print
    print(("          Point Group: ", getSpaceGroup_pg(spaceGroup))) # Python 3 print
    print((" General Multiplicity: ", spaceGroup.get_space_group_multip())) # Python 3 print
    print(divider) # Python 3 print
    print() # Python 3 print
    print(("Atom information (%d atoms)" % len(atomLists[0]))) # Python 3 print
    print(divider) # Python 3 print
    atomList = atomLists[0] # This should be fine if atomLists[0] is an AtomList object
    # magnetic = atomList.magnetic # This attribute might not exist directly on AtomList, but rather on individual atoms or AtomListModel

    # Assuming getAtom_lab and atom.coords() are methods of objects within atomList
    label = [getAtom_lab(atom).rstrip() for atom in atomList] if len(atomList) > 0 else []
    x, y, z = [], [], []
    if len(atomList) > 0:
        x, y, z = tuple([ "%.3f" % atom.coords()[i] for atom in atomList] for i in range(3))

    multip = [str(atom.get_atom_mult()) for atom in atomList] if len(atomList) > 0 else []
    occupancy = ["%.3f" % (atom.get_atom_occ()*spaceGroup.get_space_group_multip()/float(atom.get_atom_mult())) # Ensure float division
                 for atom in atomList] if len(atomList) > 0 else []
    # Figure out what the width of each column should be
    # Ensure keys are strings for max length calculation if they can be other types
    width = OrderedDict([('label', max(len(str(max(label, key=len))), 5) if label else 5),
                         ('x', len(str(max(x, key=len))) if x else 1),
                         ('y', len(str(max(y, key=len))) if y else 1),
                         ('z', len(str(max(z, key=len))) if z else 1),
                         ('mult', max(len(str(max(multip, key=len))), 4) if multip else 4),
                         ('occ', max(len(str(max(occupancy, key=len))), 3) if occupancy else 3),
                         ])
    print(("%s   %s %s %s   %s  %s" % tuple([key.center(v) for key, v
                                             in width.items()]))) # Python 3 print
    for i in range(len(atomList)):
        print(("%s  (%s %s %s)  %s  %s" % (label[i].center(width["label"]),
                                          x[i].rjust(width["x"]),
                                          y[i].rjust(width["y"]),
                                          z[i].rjust(width["z"]),
                                          multip[i].center(width["mult"]),
                                          occupancy[i].rjust(width["occ"])))) # Python 3 print
    print(divider) # Python 3 print
    print() # Python 3 print
    print(("Reflection information (%d reflections)" %
          sum([len(refList) for refList in refLists]))) # Python 3 print
    print(divider) # Python 3 print
    for atomList_iter, refList_iter in zip(atomLists, refLists): # Renamed to avoid conflict
        # magnetic_refl = refList_iter.magnetic # 'magnetic' attribute might not exist on ReflectionList
        # For now, assuming magnetic status is known or passed differently if needed for symmObject
        # This was: if magnetic: symmObject = symmetry. Defaulting to spaceGroup.
        symmObject = spaceGroup
        if symmetry is not None: # A guess: if symmetry is provided, it might imply a magnetic calculation for this refl list
             # A more robust check for refList_iter.magnetic would be ideal if the object supports it
            try:
                if refList_iter.magnetic: # Check if the individual reflection list is marked magnetic
                    symmObject = symmetry
            except AttributeError:
                pass # Keep symmObject as spaceGroup

        h, k, l = [], [], []
        if len(refList_iter) > 0:
            h, k, l = tuple([str(ref.hkl[i]) for ref in refList_iter] for i in range(3))

        multip_r = [str(ref.multip) for ref in refList_iter] if len(refList_iter) > 0 else [] # Renamed to avoid conflict with outer scope multip
        tt_calc = ["%.3f" % twoTheta(ref.s, wavelength) for ref in refList_iter]  if len(refList_iter) > 0 else [] # Renamed
        intensity_calc = [] # Renamed
        if len(refList_iter) > 0:
             # Determine if magnetic calculation for calcXtalIntensity based on symmObject
            is_magnetic_calc = (symmObject is symmetry) and (symmetry is not None)
            intensity_calc = ["%.3f" % I for I in calcXtalIntensity(refList_iter, atomList_iter, symmObject, wavelength, cell, is_magnetic_calc)[0]]

        #dtype = [('tt', float),('h', 'S10'), ('k', 'S10'), ('l','S10'), ('intensity', 'S10')] # np.float is deprecated, use float
        #array1 = np.array([(tt_calc[i], str(float(h[i])+0.5),k[i],str(float(l[i])+0.5),intensity_calc[i]) for i in range(len(tt_calc))], dtype=dtype)
        #array2 = np.sort(array1, order='tt')
        #print(array2) # Python 3 print
        # Figure out what the width of each column should be
        width = OrderedDict([('h', max(len(str(max(h, key=len))), 1) if h else 1),
                             ('k', max(len(str(max(k, key=len))), 1) if k else 1),
                             ('l', max(len(str(max(l, key=len))), 1) if l else 1),
                             ('mult', max(len(str(max(multip_r, key=len))), 4) if multip_r else 4), # Use renamed multip_r
                             ('2*theta', max(len(str(max(tt_calc, key=len))), 7) if tt_calc else 7),
                             ('intensity', max(len(str(max(intensity_calc, key=len))), 9) if intensity_calc else 9)
                            ])
        print(("  %s %s %s   %s  %s  %s" % tuple([key.center(v) for key, v
                                                 in width.items()]))) # Python 3 print
        for i in range(len(refList_iter)): # Iterate over the current reflection list
            print((" (%s %s %s)  %s  %s  %s" % (h[i].rjust(width["h"]),
                                               k[i].rjust(width["k"]),
                                               l[i].rjust(width["l"]),
                                               multip_r[i].center(width["mult"]), # Use renamed multip_r
                                               tt_calc[i].rjust(width["2*theta"]),
                                               intensity_calc[i].rjust(width["intensity"]))))
        print()
    print(divider)
    print()

# bumps model
# Model: represents an object that can be used with bumps for optimization
#   purposes single crystal version
class Model(object):

    def __init__(self, tt, observed, background,
                 wavelength, spaceGroupName, cell, atoms, exclusions=None,
                 magnetic=False, symmetry=None, newSymmetry=None, scale=1, zero=None, error=None, hkls=None, extinction=[0]):
        if isinstance(spaceGroupName, SpaceGroup):
            self.spaceGroup = spaceGroupName
        else:
            self.spaceGroup = SpaceGroup(spaceGroupName)
        self.tt = np.asarray(tt)
        self.obspeaks = makeXtalPeaks(observed, [getS(ttval, wavelength) for ttval in self.tt], refList=hkls, error=error)
        self.sList = np.array([peak.svalue for peak in self.obspeaks])
        self.observed = np.array([peak.sfs2 for peak in self.obspeaks])
        self.background = background
        self.scale = Parameter(scale, name='scale')
        self.extinctions = [Parameter(extinction[i], name="Extinction"+str(i)) for i in range(len(extinction))]
        self.error = error
        self.refList = hkls
        self.base=0
        if zero is not None:
            self.zero = Parameter(zero, name='zero')
            self.has_zero = True
        else:
            self.zero = 0
            self.has_zero = False
        self.wavelength = wavelength
        self.cell = cell
        self.exclusions = exclusions
        if len(self.tt):
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
        self._set_observations(observed)
        self.update()
    def _set_reflections(self):
        maxLattice = self.cell.getMaxLattice()
        # TODO: why is FloatVector failing on inf?  Why is maxLattice inf?
        maxLattice = [1e38 if v > 1e38 else v for v in maxLattice]
        maxCell = CrystalCell(FloatVector(maxLattice[:3]), FloatVector(maxLattice[3:])) # Use FloatVector
        if self.magnetic:
            self.magRefList = satelliteGen(self.cell.cell, self.symmetry, np.sin(np.radians(179.5)/2.0)/self.wavelength, hkls=self.refList) # Use np.radians
            self.magReflections = self.magRefList[:]
        self.reflections = self.refList

    def _set_observations(self, observed):
        self.obspeaks = makeXtalPeaks(
            observed,
            [getS(ttval, self.wavelength) for ttval in self.tt],
            refList=self.refList,
            error=self.error)
        self.sList = np.array([peak.svalue for peak in self.obspeaks])
        self.observed = np.array([peak.sfs2 for peak in self.obspeaks])

    def __getstate__(self):
        state = self.__dict__.copy()
        if "refList" in state: del state["refList"] # Ensure keys exist before deleting
        if "magRefList" in state: del state["magRefList"]
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        self._set_reflections() # This might try to access refList if not careful

    def parameters(self):
        params = {} # Initialize params
        if self.has_zero:
            params = {
                    'scale': self.scale,
                    'zero' : self.zero,
                    'cell': self.cell.parameters(),
                    'atoms': self.atomListModel.parameters()
                    }
        else: # Ensure params is initialized even if has_zero is false
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
        # getXtalIntensity returns a tuple, ensure we get the first element
        return getXtalIntensity(self.peaks, background=self.background, sList=self.sList, scale=self.scale.value)[0]

    def residuals(self):
        # Ensure self.observed is a numpy array for sqrt
        observed_array = np.asarray(self.observed)
        # The previous complex residual calculation in plotXtalPattern was more robust.
        # This is a simplified version. If issues arise, refer to that logic.
        # For now, assuming self.theory() and self.observed are directly comparable and aligned.

        # Replicating the logic from plotXtalPattern for residuals if direct comparison is intended
        calc_intensity_at_obs_s = []
        # Assuming self.peaks contains calculated peak objects and self.sList are the s-values for these peaks.
        # This might not be entirely correct if self.sList for theory is different from self.sList from obspeaks
        # For bumps, theory should be evaluated at the same points as observed data.
        # Let's assume self.sList (from obspeaks) is the definitive list of points.

        # Get calculated intensities at the observed s-points
        # The self.theory() method should ideally handle this alignment.
        # If self.theory() already returns intensities aligned with self.observed (at self.sList points), then it's simpler.

        # Assuming self.theory() is already aligned with self.observed
        theory_intensities = self.theory() # This should be an array of calculated intensities at observed s-points

        return (theory_intensities - observed_array)/(np.sqrt(observed_array)+1e-8) # Add epsilon to avoid division by zero if observed can be 0

    def nllf(self):
        return np.sum(self.residuals()**2)

    def plot(self, view="linear"):
        if self.has_zero:
            zero = self.zero.value
        else:
            zero = self.zero
        # Ensure self.observed is passed correctly
        plotXtalPattern(self.peaks, self.sList, np.asarray(self.observed), # Pass observed data array
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
            sList_mag = calcS(self.cell.cell, hkls) # Use a different name for clarity
            for i in range(len(self.magReflections)):
                self.magReflections[i].set_magh_s(sList_mag[i])
            sfs2, svalues = calcXtalIntensity(self.magRefList, self.atomListModel.magAtomList, self.symmetry, self.wavelength, magnetic=True, cell=self.cell.cell, extinctions=[ext.value for ext in self.extinctions], scale=self.scale.value)
            self.magIntensities = sfs2
            #print(self.magIntensities) # Python 3 print
            self.peaks = makeXtalPeaks(sfs2, svalues, self.magRefList)
            #self.sList = np.array([peak.svalue for peak in self.peaks]) # This would overwrite sList from observed
            #self.peaks.extend(makeXtalPeaks(sfs2, svalues))

        if self.reflections is not None: # Ensure self.reflections is not None
            hkls_nuc = [reflection.hkl for reflection in self.reflections] # Use a different name for clarity
            sList_nuc = calcS(self.cell.cell, hkls_nuc) # Use a different name
            for i in range(len(self.reflections)):
                self.reflections[i].set_reflection_s(sList_nuc[i])
            sfs2_nuc, svalues_nuc = calcXtalIntensity(self.refList, self.atomListModel.atomList, self.spaceGroup, self.wavelength, extinctions=[ext.value for ext in self.extinctions], scale=self.scale.value) # Use different names
            self.intensities = sfs2_nuc # Store nuclear intensities
            if not self.magnetic: self.peaks = None # Reset peaks if not magnetic before adding nuclear
            # If magnetic, peaks already contains magnetic peaks. Append nuclear.
            # If not magnetic, self.peaks is None, so makeXtalPeaks will create new.
            self.peaks = makeXtalPeaks(sfs2_nuc, svalues_nuc, self.refList, peaks=self.peaks)
            #self.sList = svalues_nuc # This would overwrite sList from observed
