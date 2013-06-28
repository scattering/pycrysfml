import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os

import numpy as np
import hklGen as H

np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"Al2O3 Background.BGR")
observedFile = os.path.join(DATAPATH,"Al2O3.dat")
infoFile = os.path.join(DATAPATH,"Al2O3.cif")

#(spaceGroup, cell, wavelength, sMin, sMax) = inputInfo()
#(spaceGroup, cell, wavelength, sMin, sMax, atoms) = testInfo()
(spaceGroup, crystalCell, atoms) = H.readFile(infoFile)
spaceGroup.xtalSystem = spaceGroup.xtalSystem.rstrip()
wavelength = 1.5403

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
backg = H.LinSpline(backgFile)
tt = np.linspace(3, 167.75, 3296)
# LuFeO3 exclusions:
#    exclusions = np.array([[0,3],[37.3,39.1],[43.85,45.9],[64.25,66.3],
#                           [76.15,79.8],[81.7,83.1],[89.68,99.9],[109.95,111.25],
#                           [115.25,118.45],[133.95,141.25],[156.7,180]])
exclusions = None
data = np.loadtxt(observedFile, dtype=float, skiprows=1)
observed = data.flatten()[:len(tt)]
tt, observed = H.removeRange(tt, exclusions, observed)
#    observed = data.flatten()[:len(tt)]
#    cell = HexagonalCell(4.7698, 13.0243)
#    cell = HexagonalCell(5.965, 11.702)

def fit():
    cell = H.makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
#    cell.c.pm(0.5)
    m = H.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
              atoms, exclusions)
    m.u.range(0,1)
    m.v.range(-1,0)
    m.w.range(0,10)
    m.scale.range(0,10)
    for atomModel in m.atomListModel.atomModels:
        atomModel.B.range(0, 10)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    # LuFeMnO3 occupancy:
#    m.atomListModel["Fe1"].occ.range(0, 1)
#    m.atomListModel["Mn1"].occ.range(0, 1)
#    m.atomListModel["Fe1"].occ = 1 - m.atomListModel["Mn1"].occ
#    P = m.atomListModel["Fe1"].occ
#    M = FitProblem(m, constraints=lambda:(P.value>1)*(1000+P.value**4))
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    import pylab
    for atom in H.deconstruct_dv(atoms.atoms, H.Atom): 
        print atom.label, atom.occupancy, atom.multip

    sMin, sMax = H.getS(3, wavelength), H.getS(167.8, wavelength)

    refList = H.hklGen(spaceGroup, crystalCell, wavelength, sMin, sMax, True)
    reflections = list(refList.reflections.data(H.Reflection))
    H.printReflections(reflections, spaceGroup, wavelength, sMin, sMax)
    intensities = H.calcIntensity(refList, atoms, spaceGroup, wavelength)
    g = H.makeGaussians(reflections,[.347, -.278, .166], intensities, 1, wavelength)
    H.plotPattern(g, backg, tt, observed, H.twoTheta(sMin, wavelength),
                  H.twoTheta(sMax, wavelength), .01, exclusions, labels="hkl")
    pylab.show()
    return

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()

