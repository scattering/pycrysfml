import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os

import numpy as np
import hklGen as H

np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"Al2O3 Background.BGR")
observedFile = os.path.join(DATAPATH,"Al2O3.dat")
infoFile = os.path.join(DATAPATH,"Al2O3.cif")

(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
spaceGroup.xtalSystem = spaceGroup.xtalSystem.rstrip()
wavelength = 1.5403

backg = H.LinSpline(backgFile)
ttMin = 3
ttMax = 167.75
ttStep = 0.05
exclusions = None
tt, observed = H.readData(observedFile, kind="y", skiplines=1, start=ttMin,
                          stop=ttMax, step=ttStep)

def fit():
    cell = H.makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
    cell.c.pm(0.5)
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
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    import pylab

    sMin, sMax = H.getS(ttMin, wavelength), H.getS(ttMax, wavelength)
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