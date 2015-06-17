import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = None
observedFile = os.path.join(DATAPATH,r"Mn025_4K.dat")
infoFile = os.path.join(DATAPATH,r"LFO.cfl")
(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
wavelength = 1.540000
ttMin = 3.0
ttMax = 167.750
ttStep = 0.05
exclusions = [[0.00,5.00],[38.25,39.00],[44.00,45.50],[64.75,65.75],[77.5,79.00],[81.75,83.10],[97.00,100.00],[111.50,113.00],[116.00,118.00],[135.50,139.50],[161.00,169.00]]
tt, observed = H.readData(observedFile, kind="y", skiplines=3,
                          start=ttMin, stop=ttMax, step=ttStep,
                          exclusions=exclusions, colstep=2)
base_value=min(observed)
backg = H.LinSpline(backgFile)
basisSymmetry = copy(symmetry)
def fit():
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0.226335,  -0.220125,   0.117269, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, eta=0.11935, scale=0.37163, base=base_value)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,10)
    m.base.pm(500)
    m.eta.range(0,1)
    for atomModel in m.atomListModel.atomModels:
#        atomModel.B.range(0, 10)
#        if (atomModel.atom.multip == atomModel.sgmultip):
#            # atom lies on a general position
#            atomModel.x.pm(0.1)
#            atomModel.y.pm(0.1)
#            atomModel.z.pm(0.1)
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                coeff.range(-10, 10)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    H.diffPattern(infoFile=infoFile, backgroundFile=backgFile, uvw=[ 0.226335,  -0.220125,   0.117269], wavelength=wavelength,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, exclusions=exclusions,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), scale=0.37163, base=base_value)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
