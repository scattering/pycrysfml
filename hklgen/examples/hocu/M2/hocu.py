import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
from pycrysfml import getSpaceGroup_crystalsys as xtalsys
#import bumps.parameter
#def dont(self, *args, **kw): raise Exception("don't")
#bumps.parameter.OperatorAdd.__init__ = dont

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,r"hocu.bac")
observedFile = os.path.join(DATAPATH,r"hocu.dat")
infoFile = os.path.join(DATAPATH,r"hocu.cfl")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
exclusions = []
tt, observed, error = H.readIllData(observedFile, "D1B", backgFile)
base_line = min(observed)
wavelength = 2.524000
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
backg = H.LinSpline(None)
basisSymmetry = copy(symmetry)

def fit():
    cell = Mod.makeCell(crystalCell, xtalsys(spaceGroup))
    cell.a.pm(0.5)
    cell.b.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 1.161020,  -0.658240,   0.29176, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=base_line, scale=9.6286/np.sqrt(2.0), eta=0.1031, zero=0.0197, sxtal=True, error=error)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,60)
    m.eta.range(0,1)
    m.base.pm(500)
    m.zero.pm(0.25)
    for atomModel in m.atomListModel.atomModels:
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    uvw = [ 1.161020,  -0.658240,   0.29176 ]
    cell = crystalCell
    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=9.6286/np.sqrt(2.0),
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), base=base_line, xtal=True, error=error, residuals=True)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
