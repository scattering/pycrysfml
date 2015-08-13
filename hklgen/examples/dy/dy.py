import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
import sxtal_model as S
from pycrysfml import getSpaceGroup_crystalsys as xtalsys
#import bumps.parameter
#def dont(self, *args, **kw): raise Exception("don't")
#bumps.parameter.OperatorAdd.__init__ = dont

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,r"dy.bac")
observedFile = os.path.join(DATAPATH,r"dy.dat")
infoFile = os.path.join(DATAPATH,r"dy.cfl")
exclusions = []
(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
wavelength = 1.703700
tt, observed, error = H.readIllData(observedFile, "DMC", backgFile)
base_line = min(observed)
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
backg = H.LinSpline(None)
basisSymmetry = copy(symmetry)

def fit():
    cell = Mod.makeCell(crystalCell, xtalsys(spaceGroup))
    cell.a.pm(0.5)
#    cell.b.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg,   1.814691,  -1.482098,   0.447632 , wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=base_line, scale=94.508, eta=0.0, zero=0.04101, error=error, muR=1.280)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,200)
    m.eta.range(0,1)
    m.base.pm(1000)
    m.zero.pm(0.25)
    for atomModel in m.atomListModel.atomModels:
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
                #pass
            #atomModel.phase.range(-np.pi*2, np.pi*2)
    #m.atomListModel["Mn1"].phase.range(0, np.pi*2)
    #m.atomListModel["Mn2"].phase.range(-2*np.pi, 0)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    uvw = [1.814655,  -1.482037,  0.447617 ]
    cell = crystalCell
    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=94.5,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), base=base_line, error=error, residuals=True, muR=1.280)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
