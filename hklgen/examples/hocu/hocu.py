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
tt, observed = H.readIllData(observedFile, "D1B", backgFile)
base_line = min(observed)
wavelength = 2.524000
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
backg = H.LinSpline(None)
basisSymmetry = copy(symmetry)

def makeBasis(symmetry, basisIndex):
    if (basisIndex == 2):
        symmetry.setNumBasisFunc_ind(0, 2)
        symmetry.setBasis(0, 0, 0, [1,0,0])
        symmetry.setBasis(0, 0, 1, [0,0,1])
        symmetry.setBasis(0, 1, 0, [1,0,0])
        symmetry.setBasis(0, 1, 1, [0,0,1])
def fit():
    #makeBasis(basisSymmetry, basisIndex)
    cell = Mod.makeCell(crystalCell, xtalsys(spaceGroup))
    cell.a.pm(0.5)
    cell.b.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0,0,1, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=base_line, scale=9.6286, sxtal=True)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,60)
    m.base.pm(10)
    for atomModel in m.atomListModel.atomModels:
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    #makeBasis(basisSymmetry, basisIndex)
    uvw = [ 1.161020,  -0.658240,   0.297025 ]
    cell = crystalCell
    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=9.6286,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), base=base_line, xtal=True)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
