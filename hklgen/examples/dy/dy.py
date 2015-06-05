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
backgFile = os.path.join(DATAPATH,r"dy.bac")
observedFile = os.path.join(DATAPATH,r"dy.dat")
infoFile = os.path.join(DATAPATH,r"dy.pcr")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
wavelength = 1.703700
tt, observed = H.readIllData(observedFile, "DMC", backgFile)
print tt, observed
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
exclusions = []
backg = H.LinSpline(None)
#print backg
basisSymmetry = copy(symmetry)
#basisSymmetry = symmetry
#if (basisSymmetry.centerType() == 2):
    ## change to acentric
    #basisSymmetry.setCenterType(1)
#basisSymmetry.setNumIrreps(1)

## Number of the basis from BasIreps (1-6)
#basisIndex = 2
#index = 0
#for magAtom in magAtomList:
    #magAtom.setNumkVectors(1)
    #magAtom.setIrrepNum_ind(0, 1)
    #magAtomList[index] = magAtom
    #index += 1


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
    m = Mod.Model(tt, observed, backg, 1.809863,  -1.476814,   0.446315, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=6512, scale=94.064, sxtal=True)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(11,60)
    m.base.range(6510,6514)
    for atomModel in m.atomListModel.atomModels:
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    #Ho = magAtomList[0]
    #Ni = magAtomList[1]
    #Ho.setBasis_ind(0,0,0.127)
    #Ho.setBasis_ind(0,1,8.993)
    #Ni.setBasis_ind(0,0,0.584)
    #Ni.setBasis_ind(0,1,-1.285)
    #magAtomList[0] = Ho
    #magAtomList[1] = Ni
    #makeBasis(basisSymmetry, basisIndex)
    uvw = [1.809863,  -1.476814,   0.446315]
    cell = crystalCell
    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=94.064,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), base=6512, xtal=True)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
