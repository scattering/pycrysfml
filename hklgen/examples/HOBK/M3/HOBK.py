import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
from pycrysfml import getSpaceGroup_crystalsys as xtalsys
from fswig_hklgen import satelliteGen
import bumps.parameter
#def dont(self, *args, **kw): raise Exception("don't")
#bumps.parameter.OperatorAdd.__init__ = dont

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,r"hobk_bas.bac")
observedFile = os.path.join(DATAPATH,r"hobk.dat")
infoFile = os.path.join(DATAPATH,r"hobk1.cfl")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]

wavelength = 2.524000
ttMin = 10.010000228881836
ttMax = 89.81000518798828
ttStep = 0.20000000298
exclusions = []

#print H.getMaxNumRef(H.getS(ttMax, wavelength), crystalCell.volume)

tt, observed = H.readIllData(observedFile, "D1B", backgFile)
#print H.getMaxNumRef(H.getS(ttMax, wavelength), crystalCell.volume)
#sMin, sMax = H.getS(ttMin, wavelength), H.getS(ttMax, wavelength)
#refList = H.hklGen(spaceGroup, crystalCell, sMin, sMax, True)
#print len(refList)
#print len(H.satelliteGen(crystalCell, symmetry, sMax, refList))
backg = H.LinSpline(None)
#print backg
#basisSymmetry = copy(symmetry)
basisSymmetry = symmetry
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
    #print len(H.satelliteGen(cell.cell, symmetry, float(H.getS(ttMax, wavelength))))
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=6512, scale=59.143)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(59,60)
    m.base.range(6510,6514)
    for atomModel in m.atomListModel.atomModels:
#        atomModel.B.range(0, 10)
#        if (atomModel.atom.multip == atomModel.sgmultip):
#            # atom lies on a general position
#            atomModel.x.pm(0.1)
#            atomModel.y.pm(0.1)
#            atomModel.z.pm(0.1)
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
#            atomModel.phase.range(0,1)
    # vary Fe/Mn atom positions but keep them on the special site x,0,z
#    m.atomListModel["Fe1"].x.pm(0.1)
#    m.atomListModel["Fe1"].z.pm(0.1)
#    m.atomListModel["Mn1"].x = m.atomListModel["Fe1"].x
#    m.atomListModel["Mn1"].z = m.atomListModel["Fe1"].z
    #for i in xrange(len(m.atomListModel["Fe1"].coeffs)):
        #m.atomListModel["Mn1"].coeffs[i] = m.atomListModel["Fe1"].coeffs[i]
#    m.atomListModel["Mn1"].phase = m.atomListModel["Fe1"].phase
    # Occupancy:
#    m.atomListModel["Fe1"].occ.range(0, 1)
#    m.atomListModel["Mn1"].occ.range(0, 1)
#    m.atomListModel["Mn1"].occ = 1 - m.atomListModel["Fe1"].occ
    
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    Ho = magAtomList[0]
    Ni = magAtomList[1]
    Ho.setBasis_ind(0,0,0.127)
    Ho.setBasis_ind(0,1,8.993)
    Ni.setBasis_ind(0,0,0.584)
    Ni.setBasis_ind(0,1,-1.285)
    magAtomList[0] = Ho
    magAtomList[1] = Ni
    #makeBasis(basisSymmetry, basisIndex)
    uvw = [1.548048,-0.988016,0.338780]
    cell = crystalCell
    #H.readMagInfo(infoFile)[3]
#    reflist = satelliteGen(cell, H.readMagInfo(infoFile)[3], ttMax)
    H.diffPattern(infoFile=infoFile, uvw=uvw, cell=cell, scale=59.143,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, wavelength = wavelength,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed), base=6512)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
