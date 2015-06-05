import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os
from copy import copy

import numpy as np
import fswig_hklgen as H
import hkl_model as Mod

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
print DATAPATH
backgFile = os.path.join(DATAPATH,r"background_ge.BGR")
observedFile = os.path.join(DATAPATH,r"lufmo001.gsas")
infoFile = os.path.join(DATAPATH,r"LuFeMnO3_BT1fit.cfl")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]

wavelength = 2.07750
ttMin = 3.0
ttMax = 166.3
ttStep = 0.05
exclusions = [[52,54],[60.5,62.5],[92.5,95],[115.0,117],[125.0,126]]
tt, observed = H.readData(observedFile, kind="y", skiplines=3,
                          start=ttMin, stop=ttMax, step=ttStep,
                          exclusions=exclusions, colstep=2)
backg = H.LinSpline(backgFile)

basisSymmetry = copy(symmetry)
#basisSymmetry = symmetry
#if (basisSymmetry.centerType() == 2):
    ## change to acentric
    #basisSymmetry.setCenterType(1)
#basisSymmetry.setNumIrreps(1)

# Number of the basis from BasIreps (1-6)
#basisIndex = 2
#index = 0
#for magAtom in magAtomList:
    #magAtom.setNumkVectors(1)
    #magAtom.setIrrepNum_ind(0, 1)
    #magAtomList[index] = magAtom
    #index += 1
def makeBasis(symmetry, basisIndex):
    if (basisIndex == 1):
        symmetry.setNumBasisFunc_ind(0, 1)
        symmetry.setBasis(0, 0, 0, [1,2,0])
        symmetry.setBasis(0, 1, 0, [-1,-2,0])
        symmetry.setBasis(0, 2, 0, [-2,-1,0])
        symmetry.setBasis(0, 3, 0, [2,1,0])
        symmetry.setBasis(0, 4, 0, [1,-1,0])
        symmetry.setBasis(0, 5, 0, [-1,1,0])
    elif (basisIndex == 2):
        symmetry.setNumBasisFunc_ind(0, 2)
        symmetry.setBasis(0, 0, 0, [1,0,0])
        symmetry.setBasis(0, 0, 1, [0,0,1])
        symmetry.setBasis(0, 1, 0, [-1,0,0])
        symmetry.setBasis(0, 1, 1, [0,0,1])
        symmetry.setBasis(0, 2, 0, [0,1,0])
        symmetry.setBasis(0, 2, 1, [0,0,1])
        symmetry.setBasis(0, 3, 0, [0,-1,0])
        symmetry.setBasis(0, 3, 1, [0,0,1])
        symmetry.setBasis(0, 4, 0, [-1,-1,0])
        symmetry.setBasis(0, 4, 1, [0,0,1])
        symmetry.setBasis(0, 5, 0, [1,1,0])
        symmetry.setBasis(0, 5, 1, [0,0,1])
    elif (basisIndex == 3):
        symmetry.setNumBasisFunc_ind(0, 6)
        symmetry.setBasis(0, 0, 0, [1,0,0])
        symmetry.setBasis(0, 0, 1, [0,1,0])
        symmetry.setBasis(0, 0, 2, [0,0,1])
        symmetry.setBasis(0, 0, 3, [complex(0.5,0.866),0,0])
        symmetry.setBasis(0, 0, 4, [complex(-0.5,-0.866),complex(-0.5,-0.866),0])
        symmetry.setBasis(0, 0, 5, [0,0,complex(0.5,0.866)])
        symmetry.setBasis(0, 1, 0, [-1,0,0])
        symmetry.setBasis(0, 1, 1, [0,-1,0])
        symmetry.setBasis(0, 1, 2, [0,0,1])
        symmetry.setBasis(0, 1, 3, [complex(-0.5,-0.866),0,0])
        symmetry.setBasis(0, 1, 4, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 1, 5, [0,0,complex(0.5,0.866)])
        symmetry.setBasis(0, 2, 0, [0,complex(-0.5,-0.866),0])
        symmetry.setBasis(0, 2, 1, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 2, 2, [0,0,complex(-0.5,-0.866)])
        symmetry.setBasis(0, 2, 3, [0,-1,0])
        symmetry.setBasis(0, 2, 4, [-1,0,0])
        symmetry.setBasis(0, 2, 5, [0,0,-1])
        symmetry.setBasis(0, 3, 0, [0,complex(0.5,0.866),0])
        symmetry.setBasis(0, 3, 1, [complex(-0.5,-0.866),complex(-0.5,-0.866),0])
        symmetry.setBasis(0, 3, 2, [0,0,complex(-0.5,-0.866)])
        symmetry.setBasis(0, 3, 3, [0,1,0])
        symmetry.setBasis(0, 3, 4, [1,0,0])
        symmetry.setBasis(0, 3, 5, [0,0,-1])
        symmetry.setBasis(0, 4, 0, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 4, 1, [complex(0.5,-0.866),0,0])
        symmetry.setBasis(0, 4, 2, [0,0,complex(-0.5,0.866)])
        symmetry.setBasis(0, 4, 3, [complex(-0.5,0.866),complex(-0.5,0.866),0])
        symmetry.setBasis(0, 4, 4, [0,complex(0.5,-0.866),0])
        symmetry.setBasis(0, 4, 5, [0,0,complex(0.5,-0.866)])
        symmetry.setBasis(0, 5, 0, [complex(-0.5,0.866),complex(-0.5,0.866),0])
        symmetry.setBasis(0, 5, 1, [complex(-0.5,0.866),0,0])
        symmetry.setBasis(0, 5, 2, [0,0,complex(-0.5,0.866)])
        symmetry.setBasis(0, 5, 3, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 5, 4, [0,complex(-0.5,0.866),0])
        symmetry.setBasis(0, 5, 5, [0,0,complex(0.5,-0.866)])
    elif (basisIndex == 4):
        symmetry.setNumBasisFunc_ind(0, 2)
        symmetry.setBasis(0, 0, 0, [1,0,0])
        symmetry.setBasis(0, 0, 1, [0,0,1])
        symmetry.setBasis(0, 1, 0, [1,0,0])
        symmetry.setBasis(0, 1, 1, [0,0,-1])
        symmetry.setBasis(0, 2, 0, [0,1,0])
        symmetry.setBasis(0, 2, 1, [0,0,1])
        symmetry.setBasis(0, 3, 0, [0,1,0])
        symmetry.setBasis(0, 3, 1, [0,0,-1])
        symmetry.setBasis(0, 4, 0, [-1,-1,0])
        symmetry.setBasis(0, 4, 1, [0,0,1])
        symmetry.setBasis(0, 5, 0, [-1,-1,0])
        symmetry.setBasis(0, 5, 1, [0,0,-1])
    elif (basisIndex == 5):
        symmetry.setNumBasisFunc_ind(0, 1)
        symmetry.setBasis(0, 0, 0, [1,2,0])
        symmetry.setBasis(0, 1, 0, [1,2,0])
        symmetry.setBasis(0, 2, 0, [-2,-1,0])
        symmetry.setBasis(0, 3, 0, [-2,-1,0])
        symmetry.setBasis(0, 4, 0, [1,-1,0])
        symmetry.setBasis(0, 5, 0, [1,-1,0])
    elif (basisIndex == 6):
        symmetry.setNumBasisFunc_ind(0, 6)
        symmetry.setBasis(0, 0, 0, [1,0,0])
        symmetry.setBasis(0, 0, 1, [0,1,0])
        symmetry.setBasis(0, 0, 2, [0,0,1])
        symmetry.setBasis(0, 0, 3, [complex(-0.5,-0.866),0,0])
        symmetry.setBasis(0, 0, 4, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 0, 5, [0,0,complex(-0.5,-0.866)])
        symmetry.setBasis(0, 1, 0, [1,0,0])
        symmetry.setBasis(0, 1, 1, [0,1,0])
        symmetry.setBasis(0, 1, 2, [0,0,-1])
        symmetry.setBasis(0, 1, 3, [complex(-0.5,-0.866),0,0])
        symmetry.setBasis(0, 1, 4, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 1, 5, [0,0,complex(0.5,0.866)])
        symmetry.setBasis(0, 2, 0, [0,complex(-0.5,-0.866),0])
        symmetry.setBasis(0, 2, 1, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 2, 2, [0,0,complex(-0.5,-0.866)])
        symmetry.setBasis(0, 2, 3, [0,1,0])
        symmetry.setBasis(0, 2, 4, [1,0,0])
        symmetry.setBasis(0, 2, 5, [0,0,1])
        symmetry.setBasis(0, 3, 0, [0,complex(-0.5,-0.866),0])
        symmetry.setBasis(0, 3, 1, [complex(0.5,0.866),complex(0.5,0.866),0])
        symmetry.setBasis(0, 3, 2, [0,0,complex(0.5,0.866)])
        symmetry.setBasis(0, 3, 3, [0,1,0])
        symmetry.setBasis(0, 3, 4, [1,0,0])
        symmetry.setBasis(0, 3, 5, [0,0,-1])
        symmetry.setBasis(0, 4, 0, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 4, 1, [complex(-0.5,0.866),0,0])
        symmetry.setBasis(0, 4, 2, [0,0,complex(-0.5,0.866)])
        symmetry.setBasis(0, 4, 3, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 4, 4, [0,complex(-0.5,0.866),0])
        symmetry.setBasis(0, 4, 5, [0,0,complex(-0.5,0.866)])
        symmetry.setBasis(0, 5, 0, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 5, 1, [complex(-0.5,0.866),0,0])
        symmetry.setBasis(0, 5, 2, [0,0,complex(0.5,-0.866)])
        symmetry.setBasis(0, 5, 3, [complex(0.5,-0.866),complex(0.5,-0.866),0])
        symmetry.setBasis(0, 5, 4, [0,complex(-0.5,0.866),0])
        symmetry.setBasis(0, 5, 5, [0,0,complex(0.5,-0.866)])
    print symmetry.numBasisFunc()
    print symmetry.getBasis(0, 0, 0)
    print symmetry.getBasis(0, 0, 1)
    print symmetry.getBasis(0, 1, 0)
    print symmetry.getBasis(0, 1, 1)
    print symmetry.getBasis(0, 2, 0)
    print symmetry.getBasis(0, 2, 1)
    print symmetry.getBasis(0, 3, 0)
    print symmetry.getBasis(0, 3, 1)
    print symmetry.getBasis(0, 4, 0)
    print symmetry.getBasis(0, 4, 1)
    print symmetry.getBasis(0, 5, 0)
    print symmetry.getBasis(0, 5, 1)
def fit():
    #makeBasis(basisSymmetry, basisIndex)

    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem())
    cell.a.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,10)
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
#            atomModel.phase.range(0,1)
    # vary Fe/Mn atom positions but keep them on the special site x,0,z
#    m.atomListModel["Fe1"].x.pm(0.1)
#    m.atomListModel["Fe1"].z.pm(0.1)
#    m.atomListModel["Mn1"].x = m.atomListModel["Fe1"].x
#    m.atomListModel["Mn1"].z = m.atomListModel["Fe1"].z
    for i in xrange(len(m.atomListModel["Fe1"].coeffs)):
        m.atomListModel["Mn1"].coeffs[i] = m.atomListModel["Fe1"].coeffs[i]
#    m.atomListModel["Mn1"].phase = m.atomListModel["Fe1"].phase
    # Occupancy:
#    m.atomListModel["Fe1"].occ.range(0, 1)
#    m.atomListModel["Mn1"].occ.range(0, 1)
#    m.atomListModel["Mn1"].occ = 1 - m.atomListModel["Fe1"].occ
    
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    #index = 0
    #for magAtom in magAtomList:
        #magAtom.setBasis_ind(0, 0, 3)
        #magAtom.setBasis_ind(0, 1, 0.1)
        #magAtomList[index] = magAtom
        #index += 1
    #makeBasis(basisSymmetry, basisIndex)
    H.diffPattern(infoFile=infoFile, backgroundFile=backgFile, wavelength=wavelength,
                  ttMin=ttMin, ttMax=ttMax, ttStep=ttStep, exclusions=exclusions,
                  basisSymmetry=basisSymmetry, magAtomList=magAtomList,
                  magnetic=True, info=True, plot=True,
                  observedData=(tt,observed))
    #print symmetry.numBasisFunc()
    #print symmetry.getBasis(0, 0, 0)
    #print symmetry.getBasis(0, 0, 1)
    #print symmetry.getBasis(0, 1, 0)
    #print symmetry.getBasis(0, 1, 1)
    #print symmetry.getBasis(0, 2, 0)
    #print symmetry.getBasis(0, 2, 1)
    #print symmetry.getBasis(0, 3, 0)
    #print symmetry.getBasis(0, 3, 1)
    #print symmetry.getBasis(0, 4, 0)
    #print symmetry.getBasis(0, 4, 1)
    #print symmetry.getBasis(0, 5, 0)
    #print symmetry.getBasis(0, 5, 1)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
