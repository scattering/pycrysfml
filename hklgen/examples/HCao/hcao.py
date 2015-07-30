import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
from copy import copy
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
import sxtal_model as S
#import bumps.parameter
#def dont(self, *args, **kw): raise Exception("don't")
#bumps.parameter.OperatorAdd.__init__ = dont

np.seterr(divide="ignore", invalid="ignore")

DATAPATH = os.path.dirname(os.path.abspath(__file__))
observedFile = os.path.join(DATAPATH,r"MnWO4_nuclear_5K.int")
infoFile = os.path.join(DATAPATH,r"AF4_5k.cfl")
magfile = os.path.join(DATAPATH, r"MnWO4_magneticAF4_5K.int")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
spaceGroup, crystalCell, atomList = H.readInfo(infoFile)
exclusions = []
# return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
wavelength, refList, sfs2, error = S.readIntFile(observedFile, kind="int", cell=crystalCell)
magrefList, magsfs2, magerror = S.readMagIntFile(magfile, cell=crystalCell)
tt = [H.twoTheta(H.calcS(crystalCell, ref.hkl), wavelength) for ref in refList]
sfs2 = list(sfs2)
sfs2.extend(magsfs2)
sfs2 = np.array(sfs2)
tt.extend([H.twoTheta(H.calcS(crystalCell, ref.hkl), wavelength) for ref in magrefList])
error = list(error)
error.extend(list(magerror))
error2 = [item for item in error]
backg = None
basisSymmetry = copy(symmetry)

def fit():
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)
    #cell.a.pm(5.0)
    #cell.c.pm(5.0)
    m = S.Model(tt, sfs2, backg, wavelength, spaceGroup, cell, 
                (atomList, magAtomList), exclusions, magnetic=True, 
                symmetry=symmetry, newSymmetry=basisSymmetry, scale=1.00, 
                error=error, hkls=refList, extinction=[0,0,0,0])
    m.scale.range(0,2000)
    #m.base.pm(1000)
    for ext in m.extinctions:
        ext.range(0,100.0)
    for atomModel in m.atomListModel.atomModels:
        #atomModel.x.range(0,1)
        #atomModel.y.range(0,1)
        #atomModel.z.range(0,1)
        #atomModel.B.range(0,10)
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                coeff.range(-20, 20)
                #coeff.range(0,5)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = crystalCell
    S.diffPatternXtal(infoFile=infoFile, cell=cell, scale=186.0, tt=tt, 
                      obsIntensity=sfs2, wavelength=wavelength,
                      plot=True, residuals=True, error=error2, 
                      info=True, base=0, refList=refList, extinctions=[ 3.607,2.156,4.335,0.4110,0.2721,0.0], magAtomList=magAtomList, magnetic=True)
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
