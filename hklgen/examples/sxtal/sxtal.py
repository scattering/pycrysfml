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
backgFile = None#os.path.join(DATAPATH,r"hobk_bas.bac")
observedFile = os.path.join(DATAPATH,r"NC35L.dat")
infoFile = os.path.join(DATAPATH,r"NC35L.cfl")

(spaceGroup, crystalCell, magAtomList, symmetry) = H.readMagInfo(infoFile)
atomList = H.readInfo(infoFile)[2]
exclusions = []
# return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
wavelength, refList, sfs2, error, tt, four_circle = S.readIntFile(observedFile)
backg = H.LinSpline(None)
basisSymmetry = copy(symmetry)

def fit():
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(5.0)
    cell.b.pm(5.0)
    cell.c.pm(5.0)
    m = Mod.Model(tt, observed, backg, 1.548048,-0.988016,0.338780, wavelength, spaceGroup, cell,
                (atomList, magAtomList), exclusions, magnetic=True,
                symmetry=symmetry, newSymmetry=basisSymmetry, base=6512, scale=59.08, sxtal=True, eta=0.0382, zero=0.08416, error=error)
    m.u.range(0,10)
    m.v.range(-10,0)
    m.w.range(0,10)
    m.scale.range(0,100)
    m.zero.pm(0.5)
    m.eta.range(0,1)
    m.base.pm(1000)
    for atomModel in m.atomListModel.atomModels:
        atomModel.x.range(0,1)
        atomModel.y.range(0,1)
        atomModel.z.range(0,1)
        atomModel.B.range(0,10)
        if atomModel.magnetic:
            for coeff in atomModel.coeffs:
                #coeff.range(-10, 10)
                coeff.range(-20,20)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = crystalCell
    S.diffPatternXtal(infoFile=infoFile, cell=cell, scale=2.34, tt=tt, 
                      obsIntensity=sfs2, wavelength=wavelength, basisSymmetry=basisSymmetry, 
                      magAtomList=magAtomList, plot=True, residuals=True, error=error, magnetic=False, 
                      info=True, base=0, refList=refList, extinctions=[2.974])
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
