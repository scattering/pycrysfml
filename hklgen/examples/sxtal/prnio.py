import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
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
observedFile = os.path.join(DATAPATH,r"prnio.int")
infoFile = os.path.join(DATAPATH,r"prnio.cfl")
spaceGroup, crystalCell, atomList = H.readInfo(infoFile)
exclusions = []
# return wavelength, refList, sfs2, error, two-theta, and four-circle parameters
wavelength, refList, sfs2, error = S.readIntFile(observedFile, kind="int", cell=crystalCell)
tt = [H.twoTheta(H.calcS(crystalCell, ref.hkl), wavelength) for ref in refList]
backg = None

def fit():
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem)
    #cell.a.pm(5.0)
    #cell.c.pm(5.0)
    m = S.Model(tt, sfs2, backg, wavelength, spaceGroup, cell,
                [atomList], exclusions,
                 scale=0.06298, error=error, hkls=refList, extinction=[0.0001054])
    m.scale.range(0,100)
    #m.base.pm(1000)
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
    S.diffPatternXtal(infoFile=infoFile, cell=cell, scale=5.5*0.06298, tt=tt, 
                      obsIntensity=sfs2, wavelength=wavelength, plot=True, residuals=True, error=None, magnetic=False, 
                      info=True, base=0, refList=refList, extinctions=[0.0001054])
if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
