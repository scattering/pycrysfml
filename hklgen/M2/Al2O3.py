import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
from pycrysfml import FloatVector, getSpaceGroup_crystalsys
np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"Al2O3 Background.BGR")
observedFile = os.path.join(DATAPATH,"Al2O3.dat")
infoFile = os.path.join(DATAPATH,"Al2O3.cif")

(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
#spaceGroup.xtalSystem = spaceGroup.xtalSystem.rstrip()
wavelength = 1.5403

backg = H.LinSpline(backgFile)
ttMin = 3
ttMax = 167.75
ttStep = 0.05
exclusions = None
tt, observed = H.readData(observedFile, kind="y", skiplines=1, start=ttMin,
                          stop=ttMax, step=ttStep)

def fit():
    # PYTHONPATH=. bumps Al2O3.py --fit=dream --store=M1 --burn=100 --steps=500
    cell = Mod.makeCell(crystalCell, getSpaceGroup_crystalsys(spaceGroup))
    cell.a.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                atoms, exclusions)
    m.u.range(0,2)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.scale.range(0,10)
    for atomModel in m.atomListModel.atomModels:
        atomModel.B.range(0, 10)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    m.atomListModel["Al1"].z.pm(0.1)
    m.atomListModel["O1"].x.pm(0.1)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = H.CrystalCell(FloatVector([4.761,4.761,13.000]),FloatVector([90,90,120]))
    uvw = [0.195228328354001,-0.164769183403005,0.0920158274607541]
    H.diffPattern(infoFile=infoFile, backgroundFile=backgFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=0.88765,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed))

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
