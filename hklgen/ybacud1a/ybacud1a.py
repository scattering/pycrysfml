import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"ybacud1a.bac")
observedFile = os.path.join(DATAPATH,"ybacud1a.dat")
infoFile = os.path.join(DATAPATH,"ybacud1a.pcr")

(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
wavelength = 1.907972
backg = H.LinSpline(None)
ttMin = 0.0
ttMax = 155.449996948
ttStep = 0.0499839218483
exclusions = [[0,6],[155.45,180]]
tt, observed = H.readIllData(observedFile, "D1A", backgFile)
def fit():
    # PYTHONPATH=. bumps Al2O3.py --fit=dream --store=M1 --burn=100 --steps=500
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem())
    cell.a.pm(0.5)
    cell.b.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                atoms, exclusions, base=min(observed), zero=0.19089)
    m.u.range(0,2)
    m.zero.pm(0.1)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.eta.range(0,1)
    m.scale.range(0,10)
    m.base.pm(250)
    for atomModel in m.atomListModel.atomModels:
        atomModel.B.range(0,10)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    m.atomListModel["Ba"].z.pm(1.0)
    m.atomListModel["CaB"].z.pm(1.0)
    m.atomListModel["Cu2"].z.pm(1.0)
    m.atomListModel["O1"].z.pm(1.0)
    m.atomListModel["O2"].z.pm(1.0)
    m.atomListModel["O3"].z.pm(1.0)
    #m.atomListModel["O1"].x.pm(0.1)
    #m.atomListModel["O3"].y.pm(0.1)
    #m.atomListModel["Pb"].B.range(0,10)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = H.CrystalCell([3.815149, 3.880053, 11.664815],[90,90,90])
    uvw = [0.198607,-0.322513,0.338442]
    H.diffPattern(infoFile=infoFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=2.5189,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed))

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
