import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"ncaf_3t2.bac")
observedFile = os.path.join(DATAPATH,"ncaf_3t2.dat")
infoFile = os.path.join(DATAPATH,"ncaf_3t2.pcr")

(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
wavelength = 1.225300
backg = H.LinSpline(None)
ttMin = 0.0
ttMax = 125.450004578
ttStep = 0.049980081505
exclusions = [[0,5],[125,180]]
tt, observed = H.readIllData(observedFile, "D1A", backgFile)
def fit():
    # PYTHONPATH=. bumps Al2O3.py --fit=dream --store=M1 --burn=100 --steps=500
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem())
    cell.a.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                atoms, exclusions, base=min(observed), zero=-0.09459)
    m.u.range(0,2)
    m.zero.pm(0.1)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.eta.range(0,1)
    m.scale.range(0,10)
    m.base.pm(25)
    for atomModel in m.atomListModel.atomModels:
        atomModel.x.pm(1.0)
        if H.getAtom_chemsymb(atomModel.atom).lower() != "ca":
            atomModel.y.pm(1.0)
            atomModel.z.pm(1.0)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    #m.atomListModel["Al1"].z.pm(0.1)
    #m.atomListModel["O1"].x.pm(0.1)
    #m.atomListModel["O3"].y.pm(0.1)
    #m.atomListModel["Pb"].B.range(0,10)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = H.CrystalCell([10.250222,10.250222,10.250222],[90,90,90])
    uvw = [0.270721,-0.425009,0.216537]
    H.diffPattern(infoFile=infoFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=0.34815E-02,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed), base=min(observed))

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
