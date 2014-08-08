import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
np.seterr(divide="ignore", invalid="ignore")    

DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"lamn_3t2.bac")
observedFile = os.path.join(DATAPATH,"lamn_3t2.dat")
infoFile = os.path.join(DATAPATH,"lamn_3t2.pcr")

(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
wavelength = 1.225300
backg = H.LinSpline(None)
exclusions = [[-1.00,11.00],[125.0,180.0]]
tt, observed = H.readIllData(observedFile, "D1A", backgFile)
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
def fit():
    # PYTHONPATH=. bumps Al2O3.py --fit=dream --store=M1 --burn=100 --steps=500
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem())
    cell.a.pm(0.000001)
    cell.b.pm(0.000001)
    cell.c.pm(0.000001)
    cell.a.set([5.536482,5.747128,7.692561][0])
    cell.b.set([5.536482,5.747128,7.692561][1])
    cell.c.set([5.536482,5.747128,7.692561][2])
    m = Mod.Model(tt, observed, backg, 0.176001,-0.197806,0.091452, wavelength, spaceGroup, cell,
                atoms, exclusions, base=min(observed), zero=0.00029)
    m.u.pm(0.00001)
    m.v.pm(0.00001)
    m.w.pm(0.00001)
    m.u.set([0.176001,-0.197806,0.091452][0])
    m.v.set([0.176001,-0.197806,0.091452][1])
    m.w.set([0.176001,-0.197806,0.091452][2])
    #m.u.range(0,2)
    m.zero.pm(0.1)
    #m.v.range(-2,0)
    #m.w.range(0,2)
    m.eta.range(0,1)
    m.scale.range(0,10)
    m.base.pm(250)
    for atomModel in m.atomListModel.atomModels:
        #atomModel.B.range(0,10)
        if H.getAtom_chemsymb(atomModel.atom).lower() != "mn":
            atomModel.x.range(0,1)
            atomModel.y.range(0,1)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    m.atomListModel["O2"].z.range(0,1)
    #m.atomListModel["O1"].x.pm(0.1)
    #m.atomListModel["O3"].y.pm(0.1)
    #m.atomListModel["Pb"].B.range(0,10)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = H.CrystalCell([5.536482,5.747128,7.692561],[90,90,90])
    uvw = [0.176001,-0.197806,0.091452]
    H.diffPattern(infoFile=infoFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=0.86927E-01,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed), base=min(observed))

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
