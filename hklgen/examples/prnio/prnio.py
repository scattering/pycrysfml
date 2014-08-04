import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import fswig_hklgen as H
import hkl_model as Mod
np.seterr(divide="ignore", invalid="ignore")    
def readInt(filename):
    lines = [line.strip() for line in open(filename)]
    s = []
    f = []
    for line in lines:
        split = line.split()
        if len(split) == 6:
            s.append(float(split[3]))
            f.append(float(split[4]))
    return s, f
DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"prnio.bac")
observedFile = os.path.join(DATAPATH,"prnio.int")
infoFile = os.path.join(DATAPATH,"prnio.pcr")
(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
wavelength = 0.8302
backg = H.LinSpline(None)
exclusions = None #[[0,10],[154,180]]
tt, observed = readInt(observedFile)
ttMin = min(tt)
ttMax = max(tt)
ttStep = (ttMax-ttMin)/len(tt)
#print tt , "\n", observed
def fit():
    # PYTHONPATH=. bumps Al2O3.py --fit=dream --store=M1 --burn=100 --steps=500
    cell = Mod.makeCell(crystalCell, spaceGroup.xtalSystem())
    cell.a.pm(0.5)
    cell.b.pm(0.5)
    cell.c.pm(0.5)
    m = Mod.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                atoms, exclusions, base=min(observed), zero=-0.09459, sxtal=True)
    m.u.range(0,2)
    m.zero.pm(0.1)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.eta.range(0,1)
    m.scale.range(0,10)
    m.base.pm(250)
    for atomModel in m.atomListModel.atomModels:
        atomModel.x.pm(0.1)
        atomModel.z.pm(0.1)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    #m.atomListModel["Al1"].z.pm(0.1)
    #m.atomListModel["O1"].x.pm(0.1)
    m.atomListModel["O3"].y.pm(0.1)
    m.atomListModel["Pb"].B.range(0,10)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    cell = H.CrystalCell([5.417799,5.414600,12.483399],[90,90,90])
    uvw = [0,0,1]
    H.diffPattern(infoFile=infoFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=1,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed), xtal=True)

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()
