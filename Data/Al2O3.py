import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import os

import numpy as np
import hklGen as H

# prevent numpy warnings from popping up
np.seterr(divide="ignore", invalid="ignore")    

# files should be in the same directory
DATAPATH = os.path.dirname(os.path.abspath(__file__))
backgFile = os.path.join(DATAPATH,"Al2O3 Background.BGR")
observedFile = os.path.join(DATAPATH,"Al2O3.dat")
infoFile = os.path.join(DATAPATH,"Al2O3.cif")

# read in information
(spaceGroup, crystalCell, atoms) = H.readInfo(infoFile)
wavelength = 1.5403

backg = H.LinSpline(backgFile)
ttMin = 3
ttMax = 167.75
ttStep = 0.05
exclusions = None
tt, observed = H.readData(observedFile, kind="y", skiplines=1, start=ttMin,
                          stop=ttMax, step=ttStep)

def fit():
    # construct a cell object and set ranges on parameters
    cell = H.makeCell(crystalCell, spaceGroup.xtalSystem)
    cell.a.pm(0.5)
    cell.c.pm(0.5)
    m = H.Model(tt, observed, backg, 0, 0, 1, wavelength, spaceGroup, cell,
                atoms, exclusions)
    m.u.range(0,2)
    m.v.range(-2,0)
    m.w.range(0,2)
    m.scale.range(0,10)
    for atomModel in m.atomListModel.atomModels:
        atomModel.B.range(0, 10)
        if (atomModel.atom.multip == atomModel.sgmultip):
            # atom lies on a general position; refine its location
            # (not required for Al2O3)
            atomModel.x.pm(0.1)
            atomModel.y.pm(0.1)
            atomModel.z.pm(0.1)
    # take care of special positions
    m.atomListModel["Al1"].z.pm(0.1)
    m.atomListModel["O1"].x.pm(0.1)
    M = bumps.FitProblem(m)
    M.model_update()
    return M

def main():
    # change parameters to match the results of a previous refinement
    #   (should show very good agreement with observed data)
    # do not pass these to diffPattern() to get diffraction pattern pre-refinement
    cell = H.CrystalCell([4.761,4.761,13.000],[90,90,120])
    uvw = [0.17213,-0.11884,0.07721]
    scale = 1.47113
    H.diffPattern(infoFile=infoFile, backgroundFile=backgFile, wavelength=wavelength,
                  cell=cell, uvw=uvw, scale=scale,
                  ttMin=ttMin, ttMax=ttMax, info=True, plot=True,
                  observedData=(tt,observed))

if __name__ == "__main__":
    # program run normally
    main()
else:
    # called using bumps
    import bumps.names as bumps
    problem = fit()