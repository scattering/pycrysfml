# TiO2.py
# This file demonstrates the creation of a diffraction pattern based on
#   information not read from a file

import os,sys;sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import hklGen as H

np.seterr(divide="ignore", invalid="ignore")    

spaceGroup = H.SpaceGroup("136")
crystalCell = H.CrystalCell([4.5941,4.5941,2.9595], [90,90,90])
atom1 = H.Atom("Ti1", "Ti", [0,0,0], 2, 0.125, 0)
atom2 = H.Atom("O1", "O", [0.6952,0.6952,0], 4, 0.25, 0)
atoms = H.AtomList([atom1, atom2])
wavelength = 1.5403

ttMin = 10
ttMax = 100
ttStep = 0.05
exclusions = None

def main():
    H.diffPattern(spaceGroup=spaceGroup, cell=crystalCell, atomList=atoms,
                  wavelength=wavelength, ttMin=ttMin, ttMax=ttMax, ttStep=0.05,
                  info=True, plot=True)

if __name__ == "__main__":
    # program run normally
    main()
    