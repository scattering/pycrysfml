# fswig_hklgen.py - Summer 2014
#   Generates predicted diffraction peak positions for a given unit cell and space
#   group, using the Fortran CFML library.
#   Also uses the program "bumps" to perform fitting of calculated diffraction
#   patterns to observed data.
#   
#   FortWrap/Swig Version based on original from 2013
#

from pycrysfml import *
funcs = FortFuncs()
def readInfo(filename):
    # read the file
    