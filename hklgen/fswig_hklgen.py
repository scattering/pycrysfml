# fswig_hklgen.py - Summer 2014
#   Generates predicted diffraction peak positions for a given unit cell and space
#   group, using the Fortran CFML library.
#   Also uses the program "bumps" to perform fitting of calculated diffraction
#   patterns to observed data.
#   
#   FortWrap/Swig Version based on original from 2013
import numpy as np
import pylab
from pycrysfml import FortFuncs, atom_type, crystal_cell_type, space_group_type, atom_list_type, file_list_type
funcs = FortFuncs()
