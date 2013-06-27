pycrysfml
=========

This is a library to add Python bindings to the crysfml crystallographic library

The binding uses ctypes, which requires detailed knowledge of the representation
of fortran objects in memory and the details of the fortran subroutine calling
interface.  Only gfortran is supported in this release.  Details of other compilers
can be found from the code in the [chasm-interop](http://chasm-interop.cvs.sourceforge.net/viewvc/chasm-interop/chasm/include/compilers/)
C++/F90 interoperation library from LANL.


build
=====

To build the library you need a copy of the crysfml source.  The fortran source for
crysfml should be placed in the Src subdirectory.  To generate the crysfml library
dependencies (gfortran must compile the fortran modules in the correct order), use:

	make deps crysfml

This creates Makefile.deps, which contains the dependencies, then builds the crysfml
library.  If any new modules or use statements are added to the crysfml source, then 
"make deps" will need to be rerun.

Windows builds will require tweaks to Makefile, such as selecting the correct set
of libraries to build for makedeps and renaming libcrysfml.so to crysfml.dll.
