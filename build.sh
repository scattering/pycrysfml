#!/bin/bash
# PyCrysFML build script
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
# Wraps CrysFML Library into C++ using modified FortWrap and
# auto-generates swig interface for wrapping into Python
# Compiles everything and installs pycrysfml binaries in
# bin/<Platform> and hklgen/
# Updated Jan 2017 to account for changes in CrysFML
# Upstream requires gfortran 5 or higher on all platforms due to use of IEEE library
# functions from F2003, but we patch it so that it compiles on gfortran 4.8.
set -x
wd=$(pwd)
tools=$wd/build_tools
# Mac OS Support
if [[ "$OSTYPE" == "darwin"* ]]; then
	CPPCOMP=g++-6
	SEDCOM=gsed
	LIBFLAGS='-lpython -lgfortran'
	SOFLAGS='-shared -fPIC'
	BIN_DIR='MacOS'
	PY_HEADERS='/usr/include/python3.6'
	STR_MOD='gf'
	FORTCOMP=gfortran
	SOLIB_EXT='.so'
fi
if [[ "$OSTYPE" == "linux"* ]]; then
	CPPCOMP=g++
	SEDCOM=sed
        ## Simplest option is to hardcode the python link library.
	#LIBFLAGS="-lpython3.6m -lgfortran"
        ## The -fno-lto flag was needed for the ubuntu compile. 
	#LIBFLAGS="-lgfortran -fno-lto `python3-config --ldflags`"
        ## Remove $CONDA_PREFIX/lib; it probably has a conflicting libgfortran.
        ## The symptom of this was the error:
        ##
        ##    At line 1523 of file Src/CFML_String_Util.f90
        ##    Internal Error: get_unit(): Bad internal unit KIND
        ##
	LIBFLAGS="-lgfortran -fno-lto `python3-config --ldflags | sed -e's,-L[^ ]*/lib ,,'g`"
	SOFLAGS='-shared -fPIC -rdynamic'
	BIN_DIR='Linux'
	PY_HEADERS=`python3-config --includes`
	STR_MOD='gf'
	FORTCOMP=gfortran
	SOLIB_EXT='.so'
fi
if [[ "$OSTYPE" == "msys"* ]]; then
	# running on windows
	# build requires full install of mingw and Python XY in their default locations
	# additionally command line based svn and git clients answering to standard calls are required
	# some windows svn clients do not hold the processing thread while they are checking out libraries
	# so it may be necessary to manually checkout the source first and then build from a local copy
	# In order for sed to work the user account running the build must be given "Full Control" over the
	# pycrysfml directory in the windows folder permissions window
	CPPCOMP=g++
	SEDCOM=sed
	LIBFLAGS='-L/c/Python27/libs -lgfortran -lpython27'
	SOFLAGS='-shared -fPIC'
	BIN_DIR='Windows'
	PY_HEADERS='-I/c/Python27/include'
	STR_MOD='gf'
	FORTCOMP=gfortran
	SOLIB_EXT='.pyd'
	ln -s /c/Python27/python /usr/bin/python
	mount c:/mingw /mingw
fi
if [[ "$HOSTNAME" == "darter"* ]]; then
	#running on Cray XC30
	module load python/2.7.6
	module swap PrgEnv-cray PrgEnv-gnu
	CPPCOMP=CC
	BIN_DIR='Cray_XC30'
	LIBFLAGS='-lgfortran -L/sw/xc30_cle5.2_pe2014-09/python/2.7.6/cle5.2_gnu4.9.1/lib/ -lpython3.6'
	STR_MOD='LF'
	PY_HEADERS='-I/sw/xc30_cle5.2_pe2014-09/python/2.7.6/cle5.2_gnu4.9.1/include/python3.6'
	FORTCOMP=ftn
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/xc30_cle5.2_pe2014-09/python/2.7.6/cle5.2_gnu4.9.1/lib
fi
if [[ "$HOSTNAME" == "rocks"* ]]; then
	#running on rocks cluster
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/python/lib
	LIBFLAGS='-lgfortran -L/opt/python/lib -lpython3.6'
	PY_HEADERS='-I/opt/python/include/python3.6'
	alias python=/opt/python/bin/python3.6
	STR_MOD='LF'
fi
if [ $# -lt 1 ]; then
	svn co --non-interactive --trust-server-cert https://forge.epn-campus.eu/svn/crysfml/Src
	# svn co http://forge.epn-campus.eu/svn/crysfml/SrcCpp
	# remove stale versions of String_Utilities Module
	rm Src/*_gf.f90 Src/*_LF.f90
	# Add patch to CrysFML to fix MsFac bugs
	cd Src
	patch -f < ../CFML_Msfac.patch
        # Use the following with gfortran<5.0
	patch -f < ../CFML_String_Util.patch
	# Stripping out the calls to the nexus wrapper library for ILL data becaause we do not need them.
	# The original code refers to the SrcCpp package in the crysfml directory, which depends on the
	# nexus libraries and a "bosc.h" header that isn't included in the SrcCpp, so we aren't trying
	# to build it, and instead ignoreng ILL data files.
	$SEDCOM -i 's/^ \( *call\)/!\1/' CFML_ILL_Instrm_Data_Nexus.f90
	cd $wd
	# End Patch to CrysFML
else
	cp -r $1 $wd/Src
fi
# inject python wrapper module
cp $wd/fort_methods/cfml_python/cfml_python.f90 $wd/Src/cfml_python.f90
# end injection
cd $wd
$tools/gen_list.py > $wd/list
cd $wd/Src/
mkdir $wd/Src/wrap
# Lowercase everything in the source directory
$SEDCOM -i 's/.*/\L&/' *.f90
$tools/fix_deps.py
$tools/fix_type_decl.py
# wrap library
$tools/fortwrap.py --file-list=$wd/list -d $wd/Src/wrap >& $wd/FortWrap_log || (echo "fortwrap failed" && cat $wd/FortWrap_log && exit 1)
if [ $# -lt 1 ]; then
	svn co --non-interactive --trust-server-cert https://forge.epn-campus.eu/svn/crysfml/Src
	# remove stale versions of String_Utilities Module
	rm Src/*_gf.f90 Src/*_LF.f90
	# Add patch to CrysFML to fix MsFac bugs
	cp *.patch Src/
	cd Src
	patch -f < CFML_Msfac.patch
        # Use the following with gfortran<5.0
	patch -f < CFML_String_Util.patch
	# Stripping out the calls to the nexus wrapper library for ILL data becaause we do not need them.
	# The original code refers to the SrcCpp package in the crysfml directory, which depends on the
	# nexus libraries and a "bosc.h" header that isn't included in the SrcCpp, so we aren't trying
	# to build it, and instead ignoreng ILL data files.
	mv CFML_ILL_Instrm_Data_Nexus.f90 CFML_ILL_Instrm_Data_Nexus.f90-orig
	sed -e's/^ \( *call\)/!\1/' CFML_ILL_Instrm_Data_Nexus.f90-orig > CFML_ILL_Instrm_Data_Nexus.f90
	cd $wd/Src
	# End Patch to CrysFML
	cp Src/*.f90 .
	rm -rf Src
else
	cp $1/*.f90 .
fi
# inject python wrapper module again and fix line length
cp $wd/fort_methods/cfml_python/cfml_python.f90 $wd/Src/cfml_python.f90
$wd/fort_methods/fix_line_width.py $wd/Src/cfml_python.f90
# end injection
# Build library with Makefile
cd $wd
make deps
$tools/fix_makefile_deps.py
make
make install
cd $wd/Src
$tools/add_cmplx.py
# copy C++ wrapper files
cp $wd/cpp_modules/* $wd/Src/wrap/
# auto-generate swig interface file
cd wrap
a="%module pycrysfml\n%{\n"
b="#include \""
set +x
for f in *.h; do
	a=$a$b$f"\"\n"
done
a=$a"%}\n%include \"cpointer.i\"\n%pointer_class(int, intp);\n%pointer_class(double, doublep);\n%pointer_class(float, floatp);\n%include \"std_string.i\"\n%include \"cstring.i\"\n%include \"std_vector.i\"\nnamespace std {\n\t%template(FloatVector) vector<float>;\n\t%template(FloatMatrix) vector< vector<float> >;\n\t%template(IntVector) vector<int>;\n}\n"
b="%include \""
for f in *.h; do
	a=$a$b$f"\"\n"
done
set -x
echo -e $a > pycrysfml.i
swig -python -c++ pycrysfml.i
# compile Fortran wrapper
echo "compiling wrapper"
$FORTCOMP -fPIC -o CppWrappers.o -c CppWrappers.f90 ../crysfml.so -lstdc++ -Xlinker -rpath . -I..
#compile C++
$CPPCOMP -O2 -fPIC -c *.cpp *.cxx $PY_HEADERS  ../crysfml.so $LIBFLAGS -Xlinker -rpath .
#build shared-object library
echo "making shared-object library"
##$CPPCOMP $SOFLAGS -o _pycrysfml.so -g -Wall ./*.o ../CFML_GlobalDeps_Linux.o ../CFML_Math_Gen.o ../CFML_LSQ_TypeDef.o ../CFML_Spher_Harm.o ../CFML_Random.o ../CFML_FFTs.o ../CFML_String_Util.o ../CFML_IO_Mess.o ../CFML_Profile_TOF.o ../CFML_Profile_Finger.o ../CFML_Profile_Functs.o ../CFML_Math_3D.o ../CFML_Optimization.o ../CFML_Optimization_LSQ.o ../CFML_Sym_Table.o ../CFML_Chem_Scatt.o ../CFML_Diffpatt.o ../CFML_Bonds_Table.o ../CFML_Cryst_Types.o ../CFML_Symmetry.o ../CFML_ILL_Instrm_Data.o ../CFML_Reflct_Util.o ../CFML_Atom_Mod.o ../CFML_Export_Vtk.o ../CFML_Sfac.o ../CFML_Geom_Calc.o ../CFML_Propagk.o ../CFML_Maps.o ../CFML_Molecules.o ../CFML_SXTAL_Geom.o ../CFML_Conf_Calc.o ../CFML_Form_CIF.o ../CFML_MagSymm.o ../CFML_Msfac.o ../CFML_Polar.o ../CFML_Refcodes.o ../cfml_python.o $LIBFLAGS

#$CPPCOMP $SOFLAGS -o _pycrysfml$SOLIB_EXT -g -Wall ./*.o ../CFML_GlobalDeps_Linux.o ../CFML_Math_Gen.o ../CFML_LSQ_TypeDef.o ../CFML_Spher_Harm.o ../CFML_Random.o ../CFML_FFTs.o ../CFML_String_Util.o ../CFML_IO_Mess.o ../CFML_Profile_TOF.o ../CFML_Profile_Finger.o ../CFML_Profile_Functs.o ../CFML_Math_3D.o ../CFML_Optimization.o ../CFML_Optimization_LSQ.o ../CFML_Sym_Table.o ../CFML_Chem_Scatt.o ../CFML_Diffpatt.o ../CFML_Bonds_Table.o ../CFML_Cryst_Types.o ../CFML_Symmetry.o ../CFML_ILL_Instrm_Data.o ../CFML_ILL_Instrm_Data_Nexus.o ../CFML_Reflct_Util.o ../CFML_Atom_Mod.o ../CFML_Export_Vtk_LF95.o ../CFML_Sfac.o ../CFML_Geom_Calc.o ../CFML_Propagk.o ../CFML_Maps.o ../CFML_Molecules.o ../CFML_SXTAL_Geom.o ../CFML_Conf_Calc.o ../CFML_Form_CIF.o ../CFML_MagSymm.o ../CFML_Msfac.o ../CFML_Polar.o ../CFML_Refcodes.o ../CFML_BVSpar.o ../CFML_Extinction_Correction.o ../CFML_Magnetic_Groups.o ../f2kcli.o ../cfml_python.o $LIBFLAGS

$CPPCOMP $SOFLAGS -o _pycrysfml$SOLIB_EXT -g -Wall ./*.o  \
	../CFML_Atom_Mod.o \
	../CFML_BVSpar.o \
	../CFML_Bonds_Table.o \
	../CFML_Chem_Scatt.o \
	../CFML_Conf_Calc.o \
	../CFML_Cryst_Types.o \
	../CFML_Diffpatt.o \
	../CFML_EisPack.o \
	../CFML_Export_Vtk_LF95.o \
	../CFML_Extinction_Correction.o \
	../CFML_FFTs.o \
	../CFML_Form_CIF.o \
	../CFML_Geom_Calc.o \
	../CFML_GlobalDeps_Linux.o \
	../CFML_ILL_Instrm_Data.o \
	../CFML_ILL_Instrm_Data_Nexus.o \
	../CFML_IO_Mess.o \
	../CFML_LSQ_TypeDef.o \
	../CFML_MagSymm.o \
	../CFML_Magnetic_Groups.o \
	../CFML_Maps.o \
	../CFML_Math_3D.o \
	../CFML_Math_Gen.o \
	../CFML_Molecules.o \
	../CFML_Msfac.o \
	../CFML_Optimization.o \
	../CFML_Optimization_LSQ.o \
	../CFML_Percolation.o \
	../CFML_Polar.o \
	../CFML_Profile_Finger.o \
	../CFML_Profile_Functs.o \
	../CFML_Profile_TOF.o \
	../CFML_Propagk.o \
	../CFML_Random.o \
	../CFML_Rational_Arithmetic.o \
	../CFML_Refcodes.o \
	../CFML_Reflct_Util.o \
	../CFML_SXTAL_Geom.o \
	../CFML_Sfac.o \
	../CFML_Spher_Harm.o \
	../CFML_String_Util.o \
	../CFML_Sym_Table.o \
	../CFML_Symmetry.o \
	../f2kcli.o \
        ../cfml_python.o $LIBFLAGS

	# ../CFML_Optimization_SAn.o 
	#../CFML_IO_MessWin.o 
	#../CFML_IO_MessRW.o 
	#../CFML_EoS_Mod.o 
	#../CFML_Export_Vtk.o 
	#../CFML_FlipR_Mod.o 
	#../CFML_GlobalDeps_Linux_Intel.o 
	#../CFML_GlobalDeps_MacOS.o 
	#../CFML_GlobalDeps_MacOS_Intel.o 
	#../CFML_GlobalDeps_Windows.o 
	#../CFML_GlobalDeps_Windows_intel.o 
	#../CFML_GlobalDeps_Windows_intel64.o 

#$CPPCOMP -shared -fPIC -o _pycrysfml.so -g -Wall ./*.o ../crysfml.so $LIBFLAGS -L.. -I.. -Xlinker -rpath .
cd $wd
#make clean
# install to hklgen directory
cp Src/wrap/pycrysfml.py hklgen/
cp Src/wrap/_pycrysfml$SOLIB_EXT hklgen/
#rm hklgen/pycrysfml.pyc
# update binary release
#cp Src/wrap/pycrysfml.py $wd/bin/$BIN_DIR/
#cp Src/wrap/_pycrysfml$SOLIB_EXT $wd/bin/$BIN_DIR/
# update help file
#cd $wd/bin/$BIN_DIR/
#cp $wd/gen_help_file.py .
#./gen_help_file.py
#rm ./gen_help_file.py
#cd $wd
#$wd/clean.sh
