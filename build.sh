#!/bin/bash
wd=$(pwd)
# Mac OS Support
if [[ "$OSTYPE" == "darwin"* ]]; then
CPPCOMP=g++-4.8
SEDCOM=gsed
LIBFLAGS='-lpython -lgfortran'
SOFLAGS='-shared -fPIC'
BIN_DIR='MacOS'
else
CPPCOMP=g++
SEDCOM=sed
LIBFLAGS='-lgfortran'
SOFLAGS='-shared -fPIC -rdynamic'
BIN_DIR='Linux'
fi
svn co http://forge.epn-campus.eu/svn/crysfml/Src
# inject python wrapper module
cp $wd/fort_methods/cfml_python/cfml_python.f90 $wd/Src/cfml_python.f90
# end injection
cd $wd
$wd/gen_list.py > $wd/list
cd $wd/Src/
mkdir wrap
$SEDCOM -i 's/.*/\L&/' *.f90
$wd/fix_deps.py
$wd/fix_type_decl.py
# wrap library
$wd/fortwrap.py --file-list=$wd/list -d wrap >& $wd/FortWrap_log
svn co http://forge.epn-campus.eu/svn/crysfml/Src
cp Src/*.f90 .
# inject python wrapper module again and fix line length
cp $wd/fort_methods/cfml_python/cfml_python.f90 $wd/Src/cfml_python.f90
$wd/fort_methods/fix_line_width.py $wd/Src/cfml_python.f90
# end injection
rm -r Src
# Build library with Makefile
cd $wd
make deps
$wd/fix_makefile_deps.py
make
make install
cd $wd/Src
$wd/add_cmplx.py
# auto-generate swig interface file
cd wrap
a="%module pycrysfml\n%{\n"
b="#include \""
for f in *.h
do
a=$a$b$f"\"\n"
done
a=$a"%}\n%include \"typemaps.i\"\n%include \"std_vector.i\"\nnamespace std {\n\t%template(FloatVector) vector<float>;\n}\n"
b="%include \""
for f in *.h
do
a=$a$b$f"\"\n"
done
echo -e $a > pycrysfml.i
swig -python -c++ pycrysfml.i
# compile Fortran wrapper
gfortran -fPIC -o CppWrappers.o -c CppWrappers.f90 ../crysfml.so -lstdc++ -Xlinker -rpath . -I..
#compile C++
$CPPCOMP -O2 -fPIC -c *.cpp *.cxx ../crysfml.so -I/usr/include/python2.7 $LIBFLAGS -Xlinker -rpath .
#build shared-object library
$CPPCOMP $SOFLAGS -o _pycrysfml.so -g -Wall ./*.o ../CFML_GlobalDeps_Linux.o ../CFML_Math_Gen.o ../CFML_LSQ_TypeDef.o ../CFML_Spher_Harm.o ../CFML_Random.o ../CFML_FFTs.o ../CFML_String_Util.o ../CFML_IO_Mess.o ../CFML_Profile_TOF.o ../CFML_Profile_Finger.o ../CFML_Profile_Functs.o ../CFML_Math_3D.o ../CFML_Optimization.o ../CFML_Optimization_LSQ.o ../CFML_Sym_Table.o ../CFML_Chem_Scatt.o ../CFML_Diffpatt.o ../CFML_Bonds_Table.o ../CFML_Cryst_Types.o ../CFML_Symmetry.o ../CFML_ILL_Instrm_Data.o ../CFML_Reflct_Util.o ../CFML_Atom_Mod.o ../CFML_Export_Vtk.o ../CFML_Sfac.o ../CFML_Geom_Calc.o ../CFML_Propagk.o ../CFML_Maps.o ../CFML_Molecules.o ../CFML_SXTAL_Geom.o ../CFML_Conf_Calc.o ../CFML_Form_CIF.o ../CFML_MagSymm.o ../CFML_Msfac.o ../CFML_Polar.o ../CFML_Refcodes.o ../cfml_python.o $LIBFLAGS
#$CPPCOMP -shared -fPIC -o _pycrysfml.so -g -Wall ./*.o ../crysfml.so $LIBFLAGS -L.. -I.. -Xlinker -rpath .
cd $wd
make clean
# install to hklgen directory
cp Src/wrap/pycrysfml.py hklgen/
cp Src/wrap/_pycrysfml.so hklgen/
rm hklgen/pycrysfml.pyc
# update binary release
cp Src/wrap/pycrysfml.py $wd/bin/$BIN_DIR/
cp Src/wrap/_pycrysfml.so $wd/bin/$BIN_DIR/
# update help file
cd $wd/bin/$BIN_DIR/
cp $wd/gen_help_file.py .
./gen_help_file.py
rm ./gen_help_file.py
