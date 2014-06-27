#!/bin/bash
wd=$(pwd)
svn co http://forge.epn-campus.eu/svn/crysfml/Src
# inject fortran methods
cd $wd/fort_methods
./cat.sh
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Atom_Mod.f90 $wd/fort_methods/cfml_atom_mod_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Cryst_Types.f90 $wd/fort_methods/cfml_cryst_types_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Form_CIF.f90 $wd/fort_methods/cfml_form_cif_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_MagSymm.f90 $wd/fort_methods/cfml_magsymm_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Msfac.f90 $wd/fort_methods/cfml_msfac_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Reflct_Util.f90 $wd/fort_methods/cfml_reflct_util_addns.f90
$wd/fort_methods/inject_methods.py $wd/crysfml/Src/CFML_Symmetry.f90 $wd/fort_methods/cfml_symmetry_addns.f90
./clean.sh
# end method injection
cd $wd
$wd/gen_list.py > $wd/list
cd $wd/Src/
mkdir wrap
sed -i 's/.*/\L&/' *.f90
$wd/fix_deps.py
$wd/fix_type_decl.py
# wrap library
$wd/fortwrap.py --file-list=$wd/list -d wrap >& $wd/FortWrap_log
svn co http://forge.epn-campus.eu/svn/crysfml/Src
cp Src/*.f90 .
rm -r Src
# inject fortran methods again and fix line length
cd $wd/fort_methods
./cat.sh
# fix line width for compiling
for f in *addns.f90
do
./fix_line_width.py $f
done
# inject methods
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Atom_Mod.f90 $wd/fort_methods/cfml_atom_mod_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Cryst_Types.f90 $wd/fort_methods/cfml_cryst_types_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Form_CIF.f90 $wd/fort_methods/cfml_form_cif_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_MagSymm.f90 $wd/fort_methods/cfml_magsymm_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Msfac.f90 $wd/fort_methods/cfml_msfac_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Reflct_Util.f90 $wd/fort_methods/cfml_reflct_util_addns.f90
$wd/fort_methods/inject_methods.py $wd/Src/CFML_Symmetry.f90 $wd/fort_methods/cfml_symmetry_addns.f90
./clean.sh
cd $wd
# end method injection
#$wd/build_so.sh
#update for building with makefile
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
a=$a"%}\n%include \"std_vector.i\"\nnamespace std {\n\t%template(FloatVector) vector<float>;\n}\n"
b="%include \""
for f in *.h
do
a=$a$b$f"\"\n"
done
echo -e $a > pycrysfml.i
swig -python -c++ pycrysfml.i
# compile Fortran wrapper
gfortran -fPIC -o CppWrappers.o -c CppWrappers.f90 ../crysfml.so -lstdc++ -Xlinker -rpath . -I..
# fix circular dependencies
cp $wd/fix_circular/* .
#compile c++
g++ -O2 -fPIC -c *.cpp *.cxx ../crysfml.so -I/usr/include/python2.7 -lgfortran -Xlinker -rpath .
#build shared-object library
#g++ -shared -fPIC -rdynamic -o _pycrysfml.so -g -Wall ./*.o ../CFML_GlobalDeps_Linux.o ../CFML_Math_Gen.o ../CFML_LSQ_TypeDef.o ../CFML_Spher_Harm.o ../CFML_Random.o ../CFML_FFTs.o ../CFML_String_Util.o ../CFML_IO_Mess.o ../CFML_Profile_TOF.o ../CFML_Profile_Finger.o ../CFML_Profile_Functs.o ../CFML_Math_3D.o ../CFML_Optimization.o ../CFML_Optimization_LSQ.o ../CFML_Sym_Table.o ../CFML_Chem_Scatt.o ../CFML_Diffpatt.o ../CFML_Bonds_Table.o ../CFML_Cryst_Types.o ../CFML_Symmetry.o ../CFML_ILL_Instrm_Data.o ../CFML_EoS_Mod.o ../CFML_Reflct_Util.o ../CFML_Atom_Mod.o ../CFML_Export_Vtk.o ../CFML_Sfac.o ../CFML_Geom_Calc.o ../CFML_Propagk.o ../CFML_Maps.o ../CFML_Molecules.o ../CFML_SXTAL_Geom.o ../CFML_Conf_Calc.o ../CFML_Form_CIF.o ../CFML_Optimization_SAn.o ../CFML_MagSymm.o ../CFML_Msfac.o ../CFML_Polar.o ../CFML_Refcodes.o -lgfortran -I.
g++ -shared -fPIC -rdynamic -o _pycrysfml.so -g -Wall ./*.o ../CFML_GlobalDeps_Linux.o ../CFML_Math_Gen.o ../CFML_LSQ_TypeDef.o ../CFML_Spher_Harm.o ../CFML_Random.o ../CFML_FFTs.o ../CFML_String_Util.o ../CFML_IO_Mess.o ../CFML_Profile_TOF.o ../CFML_Profile_Finger.o ../CFML_Profile_Functs.o ../CFML_Math_3D.o ../CFML_Optimization.o ../CFML_Optimization_LSQ.o ../CFML_Sym_Table.o ../CFML_Chem_Scatt.o ../CFML_Diffpatt.o ../CFML_Bonds_Table.o ../CFML_Cryst_Types.o ../CFML_Symmetry.o ../CFML_ILL_Instrm_Data.o ../CFML_EoS_Mod.o ../CFML_Reflct_Util.o ../CFML_Atom_Mod.o ../CFML_Export_Vtk.o ../CFML_Sfac.o ../CFML_Geom_Calc.o ../CFML_Propagk.o ../CFML_Maps.o ../CFML_Molecules.o ../CFML_SXTAL_Geom.o ../CFML_Conf_Calc.o ../CFML_Form_CIF.o ../CFML_MagSymm.o ../CFML_Msfac.o ../CFML_Polar.o ../CFML_Refcodes.o -lgfortran
cd $wd
make clean
