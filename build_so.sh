#!/bin/bash
#compile fortran code
COMP=gfortran
OPT1="-fPIC -g -c -Wall -rdynamic -lstdc++ -o"
$COMP CFML_GlobalDeps_Linux.f90 $OPT1 CFML_GlobalDeps_Linux.o
$COMP CFML_Math_Gen.f90 $OPT1 CFML_Math_Gen.o
$COMP CFML_LSQ_TypeDef.f90 $OPT1 CFML_LSQ_TypeDef.o
$COMP CFML_Spher_Harm.f90 $OPT1 CFML_Spher_Harm.o
$COMP CFML_Random.f90 $OPT1 CFML_Random.o
$COMP CFML_FFTs.f90 $OPT1 CFML_FFTs.o
$COMP CFML_String_Util.f90 $OPT1 CFML_String_Util.o
$COMP CFML_IO_Mess.f90 $OPT1 CFML_IO_Mess.o
$COMP CFML_Profile_TOF.f90 $OPT1 CFML_Profile_TOF.o
$COMP CFML_Profile_Finger.f90 $OPT1 CFML_Profile_Finger.o
$COMP CFML_Profile_Functs.f90 $OPT1 CFML_Profile_Functs.o
$COMP CFML_Math_3D.f90 $OPT1 CFML_Math_3D.o
$COMP CFML_Optimization.f90 $OPT1 CFML_Optimization.o
$COMP CFML_LSQ_TypeDef.f90 $OPT1 CFML_LSQ_TypeDef.o
$COMP CFML_Optimization_LSQ.f90 $OPT1 CFML_Optimization_LSQ.o
$COMP CFML_Sym_Table.f90 $OPT1 CFML_Sym_Table.o
$COMP CFML_Chem_Scatt.f90 $OPT1 CFML_Chem_Scatt.o
$COMP CFML_Diffpatt.f90 $OPT1 CFML_Diffpatt.o
$COMP CFML_Bonds_Table.f90 $OPT1 CFML_Bonds_Table.o
$COMP CFML_Cryst_Types.f90 $OPT1 CFML_Cryst_Types.o
$COMP CFML_Symmetry.f90 $OPT1 CFML_Symmetry.o
$COMP CFML_ILL_Instrm_Data.f90 $OPT1 CFML_ILL_Instrm_Data.o
$COMP CFML_EoS_Mod.f90 $OPT1 CFML_EoS_Mod.o
$COMP CFML_Reflct_Util.f90 $OPT1 CFML_Reflct_Util.o
$COMP CFML_Atom_Mod.f90 $OPT1 CFML_Atom_Mod.o
$COMP CFML_Export_Vtk.f90 $OPT1 CFML_Export_Vtk.o
$COMP CFML_Sfac.f90 $OPT1 CFML_Sfac.o
$COMP CFML_Geom_Calc.f90 $OPT1 CFML_Geom_Calc.o
$COMP CFML_Propagk.f90 $OPT1 CFML_Propagk.o
$COMP CFML_Maps.f90 $OPT1 CFML_Maps.o
$COMP CFML_Molecules.f90 $OPT1 CFML_Molecules.o
$COMP CFML_SXTAL_Geom.f90 $OPT1 CFML_SXTAL_Geom.o
$COMP CFML_Conf_Calc.f90 $OPT1 CFML_Conf_Calc.o
$COMP CFML_Form_CIF.f90 $OPT1 CFML_Form_CIF.o
#$COMP CFML_Optimization_SAn.f90 $OPT1 CFML_Optimization_SAn.o
$COMP CFML_MagSymm.f90 $OPT1 CFML_MagSymm.o
$COMP CFML_Msfac.f90 $OPT1 CFML_Msfac.o
$COMP CFML_Polar.f90 $OPT1 CFML_Polar.o
$COMP CFML_Refcodes.f90 $OPT1 CFML_Refcodes.o
#$COMP CFML_ILL_Instrm_Data_LF.f90 $OPT1 CFML_ILL_Instrm_Data_LF.o
#$COMP CFML_IO_MessRW.f90 $OPT1 CFML_IO_MessRW.o
#$COMP CFML_IO_MessWin.f90 $OPT1 CFML_IO_MessWin.o
$COMP f2kcli.f90 $OPT1 f2kcli.o


gfortran -shared -Wl,-soname,crysfml.so -o crysfml.so *.o -lc
echo "Library Compiled"
# create wrapper directory
#if [ ! -d "./wrap" ]; then
#	mkdir wrap
#fi
# run fortwrap
#for f in /home/jel/TestEnv/crysfml/Src/CFML_GlobalDeps_Linux.f90 /home/jel/TestEnv/crysfml/Src/CFML_Math_Gen.f90 /home/jel/TestEnv/crysfml/Src/CFML_LSQ_TypeDef.f90 /home/jel/TestEnv/crysfml/Src/CFML_Spher_Harm.f90 /home/jel/TestEnv/crysfml/Src/CFML_Random.f90 /home/jel/TestEnv/crysfml/Src/CFML_FFTs.f90 /home/jel/TestEnv/crysfml/Src/CFML_String_Util.f90 /home/jel/TestEnv/crysfml/Src/CFML_IO_Mess.f90 /home/jel/TestEnv/crysfml/Src/CFML_Profile_TOF.f90 /home/jel/TestEnv/crysfml/Src/CFML_Profile_Finger.f90 /home/jel/TestEnv/crysfml/Src/CFML_Profile_Functs.f90 /home/jel/TestEnv/crysfml/Src/CFML_Math_3D.f90 /home/jel/TestEnv/crysfml/Src/CFML_Optimization.f90 /home/jel/TestEnv/crysfml/Src/CFML_LSQ_TypeDef.f90 /home/jel/TestEnv/crysfml/Src/CFML_Optimization_LSQ.f90 /home/jel/TestEnv/crysfml/Src/CFML_Sym_Table.f90 /home/jel/TestEnv/crysfml/Src/CFML_Chem_Scatt.f90 /home/jel/TestEnv/crysfml/Src/CFML_Diffpatt.f90 /home/jel/TestEnv/crysfml/Src/CFML_Bonds_Table.f90 /home/jel/TestEnv/crysfml/Src/CFML_Cryst_Types.f90 /home/jel/TestEnv/crysfml/Src/CFML_Symmetry.f90 /home/jel/TestEnv/crysfml/Src/CFML_ILL_Instrm_Data.f90 /home/jel/TestEnv/crysfml/Src/CFML_EoS_Mod.f90 /home/jel/TestEnv/crysfml/Src/CFML_Reflct_Util.f90 /home/jel/TestEnv/crysfml/Src/CFML_Atom_Mod.f90 /home/jel/TestEnv/crysfml/Src/CFML_Export_Vtk.f90 /home/jel/TestEnv/crysfml/Src/CFML_Sfac.f90 /home/jel/TestEnv/crysfml/Src/CFML_Geom_Calc.f90 /home/jel/TestEnv/crysfml/Src/CFML_Propagk.f90 /home/jel/TestEnv/crysfml/Src/CFML_Maps.f90 /home/jel/TestEnv/crysfml/Src/CFML_Molecules.f90 /home/jel/TestEnv/crysfml/Src/CFML_SXTAL_Geom.f90 /home/jel/TestEnv/crysfml/Src/CFML_Conf_Calc.f90 /home/jel/TestEnv/crysfml/Src/CFML_Form_CIF.f90 /home/jel/TestEnv/crysfml/Src/CFML_Optimization_SAn.f90 /home/jel/TestEnv/crysfml/Src/CFML_MagSymm.f90 /home/jel/TestEnv/crysfml/Src/CFML_Msfac.f90 /home/jel/TestEnv/crysfml/Src/CFML_Polar.f90 /home/jel/TestEnv/crysfml/Src/CFML_Refcodes.f90 
#do
#	echo "Processing $f"
#	~/fortwrap.py -d wrap $f
#done
#~/fortwrap.py -c gfortran --file-list=./fwlist -g wrap
# compile Fortran --> C++ headers
#gfortran wrap/CppWrappers.f90 -fPIC -shared -o wrap/CppWrappers.o -lstdc++
# build swig interfaces
# cd wrap
