add_cmplx.py
========

injects complex struct into C++ wrapper to add imaginary number support

build.sh
========

Main build script. Calls forwrap related scripts, compiles library, generates swig interface, etc

clean.sh
========

Cleans up the working directory

fix_deps.py
========

fixes parsing issues with global dependencies

fix_makefile_deps.py
========

Hacks EoS and SAn modules out of makefile to avoid undefined symbol errors (SAn --> write_fst).

fix_type_decl.py
========

converts library declarations of derrived types to a FortWrap compatible form

fortwrap.py
========

Modified version of FortWrap for use with CrysFML Library

gen_list.py
========

Used to generate list of files to be wrapped

gen_help_file.py
========

outputs python help to txt file

gen_cpp_methods.py - gen_fort_methods.py
========

Used to generate manual wrappers. Not part of build script.

bfiles.txt
========

list of files in the order from the original cmake build script used by gen_list.py
TODO: combine gen_list and f90_deps to eliminate this file

ff
========

list of derived type interfaces used by gen_cpp_methods and gen_fort_methods
syntax is:

	# derived_type_name!module_name
	<Fortran declarations of type attributes>

TODO: Automatically generate this file

f90_deps.py
========

From original pycrysfml (Dylan Quintana)
used to generate Makefile.deps for compiling shared object library
