/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "atoms_cell_type.h"

// Constructor:
atoms_cell_type::atoms_cell_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_atoms_cell_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
atoms_cell_type::~atoms_cell_type() {
  __cppwrappers_MOD_deallocate1_atoms_cell_type(data_ptr); // Deallocate Fortran derived type
}

void atoms_cell_type::atoms_cell_to_list(atom_list_type* a) {
 // __cfml_atom_typedef_MOD_atoms_cell_to_list(data_ptr, a->data_ptr);
  __cfml_atom_typedef_MOD_atoms_cell_to_list(data_ptr, a);
}

void atoms_cell_type::deallocate_atoms_cell(void) {
  __cfml_atom_typedef_MOD_deallocate_atoms_cell(data_ptr);
}

