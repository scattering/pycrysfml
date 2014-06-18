/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "matom_list_type.h"

// Constructor:
matom_list_type::matom_list_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_matom_list_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
matom_list_type::~matom_list_type() {
  __cppwrappers_MOD_deallocate1_matom_list_type(data_ptr); // Deallocate Fortran derived type
}

void matom_list_type::deallocate_matom_list(void) {
  __cfml_atom_typedef_MOD_deallocate_matom_list(data_ptr);
}

void matom_list_type::mag_structure_factors(magsymm_k_type* grp, magh_list_type* reflex) {
  __cfml_magnetic_structure_factors_MOD_mag_structure_factors(data_ptr, grp->data_ptr, reflex);
}

