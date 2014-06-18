/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "molecular_crystal_type.h"

// Constructor:
molecular_crystal_type::molecular_crystal_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_molecular_crystal_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
molecular_crystal_type::~molecular_crystal_type() {
  __cppwrappers_MOD_deallocate1_molecular_crystal_type(data_ptr); // Deallocate Fortran derived type
}

void molecular_crystal_type::molcrys_to_atomlist(atom_list_type* atm) {
  __cfml_molecular_crystals_MOD_molcrys_to_atomlist(data_ptr, atm);
}

void molecular_crystal_type::write_molecular_crystal(const int* lun) {
  __cfml_molecular_crystals_MOD_write_molecular_crystal(data_ptr, lun);
}

