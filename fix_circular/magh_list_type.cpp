/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "magh_list_type.h"

// Constructor:
magh_list_type::magh_list_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_magh_list_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
magh_list_type::~magh_list_type() {
  __cppwrappers_MOD_deallocate1_magh_list_type(data_ptr); // Deallocate Fortran derived type
}

void magh_list_type::calc_mag_interaction_vector(crystal_cell_type* cell) {
  __cfml_magnetic_structure_factors_MOD_calc_mag_interaction_vector(data_ptr, cell);
}

void magh_list_type::init_mag_structure_factors(matom_list_type* atm, magsymm_k_type* grp, const int* lun) {
  __cfml_magnetic_structure_factors_MOD_init_mag_structure_factors(data_ptr, atm->data_ptr, grp->data_ptr, lun);
}

