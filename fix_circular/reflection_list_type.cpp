/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "reflection_list_type.h"

// Constructor:
reflection_list_type::reflection_list_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_reflection_list_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
reflection_list_type::~reflection_list_type() {
  __cppwrappers_MOD_deallocate1_reflection_list_type(data_ptr); // Deallocate Fortran derived type
}

void reflection_list_type::init_reflist(const int* n) {
  __cfml_reflections_utilities_MOD_init_reflist(data_ptr, n);
}

void reflection_list_type::write_reflist_info(const int* iunit, const char* mode) {
  int mode_len__ = 0;
  if (mode) mode_len__ = strlen(mode); // Protect Optional args
  __cfml_reflections_utilities_MOD_write_reflist_info(data_ptr, iunit, mode, mode_len__);
}

void reflection_list_type::init_calc_strfactors(atom_list_type* atm, space_group_type* grp, const char* mode, const float* lambda, const int* lun) {
  int mode_len__ = 0;
  if (mode) mode_len__ = strlen(mode); // Protect Optional args
  __cfml_structure_factors_MOD_init_calc_strfactors(data_ptr, atm, grp->data_ptr, mode, lambda, lun, mode_len__);
}

void reflection_list_type::init_structure_factors(atom_list_type* atm, space_group_type* grp, const char* mode, const float* lambda, const int* lun) {
  int mode_len__ = 0;
  if (mode) mode_len__ = strlen(mode); // Protect Optional args
  __cfml_structure_factors_MOD_init_structure_factors(data_ptr, atm, grp->data_ptr, mode, lambda, lun, mode_len__);
}

