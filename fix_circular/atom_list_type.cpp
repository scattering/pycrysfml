/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "atom_list_type.h"

// Constructor:
atom_list_type::atom_list_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_atom_list_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
atom_list_type::~atom_list_type() {
  __cppwrappers_MOD_deallocate1_atom_list_type(data_ptr); // Deallocate Fortran derived type
}

void atom_list_type::atom_list_to_cell(atoms_cell_type* ac) {
  __cfml_atom_typedef_MOD_atom_list_to_cell(data_ptr, ac->data_ptr);
}

void atom_list_type::copy_atom_list(atom_list_type* ac) {
  __cfml_atom_typedef_MOD_copy_atom_list(data_ptr, ac->data_ptr);
}

void atom_list_type::deallocate_atom_list(void) {
  __cfml_atom_typedef_MOD_deallocate_atom_list(data_ptr);
}

void atom_list_type::read_bin_atom_list(int lun, int* ok) {
  __cfml_atom_typedef_MOD_read_bin_atom_list(data_ptr, &lun, ok);
}

void atom_list_type::write_atom_list(const int* level, const int* lun, crystal_cell_type* cell) {
  __cfml_atom_typedef_MOD_write_atom_list(data_ptr, level, lun, cell ? cell->data_ptr : NULL);
}

void atom_list_type::write_atoms_cfl(const int* lun, crystal_cell_type* cell) {
  __cfml_atom_typedef_MOD_write_atoms_cfl(data_ptr, lun, cell ? cell->data_ptr : NULL);
}

void atom_list_type::write_bin_atom_list(int lun) {
  __cfml_atom_typedef_MOD_write_bin_atom_list(data_ptr, &lun);
}

void atom_list_type::init_calc_hkl_strfactors(const char* mode, const float* lambda, const int* lun) {
  int mode_len__ = 0;
  if (mode) mode_len__ = strlen(mode); // Protect Optional args
  __cfml_structure_factors_MOD_init_calc_hkl_strfactors(data_ptr, mode, lambda, lun, mode_len__);
}

void atom_list_type::structure_factors(space_group_type* grp, reflection_list_type* reflex, const char* mode, const float* lambda) {
  int mode_len__ = 0;
  if (mode) mode_len__ = strlen(mode); // Protect Optional args
  __cfml_structure_factors_MOD_structure_factors(data_ptr, grp->data_ptr, reflex->data_ptr, mode, lambda, mode_len__);
}

void atom_list_type::init_refcodes(matom_list_type* fmatom, magnetic_domain_type* mag_dom, molecular_crystal_type* molcrys, molecule_type* molec, nonatomic_parameter_list_type* model) {
  __cfml_keywords_code_parser_MOD_init_refcodes(data_ptr, fmatom ? fmatom->data_ptr : NULL, mag_dom ? mag_dom->data_ptr : NULL, molcrys ? molcrys->data_ptr : NULL, molec ? molec->data_ptr : NULL, model ? model->data_ptr : NULL);
}

void atom_list_type::write_restraints_obscalc(const int* iunit) {
  __cfml_keywords_code_parser_MOD_write_restraints_obscalc(data_ptr, iunit);
}

