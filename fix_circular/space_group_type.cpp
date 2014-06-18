/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "space_group_type.h"

// Constructor:
space_group_type::space_group_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_space_group_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
space_group_type::~space_group_type() {
  __cppwrappers_MOD_deallocate1_space_group_type(data_ptr); // Deallocate Fortran derived type
}

int space_group_type::spgr_equal(space_group_type* spacegroup2) {
  return __cfml_crystallographic_symmetry_MOD_spgr_equal(data_ptr, spacegroup2->data_ptr);
}

void space_group_type::get_hallsymb_from_gener(std::string *spaceh) {
  int spaceh_len__ = 0;
  if (spaceh) spaceh_len__ = spaceh->length();
  // Declare memory to store output character data
  char spaceh_c__[spaceh_len__+1];
  spaceh_c__[spaceh_len__] = '\0';
  __cfml_crystallographic_symmetry_MOD_get_hallsymb_from_gener(data_ptr, spaceh ? spaceh_c__ : NULL, spaceh_len__);
  if (spaceh) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=spaceh_len__-1; spaceh_c__[i]==' '; i--) spaceh_c__[i] = '\0';
    spaceh->assign(spaceh_c__);
  }
}

void space_group_type::get_laue_pg(std::string *laue_car, std::string *point_car) {
  int point_car_len__ = 0;
  if (point_car) point_car_len__ = point_car->length();
  // Declare memory to store output character data
  char point_car_c__[point_car_len__+1];
  point_car_c__[point_car_len__] = '\0';
  int laue_car_len__ = 0;
  if (laue_car) laue_car_len__ = laue_car->length();
  // Declare memory to store output character data
  char laue_car_c__[laue_car_len__+1];
  laue_car_c__[laue_car_len__] = '\0';
  __cfml_crystallographic_symmetry_MOD_get_laue_pg(data_ptr, laue_car ? laue_car_c__ : NULL, point_car ? point_car_c__ : NULL, laue_car_len__, point_car_len__);
  if (point_car) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=laue_car_len__-1; point_car_c__[i]==' '; i--) point_car_c__[i] = '\0';
    point_car->assign(point_car_c__);
  }
  if (laue_car) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=laue_car_len__-1; laue_car_c__[i]==' '; i--) laue_car_c__[i] = '\0';
    laue_car->assign(laue_car_c__);
  }
}

void space_group_type::read_bin_spacegroup(int lun, int* ok) {
  __cfml_crystallographic_symmetry_MOD_read_bin_spacegroup(data_ptr, &lun, ok);
}

void space_group_type::write_bin_spacegroup(int lun) {
  __cfml_crystallographic_symmetry_MOD_write_bin_spacegroup(data_ptr, &lun);
}

void space_group_type::write_spacegroup(const int* iunit, const int* full) {
  __cfml_crystallographic_symmetry_MOD_write_spacegroup(data_ptr, iunit, full);
}

void space_group_type::write_asu(const int* iunit) {
  __cfml_reflections_utilities_MOD_write_asu(data_ptr, iunit);
}

void space_group_type::atlist1_extencell_atlist2(atom_list_type* a, atom_list_type* c, int conven) {
  __cfml_atom_typedef_MOD_atlist1_extencell_atlist2(data_ptr, a, c, &conven);
}

void space_group_type::set_atom_equiv_list(crystal_cell_type* cell, atom_list_type* a, atom_equiv_list_type* ate, const int* lun) {
  __cfml_atom_typedef_MOD_set_atom_equiv_list(data_ptr, cell, a, ate, lun);
}

void space_group_type::set_new_asymunit(atom_equiv_list_type* ate, const FortranMatrix<float> *mat, const std::vector<float>* orig, atom_list_type* a_n, const char* matkind, const char* debug) {
  int debug_len__ = 0;
  if (debug) debug_len__ = strlen(debug); // Protect Optional args
  int matkind_len__ = 0;
  if (matkind) matkind_len__ = strlen(matkind); // Protect Optional args
  __cfml_geometry_calc_MOD_set_new_asymunit(data_ptr, ate->data_ptr, mat->data, &(*orig)[0], a_n, matkind, debug, matkind_len__, debug_len__);
}

void space_group_type::set_orbits_inlist(point_list_type* pl) {
  __cfml_geometry_calc_MOD_set_orbits_inlist(data_ptr, pl->data_ptr);
}

