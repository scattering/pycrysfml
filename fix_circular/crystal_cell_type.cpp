/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "crystal_cell_type.h"

// Constructor:
crystal_cell_type::crystal_cell_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_crystal_cell_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
crystal_cell_type::~crystal_cell_type() {
  __cppwrappers_MOD_deallocate1_crystal_cell_type(data_ptr); // Deallocate Fortran derived type
}

float crystal_cell_type::cell_volume_sigma(void) {
  return __cfml_crystal_metrics_MOD_cell_volume_sigma(data_ptr);
}

float crystal_cell_type::u_equiv(const std::vector<float>* th_u) {
  return __cfml_crystal_metrics_MOD_u_equiv(data_ptr, &(*th_u)[0]);
}

void crystal_cell_type::change_setting_cell(const FortranMatrix<float> *mat, crystal_cell_type* celln, const char* matkind) {
  int matkind_len__ = 0;
  if (matkind) matkind_len__ = strlen(matkind); // Protect Optional args
  __cfml_crystal_metrics_MOD_change_setting_cell(data_ptr, mat->data, celln->data_ptr, matkind, matkind_len__);
}

void crystal_cell_type::get_cryst_family(std::string *car_family, std::string *car_symbol, std::string *car_system) {
  int car_system_len__ = 0;
  if (car_system) car_system_len__ = car_system->length();
  // Declare memory to store output character data
  char car_system_c__[car_system_len__+1];
  car_system_c__[car_system_len__] = '\0';
  int car_symbol_len__ = 0;
  if (car_symbol) car_symbol_len__ = car_symbol->length();
  // Declare memory to store output character data
  char car_symbol_c__[car_symbol_len__+1];
  car_symbol_c__[car_symbol_len__] = '\0';
  int car_family_len__ = 0;
  if (car_family) car_family_len__ = car_family->length();
  // Declare memory to store output character data
  char car_family_c__[car_family_len__+1];
  car_family_c__[car_family_len__] = '\0';
  __cfml_crystal_metrics_MOD_get_cryst_family(data_ptr, car_family ? car_family_c__ : NULL, car_symbol ? car_symbol_c__ : NULL, car_system ? car_system_c__ : NULL, car_family_len__, car_symbol_len__, car_system_len__);
  if (car_system) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=car_family_len__-1; car_system_c__[i]==' '; i--) car_system_c__[i] = '\0';
    car_system->assign(car_system_c__);
  }
  if (car_symbol) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=car_family_len__-1; car_symbol_c__[i]==' '; i--) car_symbol_c__[i] = '\0';
    car_symbol->assign(car_symbol_c__);
  }
  if (car_family) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=car_family_len__-1; car_family_c__[i]==' '; i--) car_family_c__[i] = '\0';
    car_family->assign(car_family_c__);
  }
}

void crystal_cell_type::get_deriv_orth_cell(float* de_orthcell, const char* cartype) {
  // Create C array for Fortran input string data
  char cartype_c__[1+1];
  if (cartype) {
    int i;
    strncpy(cartype_c__, cartype, 1+1); cartype_c__[1+1] = 0; // strncpy protects in case cartype is too long
    for (i=strlen(cartype_c__); i<1+1; i++) cartype_c__[i] = ' '; // Add whitespace for Fortran
  }
  __cfml_crystal_metrics_MOD_get_deriv_orth_cell(data_ptr, de_orthcell, cartype ? cartype_c__ : NULL, 1);
}

void crystal_cell_type::get_transfm_matrix(crystal_cell_type* cellb, FortranMatrix<float> *trm, int* ok, const float* tol) {
  __cfml_crystal_metrics_MOD_get_transfm_matrix(data_ptr, cellb->data_ptr, trm->data, ok, tol);
}

void crystal_cell_type::get_twofold_axes(float tol, twofold_axes_type* twofold) {
  __cfml_crystal_metrics_MOD_get_twofold_axes(data_ptr, &tol, twofold);
}

void crystal_cell_type::read_bin_crystal_cell(int lun, int* ok) {
  __cfml_crystal_metrics_MOD_read_bin_crystal_cell(data_ptr, &lun, ok);
}

void crystal_cell_type::write_bin_crystal_cell(int lun) {
  __cfml_crystal_metrics_MOD_write_bin_crystal_cell(data_ptr, &lun);
}

void crystal_cell_type::write_crystal_cell(const int* lun) {
  __cfml_crystal_metrics_MOD_write_crystal_cell(data_ptr, lun);
}

void crystal_cell_type::atom_uequi_list(atom_list_type* ac) {
  __cfml_atom_typedef_MOD_atom_uequi_list(data_ptr, ac);
}

void crystal_cell_type::unitcell_to_pdbfile(space_group_type* spaceg, atom_list_type* atom_list, const char* filename) {
  int filename_len__ = 0;
  if (filename) filename_len__ = strlen(filename); // Protect Optional args
  __cfml_export_vtk_MOD_unitcell_to_pdbfile(data_ptr, spaceg->data_ptr, atom_list, filename, filename_len__);
}

void crystal_cell_type::distance_and_sigma(const float* derm, const std::vector<float>* x0, const std::vector<float>* x1, const std::vector<float>* s0, const std::vector<float>* s1, float* dis, float* s) {
  __cfml_geometry_calc_MOD_distance_and_sigma(data_ptr, derm, &(*x0)[0], &(*x1)[0], &(*s0)[0], &(*s1)[0], dis, s);
}

void crystal_cell_type::genb(FortranMatrix<float> *b) {
  __cfml_geometry_sxtal_MOD_genb(data_ptr, b->data);
}

void crystal_cell_type::calc_magnetic_strf_miv(magsymm_k_type* mgp, matom_list_type* atm, magh_type* mh) {
  __cfml_magnetic_structure_factors_MOD_calc_magnetic_strf_miv(data_ptr, mgp->data_ptr, atm->data_ptr, mh->data_ptr);
}

void crystal_cell_type::calc_magnetic_strf_miv_dom(magsymm_k_type* mgp, matom_list_type* atm, magnetic_domain_type* mag_dom, maghd_type* mh) {
  __cfml_magnetic_structure_factors_MOD_calc_magnetic_strf_miv_dom(data_ptr, mgp->data_ptr, atm->data_ptr, mag_dom->data_ptr, mh->data_ptr);
}

void crystal_cell_type::gen_satellites(magsymm_k_type* grp, float smax, magh_list_type* h, const int* ord, const int* powder, reflection_list_type* hkl) {
  __cfml_magnetic_structure_factors_MOD_gen_satellites(data_ptr, grp->data_ptr, &smax, h->data_ptr, ord, powder, hkl ? hkl->data_ptr : NULL);
}

