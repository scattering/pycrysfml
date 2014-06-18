/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#include <cstring> // For strcpy
#include "molecule_type.h"

// Constructor:
molecule_type::molecule_type() {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate1_molecule_type(&data_ptr); // Allocate Fortran derived type
  initialized = false;
}

// Destructor:
molecule_type::~molecule_type() {
  __cppwrappers_MOD_deallocate1_molecule_type(data_ptr); // Deallocate Fortran derived type
}

void molecule_type::cartesian_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_cartesian_to_fractional(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::cartesian_to_spherical(molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_cartesian_to_spherical(data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::cartesian_to_zmatrix(molecule_type* newmolecule, crystal_cell_type* cell, const float* d_min, const float* d_max) {
  __cfml_molecular_crystals_MOD_cartesian_to_zmatrix(data_ptr, newmolecule ? newmolecule->data_ptr : NULL, cell ? cell->data_ptr : NULL, d_min, d_max);
}

void molecule_type::fix_reference(molecule_type* newmolecule, const int* natom_o, const int* natom_x, const int* natom_xy) {
  __cfml_molecular_crystals_MOD_fix_reference(data_ptr, newmolecule ? newmolecule->data_ptr : NULL, natom_o, natom_x, natom_xy);
}

void molecule_type::fix_orient_cartesian(molecule_type* newmolecule, const int* natom_o, const int* natom_x, const int* natom_xy, FortranMatrix<float> *mat) {
  __cfml_molecular_crystals_MOD_fix_orient_cartesian(data_ptr, newmolecule ? newmolecule->data_ptr : NULL, natom_o, natom_x, natom_xy, mat ? mat->data : NULL);
}

void molecule_type::fractional_to_cartesian(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_fractional_to_cartesian(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::fractional_to_spherical(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_fractional_to_spherical(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::fractional_to_zmatrix(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_fractional_to_zmatrix(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::init_molecule(const int* natm) {
  __cfml_molecular_crystals_MOD_init_molecule(data_ptr, natm);
}

void molecule_type::molec_to_atomlist(atom_list_type* atm, const char* coor_type, crystal_cell_type* cell) {
  int coor_type_len__ = 0;
  if (coor_type) coor_type_len__ = strlen(coor_type); // Protect Optional args
  __cfml_molecular_crystals_MOD_molec_to_atomlist(data_ptr, atm, coor_type, cell ? cell->data_ptr : NULL, coor_type_len__);
}

void molecule_type::spherical_to_cartesian(molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_spherical_to_cartesian(data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::spherical_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_spherical_to_fractional(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::spherical_to_zmatrix(molecule_type* newmolecule, crystal_cell_type* cell) {
  __cfml_molecular_crystals_MOD_spherical_to_zmatrix(data_ptr, newmolecule ? newmolecule->data_ptr : NULL, cell ? cell->data_ptr : NULL);
}

void molecule_type::write_molecule(const int* lun) {
  __cfml_molecular_crystals_MOD_write_molecule(data_ptr, lun);
}

void molecule_type::zmatrix_to_cartesian(molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_zmatrix_to_cartesian(data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::zmatrix_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_zmatrix_to_fractional(data_ptr, cell->data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

void molecule_type::zmatrix_to_spherical(molecule_type* newmolecule) {
  __cfml_molecular_crystals_MOD_zmatrix_to_spherical(data_ptr, newmolecule ? newmolecule->data_ptr : NULL);
}

