/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef SPACE_GROUP_TYPE_H_
#define SPACE_GROUP_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "atom_equiv_list_type.h"
#include "FortranMatrix.h"
#include "point_list_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_space_group_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_space_group_type(ADDRESS caddr);
  int __cfml_crystallographic_symmetry_MOD_spgr_equal(ADDRESS spacegroup1, ADDRESS spacegroup2);
  void __cfml_crystallographic_symmetry_MOD_get_hallsymb_from_gener(ADDRESS spacegroup, char* spaceh, int spaceh_len__);
  void __cfml_crystallographic_symmetry_MOD_get_laue_pg(ADDRESS spacegroup, char* laue_car, char* point_car, int laue_car_len__, int point_car_len__);
  void __cfml_crystallographic_symmetry_MOD_read_bin_spacegroup(ADDRESS spg, int* lun, int* ok);
  void __cfml_crystallographic_symmetry_MOD_write_bin_spacegroup(ADDRESS spg, int* lun);
  void __cfml_crystallographic_symmetry_MOD_write_spacegroup(ADDRESS spacegroup, const int* iunit, const int* full);
  void __cfml_reflections_utilities_MOD_write_asu(ADDRESS spacegroup, const int* iunit);
  void __cfml_atom_typedef_MOD_atlist1_extencell_atlist2(ADDRESS spg, ADDRESS a, ADDRESS c, int* conven);
  void __cfml_atom_typedef_MOD_set_atom_equiv_list(ADDRESS spg, ADDRESS cell, ADDRESS a, ADDRESS ate, const int* lun);
  void __cfml_geometry_calc_MOD_set_new_asymunit(ADDRESS spgn, ADDRESS ate, float* mat, const float orig[], ADDRESS a_n, const char* matkind, const char* debug, int matkind_len__, int debug_len__);
  void __cfml_geometry_calc_MOD_set_orbits_inlist(ADDRESS spg, ADDRESS pl);
}
#endif // SWIG
class atom_list_type;
class crystal_cell_type;
class space_group_type {

public:
  space_group_type();
  ~space_group_type();

  int spgr_equal(space_group_type* spacegroup2);

/*! \param[out] spaceh OPTIONAL
*/
  void get_hallsymb_from_gener(std::string *spaceh=NULL);

  void get_laue_pg(std::string *laue_car, std::string *point_car);

  void read_bin_spacegroup(int lun, int* ok);

  void write_bin_spacegroup(int lun);

/*! \param[in] iunit OPTIONAL
 *
 *  \param[in] full OPTIONAL
*/
  void write_spacegroup(const int* iunit=NULL, const int* full=NULL);

/*! \param[in] iunit OPTIONAL
*/
  void write_asu(const int* iunit=NULL);

  void atlist1_extencell_atlist2(atom_list_type* a, atom_list_type* c, int conven);

/*! \param[in] lun OPTIONAL
*/
  void set_atom_equiv_list(crystal_cell_type* cell, atom_list_type* a, atom_equiv_list_type* ate, const int* lun=NULL);

/*! \param[in] orig ARRAY
 *
 *  \param[in] matkind OPTIONAL
 *
 *  \param[in] debug OPTIONAL
*/
  void set_new_asymunit(atom_equiv_list_type* ate, const FortranMatrix<float> *mat, const std::vector<float>* orig, atom_list_type* a_n, const char* matkind=NULL, const char* debug=NULL);

  void set_orbits_inlist(point_list_type* pl);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* SPACE_GROUP_TYPE_H_ */
