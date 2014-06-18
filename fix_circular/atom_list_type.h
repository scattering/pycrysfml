/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef ATOM_LIST_TYPE_H_
#define ATOM_LIST_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "reflection_list_type.h"
#include "space_group_type.h"
#include "atoms_cell_type.h"
#include "magnetic_domain_type.h"
#include "crystal_cell_type.h"
#include "matom_list_type.h"
#include "nonatomic_parameter_list_type.h"
#include "molecule_type.h"
#include "molecular_crystal_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_atom_list_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_atom_list_type(ADDRESS caddr);
  void __cfml_atom_typedef_MOD_atom_list_to_cell(ADDRESS a, ADDRESS ac);
  void __cfml_atom_typedef_MOD_copy_atom_list(ADDRESS a, ADDRESS ac);
  void __cfml_atom_typedef_MOD_deallocate_atom_list(ADDRESS a);
  void __cfml_atom_typedef_MOD_read_bin_atom_list(ADDRESS ats, int* lun, int* ok);
  void __cfml_atom_typedef_MOD_write_atom_list(ADDRESS ats, const int* level, const int* lun, ADDRESS cell);
  void __cfml_atom_typedef_MOD_write_atoms_cfl(ADDRESS ats, const int* lun, ADDRESS cell);
  void __cfml_atom_typedef_MOD_write_bin_atom_list(ADDRESS ats, int* lun);
  void __cfml_structure_factors_MOD_init_calc_hkl_strfactors(ADDRESS atm, const char* mode, const float* lambda, const int* lun, int mode_len__);
  void __cfml_structure_factors_MOD_structure_factors(ADDRESS atm, ADDRESS grp, ADDRESS reflex, const char* mode, const float* lambda, int mode_len__);
  void __cfml_keywords_code_parser_MOD_init_refcodes(ADDRESS fatom, ADDRESS fmatom, ADDRESS mag_dom, ADDRESS molcrys, ADDRESS molec, ADDRESS model);
  void __cfml_keywords_code_parser_MOD_write_restraints_obscalc(ADDRESS a, const int* iunit);
}
#endif // SWIG

class atom_list_type {

public:
  atom_list_type();
  ~atom_list_type();

  void atom_list_to_cell(atoms_cell_type* ac);

  void copy_atom_list(atom_list_type* ac);

  void deallocate_atom_list(void);

  void read_bin_atom_list(int lun, int* ok);

/*! \param[in] level OPTIONAL
 *
 *  \param[in] lun OPTIONAL
 *
 *  \param[in] cell OPTIONAL
*/
  void write_atom_list(const int* level=NULL, const int* lun=NULL, crystal_cell_type* cell=NULL);

/*! \param[in] lun OPTIONAL
 *
 *  \param[in] cell OPTIONAL
*/
  void write_atoms_cfl(const int* lun=NULL, crystal_cell_type* cell=NULL);

  void write_bin_atom_list(int lun);

/*! \param[in] mode OPTIONAL
 *
 *  \param[in] lambda OPTIONAL
 *
 *  \param[in] lun OPTIONAL
*/
  void init_calc_hkl_strfactors(const char* mode=NULL, const float* lambda=NULL, const int* lun=NULL);

/*! \param[in] mode OPTIONAL
 *
 *  \param[in] lambda OPTIONAL
*/
  void structure_factors(space_group_type* grp, reflection_list_type* reflex, const char* mode=NULL, const float* lambda=NULL);

/*! \param fatom OPTIONAL
 *
 *  \param fmatom OPTIONAL
 *
 *  \param mag_dom OPTIONAL
 *
 *  \param molcrys OPTIONAL
 *
 *  \param molec OPTIONAL
 *
 *  \param model OPTIONAL
*/
  void init_refcodes(matom_list_type* fmatom=NULL, magnetic_domain_type* mag_dom=NULL, molecular_crystal_type* molcrys=NULL, molecule_type* molec=NULL, nonatomic_parameter_list_type* model=NULL);

/*! \param[in] iunit OPTIONAL
*/
  void write_restraints_obscalc(const int* iunit=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* ATOM_LIST_TYPE_H_ */
