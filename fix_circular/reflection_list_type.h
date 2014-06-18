/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef REFLECTION_LIST_TYPE_H_
#define REFLECTION_LIST_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "space_group_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_reflection_list_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_reflection_list_type(ADDRESS caddr);
  void __cfml_reflections_utilities_MOD_init_reflist(ADDRESS reflex, const int* n);
  void __cfml_reflections_utilities_MOD_write_reflist_info(ADDRESS rfl, const int* iunit, const char* mode, int mode_len__);
  void __cfml_structure_factors_MOD_init_calc_strfactors(ADDRESS reflex, ADDRESS atm, ADDRESS grp, const char* mode, const float* lambda, const int* lun, int mode_len__);
  void __cfml_structure_factors_MOD_init_structure_factors(ADDRESS reflex, ADDRESS atm, ADDRESS grp, const char* mode, const float* lambda, const int* lun, int mode_len__);
}
#endif // SWIG
class atom_list_type;
class reflection_list_type {

public:
  reflection_list_type();
  ~reflection_list_type();

/*! \param[in] n OPTIONAL
*/
  void init_reflist(const int* n=NULL);

/*! \param[in] iunit OPTIONAL
 *
 *  \param[in] mode OPTIONAL
*/
  void write_reflist_info(const int* iunit=NULL, const char* mode=NULL);

/*! \param[in] mode OPTIONAL
 *
 *  \param[in] lambda OPTIONAL
 *
 *  \param[in] lun OPTIONAL
*/
  void init_calc_strfactors(atom_list_type* atm, space_group_type* grp, const char* mode=NULL, const float* lambda=NULL, const int* lun=NULL);

/*! \param[in] mode OPTIONAL
 *
 *  \param[in] lambda OPTIONAL
 *
 *  \param[in] lun OPTIONAL
*/
  void init_structure_factors(atom_list_type* atm, space_group_type* grp, const char* mode=NULL, const float* lambda=NULL, const int* lun=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* REFLECTION_LIST_TYPE_H_ */
