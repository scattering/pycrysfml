/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef MAGH_LIST_TYPE_H_
#define MAGH_LIST_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "magsymm_k_type.h"
#include "matom_list_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_magh_list_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_magh_list_type(ADDRESS caddr);
  void __cfml_magnetic_structure_factors_MOD_calc_mag_interaction_vector(ADDRESS reflex, ADDRESS cell);
  void __cfml_magnetic_structure_factors_MOD_init_mag_structure_factors(ADDRESS reflex, ADDRESS atm, ADDRESS grp, const int* lun);
}
#endif // SWIG
class crystal_cell_type;
class magh_list_type {

public:
  magh_list_type();
  ~magh_list_type();

  void calc_mag_interaction_vector(crystal_cell_type* cell);

/*! \param[in] lun OPTIONAL
*/
  void init_mag_structure_factors(matom_list_type* atm, magsymm_k_type* grp, const int* lun=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* MAGH_LIST_TYPE_H_ */
