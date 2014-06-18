/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef MATOM_LIST_TYPE_H_
#define MATOM_LIST_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "magsymm_k_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_matom_list_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_matom_list_type(ADDRESS caddr);
  void __cfml_atom_typedef_MOD_deallocate_matom_list(ADDRESS a);
  void __cfml_magnetic_structure_factors_MOD_mag_structure_factors(ADDRESS atm, ADDRESS grp, ADDRESS reflex);
}
#endif // SWIG
class magh_list_type;
class matom_list_type {

public:
  matom_list_type();
  ~matom_list_type();

  void deallocate_matom_list(void);

  void mag_structure_factors(magsymm_k_type* grp, magh_list_type* reflex);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* MATOM_LIST_TYPE_H_ */
