/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef MOLECULAR_CRYSTAL_TYPE_H_
#define MOLECULAR_CRYSTAL_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_molecular_crystal_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_molecular_crystal_type(ADDRESS caddr);
  void __cfml_molecular_crystals_MOD_molcrys_to_atomlist(ADDRESS molcrys, ADDRESS atm);
  void __cfml_molecular_crystals_MOD_write_molecular_crystal(ADDRESS molcrys, const int* lun);
}
#endif // SWIG
class atom_list_type;
class molecular_crystal_type {

public:
  molecular_crystal_type();
  ~molecular_crystal_type();

  void molcrys_to_atomlist(atom_list_type* atm);

/*! \param[in] lun OPTIONAL
*/
  void write_molecular_crystal(const int* lun=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* MOLECULAR_CRYSTAL_TYPE_H_ */
