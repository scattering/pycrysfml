/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef ATOMS_CELL_TYPE_H_
#define ATOMS_CELL_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_atoms_cell_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_atoms_cell_type(ADDRESS caddr);
  void __cfml_atom_typedef_MOD_atoms_cell_to_list(ADDRESS ac, ADDRESS a);
  void __cfml_atom_typedef_MOD_deallocate_atoms_cell(ADDRESS ac);
}
#endif // SWIG
class atom_list_type;
class atoms_cell_type {

public:
  atoms_cell_type();
  ~atoms_cell_type();

  void atoms_cell_to_list(atom_list_type* a);

  void deallocate_atoms_cell(void);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* ATOMS_CELL_TYPE_H_ */
