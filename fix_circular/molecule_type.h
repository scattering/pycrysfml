/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef MOLECULE_TYPE_H_
#define MOLECULE_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "FortranMatrix.h"
#include "crystal_cell_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_molecule_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_molecule_type(ADDRESS caddr);
  void __cfml_molecular_crystals_MOD_cartesian_to_fractional(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_cartesian_to_spherical(ADDRESS molecule, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_cartesian_to_zmatrix(ADDRESS molecule, ADDRESS newmolecule, ADDRESS cell, const float* d_min, const float* d_max);
  void __cfml_molecular_crystals_MOD_fix_reference(ADDRESS molecule, ADDRESS newmolecule, const int* natom_o, const int* natom_x, const int* natom_xy);
  void __cfml_molecular_crystals_MOD_fix_orient_cartesian(ADDRESS molecule, ADDRESS newmolecule, const int* natom_o, const int* natom_x, const int* natom_xy, float* mat);
  void __cfml_molecular_crystals_MOD_fractional_to_cartesian(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_fractional_to_spherical(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_fractional_to_zmatrix(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_init_molecule(ADDRESS molecule, const int* natm);
  void __cfml_molecular_crystals_MOD_molec_to_atomlist(ADDRESS molec, ADDRESS atm, const char* coor_type, ADDRESS cell, int coor_type_len__);
  void __cfml_molecular_crystals_MOD_spherical_to_cartesian(ADDRESS molecule, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_spherical_to_fractional(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_spherical_to_zmatrix(ADDRESS molecule, ADDRESS newmolecule, ADDRESS cell);
  void __cfml_molecular_crystals_MOD_write_molecule(ADDRESS molecule, const int* lun);
  void __cfml_molecular_crystals_MOD_zmatrix_to_cartesian(ADDRESS molecule, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_zmatrix_to_fractional(ADDRESS molecule, ADDRESS cell, ADDRESS newmolecule);
  void __cfml_molecular_crystals_MOD_zmatrix_to_spherical(ADDRESS molecule, ADDRESS newmolecule);
}
#endif // SWIG
class atom_list_type;
class molecule_type {

public:
  molecule_type();
  ~molecule_type();

/*! \param[out] newmolecule OPTIONAL
*/
  void cartesian_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void cartesian_to_spherical(molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
 *
 *  \param[in] cell OPTIONAL
 *
 *  \param[in] d_min OPTIONAL
 *
 *  \param[in] d_max OPTIONAL
*/
  void cartesian_to_zmatrix(molecule_type* newmolecule=NULL, crystal_cell_type* cell=NULL, const float* d_min=NULL, const float* d_max=NULL);

/*! \param[out] newmolecule OPTIONAL
 *
 *  \param[in] natom_o OPTIONAL
 *
 *  \param[in] natom_x OPTIONAL
 *
 *  \param[in] natom_xy OPTIONAL
*/
  void fix_reference(molecule_type* newmolecule=NULL, const int* natom_o=NULL, const int* natom_x=NULL, const int* natom_xy=NULL);

/*! \param[out] newmolecule OPTIONAL
 *
 *  \param[in] natom_o OPTIONAL
 *
 *  \param[in] natom_x OPTIONAL
 *
 *  \param[in] natom_xy OPTIONAL
 *
 *  \param[out] mat OPTIONAL
*/
  void fix_orient_cartesian(molecule_type* newmolecule=NULL, const int* natom_o=NULL, const int* natom_x=NULL, const int* natom_xy=NULL, FortranMatrix<float> *mat=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void fractional_to_cartesian(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void fractional_to_spherical(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void fractional_to_zmatrix(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[in] natm OPTIONAL
*/
  void init_molecule(const int* natm=NULL);

/*! \param[in] coor_type OPTIONAL
 *
 *  \param[in] cell OPTIONAL
*/
  void molec_to_atomlist(atom_list_type* atm, const char* coor_type=NULL, crystal_cell_type* cell=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void spherical_to_cartesian(molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void spherical_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
 *
 *  \param[in] cell OPTIONAL
*/
  void spherical_to_zmatrix(molecule_type* newmolecule=NULL, crystal_cell_type* cell=NULL);

/*! \param[in] lun OPTIONAL
*/
  void write_molecule(const int* lun=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void zmatrix_to_cartesian(molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void zmatrix_to_fractional(crystal_cell_type* cell, molecule_type* newmolecule=NULL);

/*! \param[out] newmolecule OPTIONAL
*/
  void zmatrix_to_spherical(molecule_type* newmolecule=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* MOLECULE_TYPE_H_ */
