/* This source file automatically generated on 2014-06-18 using 
   FortWrap wrapper generator version 1.0.4 */

#ifndef CRYSTAL_CELL_TYPE_H_
#define CRYSTAL_CELL_TYPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <string>
#include <vector>
#include "InterfaceDefs.h"
#include "reflection_list_type.h"
#include "magsymm_k_type.h"
#include "magh_list_type.h"
#include "space_group_type.h"
#include "magnetic_domain_type.h"
#include "matom_list_type.h"
#include "magh_type.h"
#include "FortranMatrix.h"
#include "maghd_type.h"

extern "C" {
  void __cppwrappers_MOD_allocate1_crystal_cell_type(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate1_crystal_cell_type(ADDRESS caddr);
  float __cfml_crystal_metrics_MOD_cell_volume_sigma(ADDRESS cell);
  float __cfml_crystal_metrics_MOD_u_equiv(ADDRESS cell, const float th_u[]);
  void __cfml_crystal_metrics_MOD_change_setting_cell(ADDRESS cell, float* mat, ADDRESS celln, const char* matkind, int matkind_len__);
  void __cfml_crystal_metrics_MOD_get_cryst_family(ADDRESS cell, char* car_family, char* car_symbol, char* car_system, int car_family_len__, int car_symbol_len__, int car_system_len__);
  void __cfml_crystal_metrics_MOD_get_deriv_orth_cell(ADDRESS cellp, float* de_orthcell, const char* cartype, int cartype_len__);
  void __cfml_crystal_metrics_MOD_get_transfm_matrix(ADDRESS cella, ADDRESS cellb, float* trm, int* ok, const float* tol);
  void __cfml_crystal_metrics_MOD_get_twofold_axes(ADDRESS celln, float* tol, ADDRESS twofold);
  void __cfml_crystal_metrics_MOD_read_bin_crystal_cell(ADDRESS celda, int* lun, int* ok);
  void __cfml_crystal_metrics_MOD_write_bin_crystal_cell(ADDRESS celda, int* lun);
  void __cfml_crystal_metrics_MOD_write_crystal_cell(ADDRESS celda, const int* lun);
  void __cfml_atom_typedef_MOD_atom_uequi_list(ADDRESS cell, ADDRESS ac);
  void __cfml_export_vtk_MOD_unitcell_to_pdbfile(ADDRESS cell, ADDRESS spaceg, ADDRESS atom_list, const char* filename, int filename_len__);
  void __cfml_geometry_calc_MOD_distance_and_sigma(ADDRESS cellp, const float* derm, const float x0[], const float x1[], const float s0[], const float s1[], float* dis, float* s);
  void __cfml_geometry_sxtal_MOD_genb(ADDRESS c, float* b);
  void __cfml_magnetic_structure_factors_MOD_calc_magnetic_strf_miv(ADDRESS cell, ADDRESS mgp, ADDRESS atm, ADDRESS mh);
  void __cfml_magnetic_structure_factors_MOD_calc_magnetic_strf_miv_dom(ADDRESS cell, ADDRESS mgp, ADDRESS atm, ADDRESS mag_dom, ADDRESS mh);
  void __cfml_magnetic_structure_factors_MOD_gen_satellites(ADDRESS cell, ADDRESS grp, float* smax, ADDRESS h, const int* ord, const int* powder, ADDRESS hkl);
}
#endif // SWIG
class atom_list_type;
class twofold_axes_type;
class crystal_cell_type {

public:
  crystal_cell_type();
  ~crystal_cell_type();

  float cell_volume_sigma(void);

/*! \param[in] th_u ARRAY
*/
  float u_equiv(const std::vector<float>* th_u);

/*! \param[in] matkind OPTIONAL
*/
  void change_setting_cell(const FortranMatrix<float> *mat, crystal_cell_type* celln, const char* matkind=NULL);

  void get_cryst_family(std::string *car_family, std::string *car_symbol, std::string *car_system);

/*! \param[in] cartype OPTIONAL
*/
  void get_deriv_orth_cell(float* de_orthcell, const char* cartype=NULL);

/*! \param[in] tol OPTIONAL
*/
  void get_transfm_matrix(crystal_cell_type* cellb, FortranMatrix<float> *trm, int* ok, const float* tol=NULL);

  void get_twofold_axes(float tol, twofold_axes_type* twofold);

  void read_bin_crystal_cell(int lun, int* ok);

  void write_bin_crystal_cell(int lun);

/*! \param[in] lun OPTIONAL
*/
  void write_crystal_cell(const int* lun=NULL);

  void atom_uequi_list(atom_list_type* ac);

  void unitcell_to_pdbfile(space_group_type* spaceg, atom_list_type* atom_list, const char* filename);

/*! \param[in] x0 ARRAY
 *
 *  \param[in] x1 ARRAY
 *
 *  \param[in] s0 ARRAY
 *
 *  \param[in] s1 ARRAY
*/
  void distance_and_sigma(const float* derm, const std::vector<float>* x0, const std::vector<float>* x1, const std::vector<float>* s0, const std::vector<float>* s1, float* dis, float* s);

  void genb(FortranMatrix<float> *b);

  void calc_magnetic_strf_miv(magsymm_k_type* mgp, matom_list_type* atm, magh_type* mh);

  void calc_magnetic_strf_miv_dom(magsymm_k_type* mgp, matom_list_type* atm, magnetic_domain_type* mag_dom, maghd_type* mh);

/*! \param[in] ord OPTIONAL
 *
 *  \param[in] powder OPTIONAL
 *
 *  \param[in] hkl OPTIONAL
*/
  void gen_satellites(magsymm_k_type* grp, float smax, magh_list_type* h, const int* ord=NULL, const int* powder=NULL, reflection_list_type* hkl=NULL);

  ADDRESS data_ptr;

private:
  bool initialized;
};

#endif /* CRYSTAL_CELL_TYPE_H_ */
