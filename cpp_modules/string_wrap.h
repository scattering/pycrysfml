#ifndef STRING_WRAP
#define STRING_WRAP
#include <iostream>
#include <string>
#include "powder_numor_type.h"
#include "crystal_cell_type.h"
#include "atom_equiv_type.h"
#include "magnetic_space_group_type.h"
#include "job_info_type.h"
#include "magsymm_k_type.h"
#include "matom_type.h"
#include "wyck_pos_type.h"
#include "atom_type.h"
#include "magnetic_group_type.h"
#include "diffraction_pattern_type.h"
#include "ns_space_group_type.h"
#include "space_group_type.h"
std::string getPowderNumor_title(powder_numor_type* obj);
std::string getPowderNumor_instrm(powder_numor_type* obj);
std::string getPowderNumor_header(powder_numor_type* obj);
std::string getPowderNumor_scantype(powder_numor_type* obj);
std::string getCrystalCell_carttype(crystal_cell_type* obj);
std::string getAtomEquiv_chemsymb(atom_equiv_type* obj);
std::string getMagneticSpaceGroup_bns_symbol(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_bns_number(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_crystalsys(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_spg_lat(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_og_number(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_spg_latsy(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_parent_spg(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_og_symbol(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_centre(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_trn_to_standard(magnetic_space_group_type* obj);
std::string getMagneticSpaceGroup_trn_from_parent(magnetic_space_group_type* obj);
std::string getJobInfo_title(job_info_type* obj);
std::string getMagsymmK_latt(magsymm_k_type* obj);
std::string getMagsymmK_bns_symbol(magsymm_k_type* obj);
std::string getMagsymmK_bns_number(magsymm_k_type* obj);
std::string getMagsymmK_magmodel(magsymm_k_type* obj);
std::string getMagsymmK_sk_type(magsymm_k_type* obj);
std::string getMagsymmK_og_number(magsymm_k_type* obj);
std::string getMagsymmK_parent_spg(magsymm_k_type* obj);
std::string getMagsymmK_og_symbol(magsymm_k_type* obj);
std::string getMatom_utype(matom_type* obj);
std::string getMatom_sfacsymb(matom_type* obj);
std::string getMatom_lab(matom_type* obj);
std::string getMatom_wyck(matom_type* obj);
std::string getMatom_chemsymb(matom_type* obj);
std::string getMatom_thtype(matom_type* obj);
std::string getMatom_atminfo(matom_type* obj);
std::string getWyckPos_str_orig(wyck_pos_type* obj);
std::string getWyckPos_site(wyck_pos_type* obj);
std::string getAtom_utype(atom_type* obj);
std::string getAtom_sfacsymb(atom_type* obj);
std::string getAtom_lab(atom_type* obj);
std::string getAtom_wyck(atom_type* obj);
std::string getAtom_chemsymb(atom_type* obj);
std::string getAtom_thtype(atom_type* obj);
std::string getAtom_atminfo(atom_type* obj);
std::string getMagneticGroup_shubnikov(magnetic_group_type* obj);
std::string getDiffractionPattern_instr(diffraction_pattern_type* obj);
std::string getDiffractionPattern_yax_text(diffraction_pattern_type* obj);
std::string getDiffractionPattern_diff_kind(diffraction_pattern_type* obj);
std::string getDiffractionPattern_filepath(diffraction_pattern_type* obj);
std::string getDiffractionPattern_title(diffraction_pattern_type* obj);
std::string getDiffractionPattern_filename(diffraction_pattern_type* obj);
std::string getDiffractionPattern_scat_var(diffraction_pattern_type* obj);
std::string getDiffractionPattern_xax_text(diffraction_pattern_type* obj);
std::string getNsSpaceGroup_crystalsys(ns_space_group_type* obj);
std::string getNsSpaceGroup_pg(ns_space_group_type* obj);
std::string getNsSpaceGroup_hall(ns_space_group_type* obj);
std::string getNsSpaceGroup_info(ns_space_group_type* obj);
std::string getNsSpaceGroup_spg_lat(ns_space_group_type* obj);
std::string getNsSpaceGroup_laue(ns_space_group_type* obj);
std::string getNsSpaceGroup_spg_latsy(ns_space_group_type* obj);
std::string getNsSpaceGroup_bravais(ns_space_group_type* obj);
std::string getNsSpaceGroup_sg_setting(ns_space_group_type* obj);
std::string getNsSpaceGroup_ghall(ns_space_group_type* obj);
std::string getNsSpaceGroup_spg_symb(ns_space_group_type* obj);
std::string getNsSpaceGroup_centre(ns_space_group_type* obj);
std::string getSpaceGroup_crystalsys(space_group_type* obj);
std::string getSpaceGroup_pg(space_group_type* obj);
std::string getSpaceGroup_hall(space_group_type* obj);
std::string getSpaceGroup_info(space_group_type* obj);
std::string getSpaceGroup_spg_lat(space_group_type* obj);
std::string getSpaceGroup_laue(space_group_type* obj);
std::string getSpaceGroup_spg_latsy(space_group_type* obj);
std::string getSpaceGroup_bravais(space_group_type* obj);
std::string getSpaceGroup_sg_setting(space_group_type* obj);
std::string getSpaceGroup_ghall(space_group_type* obj);
std::string getSpaceGroup_spg_symb(space_group_type* obj);
std::string getSpaceGroup_centre(space_group_type* obj);
#endif
