Module cfml_python
USE CFML_IO_Formats
USE CFML_Reflections_Utilities
USE CFML_Magnetic_Symmetry
Use CFML_GlobalDeps
Use CFML_Crystal_Metrics
Use CFML_Crystallographic_Symmetry
Use CFML_Atom_TypeDef
!! reflection_type array wrapper !!
TYPE reflct_array_list
	INTEGER :: current_index
	type (Reflect_Type), dimension(100) :: reflections
END TYPE reflct_array_list
CONTAINS
!!-- Wrappers for procedures called from ctypes hklGen.py --!!

!-- list methods for reflection array wrapper
SUBROUTINE reflct_array_ctor(array)
	type (reflct_array_list) :: array
	array%current_index = 1
END SUBROUTINE reflct_array_ctor
SUBROUTINE reflct_append(array, rflctn)
	type (reflct_array_list) :: array
	type (Reflect_Type) :: rflctn
	array%reflections(array%current_index) = rflctn
	array%current_index = array%current_index+1
END SUBROUTINE reflct_append
! wrapper for procedure in CFML_IO_Formats
SUBROUTINE ReadXtal_Structure_File(filenam, Cell, SpG, A, Mode, Iphase, Job_Info, file_list, CFrame)
!---- Arguments ----!
       character(len=*), intent(in) :: filenam
       Type (Crystal_Cell_Type), intent(out) :: Cell
       Type (Space_Group_Type), intent(out) :: SpG
       Type (atom_list_type), intent(out) :: A
       Character(len=*), optional, intent(in) :: Mode
       Integer, optional, intent(in) :: Iphase
       Type(Job_Info_type), optional, intent(out) :: Job_Info
       Type(file_list_type), optional, intent(in out) :: file_list
       Character(len=*), optional, intent(in) :: CFrame
       call Readn_Set_Xtal_Structure(filenam, Cell, SpG, A, Mode, Iphase, Job_Info, file_list, CFrame)
END SUBROUTINE ReadXtal_Structure_File

! wrapper for Hkl_S --> HS_R procedure from CFML_Reflect_Util
FUNCTION HklS_R(H, Crystalcell)
	real(kind=cp), dimension(3), intent(in) :: h
	type (Crystal_Cell_Type), intent(in) :: CrystalCell
	real(kind=cp) :: HklS_R
	HklS_R = Hkl_S(H, Crystalcell)
END FUNCTION HklS_R

! wrappers for HKL_GEN_SXTAL from CFML_Reflect_Util
! _reflection
SUBROUTINE hklgen_sxtal_reflection(Crystalcell, Spacegroup, stlmin, stlmax, Num_Ref, Reflex, ord, hlim)
	type (Crystal_Cell_Type), intent(in) :: crystalcell
	type (Space_Group_Type), intent(in) :: spacegroup
	real(kind=cp), intent(in) :: stlmin,stlmax
	integer, intent(out) :: num_ref
	type (reflct_array_list), intent(out) :: reflex
	Integer, dimension(3), optional, intent(in) :: ord
	Integer, dimension(3,2), optional, intent(in) :: hlim
	call HKL_GEN_SXTAL(Crystalcell, Spacegroup, stlmin, stlmax, Num_Ref, Reflex%reflections, ord, hlim)
END SUBROUTINE hklgen_sxtal_reflection

! _list
SUBROUTINE hklgen_sxtal_list(Crystalcell, Spacegroup, stlmin, stlmax, Num_Ref, Reflex, ord, hlim)
	type (Crystal_Cell_Type), intent(in) :: crystalcell
	type (Space_Group_Type), intent(in) :: spacegroup
	real(kind=cp), intent(in) :: stlmin,stlmax
	integer, intent(out) :: num_ref
	Type(Reflection_List_Type), intent(out) :: reflex
	Integer, dimension(3), optional, intent(in) :: ord
	Integer, dimension(3,2), optional, intent(in) :: hlim
	call HKL_GEN_SXTAL(Crystalcell, Spacegroup, stlmin, stlmax, Num_Ref, Reflex, ord, hlim)
END SUBROUTINE hklgen_sxtal_list

! wrappers for Hkl_Uni from CFML_Reflect_Util
! _reflection
SUBROUTINE hkluni_reflection(Crystalcell, Spacegroup, Friedel, Value1, Value2, Code, Num_Ref, Reflex, no_order)
	!---- Arguments ----!
	type (Crystal_Cell_Type), intent(in) :: crystalcell
	type (Space_Group_Type), intent(in) :: spacegroup
	Logical, intent(in) :: Friedel
	real(kind=cp), intent(in) :: value1,value2
	character(len=1), intent(in) :: code
	integer, intent(out) :: num_ref
	type (reflct_array_list), intent(out) :: reflex
	logical, optional, intent(in) :: no_order
	call Hkl_Uni(Crystalcell, Spacegroup, Friedel, Value1, Value2, Code, Num_Ref, Reflex%reflections, no_order)
END SUBROUTINE hkluni_reflection

! _ReflList
SUBROUTINE hkluni_refllist(Crystalcell, Spacegroup, Friedel, Value1, Value2, Code, Num_Ref, Reflex, no_order)
	!---- Arguments ----!
	type (Crystal_Cell_Type), intent(in) :: crystalcell
	type (Space_Group_Type), intent(in) :: spacegroup
	Logical, intent(in) :: Friedel
	real(kind=cp), intent(in) :: value1,value2
	character(len=1), intent(in) :: code
	integer, intent(in) :: Num_Ref
	type (Reflection_List_Type), intent(out) :: reflex
	logical, optional, intent(in) :: no_order
	call Hkl_Uni(Crystalcell, Spacegroup, Friedel, Value1, Value2, Code, Num_Ref, Reflex, no_order)
END SUBROUTINE hkluni_refllist

! wrappers for Readn_Set_Magnetic_Structure from CFML_Magnetic_Symmetry
! _CFL
SUBROUTINE read_mag_cfl_file(file_cfl, n_ini, n_end, MGp, Am, SGo, Mag_dom, Cell)
	!---- Arguments ----!
	type(file_list_type), intent(in) :: file_cfl
	integer, intent(in out) :: n_ini, n_end
	type(MagSymm_k_Type), intent(out) :: MGp
	type(mAtom_List_Type), intent(out) :: Am
	type(Magnetic_Group_Type), optional, intent(out) :: SGo
	type(Magnetic_Domain_type),optional, intent(out) :: Mag_dom
	type(Crystal_Cell_type), optional, intent(in) :: Cell
	call Readn_Set_Magnetic_Structure(file_cfl, n_ini, n_end, MGp, Am, SGo, Mag_dom, Cell)
END SUBROUTINE read_mag_cfl_file

! _MCIF
SUBROUTINE read_mag_mcif_file(file_mcif, mCell, MGp, Am)
	character(len=*),               intent (in)  :: file_mcif
	type(Crystal_Cell_type),        intent (out) :: mCell
	type(Magnetic_Space_Group_Type),intent (out) :: MGp
	type(mAtom_List_Type),          intent (out) :: Am
	call Readn_Set_Magnetic_Structure(file_mcif, mCell, MGp, Am)
END SUBROUTINE read_mag_mcif_file

!!-- END Procedures called from ctypes hklGen.py --!!
End Module cfml_python
