subroutine get_crystal_cell_lang(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(out) :: output_value
	output_value = obj_var%lang
end subroutine get_crystal_cell_lang

subroutine set_crystal_cell_lang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%lang = new_value
end subroutine set_crystal_cell_lang

function get_crystal_cell_RCellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: get_crystal_cell_RCellVol
	get_crystal_cell_RCellVol = obj_var%RCellVol
end function get_crystal_cell_RCellVol

subroutine set_crystal_cell_RCellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%RCellVol = new_value
end subroutine set_crystal_cell_RCellVol

subroutine get_crystal_cell_cell_std(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%cell_std
end subroutine get_crystal_cell_cell_std

subroutine set_crystal_cell_cell_std(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell_std = new_value
end subroutine set_crystal_cell_cell_std

subroutine get_crystal_cell_ang(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%ang
end subroutine get_crystal_cell_ang

subroutine set_crystal_cell_ang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%ang = new_value
end subroutine set_crystal_cell_ang

subroutine get_crystal_cell_BL_Minv(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%BL_Minv
end subroutine get_crystal_cell_BL_Minv

subroutine set_crystal_cell_BL_Minv(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_Minv = new_value
end subroutine set_crystal_cell_BL_Minv

subroutine get_crystal_cell_GR(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%GR
end subroutine get_crystal_cell_GR

subroutine set_crystal_cell_GR(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%GR = new_value
end subroutine set_crystal_cell_GR

subroutine get_crystal_cell_Cr_Orth_cel(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%Cr_Orth_cel
end subroutine get_crystal_cell_Cr_Orth_cel

subroutine set_crystal_cell_Cr_Orth_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Cr_Orth_cel = new_value
end subroutine set_crystal_cell_Cr_Orth_cel

subroutine get_crystal_cell_BL_M(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%BL_M
end subroutine get_crystal_cell_BL_M

subroutine set_crystal_cell_BL_M(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_M = new_value
end subroutine set_crystal_cell_BL_M

subroutine get_crystal_cell_Orth_Cr_cel(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%Orth_Cr_cel
end subroutine get_crystal_cell_Orth_Cr_cel

subroutine set_crystal_cell_Orth_Cr_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Orth_Cr_cel = new_value
end subroutine set_crystal_cell_Orth_Cr_cel

subroutine get_crystal_cell_CartType(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1), intent(out) :: output_value
	output_value = obj_var%CartType
end subroutine get_crystal_cell_CartType

subroutine set_crystal_cell_CartType(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1), intent(in) :: new_value
	obj_var%CartType = new_value
end subroutine set_crystal_cell_CartType

subroutine get_crystal_cell_rang(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%rang
end subroutine get_crystal_cell_rang

subroutine set_crystal_cell_rang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%rang = new_value
end subroutine set_crystal_cell_rang

subroutine get_crystal_cell_rcell(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%rcell
end subroutine get_crystal_cell_rcell

subroutine set_crystal_cell_rcell(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%rcell = new_value
end subroutine set_crystal_cell_rcell

subroutine get_crystal_cell_cell(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%cell
end subroutine get_crystal_cell_cell

subroutine set_crystal_cell_cell(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell = new_value
end subroutine set_crystal_cell_cell

subroutine get_crystal_cell_GD(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(out) :: output_value
	output_value = obj_var%GD
end subroutine get_crystal_cell_GD

subroutine set_crystal_cell_GD(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%GD = new_value
end subroutine set_crystal_cell_GD

function get_crystal_cell_CellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: get_crystal_cell_CellVol
	get_crystal_cell_CellVol = obj_var%CellVol
end function get_crystal_cell_CellVol

subroutine set_crystal_cell_CellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%CellVol = new_value
end subroutine set_crystal_cell_CellVol

subroutine get_crystal_cell_ang_std(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%ang_std
end subroutine get_crystal_cell_ang_std

subroutine set_crystal_cell_ang_std(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%ang_std = new_value
end subroutine set_crystal_cell_ang_std

subroutine get_crystal_cell_lcell(obj_var, output_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(out) :: output_value
	output_value = obj_var%lcell
end subroutine get_crystal_cell_lcell

subroutine set_crystal_cell_lcell(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%lcell = new_value
end subroutine set_crystal_cell_lcell

subroutine Crystal_Cell_Type_ctor(Crystal_Cell_Type_param, lang_param, RCellVol_param, cell_std_param, ang_param, BL_Minv_param, GR_param, Cr_Orth_cel_param, BL_M_param, Orth_Cr_cel_param, CartType_param, rang_param, rcell_param, cell_param, GD_param, CellVol_param, ang_std_param, lcell_param)
	type (Crystal_Cell_Type) :: Crystal_Cell_Type_param
	integer,      dimension(3), intent(in) :: lang_param
	real(kind=cp), intent(in) :: RCellVol_param
	real(kind=cp),dimension(3), intent(in) :: cell_std_param
	real(kind=cp),dimension(3), intent(in) :: ang_param
	real(kind=cp),dimension(3,3), intent(in) :: BL_Minv_param
	real(kind=cp),dimension(3,3), intent(in) :: GR_param
	real(kind=cp),dimension(3,3), intent(in) :: Cr_Orth_cel_param
	real(kind=cp),dimension(3,3), intent(in) :: BL_M_param
	real(kind=cp),dimension(3,3), intent(in) :: Orth_Cr_cel_param
	character (len=1), intent(in) :: CartType_param
	real(kind=cp),dimension(3), intent(in) :: rang_param
	real(kind=cp),dimension(3), intent(in) :: rcell_param
	real(kind=cp),dimension(3), intent(in) :: cell_param
	real(kind=cp),dimension(3,3), intent(in) :: GD_param
	real(kind=cp), intent(in) :: CellVol_param
	real(kind=cp),dimension(3), intent(in) :: ang_std_param
	integer,      dimension(3), intent(in) :: lcell_param
	Crystal_Cell_Type_param%lang = lang_param
	Crystal_Cell_Type_param%RCellVol = RCellVol_param
	Crystal_Cell_Type_param%cell_std = cell_std_param
	Crystal_Cell_Type_param%ang = ang_param
	Crystal_Cell_Type_param%BL_Minv = BL_Minv_param
	Crystal_Cell_Type_param%GR = GR_param
	Crystal_Cell_Type_param%Cr_Orth_cel = Cr_Orth_cel_param
	Crystal_Cell_Type_param%BL_M = BL_M_param
	Crystal_Cell_Type_param%Orth_Cr_cel = Orth_Cr_cel_param
	Crystal_Cell_Type_param%CartType = CartType_param
	Crystal_Cell_Type_param%rang = rang_param
	Crystal_Cell_Type_param%rcell = rcell_param
	Crystal_Cell_Type_param%cell = cell_param
	Crystal_Cell_Type_param%GD = GD_param
	Crystal_Cell_Type_param%CellVol = CellVol_param
	Crystal_Cell_Type_param%ang_std = ang_std_param
	Crystal_Cell_Type_param%lcell = lcell_param
end subroutine Crystal_Cell_Type_ctor
