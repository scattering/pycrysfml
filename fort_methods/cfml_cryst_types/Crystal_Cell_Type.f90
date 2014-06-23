function get_crystal_cell_lang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3) :: getlang
	getlang = obj_var%lang
end function get_crystal_cell_lang

subroutine set_crystal_cell_lang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%lang = new_value
end subroutine set_crystal_cell_lang

function get_crystal_cell_RCellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: getRCellVol
	getRCellVol = obj_var%RCellVol
end function get_crystal_cell_RCellVol

subroutine set_crystal_cell_RCellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%RCellVol = new_value
end subroutine set_crystal_cell_RCellVol

function get_crystal_cell_cell_std(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getcell_std
	getcell_std = obj_var%cell_std
end function get_crystal_cell_cell_std

subroutine set_crystal_cell_cell_std(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell_std = new_value
end subroutine set_crystal_cell_cell_std

function get_crystal_cell_ang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getang
	getang = obj_var%ang
end function get_crystal_cell_ang

subroutine set_crystal_cell_ang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%ang = new_value
end subroutine set_crystal_cell_ang

function get_crystal_cell_BL_Minv(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getBL_Minv
	getBL_Minv = obj_var%BL_Minv
end function get_crystal_cell_BL_Minv

subroutine set_crystal_cell_BL_Minv(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_Minv = new_value
end subroutine set_crystal_cell_BL_Minv

function get_crystal_cell_GR(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getGR
	getGR = obj_var%GR
end function get_crystal_cell_GR

subroutine set_crystal_cell_GR(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%GR = new_value
end subroutine set_crystal_cell_GR

function get_crystal_cell_Cr_Orth_cel(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getCr_Orth_cel
	getCr_Orth_cel = obj_var%Cr_Orth_cel
end function get_crystal_cell_Cr_Orth_cel

subroutine set_crystal_cell_Cr_Orth_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Cr_Orth_cel = new_value
end subroutine set_crystal_cell_Cr_Orth_cel

function get_crystal_cell_BL_M(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getBL_M
	getBL_M = obj_var%BL_M
end function get_crystal_cell_BL_M

subroutine set_crystal_cell_BL_M(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_M = new_value
end subroutine set_crystal_cell_BL_M

function get_crystal_cell_Orth_Cr_cel(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getOrth_Cr_cel
	getOrth_Cr_cel = obj_var%Orth_Cr_cel
end function get_crystal_cell_Orth_Cr_cel

subroutine set_crystal_cell_Orth_Cr_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Orth_Cr_cel = new_value
end subroutine set_crystal_cell_Orth_Cr_cel

function get_crystal_cell_CartType(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1) :: getCartType
	getCartType = obj_var%CartType
end function get_crystal_cell_CartType

subroutine set_crystal_cell_CartType(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1), intent(in) :: new_value
	obj_var%CartType = new_value
end subroutine set_crystal_cell_CartType

function get_crystal_cell_rang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getrang
	getrang = obj_var%rang
end function get_crystal_cell_rang

subroutine set_crystal_cell_rang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%rang = new_value
end subroutine set_crystal_cell_rang

function get_crystal_cell_rcell(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getrcell
	getrcell = obj_var%rcell
end function get_crystal_cell_rcell

subroutine set_crystal_cell_rcell(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%rcell = new_value
end subroutine set_crystal_cell_rcell

function get_crystal_cell_cell(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getcell
	getcell = obj_var%cell
end function get_crystal_cell_cell

subroutine set_crystal_cell_cell(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell = new_value
end subroutine set_crystal_cell_cell

function get_crystal_cell_GD(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getGD
	getGD = obj_var%GD
end function get_crystal_cell_GD

subroutine set_crystal_cell_GD(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%GD = new_value
end subroutine set_crystal_cell_GD

function get_crystal_cell_CellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: getCellVol
	getCellVol = obj_var%CellVol
end function get_crystal_cell_CellVol

subroutine set_crystal_cell_CellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%CellVol = new_value
end subroutine set_crystal_cell_CellVol

function get_crystal_cell_ang_std(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getang_std
	getang_std = obj_var%ang_std
end function get_crystal_cell_ang_std

subroutine set_crystal_cell_ang_std(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%ang_std = new_value
end subroutine set_crystal_cell_ang_std

function get_crystal_cell_lcell(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3) :: getlcell
	getlcell = obj_var%lcell
end function get_crystal_cell_lcell

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
