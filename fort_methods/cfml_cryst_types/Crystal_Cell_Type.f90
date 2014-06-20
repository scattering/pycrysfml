function getRCellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: getRCellVol
	getRCellVol = obj_var%RCellVol
end function getRCellVol

subroutine setRCellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%RCellVol = new_value
end subroutine setRCellVol

function getCartType(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1) :: getCartType
	getCartType = obj_var%CartType
end function getCartType

subroutine setCartType(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	character (len=1), intent(in) :: new_value
	obj_var%CartType = new_value
end subroutine setCartType

function getlcell, lang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3) :: getlcell, lang
	getlcell, lang = obj_var%lcell, lang
end function getlcell, lang

subroutine setlcell, lang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%lcell, lang = new_value
end subroutine setlcell, lang

function getBL_Minv(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getBL_Minv
	getBL_Minv = obj_var%BL_Minv
end function getBL_Minv

subroutine setBL_Minv(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_Minv = new_value
end subroutine setBL_Minv

function getCr_Orth_cel(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getCr_Orth_cel
	getCr_Orth_cel = obj_var%Cr_Orth_cel
end function getCr_Orth_cel

subroutine setCr_Orth_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Cr_Orth_cel = new_value
end subroutine setCr_Orth_cel

function getBL_M(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getBL_M
	getBL_M = obj_var%BL_M
end function getBL_M

subroutine setBL_M(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%BL_M = new_value
end subroutine setBL_M

function getOrth_Cr_cel(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getOrth_Cr_cel
	getOrth_Cr_cel = obj_var%Orth_Cr_cel
end function getOrth_Cr_cel

subroutine setOrth_Cr_cel(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%Orth_Cr_cel = new_value
end subroutine setOrth_Cr_cel

function getcell_std, ang_std(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getcell_std, ang_std
	getcell_std, ang_std = obj_var%cell_std, ang_std
end function getcell_std, ang_std

subroutine setcell_std, ang_std(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell_std, ang_std = new_value
end subroutine setcell_std, ang_std

function getGD,GR(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3) :: getGD,GR
	getGD,GR = obj_var%GD,GR
end function getGD,GR

subroutine setGD,GR(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3,3), intent(in) :: new_value
	obj_var%GD,GR = new_value
end subroutine setGD,GR

function getCellVol(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp) :: getCellVol
	getCellVol = obj_var%CellVol
end function getCellVol

subroutine setCellVol(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%CellVol = new_value
end subroutine setCellVol

function getcell, ang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getcell, ang
	getcell, ang = obj_var%cell, ang
end function getcell, ang

subroutine setcell, ang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%cell, ang = new_value
end subroutine setcell, ang

function getrcell, rang(obj_var)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3) :: getrcell, rang
	getrcell, rang = obj_var%rcell, rang
end function getrcell, rang

subroutine setrcell, rang(obj_var, new_value)
	type (Crystal_Cell_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%rcell, rang = new_value
end subroutine setrcell, rang

subroutine Crystal_Cell_Type_ctor(Crystal_Cell_Type_param, RCellVol_param, CartType_param, lcell, lang_param, BL_Minv_param, Cr_Orth_cel_param, BL_M_param, Orth_Cr_cel_param, cell_std, ang_std_param, GD,GR_param, CellVol_param, cell, ang_param, rcell, rang_param)
	type (Crystal_Cell_Type) :: Crystal_Cell_Type_param
	real(kind=cp), intent(in) :: RCellVol_param
	character (len=1), intent(in) :: CartType_param
	integer,      dimension(3), intent(in) :: lcell, lang_param
	real(kind=cp),dimension(3,3), intent(in) :: BL_Minv_param
	real(kind=cp),dimension(3,3), intent(in) :: Cr_Orth_cel_param
	real(kind=cp),dimension(3,3), intent(in) :: BL_M_param
	real(kind=cp),dimension(3,3), intent(in) :: Orth_Cr_cel_param
	real(kind=cp),dimension(3), intent(in) :: cell_std, ang_std_param
	real(kind=cp),dimension(3,3), intent(in) :: GD,GR_param
	real(kind=cp), intent(in) :: CellVol_param
	real(kind=cp),dimension(3), intent(in) :: cell, ang_param
	real(kind=cp),dimension(3), intent(in) :: rcell, rang_param
	Crystal_Cell_Type_param%RCellVol = RCellVol_param
	Crystal_Cell_Type_param%CartType = CartType_param
	Crystal_Cell_Type_param%lcell, lang = lcell, lang_param
	Crystal_Cell_Type_param%BL_Minv = BL_Minv_param
	Crystal_Cell_Type_param%Cr_Orth_cel = Cr_Orth_cel_param
	Crystal_Cell_Type_param%BL_M = BL_M_param
	Crystal_Cell_Type_param%Orth_Cr_cel = Orth_Cr_cel_param
	Crystal_Cell_Type_param%cell_std, ang_std = cell_std, ang_std_param
	Crystal_Cell_Type_param%GD,GR = GD,GR_param
	Crystal_Cell_Type_param%CellVol = CellVol_param
	Crystal_Cell_Type_param%cell, ang = cell, ang_param
	Crystal_Cell_Type_param%rcell, rang = rcell, rang_param
end subroutine Crystal_Cell_Type_ctor
