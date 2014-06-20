function getnum_k(obj_var)
	type (MagH_Type) :: obj_var
	integer :: getnum_k
	getnum_k = obj_var%num_k
end function getnum_k

subroutine setnum_k(obj_var, new_value)
	type (MagH_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_k = new_value
end subroutine setnum_k

function getH(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp), dimension(3) :: getH
	getH = obj_var%H
end function getH

subroutine setH(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine setH

function getsqMiV(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: getsqMiV
	getsqMiV = obj_var%sqMiV
end function getsqMiV

subroutine setsqMiV(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqMiV = new_value
end subroutine setsqMiV

function gets(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: gets
	gets = obj_var%s
end function gets

subroutine sets(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%s = new_value
end subroutine sets

function getMiVC(obj_var)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3) :: getMiVC
	getMiVC = obj_var%MiVC
end function getMiVC

subroutine setMiVC(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MiVC = new_value
end subroutine setMiVC

function getsignp(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: getsignp
	getsignp = obj_var%signp
end function getsignp

subroutine setsignp(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%signp = new_value
end subroutine setsignp

function getkeqv_minus(obj_var)
	type (MagH_Type) :: obj_var
	logical :: getkeqv_minus
	getkeqv_minus = obj_var%keqv_minus
end function getkeqv_minus

subroutine setkeqv_minus(obj_var, new_value)
	type (MagH_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%keqv_minus = new_value
end subroutine setkeqv_minus

function getMsF(obj_var)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3) :: getMsF
	getMsF = obj_var%MsF
end function getMsF

subroutine setMsF(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MsF = new_value
end subroutine setMsF

function getmult(obj_var)
	type (MagH_Type) :: obj_var
	integer :: getmult
	getmult = obj_var%mult
end function getmult

subroutine setmult(obj_var, new_value)
	type (MagH_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%mult = new_value
end subroutine setmult

function getMiV(obj_var)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3) :: getMiV
	getMiV = obj_var%MiV
end function getMiV

subroutine setMiV(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MiV = new_value
end subroutine setMiV

subroutine MagH_Type_ctor(MagH_Type_param, num_k_param, H_param, sqMiV_param, s_param, MiVC_param, signp_param, keqv_minus_param, MsF_param, mult_param, MiV_param)
	type (MagH_Type) :: MagH_Type_param
	integer, intent(in) :: num_k_param
	real(kind=cp), dimension(3), intent(in) :: H_param
	real(kind=cp), intent(in) :: sqMiV_param
	real(kind=cp), intent(in) :: s_param
	complex(kind=cp), dimension(3), intent(in) :: MiVC_param
	real(kind=cp), intent(in) :: signp_param
	logical, intent(in) :: keqv_minus_param
	complex(kind=cp), dimension(3), intent(in) :: MsF_param
	integer, intent(in) :: mult_param
	complex(kind=cp), dimension(3), intent(in) :: MiV_param
	MagH_Type_param%num_k = num_k_param
	MagH_Type_param%H = H_param
	MagH_Type_param%sqMiV = sqMiV_param
	MagH_Type_param%s = s_param
	MagH_Type_param%MiVC = MiVC_param
	MagH_Type_param%signp = signp_param
	MagH_Type_param%keqv_minus = keqv_minus_param
	MagH_Type_param%MsF = MsF_param
	MagH_Type_param%mult = mult_param
	MagH_Type_param%MiV = MiV_param
end subroutine MagH_Type_ctor
