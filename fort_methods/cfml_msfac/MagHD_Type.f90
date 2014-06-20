function getsqAMiV(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: getsqAMiV
	getsqAMiV = obj_var%sqAMiV
end function getsqAMiV

subroutine setsqAMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqAMiV = new_value
end subroutine setsqAMiV

function getnum_k(obj_var)
	type (MagHD_Type) :: obj_var
	integer :: getnum_k
	getnum_k = obj_var%num_k
end function getnum_k

subroutine setnum_k(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_k = new_value
end subroutine setnum_k

function getH(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp),   dimension(3) :: getH
	getH = obj_var%H
end function getH

subroutine setH(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp),   dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine setH

function getsqMiV(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: getsqMiV
	getsqMiV = obj_var%sqMiV
end function getsqMiV

subroutine setsqMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqMiV = new_value
end subroutine setsqMiV

function gets(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: gets
	gets = obj_var%s
end function gets

subroutine sets(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%s = new_value
end subroutine sets

function getMiVC(obj_var)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24) :: getMiVC
	getMiVC = obj_var%MiVC
end function getMiVC

subroutine setMiVC(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MiVC = new_value
end subroutine setMiVC

function getsignp(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: getsignp
	getsignp = obj_var%signp
end function getsignp

subroutine setsignp(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%signp = new_value
end subroutine setsignp

function getMsF(obj_var)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24) :: getMsF
	getMsF = obj_var%MsF
end function getMsF

subroutine setMsF(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MsF = new_value
end subroutine setMsF

function getAMiV(obj_var)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3) :: getAMiV
	getAMiV = obj_var%AMiV
end function getAMiV

subroutine setAMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%AMiV = new_value
end subroutine setAMiV

function getkeqv_minus(obj_var)
	type (MagHD_Type) :: obj_var
	logical :: getkeqv_minus
	getkeqv_minus = obj_var%keqv_minus
end function getkeqv_minus

subroutine setkeqv_minus(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%keqv_minus = new_value
end subroutine setkeqv_minus

function getMiV(obj_var)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24) :: getMiV
	getMiV = obj_var%MiV
end function getMiV

subroutine setMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MiV = new_value
end subroutine setMiV

subroutine MagHD_Type_ctor(MagHD_Type_param, sqAMiV_param, num_k_param, H_param, sqMiV_param, s_param, MiVC_param, signp_param, MsF_param, AMiV_param, keqv_minus_param, MiV_param)
	type (MagHD_Type) :: MagHD_Type_param
	real(kind=cp), intent(in) :: sqAMiV_param
	integer, intent(in) :: num_k_param
	real(kind=cp),   dimension(3), intent(in) :: H_param
	real(kind=cp), intent(in) :: sqMiV_param
	real(kind=cp), intent(in) :: s_param
	complex(kind=cp),dimension(3,2,24), intent(in) :: MiVC_param
	real(kind=cp), intent(in) :: signp_param
	complex(kind=cp),dimension(3,2,24), intent(in) :: MsF_param
	complex(kind=cp),dimension(3), intent(in) :: AMiV_param
	logical, intent(in) :: keqv_minus_param
	complex(kind=cp),dimension(3,2,24), intent(in) :: MiV_param
	MagHD_Type_param%sqAMiV = sqAMiV_param
	MagHD_Type_param%num_k = num_k_param
	MagHD_Type_param%H = H_param
	MagHD_Type_param%sqMiV = sqMiV_param
	MagHD_Type_param%s = s_param
	MagHD_Type_param%MiVC = MiVC_param
	MagHD_Type_param%signp = signp_param
	MagHD_Type_param%MsF = MsF_param
	MagHD_Type_param%AMiV = AMiV_param
	MagHD_Type_param%keqv_minus = keqv_minus_param
	MagHD_Type_param%MiV = MiV_param
end subroutine MagHD_Type_ctor
