function get_magh_num_k(obj_var)
	type (MagH_Type) :: obj_var
	integer :: get_magh_num_k
	get_magh_num_k = obj_var%num_k
end function get_magh_num_k

subroutine set_magh_num_k(obj_var, new_value)
	type (MagH_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_k = new_value
end subroutine set_magh_num_k

subroutine get_magh_H(obj_var, output_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%H
end subroutine get_magh_H

subroutine set_magh_H(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine set_magh_H

function get_magh_sqMiV(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: get_magh_sqMiV
	get_magh_sqMiV = obj_var%sqMiV
end function get_magh_sqMiV

subroutine set_magh_sqMiV(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqMiV = new_value
end subroutine set_magh_sqMiV

function get_magh_s(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: get_magh_s
	get_magh_s = obj_var%s
end function get_magh_s

subroutine set_magh_s(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%s = new_value
end subroutine set_magh_s

subroutine get_magh_MiVC(obj_var, output_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%MiVC
end subroutine get_magh_MiVC

subroutine set_magh_MiVC(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MiVC = new_value
end subroutine set_magh_MiVC

function get_magh_signp(obj_var)
	type (MagH_Type) :: obj_var
	real(kind=cp) :: get_magh_signp
	get_magh_signp = obj_var%signp
end function get_magh_signp

subroutine set_magh_signp(obj_var, new_value)
	type (MagH_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%signp = new_value
end subroutine set_magh_signp

function get_magh_keqv_minus(obj_var)
	type (MagH_Type) :: obj_var
	logical :: get_magh_keqv_minus
	get_magh_keqv_minus = obj_var%keqv_minus
end function get_magh_keqv_minus

subroutine set_magh_keqv_minus(obj_var, new_value)
	type (MagH_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%keqv_minus = new_value
end subroutine set_magh_keqv_minus

subroutine get_magh_MsF(obj_var, output_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%MsF
end subroutine get_magh_MsF

subroutine set_magh_MsF(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MsF = new_value
end subroutine set_magh_MsF

function get_magh_mult(obj_var)
	type (MagH_Type) :: obj_var
	integer :: get_magh_mult
	get_magh_mult = obj_var%mult
end function get_magh_mult

subroutine set_magh_mult(obj_var, new_value)
	type (MagH_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%mult = new_value
end subroutine set_magh_mult

subroutine get_magh_MiV(obj_var, output_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%MiV
end subroutine get_magh_MiV

subroutine set_magh_MiV(obj_var, new_value)
	type (MagH_Type) :: obj_var
	complex(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MiV = new_value
end subroutine set_magh_MiV

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
