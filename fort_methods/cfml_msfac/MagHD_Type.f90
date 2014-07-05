function get_maghd_sqAMiV(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: get_maghd_sqAMiV
	get_maghd_sqAMiV = obj_var%sqAMiV
end function get_maghd_sqAMiV

subroutine set_maghd_sqAMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqAMiV = new_value
end subroutine set_maghd_sqAMiV

function get_maghd_num_k(obj_var)
	type (MagHD_Type) :: obj_var
	integer :: get_maghd_num_k
	get_maghd_num_k = obj_var%num_k
end function get_maghd_num_k

subroutine set_maghd_num_k(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_k = new_value
end subroutine set_maghd_num_k

subroutine get_maghd_H(obj_var, output_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp),   dimension(3), intent(out) :: output_value
	output_value = obj_var%H
end subroutine get_maghd_H

subroutine set_maghd_H(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp),   dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine set_maghd_H

function get_maghd_sqMiV(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: get_maghd_sqMiV
	get_maghd_sqMiV = obj_var%sqMiV
end function get_maghd_sqMiV

subroutine set_maghd_sqMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%sqMiV = new_value
end subroutine set_maghd_sqMiV

function get_maghd_s(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: get_maghd_s
	get_maghd_s = obj_var%s
end function get_maghd_s

subroutine set_maghd_s(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%s = new_value
end subroutine set_maghd_s

subroutine get_maghd_MiVC(obj_var, output_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(out) :: output_value
	output_value = obj_var%MiVC
end subroutine get_maghd_MiVC

subroutine set_maghd_MiVC(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MiVC = new_value
end subroutine set_maghd_MiVC

function get_maghd_signp(obj_var)
	type (MagHD_Type) :: obj_var
	real(kind=cp) :: get_maghd_signp
	get_maghd_signp = obj_var%signp
end function get_maghd_signp

subroutine set_maghd_signp(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%signp = new_value
end subroutine set_maghd_signp

subroutine get_maghd_MsF(obj_var, output_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(out) :: output_value
	output_value = obj_var%MsF
end subroutine get_maghd_MsF

subroutine set_maghd_MsF(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MsF = new_value
end subroutine set_maghd_MsF

subroutine get_maghd_AMiV(obj_var, output_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%AMiV
end subroutine get_maghd_AMiV

subroutine set_maghd_AMiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%AMiV = new_value
end subroutine set_maghd_AMiV

function get_maghd_keqv_minus(obj_var)
	type (MagHD_Type) :: obj_var
	logical :: get_maghd_keqv_minus
	get_maghd_keqv_minus = obj_var%keqv_minus
end function get_maghd_keqv_minus

subroutine set_maghd_keqv_minus(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%keqv_minus = new_value
end subroutine set_maghd_keqv_minus

subroutine get_maghd_MiV(obj_var, output_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(out) :: output_value
	output_value = obj_var%MiV
end subroutine get_maghd_MiV

subroutine set_maghd_MiV(obj_var, new_value)
	type (MagHD_Type) :: obj_var
	complex(kind=cp),dimension(3,2,24), intent(in) :: new_value
	obj_var%MiV = new_value
end subroutine set_maghd_MiV

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
