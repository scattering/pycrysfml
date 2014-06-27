function get_interval_maxb(obj_var)
	type (interval_type) :: obj_var
	real(kind=cp) :: get_interval_maxb
	get_interval_maxb = obj_var%maxb
end function get_interval_maxb

subroutine set_interval_maxb(obj_var, new_value)
	type (interval_type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%maxb = new_value
end subroutine set_interval_maxb

function get_interval_mina(obj_var)
	type (interval_type) :: obj_var
	real(kind=cp) :: get_interval_mina
	get_interval_mina = obj_var%mina
end function get_interval_mina

subroutine set_interval_mina(obj_var, new_value)
	type (interval_type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%mina = new_value
end subroutine set_interval_mina

subroutine interval_type_ctor(interval_type_param, maxb_param, mina_param)
	type (interval_type) :: interval_type_param
	real(kind=cp), intent(in) :: maxb_param
	real(kind=cp), intent(in) :: mina_param
	interval_type_param%maxb = maxb_param
	interval_type_param%mina = mina_param
end subroutine interval_type_ctor
