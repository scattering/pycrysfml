function getmaxb(obj_var)
	type (interval_type) :: obj_var
	real(kind=cp) :: getmaxb
	getmaxb = obj_var%maxb
end function getmaxb

subroutine setmaxb(obj_var, new_value)
	type (interval_type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%maxb = new_value
end subroutine setmaxb

function getmina(obj_var)
	type (interval_type) :: obj_var
	real(kind=cp) :: getmina
	getmina = obj_var%mina
end function getmina

subroutine setmina(obj_var, new_value)
	type (interval_type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%mina = new_value
end subroutine setmina

subroutine interval_type_ctor(interval_type_param, maxb_param, mina_param)
	type (interval_type) :: interval_type_param
	real(kind=cp), intent(in) :: maxb_param
	real(kind=cp), intent(in) :: mina_param
	interval_type_param%maxb = maxb_param
	interval_type_param%mina = mina_param
end subroutine interval_type_ctor
