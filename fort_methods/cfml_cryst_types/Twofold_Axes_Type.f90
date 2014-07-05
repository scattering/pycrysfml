subroutine get_twofold_axes_a(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%a
end subroutine get_twofold_axes_a

subroutine set_twofold_axes_a(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%a = new_value
end subroutine set_twofold_axes_a

subroutine get_twofold_axes_c(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%c
end subroutine get_twofold_axes_c

subroutine set_twofold_axes_c(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%c = new_value
end subroutine set_twofold_axes_c

subroutine get_twofold_axes_b(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%b
end subroutine get_twofold_axes_b

subroutine set_twofold_axes_b(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%b = new_value
end subroutine set_twofold_axes_b

subroutine get_twofold_axes_maxes(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(out) :: output_value
	output_value = obj_var%maxes
end subroutine get_twofold_axes_maxes

subroutine set_twofold_axes_maxes(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(in) :: new_value
	obj_var%maxes = new_value
end subroutine set_twofold_axes_maxes

function get_twofold_axes_ntwo(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	integer :: get_twofold_axes_ntwo
	get_twofold_axes_ntwo = obj_var%ntwo
end function get_twofold_axes_ntwo

subroutine set_twofold_axes_ntwo(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%ntwo = new_value
end subroutine set_twofold_axes_ntwo

subroutine get_twofold_axes_cross(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(out) :: output_value
	output_value = obj_var%cross
end subroutine get_twofold_axes_cross

subroutine set_twofold_axes_cross(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(in) :: new_value
	obj_var%cross = new_value
end subroutine set_twofold_axes_cross

function get_twofold_axes_tol(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) :: get_twofold_axes_tol
	get_twofold_axes_tol = obj_var%tol
end function get_twofold_axes_tol

subroutine set_twofold_axes_tol(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%tol = new_value
end subroutine set_twofold_axes_tol

subroutine get_twofold_axes_rtwofold(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(out) :: output_value
	output_value = obj_var%rtwofold
end subroutine get_twofold_axes_rtwofold

subroutine set_twofold_axes_rtwofold(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(in) :: new_value
	obj_var%rtwofold = new_value
end subroutine set_twofold_axes_rtwofold

subroutine get_twofold_axes_caxes(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) ,dimension(3,12), intent(out) :: output_value
	output_value = obj_var%caxes
end subroutine get_twofold_axes_caxes

subroutine set_twofold_axes_caxes(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) ,dimension(3,12), intent(in) :: new_value
	obj_var%caxes = new_value
end subroutine set_twofold_axes_caxes

subroutine get_twofold_axes_dot(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(12), intent(out) :: output_value
	output_value = obj_var%dot
end subroutine get_twofold_axes_dot

subroutine set_twofold_axes_dot(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(12), intent(in) :: new_value
	obj_var%dot = new_value
end subroutine set_twofold_axes_dot

subroutine get_twofold_axes_dtwofold(obj_var, output_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(out) :: output_value
	output_value = obj_var%dtwofold
end subroutine get_twofold_axes_dtwofold

subroutine set_twofold_axes_dtwofold(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(in) :: new_value
	obj_var%dtwofold = new_value
end subroutine set_twofold_axes_dtwofold

subroutine Twofold_Axes_Type_ctor(Twofold_Axes_Type_param, a_param, c_param, b_param, maxes_param, ntwo_param, cross_param, tol_param, rtwofold_param, caxes_param, dot_param, dtwofold_param)
	type (Twofold_Axes_Type) :: Twofold_Axes_Type_param
	real(kind=cp), dimension(3), intent(in) :: a_param
	real(kind=cp), dimension(3), intent(in) :: c_param
	real(kind=cp), dimension(3), intent(in) :: b_param
	real(kind=cp), dimension(12), intent(in) :: maxes_param
	integer, intent(in) :: ntwo_param
	real(kind=cp), dimension(12), intent(in) :: cross_param
	real(kind=cp), intent(in) :: tol_param
	integer,dimension(3,12), intent(in) :: rtwofold_param
	real(kind=cp) ,dimension(3,12), intent(in) :: caxes_param
	integer,dimension(12), intent(in) :: dot_param
	integer,dimension(3,12), intent(in) :: dtwofold_param
	Twofold_Axes_Type_param%a = a_param
	Twofold_Axes_Type_param%c = c_param
	Twofold_Axes_Type_param%b = b_param
	Twofold_Axes_Type_param%maxes = maxes_param
	Twofold_Axes_Type_param%ntwo = ntwo_param
	Twofold_Axes_Type_param%cross = cross_param
	Twofold_Axes_Type_param%tol = tol_param
	Twofold_Axes_Type_param%rtwofold = rtwofold_param
	Twofold_Axes_Type_param%caxes = caxes_param
	Twofold_Axes_Type_param%dot = dot_param
	Twofold_Axes_Type_param%dtwofold = dtwofold_param
end subroutine Twofold_Axes_Type_ctor
