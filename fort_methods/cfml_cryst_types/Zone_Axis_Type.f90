function get_zone_axis_nlayer(obj_var)
	type (Zone_Axis_Type) :: obj_var
	Integer :: get_zone_axis_nlayer
	get_zone_axis_nlayer = obj_var%nlayer
end function get_zone_axis_nlayer

subroutine set_zone_axis_nlayer(obj_var, new_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%nlayer = new_value
end subroutine set_zone_axis_nlayer

subroutine get_zone_axis_rx(obj_var, output_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(out) :: output_value
	output_value = obj_var%rx
end subroutine get_zone_axis_rx

subroutine set_zone_axis_rx(obj_var, new_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(in) :: new_value
	obj_var%rx = new_value
end subroutine set_zone_axis_rx

subroutine get_zone_axis_ry(obj_var, output_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(out) :: output_value
	output_value = obj_var%ry
end subroutine get_zone_axis_ry

subroutine set_zone_axis_ry(obj_var, new_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(in) :: new_value
	obj_var%ry = new_value
end subroutine set_zone_axis_ry

subroutine get_zone_axis_uvw(obj_var, output_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(out) :: output_value
	output_value = obj_var%uvw
end subroutine get_zone_axis_uvw

subroutine set_zone_axis_uvw(obj_var, new_value)
	type (Zone_Axis_Type) :: obj_var
	Integer, dimension(3), intent(in) :: new_value
	obj_var%uvw = new_value
end subroutine set_zone_axis_uvw

subroutine Zone_Axis_Type_ctor(Zone_Axis_Type_param, nlayer_param, rx_param, ry_param, uvw_param)
	type (Zone_Axis_Type) :: Zone_Axis_Type_param
	Integer, intent(in) :: nlayer_param
	Integer, dimension(3), intent(in) :: rx_param
	Integer, dimension(3), intent(in) :: ry_param
	Integer, dimension(3), intent(in) :: uvw_param
	Zone_Axis_Type_param%nlayer = nlayer_param
	Zone_Axis_Type_param%rx = rx_param
	Zone_Axis_Type_param%ry = ry_param
	Zone_Axis_Type_param%uvw = uvw_param
end subroutine Zone_Axis_Type_ctor
