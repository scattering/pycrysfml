function get_wyckoff_num_orbit(obj_var)
	type (Wyckoff_Type) :: obj_var
	integer :: get_wyckoff_num_orbit
	get_wyckoff_num_orbit = obj_var%num_orbit
end function get_wyckoff_num_orbit

subroutine set_wyckoff_num_orbit(obj_var, new_value)
	type (Wyckoff_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_orbit = new_value
end subroutine set_wyckoff_num_orbit

subroutine get_wyckoff_orbit(obj_var, output_value)
	type (Wyckoff_Type) :: obj_var
	type(wyck_pos_type), dimension(26), intent(out) :: output_value
	output_value = obj_var%orbit
end subroutine get_wyckoff_orbit

subroutine set_wyckoff_orbit(obj_var, new_value)
	type (Wyckoff_Type) :: obj_var
	type(wyck_pos_type), dimension(26), intent(in) :: new_value
	obj_var%orbit = new_value
end subroutine set_wyckoff_orbit

subroutine Wyckoff_Type_ctor(Wyckoff_Type_param, num_orbit_param, orbit_param)
	type (Wyckoff_Type) :: Wyckoff_Type_param
	integer, intent(in) :: num_orbit_param
	type(wyck_pos_type), dimension(26), intent(in) :: orbit_param
	Wyckoff_Type_param%num_orbit = num_orbit_param
	Wyckoff_Type_param%orbit = orbit_param
end subroutine Wyckoff_Type_ctor
