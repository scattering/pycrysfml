function getnum_orbit(obj_var)
	type (Wyckoff_Type) :: obj_var
	integer :: getnum_orbit
	getnum_orbit = obj_var%num_orbit
end function getnum_orbit

subroutine setnum_orbit(obj_var, new_value)
	type (Wyckoff_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%num_orbit = new_value
end subroutine setnum_orbit

function getorbit(obj_var)
	type (Wyckoff_Type) :: obj_var
	type(wyck_pos_type), dimension(26) :: getorbit
	getorbit = obj_var%orbit
end function getorbit

subroutine setorbit(obj_var, new_value)
	type (Wyckoff_Type) :: obj_var
	type(wyck_pos_type), dimension(26), intent(in) :: new_value
	obj_var%orbit = new_value
end subroutine setorbit

subroutine Wyckoff_Type_ctor(Wyckoff_Type_param, num_orbit_param, orbit_param)
	type (Wyckoff_Type) :: Wyckoff_Type_param
	integer, intent(in) :: num_orbit_param
	type(wyck_pos_type), dimension(26), intent(in) :: orbit_param
	Wyckoff_Type_param%num_orbit = num_orbit_param
	Wyckoff_Type_param%orbit = orbit_param
end subroutine Wyckoff_Type_ctor
