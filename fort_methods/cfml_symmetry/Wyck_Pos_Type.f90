function get_wyck_pos_norb(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	integer :: get_wyck_pos_norb
	get_wyck_pos_norb = obj_var%norb
end function get_wyck_pos_norb

subroutine set_wyck_pos_norb(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%norb = new_value
end subroutine set_wyck_pos_norb

subroutine get_wyck_pos_str_orig(obj_var, output_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40), intent(out) :: output_value
	output_value = obj_var%str_orig
end subroutine get_wyck_pos_str_orig

subroutine set_wyck_pos_str_orig(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%str_orig = new_value
end subroutine set_wyck_pos_str_orig

function get_wyck_pos_multp(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	integer :: get_wyck_pos_multp
	get_wyck_pos_multp = obj_var%multp
end function get_wyck_pos_multp

subroutine set_wyck_pos_multp(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%multp = new_value
end subroutine set_wyck_pos_multp

subroutine get_wyck_pos_site(obj_var, output_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len= 6), intent(out) :: output_value
	output_value = obj_var%site
end subroutine get_wyck_pos_site

subroutine set_wyck_pos_site(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len= 6), intent(in) :: new_value
	obj_var%site = new_value
end subroutine set_wyck_pos_site

subroutine get_wyck_pos_str_orbit(obj_var, output_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40),dimension(48), intent(out) :: output_value
	output_value = obj_var%str_orbit
end subroutine get_wyck_pos_str_orbit

subroutine set_wyck_pos_str_orbit(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40),dimension(48), intent(in) :: new_value
	obj_var%str_orbit = new_value
end subroutine set_wyck_pos_str_orbit

subroutine Wyck_Pos_Type_ctor(Wyck_Pos_Type_param, norb_param, str_orig_param, multp_param, site_param, str_orbit_param)
	type (Wyck_Pos_Type) :: Wyck_Pos_Type_param
	integer, intent(in) :: norb_param
	character(len=40), intent(in) :: str_orig_param
	integer, intent(in) :: multp_param
	character(len= 6), intent(in) :: site_param
	character(len=40),dimension(48), intent(in) :: str_orbit_param
	Wyck_Pos_Type_param%norb = norb_param
	Wyck_Pos_Type_param%str_orig = str_orig_param
	Wyck_Pos_Type_param%multp = multp_param
	Wyck_Pos_Type_param%site = site_param
	Wyck_Pos_Type_param%str_orbit = str_orbit_param
end subroutine Wyck_Pos_Type_ctor
