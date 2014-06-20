function getnorb(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	integer :: getnorb
	getnorb = obj_var%norb
end function getnorb

subroutine setnorb(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%norb = new_value
end subroutine setnorb

function getstr_orig(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40) :: getstr_orig
	getstr_orig = obj_var%str_orig
end function getstr_orig

subroutine setstr_orig(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%str_orig = new_value
end subroutine setstr_orig

function getmultp(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	integer :: getmultp
	getmultp = obj_var%multp
end function getmultp

subroutine setmultp(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%multp = new_value
end subroutine setmultp

function getsite(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	character(len= 6) :: getsite
	getsite = obj_var%site
end function getsite

subroutine setsite(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len= 6), intent(in) :: new_value
	obj_var%site = new_value
end subroutine setsite

function getstr_orbit(obj_var)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40),dimension(48) :: getstr_orbit
	getstr_orbit = obj_var%str_orbit
end function getstr_orbit

subroutine setstr_orbit(obj_var, new_value)
	type (Wyck_Pos_Type) :: obj_var
	character(len=40),dimension(48), intent(in) :: new_value
	obj_var%str_orbit = new_value
end subroutine setstr_orbit

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
