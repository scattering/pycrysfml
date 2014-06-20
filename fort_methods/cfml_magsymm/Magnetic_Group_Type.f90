function gettinv(obj_var)
	type (Magnetic_Group_Type) :: obj_var
	integer, dimension(192) :: gettinv
	gettinv = obj_var%tinv
end function gettinv

subroutine settinv(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	integer, dimension(192), intent(in) :: new_value
	obj_var%tinv = new_value
end subroutine settinv

function getShubnikov(obj_var)
	type (Magnetic_Group_Type) :: obj_var
	Character(len=30) :: getShubnikov
	getShubnikov = obj_var%Shubnikov
end function getShubnikov

subroutine setShubnikov(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	Character(len=30), intent(in) :: new_value
	obj_var%Shubnikov = new_value
end subroutine setShubnikov

function getSpG(obj_var)
	type (Magnetic_Group_Type) :: obj_var
	type(Space_Group_Type) :: getSpG
	getSpG = obj_var%SpG
end function getSpG

subroutine setSpG(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	type(Space_Group_Type), intent(in) :: new_value
	obj_var%SpG = new_value
end subroutine setSpG

subroutine Magnetic_Group_Type_ctor(Magnetic_Group_Type_param, tinv_param, Shubnikov_param, SpG_param)
	type (Magnetic_Group_Type) :: Magnetic_Group_Type_param
	integer, dimension(192), intent(in) :: tinv_param
	Character(len=30), intent(in) :: Shubnikov_param
	type(Space_Group_Type), intent(in) :: SpG_param
	Magnetic_Group_Type_param%tinv = tinv_param
	Magnetic_Group_Type_param%Shubnikov = Shubnikov_param
	Magnetic_Group_Type_param%SpG = SpG_param
end subroutine Magnetic_Group_Type_ctor
