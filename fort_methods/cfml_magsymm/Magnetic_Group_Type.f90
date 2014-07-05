subroutine get_magnetic_group_tinv(obj_var, output_value)
	type (Magnetic_Group_Type) :: obj_var
	integer, dimension(192), intent(out) :: output_value
	output_value = obj_var%tinv
end subroutine get_magnetic_group_tinv

subroutine set_magnetic_group_tinv(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	integer, dimension(192), intent(in) :: new_value
	obj_var%tinv = new_value
end subroutine set_magnetic_group_tinv

subroutine get_magnetic_group_Shubnikov(obj_var, output_value)
	type (Magnetic_Group_Type) :: obj_var
	Character(len=30), intent(out) :: output_value
	output_value = obj_var%Shubnikov
end subroutine get_magnetic_group_Shubnikov

subroutine set_magnetic_group_Shubnikov(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	Character(len=30), intent(in) :: new_value
	obj_var%Shubnikov = new_value
end subroutine set_magnetic_group_Shubnikov

subroutine get_magnetic_group_SpG(obj_var, output_value)
	type (Magnetic_Group_Type) :: obj_var
	type(Space_Group_Type), intent(out) :: output_value
	output_value = obj_var%SpG
end subroutine get_magnetic_group_SpG

subroutine set_magnetic_group_SpG(obj_var, new_value)
	type (Magnetic_Group_Type) :: obj_var
	type(Space_Group_Type), intent(in) :: new_value
	obj_var%SpG = new_value
end subroutine set_magnetic_group_SpG

subroutine Magnetic_Group_Type_ctor(Magnetic_Group_Type_param, tinv_param, Shubnikov_param, SpG_param)
	type (Magnetic_Group_Type) :: Magnetic_Group_Type_param
	integer, dimension(192), intent(in) :: tinv_param
	Character(len=30), intent(in) :: Shubnikov_param
	type(Space_Group_Type), intent(in) :: SpG_param
	Magnetic_Group_Type_param%tinv = tinv_param
	Magnetic_Group_Type_param%Shubnikov = Shubnikov_param
	Magnetic_Group_Type_param%SpG = SpG_param
end subroutine Magnetic_Group_Type_ctor
