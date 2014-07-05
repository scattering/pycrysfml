function get_msym_oper_Phas(obj_var)
	type (MSym_Oper_Type) :: obj_var
	real(kind=cp) :: get_msym_oper_Phas
	get_msym_oper_Phas = obj_var%Phas
end function get_msym_oper_Phas

subroutine set_msym_oper_Phas(obj_var, new_value)
	type (MSym_Oper_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Phas = new_value
end subroutine set_msym_oper_Phas

subroutine get_msym_oper_Rot(obj_var, output_value)
	type (MSym_Oper_Type) :: obj_var
	integer, dimension(3,3), intent(out) :: output_value
	output_value = obj_var%Rot
end subroutine get_msym_oper_Rot

subroutine set_msym_oper_Rot(obj_var, new_value)
	type (MSym_Oper_Type) :: obj_var
	integer, dimension(3,3), intent(in) :: new_value
	obj_var%Rot = new_value
end subroutine set_msym_oper_Rot

subroutine MSym_Oper_Type_ctor(MSym_Oper_Type_param, Phas_param, Rot_param)
	type (MSym_Oper_Type) :: MSym_Oper_Type_param
	real(kind=cp), intent(in) :: Phas_param
	integer, dimension(3,3), intent(in) :: Rot_param
	MSym_Oper_Type_param%Phas = Phas_param
	MSym_Oper_Type_param%Rot = Rot_param
end subroutine MSym_Oper_Type_ctor
