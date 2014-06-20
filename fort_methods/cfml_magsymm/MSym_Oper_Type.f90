function getPhas(obj_var)
	type (MSym_Oper_Type) :: obj_var
	real(kind=cp) :: getPhas
	getPhas = obj_var%Phas
end function getPhas

subroutine setPhas(obj_var, new_value)
	type (MSym_Oper_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Phas = new_value
end subroutine setPhas

function getRot(obj_var)
	type (MSym_Oper_Type) :: obj_var
	integer, dimension(3,3) :: getRot
	getRot = obj_var%Rot
end function getRot

subroutine setRot(obj_var, new_value)
	type (MSym_Oper_Type) :: obj_var
	integer, dimension(3,3), intent(in) :: new_value
	obj_var%Rot = new_value
end subroutine setRot

subroutine MSym_Oper_Type_ctor(MSym_Oper_Type_param, Phas_param, Rot_param)
	type (MSym_Oper_Type) :: MSym_Oper_Type_param
	real(kind=cp), intent(in) :: Phas_param
	integer, dimension(3,3), intent(in) :: Rot_param
	MSym_Oper_Type_param%Phas = Phas_param
	MSym_Oper_Type_param%Rot = Rot_param
end subroutine MSym_Oper_Type_ctor
