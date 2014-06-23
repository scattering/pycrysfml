function get_reflect_H(obj_var)
	type (Reflect_Type) :: obj_var
	integer,dimension(3) :: getH
	getH = obj_var%H
end function get_reflect_H

subroutine set_reflect_H(obj_var, new_value)
	type (Reflect_Type) :: obj_var
	integer,dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine set_reflect_H

function get_reflect_S(obj_var)
	type (Reflect_Type) :: obj_var
	real(kind=cp) :: getS
	getS = obj_var%S
end function get_reflect_S

subroutine set_reflect_S(obj_var, new_value)
	type (Reflect_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%S = new_value
end subroutine set_reflect_S

function get_reflect_Mult(obj_var)
	type (Reflect_Type) :: obj_var
	integer :: getMult
	getMult = obj_var%Mult
end function get_reflect_Mult

subroutine set_reflect_Mult(obj_var, new_value)
	type (Reflect_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine set_reflect_Mult

subroutine Reflect_Type_ctor(Reflect_Type_param, H_param, S_param, Mult_param)
	type (Reflect_Type) :: Reflect_Type_param
	integer,dimension(3), intent(in) :: H_param
	real(kind=cp), intent(in) :: S_param
	integer, intent(in) :: Mult_param
	Reflect_Type_param%H = H_param
	Reflect_Type_param%S = S_param
	Reflect_Type_param%Mult = Mult_param
end subroutine Reflect_Type_ctor
