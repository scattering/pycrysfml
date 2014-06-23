function get_ns_sym_oper_Tr(obj_var)
	type (NS_Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3) :: getTr
	getTr = obj_var%Tr
end function get_ns_sym_oper_Tr

subroutine set_ns_sym_oper_Tr(obj_var, new_value)
	type (NS_Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Tr = new_value
end subroutine set_ns_sym_oper_Tr

function get_ns_sym_oper_Rot(obj_var)
	type (NS_Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3,3) :: getRot
	getRot = obj_var%Rot
end function get_ns_sym_oper_Rot

subroutine set_ns_sym_oper_Rot(obj_var, new_value)
	type (NS_Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3,3), intent(in) :: new_value
	obj_var%Rot = new_value
end subroutine set_ns_sym_oper_Rot

subroutine NS_Sym_Oper_Type_ctor(NS_Sym_Oper_Type_param, Tr_param, Rot_param)
	type (NS_Sym_Oper_Type) :: NS_Sym_Oper_Type_param
	real(kind=cp), dimension(3), intent(in) :: Tr_param
	real(kind=cp), dimension(3,3), intent(in) :: Rot_param
	NS_Sym_Oper_Type_param%Tr = Tr_param
	NS_Sym_Oper_Type_param%Rot = Rot_param
end subroutine NS_Sym_Oper_Type_ctor
