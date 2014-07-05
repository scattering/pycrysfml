subroutine get_sym_oper_Tr(obj_var, output_value)
	type (Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%Tr
end subroutine get_sym_oper_Tr

subroutine set_sym_oper_Tr(obj_var, new_value)
	type (Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Tr = new_value
end subroutine set_sym_oper_Tr

subroutine get_sym_oper_Rot(obj_var, output_value)
	type (Sym_Oper_Type) :: obj_var
	integer,       dimension(3,3), intent(out) :: output_value
	output_value = obj_var%Rot
end subroutine get_sym_oper_Rot

subroutine set_sym_oper_Rot(obj_var, new_value)
	type (Sym_Oper_Type) :: obj_var
	integer,       dimension(3,3), intent(in) :: new_value
	obj_var%Rot = new_value
end subroutine set_sym_oper_Rot

subroutine Sym_Oper_Type_ctor(Sym_Oper_Type_param, Tr_param, Rot_param)
	type (Sym_Oper_Type) :: Sym_Oper_Type_param
	real(kind=cp), dimension(3), intent(in) :: Tr_param
	integer,       dimension(3,3), intent(in) :: Rot_param
	Sym_Oper_Type_param%Tr = Tr_param
	Sym_Oper_Type_param%Rot = Rot_param
end subroutine Sym_Oper_Type_ctor
