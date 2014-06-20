function getTr(obj_var)
	type (Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3) :: getTr
	getTr = obj_var%Tr
end function getTr

subroutine setTr(obj_var, new_value)
	type (Sym_Oper_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Tr = new_value
end subroutine setTr

function getRot(obj_var)
	type (Sym_Oper_Type) :: obj_var
	integer,       dimension(3,3) :: getRot
	getRot = obj_var%Rot
end function getRot

subroutine setRot(obj_var, new_value)
	type (Sym_Oper_Type) :: obj_var
	integer,       dimension(3,3), intent(in) :: new_value
	obj_var%Rot = new_value
end subroutine setRot

subroutine Sym_Oper_Type_ctor(Sym_Oper_Type_param, Tr_param, Rot_param)
	type (Sym_Oper_Type) :: Sym_Oper_Type_param
	real(kind=cp), dimension(3), intent(in) :: Tr_param
	integer,       dimension(3,3), intent(in) :: Rot_param
	Sym_Oper_Type_param%Tr = Tr_param
	Sym_Oper_Type_param%Rot = Rot_param
end subroutine Sym_Oper_Type_ctor
