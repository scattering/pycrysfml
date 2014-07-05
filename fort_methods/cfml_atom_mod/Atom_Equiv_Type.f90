subroutine get_atom_equiv_x(obj_var, output_value)
	type (Atom_Equiv_Type) :: obj_var
	real(kind=cp),     dimension(:,:), intent(out) :: output_value
	output_value = obj_var%x
end subroutine get_atom_equiv_x

subroutine set_atom_equiv_x(obj_var, new_value)
	type (Atom_Equiv_Type) :: obj_var
	real(kind=cp),     dimension(:,:), intent(in) :: new_value
	obj_var%x = new_value
end subroutine set_atom_equiv_x

subroutine get_atom_equiv_ChemSymb(obj_var, output_value)
	type (Atom_Equiv_Type) :: obj_var
	character(len=2), intent(out) :: output_value
	output_value = obj_var%ChemSymb
end subroutine get_atom_equiv_ChemSymb

subroutine set_atom_equiv_ChemSymb(obj_var, new_value)
	type (Atom_Equiv_Type) :: obj_var
	character(len=2), intent(in) :: new_value
	obj_var%ChemSymb = new_value
end subroutine set_atom_equiv_ChemSymb

function get_atom_equiv_mult(obj_var)
	type (Atom_Equiv_Type) :: obj_var
	integer :: get_atom_equiv_mult
	get_atom_equiv_mult = obj_var%mult
end function get_atom_equiv_mult

subroutine set_atom_equiv_mult(obj_var, new_value)
	type (Atom_Equiv_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%mult = new_value
end subroutine set_atom_equiv_mult

subroutine get_atom_equiv_Lab(obj_var, output_value)
	type (Atom_Equiv_Type) :: obj_var
	character(len=20), dimension(:), intent(out) :: output_value
	output_value = obj_var%Lab
end subroutine get_atom_equiv_Lab

subroutine set_atom_equiv_Lab(obj_var, new_value)
	type (Atom_Equiv_Type) :: obj_var
	character(len=20), dimension(:), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine set_atom_equiv_Lab

subroutine Atom_Equiv_Type_ctor(Atom_Equiv_Type_param, x_param, ChemSymb_param, mult_param, Lab_param)
	type (Atom_Equiv_Type) :: Atom_Equiv_Type_param
	real(kind=cp),     dimension(:,:), intent(in) :: x_param
	character(len=2), intent(in) :: ChemSymb_param
	integer, intent(in) :: mult_param
	character(len=20), dimension(:), intent(in) :: Lab_param
	Atom_Equiv_Type_param%x = x_param
	Atom_Equiv_Type_param%ChemSymb = ChemSymb_param
	Atom_Equiv_Type_param%mult = mult_param
	Atom_Equiv_Type_param%Lab = Lab_param
end subroutine Atom_Equiv_Type_ctor
