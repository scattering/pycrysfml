function get_atom_list_natoms(obj_var)
	type (Atom_List_Type) :: obj_var
	integer :: get_atom_list_natoms
	get_atom_list_natoms = obj_var%natoms
end function get_atom_list_natoms

subroutine set_atom_list_natoms(obj_var, new_value)
	type (Atom_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%natoms = new_value
end subroutine set_atom_list_natoms

subroutine get_atom_list_atom(obj_var, output_value)
	type (Atom_List_Type) :: obj_var
	type(Atom_Type),dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%atom
end subroutine get_atom_list_atom

subroutine set_atom_list_atom(obj_var, new_value)
	type (Atom_List_Type) :: obj_var
	type(Atom_Type),dimension(:),allocatable, intent(in) :: new_value
	obj_var%atom = new_value
end subroutine set_atom_list_atom

subroutine Atom_List_Type_ctor(Atom_List_Type_param, natoms_param, atom_param)
	type (Atom_List_Type) :: Atom_List_Type_param
	integer, intent(in) :: natoms_param
	type(Atom_Type),dimension(:),allocatable, intent(in) :: atom_param
	Atom_List_Type_param%natoms = natoms_param
	Atom_List_Type_param%atom = atom_param
end subroutine Atom_List_Type_ctor
