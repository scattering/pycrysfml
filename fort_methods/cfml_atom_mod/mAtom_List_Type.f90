function get_matom_list_natoms(obj_var)
	type (mAtom_List_Type) :: obj_var
	integer :: get_matom_list_natoms
	get_matom_list_natoms = obj_var%natoms
end function get_matom_list_natoms

subroutine set_matom_list_natoms(obj_var, new_value)
	type (mAtom_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%natoms = new_value
end subroutine set_matom_list_natoms

subroutine get_matom_list_Atom(obj_var, output_value)
	type (mAtom_List_Type) :: obj_var
	type(mAtom_Type),dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%Atom
end subroutine get_matom_list_Atom

subroutine set_matom_list_Atom(obj_var, new_value)
	type (mAtom_List_Type) :: obj_var
	type(mAtom_Type),dimension(:),allocatable, intent(in) :: new_value
	obj_var%Atom = new_value
end subroutine set_matom_list_Atom

subroutine mAtom_List_Type_ctor(mAtom_List_Type_param, natoms_param, Atom_param)
	type (mAtom_List_Type) :: mAtom_List_Type_param
	integer, intent(in) :: natoms_param
	type(mAtom_Type),dimension(:),allocatable, intent(in) :: Atom_param
	mAtom_List_Type_param%natoms = natoms_param
	mAtom_List_Type_param%Atom = Atom_param
end subroutine mAtom_List_Type_ctor
