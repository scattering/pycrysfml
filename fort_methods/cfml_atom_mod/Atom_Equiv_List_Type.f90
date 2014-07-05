function get_atom_equiv_list_nauas(obj_var)
	type (Atom_Equiv_List_Type) :: obj_var
	integer :: get_atom_equiv_list_nauas
	get_atom_equiv_list_nauas = obj_var%nauas
end function get_atom_equiv_list_nauas

subroutine set_atom_equiv_list_nauas(obj_var, new_value)
	type (Atom_Equiv_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nauas = new_value
end subroutine set_atom_equiv_list_nauas

subroutine get_atom_equiv_list_atm(obj_var, output_value)
	type (Atom_Equiv_List_Type) :: obj_var
	type (Atom_Equiv_Type),  dimension(:), intent(out) :: output_value
	output_value = obj_var%atm
end subroutine get_atom_equiv_list_atm

subroutine set_atom_equiv_list_atm(obj_var, new_value)
	type (Atom_Equiv_List_Type) :: obj_var
	type (Atom_Equiv_Type),  dimension(:), intent(in) :: new_value
	obj_var%atm = new_value
end subroutine set_atom_equiv_list_atm

subroutine Atom_Equiv_List_Type_ctor(Atom_Equiv_List_Type_param, nauas_param, atm_param)
	type (Atom_Equiv_List_Type) :: Atom_Equiv_List_Type_param
	integer, intent(in) :: nauas_param
	type (Atom_Equiv_Type),  dimension(:), intent(in) :: atm_param
	Atom_Equiv_List_Type_param%nauas = nauas_param
	Atom_Equiv_List_Type_param%atm = atm_param
end subroutine Atom_Equiv_List_Type_ctor
