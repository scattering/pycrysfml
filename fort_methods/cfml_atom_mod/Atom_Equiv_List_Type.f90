function getnauas(obj_var)
	type (Atom_Equiv_List_Type) :: obj_var
	integer :: getnauas
	getnauas = obj_var%nauas
end function getnauas

subroutine setnauas(obj_var, new_value)
	type (Atom_Equiv_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nauas = new_value
end subroutine setnauas

function getatm(obj_var)
	type (Atom_Equiv_List_Type) :: obj_var
	type (Atom_Equiv_Type), allocatable, dimension(:) :: getatm
	getatm = obj_var%atm
end function getatm

subroutine setatm(obj_var, new_value)
	type (Atom_Equiv_List_Type) :: obj_var
	type (Atom_Equiv_Type), allocatable, dimension(:), intent(in) :: new_value
	obj_var%atm = new_value
end subroutine setatm

subroutine Atom_Equiv_List_Type_ctor(Atom_Equiv_List_Type_param, nauas_param, atm_param)
	type (Atom_Equiv_List_Type) :: Atom_Equiv_List_Type_param
	integer, intent(in) :: nauas_param
	type (Atom_Equiv_Type), allocatable, dimension(:), intent(in) :: atm_param
	Atom_Equiv_List_Type_param%nauas = nauas_param
	Atom_Equiv_List_Type_param%atm = atm_param
end subroutine Atom_Equiv_List_Type_ctor
