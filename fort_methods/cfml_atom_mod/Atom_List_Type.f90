function getnatoms(obj_var)
	type (Atom_List_Type) :: obj_var
	integer :: getnatoms
	getnatoms = obj_var%natoms
end function getnatoms

subroutine setnatoms(obj_var, new_value)
	type (Atom_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%natoms = new_value
end subroutine setnatoms

function getatom(obj_var)
	type (Atom_List_Type) :: obj_var
	type(Atom_Type),dimension(:),allocatable :: getatom
	getatom = obj_var%atom
end function getatom

subroutine setatom(obj_var, new_value)
	type (Atom_List_Type) :: obj_var
	type(Atom_Type),dimension(:),allocatable, intent(in) :: new_value
	obj_var%atom = new_value
end subroutine setatom

subroutine Atom_List_Type_ctor(Atom_List_Type_param, natoms_param, atom_param)
	type (Atom_List_Type) :: Atom_List_Type_param
	integer, intent(in) :: natoms_param
	type(Atom_Type),dimension(:),allocatable, intent(in) :: atom_param
	Atom_List_Type_param%natoms = natoms_param
	Atom_List_Type_param%atom = atom_param
end subroutine Atom_List_Type_ctor
