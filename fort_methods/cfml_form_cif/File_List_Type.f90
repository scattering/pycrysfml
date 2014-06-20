function getline(obj_var)
	type (File_List_Type) :: obj_var
	character(len=256), allocatable, dimension(:) :: getline
	getline = obj_var%line
end function getline

subroutine setline(obj_var, new_value)
	type (File_List_Type) :: obj_var
	character(len=256), allocatable, dimension(:), intent(in) :: new_value
	obj_var%line = new_value
end subroutine setline

function getnlines(obj_var)
	type (File_List_Type) :: obj_var
	integer :: getnlines
	getnlines = obj_var%nlines
end function getnlines

subroutine setnlines(obj_var, new_value)
	type (File_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nlines = new_value
end subroutine setnlines

subroutine File_List_Type_ctor(File_List_Type_param, line_param, nlines_param)
	type (File_List_Type) :: File_List_Type_param
	character(len=256), allocatable, dimension(:), intent(in) :: line_param
	integer, intent(in) :: nlines_param
	File_List_Type_param%line = line_param
	File_List_Type_param%nlines = nlines_param
end subroutine File_List_Type_ctor
