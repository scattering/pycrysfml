subroutine get_file_list_line(obj_var, output_value)
	type (File_List_Type) :: obj_var
	character(len=256),  dimension(:), intent(out) :: output_value
	output_value = obj_var%line
end subroutine get_file_list_line

subroutine set_file_list_line(obj_var, new_value)
	type (File_List_Type) :: obj_var
	character(len=256),  dimension(:), intent(in) :: new_value
	obj_var%line = new_value
end subroutine set_file_list_line

function get_file_list_nlines(obj_var)
	type (File_List_Type) :: obj_var
	integer :: get_file_list_nlines
	get_file_list_nlines = obj_var%nlines
end function get_file_list_nlines

subroutine set_file_list_nlines(obj_var, new_value)
	type (File_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nlines = new_value
end subroutine set_file_list_nlines

subroutine File_List_Type_ctor(File_List_Type_param, line_param, nlines_param)
	type (File_List_Type) :: File_List_Type_param
	character(len=256),  dimension(:), intent(in) :: line_param
	integer, intent(in) :: nlines_param
	File_List_Type_param%line = line_param
	File_List_Type_param%nlines = nlines_param
end subroutine File_List_Type_ctor
