function get_maghd_list_Nref(obj_var)
	type (MagHD_List_Type) :: obj_var
	integer :: get_maghd_list_Nref
	get_maghd_list_Nref = obj_var%Nref
end function get_maghd_list_Nref

subroutine set_maghd_list_Nref(obj_var, new_value)
	type (MagHD_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Nref = new_value
end subroutine set_maghd_list_Nref

subroutine get_maghd_list_Mh(obj_var, output_value)
	type (MagHD_List_Type) :: obj_var
	Type(MagHD_Type), dimension(:), intent(out) :: output_value
	output_value = obj_var%Mh
end subroutine get_maghd_list_Mh

subroutine set_maghd_list_Mh(obj_var, new_value)
	type (MagHD_List_Type) :: obj_var
	Type(MagHD_Type), dimension(:), intent(in) :: new_value
	obj_var%Mh = new_value
end subroutine set_maghd_list_Mh

subroutine MagHD_List_Type_ctor(MagHD_List_Type_param, Nref_param, Mh_param)
	type (MagHD_List_Type) :: MagHD_List_Type_param
	integer, intent(in) :: Nref_param
	Type(MagHD_Type), dimension(:), intent(in) :: Mh_param
	MagHD_List_Type_param%Nref = Nref_param
	MagHD_List_Type_param%Mh = Mh_param
end subroutine MagHD_List_Type_ctor
