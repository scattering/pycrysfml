function get_magh_list_Nref(obj_var)
	type (MagH_List_Type) :: obj_var
	integer :: getNref
	getNref = obj_var%Nref
end function get_magh_list_Nref

subroutine set_magh_list_Nref(obj_var, new_value)
	type (MagH_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Nref = new_value
end subroutine set_magh_list_Nref

function get_magh_list_Mh(obj_var)
	type (MagH_List_Type) :: obj_var
	Type(MagH_Type),allocatable, dimension(:) :: getMh
	getMh = obj_var%Mh
end function get_magh_list_Mh

subroutine set_magh_list_Mh(obj_var, new_value)
	type (MagH_List_Type) :: obj_var
	Type(MagH_Type),allocatable, dimension(:), intent(in) :: new_value
	obj_var%Mh = new_value
end subroutine set_magh_list_Mh

subroutine MagH_List_Type_ctor(MagH_List_Type_param, Nref_param, Mh_param)
	type (MagH_List_Type) :: MagH_List_Type_param
	integer, intent(in) :: Nref_param
	Type(MagH_Type),allocatable, dimension(:), intent(in) :: Mh_param
	MagH_List_Type_param%Nref = Nref_param
	MagH_List_Type_param%Mh = Mh_param
end subroutine MagH_List_Type_ctor
