subroutine get_reflection_list_Ref(obj_var, output_value)
	type (Reflection_List_Type) :: obj_var
	type(reflection_type), dimension(:), intent(out) :: output_value
	output_value = obj_var%Ref
end subroutine get_reflection_list_Ref

subroutine set_reflection_list_Ref(obj_var, new_value)
	type (Reflection_List_Type) :: obj_var
	type(reflection_type), dimension(:), intent(in) :: new_value
	obj_var%Ref = new_value
end subroutine set_reflection_list_Ref

function get_reflection_list_NRef(obj_var)
	type (Reflection_List_Type) :: obj_var
	integer :: get_reflection_list_NRef
	get_reflection_list_NRef = obj_var%NRef
end function get_reflection_list_NRef

subroutine set_reflection_list_NRef(obj_var, new_value)
	type (Reflection_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NRef = new_value
end subroutine set_reflection_list_NRef

subroutine Reflection_List_Type_ctor(Reflection_List_Type_param, Ref_param, NRef_param)
	type (Reflection_List_Type) :: Reflection_List_Type_param
	type(reflection_type), dimension(:), intent(in) :: Ref_param
	integer, intent(in) :: NRef_param
	Reflection_List_Type_param%Ref = Ref_param
	Reflection_List_Type_param%NRef = NRef_param
end subroutine Reflection_List_Type_ctor
