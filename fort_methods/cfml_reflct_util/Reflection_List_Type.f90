function getRef(obj_var)
	type (Reflection_List_Type) :: obj_var
	type(reflection_type),allocatable, dimension(:) :: getRef
	getRef = obj_var%Ref
end function getRef

subroutine setRef(obj_var, new_value)
	type (Reflection_List_Type) :: obj_var
	type(reflection_type),allocatable, dimension(:), intent(in) :: new_value
	obj_var%Ref = new_value
end subroutine setRef

function getNRef(obj_var)
	type (Reflection_List_Type) :: obj_var
	integer :: getNRef
	getNRef = obj_var%NRef
end function getNRef

subroutine setNRef(obj_var, new_value)
	type (Reflection_List_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NRef = new_value
end subroutine setNRef

subroutine Reflection_List_Type_ctor(Reflection_List_Type_param, Ref_param, NRef_param)
	type (Reflection_List_Type) :: Reflection_List_Type_param
	type(reflection_type),allocatable, dimension(:), intent(in) :: Ref_param
	integer, intent(in) :: NRef_param
	Reflection_List_Type_param%Ref = Ref_param
	Reflection_List_Type_param%NRef = NRef_param
end subroutine Reflection_List_Type_ctor
