function get_magnetic_domain_Chir(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: get_magnetic_domain_Chir
	get_magnetic_domain_Chir = obj_var%Chir
end function get_magnetic_domain_Chir

subroutine set_magnetic_domain_Chir(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Chir = new_value
end subroutine set_magnetic_domain_Chir

subroutine get_magnetic_domain_pop_std(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(out) :: output_value
	output_value = obj_var%pop_std
end subroutine get_magnetic_domain_pop_std

subroutine set_magnetic_domain_pop_std(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%pop_std = new_value
end subroutine set_magnetic_domain_pop_std

subroutine get_magnetic_domain_Mpop(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(out) :: output_value
	output_value = obj_var%Mpop
end subroutine get_magnetic_domain_Mpop

subroutine set_magnetic_domain_Mpop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%Mpop = new_value
end subroutine set_magnetic_domain_Mpop

function get_magnetic_domain_nd(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	integer :: get_magnetic_domain_nd
	get_magnetic_domain_nd = obj_var%nd
end function get_magnetic_domain_nd

subroutine set_magnetic_domain_nd(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nd = new_value
end subroutine set_magnetic_domain_nd

subroutine get_magnetic_domain_Lpop(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	integer      , dimension (2,24), intent(out) :: output_value
	output_value = obj_var%Lpop
end subroutine get_magnetic_domain_Lpop

subroutine set_magnetic_domain_Lpop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer      , dimension (2,24), intent(in) :: new_value
	obj_var%Lpop = new_value
end subroutine set_magnetic_domain_Lpop

subroutine get_magnetic_domain_pop(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(out) :: output_value
	output_value = obj_var%pop
end subroutine get_magnetic_domain_pop

subroutine set_magnetic_domain_pop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%pop = new_value
end subroutine set_magnetic_domain_pop

function get_magnetic_domain_Twin(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: get_magnetic_domain_Twin
	get_magnetic_domain_Twin = obj_var%Twin
end function get_magnetic_domain_Twin

subroutine set_magnetic_domain_Twin(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Twin = new_value
end subroutine set_magnetic_domain_Twin

subroutine get_magnetic_domain_Lab(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	character(len=10),dimension (2,24), intent(out) :: output_value
	output_value = obj_var%Lab
end subroutine get_magnetic_domain_Lab

subroutine set_magnetic_domain_Lab(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	character(len=10),dimension (2,24), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine set_magnetic_domain_Lab

subroutine get_magnetic_domain_DMat(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	integer,dimension(3,3,24), intent(out) :: output_value
	output_value = obj_var%DMat
end subroutine get_magnetic_domain_DMat

subroutine set_magnetic_domain_DMat(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer,dimension(3,3,24), intent(in) :: new_value
	obj_var%DMat = new_value
end subroutine set_magnetic_domain_DMat

subroutine get_magnetic_domain_Dt(obj_var, output_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (3,24), intent(out) :: output_value
	output_value = obj_var%Dt
end subroutine get_magnetic_domain_Dt

subroutine set_magnetic_domain_Dt(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (3,24), intent(in) :: new_value
	obj_var%Dt = new_value
end subroutine set_magnetic_domain_Dt

function get_magnetic_domain_trans(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: get_magnetic_domain_trans
	get_magnetic_domain_trans = obj_var%trans
end function get_magnetic_domain_trans

subroutine set_magnetic_domain_trans(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%trans = new_value
end subroutine set_magnetic_domain_trans

subroutine Magnetic_Domain_type_ctor(Magnetic_Domain_type_param, Chir_param, pop_std_param, Mpop_param, nd_param, Lpop_param, pop_param, Twin_param, Lab_param, DMat_param, Dt_param, trans_param)
	type (Magnetic_Domain_type) :: Magnetic_Domain_type_param
	logical, intent(in) :: Chir_param
	real(kind=cp), dimension (2,24), intent(in) :: pop_std_param
	real(kind=cp), dimension (2,24), intent(in) :: Mpop_param
	integer, intent(in) :: nd_param
	integer      , dimension (2,24), intent(in) :: Lpop_param
	real(kind=cp), dimension (2,24), intent(in) :: pop_param
	logical, intent(in) :: Twin_param
	character(len=10),dimension (2,24), intent(in) :: Lab_param
	integer,dimension(3,3,24), intent(in) :: DMat_param
	real(kind=cp), dimension (3,24), intent(in) :: Dt_param
	logical, intent(in) :: trans_param
	Magnetic_Domain_type_param%Chir = Chir_param
	Magnetic_Domain_type_param%pop_std = pop_std_param
	Magnetic_Domain_type_param%Mpop = Mpop_param
	Magnetic_Domain_type_param%nd = nd_param
	Magnetic_Domain_type_param%Lpop = Lpop_param
	Magnetic_Domain_type_param%pop = pop_param
	Magnetic_Domain_type_param%Twin = Twin_param
	Magnetic_Domain_type_param%Lab = Lab_param
	Magnetic_Domain_type_param%DMat = DMat_param
	Magnetic_Domain_type_param%Dt = Dt_param
	Magnetic_Domain_type_param%trans = trans_param
end subroutine Magnetic_Domain_type_ctor
