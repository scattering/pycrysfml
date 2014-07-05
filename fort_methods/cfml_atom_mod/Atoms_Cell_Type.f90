subroutine get_atoms_cell_distance(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),        dimension( :,:), allocatable, intent(out) :: output_value
	output_value = obj_var%distance
end subroutine get_atoms_cell_distance

subroutine set_atoms_cell_distance(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),        dimension( :,:), allocatable, intent(in) :: new_value
	obj_var%distance = new_value
end subroutine set_atoms_cell_distance

function get_atoms_cell_ndist(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer :: get_atoms_cell_ndist
	get_atoms_cell_ndist = obj_var%ndist
end function get_atoms_cell_ndist

subroutine set_atoms_cell_ndist(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%ndist = new_value
end subroutine set_atoms_cell_ndist

subroutine get_atoms_cell_ddist(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%ddist
end subroutine get_atoms_cell_ddist

subroutine set_atoms_cell_ddist(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%ddist = new_value
end subroutine set_atoms_cell_ddist

subroutine get_atoms_cell_xyz(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(out) :: output_value
	output_value = obj_var%xyz
end subroutine get_atoms_cell_xyz

subroutine set_atoms_cell_xyz(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%xyz = new_value
end subroutine set_atoms_cell_xyz

subroutine get_atoms_cell_charge(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%charge
end subroutine get_atoms_cell_charge

subroutine set_atoms_cell_charge(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%charge = new_value
end subroutine set_atoms_cell_charge

subroutine get_atoms_cell_ddlab(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%ddlab
end subroutine get_atoms_cell_ddlab

subroutine set_atoms_cell_ddlab(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(in) :: new_value
	obj_var%ddlab = new_value
end subroutine set_atoms_cell_ddlab

subroutine get_atoms_cell_noms(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%noms
end subroutine get_atoms_cell_noms

subroutine set_atoms_cell_noms(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(in) :: new_value
	obj_var%noms = new_value
end subroutine set_atoms_cell_noms

subroutine get_atoms_cell_moment(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%moment
end subroutine get_atoms_cell_moment

subroutine set_atoms_cell_moment(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%moment = new_value
end subroutine set_atoms_cell_moment

subroutine get_atoms_cell_neighb_atom(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	integer,              dimension( :,:), allocatable, intent(out) :: output_value
	output_value = obj_var%neighb_atom
end subroutine get_atoms_cell_neighb_atom

subroutine set_atoms_cell_neighb_atom(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer,              dimension( :,:), allocatable, intent(in) :: new_value
	obj_var%neighb_atom = new_value
end subroutine set_atoms_cell_neighb_atom

function get_atoms_cell_nat(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer :: get_atoms_cell_nat
	get_atoms_cell_nat = obj_var%nat
end function get_atoms_cell_nat

subroutine set_atoms_cell_nat(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nat = new_value
end subroutine set_atoms_cell_nat

subroutine get_atoms_cell_var_free(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(out) :: output_value
	output_value = obj_var%var_free
end subroutine get_atoms_cell_var_free

subroutine set_atoms_cell_var_free(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%var_free = new_value
end subroutine set_atoms_cell_var_free

subroutine get_atoms_cell_trans(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),      dimension(:, :,:), allocatable, intent(out) :: output_value
	output_value = obj_var%trans
end subroutine get_atoms_cell_trans

subroutine set_atoms_cell_trans(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),      dimension(:, :,:), allocatable, intent(in) :: new_value
	obj_var%trans = new_value
end subroutine set_atoms_cell_trans

subroutine get_atoms_cell_neighb(obj_var, output_value)
	type (Atoms_Cell_Type) :: obj_var
	integer,                 dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%neighb
end subroutine get_atoms_cell_neighb

subroutine set_atoms_cell_neighb(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer,                 dimension(:), allocatable, intent(in) :: new_value
	obj_var%neighb = new_value
end subroutine set_atoms_cell_neighb

subroutine Atoms_Cell_Type_ctor(Atoms_Cell_Type_param, distance_param, ndist_param, ddist_param, xyz_param, charge_param, ddlab_param, noms_param, moment_param, neighb_atom_param, nat_param, var_free_param, trans_param, neighb_param)
	type (Atoms_Cell_Type) :: Atoms_Cell_Type_param
	real(kind=cp),        dimension( :,:), allocatable, intent(in) :: distance_param
	integer, intent(in) :: ndist_param
	real(kind=cp),           dimension(:), allocatable, intent(in) :: ddist_param
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: xyz_param
	real(kind=cp),           dimension(:), allocatable, intent(in) :: charge_param
	character (len=20),      dimension(:), allocatable, intent(in) :: ddlab_param
	character (len=20),      dimension(:), allocatable, intent(in) :: noms_param
	real(kind=cp),           dimension(:), allocatable, intent(in) :: moment_param
	integer,              dimension( :,:), allocatable, intent(in) :: neighb_atom_param
	integer, intent(in) :: nat_param
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: var_free_param
	real(kind=cp),      dimension(:, :,:), allocatable, intent(in) :: trans_param
	integer,                 dimension(:), allocatable, intent(in) :: neighb_param
	Atoms_Cell_Type_param%distance = distance_param
	Atoms_Cell_Type_param%ndist = ndist_param
	Atoms_Cell_Type_param%ddist = ddist_param
	Atoms_Cell_Type_param%xyz = xyz_param
	Atoms_Cell_Type_param%charge = charge_param
	Atoms_Cell_Type_param%ddlab = ddlab_param
	Atoms_Cell_Type_param%noms = noms_param
	Atoms_Cell_Type_param%moment = moment_param
	Atoms_Cell_Type_param%neighb_atom = neighb_atom_param
	Atoms_Cell_Type_param%nat = nat_param
	Atoms_Cell_Type_param%var_free = var_free_param
	Atoms_Cell_Type_param%trans = trans_param
	Atoms_Cell_Type_param%neighb = neighb_param
end subroutine Atoms_Cell_Type_ctor
