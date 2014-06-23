function get_atoms_cell_distance(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),        dimension( :,:), allocatable :: getdistance
	getdistance = obj_var%distance
end function get_atoms_cell_distance

subroutine set_atoms_cell_distance(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),        dimension( :,:), allocatable, intent(in) :: new_value
	obj_var%distance = new_value
end subroutine set_atoms_cell_distance

function get_atoms_cell_ndist(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer :: getndist
	getndist = obj_var%ndist
end function get_atoms_cell_ndist

subroutine set_atoms_cell_ndist(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%ndist = new_value
end subroutine set_atoms_cell_ndist

function get_atoms_cell_ddist(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable :: getddist
	getddist = obj_var%ddist
end function get_atoms_cell_ddist

subroutine set_atoms_cell_ddist(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%ddist = new_value
end subroutine set_atoms_cell_ddist

function get_atoms_cell_xyz(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable :: getxyz
	getxyz = obj_var%xyz
end function get_atoms_cell_xyz

subroutine set_atoms_cell_xyz(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%xyz = new_value
end subroutine set_atoms_cell_xyz

function get_atoms_cell_charge(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable :: getcharge
	getcharge = obj_var%charge
end function get_atoms_cell_charge

subroutine set_atoms_cell_charge(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%charge = new_value
end subroutine set_atoms_cell_charge

function get_atoms_cell_ddlab(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable :: getddlab
	getddlab = obj_var%ddlab
end function get_atoms_cell_ddlab

subroutine set_atoms_cell_ddlab(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(in) :: new_value
	obj_var%ddlab = new_value
end subroutine set_atoms_cell_ddlab

function get_atoms_cell_noms(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable :: getnoms
	getnoms = obj_var%noms
end function get_atoms_cell_noms

subroutine set_atoms_cell_noms(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	character (len=20),      dimension(:), allocatable, intent(in) :: new_value
	obj_var%noms = new_value
end subroutine set_atoms_cell_noms

function get_atoms_cell_moment(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable :: getmoment
	getmoment = obj_var%moment
end function get_atoms_cell_moment

subroutine set_atoms_cell_moment(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),           dimension(:), allocatable, intent(in) :: new_value
	obj_var%moment = new_value
end subroutine set_atoms_cell_moment

function get_atoms_cell_neighb_atom(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer,              dimension( :,:), allocatable :: getneighb_atom
	getneighb_atom = obj_var%neighb_atom
end function get_atoms_cell_neighb_atom

subroutine set_atoms_cell_neighb_atom(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer,              dimension( :,:), allocatable, intent(in) :: new_value
	obj_var%neighb_atom = new_value
end subroutine set_atoms_cell_neighb_atom

function get_atoms_cell_nat(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer :: getnat
	getnat = obj_var%nat
end function get_atoms_cell_nat

subroutine set_atoms_cell_nat(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nat = new_value
end subroutine set_atoms_cell_nat

function get_atoms_cell_var_free(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable :: getvar_free
	getvar_free = obj_var%var_free
end function get_atoms_cell_var_free

subroutine set_atoms_cell_var_free(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),         dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%var_free = new_value
end subroutine set_atoms_cell_var_free

function get_atoms_cell_trans(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),      dimension(:, :,:), allocatable :: gettrans
	gettrans = obj_var%trans
end function get_atoms_cell_trans

subroutine set_atoms_cell_trans(obj_var, new_value)
	type (Atoms_Cell_Type) :: obj_var
	real(kind=cp),      dimension(:, :,:), allocatable, intent(in) :: new_value
	obj_var%trans = new_value
end subroutine set_atoms_cell_trans

function get_atoms_cell_neighb(obj_var)
	type (Atoms_Cell_Type) :: obj_var
	integer,                 dimension(:), allocatable :: getneighb
	getneighb = obj_var%neighb
end function get_atoms_cell_neighb

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
