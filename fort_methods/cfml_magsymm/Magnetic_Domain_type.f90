function getChir(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: getChir
	getChir = obj_var%Chir
end function getChir

subroutine setChir(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Chir = new_value
end subroutine setChir

function getpop_std(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24) :: getpop_std
	getpop_std = obj_var%pop_std
end function getpop_std

subroutine setpop_std(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%pop_std = new_value
end subroutine setpop_std

function getMpop(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24) :: getMpop
	getMpop = obj_var%Mpop
end function getMpop

subroutine setMpop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%Mpop = new_value
end subroutine setMpop

function getnd(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	integer :: getnd
	getnd = obj_var%nd
end function getnd

subroutine setnd(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nd = new_value
end subroutine setnd

function getLpop(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	integer      , dimension (2,24) :: getLpop
	getLpop = obj_var%Lpop
end function getLpop

subroutine setLpop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer      , dimension (2,24), intent(in) :: new_value
	obj_var%Lpop = new_value
end subroutine setLpop

function getpop(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24) :: getpop
	getpop = obj_var%pop
end function getpop

subroutine setpop(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (2,24), intent(in) :: new_value
	obj_var%pop = new_value
end subroutine setpop

function getTwin(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: getTwin
	getTwin = obj_var%Twin
end function getTwin

subroutine setTwin(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Twin = new_value
end subroutine setTwin

function getLab(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	character(len=10),dimension (2,24) :: getLab
	getLab = obj_var%Lab
end function getLab

subroutine setLab(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	character(len=10),dimension (2,24), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine setLab

function getDMat(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	integer,dimension(3,3,24) :: getDMat
	getDMat = obj_var%DMat
end function getDMat

subroutine setDMat(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	integer,dimension(3,3,24), intent(in) :: new_value
	obj_var%DMat = new_value
end subroutine setDMat

function getDt(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (3,24) :: getDt
	getDt = obj_var%Dt
end function getDt

subroutine setDt(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	real(kind=cp), dimension (3,24), intent(in) :: new_value
	obj_var%Dt = new_value
end subroutine setDt

function gettrans(obj_var)
	type (Magnetic_Domain_type) :: obj_var
	logical :: gettrans
	gettrans = obj_var%trans
end function gettrans

subroutine settrans(obj_var, new_value)
	type (Magnetic_Domain_type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%trans = new_value
end subroutine settrans

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
