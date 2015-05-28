function get_atom_LOcc(obj_var)
	type (Atom_Type) :: obj_var
	integer :: get_atom_LOcc
	get_atom_LOcc = obj_var%LOcc
end function get_atom_LOcc

subroutine set_atom_LOcc(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LOcc = new_value
end subroutine set_atom_LOcc

subroutine get_atom_LVarF(obj_var, output_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(25), intent(out) :: output_value
	output_value = obj_var%LVarF
end subroutine get_atom_LVarF

subroutine set_atom_LVarF(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(25), intent(in) :: new_value
	obj_var%LVarF = new_value
end subroutine set_atom_LVarF

subroutine get_atom_Utype(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(out) :: output_value
	output_value = obj_var%Utype
end subroutine get_atom_Utype

subroutine set_atom_Utype(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%Utype = new_value
end subroutine set_atom_Utype

function get_atom_MBiso(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_MBiso
	get_atom_MBiso = obj_var%MBiso
end function get_atom_MBiso

subroutine set_atom_MBiso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MBiso = new_value
end subroutine set_atom_MBiso

function get_atom_Occ(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Occ
	get_atom_Occ = obj_var%Occ
end function get_atom_Occ

subroutine set_atom_Occ(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ = new_value
end subroutine set_atom_Occ

function get_atom_Charge(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Charge
	get_atom_Charge = obj_var%Charge
end function get_atom_Charge

subroutine set_atom_Charge(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Charge = new_value
end subroutine set_atom_Charge

subroutine get_atom_SfacSymb(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(out) :: output_value
	output_value = obj_var%SfacSymb
end subroutine get_atom_SfacSymb

subroutine set_atom_SfacSymb(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%SfacSymb = new_value
end subroutine set_atom_SfacSymb

subroutine get_atom_Lab(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=20), intent(out) :: output_value
	output_value = obj_var%Lab
end subroutine get_atom_Lab

subroutine set_atom_Lab(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine set_atom_Lab

subroutine get_atom_MVarF(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(out) :: output_value
	output_value = obj_var%MVarF
end subroutine get_atom_MVarF

subroutine set_atom_MVarF(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(in) :: new_value
	obj_var%MVarF = new_value
end subroutine set_atom_MVarF

function get_atom_Moment(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Moment
	get_atom_Moment = obj_var%Moment
end function get_atom_Moment

subroutine set_atom_Moment(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Moment = new_value
end subroutine set_atom_Moment

subroutine get_atom_LU(obj_var, output_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(6), intent(out) :: output_value
	output_value = obj_var%LU
end subroutine get_atom_LU

subroutine set_atom_LU(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(6), intent(in) :: new_value
	obj_var%LU = new_value
end subroutine set_atom_LU

function get_atom_MOcc(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_MOcc
	get_atom_MOcc = obj_var%MOcc
end function get_atom_MOcc

subroutine set_atom_MOcc(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MOcc = new_value
end subroutine set_atom_MOcc

function get_atom_Active(obj_var)
	type (Atom_Type) :: obj_var
	logical :: get_atom_Active
	get_atom_Active = obj_var%Active
end function get_atom_Active

subroutine set_atom_Active(obj_var, new_value)
	type (Atom_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Active = new_value
end subroutine set_atom_Active

function get_atom_Mult(obj_var)
	type (Atom_Type) :: obj_var
	integer :: get_atom_Mult
	get_atom_Mult = obj_var%Mult
end function get_atom_Mult

subroutine set_atom_Mult(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine set_atom_Mult

subroutine get_atom_X_Std(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%X_Std
end subroutine get_atom_X_Std

subroutine set_atom_X_Std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X_Std = new_value
end subroutine set_atom_X_Std

subroutine get_atom_U_std(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%U_std
end subroutine get_atom_U_std

subroutine set_atom_U_std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U_std = new_value
end subroutine set_atom_U_std

function get_atom_NVar(obj_var)
	type (Atom_Type) :: obj_var
	integer :: get_atom_NVar
	get_atom_NVar = obj_var%NVar
end function get_atom_NVar

subroutine set_atom_NVar(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NVar = new_value
end subroutine set_atom_NVar

subroutine get_atom_wyck(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=1), intent(out) :: output_value
	output_value = obj_var%wyck
end subroutine get_atom_wyck

subroutine set_atom_wyck(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%wyck = new_value
end subroutine set_atom_wyck

function get_atom_Biso_std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Biso_std
	get_atom_Biso_std = obj_var%Biso_std
end function get_atom_Biso_std

subroutine set_atom_Biso_std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso_std = new_value
end subroutine set_atom_Biso_std

function get_atom_LBiso(obj_var)
	type (Atom_Type) :: obj_var
	integer :: get_atom_LBiso
	get_atom_LBiso = obj_var%LBiso
end function get_atom_LBiso

subroutine set_atom_LBiso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LBiso = new_value
end subroutine set_atom_LBiso

function get_atom_Biso(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Biso
	get_atom_Biso = obj_var%Biso
end function get_atom_Biso

subroutine set_atom_Biso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso = new_value
end subroutine set_atom_Biso

subroutine get_atom_VarF(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(out) :: output_value
	output_value = obj_var%VarF
end subroutine get_atom_VarF

subroutine set_atom_VarF(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(in) :: new_value
	obj_var%VarF = new_value
end subroutine set_atom_VarF

subroutine get_atom_U(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%U
end subroutine get_atom_U

subroutine set_atom_U(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U = new_value
end subroutine set_atom_U

function get_atom_Occ_Std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Occ_Std
	get_atom_Occ_Std = obj_var%Occ_Std
end function get_atom_Occ_Std

subroutine set_atom_Occ_Std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ_Std = new_value
end subroutine set_atom_Occ_Std

subroutine get_atom_X(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%X
end subroutine get_atom_X

subroutine set_atom_X(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X = new_value
end subroutine set_atom_X

function get_atom_Z(obj_var)
	type (Atom_Type) :: obj_var
	integer :: get_atom_Z
	get_atom_Z = obj_var%Z
end function get_atom_Z

subroutine set_atom_Z(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Z = new_value
end subroutine set_atom_Z

subroutine get_atom_MU(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%MU
end subroutine get_atom_MU

subroutine set_atom_MU(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%MU = new_value
end subroutine set_atom_MU

subroutine get_atom_LX(obj_var, output_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(3), intent(out) :: output_value
	output_value = obj_var%LX
end subroutine get_atom_LX

subroutine set_atom_LX(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%LX = new_value
end subroutine set_atom_LX

subroutine get_atom_ChemSymb(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=2), intent(out) :: output_value
	output_value = obj_var%ChemSymb
end subroutine get_atom_ChemSymb

subroutine set_atom_ChemSymb(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=2), intent(in) :: new_value
	obj_var%ChemSymb = new_value
end subroutine set_atom_ChemSymb

subroutine get_atom_Ind(obj_var, output_value)
	type (Atom_Type) :: obj_var
	integer, dimension(5), intent(out) :: output_value
	output_value = obj_var%Ind
end subroutine get_atom_Ind

subroutine set_atom_Ind(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, dimension(5), intent(in) :: new_value
	obj_var%Ind = new_value
end subroutine set_atom_Ind

subroutine get_atom_ThType(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=5), intent(out) :: output_value
	output_value = obj_var%ThType
end subroutine get_atom_ThType

subroutine set_atom_ThType(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=5), intent(in) :: new_value
	obj_var%ThType = new_value
end subroutine set_atom_ThType

subroutine get_atom_AtmInfo(obj_var, output_value)
	type (Atom_Type) :: obj_var
	character(len=40), intent(out) :: output_value
	output_value = obj_var%AtmInfo
end subroutine get_atom_AtmInfo

subroutine set_atom_AtmInfo(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%AtmInfo = new_value
end subroutine set_atom_AtmInfo

function get_atom_Ueq(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: get_atom_Ueq
	get_atom_Ueq = obj_var%Ueq
end function get_atom_Ueq

subroutine set_atom_Ueq(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Ueq = new_value
end subroutine set_atom_Ueq

subroutine get_atom_MX(obj_var, output_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%MX
end subroutine get_atom_MX

subroutine set_atom_MX(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%MX = new_value
end subroutine set_atom_MX

subroutine Atom_Type_ctor(Atom_Type_param, LOcc_param, LVarF_param, Utype_param, MBiso_param, Occ_param, Charge_param, SfacSymb_param, Lab_param, MVarF_param, Moment_param, LU_param, MOcc_param, Active_param, Mult_param, X_Std_param, U_std_param, NVar_param, wyck_param, Biso_std_param, LBiso_param, Biso_param, VarF_param, U_param, Occ_Std_param, X_param, Z_param, MU_param, LX_param, ChemSymb_param, Ind_param, ThType_param, AtmInfo_param, Ueq_param, MX_param)
	type (Atom_Type) :: Atom_Type_param
	integer, intent(in) :: LOcc_param
	integer,      dimension(25), intent(in) :: LVarF_param
	character(len=4), intent(in) :: Utype_param
	real(kind=cp), intent(in) :: MBiso_param
	real(kind=cp), intent(in) :: Occ_param
	real(kind=cp), intent(in) :: Charge_param
	character(len=4), intent(in) :: SfacSymb_param
	character(len=20), intent(in) :: Lab_param
	real(kind=cp),dimension(25), intent(in) :: MVarF_param
	real(kind=cp), intent(in) :: Moment_param
	integer,      dimension(6), intent(in) :: LU_param
	real(kind=cp), intent(in) :: MOcc_param
	logical, intent(in) :: Active_param
	integer, intent(in) :: Mult_param
	real(kind=cp),dimension(3), intent(in) :: X_Std_param
	real(kind=cp),dimension(6), intent(in) :: U_std_param
	integer, intent(in) :: NVar_param
	character(len=1), intent(in) :: wyck_param
	real(kind=cp), intent(in) :: Biso_std_param
	integer, intent(in) :: LBiso_param
	real(kind=cp), intent(in) :: Biso_param
	real(kind=cp),dimension(25), intent(in) :: VarF_param
	real(kind=cp),dimension(6), intent(in) :: U_param
	real(kind=cp), intent(in) :: Occ_Std_param
	real(kind=cp),dimension(3), intent(in) :: X_param
	integer, intent(in) :: Z_param
	real(kind=cp),dimension(6), intent(in) :: MU_param
	integer,      dimension(3), intent(in) :: LX_param
	character(len=2), intent(in) :: ChemSymb_param
	integer, dimension(5), intent(in) :: Ind_param
	character(len=5), intent(in) :: ThType_param
	character(len=40), intent(in) :: AtmInfo_param
	real(kind=cp), intent(in) :: Ueq_param
	real(kind=cp),dimension(3), intent(in) :: MX_param
	Atom_Type_param%LOcc = LOcc_param
	Atom_Type_param%LVarF = LVarF_param
	Atom_Type_param%Utype = Utype_param
	Atom_Type_param%MBiso = MBiso_param
	Atom_Type_param%Occ = Occ_param
	Atom_Type_param%Charge = Charge_param
	Atom_Type_param%SfacSymb = SfacSymb_param
	Atom_Type_param%Lab = Lab_param
	Atom_Type_param%MVarF = MVarF_param
	Atom_Type_param%Moment = Moment_param
	Atom_Type_param%LU = LU_param
	Atom_Type_param%MOcc = MOcc_param
	Atom_Type_param%Active = Active_param
	Atom_Type_param%Mult = Mult_param
	Atom_Type_param%X_Std = X_Std_param
	Atom_Type_param%U_std = U_std_param
	Atom_Type_param%NVar = NVar_param
	Atom_Type_param%wyck = wyck_param
	Atom_Type_param%Biso_std = Biso_std_param
	Atom_Type_param%LBiso = LBiso_param
	Atom_Type_param%Biso = Biso_param
	Atom_Type_param%VarF = VarF_param
	Atom_Type_param%U = U_param
	Atom_Type_param%Occ_Std = Occ_Std_param
	Atom_Type_param%X = X_param
	Atom_Type_param%Z = Z_param
	Atom_Type_param%MU = MU_param
	Atom_Type_param%LX = LX_param
	Atom_Type_param%ChemSymb = ChemSymb_param
	Atom_Type_param%Ind = Ind_param
	Atom_Type_param%ThType = ThType_param
	Atom_Type_param%AtmInfo = AtmInfo_param
	Atom_Type_param%Ueq = Ueq_param
	Atom_Type_param%MX = MX_param
end subroutine Atom_Type_ctor
