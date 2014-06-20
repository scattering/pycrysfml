function getLOcc(obj_var)
	type (Atom_Type) :: obj_var
	integer :: getLOcc
	getLOcc = obj_var%LOcc
end function getLOcc

subroutine setLOcc(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LOcc = new_value
end subroutine setLOcc

function getUtype(obj_var)
	type (Atom_Type) :: obj_var
	character(len=4) :: getUtype
	getUtype = obj_var%Utype
end function getUtype

subroutine setUtype(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%Utype = new_value
end subroutine setUtype

function getMBiso(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getMBiso
	getMBiso = obj_var%MBiso
end function getMBiso

subroutine setMBiso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MBiso = new_value
end subroutine setMBiso

function getOcc(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getOcc
	getOcc = obj_var%Occ
end function getOcc

subroutine setOcc(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ = new_value
end subroutine setOcc

function getCharge(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getCharge
	getCharge = obj_var%Charge
end function getCharge

subroutine setCharge(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Charge = new_value
end subroutine setCharge

function getSfacSymb(obj_var)
	type (Atom_Type) :: obj_var
	character(len=4) :: getSfacSymb
	getSfacSymb = obj_var%SfacSymb
end function getSfacSymb

subroutine setSfacSymb(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%SfacSymb = new_value
end subroutine setSfacSymb

function getLab(obj_var)
	type (Atom_Type) :: obj_var
	character(len=20) :: getLab
	getLab = obj_var%Lab
end function getLab

subroutine setLab(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine setLab

function getMoment(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getMoment
	getMoment = obj_var%Moment
end function getMoment

subroutine setMoment(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Moment = new_value
end subroutine setMoment

function getLU(obj_var)
	type (Atom_Type) :: obj_var
	integer,      dimension(6) :: getLU
	getLU = obj_var%LU
end function getLU

subroutine setLU(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(6), intent(in) :: new_value
	obj_var%LU = new_value
end subroutine setLU

function getMOcc(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getMOcc
	getMOcc = obj_var%MOcc
end function getMOcc

subroutine setMOcc(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MOcc = new_value
end subroutine setMOcc

function getActive(obj_var)
	type (Atom_Type) :: obj_var
	logical :: getActive
	getActive = obj_var%Active
end function getActive

subroutine setActive(obj_var, new_value)
	type (Atom_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Active = new_value
end subroutine setActive

function getMult(obj_var)
	type (Atom_Type) :: obj_var
	integer :: getMult
	getMult = obj_var%Mult
end function getMult

subroutine setMult(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine setMult

function getX_Std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3) :: getX_Std
	getX_Std = obj_var%X_Std
end function getX_Std

subroutine setX_Std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%X_Std = new_value
end subroutine setX_Std

function getU_std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6) :: getU_std
	getU_std = obj_var%U_std
end function getU_std

subroutine setU_std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6), intent(in) :: new_value
	obj_var%U_std = new_value
end subroutine setU_std

function getNVar(obj_var)
	type (Atom_Type) :: obj_var
	integer :: getNVar
	getNVar = obj_var%NVar
end function getNVar

subroutine setNVar(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NVar = new_value
end subroutine setNVar

function getwyck(obj_var)
	type (Atom_Type) :: obj_var
	character(len=1) :: getwyck
	getwyck = obj_var%wyck
end function getwyck

subroutine setwyck(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%wyck = new_value
end subroutine setwyck

function getBiso_std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getBiso_std
	getBiso_std = obj_var%Biso_std
end function getBiso_std

subroutine setBiso_std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso_std = new_value
end subroutine setBiso_std

function getLBiso(obj_var)
	type (Atom_Type) :: obj_var
	integer :: getLBiso
	getLBiso = obj_var%LBiso
end function getLBiso

subroutine setLBiso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LBiso = new_value
end subroutine setLBiso

function getBiso(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getBiso
	getBiso = obj_var%Biso
end function getBiso

subroutine setBiso(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso = new_value
end subroutine setBiso

function getVarF(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(10) :: getVarF
	getVarF = obj_var%VarF
end function getVarF

subroutine setVarF(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(10), intent(in) :: new_value
	obj_var%VarF = new_value
end subroutine setVarF

function getU(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6) :: getU
	getU = obj_var%U
end function getU

subroutine setU(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6), intent(in) :: new_value
	obj_var%U = new_value
end subroutine setU

function getOcc_Std(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getOcc_Std
	getOcc_Std = obj_var%Occ_Std
end function getOcc_Std

subroutine setOcc_Std(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ_Std = new_value
end subroutine setOcc_Std

function getX(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3) :: getX
	getX = obj_var%X
end function getX

subroutine setX(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%X = new_value
end subroutine setX

function getZ(obj_var)
	type (Atom_Type) :: obj_var
	integer :: getZ
	getZ = obj_var%Z
end function getZ

subroutine setZ(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Z = new_value
end subroutine setZ

function getMU(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6) :: getMU
	getMU = obj_var%MU
end function getMU

subroutine setMU(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(6), intent(in) :: new_value
	obj_var%MU = new_value
end subroutine setMU

function getLX(obj_var)
	type (Atom_Type) :: obj_var
	integer,      dimension(3) :: getLX
	getLX = obj_var%LX
end function getLX

subroutine setLX(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%LX = new_value
end subroutine setLX

function getChemSymb(obj_var)
	type (Atom_Type) :: obj_var
	character(len=2) :: getChemSymb
	getChemSymb = obj_var%ChemSymb
end function getChemSymb

subroutine setChemSymb(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=2), intent(in) :: new_value
	obj_var%ChemSymb = new_value
end subroutine setChemSymb

function getInd(obj_var)
	type (Atom_Type) :: obj_var
	integer, dimension(5) :: getInd
	getInd = obj_var%Ind
end function getInd

subroutine setInd(obj_var, new_value)
	type (Atom_Type) :: obj_var
	integer, dimension(5), intent(in) :: new_value
	obj_var%Ind = new_value
end subroutine setInd

function getThType(obj_var)
	type (Atom_Type) :: obj_var
	character(len=5) :: getThType
	getThType = obj_var%ThType
end function getThType

subroutine setThType(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=5), intent(in) :: new_value
	obj_var%ThType = new_value
end subroutine setThType

function getAtmInfo(obj_var)
	type (Atom_Type) :: obj_var
	character(len=40) :: getAtmInfo
	getAtmInfo = obj_var%AtmInfo
end function getAtmInfo

subroutine setAtmInfo(obj_var, new_value)
	type (Atom_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%AtmInfo = new_value
end subroutine setAtmInfo

function getUeq(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp) :: getUeq
	getUeq = obj_var%Ueq
end function getUeq

subroutine setUeq(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Ueq = new_value
end subroutine setUeq

function getMX(obj_var)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3) :: getMX
	getMX = obj_var%MX
end function getMX

subroutine setMX(obj_var, new_value)
	type (Atom_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%MX = new_value
end subroutine setMX

subroutine Atom_Type_ctor(Atom_Type_param, LOcc_param, Utype_param, MBiso_param, Occ_param, Charge_param, SfacSymb_param, Lab_param, Moment_param, LU_param, MOcc_param, Active_param, Mult_param, X_Std_param, U_std_param, NVar_param, wyck_param, Biso_std_param, LBiso_param, Biso_param, VarF_param, U_param, Occ_Std_param, X_param, Z_param, MU_param, LX_param, ChemSymb_param, Ind_param, ThType_param, AtmInfo_param, Ueq_param, MX_param)
	type (Atom_Type) :: Atom_Type_param
	integer, intent(in) :: LOcc_param
	character(len=4), intent(in) :: Utype_param
	real(kind=cp), intent(in) :: MBiso_param
	real(kind=cp), intent(in) :: Occ_param
	real(kind=cp), intent(in) :: Charge_param
	character(len=4), intent(in) :: SfacSymb_param
	character(len=20), intent(in) :: Lab_param
	real(kind=cp), intent(in) :: Moment_param
	integer,      dimension(6), intent(in) :: LU_param
	real(kind=cp), intent(in) :: MOcc_param
	logical, intent(in) :: Active_param
	integer, intent(in) :: Mult_param
	real(kind=cp), dimension(3), intent(in) :: X_Std_param
	real(kind=cp), dimension(6), intent(in) :: U_std_param
	integer, intent(in) :: NVar_param
	character(len=1), intent(in) :: wyck_param
	real(kind=cp), intent(in) :: Biso_std_param
	integer, intent(in) :: LBiso_param
	real(kind=cp), intent(in) :: Biso_param
	real(kind=cp), dimension(10), intent(in) :: VarF_param
	real(kind=cp), dimension(6), intent(in) :: U_param
	real(kind=cp), intent(in) :: Occ_Std_param
	real(kind=cp), dimension(3), intent(in) :: X_param
	integer, intent(in) :: Z_param
	real(kind=cp), dimension(6), intent(in) :: MU_param
	integer,      dimension(3), intent(in) :: LX_param
	character(len=2), intent(in) :: ChemSymb_param
	integer, dimension(5), intent(in) :: Ind_param
	character(len=5), intent(in) :: ThType_param
	character(len=40), intent(in) :: AtmInfo_param
	real(kind=cp), intent(in) :: Ueq_param
	real(kind=cp), dimension(3), intent(in) :: MX_param
	Atom_Type_param%LOcc = LOcc_param
	Atom_Type_param%Utype = Utype_param
	Atom_Type_param%MBiso = MBiso_param
	Atom_Type_param%Occ = Occ_param
	Atom_Type_param%Charge = Charge_param
	Atom_Type_param%SfacSymb = SfacSymb_param
	Atom_Type_param%Lab = Lab_param
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
