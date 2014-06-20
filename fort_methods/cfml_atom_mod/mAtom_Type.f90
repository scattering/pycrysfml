function getmmphas(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12) :: getmmphas
	getmmphas = obj_var%mmphas
end function getmmphas

subroutine setmmphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mmphas = new_value
end subroutine setmmphas

function getLOcc(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getLOcc
	getLOcc = obj_var%LOcc
end function getLOcc

subroutine setLOcc(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LOcc = new_value
end subroutine setLOcc

function getSkI_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSkI_std
	getSkI_std = obj_var%SkI_std
end function getSkI_std

subroutine setSkI_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkI_std = new_value
end subroutine setSkI_std

function getUtype(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=4) :: getUtype
	getUtype = obj_var%Utype
end function getUtype

subroutine setUtype(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%Utype = new_value
end subroutine setUtype

function getMBiso(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getMBiso
	getMBiso = obj_var%MBiso
end function getMBiso

subroutine setMBiso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MBiso = new_value
end subroutine setMBiso

function getOcc(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getOcc
	getOcc = obj_var%Occ
end function getOcc

subroutine setOcc(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ = new_value
end subroutine setOcc

function getlbas(obj_var)
	type (mAtom_Type) :: obj_var
	integer,dimension(12,12) :: getlbas
	getlbas = obj_var%lbas
end function getlbas

subroutine setlbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12,12), intent(in) :: new_value
	obj_var%lbas = new_value
end subroutine setlbas

function getCharge(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getCharge
	getCharge = obj_var%Charge
end function getCharge

subroutine setCharge(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Charge = new_value
end subroutine setCharge

function getSfacSymb(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=4) :: getSfacSymb
	getSfacSymb = obj_var%SfacSymb
end function getSfacSymb

subroutine setSfacSymb(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%SfacSymb = new_value
end subroutine setSfacSymb

function getimat(obj_var)
	type (mAtom_Type) :: obj_var
	integer,      dimension(12) :: getimat
	getimat = obj_var%imat
end function getimat

subroutine setimat(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(12), intent(in) :: new_value
	obj_var%imat = new_value
end subroutine setimat

function getlmphas(obj_var)
	type (mAtom_Type) :: obj_var
	integer,dimension(12) :: getlmphas
	getlmphas = obj_var%lmphas
end function getlmphas

subroutine setlmphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12), intent(in) :: new_value
	obj_var%lmphas = new_value
end subroutine setlmphas

function getmphas(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12) :: getmphas
	getmphas = obj_var%mphas
end function getmphas

subroutine setmphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mphas = new_value
end subroutine setmphas

function getSpher_SkI(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSpher_SkI
	getSpher_SkI = obj_var%Spher_SkI
end function getSpher_SkI

subroutine setSpher_SkI(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkI = new_value
end subroutine setSpher_SkI

function getSkR(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSkR
	getSkR = obj_var%SkR
end function getSkR

subroutine setSkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkR = new_value
end subroutine setSkR

function getlskr(obj_var)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12) :: getlskr
	getlskr = obj_var%lskr
end function getlskr

subroutine setlskr(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(in) :: new_value
	obj_var%lskr = new_value
end subroutine setlskr

function getLab(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=10) :: getLab
	getLab = obj_var%Lab
end function getLab

subroutine setLab(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=10), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine setLab

function getMoment(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getMoment
	getMoment = obj_var%Moment
end function getMoment

subroutine setMoment(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Moment = new_value
end subroutine setMoment

function getLU(obj_var)
	type (mAtom_Type) :: obj_var
	integer,      dimension(6) :: getLU
	getLU = obj_var%LU
end function getLU

subroutine setLU(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(6), intent(in) :: new_value
	obj_var%LU = new_value
end subroutine setLU

function getMOcc(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getMOcc
	getMOcc = obj_var%MOcc
end function getMOcc

subroutine setMOcc(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MOcc = new_value
end subroutine setMOcc

function getActive(obj_var)
	type (mAtom_Type) :: obj_var
	logical :: getActive
	getActive = obj_var%Active
end function getActive

subroutine setActive(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Active = new_value
end subroutine setActive

function getSkI(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSkI
	getSkI = obj_var%SkI
end function getSkI

subroutine setSkI(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkI = new_value
end subroutine setSkI

function getMult(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getMult
	getMult = obj_var%Mult
end function getMult

subroutine setMult(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine setMult

function getSpher_SkR_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSpher_SkR_std
	getSpher_SkR_std = obj_var%Spher_SkR_std
end function getSpher_SkR_std

subroutine setSpher_SkR_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkR_std = new_value
end subroutine setSpher_SkR_std

function getX_Std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3) :: getX_Std
	getX_Std = obj_var%X_Std
end function getX_Std

subroutine setX_Std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X_Std = new_value
end subroutine setX_Std

function getU_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6) :: getU_std
	getU_std = obj_var%U_std
end function getU_std

subroutine setU_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U_std = new_value
end subroutine setU_std

function getlski(obj_var)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12) :: getlski
	getlski = obj_var%lski
end function getlski

subroutine setlski(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(in) :: new_value
	obj_var%lski = new_value
end subroutine setlski

function getNVar(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getNVar
	getNVar = obj_var%NVar
end function getNVar

subroutine setNVar(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NVar = new_value
end subroutine setNVar

function getwyck(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=1) :: getwyck
	getwyck = obj_var%wyck
end function getwyck

subroutine setwyck(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%wyck = new_value
end subroutine setwyck

function getBiso_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getBiso_std
	getBiso_std = obj_var%Biso_std
end function getBiso_std

subroutine setBiso_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso_std = new_value
end subroutine setBiso_std

function getLBiso(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getLBiso
	getLBiso = obj_var%LBiso
end function getLBiso

subroutine setLBiso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LBiso = new_value
end subroutine setLBiso

function getmphas_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12) :: getmphas_std
	getmphas_std = obj_var%mphas_std
end function getmphas_std

subroutine setmphas_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mphas_std = new_value
end subroutine setmphas_std

function getBiso(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getBiso
	getBiso = obj_var%Biso
end function getBiso

subroutine setBiso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso = new_value
end subroutine setBiso

function getVarF(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(10) :: getVarF
	getVarF = obj_var%VarF
end function getVarF

subroutine setVarF(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(10), intent(in) :: new_value
	obj_var%VarF = new_value
end subroutine setVarF

function getU(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6) :: getU
	getU = obj_var%U
end function getU

subroutine setU(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U = new_value
end subroutine setU

function getOcc_Std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getOcc_Std
	getOcc_Std = obj_var%Occ_Std
end function getOcc_Std

subroutine setOcc_Std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ_Std = new_value
end subroutine setOcc_Std

function getX(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3) :: getX
	getX = obj_var%X
end function getX

subroutine setX(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X = new_value
end subroutine setX

function getZ(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getZ
	getZ = obj_var%Z
end function getZ

subroutine setZ(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Z = new_value
end subroutine setZ

function getnvk(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: getnvk
	getnvk = obj_var%nvk
end function getnvk

subroutine setnvk(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nvk = new_value
end subroutine setnvk

function getmbas(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12) :: getmbas
	getmbas = obj_var%mbas
end function getmbas

subroutine setmbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%mbas = new_value
end subroutine setmbas

function getSpher_SkI_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSpher_SkI_std
	getSpher_SkI_std = obj_var%Spher_SkI_std
end function getSpher_SkI_std

subroutine setSpher_SkI_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkI_std = new_value
end subroutine setSpher_SkI_std

function getSpher_SkR(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSpher_SkR
	getSpher_SkR = obj_var%Spher_SkR
end function getSpher_SkR

subroutine setSpher_SkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkR = new_value
end subroutine setSpher_SkR

function getmSki(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getmSki
	getmSki = obj_var%mSki
end function getmSki

subroutine setmSki(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%mSki = new_value
end subroutine setmSki

function getSkR_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getSkR_std
	getSkR_std = obj_var%SkR_std
end function getSkR_std

subroutine setSkR_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkR_std = new_value
end subroutine setSkR_std

function getMU(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6) :: getMU
	getMU = obj_var%MU
end function getMU

subroutine setMU(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%MU = new_value
end subroutine setMU

function getmSkR(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getmSkR
	getmSkR = obj_var%mSkR
end function getmSkR

subroutine setmSkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%mSkR = new_value
end subroutine setmSkR

function getLX(obj_var)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3) :: getLX
	getLX = obj_var%LX
end function getLX

subroutine setLX(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%LX = new_value
end subroutine setLX

function getChemSymb(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=2) :: getChemSymb
	getChemSymb = obj_var%ChemSymb
end function getChemSymb

subroutine setChemSymb(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=2), intent(in) :: new_value
	obj_var%ChemSymb = new_value
end subroutine setChemSymb

function getcbas(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12) :: getcbas
	getcbas = obj_var%cbas
end function getcbas

subroutine setcbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%cbas = new_value
end subroutine setcbas

function getInd(obj_var)
	type (mAtom_Type) :: obj_var
	integer, dimension(5) :: getInd
	getInd = obj_var%Ind
end function getInd

subroutine setInd(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, dimension(5), intent(in) :: new_value
	obj_var%Ind = new_value
end subroutine setInd

function getcbas_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12) :: getcbas_std
	getcbas_std = obj_var%cbas_std
end function getcbas_std

subroutine setcbas_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%cbas_std = new_value
end subroutine setcbas_std

function getThType(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=5) :: getThType
	getThType = obj_var%ThType
end function getThType

subroutine setThType(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=5), intent(in) :: new_value
	obj_var%ThType = new_value
end subroutine setThType

function getAtmInfo(obj_var)
	type (mAtom_Type) :: obj_var
	character(len=40) :: getAtmInfo
	getAtmInfo = obj_var%AtmInfo
end function getAtmInfo

subroutine setAtmInfo(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%AtmInfo = new_value
end subroutine setAtmInfo

function getUeq(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: getUeq
	getUeq = obj_var%Ueq
end function getUeq

subroutine setUeq(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Ueq = new_value
end subroutine setUeq

function getMX(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3) :: getMX
	getMX = obj_var%MX
end function getMX

subroutine setMX(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%MX = new_value
end subroutine setMX

subroutine mAtom_Type_ctor(mAtom_Type_param, mmphas_param, LOcc_param, SkI_std_param, Utype_param, MBiso_param, Occ_param, lbas_param, Charge_param, SfacSymb_param, imat_param, lmphas_param, mphas_param, Spher_SkI_param, SkR_param, lskr_param, Lab_param, Moment_param, LU_param, MOcc_param, Active_param, SkI_param, Mult_param, Spher_SkR_std_param, X_Std_param, U_std_param, lski_param, NVar_param, wyck_param, Biso_std_param, LBiso_param, mphas_std_param, Biso_param, VarF_param, U_param, Occ_Std_param, X_param, Z_param, nvk_param, mbas_param, Spher_SkI_std_param, Spher_SkR_param, mSki_param, SkR_std_param, MU_param, mSkR_param, LX_param, ChemSymb_param, cbas_param, Ind_param, cbas_std_param, ThType_param, AtmInfo_param, Ueq_param, MX_param)
	type (mAtom_Type) :: mAtom_Type_param
	real(kind=cp),dimension(12), intent(in) :: mmphas_param
	integer, intent(in) :: LOcc_param
	real(kind=cp),dimension(3,12), intent(in) :: SkI_std_param
	character(len=4), intent(in) :: Utype_param
	real(kind=cp), intent(in) :: MBiso_param
	real(kind=cp), intent(in) :: Occ_param
	integer,dimension(12,12), intent(in) :: lbas_param
	real(kind=cp), intent(in) :: Charge_param
	character(len=4), intent(in) :: SfacSymb_param
	integer,      dimension(12), intent(in) :: imat_param
	integer,dimension(12), intent(in) :: lmphas_param
	real(kind=cp),dimension(12), intent(in) :: mphas_param
	real(kind=cp),dimension(3,12), intent(in) :: Spher_SkI_param
	real(kind=cp),dimension(3,12), intent(in) :: SkR_param
	integer,      dimension(3,12), intent(in) :: lskr_param
	character(len=10), intent(in) :: Lab_param
	real(kind=cp), intent(in) :: Moment_param
	integer,      dimension(6), intent(in) :: LU_param
	real(kind=cp), intent(in) :: MOcc_param
	logical, intent(in) :: Active_param
	real(kind=cp),dimension(3,12), intent(in) :: SkI_param
	integer, intent(in) :: Mult_param
	real(kind=cp),dimension(3,12), intent(in) :: Spher_SkR_std_param
	real(kind=cp),dimension(3), intent(in) :: X_Std_param
	real(kind=cp),dimension(6), intent(in) :: U_std_param
	integer,      dimension(3,12), intent(in) :: lski_param
	integer, intent(in) :: NVar_param
	character(len=1), intent(in) :: wyck_param
	real(kind=cp), intent(in) :: Biso_std_param
	integer, intent(in) :: LBiso_param
	real(kind=cp),dimension(12), intent(in) :: mphas_std_param
	real(kind=cp), intent(in) :: Biso_param
	real(kind=cp),dimension(10), intent(in) :: VarF_param
	real(kind=cp),dimension(6), intent(in) :: U_param
	real(kind=cp), intent(in) :: Occ_Std_param
	real(kind=cp),dimension(3), intent(in) :: X_param
	integer, intent(in) :: Z_param
	integer, intent(in) :: nvk_param
	real(kind=cp),dimension(12,12), intent(in) :: mbas_param
	real(kind=cp),dimension(3,12), intent(in) :: Spher_SkI_std_param
	real(kind=cp),dimension(3,12), intent(in) :: Spher_SkR_param
	real(kind=cp),dimension(3,12), intent(in) :: mSki_param
	real(kind=cp),dimension(3,12), intent(in) :: SkR_std_param
	real(kind=cp),dimension(6), intent(in) :: MU_param
	real(kind=cp),dimension(3,12), intent(in) :: mSkR_param
	integer,      dimension(3), intent(in) :: LX_param
	character(len=2), intent(in) :: ChemSymb_param
	real(kind=cp),dimension(12,12), intent(in) :: cbas_param
	integer, dimension(5), intent(in) :: Ind_param
	real(kind=cp),dimension(12,12), intent(in) :: cbas_std_param
	character(len=5), intent(in) :: ThType_param
	character(len=40), intent(in) :: AtmInfo_param
	real(kind=cp), intent(in) :: Ueq_param
	real(kind=cp),dimension(3), intent(in) :: MX_param
	mAtom_Type_param%mmphas = mmphas_param
	mAtom_Type_param%LOcc = LOcc_param
	mAtom_Type_param%SkI_std = SkI_std_param
	mAtom_Type_param%Utype = Utype_param
	mAtom_Type_param%MBiso = MBiso_param
	mAtom_Type_param%Occ = Occ_param
	mAtom_Type_param%lbas = lbas_param
	mAtom_Type_param%Charge = Charge_param
	mAtom_Type_param%SfacSymb = SfacSymb_param
	mAtom_Type_param%imat = imat_param
	mAtom_Type_param%lmphas = lmphas_param
	mAtom_Type_param%mphas = mphas_param
	mAtom_Type_param%Spher_SkI = Spher_SkI_param
	mAtom_Type_param%SkR = SkR_param
	mAtom_Type_param%lskr = lskr_param
	mAtom_Type_param%Lab = Lab_param
	mAtom_Type_param%Moment = Moment_param
	mAtom_Type_param%LU = LU_param
	mAtom_Type_param%MOcc = MOcc_param
	mAtom_Type_param%Active = Active_param
	mAtom_Type_param%SkI = SkI_param
	mAtom_Type_param%Mult = Mult_param
	mAtom_Type_param%Spher_SkR_std = Spher_SkR_std_param
	mAtom_Type_param%X_Std = X_Std_param
	mAtom_Type_param%U_std = U_std_param
	mAtom_Type_param%lski = lski_param
	mAtom_Type_param%NVar = NVar_param
	mAtom_Type_param%wyck = wyck_param
	mAtom_Type_param%Biso_std = Biso_std_param
	mAtom_Type_param%LBiso = LBiso_param
	mAtom_Type_param%mphas_std = mphas_std_param
	mAtom_Type_param%Biso = Biso_param
	mAtom_Type_param%VarF = VarF_param
	mAtom_Type_param%U = U_param
	mAtom_Type_param%Occ_Std = Occ_Std_param
	mAtom_Type_param%X = X_param
	mAtom_Type_param%Z = Z_param
	mAtom_Type_param%nvk = nvk_param
	mAtom_Type_param%mbas = mbas_param
	mAtom_Type_param%Spher_SkI_std = Spher_SkI_std_param
	mAtom_Type_param%Spher_SkR = Spher_SkR_param
	mAtom_Type_param%mSki = mSki_param
	mAtom_Type_param%SkR_std = SkR_std_param
	mAtom_Type_param%MU = MU_param
	mAtom_Type_param%mSkR = mSkR_param
	mAtom_Type_param%LX = LX_param
	mAtom_Type_param%ChemSymb = ChemSymb_param
	mAtom_Type_param%cbas = cbas_param
	mAtom_Type_param%Ind = Ind_param
	mAtom_Type_param%cbas_std = cbas_std_param
	mAtom_Type_param%ThType = ThType_param
	mAtom_Type_param%AtmInfo = AtmInfo_param
	mAtom_Type_param%Ueq = Ueq_param
	mAtom_Type_param%MX = MX_param
end subroutine mAtom_Type_ctor
