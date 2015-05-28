subroutine get_matom_mmphas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(out) :: output_value
	output_value = obj_var%mmphas
end subroutine get_matom_mmphas

subroutine set_matom_mmphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mmphas = new_value
end subroutine set_matom_mmphas

function get_matom_LOcc(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_LOcc
	get_matom_LOcc = obj_var%LOcc
end function get_matom_LOcc

subroutine set_matom_LOcc(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LOcc = new_value
end subroutine set_matom_LOcc

subroutine get_matom_SkI_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%SkI_std
end subroutine get_matom_SkI_std

subroutine set_matom_SkI_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkI_std = new_value
end subroutine set_matom_SkI_std

subroutine get_matom_LVarF(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(25), intent(out) :: output_value
	output_value = obj_var%LVarF
end subroutine get_matom_LVarF

subroutine set_matom_LVarF(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(25), intent(in) :: new_value
	obj_var%LVarF = new_value
end subroutine set_matom_LVarF

subroutine get_matom_Utype(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(out) :: output_value
	output_value = obj_var%Utype
end subroutine get_matom_Utype

subroutine set_matom_Utype(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%Utype = new_value
end subroutine set_matom_Utype

function get_matom_MBiso(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_MBiso
	get_matom_MBiso = obj_var%MBiso
end function get_matom_MBiso

subroutine set_matom_MBiso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MBiso = new_value
end subroutine set_matom_MBiso

function get_matom_Occ(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Occ
	get_matom_Occ = obj_var%Occ
end function get_matom_Occ

subroutine set_matom_Occ(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ = new_value
end subroutine set_matom_Occ

subroutine get_matom_lbas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12,12), intent(out) :: output_value
	output_value = obj_var%lbas
end subroutine get_matom_lbas

subroutine set_matom_lbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12,12), intent(in) :: new_value
	obj_var%lbas = new_value
end subroutine set_matom_lbas

function get_matom_Charge(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Charge
	get_matom_Charge = obj_var%Charge
end function get_matom_Charge

subroutine set_matom_Charge(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Charge = new_value
end subroutine set_matom_Charge

subroutine get_matom_SfacSymb(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(out) :: output_value
	output_value = obj_var%SfacSymb
end subroutine get_matom_SfacSymb

subroutine set_matom_SfacSymb(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=4), intent(in) :: new_value
	obj_var%SfacSymb = new_value
end subroutine set_matom_SfacSymb

subroutine get_matom_imat(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(12), intent(out) :: output_value
	output_value = obj_var%imat
end subroutine get_matom_imat

subroutine set_matom_imat(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(12), intent(in) :: new_value
	obj_var%imat = new_value
end subroutine set_matom_imat

subroutine get_matom_lmphas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12), intent(out) :: output_value
	output_value = obj_var%lmphas
end subroutine get_matom_lmphas

subroutine set_matom_lmphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,dimension(12), intent(in) :: new_value
	obj_var%lmphas = new_value
end subroutine set_matom_lmphas

subroutine get_matom_mphas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(out) :: output_value
	output_value = obj_var%mphas
end subroutine get_matom_mphas

subroutine set_matom_mphas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mphas = new_value
end subroutine set_matom_mphas

subroutine get_matom_Spher_SkI(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%Spher_SkI
end subroutine get_matom_Spher_SkI

subroutine set_matom_Spher_SkI(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkI = new_value
end subroutine set_matom_Spher_SkI

subroutine get_matom_SkR(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%SkR
end subroutine get_matom_SkR

subroutine set_matom_SkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkR = new_value
end subroutine set_matom_SkR

subroutine get_matom_lskr(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(out) :: output_value
	output_value = obj_var%lskr
end subroutine get_matom_lskr

subroutine set_matom_lskr(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(in) :: new_value
	obj_var%lskr = new_value
end subroutine set_matom_lskr

subroutine get_matom_Lab(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=10), intent(out) :: output_value
	output_value = obj_var%Lab
end subroutine get_matom_Lab

subroutine set_matom_Lab(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=10), intent(in) :: new_value
	obj_var%Lab = new_value
end subroutine set_matom_Lab

function get_matom_Moment(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Moment
	get_matom_Moment = obj_var%Moment
end function get_matom_Moment

subroutine set_matom_Moment(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Moment = new_value
end subroutine set_matom_Moment

subroutine get_matom_LU(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(6), intent(out) :: output_value
	output_value = obj_var%LU
end subroutine get_matom_LU

subroutine set_matom_LU(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(6), intent(in) :: new_value
	obj_var%LU = new_value
end subroutine set_matom_LU

function get_matom_MOcc(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_MOcc
	get_matom_MOcc = obj_var%MOcc
end function get_matom_MOcc

subroutine set_matom_MOcc(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%MOcc = new_value
end subroutine set_matom_MOcc

function get_matom_Active(obj_var)
	type (mAtom_Type) :: obj_var
	logical :: get_matom_Active
	get_matom_Active = obj_var%Active
end function get_matom_Active

subroutine set_matom_Active(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Active = new_value
end subroutine set_matom_Active

subroutine get_matom_SkI(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%SkI
end subroutine get_matom_SkI

subroutine set_matom_SkI(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkI = new_value
end subroutine set_matom_SkI

function get_matom_Mult(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_Mult
	get_matom_Mult = obj_var%Mult
end function get_matom_Mult

subroutine set_matom_Mult(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine set_matom_Mult

subroutine get_matom_Spher_SkR_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%Spher_SkR_std
end subroutine get_matom_Spher_SkR_std

subroutine set_matom_Spher_SkR_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkR_std = new_value
end subroutine set_matom_Spher_SkR_std

subroutine get_matom_X_Std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%X_Std
end subroutine get_matom_X_Std

subroutine set_matom_X_Std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X_Std = new_value
end subroutine set_matom_X_Std

subroutine get_matom_U_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%U_std
end subroutine get_matom_U_std

subroutine set_matom_U_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U_std = new_value
end subroutine set_matom_U_std

subroutine get_matom_lski(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(out) :: output_value
	output_value = obj_var%lski
end subroutine get_matom_lski

subroutine set_matom_lski(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3,12), intent(in) :: new_value
	obj_var%lski = new_value
end subroutine set_matom_lski

function get_matom_NVar(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_NVar
	get_matom_NVar = obj_var%NVar
end function get_matom_NVar

subroutine set_matom_NVar(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NVar = new_value
end subroutine set_matom_NVar

subroutine get_matom_wyck(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=1), intent(out) :: output_value
	output_value = obj_var%wyck
end subroutine get_matom_wyck

subroutine set_matom_wyck(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%wyck = new_value
end subroutine set_matom_wyck

function get_matom_Biso_std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Biso_std
	get_matom_Biso_std = obj_var%Biso_std
end function get_matom_Biso_std

subroutine set_matom_Biso_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso_std = new_value
end subroutine set_matom_Biso_std

function get_matom_LBiso(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_LBiso
	get_matom_LBiso = obj_var%LBiso
end function get_matom_LBiso

subroutine set_matom_LBiso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%LBiso = new_value
end subroutine set_matom_LBiso

subroutine get_matom_mphas_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(out) :: output_value
	output_value = obj_var%mphas_std
end subroutine get_matom_mphas_std

subroutine set_matom_mphas_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12), intent(in) :: new_value
	obj_var%mphas_std = new_value
end subroutine set_matom_mphas_std

subroutine get_matom_mVarF(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(out) :: output_value
	output_value = obj_var%mVarF
end subroutine get_matom_mVarF

subroutine set_matom_mVarF(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(in) :: new_value
	obj_var%mVarF = new_value
end subroutine set_matom_mVarF

function get_matom_Biso(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Biso
	get_matom_Biso = obj_var%Biso
end function get_matom_Biso

subroutine set_matom_Biso(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Biso = new_value
end subroutine set_matom_Biso

subroutine get_matom_VarF(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(out) :: output_value
	output_value = obj_var%VarF
end subroutine get_matom_VarF

subroutine set_matom_VarF(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(25), intent(in) :: new_value
	obj_var%VarF = new_value
end subroutine set_matom_VarF

subroutine get_matom_U(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%U
end subroutine get_matom_U

subroutine set_matom_U(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%U = new_value
end subroutine set_matom_U

function get_matom_Occ_Std(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Occ_Std
	get_matom_Occ_Std = obj_var%Occ_Std
end function get_matom_Occ_Std

subroutine set_matom_Occ_Std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Occ_Std = new_value
end subroutine set_matom_Occ_Std

subroutine get_matom_X(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%X
end subroutine get_matom_X

subroutine set_matom_X(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%X = new_value
end subroutine set_matom_X

function get_matom_Z(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_Z
	get_matom_Z = obj_var%Z
end function get_matom_Z

subroutine set_matom_Z(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Z = new_value
end subroutine set_matom_Z

function get_matom_nvk(obj_var)
	type (mAtom_Type) :: obj_var
	integer :: get_matom_nvk
	get_matom_nvk = obj_var%nvk
end function get_matom_nvk

subroutine set_matom_nvk(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nvk = new_value
end subroutine set_matom_nvk

subroutine get_matom_mbas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(out) :: output_value
	output_value = obj_var%mbas
end subroutine get_matom_mbas

subroutine set_matom_mbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%mbas = new_value
end subroutine set_matom_mbas

subroutine get_matom_Spher_SkI_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%Spher_SkI_std
end subroutine get_matom_Spher_SkI_std

subroutine set_matom_Spher_SkI_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkI_std = new_value
end subroutine set_matom_Spher_SkI_std

subroutine get_matom_Spher_SkR(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%Spher_SkR
end subroutine get_matom_Spher_SkR

subroutine set_matom_Spher_SkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%Spher_SkR = new_value
end subroutine set_matom_Spher_SkR

subroutine get_matom_mSki(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%mSki
end subroutine get_matom_mSki

subroutine set_matom_mSki(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%mSki = new_value
end subroutine set_matom_mSki

subroutine get_matom_SkR_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%SkR_std
end subroutine get_matom_SkR_std

subroutine set_matom_SkR_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%SkR_std = new_value
end subroutine set_matom_SkR_std

subroutine get_matom_MU(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(out) :: output_value
	output_value = obj_var%MU
end subroutine get_matom_MU

subroutine set_matom_MU(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(6), intent(in) :: new_value
	obj_var%MU = new_value
end subroutine set_matom_MU

subroutine get_matom_mSkR(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%mSkR
end subroutine get_matom_mSkR

subroutine set_matom_mSkR(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%mSkR = new_value
end subroutine set_matom_mSkR

subroutine get_matom_LX(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3), intent(out) :: output_value
	output_value = obj_var%LX
end subroutine get_matom_LX

subroutine set_matom_LX(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer,      dimension(3), intent(in) :: new_value
	obj_var%LX = new_value
end subroutine set_matom_LX

subroutine get_matom_ChemSymb(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=2), intent(out) :: output_value
	output_value = obj_var%ChemSymb
end subroutine get_matom_ChemSymb

subroutine set_matom_ChemSymb(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=2), intent(in) :: new_value
	obj_var%ChemSymb = new_value
end subroutine set_matom_ChemSymb

subroutine get_matom_cbas(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(out) :: output_value
	output_value = obj_var%cbas
end subroutine get_matom_cbas

subroutine set_matom_cbas(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%cbas = new_value
end subroutine set_matom_cbas

subroutine get_matom_Ind(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	integer, dimension(5), intent(out) :: output_value
	output_value = obj_var%Ind
end subroutine get_matom_Ind

subroutine set_matom_Ind(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	integer, dimension(5), intent(in) :: new_value
	obj_var%Ind = new_value
end subroutine set_matom_Ind

subroutine get_matom_cbas_std(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(out) :: output_value
	output_value = obj_var%cbas_std
end subroutine get_matom_cbas_std

subroutine set_matom_cbas_std(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(12,12), intent(in) :: new_value
	obj_var%cbas_std = new_value
end subroutine set_matom_cbas_std

subroutine get_matom_ThType(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=5), intent(out) :: output_value
	output_value = obj_var%ThType
end subroutine get_matom_ThType

subroutine set_matom_ThType(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=5), intent(in) :: new_value
	obj_var%ThType = new_value
end subroutine set_matom_ThType

subroutine get_matom_AtmInfo(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	character(len=40), intent(out) :: output_value
	output_value = obj_var%AtmInfo
end subroutine get_matom_AtmInfo

subroutine set_matom_AtmInfo(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	character(len=40), intent(in) :: new_value
	obj_var%AtmInfo = new_value
end subroutine set_matom_AtmInfo

function get_matom_Ueq(obj_var)
	type (mAtom_Type) :: obj_var
	real(kind=cp) :: get_matom_Ueq
	get_matom_Ueq = obj_var%Ueq
end function get_matom_Ueq

subroutine set_matom_Ueq(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Ueq = new_value
end subroutine set_matom_Ueq

subroutine get_matom_MX(obj_var, output_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(out) :: output_value
	output_value = obj_var%MX
end subroutine get_matom_MX

subroutine set_matom_MX(obj_var, new_value)
	type (mAtom_Type) :: obj_var
	real(kind=cp),dimension(3), intent(in) :: new_value
	obj_var%MX = new_value
end subroutine set_matom_MX

subroutine mAtom_Type_ctor(mAtom_Type_param, mmphas_param, LOcc_param, SkI_std_param, LVarF_param, Utype_param, MBiso_param, Occ_param, lbas_param, Charge_param, SfacSymb_param, imat_param, lmphas_param, mphas_param, Spher_SkI_param, SkR_param, lskr_param, Lab_param, Moment_param, LU_param, MOcc_param, Active_param, SkI_param, Mult_param, Spher_SkR_std_param, X_Std_param, U_std_param, lski_param, NVar_param, wyck_param, Biso_std_param, LBiso_param, mphas_std_param, mVarF_param, Biso_param, VarF_param, U_param, Occ_Std_param, X_param, Z_param, nvk_param, mbas_param, Spher_SkI_std_param, Spher_SkR_param, mSki_param, SkR_std_param, MU_param, mSkR_param, LX_param, ChemSymb_param, cbas_param, Ind_param, cbas_std_param, ThType_param, AtmInfo_param, Ueq_param, MX_param)
	type (mAtom_Type) :: mAtom_Type_param
	real(kind=cp),dimension(12), intent(in) :: mmphas_param
	integer, intent(in) :: LOcc_param
	real(kind=cp),dimension(3,12), intent(in) :: SkI_std_param
	integer,      dimension(25), intent(in) :: LVarF_param
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
	real(kind=cp),dimension(25), intent(in) :: mVarF_param
	real(kind=cp), intent(in) :: Biso_param
	real(kind=cp),dimension(25), intent(in) :: VarF_param
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
	mAtom_Type_param%LVarF = LVarF_param
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
	mAtom_Type_param%mVarF = mVarF_param
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
