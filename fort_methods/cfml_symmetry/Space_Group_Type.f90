function getCentred(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getCentred
	getCentred = obj_var%Centred
end function getCentred

subroutine setCentred(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Centred = new_value
end subroutine setCentred

function getCentre_coord(obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3) :: getCentre_coord
	getCentre_coord = obj_var%Centre_coord
end function getCentre_coord

subroutine setCentre_coord(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Centre_coord = new_value
end subroutine setCentre_coord

function getLatt_trans(obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:) :: getLatt_trans
	getLatt_trans = obj_var%Latt_trans
end function getLatt_trans

subroutine setLatt_trans(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: new_value
	obj_var%Latt_trans = new_value
end subroutine setLatt_trans

function getHexa(obj_var)
	type (Space_Group_Type) :: obj_var
	logical :: getHexa
	getHexa = obj_var%Hexa
end function getHexa

subroutine setHexa(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Hexa = new_value
end subroutine setHexa

function getNumSpg(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumSpg
	getNumSpg = obj_var%NumSpg
end function getNumSpg

subroutine setNumSpg(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumSpg = new_value
end subroutine setNumSpg

function getSymopSymb(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=50),   allocatable,dimension(:) :: getSymopSymb
	getSymopSymb = obj_var%SymopSymb
end function getSymopSymb

subroutine setSymopSymb(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=50),   allocatable,dimension(:), intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine setSymopSymb

function getCrystalSys(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=12) :: getCrystalSys
	getCrystalSys = obj_var%CrystalSys
end function getCrystalSys

subroutine setCrystalSys(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%CrystalSys = new_value
end subroutine setCrystalSys

function getNumLat(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumLat
	getNumLat = obj_var%NumLat
end function getNumLat

subroutine setNumLat(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumLat = new_value
end subroutine setNumLat

function getPG(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getPG
	getPG = obj_var%PG
end function getPG

subroutine setPG(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%PG = new_value
end subroutine setPG

function getWyckoff(obj_var)
	type (Space_Group_Type) :: obj_var
	type(Wyckoff_Type) :: getWyckoff
	getWyckoff = obj_var%Wyckoff
end function getWyckoff

subroutine setWyckoff(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	type(Wyckoff_Type), intent(in) :: new_value
	obj_var%Wyckoff = new_value
end subroutine setWyckoff

function getHall(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=16) :: getHall
	getHall = obj_var%Hall
end function getHall

subroutine setHall(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=16), intent(in) :: new_value
	obj_var%Hall = new_value
end subroutine setHall

function getInfo(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getInfo
	getInfo = obj_var%Info
end function getInfo

subroutine setInfo(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Info = new_value
end subroutine setInfo

function getSPG_lat(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 1) :: getSPG_lat
	getSPG_lat = obj_var%SPG_lat
end function getSPG_lat

subroutine setSPG_lat(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 1), intent(in) :: new_value
	obj_var%SPG_lat = new_value
end subroutine setSPG_lat

function getLaue(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getLaue
	getLaue = obj_var%Laue
end function getLaue

subroutine setLaue(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Laue = new_value
end subroutine setLaue

function getSPG_latsy(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 2) :: getSPG_latsy
	getSPG_latsy = obj_var%SPG_latsy
end function getSPG_latsy

subroutine setSPG_latsy(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 2), intent(in) :: new_value
	obj_var%SPG_latsy = new_value
end subroutine setSPG_latsy

function getNum_gen(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNum_gen
	getNum_gen = obj_var%Num_gen
end function getNum_gen

subroutine setNum_gen(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_gen = new_value
end subroutine setNum_gen

function getBravais(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=51) :: getBravais
	getBravais = obj_var%Bravais
end function getBravais

subroutine setBravais(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=51), intent(in) :: new_value
	obj_var%Bravais = new_value
end subroutine setBravais

function getSG_setting(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=90) :: getSG_setting
	getSG_setting = obj_var%SG_setting
end function getSG_setting

subroutine setSG_setting(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%SG_setting = new_value
end subroutine setSG_setting

function getgHall(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=90) :: getgHall
	getgHall = obj_var%gHall
end function getgHall

subroutine setgHall(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%gHall = new_value
end subroutine setgHall

function getSPG_Symb(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=20) :: getSPG_Symb
	getSPG_Symb = obj_var%SPG_Symb
end function getSPG_Symb

subroutine setSPG_Symb(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%SPG_Symb = new_value
end subroutine setSPG_Symb

function getCentre(obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=80) :: getCentre
	getCentre = obj_var%Centre
end function getCentre

subroutine setCentre(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=80), intent(in) :: new_value
	obj_var%Centre = new_value
end subroutine setCentre

function getSymOp(obj_var)
	type (Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), allocatable,dimension(:) :: getSymOp
	getSymOp = obj_var%SymOp
end function getSymOp

subroutine setSymOp(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), allocatable,dimension(:), intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine setSymOp

function getNumOps(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumOps
	getNumOps = obj_var%NumOps
end function getNumOps

subroutine setNumOps(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumOps = new_value
end subroutine setNumOps

function getR_Asym_Unit(obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp),dimension(3,2) :: getR_Asym_Unit
	getR_Asym_Unit = obj_var%R_Asym_Unit
end function getR_Asym_Unit

subroutine setR_Asym_Unit(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp),dimension(3,2), intent(in) :: new_value
	obj_var%R_Asym_Unit = new_value
end subroutine setR_Asym_Unit

function getMultip(obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getMultip
	getMultip = obj_var%Multip
end function getMultip

subroutine setMultip(obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine setMultip

subroutine Space_Group_Type_ctor(Space_Group_Type_param, Centred_param, Centre_coord_param, Latt_trans_param, Hexa_param, NumSpg_param, SymopSymb_param, CrystalSys_param, NumLat_param, PG_param, Wyckoff_param, Hall_param, Info_param, SPG_lat_param, Laue_param, SPG_latsy_param, Num_gen_param, Bravais_param, SG_setting_param, gHall_param, SPG_Symb_param, Centre_param, SymOp_param, NumOps_param, R_Asym_Unit_param, Multip_param)
	type (Space_Group_Type) :: Space_Group_Type_param
	integer, intent(in) :: Centred_param
	real(kind=cp), dimension(3), intent(in) :: Centre_coord_param
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: Latt_trans_param
	logical, intent(in) :: Hexa_param
	integer, intent(in) :: NumSpg_param
	character(len=50),   allocatable,dimension(:), intent(in) :: SymopSymb_param
	character(len=12), intent(in) :: CrystalSys_param
	integer, intent(in) :: NumLat_param
	character(len= 5), intent(in) :: PG_param
	type(Wyckoff_Type), intent(in) :: Wyckoff_param
	character(len=16), intent(in) :: Hall_param
	character(len= 5), intent(in) :: Info_param
	character(len= 1), intent(in) :: SPG_lat_param
	character(len= 5), intent(in) :: Laue_param
	character(len= 2), intent(in) :: SPG_latsy_param
	integer, intent(in) :: Num_gen_param
	character(len=51), intent(in) :: Bravais_param
	character(len=90), intent(in) :: SG_setting_param
	character(len=90), intent(in) :: gHall_param
	character(len=20), intent(in) :: SPG_Symb_param
	character(len=80), intent(in) :: Centre_param
	type(Sym_Oper_Type), allocatable,dimension(:), intent(in) :: SymOp_param
	integer, intent(in) :: NumOps_param
	real(kind=cp),dimension(3,2), intent(in) :: R_Asym_Unit_param
	integer, intent(in) :: Multip_param
	Space_Group_Type_param%Centred = Centred_param
	Space_Group_Type_param%Centre_coord = Centre_coord_param
	Space_Group_Type_param%Latt_trans = Latt_trans_param
	Space_Group_Type_param%Hexa = Hexa_param
	Space_Group_Type_param%NumSpg = NumSpg_param
	Space_Group_Type_param%SymopSymb = SymopSymb_param
	Space_Group_Type_param%CrystalSys = CrystalSys_param
	Space_Group_Type_param%NumLat = NumLat_param
	Space_Group_Type_param%PG = PG_param
	Space_Group_Type_param%Wyckoff = Wyckoff_param
	Space_Group_Type_param%Hall = Hall_param
	Space_Group_Type_param%Info = Info_param
	Space_Group_Type_param%SPG_lat = SPG_lat_param
	Space_Group_Type_param%Laue = Laue_param
	Space_Group_Type_param%SPG_latsy = SPG_latsy_param
	Space_Group_Type_param%Num_gen = Num_gen_param
	Space_Group_Type_param%Bravais = Bravais_param
	Space_Group_Type_param%SG_setting = SG_setting_param
	Space_Group_Type_param%gHall = gHall_param
	Space_Group_Type_param%SPG_Symb = SPG_Symb_param
	Space_Group_Type_param%Centre = Centre_param
	Space_Group_Type_param%SymOp = SymOp_param
	Space_Group_Type_param%NumOps = NumOps_param
	Space_Group_Type_param%R_Asym_Unit = R_Asym_Unit_param
	Space_Group_Type_param%Multip = Multip_param
end subroutine Space_Group_Type_ctor
