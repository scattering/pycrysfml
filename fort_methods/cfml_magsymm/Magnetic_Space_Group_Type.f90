function getm_constr(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: getm_constr
	getm_constr = obj_var%m_constr
end function getm_constr

subroutine setm_constr(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%m_constr = new_value
end subroutine setm_constr

function getCentre_coord(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3) :: getCentre_coord
	getCentre_coord = obj_var%Centre_coord
end function getCentre_coord

subroutine setCentre_coord(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Centre_coord = new_value
end subroutine setCentre_coord

function getaLatt_trans(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:) :: getaLatt_trans
	getaLatt_trans = obj_var%aLatt_trans
end function getaLatt_trans

subroutine setaLatt_trans(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: new_value
	obj_var%aLatt_trans = new_value
end subroutine setaLatt_trans

function getMSymopSymb(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable :: getMSymopSymb
	getMSymopSymb = obj_var%MSymopSymb
end function getMSymopSymb

subroutine setMSymopSymb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%MSymopSymb = new_value
end subroutine setMSymopSymb

function getn_wyck(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getn_wyck
	getn_wyck = obj_var%n_wyck
end function getn_wyck

subroutine setn_wyck(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_wyck = new_value
end subroutine setn_wyck

function getBNS_symbol(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34) :: getBNS_symbol
	getBNS_symbol = obj_var%BNS_symbol
end function getBNS_symbol

subroutine setBNS_symbol(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%BNS_symbol = new_value
end subroutine setBNS_symbol

function getBNS_number(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15) :: getBNS_number
	getBNS_number = obj_var%BNS_number
end function getBNS_number

subroutine setBNS_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%BNS_number = new_value
end subroutine setBNS_number

function getLatt_trans(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:) :: getLatt_trans
	getLatt_trans = obj_var%Latt_trans
end function getLatt_trans

subroutine setLatt_trans(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: new_value
	obj_var%Latt_trans = new_value
end subroutine setLatt_trans

function getSh_number(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getSh_number
	getSh_number = obj_var%Sh_number
end function getSh_number

subroutine setSh_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Sh_number = new_value
end subroutine setSh_number

function getWyck_Symb(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable :: getWyck_Symb
	getWyck_Symb = obj_var%Wyck_Symb
end function getWyck_Symb

subroutine setWyck_Symb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%Wyck_Symb = new_value
end subroutine setWyck_Symb

function getCrystalSys(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=12) :: getCrystalSys
	getCrystalSys = obj_var%CrystalSys
end function getCrystalSys

subroutine setCrystalSys(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%CrystalSys = new_value
end subroutine setCrystalSys

function getm_cell(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: getm_cell
	getm_cell = obj_var%m_cell
end function getm_cell

subroutine setm_cell(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%m_cell = new_value
end subroutine setm_cell

function getirrep_id(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable :: getirrep_id
	getirrep_id = obj_var%irrep_id
end function getirrep_id

subroutine setirrep_id(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_id = new_value
end subroutine setirrep_id

function getMSymOp(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:),allocatable :: getMSymOp
	getMSymOp = obj_var%MSymOp
end function getMSymOp

subroutine setMSymOp(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:),allocatable, intent(in) :: new_value
	obj_var%MSymOp = new_value
end subroutine setMSymOp

function getNum_Lat(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: getNum_Lat
	getNum_Lat = obj_var%Num_Lat
end function getNum_Lat

subroutine setNum_Lat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Lat = new_value
end subroutine setNum_Lat

function getNum_aLat(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: getNum_aLat
	getNum_aLat = obj_var%Num_aLat
end function getNum_aLat

subroutine setNum_aLat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_aLat = new_value
end subroutine setNum_aLat

function getsmall_irrep_dim(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable :: getsmall_irrep_dim
	getsmall_irrep_dim = obj_var%small_irrep_dim
end function getsmall_irrep_dim

subroutine setsmall_irrep_dim(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%small_irrep_dim = new_value
end subroutine setsmall_irrep_dim

function getCentred(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: getCentred
	getCentred = obj_var%Centred
end function getCentred

subroutine setCentred(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Centred = new_value
end subroutine setCentred

function getSPG_lat(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 1) :: getSPG_lat
	getSPG_lat = obj_var%SPG_lat
end function getSPG_lat

subroutine setSPG_lat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 1), intent(in) :: new_value
	obj_var%SPG_lat = new_value
end subroutine setSPG_lat

function getSymopSymb(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable :: getSymopSymb
	getSymopSymb = obj_var%SymopSymb
end function getSymopSymb

subroutine setSymopSymb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine setSymopSymb

function getParent_num(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getParent_num
	getParent_num = obj_var%Parent_num
end function getParent_num

subroutine setParent_num(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Parent_num = new_value
end subroutine setParent_num

function getkv_label(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable :: getkv_label
	getkv_label = obj_var%kv_label
end function getkv_label

subroutine setkv_label(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%kv_label = new_value
end subroutine setkv_label

function getMagType(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getMagType
	getMagType = obj_var%MagType
end function getMagType

subroutine setMagType(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%MagType = new_value
end subroutine setMagType

function getirrep_action(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable :: getirrep_action
	getirrep_action = obj_var%irrep_action
end function getirrep_action

subroutine setirrep_action(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_action = new_value
end subroutine setirrep_action

function getn_kv(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getn_kv
	getn_kv = obj_var%n_kv
end function getn_kv

subroutine setn_kv(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_kv = new_value
end subroutine setn_kv

function getOG_number(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15) :: getOG_number
	getOG_number = obj_var%OG_number
end function getOG_number

subroutine setOG_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%OG_number = new_value
end subroutine setOG_number

function getn_irreps(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getn_irreps
	getn_irreps = obj_var%n_irreps
end function getn_irreps

subroutine setn_irreps(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_irreps = new_value
end subroutine setn_irreps

function getSPG_latsy(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 2) :: getSPG_latsy
	getSPG_latsy = obj_var%SPG_latsy
end function getSPG_latsy

subroutine setSPG_latsy(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 2), intent(in) :: new_value
	obj_var%SPG_latsy = new_value
end subroutine setSPG_latsy

function getNum_gen(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: getNum_gen
	getNum_gen = obj_var%Num_gen
end function getNum_gen

subroutine setNum_gen(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_gen = new_value
end subroutine setNum_gen

function getirrep_direction(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable :: getirrep_direction
	getirrep_direction = obj_var%irrep_direction
end function getirrep_direction

subroutine setirrep_direction(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_direction = new_value
end subroutine setirrep_direction

function getstandard_setting(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: getstandard_setting
	getstandard_setting = obj_var%standard_setting
end function getstandard_setting

subroutine setstandard_setting(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%standard_setting = new_value
end subroutine setstandard_setting

function getSymOp(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),allocatable :: getSymOp
	getSymOp = obj_var%SymOp
end function getSymOp

subroutine setSymOp(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),allocatable, intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine setSymOp

function getirrep_dim(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable :: getirrep_dim
	getirrep_dim = obj_var%irrep_dim
end function getirrep_dim

subroutine setirrep_dim(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_dim = new_value
end subroutine setirrep_dim

function getParent_spg(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20) :: getParent_spg
	getParent_spg = obj_var%Parent_spg
end function getParent_spg

subroutine setParent_spg(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20), intent(in) :: new_value
	obj_var%Parent_spg = new_value
end subroutine setParent_spg

function getOG_symbol(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34) :: getOG_symbol
	getOG_symbol = obj_var%OG_symbol
end function getOG_symbol

subroutine setOG_symbol(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%OG_symbol = new_value
end subroutine setOG_symbol

function getCentre(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=80) :: getCentre
	getCentre = obj_var%Centre
end function getCentre

subroutine setCentre(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=80), intent(in) :: new_value
	obj_var%Centre = new_value
end subroutine setCentre

function getmcif(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: getmcif
	getmcif = obj_var%mcif
end function getmcif

subroutine setmcif(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%mcif = new_value
end subroutine setmcif

function getNumOps(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: getNumOps
	getNumOps = obj_var%NumOps
end function getNumOps

subroutine setNumOps(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumOps = new_value
end subroutine setNumOps

function getirrep_modes_number(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable :: getirrep_modes_number
	getirrep_modes_number = obj_var%irrep_modes_number
end function getirrep_modes_number

subroutine setirrep_modes_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_modes_number = new_value
end subroutine setirrep_modes_number

function getkv(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp),     dimension(:,:),allocatable :: getkv
	getkv = obj_var%kv
end function getkv

subroutine setkv(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp),     dimension(:,:),allocatable, intent(in) :: new_value
	obj_var%kv = new_value
end subroutine setkv

function gettrn_to_standard(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40) :: gettrn_to_standard
	gettrn_to_standard = obj_var%trn_to_standard
end function gettrn_to_standard

subroutine settrn_to_standard(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(in) :: new_value
	obj_var%trn_to_standard = new_value
end subroutine settrn_to_standard

function gettrn_from_parent(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40) :: gettrn_from_parent
	gettrn_from_parent = obj_var%trn_from_parent
end function gettrn_from_parent

subroutine settrn_from_parent(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(in) :: new_value
	obj_var%trn_from_parent = new_value
end subroutine settrn_from_parent

function getMultip(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: getMultip
	getMultip = obj_var%Multip
end function getMultip

subroutine setMultip(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine setMultip

subroutine Magnetic_Space_Group_Type_ctor(Magnetic_Space_Group_Type_param, m_constr_param, Centre_coord_param, aLatt_trans_param, MSymopSymb_param, n_wyck_param, BNS_symbol_param, BNS_number_param, Latt_trans_param, Sh_number_param, Wyck_Symb_param, CrystalSys_param, m_cell_param, irrep_id_param, MSymOp_param, Num_Lat_param, Num_aLat_param, small_irrep_dim_param, Centred_param, SPG_lat_param, SymopSymb_param, Parent_num_param, kv_label_param, MagType_param, irrep_action_param, n_kv_param, OG_number_param, n_irreps_param, SPG_latsy_param, Num_gen_param, irrep_direction_param, standard_setting_param, SymOp_param, irrep_dim_param, Parent_spg_param, OG_symbol_param, Centre_param, mcif_param, NumOps_param, irrep_modes_number_param, kv_param, trn_to_standard_param, trn_from_parent_param, Multip_param)
	type (Magnetic_Space_Group_Type) :: Magnetic_Space_Group_Type_param
	logical, intent(in) :: m_constr_param
	real(kind=cp), dimension(3), intent(in) :: Centre_coord_param
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: aLatt_trans_param
	character(len=40),   dimension(:),allocatable, intent(in) :: MSymopSymb_param
	Integer, intent(in) :: n_wyck_param
	Character(len=34), intent(in) :: BNS_symbol_param
	character(len=15), intent(in) :: BNS_number_param
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: Latt_trans_param
	Integer, intent(in) :: Sh_number_param
	character(len=40),   dimension(:),allocatable, intent(in) :: Wyck_Symb_param
	character(len=12), intent(in) :: CrystalSys_param
	logical, intent(in) :: m_cell_param
	Character(len=15),   dimension(:),allocatable, intent(in) :: irrep_id_param
	type(MSym_Oper_Type),dimension(:),allocatable, intent(in) :: MSymOp_param
	integer, intent(in) :: Num_Lat_param
	integer, intent(in) :: Num_aLat_param
	Integer,             dimension(:),allocatable, intent(in) :: small_irrep_dim_param
	integer, intent(in) :: Centred_param
	character(len= 1), intent(in) :: SPG_lat_param
	character(len=40),   dimension(:),allocatable, intent(in) :: SymopSymb_param
	Integer, intent(in) :: Parent_num_param
	Character(len=15),   dimension(:),allocatable, intent(in) :: kv_label_param
	Integer, intent(in) :: MagType_param
	Character(len=20),   dimension(:),allocatable, intent(in) :: irrep_action_param
	Integer, intent(in) :: n_kv_param
	character(len=15), intent(in) :: OG_number_param
	Integer, intent(in) :: n_irreps_param
	character(len= 2), intent(in) :: SPG_latsy_param
	integer, intent(in) :: Num_gen_param
	Character(len=20),   dimension(:),allocatable, intent(in) :: irrep_direction_param
	logical, intent(in) :: standard_setting_param
	type(Sym_Oper_Type), dimension(:),allocatable, intent(in) :: SymOp_param
	Integer,             dimension(:),allocatable, intent(in) :: irrep_dim_param
	Character(len=20), intent(in) :: Parent_spg_param
	Character(len=34), intent(in) :: OG_symbol_param
	character(len=80), intent(in) :: Centre_param
	logical, intent(in) :: mcif_param
	integer, intent(in) :: NumOps_param
	Integer,             dimension(:),allocatable, intent(in) :: irrep_modes_number_param
	real(kind=cp),     dimension(:,:),allocatable, intent(in) :: kv_param
	Character(len=40), intent(in) :: trn_to_standard_param
	Character(len=40), intent(in) :: trn_from_parent_param
	Integer, intent(in) :: Multip_param
	Magnetic_Space_Group_Type_param%m_constr = m_constr_param
	Magnetic_Space_Group_Type_param%Centre_coord = Centre_coord_param
	Magnetic_Space_Group_Type_param%aLatt_trans = aLatt_trans_param
	Magnetic_Space_Group_Type_param%MSymopSymb = MSymopSymb_param
	Magnetic_Space_Group_Type_param%n_wyck = n_wyck_param
	Magnetic_Space_Group_Type_param%BNS_symbol = BNS_symbol_param
	Magnetic_Space_Group_Type_param%BNS_number = BNS_number_param
	Magnetic_Space_Group_Type_param%Latt_trans = Latt_trans_param
	Magnetic_Space_Group_Type_param%Sh_number = Sh_number_param
	Magnetic_Space_Group_Type_param%Wyck_Symb = Wyck_Symb_param
	Magnetic_Space_Group_Type_param%CrystalSys = CrystalSys_param
	Magnetic_Space_Group_Type_param%m_cell = m_cell_param
	Magnetic_Space_Group_Type_param%irrep_id = irrep_id_param
	Magnetic_Space_Group_Type_param%MSymOp = MSymOp_param
	Magnetic_Space_Group_Type_param%Num_Lat = Num_Lat_param
	Magnetic_Space_Group_Type_param%Num_aLat = Num_aLat_param
	Magnetic_Space_Group_Type_param%small_irrep_dim = small_irrep_dim_param
	Magnetic_Space_Group_Type_param%Centred = Centred_param
	Magnetic_Space_Group_Type_param%SPG_lat = SPG_lat_param
	Magnetic_Space_Group_Type_param%SymopSymb = SymopSymb_param
	Magnetic_Space_Group_Type_param%Parent_num = Parent_num_param
	Magnetic_Space_Group_Type_param%kv_label = kv_label_param
	Magnetic_Space_Group_Type_param%MagType = MagType_param
	Magnetic_Space_Group_Type_param%irrep_action = irrep_action_param
	Magnetic_Space_Group_Type_param%n_kv = n_kv_param
	Magnetic_Space_Group_Type_param%OG_number = OG_number_param
	Magnetic_Space_Group_Type_param%n_irreps = n_irreps_param
	Magnetic_Space_Group_Type_param%SPG_latsy = SPG_latsy_param
	Magnetic_Space_Group_Type_param%Num_gen = Num_gen_param
	Magnetic_Space_Group_Type_param%irrep_direction = irrep_direction_param
	Magnetic_Space_Group_Type_param%standard_setting = standard_setting_param
	Magnetic_Space_Group_Type_param%SymOp = SymOp_param
	Magnetic_Space_Group_Type_param%irrep_dim = irrep_dim_param
	Magnetic_Space_Group_Type_param%Parent_spg = Parent_spg_param
	Magnetic_Space_Group_Type_param%OG_symbol = OG_symbol_param
	Magnetic_Space_Group_Type_param%Centre = Centre_param
	Magnetic_Space_Group_Type_param%mcif = mcif_param
	Magnetic_Space_Group_Type_param%NumOps = NumOps_param
	Magnetic_Space_Group_Type_param%irrep_modes_number = irrep_modes_number_param
	Magnetic_Space_Group_Type_param%kv = kv_param
	Magnetic_Space_Group_Type_param%trn_to_standard = trn_to_standard_param
	Magnetic_Space_Group_Type_param%trn_from_parent = trn_from_parent_param
	Magnetic_Space_Group_Type_param%Multip = Multip_param
end subroutine Magnetic_Space_Group_Type_ctor
