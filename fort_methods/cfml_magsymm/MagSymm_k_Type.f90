function getLatt(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=1) :: getLatt
	getLatt = obj_var%Latt
end function getLatt

subroutine setLatt(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%Latt = new_value
end subroutine setLatt

function getSymopSymb(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:),   allocatable :: getSymopSymb
	getSymopSymb = obj_var%SymopSymb
end function getSymopSymb

subroutine setSymopSymb(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:),   allocatable, intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine setSymopSymb

function getMSymopSymb(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:,:), allocatable :: getMSymopSymb
	getMSymopSymb = obj_var%MSymopSymb
end function getMSymopSymb

subroutine setMSymopSymb(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%MSymopSymb = new_value
end subroutine setMSymopSymb

function getNumops(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getNumops
	getNumops = obj_var%Numops
end function getNumops

subroutine setNumops(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Numops = new_value
end subroutine setNumops

function getBNS_symbol(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34) :: getBNS_symbol
	getBNS_symbol = obj_var%BNS_symbol
end function getBNS_symbol

subroutine setBNS_symbol(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%BNS_symbol = new_value
end subroutine setBNS_symbol

function getBNS_number(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=15) :: getBNS_number
	getBNS_number = obj_var%BNS_number
end function getBNS_number

subroutine setBNS_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%BNS_number = new_value
end subroutine setBNS_number

function getmcentred(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getmcentred
	getmcentred = obj_var%mcentred
end function getmcentred

subroutine setmcentred(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%mcentred = new_value
end subroutine setmcentred

function getMagModel(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=31) :: getMagModel
	getMagModel = obj_var%MagModel
end function getMagModel

subroutine setMagModel(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=31), intent(in) :: new_value
	obj_var%MagModel = new_value
end subroutine setMagModel

function getkvec(obj_var)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp),dimension(3,12) :: getkvec
	getkvec = obj_var%kvec
end function getkvec

subroutine setkvec(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%kvec = new_value
end subroutine setkvec

function getnmsym(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getnmsym
	getnmsym = obj_var%nmsym
end function getnmsym

subroutine setnmsym(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nmsym = new_value
end subroutine setnmsym

function getMSymOp(obj_var)
	type (MagSymm_k_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:,:), allocatable :: getMSymOp
	getMSymOp = obj_var%MSymOp
end function getMSymOp

subroutine setMSymOp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%MSymOp = new_value
end subroutine setMSymOp

function getNum_Lat(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getNum_Lat
	getNum_Lat = obj_var%Num_Lat
end function getNum_Lat

subroutine setNum_Lat(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Lat = new_value
end subroutine setNum_Lat

function getirrep_action(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4) :: getirrep_action
	getirrep_action = obj_var%irrep_action
end function getirrep_action

subroutine setirrep_action(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(in) :: new_value
	obj_var%irrep_action = new_value
end subroutine setirrep_action

function getsmall_irrep_dim(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4) :: getsmall_irrep_dim
	getsmall_irrep_dim = obj_var%small_irrep_dim
end function getsmall_irrep_dim

subroutine setsmall_irrep_dim(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%small_irrep_dim = new_value
end subroutine setsmall_irrep_dim

function getLtr(obj_var)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp), dimension(3,4) :: getLtr
	getLtr = obj_var%Ltr
end function getLtr

subroutine setLtr(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp), dimension(3,4), intent(in) :: new_value
	obj_var%Ltr = new_value
end subroutine setLtr

function getParent_num(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer :: getParent_num
	getParent_num = obj_var%Parent_num
end function getParent_num

subroutine setParent_num(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Parent_num = new_value
end subroutine setParent_num

function getSk_type(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=15) :: getSk_type
	getSk_type = obj_var%Sk_type
end function getSk_type

subroutine setSk_type(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%Sk_type = new_value
end subroutine setSk_type

function getnirreps(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getnirreps
	getnirreps = obj_var%nirreps
end function getnirreps

subroutine setnirreps(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nirreps = new_value
end subroutine setnirreps

function getMagType(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer :: getMagType
	getMagType = obj_var%MagType
end function getMagType

subroutine setMagType(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%MagType = new_value
end subroutine setMagType

function getOG_number(obj_var)
	type (MagSymm_k_Type) :: obj_var
	character(len=15) :: getOG_number
	getOG_number = obj_var%OG_number
end function getOG_number

subroutine setOG_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%OG_number = new_value
end subroutine setOG_number

function getirrep_dim(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4) :: getirrep_dim
	getirrep_dim = obj_var%irrep_dim
end function getirrep_dim

subroutine setirrep_dim(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%irrep_dim = new_value
end subroutine setirrep_dim

function getirrep_modes_number(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4) :: getirrep_modes_number
	getirrep_modes_number = obj_var%irrep_modes_number
end function getirrep_modes_number

subroutine setirrep_modes_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%irrep_modes_number = new_value
end subroutine setirrep_modes_number

function getnkv(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getnkv
	getnkv = obj_var%nkv
end function getnkv

subroutine setnkv(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nkv = new_value
end subroutine setnkv

function getirrep_direction(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4) :: getirrep_direction
	getirrep_direction = obj_var%irrep_direction
end function getirrep_direction

subroutine setirrep_direction(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(in) :: new_value
	obj_var%irrep_direction = new_value
end subroutine setirrep_direction

function getcentred(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getcentred
	getcentred = obj_var%centred
end function getcentred

subroutine setcentred(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%centred = new_value
end subroutine setcentred

function getbasf(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Complex(kind=cp),    dimension(3,12,48,4) :: getbasf
	getbasf = obj_var%basf
end function getbasf

subroutine setbasf(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Complex(kind=cp),    dimension(3,12,48,4), intent(in) :: new_value
	obj_var%basf = new_value
end subroutine setbasf

function getParent_spg(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20) :: getParent_spg
	getParent_spg = obj_var%Parent_spg
end function getParent_spg

subroutine setParent_spg(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20), intent(in) :: new_value
	obj_var%Parent_spg = new_value
end subroutine setParent_spg

function getOG_symbol(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34) :: getOG_symbol
	getOG_symbol = obj_var%OG_symbol
end function getOG_symbol

subroutine setOG_symbol(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%OG_symbol = new_value
end subroutine setOG_symbol

function getSymOp(obj_var)
	type (MagSymm_k_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),   allocatable :: getSymOp
	getSymOp = obj_var%SymOp
end function getSymOp

subroutine setSymOp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),   allocatable, intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine setSymOp

function getirrep_id(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Character(len=15),   dimension(4) :: getirrep_id
	getirrep_id = obj_var%irrep_id
end function getirrep_id

subroutine setirrep_id(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=15),   dimension(4), intent(in) :: new_value
	obj_var%irrep_id = new_value
end subroutine setirrep_id

function getnbas(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(4) :: getnbas
	getnbas = obj_var%nbas
end function getnbas

subroutine setnbas(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(4), intent(in) :: new_value
	obj_var%nbas = new_value
end subroutine setnbas

function geticomp(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(12,4) :: geticomp
	geticomp = obj_var%icomp
end function geticomp

subroutine seticomp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(12,4), intent(in) :: new_value
	obj_var%icomp = new_value
end subroutine seticomp

function getMultip(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: getMultip
	getMultip = obj_var%Multip
end function getMultip

subroutine setMultip(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine setMultip

subroutine MagSymm_k_Type_ctor(MagSymm_k_Type_param, Latt_param, SymopSymb_param, MSymopSymb_param, Numops_param, BNS_symbol_param, BNS_number_param, mcentred_param, MagModel_param, kvec_param, nmsym_param, MSymOp_param, Num_Lat_param, irrep_action_param, small_irrep_dim_param, Ltr_param, Parent_num_param, Sk_type_param, nirreps_param, MagType_param, OG_number_param, irrep_dim_param, irrep_modes_number_param, nkv_param, irrep_direction_param, centred_param, basf_param, Parent_spg_param, OG_symbol_param, SymOp_param, irrep_id_param, nbas_param, icomp_param, Multip_param)
	type (MagSymm_k_Type) :: MagSymm_k_Type_param
	character(len=1), intent(in) :: Latt_param
	character(len=40),   dimension(:),   allocatable, intent(in) :: SymopSymb_param
	character(len=40),   dimension(:,:), allocatable, intent(in) :: MSymopSymb_param
	integer, intent(in) :: Numops_param
	Character(len=34), intent(in) :: BNS_symbol_param
	character(len=15), intent(in) :: BNS_number_param
	integer, intent(in) :: mcentred_param
	character(len=31), intent(in) :: MagModel_param
	real(kind=cp),dimension(3,12), intent(in) :: kvec_param
	integer, intent(in) :: nmsym_param
	type(MSym_Oper_Type),dimension(:,:), allocatable, intent(in) :: MSymOp_param
	integer, intent(in) :: Num_Lat_param
	Character(len=20),   dimension(4), intent(in) :: irrep_action_param
	Integer,             dimension(4), intent(in) :: small_irrep_dim_param
	real(kind=cp), dimension(3,4), intent(in) :: Ltr_param
	Integer, intent(in) :: Parent_num_param
	character(len=15), intent(in) :: Sk_type_param
	integer, intent(in) :: nirreps_param
	Integer, intent(in) :: MagType_param
	character(len=15), intent(in) :: OG_number_param
	Integer,             dimension(4), intent(in) :: irrep_dim_param
	Integer,             dimension(4), intent(in) :: irrep_modes_number_param
	integer, intent(in) :: nkv_param
	Character(len=20),   dimension(4), intent(in) :: irrep_direction_param
	integer, intent(in) :: centred_param
	Complex(kind=cp),    dimension(3,12,48,4), intent(in) :: basf_param
	Character(len=20), intent(in) :: Parent_spg_param
	Character(len=34), intent(in) :: OG_symbol_param
	type(Sym_Oper_Type), dimension(:),   allocatable, intent(in) :: SymOp_param
	Character(len=15),   dimension(4), intent(in) :: irrep_id_param
	integer,             dimension(4), intent(in) :: nbas_param
	integer,             dimension(12,4), intent(in) :: icomp_param
	integer, intent(in) :: Multip_param
	MagSymm_k_Type_param%Latt = Latt_param
	MagSymm_k_Type_param%SymopSymb = SymopSymb_param
	MagSymm_k_Type_param%MSymopSymb = MSymopSymb_param
	MagSymm_k_Type_param%Numops = Numops_param
	MagSymm_k_Type_param%BNS_symbol = BNS_symbol_param
	MagSymm_k_Type_param%BNS_number = BNS_number_param
	MagSymm_k_Type_param%mcentred = mcentred_param
	MagSymm_k_Type_param%MagModel = MagModel_param
	MagSymm_k_Type_param%kvec = kvec_param
	MagSymm_k_Type_param%nmsym = nmsym_param
	MagSymm_k_Type_param%MSymOp = MSymOp_param
	MagSymm_k_Type_param%Num_Lat = Num_Lat_param
	MagSymm_k_Type_param%irrep_action = irrep_action_param
	MagSymm_k_Type_param%small_irrep_dim = small_irrep_dim_param
	MagSymm_k_Type_param%Ltr = Ltr_param
	MagSymm_k_Type_param%Parent_num = Parent_num_param
	MagSymm_k_Type_param%Sk_type = Sk_type_param
	MagSymm_k_Type_param%nirreps = nirreps_param
	MagSymm_k_Type_param%MagType = MagType_param
	MagSymm_k_Type_param%OG_number = OG_number_param
	MagSymm_k_Type_param%irrep_dim = irrep_dim_param
	MagSymm_k_Type_param%irrep_modes_number = irrep_modes_number_param
	MagSymm_k_Type_param%nkv = nkv_param
	MagSymm_k_Type_param%irrep_direction = irrep_direction_param
	MagSymm_k_Type_param%centred = centred_param
	MagSymm_k_Type_param%basf = basf_param
	MagSymm_k_Type_param%Parent_spg = Parent_spg_param
	MagSymm_k_Type_param%OG_symbol = OG_symbol_param
	MagSymm_k_Type_param%SymOp = SymOp_param
	MagSymm_k_Type_param%irrep_id = irrep_id_param
	MagSymm_k_Type_param%nbas = nbas_param
	MagSymm_k_Type_param%icomp = icomp_param
	MagSymm_k_Type_param%Multip = Multip_param
end subroutine MagSymm_k_Type_ctor
