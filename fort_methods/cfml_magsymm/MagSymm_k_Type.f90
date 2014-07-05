subroutine get_magsymm_k_Latt(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=1), intent(out) :: output_value
	output_value = obj_var%Latt
end subroutine get_magsymm_k_Latt

subroutine set_magsymm_k_Latt(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=1), intent(in) :: new_value
	obj_var%Latt = new_value
end subroutine set_magsymm_k_Latt

subroutine get_magsymm_k_SymopSymb(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:),   allocatable, intent(out) :: output_value
	output_value = obj_var%SymopSymb
end subroutine get_magsymm_k_SymopSymb

subroutine set_magsymm_k_SymopSymb(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:),   allocatable, intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine set_magsymm_k_SymopSymb

subroutine get_magsymm_k_MSymopSymb(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:,:), allocatable, intent(out) :: output_value
	output_value = obj_var%MSymopSymb
end subroutine get_magsymm_k_MSymopSymb

subroutine set_magsymm_k_MSymopSymb(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=40),   dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%MSymopSymb = new_value
end subroutine set_magsymm_k_MSymopSymb

function get_magsymm_k_Numops(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_Numops
	get_magsymm_k_Numops = obj_var%Numops
end function get_magsymm_k_Numops

subroutine set_magsymm_k_Numops(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Numops = new_value
end subroutine set_magsymm_k_Numops

subroutine get_magsymm_k_BNS_symbol(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(out) :: output_value
	output_value = obj_var%BNS_symbol
end subroutine get_magsymm_k_BNS_symbol

subroutine set_magsymm_k_BNS_symbol(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%BNS_symbol = new_value
end subroutine set_magsymm_k_BNS_symbol

subroutine get_magsymm_k_BNS_number(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(out) :: output_value
	output_value = obj_var%BNS_number
end subroutine get_magsymm_k_BNS_number

subroutine set_magsymm_k_BNS_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%BNS_number = new_value
end subroutine set_magsymm_k_BNS_number

function get_magsymm_k_mcentred(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_mcentred
	get_magsymm_k_mcentred = obj_var%mcentred
end function get_magsymm_k_mcentred

subroutine set_magsymm_k_mcentred(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%mcentred = new_value
end subroutine set_magsymm_k_mcentred

subroutine get_magsymm_k_MagModel(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=31), intent(out) :: output_value
	output_value = obj_var%MagModel
end subroutine get_magsymm_k_MagModel

subroutine set_magsymm_k_MagModel(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=31), intent(in) :: new_value
	obj_var%MagModel = new_value
end subroutine set_magsymm_k_MagModel

subroutine get_magsymm_k_kvec(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(out) :: output_value
	output_value = obj_var%kvec
end subroutine get_magsymm_k_kvec

subroutine set_magsymm_k_kvec(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp),dimension(3,12), intent(in) :: new_value
	obj_var%kvec = new_value
end subroutine set_magsymm_k_kvec

function get_magsymm_k_nmsym(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_nmsym
	get_magsymm_k_nmsym = obj_var%nmsym
end function get_magsymm_k_nmsym

subroutine set_magsymm_k_nmsym(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nmsym = new_value
end subroutine set_magsymm_k_nmsym

subroutine get_magsymm_k_MSymOp(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:,:), allocatable, intent(out) :: output_value
	output_value = obj_var%MSymOp
end subroutine get_magsymm_k_MSymOp

subroutine set_magsymm_k_MSymOp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:,:), allocatable, intent(in) :: new_value
	obj_var%MSymOp = new_value
end subroutine set_magsymm_k_MSymOp

function get_magsymm_k_Num_Lat(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_Num_Lat
	get_magsymm_k_Num_Lat = obj_var%Num_Lat
end function get_magsymm_k_Num_Lat

subroutine set_magsymm_k_Num_Lat(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Lat = new_value
end subroutine set_magsymm_k_Num_Lat

subroutine get_magsymm_k_irrep_action(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(out) :: output_value
	output_value = obj_var%irrep_action
end subroutine get_magsymm_k_irrep_action

subroutine set_magsymm_k_irrep_action(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(in) :: new_value
	obj_var%irrep_action = new_value
end subroutine set_magsymm_k_irrep_action

subroutine get_magsymm_k_small_irrep_dim(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(out) :: output_value
	output_value = obj_var%small_irrep_dim
end subroutine get_magsymm_k_small_irrep_dim

subroutine set_magsymm_k_small_irrep_dim(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%small_irrep_dim = new_value
end subroutine set_magsymm_k_small_irrep_dim

subroutine get_magsymm_k_Ltr(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp), dimension(3,4), intent(out) :: output_value
	output_value = obj_var%Ltr
end subroutine get_magsymm_k_Ltr

subroutine set_magsymm_k_Ltr(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	real(kind=cp), dimension(3,4), intent(in) :: new_value
	obj_var%Ltr = new_value
end subroutine set_magsymm_k_Ltr

function get_magsymm_k_Parent_num(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer :: get_magsymm_k_Parent_num
	get_magsymm_k_Parent_num = obj_var%Parent_num
end function get_magsymm_k_Parent_num

subroutine set_magsymm_k_Parent_num(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Parent_num = new_value
end subroutine set_magsymm_k_Parent_num

subroutine get_magsymm_k_Sk_type(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(out) :: output_value
	output_value = obj_var%Sk_type
end subroutine get_magsymm_k_Sk_type

subroutine set_magsymm_k_Sk_type(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%Sk_type = new_value
end subroutine set_magsymm_k_Sk_type

function get_magsymm_k_nirreps(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_nirreps
	get_magsymm_k_nirreps = obj_var%nirreps
end function get_magsymm_k_nirreps

subroutine set_magsymm_k_nirreps(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nirreps = new_value
end subroutine set_magsymm_k_nirreps

function get_magsymm_k_MagType(obj_var)
	type (MagSymm_k_Type) :: obj_var
	Integer :: get_magsymm_k_MagType
	get_magsymm_k_MagType = obj_var%MagType
end function get_magsymm_k_MagType

subroutine set_magsymm_k_MagType(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%MagType = new_value
end subroutine set_magsymm_k_MagType

subroutine get_magsymm_k_OG_number(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(out) :: output_value
	output_value = obj_var%OG_number
end subroutine get_magsymm_k_OG_number

subroutine set_magsymm_k_OG_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%OG_number = new_value
end subroutine set_magsymm_k_OG_number

subroutine get_magsymm_k_irrep_dim(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(out) :: output_value
	output_value = obj_var%irrep_dim
end subroutine get_magsymm_k_irrep_dim

subroutine set_magsymm_k_irrep_dim(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%irrep_dim = new_value
end subroutine set_magsymm_k_irrep_dim

subroutine get_magsymm_k_irrep_modes_number(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(out) :: output_value
	output_value = obj_var%irrep_modes_number
end subroutine get_magsymm_k_irrep_modes_number

subroutine set_magsymm_k_irrep_modes_number(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Integer,             dimension(4), intent(in) :: new_value
	obj_var%irrep_modes_number = new_value
end subroutine set_magsymm_k_irrep_modes_number

function get_magsymm_k_nkv(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_nkv
	get_magsymm_k_nkv = obj_var%nkv
end function get_magsymm_k_nkv

subroutine set_magsymm_k_nkv(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nkv = new_value
end subroutine set_magsymm_k_nkv

subroutine get_magsymm_k_irrep_direction(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(out) :: output_value
	output_value = obj_var%irrep_direction
end subroutine get_magsymm_k_irrep_direction

subroutine set_magsymm_k_irrep_direction(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20),   dimension(4), intent(in) :: new_value
	obj_var%irrep_direction = new_value
end subroutine set_magsymm_k_irrep_direction

function get_magsymm_k_centred(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_centred
	get_magsymm_k_centred = obj_var%centred
end function get_magsymm_k_centred

subroutine set_magsymm_k_centred(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%centred = new_value
end subroutine set_magsymm_k_centred

subroutine get_magsymm_k_basf(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Complex(kind=cp),    dimension(3,12,48,4), intent(out) :: output_value
	output_value = obj_var%basf
end subroutine get_magsymm_k_basf

subroutine set_magsymm_k_basf(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Complex(kind=cp),    dimension(3,12,48,4), intent(in) :: new_value
	obj_var%basf = new_value
end subroutine set_magsymm_k_basf

subroutine get_magsymm_k_Parent_spg(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20), intent(out) :: output_value
	output_value = obj_var%Parent_spg
end subroutine get_magsymm_k_Parent_spg

subroutine set_magsymm_k_Parent_spg(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=20), intent(in) :: new_value
	obj_var%Parent_spg = new_value
end subroutine set_magsymm_k_Parent_spg

subroutine get_magsymm_k_OG_symbol(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(out) :: output_value
	output_value = obj_var%OG_symbol
end subroutine get_magsymm_k_OG_symbol

subroutine set_magsymm_k_OG_symbol(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%OG_symbol = new_value
end subroutine set_magsymm_k_OG_symbol

subroutine get_magsymm_k_SymOp(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),   allocatable, intent(out) :: output_value
	output_value = obj_var%SymOp
end subroutine get_magsymm_k_SymOp

subroutine set_magsymm_k_SymOp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),   allocatable, intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine set_magsymm_k_SymOp

subroutine get_magsymm_k_irrep_id(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=15),   dimension(4), intent(out) :: output_value
	output_value = obj_var%irrep_id
end subroutine get_magsymm_k_irrep_id

subroutine set_magsymm_k_irrep_id(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	Character(len=15),   dimension(4), intent(in) :: new_value
	obj_var%irrep_id = new_value
end subroutine set_magsymm_k_irrep_id

subroutine get_magsymm_k_nbas(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(4), intent(out) :: output_value
	output_value = obj_var%nbas
end subroutine get_magsymm_k_nbas

subroutine set_magsymm_k_nbas(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(4), intent(in) :: new_value
	obj_var%nbas = new_value
end subroutine set_magsymm_k_nbas

subroutine get_magsymm_k_icomp(obj_var, output_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(12,4), intent(out) :: output_value
	output_value = obj_var%icomp
end subroutine get_magsymm_k_icomp

subroutine set_magsymm_k_icomp(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer,             dimension(12,4), intent(in) :: new_value
	obj_var%icomp = new_value
end subroutine set_magsymm_k_icomp

function get_magsymm_k_Multip(obj_var)
	type (MagSymm_k_Type) :: obj_var
	integer :: get_magsymm_k_Multip
	get_magsymm_k_Multip = obj_var%Multip
end function get_magsymm_k_Multip

subroutine set_magsymm_k_Multip(obj_var, new_value)
	type (MagSymm_k_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine set_magsymm_k_Multip

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
