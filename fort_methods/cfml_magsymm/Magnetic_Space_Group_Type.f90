function get_magnetic_space_group_m_constr(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: get_magnetic_space_group_m_constr
	get_magnetic_space_group_m_constr = obj_var%m_constr
end function get_magnetic_space_group_m_constr

subroutine set_magnetic_space_group_m_constr(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%m_constr = new_value
end subroutine set_magnetic_space_group_m_constr

subroutine get_magnetic_space_group_Centre_coord(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%Centre_coord
end subroutine get_magnetic_space_group_Centre_coord

subroutine set_magnetic_space_group_Centre_coord(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Centre_coord = new_value
end subroutine set_magnetic_space_group_Centre_coord

subroutine get_magnetic_space_group_aLatt_trans(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(out) :: output_value
	output_value = obj_var%aLatt_trans
end subroutine get_magnetic_space_group_aLatt_trans

subroutine set_magnetic_space_group_aLatt_trans(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(in) :: new_value
	obj_var%aLatt_trans = new_value
end subroutine set_magnetic_space_group_aLatt_trans

subroutine get_magnetic_space_group_MSymopSymb(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%MSymopSymb
end subroutine get_magnetic_space_group_MSymopSymb

subroutine set_magnetic_space_group_MSymopSymb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%MSymopSymb = new_value
end subroutine set_magnetic_space_group_MSymopSymb

function get_magnetic_space_group_n_wyck(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_n_wyck
	get_magnetic_space_group_n_wyck = obj_var%n_wyck
end function get_magnetic_space_group_n_wyck

subroutine set_magnetic_space_group_n_wyck(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_wyck = new_value
end subroutine set_magnetic_space_group_n_wyck

subroutine get_magnetic_space_group_BNS_symbol(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(out) :: output_value
	output_value = obj_var%BNS_symbol
end subroutine get_magnetic_space_group_BNS_symbol

subroutine set_magnetic_space_group_BNS_symbol(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%BNS_symbol = new_value
end subroutine set_magnetic_space_group_BNS_symbol

subroutine get_magnetic_space_group_BNS_number(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(out) :: output_value
	output_value = obj_var%BNS_number
end subroutine get_magnetic_space_group_BNS_number

subroutine set_magnetic_space_group_BNS_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%BNS_number = new_value
end subroutine set_magnetic_space_group_BNS_number

subroutine get_magnetic_space_group_Latt_trans(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(out) :: output_value
	output_value = obj_var%Latt_trans
end subroutine get_magnetic_space_group_Latt_trans

subroutine set_magnetic_space_group_Latt_trans(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(in) :: new_value
	obj_var%Latt_trans = new_value
end subroutine set_magnetic_space_group_Latt_trans

function get_magnetic_space_group_Sh_number(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_Sh_number
	get_magnetic_space_group_Sh_number = obj_var%Sh_number
end function get_magnetic_space_group_Sh_number

subroutine set_magnetic_space_group_Sh_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Sh_number = new_value
end subroutine set_magnetic_space_group_Sh_number

subroutine get_magnetic_space_group_Wyck_Symb(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%Wyck_Symb
end subroutine get_magnetic_space_group_Wyck_Symb

subroutine set_magnetic_space_group_Wyck_Symb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%Wyck_Symb = new_value
end subroutine set_magnetic_space_group_Wyck_Symb

subroutine get_magnetic_space_group_CrystalSys(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=12), intent(out) :: output_value
	output_value = obj_var%CrystalSys
end subroutine get_magnetic_space_group_CrystalSys

subroutine set_magnetic_space_group_CrystalSys(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%CrystalSys = new_value
end subroutine set_magnetic_space_group_CrystalSys

function get_magnetic_space_group_m_cell(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: get_magnetic_space_group_m_cell
	get_magnetic_space_group_m_cell = obj_var%m_cell
end function get_magnetic_space_group_m_cell

subroutine set_magnetic_space_group_m_cell(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%m_cell = new_value
end subroutine set_magnetic_space_group_m_cell

subroutine get_magnetic_space_group_irrep_id(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%irrep_id
end subroutine get_magnetic_space_group_irrep_id

subroutine set_magnetic_space_group_irrep_id(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_id = new_value
end subroutine set_magnetic_space_group_irrep_id

subroutine get_magnetic_space_group_MSymOp(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%MSymOp
end subroutine get_magnetic_space_group_MSymOp

subroutine set_magnetic_space_group_MSymOp(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(MSym_Oper_Type),dimension(:),allocatable, intent(in) :: new_value
	obj_var%MSymOp = new_value
end subroutine set_magnetic_space_group_MSymOp

function get_magnetic_space_group_Num_Lat(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: get_magnetic_space_group_Num_Lat
	get_magnetic_space_group_Num_Lat = obj_var%Num_Lat
end function get_magnetic_space_group_Num_Lat

subroutine set_magnetic_space_group_Num_Lat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Lat = new_value
end subroutine set_magnetic_space_group_Num_Lat

function get_magnetic_space_group_Num_aLat(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: get_magnetic_space_group_Num_aLat
	get_magnetic_space_group_Num_aLat = obj_var%Num_aLat
end function get_magnetic_space_group_Num_aLat

subroutine set_magnetic_space_group_Num_aLat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_aLat = new_value
end subroutine set_magnetic_space_group_Num_aLat

subroutine get_magnetic_space_group_small_irrep_dim(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%small_irrep_dim
end subroutine get_magnetic_space_group_small_irrep_dim

subroutine set_magnetic_space_group_small_irrep_dim(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%small_irrep_dim = new_value
end subroutine set_magnetic_space_group_small_irrep_dim

function get_magnetic_space_group_Centred(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: get_magnetic_space_group_Centred
	get_magnetic_space_group_Centred = obj_var%Centred
end function get_magnetic_space_group_Centred

subroutine set_magnetic_space_group_Centred(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Centred = new_value
end subroutine set_magnetic_space_group_Centred

subroutine get_magnetic_space_group_SPG_lat(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 1), intent(out) :: output_value
	output_value = obj_var%SPG_lat
end subroutine get_magnetic_space_group_SPG_lat

subroutine set_magnetic_space_group_SPG_lat(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 1), intent(in) :: new_value
	obj_var%SPG_lat = new_value
end subroutine set_magnetic_space_group_SPG_lat

subroutine get_magnetic_space_group_SymopSymb(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%SymopSymb
end subroutine get_magnetic_space_group_SymopSymb

subroutine set_magnetic_space_group_SymopSymb(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=40),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine set_magnetic_space_group_SymopSymb

function get_magnetic_space_group_Parent_num(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_Parent_num
	get_magnetic_space_group_Parent_num = obj_var%Parent_num
end function get_magnetic_space_group_Parent_num

subroutine set_magnetic_space_group_Parent_num(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Parent_num = new_value
end subroutine set_magnetic_space_group_Parent_num

subroutine get_magnetic_space_group_kv_label(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%kv_label
end subroutine get_magnetic_space_group_kv_label

subroutine set_magnetic_space_group_kv_label(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=15),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%kv_label = new_value
end subroutine set_magnetic_space_group_kv_label

function get_magnetic_space_group_MagType(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_MagType
	get_magnetic_space_group_MagType = obj_var%MagType
end function get_magnetic_space_group_MagType

subroutine set_magnetic_space_group_MagType(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%MagType = new_value
end subroutine set_magnetic_space_group_MagType

subroutine get_magnetic_space_group_irrep_action(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%irrep_action
end subroutine get_magnetic_space_group_irrep_action

subroutine set_magnetic_space_group_irrep_action(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_action = new_value
end subroutine set_magnetic_space_group_irrep_action

function get_magnetic_space_group_n_kv(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_n_kv
	get_magnetic_space_group_n_kv = obj_var%n_kv
end function get_magnetic_space_group_n_kv

subroutine set_magnetic_space_group_n_kv(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_kv = new_value
end subroutine set_magnetic_space_group_n_kv

subroutine get_magnetic_space_group_OG_number(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(out) :: output_value
	output_value = obj_var%OG_number
end subroutine get_magnetic_space_group_OG_number

subroutine set_magnetic_space_group_OG_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=15), intent(in) :: new_value
	obj_var%OG_number = new_value
end subroutine set_magnetic_space_group_OG_number

function get_magnetic_space_group_n_irreps(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_n_irreps
	get_magnetic_space_group_n_irreps = obj_var%n_irreps
end function get_magnetic_space_group_n_irreps

subroutine set_magnetic_space_group_n_irreps(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%n_irreps = new_value
end subroutine set_magnetic_space_group_n_irreps

subroutine get_magnetic_space_group_SPG_latsy(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 2), intent(out) :: output_value
	output_value = obj_var%SPG_latsy
end subroutine get_magnetic_space_group_SPG_latsy

subroutine set_magnetic_space_group_SPG_latsy(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len= 2), intent(in) :: new_value
	obj_var%SPG_latsy = new_value
end subroutine set_magnetic_space_group_SPG_latsy

function get_magnetic_space_group_Num_gen(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: get_magnetic_space_group_Num_gen
	get_magnetic_space_group_Num_gen = obj_var%Num_gen
end function get_magnetic_space_group_Num_gen

subroutine set_magnetic_space_group_Num_gen(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_gen = new_value
end subroutine set_magnetic_space_group_Num_gen

subroutine get_magnetic_space_group_irrep_direction(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%irrep_direction
end subroutine get_magnetic_space_group_irrep_direction

subroutine set_magnetic_space_group_irrep_direction(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20),   dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_direction = new_value
end subroutine set_magnetic_space_group_irrep_direction

function get_magnetic_space_group_standard_setting(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: get_magnetic_space_group_standard_setting
	get_magnetic_space_group_standard_setting = obj_var%standard_setting
end function get_magnetic_space_group_standard_setting

subroutine set_magnetic_space_group_standard_setting(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%standard_setting = new_value
end subroutine set_magnetic_space_group_standard_setting

subroutine get_magnetic_space_group_SymOp(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%SymOp
end subroutine get_magnetic_space_group_SymOp

subroutine set_magnetic_space_group_SymOp(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), dimension(:),allocatable, intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine set_magnetic_space_group_SymOp

subroutine get_magnetic_space_group_irrep_dim(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%irrep_dim
end subroutine get_magnetic_space_group_irrep_dim

subroutine set_magnetic_space_group_irrep_dim(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_dim = new_value
end subroutine set_magnetic_space_group_irrep_dim

subroutine get_magnetic_space_group_Parent_spg(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20), intent(out) :: output_value
	output_value = obj_var%Parent_spg
end subroutine get_magnetic_space_group_Parent_spg

subroutine set_magnetic_space_group_Parent_spg(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=20), intent(in) :: new_value
	obj_var%Parent_spg = new_value
end subroutine set_magnetic_space_group_Parent_spg

subroutine get_magnetic_space_group_OG_symbol(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(out) :: output_value
	output_value = obj_var%OG_symbol
end subroutine get_magnetic_space_group_OG_symbol

subroutine set_magnetic_space_group_OG_symbol(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=34), intent(in) :: new_value
	obj_var%OG_symbol = new_value
end subroutine set_magnetic_space_group_OG_symbol

subroutine get_magnetic_space_group_Centre(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=80), intent(out) :: output_value
	output_value = obj_var%Centre
end subroutine get_magnetic_space_group_Centre

subroutine set_magnetic_space_group_Centre(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	character(len=80), intent(in) :: new_value
	obj_var%Centre = new_value
end subroutine set_magnetic_space_group_Centre

function get_magnetic_space_group_mcif(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical :: get_magnetic_space_group_mcif
	get_magnetic_space_group_mcif = obj_var%mcif
end function get_magnetic_space_group_mcif

subroutine set_magnetic_space_group_mcif(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%mcif = new_value
end subroutine set_magnetic_space_group_mcif

function get_magnetic_space_group_NumOps(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer :: get_magnetic_space_group_NumOps
	get_magnetic_space_group_NumOps = obj_var%NumOps
end function get_magnetic_space_group_NumOps

subroutine set_magnetic_space_group_NumOps(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumOps = new_value
end subroutine set_magnetic_space_group_NumOps

subroutine get_magnetic_space_group_irrep_modes_number(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(out) :: output_value
	output_value = obj_var%irrep_modes_number
end subroutine get_magnetic_space_group_irrep_modes_number

subroutine set_magnetic_space_group_irrep_modes_number(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer,             dimension(:),allocatable, intent(in) :: new_value
	obj_var%irrep_modes_number = new_value
end subroutine set_magnetic_space_group_irrep_modes_number

subroutine get_magnetic_space_group_kv(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp),     dimension(:,:),allocatable, intent(out) :: output_value
	output_value = obj_var%kv
end subroutine get_magnetic_space_group_kv

subroutine set_magnetic_space_group_kv(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	real(kind=cp),     dimension(:,:),allocatable, intent(in) :: new_value
	obj_var%kv = new_value
end subroutine set_magnetic_space_group_kv

subroutine get_magnetic_space_group_trn_to_standard(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(out) :: output_value
	output_value = obj_var%trn_to_standard
end subroutine get_magnetic_space_group_trn_to_standard

subroutine set_magnetic_space_group_trn_to_standard(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(in) :: new_value
	obj_var%trn_to_standard = new_value
end subroutine set_magnetic_space_group_trn_to_standard

subroutine get_magnetic_space_group_trn_from_parent(obj_var, output_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(out) :: output_value
	output_value = obj_var%trn_from_parent
end subroutine get_magnetic_space_group_trn_from_parent

subroutine set_magnetic_space_group_trn_from_parent(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Character(len=40), intent(in) :: new_value
	obj_var%trn_from_parent = new_value
end subroutine set_magnetic_space_group_trn_from_parent

function get_magnetic_space_group_Multip(obj_var)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer :: get_magnetic_space_group_Multip
	get_magnetic_space_group_Multip = obj_var%Multip
end function get_magnetic_space_group_Multip

subroutine set_magnetic_space_group_Multip(obj_var, new_value)
	type (Magnetic_Space_Group_Type) :: obj_var
	Integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine set_magnetic_space_group_Multip

subroutine Magnetic_Space_Group_Type_ctor(Magnetic_Space_Group_Type_param, m_constr_param, Centre_coord_param, aLatt_trans_param, MSymopSymb_param, n_wyck_param, BNS_symbol_param, BNS_number_param, Latt_trans_param, Sh_number_param, Wyck_Symb_param, CrystalSys_param, m_cell_param, irrep_id_param, MSymOp_param, Num_Lat_param, Num_aLat_param, small_irrep_dim_param, Centred_param, SPG_lat_param, SymopSymb_param, Parent_num_param, kv_label_param, MagType_param, irrep_action_param, n_kv_param, OG_number_param, n_irreps_param, SPG_latsy_param, Num_gen_param, irrep_direction_param, standard_setting_param, SymOp_param, irrep_dim_param, Parent_spg_param, OG_symbol_param, Centre_param, mcif_param, NumOps_param, irrep_modes_number_param, kv_param, trn_to_standard_param, trn_from_parent_param, Multip_param)
	type (Magnetic_Space_Group_Type) :: Magnetic_Space_Group_Type_param
	logical, intent(in) :: m_constr_param
	real(kind=cp), dimension(3), intent(in) :: Centre_coord_param
	real(kind=cp), dimension(:,:), intent(in) :: aLatt_trans_param
	character(len=40),   dimension(:),allocatable, intent(in) :: MSymopSymb_param
	Integer, intent(in) :: n_wyck_param
	Character(len=34), intent(in) :: BNS_symbol_param
	character(len=15), intent(in) :: BNS_number_param
	real(kind=cp), dimension(:,:), intent(in) :: Latt_trans_param
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
