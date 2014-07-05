function get_ns_space_group_Centred(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_Centred
	get_ns_space_group_Centred = obj_var%Centred
end function get_ns_space_group_Centred

subroutine set_ns_space_group_Centred(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Centred = new_value
end subroutine set_ns_space_group_Centred

subroutine get_ns_space_group_SymopSymb(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=50),      dimension(:), intent(out) :: output_value
	output_value = obj_var%SymopSymb
end subroutine get_ns_space_group_SymopSymb

subroutine set_ns_space_group_SymopSymb(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=50),      dimension(:), intent(in) :: new_value
	obj_var%SymopSymb = new_value
end subroutine set_ns_space_group_SymopSymb

function get_ns_space_group_NumSpg(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_NumSpg
	get_ns_space_group_NumSpg = obj_var%NumSpg
end function get_ns_space_group_NumSpg

subroutine set_ns_space_group_NumSpg(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumSpg = new_value
end subroutine set_ns_space_group_NumSpg

subroutine get_ns_space_group_Latt_trans(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(out) :: output_value
	output_value = obj_var%Latt_trans
end subroutine get_ns_space_group_Latt_trans

subroutine set_ns_space_group_Latt_trans(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(:,:), intent(in) :: new_value
	obj_var%Latt_trans = new_value
end subroutine set_ns_space_group_Latt_trans

subroutine get_ns_space_group_CrystalSys(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=12), intent(out) :: output_value
	output_value = obj_var%CrystalSys
end subroutine get_ns_space_group_CrystalSys

subroutine set_ns_space_group_CrystalSys(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%CrystalSys = new_value
end subroutine set_ns_space_group_CrystalSys

function get_ns_space_group_NumLat(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_NumLat
	get_ns_space_group_NumLat = obj_var%NumLat
end function get_ns_space_group_NumLat

subroutine set_ns_space_group_NumLat(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumLat = new_value
end subroutine set_ns_space_group_NumLat

subroutine get_ns_space_group_PG(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(out) :: output_value
	output_value = obj_var%PG
end subroutine get_ns_space_group_PG

subroutine set_ns_space_group_PG(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%PG = new_value
end subroutine set_ns_space_group_PG

subroutine get_ns_space_group_Hall(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=16), intent(out) :: output_value
	output_value = obj_var%Hall
end subroutine get_ns_space_group_Hall

subroutine set_ns_space_group_Hall(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=16), intent(in) :: new_value
	obj_var%Hall = new_value
end subroutine set_ns_space_group_Hall

subroutine get_ns_space_group_Info(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(out) :: output_value
	output_value = obj_var%Info
end subroutine get_ns_space_group_Info

subroutine set_ns_space_group_Info(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Info = new_value
end subroutine set_ns_space_group_Info

subroutine get_ns_space_group_SPG_lat(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 1), intent(out) :: output_value
	output_value = obj_var%SPG_lat
end subroutine get_ns_space_group_SPG_lat

subroutine set_ns_space_group_SPG_lat(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 1), intent(in) :: new_value
	obj_var%SPG_lat = new_value
end subroutine set_ns_space_group_SPG_lat

subroutine get_ns_space_group_Laue(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(out) :: output_value
	output_value = obj_var%Laue
end subroutine get_ns_space_group_Laue

subroutine set_ns_space_group_Laue(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Laue = new_value
end subroutine set_ns_space_group_Laue

subroutine get_ns_space_group_SPG_latsy(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 2), intent(out) :: output_value
	output_value = obj_var%SPG_latsy
end subroutine get_ns_space_group_SPG_latsy

subroutine set_ns_space_group_SPG_latsy(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len= 2), intent(in) :: new_value
	obj_var%SPG_latsy = new_value
end subroutine set_ns_space_group_SPG_latsy

function get_ns_space_group_Num_gen(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_Num_gen
	get_ns_space_group_Num_gen = obj_var%Num_gen
end function get_ns_space_group_Num_gen

subroutine set_ns_space_group_Num_gen(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_gen = new_value
end subroutine set_ns_space_group_Num_gen

subroutine get_ns_space_group_Bravais(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=51), intent(out) :: output_value
	output_value = obj_var%Bravais
end subroutine get_ns_space_group_Bravais

subroutine set_ns_space_group_Bravais(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=51), intent(in) :: new_value
	obj_var%Bravais = new_value
end subroutine set_ns_space_group_Bravais

subroutine get_ns_space_group_SG_setting(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=90), intent(out) :: output_value
	output_value = obj_var%SG_setting
end subroutine get_ns_space_group_SG_setting

subroutine set_ns_space_group_SG_setting(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%SG_setting = new_value
end subroutine set_ns_space_group_SG_setting

subroutine get_ns_space_group_gHall(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=90), intent(out) :: output_value
	output_value = obj_var%gHall
end subroutine get_ns_space_group_gHall

subroutine set_ns_space_group_gHall(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%gHall = new_value
end subroutine set_ns_space_group_gHall

subroutine get_ns_space_group_SPG_Symb(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=20), intent(out) :: output_value
	output_value = obj_var%SPG_Symb
end subroutine get_ns_space_group_SPG_Symb

subroutine set_ns_space_group_SPG_Symb(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%SPG_Symb = new_value
end subroutine set_ns_space_group_SPG_Symb

subroutine get_ns_space_group_Centre_coord(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%Centre_coord
end subroutine get_ns_space_group_Centre_coord

subroutine set_ns_space_group_Centre_coord(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Centre_coord = new_value
end subroutine set_ns_space_group_Centre_coord

subroutine get_ns_space_group_SymOp(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	type(NS_Sym_Oper_Type), dimension(:), intent(out) :: output_value
	output_value = obj_var%SymOp
end subroutine get_ns_space_group_SymOp

subroutine set_ns_space_group_SymOp(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	type(NS_Sym_Oper_Type), dimension(:), intent(in) :: new_value
	obj_var%SymOp = new_value
end subroutine set_ns_space_group_SymOp

subroutine get_ns_space_group_Centre(obj_var, output_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=80), intent(out) :: output_value
	output_value = obj_var%Centre
end subroutine get_ns_space_group_Centre

subroutine set_ns_space_group_Centre(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	character(len=80), intent(in) :: new_value
	obj_var%Centre = new_value
end subroutine set_ns_space_group_Centre

function get_ns_space_group_NumOps(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_NumOps
	get_ns_space_group_NumOps = obj_var%NumOps
end function get_ns_space_group_NumOps

subroutine set_ns_space_group_NumOps(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumOps = new_value
end subroutine set_ns_space_group_NumOps

function get_ns_space_group_Multip(obj_var)
	type (NS_Space_Group_Type) :: obj_var
	integer :: get_ns_space_group_Multip
	get_ns_space_group_Multip = obj_var%Multip
end function get_ns_space_group_Multip

subroutine set_ns_space_group_Multip(obj_var, new_value)
	type (NS_Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Multip = new_value
end subroutine set_ns_space_group_Multip

subroutine NS_Space_Group_Type_ctor(NS_Space_Group_Type_param, Centred_param, SymopSymb_param, NumSpg_param, Latt_trans_param, CrystalSys_param, NumLat_param, PG_param, Hall_param, Info_param, SPG_lat_param, Laue_param, SPG_latsy_param, Num_gen_param, Bravais_param, SG_setting_param, gHall_param, SPG_Symb_param, Centre_coord_param, SymOp_param, Centre_param, NumOps_param, Multip_param)
	type (NS_Space_Group_Type) :: NS_Space_Group_Type_param
	integer, intent(in) :: Centred_param
	character(len=50),      dimension(:), intent(in) :: SymopSymb_param
	integer, intent(in) :: NumSpg_param
	real(kind=cp), dimension(:,:), intent(in) :: Latt_trans_param
	character(len=12), intent(in) :: CrystalSys_param
	integer, intent(in) :: NumLat_param
	character(len= 5), intent(in) :: PG_param
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
	real(kind=cp), dimension(3), intent(in) :: Centre_coord_param
	type(NS_Sym_Oper_Type), dimension(:), intent(in) :: SymOp_param
	character(len=80), intent(in) :: Centre_param
	integer, intent(in) :: NumOps_param
	integer, intent(in) :: Multip_param
	NS_Space_Group_Type_param%Centred = Centred_param
	NS_Space_Group_Type_param%SymopSymb = SymopSymb_param
	NS_Space_Group_Type_param%NumSpg = NumSpg_param
	NS_Space_Group_Type_param%Latt_trans = Latt_trans_param
	NS_Space_Group_Type_param%CrystalSys = CrystalSys_param
	NS_Space_Group_Type_param%NumLat = NumLat_param
	NS_Space_Group_Type_param%PG = PG_param
	NS_Space_Group_Type_param%Hall = Hall_param
	NS_Space_Group_Type_param%Info = Info_param
	NS_Space_Group_Type_param%SPG_lat = SPG_lat_param
	NS_Space_Group_Type_param%Laue = Laue_param
	NS_Space_Group_Type_param%SPG_latsy = SPG_latsy_param
	NS_Space_Group_Type_param%Num_gen = Num_gen_param
	NS_Space_Group_Type_param%Bravais = Bravais_param
	NS_Space_Group_Type_param%SG_setting = SG_setting_param
	NS_Space_Group_Type_param%gHall = gHall_param
	NS_Space_Group_Type_param%SPG_Symb = SPG_Symb_param
	NS_Space_Group_Type_param%Centre_coord = Centre_coord_param
	NS_Space_Group_Type_param%SymOp = SymOp_param
	NS_Space_Group_Type_param%Centre = Centre_param
	NS_Space_Group_Type_param%NumOps = NumOps_param
	NS_Space_Group_Type_param%Multip = Multip_param
end subroutine NS_Space_Group_Type_ctor
