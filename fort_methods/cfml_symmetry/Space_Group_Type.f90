function getgHall            (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=90) :: getgHall            
	getgHall             = obj_var%gHall            
end function getgHall            

subroutine setgHall            (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%gHall             = new_value
end subroutine setgHall            

function getSPG_latsy        (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 2) :: getSPG_latsy        
	getSPG_latsy         = obj_var%SPG_latsy        
end function getSPG_latsy        

subroutine setSPG_latsy        (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 2), intent(in) :: new_value
	obj_var%SPG_latsy         = new_value
end subroutine setSPG_latsy        

function getSG_setting       (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=90) :: getSG_setting       
	getSG_setting        = obj_var%SG_setting       
end function getSG_setting       

subroutine setSG_setting       (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=90), intent(in) :: new_value
	obj_var%SG_setting        = new_value
end subroutine setSG_setting       

function getSPG_lat          (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 1) :: getSPG_lat          
	getSPG_lat           = obj_var%SPG_lat          
end function getSPG_lat          

subroutine setSPG_lat          (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 1), intent(in) :: new_value
	obj_var%SPG_lat           = new_value
end subroutine setSPG_lat          

function getNumSpg           (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumSpg           
	getNumSpg            = obj_var%NumSpg           
end function getNumSpg           

subroutine setNumSpg           (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumSpg            = new_value
end subroutine setNumSpg           

function getMultip           (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getMultip           
	getMultip            = obj_var%Multip           
end function getMultip           

subroutine setMultip           (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Multip            = new_value
end subroutine setMultip           

function getInfo             (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getInfo             
	getInfo              = obj_var%Info             
end function getInfo             

subroutine setInfo             (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Info              = new_value
end subroutine setInfo             

function getR_Asym_Unit      (obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp),dimension(3,2) :: getR_Asym_Unit      
	getR_Asym_Unit       = obj_var%R_Asym_Unit      
end function getR_Asym_Unit      

subroutine setR_Asym_Unit      (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp),dimension(3,2), intent(in) :: new_value
	obj_var%R_Asym_Unit       = new_value
end subroutine setR_Asym_Unit      

function getHall             (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=16) :: getHall             
	getHall              = obj_var%Hall             
end function getHall             

subroutine setHall             (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=16), intent(in) :: new_value
	obj_var%Hall              = new_value
end subroutine setHall             

function getPG               (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getPG               
	getPG                = obj_var%PG               
end function getPG               

subroutine setPG               (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%PG                = new_value
end subroutine setPG               

function getNumOps           (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumOps           
	getNumOps            = obj_var%NumOps           
end function getNumOps           

subroutine setNumOps           (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumOps            = new_value
end subroutine setNumOps           

function getCentred          (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getCentred          
	getCentred           = obj_var%Centred          
end function getCentred          

subroutine setCentred          (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Centred           = new_value
end subroutine setCentred          

function getWyckoff          (obj_var)
	type (Space_Group_Type) :: obj_var
	type(Wyckoff_Type) :: getWyckoff          
	getWyckoff           = obj_var%Wyckoff          
end function getWyckoff          

subroutine setWyckoff          (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	type(Wyckoff_Type), intent(in) :: new_value
	obj_var%Wyckoff           = new_value
end subroutine setWyckoff          

function getHexa             (obj_var)
	type (Space_Group_Type) :: obj_var
	logical :: getHexa             
	getHexa              = obj_var%Hexa             
end function getHexa             

subroutine setHexa             (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%Hexa              = new_value
end subroutine setHexa             

function getSPG_Symb         (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=20) :: getSPG_Symb         
	getSPG_Symb          = obj_var%SPG_Symb         
end function getSPG_Symb         

subroutine setSPG_Symb         (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%SPG_Symb          = new_value
end subroutine setSPG_Symb         

function getNumLat           (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNumLat           
	getNumLat            = obj_var%NumLat           
end function getNumLat           

subroutine setNumLat           (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%NumLat            = new_value
end subroutine setNumLat           

function getCrystalSys       (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=12) :: getCrystalSys       
	getCrystalSys        = obj_var%CrystalSys       
end function getCrystalSys       

subroutine setCrystalSys       (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%CrystalSys        = new_value
end subroutine setCrystalSys       

function getBravais          (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=51) :: getBravais          
	getBravais           = obj_var%Bravais          
end function getBravais          

subroutine setBravais          (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=51), intent(in) :: new_value
	obj_var%Bravais           = new_value
end subroutine setBravais          

function getLaue             (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len= 5) :: getLaue             
	getLaue              = obj_var%Laue             
end function getLaue             

subroutine setLaue             (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len= 5), intent(in) :: new_value
	obj_var%Laue              = new_value
end subroutine setLaue             

function getNum_gen          (obj_var)
	type (Space_Group_Type) :: obj_var
	integer :: getNum_gen          
	getNum_gen           = obj_var%Num_gen          
end function getNum_gen          

subroutine setNum_gen          (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_gen           = new_value
end subroutine setNum_gen          

function getCentre           (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=80) :: getCentre           
	getCentre            = obj_var%Centre           
end function getCentre           

subroutine setCentre           (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=80), intent(in) :: new_value
	obj_var%Centre            = new_value
end subroutine setCentre           

function getLatt_trans       (obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:) :: getLatt_trans       
	getLatt_trans        = obj_var%Latt_trans       
end function getLatt_trans       

subroutine setLatt_trans       (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: new_value
	obj_var%Latt_trans        = new_value
end subroutine setLatt_trans       

function getCentre_coord     (obj_var)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3) :: getCentre_coord     
	getCentre_coord      = obj_var%Centre_coord     
end function getCentre_coord     

subroutine setCentre_coord     (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%Centre_coord      = new_value
end subroutine setCentre_coord     

function getSymOp            (obj_var)
	type (Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), allocatable,dimension(:) :: getSymOp            
	getSymOp             = obj_var%SymOp            
end function getSymOp            

subroutine setSymOp            (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	type(Sym_Oper_Type), allocatable,dimension(:), intent(in) :: new_value
	obj_var%SymOp             = new_value
end subroutine setSymOp            

function getSymopSymb        (obj_var)
	type (Space_Group_Type) :: obj_var
	character(len=50),   allocatable,dimension(:) :: getSymopSymb        
	getSymopSymb         = obj_var%SymopSymb        
end function getSymopSymb        

subroutine setSymopSymb        (obj_var, new_value)
	type (Space_Group_Type) :: obj_var
	character(len=50),   allocatable,dimension(:), intent(in) :: new_value
	obj_var%SymopSymb         = new_value
end subroutine setSymopSymb        

subroutine Space_Group_Type_ctor(Space_Group_Type_param, gHall            _param, SPG_latsy        _param, SG_setting       _param, SPG_lat          _param, NumSpg           _param, Multip           _param, Info             _param, R_Asym_Unit      _param, Hall             _param, PG               _param, NumOps           _param, Centred          _param, Wyckoff          _param, Hexa             _param, SPG_Symb         _param, NumLat           _param, CrystalSys       _param, Bravais          _param, Laue             _param, Num_gen          _param, Centre           _param, Latt_trans       _param, Centre_coord     _param, SymOp            _param, SymopSymb        _param)
	type (Space_Group_Type) :: Space_Group_Type_param
	character(len=90), intent(in) :: gHall            _param
	character(len= 2), intent(in) :: SPG_latsy        _param
	character(len=90), intent(in) :: SG_setting       _param
	character(len= 1), intent(in) :: SPG_lat          _param
	integer, intent(in) :: NumSpg           _param
	integer, intent(in) :: Multip           _param
	character(len= 5), intent(in) :: Info             _param
	real(kind=cp),dimension(3,2), intent(in) :: R_Asym_Unit      _param
	character(len=16), intent(in) :: Hall             _param
	character(len= 5), intent(in) :: PG               _param
	integer, intent(in) :: NumOps           _param
	integer, intent(in) :: Centred          _param
	type(Wyckoff_Type), intent(in) :: Wyckoff          _param
	logical, intent(in) :: Hexa             _param
	character(len=20), intent(in) :: SPG_Symb         _param
	integer, intent(in) :: NumLat           _param
	character(len=12), intent(in) :: CrystalSys       _param
	character(len=51), intent(in) :: Bravais          _param
	character(len= 5), intent(in) :: Laue             _param
	integer, intent(in) :: Num_gen          _param
	character(len=80), intent(in) :: Centre           _param
	real(kind=cp), allocatable,dimension(:,:), intent(in) :: Latt_trans       _param
	real(kind=cp), dimension(3), intent(in) :: Centre_coord     _param
	type(Sym_Oper_Type), allocatable,dimension(:), intent(in) :: SymOp            _param
	character(len=50),   allocatable,dimension(:), intent(in) :: SymopSymb        _param
	Space_Group_Type_param%gHall             = gHall            _param
	Space_Group_Type_param%SPG_latsy         = SPG_latsy        _param
	Space_Group_Type_param%SG_setting        = SG_setting       _param
	Space_Group_Type_param%SPG_lat           = SPG_lat          _param
	Space_Group_Type_param%NumSpg            = NumSpg           _param
	Space_Group_Type_param%Multip            = Multip           _param
	Space_Group_Type_param%Info              = Info             _param
	Space_Group_Type_param%R_Asym_Unit       = R_Asym_Unit      _param
	Space_Group_Type_param%Hall              = Hall             _param
	Space_Group_Type_param%PG                = PG               _param
	Space_Group_Type_param%NumOps            = NumOps           _param
	Space_Group_Type_param%Centred           = Centred          _param
	Space_Group_Type_param%Wyckoff           = Wyckoff          _param
	Space_Group_Type_param%Hexa              = Hexa             _param
	Space_Group_Type_param%SPG_Symb          = SPG_Symb         _param
	Space_Group_Type_param%NumLat            = NumLat           _param
	Space_Group_Type_param%CrystalSys        = CrystalSys       _param
	Space_Group_Type_param%Bravais           = Bravais          _param
	Space_Group_Type_param%Laue              = Laue             _param
	Space_Group_Type_param%Num_gen           = Num_gen          _param
	Space_Group_Type_param%Centre            = Centre           _param
	Space_Group_Type_param%Latt_trans        = Latt_trans       _param
	Space_Group_Type_param%Centre_coord      = Centre_coord     _param
	Space_Group_Type_param%SymOp             = SymOp            _param
	Space_Group_Type_param%SymopSymb         = SymopSymb        _param
end subroutine Space_Group_Type_ctor
