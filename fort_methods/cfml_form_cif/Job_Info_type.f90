function getrange_2theta(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_2theta
	getrange_2theta = obj_var%range_2theta
end function getrange_2theta

subroutine setrange_2theta(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_2theta = new_value
end subroutine setrange_2theta

function getTitle(obj_var)
	type (Job_Info_type) :: obj_var
	character(len=120) :: getTitle
	getTitle = obj_var%Title
end function getTitle

subroutine setTitle(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=120), intent(in) :: new_value
	obj_var%Title = new_value
end subroutine setTitle

function getrange_tof(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_tof
	getrange_tof = obj_var%range_tof
end function getrange_tof

subroutine setrange_tof(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_tof = new_value
end subroutine setrange_tof

function getNum_Phases(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: getNum_Phases
	getNum_Phases = obj_var%Num_Phases
end function getNum_Phases

subroutine setNum_Phases(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Phases = new_value
end subroutine setNum_Phases

function getcmd(obj_var)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable :: getcmd
	getcmd = obj_var%cmd
end function getcmd

subroutine setcmd(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(in) :: new_value
	obj_var%cmd = new_value
end subroutine setcmd

function getrange_stl(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_stl
	getrange_stl = obj_var%range_stl
end function getrange_stl

subroutine setrange_stl(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_stl = new_value
end subroutine setrange_stl

function getrange_d(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_d
	getrange_d = obj_var%range_d
end function getrange_d

subroutine setrange_d(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_d = new_value
end subroutine setrange_d

function getNum_Patterns(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: getNum_Patterns
	getNum_Patterns = obj_var%Num_Patterns
end function getNum_Patterns

subroutine setNum_Patterns(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Patterns = new_value
end subroutine setNum_Patterns

function getPatt_typ(obj_var)
	type (Job_Info_type) :: obj_var
	character(len=16),  dimension(:), allocatable :: getPatt_typ
	getPatt_typ = obj_var%Patt_typ
end function getPatt_typ

subroutine setPatt_typ(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=16),  dimension(:), allocatable, intent(in) :: new_value
	obj_var%Patt_typ = new_value
end subroutine setPatt_typ

function getPhas_nam(obj_var)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable :: getPhas_nam
	getPhas_nam = obj_var%Phas_nam
end function getPhas_nam

subroutine setPhas_nam(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(in) :: new_value
	obj_var%Phas_nam = new_value
end subroutine setPhas_nam

function getNum_cmd(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: getNum_cmd
	getNum_cmd = obj_var%Num_cmd
end function getNum_cmd

subroutine setNum_cmd(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_cmd = new_value
end subroutine setNum_cmd

function getrange_q(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_q
	getrange_q = obj_var%range_q
end function getrange_q

subroutine setrange_q(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_q = new_value
end subroutine setrange_q

function getratio(obj_var)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable :: getratio
	getratio = obj_var%ratio
end function getratio

subroutine setratio(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: new_value
	obj_var%ratio = new_value
end subroutine setratio

function getdtt1,dtt2(obj_var)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable :: getdtt1,dtt2
	getdtt1,dtt2 = obj_var%dtt1,dtt2
end function getdtt1,dtt2

subroutine setdtt1,dtt2(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: new_value
	obj_var%dtt1,dtt2 = new_value
end subroutine setdtt1,dtt2

function getrange_Energy(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getrange_Energy
	getrange_Energy = obj_var%range_Energy
end function getrange_Energy

subroutine setrange_Energy(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_Energy = new_value
end subroutine setrange_Energy

function getLambda(obj_var)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable :: getLambda
	getLambda = obj_var%Lambda
end function getLambda

subroutine setLambda(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%Lambda = new_value
end subroutine setLambda

subroutine Job_Info_type_ctor(Job_Info_type_param, range_2theta_param, Title_param, range_tof_param, Num_Phases_param, cmd_param, range_stl_param, range_d_param, Num_Patterns_param, Patt_typ_param, Phas_nam_param, Num_cmd_param, range_q_param, ratio_param, dtt1,dtt2_param, range_Energy_param, Lambda_param)
	type (Job_Info_type) :: Job_Info_type_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_2theta_param
	character(len=120), intent(in) :: Title_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_tof_param
	integer, intent(in) :: Num_Phases_param
	character(len=128), dimension(:), allocatable, intent(in) :: cmd_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_stl_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_d_param
	integer, intent(in) :: Num_Patterns_param
	character(len=16),  dimension(:), allocatable, intent(in) :: Patt_typ_param
	character(len=128), dimension(:), allocatable, intent(in) :: Phas_nam_param
	integer, intent(in) :: Num_cmd_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_q_param
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: ratio_param
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: dtt1,dtt2_param
	type(interval_type),dimension(:), allocatable, intent(in) :: range_Energy_param
	type(interval_type),dimension(:), allocatable, intent(in) :: Lambda_param
	Job_Info_type_param%range_2theta = range_2theta_param
	Job_Info_type_param%Title = Title_param
	Job_Info_type_param%range_tof = range_tof_param
	Job_Info_type_param%Num_Phases = Num_Phases_param
	Job_Info_type_param%cmd = cmd_param
	Job_Info_type_param%range_stl = range_stl_param
	Job_Info_type_param%range_d = range_d_param
	Job_Info_type_param%Num_Patterns = Num_Patterns_param
	Job_Info_type_param%Patt_typ = Patt_typ_param
	Job_Info_type_param%Phas_nam = Phas_nam_param
	Job_Info_type_param%Num_cmd = Num_cmd_param
	Job_Info_type_param%range_q = range_q_param
	Job_Info_type_param%ratio = ratio_param
	Job_Info_type_param%dtt1,dtt2 = dtt1,dtt2_param
	Job_Info_type_param%range_Energy = range_Energy_param
	Job_Info_type_param%Lambda = Lambda_param
end subroutine Job_Info_type_ctor
