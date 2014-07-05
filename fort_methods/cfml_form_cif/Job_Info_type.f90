subroutine get_job_info_dtt1(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%dtt1
end subroutine get_job_info_dtt1

subroutine set_job_info_dtt1(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: new_value
	obj_var%dtt1 = new_value
end subroutine set_job_info_dtt1

subroutine get_job_info_dtt2(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%dtt2
end subroutine get_job_info_dtt2

subroutine set_job_info_dtt2(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: new_value
	obj_var%dtt2 = new_value
end subroutine set_job_info_dtt2

subroutine get_job_info_range_2theta(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_2theta
end subroutine get_job_info_range_2theta

subroutine set_job_info_range_2theta(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_2theta = new_value
end subroutine set_job_info_range_2theta

subroutine get_job_info_Title(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	character(len=120), intent(out) :: output_value
	output_value = obj_var%Title
end subroutine get_job_info_Title

subroutine set_job_info_Title(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=120), intent(in) :: new_value
	obj_var%Title = new_value
end subroutine set_job_info_Title

subroutine get_job_info_range_tof(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_tof
end subroutine get_job_info_range_tof

subroutine set_job_info_range_tof(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_tof = new_value
end subroutine set_job_info_range_tof

function get_job_info_Num_Phases(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: get_job_info_Num_Phases
	get_job_info_Num_Phases = obj_var%Num_Phases
end function get_job_info_Num_Phases

subroutine set_job_info_Num_Phases(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Phases = new_value
end subroutine set_job_info_Num_Phases

subroutine get_job_info_cmd(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%cmd
end subroutine get_job_info_cmd

subroutine set_job_info_cmd(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(in) :: new_value
	obj_var%cmd = new_value
end subroutine set_job_info_cmd

subroutine get_job_info_range_stl(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_stl
end subroutine get_job_info_range_stl

subroutine set_job_info_range_stl(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_stl = new_value
end subroutine set_job_info_range_stl

subroutine get_job_info_range_d(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_d
end subroutine get_job_info_range_d

subroutine set_job_info_range_d(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_d = new_value
end subroutine set_job_info_range_d

function get_job_info_Num_Patterns(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: get_job_info_Num_Patterns
	get_job_info_Num_Patterns = obj_var%Num_Patterns
end function get_job_info_Num_Patterns

subroutine set_job_info_Num_Patterns(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_Patterns = new_value
end subroutine set_job_info_Num_Patterns

subroutine get_job_info_Patt_typ(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	character(len=16),  dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%Patt_typ
end subroutine get_job_info_Patt_typ

subroutine set_job_info_Patt_typ(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=16),  dimension(:), allocatable, intent(in) :: new_value
	obj_var%Patt_typ = new_value
end subroutine set_job_info_Patt_typ

subroutine get_job_info_Phas_nam(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%Phas_nam
end subroutine get_job_info_Phas_nam

subroutine set_job_info_Phas_nam(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	character(len=128), dimension(:), allocatable, intent(in) :: new_value
	obj_var%Phas_nam = new_value
end subroutine set_job_info_Phas_nam

function get_job_info_Num_cmd(obj_var)
	type (Job_Info_type) :: obj_var
	integer :: get_job_info_Num_cmd
	get_job_info_Num_cmd = obj_var%Num_cmd
end function get_job_info_Num_cmd

subroutine set_job_info_Num_cmd(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Num_cmd = new_value
end subroutine set_job_info_Num_cmd

subroutine get_job_info_range_q(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_q
end subroutine get_job_info_range_q

subroutine set_job_info_range_q(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_q = new_value
end subroutine set_job_info_range_q

subroutine get_job_info_ratio(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%ratio
end subroutine get_job_info_ratio

subroutine set_job_info_ratio(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: new_value
	obj_var%ratio = new_value
end subroutine set_job_info_ratio

subroutine get_job_info_range_Energy(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%range_Energy
end subroutine get_job_info_range_Energy

subroutine set_job_info_range_Energy(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%range_Energy = new_value
end subroutine set_job_info_range_Energy

subroutine get_job_info_Lambda(obj_var, output_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(out) :: output_value
	output_value = obj_var%Lambda
end subroutine get_job_info_Lambda

subroutine set_job_info_Lambda(obj_var, new_value)
	type (Job_Info_type) :: obj_var
	type(interval_type),dimension(:), allocatable, intent(in) :: new_value
	obj_var%Lambda = new_value
end subroutine set_job_info_Lambda

subroutine Job_Info_type_ctor(Job_Info_type_param, dtt1_param, dtt2_param, range_2theta_param, Title_param, range_tof_param, Num_Phases_param, cmd_param, range_stl_param, range_d_param, Num_Patterns_param, Patt_typ_param, Phas_nam_param, Num_cmd_param, range_q_param, ratio_param, range_Energy_param, Lambda_param)
	type (Job_Info_type) :: Job_Info_type_param
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: dtt1_param
	real(kind=cp)      ,dimension(:), allocatable, intent(in) :: dtt2_param
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
	type(interval_type),dimension(:), allocatable, intent(in) :: range_Energy_param
	type(interval_type),dimension(:), allocatable, intent(in) :: Lambda_param
	Job_Info_type_param%dtt1 = dtt1_param
	Job_Info_type_param%dtt2 = dtt2_param
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
	Job_Info_type_param%range_Energy = range_Energy_param
	Job_Info_type_param%Lambda = Lambda_param
end subroutine Job_Info_type_ctor
