function get_powder_numor_nbang(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_nbang
	get_powder_numor_nbang = obj_var%nbang
end function get_powder_numor_nbang

subroutine set_powder_numor_nbang(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nbang = new_value
end subroutine set_powder_numor_nbang

function get_powder_numor_manip(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_manip
	get_powder_numor_manip = obj_var%manip
end function get_powder_numor_manip

subroutine set_powder_numor_manip(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%manip = new_value
end subroutine set_powder_numor_manip

function get_powder_numor_icalc(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_icalc
	get_powder_numor_icalc = obj_var%icalc
end function get_powder_numor_icalc

subroutine set_powder_numor_icalc(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%icalc = new_value
end subroutine set_powder_numor_icalc

function get_powder_numor_monitor(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp) :: get_powder_numor_monitor
	get_powder_numor_monitor = obj_var%monitor
end function get_powder_numor_monitor

subroutine set_powder_numor_monitor(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%monitor = new_value
end subroutine set_powder_numor_monitor

function get_powder_numor_nbdata(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_nbdata
	get_powder_numor_nbdata = obj_var%nbdata
end function get_powder_numor_nbdata

subroutine set_powder_numor_nbdata(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nbdata = new_value
end subroutine set_powder_numor_nbdata

subroutine get_powder_numor_title(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=32), intent(out) :: output_value
	output_value = obj_var%title
end subroutine get_powder_numor_title

subroutine set_powder_numor_title(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=32), intent(in) :: new_value
	obj_var%title = new_value
end subroutine set_powder_numor_title

function get_powder_numor_numor(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_numor
	get_powder_numor_numor = obj_var%numor
end function get_powder_numor_numor

subroutine set_powder_numor_numor(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%numor = new_value
end subroutine set_powder_numor_numor

subroutine get_powder_numor_Instrm(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=12), intent(out) :: output_value
	output_value = obj_var%Instrm
end subroutine get_powder_numor_Instrm

subroutine set_powder_numor_Instrm(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=12), intent(in) :: new_value
	obj_var%Instrm = new_value
end subroutine set_powder_numor_Instrm

function get_powder_numor_wave(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp) :: get_powder_numor_wave
	get_powder_numor_wave = obj_var%wave
end function get_powder_numor_wave

subroutine set_powder_numor_wave(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%wave = new_value
end subroutine set_powder_numor_wave

subroutine get_powder_numor_header(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=32), intent(out) :: output_value
	output_value = obj_var%header
end subroutine get_powder_numor_header

subroutine set_powder_numor_header(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=32), intent(in) :: new_value
	obj_var%header = new_value
end subroutine set_powder_numor_header

subroutine get_powder_numor_Scantype(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=8), intent(out) :: output_value
	output_value = obj_var%Scantype
end subroutine get_powder_numor_Scantype

subroutine set_powder_numor_Scantype(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	character(len=8), intent(in) :: new_value
	obj_var%Scantype = new_value
end subroutine set_powder_numor_Scantype

subroutine get_powder_numor_icdesc(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, dimension(11), intent(out) :: output_value
	output_value = obj_var%icdesc
end subroutine get_powder_numor_icdesc

subroutine set_powder_numor_icdesc(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, dimension(11), intent(in) :: new_value
	obj_var%icdesc = new_value
end subroutine set_powder_numor_icdesc

subroutine get_powder_numor_angles(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(5), intent(out) :: output_value
	output_value = obj_var%angles
end subroutine get_powder_numor_angles

subroutine set_powder_numor_angles(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(5), intent(in) :: new_value
	obj_var%angles = new_value
end subroutine set_powder_numor_angles

function get_powder_numor_time(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp) :: get_powder_numor_time
	get_powder_numor_time = obj_var%time
end function get_powder_numor_time

subroutine set_powder_numor_time(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%time = new_value
end subroutine set_powder_numor_time

subroutine get_powder_numor_counts(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp),dimension(:,:), intent(out) :: output_value
	output_value = obj_var%counts
end subroutine get_powder_numor_counts

subroutine set_powder_numor_counts(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp),dimension(:,:), intent(in) :: new_value
	obj_var%counts = new_value
end subroutine set_powder_numor_counts

function get_powder_numor_nframes(obj_var)
	type (POWDER_Numor_Type) :: obj_var
	integer :: get_powder_numor_nframes
	get_powder_numor_nframes = obj_var%nframes
end function get_powder_numor_nframes

subroutine set_powder_numor_nframes(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%nframes = new_value
end subroutine set_powder_numor_nframes

subroutine get_powder_numor_conditions(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(5), intent(out) :: output_value
	output_value = obj_var%conditions
end subroutine get_powder_numor_conditions

subroutine set_powder_numor_conditions(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(5), intent(in) :: new_value
	obj_var%conditions = new_value
end subroutine set_powder_numor_conditions

subroutine get_powder_numor_tmc_ang(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp),dimension(:,:), intent(out) :: output_value
	output_value = obj_var%tmc_ang
end subroutine get_powder_numor_tmc_ang

subroutine set_powder_numor_tmc_ang(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp),dimension(:,:), intent(in) :: new_value
	obj_var%tmc_ang = new_value
end subroutine set_powder_numor_tmc_ang

subroutine get_powder_numor_scans(obj_var, output_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(3), intent(out) :: output_value
	output_value = obj_var%scans
end subroutine get_powder_numor_scans

subroutine set_powder_numor_scans(obj_var, new_value)
	type (POWDER_Numor_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%scans = new_value
end subroutine set_powder_numor_scans

subroutine POWDER_Numor_Type_ctor(POWDER_Numor_Type_param, nbang_param, manip_param, icalc_param, monitor_param, nbdata_param, title_param, numor_param, Instrm_param, wave_param, header_param, Scantype_param, icdesc_param, angles_param, time_param, counts_param, nframes_param, conditions_param, tmc_ang_param, scans_param)
	type (POWDER_Numor_Type) :: POWDER_Numor_Type_param
	integer, intent(in) :: nbang_param
	integer, intent(in) :: manip_param
	integer, intent(in) :: icalc_param
	real(kind=cp), intent(in) :: monitor_param
	integer, intent(in) :: nbdata_param
	character(len=32), intent(in) :: title_param
	integer, intent(in) :: numor_param
	character(len=12), intent(in) :: Instrm_param
	real(kind=cp), intent(in) :: wave_param
	character(len=32), intent(in) :: header_param
	character(len=8), intent(in) :: Scantype_param
	integer, dimension(11), intent(in) :: icdesc_param
	real(kind=cp), dimension(5), intent(in) :: angles_param
	real(kind=cp), intent(in) :: time_param
	real(kind=cp),dimension(:,:), intent(in) :: counts_param
	integer, intent(in) :: nframes_param
	real(kind=cp), dimension(5), intent(in) :: conditions_param
	real(kind=cp),dimension(:,:), intent(in) :: tmc_ang_param
	real(kind=cp), dimension(3), intent(in) :: scans_param
	POWDER_Numor_Type_param%nbang = nbang_param
	POWDER_Numor_Type_param%manip = manip_param
	POWDER_Numor_Type_param%icalc = icalc_param
	POWDER_Numor_Type_param%monitor = monitor_param
	POWDER_Numor_Type_param%nbdata = nbdata_param
	POWDER_Numor_Type_param%title = title_param
	POWDER_Numor_Type_param%numor = numor_param
	POWDER_Numor_Type_param%Instrm = Instrm_param
	POWDER_Numor_Type_param%wave = wave_param
	POWDER_Numor_Type_param%header = header_param
	POWDER_Numor_Type_param%Scantype = Scantype_param
	POWDER_Numor_Type_param%icdesc = icdesc_param
	POWDER_Numor_Type_param%angles = angles_param
	POWDER_Numor_Type_param%time = time_param
	POWDER_Numor_Type_param%counts = counts_param
	POWDER_Numor_Type_param%nframes = nframes_param
	POWDER_Numor_Type_param%conditions = conditions_param
	POWDER_Numor_Type_param%tmc_ang = tmc_ang_param
	POWDER_Numor_Type_param%scans = scans_param
end subroutine POWDER_Numor_Type_ctor
