function get_diffraction_pattern_scal(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_scal
	get_diffraction_pattern_scal = obj_var%scal
end function get_diffraction_pattern_scal

subroutine set_diffraction_pattern_scal(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%scal = new_value
end subroutine set_diffraction_pattern_scal

subroutine get_diffraction_pattern_conv(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (5), intent(out) :: output_value
	output_value = obj_var%conv
end subroutine get_diffraction_pattern_conv

subroutine set_diffraction_pattern_conv(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (5), intent(in) :: new_value
	obj_var%conv = new_value
end subroutine set_diffraction_pattern_conv

subroutine get_diffraction_pattern_instr(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(out) :: output_value
	output_value = obj_var%instr
end subroutine get_diffraction_pattern_instr

subroutine set_diffraction_pattern_instr(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%instr = new_value
end subroutine set_diffraction_pattern_instr

subroutine get_diffraction_pattern_yax_text(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=30), intent(out) :: output_value
	output_value = obj_var%yax_text
end subroutine get_diffraction_pattern_yax_text

subroutine set_diffraction_pattern_yax_text(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=30), intent(in) :: new_value
	obj_var%yax_text = new_value
end subroutine set_diffraction_pattern_yax_text

subroutine get_diffraction_pattern_diff_kind(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(out) :: output_value
	output_value = obj_var%diff_kind
end subroutine get_diffraction_pattern_diff_kind

subroutine set_diffraction_pattern_diff_kind(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%diff_kind = new_value
end subroutine set_diffraction_pattern_diff_kind

function get_diffraction_pattern_xmin(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_xmin
	get_diffraction_pattern_xmin = obj_var%xmin
end function get_diffraction_pattern_xmin

subroutine set_diffraction_pattern_xmin(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%xmin = new_value
end subroutine set_diffraction_pattern_xmin

function get_diffraction_pattern_norm_mon(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_norm_mon
	get_diffraction_pattern_norm_mon = obj_var%norm_mon
end function get_diffraction_pattern_norm_mon

subroutine set_diffraction_pattern_norm_mon(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%norm_mon = new_value
end subroutine set_diffraction_pattern_norm_mon

function get_diffraction_pattern_ymin(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_ymin
	get_diffraction_pattern_ymin = obj_var%ymin
end function get_diffraction_pattern_ymin

subroutine set_diffraction_pattern_ymin(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%ymin = new_value
end subroutine set_diffraction_pattern_ymin

function get_diffraction_pattern_monitor(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_monitor
	get_diffraction_pattern_monitor = obj_var%monitor
end function get_diffraction_pattern_monitor

subroutine set_diffraction_pattern_monitor(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%monitor = new_value
end subroutine set_diffraction_pattern_monitor

subroutine get_diffraction_pattern_filepath(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=512), intent(out) :: output_value
	output_value = obj_var%filepath
end subroutine get_diffraction_pattern_filepath

subroutine set_diffraction_pattern_filepath(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=512), intent(in) :: new_value
	obj_var%filepath = new_value
end subroutine set_diffraction_pattern_filepath

subroutine get_diffraction_pattern_Title(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=180), intent(out) :: output_value
	output_value = obj_var%Title
end subroutine get_diffraction_pattern_Title

subroutine set_diffraction_pattern_Title(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=180), intent(in) :: new_value
	obj_var%Title = new_value
end subroutine set_diffraction_pattern_Title

subroutine get_diffraction_pattern_nd(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	integer,       dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%nd
end subroutine get_diffraction_pattern_nd

subroutine set_diffraction_pattern_nd(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	integer,       dimension (:), allocatable, intent(in) :: new_value
	obj_var%nd = new_value
end subroutine set_diffraction_pattern_nd

subroutine get_diffraction_pattern_filename(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=512), intent(out) :: output_value
	output_value = obj_var%filename
end subroutine get_diffraction_pattern_filename

subroutine set_diffraction_pattern_filename(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=512), intent(in) :: new_value
	obj_var%filename = new_value
end subroutine set_diffraction_pattern_filename

subroutine get_diffraction_pattern_ycalc(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%ycalc
end subroutine get_diffraction_pattern_ycalc

subroutine set_diffraction_pattern_ycalc(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(in) :: new_value
	obj_var%ycalc = new_value
end subroutine set_diffraction_pattern_ycalc

subroutine get_diffraction_pattern_scat_var(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(out) :: output_value
	output_value = obj_var%scat_var
end subroutine get_diffraction_pattern_scat_var

subroutine set_diffraction_pattern_scat_var(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=20), intent(in) :: new_value
	obj_var%scat_var = new_value
end subroutine set_diffraction_pattern_scat_var

function get_diffraction_pattern_step(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_step
	get_diffraction_pattern_step = obj_var%step
end function get_diffraction_pattern_step

subroutine set_diffraction_pattern_step(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%step = new_value
end subroutine set_diffraction_pattern_step

subroutine get_diffraction_pattern_x(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%x
end subroutine get_diffraction_pattern_x

subroutine set_diffraction_pattern_x(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(in) :: new_value
	obj_var%x = new_value
end subroutine set_diffraction_pattern_x

function get_diffraction_pattern_ct_step(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	logical :: get_diffraction_pattern_ct_step
	get_diffraction_pattern_ct_step = obj_var%ct_step
end function get_diffraction_pattern_ct_step

subroutine set_diffraction_pattern_ct_step(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	logical, intent(in) :: new_value
	obj_var%ct_step = new_value
end subroutine set_diffraction_pattern_ct_step

subroutine get_diffraction_pattern_xax_text(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=30), intent(out) :: output_value
	output_value = obj_var%xax_text
end subroutine get_diffraction_pattern_xax_text

subroutine set_diffraction_pattern_xax_text(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	character(len=30), intent(in) :: new_value
	obj_var%xax_text = new_value
end subroutine set_diffraction_pattern_xax_text

subroutine get_diffraction_pattern_istat(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	integer,       dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%istat
end subroutine get_diffraction_pattern_istat

subroutine set_diffraction_pattern_istat(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	integer,       dimension (:), allocatable, intent(in) :: new_value
	obj_var%istat = new_value
end subroutine set_diffraction_pattern_istat

function get_diffraction_pattern_ymax(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_ymax
	get_diffraction_pattern_ymax = obj_var%ymax
end function get_diffraction_pattern_ymax

subroutine set_diffraction_pattern_ymax(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%ymax = new_value
end subroutine set_diffraction_pattern_ymax

function get_diffraction_pattern_col_time(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_col_time
	get_diffraction_pattern_col_time = obj_var%col_time
end function get_diffraction_pattern_col_time

subroutine set_diffraction_pattern_col_time(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%col_time = new_value
end subroutine set_diffraction_pattern_col_time

function get_diffraction_pattern_Tset(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_Tset
	get_diffraction_pattern_Tset = obj_var%Tset
end function get_diffraction_pattern_Tset

subroutine set_diffraction_pattern_Tset(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Tset = new_value
end subroutine set_diffraction_pattern_Tset

function get_diffraction_pattern_Tsamp(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_Tsamp
	get_diffraction_pattern_Tsamp = obj_var%Tsamp
end function get_diffraction_pattern_Tsamp

subroutine set_diffraction_pattern_Tsamp(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Tsamp = new_value
end subroutine set_diffraction_pattern_Tsamp

subroutine get_diffraction_pattern_bgr(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%bgr
end subroutine get_diffraction_pattern_bgr

subroutine set_diffraction_pattern_bgr(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(in) :: new_value
	obj_var%bgr = new_value
end subroutine set_diffraction_pattern_bgr

function get_diffraction_pattern_xmax(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp) :: get_diffraction_pattern_xmax
	get_diffraction_pattern_xmax = obj_var%xmax
end function get_diffraction_pattern_xmax

subroutine set_diffraction_pattern_xmax(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%xmax = new_value
end subroutine set_diffraction_pattern_xmax

subroutine get_diffraction_pattern_y(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%y
end subroutine get_diffraction_pattern_y

subroutine set_diffraction_pattern_y(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(in) :: new_value
	obj_var%y = new_value
end subroutine set_diffraction_pattern_y

function get_diffraction_pattern_npts(obj_var)
	type (Diffraction_Pattern_Type) :: obj_var
	integer :: get_diffraction_pattern_npts
	get_diffraction_pattern_npts = obj_var%npts
end function get_diffraction_pattern_npts

subroutine set_diffraction_pattern_npts(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%npts = new_value
end subroutine set_diffraction_pattern_npts

subroutine get_diffraction_pattern_sigma(obj_var, output_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(out) :: output_value
	output_value = obj_var%sigma
end subroutine get_diffraction_pattern_sigma

subroutine set_diffraction_pattern_sigma(obj_var, new_value)
	type (Diffraction_Pattern_Type) :: obj_var
	real(kind=cp), dimension (:), allocatable, intent(in) :: new_value
	obj_var%sigma = new_value
end subroutine set_diffraction_pattern_sigma

subroutine Diffraction_Pattern_Type_ctor(Diffraction_Pattern_Type_param, scal_param, conv_param, instr_param, yax_text_param, diff_kind_param, xmin_param, norm_mon_param, ymin_param, monitor_param, filepath_param, Title_param, nd_param, filename_param, ycalc_param, scat_var_param, step_param, x_param, ct_step_param, xax_text_param, istat_param, ymax_param, col_time_param, Tset_param, Tsamp_param, bgr_param, xmax_param, y_param, npts_param, sigma_param)
	type (Diffraction_Pattern_Type) :: Diffraction_Pattern_Type_param
	real(kind=cp), intent(in) :: scal_param
	real(kind=cp), dimension (5), intent(in) :: conv_param
	character(len=20), intent(in) :: instr_param
	character(len=30), intent(in) :: yax_text_param
	character(len=20), intent(in) :: diff_kind_param
	real(kind=cp), intent(in) :: xmin_param
	real(kind=cp), intent(in) :: norm_mon_param
	real(kind=cp), intent(in) :: ymin_param
	real(kind=cp), intent(in) :: monitor_param
	character(len=512), intent(in) :: filepath_param
	character(len=180), intent(in) :: Title_param
	integer,       dimension (:), allocatable, intent(in) :: nd_param
	character(len=512), intent(in) :: filename_param
	real(kind=cp), dimension (:), allocatable, intent(in) :: ycalc_param
	character(len=20), intent(in) :: scat_var_param
	real(kind=cp), intent(in) :: step_param
	real(kind=cp), dimension (:), allocatable, intent(in) :: x_param
	logical, intent(in) :: ct_step_param
	character(len=30), intent(in) :: xax_text_param
	integer,       dimension (:), allocatable, intent(in) :: istat_param
	real(kind=cp), intent(in) :: ymax_param
	real(kind=cp), intent(in) :: col_time_param
	real(kind=cp), intent(in) :: Tset_param
	real(kind=cp), intent(in) :: Tsamp_param
	real(kind=cp), dimension (:), allocatable, intent(in) :: bgr_param
	real(kind=cp), intent(in) :: xmax_param
	real(kind=cp), dimension (:), allocatable, intent(in) :: y_param
	integer, intent(in) :: npts_param
	real(kind=cp), dimension (:), allocatable, intent(in) :: sigma_param
	Diffraction_Pattern_Type_param%scal = scal_param
	Diffraction_Pattern_Type_param%conv = conv_param
	Diffraction_Pattern_Type_param%instr = instr_param
	Diffraction_Pattern_Type_param%yax_text = yax_text_param
	Diffraction_Pattern_Type_param%diff_kind = diff_kind_param
	Diffraction_Pattern_Type_param%xmin = xmin_param
	Diffraction_Pattern_Type_param%norm_mon = norm_mon_param
	Diffraction_Pattern_Type_param%ymin = ymin_param
	Diffraction_Pattern_Type_param%monitor = monitor_param
	Diffraction_Pattern_Type_param%filepath = filepath_param
	Diffraction_Pattern_Type_param%Title = Title_param
	Diffraction_Pattern_Type_param%nd = nd_param
	Diffraction_Pattern_Type_param%filename = filename_param
	Diffraction_Pattern_Type_param%ycalc = ycalc_param
	Diffraction_Pattern_Type_param%scat_var = scat_var_param
	Diffraction_Pattern_Type_param%step = step_param
	Diffraction_Pattern_Type_param%x = x_param
	Diffraction_Pattern_Type_param%ct_step = ct_step_param
	Diffraction_Pattern_Type_param%xax_text = xax_text_param
	Diffraction_Pattern_Type_param%istat = istat_param
	Diffraction_Pattern_Type_param%ymax = ymax_param
	Diffraction_Pattern_Type_param%col_time = col_time_param
	Diffraction_Pattern_Type_param%Tset = Tset_param
	Diffraction_Pattern_Type_param%Tsamp = Tsamp_param
	Diffraction_Pattern_Type_param%bgr = bgr_param
	Diffraction_Pattern_Type_param%xmax = xmax_param
	Diffraction_Pattern_Type_param%y = y_param
	Diffraction_Pattern_Type_param%npts = npts_param
	Diffraction_Pattern_Type_param%sigma = sigma_param
end subroutine Diffraction_Pattern_Type_ctor
