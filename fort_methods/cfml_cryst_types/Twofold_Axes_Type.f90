function geta(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3) :: geta
	geta = obj_var%a
end function geta

subroutine seta(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(3), intent(in) :: new_value
	obj_var%a = new_value
end subroutine seta

function getmaxes(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12) :: getmaxes
	getmaxes = obj_var%maxes
end function getmaxes

subroutine setmaxes(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(in) :: new_value
	obj_var%maxes = new_value
end subroutine setmaxes

function getntwo(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	integer :: getntwo
	getntwo = obj_var%ntwo
end function getntwo

subroutine setntwo(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%ntwo = new_value
end subroutine setntwo

function getcross(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12) :: getcross
	getcross = obj_var%cross
end function getcross

subroutine setcross(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), dimension(12), intent(in) :: new_value
	obj_var%cross = new_value
end subroutine setcross

function gettol(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) :: gettol
	gettol = obj_var%tol
end function gettol

subroutine settol(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%tol = new_value
end subroutine settol

function getrtwofold(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12) :: getrtwofold
	getrtwofold = obj_var%rtwofold
end function getrtwofold

subroutine setrtwofold(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(in) :: new_value
	obj_var%rtwofold = new_value
end subroutine setrtwofold

function getcaxes(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) ,dimension(3,12) :: getcaxes
	getcaxes = obj_var%caxes
end function getcaxes

subroutine setcaxes(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	real(kind=cp) ,dimension(3,12), intent(in) :: new_value
	obj_var%caxes = new_value
end subroutine setcaxes

function getdot(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(12) :: getdot
	getdot = obj_var%dot
end function getdot

subroutine setdot(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(12), intent(in) :: new_value
	obj_var%dot = new_value
end subroutine setdot

function getdtwofold(obj_var)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12) :: getdtwofold
	getdtwofold = obj_var%dtwofold
end function getdtwofold

subroutine setdtwofold(obj_var, new_value)
	type (Twofold_Axes_Type) :: obj_var
	integer,dimension(3,12), intent(in) :: new_value
	obj_var%dtwofold = new_value
end subroutine setdtwofold

subroutine Twofold_Axes_Type_ctor(Twofold_Axes_Type_param, a_param, maxes_param, ntwo_param, cross_param, tol_param, rtwofold_param, caxes_param, dot_param, dtwofold_param)
	type (Twofold_Axes_Type) :: Twofold_Axes_Type_param
	real(kind=cp), dimension(3), intent(in) :: a_param
	real(kind=cp), dimension(12), intent(in) :: maxes_param
	integer, intent(in) :: ntwo_param
	real(kind=cp), dimension(12), intent(in) :: cross_param
	real(kind=cp), intent(in) :: tol_param
	integer,dimension(3,12), intent(in) :: rtwofold_param
	real(kind=cp) ,dimension(3,12), intent(in) :: caxes_param
	integer,dimension(12), intent(in) :: dot_param
	integer,dimension(3,12), intent(in) :: dtwofold_param
	Twofold_Axes_Type_param%a = a_param
	Twofold_Axes_Type_param%maxes = maxes_param
	Twofold_Axes_Type_param%ntwo = ntwo_param
	Twofold_Axes_Type_param%cross = cross_param
	Twofold_Axes_Type_param%tol = tol_param
	Twofold_Axes_Type_param%rtwofold = rtwofold_param
	Twofold_Axes_Type_param%caxes = caxes_param
	Twofold_Axes_Type_param%dot = dot_param
	Twofold_Axes_Type_param%dtwofold = dtwofold_param
end subroutine Twofold_Axes_Type_ctor
