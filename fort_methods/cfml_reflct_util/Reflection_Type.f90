function get_reflection_A(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getA
	getA = obj_var%A
end function get_reflection_A

subroutine set_reflection_A(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%A = new_value
end subroutine set_reflection_A

function get_reflection_AA(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getAA
	getAA = obj_var%AA
end function get_reflection_AA

subroutine set_reflection_AA(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%AA = new_value
end subroutine set_reflection_AA

function get_reflection_B(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getB
	getB = obj_var%B
end function get_reflection_B

subroutine set_reflection_B(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%B = new_value
end subroutine set_reflection_B

function get_reflection_BB(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getBB
	getBB = obj_var%BB
end function get_reflection_BB

subroutine set_reflection_BB(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%BB = new_value
end subroutine set_reflection_BB

function get_reflection_H(obj_var)
	type (Reflection_Type) :: obj_var
	integer,dimension(3) :: getH
	getH = obj_var%H
end function get_reflection_H

subroutine set_reflection_H(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	integer,dimension(3), intent(in) :: new_value
	obj_var%H = new_value
end subroutine set_reflection_H

function get_reflection_SFo(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getSFo
	getSFo = obj_var%SFo
end function get_reflection_SFo

subroutine set_reflection_SFo(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%SFo = new_value
end subroutine set_reflection_SFo

function get_reflection_S(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getS
	getS = obj_var%S
end function get_reflection_S

subroutine set_reflection_S(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%S = new_value
end subroutine set_reflection_S

function get_reflection_Fc(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getFc
	getFc = obj_var%Fc
end function get_reflection_Fc

subroutine set_reflection_Fc(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Fc = new_value
end subroutine set_reflection_Fc

function get_reflection_W(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getW
	getW = obj_var%W
end function get_reflection_W

subroutine set_reflection_W(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%W = new_value
end subroutine set_reflection_W

function get_reflection_Phase(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getPhase
	getPhase = obj_var%Phase
end function get_reflection_Phase

subroutine set_reflection_Phase(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Phase = new_value
end subroutine set_reflection_Phase

function get_reflection_Mult(obj_var)
	type (Reflection_Type) :: obj_var
	integer :: getMult
	getMult = obj_var%Mult
end function get_reflection_Mult

subroutine set_reflection_Mult(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	integer, intent(in) :: new_value
	obj_var%Mult = new_value
end subroutine set_reflection_Mult

function get_reflection_Fo(obj_var)
	type (Reflection_Type) :: obj_var
	real(kind=cp) :: getFo
	getFo = obj_var%Fo
end function get_reflection_Fo

subroutine set_reflection_Fo(obj_var, new_value)
	type (Reflection_Type) :: obj_var
	real(kind=cp), intent(in) :: new_value
	obj_var%Fo = new_value
end subroutine set_reflection_Fo

subroutine Reflection_Type_ctor(Reflection_Type_param, A_param, AA_param, B_param, BB_param, H_param, SFo_param, S_param, Fc_param, W_param, Phase_param, Mult_param, Fo_param)
	type (Reflection_Type) :: Reflection_Type_param
	real(kind=cp), intent(in) :: A_param
	real(kind=cp), intent(in) :: AA_param
	real(kind=cp), intent(in) :: B_param
	real(kind=cp), intent(in) :: BB_param
	integer,dimension(3), intent(in) :: H_param
	real(kind=cp), intent(in) :: SFo_param
	real(kind=cp), intent(in) :: S_param
	real(kind=cp), intent(in) :: Fc_param
	real(kind=cp), intent(in) :: W_param
	real(kind=cp), intent(in) :: Phase_param
	integer, intent(in) :: Mult_param
	real(kind=cp), intent(in) :: Fo_param
	Reflection_Type_param%A = A_param
	Reflection_Type_param%AA = AA_param
	Reflection_Type_param%B = B_param
	Reflection_Type_param%BB = BB_param
	Reflection_Type_param%H = H_param
	Reflection_Type_param%SFo = SFo_param
	Reflection_Type_param%S = S_param
	Reflection_Type_param%Fc = Fc_param
	Reflection_Type_param%W = W_param
	Reflection_Type_param%Phase = Phase_param
	Reflection_Type_param%Mult = Mult_param
	Reflection_Type_param%Fo = Fo_param
end subroutine Reflection_Type_ctor
