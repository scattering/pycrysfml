#!/usr/bin/python
# Joseph Lesniewski - NIST Center for Neutron Research
# Summer 2014
# global deps modifications for use with modified FortWrap
lines = [line.strip() for line in open("./CFML_GlobalDeps_Linux.f90")]
out = open("./CFML_GlobalDeps_Linux.f90", "w")
for line in lines:
	if 'integer,  parameter :: cp = sp' in line:
		out.write('integer, parameter :: cp = 4!sp')
	elif 'integer, parameter :: sp = selected_real_kind(6,30)' in line:
		out.write('integer, parameter :: sp = 4!selected_real_kind( 6, 30)')
	elif 'integer, parameter :: dp = selected_real_kind(14,150)' in line:
		out.write('integer, parameter :: dp = 8!selected_real_kind(14,150)')
	else:
		out.write(line)
	out.write('\n')
out.close()
# snb4c mod
lines = [line.strip() for line in open("./CFML_SXTAL_Geom.f90")]
out = open("./CFML_SXTAL_Geom.f90", "w")
for line in lines:
	if 'real(kind=cp), dimension(3), intent(in )  :: angl_nb(3)' in line:
		out.write('real(kind=cp), dimension(3), intent(in )  :: angl_nb !(3)')
	if 'real(kind=cp), dimension(4), intent(out)  :: angl_4c(4)' in line:
		out.write('real(kind=cp), dimension(4), intent(out)  :: angl_4c !(4)')
	else:
		out.write(line)
	out.write('\n')
out.close()
