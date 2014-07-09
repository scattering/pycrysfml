!!-------------------------------------------------------
!!---- crystallographic fortran modules library (crysfml)
!!-------------------------------------------------------
!!---- the crysfml project is distributed under lgpl. in agreement with the
!!---- intergovernmental convention of the ill,  this software cannot be used
!!---- in military applications.
!!----
!!---- copyright (c) 1999-2012  institut laue-langevin (ill),  grenoble,  france
!!----                          universidad de la laguna (ull),  tenerife,  spain
!!----                          laboratoire leon brillouin(llb),  saclay,  france
!!----
!!---- authors: juan rodriguez-carvajal (ill)
!!----          javier gonzalez-platas  (ull)
!!----
!!---- contributors: laurent chapon     (ill)
!!----               marc janoschek     (los alamos national laboratory,  usa)
!!----               oksana zaharko     (paul scherrer institute,  switzerland)
!!----               tierry roisnel     (cdifx, rennes france)
!!----               eric pellegrini    (ill)
!!----
!!---- this library is free software; you can redistribute it and/or
!!---- modify it under the terms of the gnu lesser general public
!!---- license as published by the free software foundation; either
!!---- version 3.0 of the license,  or (at your option) any later version.
!!----
!!---- this library is distributed in the hope that it will be useful,
!!---- but without any warranty; without even the implied warranty of
!!---- merchantability or fitness for a particular purpose.  see the gnu
!!---- lesser general public license for more details.
!!----
!!---- you should have received a copy of the gnu lesser general public
!!---- license along with this library; if not,  see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- module: cfml_globaldeps (windows version)
!!----   info: precision for crysfml library and operating system information
!!----         all the global variables defined in this module are implicitly public.
!!----
!!---- history
!!--..    update: january - 2009
!!--..
!!---- variables
!!--..
!!--..    operating system
!!--..
!!----    ops
!!----    ops_name
!!----    ops_sep
!!--..
!!--..    precision data
!!--..
!!----    sp
!!----    dp
!!----    cp
!!--..
!!--..    trigonometric
!!--..
!!----    pi
!!----    to_deg
!!----    to_rad
!!----    tpi
!!--..
!!--..    numeric
!!--..
!!----    deps
!!----    eps
!!--..
!!---- functions
!!--..
!!----    directory_exists
!!----
!!---- subroutines
!!--..
!!----    write_date_time
!!----
!!
module cfml_globaldeps

!---- variables ----!
implicit none

public

!------------------------------------!
!---- operating system variables ----!
!------------------------------------!

!!----
!!---- ops
!!----   integer variable 1: windows,  2: linux,  3: macos,  ....
!!----   this is a variable set by the user of the library for the case
!!----   that there is no external library with a procedure for getting
!!----   the operating system.
!!----
!!---- update: march 2009
!!
integer,  parameter :: ops= 1    ! windows

!!----
!!---- ops_name
!!----   character variable containing the name of the operating system
!!----   this is a variable set by the user of the library for the case
!!----   that there is no external library with a procedure for getting
!!----   the operating system.
!!----
!!---- update: march 2009
!!
character(len=*),  parameter :: ops_name="windows"

!!----
!!---- ops_sep
!!----   ascii code of directory separator character
!!----   here it is written explicitly as a character variable
!!----
!!---- update: march 2009
!!
character(len=*),  parameter :: ops_sep="\"

!------------------------------!
!---- precision parameters ----!
!------------------------------!

!!----
!!---- sp
!!----    sp: single precision ( sp = selected_real_kind(6, 30) )
!!----
!!---- update: january - 2009
!!
integer,  parameter :: sp = selected_real_kind(6, 30)

!!----
!!---- dp
!!----    dp: double precision ( dp = selected_real_kind(14, 150) )
!!----
!!---- update: january - 2009
!!
integer,  parameter :: dp = selected_real_kind(14, 150)

!!----
!!---- cp
!!----    cp: current precision
!!----
!!---- update: january - 2009
!!
integer,  parameter :: cp = sp

!----------------------------------!
!---- trigonometric parameters ----!
!----------------------------------!

!!----
!!---- pi
!!----    real(kind=dp),  parameter ::  pi = 3.141592653589793238463_dp
!!----
!!----    pi value
!!----
!!---- update: january - 2009
!!
real(kind=dp),  parameter ::  pi = 3.141592653589793238463_dp

!!----
!!---- to_deg
!!----    real(kind=dp),  parameter ::  to_deg = 180.0_dp/pi
!!----
!!----    conversion from radians to degrees
!!----
!!---- update: january - 2009
!!
real(kind=dp),  parameter ::  to_deg  = 180.0_dp/pi

!!----
!!---- to_rad
!!----    real(kind=dp),  parameter ::  to_rad  = pi/180.0_dp
!!----
!!----    conversion from degrees to radians
!!----
!!---- update: january - 2009
!!
real(kind=dp),  parameter ::  to_rad  = pi/180.0_dp

!!----
!!---- tpi
!!----  real(kind=dp),  parameter ::  tpi = 6.283185307179586476925_dp
!!----
!!----  2.0*pi value
!!----
!!---- update: january - 2009
!!
real(kind=dp),  parameter ::  tpi = 6.283185307179586476925_dp

!----------------------------!
!---- numeric parameters ----!
!----------------------------!

!!----
!!---- deps
!!----    real(kind=dp),  parameter :: deps=0.00000001_dp
!!----
!!----    epsilon value use for comparison of real numbers
!!----
!!---- update: january - 2009
!!
real(kind=dp),  parameter,  public :: deps=0.00000001_dp

!!----
!!----  eps
!!----     real(kind=cp),  public ::  eps=0.00001_cp
!!----
!!----     epsilon value use for comparison of real numbers
!!----
!!----  update: january - 2009
!!
real(kind=cp),   parameter,  public  ::  eps=0.00001_cp

contains

!-------------------!
!---- functions ----!
!-------------------!

!!----
!!---- function directory_exists(dirname) result(info)
!!----    character(len=*),  intent(in) :: dirname
!!----    logical                      :: info
!!----
!!---- generic function dependent of the compiler that return
!!---- a logical value if a directory exists or not.
!!----
!!---- update: april - 2009
!!
function directory_exists(dirname) result(info)
!---- argument ----!
character(len=*),  intent(in) :: dirname
logical                      :: info

!---- local variables ----!
character(len=512) :: linea
integer            :: nlong

! init value
info=.false.

linea=trim(dirname)
nlong=len_trim(linea)
if (nlong ==0) return

if (linea(nlong:nlong) /= ops_sep) linea=trim(linea)//ops_sep

! all compilers except intel
!inquire(file=trim(linea)//'.' ,  exist=info)

! intel
inquire(directory=trim(linea),  exist=info)

return
end function directory_exists

!---------------------!
!---- subroutines ----!
!---------------------!

!!----
!!---- subroutine write_date_time(lun, dtim)
!!----  integer,          optional, intent(in) :: lun
!!----  character(len=*), optional, intent(out):: dtim
!!----
!!---- generic subroutine for writing the date and time
!!---- in form   date: day/month/year  time: hour:minute:second
!!---- to a file with logical unit = lun. the output argument
!!---- can be provided to get a string with the same information
!!----
!!---- updated: january - 2014
!!
subroutine write_date_time(lun, dtim)
integer,          optional, intent(in) :: lun
character(len=*), optional, intent(out):: dtim
!--- local variables ----!
character (len=10) :: dat
character (len=10) :: tim
call date_and_time(date=dat, time=tim)
if(present(lun)) &
write(unit=lun, fmt="(/, 4a)") &
" => date: ", dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4),       &
"  time: ", tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
if(present(dtim)) &
dtim="#   date: "//dat(7:8)//"/"//dat(5:6)//"/"//dat(1:4)//      &
"  time: "//tim(1:2)//":"//tim(3:4)//":"//tim(5:10)
return
end subroutine write_date_time

end module cfml_globaldeps
