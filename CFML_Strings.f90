!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2022  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----          Nebil Ayape Katcho      (ILL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----               Ross Angel         (University of Pavia)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_String_Utilities
!!----   INFO: Manipulation of strings with alfanumeric characters
!!----
!!----
!!
 Module CFML_Strings
    !---- Use Modules ----!
    use Ieee_Arithmetic,   only: ieee_is_nan,ieee_is_finite
    use CFML_GlobalDeps,   only: cp, ops_sep, err_cfml, clear_error, set_error
    use CFML_Maths,        only: Negligible, Zbelong

    implicit none

    private

    !---- List of public functions ----!
    public :: Equal_Sets_Text, Frac_Trans_1Dig, Frac_Trans_2Dig,          &
              Get_DateTime, Get_Dirname, Get_Extension, Get_Filename,     &
              Get_Mat_From_Symb, Get_Vec_From_String,L_Case, U_Case,      &
              NumCol_from_NumFmt, Pack_String, Read_Fract,Number_Lines,   &
              Set_Symb_From_Mat, String_Count, Strip_String, String_Real, &
              String_Fraction_1Dig, String_Fraction_2Dig, String_NumStd,  &
              Reading_File, File_To_FileList, Get_Vec_from_FracStr,       &
              Num_Items


    !---- List of public subroutines ----!
    public :: Cut_string, FindFMT, &
              Get_Separator_Pos, Get_Substring_Positions, Get_Words,      &
              Get_Num, Get_NumStd, Get_Transf, Init_FindFmt, Inc_LineNum, &
              Reading_Lines, Read_Key_Str,Read_Key_StrVal,Read_Key_Value, &
              Read_Key_ValueSTD, Sort_Strings, SubString_Replace


    !!----
    !!---- TYPE :: FILE_LIST_TYPE
    !!--..
    !!---- Type,public :: File_List_Type
    !!----    integer                                       :: nlines ! Number of lines in the file
    !!----    character(len=256), allocatable, dimension(:) :: line   ! Content of the lines
    !!---- End Type file_type
    !!----
    !!---- Updated: February - 2005, November 2012, February 2020 (moved from CFML_IO_FORM)
    !!
    Type, public :: File_List_Type
       integer                                       :: nlines=0   ! Number of lines
       character(len=256), dimension(:), allocatable :: line       ! Strings containing the lines of the file
    End Type File_List_Type


    Type, public :: String_Array_Type          !Type for handling allocatable arrays of allocatable strings
      character(len=:), allocatable :: str
    End Type String_Array_Type

    !!----
    !!---- TYPE :: FILE_TYPE
    !!--..
    !!---- Type,public :: File_Type
    !!----    character(len=:),   allocatable                    :: Fname  ! Original name of the file
    !!----    integer                                            :: nlines ! Number of lines in the file
    !!----    Type(String_Array_Type), dimension(:), allocatable :: line     ! Content of the lines
    !!---- End Type file_type
    !!----
    !!---- Updated: February - 2005, November 2012, February 2020 (moved from CFML_IO_FORM)
    !!
    Type, public :: File_Type
       character(len=:),               allocatable        :: Fname      ! Name of file
       integer                                            :: nlines=0   ! Number of lines
       Type(String_Array_Type), dimension(:), allocatable :: line       ! Strings containing the lines of the file
    End Type File_Type



    !--------------------!
    !---- Parameters ----!
    !--------------------!
    character(len=*), parameter :: DIGIT         ="0123456789.-"     ! Character parameter for numbers
    character(len=*), parameter :: DIGIT_EXT     ="0123456789.-()"   ! Character parameter for numbers and parenthesis

    Interface
       Module Subroutine BuildFMT(iFld,nCar,nStr,FMTstring)
          !---- Arguments ----!
          Integer,           intent(in    ) ::   iFld       ! Format type
          Integer,           intent(in out) ::   nCar       ! integer/real field: number of characters in field
                                                            ! character field: number of characters to skip before A field
          Integer,           intent(in out) ::   nStr       ! current character number in FMTstring
          Character (len=*) ,intent(in out) ::   FMTstring  ! FORTRAN format string
       End Subroutine BuildFMT

       Pure Module Subroutine Cut_String(Str1,nlong1,Str2,nlong2)
          !---- Argument ----!
          character(len=*),           intent(in out) :: Str1     ! Input string / Out: string without the first word
          character(len=*), optional, intent(   out) :: Str2     ! The first word of String on Input
          integer,          optional, intent(   out) :: nlong1   ! Give the length of Str1 on Output
          integer,          optional, intent(   out) :: nlong2   ! Give the length of Str2 on Output
       End Subroutine Cut_String

       Pure Module Function Equal_Sets_Text(str1,n1,str2,n2) result(Equal)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in) :: str1   ! Vector of String
          character(len=*), dimension(:), intent(in) :: str2   ! Vector of String
          integer,                        intent(in) :: n1     ! Number of lines on Text1
          integer,                        intent(in) :: n2     ! Number of lines on str2
          logical                                    :: Equal  !
       End Function Equal_Sets_Text

       Module Subroutine FindFmt(IUnit,aLine,FMTfields,FMTstring,idebug)
          !---- Arguments ----!
          Integer ,           intent(in    ) ::  IUnit      ! Logical unit number
          Character (len=*) , intent(in out) ::  aLine      ! character string to be decoded
          Character (len=*) , intent(in    ) ::  FMTfields  ! description of the format fields (e.g. IIFIF)
          Character (len=*) , intent(   out) ::  FMTstring  ! format of the line (e.g. (I5,I1,F8.0,I4,F7.0,)
          Integer ,optional,  intent(in    ) ::  idebug     ! Logical unit number for writing the input file
       End Subroutine FindFmt

       Module Subroutine FindFMT_Err(aLine,nC_L)
          !---- Arguments ----!
          Character(len=*), intent(in) ::   aLine      ! Current data line
          Integer,         intent (in) ::   nC_L       ! location of last character treated
       End Subroutine FindFMT_Err

       Pure Module Function Frac_Trans_1Dig(Vec) Result(Str)
          !---- Argument ----!
          real(kind=cp), dimension(3), intent( in)   :: Vec  ! Vector
          character(:),allocatable                   :: Str  ! String with conversion to fractional
       End Function Frac_Trans_1Dig

       Pure Module Function Frac_Trans_2Dig(Vec) Result(Str)
          !---- Argument ----!
          real(kind=cp), dimension(3), intent(in) :: Vec   ! Vector
          character(:), allocatable               :: Str   ! String with conversion to fractional
       End Function Frac_Trans_2Dig

       Module Function Get_DateTime() Result(Str)
          !---- Argument ----!
          character(len=:), allocatable :: Str  ! String containing the Date and Time
       End Function Get_DateTime

       Pure Module Function Get_Dirname(Str) Result(Directory)
          !---- Argument ----!
          character(Len=*), Intent (In)  :: Str          ! String containing Path + Filename
          character(Len=:), allocatable  :: Directory    ! Path
       End Function Get_Dirname

       Pure Module Function Get_Extension(filename, dotted) Result(extension)
          !---- Arguments ----!
          character(len=*),  intent(in)  :: filename     ! Input filename
          logical, optional, intent(in)  :: dotted       ! If True, the extension will be returned with a dot
          character(len=:), allocatable  :: extension    ! Extension of the file
       End Function Get_Extension

       Pure Module Function Get_Filename(Str) Result(Filename)
          !---- Argument ----!
          character(Len=*), intent(in)  :: Str       ! String containing Path + Filename
          character(Len=:), allocatable :: Filename  ! Filename
       End Function Get_Filename

       Module Function Get_Mat_From_Symb(Symb,cod) Result(Mat)
          !---- Arguments ----!
          character(len=*),                intent(in)  :: Symb   ! String
          character(len=1), dimension(3),  intent(in)  :: cod    ! (/"u","v","w"/) or (/"x","y","z"/)
          real(kind=cp),dimension(3,3)                 :: Mat    ! Output
       End Function Get_Mat_From_Symb

       Module Subroutine Get_Num(Str,vet,ivet,iv)
          !---- Argument ----!
          character (len=*),          intent ( in) :: Str   ! Input String to convert
          real(kind=cp), dimension(:),intent (out) :: vet   ! Vector of real numbers
          integer, dimension(:),      intent (out) :: ivet  ! Vector of integer numbers
          integer,                    intent (out) :: iv    ! Number of numbers in Vet/Ivet
       End Subroutine Get_Num

       Module Subroutine Get_NumStd(Str, value, std, ic)
          !----Arguments ----!
          character(len=*),             intent( in) :: Str     ! Input String
          real(kind=cp), dimension(:),  intent(out) :: value   ! Vector of values with real numbers
          real(kind=cp), dimension(:),  intent(out) :: std     ! Vector of standard deviation values
          integer,                      intent(out) :: ic      ! Number of components of vector Value
       End Subroutine Get_NumStd

       Module Subroutine Get_Transf(str,mat,v,cod)
          !---- Arguments ----!
          character(len=*),                          intent(in)  :: str      ! Input string
          real(kind=cp),dimension(3,3),              intent(out) :: mat      ! Matrix
          real(kind=cp),dimension(3),     optional,  intent(out) :: v        ! Vector
          character(len=1), dimension(4), optional,  intent(in)  :: cod      ! Code
       End Subroutine Get_Transf

       Module Function Get_Vec_from_FracStr(Str) Result(V)
          !---- Arguments ----!
          character(len=*), intent(in) :: str
          real(kind=cp), dimension(3)  :: V
       End Function Get_Vec_from_FracStr

       Pure Module Subroutine Get_Separator_Pos(Str,car,pos,ncar)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Str   ! Inout String
          character(len=1),      intent(in)  :: car   ! Separator character
          integer, dimension(:), intent(out) :: pos   ! Vector with positions of "sep" in "Line"
          integer,               intent(out) :: ncar  ! Number of appearance of "sep" in "Line"
       End Subroutine Get_Separator_Pos

       Pure Module Subroutine Get_Substring_Positions(str,substr,pos,nsubs)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: str         ! In -> Input String
          character(len=*),      intent(in)  :: substr      ! In -> Substring
          integer, dimension(:), intent(out) :: pos         ! Out -> Vector with positions of the firs character of "substr" in "String"
          integer,               intent(out) :: nsubs       ! Out -> Number of appearance of "substr" in "String"
       End Subroutine Get_Substring_Positions

       Module Function Get_Vec_From_String(Str,Cod) Result(Vec)
          !---- Arguments ----!
          character(len=*),                intent(in)  :: str   ! Input string
          character(len=1), dimension(3),  intent(in)  :: cod   ! Code
          real(kind=cp),dimension(3)                   :: vec   ! Vector
       End Function Get_Vec_From_String

       Module Subroutine Get_Words(Str,dire,ic,sep)
          !---- Argument ----!
          character(len=*),                 intent ( in) :: Str   ! Input string
          character(len=*), dimension(:),   intent (out) :: dire  ! Vector of Words
          integer,                          intent (out) :: ic    ! Number of words
          character(len=*), optional,       intent ( in) :: sep   ! separator other than blank
       End Subroutine Get_Words

       Module Subroutine Inc_LineNum(line_n)
          !---- Argument ----!
          integer, intent(in) :: line_n
       End Subroutine Inc_LineNum

       Module Subroutine Init_FindFMT(nline)
          !---- Arguments ----!
          integer, optional, intent(in) :: nline
       End Subroutine Init_FindFMT

       Pure Module Function L_Case(Str) Result (LStr)
          !---- Argument ----!
          character (len=*), intent(in)   :: Str    ! Input String
          character (len=:), allocatable  :: LStr   ! lower case of Text
       End Function L_Case

       Pure Module Function NumCol_from_NumFmt(Str) Result(n_col)
          !---- Argument ----!
          character (len=*), intent(in)  :: Str    ! Input format string
          integer                        :: n_col  ! Integer number of columns
       End Function NumCol_from_NumFmt

       Pure Module Function Num_Items(string,separator) result(nitems)
         character(len=*), intent (in)           :: string
         character(len=1), intent (in), optional :: separator
         integer                                 :: nitems
       End Function Num_Items

       Pure Module Function Pack_String(Str) Result (Strp)
          !---- Argument ----!
          character(len=*), intent(in) :: str    ! Input String
          character(len=len_trim(str)) :: strp   ! Output string
       End Function Pack_String

       Module Function Read_Fract(str) Result(value)
          !---- Arguments ----!
          character(len=*), intent(in) :: str     ! Input String
          real(kind=cp)                :: value   ! Value
       End Function Read_Fract

       Module Subroutine Read_Key_Str(filevar,nline_ini,nline_end,keyword,string,comment)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)      :: filevar     ! Input vector of String
          integer,                        intent(in out)  :: nline_ini   ! Pointer to initial position to search
                                                                         ! Out -> Pointer to final position in search
          integer,                        intent(in)      :: nline_end   ! Pointer to final position to search
          character(len=*),               intent(in)      :: keyword     ! Word to search
          character(len=*),               intent(out)     :: string      ! Rest of the input string
          character(len=1), optional,     intent(in)      :: comment     ! Character that define a comment line
       End Subroutine Read_Key_Str

       Module Subroutine Read_Key_StrVal(filevar,nline_ini,nline_end,keyword,string,vet,ivet,iv,comment)
          !---- Arguments ----!
          character(len=*), dimension(:),           intent(in)      :: filevar       !  In -> Input vector of String
          integer,                                  intent(in out)  :: nline_ini     !  In -> Pointer to initial position to search
                                                                                     ! Out -> Pointer to final position in search
          integer,                                  intent(in)      :: nline_end     !  In -> Pointer to final position to search
          character(len=*),                         intent(in)      :: keyword       !  In -> Word to search
          character(len=*),                         intent(out)     :: string        ! Out -> Rest of the input string
          real(kind=cp),dimension(:),     optional, intent(out)     :: vet           ! Out -> Vector for real numbers
          integer,dimension(:),           optional, intent(out)     :: ivet          ! Out -> Vector for integer numbers
          integer,                        optional, intent(out)     :: iv            ! Out -> Number of numbers
          character(len=1),               optional, intent(in)      :: comment       ! Character that define a comment line
       End Subroutine Read_Key_StrVal

       Module Subroutine Read_Key_Value(filevar,nline_ini,nline_end,keyword,vet,ivet,iv,comment,line_key)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)     :: filevar            !  In -> Input vector of String
          integer,                        intent(in out) :: nline_ini          !  In -> Pointer to initial position to search
                                                                               ! Out -> Pointer to final position in search
          integer,                        intent(in)     :: nline_end          !  In -> Pointer to final position to search
          character(len=*),               intent(in)     :: keyword            !  In -> Word to search
          real(kind=cp),dimension(:),     intent(out)    :: vet                ! Out -> Vector for real numbers
          integer,dimension(:),           intent(out)    :: ivet               ! Out -> Vector for integer numbers
          integer,                        intent(out)    :: iv                 ! Out -> Number of components
          character(len=1),     optional, intent(in)     :: comment            ! Consider the character passed in comment as a comment to skip the line
          character(len=*),     optional, intent(out)    :: line_key           ! Out -> Cut line where keyword is read
       End Subroutine Read_Key_Value

       Module Subroutine Read_Key_ValueSTD(filevar,nline_ini,nline_end,keyword,vet1,vet2,iv,comment)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: filevar         !  In -> Input vector of String
          integer,                         intent(in out) :: nline_ini       !  In -> Pointer to initial position to search
                                                                             ! Out -> Pointer to final position in search
          integer,                         intent(in)     :: nline_end       !  In -> Pointer to final position to search
          character(len=*),                intent(in)     :: keyword         !  In -> Word to search
          real(kind=cp),dimension(:),      intent(out)    :: vet1            ! Out -> Vector of real numbers
          real(kind=cp),dimension(:),      intent(out)    :: vet2            ! Out -> Vector of standard deviations
          integer,                         intent(out)    :: iv              ! Out -> Number of components
          character(len=1),      optional, intent(in)     :: comment         ! Consider the character passed in comment as a comment to skip the line
       End Subroutine Read_Key_ValueSTD

       Pure Module Function Set_Symb_From_Mat(Mat,cod) Result(Symb)
          !---- Arguments ----!
          real(kind=cp),dimension(3,3),    intent(in)  :: Mat    ! Array
          character(len=1), dimension(3),  intent(in)  :: cod    ! Codes (/"u","v","w"/) or (/"x","y","z"/)
          character(len=:), allocatable                :: Symb   ! Symbol
       End Function Set_Symb_From_Mat

       Module Subroutine SGetFTMfield(GetFTMfield,FMTfields,nFld,nFldMax)
          !---- Arguments ----!
          Integer ,          intent(out)    ::  GetFTMfield
          Character (len=*) ,intent( in)    ::  FMTfields        !  -> format descriptor
          Integer ,          intent(in out) ::  nFld             ! <-> current field in format descriptor
          Integer ,          intent( in)    ::  nFldMax          !  -> max. number of fields in format descriptor
       End Subroutine SGetFTMfield

       Pure Module Subroutine Sort_PR_Partition(A, Marker)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in out) :: A
          integer,                        intent(   out) :: marker
       End Subroutine Sort_PR_Partition

       Recursive Module Subroutine Sort_Strings(Str)
          !---- Argument ----!
          character(len=*), dimension(:), intent(in out) :: Str
       End Subroutine Sort_Strings

       Pure Module Function String_Count(str,substr) Result(N)
          !---- Arguments ----!
          character(len=*), intent(in) :: str       ! Input String
          character(len=*), intent(in) :: substr    ! Substring model
          integer                      :: N         ! Number
       End Function String_Count

       Pure Module Function String_Fraction_1Dig(V) Result(Str)
          !---- Argument ----!
          real(kind=cp),    intent( in) :: V   !  Real value
          character(:), allocatable     :: Str !  Fracction in character form
       End Function String_Fraction_1Dig

       Pure Module Function String_Fraction_2Dig(V) Result(Str)
          !---- Argument ----!
          real(kind=cp),    intent( in) :: v    ! Real value
          character(:), allocatable     :: Str  ! Fraction in character form
       End Function String_Fraction_2Dig

       Pure Module Function String_NumStd(Value, Std) Result(Str)
          !---- Argument ----!
          real(kind=cp),   intent(in)  :: Value    ! Value
          real(kind=cp),   intent(in)  :: Std      ! Standard deviation
          character(len=:),allocatable :: Str      ! String containing the information
       End Function String_NumStd

       Pure Module Function String_Real(Val,W) Result(Str)
          !---- Arguments ----!
          real(kind=cp), intent(in)  :: val        ! value to be output
          integer,       intent(in)  :: w          ! Width
          character(len=w)           :: Str
       End Function String_Real

       Pure Module Function Strip_String(str, to_strip) Result(sstr)
          !---- Arguments----!
          character(len=*), intent(in)  :: str            ! Input string
          character(len=*), intent(in)  :: to_strip       ! Pattern
          character(len=len_trim(str))  :: sstr
       End Function Strip_String

       Pure Module Subroutine SubString_Replace(string, substr, repstr, warning)
          !---- Arguments ----!
          character(len=*), intent(in out) :: string   ! Input/output string
          character(len=*), intent(in)     :: substr   ! Subtring to be replaced
          character(len=*), intent(in)     :: repstr   ! String for add
          character(len=*), intent(out)    :: warning  ! Message
       End Subroutine SubString_Replace

       Module Subroutine TreatMCharField(iFld,aLine,L_Line,nC_L,nC_X)
          !---- Arguments ----!
          Integer,           intent(in out)  :: iFld      ! <-> "A" format size (1 to 9)
          Character (len=*), intent(in)      :: aLine     !  -> data line to be analysed
          Integer,           intent(in)      :: L_Line    !  -> true length of data Line
          Integer,           intent(in out)  :: nC_L      ! <-> current character in data line
          Integer,           intent(out)     :: nC_X      ! <-  number of characters in X format field (now nx -> trn)
       End Subroutine TreatMCharField

       Module Subroutine TreatNumerField(iFld,aLine,L_Line,nC_L,nCar)
          !---- Arguments ----!
          Integer ,          intent( in)    ::  iFld   ! field type
          Character (len=*), intent(in out) ::  aLine  ! data line
          Integer ,          intent( in)    ::  L_Line ! true length of the data line
          Integer ,          intent(in out) ::  nC_L   ! counts characters in data line
          Integer ,          intent(in out) ::  nCar   ! counts characters in format field
       End Subroutine TreatNumerField

       Pure Module Function U_Case(Str) Result (UStr)
          !---- Argument ----!
          character(len=*), intent(in)  :: Str   ! Input string
          character(len=:), allocatable :: UStr  ! Upper conversion
       End Function U_Case

    End Interface


 Contains
    !!----
    !!---- READING_File
    !!----    Function Reading_File(filename) result (filecont)
    !!----    character(len=*), intent( in) :: filename    ! Filename
    !!----    type(File_Type)               :: filecont    ! File_Type variable containing the lines
    !!----
    !!----
    !!----    Read the file and put the information on the File_Type object Filecont.
    !!----    This function is similar to subroutine Reading_Lines, except that it constructs
    !!----    the File_Type object Filecont. The file is opened to read the lines and closed before
    !!----    returning to the calling unit.
    !!----
    !!---- 24/02/2020
    !!
    Function Reading_File(filename) result (filecont)
       !---- Arguments ----!
       character(len=*), intent( in) :: filename    ! Filename
       type(File_Type)               :: filecont    ! File_Type object containing the lines

       !---- Local Variables ----!
       logical            :: info,opn
       integer            :: lun,i,olun,nlines,ier
       character(len=256) :: buffer

       !> Init
       info=.false.
       filecont%fname=trim(filename)
       filecont%nlines=0

       !> Exist filename ?
       inquire (file=trim(filename),exist=info)
       if (.not. info) then
          err_cfml%ierr=1
          Err_CFML%flag=.true.
          err_cfml%msg="The file: "//trim(filename)//" does not exist "
          return
       end if

       !> Is open this file?
       inquire(file=trim(filename),opened=opn, number=olun)   !Check if the file is already opened
       if (opn) then
          rewind(olun)
          lun=olun
       else
          open(newunit=lun,file=filename, status="old",action="read", position="rewind")
       end if

       !Reading the number of lines
       nlines=0
       do
         read(unit=lun,fmt="(a)",iostat=ier) buffer
         if(ier /= 0) Exit
         nlines=nlines+1
       end do
       if(nlines == 0) then
          err_cfml%ierr=1
          Err_CFML%flag=.true.
          err_cfml%msg="The file: "//trim(filename)//" contains no lines ! "
          return
       end if
       rewind(unit=lun)
       filecont%nlines=nlines
       allocate(filecont%line(nlines))
       do i=1,nlines
          read(unit=lun,fmt="(a)",iostat=ier) buffer
          filecont%line(i)%str=trim(buffer)
       end do
       if (.not. opn) close(unit=lun)

    End Function Reading_File

    !!----
    !!---- READING_LINES
    !!----    Read nlines of the file and put the information on Filevar.
    !!----    The file is opened to read the lines and closed before
    !!----    returning to the calling unit.
    !!----
    !!---- 05/04/2019
    !!
    Subroutine Reading_Lines(filename,nlines,filevar)
       !---- Arguments ----!
       character(len=*),               intent( in) :: filename    ! Filename
       integer,                        intent( in) :: nlines      ! Number of lines to be readen
       character(len=*), dimension(:), intent(out) :: filevar     ! String vector containing the lines

       !---- Local Variables ----!
       logical :: info,opn
       integer :: lun,i,olun

       !> Init
       info=.false.

       !> Exist filename ?
       inquire (file=trim(filename),exist=info)
       if (.not. info) then
          err_cfml%ierr=1
          Err_CFML%flag=.true.
          err_cfml%msg="The file: "//trim(filename)//" does not exist "
          return
       end if

       !> Is open this file?
       inquire(file=trim(filename),opened=opn, number=olun)   !Check if the file is already opened
       if (opn) then
          rewind(olun)
          lun=olun
       else
          open(newunit=lun,file=filename, status="old",action="read", position="rewind")
       end if

       !> Reading...
       do i=1,nlines
          read(unit=lun,fmt="(a)") filevar(i)
       end do

       if (.not. opn) close(unit=lun)

    End Subroutine Reading_Lines

    !!----
    !!---- NUMBER_LINES
    !!----    Return the number of lines contained in a file. The file will be opened and closed before
    !!----    returning to the calling unit. Or in the case the file is already opened the final
    !!----    status is that the pointer for reading is put at the "rewind" (first line) position.
    !!----    If 'input_string' is present, return the number of lines until 'input_string' is founded
    !!----    as first string in the line
    !!----    (example : input_string =='END' : avoid Q peaks in a SHELX file)
    !!----
    !!---- 05/04/2019
    !!
    Function Number_Lines(filename, cond_string) Result(N)
       !---- Arguments ----!
       character(len=*),           intent(in)  :: filename       ! Filename
       character(len=*), optional, intent(in)  :: cond_string    ! String to exit
       integer                                 :: n              ! Number of lines in the file

       !---- Local Variables ----!
       logical            :: info,opn
       integer            :: lun,cond,olun
       character (len=256):: read_line                             ! TR may 2013
       integer            :: lon                                  ! TR may 2013

       !> Init
       n=0

       info=.false.
       cond=0
       if (present(cond_string)) lon=len_trim(cond_string)    ! TR may 2013

       !> Exist filename ?
       inquire (file=trim(filename),exist=info)
       if (.not. info) then
          err_cfml%ierr=1
          Err_CFML%flag=.true.
          err_cfml%msg="Number_lines@STRINGS: The file: "//trim(filename)//" does not exist "
          return
       end if

       !> Is open
       inquire(file=trim(filename),opened=opn, number=olun)   !Check if the file is already opened
       if(opn) then
          rewind(olun)
          lun=olun
       else
          open(newunit=lun,file=trim(filename), status="old",action="read", position="rewind")
       end if

       !> Counting lines
       do
          read(unit=lun,fmt="(a)",iostat=cond) read_line
          if (cond /= 0) exit
          read_line=adjustl(read_line)
          if (present(cond_string)) then                                         ! TR may 2013
             if (u_case(read_line(1:lon)) == u_case(cond_string(1:lon))) exit
          end if
          n=n+1
       end do

       if (.not. opn) then
          close(unit=lun)
       else
          rewind(unit=lun)
       end if

    End Function Number_Lines

    !!----
    !!---- Function File_To_FileList(File_dat,File_list)
    !!----   character(len=*),     intent( in) :: file_dat  !Input data file
    !!----   type(file_list_type), intent(out) :: file_list !File list structure
    !!----
    !!----    Charge an external file to an object of File_List_Type.
    !!----
    !!---- Update: March - 2023
    !!
    Function File_To_FileList(File_dat) result(File_list)
       !---- Arguments ----!
       character(len=*),      intent( in) :: file_dat
       type(file_list_type)               :: file_list

       !---- Local Variables ----!
       integer                           :: nlines

       !---- Number of Lines in the input file ----!
       nlines=Number_Lines(trim(File_dat))

       if (nlines == 0) then
          err_cfml%ierr=1
          err_CFML%Flag=.true.
          err_CFML%Msg="The file "//trim(File_dat)//" contains nothing"
          return
       else
          file_list%nlines=nlines
          if (allocated(file_list%line)) deallocate(file_list%line)
          allocate(file_list%line(nlines))
          call reading_Lines(trim(File_dat),nlines,file_list%line)
       end if

    End Function File_To_FileList

 End Module CFML_Strings
