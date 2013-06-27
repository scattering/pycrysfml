      module Struct
      implicit none
      
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: qp = selected_real_kind(33, 4931)

      type :: vector
        integer :: x
        integer :: y
        integer :: z
      end type vector
      
      type :: vectorN
        integer, dimension(5) :: v
      end type vectorN
      
      Type :: Sym_Oper_Type
        integer, dimension(3,3) :: Rot
        real(kind=sp), dimension(3) :: Tr
      End Type Sym_Oper_Type
      
      type(Sym_Oper_Type) :: symop

      integer, parameter :: cp = selected_real_kind(6,30)
    
      Type, public :: Reflect_Type
          integer,dimension(3) :: H     ! H
          integer              :: Mult  ! mutiplicity
          real(kind=cp)        :: S     ! Sin(Theta)/lambda=1/2d
      End Type Reflect_Type

      Type :: reflect_list_type
          integer :: num
          type(Reflect_type), allocatable, dimension(:) :: refs
      End Type reflect_list_type

      Type :: thing
          integer :: length
          real(kind=cp),allocatable,dimension(:) :: array
      End Type thing
    
      contains

      subroutine makeThing(x)
          type(Thing) :: x
          
          integer :: length, i
          length = 10
          
          if(allocated(x%array)) then
            deallocate(x%array)
            write(*,*) "deallocated"
          else
            write(*,*) "not deallocated"
          end if
          allocate(x%array(length))
          x%length = length
   
          do i=1,length
             x%array(i) = 1.0/i
          end do
      end subroutine makeThing

      subroutine makeRefList(reflist)
          type(reflect_list_type) :: reflist
          
          integer :: length, i
          length = 6
          
          if(allocated(reflist%refs)) then
            deallocate(reflist%refs)
            write(*,*) "deallocated"
          end if
          allocate(reflist%refs(length))
          reflist%num = length
   
          do i=1,length
             reflist%refs(i)%h = (/i, i, i/)
             reflist%refs(i)%mult = 1
             reflist%refs(i)%s = i
          end do
      end subroutine makeRefList

      subroutine Hkl_Uni_Reflect(Friedel,Value1,Value2,Code,Num_Ref,Reflex,no_order)
          !---- Arguments ----!
          Logical,                              intent(in)     :: Friedel
          real(kind=cp),                        intent(in)     :: value1,value2
          character(len=1),                     intent(in)     :: code
          integer,                              intent(out)    :: num_ref
          type (Reflect_Type),    dimension(:), intent(out)    :: reflex
          logical,                   optional,  intent(in)     :: no_order
    
          Num_Ref = 1
          reflex(1)%Mult = 3
      end subroutine hkl_uni_reflect
      
      function make()
          type(Sym_Oper_Type) :: make
          symop%Rot(1,1) = 42
          make = symop
          symop%Rot(1,1) = 43
      end function make
      
      function made(sym)
          type(Sym_Oper_Type) :: sym
          integer :: made
          made = sym%Rot(1,1)
      end function made
      
      function add(u,v)
          type(vector) :: add
          type(vector) :: u,v
          !print *, "adding",u%x,u%y,u%z,"to",v%x,v%y,v%z
          add = vector(u%x + v%x, u%y + v%y, u%z + v%z)
      end function add
      
      function structret()
          type(vectorN) :: structret
          structret = vectorN(0)
          structret%v(1) = 10
      end function structret
      
      function triple(a)
          integer :: triple
          integer :: a
          triple = 3*a
      end function triple
      
      function norm(v)
          integer :: norm
          type(vector) :: v
          norm = (v%x)**2 + (v%y)**2 + (v%z)**2
      end function norm

      function length(s,n)
          character(len=*) :: s
          integer :: n
          integer :: length
          write(*,*) s
          length = len(s)*n
      end function length

      function length2(s,n)
          character(len=*) :: s
          integer, optional :: n
          integer :: length2
          write(*,*) s
          if(present(n)) then
              length2 = len(s)*n
              write(*,*) "present"
          else
              length2 = len(s)
              write(*,*) "absent"
          end if
      end function length2

      subroutine unitVector(v)
          type(vector) :: v
          v = vector(1, 0, 0)
      end subroutine unitVector

      function charIndex(c)
          character(len=1) :: c
          integer :: charIndex
          charIndex = iachar(c)
      end function charIndex

      function arraySum(a)
          integer, dimension(:) :: a
          integer :: arraySum
          arraySum = sum(a)
      end function arraySum

      function logicTest(b)
          logical :: b
          logical :: logicTest
          logicTest = .not. b
      end function logicTest

      end module
