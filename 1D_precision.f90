module precision

implicit none
!
! parameters, used in conversions: sp=single precision, hp=high (double) precision
!
integer, parameter :: ip = 4            ! precision of normal integers
integer, parameter :: rp = 4            ! precision of normal reals
integer, parameter :: sp = kind(1.0e0)  ! single precision
integer, parameter :: dp = kind(1.0d0)  ! double precision
integer, parameter :: fp = dp
integer, parameter :: hp = dp
!
! double precision integers:
!
integer, parameter :: long = SELECTED_INT_KIND(16)
!
! interfaces
!
  interface comparereal
     module procedure comparerealdouble
     module procedure comparerealsingle
  end interface
! 

contains

function comparerealdouble(val1, val2, eps)
!!--description-----------------------------------------------------------------
!
! Compares two double precision numbers
! Allow users to define the value of eps. If not, eps equals to the default machine eps
!
! Return value: -1 if val1 < val2
!                0 if val1 = val2
!               +1 if val1 > val2
!
!!--pseudo code and references--------------------------------------------------
!
! The functionality in this subroutine is copied from subroutine Ifdbl,
! written by Jaap Zeekant.
!
! eps must be machine precision dependent.
! eps may not be given by the user! See what happens when
! val1 = -666.0, val2 = -999.0, eps = 0.5
!
!!--declarations----------------------------------------------------------------
    implicit none
!
! Return value
!
integer :: comparerealdouble
!
! Global variables
!
real(hp), intent(in)           :: val1
real(hp), intent(in)           :: val2
real(hp), optional, intent(in) :: eps
!
! Local variables
!
real(hp) :: eps0
real(hp) :: value
!
!! executable statements -------------------------------------------------------
!
if (present(eps)) then
    eps0 = eps
else 
    eps0 = 2.0_hp * epsilon(val1)
endif
!
if (abs(val1)<1.0_hp .or. abs(val2)<1.0_hp) then
   value = val1 - val2
else
   value = val1/val2 - 1.0_hp
endif
!
if (abs(value)<eps0) then
   comparerealdouble = 0
elseif (val1<val2) then
   comparerealdouble = -1
else
   comparerealdouble = 1
endif
end function comparerealdouble


function comparerealsingle(val1, val2,eps)
!!--description-----------------------------------------------------------------
!
! REMARK: THE NAME OF THIS FUNCTION IS WRONG!
!         The name should be comparefp
!
! Compares two real numbers of type fp
! Allow users to define the value of eps. If not, eps equals to the default machine eps
!
! Return value: -1 if val1 < val2
!                0 if val1 = val2
!               +1 if val1 > val2
!
!!--pseudo code and references--------------------------------------------------
!
! The functionality in this subroutine is copied from subroutine Ifflt,
! written by Jaap Zeekant.
!
! eps must be machine precision dependent.
! eps may not be given by the user! See what happens when
! val1 = -666.0, val2 = -999.0, eps = 0.5
!
!!--declarations----------------------------------------------------------------
implicit none
!
! Return value
!
integer :: comparerealsingle
!
! Global variables
!
real(sp), intent(in)           :: val1
real(sp), intent(in)           :: val2
real(sp), optional, intent(in) :: eps
!
! Local variables
!
real(sp) :: eps0
real(sp) :: value
!
!! executable statements -------------------------------------------------------
!
!  
if (present(eps)) then
    eps0 = eps
else
    eps0 = 2.0_sp * epsilon(val1)
endif
!
if (abs(val1)<1.0_sp .or. abs(val2)<1.0_sp) then
   value = val1 - val2
else
   value = val1/val2 - 1.0_sp
endif
!
if (abs(value)<eps0) then
   comparerealsingle = 0
elseif (val1<val2) then
   comparerealsingle = -1
else
   comparerealsingle = 1
endif
end function comparerealsingle


end module precision
