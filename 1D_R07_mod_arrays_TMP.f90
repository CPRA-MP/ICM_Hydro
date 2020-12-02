!-------------------------------------------------------------------------------
!                                                                               
!  Developer: Vasilia Velissariou <vvelissariou@tulane.edu>                               
!                                                                                         
!  Version: 0.2
!
!    Version - 0.1 Thu Nov 21 2019
!            - Initial development
!    Version - 0.2 Jan 08 2020
!            - Modifications array sizes to accomodate code changes
!-------------------------------------------------------------------------------
!
module mod_arrays_TMP_R07

use params_TMP_R07

contains

  subroutine aloc_arrays

!    implicit none

    allocate(depth(NDx))
    allocate(width(NDx))
    allocate(area(NDx))
    allocate(flow(NDx))
    allocate(CN(NDx, Nclass))
    allocate(CN1(NDx, Nclass))

    !allocate(C(NDtUser, NDx, Nclass))
    !allocate(CLD(NDtUser, NDx, Nclass))
    !allocate(ACLD(NDtUser, NDx, Nclass))
	
	!print*,"NDtUser, NDx, Nclass", NDtUser, NDx, Nclass
	
    !allocate(depthUser(NDtUser, NDx))
    !allocate(areaUser(NDtUser, NDx))
    !allocate(flowUser(NDtUser, NDx))
    !allocate(timeUser(NDtUser))

    allocate(xx(NDx))
    allocate(time(NDt))

    allocate(C0(NDx))

    allocate(D50(Nclass))
    allocate(D90(Nclass))
    allocate(Dgr(Nclass))

    allocate(CSS(Nclass))
    allocate(SGsed(Nclass))
    allocate(rhoSed(Nclass))
    allocate(Tcrit(Nclass))
    allocate(velSet(Nclass))
    allocate(deposition(Nclass))
    allocate(resuspension(Nclass))
    allocate(SedAccumRate(Nclass))
    allocate(CSSresusOff(Nclass))

  end subroutine aloc_arrays

end module mod_arrays_TMP_R07
