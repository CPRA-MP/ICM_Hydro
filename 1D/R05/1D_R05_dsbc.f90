! Downstream boundary condition
!subroutine dsbc(n)
subroutine dsbc_R05()

    use constants_module_R05
    use arrays_module_R05
    use var_module_R05
    use sgate_module_R05

    implicit none

    ! Input
    !integer, intent(in) :: n

    ! Locals
    real(kind=4) :: ads, frds, qn, yconj

    ! compute conjugate depth at dwon stream end
    !ads=(y(n,ncomp)-z(ncomp))*bo(ncomp)
	ads=(oldY(ncomp)-z(ncomp))*bo(ncomp)
    !frds=q(n,ncomp)/sqrt(grav*ads**3.0/bo(ncomp))
	frds=oldQ(ncomp)/sqrt(grav*ads**3.0/bo(ncomp))
    !yconj=0.5*(y(n,ncomp)-z(ncomp))*(sqrt(1.0+8.0*frds*frds)-1.0)
    yconj=0.5*(oldY(ncomp)-z(ncomp))*(sqrt(1.0+8.0*frds*frds)-1.0)
    yconj=yconj+z(ncomp)
    !print *, yn,y(n,ncomp),yconj,z(ncomp)
    print *, yn,oldY(ncomp),yconj,z(ncomp)
    if(yconj < yn) then
        ! print *, 'no'
        if(option == 1) then
            ! downstream water level imposed (option 1)
            !dac(ncomp)=(yn-y(n,ncomp))*bo(ncomp)
            dac(ncomp)=(yn-oldY(ncomp))*bo(ncomp)
            !dap(ncomp)=(yn-y(n,ncomp))*bo(ncomp)
            dap(ncomp)=(yn-oldY(ncomp))*bo(ncomp)
            dqc(ncomp)=dqp(ncomp)
            ! print *, 'yes',dac(ncomp),yn,oldY(ncomp),dqc(ncomp)
        else if (option == 2) then
            ! downstream discharge imposed  (option 2)
            dqc(ncomp)=0.0
            dqp(ncomp)=0.0
            dac(ncomp)=dap(ncomp)
        elseif(option == 3) then
            ! downstream rating curve imposed (option 3)
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            qn=0.65*10*1.0*sqrt(2.0*grav*(yn-0.5))
            !dqp(ncomp)=qn-q(n,ncomp)
			dqp(ncomp)=qn-oldQ(ncomp)
            dqc(ncomp)=dqp(ncomp)
        end if
    else
        ! super critical flow exist at downstream exit
        dac(ncomp)=dap(ncomp)
        dqc(ncomp)=dqp(ncomp)
    end if

end subroutine dsbc_R05
