!subroutine weirp(i, n)
subroutine weirp_R08(i)

    use constants_module_R08
    use arrays_module_R08
    use sgate_module_R08

     !  common/var/ncomp,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,ots,
     ! 1           thes,vv,dtini
     !  common/matrix/phi,theta,alfa2,alfa4,w

    implicit none

    ! Input
    integer, intent(in) :: i !, n

    ! Locals
    real(kind=4) :: yus, ddsw, dusw, frds, qdsw, yds, ydsw, yusw

    yus=areap(i)/bo(i)+z(i)
    ! yds=(areap(i+1))/bo(i+1)+z(i+1)
    ydsn=(area(i+1)+dac(i+1))**3./bo(i+1)
    !frds=(q(n,i+1)+dqc(i+1))/sqrt(grav*ydsn)
	frds=(oldQ(i+1)+dqc(i+1))/sqrt(grav*ydsn)
    print *, 'weirp',yus,yds
    if(frds < 1) then
        ! if((yds-yw).gt.(2.*(yus-yw)/3.0)) then

        ! Flooded flow
        ddsw=(area(i+1)+dac(i+1))/bo(i+1)
        ydsw=ddsw+z(i+1)
        !qdsw=q(n,i+1)+dqc(i+1)
		qdsw=oldQ(i+1)+dqc(i+1)
        yusw=ydsw+(qdsw/(dmeu*bw*(ydsw-yw)))**2.0/(2.0*grav)
        dusw=yusw-z(i)
        dqc(i)=dqc(i+1)
        dac(i)=dusw*bo(i)-area(i)
        ! print *, 'flooded-corr',ddsw,qdsw,dusw
        ! print *, 'flooded dqp-corr dap',dqc(i),dac(i)
    else

        ! Free flowing flow
        dac(i)=dap(i)
        dqc(i)=dqp(i)
        dqc(i+1)=dqc(i)
        dac(i+1)=dap(i+1)
        ! print *, i,dqc(i+1)
        print *, 'free',dac(i),dqc(i)
    endif

end subroutine weirp_R08
