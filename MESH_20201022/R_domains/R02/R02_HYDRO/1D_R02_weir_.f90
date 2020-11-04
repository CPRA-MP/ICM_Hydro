!subroutine weir(i, n)
subroutine weir_R02(i)

    use constants_module_R02
    use arrays_module_R02
    use sgate_module_R02

     !  common/var/ncomp,cfl,f,qus,yds,rhs1,rhs2,c11,c12,c21,c22,us,ots,
     ! 1           thes,vv,dtini
     !  common/matrix/phi,theta,alfa2,alfa4,w
     !  common/sgate/ag,bg,option,yn,eps,dmeu,yw,bw,ydsn


    implicit none

    ! Input
    integer, intent(in) :: i !, n

    ! Locals
    real(kind=4), parameter :: rad = 3.1415927
    real(kind=4) :: a1, a2, a3, corssign, cortsign, dd, ddsw, dusw, qq, qusw
    real(kind=4) :: rr, s, sq, ss, st, t, tht, ver, x1, x2, x3, yds, yus, yusw

    dmeu=1.0

    ! Normal (non-reverse) flow
    !yus=y(n,i-1)
    yus=oldY(i-1)
    !yds=y(n,i)
    yds=oldY(i)
    ! print *, yus,yds,yw
    if((yds-yw) > (2.*(yus-yw)/3.0)) then
        ! print *, 'yes'
        ! Flooded flow
        dusw=(area(i-1)+dap(i-1))/bo(i-1)
        yusw=dusw+z(i-1)
        !qusw=q(n,i-1)+dqp(i-1)
        qusw=oldQ(i-1)+dqp(i-1)

        ! Solving the qubic equation and finding the proper root
        ! Note: only one root will satisfy the conditions

        a1=-2.*yw-yusw
        a2=yw*yw+2.*yusw*yw
        a3=qusw*qusw/(2.*dmeu*dmeu*bw*bw*grav)-yusw*yw*yw
        qq=(3.*a2-a1*a1)/9.
        rr=(9.*a1*a2-27.*a3-2.*a1*a1*a1)/54.
        dd=qq*qq*qq+rr*rr
        x1=0.0
        x2=0.0
        x3=0.0
        corssign=1.0
        cortsign=1.0

        ! Not clear how this worked as dd is a real number, not an integer
        ! Given that the interpretation of the original goto/if statement was
        ! always unspecified it may bear looking at the original code if this is
        ! used again
        select case(int(dd))
            case(1)
                tht=acos(rr/sqrt(-qq*qq*qq))
                x1=2*sqrt(-qq)*cos(tht/3.)-a1/3.
                x2=2*sqrt(-qq)*cos(tht/3.+120.*rad/180.)-a1/3.
                x3=2*sqrt(-qq)*cos(tht/3.+240.*rad/180.)-a1/3.

            case(2)
                ver=sqrt(dd)
                if((rr+ver) < 0.0) then
                    corssign=-1.0
                end if
                if((rr-ver) < 0.0) then
                    cortsign=-1.0
                end if
                ss=sign((rr+ver),1.0)
                st=sign((rr-ver),1.0)
                s=corssign*(ss)**(1./3.)
                t=cortsign*(st)**(1./3.)
                x1=s+t-a1/3.0
                x2=-0.5*(s+t)-a1/3.0

            case(3)
                ver=sqrt(dd)
                if((rr+ver).lt.0.0) then
                    corssign=-1.0
                end if
                if((rr-ver).lt.0.0) then
                    cortsign=-1.0
                end if
                ss=sign((rr+ver),1.0)
                st=sign((rr-ver),1.0)
                s=corssign*(ss)**(1./3.)
                t=cortsign*(st)**(1./3.)
                x1=s+t-a1/3.0

        end select

        if(x1 > 0.0 .and. (x1-yw) > (2.*(yusw-yw)/3.) .and.     &
               (x1-z(i)) > ((qusw*qusw/bo(i)/bo(i)/grav)**0.333)) then
            ddsw=x1-z(i)
        endif
        if(x2 > 0.0 .and. (x2-yw) > (2.*(yusw-yw)/3.) .and.     &
              (x2-z(i)) > ((qusw*qusw/bo(i)/bo(i)/grav)**0.333)) then
            ddsw=x2-z(i)
        endif
        if(x3 > 0.0 .and. (x3-yw) > (2.*(yusw-yw)/3.) .and.     &
              (x3-z(i)) > ((qusw*qusw/bo(i)/bo(i)/grav)**0.333)) then
            ddsw=x3-z(i)
        endif
        dqp(i)=dqp(i-1)
        dap(i)=ddsw*bo(i)-area(i)
        ! print *, x1,x2,x3,yusw,qusw,dd,ddsw,z(i)
    else
        ! Free flowing flow
        dusw=(area(i-1)+dap(i-1))/bo(i-1)
        yusw=dusw+z(i-1)
        qusw=2./3.0*dmeu*bw*sqrt(2.0*grav/3.)*(yusw-yw)**1.5
        sq=qusw/bw
        ddsw=(sq*sq/grav)**0.3333333
        !dqp(i-1)=qusw-q(n,i-1)
		dqp(i-1)=qusw-oldQ(i-1)
        dqp(i)=dqp(i-1)
        dap(i)=ddsw*bo(i)-area(i)
        ! print *, i,dqp(i-1)
        !print *, 'free',dusw,qusw,ddsw,dqp(i),dap(i),y(n,i)
		print *, 'free',dusw,qusw,ddsw,dqp(i),dap(i),oldY(i)
    endif

end subroutine weir_R02
