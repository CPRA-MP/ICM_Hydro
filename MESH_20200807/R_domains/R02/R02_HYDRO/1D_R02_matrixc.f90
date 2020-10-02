subroutine matrixc_R02()

    use constants_module_R02
    use arrays_module_R02
    use matrix_module_R02
    use var_module_R02

    implicit none

    ! Local variables
    integer :: i
    real(kind=4) :: e(2, 2), st(2,2), a(2,2), f_mat(2, 2), g(2, 2), d(2)
    real(kind=4) :: cour, crmax, crmin, di, dim1, dip1, dkda, ei, ei1, eia, ter

    do i=1,ncomp
        u(i)=qp(i)/areap(i)
        c(i)=sqrt(grav*areap(i)/bo(i))

        e(1,1)=1.0
        if(abs(u(i) - c(i)) < TOLERANCE) then
            c(i)=c(i)+0.00001
        end if
        e(1,2)=-1.0/(u(i)-c(i))
        e(2,1)=1.0
        e(2,2)=-1.0/(u(i)+c(i))

        f_mat(1,1)=-(u(i)-c(i))/(2.0*c(i))
        f_mat(1,2)=(u(i)+c(i))/(2.0*c(i))
        f_mat(2,1)=-(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))
        f_mat(2,2)=(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))

        d(1)=abs(u(i)+c(i))
        d(2)=abs(u(i)-c(i))

        a(1,1)=e(1,1)*f_mat(1,1)*d(1)+e(2,1)*f_mat(1,2)*d(2)
        a(1,2)=e(1,2)*f_mat(1,1)*d(1)+e(2,2)*f_mat(1,2)*d(2)
        a(2,1)=e(1,1)*f_mat(2,1)*d(1)+e(2,1)*f_mat(2,2)*d(2)
        a(2,2)=e(1,2)*f_mat(2,1)*d(1)+e(2,2)*f_mat(2,2)*d(2)

        if(ots == 1) then
           dt(i)=cfl*dx(i)/max(d(1),d(2))
        else
           dt(i)=dtini
        endif

        ! Calculating dK/dA (eq 15)
        !ter=bo(i)+2.0*areap(i)/bo(i)
        !Nazmul
        ter=pere(i)
        !Nazmul: dkda is corrected according to October 1, 2019 meeting at NWC
        dkda=sk(i)*((5.0/3.0*areap(i)**(2.0/3.0)*ter)     &
                   -(2.0/3.0*areap(i)**(5.0/3.0)*dpda(i)))/(ter**(5.0/3.0))

        ! Matrix S (eq 14)
        st(1, 1)=0.0
        st(1, 2)=0.0 !CHANGE
        !st(1, 2)=1.0
        st(2, 1)=grav*areap(i)/bo(i)/bo(i)*dbdx(i)+gso(i)      &
                +f*2.0*grav*areap(i)*qp(i)*abs(qp(i))/co(i)**3.0*dkda
        !Nazmul: st(2,2) term is multiplied with a gravity
        st(2,2)=-2*f*qp(i)*grav*areap(i)/co(i)/co(i)

        ! cour == sigma
        if(i == 1) then
            cour=dt(i)
        else if(dx(i-1) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i-1)
        endif
        b11(i)=0.5-phi-theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
        b12(i)=-theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
        b21(i)=-theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
        b22(i)=0.5-phi-theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

        if(dx(i) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i)
            crmax=max(crmax,cour*max(d(1),d(2)))
            crmin=min(crmin,cour*max(d(1),d(2)))
        endif

        !Nazmul: *dt(i) is multiplied to all g matrix terms
        g(1,1)=0.5+phi+theta*cour*a(1,1)+0.5*thes*st(1,1)*dt(i)
        g(1,2)=theta*cour*a(1,2)+0.5*thes*st(1,2)*dt(i)
        g(2,1)=theta*cour*a(2,1)+0.5*thes*st(2,1)*dt(i)
        g(2,2)=0.5+phi+theta*cour*a(2,2)+0.5*thes*st(2,2)*dt(i)

        g11inv(i)= g(2,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g12inv(i)=-g(1,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g21inv(i)=-g(2,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g22inv(i)= g(1,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))

        f1(i)=qp(i)
        f2(i)=qp(i)*qp(i)/areap(i)+grav*ci1(i)

        if(i >= 2 .and. i < ncomp) then
            dip1=areap(i+1)/bo(i+1)
            di=2*areap(i)/bo(i)
            dim1=areap(i-1)/bo(i-1)
            eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
        endif
    end do
    eps2(1)=eps2(2)
    eps2(ncomp)=eps2(ncomp-1)

    do i=2,ncomp-1
        if(ityp(i).ne.1) then
            eps2(i)=eps2(i-1)
            eps2(i+1)=eps2(i+2)
        endif
    end do

    do i=ncomp-1,1,-1
        eps2(i+1)=max(eps2(i+1),eps2(i))
        ! u(i+1)=(u(i+1)+u(i))/2.0
        ! c(i+1)=(c(i+1)+c(i))/2.0
        eps4(i+1)=max(0.,alfa4-eps2(i+1)/(u(i+1)+c(i+1)))
        ! print *, 'corr',i,eps2(i)
    end do
    d1(1)=0.0
    d2(1)=0.0
    d1(ncomp)=0.0
    d2(ncomp)=0.0

    do i=2,ncomp-1
        d(1)=abs(u(i)+c(i))
        d(2)=abs(u(i)-c(i))
        ei=max(d(1),d(2))
        d(1)=abs(u(i-1)+c(i-1))
        d(2)=abs(u(i-1)-c(i-1))
        ei1=max(d(1),d(2))
        eia=(ei+ei1)/2.0
        if(ityp(i-1) /= 1) then
            d1(i)=0.0
            d2(i)=0.0
        elseif(i == 2 .or. i == (ncomp-1)) then
            d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))
            d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))
        else
          d1(i)=eps2(i)*eia*(areap(i)-areap(i-1))-eps4(i)*(areap(i+1)  &
                            -3*areap(i)+3*areap(i-1)-areap(i-2))
          d2(i)=eps2(i)*eia*(qp(i)-qp(i-1))-eps4(i)*(qp(i+1)-3*qp(i)   &
                            +3*qp(i-1)-qp(i-2))
        endif
        ! print *, 'corr',i,d1(i),d2(i),eps2(i),area(i+1),area(i)
    end do
    !print*, 'eps2_corr', (eps2(i), i=1, ncomp)

end subroutine matrixc_R02
