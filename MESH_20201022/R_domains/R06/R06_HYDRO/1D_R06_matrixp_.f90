subroutine matrixp_R06()

    use constants_module_R06
    use arrays_module_R06
    use matrix_module_R06
    use var_module_R06
    use arrays_section_module_R06

    implicit none

    ! Local
    integer :: i
    real(kind=4) :: e(2, 2), f_mat(2, 2), a(2, 2), d(2), st(2, 2), g(2, 2)
    real(kind=4) :: cour, crmax, crmin, di, dim1, dip1, dkda, ei, ei1, eia, ter

    do i=1,ncomp
        u(i)=oldQ(i)/area(i)

        ! Nazmul: TODO, check if this equation is good for natural channel
        c(i)=sqrt(grav*area(i)/bo(i))

        !print*, sqrt(grav*area(i)/bo(i)), sqrt(grav*(depth(i)))

        ! This is the matrix L (left eigenvector matrix - eq 13)
        e(1,1)=1.0
        if(abs(u(i) - c(i)) < TOLERANCE) then
            c(i)=c(i)+0.00001
        end if
        e(1,2)=-1.0/(u(i)-c(i))
        e(2,1)=1.0
        e(2,2)=-1.0/(u(i)+c(i))

        ! L^{-1} (inverse of Left eigenvector matrix)
        f_mat(1,1)=-(u(i)-c(i))/(2.0*c(i))
        f_mat(1,2)=(u(i)+c(i))/(2.0*c(i))
        f_mat(2,1)=-(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))
        f_mat(2,2)=(u(i)*u(i)-c(i)*c(i))/(2.0*c(i))

        ! Diagonal wave matrix D (eq 12)
        d(1)=abs(u(i)+c(i))
        d(2)=abs(u(i)-c(i))

        ! Equation 11 (L^{-1} D L)
        a(1,1)=e(1,1)*f_mat(1,1)*d(1)+e(2,1)*f_mat(1,2)*d(2)
        a(1,2)=e(1,2)*f_mat(1,1)*d(1)+e(2,2)*f_mat(1,2)*d(2)
        a(2,1)=e(1,1)*f_mat(2,1)*d(1)+e(2,1)*f_mat(2,2)*d(2)
        a(2,2)=e(1,2)*f_mat(2,1)*d(1)+e(2,2)*f_mat(2,2)*d(2)

        if(ots == 1) then
           dt(i)=cfl*dx(i-1)/max(d(1),d(2))
        else
           dt(i)=dtini
        endif

        ! Calculating dK/dA (eq 15)
        !ter=bo(i)+2.0*area(i)/bo(i)
        !Nazmul
        ter=pere(i)

        !Nazmul: dkda is corrected according to October 1, 2019 meeting at NWC
        dkda=sk(i)*((5.0/3.0*area(i)**(2.0/3.0)*ter)     &
                   -(2.0/3.0*area(i)**(5.0/3.0)*dpda(i)))/(ter**(5.0/3.0))

        ! Matrix S (eq 14)
        st(1, 1)=0.0
        st(1, 2)=0.0 !CHANGE
        !st(1, 2)=1.0
        st(2, 1)=grav*area(i)/bo(i)/bo(i)*dbdx(i)+gso(i)      &
                +f*2.0*grav*area(i)*oldQ(i)*abs(oldQ(i))/co(i)**3.0*dkda
                !+f*2.0*grav*area(i)*q(n,i)*abs(q(n,i))/co(i)**3.0*dkda

        !Nazmul: st(2,2) term is multiplied with a gravity
        !st(2,2)=-2*f*q(n,i)*grav*area(i)/co(i)/co(i)
        st(2,2)=-2*f*oldQ(i)*grav*area(i)/co(i)/co(i)


        if(dx(i) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i)
            crmax=max(crmax,cour*max(d(1),d(2)))
            crmin=min(crmin,cour*max(d(1),d(2)))
        endif

        ! LHS of eq 7 - this is done after the next check in matrixc
        b11(i)=0.5-phi-theta*cour*a(1,1)-0.5*thes*st(1,1)*dt(i)
        b12(i)=-theta*cour*a(1,2)-0.5*thes*st(1,2)*dt(i)
        b21(i)=-theta*cour*a(2,1)-0.5*thes*st(2,1)*dt(i)
        b22(i)=0.5-phi-theta*cour*a(2,2)-0.5*thes*st(2,2)*dt(i)

        ! This is switched with the dx check about in matrixc...
        if(i == 1) then
            cour=dt(i)
        else if(dx(i-1) < TOLERANCE) then
            cour=dt(i)
        else
            cour=dt(i)/dx(i-1)
        endif

        g(1,1)=0.5+phi+theta*cour*a(1,1)-0.5*thes*st(1,1)*dt(i)
        g(1,2)=theta*cour*a(1,2)-0.5*thes*st(1,2)*dt(i)
        g(2,1)=theta*cour*a(2,1)-0.5*thes*st(2,1)*dt(i)
        g(2,2)=0.5+phi+theta*cour*a(2,2)-0.5*thes*st(2,2)*dt(i)

        g11inv(i)= g(2,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g12inv(i)=-g(1,2)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g21inv(i)=-g(2,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
        g22inv(i)= g(1,1)/(g(1,1)*g(2,2)-g(1,2)*g(2,1))

        !f1(i)=q(n,i)
        f1(i)=oldQ(i)
        !print *, ncomp, i, area(i)
        !f2(i)=q(n,i)*q(n,i)/area(i)+grav*ci1(i)
        f2(i)=oldQ(i)*oldQ(i)/area(i)+grav*ci1(i)

        if(i >= 2 .and. i < ncomp) then
            dip1=area(i+1)/bo(i+1)
!           Nazmul: WHAT IS di? Need to check how to change it for natural channel
            di=2*area(i)/bo(i)
            dim1=area(i-1)/bo(i-1)
            ! to do: CHECK LOOP
            if (abs(dip1 + di + dim1) < TOLERANCE) then
                eps2(i) = 0.0
            else
                eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
            endif
        endif

        !if(i >= 2 .and. i < ncomp) then
        !    dip1=area(i+1)/bo(i+1)
        !    di=2*area(i)/bo(i)
        !    dim1=area(i-1)/bo(i-1)
        !    eps2(i)=alfa2*abs(dip1-di+dim1)/(dip1+di+dim1)
            ! print *, i,eps2(i),area(i),area(i+1),area(i-1)
        !endif
    end do

    eps2(1)=eps2(2)
    eps2(ncomp)=eps2(ncomp-1)

    do i=2,ncomp-1
        if(ityp(i).ne.1) then
            eps2(i)=eps2(i-1)
            eps2(i+1)=eps2(i+2)
        endif
    end do

    do i=2,ncomp-1
        eps2(i)=max(eps2(i+1),eps2(i))
        ! u(i+1)=(u(i+1)+u(i))/2.0
        ! c(i+1)=(c(i+1)+c(i))/2.0
        eps4(i)=max(0.,alfa4-eps2(i)/(u(i)+c(i)))
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
        d(1)=abs(u(i+1)+c(i+1))
        d(2)=abs(u(i+1)-c(i+1))
        ei1=max(d(1),d(2))
        eia=(ei+ei1)/2.0
        if(ityp(i) /= 1) then
            d1(i)=0.0
            d2(i)=0.0
        elseif(i.eq.2.or.i.eq.(ncomp-1)) then
            d1(i)=eps2(i)*eia*(area(i+1)-area(i))
            !d2(i)=eps2(i)*eia*(q(n,i+1)-q(n,i))
            d2(i)=eps2(i)*eia*(oldQ(i+1)-oldQ(i))
            ! print *, i,d1(i),d2(i),eps2(i),area(i+1),area(i)
        else
            d1(i)=eps2(i)*eia*(area(i+1)-area(i))-eps4(i)*(area(i+2)-3*area(i+1)+3*area(i)-area(i-1))
            !d2(i)=eps2(i)*eia*(q(n,i+1)-q(n,i))-eps4(i)*(q(n,i+2)-3*q(n,i+1)+3*q(n,i)-q(n,i-1))
            d2(i)=eps2(i)*eia*(oldQ(i+1)-oldQ(i))-eps4(i)*(oldQ(i+2)-3*oldQ(i+1)+3*oldQ(i)-oldQ(i-1))
        endif
        ! print *, i,d1(i),d2(i),eps2(i),area(i+1),area(i)
    end do
    !print*, 'eps2_pred', (eps2(i), i=1, ncomp)

end subroutine matrixp_R06
