subroutine section_R03()

    use constants_module_R03
    use arrays_module_R03
    use var_module_R03
    use arrays_section_module_R03
    use xsec_attribute_module_R03

    implicit none

    ! Locals
    !integer, intent(in) :: n
    integer :: i, pp, tableLength
    real(kind=4) :: beds, fs, yyn, yyn_1, temp1, temp2, new_I2
    real(kind=4) :: xt, q_sk_multi, currentQ
    real(kind=4) :: r_interpol, r_interpo_nn

    !change 20191122
    ci1=0
    ci2=0

    do i=1,ncomp

        if(i .gt. 1) then
            upstreamI2Tablec = I2Tablec
            upstreamSquareDepth = currentSquareDepth
        end if

        depth(i)=oldY(i)-z(i)

!      Nazmul change: read all attributes from tab file
        elevTable = xsec_tab(1,:,i)
        areaTable = xsec_tab(2,:,i)
        pereTable = xsec_tab(3,:,i)
        rediTable = xsec_tab(4,:,i)
        convTable = xsec_tab(5,:,i)
        topwTable = xsec_tab(6,:,i)
        nwi1Table = xsec_tab(7,:,i)
        dPdATable = xsec_tab(8,:,i)
        currentSquareDepth=(elevTable-z(i))**2
        I2Tablep = xsec_tab(9,:,i)
        I2Tablec = xsec_tab(10,:,i)

!     interpolate the cross section attributes based on water elevation
!		pause 301
        xt=oldY(i)

        area(i)=r_interpol(elevTable,areaTable,nel,xt)
!		pause 302
        !area(i)=r_interpol(x_tab,areaTable,nel,xt)
        pere(i)=r_interpol(elevTable,pereTable,nel,xt)
!		pause 303
        hy(i)=r_interpol(elevTable,rediTable,nel,xt)
        co(i)  =r_interpol(elevTable,convTable,nel,xt)
!		pause 304
        bo(i)  =r_interpol(elevTable,topwTable,nel,xt)
!        ci1(i) =r_interpol(elevTable,nwi1Table,nel,xt)
        ci1(i) =r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2)
!		pause 305
        dpda(i)=r_interpol(elevTable,dPdATable,nel,xt)
!		pause 306

        currentQ = oldQ(i)
        q_sk_multi = 1.0
        do pp = 1, size(Q_sk_tableEntry)
            if (  ( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
				tableLength = Q_sk_tableEntry(pp)
                q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp),Q_Sk_Table(2,1:tableLength,pp),tableLength,currentQ)
            end if
        end do
        co(i) = q_sk_multi*co(i)


! ----------------------------------------------------
        if(i .gt. 1) then

! I2 calculated as interpolation start
        yyn=oldY(i)
        yyn_1=oldY(i-1)
        temp1 = r_interpol(currentSquareDepth,I2Tablep,nel,(depth(i))**2)
        temp2 = r_interpol(upstreamSquareDepth, upstreamI2Tablec, nel, (depth(i-1))**2)
        new_I2 = (temp1+temp2)/2.0
! I2 calculated as interpolation end

            if(ityp(i-1) == 1) then
                 ci2(i)=new_I2
                 beds=(z(i-1)-z(i))/dx(i-1)
                 fs=f*0.5*oldQ(i-1)*abs(oldQ(i-1))/(co(i-1)*co(i-1))+f*0.5*oldQ(i)*abs(oldQ(i))/(co(i)*co(i))
                 aso(i)=(area(i)+area(i-1))/2.0*(beds-fs)
                 gso(i)=grav*(beds-fs)
                 dbdx(i)=(bo(i)-bo(i-1))/dx(i-1)
            end if
        end if
    end do

end subroutine section_R03


