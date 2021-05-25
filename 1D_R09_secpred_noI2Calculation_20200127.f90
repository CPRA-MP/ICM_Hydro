subroutine secpred_R09()

    use constants_module_R09
    use arrays_module_R09
    use var_module_R09
    use arrays_section_module_R09
    use xsec_attribute_module_R09

    implicit none

    ! Locals
    integer :: i, pp, tableLength
    real(kind=4) :: beds, fs, yyn, yyn_1, temp1, temp2, new_I2
    real(kind=4) :: xt, q_sk_multi, currentQ 
    real(kind=4) :: r_interpol, r_interpo_nn

    !change 20191122
    ci1=0
    ci2=0


    do i=ncomp,1,-1

        if(i .lt. ncomp) then
            downstreamI2Tablep = I2Tablep
            downstreamSquareDepth = currentSquareDepth
            yyn_1 = yyn
        end if
! ----------------------------------------------------
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

!      interpolate the cross section attributes based on predicted area
        xt=areap(i)

        ! NEW CHANGE: 20191119
        ! To get rid of negative predicted area
        if (areap(i) .le. TOLERANCE) then
            yyn = elevTable(1)+(elevTable(2)-elevTable(1))/100*5
            !areap(i) = 0.01
            print*, 'At i = ', i , ' pred area is negative'
        else
            yyn    =r_interpol(areaTable,elevTable,nel,xt)
        end if

        depth(i) = yyn - z(i)
        xt = yyn

        pere(i)=r_interpol(elevTable,pereTable,nel,xt)
        hy(i)     =r_interpol(elevTable,rediTable,nel,xt)
        co(i)  =r_interpol(elevTable,convTable,nel,xt)
        bo(i)  =r_interpol(elevTable,topwTable,nel,xt)
!        ci1(i) =r_interpol(areaTable,nwi1Table,nel,xt)
        ci1(i) =r_interpol(currentSquareDepth,nwi1Table,nel,(depth(i))**2)
        dpda(i)=r_interpol(elevTable,dPdATable,nel,xt)

        currentQ = qp(i)
		!pause 22
        q_sk_multi = 1.0
        do pp = 1, size(Q_sk_tableEntry)
            if (( eachQSKtableNodeRange(1,pp) - i) * ( eachQSKtableNodeRange(2,pp) - i) .le. 0 ) then
				tableLength = Q_sk_tableEntry(pp)
				!PRINT*, tableLength
				!print*, currentQ
				!print*, Q_Sk_Table(1,1:tableLength,pp)
				!print*, Q_Sk_Table(2,1:tableLength,pp)
                q_sk_multi = r_interpo_nn(Q_Sk_Table(1,1:tableLength,pp),Q_Sk_Table(2,1:tableLength,pp),tableLength,currentQ)
            end if
        end do
        co(i) = q_sk_multi*co(i)
		!pause 33
! ----------------------------------------------------
        if(i .lt. ncomp) then

! I2 opposite direction calculated as interpolation start
        xt=areap(i+1)
        temp1 = r_interpol(currentSquareDepth,I2Tablec,nel,(yyn-z(i))**2)
        temp2 = r_interpol(downstreamSquareDepth,downstreamI2Tablep, nel,(yyn_1-z(i+1))**2)
        new_I2 = (temp1+temp2)/2.0
! I2 opposite direction calculated as interpolation end

            if(ityp(i) == 1) then
                ci2(i)=new_I2
                beds=(z(i)-z(i+1))/dx(i)
                fs=f*0.5*qp(i)*abs(qp(i))/(co(i)**2)+f*0.5*qp(i+1)*abs(qp(i+1))/(co(i+1)**2)
                aso(i)=(areap(i)+areap(i+1))/2.0*(beds-fs)
                gso(i)=grav*(beds-fs)
                dbdx(i)=(bo(i+1)-bo(i))/dx(i)
            end if
        end if
    end do

end subroutine secpred_R09
