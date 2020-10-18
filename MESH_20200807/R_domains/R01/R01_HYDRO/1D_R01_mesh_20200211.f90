!            FINITE DIFFERENCE METHOD
!
!  A program for one dimensional flow in open channel
!

!subroutine end_R01
!    close(8)
!    close(9)
!    close(51)
!	close(61)
!end subroutine end_R01

subroutine init_R01(npr, ifile, input_file, nlat, outa, outb, outc, outd, wl_lat, Q_terminal)
    use constants_module_R01
    use arrays_module_R01
    use arrays_section_module_R01
    use var_module_R01
    use matrix_module_R01
    use sgate_module_R01
    use xsec_attribute_module_R01

    implicit none
	integer, intent(in) :: npr, ifile, nlat
	character(len=128), intent(in) :: input_file
	real(kind=4), dimension(nlat), intent(out) :: wl_lat
	real(kind=4), intent(out) :: Q_terminal
	real(kind=4), dimension(npr), intent(out) :: outa, outb, outc, outd
	

    ! Local storage
    integer(kind=4) :: i, j, n, pp
    real(kind=4) :: cour, da, dq, dxini, x, saveInterval
    real(kind=4) :: qn, xt, r_interpol, maxCourant
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds, latFlowValue

    real(kind=4) :: t, tfin, t1, t2, t0 !t0 start time

    character(len=128) :: upstream_path , downstream_path
    character(len=128) :: manning_strickler_path, output_path, other_input
    character(len=128) :: QSKtablePath, lateralFlow_path
    character(len=128) :: path

    ! open file for input data
    character(len=128) :: dx_path



     open(unit=999,file=trim(input_file))
 

    print*, 'R01: Reading input file'

    ! read data
    read(999,*) dtini ; print*,dtini    ! in seconds
    read(999,*) dxini
    read(999,*) t0        ! in hours
    read(999,*) tfin      ! in hours
	ntim = floor( (tfin - t0) / dtini * 3600)
	read(999,*) ncomp

	!np_R01=ncomp
    !read(999,*) ntim
	!nt_R01=ntim
    read(999,*) phi
    read(999,*) theta
    read(999,*) thetas
    read(999,*) thesinv
    read(999,*) alfa2
    read(999,*) alfa4
    read(999,*) f
    read(999,*) skk
    read(999,*) yy
    read(999,*) qq
    read(999,*) cfl
    read(999,*) ots
    read(999,*) yw
    read(999,*) bw
    read(999,*) w
    read(999,*) option
    read(999,*) yn
    read(999,*) qn
    read(999,*) igate
    read(999,*) xSection_path
    read(999,*) manning_strickler_path
    read(999,*) upstream_path
    read(999,*) downstream_path
    read(999,*) QSKtablePath
    read(999,*) dx_path
    read(999,*) lateralFlow_path
    read(999,*) output_path
    read(999,*) option_dsbc
    read(999,*) maxTableLength
    read(999,*) nel
    read(999,*) timesDepth
    read(999,*) other_input
    read(999,*) boundaryFileMaxEntry
    read(999,*) saveInterval; saveFrequency = saveInterval / dtini
    read(999,*) noLatFlow

    allocate(latFlowLocations(noLatFlow))   ! all the first nodes where a lateral flow starts
    allocate(latFlowType(noLatFlow))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
    allocate(latFlowXsecs(noLatFlow))       ! no of x-secs at the downstream that the lateral flow is applied

    if(noLatFlow>0)read(999,*) (latFlowLocations(i), i=1, noLatFlow)
    do i=1,noLatFlow
        if ((latFlowLocations(i)-1)*(latFlowLocations(i)-ncomp) .eq. 0) then
            print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
            stop
        end if
    end do
    if(noLatFlow>0)read(999,*) (latFlowType(i), i=1, noLatFlow)
    do i=1,noLatFlow
        if (latFlowType(i) .eq. 1) then
            print*, 'Lateral flow at node = ', latFlowLocations(i), ', is a time series'
        elseif (latFlowType(i) .eq. 2) then
            print*, 'Lateral flow at node = ', latFlowLocations(i), ', is a function of upstream flow'
		elseif (latFlowType(i) .eq. 3) then
                print*, 'Lateral flow at node = ', latFlowLocations(i), ', is a coupling point with ICM model'
        else
            print*, 'Wrong lateral flow type is provided. Type ', latFlowType(i), 'is not a valid type'
            stop
        end if
    end do

    if(noLatFlow>0)read(999,*) (latFlowXsecs(i), i=1, noLatFlow)

	read(999,*) noQSKtable

	! Addition: Multiple Q-Sk table
    allocate(eachQSKtableNodeRange(2,noQSKtable))

    read(999,*) (eachQSKtableNodeRange(1,i), i=1, noQSKtable)
    read(999,*) (eachQSKtableNodeRange(2,i), i=1, noQSKtable)

    !!! Need to test so that one section does not corresponds to more than one table
    do i = 2, noQSKtable
        if ( eachQSKtableNodeRange(2,i-1) .ge. eachQSKtableNodeRange(1,i) ) then
            print*, 'Wrong range of nodes applied for Q-Sk table.'
            print*, 'Lower limit of Table ', i-1,'must be smaller than the upper limit of Table ', i
            stop
        end if
    end do

    close(999)

    ! Allocate arrays
    call setup_arrays(ntim, ncomp, maxTableLength, boundaryFileMaxEntry, noLatFlow, noQSKtable)
    call setup_arrays_section
    call setup_xsec_attribute_module(nel, ncomp)

    dt = dtini

    open(unit=999, file=trim(dx_path))
    do i=1,ncomp-1
        read(999, *) x, dx(i)
    end do
    close(999)

    ! reading Strickler's coefficient at each section
    open(unit=999,file=trim(manning_strickler_path), status='unknown') !! //'Mannings_Stricklers_coeff.txt', status='unknown')
    do i=1,ncomp
        read(999, *) x, sk(i)
        call readXsection_R01(i,(1.0/sk(i)),timesDepth)
        ! This subroutine creates attribute table for each cross sections and saves in the hdd
        ! setting initial condition
        !y(1,i) = yy ! + z(i)
        oldY(i) = yy ! + z(i)
    end do
    close(999)

    do i=1,ncomp
        call create_I2_R01(i,ncomp)
    end do

    ityp = 1


    oldQ = qq


    ! reading Q-Strickler's coefficient multiplier table
    do i=1,noQSKtable
        write(file_num,'(i4.4)')i
        open(999,file=trim(QSKtablePath)//'Q_Mannings_table_'//file_num//'.txt')
        do n=1,maxTableLength
            read(999,*,end=300) Q_Sk_Table(1, n, i), Q_Sk_Table(2, n, i)
        end do
300     close(999)
        Q_sk_tableEntry(i) = n-1
    end do
    !print*, Q_sk_tableEntry


    x = 0.0

    ! Read hydrograph input Upstream
    open(unit=999, file=upstream_path)
    do n=1,boundaryFileMaxEntry
        read(999,*,end=301) USBoundary(1, n), USBoundary(2, n)
    end do
301 close(999)
    ppp = n-1

    ! Read hydrograph input Downstream
    open(unit=999, file=downstream_path)
    do n=1,boundaryFileMaxEntry
      read(999,*,end=302)  DSBoundary(1, n), DSBoundary(2, n)
    end do
302 close(999)
    qqq = n-1

    t=t0*60.0     !! t is in minute
	!t = 315532800.
	!t=10368000.

	!oldY=0. ! Hu set init WL=downstream bnd value
	!pause 11
	call section_R01()
	!pause 22
	oldArea=area
	!pause 33
	!HU

    ! applying boundary
    ! interpolation of boundaries at the initial time step

    oldQ(1)    =r_interpol(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t) ! t is in min
    oldY(ncomp)=r_interpol(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t) ! t is in min

    ! DS Boundary treatment: from water level to area time series
    ncompElevTable = xsec_tab(1,:,ncomp)
    ncompAreaTable = xsec_tab(2,:,ncomp)

	!open(unit=81,file=trim(output_path)//'DS_area.txt', status='unknown')
    xt=oldY(ncomp)
    oldArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)
    !write(81, *) t, oldArea(ncomp)



    ! read lateral flow conditions
    do i=1,noLatFlow
        write(file_num,'(i4.4)')latFlowLocations(i)
        open(999,file=trim(lateralFlow_path)//'lateral_'//file_num//'.txt')
        do n=1,boundaryFileMaxEntry
            read(999,*,end=303) lateralFlowTable(1, n, i), lateralFlowTable(2, n, i)
        end do
303     close(999)
        dataInEachLatFlow(i) = n-1
    end do


    ! Open files for output
    path = trim(output_path) // 'output_wl.txt'
    open(unit=ifile, file=trim(path), status='unknown')
    path = trim(output_path) // 'q.txt'
    open(unit=ifile+1, file=trim(path), status='unknown')

    path = trim(output_path) // 'area.txt'
    open(unit=ifile+2, file=trim(path), status='unknown')
    path = trim(output_path) // 'hy.txt'
    open(unit=ifile+3, file=trim(path), status='unknown')

    ! Output initial condition

    write(ifile, 10)  t, (oldY(i), i=1,ncomp)
    write(ifile+1, 10)  t, (oldQ(i), i=1, ncomp)
    write(ifile+2, 10) t, (oldArea(i), i=1, ncomp)
    write(ifile+3, 10) t, (hy(i), i=1, ncomp)

	outa=oldY
	outb=oldQ
	outc=oldArea
	outd=hy
	Q_terminal=oldQ(ncomp)
	wl_lat=-99999.
	do i=1,noLatFlow
		wl_lat(i)=oldY(latFlowLocations(i))
	enddo
	
!	pause 33
!10  format(f12.2 , 1200f12.2)
10  format(f12.2 , <ncomp>f12.2)

    do i=1,ncomp  !add zw for code debugging 10/15/2020
<<<<<<< HEAD
	   if(isNan(outa(i)) .OR. isNan(outb(i)). OR. isNan(outc(i)). OR. isNan(outd(i)))then
=======
	   if(isNan(outa(i)) .OR. isNan(outb(i)) .OR. isNan(outc(i)) .OR. isNan(outd(i)))then
>>>>>>> parent of 46ac90b... Revert "Create 1D_R01_mesh_20200211.f90"
	       write(*,*) 'NaNs after init_R'
		   write(*,*) 'at 1D Cross-section ',i
		   write(*,*) 'y_R= ',outa(i)
		   write(*,*) 'q_R= ',outb(i)
		   write(*,*) 'area_R= ',outc(i)
		   write(*,*) 'hy_R=',outd(i)
		   pause
	   endif
    enddo
	
end subroutine init_R01



subroutine cal_R01(n, npr, ifile, nlat, outa, outb, outc, outd, wl_lat, WL_terminal_from_ICM, Q_upstream_from_ICM, Q_lat_from_ICM, Q_terminal)

    use constants_module_R01
    use arrays_module_R01
    use arrays_section_module_R01
    use var_module_R01
    use matrix_module_R01
    use sgate_module_R01
    use xsec_attribute_module_R01

    implicit none
	
	integer, intent(in) :: n, npr, ifile, nlat
	real(kind=4), dimension(nlat), intent(in) :: Q_lat_from_ICM
	real(kind=4), dimension(nlat), intent(out) :: wl_lat
	real(kind=4), intent(in) :: WL_terminal_from_ICM, Q_upstream_from_ICM
	real(kind=4), intent(out) :: Q_terminal
	real(kind=4), dimension(npr), intent(out) :: outa, outb, outc, outd

    ! Local storage
    integer :: i, j, pp !, areanew
    real(kind=4) :: cour, da, dq, dxini, x, tfin
    real(kind=4) :: t, qn, xt, r_interpol, t1, t2, maxCourant, t0 !t0 start time in s
    real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, s0ds !, latFlowValue

    character(len=128) :: upstream_path , downstream_path
    character(len=128) :: bed_elevation_path, output_path, other_input
    character(len=128) :: channel_width_path, lateralFlow_path
    character(len=128) :: path

    ! open file for input data
    character(len=128) :: dx_path


!	t=0
!	print*, 'dtini=',dtini
!	pause 0
    !
    ! Loop on time
    !
	t=n*dtini/60.    ! t is in min, dtini is in sec

        ! interpolation of boundaries at the desired time step
		!pause 1
        newQ(1)     =r_interpol(USBoundary(1, 1:ppp),USBoundary(2, 1:ppp),ppp,t+dtini/60.) ! t is in min, dtini is in sec
        newY(ncomp) =r_interpol(DSBoundary(1, 1:qqq),DSBoundary(2, 1:qqq),qqq,t+dtini/60.) ! t is in min, dtini is in sec
		! coupling if WL_lat_from_ICM has a valid value
		if (WL_terminal_from_ICM > -9999.) newY(ncomp) = WL_terminal_from_ICM
		! End coupling
		! coupling if Q_upstream_from_ICM has a valid value
		if (Q_upstream_from_ICM > -9999.) newQ(1) = Q_upstream_from_ICM
		! End coupling
        xt=newY(ncomp)
		newArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)
		!pause 2
        if(isNan(newQ(1)).OR.isNan(newY(ncomp)).or.isNan(newY(ncomp)).OR.isNan(newArea(ncomp)))then
		    write(*,*)'NaNs in Cal_R for boundary assignment at time ',t
			write(*,*)'newQ(1)=',newQ(1)
			write(*,*)'newY(ncomp)=',newY(ncomp)
			write(*,*)'newArea(ncomp)=',newArea(ncomp)
			pause
	    endif
		
		! applying lateral flow at the oldQ
		lateralFlow = 0
        do i=1,noLatFlow
            if (latFlowType(i) .eq. 1) then
                lateralFlow(latFlowLocations(i)) = r_interpol(lateralFlowTable(1, :, i), &
                    lateralFlowTable(2, :, i),dataInEachLatFlow(i),t)
            elseif (latFlowType(i) .eq. 2) then
                lateralFlow(latFlowLocations(i)) = r_interpol(lateralFlowTable(1, :, i), &
                    lateralFlowTable(2, :, i),dataInEachLatFlow(i),oldQ(latFlowLocations(i)))
			elseif (latFlowType(i) .eq. 3) then
                  lateralFlow(latFlowLocations(i)) = Q_lat_from_ICM(i)             ! Q coupling from ICM to MESH; (+)ve Q adds discharge to MESH
            endif
                lateralFlow(latFlowLocations(i)) = lateralFlow(latFlowLocations(i))/ &
                    sum(dx(latFlowLocations(i)-1:latFlowLocations(i)-1+latFlowXsecs(i)-1)) !update1

            do j=1,latFlowXsecs(i)-1
                lateralFlow(latFlowLocations(i)+j)=lateralFlow(latFlowLocations(i))
            end do
            !print*, 'lat flow=', latFlowValue, latFlowLocations(i), &
            !    latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5, &
            !    oldQ(latFlowLocations(i))
            ! print*, 'lat flow at', latFlowLocations(i), &
            !    'flow=',lateralFlow(latFlowLocations(i))*0.5*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i))), &
            !    dx(latFlowLocations(i)-1), dx(latFlowLocations(i))
            if(isNan(lateralFlow(latFlowLocations(i))))then
			   write(*,*)'time step=',t
			   write(*,*)'NaN lateral flow from Cal-R at 1D location ', latFlowLocations(i)
			   write(*,*)'Q_lat_from_ICM= ', Q_lat_from_ICM(i)
			   write(*,*)'LatFlowLoc, Q_lat_from_ICM'
               do j=1,noLatFlow
			       write(*,*)latFlowLocations(j),Q_lat_from_ICM(j)
			   enddo
			   pause
			endif
		end do
        !print*, 'noLatFlow',noLatFlow
        !print*, 'lateralFlow', (lateralFlow(i), i=1, ncomp)
		

        ! Set upstream discharge
        dqp(1) = newQ(1) - oldQ(1)
		!pause 3
		if(isNan(dqp(1)))then
		   write(*,*)'NaN value for dqp(1) in Cal_R at time ',t
		   write(*,*)'oldQ(1)=',oldQ(1)
		   write(*,*)'newQ(1)=',newQ(1)
		   pause
		endif
		
        call section_R01()
        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
        thes=thetas
		!pause 4

        call matrixp_R01()
		!pause 5

        do i=2,ncomp
            !cour=dt(i)/dx(i-1)
            cour=dtini/dx(i-1)
            !rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i)*dtini*dx(i)/dx(i-1) <<wrong
			rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i)*dtini !updated2
            rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dtini*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
            c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
            c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
            c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
            dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
            dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
           ! print*,'rhs1=', rhs1, rhs2, c11, c12, c21, c22
        end do

       ! print*, 'f1(i)', (f1(i), i=1, ncomp)
       ! print*, 'd1(i)', (d1(i), i=1, ncomp)
       ! print*, 'g11inv(i)', (g11inv(i), i=1, ncomp)
       ! print*, 'g12inv(i)', (g12inv(i), i=1, ncomp)
       ! print*, 'g21inv(i)', (g21inv(i), i=1, ncomp)
       ! print*, 'g22inv(i)', (g22inv(i), i=1, ncomp)
       ! print*, 'dap(i)', (dap(i), i=1, ncomp)
       ! print*, 'dqp(i)', (dqp(i), i=1, ncomp)
        ! Boundary conditions at downstream (right boundary)
        if (option_dsbc.eq.1) then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            arean=yn*bo(ncomp)
            perimn=2.0*yn+bo(ncomp)
            hyrdn=arean/perimn
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*arean*hyrdn**(2.0/3.0)*sqrt(s0ds)
            dqp(ncomp)=qn-oldQ(ncomp)
            dqc(ncomp)=dqp(ncomp)

        elseif(option_dsbc.eq.2)then
            dac(ncomp)=dap(ncomp)
            yn=(area(ncomp)+dap(ncomp))/bo(ncomp)
            areac=yn*bo(ncomp)
            perimc=2.0*yn+bo(ncomp)
            hyrdc=areac/perimc
            s0ds=-((z(ncomp)-z(ncomp-1))/dx(ncomp))
            qn=skk*areac*hyrdc**(2.0/3.0)*sqrt(s0ds)
            qcrit=1.05*(((yn**3.0)*(bo(ncomp)**2.0)*grav)**(1.0/2.0))
            write(*,*)qcrit
            dqp(ncomp)=qcrit-oldQ(ncomp)
            dqc(ncomp)=dqp(ncomp)

        else
            !dac(ncomp)=0.0
            !dap(ncomp)=0.0
            !dqc(ncomp)=dqp(ncomp)

!            dap(ncomp)=0.0	!checked email !for critical comment out
! change for unsteady flow
			dap(ncomp) = newArea(ncomp) - oldArea(ncomp)
            dac(ncomp)=dap(ncomp)	!checked email
            dqc(ncomp)=dqp(ncomp)	!checked email

        endif

        ! Update via predictor
        areap = area + dap
        qp = oldQ + dqp


        ! applying lateral flow at the qp
        !do i=1,noLatFlow
        !    latFlowValue = r_interpol(lateralFlowTable(1, :, i),lateralFlowTable(2, :, i),dataInEachLatFlow(i),t)
        !    qp(latFlowLocations(i)) = qp(latFlowLocations(i)) + &
        !        latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5
        !    print*, 'lat flow=', latFlowValue, latFlowLocations(i), &
        !        latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5, &
        !        qp(latFlowLocations(i))
        !end do

        !print*, 'areap', (areap(i), i=1, ncomp)


		!pause 6

        call secpred_R01()
        thes=thesinv
				!pause 7

        call matrixc_R01()
		!pause 8

        do i=ncomp-1,1,-1
            !cour=dt(i)/dx(i)
            cour=dtini/dx(i)
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i+1)*dtini !updated3
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dt(i)*grav*(ci2(i)+aso(i))
            c11=g11inv(i)*b11(i+1)+g12inv(i)*b21(i+1)
            c12=g11inv(i)*b12(i+1)+g12inv(i)*b22(i+1)
            c21=g21inv(i)*b11(i+1)+g22inv(i)*b21(i+1)
            c22=g21inv(i)*b12(i+1)+g22inv(i)*b22(i+1)
            dac(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dac(i+1)-c12*dqc(i+1)
            dqc(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dac(i+1)-c22*dqc(i+1)
        end do
        ! Upstream boundary condition
        ! Prescribed discharge at the upstream
        ! Area correction is calculated

        dqc(1)=dqp(1)	!checked email
        !dac(1)=dap(1) !for critical, uncomment
        dap(1)=dac(1)	!checked email

        ! Final update
        do i=1,ncomp
            da=(dap(i)+dac(i))/2.0
            dq=(dqp(i)+dqc(i))/2.0
            newArea(i)=da+area(i)
            if(newArea(i) <= 0.0) newArea(i)=0.001

!           Now calculate y based on area calculated
!-------------------------------------
            elevTable(:) = xsec_tab(1,:,i)
            areaTable(:) = xsec_tab(2,:,i)

    !       interpolate the cross section attributes based on FINAL CALCULATED area
            xt=newArea(i)
            newY(i)=r_interpol(areaTable,elevTable,nel,xt)
!-------------------------------------

            newQ(i)=oldQ(i)+dq
            froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

        end do
		!pause 9

        do i=1,ncomp-1
            courant(i)=(newQ(i)+newQ(i+1))/(newArea(i)+newArea(i+1))*dtini/dx(i)
        enddo

        if (maxCourant .lt. maxval (courant)) then
            maxCourant = maxval (courant)
        endif

        t = t0*60. + (n+1)*dtini/60.             ! t is in min, dtini is in sec
        !print "('- cycle',i9,'  completed')", n
		if(mod(n+1,24*saveFrequency) .eq. 0 .or. (n.eq.0))write(*,*)'R01: Nstep =', n+1, 'Days = ' &
		, t/60./24., real(n+1)/real(ntim)*100., '% completed'
        !print*, 'dqp', (dqp(i), i=1, ncomp)
        !print*, 'dqc', (dqc(i), i=1, ncomp)
        !print*, 'qp', (qp(i), i=1, ncomp)
        !print*, 'oldQ', (oldQ(i), i=1, ncomp)
        !print*, 'newQ', (newQ(i), i=1, ncomp)

        if (mod(n+1,saveFrequency) .eq. 0 .or. n .eq. (ntim-1)) then
        write(ifile, 10) t, (newY(i), i=1,ncomp)
        write(ifile+1, 10) t, (newQ(i), i=1,ncomp)
        write(ifile+2, 10) t, (newArea(i), i=1, ncomp)
        end if
        !  if (n .eq. 8000) pause

        !write(81, *) t, newArea(ncomp)

        ! update of Y, Q and Area vectors
        oldY   = newY
        oldQ   = newQ
        oldArea= newArea
		
!>> *********************************Adding 1d end

          !Q_1d = newQ(ncomp)                   ! 1D2D coupling 1d discharge

!          do i=1,noLatFlow
!             if (latFlowType(i) .eq. 3) then
!                  WL_lat_from_MESH = newY(latFlowLocations(i))                  ! WL lateral coupling from MESH to ICM
!              endif
!          end do

!>> -- Loop over links. Save calculated flowrates, Q as the initial condition for the next simulation timestep.
		
		
		
		call section_R01() ! Hu

	    if (mod(n+1,saveFrequency) .eq. 0 .or. n .eq. (ntim-1)) then
        write(ifile+3, 10) t, (hy(i), i=1, ncomp)
        end if

		outa=oldY
		outb=oldQ
		outc=oldArea
		outd=hy

		Q_terminal=oldQ(ncomp)
		wl_lat=-99999.
		do i=1,noLatFlow
			wl_lat(i)=oldY(latFlowLocations(i))
		enddo
		!pause 10
!10  format(f12.2 , 1200f12.2)
10  format(f12.2 , <ncomp>f12.2)

    do i=1,ncomp  !add zw for code debugging 10/15/2020
	   if(isNan(outa(i)) .OR. isNan(outb(i)). OR. isNan(outc(i)). OR. isNan(outd(i)))then
	       write(*,*) 'NaNs after Cal_R at time step ', n
		   write(*,*) 'at 1D Cross-section ',i
		   write(*,*) 'y_R= ',outa(i)
		   write(*,*) 'q_R= ',outb(i)
		   write(*,*) 'area_R= ',outc(i)
		   write(*,*) 'hy_R=',outd(i)
		   pause
	   endif
    enddo

	end subroutine cal_R01
