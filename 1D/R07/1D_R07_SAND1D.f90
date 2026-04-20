!-------------------------------------------------------------------------------
!
!  Developer: Vasilia Velissariou <vvelissariou@tulane.edu>
!
!  Version: 0.4
!
!    Version - 0.1 Thu Nov 21 2019
!            - Basic code development that uses an explicit
!              2nd order accurate (FTCS numerical scheme)
!    Version - 0.2 Mon Dec 16 2019
!            - Improvements on sediment source terms
!            - Development of sed_properties module that contains
!              many functions related to settling velocities, critical shear streses
!              van Rijn bottom erosion, etc
!    Version - 0.3 Tue Jan 28 2020
!            - General code cleanup
!            - Improved array manipulation for code speed-up and memory requirements
!            - Code modifications to account for variable (unequally spaced) cross-sections
!    Version - 0.4 Fri Jun 11 2021
!            - General code cleanup & updates to indentation/tabs
!            - Add sediment bed depth to source term to prevent over-scour
!-------------------------------------------------------------------------------
!




!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
subroutine init_SAND_R07(npr, ioutfile, input_file, rday, ndtr, nlatr_ori, nlatr, latFlowLoc, latFlowTp, latFlowXsec)
    use precision
    !use params
    use mod_arrays_SAND_R07
    
    use sed_properties_SAND_R07
    use spline_func
    
    implicit none
    integer, intent(in) :: npr, ioutfile, ndtr, nlatr_ori, nlatr
    real, intent(in) :: rday
    character(len=128), intent(in) :: input_file
    integer, dimension(nlatr), intent(in) :: latFlowLoc, latFlowTp, latFlowXsec
    integer :: nx, nt, np
    
    integer :: i, j, n, k, ilat, jlat
    
    real(fp) :: fval, fval1, fval2 !,t1,t2  !modu_user
    integer :: ival, ival1, ival2
    !integer :: cntUserDt, flgUserDt
    
    character(len=256) :: reach_dist, reach_flow, reach_cbnd, reach_depth, reach_area,  reach_qm  !out_dir,
    ! Downstream  bnd file /lateral flow
    character(len=256) :: reach_dbnd, reach_lat
    
    
    character(len=128) :: outSTR
    
    !integer :: nx, nt, np
    integer :: LUN, ioerr
    
    LUN = 999
    
    open(unit = LUN, file=trim(input_file)//'/input/sand.inp', &
    status='old', action = 'read', iostat = ioerr)
    
    if (ioerr /= 0) then
        write(*, *) 'ERROR reading the input file: sand.inp'
        stop
    endif
    
    !read(LUN, '(a)') RefDate   ! Reference date, Format: MM-DD-YYYY_HH:MM:SS or MM-DD-YYYYTHH:MM:SS (not used at the moment)
    read(LUN, *) sand_init   ! Reference date, Format: MM-DD-YYYY_HH:MM:SS or MM-DD-YYYYTHH:MM:SS (not used at the moment)
    !read(LUN, *) Tstart        ! Start time of the computation (Tstart [s]) - Relative to RefDate
    Tstart=0.
    !read(LUN, *) Tstop         ! End time of the computation (Tstop [s]) - Relative to RefDate
    Tstop=rday*86400.
    !read(LUN, *) Dt            ! Timestep to be used in the simulation (Dt [s])
    Dt=real(ndtr)               ! this is parameterized here but also used in cal_SAND_R## subroutine    
    read(LUN, *) DtUser        ! Write time interval [s]
    !read(LUN, *) Nclass        ! Number of sediment size classes (nclass)
    Nclass=1
    read(LUN, *) Ks            ! Longitudinal diffusion coefficient, reasonable starting value 10.0 m^2/s (between 10.0 and 100.0?)
    !read(LUN, *) D90x          ! representative D90/D50 ratio
    D90x=1.5
    !read(LUN, *) NDx           ! Number of cross-sections or number of spatial increments
    NDx=npr
    
    
    if ( comparereal(Dt, 0.0_fp) <= 0 ) then
        write(*, '(3x, a)') 'Dt should be a positive number'
        write(*, '(3x, a)') 'Please adjust its value in the input file'
        write(*, '(3x, a, f10.2)') 'Dt     = ', Dt
        stop
    endif
    
    if ( comparereal(DtUser, 0.0_fp) <= 0 ) then
        write(*, '(3x, a)') 'DtUser should be a positive number'
        write(*, '(3x, a)') 'Please adjust its value in the input file'
        write(*, '(3x, a, f10.2)') 'DtUser     = ', DtUser
        stop
    endif
    
    if ( comparereal(modulo(DtUser, Dt), 0.0_fp) /= 0 ) then
        write(*, '(3x, a)') 'DtUser should be an integral value of Dt: DtUser = Integer * Dt'
        write(*, '(3x, a)') 'Please adjust the values of Dt and/or DtUser in the input file'
        write(*, '(3x, a, f10.2)') 'Dt     = ', Dt
        write(*, '(3x, a, f10.2)') 'DtUser = ', DtUser
        stop
    endif
    
    
    NDt = (Tstop - Tstart) / Dt + 1
    NDtUser = floor( (Tstop - Tstart) / DtUser ) + 1
    
    
    if ( (NDt <= 1) .or. (NDtUser <= 1) ) then
        write(*, '(3x, a)') 'Tstop should be greater than Tstart'
        write(*, '(3x, a)') 'Please adjust the values of Tstart and/or Tstop in the input file'
        write(*, '(3x, a, f14.2)') 'Tstart = ', Tstart
        write(*, '(3x, a, f14.2)') 'Tstop  = ', Tstop
        stop
    endif
    
    
    nt = NDt
    nx = NDx
    np = nclass
    
    
    ! allocate memory for the arrays
    call aloc_arrays
    !allocate(mACLD(NDtUser, nx, np))
    !allocate(diffACLD(NDtUser, nx, np))
    
    !allocate(zz(nx))
    
    ! sand D50 (m)   -- USDA particle size definition: 0.05 mm < sand < 2 mm
    ! silt D50 (m)   -- USDA particle size definition: 0.002 mm < silt < 0.05 mm
    ! clay D50 (m)   -- USDA particle size definition: clay < 0.002 mm
    do i = 1, np
        read(LUN, *) D50(i)
    end do
    
    ! critical shear stress for sediment particles (N/m**2)
    do i = 1, np
        read(LUN, *) Tcrit(i)
    end do
    
    !read(LUN, *) Cf              ! friction coefficient (0.001 <= Cf <= 0.003), Cf = 0.0025 corresponds to Manning's 0.018
    read(LUN, *) Cmann            ! Manning's coefficient, a good value is 0.023
    read(LUN, *) alphaSED         ! van Rijn's suspended load coefficient (alphaSED = 0.012, dimensionless) for Sand
    read(LUN, *) SEDcalib         ! re-adaptibility coefficient (0.001 <= SEDcalib <= 0.012), SEDcalib = 0.0016 kg/m2*s for MUD
    read(LUN, *) SEDn             ! 1.0 <= SEDn <= 3.0??? (SEDn = 1.0)
    read(LUN, *) ws_fine		  ! settling velocity (m/s) for fine
    
    ! This file contains the distances of the reach increments
    read(LUN, *) reach_dist
    reach_dist=trim(input_file)//reach_dist
    
    ! This file contains the time dependent depths
    !read(LUN, *) reach_depth
    
    ! This file contains the time dependent depths, areas and flowrates at each segment
    !read(LUN, *) reach_flow
    
    ! This file contains the time dependent depths, areas and flowrates at each segment
    !read(LUN, *) reach_area
    
    ! This file contains the time dependent concentration boundary conditions for each sediment class
    read(LUN, *) reach_cbnd
    reach_cbnd=trim(input_file)//reach_cbnd
    
    ! This file contains the time dependent concentration boundary conditions for each sediment class
    read(LUN, *) reach_dbnd
    reach_dbnd=trim(input_file)//reach_dbnd
    
    ! This file contains dis_manning table
    read(LUN, *) reach_qm
    reach_qm=trim(input_file)//reach_qm
    
    ! output folder
    read(LUN, *) out_dir
    out_dir=trim(input_file)//out_dir
    
    ! This variable holds the type of the BCs to be used downstream (currently 1, 2, or 3)
    ! Recommended value is 3
    ! Use bc_option = 1 to impose the BC: Q*C - Ks*A*dC/dx = 0
    ! Use bc_option = 2 to impose the BC: Ks*A*dC/dx = 0 -> C(nx+1) = C(nx)  (no despersive transport)
    ! Use bc_option = 3 apply the 2nd order backward (upwind) in space at the downstream boundary
    !read(LUN, *) bc_option
    bc_option=3
    if( (bc_option /= 1) .and. (bc_option /= 2) .and. (bc_option /= 3) ) then
        write(*, '(3x, a)') 'Wrong value for bc_option: bc_option = 1,2 or 3'
        write(*, '(3x, a)') 'Please adjust the value of bc_option in the input file'
        write(*, '(3x, a, i4)') 'bc_option = ', bc_option
        stop
    endif
    
    ! This is used to turn the sediment source term on (>0) or off (<= 0)
    !read(LUN, *) sed_sterm
    sed_sterm=1
    if( sed_sterm > 0 ) then
        useSRC = .TRUE.
    else
        useSRC = .FALSE.
    endif
    
    ! This is used to turn bottom erosion in the source term on (>0) or off (<= 0)
    !read(LUN, *) bot_eros
    bot_eros=1
    if( bot_eros > 0 ) then
      useBOTEROS = .TRUE.
    else
      useBOTEROS = .FALSE.
    endif
    
    ! lateral flow input
    Nlat=nlatr_ori
    if( Nlat > 0 ) then
        allocate(Idlat(Nlat))
        allocate(Nsclat(Nlat))
        allocate(Ntplat(Nlat))
        !read(LUN, *) (Idlat(i),i=1, Nlat)
        !read(LUN, *) (Nsclat(i),i=1, Nlat)
        Idlat=latFlowLoc
        Ntplat=latFlowTp
        Nsclat=latFlowXsec
           ! This file contains the time dependent concentration boundary conditions for each lateral flow and con
        read(LUN, *) reach_lat 
        reach_lat=trim(input_file)//reach_lat
    endif
    !print*, Nlat
    !print*, Idlat
    !print*, Nsclat
    !print*, reach_lat
    !pause 1
    
    close (LUN)
    !!!!!!!!!!!!!!!!!!
    
    
    write(*, *) 'Sand_init    = ', sand_init
    write(*, *) 'Tstart     = ', Tstart
    write(*, *) 'Tstop      = ', Tstop
    write(*, *) 'Dt         = ', Dt
    write(*, *) 'DtUser     = ', DtUser
    write(*, *) 'Nclass     = ', Nclass
    write(*, *) 'Ks         = ', Ks
    write(*, *) 'D90x       = ', D90x
    write(*, *) 'NDt        = ', NDt
    write(*, *) 'NDtUser    = ', NDtUser
    write(*, *) 'NDx        = ', NDx
    
    do i = 1, np
        write(*, '(a, 1x, i2, 3x, a, 1x, f10.7)') ' D50 sed class', i, '=', D50(i)
    end do
    
    do i = 1, np
        write(*, '(a, 1x, i2, 1x, a, 1x, f10.7)') ' Tcrit sed class', i, '=', Tcrit(i)
    end do
    
    !write(*, *) 'CF         = ', Cf
    write(*, *) 'Cmann      = ', Cmann
    write(*, *) 'alphaSED   = ', alphaSED
    write(*, *) 'SEDcalib   = ', SEDcalib
    write(*, *) 'SEDn       = ', SEDn
    write(*, *) 'ws_fine       = ', ws_fine
    
    write(*, *) 'reach_dist = ', trim(adjustl(reach_dist))
    !write(*, *) 'reach_depth = ', trim(adjustl(reach_depth))
    !write(*, *) 'reach_flow = ', trim(adjustl(reach_flow))
    !write(*, *) 'reach_area = ', trim(adjustl(reach_area))
    write(*, *) 'reach_cbnd = ', trim(adjustl(reach_cbnd))
    write(*, *) 'out_dir = ', trim(adjustl(out_dir))
    
    
    write(*, *) 'BC_Option    = ', BC_Option
    write(*, *) 'SED_STerm    = ', SED_STerm
    write(*, *) '  useSRC     = ', useSRC
    write(*, *) 'BOT_Eros     = ', BOT_Eros
    write(*, *) '  useBOTEROS = ', useBOTEROS
    
    !pause 1
    !-------------------------------------------------------------------------------
    !----------  BEG:: READ THE INPUT DATA
    !-------------------------------------------------------------------------------
    
    
    !!!!! UNIT = 12: contains the 1D x array for the reach distance (space) = (x)
    LUN = 999
    open(unit = LUN, file = trim(reach_dist), status = 'old', action = 'read', iostat = ioerr)
    if (ioerr /= 0) then
        write(*, *) 'ERROR reading the input file: reach_dist = ' // trim(reach_dist)
        stop
    else
        read(LUN,*) !skip comment line  
        read(LUN, *) ival
        if (ival /= nx) then
            write(*, *) 'ERROR in the input file: reach_dist = ' // trim(reach_dist)
            write(*, *) '      Wrong number of spatial increments supplied'
            write(*, *) '      Please adjust the input file sed1d.inp'
            write(*, '(a, i4, a)') '       UserNDx = ', NDx, ' (from INI file)'
            write(*, '(a, i4, a)') '            nx = ', ival, ' (from DIST file)'
            close(LUN)
            stop
        end if
        do i = 1, nx
            read(LUN, *) fval
            xx(i) = fval
        end do
    end if
    close(LUN)
      
    !  print*,xx(1), xx(25)
    
    !!!!!!---------- BEG:: READ THE FLOW DATA FOR EACH CROSS SECTION
    !!!!! UNIT = 13: contains the 2D array data for the cross section flowrates (time, space) = (t, x)
      
      !print*, TimeFL(1,1), TimeFL(49,25)
      !print*, FlowFL(1,1), FlowFL(49,25)
      
      
    !!!!! UNIT = 13: contains the 2D array data for the cross section water depth (time, space) = (t, x)
    
      
      !print*, DepthFL(1,1), DepthFL(49,25)
      
    !!!!! UNIT = 13: contains the 2D array data for the cross section area (time, space) = (t, x)
     
    
      
      !print*, AreaFL(1,1), AreaFL(49,25)
      
      
      
      !pause 2  
    
    !!!!!!---------- END:: READ THE FLOW DATA FOR EACH CROSS SECTION
    
    !!!!!!---------- BEG:: READ THE UPSTREAM CONCENTRATION BOUNDARY DATA - Concentrations are imposed
    !!!!! UNIT = 14: contains the 1D array data for the concentration at the boundary (time) = (t)
    LUN = 999
    open(unit = LUN, file = trim(reach_cbnd), status = 'old', action = 'read', iostat = ioerr)
    if (ioerr /= 0) then
        write(*, *) 'ERROR reading the input file: reach_cbnd = ' // trim(reach_cbnd)
        stop
    else
        read(LUN,*) !Skip comment line
        read(LUN, *) ntcb, npcb
        if (npcb /= np) then
            write(*, *) 'ERROR in the input file: reach_cbnd = ' // trim(reach_cbnd)
            write(*, *) '      Wrong number of sediment classes were supplied'
            write(*, *) '      Please adjust the input file sed1d.inp'
            write(*, '(a, i4, a)') '       UserNclass = ', Nclass, ' (from INI file)'
            write(*, '(a, i4, a)') '           Nclass = ', npcb, ' (from CBND file)'
            close(LUN)
            stop
        end if
        ! Times [s] in the boundary files should be relative to the RefDate
        ! supplied in the sed1d.inp file. It is the user's responsibility
        ! to supply the correct boundary data (pre-process the original data)
        allocate(ConcBND(ntcb, npcb))
        allocate(TimeBND(ntcb, npcb))
        ConcBND = 0.0; TimeBND = 0.0
        
        do i=1, ntcb
            read(LUN,*)fval, (ConcBND(i,j),j=1,npcb)
            do j=1, npcb
                TimeBND(i,j)=fval*60.  
            enddo
        enddo
    end if
    close(LUN)
    
    ! Check if the times in the BC data are inside the (Tstart, Tstop) range
    ! If yes, then stop the calculations, because of possibly missing data at
    ! certain times
    do i = 1, np
        if ( (comparereal(Tstart, TimeBND(1, i)) < 0) .or. (comparereal(Tstop, TimeBND(ntcb, i)) > 0)) then
            write(*, *) 'ERROR in the input file: reach_cbnd = ' // trim(reach_cbnd)
            write(*, *) '      Data times fall inside the model''s simulation time range (Tstart, Tstop)'
            write(*, *) '      Please check the data times in the input file'
            write(*, '(a, f10.1, a)') '       Tstart = ', Tstart, ' (from INI file)'
            write(*, '(a, f10.1, a)') '        Tstop = ', Tstop, ' (from INI file)'
            stop
        endif
    end do
    !print*, TimeBND(1,1), TimeBND(2,1)
    !  print*, ConcBND(1,1), ConcBND(2,1)
    
    
    !!!Downstream bnd
    !!!!!!---------- BEG:: READ THE UPSTREAM CONCENTRATION BOUNDARY DATA - Concentrations are imposed
    !!!!! UNIT = 14: contains the 1D array data for the concentration at the boundary (time) = (t)
    LUN = 999
    open(unit = LUN, file = trim(reach_dbnd), status = 'old', action = 'read', iostat = ioerr)
    if (ioerr /= 0) then
        write(*, *) 'ERROR reading the input file: reach_dbnd = ' // trim(reach_dbnd)
        stop
    else
        read(LUN,*) !Skip comment line
        read(LUN, *) ntdb, npdb
        if (npdb /= np) then
            write(*, *) 'ERROR in the input file: reach_dbnd = ' // trim(reach_dbnd)
            write(*, *) '      Wrong number of sediment classes were supplied'
            write(*, *) '      Please adjust the input file sed1d.inp'
            write(*, '(a, i4, a)') '       UserNclass = ', Nclass, ' (from INI file)'
            write(*, '(a, i4, a)') '           Nclass = ', npcb, ' (from CBND file)'
            close(LUN)
            stop
        end if
        ! Times [s] in the boundary files should be relative to the RefDate
        ! supplied in the sed1d.inp file. It is the user's responsibility
        ! to supply the correct boundary data (pre-process the original data)
        allocate(ConcDOW(ntdb, npdb))
        allocate(TimeDOW(ntdb, npdb))
        ConcDOW = 0.0; TimeDOW = 0.0
        
        do i=1, ntdb
            read(LUN,*)fval, (ConcDOW(i,j),j=1,npdb)
            do j=1, npdb
              TimeDOW(i,j)=fval*60.  
            enddo
        enddo
    end if
    close(LUN)
    
    ! Check if the times in the BC data are inside the (Tstart, Tstop) range
    ! If yes, then stop the calculations, because of possibly missing data at
    ! certain times
    do i = 1, np
        if ( (comparereal(Tstart, TimeDOW(1, i)) < 0) .or. (comparereal(Tstop, TimeDOW(ntdb, i)) > 0)) then
            write(*, *) 'ERROR in the input file: reach_dbnd = ' // trim(reach_dbnd)
            write(*, *) '      Data times fall inside the model''s simulation time range (Tstart, Tstop)'
            write(*, *) '      Please check the data times in the input file'
            write(*, '(a, f10.1, a)') '       Tstart = ', Tstart, ' (from INI file)'
            write(*, '(a, f10.1, a)') '        Tstop = ', Tstop, ' (from INI file)'
            stop
        endif
    end do
    !  print*, TimeDOW(1,1), TimeDOW(2,1)
    !  print*, ConcDOW(1,1), ConcDOW(2,1)
    
    !!Lateral flow
    !!
    if (Nlat > 0)then
    !!!lateral flow
    !!!!! UNIT = 14: contains the 1D array data for the concentration at the boundary (time) = (t)
        LUN = 999
        open(unit = LUN, file = trim(reach_lat), status = 'old', action = 'read', iostat = ioerr)
        if (ioerr /= 0) then
            write(*, *) 'ERROR reading the input file: reach_lat = ' // trim(reach_lat)
            stop
        else
            read(LUN,*) !Skip comment line
            read(LUN, *) ntlat, nplat
            if (nplat /= Nlat) then
                write(*, *) 'ERROR in the input file: reach_lat = ' // trim(reach_lat)
                write(*, *) '      Wrong number of lateral flow were supplied'
                write(*, *) '      Please adjust the input file sed1d.inp'
                write(*, '(a, i4, a)') '       Nlat = ', Nlat, ' (from INI file)'
                write(*, '(a, i4, a)') '           nplat = ', nplat, ' (from lateral flow file)'
                close(LUN)
                stop
            end if
          ! Times [s] in the boundary files should be relative to the RefDate
          ! supplied in the sed1d.inp file. It is the user's responsibility
          ! to supply the correct boundary data (pre-process the original data)
            allocate(ConcLAT(ntlat, nplat))
            allocate(FlowLAT(ntlat, nplat))
            allocate(TimeLAT(ntlat, nplat))
    	    
            ConcLAT = 0.0; TimeLAT = 0.0
    	    FlowLAT = 0.0
            
            do i=1, ntlat
                read(LUN,*)fval, (FlowLAT(i,j), ConcLAT(i,j),j=1,nplat)
                do j=1, nplat
                    TimeLAT(i,j)=fval*60.  
                enddo
            enddo
        end if
        close(LUN)
    
        ! Check if the times in the BC data are inside the (Tstart, Tstop) range
        ! If yes, then stop the calculations, because of possibly missing data at
        ! certain times
        do i = 1, np
            if ( (comparereal(Tstart, TimeLAT(1, i)) < 0) .or. (comparereal(Tstop, TimeLAT(ntlat, i)) > 0)) then
                write(*, *) 'ERROR in the input file: reach_lat = ' // trim(reach_lat)
                write(*, *) '      Data times fall inside the model''s simulation time range (Tstart, Tstop)'
                write(*, *) '      Please check the data times in the input file'
                write(*, '(a, f10.1, a)') '       Tstart = ', Tstart, ' (from INI file)'
                write(*, '(a, f10.1, a)') '        Tstop = ', Tstop, ' (from INI file)'
                stop
            endif
        end do
    !  print*, TimeLAT(1,1), TimeLAT(2,1)
    !  print*, ConcLAT(1,1), ConcLAT(2,1)
    !  print*, ConcLAT(1,5), ConcLAT(2,5)	
    !  print*, FlowLAT(1,1), FlowLAT(2,1)
    !  print*, FlowLAT(1,5), FlowLAT(2,5)
    !
    !pause 33  
    endif
    !! Lateral flow end
      
      ! read Dis manning table
    open(unit = LUN, file = trim(reach_qm))
    do i=1,100
        read(LUN,*,end=300)Q_Sk_Table(1,i), Q_Sk_Table(2,i)
    enddo
300 close(LUN)
    Q_sk_tableEntry = i-1
    !print*, 'Q_sk_tableEntry=',Q_sk_tableEntry
    
      
      
      
      
    !  pause 3  
      
      
    !!!!!!---------- END:: READ THE UPSTREAM CONCENTRATION BOUNDARY DATA - Concentrations are imposed
    
    !-------------------------------------------------------------------------------
    !----------  END:: READ THE INPUT DATA
    !-------------------------------------------------------------------------------
    
    ! Need to make sure that input data exist at each timestep, in other words
    ! we need to interpolate in time if the input data times are not the same
    ! as the simulation times
    do n = 1, nt
        time(n) = Tstart + (n - 1) * dt
    end do
    
    !--------------------
    ! Interpolate using splines the concentration BCs to the simulation times
    ! We do this outside the main time loop because we need some BC values
    ! before entering the time loop
    allocate(Bcoef(ntcb))
    allocate(Ccoef(ntcb))
    allocate(Dcoef(ntcb))
    allocate(tmp1D(nt))
    allocate(tmp2D(nt, np))
    
    write(*, *) 'Applying spline interpolation in time on the BC sediment concentrations'
    do i = 1, np
        Bcoef = 0.0_fp; Ccoef = 0.0_fp; Dcoef = 0.0_fp
        
        call spline(ntcb, TimeBND(:, i), ConcBND(:, i), Bcoef, Ccoef, Dcoef)
        
        do n = 1, nt
            fval = ispline(time(n), ntcb, TimeBND(:, i), ConcBND(:, i), Bcoef, Ccoef, Dcoef)
            tmp1D(n) = fval
        end do
        tmp2D(:, i) = tmp1D
    end do
    
    deallocate(ConcBND); allocate(ConcBND(nt, np))
    ConcBND = tmp2D
    
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Dcoef)
    deallocate(tmp1D)
    deallocate(tmp2D)
    deallocate(TimeBND)
    !--------------------
    
    !--------------------
    ! Interpolate using splines the concentration Downstream BCs to the simulation times
    ! We do this outside the main time loop because we need some BC values
    ! before entering the time loop
    allocate(Bcoef(ntdb))
    allocate(Ccoef(ntdb))
    allocate(Dcoef(ntdb))
    allocate(tmp1D(nt))
    allocate(tmp2D(nt, np))
    
    write(*, *) 'Applying spline interpolation in time on the DWONSTREAM BC sediment concentrations'
    do i = 1, np
        Bcoef = 0.0_fp; Ccoef = 0.0_fp; Dcoef = 0.0_fp
        
        call spline(ntdb, TimeDOW(:, i), ConcDOW(:, i), Bcoef, Ccoef, Dcoef)
        
        do n = 1, nt
            fval = ispline(time(n), ntdb, TimeDOW(:, i), ConcDOW(:, i), Bcoef, Ccoef, Dcoef)
            tmp1D(n) = fval
        end do
        tmp2D(:, i) = tmp1D
    end do
    
    deallocate(ConcDOW); allocate(ConcDOW(nt, np))
    ConcDOW = tmp2D
    
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Dcoef)
    deallocate(tmp1D)
    deallocate(tmp2D)
    deallocate(TimeDOW)
    !--------------------
    
    !! Lateral flow input
    if (Nlat > 0)then
    ! First flow
    allocate(Bcoef(ntlat))
    allocate(Ccoef(ntlat))
    allocate(Dcoef(ntlat))
    allocate(tmp1D(nt))
    allocate(tmp2D(nt, nplat))
    
    write(*, *) 'Applying spline interpolation in time on the lateral flow Q'
    do i = 1, nplat
        Bcoef = 0.0_fp; Ccoef = 0.0_fp; Dcoef = 0.0_fp
        
        call spline(ntlat, TimeLAT(:, i), FlowLAT(:, i), Bcoef, Ccoef, Dcoef)
        
        do n = 1, nt
            fval = ispline(time(n), ntlat, TimeLAT(:, i), FlowLAT(:, i), Bcoef, Ccoef, Dcoef)
            tmp1D(n) = fval
        end do
        tmp2D(:, i) = tmp1D
    end do
    
    deallocate(FlowLAT); allocate(FlowLAT(nt, nplat))
    FlowLAT = tmp2D
    
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Dcoef)
    deallocate(tmp1D)
    deallocate(tmp2D)
    
    ! Second Conc
    allocate(Bcoef(ntlat))
    allocate(Ccoef(ntlat))
    allocate(Dcoef(ntlat))
    allocate(tmp1D(nt))
    allocate(tmp2D(nt, nplat))
    
    write(*, *) 'Applying spline interpolation in time on the lateral flow Conc'
    do i = 1, nplat
        Bcoef = 0.0_fp; Ccoef = 0.0_fp; Dcoef = 0.0_fp
    
        call spline(ntlat, TimeLAT(:, i), ConcLAT(:, i), Bcoef, Ccoef, Dcoef)
    
        do n = 1, nt
            fval = ispline(time(n), ntlat, TimeLAT(:, i), ConcLAT(:, i), Bcoef, Ccoef, Dcoef)
            tmp1D(n) = fval
        end do
        tmp2D(:, i) = tmp1D
    end do
    
    deallocate(ConcLAT); allocate(ConcLAT(nt, nplat))
    ConcLAT = tmp2D
    
    deallocate(Bcoef)
    deallocate(Ccoef)
    deallocate(Dcoef)
    deallocate(tmp1D)
    deallocate(tmp2D)
    deallocate(TimeLAT)
    
    endif
    !! Lateral flow input ends
    
    !pause 11
    
    
    !----------  INITIALIZE THE CONC. FIELD
    ! CN is always the concentration at the previous time
    ! This will be improved later to accomondate realistic physical conditions
    ! Initialize the field with some constant value
    CN = sand_init
    CN(1, :) = ConcBND(1, :)
    
    
    !--- Calculate the spline coefficients outside the main time loop to speed things up
    
    !----------  INITITALIZE BED SEDIMENT THICKNESS for each sediment class/cross-section
    !----------   Once a cross-section is eroded by this much, only newly deposited sediment can be resuspended
    bed_t = 0.05    ! this currently sets initial bed thicknesses to 5 cm for each sediment particle class at each cross-section
    
    
    allocate(sumConc(nx))
    !open(ioutfile,file=trim(out_dir)//'hydro_input.dat')
    open(ioutfile,file=trim(out_dir)//'SedConSand_output.dat')
    write(ioutfile,30)time(1)/60, ((CN(i,j),i=1,nx),j=1,np)
30  format(<nx*np+1>f12.2)

end subroutine init_SAND_R07


subroutine cal_SAND_R07(n, npr, ioutfile, nlatt, inp_depth, inp_area, inp_flow, out_sand, Q_lat_from_ICM, SAND_lat_from_ICM, SAND_upstream_from_ICM, SAND_terminal_from_ICM)
    use precision
    !use params
    use mod_arrays_SAND_R07

    use sed_properties_SAND_R07
    use spline_func

    implicit none
    integer, intent(in) :: n, npr, ioutfile, nlatt
    real(kind=4), intent(in) :: SAND_upstream_from_ICM, SAND_terminal_from_ICM
    real(kind=4), dimension(npr), intent(in) :: inp_depth, inp_area, inp_flow
    real(kind=4), dimension(nlatt), intent(in) :: Q_lat_from_ICM, SAND_lat_from_ICM
    real(kind=4), dimension(npr), intent(out) :: out_sand
    integer :: nx, nt, np

    integer :: i, j, k, ilat, jlat, klat

    real(fp) :: fval, fval1, fval2  !modu_user
    integer :: ival, ival1, ival2

    real(fp) :: SedMassBed
    real(fp) :: SedMassSusp
    real(fp) :: SedMassEroded
    real(fp) :: SedMassDeposit
    real(fp) :: SedDepthEroded
    real(fp) :: SedDepthDeposit
    
    !integer :: cntUserDt, flgUserDt

    !character(len=256) :: reach_dist, reach_flow, reach_cbnd, reach_depth, reach_area,  reach_qm  !out_dir,
    ! Downstream  bnd file /lateral flow
    !character(len=256) :: reach_dbnd, reach_lat
    
    
    !character(len=128) :: outSTR

    !integer :: LUN, ioerr
  
    nt = NDt
    nx = NDx
    np = nclass
  
    if(n.eq.1 .or. modulo(time(n+1), DtUser*24).eq.0. .or. n.eq.(nt-1))write(*,*),'R07_SAND Nstep =', n, 'Days = ', int(time(n+1)/3600./24.), real(n)/real(nt-1)*100., '% completed'

      ! Determine when we will store the data (not every timestep)
      ! We compare the simulation time "tt" with the supplied DtUser using
      ! the module function; if modulo(time(n), DtUser) is zero we are exactly on
      ! the required time, otherwise if modulo(time(n), DtUser) < Dt then we are close
      ! to the desired write time.
      !flgUserDt = 0
      !modu_user = abs(modulo(time(n), DtUser))
      !if ( (comparereal(modu_user, 0.0_fp) == 0) .or. (modu_user < Dt) ) then
      !  cntUserDt = cntUserDt + 1
      !  flgUserDt = 1
      !endif
  
    depth=inp_depth
    area=inp_area
    flow=inp_flow
    
    do i = 1, nx
    ! Interpolate the time/space varying water depth, cross-section area and flow rate
    ! using cubic splines (the coefficients are calculated above)
    !  depth(i) = ispline(time(n), ntfl, TimeFL(:, i), DepthFL(:, i), Bdp(:, i), Cdp(:, i), Ddp(:, i))
    !  area(i)  = ispline(time(n), ntfl, TimeFL(:, i), AreaFL(:, i), Bar(:, i), Car(:, i), Dar(:, i))
    !  flow(i)  = ispline(time(n), ntfl, TimeFL(:, i), FlowFL(:, i), Bfl(:, i), Cfl(:, i), Dfl(:, i))
        width(i) = area(i)/depth(i)
        ! Calculate the total suspended sediment concentration
        sumConc(i) = sum(CN(i, :))
    end do

20  format(<3*nx+1>f12.2)
30  format(<nx*np+1>f12.2)    
  !if ( flgUserDt > 0 ) then
  !! These are used to calculate the sediment loads outside the time loop
  !! Most likely need these calculations to be done within the main time loop
  !  flowUser(cntUserDt, :)  = flow
  !  areaUser(cntUserDt, :)  = area
  !  depthUser(cntUserDt, :) = depth
  !  timeUser(cntUserDt)     = time(n)

  !  ! Store the simulated data into the array "C" if flgUserDt > 0
  !  do k = 1, np
  !    C(cntUserDt, 1:nx, k) = CN(:, k)
  !  end do
  !endif

  ! Below we are calculating CN1 (concentration at next time step) using
  ! the discretized equation. At timestep nt -1 CN1 (concentration at nt)
  ! has been calculated and its value(s) have been stored into the CN and
  ! subsequently into the C(:,:,:) array (see above)
    if (n /= nt) then
        do k = 1, np            ! loop over sediment partical classes
            do i = 2, nx - 1    ! loop over cross-sections
                ! Calculate the coefficients of the equation first
                ! This assumes uneven spatial spacing and central differences.

                ! This is to use the 2nd order central differencing for the differential equation
                ! for unequally spaced grid points or cross-sections
                dx  = xx(i+1) - xx(i)
                dx1 = xx(i) - xx(i-1)

                ! This is the lamda parameter (see notes)
                lmd1 = dx1 / dx

                !!! The coefficients of the f' discretizetion (see notes)
                denom = lmd1 * (1.0 + lmd1) * dx

                a0 = (lmd1 ** 2.0) / denom
                a1 = (1.0 - lmd1 ** 2.0) / denom
                a2 = -1.0 / denom

                !!! The coefficients of the f'' discretizetion (see notes)
                denom = lmd1 * (1.0 + lmd1) * (dx ** 2.0)

                b0 = 2.0 * lmd1 / denom
                b1 = -2.0 * (1.0 + lmd1) / denom
                b2 = 2.0 / denom

                mu = a0 * area(i+1) + a1 * area(i) + a2 * area(i-1)
                mu = ( dt / area(i) ) * ( flow(i) - ks * mu )

                xi = ks * dt

                aa = - a0 * mu + b0 * xi
                bb = 1.0 - a1 * mu + b1 * xi
                cc = - a2 * mu + b2 * xi

                arni = area(i)
                SrcTerm = 0.0
                if( useSRC ) then
                    if (D50(k) < 0.000062) then
                      SrcTerm = srcchsv(D50(k), depth(i), width(i), flow(i), flow(i) / area(i), CN(i, k), sumConc(i), Tcrit(k)) !! Change Nazmul
                    else
                      SrcTerm = srcsand(D50(k), depth(i), width(i), flow(i) / area(i), CN(i, k), sumConc(i)) !! Change Nazmul
                    endif
                    ! SrcTerm from srchsv and srcsand is sediment flux per unit length - units: g/m/s
                    ! bed_t(i,k) = thickness of bed sediments (per particle class, k) in cross-section, i [m]
                    ! bed_bd = bulk density of bed sediments [g/m^3]                    
                    SedMassBed      = 0.0
                    SedMassSusp     = 0.0
                    SedMassEroded   = 0.0
                    SedMassDeposit  = 0.0
                    SedDepthEroded  = 0.0
                    SedDepthDeposit = 0.0
                    
                    if (SrcTerm > 0.0) then                                     ! if erosive, check that there is enough sediment on the bed to be eroded
                        SedMassBed = bed_t(i,k)*bed_bd*width(i)                 ! [g/m] mass of sediment particle class, k, in bed per unit length
                        SedMassEroded = SrcTerm*dt                              ! [g/m] mass of sediment eroded during timestep, per unit length
                        SedDepthEroded = SedMassEroded/(bed_bd*width(i))        ! [g/m]*[m3/g]*[1/m] = [m]    
                        if ( SedMassEroded >= SedMassBed ) then                
                            SrcTerm = SedMassBed/dt                             ! limit eroded mass flux to be no more than the mass of sediment in the bed
                            bed_t(i,k) = 0.0                                    ! zero out bed depth since more mass can erode than is available
                        else
                            bed_t(i,k) = bed_t(i,k) - SedDepthEroded            ! update bed depth for eroded mass
                        endif
                    else                                                        ! if depositional, update bed depth
                        SedMassSusp = CN(i,k)*area(i)                           ! [g/m] mass of sediment suspended in water column, per unit length
                        SedMassDeposit = abs(SrcTerm*dt)                        ! [g/m] mass of sediment eroded during timestep, per unit length
                        if ( SedMassDeposit >= SedMassSusp) then
                            SrcTerm = -1.0*SedMassSusp/dt                       ! limit deposited mass flux to be no more than the mass of sediment currently suspended
                            SedDepthDeposit = SedMassSusp/(bed_bd*width(i))     ! [g/m]*[m3/g]*[1/m] = [m]
                        else
                            SedDepthDeposit = SedMassDeposit/(bed_bd*width(i))  ! [g/m]*[m3/g]*[1/m] = [m]
                        endif
                        bed_t(i,k) = bed_t(i,k) + SedDepthDeposit
                    endif
                    SrcTerm = (dt / arni) * SrcTerm                             ! convert from g/s/m to g/m^3
                endif


                ! Add lateral input
                do j=1, Nlat
                    ilat=Idlat(j)
                    jlat=Nsclat(j)
                    klat=Ntplat(j)
                    if (i >= ilat .AND. i < (ilat+jlat))then
                        if (FlowLAT(n+1,j) > 0. .and. klat .eq. 1)then
                            SrcTerm = SrcTerm + FlowLAT(n+1,j)/real(jlat)*(ConcLAT(n+1,j)-CN(i,k))/(abs(flow(i))+FlowLAT(n+1,j)/real(jlat))
                        endif
                        !coupling from ICM
                        if (Q_lat_from_ICM(j) .gt. 0. .and. SAND_lat_from_ICM(j) .ge. 0. .and. klat .eq. 3)then
                            SrcTerm = SrcTerm + Q_lat_from_ICM(j)/real(jlat)*(SAND_lat_from_ICM(j)-CN(i,k))/(abs(flow(i))+Q_lat_from_ICM(j)/real(jlat))
                        endif
                        
                        if (klat.eq.2) then
                            print*, 'Rating curves for lateral SAND input have not been supported yet.'
                            stop
                        endif
                    endif
                enddo       

                !CN1(i, k) = aa * CN(i+1, k) + bb * CN(i, k) + cc * CN(i-1, k) + (dt / arni) * SrcTerm
                CN1(i, k) = aa * CN(i+1, k) + bb * CN(i, k) + cc * CN(i-1, k) + SrcTerm
                if( CN1(i, k) <= smallC ) CN1(i, k) = 0.0
            
            end do  ! end loop over cross-sections

            !!!!! Apply the boundary conditions -> (upstream: i = 1), (downstream: i = nx)
            !!!!! For upstream boundary use fixed concentration BC
            CN1(1, k) = ConcBND(n+1, k)
            if(SAND_upstream_from_ICM .ge. 0. )CN1(1, k)=SAND_upstream_from_ICM   !coupling from ICM
            if (flow(i) < 0.)then
                CN1(1, k)=CN1(2, k)
            endif
 
            !!!!! For downstream boundary use downwind flux BC or downwind differencing
            if (flow(i) > 0. )then
                if (bc_option == 1) then
                    ! This is for the downstream BC: Q*C - Ks*A*dC/dx = 0
                    ! Here we consider the ghost point at i = m+1 and second order
                    ! central differences
                    dx1 = xx(nx) - xx(nx-1)
                    dx = dx1

                    ! This is the lamda parameter (see notes)
                    lmd1 = dx1 / dx

                    !!! The coefficients of the f' discretizetion (see notes)
                    denom = lmd1 * (1.0 + lmd1) * dx

                    a0 = (lmd1 ** 2.0) / denom
                    a1 = (1.0 - lmd1 ** 2.0) / denom
                    a2 = -1.0 / denom

                    !!! The coefficients of the f'' discretizetion (see notes)
                    denom = lmd1 * (1.0 + lmd1) * (dx ** 2.0)

                    b0 = 2.0 * lmd1 / denom
                    b1 = -2.0 * (1.0 + lmd1) / denom
                    b2 = 2.0 / denom

                    ! Here we assume that area(nx+1) = area(nx)
                    mu = a0 * area(nx) + a1 * area(nx) + a2 * area(nx-1)
                    mu = ( dt / area(nx) ) * ( flow(nx) - ks * mu )

                    xi = ks * dt

                    aa = - a0 * mu + b0 * xi
                    bb = 1.0 - a1 * mu + b1 * xi
                    cc = - a2 * mu + b2 * xi

                    dd = aa * (flow(nx) / (a0 * ks * area(nx)) - a1 / a0) + bb
                    ee = - aa * (a2 / a0) + cc

                    arni = area(nx)
                    SrcTerm = 0.0
                    if( useSRC ) then
                        if (D50(k) < 0.000062) then
                            SrcTerm = srcchsv(D50(k), depth(nx), width(nx), flow(nx), flow(nx) / area(nx), CN(nx, k), sumConc(nx), Tcrit(k)) !! Change Nazmul
                        else
                            SrcTerm = srcsand(D50(k), depth(nx), width(nx), flow(nx) / area(nx), CN(nx, k), sumConc(nx)) !! Change Nazmul
                        endif
                    endif
                    CN1(nx, k) = dd * CN(nx, k) + ee * CN(nx - 1, k) + (dt / arni) * SrcTerm
                    if( CN1(nx, k) <= smallC ) CN1(nx, k) = 0.0
                
                elseif (bc_option == 2) then
                    ! This is for the downstream BC: Ks*A*dC/dx = 0, that is C(nx+1, k) = C(nx, k)
                    ! Here we consider the ghost point at i = m+1 and second order
                    ! central differences
                    dx1 = xx(nx) - xx(nx-1)
                    dx = dx1

                    ! This is the lamda parameter (see notes)
                    lmd1 = dx1 / dx

                    !!! The coefficients of the f' discretizetion (see notes)
                    denom = lmd1 * (1.0 + lmd1) * dx

                    a0 = (lmd1 ** 2.0) / denom
                    a1 = (1.0 - lmd1 ** 2.0) / denom
                    a2 = -1.0 / denom

                    !!! The coefficients of the f'' discretizetion (see notes)
                    denom = lmd1 * (1.0 + lmd1) * (dx ** 2.0)

                    b0 = 2.0 * lmd1 / denom
                    b1 = -2.0 * (1.0 + lmd1) / denom
                    b2 = 2.0 / denom

                    ! Here we assume that area(nx+1) = area(nx)
                    mu = a0 * area(nx) + a1 * area(nx) + a2 * area(nx-1)
                    mu = ( dt / area(nx) ) * ( flow(nx) - ks * mu )

                    xi = ks * dt

                    aa = - a0 * mu + b0 * xi
                    bb = 1.0 - a1 * mu + b1 * xi
                    cc = - a2 * mu + b2 * xi

                    arni = area(nx)
                    SrcTerm = 0.0
                    if( useSRC ) then
                        if (D50(k) < 0.000062) then
                            SrcTerm = srcchsv(D50(k), depth(nx), width(nx), flow(nx),flow(nx) / area(nx), CN(nx, k), sumConc(nx), Tcrit(k)) !! Change Nazmul
                        else
                            SrcTerm = srcsand(D50(k), depth(nx), width(nx), flow(nx) / area(nx), CN(nx, k), sumConc(nx)) !! Change Nazmul
                        endif
                    endif
                    CN1(nx, k) = (aa + bb) * CN(nx, k) + cc * CN(nx-1, k) + (dt / arni) * SrcTerm
                    if( CN1(nx, k) <= smallC ) CN1(nx, k) = 0.0
                
                elseif (bc_option == 3) then
                    ! This is to use the 2nd order upwind differencing for the downstream BC
                    ! for unequally spaced grid points or cross-sections
                    dx  = xx(nx) - xx(nx-1)
                    dx1 = xx(nx-1) - xx(nx-2)
                    dx2 = xx(nx-2) - xx(nx-3)

                    ! These are the lamda1 and lamda2 in the notes
                    lmd1 = dx1 / dx
                    lmd2 = dx2 / dx

                    !!! The coefficients of the f' discretizetion (see notes)
                    denom = lmd1 * (1.0 + lmd1) * dx

                    a0 = lmd1 * (2.0 + lmd1) / denom
                    a1 = - ((1.0 + lmd1) ** 2.0) / denom
                    a2 = 1.0 / denom

                    !!! The coefficients of the f'' discretizetion (see notes)
                    denom = lmd1 * lmd2 * (lmd1 + lmd2) * ((1.0 + lmd1) ** 3.0)
                    denom = denom * (1.0 + lmd1 + lmd2) * (dx ** 2.0)

                    b0 = lmd2 * (1.0 + lmd1 + lmd2) * (2.0 + 2.0 * lmd1 + lmd2) * ((1.0 + lmd1) ** 3.0 - 1.0)
                    b0 = 2.0 * ( b0 - lmd1 * (2.0 + lmd1) * ((1.0 + lmd1 + lmd2) ** 3.0 - (1.0 + lmd1) ** 3.0) )
                    b0 = b0 / denom

                    b1 = -2.0 * lmd2 * ((1.0 + lmd1) ** 3.0) * (1.0 + lmd1 + lmd2) * (2.0 + 2.0 * lmd1 + lmd2)
                    b1 = b1 / denom

                    b2 = lmd2 * (1.0 + lmd1 + lmd2) * (2.0 + 2.0 * lmd1 + lmd2)
                    b2 = 2.0 * ( b2 + lmd1 * (2.0 + lmd1) * ((1.0 + lmd1 + lmd2) ** 3.0) )
                    b2 = b2 / denom

                    b3 = -2.0 * lmd1 * (2.0 + lmd1) * (1.0 + lmd1) ** 3.0
                    b3 = b3 / denom

                  ! The coefficients of the discretized differential equation
                    mu = a0 * area(nx) + a1 * area(nx-1) + a2 * area(nx-2)
                    mu = (dt / area(nx)) * (flow(nx) - ks * mu)

                    xi = ks *dt

                    aa = 1.0 - a0 * mu + b0 * xi
                    bb = - a1 * mu + b1 * xi
                    cc = - a2 * mu + b2 * xi
                    dd = b3 * xi

                    arni = area(nx)
                    SrcTerm = 0.0
                    if( useSRC ) then
                        if (D50(k) < 0.000062) then
                            SrcTerm = srcchsv(D50(k), depth(nx), width(nx), flow(nx),flow(nx) / area(nx), CN(nx, k), sumConc(nx), Tcrit(k)) !! Change Nazmul
                        else
                            SrcTerm = srcsand(D50(k), depth(nx), width(nx), flow(nx) / area(nx), CN(nx, k), sumConc(nx)) !! Change Nazmul
                            !write(*, *) time(n), arni, SrcTerm
                        endif
                    endif
                    CN1(nx, k) = aa  * CN(nx, k) + bb * CN(nx-1, k) + cc * CN(nx - 2, k) + dd * CN(nx - 3, k) + (dt / arni) * SrcTerm
                    if( CN1(nx, k) <= smallC ) CN1(nx, k) = 0.0
                else
                    write(*,*) 'Wrong value for bc_option: bc_option = 1,2 or 3'
                    stop
                endif
      
            else
                CN1(nx, k)=ConcDOW(n+1, k)
                if(SAND_terminal_from_ICM .ge. 0. )CN1(nx, k)=SAND_terminal_from_ICM  !coupling from ICM
            endif
    
            CN = CN1

        end do ! End of nclass loop

        out_sand(:)=CN(:,1)
        if(modulo(time(n+1), DtUser).eq.0. .or. n.eq.(nt-1))then
            !write(ioutfile,20)time(n), (depth(i),i=1, nx), (area(i),i=1, nx),(flow(i),i=1, nx)
            write(ioutfile,30)time(n+1)/60, ((CN(i,j),i=1,nx),j=1,np)
        endif

    endif ! End of n /= nt

end subroutine cal_SAND_R07
