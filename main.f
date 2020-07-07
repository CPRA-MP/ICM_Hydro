
! This section contains the Doxygen comments for
!> @mainpage
!> @section Introduction
!> This is the Fortran code for the Coastal Master Plan Ecohydrology routines,
!> originally prepared by Dr. J. Alex McCorquodale at the University of New Orleans.
!> This code has been edited and updated by Eric White at The Water Institute of the Gulf.
!> @section Compiler
!> This source code was originally tested and compiled using Intel's Fortran compiler.

!BUG! EVENTUALLY ADD THIS PORTION TO THE DOXYGEN COMMENTS

! @section Participating Organizations
! @image html CPRAlogo.png
! <a href="http:\\www.coastal.LA.gov">Coastal Protection and Restoration Authority</a>
! @image html WIlogo.png
! @image html UNOlogo.png


!> @file
!> @brief This is the main program/control subroutine of the Ecohydrology model.
!> @details This is the control portion of the Ecohydrology model which controls the model I/O,
!> determines the simulation timesteps, assigns initial conditions for each timestep,
!> steps through each simulation timestep, and eventually post-processes model output to be read into
!> other ICM routines.


! Ecohydro program is set up solely to call the 'main' subroutine. The only reason 'main' is a subroutine the program is to incorporate the
! Doxygen documentation logic structure, which needs to be embedded within a subroutine, not a program.
      program Ecohydro

	!ZW 3/13/2015 move the log file open statement to here
      !open (unit=1, file='hydro_run.log', status='unknown')

      call main

      write(1,*) '----------------------------------------------------'
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) ' HYDRO MODEL RUN IS COMPLETE.'
      write(1,*) '----------------------------------------------------'

      write(*,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) ' HYDRO MODEL RUN IS COMPLETE.'
      write(*,*) '----------------------------------------------------'


!      pause
      stop
      end

!> @author J. Alex McCorquodale - University of New Orleans
!> @author Eric White - The Water Institute of the Gulf

ccccccccccccccccccLower Breton/ Lower Barataria 10/31/13
c
c         Main  MASS BALANCE MODEL FOR Lake Pontchartrain WATER QUALITY
c
c  NOTE: All distance units input as km are converted to m in code
c
c      Elements
c      i - M Internal links
c      j - N Storage cells
c      ib - Mb Boundary Links
c      +ve direction into the cell
c
c  NOTATION: follows EPA SWMM 3 and 4
c* Link no., u/s Cell no., d/s Cell No., Length, Width, Deptho, n, Km, Ken, Kx, Kb
c* Cell no., As, Eso, Bed, Ahydro,dAdz, CSS, S
c* Boundary link no., type(+1 u/s,-1 d/s), Cell no.,itype, Length, Width,depth, n, Km, Ken, Kx, Kb
c        icc(j,10) link numbers connected to cell (- in, + out)
c* Cell no., link in -, link out +                      !up to 6 numbers fo rthe links connected to a cell
c*                                                                        !connectivity
c* u/s BC Daily flows --> File=TimeSeries1.dat FORMAT free --> Day,flow m3/s, SS mg/L, S ppt
c* d/s BC Stage --> File=TimeSeries2.dat FORMAT free --> hours,stage m, SS mg/L, S ppt
c***********************************************************************************************
c* Matrix Variables:
c
c	Cells;-    j
!> @param[out] As(N,3)				surface area of cells (m2)
!> @param[out] Eso(N)					initial stage of storage cells (m)
!> @param[out] Bed(N)					bed elevation of storage cells (m)
!> @param[out] ds(N,3)				depth in storage cells         (m)
!> @param[out] Es(N,3)				stage in storage cells         (m)
!> @param[out] Ahydro(N)				area of hydrologically connected area for water balance (m2)
!> @param[out] ow						% open water in attach drainage basin
!> @param[out] dAdz(N)				change in As with depth (m)
!> @param[out] CSS(N,3)				sediment concentration	(mg/L)
!> @param[out] S(N,3)					salinity (TDS)	(ppt)
!> @param[out] acss(N)				resuspension parameter
!> @param[out] bcss(N)				resuspension parameter
!> @param[out] Vss(N)					deposition velocity		(m/d)
!> @param[out] tau(N,3)				bed shear stress		(Pa)
!> @param[out] Tauc(N)				critical shear stress	(Pa)
!> @param[out] Fetch(N,10)			Fetch EW				(km)
!> @param[out] Fetch(N,10)			Fetch NS				(km)

c	Cellsh - Properties for marsh compontent of cells
!> @param[out] Ahf(N)					Marsh area subject to flooding
!> @param[out] flood(N)				Fraction of Ahydro subject to flooding
!> @param[out] Eho(N)					Bed marsh elevation in flood prone area	(m)
!> @param[out] Eh(N,3)				stage in flooded area of marsh.
!> @param[out] CSSh(N,3)				TSS in flood area of marsh
!> @param[out] Sh(N,3)				Salinity in flood area of marsh			(ppt)
!> @param[out] acssh(N)				resupension in flooded area of marsh
!> @param[out] Vsh(N)					settling velocity in flooded area of marsh (m/d)
!> @param[out] Esho(N)				dead storage in marshes					(m)
c
cc***********************************************************************************************
c        Links:-   i
!> @param[out] Deptho(M)				initial channel depth	(m)
!> @param[out] Depth(M,3)				channel depth			(m)
!> @param[out] Length(M)				channel length			(m)
!> @param[out] Width(M)				channel width			(m)
!> @param[out] an(M)					channel Mannings n
!> @param[out] Km(M)					minor loss coefficient due to structures based on Km(Vch)^2/2g
!> @param[out] Ken(100)				entrance loss coefficient
!> @param[out] Kx(M)					exit loss coefficient based on Kx(Vch-Vrec)^2/2g
!> @param[out] Exy(M)					diffusion coefficient thru link between cells (m2/s)
!> @param[out] Q(M,3)					channel flow in m3/s
!> @param[out] Qss(M,3)				channel sediment flow	(kg/s)
!> @param[out] QSal(M,3)				channel TDS transport	(kg/s)
!> @param[out] jus(M)					u/s junction no.
!> @param[out] jds(M)					d/s junction no.
!> @param[out] itype(M)				"-1" is u/s BC, "+1" is d/s BC, "2" is two cell link.


c***********************************************************************************************
c        Forcing:-
!> @param[out] Qtrib(26,366*25)		Mean daily tributary flow			(m3/s)
!> @param[out] CssT(26,366*25)		Mean daily tributary sediment conc	(mg/L)
!> @param[out] Strib(30,366*25)		Mean daily tributary salinity		(ppt)
!> @param[out] Tcoef(M,25)			Tide generation coefficients
!> @param[out] Sal(1,1)				d/s BC on Salinity					(ppt)
!> @param[out] BCTSS(10)				d/s BC on TSS						(mg/L)	!added TSS=CSS+Algae+Detritus XX!!NOT USED JAM July 09
!> @param[out] PET(366)				daily potential evapotranspiration	(m/d)
!> @param[out] Precip(366)			daily rainfall						(m/d)
!> @param[out] QCHEM(N,15,366*25)     chemical load in tributaries (g/sec)
!> @param[out] CHEM(N,15,2)     water quality concentration
!> @param[out] QChemSUM(15)
!> @param[out] DChem(15)
!> @param[out] CSSo(N)
!> @param[out] QMult(N,40)            Multiplier on tributary-compartment lookup matrix
!> @param[out] PET(days,ETgages)
!> @param[out] Rain(days,raingages)
!> @param[out] Chem(N,15,2)
!> @param[out] Qmultdiv(N,40)         Multiplier on diversion-compartment lookup matrix
!> @param[out] QChemdiv(26,15,366*25) chemical load in diversions (g/sec)
!> @param[out] cChemdiv(1,15,366*25)  WQ concentrations in diversion input data (mg/L)
!> @param[out] QAtm(1,15,366*25)
!> @param[out] TSS(N)
!> @param[out] cChemface(15)
!> @param[out] T7(366*25)
!> @param[out] TKN(M,2)
!> @param[out] uplandNP(M,14)
c
c***********************************************************************************************
c        Control parameters:-
!> @param[out] N						number of compartments
!> @param[out] M						number of links
!> @param[out] Mus					number of u/s BCs
!> @param[out] Mds					number of d/s BCs
!> @param[out] dt						computational time step			(s)
!> @param[out] dTprint				output time step			(hours)
!> @param[out] startrun				Julian day of start of run
!> @param[out] endrun 				Julian day at end of run
!> @param[out] Specg					Specific gravity of sediment
!> @param[out] time                   elapsed model time seconds
!> @param[out] day                    elapsed model time in decimals days
!> @param[out] thour                  elapsed model time in decimal hours
!> @param[out] tmon                   elapsed model time in decimal months
!> @param[out] kthr                   integer of hourly elapsed model time +1
!> @param[out] kmon				    integer of monthly elapsed model time +1
!> @param[out] kday                   integer of daily elapsed model time +1
!> @param[out] dday                   decimal portion of day, dday=0.0 at 0:00 (midnight)
!> @param[out] hday                   decimal portion of day normalized to noon, hday = 0.0 at 12:00 (noon)

!> @param      ddayhr                 hour of day (dday/24.)
!> @param      windupdate             flag used to update wind data
!> @param      dtwind                 timestep of observed wind data file (hr)
!> @param      windrow                row of wind data array that corresponds to hour of simulation
c
c***********************************************************************************************
c	   Chemical Inputs
c        1 = NO3 + NO2
c        2 = NH4
c        3 = DIN (1+2)
c        4 = Organic N
c        5 = TIP
c        6 = TOC
c        7 = DO
c	   8 = Live Algae (ALG)
c	   9 = Dead Algae (DET)
c        10= DON
c        11= DOP
c        12= DIP   !PArtition SRP = ParP*TP
c        13= ChLa  !PArtition Chla= ParCla*LivA
c        14= POP !-EDW used to say 14=TKN !-EDW
c***********************************************************************************************

! Main subroutine is the main 'control' subroutine for the ecohydro model.
! Embedded Doxygen comments only work within subroutines, not 'programs'.
      subroutine main

      use params

      integer :: i,it,j,jj,jjj,sedclass       !iterators used in main.f

      Character*100 header

!>> *********************************Adding 1d start
      use constants_module
      use arrays_module
      use arrays_section_module
      use var_module
      use matrix_module
      use sgate_module
      use xsec_attribute_module

      integer(kind=4) :: i_1d, n_1d, ntim, igate, pp_1d, &
						boundaryFileMaxEntry, ppp, qqq, saveFrequency, noLatFlow !, areanew
      real(kind=4) :: cour, da, dq, dxini, yy, x, thetas, thesinv, tfin
      real(kind=4) :: t_1d, skk, qq, qn, xt, r_interpol, t1, t2, maxCourant
      real(kind=4) :: arean, areac, hyrdn, hyrdc, perimn, perimc, qcrit, &
					s0ds, timesDepth, latFlowValue

      character(len=128) :: upstream_path , downstream_path
      character(len=128) :: manning_strickler_path, output_path, other_input
      character(len=128) :: channel_width_path, lateralFlow_path
      character(len=128) :: path

    ! open file for input data
      character(len=128) :: dx_path

      write(*,*) '1D Code is running.'
      call cpu_time( t1 )

      open(unit=601,file="../lower_Mississippi/input/input_BC_2_SWP_2018.txt", &
       status='unknown')

    ! read data
      read(601,*) dtini
      read(601,*) dxini
      read(601,*) tfin
      read(601,*) ncomp
      read(601,*) ntim
      read(601,*) phi
      read(601,*) theta
      read(601,*) thetas
      read(601,*) thesinv
      read(601,*) alfa2
      read(601,*) alfa4
      read(601,*) f
      read(601,*) skk
      read(601,*) yy
      read(601,*) qq
      read(601,*) cfl
      read(601,*) ots
      read(601,*) yw
      read(601,*) bw
      read(601,*) w
      read(601,*) option
      read(601,*) yn
      read(601,*) qn
      read(601,*) igate
      read(601,*) xSection_path
      read(601,*) manning_strickler_path
      read(601,*) upstream_path
      read(601,*) downstream_path
      read(601,*) channel_width_path
      read(601,*) dx_path
      read(601,*) lateralFlow_path
      read(601,*) output_path
      read(601,*) option_dsbc
      read(601,*) maxTableLength
      read(601,*) nel
      read(601,*) timesDepth
      read(601,*) other_input
      read(601,*) boundaryFileMaxEntry
      read(601,*) saveFrequency
      read(601,*) noLatFlow             ! no of lateral flow coupling

      allocate(latFlowLocations(noLatFlow))   ! all the first nodes where a lateral flow starts
      allocate(latFlowType(noLatFlow))        ! Lateral flow type: Type = 1 = time series; Type 2 = flow as a function of upstream flow
      allocate(latFlowXsecs(noLatFlow))       ! no of x-secs at the downstream that the lateral flow is applied

      read(601,*) (latFlowLocations(i), i=1, noLatFlow)
      do i=1,noLatFlow
          if ((latFlowLocations(i)-1)*(latFlowLocations(i)-ncomp) .eq. 0) then
              print*, 'ERROR: Lateral flow cannot be applied at the boundaries'
              stop
          end if
      end do
      read(601,*) (latFlowType(i), i=1, noLatFlow)
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

      read(601,*) (latFlowXsecs(i), i=1, noLatFlow)

      close (601)

    ! Allocate arrays
      call setup_arrays(ntim, ncomp, maxTableLength, boundaryFileMaxEntry, noLatFlow)
      call setup_arrays_section
      call setup_xsec_attribute_module(nel, ncomp)

      dt_1d = dtini

      open(unit=690, file=trim(dx_path))
      do i=1,ncomp-1
          read(690, *) x, dx(i)
      end do
      close(690)

      ! reading Strickler's coefficient at each section
      open(unit=685,file=trim(manning_strickler_path), status='unknown') !! //'Mannings_Stricklers_coeff.txt'
      do i=1,ncomp
          read(685, *) x, sk(i)
          call readXsection(i,(1.0/sk(i)),timesDepth)
          ! This subroutine creates attribute table for each cross sections and saves in the hdd
          ! setting initial condition
          !y(1,i) = yy ! + z(i)
          oldY(i) = yy ! + z(i)
      end do
      close(685)

      do i=1,ncomp
          call create_I2(i,ncomp)
      end do

      ityp = 1

    ! setting initial condition
    ! setting initial condition from previous work
    !open(unit=91,file=trim(output_path)//'initialCondition.txt', status='unknown')
    ! read(91, *)
    !do i=1,ncomp
    !    read(91, *) oldQ(i), oldY(i)
    !end do
    !close(91)

    !q(1, :) = qq
      oldQ = qq

    ! reading Q-Strickler's coefficient multiplier table
      open(unit=686,file=trim(other_input) &
       //'Q_Mannings_table.txt', status='unknown')
      do i=1,maxTableLength
          read(686,*,end=300)Q_Sk_Table(1,i), Q_Sk_Table(2,i)
      enddo
300   close(686)
      Q_sk_tableEntry = i-1


      x = 0.0

    ! Read hydrograph input Upstream
      open(unit=687, file=upstream_path)
      do i=1,boundaryFileMaxEntry
          read(687,*,end=301) USBoundary(1,i), USBoundary(2, i)
      end do
301   close(687)
      ppp = i-1

    ! Read hydrograph input Downstream
      open(unit=688, file=downstream_path)
      do i=1,boundaryFileMaxEntry
        read(688,*,end=302)  DSBoundary(1, i), DSBoundary(2, i)
      end do
302   close(688)
      qqq = i-1

      t_1d = 0.0

    ! applying boundary
    ! interpolation of boundaries at the initial time step

      oldQ(1)    =r_interpol(USBoundary(1, :),USBoundary(2, :),ppp,t_1d)
      oldY(ncomp)=r_interpol(DSBoundary(1, :),DSBoundary(2, :),qqq,t_1d)

    ! DS Boundary treatment: from water level to area time series
      ncompElevTable = xsec_tab(1,:,ncomp)
      ncompAreaTable = xsec_tab(2,:,ncomp)

	!open(unit=81,file=trim(output_path)//'DS_area.txt', status='unknown')
      xt=oldY(ncomp)
      oldArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,xt)
    !write(81, *) t, oldArea(ncomp)
      lateralFlowTable = 0.0
      ! read lateral flow conditions
      do i=1,noLatFlow
        if ((latFlowType(i) .eq. 1) .or. (latFlowType(i) .eq. 2)) then
          write(file_num,'(i4.4)')latFlowLocations(i)
          open(689,file=trim(lateralFlow_path)//'lateral_'//file_num//'.txt')
          do n_1d=1,boundaryFileMaxEntry
              read(689,*,end=303) lateralFlowTable(1, n_1d, i), lateralFlowTable(2, n_1d, i)
          end do
303       close(689)
          dataInEachLatFlow(i) = n_1d-1
        end if
      end do

    ! Open files for output
      path = trim(output_path) // 'output_wl.txt'
      open(unit=608, file=trim(path), status='unknown')
      path = trim(output_path) // 'q.txt'
      open(unit=609, file=trim(path), status='unknown')

      path = trim(output_path) // 'area.txt'
      open(unit=651, file=trim(path), status='unknown')

    ! Output initial condition

      write(608, 10)  t_1d, (oldY(i), i=1,ncomp)
      write(609, 10)  t_1d, (oldQ(i), i=1, ncomp)
      write(651, 10) t_1d, (oldArea(i), i=1, ncomp)

!>> *********************************Adding 1d end

!>@par General Structure of Subroutine Logic:

! determine start time for calculating runtimes
      call SYSTEM_CLOCK(runtime_start,count_rate1,count_max1)
!>> Open file with model information regarding file formatting in ICM (written only before first model year - NOT updated yearly).
      open (unit=300,  file= 'ICM_info_into_EH.txt', status ='unknown')
!>> Read in names of Veg files to write output to.
      read(300,*) VegWaveAmpFile
      read(300,*) VegMeanSalFile
      read(300,*) VegSummerDepthFile
      read(300,*) VegSummerSalFile
      read(300,*) VegSummerTempFile
      read(300,*) VegTreeEstCondFile
      read(300,*) VegBIHeightFile
      read(300,*) VegPerLandFile
!>> Read in grid dimensions for 500m grid used by LAVegMod.
      read(300,*) n_500m_cells
      read(300,*) veg_matrix_rows
      read(300,*) veg_matrix_cols
      read(300,*) veg_xllcorner
      read(300,*) veg_yllcorner
!>> Read in grid dimensions for 1000m grid used by EwE.
      read(300,*) n_1000m_cells
      read(300,*) ewe_matrix_rows
      read(300,*) ewe_matrix_cols
      read(300,*) ewe_xllcorner
      read(300,*) ewe_yllcorner
!>> Read in name of log file to write output to
      read(300,*) Hydro_log
!      Hydro_log = 'ICM_NoTimeStamp_Hydro.log'
      close(300)

      open (unit=1, file=TRIM(ADJUSTL(Hydro_log)),position='append')

      write(1,*) '----------------------------------------------------'
      write(1,*) '------------- RUNNING HYDROLOGY MODEL --------------'
      write(1,*) '----------------------------------------------------'
      write(1,*)

      write(*,*) '----------------------------------------------------'
      write(*,*) '------------- RUNNING HYDROLOGY MODEL --------------'
      write(*,*) '----------------------------------------------------'
      write(*,*)

!>> Read RunControlR.dat and set simulation variables (RunControlR.dat is updated yearly by ICM.py).
      allocate(D50(4))	!sediment class parameters that need to be allocated for control data inputs (e.g. before allocate_params is called)
      allocate(Tcrit(4))	!sediment class parameters that need to be allocated for control data inputs (e.g. before allocate_params is called)

      open (unit=30, file= 'RuncontrolR.dat', status = 'unknown')
      READ(30,*) year          !  1        year of simulation (e.g. 2015)
      READ(30,*) Nyear         !  2        number of years simulated before exiting back to ICM
      READ(30,*) dt            !  3        Computational time step
      READ(30,*) M             !  4        # of links
      READ(30,*) N             !  5        # of cells
      READ(30,*) Mus           !  6        # of u/s BCs
      READ(30,*) Mds           !  7        # of d/s BCs
      READ(30,*) dTprint       !  8        Print time interval/ output time step (hr)
      READ(30,*) startrun      !  9        Julian day start of run
      READ(30,*) endrun        ! 10        Julian day end of run
      READ(30,*) Ntrib         ! 11        Number of tributaries
      READ(30,*) Ndiv          ! 12        Number of diversions.
      READ(30,*) Specg         ! 13        Specific gravity of sediment
      READ(30,*) Cp            ! 14        Calibration for wind re-suspension
      READ(30,*) dtwind        ! 15        Time step for  wind time series (hr)
      READ(30,*) windgages     ! 16        Number of windgages read in
      READ(30,*) raingages     ! 17        Number of precip gages read in
      READ(30,*) etgages       ! 18        Number of ET gages read in
      READ(30,*) tidegages     ! 19        Number of tide gages read in
!      READ(30,*) Ndttrib       !           Tributary Q & TSS u/s time interval (day)
      READ(30,*) dttide        ! 20        Time Step for tidal time series (hr)
      READ(30,*) flocA         ! 21        clay flocculation parameter (McAnally A)
      READ(30,*) flocB         ! 22        clay flocculation parameter (McAnally B)
      READ(30,*) flocN         ! 23        clay flocculation exponent (McAnally n)
      READ(30,*) flocM         ! 24        clay flocculation exponent (McAnally m)
      READ(30,*) flocC1        ! 25        clay concentration where floculation settling begins (McAnally C1)(kg/m**3)
      READ(30,*) flocC3        ! 26        clay concentration where floculation settling becominbegins to hinder settling (McAnally C3) (kg/m**3)
      READ(30,*) Csalmax       ! 27        upper limit to salinity impact on floc (salinity concentration - ppt)
      READ(30,*) Pflocmax      ! 28        upper limit on portion of fines available to floc (ratio - unitless)
      READ(30,*) D50(1)        ! 29        sand D50 (m)                  -- USDA particle size definition: 0.05 mm < sand < 2 mm
      READ(30,*) D50(2)        ! 30        silt D50 (m)                  -- USDA particle size definition: 0.002 mm < silt < 0.05 mm
      READ(30,*) D50(3)        ! 31        unflocculated clay D50 (m)    -- USDA particle size definition: clay < 0.002 mm
      READ(30,*) D50(4)        ! 32        flocculated clay D50(m)       -- USDA particle size definition: clay < 0.002 mm
      READ(30,*) D90x          ! 33        representative D90/D50 ratio
      READ(30,*) Tcrit(1)      ! 34        critical shear stress for sand particles (N/m**2)
      READ(30,*) Tcrit(2)      ! 35        critical shear stress for silt particlss (N/m**2)
      READ(30,*) Tcrit(3)      ! 36        critical shear stress for unflocculated clay particles (N/m**2)
      READ(30,*) Tcrit(4)      ! 37        critical shear stress for flocculated clay particles (N/m**2)
      READ(30,*) saltox        ! 38        Salinity concentration at which point algal growth is halved
      READ(30,*) kphy20        ! 39        temperature-dependent phytoplankton mortality rate @ 20 degC (kphy20) (1/day)
      READ(30,*) thetaphy      ! 40        temperature-dependent phytoplankton mortality rate constant (thetaphy)
      READ(30,*) kdet20        ! 41        temperature-dependent detritus hydrolysis rate @ 20 degC (kdet20) (1/day)
      READ(30,*) thetadet      ! 42        temperature-dependent detritus hydrolysis rate constant (thetaDET)
      READ(30,*) kresp20       ! 43        temperature-dependent phytoplankton respiration rate @ 20 degC (kresp20) (1/day)
      READ(30,*) thetaresp     ! 44        temperature-dependent phytoplankton respirtion rate constant (thetaresp)
      READ(30,*) kpo420        ! 45        temperature-dependent orthophosphate release rate @ 20 degC (kpo4) (1/day)
      READ(30,*) thetapo4      ! 46        temperature-dependent orthophosphate rate coefficient (thetapo4)
      READ(30,*) kdenit20      ! 47        temperature-dependent denitrification rate @ 20 degC (kdenit20) (1/day)
      READ(30,*) thetadenit    ! 48        temperature-dependent denitrification rate coefficient (thetadenit) (1/day)
      READ(30,*) kphot20       ! 49        temperature-dependent photosynthetic rate @ 20 degC (kphot20)(1/day)
      READ(30,*) thetaphot     ! 50        temperature-dependent photosynthetic rate coefficient (thetaphot)(1/day)
      READ(30,*) kdon20        ! 51        temperature-dependent dissolved organic N hydrolysis rate @ 20 degC (kdon20)(1/day)
      READ(30,*) thetadon      ! 52        temperature-dependent dissolved organic N hydrolysis rate coefficient (thetadon)(1/day)
      READ(30,*) knit20num     ! 53        depth-dependent nitrification rate numerator @ 20 degC (knit20num) (knit20 = knit20num/depth) (1/day)
      READ(30,*) thetanit      ! 54        temperature-dependent nitrification rate coefficient (thetanit) (1/day)
      READ(30,*) knitmin       ! 55        minimum threshold to apply to calculated nitrification rate @ 20-deg (knitmin) (1/day)
      READ(30,*) knitmax       ! 56        maximum threshold to apply to calculated nitrification rate @ 20-deg(knitmax) (1/day)
      READ(30,*) KKexp         ! 57        exponent on Kadlec & Knight depth term for flow into marsh (KKexp)
      READ(30,*) defbedelev    ! 58        default bed elevation of open water compartments - used to fill any missing bed elevation values (m NAVD88)
      READ(30,*) defmarelev    ! 59        default marsh elevation of open water compartments - used to fill any missing marsh elevation values (m NAVD88)
      READ(30,*) def_n         ! 60        default roughness value to use if input roughness for link is less than 0.0, or greater than 1.0
      READ(30,*) maxconnect    ! 61        maximum number of hydraulic connections per compartment (must be <= 20)
      READ(30,*) stagemax      ! 62        maximum water surface elevation value allowed to be passed to other ICM routines (m NAVD88) (999999999 if not to be used)
      READ(30,*) salmax        ! 63        maximum salinity concentration value allowed to be passed to other ICM routines (ppt) (999999999 if not to be used)
      READ(30,*) tmpmax        ! 64        maximum water temperature value allowed to be passed to other ICM routines (deg C) (999999999 if not to be used)
      READ(30,*) sedmaxow      ! 65        maximum annual open water sediment accumulation allowed to be passed to other ICM routines (kg/m2/yr) (999999999 if not to be used)
      READ(30,*) sedmaxmi      ! 66        maximum annual marsh interior sediment accumulation allowed to be passed to other ICM routines (kg/m2/yr) (999999999 if not to be used)
      READ(30,*) sedmaxme      ! 67        maximum annual marsh edge sediment accumulation allowed to be passed to other ICM routines (kg/m2/yr) (999999999 if not to be used)
      READ(30,*) prismmax      ! 68        maximum tidal prism volume allowed to be passed to other ICM routines (m3) (999999999 if not to be used)
      READ(30,*) rangemax      ! 69        maximum range in water surface elevation allowed to be passed to other ICM routines (m) (999999999 if not to be used)
      READ(30,*) depthmax      ! 70        maximum water depth value allowed to be passed to other ICM routines (m) (-9999 if not to be used)
      READ(30,*) algmax        ! 71        maximum algae concentration value allowed to be passed to other ICM routines (kg/m3) (999999999 if not to be used)
      READ(30,*) tssmax        ! 72        maximum TSS concentration value allowed to be passed to other ICM routines (kg/m3) (999999999 if not to be used)
      READ(30,*) tknmax        ! 73        maximum TKN concentration value allowed to be passed to other ICM routines (kg/m3) (999999999 if not to be used)
      READ(30,*) nlinksw       ! 74        number of links to write to FLO.out file (must match number of links in 'links_to_write.csv', set to 0 if no link flows are to be written)
      READ(30,*) nstghr        ! 75        number of compartments to write to STGhr.out file (must match number of compartments in 'hourly_stage_to_write.csv', set to 0 if no hourly stages are to be written)
      READ(30,*) modeloverland ! 76        should overland flow be modeled in link types 8 and 9? (1=yes, 0=no)
      READ(30,*) minwater      ! 77        minimum water area required in each compartment (sq. meter)
      READ(30,*) oscilflag     ! 78        vertical oscillation in water level between timesteps that will flag 'dz' warning
      READ(30,*) maxdz         ! 79        maximum change in water level allowed between between timesteps - unrealistic but set to flag instabilities (maxdz) (meters)
      READ(30,*) iTemp         ! 80        modelling Temperature (No=0, Yes=1)
!      READ(30,*) temp0         !          initial system water temp
      READ(30,*) iSed          ! 81        modelling Sediment (No=0, Yes=1)
      READ(30,*) iSal          ! 82        modelling Salinity (No=0, Yes=1)
      READ(30,*) iWQ           ! 83        modelling Water Quality chemical constituents (No=0, Yes=1)
      READ(30,*) fpet          ! 84        overwater evapotranspiration model (fpet) (1.0 = use PET input data 0.0 = Use average of PET input data)
      READ(30,*) fa_def        ! 85        upwind factor for weighting upwind/downwind contributions to WQ,salinity,TSS and temp transport calculations (fa_def) (range from 0.5 to 1) (0.5 gives equal weight to upwind and downwind, 1 will only use upwind)
      READ(30,*) upwind_vel    ! 86        channel velocity at which upwind factor, fa, is set to 1.0 (upwind_vel) (m/sec)
      READ(30,*) fe            ! 87        dispersion calibration factor (fe)
      READ(30,*) lockdays      ! 88        number of days in lock control observation timeseries data file LockControlObservedData.csv (lockdays) (days)
      READ(30,*) dtlock        ! 89        timestep of lock control observation timeseries data (dtlock) (hours)
      READ(30,*) CSSmin        ! 90        minimum CSS concentration allowed to be reached - calculated CSS below this value will be increased - this is applied to each sediment class separately (cssmin) (mg/L)
      READ(30,*) CSSmax        ! 91        maximum CSS concentration allowed to be reached - calculated CSS above this value will be throttled - this is applied to each sediment class separately (cssmax) (mg/L)
      READ(30,*) TSSresusOff   ! 92        TSS concentration threshold - used to determine maximum concentration of suspended sediment for each size class which deactivates bed resuspension (mg/L)
      READ(30,*) sal_0_1_error ! 93        error term for freshest locations (0-1 ppt) - to be used to perturb output files (sal_0_1_error)
      READ(30,*) sal_1_5_error ! 94        error term for intermediate salinity locations (1-5 ppt) - to be used to perturb output files (sal_1_5_error)
      READ(30,*) sal_5_20_error  ! 95       error term for saline locations (5-20 ppt) - to be used to perturb output files (sal_5_20_error)
      READ(30,*) sal_20_35_error ! 96       error term for most saline locations (>20 ppt) - to be used to perturb output files (sal_20_35_error)
      READ(30,*) tss_error       ! 97       error term for TSS - to be used to perturb output files -percentage adjustment (between -1 and 1) (tss_error)
      READ(30,*) stage_error     ! 98       error term for stage - to be used to perturb output files (stage_error)
      READ(30,*) stvar_error     ! 99       error term for water level variability - to be used to perturb output files (stvar_error)
      READ(30,*) nlinklimiter    ! 100      number of links to apply flow limiter (must match number of links in 'links_to_apply.csv', set to 0 if no link to apply flow limiter)
      READ(30,*) ntc             ! 101      number of terminal connection
      READ(30,*) nlc             ! 102      number of lateral connection
      close(30)
!      sal_0_1_error  = 0.0
!      sal_1_5_error  = 0.0
!      sal_5_20_error  = 0.0
!      sal_20_35_error  = 0.0
!      tss_error  = 0.0
!      stage_error  = 0.0
!      stvar_error  = 0.0

!>> Determine number of simulation timesteps
      simdays = endrun-startrun+1
      NTs =int(simdays*24*60*60/dt)
      NNN=NTs-int(24*60*60/dt)
      lastdaystep = 24*60*60/dt !dt is in seconds - number of timesteps in a day
!>> Determine number of simulation steps needed before updating tide and wind data
      lasttidestep = dttide*60*60/dt !dttide is in hours - number of timesteps before updating tide data
      lastwindstep = dtwind*60*60/dt !dtwind is in hours - number of timesteps before updating wind data
      lastlockstep = dtlock*60*60/dt !dtlock is in hours - number of timesteps before updating lock control data
!> Initialize counters for updating tide and wind data - initial values are set to 1 instead of zerob/c first day starts one timestep AFTER midnight - when midnight is hit at end of day this counter is then reset
      daystep = 0  !YW! Modified to fix the issue with skipping the first time step
      tidestep = 0 !YW!
      windstep = 0 !YW!
      lockstep = 0 !YW!

!>> Call 'allocate_params' subroutine which allocates memory for almost all arrays used by this program
!>> (some temporary arrays are allocated in ICM_InterpolateToGrid subroutine which postprocesses model results).

      call allocate_params
!>> Set initial conditions for constants and model parameters that are not included in input text files
!	Initial Conditions
	g=9.81					! Gravity (m/s2)
	pi=4.0*atan(1.0)
	TemI = 15.				! Initial Water Temperature
	KnN= 20.				! (ug/L) DIN Michaelis Constant  Thomann & Mueller
	KnP= 3.					! (ug/L) P Michaelis Constant
	KnSS= 50.				! (mg/L) SS Michaelis Constant  chged 30 to 50 JAM March 2011
	KnSal=4.				! (ppt) Salinity  Michaelis Constant

	ParP= 0.4               ! SRP/TP in Tribs   J. Day 1994 BCS trial
	ParDOP= 0.1             ! DOP/TP in Tribs   J. Day 1994 BCS trial
	PARPOP=1.-ParDOP-ParP   ! POP/TP in Tribs   J. Day 1994 BCS trial
	ParPMR= 0.2             ! SRP/TP in Diversions J. Day 1994 BCS trial
	Pardpmr= 0.05
	PPMR=1.-ParPMR-Pardpmr
	ParSand=0.05			! sand/TSS in Tribs and MR typical
	ParCLa=0.03				! Partition LivA --> ChlA

	consd=24*3600			! sec to days
	conv=0.001*consd		! (mg/L)*(m3/s) --> kg/d
	floodf(:)=0.0
	nuo=0.000001			! DEFAULT Viscosity

!>> Generate array for Day of Year value for the first day of each month
      month_DOY(1) = 1
      month_DOY(2) = 32
      month_DOY(3) = 60
      month_DOY(4) = 91
      month_DOY(5) = 121
      month_DOY(6) = 152
      month_DOY(7) = 182
      month_DOY(8) = 213
      month_DOY(9) = 244
      month_DOY(10) = 274
      month_DOY(11) = 305
      month_DOY(12) = 335
!>> Update first day of month for March through December during a leap year
      if(simdays == 366) then
          do j = 3,12
              month_DOY(j) = month_DOY(j)+1
          enddo
      endif

C> initially set GrowAlgae array equal to zero
	do j=1,N
	 do ichem=1,14
	  do me=1,14
	   GrowAlgae(j,ichem,me) = 0.
	  enddo
	 enddo
	enddo

!>> determine upper CSS threshold for each sediment class where resuspension of bed material will be deactivated
!>> These values are based on an assumption that at a TSS value equal to the threshold concentration, the particle size distribution is 10% sand, 45% silt, and 45% clay.
      CSSresusOff(1) = TSSresusOff*0.1
      CSSresusOff(2) = TSSresusOff*0.45
      CSSresusOff(3) = TSSresusOff*0.45
      CSSresusOff(4) = TSSresusOff*0.45


!>> Open input text files.

      open (unit=32, file= 'Cells.csv', status = 'unknown')
      open (unit=323, file ='Fetch.csv', status = 'unknown')
      open (unit=33, file= 'Links.csv', status = 'unknown')			! JAM Oct 2010
      open (unit=34, file= 'LinksClosedHours.dat', status = 'unknown') !-EDW !value of 1 means link is closed for the hour, zero means it is open
      open (unit=35,file='LockControlObservedData.csv',status='unknown')
!      open (unit=36, file= 'SWR.dat', status = 'unknown')
      open (unit=74, file= 'MissRToC.csv', status = 'unknown')		! JAM Oct 2010
      open (unit=39, file= 'TribQ.csv', status = 'unknown')   !Tributary flow (m3/s)
      open (unit=40, file= 'PET.csv', status = 'unknown')
      open (unit=42, file= 'Precip.csv', status='unknown')
      open (unit=45, file= 'Meteorology.csv', status = 'unknown')
      open (unit=44, file= 'AnthL.csv', status = 'unknown')			! Farm and Urban WW Loads (kg/d)
      open (unit=43, file= 'WindVectorsX.csv',status= 'unknown')
      open (unit=46, file= 'WindVectorsY.csv',status= 'unknown')
      open (unit=47, file= 'TideData.csv',status='unknown')
      open (unit=48, file= 'TideTranspose.csv',status='unknown')
      open (unit=49, file= 'TideWeight.csv',status='unknown')
      open (unit=77, file= 'QMult.csv', form= 'formatted',
     &       status = 'unknown')
      open (unit=55, file= 'TribS.csv', status = 'unknown')	! Tributary sand concentration (mg/L)
      open (unit=555,file= 'TribF.csv',status='unknown')      ! Tributary fines concentration (mg/L)
      open (unit=56, file= 'SBC.dat', status = 'unknown')
      open (unit=80, file= 'NO2NO3Data.csv', status = 'unknown')
      open (unit=81, file= 'NH4Data.csv', status = 'unknown')
      open (unit=82, file= 'OrgNData.csv', status = 'unknown')
      open (unit=83, file= 'PhosphorusData.csv', status ='unknown')
      open (unit=84, file= 'AtmChemData.csv', status ='unknown')
      open (unit=85, file= 'Decay.csv', status ='unknown')
      open (unit=86, file= 'DivQm.csv', status ='unknown')	! diversion flow multiplier on Miss Riv flow
      open (unit=87, file= 'DivWQ.csv', status ='unknown')
      open (unit=88, file= 'QMult_div.csv', form= 'formatted',
     &       status ='unknown')
      open (unit=89, file= 'DivSW.csv', status = 'unknown')
      open (unit=101, file= 'BCToC2.dat', form = 'formatted')
      open (unit=110, file= 'surge.csv', form = 'formatted')
!      open (unit=117, file= 'AsedOW.csv',form ='formatted',		! Sediment Accretion  !Status='unknown' added by Joao Pereira 5/17/2011
!     &       status ='unknown')
!      open (unit=118, file= 'UplandNP.dat', form ='formatted')
      open (unit=125, file= 'KBC.dat', status = 'unknown')		        !node numbers of open boundary
      open (unit=126, file= 'links_to_write.csv',status='unknown')	    !input csv file with the link ID numbers of links to write flowrate to output file
      open (unit=127, file='hourly_stage_to_write.csv',status='unknown')  !input csv file with the ID numbers of compartments to write hourly stage to output file

!>> Open output text files (in append mode, if needed).
      open (unit=70,file='DIN.out',form ='formatted',position='append')
      open (unit=71,file='OrgN.out',form='formatted',position='append')
      open (unit=72,file='TPH.out',form='formatted',position='append')			! TP.out
      open (unit=73,file='TOC.out',form='formatted',position='append')
      open (unit=75,file='SAL.out',form ='formatted',position='append')			! Salinity.out
      open (unit=91,file='NO3.out',form='formatted',position='append')			! NO2NO3.out
      open (unit=92,file='NH4.out',form = 'formatted',position='append')		
      open (unit=93,file='O2Sat.out',form='formatted',position='append')
      open (unit=94,file='ALG.out',form='formatted',position='append')
      open (unit=95,file='DO.out',form='formatted',position='append')
      open (unit=96,file='TSS.out',form='formatted',position='append')
      open (unit=97,file='DET.out',form='formatted',position='append')			! DeadAlgae.out
      open (unit=100,file='TMP.out',form='formatted',position='append')		   
	open (unit=103,file='SedAcc.out',form='formatted',position='append')		! Last row will be used to compute 20yr open water Acc. 
	open (unit=105,file='fflood.out',form='formatted',position='append')
      open (unit=111,file='STG.out',form='formatted',position='append')			! ESAVE.OUT
      open (unit=112,file='TRG.out',form='formatted',position='append')			! Range.out
      open (unit=113,file='DON.out',form='formatted',position='append')
      open (unit=119,file='SPH.out',form='formatted',position='append')			! SRP.out
      open (unit=121,file='NRM.out',form='formatted',position='append')			! NRAcc.out -> Denitrification
      open (unit=123,file='TKN.out',form='formatted',position='append')
      open (unit=124,file='FLOm.out',form='formatted',position='append')

! read in information for grid cells used to pass data to other ICM routines !-EDW
      open (unit=200, file='grid_lookup_500m.csv', form='formatted')              ! compartment and link lookup table for 500-m grid cells
      open (unit=201, file='grid_interp_dist_500m.csv',form='formatted')          ! distance from each 500-m grid cell centroid to the compartment and link centroids
      open (unit=202, file='grid_data_500m.csv', form='formatted')                ! mean elevation for 500 m grid cells     
      open (unit=203, file='grid_IDs_Veg_matrix.csv', form='formatted')           ! 500m grid cell names formatted in the matrix used by Vegetation ICM routine

      open (unit=204, file='grid_500m_out.csv', form='formatted')              ! output file for 500 m grid cells - in list form
      open (unit=205, file='compartment_out.csv',form='formatted')             ! output file for hydro compartments - summary values for ICM in list form

      open(unit=206,file='sal_monthly_ave_500m.out',form='formatted')
      open(unit=207,file='tmp_monthly_ave_500m.out',form='formatted')
      open(unit=208,file='tkn_monthly_ave_500m.out',form='formatted')
      open(unit=209,file='TSS_monthly_ave_500m.out',form='formatted')

      open(unit=210,file='STGhr.out',form='formatted',position='append')				!output file for hourly water level in Boundary Condition cells
      open(unit=211,file='FLO.out',form='formatted',position='append')		!output file for flowrate	
      open(unit=212,file='STGm.out',form='formatted',position='append')
! output files for use in the Vegetation ICM routine !-EDW
! these are written in append mode. ICM checks when first run as to whethere these files exist.
! If Ecohydro is run outside of the ICM these files may be erroneously appended to if they contain data and the model is re-run.
      open (unit=301, file=TRIM(ADJUSTL(VegMeanSalFile)),         ! mean salinity formatted for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=302, file=TRIM(ADJUSTL(VegSummerSalFile)),       ! mean summertime salinity for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=303, file=TRIM(ADJUSTL(VegSummerDepthFile)),     ! mean summertime water depth for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=304, file=TRIM(ADJUSTL(VegWaveAmpFile)),         ! variance in daily water depth for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=305, file=TRIM(ADJUSTL(VegSummerTempFile)),      ! mean summertime water temperature for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=306, file=TRIM(ADJUSTL(VegTreeEstCondFile)),     ! tree establishment conditions for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=307, file=TRIM(ADJUSTL(VegBIHeightFile)),      ! tree establishment conditions for input into Vegetation ICM routine
     &        form='formatted', position='append')
      open (unit=308, file=TRIM(ADJUSTL(VegPerLandFile)),      ! tree establishment conditions for input into Vegetation ICM routine
     &        form='formatted', position='append')
! hotstart files
      open(unit=400,file='hotstart_in.dat',form='formatted')
      open(unit=401,file='hotstart_out.dat',form='formatted')

!>> Read Boundary Conditions file
      Read(125,*)(KBC(jj), jj=1,mds) !AMc Oct 8 2013
      close(125)

!>> !YW! Read input link file to apply flow limiter
      open (unit=500, file= 'links_to_apply.csv',status='unknown')
      read(500,*)
      do kk = 1,nlinklimiter
          read(500,*) linkslimiter(kk)
      enddo

!>> 1D-ICM coupling input files
      if (ntc>0) then
        write(*,*) 'The number of terminal connection is ',ntc
        open (unit=402, file= '1D2Dcoupling_tc.dat', status = 'unknown')
        read(402,*)                                                           ! dump header row of compartment input file
	    do i = 1,ntc
		    read(402,*) tc2D(i)  ! read compartment list for the terminal connection
	    enddo
	    do i = 1,ntc
		    read(402,*) tc1D(i) ! read 1D CX list for the terminal connection
	    enddo
        close(402)
        
        do i = 1,ntc
            tcH(i) = 0.0
            tcQ(i) = 0.0
        enddo
      endif
      
      if (nlc>0) then
        write(*,*) 'The number of lateral connection is ',nlc
        open (unit=403, file= '1D2Dcoupling_lc.dat', status = 'unknown')
        read(403,*)                                                          ! dump header row of compartment input file
        read(403,*) 
	    do i = 1,nlc
		    read(403,*) lcr2D(i)   ! read the list of lateral flow receiving ICM compartment
        enddo
        read(403,*) 
	    do i = 1,nlc
		    read(403,*) lcf2D(i)   ! read the list of lateral flow dummy ICM compartment
        enddo
        read(403,*) 
	    do i = 1,nlc
		    read(403,*) lcl2D(i)   ! read the list of ICM link
        enddo 
        read(403,*) 
	    do i = 1,nlc
		    read(403,*) lc1D(i)   ! read the list of 1D channel cx
	    enddo
        close(403)

        do i = 1,nlc
            lcH(i) = 0.0
            lcQ(i) = 0.0
        enddo        
      endif


! Initialize some variables and arrays
      NR(:)=0.0
      stds=0.
      do j=1,N
          accsed(j)=0.0
          sal_ave(j)=0.0
          cumul_retreat(j) = 0.0
      enddo

      do i=1,M
          asedout(i)=0.0
      enddo

!>> Call 'infile' subroutine, which reads input text files and stores data into arrays
      call infile


!>> write header row for hourly stage output file (since not all compartments are printed)	
      if (nstghr > 0) then
          write(210,908) 'Compartment:',(stghrwrite(jjk),jjk=1,nstghr)
      else
          write(210,*) 'Hourly water level not saved to file.'
      endif
908	FORMAT(A,<nstghr-1>(I0,','),I0) ! first column has 'Compartment:##', followed by comma delimited list of boundary compartments

!>> write header row for flowrate output file (since not all links are printed)
	if (nlinksw > 0) then
          write(211,909) 'Link:',(linkswrite(jjk),jjk=1,nlinksw)
      else
          write(211,*) 'No links chosen to have flowrate outputs saved.'
      endif
909	format(A,<nlinksw-1>(I0,','),I0) ! first column has 'Link:##', followed by comma delimited list of links

!>> Close input files that were imported in infile subroutine
      close(32)
      close(33)
      close(34)
      close(39)
      close(40)
      close(42)
      close(44)
      close(45)
      close(43)
      close(46)
      close(47)
      close(48)
      close(49)
      close(55)
      close(56)
      close(74)
      close(77)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)
      close(90)
      close(101)
      close(110)
!      close(118)
      close(200)
      close(201)
      close(202)

!>> Take first timestep of imported wind data and save into windx and windy arrays.
!>> These arrays will be overwritten at a delta t that matches the wind data timestep (this update occurs immediately prior to calling hydrod)
!>> 'windrow' is a counter that is incrementally updated each time a wind data timestep is reached
      windrow = 1
      do jjj = 1,N
          windx(jjj) = windx_data(windrow,jwind(jjj))
          windy(jjj) = windy_data(windrow,jwind(jjj))
      enddo

!>> 'tiderow','surgerow', and 'lockrow' are counters that are incrementally updated each time a tide or lock control data timestep is reached
      tiderow = 0	!YW! Modified to match all other initialization
      surgerow = 0	!YW!
      lockrow = 0	!YW!

!>> Set initial conditions for links (from input files)
      do i=1,M
          fa(i) = fa_def*fa_mult(i) !Set array of initial upwind factor to default value
          Q(i,1)=0.0	!YW!
          Q(i,2)=Q(i,1)

          ! MP2023 zw added 04/06/2020
          EAOL(i)=0.0
          !FLO(i)=0.0 !YW! flo range nlinksw
          SlkAve(i) = 0
          SL(i,2)=0
          TL(i,1)=0
      enddo


!>> Initialize hardcoded salinity and stage control trigger flags
      SalLockStatusHNC = 1  !Open = 1 Close = -1
      SalLockTriggerHNC = 1
      StgTriggerSuperiorCanal = 1
      StgTriggerStatusSuperiorCanal = 1
      !SalLockStatusCSC = 1
      !SalLockTriggerCSC = 1
      !Atch_div_onoff = 1


!>> Read in hotstart file and set initial conditions (will overwrite some ICs set previously from input files)
      write(1,*)
      write(1,*)'-----------------------------------------------'
      write(1,*)'Reading in hotstart file and setting
     & values as initial conditons.'
      write(1,*)'-----------------------------------------------'
      write(*,*)
      write(*,*)'-----------------------------------------------'
      write(*,*)'Reading in hotstart file and setting
     & values as initial conditons.'
      write(*,*)'-----------------------------------------------'
      read(400,*)             ! ignore header row
      do j=1,N
          read(400,*) dump_int,   !no need to save compartment number - read in to a dummy integer variable
     &                    Es(j,1),
     &                    S(j,1),
     &                    Css(j,1,1),
     &                    Css(j,1,2),
     &                    Css(j,1,3),
     &                    Css(j,1,4),
     &                    Tempw(j,1),
     &                    Chem(j,1,1),
     &                    Chem(j,2,1),
     &                    Chem(j,3,1),
     &                    Chem(j,4,1),
     &                    Chem(j,5,1),
     &                    Chem(j,6,1),
     &                    Chem(j,7,1),
     &                    Chem(j,8,1),
     &                    Chem(j,9,1),
     &                    Chem(j,10,1),
     &                    Chem(j,11,1),
     &                    Chem(j,12,1),
     &                    Chem(j,13,1),
     &                    Chem(j,14,1) !Chem unit = mg/L

!>> Initialize some variables and arrays
          Sandacc(j,1) = 0.0
          Siltacc(j,1) = 0.0
          Clayacc(j,1) = 0.0

          As(j,2)= As(j,1)				! Surface area of cells (m2)
          Es(j,2)= Es(j,1)				! Stage in storage cells (m)
          ds(j,1)= Es(j,2)-Bed(j)			! Depth in storage cells (m)
          Eh(j,2)=Eh(j,1)					! Stage in Marsh storage (m)	!JAM Oct 2010
          BCnosurge(j,1) = 0.0            ! Initialize no surge BC to zero for all compartments - only BC nodes will be updated - rest of array will be 0.0
          BCnosurge(j,2) = BCnosurge(j,1) ! boundary conditions stage(m) before surge is added
          BCsurge(j,1) = 0.0              !YW! Initialize surge BC to zero for all compartments - only BC nodes will be updated - rest of array will be 0.0
          BCsurge(j,2) = BCsurge(j,1)
          Qmarsh(j,2) = Qmarsh(j,1)		! Flow in marsh cell	!JAM Oct 2010
          Qmarshmax(j) = 0.0
          S(j,2) = S(j,1)
          Tempw(j,2) = Tempw(j,1)

          do ichem=1,14    !zw 4/28/2015 added to replace above statements
              Chem(j,ichem,2)=Chem(j,ichem,1)
          enddo

          ! MP2023 zw added 04/06/2020
          SL(j,1) = 0
          SL(j,2) = 0
          do sedclass=1,4
               CSS(j,2,sedclass) = CSS(j,1,sedclass)
               CSSh(j,1,sedclass) = CSS(j,1,sedclass)
               CSSh(j,2,sedclass) = CSS(j,1,sedclass)
          enddo
          Sacc(j,1)=0.0
          Sacch_int(j,1)=0.0
          Sacch_edge(j,1)=0.0
          CSSvRs(j,1)= 0.0
          Sacc(j,2)=Sacc(j,1)
          Sacch_int(j,2)=Sacch_int(j,1)
          Sacch_edge(j,2)=Sacch_edge(j,1)
          Sandacc(j,2) = Sandacc(j,1)
          Siltacc(j,2) = Siltacc(j,1)
          Clayacc(j,2) = Clayacc(j,1)
          CSSvRs(j,2)= CSSvRs(j,1)

!>> Initialize average variable
          ESAV(j,1) = ES(j,2)*dt/(3600.*24.)
          EHAV(j,1) = EH(j,2)*dt/(3600.*24.)
          TSSave(j) = ( CSS(j,1,1) + CSS(j,1,2)
     &    + CSS(j,1,3) + CSS(j,1,4) )*dt/(3600.*24.)
          SALAV(j) = S(j,2)*dt/(3600.*24.)
          QmarshAve(j) = Qmarsh(j,1)*dt/(3600.*24.)
          TempwAve(j) = Tempw(j,1)*dt/(3600.*24.)

          ! MP2023 zw added 04/06/2020
          ESMX(j,2)=ES(j,2)
      		ESMN(j,2)=ES(j,2)
          dailyHW(j)=0.0
          dailyLW(j)=0.0
          do ichem = 1,14
               ChemAve(j,ichem) = Chem(j,ichem,2)*dt/(3600.*24.)
          enddo
          denit(j,2)=0

      enddo
      close(400)

      Emax=0.0
      Emin=0.0
      
!>> !YW! initialize BCnosurge and BCsurge with initial tide and surge      
      do jj=1,tidegages
          BCnosurge(transposed_tide(jj,1),1) = TideData(1,jj)
          do jjj=1,Mds
              if (KBC(jjj)==transposed_tide(jj,1)) then
                  BCsurge(transposed_tide(jj,1),1)= Surge(1,jjj)
              endif
          enddo
      enddo
              
      if(nlinksw > 0) then               !YW!
          do jj = 1,nlinksw
              FLO(jj) = 0.0
          enddo
      endif

 
!*****************************Start of the unsteady model run***********************************
!>> Start time stepping through model. MAIN LOOP THAT IS COMPLETED FOR EACH SIMULATION TIMESTEP.
!>> Main DO loop. Looped over each simulation timestep.
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'START MAIN HYDRODYNAMIC MODEL'
      write(1,*) '----------------------------------------------------'

      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'START MAIN HYDRODYNAMIC MODEL'
      write(*,*) '----------------------------------------------------'

      NTs_Ratio = float(Nts)/float(ntim)
      if ((NTs_Ratio < 1) .OR. (mod(NTs_Ratio,1.0) .NE. 0)) then
          write(*,*) 'Check time step for 1D and 2D model'
          write(*,*) '1D 2D time step ratio is ',NTs_Ratio
          write(*,*) '1D model number of time step is ',ntim
          write(*,*) '2D model number of time step is ',NTs
          pause
      endif
      write(*,*) '1D 2D time step ratio is ',NTs_Ratio

!	  time step of MESH n_1d=1
      n_1d=1  !!! NTs is the no of time steps in ICM while ntim is the no of time steps in MESH.

      do  mm= 1, NTs  						! ******do loop ends at ~ line 438    !YW! change NTs-1 to NTs
          t=float(mm)*dt						! lapse time in seconds  JAM 5/25/2011

!>> calculate various versions of time to be used as flags throughout program
          day = t/3600./24.
          dday=day-int(day)       ! decimal portion of day, dday=0.0 at 0:00 (midnight), 0.99 at end of day
          ddayhr = dday*24.       ! decimal portion of day converted to hours dday = 0.00 @ midnight, 23.99 at end of day

!>> update counters for current timestep in day, tide data, and wind data
!>> counters are reset to zero at end of main.f timestepping Do loop when value meets the laststep values calculated at start of main.f
          daystep = daystep + 1
          tidestep = tidestep + 1
          windstep = windstep + 1
          lockstep = lockstep + 1

!>> -- Update counter that tracks with row of the tide data to use for the current model timestep
          if (tidestep == 1) then
              tiderow = tiderow + 1
              surgerow = tiderow
          endif

!>> -- Update wind arrays if timestep matches wind data timestep
! initial wind data for timestep = 0 was read into windx and windy arrays prior to timestepping began
! all subsequent timesteps that are divisible by dtwind will update the wind data arrays
          if (windstep == 1) then
              windrow = windrow + 1
              do jjj = 1,N
                  windx(jjj) = windx_data(windrow,jwind(jjj))
                  windy(jjj) = windy_data(windrow,jwind(jjj))
              enddo
          endif

!>> -- Update lock control arrays if timesteps matches lock control data timestep
          if (lockstep == 1) then
              lockrow = lockrow + 1
          endif

!>> -- Loop over links and number of diversions and calculate sediment accumulation timeseries for each model link.
          kt=1+ifix(t/24./3600./365.25)
!		do j=1,N
!			do it=1,Ndiv
!				ASANDD(j,kt)=ParSandD(kt)*Qdiv(it,kt)*cssTdiv(it,kt)
!    &				   *Qmultdiv(j,it)
!              enddo
!		enddo


!>> -- Call 'hydrod' subroutine, which is another 'control' routine which calls all hydrodynamic subroutines at each simulation timestep.
          call hydrod(mm)

!>> Saving lateral flow           
          if (nlc>0) then
              do i = 1,nlc                          ! loop over lateral flow compartment
                  lcQ(i) = -1.0*Q(lcl2D(i),2)
              enddo
          endif          

!>> *********************************Adding 1d start

          if (mod(float(mm)/NTs_Ratio,1.0) .eq. 0) then
        ! interpolation of boundaries at the desired time step
          newQ(1)     =r_interpol(USBoundary(1, :),USBoundary(2, :),		&
            ppp,t_1d+dtini)
          !print*, newQ(1)

!>> set up terminal connection for 1D
          if (ntc>0) then
              do i = 1,ntc
			      newY(ncomp) = Es(tc2D(i),2)           ! yw need to revise for multiple tc
			  enddo
		  else
		      newY(ncomp) =r_interpol(DSBoundary(1, :),DSBoundary(2, :),		&
            qqq,t_1d+dtini)
		  endif          

          !newY(ncomp) =r_interpol(DSBoundary(1, :),DSBoundary(2, :),		&
          !  qqq,t_1d+dtini)
          xt=newY(ncomp)
          newArea(ncomp)=r_interpol(ncompElevTable,ncompAreaTable,nel,		&
            xt)

        ! applying lateral flow at the oldQ
		  lateralFlow = 0
          do i=1,noLatFlow
              if (latFlowType(i) .eq. 1) then
                  latFlowValue = r_interpol(lateralFlowTable(1, :, i), &
                    lateralFlowTable(2, :, i),dataInEachLatFlow(i),t)
              elseif (latFlowType(i) .eq. 2) then
                  latFlowValue = r_interpol(lateralFlowTable(1, :, i), &
                    lateralFlowTable(2, :, i),dataInEachLatFlow(i),oldQ(latFlowLocations(i)))

              elseif (latFlowType(i) .eq. 3) then
                  do j=1,nlc
                      if (latFlowLocations(i) == lc1D(j)) then
                          latFlowValue = lcQ(i)                                 ! Q coupling from ICM to MESH; (+)ve Q adds discharge to MESH
                          !lateralFlow(latFlowLocations(i)) = lcQ(i)            ! Q coupling from ICM to MESH; (+)ve Q adds discharge to MESH
                      endif
                  enddo
              endif


              latFlowValue = latFlowValue / &
                      sum(dx(latFlowLocations(i)-1:latFlowLocations(i)-1+latFlowXsecs(i)-1))

              do j=1,latFlowXsecs(i)
                lateralFlow(latFlowLocations(i)+j-1)=lateralFlow(latFlowLocations(i)+j-1) + latFlowValue
              enddo
            !print*, 'lat flow=', latFlowValue, latFlowLocations(i), &
            !    latFlowValue*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i)))*0.5, &
            !    oldQ(latFlowLocations(i))
            ! print*, 'lat flow at', latFlowLocations(i), &
            !    'flow=',lateralFlow(latFlowLocations(i))*0.5*(dx(latFlowLocations(i)-1)+dx(latFlowLocations(i))), &
            !    dx(latFlowLocations(i)-1), dx(latFlowLocations(i))
          end do
          !print*, 'noLatFlow',noLatFlow
          !print*, 'lateralFlow', (lateralFlow(i), i=1, ncomp)

          ! Set upstream discharge
          dqp(1) = newQ(1) - oldQ(1)
		  dap(1) = 0.0
          call section()
        ! Nazmul: The subroutine calls the attribute tables and interpolate according to the available water level
          thes=thetas

          call matrixp()

          do i=2,ncomp
            !cour=dt(i)/dx(i-1)
              cour=dtini/dx(i-1)
              rhs1=-cour*(f1(i)-f1(i-1)-d1(i)+d1(i-1))+lateralFlow(i)*dtini
              rhs2=-cour*(f2(i)-f2(i-1)-d2(i)+d2(i-1))+dtini*grav*(ci2(i)+aso(i))
              c11=g11inv(i)*b11(i-1)+g12inv(i)*b21(i-1)
              c12=g11inv(i)*b12(i-1)+g12inv(i)*b22(i-1)
              c21=g21inv(i)*b11(i-1)+g22inv(i)*b21(i-1)
              c22=g21inv(i)*b12(i-1)+g22inv(i)*b22(i-1)
              dap(i)=g11inv(i)*rhs1+g12inv(i)*rhs2-c11*dap(i-1)-c12*dqp(i-1)
              dqp(i)=g21inv(i)*rhs1+g22inv(i)*rhs2-c21*dap(i-1)-c22*dqp(i-1)
             ! print*,'rhs1=', rhs1, rhs2, c11, c12, c21, c22
          end do

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

          call secpred()
          thes=thesinv
          call matrixc()

        do i=ncomp-1,1,-1
            !cour=dt(i)/dx(i)
            cour=dtini/dx(i)
            rhs1=-cour*(f1(i+1)-f1(i)-d1(i+1)+d1(i))+lateralFlow(i+1)*dtini
            rhs2=-cour*(f2(i+1)-f2(i)-d2(i+1)+d2(i))+dtini*grav*(ci2(i)+aso(i))
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
              do pp_1d=1,nel
                  elevTable(pp_1d) = xsec_tab(1,pp_1d,i)
                  areaTable(pp_1d) = xsec_tab(2,pp_1d,i)
              enddo

    !       interpolate the cross section attributes based on FINAL CALCULATED area
              xt=newArea(i)
              newY(i)=r_interpol(areaTable,elevTable,nel,xt)
!-------------------------------------

              newQ(i)=oldQ(i)+dq
              froud(i)=abs(newQ(i))/sqrt(grav*newArea(i)**3.0/bo(i))

          end do

          do i=1,ncomp-1
              courant(i)=(newQ(i)+newQ(i+1))/(newArea(i)+newArea(i+1))*		&
             dtini/dx(i)
          enddo

          if (maxCourant .lt. maxval (courant)) then
              maxCourant = maxval (courant)
          endif

          t_1d = t_1d + dtini
          !print "('- cycle',i9,'  completed')", n_1d

          avgY = avgY + newY
          avgQ = avgQ + newQ

          if (mod(n_1d,saveFrequency) .eq. 0) then
            avgY = avgY / saveFrequency
            avgQ = avgQ / saveFrequency
            write(608, 10) t_1d, (avgY(i), i=1,ncomp)
            write(609, 10) t_1d, (avgQ(i), i=1,ncomp)
            write(651, 10) t_1d, (newArea(i), i=1, ncomp)
            avgY = 0.0
            avgQ = 0.0
          end if

        !write(81, *) t, newArea(ncomp)

        ! update of Y, Q and Area vectors
          oldY   = newY
          oldQ   = newQ
          oldArea= newArea

	  ! increasing the time step of MESH by 1
	      n_1d=n_1d+1
          endif
!>> *********************************Adding 1d end

!>> update teminal and lateral flow connections          
          if (ntc>0) then                               ! terminal connection discharge
              do i = 1,ntc
                  tcQ(i) = newQ(ncomp)					! yw need to revise for multiple tc
                  !write(*,*) 'main'
                  !write(*,*) 'i',i 
                  !write(*,*) 'tcQ(i)',tcQ(i)
                  !write(*,*) 'ncomp',ncomp
                  !write(*,*) 'newQ(ncomp)',newQ(ncomp)
                  !pause
			  enddo
		  endif
          
          if (nlc>0) then
              do i=1,nlc
                  lcH(i) = newY(lc1D(i))                  ! WL lateral coupling from MESH to ICM. need to make sure 1D 2D have the same lateral flow location. Need to modify to calculate averaged stage from 1D
              end do
              !pause
          endif                

          !do i=1,noLatFlow
          !    if (latFlowType(i) .eq. 3) then
          !        WL_lat_from_MESH = newY(latFlowLocations(i))                  ! WL lateral coupling from MESH to ICM
          !    endif
          !end do

!          Q(895,2) = newQ(1)  ! HARDCORDED FOR ASSIGNING 1D flow as control for lateral flow connection. Needs to be revised.
!          Q(896,2) = newQ(1)  ! HARDCORDED FOR ASSIGNING 1D flow as control for lateral flow connection. Needs to be revised.
          

!>> -- Loop over links. Save calculated flowrates, Q as the initial condition for the next simulation timestep.
		do i=1,M
			Q(i,1)=Q(i,2)
		enddo


!>> -- Loop over compartments. Save numerous variables that were just calculated by 'hydrod' as the initial condition for the next simulation timestep.
		do j=1,N
			Es(j,1)=Es(j,2)
			S(j,1)=S(j,2)				! resetting ICs
              SL(j,1) = SL(j,2)
              BCnosurge(j,1) = BCnosurge(j,2)
              BCsurge(j,1) = BCsurge(j,2)                !YW!
              do sedclass=1,4
                  CSS(j,1,sedclass) = CSS(j,2,sedclass)
                  CSSh(j,1,sedclass) = CSSh(j,2,sedclass)
              enddo
			Tempw(j,1)=Tempw(j,2)
!             Age(j,1)=Age(j,2)
			Sacc(j,1)=Sacc(j,2)
			Sacch_int(j,1)=Sacch_int(j,2)
		    Sacch_edge(j,1)=Sacch_edge(j,2)
              Sandacc(j,1) = Sandacc(j,2)
              Siltacc(j,1) = Siltacc(j,2)
              Clayacc(j,1) = Clayacc(j,2)

              do ichem = 1,14
			    Chem(j,ichem,1)=Chem(j,ichem,2)
		    enddo
              CSSvRs(j,1)= CSSvRs(j,2)
              Qmarsh(j,1) = Qmarsh(j,2)       ! added EDW/ZW -02/16/2015
              Eh(j,1) = Eh(j,2)               ! added EDW/ZW - 02/16/2015
!>> Check if any water surface elevation values are NaN - if so, pause model run
              if(isNAN(Es(j,2))) then
                  write(1,*)'Compartment',j,
     &                    'WSEL is NaN @ end of timestep=',mm
                  write(*,*)'Compartment',j,
     &                    'WSEL is NaN @ end of timestep=',mm
                  pause
              endif


!!>> -- Check that some  calculated water quality values do not exceed default threshold values. If they do, set equal to the threshold
!              if(Chem(j,10,2) > 0.00025) then
!                  Chem(j,10,2) = 0.00025
!              endif
!
!              if(Chem(j,11,2) > 0.00025) then
!                    Chem(j,11,2) = 0.00025
!              endif
!
!              if(Chem(j,12,2) > 0.00025) then
!                  Chem(j,12,2) = 0.00025
!              endif
!!>> -- Repeat water quality threshold check, so that the initial conditions array for the next model timestep matches the 'current' array.
!              if(Chem(j,10,1) > 0.00025) then
!                  Chem(j,10,1)=0.00025
!              endif
!
!              if(Chem(j,11,1) > 0.00025) then
!                  Chem(j,11,1)=0.00025
!              endif
!
!              if(Chem(j,12,1) > 0.00025) then
!                  Chem(j,12,1)=0.00025
!              endif
!zw 4/28/2015 remove low and high pass filters
          enddo

!>> Check if any link flowrate values are NaN - if so, pause model run
          do i=1,M
			Q(i,1)=Q(i,2)
              if(isNAN(Q(i,2))) then
                  write(1,*)'Link',i,'flow is NaN @ end of timestep=',mm
                  write(*,*)'Link',i,'flow is NaN @ end of timestep=',mm
                  write(*,*) '  Linkt=',linkt(i)
                  pause
              endif
		enddo


!>> Reset tidestep counter because end of observed tidal timestep is met
          if (tidestep == lasttidestep) then
              tidestep = 0
          endif

!>> Reset windstep counter because end of observed tidal timestep is met
          if (windstep == lastwindstep) then
              windstep = 0
          endif

!>> Reset lockstep counter because end of observed lock control timestep i smet
          if (lockstep == lastlockstep) then
              lockstep = 0
          endif


!>> Reset daystep counter because end of day is met
          if (daystep == lastdaystep) then
              write(1,3333)'Year',year,'Day',int(day),'complete.'
              write(*,3333)'Year',year,'Day',int(day),'complete.'
              daystep = 0
          endif

!>> End main model DO loop that is looped over each simulation timestep
      enddo

!>> *********************************Adding 1d start
      close(608)
      close(609)
      close(651)
    !close(81)

      print*, 'dx', (dx(i), i=1, ncomp-1)
      print*, 'Froude', (froud(i), i=1, ncomp)
      print*, 'Bed', (z(i), i=1, ncomp)
      print*, 'newArea', (newArea(i), i=1, ncomp)
      print*, 'I2_corr', (ci2(i), i=1, ncomp)
      print*, 'Courant no', (courant(i), i=1, ncomp-1)
      print*, 'Maximum Courant no', maxCourant

!10  format(f12.2 , <ncomp>f12.2)
10    format(f12.2 , 1000f12.2)

      call cpu_time( t2 )
      print '("Time = ",f10.3," seconds.")',t2 - t1
      pause 202

!>> *********************************Adding 1d end

!>> Add error terms to last timestep value before writing hotstart file
      write(1,*) 'Adding error adjustment to last timestep values.'
      write(*,*) 'Adding error adjustment to last timestep values.'

      do j=1,N
          Es(j,2) = Es(j,2) + stage_error

          if(S(j,2) < 1.0) then
              S(j,2) = max(0.0,S(j,2) + sal_0_1_error)
          elseif(S(j,2) < 5.0) then
              S(j,2) = max(0.0,S(j,2) + sal_1_5_error)
          elseif(S(j,2) < 20.0) then
              S(j,2) = max(0.0,S(j,2) + sal_5_20_error)
          else
              S(j,2) = max(0.0,S(j,2) + sal_20_35_error)
          endif

          Css(j,2,1) = max(0.0,Css(j,2,1)*(1.0+tss_error))
          Css(j,2,2) = max(0.0,Css(j,2,2)*(1.0+tss_error))
          Css(j,2,3) = max(0.0,Css(j,2,3)*(1.0+tss_error))
          Css(j,2,4) = max(0.0,Css(j,2,4)*(1.0+tss_error))

      enddo

!>> Save output values of last simulation timestep in a hotstart file
      write(1,*)'Saving values from last timestep to a hotstart file.'
      write(*,*)'Saving values from last timestep to a hotstart file.'

      write(401,11140)'Compartment',
     &                'Stage',
     &                'Salinity',
     &                'CSS_sand',
     &                'CSS_silt',
     &                'CSS_clay',
     &                'CSS_clayfloc',
     &                'WaterTemp',
     &                'NO3_NO2',
     &                'NH4',
     &                'DIN',
     &                'OrgN',
     &                'TIP',
     &                'TOC',
     &                'DO',
     &                'LiveAlg',
     &                'DeadAlg',
     &                'DON',
     &                'DOP',
     &                'DIP',
     &                'ChlA',
     &                'TKN'

      do j=1,N
	    write(401,11142) j,
     &                    Es(j,2),
     &                    S(j,2),
     &                    Css(j,2,1),
     &                    Css(j,2,2),
     &                    Css(j,2,3),
     &                    Css(j,2,4),
     &                    Tempw(j,2),
     &                    min(Chem(j,1,2),1000.0),
     &                    min(Chem(j,2,2),1000.0),
     &                    min(Chem(j,3,2),1000.0),
     &                    min(Chem(j,4,2),1000.0),
     &                    min(Chem(j,5,2),1000.0),
     &                    min(Chem(j,6,2),1000.0),
     &                    min(Chem(j,7,2),1000.0),
     &                    min(Chem(j,8,2),1000.0),
     &                    min(Chem(j,9,2),1000.0),
     &                    min(Chem(j,10,2),1000.0),
     &                    min(Chem(j,11,2),1000.0),
     &                    min(Chem(j,12,2),1000.0),
     &                    min(Chem(j,13,2),1000.0),
     &                    min(Chem(j,14,2),1000.0) ! chem unit = mg/L
      enddo
      close(401)

! Write sediment accumulation in open water output file - 1 value per year      
!      pm1 = (1.-Apctmarsh(j))
!      WRITE(117,9229)(((ASandA(j)/As(j,1)+clams
!     &			+max(0.,Sacc(j,2)))*pm1				!modified June 20, 2011 JAM for testing removed /1000  added ASANDA(
!     &			+Sacch(j,2)*Apctmarsh(j)+476.), j=1,N)		!increase fines captured
!	close (117)

!>> Calculate percentage of sand in open water bed sediments
      do j = 1,N
        if((Sandacc(j,2)+Siltacc(j,2)+Clayacc(j,2))>0)then  !zw added 04/06/2020
          pct_sand_bed(j) = max(0.0,min(100.0,Sandacc(j,2)/
     &             (Sandacc(j,2)+Siltacc(j,2)+Clayacc(j,2))))
        else
          pct_sand_bed(j) = 0
        endif
      enddo

!>> Call 'ICM_Summaries' subroutine which calculates the average output values used by the other ICM routines (Wetland Morph, Vegetation, and Ecosytems)
      call ICM_Summaries

!>> Call 'ICM_SalinityThreshold_WM' subroutine which calculates the maximum two-week salinity average for each hydro compartment (used in  Marsh Collapse calculations in Wetland Morph)
      call ICM_SalinityThreshold_WM

!>> Call 'ICM_MapToGrid' subroutine which maps (without interpolation) the hydro compartment results to the gridded structrure used by the other ICM routines.
      call ICM_MapToGrid(1,500)               ! Map annual mean salinity to grid
      call ICM_MapToGrid(2,500)               ! Map summer salinity to grid
      call ICM_MapToGrid(3,500)               ! Interpolate annual mean temperature to grid
      call ICM_MapToGrid(4,500)               ! Interpolate summer mean temperature to grid
      call ICM_MapToGrid(5,500)               ! Map annual mean water stage to grid
      call ICM_MapToGrid(6,500)               ! Map summer water stage to grid
      call ICM_MapToGrid(7,500)               ! Map summer water stage variance to grid
      call ICM_MapToGrid(8,500)               ! Map summer maximum 2-wk mean salinity to grid
      call ICM_MapToGrid(9,500)               ! Map percent sand in bed to grid

	do kk=301,312
		call ICM_MapMonthlytoGrid(kk,500)		! Map monthly TKN values to grid
      enddo

	do kk=401,412
		call ICM_MapMonthlytoGrid(kk,500)		! Map monthly TSS values to grid
	enddo

!>> Calculate water depth for 500 m-grid cells.
      do j=1,n_500m_cells
          depth_500m(j) = stage_500m(j)-bed_elev_500m(j)
          depth_summer_500m(j) = stage_summer_500m(j)-bed_elev_500m(j)
          height_500m(j) = land_elev_500m(j)-stage_500m(j)
      enddo

!>> Call 'ICM_InterpolateToGrid' subroutine which interpolates hydro compartment and link results to the gridded structure used by the other ICM routines.
      call ICM_InterpolateToGrid(1,500)       ! Interpolate annual mean salinity to grid
      call ICM_InterpolateToGrid(2,500)       ! Interpolate summer mean salinity to grid
      call ICM_InterpolateToGrid(3,500)       ! Interpolate annual mean temperature to grid
      call ICM_InterpolateToGrid(4,500)       ! Interpolate summer mean temperature to grid
      call ICM_InterpolateToGrid(5,500)       ! Interpolate summer maximum 2-wk mean salinity to grid

	do kk=101,112
		call ICM_InterpolateMonthlyToGrid(kk,500)       ! Interpolate monthly mean salinity values to grid
	enddo

	do kk=201,212
		call ICM_InterpolateMonthlyToGrid(kk,500)       ! Interpolate monthly mean temperatures to grid
	enddo

!>> Call 'ICM_TreeConditions_Veg' subroutine which determines if tree establishment conditions are met for each grid cell of the Veg model.
      call ICM_TreeConditions_Veg

!>> Call 'ICM_Formatting' subroutine which formats the gridded output into versions directly digestable by other ICM routines.
      call ICM_Formatting

!>> Calculate average tidal prism for each compartment.
      do j = 1,N
          TRsum(j) = 0.0
          do kl=1,simdays
              TRsum(j) = TRsum(j) + tidal_range_daily(kl,j)
          enddo
          TRave(j) = TRsum(j)/simdays
          tidalprism_ave(j) = TRave(j)*As(j,1)
      enddo

!>> Write compartment summary output to file in list form - one row for each compartment
!>> THESE HEADERS ARE USED BY OTHER ICM ROUTINES - DO NOT CHANGE WITHOUT UPDATING ICM.PY, HSI.PY, & WM.PY
      write(205,1117)'ICM_ID',
     &        'max_annual_stage',
     &        'ave_annual_stage',
     &        'ave_stage_summer',
     &        'var_stage_summer',
     &        'ave_annual_salinity',
     &        'ave_salinity_summer',
     &        'sal_2wk_max',
     &        'ave_tmp',
     &        'ave_tmp_summer',
     &        'openwater_sed_accum',
     &        'marsh_int_sed_accum',
     &        'marsh_edge_sed_accum',
     &        'tidal_prism_ave',
     &        'ave_sepmar_stage',
     &        'ave_octapr_stage',
     &        'marsh_edge_erosion_rate',
     &        'ave_annual_tss',
     &        'stdev_annual_tss',
     &        'totalland_m2'

      do kj=1,N
          write(205,1118) kj,
     &        max(-stagemax,min(stagemax,stage_max(kj))),
     &        max(-stagemax,min(stagemax,stage_ave(kj))),
     &        max(-stagemax,min(stagemax,stage_ave_summer(kj))),
     &        max(-rangemax,min(rangemax,stage_var_summer(kj))),
     &        min(salmax,sal_ave(kj)),
     &        min(salmax,sal_ave_summer(kj)),
     &        min(salmax,sal_2wk_ave_max(kj)),
     &        min(tmpmax,tmp_ave(kj)),
     &        min(tmpmax,tmp_ave_summer(kj)),
     &        max(-sedmaxow,min(sedmaxow,Sacc(kj,2))),
     &        max(-sedmaxmi,min(sedmaxmi,Sacch_int(kj,2))),
     &        max(-sedmaxme,min(sedmaxme,Sacch_edge(kj,2))),
     &        tidalprism_ave(kj),
     &        max(-stagemax,min(stagemax,sepmar_stage(kj))),
     &        max(-stagemax,min(stagemax,octapr_stage(kj))),
     &        MEE(kj),
     &        tss_ave(kj),
     &        tss_var_annual(kj)**0.5, !stdev = sqrt(variance)
     &        max(0.0,(Atotal(kj)-As(kj,1)))
      enddo

!>> Write gridded output to file in list form - one row for each grid cell.
!>> Write header rows in grid output files
!>> THESE HEADERS ARE USED BY OTHER ICM ROUTINES - DO NOT CHANGE WITHOUT UPDATING ICM.PY, HSI.PY, & WM.PY

      write(204,1115)'GRID_ID',
     &   'compartment_ave_salinity_ppt',
     &    'IDW_ave_salinity_ppt',
     &    'compartment_ave_summer_salinity_ppt',
     &    'IDW_ave_summer_salinity_ppt',
     &    'compartment_max_2wk_summer_salinity_ppt',
     &    'IDW_max_2wk_summer_salinity_ppt',
     &    'bed_pct_sand',
     &    'compartment_ave_temp',
     &    'IDW_ave_temp',
     &    'compartment_ave_summer_temp',
     &    'IDW_ave_summer_temp',
     &    'stage_ave',
     &    'stage_summer_ave',
     &    'var_stage_summer',
     &    'ave_depth_summer',
     &    'ave_depth'
	write(206,1119)'GRID_ID',
     &    'sal_ave_jan',
     &    'sal_ave_feb',
     &    'sal_ave_mar',
     &	'sal_ave_apr',
     &    'sal_ave_may',
     &    'sal_ave_jun',
     &    'sal_ave_jul',
     &	'sal_ave_aug',
     &    'sal_ave_sep',
     &    'sal_ave_oct',
     &    'sal_ave_nov',
     &	'sal_ave_dec'
      write(207,1119)'GRID_ID',
     &    'tmp_ave_jan',
     &    'tmp_ave_feb',
     &	'tmp_ave_mar',
     &    'tmp_ave_apr',
     &    'tmp_ave_may',
     &    'tmp_ave_jun',
     &	'tmp_ave_jul',
     &    'tmp_ave_aug',
     &    'tmp_ave_sep',
     &    'tmp_ave_oct',
     &	'tmp_ave_nov',
     &    'tmp_ave_dec'
      write(208,1119)'GRID_ID',
     &    'tkn_ave_jan',
     &    'tkn_ave_feb',
     &	'tkn_ave_mar',
     &    'tkn_ave_apr',
     &    'tkn_ave_may',
     &    'tkn_ave_jun',
     &	'tkn_ave_jul',
     &    'tkn_ave_aug',
     &    'tkn_ave_sep',
     &    'tkn_ave_oct',
     &	'tkn_ave_nov',
     &    'tkn_ave_dec'
	write(209,1119)'GRID_ID',
     &    'TSS_ave_jan',
     &    'TSS_ave_feb',
     &    'TSS_ave_mar',
     &	'TSS_ave_apr',
     &    'TSS_ave_may',
     &    'TSS_ave_jun',
     &    'TSS_ave_jul',
     &	'TSS_ave_aug',
     &    'TSS_ave_sep',
     &    'TSS_ave_oct',
     &    'TSS_ave_nov',
     &	'TSS_ave_dec'

	! write various summary results in 500m grid output file
      do k=1,n_500m_cells
          write(204,1116) grid_lookup_500m(k,1),
     &    min(salmax,salinity_500m(k)),
     &    min(salmax,salinity_IDW_500m(k)),
     &    min(salmax,salinity_summer_500m(k)),
     &    min(salmax,salinity_summer_IDW_500m(k)),
     &    min(salmax,sal_thresh_500m(k)),
     &    min(salmax,sal_thresh_IDW_500m(k)),
     &    pct_sand_bed_500m(k),
     &    min(tmpmax,tmp_500m(k)),
     &    min(tmpmax,tmp_IDW_500m(k)),
     &    min(tmpmax,tmp_summer_500m(k)),
     &    min(tmpmax,tmp_summer_IDW_500m(k)),
     &    max(-stagemax,min(stagemax,stage_500m(k))),
     &    max(-stagemax,min(stagemax,stage_summer_500m(k))),
     &    max(-rangemax,min(rangemax,stage_var_summer_500m(k))),
     &    max(-depthmax,min(depthmax,depth_summer_500m(k))),
     &    max(-depthmax,min(depthmax,depth_500m(k)))

 		write(206,1120) k,
     &	   min(salmax,sal_IDW_500m_month(1,k)),
     &	   min(salmax,sal_IDW_500m_month(2,k)),
     &       min(salmax,sal_IDW_500m_month(3,k)),
     &	   min(salmax,sal_IDW_500m_month(4,k)),
     &	   min(salmax,sal_IDW_500m_month(5,k)),
     &	   min(salmax,sal_IDW_500m_month(6,k)),
     &	   min(salmax,sal_IDW_500m_month(7,k)),
     &	   min(salmax,sal_IDW_500m_month(8,k)),
     &	   min(salmax,sal_IDW_500m_month(9,k)),
     &	   min(salmax,sal_IDW_500m_month(10,k)),
     &	   min(salmax,sal_IDW_500m_month(11,k)),
     &	   min(salmax,sal_IDW_500m_month(12,k))

 		write(207,1120) k,
     &	   min(tmpmax,tmp_IDW_500m_month(1,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(2,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(3,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(4,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(5,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(6,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(7,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(8,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(9,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(10,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(11,k)),
     &	   min(tmpmax,tmp_IDW_500m_month(12,k))

 		write(208,1120) k,
     &	   min(tknmax,tkn_500m_month(1,k)),
     &	   min(tknmax,tkn_500m_month(2,k)),
     &	   min(tknmax,tkn_500m_month(3,k)),
     &	   min(tknmax,tkn_500m_month(4,k)),
     &	   min(tknmax,tkn_500m_month(5,k)),
     &	   min(tknmax,tkn_500m_month(6,k)),
     &	   min(tknmax,tkn_500m_month(7,k)),
     &	   min(tknmax,tkn_500m_month(8,k)),
     &	   min(tknmax,tkn_500m_month(9,k)),
     &	   min(tknmax,tkn_500m_month(10,k)),
     &	   min(tknmax,tkn_500m_month(11,k)),
     &	   min(tknmax,tkn_500m_month(12,k))

 		write(209,1120) k,TSS_500m_month(1,k),
     &	   min(tssmax,TSS_500m_month(2,k)),
     &	   min(tssmax,TSS_500m_month(3,k)),
     &	   min(tssmax,TSS_500m_month(4,k)),
     &	   min(tssmax,TSS_500m_month(5,k)),
     &	   min(tssmax,TSS_500m_month(6,k)),
     &	   min(tssmax,TSS_500m_month(7,k)),
     &	   min(tssmax,TSS_500m_month(8,k)),
     &	   min(tssmax,TSS_500m_month(9,k)),
     &	   min(tssmax,TSS_500m_month(10,k)),
     &	   min(tssmax,TSS_500m_month(11,k)),
     &	   min(tssmax,TSS_500m_month(12,k))
	enddo

!>> determine end time for calculating runtimes
      call SYSTEM_CLOCK(runtime_end,count_rate2,count_max2)
      runtime_s = dble(runtime_start)
      runtime_e = dble(runtime_end)
      ! Check that system_clock parameters did not change during model run
      if (count_rate2 == count_rate1) then
          if (count_max2 == count_max1) then
              if (runtime_e <= runtime_s) then
                  runtime_e = runtime_e + dble(count_max2)
              endif
              runtime = (runtime_e-runtime_s)/(dble(count_rate2)*60.)
              write(1,1113)' Model run time = ',runtime,' minutes.'
              write(*,1113)' Model run time = ',runtime,' minutes.'
          else
              write(1,*) ' System clock was updated during run.',
     &                ' Run time not calculated.'
              write(*,*) ' System clock was updated during run.',
     &                ' Run time not calculated.'
          endif
      else
          write(1,*) ' System clock was updated during run.',
     &                ' Run time not calculated.'
          write(*,*) ' System clock was updated during run.',
     &                ' Run time not calculated.'

      endif


1111  FORMAT(A,',',A,',',A)
1112  FORMAT(I4,2(',',F))
1113  FORMAT(A,F10.2,A)
11140 format(A,21(' ',A))
11141 format(I8,21(F))
11142 format(I8,21(F20.4))
1115  format(A,16(',',A))
1116  format(I0,16(',',F0.4))
1117  format(A,19(',',A))
1118  format(I0,19(',',F0.4))
1119  format(A,12(',',A))
1120  format(I0,12(',',F0.4))
3333  FORMAT(4x,A,1x,I4,1x,A,1x,I4,1x,A)
9229	FORMAT(<cells-1>(F0.2,','),F0.2)
!9229	format(1000(F0.2,','),F0.2)

!>> End model simulation.

      return
      end