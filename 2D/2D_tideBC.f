
!> @file
!> @brief This subroutine assigns observed water level to the downstream boundary condition cells.
!> @details This subroutine assigns observed water level to the downstream boundary condition cells in two different methods.
!> If a downstream boundary cell is associated with a timeseries of data, the data is used to set the water level.
!> The observed data is shifted in time to account for travel times between the observation station and the compartment centroid.
!> The time shift is an input parameter located in the input file 'TideTranspose.csv'.
!> If a downstream cell is not associated with an observation station, it's water level is assigned by the nearest two compartments
!> that area associated with an observation station. These two associated compartments' water levels are distance weighted, using
!> the input file 'TideWeight.csv'.

!> @param[in]     Mds
!> @param[in]     dttide
!> @param[in]     tidegages
!> @param[in]     weighted_tide(,)
!> @param[in]     transposed_tide(,)
!> @param         EastWght
!> @param         WestWght
!> @param         EastComp
!> @param	        WestComp
!> @param         EastBC
!> @param         WestBC
!> @param         use_row
!> @param         row_transpose

      Subroutine TideBC
	  use params

      implicit none
	  integer :: jjk,jj,nnn
      integer :: use_row,row_transpose
      real :: EastWght,WestWght,EastComp,WestComp,EastBC,WestBC

!>> Update boundary condition water levels for compartments WITH observed water level timeseries
	  do jjk=1,Mds
		jj=KBC(jjk)
		do nnn=1,tidegages
!>> Check if compartment has an observed tide timeseries
			if (jj==transposed_tide(nnn,1)) then
!>> calculate number of tide timesteps equal to hourly adjustment applied to each gage to transpose observed tide to BC compartment
				row_transpose = transposed_tide(nnn,2)/dttide
				use_row = tiderow + row_transpose
!>> use current tide row if transposed row is outside of imported timeseries
!>> first few hours of each yearly model run will simply repeat the observed tide for a number of timesteps equal to the tranpose time
				if (use_row < 1) then
					use_row = tiderow
				endif

! MP2023 added zw-04/06/2020
! last few hours of each yearly model run will simply repeat the observed tide for a number of timesteps equal to the tranpose time
                if(use_row >= (simdays*24/dttide)) then
                    use_row = tiderow
                endif

! jj is compartment number, jjk is boundary condition number
!				BCnosurge(jj,2)=TideData(use_row,nnn)
!				ES(jj,2)=BCnosurge(jj,2)+Surge(surgerow,jjk)
!>> !YW! Interpolate tide and surge between input data time step
                BCnosurge(jj,2) = BCnosurge(jj,1)                          
     &            +(TideData(use_row+1,nnn)-TideData(use_row,nnn))
     &            /lasttidestep
                BCsurge(jj,2) = BCsurge(jj,1)
     &            +(Surge(use_row+1,jjk)-Surge(use_row,jjk))
     &            /lasttidestep  

                Es(jj,2)=BCnosurge(jj,2)+BCsurge(jj,2)
			endif
		enddo
	  end do

!>> Update boundary condition water levels for compartments WITHOUT observed water level timeseries
	  do jjk=1,Mds 
		jj=KBC(jjk)
		do nnn=1,Mds-tidegages
			if (jj==weighted_tide(nnn,1)) then
!>> import distance weighting factors to weight closest two observed water level timeseries to the BC compartment
				EastWght = weighted_tide(nnn,2)
				EastComp = weighted_tide(nnn,3)
				WestWght = weighted_tide(nnn,4)
				WestComp = weighted_tide(nnn,5)

!>> Weight east and west observed timeseries from current timestep
				EastBC = BCnosurge(EastComp,2)
				WestBC = BCnosurge(WestComp,2)

!        BCnosurge(jj,2)=EastWght*EastBC+WestWght*WestBC  !zw added 04/06/2020
				Es(jj,2)=EastWght*EastBC+WestWght*WestBC
!     &                        +Surge(surgerow,jjk)
     &                    +EastWght*BCsurge(EastComp,2)               !YW! Apply the same weighting calculation as for tide
     &                    +WestWght*BCsurge(WestComp,2)               !YW! temporary for calibration
			endif
		enddo

		Eh(jj,2) = Es(jj,2)      !assign water elevel in marsh = openwater level at offshore boundary compartments
	  enddo

      return
	  end

!***********************End Subroutine for Tidal Boundary Condition*****************************



!***********************Old Subroutine for Tidal Boundary Condition*****************************
!
! kthr and kday now global parameters - no longer needed to be passed into subroutine
!      Subroutine TideBC(jj,kthr,kday)
!
!	use params
!
!
!
!c     These parameters should be moved to input to start at any time of year ** later
!********tide information, surges, periods and phase angles
!
!      shour=0.0
!	sday=0.0
!	f1=shour*3600.					! Daily phase
!	f2=(sday/28.)
!	f2=(f2-int(f2))*28.*consd		! Lunar phase
!	f3=(sday/365.25)*12.*30.*consd	! Gulf phase
!	t1=23.5*3600.					! Daily period
!	t2=672.*3600.					! Lunar period
!	t3=4320.*3600.					! Gulf period
!      tlag=0.							! Tide lag time along the eastern boundary in the Gulf
!	aset=0.01
!c *********** chenier average
!
!      ao=0.0615						! Reduced JAM July 09 0.1 --0.061
!	amax=0.47
!      tday= t/24/3600
!	tdayj=tday-(int(tday/365.25))*365.25   ! Cal julian day JAM Nov 2010
!	t2=tday*tday
!	t3=t2*tday
!	t4=t3*tday
!
!cc _________________ JAM Nov 2010 revised
!      g0 = -0.07						!m
!      g1 = -0.000233005
!      g2 = 0.000000326479
!      g3 = -0.000000000146371
!      g4 =  0.0000000000000299401
!      g5 = -0.00000000000000000283698
!      g6 =  0.000000000000000000000100482
!	thourj=(tdayj)*24
!	tt0 = 1
!	tt1 = thourj
!	tt2 = tt1*thourj
!	tt3 = tt2*thourj
!	tt4 = tt3*thourj
!	tt5 = tt4*thourj
!	tt6 = tt5*thourj
!	agulf2= g0+g1*tt1+g2*tt2+g3*tt3+g4*tt4+g5*tt5+g6*tt6	! JAM Nov 2010
!      cdir= wspd(kday)*wspd(kday)            !new wind variables used
!      cdir = windx(jj)**2 + windy(jj)**2       !cdir = windspd**2 = (sqrt(windx**2+windy**2))**2
!      if(Tempair(10,kday).lt.T7(kday))cdir=-cdir				! cdir = !- cos(pi*wd(kday)/180.)*wspd(kday)*wspd(kday)   ! setup due to wind shear
!      setup= (aset*cdir)/10.
!	if(setup.gt.0.5)setup=0.5								! not used here
!      if(setup.lt.-0.5)setup=-0.4								! not used here
!cc       cdir=0.0  !test oct 15 2013 AMc						! JKS 10/30/13
!	agulf2=agulf2+0.05+(aset*cdir)/3.*cos(pi*kday/8.)		! JAM 2010 /3.
!      aoo=ao*float(jj)/float(113)
!
!	amp=abs(sin(2*pi*t/(24*28*3600)))*(amax-ao)+aoo
!      TP1=23.2
! only use one harmonic tide equation now - THIS NEEDS TO BE UPDATED
!	if(jj.ne.248) then
!		Es(jj,2)=Max(ao,Abs(amax*sin(2*pi*t/(24*3600*28.))))			! JAM Dec 2010
!     &		*sin(2*pi*t/(TP1*3600))										! +ao*(1.-sin(2.*pi*t/(3600*24*365.25)))   !GOM seasonal effect. JAM Sept 08
!     &		+RSSS +agulf2+surge(jj,kday)*1.5							! Sealevel shift JAM Aug 09 + JAM Nov 2010 Surge Dec 7 2010
!      else
!					StageA=Qdiv(23,kday)*0.0000541+0.29					! Stage in Atchafalaya River at Morgan City
!					Es(jj,2)=StageA+0.5
!     &					*Max(ao,Abs(amax*sin(2*pi*t/(24*3600*28.))))	! JAM Dec 2010  Atch Diversion 107
!     &					*sin(2*pi*t/(TP1*3600))							! +ao*(1.-sin(2.*pi*t/(3600*24*365.25)))   !GOM seasonal effect. JAM Sept 08
!     &					+RSSS +agulf2+surge(jj,kday)*1.5				! Sealevel shift JAM Aug 09 + JAM Nov 2010 Surge Dec 7 2010
!	endif


!c***********************End Old Subroutine for Tidal Boundary Condition*****************************