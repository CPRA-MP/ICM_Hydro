
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

                  ES(jj,2)=BCnosurge(jj,2)+BCsurge(jj,2)
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
				ES(jj,2)=EastWght*EastBC+WestWght*WestBC
!     &                        +Surge(surgerow,jjk)
     &                    +EastWght*BCsurge(EastComp,2)               !YW! Apply the same weighting calculation as for tide
     &                    +WestWght*BCsurge(WestComp,2)               !YW! temporary for calibration
			endif
		enddo
	enddo

      return
	end

!***********************End Subroutine for Tidal Boundary Condition*****************************