!> @file
!> @brief This subroutine generates the tree establishment conditions array for use in the Vegetation model.
!> @details   This subroutine determines if there is any period of time during the year in which
!> tree establishment conditions are met at each 500-m grid cell.
!> A tree establishment condition is set to 1 if at any point from March 1 throuh Aug 16, a grid cell has two weeks
!> of dry land (depth <= 0) followed by 2 weeks in which the water depth is no deeper than 10 cm.

      
      
!> @ author Eric White - The Water Institute of the Gulf

!> @param[in]     Bed(N)                  bed elevation of compartment
!> @param[in]     Es(N)                   daily compartment water level
!> @param[in]     grid_lookup_500m(n_500m_cells,16)       lookup table for 500 m grid - matches compartment and link to grid cell
!> @param[in]     month_DOY(12)           array with the starting days of each month - updated for leap/non-leap years
!> @param[in]     N                       number of compartments
!> @param[in]     n_500m_cells            number of 500-m grid cells in Veg grid
!> @param[in]     simdays                 number of days in the simulation year (either 365 or 366)
!> @param[in]     stage_daily(366,N)      daily compartment water stage, calculated in hydrod (m) 

!> @param[out]    tree_est(n_500m_cells)  tree establishment conditions array - integer, 1 or 0      

!> @param         comp                    compartment number for 500 m grid cell (from lookup table)
!> @param         dd                      number of days in moving window
!> @param         dryfuture_flag          flag (1 or 0) to determine if future 14 days have less than 10 cm of ponding
!> @param         drypast_flag            flag (1 or 0) to determine if last 14 days were dry
!> @param         firstday                first day of tree establishment window
!> @param         grid_dep_daily(n_500m_cells, simdays) array with daily water depth for each 500-m grid cell
!> @param         j                       iterator over all 365/365 simulation days
!> @param         jj                      iterator over days during tree establishment window
!> @param         jjj                     day of establishment window converted from simulation day
!> @param         k                       iterator over number of 500 m grid cells
!> @param         lastday                 last day of tree establishment window
!> @param         thresholdlength         number of days to analyze for tree establishment conditions
!> @param         tree_est_flag           combined flags to see if both conditions are met (dry past AND dry future)

      subroutine ICM_TreeConditions_Veg
      
      use params      

      implicit none
      integer :: j,k,jj,jjj,dd,comp
      integer :: firstday,lastday,thresholdlength
      real(sp), dimension(:,:), allocatable :: grid_dep_daily
      integer, dimension(:), allocatable :: drypast_flag
      integer, dimension(:), allocatable :: dryfuture_flag
      integer, dimension(:), allocatable :: tree_est_flag

    
!>@par General Structure of Subroutine Logic:
!>> Set first and last day of tree establishment window
      firstday = month_DOY(3)
      lastday = month_DOY(8)+16
      thresholdlength = lastday - firstday + 1

!>> Allocate temporary array to be of length equal to number of grid cells - these are deallocated at end of this subroutine      
      allocate(grid_dep_daily(n_500m_cells,simdays))
      allocate(drypast_flag(thresholdlength))
      allocate(dryfuture_flag(thresholdlength))
      allocate(tree_est_flag(thresholdlength))

      
      
!>> Loop through each 500 m grid cell.
      do k=1,n_500m_cells
!>> Map compartment depth values to grid cells for each day          
          comp = grid_lookup_500m(k,2)    ! second column of grid_lookup_array is compartment # for grid cell
          do j = 1,simdays
              grid_dep_daily(k,j)=stage_daily(j,comp)-land_elev_500m(k)
          enddo
          
!>> Loop through days at each grid cell and determine tree establishment criteria is met.
          do jj = firstday,lastday          
              jjj=jj-firstday+1   !convert day of year to day of tree establishment window array
!>> -- Initially set day's drypast and dryfuture flags to values of 1.
              drypast_flag(jjj) = 1
              dryfuture_flag(jjj) = 1

!>> Loop over past two weeks and determine if any past day is wet.
              do dd=0,13
                  if (grid_dep_daily(k,jj-dd) <= -0.30) then
                      drypast_flag(jjj) = drypast_flag(jjj)*1
                  else
!>> -- If any day of the past 2 weeks is wet, drypast_flag is set to 0.
                      drypast_flag(jjj) = drypast_flag(jjj)*0
                  endif
!>> Loop over next two weeks and determine if any future day is wet.                  
                  if (grid_dep_daily(k,jj+dd) <= -0.20) then
                      dryfuture_flag(jjj) = dryfuture_flag(jjj)*1
                  else
!>> -- If any day of the next 2 weeks is flooded by more than 10 cm, dryfuture_flag is set to 0.                  
                      dryfuture_flag(jjj) = dryfuture_flag(jjj)*0
                  endif
                  
              enddo
              tree_est_flag(jjj) = dryfuture_flag(jjj)*drypast_flag(jjj)
          
          enddo
          
!>> Loop over cell's timeseries of flags and set equal to 1 if any daily flags equal 1
          tree_est(k) = 0
          do jj = firstday,lastday
              jjj = jj-firstday+1
              if (tree_est_flag(jjj) > 0) then
                  tree_est(k) = 1
              endif
          enddo
          
          
      enddo

      
      deallocate(grid_dep_daily)
      deallocate(drypast_flag)
      deallocate(dryfuture_flag)
      deallocate(tree_est_flag)

      return
      end