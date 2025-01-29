!> @file
!> @brief This subroutine generates the salinity threshold data summary used by the Wetland Morphology ICM routine.
!> @details   This subroutine calculates the 2-week average salinity from March through October and reports out the max value to be used compared to a salinity threshold in the marsh collapse routine.
      
      
!> @ author Eric White - The Water Institute of the Gulf

!> @param[in]     month_DOY(12)           array with the starting days of each month - updated for leap/non-leap years
!> @param[in]     M                       number of links
!> @param[in]     N                       number of compartments
!> @param[in]     sal_daily(366,N)        daily compartment salinity, calculated in hydrod (ppt)    
!> @param[in]     sal_daily_links(366,M)  daily link salinity, calculated in hydrod (ppt)
!> @param[in]     simdays                 number of days in the simulation year (either 365 or 366)

!> @param[out]    sal_2wk_ave_max(N)      maximum of 2-week averaged salinities in  compartments from March-October (ppt)

!> @param         firstday                first day of growing season (March 1)
!> @param         lastday                 last day of growing season (Oct 31)
!> @param         thresholdlength(366,N)  number of days to examine 2-wk average salinities for salinity threshold check
!> @param         sal_2wk_ave_numer(366,N) numerator for calculating moving window 2-wk average

      
      subroutine ICM_SalinityThreshold_WM

      use params

      implicit none
      integer :: dd,j,jj,k,kk,firstday,lastday,thresholdlength
      real(dp), dimension(:,:), allocatable :: sal_2wk_ave_numer
      real(dp), dimension(:,:), allocatable :: sal_2wk_ave
      real(dp), dimension(:,:), allocatable :: sal_2wk_ave_numer_links
      real(dp), dimension(:,:), allocatable :: sal_2wk_ave_links

      write(1,*)'Preparing 2-week salinity averages from March-October.' 
      write(*,*)'Preparing 2-week salinity averages from March-October.'      
!>@par General Structure of Subroutine Logic:

!>> Take subset of daily values for days that correspond to March-October.
      firstday = month_DOY(3)
      lastday = month_DOY(11)-1
      thresholdlength = lastday - firstday + 1
      
      allocate(sal_2wk_ave_numer(thresholdlength,N))
      allocate(sal_2wk_ave(thresholdlength,N))
      allocate(sal_2wk_ave_numer_links(thresholdlength,M))
      allocate(sal_2wk_ave_links(thresholdlength,M))
      
      ! Loop over days from March through end of October      
      do j=firstday,lastday
          jj = j-firstday+1    !convert day of year to day of March-Oct for array !-EDW

          ! Loop over compartments.          
          do k=1,N
              sal_2wk_ave_numer(jj,k) = 0.
              sal_2wk_ave(jj,k) = 0.
              
              ! Sum salinity values for current day and previous 2 weeks in compartment.
              do dd=0,13
                  sal_2wk_ave_numer(jj,k) = sal_2wk_ave_numer(jj,k) +
     &                                        sal_daily(j-dd,k)
              enddo
              
              ! Calculate average salinity for previous 2 weeks for day jj.              
              sal_2wk_ave(jj,k) = sal_2wk_ave_numer(jj,k)/14.
          enddo
          ! Loop over links.          
          do kk=1,M
              sal_2wk_ave_numer_links(jj,kk) = 0.
              sal_2wk_ave_links(jj,kk) = 0.
              
              ! Sum salinity values for current day and previous 2 weeks in compartment.
              do dd=0,13
                  sal_2wk_ave_numer_links(jj,kk) = 
     &                        sal_2wk_ave_numer_links(jj,kk) + 
     &                        sal_daily_links(j-dd,kk)
              enddo
              
              ! Calculate average salinity for previous 2 weeks for day jj.              
              sal_2wk_ave_links(jj,kk) = 
     &                        sal_2wk_ave_numer_links(jj,kk)/14.
          enddo
          
      enddo
!>> Determine max 2wk salinity value for each compartment and link from March-October.
      !examining dimension 1 finds the maximum value for each compartment or link
      sal_2wk_ave_max = maxval(sal_2wk_ave, DIM=1) 
      sal_2wk_ave_max_links = maxval(sal_2wk_ave_links, DIM=1) 
      
      deallocate(sal_2wk_ave_numer)
      deallocate(sal_2wk_ave)
      deallocate(sal_2wk_ave_numer_links)
      deallocate(sal_2wk_ave_links)

      return
      end 
            
            