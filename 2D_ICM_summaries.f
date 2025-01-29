!> @file
!> @brief This subroutine generates the data summaries used by other ICM routines.
!> @details   This subroutine calculates various summary statistics on hydrom model output
!> (e.g. mean salinity, variance of stage, etc.) that are required by the other ICM routines.
      
!> @ author Eric White - The Water Institute of the Gulf

!> @param[in]     M                          number of links
!> @param[in]     N                          number of compartments
!> @param[in]     sal_daily(366,N)           daily compartment salinity, calculated in hydrod (ppt)    
!> @param[in]     sal_daily_links(366,M)     daily link salinity, calculated in hydrod (ppt)
!> @param[in]     simdays                    number of days in the simulation year (either 365 or 366)
!> @param[in]     stage_daily(366,N)         daily compartment water stage, calculated in hydrod (m)    
!> @param[in]     tmp_daily(366,N)           daily compartment water temperature, calculated in hydrod (deg C)    
!> @param[in]     tmp_daily_links(366,M)     daily link water temperature, calculated in hydrod (deg C)

!> @param[out]    sal_ave(N)                 annual mean compartment salinity (ppt)
!> @param[out]    sal_summer(180,N)          daily compartment salinity for summer days (ppt)
!> @param[out]    sal_ave_summer(N)          mean summertime compartment salinity (ppt)  
!> @param[out]    sal_ave_links(M)           annual mean link salinity (ppt)
!> @param[out]    sal_ave_summer_links(M)    mean summertime link salinity (ppt)  
!> @param[out]    stage_ave(N)               annual mean compartment stage (m)
!> @param[out]    stage_summer(180,N)        daily compartment stage for summer days (m)
!> @param[out]    stage_ave_summer(N)        mean summertime compartment stage (m)  
!> @param[out]    stage_var_summer(N)        variance in summertime compartment stage (m)  
!> @param[out]    tmp_ave(N)                 annual mean compartment temperature (deg C)
!> @param[out]    tmp_summer(180,N)          daily compartment temperature for summer days (deg C)
!> @param[out]    tmp_ave_summer(N)          mean summertime compartment temperature (deg C)
!> @param[out]    tmp_ave_links(M)           annual mean link water temperature (deg C)
!> @param[out]    tmp_ave_summer_links(M)    mean summertime link water temperature (deg C)
!> @param[out]    sal_month_ave(12,N)        monthly average salinity for compartments
!> @param[out]    tmp_month_ave(12,N)        monthly average temperature for compartments
!> @param[out]    alg_month_ave(12,N)        monthly average algae for compartments
!> @param[out]    tss_month_ave(12,N)        monthly average TSS for for compartments
!> @param[out]    sal_month_ave_links(12,M)  monthly average salinity for links
!> @param[out]    tmp_month_ave_links(12,M)  monthly average temperature for links
!> @param[out]    alg_month_ave_links(12,M)  monthly average algae for links
!> @param[out]    tss_month_ave_links(12,M)  monthly average TSS for for links


!> @param         difmean(180,N)             temporary array used in calculating the variance in summertime stage
!> @param         ndays                      number of days in month
!> @param         summerstart                first day of summer - May 1 - updated for leap years
!> @param         summerend                  last day of summer - Aug 31 - updated for leap years
!> @param         summerlength               length (days) in summer
!> @param         salsum                     local value summing the daily salinity values of the month for compartments
!> @param         salsum_links               local value summing the daily salinity values of the month for links
!> @param         tempsum                    local value summing the daily temperature values of the month for compartments
!> @param         tempsum_links              local value summing the daily temperature values of the month for links
!> @param         tsssum                     local value summing the daily TSS values of the month for compartments
!> @param         tsssum_links               local value summing the daily TSS values of the month for links
!> @param         algsum                     local value summing the daily chlorophyll A values of the month for compartments
!> @param         algsum_links               local value summing the daily chlorophyll A values of the month for links
!> @param         stgsum1                    local value summing the monthly average stages for taking subset of annual
!> @param         nmonths1                   local value used for finding the average from a subset of months
!> @param         stgsum2                    local value summing the monthly average stages for taking subset of annual
!> @param         nmonths2                   local value used for finding the average from a subset of months
      
      
      subroutine ICM_Summaries
                
      use params

      implicit none
      integer :: j,k,jj,kk,kl,mn,DOYstart,DOYend,ndays
      integer :: summerstart,summerend,summerlength
      real(dp) :: salsum,salsum_links, sa, sae
      real :: tempsum,tknsum,tsssum
      real :: tempsum_links,tknsum_links,tsssum_links
      real :: stgsum,stgsum1,nmonths1,stgsum2,nmonths2

      write(1,*)'Preparing summary values from Hydro model output.'
      write(*,*)'Preparing summary values from Hydro model output.'      
!>@par General Structure of Subroutine Logic:

      sal_summer = 0.0
      stage_summer = 0.0
      tmp_summer = 0.0
      sal_summer_links = 0.0
      tmp_summer_links = 0.0
      trg_summer = 0.0
      
!>> Take subset of daily values for days that correspond to summertime.
      ! subset daily values to generate summertime timeseries
      ! summertime defined as May 1 through Aug 31, Jan 1 = day 0 
      summerstart = month_DOY(5)
      summerend = month_DOY(9)-1
      summerlength = summerend - summerstart
      
      do j=summerstart,summerend 
          jj = j-summerstart + 1    !convert day of year to day of summer for summer array !-EDW
          do k=1,N
              sal_summer(jj,k) = sal_daily(j,k)
              stage_summer(jj,k) = stage_daily(j,k)
              tmp_summer(jj,k) = tmp_daily(j,k)
              trg_summer(jj,k) = tidal_range_daily(j,k)
          enddo
          do kk=1,M
              sal_summer_links(jj,kk) = sal_daily_links(j,kk)
              tmp_summer_links(jj,kk) = tmp_daily_links(j,kk)
          enddo
      enddo


!>> Calculate monthly averages for use in HSI equations
      do mn=1,12
          DOYstart = month_DOY(mn)
          if (mn==12) then                !update last day of month for December (since their isn't a DOY_month(13))
              DOYend = month_DOY(12)+30
          else
              DOYend = month_DOY(mn+1) - 1
          endif
          ndays = DOYend-DOYstart+1
          
          do k=1,N
              stgsum = 0.0
              salsum = 0.0
              tempsum = 0.0
              tknsum = 0.0
              tsssum = 0.0
              do jj=DOYstart,DOYend
                  stgsum = stgsum + stage_daily(jj,k)
                  salsum = salsum + sal_daily(jj,k)
                  tempsum = tempsum + tmp_daily(jj,k)
                  tknsum = tknsum + tkn_daily(jj,k)
                  tsssum = tsssum + tss_daily(jj,k)
              enddo
              stg_month_ave(mn,k) = stgsum/ndays
              sal_month_ave(mn,k) = salsum/ndays
              tmp_month_ave(mn,k) = tempsum/ndays
              tkn_month_ave(mn,k) = tknsum/ndays
              tss_month_ave(mn,k) = tsssum/ndays
          enddo
          
          
          do kk=1,M
              salsum_links = 0.0
              tempsum_links = 0.0
              tknsum_links = 0.0
              tsssum_links = 0.0
              do jj=DOYstart,DOYend
                  salsum_links = salsum_links + sal_daily_links(jj,kk)
                  tempsum_links = tempsum_links + tmp_daily_links(jj,kk)
                  tknsum_links= tknsum_links + tkn_daily_links(jj,kk)
                  tsssum_links = tsssum_links + tss_daily_links(jj,kk)
              enddo
             
              sal_month_ave_links(mn,kk) = salsum_links/ndays
              tmp_month_ave_links(mn,kk) = tempsum_links/ndays
              tkn_month_ave_links(mn,kk) = tknsum_links/ndays
              tss_month_ave_links(mn,kk)= tsssum_links/ndays
          enddo
      enddo
      
!>> Calculate annual mean values.
      ! summing over dimension 1 of sal_daily sums the daily values at each compartment or link
      sal_ave = sum(sal_daily,DIM=1)/float(simdays)
      sal_ave_links = sum(sal_daily_links,DIM=1)/float(simdays)
      stage_ave = sum(stage_daily,DIM=1)/float(simdays)
      tmp_ave = sum(tmp_daily,DIM=1)/float(simdays)
      tmp_ave_links = sum(tmp_daily_links,DIM=1)/float(simdays)
      tss_ave = sum(tss_daily,DIM=1)/float(simdays)
      
!>> Calculate summertime mean values.
      ! summing over dimension 1 sums the daily values at each compartment
      ! denominator of 123. is the range of days which define 'summer' (243-120)
      sal_ave_summer = sum(sal_summer,DIM=1)/float(summerlength)
      sal_ave_summer_links = sum(sal_summer_links,DIM=1)/
     &                         float(summerlength)
      stage_ave_summer = sum(stage_summer,DIM=1)/float(summerlength)
      tmp_ave_summer = sum(tmp_summer,DIM=1)/float(summerlength)
      tmp_ave_summer_links = sum(tmp_summer_links,DIM=1)/
     &                         float(summerlength)
      trg_ave_summer = sum(trg_summer,DIM=1)/float(summerlength)

!>> Calculate mean stage values for bird HSIs.
      do k = 1,N
          stgsum1 = 0.0
          nmonths1 = 0.0
          stgsum2 = 0.0
          nmonths2 = 0.0

!>> Calculate September-March mean stage.
          do mn = 1,3
              stgsum1 = stgsum1 + stg_month_ave(mn,k)
              nmonths1 = nmonths1 + 1.0
          enddo
          do mn = 9,12
              stgsum1 = stgsum1 + stg_month_ave(mn,k)
              nmonths1 = nmonths1 + 1.0
          enddo
      
          sepmar_stage(k) = stgsum1/nmonths1

!>> Calculate October-April mean stage.
          do mn = 1,4
              stgsum2 = stgsum2 + stg_month_ave(mn,k)
              nmonths2 = nmonths1 + 1.0
          enddo
          do mn = 10,12
              stgsum2 = stgsum2 + stg_month_ave(mn,k)
              nmonths2 = nmonths2 + 1.0
          enddo
      
          octapr_stage(k) = stgsum2/nmonths2

      enddo
        
        
!>> Calculate variance in stage height over summer
      do j=1,summerlength
          do k = 1,N
              difmean(j,k) = (stage_summer(j,k)-stage_ave_summer(k))**2
          enddo
      enddo
      stage_var_summer = sum(difmean,DIM=1)/float(summerlength)

!>> Convert variance in stage to standard deviation
      do k = 1,N
          stage_stdv_summer(k) = stage_var_summer(k)**0.5
      end do
      
!>> Set water level variabilty equal to the Hydro-calculated standard deviation in daily water levels (2017 method - replaced for 2023 model)
      !stage_wlv_summer = stage_stdv_summer
      
!>> Set water level variability for use in ICM-LAVegMod from tidal range
!>> relationship to convert from tidal range to WLV used by Veg was determined by:
!>> - calculating WLV from hourly stage standard deviation from 2006-2018 at all CRMS sites (same method used to build CRMS input tables according to Scott D-S in 2017)
!>> - comparing to the mean tidal range at each CRMS site
!>> - R-squared from above regression = 0.7674        hourly stdv WLV = 0.2647*mean tidal range + 0.0659      Rsq=0.7674
      do k = 1,N
          stage_wlv_summer(k) = 0.2647*trg_ave_summer(k) + 0.0659
      end do
      
!>> Calculate variance in TSS       
      do j=1, simdays
          do k = 1,N
              difmean_tss(j,k) = (tss_daily(j,k)-tss_ave(k))**2
          enddo
      enddo
      tss_var_annual = sum(difmean_tss,DIM=1)/float(simdays)
      
      
!>> Determine annual max stage.
      !examining dimension 1 finds the maximum daily value at each compartment
      stage_max = maxval(stage_daily, DIM=1)
      
!>> Add error terms to summary values
      write(1,*) '  Adding error terms to summary output data.'
      write(*,*) '  Adding error terms to summary output data.'
      do k = 1,N
          sepmar_stage(k) = sepmar_stage(k) + stage_error
          octapr_stage(k) = octapr_stage(k) + stage_error
          tss_ave(k) = max(0.0,tss_ave(k) + tss_error)
          stage_ave(k) = stage_ave(k) + stage_error
          stage_max(k) = stage_max(k) + stage_error
!          stage_var_summer(k) = max(0.0,stage_var_summer(k) + stvar_error)
!          stage_stdv_summer(k) = max(0.0,stage_stdv_summer(k) + stvar_error)
          stage_wlv_summer(k) = max(0.0,stage_wlv_summer(k) + stvar_error)
          
          sa = sal_ave(k)
          if(sa < 1.0) then
              sae =  sal_0_1_error
          elseif(sa < 5.0) then
              sae =   sal_1_5_error
          elseif(sa < 20.0) then
              sae = sal_5_20_error
          else
              sae = sal_20_35_error
          endif
          sal_ave(k) = min(salmax,max(0.0,sa + sae))
          
          sa = sal_ave_summer(k) 
          if(sa < 1.0) then
              sae =  sal_0_1_error
          elseif(sa < 5.0) then
              sae =   sal_1_5_error
          elseif(sa < 20.0) then
              sae =   sal_5_20_error
          else
              sae =   sal_20_35_error
          endif
          sal_ave_summer(k) = max(0.0,sa + sae)
          
          do mn = 1,12
              stg_month_ave(mn,k) = stg_month_ave(mn,k) + stage_error
              tss_month_ave(mn,k)=max(0.0,tss_month_ave(mn,k)+tss_error)
              
              sa = sal_month_ave(mn,k)
              if(sa < 1.0) then
                  sae =  sal_0_1_error
              elseif(sa < 5.0) then
                  sae =   sal_1_5_error
              elseif(sa < 20.0) then
                  sae =   sal_5_20_error
              else
                  sae =   sal_20_35_error
              endif
              sal_month_ave(mn,k) = min(salmax,max(0.0,sa + sae))
              
          enddo
      enddo
      
      do kl = 1,M
          sa = sal_ave_links(kl)
          if(sa < 1.0) then
              sae =  sal_0_1_error
          elseif(sa < 5.0) then
              sae =   sal_1_5_error
          elseif(sa < 20.0) then
              sae =   sal_5_20_error
          else
              sae =   sal_20_35_error
          endif
          sal_ave_links(kl) = min(salmax,max(0.0,sa + sae))
          
          
          sa = sal_ave_summer_links(kl)
          if(sa < 1.0) then
              sae =  sal_0_1_error
          elseif(sa < 5.0) then
              sae =   sal_1_5_error
          elseif(sa < 20.0) then
              sae =   sal_5_20_error
          else
              sae =   sal_20_35_error
          endif
          sal_ave_summer_links(kl) = min(salmax,max(0.0,sa + sae))
          
          do mn = 1,12
              tss_month_ave_links(mn,kl) = max(0.0, 
     &                         tss_month_ave_links(mn,kl) +tss_error)
              
              sa = sal_month_ave_links(mn,kl)
              if(sa < 1.0) then
                  sae =  sal_0_1_error
              elseif(sa < 5.0) then
                  sae =   sal_1_5_error
              elseif(sa < 20.0) then
                  sae =   sal_5_20_error
              else
                  sae =   sal_20_35_error
              endif
              sal_month_ave_links(mn,kl) = min(salmax,max(0.0,sa + sae))
          
          enddo
          
      enddo
      
      return
      end 
            
                    