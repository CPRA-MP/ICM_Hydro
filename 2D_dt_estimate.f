!> @file
!> @brief This is the subrouine to estimate daily dt in 2D ICM-hydro.
!> @details This is the subrouine estimate daily dt in 2D ICM-hydro 
!> based on potential highest water level within a day from offshore.
!
      Subroutine dt_estimate(isimday,dtmax)

      use params

      implicit none

      integer :: i,isimday
      real :: zfull, yfull, twave, dtmax

      dtmax=dt
      twave=dt
!>> Link updates of dtmax
      do i=1,M
!>> Link type 1 = rectangular channels
!>> Link type 3 = rectangular channels with a lock control structure
!>> Link type 11 = marsh 'composite flow' channel
!>> Link type 12 = rectangular channel with elevations not updated by ICM/morph code
          if ((linkt(i) == 1) .or. (linkt(i) == 3) .or.
     &              (linkt(i) == 11) .or. (linkt(i) == 12)) then

          !>> determine max wl elevation
              zfull = daily_maxWL(isimday) !max(ES(jus(i),2),ES(jds(i),2))

          !>> determine full depth
              yfull = zfull-Latr1(i)

          !>> determine time takes for a dynamic wave to travel the length of the conduit
              if( yfull>0) then
                  twave = Latr3(i)/sqrt(g*yfull)
              endif

!>> Link type 6 = rectangular culverts/bridges
          elseif (linkt(i) == 6)then

          !>> determine bank elevation
              zfull = Latr2(i)

          !>> determine full depth
              yfull = zfull-Latr1(i)

          !>> determine time takes for a dynamic wave to travel the length of the conduit
              if( yfull>0) then
                  twave = Latr3(i)/sqrt(g*yfull)
              endif

!>> Link type 8 = marsh over land channels
          elseif (linkt(i) == 8) then

          !>> determine bank elevation
              zfull = daily_maxWL(isimday) !max(ES(jus(i),2),ES(jds(i),2))

          !>> determine full depth
              yfull = zfull-Latr1(i)

          !>> determine time takes for a dynamic wave to travel the length of the conduit
              if( yfull>0) then
                  twave = Latr3(i)/sqrt(g*yfull)
              endif

          else
              twave=dt
          endif
          dtmax=min(dtmax,twave)
      enddo
!>> End link timestep estimates

      return
      end

!***********************End Subroutine hydrodynamic*********************************************
