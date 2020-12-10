!> @file
!> @brief This subroutine calculates waves from wind, depth, and fetch data.
!> @details The wave calculations in this subroutine were generated from the Young & Verhagen wave equations for non-dimensional frequency and energy.
!> In addition to the Young & Verhagen paper (Coastal Engineering, 29(1996)), fetch limitation conditions are determined using a methodology outlined
!> in the US Army Corps of Engineers' Coastal Engineering Manual. The CEM was also used for linear wave thoery equations.


      
!> @ author Eric White - The Water Institute of the Gulf

!> @param[in]		dtwind					timestep of wind observations
!> @param[in]     g                       gravitational acceleration = 9.81
!> @param[in]     j                       compartment number iterator from hydrod subroutine
!> @param[in]     Es(N,2)                 water surface elevation for model timestep at each compartment
!> @param[in]     Fetch(N,16)             fetch array for 16 different 22.5-deg sectors (column 1 corresponds to fetch length from 0.01-22.5 degrees)
!> @param[in]     Bed(N)                  bed elevation of compartment
!> @param[in]     S(N,2)                  salinity in compartment
!> @param[in]     windx(N)                wind X-vector for compartments at model timestep (m/s)
!> @param[in]     windy(N)                wind speed Y-vector for compartments at model timestep (m/s)

!> @param[out]	group_vel(N,2)			array of wave group velocity for compartments
!> @param[out]	Hs(N,2)					array of significant wave height for compartments (m)
!> @param[out]	Uorb(N,2)				array of orbital velocity values (at bed level) for compartments (m/s)
!> @param[out]	wave_energy(N,2)		array of wave energy values for compartments
!> @param[out]    wave_frequency(N,2)     array of wave frequency for compartments
!> @param[out]    wavelength(N,2)         array of wavelength values for compartments
!> @param[out]	wave_period(N,2)		array of wave period for compartments


!> @param         depth_value             water depth in compartment
!> @param         e                       natural log base = 2.71828
!> @param         energy_limit            asymptotic limit used in Young & Verhagen wave equation
!> @param         ff                      local do loop iterator
!> @param         fetch_value             fetch length (m) for compartment in direction wind is coming from
!> @param         frequency_limit         asymptotic limit used in Young & Verhagen wave equation
!> @param         YV_n                    constant used in Young & Verhagen wave equation
!> @param         YV_m                    constant used in Young & Verhagen wave equation
!> @param         wind_dir_rads           direction wind is blowing TO (in radians)    
        


      subroutine waves_YV(j)

      			
	use params

      implicit none
      integer :: ff,j,fetch_lookup,kday
	real :: e, wind_spd,rhow
	real :: windx_value,windy_value,wind_dir_rads,wind_dir_degs
	real :: dim_factor,fetch_value,depth_value,duration_needed
	real :: drag_coeff,wind_fric_spd,fetch_equiv
      real :: energy_limit,frequency_limit,energy_nondim,fetch_nondim
      real :: depth_nondim,frequency_nondim,waven
      real :: YV_A1,YV_A2,YV_B1,YV_B2,YV_m,YV_n,Lo,d_L,t1,t2




!>@par General Structure of Subroutine Logic:
!>> Set various constants used in Young & Verhagen wave equations
      e = 2.71828
      energy_limit = 0.00364
      frequency_limit = 0.133
      YV_n = 1.74
      YV_m = -0.37
	rhow = 1000.*(1+S(j,1)/1000.)	!adjust density for salinity

!>> Calculate wind speed and direction.      
      windx_value = max(0.01,windx(j))
      windy_value = max(0.01,windy(j))
      wind_spd = sqrt(windx_value**2 + windy_value**2)
      wind_dir_rads = atan(windy_value/windx_value)
      if (day > 360) then
          write(*,*) 'j:', j
          write(*,*) 'windx:', windx(j)
          write(*,*) 'windx:', windy(j)
          write(*,*) 'wind_spd:', wind_spd
          write(*,*) 'wind_dir_rads:', wind_dir_rads
      endif
!>> Convert wind direction to degrees and lookup fetch from 16 different directions (each separated by 22.5 degrees)      
      wind_dir_degs = wind_dir_rads*180./pi
      if (wind_dir_degs <= 0.) then
          wind_dir_degs = wind_dir_degs + 360.
      endif
      fetch_lookup = max(1, int(floor(wind_dir_degs/22.5)+1) )
      if (day > 360) then
          write(*,*) 'fetch_lookup:', fetch_lookup
      endif
!>> Pull appropriate fetch for compartment and wind direction from lookup array      
      fetch_value = max(0.1,Fetch(j,fetch_lookup))
      if (day > 360) then
          write(*,*) 'fetch_value:', fetch_value
      endif
!>> Calculate dimensionless factor for use in Young & Verhagen wave equations
      dim_factor  = g/(wind_spd**2)
 
!>> Calculate water depth for timestep      
      depth_value = max(Es(j,2) - Bed(j),0.001)
      
!>> Calculate dimensionless fetch and depth values      
      fetch_nondim = fetch_value*dim_factor
      depth_nondim = depth_value*dim_factor

!>> Calculate friction-adjusted wind speed
      drag_coeff = 0.001*(1.1+0.035*wind_spd)
      wind_fric_spd = wind_spd*sqrt(drag_coeff)
      
!>> Calculate wind duration needed (in seconds) for waves to become fetch-limited over current fetch
      duration_needed = (77.23*fetch_value**0.67/
     &                     ((wind_spd**0.34)*(g**0.33)))

     
!!!OLD DURATION LIMITED EQUIVALENT FETCH ROUTINE -START
!!!
!!!!>> Calculate equivalent fetch (from USACE CEM - eqs. II-2-36 & II-2.38) if wind interval is less than required duration
!!!      if( (duration_needed/60.) > (dtwind*60.) ) then
!!!          fetch_equiv=0.00523*sqrt(g*wind_fric_spd*duration_needed**3.0)
!!!          fetch_value = fetch_equiv
!!!!>> Update non-dimensional fetch to be used in Y&V equations
!!!          fetch_nondim = fetch_value*dim_factor
!!!      endif
!!!
!!!OLD DURATION LIMITED EQUIVALENT FETCH ROUTINE -END


!>> Check if duration of wind speed will meet fetch-limited conditions (duration_needed is in seconds, dtwind is in hours)
      if( (duration_needed/60.) > (dtwind*60.) ) then
!>> If duration-limited, calculate wave period for duration-limited, shallow water conditions (from Shore Protection Manual (1984) eq. 3-41)
          wave_period(j,1) = (wind_fric_spd/g)
     &                   *(g*dtwind*3600./(537.*wind_fric_spd))**(3./7.)
!>> If duration-limited, calculate significant wave height for duration-limited, shallow water conditions (from Shore Protection Manual (1984) eq. 3-39)
          t1 = 0.530*(g*depth_value/wind_fric_spd**2.)**0.75
          t2 = 0.00565*(g*fetch_value/wind_fric_spd**2.)**0.5
          Hs(j,1) = 0.283*(wind_fric_spd/g)*tanh(t1)*tanh(t2/tanh(t1))

!>> If fetch-limited, calculate Young & Verhagen energy and frequency          
      else
!>> If fetch-limited, A1 and B1 in Young & Verhagen (eq. 20 & eq. 21)
          YV_A1 = (0.292**(1/YV_n))*(depth_nondim**(1.3/YV_n))
          YV_B1 = (0.00004396**(1/YV_n))*(fetch_nondim**(1/YV_n))
        
!>> If fetch-limited, A2 and B2 in Young & Verhagen (eq. 23 & eq. 24)
          YV_A2 = (1.505**(1/YV_m))*(depth_nondim**(-0.375/YV_m))
          YV_B2 = (16.391**(1/YV_m))*(fetch_nondim**(-0.27/YV_m))
        
!>> If fetch-limited, calculate Young & Verhagen non-dimensional energy (eq. 19)
          energy_nondim = energy_limit*(tanh(YV_A1)*
     &                    tanh(YV_B1/tanh(YV_A1)))**YV_n
        
!>> If fetch-limited, calculate Young & Verhagen non-dimensional frequency (eq. 22)
          frequency_nondim  = frequency_limit*(tanh(YV_A2)*
     &                         tanh(YV_B2/tanh(YV_A2)))**YV_m

!>> If fetch-limited, calculate wave energy, frequency, significant wave height, and wave period from Young & Verhagen non-dimensional frequency & dimensionaless factor     
          wave_energy(j,1) = rhow*g*energy_nondim/dim_factor**2
          wave_frequency(j,1) = frequency_nondim*dim_factor*wind_spd
	    Hs(j,1)=4.*sqrt(wave_energy(j,1))
          wave_period(j,1) = 1/wave_frequency(j,1)
      endif
      
      
!>> Calculate wavelength and group velocity approximations from linear wave theory (USACE CEM Part 2.1)
      Lo = g*(wave_period(j,1)**2)/(2*pi)
      wavelength(j,1) = Lo*sqrt(tanh(2*pi*depth_value/Lo))
      d_L = depth_value/wavelength(j,1)
      Waven = (1+(4*pi*d_L)/sinh(4*pi*d_L))/2
      group_vel(j,1) = Waven*wavelength(j,1)/wave_period(j,1)      
      
!>> Calculate orbital velocity at bottom of water column
	Uorb(j,1) = g*Hs(j,1)*wave_period(j,1)/(2*wavelength(j,1)*
     &			cosh(2*pi*depth_value/wavelength(j,1)))
      

      return
      end

