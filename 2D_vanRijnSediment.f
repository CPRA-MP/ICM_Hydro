!> @file
!> @brief This subroutine calculates sediment distribution for 4 different sediment classes (sand, silt, unflocculated clay, and flocculated clay).
!> @details The sediment resuspension and deposition calculations in this subroutine were developed
!> during CPRA subtask 4.2 of the Model Improvment Plan for the 2017 Master Plan. This routine was 
!> formulated largely based on equations developed by van Rijn (2004) as well as some other supporting
!> literature. Refer to the Subtask 4.2 memorandum for full documentation of these equations.
      
!> @ author Eric White - The Water Institute of the Gulf

      
!> @param[in]     Ahydro(j)               area of compartment that is water (km**2)
!> @param[in]     Percent(j)              portion of water in compartment that is open water
!> @param[in]     g                       gravitational acceleration = 9.81
!> @param[in]     icc()                   compartment-link connectivity matrix
!> @param[in]     j                       compartment number iterator from hydrod subroutine
!> @param[in]     Es(N,2)                 water surface elevation for model timestep at each compartment
!> @param[in]     Bed(N)                  bed elevation of compartment
!> @param[in]     D50(class)			    D50 grain size for sediment class (m)
!> @param[in]     flocA                   clay flocculation parameter (McAnally A)
!> @param[in]     flocB                   clay flocculation parameter (McAnally B)
!> @param[in]     flocM                   clay flocculation exponent (McAnally m)
!> @param[in]     flocN                   clay flocculation exponent (McAnally n)
!> @param[in]     flocC1                  clay concentration where floculation settling begins (McAnally C1)(kg/m**3)
!> @param[in]     flocC3                  clay concentration where floculation settling becominbegins to hinder settling (McAnally C3) (kg/m**3)
!> @param[in]     Csalmax                 upper limit to salinity impact on floc (salinity concentration - ppt)
!> @param[in]     Pflocmax                upper limit on portion of fines available to floc (ratio - unitless)
!> @param[in]     D90x                    representative D90/D50 ratio
!> @param[in]     ka                      wind-driven currents coefficient, ranges from 0.023-0.032
!> @param[in]     cf                      bed shear stress coefficient, ranges from 0.001-0.003
!> @param[in]		group_vel(N,2)			array of wave group velocity for compartments
!> @param[in]     rhoSed(class)           density of sediment particles in each particle size class (kg/m**3)
!> @param[in]     S(N,2)                  salinity in compartment
!> @param[in]		SedDens(class)			representative particle density for each sediment class (kg/m**3)
!> @param[in]		Tempw(N,2)				array of water temperature for compartments
!> @param[in]	    Uorb(N,2)				array of orbital velocity values (at bed level) for compartments (m/s)
!> @param[in]		wave_frequency(N,2)     array of wave frequency for compartments
!> @param[in]		wavelength(N,2)         array of wavelength values for compartments
!> @param[in]		wave_period(N,2)		array of wave period for compartments
!> @param[in]     windx(N)                wind X-vector for compartments at model timestep (m/s)
!> @param[in]     windy(N)                wind speed Y-vector for compartments at model timestep (m/s)
!> @param[in]     Qsum_abs(N)             absolute value of flows into and out of compartment
!> @param[in]     Qsum_in(N)              sum of flows into compartment
!> @param[in]     Qsum_out(N)             sum of flows out of compartment
!> @param[in]     MEESedRate(class)       sediment flux (g/sec) from marsh edge erosion (negative is sediment into compartment)
!> @param[in]     QSmarsh(class)          sediment flux (g/sec) from open-water-to-marsh exchange flow (negative is sediment into compartment)
!> @param[in]     QSsum(class)            sediment flux (g/sec) from all link flow - this term for sand (class=1) is only flow into compartment  (negative is sediment into compartment)
!> @param[in]     QStrib(class)           sediment flux (g/sec) from tributary flows (negative is sediment into compartment)
!> @param[in]     QSdiv(class)            sediment flux (g/sec) from diversion flows (negative is sediment into compartment)
      
!> @param         ddy_1                   water depth in compartment at previous timestep      
!> @param         ddy_2                   water depth in compartment at current timestep
!> @param         e                       natural log base = 2.71828
!> @param         QSand                   Sum of all suspended sand flux into compartment (from links, diversions, and tribs) - calculated from current timestep flow and previous timestep CSSvRs
!> @param         vRQs                    van Rijn sand deposition flux (kg/s-m2)
!> @param         vRQse                   van Rijn equilibrium sand deposition flux (kg/s) 
!> @param         QSandTribLeave          sand load leaving via negative tributary flow
!> @param         depo_avail              sediment flux avaialble for deposition based on suspended concentration in compartment (g/sec)
!> @param         resus_avail             sediment flux available for resuspension based on amount of erodible bed remaining (g/sec)
      
      
!> @param[out]    CSSvRs(N,2)             suspended sand concentration for compartment from van Rijn calcs
!> @param[out]    kinvisc(j)		        kinematic viscosity in compartment (m**2/s)
!> @param[out]    Dr(class)				dimensionless particle diameter for sediment class
!> @param[out]    Dgr(class)		    	grain size diameter for sediment class(m)
!> @param[out]    SedDens(class)		    particle density for sediment class
!> @param[out]    SGsed(class)			specific gravity of sediment for sediment class
!> @param[out]    velF(class)			    velocity factor for sediment class
!> @param[out]    velset(N,class)			settling velocity for sediment class in compartment (m/sec)
!> @param[out]    resuspension(N,class)   resuspensed sediment (in g/sec) for sediment class in compartment
!> @param[out]    deposition(N,class)     deposited sediment (in g/sec) for sediment class in compartment
!> @param[out]    SedAccumRate(class)     net sediment accumulation rate (in g/seC) for sediment class in compartment - positve value is deposition on bed, negative value is resuspension into water column

      subroutine vanRijnSediment(j,kday)

      use params
	  implicit none
	  integer :: j,jn,jt,k,kk,kday
	
      real :: ddy_1,ddy_2,OWwidth,OWArea,Uflows
      real :: windx_value,windy_value,wind_spd,Uwind,Ubed,Tbed
      real :: aa,bb,cc,dd,ee,ff,rhow,Dr,velF,YY,ZZ
      real :: T_dimless,Tcrit_dimless,Dstar,Me,sMe
      real :: Ue,Ucr,Uterm,Ucrc,Ucrw
      real :: Csand,Csilt,CClayAll,CClayFloc,CClayUnfloc
      real :: CClayAllH,CClayFlocH,CClayUnflocH
      real :: QSand,QSandInSum,vRqs,vRQse,equilflux
      real :: depo_avail,resus_avail,sed_avail
      real :: depo_avail_rate,depo_settling
      real :: eroded_amount
!>@par General Structure of Subroutine Logic:

!>> Calculate kinematic viscosity from temperature.(eq 7d in methodology memo)
	  kinvisc(j) = 0.00000179/
     &				(1+0.03369*Tempw(j,1)+0.000221*Tempw(j,1)**2)

!>> Calculate water depths for current and previous timesteps - depth is not allowed to be less than 0.0001 in celldQ, so these values will always be non-zero
      ddy_1 = Es(j,1) - Bed(j)
      ddy_2 = Es(j,2) - Bed(j)
     
!>> Adjust water density for salinity
      rhow = 1000.*(1+S(j,1)/1000.)

      OWArea = As(j,1)


!>> Calculate portion of fines available for flocculation as function of salinity (eq. 8 in methodology memo) 
	  if (S(j,1) < Csalmax) then
		  Pfloc = Pflocmax*S(j,1)/Csalmax
	  else
		  Pfloc = Pflocmax
	  endif

!>> Read in concentrations of different sediment classes from previous timestep
      CSand = CSS(j,1,1)
      CSilt = CSS(j,1,2)
      CClayAll = CSS(j,1,3) + CSS(j,1,4)
      CClayAllH = CSSH(j,1,3) + CSSH(j,1,4)
!> Partition total clay CSS into floc and unfloc CSS's based on salinity
      CClayFloc = Pfloc*CClayAll
      CClayUnfloc = CClayAll - CClayFloc
      CClayFlocH = Pfloc*CClayAllH
      CClayUnflocH = CClayAllH - CClayFlocH
!>> Update Floc and Unfloc clay concentrations based on salinity in compartment
      CSS(j,2,3) = min(max(cssmin,CClayUnfloc),cssmax)
      CSS(j,2,4) = min(max(cssmin,CClayFloc),cssmax)
      CSSH(j,2,3) = min(max(cssmin,CClayUnflocH),cssmax)
      CSSH(j,2,4) = min(max(cssmin,CClayFlocH),cssmax)

!>> Calculate settling velocities for different sediment classes
!>> particle class size (k): 1=sand,2=silt,3=unfloc clay, 4=floc clay
      do k=1,4 
		  Dgr(k) = D50(k)         !Dgr(k) is grain size - can set to non-D50 values here in the future
          D90(k) = D90x*D50(k)
		  SGsed(k) = Specg        !ability to set different sediment specific gravities - currently use one value
          rhoSed(k) = SGsed(k)*rhow 

!>> Calculate settling velocity for sand particles (eq. 7a & 7b in methodology  memo)		
		  if (k == 1) then
		 	  Dr = Dgr(k)*sqrt((g*(SGsed(k)-1.)/kinvisc(j)**2.))
			  velF = max(sqrt((36./Dr**3.) + 2./3.) - sqrt(36./(Dr**3.)),0.0)
			  velset(j,k) = velF*sqrt((SGsed(k)-1.)*Dgr(k)*g)
		  else
!>> Calculate settling velocity for silt & clay particles (eq. 6 in methodology memo)
			  velset(j,k) = (g*(SGsed(k)-1.0)*D50(k)**2.0)/(18.0*kinvisc(j))
		  endif
	  enddo

!>> Update settling velocity for flocculated clay particles (sediment class (k) = 4) 
! BUG - flocC1 and flocC3 unit is kg/m3 in the inputs whereas CSS is in g/m3!!
! BUG - flocA/flocB/flocN/flocM in eq.9 of methodology memo are for unit of kg/m3!!
! BUG fix - need to convert CClayFloc from g/m3 to kg/m3 in the following codes - zw 1/4/2024
      CClayFloc = CClayFloc/1000.0  !convert CClayFloc from g/m3 to kg/m3
	  if (CClayFloc < flocC1) then
		  flocset = velset(j,4)
	  elseif (CClayFloc < flocC3) then
		  flocset = flocA*(CClayFloc**flocN)/(CClayFloc**2+flocB**2)**flocM
	  else
		  flocset = 0.000001
	  endif
	  velset(j,4) = flocset
	
!>> Calculate wind-driven currents.
      windx_value = windx(j)
      windy_value = windy(j)
      wind_spd = sqrt(windx_value**2 + windy_value**2)
      Uwind = ka(j)*wind_spd

!>> Determine cumulative flow into and out of compartment
!>> Convert cumulative flow into magnitude of velocities into and out of compartment (assume half of flows leave compartment, half enter - hence divide by 2)
      OWwidth = sqrt(As(j,1))        ! assume open water portion of compartment is square
!      Uflows = Qsum_abs(j)/(2*ddy_2*OWwidth) 
      Uflows = Uwind + ave_vel(j)  !Uflows should exclude tributary flows
      
!>> Calculate total velocity at bed of open water compartment (winds, flows, waves).
      Ubed = Uwind + Uorb(j,1) + ave_vel(j) !Uflows
!      write(*,*) j,Uflows,ave_vel(j),min_vel(j),max_vel(j)

!>> Calculate bed shear stress from velocity at bed.
      Tbed = cf(j)*rhow*Ubed**2
      !write(*,*) j,Uwind,Uflows,Uorb(j,1),Tbed

!>> Loop over 4 sediment classes and determine maximum possible sediment accumulation rates based on flow conditions and an infinite sediment source    
      do k=1,4
          if (k /= 1) then
!>> -- Determine deposition fluxes (g/sec) for non-sand particles (via Krone methodology - eq. 10 in methodology memo)	
              if (Tbed<Tcrit(k)) then
                  deposition(j,k) = velset(j,k)*(1.0-Tbed/Tcrit(k))
     &                                *Css(j,1,k)*OWArea
                  resuspension(j,k) = 0.0
!>> -- Determine resuspension fluxes for non-sand particles (g/sec) (eq. 11 in methodology memo is for kg/sec - convert to g/sec).
              else
                  deposition(j,k) = 0.0
                  resuspension(j,k)=(((Tbed/Tcrit(k))-1.0)**sedn(j))
     &                                *sedcalib(j)*OWArea*1000.0     !sedcalib unit should be kg/m2/s
                         
              endif
!>> -- calculate maximum possible non-sand sediment accumulation rate (g/sec) based on sediment fluxes (positive rate is deposition on open water bed)
              SedAccumRate(k) = deposition(j,k)
     &                         - resuspension(j,k)
!     &                         - MEESedRate(k)
!     &                         - QSmarsh(k)
!     &                         - QSsum(k)
!     &                         - QStrib(k)
!     &                         - QSdiv(k)

!>> -- Determine resuspension and deposition fluxes for sand particles (eqs. 32-40 in methodology memo)      
          else	
!>> -- set van Rijn constants based on D50 size
              if (D50(k) < 0.0005) then
                  aa = 0.19
                  bb = 0.1
                  cc = 0.24
                  dd = 0.66
                  ee = 0.33
                  ff = 0.33
              elseif (D50(k) <= 0.002) then
                  aa = 8.5
                  bb = 0.6
!                  cc = 8.5
!BUG fix - 0.95 in van Rijn 2013 as cited in methodology memo - zw 1/4/2024
                  cc = 0.95  
                  dd = 0.57
                  ee = 0.43
                  ff = 0.14
              else
				  write(1,*) '*************ERROR**************'
				  write(1,*) 'D50 for sand is greater than 2 mm.'
				  write(1,*) 'The van Rijn equation should not be used.'
				  write(1,*) 'Correct grain size definitions and re-run.' 
                            
                  write(*,*) '*************ERROR**************'
				  write(*,*) 'D50 for sand is greater than 2 mm.'
				  write(*,*) 'The van Rijn equation should not be used.'
				  write(*,*) 'Correct grain size definitions and re-run.' 
				  stop !pause
                  
              endif
                  
!>> -- calculate van Rijn terms
              if (ddy_2>dry_threshold) then
                  YY = log10( 12*ddy_2/(3*D90(k)) )
                  ZZ = (SGsed(k)-1)*g
                  Ucrc = aa*(D50(k)**bb)*YY          
                  Ucrw = cc*ZZ**dd*(D50(k)**ee)*(wave_period(j,1)**ff)
                  if ((Uflows+Uorb(j,1))>0) then
                      Uterm = Uflows/(Uflows+Uorb(j,1))
                  else
                      Uterm=0
                  endif
!zw 1/23/2024 Not sure why Ucrw used to calculate critical velocity. In the 1D codes, only Ucrc is used for sand
!                  Ucr=Ucrc  !simialr to 1D sand 
                  Ucr = Uterm*Ucrc + (1-Uterm)*Ucrw
                  Ue = Uflows + 0.4*Uorb(j,1)   !0.4 for irregular waves, 0.8 for regular waves
!              Me = abs((Ue-Ucr)/sqrt(g*D50(k)*(SGsed(k)-1))) ! need absolute value since this value is raised to 2.4 power - removes sign of Ucr
!              sMe = Me/((Ue-Ucr)/sqrt(g*D50(k)*(SGsed(k)-1)))
! zw 1/22/2024 revised to avoid dividing by 0 condition when Ue=Ucr
                  Me = (Ue-Ucr)/sqrt(g*D50(k)*(SGsed(k)-1))
                  sMe=1.0
                  if (Me<0) sMe=-1.0
!>> -- dimensionless shear stress
!              T_dimless = Tcrit(k)/((rhoSed(k)-rhow)*g*D50(k))
!  BUG fix - Tbed should be used instead of Tcrit - zw 1/4/2024
                  T_dimless = Tbed/((rhoSed(k)-rhow)*g*D50(k))
!>> -- reference particle diameter, Dstar
                  Dstar = D50(k)*(g*(SGsed(k)-1)/(kinvisc(j)**2.))**(1./3.)
!>> -- dimensionless CRITICAL shear stress for sand
                  if (Dstar<4) then
                      Tcrit_dimless = 0.115*(Dstar**(-0.5))
                  else
                      Tcrit_dimless = 0.14*(Dstar**(-0.64))
                  endif
!>> -- van Rijn suspended sediment flux (kg/m/sec) - represents the sediment carrying capacity of the flowing water - can be negative, this would represent periods of deposition b/c flow does not have any carrying capacity for sand                 
                  if (T_dimless <= Tcrit_dimless) then
                      vRqs = 0.0
                  else
                      vRqs = sMe*alphaSed(j)*rhoSed(k)*Uflows*
     &                     D50(k)*(abs(Me)**2.4)*(Dstar**(-0.6))  !qs in 1D code qrijn07 (eroflux), unit=kg/m/s
                  endif
! !>> -- van Rijn equilibrium sediment flux out of compartment, Qse        
!              vRQse = vRqs*OWwidth        !equilibrium sand load (in kg/sec)
! !>> -- van Rijn suspended sand concentration in compartment for current timestep - CSS value of flows leaving compartment - in mg/L, same as CSS(:)
! !>> -- if effective velocity is less than critical velocity, van Rijn sediment flux calculation is negative - deposition will occur and van Rijn suspended sand concentration should be set to zero
!              if (Qsum_out(j) /= 0.0) then
!                  CSSvRs(j,2) = min(CSSresusOff(k), max(
!     &                           1000.0*vRQse*dt/(As(j,1)*ddy_2),0.0))
!              else
!                  CSSvRs(j,2) = 0.0
!              endif
!>> -- var Rijn equilibrium concentration (kg/m3 -> g/m3) - zw 1/4/2024
			      CSSvRs(j,2) = vRqs/(Uflows*ddy_2)*1000.0 
!>> -- calculate maximum possible sand sediment accumulation rate (g/sec) based on sediment fluxes (positive rate is deposition on open water bed)              
!              SedAccumRate(k) = (CSS(j,1,k)*As(j,1)*ddy_1/dt)
!     &                         - (CSSvRs(j,2)*As(j,1)*ddy_2/dt)   
!     &                         - MEESedRate(k)
!     &                         - QSmarsh(k)
!     &                         - QSsum(k)
!     &                         - QStrib(k)
!     &                         - QSdiv(k)
! Source/sink term for sand (g/sec) = alpha*settling velocity *(ceq-c) where ceq=CSSvRs(j,2), c=CSS(j,1,k) - zw 1/4/2024
! alpha is non equilibrium adaptation coefficient (constant or varying)
! The nonequilibrium adaptation coefficient is assigned different values by different
! investigators in their studies. Han et al. (1980) and Wu and Li (1992) used a=1 for strong erosion,
! a = 0:25 for strong deposition and a = 0.5 for weak erosion and deposition. Yang (1998) used a
! very small value 0.001. Thus, the nonequilibrium adaptation coefficient is defined as
! a user-defined parameter in the model as a compartment attribute in cells.csv
                  SedAccumRate(k) = (CSS(j,1,k)-CSSvRs(j,2))*velset(j,k)
     &                              *As(j,1)*adaption_coeff(j)   
              else  !depth considered as dry condition, only deposition
                  SedAccumRate(k) = CSS(j,1,k)*velset(j,k)*As(j,1)  
              endif
          endif

!>> Reduce maximum sediment accumulation rates (g/sec) based on available sediment in water column or in erodible bed
!>> -- if deposition occurs, determine maximum sediment availble to deposit (g/sec)
          if (SedAccumRate(k) > 0) then
!>> -- Calculate sediment suspended during during timestep (g)
              sed_avail = CSS(j,1,k)*As(j,1)*ddy_1
!>> -- Calculate deposition rate if all sediment were to settle out (g/sec)
              depo_avail = max(sed_avail,0.0)/dt
!>> -- Calculate deposition rate based on current settling velocity (g/sec)
              depo_settling = CSS(j,1,k)*velset(j,k)*As(j,1)
!>> -- Determine deposition rate 
              depo_avail_rate = min(depo_settling,depo_avail)
!>> -- Set sediment accumulation rate to smaller of two values and multiply by flag for compartments where settling is turned off              
              SedAccumRate(k) = depo_on_off(j)*
     &                            min(SedAccumRate(k),depo_avail_rate)

!>> -- if resuspension occurs, determine if enough material remains in the erodible bed - material remains in bed as long as cumulative sediment accumulation value is greater than the mass that would be removed if entire erodible bed were eroded
          else
!>> -- determine amount of sediment in erodible bed available for resuspension
              if (k == 1) then
                  eroded_amount = Sandacc(j,1)
              elseif (k == 2) then
                  eroded_amount = Siltacc(j,1)
              elseif (k == 3) then
                  eroded_amount = Clayacc(j,1)/2.0
              else
                  eroded_amount = Clayacc(j,1)/2.0
              endif
              
              resus_avail = (eroded_amount + erBedAvail(j,k))*As(j,1)/dt           !if sacc > 0 , total sediment available bed + deposited amount

              resus_avail = max(0.0,resus_avail) ! if resus_avail < 0, bed has already been eroded

!>> -- limit resupsension rate to only erode the amount of erodible bed available for the sediment class
              SedAccumRate(k) = depo_on_off(j)*max(SedAccumRate(k),-resus_avail) ! SedAccumRate and available sediment for resuspension are negative here, since resuspension is occurring)
              
!>> -- if CSS is greater than threshold value (in RuncontrolR.dat) do not allow any resuspension
              if (css(j,1,k) >= CSSresusOff(k)) then
                  SedAccumRate(k) = 0.0
              endif
              
          endif
          
          
      enddo        
     	
      return
      end

      
