!> @file
!> @brief This subroutine calculates phytoplankton photosynthesis rates.
!> @details This subroutine calculates photosynthesis rates, as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @author  Eric White - The Water Institute of the Gulf

!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     Es(N,2)         current water surface elevation in compartment (m)
!> @param[in]     Bed(N)          bed elevation of compartment (m)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]     css(N,2,4)      array with current suspended sediment concentrations for 4 particle classes (mg/L)
!> @param[in]		saltox			salinity concentration where algal growth is halved
!> @param[in]     S(N,2)          array with current salinity concentration in compartment (mg/L)
!> @param[in]     totalsolrad()   array with daily solar radiation (langly/day)

!> @param[out]     muph           phytoplankton photosynthesis rate (1/day)

!> @param 		e				constant e
!> @param 		perOW			percent open water
!> @param 		dd				water depth
!> @param 		daylight		hours of daylight in day
!> @param 		solrad			incident solar radiation during daylight (langly/day)
!> @param 		dailysolrad		total daily incident solar radiation (langly/day)
!> @param 		alg				algae concentration at previous timestep
!> @param 		det				detritus concentration at previous timestep
!> @param 		din				dissolved inorganic N concentration at previous timestep
!> @param 		tip				total inorganic P concentration at previous timestep
!> @param 		iss				inorganic suspended sediment concentration at previous timestep
!> @param 		khn				half saturation concentration for algal uptake of N (mg/L)
!> @param 		khp				half saturation concentration for algal uptake of P (mg/L)				

!> @param 		fnfix			fraction of algae that fixes N
!> @param 		kd				adsorption distribution coefficient (L/kg)
!> @param 		fpp				fraction of TIP in particulate form
!> @param 		dip_up			portion of TIP that is dissolved and available for uptake (1-fpp)*TIP
!> @param 		kphot20			reaction rate @ 20 deg-C for photosynthesis
!> @param 		thetaphot		reaction rate coefficient for photosnythesis
!> @param 		kphot			temperature dependent reaction rate for photosynthesis
!> @param 		limnp1			nutrient limitation - intermdediate calculation
!> @param 		limnp2			nutrient limitation - intermdediate calculation
!> @param 		limnp			nutrient limitation factor
!> @param 		keb				background light extinction coefficient due to water and coloer (1/m)
!> @param 		alphi			light extinction coefficient due to suspend solids (L/mg/m)
!> @param 		alphao			light extinction coefficient due to detritus (L/mg/m)
!> @param 		alphap			light extinction coefficient due to Chlorophyll A
!> @param 		alphapn			light extinction coefficient due to Chlorophyll A
!> @param 		ke				light extinction coefficient (1/m)
!> @param 		klp				light where phytoplankton growth is optimal
!> @param 		par0ave			average light intensity over daylight hours at water surface
!> @param 		limlight1		light limitation - intermdediate calculation
!> @param 		limlight2		light limitation - intermdediate calculation
!> @param 		limlight3		light limitation - intermdediate calculation
!> @param 		limlight		light limitation factor
!> @param 		limsal			salinity limitation factor

      Subroutine photo(j)
	

      use params
      
      implicit none
      integer :: j
      real :: e,perOW,dd,daylight,solrad,dailysolrad
      real :: alg,det,din,tip,iss
      real :: khn,khp,fnfix,kd,kdc,dip_up,kphot
      real :: limnp1,limnp2,limnp
      real :: keb,alphai,alphao,alphap,alphapn,ke,klp
      real :: par0ave,limlight1,limlight2,limlight3,limlight
      real :: limsal    

!>>  photosynthesis rate routine (eq. 27 of 2012 Master Plan Appendix D-1)
	e = 2.71828	

!>> percent open water of compartment (used in light limitation factor)
	perOW = Apctwater(j)

!>> current depth in compartment
!	dd = max(Es(j,2) - Bed(j),0.01)
	dd = max(Es(j,2) - Bed(j),dry_threshold)

!>> fraction of daylight hours
	daylight = 0.5

! average daily solar radiation at NO Int Airport ~ 4000 Wh/m2 (NCDC National Solar Radiation Database - MSY is NSRD station # 722310) 
      dailysolrad=4000.0/24.0
!>> incident solar radiation @ surface (langly/day) (0.484 converts W/m2 to langley/day)
	solrad = dailysolrad/0.484

!>> previous time step WQ concentrations
	alg = chem(j,8,1)		
	det = chem(j,9,1)
	din = chem(j,3,1)
      tip = chem(j,5,1)
	iss = css(j,2,1) + css(j,2,2) + css(j,2,3) + css(j,2,4)

!>> temperature-dependent max photosynthetic rate
	kphot = kphot20*thetaphot**(Tempw(j,2)-20.)
      
!>> half saturation concentration  for algal uptake of N (mg/L)
	khn = 0.02

!>> half saturation concentration for algal uptake of P (mg/L)
	khp = 0.001

!>> calculate salinity-dependent algae coefficients: phythoplankton preference for NH4 uptake and fraction of algae that is nitrogen-fixing
	if (S(j,2) > 2.0) then
		fnfix = 0.0
	else
		fnfix = 0.1
	endif 

!>> convert total inorganic phosphorus to dissolved inorganic phosphorus readily available for uptake (eq 37 in App D1)
      kd = 500.0              ! TIP sorption distribution coefficient (L/kg)
      kdc = kd/1000000.0
      fpp = kdc*iss/(1.0+kdc*iss)
      dip_up = (1 - fpp)*tip      

!>> nutrient limitation factor
	limnp1 = (1-fnfix)*din/(khn+din)+fnfix
	limnp2 = dip_up/(khp+dip_up)
	limnp = min(limnp1,limnp2)

!>> light extinction factor
	keb = 0.3
	alphai = 0.052
	alphao = 0.174
	alphap = 8.8		!units = 8.8 L/mg/m 
	alphapn = 5.412		!units = 5.412 (L/mg)^0.667/m.
	ke = keb + alphai*iss + alphao*det + alphap*alg
     &		+ alphapn*alg**(2./3.)
      
!>> photosynthetically available radiation (langly/day) - PAR for optimal growth = 300.0
	klp = 300.0
	par0ave = 0.47*solrad

!>> light limitation factor
	limlight1 = (-par0ave*e**(-ke*dd))/klp
	limlight2 = -par0ave/klp
	limlight3 = e**limlight1 - e**limlight2
      !	limlight = 2.718*daylight*limlight3/(ke*dd)
      limlight = limlight3
      
!>> salinity limitation factor
	limsal = saltox**2/(saltox**2+S(j,2)**2)

!>> photosynthesis rate
	muph = kphot*limnp*limlight*limsal
      
      if(isNan(ke)) then
          write(*,*) alphai,alphao,alphap,alphapn,iss,det,alg
          write(*,*) 'NaN encountered in Photosynthesis calculations'
          stop!pause
      endif
      
      
      return
      end