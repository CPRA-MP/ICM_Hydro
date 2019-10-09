!> @file
!> @brief This subroutine calculates change in ammonium concentration.
!> @details This subroutine calculates change in ammonium concentration as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @author  Eric White - The Water Institute of the Gulf


!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]     S(N,2)          array with current salinity concentration in compartment (mg/L)
!> @param[in]     muph            phytoplankton photosynthesis rate (1/day)
!> @param[in]		kresp20	    	reaction rate @ 20 deg-C for phytoplankton respiration
!> @param[in]		thetaresp		reaction rate coefficient for phytoplankton respiration

!> @param[out]     DChemSum       change in concentration of water quality constituent during timestep

!> @param 		no3				NO3 concentration at previous timestep
!> @param 		nh4				NH4 concentration at previous timestep
!> @param 		alg				algae concentration at previous timestep
!> @param 		din				DIN concentration at previous timestep
!> @param 		kdon20			reaction rate @ 20 deg-C for DON hydrolysis
!> @param 		thetadon		reaction rate coefficient for DON hydrolysis
!> @param 		kdon			temperature dependent reaction rate for DON hydrolysis
!> @param 		kresp			temperature dependent reaction rate for phytoplankton respiration
!> @param 		khn				half saturation concentration for algal uptake of N (mg/L)
!> @param 		khp				half saturation concentration for algal uptake of P (mg/L)
!> @param 		fnfix			fraction of algae that fixes N
!> @param         pap             phytoplankton preference for ammonium uptake
!> @param         rca             carbon-to-chlorophyll A stoichiometric mass ratio
!> @param         rna             nitrogen-to-chlorphyll A stoichiometric mass ratio


      
      !      Subroutine NH4(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine dNH4(DChemSum,ichem,j)
	!JAM Oct 2010 Chem #2

      use params
      
      implicit none
      integer :: ichem,j
      real :: dd,e
      real :: no3,nh4,alg,don,kdon
      real :: kresp,knit,knitf
      real :: khn,khp,pap,fnfix
      real :: rca,rna
      real :: dChemSUM 
      
!>> NH4 routine (eq. 20 of 2012 Master Plan Appendix D-1)
	e = 2.71828	

!>> previous time step WQ concentrations
	no3 = chem(j,1,1)
	nh4 = chem(j,2,1)
	alg = chem(j,8,1)		
	don = chem(j,10,1)

!>> temperature-dependent DON hydrolysis rate coefficient 
	kdon = kdon20*thetadon**(Tempw(j,2)-20.)
	
!>> temperature-dependent phytoplankton respiration rate coefficient
	kresp = kresp20*thetaresp**(Tempw(j,2)-20.)
      
!>> current depth in compartment
	dd = max(Es(j,2) - Bed(j),0.01)
            
!>> temperature-dependent nitrification rate coefficient 
	knit20 = max(min(knit20num/dd,knitmax),knitmin)
      knit = knit20*thetanit**(Tempw(j,2)-20.)
	      

!>> half saturation concentration  for algal uptake of N (mg/L)
	khn = 0.02

!>> calculate salinity-dependent fraction of algae that is nitrogen-fixing
	if (S(j,2) > 2.0) then
		fnfix = 0.0
	else
		fnfix = 0.1
	endif 

!>> calculate phythoplankton preference for N uptake
      pap = (nh4*no3/((khn+nh4)*(khn+no3))) 
     &			+ (nh4*khn/((nh4+no3+0.001)*(khn+no3)))

!>> carbon-to-chlorophyll ratio
	rca = 75.0

!>> nitrotgen-chlorophyll A stoichometric mass ratio
	rna = 0.176*rca

!>> change in ammonium concentration
	dChemSUM = kdon*don - knit*nh4 + 
     &                (rna*kresp - rna*pap*muph*(1-fnfix))*alg
      
      return
      end
      
!	rca=0.004
!	fixN=0.000001
!	FNoC=1./5.7
!	fixNO3=rca*Chem(j,8,1)*fixN
!	Tcorn=1.085**(Tempw(j,2)-20.)					! Temperature Correction  Nitrification
!      frr=0.85
!	fmm=1
!	if(j.eq.243) then						! Changed from 	"if(j.eq.71.or.j.eq.39) then" Oct 16 2013 JKS	! JAM March 2001	if(j.eq.67.or.j.eq.71.or.j.eq.39) then
!		fmm=2.
!		frr=0.5
!      endif
!
!	fmarsh=fmm+0.5*(Ahf(j)/As(j,1))*floodf(j)  
!	FNoC=1./5.681
!	fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
!	DChemSUM=DChemSUM+decay(ichem,ichem)*Chem(j,ichem,1)*Tcorn	! Decrease due to nitrification NH4--> NO2-->NO3
!     DChemSUM=DChemSUM+decay(ichem,10)*Chem(j,10,1)*Tcorn		! Anaerobic decomposition of organic material   
!    DChemSUM=DChemSUM+decay(ichem,8)*GrowAlgae(j,8,8)*			! Loss due to uptake by live algae
!     &	Chem(j,8,1)*(1-fnodin)*FNoC*Tcorn*fmarsh
	
!	return
!	end

c______________________________________________________________________________________
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
cc   1	2	3	4	 5	 6	7	  8	     9	   10	11	12	13	  14 ***********
c*********************************************************************************
cc-0.0815	0.1	0	0	 0	 0	0    -0.5	 0	   0	 0	 0	 0	  0		!1 NO2+NO3
cc  0 -0.21	0	0    0	 0	0    -0.5	 0	   0.01	 0	 0	 0	  0		!2 NH4
c    0	0	0	0	 0	 0	0	  0	     0	   0	 0	 0	 0	  0		!3 DIN
cc 0-0.00001 0 -0.000001 0 0	0     1	-0.00001   0	 0	 0	 0	  0		!4 ON
cc 0	    0	0	0 -0.002 0	0     -1	  0    0     0	 0	 0	  0		!5 TP
c   0	   0    0	0	0	0 -0.005 0.265 -0.055  0	 0	 0	 0	  0		!8 LiveAlgae
c   0	   0    0	0	0	0 -0.005-0.055  0.055 -0.01-0.005 0	 0	  0		!9 DeadAlgae
cc  0 -0.01	0	0	0	0	0    0.0	0.005	0	 0	 0	 0	  0		!10 DON
cc  0	   0	0	0	0	0	0   0.01	  0	   0  0.01 -0.01 0	  0		!11 DOP
cc  0	   0	0	0	0	0	0	  0	      0	   0  -0.01	0.01 0	  0		!12 SRP
cc  0	   0	0	0	0	0	0	0.06667	  0	   0	0	0  0.01	  0		!13 Chla
cc  0	   0	0	0	0	0	0	  1	      0	-0.01	 0	 0	 0	  0		!14 POP
c**********************************************************************************
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
c______________________________________________________________________________________

c***********************End Subroutine for chemical NH4*****************************************