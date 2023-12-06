!> @file
!> @brief This subroutine calculates change in nitrate-nitrite concentrations.
!> @details This subroutine calculates change in nitrate-nitrite concentrations as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @author  Eric White - The Water Institute of the Gulf


!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]     S(N,2)          array with current salinity concentration in compartment (mg/L)
!> @param[in]     muph            phytoplankton photosynthesis rate (1/day)
!> @param[in]		kdenit20		reaction rate @ 20 deg-C for denitrification
!> @param[in]		thetadenit		reaction rate coefficient for denitrification
!> @param[in]		knit20			reaction rate @ 20 deg-C for nitrification
!> @param[in]		thetanit		reaction rate coefficient for nitrification

!> @param[out]     DChemSum       change in concentration of water quality constituent during timestep

!> @param 		no3				NO3 concentration at previous timestep
!> @param 		nh4				NH4 concentration at previous timestep
!> @param 		alg				algae concentration at previous timestep
!> @param 		knit			temperature dependent reaction rate for nitrification
!> @param         knitf           knit filtered with max and min thresholds on nitrification rate
!> @param 		kdenit			temperature dependent reaction rate for denitrification
!> @param 		khn				half saturation concentration for algal uptake of N (mg/L)
!> @param 		fnfix			fraction of algae that fixes N
!> @param         pap             phytoplankton preference for ammonium uptake
!> @param         rca             carbon-to-chlorophyll A stoichiometric mass ratio
!> @param         rna             nitrogen-to-chlorphyll A stoichiometric mass ratio



!      Subroutine NO3(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine dNO3(DChemSum,ichem,j)
	!JAM Chem # 1

	use params

      implicit none
      integer :: ichem,j
      real :: dd,e
      real :: no3,nh4,alg
      real :: knit,knitf
      real :: kdenit
      real :: khn,pap,fnfix
      real :: rca,rna
      real :: dChemSUM 
      
!>> NO3 routine (eq. 19 of 2012 Master Plan Appendix D-1)
	e = 2.71828	

      !>> previous time step WQ concentrations
	no3 = chem(j,1,1)
	nh4 = chem(j,2,1)
	alg = chem(j,8,1)		

!>> current depth in compartment
!	dd = max(Es(j,2) - Bed(j),0.01)
	dd = max(Es(j,2) - Bed(j),dry_threshold)
            
!>> temperature-dependent nitrification rate coefficient 
	knit20 = max(min(knit20num/dd,knitmax),knitmin)
      knit = knit20*thetanit**(Tempw(j,2)-20.)
	      
!>> temperature-dependent denitrifcation rate coefficient
	kdenit = kdenit20*thetadenit**(Tempw(j,2)-20.)

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

!>> nitrogen-chlorophyll A stoichiometric mass ratio
	rna = 0.176*rca

!>> change in nitrate/nitrite concentration
	dChemSUM = knit*nh4-(kdenit*no3/dd)-rna*(1-pap)*muph*(1-fnfix)*alg
      
      return
      end
      
!	rca=0.004
!	fixN=0.00000001
!	fixNO3=rca*Chem(j,8,1)*fixN
!      Tcorn=1.085**(Tempw(j,2)-20.)						! Temperature Correction  Nitrification
!	Tcord=1.05**(Tempw(j,2)-20.)						! Temperature Correction Denitrification
!	ALAG=0.2
!	PP=1-0.2*cos(2*PI*(Float(kday)/365.25-ALAG)) 
!	fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)				! Assumes proportionate uptake of NO3 and NH4
!	FNoC=1./5.681
!	cxx=1/3600/24/1000000.
!	fcN=1.
!      FSEASON2=1.0+0.2*cos(2*PI*(Float(kday)/365.25+0.25))			! JAM Jan 09 11
!      frr=0.85
!	fmm=1.
!	fgrow=1.0
!
!cc	if(j.eq.38.or.j.eq.39.or.j.eq.40.or.j.eq.41.or.j.eq.47) then	! JKS Oct 16 2013 temporary!	! Need to increase NO3 uptake in these cells
!cc		frr=0.10													! JKS Oct 16 2013 temporary!
!cc		fmm=6.														! JKS Oct 16 2013 temporary!
!cc		fgrow=2.													! JKS Oct 16 2013 temporary!
!cc		Tcorn=1.085**(Tempw(j,2)+1.-20.)							! JKS Oct 16 2013 temporary!	! Temperature Correction Nitrification
!cc		Tcord=1.05**(Tempw(j,2)+1.-20.)								! JKS Oct 16 2013 temporary!	! Temperature Correction deNitrification
!cc    endif															! JKS Oct 16 2013 temporary!
!
!	fmarsh=fmm+0.5*(Ahf(j)/As(j,1))*floodf(j)						! Loss from marsh uptake JAM March 9 2011                                  !
!	                                      
!      DChemSUM=-(1.-flood(j))*Ahydro(j)*cxx*
!     &	Max(0.,(Rain(kday,jrain(j))-PET(kday)))*UplandNP(j,1)*
!     &	FSEASON2*frr												! JAM Jan 09 11  !ungauged upland contribution
!      DChemSUM=DChemSUM												! +fixNO3*Tcorn	!Nitrogen Fixation ~Chla*kfx
!      dnx=1.0
!      denit(j,2)=decay(1,1)*Chem(j,1,1)*dnx/ds(j,1)*Tcorn*PP*fmm		! JAM April 2, 2011
!      DChemSUM=DChemSUM+ denit(j,2)									! Decay(1,1)*Chem(j,1,1)*dnx/ds(j,1)*Tcorn*PP*fmm        !denitrification
!      DChemSUM=DChemSUM+decay(1,2)*Chem(j,2,1)*dnx/ds(j,1)*Tcorn		! Increase due to nitrification NH4--> NO2-->NO3
!      DChemSUM=DChemSUM+decay(ichem,8)*GrowAlgae(j,8,8)*				! Loss due to uptake by live algae
!     &	Chem(j,8,1)*fnodin*FNoC*fmarsh*fgrow
!	return
!	end

c______________________________________________________________________________________
cc
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


c***********************End Subroutine for Chemical chemical NO3+NO2****************************