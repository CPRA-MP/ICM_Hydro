!> @file
!> @brief This subroutine calculates change in live algae concentrations.
!> @details This subroutine calculates change in live algae concentrations as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     Es(N,2)         current water surface elevation in compartment (m)
!> @param[in]     Bed(N)          bed elevation of compartment (m)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]     muph            phytoplankton photosynthesis rate (1/day)
!> @param[in]		kphy20			reaction rate @ 20 deg-C for phytoplankton mortality
!> @param[in]		thetaphy		reaction rate coefficient for phytoplankton mortality
!> @param[in] 	kresp20	    	reaction rate @ 20 deg-C for phytoplankton respiration
!> @param[in]		thetaresp		reaction rate coefficient for phytoplankton respiration

!> @param[out]     DChemSum       change in concentration of water quality constituent during timestep

!> @param 		kphy			temperature dependent reaction rate for phytoplankton mortality
!> @param 		kresp			temperature dependent reaction rate for phytoplankton respiration
	
	Subroutine dLivA(DChemSum,ichem,j)
!	dLivA(DChemSum,ichem,mex,j,k,kday)
      use params
      
	implicit none

	integer :: j, ichem
	real :: dd
	real :: alg
	real :: kresp
      real :: kphy
	real :: va
      real :: DChemSUM

!>> live algae calculation (eq. 25 of 2012 Master Plan Appendix D-1)
!>> current depth in compartment
!	dd = Es(j,2) - Bed(j)
	dd = max(Es(j,2) - Bed(j),0.01)	 !zw 4/28/2015 be consistent with NO3 and NH4

!>> previous time step WQ concentrations
	alg = chem(j,8,1)
      
!>> temperature-dependent phytoplankton respiration rate coefficient
	kresp = kresp20*thetaresp**(Tempw(j,2)-20.)

!>> temperature-dependent phytoplankton mortality rate coefficient
	kphy = kphy20*thetaphy**(Tempw(j,2)-20.)      

!>> algae settling rate in compartments (m/day)
	va = 0.01

!>> change in live algae concentration (mg/L)
	DChemSUM = (muph-kresp-kphy-va/dd)*alg

	return
	end


!      Subroutine LivA(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
!      Subroutine LivA(DChemSum,ichem,mex,j,k,kday)			! Live Algae
!	
!      use params
!      
!      !JAM Oct 2010 Chem #8
!
!
!	xsqrt=2.01
!	Tcorn=1.0667**(Tempw(j,2)-20.) 
!      baseLG=1.0								!0.85
!	FLG=0.25
!	ALAG=0.25
!	KnSal=1.
!	FLGT=baseLG-FLG*cos(2*PI*(Float(kday)/365.25+ALAG))		! JAM Feb 2011   !Light 
!	PP=1-0.2*cos(2*PI*(Float(kday)/365.25+ALAG))			! Relative Photo period 
!      CSSsum = 0.0
!      do sedclass=1,4
!          CSSsum = CSSsum + CSS(j,1,sedclass)
!      enddo
!      FCSS= 1.- CSSsum/(KnSS+CSSsum)      !replaced CSS(j,1) with CSSsum calculated in do loop now
!	Fsalinity= 1.-S(j,1)/(KnSal/2.+S(j,1))
!	GM6=(1-Chem(j,8,1)/(0.00251+Chem(j,8,1)))
!      fdead=1.0
!	GM7=1+Chem(j,8,1)/(0.00251+Chem(j,8,1))		!Limits
!
!      if(j.eq.14)then						! Changed from 2012 SMP Cell 6 to Lower B&B Cell 14 Oct 16 2013 JKS
!		fdead=0.9
!		gm6=4.   
!	endif								! JAM April 8, 2011
!
!      if(j.eq.22)gm6=2.					! Changed from 2012 SMP Cell 5 to Lower B&B Cell 22 Oct 16 2013 JKS		! JAM April 8, 2011
!      if(j.eq.33)then						! Changed from 2012 SMP Cell 1 to Lower B&B Cell 33 Oct 16 2013 JKS		! JAM April 9, 2011
!      fdead=2.75
!	endif
!      
!      if(j.eq.51)then						! Changed from 2012 SMP Cell 67 to Lower B&B Cell 51 Oct 16 2013 JKS		! JAM April 9, 2011
!	 gm6=2.2							! 0.000000001  April 22, 2011
!       fdead=0.6
!	endif
!
!     	Fnutrient=MIN(1.,Chem(j,3,1)/(1.1*KnN/1000000.+Chem(j,3,1)),	! 2/1000000-->1.1/1000000
!     &	Chem(j,12,1)/(KnP/1000000./2.+Chem(j,12,1)))				! *Chem(j,8,1)*1000./0.5	! Chem(j,12 instead of chem(j,5  JAM April 2011 BC effect
!
!cccccccccccccccc      
!	GrowAlgae(j,ichem,ichem)=decay(ichem,ichem)*FLGT*xSQRT*FCSS		!(abs(1.-
!    &	*FSalinity*Fnutrient*Tcorn/20.*gm6							!/24	!Respiration and elimination JAM April 18 2011
!	Respir=decay(ichem,7)*Chem(j,ichem,1)*Tcorn 
!	DChemSUM=DChemSUM+GrowAlgae(j,ichem,ichem)*Chem(j,ichem,1)      ! Increase due to uptake
!     &	+decay(ichem,9)*Chem(j,ichem,1)*Tcorn/1.5*gm7*fdead+Respir	! * Chem(j,9,1)/			! net loss ~ death rate & elimination rate  = - DA
!	
!	return
!	end     
	
!********apply growth restraint to Dead Algae based on Redfield
!c______________________________________________________________________________________
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
cc   1	2	3	4	 5	 6	7	  8	     9	   10	11	12	13	  14 ***********
c*********************************************************************************
cc-0.0815	0.1	0	0	 0	 0	0    -0.5	 0	   0	 0	 0	 0	  0		!1 NO2+NO3
cc  0 -0.21	0	0    0	 0	0    -0.5	 0	   0.01	 0	 0	 0	  0		!2 NH4
c    0	0	0	0	 0	 0	0	  0	     0	   0	 0	 0	 0	  0		!3 DIN
cc 0-0.00001 0 -0.000001 0 0	0     1	-0.00001   0	 0	 0	 0	  0		!4 ON
cc 0	    0	0	0 -0.002 0	0     -1	  0    0     0	 0	 0	  0		!5 TP
c   0	   0    0	0	0	0 -0.005* 0.265* -0.055*  0	 0	 0	 0	  0		!8 LiveAlgae
c   0	   0    0	0	0	0 -0.005-0.055  0.055 -0.01-0.005 0	 0	  0		!9 DeadAlgae
cc  0 -0.01	0	0	0	0	0    0.0	0.005	0	 0	 0	 0	  0		!10 DON
cc  0	   0	0	0	0	0	0   0.01	  0	   0  0.01 -0.01 0	  0		!11 DOP
cc  0	   0	0	0	0	0	0	  0	      0	   0  -0.01	0.01 0	  0		!12 SRP
cc  0	   0	0	0	0	0	0	0.06667	  0	   0	0	0  0.01	  0		!13 Chla
cc  0	   0	0	0	0	0	0	  1	      0	-0.01	 0	 0	 0	  0		!14 POP
c**********************************************************************************
cc! NO3	NH4	DIN	ON	 TP	TOC	DO	  LA-C	DA-C	DON	DOP	SRP	ChLa POP
c______________________________________________________________________________________

c***********************End Subroutine for Live Algae ******************************************