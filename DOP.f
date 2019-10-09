!> @file
!> @brief This subroutine calculates change in dissolved organic phosphorus concentration.
!> @details This subroutine calculates change in dissolved organic phosphorus concentration as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @author  Eric White - The Water Institute of the Gulf


!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]		kdet20			reaction rate @ 20 deg-C for detritus dissolution
!> @param[in]		thetadet		reaction rate coefficient for detritus dissolution

!> @param[out]     DChemSum       change in concentration of water quality constituent during timestep

!> @param 		det				detritus concentration at previous timestep
!> @param 		don				DON concentration at previous timestep
!> @param 		kdon20			reaction rate @ 20 deg-C for DON hydrolysis
!> @param 		thetadon		reaction rate coefficient for DON hydrolysis
!> @param 		kdon			temperature dependent reaction rate for DON hydrolysis
!> @param 		kdop20			reaction rate @ 20 deg-C for DOP hydrolysis
!> @param 		thetadop		reaction rate coefficient for DOP hydrolysis
!> @param 		kdop			temperature dependent reaction rate for DOP hydrolysis
!> @param			rpd				phosphorus-to-detritus stoichiometric mass ratio	

      
!      Subroutine DON(DChemSum,ichem,mex,j,k)
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine dDOP(DChemSum,ichem,j)

	use params
      
      implicit none
      integer :: ichem,j
      real :: det,dop
      real :: kdet
      real :: kdop20,thetadop,kdop
      real :: rpd
      real :: dChemSUM
      
!>> DOP routine (eq. 23 of 2012 Master Plan Appendix D-1)

!>> previous time step WQ concentrations
	det = chem(j,9,1)
	dop = chem(j,11,1)
      
!>> temperature-dependent detritus dissolution hydrolysis rate coefficient 
	kdet = kdet20*thetadet**(Tempw(j,2)-20.)
      
!>> temperature-dependent DOP hydrolysis rate coefficient 
	kdop20 = 0.1
      thetadop = 1.047
	kdop = kdop20*thetadop**(Tempw(j,2)-20.)
      
!>> phosphorus-to-detritus stoichometric mass ration
      rpd = 0.01
      
!>> change in dissolved organic phosphorus concentration
      dChemSUM = rpd*kdet*det-kdop*dop
 
      return
      end




!      Subroutine DOP(DChemSum,ichem,mex,j,k)
! kday now global parameter - no longer needed to be passed into subroutine   
!      Subroutine DOP(DChemSum,ichem,mex,j,k,kday)				!Ammonium as N	
!	!JAM Oct 2010 Chem #11
!
!      use params
!
!
!	FNoC=1./5.681
!	FPoC=1/41.
!	fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
!      DChemSUM=DChemSUM+decay(ichem,11)*Chem(j,9,1)*FPoC		!Dissolution = anaerobic decomposition of organic material   
!      DChemSUM=DChemSUM+decay(12,11)*Chem(j,11,1)				!Dissolution = anaerobic decomposition of organic material   
!
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

c***********************End Subroutine for chemical DeadA*StchP --> DOP*************************