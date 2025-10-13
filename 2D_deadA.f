!> @file
!> @brief This subroutine calculates change in dead algae (detritus) concentrations.
!> @details This subroutine calculates change in dead algae (detritus) concentrations as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     Es(N,2)         current water surface elevation in compartment (m)
!> @param[in]     Bed(N)          bed elevation of compartment (m)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]		kphy20			reaction rate @ 20 deg-C for phytoplankton mortality
!> @param[in]		thetaphy		reaction rate coefficient for phytoplankton mortality
!> @param[in]		kdet20			reaction rate @ 20 deg-C for detritus dissolution
!> @param[in]		thetadet		reaction rate coefficient for detritus dissolution

!> @param[out]     DChemSum       change in concentration of water quality constituent during timestep
     
!> @param 		dd				water depth
!> @param 		alg				total algae concentration at previous timestep
!> @param 		det				detritus/dead algae concentration at previous timestep      
!> @param 		kphy			temperature dependent reaction rate for phytoplankton mortality
!> @param 		kdet			temperature dependent reaction rate for detritus dissolution
!> @param         rca             carbon-to-chlorophyll A stoichiometric mass ratio
!> @param 		rda  			DET-to-chlorophyll A stoichiometric mass ratio
!> @param			vd				settling rate of detritus (m/day)
      
	  Subroutine dDeadA(DChemSum,ichem,j)
	
	  use params      

      implicit none
      
      integer :: ichem,j
      real :: dd
      real :: alg,det
      real :: kphy
      real :: kdet
	  real :: rca, rda, vd
      real :: DChemSUM


!>> dead algae (detritus) calculation (eq. 24 of 2012 Master Plan Appendix D-1)

!>> current depth in compartment
!	dd = Es(j,2) - Bed(j)
!	dd = max(Es(j,2) - Bed(j),0.01)	 !zw 4/28/2015 be consistent with NO3 and NH4
	  dd = max(Es(j,2) - Bed(j),dry_threshold)	 !zw 4/28/2015 be consistent with NO3 and NH4

!>> previous time step WQ concentrations
	  alg = chem(j,8,1)
	  det = chem(j,9,1)

!>> temperature-dependent phytoplankton mortality rate coefficient
	  kphy = kphy20*thetaphy**(Tempw(j,2)-20.)      

!>> temperature-dependent detritus dissolution hydrolysis rate coefficient 

	  kdet = kdet20*thetadet**(Tempw(j,2)-20.)

!>> detritus settling rate in compartments (m/day) - weighted by marsh area - if compartment is all marsh, settling is at a maximum rate of 0.05 m/day, if there is no marsh (e.g. all open water), there is no settling of detritus
	  vd = 0.05*apctmarsh(j)

!>> carbon-to-chlorophyll ratio
	  rca = 75.0
      
!>> detritus-to-chlorophyll A stoichometric mass ratio
	  rda = 2.44*rca

!>> change in dead algae (detritus) concentration
	  DChemSUM = rda*kphy*alg - (kdet+vd/dd)*det

	  return
	  end


!      Subroutine DeadA(DChemSum,ichem,mex,j,k)
      
! kday now global parameter - no longer needed to be passed into subroutine   
!	Subroutine DeadA(DChemSum,ichem,mex,j,k,kday)  
!	!JAM Oct 2010 Chem #9
!
!	use params      
!
!     Elim=-decay(8,7)*Chem(j,ichem-1,1)*Tcorn/2. 
!	DChemSUM=DChemSUM+decay(ichem,7)*Chem(j,ichem,1)			! Loss to CO2          
!     &	+decay(ichem,2)*Chem(j,ichem,1)+Elim					! Loss to CH4   
!     &	-decay(8,9)*Chem(j,8,1)*Tcorn/1.5						! *Chem(j,9,1)/		!50-->1.5 JAM March 2011 net loss ~ death rate & elimination rate  = - DA
!     &	+decay(ichem,ichem)*Chem(j,ichem-1,1)*Chem(j,ichem,1)	! /     !death of algae = +DA ccc JAM April 18, 2011     &    (Chem(j,ichem,1)+Chem(j,ichem-1,1)+0.0000000001)*
!     &	*Tcorn/50.
!     &	+decay(ichem,10)*Chem(j,ichem,1)	
!     &	+decay(ichem,7)*Chem(j,ichem,1)
!       
!cc	! Loss to CaCO3 mineralization 
!
!	return
!	end
!c***********************End Subroutine for Dead Algae ******************************************