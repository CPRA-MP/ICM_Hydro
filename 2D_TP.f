!> @file
!> @brief This subroutine calculates change in total inorganic phosphorus concentrations.
!> @details This subroutine calculates change in total phosphorus concentrations as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.
!> This subroutine relies on settling velocites of sediments, which is calculated in subroutine vanRijnSediments.
!> Therefore, this dTP subroutine should be called after vanRijnSediments is called.
!> Currently hydrod calls celldSS (which calls vanRijnSediment)  before hydrod calls celldChem (which calls dTP).
!> If celldChem is called before celldSS, this subroutine will use the compartment's settling velocities from the previous timestep.

!> @param[in]     Ahf(N)          portion of compartment that is marsh area(0-1)
!> @param[in]     TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]     Es(N,2)         current water surface elevation in compartment (m)
!> @param[in]     Bed(N)          bed elevation of compartment (m)
!> @param[in]     chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]     s(N,2)          array with current salinity concentration in compartment (mg/L)
!> @param[in]     velset(N,class) array with settling velocities for compartments at current timestep (m/sec)
!> @param[in]     fpp             fraction of TIP in particulate form
!> @param[in]     muph            photosynthesis rate
!> @param[in]     kresp20         reaction rate @ 20 deg-C for phytoplankton respiration
!> @param[in]     thetaresp       reaction rate coefficient for phytoplankton respiration
!> @param[in]     kpo420          temperature-dependent orthophosphate release rate @ 20 degC (kpo4) (1/day) 
!> @param[in]     thetapo4        temperature-dependent orthophosphate rate coefficient (thetapo4)
      

!> @param[out]    DChemSum        change in concentration of water quality constituent during timestep
     
!> @param         dd              water depth
!> @param         alg             total algae concentration at previous timestep
!> @param         dop             dissolved organic P concentration at previous timestep      
!> @param         tip             total inorganic P concentration at previous timestep
!> @param         kdop20          reaction rate @ 20 deg-C for DOP hydrolysis
!> @param         thetadop        reaction rate coefficient for DOP hydrolysis
!> @param         kdop            temperature dependent reaction rate for DOP hydrolysis
!> @param         kresp           temperature dependent reaction rate for phytoplankton respiration
!> @param         rca             carbon-to-chlorophyll A stoichiometric mass ratio
!> @param         rp              temperature dependent orthophosphate release rate from marsh sediments
!> @param         po4release      orthophosphate released from salt marshes - conditional on salinity concentration      
      
!      Subroutine TP(DChemSum,ichem,mex,j,k)
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine dTP(DChemSum,ichem,j) ! THIS CALCULATED TOTAL INORGANIC PHOSPHORUS
    !JAM Oct 2010 Chem #5

      use params      

      implicit none
      
      integer :: ichem,j
      real :: dd
      real :: alg,tip,dop
      real :: kresp
      real :: kdop20,thetadop,kdop
      real :: kpo4
      real :: rca,rp,po4release
      real :: vs,fmarsh
      real :: DChemSUM

!>> total inorganic phosphorus calculation (eq. 22 of 2012 Master Plan Appendix D-1)

!>> current depth in compartment
!   dd = Es(j,2) - Bed(j)
!   dd = max(Es(j,2) - Bed(j),0.01)  !zw 4/28/2015 be consistent with NO3 and NH4
      dd = max(Es(j,2) - Bed(j),dry_threshold)   !zw 4/28/2015 be consistent with NO3 and NH4
      
!>> portion of compartment that is marsh
      fmarsh = Apctmarsh(j)
      
!>> previous time step WQ concentrations
      alg = chem(j,8,1)     
      tip = chem(j,5,1)
      dop = chem(j,11,1)
      
!>> temperature-dependent phytoplankton respiration rate coefficient
      kresp = kresp20*thetaresp**(Tempw(j,2)-20.)      
      
!>> temperature-dependent DOP hydrolysis rate coefficient 
      kdop20 = 0.1
      thetadop = 1.047
      kdop = kdop20*thetadop**(Tempw(j,2)-20.)

!>> temperature-dependent orthophosphate release rate coefficient 
      kpo4 = kpo420*thetapo4**(Tempw(j,2)-20.)

!>> carbon-to-chlorophyll ratio
      rca = 75.0
      
!>> phosphorus-to-chlorophyll A stoichometric mass ratio
      rp = 0.0244*rca
      
!>> settling rate for phosphorus calculations - average settling velocity of four particle classes converted to m/day
!      vs=3600.*24.*(velset(j,1)+velset(j,2)+velset(j,3)+velset(j,4))/4.
      vs=3600.*24.*velset(j,4)
      
!>> calculate salinity-dependent orthophosphate release from marsh sediments
      if (S(j,2) > 1.0) then
          po4release = fmarsh*kpo4/dd
      else
          po4release = 0.0
      endif
      
!>> change in total phosphorus concentration      
      DChemSUM = kdop*dop + (kresp - muph)*rp*alg
     &             - vs*fpp*tip/dd + po4release
    
     
      return
      end
      
      
      
! Nitrogen to carbon ratio in algae
!ccc using DOP and DIP
!      XTP=0.0                                                      ! Use SRP not TP
!   DChemSUM=DChemSUM+decay(ichem,ichem)*Chem(j,ichem,1)*XTP    ! Loss to mineralization
!   DChemSUM=DChemSUM+decay(ichem,8)*GrowAlgae(j,8,8)*FPoC      ! Addition by uptake by live algae
!     & *chem(j,8,1)*XTP 

!   return
!   end

!c***********************End Subroutine for chemical TP******************************************