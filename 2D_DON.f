!> @file
!> @brief This subroutine calculates change in dissolved organic nitrogen concentration.
!> @details This subroutine calculates change in dissolved organic nitrogen concentration as developed for the 2012 Master Plan models.
!> Refer to Appendix D-1 'Eco-hydrology Model Technical Report' from the 2012 Master Plan report for detailed
!> descriptions and development of these water quality routines.

!> @author  Eric White - The Water Institute of the Gulf


!> @param[in]   TempW(N,2)      current water temperature in compartment (deg Celsius)
!> @param[in]   chem(N,ichem,1) array with previous timestep water quality constituent concentrations (mg/L)
!> @param[in]   muph            phytoplankton photosynthesis rate (1/day)
!> @param[in]   kdet20          reaction rate @ 20 deg-C for detritus dissolution
!> @param[in]   thetadet        reaction rate coefficient for detritus dissolution

!> @param[out]  DChemSum       change in concentration of water quality constituent during timestep

!> @param       det             detritus concentration at previous timestep
!> @param       don             DON concentration at previous timestep
!> @param       kdet            temperature dependent reaction rate for detritus dissolution
!> @param       kdon20          reaction rate @ 20 deg-C for DON hydrolysis
!> @param       thetadon        reaction rate coefficient for DON hydrolysis
!> @param       kdon            temperature dependent reaction rate for DON hydrolysis
!> @param       kresp20         reaction rate @ 20 deg-C for phytoplankton respiration
!> @param       thetaresp       reaction rate coefficient for phytoplankton respiration
!> @param       kresp           temperature dependent reaction rate for phytoplankton respiration
!> @param       rnd             nitrogen-to-chlorphyll A stoichiometric mass ratio


      
!      Subroutine DON(DChemSum,ichem,mex,j,k)
! kday now global parameter - no longer needed to be passed into subroutine   
      Subroutine dDON(DChemSum,ichem,j)

      use params
      
      implicit none
      integer :: ichem,j
      real :: det,don
      real :: kdet
      real :: kdon
      real :: rnd
      real :: dChemSUM
      
!>> DON routine (eq. 20 of 2012 Master Plan Appendix D-1)

!>> previous time step WQ concentrations
      det = chem(j,9,1)
      don = chem(j,10,1)
      
!>> temperature-dependent detritus dissolution hydrolysis rate coefficient 
      kdet = kdet20*thetadet**(Tempw(j,2)-20.)
      
!>> temperature-dependent DON hydrolysis rate coefficient 
      kdon20 = 0.015
      thetadon = 1.047
      kdon = kdon20*thetadon**(Tempw(j,2)-20.)
      
!>> nitrogen-detritus stoichometric mass ratio
      rnd = 0.072
      
!>> change in dissolved organic nitrogen concentration
      dChemSUM = rnd*kdet*det-kdon*don

      return
      end
      
      !JAM Oct 2010 Chem #10

!   FNoC=1./5.681
!   fnodin=Chem(j,1,1)/(Chem(j,3,1)+10E-08)
!      DChemSUM=DChemSUM+decay(ichem,9)*                    ! Hydrolysis
!     & Chem(j,9,1)*FNoC+decay(ichem,2)*Chem(j,ichem,1)

!   return
!   end
    
!c-------------------------------
!c        1 = NO3 + NO2
!c        2 = NH4 
!c        3 = DIN (1+2)
!c        4 = Organic N
!c        5 = TP
!c        6 = TOC
!c        7 = DO
!c     8 = Live Algae
!c     9 = Dead Algae = Detritus
!c        10= DON
!c        11= DOP
!c        12= DIP   !Partition SRP = ParP*TP
!c        13= ChLa  !Partition Chla= ParCla*LivA
!c        14= POP !-EDW used to say 14=TKN !-EDW
!c-------------------------------
    
!c***********************End Subroutine for DeadA --> DON****************************************