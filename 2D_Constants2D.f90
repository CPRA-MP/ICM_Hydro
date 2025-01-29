!> @file
!> @brief This is the subroutine to initialize constant variables.
!> @details This is the subroutine to initialize constant variables

      subroutine Constants2D

      use params

      implicit none
      integer :: j            !iterators
      real :: Pardpmr
      
!>> Set initial conditions for constants and model parameters that are not included in input text files
!   Initial Conditions
      g=9.81                  ! Gravity (m/s2)
      pi=4.0*atan(1.0)
!      TemI = 15.                ! Initial Water Temperature
      KnN= 20.                ! (ug/L) DIN Michaelis Constant  Thomann & Mueller
      KnP= 3.                 ! (ug/L) P Michaelis Constant
      KnSS= 50.               ! (mg/L) SS Michaelis Constant  chged 30 to 50 JAM March 2011
      KnSal=4.                ! (ppt) Salinity  Michaelis Constant
    
      ParP= 0.4               ! SRP/TP in Tribs   J. Day 1994 BCS trial
      ParDOP= 0.1             ! DOP/TP in Tribs   J. Day 1994 BCS trial
      PARPOP=1.-ParDOP-ParP   ! POP/TP in Tribs   J. Day 1994 BCS trial
      ParPMR= 0.2             ! SRP/TP in Diversions J. Day 1994 BCS trial
      Pardpmr= 0.05
      PPMR=1.-ParPMR-Pardpmr
      ParSand=0.05            ! sand/TSS in Tribs and MR typical
      ParCLa=0.03             ! Partition LivA --> ChlA
    
      consd=24*3600           ! sec to days
      conv=0.001*consd        ! (mg/L)*(m3/s) --> kg/d
      stds=0.
!      nuo=0.000001           ! DEFAULT Viscosity
      fa_def=1                !fa is an link attribute now and a calibration factor, so fa_def should be always 1 
      cden=1./1000./24./3600.       ! mm/d to m/s conversion

!>> Generate array for Day of Year value for the first day of each month
      month_DOY(1) = 1
      month_DOY(2) = 32
      month_DOY(3) = 60
      month_DOY(4) = 91
      month_DOY(5) = 121
      month_DOY(6) = 152
      month_DOY(7) = 182
      month_DOY(8) = 213
      month_DOY(9) = 244
      month_DOY(10) = 274
      month_DOY(11) = 305
      month_DOY(12) = 335
!>> Update first day of month for March through December during a leap year
      if(simdays == 366) then
          do j = 3,12
              month_DOY(j) = month_DOY(j)+1
          enddo
      endif


!>> determine upper CSS threshold for each sediment class where resuspension of bed material will be deactivated
!>> These values are based on an assumption that at a TSS value equal to the threshold concentration, the particle size distribution is 10% sand, 45% silt, and 45% clay.
      CSSresusOff(1) = TSSresusOff*0.1
      CSSresusOff(2) = TSSresusOff*0.45
      CSSresusOff(3) = TSSresusOff*0.45
      CSSresusOff(4) = TSSresusOff*0.45

!>> Initialize hardcoded salinity and stage control trigger flags
      SalLockStatusHNC = 1  !Open = 1 Close = -1
      SalLockTriggerHNC = 1
      StgTriggerSuperiorCanal = 1
      StgTriggerStatusSuperiorCanal = 1
      !SalLockStatusCSC = 1
      !SalLockTriggerCSC = 1
      !Atch_div_onoff = 1

      return
      end
