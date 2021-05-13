      Subroutine CelldSal(j,kday,QSalSum,iab,jnb)


      use params
      implicit none
      
      integer, intent(in) :: j, kday              ! local variables declared in hydrod.f
      real, intent(inout) :: QSalSum              ! summation of salinity masses from all tributaries, diversions, and links
      integer, intent(out) :: iab, jnb
      integer :: k, ktrib, kdiv
      real :: dry_depth, ddy1, ddy2, ddym1, ddym2
      real ::  Saltrib, vol1, vol2
!      real :: dddy, dddym                         ! filtered water and marsh depths to use in denominators to avoid div by 0 error

      !>> Set minimum depth value (avoids div-by-zero errors and sets dry cells to turn off salinity change calculations)
      dry_depth = 0.1             ! depth value underwhich the open water (or marsh area) will be considered dry

      ddy1  = Es(j,1) - Bed(j)    ! previous timestep water area depth
      ddy2  = Es(j,2) - Bed(j)    ! current timestep water area depth
      ddym1 = Eh(j,1) - BedM(j)   ! previous timestep marsh area depth
      ddym2 = Eh(j,2) - BedM(j)   ! current timestep marsh area depth
      
!      if(ddy1 <= 0.01) then       ! set non-zero water depth for any use as a denominator (avoid div by 0)
!          dddy = 0.01
!      else
!          dddy = ddy1
!      endif

!      if(ddym1 <= 0.01) then
!          dddym = 0.01            ! set non-zero marsh depth for any use as a denominator (avoid div by 0)
!      else
!          dddym = ddym1
!      endif
      
      QSalsum = 0 
!>> update salinity mass flux (Qsalsum) for tributary flows into compartment      
      do ktrib=1,Ntrib
!>> set salinity in tributary to default freshwater salinity value
          Saltrib = 0.0
!>> if tributary flow is negative, use compartment salinity concentration instead of default tributary salinity concentration
          if (Qtrib(ktrib,kday) < 0.0) then
              Saltrib = S(j,1)
          endif             
          QSalsum = QSalsum - Qtrib(ktrib,kday)*Saltrib*Qmult(j,ktrib)
      enddo

!>> update salinity mass flux (Qsalsum) for diversion flows into compartment (diversions no longer modeled separately, but instead are treated as tributaries)
      do kdiv=1,Ndiv
          QSalsum = QSalsum - Qdiv(kdiv,kday)*0.15*Qmultdiv(j,kdiv)
      enddo

!>> update salinity mass flux (Qsalsum) for link flows into/out of compartment            
      do k=1,nlink2cell(j)
          if(icc(j,k) /= 0) then
              if(icc(j,k) < 0) then
                  jnb=jus(abs(icc(j,k)))
              elseif(icc(j,k) > 0) then
                  jnb=jds(abs(icc(j,k)))
              endif  
          endif
          
          iab=abs(icc(j,k))
    
          call salinity(iab,jnb,j,k,Qsalsum)
      enddo
 
      
      !if (Qmarsh(j,2) > 0.0) then                ! combining marsh and OW volumne
      !QSalSum = QSalSum + Qmarsh(j,2)*S(j,1)      ! add marsh-openwater flow exchange to salinity flux
      !endif

      
!>> original salinity change equation 
!      QRain = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
!     &        -fpet*PET(kday,Jet(j)))*As(j,1)*cden      
!      
!      if (dddy <= 0.1) then  
!          QRain = max(QRain,0.0)
!      endif
!
!      QRainm = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
!!     &        -fpet*PET(kday,Jet(j)))*Ahf(j)*cden      
!     &        -fpet*PET(kday,Jet(j)))*Ahydro(j)*cden
!      
!      if (dddym <= 0.1) then  
!          QRainm = max(QRainm,0.0)
!      endif
!     
!      S(j,2) = (S(j,1)*As(j,1)*dddy-QSalsum*dt)
!     &      /(As(j,1)*dddy+(Qsum_in(j)-Qsum_out(j)+QRain)*dt)

!>> updated salinity mass balance equation - with treatment for dry cells - Jan 21
!      if (Es(j,2) - Bed(j) <= 0.01) then 
!          S(j,2) = S(j,1) ! 0.10      ! if dry, don't update salinity, previously this set salinity to min value
!      else
!          S(j,2) = ( S(j,1)*As(j,1)*dddy - QSalsum*dt ) / ( As(j,1)*( Es(j,2)-Bed(j) ) )
!      endif

      !>> new sal concentration = [ (old salinity concentration * old volume) + salinity mass flux at timestep ] / current volume 
      !>>       [kg/m3]         = [ (    [kg/m3]               *     [m3]  ) +      [kg/sec]*[sec]            ] /    [m3]
          
      !>> QSalsum is salinity mass flux at timestep into/out of open water from all flow mechanisms that change the water surface elevation (Es):
      !>>       - tributary flows into compartment
      !>>       - diversion flow into compartment (not used anymore, diversions treated same as tribs now)
      !>>       - all link flows into/out of compartment - tabulated via salinity.f subroutine
      !>>       - marsh-openwater exchange flow
      !>>       - changes in water level due to rain/ET result in changes to total volume (denominator), but does not impact salinity mass within compartment, therefore concentration will be updated correctly      
      
!>> YW testing equation for MP 2023. combining marsh and OW volumne

!      if (dddy > 0.1) then 
!          S(j,2) = (S(j,1)*(As(j,1)*dddy+Ahf(j)*dddym)-QSalsum*dt)
!     &    /(As(j,1)*dddy+Ahf(j)*dddym
!     &    +(Qsum_in(j)-Qsum_out(j)+QRain+QRainm)*dt)
!      else
!          S(j,2) = S(j,1)
!      endif
      
!>> YW TEST
!      if ((dddy > 0.1) .AND. (dddym > 0.1)) then
!         S(j,2) = (S(j,1)*(As(j,1)*dddy+Ahf(j)*dddym)-QSalsum*dt)
!     &    /(As(j,1)*dddy+Ahf(j)*dddym
!     &    +(Qsum_in(j)-Qsum_out(j)+QRain+QRainm)*dt)
!      else if ((dddy > 0.1) .AND. (dddym <= 0.1)) then
!         S(j,2) = (S(j,1)*(As(j,1)*dddy)-QSalsum*dt)
!     &    /(As(j,1)*dddy
!     &    +(Qsum_in(j)-Qsum_out(j)+QRain)*dt)
!      else
!          S(j,2) = S(j,1)
!      endif

!>> ZW TEST 07/22/2020 - better than YW TEST
      
      
      
      if (Ahf(j) > 0) then                                                    ! if there is marsh area
          if (ddym1 > dry_depth) then                                         ! and marsh was not dry in previous timestep
              if (ddym2 > dry_depth) then                                     ! and marsh is not dry in current timestep
                  vol1=(Es(j,1)-Bed(j))*As(j,1)+(Eh(j,1)-BedM(j))*Ahf(j)      ! then use marsh area AND water area for water volume in previous timestep
                  vol2=(Es(j,2)-Bed(j))*As(j,1)+(Eh(j,2)-BedM(j))*Ahf(j)      ! then use marsh area AND water area for water volume in current timestep
              endif
          endif
      else
          vol1=(Es(j,1)-Bed(j))*As(j,1)                                       ! otherwise, only use water area for water volume in previous timestep
          vol2=(Es(j,2)-Bed(j))*As(j,1)                                       ! otherwise, only use water area for water volume in current timestep
      endif
      
      if(ddy2 > dry_depth) then                                 ! if oepn water area is not dry
          S(j,2)= (S(j,1)*vol1-QSalsum*dt)/vol2                               ! calculate salinity mass
      else                                                                    ! otherwise, water area is dry
          S(j,2) = 0.0                                                        ! set salinity of dry cell to zero
!          S(j,2) = S(j,1)                                                     ! set salinity of dry cell to be the same as previous timestep
      endif
      
!>> equation for MP2017 to avoid salinity spike
!      if (dddy > 0.1) then   
!          DSal =  -QSalsum/(As(j,1)*dddy)*dt
!     &       -Dz*S(j,1)/dddy
!      else
!          DSal = 0.0
!      endif
!      
!     S(j,2)=S(j,1)+DSal
         
!      if (S(j,2) > 100) then
!          write(*,*)'comp = ',j
!          write(*,*) 'As =',As(j,2)
!          write(*,*)'sal(t-1) = ',s(j,1)
!          write(*,*)'sal(t) = ',s(j,2)
!          write(*,*)'qsalsum =',qsalsum
!          write(*,*)'depth = ',Es(j,2)-Bed(j)
!          write(*,*)'Es(t-1) =', Es(j,1)
!          write(*,*)'Es(t) =', Es(j,2)
!          write(*,*) 'Dz =',Dz
!      endif

!>> High-pass and low-pass filters on salinity calculation
!      if(S(j,2) < 0.10) then
!          S(j,2) = 0.10
      if(S(j,2) < 0.0) then           ! try a lower floor on allowed salinity value
          S(j,2) = 0.0
      elseif (S(j,2) > 36.) then
          S(j,2)=36.
      endif

!>> Set marsh salinity equal to open water salinity
      Sh(j,2)=S(j,2)  


      return 
	end
c***********************End Subroutine for Change in Cell Salinity*******JAM Oct 2010***********