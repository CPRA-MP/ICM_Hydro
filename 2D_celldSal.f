!      Subroutine CelldSal(QSalSUM,j,kday,k,SalTRIBj,dref,Tres)
      Subroutine CelldSal(j,kday)
      
      !>> QSalsum is salinity mass flux at timestep into/out of open water from all flow mechanisms that change the water surface elevation (Es):
      !>>       - tributary flows into compartment
      !>>       - diversion flow into compartment (not used anymore, diversions treated same as tribs now)
      !>>       - all link flows into/out of compartment - tabulated via salinity.f subroutine
      !>>       - marsh-openwater exchange flow
      !>>       - changes in water level due to rain/ET result in changes to total volume (denominator), but does not impact salinity mass within compartment, therefore concentration will be updated correctly      

      !*************************************************************************************************************************************
      !>> new sal concentration = [ (old salinity concentration * old volume) + salinity mass flux at timestep ] / current volume 
      !>>       [kg/m3]         = [ (    [kg/m3]               *     [m3]  ) +      [kg/sec]*[sec]            ] /    [m3]
      !*************************************************************************************************************************************    


      use params
      
!      implicit none
      integer :: j, kday, ktrib, kdiv, marsh_link_flow, k, iab, jnb
      real :: QSalsum, Saltrib, Qsalsum_b4link
      real :: dry_depth, dry_salinity
      real :: vol1, vol2, marsh_vol1, marsh_vol2
      real :: ddy1, ddy2, dddy, ddym1, ddym2, dddym
      real :: salmaxcon, ds, Qlink, Csalface
      
      !cden=1./1000./24./3600.		! mm/d to m/s conversion

      !>> Define depth, in meters, for dry cells that will turn off salinity change calculations 
      !      this is used in other celldXXX subroutines but each subroutine may have a separate dry depth value assigned - double check for consistency
!      dry_depth = 0.05
      dry_depth = dry_threshold
!      dry_salinity = 0.1         ! set dry cell salinity to 0.1 ppt
      dry_salinity = S(j,1)      ! set dry cell salinity to value from previous timestep
      
!>> Calculate water and marsh depths for current and previous timesteps
      ddy1 = Es(j,1)-Bed(j)
      ddy2 = Es(j,2)-Bed(j)
      ddym1 = Eh(j,1)-BedM(j)
      ddym2 = Eh(j,2)-BedM(j)
      
!>> Set minimum depth value to avoid div-by-zero errors
!      if(ddy1 <= 0.01) then
!          dddy = 0.01
      if(ddy1 <= dry_threshold) then
          dddy = dry_threshold
      else
          dddy = ddy1
      endif

!      if(ddym1 <= 0.01) then
!          dddym = 0.01
      if(ddym1 <= dry_threshold) then
          dddym = dry_threshold
      else
          dddym = ddym1
      endif
      
      QSalsum = 0
      salmaxcon = 0.0      
          
!>> update salinity mass flux (Qsalsum) for tributary flows into compartment      
      do ktrib=1,Ntrib
!>> set salinity in tributary to default freshwater salinity value (0.1 ppt)
          Saltrib = 0.1
!>> if tributary flow is negative, use compartment salinity concentration instead of default tributary salinity concentration
          if (Qtrib(ktrib,kday) < 0.0) then
              Saltrib = S(j,1)
          endif             
          QSalsum=QSalsum-Qtrib(ktrib,kday)*Saltrib*Qmult(j,ktrib)
      enddo
      
      if (Qsalsum .ne. 0) then
          salmaxcon = Saltrib     ! set first value of maximum salinity concentration to salinity of tributaries if there was any tributary flow
      endif
      
!>> update salinity mass flux (Qsalsum) for diversion flows into compartment (diversions no longer modeled separately, but instead are treated as tributaries)
      do kdiv=1,Ndiv
          QSalsum = QSalsum-Qdiv(kdiv,kday)*0.15*Qmultdiv(j,kdiv)     ! diversion salinity assumed here to be 0.15 - but this is not used anymore since Ndiv is always set to 0
      enddo

!>> flag that will be set if overland marsh links have flow
      marsh_link_flow = 0          
      Qsalsum_b4link = 0
      
!>> update salinity mass flux (Qsalsum) for link flows into/out of compartment            
      do k=1,nlink2cell(j)
          iab=abs(icc(j,k))
          if(icc(j,k) /= 0) then
              if(icc(j,k) < 0) then
                  jnb=jus(iab)
              elseif(icc(j,k) > 0) then
                  jnb=jds(iab)
              endif  
          endif

          Qsalsum_b4link = Qsalsum

!          call salinity(mm,iab,jnb,j,k,Qsalsum)
          if(iab > 0) call salinity(iab,jnb,j,k,Qsalsum)
          
!>> check if current link has flow during timestep
          if (Q(iab,2) .ne. 0.0) then
              salmaxcon = max( salmaxcon,SL(iab,2) )          ! SL() is the updated face salinity concentration for link iab calculated in salinity()
          endif
          
!>> check if marsh overland links have salinity convection and/or dispersion through link flow
          if (linkt(iab) == 8) then
              if ( Qsalsum_b4link .ne. Qsalsum ) then
                  marsh_link_flow = 1
              endif
          endif
      
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


!>> Check if the marsh area of the compartment should be used for calculating total water volume in salinity concentration calculations
!>> if marsh depth is below dry depth threshold, there still may be shallow flow into the marsh area bringing salinity mass
!>> check if there was flow into the marsh that could potentially add salinity mass
!>> if so, then include volume of shallow water over marsh surface in total volume
      marsh_vol1 = 0.0
      marsh_vol2 = 0.0
      if( Ahf(j) > 0 ) then                           ! check if there is marsh area
!          if ( ddym1 > dry_depth ) then               ! check if marsh was dry in previous timestep
!              if ( ddym2 > dry_depth ) then           ! check if marsh is dry in current timestep
!                  marsh_vol1 = ddym1*Ahf(j)
!                  marsh_vol2 = ddym2*Ahf(j)
!              else
!                  marsh_vol1 = 0.0
!                  marsh_vol2 = 0.0
!                  if (marsh_link_flow == 1) then      ! marsh is dry but there was overland marsh flow
!                      marsh_vol1 = ddym1*Ahf(j)
!                      marsh_vol2 = ddym2*Ahf(j)
!                  endif
!              endif
!          else
!              marsh_vol1 = 0.0
!              marsh_vol2 = 0.0
!          endif
!      else
!          marsh_vol1 = 0.0
!          marsh_vol2 = 0.0
          if ( ddym1 > dry_depth ) then               ! check if marsh was dry in previous timestep
              marsh_vol1 = max(ddym1*Ahf(j),0.0)
              marsh_vol2 = max(ddym2*Ahf(j),0.0)
          
          elseif ( ddym2 > dry_depth ) then               ! check if marsh was dry in current timestep
              marsh_vol1 = max(ddym1*Ahf(j),0.0)
              marsh_vol2 = max(ddym2*Ahf(j),0.0)
          endif
      endif

!     openwater volume
      vol1 = 0.0
      vol2 = 0.0
      if ( ddy1 > dry_depth ) then               ! check if openwater was dry in previous timestep
          vol1 = max(ddy1*As(j,1),0.0)
          vol2 = max(ddy2*As(j,1),0.0)
      
      elseif ( ddy2 > dry_depth ) then               ! check if openwater was dry in current timestep
          vol1 = max(ddy1*As(j,1),0.0)
          vol2 = max(ddy2*As(j,1),0.0)
      endif

!      vol1 = ddy1*As(j,1) + marsh_vol1
!      vol2 = ddy2*As(j,1) + marsh_vol2
      vol1 = vol1 + marsh_vol1
      vol2 = vol2 + marsh_vol2


      if(ddy2 > dry_depth) then
!1/15/2024      if(vol2 > 0) then
!          S(j,2)= ( S(j,1)*vol1 - QSalsum*dt ) / max(0.01,vol2)   
          S(j,2)= ( S(j,1)*vol1 - QSalsum*dt ) / vol2   
          ds = S(j,2) - S(j,1)
          
          !>> vol2 includes changes to water volume from precip and ET (since it is calculated from depth, ddy2) 
          !>> this net precip volumetric change does not need to be included in Qsalsum since precip and ET are 
          !>> always assumed to have zero salinity mass in those water volumes
          
!          !>> set default value that will cap the allowable change in salinity for the current timestep
!          maxDSal = 0.1           ! this dS cap will cap salinity changes to 0.1 ppt but is only applied to forested compartments below
!          
!          ! maxDSal = 0.1*SALAV(j)  ! this dS cap will cap a change in salinity to 10% of the background mean salinity of the compartment
!                                  ! (mean salinity is updated every day using a rolling window that will equal the mean one-day salinity at the end of Jan 1, 
!                                  ! will equal the 6 month salinity on June 30, and so on      
!
!          !>> check to see if the compartment is 10% open water or less (indicates possibly location in swamp/forested region
!          !>> if it is small water-to-land ratio, then check to see if it is fresh
!          !>> if both criteria are met - do not allow the salinity to increase a large amount for the current timestep
!          !>> this is a filter to prevent spiking of salinity due to swamp forest compartment geographic representation issues
!          if ( Apctwater(j) < 0.1) then
!              if (SALAV(j) < 2.0) then
!                  if ( abs(DSal) > maxDSal ) then
!                      S(J,2) = S(j,1) + maxDSal*DSal/abs(Dsal)      ! ds/abs(ds) gets the directionality of the ds vector and applies the max dS filter to the current timestep
!                      write(*,'(A,x,I)') 'max dsal exceeded. comp:',j
!                  endif
!              endif
!          endif
      else
          S(j,2) = dry_salinity
!          write(*,'(A,I,A,F,A,F,A,I)') 'dry sal in comp: ',j,' vol:',vol2,' m_vol: ',marsh_vol2,' marsh_flow: ',marsh_link_flow
      endif
      
!>> equation for MP2017 to avoid salinity spike
!      if (dddy > 0.1) then   
!          DSal =  -QSalsum/(As(j,1)*dddy)*dt
!     &	    -Dz*S(j,1)/dddy
!      else
!          DSal = 0.0
!      endif
!      
!     S(j,2)=S(j,1)+DSal
         
! for code debugging
      if ((S(j,2) < 0) .or. (S(j,2) > salmax)) then
          write(1,*)'comp = ',j
          write(1,*)'As =',As(j,1)
          write(1,*)'sal(t-1) = ',S(j,1)
          write(1,*)'sal(t) = ',S(j,2)
          write(1,*)'depth(t-1) = ',Es(j,1)-Bed(j)
          write(1,*)'depth(t) =', Es(j,2)-Bed(j)
          write(1,*)'Dz =',Es(j,2)-Es(j,1)
          write(1,*)'vol(t-1) =', vol1,marsh_vol1
          write(1,*)'vol(t) =', vol2,marsh_vol2
          write(1,*)'qsalsum =',qsalsum
          do k=1,nlink2cell(j)
              iab=abs(icc(j,k))
              if(icc(j,k) /= 0) then
                  if(icc(j,k) < 0) then
                      jnb=jus(iab)
                  elseif(icc(j,k) > 0) then
                      jnb=jds(iab)
                  endif  
              endif
              Qlink = sicc(j,k)*Q(iab,2)
              if(abs(Qlink)>0) then
                  write(1,*)'LinkID=',iab,'Type=',linkt(iab),'Q=',Qlink
                  if(Q(iab,2) >= 0.0) then
                      Csalface= ((fa(iab)*S(jus(iab),1)
     &                  +fb(iab)*S(jds(iab),1)))
                  else
                      Csalface= ((fa(iab)*S(jds(iab),1)
     &                  +fb(iab)*S(jus(iab),1)))
                  endif
                  write(1,*)'S(jus)=',S(jus(iab),1),
     &                      'S(jds)=',S(jds(iab),1), 'Csalface=',Csalface
                  write(1,*)'QSALadvec=',sicc(j,k)*(Q(iab,2))*Csalface
                  write(1,*)'QSALdiffu=',fe*EAOL(iab)*(S(j,1)-S(jnb,1))
              endif
          enddo
          stop
      endif

!>> High-pass and low-pass filters on salinity calculation
!      if(S(j,2) < 0.10) then
!          S(j,2) = 0.10
      if(S(j,2) < 0.0) then
          S(j,2) = 0.0
!      elseif (salmaxcon > 0.0) then           ! if salmaxcon is zero then there were no tributary or link flows into compartment for timestep - so salinity does not need to be capped by surrounding concentrations
!          if ( S(j,2) > S(j,1) ) then         ! only filter by max salinity concentration from flows to compartments that have increased in salinity - otherwise a higher saline waterbody would be reduced to match fresh inflows
!              if (S(j,2) > salmaxcon ) then   ! salinity concentration in compartment cannot be greater than the maximum salinity concentration of all connecting links for the timestep
!                  S(j,2) = salmaxcon
!              endif
!          endif
!      elseif (S(j,2) > 36.) then
!          S(j,2)=36.
      elseif (S(j,2) > salmax) then
          S(j,2) = salmax
      endif

!>> Set marsh salinity equal to open water salinity
      Sh(j,2)=S(j,2)  

!!>> Calculate daily average values reported out to output files ! moved to hydro
!     if(daystep == 1) then
!     !>> reset average salinity value at start of day
!         SALAV(j) = S(j,2)*dt/(3600.*24.)
!      else
!          !>> Update average salinity by timestep's contribution to daily average 
!         SALAV(j)=SALAV(j) + S(j,2)*dt/(3600.*24.)
!     endif

      return 
      end
!***********************End Subroutine for Change in Cell Salinity*******JAM Oct 2010***********
