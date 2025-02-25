!     Subroutine CelldQ(QSUM,Dz,j,fcrop,mm)
	
! day, dday, and kday now global parameters - no longer needed to be passed into subroutine      
      Subroutine CelldQ(j,kday,fcrop,mm,dday)			!Solves Continuity for cell j 

      
!> @param[out]        Qsum_out(j)     sum of all flows leaving compartment
!> @param[out]        Qsum_in(j)      sum of all flows entering compartment
!> @param[out]        Qsum_abs(j)     magnitude of all flows entering AND leaving compartment
      
      use params      

      implicit none
      real :: Qlink,Dzhlim,Elevel,flo_trib,flo_div,mindz
      integer :: j,kday,mm,jn,jjn,k,iab,jnb
      real :: Qsum,Qsumh,dday,fcrop,fpc,Qhhf,Ahmf,Qupld,ddy1,ddym1,Qavail
      real :: Qow,sndz,Dzh,sndzh,PETuse,ETmin,Het,fET

      Qsum=0.0						!JAM Oct 2010
      Qsumh=0.0
      cden=1./1000./24./3600.		! mm/d to m/s conversion

!>> initialize directional components of cumulative flowrates that are used for sand transport equations in vanRijn
      Qsum_out(j) = 0.0
      Qsum_in(j) = 0.0
      Qsum_abs(j) = 0.0
!      Qsum_out_link(j) = 0.0

      ddy1=max(Es(j,1)-Bed(j),0.0)
      ddym1=max(Eh(j,1)-BedM(j),0.0)

!>> Update cumulative open water flow rates in compartment from input boundary tributary flows      
!>> sign convention on open water flow = positive is flow out of compartment
      do jn = 1,ntrib
!          if (isnan(Qmult(j,jn))) then  !move to infile.f - ZW 12/18/2023
!              Qmult(j,jn) = 0.0
!          endif
          flo_trib = Qtrib(jn,kday)*Qmult(j,jn)
          if((ddy1==0) .and. (flo_trib<0)) flo_trib = 0   !no outflow when compartment is dry in previous time step
          Qsum = Qsum - flo_trib
          Qsum_abs(j) = Qsum_abs(j) + abs(flo_trib)
          if (flo_trib >= 0.0) then
              Qsum_in(j) = Qsum_in(j) + flo_trib
          else
              Qsum_out(j) = Qsum_out(j) + abs(flo_trib)
          endif
      enddo

!>> Update cumulative open water flow rates in compartment from input boundary diversion flows
!>> sign convention on open water flow = positive is flow out of compartment
      do jjn=1,Ndiv
          flo_div = Qdiv(jjn,kday)*Qmultdiv(j,jjn)
          if((ddy1==0) .and. (flo_div<0)) flo_div = 0    !no outflow when compartment is dry in previous time step
          Qsum = Qsum - flo_div
          Qsum_abs(j) = Qsum_abs(j) + abs(flo_div)
          if (flo_div >= 0.0) then
              Qsum_in(j) = Qsum_in(j) + flo_div
          else
              Qsum_out(j) = Qsum_out(j) + abs(flo_div)
          endif
      enddo

!=====Rainfall runoff (precip-evap)
      ! cden: mm/day to m/s conversion factor
      ! Rain: mm/day rainfall
      ! PET:  mm/day potential ET
      ! fpet: 1=use PET input data; 0: use average ET value
      ! ETA:  average ET value to use if fpet = 0
!c correction JAM March 26 2007  ---correction JAM Aug 10 090.50093
      fpc= percent(j)*(1.-fcrop)/100.+fcrop               ! multiplier on PET for marsh areas ; if percent(j)=10 and fcrop=0.5, then the marsh area will evaporate 0.55*PET for the day
!      soilm = max(0.0001, esho(j)-BedM(j))				!JAM Oct 2010
!      shh   = max(0.0001, Eh(j,1)-BedM(j))				!JAM Oct 2010
!      rhh   = max(0.0001, shh/soilm)			            !JAM Oct 2010

      PETuse=(1-fpet)*ETA(Jet(j))-fpet*PET(kday,Jet(j))
!>> Excess rainfall runoff on marsh
!      Qhhf=Ahf(j)*(Rain(kday,jrain(j))-PET(kday,Jet(j))*fpc)	!include fpc for marsh area PET - ZW 12/18/2023
!zw 1/18/2024 add ET correction in marsh area following MP2012 AA code
      ETmin=0.20                 !minimum ET reduction factor
      Het=0.25                   !depth below which ET is reduced
      fET=max(ETmin,min(1.0,ddym1/Het)) !reduction factor for reduced sunlight through marsh
      if (ddym1<=dry_threshold) then
          Qhhf=Ahf(j)*Rain(kday,jrain(j))*cden 
      else
          Qhhf=Ahf(j)*(Rain(kday,jrain(j))-PETuse*fET)*cden ! in m^3/s !*Max(1.,rhh*rhh))	!JAM Oct 2010
          Qavail=ddym1*Ahf(j)/dt
          Qhhf=max(Qhhf,-Qavail)                            !prevent excessive evap over marsh
      endif
!>> !YW! ignore PET at low marsh water level to avoid Eh too low
!      if((Eh(j,1)-BedM(j))>0.1) then
!          Qhhf=Ahf(j)*(Rain(kday,jrain(j))
!     &    -PET(kday,Jet(j))*Max(1.,rhh*rhh))	!JAM Oct 2010
!      else
!          Qhhf=Ahf(j)*Rain(kday,jrain(j))
!      endif

      Ahmf=Ahydro(j)-Ahf(j)
!>> Update cumulative flow rate in marsh based on excess rainfall runoff on upland area
!>> sign convention on marsh flow = positive flow is from marsh to open water
!      Qupld=Qhhf+(max(0.0,(Rain(kday,jrain(j))-PETuse*fpc))*Ahmf)*cden	 
      Qupld=Qhhf+(max(0.0,0.05*Rain(kday,jrain(j)))*Ahmf)*cden	 !zw 02/23/2025 Rational Method = CiA (assuming C=0.05 in all upland types)

!>> Update cumulative flow rate in open water based on excess rainfall runoff on open water area
!>> sign convention on open water flow = positive is flow out of compartment
      if (ddy1<=dry_threshold) then
          Qow=Rain(kday,jrain(j))*As(j,1)*cden
      else
          Qow=(Rain(kday,jrain(j))-PETuse)*As(j,1)*cden
          Qavail=ddy1*As(j,1)/dt
          Qsum=Qsum-max(Qow,-Qavail)                      !prevent excessive evap over openwater
      endif
      QRain(j)=Qupld+Qow  !ZW 1/31/2024: QRain is the total rainfall runoff within compartment j

!      Qsumh=Qsumh-(Qhhf+max(0.0,(Rain(kday,jrain(j))
!     &	 -PET(kday,Jet(j))*fpc))*Ahmf)*cden								!Runoff>0    !JAM Oct 2010          
!>> modified zw 04/29/2020 to deal with upland compartments w/o marsh area
      if(Ahf(j) > 0.0) then	
          Qsumh=Qsumh-Qupld
      else
          Qsum=Qsum-Qupld
      endif
!=====end rainfall runoff

!>> Update cumulative marsh flowrate for marsh-to-open water exchange flow (calculated via Kadlec Knight in hydrod)      
!>> sign convention on marsh flow Qmarsh = positive flow is from open water-to-marsh flow
      Qsumh=Qsumh-Qmarsh(j,2)											!Qsumh is net of marsh flows 
!>> add marsh exchange flow to open water cumulative flow
      Qsum=Qsum+Qmarsh(j,2)

!>> collect flow from connecting links
!>> sign convention on link flows: positive is flow out of compartment
      do k=1,nlink2cell(j)
          iab=abs(icc(j,k))
          Qlink = 0
          if(iab > 0) Qlink = sicc(j,k)*Q(iab,2)
          if(abs(Qlink)>0) then
!>> if link type is marsh overland flow type, add flow to marsh flow sum
!>> if compartment has no marsh area, add marsh overland flow in compartment to water area
!>> if  marsh link flow is negative, flow is entering marsh from neighboring marsh
              if (linkt(iab) == 8) then
                  if (Ahf(j) > 0) then
                      Qsumh = Qsumh + Qlink
                  else
                      Qsum = Qsum + Qlink
                  endif
! distribute ridge link flows to both marsh & OW proportionly based on the percentage area
              elseif (linkt(iab) == 9) then
                  Qsumh = Qsumh + Qlink*Ahf(j)/(Ahf(j)+As(j,1))
                  Qsum = Qsum + Qlink*As(j,1)/(Ahf(j)+As(j,1))
!>> if link is not marsh overland flow type, then add flow to water flow sum
              else    
                  Qsum = Qsum + Qlink
              endif

!>> calculate magnitude of all flow and only in/out flows for use in sediment routing equations            
              Qsum_abs(j) = Qsum_abs(j) + abs(Qlink)
            
              if (Qlink >= 0.0) then
                  Qsum_out(j) = Qsum_out(j) + Qlink
!                  Qsum_out_link(j) = Qsum_out_link(j) + Qlink
              else
                  Qsum_in(j) = Qsum_in(j) + abs(Qlink)
              endif
              
          endif
      enddo

!>> Reduce calculated flow rates out of compartment based on available volume FOR SEDIMENT ROUTING
!>> ***NOTE**** this assumes that tributary and diversion flowrates are by and large only flowing INTO the compartment
!>> ***NOTE**** if large flowrates out of the compartment are due to the observed tributary/diversion flowrates in the input files
!>> ***NOTE**** then mass will not be conserved correctly for sediment routing
!>> Determine maximum flowrate out if all volume in compartment drains      
!      Qsum_out_max = ((Es(j,2) - Bed(j))*As(j,1))/dt
!!>> If compartment will drain completely, reduce link flowrates
!      if (Qsum_out(j) > Qsum_out_max) then
!          do k=1,nlink2cell(j)
!              if(abs(icc(j,k)) /= 0) then
!                  Qlink = sicc(j,k)*Q(abs(icc(j,k)),1)
!          !>> determine portion of flow that is fixed due to tributary and diversion flows
!                  Qfixed = Qsum_out(j) - Qsum_out_link(j)
!          !>> determine reduction in flowrated needed to meet available
!                  Qsum_red = Qsum_out_max - Qfixed
!                  
!                  if (Qlink >= 0.0) then
!                  !>> Determine current link's contribution to total flow out
!                      Qwt = Qlink/Qsum_out_link(j)
!                  !>> Reduce link flow rate
!                      Qred(abs(icc(j,k)),1) = Qwt*Qsum_red
!                  endif
!              endif
!          enddo
!      endif


!>> Calculate change in open water stage
!>> negative flows are into compartment and will result in positive deltaZ
      Dz=((-Qsum)/As(j,1))*dt				                !Euler method,3.
      sndz = 1.0
      if (Dz < 0) sndz = -1.0

! for code debugging
      if (abs(Dz) >= oscilflag) then
          write(*,*) 'Comp=',j,'timestep=',mm,'dz=',dz
          write(*,*) 'Es(t-1)=',Es(j,1),'Eh(t-1)=',Eh(j,1)
          write(*,*) 'Qmarsh=',Qmarsh(j,2)
          write(*,*) 'Es(t)=',Es(j,2)
          write(*,*) 'Qsum=',Qsum
          write(1,*) 'Comp=',j,'timestep=',mm,'dz=',dz
          write(1,*) 'Es(t-1)=',Es(j,1),'Eh(t-1)=',Eh(j,1)
          write(1,*) 'Qmarsh=',Qmarsh(j,2)
          write(1,*) 'Es(t)=',Es(j,2)
          write(1,*) 'Qsum=',Qsum
          do k=1,nlink2cell(j)
              iab=abs(icc(j,k))
              if(iab /= 0) then
                  Qlink = sicc(j,k)*Q(iab,2)
                  if(abs(Qlink)>0) then
                     write(1,*)'LinkID=',iab,'Type=',linkt(iab),'Q=',Qlink
                     if(linkt(iab)==8) then
                         write(1,*)'Eh(jus)=',Eh(jus(iab),1),
     &                         'Eh(jds)=',Eh(jds(iab),1)
                     else
					     write(1,*)'Es(jus)=',Es(jus(iab),1),
     &                         'Es(jds)=',Es(jds(iab),1)
                     endif
                  endif
              endif
          enddo
          stop
      endif

      if(abs(Dz) > maxdz) then
          Dz = sndz*oscilflag   ! Dz = sndz*maxdz
          write(*,898) 'Large deltaZ in compartment ',j,'timestep=',mm,
     & 'deltaZ set to max value allowed in RunControlR.dat'
          write(1,898) 'Large deltaZ in compartment ',j,'timestep=',mm,
     & 'deltaZ set to max value allowed in RunControlR.dat'
          !pause
898   Format(x,A,I0,x,A,x,I0,x,A)
      endif

!>> Update water elevation of open water portion - don't let water level drop below bed elevation
      Es(j,2)=max(Es(j,1) + Dz,Bed(j))				                !Stage in storage cells (m)
      
!>> Calculate change in marsh stage
!>> Negative marsh flow is into marsh and will result in positive deltaZmarsh
      if (Ahf(j) > 0.0) then
          Dzh=(-Qsumh)/Ahf(j)*dt
          sndzh = 1.0
          if (Dzh < 0) sndzh = -1.0

! for code debugging
          if (abs(Dzh) >= oscilflag) then
              write(*,*) 'Comp=',j,'timestep=',mm,'dzh=',dzh
              write(*,*) 'Es(t-1)=',Es(j,1),'Eh(t-1)=',Eh(j,1)
              write(*,*) 'Qmarsh=',Qmarsh(j,2)
              write(*,*) 'Eh(t)=',Eh(j,2)
              write(*,*) 'Qsumh=',Qsumh
              write(1,*) 'Comp=',j,'timestep=',mm,'dz=',dz
              write(1,*) 'Es(t-1)=',Es(j,1),'Eh(t-1)=',Eh(j,1)
              write(1,*) 'Qmarsh=',Qmarsh(j,2)
              write(1,*) 'Eh(t)=',Eh(j,2)
              write(1,*) 'Qsumh=',Qsumh
              do k=1,nlink2cell(j)
                  iab=abs(icc(j,k))
                  if((iab>0) .and. (linkt(iab)==8)) then
                      Qlink = sicc(j,k)*Q(iab,2)
                      if(abs(Qlink)>0) then
                         write(1,*)'LinkID=',iab,'Type=',linkt(iab),'Q=',Qlink
                         write(1,*)'Eh(jus)=',Eh(jus(iab),1),
     &                         'Eh(jds)=',Eh(jds(iab),1)
                      endif
                  endif
              enddo
              stop
          endif

          if(abs(Dzh) > maxdz) then
            Dzh = sndzh*oscilflag   ! Dz = sndz*maxdz
            write(*,898) 'Large deltaZh in compartment ',j,'timestep=',mm,
     & 'deltaZh set to max value allowed in RunControlR.dat'
            write(1,898) 'Large deltaZh in compartment ',j,'timestep=',mm,
     & 'deltaZh set to max value allowed in RunControlR.dat'
          endif

!>> Update marsh water elevation
          Eh(j,2)= Max(Eh(j,1)+Dzh, BedM(j))		!JAM Oct 2010 Marsh stage (m) 
! move Qmarshmax to hydrod.f - ZW 12/19/2023
! !>> Determine level where open water stage and marsh stage will be equal in compartment
!          Elevel = (As(j,1)*Es(j,2)+Ahf(j)*Eh(j,2)) / (As(j,1)+Ahf(j))
!!>> Determine change in marsh level to reach equilbrium with open water marsh stage
!          Dzhlim = max(Elevel,BedM(j)) - Eh(j,2)
!!>> Calculate magnitude of flowrate needed into/out of marsh - directionality assigned in hydrod
!          Qmarshmax(j) = abs(Dzhlim)*Ahf(j)/dt
              
      else
!>> If no marsh area in cell, set marsh water elevation to open water stage
          Eh(j,2) = Es(j,2)
!          Qmarshmax(j) = 0.0
      endif      

!!>> Calculate daily values reported out to output files ! moved to hydro
!!	if(kday <= mmmd) then 
!	if(daystep == 1) then
!		ESMX(j,2)=ES(j,2)
!		ESMN(j,2)=ES(j,2)
!	!>> reset average water elevation value at start of day
!		ESAV(j,1) = ES(j,2)*dt/(3600.*24.)
!          EHAV(j,1) = EH(j,2)*dt/(3600.*24.)
!      else
!          ESMX(j,2)=max(ESMX(j,2),ES(j,2))
!          ESMN(j,2)=min(ESMN(j,2),ES(j,2))
!          !>> Update average elevation term by timestep's contribution to daily average 
!		ESAV(j,1)=ESAV(j,1) + ES(j,2)*dt/(3600.*24.)
!          EHAV(j,1) = EHAV(j,1) + EH(j,2)*dt/(3600.*24.)
!          dailyHW(j) = max(dailyHW(j),ES(j,2))
!          dailyLW(j) = min(dailyLW(j),ES(j,2))
!          EsRange(j,1)=ESMX(j,2)-ESMN(j,2)        
!	endif

!>> Calculate cumulative time marsh is flooded
!      if (ES(j,2) > (BedM(j)+dry_threshold)) then
      if (Eh(j,2) > (BedM(j)+dry_threshold)) then
          floodf(j)=floodf(j)+dt/(24.*3600.)	!accum days of flooding JAM Nov 13, 2010
      endif 


2222  FORMAT(<cells-1>(F20.2,','),F20.2)

3333  Format(5(1x,f15.4))
      return
      end

!c***********************End Subroutine for CelldQ***********************************************