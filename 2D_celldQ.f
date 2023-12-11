!     Subroutine CelldQ(QSUM,Dz,j,fcrop,mm)
	
! day, dday, and kday now global parameters - no longer needed to be passed into subroutine      
 	Subroutine CelldQ(j,kday,fcrop,mm,dday)			!Solves Continuity for cell j 

      
!> @param[out]        Qsum_out(j)     sum of all flows leaving compartment
!> @param[out]        Qsum_in(j)      sum of all flows entering compartment
!> @param[out]        Qsum_abs(j)     magnitude of all flows entering AND leaving compartment
      
	use params      
      real(sp) :: Qlink,Dzhlim,Elevel,flo_trib,flo_div,mindz

	Qsum=0.0						!JAM Oct 2010
	Qsumh=0.0
      Qlink = 0.0
	cden=1./1000./24./3600.		! mm/d to m/s conversion

!	if((mm>=626160).AND.(j==458))then
!		write(1,*)'time step = ',mm, 'compartment = ',j,'At start: ',Qsumh
!	endif

!>> initialize directional components of cumulative flowrates that are used for sand transport equations in vanRijn
      Qsum_out(j) = 0.0
      Qsum_in(j) = 0.0
      Qsum_abs(j) = 0.0
!      Qsum_out_link(j) = 0.0
!>> Update cumulative open water flow rates in compartment from input boundary tributary flows      
!>> sign convention on open water flow = positive is flow out of compartment
      do jn = 1,ntrib
		if (isnan(Qmult(j,jn))) then
              Qmult(j,jn) = 0.0
          endif
          flo_trib = Qtrib(jn,kday)*Qmult(j,jn)
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
          Qsum = Qsum - flo_div
          Qsum_abs(j) = Qsum_abs(j) + abs(flo_div)
          if (flo_div >= 0.0) then
              Qsum_in(j) = Qsum_in(j) + flo_div
          else
              Qsum_out(j) = Qsum_out(j) + abs(flo_div)
          endif
      enddo

c correction JAM March 26 2007  ---correction JAM Aug 10 090.50093
	fpc= percent(j)*(1.-fcrop)/100.+fcrop               ! multiplier on PET for marsh areas ; if percent(j)=10 and fcrop=0.5, then the marsh area will evaporate 0.55*PET for the day
	soilm = max(0.0001, esho(j)-BedM(j))				!JAM Oct 2010
	shh   = max(0.0001, Eh(j,1)-BedM(j))				!JAM Oct 2010
	rhh   = max(0.0001, shh/soilm)			            !JAM Oct 2010

!>> Excess rainfall runoff on marsh
      Qhhf=Ahf(j)*(Rain(kday,jrain(j))-PET(kday,Jet(j)))!*Max(1.,rhh*rhh))	!JAM Oct 2010

!>> !YW! ignore PET at low marsh water level to avoid Eh too low
!      if((Eh(j,1)-BedM(j))>0.1) then
!          Qhhf=Ahf(j)*(Rain(kday,jrain(j))
!     &    -PET(kday,Jet(j))*Max(1.,rhh*rhh))	!JAM Oct 2010
!      else
!          Qhhf=Ahf(j)*Rain(kday,jrain(j))
!      endif

      Ahmf=Ahydro(j)-Ahf(j)

!>> Update cumulative flow rate in open water based on excess rainfall runoff on open water area
!>> sign convention on open water flow = positive is flow out of compartment
      ! cden: mm/day to m/s conversion factor
      ! Rain: mm/day rainfall
      ! PET:  mm/day potential ET
      ! fpet: 1=use PET input data; 0: use average ET value
      ! ETA:  average ET value to use if fpet = 0
      Qsum=Qsum-(Rain(kday,jrain(j))-(1-fpet)*ETA(Jet(j))		!openwater As 
     &	  -fpet*PET(kday,Jet(j)))*As(j,1)*cden
!>> Update cumulative flow rate in marsh based on excess rainfall runoff on marsh area
!>> sign convention on marsh flow = positive flow is from marsh to open water
	Qupld=(Qhhf+max(0.0,(Rain(kday,jrain(j))
     &	 -PET(kday,Jet(j))*fpc))*Ahmf)*cden	 
!      Qsumh=Qsumh-(Qhhf+max(0.0,(Rain(kday,jrain(j))
!     &	 -PET(kday,Jet(j))*fpc))*Ahmf)*cden								!Runoff>0    !JAM Oct 2010          

!>> modified zw 04/29/2020 to deal with upland compartments w/o marsh area
      if(Ahf(j) > 0.0) then	
          Qsumh=Qsumh-Qupld
      else
          Qsum=Qsum-Qupld
      endif

!	if((mm>=626160).AND.(j==458))then
!		write(1,*)'time step = ',mm, 'compartment = ',j,'Add runoff:',Qsumh
!	endif
!>> Update cumulative marsh flowrate for marsh-to-open water exchange flow (calculated via Kadlec Knight in hydrod)      
!>> sign convention on marsh flow = positive flow is from marsh to open water
      Qsumh=Qsumh-Qmarsh(j,2)											!Qsumh is net of marsh flows 
!	if((mm>=626160).AND.(j==458))then
!		write(1,*)'time step = ',mm, 'compartment = ',j,'Add Qmarsh:',Qsumh
!	endif
!	 (sign convention outward normal is positive and inward is negative. 
cccccccccc cjam collects flow from connecting links

!     do k=1,13      
!      do k=1,maxconnect						!more connection: note max number of connected links is 11 to one cell (non-marsh) ***

!>> collect flow from connecting links
!>> sign convention on link flows: positive is flow out of compartment
      do k=1,nlink2cell(j)
          if(abs(icc(j,k)) /= 0) then
              Qlink = sicc(j,k)*Q(abs(icc(j,k)),2)
!>> if link type is marsh overland flow type, add flow to marsh flow sum
!>> if compartment has no marsh area, add marsh overland flow in compartment to water area
!>> if  marsh link flow is negative, flow is entering marsh from neighboring marsh
              if (linkt(abs(icc(j,k))) == 8) then
                  if (Ahf(j) > 0) then
                      Qsumh = Qsumh + Qlink
                  else
                      Qsum = Qsum + Qlink
                  endif
!				if((mm>=626160).AND.(j==458))then
!					write(1,*)'time step = ',mm, 'compartment = ',j,abs(icc(j,k)),Qlink
!				endif
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
!                      
!                  
                  
                  
                  
!	if((mm>=626160).AND.(j==458))then
!		write(1,*)'time step = ',mm, 'compartment = ',j,'Add Qlink:',Qsumh
!	endif
      
!      if (isNan(Qmarsh(j,1))) then
!          Qmarsh(j,1) = 0.0
!          write(*,*) 'Qmarsh in ',j,' set to 0.0 @ step ',mm
!      endif
!>> add marsh exchange flow to open water cumulative flow
!>> if positive, flow is open water-to-marsh flow	
      Qsum=Qsum+Qmarsh(j,2)				!JAM Oct 2010 Transfer from marsh runoff or ebb flow NOTE +ve means outflow from open water cell.

!JAM   change in elevation for cell j

!>> Calculate change in open water stage
!>> negative flows are into compartment and will result in positive deltaZ
      Dz=((-Qsum)/As(j,1))*dt				                !Euler method,3.

      sndz = 1.0
      if (Dz < 0) then
          sndz = -1.0
      endif

!>> do not allow Dz to be smaller than a value that would result in a change in water level of less than 1 mm over an entire day
      !mindz = 0.00000035
      !if (abs(Dz) < mindz) then
      !    Dz = 0.0
      !endif
      !
      if(abs(Dz) > maxdz) then
          Dz = sndz*oscilflag   ! Dz = sndz*maxdz
          write(*,898) 'Large deltaZ in compartment ',j,'timestep=',mm,'
     & deltaZ set to max value allowed in RunControlR.dat'
          write(1,898) 'Large deltaZ in compartment ',j,'timestep=',mm,'
     & deltaZ set to max value allowed in RunControlR.dat'
          !pause
898   Format(x,A,I0,x,A,x,I0,x,A)
      endif

!>> Update water elevation of open water portion - don't let water level drop below bed elevation
      Es(j,2)=max(Es(j,1) + Dz,Bed(j)+0.0001)				                !Stage in storage cells (m)
      
      mmmd=float(mm+1)*dt/3600./24.                       !mmmd is 
      Dzz = sndz*min(abs(Dz),oscilflag)
      if (abs(dzz) == oscilflag) then
          write(*,*) 'Comp=',j,'timestep=',mm,'dz=',dz
          write(*,*) 'Qsum=',Qsum,'Qmarsh=',Qmarsh(j,2)
          write(*,*) 'OWstg(t-1)=',Es(j,1),'OWstg(t)=',Es(j,2)
          write(*,*) 'Mstg(t-1)=',Eh(j,1),'Mstg(t)=',Eh(j,2)
          write(1,*) 'Comp=',j,'timestep=',mm,'dz=',dz
          write(1,*) 'Qsum=',Qsum,'Qmarsh=',Qmarsh(j,2)
          write(1,*) 'OWstg(t-1)=',Es(j,1),'OWstg(t)=',Es(j,2)
          write(1,*) 'Mstg(t-1)=',Eh(j,1),'Mstg(t)=',Eh(j,2)
      endif

!>> Calculate change in marsh stage
!>> Negative marsh flow is into marsh and will result in positive deltaZmarsh
      if (Ahf(j) > 0.0) then
          Dzh=(-Qsumh)/Ahf(j)*dt

!>> Update marsh water elevation
          Eh(j,2)= Max(Eh(j,1)+Dzh, BedM(j))		!JAM Oct 2010 Marsh stage (m) 
!>> Determine level where open water stage and marsh stage will be equal in compartment
          Elevel = (As(j,1)*Es(j,2)+Ahf(j)*Eh(j,2)) / (As(j,1)+Ahf(j))
!>> Determine change in marsh level to reach equilbrium with open water marsh stage
          Dzhlim = max(Elevel,BedM(j)) - Eh(j,2)
!>> Calculate magnitude of flowrate needed into/out of marsh - directionality assigned in hydrod
          Qmarshmax(j) = abs(Dzhlim)*Ahf(j)/dt
              
      else
!>> If no marsh area in cell, set marsh water elevation to open water stage
          Eh(j,2) = Es(j,2)
          Qmarshmax(j) = 0.0
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
      if (ES(j,2) > (BedM(j)+dry_threshold)) then
	    floodf(j)=Floodf(j)+dt/(24.*3600.)	!accum days of flooding JAM Nov 13, 2010
	endif 


2222	FORMAT(<cells-1>(F20.2,','),F20.2)

	!Zw added 3/13/2015 for compartment 458 instability diagonis at dt=30s
	!need to be deleted or commented out after the diagonis is done
!	if((mm>=626160).AND.(j==458))then
!		write(1,*)'At time step = ',mm, 'for compartment = ',j
!		write(1,3333)Q(2203,1),Q(2207,1),Q(2209,1),Q(2208,1),Qmarsh(j,1)
!		write(1,3333)Qsumh,Qhhf*cden,Qupld,Rain(kday,jrain(j)),rhh
!	endif
3333  Format(5(1x,f15.4))
	return
	end

c***********************End Subroutine for CelldQ***********************************************