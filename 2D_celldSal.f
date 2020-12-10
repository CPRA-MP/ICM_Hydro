!      Subroutine CelldSal(QSalSUM,Dz,j,SalTRIBj,dref,Tres)

! kthr and kday now global parameters - no longer needed to be passed into subroutine      
	Subroutine CelldSal(QSalSUM,j,kday,k	,SalTRIBj,dref,Tres)
! kthr was passed into this subroutine as k which was then overwritten during DO loop
cjam     c Salinity  computations ****************************

	use params
      
      real :: Saltrib

!>> Set minimum depth value (avoids div-by-zero errors)
      ddy= Es(j,1)-Bed(j)
	
      if(ddy <= 0.01) then
          dddy = 0.01
      else
          dddy = ddy
      endif
      
      ddym= Eh(j,1)-BedM(j)
	
      if(ddym <= 0.01) then
          dddym = 0.01
      else
          dddym = ddym
      endif
      
      QSalsum=0

      do ktrib=1,Ntrib
			
!>> set salinity in tributary to default freshwater salinity value (assigned in hydrod)
          Saltrib = Saltribj
!>> if tributary flow is negative, use compartment salinity concentration instead of default tributary salinity concentration
          if (Qtrib(ktrib,kday) < 0.0) then
              Saltrib = S(j,1)
          endif             
          QSalsum=QSalsum-Qtrib(ktrib,kday)*Saltrib*
     &			Qmult(j,ktrib)
      enddo

      do kdiv=1,Ndiv
          QSalsum=QSalsum-Qdiv(kdiv,kday)*0.15*	!JAM jan 09 11
     &			Qmultdiv(j,kdiv)
      enddo
!		do k=1,13 									! note max number of connected links is 5 ***
!		do k=1,maxconnect

      do k=1,nlink2cell(j)
          if(icc(j,k) /= 0) then
              if(icc(j,k) < 0) then
		        jnb=jus(abs(icc(j,k)))
		    elseif(icc(j,k) > 0) then
                  jnb=jds(abs(icc(j,k)))
	        endif  
          endif
          
	    iab=abs(icc(j,k))
    
	    call salinity(mm,iab,jnb,j,k,Qsalsum)
      enddo										!k do loop   

      
      if (Qmarsh(j,2) > 0.0) then                !YW! combining marsh and OW volumne
          QSalSum=QSalSum + Qmarsh(j,2)*S(j,1)
      endif

c  salinity change computations  *********************************
      QRain = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
     &        -fpet*PET(kday,Jet(j)))*As(j,1)*cden      
      
      if (dddy <= 0.1) then  
          QRain = max(QRain,0.0)
      endif

      QRainm = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
!     &        -fpet*PET(kday,Jet(j)))*Ahf(j)*cden      
     &        -fpet*PET(kday,Jet(j)))*Ahydro(j)*cden
      
      if (dddym <= 0.1) then  
          QRainm = max(QRainm,0.0)
      endif

!>> original equation
      S(j,2) = (S(j,1)*As(j,1)*dddy-QSalsum*dt)
     &      /(As(j,1)*dddy+(Qsum_in(j)-Qsum_out(j)+QRain)*dt)


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
!      ddym1=Eh(j,1)-BedM(j)
!      ddym2=Eh(j,2)-BedM(j)
!      if(Ahf(j)>0 .AND. (ddym1>0.1 .AND. ddym2>0.1))then
!	      vol1=(Es(j,1)-Bed(j))*As(j,1)+(Eh(j,1)-BedM(j))*Ahf(j)
!	      vol2=(Es(j,2)-Bed(j))*As(j,1)+(Eh(j,2)-BedM(j))*Ahf(j)
!	  else
!	      vol1=(Es(j,1)-Bed(j))*As(j,1)
!	      vol2=(Es(j,2)-Bed(j))*As(j,1)
!	  endif
	  
!	  if((Es(j,2)-Bed(j))>0.1)then
!	      S(j,2)= (S(j,1)*vol1-QSalsum*dt)/vol2
!     else
!          S(j,2) = S(j,1)
!      endif
      
!>> equation for MP2017 to avoid salinity spike
!      if (dddy > 0.1) then   
!          DSal =  -QSalsum/(As(j,1)*dddy)*dt
!     &	    -Dz*S(j,1)/dddy
!      else
!          DSal = 0.0
!      endif
!      
!	S(j,2)=S(j,1)+DSal
         
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
      if(S(j,2) < 0.10) then
          S(j,2) = 0.10
      elseif (S(j,2) > 36.) then
	    S(j,2)=36.
      endif

 !>> Test setting salinity to 0.1 when the compartment is dry   Nov 2020  YW
      if (dddy == 0.01) then
          S(j,2) = 0.10
      endif

	Sh(j,2)=S(j,2)  


!!>> Calculate daily average values reported out to output files ! moved to hydro
!	if(daystep == 1) then
!	!>> reset average salinity value at start of day
!		SALAV(j) = S(j,2)*dt/(3600.*24.)
!      else
!          !>> Update average salinity by timestep's contribution to daily average 
!		SALAV(j)=SALAV(j) + S(j,2)*dt/(3600.*24.)
!	endif

      return 
	end
c***********************End Subroutine for Change in Cell Salinity*******JAM Oct 2010***********