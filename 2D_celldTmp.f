!	Subroutine CelldTmp(QSalSUM,Dz,j,SalTRIBj,dref,Tres)
! kthr and kday now global parameters - no longer needed to be passed into subroutine      
	Subroutine CelldTmp(QTmpSUM,j,kday,kthr,SalTRIBj,dref,Tres)
!JAM     c Salinity  computations ****************************
	
      use params

      real :: QRain, CSHEAT,rhoj,ddy,dddy,ake,aktmp,DTempw2,QTMPsum
      integer :: kdiv, ktrib,k,iab,jnb,j

	CSHEAT=4182.					!Specific heat
	rhoj=1000.*(1+S(j,1)/1000.)		!density

!>> Set minimum depth value (avoids div-by-zero errors)
      ddy= Es(j,1)-Bed(j)
	
!      if(ddy <= 0.1) then
!          dddy = 0.1
      if(ddy <= dry_threshold) then
          dddy = dry_threshold
      else
          dddy = ddy
      endif
      
      ake=5.0
      cden=1./1000./24./3600.		!JAM Oct 2010 mm/d to m/s
	aktmp=ake/rhoj/CSHEAT/dddy !JAM Oct 2010 temperature heat conduction  *** will need calibration

      QTMPsum=0

      do kdiv=1,Ndiv
	    QTMPsum=QTMPsum  -Qdiv(kdiv,kday)*TempMR(kday)*
     &	    Qmultdiv(j,kdiv)										!!!JAM Oct 2010
      enddo

	do ktrib=1,Ntrib
!>> If Qtrib is negative, flow is leaving system via tributary, use compartment temp if this is the case
          if (Qtrib(ktrib,kday) > 0.0) then
		    QTMPsum=QTMPsum-Qtrib(ktrib,kday)*TMtrib(ktrib,kday)*
     &			            Qmult(j,ktrib)
          else
              QTMPsum = QTMPsum-Qtrib(ktrib,kday)*TempW(j,1)*
     &                    Qmult(j,ktrib)
          endif
      enddo
      QRain = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
     &        -fpet*PET(kday,Jet(j)))*As(j,1)*cden
     
	QTMPsum=QTMPsum-QRain*ta(kday)            !openwater As 
      
      
      do k=1,nlink2cell(j)
          if(icc(j,k) /= 0.0) then
		    if(icc(j,k) < 0.0)then
		        jnb=jus(abs(icc(j,k)))
              else 
		        jnb=jds(abs(icc(j,k)))
              endif
          endif  
          iab=abs(icc(j,k))
		call Temperature(mm,iab,jnb,j,k,QTMPsum)	!Temperature change computations
      enddo


!      DTempw1=  -QTMPsum/(As(j,1)*dddy)*dt
!     &	-Dz*Tempw(j,1)/dddy

!>> New temperature equations:
! T2*V2 = T1*V1 + Tin*Vin + Tsource

! V2 = As(j,1)*dddy - QSum_out(j)*dt + Qsum_in(j)*dt + Qrain*dt
! T1 = Tempw(j,1)
! V1 = As(j,1)*dddy
! Tin*Vin = QTMPsum*dt
! Tsource = DTempw2

!>> Calculate temperature source from equilbrium tempertaure at air-water interface

      DTempw2=aktmp*(Tempe(j,kday)-Tempw(j,1))*dt			!JAM Oct 2010
      
      ! debug for v23.4.2
      !if(isNan(QTMPsum)) then
      !	  QTMPsum = 0.0
      !endif

      Tempw(j,2) = ((Tempw(j,1)*As(j,1)*dddy - QTMPsum*dt)
     &   / (As(j,1)*dddy +(Qsum_in(j)-Qsum_out(j)+QRain)*dt) ) + DTempw2

!	Tempw(j,2)=Tempw(j,1)+DTempw1+DTempw2
	
     
      
!		if(Tempw(j,2).lt.TempMR(kday))Tempw(j,2)=TempMR(kday)

! low high pass filter
      if(Tempw(j,2).lt.2.0)Tempw(j,2)=2.0
      if(Tempw(j,2).gt.36.)Tempw(j,2)=36.

      if(isNan(Tempw(j,2))) then      ! YW
          write(*,*) 'j:',j
          write(*,*) 'Tempw(j,2):',Tempw(j,2)
          write(*,*) 'Tempw(j,1):',Tempw(j,1)
          write(*,*) 'As(j,1):',As(j,1)
          write(*,*) 'dddy:',dddy
          write(*,*) 'QTMPsum:',QTMPsum
          write(*,*) 'dt:',dt
          write(*,*) 'Qsum_in(j):',Qsum_in(j)
          write(*,*) 'Qsum_out(j):',Qsum_out(j)
          write(*,*) 'QRain:',QRain
          write(*,*) 'DTempw2:',DTempw2
          write(*,*) 'aktmp:',aktmp
          write(*,*) 'kday:',kday
          write(*,*) 'Tempe(j,kday):',Tempe(j,kday)   
          stop
	  !write(*,*) 'setting water temp to equil temp due to NaN'
	  !Tempw(j,2) = Tempe(j,kday)
      endif

!	fsal=(1+S(j,2)/35) !-EDW not used anywhere									!salinity correction on Vs
!	Temph(j,2)=Tempw(j,2)+0.5							!JKS 10/31/13

      return 
	end
!***********************End Subroutine for Change in Cell Temperature*******JAM Oct 2010********
