!	Subroutine CelldTmp(QSalSUM,Dz,j,SalTRIBj,dref,Tres)
! kthr and kday now global parameters - no longer needed to be passed into subroutine      
!	Subroutine CelldTmp(QTmpSUM,j,kday,kthr,SalTRIBj,dref,Tres)
	  Subroutine CelldTmp(j,kday)
!JAM     c Salinity  computations ****************************
	
        use params

        real :: QRain, CSHEAT,rhoj,ddy1,ddy2,dddy,ake,aktmp,DTempw2,QTMPsum
        integer :: kdiv, ktrib,k,iab,jnb,j

        CSHEAT=4182.					!Specific heat
        rhoj=1000.*(1+S(j,1)/1000.)		!density

!>> Set minimum depth value (avoids div-by-zero errors)
        ddy1 = Es(j,1)-Bed(j)
        ddy2 = Es(j,2)-Bed(j)
	
!      if(ddy1 <= 0.1) then
!          dddy = 0.1
        if(ddy1 <= dry_threshold) then
            dddy = dry_threshold
        else
            dddy = ddy1
        endif
      
!        ake=5.0
        ake=26.5  !zw - 1/12/2024 based on MP2012 Ke averages
        cden=1./1000./24./3600.		!JAM Oct 2010 mm/d to m/s
        if(ddy1 <= dry_threshold) then
            aktmp=0
        else
            aktmp=ake/rhoj/CSHEAT/dddy !JAM Oct 2010 temperature heat conduction  *** will need calibration
        endif

        QTMPsum=0

!>> update  mass flux (QTMPsum) for diversion flows into compartment (diversions no longer modeled separately, but instead are treated as tributaries)
        do kdiv=1,Ndiv
            QTMPsum=QTMPsum - Qdiv(kdiv,kday)*TempMR(kday)*
     &	    Qmultdiv(j,kdiv)										!!!JAM Oct 2010
        enddo

!>> update mass flux (QTMPsum) for tributary flows into compartment      
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

!>> temperature source term from rainfall
  !ZW 1/13/2024 adding marsh area rainfall
        QRain = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
     &        -fpet*PET(kday,Jet(j)))*(As(j,1)+Ahf(j))*cden
!     &        -fpet*PET(kday,Jet(j)))*As(j,1)*cden
     
        QTMPsum=QTMPsum-QRain*ta(kday)            !openwater As 
      
!>> update mass flux (QTMPsum) for link flows into/out of compartment            
        do k=1,nlink2cell(j)
            if(icc(j,k) /= 0.0) then
                if(icc(j,k) < 0.0)then
                    jnb=jus(abs(icc(j,k)))
                else 
                    jnb=jds(abs(icc(j,k)))
                endif
            endif  
            iab=abs(icc(j,k))
!          call Temperature(mm,iab,jnb,j,k,QTMPsum)	!Temperature change computations
            call Temperature(iab,jnb,j,k,QTMPsum)	!Temperature change computations
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

!     openwater volume
        vol1 = 0.0
        vol2 = 0.0
        if( As(j,1) > 0 ) then                           ! check if there is openwater area
            if ( ddy1 > dry_threshold ) then               ! check if openwater was dry in previous timestep
                vol1 = ddy1*As(j,1)
            endif
            if ( ddy2 > dry_threshold ) then               ! check if openwater was dry in current timestep
                vol2 = ddy2*As(j,1)
            endif
        endif

! ZW-1/13/2024 adding marsh water volume
        ddym1 = Eh(j,1)-BedM(j)
        ddym2 = Eh(j,2)-BedM(j)

        marsh_vol1 = 0.0
        marsh_vol2 = 0.0
        if( Ahf(j) > 0 ) then                           ! check if there is marsh area
            if ( ddym1 > dry_threshold ) then               ! check if marsh was dry in previous timestep
                marsh_vol1 = ddym1*Ahf(j)
            endif
            if ( ddym2 > dry_threshold ) then               ! check if marsh was dry in current timestep
                marsh_vol2 = ddym2*Ahf(j)
            endif
        endif
        vol1 = vol1 + marsh_vol1
        vol2 = vol2 + marsh_vol2

!      Tempw(j,2) = ((Tempw(j,1)*As(j,1)*dddy - QTMPsum*dt)
!     &   / (As(j,1)*dddy +(Qsum_in(j)-Qsum_out(j)+QRain)*dt) ) + DTempw2

        if(vol2 > 0) then
            Tempw(j,2) = (Tempw(j,1)*vol1 - QTMPsum*dt)/vol2 + DTempw2
        else
            Tempw(j,2) = ta(kday)
        endif
!	Tempw(j,2)=Tempw(j,1)+DTempw1+DTempw2
	
     
      
!		if(Tempw(j,2).lt.TempMR(kday))Tempw(j,2)=TempMR(kday)

! low high pass filter
        if(Tempw(j,2).lt.2.0)Tempw(j,2)=2.0
!      if(Tempw(j,2).gt.36.)Tempw(j,2)=36.
!        if(Tempw(j,2).gt.tmpmax)Tempw(j,2)=tmpmax  !temperature instable
        if(Tempw(j,2).gt.ta(kday))Tempw(j,2)=ta(kday)

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
