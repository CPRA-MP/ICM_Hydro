      Subroutine CelldChem(j,kday,ichem)		!Chemical  computations

	use params      
      

      real :: QRain
      
	QChemSUM(ichem) = 0.0
      QChemSUManth(ichem) = 0.0
      QChemSUMtrib(ichem) = 0.0
      QChemSUMdiv(ichem) = 0.0
      QChemSUMflows(ichem) = 0.0
      QChemSUMatm(ichem) = 0.0
      
!>> water depth in compartment      
      dyy=max((Es(j,1)-Bed(j)),0.01)

!>> Anthropogenic loads 
      if(ichem == 1) then
          if(dyy > 0.01) then
!          QChemSUM(ichem)=QChemSUM(ichem)-AnthL(j)								! kg/d N-NO3  JAM April 3, 2011
              QChemSUManth(ichem) = -AnthL(j)*1000.0/(24.0*3600.0)			!zw 4/28/2015 AnthL kg/d to g/s
          endif
      endif
      
      if(ichem == 5) then
          if (dyy > 0.01) then
!          QChemSUM(ichem)=QChemSUM(ichem)-AnthL(j)/7.1                          
              QChemSUManth(ichem) = -AnthL(j)*1000.0/(24.0*3600.0)/7.1     !zw 4/28/2015 AnthL kg/d to g/s                     
          endif
      endif

!>> Determine tributary contributions to cumulative WQ load	
      do ktrib=1,Ntrib
!>> If Qchem is negative, flow is leaving system via tributary, use compartment WQ concentration if this is the case
          if (QChem(ktrib,ichem,kday) > 0.0) then
!			QChemSUM(ichem)=QChemSUM(ichem)-QChem(ktrib,ichem,kday)*
!     &			Qmult(j,ktrib)
              QChemSUMtrib(ichem) = QChemSUMtrib(ichem)
     &			-QChem(ktrib,ichem,kday)*Qmult(j,ktrib)		!-EDW 6/4/2015 QChem is now in g/s!zw 4/28/2015 QChem kg/s to g/s	(QChem is converted from kg/d to kg/s in infile.f)
          else
              QChemSUMtrib(ichem) = QChemSUMtrib(ichem)
     &            - Qtrib(ktrib,kday)*Qmult(j,ktrib)*Chem(j,ichem,1)
          endif
	enddo															! end trib chem   g/s
!>> Determine diversion contributions to cumulative WQ load      
      do kdiv=1,Ndiv
!			QChemSUM(ichem)=QChemSUM(ichem)-QChemdiv(kdiv,ichem,kday)*	! JAM Feb 22, 2010
!     &			Qmultdiv(j,kdiv)
			QChemSUMdiv(ichem) = QChemSUMdiv(ichem)
     &			    -QChemdiv(kdiv,ichem,kday)*Qmultdiv(j,kdiv)     !-EDW 6/4/2015 QChemdiv is now in g/s !zw 4/28/2015 QChem kg/s to g/s	(QChemdiv is converted from kg/d to kg/s in infile.f)
	enddo															
													! note max number of connected links is 11 *** JAM Oct 2010
      do k=1,nlink2cell(j)
          if(icc(j,k) /= 0) then
	    	if (icc(j,k) < 0) then
			    jnb=jus(abs(icc(j,k)))
		    else
			    jnb=jds(abs(icc(j,k)))
		    endif  
          endif
 		iab=abs(icc(j,k))

!>> call chemical subroutine (mass balance)
	    call chemical(mm,iab,jnb,j,k,ichem)
      enddo															! k do loop neighbouring cell contributions
      
      if (dyy > 0.01) then
          QChemSUMatm(ichem) = -QAtm(1,ichem,kday)
     &                     *As(j,1)/(1000000.)*1000.0	  ! zw 4/28/2015 QAtm kg/s/km2 to g/s/km2 (QAtm is converted from kg/d/km2 to kg/s/km2 in infile.f)
      endif
!      QChemSUM(ichem)=QChemSUM(ichem) 
!     &	- QAtm(1,ichem,kday)*(As(j,1)+AHydro(j))/(1000000.)			! kg/day

      QChemSUM(ichem) = QChemSUManth(ichem)
     &                + QChemSUMtrib(ichem)
     &                + QChemSUMdiv(ichem)
     &                + QChemSUMflows(ichem)
     &                + QChemSUMatm(ichem)

     
     
c  chemical change computations  *********************************
      mex=9
      DChemSUM = 0.

! WQ subroutines that start with 'd' (e.g. dNO3) were updated to match the 2012 equations used in AA and CP
! these subroutine names do not necessarily match their .f filename (e.g. dTP is located in TP.f)
      if (dyy > 0.01) then
          if(ichem == 1) then
              call dNO3(DChemSum,ichem,j)
          elseif(ichem == 2) then
              call dNH4(DChemSum,ichem,j)
            elseif(ichem == 5) then
              call dTP(DChemSum,ichem,j)      ! this used to be TP, now TIP since updated WQ equations
            elseif(ichem == 7) then
              call DOx(DChemSum,ichem,mex,j,k,kday)
          elseif(ichem == 8) then 
              call dLivA(DChemSum,ichem,j)
          elseif(ichem == 9) then
                call dDeadA(DChemSum,ichem,j)      
          elseif(ichem == 10) then
              call dDON(DChemSum,ichem,j)
          elseif (ichem == 11) then
              call dDOP(DChemSum,ichem,j)     
          endif
      endif

!>> Calculate contribution to flow from precip and evapotranspiration      
      QRain = (Rain(kday,Jrain(j))-(1-fpet)*ETA(Jet(j))
     &        -fpet*PET(kday,Jet(j)))*As(j,1)*cden      
      if (dyy <= 0.01) then  
          QRain = max(QRain,0.0)
      endif

      
!>> Transport equation for WQ constituents, dChem      
      Chem(j,ichem,2) = (Chem(j,ichem,1)*As(j,1)*dyy-QChemSUM(ichem)*(dt/max(1,NTs2_ICM))
     &                + DChemSUM*(dt/max(1,NTs2_ICM))/(24.0*3600.0)*As(j,1)*dyy) 
     &                / (As(j,1)*dyy +(Qsum_in(j)-Qsum_out(j)+QRain)*(dt/max(1,NTs2_ICM)))      
      
!      DChem(ichem)=  -1.0*QChemSUM(ichem)/(As(j,1)*dyy)*dt
!     &	-Dz*Chem(j,ichem,1)/dyy+DChemSUM*(dt/(24.0*3600.0))   !zw 4/28/2015 change dt from sec to day in DChemSUM term because the paramters in the source/sink terms are all in 1/day
!     &	-Dz*Chem(j,ichem,1)/dyy+DChemSUM*dt

!>> Calculate new WQ concentration from dChem
      Chem(j,ichem,2)=min(1000.0,max(Chem(j,ichem,2),0.0))
!      Chem(j,ichem,2)=min(1000.0,max(Chem(j,ichem,1)+DChem(ichem),0.0))			  !zw 4/28/2015 if Chem unit is mg/L (=g/m3), then zw changes are correct; otherwise, need to check the WQ units 
      !
      !if (chem(j,ichem,2)>1000.0) then
      !    write(*,*) 'ICMID:',j,'ichem:',ichem
      !    write(*,*) '1:',chem(j,ichem,1),dchem(ichem),chem(j,ichem,2)
      !    write(*,*) '2:',QChemSUM(ichem)
      !    write(*,*) '3:',QChemSUManth(ichem)
      !    write(*,*) '4:',QChemSUMtrib(ichem)
      !    write(*,*) '5:',QChemSUMdiv(ichem)
      !    write(*,*) '6:',QChemSUMflows(ichem)
      !    write(*,*) '7:',QChemSUMatm(ichem)
      !    pause
      !
      !endif

	return 
	end

c***********************End Subroutine for Change in Cell Chemistry*******JAM Oct 2010**********