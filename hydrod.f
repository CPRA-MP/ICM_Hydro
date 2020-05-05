	Subroutine hydrod(mm)			!face densities from node densities

      use params

      implicit none

      integer :: ichem
      integer :: day_int,i,j,jj,jjk,k,kk,kl,klk,mm
      integer :: sedclass,kthr,kmon,kday,simhr,khr
      real :: dday,ttday
      real :: Rh,AK,Res,Resist,Deta,Detah,Dmarsh,Sf,sn,dkd,linkdepth,dhr
      real :: time,thour,tmon,tday,tdayj,hday
      real :: g0,g1,g2,g3,g4,g5,g6
      real :: thourj
      real :: tt0,tt1,tt2,tt3,tt4,tt5,tt6,agulf3,Sseason
      real :: shour,sday,f1,f2,f3,t1,t2,t3,tlag,aset
      real :: fcrop,Tres,Vettl,dref,CSSTRIBj,Saltribj,dmon,dmod
      real :: akL,akns
      real :: QSalSum,QTmpSum,QMarshKK
      real :: rca,rna,rnd,rcd,rcn
      real :: Q_filter1,Q_filter2                                               ! yw flow filter
      integer :: flag_skip                                                      ! yw skip flag
	  integer :: flag_offbc(Cells)  !zw offshore bc cells flag 04/07/2020
      real :: dkd_h,Marsh_ruf,MarshRh,MarshAch                                  ! edw new parameters for replacing Kadlec-Knight with Manning's
      real :: MarshRes,MarshResist,QMarshMann                                   ! edw new parameters for replacing Kadlec-Knight with Manning's

!parameters for flow by linktype
      real :: ach,AET,avdep,delp,deflength,delh,dv
      real :: beta,eqC,orarea,orc
      real :: hf,htdn,htupf,htup
      real :: infilrate,invdown,invup
      real :: p,subp,pexcess,marshedge,marshl
      real :: Qmax,Qpump,Qrunoff,reg_r,reg_s,reg_p,reg_fs
      real :: ruf,w_k,w_ksub,wid
      real :: sDetah,volavailable
      integer :: downN,upN,pumpon

! parameters for Atchafalya River diversion for MP project runs
      integer :: Atch_US_link, Atch_DS_link,BayouShaffer_link,Div_link
      real :: Atch_US_Q,BayouShaffer_Qold,Div_Q

			!c     time in seconds
			      time=float(mm)*dt           ! elapsed time
			      thour=mm*dt/3600            ! elapsed time in hours
				tmon=thour/730.				!!added JAM June 23, 2009  -- 730.5--> 730 June 26, 2009
				kthr=int(thour+1)
				kmon=ifix(tmon)				!+1    !!added JAM June 23, 2009
			!
			!! calculate various versions of time to be used as flags throughout program
				kday = ifix(time/3600./24.)+1   !convert time to integer days !BUG! Why is kday day+1?
			      day = time/3600./24.
			      dday=day-int(day)           ! decimal portion of day, dday=0.0 at 0:00 (midnight)
				hday=day-int(day+0.5)       ! decimal portion of day normalized to noon, hday=0.0 at 12:00 (noon)
			      simhr = floor(dday*24.)     ! hour of the simulation day, in integer
			      dhr = thour -int(thour)     !decimal portion of hour , dhr - 0.0 at XX:00

						dmon=kmon-int(tmon)
			      dmod=tmon-int(tmon)

			! Moved to start of main.f time looping and added variables to global parameters list      !-EDW

c Temporary values
				Cp =0.5
				fcrop=0.5          !0.1  !0.59                        !potential coef
				Tres = 3600.
			! 	Vsettl=8/24/3600 !-EDW not used
				dref=4.0
			!>> default CSS and salinity concentrations in tributary flow
			      CSSTRIBj=25.
				Saltribj=0.205


ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc OPEN WATER BOUNDARY CONDITIONS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc		do jj=101,mds+101-1
      do jjk=1,mds  !AMc 8 oct 2013
	    jj=KBC(jjk) !AMc 8 oct 2013
		tday= t/24/3600
		tdayj=tday-(int(tday/365.25))*365.25					!cal julian day JAM Nov 2010

cc _________________ JAM Nov 2010 revised
		g0 = -0.07  !m
		g1 = -0.000233005
		g2 = 0.000000326479
		g3 = -0.000000000146371
		g4 =  0.0000000000000299401
		g5 = -0.00000000000000000283698
		g6 =  0.000000000000000000000100482
		thourj=(tdayj)*24
		tt0 = 1
		tt1 = thourj
		tt2 = tt1*thourj
		tt3 = tt2*thourj
		tt4 = tt3*thourj
		tt5 = tt4*thourj
		tt6 = tt5*thourj
		agulf3= g0+g1*tt1+g2*tt2+g3*tt3+g4*tt4+g5*tt5+g6*tt6	!JAM Nov 2010
		Sseason=agulf3*10.
		S(jj,1) = SBC(jj)+SSeason                       !  0.6*SEARD*t/0.8+   Dec 2011
!		age(jj,1)=0.0											!JAM Oct 2010

          Chem(jj,1,1) = BCNO3(jj) !YW removing dividing by 1000. assuming input in mg/l
          Chem(jj,2,1) = BCNH4(jj)
          Chem(jj,3,1) = Chem(jj,1,1) + Chem(jj,2,1)  !zw 4/28/2015 add DIN=NO3+NH4 at offshore BCs
          Chem(jj,4,1) = BCON(jj)
          Chem(jj,5,1) = BCTP(jj)
          Chem(jj,6,1) = BCTOC(jj)
          Chem(jj,7,1) = BCDO(jj)
          Chem(jj,8,1) = BCLA(jj)
          Chem(jj,9,1) = BCDA(jj)
          Chem(jj,12,1) = Chem(jj,5,1)*0.9
          Chem(jj,10,1) = 0.2               !YW originally 0.0, average calculation see infile New 0.4
          Chem(jj,11,1) = 0.012             !YW  TP*ParDOP originally 0.0 New 0.012 0.05
          Chem(jj,13,1) = BCTOC(jj)*0.04
          Chem(jj,14,1) = BCTOC(jj)*0.025   !!marsh POP JAM April 16, 2011


!		Age(jj,1)=BCage(jj)
!		Chem(jj,1,2)=BCNO3(jj)/1000.
!		Chem(jj,2,2)=BCNH4(jj)/1000.
!		Chem(jj,4,2)=BCON(jj)/1000.
!		Chem(jj,5,2)=BCTP(jj)/1000.
!		Chem(jj,6,2)=BCDO(jj)/1000.
!		Chem(jj,7,2)=BCTOC(jj)/1000.
!		Chem(jj,6,2)=BCTOC(jj)/1000.
!          Chem(jj,7,2)=BCDO(jj)/1000.
!          Chem(jj,8,2)=BCLA(jj)/1000.
!		Chem(jj,9,2)=BCDA(jj)/1000.
!		Chem(jj,12,2)=Chem(jj,5,1)*0.9
!		Chem(jj,10,2)=0.0
!		Chem(jj,11,2)=0.0
!		Chem(jj,13,2)=BCTOC(jj)/1000.*0.04
!		Chem(jj,14,2)=BCTOC(jj)/1000.*0.025   !!marsh POP JAM April 16, 2011
!         Age(jj,2)=BCage(jj)
		do ichem=1,14	 !zw 4/28/2015 added to replace above statements
	       Chem(jj,ichem,2)=Chem(jj,ichem,1)
	    enddo
			S(jj,2)=S(jj,1) !zw added 04/07/2020
	enddo

!MP2023 zw moved tide bc and css bc update to here 04/07/2020

c*********************TIDE BC*******************************************************************

c     These parameters should be moved to input to start at any time of year ** later
********tide information, surges, periods and phase angles

	      shour=0.0
		sday=0.0
		f1=shour*3600.						!daily phase
		f2=(sday/28.)
		f2=(f2-int(f2))*28.*24.*3600.		!lunar phase
		f3=(sday/365.25)*12.*30.*24.*3600.	!Gulf phase
		t1=23.5*3600.						!daily period
		t2=672.*3600.						!lunar period
		t3=4320.*3600.						!Gulf period
	      tlag=0.								!Tide lag time along the eastern boundary in the Gulf
		aset=0.01

c***********************Open Boundary Conditions I.GEORGIOU/JAMc********************************
c***********************************************************************************************

	!      do jjk=1,Mds !AMc 8 oct 2013 revised Boundary
	!	    jj=KBC(jjk)
	!          Call TideBC(jj,jjk)
	!! kthr and kday now global parameters - no longer needed to be passed into subroutine
	!!          Call TideBC(jj,kthr,kday)
	!      enddo
	!>> Update water level for boundary condition cells - this subroutine will loop through all boundary condition compartments
	      Call TideBC
c***********************END TIDE****************************************************************

c Open Boundary Sed Conc.           !needs revision to reflect MR TSS JAM Oct 2010
CCCCCCCCCCCCCCCCCCCC AMc
cc		Css(101,2)=(BCTSS(101)+fcbc*50.*sin(pi*wd(kday)/180.))/3.
cc		Css(102,1)=(BCTSS(101)+fcbc*40.*sin(pi*wd(kday)/180.))/3.
cc		if(mds.gt.2) then
cc			do ii=103,110
cc				Css(ii,2)=75.+ fcbc*(1-float(ii-103)/float(110-103))
cc     &				*150.*sin(pi*wd(kday)/180.)
cc			enddo
cc		endif
cccccc  AMc 8 oct 2013  change TSS BC
cc		Css(101,2)=(BCTSS(101)+fcbc*50.*sin(pi*wd(kday)/180.))/3.
cc		Css(102,1)=(BCTSS(101)+fcbc*40.*sin(pi*wd(kday)/180.))/3.
cc		if(mds.gt.2) then


	!>> seasonal adjustment of boundary condition sediment data
	!>> Assume Boundary condition sediment is evenly distributed amongst clay and silt size classes and half of the clay is flocced.
	      BCSedRatio(1) = 0.0
	      BCSedRatio(2) = 1./2.
	      BCSedRatio(3) = 1./4.
	      BCSedRatio(4) = 1./4.

	 do jjk=1,mds
	    jj=KBC(jjk)
			do kk=1,4
			  !Css(jj,2,kk)=BCTSS(jj)*(1.+ 0.5*sin(pi*wd(kday)/180.))*BCSedRatio(kk)         !BCSedRatio(k) is multipler on BCTss that separates into different classes
			  Css(jj,2,kk)=BCTSS(jj)*BCSedRatio(kk)         !BUG wd not defined zw 04/07/2020

        !zw added 04/07/2020
			  Css(jj,1,kk)=Css(jj,2,kk)
			  Es(jj,1)=Es(jj,2)
			  Eh(jj,2)=Es(jj,2)
			  Eh(jj,1)=Eh(jj,2)
			  BCnosurge(jj,1)=BCnosurge(jj,2)
			  Qmarshmax(jj) = 0.0
        !zw added 04/15/2020 Water Temperature at offshore boundary cells
		      Tempw(jj,1)=TempwBC(jjk,kday)
			  Tempw(jj,2)=Tempw(jj,1)
	    enddo
	 enddo    !AMc 8 oct 2013
cc		endif


c****************************Central Difference -Hybrid Upwind for scalars**********************
c        fa=0.750  !1.0 upwind !0.5 central
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>> update upwind & downwind factors for dispersion based on previous timestep's upwind factor
      do i = 1,M
          fb(i)=(1.0-fa(i))				!weighting
      enddo

      fbc=0.705				!ch JAM July 13, 2009  was 0.75
	akL = 0.01				!wind induce currents
	akns=0.01
c      beginning of cell loop (flow, SS, Salinity, chem)

	do ichem = 1,14
		QChemSUM(ichem) = 0.0
	enddo

	!zw added 04/07/2020 to determine whether a cell is offbc or not
	flag_offbc(:)=0
	do jjk=1,mds
			jj=KBC(jjk)
			flag_offbc(jj)=1
	enddo


!============================cell updating
	do j=1,N
!============================cell continuity
! day, dday, kthr, and kday now global parameters - no longer needed to be passed into subroutines
      if(flag_offbc(j)==0) then !zw added 4/07/2020 for only non-offbc cells

          call CelldQ(j,kday,fcrop,mm,dday)

          if (iSed == 1) then
              call waves_YV(j)                                                            !-EDW
              call CelldSS(j,kday,kthr,CSSTRIBj,dref,Tres)               !JAM Oct 2010!
          endif

          if (iSal == 1) then
              call CelldSal(QSalSUM,j,kday,kthr,SalTRIBj,dref,Tres)
          endif

          if (iTemp == 1) then
              call CelldTmp(QTmpSUM,j,kday,kthr,SalTRIBj,dref,Tres)
          endif

!          call Celldage(QageSUM,Dz,j,kday,kthr,ageTRIBj,dref,Tres)

!          call CelldQ(QSUM,Dz,j,fcrop,mm)
!          call waves_YV(j)                                                            !-EDW
! QSSUM is now global parameter - doesn't need to be passed into subroutine
! QSsumh is calculated in CelldSS - therefore it is no longer passed into subroutine but is calculated only locally
!          call CelldQ(QSUM,Dz,j,fcrop,mm)
!          call waves_YV(j)

!          call CelldSS(Dz,j,CSSTRIBj,dref,Tres)               !JAM Oct 2010
!	    call CelldSal(QSalSUM,Dz,j,SalTRIBj,dref,Tres)
!          call CelldTmp(QSalSUM,Dz,j,SalTRIBj,dref,Tres)


!>> reset photosynthesis parameters so they can be reset for current compartment at timestep in photo.f
          muph = 0.0
          fpp = 0.0

!>> call photosnythesis rate calculation to pass into water quality routines that require it, this updates muph and fpp values
          if (iWQ == 1) then
              call photo(j)

!>> loop through WQ constituents for each WQ type that utilizes transport equation
              do ichem = 1,2
                  call CelldChem(j,kday,ichem)
		    enddo
!              write(*,*) chem(j,kday,1),chem(j,kday,2)
              ichem = 5
              call CelldChem(j,kday,ichem)

              do ichem = 7,11
                  call CelldChem(j,kday,ichem)
		    enddo

              !ichem = 14
              !call CelldChem(j,kday,ichem)  !zw 4/28/2015 POP calculated in MP2012 by ALG and DET see below

!>> Update WQ concentrations that are calculated from other constituents and did not rely on the transport equation
!>>-- set carbon-to-chlorophyll A ratio
	        rca = 75.0
!>>-- set nitrogen-chlorophyll A stoichiometric mass ratio for Organic N calculation
	        rna = 0.176*rca
!>>-- set nitrogen-detritus stoichiometric mass ratio for Organic N calculation
              rnd = 0.072
!>>-- set carbron-detritus stoichiometric mass ratio for Organic C calculation
              rcd = 0.41
!>>-- set carbon-nitrogen stoichiometric mass ratio for Organic C calculation
              rcn = 5.69

!>>-- Dissolved Inorganic N is function of NO3 and NH4
              Chem(j,3,2) = Chem(j,1,2) + Chem(j,2,2)

!>>-- Organic N is function of: dissovled organic N (ichem = 10) Algae (ichem = 8) and Detritus (ichem = 9) and stoichiometric ratios
              Chem(j,4,2) = Chem(j,10,2)+rna*Chem(j,8,2)+rnd*Chem(j,9,2)

!>>-- Total Organic Carbon is function of:dissovled organic N (ichem = 10) Algae (ichem = 8) and Detritus (ichem = 9) and stoichiometric ratios
              Chem(j,6,2) = rcn*Chem(j,10,2)+
     &                        rca*Chem(j,8,2)+rcd*Chem(j,9,2)

!>>-- Dissolved Inorganic P is function of Total Inorganic P (ichem = 5) and fraction of TIP that is in particulate form (calculated in photo.f)
              Chem(j,12,2) = (1-fpp)*chem(j,5,2)

!>>-- Chlorophyll A not calculated - set to -9999
              Chem(j,13,2) = -9999

! used to be located at end of hydrod
		    Chem(j,7,2) = O2Sat(int(day)+1)		!zw 4/28/2015 change Chem(j,7,1) to Chem(j,7,2)

! Particulate Organic P (POP) is function of ALG and DET
			Chem(j,14,2) = 0.0244*rca*Chem(j,8,2)+0.01*Chem(j,9,2)  !zw 4/28/2015 added POP calculation from MP2012

          endif
      endif

!>> Apply low and high pass filters to newly calculated WQ concentrations
	    !zw 4/28/2015 rove the WQ low and high pass filters
          !do ichem = 1,14
          !    if(Chem(j,ichem,2) < 0.0) Chem(j,ichem,2)=0.0
	    !enddo

          !if(Chem(j,8,2) <= 0.0001) Chem(j,8,2)=0.0001
	    !if(Chem(j,5,2) > 0.001) Chem(j,5,2)=0.001					! JKS Sept 2011
	    !if(Chem(j,12,2) > 0.001) Chem(j,12,2)=0.001				! JKS Sept 2011
	    !if(Chem(j,10,2) > 0.0005) Chem(j,10,2)=0.0005

          !if(Chem(j,1,2) > 0.003) Chem(j,1,2)=0.003					! JAM June 2011
	    !if(Chem(j,2,2) > 0.001) Chem(j,2,2)=0.001
	    !if(Chem(j,3,2) > 10.) chem(j,3,2)=10.0
	    !if(Chem(j,4,2) > 10.) chem(j,4,2)=10.0
	    !if(Chem(j,5,2) > 10.) chem(j,5,2)=10.0
	    !if(Chem(j,6,2) > 10.) chem(j,6,2)=10.0
	    !if(Chem(j,7,2) > 10.) chem(j,7,2)=10.0
	    !if(Chem(j,8,2) > 10.) chem(j,8,2)=10.0
	    !if(Chem(j,9,2) > 10.) chem(j,9,2)=10.0
	    !if(Chem(j,10,2) > 10.) chem(j,10,2)=10.0
	    !if(Chem(j,11,2) > 10.) chem(j,11,2)=10.0
	    !if(Chem(j,12,2) > 10.) chem(j,12,2)=10.0
!	    if(Chem(j,13,2) > 10.) chem(j,13,2)=10.0
	    !if(Chem(j,14,2) > 10.) chem(j,14,2)=10.0

      enddo										!j do loop



!>> Link updates of Q
      do i=1,M
!>> Reset upwind dispersion factor to default value - this will be updated to 1.0 if flows in each link are above threshold value
          fa(i) = fa_def*fa_mult(i)

!>> Link type < 0 = dummy link that is not active in current run
          if (linkt(i) <= 0) then
              Q(i,2) = 0.0

!>> Link type 1 = rectangular channels
!>> Link type 6 = rectangular culverts/bridges
!>> Link type 11 = marsh 'composite flow' channel
!>> Link type 12 = rectangular channel with elevations not updated by ICM/morph code
          elseif ((linkt(i) == 1) .or. (linkt(i) == 6) .or.
     &              (linkt(i) == 11) .or. (linkt(i) == 12)) then
              !! Channel attributes for channels (types 1,11,12)
              !! Latr1 = Channel invert elevation
              !! Latr2 = Channel bank elevation
              !! Latr3 = length
              !! Latr4 = width
              !! Latr5 = Roughness, n
              !! Latr6 = entrance loss, Ken
              !! Latr7 = exit loss, Kex
              !! Latr8 = Structure loss, Kstr
              !! Latr9 = n/a
              !! Latr10 = n/a

			!! Channel attributes for bridges/culverts (linkt=6)
              !! Latr1 = culvert invert elevation
              !! Latr2 = culvert crown elevation
              !! Latr3 = length
              !! Latr4 = width
              !! Latr5 = Roughness, n
              !! Latr6 = entrance loss, Ken
              !! Latr7 = exit loss, Kex
              !! Latr8 = Structure loss, Kstr
              !! Latr9 = n/a
              !! Latr10 = n/a
          !>> Set zero multiplier for infinitely high roughness values
              dkd=1.0

          !>> Determine flow direction
              if ( ES(jus(i),2) >= ES(jds(i),2) ) then
                  upN = jus(i)
                  downN = jds(i)
                  sn = 1.0
              else
                  upN = jds(i)
                  downN = jus(i)
                  sn = -1.0
              endif

          !>> Check for missing roughness attribute values and reassign to default values if missing
              if (Latr5(i) < 0.0) then
                  Latr5(i) = def_n
              elseif (Latr5(i) >= 1) then
          !>> Set zero multiplier for infinitely high roughness values
                  if (Latr5(i) >= 999.) then
                      dkd = 0.0
                  else
                      Latr5(i) = def_n
                  endif
              endif

          !>> If upstream water stage is lower than channel invert, flow is zero (downstream water stage is lower than upstream, so it is also lower than invert)
!              if (Es(upN,2) < latr1(i)) then
              if (Es(upN,2)<= latr1(i))then
                  Q(i,2) = 0.0
                  EAOL(i) = 0.0
          !>> If upstream water stage is lower than bed elevation, flow is zero (downstream water stage is lower than upstream, so it is also lower than bed)
              elseif( (Es(upN,2)-Bed(upN)) <=0.01)then	!zw 3/16/2015
                  Q(i,2) = 0.0
                  EAOL(i) = 0.0
          !>> If only downstream water stage is lower than channel invert, calculate average depth in channel with a zero depth at downstream end of channel
              else
                  !if(ES(downN,2) < latr1(i)) then
                  !    !linkdepth=max(0.5*(ES(upN,2)-latr1(i)),0.01)
                  !    linkdepth=0.5*(ES(upN,2)-latr1(i))	   !zw 3/16/2015

          !>> If both upstream and downstream water stages are above channel invert, calculate average depth from both upstream and downstream channel depths
                  !else
                  !    !linkdepth=max(0.5*(Es(upN,2)+ES(downN,2))-latr1(i),0.01)
                  !    linkdepth=0.5*(Es(upN,2)+ES(downN,2))-latr1(i)	 !zw 3/16/2015
                  !endif
				  
				  !zw added 5/5/2020 testing upwind approach
				  linkdepth=Es(upN,2)-latr1(i)
          !>> Calculate hydraulic radius from culvert dimensions (if bridge/culvert)
                  if (linkt(i) == 6) then
                      Rh=Latr4(i)*(Latr2(i)-Latr1(i))/
     &                 (2.*(Latr2(i)-Latr1(i))+Latr4(i))
                  else
		!>> Calculate hydraulic radius from average flow depth (if open channel)
		      Rh=linkdepth*Latr4(i)/(linkdepth*2.+Latr4(i))
	          endif

                  Ach = linkdepth*Latr4(i)

          !>> Calculate channel 'resistance'
		        AK = Latr8(i)+Latr6(i)+Latr7(i) ! Minor loss coefficients
                  ruf=Latr3(i)*Latr5(i)*Latr5(i)

                  Res=(AK/2./g+ruf/Rh**(4./3.))/(Ach*Ach)	!Resistance in links
                  Resist=1.0/sqrt(abs(Res))		!Resistance in links
                  if(isNaN(Resist)) then
                      Resist = 0.0
                  endif
                  Deta=Es(jus(i),2)-Es(jds(i),2)

          !>> Determine direction of flow
                  !sn=1.0		 !zw 3/14/2015 comment: sn already determine above
                  !if(Deta <= 0) then
                  !    sn=-1.0
                  !endif

                  Q(i,2)=sqrt(abs(Deta))*Resist*sn*dkd

                  Q_filter1 = abs(Deta)*As(downN,1)/dt                 ! yw flow filter
                  Q_filter2 = sqrt(g*linkdepth)*Latr4(i)*abs(Deta)     ! yw flow filter

                  flag_skip = 0
                  if (nlinkskip > 0) then
                      do jj = 1,nlinkskip
                          if (linkskip(jj) == i) then
                              flag_skip = 1
                          endif
                      enddo
                  endif

                  if (flag_skip == 0) then
                      Q(i,2) = sn*min(abs(Q(i,2)),Q_filter1,Q_filter2)     ! yw flow filter
                  endif

                  if(isNan(Q(i,2))) then
                      write (*,*) 'Link',i,'flow is NaN'
                      write(*,*) 'Deta=', Deta
                      write(*,*) 'Res=',Res
                      write(*,*) 'Resist=',Resist
                      write(*,*) 'dkd=',dkd
                      write(*,*) 'ruf=',ruf
                      write(*,*) 'AK=',AK
                      write(*,*) 'Ach=',Ach
                      write(*,*) 'linkdepth=',linkdepth
                      write(*,*) 'invert=',latr1(i)
                      write(*,*) 'Es(jus)=',Es(jus(i),2)
                      write(*,*) 'Es(jds)=',Es(jds(i),2)
                      pause
                  endif
                  EAOL(i)=Exy(i)*Ach/Latr3(i)*dkd  !zw 3/14/2015 add *dkd for high roughness no flow case
              endif
          !>> update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
              if (abs(Q(i,2)/Ach) >= upwind_vel) then
                  fa(i) = 1.0
              endif

          !>> Link type 2 = weirs
          elseif( linkt(i) == 2) then
              !! Channel attributes for weirs
              !! Latr1 = weir crest elevation
              !! Latr2 = invert elevation upstream of weir
              !! Latr3 = invert elevation downstream of weir
              !! Latr4 = weir width (crest width - perpendicular to flow)
              !! Latr5 = n/a
              !! Latr6 = n/a
              !! Latr7 = n/a
              !! Latr8 = weir coefficient
              !! Latr9 = n/a
              !! Latr10 = n/a

          !>> set weir coefficient to default value if outside of standard range
              if (Latr8(i) < 0.55) then
                  Latr8(i) = 0.62
              elseif (Latr8(i) > 0.65) then
                  Latr8(i) = 0.65
              endif
           !>> check for missing invert attributes and assign upstream invert to downstream value
              if (Latr2(i) < -9990.) then
                  Latr2(i) = Latr3(i)
              elseif (Latr2(i) > 100.) then
                  Latr2(i) = Latr3(i)
              endif

            !>> determine direction of flow through link and re-assign upstream/downstream attributes if necessary

              if ( ES(jus(i),2) >= ES(jds(i),2) ) then
                  upN = jus(i)
                  downN = jds(i)
                  sn = 1.0
                  invup = Latr2(i)
                  invdown = Latr3(i)
              else
                  upN = jds(i)
                  downN = jus(i)
                  sn = -1.0
                  invup = Latr3(i)
                  invdown = Latr2(i)
              endif

              htup = ES(upN,2) - Latr1(i)
          !>> check if upstream water elevation is above weir crest - if not, set flow to 0
              !if (htup <= 0.0) then
              if ((htup <= 0.0).OR.((Es(upN,2)-Bed(upN))<=0.01)) then	   !zw 3/16/2015
                  Q(i,2) = 0.0
          !>> otherwise, if water is above weir crest calculate flow
              else
              !>> check if downstream water is below ridge elevation - if so, set downstream height to 0
                  if (ES(downN,2) > Latr1(i)) then
                      htdn = ES(downN,2) - Latr1(i)
                  else
                      htdn = 0.0
                  endif
                  p = Latr1(i) - invup

                  if(htup /= 0.0) then
                      subp = htdn/htup
                  else
                      subp = 1.0
                  endif


          !>> If weir is fully submerged, calculate flow with submerged weir equation
                  if (subp > 0.95) then
                      Q(i,2) = sn*0.6*Latr4(i)*htup*
     &                        sqrt(abs(2*g*(htup-htdn)))
                  else
          !>> If weir is partially submerged calculate submergence coefficient, ksub
                      if (subp >= 0.85) then
                          w_ksub = -14.167*subp**2.+23.567*subp-8.815
                      else
                          w_ksub = 1.0
                      endif
          !>> If not fully submerged calculate flow with standard homogeneous weir equation Ven Te Chow (1959)
                      w_k = w_ksub*(Latr8(i)+0.05*htup/p)*
     &                          max(0.6,(1-0.2*htup/p))
                      Q(i,2) = sn*w_k*Latr4(i)*sqrt(2*g)*htup**1.5
                  endif
              endif
              ! Assume flow length through weir is 10 m for EAOL calculations
              if (htup > 0.0) then
                  !if(Exy(i) > -9990.0) then
                  if(Exy(i) < 0) then	!zw 3/14/2015 revised
                      if (mm == 1) then

                          write(1,*)'Exy value for link',i,'is missing.'
                          write(1,*) 'Default value of 100 is assigned
     &    but this should be fixed in the input file.'

                          write(*,*)'Exy value for link',i,'is missing.'
                          write(*,*) 'Default value of 100 is assigned
     &    but this should be fixed in the input file.'

                      endif
                      Exy(i) = 100.0
                  endif
                  EAOL(i)=Exy(i)*htup*Latr4(i)/10.

          !>> update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
          !        if (abs(Q(i,2)/(htup*Latr4(i))) >= upwind_vel) then
          !            fa(i) = 1.0
          !        endif
          !    else
                  EAOL(i) = 0.0
              endif

!>> Link type 3 = rectangular channels with a lock control structure
          elseif (linkt(i) == 3) then

              !! Channel attributes for channels with locks
              !! Latr1 = Channel invert elevation
              !! Latr2 = second control scheme threshold that shuts lock
                      ! if Latr9 = 1, Latr2 = -9999
                      ! if Latr9 = 2, Latr2 = -9999
                      ! if Latr9 = 3, Latr2 = -9999
                      ! if Latr9 = 4, Latr2 = -9999
                      ! if Latr9 = 5, Latr2 = d/s salinty (ppt)
                      ! if Latr9 = 6, Latr2 = -9999
              !! Latr3 = channel length
              !! Latr4 = channel width
              !! Latr5 = channel Roughness, n
              !! Latr6 = channel entrance loss, Ken
              !! Latr7 = channel exit loss, Kex
              !! Latr8 = Structure loss through lock (when open), Kstr
              !! Latr9 = control scheme (1=diff stage, 2=hour, 3=d/s WSEL, 4=d/s Sal, 5=d/s WSEL or Sal,6=observed timeseries)
              !! Latr10 = control scheme threshold value that shuts lock
                      ! if Latr9 = 1, Latr10 = diff stage (m)
                      ! if Latr9 = 2, Latr10 = -9999
                      ! if Latr9 = 3, Latr10 = d/s WSEL (m)
                      ! if Latr9 = 4, Latr10 = d/s salinty (ppt)
                      ! if Latr9 = 5, Latr10 = d/s WSEL (m)
                      ! if Latr9 = 6, Latr10 = corresponding column in LockControlObservedData.csv
          !>> Initialize link's on/off flag to on
              dkd=1.0	  !zw 3/14/2015 moved to here from below

          !>> Check for missing roughness attribute values and reassign to default values if missing
              if (Latr5(i) < 0.0) then
                  Latr5(i) = def_n
              elseif (Latr5(i) >= 1) then
          !>> Set zero multiplier for infinitely high roughness values
                  if (Latr5(i) >= 999.) then
                      dkd = 0.0
                  else
                      Latr5(i) = def_n
                  endif
              endif

          !>> Set channel length to default value if missing
              if (Latr3(i) <= 0.0) then
                  if (mm == 1) then

                      write(1,*)'Length (Latr3) for link',i,'is missing'
                      write(1,*) 'Default value of 1000 is assigned
     & but this should be fixed in the input file.'

                      write(*,*)'Length (Latr3) for link',i,'is missing'
                      write(*,*) 'Default value of 1000 is assigned
     & but this should be fixed in the input file.'

                  endif
                  Latr3(i) = 1000.0
              endif

              if ( ES(jus(i),2) >= ES(jds(i),2) ) then
                  upN = jus(i)
                  downN = jds(i)
                  sn = 1.0
              else
                  upN = jds(i)
                  downN = jus(i)
                  sn = -1.0
              endif

          !>> Lock control rules

          !>> Set zero multiplier if differential head is greater than differential stage trigger
          !>> Lock is closed if differential stage is greater than trigger
          !>> This does not take into account which stage (u/s or d/s) is greater
              if (Latr9(i) == 1) then
                  if (abs(ES(upN,2)-ES(downN,2)) > Latr10(i)) then
                      dkd = 0.0
                  endif
          !>> Set zero multiplier if lock is closed due to time of day
              elseif (Latr9(i) == 2) then
                  if (hourclosed(i,simhr+1) == 1) then
                      dkd = 0.0
                  endif
          !>> Set zero multiplier if lock should be closed due to high downstream water level
              elseif (Latr9(i) == 3) then
                  if(Es(downN,2) > Latr10(i)) then
                      dkd = 0.0
                  endif
          !>> Set zero multiplier if lock should be closed due to high downstream salinity
              elseif (Latr9(i) == 4) then
                  if(S(downN,2) > Latr10(i)) then
                      dkd = 0.0
                  endif

          !>> Set zero multiplier if lock should be closed due to high downstream water level OR salinity
              elseif (Latr9(i) == 5) then
                  if(S(downN,2) > Latr2(i)) then
                      dkd = 0.0
                  elseif(Es(downN,2) > Latr10(i)) then
                      dkd = 0.0
                  endif

          !>> Set zero multiplier if lock should be closed based on observed daily timeseries
              elseif (Latr9(i) == 6) then
                  if (lockhours(lockrow,Latr10(i)) <= 0.0) then    ! lockrow is data timestep counter updated in main.f
                      dkd = 0.0
                  endif

              endif
           !>> If both upstream and downstream water stages are lower than channel invert, flow is zero
              !if (Es(upN,2) < latr1(i)) then
              if (Es(upN,2) <= latr1(i) ) then	!zw 3/16/2015
                  Q(i,2) = 0.0
                  EAOL(i) = 0.0
              elseif ( ( Es(upN,2)-Bed(upN) ) <=0.01) then
                  Q(i,2) = 0.0
                  EAOL(i) = 0.0
          !>> If only downstream water stage is lower than channel invert, calculate average depth in channel with a zero depth at downstream end of channel
              else
                  if(ES(downN,2) < latr1(i)) then
                      !linkdepth=max(0.5*(ES(upN,2)-latr1(i)),0.01)
                      linkdepth=0.5*(ES(upN,2)-latr1(i))	   !zw 3/16/2015

          !>> If both upstream and downstream water stages are above channel invert, calculate average depth from both upstream and downstream channel depths
                  else
                      !linkdepth=max(0.5*(Es(upN,2)+ES(downN,2))-latr1(i),0.01)
                      linkdepth=0.5*(Es(upN,2)+ES(downN,2))-latr1(i)	  !zw 3/16/2015
                  endif


			    Rh=linkdepth*Latr4(i)/(linkdepth*2.0+Latr4(i))

                  Ach = linkdepth*Latr4(i)


          !>> Calculate channel 'resistance'
		        AK = Latr8(i)+Latr6(i)+Latr7(i) ! Minor loss coefficients
                  ruf=Latr3(i)*Latr5(i)*Latr5(i)

                  Res=(AK/2./g+ruf/Rh**(4./3.))/(Ach*Ach)	!Resistance in links
                  Resist=1.0/sqrt(abs(Res))		!Resistance in links
                  Deta=Es(upN,2)-Es(downN,2)

          !>> Put a low pass filter on Deta to avoid infinity and NaN values
                  !Deta = min(abs(Deta),100.)	  !zw 3/14/2015 comments out:Es is already checked in celldQ

                  Q(i,2)=sqrt(abs(Deta))*Resist*sn*dkd

                  EAOL(i)=Exy(i)*Ach/Latr3(i)*dkd  !zw 3/14/2015 add *dkd for no flow condition under lock control rules
              endif

          !>> update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
              if (abs(Q(i,2)/Ach) >= upwind_vel) then
                  fa(i) = 1.0
              endif


!>> Link type 4 = tide gate, orifice flow in only one direction
!>> Link type 5 = orifice flow in two directions
          elseif ((linkt(i) == 4) .or. (linkt(i) == 5)) then

              !! Channel attributes for orifices
              !! Latr1 = Invert elevation of orifice opening
              !! Latr2 = Crown elevation of orifice opening
              !! Latr3 = Ground elevation upstream
              !! Latr4 = average width of opening
              !! Latr5 = Ground elevation downstream
              !! Latr6 = n/a
              !! Latr7 = n/a
              !! Latr8 = orifice coefficient
              !! Latr9 = n/a
              !! Latr10 = n/a

			if ( ES(jus(i),2) >= ES(jds(i),2) ) then
                  upN = jus(i)
                  downN = jds(i)
                  invup = Latr3(i)
				sn = 1.0
              else
                  upN = jds(i)
                  downN = jus(i)
                  invup = Latr5(i)

	!>> if link type is orifice, set flow sign to negative if downstream water surface is greater than upstream
				if(linkt(i) == 5) then
					sn = -1.0
				else
	!>> if link type is flood gate and downstream water elevation is greater, set flow sign to 0, so backwards flow is set equal to zero
					sn = 0.0
				endif
			endif
			delp = ES(upN,2) - ES(downN,2)
			htup = ES(upN,2) - Latr1(i)
              htdn = ES(downN,2) - Latr1(i)

	!>> if water level is below orifice invert, there is no flow
			!if (htup <= 0.0) then
			if ((htup <= 0.0).OR.((Es(upN,2)-Bed(upN))<=0.01)) then	  !zw 3/16/2015
				Q(i,2) = 0.0
                  ! set flow area for no-flow orifice
                  Ach = 0.0
	!>> If water water level is above tide gate invert, but below crown, flow is calculated with weir equation
			elseif (Es(upN,2) <= Latr2(i)) then
	            if(htup /= 0.0) then
                      subp = htdn/htup
                  else
                      subp = 1.0
                  endif
                  p = Latr1(i) - invup
                  eqC = 0.62
		!>> If orifice-as-equivalent-weir is fully submerged, calculate flow with submerged weir equation
				if (subp > 0.95) then
					Q(i,2) = sn*0.6*Latr4(i)*htup*sqrt(2*g*(htup-htdn))
				else
          !>> If orifice-as-equivalent-weir is partially submerged calculate submergence coefficient, ksub
					if (subp >= 0.85) then
						w_ksub = -14.167*subp**2.+23.567*subp-8.815
					elseif (subp < 0.85) then
						w_ksub = 1.0
					endif
          !>> If orifice-as-equivalent-weir is not fully submerged calculate flow with standard homogeneous weir equation Ven Te Chow (1959)
					w_k = w_ksub*(eqC+0.05*htup/p)*max(0.6,(1-0.2*htup/p))
					Q(i,2) = sn*w_k*Latr4(i)*sqrt(2*g)*htup**1.5
				endif
                  ! set flow area for non-full orifice
				Ach = Latr4(i)*htup
	!>> If orifice crown is submerged, use orifice equation
			else
				orarea = (Latr2(i)-Latr1(i))*Latr4(i)
				beta = (Latr2(i)-Latr1(i))/(Es(upN,2)-invup)
				orC = Latr8(i)/(1-beta**4.0)

				Q(i,2) = sn*orC*orarea*sqrt(abs(2*delp*g))

                  ! set flow area for full orifice
                  Ach = orarea
			endif
! update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
!             if (abs(Q(i,2)/Ach) >= upwind_vel) then
!                  fa(i) = 1.0
!              endif
      ! Assume flow length through weir is 10 m for EAOL calculations
              EAOL(i)=abs(sn)*Exy(i)*Ach/10. ! multiply EAOL here by abs(sn) to zero out EAOL for backwards flow in tide gates (where sn=0)

!>> Link type 7 = pumps
          elseif (linkt(i) == 7) then

              !! Channel attributes for over pumps
              !! Latr1 = upstream water elevation stage which triggers pump to turn on
              !! Latr2 = upstream water elevation stage which triggers pump to turn off
              !! Latr3 = n/a
              !! Latr4 = n/a
              !! Latr5 = n/a
              !! Latr6 = n/a
              !! Latr7 = n/a
              !! Latr8 = n/a
              !! Latr9 = pump capacity (m**3/s)
              !! Latr10 = n/a

          !>> Set pump to off by default - will turn on if certain conditions are met
              pumpon = 0
          !>> If pump was off during previous timestep, check if water is now above turn-on level, this will keep pump off if water is above shut-off level but below turn-on level
              if (Q(i,1) == 0.0) then
                  if (Es(jus(i),2) > Latr1(i)) then
                      pumpon = 1
                  endif
          !>> If pump was on during previous timestep, check if water is still above shut-off level
              else
                  if (Es(jus(i),2) > Latr2(i)) then
                      pumpon = 1
                  endif
              endif

              if (pumpon == 0) then
                  Q(i,2) = 0.0
              else
                  cden=1/1000./24/3600.		!JAM Oct 2010 mm/d to m/s
                  infilrate= 2.54              ! urban infiltration rate (mm/hr)
                  AET = PET(kday,Jet(jus(i)))
                  Pexcess =(Rain(kday,jrain(jus(i)))-
     &                         AET-infilrate*24.0)*cden
                  Qrunoff = Ahydro(jus(i))*max(0.0,Pexcess)
                  Qpump = (Latr9(i) + Qrunoff)/2.
                  Q(i,2) = min(Latr9(i),Qpump)
              endif

              ! EAOL term is inactive for pumps due to no flow area/flow length properties
              EAOL(i) = 0.0

!>> Link type 8 = marsh over land channels
          elseif (linkt(i) == 8) then

              !! Channel attributes for over land marsh channels
              !! Latr1 = Marsh elevation along channel
              !! Latr2 = Marsh elevation in upstream compartment
              !! Latr3 = length
              !! Latr4 = width
              !! Latr5 = Roughness, Manning's n
              !! Latr6 = n/a
              !! Latr7 = n/a
              !! Latr8 = n/a
              !! Latr9 = n/a
              !! Latr10 = Marsh elevation in downstream compartment

          !>> Set zero multiplier for infinitely high roughness values
              dkd=1.0
          !>> Check for missing roughness attribute values and reassign to default values if missing
              if (Latr5(i) < 0.0) then
                  Latr5(i) = def_n
              elseif (Latr5(i) >= 1) then
          !>> Set zero multiplier for infinitely high roughness values
                  if (Latr5(i) >= 999.) then
                      dkd = 0.0
                  else
                      Latr5(i) = def_n
                  endif
              endif

          !>> determine direction of flow through link and re-assign upstream/downstream attributes if necessary
              if ( EH(jus(i),2) >= EH(jds(i),2) ) then
                  upN = jus(i)
                  downN = jds(i)
                  sn = 1.0
              else
                  upN = jds(i)
                  downN = jus(i)
                  sn = -1.0
              endif

          !>> if upstream water level is at or below marsh elevation, flow is equal to zero
              !if (EH(upN,2) <= Latr1(i)) then
              if (EH(upN,2) <= Latr1(i))then
                  Q(i,2) = 0.0
                  Ach = 0.0
              elseif ((EH(upN,2)-BedM(upN))<=0.1) then	!zw 3/16/2015
                  Q(i,2) = 0.0
                  Ach = 0.0
          !>> calculate flow if water elevation is above marsh invert
              else
          !>> if both downstream and upstream water levels are above marsh, calculate average depth
                  if (EH(downN,2) > Latr1(i)) then
                      !avdep  =max(0.5*(EH(upN,2)+EH(downN,2))-Latr1(i),0.01)
                      avdep  =0.5*(EH(upN,2)+EH(downN,2))-Latr1(i)  !zw 3/16/2015
                      delh = abs(EH(upN,2)-EH(downN,2))
          !>> if downstream water elevation is below marsh, calculate depth from marsh elevation
                  else
                     !delh = max(EH(upN,2)-Latr1(i),0.01)
                     delh = EH(upN,2)-Latr1(i)
                     avdep = .5*delh
                  endif

                  Q(i,2) = sn*dkd*avdep**(5./3.)*(Latr4(i)/Latr5(i))*
     &                            (delh/Latr3(i))**(1./2.)

              ! set flow area for EAOL calculation
                  Ach = avdep*latr4(i)

              endif
          ! update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
          !    if (abs(Q(i,2)/Ach) >= upwind_vel) then
          !        fa(i) = 1.0
          !    endif

              EAOL(i)=Exy(i)*Ach/Latr3(i)*dkd  !zw 3/14/2015 add *dkd for high roughness no flow conditions

!>> Link type 9 = ridges/levees
          elseif( linkt(i) == 9) then
              !! Channel attributes for ridges/levees
              !! Latr1 = ridge crest elevation
              !! Latr2 = invert elevation upstream of ridge
              !! Latr3 = length of ridge crest (flowlength across top of ridge - parallel to flow)
              !! Latr4 = width of ridge crest (flowpath width - perpendicular to flow)
              !! Latr5 = Roughness, Manning's n
              !! Latr6 = entrance loss, Ken
              !! Latr7 = exit loss, Kex
              !! Latr8 = Structure loss, Kstr
              !! Latr9 = Initial Flowrate, Q
              !! Latr10 = invert elevation downstream of ridge


!>> determine direction of flow through link and re-assign upstream/downstream attributes if necessary
              if ( ES(jus(i),2) >= ES(jds(i),2) ) then
                   upN = jus(i)
                   downN = jds(i)
                   sn = 1.0
                   invup = Latr2(i)
                   invdown = Latr10(i)
              else
                   upN = jds(i)
                   downN = jus(i)
                   sn = -1.0
                   invup = Latr10(i)
                   invdown = Latr2(i)
              endif
!>> check for missing invert elevations
              if (invup < -9990.0) then          ! invup is missing
                  if (invdown > -9990.0) then     ! invdown exists
                      invup = invdown + 0.1
                  else                            ! both are missing
                      invup = -3.0
                      invdown = -3.1
                  endif
              elseif (invdown < -9990.0) then    ! invup exists and invdown is missing
                  invdown = invup - 0.1
              endif

              htup = ES(upN,2) - Latr1(i)

!>> check if upstream water elevation is above ridge elevation - set flow to zero if it is lower than ridge
              if (htup <= 0.0) then
                  Q(i,2) = 0.0
				Ach = 0.0    !zw 3/14/2015 add to ensure diffusion is zero for no flow condtion
!>> otherwise, calculate flow across ridge
              else
!>> check if downstream water is below ridge elevation - if so, set downstream height to 0
                  if (ES(downN,2) > Latr1(i)) then
                      htdn = ES(downN,2) - Latr1(i)
                  else
                      htdn = 0.0
                  endif

                  if(htup /= 0.0) then
                      subp = htdn/htup
                  else
                      subp = 1.0
                  endif

                  p = Latr1(i) - invup

           !>> calculate average flow depth across ridge crest
                  dv = max(0.5*(htup+htdn),0.01)
           !>> calculate headloss due to friction as flow passes over ridge crest (function of Q@t-1)
                  hf=Latr3(i)*(Latr5(i)*Q(i,1)/
     &                          (Latr4(i)*dv**(5./3.)))**2.0
           !>> correct upstream water elevation for headloss (only used in flow calculations - NOT in submerged coefficient calculations)
           !>> if headloss is too great, set upstream water elevation to value JUST above downstream water elevation
                  htupf = max((htup-hf),(htdn+0.01))

               ! set flow area for use in EAOL calculation
                  Ach = dv*Latr4(i)

           !>> If weir is fully submerged, calculate flow with submerged weir equation (with headloss corrected upstream water elevation)
                  if (subp > 0.95) then
                      Q(i,2)=sn*0.6*Latr4(i)*htupf
     &                        *sqrt(abs(2*g*(htupf-htdn)))
                  else
           !>> If weir is partially submerged calculate submergence coefficient
                      if (subp >= 0.85) then
                           w_ksub = -14.167*subp**2.+23.567*subp-8.815
                      elseif (subp < 0.85) then
                           w_ksub = 1.0
                      endif
           !>> If not fully submerged calculate flow with standard homogeneous weir equation Ven Te Chow (1959) (with headloss corrected upstream water elevation)
                      w_k=w_ksub*(Latr8(i)+0.05*htup/p)*
     &                          max(0.6,(1-0.2*htup/p))
                      Q(i,2) = sn*w_k*Latr4(i)*sqrt(2*g)*abs(htupf)**1.5

                  endif
              endif
          ! update upwind factor for salinity/WQ dispersion if channel velocity is greater than threshold value
          !    if (abs(Q(i,2)/Ach) >= upwind_vel) then
          !        fa(i) = 1.0
          !    endif

              EAOL(i)=Exy(i)*Ach/Latr3(i)

!>> Link type 10 = regime channels for sediment distribution
          elseif (linkt(i) == 10) then

              !! Channel attributes for regime channels
              !! Latr1 = not used
              !! Latr2 = corresponding channel link number, M
              !! Latr3 = length of corresponding channel link
              !! Latr4 = width
              !! Latr5 = Roughness, n = 0.025
              !! Latr6 = Entrance loss, ken
              !! Latr7 = Exit loss, kex
              !! Latr8 = Structure loss, kstr
              !! Latr9 = Regime flowrate (cms)
              !! Latr10 = Regime D50 (mm)

              dkd=1.0
          !>> Check for missing roughness attribute values and reassign to default values if missing
              if (Latr5(i) < 0.0) then
                  Latr5(i) = 0.025 !0.025 is default n for regime channels
              elseif (Latr5(i) >= 1) then
          !>> Set zero multiplier for infinitely high roughness values
                  if (Latr5(i) >= 999.) then
                      dkd = 0.0
                  else
                      Latr5(i) = 0.025 !0.025 is default n for regime channels
                  endif
              endif
           !>> calculate flow from Lacey Regime channel equations
              reg_p = 4.84*Latr9(i)**(1./2.)
              reg_fs = 1.587*Latr10(i)**(1./2.)
              reg_r = 0.468*(Latr9(i)/reg_fs)**(1./3.)
              reg_s = 0.000316*reg_fs**(5./3.)/Latr9(i)**(1./6.)
              sn = Q(Latr2(i),2)/abs(Q(Latr2(i),2))
              Q(i,2) = sn*reg_p*reg_r**(5./3.)*reg_s**(1./2.)/Latr5(i)

!>> Check to see if regime channel link flow is greater than original channel flow
!>> If so, set original channel to zero and use regime channel flow
!>> If not, set regime channel flow to zero and use original channel flow
              if(abs(Q(i,2)) >= abs(Q(Latr2(i),2))) then
                  Q(Latr2(i),2) = 0.0
              else
                  Q(i,2) = 0.0
              endif
              EAOL(i)=EAOL(Latr2(i))

!>> end flowrate calculations based on linktype
          endif



!>> Calculate volume of water in upstream compartment available for exchange during timestep
          sn =(Es(jus(i),2)-Es(jds(i),2))/abs(Es(jus(i),2)-Es(jds(i),2))
          if (sn > 0) then
!              volavailable = abs(Es(jus(i),2)-Es(jds(i),2))*As(jus(i),1)
              volavailable = Max(Es(jus(i),2)-Bed(jus(i)),0.01)*
     &                        As(jus(i),1)

          else
!              volavailable = abs(Es(jus(i),2)-Es(jds(i),2))*As(jds(i),1)
              volavailable = Max(Es(jds(i),2)-Bed(jds(i)),0.01)*
     &                        As(jds(i),1)
          endif
!>> Calculate separate volume of water available for exchange for overland marsh flow links
          if (linkt(i) == 8) then
              if (Eh(jus(i),2) > Eh(jds(i),2)) then
                  volavailable=Max(Eh(jus(i),2)-Bedm(jus(i)),0.01)*
     &                            Ahf(jus(i))
              else
                  volavailable=Max(Eh(jds(i),2)-Bedm(jds(i)),0.01)*
     &                            Ahf(jds(i))
              endif
          endif

!>> Determine flowrate required to exchange all available water during timestep
          Qmax = volavailable/dt
!>> If calculated flowrate is greater than needed to exchange all available water, set flowrate to this max rate
          if (linkt(i) >= 0.0) then
              if (Qmax < abs(Q(i,2))) then
                  Q(i,2) = Qmax*sn
                  if (Qmax /= 0.0) then
                      if (daystep == lastdaystep) then
                          Write(1,*) 'Max flowrate reached. Link:',i
                          Write(*,*) 'Max flowrate reached. Link:',i
                      endif
                  endif
              endif
          endif

!>> End link updates of Q
      enddo

!>> ****HARDCODED FLOW CONTROL TRIGGERS - START*****
!>> MP project - CSC - set flow to zero for Calcesieu Salinity Control measures if salinity threshold @ Calcesieu Lake is met
!! Turner's Bay and East Pass structures are links 3894 and 3891 - MP project # 004.HR.06
!      if (S(860,2) > 12 )  then
!          Q(3884,2) = 0.0
!          Q(3891,2) = 0.0
!          Q(4023,2) = 0.0
!          Q(4030,2) = 0.0
!
!          SalLockTriggerCSC = -1
!      else
!          SalLockTriggerCSC = 1
!      endif
!>> print statements to indicate lock opening/closing
!      if (SalLockStatusCSC /= SalLockTriggerCSC) then
!          if (SalLockStatusCSC == 1) then
!               write(*,*) ' Calcesieu Lake salinity > 12:'
!               write(*,*)'  - Control structures closed at timestep ',mm
!               write(1,*) ' Calcesieu Lake salinity > 12:'
!               write(1,*)'  - Control structures closed at timestep ',mm
!          elseif (SalLockStatusCSC == -1) then
!               write(*,*) ' Calcesieu Lake salinity > 12:'
!               write(*,*)'  - Control structures opened at timestep ',mm
!               write(1,*) ' Calcesieu Lake salinity > 12:'
!               write(1,*)'  - Control structures opened at timestep ',mm
!          endif
!          SalLockStatusCSC = -1*SalLockStatusCSC
!      endif



!!>> MP project - Increase Atch Flow to Terrebonne - set flowrates for Atchafalaya diversion projects - MP project # 03b.DI.04
!!      link numbers to use
!      Atch_US_link = 1272
!!      Div_link = 3859
!      Div_link = 4023
!!      Div_link = 4035
!!      flowrates to use
!      Atch_US_Q = Q(Atch_US_link,2)
!
!      if (linkt(Div_link) > 0) then
!
!            ! if diversion is open, check stage on either side of lock to determine if diversion should be closed
!          if (Atch_div_onoff == 1) then
!              if (Es(509,2) > 0.9144) then
!                  Atch_div_onoff = 0
!                  write(*,*) ' Stage EAST of Bayou Boeuf lock > 3ft'
!                  write(*,*) ' Atch diversion closed at timestep',mm
!                  write(1,*) ' Stage EAST of Bayou Boeuf lock > 3 ft'
!                  write(1,*) ' Atch diversion closed at timestep',mm
!              endif
!          ! if diversion is closed, check stage on either side of lock to determine if diversion should be re-opened
!          else
!              if (Es(509,2) <= 0.6096) then
!                  Atch_div_onoff = 1
!                  write(*,*)' Stage EAST of Bayou Boeuf lock <= 2ft'
!                  write(*,*)' Atch. diversion opened at timestep',mm
!                  write(1,*)' Stage EAST of Bayou Boeuf lock <= 2ft'
!                  write(1,*)' Atch. diversion opened at timestep',mm
!              endif
!          endif
!! flow in diversion is 11% of river if operated separately
!!          Div_Q = 0.11*Atch_US_Q*Atch_div_onoff
!! flow in diversion is 8% of river if operated jointly with 03a.DI.05
!          Div_Q = 0.08*Atch_US_Q*Atch_div_onoff
!      else
!          Div_Q = 0.0
!      endif
!
!
!
!!>> MP project - Atchafalaya River Diversion - set flowrates for Atchafalaya diversion projects - MP project # 03a.DI.05
!      !! link numbers to use
!      Atch_US_link = 1274
!      !Atch_DS_link = 1275
!!      Div_link = 3859
!      Div_link = 4031
!
!      !! flowrates to use
!      Atch_US_Q = Q(Atch_US_link,2)
!
!      if (linkt(Div_link) > 0) then
!! flow in diversion is 26% of river if operated separately
!!          Div_Q = 0.26*Atch_US_Q
!! flow in diversion is 17% of river if operated jointly with 03a.DI.04
!          Div_Q = 0.17*Atch_US_Q
!      else
!          Div_Q = 0.0
!      endif
!
!!! update diversion link with diversion flowrate
!      volavailable = Max(Es(jus(Div_link),2)-Bed(jus(Div_link)),0.01)*
!     &                        As(jus(Div_link),1)
!      Qmax = volavailable/dt
!
!      Q(Div_link,2) = min(Div_Q,Qmax)
!
!
!

!>> FWOA project - HNC - set flow to zero for Bubba Dove floodgate on the Houma Navigation Canal if salinity threshold @ LUMCON is met
! Bubba Dove flood gate is link #1351, LUMCON is located in compartment #494
!      if (S(494,2) > 7.5 )  then
!          Q(1351,2) = 0.0
!          SalLockTriggerHNC = -1
!      else
!          SalLockTriggerHNC = 1
!      endif
!>> print statements to indicate lock opening/closing
!      if (SalLockStatusHNC /= SalLockTriggerHNC) then
!          if (SalLockStatusHNC == 1) then
!               write(*,*) ' LUCMON salinity > 7.5:'
!               write(*,*) '  - HNC @ Bubba Dove closed at timestep ',mm
!               write(1,*) ' LUCMON salinity > 7.5:'
!               write(1,*) '  - HNC @ Bubba Dove closed at timestep ',mm
!          elseif (SalLockStatusHNC == -1) then
!               write(*,*) ' LUCMON salinity < 7.5:'
!               write(*,*) '  - HNC @ Bubba Dove opened at timestep ',mm
!               write(1,*) ' LUCMON salinity < 7.5:'
!               write(1,*) '  - HNC @ Bubba Dove opened at timestep ',mm
!          endif
!          SalLockStatusHNC = -1*SalLockStatusHNC
!      endif

!>> FWOA project - Superior Canal - set flow to zero for floodgates that are controlled by stage in Superior Canal @ the Hwy 82 Gauge
!      if (max(ES(jus(2807),2),ES(jds(2807),2)) <= 0.229) then
!          Q(2870,2) = 0.0
!          Q(3006,2) = 0.0
!          Q(3844,2) = 0.0
!          Q(3845,2) = 0.0
!          StgTriggerSuperiorCanal = -1
!      else
!          StgTriggerSuperiorCanal = 1
!      endif
!>> print statements to indicate tidegate opening/closing
!      if (StgTriggerStatusSuperiorCanal /= StgTriggerSuperiorCanal) then
!          if (StgTriggerStatusSuperiorCanal == 1) then
!               write(*,*) ' Water stage in Superior Canal <= 0.229 m'
!               write(*,*) '  - tidegates closed at timestep =',mm
!               write(1,*) ' Water stage in Superior Canal <= 0.229 m'
!               write(1,*) '  - tidegates closed at timestep =',mm
!          elseif (StgTriggerStatusSuperiorCanal == -1) then
!               write(*,*) ' Water stage in Superior Canal > 0.229 m'
!               write(*,*) '  - tidegates opened at timestep =',mm
!               write(1,*) ' Water stage in Superior Canal > 0.229 m'
!               write(1,*) '  - tidegates opened at timestep =',mm
!          endif
!          StgTriggerStatusSuperiorCanal=-1*StgTriggerStatusSuperiorCanal
!      endif


!>> ****HARDCODED FLOW TRIGGERS - END*****




!>> Kadlec-Knight flow into/out of marsh.
	do j=1,N
!		dkd=1.0
!		drefm=0.1										!Manning's n based on reference depth of 1/3 foot.
!		if(han(j) >= 999.) then
!              dkd=0.0
!         endif
!		Rh=max(Eh(j,1)-BedM(j),0.01)					!hydraulic radius for exchange flow
!		Ach = Rh*hWidth(j)+100.*RSSS					!exchange X-sectional area & 1% slope of topography
!		AK=hkm(j)										!minor loss

          !>> sign of flow direction (must match celldQ), positive is open water-to-marsh flow, negative is marsh-to-open water flow
          if (Es(j,2) > Eh(j,2)) then
              sDetah = 1.0
          else
              sDetah = -1.0
          endif


          !>> Stage differential between open water and marsh (with sign, negative is marsh-to-open water flow direction)
		Detah=sDetah*max(0.0,abs(Es(j,2)-Eh(j,2)))							!driving head into marsh

          !>> Depth in marsh, used in Kadlec-Knight equation - minimum allowed is input parameter to each compartment
          Dmarsh=max(KKdepth(j),Eh(j,2)-BedM(j))

!          Dmarsh=Max(0.01,Eh(j,1)-BedM(j))				!marsh depth
!		Sf = abs(Deta)/hLength(j)						!Gradient towards marsh
!		hann=han(j)*drefm/Dmarsh						!increased resistance as marsh depth-->0
!		Res= (AK/2./g+hLength(j)*hann*hann/Rh**(4/3))/(Ach*Ach)
!		hResist=1.0/sqrt(abs(Res))
!		if(Eh(j,2) <= BedM(j)) then
!              Eh(j,2)=BedM(j)+0.01
!          endif
!		sn=1.0											!positive into marsh
!		if(Deta <= 0) then
!              sn=-1.0							!negative leaving marsh
!		endif
!          QMarsh(j,2)=sqrt(abs(Deta))*hResist*sn*dkd
          MarshEdge = max(hwidth(j),0.01)
          MarshL = hLength(j)

! 2023test ! testing replacing Kadlec Knight marsh flow with Manning's equation
! 2023test ! this will use KKa(j) parameter from cells.csv as Manning's n
          
! 2023test !         QMarshKK=KKa(j)*MarshEdge*(Dmarsh**KKexp)*Detah/MarshL
          
          !>> Calculate marsh flow 'resistance'
		
          !>> Check for missing or negative roughness attribute values (here using the KKa input value for compartment and reassign to default values if missing
          dkd_h = 1.0
          if (KKa(j) < 0.0) then
              KKa(j) = def_n
          elseif (KKa(j) >= 1) then
              if (KKa(j) >= 999.) then
          !>> Set zero multiplier for infinitely high roughness values
                  dkd_h = 0.0
              else
                  KKa(j) = def_n
              endif
          endif

          Marsh_ruf=MarshL*KKa(j)*KKa(j)
          MarshRh=abs(Detah)   ! hyd. radius of marsh is approx. equal to depth
          MarshAch=MarshEdge*abs(Detah)
          MarshRes=(Marsh_ruf/MarshRh**(4./3.))/(MarshAch*MarshAch)	!Resistance in marsh-open-water exchange links
          MarshResist=1.0/sqrt(abs(MarshRes))		!Resistance in marsh-open-water exchange link
          
          if(isNaN(MarshResist)) then
              MarshResist = 0.0
          endif
          
          QMarshMann = sqrt(abs(Detah))*MarshResist*sDetah*dkd_h
      

!>> Determine volume of water available for exchange with marsh during timestep
!          if (sDetaH < 0) then
!              volavailable = (Es(j,2) - Eh(j,2))*As(j,1)
!          else
!              volavailable = (Eh(j,2) - BedM(j))*Ahf(j)
!          endif
!>> Determine maximum flowrate from available exchange water
!          Qmax = volavailable/dt
!>> If calculated flowrate is greater than needed to exchange all available water, set flowrate to the max flowrate calculated in celldQ
          if (Ahf(j) > 0.0) then
! 2023test !              QMarsh(j,2) = sDetaH*min(abs(QMarshKK),Qmarshmax(j))
              QMarsh(j,2) = sDetah*min(abs(QMarshMann),Qmarshmax(j))
          else
              QMarsh(j,2) = 0.0
          endif
	enddo												!JAM Oct 2010

      dmod=thour-int(thour)

      if (dmod == 0) then
		khr=int(thour)
		dmod=khr-3*int(float(khr)/3.)

	endif


! Save salinity values at mid-day to calculate average daily salinities
!	if (hday == 0.0) then
!          ! save mid-day compartment salinities
!          do j=1,N						!JAM Aug 1, 2009 time averaging
!			STEMP(j,2)=S(j,2)
!		enddo
!          ! save mid-day link salinities
!          do kk=1,M
!              SLTEMP(kk,2)=SL(kk,2)
!          enddo
!	endif

!>> generate array of hourly water level in boundary condition compartments (if flag is set to write hourly output file)
      if (dhr == 0.0) then

!>> Write hourly water level in boundary compartments
          if(nstghr > 0) then
		    do jj = 1,nstghr
			    EShrly(jj)=ES(stghrwrite(jj),2)
		    enddo
		    WRITE(210,1111) (EShrly(jj),jj=1,nstghr)
          endif
      endif

!>> Write flow output for select links
      if(nlinksw > 0) then
!>> If first timestep of day, reset average flow term
          if (daystep == 1) then
              do jj = 1,nlinksw
			    FLO(jj) = Q(linkswrite(jj),2)*dt/(3600.*24.)
              enddo
          else
              do jj = 1,nlinksw
                  FLO(jj) = FLO(jj) + Q(linkswrite(jj),2)*dt/(3600.*24.)
              enddo
          endif
!>> Write flow output for links named in 'links_to_write.csv'
          if (daystep == lastdaystep) then
              if (nlinksw == 1) then
                  WRITE(211,1113) (FLO(jj),jj=1,nlinksw)
              else
                  WRITE(211,1112) (FLO(jj),jj=1,nlinksw)
              endif
          endif
      endif


      do kl = 1,N
!>> Reset average WQ values for compartments at start of day
          if(daystep == 1) then
              do klk = 1,14   !zw 4/28/2015 replace ichem with 14
                  ChemAve(kl,klk) = Chem(kl,klk,2)*dt/(3600.*24.)
              enddo
              TempwAve(kl) = Tempw(kl,2)*dt/(3600.*24.)
              ! TSS is sum of 4 different CSS concentrations
              TSSave(kl) = ( CSS(kl,2,1) + CSS(kl,2,2)
     &                    + CSS(kl,2,3) + CSS(kl,2,4) )*dt/(3600.*24.)
              ! Average flow in marsh exchange link
              QmarshAve(kl) = Qmarsh(kl,2)*dt/(3600.*24.)
              !if (kl == 356) then
              !    write(*,*) 'CSS in hydrod:'
              ! write(*,*)CSS(kl,2,1),CSS(kl,2,2),CSS(kl,2,3),CSS(kl,2,4)
              !    write(*,*) 'daystep:',daystep
              !    write(*,*) 'TSSave:',TSSave(kl)
              !    pause
              !endif
			  
			  !zw added 04/15/2020 - move daily values calculations to here from subroutines celldQ/celldSal
		      ESMX(kl,2) = ES(kl,2)
		      ESMN(kl,2) = ES(kl,2)
		      ESAV(kl,1) = ES(kl,2)*dt/(3600.*24.)
              EHAV(kl,1) = EH(kl,2)*dt/(3600.*24.)
              dailyHW(kl) = 0.0
              dailyLW(kl) = 0.0
		      SALAV(kl) = S(kl,2)*dt/(3600.*24.)
!>> Update average WQ values for compartments by timestep's contribution to daily average
          else
              do klk = 1,14  !zw 4/28/2015 replace ichem with 14
                  ChemAve(kl,klk) = ChemAve(kl,klk)
     &                                + Chem(kl,klk,2)*dt/(3600.*24.)
              enddo
              TempwAve(kl) = TempwAve(kl) + Tempw(kl,2)*dt/(3600.*24.)
              TSSAve(kl) = TSSAve(kl) + ( CSS(kl,2,1)+CSS(kl,2,2)
     &                  + CSS(kl,2,3)+CSS(kl,2,4) )*dt/(3600.*24.)
              QmarshAve(kl) = QmarshAve(kl)+Qmarsh(kl,2)*dt/(3600.*24.)
			  
			  !zw added 04/15/2020 - move daily values calculations to here from subroutines celldQ/celldSal
              ESMX(kl,2) = max(ESMX(kl,2),ES(kl,2))
              ESMN(kl,2) = min(ESMN(kl,2),ES(kl,2))
		      ESAV(kl,1) = ESAV(kl,1) + ES(kl,2)*dt/(3600.*24.)
              EHAV(kl,1) = EHAV(kl,1) + EH(kl,2)*dt/(3600.*24.)
              dailyHW(kl) = max(dailyHW(kl),ES(kl,2))
              dailyLW(kl) = min(dailyLW(kl),ES(kl,2))
              EsRange(kl,1)=ESMX(kl,2)-ESMN(kl,2)
 		      SALAV(kl)=SALAV(kl) + S(kl,2)*dt/(3600.*24.)
          endif
      enddo



      do kk=1,M
          if(daystep == 1) then
!>> reset average salinity value for links at start of day
		    SlkAve(kk) = SL(kk,2)*dt/(3600.*24.)
          else
!>> Update average salinity value for links by timestep's contribution to daily average
		    SlkAve(kk)=SlkAve(kk) + SL(kk,2)*dt/(3600.*24.)
	    endif
      enddo



!>> Write daily averages to output files
	if (daystep == lastdaystep) then
!>-- Write daily average water stage elevation in open water to STG.out
          WRITE(111,2222)  (ESAV(j,1),j=1,N)
!>-- Write daily average water stage elevation in marsh to STGm.out
          WRITE(212,2222) (EHAV(j,1),j=1,N)
!>-- Write daily tidal range to TRG.out
          WRITE(112,2222)  (EsRange(j,1),j=1,N)
!>-- Write daily average suspended solids to TSS.out
          WRITE(96,2222) (TSSave(j),j=1,N)
!>-- Write daily average nitrate-N to NO3.out
          WRITE(91,1324)  (ChemAve(j,1),j=1,N)		 !zw 4/28/2015 why multiple the ChemAve by 1000.0? What's the output unit then?
!>-- Write daily average ammonium-N to NH4.out
          WRITE(92,1324)  (ChemAve(j,2),j=1,N)
!>-- Write daily average dissolved organic N to DON.out
          WRITE(113,1324)  (ChemAve(j,10),j=1,N)
!>-- Write daily average dissolved inorganic N to DIN.out
          WRITE(70,1324)  (ChemAve(j,3),j=1,N)
!>-- Write daily average total organic N to OrgN.out
          WRITE(71,1324)  (ChemAve(j,4),j=1,N)
!>-- Write daily total Kjedahl N to TKN.out (TKN = NH4 + OrgN)
          WRITE(123,1324) ((ChemAve(j,2)+ChemAve(j,4)),j=1,N)
!>-- Write daily average total P to TPH.out (TP = TIP + DOP + POP)
          WRITE(72,1324)  ((ChemAve(j,5)+ChemAve(j,11)+ChemAve(j,14))
     &                                ,j=1,N)
!>-- Write daily average total organic C to TOC.out
          !WRITE(73,2222)  (ChemAve(j,10)*1000.*2.,j=1,N) !BUG!? why is Carbon.out calculated as 2*DON?
          WRITE(73,1324)  (ChemAve(j,6),j=1,N) !zw 4/28/2015 TOC as Chem(j,6)
!>-- Write daily average dissolved oxygen to DO.out
          WRITE(95,1324)  (ChemAve(j,7),j=1,N)
!>-- Write daily average live algae to ALG.out
	    WRITE(94,1324)  (ChemAve(j,8),j=1,N)
!>-- Write daily average detritus (dead algae) to DET.out
          WRITE(97,1324)  (ChemAve(j,9),j=1,N)
!	    WRITE(116,2222)  (min(Chem(j,13,2)*1000000.,60.),j=1,N)		!JAM May 21, 2011
!>-- Write daily average to SPH.out
          WRITE(119,1324)  (min(ChemAve(j,1)/3.,max(ChemAve(j,12),0.2
     &                            *ChemAve(j,5))),j=1,N)

          cden=86400000.*365.25
!>-- Write to NRM.out
          WRITE(121,2222)(max(cden*denit(j,2),-39.9),j=1,N)
	    !WRITE(75,2222)(max(min((S(j,2)+STEMP(j,2))/2.,35.),0.2),j=1,N)         !Nov 2010 JAM Aug 1, 2009 time averaging
!>-- Write daily average salinity to SAL.out
          WRITE(75,2222)(SALAV(j),j=1,N)         !zw 3/22/2015 daily average salinity from each time step

!          WRITE(96,2222)  (TSS(j),j=1,N)                                 ! removed - not used !-EDW

!>-- Write daily average marsh-open-water exchange flow to FLOm.out
          WRITE(124,2222) (QmarshAve(j),j=1,N)

!>-- Write daily average water temperature to TMP.out
          WRITE(100,2222)	(TempwAve(j),j=1,N)								!!JAM Oct 24, 2010 should be tempw
!     	WRITE(102,2222)	(age(j,2),j=1,N)								!!JAM Oct 24, 2010 should be age

!>-- Write cumulative sediment accumulation in open water to SedAcc.out (Sacc + AsandA)
          WRITE(103,2227)	(Sacc(j,2),j=1,N)

!>-- Write time marsh is flooded to fflood.out
          WRITE(105,2222)	(floodf(j),j=1,N)								!!JAM Oct 24, 2010 should be Sed accumulation
          ttday=t/(24*3600)

          if(kday == int(ttday+1)) then
              floodf(:)=0.0
          endif

! save daily timeseries for use in developing summary statistics to be used in other ICM routines
          day_int = int(day)    ! !BUG! this will need to be equal to Day of Year, not day (in case multiple years are run at once).
! save daily values at compartments
          do j=1,N
              stage_daily(day_int,j)= ESAV(j,1)
              tidal_range_daily(day_int,j) = dailyHW(j)-dailyLW(j)
              sal_daily(day_int,j) = SALAV(j)
! SALAV is already filtered by filter in celldSal (all filters removed here)
!              sal_daily(day_int,j) = max(min(SALAV(j),35.0),0.2)
!              sal_daily(day_int,j) =
!     &                max(min((S(j,2)+STEMP(j,2))/2.,35.),0.2)        ! salinity equation same as used to write SAL.out file (unit=75)
              tmp_daily(day_int,j) = Tempw(j,2)
              tkn_daily(day_int,j) = ChemAve(j,2)+ChemAve(j,4)
              tss_daily(day_int,j) = TSSAve(j)
              ! Reset daily high and low waters to zero
              dailyHW(j) = 0.0
              dailyLW(j) = 0.0
          enddo



! save daily values at links
          do kk=1,M
              sal_daily_links(day_int,kk) = SlkAve(kk)
          ! SALAV is already filtered by filter in celldSal (all filters removed here)
!          sal_daily_links(day_int,kk)=max(min(SlkAve(kk),35.),0.2)
!              sal_daily_links(day_int,kk) =
!     &                max(min((SL(kk,2)+SLTEMP(kk,2))/2.,35.),0.2)
              tmp_daily_links(day_int,kk) = TL(kk,1)
	    enddo

!          kyear=int(float(kday)/365.25)+1
!          WRITE(*,*)kyear,tmon											!,S(10,2)
! changed the WRITE statement to print the year and day, rather than fractional month values !-EDW
! now read year from input data do not calculate within model


!  This summary output file moved to the end of main.f
!!BUG! this kday flag is never met...therefore ASedOW.csv is never written
!          if (kday >= Nyear*365.2) then
!
!!              pm1= (1.-Pmsh(j))
!          	pm1 = (1.-Apctmarsh(j))
!              WRITE(117,9229)(((ASandA(j)/float(Nyear)/As(j,1)+clams
!     &			+max(0.,Sacc(j,2))/float(Nyear))*pm1				!modified June 20, 2011 JAM for testing removed /1000  added ASANDA(
!     &			+Sacch(j,2)/float(Nyear)*Apctmarsh(j)+476.), j=1,N)		!increase fines captured
!			close (117)
!
!			WRITE(*,*)' day =',kday
!          endif

!          WRITE(103,2227)	(Sacc(j,2),j=1,N)           !!JAM Oct 24, 2010 should be Sed accumulation

      endif               !endif XX

 9229	FORMAT(<cells-1>(F0.4,','),F0.4)
 2227	FORMAT(<cells-1>(F0.4,','),F0.4)
 2222	FORMAT(<cells-1>(F0.4,','),F0.4)
 1324 FORMAT(<cells-1>(F0.6,','),F0.6)
 1111 FORMAT(<nstghr-1>(F0.4,','),F0.4)
 1112 FORMAT(<nlinksw-1>(F0.4,','),F0.4)
 1113 FORMAT(F0.4)
      return
	end

c***********************End Subroutine hydrodynamic*********************************************
