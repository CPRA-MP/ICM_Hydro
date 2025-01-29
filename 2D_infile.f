!> @file
!> @brief This subroutine reads the input text files into Fortran arrays.
!> @details This subroutine reads the various input files from the project directory
!> and saves them into arrays. Some processing of results is also performed.

!> @author J. Alex McCorquodale - University of New Orleans
!> @author Eric White - The Water Institute of the Gulf

      Subroutine infile

      use params

      implicit none
      integer :: i,ichem,ichem2,it,itrib,j,jj,jjn,jn,jt,jkk,mk,kd,kt,etg,rg
      integer :: bedcount,marshcount,sedclass
      integer :: gridbedcount,gridmarcount,gridplcount
      integer :: windstartrun,tidestartrun
      integer :: ETzero,Rzero
      integer :: node,lnkid,jmds,k,jk,ndtt
      real :: Athresh,Area_change,Area_change2,upl,mr,maxmarel,edge_pct
      real :: faN,FSEASON,FSEASON2,FSEASON3,CtoN,fctrib,Qmax,tadd,dlow,a
      integer,dimension(:),allocatable :: jqtrib
      real :: Asum,phi_us,Aus,Ads

!>@par General Structure of Subroutine Logic:
!>> Input junction geometry and properties.
      write(1,*) '----------------------------------------------------'
      write(1,*) 'IMPORTING INPUT DATA'
      write(1,*) '----------------------------------------------------'

      write(*,*) '----------------------------------------------------'
      write(*,*) 'IMPORTING INPUT DATA'
      write(*,*) '----------------------------------------------------'


!>> counters used to report out missing data to console screen
      bedcount = 0
      marshcount = 0
      gridbedcount = 0
      gridmarcount = 0
      gridplcount = 0

      Es(:,:)=0
      S(:,:)=0
      Css(:,:,:)=0
      Tempw(:,:)=0
      Chem(:,:,:)=0
      Eh(:,:)=0

!>> Loop over compartments and save attributes in appropriate arrays.
! change array(j,*)=array(node,*) to ensure the compartment atributes are assigned correctly if the records in cells.csv are not in ascending order - zw 12/11/2023 
      read(32,*)                  ! dump header row of compartment input file
      do j=1,N
          READ(32,*) node,        ! node number - not saved - just a placeholder so input file can have ID as first column
!     &       As(node,1),             ! open water surface area of cells            (km2)
     &        Atotal(node),           ! total area of cells (m2)
     &        Apctwater(node),        ! portion of cell that is water (0-1)
     &        Apctupland(node),       ! portion of cell that is upland (0-1)
     &        Apctmarsh(node),        ! portion of cell that is marsh (0-1)
!     &       Ahf(node),              ! marsh area of cells (m2)
     &        ar_ed(node),            ! marsh edge area of cell (m2)
     &        dump_float,             ! initial stage of storage cells  (m)     !BUG! Eso(node) is not used !BUG! using values from hotstar_in.dat
     &        Bed(node),              ! bed elevation of storage cells  (m)
     &        erBedDepth(node),       ! depth of erodible bed in open water area (m)
     &        erBedBD(node),          ! bulk density of erodible bed in open water area (g/cm3)
!     &       Ahydro(node),           ! area of hydrologically connected area for water balance (marsh+upland) (m2)
     &        Percent(node),          ! percentage open water in Ahydro for PET-to-AET conversion(0-100%)
!     &       dAdz(node),             ! change in As with depth (m)
!     &       CSS(node,1,1),          ! sediment concentration  - sand only (mg/L)  !BUG! CSS(node,1,1) is not used, using values from hotstar_in.dat
     &        adaption_coeff(node),   ! non-equilibrium adaption coefficient for sand resuspension/deposition terms  
     &        dump_float,             ! initial salinity (TDS) (ppt)       !BUG! S(node,1) is not used, using values from hotstar_in.dat
!     &       acss(node),             ! resuspension parameters
!     &       bcss(node),             ! resuspension parameters
!     &       Vss(node),              ! deposition velocity     (m/d)
!     &       tau(node,1),            ! bed shear stress        (Pa)
!     &       tauc(node),             ! critical shear stress   (Pa)
!     &       Fetch(node,1),          ! Fetch E-W               (km)
!     &       Fetch(node,2),          ! Fetch N-S               (km)
!     &       SSource(node),          ! Avg annual source due to coastal erosion (m3/yr/m)
     &        Jrain(node),            ! rain gage number
     &        Jwind(node),            ! wind gage number
     &        Jet(node),              ! ET gage number
     &        ka(node),               ! Wind-current surface coefficient
     &        cf(node),               ! Wind-current bed coefficient
     &        sedn(node),             ! sediment resuspension exponent (non-sand particles)
     &        sedcalib(node),         ! sediment resuspension coefficient (non-sand particles)
     &        alphased(node),         ! sediment vanRijn resuspension coefficient (sand particles)
     &        KKa(node),              ! Kadlec-Knight marsh flow calibration coefficient, A
     &        KKdepth(node),          ! Minimum marsh water depth allowed in Kadlec-Knight marsh flow (m)
     &        MEE(node),              ! Annual marsh edge erosion retreat rate for compartment (m/yr)
!          READ(90,*) nodem,          ! Marsh node number                                       !JAM Oct 2010
!     &       Flood(node),            ! portion of cells that is  (-)
     &        CSSos(node),            ! TSS concentration to be applied to offshore cells - set to -9999 if not used        (mg/L)
     &        BedMOrig(node),         ! bed elevation of marsh storage cells        (m)
     &        BedMSD(node),           ! marsh elevation standard deviation          (m)         ! yw Nov 2020     
     &        Esho(node),             ! soil moisture depth in marsh storage cells  (m)         !JAM Oct 2010
     &        depo_on_off(node),      ! sediment deposition on/off flag - 1=deposition allowed in compartment, 0=no deposition in cell (e.g. main channels) (0,1)
     &        BedMAdj(node)           ! vertical adjustment to marsh elevation (m)
!     &       acssh(node),            ! resuspension parameters                                 !JAM Oct 2010
!     &       Vsh(node),              ! deposition velocity                         (m/d)       !JAM Oct 2010
!     &       SourceBM(node)          ! Avg annual source due to biomass detritus.  (g/yr/m2)   !JAM Oct 2010

          ! Eh(node,1)=Eho(node)  ! BUG! Eh is assigned as BedM+0.1 later on - zw 04/04/2020
          BedM(node) = BedMOrig(node) + BedMAdj(node)
          if ((adaption_coeff(node)<=0) .or.(adaption_coeff(node)>1))  adaption_coeff(node)=0.25

!>> Calculate amount of sediment available in erodible bed (g/m2) for each sediment class (10^6 is unit conversion for bulk density from cm3 to m3)
!>> - once compartment is eroded by this much, only newly deposited sediment can be resuspended
!>> Assume that 10% of bed is sand, 45% is silt, 22.5% is unflocculated clay, and 22.5% is flocculated clay
! This assumes that the erodible bed depth and bulk density are representing all inorganic matter
          erBedAvail(node,1)=erBedDepth(node)*erBedBD(node)*1000000.0*0.1
          erBedAvail(node,2)=erBedDepth(node)*erBedBD(node)*1000000.0*0.45
          erBedAvail(node,3)=erBedDepth(node)*erBedBD(node)*1000000.0*0.225
          erBedAvail(node,4)=erBedDepth(node)*erBedBD(node)*1000000.0*0.225

          !>> Check for missing or negative roughness attribute values (here using the KKa input value for compartment and reassign to default values if missing
          if (KKa(node) < 0.0) then
              KKa(node) = def_n
          elseif ((KKa(node) > 1)) then
              KKa(node) = def_n
          endif

      enddo
      close(32)
          
!>>Loop over compartments and save fetch data in appropriate arrays.          
! change array(j,*)=array(node,*) to ensure the compartment atributes are assigned correctly if the records in fetch.csv are not in ascending order - zw 12/11/2023 
      Fetch(:,:)=0
      do j=1,N
          READ(323,*) node,      ! node number
     &        Fetch(node,1),         ! Fetch length (m) in 0-deg sector
     &        Fetch(node,2),         ! Fetch length (m) in 22.5-deg sector
     &        Fetch(node,3),         ! Fetch length (m) in 45-deg sector
     &        Fetch(node,4),         ! Fetch length (m) in 67.5-deg sector
     &        Fetch(node,5),         ! Fetch length (m) in 90-deg sector
     &        Fetch(node,6),         ! Fetch length (m) in 112.5-deg sector
     &        Fetch(node,7),         ! Fetch length (m) in 135-deg sector
     &        Fetch(node,8),         ! Fetch length (m) in 157.5-deg sector
     &        Fetch(node,9),         ! Fetch length (m) in 180-deg sector
     &        Fetch(node,10),        ! Fetch length (m) in 202.5-deg sector
     &        Fetch(node,11),        ! Fetch length (m) in 225-deg sector
     &        Fetch(node,12),        ! Fetch length (m) in 247.5-deg sector
     &        Fetch(node,13),        ! Fetch length (m) in 270-deg sector
     &        Fetch(node,14),        ! Fetch length (m) in 292.5-deg sector
     &        Fetch(node,15),        ! Fetch length (m) in 315-deg sector
     &        Fetch(node,16)         ! Fetch length (m) in 337.5-deg sector
      enddo
      close(323)
      
!>> Add minimal area of water portion - each cell must have at least some water
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Updating area attributes to meet minimum water area:'
      write(1,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Updating area attributes to meet minimum water area:'
      write(*,*) '----------------------------------------------------'

      do j = 1,N
      !>> check for incorrect area percentages
          if (Apctwater(j) < 0.0) then
              Apctwater(j) = 0.0
          endif

          if (Apctmarsh(j) < 0.0) then
              Apctmarsh(j) = 0.0
          endif

          if (Apctupland(j) < 0.0) then
              Apctupland(j) = 0.0
          endif

          if (Atotal(j) < minwater) then
              Atotal(j) = minwater
              write(1,919)'Total area for compartment',j,'increased to
     & ',minwater,'sq m.'
              write(*,919)'Total area for compartment',j,'increased to
     & ',minwater,'sq m.'
          endif

          Athresh = minwater/Atotal(j)
!>> If total area is less than or equal to the mimimum water area required, set entire compartment to water
          if(Athresh >= 1) then
              Apctwater(j) = 1.0
              Apctmarsh(j) = 0.0
              ar_ed(j) = 0.0
              Apctupland(j) = 0.0
              BedM(j) = Bed(j)    !change marsh elev to bed elevation if converted entire marsh area to water
!>> Otherwise, determine amount of upland and marsh needed to be converted to water to meet minimum requirement
          else
              if (Apctwater(j) <= Athresh) then
                  Area_change = Athresh - Apctwater(j)
                  upl = Apctupland(j)
                  mr = Apctmarsh(j)
                  edge_pct = ar_ed(j)/(Apctmarsh(j)*Atotal(j))
                  Apctwater(j) = Athresh

!>> Determine if there is enough upland area to only reduce upland area for increase to water area
                  if (Area_change <= Apctupland(j)) then
                      Apctupland(j) = upl - Area_change
                      write(1,919) 'Compartment',j,' water area
     & increased by',Area_change*Atotal(j),'sq m of converted upland.'
                      write(*,919)  'Compartment',j,' water area
     & increased by',Area_change*Atotal(j),'sq m of converted upland.'
                  else
!>> If there is not enough upland, set upland to zero and reduce marsh area by remainder

                      write(1,919) 'Compartment',j,' water area
     & increased by',Apctupland(j)*Atotal(j),'sq m of converted upland.'
                      write(*,919)  'Compartment',j,' water area
     & increased by',Apctupland(j)*Atotal(j),'sq m of converted upland.'

                      Apctupland(j) = 0.0
                      Area_change2 = Area_change - upl
                      Apctmarsh(j) = mr - Area_change2
                      ar_ed(j) = edge_pct*Apctmarsh(j)

                      write(1,919) 'Compartment',j,' water area
     & increased by',Area_change2*Atotal(j),'sq m of converted marsh.'
                      write(*,919)  'Compartment',j,' water area
     & increased by',Area_change2*Atotal(j),'sq m of converted marsh.'

                   endif
              endif
          endif
      enddo
919   Format(7x,a,x,I0,x,a,x,F0.0,x,a)
920   Format(7x,F0.2,x,a,x,I0)
921   Format(7x,F0.2,x,a,x,F0.2,x,a,x,I0)

!>> check that compartment bed and marsh elevations have values, if not fill with default
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Updating bed and marsh elevations if needed:'
      write(1,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Updating bed and marsh elevations if needed:'
      write(*,*) '----------------------------------------------------'
      do j=1,N
          if (bed(j) < -9990) then
              bed(j) = defbedelev
              bedcount = bedcount + 1
          endif
          if (bedM(j) < -9990) then
              bedM(j) = defmarelev
              marshcount = marshcount + 1
          endif

          if (BedM(j) <= Bed(j)) then
              BedM(j) = Bed(j)
              write(1,922) 'Compartment',j,'had an input marsh
     & elevation below the open water bed elevation. Set to bed elev.'
              write(*,922) 'Compartment',j,'had an input marsh
     & elevation below the open water bed elevation. Set to bed elev.'
          endif
      enddo

922   Format(7x,a,x,I0,x,a)
!>> Check for areas of very small marsh.
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Removing small marsh areas, if needed:'
      write(1,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Removing small marsh areas, if needed:'
      write(*,*) '----------------------------------------------------'

      do j = 1,N
          if (Apctmarsh(j) > 0.0) then
              if (Apctmarsh(j) < 0.01) then
                  Apctwater(j) = Apctwater(j) + Apctmarsh(j)
                  Apctmarsh(j) = 0.0
                  ar_ed(j) = 0.0
                  BedM(j) = Bed(j)    !change marsh elev to bed elevation if converted entire marsh area to water
                  write(1,923) 'Compartment',j,'had less than 1% marsh.
     & All marsh converted to open water.'
                  write(*,923) 'Compartment',j,'had less than 1% marsh.
     & All marsh converted to open water.'
              endif
          endif

!>> Calculate water/marsh/upland areas from input percentages.
          As(j,1) = Apctwater(j)*Atotal(j)
          As(j,2) = As(j,1)
          Ahf(j) = Apctmarsh(j)*Atotal(j)
          Ahydro(j) = (Apctupland(j)+Apctmarsh(j))*Atotal(j)
!          Ahydro(j)=Ahydro(j)*1000000.
!          Ahf(j)=flood(j)*Ahydro(j)        ! JAM Oct 2010
!          Vss(j)=Vss(j)/consd              ! JAM Oct 2010
!          Vsh(j)=Vsh(j)/consd              ! JAM Oct 2010
!          ES(j,1)=ES(j,1)+RSSSo  !BUG RSSSo is undefined - zw 04/06/2020


!>> Calculate length of marsh edge from edge area (for KadlecKnight equation)
!>> Set edge length to perimeter of idealized (e.g. square)  marsh or water area, whichever is smaller
          hwidth(j) = min(sqrt(As(j,2)),sqrt(Ahf(j)))
!>> The raster used to calculate edge area is 30.0 m square, therefore the length of marsh is
!>> assumed to be equal to the area divided by this dimension
!         hwidth(j) = ar_ed(j)/30.0

!>> Calculate distance from marsh edge to center of marsh (for KadlecKnight equation)
!>> this assumes that the marsh area in cells.csv is a square border surrounding a square water compartment
!>> Set minimum value of marsh length equal to 30.0 m (corresponds to pixel size of land/water raster in ICM)
          hLength(j) = max(30.0,(sqrt(Ahf(j)+ As(j,1))-Hwidth(j))/4.0)
          ar_int(j) = max(0.0,Ahf(j) - ar_ed(j))
!          rmo=sqrt(As(j,1)/pi)             ! JAM Oct 2010
!          rmm=sqrt((As(j,1)+Ahf(j))/pi)    ! JAM Oct 2010
!          hLength(j)=(rmm-rmo)*0.5         ! JAM Oct 2010
!          hwidth(j)=2.*pi*rmo              ! JAM Oct 2010
      enddo

923   Format(7x,a,x,I0,x,a)

!      clams=40.            !   clam grow rate          (g/m2/yr) !it is not used anywhere

! input link geometry   and properties of OPEN WATER

!>> Loop over links and save input data in appropriate arrays.
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Reading link attributes:'
      write(1,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Reading link attributes:'
      write(*,*) '----------------------------------------------------'
      read(33,*)          ! dump header row of links input file
! change array(i,*)=array(lnkid,*) to ensure the compartment atributes are assigned correctly if the records in fetch.csv are not in ascending order - zw 12/11/2023 
      do i=1,M
          READ(33,*) lnkid,         ! link number - not saved - just a placeholder so input file can have ID as first column
     &        jus(lnkid),                  ! upstream compartment number
     &        jds(lnkid),                  ! downstream compartment number
!     &       itype(lnkid),                ! "-1" is u/s BC, "+1" is d/s BC, "2" is two cell link.
     &        USx(lnkid),                  ! x coordinate (UTM meters) of upstream end of link
     &        USy(lnkid),                  ! x coordinate (UTM meters) of upstream end of link
     &        DSx(lnkid),                  ! x coordinate (UTM meters) of downstream end of link
     &        DSy(lnkid),                  ! x coordinate (UTM meters) of downstream end of link
     &        linkt(lnkid),                ! link type (1=channel,2=weir,3=lock,4=lock,5=orifice,6=culvert/bridge,7=pump,8=marsh overland flow,9=ridge/levee,10=regime channel for sediment)
     &        Latr1(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr2(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr3(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr4(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr5(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr6(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr7(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr8(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Exy(lnkid),                  ! diffusion coefficient thru link between cells (m2/s)
     &        Latr9(lnkid),                ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        Latr10(lnkid),               ! link attribute (varies by link type - see Hydrod.f for descriptions)
     &        fa_mult(lnkid)               ! link multiplier to apply to default upwind factor for diffusivity calculations

!>> flag overland marsh and overland ridge links to turn off if model overland flow flag is 0
          if (modeloverland == 0) then
              if (linkt(lnkid) == 8) then
                  linkt(lnkid) = -8
              elseif (linkt(lnkid) == 9) then
                  linkt(lnkid) = -9
              endif
          write(1,*)
          write(1,*) '       Overland flow links turned off.'
          write(1,*)
          write(*,*)
          write(*,*) '       Overland flow links turned off.'
          write(*,*)
          endif
      enddo
      close(33)

!============ check link attributes & assign default values if missing or violeting possible ranges 
      do lnkid=1,M

!>> channel/culvert link attribute checks
          if ((linkt(lnkid) == 1) .or. (linkt(lnkid) == 3) 
     &        .or. (linkt(lnkid) == 6) .or. (linkt(lnkid) == 11)  
     &        .or. (linkt(lnkid) == 12)) then
              
!zw 1/30/2024              Latr1(lnkid)=max(Bed(jus(lnkid)),Bed(jds(lnkid)))

              if((linkt(lnkid)==6) .and. (Latr2(lnkid)<=Latr1(lnkid)))then
                  write(1,925) 'Culvert link',lnkid,'has crown elevation lower than 
     & invert elevation'
                  write(*,925) 'Culvert link',lnkid,'has crown elevation lower than 
     & invert elevation'
                  Latr2(lnkid)=Latr1(lnkid)+3.0
                  write(1,*) 'Crown elevation is set to be 3m 
     & above crest elevation: ',Latr2(lnkid)
                  write(*,*) 'Crown elevation is set to be 3m 
     & above crest elevation: ',Latr2(lnkid)
              endif

              if(linkt(lnkid)==3)then
                  if (Latr9(lnkid)<=0) then
                      write(1,925)'Lock control scheme of link',lnkid,'is missing'
                      write(*,925)'Lock control scheme of link',lnkid,'is missing'
                      stop
                  endif
                  
                  if(Latr2(lnkid)<0) then
                      if((Latr9(lnkid)==5) .or. (Latr9(lnkid)==10)
     &                    .or. (Latr9(lnkid)==11)) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing salinity (Latr2)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing salinity (Latr2)'
                          stop
                      elseif (Latr9(lnkid)==7) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing control link number (Latr2)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing control link number (Latr2)'
                          stop
                      endif
                  endif
  
                  if(Latr10(lnkid)<0) then
                      if((Latr9(lnkid)==1) .or. (Latr9(lnkid)==3)
     &                    .or. (Latr9(lnkid)==5) .or. (Latr9(lnkid)==8)
     &                    .or. (Latr9(lnkid)==9) .or. (Latr9(lnkid)==10)
     &                    .or. (Latr9(lnkid)==11)) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing stage (Latr10)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing stage (Latr10)'
                          stop
                      elseif(Latr9(lnkid)==4) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing salinity (Latr10)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing salinity (Latr10)'
                          stop
                      elseif (Latr9(lnkid)==6) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing corresponding column in LockControlObservedData.csv (Latr10)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing corresponding column in LockControlObservedData.csv (Latr10)'
                          stop
                      elseif (Latr9(lnkid)==7) then
                          write(1,925)'Lock control scheme of link',lnkid,
     & 'missing discharge (Latr10)'
                          write(*,925)'Lock control scheme of link',lnkid,
     & 'missing discharge (Latr10)'
                          stop
                      endif
                  endif
              endif

          !>> Set channel length to default value if missing
              if (Latr3(lnkid) <= 0.0) then
                  write(1,*)'Length (Latr3) for link',lnkid,'is missing'
                  write(1,*) 'Default value of 1000m is assigned
     & but this should be fixed in the input file.'

                  write(*,*)'Length (Latr3) for link',lnkid,'is missing'
                  write(*,*) 'Default value of 1000m is assigned
     & but this should be fixed in the input file.'
                  Latr3(lnkid) = 1000.0
              endif

              if (Latr4(lnkid) <= 0.0) then
                  write(1,*)'Width (Latr4) for link',lnkid,'is missing'
                  write(1,*) 'Default value of 3m is assigned
     & but this should be fixed in the input file.'

                  write(*,*)'Width (Latr4) for link',lnkid,'is missing'
                  write(*,*) 'Default value of 3m is assigned
     & but this should be fixed in the input file.'
                  Latr4(lnkid) = 3.0
              endif

          !>> Check for missing roughness attribute values and reassign to default values if missing
              if (Latr5(lnkid) <= 0.0) then
                  Latr5(lnkid) = def_n
              elseif (Latr5(lnkid) >= 1) then
                  Latr5(lnkid) = def_n
              endif

              Latr6(lnkid)=max(Latr6(lnkid),0.0)
              Latr7(lnkid)=max(Latr7(lnkid),0.0)
              Latr8(lnkid)=max(Latr8(lnkid),0.0)

              !ZW 2/1/2024 calculate fa as fa(us) for blended differencing (BD) scheme to determine link face salinity
              if (iAdvTrans==1) then
                  Asum=As(jus(lnkid),1)+As(jds(lnkid),1)
                  phi_us=As(jus(lnkid),1)/Asum
                  fa_mult(lnkid)=1.0-r_BD*phi_us  !this is fa for flow from US to DS (Q>0)
              endif
!>> weir link attribute checks
          elseif(linkt(lnkid) == 2) then
!>> set weir coefficient to default value if outside of standard range
              if (Latr8(lnkid) < 0.55) then
                  Latr8(lnkid) = 0.62
              elseif (Latr8(lnkid) > 0.65) then
                  Latr8(lnkid) = 0.65
              endif
!>> check for missing weir invert attributes and assign to compartment bed elev
!              if (Latr2(lnkid) < -9990.) then
!                  Latr2(lnkid) = Latr3(lnkid)
!              elseif (Latr2(lnkid) > 100.) then
!                  Latr2(lnkid) = Latr3(lnkid)
!              endif
!   weir upstream & downstream ground elevation = bed elevation of corresponding us/ds compartment
              Latr2(lnkid)=Bed(jus(lnkid))
              Latr3(lnkid)=Bed(jds(lnkid))
!   weir upstream & downstream ground elevation can not be higher than crest elevation
              if(Latr1(lnkid)<=max(Latr2(lnkid),Latr3(lnkid)))then
                  write(1,925) 'Weir link',lnkid,'has crest elevation lower than 
     & bed elevation of connecting compartments'
                  write(*,925) 'Weir link',lnkid,'has crest elevation lower than 
     & bed elevation of connecting compartments'
                  Latr1(lnkid)=max(Latr2(lnkid),Latr3(lnkid))+0.5
                  write(1,*) 'Weir crest elevation is set to be 0.5m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
                  write(*,*) 'Weir crest elevation is set to be 0.5m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
              endif
              if(Latr4(lnkid) <=0)then
                  write(1,925) 'Weir link',lnkid,'has crest length lower than 0'
                  write(*,925) 'Weir link',lnkid,'has crest length lower than 0'
                  Latr4(lnkid)=10.0
                  write(1,*) 'Weir crest length is set to be 10m'
                  write(*,*) 'Weir crest length is set to be 10m'
              endif

!!!ZW 12/15/2023 upwind factor for weir link should always be 1???
              fa_mult(lnkid) = 1.0  

!>> Tidal gate/orifice link attribute checks
          elseif((linkt(lnkid) == 4) .or. (linkt(lnkid) == 5)) then
!   orifice upstream & downstream ground elevation = bed elevation of corresponding us/ds compartment
!             Latr3(lnkid)=bed(jus(lnkid))
!             Latr5(lnkid)=bed(jds(lnkid))

!   orifice upstream & downstream ground elevation can not be higher than structure invert elevation
              if(Latr1(lnkid)<=max(Latr3(lnkid),Latr5(lnkid)))then
                  write(1,925) 'Orifice/tidal gate link',lnkid,'has invert elevation lower than 
     & bed elevation of connecting compartments'
                  write(*,925) 'Orifice/tidal gate link',lnkid,'has invert elevation lower than 
     & bed elevation of connecting compartments'
                  Latr1(lnkid)=max(Latr3(lnkid),Latr5(lnkid))+0.5
                  write(1,*) 'Invert elevation is set to be 0.5m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
                  write(*,*) 'Invert elevation is set to be 0.5m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
              endif

              if(Latr2(lnkid)<=Latr1(lnkid))then
                  write(1,925) 'Orifice/tidal gate link',lnkid,'has crown elevation lower than 
     & invert elevation'
                  write(*,925) 'Orifice/tidal gate link',lnkid,'has crown elevation lower than 
     & invert elevation'
                  Latr2(lnkid)=Latr1(lnkid)+1.0
                  write(1,*) 'Crown elevation is set to be 1m 
     & above crest elevation: ',Latr2(lnkid)
                  write(*,*) 'Crown elevation is set to be 1m 
     & above crest elevation: ',Latr2(lnkid)
              endif

              if(Latr4(lnkid) <=0)then
                  write(1,925) 'Orifice/tidal gate link',lnkid,'has mean width lower than 0'
                  write(*,925) 'Orifice/tidal gate link',lnkid,'has mean width lower than 0'
                  Latr4(lnkid)=5.0
                  write(1,*) 'Tidegate/Orifice width is set to be 5m'
                  write(*,*) 'Tidegate/Orifice width is set to be 5m'
              endif
              if (Latr8(lnkid)<0) Latr8(lnkid)=0.4
!!!ZW 12/15/2023 upwind factor for orifice link should always be 1???
              fa_mult(lnkid) = 1.0  

!>> pump link attribute checks
          elseif(linkt(lnkid) == 7) then
!   pump upstream stage threshold can not be lower than bed/marsh elevation
              if(Latr1(lnkid)<=bed(jus(lnkid)))then
                  write(1,925) 'Pump link',lnkid,'has upstream pumpon 
     & stage threshold lower than bed elevation'
                  write(*,925) 'Pump link',lnkid,'has upstream pumpon 
     & stage threshold lower than bed elevation'
                  Latr1(lnkid)=bed(jus(lnkid))+dry_threshold
              endif
              if(Latr2(lnkid)<=bed(jus(lnkid)))then
                  write(1,925) 'Pump link',lnkid,'has upstream pumpoff 
     & stage threshold lower than bed elevation'
                  write(*,925) 'Pump link',lnkid,'has upstream pumpoff 
     & stage threshold lower than bed elevation'
                  Latr2(lnkid)=bed(jus(lnkid))+dry_threshold
              endif
              if(Latr1(lnkid)<Latr2(lnkid))then
                  write(1,925) 'Pump link',lnkid,'has pumpon stage 
     & threshold lower than pumpoff stage threshold'
                  write(*,925) 'Pump link',lnkid,'has pumpon stage 
     & threshold lower than pumpoff stage threshold'
                  Latr1(lnkid)=Latr2(lnkid)
              endif
              if(Latr9(lnkid) <=0)then
                  write(1,925) 'pump link',lnkid,'has pump capacity lower than 0'
                  write(*,925) 'pump link',lnkid,'has pump capacity lower than 0'
                  Latr9(lnkid) = 20.0
                  write(1,*) 'Pump capacity is set to 20 m3/s'
              endif

!!!ZW 12/15/2023 upwind factor for weir link should always be 1 & no diffusion term
              fa_mult(lnkid) = 1.0  
              Exy(lnkid) = 0.0

!>> Marsh link attribute checks
          elseif (linkt(lnkid) == 8) then
              Latr2(lnkid)=BedM(jus(lnkid))
              Latr10(lnkid)=BedM(jds(lnkid))

!>> if marsh overland links connect to a compartment that has zero marsh area, update that compartment's marsh elevation (in the link attributes) to the bed elevation of the open water
! check if upstream is now water
              if (Ahf(jus(lnkid)) == 0.0) then
                  Latr2(lnkid) = Bed(jus(lnkid))
                  write(1,924) 'Overland link',lnkid,'no longer connects
     & to marsh upstream. Compartment',jus(lnkid),'is now all water.'
                  write(*,924) 'Overland link',lnkid,'no longer connects
     & to marsh in upstream. Compartment',jus(lnkid),'is now all water.'
              endif

! check if downstream is now water
! this condition is met if BOTH us/ds are now water
              if (Ahf(jds(lnkid)) == 0.0) then
                  Latr10(lnkid) = Bed(jds(lnkid))
                  write(1,924) 'Overland link',lnkid,'no longer connects
     & to marsh downstream. Compartment',jds(lnkid),'is now all water.'
                  write(*,924) 'Overland link',lnkid,'no longer connects
     & to marsh in downstream. Compartment',jds(lnkid),'is now all water.'
              endif
! Set the elevation of the marsh overland link to be the maximum of the input link elevation, and the us and downstream marsh elevations (that may or may not have been updated for non-marsh areas)
! if both us/ds are now water, Latr1 was reset, otherwise the input value for latr1 will still be used in comparison of marsh elev

! both ends are water, update roughness value to low (water) value of 0.03
              if((Ahf(jus(lnkid)) == 0.0) .and. (Ahf(jds(lnkid)) == 0.0)) then 
                      Latr1(i) = max(Latr2(i),Latr10(i))
                      Latr5(lnkid) = 0.03
                      write(1,925) 'Overland link',lnkid,'no longer connects
     & to marsh at either end. Roughness is set to 0.03.'
                      write(*,925) 'Overland link',lnkid,'no longer connects
     & to marsh at either end. Roughness is set to 0.03.'
              endif
!              maxmarel = max(Latr1(lnkid),Latr2(lnkid),Latr10(lnkid))
!              Latr1(lnkid) = maxmarel
              Latr1(lnkid) = max(Latr2(lnkid),Latr10(lnkid))  !ZW 12/12/2023 marsh link invert should not higher than us & ds marsh/bed elevations

              if(Latr3(lnkid) <=0)then
                  write(1,925) 'Marsh link',lnkid,'has length lower than 0'
                  write(*,925) 'Marsh link',lnkid,'has length lower than 0'
                  Latr3(lnkid) = 1000.0
                  write(1,*) 'Marsh link length is set to be 1000m'
                  write(*,*) 'Marsh link length is set to be 1000m'
              endif

              if(Latr4(lnkid) <=0)then
                  write(1,925) 'Marsh link',lnkid,'has width lower than 0'
                  write(*,925) 'Marsh link',lnkid,'has width lower than 0'
                  Latr4(lnkid) = 1000.0
                  write(1,*) 'Marsh link width is set to be 1000m'
                  write(*,*) 'Marsh link width is set to be 1000m'
              endif

              if (Latr5(lnkid) <= 0.0) then
                  Latr5(lnkid) = 0.05
              elseif (Latr5(lnkid) >= 1) then
                  Latr5(lnkid) = 0.05
              endif
              !ZW 2/1/2024 calculate fa as fa(us) for blended differencing (BD) scheme to determine link face salinity
              if (iAdvTrans==1) then
                  Aus=Ahf(jus(lnkid))
                  Ads=Ahf(jds(lnkid))
                  If(Ahf(jus(lnkid))==0) Aus=As(jus(lnkid),1)
                  If(Ahf(jds(lnkid))==0) Ads=As(jds(lnkid),1)
                  Asum=Aus+Ads
                  phi_us=Aus/Asum
                  fa_mult(lnkid)=1.0-r_BD*phi_us  !this is fa for flow from US to DS (Q>0)
              endif

!>> ridge link attribute checks
          elseif(linkt(lnkid) == 9) then
!   ridge upstream & downstream ground elevation = bed elevation of corresponding us/ds compartment
              Latr2(lnkid)=BedM(jus(lnkid))
              Latr10(lnkid)=BedM(jds(lnkid))
!   weir upstream & downstream ground elevation can not be higher than crest elevation
              if(Latr1(lnkid)<=max(Latr2(lnkid),Latr10(lnkid)))then
                  write(1,925) 'Ridge link',lnkid,'has crest elevation lower than 
     & bed elevation of connecting compartments'
                  write(*,925) 'Ridge link',lnkid,'has crest elevation lower than 
     & bed elevation of connecting compartments'
                  Latr1(lnkid)=max(Latr2(lnkid),Latr10(lnkid))+1.0
                  write(1,*) 'Ridge crest elevation is set to be 1m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
                  write(*,*) 'Ridge crest elevation is set to be 1m 
     & above higher bed elevation of us/ds compartments: ',Latr1(lnkid)
              endif

              if(Latr3(lnkid) <=0)then
                  write(1,925) 'Ridge link',lnkid,'has crest width 
     & (parallel to flow) lower than 0'
                  write(*,925) 'Ridge link',lnkid,'has crest width 
     & (parallel to flow) lower than 0'
                  Latr3(lnkid) = 10.0  !default 10m ridge crest width
                  write(1,*) 'Ridge link crest width is set to be 10m'
                  write(*,*) 'Ridge link crest width is set to be 10m'
              endif
!              if(Latr3(lnkid) >30)then
!                  write(1,925) 'Ridge link',lnkid,'has crest width longer than 30m'
!                  write(*,925) 'Ridge link',lnkid,'has crest width longer than 30m'
!                 Latr3(lnkid) = 10.0  !default 10m ridge crest width
!                  write(1,*) 'Ridge link crest width is set to be 10m'
!                  write(*,*) 'Ridge link crest width is set to be 10m'
!             endif

              if(Latr4(lnkid) <=0)then
                  write(1,925) 'Ridge link',lnkid,'has crest length 
     & (perpendicular to flow) lower than 0'
                  write(*,925) 'Ridge link',lnkid,'has crest length 
     & (perpendicular to flow) lower than 0'
                  Latr4(lnkid)=100.0
                  write(1,*) 'Ridge link crest length is set to be 100m'
                  write(*,*) 'Ridge link crest length is set to be 100m'
              endif

              if (Latr5(lnkid) <= 0.0) then
                  Latr5(lnkid) = 0.1
              elseif (Latr5(lnkid) >= 1) then
                  Latr5(lnkid) = 0.1
              endif

              if(Latr8(lnkid)<0) Latr8(lnkid) = 0.37  !default for long broad-crested weir
!!!ZW 12/15/2023 upwind factor for weir link should always be 1???
              fa_mult(lnkid) = 1.0  
          endif

!>> Exy check
          if(Exy(lnkid) < 0) then   !zw 3/14/2015 revised
              write(1,925)'Exy value for link',lnkid,'is missing.'
              write(1,*) 'Default value of 10 is assigned
     &    but this should be fixed in the input file.'

              write(*,925)'Exy value for link',lnkid,'is missing.'
              write(*,*) 'Default value of 10 is assigned
     &    but this should be fixed in the input file.'

              Exy(lnkid) = 10.0
          endif

!>> fa check
          if(fa_mult(lnkid) < 0.5) then !zw 3/14/2015 revised
              write(1,925)'fa value for link',lnkid,'is less than 0.5'
              write(1,*) 'Default value of 0.5 is assigned'

              write(*,925)'fa value for link',lnkid,'is missing.'
              write(*,*) 'Default value of 0.5 is assigned'

              fa_mult(lnkid) = 0.56
          endif

! Latr11 is for type 3 Lock Links, not in the attribute inputs.  Used in hydrod.f         
          Latr11(lnkid) = 0.0            
      enddo
!================end link attribute checks
      
924   Format(7x,a,x,I0,x,a,x,I0,x,a)
925   Format(7x,a,x,I0,x,a)
      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Reading various model configuration settings:'
      write(1,*) '----------------------------------------------------'
      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Reading various model configuration settings:'
      write(*,*) '----------------------------------------------------'

!>> Dump header row of LinksClosedHours file
      hourclosed(:,:)=0
      read(34,*)
!>> Read LinksClosedHours file to create array of hours where links are closed (only used for lock-type links)
!>> This hourly control pattern is repeated for every day of the simulation
! change array(i,*)=array(lnkid,*) to ensure the compartment atributes are assigned correctly if the records in fetch.csv are not in ascending order - zw 12/11/2023 
      do i=1,M
          READ(34,*) lnkid,        ! read link number - don't save
     &        hourclosed(lnkid,1),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,2),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,3),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,4),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,5),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,6),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,7),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,8),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,9),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,10),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,11),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,12),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,13),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,14),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,15),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,16),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,17),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,18),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,19),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,21),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,22),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,23),        ! flag value 1 if link is closed during hour
     &        hourclosed(lnkid,24)        ! flag value 1 if link is closed during hour
      enddo
      close(34)

!>> Read in link numbers for locks that have a  timeseries in the hourly timeseries control rule file
!>> This hourly control timeseries will override any other lock control rules
!>> count number of locks that have observed lock control data
      nlockobs_links = 0
      nlockobs = 0
      do i = 1,M
          if (linkt(i) == 3) then
              if (Latr9(i) == 6) then
                  nlockobs_links = nlockobs_links + 1
                  nlockobs = max(nlockobs,int(Latr10(i)))      ! read in the observation record that the link will be looking for set the max value to nlockobs to be used for array allocation
              endif
          endif
      enddo

      if (nlockobs_links > 0) then
          write(*,4123) nlockobs_links,
     &       'links are classified as locks with observed control data.'
          write(1,4123) nlockobs_links,
     &       'links are classified as locks with observed control data.'

          if (nlockobs > 0) then
              write(*,4123) nlockobs,
     &            'lock observation records were found.'
              write(1,4123) nlockobs,
     &            'lock observation records were found.'
          endif
            
          allocate(lockhours(simdays*24/dtlock,nlockobs))
          lockhours(:,:)=0 !zw added 04/06/2020

          read(35,*) dump_text ! dump header row

          !>> Skip Lock operation observations that occur in years prior to current model year
          do jt=1,startrun*24/dtlock
              read(35,*)
              if (jt == startrun) then
                  write(1,*)'Lock observation data starts at:'
                  write(1,66) '      line ',int(startrun*24/dtlock+2)
                  write(*,*)'Lock observation data starts at:'
                  write(*,66) '      line ',int(startrun*24/dtlock+2)
              endif
          enddo

          !>> Read lock control operation observations for current model year
          do jt = 1,simdays*24/dtlock
              READ(35,*) dump_int,dump_text,lockhours(jt,1:nlockobs)
          enddo
      endif
      close(35)

      !>> Read in links that will have flowrates printed to FLO.out file
      linkswrite(:)=0    !zw added 04/06/2020
      if (nlinksw>0) then
          read(126,*)
          do i = 1,nlinksw
              read(126,*) linkswrite(i)
          enddo
          close(126)
    
          write(*,4123) nlinksw,
     &    'links will have average daily flow output written to FLO.out'
          write(1,4123) nlinksw,
     &    'links will have average daily flow output written to FLO.out'
      endif

!>> Read in compartments that will have hourly stage printed to STGhr.out file
      stghrwrite(:)=0   !zw added 04/06/2020
      if (nstghr>0) then
          read(127,*)
          do i = 1,nstghr
              read(127,*) stghrwrite(i)
          enddo
          close(127)

          write(*,4123) nstghr,'compartments will
     &     have average hourly flow output written to STGhr.out'
          write(1,4123) nstghr,'compartments will
     &     have average hourly flow output written to STGhr.out'
      endif
      
4123  Format(7X,I0,X,A)

      write(1,*)
      write(1,*) '----------------------------------------------------'
      write(1,*) 'Reading input data timeseries:'
      write(1,*) '----------------------------------------------------'
      write(1,*)'Only data for current year of ICM run will be imported'
      write(1,*) 'Current ICM model year =', year
      write(1,*)

      write(*,*)
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Reading input data timeseries:'
      write(*,*) '----------------------------------------------------'
      write(*,*)'Only data for current year of ICM run will be imported'
      write(*,*) 'Current ICM model year =', year
      write(*,*)


      Qtrib(:,:)=0.0
      Qmult(:,:)=0
      jtrib(:)=0
!>> Read Tributary and Diversion Numbers from header row of input file
      READ(39,*) (jtrib(jt), jt=1,Ntrib)
      do jt=1,Ntrib
          if((jtrib(jt)>Ntrib) .or. (jtrib(jt)<1))then
              write(*,*)'The tributary number of column in TribQ ',jt,
     &             'is ',jtrib(jt), '. It should be between 1 & ',Ntrib
              stop
          endif
      enddo

!>> Skip Tributary and Diversion Flows that occur in years prior to current model year
      do kt=1,startrun
          read(39,*)
          if (kt == startrun) then
              write(1,*)'Tributary flow input data starts at:'
              write(1,66) '      line ',int(startrun+2)
              write(*,*)'Tributary flow input data starts at:'
              write(*,66) '      line ',int(startrun+2)
          endif
      enddo

!>> Read Tributary and Diversion Flows for current model year
      do kt=1,simdays
          READ(39,*)(Qtrib(jtrib(jt),kt),jt=1,Ntrib)
      enddo
      close(39)
!>> Read Tributary Matrix
    ! dump header row of tributary and diversion multiplier matrix files (Qmult, Qmult_div)
      read(77,*)

      do j=1,N
          READ(77,*) node,(Qmult(node,jn), jn=1,Ntrib)   !dumps first column (which is compartment number)
          do jn = 1,Ntrib
              if (isnan(Qmult(node,jn))) then
                  Qmult(node,jn) = 0.0
              endif
          enddo
      enddo
      close(77)


!>> Read multiplier for diversion flows as function of Mississippi River flow
      DivMult(:)=0  !zw added 04/06/2020
      Qdiv(:,:)=0
      Qmultdiv(:,:)=0  !zw added 04/06/2020
      if (Ndiv>0) then
          read(86,*) !skip header row
          write(1,*) 'Reading multiplier on flow for each diversion.'
          write(*,*) 'Reading multipliers on flow for each diversion.'
          do it = 1,Ndiv
              READ(86,*) dump_int,DivMult(it)
          enddo

!>> Develop timeseries of diversion flows
          do kt = 1,simdays
          do it = 1,Ndiv
              Qdiv(it,kt) = Qtrib(11,kt)*DivMult(it)
          enddo
          enddo
          close(86)
          
          ! dump header row of tributary and diversion multiplier matrix files (Qmult, Qmult_div)
          read(88,*)
          do j=1,N
              READ(88,*) node,(Qmultdiv(node,jjn), jjn=1,Ndiv)   !dumps first column (which is compartment number)
          do jjn = 1,Ndiv
              if (isnan(Qmultdiv(node,jjn))) then
                  Qmultdiv(node,jjn) = 0.0
              endif
          enddo
          enddo
          close(88)
      endif


 3377 Format(I6,x,26(F7.1,x))


!>> Read in decay constants for water quality constituents
!---------------------------------
! ichem numbers of WQ constituents
!---------------------------------
!        1 = NO3 + NO2
!        2 = NH4
!        3 = DIN (1+2)
!        4 = Organic N
!        5 = TIP      !zw 4/28/2015 used to say TP
!        6 = TOC
!        7 = DO
!        8 = Live Algae
!      9 = Dead Algae = Detritus
!        10= DON
!        11= DOP
!        12= DIP   !Partition SRP = ParP*TP
!        13= ChLa  !Partition Chla= ParCla*LivA
!        14= POP !-EDW used to say 14=TKN !-EDW
!---------------------------------
      decay(:,:)=0  !zw added 04/06/2020
      if (iWQ >0) then
          do ichem = 1,14
             READ(85,*) (decay(ichem,i),i=1,14)       !zw 4/28/2015 decay unit=1/day
          enddo
          close(85)

          do ichem = 1,14
             do ichem2=1,14
                decay(ichem,ichem2) = decay(ichem,ichem2)/consd    !zw 4/28/2015 decay unit convert to 1/s
             enddo
          enddo
      endif

 1212 Format(13(F9.4,2x))

!>> Read in timeseries of water quality chemical loads (in mg/L)
!   do jkk=1,N
!!      do kd=1,365.25*Nyear        !-EDW replaced with loop over simdays
!          do kd=1,simdays
!!              do mk=1,13
!              do mk=1,14
!               QChem(jkk,mk,kd)=0.0
!           enddo
!       enddo
!   enddo
    !zw 4/28/2015 remove QChem initialization for all compartments, it is used to store the WQ loads from tributary only, see below

! Tributary Nutrients:  NO3, NH4, ON, TP

!>> Read and advance past the Tributary and Diversion Numbers from header row of input Water Quality files
      cChem(:,:,:)=0  !zw added 04/06/2020
      QChem(:,:,:)=0  !zw added 12/05/2023
      if (iWQ>0) then
          allocate(jqtrib(Ntrib))
          READ(80,*) (jqtrib(jt), jt=1,Ntrib) !changed in case trib number sequence different from TribQ
          READ(81,*)
          READ(82,*)
          READ(83,*)

!>> Skip Water Quality data that occurs in years prior to current model year
          do kt=1,startrun
              read(80,*)
              read(81,*)
              read(82,*)
              read(83,*)
              if (kt == startrun) then
                  write(1,*)'Nutrient loading input data starts at:'
                  write(1,66) '      line ',int(startrun+2)
                  write(*,*)'Nutrient loading input data starts at:'
                  write(*,66) '      line ',int(startrun+2)
              endif
          enddo

!>> Read Water Quality input data for current model year
          faN=0.1
          do kt=1,simdays
              FSEASON=1+0.5*cos(2*PI*(Float(kt)/365.25+0.8))
              FSEASON2=1.0+faN*cos(2*PI*(Float(kt)/365.25+0.25))     ! JAM Jan 09, 2011
!   Set loads to zero
!             do ichem = 1,14
!            do itrb =1,Ntrib
!                QChem(itrb,ichem,kt)=0.0
!            enddo
!            enddo

              READ(80,*) (cChem(jqtrib(jt),1,kt),jt=1,Ntrib)         !NO3 (mg/L) daily values    !JAM June 2008
              READ(81,*) (cChem(jqtrib(jt),2,kt),jt=1,Ntrib)         !NH4    ! JAM June 2008
              READ(82,*) (cChem(jqtrib(jt),4,kt),jt=1,Ntrib)         !ON     ! JAM June 2008
              READ(83,*) (cChem(jqtrib(jt),5,kt),jt=1,Ntrib)         !TP     ! JAM June 2008

              CtoN=5.7

              do jt=1,Ntrib
                  fctrib=0.0
!>> Add adjustment to Pearl River tributary flow !HARDCODED ADJUSTMENT
!                   if(jt == 5) then
                  if(jt == 19) then
                      fctrib=0.5
                  endif
!>> Convert Water Quality input data from concentration to loading rate (input data in mg/L, converted to g/sec here)
                  QChem(jtrib(jt),1,kt)=cChem(jtrib(jt),1,kt)
     &             *Qtrib(jtrib(jt),kt)*FSEASON                                      !NO3     ! JAM June 2008 & Jan 09, 2011
                  QChem(jtrib(jt),2,kt)=cChem(jtrib(jt),2,kt)
     &             *Qtrib(jtrib(jt),kt)*FSEASON                                      !NH4      ! JAM June 2008
                  QChem(jtrib(jt),3,kt)=QChem(jtrib(jt),1,kt)
     &             +QChem(jtrib(jt),2,kt)
                  QChem(jtrib(jt),4,kt)=cChem(jtrib(jt),4,kt)
     &             *Qtrib(jtrib(jt),kt)*FSEASON2 !*1.5 zw 4/28/2015 remove *1.5, adjust input file instead                 !ON                        ! JAM June 2008  *****increased to get TOC and ON to calibrate
                  QChem(jtrib(jt),6,kt)=QChem(jtrib(jt),4,kt)*5.7                        !*****increased to get TOC and ON to calibrate    JAM May 2011
                  QChem(jtrib(jt),5,kt)=cChem(jtrib(jt),5,kt)
     &             *Qtrib(jtrib(jt),kt)*FSEASON2*0.4                                           !TP                         ! JAM June 2008 !YW 0.4=1-PARDOP-PARPOP
                  QChem(jtrib(jt),12,kt)=cChem(jtrib(jt),5,kt)
     &             *Qtrib(jtrib(jt),kt)*ParP                                                                !TP                         ! JAM Dec 12, 2010
!                  QChem(jtrib(jt),10,kt)=0.0                                                                                                                           ! JAM March 5, 2011
                  QChem(jtrib(jt),11,kt)=cChem(jtrib(jt),5,kt)
     &             *Qtrib(jtrib(jt),kt)*ParDOP                                         !TP                         ! JAM Dec 12, 2010 !0.0      !JAM March 5 2011
!                  QChem(jtrib(jt),8,kt)=cChem(jtrib(jt),4,kt)
!     &           *Qtrib(jtrib(jt),kt)*FSEASON2
!     &           *CtoN*(1-fctrib)/2                                                                                                                                          ! JAM March 5, 2011
!                  QChem(jtrib(jt),9,kt)=cChem(jtrib(jt),4,kt)
!     &           *Qtrib(jtrib(jt),kt)*FSEASON2*CtoN
!     &           *(1+fctrib)/2.                                                                                                                                                     ! JAM March 5, 2011 /// April 8, 2011
                  QChem(jtrib(jt),8,kt)=cChem(jtrib(jt),4,kt)        !YW
     &              *Qtrib(jtrib(jt),kt)*FSEASON2/0.176/75.0/10.0  !YW rna*rca*calibration factor

                  QChem(jtrib(jt),9,kt)=cChem(jtrib(jt),4,kt)        !YW
     &             *Qtrib(jtrib(jt),kt)*FSEASON2/0.072/0.9        !YW        rnd*calibration factor

                  QChem(jtrib(jt),10,kt)=0.4*Qtrib(jtrib(jt),kt)*FSEASON2                              !YW average DON calculated from TOC = POC + DOC, POC = rca*ALG+rcd*DET, DOC = 5.69*DON, ALG and DET from table 77 in the 2012 report.
                  QChem(jtrib(jt),14,kt)=cChem(jtrib(jt),5,kt)
     &              *Qtrib(jtrib(jt),kt)*ParPOP                                          !TP                         ! JAM Dec 12, 2010      !JAM March 5 2011

              enddo
        !-EDW 6/4/2015 - Qchem units now in g/sec !!!changed in celldChem.f too!!!
          ! zw 4/28/2015 QChem unit=kg/day in here

    !zw 4/28/2015 remove the following QChem adjustments for high flows, should adjust the tributary WQ inputs instead during calibration
!**** seem to be underestimating nitrogen load for high flows, therefore mult by 2 if greater than average LB increased JAM April 4 08

!       If(QChem(1,1,kt).ge.995.6982)QChem(1,1,kt)=QChem(1,1,kt)
!     &         *2.02502                                                    ! increased to match observations
!       If(QChem(2,1,kt).ge.1314.888)QChem(2,1,kt)=QChem(2,1,kt)
!     &         *1.202                                                      ! increased to match observations
!       If(QChem(5,1,kt).ge.849.0132)QChem(5,1,kt)=QChem(5,1,kt)
!     &         *1.022
!       If(QChem(6,1,kt).ge.5.66589)QChem(6,1,kt)=QChem(6,1,kt)*2.0502      ! increased to match observations QChem(6,1,kt)*2.502
!       If(QChem(7,1,kt).ge.197.6641)QChem(7,1,kt)=QChem(7,1,kt)*1.02
!       If(QChem(9,1,kt).ge.1275.069)QChem(9,1,kt)=QChem(9,1,kt)*1.102      ! increased to match observations
!
!       QChem(8,5,kt) = QChem(8,5,kt)*1.05
!       QChem(8,4,kt) = QChem(8,4,kt)*1.82          !2.05
!       QChem(8,1,kt) = QChem(8,1,kt)*1.82          !reduced 2-->1.8        ! JAM Feb 27, 2011
!       QChem(8,2,kt) = QChem(8,2,kt)*1.82
          enddo
          close(80)
          close(81)
          close(82)
          close(83)
          deallocate(jqtrib)
      endif
!>> Read in anthropogenic N-NO3 loads to compartments (in kg/d) - loads from farms and MWWTPs
      AnthL(:)=0  !zw added 04/06/2020
      if (iWQ>0) then
          do j=1,N
              READ(44,*)AnthL(j)
          enddo
          close(44)
      endif

!>> Initialize Atmospheric Load arrays to zero
!      do kt=1,simdays
!
!       do ichem=1,14
!           Qatm(1,ichem,kt)=0.0
!       enddo
!     enddo
      QAtm(:,:,:)=0
      if (iWQ>0) then
!>> Skip header row of Atmospheric Load file
          READ(84,*)

!>> Skip Atmospheric Load data that occurs in years prior to current model year
          do kt=1,startrun
              read(84,*)
              if (kt == startrun) then
                  write(1,*)'Atmospheric loading input data starts at:'
                  write(1,66) '      line ',int(startrun+2)
                  write(*,*)'Atmospheric loading input data starts at:'
                  write(*,66) '      line ',int(startrun+2)
              endif
          enddo

!>> Read Atmospheric Load data that occurs in current model year
          do kt=1,simdays
              READ(84,*)QAtm(1,5,kt),QAtm(1,1,kt),QAtm(1,2,kt),QAtm(1,4,kt) !TP, NO3  NH4  ON (kg/km2/day)      !JAM April 17, 2011
              QAtm(1,3,kt) = QAtm(1,1,kt)+QAtm(1,2,kt)
              !QAtm(1,4,kt) = QAtm(1,4,kt)*1.5  !zw 4/28/2015 remove, adjust the input file instead ! JAM TKN & TOC too low May 2011
              QAtm(1,6,kt)=QAtm(1,4,kt)/3.                          ! JAM April 17, 2011
              QAtm(1,7,kt)=0.1
              QAtm(1,8,kt) = 0.
              QAtm(1,9,kt)=QAtm(1,4,kt)*2/3.                            ! JAM April 17, 2011
              QAtm(1,10,kt) =0.1*QAtm(1,4,kt)                           ! JAM April 17, 2011
              QAtm(1,11,kt)=0.0
              QAtm(1,12,kt) = 0.5*QAtm(1,5,kt)                      ! JAM April 17, 2011
              QAtm(1,13,kt)=0.
              QAtm(1,14,kt)=0.5*QAtm(1,5,kt)                            ! JAM April 17, 2011
          enddo
 
!>>  Calculate DIN values from other WQ input data
      !zw 4/28/2015 remove QAtm ajustments by Faatm, should adjust the input file instead
      !Faatm=0.91
    !if(j>10) then
      !    Faatm=0.75    !zw 4/28/2015 j=N so Faatm=0.75 always
      !endif

      !do kt=1,simdays
!     ! do kt=1,365.25*Nyear                       !-EDW replaced with loop over simdays
    !   QAtm(1,1,kt) = QAtm(1,1,kt)*Faatm
    !   QAtm(1,2,kt) = QAtm(1,2,kt)*Faatm
    !   QAtm(1,3,kt) = QAtm(1,1,kt) + QAtm(1,2,kt)
!       do jt=1,Ntrib    !zw 4/28/2015 remove do loop for DIN, already done in the tributary load inputs section above
!           QChem(jtrib(jt),3,kt)=(QChem(jtrib(jt),1,kt)
!     &             +QChem(jtrib(jt),2,kt))
!       enddo
    !enddo
          close(84)
      endif

!zw 4/28/2015 UplandNP never used during WQ calculations.
!Read in Upland nutrient loads (constants for each compartment, not timeseries)
!      do j=1,N
!       READ(118,*) UplandNP(j,1),UplandNP(j,2),UplandNP(j,4),UplandNP(j,5)
!
! Calculate various upland nutrient loads for each compartment based on imported upland data (constants ofr each compartment, not timeseries)
!          UplandNP(j,4)=UplandNP(j,4)*1.5
!       UplandNP(j,3)=UplandNP(j,2)+UplandNP(j,1)
!       UplandNP(j,8)=UplandNP(j,4)*0.25*5.7            ! JAM April 17, 2011
!       UplandNP(j,9)=UplandNP(j,4)*0.75*5.7
!       UplandNP(j,12)=UplandNP(j,5)*0.5
!       UplandNP(j,13)=UplandNP(j,8)*0.5*.03            ! JAM April 17, 2011
!       UplandNP(j,14)=UplandNP(j,5)*0.5                ! JAM April 17, 2011
!   enddo

      cChemdiv(:,:,:) = 0.0
      QChemdiv(:,:,:) = 0.0
      if (iWQ>0) then
!>> Skip header row of Diversion Nutrient Load file
          read(87,*)

!>> Skip Diversion Nutrient Load data that occurs in years prior to current model year
          do kt=1,startrun
              read(87,*)
              if (kt == startrun) then
                  write(1,*)'Diversion nutrient input data starts at:'
                  write(1,66) '      line ',int(startrun+2)
                  write(*,*)'Diversion nutrient input data starts at:'
                  write(*,66) '      line ',int(startrun+2)
              endif
          enddo

!>> Read Diversion Nutrient Load data for current model year
          do kt=1,simdays
!>> initialize CchemDiv values to 0.0 before reading in input data
!             do ktk = 1,14
!                 cChemdiv(1,ktk,kt) = 0.0
!             enddo

              READ(87,*) cChemdiv(1,5,kt),cChemdiv(1,1,kt),
     &          cChemdiv(1,2,kt), cChemdiv(1,4,kt)      !zw 4/28/2015 unit=mg/L

              cChemdiv(1,8,kt)=cChemdiv(1,4,kt)*5.7/2.        ! JAM March 2011
              cChemdiv(1,9,kt)=cChemdiv(1,4,kt)*5.7/2.

!>> Calculate DIN values for Diversions from other WQ input data
              cChemdiv(1,3,kt) = cChemdiv(1,1,kt)+cChemdiv(1,2,kt)
!>> Convert loads to kg/day
              do it=1,Ndiv
                  do mk=1,14                                    ! check limit on mk JAM Aug 10, 2009 increased 7-->9?--11-->13
! add in reductions for low flow from Lane et. al
                      QChemdiv(it,mk,kt)=(cChemdiv(1,mk,kt))*
     &                  (Qdiv(it,kt))                   !JAM Jan 2, 2011
                  enddo
                  QChemdiv(it,12,kt)=QChemdiv(1,4,kt)*ParPMR     !zw 4/28/2015 move it out of the above do loop
              enddo
          enddo
      !-EDW 6/4/2015 - Qchemdiv units now in g/sec !!!changed in celldChem.f too!!!
    !zw 4/28/2015 QChemdiv unit is kg/day in here

!>> Input algae conentrations for each tributary (assume 1 mg/L/trib) and convert to a load
      !zw 4/28/2015 remove the following do loop for tributary ALG and DET loads, they have already ben calculated in the tributary WQ load section
      !do kt=1,simdays
!     ! do kt=1,365.25*Nyear                            !-EDW replaced with loop over simdays
    !   do jt=1,Ntrib
    !       cChem(jtrib(jt),8,kt)=cChem(jtrib(jt),4,kt)*5.7
    !       cChem(jtrib(jt),9,kt)=cChem(jtrib(jt),4,kt)*5.7     !JAM Dec 12, 2010 increased detritus 0.01-->1
    !       QChem(jtrib(jt),8,kt)=(cChem(jtrib(jt),8,kt))
      !&                *(Qtrib(jtrib(jt),kt)*conv)
    !       QChem(jtrib(jt),9,kt)=(cChem(jtrib(jt),9,kt))
      !&                *(Qtrib(jtrib(jt),kt)*conv)
    !   enddo
    !enddo
          close(87)
      endif
!>> Change units of WQ chemicals from kg/day to kg/s  ===/consd=24*3600
      do mk = 1,numChem
          do kt=1,simdays
!         do kt=1,365.25*Nyear                            !-EDW replaced with loop over simdays
            QAtm(1,mk,kt) = QAtm(1,mk,kt)/consd
!           do jt =1,Ntrib
!               QChem(jtrib(jt),mk,kt)=QChem(jtrib(jt),mk,kt)/consd
!           enddo
!           do it=1,Ndiv
!               QChemdiv(it,mk,kt)=QChemdiv(it,mk,kt)/consd     ! JAM Feb 2010
!           enddo
          enddo
      enddo

!>> Calculate Carbon loads in tributaries and diversions
    !zw 4/28/2015 remove the following DO loop,these are all hardcoded ajustments for WQ inputs, should adjust the input files instead
      !do kt=1,simdays
!     ! do kt=1,365.25*Nyear                            !-EDW replaced with loop over simdays
    !   do jt=1,Ntrib
    !       cChem(jtrib(jt),6,kt)=MIN((cChem(jtrib(jt),3,kt))
      !&                *5.68,(cChem(jtrib(jt),5,kt)*41.11))
    !       QChem(jtrib(jt),6,kt)=cChem(jtrib(jt),6,kt)
      !&                *(Qtrib(jtrib(jt),kt)*.001)                     !!?? JAM Feb 2011
    !   enddo
    !
    !   do it=1,Ndiv
    !       if(Qdiv(jdiv(it),kt).le.500.) then
    !
      !            cChemdiv(1,6,kt)=MIN((cChemdiv(1,1,kt)*0.65+
      !&                cChemdiv(1,2,kt)*0.72)*5.68,
      !&                (cChemdiv(1,5,kt)*0.75)*41.11)
    !
      !            cChemdiv(1,1,kt)=    cChemdiv(1,6,kt)*0.5
    !
      !            cChemdiv(1,4,kt)=    cChemdiv(1,6,kt)*0.5
    !
      !            QChemdiv(jdiv(it),6,kt)=cChemdiv(1,6,kt)
      !&                *(Qdiv(jdiv(it),kt)*.001)
    !
      !            QChemdiv(jdiv(it),1,kt)=cChemdiv(1,1,kt)
      !&                *(Qdiv(jdiv(it),kt)*.001)
    !
      !            QChemdiv(jdiv(it),4,kt)=cChemdiv(1,4,kt)
      !&                *(Qdiv(jdiv(it),kt)*.001)
      !        else
      !            cChemdiv(1,6,kt)=MIN((cChemdiv(1,3,kt))*5.68,
      !&                    (cChemdiv(1,5,kt)*41.11))
    !
      !            QChemdiv(jdiv(it),6,kt)=cChemdiv(1,6,kt)
      !&                    *(Qdiv(jdiv(it),kt)*.001)
    !       endif
    !   enddo
    !enddo

!>> Initialize Sediment arrays to zero
!   Qss(:,:)=0.0
      QssT(:,:)=0.0
      QssTdiv(:,:)=0.0
      ASandT(:,:)=0.0
      ASandD(:,:)=0.0
      ACCSEDj(:)=0.0

  !MP2023 added zw 04/06/2020
      cssT(:,:,:)=0
      cssFines(:)=0

!>> Read Tributary Numbers from header row of Tributary sand concentration input file
      allocate(jqtrib(Ntrib))
      READ(55,*) (jqtrib(jt), jt=1,Ntrib)  !changed in case trib number sequence different from TribQ
!>> Skip Tributary Sediment data for years prior to current model year
      do kt=1,startrun
          read(55,*)
          if (kt == startrun) then
              write(1,*)'Tributary sand concentration input data'
              write(1,66) ' starts at line : ',int(startrun+2)
              write(*,*)'Tributary sand concentration input data'
              write(1,66)' starts at line : ',int(startrun+2)
          endif
      enddo
!>> Read Tributary Sediment data for current model year (in mg/L)
      do kt=1,simdays
          READ(55,*) (cssT(jqtrib(jt),kt,1), jt=1,Ntrib)    !READ in tributary sediment conc cssT   !JAM correction April 8, 2007
      enddo
      close(55)

      read(555,*)(jqtrib(jt), jt=1,Ntrib)
      do kt=1,startrun
          read(555,*)
          if (kt == startrun) then
              write(1,*)'Tributary fine concentration input data'
              write(1,66)' starts at line : ',int(startrun+2)
              write(*,*)'Tributary fine concentration input data'
              write(*,66)' starts at line : ',int(startrun+2)
          endif
      enddo
!>> Read Tributary Sediment data for current model year (in mg/L)
      do kt=1,simdays
      !>> Read each row of tributary fines data, then divide by three to partition into silt, clay, and floc
          READ(555,*) (cssFines(jqtrib(jt)), jt=1,Ntrib)    !READ in tributary sediment conc cssT   !JAM correction April 8, 2007
          do jt = 1, Ntrib
              cssT(jt,kt,2) = cssFines(jt)/3.0
              cssT(jt,kt,3) = cssFines(jt)/3.0
              cssT(jt,kt,4) = cssFines(jt)/3.0
          enddo
      enddo
      close(555)
      deallocate(jqtrib)


!>> Create tributary sediment concentrations and convert loads to kg/d
      Parsand=0.05

!      do kt=1,simdays
!       do jt=1,Ntrib
!           QssT(jt,kt)=cssT(jt,kt)*(QTrib(jt,kt)*conv)                 ! kg/d
!              ACCSED(jtrib(jt))=QssT(jtrib(jt),kt)+ACCSED(jtrib(jt))       !*24.*3600  (kg/d)  May 17, 2011 JAM ACCSED(jtrib(jt))=QssT(jtrib(jt),kt)*24.*3600+ACCSED(jtrib(jt))
!       enddo
!
!       do j=1,N
!           do jt=1,Ntrib
!               ASANDT(j,kt)=ParSand*2.*cssT(jt,kt)*QTrib(jt,kt)
!     &                     *QMult(j,jt)+ASANDT(j,kt)
!           enddo
!       enddo
!   enddo

!      do j=1,N
!       do jt=1,Ntrib
!           ACCSEDj(j)=ACCSED(jt)*Qmult(j,jt)
!       enddo
!   enddo

!>> Read Sediment-water-ratios for diversions
      write(1,*) 'Reading in sediment-water-ratios for diversions.'
      write(*,*) 'Reading in sediment-water-ratios for diversions.'
      SWRsand(:)=0  !zw added 04/06/2020
      SWRfines(:)=0
      ParSandD(:)=0
      CSSTdiv(:,:,:)=0
      if (Ndiv>0) then
          read(89,*) !skip header row
          do it = 1,Ndiv
              read(89,*) dump_int,SWRsand(it),SWRfines(it)
          enddo
          close(89)

          write(1,*) 'Calculating sediment concentration in diversions
     & based on Mississippi River sediment downstream of Bonnet Carre.'

          write(*,*) 'Calculating sediment concentration in diversions
     & based on Mississippi River sediment downstream of Bonnet Carre.'

          Qmax=35400.   !flowrate allowed past Bonnet Carre (m3/s)
          do kt=1,simdays
          do it = 1,Ndiv
!>> Sand partition is function of Mississippi River flow d/s Bonnet Carre --->Sand Partition factor varies with Q^2
              ParSandD(kt)=SWRsand(it)*6.*(Qtrib(11,kt)/Qmax)
     &              *(Qtrib(11,kt)/Qmax)
              CSSTdiv(it,kt,1) = ParSandD(kt)*CSST(11,kt,1) !tributary 11 is Mississippi River CSS

!>> Split Miss Riv fines evenly between silt, clay and floc
              CSSTdiv(it,kt,2) = SWRfines(it)*CSST(11,kt,2)
              CSSTdiv(it,kt,3) = SWRfines(it)*CSST(11,kt,3)
              CSSTdiv(it,kt,4) = SWRfines(it)*CSST(11,kt,4)
          enddo
          enddo
      endif
!>> Calculate Diversion sediment load (kg/day)


!      do kt=1,simdays
!       do it=1,Ndiv
! Apply equations to  to diversion sediment loads to account for increased sand percentage in diversions
!           ParSandD(kt)=ParSand*6.*(Qdiv(23,kt)/Qmax)   !Qdiv is flow d/s Bonnet Carre --->Sand Partition factor varies with Q^2
!     &             *(Qdiv(23,kt)/Qmax)             ! near outfall sand and coarse material sedimentation
!           if (Qdiv(it,kt).le.10.) then            ! removal near diversion oufall
!               QssTdiv(it,kt)=(cssTdiv(it,kt)*0.97)*
!     &                 (Qdiv(it,kt)*conv)
!              else
!               QssTdiv(it,kt)=cssTdiv(it,kt)*
!     &                 (Qdiv(jdiv(it),kt)*conv)
!                  ACCSED(jdiv(it))=QssTdiv(jdiv(it),kt)
!     &                     +ACCSED(jdiv(it))           ! (kg/d)  *24.*3600  JAM June 21, 2011
!           endif
!       enddo
!   enddo

!      do kt=1,simdays
!       do j=1,N
!           do it=1,Ndiv
!               ACCSEDj(j)=ACCSED(it)*Qmultdiv(j,it)+ ACCSEDj(j)        !kg accum
!               ASANDD(j,kt)=ParSandD(kt)*Qdiv(it,kt)*cssTdiv(it,kt)        !ParSandD fraction of sand in total River sediment load
!     &                *SWR(it)                                             !Sediment water Ratio for sand capture by the outlet; FWOP SWR = 1.0 AMc Jan 8 2014
!     &                *Qmultdiv(j,it)+(1.-ParSandD(kt))*Qdiv(it,kt)            !trap 10% of fines with sand JAM June 25, 2001
!     &                *cssTdiv(it,kt)*Qmultdiv(j,it)*0.10+ASANDD(j,kt)     !g/s QssTdiv(it,kt) Assumes 10% of fines are rapidly settleable silt
!           enddo
!       enddo
!   enddo

! Convert diversion sediment loads to kg/s
!      do kt=1,simdays
!       do it=1,Ndiv
!           QssTdiv(jdiv(it),kt)=QssTdiv(jdiv(it),kt)/consd
!       enddo
!   enddo

! Convert tributary sediment loads to kg/s
!      do kt=1,simdays
!       do jt=1,NTrib
!           QssT(jt,kt)=QssT(jt,kt)/consd           ! (kg/s)        ! JAM May 17, 2011
!       enddo
!   enddo

!>> Skip header row of Precip input file
      Rain(:,:)=0  !zw added 04/06/2020
      read(42,*)                                        !dump header

!>> Skip Precip data for years prior to current model year
      do kt=1,startrun
          read(42,*)
          if (kt == startrun) then
              write(1,*)'Precip input data starts at:'
              write(1,66) '      line ',int(startrun+2)
              write(*,*)'Precip input data starts at:'
              write(*,66) '      line ',int(startrun+2)
          endif
      enddo

!>> Read Precip data for current model year
      do kt=1,simdays
          READ(42,*) dump_int,Rain(kt,1:raingages)
      enddo
      close(42)

!>> Skip header row of ET input file
      PET(:,:)=0  !zw added 04/06/2020
      read(40,*)                                          !dump header

!>> Skip ET data for years prior to current model year
      do kt=1,startrun
          read(40,*)
          if (kt == startrun) then
              write(1,*)'ET input data starts at:'
              write(1,66) '      line ',int(startrun+2)
              write(*,*)'ET input data starts at:'
              write(*,66) '      line ',int(startrun+2)
          endif
      enddo

!>> Read ET data for current model year
      do kt=1,simdays
          READ(40,*) dump_int,PET(kt,1:etgages)
      enddo
      close(40)

!>> Filter for negative precip and ET input data and set to zero

      write(1,*)'Checking Precip and ET input data for negative values.'
      write(*,*)'Checking Precip and ET input data for negative values.'
      ETzero=0  !zw added 04/06/2020
      Rzero=0
      do kt = 1,simdays
          do etg = 1,etgages
              if (PET(kt,etg)<0.0) then
                  PET(kt,etg) = 0.0
                  ETzero = ETzero + 1
              endif
          enddo

          do rg = 1,raingages
              if (Rain(kt,rg)<0.0) then
                  Rain(kt,rg) = 0.0
                  Rzero = Rzero + 1
              endif
          enddo
      enddo

      if(Rzero > 0) then
          write(1,2) Rzero, 'precip records were negative and set to 0.'
          write(*,2) Rzero, 'precip records were negative and set to 0.'
      endif

      if(ETzero > 0) then
          write(1,2) ETzero,'ET records were negative and set to 0.'
          write(*,2) ETzero,'ET records were negative and set to 0.'
      endif

!>> Calculate average annual ET (mm/day) at ET gage
      ! this is only used if fpet flag is set to 0,
      do etg=1,etgages
          ETA(etg) = 0.0
          do kt=1,simdays
              ETA(etg)=ETA(etg)+PET(kt,etg)
          enddo
          ETA(etg)=ETA(etg)/simdays
      enddo

!>> Skip header row of wind vector input files (both X & Y vectors)
      windx_data(:,:)=0.0  !zw added 04/06/2020
      windy_data(:,:)=0.0
      read(43,*)
      read(46,*)                                          !dump header

!>> Skip wind data for years prior to current model year
      windstartrun=startrun*24/dtwind

      do kt=1,windstartrun
          read(43,*)
          read(46,*)
          if (kt == windstartrun) then
              write(1,*)'Wind input data starts at:'
              write(1,66) '      line ',int(windstartrun+2) !added zw 2/9/2015
              write(*,*)'Wind input data starts at:'
              write(*,66) '      line ',int(windstartrun+2) !added zw 2/9/2015
          endif
      enddo

!>> Read wind data for current model year
      do kt=1,simdays*24/dtwind
          READ(43,*) dump_int,dump_text,windx_data(kt,1:windgages)
          READ(46,*) dump_int,dump_text,windy_data(kt,1:windgages)
      enddo
      close(43)
      close(46)

!>> Read Boundary Conditions file
      KBC(:)=0
      Read(125,*)(KBC(jj), jj=1,mds) !AMc Oct 8 2013
      close(125)
    !zw added 04/07/2020 to determine whether a cell is offbc or not
      flag_offbc(:)=0
      do i=1,mds
          jj=KBC(i)
          flag_offbc(jj)=1
      enddo    

!>> Read in data to transpose near-shore observed water level timeseries to off-shore water levels
      transposed_tide(:,:)=0  !zw added 04/06/2020
      read(48,*)
      do i = 1,tidegages
          read(48,*)dump_int,transposed_tide(dump_int,1),transposed_tide(dump_int,2)
      enddo
      close(48)

!>> Skip header row of Tide Gage and Surge data input files
      TideData(:,:)=0.0  !zw added 04/06/2020
      Surge(:,:)=0.0
      read(47,*)
      read(110,*)

!>> Skip Tide Gage data for years prior to current model year
      tidestartrun = startrun*24/dttide  !added zw 2/9/2015
      do kt=1,tidestartrun
          read(47,*)
          read(110,*)

          if (kt == tidestartrun) then  !added zw 2/9/2015
              write(1,*)'Tide gage and surge input data start at:'
              write(1,66) '      line ',int(tidestartrun+2)   !added zw 2/9/2015
              write(*,*)'Tide gage and surge input data start at:'
              write(*,66) '      line ',int(tidestartrun+2)   !added zw 2/9/2015
          endif
      enddo

!>> Read Tide Gage and Surge data for current model year
      do kt=1,(simdays*24/dttide+1)  !zw modififed 04/06/2020   !YW! +1 to include the final row
          read(47,*)dump_int,dump_text,TideData(kt,1:tidegages)
          read(110,*)dump_int,dump_text,Surge(kt,1:mds)
      enddo
      close(47)
      close(110)

!>> Read in data to weight (by distance) the nearest observed water level timeseries to off-shore boundary compartments that do not have an observed WSEL timeseries
      weighted_tide(:,:)=0  !zw added 04/06/2020
      read(49,*)
      do i = 1,(mds-tidegages)
          read(49,*) weighted_tide(i,1), ! compartment ICM ID number for BC comparments WITHOUT observed water level timeseries
     &            weighted_tide(i,2),    ! distance weighting factor for nearest BC comparment WITH observed water level timeseries to the East
     &            weighted_tide(i,3),    ! nearest BC compartment WITH observed water level timeseries to the East
     &            weighted_tide(i,4),    ! distance weighting factor for nearest BC comparment WITH observed water level timeseries to the West
     &            weighted_tide(i,5)     ! nearest BC comparment WITH observed water level timeseries to the West
      enddo
      close(49)

!

!>> Skip header row of Meteorology and Mississippi River Temperature input files
      ta(:)=0 !zw added 04/06/2020
      tw(:)=0
      TempMR(:)=0
      read(45,*)

!>> Skip Meteorology and Mississippi River Temperature data for years prior to current model year
      do kt=1,startrun
          read(45,*)
          read(74,*)
          if (kt == startrun) then
              write(1,*)'Met. tmp input data starts at:'
              write(1,66) '      line ',int(startrun+2)
              write(*,*)'Met. tmp input data starts at:'
              write(*,66) '      line ',int(startrun+2)
              write(1,*)'River tmp input data starts at:'
              write(1,66) '      line ',int(startrun+1)
              write(*,*)'River tmp input data starts at:'
              write(*,66) '      line ',int(startrun+1)
          endif
      enddo

!>> Read Meteorology (temp) and Mississippi River Temperature data for current model year
      do kt=1,simdays
          read(45,*) ta(kt),tw(kt)
          read(74,*) TempMR(kt)
      enddo
      close(45)
      close(74)

      do kt=1,simdays
          FSEASON3=0.1-2.0*cos(2*PI*(Float(kt)/365.25-0.05))            ! JAM April 2011 Jan 09 11 **March 2011**April 2011

          do j=1,N
              tadd=0.0                                                ! *3.0*(ta(kt)-20.)/10.
              Tempair(j,kt)=ta(kt)
              TMtrib(j,kt)=tw(kt)
              dlow=0.25
              if(Tempair(j,kt).lt.4.)dlow=(1.-Tempair(j,kt)/4.)       ! April 2011 JAM
              T7(kt)=(Tempair(j,kt)+Tempair(j,max(1,kt-1))+
     &          Tempair(j,max(1,kt-2))+Tempair(j,max(1,kt-3))
     &          +Tempair(j,max(1,kt-4))
     &          +Tempair(j,max(1,kt-5))
     &          +Tempair(j,max(1,kt-6)))/7.
              Tempe(j,kt)=0.98*(T7(kt)-20.)+21.5+0.5*2.+dlow
     &          +FSEASON3+Tadd                              ! April 2011 JAM Oct 2010/March 2011/April 2, 2011   0.99*T7(kt)+1.48+dlow+FSEASON3
          enddo

          ta_k(kt) = tw(kt) + 273.15

!***********Calculate Oxygen Saturation from QUAL2K eq. 140 - ALPHA 1992
          lnO2Sat(kt) = 0.0
          lnO2Sat(kt) = -139.34411+157570.1/ta_k(kt)-66423080./
     &          (ta_k(kt)**2)+12438000000./(ta_k(kt)**3)-862194900000.
     &          /(ta_k(kt)**4)
          a = lnO2Sat(kt)
          O2Sat(kt) = exp(a)
          WRITE(93,*) O2Sat(kt)

      enddo

!!******************* Initial Conditions
!      FSEASON3=0.1-2.0*cos(2*PI*(1.0/365.25-0.05))  !added zw 04/06/2020
!     do j=1,N
!         Tempw(j,1)=0.98*(Tempair(j,1)-20.)+21.5+0.01
!     &             *2.+dlow/2.+FSEASON3            ! JAM Oct 2010/March 2011 0.9855*Tempair(j,1)+1.38; this is replaced by initial condition file
!         Sacc(j,1)=0.0                             !in main.f
!         Sacch_edge(j,1)=0.0                       !in main.f
!         Sacch_int(j,1)=0.0                        !in main.f
!     enddo

  !>>boundary conditions data for salinity and WQ
      SBC(:)=0  !zw added 04/06/2020
      BCTSS(:)=0
      BCNO3(:)=0
      BCNH4(:)=0
      BCON(:)=0
      BCTP(:)=0
      BCDO(:)=0
      BCTOC(:)=0
      BCLA(:)=0
      BCDA(:)=0
      BCage(:)=0


      do i=1,mds   !AMc Oct 8 2013
!         jj=kbc(i)   !AMc Oct 8 2013
          READ(56,*) jmds,SBC(jmds),BCTSS(jmds),BCNO3(jmds),BCNH4(jmds),
     &          BCON(jmds),BCTP(jmds),BCDO(jmds),BCTOC(jmds),BCLA(jmds),
     &          BCDA(jmds),BCage(jmds)                                  !added age

      enddo
      close(56)



!>> Skip Boundary Condition Temperature data for years prior to current model year (no header row to skip)
      TempwBC(:,:)=0  !zw added 04/06/2020
      do kt=1,startrun
          read(101,*)
          if (kt == startrun) then
              write(1,*)'Bndry. cond. tmp. input data starts at:'
              write(1,66) '      line ',int(startrun+1)
              write(*,*)'Bndry. cond. tmp. input data starts at:'
              write(*,66) '      line ',int(startrun+2)
          endif
      enddo

!>> Read Bounday Condition Temperature data for current model year
      do kt=1,simdays
          READ(101,*)(TempwBC(jj,kt), jj=1,mds)

!>> If compartment has a boundary condition, replace temperature values with downstream boundary conditions that were just read in for BC locations
          do jkk=1,mds   !AMc Oct 8 2013
              jj=kbc(jkk)   !AMc Oct 8 2013
              Tempe(jj,kt)=TempwBC(jkk,kt)

          enddo
      enddo
      close(101)

!>> Read input link file to apply flow limiter  !YW
      linkslimiter(:)=0
      flag_apply(:)=0
      if (nlinklimiter>0) then
          read(500,*)
          do i = 1,nlinklimiter
              read(500,*) linkslimiter(i)
              flag_apply(linkslimiter(i))=1
          enddo
          close(500)
      endif
      
!>> Generate Compartment/Link Connectivity Table
      write(1,*)
      write(1,*)'-----------------------------------------'
      write(1,*)'Generating link connectivity table.'
      write(1,*)'-----------------------------------------'
      write(1,*)
      write(*,*)
      write(*,*)'-----------------------------------------'
      write(*,*)'Generating link connectivity table.'
      write(*,*)'-----------------------------------------'
      write(*,*)

      do  j=1,N
          nlink2cell(j) = 0
          do k=1,maxconnect
              icc(j,k) = 0
          enddo
      enddo

      do j=1,N
          k=0
          do i=1,M
              if(jus(i).eq.j)then
                  k=k+1
                  icc(j,k) = i
              endif
              if(jds(i).eq.j) then
                  k=k+1
                  icc(j,k)=-i
              endif
          enddo
          nlink2cell(j) = k   !> set number of links in compartment value in array

          if (k > maxconnect) then
              write(1,1121)'Number of links in compartment',j,
     &        ' exceeds maxconnect value in RunControlR.dat.'
              write(1,*)'***** Please increase maxconnect in RunControlR.dat!!!****'

              write(*,1121)' ****** Number of links in compartment ',j,
     &        ' exceeds maxconnect value in RunControlR.dat.'
              write(*,*)'***** Please increase maxconnect in RunControlR.dat!!!****'
              stop
          endif

1121  format(7x,A,x,I0,x,A,x,I0,x,A)
1122  format(7x,A,x,I0,x,A,x,I0,x,A)
      enddo

      do j=1,N                  !chg JAM Aug 10, 2009  sicc is the sign of the link k connecting node j
          do k=1,maxconnect
              if(icc(j,k).ne.0) then
                  if(icc(j,k).lt.0.0)sicc(j,k)=-1.0
                  if(icc(j,k).gt.0.0)sicc(j,k)=1.0
              else
                  sicc(j,k)=0.0
              endif
          enddo
      enddo

!>> Lookup each corresponding hydro compartment and links for interpolating data to each 500 m grid cell
!>> Input grid lookup tables MUST have 22 columns - Col1= gridID, Col2=compartmentID,Col3-Col22=links with -9999 for no connections
      read(200,*) dump_text !dump header row
      do jk=1,n_500m_cells    !j = number of grid cells in lookup table + 1 (for header row)
          read(200,*) grid_lookup_500m(jk,1),     ! 500 m grid ID
     &                    grid_lookup_500m(jk,2),     ! hydro compartment associated with 500 m grid cell
     &                    grid_lookup_500m(jk,3),     ! hydro link #1 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,4),     ! hydro link #2 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,5),     ! hydro link #3 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,6),     ! hydro link #4 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,7),     ! hydro link #5 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,8),     ! hydro link #6 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,9),     ! hydro link #7 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,10),    ! hydro link #8 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,11),    ! hydro link #9 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,12),    ! hydro link #10 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,13),    ! hydro link #11 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,14),    ! hydro link #12 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,15),    ! hydro link #13 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,16),    ! hydro link #14 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,17),    ! hydro link #15 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,18),    ! hydro link #16 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,19),    ! hydro link #17 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,20),    ! hydro link #18 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,21),    ! hydro link #19 associated with 500 m grid cell
     &                    grid_lookup_500m(jk,22)     ! hydro link #20 associated with 500 m grid cell
      enddo
      close(200)

!>> Lookup distances to centroids for each corresponding hydro compartment and link for interpolating data to each 500 m grid cell
!>> Input grid lookup tables MUST have 22 columns - Col1= gridID, Col2=compartmentID,Col3-Col22=links with -9999 for no connections

      read(201,*) dump_text !dump header row
      do jk=1,n_500m_cells    !j = number of grid cells in lookup table + 1 (for header row)
          read(201,*) grid_interp_dist_500m(jk,1),     ! 500 m grid ID
     &                    grid_interp_dist_500m(jk,2),     ! distance to centroid of hydro compartment associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,3),     ! distance to centroid of hydro link #1 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,4),     ! distance to centroid of hydro link #2 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,5),     ! distance to centroid of hydro link #3 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,6),     ! distance to centroid of hydro link #4 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,7),     ! distance to centroid of hydro link #5 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,8),     ! distance to centroid of hydro link #6 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,9),     ! distance to centroid of hydro link #7 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,10),    ! distance to centroid of hydro link #8 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,11),    ! distance to centroid of hydro link #9 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,12),    ! distance to centroid of hydro link #10 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,13),    ! distance to centroid of hydro link #11 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,14),    ! distance to centroid of hydro link #12 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,15),    ! distance to centroid of hydro link #13 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,16),    ! distance to centroid of hydro link #14 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,17),    ! distance to centroid of hydro link #15 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,18),    ! distance to centroid of hydro link #16 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,19),    ! distance to centroid of hydro link #17 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,20),    ! distance to centroid of hydro link #18 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,21),    ! distance to centroid of hydro link #19 associated with 500 m grid cell
     &                    grid_interp_dist_500m(jk,22)     ! distance to centroid of hydro link #20 associated with 500 m grid cell
      enddo
      close(201)

!>> Lookup elevation from topo/bathy data for each 500 m grid cell (updated in the Wetland Morph ICM routine)
      read(202,*) dump_text !dump header row
      do jk=1,n_500m_cells    !j = number of grid cells in lookup table + 1 (for header row)
          read(202,*)dump_int,bed_elev_500m(jk),land_elev_500m(jk),
     &                per_land_500m(jk),dump_float,per_water_500m(jk)
!>> if grid cell bed elevation is missing (flagged as -9999) then set it equal to the compartment's bed elevation
          if (bed_elev_500m(jk) < -9990) then
              bed_elev_500m(jk) = bed(grid_lookup_500m(jk,2))
              gridbedcount = gridbedcount + 1
          endif
!>> if grid cell land elevation is missing (flagged as -9999) then set it equal to the compartment's marsh elevation
          if (land_elev_500m(jk) < -9990) then
              land_elev_500m(jk) = BedM(grid_lookup_500m(jk,2))
              gridmarcount = gridmarcount + 1
          endif
!>> if grid cell is missing percent land value - assume that it is 100% land
          if (per_land_500m(jk) < -9990) then
              per_land_500m(jk) = 100.0
              gridplcount = gridplcount + 1
!>> If a grid cell is not 100% land, but it does not have any water associated with it, then reclassify thegrid as being 100% land
!>> This is necessary due to mismatch of resolution in the 30-m land/water data and the 500-m veg grid
          elseif(per_land_500m(jk) < 100) then
              if(per_water_500m(jk) <= 0.0) then
                  per_land_500m(jk) = 100.0
              endif
          endif
      enddo
      close(202)

!>> Report out how many compartments and grid cells had missing elevation data and were assigned default values.
      write(1,*)
      write(1,*)'-----------------------------------------------'
      write(1,*)'Missing input elevation data:'
      write(1,*)'-----------------------------------------------'
      write(*,*)
      write(*,*)'-----------------------------------------------'
      write(*,*)'Missing input elevation data:'
      write(*,*)'-----------------------------------------------'
      if (bedcount > 0) then
          write(1,412) 'Compartments without input bed elev:',
     &    bedcount
          write(*,412) 'Compartments without input bed elev:',
     &    bedcount
      endif

      if (marshcount > 0) then
          write(1,412) 'Compartments without input marsh elev:',
     &     marshcount
          write(*,412) 'Compartments without input marsh elev:',
     &     marshcount
      endif

      if (gridbedcount > 0) then
          write(1,412) '500-m grid cells without input bed elev:',
     &    gridbedcount
          write(*,412) '500-m grid cells without input bed elev:',
     &    gridbedcount
      endif
      if (gridmarcount > 0) then
          write(1,412) '500-m grid cells without input marsh elev:',
     &    gridmarcount
          write(*,412) '500-m grid cells without input marsh elev:',
     &    gridmarcount
      endif
      if (gridplcount > 0) then
          write(1,412) '500-m grid cells without percent land values:',
     &    gridplcount
          write(*,412) '500-m grid cells without percent land values:',
     &    gridplcount
      endif

412   Format(7x,A,3x,I0)
!THESE AREN'T USED ANYWHERE
! Read in constants pertaining to settling data
!      write(*,*) '----------------------------------------------------'
!      write(*,*) 'READING IN PARTICLE SETTLING CONSTANTS'
!      write(*,*) '----------------------------------------------------'
!      read(31,*) Vmax                                                      !Zone Settling Css>600 mg/L
!   write(*,500) 'Maximum Settling Velocity (m/h) =  ' ,Vmax
!      Vmax = Vmax/3600.0
!
!   READ(31,*) Fsp
!   WRITE(*,500) 'Floc Settling Parameter (m^3/Kg) =    ' ,Fsp
!
!   READ(31,*) Csp
!   WRITE(*,500) 'Colloids Settling Parameter (m^3/Kg) =    ',Csp
!
!   READ(31,*) Cmin
!   WRITE(*,501) 'Concentration of nonsettling floc (Kg/m3) =',Cmin

  2   Format(x,I0,x,A)

  66  Format(A,12I)
  93  Format(x,I5, 2x,10(F8.3,x))

  500 Format (1X,A,2X,1F5.2)
  501 Format (1X,A,2X,1F6.3)

  811 Format(2x,I3,2x,18(F9.3,x))
  812 Format(2x,4(I3,x), 12(F10.3,x))

  913 Format(12(I4,2x))
  914 Format(I3,2x,11(F4.1,2x))

!>> 1D-ICM coupling input files
      if (n1D > 0) then
          write(*,*) 'Reading input files to couple 1D and 2D models'
          write(1,*) 'Reading input files to couple 1D and 2D models'
          
          write(*,*) '  - the number of terminal connections is ',ntc
          write(1,*) '  - the number of terminal connections is ',ntc          
          if (ntc>0) then
              read(402,*)                                                           ! dump header row of compartment input file
              read(402,*)                                                           ! dump header row of compartment input file
              do i = 1,ntc
                  read(402,*) tcr1D(i),                                            ! 1D region
     &                        tcn1D(i),                                                    ! 1D node            
     &                        tcr2D(i),                                                    ! ICM receiving compartment
     &                        tcf2D(i),                                                    ! ICM connecting compartment
     &                        tcl2D(i)                                                     ! ICM connecting link
              enddo
              close(402)
          endif

          write(*,*) '  - the number of lateral connections is ',nlc
          write(1,*) '  - the number of lateral connections is ',nlc
          if (nlc>0) then
              read(403,*)                                                          ! dump header row of compartment input file
              read(403,*)
              do i = 1,nlc
                  read(403,*) lcr1D(i), 
     &                        lcn1D(i),                 
     &                        lcr2D(i),
     &                        lcf2D(i),
     &                        lcl2D(i)
              enddo        
              close(403)
          endif

          write(*,*) '  - the number of upstream connections is ',nuc
          write(1,*) '  - the number of upstream connections is ',nuc
          if (nuc>0) then
              read(404,*)                                                          ! dump header row of compartment input file
              read(404,*)
              do i = 1,nuc
                  read(404,*) ucr1D(i), 
     &                        ucn1D(i),               
     &                        ucr2D(i), 
     &                        ucf2D(i), 
     &                        ucl2D(i)
              enddo        
              close(404)
          endif
      endif


!>> Read in hotstart file and set initial conditions (will overwrite some ICs set previously from input files)
      write(1,*)
      write(1,*)'-----------------------------------------------'
      write(1,*)'Reading in hotstart file and setting values as initial conditons.'
      write(1,*)'-----------------------------------------------'
      write(*,*)
      write(*,*)'-----------------------------------------------'
      write(*,*)'Reading in hotstart file and setting values as initial conditons.'
      write(*,*)'-----------------------------------------------'
      read(400,*)                       ! ignore header row
      do j=1,N
          read(400,*) node,        ! no need to save compartment number - read in to a dummy integer variable
     &                Es(node,1),       
     &                S(node,1),        
     &                Css(node,1,1),    
     &                Css(node,1,2),    
     &                Css(node,1,3),    
     &                Css(node,1,4),    
     &                Tempw(node,1),    
     &                Chem(node,1,1),   
     &                Chem(node,2,1),   
     &                Chem(node,3,1),   
     &                Chem(node,4,1),   
     &                Chem(node,5,1),   
     &                Chem(node,6,1),   
     &                Chem(node,7,1),   
     &                Chem(node,8,1),   
     &                Chem(node,9,1),   
     &                Chem(node,10,1),  
     &                Chem(node,11,1),  
     &                Chem(node,12,1),  
     &                Chem(node,13,1),  
     &                Chem(node,14,1),  !Chem unit = mg/L
     &                Eh(node,1)

      enddo
      close(400)

!>> Update the offshore boundary compartments' initial conditions 
!   - overwrite values from hotstart_in.dat to make sure they match the BC inputs from SBC.dat
!>> Assume Boundary condition sediment is evenly distributed amongst clay and silt size classes and half of the clay is flocced.
      BCSedRatio(1) = 0.0
      BCSedRatio(2) = 1./2.
      BCSedRatio(3) = 1./4.
      BCSedRatio(4) = 1./4.
      do jkk=1,mds
          jj=KBC(jkk)
          S(jj,1) = SBC(jj)
          Tempw(jj,1)=TempwBC(jkk,1)
          do i=1,4
              !Css(jj,2,i)=BCTSS(jj)*(1.+ 0.5*sin(pi*wd(kday)/180.))*BCSedRatio(i)         !BCSedRatio(k) is multipler on BCTss that separates into different classes
              Css(jj,1,i)=BCTSS(jj)*BCSedRatio(i)         !BUG wd not defined zw 04/07/2020
          enddo
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
      enddo

!>> user-specified varying time step option - ZW 1/27/2025
!     If idt_schm = 2 user specified varying time step input file
      if(idt_schm == 2) then
          dt_var_user(:)=0

!>> Skip timestep data for years prior to current model year
          do kt=1,startrun
              read(900,*)
              if (kt == startrun) then
                  write(1,*)'variable timestep input data starts at:'
                  write(1,66) '      line ',int(startrun+1)
              endif
          enddo

!>> Read variable tiemstep data for current model year
          do kt=1,simdays
              read(900,*) dt_var_user(kt)
              if (dt_var_user(kt) > dt) then  !varying time step should not greater than dt specified in RunControlR.dat
                  dt_var_user(kt) = dt
              else
                  if (mod(dt,dt_var_user(kt)) > 0) then  !variable timestep should be an integer fraction of dt
                     ndtt = int(dt/dt_var_user(kt))
                     dt_var_user(kt) = dt/ndtt
                  endif
              endif
          enddo
          close(900)
      endif

      return
      end
