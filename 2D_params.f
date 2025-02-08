      module params
      
      !determine single precision kind value
      integer,parameter :: sp=selected_real_kind(p=6)     
      
      !determine double precision kind value
      integer,parameter :: dp=selected_real_kind(p=13)

      ! diversion on/off flag used for project run - MP project # 03b.DI.04 
      integer :: Atch_div_onoff
      
      integer :: cells,links,maxconnect,simdays,lastdaystep,daystep
      integer :: lasttidestep,tidestep,tiderow,surgerow
      integer :: lastwindstep,windstep,windrow
      integer :: lastlockstep,lockstep,lockrow,lockOPstep
      integer :: windgages,raingages,etgages,tidegages 
      integer :: writehourly,modeloverland
      integer :: iSed,iSal,iTemp,iWQ

! zw 2/1/2024 added for option to select scalar advection transport scheme: 
!     1-Blended Differencing (BD); 2-User defined fa; 3-SWMM5 WQ scheme (w/o fa)
      integer :: iAdvTrans    
      real(sp):: r_BD

! zw 11/6/2024 added for option to select time step options: 
!     1-fixed; 2-Variable User Defined; 3-Varibale SWMM Scheme
      integer :: idt_schem
      real(sp),dimension(:),allocatable :: dt_var_user
      real(sp),dimension(:),allocatable :: daily_maxWL
! hardcoded specific salinity and stage trigger flags      
      integer :: SalLockTriggerHNC,SalLockStatusHNC
      integer :: StgTriggerSuperiorCanal,StgTriggerStatusSuperiorCanal
      integer :: SalLockTriggerCSC,SalLockStatusCSC
      
! filters applied to output data used by other ICM routines
      real(sp) :: stagemax
      real(sp) :: salmax
      real(sp) :: tmpmax
      real(sp) :: sedmaxow
      real(sp) :: sedmaxmi
      real(sp) :: sedmaxme
      real(sp) :: prismmax
      real(sp) :: rangemax
      real(sp) :: depthmax
      real(sp) :: algmax
      real(sp) :: tssmax
      real(sp) :: tknmax
      real(sp) :: TSSresusOff

      ! input variables to check for instabilities
      real(sp) :: oscilflag
      real(sp) :: minwater
      real(sp) :: maxdz

      
      ! input variables with error terms to perturb to output files
      real(dp) :: sal_0_1_error
      real(dp) :: sal_1_5_error
      real(dp) :: sal_5_20_error
      real(dp) :: sal_20_35_error
      real(dp) :: tss_error
      real(dp) :: stage_error
      real(dp) :: stvar_error
      
      
! new global parameters - formally local to hydrod or main
      !real(sp) :: time
      !real(sp) :: thour
      !real(sp) :: tmon
      real(sp) :: day
      !real(sp) :: dday
      !real(sp) :: hday
      !integer :: kday
      !integer :: kthr
      integer :: NTs
      real(sp) :: dtwind
      real(sp) :: dttide
      real(sp) :: dz
      real(sp) :: dry_threshold  !dry water depth threshold >0

! global water quality parameters
      real(sp) :: saltox
      real(sp) :: kphy20
      real(sp) :: thetaphy
      real(sp) :: kdet20
      real(sp) :: thetadet
      real(sp) :: kresp20
      real(sp) :: thetaresp
      real(sp) :: kpo420
      real(sp) :: thetapo4
      real(sp) :: knit20num
      real(sp) :: knit20
      real(sp) :: thetanit
      real(sp) :: knitmax
      real(sp) :: knitmin
      real(sp) :: kdenit20
      real(sp) :: thetadenit
      real(sp) :: kphot20
      real(sp) :: thetaphot
      real(sp) :: kdon20
      real(sp) :: thetadon
      

      integer :: nlockobs
      integer :: nlockobs_links
      integer :: dtlock
      real(sp),dimension(:,:),allocatable:: lockhours
      
! parameters used in hourly stage output file
      integer :: nstghr
      integer,dimension(:),allocatable :: stghrwrite
      real(sp),dimension(:),allocatable :: EShrly
      
!parameters used in flowrate output file
      integer :: nlinksw
      integer,dimension(:),allocatable :: linkswrite
      real(sp),dimension(:),allocatable :: FLO
      
! wind data arrays
      real(sp),dimension(:,:),allocatable :: windx_data
      real(sp),dimension(:,:),allocatable :: windy_data
      real(sp),dimension(:),allocatable :: windx
      real(sp),dimension(:),allocatable :: windy
! tide data arrays      
      real(sp),dimension(:,:),allocatable :: TideData
      integer,dimension(:,:),allocatable :: transposed_tide
      real(sp),dimension(:,:),allocatable :: weighted_tide
    

! Young&Verhagen wave arrays
      real(sp),dimension(:,:),allocatable :: group_vel
      real(sp),dimension(:,:),allocatable :: Hs
      real(sp),dimension(:,:),allocatable :: Uorb
      real(sp),dimension(:,:),allocatable :: wave_energy
      real(sp),dimension(:,:),allocatable :: wave_frequency
      real(sp),dimension(:,:),allocatable :: wavelength
      real(sp),dimension(:,:),allocatable :: wave_period

! Sediment distribution variables and arrays
      real(sp) :: flocA
      real(sp) :: flocB
      real(sp) :: flocN
      real(sp) :: flocM
      real(sp) :: flocC1
      real(sp) :: flocC3
      real(sp) :: Pfloc
      real(sp) :: flocset
      real(sp) :: Csalmax
      real(sp) :: Pflocmax
      real(sp) :: D90x
      real(sp) :: CSSmin
      real(sp) :: CSSmax
      real(sp),dimension(:),allocatable :: ka
      real(sp),dimension(:),allocatable :: cf
      real(sp),dimension(:),allocatable :: sedn
      real(sp),dimension(:),allocatable :: sedcalib
      real(sp),dimension(:),allocatable :: alphaSed
      real(sp),dimension(:),allocatable :: Dgr
      real(sp),dimension(:),allocatable :: D50
      real(sp),dimension(:),allocatable :: D90
      real(sp),dimension(:),allocatable :: rhoSed
      real(sp),dimension(:),allocatable :: SGsed
      real(sp),dimension(:),allocatable :: kinvisc
      real(sp),dimension(:),allocatable :: Tcrit
      real(sp),dimension(:),allocatable :: dCSS
      real(sp),dimension(:),allocatable :: dCSSh
      real(sp),dimension(:),allocatable :: dSacc
      real(sp),dimension(:),allocatable :: dSacch
      real(sp),dimension(:),allocatable :: depo_on_off
      real(sp),dimension(:),allocatable :: SedAccumRate
      real(sp),dimension(:),allocatable :: SedAccumRate_h
      real(sp),dimension(:),allocatable :: CSSresusOff
      real(sp),dimension(:),allocatable :: erBedDepth
      real(sp),dimension(:),allocatable :: erBedBD
      real(sp),dimension(:,:),allocatable :: erBedAvail
      real(sp),dimension(:,:),allocatable :: deposition
      real(sp),dimension(:,:),allocatable :: resuspension
      real(sp),dimension(:,:),allocatable :: velSet
      real(sp),dimension(:,:),allocatable :: CSSvRs
      real(sp),dimension(:),allocatable :: Qsum_in
      real(sp),dimension(:),allocatable :: Qsum_out
      real(sp),dimension(:),allocatable :: Qsum_abs
      real(sp),dimension(:),allocatable :: QSsum
      real(sp),dimension(:),allocatable :: KKa
      real(sp),dimension(:),allocatable :: KKdepth
      real(sp),dimension(:),allocatable :: MEE
      real(sp),dimension(:),allocatable :: QSmarsh
      real(sp),dimension(:),allocatable :: daccm
      real(sp),dimension(:),allocatable :: QSsumh
      real(sp),dimension(:),allocatable :: QStrib
      real(sp),dimension(:),allocatable :: QSdiv
      real(sp),dimension(:),allocatable :: netflux
      real(sp),dimension(:,:),allocatable :: Sacch_int
      real(sp),dimension(:,:),allocatable :: Sacch_edge
      real(sp),dimension(:,:),allocatable :: Sandacc
      real(sp),dimension(:,:),allocatable ::Siltacc
      real(sp),dimension(:),allocatable :: pct_sand_bed
      real(sp),dimension(:),allocatable :: pct_sand_bed_500m
      real(sp),dimension(:,:),allocatable :: Clayacc
      real(sp),dimension(:),allocatable :: BCSedRatio
      real(sp),dimension(:),allocatable :: MEESedRatio
      real(sp),dimension(:),allocatable :: cumul_retreat
      real(sp),dimension(:),allocatable :: MEESedRate
      real(sp),dimension(:),allocatable :: DivMult
      real(sp),dimension(:),allocatable :: SWRsand
      real(sp),dimension(:),allocatable :: SWRfines
      real(sp),dimension(:),allocatable :: cssFines
      real(sp),dimension(:),allocatable :: adaption_coeff   !non-equilibrium adaption coefficient for sand resuspension/deposition source term
! arrays used for BIMODE tidal prism
      real(sp),dimension(:,:),allocatable :: tidal_range_daily
      real(sp),dimension(:),allocatable :: dailyLW
      real(sp),dimension(:),allocatable :: dailyHW
      real(sp),dimension(:),allocatable :: TRsum
      real(sp),dimension(:),allocatable :: TRave
      real(sp),dimension(:),allocatable :: tidalprism_ave
      
! arrays used in ICM_formatting and InterpolateToGrid      
      real(dp) :: runtime,runtime_s,runtime_e
      integer :: runtime_start,runtime_end
      integer :: count_rate1,count_rate2,count_max1,count_max2
      character*1000 :: dump_text,dump_text2
      character*300 :: VegWaveAmpFile
      character*300 :: VegMeanSalFile  
      character*300 :: VegSummerDepthFile
      character*300 :: VegSummerSalFile 
      character*300 :: VegSummerTempFile
      character*300 :: VegTreeEstCondFile
      character*300 :: VegBIHeightFile
      character*300 :: VegPerLandFile
      character*300 :: Hydro_log
      
      integer :: dump_int
      real :: dump_float
      
      integer :: n_500m_cells
      integer :: veg_matrix_rows,veg_matrix_cols
      integer :: veg_xllcorner,veg_yllcorner
      
      integer :: n_1000m_cells
      integer :: ewe_matrix_rows,ewe_matrix_cols
      integer :: ewe_xllcorner,ewe_yllcorner
      
      integer :: year
      integer :: output_interp_flag       ! flag value corresponds to data to be interpolated to grid
      
      real(sp) :: IDW_exp
      real(sp) :: defbedelev
      real(sp) :: defmarelev
      real(sp) :: def_n
      
      integer, dimension(:), allocatable :: month_DOY
      integer, dimension(:), allocatable :: nlink2cell
      
      integer, dimension(:,:), allocatable :: grid_lookup_500m        !lookup table - relates 500 m grid cell to hydro compartment and 14 links
      real(sp), dimension(:,:), allocatable :: grid_interp_dist_500m  !lookup table - relates 500 m grid cell to distance to centorid of hydro compartment and 14 links
      real(sp), dimension(:), allocatable :: bed_elev_500m            !lookup table - relates 500 m grid cell to mean elevation of open water bed (as calculated in Wetland Morph routine)
      real(sp), dimension(:), allocatable :: land_elev_500m           !lookup table - relates 500 m grid cell to mean elevation of land (as calculated in Wetland Morph routine)
      real(sp), dimension(:), allocatable :: per_land_500m            !lookup table - relates 500 m grid cell to percentage of grid that is  land (as calculated in Wetland Morph routine)
      real(sp), dimension(:), allocatable :: per_water_500m            !lookup table - relates 500 m grid cell to percentage of grid that is  water (as calculated in Wetland Morph routine)
      real(sp), dimension(:),allocatable :: height_500m               !lookup table - relates 500 m grid cell to mean height of land above mean water level
      real(dp), dimension(:), allocatable :: SlkAve                 !average daily salinity values for links
      real(dp), dimension(:), allocatable :: TlkAve                 !average daily salinity values for links
      
      real(dp), dimension(:,:), allocatable :: sal_daily              ! daily salinity values for compartments
      real(dp), dimension(:,:), allocatable :: sal_daily_links        ! daily salinity values for links
      real(dp), dimension(:), allocatable :: sal_ave                  ! average salinity values for compartments
      real(dp), dimension(:), allocatable :: sal_ave_links            ! average salinity values for links
      real(dp), dimension(:), allocatable :: salinity_500m            ! average salinity values mapped to 500m grid cells - no smoothing/interpolation
      real(dp), dimension(:), allocatable :: salinity_IDW_500m        ! average salinity values mapped to 500m grid cells - inverse distance weighting interpolation
      real(dp), dimension(:), allocatable :: sal_2wk_ave_max          ! maximum of 2-week average salinity from March to October in compartments
      real(dp), dimension(:), allocatable :: sal_2wk_ave_max_links    ! maximum of 2-week average salinity from March to October in links
      real(dp), dimension(:), allocatable :: sal_thresh_500m          ! maximum of 2-week average salinity from March to October mapped to 500m grid cells - no smoothing/interpolation
      real(dp), dimension(:), allocatable :: sal_thresh_IDW_500m      ! maximum of 2-week average salinity from March to October mapped to 500m grid cells - inverse distance weighting interpolation
      
      real(dp), dimension(:,:), allocatable :: sal_summer             ! daily summertime salinity values for compartments
      real(dp), dimension(:,:), allocatable :: sal_summer_links       ! daily summertime salinity values for links
      real(dp), dimension(:), allocatable :: sal_ave_summer           ! average summertime salinity values for compartments
      real(dp), dimension(:), allocatable :: sal_ave_summer_links     ! average summertime salinity values for links
      real(dp),dimension(:),allocatable :: salinity_summer_500m       ! average salinity values mapped to 500m grid cells - no smoothing/interpolation
      real(dp),dimension(:),allocatable :: salinity_summer_IDW_500m   ! average salinity values mapped to 500m grid cells - inverse distance weighting interpolation
       
      real(sp), dimension(:,:), allocatable :: tmp_daily              ! daily temperature values for compartments
      real(sp), dimension(:,:), allocatable :: tmp_daily_links        ! daily temperature values for links
      real(sp), dimension(:), allocatable :: tmp_ave                  ! average temperature values for compartments
      real(sp), dimension(:), allocatable :: tmp_ave_links            ! average temperature values for links
      real(sp),dimension(:),allocatable :: tmp_500m                   ! average temperature values mapped to 500m grid cells - no smoothing/interpolation
      real(sp),dimension(:),allocatable :: tmp_IDW_500m               ! average temperature values mapped to 500m grid cells - inverse distance weighting interpolation
      
      real(sp), dimension(:,:), allocatable :: tkn_daily              ! daily TKN values for compartments
      real(sp), dimension(:,:), allocatable :: tkn_daily_links        ! daily TKN values for links
      real(sp), dimension(:,:), allocatable :: tss_daily              ! daily TSS values for compartments
      real(sp), dimension(:,:), allocatable :: tss_daily_links        ! daily TSS values for links 
      real(sp), dimension(:), allocatable :: tss_ave                  ! average TSS values for compartments
      real(sp), dimension(:), allocatable :: tss_var_annual
      real(sp), dimension(:,:), allocatable :: difmean_tss

      
      
      ! arrays for HSI monthly values
      real(sp), dimension(:,:), allocatable :: stg_month_ave            ! monthly average stage values for compartments
      real(dp), dimension(:,:), allocatable :: sal_month_ave          ! monthly average salinity values for compartments
      real(dp), dimension(:,:), allocatable :: sal_month_ave_links    ! monthly average salinity values for links
      real(sp), dimension(:,:), allocatable :: tmp_month_ave          ! monthly average temperature values for compartments
      real(sp), dimension(:,:), allocatable :: tmp_month_ave_links    ! monthly average temperature values for links
      real(sp), dimension(:,:), allocatable :: tkn_month_ave          ! monthly average TKN values for compartments
      real(sp), dimension(:,:), allocatable :: tkn_month_ave_links    ! monthly average TKN values for links
      real(sp), dimension(:,:), allocatable :: tss_month_ave          ! monthly average TSS values for compartments
      real(sp), dimension(:,:), allocatable :: tss_month_ave_links    ! monthly average TSS values for links
      real(dp), dimension(:,:), allocatable :: sal_IDW_500m_month     ! average monthly salinity values mapped to 500m grid cells - inverse distance weighting interpolation
      real(sp),dimension(:,:),allocatable :: tmp_IDW_500m_month       ! average monthly temperature values mapped to 500m grid cells - inverse distance weighting interpolation
      real(sp),dimension(:,:),allocatable :: tkn_500m_month           ! average monthly algae values mapped to 500m grid cells - no smoothing/interpolation
      real(sp),dimension(:,:),allocatable :: tss_500m_month           ! average monthly TSS values mapped to 500m grid cells - no smoothing/interpolation
            
      real(sp), dimension(:,:), allocatable :: tmp_summer             ! daily temperature values for compartments
      real(sp), dimension(:,:), allocatable :: tmp_summer_links       ! daily temperature values for links
      real(sp), dimension(:), allocatable :: tmp_ave_summer           ! average temperature values for compartments
      real(sp), dimension(:), allocatable :: tmp_ave_summer_links     ! average temperature values for links
      real(sp),dimension(:),allocatable :: tmp_summer_500m            ! average temperature values mapped to 500m grid cells - no smoothing/interpolation
      real(sp),dimension(:),allocatable :: tmp_summer_IDW_500m        ! average temperature values mapped to 500m grid cells - inverse distance weighting interpolation
          
      real(sp), dimension(:,:), allocatable :: stage_daily            ! daily stage values for compartments
      real(sp), dimension(:), allocatable :: stage_ave                ! average stage values for compartments
      real(sp), dimension(:), allocatable :: stage_max                ! maximum stage values for compartments
      real(sp), dimension(:), allocatable :: sepmar_stage               ! average stage for September-March (Green-winged Teel HSI)
      real(sp), dimension(:), allocatable :: octapr_stage             ! average stage for October-April (Gadwall HSI)
      real(sp), dimension(:), allocatable :: stage_500m               ! average stage values mapped to 500m grid cells - no smoothing/interpolation
      real(sp), dimension(:), allocatable :: depth_500m               ! average depth values mapped to 500m grid cells - no smoothing/interpolation
             
      real(sp), dimension(:,:), allocatable :: stage_summer           ! daily summertime stage values for compartments
      real(sp), dimension(:), allocatable :: stage_ave_summer         ! average summertime stage values for compartments
      real(sp), dimension(:), allocatable :: stage_summer_500m        ! average summertime stage values mapped to 500m grid cells - no smoothing/interpolation
      real(sp), dimension(:), allocatable :: depth_summer_500m        ! average summertime depth values mapped to 500m grid cells - no smoothing/interpolation
      integer, dimension(:), allocatable :: tree_est                  ! tree establishment condition for 500 m grid cells
      
      real(sp), dimension(:,:), allocatable :: difmean                ! temporary array used in calculating variance in stage
      real(sp), dimension(:), allocatable :: stage_var_summer         ! variance in summertime stage for compartments
      real(sp), dimension(:), allocatable :: stage_stdv_summer        ! standard deviation of summertime stage for compartments

      real(sp), dimension(:,:), allocatable :: trg_summer             ! daily summertime tidal range values for compartments
      real(sp), dimension(:), allocatable :: trg_ave_summer           ! average summertime tidal range for compartments
      
      real(sp), dimension(:), allocatable :: stage_wlv_summer         ! water level variablilty in summertime stage
      real(sp), dimension(:), allocatable :: stage_wlv_summer_500m    ! water level variablilty in summertime stage mapped to 500m grid cells (as used by ICM-LAVegMod)
      
      real(sp), dimension(:), allocatable :: SedOW                    ! annual sediment accumulation in Open Water portion of compartment
      real(sp), dimension(:), allocatable :: SedMarshInt              ! annual sediment accumulation in interior of Marsh portion of compartment
      real(sp), dimension(:), allocatable :: SedMarshEdge             ! annual sediment accumulation in edge of Marsh portion of compartment
      
      integer,dimension(:,:),allocatable :: veg_grid_IDs              ! lookup matrix of grid IDs for Veg model formatting
      real(sp),dimension(:,:),allocatable :: depth_summer_forVeg      ! mapped summertime depth values formatted for Veg model
      real(dp),dimension(:,:),allocatable :: salinity_forVeg          ! interpolated salinity values formatted for Veg model
      real(dp),dimension(:,:),allocatable :: salinity_summer_forVeg   ! interpolated summertime salinity values formatted for Veg model
      real(sp),dimension(:,:),allocatable :: stage_var_forVeg         ! mapped stage variance values formatted for Veg model
      real(sp),dimension(:,:),allocatable :: tmp_summer_forVeg        ! interpolated summertime temp values formatted for Veg model
      integer,dimension(:,:),allocatable :: tree_est_forVeg           ! mapped tree establishment conditions formatted for Veg model
      real(sp),dimension(:,:),allocatable :: ht_abv_water_forVeg      ! mapped height above mean water level formatted for Veg model
      real(sp),dimension(:,:),allocatable :: per_land_forVeg          ! mapped percentage land of grid cell formatted for Veg model
          
      ! new link attributes
      real(sp),dimension(:),allocatable :: USx
      real(sp),dimension(:),allocatable :: USy
      real(sp),dimension(:),allocatable :: DSx
      real(sp),dimension(:),allocatable :: DSy
      integer, dimension(:),allocatable :: linkt
      real(sp),dimension(:),allocatable :: Latr1
      real(sp),dimension(:),allocatable :: Latr2
      real(sp),dimension(:),allocatable :: Latr3
      real(sp),dimension(:),allocatable :: Latr4
      real(sp),dimension(:),allocatable :: Latr5
      real(sp),dimension(:),allocatable :: Latr6
      real(sp),dimension(:),allocatable :: Latr7
      real(sp),dimension(:),allocatable :: Latr8
      real(sp),dimension(:),allocatable :: Latr9
      real(sp),dimension(:),allocatable :: Latr10
      real(sp),dimension(:),allocatable :: Latr11      
      real(sp),dimension(:),allocatable :: fa_mult

      integer,dimension(:,:),allocatable :: hourclosed
      
      ! new variables calculatedin photo.f and passed to various WQ routines
      real(sp) :: fpp
      real(sp) :: muph
      
! old comdeck.h memory block: COMMON/hyd1/
!      integer :: kyear                                       !replaced with reading year into model run
      
      real(sp) :: Clams
      real(sp) :: fFetch
      real(sp),dimension(:),allocatable :: ParSandD
      real(sp) :: TP1
!      real(sp) :: ESLR                                        !variable not used
!      real(sp) :: ESLRd                                       !variable not used
!      real(sp) :: NBND                                        !variable not used
!      real(sp) :: rslr                                        !variable not used

      integer, dimension(:), allocatable :: Jrain     !should this be integer or real? implicitly defined in original
      integer, dimension(:), allocatable :: KBC       !should this be integer or real? implicitly defined in original
      integer, dimension(:), allocatable :: Jwind
      integer, dimension(:), allocatable :: Jet
      integer, dimension(:), allocatable :: flag_offbc !zw offshore bc cells flag 04/07/2020
      
!      real(sp), dimension(:), allocatable :: acss
      real(sp), dimension(:), allocatable :: Ahydro
      real(sp), dimension(:), allocatable :: Apctmarsh
      real(sp), dimension(:), allocatable :: Apctwater
      real(sp), dimension(:), allocatable :: Apctupland
      real(sp), dimension(:), allocatable :: Atotal
!      real(sp), dimension(:), allocatable :: bcss
      real(sp), dimension(:), allocatable :: Bed
      real(sp), dimension(:), allocatable :: dAdz
      real(sp), dimension(:), allocatable :: Eso
      real(sp), dimension(:), allocatable :: Percent
      real(sp), dimension(:), allocatable :: SBC
      real(sp), dimension(:), allocatable :: SSource      
!      real(sp), dimension(:), allocatable :: tauc
!      real(sp), dimension(:), allocatable :: Vss
            
             
      integer, dimension(:,:), allocatable :: icc            
      
      real(sp), dimension(:,:), allocatable :: Aas
      real(dp), dimension(:,:), allocatable :: As
      real(sp), dimension(:,:,:), allocatable :: CSS !third array dimension for 4 different sediment classes
      real(sp), dimension(:), allocatable :: CSSos 
      real(sp), dimension(:,:), allocatable :: daymo
      real(sp), dimension(:,:), allocatable :: ds
      real(sp), dimension(:,:), allocatable :: ESAV
      real(sp), dimension(:,:), allocatable :: ESMN
      real(sp), dimension(:,:), allocatable :: ESMX
      real(sp), dimension(:,:), allocatable :: Fetch
      real(sp), dimension(:,:), allocatable :: sicc
!      real(sp), dimension(:,:), allocatable :: SA             !variable not used
!      real(sp), dimension(:,:), allocatable :: tau
      real(sp), dimension(:,:), allocatable :: EHAV
      
      real(dp), dimension(:,:), allocatable :: Es
      real(dp), dimension(:,:), allocatable :: BCnosurge      ! ES value for BC comparments before surge is added
      real(dp), dimension(:,:), allocatable :: BCsurge        !YW! added for tide calculation
      real(dp), dimension(:,:), allocatable :: S 
      real(dp), dimension(:,:), allocatable :: SL             !changed to double precision !-EDW
      real(dp), dimension(:,:), allocatable :: STEMP          !changed to double precision !-EDW
      real(dp), dimension(:), allocatable :: SALAV  !zw added 3/22/2015 Salinity daily average from each time step

      
! old comdeck.h memory block: COMMON/hyd2/
      real(sp) :: cden
      real(sp) :: fa_def
      real(sp), dimension(:), allocatable :: fa
      real(sp), dimension(:), allocatable :: fb
      real(sp) :: upwind_vel
      real(sp) :: fbc
      real(sp) :: fe
      real(sp) :: KKexp
      
!      real(sp), dimension(:), allocatable :: Achan
      real(sp), dimension(:), allocatable :: an
      real(sp), dimension(:), allocatable :: AnthL
!      real(sp), dimension(:), allocatable :: Cinvert
!      real(sp), dimension(:), allocatable :: Deptho
      real(sp), dimension(:), allocatable :: EAOL
      real(sp), dimension(:), allocatable :: Exy
!      real(sp), dimension(:), allocatable :: Ken
!      real(sp), dimension(:), allocatable :: Km
!      real(sp), dimension(:), allocatable :: Kx
!      real(sp), dimension(:), allocatable :: Length
!      real(sp), dimension(:,:), allocatable :: Qss
!      real(sp), dimension(:,:), allocatable :: Qsal
!      real(sp), dimension(:), allocatable :: Resist ! no longer an array - now variable local to hydrod
!      real(sp), dimension(:), allocatable :: Width

!      integer, dimension(:), allocatable :: itype
      integer, dimension(:), allocatable :: jus
      integer, dimension(:), allocatable :: jds

      real(sp), dimension(:,:), allocatable :: denit
      real(sp), dimension(:,:), allocatable :: Depth

      real(dp), dimension(:,:), allocatable :: Q                  ! discharge [m3/s] in each link - negative rate indicates flow into downstream compartment from upstream for the link
      real(dp), dimension(:), allocatable :: link_vel             ! calculated velocity [m/s] in each link - 
      real(dp), dimension(:), allocatable :: ave_vel_sum          ! array used to accumulate the numerator when calculating the average velocity magnitude of all links flowing into/out of each compartment
      real(dp), dimension(:), allocatable :: ave_vel_cnt          ! array used to accumulate the denominator when calculating the average velocity magnitude of all links flowing into/out of each compartment
      real(dp), dimension(:), allocatable :: ave_vel              ! average velocity magnitude [m/s] in each compartment during current timestep - determined by looping over all connecting links and averaging velocities together, neglecting directionality
      real(dp), dimension(:), allocatable :: max_vel              ! maximum velocity magnitude [m/s] for each compartment during current timestep
      real(dp), dimension(:), allocatable :: min_vel              ! minimum velocity magnitude [m/s] for each compartment during current timestep
      
      real(sp), dimension(:,:,:), allocatable :: GrowAlgae
      real(sp), dimension(:,:,:), allocatable :: GrowChlA
      
! old comdeck.h memory block: COMMON/sed/
!      real(sp) :: fcbc !-EDW not used
!      real(sp) :: fcss !-EDW no longer used
!      real(sp) :: fvs  !-EDW no longer used
      real(sp) :: PP
!      real(sp) :: VsettL !-EDW not used
      
      real(sp), dimension(:), allocatable :: accsed
      real(sp), dimension(:), allocatable :: ACCSEDj
      real(sp), dimension(:), allocatable :: ASandA
      real(sp), dimension(:), allocatable :: asedout
      real(sp), dimension(:), allocatable :: CssL
!      real(sp), dimension(:), allocatable :: NR
!      real(sp), dimension(:), allocatable :: Pmsh
      real(sp), dimension(:), allocatable :: por
      real(sp), dimension(:), allocatable :: Qsed
      real(sp), dimension(:), allocatable :: Qsedsm
      real(sp), dimension(:), allocatable :: SWR
      real(sp), dimension(:), allocatable :: Ws
      
      real(sp), dimension(:,:), allocatable :: AsandD
      real(sp), dimension(:,:), allocatable :: AsandT
      real(sp), dimension(:,:), allocatable :: QssT      
      
! old comdeck.h memory block: COMMON/Forcing/
!      real(sp) :: ETA
      real(sp),dimension(:),allocatable :: ETA
      real(sp) :: fpet
      real(sp) :: g
      real(sp) :: ParChlA
      real(sp) :: ParSand
!      real(sp) :: RSSSo
      real(dp) :: t

!      real(dp) :: RSSS !-EDW no longer used
!      real(dp) :: SEARD !-EDW no longer used
      
      real(sp), dimension(:), allocatable :: cChemface
      real(sp), dimension(:), allocatable :: DChem
      real(sp), dimension(:), allocatable :: QChemSUM
      real(sp), dimension(:), allocatable :: QChemSUManth
      real(sp), dimension(:), allocatable :: QChemSUMtrib
      real(sp), dimension(:), allocatable :: QChemSUMdiv
      real(sp), dimension(:), allocatable :: QChemSUMflows
      real(sp), dimension(:), allocatable :: QChemSUMatm
      real(sp), dimension(:), allocatable :: CSSo
      real(sp), dimension(:), allocatable :: TSS
      real(sp), dimension(:), allocatable :: TSSAve !-EDW new array to save daily average WQ values

      real(sp), dimension(:,:), allocatable :: PET
      real(sp), dimension(:), allocatable :: T7
      
      real(sp), dimension(:,:), allocatable :: Sal
      real(sp), dimension(:,:), allocatable :: Strib
      real(sp), dimension(:,:), allocatable :: Qtrib
      real(sp), dimension(:,:), allocatable :: QMult
      real(sp), dimension(:,:), allocatable :: QMultdiv
      real(sp), dimension(:,:), allocatable :: Rain
!      real(sp), dimension(:,:), allocatable :: uplandNP
      real(sp), dimension(:,:), allocatable :: Tcoef
      real(sp), dimension(:), allocatable :: QRain  !ZW 1/31/2024 save total runoff volume (m3/s) at each time step
      
      real(sp), dimension(:,:,:), allocatable :: cChemdiv
      real(sp), dimension(:,:,:), allocatable :: QAtm
      real(sp), dimension(:,:,:), allocatable :: cCHEM
      real(sp), dimension(:,:,:), allocatable :: QChemdiv
      real(sp), dimension(:,:,:), allocatable :: Chem
      real(sp), dimension(:,:), allocatable :: ChemAve !-EDW new array to save daily average WQ values
      real(sp), dimension(:,:,:), allocatable :: QCHEM
      
! old comdeck.h memory block: COMMON/Control/
      integer :: M
      integer :: Mds
      integer :: Mus
      integer :: N
      integer :: Ntrib
      
      real(sp) :: consd
      real(sp) :: conv
      real(sp) :: dt
      real(sp) :: dTprint
      real(sp) :: endrun
      real(sp) :: pi
      real(sp) :: Specg
      real(sp) :: startrun
      real(sp) :: STDS
      real(sp) :: Tcorn

      integer,dimension(:),allocatable :: jdiv
      integer,dimension(:),allocatable :: jtrib
      
      real(sp),dimension(:,:,:),allocatable :: cssT
      real(sp),dimension(:,:,:),allocatable :: cSSTdiv
      real(sp),dimension(:,:),allocatable :: QSSTdiv
      real(sp),dimension(:,:),allocatable :: EsRange

! old comdeck.h memory block: COMMON/Ftide/
      real(sp),dimension(:),allocatable :: AM1     
      real(sp),dimension(:),allocatable :: AM2
      real(sp),dimension(:),allocatable :: phz
      real(sp),dimension(:),allocatable :: TM1
      real(sp),dimension(:),allocatable :: TM2
      
! old comdeck.h memory block: COMMON/seas/
!     real(sp) :: RSLS !-EDW not used anymore
      real(sp),dimension(:,:),allocatable :: BL
      real(sp),dimension(:,:),allocatable :: surge

! old comdeck.h memory block: COMMON/MET/
      integer :: nbar
      
      real(sp) :: Cp
      
      integer,dimension(:),allocatable :: jbar
      real(sp),dimension(:),allocatable :: lnO2Sat
      real(sp),dimension(:),allocatable :: O2Sat
      real(sp),dimension(:),allocatable :: ta
      real(sp),dimension(:),allocatable :: ta_k
      real(sp),dimension(:),allocatable :: tw
      real(sp),dimension(:),allocatable :: wd
!      real(sp),dimension(:),allocatable :: wspd
      
      real(sp),dimension(:,:),allocatable :: barp
      real(sp),dimension(:,:),allocatable :: Tmp
      real(sp),dimension(:,:),allocatable :: Tmpe
      real(sp),dimension(:,:),allocatable :: Tmpm
      
! old comdeck.h memory block: COMMON/Qdivin
      integer :: NDIV
      
! old comdeck.h memory block: COMMON/Chem
      integer :: Nyear
      integer :: NCCC
      
      real(sp) :: Vsettoc
      
      real(sp),dimension(:),allocatable :: SOD
      
      real(sp),dimension(:,:),allocatable :: decay
      real(sp),dimension(:,:),allocatable :: PAR
      real(sp),dimension(:,:),allocatable :: Qdiv
      
      real(sp),dimension(:,:,:),allocatable :: dccc1
      real(sp),dimension(:,:,:),allocatable :: dccc2
      
! old comdeck.h memory block: COMMON/Michaelis
      real(sp) :: KnN
      real(sp) :: KnP
      real(sp) :: KnSS
      real(sp) :: KnSal

! old comdeck.h memory block: COMMON/NutroOB
      real(sp),dimension(:),allocatable :: BCDIN      
      real(sp),dimension(:),allocatable :: BCNH4
      real(sp),dimension(:),allocatable :: BCNO3
      real(sp),dimension(:),allocatable :: BCON
      
! old comdeck.h memory block: COMMON/OBC1
      real(sp),dimension(:),allocatable :: BCage
      real(sp),dimension(:),allocatable :: BCDO
      real(sp),dimension(:),allocatable :: BCTP

! old comdeck.h memory block: COMMON/OBC
      real(sp),dimension(:),allocatable :: BCDA
      real(sp),dimension(:),allocatable :: BCLA
      real(sp),dimension(:),allocatable :: BCTOC
      real(sp),dimension(:),allocatable :: BCTSS

      
! old comdeck.h memory block: COMMON/Marsh
      real(sp),dimension(:),allocatable :: Ahf
!      real(sp),dimension(:),allocatable :: Eho
!      real(sp),dimension(:),allocatable :: flood
      real(sp),dimension(:),allocatable :: floodf

      real(dp),dimension(:,:),allocatable :: Eh
      
! old comdeck.h memory block: COMMON/Marsh2
      real(sp),dimension(:),allocatable :: acssh
!      real(sp),dimension(:),allocatable :: Vsh

      real(sp),dimension(:,:),allocatable :: Cssf
      real(sp),dimension(:,:),allocatable :: Sh
      
! old comdeck.h memory block: COMMON/Marsh3
      real(sp),dimension(:),allocatable :: BedM
      real(sp),dimension(:),allocatable :: BedMOrig
      real(sp),dimension(:),allocatable :: BedMSD
      real(sp),dimension(:),allocatable :: BedMAdj
      real(sp),dimension(:),allocatable :: Esho
      
      real(sp),dimension(:,:,:),allocatable :: CSSh !third array dimension is for 4 different sediment classes
      real(sp),dimension(:,:),allocatable :: Shg

! old comdeck.h memory block: COMMON/Marsh4
      real(sp),dimension(:),allocatable :: hLength
      real(sp),dimension(:),allocatable :: SourceBM
      
      real(sp),dimension(:,:),allocatable :: hmarsho

! old comdeck.h memory block: COMMON/Marsh5
      real(sp),dimension(:),allocatable :: hwidth
      real(sp),dimension(:),allocatable :: ar_ed
      real(sp),dimension(:),allocatable :: ar_int
!      real(sp),dimension(:),allocatable :: han
!      real(sp),dimension(:),allocatable :: hkm
      
      real(sp),dimension(:,:),allocatable :: Qmarsh
      real(sp),dimension(:),allocatable :: Qmarshmax
      real(sp),dimension(:),allocatable :: QmarshAve
      
! old comdeck.h memory block: COMMON/Marsh6
      real(sp),dimension(:,:),allocatable :: Sacc
      real(sp),dimension(:,:),allocatable :: Sacch
      
      
! old comdeck.h memory block: COMMON/Tempr
      real(sp),dimension(:),allocatable :: TempMR
      
      real(sp),dimension(:,:),allocatable :: Tempair
      real(sp),dimension(:,:),allocatable :: Tempe
      real(sp),dimension(:,:),allocatable :: Tempw
      real(sp),dimension(:),allocatable :: TempwAve !-EDW new array for saving daily average temp
      real(sp),dimension(:,:),allocatable :: TempwBC  !-EDW new array for reading in BCs
      real(sp),dimension(:,:),allocatable :: TL
      
! old comdeck.h memory block: COMMON/Waterage
      real(sp),dimension(:),allocatable :: HDTo
      
      real(sp),dimension(:,:),allocatable :: HDT
      real(sp),dimension(:,:),allocatable :: TMtrib
      
! old comdeck.h memory block: COMMON/ChLA1
      real(sp) :: ParCla
      
      real(sp),dimension(:),allocatable :: ChLAo
      
      real(sp),dimension(:,:),allocatable :: ChLa
      
! old comdeck.h memory block: COMMON/SRP
      real(sp) :: ParP
      real(sp) :: ParPMR
      
      real(sp),dimension(:),allocatable :: SRPo
      
      real(sp),dimension(:,:),allocatable :: SRP
      
! old comdeck.h memory block: COMMON/DON1
      real(sp) :: STchN
      
      real(sp),dimension(:),allocatable :: ChemDONo
      
      real(sp),dimension(:,:),allocatable :: ChemDON
      
! old comdeck.h memory block: COMMON/DOP1
      real(sp) :: ParDOP
      real(sp) :: PARPOP
      real(sp) :: PPMR
      real(sp) :: STchP
      
      real(sp),dimension(:),allocatable :: ChemDOPo
      
      real(sp),dimension(:,:),allocatable :: ChemDOP
      real(sp),dimension(:,:),allocatable :: ChemPOP

      integer :: nlinklimiter                            !YW! number of link to apply flow limiter   
      integer,dimension(:),allocatable :: linkslimiter   !YW! list of link number to apply flow limiter      
      integer,dimension(:),allocatable :: flag_apply     !ZW! flag to apply flow limiter in links (1-apply; 0-no)     
      
      integer :: numChem  !add zw 04/08/2020      
      
!1D-ICM coupling variables
      integer :: ntim_all_ICM, ndt_all_ICM!, ndt_ICM       ! time stepping parameters used across 1D and 2D codes
      integer :: ntc,nlc,nuc                      ! number of terminal, lateral, and upstream connection
      integer :: n_region,n1D                     ! number of 1D rivers
    
      integer,dimension(:),allocatable :: tcr2D   ! terminal connection ICM receving compartment
      integer,dimension(:),allocatable :: tcf2D   ! terminal connection ICM connecting compartment      
      integer,dimension(:),allocatable :: tcl2D   ! terminal connection ICM link      
      integer,dimension(:),allocatable :: tcn1D   ! terminal connection 1D connection node
      integer,dimension(:),allocatable :: tcr1D   ! terminal connection 1D region  

      integer,dimension(:),allocatable :: lcr2D   ! lateral connection ICM receving compartment
      integer,dimension(:),allocatable :: lcf2D   ! lateral connection ICM connecting compartment      
      integer,dimension(:),allocatable :: lcl2D   ! lateral connection ICM link
      integer,dimension(:),allocatable :: lcn1D   ! lateral connection 1D connection node
      integer,dimension(:),allocatable :: lcr1D   ! lateral connection 1D region        
      
      integer,dimension(:),allocatable :: ucr2D   ! upstream connection ICM upstream compartment
      integer,dimension(:),allocatable :: ucf2D   ! upstream connection ICM connecting compartment      
      integer,dimension(:),allocatable :: ucl2D   ! upstream connection ICM link
      integer,dimension(:),allocatable :: ucn1D    ! upstream connection 1D connection node
      integer,dimension(:),allocatable :: ucr1D   ! upstream connection 1D region         
      
      real(dp),dimension(:),allocatable :: tcH    ! terminal connection stage (applied to the connecting compartment from ICM)
      real(dp),dimension(:),allocatable :: tcQ    ! terminal connection discharge (applied to the connecting link from 1D)
      real(dp),dimension(:),allocatable :: lcH    ! lateral connection stage (applied to the connecting compartment from 1D)
      real(dp),dimension(:),allocatable :: lcQ    ! lateral connection discharge (calculated)
      real(dp),dimension(:),allocatable :: ucH    ! upstream connection stage (applied to the connecting compartment from 1D)
      real(dp),dimension(:),allocatable :: ucQ    ! upstream connection discharge (calculated)      
      
      real(dp),dimension(:),allocatable :: NTs_Ratio  ! time step ratio between 1D river and ICM
      
      end module params