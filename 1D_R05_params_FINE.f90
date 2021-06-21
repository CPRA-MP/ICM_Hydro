!-------------------------------------------------------------------------------
!                                                                               
!  Developer: Vasilia Velissariou <vvelissariou@tulane.edu>                               
!                                                                                         
!  Version: 0.2
!
!    Version - 0.1 Thu Nov 21 2019
!            - Basic setup
!    Version - 0.2 Jan 08 2020
!            - Modifications to inlude all latest options and arrays
!-------------------------------------------------------------------------------
!
module params_SAND_R05

  use precision

  implicit none

  real(fp), parameter :: Grav = 9.80665      ! gravitational acceleration (m/s2)

  real(fp), parameter :: DensW0 = 999.89     ! reference density of water at T = 0 degrees Celcius (kg/m3 = g/L)
  real(fp), parameter :: DensW4 = 999.95     ! reference density of water at T = 4 degrees Celcius (kg/m3 = g/L)
  real(fp), parameter :: DensW6 = 999.99     ! reference density of water at T = 6 degrees Celcius (kg/m3 = g/L)
  real(fp), parameter :: DensW = DensW6      ! default value

  real(hp), parameter :: ViscW0 = 1.7919E-6  ! kinematic viscosity for clear water at T = 0 degrees Celcius (m2/s)
  real(hp), parameter :: ViscW4 = 1.5707E-6  ! kinematic viscosity for clear water at T = 4 degrees Celcius (m2/s)
  real(hp), parameter :: ViscW6 = 1.4749E-6  ! kinematic viscosity for clear water at T = 6 degrees Celcius (m2/s)
  real(hp), parameter :: ViscW = ViscW6      ! default value

  real(fp), parameter :: SpSed = 2.65        ! specific gravity of quartz sediments (dimensionless)

  real(fp), parameter :: TSSresusOff = 250.0   ! TSSresusOff: CSS concentration threshold for bed resuspension -
                                               ! if sediment class CSS is greater than this threshold value,
                                               ! bed resuspension will be turned off for sediment class
                                               ! (mg/L = g/m3)

  ! This is the reference date for the simulation. Tstart and Tstop are times [s]
  ! from the reference date. - (Currently RefDate is not used)
  !character(len=19) :: RefDate
  real(fp) :: sand_init, Tstart, Tstop, Dt
  integer  :: Nclass                            ! number of sediment classes to be considered
  integer  :: NDt                               ! number of time steps in the simulation
  integer  :: NDx                               ! number of spatial increments (cross-sections) in the simulation

  real(fp) :: DtUser                            ! write time interval [s]
  integer  :: NDtUser                           ! number of write time intervals (simulation time/DtUser)
  
  ! for lateral flow khu
  integer :: Nlat                                ! No 0f lateral flows
  integer, allocatable :: Idlat(:), Nsclat(:), Ntplat(:)    ! Lateral flow id and # of CSs

  ! Longitudinal dispersion coefficient for the sediments
  real(fp) :: Ks

  ! Do we want a variable Manning's coefficient in the future?
  real(fp) :: Cmann                             ! Manning's coefficient, a good value is 0.023

  ! sand D50(1) = 0.001      -- USDA particle size definition: 0.05 mm < sand < 2 mm
  ! silt D50(2) = 0.00003    -- USDA particle size definition: 0.002 mm < silt < 0.05 mm
  ! clay D50(3) = 0.000001   -- USDA particle size definition: clay < 0.002 mm
  real(fp), dimension(:), allocatable :: D50, D90, Dgr
  real(fp) :: D90x                             ! representative D90/D50 ratio (a good value is: D90x = 1.5)

  ! sand Tcrit(1) = 0.47      -- critical shear stress for sand particlss (N/m**2)
  ! silt Tcrit(2) = 0.08      -- critical shear stress for silt particlss (N/m**2)
  ! silt Tcrit(2) = 0.045     -- critical shear stress for clay particles (N/m**2)
  real(fp), dimension(:), allocatable :: Tcrit

  real(fp), dimension(:), allocatable :: CSS ! array dimension for "nclass" different sediment classes
  real(fp), dimension(:), allocatable :: SGsed, rhoSed, velSet

  real(fp), dimension(:), allocatable :: deposition, resuspension, SedAccumRate

  real(fp) :: Cf, alphaSED, SEDn, SEDcalib, ws_fine
 
  real(fp), allocatable :: CSSresusOff(:)

  real(fp), allocatable :: xx(:), time(:), depth(:), area(:), flow(:), width(:)
  real(fp), allocatable :: CN(:, :), CN1(:, :)
  
  !real(fp), allocatable :: C(:, :, :), CLD(:, :, :), ACLD(:, :, :)
  !real(fp), allocatable :: timeUser(:), depthUser(:, :), areaUser(:, :), flowUser(:, :)

  ! C0 = initial condition array for the concentration
  real(fp), allocatable :: C0(:)

  logical :: useSRC, useBOTEROS
  integer :: BC_Option, SED_STerm, BOT_Eros
  
  ! for changing manning
  real(kind=4) ::  Q_Sk_Table(2,100)
  integer :: Q_sk_tableEntry
  
  ! from main
    ! For upstream BCs with the concentration is specified
  integer  :: nxfl, ntfl
  real(fp), allocatable :: DepthFL(:, :), AreaFL(:, :), FlowFL(:, :), TimeFL(:, :)

  real(fp), allocatable :: tmp1D(:), tmp2D(:, :)

  ! For upstream BCs with the concentration is specified
  integer  :: ntcb, npcb
  real(fp), allocatable :: ConcBND(:, :), TimeBND(:, :)

  ! For downstream BCs with the concentration is specified
  integer  :: ntdb, npdb
  real(fp), allocatable :: ConcDOW(:, :), TimeDOW(:, :)

! For lateral flow with the concentration is specified
  integer  :: ntlat, nplat
  real(fp), allocatable :: ConcLAT(:, :), TimeLAT(:, :), FlowLAT(:, :)
  
  
  ! input bed level
  !real(fp), allocatable :: zz(:)

  ! These are for the spline interpolations
  real(fp), allocatable :: Bcoef(:), Ccoef(:), Dcoef(:)
  real(fp), allocatable :: Bdp(:, :), Cdp(:, :), Ddp(:, :)
  real(fp), allocatable :: Bar(:, :), Car(:, :), Dar(:, :)
  real(fp), allocatable :: Bfl(:, :), Cfl(:, :), Dfl(:, :)

  real(fp) :: smallC = 0.0_fp

  real(fp) :: dx, dx1, dx2, lmd1, lmd2
  real(fp) :: a0, a1, a2, b0, b1, b2, b3, denom
  real(fp) :: mu, xi
  real(fp) :: aa, bb, cc, dd, ee, arni

  real(fp) :: diffVAL, accSL, bdjSL, trapEFF
  real(fp) :: SrcTerm

  real(fp), allocatable :: sumConc(:)  !, mACLD(:, :, :), diffACLD(:, :, :)
  
  real(fp), allocatable :: bed_t(:,:)   ! sediment bed thickness at cross section for particle size class
  real(fp) :: bed_bd                    ! bulk density of bed sediments
  
  character(len=256) :: out_dir


end module params_SAND_R05
