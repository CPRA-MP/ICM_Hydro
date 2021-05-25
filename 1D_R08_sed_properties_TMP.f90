!-------------------------------------------------------------------------------
!
!  Developer: Vasilia Velissariou <vvelissariou@tulane.edu>
!
!  Version: 0.1
!
!    Version - 0.1 Mon Dec 16 2019
!            - Development of the module that contains
!              many functions related to settling velocities, critical shear streeses
!              van Rijn bottom erosion, etc
!-------------------------------------------------------------------------------
!
module sed_properties_TMP_R08

use params_TMP_R08

contains


!!!========================================!!!
function svel0(sdiam, sgrav)

!  use params

  implicit none

!!!--- Return value
  real(hp) :: svel0

!!!--- Global variables
  real(fp), intent(in)       :: sdiam
  real, optional, intent(in) :: sgrav

!!!--- Local variables
  real     :: sg
  real(hp) :: vis, gsg
  real(hp) :: DSED0, DSED1, DSED100, DSED1000
  real(hp) :: velo

!
!!!--- executable statements ---
!
  !!! default kinematic viscosity of water (m^2/s)
  vis  = ViscW

  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! default
  endif

  gsg = Grav * ( sg - 1.0_hp )

  !!! The following equations follow after van Rijn:
  !!! Principles of Sediment Transport In Rivers, Estuaries
  !!! And Coastal Seas, Leo C. van Rijn, 1993
  !!! Equations 3.2.21, 3.2.22 and 3.2.23
  DSED0    = 0.0_hp     ! meters
  DSED1    = 1.0D-6     ! meters
  DSED100  = 100.0D-6   ! meters
  DSED1000 = 1000.0D-6  ! meters
  if ((DSED0 <= sdiam) .and. (sdiam <= DSED1)) then
    velo = 0.0_hp
  elseif ((DSED1 < sdiam) .and. (sdiam <= DSED100)) then
    velo = sdiam ** 2.0
    velo = (1.0_hp / 18.0_hp) * (gsg / vis) * velo
  elseif ((DSED100 < sdiam) .and. (sdiam < DSED1000)) then
    velo = 1.0_hp / vis
    velo = velo * velo * (sdiam ** 3.0)
    velo = 0.01_hp * gsg * velo
    velo = sqrt(1.0_hp + velo) - 1.0_hp
    velo = (10.0_hp * vis / sdiam) * velo
  elseif (DSED1000 <= sdiam) then
    velo = 1.1_hp * sqrt(gsg * sdiam)
  else
    velo = -1.0_hp
  endif

  svel0 = velo

  if (velo < 0) then
    write(*, '(a, f15.10)') 'ERROR in "svel0": illegal sediment diameter supplied -> ', sdiam
    stop
  endif

end function svel0


!!!========================================!!!
!!! In the future we might need to involve move advanced
!!! formulations for the calculation of the Cf grag coefficient.
!!! Keep this function as a place holder for future adjustments.
function cfdrag(wdepth,wdis)

!  use params

  implicit none

!!!--- Return value
  real(fp) :: cfdrag

!!!--- Global variables
  real(fp), intent(in) :: wdepth, wdis

!!!--- Local variables
  real(fp) :: an, gr, q_sk_multi

  ! Manning's coefficient (dimensionless)
  q_sk_multi = r_interpo_nn(Q_Sk_Table(1,:),Q_Sk_Table(2,:),Q_sk_tableEntry,wdis)

  an  = Cmann/q_sk_multi

  ! Gravitational acceleration
  gr = Grav

  ! Cf = g / C^2, where C is Chezy's coefficient
  ! and C ~ h^(1/6) / n where h = water depth and n = Manning's coefficient.
  ! Therefore: Cf = g * n^2 * h^(-1/3)
  cfdrag = gr * (an ** 2.0_fp) * (wdepth **(-1.0_fp / 3.0_fp))

  ! According to ICM code the Cf values range: 0.001 - 0.003
  ! Need to research on these for channel flows. Most likely we need
  ! to implement a more advance methodology for calculating Cf
  
  !if (cfdrag < 0.001) cfdrag = 0.001_fp   ! HU turned off this limit
  !if (cfdrag > 0.003) cfdrag = 0.003_fp

end function cfdrag


!!!========================================!!!
function dstar(sdiam, sgrav)

!  use params

  implicit none

!!!--- Return value
  real(hp) :: dstar   ! units: dimensionless

!!!--- Global variables
  real(fp), intent(in)       :: sdiam   ! meters
  real, optional, intent(in) :: sgrav   ! dimensionless

!!!--- Local variables
  real     :: sg
  real(hp) :: vis, gsg, sdiam_star

!
!!!--- executable statements ---
!
  ! default kinematic viscosity of water (m^2/s)
  ! can be a function of temperature
  vis  = ViscW

  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! default
  endif

  gsg = Grav * ( sg - 1.0_hp )

  ! Shields parameter (dimensionless particle size dstar)
  sdiam_star = 1.0_hp / vis
  sdiam_star = (sdiam_star * sdiam_star * gsg) ** (1.0_hp / 3.0_hp)
  sdiam_star = sdiam * sdiam_star

  dstar = sdiam_star
end function dstar

!!!========================================!!!
function shlds(sdiam_star, vanrijn)

!  use precision

  implicit none

!!!--- Return value
  real(hp) :: shlds   ! units: dimensionless

!!!--- Global variables
  real(hp), intent(in)          :: sdiam_star  ! This is Shields dstar (dimensionless)
  integer, optional, intent(in) :: vanrijn

!!!--- Local variables
  real(hp) :: thcrs
  integer  :: use_vanrijn

!
!!!--- executable statements ---
!
  if (present(vanrijn)) then
      use_vanrijn = vanrijn
  else
      use_vanrijn = 1   ! default
  endif

  !!!===== Shields critical shear stress parameter =====!!!
  ! "thcrs" is the Shields parameter determined from the Shields diagram
  ! and computationally can be estimated using eithshieldser the Soulsby-Whitehouse equation
  ! ( R. Soulsby. Dynamics of Marine Sands: A Manual for Practical Applications.
  !   Thomas Telford Publications, London, 1997. Pages 249, ISBN 0-72-772584-X.)
  ! or the van Rijn equations
  ! ( L. C. van Rijn. Unified View of Sediment Transport by Currents and Waves. I: Initiation of Motion,
  ! Bed Roughness, and Bed-Load Transport. Journal of Hydraulic Engineering, ASCE, 133(6):649-667, 2007.

  select case(use_vanrijn)
    case( : 0)
      !!! Use the Soulsby-Whitehouse equation
      thcrs = 0.055_hp * ( 1.0_hp - exp(-0.02_hp * sdiam_star) )
      thcrs = ( 0.30_hp / (1.0_hp + 1.2_hp * sdiam_star) ) + thcrs
    case(1 : )
      !!! Use the van Rijn equations
      if (sdiam_star < 4.0) then
        thcrs = 0.115_hp * ( sdiam_star ** (-0.5_hp) )
      elseif ((4.0 <= sdiam_star) .and. (sdiam_star < 10.0)) then
        thcrs = 0.14_hp * ( sdiam_star ** (-0.64_hp) )
      elseif ((10.0 <= sdiam_star) .and. (sdiam_star < 20.0)) then
        thcrs = 0.04_hp * ( sdiam_star ** (-0.10_hp) )
      elseif ((20.0 <= sdiam_star) .and. (sdiam_star < 150.0)) then
        thcrs = 0.013_hp * ( sdiam_star ** (0.29_hp) )
      elseif (150.0 <= sdiam_star) then
        thcrs = 0.055_hp
      endif
  end select

  shlds = thcrs

end function shlds


!!!========================================!!!
function taucrs(sdiam, vanrijn, sgrav)

!  use params

  implicit none

!!!--- Return value
  real(hp) :: taucrs   ! units: kg / m*s^2

!!!--- Global variables
  real(fp), intent(in)       :: sdiam
  integer, intent(in)        :: vanrijn
  real, optional, intent(in) :: sgrav

!!!--- Local variables
  real     :: sg
  real(hp) :: vis, gsg
  real(hp) :: sdiam_dstar, thcrs

!
!!!--- executable statements ---
!
  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! defaultgfortran precision.f90 params.f90 mod_arrays.f90 sed_properties.f90 test.f90 -o test_sed
  endif

  gsg = Grav * ( sg - 1.0_hp )

  ! Shields parameter (dimensionless particle size dstar)
  sdiam_dstar = dstar(sdiam, sgrav)

  ! Shields critical shear stress parameter
  thcrs = shlds(sdiam_dstar, vanrijn)

  ! Calculate the bed critical shear stress for the initiation
  ! of the motion of cohesionless sediment particles (tcr,0 in van Rijn)
  taucrs = thcrs * DensW * gsg * sdiam

end function taucrs


!!!========================================!!!
function qrijn07(sdiam, wdepth, bedflowvel, sgrav)

!  use params

  implicit none

!!!--- Return value
  real(hp) :: qrijn07   ! units: kg/m/s

!!!--- Global variables
  real(fp), intent(in) :: sdiam
  ! "bedflowvel" is the effective near bed flow velocity a function
  ! of wind, current and wave velocities. It is calculated outside this function.
  real(fp), intent(in) :: wdepth, bedflowvel

  real, optional, intent(in) :: sgrav

!!!--- Local variables
  real     :: sg
  real(fp) :: gsg, ds90, sdiam_star
  real(fp) :: Me, sMe, Ue, Ucr, qs
  real(fp) :: DSED0, DSED1, DSED2

!
!!!--- executable statements ---
!
  DSED0 = 0.00005_fp  ! meters
  DSED1 = 0.0005_fp   ! meters
  DSED2 = 0.002_fp    ! meters

  if ((sdiam <= DSED0) .or. (sdiam >= DSED2)) then
    write(*, '(a, f15.10)') 'ERROR in "qrijn07": illegal sediment diameter supplied -> ', sdiam
    write(*, '(a)')         '                    valid range 0.00005 m < d < 0.002 m'
    stop
  endif

  if (D90x <= 0.0) then
    write(*, '(a, f8.5)') 'ERROR in "qrijn07": illegal/undefined D90x -> ', D90x
    stop
  endif

  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! default
  endif

  gsg = Grav * ( sg - 1.0_hp )

  ds90 = D90x * sdiam
  Ue = bedflowvel

  !!! The following equations follow after van Rijn:
  !!! Unified View of Sediment Transport by Currents and Waves.
  !!! I: Initiation of Motion, Bed Roughness, and Bed-Load Transport
  !!! Journal of Hydraulic Engineering 2007
  if ((DSED0 < sdiam) .and. (sdiam <= DSED1)) then
    Ucr = 0.19_fp * (sdiam ** 0.1_fp) * log10(12.0_fp * wdepth / (3.0_fp * ds90))
  elseif ((DSED1 < sdiam) .and. (sdiam < DSED2)) then
    Ucr = 8.50_fp * (sdiam ** 0.6_fp) * log10(12.0_fp * wdepth / (3.0_fp * ds90))
  endif

  ! Mobility parameter
  ! Need absolute value since this value is raised to 2.4 power - removes sign of Ue - Ucr
  Me = abs(Ue - Ucr) / sqrt(sdiam * gsg) ! mobility parameter
  sMe = abs(Ue - Ucr) / (Ue - Ucr)

  ! Suspended load under currents only
  sdiam_star = dstar(sdiam, sg)
  qs = sMe * (Me ** 2.4)* (sdiam_star ** (-0.6))
  qs = alphaSED * sg * DensW * bedflowvel * sdiam * qs

  ! If Ue <= Ucr then there is no resuspension
  if (qs <= 0.0) qs = 0.0
   ! print*, 'Ucr = ', Ucr, Ue
  qrijn07 = qs   ! in kg/m/s

end function qrijn07


!!!========================================!!!
function srcsand(sdiam, wdepth, wwidth, bedflowvel, sconc, tsconc, sgrav)
! This function calculates the source/sink term in the
! sediment transport equation: R = Erosion - Deposition
! for sand size sediments according to van Rijn (2007).

!  use params

  implicit none

!!!--- Return value
  real(fp) :: srcsand   ! units: g/m/s : it is per unit length

!!!--- Global variables
  real(fp), intent(in) :: sdiam
  ! "bedflowvel" is the effective near bed flow velocity a function
  ! of wind, current and wave velocities. It is calculated outside this function.
  ! sconc  : the sediment concentration for the particular sediment size
  ! tsconc : total sediment concentartion (all particle sizes in suspension)
  real(fp), intent(in) :: wdepth, wwidth, bedflowvel, sconc, tsconc   ! m, m, m/s, mg/L = g/m^3, mg/L = g/m^3

  real, optional, intent(in) :: sgrav

!!!--- Local variables
  real     :: sg, ff_hu
  real(fp) :: ws, tau_bed, tau_crs
  real(fp) :: depflux, eroflux, readapt, nexp
  real(fp) :: dsand, alpha, Cgels, Ffloc, Fhs

!
!!!--- executable statements ---
!
  dsand = 0.000062

  if (sdiam < dsand) then
    write(*, '(a, f15.10)') 'ERROR in "srcsand": illegal sediment diameter supplied -> ', sdiam
    write(*, '(a)')         '                    valid range d >= 0.000062 m'
    stop
  endif

  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! default
  endif

  !!!--- Settling velocity of sediment particles ---
  ! Adjust the settling velocity for flocculation and hindered settling effects
  ! ws = Ffloc * Fhs * ws
  Cgels = 0.65_fp * sg * DensW * 1000.0_fp ! in g/m^3
  !--- First flocculation effects
  alpha = dsand / sdiam - 1.0_fp
  if (alpha <= 0.0) then
    alpha = 0.0
    Ffloc = 1.0
  else
    if (alpha >= 3.0) alpha = 3.0
    Ffloc = 4.0_fp + log10( 2.0_fp * (dsand / sdiam) * (sconc / Cgels) )
    Ffloc = Ffloc ** alpha
  endif
  !--- Next the hindered settling effects
  Fhs = 1.0_fp - 0.65_fp * (tsconc / Cgels) ** 5.0_fp

  ! Settling velocity in clear water
  ws = svel0(sdiam, sgrav)
  ! Adjusted settling velocity
  !ws = Ffloc * Fhs * ws

  ! Deposition of suspended sediments
  depflux = ws * sconc * wwidth

  ! Erosion of bottom sediments according the van Rijn formulation (units: kg/m/s)
  if( useBOTEROS ) then
    eroflux = qrijn07(sdiam, wdepth, bedflowvel, sg)
    eroflux = 1000.0_fp * eroflux   ! in g/m/s
  else
    eroflux = 0.0_fp
  endif

  ! positive value of srcsand means that sediments enter the water column
  
  !if(sconc .lt. 30.)ff_hu=1.
  !if(sconc .ge. 30. .and. sconc .lt. 60.)ff_hu=1.-(sconc-30.)/(60.-30.)*(1.-0.25)
  !if(sconc .gt. 60.)ff_hu=0.25
  if(sconc .lt. 20.)ff_hu=1.
  if(sconc .ge. 20. .and. sconc .lt. 60.)ff_hu=1.-(sconc-20.)/(60.-20.)*(1.-0.2)
  if(sconc .gt. 60.)ff_hu=0.2
  
  srcsand = eroflux * ff_hu - depflux

end function srcsand


!!!========================================!!!
function srcchsv(sdiam, wdepth, wwidth, wdis, bedflowvel, sconc, tsconc, userCriticalShear, sgrav)
! This function calculates the source/sink term in the
! sediment transport equation: R = Erosion - Deposition
! according to Krone.

!  use params

  implicit none

!!!--- Return value
  real(fp) :: srcchsv   ! units: g/m/s : it is per unit length

!!!--- Global variables
  real(fp), intent(in) :: sdiam
  ! "bedflowvel" is the effective near bed flow velocity a function
  ! of wind, current and wave velocities. It is calculated outside this function.
  ! sconc  : the sediment concentration for the particular sediment size
  ! tsconc : total sediment concentartion (all particle sizes in suspension)
  real(fp), intent(in) :: wdepth, wwidth, wdis, bedflowvel, sconc, tsconc   ! m, m, m3/s, m/s, mg/L = g/m^3, mg/L = g/m^3

  real(fp), intent(in) :: userCriticalShear   !! Change Nazmul

  real, optional, intent(in) :: sgrav

!!!--- Local variables
  real     :: sg
  real(fp) :: ws, tau_bed, tau_crs
  real(fp) :: depflux, eroflux, readapt, nexp
  real(fp) :: dsand, alpha, Cgels, Ffloc, Fhs

!
!!!--- executable statements ---
!
  dsand = 0.000062

  if (sdiam >= dsand) then
    write(*, '(a, f15.10)') 'ERROR in "srcchsv": illegal sediment diameter supplied -> ', sdiam
    write(*, '(a)')         '                    valid range d < 0.000062 m'
    stop
  endif

  !!! specific gravity of sediments (dimensionless)
  if (present(sgrav)) then
      sg = sgrav
  else
      sg = SpSed   ! default
  endif

  readapt = SEDcalib * 1000.0_fp   ! readaptability coefficient; convert it to g/m^2/s from kg/m^2/s
  nexp = SEDn
  if (nexp < 1.0) nexp = 1.0
  if (nexp > 3.0) nexp = 3.0

  !!!--- Settling velocity of sediment particles ---
  ! Adjust the settling velocity for flocculation and hindered settling effects
  ! ws = Ffloc * Fhs * ws
  Cgels = 0.65_fp * sg * DensW * 1000.0_fp ! in g/m^3
  !--- First flocculation effects
  alpha = dsand / sdiam - 1.0_fp
  if (alpha <= 0.0) then
    alpha = 0.0
    Ffloc = 1.0
  else
    if (alpha >= 3.0) alpha = 3.0
    Ffloc = 4.0_fp + log10( 2.0_fp * (dsand / sdiam) * (sconc / Cgels) )
    Ffloc = Ffloc ** alpha
  endif
  !--- Next the hindered settling effects
  Fhs = 1.0_fp - 0.65_fp * (tsconc / Cgels) ** 5.0_fp

  ! Settling velocity in clear water
  !ws = svel0(sdiam, sg)
  ws=0.05e-3
  ! Adjusted settling velocity
  !ws = Ffloc * Fhs * ws

  ! Calculate bed shear stress
  tau_bed = cfdrag(wdepth,wdis) * DensW * (abs(bedflowvel) ** 2.0)
  !print*, "fine tau_bed and userCriticalShear", tau_bed, userCriticalShear
  ! Calculate bed critical shear stress for initiation of sediment motion
  !! tau_crs = taucrs(sdiam, 1, sg)  !! !! Change Nazmul : tau_crit is taken from user defined value to use as a calibration parameter
  tau_crs = userCriticalShear
  
  depflux=ws * sconc * wwidth

  if( useBOTEROS ) then

    if (tau_bed < tau_crs) then
      eroflux = 0.
    else
      eroflux = readapt * ( (tau_bed / tau_crs - 1.0_fp) ** nexp ) * wwidth
    endif
  else
     eroflux = 0. 
  endif
  
  ! positive value of srcchsv means that sediments enter the water column
  srcchsv = eroflux - depflux

end function srcchsv

function r_interpo_nn(x,y,jj,xt)

    integer, intent(in) :: jj
    real(fp), intent(in) :: xt
    real, intent(in) :: x(jj), y(jj)

    real :: yt

    ! nn means nearest neighbour

    if (xt.le. x(1)) then
        yt=y(1)
    elseif (xt.ge. x(jj)) then
        yt=y(jj)
    else
        do j=1,jj-1
            if((x(j)-xt)*(x(j+1)-xt).le.0)then

                yt=(xt-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j))+y(j)

                EXIT
            endif
        end do
    end if
    r_interpo_nn = yt
    return
end function

! add temperature functions
! f(U10)
function wind_temp_coef(wind10)

!  use params

  implicit none

!!!--- Return value
  real(hp) :: wind_temp_coef

!!!--- Global variables
  real(fp), intent(in)       :: wind10

!
!!!--- executable statements ---
!
 wind_temp_coef=(3.5+2.0*wind10)*(5.e6/Total_area)**0.05

end function wind_temp_coef


function srctemp(wind10,temp_back,temp_water,wwidth)

!  use params

  implicit none

!!!--- Return value
  real(fp) :: srctemp

!!!--- Global variables
  real(fp), intent(in)       :: wind10, temp_back, temp_water, wwidth

!!!--- Local variables
  real(fp) :: coef
!
!!!--- executable statements ---
!
  coef=4.48+0.049*temp_water+wind_temp_coef(wind10)*(1.12+0.018*temp_water+0.00158*temp_water**2)
  srctemp=-1.*coef*(temp_water-temp_back)*wwidth/DensW/3930.  ! Cp=3930 specific sea water capacity

end function srctemp						   

end module sed_properties_TMP_R08
