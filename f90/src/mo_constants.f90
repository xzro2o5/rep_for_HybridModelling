MODULE constants

!- Description:
!
!  This module contains basic and derived constants
!
!  Written Jan 2011, Matthias Cuntz - Ported C-Code

  USE kinds, ONLY: i4, wp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: version       = '4.0'        ! model version number
  CHARACTER(len=*), PARAMETER :: main_file     = 'canveg.f90' ! main file name
  CHARACTER(len=*), PARAMETER :: namelist_file = 'canveg.nml' ! namelist file name

  ! maximum model parameters
  INTEGER(i4), PARAMETER :: nsoilmax = 10       ! max # of soil layers
  INTEGER(i4), PARAMETER :: nbetamax = 5        ! max # of levels for beta distribution of e.g. lai
  INTEGER(i4), PARAMETER :: nleafopticalmax = 4 ! max # of leaf optical properties

  REAL(wp), PARAMETER :: undef = 9999._wp      ! undefinde values

  ! file units
  INTEGER, PARAMETER :: nin            = 5   ! standard input stream
  INTEGER, PARAMETER :: nout           = 6   ! standard output stream
  INTEGER, PARAMETER :: nerr           = 0   ! error output stream
  INTEGER, PARAMETER :: ninnml         = 100 ! namelist unit

  INTEGER, PARAMETER :: ninmet         = 101 ! fptr1   met in
  INTEGER, PARAMETER :: nindisp        = 104 ! fptr4   dispersion matrix
  INTEGER, PARAMETER :: ninlai         = 115 ! fptr15  lai in

  INTEGER, PARAMETER :: noutseas       = 106 ! fptr6   season out
  INTEGER, PARAMETER :: noutsoil       = 107 ! fptr7   soil out
  INTEGER, PARAMETER :: noutdaily      = 108 ! fptr8   daily out
  INTEGER, PARAMETER :: noutprof       = 109 ! fptr9   profile out
  INTEGER, PARAMETER :: noutflux       = 110 ! fptr10  flux profile out
  INTEGER, PARAMETER :: noutopti       = 111 ! fptr11  optimise
  INTEGER, PARAMETER :: noutdum        = 112 ! fptr12  dummy out
  INTEGER, PARAMETER :: nouth2osoil    = 113 ! fptr13  h2o soil out

  INTEGER, PARAMETER :: ninwiso        = 114 ! fptr14  water iso in
  INTEGER, PARAMETER :: noutwisoleaf   = 116 ! fptr16  wiso h2o leaf out
  INTEGER, PARAMETER :: noutwisoprof   = 117 ! fptr17  wiso profile out
  INTEGER, PARAMETER :: noutwisoflux   = 118 ! fptr18  wiso flux profile out
  INTEGER, PARAMETER :: noutwisosoil   = 119 ! fptr19  wiso h2o soil out

  INTEGER, PARAMETER :: noutcisodaily  = 120 ! fptr20  13c daily out
  INTEGER, PARAMETER :: noutcisoseason = 121 ! fptr21  13c season out
  INTEGER, PARAMETER :: noutcisoprof   = 122 ! fptr22  13c profile air out
  INTEGER, PARAMETER :: noutcisoflux   = 123 ! fptr23  13c profile flux out

  CHARACTER(len=*), PARAMETER :: dailyfile          = "daily_ave"
  CHARACTER(len=*), PARAMETER :: daily13cfile       = "daily_ave_13c"
  CHARACTER(len=*), PARAMETER :: hourlyfile         = "season"
  CHARACTER(len=*), PARAMETER :: hourly13cfile      = "season_13c"
  CHARACTER(len=*), PARAMETER :: optimisefile       = "optimise"
  CHARACTER(len=*), PARAMETER :: profilefile        = "profile_air"
  CHARACTER(len=*), PARAMETER :: profile13cfile     = "profile_air_13c"
  CHARACTER(len=*), PARAMETER :: profileisofile     = "profile_air_wiso"
  CHARACTER(len=*), PARAMETER :: fluxprofilefile    = "profile_fluxes"
  CHARACTER(len=*), PARAMETER :: fluxprofile13cfile = "profile_fluxes_13c"
  CHARACTER(len=*), PARAMETER :: fluxprofileisofile = "profile_fluxes_wiso"
  CHARACTER(len=*), PARAMETER :: soilfile           = "soil"
  CHARACTER(len=*), PARAMETER :: h2osoilfile        = "h2osoil"
  CHARACTER(len=*), PARAMETER :: h2osoilisofile     = "h2osoil_wiso"
  CHARACTER(len=*), PARAMETER :: h2oleafisofile     = "h2oleaf_wiso"

  ! computational
  REAL(wp), PARAMETER :: zero     = 0.0_wp      ! 0
  REAL(wp), PARAMETER :: half     = 0.5_wp      ! 1/2
  REAL(wp), PARAMETER :: twothird = 2._wp/3._wp ! 1/2
  REAL(wp), PARAMETER :: one      = 1.0_wp      ! 1
  REAL(wp), PARAMETER :: two      = 2.0_wp      ! 2
  REAL(wp), PARAMETER :: e1       = 1e-1_wp     ! 1e-1
  REAL(wp), PARAMETER :: e2       = 1e-2_wp     ! 1e-2
  REAL(wp), PARAMETER :: e3       = 1e-3_wp     ! 1e-3
  REAL(wp), PARAMETER :: e6       = 1e-6_wp     ! 1e-6
  REAL(wp), PARAMETER :: e12      = 1e-12_wp    ! 1e-12
  REAL(wp), PARAMETER :: e15      = 1e-15_wp    ! 1e-15
  REAL(wp), PARAMETER :: s2h      = one / 3600._wp ! s -> h

  ! mathematical
  REAL(wp), PARAMETER :: pi    = 3.14159265358979323846_wp ! pi
  REAL(wp), PARAMETER :: pi2   = 2.0_wp * pi               ! 2*pi
  REAL(wp), PARAMETER :: pi4   = 4.0_wp * pi               ! 4*pi
  REAL(wp), PARAMETER :: pi9   = 9.0_wp / pi               ! 9/pi
  REAL(wp), PARAMETER :: pi180 = pi / 180.0_wp             ! pi/180
  ! REAL(wp), PARAMETER :: pi    = 3.1415926536_wp
  ! REAL(wp), PARAMETER :: pi2   = 6.2831853072_wp
  ! REAL(wp), PARAMETER :: pi4   = 12.5663706_wp
  ! REAL(wp), PARAMETER :: pi180 = 0.0174532925_wp
  ! REAL(wp), PARAMETER :: pi9   = 2.86478898_wp

  ! physical
  REAL(wp), PARAMETER :: TN0          = 273.15_wp    ! Celcius <-> Kelvin [K]
  REAL(wp), PARAMETER :: P0           = 101325._wp   ! Standard pressure [Pa]
  REAL(wp), PARAMETER :: tk_25        = TN0 + 25._wp ! 25 Celcius <-> Kelvin [K]
  REAL(wp), PARAMETER :: Gravity      = 9.81_wp      ! Gravity acceleration [m2/s]
  REAL(wp), PARAMETER :: vonKarman    = 0.41_wp      ! von Karman constant
  REAL(wp), PARAMETER :: waterdensity = 1000._wp     ! density of water [kg m-3]
  REAL(wp), PARAMETER :: rugc         = 8.3144_wp    ! Universal gas constant [J mole-1 K-1]
  REAL(wp), PARAMETER :: rgc1000      = 8314.4_wp    ! gas constant times 1000. [mJ mole-1 K-1]
  REAL(wp), PARAMETER :: Rw           = 461.5_wp     ! gas constant for water vapour [J kg-1 K-1] (was 461.89 in v2.0)

  ! Constants for leaf energy balance
  REAL(wp), PARAMETER ::  sigma      = 5.67e-08_wp ! Stefan-Boltzmann constant W M-2 K-4
  REAL(wp), PARAMETER ::  cp         = 1005._wp    ! Specific heat of air, J KG-1 K-1
  REAL(wp), PARAMETER ::  mass_air   = 28.97_wp    ! Molecular mass of dry air, g mole-1
  REAL(wp), PARAMETER ::  mass_CO2   = 44._wp      ! molecular mass of CO2, g mole-1
  REAL(wp), PARAMETER ::  mass_13CO2 = 45._wp      ! molecular mass of 13CO2, g mole-1
  REAL(wp), PARAMETER ::  dldt       = -2370._wp   ! Derivative of the latent heat of vaporization

  ! constants for the polynomial equation for saturation vapor pressure-T function, es=f(t)
  REAL(wp), PARAMETER :: a1en = 617.4_wp     !
  REAL(wp), PARAMETER :: a2en = 42.22_wp     !
  REAL(wp), PARAMETER :: a3en = 1.675_wp     !
  REAL(wp), PARAMETER :: a4en = 0.01408_wp   !
  REAL(wp), PARAMETER :: a5en = 0.0005818_wp !

  ! Diffusivity values for 273 K and 1013 mb (STP) using values from Massman (1998) Atmos Environment
  ! These values are for diffusion in air. When used these values must be adjusted for
  ! temperature and pressure
  ! nu, Molecular viscosity
  REAL(wp), PARAMETER :: nuvisc = 13.27_wp ! mm2 s-1
  REAL(wp), PARAMETER :: nnu = 0.00001327_wp ! m2 s-1
  ! Diffusivity of CO2
  REAL(wp), PARAMETER :: dc = 13.81_wp ! mm2 s-1
  REAL(wp), PARAMETER :: ddc = 0.00001381_wp ! m2 s-1
  ! Diffusivity of heat
  REAL(wp), PARAMETER :: dh = 18.69_wp ! mm2 s-1
  REAL(wp), PARAMETER :: ddh = 0.00001869_wp ! m2 s-1
  ! Diffusivity of water vapor
  REAL(wp), PARAMETER :: dv = 21.78_wp ! mm2 s-1
  REAL(wp), PARAMETER :: ddv = 0.00002178_wp ! m2 s-1
  ! Diffusivity of ozone
  REAL(wp), PARAMETER :: do3 = 14.44_wp ! mm2 s-1
  REAL(wp), PARAMETER :: ddo3 = 0.00001444_wp ! m2 s-1

  ! 13CO2
  ! Isotope ratio of PeeDee Belimdite standard (PDB)
  ! redefined as 13C/(12C+13C) for PDB, Tans et al 93
  REAL(wp), PARAMETER :: Rpdb_CO2 = 0.01115_wp
  ! Isotope ratio of PeeDee Belimdite standard (PDB)
  ! defined as 13C/12C, from Farquhar
  REAL(wp), PARAMETER :: Rpdb_12C = 0.01124_wp

  ! WISO
  ! Vienna Standard Mean Ocean Water (VSMOW) for water isotopes
  REAL(wp), PARAMETER :: RVSMOW_18O = 2005.2e-6_wp  ! 18O/16O
  REAL(wp), PARAMETER :: RVSMOW_D   = 155.76e-6_wp  ! 2H/1H

  ! numerical
  REAL(wp), PARAMETER :: isotolerance = 1e-15_wp ! tolerance in isotope variables so that two values are the same
  REAL(wp), PARAMETER :: isnight      = 0.05_wp  ! if (solar%sine_beta <= isnight) then it is dark

END MODULE constants
