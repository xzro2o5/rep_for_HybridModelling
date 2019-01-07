MODULE parameters

  ! Stores parameters and namelist variables.
  ! Provides routines to map namelist inputs to parameters

  ! Written Jan 2011, Matthias Cuntz - Ported C-Code

  USE kinds,       ONLY: wp, i4, i8
  USE constants,   ONLY: nbetamax, nsoilmax, nleafopticalmax, & ! max dims
                         zero, half, twothird, one, &           ! numeric
                         ninnml, namelist_file, &                  ! namelist
                         rugc, sigma, TN0                       ! physical
  USE setup,       ONLY: ncl, ntl, nsky, nl, nsoil, nbeta, ndaysc13, nwiso, nlop

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_namelist
  PUBLIC :: zd, z0, izref, delz, zh65, epsigma, epsigma2, epsigma4, epsigma6, epsigma8, epsigma12, qalpha2, &
       vc25, jm_vc, rd_vc, g0, a1, D0, kball, bprime, vcopt, jmopt, ht_midpt, lai_freq, pai_freq, &
       par_reflect, par_trans, par_soil_refl_dry, nir_reflect, nir_trans, nir_soil_refl_dry, workdir, &
       indir, metinfile, laiinfile, wisoinfile, outdir, outsuffix, dispfile, netcdf_in, netcdf_out, netcdf_disp, &
       start_run, end_run, start_profiles, end_profiles, latitude, longitude, perc_up, perc_dn, zone, ht, pai, lai, ustar_ref, &
       zm, hkin, skin, ejm, evc, kc25, ko25, o2, tau25, ekc, eko, erd, ektau, toptvc, &
       toptjm, curvature, qalpha, gm_vc, rsm, brs, ep, n_stomata_sides, betfact, markov, lleaf, leaf_out, leaf_full, &
       leaf_fall, leaf_fall_complete, attfac, eabole, epsoil, water_film_thickness, tau_water, extra_nate, nup ! add perc_up and down for sensitivity analysis Yuan 2017.10.26

  ! Derived parameters
  INTEGER(i4) :: izref     ! array value of reference height = measurement height*ncl/ht
  REAL(wp)    :: zd        ! displacement height [m]
  REAL(wp)    :: z0        ! rougness lenght [m]
  REAL(wp)    :: delz      ! height of each layer, ht/ncl [m]
  REAL(wp)    :: zh65      ! 0.65/ht [m]
  REAL(wp)    :: epsigma   ! ep*sigma
  REAL(wp)    :: epsigma2  ! 2*ep*sigma
  REAL(wp)    :: epsigma4  ! 4.0 * ep * sigma
  REAL(wp)    :: epsigma6  ! 6.0 * ep * sigma
  REAL(wp)    :: epsigma8  ! 8.0 * ep * sigma
  REAL(wp)    :: epsigma12 ! 12.0 * ep * sigma
  REAL(wp)    :: qalpha2   ! qalpha**2

  ! Parameter arrays
  REAL(wp), DIMENSION(:), ALLOCATABLE :: vc25   ! carboxylation rate at 25 deg C, [umol m-2 s-1]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: jm_vc  ! jmax/vcmax, jmax = electron transport rate at 25 deg C, [umol m-2 s-1]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: rd_vc  ! rd/vcmax, rd = dark respiration at 25 C, [umol m-2 s-1]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: g0     ! g0, intercept of Leuning stomata model for water vapor [mol(H2O) m-2 s-1]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: a1     ! a1, slope of Leuning stomata model for water vapor [-]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: D0     ! empirical coefficient for humidity sensitivity of Leuning stomata model [Pa]
  REAL(wp), DIMENSION(:), ALLOCATABLE :: kball  ! Ball-Berry stomatal coefficient for stomatal conductance
  REAL(wp), DIMENSION(:), ALLOCATABLE :: bprime ! intercept of Ball-Berry model, mol m-2 s-1
  REAL(wp), DIMENSION(:), ALLOCATABLE :: vcopt  ! maximum rate of RuBP Carboxylase/oxygenase
  REAL(wp), DIMENSION(:), ALLOCATABLE :: jmopt  ! optimal rate of electron transport
  REAL(wp), DIMENSION(:), ALLOCATABLE :: ht_midpt ! height of mid point of layer of lai_freq
  REAL(wp), DIMENSION(:), ALLOCATABLE :: lai_freq ! fraction of total LAI per layer
  REAL(wp), DIMENSION(:), ALLOCATABLE :: pai_freq ! fraction of total PAI per layer
  REAL(wp), DIMENSION(:), ALLOCATABLE :: par_reflect       ! PAR leaf reflectivity, BARK REFLECTIVITY, AVG TOP AND BOTTOM
  REAL(wp), DIMENSION(:), ALLOCATABLE :: par_trans         ! PAR leaf transmissivity, LEAF TRANSMISSIVITY
  REAL(wp), DIMENSION(:), ALLOCATABLE :: par_soil_refl_dry ! PAR soil reflectivitiy
  REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_reflect       ! NIR leaf reflectivity
  REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_trans         ! NIR leaf transmissivity
  REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_soil_refl_dry ! NIR soil reflectivitiy

  ! Namelist parameters
  CHARACTER(LEN=256) :: workdir           ! working directory
  CHARACTER(LEN=256) :: indir             ! location of input file
  CHARACTER(LEN=256) :: metinfile         ! met input file
  CHARACTER(LEN=256) :: laiinfile         ! lai input file
  CHARACTER(LEN=256) :: wisoinfile        ! wiso input fie
  CHARACTER(LEN=256) :: outdir            ! location of output file
  CHARACTER(LEN=256) :: outsuffix         ! output suffix
  CHARACTER(LEN=256) :: dispfile          ! location and name of disperion matrix
  INTEGER(i4)        :: netcdf_in         ! 1: input in netcdf; 0: input in ascii format
  INTEGER(i4)        :: netcdf_out        ! 1: output in netcdf; 0: output in ascii format
  INTEGER(i4)        :: netcdf_disp       ! 1: dispersion matrix in netcdf; 0: dispersion in ascii format
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  INTEGER(i4)        :: year0             ! Year
  INTEGER(i8)        :: start_run         ! start of model run, [dddhhmm] > 0010100
  INTEGER(i8)        :: end_run           ! end of model run, [dddhhmm] < 3662400
  INTEGER(i8)        :: start_profiles    ! start of profile output, [dddhhmm]
  INTEGER(i8)        :: end_profiles      ! end of profile output, [dddhhmm]
  REAL(wp)           :: time_step         ! time step between consecutive computation times
  !  ! sets soil respiration reference temperature to
  !   (1) input soil temperature or
  !   (0) modeled soil temperature at 5 cm
  INTEGER(i4)        :: switch_soil_resp_temp
  ! set coupled stomata-photosynthesis model to
  !   (0) Ball Berry (Baldocchi analytical solution) or
  !   (1) Leuning and mesophyll conductance (Knohl analytical solution)
  INTEGER(i4)        :: switch_ball       !
  INTEGER(i4)        :: switch_isoprene   ! (1) calc isoprene ! (0) no isoprene
  INTEGER(i4)        :: switch_d13c       ! (1) calc d13C ! (0) no d13C
  INTEGER(i4)        :: switch_wiso       ! (1) calc water isotopes ! (0) no water isotopes
  ! How to determine autotrophic respiration
  !   (0) total = 50% auto + 50% hetero
  !   (1) from BETHY (Knorr 1997)
  INTEGER(i4)        :: switch_bethy_resp !
  ! Debug by not allowing condensation
  INTEGER(i4)        :: switch_no_negative_water_flux ! (1) all water fluxes >= 0 ! (0) no restriction, <0 possible
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Ecosystem Set-up --- Basic Parameter
  ! Site location = Mesita del Buey, NM, USA
  REAL(wp)           :: latitude          ! latitude  N
  REAL(wp)           :: longitude         ! longitude E
  REAL(wp)           :: perc_up           ! range of parameter variation for sensitivity analysis Yuan 2017.10.26
  REAL(wp)           :: perc_dn           ! range of parameter variation for sensitivity analysis Yuan 2017.10.26

  ! Eastern Standard Time
  REAL(wp)           :: zone              ! delay from GMT
  REAL(wp)           :: ht                ! Canopy height [m]
  REAL(wp)           :: pai               ! Plant area index [m2 m-2]
  REAL(wp)           :: lai               ! Maximum leaf area index [m2 m-2]
  REAL(wp)           :: vc25_in              ! carboxylation rate at 25 deg C, [umol m-2 s-1], old = 66.3
  ! ratios of x to vcmax
  !   jmax = electron transport rate at 25 deg C, [umol m-2 s-1] old = 127.5
  REAL(wp)           :: jm_vc_in             !
  !  ratio of rd to vcmax, rd = dark respiration at 25 C, [umol m-2 s-1]
  REAL(wp)           :: rd_vc_in             !

  REAL(wp)           :: ustar_ref         ! reference ustar for Dij [umol m-2 s-1]
  !INTEGER(i4) :: ncl   ! canopy layers
  !INTEGER(i4) :: nsky  ! # of sky angle classes
  !INTEGER(i4) :: nsoil ! # of soil layers (0: old Canoak formulation)
  !INTEGER(i4) :: nbeta ! # of levels for beta distribution of e.g. lai
  !INTEGER(i4) :: ndaysc13 ! # of days to remember for mean 13C discrimination
  !INTEGER(i4) :: nwiso    ! # of number of (isotopic) waters
  REAL(wp)           :: zm                ! measurement height
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Ecosystem Set-up --- Photosynthesis Parameter
  ! Taken from Harley and Baldocchi (1995, PCE)
  REAL(wp)           :: hkin              ! hkin, enthalpy term, old 200000
  REAL(wp)           :: skin              ! skin, entropy term
  REAL(wp)           :: ejm               ! ejm, activation energy for electron transport, old 55000
  REAL(wp)           :: evc               ! evc, activation energy for carboxylation, old 55000
  !  Enzyme constants & partial pressure of O2 and CO2
  !  Michaelis-Menten K values. From survey of literature.
  REAL(wp)           :: kc25              ! kc25, ,kinetic coef for CO2 at 25 C,   microbars old = 274.6
  REAL(wp)           :: ko25              ! ko25, kinetic coef for O2 at 25C,  millibars old = 419.8
  REAL(wp)           :: o2                ! o2, oxygen concentration  umol mol-1
  ! tau is computed on the basis of the Specificity factor (102.33)
  ! times Kco2/Kh2o (28.38) to convert for value in solution to that based in air/
  ! The old value was 2321.1. New value for Quercus robor from Balaguer et al. (1996)
  REAL(wp)           :: tau25             ! tau25, tau coefficient
  ! Arrhenius constants
  ! Eact for Michaelis-Menten const. for KC, KO and dark respiration
  ! These values are from Harley
  REAL(wp)           :: ekc               ! ekc, Activation energy for K of CO2 ! J mol-1, old = 80500.0
  REAL(wp)           :: eko               ! eko, Activation energy for K of O2, old = 14500.0
  REAL(wp)           :: erd               ! erd, activation energy for dark respiration, eg Q10=2
  REAL(wp)           :: ektau             ! ektau, old = -29000
  REAL(wp)           :: toptvc            ! toptvc, optimum temperature for maximum carboxylation
  REAL(wp)           :: toptjm            ! toptjm, optimum temperature for maximum electron transport
  ! curvature for light response function (ranges from 0 to 1)
  REAL(wp)           :: curvature
  ! leaf quantum yield, electrons
  REAL(wp)           :: qalpha
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Ecosystem Set-up --- Stomata Parameter
  REAL(wp)           :: g0_in ! g0, intercept of Leuning stomata model for water vapor [mol(H2O) m-2 s-1]
  REAL(wp)           :: a1_in ! a1, slope of Leuning stomata model for water vapor [-]
  REAL(wp)           :: D0_in ! Do, empirical coefficient for humidity sensitivity of Leuning stomata model [Pa]
  REAL(wp)           :: gm_vc ! ratio of mesophyll conductance to Vcmax [mol m-2 s-1]
  REAL(wp)           :: fthreshold ! threshold of plant available water below which stomata response to water stress
  REAL(wp)           :: fslope ! slope of stomata response to water stress (should be: slope = 1/threshold)
  ! Ball-Berry stomatal coefficient for stomatal conductance
  REAL(wp)           :: kball_in
  ! intercept of Ball-Berry model, mol m-2 s-1
  REAL(wp)           :: bprime_in ! intercept for H2O 0.0175
  ! Minimum stomatal resistance, s m-1.
  REAL(wp)           :: rsm
  REAL(wp)           :: brs ! curvature coeffient for light response
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Ecosystem Set-up --- Leaf Parameter
  REAL(wp)           :: ep ! emissivity of leaves
  ! number of leaf sides with stomata (1) hypostomatous, (2) amphistomatous
  INTEGER(i4)        :: n_stomata_sides
  ! multiplication factor for aerodynamic sheltering
  !   based on work by Grace and Wilson
  REAL(wp)           :: betfact
  ! leaf clumping factor
  REAL(wp)           :: markov
  ! Leaf dimension. geometric mean of length and width (m)
  REAL(wp)           :: lleaf ! leaf length, m
  ! leaf development
  INTEGER(i4)        :: leaf_out ! start of leaf out
  INTEGER(i4)        :: leaf_full ! day of full leaves
  INTEGER(i4)        :: leaf_fall ! day of leaf fall
  INTEGER(i4)        :: leaf_fall_complete ! day of leaf fall end
  ! Distribution
  ! height of mid point of layer of lai_freq
  REAL(wp), DIMENSION(nbetamax) :: ht_midpt_in ! ht_midpt
  ! fraction of total LAI per layer
  REAL(wp), DIMENSION(nbetamax) :: lai_freq_in ! lai_freq
  ! fraction of total PAI per layer
  REAL(wp), DIMENSION(nbetamax) :: pai_freq_in
  ! optical properties of leaves and soil
  !   1. date <= leafout || date >= leaf_fall_complete
  !   2. date > leafout and date < fulleaf
  !   3. date >= fulleaf and date <= leaf_fall
  !   4. date > leaf_fall and date < leaf_fall_complete
  REAL(wp), DIMENSION(nleafopticalmax) :: par_reflect_in       ! PAR leaf reflectivity, BARK REFLECTIVITY, AVG TOP AND BOTTOM
  REAL(wp), DIMENSION(nleafopticalmax) :: par_trans_in         ! PAR leaf transmissivity, LEAF TRANSMISSIVITY
  REAL(wp), DIMENSION(nleafopticalmax) :: par_soil_refl_dry_in ! PAR soil reflectivitiy
  REAL(wp), DIMENSION(nleafopticalmax) :: nir_reflect_in       ! NIR leaf reflectivity
  REAL(wp), DIMENSION(nleafopticalmax) :: nir_trans_in         ! NIR leaf transmissivity
  REAL(wp), DIMENSION(nleafopticalmax) :: nir_soil_refl_dry_in ! NIR soil reflectivitiy
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Ecosystem Set-up --- Miscellanous Parameter
  REAL(wp)           :: attfac ! attenuation factor for wind speed inside canopy
  REAL(wp)           :: eabole ! activation energy for bole respiration for Q10 = 2.02
  ! Constants for leaf energy balance
  REAL(wp)           :: epsoil ! epsoil, Emissivity of soil
  ! Interception reservoir
  REAL(wp)           :: water_film_thickness ! leaf water film thickness [mm m-2 LAI]
  REAL(wp)           :: tau_water ! rain interception efficiency per m-2 LAI [0-1]
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- 13CO2 ---
  !  d13C of heterotrophic soil respiration (= longterm) [per mil]
  REAL(wp)           :: delta_soil
  REAL(wp)           :: da_m_day ! slope of daytime regression deltaCa*Ca=m*Ca+b
  REAL(wp)           :: da_b_day ! intercept of daytime regression deltaCa*Ca=m*Ca+b
  REAL(wp)           :: da_m_night ! slope of nighttime regression deltaCa*Ca=m*Ca+b
  REAL(wp)           :: da_b_night ! intercept of nighttime regression deltaCa*Ca=m*Ca+b
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Litter ---
  REAL(wp)           :: z_litter ! depth of litter layer [m] ! 0 = no litter
  ! z_litter = 0 !
  ! litter density
  !   Hainich: annual litter=200 gC/m2 =400 gTG/m2 = 0.04 g/cm3
  !   for 1 cm litter height = 40 kg/m3
  REAL(wp)           :: rho_l
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Soil ---
  INTEGER(i4)        :: saxton ! (0) Clapp & Hornberger, (1) Saxton et al. for soil texture
  ! rel. root distribution calculated with Jackson et al. (Oecologia, 1996)
  REAL(wp)           :: root_factor ! Jackson''s gamma
  REAL(wp)           :: clay_factor ! clay factor
  REAL(wp)           :: sand_factor ! sand factor
  ! soil layer depth (layer boundaries) [m]
  REAL(wp), DIMENSION(0:nsoilmax) :: z_soil_in
  ! soil bulk density profile
  REAL(wp), DIMENSION(nsoilmax) :: bulk_density_in
  ! clay content in each layer
  REAL(wp), DIMENSION(nsoilmax) :: clay_in
  ! sand content in each layer
  REAL(wp), DIMENSION(nsoilmax) :: sand_in
  ! Organic matter [%]
  REAL(wp), DIMENSION(nsoilmax) :: om_in
  ! Volumetric Gravel fraction
  REAL(wp), DIMENSION(nsoilmax) :: gravel_in
  ! initial values for soil water content [fraction]
  REAL(wp), DIMENSION(nsoilmax) :: theta_in
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Water Isotopes ---
  ! For water isotope test purposes
  ! (0) normal soil water isotope fractionation ! (1) test code with no fractionation during soil evaporation
  INTEGER(i4)        :: wiso_nofracsoil !
  ! (0) normal litter water isotope fractionation ! (1) test code with no fractionation during litter evaporation
  INTEGER(i4)        :: wiso_nofraclitter !
  ! (0) normal leaf water isotope fractionation ! (1) test code with no fractionation during leaf transpiration
  INTEGER(i4)        :: wiso_nofracleaf !
  ! (0) normal rain water isotopes ! (1) test code with rain water isotopes same as initial soil water isotopes
  INTEGER(i4)        :: wiso_nofracin !
  ! Calc water isotopes iteratively in loop or just once after the
  ! normal water iteration loop: (0) Once ! (1) In loop
  INTEGER(i4)        :: wiso_implicit !
  INTEGER(i4)        :: merlivat ! (0) Cappa et al., (1) Merlivat kinetic fractionation factors
  ! delta-18O initial values
  REAL(wp), DIMENSION(nsoilmax) :: theta1_in
  ! delta-2H initial values
  REAL(wp), DIMENSION(nsoilmax) :: theta2_in
  ! delta-16O initial values
  REAL(wp), DIMENSION(nsoilmax) :: theta3_in
  ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  ! --- Special Nate McDowell Section ---
  ! Physiological parameters for upper canopy and understory
  ! On/Off: (0) 'normal' CanVeg, (1) special understory etc. treatment for Nate McDowell''s New Mexico site
  INTEGER(i4)        :: extra_nate
  ! Level at which starts upper canopy, i.e. below is understory
  INTEGER(i4)        :: nup
  ! Photosynthesis
  ! Vcmax(25C) carboxylation rate at 25°C, [umol m-2 s-1]
  REAL(wp)           :: vc25_up ! upper canopy
  REAL(wp)           :: vc25_down ! understory
  ! Jmax/Vcmax(25C) ratio of jmax to vcmax at 25 deg C, jmax = electron transport rate at 25 deg C, [umol m-2 s-1]
  REAL(wp)           :: jm_vc_up ! upper canopy
  REAL(wp)           :: jm_vc_down ! understory
  ! Rd/Vcmax ratio of rd to vcmax, rd = dark respiration at 25 C, [umol m-2 s-1]
  REAL(wp)           :: rd_vc_up ! upper canopy
  REAL(wp)           :: rd_vc_down ! understory
  ! Stomata - Leunig
  ! g0, intercept of Leuning stomata model for water vapor [mol(H2O) m-2 s-1]
  REAL(wp)           :: g0_up ! upper canopy
  REAL(wp)           :: g0_down ! understory
  ! a1, slope of Leuning stomata model for water vapor [-]
  REAL(wp)           :: a1_up ! upper canopy
  REAL(wp)           :: a1_down ! understory
  ! Do, empirical coefficient for humidity sensitivity of Leuning stomata model [Pa]
  REAL(wp)           :: D0_up ! upper canopy
  REAL(wp)           :: D0_down ! understory
  ! Stomata - Ball-Berry
  ! Ball-Berry stomatal coefficient for stomatal conductance, kball [-]
  REAL(wp)           :: kball_up ! upper canopy
  REAL(wp)           :: kball_down ! understory
  ! intercept of Ball-Berry model, (mol(H2O) m-2 s-1), bprime, intercept for H2O
  REAL(wp)           :: bprime_up ! upper canopy
  REAL(wp)           :: bprime_down ! understory

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE ini_namelist()

    IMPLICIT NONE

    INTEGER(i4) :: i

    workdir     = "./"                     ! working directory
    indir       = "../input/"              ! location of input file
    metinfile   = "2006input_eddy.dat"     ! met input file
    laiinfile   = "2006laiinput.dat"       ! lai input file
    wisoinfile  = "2006wisoinput.dat"      ! wiso input fie
    outdir      = "../output/"             ! location of output files
    outsuffix   = ".csv"                   ! output suffix
    dispfile    = "../input/DIJHainich._C" ! location and name of disperion matrix
    netcdf_in   = 0                        ! 1: input in netcdf; 0: input in ascii format
    netcdf_out  = 0                        ! 1: output in netcdf; 0: output in ascii format
    netcdf_disp = 0                        ! 1: dispersion matrix in netcdf; 0: dispersion in ascii format
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Model Set-up
    year0          = 2006 ! Year
    start_run      = 0010000_i8 ! start of model run, [dddhhmm] > 0010100
    end_run        = 0050100_i8 ! end of model run, [dddhhmm] < 3662400
    start_profiles = 0020000_i8 ! start of profile output, [dddhhmm]
    end_profiles   = 0050000_i8 ! end of profile output, [dddhhmm]
    time_step      = 3600._wp   ! time step between consecutive computation times
    ! sets soil respiration reference temperature to
    !   (1) input soil temperature or
    !   (0) modeled soil temperature at 5 cm
    switch_soil_resp_temp = 0   !
    ! set coupled stomata-photosynthesis model to
    !   (0) Ball Berry (Baldocchi analytical solution) or
    !   (1) Leuning and mesophyll conductance (Knohl analytical solution)
    switch_ball = 1             !
    ! Isoprene
    switch_isoprene = 1         ! (1) calc isoprene ! (0) no isoprene
    ! 13CO2
    switch_d13c = 1             ! (1) calc d13C ! (0) no d13C
    ! Water isotopes
    switch_wiso = 1             ! (1) calc water isotopes ! (0) no water isotopes
    ! How to determine autotrophic respiration
    !   (0) total = 50% auto + 50% hetero
    !   (1) from BETHY (Knorr 1997)
    switch_bethy_resp = 0       !
    ! Debug by not allowing condensation
    switch_no_negative_water_flux = 0 ! (1) all water fluxes >= 0 ! (0) no restriction, <0 possible
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Ecosystem Set-up --- Basic Parameter
    ! Site location = Mesita del Buey, NM, USA
    latitude = 35.85_wp   ! latitude  N
    longitude = 106.27_wp ! longitude E
    ! Eastern Standard Time
    perc_up = 1.1_wp          ! range of parameter variation for sensitivity analysis Yuan 2017.10.26
    perc_dn = 0.9_wp
    zone = 7.0_wp         ! delay from GMT
    ht = 2.66_wp          ! Canopy height [m]
    pai = 0.1_wp          ! Plant area index [m2 m-2]
    lai = 2.99_wp         ! Maximum leaf area index [m2 m-2]
    vc25_in = 45_wp          ! carboxylation rate at 25 deg C, [umol m-2 s-1], old = 66.3
    ! ratios of x to vcmax
    !   jmax = electron transport rate at 25 deg C, [umol m-2 s-1] old = 127.5
    jm_vc_in = 1.6_wp        !
    !  ratio of rd to vcmax, rd = dark respiration at 25 C, [umol m-2 s-1]
    rd_vc_in = 0.011_wp      !
    ustar_ref = 1._wp     ! reference ustar for Dij [umol m-2 s-1]
    ncl   = 40      ! canopy layers
    nsky  = 16      ! # of sky angle classes
    nl    = 9       ! # leaf angle classes
    nsoil = 10      ! # of soil layers (0: old Canoak formulation)
    nbeta = 5       ! # of levels for beta distribution of e.g. lai
    ndaysc13 = 20   ! # of days to remember for mean 13C discrimination
    nwiso = 4    ! # of (isotopic) waters >= 1
    zm = 5._wp            ! measurement height
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Ecosystem Set-up --- Photosynthesis Parameter
    ! Taken from Harley and Baldocchi (1995, PCE)
    hkin = 220000.0_wp   ! hkin, enthalpy term, old 200000
    skin = 710.0_wp      ! skin, entropy term
    ejm = 37000.0_wp     ! ejm, activation energy for electron transport, old 55000
    evc = 37000.0_wp     ! evc, activation energy for carboxylation, old 55000
    !  Enzyme constants & partial pressure of O2 and CO2
    !  Michaelis-Menten K values. From survey of literature.
    kc25 = 260.0_wp      ! kc25, ,kinetic coef for CO2 at 25 C,   microbars old = 274.6
    ko25 = 179.0_wp      ! ko25, kinetic coef for O2 at 25C,  millibars old = 419.8
    o2 = 210.0_wp        ! o2, oxygen concentration  umol mol-1
    ! tau is computed on the basis of the Specificity factor (102.33)
    ! times Kco2/Kh2o (28.38) to convert for value in solution to that based in air/
    ! The old value was 2321.1. New value for Quercus robor from Balaguer et al. (1996)
    tau25 = 2904.12_wp   ! tau25, tau coefficient
    ! Arrhenius constants
    ! Eact for Michaelis-Menten const. for KC, KO and dark respiration
    ! These values are from Harley
    ekc = 59356.0_wp    ! ekc, Activation energy for K of CO2 ! J mol-1, old = 80500.0
    eko = 35948.0_wp    ! eko, Activation energy for K of O2, old = 14500.0
    erd = 66400.0_wp    ! erd, activation energy for dark respiration, eg Q10=2
    ektau = -23408.0_wp ! ektau, old = -29000
    toptvc = 311._wp    ! toptvc, optimum temperature for maximum carboxylation
    toptjm = 311._wp    ! toptjm, optimum temperature for maximum electron transport
    ! curvature for light response function (ranges from 0 to 1)
    curvature = 0.9_wp
    ! leaf quantum yield, electrons
    qalpha = 0.22_wp
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Ecosystem Set-up --- Stomata Parameter
    g0_in = 0.001_wp      ! g0, intercept of Leuning stomata model for water vapor [mol(H2O) m-2 s-1]
    a1_in = 7._wp         ! a1, slope of Leuning stomata model for water vapor [-]
    D0_in = 20._wp        ! Do, empirical coefficient for humidity sensitivity of Leuning stomata model [Pa]
    gm_vc = 0.004_wp   ! ratio of mesophyll conductance to Vcmax [mol m-2 s-1]
    fthreshold = 0._wp ! threshold of plant available water below which stomata response to water stress
    fslope = 0._wp     ! slope of stomata response to water stress (should be: slope = 1/threshold)
    ! Ball-Berry stomatal coefficient for stomatal conductance
    kball_in = 10.0_wp
    ! intercept of Ball-Berry model, mol m-2 s-1
    bprime_in = 0.001_wp  ! intercept for H2O 0.0175
    ! Minimum stomatal resistance, s m-1.
    rsm = 145.0_wp
    brs = 60.0_wp      ! curvature coeffient for light response
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Ecosystem Set-up --- Leaf Parameter
    ep = 0.98_wp ! emissivity of leaves
    ! number of leaf sides with stomata (1) hypostomatous, (2) amphistomatous
    n_stomata_sides = 1
    ! multiplication factor for aerodynamic sheltering
    !   based on work by Grace and Wilson
    betfact = 1.5_wp
    ! leaf clumping factor
    markov = 0.8_wp
    ! Leaf dimension. geometric mean of length and width (m)
    lleaf = 0.1_wp                ! leaf length, m
    ! leaf development
    leaf_out = 0                  ! start of leaf out
    leaf_full = 0                 ! day of full leaves
    leaf_fall = 366               ! day of leaf fall
    leaf_fall_complete = 366      ! day of leaf fall end
    ! Distribution
    ! height of mid point of layer of lai_freq
    forall(i=1:nbetamax) ht_midpt_in(i) = real(i,wp)/0.5_wp
    ! fraction of total LAI per layer
    lai_freq_in(:) = zero
    lai_freq_in(1:5) = (/ 0.30_wp, 0.25_wp, 0.21_wp, 0.14_wp, 0.10_wp /)
    ! fraction of total PAI per layer
    pai_freq_in(:) = zero
    pai_freq_in(1:5) = (/ 0.15_wp, 0.20_wp, 0.21_wp, 0.24_wp, 0.20_wp /)
    ! optical properties of leaves and soil
    !   1. date <= leafout || date >= leaf_fall_complete
    !   2. date > leafout and date < fulleaf
    !   3. date >= fulleaf and date <= leaf_fall
    !   4. date > leaf_fall and date < leaf_fall_complete
    par_reflect_in(:)       = zero
    par_trans_in(:)         = zero
    par_soil_refl_dry_in(:) = zero
    nir_reflect_in(:)       = zero
    nir_trans_in(:)         = zero
    nir_soil_refl_dry_in(:) = zero
    par_reflect_in(1:4)       = (/ 0.30_wp, 0.09_wp, 0.15_wp, 0.09_wp /)
    par_trans_in(1:4)         = (/ 0.00_wp, 0.06_wp, 0.06_wp, 0.06_wp /)
    par_soil_refl_dry_in(1:4) = (/ 0.10_wp, 0.10_wp, 0.31_wp, 0.10_wp /)
    nir_reflect_in(1:4)       = (/ 0.50_wp, 0.43_wp, 0.35_wp, 0.43_wp /)
    nir_trans_in(1:4)         = (/ 0.00_wp, 0.26_wp, 0.26_wp, 0.26_wp /)
    nir_soil_refl_dry_in(1:4) = (/ 0.10_wp, 0.10_wp, 0.35_wp, 0.10_wp /)
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Ecosystem Set-up --- Miscellanous Parameter
    attfac = 2.5_wp               ! attenuation factor for wind speed inside canopy
    eabole = 45162_wp             ! activation energy for bole respiration for Q10 = 2.02
    ! Constants for leaf energy balance
    epsoil = 0.98_wp              ! epsoil, Emissivity of soil
    ! Interception reservoir
    water_film_thickness = 0.1_wp ! leaf water film thickness [mm m-2 LAI]
    tau_water = 0.5_wp            ! rain interception efficiency per m-2 LAI [0-1]
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- 13CO2 ---
    !  d13C of heterotrophic soil respiration (= longterm) [per mil]
    delta_soil = -26.6_wp
    da_m_day = -24.85_wp   ! slope of daytime regression deltaCa*Ca=m*Ca+b
    da_b_day = 6254.1_wp   ! intercept of daytime regression deltaCa*Ca=m*Ca+b
    da_m_night = -23.65_wp ! slope of nighttime regression deltaCa*Ca=m*Ca+b
    da_b_night = 5778.1_wp ! intercept of nighttime regression deltaCa*Ca=m*Ca+b
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Litter ---
    z_litter = 0.037_wp ! depth of litter layer [m] ! 0 = no litter
    ! z_litter = 0 !
    ! litter density
    !   Hainich: annual litter=200 gC/m2 =400 gTG/m2 = 0.04 g/cm3
    !   for 1 cm litter height = 40 kg/m3
    rho_l = 40._wp
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Soil ---
    saxton = 1 ! (0) Clapp & Hornberger, (1) Saxton et al. for soil texture
    ! rel. root distribution calculated with Jackson et al. (Oecologia, 1996)
    root_factor = 0.964_wp ! Jackson''s gamma
    clay_factor = 1._wp ! clay factor
    sand_factor = 1._wp ! sand factor
    ! soil layer depth [m]
    forall(i=1:nsoilmax) z_soil_in(i) = real(i,wp)
    z_soil_in(0:10) = (/ 0._wp, 0.037_wp, 0.086_wp, 0.153_wp, 0.243_wp, 0.365_wp, &
         0.529_wp, 0.751_wp, 1.050_wp, 1.454_wp, 2.000_wp /)
    ! soil bulk density profile
    bulk_density_in(:) = one
    bulk_density_in(1:10) = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.1_wp, &
         1.1_wp, 1.2_wp, 1.3_wp, 1.4_wp, 1.4_wp /)
    ! clay content in each layer
    clay_in(:) = 15._wp
    ! sand content in each layer
    sand_in(:) = 35._wp
    ! Organic matter [%]
    om_in(:) = half
    ! Volumetric Gravel fraction
    gravel_in(:) = zero
    ! initial values for soil water content [fraction]
    theta_in(:) = zero
    theta_in(1:10) = (/ 0.12_wp, 0.12_wp, 0.12_wp, 0.12_wp, 0.12_wp, &
         0.12_wp, 0.13_wp, 0.13_wp, 0.13_wp, 0.14_wp /)
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Water Isotopes ---
    ! For water isotope test purposes
    ! (0) normal soil water isotope fractionation ! (1) test code with no fractionation during soil evaporation
    wiso_nofracsoil   = 0 !
    ! (0) normal litter water isotope fractionation ! (1) test code with no fractionation during litter evaporation
    wiso_nofraclitter = 0 !
    ! (0) normal leaf water isotope fractionation ! (1) test code with no fractionation during leaf transpiration
    wiso_nofracleaf   = 0 !
    ! (0) normal rain water isotopes ! (1) test code with rain water isotopes same as initial soil water isotopes
    wiso_nofracin     = 0 !
    ! Calc water isotopes iteratively in loop or just once after the
    ! normal water iteration loop: (0) Once ! (1) In loop
    wiso_implicit = 1 !
    merlivat = 0 ! (0) Cappa et al., (1) Merlivat kinetic fractionation factors
    ! delta-18O initial values
    theta1_in(:) = -2._wp
    ! delta-2H initial values
    theta2_in(:) = -10._wp
    ! delta-16O initial values
    theta3_in(:) = zero
    ! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    ! --- Special Nate McDowell Section ---
    ! Physiological parameters for upper canopy and understory
    ! On/Off: (0) 'normal' CanVeg, (1) special understory etc. treatment for Nate McDowell''s New Mexico site
    extra_nate = 1
    ! Level at which starts upper canopy, i.e. below is understory
    nup = 9
    ! Photosynthesis
    ! Vcmax(25C) carboxylation rate at 25°C, [umol m-2 s-1]
    vc25_up   = 35._wp ! upper canopy
    vc25_down = 35._wp ! understory
    ! Jmax/Vcmax(25C) ratio of jmax to vcmax at 25 deg C, jmax = electron transport rate at 25 deg C, [umol m-2 s-1]
    jm_vc_up   = 1.6_wp ! upper canopy
    jm_vc_down = 1.6_wp ! understory
    ! Rd/Vcmax ratio of rd to vcmax, rd = dark respiration at 25 C, [umol m-2 s-1]
    rd_vc_up   = 0.020_wp ! upper canopy
    rd_vc_down = 0.020_wp ! understory
    ! Stomata - Leunig
    ! g0, intercept of Leuning stomata model for water vapor [mol(H2O) m-2 s-1]
    g0_up   = 0.001_wp ! upper canopy
    g0_down = 0.001_wp ! understory
    ! a1, slope of Leuning stomata model for water vapor [-]
    a1_up   = 4._wp ! upper canopy
    a1_down = 7._wp ! understory
    ! Do, empirical coefficient for humidity sensitivity of Leuning stomata model [Pa]
    D0_up   = 20._wp ! upper canopy
    D0_down = 20._wp ! understory
    ! Stomata - Ball-Berry
    ! Ball-Berry stomatal coefficient for stomatal conductance, kball [-]
    kball_up   = 7._wp ! upper canopy
    kball_down = 10._wp ! understory
    ! intercept of Ball-Berry model, (mol(H2O) m-2 s-1), bprime, intercept for H2O
    bprime_up   = 0.001_wp ! upper canopy
    bprime_down = 0.001_wp ! understory

  END SUBROUTINE ini_namelist


  ! ------------------------------------------------------------------
  SUBROUTINE read_namelist()

    USE nml,          ONLY: open_nml, position_nml, close_nml, nnml, POSITIONED
    USE types,        ONLY: time, iswitch, srf_res, soil, ciso, wiso, &
                            alloc_type_vars, zero_new_timestep, zero_initial
    USE utils,        ONLY: inv_boltz
#ifdef DEBUG
    USE messages,     ONLY: message
#endif

    IMPLICIT NONE

    INTEGER(i4) :: ierr
    REAL(wp), DIMENSION(:), ALLOCATABLE :: jm25 ! jmax @ 25C

    NAMELIST /canctl/  &
         workdir, indir, metinfile, laiinfile, wisoinfile, outdir, outsuffix, dispfile, &
         netcdf_in, netcdf_out, netcdf_disp, year0, start_run, end_run, start_profiles, &
         end_profiles, time_step, switch_soil_resp_temp, switch_ball, switch_isoprene, switch_d13c, &
         switch_wiso, switch_bethy_resp, switch_no_negative_water_flux, latitude, longitude, &
         perc_up, perc_dn, zone, ht, pai, lai, vc25_in, &
         jm_vc_in, rd_vc_in, ustar_ref, ncl, nsky, nl, nsoil, nbeta, ndaysc13, nwiso, zm, hkin, skin, ejm, evc, kc25, ko25, &
         o2, tau25, ekc, eko, erd, ektau, toptvc, toptjm, curvature, qalpha, g0_in, a1_in, D0_in, gm_vc, &
         fthreshold, fslope, kball_in, bprime_in, rsm, brs, ep, n_stomata_sides, betfact, markov, lleaf, leaf_out, &
         leaf_full, leaf_fall, leaf_fall_complete, ht_midpt_in, lai_freq_in, pai_freq_in, par_reflect_in, par_trans_in, &
         par_soil_refl_dry_in, nir_reflect_in, nir_trans_in, nir_soil_refl_dry_in, attfac, eabole, epsoil, &
         water_film_thickness, tau_water, delta_soil, da_m_day, da_b_day, da_m_night, da_b_night, z_litter, rho_l, &
         saxton, root_factor, clay_factor, sand_factor, z_soil_in, bulk_density_in, clay_in, sand_in, om_in, &
         gravel_in, theta_in, wiso_nofracsoil, wiso_nofraclitter, wiso_nofracleaf, wiso_nofracin, wiso_implicit, merlivat, &
         theta1_in, theta2_in, theta3_in, extra_nate, nup, vc25_up, vc25_down, jm_vc_up, jm_vc_down, rd_vc_up, rd_vc_down, &
         g0_up, g0_down, a1_up, a1_down, D0_up, D0_down, kball_up, kball_down, bprime_up, bprime_down
    call ini_namelist()
    call open_nml(namelist_file, ninnml, quiet=1)
    call position_nml('canctl', status=ierr)
    select case (ierr)
     case (POSITIONED)
       read(nnml, canctl)
    end select
    call close_nml()

    ! Setup model
    ntl  = 3*ncl
    if (switch_wiso == 0) nwiso=0 !originally 1, Yuan changed based on C4.0 line 1619. 2017.08.23
    if (nwiso < 1) then !2017.10.04
       nwiso = 1
#ifdef DEBUG
       call message('READ_NAMELIST: ','nwiso set to 1 because it includes normal water.')
#endif
    end if

    call alloc_type_vars()
    call zero_new_timestep()
    call zero_initial()

    ! Set parameters in types
    time%year0     = year0
    time%time_step = time_step
    iswitch%soil_resp_temp    = switch_soil_resp_temp
    iswitch%ball              = switch_ball
    iswitch%isoprene          = switch_isoprene
    iswitch%d13c              = switch_d13c
    iswitch%wiso              = switch_wiso
    iswitch%bethy_resp        = switch_bethy_resp
    iswitch%no_neg_water_flux = switch_no_negative_water_flux
    srf_res%fthreshold = fthreshold
    srf_res%fslope     = fslope
    ciso%delta_soil = delta_soil
    ciso%da_m_day   = da_m_day
    ciso%da_b_day   = da_b_day
    ciso%da_m_night = da_m_night
    ciso%da_b_night = da_b_night
    soil%z_litter    = z_litter
    soil%rho_l       = rho_l
    soil%saxton      = saxton
    soil%root_factor = root_factor
    soil%clay_factor = clay_factor
    soil%sand_factor = sand_factor
    wiso%nofracsoil   = wiso_nofracsoil
    wiso%nofraclitter = wiso_nofraclitter
    wiso%nofracleaf   = wiso_nofracleaf
    wiso%nofracin     = wiso_nofracin
    wiso%implicit     = wiso_implicit
    wiso%merlivat     = merlivat

    ! Allocate and set arrays
    if (.not. allocated(vc25)) allocate(vc25(ncl))
    if (.not. allocated(jm_vc)) allocate(jm_vc(ncl))
    if (.not. allocated(rd_vc)) allocate(rd_vc(ncl))
    if (.not. allocated(g0)) allocate(g0(ncl))
    if (.not. allocated(a1)) allocate(a1(ncl))
    if (.not. allocated(D0)) allocate(D0(ncl))
    if (.not. allocated(kball)) allocate(kball(ncl))
    if (.not. allocated(bprime)) allocate(bprime(ncl))
    vc25(:)   = vc25_in*perc_up ! 10% up for sensitivity analysis
 !   print *, "perc_up = ",perc_up
!    print *, "vc25 = ",vc25
    jm_vc(:)  = jm_vc_in*perc_up
    rd_vc(:)  = rd_vc_in*perc_up
    g0(:)     = g0_in*perc_up
    a1(:)     = a1_in*perc_up
    D0(:)     = D0_in*perc_up
    kball(:)  = kball_in*perc_up
    bprime(:) = bprime_in*perc_up

    ! Calc derived variables
    zd        = twothird*ht ! displacement height [m]
    z0        = 0.1_wp*ht     ! rougness lenght [m]
    izref     = ceiling(zm*real(ncl,kind=wp)/ht,kind=i4) ! array value of reference height = measurement height*ncl/ht
    delz      = ht/ncl ! height of each layer, ht/ncl
    zh65      = 0.65_wp/ht
    epsigma   = ep * sigma !2017.10.04
    epsigma2  = 2.0_wp * ep * sigma !2017.10.04
    epsigma4  = 4.0_wp * ep * sigma
    epsigma6  = 6.0_wp * ep * sigma
    epsigma8  = 8.0_wp * ep * sigma
    epsigma12 = 12.0_wp * ep * sigma
    qalpha2   = qalpha*qalpha*perc_up*perc_up

    ! Set parameters in type arrays
    soil%z_soil(0:nsoil)       = z_soil_in(0:nsoil) *perc_up!2017.10.04
!    print *, "z_soil_in ",z_soil_in
!    print *, "z_soil", soil%z_soil
    soil%bulk_density(1:nsoil) = bulk_density_in(1:nsoil)*perc_up
    soil%clay_in(1:nsoil)      = clay_in(1:nsoil)*perc_up
    soil%sand_in(1:nsoil)      = sand_in(1:nsoil)*perc_up
    soil%om(1:nsoil)           = om_in(1:nsoil)*perc_up
    soil%gravel(1:nsoil)       = gravel_in(1:nsoil)*perc_up
    soil%theta(1:nsoil,1)      = theta_in(1:nsoil)*perc_up
    wiso%dtheta(1:nsoil,1)     = theta_in(1:nsoil)*perc_up
    if (nwiso >= 2) wiso%dtheta(1:nsoil,2) = theta1_in(1:nsoil)*perc_up ! 18O
    if (nwiso >= 3) wiso%dtheta(1:nsoil,3) = theta2_in(1:nsoil)*perc_up ! 2H
    if (nwiso >= 4) wiso%dtheta(1:nsoil,4) = theta3_in(1:nsoil)*perc_up ! Normal water

    ! LAI beta distribution
    if (.not. allocated(ht_midpt)) allocate(ht_midpt(nbeta))
    if (.not. allocated(lai_freq)) allocate(lai_freq(nbeta))
    if (.not. allocated(pai_freq)) allocate(pai_freq(nbeta))
    ht_midpt(1:nbeta) = ht_midpt_in(1:nbeta)
    lai_freq(1:nbeta) = lai_freq_in(1:nbeta)
    pai_freq(1:nbeta) = pai_freq_in(1:nbeta)

    ! Leaf optical properties
    nlop = nleafopticalmax
    if (.not. allocated(par_reflect)) allocate(par_reflect(nlop))
    if (.not. allocated(par_trans)) allocate(par_trans(nlop))
    if (.not. allocated(par_soil_refl_dry)) allocate(par_soil_refl_dry(nlop))
    if (.not. allocated(nir_reflect)) allocate(nir_reflect(nlop))
    if (.not. allocated(nir_trans)) allocate(nir_trans(nlop))
    if (.not. allocated(nir_soil_refl_dry)) allocate(nir_soil_refl_dry(nlop))
    par_reflect(1:nlop)       = par_reflect_in(1:nlop)*perc_up
 !   print *, par_reflect
    par_trans(1:nlop)         = par_trans_in(1:nlop)*perc_up
 !   print *, par_trans
    par_soil_refl_dry(1:nlop) = par_soil_refl_dry_in(1:nlop)*perc_up
    nir_reflect(1:nlop)       = nir_reflect_in(1:nlop)*perc_up
    nir_trans(1:nlop)         = nir_trans_in(1:nlop)*perc_up
    nir_soil_refl_dry(1:nlop) = nir_soil_refl_dry_in(1:nlop)*perc_up

    ! Extra Nate has two layers of vegetation: 1:nup-1 and nup:ncl
    if (extra_nate == 1) then
       vc25(1:nup-1)   = vc25_down
       jm_vc(1:nup-1)  = jm_vc_down
       rd_vc(1:nup-1)  = rd_vc_down
       g0(1:nup-1)     = g0_down
       a1(1:nup-1)     = a1_down
       D0(1:nup-1)     = D0_down
       kball(1:nup-1)  = kball_down
       bprime(1:nup-1) = bprime_down
       vc25(nup:ncl)   = vc25_up
       jm_vc(nup:ncl)  = jm_vc_up
       rd_vc(nup:ncl)  = rd_vc_up
       g0(nup:ncl)     = g0_up
       a1(nup:ncl)     = a1_up
       D0(nup:ncl)     = D0_up
       kball(nup:ncl)  = kball_up
       bprime(nup:ncl) = bprime_up
    end if

    ! Variables per layer
    if (.not. allocated(jm25)) allocate(jm25(ncl))
    if (.not. allocated(vcopt)) allocate(vcopt(ncl))
    if (.not. allocated(jmopt)) allocate(jmopt(ncl))
    jm25(:)  = jm_vc(:)*vc25(:) ! jmax from vcmax @ 25C
    vcopt(:) = inv_boltz(vc25(:), evc, toptvc, TN0+25._wp, hkin)
    jmopt(:) = inv_boltz(jm25(:), ejm, toptjm, TN0+25._wp, hkin)

  END SUBROUTINE read_namelist


  ! ------------------------------------------------------------------
  FUNCTION inv_boltz(rate, eakin, topt, tl, hkin_in)
    ! calculates the inverse of Boltzmann temperature distribution for photosynthesis
    ! used to get Vcmax and Jmax at Toptimum from Vcmax and Jmax at 25 deg C
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN) :: rate
    REAL(wp),               INTENT(IN) :: eakin
    REAL(wp),               INTENT(IN) :: topt
    REAL(wp),               INTENT(IN) :: tl
    REAL(wp),               INTENT(IN) :: hkin_in
    REAL(wp), DIMENSION(1:size(rate))  :: inv_boltz

    REAL(wp) :: dtlopt, prodt, denom, numm

    dtlopt    = tl - topt
    prodt     = rugc * topt * tl
    denom     = hkin_in * exp(eakin * (dtlopt) / (prodt))
    numm      = hkin_in - eakin * (one - exp(hkin_in * (dtlopt) / (prodt)))
    inv_boltz = rate * (numm / denom)

  END FUNCTION inv_boltz

END MODULE parameters
