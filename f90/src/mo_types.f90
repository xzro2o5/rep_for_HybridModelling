MODULE types

  !
  ! Module defines types, assigns variables of these types and allocates memory
  ! to members of the assigned variables.
  !
  ! Written Jan 2011, Matthias Cuntz - Ported C-Code
  !

  USE kinds, ONLY: wp, i4, i8
  USE setup, ONLY: ncl, ntl, nsoil, nl, nwiso, ndaysc13!, nsky

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: output_variables, flux_variables, input_variables, time_variables, &          ! Types
            meteorology, surface_resistances, factors, bole_respiration_structure, &
            canopy_architecture, non_dimensional_variables, boundary_layer_resistances, &
            solar_radiation_variables, soil_variables, profile, switch_variables, &
            isotope_variables, water_isotope_variables, debug_variables, nitrogen_variables ! Yuan added 2018.05.09, store variables to be debug
  PUBLIC :: debug, output, flux, input, time, met, srf_res, fact, bole, canopy, &                ! Type variables
            non_dim, bound_lay_res, solar, soil, prof, iswitch, ciso, wiso, nitrogen
  PUBLIC :: alloc_type_vars, zero_new_timestep, zero_initial                              ! Subroutines

  ! -----------------------------------------------------------------------------------------------
  ! Type definitions
  ! debug variables: to store variables for debug use
    TYPE debug_variables
     REAL(wp) :: R1
     REAL(wp) :: R2
     REAL(wp) :: R3
     REAL(wp) :: R4
     REAL(wp) :: R5
     REAL(wp) :: R6
     REAL(wp) :: R7
     REAL(wp) :: R8
     REAL(wp) :: R9
     INTEGER(i8) :: in1
     INTEGER(i8) :: in2
     INTEGER(i8) :: in3
     INTEGER(i8) :: in4
     INTEGER(i8) :: in5
     INTEGER(i8) :: in6
     INTEGER(i8) :: in7
     INTEGER(i8) :: in8
     INTEGER(i8) :: in9
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D1
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D3
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D4
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D5
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D6
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D7
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D8
     REAL(wp), DIMENSION(:), ALLOCATABLE :: D9

  END TYPE debug_variables

  TYPE output_variables
     REAL(wp) :: c1
     REAL(wp) :: c2
     REAL(wp) :: c3
     REAL(wp) :: c4
     REAL(wp) :: c5
     REAL(wp) :: c6
     REAL(wp) :: c7
     REAL(wp) :: c8
     REAL(wp) :: c9
     REAL(wp) :: c10
     REAL(wp) :: netrad
     REAL(wp) :: sumrn
     REAL(wp) :: sumlout
     REAL(wp) :: sumh
     REAL(wp) :: sumle
     REAL(wp) :: can_ps_mol
     REAL(wp) :: can_gpp
     REAL(wp) :: c_transpiration_mole
     REAL(wp) :: canresp
     REAL(wp) :: sumksi
     REAL(wp) :: tleaf_mean
     REAL(wp) :: tavg_sun
     REAL(wp) :: tavg_shd
     REAL(wp) :: ave_disc13
     REAL(wp) :: ave_disc13_long
     REAL(wp) :: ave_cica
     REAL(wp) :: ave_gs
     REAL(wp) :: ave_daC13
     REAL(wp) :: ave_disc13_a
     REAL(wp) :: ave_disc13_ab
     REAL(wp) :: ave_disc13_asal
     REAL(wp) :: ave_disc13_b
     REAL(wp) :: ave_csca
     REAL(wp) :: ave_ccca
     REAL(wp) :: ave_gm
     REAL(wp) :: isoprene_efflux
     REAL(wp) :: diff_par
     REAL(wp) :: rnet_soil
     REAL(wp) :: sumfc
     REAL(wp) :: sumevap
     REAL(wp) :: sumsens
     REAL(wp) :: sumps
     REAL(wp) :: sumgpp
     REAL(wp) :: sumpar
     REAL(wp) :: sumnet
     REAL(wp) :: sumbole
     REAL(wp) :: sumsoil
     REAL(wp) :: sumta
     REAL(wp) :: sumgs
     REAL(wp) :: sumresp
     REAL(wp) :: sumTleaf
     REAL(wp) :: sumF13C
     REAL(wp) :: sumdisc13C
     REAL(wp) :: sumdisc13C_long_day
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dLAIdz ! ! [ncl]
     ! oxygen output Yuan 2018.01.30
     REAL(wp) :: sumo ! sum*** means daily output
     REAL(wp) :: sumneto
     REAL(wp) :: sumROC
     REAL(wp) :: sumresp_o
     REAL(wp) :: hour_canrespo ! O2 via leaf respiration
     REAL(wp) :: houro
     REAL(wp) :: hourneto
     REAL(wp) :: hourROC ! Yuan added hourly ROC 2018.05.07
  END TYPE output_variables


  TYPE flux_variables
     REAL(wp) :: photosyn ! photosynthesis, um m-2 s-1
     REAL(wp) :: nee ! net ecosystem exchange, umol m-2 s-1
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c_evaporation ! canopy evaporation from interception reservoir, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c_transpiration ! canopy transpiration through stomata, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c_evapotranspiration ! canopy transpiration + evaporation, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: soilevap ! soil evaporation, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: litterevap ! litter evaporation, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: s_evap ! soil+litter evaporation, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: evapotranspiration ! canopy+soil+litter evaporationtranspiration, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: soilinfl ! soil infiltration, mm m-2 s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: surfrun ! surface runoff = water flux lateral to surface, mm m-2 s-1 ! [nwiso]
     REAL(wp) :: albedo                             ! surface albedo = (PARout + NIRout)/(PARin + NIRin)
  END TYPE flux_variables


  TYPE input_variables
     INTEGER(i4) :: dayy ! day
     REAL(wp) :: hhrr ! hour
     REAL(wp) :: ta ! air temperature, C
     REAL(wp) :: rglobal ! global radiation, W m-2
     REAL(wp) :: parin ! photosynthetically active radiation, micromole m-2 s-1
     REAL(wp) :: pardif ! diffuse PAR, micromol m-2 s-1
     REAL(wp) :: lai ! for Nate McDowell''s juniper site read lai instead of pardif
     REAL(wp) :: lai_up ! for Nate McDowell''s juniper site read lai in upper canopy
     REAL(wp) :: lai_down ! for Nate McDowell''s juniper site read lai in understory
     REAL(wp) :: ea ! vapor pressure, kPa
     REAL(wp) :: wnd ! wind speed, m s-1
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ppt ! precipitation, mm per hour ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dppt ! delta value of precipitation [per mi] ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dvapour ! delta value of vapour [per mi] ! [nwiso]
     REAL(wp) :: co2air       ! CO2 concentration, ppm
     REAL(wp) :: press_mb     ! air pressure, mb
     REAL(wp) :: tsoil        ! soil temperature in 50 cm, degree C ! Yuan's correction 2018.07.30. soil T at 5cm depth!!! refer to paper Astrid R. B. S�e  Nina Buchmann
     REAL(wp) :: soilmoisture ! soil moisture in 15 cm, %
     INTEGER(i4) :: flag      ! input coding
     REAL(wp) :: longwave     ! long wave irradiannce, Wm-2
     REAL(wp) :: d13CO2       ! d13C of atmospheric CO2 [permille]
     REAL(wp) :: d18CO2       ! d18O of atmospheric CO2 [permille]
     REAL(wp) :: o2air        ! O2 concentration, ppm Yuan 2018.02.14
     REAL(wp) :: ER           ! net Ass O2:CO2 chamber measurements Yuan 2022.07.14
  END TYPE input_variables


  ! structure for time variables
  TYPE time_variables
     REAL(wp) :: local_time
     INTEGER(i8) :: daytime ! day+hour coding, e.g. 1230830
     INTEGER(i4) :: doy  ! day of year
     INTEGER(i4) :: days ! day
     INTEGER(i4) :: jdold ! previous day
     INTEGER(i4) :: count ! number of iterations

     INTEGER(i4) :: year0 ! always first year of data
     INTEGER(i4) :: year ! actual year
     INTEGER(i4) :: leafout ! day of leaf out
     INTEGER(i4) :: leaffull ! date of full leaf
     INTEGER(i4) :: leaffall ! date of leaf shedding
     INTEGER(i4) :: leaffallcomplete ! date of leaf shedding completed
     REAL(wp) :: lai ! lai as a function of time
     REAL(wp) :: pai ! pai as wai+lai Yuan 2018.03.02
     REAL(wp) :: wai ! daily wai as a constant 1.1
     REAL(wp) :: time_step ! time between subsequent computation times
  END TYPE time_variables


  ! structure for meteorological variables
  TYPE meteorology
     REAL(wp) :: ustar ! friction velocity, m s-1
     REAL(wp) :: ustar_filter ! updated friction velocity with new H, m s-1
     REAL(wp) :: K ! momentum diffusivity
     REAL(wp) :: K_filter ! momentum diffusivity
     REAL(wp) :: phim ! stability correction of ustar
     REAL(wp) :: rhova_g ! absolute humidity, g m-3
     REAL(wp) :: rhova_kg ! absolute humidity, kg m-3
     REAL(wp) :: H ! sensible heat flux, W M-2
     REAL(wp) :: H_filter ! old sensible heat flux, W m-2
     REAL(wp) :: air_density ! air density, kg m-3
     REAL(wp) :: T_Kelvin ! absolute air temperature, K
     REAL(wp) :: zl ! z over L, height normalized by the Obukhov length
     REAL(wp) :: press_kpa ! station pressure, kPa
     REAL(wp) :: press_bars ! station pressure, bars
     REAL(wp) :: press_Pa ! pressure, Pa
     REAL(wp) :: pstat273 ! gas constant computations
     REAL(wp) :: air_density_mole ! air density, mole m-3
     REAL(wp) :: relative_humidity ! relative humidity

     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dispersion ! Lagrangian dispersion matrix, s m-1 ! [ntl,ncl]
  END TYPE meteorology


  ! structure for surface resistances and conductances
  TYPE surface_resistances
     REAL(wp) :: fdrought ! drought response factor of stomata conducatance [-]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rcuticle ! cuticle resistance, m2 s mol-1 ! [ncl]
     REAL(wp) :: fthreshold ! threshold of relative soil water below which stomata show response to water stress [mm]
     REAL(wp) :: fslope ! slope of water stress response
  END TYPE surface_resistances


  ! structure for plant and physical factors
  TYPE factors
     REAL(wp) :: latent ! latent heat of vaporization, J kg-1
     REAL(wp) :: heatcoef ! factor for sensible heat flux density
     REAL(wp) :: a_filt ! filter coefficients
     REAL(wp) :: co2 ! CO2 factor, ma/mc * rhoa (mole m-3)
     REAL(wp) :: o2  ! O2 factor, ma/mc * rhoa (mole m-3)
  END TYPE factors


  ! structure for bole respiration and structure
  TYPE bole_respiration_structure
     REAL(wp) :: factor ! base respiration, micromoles m-2 s-1, data from Edwards
     REAL(wp) :: respiration_mole ! bole respiration, micromol m-2 s-1
     REAL(wp) :: respiration_mg ! bole respiration, mg CO2 m-2 s-1
     REAL(wp) :: calc ! calculation factor
     REAL(wp), DIMENSION(:), ALLOCATABLE :: layer ! bole pai per layer ! [ncl]
  END TYPE bole_respiration_structure


  ! structure for canopy architecture
  TYPE canopy_architecture
     REAL(wp), DIMENSION(:), ALLOCATABLE :: bdens ! probability density of leaf angle ! [10]
  END TYPE canopy_architecture


  ! structure for non dimensional variables
  TYPE non_dimensional_variables
     ! Prandtl Number
     REAL(wp) :: pr
     REAL(wp) :: pr33
     ! Schmidt number for vapor
     REAL(wp) :: sc
     REAL(wp) :: sc33
     ! Schmidt number for CO2
     REAL(wp) :: scc
     REAL(wp) :: scc33
     ! Schmidt number for ozone
     REAL(wp) :: sco3
     REAL(wp) :: sco333
     ! Grasshof number
     REAL(wp) :: grasshof
     ! multiplication factors with leaf length and diffusivity
     REAL(wp) :: lfddv
     REAL(wp) :: lfddh
  END TYPE non_dimensional_variables


  ! boundary layer resistances
  TYPE boundary_layer_resistances
     REAL(wp) :: vapor ! resistance for water vapor, s/m
     REAL(wp) :: heat ! resistance for heat, s/m
     REAL(wp) :: co2 ! resistance for CO2, s/m
  END TYPE boundary_layer_resistances


  ! radiation variables, visible, near infrared and infrared
  TYPE solar_radiation_variables
     ! profiles of the probabilities of sun and shade leaves
     REAL(wp), DIMENSION(:), ALLOCATABLE :: prob_beam ! probability of beam or sunlit fraction ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: prob_shd ! probability of shade ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ir_dn ! downward directed infrared radiation, W m-2 ! [ncl+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ir_up ! upward directed infrared radiation. W m-2 ! [ncl+1]
     ! inputs of near infrared radiation and components, W m-2
     REAL(wp) :: nir_beam ! incoming beam component near infrared radiation, W m-2
     REAL(wp) :: nir_diffuse ! incoming diffuse component near infrared radiation, W m-2
     REAL(wp) :: nir_total ! incoming total near infrared radiaion, W m-2
     ! computed profiles of near infrared radiation, W m-2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_dn ! downward scattered near infrared radiation ! [ncl+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_up ! upward scattered near infrared radiation ! [ncl+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_sun ! near infrared radiation on sunlit fraction of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: nir_shd ! near infrared radiation on shaded fraction of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: beam_flux_nir ! flux density of direct near infrared radiation ! [ncl+1]
     ! inputs of visible light, PAR, W m-2
     REAL(wp) :: par_diffuse ! diffuse component of incoming PAR, parin
     REAL(wp) :: par_beam ! beam component of incoming PAR, parin
     ! computed profiles of visible radiation, PAR, W m-2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: par_shd ! PAR on shaded fraction of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: par_sun ! PAR on sunlit fraction of layer, beam and diffuse ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: beam_flux_par ! PAR in the beam component ! [ncl+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: par_down ! upward scattered PAR ! [ncl+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: par_up ! downward scattered PAR ! [ncl+1]
     ! flux densities of visible quanta on sun and shade leaves for photosynthesis
     ! calculations, micromoles m-2 s-1
     REAL(wp), DIMENSION(:), ALLOCATABLE :: quantum_sun ! par on sunlit leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: quantum_shd ! par on shaded leaves ! [ncl]
     ! Net radiation profiles, W m-2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rnet_sun ! net radiation flux density on sunlit fraction of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rnet_shd ! net radiation flux density on shade fraction of layer ! [ncl]
     REAL(wp) :: beta_rad ! solar elevation angle, radians
     REAL(wp) :: sine_beta ! sine of beta
     REAL(wp) :: beta_deg ! solar elevation angle, degrees
     REAL(wp) :: ratrad ! radiation ratio to detect cloud amount
     ! leaf and soil optical properities of near infrared radiation
     REAL(wp) :: nir_reflect ! leaf reflectance in the near infrared
     REAL(wp) :: nir_trans ! leaf transmittance in the near infrared
     REAL(wp) :: nir_soil_refl_dry ! soil reflectance in the near infrared of dry soil
     REAL(wp) :: nir_soil_refl ! soil reflectance in the near infrared
     REAL(wp) :: nir_absorbed ! leaf absorptance in the near infrared
     ! optical properties of leaves and soil for PAR
     REAL(wp) :: par_absorbed ! PAR leaf absorptance
     REAL(wp) :: par_reflect ! PAR leaf reflectance
     REAL(wp) :: par_trans ! PAR leaf transmittance
     REAL(wp) :: par_soil_refl_dry ! PAR soil reflectance of dry soil
     REAL(wp) :: par_soil_refl ! PAR soil reflectance
     ! Net radiation profiles, W m-2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: exxpdir ! exponential transmittance of diffuse radiation through a layer ! [ncl]
     REAL(wp) :: ratradnoon ! radiation ratio at noon for guestimating cloud amount at night
  END TYPE solar_radiation_variables


  ! physical properties of soil abd soil energy balance variables
  TYPE soil_variables
     ! soil properties
     REAL(wp), DIMENSION(:), ALLOCATABLE :: T_soil ! soil temperature [deg C] ! [0:nsoil+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: T_soil_filter ! soil temperature [deg C] ! [0:nsoil+1]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: k_conductivity_soil ! thermal conductivity of soil *dz ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: cp_soil ! specific heat of soil, f(texture, moisture) ! [nsoil]
     REAL(wp) :: T_Kelvin ! soil surface temperature in Kelvin
     REAL(wp) :: T_air ! air temperature above soil, C
     REAL(wp) :: tsrf ! soil surface temperature in C
     REAL(wp) :: tsrf_filter ! filtered soil surface temperature in C
     REAL(wp) :: rs ! soil resistances
     REAL(wp) :: rb ! soil boundary layer resistances
     REAL(wp) :: T_15cm ! soil temperature at 15 cm
     ! soil energy flux densities, W m-2
     REAL(wp) :: lout ! longwave efflux from soil
     REAL(wp) :: soilevap ! soil evaporation
     REAL(wp) :: soilevap_filter ! filtered soil evaporation
     REAL(wp) :: litterevap ! litter evaporation
     REAL(wp) :: litterevap_filter ! filtered litter evaporation
     REAL(wp) :: litterevap_save ! save litter evaporation for wiso
     REAL(wp) :: maxlitterevap ! maximum possible litter evaporation
     REAL(wp) :: evap ! total forest floor evaporation (soil+litter)
     REAL(wp) :: heat ! soil sensible heat flux density
     REAL(wp) :: rnet ! net radiation budget of the soil
     REAL(wp) :: gsoil ! soil heat flux density
     ! soil CO2 respiratory efflux
     REAL(wp) :: respiration_mole ! soil respiration, micromol m-2 s-1
     REAL(wp) :: respiration_mg ! soil respiration, mg CO2 m-2 s-1
     REAL(wp) :: base_respiration ! base rate of soil respiration, micromol m-2 s-1
     REAL(wp) :: respiration_auto ! autotrohic soil respiration, micromol m-2 s-1
     REAL(wp) :: respiration_hetero ! heterotrophic soil respiration, micromol m-2 s-1
     REAL(wp) :: resp_13 ! respiration of 13C micromole m-2 s-1
     REAL(wp) :: SR_ref_temp ! refernce temperature for soil respiration, degC
     ! from ISOLSM
     REAL(wp), DIMENSION(:), ALLOCATABLE :: swp ! soil water potential (kPa) ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: swp_mm ! soil water potential (mm) ! [nsoil]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: r ! ! [nsoil,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: a ! ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: b ! ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c ! ! [nsoil]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dwat ! ! [nsoil,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: theta ! volumetric soil water content ! [nsoil,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: qdrai ! sub-surface runoff (mm h2o s-1) ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: qinfl ! infiltration rate, mm h2o s-1 ! [nwiso]
     REAL(wp) :: soil_mm ! water content in the soil, mm m-2
     REAL(wp) :: soil_mm_root ! total water content in the soil weighted by roots, mm
     REAL(wp) :: soil_mm_50 ! water content in the top 50 cm of soil, mm m-2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: qtran ! plant transpiration, mm s-1 ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: qseva ! plant transpiration, mm s-1 ! [nwiso]
     ! Saxton & Rawls (2006)
     REAL(wp), DIMENSION(:), ALLOCATABLE :: psi ! soil water potential [Pa] ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: k_theta ! moisture conductivity [mm s-1] ! [nsoil]
     ! Litter (Ogee & Brunet 2002)
     REAL(wp), DIMENSION(:), ALLOCATABLE :: theta_l ! volumetric litter water content ! [nwiso]
     REAL(wp) :: T_l ! litter temperature [K]
     REAL(wp) :: T_l_filter ! litter temperature [K]
     REAL(wp) :: rl ! litter resistances
     REAL(wp) :: c_litterevap ! factor regulating that litter evaporation < litter moisture
     REAL(wp) :: c_litterevap_save ! save factor regulating that litter evaporation < litter moisture
     INTEGER(i4) :: maxlitter ! switch if litterevap was restricted to maximum litter evaporation
     REAL(wp) :: latent ! latent heat coefficient before soil-air vapour pressure deficit
     REAL(wp) :: lecoef ! latent heat coefficient before soil-air vapour pressure deficit
     REAL(wp) :: rh_soil ! relative humidity in first soil layer air space
     REAL(wp) :: rh_litter ! relative humidity in litter layer air space
     REAL(wp), DIMENSION(:), ALLOCATABLE :: xylem ! xylem water isotopic composition ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rxylem ! xylem water isotope ratio ! [nwiso]
     ! Fixed variables
     ! soil properties
     REAL(wp), DIMENSION(:), ALLOCATABLE :: z_soil ! depth of soil layer boundaries ! [0:nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: bulk_density ! soil bulk density ! [nsoil]
     REAL(wp) :: d ! dispacement height for soil [m]
     REAL(wp) :: z0 ! soil roughness length [m]
     REAL(wp) :: T_base ! base soil temperature
     ! For numerical soil temperature
     REAL(wp) :: temperature_dt ! time step, s
     INTEGER(i4) :: temperature_mtime ! number of time steps per dt
     ! from ISOLSM
     REAL(wp) :: moisture_dt ! time step, s
     INTEGER(i4) :: moisture_mtime ! number of time steps per dt
     REAL(wp) :: qinfl_max ! maximal infiltration rate into soil, mm m-2 s-1
     REAL(wp), DIMENSION(:), ALLOCATABLE :: clay_in ! clay fraction in each layer input ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sand_in ! sand fraction in each layer input ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: root ! actual (optimised) relative root abundance (0 to 1) ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: clay ! actual (optimised) clay fraction in each layer ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sand ! actual (optimised) sand fraction in each layer ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: om ! soil organic matter content in each layer ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: gravel ! volumetric gravel content in each layer ! [nsoil]
     REAL(wp) :: root_factor ! gamma of Jackson et al. (Oecologia 1996)
     REAL(wp) :: sand_factor ! gamma of Jackson et al. (Oecologia 1996)
     REAL(wp) :: clay_factor ! gamma of Jackson et al. (Oecologia 1996)
     REAL(wp) :: soil_mm_33_root ! total water content in the soil at field capacity (-33 kPa) weighted by roots, mm
     REAL(wp) :: soil_mm_1500_root ! total water content in the soil at wilting point (-1500 kPa) weighted by roots, mm
     ! Saxton & Rawls (2006)
     INTEGER(i4) :: saxton ! 0: Clapp & Hornberger (1978) 1: Saxton & Rawls (2006)
     REAL(wp), DIMENSION(:), ALLOCATABLE :: theta_1500 ! volumetric soil water content at wilting point ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: theta_33 ! volumetric soil water content at field capacity ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: theta_s33 ! saturation - field capacity volumetric soil water content ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: theta_s ! saturated volumetric soil water content ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: psi_e ! soil water potential at air entry (bubbling pressure) [Pa] ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rho ! normal soil density (without gravel) [kg m-3] ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: k_s ! saturated conductivity [mm s-1] ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: lambda ! coefficient for moisture conductivity ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: big_a ! coefficient for moisture conductivity ! [nsoil]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: big_b ! coefficient for moisture conductivity ! [nsoil]
     ! Litter (Ogee & Brunet 2002)
     REAL(wp) :: z_litter ! litter layer depth [m]
     REAL(wp) :: theta_ls ! saturated volumetric litter water content
     REAL(wp) :: theta_l33 ! volumetric litter water content at field capacity
     REAL(wp) :: n_l ! exponent in litter resistance calculation
     REAL(wp) :: rho_l ! litter bulk density [kg m-3]
     ! Original Dennis code
     INTEGER(i4) :: camillo ! 0: Passerat (1986) 1: Camillo & Guerney, i.e. no soil water
     ! Water
     REAL(wp), DIMENSION(:), ALLOCATABLE :: watmin ! residual vol. soil water content depending on texture ! [nsoil]
     ! Drainage
     INTEGER(i4) :: drain0 ! Switch to write out drainage <0 exactly once
     ! Lose water
     INTEGER(i4) :: lost0 ! Switch to write out when losing water starts
     ! oxygen module Yuan 2018.01.18
     REAL(wp) :: ROC_soil_air ! ROC during soil-air O2 exchange
     REAL(wp) :: o2 ! soil O2 sink
  END TYPE soil_variables


  ! Structure for Profile information,fluxes and concentrations
  TYPE profile
     ! microclimate profiles
     REAL(wp), DIMENSION(:), ALLOCATABLE :: tair ! air temp (C) ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: tair_filter ! numerical filter of Tair ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: tair_filter_save ! save for later calculations ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: u ! wind speed (m/s) ! [ntl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rhov_air ! water vapor density for all water isotopes [kg m-3] ! [ntl,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rhov_air_save ! save water vapor density [kg m-3] ! [ntl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rhov_air_filter ! numerical filter of rhov_air ! [ntl,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rhov_air_filter_save ! save numerical filter of rhov_air ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: co2_air ! co2 concentration (ppm) ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: co2_air_filter ! co2 concentration (ppm) ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: vpd_air ! vapor pressure deficit[kPa] ! [ntl]
     ! canopy structure profiles
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Gfunc_solar ! leaf-sun direction cosine function ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: tleaf ! leaf temperature per layer (sun and shade) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: isopreneflux ! isoprene flux per layer (sun and shade) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: vcmax ! vcmax in per layer (temperature corrected) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: jmax ! jmax per layer (temperature corrected) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: tp
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rd ! rd per layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: vcmaxz ! vcmaxz per layer ! [ncl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: throughfall ! throughfall in layer [mm], one extra layer for input as precipitation ! [ncl+1,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: cws ! canopy water storage in layer [mm] ! [ncl,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: wet_coef ! coefficient to regulate evaporation from wet leaf surface [0-1] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: wet_coef_filter ! filter coefficient to regulate evaporation from wet leaf surface [0-1] ! [ncl]
     ! variables for 13C isotopes
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c13cnc ! concentration of 13C ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: c12cnc ! concentration of 12C ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sour13co2 ! source/sink strength 13C ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: d13C ! del 13C ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C ! photosynthetic weighted discrimination, 13D ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C_long ! photosynthetic weighted discrimination, 13D ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C_ab    ! boundary layer discrimination ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C_a     ! stomata discrimination ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C_asal    ! mesophyll discrimination ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: disc13C_b        ! rubisco discrimination ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: cs      ! CO2 concentration at leaf surface ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ci      ! CO2 concentration inside stomatal cavity ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: cc      ! CO2 concentration at site of carboxlylation ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: csca ! ratio of cs and ca ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: cica ! ratio of ci and ca ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ccca ! ratio of cc and ca ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: gm      ! mesophyll conductance ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: d13Cair ! del 13C of the air ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: d13Cplant ! del 13C of the plant ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: R13_12_air ! 13C/12C ratio of the air ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Rplant_sun ! ratio of discriminated 13C in sunlit leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Rplant_shd ! ratio of discriminated 13C in shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Rresp ! discriminated 13C for later respiration, Ps weight ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Rresp_sum ! summing of Ps weighted values for daily ave ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Rresp_ave ! previous days discriminated 13C for later respiration ! [ncl]
     INTEGER(i4), DIMENSION(:), ALLOCATABLE :: cnt_Rresp ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: recycle ! fraction of recycled CO2, after Yakir and Sternberg ! [ntl]
     ! source/sink strengths
     REAL(wp), DIMENSION(:), ALLOCATABLE :: source_co2 ! source/sink strength of CO2 ! [ncl]
     ! fluxes for total layer, with sun/shade fractions
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dPsdz ! layer photosynthesis (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dPsdz_sun ! layer photosynthesis of sunlit area (micromol m-2 s-1) = Aphoto sun ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dPsdz_shd ! layer photosynthesis of shaded area (micromol m-2 s-1) = Aphoto shd ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGPPdz       ! layer gross primary productivity (Ps + Resp) (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGPPdz_sun   ! layer gross primary productivity (Ps + Resp) of sunlit area (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGPPdz_shd ! layer gross primary productivity (Ps + Resp) of shaded area (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGOPdz       ! layer gross primary productivity (Ps + Resp) (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGOPdz_sun   ! layer gross primary productivity (Ps + Resp) of sunlit area (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dGOPdz_shd
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dHdz ! layer sensible heat flux (W m-2) ! [ncl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dLEdz ! layer latent heat flux ! [ncl,nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dLEdz_sun ! layer latent heat flux of sunlit area (W m-2) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dLEdz_shd ! layer latent heat flux of shaded area (W m-2) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dRNdz ! layer net radiation flux (W m-2) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dLoutdz ! layer radiation flux ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dboledz ! layer bole respiration ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dRESPdz ! layer respiration (micromol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dRESPdz_sun ! layer respiration of sunlit area (micromol m-2 s-1) = rd sun ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dRESPdz_shd ! layer respiration of shaded area (micromol m-2 s-1) = rd shd ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dStomCondz ! layer stomatal conductance (mol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dStomCondz_sun ! layer stomatal conductance (mol m-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dStomCondz_shd ! layer stomatal conductance (mol m-2 s-1) ! [ncl]
     ! sun leaf variables
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_frac ! sun leaf fraction ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_tleaf ! leaf temp (C) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_GPP ! layer GPP flux for sun only (micromol mn-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_GOP
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_A ! layer A flux for sun only (micromol mn-2 s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_gs ! stomatal conductance (m s-1) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_gs_mol ! stomatal conductance of sun leaves mol m-2 s-1 ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rs ! stomatal resistance to H2O (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rs_filter ! filter stomatal resistance to H2O (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rs_save ! stomatal resistance to H2O (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rbh ! boundary layer resistance to heat (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rbv ! boundary layer resistance to H2O (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rbco2 ! boundary layer resistance to CO2 (s/m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_cs ! Cs on sun leaves (CO2 mixing ratio on leaf surface ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_ci ! Ci on sun leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_cc ! Cc (CO2 mixing ratio at site of carboxylation) on sun leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13 ! discrimination 13C ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13_ab ! del 13C boundary resistence on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13_a ! del 13C stomatal resistance on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13_asal ! del 13C mesophyll resistance on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13_b ! del 13C rubisco on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_disc13_long ! discrimination 13C ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_csca ! Cs/Ca on sunlit leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_cica ! Ci/Ca on sunlit leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_ccca ! Cc/Ca on sunlit leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_lai ! sunlit lai of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_tleaf_filter ! filtered sunlit temperature ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_wj ! electron transport rate of Ps for sun leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_wc ! carboxylatio velocity for sun leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_wp !
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_alphag ! N assimilation in photorespiration to glycine
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_alphas ! N assimilation in photorespiration to serine
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_tpu_coeff ! tpu occurs = 1 else = 0
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_resp ! respiration ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_isopreneflux ! isoprene flux per layer for sunleaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: iso_sun ! isoprene flux per leaf area in the sun ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_vpd ! vapor pressure deficit at leaf surface [hPa] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_rh ! relative humidity at leaf surface [hPa] ! [ncl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_LEstoma ! LE flux through stomata [W m-2 leaf] ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_LEwet ! LE flux from wet leaves [W m-2 leaf] ! [ncl,nwiso]
     ! shade leaf variables
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_frac ! shade leaf fraction ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_tleaf ! temperature of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_GPP ! photosynthesis (GPP) of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_GOP
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_A ! photosynthesis of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_gs ! stomatal conductance of shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_gs_mol ! stomatal conductance of shade leaves mol m-2 s-1 ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rs ! stomatal resistance of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rs_filter ! filter stomatal resistance of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rs_save ! stomatal resistance of shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rbh ! boundary layer resistance for heat on shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rbv ! boundary layer resistance for vapor on shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rbco2 ! boundary layer resistance for CO2 on shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_cs ! Cs on shade leaves (CO2 mixing ratio on leaf surface ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_ci ! Ci on shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_cc ! Cc (CO2 mixing ratio at site of carboxylation) on sun leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13 ! del 13C on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13_ab ! del 13C boundary resistence on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13_a ! del 13C stomatal resistance on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13_asal ! del 13C mesophyll resistance on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13_b ! del 13C rubisco on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_disc13_long ! discrimination 13C ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_csca ! Cs/Ca ratio on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_cica ! Ci/Ca ratio on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_ccca ! Cc/Ca ratio on shaded leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_lai ! shaded lai of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_tleaf_filter ! previous temperature ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_wc ! carboxylation rate for shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_wj ! electron transport rate for shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_wp
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_alphag
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_alphas
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_tpu_coeff ! tpu occurs = 1 else = 0
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_resp ! respiration ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_isopreneflux ! isoprene flux per layer for shade leaves ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: iso_shd ! isoprene flux per leaf area in the shade ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_vpd ! vapor pressure deficit at leaf surface [hPa] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_rh ! relative humidity at leaf surface [hPa] ! [ncl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_LEstoma ! shd leaf LE flux through stomata [W m-2 leaf] ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_LEwet ! LE flux from wet shaded leaves [W m-2 leaf] ! [ncl,nwiso]
     ! water isotope variables
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_wi !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_wi !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: wa !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_h !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_h !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: rs_fact !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_gross !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_gross !  ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_LEstoma_new ! sun leaf LE flux through stomata in accordance with water isotope [mol(H2O)/m2(leaf)s] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_LEstoma_new ! shd leaf LE flux through stomata in accordance with water isotope [mol(H2O)/m2(leaf)s] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_LEstoma_save ! save LE flux through stomata [W m-2 leaf] ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_LEstoma_save ! save LE flux from wet leaves [W m-2 leaf] ! [ncl]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_alpha_k !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_alpha_k !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_alpha_equ !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_alpha_equ !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: rvapour !  ! [ntl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dvapour !  ! [ntl,nwiso] ! Yuan 2015.05.28
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_peclet !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_peclet !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_fem !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_fem !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_craig !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_craig !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_leafwater_e !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_leafwater_e !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_leafwater_e_old !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_leafwater_e_old !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_leafwater !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_leafwater !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_trans_rtrans !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_trans_rtrans !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sun_rtrans !  ! [ncl,nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: shd_rtrans !  ! [ncl,nwiso]
     ! canopy structure profiles
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ht ! layer height (m) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dLAIdz ! leaf area index of layer (m2/m2) ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dPAIdz ! plant area index of layer ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: dWAIdz ! plant wood index of layer ! [ncl] Yuan 2018.03.01

     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: Gfunc_sky ! leaf-sky sector direction cosine function ! [ncl,nl]
     ! oxygen profiles Yuan 2018.01.17
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     sun_psn_O2! O2 in net psn, per leaf level
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     shd_psn_O2!
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dPsdz_O2_sun! per ground level
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dPsdz_O2_shd!
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dPsdz_O2!
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     sun_resp_O2! O2 in net dark resp, per leaf level
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     shd_resp_O2!
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dRESPdz_O2_sun! per ground level
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dRESPdz_O2_shd!
     REAL(wp), DIMENSION(:), ALLOCATABLE ::     dRESPdz_O2!

     REAL(wp), DIMENSION(:), ALLOCATABLE :: O2_air ! O2 concentration (ppm) ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: O2_soil ! dispersion of soil O2 flux ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: CO2_soil ! dispersion of soil CO2 flux ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: O2_disp ! dispersion of canopy O2 flux ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: CO2_disp ! dispersion of canopy CO2 flux ! [ntl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: O2_air_filter ! filtered O2 concentration (ppm) ! [ntl] Yuan 2018.07.02
     REAL(wp), DIMENSION(:), ALLOCATABLE :: source_O2 ! source/sink strength of O2 ! [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: gpp_O2 ! gross O2 flux
!     REAL(wp), DIMENSION(:), ALLOCATABLE :: rd_O2 ! dark leaf respiration [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ROC_layer ! ROC  of the layer Yuan 2018.02.26
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ROC_leaf_air !ROC during leaf-air O2 exchange [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ROC_bole_air !ROC during bole-air O2 exchange [ncl]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: RQ_sun    ! respiratory quotient of sunlit leaves
     REAL(wp), DIMENSION(:), ALLOCATABLE :: RQ_shd    ! respiratory quotient of shaded leaves
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ROC_sun    ! respiratory quotient of sunlit leaves
     REAL(wp), DIMENSION(:), ALLOCATABLE :: ROC_shd    ! respiratory quotient of shaded leaves
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_NBusch
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_NBusch
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_NO3   ! assimilated nitrate, nitrite and ammonia
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_NO3
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_NO2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_NO2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_NH4
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_NH4
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Ja_sun    ! electron for CO2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Jglu_sun    ! electron for N->glutamate
     REAL(wp), DIMENSION(:), ALLOCATABLE :: JBusch_sun    ! electron for N -> gly and serine
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Ja_shd
     REAL(wp), DIMENSION(:), ALLOCATABLE :: Jglu_shd
     REAL(wp), DIMENSION(:), ALLOCATABLE :: JBusch_shd
     REAL(wp), DIMENSION(:), ALLOCATABLE :: jphoton_sun    ! electron for CO2
     REAL(wp), DIMENSION(:), ALLOCATABLE :: jphoton_shd    ! electron for N
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_quad
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_quad
     REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_ABusch
     REAL(wp), DIMENSION(:), ALLOCATABLE :: shd_ABusch
  END TYPE profile


  TYPE switch_variables
     INTEGER(i4) :: soil_resp_temp ! 1: temp for hetero resp = input soil 0: modelled soil at 5cm
     INTEGER(i4) :: ball ! 0: Ball-Berry 1: Leuning
     INTEGER(i4) :: isoprene ! 1: calc isoprene
     INTEGER(i4) :: d13c ! 1: calc 13CO2
     INTEGER(i4) :: wiso ! 1: calc water isotopes
     INTEGER(i4) :: bethy_resp ! 0: auto/hetero=50/50 1: Bethy
     INTEGER(i4) :: no_neg_water_flux ! 1: restrict water fluxes >= 0
     INTEGER(i4) :: oxygen ! 1: calc oxygen flux
     INTEGER(i4) :: wai_new ! to switch between default wai distribution and Yuan's field measurement. 2018.07.16
     INTEGER(i4) :: tpu    ! add tpu limits
     INTEGER(i4) :: n_limit ! consider Nitrogen limit
     INTEGER(i4) :: n_random ! random mixed N source
     INTEGER(i4) :: ER ! if use chamber ER measurements

  END TYPE switch_variables


  TYPE isotope_variables
     REAL(wp) :: delta_soil ! delta13C of heterotrophic soil respiration
     REAL(wp) :: da_m_day ! slope of daytime regression deltaCa*Ca=m*Ca+b
     REAL(wp) :: da_b_day ! intercept of daytime regression deltaCa*Ca=m*Ca+b
     REAL(wp) :: da_m_night ! slope of nighttime regression deltaCa*Ca=m*Ca+b
     REAL(wp) :: da_b_night ! intercept of nighttime regression deltaCa*Ca=m*Ca+b
     REAL(wp), DIMENSION(:), ALLOCATABLE :: bigdelta ! daily flux weighted canopy discrimination after simple equation for 20 days ! [ndaysc13]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: bigdelta_long ! daily flux weighted canopy discrimination after extended equation for 20 days ! [ndaysc13]
  END TYPE isotope_variables


  TYPE water_isotope_variables
     INTEGER(i4) :: merlivat ! 0: Cappa et al. 1: Merlivat
     INTEGER(i4) :: nofracsoil ! 0: normal 1: no fractionation at soil evaporation
     INTEGER(i4) :: nofraclitter ! 0: normal 1: no fractionation at litter evaporation
     INTEGER(i4) :: nofracleaf ! 0: normal 1: no fractionation at leaf transpiration
     INTEGER(i4) :: nofracin ! 0: normal 1: rain same as initial
     INTEGER(i4) :: implicit ! 0: explicit, i.e. out of loop 1:implicit, i.e. in loop
     REAL(wp), DIMENSION(:), ALLOCATABLE :: vsmow ! ! [nwiso]
     REAL(wp), DIMENSION(:), ALLOCATABLE :: lost ! stores water lost by rounding errors during ratio calculation ! [nwiso]
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: dtheta ! initial soil delta-values ! [nsoil,nwiso]
  END TYPE water_isotope_variables

  ! structure for nitrogen variables
  TYPE nitrogen_variables
     REAL(wp) :: Ja ! extra electrons requirements for N assimilation, umol m-2 s-1
     REAL(wp) :: J_glu
     REAL(wp) :: J_Busch
     REAL(wp) :: cn_bulk ! bulk C:N ratio
     REAL(wp) :: n_supply ! field N supply, ambient for fertilization, umol m-2 s-1
     REAL(wp) :: n_mult ! set nassimilation as a multiple of glycine and serine
     REAL(wp) :: nitrate_per !fraction of nitrate in N supply
     REAL(wp) :: nitrite_per !fraction of nitrite in N supply
     REAL(wp) :: ammonia_per !fraction of ammonia in N supply
     REAL(wp) :: Busch_mol
     REAL(wp) :: nitrate_mol !mole of nitrate in N supply
     REAL(wp) :: nitrite_mol !mole of nitrite in N supply
     REAL(wp) :: ammonia_mol !mole of ammonia in N supply

     !REAL(wp) :: ammonia !fraction of ammonia in N supply


  END TYPE nitrogen_variables
  ! -----------------------------------------------------------------------------------------------
  ! Assign type variables
  TYPE(debug_variables)            :: debug
  TYPE(output_variables)           :: output
  TYPE(flux_variables)             :: flux
  TYPE(input_variables)            :: input
  TYPE(time_variables)             :: time
  TYPE(meteorology)                :: met
  TYPE(surface_resistances)        :: srf_res
  TYPE(factors)                    :: fact
  TYPE(bole_respiration_structure) :: bole
  TYPE(canopy_architecture)        :: canopy
  TYPE(non_dimensional_variables)  :: non_dim
  TYPE(boundary_layer_resistances) :: bound_lay_res
  TYPE(solar_radiation_variables)  :: solar
  TYPE(soil_variables)             :: soil
  TYPE(profile)                    :: prof
  TYPE(switch_variables)           :: iswitch
  TYPE(isotope_variables)          :: ciso
  TYPE(water_isotope_variables)    :: wiso
  TYPE(nitrogen_variables)         :: nitrogen

  ! -----------------------------------------------------------------------------------------------

CONTAINS

  ! -----------------------------------------------------------------------------------------------
  ! Allocate arrays in type variables - wrapper routine
  SUBROUTINE alloc_type_vars()

    IMPLICIT NONE

    call alloc_debug()
    call alloc_output()
    call alloc_flux()
    call alloc_input()
    call alloc_met()
    call alloc_srf_res()
    call alloc_bole()
    call alloc_canopy()
    call alloc_solar()
    call alloc_soil()
    call alloc_prof()
    call alloc_ciso()
    call alloc_wiso()

  END SUBROUTINE alloc_type_vars


  ! -----------------------------------------------------------------------------------------------
  ! Allocate arrays in type variables - individual routines
    SUBROUTINE alloc_debug()

    IMPLICIT NONE

    if (.not. allocated(debug%D1)) allocate(debug%D1(ncl))
    if (.not. allocated(debug%D2)) allocate(debug%D2(ncl))
    if (.not. allocated(debug%D3)) allocate(debug%D3(ncl))
    if (.not. allocated(debug%D4)) allocate(debug%D4(ncl))
    if (.not. allocated(debug%D5)) allocate(debug%D5(ncl))
    if (.not. allocated(debug%D6)) allocate(debug%D6(ncl))
    if (.not. allocated(debug%D7)) allocate(debug%D7(ncl))
    if (.not. allocated(debug%D8)) allocate(debug%D8(ncl))
    if (.not. allocated(debug%D9)) allocate(debug%D9(ncl))

  END SUBROUTINE alloc_debug
  SUBROUTINE alloc_output()

    IMPLICIT NONE

    if (.not. allocated(output%dLAIdz)) allocate(output%dLAIdz(ncl))

  END SUBROUTINE alloc_output


  SUBROUTINE alloc_flux()

    IMPLICIT NONE

    if (.not. allocated(flux%c_evaporation)) allocate(flux%c_evaporation(nwiso))
    if (.not. allocated(flux%c_transpiration)) allocate(flux%c_transpiration(nwiso))
    if (.not. allocated(flux%c_evapotranspiration)) allocate(flux%c_evapotranspiration(nwiso))
    if (.not. allocated(flux%soilevap)) allocate(flux%soilevap(nwiso))
    if (.not. allocated(flux%litterevap)) allocate(flux%litterevap(nwiso))
    if (.not. allocated(flux%s_evap)) allocate(flux%s_evap(nwiso))
    if (.not. allocated(flux%evapotranspiration)) allocate(flux%evapotranspiration(nwiso))
    if (.not. allocated(flux%soilinfl)) allocate(flux%soilinfl(nwiso))
    if (.not. allocated(flux%surfrun)) allocate(flux%surfrun(nwiso))

  END SUBROUTINE alloc_flux


  SUBROUTINE alloc_input()

    IMPLICIT NONE

    if (.not. allocated(input%ppt)) allocate(input%ppt(nwiso))
    if (.not. allocated(input%dppt)) allocate(input%dppt(nwiso))
    if (.not. allocated(input%dvapour)) allocate(input%dvapour(nwiso))

  END SUBROUTINE alloc_input


  SUBROUTINE alloc_met()

    IMPLICIT NONE

    if (.not. allocated(met%dispersion)) allocate(met%dispersion(ntl,ncl))

  END SUBROUTINE alloc_met


  SUBROUTINE alloc_srf_res()

    IMPLICIT NONE

    if (.not. allocated(srf_res%rcuticle)) allocate(srf_res%rcuticle(ncl))

  END SUBROUTINE alloc_srf_res


  SUBROUTINE alloc_bole()

    IMPLICIT NONE

    if (.not. allocated(bole%layer)) allocate(bole%layer(ncl))

  END SUBROUTINE alloc_bole


  SUBROUTINE alloc_canopy()

    IMPLICIT NONE

    if (.not. allocated(canopy%bdens)) allocate(canopy%bdens(10))

  END SUBROUTINE alloc_canopy


  SUBROUTINE alloc_solar()

    IMPLICIT NONE

    if (.not. allocated(solar%prob_beam)) allocate(solar%prob_beam(ncl))
    if (.not. allocated(solar%prob_shd)) allocate(solar%prob_shd(ncl))

    if (.not. allocated(solar%ir_dn)) allocate(solar%ir_dn(ncl+1))
    if (.not. allocated(solar%ir_up)) allocate(solar%ir_up(ncl+1))
    if (.not. allocated(solar%nir_dn)) allocate(solar%nir_dn(ncl+1))
    if (.not. allocated(solar%nir_up)) allocate(solar%nir_up(ncl+1))
    if (.not. allocated(solar%nir_sun)) allocate(solar%nir_sun(ncl))
    if (.not. allocated(solar%nir_shd)) allocate(solar%nir_shd(ncl))
    if (.not. allocated(solar%beam_flux_nir)) allocate(solar%beam_flux_nir(ncl+1))
    if (.not. allocated(solar%par_shd)) allocate(solar%par_shd(ncl))
    if (.not. allocated(solar%par_sun)) allocate(solar%par_sun(ncl))
    if (.not. allocated(solar%beam_flux_par)) allocate(solar%beam_flux_par(ncl+1))
    if (.not. allocated(solar%par_down)) allocate(solar%par_down(ncl+1))
    if (.not. allocated(solar%par_up)) allocate(solar%par_up(ncl+1))
    if (.not. allocated(solar%quantum_sun)) allocate(solar%quantum_sun(ncl))
    if (.not. allocated(solar%quantum_shd)) allocate(solar%quantum_shd(ncl))
    if (.not. allocated(solar%rnet_sun)) allocate(solar%rnet_sun(ncl))
    if (.not. allocated(solar%rnet_shd)) allocate(solar%rnet_shd(ncl))
    if (.not. allocated(solar%exxpdir)) allocate(solar%exxpdir(ncl))

  END SUBROUTINE alloc_solar


  SUBROUTINE alloc_soil()

    IMPLICIT NONE

    if (.not. allocated(soil%T_soil)) allocate(soil%T_soil(0:nsoil+1))
    if (.not. allocated(soil%T_soil_filter)) allocate(soil%T_soil_filter(0:nsoil+1))
    if (.not. allocated(soil%k_conductivity_soil)) allocate(soil%k_conductivity_soil(0:nsoil+1))
    if (.not. allocated(soil%cp_soil)) allocate(soil%cp_soil(0:nsoil+1))
    if (.not. allocated(soil%swp)) allocate(soil%swp(nsoil))
    if (.not. allocated(soil%swp_mm)) allocate(soil%swp_mm(nsoil))
    if (.not. allocated(soil%r)) allocate(soil%r(nsoil,nwiso))
    if (.not. allocated(soil%a)) allocate(soil%a(nsoil))
    if (.not. allocated(soil%b)) allocate(soil%b(nsoil))
    if (.not. allocated(soil%c)) allocate(soil%c(nsoil))
    if (.not. allocated(soil%dwat)) allocate(soil%dwat(nsoil,nwiso))
    if (.not. allocated(soil%theta)) allocate(soil%theta(nsoil,nwiso))
    if (.not. allocated(soil%qdrai)) allocate(soil%qdrai(nwiso))
    if (.not. allocated(soil%qinfl)) allocate(soil%qinfl(nwiso))
    if (.not. allocated(soil%qtran)) allocate(soil%qtran(nwiso))
    if (.not. allocated(soil%qseva)) allocate(soil%qseva(nwiso))
    if (.not. allocated(soil%psi)) allocate(soil%psi(nsoil))
    if (.not. allocated(soil%k_theta)) allocate(soil%k_theta(nsoil))
    if (.not. allocated(soil%theta_l)) allocate(soil%theta_l(nwiso))
    if (.not. allocated(soil%xylem)) allocate(soil%xylem(nwiso))
    if (.not. allocated(soil%rxylem)) allocate(soil%rxylem(nwiso))
    if (.not. allocated(soil%z_soil)) allocate(soil%z_soil(0:nsoil))
    if (.not. allocated(soil%bulk_density)) allocate(soil%bulk_density(nsoil))
    if (.not. allocated(soil%clay_in)) allocate(soil%clay_in(nsoil))
    if (.not. allocated(soil%sand_in)) allocate(soil%sand_in(nsoil))
    if (.not. allocated(soil%root)) allocate(soil%root(nsoil))
    if (.not. allocated(soil%clay)) allocate(soil%clay(nsoil))
    if (.not. allocated(soil%sand)) allocate(soil%sand(nsoil))
    if (.not. allocated(soil%om)) allocate(soil%om(nsoil))
    if (.not. allocated(soil%gravel)) allocate(soil%gravel(nsoil))
    if (.not. allocated(soil%theta_1500)) allocate(soil%theta_1500(nsoil))
    if (.not. allocated(soil%theta_33)) allocate(soil%theta_33(nsoil))
    if (.not. allocated(soil%theta_s33)) allocate(soil%theta_s33(nsoil))
    if (.not. allocated(soil%theta_s)) allocate(soil%theta_s(nsoil))
    if (.not. allocated(soil%psi_e)) allocate(soil%psi_e(nsoil))
    if (.not. allocated(soil%rho)) allocate(soil%rho(nsoil))
    if (.not. allocated(soil%k_s)) allocate(soil%k_s(nsoil))
    if (.not. allocated(soil%lambda)) allocate(soil%lambda(nsoil))
    if (.not. allocated(soil%big_a)) allocate(soil%big_a(nsoil))
    if (.not. allocated(soil%big_b)) allocate(soil%big_b(nsoil))
    if (.not. allocated(soil%watmin)) allocate(soil%watmin(nsoil))

  END SUBROUTINE alloc_soil


  SUBROUTINE alloc_prof()

    IMPLICIT NONE

    if (.not. allocated(prof%tair)) allocate(prof%tair(ntl))
    if (.not. allocated(prof%tair_filter)) allocate(prof%tair_filter(ntl))
    if (.not. allocated(prof%tair_filter_save)) allocate(prof%tair_filter_save(ntl))
    if (.not. allocated(prof%u)) allocate(prof%u(ntl))
    if (.not. allocated(prof%rhov_air)) allocate(prof%rhov_air(ntl,nwiso))
    if (.not. allocated(prof%rhov_air_save)) allocate(prof%rhov_air_save(ntl))
    if (.not. allocated(prof%rhov_air_filter)) allocate(prof%rhov_air_filter(ntl,nwiso))
    if (.not. allocated(prof%rhov_air_filter_save)) allocate(prof%rhov_air_filter_save(ntl))
    if (.not. allocated(prof%co2_air)) allocate(prof%co2_air(ntl))
    if (.not. allocated(prof%co2_air_filter)) allocate(prof%co2_air_filter(ntl))
    if (.not. allocated(prof%vpd_air)) allocate(prof%vpd_air(ntl))
    if (.not. allocated(prof%Gfunc_solar)) allocate(prof%Gfunc_solar(ncl))
    if (.not. allocated(prof%tleaf)) allocate(prof%tleaf(ncl))
    if (.not. allocated(prof%isopreneflux)) allocate(prof%isopreneflux(ncl))
    if (.not. allocated(prof%vcmax)) allocate(prof%vcmax(ncl))
    if (.not. allocated(prof%jmax)) allocate(prof%jmax(ncl))
    if (.not. allocated(prof%tp)) allocate(prof%tp(ncl))
    if (.not. allocated(prof%rd)) allocate(prof%rd(ncl))
    if (.not. allocated(prof%vcmaxz)) allocate(prof%vcmaxz(ncl))
    if (.not. allocated(prof%throughfall)) allocate(prof%throughfall(ncl+1,nwiso))
    if (.not. allocated(prof%cws)) allocate(prof%cws(ncl,nwiso))
    if (.not. allocated(prof%wet_coef)) allocate(prof%wet_coef(ncl))
    if (.not. allocated(prof%wet_coef_filter)) allocate(prof%wet_coef_filter(ncl))
    if (.not. allocated(prof%c13cnc)) allocate(prof%c13cnc(ntl))
    if (.not. allocated(prof%c12cnc)) allocate(prof%c12cnc(ntl))
    if (.not. allocated(prof%sour13co2)) allocate(prof%sour13co2(ncl))
    if (.not. allocated(prof%d13C)) allocate(prof%d13C(ntl))
    if (.not. allocated(prof%disc13C)) allocate(prof%disc13C(ncl))
    if (.not. allocated(prof%disc13C_long)) allocate(prof%disc13C_long(ncl))
    if (.not. allocated(prof%disc13C_ab)) allocate(prof%disc13C_ab(ncl))
    if (.not. allocated(prof%disc13C_a)) allocate(prof%disc13C_a(ncl))
    if (.not. allocated(prof%disc13C_asal)) allocate(prof%disc13C_asal(ncl))
    if (.not. allocated(prof%disc13C_b)) allocate(prof%disc13C_b(ncl))
    if (.not. allocated(prof%cs)) allocate(prof%cs(ncl))
    if (.not. allocated(prof%ci)) allocate(prof%ci(ncl))
    if (.not. allocated(prof%cc)) allocate(prof%cc(ncl))
    if (.not. allocated(prof%csca)) allocate(prof%csca(ncl))
    if (.not. allocated(prof%cica)) allocate(prof%cica(ncl))
    if (.not. allocated(prof%ccca)) allocate(prof%ccca(ncl))
    if (.not. allocated(prof%gm)) allocate(prof%gm(ncl))
    if (.not. allocated(prof%d13Cair)) allocate(prof%d13Cair(ntl))
    if (.not. allocated(prof%d13Cplant)) allocate(prof%d13Cplant(ncl))
    if (.not. allocated(prof%R13_12_air)) allocate(prof%R13_12_air(ntl))
    if (.not. allocated(prof%Rplant_sun)) allocate(prof%Rplant_sun(ncl))
    if (.not. allocated(prof%Rplant_shd)) allocate(prof%Rplant_shd(ncl))
    if (.not. allocated(prof%Rresp)) allocate(prof%Rresp(ncl))
    if (.not. allocated(prof%Rresp_sum)) allocate(prof%Rresp_sum(ncl))
    if (.not. allocated(prof%Rresp_ave)) allocate(prof%Rresp_ave(ncl))
    if (.not. allocated(prof%cnt_Rresp)) allocate(prof%cnt_Rresp(ncl))
    if (.not. allocated(prof%recycle)) allocate(prof%recycle(ntl))
    if (.not. allocated(prof%source_co2)) allocate(prof%source_co2(ncl))
    if (.not. allocated(prof%dPsdz)) allocate(prof%dPsdz(ncl))
    if (.not. allocated(prof%dPsdz_sun)) allocate(prof%dPsdz_sun(ncl))
    if (.not. allocated(prof%dPsdz_shd)) allocate(prof%dPsdz_shd(ncl))
    if (.not. allocated(prof%dGPPdz)) allocate(prof%dGPPdz(ncl))
    if (.not. allocated(prof%dGPPdz_sun)) allocate(prof%dGPPdz_sun(ncl))
    if (.not. allocated(prof%dGPPdz_shd)) allocate(prof%dGPPdz_shd(ncl))
    if (.not. allocated(prof%dGOPdz)) allocate(prof%dGOPdz(ncl))
    if (.not. allocated(prof%dGOPdz_sun)) allocate(prof%dGOPdz_sun(ncl))
    if (.not. allocated(prof%dGOPdz_shd)) allocate(prof%dGOPdz_shd(ncl))
    if (.not. allocated(prof%dHdz)) allocate(prof%dHdz(ncl))
    if (.not. allocated(prof%dLEdz)) allocate(prof%dLEdz(ncl,nwiso))
    if (.not. allocated(prof%dLEdz_sun)) allocate(prof%dLEdz_sun(ncl))
    if (.not. allocated(prof%dLEdz_shd)) allocate(prof%dLEdz_shd(ncl))
    if (.not. allocated(prof%dRNdz)) allocate(prof%dRNdz(ncl))
    if (.not. allocated(prof%dLoutdz)) allocate(prof%dLoutdz(ncl))
    if (.not. allocated(prof%dboledz)) allocate(prof%dboledz(ncl)) ! bole resp, Yuan 2018.02.12
    if (.not. allocated(prof%dRESPdz)) allocate(prof%dRESPdz(ncl))
    if (.not. allocated(prof%dRESPdz_sun)) allocate(prof%dRESPdz_sun(ncl))
    if (.not. allocated(prof%dRESPdz_shd)) allocate(prof%dRESPdz_shd(ncl))
    if (.not. allocated(prof%dStomCondz)) allocate(prof%dStomCondz(ncl))
    if (.not. allocated(prof%dStomCondz_sun)) allocate(prof%dStomCondz_sun(ncl))
    if (.not. allocated(prof%dStomCondz_shd)) allocate(prof%dStomCondz_shd(ncl))
    if (.not. allocated(prof%sun_frac)) allocate(prof%sun_frac(ncl))
    if (.not. allocated(prof%sun_tleaf)) allocate(prof%sun_tleaf(ncl))
    if (.not. allocated(prof%sun_GPP)) allocate(prof%sun_GPP(ncl))
    if (.not. allocated(prof%sun_GOP)) allocate(prof%sun_GOP(ncl))
    if (.not. allocated(prof%sun_A)) allocate(prof%sun_A(ncl))
    if (.not. allocated(prof%sun_gs)) allocate(prof%sun_gs(ncl))
    if (.not. allocated(prof%sun_gs_mol)) allocate(prof%sun_gs_mol(ncl))
    if (.not. allocated(prof%sun_rs)) allocate(prof%sun_rs(ncl))
    if (.not. allocated(prof%sun_rs_filter)) allocate(prof%sun_rs_filter(ncl))
    if (.not. allocated(prof%sun_rs_save)) allocate(prof%sun_rs_save(ncl))
    if (.not. allocated(prof%sun_rbh)) allocate(prof%sun_rbh(ncl))
    if (.not. allocated(prof%sun_rbv)) allocate(prof%sun_rbv(ncl))
    if (.not. allocated(prof%sun_rbco2)) allocate(prof%sun_rbco2(ncl))
    if (.not. allocated(prof%sun_cs)) allocate(prof%sun_cs(ncl))
    if (.not. allocated(prof%sun_ci)) allocate(prof%sun_ci(ncl))
    if (.not. allocated(prof%sun_cc)) allocate(prof%sun_cc(ncl))
    if (.not. allocated(prof%sun_disc13)) allocate(prof%sun_disc13(ncl))
    if (.not. allocated(prof%sun_disc13_ab)) allocate(prof%sun_disc13_ab(ncl))
    if (.not. allocated(prof%sun_disc13_a)) allocate(prof%sun_disc13_a(ncl))
    if (.not. allocated(prof%sun_disc13_asal)) allocate(prof%sun_disc13_asal(ncl))
    if (.not. allocated(prof%sun_disc13_b)) allocate(prof%sun_disc13_b(ncl))
    if (.not. allocated(prof%sun_disc13_long)) allocate(prof%sun_disc13_long(ncl))
    if (.not. allocated(prof%sun_csca)) allocate(prof%sun_csca(ncl))
    if (.not. allocated(prof%sun_cica)) allocate(prof%sun_cica(ncl))
    if (.not. allocated(prof%sun_ccca)) allocate(prof%sun_ccca(ncl))
    if (.not. allocated(prof%sun_lai)) allocate(prof%sun_lai(ncl))
    if (.not. allocated(prof%sun_tleaf_filter)) allocate(prof%sun_tleaf_filter(ncl))
    if (.not. allocated(prof%sun_wj)) allocate(prof%sun_wj(ncl))
    if (.not. allocated(prof%sun_wc)) allocate(prof%sun_wc(ncl))
    if (.not. allocated(prof%sun_wp)) allocate(prof%sun_wp(ncl))
    if (.not. allocated(prof%sun_alphag)) allocate(prof%sun_alphag(ncl))
    if (.not. allocated(prof%sun_alphas)) allocate(prof%sun_alphas(ncl))
    if (.not. allocated(prof%sun_tpu_coeff)) allocate(prof%sun_tpu_coeff(ncl))
    if (.not. allocated(prof%sun_resp)) allocate(prof%sun_resp(ncl))
    if (.not. allocated(prof%sun_isopreneflux)) allocate(prof%sun_isopreneflux(ncl))
    if (.not. allocated(prof%iso_sun)) allocate(prof%iso_sun(ncl))
    if (.not. allocated(prof%sun_vpd)) allocate(prof%sun_vpd(ncl))
    if (.not. allocated(prof%sun_rh)) allocate(prof%sun_rh(ncl))
    if (.not. allocated(prof%sun_LEstoma)) allocate(prof%sun_LEstoma(ncl,nwiso))
    if (.not. allocated(prof%sun_LEwet)) allocate(prof%sun_LEwet(ncl,nwiso))
    if (.not. allocated(prof%shd_frac)) allocate(prof%shd_frac(ncl))
    if (.not. allocated(prof%shd_tleaf)) allocate(prof%shd_tleaf(ncl))
    if (.not. allocated(prof%shd_GPP)) allocate(prof%shd_GPP(ncl))
    if (.not. allocated(prof%shd_GOP)) allocate(prof%shd_GOP(ncl))
    if (.not. allocated(prof%shd_A)) allocate(prof%shd_A(ncl))
    if (.not. allocated(prof%shd_gs)) allocate(prof%shd_gs(ncl))
    if (.not. allocated(prof%shd_gs_mol)) allocate(prof%shd_gs_mol(ncl))
    if (.not. allocated(prof%shd_rs)) allocate(prof%shd_rs(ncl))
    if (.not. allocated(prof%shd_rs_filter)) allocate(prof%shd_rs_filter(ncl))
    if (.not. allocated(prof%shd_rs_save)) allocate(prof%shd_rs_save(ncl))
    if (.not. allocated(prof%shd_rbh)) allocate(prof%shd_rbh(ncl))
    if (.not. allocated(prof%shd_rbv)) allocate(prof%shd_rbv(ncl))
    if (.not. allocated(prof%shd_rbco2)) allocate(prof%shd_rbco2(ncl))
    if (.not. allocated(prof%shd_cs)) allocate(prof%shd_cs(ncl))
    if (.not. allocated(prof%shd_ci)) allocate(prof%shd_ci(ncl))
    if (.not. allocated(prof%shd_cc)) allocate(prof%shd_cc(ncl))
    if (.not. allocated(prof%shd_disc13)) allocate(prof%shd_disc13(ncl))
    if (.not. allocated(prof%shd_disc13_ab)) allocate(prof%shd_disc13_ab(ncl))
    if (.not. allocated(prof%shd_disc13_a)) allocate(prof%shd_disc13_a(ncl))
    if (.not. allocated(prof%shd_disc13_asal)) allocate(prof%shd_disc13_asal(ncl))
    if (.not. allocated(prof%shd_disc13_b)) allocate(prof%shd_disc13_b(ncl))
    if (.not. allocated(prof%shd_disc13_long)) allocate(prof%shd_disc13_long(ncl))
    if (.not. allocated(prof%shd_csca)) allocate(prof%shd_csca(ncl))
    if (.not. allocated(prof%shd_cica)) allocate(prof%shd_cica(ncl))
    if (.not. allocated(prof%shd_ccca)) allocate(prof%shd_ccca(ncl))
    if (.not. allocated(prof%shd_lai)) allocate(prof%shd_lai(ncl))
    if (.not. allocated(prof%shd_tleaf_filter)) allocate(prof%shd_tleaf_filter(ncl))
    if (.not. allocated(prof%shd_wc)) allocate(prof%shd_wc(ncl))
    if (.not. allocated(prof%shd_wj)) allocate(prof%shd_wj(ncl))
    if (.not. allocated(prof%shd_wp)) allocate(prof%shd_wp(ncl))
    if (.not. allocated(prof%shd_alphag)) allocate(prof%shd_alphag(ncl))
    if (.not. allocated(prof%shd_alphas)) allocate(prof%shd_alphas(ncl))
    if (.not. allocated(prof%shd_tpu_coeff)) allocate(prof%shd_tpu_coeff(ncl))
    if (.not. allocated(prof%shd_resp)) allocate(prof%shd_resp(ncl))
    if (.not. allocated(prof%shd_isopreneflux)) allocate(prof%shd_isopreneflux(ncl))
    if (.not. allocated(prof%iso_shd)) allocate(prof%iso_shd(ncl))
    if (.not. allocated(prof%shd_vpd)) allocate(prof%shd_vpd(ncl))
    if (.not. allocated(prof%shd_rh)) allocate(prof%shd_rh(ncl))
    if (.not. allocated(prof%shd_LEstoma)) allocate(prof%shd_LEstoma(ncl,nwiso))
    if (.not. allocated(prof%shd_LEwet)) allocate(prof%shd_LEwet(ncl,nwiso))
    if (.not. allocated(prof%sun_wi)) allocate(prof%sun_wi(ncl))
    if (.not. allocated(prof%shd_wi)) allocate(prof%shd_wi(ncl))
    if (.not. allocated(prof%wa)) allocate(prof%wa(ncl))
    if (.not. allocated(prof%sun_h)) allocate(prof%sun_h(ncl))
    if (.not. allocated(prof%shd_h)) allocate(prof%shd_h(ncl))
    if (.not. allocated(prof%rs_fact)) allocate(prof%rs_fact(ncl))
    if (.not. allocated(prof%sun_gross)) allocate(prof%sun_gross(ncl))
    if (.not. allocated(prof%shd_gross)) allocate(prof%shd_gross(ncl))
    if (.not. allocated(prof%sun_LEstoma_new)) allocate(prof%sun_LEstoma_new(ncl))
    if (.not. allocated(prof%shd_LEstoma_new)) allocate(prof%shd_LEstoma_new(ncl))
    if (.not. allocated(prof%sun_LEstoma_save)) allocate(prof%sun_LEstoma_save(ncl))
    if (.not. allocated(prof%shd_LEstoma_save)) allocate(prof%shd_LEstoma_save(ncl))
    if (.not. allocated(prof%sun_alpha_k)) allocate(prof%sun_alpha_k(ncl,nwiso))
    if (.not. allocated(prof%shd_alpha_k)) allocate(prof%shd_alpha_k(ncl,nwiso))
    if (.not. allocated(prof%sun_alpha_equ)) allocate(prof%sun_alpha_equ(ncl,nwiso))
    if (.not. allocated(prof%shd_alpha_equ)) allocate(prof%shd_alpha_equ(ncl,nwiso))
    if (.not. allocated(prof%rvapour)) allocate(prof%rvapour(ntl,nwiso))
    if (.not. allocated(prof%dvapour)) allocate(prof%dvapour(ntl,nwiso))
    if (.not. allocated(prof%sun_peclet)) allocate(prof%sun_peclet(ncl,nwiso))
    if (.not. allocated(prof%shd_peclet)) allocate(prof%shd_peclet(ncl,nwiso))
    if (.not. allocated(prof%sun_fem)) allocate(prof%sun_fem(ncl,nwiso))
    if (.not. allocated(prof%shd_fem)) allocate(prof%shd_fem(ncl,nwiso))
    if (.not. allocated(prof%sun_craig)) allocate(prof%sun_craig(ncl,nwiso))
    if (.not. allocated(prof%shd_craig)) allocate(prof%shd_craig(ncl,nwiso))
    if (.not. allocated(prof%sun_leafwater_e)) allocate(prof%sun_leafwater_e(ncl,nwiso))
    if (.not. allocated(prof%shd_leafwater_e)) allocate(prof%shd_leafwater_e(ncl,nwiso))
    if (.not. allocated(prof%sun_leafwater_e_old)) allocate(prof%sun_leafwater_e_old(ncl,nwiso))
    if (.not. allocated(prof%shd_leafwater_e_old)) allocate(prof%shd_leafwater_e_old(ncl,nwiso))
    if (.not. allocated(prof%sun_leafwater)) allocate(prof%sun_leafwater(ncl,nwiso))
    if (.not. allocated(prof%shd_leafwater)) allocate(prof%shd_leafwater(ncl,nwiso))
    if (.not. allocated(prof%sun_trans_rtrans)) allocate(prof%sun_trans_rtrans(ncl,nwiso))
    if (.not. allocated(prof%shd_trans_rtrans)) allocate(prof%shd_trans_rtrans(ncl,nwiso))
    if (.not. allocated(prof%sun_rtrans)) allocate(prof%sun_rtrans(ncl,nwiso))
    if (.not. allocated(prof%shd_rtrans)) allocate(prof%shd_rtrans(ncl,nwiso))
    if (.not. allocated(prof%ht)) allocate(prof%ht(ncl))
    if (.not. allocated(prof%dLAIdz)) allocate(prof%dLAIdz(ncl))
    if (.not. allocated(prof%dPAIdz)) allocate(prof%dPAIdz(ncl))
    if (.not. allocated(prof%dWAIdz)) allocate(prof%dWAIdz(ncl)) ! YUAN 2018.03.01
    if (.not. allocated(prof%Gfunc_sky)) allocate(prof%Gfunc_sky(ncl,nl))
 ! oxygen module Yuan 2018.01.17
    if (.not. allocated(prof%O2_air)) allocate(prof%O2_air(ntl))
    if (.not. allocated(prof%O2_air_filter)) allocate(prof%O2_air_filter(ntl))
    if (.not. allocated(prof%O2_soil)) allocate(prof%O2_soil(ntl))
    if (.not. allocated(prof%CO2_soil)) allocate(prof%CO2_soil(ntl))
    if (.not. allocated(prof%O2_disp)) allocate(prof%O2_disp(ntl))
    if (.not. allocated(prof%CO2_disp)) allocate(prof%CO2_disp(ntl))
    if (.not. allocated(prof%source_O2)) allocate(prof%source_O2(ncl))
    if (.not. allocated(prof%sun_psn_O2)) allocate(prof%sun_psn_O2(ncl))
    if (.not. allocated(prof%shd_psn_O2)) allocate(prof%shd_psn_O2(ncl))
    if (.not. allocated(prof%dPsdz_O2_sun)) allocate(prof%dPsdz_O2_sun(ncl))
    if (.not. allocated(prof%dPsdz_O2_shd)) allocate(prof%dPsdz_O2_shd(ncl))
    if (.not. allocated(prof%dPsdz_O2)) allocate(prof%dPsdz_O2(ncl))
    if (.not. allocated(prof%sun_resp_O2)) allocate(prof%sun_resp_O2(ncl))
    if (.not. allocated(prof%shd_resp_O2)) allocate(prof%shd_resp_O2(ncl))
    if (.not. allocated(prof%dRESPdz_O2_sun)) allocate(prof%dRESPdz_O2_sun(ncl))
    if (.not. allocated(prof%dRESPdz_O2_shd)) allocate(prof%dRESPdz_O2_shd(ncl))
    if (.not. allocated(prof%dRESPdz_O2)) allocate(prof%dRESPdz_O2(ncl))
    if (.not. allocated(prof%sun_NO3)) allocate(prof%sun_NO3(ncl))
    if (.not. allocated(prof%shd_NO3)) allocate(prof%shd_NO3(ncl))
    if (.not. allocated(prof%sun_NO2)) allocate(prof%sun_NO2(ncl))
    if (.not. allocated(prof%shd_NO2)) allocate(prof%shd_NO2(ncl))
    if (.not. allocated(prof%sun_NH4)) allocate(prof%sun_NH4(ncl))
    if (.not. allocated(prof%shd_NH4)) allocate(prof%shd_NH4(ncl))
    if (.not. allocated(prof%Ja_sun)) allocate(prof%Ja_sun(ncl))
    if (.not. allocated(prof%Jglu_sun)) allocate(prof%Jglu_sun(ncl))
    if (.not. allocated(prof%JBusch_sun)) allocate(prof%JBusch_sun(ncl))
    if (.not. allocated(prof%Ja_shd)) allocate(prof%Ja_shd(ncl))
    if (.not. allocated(prof%Jglu_shd)) allocate(prof%Jglu_shd(ncl))
    if (.not. allocated(prof%JBusch_shd)) allocate(prof%JBusch_shd(ncl))
    if (.not. allocated(prof%jphoton_sun)) allocate(prof%jphoton_sun(ncl))
    if (.not. allocated(prof%jphoton_shd)) allocate(prof%jphoton_shd(ncl))
    if (.not. allocated(prof%sun_quad)) allocate(prof%sun_quad(ncl))
    if (.not. allocated(prof%shd_quad)) allocate(prof%shd_quad(ncl))
    if (.not. allocated(prof%sun_ABusch)) allocate(prof%sun_ABusch(ncl))
    if (.not. allocated(prof%shd_ABusch)) allocate(prof%shd_ABusch(ncl))


    if (.not. allocated(prof%gpp_O2)) allocate(prof%gpp_O2(ncl))
!    if (.not. allocated(prof%rd_O2)) allocate(prof%rd_O2(ncl))
    if (.not. allocated(prof%RQ_sun)) allocate(prof%RQ_sun(ncl))
    if (.not. allocated(prof%RQ_shd)) allocate(prof%RQ_shd(ncl))
    if (.not. allocated(prof%ROC_sun)) allocate(prof%ROC_sun(ncl))
    if (.not. allocated(prof%ROC_shd)) allocate(prof%ROC_shd(ncl))
    if (.not. allocated(prof%ROC_layer)) allocate(prof%ROC_layer(ncl))
    if (.not. allocated(prof%ROC_leaf_air)) allocate(prof%ROC_leaf_air(ncl))
    if (.not. allocated(prof%ROC_bole_air)) allocate(prof%ROC_bole_air(ncl))

  END SUBROUTINE alloc_prof


  SUBROUTINE alloc_ciso()

    IMPLICIT NONE

    if (.not. allocated(ciso%bigdelta)) allocate(ciso%bigdelta(ndaysc13))
    if (.not. allocated(ciso%bigdelta_long)) allocate(ciso%bigdelta_long(ndaysc13))

  END SUBROUTINE alloc_ciso


  SUBROUTINE alloc_wiso()

    IMPLICIT NONE

    if (.not. allocated(wiso%vsmow)) allocate(wiso%vsmow(nwiso))
    if (.not. allocated(wiso%lost)) allocate(wiso%lost(nwiso))
    if (.not. allocated(wiso%dtheta)) allocate(wiso%dtheta(nsoil,nwiso))

  END SUBROUTINE alloc_wiso



  ! -----------------------------------------------------------------------------------------------
  ! Zero all type variables - wrapper routines
  SUBROUTINE zero_new_timestep()
    ! This calls the zeroing routines after each time step.
    ! The individual routines leave out the variables
    ! that are passed on to the next timestep.
    ! If you want to zero all variables, call zero_initial() as well.

    IMPLICIT NONE

    call zero_new_timestep_debug()
    call zero_new_timestep_output()
    call zero_new_timestep_flux()
    call zero_new_timestep_input()
    call zero_new_timestep_time()
    call zero_new_timestep_met()
    call zero_new_timestep_srf_res()
    call zero_new_timestep_fact()
    call zero_new_timestep_bole()
    call zero_new_timestep_bound_lay_res()
    call zero_new_timestep_solar()
    call zero_new_timestep_soil()
    call zero_new_timestep_prof()

  END SUBROUTINE zero_new_timestep


  SUBROUTINE zero_initial()
    ! This calls the routines that zero the variables that are left out
    ! in the individual zero_new_timetsp routines.

    IMPLICIT NONE

    CALL zero_initial_output()
    CALL zero_initial_input()
    CALL zero_initial_met()
    CALL zero_initial_srf_res()
    CALL zero_initial_canopy()
    CALL zero_initial_time()
    CALL zero_initial_solar()
    CALL zero_initial_soil()
    CALL zero_initial_prof()
    CALL zero_initial_ciso()
    CALL zero_initial_wiso()
    CALL zero_initial_iswitch()
    CALL zero_initial_non_dim()
    CALL zero_initial_nitrogen()

  END SUBROUTINE zero_initial

  ! -----------------------------------------------------------------------------------------------
  ! Zero all type variables for new timestep - individual routines
  SUBROUTINE zero_new_timestep_debug()

  USE constants, ONLY: undef, zero

  IMPLICIT NONE

  debug%R1 = undef
  debug%R2 = undef
  debug%R3 = undef
  debug%R4 = undef
  debug%R5 = undef
  debug%R6 = undef
  debug%R7 = undef
  debug%R8 = undef
  debug%R9 = undef

  debug%in1 = int(zero)
  debug%in2 = int(zero)
  debug%in3 = int(zero)
  debug%in4 = int(zero)
  debug%in5 = int(zero)
  debug%in6 = int(zero)
  debug%in7 = int(zero)
  debug%in8 = int(zero)
  debug%in9 = int(zero)

  END SUBROUTINE zero_new_timestep_debug


  SUBROUTINE zero_new_timestep_output()

    USE constants, ONLY: zero

    IMPLICIT NONE

    output%c1 = zero
    output%c2 = zero
    output%c3 = zero
    output%c4 = zero
    output%c5 = zero
    output%c6 = zero
    output%c7 = zero
    output%c8 = zero
    output%c9 = zero
    output%c10 = zero
    output%netrad = zero
    output%sumrn = zero
    output%sumlout = zero
    output%sumh = zero
    output%sumle = zero
    output%can_ps_mol = zero
    output%can_gpp = zero
    output%c_transpiration_mole = zero
    output%canresp = zero
    output%sumksi = zero
    output%tleaf_mean = zero
    output%tavg_sun = zero
    output%tavg_shd = zero
    output%ave_disc13 = zero
    output%ave_disc13_long = zero
    output%ave_cica = zero
    output%ave_gs = zero
    output%ave_daC13 = zero
    output%ave_disc13_a = zero
    output%ave_disc13_ab = zero
    output%ave_disc13_asal = zero
    output%ave_disc13_b = zero
    output%ave_csca = zero
    output%ave_ccca = zero
    output%ave_gm = zero
    output%isoprene_efflux = zero
    output%diff_par = zero
    output%rnet_soil = zero
    ! output%sumfc = zero
    ! output%sumevap = zero
    ! output%sumsens = zero
    ! output%sumps = zero
    ! output%sumgpp = zero
    ! output%sumpar = zero
    ! output%sumnet = zero
    ! output%sumbole = zero
    ! output%sumsoil = zero
    ! output%sumta = zero
    ! output%sumgs = zero
    ! output%sumresp = zero
    ! output%sumTleaf = zero
    ! output%sumF13C = zero
    ! output%sumdisc13C = zero
    ! output%sumdisc13C_long_day = zero
    output%dLAIdz = zero

  END SUBROUTINE zero_new_timestep_output


  SUBROUTINE zero_new_timestep_flux()

    USE constants, ONLY: zero

    IMPLICIT NONE

    flux%photosyn = zero
    flux%nee = zero
    flux%c_evaporation = zero
    flux%c_transpiration = zero
    flux%c_evapotranspiration = zero
    flux%soilevap = zero
    flux%litterevap = zero
    flux%s_evap = zero
    flux%evapotranspiration = zero
    flux%soilinfl = zero
    flux%surfrun = zero
    flux%albedo = zero

  END SUBROUTINE zero_new_timestep_flux


  SUBROUTINE zero_new_timestep_input()

    USE constants, ONLY: zero

    IMPLICIT NONE

    input%dayy = 0
    input%hhrr = zero
    input%ta = zero
    input%rglobal = zero
    input%parin = zero
    input%pardif = zero
    ! input%lai = zero
    ! input%lai_up = zero
    ! input%lai_down = zero
    input%ea = zero
    input%wnd = zero
    input%ppt = zero
    input%dppt = zero
    input%dvapour = zero
    input%co2air = zero
    input%press_mb = zero
    input%tsoil = zero
    input%soilmoisture = zero
    input%flag = 0
    input%longwave = zero
    input%d13CO2 = zero
    input%d18CO2 = zero
    input%o2air = zero ! atom o2 ppm Yuan 2018.02.14
    input%ER = zero

  END SUBROUTINE zero_new_timestep_input


  SUBROUTINE zero_new_timestep_time()

    USE constants, ONLY: zero

    IMPLICIT NONE

    time%local_time = zero
    time%daytime = 0
    time%doy = 0
    time%days = 0
    time%jdold = 0
    time%count = 0
    !time%year0 = 0
    time%year = 0
    !time%leafout = 0
    !time%leaffull = 0
    !time%leaffall = 0
    !time%leaffallcomplete = 0
    !time%lai = zero
    !time%time_step = zero

  END SUBROUTINE zero_new_timestep_time


  SUBROUTINE zero_new_timestep_met()

    USE constants, ONLY: zero

    IMPLICIT NONE

    met%ustar = zero
    met%ustar_filter = zero
    met%K = zero
    met%K_filter = zero
    met%phim = zero
    met%rhova_g = zero
    met%rhova_kg = zero
    ! met%H = zero
    ! met%H_filter = zero
    met%air_density = zero
    met%T_Kelvin = zero
    met%zl = zero
    met%press_kpa = zero
    met%press_bars = zero
    met%press_Pa = zero
    met%pstat273 = zero
    met%air_density_mole = zero
    met%relative_humidity = zero
    ! met%dispersion = zero

  END SUBROUTINE zero_new_timestep_met


  SUBROUTINE zero_new_timestep_srf_res()

    USE constants, ONLY: zero

    IMPLICIT NONE

    srf_res%fdrought   = zero
    srf_res%rcuticle   = zero
    !srf_res%fthreshold = zero
    !srf_res%fslope     = zero

  END SUBROUTINE zero_new_timestep_srf_res


  SUBROUTINE zero_new_timestep_fact()

    USE constants, ONLY: zero

    IMPLICIT NONE

    fact%latent   = zero
    fact%heatcoef = zero
    fact%a_filt   = zero
    fact%co2      = zero

  END SUBROUTINE zero_new_timestep_fact


  SUBROUTINE zero_new_timestep_bole()

    USE constants, ONLY: zero

    IMPLICIT NONE

    bole%factor           = zero
    bole%respiration_mole = zero
    bole%respiration_mg   = zero
    bole%calc             = zero
    bole%layer            = zero

  END SUBROUTINE zero_new_timestep_bole


  SUBROUTINE zero_new_timestep_bound_lay_res()

    USE constants, ONLY: zero

    IMPLICIT NONE

    bound_lay_res%vapor = zero
    bound_lay_res%heat  = zero
    bound_lay_res%co2   = zero

  END SUBROUTINE zero_new_timestep_bound_lay_res


  SUBROUTINE zero_new_timestep_solar()

    USE constants, ONLY: zero

    IMPLICIT NONE

    solar%prob_beam = zero
    solar%prob_shd = zero
    solar%ir_dn = zero
    solar%ir_up = zero
    solar%nir_beam = zero
    solar%nir_diffuse = zero
    solar%nir_total = zero
    solar%nir_dn = zero
    solar%nir_up = zero
    solar%nir_sun = zero
    solar%nir_shd = zero
    solar%beam_flux_nir = zero
    solar%par_diffuse = zero
    solar%par_beam = zero
    solar%par_shd = zero
    solar%par_sun = zero
    solar%beam_flux_par = zero
    solar%par_down = zero
    solar%par_up = zero
    solar%quantum_sun = zero
    solar%quantum_shd = zero
    solar%rnet_sun = zero
    solar%rnet_shd = zero
    solar%beta_rad = zero
    solar%sine_beta = zero
    solar%beta_deg = zero
    solar%ratrad = zero
    !solar%nir_reflect = zero
    !solar%nir_trans = zero
    !solar%nir_soil_refl_dry = zero
    !solar%nir_soil_refl = zero
    !solar%nir_absorbed = zero
    !solar%par_absorbed = zero
    !solar%par_reflect = zero
    !solar%par_trans = zero
    !solar%par_soil_refl_dry = zero
    !solar%par_soil_refl = zero
    !solar%exxpdir = zero
    !solar%ratradnoon = zero

  END SUBROUTINE zero_new_timestep_solar


  SUBROUTINE zero_new_timestep_soil()

    USE constants, ONLY: zero

    IMPLICIT NONE

    !soil%T_soil = zero
    !soil%T_soil_filter = zero
    soil%k_conductivity_soil = zero
    soil%cp_soil = zero
    soil%T_Kelvin = zero
    soil%T_air = zero
    !soil%tsrf = zero
    !soil%tsrf_filter = zero
    soil%rs = zero
    soil%rb = zero
    soil%T_15cm = zero
    soil%lout = zero
    !soil%soilevap = zero
    !soil%soilevap_filter = zero
    !soil%litterevap = zero
    !soil%litterevap_filter = zero
    soil%litterevap_save = zero
    soil%maxlitterevap = zero
    soil%evap = zero
    soil%heat = zero
    soil%rnet = zero
    soil%gsoil = zero
    soil%respiration_mole = zero
    soil%respiration_mg = zero
    soil%base_respiration = zero
    soil%respiration_auto = zero
    soil%respiration_hetero = zero
    soil%resp_13 = zero
    soil%SR_ref_temp = zero
    soil%swp = zero
    soil%swp_mm = zero
    soil%r = zero
    soil%a = zero
    soil%b = zero
    soil%c = zero
    soil%dwat = zero
    !soil%theta = zero
    soil%qdrai = zero
    soil%qinfl = zero
    soil%soil_mm = zero
    !soil%soil_mm_root = zero
    soil%soil_mm_50 = zero
    soil%qtran = zero
    soil%qseva = zero
    ! soil%psi = zero
    soil%k_theta = zero
    !soil%theta_l = zero
    !soil%T_l = zero
    !soil%T_l_filter = zero
    soil%rl = zero
    soil%c_litterevap = zero
    soil%c_litterevap_save = zero
    soil%maxlitter = 0
    soil%latent = zero
    soil%lecoef = zero
    soil%rh_soil = zero
    soil%rh_litter = zero
    soil%xylem = zero
    soil%rxylem = zero
    ! soil%z_soil = zero
    ! soil%bulk_density = zero
    ! soil%d = zero
    ! soil%z0 = zero
    ! soil%T_base = zero
    ! soil%temperature_dt = zero
    ! soil%temperature_mtime = 0
    ! soil%moisture_dt = zero
    ! soil%moisture_mtime = 0
    ! soil%qinfl_max = zero
    ! soil%clay_in = zero
    ! soil%sand_in = zero
    ! soil%root = zero
    ! soil%clay = zero
    ! soil%sand = zero
    ! soil%om = zero
    ! soil%gravel = zero
    ! soil%root_factor = zero
    ! soil%sand_factor = zero
    ! soil%clay_factor = zero
    ! soil%soil_mm_33_root = zero
    ! soil%soil_mm_1500_root = zero
    ! soil%saxton = 0
    ! soil%theta_1500 = zero
    ! soil%theta_33 = zero
    ! soil%theta_s33 = zero
    ! soil%theta_s = zero
    ! soil%psi_e = zero
    ! soil%rho = zero
    ! soil%k_s = zero
    ! soil%lambda = zero
    ! soil%big_a = zero
    ! soil%big_b = zero
    ! soil%z_litter = zero
    ! soil%theta_ls = zero
    ! soil%theta_l33 = zero
    ! soil%n_l = zero
    ! soil%rho_l = zero
    ! soil%camillo = 0
    ! soil%watmin = zero
    ! soil%drain0 = 0
    ! soil%lost0 = 0

  END SUBROUTINE zero_new_timestep_soil


  SUBROUTINE zero_new_timestep_prof()

    USE constants, ONLY: zero

    IMPLICIT NONE

    !prof%tair = zero
    !prof%tair_filter = zero
    prof%tair_filter_save = zero
    prof%u = zero
    !prof%rhov_air = zero
    prof%rhov_air_save = zero
    !prof%rhov_air_filter = zero
    prof%rhov_air_filter_save = zero
    !prof%co2_air = zero
    !prof%co2_air_filter = zero
    prof%vpd_air = zero
    prof%Gfunc_solar = zero
    prof%tleaf = zero
    prof%isopreneflux = zero
    prof%vcmax = zero
    prof%jmax = zero
    prof%tp = zero
    prof%rd = zero
    prof%vcmaxz = zero
    prof%throughfall = zero
    !prof%cws = zero
    prof%wet_coef = zero
    prof%wet_coef_filter = zero
    prof%c13cnc = zero
    prof%c12cnc = zero
    prof%sour13co2 = zero
    prof%d13C = zero
    prof%disc13C = zero
    prof%disc13C_long = zero
    prof%disc13C_ab = zero
    prof%disc13C_a = zero
    prof%disc13C_asal = zero
    prof%disc13C_b = zero
    prof%cs = zero
    prof%ci = zero
    prof%cc = zero
    prof%csca = zero
    prof%cica = zero
    prof%ccca = zero
    prof%gm = zero
    !prof%d13Cair = zero
    prof%d13Cplant = zero
    !prof%R13_12_air = zero
    prof%Rplant_sun = zero
    prof%Rplant_shd = zero
    prof%Rresp = zero
    prof%Rresp_sum = zero
    prof%Rresp_ave = zero
    prof%cnt_Rresp = 0
    prof%recycle = zero
    prof%source_co2 = zero
    prof%dPsdz = zero
    prof%dPsdz_sun = zero
    prof%dPsdz_shd = zero
    prof%dPsdz_O2 = zero
    prof%dPsdz_O2_sun = zero
    prof%dPsdz_O2_shd = zero
    prof%dGPPdz = zero
    prof%dGPPdz_sun = zero
    prof%dGPPdz_shd = zero
    prof%dGOPdz = zero
    prof%dGOPdz_sun = zero
    prof%dGOPdz_shd = zero
    prof%dHdz = zero
    prof%dLEdz = zero
    prof%dLEdz_sun = zero
    prof%dLEdz_shd = zero
    prof%dRNdz = zero
    prof%dLoutdz = zero
    prof%dboledz = zero ! Yuan 2018.02.12
    prof%dRESPdz = zero
    prof%dRESPdz_sun = zero
    prof%dRESPdz_shd = zero
    prof%dRESPdz_O2 = zero
    prof%dRESPdz_O2_sun = zero
    prof%dRESPdz_O2_shd = zero
    prof%dStomCondz = zero
    prof%dStomCondz_sun = zero
    prof%dStomCondz_shd = zero
    prof%sun_frac = zero
    !prof%sun_tleaf = zero
    prof%sun_GPP = zero
    prof%sun_GOP = zero
    prof%sun_A = zero
    prof%sun_gs = zero
    prof%sun_gs_mol = zero
    !prof%sun_rs = zero
    !prof%sun_rs_filter = zero
    prof%sun_rs_save = zero !Yuan swich off 2018.05.15
    prof%sun_rbh = zero
    prof%sun_rbv = zero
    prof%sun_rbco2 = zero
    prof%sun_cs = zero
    prof%sun_ci = zero
    prof%sun_cc = zero
    prof%sun_disc13 = zero
    prof%sun_disc13_ab = zero
    prof%sun_disc13_a = zero
    prof%sun_disc13_asal = zero
    prof%sun_disc13_b = zero
    prof%sun_disc13_long = zero
    prof%sun_csca = zero
    prof%sun_cica = zero
    prof%sun_ccca = zero
    prof%sun_lai = zero
    !prof%sun_tleaf_filter = zero
    prof%sun_wj = zero
    prof%sun_wc = zero
    prof%sun_wp = zero
    prof%sun_alphag = zero
    prof%sun_alphas = zero
    prof%sun_tpu_coeff = zero
    prof%sun_resp = zero
    prof%sun_isopreneflux = zero
    prof%iso_sun = zero
    prof%sun_vpd = zero
    prof%sun_rh = zero
    prof%sun_LEstoma = zero
    prof%sun_LEwet = zero
    prof%shd_frac = zero
    !prof%shd_tleaf = zero
    prof%shd_GPP = zero
    prof%shd_GOP = zero
    prof%shd_A = zero
    prof%shd_gs = zero
    prof%shd_gs_mol = zero
    !prof%shd_rs = zero
    !prof%shd_rs_filter = zero
    prof%shd_rs_save = zero !Yuan swich off 2018.05.15
    prof%shd_rbh = zero
    prof%shd_rbv = zero
    prof%shd_rbco2 = zero
    prof%shd_cs = zero
    prof%shd_ci = zero
    prof%shd_cc = zero
    prof%shd_disc13 = zero
    prof%shd_disc13_ab = zero
    prof%shd_disc13_a = zero
    prof%shd_disc13_asal = zero
    prof%shd_disc13_b = zero
    prof%shd_disc13_long = zero
    prof%shd_csca = zero
    prof%shd_cica = zero
    prof%shd_ccca = zero
    prof%shd_lai = zero
    !prof%shd_tleaf_filter = zero
    prof%shd_wc = zero
    prof%shd_wj = zero
    prof%shd_wp = zero
    prof%shd_alphag = zero
    prof%shd_alphas = zero
    prof%shd_tpu_coeff = zero
    prof%shd_resp = zero
    prof%shd_isopreneflux = zero
    prof%iso_shd = zero
    prof%shd_vpd = zero
    prof%shd_rh = zero
    prof%shd_LEstoma = zero
    prof%shd_LEwet = zero
    prof%sun_wi = zero
    prof%shd_wi = zero
    prof%wa = zero
    prof%sun_h = zero
    prof%shd_h = zero
    prof%rs_fact = zero
    prof%sun_gross = zero
    prof%shd_gross = zero
    prof%sun_LEstoma_new = zero
    prof%shd_LEstoma_new = zero
    prof%sun_LEstoma_save = zero
    prof%shd_LEstoma_save = zero
    prof%sun_alpha_k = zero
    prof%shd_alpha_k = zero
    prof%sun_alpha_equ = zero
    prof%shd_alpha_equ = zero
    !prof%rvapour = zero
    prof%sun_peclet = zero
    prof%shd_peclet = zero
    prof%sun_fem = zero
    prof%shd_fem = zero
    prof%sun_craig = zero
    prof%shd_craig = zero
    prof%sun_leafwater_e = zero
    prof%shd_leafwater_e = zero
    !prof%sun_leafwater_e_old = zero
    !prof%shd_leafwater_e_old = zero
    prof%sun_leafwater = zero
    prof%shd_leafwater = zero
    prof%sun_trans_rtrans = zero
    prof%shd_trans_rtrans = zero
    prof%sun_rtrans = zero
    prof%shd_rtrans = zero
    !prof%ht = zero
    !prof%dLAIdz = zero
    !prof%dPAIdz = zero
    !prof%Gfunc_sky = zero
    ! oxygen module: Yuan 2018.01.17
    !prof%O2_air = zero
    prof%source_O2 = zero
    prof%O2_soil = zero
    prof%CO2_soil = zero
    prof%O2_disp = zero
    prof%CO2_disp = zero
    prof%gpp_O2 = zero
!    prof%rd_O2 = zero
    prof%ROC_layer = zero
    prof%RQ_sun = zero
    prof%RQ_shd = zero
    !prof%ROC_leaf_air = zero
    !prof%ROC_bole_air = zero

  END SUBROUTINE zero_new_timestep_prof


  ! -----------------------------------------------------------------------------------------------
  ! Zero all type variables for initial time step - individual routines
  SUBROUTINE zero_initial_output()

    USE constants, ONLY: zero

    IMPLICIT NONE

    output%sumfc               = zero
    output%sumevap             = zero
    output%sumsens             = zero
    output%sumps               = zero
    output%sumgpp              = zero
    output%sumpar              = zero
    output%sumnet              = zero
    output%sumbole             = zero
    output%sumsoil             = zero
    output%sumta               = zero
    output%sumgs               = zero
    output%sumresp             = zero
    output%sumTleaf            = zero
    output%sumF13C             = zero
    output%sumdisc13C          = zero
    output%sumdisc13C_long_day = zero
    ! total and net oxygen flux Yuan 2018.01.30
    output%sumo                = zero
    output%sumneto             = zero
    output%sumresp_o           = zero
    output%hour_canrespo       = zero
    output%houro               = zero
    output%hourneto            = zero
  END SUBROUTINE zero_initial_output


  SUBROUTINE zero_initial_input()

    USE constants, ONLY: zero

    IMPLICIT NONE

    input%lai = zero
    input%lai_up = zero
    input%lai_down = zero

  END SUBROUTINE zero_initial_input


  SUBROUTINE zero_initial_met()

    USE constants, ONLY: zero

    IMPLICIT NONE

    met%H        = zero
    met%H_filter = zero

  END SUBROUTINE zero_initial_met


  SUBROUTINE zero_initial_srf_res()

    USE constants, ONLY: zero

    IMPLICIT NONE

    srf_res%fthreshold = zero
    srf_res%fslope     = zero

  END SUBROUTINE zero_initial_srf_res


  SUBROUTINE zero_initial_canopy()

    USE constants, ONLY: zero

    IMPLICIT NONE

    canopy%bdens = zero

  END SUBROUTINE zero_initial_canopy


  SUBROUTINE zero_initial_time()

    USE constants, ONLY: zero

    IMPLICIT NONE

    time%year0            = 0
    time%leafout          = 0
    time%leaffull         = 0
    time%leaffall         = 0
    time%leaffallcomplete = 0
    time%lai              = zero
    time%time_step        = zero

  END SUBROUTINE zero_initial_time


  SUBROUTINE zero_initial_solar()

    USE constants, ONLY: zero

    IMPLICIT NONE

    solar%nir_reflect = zero
    solar%nir_trans = zero
    solar%nir_soil_refl_dry = zero
    solar%nir_soil_refl = zero
    solar%nir_absorbed = zero
    solar%par_absorbed = zero
    solar%par_reflect = zero
    solar%par_trans = zero
    solar%par_soil_refl_dry = zero
    solar%par_soil_refl = zero
    solar%exxpdir = zero
    solar%ratradnoon = zero

  END SUBROUTINE zero_initial_solar


  SUBROUTINE zero_initial_soil()

    USE constants, ONLY: zero

    IMPLICIT NONE

    soil%T_soil = zero
    soil%T_soil_filter = zero
    soil%tsrf = zero
    soil%tsrf_filter = zero
    soil%soilevap = zero
    soil%soilevap_filter = zero
    soil%litterevap = zero
    soil%litterevap_filter = zero
    soil%theta = zero
    soil%soil_mm_root = zero
    soil%psi = zero
    soil%theta_l = zero
    soil%T_l = zero
    soil%T_l_filter = zero
    soil%z_soil = zero
    soil%bulk_density = zero
    soil%d = zero
    soil%z0 = zero
    soil%T_base = zero
    soil%temperature_dt = zero
    soil%temperature_mtime = 0
    soil%moisture_dt = zero
    soil%moisture_mtime = 0
    soil%qinfl_max = zero
    soil%clay_in = zero
    soil%sand_in = zero
    soil%root = zero
    soil%clay = zero
    soil%sand = zero
    soil%om = zero
    soil%gravel = zero
    soil%root_factor = zero
    soil%sand_factor = zero
    soil%clay_factor = zero
    soil%soil_mm_33_root = zero
    soil%soil_mm_1500_root = zero
    soil%saxton = 0
    soil%theta_1500 = zero
    soil%theta_33 = zero
    soil%theta_s33 = zero
    soil%theta_s = zero
    soil%psi_e = zero
    soil%rho = zero
    soil%k_s = zero
    soil%lambda = zero
    soil%big_a = zero
    soil%big_b = zero
    soil%z_litter = zero
    soil%theta_ls = zero
    soil%theta_l33 = zero
    soil%n_l = zero
    soil%rho_l = zero
    soil%camillo = 0
    soil%watmin = zero
    soil%drain0 = 0
    soil%lost0 = 0

  END SUBROUTINE zero_initial_soil


  SUBROUTINE zero_initial_prof()

    USE constants, ONLY: zero

    IMPLICIT NONE

    prof%tair = zero
    prof%tair_filter = zero
    prof%rhov_air = zero
    prof%rhov_air_filter = zero
    prof%co2_air = zero
    prof%co2_air_filter = zero
    prof%O2_air = zero
    prof%O2_air_filter = zero
    prof%cws = zero
    prof%d13Cair = zero
    prof%R13_12_air = zero
    prof%sun_tleaf = zero
    prof%sun_rs = zero
    prof%sun_rs_filter = zero
    prof%sun_tleaf_filter = zero
    prof%shd_tleaf = zero
    prof%shd_rs = zero
    prof%shd_rs_filter = zero
    prof%shd_tleaf_filter = zero
    prof%rvapour = zero
    prof%dvapour = zero
    prof%sun_leafwater_e_old = zero
    prof%shd_leafwater_e_old = zero
    prof%ht = zero
    prof%dLAIdz = zero
    prof%dPAIdz = zero
    prof%dWAIdz = zero
    prof%Gfunc_sky = zero

  END SUBROUTINE zero_initial_prof


  SUBROUTINE zero_initial_ciso()

    USE constants, ONLY: zero

    IMPLICIT NONE

    ciso%delta_soil    = zero
    ciso%da_m_day      = zero
    ciso%da_b_day      = zero
    ciso%da_m_night    = zero
    ciso%da_b_night    = zero
    ciso%bigdelta      = zero
    ciso%bigdelta_long = zero

  END SUBROUTINE zero_initial_ciso


  SUBROUTINE zero_initial_wiso()

    USE constants, ONLY: zero

    IMPLICIT NONE

    wiso%merlivat     = 0
    wiso%nofracsoil   = 0
    wiso%nofraclitter = 0
    wiso%nofracleaf   = 0
    wiso%nofracin     = 0
    wiso%implicit     = 0
    wiso%vsmow        = zero
    wiso%lost         = zero
    wiso%dtheta       = zero

  END SUBROUTINE zero_initial_wiso


  SUBROUTINE zero_initial_iswitch()

    IMPLICIT NONE

    iswitch%soil_resp_temp    = 0
    iswitch%ball              = 0
    iswitch%isoprene          = 0
    iswitch%d13c              = 0
    iswitch%wiso              = 0
    iswitch%bethy_resp        = 0
    iswitch%no_neg_water_flux = 0

  END SUBROUTINE zero_initial_iswitch


  SUBROUTINE zero_initial_non_dim()

    USE constants, ONLY: zero

    IMPLICIT NONE

    non_dim%pr       = zero
    non_dim%pr33     = zero
    non_dim%sc       = zero
    non_dim%sc33     = zero
    non_dim%scc      = zero
    non_dim%scc33    = zero
    non_dim%sco3     = zero
    non_dim%sco333   = zero
    non_dim%grasshof = zero
    non_dim%lfddv    = zero
    non_dim%lfddh    = zero

  END SUBROUTINE zero_initial_non_dim

  SUBROUTINE zero_initial_nitrogen()

    USE constants, ONLY: zero

    IMPLICIT NONE

    nitrogen%Ja   = zero
    nitrogen%J_glu   = zero
    nitrogen%J_Busch   = zero
    nitrogen%cn_bulk   = zero
    nitrogen%n_supply  = zero
    nitrogen%n_mult    = zero
    nitrogen%nitrate_per   = zero
    nitrogen%nitrite_per   = zero
    nitrogen%ammonia_per   = zero
    nitrogen%Busch_mol     = zero
    nitrogen%nitrate_mol   = zero
    nitrogen%nitrite_mol   = zero
    nitrogen%ammonia_mol   = zero

  END SUBROUTINE zero_initial_nitrogen

END MODULE types


