MODULE soils

  ! This module contains the soil water and energy routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4, i8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: litter_h2o
  PUBLIC :: litter_rain
  PUBLIC :: litter_resistance
  PUBLIC :: litter_rh
  PUBLIC :: set_litter_moisture      ! set humidity per soil layer
  PUBLIC :: set_litter_temp
  PUBLIC :: set_litter_texture
  PUBLIC :: set_soil_clapp           ! set hydraulic soil parameters after Clapp & Hornberger (1978)
  PUBLIC :: set_soil_litter_capacity ! compute soil heat capacity and conductivity
  PUBLIC :: set_soil_moisture        ! set humidity per soil layer
  PUBLIC :: set_soil_root            ! set rooting profile
  PUBLIC :: set_soil_saxton          ! set hydraulic soil parameters after Saxton & Rawls (2006)
  PUBLIC :: set_soil_temp            ! initialize deep soil temperature
  PUBLIC :: set_soil_texture         ! set soil texture per soil layer
  PUBLIC :: set_soil_time            ! initialize soil heat conduction equation time step
  PUBLIC :: soil_energy_balance      ! soil energy balance
  PUBLIC :: soil_h2o                 ! calculate soil moisture profil
  PUBLIC :: soil_resistance
  PUBLIC :: soil_rh
  PUBLIC :: throughfall              ! calculates throughfall and interception for all layers
  PUBLIC :: tridia                   ! solves tridagonal matrix

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE litter_h2o()

    USE types,         ONLY: soil, flux, time, wiso, prof
    USE constants,     ONLY: zero, one, e3, isotolerance
    USE setup,         ONLY: nwiso
    USE isotope_utils, ONLY: isorat
#ifdef DEBUG
    USE string_utils,  ONLY: num2str
    USE messages,      ONLY: message
#endif

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nwiso) :: rsoil
    REAL(wp), DIMENSION(1:nwiso) :: dw_l
    REAL(wp), DIMENSION(1:nwiso) :: rqinfl
    REAL(wp), DIMENSION(1:nwiso) :: tmp1
    REAL(wp) :: tmp2

    if (soil%z_litter < 1e-6_wp) then
       tmp2 = one / time%time_step
       flux%soilinfl(1:nwiso) = prof%throughfall(1,1:nwiso) * tmp2 ! convert from mm to mm/s
       ! if z_litter = 0 then theta_l is something like a rain buffer or a puddle or similar (mm)
       tmp1(1:nwiso)       = flux%soilinfl(1:nwiso) + soil%theta_l(1:nwiso) * tmp2
       rqinfl(1:nwiso)     = isorat(tmp1(1:nwiso), tmp1(1), wiso%lost(1:nwiso), wiso%lost(1))
       soil%qinfl(1)       = min(tmp1(1), soil%qinfl_max) ! (kg/m2/s = mm/s)
       soil%qinfl(1:nwiso) = soil%qinfl(1) * rqinfl(1:nwiso)
#ifdef DEBUG
       if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('LITTER_H2O: ', 'Lost 09 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       soil%theta_l(1:nwiso) = soil%theta_l(1:nwiso) + flux%soilinfl(1:nwiso)*time%time_step - &
            soil%qinfl(1:nwiso)*time%time_step
#ifdef DEBUG
       ! Check: if this happens, we have to code some thresholds here
       if (any(soil%theta_l(1:nwiso) < (-isotolerance))) then
          call message('LITTER_H2O: ', 'L1<0 @ ', num2str(time%daytime))
       end if
#endif
       where (soil%theta_l(1:nwiso) < isotolerance)
          soil%qinfl(1:nwiso)   = soil%qinfl(1:nwiso) + soil%theta_l(1:nwiso)
          soil%theta_l(1:nwiso) = zero
       end where
    else
       !  Litter water drainage into soil, from Ogee & Brunet (2002)
       soil%qinfl(1)       = min(max(zero, 1000._wp*soil%z_litter*(soil%theta_l(1)-soil%theta_l33)/ &
            time%time_step), soil%qinfl_max) !  kg/m2/s = mm/s
       rsoil(1:nwiso)      = isorat(soil%theta_l(1:nwiso), soil%theta_l(1), wiso%lost(1:nwiso), wiso%lost(1))
       soil%qinfl(1:nwiso) = rsoil * soil%qinfl(1)
#ifdef DEBUG
       if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('LITTER_H2O: ', 'Lost 10 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       ! new litter water status, mass balance
       dw_l(1:nwiso)  = (zero - flux%litterevap(1:nwiso)  - soil%qinfl(1:nwiso)) * (e3 / soil%z_litter)
#ifdef DEBUG
       if (any((soil%theta_l(1:nwiso)+dw_l(1:nwiso)*time%time_step) < zero)) then
          call message('LITTER_H2O: ', 'L2<0 @  ', num2str(time%daytime))
       end if
#endif
       soil%theta_l(1:nwiso) = soil%theta_l(1:nwiso) + dw_l(1:nwiso)*time%time_step
#ifdef DEBUG
       ! Check: if this happens, we have to code some thresholds here
       if (any(soil%theta_l(1:nwiso) < (-isotolerance))) then
          call message('LITTER_H2O: ', 'L3<0 @ ', num2str(time%daytime))
       end if
#endif
    end if

  END SUBROUTINE litter_h2o


  ! ------------------------------------------------------------------
  SUBROUTINE litter_rain()

    USE types,     ONLY: prof, time, flux, soil
    USE constants, ONLY: one, e3
    USE setup,     ONLY: nwiso

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nwiso) :: dw_l
    REAL(wp) :: ztmp

    ztmp = one / time%time_step
    flux%soilinfl(1:nwiso) = prof%throughfall(1,1:nwiso) * ztmp  ! convert from mm to mm/s
    ! new litter water status, mass balance
    dw_l(1:nwiso) = flux%soilinfl(1:nwiso) * (e3 / soil%z_litter) ! mm/s
    soil%theta_l(1:nwiso) = soil%theta_l(1:nwiso) + dw_l(1:nwiso)*time%time_step ! mm

  END SUBROUTINE litter_rain


  ! ------------------------------------------------------------------
  FUNCTION litter_resistance()

    USE types,      ONLY: soil
    USE constants,  ONLY: zero, twothird, one, TN0
    USE parameters, ONLY: extra_nate

    IMPLICIT NONE

    REAL(wp) :: litter_resistance

    REAL(wp), PARAMETER :: dw_0 = 0.229e-4_wp ! vapour diffusion in air @ TN0
    REAL(wp) :: dw ! vapour diffusion in air @ TN0
    REAL(wp) :: rlmax

    dw = dw_0 * (soil%T_l_filter/TN0)**1.75_wp ! from Kondo et al. (1990)
    if (extra_nate == 1) then
       ! for Nate McDowell''s juniper site, fix maximum litter resistance
       rlmax = 5050._wp    !soil%z_litter/(dw* pow(soil%theta_ls, 2./3.)) ! after Schaap & Bouten (1997)
       rlmax = 4200._wp
    else
       rlmax = soil%z_litter/(dw*soil%theta_ls**twothird) ! after Schaap & Bouten (1997)
       ! with Milly model (2/3) after Saravanapavan et al. (2000).
    end if
    litter_resistance = zero
    if (soil%theta_l(1) < soil%theta_l33) then
       litter_resistance = rlmax * (one-soil%theta_l(1)/soil%theta_l33)**soil%n_l
    end if

  END FUNCTION litter_resistance


  ! ------------------------------------------------------------------
  FUNCTION litter_rh()

    USE types,      ONLY: soil
    USE constants,  ONLY: zero, one, Gravity, Rw

    IMPLICIT NONE

    REAL(wp) :: litter_rh

    REAL(wp), PARAMETER :: psi_l1 = 35.3_wp  ! Ogee & Brunet (2002)
    REAL(wp), PARAMETER :: bl     = 2.42_wp
    REAL(wp), PARAMETER :: rho_w  = 1000._wp ! density of water
    REAL(wp) :: psi_l

    if (abs(soil%theta_l(1)) > tiny(one)) then
       psi_l = psi_l1 * (rho_w*soil%theta_l(1)/soil%rho_l)**(-bl)
    else
       psi_l = zero
    end if
    litter_rh = exp(-Gravity*psi_l/(Rw*soil%T_l_filter))
    litter_rh = max(min(litter_rh, one), zero)

  END FUNCTION litter_rh


  ! ------------------------------------------------------------------
  SUBROUTINE set_litter_moisture()

    USE types,      ONLY: soil
    USE constants,  ONLY: zero
    USE setup,      ONLY: nwiso

    IMPLICIT NONE

    if (soil%z_litter < 1e-6_wp) then
       soil%theta_l(1:nwiso) = zero !m3 m-3
    else
       soil%theta_l(1:nwiso) = soil%theta(1,1:nwiso) !m3 m-3
    end if

  END SUBROUTINE set_litter_moisture


  ! ------------------------------------------------------------------
  SUBROUTINE set_litter_temp()
    ! Assign litter temperature
    USE types,      ONLY: soil
    USE constants,  ONLY: TN0

    IMPLICIT NONE

    soil%T_l        = soil%T_base+TN0
    soil%T_l_filter = soil%T_base+TN0

  END SUBROUTINE set_litter_temp


  ! ------------------------------------------------------------------
  SUBROUTINE set_litter_texture()
    ! Assign litter temperature
    USE types,      ONLY: soil
    USE constants,  ONLY: one

    IMPLICIT NONE

    ! dry density of std organic matter (kg m-3)
    ! average of table 2 of Redding et al. (Soil Science Society of America Journal, 2005)
    REAL(wp), PARAMETER :: rho_0 = 80._wp

    ! same as Saxton & Rawls (2006), i.e. theta_s = porosity
    soil%theta_ls  = one - soil%rho_l/rho_0
    ! fit of rmax*(1-theta/theta_field_capacity)^n to data of Schaap & Bouten (1997)
    soil%theta_l33 = 0.25_wp
    ! dito, see fit_litter_resistance.pro
    soil%n_l       = 2.5_wp

  END SUBROUTINE set_litter_texture


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_clapp()
    ! Soil characteristics after Clapp & Hornberger (1978)
    USE types,      ONLY: soil
    USE constants,  ONLY: one, e2
    USE parameters, ONLY: extra_nate
    USE setup,      ONLY: nsoil

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nsoil) :: sand, clay
    REAL(wp), DIMENSION(1:nsoil) :: alpha, gravelw ! Matric soil density/gravel density, Weight fraction of gravel
    REAL(wp) :: ztmp

    ! Sand, Clay, OM in % in parameter file
    sand(1:nsoil) = soil%sand(1:nsoil)
    clay(1:nsoil) = soil%clay(1:nsoil)
    ! Clapp & Hornberger
    soil%theta_s(1:nsoil)    = 0.489_wp - 0.00126_wp*soil%sand(1:nsoil)
    soil%k_s(1:nsoil)        = 0.0070556_wp * 10._wp**(-0.884_wp+0.0153_wp*soil%sand(1:nsoil))
    soil%psi_e(1:nsoil)      = -10._wp * 10._wp**(1.88_wp-0.0131_wp*soil%sand(1:nsoil))
    soil%big_b(1:nsoil)      = 2.91_wp + 0.159_wp*soil%clay(1:nsoil)
    soil%theta_1500(1:nsoil) = soil%theta_s(1:nsoil) * (-316230._wp/soil%psi_e(1:nsoil))**(-one/soil%big_b(1:nsoil))
    soil%theta_33(1:nsoil)   = soil%theta_s(1:nsoil) * (-158490._wp/soil%psi_e(1:nsoil))**(-one/soil%big_b(1:nsoil))
    soil%theta_s33(1:nsoil)  = soil%theta_s(1:nsoil) - soil%theta_33(1:nsoil)
    ! limit theta_s33 to 0.5 %v and recalc theta_33
    where (soil%theta_s33(1:nsoil) < 0.005_wp) soil%theta_s33(1:nsoil) = 0.005_wp
    soil%theta_33(1:nsoil) = soil%theta_s(1:nsoil) - soil%theta_s33(1:nsoil)
    ! Calc density for Nate McDowell''s Juniper site
    if (extra_nate == 1) then
       soil%rho(1:nsoil) = (one-soil%theta_s(1:nsoil)) * 2.65_wp
    else
       soil%rho(1:nsoil) = soil%bulk_density(1:nsoil)
    end if
    ! Correct for gravel
    soil%rho(1:nsoil) = soil%rho(1:nsoil)*(one-soil%gravel(1:nsoil)*e2) + soil%gravel(1:nsoil)*e2*2.65_wp
    ztmp = one / 2.65_wp
    alpha(1:nsoil)    = soil%rho(1:nsoil) * ztmp
    gravelw(1:nsoil)  = soil%gravel(1:nsoil)*e2 / (alpha(1:nsoil) + soil%gravel(1:nsoil)*e2*(one-alpha(1:nsoil)))
    soil%k_s(1:nsoil) = soil%k_s(1:nsoil) * (one-gravelw(1:nsoil)) / (one-gravelw(1:nsoil)*(one-1.5_wp*alpha(1:nsoil)))
    ! Calculate root weighted field capacity and wilting point
    ! soil%soil_mm_33_root is used for soil water stress
    soil%soil_mm_33_root   = sum(soil%root(1:nsoil)*soil%theta_33(1:nsoil)*1000._wp * &
         (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))
    ! soil%soil_mm_1500_root is used for soil water stress
    soil%soil_mm_1500_root = sum(soil%root(1:nsoil)*soil%theta_1500(1:nsoil)*1000._wp * &
         (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))

  END SUBROUTINE set_soil_clapp


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_litter_capacity()
    ! Compute layer heat capacities and conductivities
    ! as a function of texture, bulk density and soil moisture
    ! Reference: Campbell, Soil physics with Basic (1985, Chapter 4)
    USE types,      ONLY: soil
    USE constants,  ONLY: one, two, e2
    USE setup,      ONLY: nsoil

    IMPLICIT NONE

    REAL(wp), PARAMETER :: rho_0 = 1300._wp ! wet density of std organic matter (kg m-3)
    REAL(wp), PARAMETER :: c_0   = 2511600._wp ! heat capacity of std organic matter (J m-3 K-1)
    REAL(wp), PARAMETER :: c_w   = 4185000._wp ! heat capacity of water (J m-3 K-1)
    REAL(wp), PARAMETER :: c_s   = 2000000._wp ! heat capacity of soil matter (J m-3 K-1)
    ! soil
    REAL(wp), DIMENSION(0:nsoil) :: C1, C2, C3, C4, C5
    REAL(wp), DIMENSION(1:nsoil) :: rho_local
    REAL(wp) :: ztmp1

    ! Soil
    ! as in Campbell (1985)
    rho_local(1:nsoil-1)    = (one-soil%gravel(1:nsoil-1)*e2)*soil%rho(1:nsoil-1) + soil%gravel(1:nsoil-1)*(e2*2.65_wp)
    soil%cp_soil(1:nsoil-1) = (((c_s/2.65_wp)*rho_local(1:nsoil-1)) + &
         c_w*soil%theta(1:nsoil-1,1)*(one-soil%gravel(1:nsoil-1)*e2)) * &
         (soil%z_soil(2:nsoil) - soil%z_soil(0:nsoil-2)) / (two*soil%temperature_dt)
    !old=0.6 instead of 0.9
    C1(1:nsoil-1) = 0.65_wp - 0.78_wp * rho_local(1:nsoil-1) + 0.6_wp * rho_local(1:nsoil-1)*rho_local(1:nsoil-1)
    C2(1:nsoil-1) = 1.06_wp * rho_local(1:nsoil-1) !old =1.06
    C3(1:nsoil-1) = one + 2.6_wp / sqrt(soil%clay(1:nsoil-1))
    C4(1:nsoil-1) = 0.03_wp + 0.1_wp * rho_local(1:nsoil-1)*rho_local(1:nsoil-1)
    C5(1:nsoil-1) = 4._wp
    soil%k_conductivity_soil(1:nsoil-1) = (C1(1:nsoil-1) + &
         C2(1:nsoil-1) * soil%theta(1:nsoil-1,1)*(one-soil%gravel(1:nsoil-1)*e2) - &
         (C1(1:nsoil-1) - C4(1:nsoil-1)) * &
         exp(-(C3(1:nsoil-1)*soil%theta(1:nsoil-1,1)*(one-soil%gravel(1:nsoil-1)*e2))**C5(1:nsoil-1))) / &
         (soil%z_soil(2:nsoil) - soil%z_soil(1:nsoil-1))
    ! ! ISOLSM formulation
    ! C1 = (8.80*soil%sand(i)+2.92*soil%clay(i))/(soil%sand(i)+soil%clay(i))
    ! C2 = 0.15
    ! C3 = (pow(C1,1.-soil%theta_s(i)) * pow(0.6,soil%theta(i,1)) - C2 )
    !   * soil%theta(i,1)/soil%theta_s(i) + C2
    ! soil%k_conductivity_soil(i) = C3 / (soil%z_soil(i+1) - soil%z_soil(i))
    ! Last layer is different (not in Campbell (1985) who has an extra layer underneath the last soil layer)
    ! Take thickness of last soil layer instead of mean thickness of two layers
    ztmp1 = one / 2.65_wp
    rho_local(nsoil)    = (one-soil%gravel(nsoil)*e2)*soil%rho(nsoil) + soil%gravel(nsoil)*e2*2.65_wp
    soil%cp_soil(nsoil) =  ((c_s*rho_local(nsoil)*ztmp1) + c_w*soil%theta(nsoil,1)*(one-soil%gravel(nsoil)*e2)) * &
         (soil%z_soil(nsoil) - soil%z_soil(nsoil-1)) / soil%temperature_dt
    C1(nsoil) = 0.65_wp - 0.78_wp * rho_local(nsoil) + 0.6_wp * rho_local(nsoil) * rho_local(nsoil)
    C2(nsoil) = 1.06_wp * rho_local(nsoil)
    C3(nsoil) = one + 2.6_wp / sqrt(soil%clay(nsoil))
    C4(nsoil) = 0.03_wp + 0.1_wp * rho_local(nsoil) * rho_local(nsoil)
    C5(nsoil) = 4._wp
    soil%k_conductivity_soil(nsoil) = (C1(nsoil) + &
         C2(nsoil) * soil%theta(nsoil,1)*(one-soil%gravel(nsoil)*e2) - &
         (C1(nsoil) - C4(nsoil)) * exp(-(C3(nsoil)*soil%theta(nsoil,1)*(one-soil%gravel(nsoil)*e2))**C5(nsoil))) / &
         (soil%z_soil(nsoil) - soil%z_soil(nsoil-1))
    ! ! ISOLSM formulation
    ! C1 = (8.80*soil%sand(i)+2.92*soil%clay(i))/(soil%sand(i)+soil%clay(i))
    ! C2 = 0.15
    ! C3 = (pow(C1,1.-soil%theta_s(i)) * pow(0.6,soil%theta(i,1)) - C2 )
    !   * soil%theta(i,1)/soil%theta_s(i) + C2
    ! soil%k_conductivity_soil(i) = C3  / (soil%z_soil(i) - soil%z_soil(i-1))
    if (soil%z_litter >= 1e-6_wp) then
       ! Layer 0 = litter
       !soil%cp_soil(0) = (c_0/rho_0*soil%rho_l + c_w*soil%theta_l(1)) * soil%z_litter / soil%temperature_dt
       soil%cp_soil(0) = (c_0*soil%rho_l/rho_0 + c_w*soil%theta_l(1)) * soil%z_litter / (two*soil%temperature_dt)
       ! Campbell (1985, Chapter 4)
       C1(0) = 0.4_wp
       C2(0) = 0.5_wp
       C3(0) = 1._wp
       C4(0) = 0.06_wp
       C5(0) = 4._wp
       soil%k_conductivity_soil(0) = (C1(0) + C2(0) * soil%theta_l(1) - &
            (C1(0)-C4(0)) * exp(-(C3(0)*soil%theta_l(1))**C5(0))) / soil%z_litter
    end if

  END SUBROUTINE set_soil_litter_capacity


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_moisture()
    USE types,         ONLY: soil, wiso
    USE setup,         ONLY: nsoil, nwiso
    USE isotope_utils, ONLY: invdelta1000_h2o

    IMPLICIT NONE

    INTEGER(i4) :: i
    INTEGER(i4), DIMENSION(nwiso-1) :: mc

    if (nwiso > 1) then
       forall(i=1:nwiso-1) mc(i) = i+1
       ! convert delta values read from parameter file to m3 m-3
       soil%theta(1:nsoil,2:nwiso) = spread(soil%theta(1:nsoil,1),dim=2,ncopies=nwiso-1) * &
            invdelta1000_h2o(wiso%dtheta(1:nsoil,2:nwiso), mc(1:nwiso-1)) !m3 m-3
    end if

  END SUBROUTINE set_soil_moisture


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_root()
    USE types,      ONLY: soil
    USE constants,  ONLY: one
    USE setup,      ONLY: nsoil
    USE parameters, ONLY: extra_nate

    IMPLICIT NONE

    REAL(wp) :: root_total ! total has to sum up to 1

    soil%root(1:nsoil) = soil%root_factor**(100._wp*soil%z_soil(1:nsoil)) - &
         soil%root_factor**(100._wp*soil%z_soil(0:nsoil-1))
    ! for Nate McDowell''s juniper site give root distribution
    if (extra_nate == 1) then
       soil%root(1)  = 0.019927536_wp
       soil%root(2)  = 0.027173913_wp
       soil%root(3)  = 0.038043478_wp
       soil%root(4)  = 0.052536232_wp
       soil%root(5)  = 0.066847826_wp
       soil%root(6)  = 0.084057971_wp
       soil%root(7)  = 0.102717391_wp
       soil%root(8)  = 0.122826087_wp
       soil%root(9)  = 0.143115942_wp
       soil%root(10) = 0.342753623_wp
    end if
    ! normalise to total=1 ! for comparison with v3.3.8 move above fixed root distribution
    root_total = one / sum(soil%root(1:nsoil))
    soil%root(1:nsoil) = soil%root(1:nsoil) * root_total

  END SUBROUTINE set_soil_root


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_saxton()
    ! Soil characteristics after Saxton & Rawls (2006)
    USE types,      ONLY: soil
    USE constants,  ONLY: one, e2
#ifdef DEBUG
    USE constants,  ONLY: nerr
#endif
    USE setup,      ONLY: nsoil
    USE parameters, ONLY: extra_nate
#ifdef DEBUG
    USE messages,   ONLY: message
#endif

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nsoil) :: sand, clay, om
    REAL(wp), DIMENSION(1:nsoil) :: theta_1500t
    REAL(wp), DIMENSION(1:nsoil) :: theta_33t
    REAL(wp), DIMENSION(1:nsoil) :: theta_s33t
    REAL(wp), DIMENSION(1:nsoil) :: alpha, gravelw ! Matric soil density/gravel density, Weight fraction of gravel
    REAL(wp), DIMENSION(1:nsoil) :: psi_et
    REAL(wp), DIMENSION(1:nsoil) :: df ! used in density correction for compacted soil
    REAL(wp), DIMENSION(1:nsoil) :: theta_temp  ! used in density correction for compacted soil
    REAL(wp) :: z265

    ! Sand, Clay, OM in % in parameter file
    sand(1:nsoil) = soil%sand(1:nsoil)*e2
    clay(1:nsoil) = soil%clay(1:nsoil)*e2
    om(1:nsoil)   = soil%om(1:nsoil)*e2
    ! Set normal density wilting point (1500 kPa), field capacity (33 kPa), excess (Sat-field)
    ! air entry tension (potential), saturation moisture, density
    theta_1500t(1:nsoil) = -0.024_wp*sand(1:nsoil) + 0.487_wp*clay(1:nsoil) + 0.006_wp*om(1:nsoil) + &
         0.005_wp*sand(1:nsoil)*om(1:nsoil) - 0.013_wp*clay(1:nsoil)*om(1:nsoil) + 0.068_wp*sand(1:nsoil)*clay(1:nsoil) + &
         0.031_wp
    soil%theta_1500(1:nsoil) = theta_1500t(1:nsoil) + (0.14_wp*theta_1500t - 0.02_wp)
    theta_33t(1:nsoil) = -0.251_wp*sand(1:nsoil) + 0.195_wp*clay(1:nsoil) + 0.011_wp*om(1:nsoil) + &
         0.006_wp*sand(1:nsoil)*om(1:nsoil) - 0.027_wp*clay(1:nsoil)*om(1:nsoil) + 0.452_wp*sand(1:nsoil)*clay(1:nsoil) + &
         0.299_wp
    soil%theta_33(1:nsoil) = theta_33t(1:nsoil) + (1.283_wp*theta_33t*theta_33t - 0.374_wp*theta_33t - 0.015_wp)
    theta_s33t(1:nsoil) = 0.278_wp*sand(1:nsoil) + 0.034_wp*clay(1:nsoil) + 0.022_wp*om(1:nsoil) - &
         0.018_wp*sand(1:nsoil)*om(1:nsoil) - 0.027_wp*clay(1:nsoil)*om(1:nsoil) - 0.584_wp*sand(1:nsoil)*clay(1:nsoil) + &
         0.078_wp
    soil%theta_s33(1:nsoil) = theta_s33t(1:nsoil) + (0.636_wp*theta_s33t - 0.107_wp)
    psi_et(1:nsoil) = -21.674_wp*sand(1:nsoil) - 27.932_wp*clay(1:nsoil) - 81.975_wp*soil%theta_s33(1:nsoil) + &
         71.121_wp*sand(1:nsoil)*soil%theta_s33(1:nsoil) + 8.294_wp*clay(1:nsoil)*soil%theta_s33(1:nsoil) + &
         14.05_wp*sand(1:nsoil)*clay(1:nsoil) + 27.161_wp
    soil%psi_e(1:nsoil) = psi_et(1:nsoil) + (0.02_wp*psi_et*psi_et - 0.113_wp*psi_et - 0.70_wp)
    soil%theta_s(1:nsoil) = soil%theta_33(1:nsoil) + soil%theta_s33(1:nsoil) - 0.097_wp*sand(1:nsoil) + 0.043_wp
    soil%rho(1:nsoil) = (one-soil%theta_s(1:nsoil)) * 2.65_wp
    ! exclude for Nate McDowell''s Juniper site
    z265 = one / 2.65_wp
    if (extra_nate == 0) then
       ! Correct for compacted soil
       df(1:nsoil) = soil%bulk_density(1:nsoil) / soil%rho(1:nsoil)
#ifdef DEBUG
       if (any(df < 0.9_wp) .or. any(df > 1.3_wp)) then
          call message('SET_SOIL_SAXTON: ','Warning: density factor not in range (0.9,1.3) of Saxton & Rawls (2006) ')
          write(nerr,*) 'df:            ', df
          write(nerr,*) 'Given density: ', soil%bulk_density(1:nsoil)
          write(nerr,*) 'Calc density:  ', soil%rho(1:nsoil)
       end if
#endif
       soil%rho(1:nsoil)       = soil%bulk_density(1:nsoil)
       theta_temp(1:nsoil)     = one - soil%rho(1:nsoil) * z265
       soil%theta_33(1:nsoil)  = soil%theta_33(1:nsoil) - 0.2_wp*(soil%theta_s(1:nsoil) - theta_temp(1:nsoil))
       soil%theta_s(1:nsoil)   = theta_temp(1:nsoil)
       soil%theta_s33(1:nsoil) = soil%theta_s(1:nsoil) - soil%theta_33(1:nsoil)
       ! limit theta_s33 to 0.5 %v and recalc theta_33
       where (soil%theta_s33(1:nsoil) < 0.005_wp) soil%theta_s33(1:nsoil) = 0.005_wp
       soil%theta_33(1:nsoil)  = soil%theta_s(1:nsoil) - soil%theta_s33(1:nsoil)
    end if
    ! Moisture parameters
    soil%big_b(1:nsoil)  = (log(1500._wp)-log(33._wp)) / (log(soil%theta_33(1:nsoil))-log(soil%theta_1500(1:nsoil)))
    soil%big_a(1:nsoil)  = exp(log(33._wp)+soil%big_b(1:nsoil)*log(soil%theta_33(1:nsoil)))
    soil%lambda(1:nsoil) = one / soil%big_b(1:nsoil)
    soil%k_s(1:nsoil)    = 1930._wp * (soil%theta_s(1:nsoil)-soil%theta_33(1:nsoil))**(3._wp-soil%lambda(1:nsoil))
    ! Correct for gravel
    soil%rho(1:nsoil) = soil%rho(1:nsoil)*(one-soil%gravel(1:nsoil)*e2) + soil%gravel(1:nsoil)*e2*2.65_wp
    alpha(1:nsoil)    = soil%rho(1:nsoil) * z265
    gravelw(1:nsoil)  = soil%gravel(1:nsoil)*e2 / (alpha(1:nsoil) + soil%gravel(1:nsoil)*e2*(one-alpha(1:nsoil)))
    soil%k_s(1:nsoil) = soil%k_s(1:nsoil) * (one-gravelw(1:nsoil)) / (one-gravelw(1:nsoil)*(one-1.5_wp*alpha(1:nsoil)))
    ! Calculate root weighted field capacity and wilting point
    ! soil%soil_mm_33_root is used for soil water stress
    soil%soil_mm_33_root   = sum(soil%root(1:nsoil)*soil%theta_33(1:nsoil)*1000._wp * &
         (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))
    ! soil%soil_mm_1500_root is used for soil water stress
    soil%soil_mm_1500_root = sum(soil%root(1:nsoil)*soil%theta_1500(1:nsoil)*1000._wp * &
         (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))

  END SUBROUTINE set_soil_saxton


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_temp()
    ! Assign soil temperature
    USE types, ONLY: soil
    USE setup, ONLY: nsoil

    IMPLICIT NONE

    soil%T_soil(0:nsoil+1)        = soil%T_base
    soil%T_soil_filter(0:nsoil+1) = soil%T_base

  END SUBROUTINE set_soil_temp


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_texture()
    USE types,      ONLY: soil
    USE constants,  ONLY: one, e2, s2h
    USE parameters, ONLY: extra_nate
    USE setup,      ONLY: nsoil

    IMPLICIT NONE

    soil%clay(1:nsoil) = soil%clay_in(1:nsoil)*soil%clay_factor !clay fraction
    soil%sand(1:nsoil) = soil%sand_in(1:nsoil)*soil%sand_factor !sand fraction
    ! watmin = 0.001  !(m3 m-3), based on Konukcu et al. (AJSR, 2004)
    if (extra_nate == 1) then
       ! for Nate McDowell''s juniper site
       soil%watmin(1:nsoil) = 0.1_wp ! from measurements
    else
       ! estimate residual volumetric water content based on data from EPA report and selfmade regression
       soil%watmin(1:nsoil) = 0.211_wp*soil%clay(1:nsoil)*e2 + 0.009_wp*soil%sand(1:nsoil)*e2 + 0.028_wp  !(m3 m-3)
    end if
    ! soil%qinfl_max = 5./3600. ! (mm) 5 mm h-1 is typical values from literature for loamy-clay soil
    ! set maximal infiltration rate depending on soil texture at soil surface
    !   (based on FAO report, selfmade regression)
    soil%qinfl_max = max(-26.94_wp*soil%clay(1) + 19.27_wp*soil%sand(1) + 13.26_wp, one) * s2h

  END SUBROUTINE set_soil_texture


  ! ------------------------------------------------------------------
  SUBROUTINE set_soil_time()
    USE types, ONLY: soil, time

    IMPLICIT NONE

    soil%temperature_dt = 15._wp ! Time step in seconds
    soil%temperature_mtime = nint(time%time_step/soil%temperature_dt, kind=i4) ! time steps per hour

    soil%moisture_dt = 300._wp ! Time step in seconds
    !soil%moisture_dt = 600._wp ! Time step in seconds
    soil%moisture_mtime = nint(time%time_step/soil%moisture_dt, kind=i4) ! time steps per hour

  END SUBROUTINE set_soil_time


  ! ------------------------------------------------------------------
  SUBROUTINE soil_energy_balance()
    ! The soil energy balance model of Campbell has been adapted to
    ! compute soil energy fluxes and temperature profiles at the soil
    ! surface. The model has been converted from BASIC to C. We
    ! also use an analytical version of the soil surface energy
    ! balance to solve for LE, H and G.

    ! Campbell GS (1985) Soil physics with Basic: Transport models for soil-plant systems.
    ! Elsevier, New York, p 150

    ! Combine surface energy balance calculations with soil heat
    ! transfer model to calculate soil conductive and convective heat
    ! transfer and evaporation rates. Here, only the deep temperature
    ! is needed and G, Hs and LEs can be derived from air temperature
    ! and energy inputs.
    ! Soil evaporation models by Kondo, Mafouf et al. and
    ! Dammond and Simmonds are used. Dammond and Simmonds for example
    ! have a convective adjustment to the resistance to heat transfer.
    ! Our research in Oregon and Canada have shown that this consideration
    ! is extremely important to compute G and Rn_soil correctly.
    USE types,        ONLY: soil, flux, output, solar, prof, met, fact, time, input, iswitch
    USE constants,    ONLY: zero, half, one, two, &
         TN0, sigma, Gravity, vonKarman, Rw, cp
    USE setup,        ONLY: nsoil
    USE parameters,   ONLY: epsoil
    USE transport,    ONLY: uz
    USE utils,        ONLY: es, desdt, des2dt, lambda
#ifdef DEBUG
    USE messages,     ONLY: message
    USE string_utils, ONLY: num2str
#endif

    IMPLICIT NONE

    INTEGER(i4) :: j, mm1, i, ip1
    INTEGER(i4) :: n_soil_end
    REAL(wp) :: Fst, Gst
    REAL(wp) :: soil_par, soil_nir
    REAL(wp) :: u_soil, Rh_soil, kv_soil, kcsoil
    REAL(wp), DIMENSION(0:nsoil+2) :: a_soil, b_soil, c_soil, d_soil
    REAL(wp) :: est, dest, d2est, tk2, tk3, tk4, llout, lecoef
    REAL(wp) :: product
    REAL(wp) :: att, btt, ctt
    REAL(wp) :: storage
    REAL(wp) :: facstab, stabdel, ea
    REAL(wp), DIMENSION(0:nsoil+2) :: T_new_soil
    REAL(wp), DIMENSION(0:nsoil+2) :: k_soil, cp_soil, T_soil
    REAL(wp) :: kv_litter, temp
    REAL(wp) :: firstlayerevap, secondlayerevap

    ! Compute soilevap as a function of energy balance at the soil
    ! surface. Net incoming short and longwave energy

    ! radiation balance at soil in PAR band, W m-2
    !soil_par = (solar%beam_flux_par(1) + solar%par_down(1) - solar%par_up(1)) / 4.6
    soil_par   = (solar%beam_flux_par(1) + solar%par_down(1)) / 4.6_wp * (one-solar%par_soil_refl)
    ! radiation balance at soil in NIR band, W m-2
    !soil_nir = (solar%beam_flux_nir(1) + solar%nir_dn(1) - solar%nir_up(1))
    soil_nir   = (solar%beam_flux_nir(1) + solar%nir_dn(1)) * (one-solar%nir_soil_refl)
    ! incoming radiation balance at soil, solar and terrestrial, W m-2
    soil%rnet  = soil_par + soil_nir + solar%ir_dn(1)*epsoil
    ! set air temperature over soil with lowest air layer, filtered
    soil%T_air = prof%tair_filter(1)
    ! Compute Rh_soil and rv_soil from wind log profile for lowest layer
    prof%u(1)  = uz(prof%ht(1)*half)
    u_soil     = prof%u(1) ! wind speed one layer above soil
    ! Stability factor from Daamen and Simmonds (1996)
    stabdel = 5._wp*Gravity*(prof%ht(1)*half)*(soil%tsrf_filter-soil%T_air)/((soil%T_air+TN0)*u_soil*u_soil)
    if (stabdel > zero) then
       facstab = one/((one+stabdel)**0.75_wp)
    else
       facstab = one/((one+stabdel)*(one+stabdel))
    end if
    if (time%count <= 1) facstab = one
    facstab = min(max(facstab, 0.1_wp), 5._wp)
    Rh_soil = log((half*prof%ht(1)-soil%d)/soil%z0)*log((half*prof%ht(1)-soil%d)/soil%z0) / &
         (vonKarman*vonKarman*u_soil) * facstab
    Rh_soil = min(max(Rh_soil, 50._wp), 3000._wp)
    ! kcsoil is the convective transfer coeff for the soil% (W m-2 K-1)
    kcsoil         = (cp * met%air_density) / Rh_soil
    ! soil boundary layer resistance
    soil%rb        = Rh_soil
    ! soil resistance
    soil%rs        = soil_resistance()
    ! litter resistance
    soil%rl        = litter_resistance()
    ! kvsoil is soil surface conductance to water vapor transfer (s m-1)
    kv_soil        = one / (soil%rb+soil%rs)
    ! kvsoil is litter surface conductance to water vapor transfer (s m-1)
    kv_litter      = one / (soil%rb+soil%rl)
    ! soil relative humidty
    soil%rh_soil   = soil_rh()
    ! litter relative humidty
    soil%rh_litter = litter_rh()
    ! initialise soil temperatures
    ! Note: the first internal soil layer becomes the litter layer
    ! all subsequent soil layers are shifted one layer down
    if (soil%z_litter < 1e-6_wp) then
       n_soil_end = nsoil
       ! does not include litter
       T_soil(1:nsoil)  = soil%T_soil_filter(1:nsoil)
       k_soil(1:nsoil)  = soil%k_conductivity_soil(1:nsoil)
       cp_soil(1:nsoil) = soil%cp_soil(1:nsoil)
    else
       n_soil_end = nsoil+1
       ! includes litter (0) at soil%T_soil, (1) at T_soil
       T_soil(1:nsoil+1) = soil%T_soil_filter(0:nsoil)
       k_soil(1:nsoil+1) = soil%k_conductivity_soil(0:nsoil)
       cp_soil(1:nsoil+1) = soil%cp_soil(0:nsoil)
       T_soil(1) = soil%T_l_filter-TN0 ! redundant
    end if
    T_soil(0)                = soil%T_air
    T_new_soil(0)            = soil%T_air
    T_soil(n_soil_end+1)     = soil%T_base
    T_new_soil(n_soil_end+1) = soil%T_base
    ! Boundary layer conductance at soil surface, W m-2 K-1
    k_soil(0)     = kcsoil
    ! initialize absolute soil surface temperature
    soil%T_Kelvin = soil%T_air+TN0
    ! evaluate latent heat of evaporation at T_soil=T_air
    soil%latent   = lambda(soil%T_Kelvin)
    !soil%latent = LAMBDA(soil%tsrf_filter+TN0)
    ! evaluate saturation vapor pressure at T_air
    est          = es(soil%T_Kelvin) * 100._wp ! 100 converts es(T) from mb to Pa
    ! Slope of the vapor pressure-temperature curve, Pa/C, evaluate as function of Tk=T_air
    fact%latent  = soil%latent
    dest         = desdt(soil%T_Kelvin, fact%latent)
    ! Second derivative of the vapor pressure-temperature curve, Pa/C, evaluate as function of Tk=T_air
    d2est        = des2dt(soil%T_Kelvin)
    ! Atmospheric vapour pressure (Pa)
    ea           = prof%rhov_air_filter(1,1) * soil%T_Kelvin * Rw
    ! output
    output%c8    = est - ea ! Vapor pressure deficit, Pa
    ! Compute products of absolute air temperature
    tk2          = soil%T_Kelvin * soil%T_Kelvin
    tk3          = tk2 * soil%T_Kelvin
    tk4          = tk3 * soil%T_Kelvin
    ! Longwave emission at air temperature, W m-2
    llout        = epsoil * sigma * tk4
    ! IR emissive flux density from soil (W m-2)
    soil%lout    = llout
    ! coefficients for latent heat flux density
    lecoef       = met%air_density * 0.622_wp * soil%latent / met%press_Pa
    ! Weighting factors for solving diff eq.
    Fst          = 0.6_wp
    Gst          = one - Fst
    ! solve by looping through the d()/dt term of the Fourier
    ! heat transfer equation
    soil%soilevap   = soil%soilevap_filter
    soil%litterevap = soil%litterevap_filter
    do j=1, soil%temperature_mtime
       if (soil%z_litter < 1e-6_wp) then
          firstlayerevap  = soil%soilevap
          secondlayerevap = zero
       else
          firstlayerevap  = soil%litterevap
          secondlayerevap = soil%soilevap
       end if
       ! define coef for each soil layer
       c_soil(1:n_soil_end)   = -k_soil(1:n_soil_end)*Fst
       a_soil(2:n_soil_end+1) = c_soil(1:n_soil_end)
       b_soil(1:n_soil_end)   = Fst*(k_soil(1:n_soil_end)+k_soil(0:n_soil_end-1)) + cp_soil(1:n_soil_end)
       d_soil(1:n_soil_end)   = Gst*k_soil(0:n_soil_end-1)*T_soil(0:n_soil_end-1) + &
            (cp_soil(1:n_soil_end)-Gst*(k_soil(1:n_soil_end) + k_soil(0:n_soil_end-1)))*T_soil(1:n_soil_end) + &
            Gst*k_soil(1:n_soil_end)*T_soil(2:n_soil_end+1)
       d_soil(1)              = d_soil(1) + k_soil(0)*T_new_soil(0)*Fst + soil%rnet - soil%lout - firstlayerevap
       d_soil(2)              = d_soil(2) - secondlayerevap
       d_soil(n_soil_end)     = d_soil(n_soil_end) + k_soil(n_soil_end)*Fst*T_new_soil(n_soil_end+1)
       mm1 = n_soil_end-1
       do i=1, mm1
          ip1 = i+1
          c_soil(i)   = c_soil(i)/b_soil(i)
          d_soil(i)   = d_soil(i)/b_soil(i)
          b_soil(ip1) = b_soil(ip1) - a_soil(ip1)*c_soil(i)
          d_soil(ip1) = d_soil(ip1) - a_soil(ip1)*d_soil(i)
       end do
       T_new_soil(n_soil_end) = d_soil(n_soil_end)/b_soil(n_soil_end)
       do i=mm1, 1, -1
          T_new_soil(i) = d_soil(i) - c_soil(i)*T_new_soil(i+1)
       end do
       ! soil temperature at 15 cm
       soil%T_15cm = T_new_soil(7)
       ! compute soil conductive heat flux density, W m-2
       if (soil%z_litter < 1e-6_wp) then
          soil%gsoil = k_soil(1)*(T_new_soil(1)-T_new_soil(2)) ! Matthias ??? why 1 and 2
          ! Does not work without litter, do not know why
          !storage = 0.
          !for (k=1 k<=n_soil_end k++)
          ! storage += cp_soil(k)*(T_new_soil(k)-T_soil(k))
          !soil%gsoil += storage
       else
          ! take 2 and 3 because 1 is litter
          storage    = sum(cp_soil(2:n_soil_end)*(T_new_soil(2:n_soil_end)-T_soil(2:n_soil_end)))
          soil%gsoil = k_soil(2)*(T_new_soil(2)-T_new_soil(3)) + storage
       end if
       !output soil heat flux (W m-1)
       output%c7 = soil%gsoil
       ! solve for T_litter using quadratic solution x1/2 = (-b +- sqrt(b^2-4ac) )/2a
       att = 6._wp*epsoil*sigma*tk2 + half*d2est*lecoef*(kv_soil*soil%rh_soil+soil%c_litterevap*kv_litter*soil%rh_litter)
       btt = 4._wp*epsoil*sigma*tk3 + kcsoil + &
            lecoef*(dest*(kv_soil*soil%rh_soil+soil%c_litterevap*kv_litter*soil%rh_litter) - &
            kv_soil*soil%rh_soil*d2est*soil%gsoil/k_soil(1))
       !       ctt = llout - soil%rnet + soil%gsoil
       !         + lecoef*(__max(soil%rh_soil*est-ea, 0.)*kv_soil+soil%c_litterevap*__max(soil%rh_litter*est-ea, 0.)*kv_litter)
       !         + lecoef*kv_soil*soil%rh_soil*(d2est/2.*pow(soil%gsoil/k_soil(1), 2) - dest*soil%gsoil/k_soil(1))
       ctt = llout - soil%rnet + soil%gsoil + &
            lecoef*((soil%rh_soil*est-ea)*kv_soil+soil%c_litterevap*(soil%rh_litter*est-ea)*kv_litter) + &
            lecoef*kv_soil*soil%rh_soil*(half*d2est*soil%gsoil/k_soil(1)*soil%gsoil/k_soil(1) - &
            dest*soil%gsoil/k_soil(1))
       product = btt*btt - 4._wp*att*ctt
       if (product >= zero) then
          soil%tsrf = soil%T_air + (-btt+sqrt(product))/(two*att)
       else
          soil%tsrf = soil%T_air
       end if
       ! ToDo for isotopes ! litter evaporation from second Taylor expansion
       ! soil%litterevap = soil%c_litterevap*lecoef*kv_litter*(__max(soil%rh_litter*est-ea, 0.)
       !                                                                     + soil%rh_litter*dest*(soil%tsrf-soil%T_air)
       !                                                                     + soil%rh_litter*d2est/2.*pow(soil%tsrf-soil%T_air, 2)
       !                                                                     )
       ! litter evaporation from energy budget
       if (1 == 1) then
          soil%litterevap = soil%c_litterevap*lecoef*kv_litter * &
               (soil%rh_litter*ES(soil%tsrf+TN0)*100._wp - &
               prof%rhov_air_filter(1,1)*(soil%T_air+TN0)*Rw)
       else
          soil%litterevap = soil%c_litterevap*lecoef*kv_litter * &
               ((soil%rh_litter*est-ea) + soil%rh_litter*dest*(soil%tsrf-soil%T_air) + &
               soil%rh_litter*half*d2est*(soil%tsrf-soil%T_air)*(soil%tsrf-soil%T_air))
       end if
       if (iswitch%no_neg_water_flux == 1) soil%litterevap = max(soil%litterevap, zero)
       ! ToDo for isotopes ! soil evaporation from second Taylor expansion
       ! soil%soilevap = lecoef*kv_soil*(__max(soil%rh_soil*est-ea, 0.)
       !                                        + soil%rh_soil*dest*(temp-soil%T_air)
       !                                        + soil%rh_soil*d2est/2.*pow(temp-soil%T_air, 2)
       !                                        )
       ! soil evaporation from energy budget
       if (1 == 1) then
          soil%soilevap = lecoef*kv_soil * (soil%rh_soil*ES(soil%tsrf+TN0)*100._wp - &
               prof%rhov_air_filter(1,1)*(soil%T_air+TN0)*Rw)
       else
          temp = soil%tsrf - soil%gsoil/k_soil(1)
          soil%soilevap = lecoef*kv_soil * ((soil%rh_soil*est-ea) + soil%rh_soil*dest*(temp-soil%T_air) + &
               soil%rh_soil*half*d2est*(temp-soil%T_air)*(temp-soil%T_air))
       end if
       if (iswitch%no_neg_water_flux == 1) soil%soilevap = max(soil%soilevap, zero)
       soil%evap = soil%litterevap + soil%soilevap
       ! IR emissive flux density from soil (W m-2)
       !soil%lout = epsoil * sigma * pow(soil%tsrf+TN0, 4)
       soil%lout = llout + 4._wp*epsoil*sigma*tk3*(soil%tsrf-soil%T_air) + 6._wp*epsoil*sigma*tk2 * &
            (soil%tsrf-soil%T_air)*(soil%tsrf-soil%T_air)
       ! Sensible heat flux density over soil (W m-2)
       soil%heat = kcsoil*(soil%tsrf-soil%T_air)
       ! update latent heat
       !soil%latent = LAMBDA(soil%tsrf_filter+TN0)
       !lecoef = met%air_density * 0.622 * soil%latent / met%press_Pa
       ! set new temperature profile, check for extremes
       where (T_new_soil(0:n_soil_end) < -30._wp .or. T_new_soil(0:n_soil_end) > 60._wp) &
            T_new_soil(0:n_soil_end) = input%ta
       T_soil(0:n_soil_end) = T_new_soil(0:n_soil_end)
    end do ! j of soil%temperature_mtime

    ! set 'real' soil & litter temperatures
    if (soil%z_litter < 1e-6_wp) then
       soil%T_l = soil%tsrf + TN0
    else
       soil%T_l = T_new_soil(1)+TN0
    end if
    soil%T_soil(0:n_soil_end) = T_new_soil(1:n_soil_end+1)
    ! store coefficient of latent heat
    soil%lecoef = lecoef
    ! for output
    output%c1 = kcsoil
    output%c2 = kv_soil
    output%c3 = Rh_soil
    output%c4 = soil%rb
    output%c5 = met%air_density
    ! Restrict litterevap to maximum available litter water content
    if (soil%z_litter < 1e-6_wp) then
       soil%c_litterevap_save = zero
       soil%c_litterevap      = zero
    else
       soil%c_litterevap_save = soil%c_litterevap
       if (soil%litterevap > zero) then
          !* soil%latent ! kg/m2/s = mm/s -> W/m2
          soil%maxlitterevap = (1000._wp*soil%z_litter/time%time_step * soil%theta_l(1) - & ! water in litter
               min(max(zero, 1000._wp*soil%z_litter/time%time_step*(soil%theta_l(1)-soil%theta_l33)), & ! - drainage in soil
               soil%qinfl_max))* soil%latent ! kg/m2/s = mm/s -> W/m2
          if (abs(soil%c_litterevap) > tiny(one)) then
             temp = soil%litterevap / soil%c_litterevap
          else
             temp = zero
          end if
          soil%c_litterevap = min(one, soil%maxlitterevap/max(temp, 1e-20_wp))
          soil%maxlitter = 0
          if (soil%litterevap > soil%maxlitterevap) then
             soil%litterevap_save = soil%litterevap
             soil%litterevap      = soil%maxlitterevap
             soil%maxlitter       = 1
#ifdef DEBUG
             call message('SOIL_ENERGY_BALANCE: ','Litterevap', num2str(soil%litterevap), &
                  ' > Maxlitterevap ', num2str(soil%maxlitterevap), ' at ', num2str(time%daytime))
#endif
          end if
       end if
    end if
    ! W/m2 -> kg/m2/s = mm/s
    flux%soilevap(1)   = soil%soilevap / soil%latent
    flux%litterevap(1) = soil%litterevap / soil%latent
    flux%s_evap(1)     = flux%litterevap(1) + flux%soilevap(1)

  END SUBROUTINE soil_energy_balance


  ! ------------------------------------------------------------------
  SUBROUTINE soil_h2o()
    ! * ------------------------ code history ------------------------------
    ! * source file: soih2o.F
    ! * purpose: soil hydrology and sub-surface drainage
    ! * date last revised: March 1996 - lsm version 1
    ! * author: Gordon Bonan
    ! * standardized: J. Truesdale, Feb. 1996
    ! * reviewed: G. Bonan, Feb. 1996
    ! * --------------------------------------------------------------------

    ! * ------------------------ notes -------------------------------------
    ! * use tridiagonal system of equations to solve one-dimensional water balance:

    ! * d wat
    ! * dz ----- = -qi + qo - s
    ! * dt

    ! * with q = -k d(psi+z)/dz = -k (d psi/dz + 1) and with s=0
    ! * this is the Richards equation for vertical water flow

    ! * d wat d d wat d psi
    ! * ----- = -- [ k(----- ----- + 1) ]
    ! * dt dz dz d wat

    ! * where: wat = volume of water per volume of soil (mm**3/mm**3)
    ! * psi = soil matrix potential (mm)
    ! * dt = time step (s)
    ! * z = depth (mm)
    ! * dz = thickness (mm)
    ! * qi = inflow at top (mm h2o /s) (+ = up, - = down)
    ! * qo = outflow at bottom (mm h2o /s) (+ = up, - = down)
    ! * s = source/sink flux (mm h2o /s) (+ = loss, - = gain)
    ! * k = hydraulic conductivity (mm h2o /s)

    ! * solution: linearize k and psi about d wat and use tridiagonal system
    ! * of equations to solve for d wat, where for layer i
    ! * r_i = a_i [d wat_i-1] + b_i [d wat_i] + c_i [d wat_i+1]

    ! * the solution conserves water as:
    ! * [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n+1 =
    ! * [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n +
    ! * (qinfl - qseva - qtran - qdrai)*dtlsm

    ! * --------------------------------------------------------------------
    ! * ------------------------ input/output variables --------------------
    ! * input
    USE constants,     ONLY: zero, one, e2, e3, isotolerance, s2h, Gravity
    USE setup,         ONLY: nsoil, nwiso
    USE types,         ONLY: flux, soil, wiso, time
    USE isotope_utils, ONLY: isorat
    USE messages,      ONLY: message
    USE string_utils,  ONLY: num2str


    IMPLICIT NONE

    REAL(wp), DIMENSION(1:nwiso) :: qtran !transpiration water flux (mm h2o/s)
    REAL(wp), DIMENSION(1:nwiso) :: qinfl !infiltration rate (mm h2o/s)
    REAL(wp), DIMENSION(1:nwiso) :: qseva !ground surface evaporation rate (mm h2o/s)
    REAL(wp), DIMENSION(1:nsoil) :: dzsoi !soil layer thickness (m)
    REAL(wp), DIMENSION(1:nsoil) :: root !relative root abundance (0 to 1)
    !* input/output
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: h2osoi !volumetric soil water content
    REAL(wp), DIMENSION(1:nwiso) :: lost_h2osoi ! water conservation check variable
    !* output
    REAL(wp), DIMENSION(1:nwiso) :: qdrai !sub-surface runoff (mm h2o/s)
    !* ------------------------ local variables ---------------------------
    INTEGER(i4) :: j !do loop/array indices
    INTEGER(i4) :: iter !loop
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: r !solution matrix
    REAL(wp), DIMENSION(1:nsoil) :: a !"a" vector for tridiagonal matrix
    REAL(wp), DIMENSION(1:nsoil) :: b !"b" vector for tridiagonal matrix
    REAL(wp), DIMENSION(1:nsoil) :: c !"c" vector for tridiagonal matrix
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: dwat !change in soil water
    REAL(wp), DIMENSION(1:nsoil) :: smp !soil matrix potential (mm)
    REAL(wp), DIMENSION(1:nsoil) :: hk !hydraulic conductivity (mm h2o/s)
    REAL(wp), DIMENSION(1:nsoil) :: hk2 !hk**2
    REAL(wp), DIMENSION(1:nsoil) :: dsmpdw !d(smp)/d(h2osoi)
    REAL(wp), DIMENSION(1:nsoil) :: dhkdw !d(hk)/d(h2osoi)
    REAL(wp), DIMENSION(1:nsoil) :: s ! relative soil water content in Clapp & Hornberger
    REAL(wp), DIMENSION(1:nsoil) :: hydcon ! temp hydraulic conductivity
    REAL(wp) :: qin !flux of water into soil layer (mm h2o/s)
    REAL(wp) :: qout !flux of water out of soil layer (mm h2o/s)
    REAL(wp) :: num !used in calculating qi, qout
    REAL(wp) :: den !used in calculating qi, qout
    REAL(wp) :: den2 !den**2 used in calculating qi, qout
    REAL(wp) :: dqidw1 !d(qin)/d(h2osoi(i-1))
    REAL(wp) :: dqidw2 !d(qin)/d(h2osoi(i))
    REAL(wp) :: dqodw1 !d(qout)/d(h2osoi(i))
    REAL(wp) :: dqodw2 !d(qout)/d(h2osoi(i+1))
    REAL(wp) :: xs !soil water > sat or < some minimum (mm h2o)
    REAL(wp) :: thisxs
    REAL(wp) :: hkdrai
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: rsoil ! ratio of levels
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: xstmp1
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: xstmp2
    REAL(wp), DIMENSION(1:nsoil) :: ztmp
    REAL(wp) :: ztmp1
    INTEGER(i4), DIMENSION(1) :: j55

    ! Assign local variables
    qtran(1:nwiso)          = flux%c_transpiration(1:nwiso) !transpiration water flux (kg H2O m-2 s-1) = mm h2o/s
    soil%qtran(1:nwiso)     = flux%c_transpiration(1:nwiso) !transpiration water flux (kg H2O m-2 s-1) = mm h2o/s
    qinfl(1:nwiso)          = soil%qinfl(1:nwiso)
    soil%qseva(1:nwiso)     = flux%soilevap(1:nwiso) !ground surface evaporation rate (mm h2o/s) = (kg H2O m-2 s-1)
    qseva(1:nwiso)          = flux%soilevap(1:nwiso) !ground surface evaporation rate (mm h2o/s) = (kg H2O m-2 s-1)
    soil%qdrai(1:nwiso)     = zero
    qdrai(1:nwiso)          = zero
    ! set effective soil layer thickness
    dzsoi(1:nsoil)          = (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1)) * (one-soil%gravel(1:nsoil)*e2)
    root(1:nsoil)           = soil%root(1:nsoil) !relative root abundance (0 to 1)
    ! volumetric soil water content, get soil water content from last time step
    h2osoi(1:nsoil,1:nwiso) = soil%theta(1:nsoil,1:nwiso)
    lost_h2osoi(1:nwiso)    = wiso%lost(1:nwiso)

    do iter=1, soil%moisture_mtime
       rsoil(1:nsoil,1:nwiso) = isorat(h2osoi(1:nsoil,1:nwiso), h2osoi(1:nsoil,1), &
            lost_h2osoi(1:nwiso), lost_h2osoi(1))
#ifdef DEBUG
       if (any(abs(lost_h2osoi(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('SOIL_H2O: ', 'Lost 11.1 @  ', num2str(iter), num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       ! * evaluate hydraulic conductivity, soil matrix potential,
       ! * d(smp)/d(h2osoi), and d(hk)/d(h2osoi). limit s >= 0.05
       ! * when evaluating these terms. this helps prevent numerical
       ! * problems for very small h2osoi. also limit hk >= some very
       ! * small number for same reason.
       if (soil%saxton == 1) then
          where (h2osoi(1:nsoil,1) <= soil%theta_33(1:nsoil)) ! swp in kPa
             smp(1:nsoil) = soil%big_a(1:nsoil) * &
                  max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))**(-soil%big_b(1:nsoil))
          elsewhere
             smp(1:nsoil) = 33._wp - (33._wp-soil%psi_e(1:nsoil)) * &
                  (min(max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil)), soil%theta_s(1:nsoil)) - &
                  soil%theta_33(1:nsoil)) / (soil%theta_s(1:nsoil)-soil%theta_33(1:nsoil))
          end where
          soil%swp(1:nsoil) = smp(1:nsoil) !kPa
          where (h2osoi(1:nsoil,1) <= soil%theta_33(1:nsoil)) ! derivative of swp to sm
             ! same as Clapp & Hornberger
             dsmpdw(1:nsoil) = -soil%big_b(1:nsoil) * smp(1:nsoil) / &
                  max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))
          elsewhere
             dsmpdw(1:nsoil) = -(33._wp-soil%psi_e(1:nsoil)) / (soil%theta_s(1:nsoil)-soil%theta_33(1:nsoil))
          end where
          smp(1:nsoil)          = -smp(1:nsoil)*(1000._wp/Gravity) ! kPa -> mm
          soil%swp_mm(1:nsoil)  = smp(1:nsoil)
          dsmpdw(1:nsoil)       = -dsmpdw(1:nsoil)*(1000._wp/Gravity) ! kPa -> mm
          ztmp(1:nsoil)         = 3._wp + 2._wp / soil%lambda(1:nsoil)
          hk(1:nsoil)           = soil%k_s(1:nsoil) * &
               (max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))/soil%theta_s(1:nsoil))**ztmp(1:nsoil)
          soil%k_theta(1:nsoil) = hk(1:nsoil)  ! mm h-1
          ! derivative of hk on sm
          dhkdw(1:nsoil)        = ztmp(1:nsoil) * hk(1:nsoil) / max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))
          hk(1:nsoil)           = hk(1:nsoil) * s2h ! mm s-1
          soil%k_theta(1:nsoil) = hk(1:nsoil)
          dhkdw(1:nsoil)        = dhkdw(1:nsoil) * s2h
          hk2(1:nsoil)          = hk(1:nsoil)*hk(1:nsoil) ! hk^2
       else
          ! Clapp & Hornberger
          s(1:nsoil)           = max(h2osoi(1:nsoil,1)/soil%theta_s(1:nsoil), 0.05_wp) ! % de saturation
          ! soil water potential swp (mm water)
          smp(1:nsoil)         = soil%psi_e(1:nsoil) * s(1:nsoil)**(-soil%big_b(1:nsoil))
          soil%swp_mm(1:nsoil) = smp(1:nsoil) ! mm
          soil%swp(1:nsoil)    = -smp(1:nsoil)*Gravity*e3 ! kPa
          ! derivative of swp on sm
          dsmpdw(1:nsoil)      = -soil%big_b(1:nsoil)*smp(1:nsoil)/max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))
          ztmp(1:nsoil)        = 2._wp * soil%big_b(1:nsoil) + 3._wp
          hydcon(1:nsoil)      = soil%k_s(1:nsoil) * s(1:nsoil)**ztmp(1:nsoil) ! hydraulic conductivity (mm h2o/s)
          hk(1:nsoil)          = max(hydcon(1:nsoil), 1.e-10_wp) ! hydraulic conductivity hk is not < 0.0
          ! derivative of hk on sm
          dhkdw(1:nsoil)       = hydcon(1:nsoil)*ztmp(1:nsoil)/max(h2osoi(1:nsoil,1),soil%watmin(1:nsoil))
          hk2(1:nsoil)         = hk(1:nsoil)*hk(1:nsoil) ! hk^2
       end if
       !
       !* set up r, a, b, and c vectors for tridiagonal solution
       !* node j=1
       j = 1
       num     = -2._wp*(smp(j)-smp(j+1)) - 1000._wp*(dzsoi(j)+dzsoi(j+1))
       den     = 1000._wp*(dzsoi(j)/hk(j) + dzsoi(j+1)/hk(j+1))
       den2    = den*den
       qout    = num/den
       dqodw1  = (-2._wp*den*dsmpdw(j ) + num*1000._wp*dzsoi(j )/hk2(j )*dhkdw(j ))/ den2
       dqodw2  = ( 2._wp*den*dsmpdw(j+1) + num*1000._wp*dzsoi(j+1)/hk2(j+1)*dhkdw(j+1))/ den2
       r(j,1:nwiso) = qseva(1:nwiso) + rsoil(j,1:nwiso)*qtran(1)*root(j) - &
            qinfl(1:nwiso) - qout*rsoil(merge(j,j+1,qout<=zero),1:nwiso)
       a(j) = zero
       b(j) = dqodw1 - 1000._wp*dzsoi(j) / soil%moisture_dt
       c(j) = dqodw2
       !* node j=n_soil
       j = nsoil
       num    = -2._wp*(smp(j-1)-smp(j)) -1000._wp*(dzsoi(j-1)+dzsoi(j))
       den    = 1000._wp*(dzsoi(j-1)/hk(j-1) + dzsoi(j)/hk(j))
       den2   = den*den
       qin    = num/den
       dqidw1 = (-2._wp*den*dsmpdw(j-1) +num*1000._wp*dzsoi(j-1)/hk2(j-1)*dhkdw(j-1))/ den2
       dqidw2 = ( 2._wp*den*dsmpdw(j ) +num*1000._wp*dzsoi(j )/hk2(j )*dhkdw(j ))/ den2
       qout   = -hk(j)
       dqodw1 = -dhkdw(j)
       r(j,1:nwiso) = rsoil(j,1:nwiso)*qtran(1)*root(j) + &
            qin*rsoil(merge(j-1,j,qin<=zero),1:nwiso) - qout*rsoil(nsoil,1:nwiso)
       a(j) = -dqidw1
       b(j) = dqodw1 - dqidw2 - 1000._wp*dzsoi(j)/soil%moisture_dt
       c(j) = zero
       !* nodes j=2 to j=n_soil-1
       do j=2, nsoil-1
          num    = -2._wp*(smp(j-1)-smp(j)) -1000._wp*(dzsoi(j-1)+dzsoi(j))
          den    = 1000._wp*(dzsoi(j-1)/hk(j-1) + dzsoi(j)/hk(j))
          den2   = den*den
          qin    = num/den
          dqidw1 = (-2._wp*den*dsmpdw(j-1) +num*1000._wp*dzsoi(j-1)/hk2(j-1)*dhkdw(j-1))/ den2
          dqidw2 = ( 2._wp*den*dsmpdw(j ) +num*1000._wp*dzsoi(j )/hk2(j )*dhkdw(j ))/ den2
          num    = -2._wp*(smp(j)-smp(j+1)) -1000._wp*(dzsoi(j)+dzsoi(j+1))
          den    = 1000._wp*(dzsoi(j)/hk(j) + dzsoi(j+1)/hk(j+1))
          den2   = den*den
          qout   = num/den
          dqodw1 = (-2._wp*den*dsmpdw(j ) +num*1000._wp*dzsoi(j )/hk2(j )*dhkdw(j ))/ den2
          dqodw2 = ( 2._wp*den*dsmpdw(j+1) +num*1000._wp*dzsoi(j+1)/hk2(j+1)*dhkdw(j+1))/ den2
          r(j,1:nwiso) = rsoil(j,1:nwiso)*qtran(1)*root(j) + &
               qin*rsoil(merge(j-1,j,qin<=zero),1:nwiso) - qout*rsoil(merge(j,j+1,qout<=zero),1:nwiso)
          a(j) = -dqidw1
          b(j) = dqodw1 - dqidw2 - 1000._wp*dzsoi(j)/soil%moisture_dt
          c(j) = dqodw2
       end do
       !* solve for dwat: a, b, c, r, dwat
       soil%a(1:nsoil)         = a(1:nsoil)
       soil%b(1:nsoil)         = b(1:nsoil)
       soil%c(1:nsoil)         = c(1:nsoil)
       soil%r(1:nsoil,1:nwiso) = r(1:nsoil,1:nwiso)
       call tridia()
       ! * could now update h2osoi = h2osoi + dwat except for one problem:
       ! * need to make sure h2osoi <= watsat. if not, what to do with
       ! * excess water?
       ! * We percolate any excess water down the soil% This is different from ISOLSM.
       ! * Any excess water in the last soil layer is sub-surface runoff
       dwat(1:nsoil,1:nwiso)   = soil%dwat(1:nsoil,1:nwiso)
       h2osoi(1:nsoil,1:nwiso) = h2osoi(1:nsoil,1:nwiso) + dwat(1:nsoil,1:nwiso)
       rsoil(1:nsoil,1:nwiso)  = isorat(h2osoi(1:nsoil,1:nwiso), h2osoi(1:nsoil,1), &
            lost_h2osoi(1:nwiso), lost_h2osoi(1))
#ifdef DEBUG
       if (any(abs(lost_h2osoi(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('SOIL_H2O: ', 'Lost 11.3 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       ! from 1 to n_soil-1
       do j=1, nsoil-1
          thisxs              = max(h2osoi(j,1)-soil%theta_s(j), zero) * 1000._wp * dzsoi(j)
          xstmp1(j,1:nwiso)   = thisxs * rsoil(j,1:nwiso)
          h2osoi(j  ,1:nwiso) = h2osoi(j  ,1:nwiso) - xstmp1(j,1:nwiso)/(1000._wp*dzsoi(j  ))
          h2osoi(j+1,1:nwiso) = h2osoi(j+1,1:nwiso) + xstmp1(j,1:nwiso)/(1000._wp*dzsoi(j+1))
       end do
       ! nsoil
       thisxs                = max(h2osoi(nsoil,1)-soil%theta_s(nsoil), zero) * 1000._wp * dzsoi(nsoil)
       xstmp1(nsoil,1:nwiso) = thisxs * rsoil(nsoil,1:nwiso)
       h2osoi(nsoil,1:nwiso) = h2osoi(nsoil,1:nwiso) - xstmp1(nsoil,1:nwiso)/(1000._wp*dzsoi(nsoil))
       qdrai(1:nwiso)        = qdrai(1:nwiso) + xstmp1(nsoil,1:nwiso)
       ! sub-surface drainage (accumulate over time%time_step/soil%moisture_dt iterations)
       ! isotopes assume ratio of deepest soil
       hkdrai         = (hk(nsoil)+dhkdw(nsoil)*dwat(nsoil,1))*soil%moisture_dt
       qdrai(1:nwiso) = qdrai(1:nwiso) + hkdrai * rsoil(nsoil,1:nwiso)
    end do ! iter=1, moisture_mtime
    ! Quick check for numerics
    if (any(h2osoi(1:nsoil,1) < zero)) then
       call message('SOIL_H2O: ', 'Soil water < zero @ ', num2str(time%daytime))
    end if
    !* sub-surface drainage over time step = time%time_step
    ztmp1 = one / time%time_step
    qdrai(1:nwiso) = qdrai(1:nwiso) * ztmp1
    ! limit h2osoi >= watmin. get water needed to bring h2osoi = watmin
    ! from lower layer. do for all points so inner loop vectorizes
    ! This is a "recharge" hack. xs is now a deficit if you will, and > 0.
    ! For isotopes assume this has the isotopic ratio of the lower layer.
    rsoil(1:nsoil,1:nwiso) = isorat(h2osoi(1:nsoil,1:nwiso), h2osoi(1:nsoil,1), &
         lost_h2osoi(1:nwiso), lost_h2osoi(1))
#ifdef DEBUG
    if (any(abs(lost_h2osoi(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('SOIL_H2O: ', 'Lost 11.4 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif
    ! 1 to n_soil-1
    do j=1, nsoil-1
       xs                  = max(soil%watmin(j)-h2osoi(j,1), zero)*1000._wp*dzsoi(j)
       xstmp2(j,1:nwiso)   = xs * rsoil(j,1:nwiso)
       h2osoi(j,1:nwiso)   = h2osoi(j,1:nwiso)   + xstmp2(j,1:nwiso)/(1000._wp*dzsoi(j))
       h2osoi(j+1,1:nwiso) = h2osoi(j+1,1:nwiso) - xstmp2(j,1:nwiso)/(1000._wp*dzsoi(j+1))
    end do
    ! n_soil
    xs                    = max(soil%watmin(nsoil)-h2osoi(nsoil,1), zero)*1000._wp*dzsoi(nsoil)
    xstmp2(nsoil,1:nwiso) = xs * rsoil(nsoil,1:nwiso)
    h2osoi(nsoil,1:nwiso) = h2osoi(nsoil,1:nwiso) + xstmp2(nsoil,1:nwiso)/(1000._wp*dzsoi(nsoil))
    qdrai(1:nwiso)        = qdrai(1:nwiso) - xstmp2(nsoil,1:nwiso)/time%time_step
    if (any(qdrai(1:nwiso) < zero) .and. (soil%drain0 == 0)) then
       soil%drain0 = 1
#ifdef DEBUG
       call message('SOIL_H2O: ', 'Drainage <0 started at time step: ', num2str(time%daytime))
#endif
    end if
    ! Check
    if (any((soil%watmin(1:nsoil)-h2osoi(1:nsoil,1)) > isotolerance)) then
       ! should not come to this
       rsoil(1:nsoil,1:nwiso) = isorat(h2osoi(1:nsoil,1:nwiso), h2osoi(1:nsoil,1), &
            lost_h2osoi(1:nwiso), lost_h2osoi(1))
#ifdef DEBUG
       if (any(abs(lost_h2osoi(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('SOIL_H2O: ', 'Lost 11.4 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       call message('SOIL_H2O: ', 'Soil water in layer < watmin at time step: ', num2str(time%daytime))
       h2osoi(1:nsoil,1) = max(h2osoi(1:nsoil,1), soil%watmin(1:nsoil))
       h2osoi(1:nsoil,1:nwiso) = rsoil(1:nsoil,1:nwiso)*spread(h2osoi(1:nsoil,1),2,nwiso)
    end if
    ! Calc total soil water and plant available water
    soil%soil_mm      = sum(h2osoi(1:nsoil,1) *1000._wp *dzsoi(1:nsoil)) ! gravel included in dzsoi
    soil%soil_mm_root = sum(soil%root(1:nsoil) *h2osoi(1:nsoil,1) *1000._wp *dzsoi(1:nsoil))
    j55 = minloc(abs(soil%z_soil(1:nsoil)-0.55_wp))
    soil%soil_mm_50   = sum(h2osoi(1:j55(1),1) *1000._wp *dzsoi(1:j55(1)))
    ! copy back arrays
    soil%qdrai(1:nwiso) = qdrai(1:nwiso)
    wiso%lost(1:nwiso)  = lost_h2osoi(1:nwiso)

  END SUBROUTINE soil_h2o


  ! ------------------------------------------------------------------
  FUNCTION soil_resistance()
    USE types,      ONLY: soil
    USE constants,  ONLY: zero, TN0

    IMPLICIT NONE

    REAL(wp) :: soil_resistance

    REAL(wp), PARAMETER :: c = 13.515_wp ! soil resistance of Passerat (1986), from Mahfouf & Noilhan (1991)
    REAL(wp), PARAMETER :: rsmax = 3.8113e4_wp
    REAL(wp), PARAMETER :: dw_0 = 0.229e-4_wp ! vapour diffusion in air @ TN0
    REAL(wp) :: dw, rs, rl

    rs = rsmax * exp(-c*(soil%theta(1,1)-soil%watmin(1)) / (soil%theta_33(1)-soil%watmin(1)))
    ! Camillo and Gurney model for soil resistance
    ! Rsoil= 4104 (ws-wg)-805, ws=.395, wg=0
    ! ws = max soil moisture = 0.55 at the Hainich site, wg=0
    if (soil%camillo == 1) rs = 1452.2_wp
    if (soil%z_litter < 1e-6_wp) then
       rl = zero
    else
       dw = dw_0 * (soil%T_l_filter/TN0)**1.75_wp ! from Kondo et al. (1990)
       ! with Milly model (5/3) after Saravanapavan et al. (2000).
       rl = soil%z_litter/(dw*max(soil%theta_ls-soil%theta_l(1),1e-12_wp)**(5._wp/3._wp))
    end if
    soil_resistance = rs + rl

  END FUNCTION soil_resistance


  ! ------------------------------------------------------------------
  FUNCTION soil_rh()
    ! Soil characteristics after Clapp & Hornberger (1978)
    USE types,      ONLY: soil
    USE constants,  ONLY: zero, one, TN0, Gravity, Rw

    IMPLICIT NONE

    REAL(wp) :: soil_rh

    REAL(wp) :: s, smp

    if (soil%camillo == 1) then
       soil_rh = one
    else
       if (soil%saxton == 1) then
          if (soil%theta(1,1) <= soil%theta_33(1)) then ! swp in kPa
             smp = soil%big_a(1) * soil%theta(1,1)**(-soil%big_b(1))
          else
             smp = 33._wp - (33._wp-soil%psi_e(1)) * &
                  (min(soil%theta(1,1),soil%theta_s(1))-soil%theta_33(1)) / &
                  (soil%theta_s(1)-soil%theta_33(1))
          end if
          smp = -smp*1000._wp/Gravity ! kPa -> mm
       else
          s = max(soil%theta(1,1)/soil%theta_s(1), 0.05_wp) ! % saturation
          smp = soil%psi_e(1) * s**(-soil%big_b(1)) ! soil water potential swp (mm water)
       end if
       soil_rh = exp(-Gravity*smp/(1000._wp*Rw*(soil%T_soil_filter(1)+TN0)))
    end if
    if (soil_rh > one)  soil_rh = one
    if (soil_rh < zero) soil_rh = zero

  END FUNCTION soil_rh


  ! ------------------------------------------------------------------
  SUBROUTINE throughfall()
    ! function that calculates throughfall based on
    ! input.ppt : precipitation above canopy (mm)
    ! prof%dLAIdz: LAI per layer (m2 m-2)
    ! markov : clumping facotr ()
    ! tau_water : interception efficiency
    ! cws_max : maximum canopy water storage per m2 LAI
    !
    ! output
    ! prof%throughfall : throughfall (mm)  that leaves this layer
    ! prof%intercept : interception (mm)   that is intercepted in this layer
    ! prof%cws : canopy water storage (mm) that is stored in this layers leaves
    USE constants,     ONLY: zero, one, two
    USE types,         ONLY: prof, input, wiso
    USE setup,         ONLY: ncl, nwiso
    USE parameters,    ONLY: water_film_thickness, tau_water, markov
    USE isotope_utils, ONLY: isorat
#ifdef DEBUG
    USE types,         ONLY: soil, time
    USE string_utils,  ONLY: num2str
    USE messages,      ONLY: message
#endif

    IMPLICIT NONE

    INTEGER(i4) :: i         ! counting variable for layers
    REAL(wp) :: cws_max      ! maximum canopy water storage per m2 LAI (mm m-2)
    REAL(wp) :: drip         ! water dripping of leaves after canopy water storage has reached max (mm)
    REAL(wp) :: intercept    ! interception of rain water in layer (mm)
    REAL(wp), DIMENSION(nwiso) :: rthrough ! isotope ratios
    REAL(wp), DIMENSION(nwiso) :: rcws

    cws_max   = zero   ! maximum canopy water storage per m2 LAI (mm m-2)
    drip      = zero      ! water dripping of leaves after canopy water storage has reached max (mm)
    intercept = zero ! interception of rain water in layer (mm)
    cws_max   = water_film_thickness*two   ! leaves can be wet from both sides (mm m-2 LAI)
    ! at canopy top: throughfall = precipitation, drip = 0
    prof%throughfall(1:ncl+1,1:nwiso) = spread(input%ppt(1:nwiso),1,ncl+1)
    ! calculate for each layer starting from top layer (i= ncl, typically 40)
    do i=ncl, 1, -1
       ! intercept per layer (mm)
       intercept = min(tau_water*markov*prof%dLAIdz(i)*prof%throughfall(i+1,1), &
            prof%throughfall(i+1,1))
       rthrough(1:nwiso) = isorat(prof%throughfall(i+1,1:nwiso), &
            prof%throughfall(i+1,1), wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
       if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('THROUGHFALL: ', 'Lost 07 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       ! add to canopy water storage in layer (sun and shade) from previous time step
       prof%cws(i,1:nwiso) = prof%cws(i,1:nwiso) + rthrough(1:nwiso) * intercept
       rcws(1:nwiso) = isorat(prof%cws(i,1:nwiso), prof%cws(i,1), &
            wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
       if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
          call message('THROUGHFALL: ', 'Lost 08 @  ', num2str(time%daytime))
          soil%lost0 = 1
       end if
#endif
       ! check if canopy water storage is at maximum
       if (prof%cws(i,1) > cws_max*prof%dLAIdz(i)) then
          drip                = prof%cws(i,1) -  cws_max * prof%dLAIdz(i)
          prof%cws(i,1:nwiso) = rcws(1:nwiso) * (cws_max * prof%dLAIdz(i))
       else
          drip = zero
       end if
       ! throughfall from this layer to next layer
       prof%throughfall(i,1:nwiso) = prof%throughfall(i+1,1:nwiso) - &
            rthrough(1:nwiso)*intercept + rcws(1:nwiso)*drip
#ifdef DEBUG
       ! Check: if this happens, we have to code some thresholds here
       if (any(prof%throughfall(i,1:nwiso) < zero)) then
          call message('THROUGHFALL: ', 'T<0 @ ', num2str(time%daytime))
       end if
       if (any(prof%cws(i,1:nwiso) < zero)) then
          call message('THROUGHFALL: ', 'C2<0 @ ', num2str(time%daytime))
       end if
#endif
       ! set coefficient for evaporation from wet leaf surface
       if (prof%cws(i,1) > zero) then
          prof%wet_coef_filter(i) = one
       else
          prof%wet_coef_filter(i) = zero
       end if
    end do

  END SUBROUTINE throughfall


  ! ------------------------------------------------------------------
  SUBROUTINE tridia()
    ! * ------------------------ code history ---------------------------
    ! * source file: tridia.F
    ! * purpose: solve tridiagonal system of equations
    ! * date last revised: March 1996 - lsm version 1
    ! * author: Gordon Bonan
    ! * standardized: J. Truesdale, Feb. 1996
    ! * reviewed: G. Bonan, Feb. 1996
    ! * ------------------------ notes ----------------------------------
    ! * solve for u given the set of equations f * u = r, where u is a
    ! * vector of length n, r is a vector of length n, and f is an n x n
    ! * tridiagonal matrix defined by the vectors a, b, c [each of length n].
    ! * a(1) and c(n) are undefined and are not referenced by the subroutine.

    ! * |b(1) c(1) 0 ... | |u(1 )| |r(1 )|
    ! * |a(2) b(2) c(2) ... | |u(2 )| |r(2 )|
    ! * | ... | * | ... | = | ... |
    ! * | ... a(n-1) b(n-1) c(n-1)| |u(n-1)| |r(n-1)|
    ! * | ... 0 a(n ) b(n )| |u(n )| |r(n )|
    USE types, ONLY: soil
    USE setup, ONLY: nsoil, nwiso

    IMPLICIT NONE

    ! * input/output variables
    REAL(wp), DIMENSION(1:nsoil) :: a ! input vector
    REAL(wp), DIMENSION(1:nsoil) :: b ! input vector
    REAL(wp), DIMENSION(1:nsoil) :: c ! input vector
    REAL(wp), DIMENSION(1:nsoil,1:nwiso) :: r ! input vector
    REAL(wp), DIMENSION(1:nsoil) :: u ! solution vector
    ! * local variables
    REAL(wp), DIMENSION(1:nsoil) :: gam
    REAL(wp) :: bet
    INTEGER(i4) :: j, mc ! loop index

    r(1:nsoil,1:nwiso) = soil%r(1:nsoil,1:nwiso)
    a(1:nsoil)         = soil%a(1:nsoil)
    b(1:nsoil)         = soil%b(1:nsoil)
    c(1:nsoil)         = soil%c(1:nsoil)
    do mc=1, nwiso
       bet  = b(1)
       u(1) = r(1,mc) / bet
       do j=2, nsoil
          gam(j) = c(j-1) / bet
          bet    = b(j) - a(j) * gam(j)
          u(j)   = (r(j,mc) - a(j)*u(j-1)) / bet
       end do
       do j=nsoil-1, 1, -1
          u(j) = u(j) - gam(j+1) * u(j+1)
       end do
       soil%dwat(1:nsoil,mc) = u(1:nsoil)
    end do

  END SUBROUTINE tridia


END MODULE soils
