MODULE transport

  ! This module contains the turbulence and boundary layer routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, rp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: boundary_resistance ! leaf boundary layer resistances for heat, water, CO2
  PUBLIC :: friction_velocity   ! updates friction velocity with new z/L
  PUBLIC :: conc                ! Concentration calculation routines for q, T and CO2
  PUBLIC :: uz                  ! wind speed computation as a function of z

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE boundary_resistance(zzz, tlf, cws, jlay)
    ! BOUNDARY_RESISTANCE

    ! This subroutine computes the leaf boundary layer
    ! resistances for heat, vapor and CO2 (s/m).

    ! Flat plate theory is used, as discussed in Schuepp (1993) and
    ! Grace and Wilson (1981).

    ! We consider the effects of turbulent boundary layers and sheltering.
    ! Schuepp''s review shows a beta factor multiplier is necessary for SH in
    ! flows with high turbulence. The concepts and theories used have been
    ! validated with our work on HNO3 transfer to the forest.

    ! Schuepp. 1993 New Phytologist 125: 477-507

    ! Diffusivities have been corrected using the temperature/Pressure algorithm in Massman (1998)
    USE types,      ONLY: prof, non_dim, input, bound_lay_res
    USE constants,  ONLY: zero, half, one, ddc, ddv, ddh, nnu, TN0, e15
    USE parameters, ONLY: lleaf, betfact
    USE messages,   ONLY: message

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: zzz  ! layer height
    REAL(wp),    INTENT(IN) :: tlf  ! leaf temp
    REAL(wp),    INTENT(IN) :: cws
    INTEGER(i4), INTENT(IN) :: jlay ! layer number

    REAL(wp) :: Re, Re5, Re8 ! Reynolds numbers
    REAL(wp) :: Sh_heat, Sh_vapor, Sh_CO2 ! Sherwood numbers
    REAL(wp) :: graf, GR25 ! Grasshof numbers
    REAL(wp) :: deltlf
    REAL(wp) :: Res_factor!, ddia, dRe
    REAL(wp) :: nnu_T_P, ddh_T_P, ddv_T_P, ddc_T_P, T_kelvin

    ! Difference between leaf and air temperature
    deltlf   = (tlf - prof%tair_filter(jlay))
    T_kelvin = prof%tair_filter(jlay) + TN0
    if (deltlf > zero) then
       graf = non_dim%grasshof * deltlf / T_kelvin
    else
       graf = zero
    end if
    nnu_T_P      = nnu * (1013._wp/input%press_mb) * (T_kelvin/TN0)**1.81_wp
    prof%u(jlay) = uz(zzz)
!    print *, prof%u(jlay)
    Re           = lleaf * prof%u(jlay) / nnu_T_P
    if (Re > zero) then
       Re5 = sqrt(Re)
    else
       Re5 = 100._wp
       call message('BOUNDARY_RESISTANCE: ','bad RE in RESHEAT')
    end if
    Re8 = Re**0.8_wp
    if (Re > 14000._wp) then
       Res_factor = 0.036_wp * Re8 * betfact
       ! turbulent boundary layer
       ! SH = .036 * Re8 * pr33*betfact
       ! SHV = .036 * Re8 * sc33*betfact
       ! SHCO2 = .036 * Re8 * scc33*betfact
       Sh_heat    = Res_factor * non_dim%pr33
       Sh_vapor   = Res_factor * non_dim%sc33
       Sh_CO2     = Res_factor * non_dim%scc33
    else
       Res_factor = 0.66_wp * Re5 * betfact
       ! laminar sublayer
       ! SH = .66 * Re5 * pr33*betfact
       ! SHV = .66 * Re5 * sc33*betfact
       ! SHCO2 = .66 * Re5 * scc33*betfact
       Sh_heat    = Res_factor * non_dim%pr33
       Sh_vapor   = Res_factor * non_dim%sc33
       Sh_CO2     = Res_factor * non_dim%scc33
       if (cws > zero) Sh_vapor = 0.66_wp * non_dim%sc33 * Re**0.4_wp ! originally cws > 0, Yuan changed 2017.0901 1e-15_wp
       !   Sh_vapor = 0.66 * non_dim%sc33 * pow(Re, 0.4) * betfact
    end if
    ! If there is free convection
    if (graf/(Re*Re) > one) then
       ! Compute Grashof number for free convection
       if (graf < 100000._wp) then
          GR25 = half * graf**0.25_wp
       else
          GR25 = 0.13_wp * graf**0.33_wp
       end if
       Sh_heat  = non_dim%pr33 * GR25
       Sh_vapor = non_dim%sc33 * GR25
       Sh_CO2   = non_dim%scc33 * GR25
    end if
    ! Correct diffusivities for temperature and pressure
    ddh_T_P = ddh * (1013._wp/input%press_mb) * (T_kelvin/TN0)**1.81_wp
    ddv_T_P = ddv * (1013._wp/input%press_mb) * (T_kelvin/TN0)**1.81_wp
    ddc_T_P = ddc * (1013._wp/input%press_mb) * (T_kelvin/TN0)**1.81_wp
    bound_lay_res%heat  = lleaf / (ddh_T_P * Sh_heat)
    bound_lay_res%vapor = lleaf / (ddv_T_P * Sh_vapor)
    bound_lay_res%co2   = lleaf / (ddc_T_P * Sh_CO2)

  END SUBROUTINE boundary_resistance


  ! ------------------------------------------------------------------
  SUBROUTINE conc(source, cncc, cref, soilflux, factor)
    ! Subroutine to compute scalar concentrations from source
    ! estimates and the Lagrangian dispersion matrix
    USE types,      ONLY: met
    USE constants,  ONLY: zero
    USE setup,      ONLY: ncl, ntl
    USE parameters, ONLY: delz, izref, ustar_ref

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN)  :: source   ! canopy sources
    REAL(wp), DIMENSION(:), INTENT(OUT) :: cncc     !
    REAL(wp),               INTENT(IN)  :: cref     ! ref concentration = upper boundary
    REAL(wp),               INTENT(IN)  :: soilflux ! soil efflux = lower boundary
    REAL(wp),               INTENT(IN)  :: factor   !

    REAL(wp), DIMENSION(1:ntl,1:ncl) :: disper, disperzl
    REAL(wp), DIMENSION(1:ntl) :: disper1, disperzl1
    REAL(wp), DIMENSION(1:ntl) :: sumcc, cc, soilbnd
    REAL(wp) :: ustfact

    ! Compute concentration profiles from Dispersion matrix
    ustfact = ustar_ref / met%ustar ! factor to adjust Dij with alternative u* values
    ! Note that disperion matrix was computed using u* = 0.405
    ! CC is the differential concentration (Ci-Cref)
    ! Ci-Cref = SUM (Dij S DELZ), units mg m-3,mole m-3, or heat m-3
    ! S = dfluxdz/DELZ
    ! note delz values cancel
    ! scale dispersion matrix according to friction velocity
    disper(1:ntl,1:ncl) = ustfact * met%dispersion(1:ntl,1:ncl) ! units s/m
    ! scale dispersion matrix according to Z/L
    ! if (met%zl < 0)
    ! disperzl = disper / (1.- .812 * met%zl)
    ! else
    ! disperzl=disper
    ! updated Dispersion matrix (Dec, 2002). New Tl and 1,000,000 particles Hainich
    if (met%zl < zero) then
       disperzl(1:ntl,1:ncl) = disper(1:ntl,1:ncl) * (0.2469_wp* met%zl -0.4701_wp)/(met%zl -0.3716_wp) !changed
    else
       disperzl(1:ntl,1:ncl) = disper(1:ntl,1:ncl)
    end if
    sumcc(1:ntl) = sum(delz * disperzl(1:ntl,1:ncl) * spread(source(1:ncl),dim=1,ncopies=ntl), dim=2)
    ! scale dispersion matrix according to Z/L
    ! case for j=1, soil level
    disper1(1:ntl) = ustfact * met%dispersion(1:ntl,1)
    if (met%zl < zero) then
       disperzl1(1:ntl) = disper1(1:ntl) * (0.2469_wp* met%zl -0.4701_wp)/(met%zl -0.3716_wp) !changed
    else
       disperzl1(1:ntl) = disper1(1:ntl)
    end if
    ! add soil flux to the lowest boundary condition
    ! convert to the units we need
    soilbnd(1:ntl) = soilflux * disperzl1(1:ntl) / factor
    cc(1:ntl)      = sumcc(1:ntl) / factor + soilbnd(1:ntl)
    ! Compute scalar profile below reference
    cncc(1:ntl) = cc(1:ntl) + (cref - cc(izref))

  END SUBROUTINE conc


  ! ------------------------------------------------------------------
  SUBROUTINE friction_velocity()
    ! this subroutine updates ustar and stability corrections
    ! based on the most recent H and z/L values
    USE types,      ONLY: met, input
    USE constants,  ONLY: zero, half, one, two, pi, vonKarman, Gravity
    USE parameters, ONLY: zm, zd, z0

    IMPLICIT NONE

    REAL(wp) :: xzl, logprod, phim

    met%zl = -(vonKarman*Gravity*met%H_filter*(zm-zd)) / &
         (met%air_density*1005._wp*met%T_Kelvin*met%ustar_filter*met%ustar_filter*met%ustar_filter) ! z/L
    ! restrict z/L within reasonable bounds to minimize numerical errors
    if (met%zl > 0.25_wp) met%zl = 0.25_wp
    if (met%zl < -3._wp)  met%zl = -3._wp
    ! calculation based on Businger et a. 1977 and Hoegstroem 1988 (see Foken 2003)
    if (met%zl < zero) then
       xzl = (one-met%zl*16._wp)**0.25_wp
       logprod = half*(one+xzl*xzl) * half*(one+xzl) *half*(one+xzl)
       if (logprod > zero .or. logprod < 100000._wp) then
          logprod = log(logprod)
       else
          logprod = -one
       end if
       phim = logprod - two*atan(xzl) + half*pi
    else
       phim = -5._wp*met%zl
    end if
    if ((log((zm-zd)/z0)-phim) > zero) then
       met%ustar = vonKarman*input%wnd/(log((zm-zd)/z0)-phim)
    else
       met%ustar = met%ustar
    end if
    if (met%ustar > two)     met%ustar = vonKarman*input%wnd/log((zm-zd)/z0)
    if (met%ustar < 0.02_wp) met%ustar = 0.02_wp

  END SUBROUTINE friction_velocity


  ! ------------------------------------------------------------------
  RECURSIVE FUNCTION uz(zzz) RESULT(speed)
    ! Wind speed inside canopy is calculated based on Campbell and Norman, 1998
    USE types,      ONLY: met, time
    USE constants,  ONLY: one, vonKarman
    USE parameters, ONLY: ht, zd, z0, lai, attfac

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: zzz ! height
    REAL(wp) :: speed

    REAL(wp) :: zh, factor, uh, aa, ulow, ustarlow

    zh = zzz/ht
    if (zh >= 0.11_wp) then
       factor = log((ht-zd)/z0)/vonKarman
       uh     = factor*met%ustar ! wind speed at canopy top
       ! calculating attenuation factor results in very large values (>8) due to the height canopy.
       ! Measured values are typically around 2.5
       aa     = (attfac-1.5_wp)*time%lai/lai + 1.5_wp ! vary attenuation factor between 1.5 and attfac (from parameterfile)
       speed  = uh * exp(aa*(zh-one)) ! calculate wind speed at height zh with attenuation
    else
       ulow     = uz(0.11_wp*ht) ! calculate wind speed at 0.11 * canopy height
       ustarlow = ulow * vonKarman / log((0.11_wp*ht*100._wp)) ! calculate ustar at 0.11*canopy height
       speed    = ustarlow / vonKarman * log((zh*ht)*100._wp) ! calculate wind speed at height zh with logarithmic wind profile
    end if
    ! keep routines from blowing up
    if (speed < 0.1_wp) speed = 0.1_wp

  END FUNCTION uz


END MODULE transport
