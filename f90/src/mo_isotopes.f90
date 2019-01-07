MODULE isotopes

  ! This module contains the carbon and water isotope routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, rp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: canopy_flux_wiso    ! H2O-isotope transpiration and evaporation from interception
  PUBLIC :: carbon_isotopes     ! C-isotope computations
  PUBLIC :: le_wiso             ! Latent heat so that it conforms with E=gt*(ei-ea)
  PUBLIC :: leaf_wiso           ! Transpiration and leaf water isotopes
  PUBLIC :: soil_flux_wiso      ! Soil and litter evaporation isotopes

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE canopy_flux_wiso()
    ! Calculates canopy water isotope transpiration and evaporation
    !   from leaf level fluxes
    USE constants,     ONLY: zero, one, TN0
    USE types,         ONLY: flux, wiso, prof, time, solar
#ifdef DEBUG
    USE types,         ONLY: soil
#endif
    USE setup,         ONLY: ncl, nwiso
    USE utils,         ONLY: lambda
    USE isotope_utils, ONLY: isorat
    USE string_utils,  ONLY: num2str
    USE messages,      ONLY: message

    IMPLICIT NONE

    INTEGER(i4) :: mc
    REAL(wp), DIMENSION(1:ncl,1:nwiso) :: rrcws, lam
    REAL(wp), DIMENSION(1:ncl)         :: mask

    ! isotope canopy transpiration
    lam = spread(lambda(prof%tair_filter_save(1:ncl)+TN0),2,nwiso) ! m/s -> J/m2s = W/m2
    prof%sun_LEstoma(1:ncl,1:nwiso) = prof%sun_trans_rtrans(1:ncl,1:nwiso) * lam(1:ncl,1:nwiso)
    prof%shd_LEstoma(1:ncl,1:nwiso) = prof%shd_trans_rtrans(1:ncl,1:nwiso) * lam(1:ncl,1:nwiso)
    rrcws(1:ncl,1:nwiso) = isorat(prof%cws(1:ncl,1:nwiso), prof%cws(1:ncl,1), &
         wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
    if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('CANOPY_FLUX_WISO: ', 'Lost 01 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif
    ! isotope canopy evaporation
    ! sun
    if (any(abs(prof%sun_LEwet(1:ncl,1)) > epsilon(one) .and. abs(prof%cws(1:ncl,1)) < tiny(one))) then
       call message('CANOPY_FLUX_WISO: ', 'Problem 01 @  ', num2str(time%daytime))
    end if
    mask(1:ncl) = zero
    where (abs(prof%sun_LEwet(1:ncl,1)) > epsilon(one) .and. abs(prof%cws(1:ncl,1)) > epsilon(one))
       mask(1:ncl) = prof%sun_LEwet(1:ncl,1)
    end where
    prof%sun_LEwet(1:ncl,1:nwiso) = rrcws(1:ncl,1:nwiso) * spread(mask(1:ncl),2,nwiso)
    ! shd
    if (any(abs(prof%shd_LEwet(1:ncl,1)) > epsilon(one) .and. abs(prof%cws(1:ncl,1)) < tiny(one))) then
       call message('CANOPY_FLUX_WISO: ', 'Problem 01 @  ', num2str(time%daytime))
    end if
    mask(1:ncl) = zero
    where (abs(prof%shd_LEwet(1:ncl,1)) > epsilon(one) .and. abs(prof%cws(1:ncl,1)) > epsilon(one))
       mask(1:ncl) = prof%shd_LEwet(1:ncl,1)
    end where
    prof%shd_LEwet(1:ncl,1:nwiso) = rrcws(1:ncl,1:nwiso) * spread(mask(1:ncl),2,nwiso)
    ! Fluxes
    do mc=2, nwiso
       ! isotope canopy evapotranspiration
       prof%dLEdz(1:ncl,mc) = prof%dLAIdz(1:ncl) * &
            (solar%prob_beam(1:ncl) * (prof%sun_LEstoma(1:ncl,mc)+prof%sun_LEwet(1:ncl,mc)) + &
            solar%prob_shd(1:ncl) * (prof%shd_LEstoma(1:ncl,mc)+prof%shd_LEwet(1:ncl,mc)))
       ! total isotope canopy transpiration and evaporation
       flux%c_evaporation(mc) = sum(prof%dLAIdz(1:ncl) * &
            (prof%sun_LEwet(1:ncl,mc)*solar%prob_beam(1:ncl) + &
            prof%shd_LEwet(1:ncl,mc)*solar%prob_shd(1:ncl)) / lam(1:ncl,mc))
       flux%c_transpiration(mc) = sum(prof%dLAIdz(1:ncl) * &
            (prof%sun_LEstoma(1:ncl,mc)*solar%prob_beam(1:ncl) + &
            prof%shd_LEstoma(1:ncl,mc)*solar%prob_shd(1:ncl)) / lam(1:ncl,mc))
    end do
    ! total flux
    flux%c_evapotranspiration(1:nwiso) = flux%c_evaporation(1:nwiso) + flux%c_transpiration(1:nwiso)
    flux%evapotranspiration(1:nwiso)   = flux%c_evapotranspiration(1:nwiso) + flux%s_evap(1:nwiso)

  END SUBROUTINE canopy_flux_wiso


  ! ------------------------------------------------------------------
  SUBROUTINE carbon_isotopes()
    ! This subroutine, ISOTOPE, computes concentrations of 13C
    ! It was developed in collaboration with David Bowling (circa Dec, 2000).

    ! This program computes2 fluxes as 13C/(13C+12C)

    ! We apply the algorithm of Farquhar et al that computes disc13 as a function of Ci/Ca.
    ! Ci/Ca is computed using the photosynthesis, respiration and leaf energy balance algorithms
    ! of CANOAK. The Farquhar photosynthesis model is coupled to the Ball-Berry stomatal conductance
    ! model.

    ! A Lagrangian random walk model is used to compute the turbulence Dispersion matrix
    ! and compute concentration profiles in and above the canopy. This scheme accounts
    ! for counter gradient transfer.

    ! Summary of variables computed:
    !   Isotopic discrimination due to photosynthesis and diffusion
    !   disc13(big Delta) = a + (b-a) ci/ca
    !   a is fractionation by diffusion, 4.4 per mill
    !   b is net fractionation by carboxylation, typically 27 per mill
    !   Isotopic discrimination in term of the isotopic content of the air and plant
    !   disc13 = (d13_air - d13_plant)/(1+d13_plant)
    !   disc13 = (Rair/Rplant-1) 1000
    !   disc13 = alpha -1
    !   Rplant = Rair/(disc13+1)
    !   Isotopic content relative to the PeeDee Standard
    !   d13_plant(little delta)= (Rplant/Rpdb-1) * 1000 (per mill)
    !   d13_air(little delta)= (Rair/Rpdb-1) * 1000 (per mill)
    !   Rpdb_13 = 0.01124 (13C/12C Farquhar et al. 1989) or
    !   Rpdb_12_13 = 0.01115 (13C/(12C+13C), Tans et al)
    !   The ratio Rpdb_13/Rpdb_12_13 = 1.00807. It can be used as a multiplier to convert the
    !     isotopic ratio.
    !   alpha= Rair/Rplant
    !   photosynthetic flux of 13C= A*Ra/(1+disc13) = A * Rplant
    !     This relation is in terms of net Ps, A, not gross Ps
    USE constants,     ONLY: zero, one, e3, mass_air, mass_13CO2, Rpdb_12C, Rpdb_CO2, undef
    USE setup,         ONLY: ncl, ntl
    USE types,         ONLY: ciso, bole, prof, soil, met, input, solar, fact
    USE isotope_utils, ONLY: invdelta1000
    USE transport,     ONLY: conc

    IMPLICIT NONE

    REAL(wp) :: R_soil_12C, R_air_12C, co2_13 !, da13, R_soil_CO2, R_air_CO2
    REAL(wp) :: b_a, al, as, ab, a, b, source_sun, source_shade
    REAL(wp) :: sun_A, shd_A, sumA, ztmp
    INTEGER(i4) :: j, timelag

    source_sun   = zero
    source_shade = zero
    sun_A        = zero
    shd_A        = zero
    sumA         = zero
    ! R_soil_C = 13C/(12C) is the isotopic ratio of the soil C flux and is calculated as
    ! the sum of heterotropic signal (constant) and autotrophic signal from prevoius day with time lag (=canopy discrimination)
    timelag = 1  ! for now we set timelag to 1 day,but depending on site could be longer, Hainich=4 days
    R_soil_12C = invdelta1000(ciso%delta_soil)*Rpdb_12C*soil%respiration_hetero + &
         invdelta1000(ciso%bigdelta_long(timelag))*Rpdb_12C*soil%respiration_auto
    ! soil respiration is for total CO2, so we need to use the
    ! Tans version of Pee Dee Bee, 13C/(12C + 13C) to compute the flux
    ! of 13 C for respiration
    !R_soil_CO2 = (ciso%delta_soil/1000.+1.)*Rpdb_CO2
    !R_soil_CO2 = R_soil_C12/(1+R_soil_C12)
    ! co2_13 is upper boundary condition for (13CO2), in ppm
    R_air_12C = invdelta1000(input%d13CO2) * Rpdb_12C ! ratio of 13C relative to 12C
    !R_air_CO2 = R_air_C12/(1+R_air_C12)  ! convert to 13C relatve to (12C + 13C)
    !   Or calculate it with Bowling method
    !   Calculate calculate d13Cair from regression of deltaCa*Ca=m*Ca+b
    !   and then solve for delta
    !   Data are coming from field measurements and coefficients are provided in parameter file
    !   Note different functions exist for day and night
    ! if(input%parin > 0)
    !     da13 = (ciso%da_m_day * input%co2air + ciso%da_b_day)/input%co2air  ! daytime
    ! else
    !     da13= (ciso%da_m_night * input%co2air + ciso%da_b_night)/input%co2air   !  nighttime
    ! R_air_12C = (da13/1000.+1.)*Rpdb_12C  ! ratio of 13C relative to CO2
    ! We have to be careful here as CO2air is 12C + 13C. We need to multiply
    ! R_13_12 times 12C to compute the concentration of (13CO2)
    ! R_13_12 * 12C = 13C
    ! C = 12C + 13C
    ! 12C = C - 13C
    ! 12C= C/(1+R_13_12)
    ! 13C = R_13_12 * 12C = R_13_12 * (C - 13C)
    ! 13C + R_13_12 * 13C = R_13_12 * C
    ! 13C = R_13_12 * C/(1+ R_13_12)
    co2_13 = R_air_12C * input%co2air ! ppm concentration of 13C
    ! Respiration flux of 13C = Rr*Resp, where Resp = soil respiration
    soil%resp_13 = R_soil_12C * soil%respiration_mole
    b_a = 27._wp - 4.4_wp
    ab  = 2.9_wp !fractionation due to diffusion in the laminar boundary layer, permille
    a   = 4.4_wp !fractionation due to diffusion from the leaves to stomata cavity, permille
    as  = 1.1_wp !fractionation as CO2 enters solution at 25 deg C
    al  = 0.7_wp !fractionation by diffusion within cell
    b   = 29._wp !fractionation during carboxylation
    ! calculate canopy discrimination (only if there is photosynthesis)
    !  if (flux.photosyn > zero) {
    ! compute cc (CO2 mixing ratio at site of carboxylation) on sun and shade fractions for layers of the canopy
    ! and use this information to compute the source/sink function of 13C
    do j=1, ncl
       ! ! get resistances and convert from m s-1 to mol m-2 s-1
       ! fact_rb = rugc * (prof%tair(j)+TN0)/ met%press_Pa !conversion factor
       ! ! J mol-1 K-1 * K * Pa-1 = J mol-1 Pa-1 = (Kg m2 s-2) mol-1 (Kg m-1 s-2)-1 = m3 mol-1
       ! fact_rs_sun = rugc * (prof%sun_tleaf_filter(j)+TN0)/ met%press_Pa !conversion factor
       ! fact_rs_shd = rugc * (prof%shd_tleaf_filter(j)+TN0)/ met%press_Pa !conversion factor
       ! rb = prof%sun_rbco2(j)*fact_rs_sun ! s m-1 to mol-1 m2 s1
       ! rs = prof%sun_rs(j)*fact_rs_sun*1.577 ! s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2)
       ! cs = prof%co2_air(j) - rb*prof%sun_A(j) !unit ppm
       ! prof%sun_ci(j) = prof%shd_cs(j) - rs*prof%sun_A(j) !unit ppm
       ! ci = cs - rs*prof%sun_A(j) !unit ppm
       ! DONE IN PHOTOSYNTHESIS NOW
       ! ! mesophyll conductance
       ! rm = 1/(0.0045*prof%vcmax(j)) ! mol-1 m2 s1 similar to Suits et al. 2005
       ! rm_epron = 80/(prof%sun_A(j)-0.05) ! according to Epron et al. 1995
       ! rm_manter = 1/(0.0028*prof%vcmax(j) - 0.0176) ! accornding to data from Manter 2004
       ! fact_rs_sun = rugc * (prof%sun_tleaf_filter(j)+TN0)/ met%press_Pa !conversion factor
       ! fact_rs_shd = rugc * (prof%shd_tleaf_filter(j)+TN0)/ met%press_Pa !conversion factor
       ! rm_gs_sun = (1/0.9)*prof%sun_rs(j)*fact_rs_sun*1.577 ! s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2), according to Manter (data, pers. comm)
       ! rm_gs_shd = (1/0.9)*prof%shd_rs(j)*fact_rs_shd*1.577 ! s m-1 to mol-1 m2 s1 and convert from rs(H2O) to rs(CO2), according to Manter (data, pers. comm)
       ! rm_sun = rm
       ! rm_shd = rm
       ! !sun leaves
       ! prof%sun_cc(j) = prof%sun_ci(j) - rm_sun*prof%sun_A(j) !unit ppm
       ! !shade leaves
       ! prof%shd_cc(j) = prof%shd_ci(j) - rm_shd*prof%shd_A(j) !unit ppm
       ! ! Cc/Ca on sun leaves
       ! prof%sun_cica(j)= prof%sun_ci(j)/prof%co2_air(j)
       ! Compute Del 13C using Farquhar eq. disc13C= a + (b-a)Ci/Ca
       ! a = 4.4 per mill, b = 27.5 per mill
       prof%sun_disc13(j) = a + b_a*min(prof%sun_cica(j), one) ! parts per mill
       ! compute Del 13C using Farquhar eq Appendix
       prof%sun_disc13_ab(j)  = ab * max((prof%co2_air_filter(j)-prof%sun_cs(j))/prof%co2_air_filter(j), zero)
       prof%sun_disc13_a(j)    = a * max((prof%sun_cs(j)-prof%sun_ci(j))/prof%co2_air_filter(j), zero)
       prof%sun_disc13_asal(j) = (as+al) * max((prof%sun_ci(j)-prof%sun_cc(j))/prof%co2_air_filter(j), zero)
       prof%sun_disc13_b(j)    = b * min((prof%sun_cc(j))/prof%co2_air_filter(j), one)
       prof%sun_disc13_long(j) = prof%sun_disc13_ab(j) + prof%sun_disc13_a(j) + &
            prof%sun_disc13_asal(j) + prof%sun_disc13_b(j)
       ! Compute Del 13C using Farquhar eq. disc13C= a + (b-a)Ci/Ca
       prof%shd_disc13(j)      = a + b_a * min(prof%shd_cica(j), one)
       prof%shd_disc13_ab(j)   = ab * max((prof%co2_air_filter(j)-prof%shd_cs(j))/prof%co2_air_filter(j), zero)
       prof%shd_disc13_a(j)    = a * max((prof%shd_cs(j)-prof%shd_ci(j))/prof%co2_air_filter(j), zero)
       prof%shd_disc13_asal(j) = (as+al) * max((prof%shd_ci(j)-prof%shd_cc(j))/prof%co2_air_filter(j), zero)
       prof%shd_disc13_b(j)    = b * min((prof%shd_cc(j))/prof%co2_air_filter(j), one)
       prof%shd_disc13_long(j) = prof%shd_disc13_ab(j) + prof%shd_disc13_a(j) + &
            prof%shd_disc13_asal(j) + prof%shd_disc13_b(j)
       ! sour13co2 is photsynthetic flux of 13C
       ! disc13 is integrated (sun+shade) discrimination for this layer
       ! need to work up disc13 from prof%sun_disc13 and prof%sh_disc13!
       ! - may need to iterate between sour13CO2 and Rair
       ! need to divide disc13 by 1000
       prof%Rplant_sun(j) = prof%R13_12_air(j) / (one+prof%sun_disc13_long(j)*e3)
       prof%Rplant_shd(j) = prof%R13_12_air(j) / (one+prof%shd_disc13_long(j)*e3)
       ! sour13CO2 source/sink strength of 13C in micromoles
       ! apply minus sign as Ps is an atmospheric sink for C
       ! fluxes on the sun/shade fractions are weighted and considered
       sun_A = max(prof%sun_A(j)*solar%prob_beam(j), zero)
       shd_A = max(prof%shd_A(j)*solar%prob_shd(j), zero)
       sumA  = sun_A + shd_A
       if (sun_A > zero .or. shd_A > zero) then
          prof%disc13C_long(j) = (shd_A*prof%shd_disc13_long(j) + sun_A*prof%sun_disc13_long(j)) / sumA
          prof%disc13C(j)      = (shd_A*prof%shd_disc13(j) + sun_A*prof%sun_disc13(j)) / sumA
          ! step of discrimination
          prof%disc13C_ab(j)   = (shd_A*prof%shd_disc13_ab(j) + sun_A*prof%sun_disc13_ab(j)) / sumA
          prof%disc13C_a(j)    = (shd_A*prof%shd_disc13_a(j) + sun_A*prof%sun_disc13_a(j)) / sumA
          prof%disc13C_asal(j) = (shd_A*prof%shd_disc13_asal(j) + sun_A*prof%sun_disc13_asal(j)) / sumA
          prof%disc13C_b(j)    = (shd_A*prof%shd_disc13_b(j) + sun_A*prof%sun_disc13_b(j)) / sumA
          ! cs, ci, cc
          prof%cs(j)           = (shd_A*prof%shd_cs(j) + sun_A*prof%sun_cs(j)) / sumA
          prof%ci(j)           = (shd_A*prof%shd_ci(j) + sun_A*prof%sun_ci(j)) / sumA
          prof%cc(j)           = (shd_A*prof%shd_cc(j) + sun_A*prof%sun_cc(j)) / sumA
          prof%csca(j)         = prof%cs(j)/prof%co2_air_filter(j)
          prof%cica(j)         = prof%ci(j)/prof%co2_air_filter(j)
          prof%ccca(j)         = prof%cc(j)/prof%co2_air_filter(j)
          source_shade         = shd_A*prof%Rplant_shd(j)
          source_sun           = sun_A*prof%Rplant_sun(j)
       else
          prof%disc13C_long(j) = undef
          prof%disc13C(j)      = undef
          ! step of discrimination
          prof%disc13C_ab(j)   = undef
          prof%disc13C_a(j)    = undef
          prof%disc13C_asal(j) = undef
          prof%disc13C_b(j)    = undef
          ! cs, ci, cc
          prof%cs(j)           = undef
          prof%ci(j)           = undef
          prof%cc(j)           = undef
          prof%csca(j)         = undef
          prof%cica(j)         = undef
          prof%ccca(j)         = undef
       end if
       ! if dark, respire only fixed 13C assume signature of soil for now
       ! latter use data assimilated from day
       if (prof%shd_A(j) <= zero) then
          source_shade = solar%prob_shd(j) * prof%shd_A(j) * prof%Rresp_ave(j)
       end if
       if (prof%sun_A(j) <= zero) then
          source_sun = solar%prob_beam(j) * prof%sun_A(j) * prof%Rresp_ave(j)
       end if
       ! micromol C m-2 s-1
       prof%sour13co2(j) = -prof%dLAIdz(j) * (source_shade+source_sun)
    end do ! j=1, ncl
    ! add bole respiration to 13CO2 sources assuming it has the canopy discrimination signal
    ! shifted by one day
    prof%sour13co2(1:ncl) = prof%sour13co2(1:ncl) + &
         bole%layer(1:ncl)*(invdelta1000(ciso%bigdelta_long(timelag))*Rpdb_12C)
    ! generate 13CO2 concentration profiles using information on 13C leaf
    ! sources and sinks, turbulent mixing of the air and source of soil
    ! respiration, each with distinct signatures
    ! need to convert umol m-3 to umol/mol we have to consider
    ! Pc/Pa = (CO2)ppm = rhoc ma/ rhoa mc
    fact%co2 = (mass_air/mass_13CO2)*met%air_density_mole
    call conc(prof%sour13co2, prof%c13cnc, co2_13, soil%resp_13, fact%co2)
    ! this is the 13C/(12C+13C) ratio of air in the canopy
    prof%R13_12_air(1:ntl) = prof%c13cnc(1:ntl)/prof%co2_air_filter(1:ntl)
    ztmp = one / Rpdb_CO2
    prof%d13Cair(1:ntl) = (prof%R13_12_air(1:ntl)*ztmp-one)*1000._wp
 !   print*, 'd13a: ', prof%R13_12_air(1), prof%c13cnc(1), prof%co2_air_filter(1), &
 !       prof%sour13co2(1), bole%layer(1), (invdelta1000(ciso%bigdelta_long(timelag))*Rpdb_12C)

  END SUBROUTINE carbon_isotopes


  ! ------------------------------------------------------------------
  SUBROUTINE le_wiso()
    ! Calculates latent heat of transpiration in accordance with leaf water isotopes
    !   cf. Cuntz et al. (2007)
    USE constants,    ONLY: zero, one, TN0, Rw
    USE types,        ONLY: prof, met, time, iswitch
    USE setup,        ONLY: ncl
    USE parameters,   ONLY: n_stomata_sides
    USE utils,        ONLY: es
    USE string_utils, ONLY: num2str
    USE messages,     ONLY: message

    IMPLICIT NONE

    INTEGER(i4) :: j

    ! isotope independent variables
    do j=1, ncl
       ! wi
       prof%sun_wi(j) = es(prof%sun_tleaf(j)+TN0)*100._wp/met%press_Pa ! (Pa)
       prof%shd_wi(j) = es(prof%shd_tleaf(j)+TN0)*100._wp/met%press_Pa ! (Pa)
       ! wa
       prof%wa(j) = prof%rhov_air_filter(j,1)*(prof%tair_filter(j)+TN0)*Rw/met%press_Pa ! (Pa)
       ! h
       if (prof%wa(j) < zero) then
          call message('LE_WISO: ','Abs. humidity wa < 0. @ ', num2str(j), num2str(time%daytime))
       end if
       if (iswitch%no_neg_water_flux == 1) then
          prof%sun_h(j) = max(min(prof%wa(j)/prof%sun_wi(j), one), zero)
          prof%shd_h(j) = max(min(prof%wa(j)/prof%shd_wi(j), one), zero)
       else
          prof%sun_h(j) = prof%wa(j)/prof%sun_wi(j)
          prof%shd_h(j) = prof%wa(j)/prof%shd_wi(j)
       end if
       ! (m/s) -> (mol(H2O)/m2s)
       prof%rs_fact(j) = n_stomata_sides * met%air_density * 0.622_wp ! m->mol
       ! gross transpiration flux (mol(H2O)/m2s)
       prof%sun_gross(j) = prof%rs_fact(j)/(prof%sun_rs_save(j)+prof%sun_rbv(j)) * prof%sun_wi(j)
       prof%shd_gross(j) = prof%rs_fact(j)/(prof%shd_rs_save(j)+prof%shd_rbv(j)) * prof%shd_wi(j)
       ! E (mol(H2O)/m2s)
       prof%sun_LEstoma_new(j) = prof%sun_gross(j)*(one-prof%sun_h(j))
       prof%shd_LEstoma_new(j) = prof%shd_gross(j)*(one-prof%shd_h(j))
    end do

  END SUBROUTINE le_wiso


  ! ------------------------------------------------------------------
  SUBROUTINE leaf_wiso()
    ! Calculates leaf water isotopic composition after Farquhar & Cernusak (2005)
    !   with fixed leaf water volume. It uses a Dongmann-style solution
    !   (Dongmann et al. 1974) for enrichment at the evaporating sites and does the
    !   'brave assumption' for bulk mesophyll water.
    ! Cf. Cuntz et al. (2007)
    USE constants,     ONLY: zero, one, TN0
    USE types,         ONLY: soil, prof, wiso, time
    USE setup,         ONLY: ncl, nwiso, nsoil
    USE isotope_utils, ONLY: alpha_kin_h2o, alpha_equ_h2o, isorat, vtf
#ifdef DEBUG
    USE string_utils,  ONLY: num2str
    USE messages,      ONLY: message
#endif

    IMPLICIT NONE

    INTEGER(i4) :: mc
    REAL(wp), PARAMETER :: Vm     = 18._wp
    REAL(wp), PARAMETER :: Leff   = 0.020_wp
    REAL(wp), PARAMETER :: waterc = 55555.55555_wp

    ! xylem
    soil%xylem(1:nwiso)  = sum(soil%theta(1:nsoil,1:nwiso)*spread(soil%root(1:nsoil),2,nwiso),1)
    soil%rxylem(1:nwiso) = isorat(soil%xylem(1:nwiso), soil%xylem(1), wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
    if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('LEAF_WISO: ', 'Lost 01 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif
    ! Variables - isotope dependent
    do mc=2, nwiso
       ! kinetic fractionation
       prof%sun_alpha_k(1:ncl,mc) = alpha_kin_h2o(prof%sun_rs_save(1:ncl), prof%sun_rbv(1:ncl), mc, wiso%merlivat)
       prof%shd_alpha_k(1:ncl,mc) = alpha_kin_h2o(prof%shd_rs_save(1:ncl), prof%shd_rbv(1:ncl), mc, wiso%merlivat)
       ! equilibrium fractionation
       prof%sun_alpha_equ(1:ncl,mc) = alpha_equ_h2o(prof%sun_tleaf(1:ncl)+TN0, mc)
       prof%shd_alpha_equ(1:ncl,mc) = alpha_equ_h2o(prof%shd_tleaf(1:ncl)+TN0, mc)
       ! Peclet, P
       prof%sun_peclet(1:ncl,mc) = prof%sun_LEstoma_new(1:ncl) * Leff / waterc / vtf(prof%sun_tleaf(1:ncl)+TN0, mc)
       prof%shd_peclet(1:ncl,mc) = prof%shd_LEstoma_new(1:ncl) * Leff / waterc / vtf(prof%shd_tleaf(1:ncl)+TN0, mc)
       ! (1-exp(-P))/P
       prof%sun_fem(1:ncl,mc) = (one-exp(-prof%sun_peclet(1:ncl,mc)))/prof%sun_peclet(1:ncl,mc)
       prof%shd_fem(1:ncl,mc) = (one-exp(-prof%shd_peclet(1:ncl,mc)))/prof%shd_peclet(1:ncl,mc)
       ! Craig-Gordon
       prof%sun_craig(1:ncl,mc) = ((one-prof%sun_h(1:ncl))*soil%rxylem(mc)/prof%sun_alpha_k(1:ncl,mc) + &
            prof%sun_h(1:ncl)*prof%rvapour(1:ncl,mc)) / prof%sun_alpha_equ(1:ncl,mc)
       prof%shd_craig(1:ncl,mc) = ((one-prof%shd_h(1:ncl))*soil%rxylem(mc)/prof%shd_alpha_k(1:ncl,mc) + &
            prof%shd_h(1:ncl)*prof%rvapour(1:ncl,mc)) / prof%shd_alpha_equ(1:ncl,mc)
       ! Dongmann
       where (prof%sun_LEstoma_new(1:ncl) /= zero)
          prof%sun_leafwater_e(1:ncl,mc) = prof%sun_craig(1:ncl,mc) - &
               (prof%sun_craig(1:ncl,mc) - prof%sun_leafwater_e_old(1:ncl,mc)) * &
               exp(-prof%sun_alpha_k(1:ncl,mc)*prof%sun_alpha_equ(1:ncl,mc) * &
               prof%sun_gross(1:ncl)*time%time_step/Vm)
       elsewhere
          prof%sun_leafwater_e(1:ncl,mc) = prof%sun_leafwater_e_old(1:ncl,mc)
       end where
       ! evaporative site
       where (prof%shd_LEstoma_new(1:ncl) /= zero)
          prof%shd_leafwater_e(1:ncl,mc) = prof%shd_craig(1:ncl,mc) - &
               (prof%shd_craig(1:ncl,mc) - prof%shd_leafwater_e_old(1:ncl,mc)) * &
               exp(-prof%shd_alpha_k(1:ncl,mc)*prof%shd_alpha_equ(1:ncl,mc) * &
               prof%shd_gross(1:ncl)*time%time_step/Vm)
       elsewhere
          prof%shd_leafwater_e(1:ncl,mc) = prof%shd_leafwater_e_old(1:ncl,mc)
       end where
       ! nofrac = steady-state
       if (wiso%nofracleaf == 1) then
          prof%sun_leafwater_e(1:ncl,mc) = prof%sun_craig(1:ncl,mc)
          prof%shd_leafwater_e(1:ncl,mc) = prof%shd_craig(1:ncl,mc)
       end if
       ! bulk mesophyll
       prof%sun_leafwater(1:ncl,mc) = prof%sun_fem(1:ncl,mc) * prof%sun_leafwater_e(1:ncl,mc) + &
            (one-prof%sun_fem(1:ncl,mc)) * soil%rxylem(mc)
       prof%shd_leafwater(1:ncl,mc) = prof%shd_fem(1:ncl,mc) * prof%shd_leafwater_e(1:ncl,mc) + &
            (one-prof%shd_fem(1:ncl,mc)) * soil%rxylem(mc)
       ! E*R_E
       where (prof%sun_LEstoma_new(1:ncl) /= zero)
          prof%sun_trans_rtrans(1:ncl,mc) = prof%sun_alpha_k(1:ncl,mc) * prof%sun_gross(1:ncl) * &
               (prof%sun_alpha_equ(1:ncl,mc)*prof%sun_leafwater_e(1:ncl,mc) - &
               prof%sun_h(1:ncl)*prof%rvapour(1:ncl,mc))
       elsewhere
          prof%sun_trans_rtrans(1:ncl,mc) = zero
       end where
       where (prof%shd_LEstoma_new(1:ncl) /= zero)
          prof%shd_trans_rtrans(1:ncl,mc) = prof%shd_alpha_k(1:ncl,mc) * prof%shd_gross(1:ncl) * &
               (prof%shd_alpha_equ(1:ncl,mc)*prof%shd_leafwater_e(1:ncl,mc) - &
               prof%shd_h(1:ncl)*prof%rvapour(1:ncl,mc))
       elsewhere
          prof%shd_trans_rtrans(1:ncl,mc) = zero
       end where
       ! R_E
       where (prof%sun_LEstoma_new(1:ncl) /= zero)
          prof%sun_rtrans(1:ncl,mc) = prof%sun_trans_rtrans(1:ncl,mc)/prof%sun_LEstoma_new(1:ncl)
       elsewhere
          prof%sun_rtrans(1:ncl,mc) = zero
       end where
       where (prof%shd_LEstoma_new(1:ncl) /= zero)
          prof%shd_rtrans(1:ncl,mc) = prof%shd_trans_rtrans(1:ncl,mc)/prof%shd_LEstoma_new(1:ncl)
       elsewhere
          prof%shd_rtrans(1:ncl,mc) = zero
       end where
    end do

  END SUBROUTINE leaf_wiso


  ! ------------------------------------------------------------------
  SUBROUTINE soil_flux_wiso()
    USE constants,     ONLY: zero, one, e1, TN0, Rw
    USE types,         ONLY: flux, soil, wiso, prof, time
    USE setup,         ONLY: nwiso
    USE utils,         ONLY: es
    USE isotope_utils, ONLY: alpha_kin_h2o, alpha_equ_h2o, isorat
    USE string_utils,  ONLY: num2str
    USE messages,      ONLY: message

    IMPLICIT NONE

    INTEGER(i4) :: mc
    REAL(wp), DIMENSION(1:nwiso) :: litterevap, soilevap, rl, rs, ra, ztmp
    REAL(wp) :: kv_soil, kv_litter
    REAL(wp) :: alpha_equ, alpha_l, alpha_s

    ! water
    soilevap(1)   = flux%soilevap(1)   * soil%latent
    litterevap(1) = flux%litterevap(1) * soil%latent
    ! kvsoil is soil surface conductance to water vapor transfer (s m-1)
    kv_soil = one / (soil%rb+soil%rs)
    ! kvsoil is litter surface conductance to water vapor transfer (s m-1)
    kv_litter = one / (soil%rb+soil%rl)
    ! litter variables
    rl(1:nwiso) = isorat(soil%theta_l(1:nwiso), soil%theta_l(1), wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
    if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('SOIL_FLUX_WISO: ', 'Lost 01 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif
    rs(1:nwiso) = isorat(soil%theta(1,1:nwiso), soil%theta(1,1), wiso%lost(1:nwiso), wiso%lost(1))
    ! litter & soil
#ifdef DEBUG
    if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('SOIL_FLUX_WISO: ', 'Lost 02 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif
    ra(1:nwiso) = isorat(prof%rhov_air_filter(1,1:nwiso), prof%rhov_air_filter_save(1), &
         wiso%lost(1:nwiso), wiso%lost(1))
#ifdef DEBUG
    if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
       call message('SOIL_FLUX_WISO: ', 'Lost 02 @  ', num2str(time%daytime))
       soil%lost0 = 1
    end if
#endif

    do mc=2, nwiso
       alpha_equ = alpha_equ_h2o(soil%tsrf+TN0, mc)
       alpha_l   = alpha_kin_h2o(soil%rl, soil%rb, mc, wiso%merlivat)
       alpha_s   = alpha_kin_h2o(soil%rs, soil%rb, mc, wiso%merlivat)
       ! litter evap
       if (wiso%nofraclitter == 1) then
          if (litterevap(1) >= zero) then
             litterevap(mc) = rl(mc) * litterevap(1)
          else
             litterevap(mc) = ra(mc) * litterevap(1)
          end if
       else
          if (litterevap(1) /= zero) then
             litterevap(mc) = soil%c_litterevap_save * soil%lecoef * alpha_l * kv_litter * &
                  (soil%rh_litter * alpha_equ * rl(mc) * es(soil%tsrf+TN0)*100._wp - &
                  prof%rhov_air_filter(1,mc) * (soil%T_air+TN0) * Rw)
          else
             litterevap(mc) = zero
          end if
       end if
       ! soil evap
       if (wiso%nofracsoil == 1) then
          if (soilevap(1) >= zero) then
             soilevap(mc) = rs(mc) * soilevap(1)
          else
             soilevap(mc) = ra(mc) * soilevap(1)
          end if
       else
          if (soilevap(1) /= zero) then
             soilevap(mc) = soil%lecoef * alpha_s * kv_soil * &
                  (soil%rh_soil * alpha_equ * rs(mc) * es(soil%tsrf+TN0)*100._wp - &
                  prof%rhov_air_filter(1,mc) * (soil%T_air+TN0) * Rw)
          else
             soilevap(mc) = zero
          end if
       end if
    end do
    if ((soil%maxlitter == 1) .and. (nwiso > 1)) then
       litterevap(2:nwiso) = soil%maxlitterevap * rl(2:nwiso)
    end if
    ! Checks
    ztmp(1:nwiso) = one
    if (soilevap(1) /= zero) then
       ztmp(1:nwiso) = (one-(abs(soilevap(1))-abs(soilevap(1:nwiso))/wiso%vsmow(1:nwiso)) / &
            abs(soilevap(1)))*100._wp
    end if
    if (any(ztmp > e1) .and. (abs(soilevap(1)) <= tiny(one))) then
       if (any(abs(soilevap(1:nwiso)) > tiny(one))) then
          call message('SOIL_FLUX_WISO: ', 'bad s1=0 @ ', num2str(time%daytime))
       end if
    end if
    ztmp(1:nwiso) = one
    if (litterevap(1) /= zero) then
       ztmp(1:nwiso) = (one-(abs(litterevap(1))-abs(litterevap(1:nwiso))/wiso%vsmow(1:nwiso)) / &
            abs(litterevap(1)))*100._wp
    end if
    if (any(ztmp > e1) .and. (abs(litterevap(1)) <= tiny(one))) then
       if (any(abs(litterevap(1:nwiso)) > tiny(one))) then
          call message('SOIL_FLUX_WISO: ', 'bad l1=0 @ ', num2str(time%daytime))
       end if
    end if
    ! mm -> W/m2
    flux%soilevap(1:nwiso)   = soilevap(1:nwiso) / soil%latent
    flux%litterevap(1:nwiso) = litterevap(1:nwiso) / soil%latent
    flux%s_evap(1:nwiso)     = flux%litterevap(1:nwiso) + flux%soilevap(1:nwiso)

  END SUBROUTINE soil_flux_wiso


END MODULE isotopes
