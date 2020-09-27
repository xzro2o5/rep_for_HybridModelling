MODULE energy_carbon

  ! This module contains the energy and carbon balance routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bole_respiration         ! computes bole respiration
  PUBLIC :: energy_and_carbon_fluxes ! profiles of energy, stomatal conductance and photosynthesis
  PUBLIC :: energy_balance           ! computes leaf energy balance
  PUBLIC :: photosynthesis           ! computes leaf photosynthesis with Farquhar model
  PUBLIC :: soil_respiration         ! computes soil respiration
  PUBLIC :: srf_vpd                  ! humidity at leaf surface

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE bole_respiration()
    ! Compute bole respiration
    ! Using data of Edwards and Hanson from WBW to
    ! estimate bole respiration on ground area
    ! SAI is 1.5. Rbole is the sum of growth and maintenance
    ! respiration. Growth respiration is
    ! about 30% of maintenance respiration. Q10 was about 2.4
    ! for oak and 1.7 for maple. Geometric mean is 2.02.
    USE constants,  ONLY: half, e3, rugc, TN0, mass_CO2, mass_air
    USE types,      ONLY: soil, input, bole, time, met, prof, iswitch
    USE parameters, ONLY: extra_nate, eabole, pai, wai, R_base1, R_base2!R_base is base stem resp rate, Yuan 2018.02.13
    USE setup,      ONLY: ncl

    IMPLICIT NONE

    REAL(wp) :: tempbole, air_ratio
!    integer(i4) :: j
    tempbole  = half*soil%T_15cm + half*input%ta
    bole%calc = (tempbole-10._wp) / (rugc*283._wp*(tempbole+TN0))
    if (time%days >= time%leafout .and. time%days <= time%leaffallcomplete) then
       bole%factor = R_base1!0.43_wp
    else
       bole%factor = R_base2!0.259_wp
    end if

    ! for Nate McDowell''s Juniper site
    if (extra_nate == 1) bole%factor = bole%factor * half
    ! bole respiration flux density, micromol m-2 s-1
    bole%respiration_mole = bole%factor * exp(eabole*bole%calc)
    ! convert bole%respiration_mole to mg m-2 s-1
    bole%respiration_mg = bole%respiration_mole * mass_CO2 * e3
    ! Divide the bole resp into layers in the stem space and
    ! subtract from layer ps, mg CO2 m-3 s-1
    air_ratio = mass_air/(met%air_density*mass_CO2)
     if (iswitch%wai_new==0) then
!    bole%layer(1:ncl) = prof%dPAIdz(1:ncl) * (bole%respiration_mole / pai)
        else
     bole%layer(1:ncl) = prof%dWAIdz(1:ncl) * (bole%respiration_mole/wai)
     end if


!  		do  j=1,ncl
! 			if (prof%dRESPdz(j)>0) then
!                bole%layer(1:ncl) = prof%dRESPdz*3.25_wp
!            end if

!  		end  do


!    print *, prof%dRESPdz
!    print *, bole%layer
    prof%dboledz(1:ncl) = bole%layer ! add bole resp to the profile Yuan 2018.02.12
 !print *,   "gpp: ", prof%source_co2
 ! source gpp is negetive, bole resp is positive !!!!
    prof%source_co2(1:ncl) = prof%source_co2(1:ncl) + bole%layer(1:ncl)
!print *, "bole: ",bole%layer
  END SUBROUTINE bole_respiration

  ! ------------------------------------------------------------------

  SUBROUTINE energy_and_carbon_fluxes()
    ! The ENERGY_AND_CARBON_FLUXES routine to computes coupled fluxes
    ! of energy, water and CO2 exchange, as well as leaf temperature. Computataions
    ! are performed for each layer in the canopy and on the sunlit and shaded fractions.

    ! Analytical solution for leaf energy balance and leaf temperature is used. The program
    ! is derived from work by Paw U (1986) and was corrected for errors with a re-derivation
    ! of the equations. The quadratic form of the solution is used, rather than the quartic
    ! version that Paw U prefers.

    ! Estimates of leaf temperature are used to adjust kinetic coefficients for enzyme kinetics,
    ! respiration and photosynthesis, as well as the saturation vapor pressure at the leaf surface.

    ! The Analytical solution of the coupled set of equations for photosynthesis and stomatal
    ! conductance by Baldocchi (1994, Tree Physiology) is used. This equation is a solution to
    ! a cubic form of the photosynthesis equation. The photosynthesis algorithms are from the
    ! model of Farquhar. Stomatal conductance is based on the algorithms of Ball-
    ! Berry and Collatz et al., which couple gs to A

    ! Layer 1 is the soil, Layer 40 is top of the canopy
    USE constants,  ONLY: zero, one, TN0, Rw, rugc, mass_CO2
    USE types,      ONLY: prof, solar, fact, bound_lay_res, met, time, debug, iswitch
    USE parameters, ONLY: pai
    USE setup,      ONLY: ncl
    USE utils,      ONLY: es, lambda
    USE transport,  ONLY: boundary_resistance
    USE messages,   ONLY: message
    USE string_utils,  ONLY: num2str

    IMPLICIT NONE

    INTEGER(i4) :: j
    REAL(wp) :: Tair_K_filtered ! temporary absolute air temperature
    REAL(wp) :: T_srf_K, T_srf_C ! surface temperatures in Kelvin and Centigrade
    REAL(wp) :: H_sun, loutsun, Rn_sun, A_sun! energy fluxes on sun leaves
    REAL(wp) :: H_shade, loutsh, Rn_shade, A_shade ! energy fluxes on shaded leaves
    REAL(wp) :: LE_leaf, LE_wet, H_leaf, lout_leaf
    REAL(wp) :: wj_leaf, wc_leaf, wp_leaf, surface_rh, surface_vpd
    REAL(wp) :: rs_sun, rs_shade, A_mg, GPP, resp, internal_CO2, surface_CO2, chloroplast_CO2
    REAL(wp) :: GOP, O_sun, O_shade, A_O2, resp_O2, resp_ROC ! leaf level oxygen flux in photosynthesis and dark respirations
    REAL(wp) :: csca, cica, ccca
    REAL(wp) :: fact_rs_sun, fact_rs_shd

    wj_leaf = zero
    wc_leaf = zero
    wp_leaf = zero

    do j=1, ncl
       ! zero summing values
       H_sun   = zero
       LE_leaf = zero
       LE_wet  = zero
       Rn_sun  = zero
       loutsun = zero
!       call message("a new timestep:", trim(num2str(time%daytime)))
!       call message("a new layer", trim(num2str(j)))
!       call message("load rs_sun from previous hour:",trim(num2str(prof%sun_rs_filter(j))))
    !   print *, "load rs_sun from previous hour:"
    !   print *, prof%sun_rs_filter(j)
       rs_sun  = prof%sun_rs_filter(j)
!       prof%sun_rs_save(j) = rs_sun
!       prof%sun_gs(j)      = one / rs_sun
       A_sun   = zero
       O_sun   = zero
       A_mg    = zero
       A_O2    = zero
       GPP     = zero
       GOP     = zero
       resp    = zero
       resp_O2 = zero
       internal_CO2    = prof%co2_air_filter(j)
       surface_CO2     = prof%co2_air_filter(j)
       chloroplast_CO2 = prof%co2_air_filter(j)
       csca    = one
       cica    = one
       ccca    = one
       surface_rh  = prof%rhov_air_filter(j,1)*(prof%tair_filter(j)+TN0)*Rw / &
            (es(prof%tair_filter(j)+TN0)*100._wp)
       surface_vpd = (one-surface_rh)*es(prof%tair_filter(j)+TN0)*100._wp
       ! First compute energy balance of sunlit leaf, then
       ! repeat and compute algorithm for shaded leaves.
       ! Remember layer is two-sided so include upper and lower streams
       ! are considered.
       ! KC is the convective transfer coeff. (W M-2 K-1). A factor
       ! of 2 is applied since convective heat transfer occurs
       ! on both sides of leaves.
       ! To calculate the total leaf resistance we must combine stomatal
       ! and the leaf boundary layer resistance. Many crops are amphistomatous
       ! so KE must be multiplied by 2. Deciduous forest, on the other hand
       ! has stomata on one side of the leaf.
       Tair_K_filtered = prof%tair_filter(j) + TN0 ! absolute air temperature
       ! Initialize surface temperature with air temperature
       T_srf_K         = Tair_K_filtered
       ! Energy balance on sunlit leaves
       ! update latent heat with new temperature during each call of this routine
       fact%latent     = lambda(Tair_K_filtered)
       ! if (j==1 .or. j==40) print*, 'EC01 ', j, surface_rh, T_srf_K, fact%latent
       if (solar%prob_beam(j) > zero) then
          ! Compute the resistances for heat and vapor transfer, rh and rv,
          ! for each layer, s/m
          call boundary_resistance(prof%ht(j), prof%sun_tleaf_filter(j), prof%cws(j,1), j)
          ! if (j==1 .or. j==40) print*, 'EC02 ', prof%ht(j), prof%sun_tleaf_filter(j), prof%cws(j,1)
          ! compute energy balance of sunlit leaves
!          call message("rs_sun used for LE:",trim(num2str(rs_sun)))
 !print *, "rs_sun used for LE:"
 !print *, rs_sun
          call energy_balance(solar%rnet_sun(j), T_srf_K, Tair_K_filtered, prof%rhov_air_filter(j,1), &
               bound_lay_res%vapor, rs_sun, LE_leaf, LE_wet, H_leaf, lout_leaf, &
               prof%wet_coef_filter(j))
          ! if (j==1 .or. j==40) print*, 'EC03 ', T_srf_K, LE_leaf, LE_wet
          ! if (j==1 .or. j==40) print*, 'EC04 ', H_leaf, lout_leaf
          ! compute photosynthesis of sunlit leaves if leaves have emerged
        !  if (prof%dLAIdz(j) > pai/ncl) then
          if (prof%dLAIdz(j) > zero) then
!            if (j== 38) then
!                print *, prof%dLAIdz(j),solar%quantum_sun(j)
!            end if
! call message("rs_sun before ps:",trim(num2str(rs_sun)))
 !print *, "rs_sun before ps:"
 !print *, rs_sun
 !print *, "sunlit leaves:"
 !print *, j
             call photosynthesis(solar%quantum_sun(j), rs_sun, prof%ht(j), &
                  prof%co2_air_filter(j), T_srf_K, LE_leaf, A_mg, A_O2, GPP, GOP,&
                  resp, resp_O2, resp_ROC, internal_CO2, surface_CO2, chloroplast_CO2, cica, ccca, &
                  surface_rh, surface_vpd, wj_leaf, wc_leaf, wp_leaf, prof%sun_tpu_coeff(j),j)
          end if
          ! Assign values of function to the LE and H source/sink strengths
          T_srf_C = T_srf_K - TN0 ! surface temperature, Centigrade
          H_sun   = H_leaf ! sensible heat flux
          prof%sun_LEstoma(j,1) = LE_leaf    ! latent heat flux from stomata
          prof%sun_LEwet(j,1)   = LE_wet     ! latent heat flux from wet leaf surface
          prof%sun_tleaf(j)     = T_srf_C
          loutsun               = lout_leaf ! long wave out
          Rn_sun                = solar%rnet_sun(j) - lout_leaf ! net radiation
!          if (j==40) then
!          print *, "Rn_leaf,  net_leaf,    longwave_leaf",Rn_sun, solar%rnet_sun(j) , lout_leaf
!          endif
          A_sun                 = A_mg ! leaf photosynthesis, mg CO2 m-2 s-1
          O_sun                 = A_O2
          prof%sun_resp(j)      = resp ! respiration on sun leaves
          prof%sun_resp_O2(j)   = resp_O2
          prof%RQ_sun(j)        = 1/resp_ROC
          prof%ROC_sun(j)       = resp_ROC
          prof%sun_gs(j)        = one / rs_sun ! stomatal conductance
          prof%sun_rs(j)        = rs_sun
          prof%sun_rs_filter(j) = prof%sun_rs(j)! filter should be updated because it is used in every iteration. Yuan 2018.09.13
 !print *, "saved rs_sun for next hour:"
 !print *, prof%sun_rs_filter(j)
!  call message("saved rs_sun for next hour:",trim(num2str(prof%sun_rs_filter(j))))
 !          prof%sun_rs_save(j) = rs_sun !!Yuan test 2018.05.15
          fact_rs_sun           = rugc * (prof%sun_tleaf(j)+TN0)/met%press_Pa !conversion factor
          prof%sun_gs_mol(j)    = one / (prof%sun_rs(j)*fact_rs_sun*1.577_wp) ! convert to mol m-2 s-1 CO2
          prof%sun_wj(j)        = wj_leaf
          prof%sun_wc(j)        = wc_leaf
          prof%sun_wp(j)        = wp_leaf
          prof%sun_GPP(j)       = GPP ! micromolC m-2 s-1
          prof%sun_GOP(j)       = GOP
          prof%sun_A(j)         = A_sun * 1000._wp / mass_CO2 ! micromolC m-2 s-1
          prof%sun_psn_O2(j)    = O_sun
          prof%sun_rbh(j)       = bound_lay_res%heat
          prof%sun_rbv(j)       = bound_lay_res%vapor
          prof%sun_rbco2(j)     = bound_lay_res%co2
          prof%sun_ci(j)        = internal_CO2
          prof%sun_cs(j)        = surface_CO2
          prof%sun_cc(j)        = chloroplast_CO2
          prof%sun_csca(j)      = surface_CO2/prof%co2_air_filter(j)
          prof%sun_cica(j)      = cica
          prof%sun_ccca(j)      = ccca
          prof%sun_rh(j)        = surface_rh ! relative humidity at leaf surface (0 to 1)
          prof%sun_vpd(j)       = surface_vpd ! vapor pressure deficit at leaf surface (hPa)
       end if
       ! Energy balance on shaded leaves
       T_srf_K = Tair_K_filtered
       ! initial value of stomatal resistance based on light
       rs_shade = prof%shd_rs_filter(j) ! should switch off as v3.4. Yuan 2018.06.04
!       prof%shd_rs_save(j) = rs_shade ! stomatal resistance, shaded leaf ! should switch off as v3.4. Yuan 2018.06.04
       ! boundary layer resistances on shaded leaves. With different
       ! surface temperature, the convective effect may differ from that
       ! computed on sunlit leaves
       call boundary_resistance(prof%ht(j), prof%shd_tleaf_filter(j), prof%cws(j,1), j)
       ! if (j==1 .or. j==40) print*, 'EC05 ', prof%ht(j), prof%shd_tleaf_filter(j), prof%cws(j,1)
       ! Energy balance of shaded leaves
       call energy_balance(solar%rnet_shd(j), T_srf_K, Tair_K_filtered, prof%rhov_air_filter(j,1), &
            bound_lay_res%vapor, rs_shade, LE_leaf, LE_wet, H_leaf, lout_leaf, &
            prof%wet_coef_filter(j))
       ! if (j==1 .or. j==40) print*, 'EC06.0 ', solar%rnet_shd(j), Tair_K_filtered, prof%rhov_air_filter(j,1)
       ! if (j==1 .or. j==40) print*, 'EC06.1 ', bound_lay_res%vapor, rs_shade, prof%wet_coef_filter(j)
       ! if (j==1 .or. j==40) print*, 'EC06 ', T_srf_K, LE_leaf, LE_wet
       ! if (j==1 .or. j==40) print*, 'EC07 ', H_leaf, lout_leaf

       ! compute photosynthesis and stomatal conductance of shaded leaves
   !    if (prof%dLAIdz(j) > pai/ncl) then
        if (prof%dLAIdz(j) > zero) then
!            print *, "shaded leaves:"
!            print *, j
          call photosynthesis(solar%quantum_shd(j), rs_shade,prof%ht(j), prof%co2_air_filter(j), &
               T_srf_K, LE_leaf, A_mg, A_O2, GPP, GOP, resp, resp_O2, resp_ROC, internal_CO2, surface_CO2, &
               chloroplast_CO2, cica, ccca, surface_rh, surface_vpd, wj_leaf, wc_leaf, wp_leaf, prof%shd_tpu_coeff(j),j)
      ! call energy_balance(solar%rnet_shd(j), T_srf_K, Tair_K_filtered, prof%rhov_air_filter(j,1), &
      !      bound_lay_res%vapor, rs_shade, LE_leaf, LE_wet, H_leaf, lout_leaf, &
      !      prof%wet_coef_filter(j))
       end if
       ! re-assign variable names from functions output
       T_srf_C               = T_srf_K - TN0 ! surface temperature, C
       prof%shd_LEstoma(j,1) = LE_leaf    ! latent heat flux from stomata
       prof%shd_LEwet(j,1)   = LE_wet     ! latent heat flux from wet leaf surface
       H_shade               = H_leaf ! sensible heat flux density, shaded leaves, W m-2
       loutsh                = lout_leaf ! long wave energy emissive flux density, shaded leaves, W m-2
       Rn_shade              = solar%rnet_shd(j) - lout_leaf ! net radiation balance, shaded leaves, W m-2
       prof%shd_wj(j)        = wj_leaf ! electron transport velocity, shaded leaves, micromol m-2 s-1
       prof%shd_wc(j)        = wc_leaf ! carboxylation velocity, shaded leaves, micromol m-2 s-1
       prof%shd_wp(j)        = wp_leaf
       A_shade               = A_mg ! photosynthesis, shaded leaves, mgCO2 m-2 s-1
       O_shade               = A_O2
       prof%shd_rh(j)        = surface_rh ! relative humidity at leaf surface (0 to 1)
       prof%shd_vpd(j)       = surface_vpd ! vapor pressure deficit at leaf surface (hPa)
       ! compute profiles
       prof%shd_resp(j)   = resp
       prof%shd_resp_O2(j)   = resp_O2
       prof%RQ_shd(j)        = 1/resp_ROC
       prof%ROC_shd(j)       = resp_ROC
       prof%shd_tleaf(j)  = T_srf_C
       prof%shd_gs(j)     = one/rs_shade ! stomatal conductance, shaded leaf
       prof%shd_rs(j)     = rs_shade ! stomatal resistance, shaded leaf
       prof%shd_rs_filter(j)=prof%shd_rs(j)! Yuan updated 2018.09.13
       prof%shd_rs_save(j) = rs_shade !!Yuan test 2018.05.15
       fact_rs_shd        = rugc * (prof%shd_tleaf(j)+TN0)/ met%press_Pa !conversion factor
       prof%shd_gs_mol(j) = one/(prof%shd_rs(j)*fact_rs_shd*1.577_wp) ! convert to mol m-2 s-1 CO2
       prof%shd_GPP(j)    = GPP ! micromolC m-2 s-1
       prof%shd_GOP(j)    = GOP
       prof%shd_A(j)      = A_shade*1000._wp/mass_CO2 ! micromolC m-2 s-1
       prof%shd_psn_O2(j) = O_shade! net O2 flux in photosynthesis
       prof%shd_rbh(j)    = bound_lay_res%heat
       prof%shd_rbv(j)    = bound_lay_res%vapor
       prof%shd_rbco2(j)  = bound_lay_res%co2
       prof%shd_ci(j)     = internal_CO2
       prof%shd_cs(j)     = surface_CO2
       prof%shd_cc(j)     = chloroplast_CO2
       prof%shd_csca(j)   = surface_CO2/prof%co2_air_filter(j)
       prof%shd_cica(j)   = cica
       prof%shd_ccca(j)   = ccca
       ! compute layer energy fluxes, weighted by leaf area and sun and shaded fractions
       ! prof%dLEdz_sun(j) = prof%dPAIdz(j) * solar%prob_beam(j) * &
       !     (prof%sun_LEstoma(j,1)+prof%sun_LEwet(j,1))
       ! scale transpiration with LAI but wet leaf evaperation with WAI, Yuan 2018.08.28
      ! prof%dLEdz_sun(j) = prof%dLAIdz(j) * solar%prob_beam(j) * &
      !      prof%sun_LEstoma(j,1)+ prof%dWAIdz(j) * solar%prob_beam(j)* prof%sun_LEwet(j,1)
        prof%dLEdz_sun(j) = prof%dLAIdz(j) * solar%prob_beam(j) * prof%sun_LEstoma(j,1)+ &
                            prof%dLAIdz(j) * solar%prob_beam(j)* prof%sun_LEwet(j,1)

       !prof%dLEdz_shd(j) = prof%dPAIdz(j) * solar%prob_shd(j) * &
       !     (prof%shd_LEstoma(j,1)+prof%shd_LEwet(j,1))
  !     prof%dLEdz_shd(j) = prof%dLAIdz(j) * solar%prob_shd(j) * &
  !          prof%shd_LEstoma(j,1)+ prof%dWAIdz(j) * solar%prob_shd(j)* prof%shd_LEwet(j,1)
       prof%dLEdz_shd(j) = prof%dLAIdz(j) * solar%prob_shd(j) * prof%shd_LEstoma(j,1)+ &
                           prof%dLAIdz(j) * solar%prob_shd(j)* prof%shd_LEwet(j,1)
       prof%dLEdz(j,1)   = prof%dLEdz_sun(j) + prof%dLEdz_shd(j)
               if ((prof%dLEdz(j,1)==0)) then
                !debug%R1=rvsrf
                !debug%R2=stomsrf
                !debug%R3=tsrfkpt
                !debug%R4=es(tsrfkpt)*100._wp-ea
!                print *, prof%dLEdz_sun(j)
!                print *, prof%dLEdz_shd(j)
!                print *, prof%dLAIdz(j)
!                print *, solar%prob_beam(j)
!                print *, prof%sun_LEstoma(j,1)
!                print *, prof%sun_LEwet(j,1)

               end if
               ! H, radiation should be scaled with PAI
               if (iswitch%wai_new==0) then
                    prof%dHdz(j)      = prof%dLAIdz(j) * (solar%prob_beam(j) * H_sun + solar%prob_shd(j) * H_shade)
                    prof%dRNdz(j)     = prof%dLAIdz(j) * (solar%prob_beam(j) * Rn_sun + solar%prob_shd(j) * Rn_shade)
                    prof%dLoutdz(j)   = prof%dLAIdz(j) * (solar%prob_beam(j) * loutsun + solar%prob_shd(j) * loutsh)
                else
                    prof%dHdz(j)      = prof%dPAIdz(j) * (solar%prob_beam(j) * H_sun + solar%prob_shd(j) * H_shade)
                    prof%dRNdz(j)     = prof%dPAIdz(j) * (solar%prob_beam(j) * Rn_sun + solar%prob_shd(j) * Rn_shade)
                    prof%dLoutdz(j)   = prof%dPAIdz(j) * (solar%prob_beam(j) * loutsun + solar%prob_shd(j) * loutsh)


               end if
 !      prof%dHdz(j)      = prof%dLAIdz(j) * (solar%prob_beam(j) * H_sun + solar%prob_shd(j) * H_shade)
       !print *, prof%dHdz(j)
 !      if (ISNAN(prof%dHdz(j)) ) then
 !       print *, "H=NAN"
 !       print *, prof%dLAIdz(j), prof%dHdz(j), solar%prob_beam(j), solar%prob_shd(j), H_sun, H_shade
 !      end if
 !      prof%dRNdz(j)     = prof%dLAIdz(j) * (solar%prob_beam(j) * Rn_sun + solar%prob_shd(j) * Rn_shade)
       !print *, prof%dPAIdz(j), prof%dRNdz(j)
       !print *, prof%dLAIdz(j), solar%prob_beam(j), Rn_sun, solar%prob_shd(j) , Rn_shade

 !      prof%dLoutdz(j)   = prof%dLAIdz(j) * (solar%prob_beam(j) * loutsun + solar%prob_shd(j) * loutsh)
!print *, prof%dLoutdz(j)
       ! photosynthesis of the layer, prof%dPsdz has units mg m-3 s-1
       !prof%dPsdz(j) = prof%dLAIdz(j) * (A_sun * solar%prob_beam(j) + A_shade * solar%prob_shd(j))
       ! photosynthesis of layer, prof%dPsdz has units of micromoles m-2 s-1
       prof%dPsdz_sun(j) = prof%dLAIdz(j) * prof%sun_A(j) * solar%prob_beam(j)
       prof%dPsdz_shd(j) = prof%dLAIdz(j) * prof%shd_A(j) * solar%prob_shd(j)
       prof%dPsdz(j)     = prof%dPsdz_sun(j) + prof%dPsdz_shd(j)
       ! photosynthetic O2 per ground level:
       prof%dPsdz_O2_sun(j) = prof%dLAIdz(j) * prof%sun_psn_O2(j) * solar%prob_beam(j)
       prof%dPsdz_O2_shd(j) = prof%dLAIdz(j) * prof%shd_psn_O2(j) * solar%prob_shd(j)
       prof%dPsdz_O2(j)     = prof%dPsdz_O2_sun(j) + prof%dPsdz_O2_shd(j)
       ! GPP of layer, prof%dGPPdz has units of micromoles m-2 s-1
       prof%dGPPdz_sun(j) = prof%dLAIdz(j) * prof%sun_GPP(j) * solar%prob_beam(j)
       prof%dGPPdz_shd(j) = prof%dLAIdz(j) * prof%shd_GPP(j) * solar%prob_shd(j)
       prof%dGPPdz(j)     = prof%dGPPdz_sun(j) + prof%dGPPdz_shd(j)
       ! respiration of the layer, micromol m-2 s-1
       prof%dRESPdz_sun(j) = prof%dLAIdz(j) * prof%sun_resp(j) * solar%prob_beam(j)
       prof%dRESPdz_shd(j) = prof%dLAIdz(j) * prof%shd_resp(j) * solar%prob_shd(j)
       prof%dRESPdz(j)     = prof%dRESPdz_sun(j) + prof%dRESPdz_shd(j)
       ! O2 in dark respiration:
       prof%dRESPdz_O2_sun(j) = prof%dLAIdz(j) * prof%sun_resp_O2(j) * solar%prob_beam(j)
       prof%dRESPdz_O2_shd(j) = prof%dLAIdz(j) * prof%shd_resp_O2(j) * solar%prob_shd(j)
       prof%dRESPdz_O2(j)     = prof%dRESPdz_O2_sun(j) + prof%dRESPdz_O2_shd(j)
 !      print *, 'leaf resp:'
 !      print *, prof%dRESPdz_sun(j)
 !      print *, prof%dRESPdz_shd(j)
       ! prof%sun_gs has units of: m s-1
       ! prof%dStomCondz has units of: mol m-2 s-1
       prof%dStomCondz_sun(j) = prof%dLAIdz(j) * solar%prob_beam(j)*prof%sun_gs(j)
       prof%dStomCondz_shd(j) = prof%dLAIdz(j) * solar%prob_shd(j)*prof%shd_gs(j)
       prof%dStomCondz(j)     = prof%dStomCondz_sun(j) + prof%dStomCondz_shd(j)
       if ((prof%sun_LEwet(j,1)*solar%prob_beam(j) + prof%shd_LEwet(j,1)*solar%prob_shd(j)) /= zero) then
          prof%wet_coef(j) = max(min(prof%cws(j,1) * fact%latent / (time%time_step * prof%dLAIdz(j) * &
               (prof%sun_LEwet(j,1)*solar%prob_beam(j) + prof%shd_LEwet(j,1)*solar%prob_shd(j))) , one), zero)
        !   prof%wet_coef(j) = max(min((time%time_step * prof%dLAIdz(j) * &
        !       (prof%sun_LEwet(j,1)*solar%prob_beam(j) + prof%shd_LEwet(j,1)*solar%prob_shd(j))/prof%cws(j,1) * fact%latent) , &
        !       one), zero)
!               print *, prof%wet_coef(j)

       else
          prof%wet_coef(j) = zero
       end if
    end do

  END SUBROUTINE energy_and_carbon_fluxes


  ! ------------------------------------------------------------------
  SUBROUTINE energy_balance(qrad, tsrfkpt, taa, rhovva, rvsrf, stomsrf, &
       lept, lewet, H_leafpt, lout_leafpt, wet_coef)
    ! ENERGY BALANCE COMPUTATION
    ! A revised version of the quadratic solution to the leaf energy balance relationship is used.
    ! Paw U, KT. 1987. J. Thermal Biology. 3: 227-233
    ! H is sensible heat flux density on the basis of both sides of a leaf
    ! J m-2 s-1 (W m-2). Note KC includes a factor of 2 here for heat flux
    ! because it occurs from both sides of a leaf.
    ! LE is latent heat flux density through stomata. For hypostomatous leaves (n_stomata_sides = 1) LE
    ! occurs on only one side. For amphistomatous leaves (n_stomata_sides = 2) LE occurs
    ! on both sides.
    ! LEw is latent heat flux density from wet leave surfaces. It occurs on both sides.
    USE constants,  ONLY: zero, half, one, two, cp, Rw
    USE types,      ONLY: fact, met, bound_lay_res, iswitch, debug
    USE parameters, ONLY: epsigma2, epsigma8, epsigma12, n_stomata_sides
    USE utils,      ONLY: es, desdt, des2dt

    IMPLICIT NONE

    REAL(wp), INTENT(IN)  :: qrad
    REAL(wp), INTENT(OUT) :: tsrfkpt
    REAL(wp), INTENT(IN)  :: taa
    REAL(wp), INTENT(IN)  :: rhovva
    REAL(wp), INTENT(IN)  :: rvsrf
    REAL(wp), INTENT(IN)  :: stomsrf
    REAL(wp), INTENT(OUT) :: lept
    REAL(wp), INTENT(OUT) :: lewet
    REAL(wp), INTENT(OUT) :: H_leafpt
    REAL(wp), INTENT(OUT) :: lout_leafpt
    REAL(wp), INTENT(IN)  :: wet_coef

    REAL(wp) :: est, ea, tkta
    REAL(wp) :: tk2, tk3, tk4
    REAL(wp) :: dest, d2est
    REAL(wp) :: lecoef, hcoef, hcoef2, product
    !REAL(wp) :: bcoef, ccoef, repeat,acoef, acoeff, le2
    REAL(wp) :: atlf, btlf, ctlf,vpd_leaf,llout
    REAL(wp) :: ke

    tkta     = taa               ![K]
    est      = es(tkta) * 100._wp   ! converts es(T) from mb to Pa
    ! ea = RHOA * TAA * 1000 / 2.165
    ea       = rhovva * taa * Rw ! vapor pressure above leaf
    ! Vapor pressure deficit, Pa
    vpd_leaf = est - ea
    ! Slope of the vapor pressure-temperature curve, Pa/C
    ! evaluate as function of Tk
    dest     = desdt(tkta, fact%latent)
    ! print*, 'EB01 ', tkta, fact%latent, dest
    ! Second derivative of the vapor pressure-temperature curve, Pa/C
    ! Evaluate as function of Tk
    d2est    = des2dt(tkta)
    ! print*, 'EB02 ', d2est
    ! Compute products of air temperature, K
    tk2      = tkta * tkta
    tk3      = tk2 * tkta
    tk4      = tk3 * tkta
    ! Longwave emission at air temperature, W m-2
    llout    = epsigma2 * tk4
    ! print*, 'EB03 ', epsigma2, tk4, llout
    ! Check if leaves have stomata on one or both sides
    ! Cuticle resistance is included in STOM.
    ! hypostomatous  n_stomata_sides = 1
    ! amphistomatous  n_stomata_sides = 2
    ! read in from parameter file in set_parameter()
    ! Interception led to wet leaf surface (both sides of leaves), which can contribute
    ! to evaporation from leaf surfaces. Therefore we have to include an evaporaton term
    ! without stomata resistance (only boundary layer resistance) We introduce
    ! wet_coef[0-1] to make sure that LE from wet surface does not exceed the amount
    ! of water on the leaves.
    !  Coefficient for latent heat flux
    ke     = one/ (rvsrf + stomsrf) + wet_coef/rvsrf
!    print *,"wet_coef",wet_coef
 !   print *, "rvsrf + stomsrf:         ", rvsrf , stomsrf
    lecoef = met%air_density * 0.622_wp * fact%latent * ke / met%press_Pa
 !   print *, "air_density, fact%latent,ke,press"
 !   print *, met%air_density, fact%latent,ke,met%press_Pa
    ! Coefficients for sensible heat flux
    hcoef  = met%air_density*cp/bound_lay_res%heat
!    print *, cp
    hcoef2 = two * hcoef
    ! now LE is not directly calculated with the quadratic solution anymore,
    ! but we first solve for leaf temperature with a quadratic equation. Then we use
    ! the outcome of leaf temperature to calculate LE
    ! ! The quadratic coefficients for the a LE^2 + b LE +c =0
    ! repeat = hcoef + epsigma4 * tk3
    ! acoeff = lecoef * d2est / (2. * repeat)
    ! acoef = acoeff / 4.
    ! bcoef = -(repeat) - lecoef * dest / 2. + acoeff * (-qrad / 2. + llout)
    ! ccoef = repeat * lecoef * vpd_leaf + lecoef * dest * (qrad / 2. - llout) + acoeff * ((qrad * qrad) / 4. + llout * llout - qrad * llout)
    ! ! LE1 = (-BCOEF + (BCOEF ^ 2 - 4 * ACOEF * CCOEF) ^ .5) / (2 * ACOEF)
    ! product = bcoef * bcoef - 4. * acoef * ccoef
    ! ! LE2 = (-BCOEF - (BCOEF * BCOEF - 4 * acoef * CCOEF) ^ .5) / (2. * acoef)
    ! le2= (-bcoef - sqrt(product)) / (2. * acoef)
    ! *lept=le2 ! need to pass pointer out of subroutine
    ! solve for Ts using quadratic solution
    ! coefficients to the quadratic solution
    atlf    = epsigma12 * tk2 + n_stomata_sides * d2est * lecoef * half
    btlf    = epsigma8 * tk3 + hcoef2 + n_stomata_sides * lecoef * dest
    ctlf    = -qrad + llout + n_stomata_sides * lecoef * vpd_leaf
    product = btlf * btlf - 4._wp * atlf * ctlf
    ! print*, 'EB04 ', lecoef, hcoef
   ! print *, "c=:" ,ctlf
!    print *, "lecoef , vpd_leaf:  ", lecoef , vpd_leaf
!    print *, "-b:      ", -btlf,"sqrt:   ",sqrt(product)
    if (product >= zero) then
       tsrfkpt = tkta + (-btlf + sqrt(product)) / (two * atlf) ! [K]
    else
       tsrfkpt = tkta ! [K]
    end if
    if (tsrfkpt < 230._wp .or. tsrfkpt > 335._wp) then
       tsrfkpt = tkta ! [K]
    end if
    ! long wave emission of energy
    !*lout_leafpt = epsigma2 * pow(*tsrfkpt, 4)
    lout_leafpt = llout + epsigma8 * tkta*tkta*tkta * (tsrfkpt-tkta) + &
         epsigma12 * tkta*tkta * (tsrfkpt-tkta)*(tsrfkpt-tkta)
    ! H is sensible heat flux
    H_leafpt    = hcoef2 * (tsrfkpt-tkta)
    !print *, "H_coef,    tleaf-tair,      sensible heat"
    !print *, hcoef2,(tsrfkpt-tkta),H_leafpt
    ! lept is latent heat flux through stomata
    ! ToDo for isotopes ! transpiration from second Taylor expansion
    !*lept = n_stomata_sides * met%air_density * 0.622 * fact%latent
    !  / (met%press_Pa * (rvsrf + stomsrf))
    !  * (vpd_leaf + dest*(*tsrfkpt-tkta) + d2est/2.*(*tsrfkpt-tkta)*(*tsrfkpt-tkta))
    lept = n_stomata_sides * met%air_density * 0.622_wp * fact%latent / &
         (met%press_Pa * (rvsrf + stomsrf)) * (es(tsrfkpt)*100._wp-ea)
         ! add debug output:
debug%R1=rvsrf
debug%R2=stomsrf
debug%R3=tsrfkpt
debug%R4=es(tsrfkpt)*100._wp-ea
    if (iswitch%no_neg_water_flux == 1) lept = max(lept, zero)
    ! lewet is latent heat flux from wet leaves
    ! ToDo for isotopes ! soil evaporation from second Taylor expansion
    !*lewet = wet_coef * n_stomata_sides*met%air_density * 0.622 * fact%latent
    !  / (met%press_Pa * rvsrf) * (vpd_leaf + dest*(*tsrfkpt-tkta) + d2est/2.*(*tsrfkpt-tkta)*(*tsrfkpt-tkta))
    lewet = wet_coef * n_stomata_sides*met%air_density * 0.622_wp * fact%latent / &
         (met%press_Pa * rvsrf) * (es(tsrfkpt)*100._wp-ea)
    lewet = max(lewet, zero)

  END SUBROUTINE energy_balance


  ! ------------------------------------------------------------------
  SUBROUTINE photosynthesis(Iphoton, rstompt, zzz, cca, tlk, leleaf, A_mgpt, O_pt, &
       GPPpt, GOPpt, resppt, resOppt, ROC_rd, cipnt, cspnt, ccpnt, cicapnt, cccapnt, rh_leafpnt, vpd_leafpnt, &
       wjpnt, wcpnt, wppnt, tpupnt, JJ)
    ! This program solves a cubic equation to calculate
    ! leaf photosynthesis. This cubic expression is derived from solving
    ! five simultaneous equations for A, PG, cs, CI and GS.
    ! Stomatal conductance is computed with the Ball-Berry model.
    ! The cubic derivation assumes that b'', the intercept of the Ball-Berry
    ! stomatal conductance model, is non-zero.

    ! Gs = k A rh/cs + b''

    ! We also found that the solution for A can be obtained by a quadratic equation
    ! when Gs is constant or b'' is zero.

    ! The derivation is published in:

    ! Baldocchi, D.D. 1994. An analytical solution for coupled leaf photosynthesis
    ! and stomatal conductance models. Tree Physiology 14: 1069-1079.

    ! -----------------------------------------------------------------------

    ! A Biochemical Model of C3 Photosynthesis

    ! After Farquhar, von Caemmerer and Berry (1980) Planta.
    ! 149: 78-90.

    ! The original program was modified to incorporate functions and parameters
    ! derived from gas exchange experiments of Harley, who paramertized Vc and J in
    ! terms of optimal temperature, rather than some reference temperature, eg 25C.

    ! Program calculates leaf photosynthesis from biochemical parameters

    ! rd25 - Dark respiration at 25 degrees C (umol m-2 s-1)
    ! tlk - leaf temperature, Kelvin
    ! jmax - optimal rate of electron transport
    ! vcopt - maximum rate of RuBP Carboxylase/oxygenase
    ! iphoton - incident photosynthetically active photon flux (mmols m-2 s-1)

    ! note: Harley parameterized the model on the basis of incident PAR

    ! gs - stomatal conductance (mols m-2 s-1), typically 0.01-0.20
    ! pstat-station pressure, bars
    ! aphoto - net photosynthesis (umol m-2 s-1)
    ! ps - gross photosynthesis (umol m-2 s-1)
    ! aps - net photosynthesis (mg m-2 s-1)
    ! aphoto (umol m-2 s-1)

    ! --------------------------------------------------

    ! iphoton is radiation incident on leaves

    ! The temperature dependency of the kinetic properties of
    ! RUBISCO are compensated for using the Arrhenius and
    ! Boltzmann equations. From biochemistry, one observes that
    ! at moderate temperatures enzyme kinetic rates increase
    ! with temperature. At extreme temperatures enzyme
    ! denaturization occurs and rates must decrease.

    ! Arrhenius Eq.

    ! f(T)=f(tk_25) exp(tk -298)eact/(298 R tk)), where eact is the
    ! activation energy.

    ! Boltzmann distribution

    ! F(T)=tboltz()

    ! Define terms for calculation of gross photosynthesis, PG

    ! PG is a function of the minimum of RuBP saturated rate of
    ! carboxylation, Wc, and the RuBP limited rate of carboxylation, Wj.
    ! Wj is limiting when light is low and electron transport, which
    ! re-generates RuBP, is limiting. Wc is limiting when plenty of RuBP is
    ! available compared to the CO2 that is needed for carboxylation.

    ! Both equations take the form:

    ! PG-photorespiration= (a CI-a d)/(e CI + b)

    ! PG-photorespiration=min(Wj,Wc) (1-gamma/Ci)

    ! Wc=Vcmax Ci/(Ci + Kc(1+O2/Ko))

    ! Wj=J Ci/(4 Ci + 8 gamma)

    ! Ps kinetic coefficients from Harley at WBW.

    ! Gamma is the CO2 compensation point

    ! Jan 14, 1999 Updated the cubic solutions for photosynthesis. There are
    ! times when the restriction that R^2 < Q^3 is violated. I therefore need
    ! alternative algorithms to solve for the correct root.
    USE constants,    ONLY: zero, half, one, two, e3, TN0, rugc, tk_25, pi2, mass_co2
    USE parameters,   ONLY: kc25, ko25, tau25, o2, extra_nate, hkin, skin, &
         ekc, eko, ektau, jmopt, vcopt, htFrac, zh65, lai, evc, &
         toptvc, rd_vc, ejm, toptjm, kball, g0, a1, erd, &
         D0, gm_vc, qalpha, curvature, bprime, g0_mly_in, g1_mly_in, &
         tp_vc, alpha_g_max, alpha_s_max, n_max
    USE types,        ONLY: time, prof, met, bound_lay_res, srf_res, soil, &
         iswitch, output, input
    USE utils,        ONLY: temp_func, tboltz, es
    USE messages,     ONLY: message
    USE string_utils, ONLY: num2str

    IMPLICIT NONE

    REAL(wp),    INTENT(INOUT) :: Iphoton ! Hallo
    REAL(wp),    INTENT(INOUT) :: rstompt
    REAL(wp),    INTENT(IN) :: zzz
    REAL(wp),    INTENT(IN) :: cca
    REAL(wp),    INTENT(IN) :: tlk
    REAL(wp),    INTENT(IN) :: leleaf
    REAL(wp),    INTENT(OUT) :: A_mgpt
    REAL(wp),    INTENT(OUT) :: O_pt ! net photosynthetic O2 in umol m-2 leaf area s-2
    REAL(wp),    INTENT(OUT) :: GPPpt
    REAL(wp),    INTENT(OUT) :: GOPpt
    REAL(wp),    INTENT(OUT) :: resppt
    REAL(wp),    INTENT(OUT) :: resOppt ! leaf dark respiration!
    REAL(wp),    INTENT(OUT) :: ROC_rd ! ROC of leaf dark respiration!
    REAL(wp),    INTENT(OUT) :: cipnt
    REAL(wp),    INTENT(OUT) :: cspnt
    REAL(wp),    INTENT(OUT) :: ccpnt
    REAL(wp),    INTENT(OUT) :: cicapnt
    REAL(wp),    INTENT(OUT) :: cccapnt
    REAL(wp),    INTENT(OUT) :: rh_leafpnt
    REAL(wp),    INTENT(OUT) :: vpd_leafpnt
    REAL(wp),    INTENT(OUT) :: wjpnt
    REAL(wp),    INTENT(OUT) :: wcpnt
    REAL(wp),    INTENT(OUT) :: wppnt
    REAL(wp),    INTENT(INOUT) :: tpupnt
    INTEGER(i4), INTENT(IN) :: JJ

    REAL(wp) :: tprime25, bc, ttemp, gammac
    REAL(wp) :: jmax, vcmax, jmaxz, vcmaxz, cs, ci, cc, vc25z
    REAL(wp) :: kct, ko, tau
    REAL(wp) :: rd, rdz, rd_O2, RQ_rd!O2 in dark respiration
    REAL(wp) :: rb_mole, gb_mole, dd, b8_dd
    REAL(wp) :: rh_leaf, vpd_leaf, es_leaf, k_rh, ci_guess
    REAL(wp) :: j_photon, alpha_ps, bbeta, ggamma
    REAL(wp) :: denom, Pcube, Qcube, Rcube
    REAL(wp) :: P2, P3, Q, R
    REAL(wp) :: root1, root2
    REAL(wp) :: root3, arg_U, ang_L
    REAL(wp) :: aphoto, Ophoto, gpp, gpp_o2, j_sucrose, wj, wj2 ! Ophoto:net phptosynthetic O2
    REAL(wp) :: phi, vo, alpha_g, alpha_g_guess, alpha_s, g_tmp, s_tmp, beta_tpu, tpu_coeff ! add TPU limits to photosynthesis. Yuan 2019.12.20
    REAL(wp) :: tp, wp_tpu
    REAL(wp) :: gs_leaf_mole, gs_co2, gs_m_s
    REAL(wp) :: ps_1,delta_1, Aquad1, Bquad1, Cquad1
    REAL(wp) :: theta_ps, wc, b_ps, a_ps, e_ps, psguess
    REAL(wp) :: sqrprod, product ! delday
    REAL(wp) :: rt
    REAL(wp) :: gm
    REAL(wp) :: beta_ps, gamma_ps, delta_ps, epsilon_ps, zeta_ps, eta_ps
    REAL(wp) :: g1_local !alpha_ps1, alpha_ps2,  ! Medlyn's stomatal model Yuan 2020.02.27
    REAL(wp) :: theta_soil
    REAL(wp) :: g0_local, a1_local, D0_local

    ! double a_cubic, b_cubic, rootprod
    REAL(wp) :: rr, qqq, minroot, maxroot, midroot
    REAL(wp) :: bprime_local, bprime16_local
    INTEGER(i4) :: quad

    jmaxz        = zero
    vcmaxz       = zero
    alpha_ps     = zero
    bbeta        = zero
    ggamma       = zero
    Pcube        = zero
    Qcube        = zero
    Rcube        = zero
    aphoto       = zero
    Ophoto       = zero
    RQ_rd        = zero
    gpp          = zero
    gs_leaf_mole = zero
    theta_ps     = zero
    sqrprod      = zero
    beta_ps      = zero
    gamma_ps     = zero
    delta_ps     = zero
    epsilon_ps   = zero
    zeta_ps      = zero
    eta_ps       = zero
    g0_local     = zero
    a1_local     = zero
    D0_local     = zero
    minroot      = zero
    maxroot      = zero
    midroot      = zero
    bprime_local = zero
    quad         = 0
    cs           = one
    ci           = one
    cc           = one
    tpu_coeff    = zero
   ! alpha_ps1    = zero
   ! alpha_ps2    = zero
    g1_local     = zero

    rt       = rugc * tlk ! product of universal gas constant and abs temperature
    tprime25 = tlk - tk_25 ! temperature difference
    ttemp    = exp((skin * tlk - hkin) / rt) + one ! denominator term
    if (Iphoton < one) Iphoton = zero
    ! KC and KO are solely a function of the Arrhenius Eq.
    kct = temp_func(kc25, ekc, tprime25, tk_25, tlk) ! mubar
    ko  = temp_func(ko25, eko, tprime25, tk_25, tlk) ! mbar
    tau = temp_func(tau25, ektau, tprime25, tk_25, tlk) ! dimensonless
    ! leaf dark RQ as a function of leaf temperature
    RQ_rd = -0.0147*(tlk-TN0)+1.24 !RQ=CO2/O2
    bc  = kct * (one + o2 / ko) ! mubar*(1+mbar/mbar) = mubar
    ! gammac is the CO2 compensation point due to photorespiration, umol mol-1
    ! Recalculate gammac with the new temperature dependent KO and KC
    ! coefficients
    ! gammac = .5 * O2*1000/TAU
    gammac = 500.0_wp * o2 / tau !0.5*mmol mol-1 * 100/tau = mumol mol-1 = ppm
    gammac = 500.0_wp * input%o2air / tau ! use dynamic atom o2 instead of fixed 210000. Yuan 2018.01.31

    ! temperature corrections for Jmax and Vcmax
    ! Scale jmopt and VCOPT with a surrogate for leaf nitrogen
    ! specific leaf weight (Gutschick and Weigel).
    ! normalized leaf wt is 1 at top of canopy and is 0.35
    ! at forest floor. Leaf weight scales linearly with height
    ! and so does jmopt and vcmax
    ! zoverh=0.65/HT=zh65
    if (extra_nate == 1) then
       ! for Nate McDowell''s juniper site, no scaling of Vcmax with height and LAI
       jmaxz  = jmopt(JJ)
       vcmaxz = vcopt(JJ)
    else
       ! time before leaf out
       if (time%days < time%leafout) then
          jmaxz  = zero
          vcmaxz = zero
       end if
       ! spring, increase Ps capacity with leaf expansion as a function of leaf area changes
       if (time%days >= time%leafout .and. time%days < time%leaffull) then
          jmaxz  = jmopt(JJ) * (zh65 * zzz + (1-htFrac)) * time%lai/lai! use htFrac in parameter file instead of a fixed value Yuan 2018.02.21
          vcmaxz = vcopt(JJ) * (zh65 * zzz + (1-htFrac)) * time%lai/lai
!          print *, htFrac, (1-htFrac)
       end if
       ! growing season, full Ps capacity (note newer data by Wilson et al shows more
       ! dynamics
       if (time%days >= time%leaffull .and. time%days < time%leaffall) then
          jmaxz  = jmopt(JJ) * (zh65 * zzz + (1-htFrac))
          vcmaxz = vcopt(JJ) * (zh65 * zzz + (1-htFrac))
       end if
       ! gradual decline in fall
       if (time%days >= time%leaffall .and. time%days <= time%leaffallcomplete) then
          !delday=1-(time%days-270)/30
          jmaxz  = jmopt(JJ) * (zh65 * zzz + (1-htFrac)) * time%lai/lai
          vcmaxz = vcopt(JJ) * (zh65 * zzz + (1-htFrac)) * time%lai/lai
       end if
       if (time%days > time%leaffallcomplete) then
          jmaxz  = zero
          vcmaxz = zero
       end if
    end if
    prof%vcmaxz(JJ) = vcmaxz
    ! Scale rd with height via vcmax and apply temperature
    ! correction for dark respiration
    ! get vcmax at 25 deg C and then take fraction of it as dark respiration (rd_vc)
    vc25z = tboltz(vcmaxz, evc, toptvc, TN0+25._wp, hkin)
    rdz   = vc25z * rd_vc(JJ)
    !? from Harley 1995, sun leaves: Rd(25 deg C)/Vcmax(Topt)=0.34/73=0.004657?
    ! but if we use Rd(25 deg C)/Vcmax(25 deg C)=0.34/34=0.01 (Harley 1995 data)
    ! collatz 1991 gives rd=0.015*vcmax
    ! Farqhuar 1980 gives rd=0.011*vcmax
    ! reduce respiration by 50% in light according to Amthor (might be less, see Pinelli and Loreto, 2003)
    if (Iphoton > 1.0_wp) rdz = rdz * half !changed to 40% reduction ! >10 or >1?
    ! apply temperature correction for rd at 25 deg C to leaf level temperature
    rd          = temp_func(rdz, erd, tprime25, tk_25, tlk)
    ! Apply temperature correction to JMAX and vcmax
    jmax  = tboltz(jmaxz, ejm, toptjm, tlk, hkin)
    vcmax = tboltz(vcmaxz, evc, toptvc, tlk, hkin)
    prof%jmax(JJ)  = jmax !store jmax in gobal structure
    prof%vcmax(JJ) = vcmax !store vcmax in gobal structure
    ci_guess    = cca * 0.7_wp ! initial guess of internal CO2 to estimate Wc and Wj
    if (iswitch%tpu == 1 ) alpha_g_guess  = alpha_g_max*0.3_wp ! initial guess of glycine reduction
 !   if (JJ == 34)then
 !       print *, wj
 !   end if
    !tp~~beta_tpu are related to tpu limits Yuan 2019.12.20
    if (iswitch%tpu == 1) then
       !tp           = 6_wp
        tp = vcmax * tp_vc
        prof%tp(JJ)    = tp ! TP is 1/12 of vcmax
   !     alpha_g_max  = 0.09
   !     alpha_s_max  = 0.38
   !     n_max        = 1.21
        beta_tpu     = 3*alpha_g_max/(3*alpha_g_max+2*alpha_s_max)
        gammac       = 500.0_wp * input%o2air * (1-alpha_g_max) / tau
        if (alpha_g_max==0 .and. alpha_s_max==0) then
          beta_tpu = 0
        end if
   !     print*, beta_tpu
    end if
    gammac = gammac/1000_wp ! Because input o2 is around 210000 in ppm instead of 210. Yuan 2018.01.31
    if (iswitch%tpu == 1) phi = 2*gammac/ci_guess! ratio of oxygenation to carboxylation
    ! Compute the leaf boundary layer resistance
    ! gb_mole leaf boundary layer conductance for CO2 exchange,
    ! mol m-2 s-1
    ! RB has units of s/m, convert to mol-1 m2 s1 to be
    ! consistant with R.
    ! rb_mole = RBCO2 * .0224 * 1.01 * tlk / (met%pstat * TN0)
    rb_mole = bound_lay_res%co2 * tlk * (met%pstat273) ! met%pstat273 = rugc / (100000 * met%press_bars)
    gb_mole = one / rb_mole
    dd = gammac
    b8_dd = 8._wp * dd
    ! **************************************
    ! aphoto = PG - rd, net photosynthesis is the difference
    ! between gross photosynthesis and dark respiration. Note
    ! photorespiration is already factored into PG.
    ! ****************************************
    ! coefficients for Ball-Berry stomatal conductance model
    ! Gs = k A rh/cs + b''
    ! rh is relative humidity, which comes from a coupled
    ! leaf energy balance model
    rh_leaf = srf_vpd(tlk, zzz, leleaf) ! calculate relative humidity at leaf surface
    es_leaf = es(tlk) ! saturation vapor pressure at leaf temperature
    vpd_leaf = es_leaf - rh_leaf*es_leaf ! calculate vapor pressure deficit at leaf surface
    if (soil%camillo == 0) then
       ! sets response of stomata conductance to soil water stress
       ! calculate relative plant available water
       ! normalise total soil water by total soil water at field capacity (-33 kPa)
       theta_soil = min(max((soil%soil_mm_root-soil%soil_mm_1500_root) / &
            (soil%soil_mm_33_root-soil%soil_mm_1500_root), zero), one)
       output%c9 = theta_soil
       if (theta_soil > srf_res%fthreshold) then
          srf_res%fdrought = one
       else ! prevent instability by setting minimum to 0
          srf_res%fdrought = min(max(one-srf_res%fslope*(srf_res%fthreshold-theta_soil), zero), one)
       end if
    else
       srf_res%fdrought = one
    end if
    ! set Ball-Berry stomatal factor
    ! combine product of rh and K ball-berry times drought response factor
    k_rh = srf_res%fdrought * rh_leaf * kball(JJ)
!    print *, kball(JJ)
    ! Gs from Ball-Berry is for water vapor. It must be divided
    ! by the ratio of the molecular diffusivities to be valid
    ! for A
    k_rh = k_rh / 1.577_wp ! adjust the coefficient for the diffusion of CO2 rather than H2O
    ! parameter for Leuning mesophyll model
    if (extra_nate == 1) then
       ! !!! ASK Alex about 1.577 etc. !!!
       g0_local = g0(JJ) !empirical coefficient, intercept, converted from water to CO2
       ! slope of stomata function (-) converted from water to CO2 times drought response factor
       ! a1_local = srf_res%fdrought*a1(JJ)
       a1_local = a1(JJ)
       D0_local = D0(JJ) !empirical coefficient, intercept, converted from water to CO2
    else
       g0_local = g0(JJ)/1.577_wp !empirical coefficient, intercept, converted from water to CO2
       ! slope of stomata function (-) converted from water to CO2 times drought response factor
       a1_local = srf_res%fdrought*a1(JJ)/1.577_wp
       D0_local = D0(JJ) !empirical coefficient, intercept, converted from water to CO2
    end if
    gm          = gm_vc*vcmax !mesophyll conductance (mol m-2 s-1)
    prof%gm(JJ) = gm
    if (Iphoton < one) quad = 1 ! shortcut to dark
    ! Test for the minimum of Wc and Wj. Both have the form:
    ! W = (a ci - ad)/(e ci + b)
    ! after the minimum is chosen set a, b, e and d for the cubic solution.
    ! estimate of J according to Farquhar and von Cammerer (1981)
    ! J photon from Harley
    if (jmax > 0) then
       ! Michaelis Menten type (originally in CANOAK)
       ! j_photon = qalpha * Iphoton / sqrt(1. +(qalpha2 * Iphoton * Iphoton / (jmax * jmax)))
       ! non rectangular hyperbola with curvature parameter (see Von Cammerer 2000)
       ! j_photon = (qalpha * Iphoton + jmax - sqrt(pow(qalpha*Iphoton+jmax,2)-4.*curvature*qalpha*jmax*Iphoton)) / (2.*curvature)
       j_photon = (qalpha * Iphoton + jmax - &
            sqrt((qalpha*Iphoton+jmax)*(qalpha*Iphoton+jmax)-4._wp*curvature*qalpha*jmax*Iphoton)) / &
            (two*curvature)
    else
       j_photon = zero
    end if
    wj = j_photon * (ci_guess - dd) / (4._wp * ci_guess + b8_dd)
    wc = vcmax * (ci_guess - dd) / (ci_guess + bc)
    ! Vc under TPU limits:
    if (iswitch%tpu == 1) then
      wj2 = j_photon *(1-alpha_g_max)* (ci_guess - dd) / &
      (4._wp * ci_guess*(1-alpha_g_max) + (8_wp+16_wp*alpha_g_max+8_wp*alpha_s_max)*dd)
      wp_tpu = 3*tp*(1-alpha_g_max) * (ci_guess - dd) / &
      (ci_guess*(1-alpha_g_max)-(1_wp+3_wp*alpha_g_max+4_wp*alpha_s_max)*dd)
    end if
    ! frost and end of leaf photosynthesis and respiration
    if (time%days > time%leaffallcomplete) then ! old: 300
       wj = zero
       wj2 = zero
       j_photon = zero
       wc = zero
       wp_tpu = zero
       rd = zero
    end if
    if (iswitch%tpu == 0) then
        if (wj < wc) then
         ! for Harley and Farquhar type model for Wj
           tpu_coeff = 2
           psguess = wj
           b_ps = b8_dd
           a_ps = j_photon
           e_ps = 4._wp
        else
           tpu_coeff = 1
           psguess = wc
           b_ps = bc
           a_ps = vcmax
           e_ps = one
        end if
    else if (iswitch%tpu == 1) then! include TPU
        if (wj2 < min(wc,wp_tpu)) then
         ! for Harley and Farquhar type model for Wj
           tpu_coeff = 2
           psguess = wj2
           if (j_photon > n_max*(2*beta_tpu+6)) then
            g_tmp = (4*n_max*beta_tpu*(1+1/phi))/(j_photon - n_max*(2*beta_tpu+6))
            s_tmp = (6*n_max*(1-beta_tpu)*(1+2/phi))/(j_photon - n_max*(2*beta_tpu+6))
            alpha_g = min(g_tmp, alpha_g_max)
            alpha_s = min(s_tmp, alpha_s_max)
           else
            alpha_g = alpha_g_max
            alpha_s = alpha_s_max
           end if
           b_ps = (8_wp+16_wp*alpha_g_max+8_wp*alpha_s_max)*dd
           a_ps = j_photon*(1-alpha_g_max)
           e_ps = 4._wp*(1-alpha_g_max)

        else if (wc < min(wj2,wp_tpu)) then
           tpu_coeff = 1
           psguess = wc
           vo = 2*gammac*psguess/ci_guess
           alpha_g = min(alpha_g_max,n_max*beta_tpu/vo)
           alpha_s = min(alpha_s_max,3*n_max*(1-beta_tpu)/2*vo)
           b_ps = bc
           a_ps = vcmax
           e_ps = one
        else
           tpu_coeff = 3
           psguess=wp_tpu
           g_tmp = n_max*beta_tpu*(2/phi-1)/(6*tp+3*n_max*(2-beta_tpu))
           alpha_g = min(alpha_g_max,g_tmp)
           s_tmp = 1.5*n_max*(1-beta_tpu)*(2/phi-1)/(6*tp+3*n_max*(2-beta_tpu))
           alpha_s = min(alpha_s_max, s_tmp)
           b_ps = -(1_wp+3_wp*alpha_g_max+4_wp*alpha_s_max)*dd
           a_ps = 3*tp*(1-alpha_g_max)
           e_ps = (1-alpha_g_max)
        end if

    end if
    ! cubic coefficients that are only dependent on CO2 levels
    ! for the Ball Berry Farquhar model based on Baldocchis analytical solution
    if (iswitch%ball == 0) then
       bprime_local   = bprime(JJ)
       bprime16_local = bprime(JJ)/1.577_wp
       alpha_ps       = one + (bprime16_local / gb_mole) - k_rh
       bbeta          = cca * (gb_mole * k_rh - two * bprime16_local - gb_mole)
       ggamma         = cca * cca * gb_mole * bprime16_local
       theta_ps       = gb_mole * k_rh - bprime16_local
    else if (iswitch%ball == 1) then ! for the Leuning Farquhar model based on Knohls analytical solution
       alpha_ps   = 1/gm + 1/gb_mole
       beta_ps    = (D0_local+vpd_leaf)*gb_mole*(cca-gammac)
       gamma_ps   = a1_local*gb_mole*vpd_leaf
       delta_ps   = b_ps+e_ps*cca
       epsilon_ps = a_ps-e_ps*rd
       zeta_ps    = D0_local+vpd_leaf
       eta_ps     = a_ps*(gammac-cca)
       ! alpha_ps = (gb_mole*gm+g0*(gm+gb_mole))*(D0+vpd_leaf)
       ! beta_ps = gb_mole*gm*g0*(D0+vpd_leaf)
       ! gamma_ps = gb_mole*D0*a1
       ! delta_ps = a_ps*dd-a_ps*cca+rd*e_ps*cca+rd*e_ps
       ! epsilon_ps = gb_mole*gm*gamma_ps-beta_ps
       ! zeta_ps = gm*gamma_ps+gb_mole*gamma_ps-alpha_ps
       bprime_local   = bprime(JJ)
       bprime16_local = bprime(JJ)/1.577_wp
    else if (iswitch%ball == 2) then ! for Medlyn Farquar model without mesophyll conductance
        ! gs = g0+1.6*(1+g1/sqrt(D))*A/Cs
        g1_local = (1+g1_mly_in/sqrt(vpd_leaf))*1.6_wp
        g1_local = g1_local/1.6_wp
        g0_local = g0_mly_in*1.6_wp
!        print *, g0_mly_in
   !     alpha_ps1      = g0_local + gb_mole - g1_local * gb_mole
   !     alpha_ps2      = g0_local + gb_mole - g1_local * g0_local
        alpha_ps       = g0_local + gb_mole - g1_local * gb_mole
        bbeta          = cca * gb_mole * ( gb_mole * g1_local - two * g0_local - gb_mole)
        ggamma         = cca * cca * gb_mole * gb_mole* g0_local
        theta_ps       = gb_mole * gb_mole * g1_local - g0_local * gb_mole
        bprime_local   = g0_local
        bprime16_local = g0_local/1.6_wp
    else if (iswitch%ball == 3) then !Medlyn's model with mesophyll conductance
    ! gs = g0+1.6*(1+g1/sqrt(D))*A/Cs
        g0_local = g0_mly_in*1.6!/1.6_wp
!        print *, g0_mly_in
        g1_local = (1+g1_mly_in/sqrt(vpd_leaf))*1.6_wp
        g1_local = g1_local/1.6_wp
        alpha_ps       = g0_local + gm - g1_local*gb_mole - g1_local*gm + g0_local * gm/gb_mole
        bbeta          = cca * (gm*gb_mole*g1_local-2*g0_local*gm-g0_local*gb_mole-gb_mole*gm)
        ggamma         = cca * cca * gb_mole * gm* g0_local
        theta_ps       = gb_mole * gm * g1_local - g0_local * gm
        bprime_local   = g0_local
        bprime16_local = g0_local/1.6_wp

    end if
    rd_O2       = rd/RQ_rd
    prof%rd(JJ) = rd !store rd in gobal structure
!    prof%rd_O2(JJ) =  rd_O2!store rd_O2 in gobal structure
    ! if wj or wc are less than rd then A would probably be less than zero. This would yield a
    ! negative stomatal conductance. In this case, assume gs equals the cuticular value. This
    ! assumptions yields a quadratic rather than cubic solution for A
    if (wj <= rd) quad = 1
    if (iswitch%tpu == 1 .and. wj2 <= rd) quad = 1
    if (wc <= rd) quad = 1
    if (iswitch%tpu == 1 .and. wp_tpu <= rd) quad = 1
    if (quad == 0) then
       ! cubic solution: A^3 + p A^2 + q A + r = 0
       ! for the Ball Berry Farquhar model based on Baldocchis analytical solution
       if (iswitch%ball == 0) then
          denom = e_ps * alpha_ps
          Pcube = (e_ps * bbeta + b_ps * theta_ps - a_ps * alpha_ps + e_ps * rd * alpha_ps)
          Pcube = Pcube/denom
          Qcube = e_ps * ggamma + (b_ps * ggamma / cca) - a_ps * bbeta + a_ps * dd * theta_ps + &
               e_ps * rd * bbeta + rd * b_ps * theta_ps
          Qcube = Qcube/denom
          Rcube = -a_ps * ggamma + a_ps * dd * (ggamma / cca) + e_ps * rd * ggamma + &
               rd * b_ps * ggamma / cca
          Rcube = Rcube/denom
       else if (iswitch%ball == 1) then! for the Leuning Farquhar model based on Knohls analytical solution
          denom = -alpha_ps*e_ps*gamma_ps + alpha_ps*e_ps*g0_local*zeta_ps + e_ps*zeta_ps
          Pcube = alpha_ps*epsilon_ps*gamma_ps - alpha_ps*epsilon_ps*g0_local*zeta_ps - &
               alpha_ps*e_ps*g0_local*beta_ps - zeta_ps*epsilon_ps - &
               beta_ps*e_ps+delta_ps*gamma_ps - delta_ps*zeta_ps*g0_local
          Pcube = Pcube/denom
          Qcube = alpha_ps*epsilon_ps*g0_local*beta_ps + beta_ps*epsilon_ps + &
               delta_ps*beta_ps*g0_local + eta_ps*gamma_ps-eta_ps*zeta_ps*g0_local + &
               delta_ps*gamma_ps*rd - delta_ps*zeta_ps*g0_local*rd
          Qcube = Qcube/denom
          Rcube = eta_ps*beta_ps*g0_local + delta_ps*beta_ps*g0_local*rd
          Rcube = Rcube/denom
          ! denom = e_ps * (-zeta_ps)
          ! Pcube = (a_ps-e_ps*rd)*zeta_ps-e_ps*alpha_ps*cca*gb_mole+(b_ps+e_ps*cca)*epsilon_ps
          ! Pcube = Pcube/denom
          ! Qcube = (b_ps+e_ps*cca)*beta_ps*cca*gb_mole+delta_ps*epsilon_ps
          ! Qcube = Qcube/denom
          ! Rcube = delta_ps*beta_ps*cca*gb_mole
          ! Rcube = Rcube/denom
       else if (iswitch%ball == 2) then ! Medlyn's stomatal model without Mesophyll conductance
          denom = e_ps * alpha_ps
          Pcube = (e_ps * bbeta + b_ps * theta_ps - a_ps * alpha_ps + e_ps * rd * alpha_ps)
          Pcube = Pcube/denom
          Qcube = e_ps * ggamma + (b_ps * ggamma / cca) - a_ps * bbeta + a_ps * dd * theta_ps + &
                  e_ps * rd * bbeta + rd * b_ps * theta_ps
          Qcube = Qcube/denom
          Rcube = e_ps * rd * ggamma + a_ps * dd * (ggamma / cca) + b_ps * rd * (ggamma / cca) - &
                  a_ps * ggamma
          Rcube = Rcube/denom
       else if (iswitch%ball == 3) then ! Medlyn's stomatal model with Mesophyll conductance
          denom = e_ps * alpha_ps
          Pcube = (e_ps * bbeta + b_ps * theta_ps - a_ps * alpha_ps + e_ps * rd * alpha_ps)
          Pcube = Pcube/denom
          Qcube = e_ps * ggamma + (b_ps * ggamma / cca) - a_ps * bbeta + a_ps * dd * theta_ps + &
                  e_ps * rd * bbeta + rd * b_ps * theta_ps
          Qcube = Qcube/denom
          Rcube = e_ps * rd * ggamma + a_ps * dd * (ggamma / cca) + b_ps * rd * (ggamma / cca) - &
                  a_ps * ggamma
          Rcube = Rcube/denom
       end if
       ! Use solution from Numerical Recipes from Press
       P2 = Pcube * Pcube
       P3 = P2 * Pcube
       Q = (P2 - 3.0_wp * Qcube) / 9.0_wp
       R = (2.0_wp * P3 - 9.0_wp * Pcube * Qcube + 27.0_wp * Rcube) / 54.0_wp
       ! Test = Q ^ 3 - R ^ 2
       ! if test >= O then all roots are real
       rr = R*R
       qqq = Q*Q*Q
       ! real roots
       arg_U = R / sqrt(qqq)
       ang_L = acos(arg_U)
       root1 = -two * sqrt(Q) * cos(ang_L / 3.0_wp) - Pcube / 3.0_wp
       root2 = -two * sqrt(Q) * cos((ang_L + pi2) / 3.0_wp) - Pcube / 3.0_wp
       root3 = -two * sqrt(Q) * cos((ang_L - pi2) / 3.0_wp) - Pcube / 3.0_wp
       ! rank roots #1,#2 and #3 according to the minimum, intermediate and maximum
       ! value
       if (root1 <= root2 .and. root1 <= root3) then
          minroot=root1
          if (root2 <= root3) then
             midroot=root2
             maxroot=root3
          else
             midroot=root3
             maxroot=root2
          end if
       end if
       if (root2 <= root1 .and. root2 <= root3) then
          minroot=root2
          if (root1 <= root3) then
             midroot=root1
             maxroot=root3
          else
             midroot=root3
             maxroot=root1
          end if
       end if
       if (root3 <= root1 .and. root3 <= root2) then
          minroot=root3
          if (root1 < root2) then
             midroot=root1
             maxroot=root2
          else
             midroot=root2
             maxroot=root1
          end if
       end if ! end of the loop for real roots
       ! find out where roots plop down relative to the x-y axis
       if (minroot > zero .and. midroot > zero .and. maxroot > zero) aphoto = minroot
       if (minroot < zero .and. midroot < zero .and. maxroot > zero) aphoto = maxroot
       if (minroot < zero .and. midroot > zero .and. maxroot > zero) aphoto = midroot
       ! Ophoto = aphoto * prof%ROC_leaf_air(JJ)
       ! Here A = x - p / 3, allowing the cubic expression to be expressed
       ! as: x^3 + ax + b = 0
       ! aphoto=root3 ! back to original assumption
       ! also test for sucrose limitation of photosynthesis, as suggested by
       ! Collatz. Js=Vmax/2
       j_sucrose = vcmax * half - rd
       if (j_sucrose < aphoto) aphoto = j_sucrose
       cs = cca - aphoto / gb_mole
       if (cs > 3._wp*cca) cs = input%co2air
       gpp = aphoto + rd
       gpp_o2 = gpp * prof%ROC_leaf_air(JJ)
       if (iswitch%tpu == 1) then
        gpp_o2 = gpp_o2 + 2*alpha_g + alpha_s
       end if
       Ophoto = gpp_o2 - rd_O2
       ! Stomatal conductance for water vapor
       ! forest are hypostomatous.
       ! Hence we don''t divide the total resistance
       ! by 2 since transfer is going on only one side of a leaf.
       if (iswitch%ball == 0) then
          gs_leaf_mole = (k_rh*1.577_wp * aphoto / cs) + bprime_local
       else if(iswitch%ball == 1) then
          gs_leaf_mole = a1_local*1.577_wp*aphoto/((cs-dd)*(1+vpd_leaf/D0_local)) + g0_local*1.577_wp
       else if (iswitch%ball == 2) then
          gs_leaf_mole = g1_local*1.6_wp*aphoto/cs + g0_local
       else if (iswitch%ball == 3) then
          gs_leaf_mole = g1_local*1.6_wp*aphoto/cs + g0_local
       end if
       ! convert Gs from vapor to CO2 diffusion coefficient based on diffusivities in Massman (1998)
       gs_co2 = gs_leaf_mole /1.577_wp
       if (aphoto < zero) then
          call message('PHOTOSYNTHESIS: ','aphoto<0 should not be here: ', num2str(JJ), num2str(time%daytime))
       end if
       ci = cs - aphoto / gs_co2
       cc = ci - aphoto / gm
       if (iswitch%ball == 0) then
          wj = j_photon * (ci - dd) / (4._wp * ci + b8_dd)
          wc = vcmax * (ci - dd) / (ci + bc)
       else if (iswitch%ball == 1) then
          wj = j_photon * (cc - dd) / (4._wp * cc + b8_dd)
          wc = vcmax * (cc - dd) / (cc + bc)
       else if (iswitch%ball == 2) then
          wj = j_photon * (ci - dd) / (4._wp * ci + b8_dd)
          wc = vcmax * (ci - dd) / (ci + bc)
       else if (iswitch%ball == 3) then
          wj = j_photon * (cc - dd) / (4._wp * cc + b8_dd)
          wc = vcmax * (cc - dd) / (cc + bc)
       end if

       if (iswitch%tpu == 1) then
        wj2 = j_photon * (ci - dd) / (4._wp * ci + (8_wp+16_wp*alpha_g+8_wp*alpha_s)*dd)
        wp_tpu = 3*tp * (ci - dd) / (ci-(1_wp+3_wp*alpha_g+4_wp*alpha_s)*dd)
        if (iswitch%ball == 1 .OR. iswitch%ball == 3) then
            wc = vcmax * (cc - dd) / (cc + bc)
            wj2 = j_photon * (cc - dd) / (4._wp * cc + (8_wp+16_wp*alpha_g+8_wp*alpha_s)*dd)
            wp_tpu = 3*tp * (cc - dd) / (cc-(1_wp+3_wp*alpha_g+4_wp*alpha_s)*dd)
        end if
       end if
       ! stomatal conductance is mol m-2 s-1
       ! convert back to resistance (s/m) for energy balance routine
       gs_m_s = gs_leaf_mole * tlk * met%pstat273
       ! need point to pass rstom out of subroutine
       rstompt = one / gs_m_s
       !   ! to compute ci, Gs must be in terms for CO2 transfer
       !   ci = cs - aphoto / gs_co2
       !   ! compute cc with mesophyll conductance
       !   cc = ci - aphoto / gm
       ! if A < 0 then gs should go to cuticular value and recalculate A
       ! using quadratic solution
       if (aphoto <= zero) quad = 1
    end if ! quad == 0
    if (quad == 1) then
       ! if aphoto < 0 set stomatal conductance to cuticle value
       gs_leaf_mole = bprime_local
       gs_co2       = gs_leaf_mole /1.577_wp
       ! stomatal conductance is mol m-2 s-1
       ! convert back to resistance (s/m) for energy balance routine
       gs_m_s = gs_leaf_mole * tlk * (met%pstat273)
       ! need point to pass rstom out of subroutine as a pointer
       rstompt = one / gs_m_s
       ! a quadratic solution of A is derived if gs=ax, but a cubic form occurs
       ! if gs =ax + b. Use quadratic case when A is less than zero because gs will be
       ! negative, which is nonsense
       if (Iphoton < one) then ! shortcut dark
          gpp    = zero
          aphoto = -rd

       else
          ps_1 = cca * gb_mole * gs_co2
          delta_1 = gs_co2 + gb_mole
          denom = gb_mole * gs_co2
          Aquad1 = delta_1 * e_ps
          Bquad1 = -ps_1 * e_ps - a_ps * delta_1 + e_ps * rd * delta_1 - b_ps * denom
          Cquad1 = a_ps * ps_1 - a_ps * dd * denom - e_ps * rd * ps_1 - rd * b_ps * denom
          product = Bquad1 * Bquad1 - 4._wp * Aquad1 * Cquad1
          if (product >= zero) then
             sqrprod= sqrt(product)
             aphoto = (-Bquad1 - sqrprod) / (two * Aquad1)
             ! Tests suggest that APHOTO2 is the correct photosynthetic root when
             ! light is zero because root 2, not root 1 yields the dark respiration
             ! value rd.)
          else
             aphoto = zero
          end if
          ! correct for gpp>0 - should only be numerical adjustments
          gpp = max(aphoto+rd, zero)
          aphoto = gpp-rd
       end if

       gpp_o2 = gpp * prof%ROC_leaf_air(JJ)
       if (iswitch%tpu == 1) then
           gpp_o2 = gpp_o2 + 2*alpha_g + alpha_s
       end if
       Ophoto = gpp_o2 - rd_O2
       cs = cca - aphoto / gb_mole
       ci = cs - aphoto / gs_co2
       cc = ci - aphoto / gm
    end if
    if (ci < 0) then
        print *, cs
        print *, aphoto
        print *, gs_co2
    end if
    if (cc < 0) then
        print *, gm
    end if
    if (gpp < zero) then
       call message('PHOTOSYNTHESIS: ','gpp<0 should not be here: ', num2str(JJ), num2str(time%daytime))
    end if
    ! compute photosynthesis with units of mg m-2 s-1 and pass out as pointers
    ! A_mg = APHOTO * 44 / 1000
    A_mgpt      = aphoto * mass_CO2 * e3
    O_pt        = Ophoto ! units of umol m-2 s-1
    GPPpt       = gpp
    GOPpt       = gpp_o2
    resppt      = rd ! leaf dark respiration!
    resOppt     = rd_O2
    ROC_rd      = 1/RQ_rd
    cipnt       = ci
    cspnt       = cs
    ccpnt       = cc
    cicapnt     = ci/cca
    cccapnt     = cc/cca
    wcpnt       = wc
    if (iswitch%tpu == 0) then
        wjpnt   = wj
    else
        wjpnt   = wj2
    end if
    wppnt       = wp_tpu
    tpupnt      = tpu_coeff
    rh_leafpnt  = rh_leaf
    vpd_leafpnt = vpd_leaf
 !   if (JJ == 34)then
 !       print *, wj, wj2
 !   end if

  END SUBROUTINE photosynthesis


  ! ------------------------------------------------------------------
  SUBROUTINE soil_respiration()
    ! Computes soil respiration
    USE constants,  ONLY: zero, half, one, e3, mass_co2
    USE types,      ONLY: soil, iswitch, input, bole, prof
    USE parameters, ONLY: extra_nate, rsoil1, rsoil2
    USE setup,      ONLY: ncl

    IMPLICIT NONE

    REAL(wp), PARAMETER :: fautleaf = 0.40_wp
    REAL(wp), PARAMETER :: ccost = 0.1_wp
    REAL(wp) :: zass   ! GPP of entire canopy layer
    REAL(wp) :: zrd    ! dark respiration of leaves
    REAL(wp) :: zauto  ! autotrophic respiration of entire canopy
    REAL(wp) :: ztmp2d ! maintenace respiration

    zass   = zero   ! GPP of entire canopy layer
    zrd    = zero   ! dark respiration of leaves
    zauto  = zero  ! autotrophic respiration of entire canopy
    ztmp2d = zero ! maintenace respiration
    ! OLD
    ! soil.base_respiration=0.71 ! new value Recosystem=1.34
    ! soil.respiration_mole = soil.base_respiration * exp((51609. / 8.314) * ((1. / 273.) - 1. / (soil.T_15cm + TN0)))
    ! ! soil wetness factor from the Hanson model, assuming constant and wet soils
    ! soil.respiration_mole *= 0.86
    if (extra_nate == 1) then
       ! Equation of Nate McDowell for Juniper site
       ! Take 5cm instead of Nate''s 2cm because Canoak''s 2cm is too variable
       soil%respiration_mole = 0.096_wp * soil%theta(4,1) * 100._wp + 0.5089_wp
       ! print*, "R: ", soil%respiration_mole, soil%theta(4,1)
    else
       ! calculate soil respiration based on chamber measurements by Astrd Soe for Hainich site (Soe 2003) Astrid R. B. S�e  Nina Buchmann!! Yuan 2018.07.30
       ! soil respiraton = f(Tsoil in 5 cm)
       if (iswitch%soil_resp_temp == 1) then
          ! set soil reference temperature to input soil temperature
          soil%SR_ref_temp = input%tsoil
       else
          ! set soil reference temperature to modeled soil temperature in 5 cm, depending on switch set in parameter file
          soil%SR_ref_temp = soil%T_soil(4)
       end if
       ! canisotope v3.1
       ! soil.respiration_mole = 0.71*exp(0.09*soil.SR_ref_temp)    !factor changed, orignal 0.11
       ! from Hainich soil respiration measurements
       ! print *, rsoil1,rsoil2
       soil%respiration_mole = rsoil1 * exp(rsoil2*soil%SR_ref_temp)!0.11,0.09,0.06
       if (iswitch%oxygen == 1 .and. soil%SR_ref_temp>=0) then
        ! add temperature function of soil ROC of Quercus calliprinos
        ! by Angert et al JGR Biogeosciences 2020
        soil%ROC_soil_air = 0.15_wp*log(soil%SR_ref_temp)+0.34_wp
        soil%ROC_soil_air = 1/soil%ROC_soil_air
       end if
    end if

    if (iswitch%bethy_resp == 1) then
       ! Compute soil respiration as autotrophic and heterotrophic respiration based on BETHY (Knorr PhD thesis 1997)
       ! Autotrophic is maintenance and growth respiration
       ! Rleaf = 0.4 Rmaintenance => Rmaintenance = Rleaf/0.4
       ! Rgrowth = 0.25 * NPP = 0.25*(GPP - Rmaintenence - Rgrowth) => Rgrowth = 0.25/(1+0.25) * (GPP - Rmaintenance)
       zass   = sum(prof%dPsdz(1:ncl) + prof%dRESPdz(1:ncl))
       zrd    = sum(prof%dRESPdz(1:ncl))
       ztmp2d = zrd/fautleaf  ! maintenance respiration of entire plants
       zauto  = ztmp2d + max((zass-ztmp2d)*(ccost/(one+ccost)), zero) !maintenance + growth respiration
       ! soil autotroph respiration = total autotroph - leaf - bole
       soil%respiration_auto   = zauto - zrd - bole%respiration_mole
       soil%respiration_hetero = soil%respiration_mole - soil%respiration_auto
    else
       ! assume autotroph and heterotroph are both 50% of soil respiration
       soil%respiration_auto   = soil%respiration_mole * half
       soil%respiration_hetero = soil%respiration_mole * half
    end if
    ! convert soilresp to mg m-2 s-1 from umol m-2 s-1
    soil%respiration_mg = soil%respiration_mole * mass_CO2 * e3

  END SUBROUTINE soil_respiration


  ! ------------------------------------------------------------------
  FUNCTION srf_vpd(tlk, Z, leleafpt)
    ! this function computes the relative humidity at the leaf surface for
    ! application in the Ball Berry Equation
    ! latent heat flux, LE, is passed through the function, mol m-2 s-1
    ! and it solves for the humidity at leaf surface
    USE constants,  ONLY: one
    USE types,      ONLY: fact, bound_lay_res, prof
    USE parameters, ONLY: delz
    USE utils,      ONLY: es

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: tlk
    REAL(wp), INTENT(IN) :: Z
    REAL(wp), INTENT(IN) :: leleafpt
    REAL(wp)             :: srf_vpd

    INTEGER(i4) :: J
    REAL(wp) :: rhov_srf, e_srf, vpd_srf
    REAL(wp) :: es_leaf

    es_leaf   = es(tlk)       ! saturation vapor pressure at leaf temperature
    J         = nint(Z/delz,kind=i4) ! layer number
    rhov_srf  = (leleafpt / (fact%latent)) * bound_lay_res%vapor + prof%rhov_air(J,1) ! kg m-3
    e_srf     = rhov_srf * tlk / 0.2165_wp ! mb
    vpd_srf   = es_leaf - e_srf ! mb
    srf_vpd   = one - vpd_srf / es_leaf ! 0 to 1

  END FUNCTION srf_vpd


END MODULE energy_carbon
