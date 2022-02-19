MODULE oxygen

  ! This module contains the oxygen release and uptake of Canveg

  ! Written June 2021, Yuan Yan
    USE kinds, ONLY: wp, i4
    USE constants,     ONLY: zero, one, undef, mass_air, mass_O2
    USE setup,         ONLY: ncl, ntl
    USE types,         ONLY: bole, prof, soil, met, input, solar, fact
!    USE isotope_utils, ONLY: invdelta1000
    USE transport,     ONLY: conc, conc_seperate
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: gross_emission
  PUBLIC  :: uptake
  PUBLIC  :: OXYFLUX

  CONTAINS

!  FUNCTION gross_emission (J_CO2,J_extra,I_photon,ETRmax)
!
!    ! This subroutine computes primary O2 release via water split reaction,
!    ! where 4 mole electrons are provided while 1 mole O2 released.
!
!  IMPLICIT NONE
!        REAL(wp), INTENT(IN) :: J_CO2, J_extra, I_photon, ETRmax
!        REAL(wp) :: gross_o2, o2_test
!        REAL(wp) :: ETR1, ETR2 ! electron transport rate
!
!        ETR1 = I_photon*0.35*0.95
!        ETR2 = ETRmax*tanh(0.35*I_photon/ETRmax)
!        o2_test = ETR2/4
!        gross_o2 = (J_CO2+J_extra)/4
!
!  END FUNCTION gross_emission

  SUBROUTINE gross_emission (carboxlation,vo_vc,o2_emit,Ja,Jn)


    ! This subroutine computes primary O2 release via water split reaction,
    ! where 4 mole electrons are provided while 1 mole O2 released.
    ! The total O2 evolution (in support of PCR and PCO cycle activity) at PS II derived from the
    ! splitting of H2O is therefore:
    ! gross emission = JA/4 = Vc + Vo = (1+phi)*vc
    ! where, vc (carboxlation)= min{wc,wj,wp}

    USE types, ONLY: nitrogen
    IMPLICIT NONE
        REAL(wp), INTENT(IN) :: carboxlation,vo_vc
        REAL(wp), INTENT(OUT) :: o2_emit, Ja, Jn

        Jn = nitrogen%J_extra
        Ja = 4*(1+vo_vc)*carboxlation
        o2_emit = (Ja+Jn)/4
!        print *, Jn/Ja

  END SUBROUTINE gross_emission

  FUNCTION uptake (oxygenation,dark_resp)

    ! Rubisco oxygenase activity is a major component of the leaf�s oxygen uptake processes. One
    !mol of O2 is consumed per mol of RuBP oxygenated and a further half mol of O2 is consumed
    !in the PCO cycle by glycolate oxidase in the peroxisomes (Badger 1985). Thus the
    !rubisco-linked O2 uptake processes are given by 1.5Vo = 1.5*phi*vc, where vc = min{wc,wj,wp}

    ! Uo = 1.5Vo + Rd + MAP

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: oxygenation
    REAL(wp), INTENT(IN) :: dark_resp
    REAL(wp)             :: MAP, Ro
    REAL(wp)             :: uptake

    Ro = dark_resp
    MAP = 0
    uptake = 1.5*oxygenation+Ro+MAP

  END FUNCTION uptake

  SUBROUTINE OXYFLUX
    ! This subroutine, LEAF_O2, computes O2 flux contributed by the leaf.

    ! This subroutine computes O2: CO2 ratio (ROC) as deltaO2/deltaCO2,
    ! BUT firstly we use fixed ROC from the publications:
    ! ROC leaf =1 by assuming the photosynthesis products are all sugar

    ! Later there will be O2 flux from bole and soil respirations by 2 other subroutines.
    ! Ci/Ca is computed using the photosynthesis, respiration and leaf energy balance algorithms
    ! of CANOAK. The Farquhar photosynthesis model is coupled to the Ball-Berry stomatal conductance
    ! model.

    ! A Lagrangian random walk model is used to compute the turbulence Dispersion matrix
    ! and compute concentration profiles in and above the canopy. This scheme accounts
    ! for counter gradient transfer.

    ! Summary of variables computed:

    !   photosynthetic flux of O2_leaf = A * ROC_leaf
    !     This relation is in terms of net Ps, A, not gross Ps
    USE parameters,       only: o2
  IMPLICIT NONE
  ! respiratory quotient of leaf respiration as a function of leaf temperature. Yuan 2019.03.11
        REAL(wp) :: source_sun, source_shade!, RQ_sun, RQ_shade
        REAL(wp) :: o2air
        INTEGER(i4) :: j
        do j=1, ncl

 !           sun_A = max(prof%sun_A(j)*solar%prob_beam(j), zero) ! photosynthetic CO2 of sunlit&shade leaves
 !           shd_A = max(prof%shd_A(j)*solar%prob_shd(j), zero)
 !           if (sun_A > zero .or. shd_A > zero) then
 !               source_sun         = sun_A*prof%ROC_leaf_air(j)
 !               source_shade       = shd_A*prof%ROC_leaf_air(j)
 !           else
 !               source_sun         = zero
 !               source_shade       = zero

 !           end if

       !scaling from leaf area to ground area for each layer
       ! source_O2 is positive when emitted, negetive when consumed
       ! prof%source_O2(j) = prof%dLAIdz(j) * (source_shade+source_sun) ! o2 source/sink from photosynthesis
     ! RQ = -0.0147*(T) +1.24
     !RQ_sun = -0.0147_wp*prof%sun_tleaf(j)+1.24_wp
     !RQ_shade = -0.0147_wp*prof%shd_tleaf(j) +1.24_wp
     !prof%RQ_sun(j) = RQ_sun
     !prof%RQ_shade(j) = RQ_shade
        end do
!     rd_o2 updated in photosynthesis()
!     prof%rd_O2(1:ncl) = prof%dRESPdz_sun(1:ncl)/prof%RQ_sun(1:ncl)+prof%dRESPdz_shd(1:ncl)/prof%RQ_shade(1:ncl)
    ! prof%gpp_O2(1:ncl) = prof%dPsdz(1:ncl)*prof%ROC_leaf_air(1:ncl) + prof%rd_O2(1:ncl)! o2 flux of gpp = psn + rd
    ! prof%gpp_O2(1:ncl) = prof%dPsdz_O2(1:ncl) + prof%dRESPdz_O2(1:ncl)! o2 flux of gpp = psn + rd
     prof%dGOPdz_sun(1:ncl) = prof%dLAIdz(1:ncl) * prof%sun_GOP(1:ncl) * solar%prob_beam(1:ncl)
     prof%dGOPdz_shd(1:ncl) = prof%dLAIdz(1:ncl) * prof%shd_GOP(1:ncl) * solar%prob_shd(1:ncl)
     prof%dGOPdz(1:ncl)     = prof%dGOPdz_sun(1:ncl) + prof%dGOPdz_shd(1:ncl)
     prof%gpp_O2(1:ncl) = prof%dGOPdz(1:ncl)
!     print *, prof%shd_GOP
!print *, sum(prof%gpp_O2(1:ncl))
!print *, sum(prof%dGPPdz(1:ncl))
!print *, sum(prof%dPsdz(1:ncl)*prof%ROC_leaf_air(j))
!     prof%gpp_O2(1:ncl) = prof%dGPPdz(1:ncl)* prof%ROC_leaf_air(1:ncl)
!    prof%source_co2(1:ncl) = -prof%dPsdz(1:ncl) + bole%layer(1:ncl) (negative + positive)
!    prof%source_O2(1:ncl) = prof%dPsdz(1:ncl)* prof%ROC_leaf_air(1:ncl)- &! add o2 source/sink from bole respiration
!                            bole%layer(1:ncl)* prof%ROC_bole_air(1:ncl)
     prof%source_O2(1:ncl) = prof%gpp_O2(1:ncl) - prof%dRESPdz_O2(1:ncl) - bole%layer(1:ncl)* prof%ROC_bole_air(1:ncl)
   !  print *, prof%source_O2(33)
   !  print *, prof%dPsdz(33)
   !  print *, bole%layer(33)
     soil%o2 = -soil%respiration_mole*soil%ROC_soil_air ! negative!!!
   !  print *, prof%dGPPdz(33), bole%layer(33), soil%respiration_mole
   !  print *, prof%gpp_O2(33), bole%layer(33)*prof%ROC_bole_air(33), soil%o2
!     o2air = o2*1000 ! convert to ppm
     o2air = o2
     fact%o2 = (mass_air/mass_O2)*met%air_density_mole ! convert from density to mole concentration
!print *,fact%o2
!print *,fact%co2
!     call conc(prof%source_O2, prof%O2_air, input%o2air, soil%o2, fact%o2) ! use input o2 con. Yuan 2018.02.14
     call conc_seperate(prof%source_O2, prof%O2_air, input%o2air, soil%o2, prof%O2_soil, prof%O2_disp, fact%o2)
!    calculate ROC of each layer
     prof%ROC_layer(1:ncl) = abs(prof%source_O2(1:ncl)/prof%source_co2(1:ncl))
!    print *, prof%ROC_layer
!    call conc(prof%source_co2, prof%co2_air, input%co2air, soil%respiration_mole, fact%co2)

!    call conc(prof%sour13co2, prof%c13cnc, co2_13, soil%resp_13, fact%co2)
    ! this is the o2/(co2) ratio of air in the canopy
    ! prof%R13_12_air(1:ntl) = prof%c13cnc(1:ntl)/prof%co2_air_filter(1:ntl)
  END SUBROUTINE OXYFLUX

END MODULE oxygen
