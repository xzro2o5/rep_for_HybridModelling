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
  PUBLIC  :: uptake, O_to_N, N_to_O
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

  SUBROUTINE gross_emission(carboxlation,vo_vc,o2_emit,Ja,Jn)


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

        Jn = nitrogen%J_glu
        Ja = 4*(1+vo_vc)*carboxlation
        o2_emit = (Ja+Jn)/4
!        print *, Jn/Ja

  END SUBROUTINE gross_emission



  SUBROUTINE O_to_N (An,carboxylation,Vo_Vc,ER_An,Rd,leaf_T,gly,serine,Etot,En, Uo, dark_resp_O, &
    Ja, J_glu, J_Busch, N_demand, N_tot, source_Busch, source_NO3,source_NO2,source_NH4)

    USE utils,        ONLY: ER_rd_func
    USE nitrogen_assimilation, ONLY: N_fraction
    USE parameters, ONLY: cn_bulk

    REAL(wp), INTENT(IN)  :: An,carboxylation,Vo_Vc,ER_An,Rd,leaf_T,gly,serine
    REAL(wp), INTENT(OUT) :: Etot, En, Uo, dark_resp_O, J_glu, Ja, J_Busch, &
    N_demand, N_tot, source_Busch, source_NO3, source_NO2, source_NH4
    REAL(wp)              :: ER_rd, MAP, source_glu
    REAL(wp)              :: Jtot, f1, f2, f3

    ! oxygen uptake by leaf dark resp:
    MAP = 0
    ER_rd = ER_rd_func(leaf_T)
    dark_resp_O = Rd * ER_rd
    ! net Ass oxygen
    En = An*ER_An

    if (En<An) MAP = An-En
    Uo = 1.5*Vo_Vc*carboxylation+dark_resp_O + MAP
    Etot = En + Uo

    ! calculate N assimilation:
    ! total electrons
    Jtot = Etot * 4
    ! electrons for CO2 assimilation
    Ja = 4*(1+Vo_Vc)*carboxylation
    ! electrons by Busch's photorespiration
    J_Busch = (8*gly+4*serine)*Vo_Vc*carboxylation
    ! e- required for B assimilation
    J_glu = Jtot-Ja-J_Busch
    J_glu = max(J_glu,zero)
    call N_fraction (f1,f2,f3)
    source_glu = J_glu/(f1*10+f2*8+f3*2)
    source_Busch = (gly+2*serine/3)*Vo_Vc*carboxylation
    source_NO3 = source_glu*f1
    source_NO2 = source_glu*f2
    source_NH4 = source_glu*f3
    N_tot = source_Busch + source_glu
    N_demand = carboxylation/cn_bulk


! save N source in structures:
!    nitrogen%Ja = Ja
!    nitrogen%J_glu = J_glu
!    nitrogen%J_Busch = J_Busch
!    nitrogen%Busch_mol   = source_Busch
!    nitrogen%nitrate_mol = source_NO3
!    nitrogen%nitrite_mol = source_NO2
!    nitrogen%ammonia_mol = source_NH4


  END SUBROUTINE O_to_N

  SUBROUTINE N_to_O (carboxylation,gross_CO2, Vo_Vc,Rd,leaf_T,gly,serine, cnleaf,nmax_extra,Etot, En, Uo, dark_resp_O, &
    Ja, J_glu, J_Busch, N_demand, N_tot, source_Busch, source_NO3, source_NO2, source_NH4)

    USE utils,        ONLY: ER_rd_func
!    USE OXYGEN,       ONLY: uptake
    USE nitrogen_assimilation, ONLY: N_fraction
    USE constants,  ONLY: one
    USE types,      ONLY: iswitch, nitrogen, time
    USE parameters, ONLY: cn_bulk,  n_mult
    USE setup,      ONLY: ncl


    REAL(wp), INTENT(IN)  :: carboxylation,gross_CO2, Vo_Vc,Rd,leaf_T,gly,serine,cnleaf,nmax_extra
    REAL(wp), INTENT(OUT) :: Etot, En, Uo, dark_resp_O, J_glu, Ja, J_Busch, &
    N_demand, N_tot, source_Busch, source_NO3, source_NO2, source_NH4
    REAL(wp)              :: ER_rd, MAP, source_glu
    REAL(wp)              :: J1,J2,J3, JCO2, Jtot, f1, f2, f3
    ! INTEGER(i4), INTENT(IN) :: layer

        ! to calculate N assimilation accompany with carbon assimilation
    ! electron requirements
    ! assuming N reducted to glutamate, CO2 reducted to carbohydrate,
    ! co2~4e, nitrate~10e, nitrite~8e and ammonia~2e
    ! J_extra/J_a =(e_source)*ncbulk/e_co2 without N limit
    ! J_extra/J_a =(e_source)*N_supply/(C_ass*e_co2) under N limit

    ! determine N assimilation amount:
    source_Busch = (gly+2*serine/3)*Vo_Vc*carboxylation
    N_demand = carboxylation/cnleaf
    !N_demand = gross_CO2/cn_bulk

    SELECT CASE (iswitch%n_limit)

    CASE (0)


        !n_ass = gpp/cn_bulk
        ! use the vertical N:C profile (Bachofen et al. 2020)

        ! spring
!        if (time%days >= time%leafout .and. time%days < time%leaffull) then
!            if (layer<=10) then
!                N_C = 0.065
!            else if (layer >10 .and. layer <= 20) then
!                N_C = 0.066
!            else if (layer >20 .and. layer <= 30) then
!                N_C = 0.071
!            else if (layer >30 .and. layer <= 40) then
!                N_C = 0.07
!            end if
!
!        end if
!
!!print *, time%days
!        ! summer
!        if (time%days >= time%leaffull .and. time%days <= time%leaffall) then
!            if (layer<=10) then
!                N_C = 0.0547
!            else if (layer >10 .and. layer <= 20) then
!                N_C = 0.0531
!            else if (layer >20 .and. layer <= 30) then
!                N_C = 0.0502
!            else if (layer >30 .and. layer <= 40) then
!                N_C = 0.0498
!            end if
!
!        end if
        source_glu = N_demand
    CASE (1)
        source_glu = nmax_extra
        print *, 'extra N supply=', nmax_extra
    CASE (2)
        source_glu = min(nmax_extra, N_demand) ! or n_supply/ncl, per layer
        !source_glu = min(min(n_supply,n_max),N_demand) ! or n_supply/ncl, per layer
    CASE (3)
      source_glu = n_mult*source_Busch

    END SELECT

! calculate N source fractions:
    call N_fraction (f1,f2,f3)
    source_NO3 = source_glu*f1
    source_NO2 = source_glu*f2
    source_NH4 = source_glu*f3

    N_tot = source_glu+source_Busch
!    J_nitrate = J_co2*(10*source_nitrate/gpp)/4
!    J_nitrite = J_co2*(8*source_nitrite/gpp)/4
!    J_ammonia = J_co2*(2*source_ammonia/gpp)/4

! from nitrogen to electron
    J1 = 10*source_NO3
    J2 = 8*source_NO2
    J3 = 2*source_NH4
    J_glu = J1 + J2 + J3

! e- by carboxylation and oxygenation:
    MAP = 0
    ER_rd = ER_rd_func(leaf_T)
    dark_resp_O = Rd * ER_rd

    ! electrons for CO2 assimilation
    Ja = 4*(1+Vo_Vc)*carboxylation
    JCO2 = 4*gross_CO2
    ! electrons by Busch's photorespiration
    J_Busch = (8*gly+4*serine)*Vo_Vc*carboxylation
    ! total electrons
    Jtot = J_glu+JCO2+J_Busch
    ! total O2 emission:
    Etot = Jtot/4
    Uo = 1.5*Vo_Vc*carboxylation+dark_resp_O + MAP
    En = Etot-Uo
!print *, "GOP and GPP in N to O:",Etot, gross_CO2
    ! save N source in structures:
!    nitrogen%Ja = Ja
!    nitrogen%J_glu = J_glu
!    nitrogen%J_Busch = J_Busch
!    nitrogen%Busch_mol   = source_Busch
!    nitrogen%nitrate_mol = source_NO3
!    nitrogen%nitrite_mol = source_NO2
!    nitrogen%ammonia_mol = source_NH4

  END SUBROUTINE N_to_O

  FUNCTION uptake (oxygenation,Ro)

    ! Rubisco oxygenase activity is a major component of the leaf’s oxygen uptake processes. One
    !mol of O2 is consumed per mol of RuBP oxygenated and a further half mol of O2 is consumed
    !in the PCO cycle by glycolate oxidase in the peroxisomes (Badger 1985). Thus the
    !rubisco-linked O2 uptake processes are given by 1.5Vo = 1.5*phi*vc, where vc = min{wc,wj,wp}

    ! Uo = 1.5Vo + Rd + MAP

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: oxygenation
    REAL(wp), INTENT(IN) :: Ro
    REAL(wp)             :: uptake

    uptake = 1.5*oxygenation+Ro

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
 !               source_sun         = sun_A*prof%ROC_leaf_G(j)
 !               source_shade       = shd_A*prof%ROC_leaf_G(j)
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
    ! prof%gpp_O2(1:ncl) = prof%dPsdz(1:ncl)*prof%ROC_leaf_G(1:ncl) + prof%rd_O2(1:ncl)! o2 flux of gpp = psn + rd
    ! prof%gpp_O2(1:ncl) = prof%dPsdz_O2(1:ncl) + prof%dRESPdz_O2(1:ncl)! o2 flux of gpp = psn + rd
     prof%dGOPdz_sun(1:ncl) = prof%dLAIdz(1:ncl) * prof%sun_GOP(1:ncl) * solar%prob_beam(1:ncl)
     prof%dGOPdz_shd(1:ncl) = prof%dLAIdz(1:ncl) * prof%shd_GOP(1:ncl) * solar%prob_shd(1:ncl)
     prof%dGOPdz(1:ncl)     = prof%dGOPdz_sun(1:ncl) + prof%dGOPdz_shd(1:ncl)
     prof%ROC_leaf_G(1:ncl) = prof%dGOPdz(1:ncl)/prof%dGPPdz(1:ncl)

    ! photosynthetic O2 per ground level:
     prof%dPsdz_O2_sun(1:ncl) = prof%dLAIdz(1:ncl) * prof%sun_psn_O2(1:ncl) * solar%prob_beam(1:ncl)
     prof%dPsdz_O2_shd(1:ncl) = prof%dLAIdz(1:ncl) * prof%shd_psn_O2(1:ncl) * solar%prob_shd(1:ncl)
     prof%dPsdz_O2(1:ncl)     = prof%dPsdz_O2_sun(1:ncl) + prof%dPsdz_O2_shd(1:ncl)
     prof%ROC_leaf_net(1:ncl) = prof%dPsdz_O2(1:ncl)/prof%dPsdz(1:ncl)
!     print *, prof%shd_GOP
!print *, sum(prof%gpp_O2(1:ncl))
!print *, sum(prof%dGPPdz(1:ncl))
!print *, sum(prof%dPsdz(1:ncl)*prof%ROC_leaf_G(j))
!     prof%gpp_O2(1:ncl) = prof%dGPPdz(1:ncl)* prof%ROC_leaf_G(1:ncl)
!    prof%source_co2(1:ncl) = -prof%dPsdz(1:ncl) + bole%layer(1:ncl) (negative + positive)
!    prof%source_O2(1:ncl) = prof%dPsdz(1:ncl)* prof%ROC_leaf_G(1:ncl)- &! add o2 source/sink from bole respiration
!                            bole%layer(1:ncl)* prof%ROC_bole_air(1:ncl)
     prof%source_O2(1:ncl) = prof%gpp_O2(1:ncl) - prof%dRESPdz_O2(1:ncl) - bole%layer(1:ncl)* prof%ROC_bole_air(1:ncl)
   !  print *, prof%source_O2(33)
   !  print *, prof%dPsdz(33)
   !  print *, bole%layer(33)
     soil%o2 = -soil%respiration_mole*soil%ROC_soil_air ! negative!!!
!     print *, "ground level GPP and GOP:", prof%dGPPdz(23), prof%gpp_O2(23) , prof%dGOPdz(23)!,bole%layer(33), soil%respiration_mole
    ! print *, prof%gpp_O2(33), bole%layer(33)*prof%ROC_bole_air(33), soil%o2
!     o2air = o2*1000 ! convert to ppm
     o2air = o2
     fact%o2 = (mass_air/mass_O2)*met%air_density_mole ! convert from density to mole concentration
!print *,fact%o2
!print *,fact%co2
!     call conc(prof%source_O2, prof%O2_air, input%o2air, soil%o2, fact%o2) ! use input o2 con. Yuan 2018.02.14
     call conc_seperate(prof%source_O2, prof%O2_air, input%o2air, soil%o2, prof%O2_soil, prof%O2_disp, fact%o2)
!    calculate ROC of each layer
     prof%ROC_layer(1:ncl) = abs(prof%source_O2(1:ncl)/prof%source_co2(1:ncl))
    print *, 'ROC = ',prof%ROC_layer(23)
    print *, 'source O2 and CO2',prof%source_O2(23), prof%source_co2(23)
    print *, 'GOP, resp O2 and bole resp',prof%gpp_O2(23), prof%dRESPdz_O2(23), bole%layer(23)
    print *, 'GPP, resp CO2 and bole resp',prof%dGOPdz(23), prof%dRESPdz(23), bole%layer(23)
    print *, 'Psn O2:', prof%dPsdz_O2(23)
    print *, 'Psn CO2', prof%dPsdz(23)
!    call conc(prof%source_co2, prof%co2_air, input%co2air, soil%respiration_mole, fact%co2)

!    call conc(prof%sour13co2, prof%c13cnc, co2_13, soil%resp_13, fact%co2)
    ! this is the o2/(co2) ratio of air in the canopy
    ! prof%R13_12_air(1:ntl) = prof%c13cnc(1:ntl)/prof%co2_air_filter(1:ntl)
  END SUBROUTINE OXYFLUX

END MODULE oxygen
