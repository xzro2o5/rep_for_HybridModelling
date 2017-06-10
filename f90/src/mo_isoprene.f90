MODULE isoprene

  ! This module contains the isoprene flux routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: isoprene_leaf_flux    ! isoprene emission from leaves
  PUBLIC :: isoprene_canopy_flux

  ! ------------------------------------------------------------------

CONTAINS
  
  ! ------------------------------------------------------------------
  FUNCTION isoprene_canopy_flux()

    USE types,      ONLY: prof, solar
    USE setup,      ONLY: ncl
    USE parameters, ONLY: delz, ht

    IMPLICIT NONE

    REAL(wp) :: isoprene_canopy_flux

    INTEGER(i4) :: i
    REAL(wp), DIMENSION(1:ncl) :: zz, hh

    ! leaf level flux per leaf area
    ! This gave always 0 in C-Code: zz = (int)(delz * j);
    ! TODO: forall(i=1:ncl) zz(i) = delz*real(i,wp)
    forall(i=1:ncl) zz(i) = floor(delz*real(i,wp))
    hh(1:ncl) = ht
    prof%iso_sun(1:ncl) = isoprene_leaf_flux(solar%quantum_sun(1:ncl), prof%sun_tleaf_filter(1:ncl), zz(1:ncl), hh(1:ncl))
    prof%iso_shd(1:ncl) = isoprene_leaf_flux(solar%quantum_shd(1:ncl), prof%shd_tleaf_filter(1:ncl), zz(1:ncl), hh(1:ncl))
    ! leaf level flux per ground area
    prof%sun_isopreneflux(1:ncl) = prof%dLAIdz(1:ncl) * (solar%prob_beam(1:ncl) * prof%iso_sun(1:ncl))
    prof%shd_isopreneflux(1:ncl) = prof%dLAIdz(1:ncl) * (solar%prob_shd(1:ncl)  * prof%iso_shd(1:ncl))
    prof%isopreneflux(1:ncl)     = prof%sun_isopreneflux(1:ncl) + prof%shd_isopreneflux(1:ncl)
    ! canopy flux
    isoprene_canopy_flux = sum(prof%isopreneflux(1:ncl))

  END FUNCTION isoprene_canopy_flux


  ! ------------------------------------------------------------------
  ELEMENTAL PURE FUNCTION isoprene_leaf_flux(ppfd, tlf, zzz, hhh)
    ! leaf level algorithms of Guenther/Harley for oak leaves
    USE constants, ONLY: one, TN0, rugc

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: ppfd
    REAL(wp), INTENT(IN) :: tlf
    REAL(wp), INTENT(IN) :: zzz
    REAL(wp), INTENT(IN) :: hhh
    REAL(wp)             :: isoprene_leaf_flux

    REAL(wp) :: t_lfk, light, c_t1, c_t2, isoref

    t_lfk = tlf + TN0
    ! equation from original CANOAK
    light = 0.0027_wp * ppfd * 1.066_wp / sqrt(one + 0.0027_wp * 0.0027_wp * ppfd * ppfd)
    c_t1 = exp(95000._wp * (t_lfk - 303._wp) / (rugc * t_lfk * 303._wp))
    c_t2 = one + exp(230000._wp * (303._wp - 314._wp) / (rugc * 303._wp * 314._wp))
    ! nmol m-2 s-1
    isoref = 45.8_wp * (0.65_wp * zzz / hhh + 0.35_wp)
    ! ! equation after adaptation from Jenifer Funk (Stanford)
    ! light = 0.0015_wp * ppfd * 1.1032_wp / sqrt(one + 0.0015_wp * 0.0015_wp * ppfd * ppfd)
    ! c_t1 = exp(93960._wp * (t_lfk - 301.15_wp) / (rugc * t_lfk * 301.15_wp))
    ! c_t2 = one + exp(540600._wp * (t_lfk - 311.3_wp) / (rugc * t_lfk * 301.15_wp))
    ! ! nmol m-2 s-1
    ! isoref = 32.1_wp * (0.65_wp * zzz / hhh + 0.35_wp)
    isoprene_leaf_flux = isoref * light * c_t1 / c_t2
    
  END FUNCTION isoprene_leaf_flux


END MODULE isoprene
