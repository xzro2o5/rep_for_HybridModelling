MODULE isotope_utils

  ! This module contains isotope helper routines such as calculating delta values

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: alpha_equ_co2_h2o    ! CO2-H2O exchange equilibrium fractionation factor
  PUBLIC :: alpha_equ_h2o        ! H2O liquid-vapour equilibrium fractionation factor
  PUBLIC :: alpha_kin_co18o      ! CO18O kinetic fractionation factor
  PUBLIC :: alpha_kin_h2o        ! Vapour H218O and HDO kinetic fractionation factors
  PUBLIC :: delta                ! delta value (rare, abundant, standard)
  PUBLIC :: delta1000            ! delta value in permil
  PUBLIC :: delta_h2o            ! H2O delta value (rare, abundant, species)
  PUBLIC :: delta1000_h2o        ! H2O delta value in permil
  PUBLIC :: invdelta1000         ! Isotope ratio over standard from delta value in permil
  PUBLIC :: invdelta1000_h2o     ! Isotope ratio of H2O from delta value in permil (iso, species)
  PUBLIC :: isorat               ! Isotope ratio with threshold (rare, abundant, collectrare, collectabundant, species)
  PUBLIC :: vtf

  INTERFACE delta_h2o
     MODULE PROCEDURE delta_h2o_00d, delta_h2o_10d, delta_h2o_20d, delta_h2o_11d, delta_h2o_22d
  END INTERFACE
  INTERFACE delta1000_h2o
     MODULE PROCEDURE delta1000_h2o_00d, delta1000_h2o_10d, delta1000_h2o_20d, delta1000_h2o_11d, delta1000_h2o_22d
  END INTERFACE
  INTERFACE invdelta1000_h2o
     MODULE PROCEDURE invdelta1000_h2o_00d, invdelta1000_h2o_10d, invdelta1000_h2o_11d, &
          invdelta1000_h2o_20d, invdelta1000_h2o_21d, invdelta1000_h2o_22d
  END INTERFACE
  INTERFACE isorat
     MODULE PROCEDURE isorat_00d, isorat_10d, isorat_20d
  END INTERFACE

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION alpha_equ_co2_h2o(T, species)
    ! Calculates the equilibrium fractionation factor 
    !   of CO2 equilibrated with water
    ! Brenninkmeier et al., Isotope Geoscience, 1, 181-190, 1983
    !  T in [KELVIN]
    !  species: 2: CO18O; else: no fractionation, return 1.
    USE constants, ONLY: one, e3, TN0

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: T
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp)                :: alpha_equ_co2_h2o

    if (species == 2) then
       if (T > TN0) then
          alpha_equ_co2_h2o = one + (17604._wp/T - 17.93_wp)*e3
       else ! ice
          alpha_equ_co2_h2o = one
       end if
    else
       alpha_equ_co2_h2o = one
    end if

  END FUNCTION alpha_equ_co2_h2o

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION alpha_equ_h2o(temp, species)
    ! Calculates fractionation factors during liquid-water vapour
    ! equilibration at temperature temp [K]
    ! species: 2: H218O 3: HDO else: no fractionation, return 1
    ! Reference: Majoube (1971)
    USE constants, ONLY: zero, one

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: temp
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp)                :: alpha_equ_h2o

    REAL(wp) :: a, b, c

    select case (species)
    case (2)
       a = +1.137e+3_wp
       b = -4.156e-1_wp
       c = -2.067e-3_wp
    case (3)
       a = +2.4844e+4_wp
       b = -7.6248e+1_wp
       c = +5.261e-2_wp
    case default
       a = zero
       b = zero
       c = zero
    end select
    alpha_equ_h2o = one/exp((a/temp + b)/temp +c)

  END FUNCTION alpha_equ_h2o

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION alpha_kin_co18o(rs, rb, species)
    ! Calculates kinetic fractionation factor for co18o diffusion
    !   through a pure molecular diffusion resistance (rs)
    !   and a boundary layer resistance (rb).
    ! Assume 2/3 power law for boundary layer.
    ! species: 2: CO18O; else: no fractionation, return 1.
    ! Reference: Farquhar & Lloyd (1993)
    USE constants, ONLY: zero, twothird, one

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: rs
    REAL(wp),    INTENT(IN) :: rb
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp)                :: alpha_kin_co18o

    REAL(wp) :: alpha_k

    select case (species)
    case (2)
       alpha_k = 0.9912_wp
    case default
       alpha_k = one
    end select
    !
    if (abs(rs+rb) > tiny(one)) then
       alpha_kin_co18o = (rs*alpha_k + rb*alpha_k**twothird) / (rs+rb)
    else
       alpha_kin_co18o = zero
    end if

  END FUNCTION alpha_kin_co18o

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION alpha_kin_h2o(rs, rb, species, merlivat)
    ! Calculates kinetic fractionation factor for water vapour diffusion
    ! through a pure molecular diffusion resistance (rs)
    ! and a boundary layer resistance (rb).
    ! Assume 2/3 power law for boundary layer.
    ! species: 2: H218O 3: HDO else: no fractionation, return 1
    ! Reference: Farquhar & Lloyd (1993) - weighting
    !            Merlivat (1978) - molecular diffusision fractionation factors
    !            Cappa et al. (2003) - dito
    USE constants, ONLY: zero, twothird, one

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: rs
    REAL(wp),    INTENT(IN) :: rb
    INTEGER(i4), INTENT(IN) :: species
    INTEGER(i4), INTENT(IN) :: merlivat
    REAL(wp)                :: alpha_kin_h2o

    REAL(wp) :: alpha_k

    if (merlivat == 1) then
       select case (species)
       case (2)
          alpha_k = 0.9727_wp
       case (3)
          alpha_k = 0.9755_wp
       case default
          alpha_k = one
       end select
    else
       select case (species)
       case (2)
          alpha_k = 0.9691_wp
       case (3)
          alpha_k = 0.9839_wp
       case default
          alpha_k = one
       end select
    end if
    !
    if (abs(rs+rb) > tiny(one)) then
       alpha_kin_h2o = (alpha_k*rs + rb*alpha_k**twothird) / (rs+rb)
    else
       alpha_kin_h2o = zero
    end if

  END FUNCTION alpha_kin_h2o

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION delta(rare, abundant, standard)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: rare
    REAL(wp), INTENT(IN) :: abundant
    REAL(wp), INTENT(IN) :: standard
    REAL(wp)             :: delta

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*standard)
       delta = rare*ztmp - one
    else
       delta = -one
    end if

  END FUNCTION delta

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION delta1000(rare, abundant, standard)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: rare
    REAL(wp), INTENT(IN) :: abundant
    REAL(wp), INTENT(IN) :: standard
    REAL(wp)             :: delta1000

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*standard)
       delta1000 = rare*ztmp - one
    else
       delta1000 = -one
    end if
    delta1000 = delta1000 * 1000._wp

  END FUNCTION delta1000

  ! ------------------------------------------------------------------

  FUNCTION delta_h2o_00d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: rare
    REAL(wp),    INTENT(IN) :: abundant
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp) :: delta_h2o_00d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta_h2o_00d = rare*ztmp - one
    else
       delta_h2o_00d = -one
    end if

  END FUNCTION delta_h2o_00d


  FUNCTION delta_h2o_10d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: rare
    REAL(wp),                  INTENT(IN) :: abundant
    INTEGER(i4),               INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare)) :: delta_h2o_10d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta_h2o_10d = rare*ztmp - one
    else
       delta_h2o_10d = -one
    end if

  END FUNCTION delta_h2o_10d


  FUNCTION delta_h2o_20d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: rare
    REAL(wp),                    INTENT(IN) :: abundant
    INTEGER(i4),                 INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: delta_h2o_20d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta_h2o_20d = rare*ztmp - one
    else
       delta_h2o_20d = -one
    end if

  END FUNCTION delta_h2o_20d


  FUNCTION delta_h2o_11d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: rare
    REAL(wp),    DIMENSION(:), INTENT(IN) :: abundant
    INTEGER(i4),               INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare)) :: delta_h2o_11d

    REAL(wp), DIMENSION(1:size(rare)) :: ztmp

    where (abs(abundant) > tiny(one))
       ztmp = one / (abundant*wiso%vsmow(species))
       delta_h2o_11d = rare*ztmp - one
    elsewhere
       delta_h2o_11d = -one
    end where

  END FUNCTION delta_h2o_11d


  FUNCTION delta_h2o_22d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: rare
    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: abundant
    INTEGER(i4),                 INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: delta_h2o_22d

    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: ztmp

    where (abs(abundant) > tiny(one))
       ztmp = one / (abundant*wiso%vsmow(species))
       delta_h2o_22d = rare*ztmp - one
    elsewhere
       delta_h2o_22d = -one
    end where

  END FUNCTION delta_h2o_22d

  ! ------------------------------------------------------------------

  FUNCTION delta1000_h2o_00d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: rare
    REAL(wp),    INTENT(IN) :: abundant
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp) :: delta1000_h2o_00d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta1000_h2o_00d = rare*ztmp - one
    else
       delta1000_h2o_00d = -one
    end if
    delta1000_h2o_00d = delta1000_h2o_00d * 1000._wp

  END FUNCTION delta1000_h2o_00d


  FUNCTION delta1000_h2o_10d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: rare
    REAL(wp),                  INTENT(IN) :: abundant
    INTEGER(i4),               INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare)) :: delta1000_h2o_10d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta1000_h2o_10d = rare*ztmp - one
    else
       delta1000_h2o_10d = -one
    end if
    delta1000_h2o_10d = delta1000_h2o_10d * 1000._wp

  END FUNCTION delta1000_h2o_10d


  FUNCTION delta1000_h2o_20d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: rare
    REAL(wp),                    INTENT(IN) :: abundant
    INTEGER(i4),                 INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: delta1000_h2o_20d

    REAL(wp) :: ztmp

    if (abs(abundant) > tiny(one)) then
       ztmp = one / (abundant*wiso%vsmow(species))
       delta1000_h2o_20d = rare*ztmp - one
    else
       delta1000_h2o_20d = -one
    end if
    delta1000_h2o_20d = delta1000_h2o_20d * 1000._wp

  END FUNCTION delta1000_h2o_20d


  FUNCTION delta1000_h2o_11d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: rare
    REAL(wp),    DIMENSION(:), INTENT(IN) :: abundant
    INTEGER(i4),               INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare)) :: delta1000_h2o_11d

    REAL(wp), DIMENSION(1:size(rare)) :: ztmp

    where (abs(abundant) > tiny(one))
       ztmp = one / (abundant*wiso%vsmow(species))
       delta1000_h2o_11d = rare*ztmp - one
    elsewhere
       delta1000_h2o_11d = -one
    end where
    delta1000_h2o_11d = delta1000_h2o_11d * 1000._wp

  END FUNCTION delta1000_h2o_11d


  FUNCTION delta1000_h2o_22d(rare, abundant, species)
    ! Calculates delta value for isotopes expressed in per mil
    USE constants, ONLY: one
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: rare
    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: abundant
    INTEGER(i4),                 INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: delta1000_h2o_22d

    REAL(wp), DIMENSION(1:size(rare,1),1:size(rare,2)) :: ztmp

    where (abs(abundant) > tiny(one))
       ztmp = one / (abundant*wiso%vsmow(species))
       delta1000_h2o_22d = rare*ztmp - one
    elsewhere
       delta1000_h2o_22d = -one
    end where
    delta1000_h2o_22d = delta1000_h2o_22d * 1000._wp

  END FUNCTION delta1000_h2o_22d

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION invdelta1000(iso)
    ! Calculates ratio of isotope ratio to standard from delta value expressed in per mil
    USE constants, ONLY: one, e3

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: iso
    REAL(wp) :: invdelta1000

    invdelta1000 = iso*e3 + one

  END FUNCTION invdelta1000

  ! ------------------------------------------------------------------

  FUNCTION invdelta1000_h2o_00d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: iso
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp)                :: invdelta1000_h2o_00d

    invdelta1000_h2o_00d = (iso*e3+one)*wiso%vsmow(species)

  END FUNCTION invdelta1000_h2o_00d


  FUNCTION invdelta1000_h2o_10d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: iso
    INTEGER(i4),               INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(iso))      :: invdelta1000_h2o_10d

    invdelta1000_h2o_10d = (iso*e3+one)*wiso%vsmow(species)

  END FUNCTION invdelta1000_h2o_10d


  FUNCTION invdelta1000_h2o_11d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso
    USE finishes,  ONLY: finish

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:), INTENT(IN) :: iso
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(iso))      :: invdelta1000_h2o_11d

    if (size(species) /= size(iso)) call finish('invdelta1000_h2o_11d','size(spezies) /= size(iso)')
    invdelta1000_h2o_11d = (iso*e3+one)*wiso%vsmow(species)

  END FUNCTION invdelta1000_h2o_11d


  FUNCTION invdelta1000_h2o_20d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: iso
    INTEGER(i4),                 INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(iso,1),1:size(iso,2)) :: invdelta1000_h2o_20d

    invdelta1000_h2o_20d = (iso*e3+one)*wiso%vsmow(species)

  END FUNCTION invdelta1000_h2o_20d


  FUNCTION invdelta1000_h2o_21d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso
    USE finishes,  ONLY: finish

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: iso
    INTEGER(i4), DIMENSION(:),   INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(iso,1),1:size(iso,2)) :: invdelta1000_h2o_21d

    if (size(species,1) /= size(iso,2)) call finish('invdelta1000_h2o_21d','size(spezies) /= size(iso,2)')
    invdelta1000_h2o_21d = (iso*e3+one)*spread(wiso%vsmow(species),1,size(iso,1))

  END FUNCTION invdelta1000_h2o_21d


  FUNCTION invdelta1000_h2o_22d(iso, species)
    ! Calculates ratio from water isotope value expressed in per mil
    USE constants, ONLY: one, e3
    USE types,     ONLY: wiso
    USE finishes,  ONLY: finish

    IMPLICIT NONE

    REAL(wp),    DIMENSION(:,:), INTENT(IN) :: iso
    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: species
    REAL(wp), DIMENSION(1:size(iso,1),1:size(iso,2)) :: invdelta1000_h2o_22d

    INTEGER(i4) :: i, n1, n2
    REAL(wp), DIMENSION(1:size(iso,1),1:size(iso,2)) :: ztmp

    n1 = size(iso,1)
    n2 = size(iso,2)
    if (size(species,1) /= n1) call finish('invdelta1000_h2o_22d','size(spezies,1) /= size(iso,1)')
    if (size(species,2) /= n2) call finish('invdelta1000_h2o_22d','size(spezies,2) /= size(iso,2)')
    do i=1, n1
       ztmp(i,1:n2) = wiso%vsmow(species(i,1:n2))
    end do
    invdelta1000_h2o_22d = (iso*e3+one) * ztmp

  END FUNCTION invdelta1000_h2o_22d

  ! ------------------------------------------------------------------

  FUNCTION isorat_00d(qi, q, lostqi, lostq)
    ! returns isotope ratio qi/q
    ! checks for underflow. If abs(q)<1e-15 -> lostqi=lostqi+qi and qi=q=0
    ! double precision allows for about 15 digits of precision -> 1e-15
    USE constants,    ONLY: zero, one
#ifdef DEBUG
    USE messages,     ONLY: message
    USE types,        ONLY: time
    USE string_utils, ONLY: num2str
#endif
    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: qi
    REAL(wp), INTENT(INOUT) :: q
    REAL(wp), INTENT(INOUT) :: lostqi
    REAL(wp), INTENT(INOUT) :: lostq
    REAL(wp)                :: isorat_00d

    if (abs(q) >= epsilon(one)) then
       isorat_00d = qi / q
    else
#ifdef DEBUG
       if (abs(q) > tiny(one) .or. abs(qi) > tiny(one)) then
          call message('ISORAT: ', num2str(qi), num2str(q), num2str(time%daytime))
       end if
#endif
       lostqi     = lostqi + qi
       lostq      = lostq + q
       qi         = zero
       q          = zero
       isorat_00d = zero
    end if

  END FUNCTION isorat_00d


  FUNCTION isorat_10d(qi, q, lostqi, lostq)
    ! returns isotope ratio qi/q
    ! checks for underflow. If abs(q)<1e-15 -> lostqi=lostqi+qi and qi=q=0
    ! double precision allows for about 15 digits of precision -> 1e-15
    USE constants,    ONLY: zero, one

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(INOUT) :: qi     ! nwiso
    REAL(wp),               INTENT(INOUT) :: q
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: lostqi ! nwiso
    REAL(wp),               INTENT(INOUT) :: lostq
    REAL(wp), DIMENSION(1:size(qi))       :: isorat_10d

    REAL(wp) :: ztmp

    if (abs(q) >= epsilon(one)) then
       ztmp       = one / q
       isorat_10d = qi * ztmp
    else
       lostqi     = lostqi + qi
       lostq      = lostq + q
       qi         = zero
       q          = zero
       isorat_10d = zero
    end if

  END FUNCTION isorat_10d

  FUNCTION isorat_20d(qi, q, lostqi, lostq)
    ! returns isotope ratio qi/q
    ! checks for underflow. If abs(q)<1e-15 -> lostqi=lostqi+qi and qi=q=0
    ! double precision allows for about 15 digits of precision -> 1e-15
    USE constants,    ONLY: zero, one

    IMPLICIT NONE

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) :: qi     ! nsoil,nwiso
    REAL(wp), DIMENSION(:),   INTENT(INOUT) :: q      ! nsoil
    REAL(wp), DIMENSION(:),   INTENT(INOUT) :: lostqi !       nwiso
    REAL(wp),                 INTENT(INOUT) :: lostq  ! 
    REAL(wp), DIMENSION(1:size(qi,1),1:size(qi,2)) :: isorat_20d

    INTEGER(i4) :: n1, n2
    REAL(wp), DIMENSION(1:size(qi,1),1:size(qi,2)) :: ztmp
    REAL(wp), DIMENSION(1:size(qi,1),1:size(qi,2)) :: zmask

    n1 = size(qi,1)
    n2 = size(qi,2)
    ztmp  = spread(q,2,n2)
    where (abs(ztmp) >= epsilon(one))
       isorat_20d = qi / ztmp
       zmask      = one
    elsewhere
       isorat_20d = zero
       zmask      = zero
    end where
    lostqi = lostqi + sum(qi*(one-zmask),1)
    lostq  = lostq  + sum(q*(one-zmask(1:n1,1)))
    qi = qi * zmask
    q  = q  * zmask(1:n1,1)
    isorat_20d = isorat_20d * zmask

  END FUNCTION isorat_20d

  ! ------------------------------------------------------------------

  ELEMENTAL PURE FUNCTION vtf(T, species)
    ! Vogel-Tamman-Fulcher relationship: D(T)=D0*exp(-B0/(T-T0))
    ! T [K]
    USE constants, ONLY: one

    IMPLICIT NONE

    REAL(wp),    INTENT(IN) :: T ! [K]
    INTEGER(i4), INTENT(IN) :: species
    REAL(wp)                :: vtf

    REAL(wp), PARAMETER :: A0 = 100e-09_wp
    REAL(wp), PARAMETER :: A1 = 577._wp
    REAL(wp), PARAMETER :: A2 = 145._wp
    REAL(wp) :: c1

    select case (species)
    case (2)
       c1 = sqrt(18._wp/20._wp*(20._wp+18._wp)/(18._wp+18._wp))
    case (3)
       c1 = sqrt(18._wp/19._wp*(19._wp+18._wp)/(18._wp+18._wp))
    case default
       c1 = one
    end select

    vtf = c1*A0*exp(-A1/(T-A2))

  END FUNCTION vtf

END MODULE isotope_utils
