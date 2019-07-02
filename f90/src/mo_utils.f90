MODULE utils

  ! Helper routines such as calculate correlation coefficient, etc.

  ! Written Jan 2011, Matthias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: corr       ! Pearson''s correlation coefficient
  PUBLIC :: desdt      ! first derivative of saturation vapor pressure function
  PUBLIC :: des2dt     ! second derivative of saturation vapor pressure function
  PUBLIC :: es         ! saturation vapor pressure function
  PUBLIC :: gammaf     ! gamma function
  PUBLIC :: inv_boltz  ! inverse of Boltmann function
  PUBLIC :: lambda     ! latent heat of vaporization
  PUBLIC :: tboltz     ! Boltzmann function for temperature
  PUBLIC :: temp_func  ! Arrenhius function for temperature

  INTERFACE tboltz
     MODULE PROCEDURE tboltz0, tboltz1
  END INTERFACE tboltz

  ! ----------------------------------------------------

CONTAINS

  ! ----------------------------------------------------

  FUNCTION corr(x,y)
    ! modified from pearsn.f90 of numerical recipes in fortran, 2nd edition, 1992
    ! calcs pearsons correlation coefficient
    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN) :: x, y
    REAL(wp) :: corr

    REAL(wp), PARAMETER          :: TINY=1.0e-20_wp
    REAL(wp), DIMENSION(size(x)) :: xt, yt
    REAL(wp)    :: ax, ay, sxx, sxy, syy
    INTEGER(i4) :: n

    n     = size(x)
    ax    = sum(x)/n
    ay    = sum(y)/n
    xt(:) = x(:)-ax
    yt(:) = y(:)-ay
    sxx   = dot_product(xt,xt)
    syy   = dot_product(yt,yt)
    sxy   = dot_product(xt,yt)
    corr  = sxy/(sqrt(sxx*syy)+TINY)

  END FUNCTION corr

  ! ----------------------------------------------------

  ELEMENTAL PURE FUNCTION lambda(tak)
    ! Latent heat of Vaporisation, J kg-1
    USE constants, ONLY: TN0

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: tak
    REAL(wp) :: lambda

    lambda = 3149000._wp - 2370._wp * tak
    ! add heat of fusion for melting ice
    if (tak < TN0) lambda = lambda + 333_wp

  END FUNCTION lambda

  ! ----------------------------------------------------

#ifdef DEBUG
  FUNCTION es(t)
#else
  ELEMENTAL PURE FUNCTION es(t)
#endif
    ! saturation vapor pressure function (mb)
    ! T is temperature in Kelvin
    USE constants, ONLY: zero
#ifdef DEBUG
    USE messages,  ONLY: message
#endif
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: t
    REAL(wp)             :: es

    REAL(wp) :: ess

    ess = zero
    if (t > zero) then
       ess = 54.8781919_wp - 6790.4985_wp / t - 5.02808_wp * log(t)
       ess = exp(ess)
    end if
#ifdef DEBUG
    if (t <= zero) call message('ES: ','input must be > 0')
#endif
    es = ess
    
  END FUNCTION es

  ! ----------------------------------------------------

  FUNCTION desdt(t, latent)
    ! first derivative of es with respect to tk
    ! routine needs to convert es(t) (mb) to Pa
    USE constants, ONLY: rgc1000

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: t
    REAL(wp), INTENT(IN) :: latent
    REAL(wp) :: desdt

    desdt = es(t)*100._wp * latent*18._wp / (rgc1000*t*t)

  END FUNCTION desdt

  ! ----------------------------------------------------

  ELEMENTAL PURE FUNCTION des2dt(t)
    ! The second derivative of the saturation vapor pressure
    ! temperature curve, using the polynomial equation of Paw U
    USE constants, ONLY: TN0, a3en, a4en, a5en

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: t
    REAL(wp) :: des2dt

    REAL(wp) :: tcel

    tcel   = t-TN0
    des2dt = 2._wp*a3en + 6._wp*a4en*tcel + 12._wp*a5en*tcel*tcel

  END FUNCTION des2dt

  ! ----------------------------------------------------

#ifdef DEBUG
  FUNCTION gammaf(x) ! use message for the start
#else
  ELEMENTAL PURE FUNCTION gammaf(x)
#endif
    ! gamma function
    USE constants, ONLY: zero, one, pi2
#ifdef DEBUG
    USE messages,  ONLY: message
#endif

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: gammaf

    REAL(wp) :: gam

    if (x <= zero) then
#ifdef DEBUG
       call message('GAMMAF: ','input must be > 0')
#endif
       gammaf   = zero
    else
       !  gam= (1.0 / (12.0 * x)) + (1.0 / (288.0 * x*x)) - (139.0 / (51840.0 * pow(x,3.0)))
       gam = one + (one/(12._wp*x)) + (one/(288._wp*x*x)) - (139._wp/(51840._wp*x*x*x))
       gammaf = sqrt(pi2/x) * (x**x) * exp(-x) * gam
    end if

  END FUNCTION gammaf

  ! ------------------------------------------------------------------

  FUNCTION inv_boltz(rate, eakin, topt, tl, hkin_in)
    ! calculates the inverse of Boltzmann temperature distribution for photosynthesis
    ! used to get Vcmax and Jmax at Toptimum from Vcmax and Jmax at 25 deg C
    USE constants, ONLY: one, rugc

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

  ! ----------------------------------------------------

  FUNCTION tboltz0(rate, eakin, topt, tl, hkin_in)
    ! Boltzmann temperature distribution for photosynthesis
    USE constants, ONLY: one, rugc

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: rate
    REAL(wp), INTENT(IN) :: eakin
    REAL(wp), INTENT(IN) :: topt
    REAL(wp), INTENT(IN) :: tl
    REAL(wp), INTENT(IN) :: hkin_in
    REAL(wp) :: tboltz0

    REAL(wp) :: dtlopt, prodt, denom, numm

    dtlopt  = tl - topt
    prodt   = rugc * topt * tl
    numm    = rate * hkin_in * exp(eakin * (dtlopt) / (prodt))
    denom   = hkin_in - eakin * (one - exp(hkin_in * (dtlopt) / (prodt)))
    tboltz0 = numm / denom

  END FUNCTION tboltz0


  FUNCTION tboltz1(rate, eakin, topt, tl, hkin_in)
    ! Boltzmann temperature distribution for photosynthesis
    USE constants, ONLY: one, rugc

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN) :: rate
    REAL(wp),               INTENT(IN) :: eakin
    REAL(wp),               INTENT(IN) :: topt
    REAL(wp),               INTENT(IN) :: tl
    REAL(wp),               INTENT(IN) :: hkin_in
    REAL(wp), DIMENSION(1:size(rate))  :: tboltz1

    REAL(wp), DIMENSION(1:size(rate))  :: numm
    REAL(wp) :: dtlopt, prodt, denom, ztmp

    dtlopt  = tl - topt
    prodt   = rugc * topt * tl
    numm    = rate * (hkin_in * exp(eakin * (dtlopt) / (prodt)))
    denom   = hkin_in - eakin * (one - exp(hkin_in * (dtlopt) / (prodt)))
    ztmp    = one / denom
    tboltz1 = numm * ztmp

  END FUNCTION tboltz1

  ! ----------------------------------------------------

  FUNCTION temp_func(rate, eact, tprime, tref, t_lk)
    ! Arhennius temperature function
    USE constants, ONLY: rugc

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: rate
    REAL(wp), INTENT(IN) :: eact
    REAL(wp), INTENT(IN) :: tprime
    REAL(wp), INTENT(IN) :: tref
    REAL(wp), INTENT(IN) :: t_lk
    REAL(wp) :: temp_func

    temp_func = rate * exp(tprime*eact/(tref*rugc*t_lk))

  END FUNCTION temp_func

  ! ----------------------------------------------------

END MODULE utils
