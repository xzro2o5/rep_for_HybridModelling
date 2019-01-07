MODULE kinds

  ! From Echam5, (C) MPI-MET, Hamburg, Germany

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15
  !                   exponent  = 37   exponent  =  307
  !
  ! Most likely these are the only possible models.

  ! Author
  !     L. Kornblueh, MPI, August 2001, added working precision and comments

  ! Modified Jan 2011, Matthias Cuntz - Adapted for Canveg
  !                                   - Private/Public

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sp, dp, rp, wp, i4, i8

  ! Floating point section
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER :: rp = SELECTED_REAL_KIND(24,307)
  INTEGER, PARAMETER :: wp = dp   ! working precision

  ! Integer section
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

END MODULE kinds
