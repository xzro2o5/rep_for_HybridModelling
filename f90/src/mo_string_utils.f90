MODULE string_utils

  ! From Echam5, (C) MPI-MET, Hamburg, Germany

  ! This module holds string conversion utilities

  ! Modified Jan 2011, Matthias Cuntz - Adapted for Canveg
  !                                     renamed from util_string, rename __PGI
  !                                     moved int2string etc. from mo_exception.f90, renamed them
  !                                     format for num2str

  USE kinds, ONLY: i4, i8, sp, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tolower   ! Conversion   : 'ABCXYZ' -> 'abcxyz'   
  PUBLIC :: toupper   ! Conversion   : 'abcxyz' -> 'ABCXYZ'
  PUBLIC :: separator ! Format string: '-----...-----'
  PUBLIC :: num2str

  INTERFACE num2str
     MODULE PROCEDURE i42str, i82str, sp2str, dp2str, log2str
  END INTERFACE

  !CHARACTER(len=*), PARAMETER :: separator = '(70("-"))'
  CHARACTER(len=*), PARAMETER :: separator = repeat('-',70)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  FUNCTION tolower (upper)
    ! Conversion: Uppercase -> Lowercase
    CHARACTER(LEN=*)              ,INTENT(in) :: upper
    CHARACTER(LEN=LEN_TRIM(upper))            :: tolower

    INTEGER            :: i
    INTEGER ,PARAMETER :: idel = ICHAR('a')-ICHAR('A')

    DO i=1,LEN_TRIM(upper)
      IF (ICHAR(upper(i:i)) >= ICHAR('A') .AND. &
          ICHAR(upper(i:i)) <= ICHAR('Z')) THEN
        tolower(i:i) = CHAR( ICHAR(upper(i:i)) + idel )
      ELSE
        tolower(i:i) = upper(i:i)
      END IF
    END DO

  END FUNCTION tolower

  ! ------------------------------------------------------------------

  FUNCTION toupper (lower)
    ! Conversion: Lowercase -> Uppercase
    CHARACTER(LEN=*)              ,INTENT(in) :: lower
    CHARACTER(LEN=LEN_TRIM(lower))            :: toupper

    INTEGER            :: i
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')

    DO i=1,LEN_TRIM(lower)
      IF (ICHAR(lower(i:i)) >= ICHAR('a') .AND. &
          ICHAR(lower(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lower(i:i)) + idel )
      ELSE
        toupper(i:i) = lower(i:i)
      END IF
    END DO

  END FUNCTION toupper

  ! ------------------------------------------------------------------

  PURE FUNCTION i42str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i4),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: i42str

    if (present(form)) then
       write(i42str,form) nn
    else
       write(i42str,'(I10)') nn
    end if
    !i42str = adjustl(i42str)

  END FUNCTION i42str


  PURE FUNCTION i82str(nn,form)
    ! returns integer nn as a string (often needed in printing messages)
    IMPLICIT NONE
    INTEGER(i8),      INTENT(IN)           :: nn
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=20) :: i82str

    if (present(form)) then
       write(i82str,form) nn
    else
       write(i82str,'(I20)') nn
    end if
    !i82str = adjustl(i82str)

  END FUNCTION i82str


  PURE FUNCTION sp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(sp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: sp2str

    if (present(form)) then
       write(sp2str,form) rr
    else
       write(sp2str,'(G32.5)') rr
    end if
    !sp2str = adjustl(sp2str)

  END FUNCTION sp2str


  PURE FUNCTION dp2str(rr,form)
    ! returns real rr as a string (often needed in printing messages)
    IMPLICIT NONE
    REAL(dp),         INTENT(IN)           :: rr
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=32) :: dp2str

    if (present(form)) then
       write(dp2str,form) rr
    else
       write(dp2str,'(G32.5)') rr
    end if
    !dp2str = adjustl(dp2str)

  END FUNCTION dp2str


  PURE FUNCTION log2str(ll,form)
    ! returns logical ll as a string (often needed in printing messages)
    IMPLICIT NONE
    LOGICAL,          INTENT(in)           :: ll
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: form
    CHARACTER(len=10) :: log2str

    if (present(form)) then
       write(log2str,form) ll
    else
       write(log2str,'(L10)') ll
    end if
    !log2str = adjustl(log2str)

  END FUNCTION log2str

END MODULE string_utils
 
