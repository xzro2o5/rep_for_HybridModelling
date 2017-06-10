MODULE finishes

  ! From Echam5, (C) MPI-MET, Hamburg, Germany

  ! This module supplies routines to finish gracefully

  ! Modified Jan 2011, Matthias Cuntz - Adapted for Canveg

  USE constants,    ONLY: nerr
  USE string_utils, ONLY: separator

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finish     ! Write out error message and stop program

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE finish(name, text)

    CHARACTER(len=*), INTENT(IN)           :: name
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: text

    WRITE (nerr,'(a)') separator

    IF (PRESENT(text)) THEN
      WRITE (nerr,'(a,a,a)') name, ': ', text
    ELSE
      WRITE (nerr,'(a)') name
    END IF

    WRITE (nerr,'(a)') separator

    STOP

  END SUBROUTINE finish

END MODULE finishes
