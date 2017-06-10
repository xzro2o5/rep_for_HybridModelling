MODULE messages

  ! This module supplies routines to write out text

  ! Written Jul 2011, Matthias Cuntz - Inspired from Echam5 mo_exception.f90

  USE constants, ONLY: nout

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text    ! dummy string to use in subroutines
  PUBLIC :: message         ! versatile routine to write out strings in file or on screen

  CHARACTER(len=1024) :: message_text = ''

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  SUBROUTINE message(t01, t02, t03, t04, t05, t06, t07, t08, t09, t10, uni, advance)

    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t01
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t02
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t03
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t04
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t05
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t06
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t07
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t08
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t09
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t10
    INTEGER,          INTENT(IN), OPTIONAL :: uni
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: advance

    INTEGER              :: iout
    CHARACTER(len=32000) :: out
    CHARACTER(len=3)     :: iadv


    if (present(uni)) then
       iout = uni
    else
       iout = nout
    end if
    if (present(advance)) then
       iadv = ''
       iadv(1:min(len(advance),3)) = advance(1:min(len(advance),3))
    else
       iadv = 'yes'
    end if

    out = ''
    ! start from back so that trim does not remove user desired blanks
    if (present(t10)) out = t10//trim(out)
    if (present(t09)) out = t09//trim(out)
    if (present(t08)) out = t08//trim(out)
    if (present(t07)) out = t07//trim(out)
    if (present(t06)) out = t06//trim(out)
    if (present(t05)) out = t05//trim(out)
    if (present(t04)) out = t04//trim(out)
    if (present(t03)) out = t03//trim(out)
    if (present(t02)) out = t02//trim(out)
    if (present(t01)) out = t01//trim(out)

    ! output at least one space otherwise some compilers get confused on Mac (empty assembler statement)
    write(iout,'(a)',advance=iadv) trim(out)//' '

  END SUBROUTINE message

END MODULE messages
