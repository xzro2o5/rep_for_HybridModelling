MODULE nml

  ! Adapted from Echam5, (C) MPI-MET, Hamburg, Germany

  ! Reading and positioning namelist file
  ! Author:
  !     L. Kornblueh, MPI, March 2001, original source

  ! Modified Jan 2011, Matthias Cuntz - Adapted for Canveg
  !                                     compatible with gfortran <= version 4.3
  !                                     all integer(i4)
  !                                     quiet

  USE kinds,        ONLY: i4
  USE string_utils, ONLY: tolower
  USE messages,     ONLY: message, message_text
  USE finishes,     ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: open_nml                                      ! subroutine: open namelist file
  PUBLIC :: close_nml                                     ! subroutine: close namelist file
  PUBLIC :: position_nml                                  ! subroutine: position namelist file
  PUBLIC :: nnml                                          ! namelist unit
  PUBLIC :: POSITIONED, MISSING, LENGTH_ERROR, READ_ERROR ! return values from position_nml

  ! return values of function 'position_nml'
  INTEGER(i4), PARAMETER :: POSITIONED   =  0 ! file pointer set to namelist group
  INTEGER(i4), PARAMETER :: MISSING      =  1 ! namelist group is missing
  INTEGER(i4), PARAMETER :: LENGTH_ERROR =  2 !
  INTEGER(i4), PARAMETER :: READ_ERROR   =  3 !
  ! default namelist unit
  INTEGER,     SAVE      :: nnml         = -1

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE open_nml(file, unit, quiet)

    CHARACTER(len=*), INTENT(IN) :: file
    INTEGER         , INTENT(IN) :: unit
    INTEGER         , INTENT(IN), OPTIONAL :: quiet
    INTEGER :: istat

    nnml = unit
    if (.not. present(quiet)) then
       message_text = '    This is namelist '//file
       CALL message(trim(message_text))
    end if
    OPEN (nnml, file=file, iostat=istat, status='old', action='read', &
         delim='apostrophe')

    IF (istat /= 0) THEN
       message_text = 'Could not open namelist file '//TRIM(file)
      CALL finish('OPEN_NML',trim(message_text))
    END IF

  END SUBROUTINE open_nml  

  ! ------------------------------------------------------------------
  SUBROUTINE close_nml ()

    INTEGER :: istat

    IF (nnml == -1) THEN
      CALL finish('CLOSE_NML','No namelist file opened.')
    END IF

    CLOSE(nnml, IOSTAT=istat)

    IF (istat /= 0) THEN
      CALL finish('CLOSE_NML','Could not close namelist file.')
    END IF

    nnml = -1

  END SUBROUTINE close_nml

  ! ------------------------------------------------------------------
  SUBROUTINE position_nml (name, unit, REWIND, status)

    ! position_nml - position namelist file for reading
    !
    ! Purpose:
    !
    ! To position namelist file at correct place for reading
    ! namelist /name/ (case independent). 
    !

    CHARACTER(len=*), INTENT(in)            :: name   ! namelist group name
    INTEGER,          INTENT(in)  ,OPTIONAL :: unit   ! file unit number
    LOGICAL,          INTENT(in)  ,OPTIONAL :: REWIND ! default: true
    INTEGER(i4),      INTENT(out) ,OPTIONAL :: status ! error return value

    CHARACTER(len=256) :: yline    ! line read
    CHARACTER(len=256) :: test     ! uppercase namelist group name
    INTEGER(i4)        :: stat     ! local copy of status variable
    INTEGER            :: ios      ! status variable from read operation
    LOGICAL            :: lrew     ! local copy of rewind flag
    INTEGER(i4)        :: iunit    ! local copy of unit number
    INTEGER(i4)        :: len_name ! length of requested namelist group name
    CHARACTER          :: ytest    ! character to test for delimiter
    CHARACTER(len=12)  :: code     ! error code printed
    INTEGER(i4)        :: ind      ! index from index routine
    INTEGER(i4)        :: indc     ! index of comment character (!)

    lrew  = .TRUE.
    IF (PRESENT(REWIND)) lrew  = REWIND
    iunit =  nnml
    IF (PRESENT(unit)) iunit = unit   
    stat  =  MISSING
    code  = 'MISSING'

    len_name = LEN_TRIM(name)

    IF (len_name > LEN(test)) THEN
       stat =  LENGTH_ERROR
       code = 'LENGTH_ERROR'
    END IF

    test = '&'//tolower(name)

    ! Reposition file at beginning:

    IF (lrew) REWIND(iunit)

    ! Search start of namelist

    DO
       IF (stat /= MISSING) EXIT

       yline = ' '

       READ (iunit,'(a)',IOSTAT=ios) yline
       IF (ios < 0) THEN
          EXIT  ! MISSING
       ELSE IF (ios > 0) THEN
          stat =  READ_ERROR
          code = 'READ_ERROR'
          EXIT
       END IF

       yline = tolower(yline)

       ind = INDEX(yline,TRIM(test))

       IF (ind == 0) CYCLE

       indc = INDEX(yline,'!')

       IF (indc > 0 .AND. indc < ind) CYCLE

       ! test for delimiter

       ytest = yline(ind+len_name+1:ind+len_name+1)

       IF ( (LGE(ytest,'0') .AND. LLE(ytest,'9')) .OR. &
            (LGE(ytest,'a') .AND. LLE(ytest,'z')) .OR. &
            ytest == '_'                         .OR. &
            (LGE(ytest,'A') .AND. LLE(ytest,'Z'))) THEN
          CYCLE
       ELSE 
          stat = POSITIONED
          BACKSPACE(iunit)
          EXIT
       END IF
    END DO

    IF (PRESENT(status)) status = stat
    SELECT CASE (stat)
    CASE (POSITIONED)
       RETURN
    CASE (MISSING)
       IF (PRESENT(status)) RETURN
    END SELECT

    ! Error if it reaches here
    message_text = 'namelist /'//TRIM(name)//'/ '//code
    CALL finish('POSITION_NML',message_text)

  END SUBROUTINE position_nml

END MODULE nml
