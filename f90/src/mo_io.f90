MODULE io

  ! This module contains the wrapper routines for text and netcdf input/output of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE io_text, ONLY: lastin

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lastin

  PUBLIC :: close_files
  PUBLIC :: open_files
  PUBLIC :: read_disp
  PUBLIC :: read_in
  PUBLIC :: skip_input
  PUBLIC :: write_daily
  PUBLIC :: write_profiles
  PUBLIC :: write_output

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE close_files()
    ! Calls the routines to open input and output files
    USE parameters, ONLY: netcdf_in, netcdf_out
    USE io_text,    ONLY: close_in_text, close_out_text

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call close_in_text()
    else
       continue
    end if

    if (netcdf_out == 0) then
       call close_out_text()
    else
       continue
    end if

  END SUBROUTINE close_files


  ! ------------------------------------------------------------------
  SUBROUTINE open_files()
    ! Calls the routines to open input and output files
    USE parameters, ONLY: netcdf_in, netcdf_out
    USE io_text,    ONLY: open_in_text, open_out_text

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call open_in_text()
    else
       continue
    end if

    if (netcdf_out == 0) then
       call open_out_text()
    else
       continue
    end if

  END SUBROUTINE open_files


  ! ------------------------------------------------------------------
  SUBROUTINE read_disp()
    ! Input data on Thomson dispersion matrix that was computed offline with
    ! MOVOAK.C, Dij (s m-1)
    USE parameters, ONLY: netcdf_disp
    USE io_text,    ONLY: read_disp_text

    IMPLICIT NONE

    if (netcdf_disp == 0) then
       call read_disp_text()
    else
       continue
    end if

  END SUBROUTINE read_disp


  ! ------------------------------------------------------------------
  SUBROUTINE read_in()
    ! input data and check for bad data
    USE parameters, ONLY: netcdf_in
    USE io_text,    ONLY: read_in_text

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call read_in_text()
    else
       continue
    end if

  END SUBROUTINE read_in


  ! ------------------------------------------------------------------
  SUBROUTINE skip_input()
    ! Skip input lines or set netcdf data record
    USE parameters, ONLY: netcdf_in
    USE io_text,    ONLY: skip_in_text

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call skip_in_text()
    else
       continue
    end if

  END SUBROUTINE skip_input


  ! ------------------------------------------------------------------
  SUBROUTINE write_daily()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_daily_text

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_daily_text()
    else
       continue
    end if
  !  print *, "called"

  END SUBROUTINE write_daily


  ! ------------------------------------------------------------------
  SUBROUTINE write_profiles()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_prof_text

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_prof_text()
    else
       continue
    end if

  END SUBROUTINE write_profiles


  ! ------------------------------------------------------------------
  SUBROUTINE write_output()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_out_text

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_out_text()
    else
       continue
    end if

  END SUBROUTINE write_output


END MODULE io
