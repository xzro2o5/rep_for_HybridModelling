MODULE io

  ! This module contains the wrapper routines for text and netcdf input/output of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE io_text, ONLY: lastin

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lastin

  PUBLIC :: create_dir ! Yuan 2018.01.22
  PUBLIC :: close_files
  PUBLIC :: open_files
  PUBLIC :: read_disp
  PUBLIC :: read_in
  PUBLIC :: skip_input
  PUBLIC :: write_daily
  PUBLIC :: write_profiles
  PUBLIC :: write_output
  PUBLIC :: copy_code ! Yuan 2018.05.07

  ! ------------------------------------------------------------------

CONTAINS
  ! ------------------------------------------------------------------
  SUBROUTINE create_dir()
    ! create
    USE io_text, ONLY: create_text_dir

    IMPLICIT NONE

    call create_text_dir()

  END SUBROUTINE create_dir
  ! ------------------------------------------------------------------
  SUBROUTINE copy_code()
   USE io_text, ONLY: copy_text_code

   IMPLICIT NONE

   call copy_text_code()

  END SUBROUTINE copy_code
  ! ------------------------------------------------------------------
  SUBROUTINE close_files()
    ! Calls the routines to open input and output files
    USE parameters, ONLY: netcdf_in, netcdf_out
    USE io_text,    ONLY: close_text_in, close_text_out

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call close_text_in()
    else
       continue
    end if

    if (netcdf_out == 0) then
       call close_text_out()
    else
       continue
    end if

  END SUBROUTINE close_files


  ! ------------------------------------------------------------------
  SUBROUTINE open_files()
    ! Calls the routines to open input and output files
    USE parameters, ONLY: netcdf_in, netcdf_out
    USE io_text,    ONLY: open_text_in, open_text_out

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call open_text_in()
    else
       continue
    end if

    if (netcdf_out == 0) then
       call open_text_out()
    else
       continue
    end if

  END SUBROUTINE open_files


  ! ------------------------------------------------------------------
  SUBROUTINE read_disp()
    ! Input data on Thomson dispersion matrix that was computed offline with
    ! MOVOAK.C, Dij (s m-1)
    USE parameters, ONLY: netcdf_disp
    USE io_text,    ONLY: read_text_disp

    IMPLICIT NONE

    if (netcdf_disp == 0) then
       call read_text_disp()
    else
       continue
    end if

  END SUBROUTINE read_disp


  ! ------------------------------------------------------------------
  SUBROUTINE read_in()
    ! input data and check for bad data
    USE parameters, ONLY: netcdf_in
    USE io_text,    ONLY: read_text_in

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call read_text_in()
    else
       continue
    end if

  END SUBROUTINE read_in


  ! ------------------------------------------------------------------
  SUBROUTINE skip_input()
    ! Skip input lines or set netcdf data record
    USE parameters, ONLY: netcdf_in
    USE io_text,    ONLY: skip_text_in

    IMPLICIT NONE

    if (netcdf_in == 0) then
       call skip_text_in()
    else
       continue
    end if

  END SUBROUTINE skip_input


  ! ------------------------------------------------------------------
  SUBROUTINE write_daily()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_text_daily

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_text_daily()
    else
       continue
    end if

  END SUBROUTINE write_daily


  ! ------------------------------------------------------------------
  SUBROUTINE write_profiles()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_text_prof

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_text_prof()
    else
       continue
    end if

  END SUBROUTINE write_profiles


  ! ------------------------------------------------------------------
  SUBROUTINE write_output()
    ! Calls the routines to open input, output and dispersion file
    USE parameters, ONLY: netcdf_out
    USE io_text,    ONLY: write_text_out

    IMPLICIT NONE

    if (netcdf_out == 0) then
       call write_text_out()
    else
       continue
    end if

  END SUBROUTINE write_output


END MODULE io
