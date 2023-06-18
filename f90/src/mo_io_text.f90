MODULE io_text

  ! This module contains routines for ascii input/output of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4, i8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lastin

  !PRIVATE :: close_text_disp
  PUBLIC :: create_text_dir !Yuan 2018.01.22
  PUBLIC :: close_in_text
  PUBLIC :: close_out_text
  !PRIVATE :: error_opening
  !PRIVATE :: error_reading
  !PRIVATE :: error_writing
  !PRIVATE :: open_disp_text
  PUBLIC :: open_in_text
  PUBLIC :: open_out_text
  PUBLIC :: read_disp_text
  PUBLIC :: read_in_text
  PUBLIC :: skip_in_text
  PUBLIC :: write_daily_text
  PUBLIC :: write_out_text
  PUBLIC :: write_prof_text
  PUBLIC :: copy_text_code!Yuan 2018.05.07


  INTEGER(i4) :: lastin

  ! ------------------------------------------------------------------

CONTAINS

  ! create reset output folder as specified in parameter file. Yuan 2018.01.22
  SUBROUTINE create_text_dir()
    USE parameters, ONLY : outdir
    IMPLICIT NONE
    CHARACTER(LEN=256) :: new_commend, new_dir, cwd, full_dir
    CHARACTER(LEN=256) :: backslash = "\"
    LOGICAL :: ex

    call getcwd(cwd)

    new_dir = trim(cwd)//trim(backslash)//trim(outdir)
    full_dir = new_dir
    ! test if a directory exists
    INQUIRE(FILE=full_dir,EXIST=ex)
    if (.not.ex) then

        new_commend = "mkdir "//trim(full_dir)
!        print *, new_commend
        call system(new_commend)

    end if

  END SUBROUTINE create_text_dir

  ! ------------------------------------------------------------------
  SUBROUTINE copy_text_code()

    USE parameters, ONLY : outdir
    USE constants, ONLY : namelist_file

! source code folder:
  !CHARACTER(LEN=256) :: folder_from = "D:\Fortran\canveg\f90\src-Copy2"
  CHARACTER(LEN=256) :: folder_from = "D:\CANVEG_compare\canveg_git\f90\src"
  CHARACTER(LEN=256) :: folder_to
  CHARACTER(LEN=256) :: cmd_tmp
  LOGICAL :: ex
  ! create a code folder inside outdir
  write(folder_to,'(a,a,a)') trim(outdir),"\","src"
!  print *, outdir
!  print *, folder_to
  INQUIRE(FILE=folder_to,EXIST=ex)
  if (.not.ex) then
    cmd_tmp = "mkdir "//trim(folder_to)
!    print *, cmd_tmp
    call system(cmd_tmp)
  end if
  write(cmd_tmp,'(a,a,a,a)') "copy ", trim(folder_from), " ", trim(folder_to)
!  print *, cmd_tmp
  call system(cmd_tmp)

  ! copy parameter file
  write(cmd_tmp,'(a,a,a,a,a,a)') "copy ", trim(namelist_file), " ", trim(folder_to),"\",trim(namelist_file)
 ! call getcwd(cmd_tmp)
 ! print *, cmd_tmp
  call system(cmd_tmp)


  END SUBROUTINE copy_text_code
  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
  SUBROUTINE close_disp_text()
    ! Close dispersion file
    USE constants,  ONLY: nindisp

    IMPLICIT NONE

    close(nindisp)

  END SUBROUTINE close_disp_text


  ! ------------------------------------------------------------------
  SUBROUTINE close_in_text()
    ! Close ascii input files
    USE constants,  ONLY: ninmet, ninlai, ninwiso
    USE parameters, ONLY: extra_nate
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    close(unit=ninmet)
    ! water isotopes
    if (iswitch%wiso == 1) close(unit=ninwiso)
    ! LAI
    if (extra_nate == 1) close(unit=ninlai)

  END SUBROUTINE close_in_text


  ! ------------------------------------------------------------------
  SUBROUTINE close_out_text()
    ! Calls the routines to open input and output files
    USE constants, ONLY: noutseas, noutopti, noutsoil, noutdaily, noutprof, &
         noutflux, nouth2osoil, &
         noutcisodaily, noutcisoseason, noutcisoprof, &
         noutwisoleaf, noutwisoprof, noutwisoflux, noutwisosoil, noutdebug, nouttpu, &
         noutnit
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    close(unit=noutseas)
    close(unit=noutopti)
    close(unit=noutsoil)
    close(unit=noutdaily)
    close(unit=noutprof)
    close(unit=noutflux)
    close(unit=nouth2osoil)
    close(unit=noutdebug)
    close(unit=nouttpu)
    close(unit=noutnit)
    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       close(unit=noutcisodaily)
       close(unit=noutcisoseason)
       close(unit=noutcisoprof)
    end if ! 13C
    ! Water isotope files
    if (iswitch%wiso == 1) then
       close(unit=noutwisoleaf)
       close(unit=noutwisoprof)
       close(unit=noutwisoflux)
       close(unit=noutwisosoil)
    end if

  END SUBROUTINE close_out_text


  ! ------------------------------------------------------------------
  SUBROUTINE error_opening(routine, stmp)
    ! writes out error message and finishes run
    USE messages, ONLY: message
    USE finishes, ONLY: finish

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: routine
    CHARACTER(len=*), INTENT(IN) :: stmp

    call message(routine, ': ', 'Error opening ', trim(stmp))
    call finish(routine, 'IO error')

  END SUBROUTINE error_opening


  ! ------------------------------------------------------------------
  SUBROUTINE error_reading(routine, stmp, text)
    ! writes out error message and finishes run
    USE string_utils, ONLY: num2str
    USE messages,     ONLY: message
    USE finishes,     ONLY: finish

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: routine
    INTEGER, INTENT(IN) :: stmp
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: text

    if (present(text)) then
       call message(routine, ': ', 'Error reading ', num2str(stmp), ' - ', text)
    else
       call message(routine, ': ', 'Error reading ', num2str(stmp))
    end if
    call finish(routine, 'IO error')

  END SUBROUTINE error_reading


  ! ------------------------------------------------------------------
  SUBROUTINE error_writing(routine, stmp, text)
    ! writes out error message and finishes run
    USE string_utils, ONLY: num2str
    USE messages,     ONLY: message
    USE finishes,     ONLY: finish

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: routine
    INTEGER, INTENT(IN) :: stmp
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: text

    if (present(text)) then
       call message(routine, ': ', 'Error writing ', num2str(stmp), ' - ', text)
    else
       call message(routine, ': ', 'Error writing ', num2str(stmp))
    end if
    call finish(routine, 'IO error')

  END SUBROUTINE error_writing


  ! ------------------------------------------------------------------
  SUBROUTINE open_disp_text()
    ! Opens the ascii input files
    USE constants,  ONLY: nindisp
    USE parameters, ONLY: indir, dispfile

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_DISP_TEXT'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr

    ! Met
    write(stmp,'(a,a,a)') trim(indir), '/', trim(dispfile)
    open(unit=nindisp, file=stmp,action="read", status="old", &
         form="formatted", iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)

  END SUBROUTINE open_disp_text


  ! ------------------------------------------------------------------
  SUBROUTINE open_in_text()
    ! Opens the ascii input files
    USE constants,  ONLY: ninmet, ninlai, ninwiso
    USE parameters, ONLY: indir, metinfile, laiinfile, wisoinfile, extra_nate, &
        scenariodir, scenariofile
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_IN_TEXT'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr

    ! Met
    select case (iswitch%scenario)
    case (0)
    write(stmp,'(a,a,a)') trim(indir), '/', trim(metinfile)
    case (1)
    write(stmp,'(a,a,a)') trim(scenariodir), '/', trim(scenariofile)
    end select

    open(unit=ninmet, file=stmp,action="read", status="old", &
         form="formatted", iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    read(ninmet,*,iostat=ierr) stmp ! read header
    if (ierr > 0) call error_reading(isroutine, ninmet)
    if (ierr < 0) call error_reading(isroutine, ninmet, 'reached EOF.')

    ! water isotopes
    if (iswitch%wiso == 1) then
       write(stmp,'(a,a,a)') trim(indir), '/', trim(wisoinfile)
       open(unit=ninwiso, file=stmp,action="read", status="old", &
            form="formatted", iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       read(ninwiso,*,iostat=ierr) stmp ! read header
       if (ierr > 0) call error_reading(isroutine, ninwiso)
       if (ierr < 0) call error_reading(isroutine, ninwiso, 'reached EOF.')
    end if

    ! LAI
    if (extra_nate == 1) then
       write(stmp,'(a,a,a)') trim(indir), '/', trim(laiinfile)
       open(unit=ninlai, file=stmp,action="read", status="old", &
            form="formatted", iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       read(ninlai,*,iostat=ierr) stmp ! read header
       if (ierr > 0) call error_reading(isroutine, ninlai)
       if (ierr < 0) call error_reading(isroutine, ninlai, 'reached EOF.')
    end if

  END SUBROUTINE open_in_text


  ! ------------------------------------------------------------------
  SUBROUTINE open_out_text()
    ! Calls the routines to open input and output files
    USE constants, ONLY: &
         noutseas, noutopti, noutsoil, noutdaily, noutprof, &                ! units
         noutflux, nouth2osoil, nouttpu, noutnit, &
         noutcisodaily, noutcisoseason, noutcisoprof, &
         noutwisoleaf, noutwisoprof, noutwisoflux, noutwisosoil, noutdebug, &
         dailyfile, daily13cfile, hourlyfile, hourly13cfile, optimisefile, & ! filenames
         profilefile, profile13cfile, profileisofile, fluxprofilefile, &
         fluxprofileisofile, soilfile, h2osoilfile, &
         h2osoilisofile, h2oleafisofile, debugfile, tpufile, nitfile
    USE parameters, ONLY: outdir, outsuffix
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_OUT_TEXT'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr
    CHARACTER(LEN=30) :: form1

    ierr = 0
    ! Hourly/Season
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(hourlyfile), trim(outsuffix)
    open(unit=noutseas, file=stmp,action="write", status="replace", &
         form="formatted", recl=31*25, iostat=ierr) ! add hourlyROC Yuan 2018.05.07
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 31, '(",",a))'
!    print *, form1
    write(noutseas,form1) "daytime ", "netrn ", "sumrn ", "sumh ", "sumle ", &
    "canps ", "gpp ", "canresp ", "soilresp ", "boleresp ", &
    "sumOXY ", "netOXY ", "netROC ", "ROC_leaf ", "ROC_bole ", "ROC_soil ", &
    "ustar ", "Kdiff ","phim ","canresp_o ", &
    "inputPAR ","PAR_dir ","PAR_dn ", "PAR_up ", &
    "NIR_dir ","NIR_dn ","NIR_up ","IR_dn ","IR_up ", &
    "Nsupply ", "Ndemand"
    !write(noutseas,form1) &
    !     "daytime", "netrn", "sumrn", "sumh", "sumle"

    ! Optimise
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(optimisefile), trim(outsuffix)
    open(unit=noutopti, file=stmp,action="write", status="replace", &
         form="formatted", recl=2*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 2-1, '(",",a))'
    write(noutopti,form1) "daytime ", "soil_mm_50"

    ! Soil
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(soilfile), trim(outsuffix)
    open(unit=noutsoil, file=stmp,action="write", status="replace", &
         form="formatted", recl=26*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 26-1, '(",",a))'
    write(noutsoil,form1) "daytime", "netrn", "soilh", "soille", "soil%tsrf", &
         "soil%T%base", "soil%T%15cm", "flux%transp", "flux%evap", "prof%throughfall", &
         "soil%soil_mm", "soil%qinfl", "soil%qdrai", "soil%gsoil", "soil%qtran", &
         "soil%surfrun", &
         "Tsoil%1", "Tsoil%2", "Tsoil%3", "Tsoil%4", "Tsoil%5", &
         "Tsoil%6", "Tsoil%7", "Tsoil%8", &
         "Tsoil%9", "Tsoil%10"

    ! Daily average
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(dailyfile), trim(outsuffix)
    open(unit=noutdaily, file=stmp,action="write", status="replace", &
         form="formatted", recl=22*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 22-1, '(",",a))'
    write(noutdaily,form1) "Day ", "Avg_NEE ", "Avg_EVAP ", &
         "Avg_H ", "Avg_PAR ", "Avg_RNET ", "lai ", "wai ", "pai ", "Avg_PS ", &
         "Ave_Resp ", "Avg_BOLE ", "Avg_SOIL ", "Avg_TLeaf ", "Avg_Gs ", "Tleaf_day ", &
         "Avg_GPP ", "Avg_OXY ", "Net_OXY ", "Avg_ROC ", "Avg_respo" ! add daily mean ROC
    !     form="formatted", recl=14*25, iostat=ierr)
   ! if (ierr > 0) call error_opening(isroutine, stmp)
   ! write(form1,'(A,I3,A)') '(a,', 14-1, '(",",a))'
   ! write(noutdaily,form1) "Day", "Avg_FC", "Avg_EVAP", "AVG_H", &
   !      "Avg_PAR", "Avg_RNET", "lai", "Avg_PS", &
   !      "Ave_Resp", "Avg_BOLE", "Avg_SOIL", "Avg_TLeaf", "Avg_Gs", "Tleaf_day"

    ! Profile air
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profilefile), trim(outsuffix)
    open(unit=noutprof, file=stmp,action="write", status="replace", &
         form="formatted", recl=13*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 14-1, '(",",a))'
   ! print *, form1
    write(noutprof,form1) "Time ", "i ", "tair ", "tair_f ", "qair ", "co2 ", "o2 ", "co2_f ", "o2_f ", &
          "co2_soil ", "o2_soil ", "wnd ", "vpd"!"co2_disp ", "o2_disp" !
    !write(form1,'(A,I3,A)') '(a,', 5-1, '(",",a))'
   ! write(noutprof,form1) "Time", "i", "tair", "qair", "co2"

    ! Profile fluxes
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(fluxprofilefile), trim(outsuffix)
    open(unit=noutflux, file=stmp,action="write", status="replace", &
         form="formatted", recl=40*86, iostat=ierr) ! original 40*25
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 86, '(",",a))'
    write(noutflux,form1) "daytime ", "i ", "dHdz ", "dLEdz ", "dLEdz%sun ", "dLEdz%shd ", &
         "sun_A ", "shd_A ", "dGPPdz ", "dGPPdz%sun ", "dGPPdz%shd ", "dPsdz ", &
         "dPsdz%sun ", "dPsdz%shd ", "dRESPdz ", "dRESPdz%sun ", "dRESPdz%shd ", &
         "PARdirect ", "PARdiffuse ", "Tleaf ", "Tleaf_sun ", "Tleaf_shd ", &
         "Beam ", "Nonbeam ", "LAI ", "gs ", "gs%sun ", "gs%shd ", &
         "PARin ", "CO2_source ", "O2_source ", "RQ_sun ", "RQ_shade ", "RESP_O2 ", &
         "bole_resp ", "ROC_layer ", "PAI ", "WAI ", "RN ", &
         "netrn_sun ", "netrn_shd ","par_sun ", "nir_sun ", "par_shd ", "nir_shd ", &
         "par_up ","par_dn ","nir_up ","nir_dn ","ir_dn ", "ir_up ", &
         "par_beam_dn ", "nir_beam_dn ", &
         "T_sun_f ", "T_shd_f ", "filt_T ", "filt_H ", "filt_LE ", &
         "sun_rs ", "shd_rs ", "LAI_sun ", "LAI_shd ", "par_diff ", "nir_diff ", &
         "PSN_O2 ", "cws ", "wet coef ", "sun_tpucoef ", "shd_tpucoef ","GPO ", "GPO_sun ","GPO_shd ", &
         "Ja_sun ", "Jglu_sun ", "JBusch_sun ", "Ja_shd ", "Jglu_shd ", "JBusch_shd ", &
         "iphoton_sun ", "iphoton_shd ", "sun_quad ", "shd_quad ", "ci_sun ", "ci_shd "
        ! write(form1,'(A,I3,A)') '(a,', 44-1, '(",",a))'
    !write(noutflux,form1) "daytime",

    ! Profile N ass
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(nitfile), trim(outsuffix)
    open(unit=noutnit, file=stmp,action="write", status="replace", &
         form="formatted", recl=40*11, iostat=ierr) ! original 40*25
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 11, '(",",a))'
    write(noutnit,form1) "daytime ", "i ", &
                       "dNdemanddz ", "dNdemand_sundz ", "dNdemand_shddz ", &
                       "dNtotdz ", "dNtot_sundz ", "dNtot_shddz ",&
                       "dNBuschdz ", "dNBusch_sundz ", "dNBusch_shddz"


    ! Profile tpu limitation
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(tpufile), trim(outsuffix)
    open(unit=nouttpu, file=stmp,action="write", status="replace", &
         form="formatted", recl=40*41, iostat=ierr) ! original 40*25
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 43, '(",",a))'
    write(nouttpu,form1) "daytime ", "i ", &
                       "lai ", "air_co2 ", &
                       "sun_lai_f ", "shd_lai_f ", &
                       "sun_ci ", "shd_ci ", &
                       "sun_cs ", "shd_cs ", &
                       "sun_cc ", "shd_cc ", &
                       "sun_wc ", "shd_wc ", &
                       "sun_wj ", "shd_wj ", &
                       "sun_wp ", "shd_wp ", &
                       "sun_alphag ", "shd_alphag ", &
                       "sun_alphas ", "shd_alphas ", &
                       "sun_Ja ", "shd_Ja ", &
                       "sun_Jglu ", "shd_Jglu ", &
                       "sun_JBusch ", "shd_JBusch ", &
                       "sun_Ndemand ", "shd_Ndemand ", &
                       "sun_Ntot ", "shd_Ntot ", &
                       "sun_NBusch ", "shd_NBusch ", &
                       "sun_NO3 "   , "shd_NO3 "   , &
                       "sun_NO2 "   , "shd_NO2 "   , &
                       "sun_NH4 "   , "shd_NH4 "   , &
                       "sun_tpu_coeff ", "shd_tpu_coeff"

    ! Soil water
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2osoilfile), trim(outsuffix)
    open(unit=nouth2osoil, file=stmp,action="write", status="replace", &
         form="formatted", recl=26*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(form1,'(A,I3,A)') '(a,', 26-1, '(",",a))'
    write(nouth2osoil,form1) "Time ", "wiso%lost1 ", "flux%litterevap1 ", &
         "flux%soilevap1 ", "flux%evap1 ", "flux%c_evaporation1 ", &
         "flux%c_transpiration1 ", "flux%c_evapotranspiration1 ", &
         "flux%evapotranspiration1 ", "prof%rhov_air11 ", "prof%rhov_air1201 ", &
         "input%ppt1 ", "prof%throughfall1 ", "flux%surfrun1 ", "soil%qdrai1 ", &
         "soil%thetal1 ", "soil%theta01 ", "soil%theta11 ", "soil%theta21 ", &
         "soil%theta31 ", "soil%theta41 ", "soil%theta51 ", "soil%theta61 ", &
         "soil%theta71 ", "soil%theta81 ", "soil%theta91"
    ! debugfile:
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(debugfile), trim(".txt")
    open(unit=noutdebug, file=stmp,action="write", status="replace", &
         form="formatted", iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutdebug,*) "This is a run log or debug!"

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Daily average
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(daily13cfile), trim(outsuffix)
       open(unit=noutcisodaily, file=stmp,action="write", status="replace", &
            form="formatted", recl=3*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 3-1, '(",",a))'
       write(noutcisodaily,form1) "Day", "Avg_Fc_13C", "Avg_disc13_day"
       ! Hourly/Season
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(hourly13cfile), trim(outsuffix)
       open(unit=noutcisoseason, file=stmp,action="write", status="replace", &
            form="formatted", recl=7*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 6-1, '(",",a))'
       write(noutcisoseason,form1) "daytime", "disc13", "disc13_long", "ave_dair_13C", &
            "disc13Clong40", "disc13C40"
       ! Profile air
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profile13cfile), trim(outsuffix)
       open(unit=noutcisoprof, file=stmp,action="write", status="replace", &
            form="formatted", recl=5*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 5-1, '(",",a))'
       write(noutcisoprof,form1) "Time", "i", "R13C", "del13C", "13C"
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso leaf
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2oleafisofile), trim(outsuffix)
       open(unit=noutwisoleaf, file=stmp,action="write", status="replace", &
            form="formatted", recl=20*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 7-1, '(",",a))'
       write(noutwisoleaf,form1) "Time", "soil%rxylem2", "soil%rxylem3", "soil%rxylem4", &
            "prof%rvapour12", "prof%rvapour13", "prof%rvapour14"
       ! Iso profile air
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profileisofile), trim(outsuffix)
       open(unit=noutwisoprof, file=stmp,action="write", status="replace", &
            form="formatted", recl=9*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 5-1, '(",",a))'
       write(noutwisoprof,form1) "Time", "i", "qair2", "qair3", "qair4"
       ! Iso profile flux
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(fluxprofileisofile), trim(outsuffix)
       open(unit=noutwisoflux, file=stmp,action="write", status="replace", &
            form="formatted", recl=29*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 29-1, '(",",a))'
       write(noutwisoflux,form1) "daytime", "i", "dLEdz2", "dLEdz3", "dLEdz4", &
            "sun_craig2", "sun_craig3", "sun_craig4", "shd_craig2", "shd_craig3", &
            "shd_craig4", "sun_leafwater_e2", "sun_leafwater_e3", "sun_leafwater_e4", &
            "shd_leafwater_e2", "shd_leafwater_e3", "shd_leafwater_e4", &
            "sun_leafwater2", "sun_leafwater3", "sun_leafwater4", "shd_leafwater2", &
            "shd_leafwater3", "shd_leafwater4", "sun_trans2", "sun_trans3", &
            "sun_trans4", "shd_trans2", "shd_trans3", "shd_trans4"
       ! Iso soil
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2osoilisofile), trim(outsuffix)
       open(unit=noutwisosoil, file=stmp,action="write", status="replace", &
            form="formatted", recl=101*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(form1,'(A,I3,A)') '(a,', 101-1, '(",",a))'
       write(noutwisosoil,form1) "Time", "wiso%lost1", "wiso%lost2", "wiso%lost3", &
            "wiso%lost4", "flux%litterevap1", "flux%litterevap2", "flux%litterevap3", &
            "flux%litterevap4", "flux%soilevap1", "flux%soilevap2", "flux%soilevap3", &
            "flux%soilevap4", "flux%s_evap1", "flux%s_evap2", "flux%s_evap3", "flux%s_evap4", &
            "flux%c_evaporation1", "flux%c_evaporation2", "flux%c_evaporation3", &
            "flux%c_evaporation4", "flux%c_transpiration1", "flux%c_transpiration2", &
            "flux%c_transpiration3", "flux%c_transpiration4", "flux%c_evapotranspiration1", &
            "flux%c_evapotranspiration2", "flux%c_evapotranspiration3", "flux%c_evapotranspiration4", &
            "flux%evapotranspiration1", "flux%evapotranspiration2", "flux%evapotranspiration3", &
            "flux%evapotranspiration4", "prof%rhov_air11", "prof%rhov_air12", "prof%rhov_air13", &
            "prof%rhov_air14", "prof%rhov_air1201", "prof%rhov_air1202", "prof%rhov_air1203", &
            "prof%rhov_air1204", "input%ppt1", "input%ppt2", "input%ppt3", "input%ppt4", &
            "prof%throughfall1", "prof%throughfall2", "prof%throughfall3", "prof%throughfall4", &
            "flux%surfrun1", "flux%surfrun2", "flux%surfrun3", "flux%surfrun4", "soil%qdrai1", &
            "soil%qdrai2", "soil%qdrai3", "soil%qdrai4", "soil%thetal1", "soil%thetal2", &
            "soil%thetal3", "soil%thetal4", "soil%theta01", "soil%theta02", "soil%theta03", &
            "soil%theta04", "soil%theta11", "soil%theta12", "soil%theta13", "soil%theta14", &
            "soil%theta21", "soil%theta22", "soil%theta23", "soil%theta24", "soil%theta31", &
            "soil%theta32", "soil%theta33", "soil%theta34", "soil%theta41", "soil%theta42", &
            "soil%theta43", "soil%theta44", "soil%theta51", "soil%theta52", "soil%theta53", &
            "soil%theta54", "soil%theta61", "soil%theta62", "soil%theta63", "soil%theta64", &
            "soil%theta71", "soil%theta72", "soil%theta73", "soil%theta74", "soil%theta81", &
            "soil%theta82", "soil%theta83", "soil%theta84", "soil%theta91", "soil%theta92", &
            "soil%theta93", "soil%theta94"
    end if

  END SUBROUTINE open_out_text


  ! ------------------------------------------------------------------
  SUBROUTINE read_disp_text()
    ! Input data on Thomson dispersion matrix that was computed offline with
    ! MOVOAK.C, Dij (s m-1)
    USE constants, ONLY: nindisp
    USE types,     ONLY: met
    USE setup,     ONLY: ncl, ntl

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'READ_DISP_TEXT'
    INTEGER :: ierr
    INTEGER(i4) :: i, j
    REAL(wp) :: in1, in2

    call open_disp_text()

    ierr = 0
    do j=1, ncl
       do i=1, ntl
          read(nindisp,*,iostat=ierr) in1, in2
          met%dispersion(i,j) = in1
       end do
    end do
    if (ierr > 0) call error_reading(isroutine, nindisp)
    if (ierr < 0) call error_reading(isroutine, nindisp, 'reached EOF.')

    call close_disp_text()

  END SUBROUTINE read_disp_text


  ! ------------------------------------------------------------------
  SUBROUTINE read_in_text()
    ! input data and check for bad data
    ! note that the data were produced in single precision (float)
    ! so I had to read them as single precision, otherwise I ingested garbage
    USE constants,     ONLY: ninmet, ninwiso, ninlai, zero, one, two, e1, e3, &
         TN0, rugc, vonKarman, mass_air, o2_ref
    USE types,         ONLY: met, input, solar, time, srf_res, wiso, iswitch
    USE setup,         ONLY: ncl, nwiso
    USE parameters,    ONLY: extra_nate, bprime, zm, zd, z0, end_run, lai, scenario_c, scenario_temp
    USE utils,         ONLY: es
    USE isotope_utils, ONLY: invdelta1000_h2o

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'READ_IN_TEXT'
    INTEGER :: ierr
    INTEGER(i4) :: i, dayy, flag
    REAL(wp) :: hhrr, in01, in02, in03, in04, in05, in06, in07, in08
    REAL(wp) :: in09, in10, in11, in13, in14
    REAL(wp) :: in15 ! in15 for o2 con ppm Yuan 2018.02.14
    REAL(wp) :: in16 ! in16 for chamber ER_An measurements Yuan 2012.07.14
    INTEGER(i8) :: dt
    INTEGER(i4), DIMENSION(nwiso-1) :: mc

    dt = int(time%time_step,kind=i8) * 100_i8 / 3600_i8
    ! use input%dayy instead of time%days because of several years at once
    time%jdold = input%dayy ! identify previous day
    read(ninmet,*,iostat=ierr) dayy, hhrr, in01, in02, in03, in04, &
         in05, in06, in07, in08, in09, in10, in11, flag, in13, in14, in15, in16

    ! Although O2 is included as "in15" in input file, it is more convenient to calculate atm O2 directly in the model,
    ! especially when RCP CO2 concentration is implemented.

    in15=o2_ref-1.15_wp*in08
    !the following 3 lines are not needed when RCP scenario files are used,
    !cause scenario_temp and scenario_c are only a fixed value of annual T and CO2 variation.
    in01=in01+scenario_temp
    in08=in08+scenario_c ! input CO2 ppm
    in15=in15-1.15_wp*scenario_c ! deltaO2+refO2=real O2 ppm Yuan 2018.02.14
    if (ierr > 0) call error_reading(isroutine, ninmet)
    if (ierr < 0) call error_reading(isroutine, ninmet, 'reached EOF.')
    input%dayy         = dayy
    input%hhrr         = hhrr
    input%ta           = in01
    input%rglobal      = in02
    input%parin        = in03
    input%pardif       = in04
    input%ea           = in05
    input%wnd          = in06
    input%ppt(1)       = in07
    input%co2air       = in08
    input%press_mb     = in09
    input%tsoil        = in10
    input%soilmoisture = in11
    input%flag         = flag
    input%d13CO2       = in13
    input%d18CO2       = in14
    input%o2air        = in15
    input%ER           = in16
   !print *, "check scenario file:    ", input%co2air
    ! write(*,'(a,i10,3f20.14)') 'RI01.01 ', input%dayy, input%hhrr, input%ta
    ! write(*,'(a,3f20.14)') 'RI01.02 ', input%rglobal, input%parin, input%pardif
    ! write(*,'(a,3f20.14)') 'RI01.03 ', input%ea, input%wnd, input%ppt(1)
    ! write(*,'(a,3f20.14)') 'RI01.04 ', input%co2air, input%press_mb, input%tsoil
    ! write(*,'(a,f20.14,i10,3f20.14)') 'RI01.05 ', input%soilmoisture, input%flag, input%d13CO2
    ! write(*,'(a,3f20.14)') 'RI01.06 ', input%d18CO2

    ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
    if (extra_nate == 1) input%lai = input%pardif ! redundant with fptr15
    ! preliminary approach to calculate consecutive year
    time%year = time%year0 + (input%dayy-1)/365
    time%daytime = int(input%dayy,kind=i8)*10000_i8 + int(input%hhrr*100._wp,kind=i8) ! define daytime
    time%doy = mod(input%dayy,365)
    if (time%doy == 0) time%doy = 365
    time%local_time = input%hhrr
    !print *, input%dayy
    time%days = time%doy
    !print *, time%days
    ! compute derived quantities for the model
    met%T_Kelvin = input%ta + TN0 ! compute absolute air temperature
!   print *, input%ea
    select case (iswitch%scenario) ! in RCP file, column ea is actually relative humidity
    case (0)
        met%relative_humidity = input%ea*10._wp/es(met%T_Kelvin) ! relative humidity

    case (1)
        met%relative_humidity = input%ea
        input%ea=met%relative_humidity*es(met%T_Kelvin)/10._wp
    end select
!print *, iswitch%scenario,input%ea,met%relative_humidity
    met%rhova_g = input%ea * 2165._wp/met%T_Kelvin ! compute absolute humidity, g m-3
    ! limit humidity
    if (met%rhova_g < zero) met%rhova_g = zero
    if (met%rhova_g > 30._wp) met%rhova_g = 30._wp
    met%rhova_kg = met%rhova_g *e3 ! absolute humidity, kg m-3
    met%press_kpa = input%press_mb *e1 ! air pressure, kPa
    met%press_bars = input%press_mb * e3 ! air pressure, bars
    met%press_Pa = met%press_kpa*1000._wp ! pressure, Pa

    ! combining gas law constants
    met%pstat273 = rugc / (100000._wp * met%press_bars)
    ! cuticular conductance adjusted for pressure and T, mol m-2 s-1
    ! cuticular resistance
    srf_res%rcuticle(1:ncl) = one / (bprime(1:ncl) * met%T_Kelvin * met%pstat273)
    ! check for bad CO2 data
    !if (abs(input%co2air) >= 998._wp) input%co2air = 370._wp
    ! write(*,'(a,3f20.14)') 'RI01.07 ', input%parin
    if (input%parin < zero) input%parin = zero ! check for bad par data
 !   print *, "inputPAR:    ",input%parin
    ! write(*,'(a,3f20.14)') 'RI01.08 ', input%parin
    ! if (input%parin <= zero) then ! check for night
    !    solar%par_beam = zero
    !    solar%par_diffuse = zero
    !    solar%nir_beam = zero
    !    solar%nir_diffuse = zero
    ! end if
    ! write(*,'(a,3f20.14)') 'RI01.09 ', solar%par_diffuse, solar%nir_beam, solar%nir_diffuse
    ! write(*,'(a,3f20.14)') 'RI01.10 ', solar%ratrad
    ! set some limits on bad input data to minimize the model from blowing up
    if (solar%ratrad > 0.9_wp .or. solar%ratrad < 0.2_wp) solar%ratrad = 0.5_wp
    ! write(*,'(a,3f20.14)') 'RI01.11 ', solar%ratrad
    ! limit wind speed
    if (input%wnd < one) input%wnd = one
    if (input%wnd > 10._wp) input%wnd = 10._wp
    met%ustar_filter = vonKarman/log((zm-zd)/z0)*input%wnd
    ! limit air temperature
    if (input%ta < -30._wp) input%ta = -30._wp
    if (input%ta > 60._wp)  input%ta = 60._wp
    ! air density, kg m-3
    met%air_density = met%press_kpa * mass_air / (rugc * met%T_Kelvin)
    ! air density, mole m-3
    met%air_density_mole = met%press_kpa/ (rugc * met%T_Kelvin) * 1000._wp

    ! water isotopes
    if (iswitch%wiso == 1) then
       read(ninwiso,*,iostat=ierr) dayy, hhrr, in01, in02, in03
       if (ierr > 0) call error_reading(isroutine, ninwiso)
       if (ierr < 0) call error_reading(isroutine, ninwiso, 'reached EOF.')
       input%dppt(1)    = zero
       input%dppt(2)    = real(in01,kind=wp)
       input%dppt(3)    = real(in02,kind=wp)
       input%dppt(4)    = zero
       input%dvapour(2) = real(in03,kind=wp)
       ! Test Meteoric water line for deuterium
       input%dppt(3) = 8._wp * input%dppt(2) + 10._wp
       if (wiso%nofracin == 1) input%dppt(1:nwiso) = wiso%dtheta(1,2:nwiso)
       forall(i=1:nwiso-1) mc(i) = i+1
       input%ppt(2:nwiso) = input%ppt(1) * invdelta1000_h2o(input%dppt(2:nwiso), mc(1:nwiso-1))
    end if

    ! LAI
    if (extra_nate == 1) then
       read(ninlai,*,iostat=ierr) dayy, hhrr, in01, in02
       if (ierr > 0) call error_reading(isroutine, ninlai)
       if (ierr < 0) call error_reading(isroutine, ninlai, 'reached EOF.')
       input%lai_up   = real(in01,kind=wp)
       input%lai_down = real(in02,kind=wp)
       ! Test: set LAI to fix value of 1 + 1 m2 m-2
       input%lai_up   = one
       input%lai_down = one
       input%lai      = input%lai_up + input%lai_down
       lai            = two
    end if

    ! Check for end of run or end-of-file
    if (ierr < 0) lastin = 1
    if ((time%daytime+dt) > end_run) lastin = 1

  END SUBROUTINE read_in_text


  ! ------------------------------------------------------------------
  SUBROUTINE skip_in_text()
    ! Skip input lines
    USE constants,  ONLY: ninmet, ninlai, ninwiso
    USE parameters, ONLY: extra_nate, start_run
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'SKIP_IN_TEXT'
    INTEGER :: ierr
    LOGICAL :: notreached
    REAL :: dayy, hhrr, in01, in02, in03, in04, in05, in06, in07, in08
    REAL :: in09, in10, in11, in12, in13, in14, in15
    INTEGER(i8) :: daytime

    lastin = 0

    ierr = 0
    notreached = .true.
    do while ((ierr==0) .and. notreached)
       read(ninmet,*,iostat=ierr) dayy, hhrr, in01, in02, in03, in04, &
            in05, in06, in07, in08, in09, in10, in11, in12, in13, in14, in15
       daytime = int(dayy,kind=i8)*10000_i8 + int(hhrr*100., kind=i8)
       if (daytime >= start_run) notreached = .false.
    end do
    if (ierr > 0) call error_reading(isroutine, ninmet)
    if (ierr < 0) call error_reading(isroutine, ninmet, 'reached EOF.')
    backspace ninmet

    ! water isotopes
    if (iswitch%wiso == 1) then
       notreached = .true.
       do while ((ierr==0) .and. notreached)
          read(ninwiso,*,iostat=ierr) dayy, hhrr, in01, in02, in03
          daytime = int(dayy,kind=i8)*10000_i8 + int(hhrr*100., kind=i8)
          if (daytime >= start_run) notreached = .false.
       end do
       if (ierr > 0) call error_reading(isroutine, ninwiso)
       if (ierr < 0) call error_reading(isroutine, ninwiso, 'reached EOF.')
       backspace ninwiso
    end if

    ! LAI
    if (extra_nate == 1) then
       notreached = .true.
       do while ((ierr==0) .and. notreached)
          read(ninlai,*,iostat=ierr) dayy, hhrr, in01, in02
          daytime = int(dayy,kind=i8)*10000_i8 + int(hhrr*100., kind=i8)
          if (daytime >= start_run) notreached = .false.
       end do
       if (ierr > 0) call error_reading(isroutine, ninlai)
       if (ierr < 0) call error_reading(isroutine, ninlai, 'reached EOF.')
       backspace ninlai
    end if

  END SUBROUTINE skip_in_text


  ! ------------------------------------------------------------------
  SUBROUTINE write_daily_text()
    ! Write out daily means
    USE constants, ONLY: noutdaily, noutcisodaily
    USE types,     ONLY: iswitch, time, output, input

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_DAILY_TEXT'
    INTEGER :: ierr
    CHARACTER(LEN=30) :: form1

    ierr = 0
    ! Daily average!output%sumo, output%sumneto,
    write(form1,'(A,I3,A)') '(i03,', 21-1, '(",",es22.14))'
    write(noutdaily,form1,iostat=ierr) &
         input%dayy, output%sumfc, &
         output%sumevap, output%sumsens, &
         output%sumpar, output%sumnet, time%lai, time%wai, time%pai, &
         output%sumps, output%sumresp, output%sumbole, output%sumsoil, &
         output%sumta, output%sumgs, output%sumTleaf, &
         output%sumgpp, output%sumo, output%sumneto, output%sumROC, output%sumresp_o
!         print *, output%sumgpp, output%sumfc, output%sumo, output%sumneto
    if (ierr > 0) call error_writing(isroutine, noutdaily)

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Daily average
       write(form1,'(A,I3,A)') '(i03,', 3-1, '(",",es22.14))'
       write(noutcisodaily,form1,iostat=ierr) &
            input%dayy, output%sumF13C, output%sumdisc13C
       if (ierr > 0) call error_writing(isroutine, noutcisodaily)
    end if ! 13C

  END SUBROUTINE write_daily_text


  ! ------------------------------------------------------------------
  SUBROUTINE write_out_text()
    ! Writes output except the profile
    USE constants, ONLY: noutseas, noutopti, noutsoil, nouth2osoil, &
         noutcisoseason, noutwisoleaf, noutwisosoil, noutdebug
    USE types,     ONLY: iswitch, time, prof, soil, bole, &
         input, output, flux, met, solar, wiso, nitrogen
    USE setup,     ONLY: ncl, ntl, nwiso, nsoil
    USE parameters,    ONLY: ROC_leaf_in, ROC_bole_in, ROC_soil_in

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_OUT_TEXT'
    INTEGER :: ierr
    CHARACTER(LEN=30) :: form1

    ierr = 0
    ! Hourly/Season
!write(form1,'(A,I3,A)') '(i07,",",i03,', 19-1+4, '(",",es22.14))'
    write(form1,'(A,I3,A)') '(i07,', 31, '(",",es22.14))'
    !print *, form1
    write(noutseas,form1,iostat=ierr) &
         time%daytime, output%netrad, output%sumrn, output%sumh, output%sumle, &
         output%can_ps_mol, output%can_gpp, output%canresp, soil%respiration_mole, bole%respiration_mole, &
         output%houro, output%hourneto, output%hourROC, ROC_leaf_in, ROC_bole_in, ROC_soil_in, &
         met%ustar_filter, met%K, met%phim, output%hour_canrespo, &
         input%parin,solar%beam_flux_par(ncl+1)/ 4.6_wp,solar%par_down(ncl+1)/ 4.6_wp,solar%par_up(ncl+1)/ 4.6_wp,&
         solar%beam_flux_nir(ncl+1),solar%nir_dn(ncl+1),solar%nir_up(ncl+1),&
         solar%ir_dn(ncl+1),solar%ir_up(ncl+1),nitrogen%Nsupply, nitrogen%Ndemand
    if (ierr > 0) call error_writing(isroutine, noutseas)

    ! Optimise
    write(form1,'(A,I3,A)') '(i07,', 2-1, '(",",es22.14))'
    write(noutopti,form1,iostat=ierr) time%daytime, soil%soil_mm_50
    if (ierr > 0) call error_writing(isroutine, noutopti)

    ! Soil
    write(form1,'(A,I3,A)') '(i07,', nsoil+16-1, '(",",es22.14))'
    write(noutsoil,form1,iostat=ierr) &
         time%daytime, output%rnet_soil, soil%heat, soil%evap, &
         soil%tsrf, soil%T_base, soil%T_15cm, flux%c_transpiration(1), &
         flux%s_evap(1), prof%throughfall(1,1), soil%soil_mm, &
         soil%qinfl(1), soil%qdrai(1), output%c7,soil%qtran(1), &
         flux%surfrun(1), soil%T_soil(1:nsoil)
    if (ierr > 0) call error_writing(isroutine, noutsoil)

    ! Soil water
    write(form1,'(A,I3,A)') '(i07,', nsoil+16-1, '(",",es22.14))'
    write(nouth2osoil,form1,iostat=ierr) &
         time%daytime, wiso%lost(1), flux%litterevap(1), &
         flux%soilevap(1), flux%s_evap(1), flux%c_evaporation(1), &
         flux%c_transpiration(1), flux%c_evapotranspiration(1), &
         flux%evapotranspiration(1), prof%rhov_air(1,1), &
         prof%rhov_air(ntl,1), input%ppt(1), prof%throughfall(1,1), &
         flux%surfrun(1), soil%qdrai(1), soil%theta_l(1), &
         soil%theta(1:nsoil,1)
    if (ierr > 0) call error_writing(isroutine, nouth2osoil)

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Hourly/Season
       write(form1,'(A,I3,A)') '(i07,', 6-1, '(",",es22.14))'
       write(noutcisoseason,form1,iostat=ierr) &
            time%daytime, output%ave_disc13, output%ave_disc13_long, &
            output%ave_daC13, prof%disc13C_long(ncl), prof%disc13C(ncl)
       if (ierr > 0) call error_writing(isroutine, noutcisoseason)
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso leaf
       write(form1,'(A,I3,A)') '(i07,', 7-1, '(",",es22.14))'
       write(noutwisoleaf,form1,iostat=ierr) &
            time%daytime, soil%rxylem(2), soil%rxylem(3), soil%rxylem(4), &
            prof%rvapour(35,2), prof%rvapour(35,3), prof%rvapour(35,4), &
            prof%rhov_air_filter(35,1:nwiso), wiso%lost(1:nwiso), &
            prof%dvapour(35,1:nwiso)
       if (ierr > 0) call error_writing(isroutine, noutwisoleaf)
       ! Iso soil
       write(form1,'(A,I3,A)') '(i07,', 1+15*nwiso+nsoil*nwiso-1, '(",",es22.14))'
       write(noutwisosoil,form1,iostat=ierr) &
            time%daytime, wiso%lost(1:nwiso), &
            flux%litterevap(1:nwiso), flux%soilevap(1:nwiso), &
            flux%s_evap(1:nwiso), flux%c_evaporation(1:nwiso), &
            flux%c_transpiration(1:nwiso), &
            flux%c_evapotranspiration(1:nwiso), &
            flux%evapotranspiration(1:nwiso), &
            prof%rhov_air(1,1:nwiso), prof%rhov_air(ntl,1:nwiso), &
            input%ppt(1:nwiso), prof%throughfall(1,1:nwiso), &
            flux%surfrun(1:nwiso), soil%qdrai(1:nwiso), soil%theta_l(1:nwiso), &
            transpose(soil%theta(1:nsoil,1:nwiso))
       if (ierr > 0) call error_writing(isroutine, noutwisosoil)
    end if

  END SUBROUTINE write_out_text


  ! ------------------------------------------------------------------
  SUBROUTINE write_prof_text()
    ! Calls the routines to open input and output files
    USE constants, ONLY: noutprof, noutflux, noutcisoprof, noutwisoprof, noutwisoflux, nouttpu, noutnit
    USE types,     ONLY: iswitch, prof, solar, time, input, fact
    USE setup,     ONLY: ncl, ntl, nwiso

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_PROF_TEXT'
    INTEGER :: ierr
    INTEGER(i4) :: j
    REAL(wp), DIMENSION(nwiso-1) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    CHARACTER(LEN=30) :: form1

    ierr = 0
    ! Profile air
    write(form1,'(A,I3,A)') '(i07,",",i03,', 13-1, '(",",es22.14))'
    do j=1, ntl
       write(noutprof,form1,iostat=ierr) &
            time%daytime, j, prof%tair(j), prof%tair_filter(j), prof%rhov_air(j,1), &
            prof%co2_air(j), prof%O2_air(j),  prof%co2_air_filter(j), prof%O2_air_filter(j), &! oxygen module, Yuan 2018.01.19
            prof%CO2_soil(j), prof%O2_soil(j), prof%u(j) , prof%vpd_air(j) !prof%CO2_disp(j), prof%O2_disp(j)
    end do
    if (ierr > 0) call error_writing(isroutine, noutprof, ' - 1')

    ! Profile fluxes
    write(form1,'(A,I3,A)') '(i07,",",i03,', 86, '(",",es22.14))'
    do j=1, ncl
       write(noutflux,form1,iostat=ierr) &
            time%daytime,j,prof%dHdz(j), prof%dLEdz(j,1), &
            prof%dLEdz_sun(j), prof%dLEdz_shd(j), prof%sun_A(j), &
            prof%shd_A(j), prof%dGPPdz(j), prof%dGPPdz_sun(j), &
            prof%dGPPdz_shd(j), prof%dPsdz(j), prof%dPsdz_sun(j), &
            prof%dPsdz_shd(j), prof%dRESPdz(j), prof%dRESPdz_sun(j), &
            prof%dRESPdz_shd(j), solar%quantum_sun(j), &
            solar%quantum_shd(j), prof%tleaf(j), prof%sun_tleaf(j), &
            prof%shd_tleaf(j), solar%prob_beam(j), solar%prob_shd(j), &
            prof%dLAIdz(j), &
            prof%dStomCondz(j)*1e3_wp,prof%dStomCondz_sun(j)*1e3_wp,prof%dStomCondz_shd(j)*1e3_wp, &
            input%parin, prof%source_co2(j), prof%source_O2(j), prof%RQ_sun(j), prof%RQ_shd(j), prof%dRESPdz_O2(j),&
            prof%dboledz(j),prof%ROC_layer(j),prof%dPAIdz(j),prof%dWAIdz(j), &
            prof%dRNdz(j), solar%rnet_sun(j), solar%rnet_shd(j), solar%par_sun(j), solar%nir_sun(j),  &
            solar%par_shd(j), solar%nir_shd(j), &
            solar%par_up(j), solar%par_down(j), &
            solar%nir_up(j), solar%nir_dn(j), &
            solar%ir_dn(j), solar%ir_up(j), &
            solar%beam_flux_par(j), solar%beam_flux_nir(j), &
            prof%sun_tleaf_filter(j), prof%shd_tleaf_filter(j), fact%a_filt, fact%heatcoef, fact%latent, &
            prof%sun_rs_filter(j), prof%shd_rs_filter(j), prof%sun_lai(j),prof%shd_lai(j),solar%par_diffuse,solar%nir_diffuse, &
            prof%dPsdz_O2(j), prof%cws(j,1), prof%wet_coef(j), real(prof%sun_tpu_coeff(j)), real(prof%shd_tpu_coeff(j)), &
            prof%gpp_O2(j),prof%dGOPdz_sun(j), prof%dGOPdz_shd(j), prof%Ja_sun(j),prof%Jglu_sun(j),prof%JBusch_sun(j),&
            prof%Ja_shd(j),prof%Jglu_shd(j),prof%JBusch_shd(j), &
            prof%jphoton_sun(j), prof%jphoton_shd(j),prof%sun_quad(j), prof%shd_quad(j), prof%sun_ci(j), prof%shd_ci(j)

!if (j==40) then
!    print *, 'output:', prof%sun_quad(j)
!end if


!print *, prof%dLEdz(j,1)
if (ISNAN(prof%dLEdz(j,1))) then
 !print *, prof%dLEdz(j,1)
 else
 !   print *, "normal"
end if
    end do
    if (ierr > 0) call error_writing(isroutine, noutflux, ' - 2')

     ! Profile tpu
    write(form1,'(A,I3,A)') '(i07,",",i03,', 43, '(",",es22.14))'
    do j=1, ncl
       write(nouttpu,form1,iostat=ierr) &
            time%daytime,j,&
            prof%dLAIdz(j), prof%co2_air_filter(j), &
            solar%prob_beam(j), solar%prob_shd(j), &
            prof%sun_ci(j), prof%shd_ci(j), &
            prof%sun_cs(j), prof%shd_cs(j), &
            prof%sun_cc(j), prof%shd_cc(j), &
            prof%sun_wc(j), prof%shd_wc(j), &
            prof%sun_wj(j), prof%shd_wj(j), &
            prof%sun_wp(j), prof%shd_wp(j), &
            prof%sun_alphag(j), prof%shd_alphag(j), &
            prof%sun_alphas(j), prof%shd_alphas(j), &
            prof%Ja_sun(j),     prof%Ja_shd(j), &
            prof%Jglu_sun(j),   prof%Jglu_shd(j), &
            prof%JBusch_sun(j), prof%JBusch_shd(j), &
            prof%sun_Ndemand(j),prof%shd_Ndemand(j), &
            prof%sun_Ntot(j), prof%shd_Ntot(j), &
            prof%sun_ABusch(j), prof%shd_ABusch(j), &
            prof%sun_NO3(j), prof%shd_NO3(j), &
            prof%sun_NO2(j), prof%shd_NO2(j), &
            prof%sun_NH4(j), prof%shd_NH4(j), &
            prof%sun_tpu_coeff(j), prof%shd_tpu_coeff(j)
            !prof%sun_wp(j), prof%shd_wp(j), &

!print *, prof%sun_tpu_coeff, prof%shd_tpu_coeff
    end do
    if (ierr > 0) call error_writing(isroutine, nouttpu, ' - 2')

     ! Profile nitrogen ass
    write(form1,'(A,I3,A)') '(i07,",",i03,', 11, '(",",es22.14))'
    do j=1, ncl
       write(noutnit,form1,iostat=ierr) &
            time%daytime,j,&
            prof%dNdemanddz(j),prof%dNdemanddz_sun(j),prof%dNdemanddz_shd(j),&
            prof%dNtotdz(j),prof%dNtotdz_sun(j),prof%dNtotdz_shd(j),&
            prof%dNBuschdz(j), prof%dNBuschdz_sun(j), prof%dNBuschdz_shd(j)

    end do
    if (ierr > 0) call error_writing(isroutine, noutnit, ' - 2')
    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Profile air
       write(form1,'(A,I3,A)') '(i07,",",i03,', 3, '(",",es22.14))'
       do j=1, ntl
          write(noutcisoprof,form1,iostat=ierr) &
               time%daytime, j, prof%R13_12_air(j), &
               prof%d13Cair(j), prof%c13cnc(j)
       end do
       if (ierr > 0) call error_writing(isroutine, noutcisoprof, ' - 3')
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso profile air
       write(form1,'(A,I3,A)') '(i07,",",i03,', 1*(nwiso-1), '(",",es22.14))'
       do j=1, ntl
          tmp1 = prof%rhov_air(j,2:nwiso)
          write(noutwisoprof,form1,iostat=ierr) &
               time%daytime, j, tmp1
       end do
       if (ierr > 0) call error_writing(isroutine, noutwisoprof, ' - 4')
       ! Iso profile flux
       write(form1,'(A,I3,A)') '(i07,",",i03,', 9*(nwiso-1), '(",",es22.14))'
       do j=1, ncl
          tmp1 = prof%dLEdz(j,2:nwiso)
          tmp2 = prof%sun_craig(j,2:nwiso)
          tmp3 = prof%shd_craig(j,2:nwiso)
          tmp4 = prof%sun_leafwater_e(j,2:nwiso)
          tmp5 = prof%shd_leafwater_e(j,2:nwiso)
          tmp6 = prof%sun_leafwater(j,2:nwiso)
          tmp7 = prof%shd_leafwater(j,2:nwiso)
          tmp8 = prof%sun_rtrans(j,2:nwiso)
          tmp9 = prof%shd_rtrans(j,2:nwiso)
          write(noutwisoflux,form1,iostat=ierr) &
               time%daytime, j, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, &
               tmp7, tmp8, tmp9
       end do
       if (ierr > 0) call error_writing(isroutine, noutwisoflux, ' - 5')
    end if

  END SUBROUTINE write_prof_text

  ! ------------------------------------------------------------------

END MODULE io_text
