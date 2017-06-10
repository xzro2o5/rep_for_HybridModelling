MODULE io_text

  ! This module contains routines for ascii input/output of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4, i8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lastin

  !PRIVATE :: close_text_disp
  PUBLIC :: close_text_in
  PUBLIC :: close_text_out
  !PRIVATE :: error_opening
  !PRIVATE :: error_reading
  !PRIVATE :: error_writing
  !PRIVATE :: open_text_disp
  PUBLIC :: open_text_in
  PUBLIC :: open_text_out
  PUBLIC :: read_text_disp
  PUBLIC :: read_text_in
  PUBLIC :: skip_text_in
  PUBLIC :: write_text_daily
  PUBLIC :: write_text_out
  PUBLIC :: write_text_prof

  INTEGER(i4) :: lastin

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE close_text_disp()
    ! Close dispersion file
    USE constants,  ONLY: nindisp

    IMPLICIT NONE

    close(nindisp)

  END SUBROUTINE close_text_disp


  ! ------------------------------------------------------------------
  SUBROUTINE close_text_in()
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

  END SUBROUTINE close_text_in


  ! ------------------------------------------------------------------
  SUBROUTINE close_text_out()
    ! Calls the routines to open input and output files
    USE constants, ONLY: noutseas, noutopti, noutsoil, noutdaily, noutprof, &
         noutflux, nouth2osoil, &
         noutcisodaily, noutcisoseason, noutcisoprof, &
         noutwisoleaf, noutwisoprof, noutwisoflux, noutwisosoil
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    close(unit=noutseas)
    close(unit=noutopti)
    close(unit=noutsoil)
    close(unit=noutdaily)
    close(unit=noutprof)
    close(unit=noutflux)
    close(unit=nouth2osoil)
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

  END SUBROUTINE close_text_out


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
  SUBROUTINE open_text_disp()
    ! Opens the ascii input files
    USE constants,  ONLY: nindisp
    USE parameters, ONLY: indir, dispfile

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_TEXT_DISP'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr

    ! Met
    write(stmp,'(a,a,a)') trim(indir), '/', trim(dispfile)
    open(unit=nindisp, file=stmp,action="read", status="old", &
         form="formatted", iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)

  END SUBROUTINE open_text_disp


  ! ------------------------------------------------------------------
  SUBROUTINE open_text_in()
    ! Opens the ascii input files
    USE constants,  ONLY: ninmet, ninlai, ninwiso
    USE parameters, ONLY: indir, metinfile, laiinfile, wisoinfile, extra_nate
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_TEXT_IN'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr

    ! Met
    write(stmp,'(a,a,a)') trim(indir), '/', trim(metinfile)
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

  END SUBROUTINE open_text_in


  ! ------------------------------------------------------------------
  SUBROUTINE open_text_out()
    ! Calls the routines to open input and output files
    USE constants, ONLY: &
         noutseas, noutopti, noutsoil, noutdaily, noutprof, &                ! units
         noutflux, nouth2osoil, &
         noutcisodaily, noutcisoseason, noutcisoprof, &
         noutwisoleaf, noutwisoprof, noutwisoflux, noutwisosoil, &
         dailyfile, daily13cfile, hourlyfile, hourly13cfile, optimisefile, & ! filenames
         profilefile, profile13cfile, profileisofile, fluxprofilefile, &
         fluxprofileisofile, soilfile, h2osoilfile, &
         h2osoilisofile, h2oleafisofile
    USE parameters, ONLY: outdir, outsuffix
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'OPEN_TEXT_OUT'
    CHARACTER(len=256) :: stmp
    INTEGER :: ierr

    ierr = 0
    ! Hourly/Season
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(hourlyfile), trim(outsuffix)
    open(unit=noutseas, file=stmp,action="write", status="replace", &
         form="formatted", recl=150*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutseas,*) "daytime ", "netrn ", "sumrn ", "sumh ", "sumle ", "canps ", &
         "gpp ", "transpiration ", "canresp ", "soilresp ", "boleresp ", "gsoil ", &
         "ksi ", "Tleaf ", "Tlsun ", "Tlshade ", "disc13 ", "disc13_long ", "disc13_a ", &
         "disc13_ab ", "disc13_asal ", "disc13_b ", "CsCa ", "CiCa ", "CcCa ", "gm ", "gs ", &
         "ave_dair_13C ", "isoprene%flux ", "diff_par ", "par ", "ta ", "u ", "ustar ", "rhov ", &
         "zL ", "press%Pa ", "relative%humidity ", "CO2air ", "Tair40 ", "CO2air40 ", "Vpd40 ", &
         "Tleaf40 ", "Sun_A40 ", "Sun_gs_mol40 ", "Sun_rbCO240 ", "Sun_cs40 ", "Sun_ci40 ", &
         "Sun_cc40 ", "Sun_wj40 ", "Sun_wc40 ", "Quantum_sun40 ", "Sun_rh40 ", "Sun_vpd40 ", &
         "H_40 ", "Trans_40 ", "canps_40 ", "Tleaf_40 ", "disc13_40 ", "disc13C_long_40 ", &
         "disc13C_a_40 ", "disc13C_ab_40 ", "disc13C_asl_40 ", "disc13C_b_40 ", "cica_40 ", &
         "gs_40 ", "Fisoprene_40 ", "PARdirect_40 ", "PARdiffuse_40 ", "NIRdirect_40 ", &
         "NIRdiffuse_40 ", "H_30 ", "Trans_30 ", "canps_30 ", "Tleaf_30 ", "disc13_30 ", &
         "disc13C_long_30 ", "disc13C_a_30 ", "disc13C_ab_30 ", "disc13C_asl_30 ", &
         "disc13C_b_30 ", "cica_30 ", "gs_30 ", "Fisoprene_30 ", "PARdirect_30 ", &
         "PARdiffuse_30 ", "NIRdirect_30 ", "NIRdiffuse_30 ", "H_20 ", "Trans_20 ", &
         "canps_20 ", "Tleaf_20 ", "disc13_20 ", "disc13C_long_20 ", "disc13C_a_20 ", &
         "disc13C_ab_20 ", "disc13C_asl_20 ", "disc13C_b_20 ", "cica_20 ", "gs_20 ", &
         "Fisoprene_20 ", "PARdirect_20 ", "PARdiffuse_20 ", "NIRdirect_20 ", "NIRdiffuse_20 ", &
         "H_10 ", "Trans_10 ", "canps_10 ", "Tleaf_10 ", "disc13_10 ", "disc13C_long_10 ", &
         "disc13C_a_10 ", "disc13C_ab_10 ", "disc13C_asl_10 ", "disc13C_b_10 ", "cica_10 ", &
         "gs_10 ", "Fisoprene_10 ", "PARdirect_10 ", "PARdiffuse_10 ", "NIRdirect_10 ", &
         "NIRdiffuse_10 ", "H_1 ", "Trans_1 ", "canps_1 ", "Tleaf_1 ", "disc13_1 ", &
         "disc13C_long_1 ", "disc13C_a_1 ", "disc13C_ab_1 ", "disc13C_asl_1 ", "disc13C_b_1 ", &
         "cica_1 ", "gs_1 ", "Fisoprene_1 ", "PARdirect_1 ", "PARdiffuse_1 ", "NIRdirect_1 ", &
         "NIRdiffuse_1"

    ! Optimise
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(optimisefile), trim(outsuffix)
    open(unit=noutopti, file=stmp,action="write", status="replace", &
         form="formatted", recl=2*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutopti,*) "daytime ", "soil_mm_50"

    ! Soil
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(soilfile), trim(outsuffix)
    open(unit=noutsoil, file=stmp,action="write", status="replace", &
         form="formatted", recl=25*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutsoil,*) "daytime ", "netrn ", "soilh ", "soille ", "soil%tsrf ", &
         "soil%T%base ", "soil%T%15cm ", "flux%transp ", "flux%evap ", "prof%throughfall ", &
         "soil%soil_mm ", "soil%qinfl ", "soil%qdrai ", "soil%gsoil ", "soil%qtran ", "soil%surfrun ", &
         "Tsoil%1 ", "Tsoil%2 ", "Tsoil%3 ", "Tsoil%4 ", "Tsoil%5 ", "Tsoil%6 ", "Tsoil%7 ", "Tsoil%8 ", &
         "Tsoil%9 ", "Tsoil%10"

    ! Daily average
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(dailyfile), trim(outsuffix)
    open(unit=noutdaily, file=stmp,action="write", status="replace", &
         form="formatted", recl=14*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutdaily,*) "Day ", "Avg_FC ", "Avg_EVAP ", "AVG_H ", "Avg_PAR ", "Avg_RNET ", "lai ", "Avg_PS ", &
         "Ave_Resp ", "Avg_BOLE ", "Avg_SOIL ", "Avg_TLeaf ", "Avg_Gs ", "Tleaf_day"

    ! Profile air
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profilefile), trim(outsuffix)
    open(unit=noutprof, file=stmp,action="write", status="replace", &
         form="formatted", recl=5*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutprof,*) "Time ", "i ", "tair ", "qair ", "co2"

    ! Profile fluxes
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(fluxprofilefile), trim(outsuffix)
    open(unit=noutflux, file=stmp,action="write", status="replace", &
         form="formatted", recl=40*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(noutflux,*) "daytime ", "i ", "dHdz ", "dLEdz ", "dLEdz%sun ", "dLEdz%shd ", &
         "sun_A ", "shd_A ", "dGPPdz ", "dGPPdz%sun ", "dGPPdz%shd ", "dPsdz ", &
         "dPsdz%sun ", "dPsdz%shd ", "dRESPdz ", "dRESPdz%sun ", "dRESPdz%shd ", &
         "PARdirect ", "PARdiffuse ", "Tleaf ", "Tleaf_sun ", "Tleaf_shd ", &
         "Beam ", "Nonbeam ", "LAI ", "Isopreneflux ", "gs ", "gs%sun ", "gs%shd ", &
         "disc13Clong ", "disc13Clong%sun ", "disc13Clong%shd ", "disc13C ", &
         "disc13C%sun ", "disc13C%shd ", "disc13Ca%sun ", "disc13Cab%sun ", &
         "disc13Casal%sun ", "disc13Cb%sun ", "disc13Ca%shd ", "disc13Cab%shd ", &
         "disc13Casal%shd ", "disc13Cb%shd ", "PARin"

    ! Soil water
    write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2osoilfile), trim(outsuffix)
    open(unit=nouth2osoil, file=stmp,action="write", status="replace", &
         form="formatted", recl=40*25, iostat=ierr)
    if (ierr > 0) call error_opening(isroutine, stmp)
    write(nouth2osoil,*) "Time ", "wiso%lost1 ", "flux%litterevap1 ", &
         "flux%soilevap1 ", "flux%evap1 ", "flux%c_evaporation1 ", &
         "flux%c_transpiration1 ", "flux%c_evapotranspiration1 ", &
         "flux%evapotranspiration1 ", "prof%rhov_air11 ", "prof%rhov_air1201 ", &
         "input%ppt1 ", "prof%throughfall1 ", "flux%surfrun1 ", "soil%qdrai1 ", &
         "soil%thetal1 ", "soil%theta01 ", "soil%theta11 ", "soil%theta21 ", &
         "soil%theta31 ", "soil%theta41 ", "soil%theta51 ", "soil%theta61 ", &
         "soil%theta71 ", "soil%theta81 ", "soil%theta91"

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Daily average
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(daily13cfile), trim(outsuffix)
       open(unit=noutcisodaily, file=stmp,action="write", status="replace", &
            form="formatted", recl=3*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutcisodaily,*) "Day ", "Avg_Fc_13C ", "Avg_disc13_day"
       ! Hourly/Season
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(hourly13cfile), trim(outsuffix)
       open(unit=noutcisoseason, file=stmp,action="write", status="replace", &
            form="formatted", recl=6*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutcisoseason,*) "daytime ", "disc13 ", "disc13_long ", "ave_dair_13C ", &
            "disc13Clong40 ", "disc13C40"
       ! Profile air
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profile13cfile), trim(outsuffix)
       open(unit=noutcisoprof, file=stmp,action="write", status="replace", &
            form="formatted", recl=5*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutcisoprof,*) "Time ", "i ", "R13C ", "del13C ", "13C"
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso leaf
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2oleafisofile), trim(outsuffix)
       open(unit=noutwisoleaf, file=stmp,action="write", status="replace", &
            form="formatted", recl=7*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutwisoleaf,*) "Time ", "soil%rxylem2 ", "soil%rxylem3 ", "soil%rxylem4 ", &
            "prof%rvapour12 ", "prof%rvapour13 ", "prof%rvapour14"
       ! Iso profile air
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(profileisofile), trim(outsuffix)
       open(unit=noutwisoprof, file=stmp,action="write", status="replace", &
            form="formatted", recl=5*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutwisoprof,*) "Time ", "i ", "qair2 ", "qair3 ", "qair4"
       ! Iso profile flux
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(fluxprofileisofile), trim(outsuffix)
       open(unit=noutwisoflux, file=stmp,action="write", status="replace", &
            form="formatted", recl=28*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutwisoflux,*) "daytime ", "i ", "dLEdz2 ", "dLEdz3 ", "dLEdz4 ", &
            "sun_craig2 ", "sun_craig3 ", "sun_craig4 ", "shd_craig2 ", "shd_craig3 ", &
            "shd_craig4 ", "sun_leafwater_e2 ", "sun_leafwater_e3 ", "sun_leafwater_e4 ", &
            "shd_leafwater_e2 ", "shd_leafwater_e3 ", "shd_leafwater_e4 ", &
            "sun_leafwater2 ", "sun_leafwater3 ", "sun_leafwater4 ", "shd_leafwater2 ", &
            "shd_leafwater3 ", "shd_leafwater4 ", "sun_trans2 ", "sun_trans3 ", &
            "sun_trans4 ", "shd_trans2 ", "shd_trans3 ", "shd_trans4"
       ! Iso soil
       write(stmp,'(a,a,a,a)') trim(outdir), '/', trim(h2osoilisofile), trim(outsuffix)
       open(unit=noutwisosoil, file=stmp,action="write", status="replace", &
            form="formatted", recl=128*25, iostat=ierr)
       if (ierr > 0) call error_opening(isroutine, stmp)
       write(noutwisosoil,*) "Time ", "wiso%lost1 ", "wiso%lost2 ", "wiso%lost3 ", &
            "wiso%lost4 ", "flux%litterevap1 ", "flux%litterevap2 ", "flux%litterevap3 ", &
            "flux%litterevap4 ", "flux%soilevap1 ", "flux%soilevap2 ", "flux%soilevap3 ", &
            "flux%soilevap4 ", "flux%s_evap1 ", "flux%s_evap2 ", "flux%s_evap3 ", "flux%s_evap4 ", &
            "flux%c_evaporation1 ", "flux%c_evaporation2 ", "flux%c_evaporation3 ", &
            "flux%c_evaporation4 ", "flux%c_transpiration1 ", "flux%c_transpiration2 ", &
            "flux%c_transpiration3 ", "flux%c_transpiration4 ", "flux%c_evapotranspiration1 ", &
            "flux%c_evapotranspiration2 ", "flux%c_evapotranspiration3 ", "flux%c_evapotranspiration4 ", &
            "flux%evapotranspiration1 ", "flux%evapotranspiration2 ", "flux%evapotranspiration3 ", &
            "flux%evapotranspiration4 ", "prof%rhov_air11 ", "prof%rhov_air12 ", "prof%rhov_air13 ", &
            "prof%rhov_air14 ", "prof%rhov_air1201 ", "prof%rhov_air1202 ", "prof%rhov_air1203 ", &
            "prof%rhov_air1204 ", "input%ppt1 ", "input%ppt2 ", "input%ppt3 ", "input%ppt4 ", &
            "prof%throughfall1 ", "prof%throughfall2 ", "prof%throughfall3 ", "prof%throughfall4 ", &
            "flux%surfrun1 ", "flux%surfrun2 ", "flux%surfrun3 ", "flux%surfrun4 ", "soil%qdrai1 ", &
            "soil%qdrai2 ", "soil%qdrai3 ", "soil%qdrai4 ", "soil%thetal1 ", "soil%thetal2 ", &
            "soil%thetal3 ", "soil%thetal4 ", "soil%theta01 ", "soil%theta02 ", "soil%theta03 ", &
            "soil%theta04 ", "soil%theta11 ", "soil%theta12 ", "soil%theta13 ", "soil%theta14 ", &
            "soil%theta21 ", "soil%theta22 ", "soil%theta23 ", "soil%theta24 ", "soil%theta31 ", &
            "soil%theta32 ", "soil%theta33 ", "soil%theta34 ", "soil%theta41 ", "soil%theta42 ", &
            "soil%theta43 ", "soil%theta44 ", "soil%theta51 ", "soil%theta52 ", "soil%theta53 ", &
            "soil%theta54 ", "soil%theta61 ", "soil%theta62 ", "soil%theta63 ", "soil%theta64 ", &
            "soil%theta71 ", "soil%theta72 ", "soil%theta73 ", "soil%theta74 ", "soil%theta81 ", &
            "soil%theta82 ", "soil%theta83 ", "soil%theta84 ", "soil%theta91 ", "soil%theta92 ", &
            "soil%theta93 ", "soil%theta94"
    end if

  END SUBROUTINE open_text_out


  ! ------------------------------------------------------------------
  SUBROUTINE read_text_disp()
    ! Input data on Thomson dispersion matrix that was computed offline with
    ! MOVOAK.C, Dij (s m-1)
    USE constants, ONLY: nindisp
    USE types,     ONLY: met
    USE setup,     ONLY: ncl, ntl

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'READ_TEXT_DISP'
    INTEGER :: ierr
    INTEGER(i4) :: i, j
    REAL :: in1, in2

    call open_text_disp()

    ierr = 0
    do j=1, ncl
       do i=1, ntl
          read(nindisp,*,iostat=ierr) in1, in2
          met%dispersion(i,j) = real(in1,kind=wp)
       end do
    end do
    if (ierr > 0) call error_reading(isroutine, nindisp)
    if (ierr < 0) call error_reading(isroutine, nindisp, 'reached EOF.')

    call close_text_disp()

  END SUBROUTINE read_text_disp


  ! ------------------------------------------------------------------
  SUBROUTINE read_text_in()
    ! input data and check for bad data
    ! note that the data were produced in single precision (float)
    ! so I had to read them as single precision, otherwise I ingested garbage
    USE constants,     ONLY: ninmet, ninwiso, ninlai, zero, one, two, e1, e3, &
         TN0, rugc, vonKarman, mass_air
    USE types,         ONLY: met, input, solar, time, srf_res, wiso, iswitch
    USE setup,         ONLY: ncl, nwiso
    USE parameters,    ONLY: extra_nate, bprime, zm, zd, z0, end_run, lai
    USE utils,         ONLY: es
    USE isotope_utils, ONLY: invdelta1000_h2o

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'READ_TEXT_IN'
    INTEGER :: ierr
    INTEGER(i4) :: i
    REAL :: dayy, hhrr, in01, in02, in03, in04, in05, in06, in07, in08
    REAL :: in09, in10, in11, in12, in13, in14
    INTEGER(i8) :: dt
    INTEGER(i4), DIMENSION(nwiso-1) :: mc

    dt = int(time%time_step,kind=i8) * 100_i8 / 3600_i8
    ! use input%dayy instead of time%days because of several years at once
    time%jdold = input%dayy ! identify previous day
    read(ninmet,*,iostat=ierr) dayy, hhrr, in01, in02, in03, in04, &
         in05, in06, in07, in08, in09, in10, in11, in12, in13, in14
    if (ierr > 0) call error_reading(isroutine, ninmet)
    if (ierr < 0) call error_reading(isroutine, ninmet, 'reached EOF.')
    input%dayy         = int(dayy,kind=i4)
    input%hhrr         = real(hhrr,kind=wp)
    input%ta           = real(in01,kind=wp)
    input%rglobal      = real(in02,kind=wp)
    input%parin        = real(in03,kind=wp)
    input%pardif       = real(in04,kind=wp)
    input%ea           = real(in05,kind=wp)
    input%wnd          = real(in06,kind=wp)
    input%ppt(1)       = real(in07,kind=wp)
    input%co2air       = real(in08,kind=wp)
    input%press_mb     = real(in09,kind=wp)
    input%tsoil        = real(in10,kind=wp)
    input%soilmoisture = real(in11,kind=wp)
    input%flag         = nint(in12,kind=i4)
    input%d13CO2       = real(in13,kind=wp)
    input%d18CO2       = real(in14,kind=wp)

    ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
    if (extra_nate == 1) input%lai = input%pardif ! redundant with fptr15
    ! preliminary approach to calculate consecutive year
    time%year = time%year0 + (input%dayy-1)/365
    time%daytime = int(input%dayy,kind=i8)*10000_i8 + int(input%hhrr*100._wp,kind=i8) ! define daytime
    time%doy = mod(input%dayy,365)
    if (time%doy == 0) time%doy = 365
    time%local_time = input%hhrr
    time%days = time%doy
    ! compute derived quantities for the model
    met%T_Kelvin = input%ta + TN0 ! compute absolute air temperature
    met%rhova_g = input%ea * 2165._wp/met%T_Kelvin ! compute absolute humidity, g m-3
    ! limit humidity
    if (met%rhova_g < zero) met%rhova_g = zero
    if (met%rhova_g > 30._wp) met%rhova_g = 30._wp
    met%rhova_kg = met%rhova_g *e3 ! absolute humidity, kg m-3
    met%press_kpa = input%press_mb *e1 ! air pressure, kPa
    met%press_bars = input%press_mb * e3 ! air pressure, bars
    met%press_Pa = met%press_kpa*1000._wp ! pressure, Pa
    met%relative_humidity = input%ea*10._wp/es(met%T_Kelvin) ! relative humidity
    ! combining gas law constants
    met%pstat273 = rugc / (100000._wp * met%press_bars)
    ! cuticular conductance adjusted for pressure and T, mol m-2 s-1
    ! cuticular resistance
    srf_res%rcuticle(1:ncl) = one / (bprime(1:ncl) * met%T_Kelvin * met%pstat273)
    ! check for bad CO2 data
    if (abs(input%co2air) >= 998._wp) input%co2air = 370._wp
    if (input%parin < zero) input%parin = zero ! check for bad par data
    if (input%parin <= zero) then ! check for night
       solar%par_beam = zero
       solar%par_diffuse = zero
       solar%nir_beam = zero
       solar%nir_diffuse = zero
    end if
    ! set some limits on bad input data to minimize the model from blowing up
    if (solar%ratrad > 0.9_wp .or. solar%ratrad < 0.2_wp) solar%ratrad = 0.5_wp
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
      input%ppt(2:nwiso) = input%ppt(2:nwiso) * invdelta1000_h2o(input%dppt(2:nwiso), mc(1:nwiso-1))
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

 END SUBROUTINE read_text_in


  ! ------------------------------------------------------------------
  SUBROUTINE skip_text_in()
    ! Skip input lines
    USE constants,  ONLY: ninmet, ninlai, ninwiso
    USE parameters, ONLY: extra_nate, start_run
    USE types,      ONLY: iswitch

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'SKIP_TEXT_IN'
    INTEGER :: ierr
    LOGICAL :: notreached
    REAL :: dayy, hhrr, in01, in02, in03, in04, in05, in06, in07, in08
    REAL :: in09, in10, in11, in12, in13, in14
    INTEGER(i8) :: daytime

    lastin = 0

    ierr = 0
    notreached = .true.
    do while ((ierr==0) .and. notreached)
       read(ninmet,*,iostat=ierr) dayy, hhrr, in01, in02, in03, in04, &
            in05, in06, in07, in08, in09, in10, in11, in12, in13, in14
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

  END SUBROUTINE skip_text_in


  ! ------------------------------------------------------------------
  SUBROUTINE write_text_daily()
    ! Write out daily means
    USE constants, ONLY: noutdaily, noutcisodaily
    USE types,     ONLY: iswitch, time, output, input

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_TEXT_DAILY'
    INTEGER :: ierr

    ierr = 0
    ! Daily average
    write(noutdaily,*,iostat=ierr) &
         input%dayy, output%sumfc, output%sumevap, output%sumsens, &
         output%sumpar, output%sumnet, time%lai, &
         output%sumps, output%sumresp, output%sumbole, output%sumsoil, &
         output%sumta, output%sumgs, output%sumTleaf
    if (ierr > 0) call error_writing(isroutine, noutdaily)

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Daily average
       write(noutcisodaily,*,iostat=ierr) &
            input%dayy, output%sumF13C, output%sumdisc13C
       if (ierr > 0) call error_writing(isroutine, noutcisodaily)
    end if ! 13C

  END SUBROUTINE write_text_daily


  ! ------------------------------------------------------------------
  SUBROUTINE write_text_out()
    ! Writes output except the profile
    USE constants, ONLY: noutseas, noutopti, noutsoil, nouth2osoil, &
         noutcisoseason, noutwisoleaf, noutwisosoil
    USE types,     ONLY: iswitch, time, prof, soil, bole, &
         input, output, flux, met, solar, wiso
    USE setup,     ONLY: ncl, ntl, nwiso, nsoil

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_TEXT_OUT'
    INTEGER :: ierr

    ierr = 0
    ! Hourly/Season
    write(noutseas,*,iostat=ierr) &
         time%daytime, output%netrad, output%sumrn, output%sumh, output%sumle, output%can_ps_mol, &
         output%can_gpp, output%c_transpiration_mole, output%canresp, soil%respiration_mole, &
         bole%respiration_mole, soil%gsoil, output%sumksi, output%tleaf_mean, output%tavg_sun, &
         output%tavg_shd, output%ave_disc13, output%ave_disc13_long, output%ave_disc13_a, &
         output%ave_disc13_ab, output%ave_disc13_asal, output%ave_disc13_b, output%ave_csca, &
         output%ave_cica, output%ave_ccca, output%ave_gm, output%ave_gs, output%ave_daC13, &
         output%isoprene_efflux, output%diff_par, input%parin,input%ta, input%wnd, met%ustar, &
         met%rhova_g, met%zl, met%press_Pa, met%relative_humidity, input%co2air, prof%tair(ncl), &
         prof%co2_air(ncl), prof%vpd_air(ncl), prof%tleaf(ncl), prof%sun_A(ncl), &
         prof%sun_gs_mol(ncl), prof%sun_rbco2(ncl), prof%sun_cs(ncl), prof%sun_ci(ncl), &
         prof%sun_cc(ncl), prof%sun_wj(ncl), prof%sun_wc(ncl), solar%quantum_sun(ncl), &
         prof%sun_rh(ncl), prof%sun_vpd(ncl), prof%dHdz(ncl), prof%dLEdz(ncl,1), prof%dPsdz(ncl), &
         prof%tleaf(ncl), prof%disc13C(ncl), prof%disc13C_long(ncl), prof%disc13C_a(ncl), &
         prof%disc13C_ab(ncl), prof%disc13C_asal(ncl), prof%disc13C_b(ncl), prof%cica(ncl), &
         prof%dStomCondz(ncl), prof%isopreneflux(ncl), solar%beam_flux_par(ncl), &
         solar%par_down(ncl), solar%beam_flux_nir(ncl), solar%nir_dn(ncl), prof%dHdz(3*ncl/4), &
         prof%dLEdz(3*ncl/4,1), prof%dPsdz(3*ncl/4), prof%tleaf(3*ncl/4), prof%disc13C(3*ncl/4), &
         prof%disc13C_long(3*ncl/4), prof%disc13C_a(3*ncl/4), prof%disc13C_ab(3*ncl/4), &
         prof%disc13C_asal(3*ncl/4), prof%disc13C_b(3*ncl/4), prof%cica(3*ncl/4), &
         prof%dStomCondz(3*ncl/4), prof%isopreneflux(3*ncl/4), solar%beam_flux_par(3*ncl/4), &
         solar%par_down(3*ncl/4), solar%beam_flux_nir(3*ncl/4), solar%nir_dn(3*ncl/4), &
         prof%dHdz(ncl/2), prof%dLEdz(ncl/2,1), prof%dPsdz(ncl/2), prof%tleaf(ncl/2), &
         prof%disc13C(ncl/2), prof%disc13C_long(ncl/2), prof%disc13C_a(ncl/2), &
         prof%disc13C_ab(ncl/2), prof%disc13C_asal(ncl/2), prof%disc13C_b(ncl/2), prof%cica(ncl/2), &
         prof%dStomCondz(ncl/2), prof%isopreneflux(ncl/2), solar%beam_flux_par(ncl/2), &
         solar%par_down(ncl/2), solar%beam_flux_nir(ncl/2), solar%nir_dn(ncl/2), prof%dHdz(ncl/4), &
         prof%dLEdz(ncl/4,1), prof%dPsdz(ncl/4), prof%tleaf(ncl/4), prof%disc13C(ncl/4), &
         prof%disc13C_long(ncl/4), prof%disc13C_a(ncl/4), prof%disc13C_ab(ncl/4), &
         prof%disc13C_asal(ncl/4), prof%disc13C_b(ncl/4), prof%cica(ncl/4), prof%dStomCondz(ncl/4), &
         prof%isopreneflux(ncl/4), solar%beam_flux_par(ncl/4), solar%par_down(ncl/4), &
         solar%beam_flux_nir(ncl/4), solar%nir_dn(ncl/4), prof%dHdz(1), prof%dLEdz(1,1), &
         prof%dPsdz(1), prof%tleaf(1), prof%disc13C(1), prof%disc13C_long(1), prof%disc13C_a(1), &
         prof%disc13C_ab(1), prof%disc13C_asal(1), prof%disc13C_b(1), prof%cica(1), &
         prof%dStomCondz(1), prof%isopreneflux(1), solar%beam_flux_par(1), solar%par_down(1), &
         solar%beam_flux_nir(1), solar%nir_dn(1)

    if (ierr > 0) call error_writing(isroutine, noutseas)

    ! Optimise
    write(noutopti,*,iostat=ierr) time%daytime, soil%soil_mm_50
    if (ierr > 0) call error_writing(isroutine, noutopti)

    ! Soil
    write(noutsoil,*,iostat=ierr) &
         time%daytime, output%rnet_soil, soil%heat, soil%evap, &
         soil%tsrf, soil%T_base, soil%T_15cm, flux%c_transpiration(1), &
         flux%s_evap(1), prof%throughfall(1,1), soil%soil_mm, &
         soil%qinfl(1), soil%qdrai(1), output%c7,soil%qtran(1), &
         flux%surfrun(1), soil%T_soil(1:nsoil)
    if (ierr > 0) call error_writing(isroutine, noutsoil)

    ! Soil water
    write(nouth2osoil,*,iostat=ierr) &
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
       write(noutcisoseason,*,iostat=ierr) &
            time%daytime, output%ave_disc13, output%ave_disc13_long, &
            output%ave_daC13, prof%disc13C_long(ncl), prof%disc13C(ncl)
       if (ierr > 0) call error_writing(isroutine, noutcisoseason)
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso leaf
       write(noutwisoleaf,*,iostat=ierr) &
            time%daytime, soil%rxylem(2), soil%rxylem(3), soil%rxylem(4), &
            prof%rvapour(1,2), prof%rvapour(1,3), prof%rvapour(1,4)
       if (ierr > 0) call error_writing(isroutine, noutwisoleaf)
       ! Iso soil
       write(noutwisosoil,*,iostat=ierr) &
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

  END SUBROUTINE write_text_out


  ! ------------------------------------------------------------------
  SUBROUTINE write_text_prof()
    ! Calls the routines to open input and output files
    USE constants, ONLY: noutprof, noutflux, noutcisoprof, noutwisoprof, noutwisoflux
    USE types,     ONLY: iswitch, prof, solar, time, input
    USE setup,     ONLY: ncl, ntl, nwiso

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: isroutine = 'WRITE_TEXT_PROF'
    INTEGER :: ierr
    INTEGER(i4) :: j
    REAL(wp), DIMENSION(nwiso-1) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9

    ierr = 0
    ! Profile air
    do j=1, ntl
       write(noutprof,*,iostat=ierr) &
            time%daytime, j, prof%tair(j), prof%rhov_air(j,1), &
            prof%co2_air(j)
    end do
    if (ierr > 0) call error_writing(isroutine, noutprof)

    ! Profile fluxes
    do j=1, ncl
       write(noutflux,*,iostat=ierr) &
            time%daytime,j,prof%dHdz(j), prof%dLEdz(j,1), &
            prof%dLEdz_sun(j), prof%dLEdz_shd(j), prof%sun_A(j), &
            prof%shd_A(j), prof%dGPPdz(j), prof%dGPPdz_sun(j), &
            prof%dGPPdz_shd(j), prof%dPsdz(j), prof%dPsdz_sun(j), &
            prof%dPsdz_shd(j), prof%dRESPdz(j), prof%dRESPdz_sun(j), &
            prof%dRESPdz_shd(j), solar%quantum_sun(j), &
            solar%quantum_shd(j), prof%tleaf(j), prof%sun_tleaf(j), &
            prof%shd_tleaf(j), solar%prob_beam(j), solar%prob_shd(j), &
            prof%dLAIdz(j), prof%isopreneflux(j), &
            prof%dStomCondz(j)*1e3_wp,prof%dStomCondz_sun(j)*1e3_wp,prof%dStomCondz_shd(j)*1e3_wp, &
            prof%disc13C_long(j), &
            prof%sun_disc13_long(j), prof%shd_disc13_long(j), &
            prof%disc13C(j), prof%sun_disc13(j), prof%shd_disc13(j), &
            prof%sun_disc13_ab(j), prof%sun_disc13_a(j), &
            prof%sun_disc13_asal(j), prof%sun_disc13_b(j), &
            prof%shd_disc13_ab(j), prof%shd_disc13_a(j), &
            prof%shd_disc13_asal(j), prof%shd_disc13_b(j), input%parin
    end do
    if (ierr > 0) call error_writing(isroutine, noutflux)

    ! 13CO2 isotope files
    if (iswitch%d13c == 1) then
       ! Profile air
       do j=1, ntl
          write(noutcisoprof,*,iostat=ierr) &
               time%daytime, j, prof%R13_12_air(j), &
               prof%d13Cair(j), prof%c13cnc(j)
       end do
       if (ierr > 0) call error_writing(isroutine, noutcisoprof)
    end if ! 13C

    ! Water isotope files
    if (iswitch%wiso == 1) then
       ! Iso profile air
       do j=1, ntl
          tmp1 = prof%rhov_air(j,2:nwiso)
          write(noutwisoprof,*,iostat=ierr) &
               time%daytime, j, tmp1
       end do
       if (ierr > 0) call error_writing(isroutine, noutwisoprof)
       ! Iso profile flux
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
          write(noutwisoflux,*,iostat=ierr) &
               time%daytime, j, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, &
               tmp7, tmp8, tmp9
       end do
       if (ierr > 0) call error_writing(isroutine, noutwisoflux)
    end if

  END SUBROUTINE write_text_prof

END MODULE io_text
