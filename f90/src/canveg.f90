PROGRAM canveg
! start adding ROC to the model!
  ! Parameters / Types / Switches / etc.
  USE kinds,         ONLY: wp, i4, i8                          ! Precision
  USE constants,     ONLY: zero, half, one, e2, e3, e6, &      ! numerical
       version, main_file, namelist_file, &                    ! info
       TN0, P0, Gravity, rugc, cp, &                           ! physical
       mass_air, mass_CO2, o2_ref, &
       RVSMOW_18O, RVSMOW_D, Rpdb_12C, &
       nnu, nuvisc, dc, ddh, dh, ddv, dv, &
       isnight, undef                                          ! computational
  USE setup,         ONLY: ncl, ntl, nsoil, nwiso, ndaysc13    ! dimensions
  USE parameters,    ONLY: read_namelist, &                    ! namelist parameters
       longitude, latitude, ht, lai, &
       indir, outdir, workdir, start_run, end_run, &
       start_profiles, end_profiles, &
       extra_nate, lleaf, delz, bprime, &
       ROC_leaf_in, ROC_bole_in, ROC_soil_in ! Yuan 2018.01.17
  USE types,         ONLY: time, iswitch, soil, wiso, &        ! Canveg''s types
       non_dim, bound_lay_res, input, prof, solar, met, &
       ciso, fact, flux, output, bole, srf_res, debug, &
       ! Routines
       zero_new_timestep, nitrogen
  USE io,            ONLY: create_dir, open_files, close_files, &          ! I/O wrappers
       read_disp, skip_input, lastin, write_daily, read_in, &              ! Yuan 2018.01.22 to create a new directory
       write_profiles, write_output, copy_code ! copy_code, Yuan 2018.05.07

  USE soils,         ONLY: set_litter_texture, &               ! Soil & Litter
       set_soil_moisture, set_soil_temp, set_soil_root, &      ! water & energy
       set_soil_clapp, set_soil_texture, set_soil_saxton, &
       set_litter_moisture, set_litter_temp, set_soil_time, &
       set_soil_litter_capacity, soil_h2o, litter_h2o, &
       litter_rain, throughfall, soil_energy_balance
  USE energy_carbon, ONLY: bole_respiration, &                 ! Carbon & Energy
       soil_respiration, energy_and_carbon_fluxes
  USE radiation,     ONLY: set_leaf_phenology, lai_time, &     ! Phenology & Radiation
       irflux, angle, nir, par, rnet, &
       diffuse_direct_radiation, gfunc
  USE transport,     ONLY: conc, conc_seperate, friction_velocity, uz             ! Dispersion
  USE isotopes,      ONLY: leaf_wiso, soil_flux_wiso, &        ! Isotopes
       le_wiso, canopy_flux_wiso, carbon_isotopes
  USE  OXYGEN,       ONLY: OXYFLUX                             ! oxygen module Yuan 2018.01.19
  USE isoprene,      ONLY: isoprene_canopy_flux                ! Isoprene
  USE utils,         ONLY: es, lambda                          ! Utilities
  USE isotope_utils, ONLY: alpha_equ_h2o, delta1000, &
       delta1000_h2o, invdelta1000, invdelta1000_h2o, isorat
#ifndef QUIET
  USE string_utils,  ONLY: num2str, separator                  ! Print out
  USE messages,      ONLY: message, message_text
  USE finishes,      ONLY: finish
#endif

  IMPLICIT NONE

  LOGICAL judgenight ! judge if night, more complex than isnight. Yuan 2018.01.16
  COMMON  judgenight

  INTEGER(i4) :: i_count=0, i_max=0, j=0
  INTEGER(i8) :: daycnt=0, gppcnt=0
  INTEGER(i4) :: mc=0
  INTEGER(i4) :: count40=0, totalcount=0, dcount40=0, dtotalcount=0, ncount40=0, ntotalcount=0
  INTEGER(i4) :: print_balance=0
  REAL(wp) :: rnet_soil=zero, netrad=zero
  REAL(wp) :: fc_mg=zero, fc_mol=zero, fc_13C=zero
  !REAL(wp) :: wue=zero
  REAL(wp) :: sumA=zero
  REAL(wp) :: sumlai=zero ! sums for daily averages and totals
  REAL(wp) :: sumpai=zero ! YUAN 2018.03.04
  REAL(wp) :: sumh=zero, sumle=zero, sumrn=zero, sumlout=zero, sumksi=zero, sum13C=zero
  REAL(wp) :: ebalance=zero, sbalance=zero, tbalance=zero
  REAL(wp) :: sumdisc13=zero, sumdisc13_long=zero, sumdisc13_long_isotope=zero
  REAL(wp) :: ave_disc13=zero, ave_disc13_long=zero, ave_disc13_long_isotope=zero, ave_daC13=zero
  REAL(wp) :: sum_csca=zero, sum_cica=zero, sum_ccca=zero, sum_gs=zero, sum_gm=zero
  REAL(wp) :: ave_csca=zero, ave_cica=zero, ave_ccca=zero, ave_gs=zero, ave_gm=zero
  REAL(wp) :: sumdisc13_a=zero, sumdisc13_ab=zero, sumdisc13_asal=zero, sumdisc13_b=zero
  REAL(wp) :: ave_disc13_a=zero, ave_disc13_ab=zero, ave_disc13_asal=zero, ave_disc13_b=zero
  REAL(wp) :: Rlongterm_12C=zero
  REAL(wp) :: can_ps_mol=zero, can_gpp=zero, can_ps_mg=zero, canresp=zero, c_transpiration_mole=zero
  REAL(wp) :: isoprene_efflux=zero !, recycle=zero
  REAL(wp) :: etest=zero, etest_old1=zero, etest_old2=zero, etest_diff1=zero, etest_diff2=zero ! test value for energy closure
  REAL(wp) :: itest=zero, itest_old1=zero, itest_old2=zero, itest_diff1=zero, itest_diff2=zero ! test for wiso convergence
  REAL(wp) :: tleaf_mean=zero, tavg_sun=zero, tavg_shd=zero ! integrated leaf temperatures
  REAL(wp) :: ztmp=zero, temp3=zero
  REAL(wp), DIMENSION(:), ALLOCATABLE :: sun_A, shd_A ! ncl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tmpntl       ! ntl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tmpncl       ! ncl
  REAL(wp), DIMENSION(:), ALLOCATABLE :: rcws         ! nwiso

!  CHARACTER(LEN=256),PARAMETER :: stmp

  ! output of oxygen module Yuan 2018.01.30
  REAL(wp) :: can_gpp_o=zero, can_ps_o=zero, canresp_o = zero, sumo=zero, sumneto=zero, sumcanresp_o=zero, &
  sumNsupply=zero, sumNdemand=zero


#ifndef QUIET
  REAL :: ctime1=zero, ctime2=zero
  INTEGER, DIMENSION(8) :: dtime
#endif

  !
  ! Start
  !

#ifndef QUIET
  call cpu_time(ctime1)
  call date_and_time(values=dtime)
  call message(separator)
  message_text = trim(num2str(dtime(3),'(I2.2)'))//"."//trim(num2str(dtime(2),'(I2.2)')) &
       //"."//trim(num2str(dtime(1),'(I4.4)'))//" "//trim(num2str(dtime(5),'(I2.2)')) &
       //":"//trim(num2str(dtime(6),'(I2.2)'))//":"//trim(num2str(dtime(7),'(I2.2)'))
  call message('Canveg: ', 'Start at ', trim(message_text))
  call message()
  call message('    This is Canveg version ', trim(version))
  call message('    Using main file ', trim(main_file), ' and namelist ', trim(namelist_file))
#endif
! copy source code:
! Copy_Folder (trim(old_file), trim(new_file))
  !
  ! Read namelist
  !

#ifndef QUIET
  call message()
  call message('    Read namelist file: ', trim(namelist_file))
#endif
  call read_namelist() ! C: L1606-1630 READ_PARAMETER() & SET_PARAMETER()

  ! Show model setup
#ifndef QUIET
  ! The following is not possible in Fortran90.
  !    write(nout,*) '    # canopy layers:     ', num2str(ncl)
  ! because the write statement calls num2str which itself has a write statement in it.
  ! e.g. NAG complains with: Recursive I/O reference.
  ! It should be possible in Fortran2003.
  ! However, MC''s NAG compiler is not allowing it even with the -f2003 compiler flag, i.e. use call message() instead.
  call message()
  call message('    Model setup')
  call message('    Extra Nate:          ', trim(num2str(extra_nate)))
  call message('    Working directory:   ', trim(workdir))
  call message('    Input directory:     ', trim(indir))
  call message('    Output directory:    ', trim(outdir))
  call message('    # canopy layers:     ', trim(num2str(ncl)))
  call message('    # soil layers:       ', trim(num2str(nsoil)))
  call message('    # water isotopes:    ', trim(num2str(nwiso-1)))
  call message('    Longitude (N):       ', trim(num2str(longitude,'(F10.2)')))
  call message('    Latitude (E):        ', trim(num2str(latitude,'(F10.2)')))
  call message('    Canopy height (m):   ', trim(num2str(ht,'(F10.2)')))
  call message('    LAI:                 ', trim(num2str(lai,'(F10.2)')))
  call message('    Litter depth (m):    ', trim(num2str(soil%z_litter,'(F10.2)')))
  call message('    Year0:               ', trim(num2str(time%year0)))
  call message('    Start run (dddhhmm): ', trim(num2str(start_run,'(I10.8)')))
  call message('    End run (dddhhmm):   ', trim(num2str(end_run,'(I10.8)')))
  call message('    Time step (s):       ', trim(num2str(nint(time%time_step,i4))))
  call message('    Switch soil resp:    ', trim(num2str(iswitch%soil_resp_temp)))
  call message('    Switch Ball-Berry:   ', trim(num2str(iswitch%ball)))
  call message('    Switch isoprene:     ', trim(num2str(iswitch%isoprene)))
  call message('    Switch 13C:          ', trim(num2str(iswitch%d13c)))
  call message('    Switch water iso:    ', trim(num2str(iswitch%wiso)))
  call message('    Switch Bethy resp:   ', trim(num2str(iswitch%bethy_resp)))
  call message('    Extra Nate:          ', trim(num2str(extra_nate)))
#endif
  !
  ! Further model setup
  !

  ! Vienna Standard Mean Ocean Water (VSMOW) for water isotopes
  wiso%vsmow(:) = one ! for vsmow(4) = 1H1H16O; vsmow(1) not used
  if (nwiso >= 2) wiso%vsmow(2) = RVSMOW_18O ! 18O/16O
  if (nwiso >= 3) wiso%vsmow(3) = RVSMOW_D   ! 2H/1H

  ! Soil setup
  soil%d       = e6       ! displacement height for soil ~ 0.
  soil%z0      = 0.015_wp ! roughness length of soil = 1.5 cm
  soil%drain0  = 0
  soil%lost0   = 0
  ! soil%T_base is set to mean annual air temperature
  if (extra_nate==1) then
     soil%T_base = 10.8_wp ! Nate McDowells Juniper site
  else
     soil%T_base = 7.5_wp ! Hainich
  end if
  ! 'secure' T profile
  soil%T_base = min(max(soil%T_base, -10._wp), 50._wp)
  ! 0: layered soil water 1: no soil water model
  if (nsoil==0) soil%camillo = 1

  ! Misc
  i_max = 501            ! maximum number of iterations for energy balance
  time%year = time%year0 ! Save first year
  if (.not. allocated(sun_A))  allocate(sun_A(ncl))
  if (.not. allocated(shd_A))  allocate(shd_A(ncl))
  if (.not. allocated(tmpntl))   allocate(tmpntl(ntl))
  if (.not. allocated(tmpncl)) allocate(tmpncl(ncl))
  if (.not. allocated(rcws))   allocate(rcws(nwiso))

  ! create a new directory for output Yuan 2018.01.22
  call create_dir()
!  PRINT *, outdir
! copy source code
  call copy_code()
  !
  ! Open files
  call open_files()
  !
  ! Setup leaves
  call set_leaf_phenology() ! set leaf onset and full
  ! write(*,'(a,2i5)') 'CV01 ', time%leafout, time%leaffull
  ! write(*,'(a,2i5)') 'CV02 ', time%leaffall, time%leaffallcomplete
  !
  ! Setup soil
  call set_soil_texture()   ! set soil texture per soil layer
  call set_soil_root()      ! set rooting profile
  call set_soil_moisture()  ! set humidity per soil layer
  call set_soil_temp()      ! set deep soil temperature
  ! write(*,'(a,3f20.14)') 'CV03 ', soil%theta_ls, soil%theta_l33, soil%n_l
  ! write(*,'(a,2f20.14)') 'CV04 ', soil%root(1), soil%root(nsoil)
  ! write(*,'(a,2f20.14)') 'CV05 ', soil%theta(1,1), soil%theta(nsoil,1)
  ! write(*,'(a,2f20.14)') 'CV06 ', soil%T_soil(0), soil%T_soil(nsoil+1)
  ! write(*,'(a,2f20.14)') 'CV07 ', soil%T_soil_filter(0), soil%T_soil_filter(nsoil+1)

  ! Setup litter
  call set_litter_texture()
  call set_litter_temp()
  call set_litter_moisture()
  ! write(*,'(a,3f20.14)') 'CV08 ', soil%theta_ls, soil%theta_l33, soil%n_l
  ! write(*,'(a,2f20.14)') 'CV09 ', soil%T_l, soil%T_l_filter
  ! write(*,'(a,1f20.14)') 'CV10 ', soil%theta_l(1)
  if (soil%saxton==1) then  ! set hydraulic soil parameters
     call set_soil_saxton()
     ! write(*,'(a,2f20.14)') 'CV11.01 ', soil%theta_1500(1), soil%theta_1500(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.02 ', soil%theta_33(1), soil%theta_33(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.03 ', soil%theta_s33(1), soil%theta_s33(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.04 ', soil%psi_e(1), soil%psi_e(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.05 ', soil%theta_s(1), soil%theta_s(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.06 ', soil%rho(1), soil%rho(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.07 ', soil%big_b(1), soil%big_b(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.08 ', soil%big_a(1), soil%big_a(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.09 ', soil%lambda(1), soil%lambda(nsoil)
     ! write(*,'(a,2f20.14)') 'CV11.10 ', soil%k_s(1), soil%k_s(nsoil)
     ! write(*,'(a,1f20.14)') 'CV11.11 ', soil%soil_mm_33_root
     ! write(*,'(a,1f20.14)') 'CV11.12 ', soil%soil_mm_1500_root
  else
     call set_soil_clapp()
     ! write(*,'(a,2f20.14)') 'CV12.01 ', soil%theta_s(1), soil%theta_s(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.02 ', soil%k_s(1), soil%k_s(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.03 ', soil%psi_e(1), soil%psi_e(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.04 ', soil%big_b(1), soil%big_b(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.05 ', soil%theta_1500(1), soil%theta_1500(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.06 ', soil%theta_33(1), soil%theta_33(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.07 ', soil%theta_s33(1), soil%theta_s33(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.08 ', soil%rho(1), soil%rho(nsoil)
     ! write(*,'(a,2f20.14)') 'CV12.09 ', soil%k_s(1), soil%k_s(nsoil)
     ! write(*,'(a,1f20.14)') 'CV12.10 ', soil%soil_mm_33_root
     ! write(*,'(a,1f20.14)') 'CV12.11 ', soil%soil_mm_1500_root
  end if
  !
  !
  ! Relative plant available water
  soil%soil_mm_root = sum(soil%root(1:nsoil)*soil%theta(1:nsoil,1)*1000._wp * &
       (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))
  ! write(*,'(a,1f20.14)') 'CV13 ', soil%soil_mm_root
  !
  ! Constants for leaf boundary layers
  non_dim%lfddh = lleaf / ddh
  non_dim%pr    = nuvisc / dh                            ! Prandtl Number
  non_dim%pr33  = non_dim%pr**0.33_wp
  non_dim%lfddv = lleaf / ddv                            ! Diffusivity of water vapor, m2 s-1
  non_dim%sc    = nuvisc / dv                            ! Schmidt number for vapor
  non_dim%sc33  = non_dim%sc**0.33_wp
  non_dim%scc   = nuvisc / dc                            ! Schmidt number for CO2
  non_dim%scc33 = non_dim%scc**0.33_wp
  non_dim%grasshof = Gravity*lleaf*lleaf*lleaf/(nnu*nnu) ! Grasshof Number
  ! boundary layer resistances
  bound_lay_res%heat  = zero
  bound_lay_res%vapor = zero
  bound_lay_res%co2   = zero
  tmpntl(:) = zero
  tmpncl(:) = zero
  ! assign heights to array
  forall(j=1:ncl) prof%ht(j) = delz * real(j,kind=wp)
  ! write(*,'(a,3f20.14)') 'CV14 ', non_dim%lfddh, non_dim%pr, non_dim%pr33
  ! write(*,'(a,3f20.14)') 'CV15 ', non_dim%lfddv, non_dim%sc, non_dim%sc33
  ! write(*,'(a,3f26.14)') 'CV16 ', non_dim%scc, non_dim%scc33, non_dim%grasshof
  ! write(*,'(a,2f20.14)') 'CV17 ', prof%ht(1), prof%ht(ncl)

  ! Define parameters for soil energy balance model
  ! Preliminary tests show that the soil model can be sensitive to the depth of the litter layer
  call set_soil_time()
  ! write(*,'(a,1f20.14,i5)') 'CV18 ', soil%temperature_dt, soil%temperature_mtime
  ! write(*,'(a,2f20.14,i5)') 'CV19 ', soil%moisture_dt, soil%moisture_dt, soil%moisture_mtime

  ! Input data on Thomson dispersion matrix that was computed offline with
  ! MOVOAK%C, Dij (s m-1)
  call read_disp()
  ! write(*,'(a,3f20.14)') 'CV20 ', met%dispersion(1,1), met%dispersion(ntl/2,ncl/2), met%dispersion(ntl,ncl)

  ! Initialize respired 13C to the heterotrophic soil respiration value
  ! (assumed to equal longterm average), then update each day
  if (iswitch%d13c==1) then
     Rlongterm_12C         = invdelta1000(ciso%delta_soil)*Rpdb_12C
     prof%Rresp_ave(1:ncl) = Rlongterm_12C
     prof%Rresp_sum(1:ncl) = zero
     prof%cnt_Rresp(1:ncl) = 0
     ! arbitary values are taken for initialisation
     ciso%bigdelta(1:ndaysc13)      = 18._wp
     ciso%bigdelta_long(1:ndaysc13) = 18._wp
     ! write(*,'(a,3f20.14)') 'CV21 ', prof%Rresp_ave(1)
  end if

  ! Loop through the input met file. There should be a line of data for each hour of the year.
  call skip_input()
  time%days  = int(start_run / 10000_i8, i4)
  time%jdold = int(start_run / 10000_i8, i4)
  ! set input%dayy because of several years at once
  input%dayy = time%days ! for more then 1 year
  ! print *, time%days
  ! write(*,'(a,3i10)') 'CV22 ', time%days, time%jdold, input%dayy
  ! initialise arbitrary value
  solar%ratradnoon = 0.5_wp
  ! initialise interception reservoir
  prof%cws(:,:) = zero
  ! initialise leaf water at evaporating site
  forall(mc=2:nwiso)
     prof%sun_leafwater_e(1:ncl,mc)     = wiso%vsmow(mc)
     prof%shd_leafwater_e(1:ncl,mc)     = wiso%vsmow(mc)
     prof%sun_leafwater_e_old(1:ncl,mc) = wiso%vsmow(mc)
     prof%shd_leafwater_e_old(1:ncl,mc) = wiso%vsmow(mc)
  end forall
!  print *, wiso%vsmow
!  print *,prof%sun_leafwater_e(1:ncl,2:nwiso)
  ! initialise leaf area index
  input%lai = lai
  if (extra_nate==1) then
     ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
     input%lai_up   = lai
     input%lai_down = lai
  end if
  call lai_time() ! define leaf area and canopy structure
  ! write(*,'(a,1f20.14)') 'CV23 ', time%lai
  ! write(*,'(a,3f20.14)') 'CV24 ', solar%par_reflect, solar%par_trans, solar%par_soil_refl_dry
  ! write(*,'(a,3f20.14)') 'CV25 ', solar%par_absorbed, solar%nir_reflect, solar%nir_trans
  ! write(*,'(a,3f20.14)') 'CV26 ', solar%nir_soil_refl_dry, solar%nir_absorbed
  ! write(*,'(a,2f20.14)') 'CV27 ', prof%dLAIdz(1), prof%dLAIdz(ncl)
  ! write(*,'(a,2f20.14)') 'CV28 ', prof%dPAIdz(1), prof%dPAIdz(ncl)
  ! write(*,'(a,2f20.14)') 'CV29 ', solar%exxpdir(1), solar%exxpdir(ncl)

  ! Initialize humidity and temperature profiles, arbitrary values
  prof%sun_tleaf(1:ncl)        = 25._wp
  prof%shd_tleaf(1:ncl)        = 25._wp
  prof%sun_tleaf_filter(1:ncl) = 25._wp
  prof%shd_tleaf_filter(1:ncl) = 25._wp
  tmpncl(1:ncl)                = one/(bprime(1:ncl)*(rugc*TN0/P0))
  prof%sun_rs(1:ncl)           = tmpncl(1:ncl)
  prof%sun_rs_filter(1:ncl)    = tmpncl(1:ncl)
!  print *, prof%sun_rs_filter
  prof%shd_rs(1:ncl)           = tmpncl(1:ncl)
  prof%shd_rs_filter(1:ncl)    = tmpncl(1:ncl)
  ! write(*,'(a,2f23.14)') 'CV30 ', prof%sun_rs(1), prof%sun_rs(ncl)

  prof%tair(1:ntl)                 = 25._wp
!  print *, "1\n"
!  print *, prof%tair
  prof%tair_filter(1:ntl)          = 25._wp
  prof%tair_filter_save(1:ntl)     = 25._wp
  prof%rhov_air(1:ntl,1)           = 0.02_wp
  prof%rhov_air_save(1:ntl)        = 0.02_wp
  prof%rhov_air_filter(1:ntl,1)    = 0.02_wp
  prof%rhov_air_filter_save(1:ntl) = 0.02_wp
  prof%co2_air(1:ntl)              = 380._wp
  prof%co2_air_filter(1:ntl)       = 380._wp
  prof%O2_air(1:ntl)              = o2_ref
  prof%O2_air_filter(1:ntl)       = o2_ref

  ! Initialise atmospheric d13C
  if (iswitch%d13c==1) then
     prof%R13_12_air(1:ntl) = invdelta1000(-8._wp)*Rpdb_12C ! ratio of 13C relative to 12C
     prof%d13Cair(1:ntl)    = -8.0_wp
  end if
  ! write(*,'(a,2f20.14)') 'CV31 ', prof%R13_12_air(1), prof%R13_12_air(ntl)

  ! Initialise atmospheric water isotopes
  forall(mc=2:nwiso) ! 2017.10.04
     prof%rhov_air(1:ntl,mc)        = 0.02_wp*wiso%vsmow(mc)
     prof%rhov_air_filter(1:ntl,mc) = 0.02_wp*wiso%vsmow(mc)
     prof%rvapour(1:ntl,mc)         = wiso%vsmow(mc)
  end forall
!print *, prof%rhov_air(1:ntl,2:nwiso)
!print *, prof%rvapour(1,2:nwiso)
  ! initialize O2: CO2 Yuan 2018.01.17
  if (iswitch%oxygen==1) then
     !prof%ROC_leaf_G(1:ncl) =  ROC_leaf_in!
     prof%ROC_bole_air(1:ncl) =  ROC_bole_in!
     soil%ROC_soil_air =  ROC_soil_in!
!     print *, ROC_leaf_in, ROC_bole_in, ROC_soil_in
 !    print *, prof%ROC_leaf_G, prof%ROC_bole_air, soil%ROC_soil_air
  end if
  ! write(*,'(a,2f20.14)') 'CV32 ', prof%rhov_air(1,2), prof%rhov_air(ntl,3)
  ! initialize soil surface temperature with air temperature
  soil%tsrf        = 25._wp
  soil%tsrf_filter = 25._wp

  ! initialise evapotranspiration
  soil%soilevap          = zero
  soil%soilevap_filter   = zero
  soil%litterevap        = zero
  soil%litterevap_filter = zero

  ! initialise sensible heat flux
  met%H        = zero
  met%H_filter = zero

  ! -------------------------------------------------------------------
  ! MAIN PROGRAM

  ! Describe canopy attributes

  ! Flow chart of the Main Program:

  ! 1a) start run at beginning of year
  ! 1b) establish new LAI
  ! 1c) Compute continuous distrib of height at a function of lai
  ! 1d) Set soil parameters

  ! 2) input met values

  ! 3) Compute solar elevation angle

  ! 4) Compute PAR and NIR profiles and their flux densities on
  ! the sunlit and shaded leaf fractions

  ! 5) Compute sunlit and shaded leaf fractions

  ! 6) Compute first estimate of stomatal conductance

  ! 7) Compute first estimate of IRFLUX assuming leaf temperature
  ! equals air temperature.
  ! 8) Compute first estimate of leaf energy balance and leaf
  ! temperature

  ! 9) Compute photosynthesis, transpiration

  ! 10) Compute soil energy balance

  ! 11) update computation of friction velocity with new H and z/L

  ! 12) Compute new scalar and source/sink profiles for CO2, T and e

  ! 13) Iterate among 7 through 12 until convergence

  ! 14) compute fluxes of isoprene, 13C isotopes or whatever

  ! 15) Thats all Folks!!!
  ! -------------------------------------------------------------------

  ! Start time steps
#ifndef QUIET
  call message()
  if (print_balance==0) call message("Day-oy-year")
#endif
  do ! loop over time steps until lastin==1

     ! init new time step
     if (nwiso > 1) then
        prof%sun_leafwater_e_old(1:ncl,2:nwiso) = prof%sun_leafwater_e(1:ncl,2:nwiso)
        prof%shd_leafwater_e_old(1:ncl,2:nwiso) = prof%shd_leafwater_e(1:ncl,2:nwiso)
     end if
     call zero_new_timestep() ! zero variables
     call read_in()           ! input new time step
       ! print doy and hours"
     call message('    day and time:          ', trim(num2str(time%daytime)))
!    call message('    hours:          ', trim(num2str(time%daytime)))
!     print *, time%daytime
!     print *, input%co2air
!     print *, input%o2air
!    if (time%daytime==1771200) then  !Yuan 2017.09.20

!	      print *, "what an hour!!!\n"
!	      print *, "daytime = ", time%daytime
!		  print *, "\n"
!	 end if
     ! write(*,'(a,i10,2f20.14)') 'CV33.01 ', input%dayy, input%hhrr, input%ta
     ! write(*,'(a,3f20.14)') 'CV33.02 ', input%ea, input%wnd, input%ppt(1)
     ! write(*,'(a,3f20.14)') 'CV33.03 ', input%co2air, input%press_mb, input%tsoil
     ! write(*,'(a,f20.14,i10,f20.14)') 'CV33.04 ', input%soilmoisture, input%flag, input%d13CO2
     ! write(*,'(a,f20.14,2i10)') 'CV33.05 ', input%d18CO2, time%year, time%daytime
     ! write(*,'(a,i10,f20.14,i10)') 'CV33.06 ', time%doy, time%local_time, time%days
     ! write(*,'(a,3f20.14)') 'CV33.07 ', met%T_Kelvin, met%rhova_g, met%rhova_kg
     ! write(*,'(a,3f23.14)') 'CV33.08 ', met%press_kpa, met%press_bars, met%press_Pa
     ! write(*,'(a,3f20.14)') 'CV33.09 ', met%relative_humidity, met%pstat273
     ! write(*,'(a,2f23.14)') 'CV33.10 ', srf_res%rcuticle(1), srf_res%rcuticle(ncl)
     ! write(*,'(a,3f20.14)') 'CV33.11 ', input%co2air, input%parin
     ! write(*,'(a,3f20.14)') 'CV33.13 ', input%wnd, met%ustar_filter
     ! write(*,'(a,3f20.14)') 'CV33.14 ', input%ta, met%air_density, met%air_density_mole
     ! write(*,'(a,3f20.14)') 'CV33.15 ', input%dppt(1), input%dppt(2), input%dppt(3)
     ! write(*,'(a,3f20.14)') 'CV33.16 ', input%dppt(4), input%dvapour(2), input%dppt(3)
     ! write(*,'(a,2f20.14)') 'CV33.17 ', wiso%dtheta(1,2), wiso%dtheta(1,nwiso-1)
     ! write(*,'(a,2f20.14)') 'CV33.18 ', input%ppt(2), input%ppt(nwiso-1)
     ! write(*,'(a,3f20.14)') 'CV33.19 ', input%lai_up, input%lai_down, input%lai
     ! write(*,'(a,3f20.14)') 'CV33.20 ', input%rglobal, input%parin, input%pardif
     ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
     ! update LAI in each time step
     if (extra_nate==1) then
        call lai_time()
     else
        ! start of new day
        ! if ((time%daytime+(nint(time%time_step,kind=i8)*100_i8/3600_i8)) & ! cf. mo_io_text.f90:read_text_in
        !      >= (int(input%dayy,kind=i8)*10000_i8+2400)) then
        ! Matthias: <= or <
        ! Yuan 2017: <=
        if ((time%daytime-(int(time%time_step,kind=i8)*100_i8/3600_i8)) &
             <= (int(input%dayy,kind=i8)*10000_i8)) call lai_time()
     endif
     ! write(*,'(a,1f20.14)') 'CV34.01 ', time%lai
     ! write(*,'(a,3f20.14)') 'CV34.02 ', solar%par_reflect, solar%par_trans, solar%par_soil_refl_dry
     ! write(*,'(a,3f20.14)') 'CV34.03 ', solar%par_absorbed, solar%nir_reflect, solar%nir_trans
     ! write(*,'(a,3f20.14)') 'CV34.04 ', solar%nir_soil_refl_dry, solar%nir_absorbed
     ! write(*,'(a,2f20.14)') 'CV34.05 ', prof%dLAIdz(1), prof%dLAIdz(ncl)
     ! write(*,'(a,2f20.14)') 'CV34.06 ', prof%dPAIdz(1), prof%dPAIdz(ncl)
     ! write(*,'(a,2f20.14)') 'CV34.07 ', solar%exxpdir(1), solar%exxpdir(ncl)

     ! Compute solar elevation angle
     call angle()
     ! write(*,'(a,3f20.14)') 'CV35 ', solar%beta_rad, solar%sine_beta, solar%beta_deg
     ! make sure PAR is zero at night. Some data have negative offset, which causes numerical problems
     ! Yuan 2017-07-01: ((solar%sine_beta <= isnight) .OR. (input%parin <=0))
     ! Matthias 2019: Not done in mo_radiation so inconsistent. Set parin >= 0 instead.
     if (solar%sine_beta <= isnight) input%parin = zero
     if (input%parin < zero) input%parin = zero
     ! Compute the fractions of beam and diffuse radiation from incoming measurements
     ! Set the radiation factor to the day before for night calculations. This way if the day
     !   was cloudy, so will IR calculations for the night reflect this.
     if (.not. (solar%sine_beta <= isnight)) then
        call diffuse_direct_radiation()
     else
        solar%ratrad      = solar%ratradnoon
        ! moved from read_data
        solar%par_beam    = zero
        solar%par_diffuse = zero
        solar%nir_beam    = zero
        solar%nir_diffuse = zero
     end if
     ! write(*,'(a,3f20.14)') 'CV36.01 ', solar%ratrad, output%c10, solar%ratradnoon
     ! write(*,'(a,3f20.14)') 'CV36.02 ', solar%par_beam, solar%par_diffuse
     ! write(*,'(a,3f20.14)') 'CV36.03 ', solar%nir_beam, solar%nir_diffuse
     ! computes leaf inclination angle distribution function, the mean direction cosine
     ! between the sun zenith angle and the angle normal to the mean leaf
     ! for CANOAK we use the leaf inclination angle data of Hutchison et al. 1986, J Ecology
     ! for other canopies we assume the leaf angle distribution is spherical
     if (.not. (solar%sine_beta <= isnight)) then
        call gfunc()
     else
        prof%Gfunc_solar(1:ncl) = 0.01_wp
     end if
     ! write(*,'(a,3f20.14)') 'CV37 ', prof%Gfunc_solar(1), prof%Gfunc_solar(ncl)
     ! Set soil reflectivity depending on soil moisture of first layer
     !   i.e. wet soil seems half as bright as dry soil
     !   after Wilson & Henderson-Sellers (1985)
     !   cited in Knorr PhD (1997, Table 2.4, p. 34)
     ztmp = min(max((soil%theta(1,1)-soil%watmin(1))/(soil%theta_s(1)-soil%watmin(1)), zero), one)
     solar%par_soil_refl = solar%par_soil_refl_dry * (one-half*ztmp)
     solar%nir_soil_refl = solar%nir_soil_refl_dry * (one-half*ztmp)
     ! write(*,'(a,3f20.14)') 'CV38.01 ', solar%par_soil_refl_dry, solar%nir_soil_refl_dry
     ! write(*,'(a,3f20.14)') 'CV38.02 ', solar%par_soil_refl, solar%nir_soil_refl
     ! write(*,'(a,3f20.14)') 'CV38.03 ', soil%theta(1,1), soil%watmin(1), soil%theta_s(1)
     ! Compute PAR profiles
     call par()
     ! write(*,'(a,3f20.14)') 'CV39.01 ', prof%sun_lai(1), prof%sun_lai(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.02 ', prof%shd_lai(1), prof%shd_lai(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.03 ', solar%prob_beam(1), solar%prob_beam(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.04 ', solar%prob_shd(1), solar%prob_shd(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.05 ', solar%par_down(1), solar%par_down(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.06 ', solar%par_up(1), solar%par_up(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.07 ', solar%beam_flux_par(1), solar%beam_flux_par(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.08 ', solar%quantum_shd(1), solar%quantum_shd(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.09 ', solar%quantum_sun(1), solar%quantum_sun(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.10 ', solar%par_shd(1), solar%par_shd(ncl)
     ! write(*,'(a,3f20.14)') 'CV39.11 ', solar%par_sun(1), solar%par_sun(ncl)
     call nir()
     ! write(*,'(a,3f20.14)') 'CV40.01 ', solar%nir_up(1), solar%nir_up(ncl)
     ! write(*,'(a,3f20.14)') 'CV40.02 ', solar%nir_dn(1), solar%nir_dn(ncl)
     ! write(*,'(a,3f20.14)') 'CV40.03 ', solar%nir_shd(1), solar%nir_shd(ncl)
     ! write(*,'(a,3f20.14)') 'CV40.04 ', solar%nir_sun(1), solar%nir_sun(ncl)
     ! write(*,'(a,3f20.14)') 'CV40.05 ', solar%beam_flux_nir(1), solar%beam_flux_nir(ncl)
     ! write(*,'(a,3f20.14)') 'CV40.06 ', solar%nir_total
     ! Interception reservoir
     call throughfall()
     ! write(*,'(a,3f20.14)') 'CV41.01 ', prof%throughfall(1,1), prof%throughfall(ncl,nwiso-1)
     ! write(*,'(a,3f20.14)') 'CV41.02 ', prof%cws(1,1), prof%cws(ncl,nwiso-1)
     ! write(*,'(a,3f20.14)') 'CV41.03 ', prof%wet_coef_filter(1), prof%wet_coef_filter(ncl)
     ! Update litter with rain
     if (soil%camillo==0 .and. soil%z_litter >= e6) call litter_rain()
     ! calculate soil heat capacity and conductivity based on soil moisture, texture and bulk density
     call set_soil_litter_capacity()
     fact%heatcoef = met%air_density * cp

     ! any recursive looping would occur here, below TL initialization
     time%count  = 0
     i_count     = 0
     etest_old1  = zero
     etest_old2  = -one
     itest_old1  = zero
     itest_old2  = -one
     fact%a_filt = 0.85_wp
     met%ustar   = met%ustar_filter
     if (soil%z_litter < e6) then
        soil%c_litterevap = zero
     else
        soil%c_litterevap = one
     end if
     ! write(*,'(a,3f23.14)') 'CV42.01 ', soil%cp_soil(1), soil%cp_soil(nsoil)
     ! write(*,'(a,3f20.14)') 'CV42.02 ', soil%k_conductivity_soil(1), soil%k_conductivity_soil(nsoil)
     ! write(*,'(a,3f20.14)') 'CV42.03 ', fact%heatcoef, met%ustar

     ! iteration looping for energy fluxes and scalar fields
     ! iterate until energy balance closure occurs or i_max iterations are reached
     do
        ! save input filtered variables for use in wiso routines
        prof%tair_filter_save(1:ntl) = prof%tair_filter(1:ntl)
        prof%rhov_air_filter_save(1:ntl) = prof%rhov_air_filter(1:ntl,1)
        prof%rhov_air_save(1:ntl) = prof%rhov_air(1:ntl,1)
        ! write(*,'(a,3f23.14)') 'CV43.01 ', prof%tair_filter_save(1), prof%tair_filter_save(ntl)
        ! write(*,'(a,3f23.14)') 'CV43.01 ', prof%rhov_air_filter_save(1), prof%rhov_air_filter_save(ntl)
        ! write(*,'(a,3f23.14)') 'CV43.01 ', prof%rhov_air_save(1), prof%rhov_air_save(ntl)
        !print *, prof%rhov_air_save(1:ntl)
        call irflux() ! first guess
        ! write(*,'(a,3f23.14)') 'CV44.01 ', solar%ir_up(1), solar%ir_up(ncl)
        ! write(*,'(a,3f23.14)') 'CV44.02 ', solar%ir_dn(1), solar%ir_dn(ncl)
        call friction_velocity()
        ! write(*,'(a,3f23.14)') 'CV45 ', met%zl, met%ustar
        ! compute net radiation balance on sunlit and shaded leaves
        call rnet()
        ! write(*,'(a,i5,2f20.14)') 'CV46.01 ', i_count, solar%par_sun(1), solar%nir_sun(1)
        ! write(*,'(a,2f20.14)') 'CV46.02 ', solar%par_sun(40), solar%nir_sun(40)
        ! write(*,'(a,2f20.14)') 'CV46.03 ', solar%par_shd(1), solar%nir_shd(1)
        ! write(*,'(a,2f20.14)') 'CV46.04 ', solar%par_shd(40), solar%nir_shd(40)
        ! write(*,'(a,2f20.14)') 'CV46.05 ', solar%ir_dn(1), solar%ir_up(1)
        ! write(*,'(a,2f20.14)') 'CV46.06 ', solar%ir_dn(40), solar%ir_up(40)
        ! write(*,'(a,2f20.14)') 'CV46.07 ', solar%rnet_sun(1), solar%rnet_shd(1)
        ! write(*,'(a,2f20.14)') 'CV46.08 ', solar%rnet_sun(40), solar%rnet_shd(40)

        ! Compute leaf energy balance, leaf temperature, photosynthesis and stomatal conductance
        call energy_and_carbon_fluxes()
        ! write(*,'(a,i5,2f20.14)') 'CV47.01 ', i_count, prof%dLEdz(1,1), prof%dLEdz(40,1)
        ! write(*,'(a,2f20.14)') 'CV47.02 ', prof%dHdz(1), prof%dHdz(40)
        ! write(*,'(a,2f20.14)') 'CV47.03 ', prof%dRNdz(1), prof%dRNdz(40)
        ! write(*,'(a,2f20.14)') 'CV47.04 ', prof%dLoutdz(1), prof%dLoutdz(40)
        ! write(*,'(a,2f20.14)') 'CV47.00 ', prof%sun_tleaf(1), prof%sun_tleaf(40)
        ! write(*,'(a,2f20.14)') 'CV47.05 ', prof%shd_tleaf(1), prof%shd_tleaf(40)
        ! Soil energy balance
        call soil_energy_balance()
        ! write(*,'(a,3f20.14)') 'CV47.06 ', flux%soilevap(1), flux%litterevap(1), flux%s_evap(1)
        ! write(*,'(a,3f20.14)') 'CV47.07 ', soil%tsrf, soil%T_l, soil%litterevap
        ! write(*,'(a,3f20.14)') 'CV47.08 ', output%c1, output%c2, output%c3
        ! write(*,'(a,2f20.14)') 'CV47.09 ', output%c4, output%c7
        ! write(*,'(a,3f20.14)') 'CV47.10 ', soil%theta(1:3,1)

        ! --- Water ---
        if (iswitch%wiso==1) then
           ! calculates latent heat in accordance with leaf water isotope calculations
           call le_wiso()
           ! re-define/re-calc LEstoma in accordance with leaf water isotope calculations
           prof%sun_LEstoma_save(1:ncl) = prof%sun_LEstoma(1:ncl,1)
           prof%shd_LEstoma_save(1:ncl) = prof%shd_LEstoma(1:ncl,1)
           prof%sun_LEstoma(1:ncl,1)    = prof%sun_LEstoma_new(1:ncl) &
                * lambda(prof%tair_filter(1:ncl)+TN0)
!                print *, prof%sun_LEstoma_new(1:ncl)
 !               print *, prof%tair_filter(1:ncl)
           prof%shd_LEstoma(1:ncl,1)    = prof%shd_LEstoma_new(1:ncl) &
                * lambda(prof%tair_filter(1:ncl)+TN0)
           prof%dLEdz(1:ncl,1)          = &
                (solar%prob_beam(1:ncl) * (prof%sun_LEstoma(1:ncl,1)+prof%sun_LEwet(1:ncl,1)) &
                + solar%prob_shd(1:ncl) * (prof%shd_LEstoma(1:ncl,1)+prof%shd_LEwet(1:ncl,1))) &
                *  prof%dLAIdz(1:ncl)
           ! write(*,'(a,2f20.14)') 'CV48.01 ', prof%sun_LEstoma(1,1), prof%sun_LEstoma(ncl,1)
           ! write(*,'(a,2f20.14)') 'CV48.02 ', prof%dLEdz(1,1), prof%dLEdz(ncl,1)
        end if

        ! compute canopy transpiration and evaporation
        flux%c_evaporation(1) = sum((prof%sun_LEwet(1:ncl,1)*solar%prob_beam(1:ncl) &
             + prof%shd_LEwet(1:ncl,1)*solar%prob_shd(1:ncl)) &          !convert from W m-2
             * prof%dLAIdz(1:ncl) / LAMBDA(prof%tair_filter(1:ncl)+TN0)) !        to kg H2O m-2 s-1
        flux%c_transpiration(1) = sum((prof%sun_LEstoma(1:ncl,1)*solar%prob_beam(1:ncl) &
             + prof%shd_LEstoma(1:ncl,1)*solar%prob_shd(1:ncl)) &
             * prof%dLAIdz(1:ncl) / LAMBDA(prof%tair_filter(1:ncl)+TN0))
        temp3 = sum((solar%rnet_sun(1:ncl)*solar%prob_beam(1:ncl) &
             + solar%rnet_shd(1:ncl)*solar%prob_shd(1:ncl)) * prof%dLAIdz(1:ncl))

        ! total flux
        flux%c_evapotranspiration(1) = flux%c_evaporation(1) + flux%c_transpiration(1)
        flux%evapotranspiration(1)   = flux%c_evapotranspiration(1) + flux%s_evap(1)
        ! write(*,'(a,2f20.14)') 'CV48.03 ', flux%c_evaporation(1), flux%c_evapotranspiration(1)

        ! Compute air temperature profiles from dispersion matrix
        ! Adjust dHdz(1) for soil heat flux
        ! sign convention used: fluxes from surface are positive
        ! those toward the surface are negative
        ! dHdz is net sensible heat exchange per unit leaf area
        ! Using filtered temperatures to minimize the system from being mathematically unstable.

        ! filter temperatures with each interation to minimize numerical instability

        ! Original filter
        !   if (i_count < 10) then
        !     fact%a_filt = 0.85_wp
        !   else
        !     fact%a_filt = half
        !   end if

        ! Matthias, constant filter
        !   fact%a_filt = 0.85_wp
        ! Matthias, force conversion with narrowing filter
        !   fact%a_filt = 0.99_wp - 0.98_wp*real(i_count,wp)/real(i_max-1,wp)
        ! Matthias, narrowing filter to 0.5 at i_max/2
        fact%a_filt = 0.85_wp - 0.7_wp*real(i_count,wp)/real(i_max-1,wp)
!       fact%a_filt = fact%a_filt*0.7_wp
        ! write(*,'(a,2f20.14)') 'CV48.04 ', fact%a_filt
        ! Matthias, a_filt=0 from 10 steps before max iteration
        ! the noise one sees should be filtered out as well
        !   fact%a_filt = max(0.99_wp - one*real(i_count,wp)/real(i_max-10,wp), zero)

        ! conc, for temperature profiles using source/sinks
        ! inputs are source/sink(), scalar_profile(), ref_val, boundary_flux, unit conversions
! print *, "1.5\n"
!         print *, prof%tair
!        call message('conc Tair')
        call conc(prof%dHdz, prof%tair, input%ta, soil%heat, fact%heatcoef)
        ! write(*,'(a,3f20.14)') 'CV49.01 ', prof%dHdz(1), prof%dHdz(ncl)
        ! write(*,'(a,3f20.14)') 'CV49.02 ', input%ta, soil%heat, fact%heatcoef
        ! write(*,'(a,3f20.14)') 'CV49.03 ', prof%tair(1), prof%tair(ntl-1)

        ! filter temperatures to remove numerical instabilities for each iteration
        where (prof%tair(1:ntl) < -30._wp .or. prof%tair(1:ntl) > 60._wp) prof%tair(1:ntl) = input%ta
        ! write(*,'(a,3f20.14)') 'CV49.04 ', prof%tair(1), prof%tair(ntl-1)

        ! Compute vapor density profiles from Dispersion matrix
        ! sign convention used: fluxes from surface are positive
        ! those toward the surface are negative
        ! dLEdZ is net latent heat exchange per unit leaf area
        ! turbulent transport of H2O
        tmpncl(1:ncl) = prof%dLEdz(1:ncl,1) / lambda(prof%tair_filter_save(1:ncl)+TN0)
!print *,  prof%dLEdz(1:ncl,1)
!print *, prof%tair_filter_save
!print *, tmpncl(1:ncl)
        call conc(tmpncl, tmpntl, met%rhova_kg, flux%s_evap(1), one) ! negetive flux%s_evap! Yuan 2018.05.29
        prof%rhov_air(1:ntl,1) = tmpntl(1:ntl)
        ! write(*,'(a,3f20.14)') 'CV49.05 ', tmpncl(1), tmpncl(ncl)
        ! write(*,'(a,3f20.14)') 'CV49.06 ', met%rhova_kg, flux%s_evap(1), one
        ! write(*,'(a,3f20.14)') 'CV49.07 ', prof%rhov_air(1,1), prof%rhov_air(ntl-1,1)

        ! filter humidity computations
        where ((prof%rhov_air(1:ntl,1) < zero) .or. (prof%rhov_air(1:ntl,1) > 0.03_wp)) &
             prof%rhov_air(1:ntl,1) = met%rhova_kg
        ! write(*,'(a,3f20.14)') 'CV49.08 ', prof%rhov_air(1,1), prof%rhov_air(ntl-1,1)
        ! ea
        ztmp = one / 2.165_wp
        tmpntl(1:ntl) = ztmp * prof%rhov_air(1:ntl,1) * (prof%tair(1:ntl)+TN0)
        ! vapor pressure deficit in (kPa)
#ifdef DEBUG
        do j=1, ntl
           prof%vpd_air(j) = es(prof%tair(j)+TN0) - tmpntl(j)
        end do
#else
        prof%vpd_air(1:ntl) = es(prof%tair(1:ntl)+TN0) - tmpntl(1:ntl)
#endif
!print *, input%ea
!print *, prof%rhov_air
!print *, prof%vpd_air
!        do j=1, ntl
!           prof%u(j) = uz(j*ht)
!        end do
        ! write(*,'(a,3f20.14)') 'CV49.09 ', prof%vpd_air(1), prof%vpd_air(ntl-1)
        ! Implicit water isotopes diagnostics
        if (iswitch%wiso==1 .and. wiso%implicit==1) then
           ! isotope soil water flux
           call soil_flux_wiso()
           ! write(*,'(a,3f20.14)') 'CV52.01 ', flux%soilevap(1), flux%soilevap(nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.02 ', flux%litterevap(1), flux%litterevap(nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.03 ', flux%s_evap(1), flux%s_evap(nwiso-1)
           ! leaf water enrichment
           call leaf_wiso()
           ! write(*,'(a,3f20.14)') 'CV52.04 ', prof%sun_wi(1), prof%sun_wi(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.05 ', prof%shd_wi(1), prof%shd_wi(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.06 ', prof%wa(1), prof%wa(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.07 ', prof%sun_h(1), prof%sun_h(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.08 ', prof%rs_fact(1), prof%rs_fact(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.09 ', prof%sun_gross(1), prof%sun_gross(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.10 ', prof%shd_gross(1), prof%shd_gross(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.11 ', prof%sun_LEstoma_new(1), prof%sun_LEstoma_new(ncl)
           ! write(*,'(a,3f20.14)') 'CV52.12 ', prof%shd_LEstoma_new(1), prof%shd_LEstoma_new(ncl)
           ! isotope canopy transpiration and evaporation
           call canopy_flux_wiso()
           ! write(*,'(a,3f20.14)') 'CV52.13 ', prof%sun_alpha_k(1,1), prof%sun_alpha_k(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.14 ', prof%shd_alpha_k(1,1), prof%shd_alpha_k(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.15 ', prof%sun_alpha_equ(1,1), prof%sun_alpha_equ(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.16 ', prof%shd_alpha_equ(1,1), prof%shd_alpha_equ(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.17 ', prof%sun_peclet(1,1), prof%sun_peclet(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.18 ', prof%shd_peclet(1,1), prof%shd_peclet(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.19 ', prof%sun_fem(1,1), prof%sun_fem(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.20 ', prof%shd_fem(1,1), prof%shd_fem(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.21 ', prof%sun_craig(1,1), prof%sun_craig(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.22 ', prof%shd_craig(1,1), prof%shd_craig(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.23 ', prof%sun_leafwater_e(1,1), prof%sun_leafwater_e(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.24 ', prof%shd_leafwater_e(1,1), prof%shd_leafwater_e(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.25 ', prof%sun_leafwater(1,1), prof%sun_leafwater(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.26 ', prof%shd_leafwater(1,1), prof%shd_leafwater(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.27 ', prof%sun_trans_rtrans(1,1), prof%sun_trans_rtrans(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.28 ', prof%shd_trans_rtrans(1,1), prof%shd_trans_rtrans(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.29 ', prof%sun_rtrans(1,1), prof%sun_rtrans(ncl,nwiso-1)
           ! write(*,'(a,3f20.14)') 'CV52.30 ', prof%shd_rtrans(1,1), prof%shd_rtrans(ncl,nwiso-1)
           ! turbulent transport of H2O18, DHO and H216O
           do mc=2, nwiso
              ! Vapour input in layers
              tmpncl(1:ncl) = prof%dLEdz(1:ncl,mc) / lambda(prof%tair_filter_save(1:ncl)+TN0)
              ! Vapour input at top
              ztmp = alpha_equ_h2o(input%ta+TN0, mc) * invdelta1000_h2o(input%dppt(mc), mc) * met%rhova_kg ! in equi with rain
              call conc(tmpncl, tmpntl, ztmp, flux%s_evap(mc), one)
              prof%rhov_air(1:ntl,mc) = tmpntl(1:ntl)
!              if (prof%rhov_air(1,mc) < zero ) then
!                print *, prof%rhov_air(1,mc)
!              end if
           end do
           ! filter humidity computations
           do mc=2, nwiso
!            print *, prof%rhov_air(1:ntl-1,1)
! .or. prof%rhov_air(1:ntl-1,mc) < zero
              where ((prof%rhov_air(1:ntl-1,1) == met%rhova_kg)) &
                   prof%rhov_air(1:ntl-1,mc) = prof%rhov_air(1:ntl-1,1) * prof%rvapour(1:ntl-1,mc)

           end do
        end if ! implicit water isotope diagnostics

        ! sign convention used: photosynthetic uptake is positive
        ! respiration is negative
        ! prof%dPsdz is net photosynthesis per unit ground area,
        ! the profile was converted to units of mg m-2 s-1 to be
        ! consistent with inputs to CONC
        ! change sign of dPsdz
        prof%source_co2(1:ncl) = -prof%dPsdz(1:ncl)
        ! write(*,'(a,3f20.14)') 'CV50.01 ', prof%source_co2(1), prof%source_co2(ncl)

        ! compute bole respiration
        ! It is added to the prof%source_CO2 array.
        call bole_respiration()
        ! write(*,'(a,3f20.14)') 'CV50.02 ', bole%calc, bole%factor
        ! write(*,'(a,3f20.14)') 'CV50.03 ', bole%respiration_mole, bole%respiration_mg
        ! write(*,'(a,3f20.14)') 'CV50.04 ', bole%layer(1), bole%layer(ncl)
        ! write(*,'(a,3f20.14)') 'CV50.05 ', prof%source_co2(1), prof%source_co2(ncl)
        ! compute soil respiration
        call soil_respiration()
        ! write(*,'(a,3f20.14)') 'CV50.06 ', soil%respiration_mole, soil%respiration_mg
        ! write(*,'(a,3f20.14)') 'CV50.07 ', soil%respiration_auto, soil%respiration_hetero

        ! to convert umol m-3 to umol/mol we have to consider
        ! Pc/Pa (CO2 density/ air density) = (CO2)ppm = rhoc ma/ rhoa mc
        fact%co2 = (mass_air/mass_CO2)*met%air_density_mole
!        call conc(prof%source_co2, prof%co2_air, input%co2air, soil%respiration_mole, fact%co2)
     call conc_seperate(prof%source_CO2, prof%co2_air, input%co2air, soil%respiration_mole, prof%CO2_soil, prof%CO2_disp, fact%co2)
!print *,  prof%co2_air
!print *,  soil%respiration_mole
!print *,  fact%co2
        if (iswitch%oxygen==1) then
            call OXYFLUX ()
        end if
!        call conc(prof%source_co2, prof%co2_air, input%co2air, soil%respiration_mole, fact%co2)
        ! write(*,'(a,3f20.14)') 'CV50.08 ', prof%source_co2(1), prof%source_co2(ncl)
        ! write(*,'(a,3f20.14)') 'CV50.09 ', input%co2air, soil%respiration_mole, fact%co2
        ! write(*,'(a,3f20.14)') 'CV50.10 ', prof%co2_air(1), prof%co2_air(ntl-1)

        ! Integrate source-sink strengths to estimate canopy flux
        sumh       = sum(prof%dHdz(1:ncl))       ! sensible heat
        sumle      = sum(prof%dLEdz(1:ncl,1))    ! latent heat
!        print *, "canopy H and LE:"
!        print *, sumh, sumle
        sumrn      = sum(prof%dRNdz(1:ncl))      ! net radiation
        sumlout    = sum(prof%dLoutdz(1:ncl))    ! radiation
        can_ps_mol = sum(prof%dPsdz(1:ncl))      ! canopy photosynthesis
        can_ps_o   = sum(prof%dPsdz_O2(1:ncl))   ! net O2 emission
        can_gpp    = sum(prof%dGPPdz(1:ncl))     ! canopy GPP = can_ps_mol + dark respiration
        can_gpp_o  = sum(prof%gpp_O2(1:ncl))     ! O2 emmision via gpp
        canresp    = sum(prof%dRESPdz(1:ncl))    ! canopy respiration
        canresp_o  = sum(prof%dRESPdz_O2(1:ncl)) ! O2 via canopy respiration
        sumksi     = sum(prof%dStomCondz(1:ncl)) ! canopy stomatal conductance
        sumlai     = sum(prof%dLAIdz(1:ncl))     ! leaf area
        sumpai     = sum(prof%dPAIdz(1:ncl))     ! use plant area Yuan 2018.03.04
        sumNsupply = sum(prof%dNsupplydz(1:ncl)) ! total N supply per hour
        sumNdemand = sum(prof%dNdemanddz(1:ncl))    ! total N demand per hour
 !       print*, 'T: ', i_count, prof%shd_tleaf(1), prof%shd_tleaf(40)
! print *, i_count
        tleaf_mean = sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
             + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))    ! mean leaf temperature
        ztmp = one / real(ncl,wp)
!        prof%tleaf(1:ncl) = (sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
!             + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))) * ztmp    ! Tleaf per layer (sun and shade)
        prof%tleaf(1:ncl) = ((prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
             + (prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl)))    ! Tleaf per layer (sun and shade)

        ! need to weight by sun and shaded leaf areas then divide by LAI
        tavg_sun   = sum(prof%sun_tleaf(1:ncl)*prof%dPAIdz(1:ncl)) ! avg sunlit temperature USE PAI YUAN 2018.03.04
        tavg_shd   = sum(prof%shd_tleaf(1:ncl)*prof%dPAIdz(1:ncl)) ! avg shaded temperature
        ! write(*,'(a,3f20.14)') 'CV50.11 ', sumh, sumle, sumrn
        ! write(*,'(a,3f20.14)') 'CV50.12 ', sumlout, can_ps_mol, can_gpp
        ! write(*,'(a,3f20.14)') 'CV50.13 ', canresp, sumksi, sumlai
        ! write(*,'(a,i5,2f20.14)') 'T: ', i_count, prof%shd_tleaf(1), prof%shd_tleaf(40)
        !tleaf_mean = sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
        !     + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))    ! mean leaf temperature
        ! ztmp = one / real(ncl,wp)
        ! prof%tleaf(1:ncl) = (sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
        !      + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))) * ztmp    ! Tleaf per layer (sun and shade)
        !prof%tleaf(1:ncl) = prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl) &
        !     + prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl) ! Tleaf per layer (sun and shade)
        ! write(*,'(a,3f20.14)') 'CV50.14 ', prof%tleaf(1), prof%tleaf(ncl)
        ! need to weight by sun and shaded leaf areas then divide by LAI
        !tavg_sun   = sum(prof%sun_tleaf(1:ncl)*prof%dLAIdz(1:ncl)) ! avg sunlit temperature
        !tavg_shd   = sum(prof%shd_tleaf(1:ncl)*prof%dLAIdz(1:ncl)) ! avg shaded temperature
        ! write(*,'(a,3f20.14)') 'CV50.15 ', tavg_sun, tavg_shd

        ebalance      = sumrn - sumle - sumh
!        print *, 'energy balance:'
!        print *, ebalance
        flux%photosyn = can_ps_mol
        ! write(*,'(a,3f20.14)') 'CV50.16 ', ebalance, flux%photosyn
        ! calculate gpp : can_ps_mol = photosynthesis - photorespiration
        !                              - dark respiration or can_ps_mol
        !                            = GPP - leaf respiration
        ! can_gpp = can_ps_mol + canresp

        ! mean canopy leaf temperature
        tleaf_mean = tleaf_mean / real(ncl,wp)
        ! leaf area weighted temperatures
        tavg_sun = tavg_sun / sumpai
        tavg_shd = tavg_shd / sumpai
        !tavg_sun = tavg_sun / sumlai
        !tavg_shd = tavg_shd / sumlai
        ! write(*,'(a,3f20.14)') 'CV50.017 ', tleaf_mean, tavg_sun, tavg_shd
        ! Energy exchanges at the soil
        rnet_soil = soil%rnet - soil%lout
        sbalance  = rnet_soil - soil%soilevap - soil%litterevap - soil%heat - soil%gsoil
        ! write(*,'(a,3f20.14)') 'CV50.18 ', rnet_soil, sbalance
        ! canopy scale flux densities, vegetation plus soil
        sumh     = sumh + soil%heat
        sumle    = sumle + soil%evap
!        print *, "soil heat:"
!        print *, soil%heat, soil%evap, rnet_soil
        sumrn    = sumrn + rnet_soil
        sumlout  = sumlout + soil%lout
        temp3    = temp3 + soil%rnet
        tbalance = sbalance + ebalance
        met%H    = sumh
        ! write(*,'(a,3f20.14)') 'CV50.19 ', sumh, sumle, sumrn
        ! write(*,'(a,3f20.14)') 'CV50.20 ', sumlout, temp3, tbalance
        ! write(*,'(a,3f20.14)') 'CV50.21 ', met%H

        ! filter iterative variables
        ! soil temp and fluxes
        soil%tsrf_filter              = fact%a_filt * soil%tsrf &
             + (one-fact%a_filt) * soil%tsrf_filter
        soil%soilevap_filter          = fact%a_filt * soil%soilevap &
             + (one-fact%a_filt) * soil%soilevap_filter
        soil%litterevap_filter        = fact%a_filt * soil%litterevap &
             + (one-fact%a_filt) * soil%litterevap_filter
        soil%T_l_filter               = fact%a_filt * soil%T_l &
             + (one-fact%a_filt) * soil%T_l_filter
        soil%T_soil_filter(0:nsoil+1) = fact%a_filt * soil%T_soil(0:nsoil+1) &
             + (one-fact%a_filt) * soil%T_soil_filter(0:nsoil+1)
        ! write(*,'(a,3f20.14)') 'CV50.22 ', soil%tsrf_filter, soil%soilevap_filter, soil%litterevap_filter
        ! write(*,'(a,3f20.14)') 'CV50.23 ', soil%T_l_filter
        ! write(*,'(a,3f20.14)') 'CV50.231 ', soil%T_soil_filter(0), soil%T_soil_filter(nsoil+1)
        ! met variables
        met%H_filter     = fact%a_filt * met%H + (one-fact%a_filt) * met%H_filter
        met%ustar_filter = fact%a_filt * met%ustar + (one-fact%a_filt) * met%ustar_filter
        met%K_filter = fact%a_filt * met%K + (one-fact%a_filt) * met%K_filter
        ! write(*,'(a,3f20.14)') 'CV50.24 ', met%H_filter, met%ustar_filter
        ! air variables
!        print *, "3\n"
!        print *, prof%tair
        prof%tair_filter(1:ntl)       = fact%a_filt * prof%tair(1:ntl) &
             + (one-fact%a_filt) * prof%tair_filter(1:ntl)
!             print *, prof%rhov_air_filter(1,1)
        prof%rhov_air_filter(1:ntl,1) = fact%a_filt * prof%rhov_air(1:ntl,1) &
             + (one-fact%a_filt) * prof%rhov_air_filter(1:ntl,1)
 !            print *, prof%rhov_air_filter(1,1)
             if (prof%rhov_air_filter(1,1)==zero) then
 !            print *, prof%rhov_air
 !            print *, fact%a_filt
             end if
        prof%co2_air_filter(1:ntl)    = fact%a_filt * prof%co2_air(1:ntl) &
             + (one-fact%a_filt) * prof%co2_air_filter(1:ntl)
        prof%O2_air_filter(1:ntl)    = fact%a_filt * prof%O2_air(1:ntl) &
             + (one-fact%a_filt) * prof%O2_air_filter(1:ntl) ! Yuan 2018.07.02
        ! write(*,'(a,3f20.14)') 'CV50.25 ', prof%tair_filter(1), prof%tair_filter(ntl-1)
        ! write(*,'(a,3f20.14)') 'CV50.26 ', prof%rhov_air_filter(1,1), prof%rhov_air_filter(ntl-1,1)
        ! write(*,'(a,3f20.14)') 'CV50.27 ', prof%co2_air_filter(1), prof%co2_air_filter(ntl-1)
        ! leaf variables
        prof%sun_tleaf_filter(1:ncl) = fact%a_filt * prof%sun_tleaf(1:ncl) &
             + (one-fact%a_filt) * prof%sun_tleaf_filter(1:ncl)
        prof%shd_tleaf_filter(1:ncl) = fact%a_filt * prof%shd_tleaf(1:ncl) &
             + (one-fact%a_filt) * prof%shd_tleaf_filter(1:ncl)
        prof%wet_coef_filter(1:ncl)  = fact%a_filt * prof%wet_coef(1:ncl) &
             + (one-fact%a_filt) * prof%wet_coef_filter(1:ncl)
        ! write(*,'(a,3f20.14)') 'CV50.28 ', prof%sun_tleaf_filter(1), prof%sun_tleaf_filter(ncl)
        ! write(*,'(a,3f20.14)') 'CV50.29 ', prof%shd_tleaf_filter(1), prof%shd_tleaf_filter(ncl)
        ! write(*,'(a,3f20.14)') 'CV50.30 ', prof%wet_coef_filter(1), prof%wet_coef_filter(ncl)

        ! Implicit water isotopes diagnostics
        if (iswitch%wiso==1 .and. wiso%implicit==1) then
           prof%rhov_air_filter(1:ntl,2:nwiso) = fact%a_filt*prof%rhov_air(1:ntl,2:nwiso) &
                + (one-fact%a_filt)*prof%rhov_air_filter(1:ntl,2:nwiso)
           ! isotopes in vapour
           prof%rvapour(1:ntl,2:nwiso) = isorat(prof%rhov_air_filter(1:ntl,2:nwiso), &
                prof%rhov_air_filter(1:ntl,1), wiso%lost(2:nwiso), wiso%lost(1))
 !               if ( abs(prof%rvapour(1,2)) < 0.0016_wp) then
 !               print *, prof%rvapour(1,1:nwiso)
 !               print *, prof%rhov_air_filter(1,1:nwiso)
 !               print *, prof%rhov_air(1,1:nwiso)
 !               print *, wiso%lost(1:nwiso)
 !               end if

#ifdef DEBUG
           if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
              call message('CANVEG', 'Lost 01 @, count ', num2str(time%daytime), ',', num2str(i_count))
              soil%lost0 = 1
           end if
#endif
        end if

        ! !!
        ! Net radiation balance at the top of the canopy
        ! PAR + NIR + IR
        netrad = (solar%beam_flux_par(ncl+1) + solar%par_down(ncl+1) - solar%par_up(ncl+1)) / 4.6_wp &
             + solar%beam_flux_nir(ncl+1) + solar%nir_dn(ncl+1) - solar%nir_up(ncl+1) &
             + solar%ir_dn(ncl+1) - solar%ir_up(ncl+1)
        ! write(*,'(a,3f20.14)') 'CV50.31 ', netrad
!print *, "direct par:   ", solar%beam_flux_par(ncl+1)/4.6_wp
!print *, "downward par:   ", solar%par_down(ncl+1)/4.6_wp
!print *, "upward par:   ", solar%par_down(ncl+1)/4.6_wp
        ! test for convergence between the sum of the net radiation flux profile and the
        ! net flux exiting the canopy
        ! etest = abs((sumrn-netrad)/sumrn)
        etest       = (sumrn-netrad) / sumrn
        etest_diff1 = etest_old1 - etest
        etest_diff2 = etest_old2 - etest_old1
        etest_old2  = etest_old1
        etest_old1  = etest
!        print *, "sumrn, netrd"
!        print *, sumrn,netrad
        itest = zero
        if (iswitch%wiso==1 .and. wiso%implicit==1) &
             itest = delta1000_h2o(flux%evapotranspiration(2), flux%evapotranspiration(1), 2)
        itest_diff1 = itest_old1 - itest
        itest_diff2 = itest_old2 - itest_old1
        itest_old2  = itest_old1
        itest_old1  = itest
        ! write(*,'(a,3f20.14)') 'CV50.32 ', etest, etest_diff1, etest_diff2
        ! write(*,'(a,3f20.14)') 'CV50.33 ', etest_old2, etest_old1
        ! write(*,'(a,3f20.14)') 'CV50.34 ', itest, itest_diff1, itest_diff2
        ! write(*,'(a,3f20.14)') 'CV50.35 ', itest_old2, itest_old1

        i_count    = i_count + 1
        time%count = i_count
!!        print *,time%count
        ! write(*,'(a,3i10)') 'CV50.36 ', time%count
!write (*,*) i_count
        ! check for convergence
        ! if (etest <= 0.005_wp .or. i_count >= i_max) exit
        !MC if ((abs(etest_diff1) <= 0.001_wp .and. abs(etest_diff2) <= 0.001_wp &
        !MC     .and. abs(itest_diff1) <= 0.01_wp .and. abs(itest_diff2) <= 0.01_wp) &
        !MC     .or. i_count >= i_max) exit
        if ((abs(etest_diff1) <= 0.005_wp .and. abs(etest_diff2) <= 0.005_wp &
             .and. abs(itest_diff1) <= 0.03_wp .and. abs(itest_diff2) <= 0.03_wp) &
             .or. i_count >= i_max) exit

        ! if ((abs(etest_diff1) <= 1e-15_wp .and. abs(etest_diff2) <= 1e-15_wp &
        !      .and. abs(itest_diff1) <= 1e-15_wp .and. abs(itest_diff2) <= 1e-15_wp) &
        !      .or. i_count >= i_max) exit
     end do ! until energy balance closure in this time step (or i_max)
!     print *, etest, etest_diff1,etest_diff2,itest_diff1,itest_diff2
     totalcount = totalcount + 1
     if (solar%sine_beta <= isnight) then
        ntotalcount = ntotalcount + 1
     else
        dtotalcount = dtotalcount + 1
     end if
     if (i_count == i_max) then
        count40 = count40 + 1
        if (solar%sine_beta <= isnight) then
           ncount40 = ncount40 + 1
        else
           dcount40 = dcount40 + 1
        end if
     end if
     if (print_balance==1) then
        message_text = '    Day '//trim(num2str(time%days,'(I5.5)'))//' Time ' &
             //trim(num2str(time%local_time,'(F4.1)'))//' i_c ' &
             //trim(num2str(i_count,'(I3.3)'))//': diff '//trim(num2str(abs(etest)*100._wp,'(F6.2)')) &
             //'% net '//trim(num2str(netrad,'(F7.2)'))//' sum '//trim(num2str(sumrn,'(F7.2)')) &
             //', soil '//trim(num2str(sbalance,'(F7.2)'))//' leaf ' &
             //trim(num2str(ebalance,'(F7.2)'))//' total '//trim(num2str(tbalance,'(F7.2)'))
        call message(trim(message_text))
     endif

     ! Explicit water isotope diagnostics
     if (iswitch%wiso == 1 .and. wiso%implicit == 0) then
        ! isotopes in vapour
        ! like that, we do not have to transport isotopes all the time
        !   but take the isotopes from one time step before
!        print *, prof%rhov_air_filter(1,1)
        ! but take the isotopes from one time step before
        prof%rhov_air_filter(1:ntl,2:nwiso) = spread(prof%rhov_air_filter_save(1:ntl),2,nwiso) * prof%rvapour(1:ntl,2:nwiso)
if (prof%rhov_air_filter(1,1)==zero) then
!    print *, prof%rhov_air_filter_save(1)
end if

        ! isotope soil water flux
!        print *, wiso%lost(1:nwiso)
        call soil_flux_wiso()
!                print *, wiso%lost(1:nwiso)
        ! write(*,'(a,3f20.14)') 'CV51.01 ', flux%soilevap(1), flux%soilevap(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.02 ', flux%litterevap(1), flux%litterevap(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.03 ', flux%s_evap(1), flux%s_evap(nwiso-1)
        ! leaf water enrichment
        call leaf_wiso()
        ! write(*,'(a,3f20.14)') 'CV51.04 ', prof%sun_wi(1), prof%sun_wi(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.05 ', prof%shd_wi(1), prof%shd_wi(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.06 ', prof%wa(1), prof%wa(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.07 ', prof%sun_h(1), prof%sun_h(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.08 ', prof%rs_fact(1), prof%rs_fact(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.09 ', prof%sun_gross(1), prof%sun_gross(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.10 ', prof%shd_gross(1), prof%shd_gross(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.11 ', prof%sun_LEstoma_new(1), prof%sun_LEstoma_new(ncl)
        ! write(*,'(a,3f20.14)') 'CV51.12 ', prof%shd_LEstoma_new(1), prof%shd_LEstoma_new(ncl)
        ! isotope canopy transpiration and evaporation
        call canopy_flux_wiso()
        ! write(*,'(a,3f20.14)') 'CV51.13 ', prof%sun_alpha_k(1,1), prof%sun_alpha_k(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.14 ', prof%shd_alpha_k(1,1), prof%shd_alpha_k(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.15 ', prof%sun_alpha_equ(1,1), prof%sun_alpha_equ(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.16 ', prof%shd_alpha_equ(1,1), prof%shd_alpha_equ(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.17 ', prof%sun_peclet(1,1), prof%sun_peclet(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.18 ', prof%shd_peclet(1,1), prof%shd_peclet(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.19 ', prof%sun_fem(1,1), prof%sun_fem(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.20 ', prof%shd_fem(1,1), prof%shd_fem(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.21 ', prof%sun_craig(1,1), prof%sun_craig(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.22 ', prof%shd_craig(1,1), prof%shd_craig(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.23 ', prof%sun_leafwater_e(1,1), prof%sun_leafwater_e(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.24 ', prof%shd_leafwater_e(1,1), prof%shd_leafwater_e(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.25 ', prof%sun_leafwater(1,1), prof%sun_leafwater(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.26 ', prof%shd_leafwater(1,1), prof%shd_leafwater(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.27 ', prof%sun_trans_rtrans(1,1), prof%sun_trans_rtrans(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.28 ', prof%shd_trans_rtrans(1,1), prof%shd_trans_rtrans(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.29 ', prof%sun_rtrans(1,1), prof%sun_rtrans(ncl,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV51.30 ', prof%shd_rtrans(1,1), prof%shd_rtrans(ncl,nwiso-1)

        ! turbulent transport of H2O18, DHO and H216O
        do mc=2, nwiso
           ! Vapour input in layers
           tmpncl(1:ncl) = prof%dLEdz(1:ncl,mc) / lambda(prof%tair_filter_save(1:ncl)+TN0)
           ! Vapour input at top
           ztmp = alpha_equ_h2o(input%ta+TN0, mc) * invdelta1000_h2o(input%dppt(mc), mc) * met%rhova_kg ! in equi with rain
           call conc(tmpncl, tmpntl, ztmp, flux%s_evap(mc), one)
           prof%rhov_air(1:ntl,mc) = tmpntl(1:ntl)
!           print *, ztmp, flux%s_evap(mc)
        end do
        ! TODO: check order of statements
        ! filter humidity computations
        do mc=2, nwiso
            ! Yuan added filter prof%rhov_air(1:ntl-1,mc) < zero
           where ((prof%rhov_air(1:ntl-1,1) == met%rhova_kg)) &
                prof%rhov_air(1:ntl-1,mc) = prof%rhov_air(1:ntl-1,1) * prof%rvapour(1:ntl-1,mc)
        end do
        ! isotopes in vapour
        prof%rvapour(1:ntl,2:nwiso) = isorat(prof%rhov_air(1:ntl,2:nwiso), &
             prof%rhov_air(1:ntl,1), wiso%lost(2:nwiso), wiso%lost(1))
 !            if (wiso%lost(1)>zero) then
 !              print *, wiso%lost(1:nwiso)
 !              print *, prof%rhov_air(1:ntl,2:nwiso)
 !            end if

        prof%dvapour(1:ntl,2) = delta1000_h2o(prof%rhov_air(1:ntl,2),prof%rhov_air(1:ntl,1),2)
        prof%dvapour(1:ntl,3) = delta1000_h2o(prof%rhov_air(1:ntl,3),prof%rhov_air(1:ntl,1),3)
        ! range d18O
        where (prof%dvapour(1:ntl,2)<-30)
            prof%dvapour(1:ntl,2) = input%dvapour
        end where



#ifdef DEBUG
        if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
           call message('CANVEG', 'Lost 02 @  ', num2str(time%daytime))
           soil%lost0 = 1
        end if
#endif
     end if ! explicit water isotope diagnostics

     ! subtract evaporation from wet leaf surfaces from canopy water storage
     do j=1, ncl
        if (prof%cws(j,1) > zero) then
           if (iswitch%wiso == 1) then
              rcws(2:nwiso) = isorat(prof%cws(j,2:nwiso), prof%cws(j,1), wiso%lost(2:nwiso), wiso%lost(1))
             if (wiso%lost(1)>zero) then
!               print *, wiso%lost(1:nwiso)
!               print *, prof%cws(j,2:nwiso)
             end if
#ifdef DEBUG
              if (any(abs(wiso%lost(1:nwiso)) > epsilon(one)) .and. soil%lost0 == 0) then
                 call message('CANVEG', 'Lost 03 @  ', num2str(time%daytime))
                 soil%lost0 = 1
              end if
#endif
           end if
           prof%cws(j,1) = max(zero, &
                prof%cws(j,1) - (prof%sun_LEwet(j,1)*solar%prob_beam(j) &
                + prof%shd_LEwet(j,1)*solar%prob_shd(j)) * prof%dLAIdz(j) &
                / lambda(prof%tair_filter_save(j) + TN0) * time%time_step)
!                print *, "canopy water storage", j, prof%cws(j,1)
           if (iswitch%wiso == 1) then
              ! No fractionation for interception evaporation
              prof%cws(j,2:nwiso) = rcws(2:nwiso) * prof%cws(j,1)
              ! Check: if this happens, we have to code some threshold here
              if (any(prof%cws(j,1:nwiso) < zero)) then
                 call message('CANVEG', 'prof%cws < 0 @ ', num2str(j))
              end if
           end if
        end if
     end do

     c_transpiration_mole = 1e6_wp * flux%c_transpiration(1) / 18._wp ! mmol m-2 s-1

     if (soil%camillo == 0) then
        ! compute litter moisture and litter drainage
        call litter_h2o()
        ! write(*,'(a,3f20.14)') 'CV53.01 ', flux%soilinfl(1), flux%soilinfl(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV53.02 ', soil%qinfl(1), soil%qinfl(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV53.03 ', soil%theta_l(1), soil%theta_l(nwiso-1)
        ! compute soil moisture in different layers
        call soil_h2o()
        ! write(*,'(a,3f20.14)') 'CV54.01 ', soil%qseva(1), soil%qseva(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV54.03 ', soil%qdrai(1), soil%qdrai(nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV54.04 ', soil%swp(1), soil%swp(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.05 ', soil%swp_mm(1), soil%swp_mm(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.06 ', soil%k_theta(1), soil%k_theta(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.07 ', soil%a(1), soil%a(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.08 ', soil%b(1), soil%b(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.09 ', soil%c(1), soil%c(nsoil)
        ! write(*,'(a,3f20.14)') 'CV54.10 ', soil%r(1,1), soil%r(nsoil,nwiso-1)
        ! write(*,'(a,3f20.14)') 'CV54.11 ', soil%soil_mm, soil%soil_mm_root, soil%soil_mm_50
        ! write(*,'(a,3f20.14)') 'CV54.13 ', soil%theta(1,1), soil%theta(nsoil,nwiso-1)
     endif

     ! check and convert units of all components to desirable values
     ! Convert to umol m-2 s-1
     can_ps_mg    = can_ps_mol * mass_CO2 *e3
     !wue          = can_ps_mg / (1000._wp * flux%c_transpiration(1)) ! mg co2 / g h20
     fc_mg        = -(can_ps_mg - soil%respiration_mg - bole%respiration_mg)
     fc_mol       = fc_mg * 1000._wp / mass_CO2

     sun_A(1:ncl) = max(prof%dLAIdz(1:ncl)*prof%sun_A(1:ncl)*solar%prob_beam(1:ncl), zero)
     shd_A(1:ncl) = max(prof%dLAIdz(1:ncl)*prof%shd_A(1:ncl)*solar%prob_shd(1:ncl),  zero)
     sumA         = sum(sun_A(1:ncl)) + sum(shd_A(1:ncl))
     sum_cica     = sum(sun_A(1:ncl)*prof%sun_cica(1:ncl)) + sum(shd_A(1:ncl)*prof%shd_cica(1:ncl))
     sum_gs       = sum(sun_A(1:ncl)*prof%sun_gs_mol(1:ncl)) + sum(shd_A(1:ncl)*prof%shd_gs_mol(1:ncl))
     if (sumA > zero) then
        ave_cica  = sum_cica / sumA
        ave_gs    = sum_gs / sumA
     else
        ave_cica  = undef
        ave_gs    = undef
     end if
     ! calculate total hourly O2 flux from the profile
     if (iswitch%oxygen == 1) then
        !sumo = can_gpp*ROC_leaf_in
        sumo = can_gpp_o
        sumcanresp_o = canresp_o
        sumneto = can_ps_o-soil%respiration_mole*ROC_soil_in-bole%respiration_mole*ROC_bole_in
        output%houro = sumo ! save hourly o flux in a global variable
        output%hourneto = sumneto
        nitrogen%Nsupply = sumNsupply
        nitrogen%Ndemand = sumNdemand
        output%hour_canrespo = sumcanresp_o
        if (fc_mol == zero) then
            output%hourROC = zero ! Yuan added hourly ROC output 2018.05.07
        else
            output%hourROC = abs(sumneto/fc_mol)
        end if
     end if
     ! 13C calculations
     if (iswitch%d13c == 1) then
        ! call carbon isotope subroutine
        ! it produces outputs based on inputs of Ci/Ca, gs, A, Tl
        call carbon_isotopes()
        ! integrated source fluxe of 13C
        sum13C = sum(prof%sour13co2(1:ncl))
        ! avg 13C discrimination, disc13, weighted by photosynthesis
        sumdisc13      = sum(sun_A(1:ncl)*prof%sun_disc13(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13(1:ncl))
        sumdisc13_long = sum(sun_A(1:ncl)*prof%sun_disc13_long(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_long(1:ncl))
        sumdisc13      = sum(sun_A(1:ncl)*prof%sun_disc13(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13(1:ncl))
        sumdisc13_long = sum(sun_A(1:ncl)*prof%sun_disc13_long(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_long(1:ncl))
        sumdisc13_a    = sum(sun_A(1:ncl)*prof%sun_disc13_a(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_a(1:ncl))
        sumdisc13_ab   = sum(sun_A(1:ncl)*prof%sun_disc13_ab(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_ab(1:ncl))
        sumdisc13_asal = sum(sun_A(1:ncl)*prof%sun_disc13_asal(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_asal(1:ncl))
        sumdisc13_b    = sum(sun_A(1:ncl)*prof%sun_disc13_b(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_disc13_b(1:ncl))
        sumdisc13_long_isotope = sum(prof%disc13C_long(1:ncl)*prof%dPsdz(1:ncl)) ! ???
        sum_csca       = sum((sun_A(1:ncl) + shd_A(1:ncl))*prof%csca(1:ncl))
        sum_cica       = sum(sun_A(1:ncl)*prof%sun_cica(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_cica(1:ncl))
        sum_ccca       = sum((sun_A(1:ncl) + shd_A(1:ncl))*prof%ccca(1:ncl))
        sum_gs         = sum(sun_A(1:ncl)*prof%sun_gs_mol(1:ncl)) &
             + sum(shd_A(1:ncl)*prof%shd_gs_mol(1:ncl))
        sum_gm         = sum((sun_A(1:ncl)+shd_A(1:ncl))*prof%gm(1:ncl))
        if (sumA > zero) then
           ztmp = one / sumA
           ave_disc13              = sumdisc13 * ztmp
           ave_disc13_long         = sumdisc13_long * ztmp
           ave_disc13_a            = sumdisc13_a * ztmp
           ave_disc13_ab           = sumdisc13_ab * ztmp
           ave_disc13_asal         = sumdisc13_asal * ztmp
           ave_disc13_b            = sumdisc13_b * ztmp
           ave_disc13_long_isotope = sumdisc13_long_isotope * ztmp
           ave_csca                = sum_csca * ztmp
           ave_cica                = sum_cica * ztmp
           ave_ccca                = sum_ccca * ztmp
           ave_gs                  = sum_gs * ztmp
           ave_gm                  = sum_gm * ztmp
        else
           ave_disc13              = undef
           ave_disc13_long         = undef
           ave_disc13_a            = zero
           ave_disc13_ab           = zero
           ave_disc13_asal         = zero
           ave_disc13_b            = zero
           ave_disc13_long_isotope = undef
           ave_csca                = undef
           ave_cica                = undef
           ave_ccca                = undef
           ave_gs                  = undef
           ave_gm                  = undef
        end if

        ! Flux of 13C in micromoles of 13C m-2 s-1
        ! it should relate to Rresp Reco + A (Rair/(1 + disc13)
        fc_13C = sum13C + soil%resp_13

        ! compute various integrated isotope factors
        if (can_ps_mol > zero) then
           do j=1, ncl
              ! compute flux weighted discrimination ratio for sun and shaded leaves
              sun_A(j) = max(prof%sun_A(j), zero)
              shd_A(j) = max(prof%shd_A(j), zero)
              prof%Rresp(j) = sun_A(j)*solar%prob_beam(j)*prof%Rplant_sun(j) &
                   + shd_A(j)*solar%prob_shd(j)*prof%Rplant_shd(j)
              prof%Rresp(j) = prof%Rresp(j) &
                   / (sun_A(j)*solar%prob_beam(j) + shd_A(j)*solar%prob_shd(j))
              ! sum the isotope ratio to define the respiratory signal
              ! based on the carbon that is assimilated
              if ((shd_A(j) > zero) .and. (sun_A(j) > zero)) then
                 prof%Rresp_sum(j) = prof%Rresp_sum(j) + prof%Rresp(j)
                 prof%cnt_Rresp(j) = prof%cnt_Rresp(j) + 1
              end if
           end do
        else
           prof%Rresp(1:ncl) = zero
        end if

        ! compute isotope flux
        ! computing profile of recycled CO2
        ! theta=(Cf (df-dR)- Ca(da-dR))/(D Cf - (da-dR)Ca)
        ave_daC13 = zero ! average del 13C of the air in the canopy
        !MC ave_Ca_f = zero ! average CO2 concentration in the forest
        !MC ave_d13resp = zero ! ave del 13C of respiring plant
        ztmp = one / real(ncl,wp)
        ave_daC13 = sum(prof%d13Cair(1:ncl)) * ztmp
        !print*, 'd13a ', prof%d13Cair(1:ncl)
        !MC ave_Ca_f = sum(prof%co2_air(1:ncl)) / ncl
        prof%d13Cplant(1:ncl) = delta1000(prof%Rresp(1:ncl), 1._wp, Rpdb_12C)
        !MC ave_d13resp = sum(prof%d13Cplant(1:ncl)) / ncl

        ! my only concern is weighting plants and soil correctly.
        ! compute ave_d13resp from a simple Keeling plot
        ! Feb 28, 2002 having problems with recyling, getting negative
        ! March 28, 2002 found algebraic error
        !MC recycle = ave_Ca_f*(ave_daC13-ave_d13resp) - input%co2air*(input%d13CO2-ave_d13resp)
        !MC recycle = recycle / ave_Ca_f*ave_disc13 - input%co2air*(input%d13CO2-ave_d13resp)
        ! Reco = Rsoil + Rbole
        ! Bowling et al. 2001 Global Change Biology 7, 127-145
        ! Fisotope = d13resp Reco +(d13air-disc13)A
        ! Need to adjust the signs with the direction of the fluxes to do the accounting correctly
        ! plus up like respiration, negative down like Ps
        ! Ecosystem respiration contains soil and bole respiration.
        !MC Fisotope = -26.32 * (soil%respiration_mole+bole%respiration_mole) + &
        !              (ave_daC13-ave_disc13)*(-can_ps_mol)
     end if ! 13C

     ! calculate canopy albedo
     ztmp = (solar%beam_flux_par(ncl+1) + solar%par_down(ncl+1))/4.6_wp &
          + solar%beam_flux_nir(ncl+1) + solar%nir_dn(ncl+1)
     if (ztmp /= zero) then
        flux%albedo = (solar%par_up(ncl+1)/4.6_wp + solar%nir_up(ncl+1)) / ztmp
     else
        flux%albedo = undef
     end if

     ! Isoprene flux density
     if (iswitch%isoprene == 1) isoprene_efflux = isoprene_canopy_flux()

     ! Daily sums of hourly data
     output%sumevap = output%sumevap + sumle ! latent heat flux
     output%sumsens = output%sumsens + sumh ! sensible heat flux
     output%sumfc   = output%sumfc + fc_mol ! CO2 flux, micromoles m-2 s-1
     output%sumpar  = output%sumpar + input%parin ! PAR
     output%sumnet  = output%sumnet + netrad ! Net radiation
     output%sumbole = output%sumbole + bole%respiration_mole ! bole respiration
     output%sumsoil = output%sumsoil + soil%respiration_mole
     output%sumresp = output%sumresp + canresp ! canopy respiration
     output%sumta   = output%sumta + tleaf_mean ! mean leaf temperature
     output%sumgs   = output%sumgs + sumksi ! canopy stomatal conductance
!     daily sums of hourly oxygen flux
      if (iswitch%oxygen == 1) then
         output%sumo = output%sumo + sumo
         output%sumneto = output%sumneto + sumneto
         output%sumresp_o = output%sumresp_o + canresp_o ! O2 via canopy respiration
      end if
     if (iswitch%d13c == 1) output%sumF13C = output%sumF13C + fc_13C ! net 13C flux
     daycnt = daycnt + 1
     output%sumps = output%sumps + can_ps_mol ! canopy photosynthesis
     if (can_gpp > 0) then
        output%sumgpp = output%sumgpp + can_gpp ! GPP
        output%sumTleaf = output%sumTleaf + tleaf_mean ! daytime leaf temperature
        gppcnt = gppcnt + 1
     end if
     if (iswitch%d13c == 1 .and. can_gpp > zero) then
        output%sumdisc13C = output%sumdisc13C + ave_disc13*can_ps_mol ! isotopic discrimination
        output%sumdisc13C_long_day = output%sumdisc13C_long_day + ave_disc13_long*can_ps_mol
     end if

     ! write to output structure
     output%netrad     = netrad
     output%sumrn      = sumrn
     output%sumlout    = sumlout
     output%sumh       = sumh
     output%sumle      = sumle
     output%can_ps_mol = can_ps_mol
     output%can_gpp    = can_gpp
     output%c_transpiration_mole = c_transpiration_mole
     output%canresp    = canresp
     output%sumksi     = sumksi
     output%tleaf_mean = tleaf_mean
     output%tavg_sun   = tavg_sun
     output%tavg_shd   = tavg_shd
     output%ave_cica   = ave_cica
     output%ave_gs     = ave_gs
     output%ave_disc13_a    = ave_disc13_a
     output%ave_disc13_ab   = ave_disc13_ab
     output%ave_disc13_asal = ave_disc13_asal
     output%ave_disc13_b    = ave_disc13_b
     output%ave_csca        = ave_csca
     output%ave_ccca        = ave_ccca
     output%ave_gm          = ave_gm
     output%isoprene_efflux = isoprene_efflux
     if (input%parin > 5._wp) then
        output%diff_par = solar%par_diffuse/input%parin
     else
        output%diff_par = undef
     end if
     output%rnet_soil = rnet_soil
     output%dLAIdz(1:ncl) = prof%dLAIdz(1:ncl)
     if (iswitch%d13c == 1) then
        output%ave_disc13      = ave_disc13
        output%ave_disc13_long = ave_disc13_long
        output%ave_daC13       = ave_daC13
     endif

     ! hourly output
     call write_output()
     ! output profiles for only specified periods
     if ((time%daytime >= start_profiles) .and. (time%daytime <= end_profiles)) call write_profiles()

     ! End of day
     ! Matthias: > or >= ?
     ! Yuan 2017: >
     if ((time%daytime+(nint(time%time_step,kind=i8)*100_i8/3600_i8)) & ! cf. mo_io_text.f90:read_text_in
          > (int(input%dayy,kind=i8)*10000_i8+2400)) then ! Yuan 20180104 >= to >
        ! compute daily averages
        if (gppcnt==0) then
           output%sumgpp   = zero
           output%sumTleaf = undef
        else
           output%sumgpp   = output%sumgpp   / gppcnt
           output%sumTleaf = output%sumTleaf / gppcnt
        end if
        if (iswitch%d13c==1) then
           if (gppcnt==0) then
              output%sumdisc13C          = undef
              output%sumdisc13C_long_day = undef
           else
              output%sumdisc13C          = output%sumdisc13C / output%sumps
              output%sumdisc13C_long_day = output%sumdisc13C_long_day / output%sumps
           end if
        end if
        if (daycnt==0) call finish('CANVEG','Unknown daycnt error')
        ztmp = one / real(daycnt,wp)
        output%sumfc   = output%sumfc   * ztmp
        output%sumevap = output%sumevap * ztmp
        output%sumsens = output%sumsens * ztmp
        output%sumpar  = output%sumpar  * ztmp
        output%sumnet  = output%sumnet  * ztmp
        output%sumbole = output%sumbole * ztmp
        output%sumsoil = output%sumsoil * ztmp
        output%sumta   = output%sumta   * ztmp
        output%sumgs   = output%sumgs   * ztmp
        output%sumresp = output%sumresp * ztmp
        output%sumps   = output%sumps   * ztmp
        ! oxyflux
        if (iswitch%oxygen==1) then ! write daily oxyflux and ROC
           output%sumo               = output%sumo   * ztmp
           output%sumneto            = output%sumneto   * ztmp
           output%sumresp_o          = output%sumresp_o * ztmp
           if (output%sumfc==zero) then
            output%sumROC = zero
           else
            output%sumROC = abs(output%sumneto/output%sumfc)
           end if

        end if
        ! 13C
        if (iswitch%d13c==1) then
           output%sumF13C = output%sumF13C * ztmp
           ! compute the previous days isotope ratio based on its photosynthesis
           ! TODO: sumps > 0 -> prof%cnt_Rresp(1:ncl) > 0
           if (output%sumps > 0) then
              prof%Rresp_ave(1:ncl) = prof%Rresp_sum(1:ncl)/prof%cnt_Rresp(1:ncl)
           else
              prof%Rresp_ave(1:ncl) = Rlongterm_12C
           end if
           prof%Rresp_sum(1:ncl) = zero
           prof%cnt_Rresp(1:ncl) = zero
           ! shift memory bigdelta one day back
           ciso%bigdelta(2:ndaysc13)      = ciso%bigdelta(1:ndaysc13-1)
           ciso%bigdelta_long(2:ndaysc13) = ciso%bigdelta_long(1:ndaysc13-1)
           if (output%sumps > zero) then
              ciso%bigdelta(1)      = output%sumdisc13C
              ciso%bigdelta_long(1) = output%sumdisc13C_long_day
           else  ! in case of no photosynthesis keep bigdelta from previous day
              ! TODO: replace ndaysc13 with 2
              ciso%bigdelta(1)      = ciso%bigdelta(ndaysc13)
              ciso%bigdelta_long(1) = ciso%bigdelta_long(ndaysc13)
           end if
        end if
        ! Daily output file
        call write_daily()
        ! Re-zero summing variables
        output%sumfc               = zero
        output%sumevap             = zero
        output%sumsens             = zero
        output%sumps               = zero
        output%sumgpp              = zero
        output%sumpar              = zero
        output%sumnet              = zero
        output%sumbole             = zero
        output%sumsoil             = zero
        output%sumta               = zero
        output%sumgs               = zero
        output%sumresp             = zero
        output%sumTleaf            = zero
        if (iswitch%d13c==1) then
           output%sumF13C             = zero
           output%sumdisc13C          = zero
           output%sumdisc13C_long_day = zero
        end if
        ! zero daily oxygen varibles
        if (iswitch%oxygen==1) then
           output%sumo               = zero
           output%sumresp_o          = zero
           output%sumneto             = zero
        end if
        gppcnt = 0
        daycnt = 0
#ifndef QUIET
        if (print_balance==0) then
           call message(trim(num2str(time%days,"(i5.5)")),advance='no')
           !write(nout,"(i5.5)",advance="no") time%days
           if (mod(time%days,24)==0) call message(advance='yes')
        end if
#endif
     end if ! end if new day

!     print*, 'T: ', prof%shd_tleaf_filter(1), prof%shd_tleaf_filter(40)
     ! write(*,'(a,3f20.14)') 'T: ', time%daytime, prof%shd_tleaf_filter(1), prof%shd_tleaf_filter(40)

     if (lastin == 1) exit ! interrupt timestep loop

  end do ! time steps
  if (print_balance==0) call message(advance='yes')

  !
  ! Finish
  !

  call close_files()

  deallocate(sun_A)
  deallocate(shd_A)
  deallocate(tmpntl)
  deallocate(tmpncl)
  deallocate(rcws)

#ifndef QUIET
  call cpu_time(ctime2)
  call date_and_time(values=dtime)
  call message()
  message_text = 'Done '//trim(num2str(totalcount,'(I6)'))//" time steps."
  call message(trim(message_text))
  message_text = 'Energy and carbon balance converged'//' ' &
       //trim(num2str(totalcount-count40,'(I6)'))//' times = ' &
       //trim(num2str(100._wp*real(totalcount-count40,wp)/real(totalcount,wp),'(F6.2)'))//"%."
  call message(trim(message_text))
  if (dtotalcount > 0 .and. ntotalcount > 0) then
     message_text = 'Converged '//trim(num2str(dtotalcount-dcount40,'(I6)')) &
          //' times during daytime = ' &
          //trim(num2str(100._wp*real(dtotalcount-dcount40,wp)/real(dtotalcount,wp),'(F6.2)')) &
          //"% and "//trim(num2str(ntotalcount-ncount40,'(I6)'))//' times during night time = ' &
          //trim(num2str(100._wp*real(ntotalcount-ncount40,wp)/real(ntotalcount,wp),'(F6.2)'))//"%."
     call message(trim(message_text))
  endif
  message_text = trim(num2str(dtime(3),'(I2.2)'))//"."//trim(num2str(dtime(2),'(I2.2)')) &
       //"."//trim(num2str(dtime(1),'(I4.4)'))//" "//trim(num2str(dtime(5),'(I2.2)')) &
       //":"//trim(num2str(dtime(6),'(I2.2)'))//":"//trim(num2str(dtime(7),'(I2.2)'))
  call message('Finished at ', trim(message_text), ' in ', trim(num2str(ctime2-ctime1,'(F12.3)')), ' seconds.')
  call message()
  call finish('Canveg','Finished!')
#endif

END PROGRAM canveg
