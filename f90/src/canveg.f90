PROGRAM canveg

  ! Parameters / Types / Switches / etc.
  USE kinds,         ONLY: wp, i4, i8                          ! Precision
  USE constants,     ONLY: zero, half, one, e2, e3, e6, &      ! numerical
       version, main_file, namelist_file, &                    ! info
       TN0, P0, Gravity, rugc, cp, &                           ! physical
       mass_air, mass_CO2, &
       RVSMOW_18O, RVSMOW_D, Rpdb_12C, &
       nnu, nuvisc, dc, ddh, dh, ddv, dv, &
       isnight, undef, judgenight                                          ! computational
  USE setup,         ONLY: ncl, ntl, nsoil, nwiso, ndaysc13    ! dimensions
  USE parameters,    ONLY: read_namelist, &                    ! namelist parameters
       longitude, latitude, ht, lai, &
       indir, outdir, workdir, start_run, end_run, &
       start_profiles, end_profiles, &
       extra_nate, lleaf, delz, bprime
  USE types,         ONLY: time, iswitch, soil, wiso, &        ! Canveg''s types
       non_dim, bound_lay_res, input, prof, solar, met, &
       ciso, fact, flux, output, bole, &
       ! Routines
       zero_new_timestep
  USE io,            ONLY: open_files, close_files, &          ! I/O wrappers
       read_disp, skip_input, lastin, write_daily, read_in, &
       write_profiles, write_output,write_debug
  USE soils,         ONLY: set_litter_texture, &               ! Soil & Litter
       set_soil_moisture, set_soil_temp, set_soil_root, &      !   water & energy
       set_soil_clapp, set_soil_texture, set_soil_saxton, &
       set_litter_moisture, set_litter_temp, set_soil_time, &
       set_soil_litter_capacity, soil_h2o, litter_h2o, &
       litter_rain, throughfall, soil_energy_balance
  USE energy_carbon, ONLY: bole_respiration, &                 ! Carbon & Energy
       soil_respiration, energy_and_carbon_fluxes
  USE radiation,     ONLY: set_leaf_phenology, lai_time, &     ! Phenology & Radiation
       irflux, angle, nir, par, rnet, &
       diffuse_direct_radiation, gfunc
  USE transport,     ONLY: conc, friction_velocity             ! Dispersion
  USE isotopes,      ONLY: leaf_wiso, soil_flux_wiso, &        ! Isotopes
       le_wiso, canopy_flux_wiso, carbon_isotopes
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
  logical(kind=4) ll ! Yuan 2017.08.14
  INTEGER(i8) :: d1=0
  INTEGER(i8) :: d2=0
  INTEGER(i8) :: d3=0
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

  !
  ! Read namelist
  !

#ifndef QUIET
  call message()
  call message('    Read namelist file: ', trim(namelist_file))
#endif
  call read_namelist()!key 1 = C 1606-1630 READ_PARAMETER()& SET_PARAMETER()

  ! Show model setup
#ifndef QUIET
  ! The following is not possible in Fortran90.
  !    write(nout,*) '    # canopy layers:     ', num2str(ncl)
  ! because the write statement calls num2str which itself has a write statement in it.
  ! e.g. NAG complains with: Recursive I/O reference.
  ! It should be possible in Fortran2003.
  ! However, MC''s NAG compiler is not allowing it even with the -f2003 compiler flag.
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
   wiso%vsmow(:) = one !in C model it is set to be zero Yuan 2017.09.25
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
  i_max = 101            ! maximum number of iterations for energy balance
  time%year = time%year0 ! Save first year
  if (.not. allocated(sun_A))  allocate(sun_A(ncl))
  if (.not. allocated(shd_A))  allocate(shd_A(ncl))
  if (.not. allocated(tmpntl))   allocate(tmpntl(ntl))
  if (.not. allocated(tmpncl)) allocate(tmpncl(ncl))
  if (.not. allocated(rcws))   allocate(rcws(nwiso))

  !
  ! Open files
  call open_files()
 ! call write_debug()
  !
  ! Setup leaves
  call set_leaf_phenology() ! set leaf onset and full
  !
  ! Setup soil
  call set_soil_texture()   ! set soil texture per soil layer
  call set_soil_root()      ! set rooting profile
  call set_soil_moisture()  ! set humidity per soil layer
  call set_soil_temp()      ! set deep soil temperature

  ! Setup litter
  call set_litter_texture()
  call set_litter_temp()
  call set_litter_moisture()
  if (soil%saxton==1) then  ! set hydraulic soil parameters
     call set_soil_saxton()
  else
     call set_soil_clapp()
  end if
  !
  !
  ! Relative plant available water
  soil%soil_mm_root = sum(soil%root(1:nsoil)*soil%theta(1:nsoil,1)*1000._wp * &
       (soil%z_soil(1:nsoil)-soil%z_soil(0:nsoil-1))*(one-soil%gravel(1:nsoil)*e2))
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

  ! Define parameters for soil energy balance model
  ! Preliminary tests show that the soil model can be sensitive to the depth of the litter layer
  call set_soil_time()

  ! Input data on Thomson dispersion matrix that was computed offline with
  ! MOVOAK%C, Dij (s m-1)
  call read_disp()

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
  end if

  ! Loop through the input met file. There should be a line of data for each hour of the year.
  call skip_input()
  time%days  = int(start_run / 10000_i8, i4)
  time%jdold = int(start_run / 10000_i8, i4)
  ! set input%dayy because of several years at once
  input%dayy = time%days ! for more then 1 year
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
!  print *, wiso%vsmow ! study later. Yuan 2017.09.24
  ! initialise leaf area index
  input%lai = lai
  if (extra_nate==1) then
     ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
     input%lai_up   = lai
     input%lai_down = lai
  end if
  call lai_time() ! define leaf area and canopy structure

  ! Initialize humidity and temperature profiles, arbitrary values
  prof%sun_tleaf(1:ncl)        = 25._wp
  prof%shd_tleaf(1:ncl)        = 25._wp
  prof%sun_tleaf_filter(1:ncl) = 25._wp
  prof%shd_tleaf_filter(1:ncl) = 25._wp
  tmpncl(1:ncl)                = one/(bprime(1:ncl)*(rugc*TN0/P0))
  prof%sun_rs(1:ncl)           = tmpncl(1:ncl)
  prof%sun_rs_filter(1:ncl)    = tmpncl(1:ncl)
  prof%shd_rs(1:ncl)           = tmpncl(1:ncl)
  prof%shd_rs_filter(1:ncl)    = tmpncl(1:ncl)

  prof%tair(1:ntl)                 = 25._wp
  prof%tair_filter(1:ntl)          = 25._wp
  prof%tair_filter_save(1:ntl)     = 25._wp
  prof%rhov_air(1:ntl,1)           = 0.02_wp
  prof%rhov_air_save(1:ntl)        = 0.02_wp
  prof%rhov_air_filter(1:ntl,1)    = 0.02_wp
  prof%rhov_air_filter_save(1:ntl) = 0.02_wp
  prof%co2_air(1:ntl)              = 380._wp
  prof%co2_air_filter(1:ntl)       = 380._wp

  ! Initialise atmospheric d13C
  if (iswitch%d13c==1) then
     prof%R13_12_air(1:ntl) = invdelta1000(-8._wp)*Rpdb_12C ! ratio of 13C relative to 12C
     prof%d13Cair(1:ntl)    = -8.0_wp
  end if

  ! Initialise atmospheric water isotopes
  forall(mc=2:nwiso)!2017.10.04
     prof%rhov_air(1:ntl,mc)        = 0.02_wp*wiso%vsmow(mc)
     prof%rhov_air_filter(1:ntl,mc) = 0.02_wp*wiso%vsmow(mc)
     prof%rvapour(1:ntl,mc)         = wiso%vsmow(mc)
  end forall
  !print *, "wiso%vsmow", wiso%vsmow
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
#endif
  do ! loop over time steps until lastin==1

     ! init new time step
     if (nwiso > 1) then
        prof%sun_leafwater_e_old(1:ncl,2:nwiso) = prof%sun_leafwater_e(1:ncl,2:nwiso)
        prof%shd_leafwater_e_old(1:ncl,2:nwiso) = prof%shd_leafwater_e(1:ncl,2:nwiso)
     end if
     call zero_new_timestep() ! zero variables
     call read_in()           ! input new time step
  !  if (time%daytime==1241200) then  !Yuan 2017.09.20

!	      print *, "what an hour!!!\n"
!	      print *, "daytime = ", time%daytime
!		  print *, "\n"
!	 end if
     ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
     ! update LAI in each time step
 !    d3 = time%daytime-(int(time%time_step,kind=i8)*100_i8/3600_i8)
     d3 = int(time%time_step,kind=i8)*100_i8/3600_i8
     if (extra_nate==1) then
        call lai_time()
     else
        ! start of new day
        ! if ((time%daytime+(nint(time%time_step,kind=i8)*100_i8/3600_i8)) & ! cf. mo_io_text.f90:read_text_in
        !      >= (int(input%dayy,kind=i8)*10000_i8+2400)) then
        if ((time%daytime-(int(time%time_step,kind=i8)*100_i8/3600_i8)) &
             <= (int(input%dayy,kind=i8)*10000_i8)) call lai_time()
     endif

     ! Compute solar elevation angle
     call angle()
     judgenight =((solar%sine_beta <= isnight) .OR. (input%parin <=0))
     ! make sure PAR is zero at night. Some data have negative offset, which causes numerical problems
     if (judgenight) input%parin = zero
     ! Compute the fractions of beam and diffuse radiation from incoming measurements
     ! Set the radiation factor to the day before for night calculations. This way if the day
     !   was cloudy, so will IR calculations for the night reflect this.
     if (.NOT.judgenight) then
        call diffuse_direct_radiation()
     else
        solar%ratrad      = solar%ratradnoon
  !      solar%nir_beam    = 0.1_wp ! Should not be here, to agree with C version. Yuan 2017.09.24
  !      solar%nir_diffuse = 0.1_wp
     end if
     ! computes leaf inclination angle distribution function, the mean direction cosine
     ! between the sun zenith angle and the angle normal to the mean leaf
     ! for CANOAK we use the leaf inclination angle data of Hutchison et al. 1983, J Ecology
     ! for other canopies we assume the leaf angle distribution is spherical
     if (.NOT.judgenight) then
        call gfunc()
     else
        prof%Gfunc_solar(1:ncl) = 0.01_wp
     end if
     ! Set soil reflectivity depending on soil moisture of first layer
     !   i.e. wet soil seems half as bright as dry soil
     !   after Wilson & Henderson-Sellers (1985)
     !   cited in Knorr PhD (1997, Table 2.4, p. 34)
     ztmp = min(max((soil%theta(1,1)-soil%watmin(1))/(soil%theta_s(1)-soil%watmin(1)), zero), one)
     solar%par_soil_refl = solar%par_soil_refl_dry * (one-half*ztmp)
     solar%nir_soil_refl = solar%nir_soil_refl_dry * (one-half*ztmp)
     ! Compute PAR profiles
     call par()
     ! Compute NIR profiles
     call nir()
     ! Interception reservoir
     call throughfall()
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

     ! iteration looping for energy fluxes and scalar fields
     ! iterate until energy balance closure occurs or i_max iterations are reached
     do
        ! save input filtered variables for use in wiso routines
        prof%tair_filter_save(1:ntl) = prof%tair_filter(1:ntl)
        prof%rhov_air_filter_save(1:ntl) = prof%rhov_air_filter(1:ntl,1)
        prof%rhov_air_save(1:ntl) = prof%rhov_air(1:ntl,1)
!print *, prof%rhov_air_save(1:ntl)
        call irflux() ! first guess
        call friction_velocity()
        ! compute net radiation balance on sunlit and shaded leaves
        call rnet()
        ! print*, 'CA11 ', i_count, solar%par_sun(1), solar%nir_sun(1)
        ! print*, 'CA12 ', i_count, solar%par_sun(40), solar%nir_sun(40)
        ! print*, 'CA13 ', i_count, solar%par_shd(1), solar%nir_shd(1)
        ! print*, 'CA14 ', i_count, solar%par_shd(40), solar%nir_shd(40)
        ! print*, 'CA15 ', i_count, solar%ir_dn(1), solar%ir_up(1)
        ! print*, 'CA16 ', i_count, solar%ir_dn(40), solar%ir_up(40)
        ! print*, 'CA17 ', i_count, solar%rnet_sun(1), solar%rnet_shd(1)
        ! print*, 'CA18 ', i_count, solar%rnet_sun(40), solar%rnet_shd(40)
        call energy_and_carbon_fluxes()
        ! Compute leaf energy balance, leaf temperature, photosynthesis and stomatal conductance        call energy_and_carbon_fluxes()
        ! print*, 'CA01 ', i_count, prof%dLEdz(1,1), prof%dLEdz(40,1)
        ! print*, 'CA02 ', i_count, prof%dHdz(1), prof%dHdz(40)
        ! print*, 'CA03 ', i_count, prof%dRNdz(1), prof%dRNdz(40)
        ! print*, 'CA04 ', i_count, prof%dLoutdz(1), prof%dLoutdz(40)
        ! print*, 'CA05 ', i_count, prof%shd_tleaf(1), prof%shd_tleaf(40)
        ! Soil energy balance
        call soil_energy_balance()
        ! print*, 'CA06 ', flux%soilevap(1), flux%litterevap(1), flux%s_evap(1)
        ! print*, 'CA07 ', soil%tsrf, soil%T_l, soil%litterevap
        ! print*, 'CA08 ', output%c1, output%c2, output%c3
        ! print*, 'CA09 ', output%c4, output%c7
        ! print*, 'CA10 ', soil%theta(1:3,1)

        ! --- Water ---
        if (iswitch%wiso==1) then
           ! calculates latent heat in accordance with leaf water isotope calculations
           call le_wiso()
           ! re-define/re-calc LEstoma in accordance with leaf water isotope calculations
           prof%sun_LEstoma_save(1:ncl) = prof%sun_LEstoma(1:ncl,1)
           prof%shd_LEstoma_save(1:ncl) = prof%shd_LEstoma(1:ncl,1)
           prof%sun_LEstoma(1:ncl,1)    = prof%sun_LEstoma_new(1:ncl) &
                * lambda(prof%tair_filter(1:ncl)+TN0)
           prof%shd_LEstoma(1:ncl,1)    = prof%shd_LEstoma_new(1:ncl) &
                * lambda(prof%tair_filter(1:ncl)+TN0)
           prof%dLEdz(1:ncl,1)          = &
                (solar%prob_beam(1:ncl) * (prof%sun_LEstoma(1:ncl,1)+prof%sun_LEwet(1:ncl,1)) &
                + solar%prob_shd(1:ncl) * (prof%shd_LEstoma(1:ncl,1)+prof%shd_LEwet(1:ncl,1))) &
                *  prof%dLAIdz(1:ncl)
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

        ! Compute air temperature profiles from dispersion matrix
        ! Adjust dHdz(1) for soil heat flux
        ! sign convention used: fluxes from surface are positive
        ! those toward the surface are negative
        ! dHdz is net sensible heat exchange per unit leaf area
        ! Using filtered temperatures to minimize the system from being mathematically unstable.

        ! filter temperatures with each interation to minimize numerical instability

        ! Matthias, constant filter
          fact%a_filt = 0.85_wp !Yuan2017.10.05
        ! Matthias, force conversion with narrowing filter
        !   fact%a_filt = 0.99_wp - 0.98_wp*real(i_count,wp)/real(i_max-1,wp)
        ! Matthias, narrowing filter to 0.5 at i_max/2
         fact%a_filt = 0.85_wp - 0.7_wp*real(i_count,wp)/real(i_max-1,wp)!Yuan 2017.10.05
        !fact%a_filt = 0.85_wp - 0.7_wp*real(i_count,wp)/real(i_max-1,wp)!replace 0.85 and 0.7 with two parameters that calibrated by Yuan 2017.11.02
        !fact%a_filt = fact%a_filt*0.7_wp
        ! Matthias, a_filt=0 from 10 steps before max iteration
        ! the noise one sees should be filtered off as well
        !   fact%a_filt = max(0.99_wp - one*real(i_count,wp)/real(i_max-10,wp), zero)
        ! Matthias, original filter
        !   if (i_count < 10) then
        !     fact%a_filt = 0.85_wp
        !   else
        !     fact%a_filt = half
        !   end if

        ! conc, for temperature profiles using source/sinks
        ! inputs are source/sink(), scalar_profile(), ref_val, boundary_flux, unit conversions
        call conc(prof%dHdz, prof%tair, input%ta, soil%heat, fact%heatcoef)


        ! filter temperatures to remove numerical instabilities for each iteration
        where (prof%tair(1:ntl) < -30._wp .or. prof%tair(1:ntl) > 60._wp) prof%tair(1:ntl) = input%ta

        ! Compute vapor density profiles from Dispersion matrix
        ! sign convention used: fluxes from surface are positive
        ! those toward the surface are negative
        ! dLEdZ is net latent heat exchange per unit leaf area
        ! turbulent transport of H2O
        tmpncl(1:ncl) = prof%dLEdz(1:ncl,1) / lambda(prof%tair_filter_save(1:ncl)+TN0)
        call conc(tmpncl, tmpntl, met%rhova_kg, flux%s_evap(1), one)
        prof%rhov_air(1:ntl,1) = tmpntl(1:ntl)

        ! filter humidity computations
        where (prof%rhov_air(1:ntl,1) < zero .or. prof%rhov_air(1:ntl,1) > e2) &
             prof%rhov_air(1:ntl,1) = met%rhova_kg
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
        ! Implicit water isotopes diagnostics
        if (iswitch%wiso==1 .and. wiso%implicit==1) then
           ! isotope soil water flux
           call soil_flux_wiso()
           ! leaf water enrichment
           call leaf_wiso()
           ! isotope canopy transpiration and evaporation
           call canopy_flux_wiso()
           ! turbulent transport of H2O18, DHO and H216O
           do mc=2, nwiso
              ! Vapour input in layers
              tmpncl(1:ncl) = prof%dLEdz(1:ncl,mc) / lambda(prof%tair_filter_save(1:ncl)+TN0)
              ! Vapour input at top
              ztmp = alpha_equ_h2o(input%ta+TN0, mc) * invdelta1000_h2o(input%dppt(mc), mc) * met%rhova_kg ! in equi with rain
              call conc(tmpncl, tmpntl, ztmp, flux%s_evap(mc), one)
              prof%rhov_air(1:ntl,mc) = tmpntl(1:ntl)
           end do
           ! filter humidity computations
           do mc=2, nwiso
              where (prof%rhov_air(1:ntl-1,1) == met%rhova_kg) &
                   prof%rhov_air(1:ntl-1,mc) = prof%rhov_air(1:ntl-1,1) * prof%rvapour(1:ntl-1,mc)
           end do
        end if ! implicit water isotope diagnostics

        ! sign convention used: photosynthetic uptake is positive
        ! respiration is negative
        ! prof%dPsdz is net photosynthesis per unit leaf area,
        ! the profile was converted to units of mg m-2 s-1 to be
        ! consistent with inputs to CONC
        ! change sign of dPsdz
        prof%source_co2(1:ncl) = -prof%dPsdz(1:ncl)

        ! compute bole respiration
        ! It is added to the prof%source_CO2 array.
        call bole_respiration()
        ! compute soil respiration
        call soil_respiration()

        ! to convert umol m-3 to umol/mol we have to consider
        ! Pc/Pa = (CO2)ppm = rhoc ma/ rhoa mc
        fact%co2 = (mass_air/mass_CO2)*met%air_density_mole
        call conc(prof%source_co2, prof%co2_air, input%co2air, soil%respiration_mole, fact%co2)

        ! Integrate source-sink strengths to estimate canopy flux
        sumh       = sum(prof%dHdz(1:ncl))       ! sensible heat
        sumle      = sum(prof%dLEdz(1:ncl,1))    ! latent heat
        sumrn      = sum(prof%dRNdz(1:ncl))      ! net radiation
        sumlout    = sum(prof%dLoutdz(1:ncl))    ! radiation
        can_ps_mol = sum(prof%dPsdz(1:ncl))      ! canopy photosynthesis
        can_gpp    = sum(prof%dGPPdz(1:ncl))     ! canopy GPP = can_ps_mol + dark respiration
        canresp    = sum(prof%dRESPdz(1:ncl))    ! canopy respiration
        sumksi     = sum(prof%dStomCondz(1:ncl)) ! canopy stomatal conductance
        sumlai     = sum(prof%dLAIdz(1:ncl))     ! leaf area
!        print*, 'T: ', i_count, prof%shd_tleaf(1), prof%shd_tleaf(40)
        tleaf_mean = sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
             + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))    ! mean leaf temperature
        ztmp = one / real(ncl,wp)
        prof%tleaf(1:ncl) = (sum(prof%sun_tleaf(1:ncl)*solar%prob_beam(1:ncl)) &
             + sum(prof%shd_tleaf(1:ncl)*solar%prob_shd(1:ncl))) * ztmp    ! Tleaf per layer (sun and shade)
        ! need to weight by sun and shaded leaf areas then divide by LAI
        tavg_sun   = sum(prof%sun_tleaf(1:ncl)*prof%dLAIdz(1:ncl)) ! avg sunlit temperature
        tavg_shd   = sum(prof%shd_tleaf(1:ncl)*prof%dLAIdz(1:ncl)) ! avg shaded temperature

        ebalance      = sumrn - sumle - sumh
        flux%photosyn = can_ps_mol
        ! calculate gpp : can_ps_mol = photosynthesis - photorespiration
        !                              - dark respiration or can_ps_mol
        !                            = GPP - leaf respiration
        ! can_gpp = can_ps_mol + canresp

        ! mean canopy leaf temperature
        tleaf_mean = tleaf_mean / real(ncl,wp)
        ! leaf area weighted temperatures
        tavg_sun = tavg_sun / sumlai
        tavg_shd = tavg_shd / sumlai
        ! Energy exchanges at the soil
        rnet_soil = soil%rnet - soil%lout
        sbalance  = rnet_soil - soil%soilevap - soil%litterevap - soil%heat - soil%gsoil
        ! canopy scale flux densities, vegetation plus soil
        sumh     = sumh + soil%heat
        sumle    = sumle + soil%evap
        sumrn    = sumrn + rnet_soil
        sumlout  = sumlout + soil%lout
        temp3    = temp3 + soil%rnet
        tbalance = sbalance + ebalance
        met%H    = sumh

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
        ! met variables
        met%H_filter     = fact%a_filt * met%H + (one-fact%a_filt) * met%H_filter
        met%ustar_filter = fact%a_filt * met%ustar + (one-fact%a_filt) * met%ustar_filter
        ! air variables
        prof%tair_filter(1:ntl)       = fact%a_filt * prof%tair(1:ntl) & !2017.10.04
             + (one-fact%a_filt) * prof%tair_filter(1:ntl)
        prof%rhov_air_filter(1:ntl,1) = fact%a_filt * prof%rhov_air(1:ntl,1) &
             + (one-fact%a_filt) * prof%rhov_air_filter(1:ntl,1)
        prof%co2_air_filter(1:ntl)    = fact%a_filt * prof%co2_air(1:ntl) &
             + (one-fact%a_filt) * prof%co2_air_filter(1:ntl)
        ! leaf variables
        prof%sun_tleaf_filter(1:ncl) = fact%a_filt * prof%sun_tleaf(1:ncl) &
             + (one-fact%a_filt) * prof%sun_tleaf_filter(1:ncl)
        prof%shd_tleaf_filter(1:ncl) = fact%a_filt * prof%shd_tleaf(1:ncl) &
             + (one-fact%a_filt) * prof%shd_tleaf_filter(1:ncl)
        prof%wet_coef_filter(1:ncl)  = fact%a_filt * prof%wet_coef(1:ncl) &
             + (one-fact%a_filt) * prof%wet_coef_filter(1:ncl)

        ! Implicit water isotopes diagnostics
        if (iswitch%wiso==1 .and. wiso%implicit==1) then
           prof%rhov_air_filter(1:ntl,2:nwiso) = fact%a_filt*prof%rhov_air(1:ntl,2:nwiso) &
                + (one-fact%a_filt)*prof%rhov_air_filter(1:ntl,2:nwiso)
           ! isotopes in vapour
           prof%rvapour(1:ntl,2:nwiso) = isorat(prof%rhov_air_filter(1:ntl,2:nwiso), &
                prof%rhov_air_filter(1:ntl,1), wiso%lost(2:nwiso), wiso%lost(1))
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

        ! test for convergence between the sum of the net radiation flux profile and the
        ! net flux exiting the canopy
        ! etest = abs((sumrn-netrad)/sumrn)
        etest       = (sumrn-netrad) / sumrn
        etest_diff1 = etest_old1 - etest
        etest_diff2 = etest_old2 - etest_old1
        etest_old2  = etest_old1
        etest_old1  = etest

        itest = zero
        if (iswitch%wiso==1 .and. wiso%implicit==1) &
             itest = delta1000_h2o(flux%evapotranspiration(2), flux%evapotranspiration(1), 2)
        itest_diff1 = itest_old1 - itest
        itest_diff2 = itest_old2 - itest_old1
        itest_old2  = itest_old1
        itest_old1  = itest

        i_count    = i_count + 1
        time%count = i_count

        ! check for convergence
        ! if (etest <= 0.005_wp .or. i_count >= i_max) exit
        if ((abs(etest_diff1) <= 0.001_wp .and. abs(etest_diff2) <= 0.001_wp &
             .and. abs(itest_diff1) <= 0.01_wp .and. abs(itest_diff2) <= 0.01_wp) &
             .or. i_count >= i_max) exit
        ! if ((abs(etest_diff1) <= 1e-15_wp .and. abs(etest_diff2) <= 1e-15_wp &
        !      .and. abs(itest_diff1) <= 1e-15_wp .and. abs(itest_diff2) <= 1e-15_wp) &
        !      .or. i_count >= i_max) exit
     end do ! until energy balance closure in this time step (or i_max)

     totalcount = totalcount + 1
     if (judgenight) then
        ntotalcount = ntotalcount + 1
     else
        dtotalcount = dtotalcount + 1
     end if
     if (i_count == i_max) then
        count40 = count40 + 1
        if (judgenight) then
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
        prof%rhov_air_filter(1:ntl,2:nwiso) = spread(prof%rhov_air_filter_save(1:ntl),2,nwiso) * prof%rvapour(1:ntl,2:nwiso)
        ! isotope soil water flux
        call soil_flux_wiso()
        ! leaf water enrichment
        call leaf_wiso()
        ! isotope canopy transpiration and evaporation
        call canopy_flux_wiso()
        ! turbulent transport of H2O18, DHO and H216O
        do mc=2, nwiso
           ! Vapour input in layers
           tmpncl(1:ncl) = prof%dLEdz(1:ncl,mc) / lambda(prof%tair_filter_save(1:ncl)+TN0)
           ! Vapour input at top
           ztmp = alpha_equ_h2o(input%ta+TN0, mc) * invdelta1000_h2o(input%dppt(mc), mc) * met%rhova_kg ! in equi with rain
           call conc(tmpncl, tmpntl, ztmp, flux%s_evap(mc), one)
           prof%rhov_air(1:ntl,mc) = tmpntl(1:ntl)
        end do
        ! TODO: check order of statements
        ! filter humidity computations
        do mc=2, nwiso
           where (prof%rhov_air(1:ntl-1,1) == met%rhova_kg) &
                prof%rhov_air(1:ntl-1,mc) = prof%rhov_air(1:ntl-1,1) * prof%rvapour(1:ntl-1,mc)
        end do
        ! isotopes in vapour
        prof%rvapour(1:ntl,2:nwiso) = isorat(prof%rhov_air(1:ntl,2:nwiso), &
             prof%rhov_air(1:ntl,1), wiso%lost(2:nwiso), wiso%lost(1))
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
        ! compute soil moisture in different layers
        call soil_h2o()
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

     ! output
     call write_output()
!     call write_debug()
     ! output profiles for only specified periods
     ll = time%daytime >= start_profiles .and. time%daytime <= end_profiles
     d1 = time%daytime-start_profiles
     d2 = time%daytime-end_profiles
     if (ll) then
     call write_profiles()
!     call write_debug()
 !    else
!        call write_debug()
    endif

     ! End of day
     if ((time%daytime+(nint(time%time_step,kind=i8)*100_i8/3600_i8)) & ! cf. mo_io_text.f90:read_text_in
          > (int(input%dayy,kind=i8)*10000_i8+2400)) then
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

!     print*, 'T: ', time%daytime, prof%shd_tleaf_filter(1), prof%shd_tleaf_filter(40)

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
  message_text = 'Done '//trim(num2str(totalcount,'(I6.6)'))//" time steps."
  call message(trim(message_text))
  message_text = 'Energy and carbon balance converged'//' ' &
       //trim(num2str(totalcount-count40,'(I6.6)'))//' times = ' &
       //trim(num2str(100._wp*real(totalcount-count40,wp)/real(totalcount,wp),'(F6.2)'))//"%."
  call message(trim(message_text))
  if (dtotalcount > 0 .and. ntotalcount > 0) then
     message_text = 'Converged '//trim(num2str(dtotalcount-dcount40,'(I6.6)')) &
          //' times during daytime = ' &
          //trim(num2str(100._wp*real(dtotalcount-dcount40,wp)/real(dtotalcount,wp),'(F6.2)')) &
          //"% and "//trim(num2str(ntotalcount-ncount40,'(I6.6)'))//' times during night time = ' &
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
