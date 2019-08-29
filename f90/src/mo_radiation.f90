MODULE radiation

  ! This module contains the radiation transfer and canopy structure routines of Canveg

  ! Written Jan 2011, Mathias Cuntz - Ported C-Code

  USE kinds, ONLY: wp, i4

  IMPLICIT NONE
  LOGICAL judgenight ! judge if night, more complex than isnight. Yuan 2018.01.16
  COMMON  judgenight
  PRIVATE

  PUBLIC :: angle                     ! computes sun angles
  PUBLIC :: diffuse_direct_radiation  ! computes direct and diffuse radiation from radiation inputs
  PUBLIC :: freq                      ! leaf angle frequency distribution for G
  PUBLIC :: gfunc                     ! computes leaf orientation angle for the given sun angle and a spherical leaf inclination
  PUBLIC :: gfunc_diffuse             ! computes the leaf orientation function for diffuse radiation
  PUBLIC :: irflux                    ! computes longwave radiation flux density
  PUBLIC :: lai_time                  ! updates LAI with time
  PUBLIC :: nir                       ! computes near infrared radiation profiles
  PUBLIC :: par                       ! computes visible light profiles
  PUBLIC :: rnet                      ! computes net radiation profile
  PUBLIC :: set_leaf_phenology        ! initialize leaf start and leaf full
  PUBLIC :: sky_ir                    ! computes the incoming longwave radiation flux density

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------
  SUBROUTINE angle()
    ! ANGLE computes solar elevation angles,
    ! This subroutine is based on algorithms in Walraven. 1978. Solar Energy. 20: 393-397
    USE types,      ONLY: time, solar
    USE constants,  ONLY: zero, one, two, pi, pi2, pi180
    USE parameters, ONLY: longitude, latitude, zone
    USE messages,   ONLY: message

    IMPLICIT NONE

    REAL(wp) :: delyr, leap_yr, delyr4, day_local, leaf_yr_4
    REAL(wp) :: theta_angle, G, EL, EPS, sin_el, A1, A2, RA
    REAL(wp) :: T_local, time_1980
    REAL(wp) :: day_savings_time
    REAL(wp) :: S, HS, phi_lat_radians, val, declination_ang, ST, SSAS
    REAL(wp) :: E_ang, zenith, elev_ang_deg, cos_zenith

    E_ang            = zero
    declination_ang  = zero
    delyr            = time%year - 1980._wp
    delyr4           = delyr/4._wp
    leap_yr          = mod(delyr4,4._wp)
    ! Daylight Savings Time, Dasvtm =1
    ! Standard time, Dasvtm= 0
    day_savings_time = zero
    T_local          = time%local_time
    day_local        = real(time%days, kind=wp)
    time_1980        = delyr*365._wp + leap_yr + (day_local-one) + T_local / 24._wp
    leaf_yr_4        = leap_yr*4
    if (delyr == leaf_yr_4) time_1980 = time_1980 - one
    if ((delyr < 0) .and. (delyr /= leaf_yr_4)) time_1980 = time_1980 - one
    theta_angle = (360.0_wp * time_1980 / 365.25_wp) * pi180
    G           = -0.031272_wp - 4.53963E-7_wp * time_1980 + theta_angle
    EL          = 4.900968_wp + 3.6747E-7_wp * time_1980 + &
         (0.033434_wp - 2.3E-9_wp * time_1980) * sin(G) + &
         0.000349_wp * sin(two*G) + theta_angle
    EPS         = 0.40914_wp - 6.2149E-9_wp * time_1980
    sin_el      = sin(EL)
    A1          = sin_el * cos(EPS)
    A2          = cos(EL)
    !  for ATAN2
    ! RA = atan2(A1,A2)
    ! if (RA < zero) RA = RA + pi2
    RA = atan(A1/A2)
    ! The original program was in FORTRAN. It used the ATAN2 function.
    ! In C we must find the correct quadrant and evaluate the arctan correctly.
    ! Note ATAN2 is -PI TO PI, while ATAN is from PI/2 TO -PI/2
    ! QUAD II, TAN theta_angle < 0
    if (A1 > zero .and. A2 <= zero) RA = RA + pi
    ! QUAD III, TAN theta_angle > 0
    if (A1 <= zero .and. A2 <= zero) RA = Ra + pi
    val = sin_el * sin(EPS)
    if ((one-val*val) > zero) then
       declination_ang = atan(val/sqrt(one-val*val))
    else
       call message('ANGLE: ','bad declination_ang')
    end if
    ! declination_ang=asin(val)
    ST              = 1.759335_wp + pi2 * (time_1980 / 365.25_wp - real(delyr,kind=wp)) + 3.694E-7_wp * time_1980
    if (ST >= pi2) ST = ST - pi2
    S               = ST - longitude * pi180 + 1.0027379_wp * (zone - day_savings_time + T_local) * 15._wp * pi180
    if (S >= pi2) S = S - pi2
    HS              = RA - S
    phi_lat_radians = latitude * pi180
    ! DIRECTION COSINE
    SSAS = sin(phi_lat_radians) * sin(declination_ang) + cos(phi_lat_radians) * cos(declination_ang) * cos(HS)
    if ((one-SSAS*SSAS) > zero) then
       E_ang = atan(sqrt(one-SSAS*SSAS)/SSAS)
    else
       call message('ANGLE: ','bad SSAS')
    end if
    if (SSAS < zero) E_ang = E_ang + pi
    ! E=asin(SSAS)
    if (E_ang < zero) E_ang = pi/two
    zenith          = E_ang / pi180
    elev_ang_deg    = 90._wp - zenith
    solar%beta_rad  = elev_ang_deg * pi180
    solar%sine_beta = sin(solar%beta_rad)
    cos_zenith      = cos(zenith)
    solar%beta_deg  = solar%beta_rad / pi180

  END SUBROUTINE angle


  ! ------------------------------------------------------------------
  SUBROUTINE diffuse_direct_radiation()
    ! This subroutine uses the Weiss-Norman (1985, Agric. forest Meteorol. 34: 205-213)
    ! routine to compute direct and diffuse PAR from total par
    USE types,      ONLY: solar, input, time, output
    USE constants,  ONLY: zero, one, two, e3
    USE parameters, ONLY: extra_nate

    IMPLICIT NONE

    REAL(wp) :: fAND
    REAL(wp) :: rdir, rdvis, rsvis, wa
    REAL(wp) :: ru, rsdir, rvt, rit, nirx
    REAL(wp) :: xvalue, fvsb, fvd, fansb

    ! fractions of NIR and PAR (visible)
    !fir   = 0.54_wp
    !fv    = 0.46_wp
    ru    = 98.5_wp / (101.3_wp * solar%sine_beta)
    ! visible direct PAR
    rdvis = 600._wp * exp(-0.185_wp * ru) * solar%sine_beta
    ! potential diffuse PAR
    ! rsvis = .4 * (600.0 - rdvis) * solar%sine_beta
    rsvis = 0.4_wp * (600._wp * solar%sine_beta - rdvis) !new equation after Al Weiss, 1/14/05
    !solar constant: 1320 W m-2
    !water absorption in NIR for 10 mm precip water
    wa    = 1320._wp * 0.077_wp * (two*ru)**0.3_wp
    ! direct beam NIR
    rdir  = (720._wp * exp(-0.06_wp * ru) - wa) * solar%sine_beta
    if (rdir < zero) rdir = zero
    !potential diffuse NIR
    !rsdir = .6 * (720 - wa - rdir) * solar%sine_beta
    rsdir = 0.6_wp * (720._wp - wa - rdir/solar%sine_beta) * solar%sine_beta !new equation after Al Weiss, 1/14/05
    if (rsdir < zero) rsdir = zero
    rvt  = rdvis + rsvis
    rit  = rdir + rsdir
    if (rit < 0.1_wp) rit = 0.1_wp
    if (rvt < 0.1_wp) rvt = 0.1_wp
    solar%ratrad = input%rglobal / (rvt + rit)
    output%c10   = (rvt + rit)
    if (time%local_time >= 12._wp .and. time%local_time < 13._wp) solar%ratradnoon = solar%ratrad
    ! ratio is the ratio between observed and potential radiation
    ! NIR flux density as a function of PAR
    ! since NIR is used in energy balance calculations
    ! convert it to W m-2: divide PAR by 4.6
    nirx = input%rglobal - (input%parin / 4.6_wp)
    ! ratio = (PARIN / 4.6 + NIRX) / (rvt + rit)
    if (solar%ratrad > 0.89_wp) solar%ratrad = 0.89_wp
    if (solar%ratrad < 0.22_wp) solar%ratrad = 0.22_wp
    ! fraction PAR direct and diffuse
    xvalue = (0.9_wp-solar%ratrad)/0.70_wp
    fvsb   = rdvis / rvt * (one - xvalue**0.67_wp)
    if (fvsb < zero) fvsb = zero
    if (fvsb > one) fvsb = one
    fvd = one - fvsb
    ! note PAR has been entered in units of umol m-2 s-1
    ! use modeled direct and diffuse relationship
    solar%par_beam    = fvsb * input%parin
    solar%par_diffuse = fvd * input%parin
    ! use calculated PAR for Nate McDowell''s Juniper data
    if (extra_nate == 0) then ! use input%data
       solar%par_beam    = input%parin-input%pardif
       solar%par_diffuse = input%pardif
    end if
    if (solar%par_beam <= zero) then
       solar%par_beam    = zero
       solar%par_diffuse = input%parin
    end if
    if (input%parin <= e3) then
       solar%par_beam    = e3
       solar%par_diffuse = e3
    end if
    xvalue=(0.9_wp-solar%ratrad)/0.68_wp
    fansb = rdir / rit * (one - xvalue**0.67_wp)
    if (fansb < zero) fansb = zero
    if (fansb > one) fansb = one
    fAND = one - fansb
    !NIR beam and diffuse flux densities
    solar%nir_beam    = fansb * nirx
    solar%nir_diffuse = fAND * nirx
    if (solar%nir_beam <= zero) then
       solar%nir_beam    = zero
       solar%nir_diffuse = nirx
    end if
    if (nirx <= 0.1_wp) then
       solar%nir_beam    = 0.1_wp
       solar%nir_diffuse = 0.1_wp
    end if

  END SUBROUTINE diffuse_direct_radiation


  ! ------------------------------------------------------------------
  SUBROUTINE freq(lflai)
    ! THIS PROGRAM USES THE BETA DISTRIBUTION
    ! TO COMPUTE THE PROBABILITY FREQUENCY
    ! DISTRIBUTION FOR A KNOWN MEAN LEAF INCLINATION ANGLE
    ! STARTING FROM THE TOP OF THE CANOPY, WHERE llai=0

    ! AFTER GOEL AND STREBEL (1984)

    ! Note: there is a typo in equation [8] in Goel and Strebel (1984)
    ! correct form is:

    ! nuu = (1. - <theta>^2 / (90. * <theta>)) / (<theta^2> / <theta>^2 - 1.);

    ! theta : leaf angle
    ! <theta> : mean of leaf angle

    ! Canopy type MEAN THETA2 STD
    ! Planophile 26.8 1058.6 18.5
    ! Erectophile 63.2 4341.4 18.5
    ! Plagiophile 45.0 2289.7 16.3
    ! Extremophile 45.0 3110.4 32.9
    ! Uniform 45.0 2700.0 26.0
    ! Sperical 57.3 3747.6 21.5


    ! ! old
    ! if (lflai < 2.4)    !till here plagiophile
    ! MEAN = 45.;
    ! else
    ! MEAN = 78.99 - 14.84 * lflai + .244 * lflai * lflai;

    ! if(MEAN < 26)        !from here planophile
    ! MEAN=26;

    ! ! MEAN=70.;            ! mean leaf angle for eucalyptus from King 1997
    ! STD = 18;
    USE constants,  ONLY: one
    USE types,      ONLY: canopy
    USE parameters, ONLY: extra_nate
    USE utils,      ONLY: gammaf

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: lflai

    INTEGER, PARAMETER :: nc = 9
    REAL(wp) :: STD, MEAN, CONS
    REAL(wp) :: THETA2,nuu,SUM,MU,FL1,MU1,nu1
    REAL(wp) :: ANG, FL2, FL3, ztmp
    REAL(wp) :: leaf_angle_top, leaf_angle_bottom, leaf_angle_SD_top, leaf_angle_SD_bottom
    INTEGER(i4) :: i

    leaf_angle_top       = 55._wp    ! Mean leaf inclination angle at canopy top
    leaf_angle_bottom    = 25._wp    ! Mean leaf inclination angle at canopy bottom, old=25
    leaf_angle_SD_top    = 21.5_wp  ! Standard deviation of leaf inclination angle at canopy top
    leaf_angle_SD_bottom = 18._wp    ! Standard deviation of leaf inclination angle at canopy bottom, old =18

    MEAN = leaf_angle_top
    STD = leaf_angle_SD_top

    if (lflai < 2.4_wp) then
       MEAN = leaf_angle_top
       STD  = leaf_angle_SD_top
    end if
    if ((lflai >= 2.4_wp) .and. (lflai < 3.5_wp)) then
       MEAN = leaf_angle_top - (lflai-2.4_wp)/(3.5_wp-2.4_wp)*(leaf_angle_top-leaf_angle_bottom)
       STD = leaf_angle_SD_top - (lflai-2.4_wp)/(3.5_wp-2.4_wp)*(leaf_angle_SD_top-leaf_angle_SD_bottom)
    end if
    if (lflai >= 3.5_wp) then
       MEAN = leaf_angle_bottom
       STD = leaf_angle_SD_bottom
    end if

    if (extra_nate == 1) then
       ! !MC - v1
       ! MEAN = 57.3_wp
       ! STD = 21.5_wp
       ! !MC - v2
       ! MEAN = 45._wp
       ! STD = 18._wp
       !MC - smallest netrad error
       MEAN = 35._wp
       STD = 18._wp
    end if

    THETA2 = STD * STD + MEAN * MEAN ! = <theta^2> since Variance = STD*STD = <theta^2> - <theta>^2
    nuu    = (one - THETA2 / (90._wp * MEAN)) / (THETA2 / (MEAN * MEAN) - one)
    MU     = nuu * ((90._wp / MEAN) - one)
    SUM    = nuu + MU

    FL1 = GAMMAF(SUM) / (GAMMAF(nuu) * GAMMAF(MU))
    MU1 = MU - one
    nu1 = nuu - one

    CONS = one / real(nc,kind=wp)

    ! COMPUTE PROBABILITY DISTRIBUTION FOR 9 ANGLE CLASSES
    ! BETWEEN 5 AND 85 DEGREES, WITH INCREMENTS OF 10 DEGREES
    ztmp = one/90._wp
    do i=1, nc
       ANG = real(10*I-5,kind=wp)
       FL2 = (one-ANG*ztmp)**MU1
       FL3 = (ANG*ztmp)**nu1
       canopy%bdens(i) = CONS * FL1 * FL2 * FL3
    end do

  END SUBROUTINE freq


  ! ------------------------------------------------------------------
  SUBROUTINE gfunc()
    ! This subroutine computes the G function according to the
    ! algorithms of Lemeur(1973). This progrom computes G for a given
    ! sun angle. G changes with height due to change leaf angles
    USE setup,      ONLY: ncl, nsky, nl
    USE constants,  ONLY: zero, one, two, e2, pi, pi180
    USE types,      ONLY: solar, prof, canopy
    USE messages,   ONLY: message

    IMPLICIT NONE

    !INTEGER, PARAMETER :: nl  = 9 ! Lemeur defines bdens as delta F/(PI/N), WHERE HERE N=9
    REAL(wp), DIMENSION(nsky)  :: aden, pgg
    REAL(wp), DIMENSION(nsky+1) :: TT, sin_TT, del_TT,del_sin
    REAL(wp) :: PPP, PP, aang
    REAL(wp) :: cos_A, cos_B, sin_A, sin_B, X, Y, sin_TT0, sin_TT1
    REAL(wp) :: T0, TII, TT0, TT1
    REAL(wp) :: R, S, square
    REAL(wp) :: llai
    INTEGER(i4) :: i, j, k, kp1

    llai   = zero
    square = zero
    ! Midpoint of azimuthal intervals
    aden(1:nsky) = 0.0625_wp
    do j=1, nsky+1
       k = 2*j-3
       TT(j) = pi / real(nsky,kind=wp) * real(k,kind=wp)
       sin_TT(j) = sin(TT(j))
    end do
    del_TT(1:nsky)  = TT(2:nsky+1) - TT(1:nsky)
    del_sin(1:nsky) = sin_TT(2:nsky+1) - sin_TT(1:nsky)
    ! Compute the G function for each layer
    do j=ncl, 1, -1
       ! need LAI from LAIZ
       llai  = llai  + prof%dLAIdz(j)
       ! Calculate the leaf angle probabilty distribution, bdens
       call freq(llai)
       ! Lemeur defines bdens as delta F/(PI/N), WHERE HERE N=9
       PPP = zero
       do i=1, nl
          aang = real((i-1)*10+5,kind=wp) * pi180 ! 5, 15, ..., 85 degree
          cos_A = cos(aang)
          cos_B = cos(solar%beta_rad)
          sin_A = sin(aang)
          sin_B = sin(solar%beta_rad)
          X = cos_A * sin_B
          Y = sin_A * cos_B
          if ((aang-solar%beta_rad) <= zero) then
             do k=1, nsky
                pgg(k) = X * del_TT(k) + Y * del_sin(k)
             end do
          else
             T0  = (one + X/Y)
             TII = (one - X/Y)
             if (T0/TII > zero) then
                square = sqrt(T0/TII)
             else
                call message('GFUNC: ','bad T0/TII')
             end if
             TT0 = two * atan(square)
             TT1 = two * pi - TT0
             sin_TT0 = sin(TT0)
             sin_TT1 = sin(TT1)
             do k=1, nsky
                kp1 = k+1
                if ((TT(kp1) - TT0) <= zero) then
                   pgg(k) = X * del_TT(k) + Y *del_sin(k)
                else
                   if ((TT(kp1) - TT1) <= zero) then
                      if ((TT0 - TT(k)) <= zero) then
                         pgg(k) = -X * del_TT(k) - Y * del_sin(k)
                      else
                         R = X * (TT0 - TT(k)) + Y * (sin_TT0 - sin_TT(k))
                         S = X * (TT(kp1) - TT0) + Y * (sin_TT(kp1) - sin_TT0)
                         pgg(k) = R - S
                      end if
                   else
                      if ((TT1 - TT(k)) <= zero) then
                         pgg(k) = X * del_TT(k) + Y * del_sin(k)
                      else
                         R = X * (TT1 - TT(k)) + Y * (sin_TT1 - sin_TT(k))
                         S = X * (TT(kp1) - TT1) + Y * (sin_TT(kp1) - sin_TT1)
                         pgg(k) = S - R
                      end if
                   end if
                end if
             end do
          end if
          ! Compute the integrated leaf orientation function
          PP = sum(pgg(1:nsky) * aden(1:nsky))
          PPP  = PPP  + (PP * canopy%bdens(i) * (real(nl,kind=wp)/pi))
       end do
       prof%Gfunc_solar(j) = PPP
       if (prof%Gfunc_solar(j) < e2)  prof%Gfunc_solar(j) = e2
       if (prof%Gfunc_solar(j) > one) prof%Gfunc_solar(j) = one
    end do

  END SUBROUTINE gfunc


  ! ------------------------------------------------------------------
  SUBROUTINE gfunc_diffuse()
    ! This subroutine computes the G Function according to
    ! the algorithms of Lemeur (1973, Agric. Meteorol. 12: 229-247).
    ! The original code was in Fortran and was converted to C
    ! This program computs G for each sky sector, as
    ! needed to compute the transmission of diffuse light
    ! G varies with height to account for vertical variations
    ! in leaf angles
    USE setup,     ONLY: ncl, nsky, nl
    USE constants, ONLY: zero, one, two, pi180, pi, pi2, pi9
    USE types,     ONLY: prof, canopy
    USE messages,  ONLY: message

    IMPLICIT NONE

    !INTEGER, PARAMETER :: nl  = 9 ! Lemeur defines bdens as delta F/(PI/N), WHERE HERE N=9
    REAL(wp), DIMENSION(nsky)  :: aden, PGF
    REAL(wp), DIMENSION(nsky+1) :: TT, sin_TT(18), del_TT(18), del_sin(18)
    REAL(wp) :: PPP
    REAL(wp) :: ang, dang, aang
    REAL(wp) :: cos_A, cos_B, sin_A, sin_B, X, Y
    REAL(wp) :: T0, TII, TT0, TT1
    REAL(wp) :: R,S, PP, bang
    REAL(wp) :: sin_TT1, sin_TT0, square
    REAL(wp) :: llai, ztmp
    INTEGER(i4) :: i, j, k, l, m, mp1

    square = zero
    llai   = zero
    ang    = 5._wp*pi180
    dang   = two*ang
    ! Midpoints for azimuth intervals
    aden(1:nsky) = 0.0625_wp
    ztmp = pi/real(nsky,kind=wp)
    ! 3.1415/16 =0.1963
    do j=1, nsky+1
       k = 2*j-3
       !TODO: TT(j)    = ztmp*real(k,kind=wp)
       TT(j)    = 0.1963_wp*real(k,kind=wp)
       sin_TT(j) = sin(TT(j))
    end do
    del_TT(1:nsky)  = TT(2:nsky+1) - TT(1:nsky)
    del_sin(1:nsky) = sin_TT(2:nsky+1) - sin_TT(1:nsky)

    do l=1, nl
       bang = ang
       ! COMPUTE G FUNCTION FOR EACH LAYER IN THE CANOPY
       llai = zero
       do k=ncl, 1, -1
          ! top of the canopy is ncl. Its cumulative LAI starts with 0.
          llai = llai + prof%dLAIdz(k)
          ! CALCULATE PROBABILITY FREQUENCY DISTRIBUTION, BDENS
          call freq(llai)
          ! LEMEUR DEFINES BDENS AS DELTA F/(PI/N), WHERE HERE N=9
          PPP = zero
          do i=1, nl
             aang = real((i-1)*10+5,kind=wp) * pi180
             cos_B = cos(bang)
             sin_B = sin(bang)
             cos_A = cos(aang)
             sin_A = sin(aang)
             X = cos_A * sin_B
             Y = sin_A * cos_B
             if ((aang - bang) <= zero) then
                do m=1, nsky
                   PGF(m) = X * del_TT(m) + Y * del_sin(m)
                end do
             else
                T0 = one + X /Y
                TII = one - X / Y
                if (T0/TII > 0) then
                   square = sqrt(T0/TII)
                else
                   call message('GFUNC_DIFFUSE: ','bad sqrt in ggfun')
                end if
                TT0     = two * atan(square)
                sin_TT0 = sin(TT0)
                TT1     = pi2 - TT0
                sin_TT1 = sin(TT1)
                do m=1, nsky
                   mp1 = m+1
                   if ((TT(mp1) - TT0) <= zero) then
                      PGF(m) = X * del_TT(m) + Y * del_sin(m)
                   else
                      if ((TT(mp1) - TT1) <= zero) then
                         if ((TT0 - TT(m)) <= zero) then
                            PGF(m) = -X * del_TT(m) - Y * del_sin(m)
                         else
                            R = X * (TT0 - TT(m)) + Y * (sin_TT0 - sin_TT(m))
                            S = X * (TT(mp1) - TT0) + Y * (sin_TT(mp1) - sin_TT0)
                            PGF(m) = R - S
                         end if
                      else
                         if ((TT1 - TT(m)) <= zero) then
                            PGF(m) = X * del_TT(m) + Y * del_sin(m)
                         else
                            R = X * (TT1 - TT(m)) + Y * (sin_TT1 - sin_TT(m))
                            S = X * (TT(mp1) - TT1) + Y * (sin_TT(mp1) - sin_TT1)
                            PGF(m) = S - R
                         end if
                      end if
                   end if
                end do
             end if
             ! Compute the integrated leaf orientation function, Gfun
             PP = zero
             do m=1, nsky
                PP = PP + (PGF(m) * aden(m))
             end do
             PPP = PPP + (PP * canopy%bdens(i) * pi9)
          end do
          prof%Gfunc_sky(k,l) = PPP
       end do
       ang = ang + dang
    end do

  END SUBROUTINE gfunc_diffuse


  ! ------------------------------------------------------------------
  SUBROUTINE irflux()
    ! This subroutine is adapted from:
    ! Norman, J.M. 1979. Modeling the complete crop canopy.
    ! Modification of the Aerial Environment of Crops.
    ! B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.

    ! Compute probability of penetration for diffuse radiation for each layer in the canopy.
    ! IR radiation is isotropic.
    USE types,      ONLY: solar, met, prof, soil
    USE constants,  ONLY: nerr, one, e12, TN0, sigma
    USE setup,      ONLY: ncl
    USE parameters, ONLY: ep, epsigma, epsoil
    USE utils,      ONLY: corr

    IMPLICIT NONE

    INTEGER  :: j
    REAL(wp) :: reflc_lay_IR
    REAL(wp), DIMENSION(1:ncl)  :: Tk_sun_filt, Tk_shade_filt, IR_source_sun, IR_source_shade, IR_source
    REAL(wp), DIMENSION(1:ncl+1) :: SUP, ir_dn, ir_up
    REAL(wp), DIMENSION(1:ncl)  :: SDN
    REAL(wp) :: emiss_IR_soil
    REAL(wp) :: r2, r22
    INTEGER  :: i_check

    ! IR down flux at top of canopy
    solar%ir_dn(ncl+1) = sky_ir(met%T_Kelvin)

    ! Integrated probability of diffuse sky radiation penetration
    ! EXPDIF(JJ) is computed in RAD

    ! compute IR radiative source flux as a function of
    ! leaf temperature weighted according to sunlit and shaded fractions

    ! source=ep*sigma*(laisun*tksun^4 + laish*tksh^4)
    ! remember energy balance is done on layers not levels.
    ! so level ncl+1 must use tl from layer ncl
    Tk_sun_filt(:)     = prof%sun_tleaf_filter(:)+TN0
    Tk_shade_filt(:)   = prof%shd_tleaf_filter(:)+TN0
    IR_source_sun(:)   = solar%prob_beam(:) * Tk_sun_filt(:)*Tk_sun_filt(:)*Tk_sun_filt(:)*Tk_sun_filt(:)
    IR_source_shade(:) = solar%prob_shd(:)  * Tk_shade_filt(:)*Tk_shade_filt(:)*Tk_shade_filt(:)*Tk_shade_filt(:)
    IR_source(:)       = epsigma * (IR_source_sun(:) + IR_source_shade(:))
    ! Intercepted IR that is radiated up
    SUP(2:ncl+1)       = IR_source(1:ncl) * (one-solar%exxpdir(1:ncl))
    ! Intercepted IR that is radiated downward
    SDN(:)             = IR_source(:) * (one-solar%exxpdir(:))

    ! First guess of up and down values

    ! Downward IR radiation, sum of that from upper layer that is transmitted
    !   and the downward source generated in the upper layer.
    do j=ncl, 1, -1
       solar%ir_dn(j) = solar%exxpdir(j) * solar%ir_dn(j+1) + SDN(j)
    end do
    ! same for upward
    do j=2, ncl+1
       solar%ir_up(j) = solar%exxpdir(j-1) * solar%ir_up(j-1) + SUP(j)
    end do
    ! ground emission
    emiss_IR_soil  = epsoil*sigma &
         * (soil%tsrf_filter+TN0)*(soil%tsrf_filter+TN0) &
         * (soil%tsrf_filter+TN0)*(soil%tsrf_filter+TN0)
    SUP(1)         = solar%ir_dn(1) * (one-epsoil)
    solar%ir_up(1) = emiss_IR_soil + SUP(1)

    i_check = 0
    r2      = one
    r22     = one
    do while (r2 > e12 .or. r22 > e12)
       ! save result
       ir_up(1:ncl+1) = solar%ir_up(1:ncl+1)
       ir_dn(1:ncl+1) = solar%ir_dn(1:ncl+1)
       ! go down
       do j=ncl, 1, -1
          reflc_lay_IR   = (one-solar%exxpdir(j)) * (one-ep)
          solar%ir_dn(j) = solar%exxpdir(j) * solar%ir_dn(j+1) + solar%ir_up(j) * reflc_lay_IR + SDN(j)
       end do
       ! go up
       SUP(1)         = solar%ir_dn(1) * (one-epsoil)
       solar%ir_up(1) = emiss_IR_soil + SUP(1)
       do j=2, ncl+1
          reflc_lay_IR   = (one-solar%exxpdir(j-1)) * (one-ep)
          solar%ir_up(j) = reflc_lay_IR * solar%ir_dn(j) + solar%ir_up(j-1) * solar%exxpdir(j-1) + SUP(j)
       end do
       ! check if stabilised
       r2  = one - corr(ir_up(1:ncl+1),solar%ir_up(1:ncl+1))**2
       r22 = one - corr(ir_dn(1:ncl+1),solar%ir_dn(1:ncl+1))**2
       i_check = i_check+1
       if (i_check > 10) write(nerr,*) "Bad irflux"
    end do

  END SUBROUTINE irflux


  SUBROUTINE Beta_distr(HHH,FQ,BF,III,NNN)

    ! add individual B_Beta_distr subroutine to calculate beta distribution for WAI and LAI distribution
    ! Yuan 2018.04.23
        ! Beta distribution
    ! f(x) = x^(p-1) (1-x)^(q-1) / B(v,w)
    ! B(v,w) = int from 0 to 1 x^(p-1) (1-x)^(q-1) dx
    ! p = mean{(mean(1-mean)/var)-1}
    ! q =(1-mean){(mean(1-mean)/var)-1}
    ! ENTER THE HEIGHT AT THE MIDPOINT OF THE LAYER
    ! Normalize height
    USE constants,  ONLY: zero, half, one, two
    USE setup,      ONLY: ncl
    USE parameters, ONLY: ht

    IMPLICIT NONE
    INTEGER(i4) :: i, NNN! number of elements in midpoint array
    REAL(wp) :: ztmp
    REAL(wp) :: TF, MU1, MU2, III ! III is beta integer
    REAL(wp) :: dx, DX2, DX4, X, P_beta, Q_beta, F1, F2, F3
    REAL(wp), DIMENSION(1:NNN)   :: HHH, Htmp, FQ ! midpt height and midpt frequency
    REAL(wp), DIMENSION(1:ncl)     :: BF ! beta function
!    print *, HHH
!    print *, FQ
    ztmp = one / ht
!    print *, "HHH before: "
!    print *, HHH
    Htmp(1:NNN) = HHH(1:NNN) * ztmp
!    print *, "HHH after: "
!    print *, Htmp
    ! TOTAL F IN EACH LAYER, TF SHOULD SUM TO EQUAL lai
    TF  = sum(FQ(1:NNN))
    ! weighted mean lai
    MU1 = sum(Htmp(1:NNN)*FQ(1:NNN))
    ! weighted variance
    MU2 = sum(Htmp(1:NNN)*Htmp(1:NNN)*FQ(1:NNN))
    ! normalize mu by lai
    MU1 = MU1 / TF
    MU2 = MU2 / TF
    ! compute Beta parameters
    P_beta = MU1 * (MU1 - MU2) / (MU2 - MU1 * MU1)
    Q_beta = (one - MU1) * (MU1 - MU2) / (MU2 - MU1 * MU1)
    P_beta = P_beta - one
    Q_beta = Q_beta - one
    ! integrate Beta function, with Simpson''s Approx.
    ! The boundary conditions are level 1 is height of ground
    ! and level ncl+1 is height of canopy. Layer 1 is between
    ! height levels 1 and 2. Layer ncl is between levels
    ! ncl and ncl+1
    ! Thickness of layer
    dx  = one / ncl
    DX2 = half * dx
    DX4 = dx / 4._wp
    X   = DX4
    F2  = X**P_beta * (one-X)**Q_beta
    X   = X + DX4
    F3  = X**P_beta * (one-X)**Q_beta
    ! start integration at lowest boundary
    BF(1) = DX4 * (4._wp * F2 + F3) / 3._wp
    III = BF(1)
    do i=2, ncl-1
       F1 = F3
       X  = X + DX2
       F2 = X**P_beta * (one-X)**Q_beta
       X  = X + DX2
       F3 = X**P_beta * (one-X)**Q_beta
       BF(i) = DX2 * (F1 + 4._wp * F2 + F3) / 3._wp
       III = III + BF(i)
    end do
    F1 = F3
    X  = X + DX4
    F2 = X**P_beta * (one-X)**Q_beta
    ! compute integrand at highest boundary
    BF(ncl) = DX4 * (F1 + 4._wp * F2) / 3._wp
    III   = III + BF(ncl)

!    print *, HHH
!    print *, FQ
!    print *, (BF)/III

  END SUBROUTINE Beta_distr


  SUBROUTINE wai_distr(wai,B_func,wai_z)
  ! add WAI distribution Yuan 2018.03.15:

  ! 1. seperate PAI into WAI and LAI: PAI = WAI + LAI
  ! 2. WAI = stem + branch, the fraction can be set mannually.
  ! 3. WAI% of stem layer = volume% layer to the total stem
  ! 4. volume% layer is a function of diameter per layer (D).
  ! 5. diameter per layer is up to the D/DBH, which is a function of total tree height and the layer height
  ! 6. WAI% of branch layer follows beta distribution, but indipendent from beta_lai
    USE setup,      ONLY: ncl
    USE parameters, ONLY: delz
    USE constants,  ONLY: zero, pi
    USE types, ONLY: debug
    IMPLICIT NONE

    REAL(wp), DIMENSION(1:ncl)   :: B_func
    REAL(wp), DIMENSION(1:ncl)   :: vol_perc
    REAL(wp), DIMENSION(1:ncl)   :: wai_z
    REAL(wp) :: h_stem, stem_wai, branch_wai, wai, htmp
    REAL(wp) :: d_frac! diameter fraction against DBH
    REAL(wp) :: p0, p1, p2, p3
    INTEGER(i4) :: j, j_branch ! j_stem = layer on which stem ends and branch extension starts

  h_stem = 17 ! stem height = 17m
  stem_wai = 0.2 ! 0.2 and 0.8 ratio, from Yuan's measurement. 2018.04.23
  branch_wai = 0.8
  p0 = 102
  p1 = -2.6
  p2 = 0.08
  p3 = -0.0023

  !QD <- 0.55 ! quadratic mean diameter
  j_branch = nint(h_stem/delz,i4)-1 ! branch starts at layer 17
  ! calculate stem volume fraction on each layer
  vol_perc = zero
   ! lower layers, stems
  do j=1, ncl

    htmp = j*delz
    ! diameter function with height
    d_frac =p0+p1*(htmp)+p2*(htmp**2)+p3*(htmp**3)
    vol_perc(j) = d_frac**3
  end do
  ! volume frac per layer stem
  vol_perc =  vol_perc/sum(vol_perc)
!  print *, sum(vol_perc)
 ! print *, j_branch
  ! wai distr of the stem based on vol%
  j_branch = 20 ! j_branch = 18 from Yuan's measurements 2018.04.23
  wai_z(1:j_branch) = stem_wai*wai*vol_perc(1:j_branch)/sum(vol_perc(1:j_branch))
  ! wai of branch is scaled with lai
  wai_z((j_branch+1):ncl)= branch_wai*wai*B_func((j_branch+1):ncl)/sum(B_func((j_branch+1):ncl))
!  debug%R3 =sum(B_func((j_branch+1):ncl))
!  wai_z = wai*B_func/sum(B_func)
  END SUBROUTINE wai_distr

  ! ------------------------------------------------------------------
  SUBROUTINE lai_time()
    ! Evaluate how LAI and other canopy structural variables vary with time
    ! ATTENTION: We separate original LAI into PAI - WAI; old LAI = new PAI, old PAI = new WAI; New PAI = New LAI + New WAI
    ! new lai for photosyn, new PAI for radiation, new WAI for bole resp
    USE types,      ONLY: time, solar, prof, input, debug, iswitch
    USE constants,  ONLY: zero, half, one, two, nbetamax
    USE setup,      ONLY: ncl, nbeta, nbeta_w
    USE parameters, ONLY: delz, extra_nate, ht, lai, pai, wai, markov, nup, ht_midpt, lai_freq, pai_freq, wai_freq, &
         wai_midpt, par_reflect, par_trans, par_soil_refl_dry, nir_reflect, nir_trans, nir_soil_refl_dry

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:ncl)   :: lai_z
    REAL(wp), DIMENSION(1:ncl)   :: wai_z ! add wai distribution Yuan 2018.03.04
    REAL(wp), DIMENSION(1:ncl)   :: beta_fnc, beta_fnc_w
    REAL(wp), DIMENSION(1:nbetamax) :: lai_freq_local
    REAL(wp), DIMENSION(1:nbetamax) :: ht_midpt_local
    REAL(wp), DIMENSION(1:ncl)   :: dff, XX, AA, DA, dff_Markov
    REAL(wp), DIMENSION(1:ncl)   :: cos_AA, sin_AA, exp_diffuse
!    REAL(wp) :: TF, MU1, MU2,
    REAL(wp) :: integr_beta, integr_beta_w
!    REAL(wp) :: dx, DX2, DX4, X, P_beta, Q_beta, F1, F2, F3
    REAL(wp) :: cum_lai, sumlai, cum_ht! sumpai Yuan 2018.03.02
    INTEGER(i4) :: i ! number of elements in wai midpoint array. Yuan 2018.04.25

    ! seasonal update of model parameters
    time%wai = wai ! Yuan 2018.03.02
    ! Winter 1
    if (time%days < time%leafout) then
        if (iswitch%wai_new==0) then
            time%lai                = pai
            else
            time%lai                = zero!time%wai!old pai = new wai
        end if


       ! PAR wave band: after Norman (1979) and NASA report
       solar%par_reflect       = par_reflect(1)
       solar%par_trans         = par_trans(1)
       solar%par_soil_refl_dry = par_soil_refl_dry(1)
       solar%par_absorbed      = (one - solar%par_reflect - solar%par_trans)
       ! NIR wave band: after Norman (1979) and NASA report
       solar%nir_reflect       = nir_reflect(1)
       solar%nir_trans         = nir_trans(1)
       solar%nir_soil_refl_dry = nir_soil_refl_dry(1)
       solar%nir_absorbed      = (one - solar%nir_reflect - solar%nir_trans)
    end if
    ! Spring
    if (time%days >= time%leafout .and. time%days < time%leaffull) then
        if (iswitch%wai_new==0) then
            time%lai                = pai + (time%days - time%leafout) * (lai - pai) / (time%leaffull-time%leafout) ! old LAI
        else
                  ! time%lai                = wai + (time%days - time%leafout) * (lai-wai) / (time%leaffull-time%leafout) !new LAI
            time%lai                = (time%days - time%leafout) * (lai) / (time%leaffull-time%leafout) !new LAI
        end if
      ! PAR wave band: after Norman (1979) and NASA report
       solar%par_reflect       = par_reflect(2)
       solar%par_trans         = par_trans(2)
       solar%par_soil_refl_dry = par_soil_refl_dry(2)
       solar%par_absorbed      = (one - solar%par_reflect - solar%par_trans)
       ! NIR wave band: after Norman (1979) and NASA report
       solar%nir_reflect       = nir_reflect(2)
       solar%nir_trans         = nir_trans(2)
       solar%nir_soil_refl_dry = nir_soil_refl_dry(2)
       solar%nir_absorbed      = (one - solar%nir_reflect - solar%nir_trans)
    end if
    ! Summer
    if (time%days >= time%leaffull .and. time%days < time%leaffall) then
       time%lai                = lai ! NEW lai
 !      time%pai                = time%lai + time%wai ! both leaf and wood
       ! PAR wave band: after Norman (1979) and NASA report
       solar%par_reflect       = par_reflect(3)
       solar%par_trans         = par_trans(3)
       solar%par_soil_refl_dry = par_soil_refl_dry(3)
       solar%par_absorbed      = (one - solar%par_reflect - solar%par_trans)
       ! NIR wave band: after Norman (1979) and NASA report
       solar%nir_reflect       = nir_reflect(3)
       solar%nir_trans         = nir_trans(3)
       solar%nir_soil_refl_dry = nir_soil_refl_dry(3)
       solar%nir_absorbed      = (one - solar%nir_reflect - solar%nir_trans)
    end if
    ! Fall
    if (time%days >= time%leaffall .and. time%days < time%leaffallcomplete) then
        if (iswitch%wai_new==0) then
                  time%lai   = lai - (time%days - time%leaffall) * (lai-pai) / (time%leaffallcomplete-time%leaffall)
       else
                ! new lai , new wai
!      time%lai                = lai - (time%days - time%leaffall) * (lai-wai) / (time%leaffallcomplete-time%leaffall) ! Yuan 2018.03.01
       time%lai                = lai - (time%days - time%leaffall) * (lai) / (time%leaffallcomplete-time%leaffall) ! Yuan 2018.03.02
       !time%pai                = time%lai + time%wai

        end if
       ! PAR wave band: after Norman (1979) and NASA report
       solar%par_reflect       = par_reflect(4)
       solar%par_trans         = par_trans(4)
       solar%par_soil_refl_dry = par_soil_refl_dry(4)
       solar%par_absorbed      = (one - solar%par_reflect - solar%par_trans)
       ! NIR wave band: after Norman (1979) and NASA report
       solar%nir_reflect       = nir_reflect(4)
       solar%nir_trans         = nir_trans(4)
       solar%nir_soil_refl_dry = nir_soil_refl_dry(4)
       solar%nir_absorbed      = (one - solar%nir_reflect - solar%nir_trans)
    end if
    ! Winter 2
    if (time%days >= time%leaffallcomplete) then
        if (iswitch%wai_new==0) then
         time%lai                = pai
            else
       ! time%lai                = wai ! YUAN 2018.03.01
       time%lai                = zero ! YUAN 2018.03.02
        end if



       ! PAR wave band: after Norman (1979) and NASA report
       solar%par_reflect       = par_reflect(1)
       solar%par_trans         = par_trans(1)
       solar%par_soil_refl_dry = par_soil_refl_dry(1)
       solar%par_absorbed      = (one - solar%par_reflect - solar%par_trans)
       ! NIR wave band: after Norman (1979) and NASA report
       solar%nir_reflect       = nir_reflect(1)
       solar%nir_trans         = nir_trans(1)
       solar%nir_soil_refl_dry = nir_soil_refl_dry(1)
       solar%nir_absorbed      = (one - solar%nir_reflect - solar%nir_trans)
    end if
    ! for Nate McDowell''s juniper site read LAI instead of diffuse PAR
    if (extra_nate == 1) time%lai = max(input%lai_up+input%lai_down, pai)
    ! height of mid point of layer scaled to height of forest
    time%pai                = time%lai + time%wai
    ht_midpt_local(1:nbeta) = ht_midpt(1:nbeta)
    ! use LAI distribution
!    if (time%days >= time%leafout .and. time%days <= time%leaffallcomplete) then
       ! lai of the layers at the midpoint of height
!       lai_freq_local(1:nbeta) = lai_freq(1:nbeta)*lai
       lai_freq_local(1:nbeta) = lai_freq(1:nbeta)*lai
    ! get beta function for LAI:
    call Beta_distr(ht_midpt_local,lai_freq_local,beta_fnc,integr_beta,nbeta)
    lai_z(1)       = beta_fnc(1) * time%lai / integr_beta
    lai_z(1:ncl-1) = beta_fnc(1:ncl-1) * (time%lai / integr_beta)
    lai_z(ncl)     = beta_fnc(ncl) * time%lai / integr_beta
        ! re-index layers of lai_z.
    ! layer 1 is between ground and 1st level
    ! layer jtot is between level ncl and top of canopy (ncl+1)
    cum_ht  = ncl*delz
    cum_lai = sum(lai_z(1:ncl))
    ! use prof%dLAIdz for radiative transfer model
    prof%dLAIdz(1:ncl) = lai_z(1:ncl)
    ! Beta distribution
    ! f(x) = x^(p-1) (1-x)^(q-1) / B(v,w)
    ! B(v,w) = int from 0 to 1 x^(p-1) (1-x)^(q-1) dx
    ! p = mean{(mean(1-mean)/var)-1}
    ! q =(1-mean){(mean(1-mean)/var)-1}
    ! ENTER THE HEIGHT AT THE MIDPOINT OF THE LAYER
    ! Normalize height
!    ztmp = one / ht
!    ht_midpt_local(1:nbeta) = ht_midpt_local(1:nbeta) * ztmp
!    ! TOTAL F IN EACH LAYER, TF SHOULD SUM TO EQUAL lai
!    TF  = sum(lai_freq_local(1:nbeta))
!    ! weighted mean lai
!    MU1 = sum(ht_midpt_local(1:nbeta)*lai_freq_local(1:nbeta))
!    ! weighted variance
!    MU2 = sum(ht_midpt_local(1:nbeta)*ht_midpt_local(1:nbeta)*lai_freq_local(1:nbeta))
!    ! normalize mu by lai
!    MU1 = MU1 / TF
!    MU2 = MU2 / TF
!    ! compute Beta parameters
!    P_beta = MU1 * (MU1 - MU2) / (MU2 - MU1 * MU1)
!    Q_beta = (one - MU1) * (MU1 - MU2) / (MU2 - MU1 * MU1)
!    P_beta = P_beta - one
!    Q_beta = Q_beta - one
!    ! integrate Beta function, with Simpson''s Approx.
!    ! The boundary conditions are level 1 is height of ground
!    ! and level ncl+1 is height of canopy. Layer 1 is between
!    ! height levels 1 and 2. Layer ncl is between levels
!    ! ncl and ncl+1
!    ! Thickness of layer
!    dx  = one / ncl
!    DX2 = half * dx
!    DX4 = dx / 4._wp
!    X   = DX4
!    F2  = X**P_beta * (one-X)**Q_beta
!    X   = X + DX4
!    F3  = X**P_beta * (one-X)**Q_beta
!    ! start integration at lowest boundary
!    beta_fnc(1) = DX4 * (4._wp * F2 + F3) / 3._wp
!    integr_beta = beta_fnc(1)
!    do i=2, ncl-1
!       F1 = F3
!       X  = X + DX2
!       F2 = X**P_beta * (one-X)**Q_beta
!       X  = X + DX2
!       F3 = X**P_beta * (one-X)**Q_beta
!       beta_fnc(i) = DX2 * (F1 + 4._wp * F2 + F3) / 3._wp
!       integr_beta = integr_beta + beta_fnc(i)
!    end do
!    F1 = F3
!    X  = X + DX4
!    F2 = X**P_beta * (one-X)**Q_beta
!    ! compute integrand at highest boundary
!    beta_fnc(ncl) = DX4 * (F1 + 4._wp * F2) / 3._wp
!    integr_beta   = integr_beta + beta_fnc(ncl)
    ! lai_z IS THE LEAF AREA AS A FUNCTION OF Z
    !
    ! beta_fnc is the pdf for the interval dx


    !!add wai distribution, Yuan 2018.03.04
! if to switch on Yuan's wai distribution
     if (iswitch%wai_new==0) then
       prof%dWAIdz(1:ncl) = prof%dLAIdz(1:ncl)*(time%wai/time%lai)
       else
    call Beta_distr(wai_midpt,wai_freq,beta_fnc_w,integr_beta_w,nbeta_w)
!    print *, "wai_mdpt="
!    print *, wai_midpt
!    print *, "wai_freq="
!    print *, wai_freq
!    print *, "beta_fnc_w="
!    print *, beta_fnc_w
!    print *, "int_beta_w="
!    print *, integr_beta_w
!    print *, "nbeta_w="
!    print *, nbeta_w

 !   wai_z(1)       = beta_fnc(1) * time%wai / integr_beta
 !   wai_z(1:ncl-1) = beta_fnc(1:ncl-1) * (time%wai / integr_beta)
 !   wai_z(ncl)     = beta_fnc(ncl) * time%wai / integr_beta
 !   wai_z (1:ncl)        = time%wai/ncl

    call wai_distr(time%wai,beta_fnc_w,wai_z)
        prof%dWAIdz(1:ncl) = wai_z(1:ncl)
     end if

    ! get beta function for WAI:

!    print *, "wai="
!    print *, wai_z




!debug%D1 = wai_z

!debug%D2 = beta_fnc_w

!debug%D3 = wai_midpt

!debug%D4 = wai_freq

!debug%R1 = time%wai

!debug%R2 = integr_beta_w

!debug%in1 = nbeta_w

!print *, prof%dLAIdz
    ! for Nate McDowell''s juniper site, hardcode lai distribution for 40 layers
    if (extra_nate == 1) then
       ! equall distribution in lowest 50cm = 8 layers -> nup=9
       prof%dLAIdz(1:nup-1) = max(input%lai_down/real(nup-1,kind=wp), pai/real(ncl,kind=wp))
       ! pear shaped profile of upper canopy
       ! start tapering from ca. 2m = layer 31
       prof%dLAIdz(nup:ncl-nup-1)  = one
       do i=ncl-nup, ncl
          prof%dLAIdz(i) = one - one/real(nup+1,kind=wp)*real(nup-(ncl-i),kind=wp)
       end do
       ! normalise pear shape
       sumlai = sum(prof%dLAIdz(nup:ncl))
       prof%dLAIdz(nup:ncl) = max(prof%dLAIdz(nup:ncl)/sumlai*input%lai_up, pai/real(ncl,kind=wp))
       ! sync time%lai with dLAIdz
       time%lai = sum(prof%dLAIdz(1:ncl))
    else
       ! normalise total LAI including layers with PAI
       sumlai = sum(prof%dLAIdz(1:ncl)) ! only leaves
!       prof%dLAIdz(1:ncl) = prof%dLAIdz(1:ncl)*(time%lai/sumlai)! sumlai=0
    end if
    ! PAI
    if (iswitch%wai_new==0) then
       prof%dPAIdz(1:ncl) = prof%dLAIdz(1:ncl)*(pai/time%lai)
       else
            !!own distribution!!!! Yuan 2018.03.01
     prof%dPAIdz(1:ncl) = prof%dWAIdz(1:ncl) + prof%dLAIdz(1:ncl)! use PAI for radiation
    end if

 !   print *, prof%dPAIdz(1:ncl)
 !   print *, prof%dLAIdz(1:ncl)


    ! Direction cosine for the normal between the mean
    call gfunc_diffuse()
    ! compute the probability of diffuse radiation penetration through the
    ! hemisphere. This computation is not affected by penubra
    ! since we are dealing only with diffuse radiation from a sky
    ! sector.
    !
    ! The probability of beam penetration is computed with a
    ! Markov distribution.
    ! both use new PAI (leaf and wood) instead of LAI only Yuan 2018.03.01
    if (iswitch%wai_new==0) then
        dff(1:ncl) = prof%dLAIdz(1:ncl) !+ prof%dPAIdz(1:ncl)
        else
        dff(1:ncl) = prof%dPAIdz(1:ncl)
    end if


    XX(1:ncl) = zero
    AA(1:ncl) = 0.087_wp
    DA(1:ncl) = 0.1745_wp
    ! The leaf clumping coefficient. From Chason et al. 1990 and studies on WBW
    dff_Markov(1:ncl) = dff(1:ncl)*markov
    do i=1, 9
       cos_AA(1:ncl) = cos(AA(1:ncl))
       sin_AA(1:ncl) = sin(AA(1:ncl))
       ! probability of photon transfer through a sky section
       ! for clumped foliage and Markov model
       ! for spherical distribution
       ! exp_diffuse = exp(-DFF * prof%Gfunc_sky(J, II) / cos_AA)
       exp_diffuse(1:ncl) = exp(-dff_Markov(1:ncl) * prof%Gfunc_sky(1:ncl,i) / cos_AA(1:ncl))
       XX(1:ncl) = XX(1:ncl) + (cos_AA(1:ncl) * sin_AA(1:ncl) * exp_diffuse(1:ncl))
       AA(1:ncl) = AA(1:ncl) + DA(1:ncl)
    end do
    ! Itegrated probability of diffuse sky radiation penetration for each layer
    solar%exxpdir(1:ncl) = two * XX(1:ncl) * DA(1:ncl)
    where (solar%exxpdir(1:ncl) > one) solar%exxpdir(1:ncl) = one

  END SUBROUTINE lai_time


  ! ------------------------------------------------------------------
  SUBROUTINE nir()
    ! This subroutine computes the flux density of direct and diffuse
    ! radiation in the near infrared waveband. The Markov model is used
    ! to compute the probability of beam penetration through clumped foliage.

    ! The algorithms of Norman (1979) are used.

    ! Norman, J.M. 1979. Modeling the complete crop canopy.
    ! Modification of the Aerial Environment of Crops.
    ! B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
    USE setup,      ONLY: ncl
    USE constants,  ONLY: zero, one, e1, e2, isnight
    USE types,      ONLY: solar, prof, iswitch, input
    USE parameters, ONLY: markov

    IMPLICIT NONE

    INTEGER(i4) :: j, jp1, jm1, ITER, IREP
    REAL(wp), DIMENSION(ncl)   :: SDN, transmission_layer, reflectance_layer
    REAL(wp), DIMENSION(ncl)   :: nir_normal, NSUNEN
    REAL(wp), DIMENSION(ncl+1) :: SUP, TBEAM, ADUM, beam
    REAL(wp) :: fraction_beam
    REAL(wp) :: exp_direct, dff
    REAL(wp) :: TLAY2
    REAL(wp) :: DOWN, UP, ztmp

    solar%nir_total     = solar%nir_beam + solar%nir_diffuse
    fraction_beam       = solar%nir_beam / solar%nir_total
    beam(ncl+1)         = fraction_beam
    TBEAM(ncl+1)        = fraction_beam
    solar%nir_dn(ncl+1) = one - fraction_beam

    if (solar%nir_total > one .and. solar%sine_beta > isnight .and. input%parin > zero) then
       SDN(1) = 0
       ! Compute probability of penetration for direct and
       ! diffuse radiation for each layer in the canopy
       ! Level 1 is the soil surface and level ncl+1 is the
       ! top of the canopy. layer 1 is the layer above
       ! the soil and layer ncl is the top layer.
       ! diffuse NIR reflected by each layer
       ! compute the probability of diffuse radiation penetration through the
       ! hemisphere. this computation is not affected by penubra
       ! since we are dealing only with diffuse radiation from a sky
       ! sector.
       ! The probability of beam penetration is computed with a
       ! negative binomial distriubution for LAI.
       ! Radiation attenuation is a function of leaf and woody
       ! biomass
       ! Itegrated probability of diffuse sky radiation penetration
       ! for each layer
       !
       ! EXPDIF(JJ) is computed in PAR and can be used in NIR and IRFLUX
       reflectance_layer(1:ncl) = (one - solar%exxpdir(1:ncl)) * solar%nir_reflect
       ! DIFFUSE RADIATION TRANSMITTED THROUGH LAYER
       transmission_layer(1:ncl) = (one - solar%exxpdir(1:ncl)) * solar%nir_trans + solar%exxpdir(1:ncl)
       ! COMPUTE THE PROBABILITY OF beam PENETRATION
       ztmp = one / solar%sine_beta
       do j=ncl, 1, -1
          jp1 = j+1
          ! Probability of beam penetration.
          if (iswitch%wai_new == 0) then
          dff = prof%dLAIdz(j) ! + prof%dPAIdz(j)
            else
          dff = prof%dPAIdz(j) ! lai+wai Yuan 2018.03.01
          end if


          exp_direct = exp(-dff * markov*prof%Gfunc_solar(j) * ztmp)
          ! PEN1 = exp(-llai * prof%Gfunc_solar(JJ) / solar%sine_beta)
          ! exp_direct = exp(-DFF * prof%Gfunc_solar(JJ) / solar%sine_beta)
          ! Beam transmission
          beam(j)  = beam(jp1) * exp_direct
          TBEAM(j) = beam(j)
          SUP(jp1) = (TBEAM(jp1) - TBEAM(j)) * solar%nir_reflect
          SDN(j)   = (TBEAM(jp1) - TBEAM(j)) * solar%nir_trans
       end do
       ! initiate scattering using the technique of NORMAN (1979).
       ! scattering is computed using an iterative technique
       SUP(1)  = TBEAM(1) * solar%nir_soil_refl
       ADUM(1) = solar%nir_soil_refl
       do j=2, ncl+1
          jm1 = j-1
          TLAY2 = transmission_layer(jm1) * transmission_layer(jm1)
          ADUM(j) = ADUM(jm1) * TLAY2 / (one - ADUM(jm1) * reflectance_layer(jm1)) + reflectance_layer(jm1)
       end do
       do j=ncl, 1, -1
          jp1 = j+1
          solar%nir_dn(j)   = solar%nir_dn(jp1) * transmission_layer(j) / (one - ADUM(jp1) * reflectance_layer(j)) + SDN(j)
          solar%nir_up(jp1) = ADUM(jp1) * solar%nir_dn(jp1) + SUP(jp1)
       end do
       ! lower boundary: upward radiation from soil
       solar%nir_up(1) = solar%nir_soil_refl * solar%nir_dn(1) + SUP(1)
       ! Iterative calculation of upward diffuse and downward beam +
       ! diffuse NIR to compute scattering
       ITER  = 0
       IREP  = 1
       ITER  = ITER  + 1
       do while(IREP == 1)
          IREP=0
          do j=ncl, 1, -1
             jp1 = j+1
             DOWN = transmission_layer(j) * solar%nir_dn(jp1) + solar%nir_up(j) * reflectance_layer(j) + SDN(j)
             if (abs(DOWN-solar%nir_dn(j)) > e2) IREP = 1
             solar%nir_dn(j) = DOWN
          end do
          ! upward radiation at soil is reflected beam and downward diffuse
          solar%nir_up(1) = (solar%nir_dn(1) + TBEAM(1)) * solar%nir_soil_refl
          do j=2, ncl+1
             jm1 = j-1
             UP = reflectance_layer(jm1) * solar%nir_dn(j) + solar%nir_up(jm1) * transmission_layer(jm1) + SUP(j)
             if (abs(UP-solar%nir_up(j)) > e2) IREP = 1
             solar%nir_up(j) = UP
          end do
       end do
       ! Compute NIR flux densities
       ! upward diffuse NIR flux density, on the horizontal
       solar%nir_up(1:ncl+1)  = solar%nir_up(1:ncl+1)  * solar%nir_total
       where (solar%nir_up(1:ncl+1) < e1) solar%nir_up(1:ncl+1) = e1
       ! downward beam NIR flux density, incident on the horizontal
       solar%beam_flux_nir(1:ncl+1) = beam(1:ncl+1) * solar%nir_total
       where (solar%beam_flux_nir(1:ncl+1) < e1) solar%beam_flux_nir(1:ncl+1) = e1
       ! downward diffuse radiaiton flux density on the horizontal
       solar%nir_dn(1:ncl+1)  = solar%nir_dn(1:ncl+1)  * solar%nir_total
       where (solar%nir_dn(1:ncl+1) < e1) solar%nir_dn(1:ncl+1) = e1
       ! normal radiation on sunlit leaves
       if (solar%sine_beta > isnight .and. input%parin > zero) then
          nir_normal(1:ncl) = solar%nir_beam * prof%Gfunc_solar(1:ncl) * ztmp
       else
          nir_normal(1:ncl) = zero
       end if
       NSUNEN(1:ncl) = nir_normal(1:ncl) * solar%nir_absorbed
       ! Diffuse radiation received on top and bottom of leaves
       ! drive photosynthesis and energy exchanges
       solar%nir_shd(1:ncl) = (solar%nir_dn(1:ncl) + solar%nir_up(1:ncl))
       ! absorbed radiation, shaded
       solar%nir_shd(1:ncl) = solar%nir_shd(1:ncl)  * solar%nir_absorbed
       ! plus diffuse component
       solar%nir_sun(1:ncl) = NSUNEN(1:ncl) + solar%nir_shd(1:ncl)
    else ! solar%nir_total <= one .or. solar%sine_beta <= isnight .or. input%parin <= zero
       solar%nir_up(1:ncl)        = zero
       solar%nir_dn(1:ncl)        = zero
       solar%nir_shd(1:ncl)       = zero
       solar%nir_sun(1:ncl)       = zero
       solar%beam_flux_nir(1:ncl) = zero
       solar%nir_total            = zero
       solar%nir_up(ncl+1)        = zero
       solar%nir_dn(ncl+1)        = zero
       solar%beam_flux_nir(ncl+1) = zero
    end if

  END SUBROUTINE nir


  ! ------------------------------------------------------------------
  SUBROUTINE par()
    ! This subroutine computes the flux densities of direct and
    ! diffuse radiation using the measured leaf distrib.
    ! We apply the Markov model to account for clumping of leaves.

    ! The model is based on the scheme of Norman.

    ! Norman, J.M. 1979. Modeling the complete crop canopy.
    ! Modification of the Aerial Environment of Crops.
    ! B. Barfield and J. Gerber, Eds. American Society of Agricultural Engineers, 249-280.
    USE setup,      ONLY: ncl
    USE constants,  ONLY: zero, one, e2, e3, isnight
    USE types,      ONLY: solar, input, prof, debug, iswitch
    USE parameters, ONLY: markov

    IMPLICIT NONE

    REAL(wp), DIMENSION(ncl)   :: SDN, transmission_layer, reflectance_layer
    REAL(wp), DIMENSION(ncl)   :: par_normal_abs_energy, par_normal_abs_quanta, par_normal_quanta
    REAL(wp), DIMENSION(ncl+1) :: SUP, ADUM, TBEAM, beam
    REAL(wp) :: fraction_beam, sumlai
    REAL(wp) :: PEN2, dff
    REAL(wp) :: exp_direct, QU, TLAY2
    REAL(wp) :: DOWN, UP, ztmp

    INTEGER(i4) :: j, jp1, jm1
    INTEGER(i4) :: IREP, ITER

    if (solar%sine_beta > isnight .and. input%parin > zero) then
       fraction_beam = solar%par_beam / input%parin
       beam(ncl+1)   = fraction_beam
       TBEAM(ncl+1)  = fraction_beam
       SDN(1)        = zero
       ! Compute probability of penetration for direct and
       ! diffuse radiation for each layer in the canopy
       ! Level 1 is the soil surface and level ncl+1 is the
       ! top of the canopy. layer 1 is the layer above
       ! the soil and layer ncl is the top layer.
       ! Diffuse PAR reflected by each layer
       ! compute the probability of diffuse radiation penetration through the
       ! hemisphere. this computation is not affected by penubra
       ! since we are dealing only with diffuse radiation from a sky
       ! sector.
       ! The probability of beam penetration is computed with a
       ! negative binomial distriubution.
       reflectance_layer(1:ncl)  = (one - solar%exxpdir(1:ncl)) * solar%par_reflect

       ! DIFFUSE RADIATION TRANSMITTED THROUGH LAYER
       transmission_layer(1:ncl) = (one - solar%exxpdir(1:ncl)) * solar%par_trans + solar%exxpdir(1:ncl)

       ! COMPUTE THE PROBABILITY OF beam PENETRATION
       sumlai = zero
       do j=ncl, 1, -1
          jp1 = j+1
          ! Probability of beam penetration. This is that which
          ! is not umbral. Hence, this radiation
          ! is attenuated by the augmented leaf area: DF+PA.
          ! leaf and wood Yuan 2018.03.01
          if (iswitch%wai_new==0) then
           dff        = prof%dLAIdz(j) ! + prof%dPAIdz(JJ)
            else
           dff        = prof%dPAIdz(j)
          end if


!          print *, "sumlai"
!          print *, sumlai
          sumlai     = sumlai + dff ! sumlai includes wood and leaves
          exp_direct = exp(-dff*markov*prof%Gfunc_solar(j)/ solar%sine_beta)
          PEN2       = exp(-sumlai*markov*prof%Gfunc_solar(j)/ solar%sine_beta)
          ! lai Sunlit and shaded
          prof%sun_lai(j) = solar%sine_beta * (one-PEN2)/(markov*prof%Gfunc_solar(j))
          prof%shd_lai(j) = sumlai - prof%sun_lai(j)
!          print *, prof%sun_lai(j),prof%shd_lai(j),sumlai
          ! note that the integration of the source term time solar%prob_beam with respect to
          !   leaf area will yield the sunlit leaf area, and with respect to solar%prob_shd the
          !   shaded leaf area.
          !   In terms of evaluating fluxes for each layer
          !   Fcanopy = sum {fsun psun + fshade pshade} (see Leuning et al. Spitters et al.)
          !   psun is equal to exp(-lai G markov/sinbet)
          !   pshade = 1 - psun
          solar%prob_beam(j) = markov*PEN2
 !      if (ISNAN(solar%prob_beam(j))) then
 !       print *, solar%prob_beam(j), prof%Gfunc_solar(j), solar%sine_beta, markov, PEN2, exp_direct, sumlai , dff

 !       print *, prof%dWAIdz(j)

!debug%D1 = wai_z
!debug%D2 = beta_fnc_w
!debug%D3 = wai_midpt
!debug%D4 = wai_freq
!debug%R1 = time%wai
!debug%R2 = integr_beta_w
!debug%in1 = nbeta_w
!    print *, "debug:"
!    print *, "wai_z="
!        print *, debug%D1
!print *, "beta_fnc_w_z="
!        print *, debug%D2
!print *, "wai_midpt="
!        print *, debug%D3
!print *, "wai_freq="
!        print *, debug%D4
!print *, "time%wai="
!        print *, debug%R1
!print *, "integr_beta_w="
!        print *, debug%R2
!print *, "sum_beta="
!print *, debug%R3

!print *, "nbeta_w="
!        print *, debug%in1







!       end if
          ! probability of beam
          beam(j) = beam(jp1) * exp_direct
          QU = one - solar%prob_beam(j)
          if (QU > one)  QU = one
          if (QU < zero) QU = zero
          ! probability of umbra
          solar%prob_shd(j) = QU
!       if (ISNAN(solar%prob_shd(j)) ) then
!        print *, debug%D1
!        print *, debug%D2
!        print *, debug%R1
!        print *, debug%R2
!       end if
          TBEAM(j) = beam(j)
          ! beam PAR that is reflected upward by a layer
          SUP(jp1) = (TBEAM(jp1) - TBEAM(j)) * solar%par_reflect
          ! beam PAR that is transmitted downward
          SDN(j) = (TBEAM(jp1) - TBEAM(j)) * solar%par_trans
       end do

       ! initiate scattering using the technique of NORMAN (1979).
       ! scattering is computed using an iterative technique.
       ! Here Adum is the ratio up/down diffuse radiation.
       SUP(1)                = TBEAM(1) * solar%par_soil_refl
       solar%par_down(ncl+1) = one - fraction_beam
       ADUM(1)               = solar%par_soil_refl
       do j=1, ncl
          jp1 = j+1
          TLAY2 = transmission_layer(j) * transmission_layer(j)
          ADUM(jp1) = ADUM(j) * TLAY2 / (one - ADUM(j) * reflectance_layer(j)) + reflectance_layer(j)
       end do
       do j=ncl, 1, -1
          jp1 = j+1
          solar%par_down(j) = solar%par_down(jp1) * transmission_layer(j) / (one - ADUM(jp1) * reflectance_layer(j)) + SDN(j)
          solar%par_up(jp1) = ADUM(jp1) * solar%par_down(jp1) + SUP(jp1)
       end do
       ! lower boundary: upward radiation from soil
       solar%par_up(1) = solar%par_soil_refl * solar%par_down(1) + SUP(1)
       ! Iterative calculation of upward diffuse and downward beam +
       ! diffuse PAR.
       ! This section has been commented out for the negative binomial
       ! model. It seems not to apply and incorrectly calculates
       ! scattering. When I ignore this section, I get perfect
       ! agreement between measured and calculated Rn.
       ! Scattering
       ITER = 0
       IREP = 1
       do while(IREP == 1)
          IREP = 0
          ITER = ITER + 1
          do j=ncl, 1, -1
             jp1 = j+1
             DOWN = transmission_layer(j) * solar%par_down(jp1) + solar%par_up(j) * reflectance_layer(j) + SDN(j)
             if (abs(DOWN-solar%par_down(j)) > e2) IREP = 1
             solar%par_down(j) = DOWN
          end do
          ! upward radiation at soil is reflected beam and downward diffuse
          solar%par_up(1) = (solar%par_down(1) + TBEAM(1)) * solar%par_soil_refl
          do j=2, ncl+1
             jm1 = j-1
             UP = reflectance_layer(jm1) * solar%par_down(j) + solar%par_up(jm1) * transmission_layer(jm1) + SUP(j)
             if (abs(UP-solar%par_up(j)) > e2) IREP = 1
             solar%par_up(j) = UP
          end do
       end do
       ! Compute flux density of PAR
       solar%par_up(1:ncl+1) = solar%par_up(1:ncl+1) * input%parin
       where (solar%par_up(1:ncl+1) < e3) solar%par_up(1:ncl+1) = e3
       ! downward beam PAR flux density, incident on the horizontal
       solar%beam_flux_par(1:ncl+1) = beam(1:ncl+1) * input%parin
       !print *, solar%beam_flux_par
      ! print *, input%parin
      ! print *, beam(1:ncl+1)
      ! print *, sum(beam(1:ncl))
       where (solar%beam_flux_par(1:ncl+1) < e3) solar%beam_flux_par(1:ncl+1) = e3
       ! Downward diffuse radiatIon flux density on the horizontal
       solar%par_down(1:ncl+1) = solar%par_down(1:ncl+1) * input%parin
       where (solar%par_down(1:ncl+1) < e3) solar%par_down(1:ncl+1) = e3

       if (solar%par_beam < e3) solar%par_beam = e3
       ! PSUN is the radiation incident on the mean leaf normal

       if (solar%sine_beta > isnight .and. input%parin > zero) then
          ztmp = one / solar%sine_beta
          ! PAR received normal to a leaf on a sunlit spot
          par_normal_quanta(1:ncl) = solar%par_beam * prof%Gfunc_solar(1:ncl) * ztmp
       else
          par_normal_quanta(1:ncl) = zero
       end if
       ! amount of energy absorbed by sunlit leaf
       par_normal_abs_energy(1:ncl) = par_normal_quanta(1:ncl) * solar%par_absorbed / 4.6_wp ! W m-2
       par_normal_abs_quanta(1:ncl) = par_normal_quanta(1:ncl) * solar%par_absorbed ! umol m-2 s-1
       ! Convert PAR to W m-2 for energy balance computations
       ! but remember umol m-2 s-1 is needed for photosynthesis.
       ! Average fluxes for the values on the top and bottom of
       ! each layer that drives energy fluxes.

       ! Energy balance computations are on the basis of
       ! absorbed energy. Harley''s photosynthesis
       ! parameterizations are on the basis of incident
       ! PAR
       solar%quantum_shd(1:ncl) = (solar%par_down(1:ncl) + solar%par_up(1:ncl)) * solar%par_absorbed ! umol m-2 s-1
       solar%quantum_sun(1:ncl) = solar%quantum_shd(1:ncl) + par_normal_abs_quanta(1:ncl)
       ! calculate absorbed par
       solar%par_shd(1:ncl)     = solar%quantum_shd(1:ncl) / 4.6_wp ! W m-2
       ! solar%par_sun is the total absorbed radiation on a sunlit leaf,
       ! which consists of direct and diffuse radiation
       solar%par_sun(1:ncl) = par_normal_abs_energy(1:ncl) + solar%par_shd(1:ncl)
    else ! solar%sine_beta <= isnight .or. input%parin <= zero
       solar%prob_shd(1:ncl)      = one
       solar%prob_beam(1:ncl)     = zero
       solar%par_up(1:ncl)        = zero
       solar%par_down(1:ncl)      = zero
       solar%par_sun(1:ncl)       = zero
       solar%par_shd(1:ncl)       = zero
       solar%beam_flux_par(1:ncl) = zero
       solar%quantum_shd(1:ncl)   = zero
       solar%quantum_sun(1:ncl)   = zero
       prof%sun_lai(1:ncl)        = zero
       prof%shd_lai(1:ncl)        = prof%dLAIdz(1:ncl)
       solar%par_beam             = zero
       solar%par_up(ncl+1)        = zero
       solar%par_down(ncl+1)      = zero
       solar%beam_flux_par(ncl+1) = zero
    end if

  END SUBROUTINE par


  ! ------------------------------------------------------------------
  SUBROUTINE rnet()
    ! Radiation layers go from jtot+1 to 1. jtot+1 is the top of the
    ! canopy and level 1 is at the soil surface.
    !
    ! Energy balance and photosynthesis are performed for vegetation
    ! between levels and based on the energy incident to that level
    USE setup,      ONLY: ncl
    USE types,      ONLY: solar
    USE parameters, ONLY: ep

    IMPLICIT NONE

    REAL(wp), DIMENSION(1:ncl) :: ir_shade

    ! Infrared radiation on leaves
    ir_shade(1:ncl) = (solar%ir_dn(1:ncl) + solar%ir_up(1:ncl)) * ep
    ! Available energy on leaves for evaporation.
    ! Values are average of top and bottom levels of a layer.
    ! The index refers to the layer. So layer 3 is the average
    ! of fluxes at level 3 and 4. Level 1 is soil and level
    ! j+1 is the top of the canopy. layer ncl is the top layer

    ! Sunlit, shaded values
    solar%rnet_sun(1:ncl) = solar%par_sun(1:ncl) + solar%nir_sun(1:ncl) + ir_shade(1:ncl)
 !print *, solar%rnet_sun(ncl)
 !print *, solar%par_sun(ncl)
 !print *, solar%nir_sun(ncl)
 !print *, ir_shade(ncl)
 !print *, solar%ratrad
 !print *, solar%ir_up(ncl+1)
 !print *, solar%ir_dn(ncl+1)
    solar%rnet_shd(1:ncl) = solar%par_shd(1:ncl) + solar%nir_shd(1:ncl) + ir_shade(1:ncl)

  END SUBROUTINE rnet


  ! ------------------------------------------------------------------
  SUBROUTINE set_leaf_phenology()

    USE types,      ONLY: time
    USE parameters, ONLY: leaf_out, leaf_full, leaf_fall, leaf_fall_complete

    IMPLICIT NONE

    !changed: now fixed to input file, old: assume leaf out starts when soil temperature reaches 13 C at 32 cm
    time%leafout          = leaf_out !changed: (int)(365.*(1.5707+asin(14.-soil.Temp_ref)/soil.amplitude)/6.283)
    ! assume 'leaf_full' days till full leaf
    time%leaffull         = leaf_full ! changed: old=+30
    ! leaf shedding
    time%leaffall         = leaf_fall
    time%leaffallcomplete = leaf_fall_complete

  END SUBROUTINE set_leaf_phenology


  ! ------------------------------------------------------------------
  FUNCTION sky_ir(T)
    ! Infrared radiation from sky, W m-2, using algorithm from Norman
    USE types,     ONLY: solar
    USE constants, ONLY: one, sigma, TN0

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: T ! Temp [K]
    REAL(wp) :: sky_ir

    ! sigma * T**4 * ((1.-0.261*exp(-.000777 * (TN0-T)**))*solar.ratrad + 1.-solar.ratrad)
    sky_ir = sigma * T**4 * ((one -0.261_wp*exp(-0.000777_wp * (T-TN0)**2))*solar%ratrad + one-solar%ratrad)

  END FUNCTION sky_ir


END MODULE radiation
