MODULE nitrogen_assimilation


  USE kinds, ONLY: wp, i4
  USE constants,  ONLY: one, zero

 IMPLICIT NONE

  PRIVATE

  PUBLIC :: N_assimilation


  CONTAINS

  SUBROUTINE N_assimilation(gpp,J_co2,gly,serine,source_NO3,source_NO2,source_NH4,layer)
    ! to calculate N assimilation accompany with carbon assimilation
    ! electron requirements
    ! assuming N reducted to glutamate, CO2 reducted to carbohydrate,
    ! co2~4e, nitrate~10e, nitrite~8e and ammonia~2e
    ! J_extra/J_a =(e_source)*ncbulk/e_co2 without N limit
    ! J_extra/J_a =(e_source)*N_supply/(C_ass*e_co2) under N limit

    USE constants,  ONLY: one
    USE types,      ONLY: iswitch, nitrogen, time
    USE parameters, ONLY: nc_bulk, n_supply, n_mult, nitrate, nitrite, ammonia, n_max
    USE setup,      ONLY: ncl

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: gpp
    REAL(wp), INTENT(IN) :: J_co2
    REAL(wp), INTENT(IN) :: gly
    REAL(wp), INTENT(IN) :: serine
    REAL(wp), INTENT(OUT) :: source_NO3,source_NO2,source_NH4
    INTEGER(i4), INTENT(IN) :: layer

    !REAL(wp) :: J_co2
    REAL(wp) :: J_extra ! total electrons for N reduction
    REAL(wp) :: n_ass ! total assimilated N in umol m-2 s-1
    REAL(wp) :: source_nitrate ! supplied nitrate as N source, umol m-2 s-1
    REAL(wp) :: source_nitrite ! supplied nitrite as N source, umol m-2 s-1
    REAL(wp) :: source_ammonia ! supplied ammonia as N source, umol m-2 s-1
    REAL(wp) :: J_nitrate ! electrons for nitrate reduction, umol m-2 s-1
    REAL(wp) :: J_nitrite ! electrons for nitrite reduction, umol m-2 s-1
    REAL(wp) :: J_ammonia ! electrons for ammonia reduction, umol m-2 s-1
    REAL(wp) :: N_C ! leaf N:C ratio
    REAL(wp) :: test

    if (iswitch%n_limit==0) then  ! 0 without N limit; 1 with N limit; 2 multiply of gly + serine

        !
        N_C = nc_bulk! default N:C ratio
        !n_ass = nc_bulk*gpp
        ! use the vertical N:C profile (Bachofen et al. 2020)

        ! spring
!        if (time%days >= time%leafout .and. time%days < time%leaffull) then
!            if (layer<=10) then
!                N_C = 0.065
!            else if (layer >10 .and. layer <= 20) then
!                N_C = 0.066
!            else if (layer >20 .and. layer <= 30) then
!                N_C = 0.071
!            else if (layer >30 .and. layer <= 40) then
!                N_C = 0.07
!            end if
!
!        end if
!
!!print *, time%days
!        ! summer
!        if (time%days >= time%leaffull .and. time%days <= time%leaffall) then
!            if (layer<=10) then
!                N_C = 0.0547
!            else if (layer >10 .and. layer <= 20) then
!                N_C = 0.0531
!            else if (layer >20 .and. layer <= 30) then
!                N_C = 0.0502
!            else if (layer >30 .and. layer <= 40) then
!                N_C = 0.0498
!            end if
!
!        end if
        n_ass = N_C*gpp


    else if (iswitch%n_limit==1) then

        n_ass = min(n_supply,nc_bulk*gpp) ! or n_supply/ncl, per layer

    else if (iswitch%n_limit==2) then

        n_ass = n_mult*(gly+serine)

    end if

    if (iswitch%n_random==1) then

       ! create random N source fraction:
       !print *, nitrate
       !call random_number(nitrate)
       call random_uniform(0.2_wp,0.3_wp,nitrate)
       call random_uniform(zero,(1-nitrate),nitrite)
       ammonia = 1-nitrate-nitrite
       !print *, nitrate, nitrite, ammonia
    end if


    source_NO3 = n_ass*nitrate
    source_NO2 = n_ass*nitrite
    source_NH4 = n_ass*ammonia

!    J_nitrate = J_co2*(10*source_nitrate/gpp)/4
!    J_nitrite = J_co2*(8*source_nitrite/gpp)/4
!    J_ammonia = J_co2*(2*source_ammonia/gpp)/4

! from nitrogen to electron
    J_nitrate = 10*source_NO3
    J_nitrite = 8*source_NO2
    J_ammonia = 2*source_NH4
    J_extra = J_nitrate + J_nitrite + J_ammonia

    ! save N source in structures:
    nitrogen%J_extra = J_extra
    nitrogen%nitrate_mol = source_NO3
    nitrogen%nitrite_mol = source_NO2
    nitrogen%ammonia_mol = source_NH4

!call random_uniform(0.4,0.8,test)
!print *, test

  END SUBROUTINE N_assimilation

! assuming a<b
! generate random nitrate, nitrite and ammonia fraction
SUBROUTINE random_uniform(a,b,x)
   implicit none
   REAL(wp),intent(in) :: a,b
   !real,intent(out) :: x
   REAL(wp),intent(out) :: x
   REAL(wp) :: u
   call random_stduniform(u)
   x = (b-a)*u + a
END SUBROUTINE random_uniform

SUBROUTINE random_stduniform(u)
   implicit none
   real(wp),intent(out) :: u
   real(wp) :: r
   call random_number(r)
   u = 1 - r
end SUBROUTINE random_stduniform


END MODULE nitrogen_assimilation
