MODULE nitrogen_assimilation


  USE kinds, ONLY: wp, i4
  USE constants,  ONLY: one, zero

 IMPLICIT NONE

  PRIVATE

  PUBLIC :: N_fraction, NC_canopy, N_top


  CONTAINS
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
END SUBROUTINE random_stduniform

SUBROUTINE random_stdnormal(x)
    use constants, ONLY: pi
   implicit none
   real(wp),intent(out) :: x
   real(wp) :: u1,u2
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)
END SUBROUTINE random_stdnormal

SUBROUTINE N_fraction (x1,x2,x3)

    ! calculate fractions of N sources
    USE parameters, ONLY: nitrate, nitrite, ammonia
    USE types,      ONLY: iswitch
    implicit none

    real(wp),intent(out) :: x1,x2,x3
    real(wp) :: n1, n2, n3, tmp

    SELECT CASE (iswitch%n_random)
   CASE (0)
      x1 = nitrate
      x2 = nitrite
      x3 = ammonia
   CASE (1)
      call random_stdnormal (n1)
      call random_stdnormal (n2)
      call random_stdnormal (n3)
      n1 = abs(n1)
      n2 = abs(n2)
      n3 = abs(n3)
      x1 = n1/(n1+n2+n3)
      x2 = n2/(n1+n2+n3)
      x3 = n3/(n1+n2+n3)
    CASE (2)
      call random_uniform(0.8_wp,0.9_wp,n1)! determine the fraction of nitrate
      tmp = one-n1
      call random_uniform(zero,tmp,n2)
      n3 = tmp-n2
      x1 = n1
      x2 = n2
      x3 = n3
END SELECT


END SUBROUTINE N_fraction

FUNCTION NC_canopy (JJ)

    USE types,      ONLY: iswitch
    USE parameters, ONLY: cn_bulk

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: JJ
    REAL(wp)             :: tmp, NC_canopy

    if (iswitch%cn==0) then
        tmp = cn_bulk
    else
        tmp = 14.55608 + 11.53367* (JJ/40)
    end if
    NC_canopy = one/tmp

END FUNCTION NC_canopy

FUNCTION N_top (J_ref,N_ref)

    ! derive N supply of top canopy layer from actual N assimilation of our measurement layer (J_ref),&
    ! which is deduced from O2 flux in subroutine "O_to_N"
    USE types,      ONLY: iswitch, time, prof
    USE parameters, ONLY: cn_bulk, htFrac, zh65, lai

    IMPLICIT NONE

    INTEGER(i4), INTENT(IN) :: J_ref
    REAL(wp), INTENT(IN)    :: N_ref
    REAL(wp)                :: N_top
    ! N_ref = N_top * (zh65 * zzz + (1-htFrac)) * time%lai/lai
    N_top = N_ref/(zh65 * prof%ht(J_ref) + (1-htFrac))
    !print *, zh65, prof%ht(J_ref)

END FUNCTION N_top

END MODULE nitrogen_assimilation
