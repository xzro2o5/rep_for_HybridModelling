MODULE setup

  ! This module contains the model''s dimensions such as number of layers, etc.

  ! Written Jan 2011, Matthias Cuntz - Ported C-Code

  USE kinds, ONLY: i4

  IMPLICIT NONE

  INTEGER(i4) :: ncl         ! # of canopy layers
  INTEGER(i4) :: ntl         ! total # of layers = 3*ncl
  INTEGER(i4) :: nsky        ! # of sky angle classes
  INTEGER(i4) :: nl          ! # of leaf angle classes
  INTEGER(i4) :: nsoil       ! # of soil layers
  INTEGER(i4) :: nbeta       ! # of levels for beta distribution of e.g. lai
  INTEGER(i4) :: ndaysc13    ! # of days to remember for mean 13C discrimination
  INTEGER(i4) :: nwiso       ! # of maximum number of isotopic waters
  INTEGER(i4) :: nlop        ! # of leaf optical properties
  
END MODULE setup
