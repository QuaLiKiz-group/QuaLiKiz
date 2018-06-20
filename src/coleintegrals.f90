MODULE coleintegrals
  USE kind

  IMPLICIT NONE

CONTAINS
  SUBROUTINE coleint(ndim,aa,bb,minpts,maxpts,integrand,relacc2,acc,lenwrk,wrkstr,COLEresult,ifailloc)


    INTEGER, INTENT(IN) :: ndim,minpts,maxpts,ifailloc,lenwrk
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: relacc2,acc
    REAL(KIND=DBL), INTENT(OUT) :: COLEresult
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: wrkstr
    REAL(KIND=DBL) :: integrand

    REAL(KIND=DBL), DIMENSION(ndim) :: xy !local position
    REAL(KIND=DBL) :: count

    !just an example
    xy(1) = aa(1)+0.1
    xy(2) = bb(1)+0.1
    count = integrand(ndim,xy)

    COLEresult = 42./4.
    
  END SUBROUTINE coleint

END MODULE coleintegrals
