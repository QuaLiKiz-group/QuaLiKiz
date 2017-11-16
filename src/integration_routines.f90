!
MODULE integration_routines
  USE kind
  USE datcal
  USE datmat
  USE calltrapints
  USE callpassints

  IMPLICIT NONE
INTERFACE integrate_1d
  MODULE PROCEDURE integrate_1d_nag
END INTERFACE

CONTAINS
  REAL(KIND=DBL) FUNCTION integrate_1d_nag(cc, dd, relacc1, npts, relerr, lw, integrand, verbose, p, nu, omega, intname)
!     INPUT ARGUMENTS
!     ----- ----------
!     (A,B) cc,dd     -  LOWER AND UPPER INTEGRATION LIMITS.
!     (EPR) relacc1     -  REQUIRED RELATIVE ACCURACY.
!     (NL) lw      -  APPROXIMATE LIMIT ON NUMBER OF INTEGRAND
!                 EVALUATIONS. IF SET NEGATIVE OR ZERO THE
!                 DEFAULT IS 10000.
!     (F) integrand       -  THE USER NAMED AND PREPARED FUNCTION  F(X)
!                 GIVES THE VALUE OF THE INTEGRAND AT X.
!     (IFAIL) ifailloc      INTEGER VARIABLE
!             - 0  FOR HARD FAIL REPORT
!             - 1  FOR SOFT FAIL REPORT
!
!     OUTPUT ARGUMENTS
!     ------ ----------
!     (NPTS) npts    -  NUMBER OF INTEGRAND EVALUATIONS USED IN OBTAINING
!                THE RESULT.
!     (RELERR) relerr  -  ROUGH ESTIMATE OF RELATIVE ACCURACY ACHIEVED.
!     IFAIL   -  VALUE INDICATES THE OUTCOME OF THE INTEGRATION -
    REAL(KIND=DBL), INTENT(IN)     :: cc, dd, relacc1
    INTEGER, INTENT(IN)            :: npts, lw  !output number of integrand evaluations
    EXTERNAL :: integrand

    LOGICAL, INTENT(IN) :: verbose
    INTEGER, INTENT(IN)  :: p, nu
    CHARACTER(*), INTENT(IN) :: intname
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega

    REAL(KIND=DBL) :: relerr
    INTEGER :: ifailloc

    REAL(KIND=DBL), DIMENSION(8) :: K

    !WRITE(stderr, "(A)") intname
    ifailloc = 1
    !integrate_1d_nag = d01ahf(cc, dd, relacc1, npts, relerr, integrand, lw, ifailloc)
    CALL QUAD(cc, dd, K, integrate_1d_nag, relacc1, npts, ifailloc, integrand)
    !WRITE(stderr, "(G10.3)") K
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,A,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of ', intname, ' integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF
  END FUNCTION

END MODULE integration_routines
