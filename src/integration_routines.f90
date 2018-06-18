!
MODULE integration_routines
  USE kind
  USE datcal

  IMPLICIT NONE

CONTAINS
  REAL(KIND=DBL) FUNCTION integrate_1d(cc, dd, relacc1, npts, relerr, lw, integrand, verbose, p, nu, omega, intname)

    REAL(KIND=DBL), EXTERNAL        :: integrand
    REAL(KIND=DBL), EXTERNAL        :: d01ahf

    REAL(KIND=DBL), INTENT(IN)      :: cc, dd, relacc1
    INTEGER, INTENT(IN)             :: lw  !output number of integrand evaluations

    REAL(KIND=DBL), INTENT(OUT)     :: relerr


    LOGICAL, INTENT(IN)             :: verbose
    INTEGER, INTENT(IN)             :: p, nu
    CHARACTER(*), INTENT(IN)        :: intname
    COMPLEX(KIND=DBL), INTENT(IN)   :: omega

    INTEGER                         :: ifailloc, npts

    ifailloc = 1
    CALL integrate_1d_nag(cc, dd, relacc1, npts, relerr, integrand, lw, integrate_1d, ifailloc)

    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,A,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of ', intname, ' integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

  END FUNCTION

  SUBROUTINE integrate_1d_nag(cc, dd, relacc1, npts, relerr, integrand, lw, res, ifailloc)
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
    REAL(KIND=DBL), EXTERNAL        :: integrand
    REAL(KIND=DBL), EXTERNAL        :: d01ahf

    REAL(KIND=DBL), INTENT(IN)      :: cc, dd, relacc1
    INTEGER, INTENT(IN)             :: lw  !output number of integrand evaluations

    REAL(KIND=DBL), INTENT(OUT)     :: relerr, res
    INTEGER, INTENT(OUT)            :: npts, ifailloc

    ifailloc = 1
    res = d01ahf(cc, dd, relacc1, npts, relerr, integrand, lw, ifailloc)

    !Netlib trial
    !CALL QUAD(cc, dd, K, integrate_1d_nag, relacc1, npts, ifailloc, integrand)
  END SUBROUTINE

END MODULE integration_routines
