!
MODULE integration_routines
  USE kind
  USE datcal

  IMPLICIT NONE

CONTAINS
  REAL(KIND=DBL) FUNCTION integrate_1d(cc, dd, relacc, npts, relerr, integrand, lw, verbose, p, nu, omega, intname)

    REAL(KIND=DBL), EXTERNAL        :: integrand
    REAL(KIND=DBL), EXTERNAL        :: d01ahf

    REAL(KIND=DBL), INTENT(IN)      :: cc, dd, relacc
    INTEGER, INTENT(IN)             :: lw  !output number of integrand evaluations

    REAL(KIND=DBL), INTENT(OUT)     :: relerr

    ! For debugging printout
    LOGICAL, INTENT(IN)             :: verbose
    INTEGER, INTENT(IN)             :: p, nu
    CHARACTER(*), INTENT(IN)        :: intname
    COMPLEX(KIND=DBL), INTENT(IN)   :: omega

    INTEGER                         :: ifailloc, npts

    ifailloc = 1
    CALL integrate_1d_nag(cc, dd, relacc, npts, relerr, integrand, lw, integrate_1d, ifailloc)

    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,A,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of ', intname, ' integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

  END FUNCTION

  SUBROUTINE integrate_1d_nag(cc, dd, relacc, npts, relerr, integrand, lw, res, ifailloc)
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

    REAL(KIND=DBL), INTENT(IN)      :: cc, dd, relacc
    INTEGER, INTENT(IN)             :: lw  !output number of integrand evaluations

    REAL(KIND=DBL), INTENT(OUT)     :: relerr, res
    INTEGER, INTENT(OUT)            :: npts, ifailloc

    ifailloc = 1
    res = 0
    !res = d01ahf(cc, dd, relacc, npts, relerr, integrand, lw, ifailloc)

    !Netlib trial
    !CALL QUAD(cc, dd, K, integrate_1d_nag, relacc1, npts, ifailloc, integrand)
  END SUBROUTINE

  REAL(KIND=DBL) FUNCTION integrate_2d(ndim, a, b, minpts, maxpts, integrand, relacc, acc, lenwrk, wrkstr, verbose, p, nu, omega, intname)
    REAL(KIND=DBL), EXTERNAL        :: integrand

    REAL(KIND=DBL), DIMENSION(:), INTENT(IN)      :: a, b, wrkstr
    REAL(KIND=DBL), INTENT(IN)      :: relacc
    INTEGER, INTENT(IN)             :: ndim, minpts, maxpts, lenwrk

    REAL(KIND=DBL), INTENT(OUT)     :: acc
    INTEGER                         :: ifailloc

    ! For debugging printout
    LOGICAL, INTENT(IN)             :: verbose
    INTEGER, INTENT(IN)             :: p, nu
    CHARACTER(*), INTENT(IN)        :: intname
    COMPLEX(KIND=DBL), INTENT(IN)   :: omega

    ifailloc = 1
    CALL integrate_2d_nag(ndim, a, b, minpts, maxpts, integrand, relacc, acc, lenwrk, wrkstr, integrate_2d, ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 2DNAG rFFke integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF
  END FUNCTION

  SUBROUTINE integrate_2d_nag(ndim, a, b, minpts, maxpts, integrand, relacc, acc, lenwrk, wrkstr, res, ifailloc)
    !      INPUT PARAMETERS
    !
    !     NDIM    INTEGER NUMBER OF VARIABLES, MUST EXCEED 1 BUT
    !             NOT EXCEED 15.
    !
    !     A       REAL ARRAY OF LOWER LIMITS, WITH DIMENSION NDIM
    !
    !     B       REAL ARRAY OF UPPER LIMITS, WITH DIMENSION NDIM
    !
    !     MINPTS  INTEGER MINIMUM NUMBER OF INTEGRAND VALUES TO BE
    !             ALLOWED, WHICH MUST NOT EXCEED MAXPTS.
    !
    !     MAXPTS  INTEGER MAXIMUM NUMBER OF INTEGRAND VALUES TO BE
    !             ALLOWED, WHICH MUST BE AT LEAST
    !             2**NDIM+2*NDIM**2+2*NDIM+1.
    !
    !     FUNCTN  EXTERNALLY DECLARED USER DEFINED REAL FUNCTION
    !             INTEGRAND. IT MUST HAVE PARAMETERS (NDIM,Z),
    !             WHERE Z IS A REAL ARRAY OF DIMENSION NDIM.
    !
    !     EPS     REAL REQUIRED RELATIVE ACCURACY, MUST BE GREATER
    !             THAN ZERO
    !
    !     LENWRK  INTEGER LENGTH OF ARRAY WRKSTR, MUST BE AT LEAST
    !             2*NDIM+4.
    !
    !     IFAIL   INTEGER NAG FAILURE PARAMETER
    !             IFAIL=0 FOR HARD FAIL
    !             IFAIL=1 FOR SOFT FAIL
    !
    !      OUTPUT PARAMETERS
    !
    !     MINPTS  INTEGER NUMBER OF INTEGRAND VALUES USED BY THE
    !             ROUTINE
    !
    !     WRKSTR  REAL ARRAY OF WORKING STORAGE OF DIMENSION (LENWRK).
    !
    !     ACC     REAL ESTIMATED RELATIVE ACCURACY OF FINVAL
    !
    !     FINVAL  REAL ESTIMATED VALUE OF INTEGRAL
    !
    !     IFAIL   IFAIL=0 FOR NORMAL EXIT, WHEN ESTIMATED RELATIVE
    !                  LESS INTEGACCURACY RAND VALUES USED.
    REAL(KIND=DBL), EXTERNAL        :: integrand
    REAL(KIND=DBL), EXTERNAL        :: d01fcf

    REAL(KIND=DBL), DIMENSION(:), INTENT(IN)      :: a, b, wrkstr
    REAL(KIND=DBL), INTENT(IN)      :: relacc
    INTEGER, INTENT(IN)             :: ndim, minpts, maxpts, lenwrk

    REAL(KIND=DBL), INTENT(OUT)     :: acc, res
    INTEGER, INTENT(OUT)            :: ifailloc

    ifailloc = 1
    res = 0 
    !res = d01fcf(ndim, a, b, minpts, maxpts, integrand, relacc, acc, lenwrk, wrkstr, res, ifailloc)

  END SUBROUTINE

END MODULE integration_routines
