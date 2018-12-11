MODULE datcal
  USE kind
  IMPLICIT NONE

  !Numerical and physical constants
  REAL(KIND=DBL), PARAMETER :: pi   = 3.14159265358979
  REAL(KIND=DBL), PARAMETER :: twopi   = 6.28318530718
  REAL(KIND=DBL), PARAMETER :: sqrtpi = 1.77245385090552
  REAL(KIND=DBL), PARAMETER :: qe = 1.602176565d-19 !SI electron charge
  REAL(KIND=DBL), PARAMETER :: me = 9.10938291d-31 !SI electron mass
  REAL(KIND=DBL), PARAMETER :: mp = 1.672621777d-27 !SI proton mass
  REAL(KIND=DBL), PARAMETER :: Ze = -1._DBL  !electron charge in qe units
  COMPLEX(KIND=DBL),PARAMETER :: ci   = (0.,1.) !i

  !Code definitions
  REAL(KIND=DBL), PARAMETER :: eps=2.**(-52.) !small number epsilon
  REAL(KIND=DBL), PARAMETER :: epsD=1.d-14 !double precision epsilon
  REAL(KIND=DBL), PARAMETER :: epsS=1.d-7 !single precision epsilon

  !FLR parameters
  REAL(KIND=DBL), PARAMETER :: epsFLR = 2.0d-3 !Relative accuracy for FLR integrations
  INTEGER, PARAMETER :: nlimit=10000    !limit of integrand evaluations in integral functions
  REAL(KIND=DBL), PARAMETER :: normkrfac = 20 !Scale factor to set kr integration boundaries in FLR. Arbitrary but can optimize integration

  !Input and output variables
  INTEGER, PARAMETER :: n_order_max = 3
  INTEGER, PARAMETER :: ktype = 1 ! BINARY FILES
  LOGICAL, PARAMETER :: print_binary = .FALSE. !toggle printing of bin BUFF files
  LOGICAL, PARAMETER :: print_ascii = .TRUE. !toggle printing of output in ASCII files

  !Sets if we calculate the trapped and circulated particles
  LOGICAL, PARAMETER :: calctrap = .TRUE.
  LOGICAL, PARAMETER :: calccirc = .TRUE.

  !Sets if we calculate only positive or only negative frequency modes
  LOGICAL, PARAMETER :: onlyion = .FALSE.
  LOGICAL, PARAMETER :: onlyelec = .FALSE.
  REAL(KIND=DBL), PARAMETER :: impfac=0.05

  !Parameter for setting min k of ETG regime
  REAL(KIND=DBL), PARAMETER :: ETGk = 2._DBL

  !Contour parameters
  INTEGER, PARAMETER :: MM=17 ! Initial number of segments the contour is split into. Captures rectangle corners
  INTEGER, PARAMETER :: maxM=MM*(2**4) !Highest allowed  number of segments in the contour

  REAL(KIND=DBL), PARAMETER :: maxangle = pi/3. !Maximum angle above when neighbouring fonct points cannot go above 
  REAL(KIND=DBL), PARAMETER :: mincont = 0.05 !Minimum imaginary value in contour
  REAL(KIND=DBL), PARAMETER :: Elli = 0.1 !Ellipicity of contour parameter
  REAL(KIND=DBL), PARAMETER :: soldel=0.05 !Parameter for removing duplicate solutions
  REAL(KIND=DBL), PARAMETER :: nearlysol=10. !Determines if a solution is close enough for the Newton step
  REAL(KIND=DBL), PARAMETER :: overlapfac = 0.5 !percentage of contour overlap
  REAL(KIND=DBL), PARAMETER :: squirclecoef = 0.975 !default is 0.975
  REAL(KIND=DBL), PARAMETER :: centeroverlapfac = 0.1 !percentage of center contour overlap
  REAL(KIND=DBL), PARAMETER :: centerwidth = 2. !for narrow center contour

  !Integration parameters
  INTEGER, PARAMETER :: ndim=2 !set dimension in multidimensional ints
  INTEGER, PARAMETER :: nf=2 !sets number of components in multidimensional ints

  !Integration limits
  REAL(KIND=DBL), PARAMETER :: vuplim = 5. ! upper limit for v integration trapped particles
  REAL(KIND=DBL), PARAMETER :: rkuplim = 5.d0 * 1. ! limit for k and r passing integrations

!  REAL(KIND=DBL), PARAMETER :: barelyavoid = 0.01d0 ! Cutoff value below 1 to remove barely trapped singularity

  REAL(KIND=DBL), PARAMETER :: barelyavoid = epsD ! Cutoff value below 1 to remove barely trapped singularity
!  REAL(KIND=DBL), PARAMETER :: barelyavoid = 0.04 ! Cutoff value below 1 to remove barely trapped singularity
  REAL(KIND=DBL), PARAMETER :: minfki = 1d-3 ! Minimum allowed absolute vertical drift frequency
  REAL(KIND=DBL), PARAMETER :: min_ninorm = 1d-2 ! ni/ne below which the energy QL integral isn't carried out

  LOGICAL, PARAMETER :: traporder1 = .FALSE.

  INTEGER, PARAMETER :: inttype = 2 !1 for CUBART and DQAGSE_QLK, 2 for NAG proxy, 3 for Gauss-Hermite

  ! For inttype = 1 (CUBATR and DQ)
  LOGICAL, PARAMETER :: restar = .FALSE. !restart parameter for 2D integrals
  INTEGER, PARAMETER :: key = 0 !Algorithm parameter choice 
  INTEGER, PARAMETER :: job = 1 !Algorithm choice
  REAL(KIND=DBL), PARAMETER :: tune = 1. !Internal 
  !  relative accuracies for integrations
  !  REAL(KIND=DBL) , PARAMETER :: relacc1 = 1.0d-3 !default for 1D integrals
  !  REAL(KIND=DBL) , PARAMETER :: relacc2 = 1.0d-2 !default for 2D integrals
  INTEGER, PARAMETER :: limit = 5000 !for DQAGE and DQAGSE_QLK routines

  REAL(KIND=DBL) , PARAMETER :: absacc1 = 0.0d-0 !absolute accuracy for 1D integrals
  REAL(KIND=DBL) , PARAMETER :: abacc = 0.0d-0 !5.0d-1 / 5.0d-1 !absolute accuracy for 2D CUBATR int
  REAL(KIND=DBL) , PARAMETER :: abaccQL1 = 0.0d-0  !higher accuracies for 1D QL integrals
  REAL(KIND=DBL) , PARAMETER :: abaccQL2 = 0.0d-0  !higher accuracies for 2D QL integrals

  !For inttype = 2 (nag proxies)
  INTEGER, PARAMETER :: lw=10000 !Limit of integrand evaluations in 1D integrals
  INTEGER, PARAMETER :: lw2=100000 !Limit of integrand evaluations in 1D integrals in makecoefsgau

  !Parameters for asymmetry module
  REAL(KIND=DBL), PARAMETER :: B=-1d5 !lower limit of root solver
  REAL(KIND=DBL), PARAMETER :: C= 1d5 !upper limit of root solver
  REAL(KIND=DBL), PARAMETER :: RE= 1d-6 !Relative error of root solver
  REAL(KIND=DBL), PARAMETER :: AE= 1d-6 !Absolute error of root solver
  INTEGER, PARAMETER :: ntheta = 64 !resolution of field line coordinate
  INTEGER, PARAMETER :: numecoefs = 13 !number of coefficients in coef (e0-9, <R/Ln>, <n>, in-out asym)
  INTEGER, PARAMETER :: numicoefs = 7 !number of coefficients in coef
  REAL(KIND=DBL), PARAMETER :: epsasym = 1.0d-8 !Relative and absolute accuracy for e# integrations
  REAL(KIND=DBL), PARAMETER :: dtheta = 1.0d-5  !dtheta for numerical differentiation

  !Parameters for Newton solver
  REAL(KIND=DBL) , PARAMETER :: maxnerr = 5.d-3 !5.d-2 is default
  REAL(KIND=DBL) , PARAMETER :: maxnerr2 = 5.d-3 !5.d-3 is default (used when going directly into newton for more precision)
  REAL(KIND=DBL) , PARAMETER :: maxiter = 20
  REAL(KIND=DBL) , PARAMETER :: ndif = 5.d-2  !2.5d-2 is default. For function differentiation

  !Used in dispfuncs
  INTEGER, PARAMETER :: nerr=10 !Number of terms in Z# asymptotic limit expansion
  REAL(KIND=DBL), PARAMETER :: abslim = 50. !Value above which asymptotic limit is carried out
  REAL(KIND=DBL), SAVE, DIMENSION(nerr) :: pduittab, puiss2   

  !Used in mod_fluidsol
  INTEGER, PARAMETER :: ndegpoly = 3
  INTEGER, PARAMETER :: ndegx0 = 2
  
  !Integration parameters
  INTEGER, PARAMETER :: int_method = 2 !Uses pcubature in contour integral
  INTEGER, PARAMETER :: newt_method = 2 !Uses pcubature in root finding
  INTEGER, PARAMETER :: QL_method = 2 !Uses pcubature in flux calculation
  INTEGER, PARAMETER :: fluid_method = 1 !Uses hcubature in fluid solver
  INTEGER, PARAMETER :: newt_conv = 2 !Determines convergence criterion for root finding
  INTEGER, PARAMETER :: int_split = 1 !Determines how the integrals are split in the contour integral and root finder
  INTEGER, PARAMETER :: norm = 2 !Determines norm used in cubature routines

CONTAINS  
  SUBROUTINE init_asym()
    INTEGER :: j
    !Initialize data for asymptotic expansion in dispfuncs
    pduittab(1) = 1.
    puiss2(1) = 2.0
    DO j=2,nerr
       pduittab(j) = pduittab(j-1) * (2*j-1)
       puiss2(j)    = puiss2(j-1)*2.0 
    ENDDO
    DO j=1,nerr
       pduittab(j) = pduittab(j)/puiss2(j)
    ENDDO
  END SUBROUTINE init_asym

END MODULE datcal
