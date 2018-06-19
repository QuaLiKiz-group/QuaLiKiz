PROGRAM integration_driver

    USE kind
    USE datcal
    USE datmat
    USE integration_routines
    USE callpassints
    USE calltrapints

    IMPLICIT NONE

    INTEGER  :: p, nu
    COMPLEX(KIND=DBL) :: omega

    REAL(KIND=DBL),  DIMENSION(:), ALLOCATABLE :: aa, bb
    REAL(KIND=DBL),  DIMENSION(:), ALLOCATABLE :: wrkstr
    REAL(KIND=DBL),  DIMENSION(:), ALLOCATABLE :: intout
    !REAL(KIND=DBL), DIMENSION(ndim) :: a, b
    !REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    !REAL(KIND=DBL), DIMENSION(nf) :: intout

    REAL(KIND=DBL)     :: acc
    INTEGER            :: minpts
    !REAL(KIND=DBL)     :: relacc2
    !INTEGER            :: maxpts, lenwrk, ndim !Defined somewhere else
    !LOGICAL            :: verbose

    PRINT *, "Hello World!"
    ! For generating the debug printout
    p = 0
    nu = 0
    omega = 0
    verbose = .TRUE.

    ! Integation hyperpars
    minpts = 0
    maxpts = 5e5
    relacc2 = 0.02
    !ndim = 2

    ! Routine temp arrays
    !lenwrk = 1e5

    ALLOCATE(aa(ndim))
    ALLOCATE(bb(ndim))
    ALLOCATE(wrkstr(lenwrk))
    ALLOCATE(intout(2))

    ! Integration bounds
    aa(1) = 0.0d0
    aa(2) = 0.0d0
    bb(:) = 5.d0 * 1.


    intout(:) = 0
    ! iFkstarrstare(nf, ny) calls Fkstarrstare(ndim,xy,1)
    ! Fkstarrstare(ndim, xx, caseflag) depends on:
    !       INTEGER, INTENT(IN) :: ndim
    !       REAL(KIND=DBL)   , INTENT(IN) :: xx(ndim)
    !       INTEGER , INTENT(IN) :: caseflag

    !       REAL(KIND=DBL)    :: Ather⋅
    !       COMPLEX(KIND=DBL) :: aae, bbe, cce, dde, sqrtde, Vme, Vpe, Zae
    !       COMPLEX(KIND=DBL) :: alphae
    !       COMPLEX(KIND=DBL) :: face
    !       REAL(KIND=DBL)    :: nwge,var2,var3,bessm2
    !       REAL(KIND=DBL)    :: rstar, kstar, teta, fkstar
    !       COMPLEX(KIND=DBL) :: inte3, inte5
    !       COMPLEX(KIND=DBL) :: Fekstarrstar
    !       COMPLEX(KIND=DBL) :: Febkstarrstar

    !!!!!!!!!!!!!!!!!
    ! List 'o hidden dependencies
    !!!!!!!!!!!!!!!!!
    ! input:
    ! - smag(pFkr)
    ! - alphax(pFkr)
    ! - Ane(pFkr)
    ! - Ate(pFkr)
    ! - Nex(pFkr)
    ! - Tex(pFkr)
    ! - Rhoe(pFkr)

    ! calcroutines
    ! - mwidth = modewidth(p,nu) [calcroutines.f90:154]
    !   - modewidth(p,nu) = jon_modewidth(p,nu) [calcroutines.f90:153]
    !     - jon_modewidth(p,nu) = mwidth [calcroutines.f90:137]
    ! - nwg = ntor(p,nu)*wg(p)
    !   - ntor(p,nu) = kthetarhos(nu)*x(p)/(qx(p)*rhostar(p)) [mod_make_io]
    !     - rhostar(:) = 1./Rmin(:)*SQRT(Tex(:)*1.d3*qe/mi(:,1))/(Zi(:,1)*qe*Bo(:)/mi(:,1)) !With respect to main ion
    !   - wg(:) = qx(:)*1.d3/(Ro(:)*Bo(:)*Rmin(:)*(x(:)+eps)) [mod_make_io]
    ! - ktetaRhoe = kteta*Rhoe(p)
    !   - kteta     = ntor(p,nu)*ktetasn(p)
    !     - ktetasn(:) = qx(:)/(Rmin(:)*(x(:)+eps)) [mod_make_io]
    !   - Rhoe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3*me)/(qe*Bo(:)) [mod_make_io]
    ! - Athe=widthhat*cthe(p)/qRd
    !   - widthhat = ABS(mwidth)**2 / SQRT(REAL(mwidth**2)) [calcroutines]
    !   - cthe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3/me) [mod_make_io]
    !   - qRd = qx(p)*Ro(p)*d*1./alam1 [calcroutines]
    !     - alam1=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam1int,lw,ifailloc)/alamnorm !pitch angle average of sqrt(1-lambda*b) [calcroutines]
    !       - alam1int = Rlam2/Rlam1*1./(4.*pi)*Tlam [calcroutines]
    !         - CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,5)
    !           - Rint2 = Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**something !something depends on where it is called
    !         - CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,4)
    !           - Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    !             - lamin ! input of function
    !             - plam = p
    !             - theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    !         - CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,3)
    !           - Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))
    !       - alamnorm = fc(p) !to be consistent with passing particle fraction

    ! Constants of nature:
    ! - minfki = 1d-3 (datcal)
    ! - Ze = -1._DBL
    ! - twopi   = 6.28318530718
    ! - epsD=1.d-14

    ! Always the case all over the place
    ! - omFkr = omega ! Just a placeholder to store the result in?

    ! Fried-Conte functions (globals)
    ! - Z1(aae)
    ! - Z1(Vpe)
    ! - Z2
    ! - Z3

    ! Some other guys:
    ! - d ?= ABS(1./(ntor(p,nu)*(epsD+ABS(qprim(p))))) ! distance between rational surfaces
    !   - qprim(:)=smag(:)*qx(:)/(Rmin(:)*x(:)) [mod_make_io]
    ! - fc(pFkr)
    !   - fc(:) = 1. - ft(:) !passing particle fraction
    !     - ft(:) = 2.*(2.*epsilon(:))**0.5/pi; !trapped particle fractio

    intout(1) = integrate_2d(ndim,aa,bb,minpts,maxpts,iFkstarrstar,relacc2,acc,lenwrk,wrkstr,verbose,p,nu,omega,'iFkstarrstar')
    PRINT *, intout(1)
END PROGRAM integration_driver
