!Calculate integrals
MODULE mod_fonct
  USE kind
  USE datcal
  USE datmat
  USE dispfuncs
  USE calltrapints
  USE callpassints

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calcfonctp( p, nu, omega, fonctp )
    !-----------------------------------------------------------
    ! Calculate the trapped particle integrals
    ! Includes collisions for the electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp

    REAL(KIND=DBL)     :: acc, relerr, cc, dd
    REAL(KIND=DBL), DIMENSION(2) :: a,b
    INTEGER            :: minpts, neval, npts  !output number of integrand evaluations
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    INTEGER :: ifailloc
    REAL(KIND=DBL)    :: rfonctpe, ifonctpe, rfonctpiz, ifonctpiz

    REAL(KIND=DBL), DIMENSION(nf) :: intout

    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu

    !Set integration limits. (1) for kappa and (2) for v (for the electrons)
    a(1) = 0.0d0
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    ! The kappa integral is calculated. The real component with rFFkiz,
    ! and the imaginary component with iFFkiz, for all ions

    ifailloc = 1
    rfonctpiz = d01ahf(cc,dd,relacc1,npts,relerr,rFFkiz,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 1DNAG rFFkiz integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF
    ifailloc = 1
    ifonctpiz = d01ahf(cc,dd,relacc1,npts,relerr,iFFkiz,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 1DNAG iFFkiz integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

    minpts = 0; ifailloc=-1;
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
          minpts=0; ifailloc = 1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFke,relacc2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 2DNAG rFFke integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
          minpts=0; ifailloc = 1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFke,relacc2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 2DNAG iFFke integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
       ELSE !Collisionless simulation, revert to faster single integral
          ifailloc = 1
          intout(1) = d01ahf(cc,dd,relacc1,npts,relerr,rFFke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 1DNAG rFFke_nocoll integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
          ifailloc = 1
          intout(2) = d01ahf(cc,dd,relacc1,npts,relerr,iFFke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 1DNAG iFFke_nocoll integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF

    fonctp = intout(1) + rfonctpiz + ci * intout(2) + ci * ifonctpiz

  END SUBROUTINE calcfonctp

  SUBROUTINE calcfonctc( p, nu, omega, fonctc )
    !-----------------------------------------------------------
    ! Calculate the passing particle integrands
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc

    REAL(KIND=DBL), DIMENSION(ndim) :: a, b, c, dd , xy
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    REAL(KIND=DBL), DIMENSION(nf) :: intout

    REAL(KIND=DBL)     :: acc
    INTEGER            :: npts, minpts, neval, iii,jjj,rc !, lims=100
    INTEGER :: ifailloc

    REAL(KIND=DBL)     :: rfonctc, ifonctc !outputs of integrals

    omFkr = omega
    pFkr = p
    nuFkr = nu

    a(1) = 0.0d0 
    a(2) = 0.0d0 
    !Note the 0 lower limit, different from callpassQLints. 
    !This leads to a factor 4 (even symmetry and 2*2) prefactor for the integrals here
    b(:) = rkuplim

    !The k* and rho* double integral is calculated
    !The real and imaginary parts are calculated separately due to NAG function limitations
    ccount=0
    minpts = 0

    minpts=0; ifailloc = 1
    CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstar,relacc2,acc,lenwrk,wrkstr,intout(1),ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 2DNAG rFkstarrstar integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
       IF (ifailloc == -399) THEN
          WRITE(stderr,"(A)") 'NAG license error! Exiting'
          STOP
       ENDIF
    ENDIF

    minpts=0; ifailloc = 1
    CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstar,relacc2,acc,lenwrk,wrkstr,intout(2),ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 2DNAG iFkstarrstar integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

    fonctc = intout(1) + ci * intout(2)

  END SUBROUTINE calcfonctc

  !*********************************************************************************************************************************************
  ! same functional calculations but with rotation
  !*********************************************************************************************************************************************

  SUBROUTINE calcfonctrotp( p, nu, omega, fonctp )
    !-----------------------------------------------------------
    ! Calculate the trapped particle integrals
    ! Includes collisions for the electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp

    REAL(KIND=DBL)     :: acc, relerr, cc, dd
    REAL(KIND=DBL), DIMENSION(2) :: a,b
    INTEGER            :: minpts, neval, npts !output number of integrand evaluations
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    INTEGER :: ifailloc
    REAL(KIND=DBL)    :: rfonctpe, ifonctpe, rfonctpiz, ifonctpiz

    REAL(KIND=DBL), DIMENSION(nf) :: intout

    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu

    !Set integration limits. (1) for kappa and (2) for v (for the electrons)
    a(1) = 0.0d0
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    ! The kappa integral is calculated. The real component with rFFkizrot,
    ! and the imaginary component with iFFkizrot, for all ions

    ifailloc = 1
    rfonctpiz = d01ahf(cc,dd,relacc1,npts,relerr,rFFkizrot,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 1DNAG rFFkizrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF
    ifailloc = 1
    ifonctpiz = d01ahf(cc,dd,relacc1,npts,relerr,iFFkizrot,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 1DNAG iFFkizrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
          minpts=0; ifailloc = 1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkerot,relacc2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 2DNAG rFFkerot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
          minpts=0; ifailloc = 1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkerot,relacc2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 2DNAG iFFkerot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
       ELSE !collisionless integral, do single integral
          ifailloc = 1
          intout(1) = d01ahf(cc,dd,relacc1,npts,relerr,rFFke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 1DNAG rFFke_nocollrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
          ifailloc = 1
          intout(2) = d01ahf(cc,dd,relacc1,npts,relerr,iFFke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                  &'. Abnormal termination of 1DNAG iFFke_nocollrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
          ENDIF
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF

    fonctp = intout(1) + rfonctpiz + ci * intout(2) + ci * ifonctpiz

  END SUBROUTINE calcfonctrotp


  !**************************************************************************************************************


  SUBROUTINE calcfonctrotc( p, nu, omega, fonctc )
    !-----------------------------------------------------------
    ! Calculate the passing particle integrands
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    REAL(KIND=DBL), DIMENSION(ndim) :: a, b, c, dd , xy 
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    REAL(KIND=DBL)     :: acc
    INTEGER            :: maxpts2, npts, minpts, neval, iii,jjj,rc,lenwrk2 
    INTEGER :: ifailloc
    REAL(KIND=DBL)     :: rfonctc, ifonctc !outputs of integrals

    omFkr = omega
    pFkr = p
    nuFkr = nu

    a(:) = -rkuplim ! for passing with rotation need to integrate from -inf to +inf 
    !hence get rid of factor 2x2 in front of integrands in call_passints, also inside QLKfluxes change intmult from 4 to 1
    b(:) = rkuplim

    !The k* and rho* double integral is calculated
    !The real and imaginary parts are calculated separately due to NAG function limitations
    ccount=0
    minpts = 0

    minpts=0; ifailloc = 1
    CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarrot,relacc2,acc,lenwrk,wrkstr,intout(1),ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 2DNAG rFkstarrstarrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
       IF (ifailloc == -399) THEN
          WRITE(stderr,"(A)") 'NAG license error! Exiting'
          STOP
       ENDIF
    ENDIF

    minpts=0; ifailloc = 1  
    CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarrot,relacc2,acc,lenwrk,wrkstr,intout(2),ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of 2DNAG iFkstarrstarrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    ENDIF

    fonctc = intout(1) + ci * intout(2)

  END SUBROUTINE calcfonctrotc

END MODULE mod_fonct
