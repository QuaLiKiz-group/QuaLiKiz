MODULE calltrapQLints

  USE kind
  USE datcal
  USE datmat
  USE calltrapints

  IMPLICIT NONE

CONTAINS

  SUBROUTINE trapQLints( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne,&
       fonctpgni, fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni, fonctepce, fonctepci)
    !-----------------------------------------------------------
    ! Trapped particle integrals at the omega of the linear solutions
    ! Double integral done for trapped electrons due to collisions
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! To save some computational time, only the imaginary component is kept - relevant for transport
    !----------------------------------------------------------- 
    INTEGER, INTENT(IN)  :: p,nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctpe, fonctpgte, fonctpgne, fonctpce, fonctepe
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctepgte, fonctepgne, fonctepce
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctpi, fonctpgti, fonctpgni, fonctpci, fonctepi
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctepgti, fonctepgni, fonctepci

    INTEGER, PARAMETER :: lwloc=50000 !default 10000
    INTEGER :: npts, minpts, neval

    REAL(KIND=DBL), DIMENSION(ndim) :: a,b
    REAL(KIND=DBL) :: acc, cc, dd, relerr 
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL)    :: rfonctpe, rfonctpgte, rfonctpgne, rfonctpce, rfonctepe
    REAL(KIND=DBL)    :: rfonctepgte, rfonctepgne, rfonctepce
    REAL(KIND=DBL)    :: ifonctpe, ifonctpgte, ifonctpgne, ifonctpce, ifonctepe
    REAL(KIND=DBL)    :: ifonctepgte, ifonctepgne, ifonctepce

    REAL(KIND=DBL), DIMENSION(nions) :: rfonctpi, rfonctpgti, rfonctpgni, rfonctpci, rfonctepi
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctepgti, rfonctepgni, rfonctepci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctpi, ifonctpgti, ifonctpgni, ifonctpci, ifonctepi
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctepgti, ifonctepgni, ifonctepci

    REAL(KIND=DBL), DIMENSION(nf) :: intout
    INTEGER :: ifailloc

    omFFk = omega
    pFFk = p

    !! Integration bounds
    !! a,b are for the 1D ions. 
    a(1) = 0.0d0
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    DO ion=1,nions

       !! ION PARTICLE FLUX 
       !ifailloc=1
       !rfonctpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFki,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,*) 'Abnormal termination of rFFki integration at p= ',p,', ion=',ion
       !ENDIF
       ifailloc=1     
       ifonctpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFki,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFki integration at p=',p,', nu=',nu,', ion=',ion
       ENDIF
       rfonctpi(ion) = 0.

       ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
       IF (phys_meth .NE. 0.0) THEN

          !ifailloc=1
          !rfonctpgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgti,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgti integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgti,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgti integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpgti(ion)=0.

          !ifailloc=1
          !rfonctpgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgni,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgni integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgni,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgni integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpgni(ion)=0.

          !ifailloc=1
          !rfonctpci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkci,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkci integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkci,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkci integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpci(ion)=0.
!!!
          IF (phys_meth == 2) THEN
             !ifailloc=1
             !rfonctepgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgti,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgti integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgti,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgti integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepgti(ion) = 0.
             ENDIF
             rfonctepgti(ion)=0.

             !ifailloc=1
             !rfonctepgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgni,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgni integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgni,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgni integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepgni(ion) = 0.
             ENDIF
             rfonctepgni(ion)=0.

             !ifailloc=1
             !rfonctepci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekci,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekci integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekci,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekci integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepci(ion) = 0.
             ENDIF
             rfonctepci(ion)=0.
          ELSE
             rfonctepgti(ion) = 0.
             ifonctepgti(ion) = 0.
             rfonctepgni(ion) = 0.
             ifonctepgni(ion) = 0.
             rfonctepci(ion) = 0.
             ifonctepci(ion) = 0.
          ENDIF
       ELSE
          rfonctpgti(ion) = 0.
          ifonctpgti(ion) = 0.
          rfonctpgni(ion) = 0.
          ifonctpgni(ion) = 0.
          rfonctpci(ion) = 0.
          ifonctpci(ion) = 0.
          rfonctepgti(ion) = 0.
          ifonctepgti(ion) = 0.
          rfonctepgni(ion) = 0.
          ifonctepgni(ion) = 0.
          rfonctepci(ion) = 0.
          ifonctepci(ion) = 0.

       ENDIF

       ! ION ENERGY FLUX

       !ifailloc=1
       !rfonctepi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFeki,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFeki integration at p=',p,', nu=',nu,', ion=',ion
       !ENDIF

       ifailloc=1
       IF (ninorm(p,ion) > min_ninorm) THEN
          ifonctepi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFeki,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFeki integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
       ELSE
          ifonctepi(ion) = 0.
       ENDIF
       rfonctepi(ion)=0.
    ENDDO


    ! 2D INTEGRALS FOR ELECTRONS

    minpts = 0; ifailloc=1;
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 

       IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFke,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFke integration at p=',p,' nu=',nu
          ENDIF

          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFke,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFke integration at p=',p,' nu=',nu
          ENDIF
       ELSE ! Collisionless simulation, revert to faster single integral
          ifailloc=1
          intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFke_nocoll integration at p=',p,' nu=',nu
          ENDIF

          ifailloc=1
          intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFke_nocoll integration at p=',p,' nu=',nu
          ENDIF
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF
    rfonctpe = intout(1)
    ifonctpe = intout(2)

    !   !!ADDITIONAL ELECTRON INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
    IF (phys_meth .NE. 0.0) THEN

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkgte,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkgte integration at p=',p,' nu=',nu
             ENDIF

             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkgte,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkgte integration at p=',p,' nu=',nu
             ENDIF

          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgte_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgte_nocoll integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgte_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgte_nocoll integration at p=',p,' nu=',nu
             ENDIF

          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpgte = intout(1)
       ifonctpgte = intout(2)

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkgne,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkgne integration at p=',p,' nu=',nu
             ENDIF
             minpts=0;  ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkgne,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkgne integration at p=',p,' nu=',nu
             ENDIF
          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgne_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgne_nocoll integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgne_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgne_nocoll integration at p=',p,' nu=',nu
             ENDIF
          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpgne = intout(1)
       ifonctpgne = intout(2)

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkce,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkce integration at p=',p,' nu=',nu
             ENDIF
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkce,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkce integration at p=',p,' nu=',nu
             ENDIF
          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkce_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkce_nocoll integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkce_nocoll,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkce_nocoll integration at p=',p,' nu=',nu
             ENDIF
          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpce = intout(1)
       ifonctpce = intout(2)
!!!
       IF (phys_meth == 2) THEN
          minpts = 0; ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekgte,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekgte integration at p=',p,' nu=',nu
                ENDIF

                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekgte,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekgte integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgte_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgte_nocoll integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgte_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgte_nocoll integration at p=',p,' nu=',nu
                ENDIF

             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepgte = intout(1)
          ifonctepgte = intout(2)

          minpts = 0; ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekgne,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekgne integration at p=',p,' nu=',nu
                ENDIF
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekgne,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekgne integration at p=',p,' nu=',nu
                ENDIF
             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgne_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgne_nocoll integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgne_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgne_nocoll integration at p=',p,' nu=',nu
                ENDIF
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepgne = intout(1)
          ifonctepgne = intout(2)

          minpts = 0; ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekce,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekce integration at p=',p,' nu=',nu
                ENDIF
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekce,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekce integration at p=',p,' nu=',nu
                ENDIF
             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekce_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekce_nocoll integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekce_nocoll,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekce_nocoll integration at p=',p,' nu=',nu
                ENDIF
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepce = intout(1)
          ifonctepce = intout(2)
       ELSE
          rfonctepgte = 0.
          ifonctepgte = 0.
          rfonctepgne = 0.
          ifonctepgne = 0.
          rfonctepce = 0.
          ifonctepce = 0.
       ENDIF
    ELSE
       rfonctpgte = 0.
       ifonctpgte = 0.
       rfonctpgne = 0.
       ifonctpgne = 0.
       rfonctpce = 0.
       ifonctpce = 0.
       rfonctepgte = 0.
       ifonctepgte = 0.
       rfonctepgne = 0.
       ifonctepgne = 0.
       rfonctepce = 0.
       ifonctepce = 0.
    ENDIF

    !ELECTRON ENERGY FLUX

    minpts = 0; ifailloc=1;
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFeke,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFeke integration at p=',p,' nu=',nu
          ENDIF
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFeke,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFeke integration at p=',p,' nu=',nu
          ENDIF

       ELSE ! Collisionless simulation, revert to faster single integral
          ifailloc=1
          intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFeke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFeke_nocoll integration at p=',p,' nu=',nu
          ENDIF
          ifailloc=1
          intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFeke_nocoll,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFeke_nocoll integration at p=',p,' nu=',nu
          ENDIF
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF
    rfonctepe = intout(1)
    ifonctepe = intout(2)

    !The complex forms are reconstructed
    fonctpe = rfonctpe + ci * ifonctpe
    fonctpi(:) = rfonctpi(:) + ci * ifonctpi(:)

    fonctpgte = rfonctpgte + ci * ifonctpgte
    fonctpgti(:) = rfonctpgti(:) + ci * ifonctpgti(:)

    fonctpgne = rfonctpgne + ci * ifonctpgne
    fonctpgni(:) = rfonctpgni(:) + ci * ifonctpgni(:)

    fonctpce = rfonctpce + ci * ifonctpce
    fonctpci(:) = rfonctpci(:) + ci * ifonctpci(:)

    fonctepe = rfonctepe + ci * ifonctepe
    fonctepi(:) = rfonctepi(:) + ci * ifonctepi(:)

    fonctepgte = rfonctepgte + ci * ifonctepgte
    fonctepgti(:) = rfonctepgti(:) + ci * ifonctepgti(:)

    fonctepgne = rfonctepgne + ci * ifonctepgne
    fonctepgni(:) = rfonctepgni(:) + ci * ifonctepgni(:)

    fonctepce = rfonctepce + ci * ifonctepce
    fonctepci(:) = rfonctepci(:) + ci * ifonctepci(:)


  END SUBROUTINE trapQLints

  
!**************************************************************************************************************************************
! with rotation
!**************************************************************************************************************************************

  SUBROUTINE trapQLintsrot( p, nu, omega, fonctpe, fonctpi, fonctpgte, fonctpgti, fonctpgne, fonctpgni, fonctpgui, &
       & fonctpce, fonctpci, fonctepe, fonctepi, fonctepgte, fonctepgti, fonctepgne, fonctepgni, fonctepgui, fonctepce, fonctepci, fonctvpi)
    !-----------------------------------------------------------
    ! Trapped particle integrals at the omega of the linear solutions
    ! Double integral done for trapped electrons due to collisions
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! To save some computational time, only the imaginary component is kept - relevant for transport
    !----------------------------------------------------------- 
    INTEGER, INTENT(IN)  :: p,nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctpe, fonctpgte, fonctpgne, fonctpce, fonctepe
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctepgte, fonctepgne, fonctepce
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctpi, fonctpgti, fonctpgni, fonctpgui, fonctpci, fonctepi, fonctvpi
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctepgti, fonctepgni, fonctepgui, fonctepci

    INTEGER, PARAMETER :: lwloc=50000 !default 10000
    INTEGER :: npts, minpts, neval

    REAL(KIND=DBL), DIMENSION(ndim) :: a,b
    REAL(KIND=DBL) :: acc, cc, dd, relerr 
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL)    :: rfonctpe, rfonctpgte, rfonctpgne, rfonctpce, rfonctepe
    REAL(KIND=DBL)    :: rfonctepgte, rfonctepgne, rfonctepce
    REAL(KIND=DBL)    :: ifonctpe, ifonctpgte, ifonctpgne, ifonctpce, ifonctepe
    REAL(KIND=DBL)    :: ifonctepgte, ifonctepgne, ifonctepce

    REAL(KIND=DBL), DIMENSION(nions) :: rfonctpi, rfonctpgti, rfonctpgni, rfonctpgui, rfonctpci, rfonctepi, rfonctvpi
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctepgti, rfonctepgni, rfonctepgui, rfonctepci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctpi, ifonctpgti, ifonctpgni, ifonctpgui, ifonctpci, ifonctepi, ifonctvpi
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctepgti, ifonctepgni, ifonctepgui, ifonctepci

    REAL(KIND=DBL), DIMENSION(nf) :: intout
    INTEGER :: ifailloc

    omFFk = omega
    pFFk = p

    !! Integration bounds
    !! a,b are for the 1D ions.
    a(1) = 0.0d0
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    DO ion=1,nions

       !! ION PARTICLE FLUX 
       !ifailloc=1
       !rfonctpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkirot,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,*) 'Abnormal termination of rFFkirot integration at p= ',p,', ion=',ion
       !ENDIF
       ifailloc=1     
       ifonctpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkirot,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkirot integration at p=',p,', nu=',nu,', ion=',ion
       ENDIF
       rfonctpi(ion) = 0.

       ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
       IF (phys_meth .NE. 0.0) THEN

          !ifailloc=1
          !rfonctpgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgtirot,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgtirot integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgtirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgtirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpgti(ion)=0.

          !ifailloc=1
          !rfonctpgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgnirot,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgnirot integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgnirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgnirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpgni(ion)=0.

          !ifailloc=1
          !rfonctpgui(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkguirot,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkguirot integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpgui(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkguirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkguirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpgui(ion)=0.

          !ifailloc=1
          !rfonctpci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkcirot,lw,ifailloc)
          !IF (ifailloc /= 0) THEN
          !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkcirot integration at p=',p,', nu=',nu,', ion=',ion
          !ENDIF

          ifailloc=1
          ifonctpci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkcirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkcirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpci(ion)=0.
!!!
          IF (phys_meth == 2) THEN
             !ifailloc=1
             !rfonctepgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgtirot,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgtirot integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepgti(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgtirot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgtirot integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepgti(ion) = 0.
             ENDIF
             rfonctepgti(ion)=0.

             !ifailloc=1
             !rfonctepgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgnirot,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgnirot integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepgni(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgnirot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgnirot integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepgni(ion) = 0.
             ENDIF
             rfonctepgni(ion)=0.

             !ifailloc=1
             !rfonctepgui(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekguirot,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekguirot integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepgui(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekguirot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekguirot integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepgui(ion) = 0.
             ENDIF
             rfonctepgui(ion)=0.

             !ifailloc=1
             !rfonctepci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekcirot,lw,ifailloc)
             !IF (ifailloc /= 0) THEN
             !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekcirot integration at p=',p,', nu=',nu,', ion=',ion
             !ENDIF

             ifailloc=1
             IF (ninorm(p,ion) > min_ninorm) THEN
                ifonctepci(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekcirot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekcirot integration at p=',p,', nu=',nu,', ion=',ion
                ENDIF
             ELSE
                ifonctepci(ion) = 0.
             ENDIF
             rfonctepci(ion)=0.
          ELSE
             rfonctepgti(ion) = 0.
             ifonctepgti(ion) = 0.
             rfonctepgni(ion) = 0.
             ifonctepgni(ion) = 0.
             rfonctepgui(ion) = 0.
             ifonctepgui(ion) = 0.
             rfonctepci(ion) = 0.
             ifonctepci(ion) = 0.
          ENDIF
       ELSE
          rfonctpgti(ion) = 0.
          ifonctpgti(ion) = 0.
          rfonctpgni(ion) = 0.
          ifonctpgni(ion) = 0.
          rfonctpgui(ion) = 0.
          ifonctpgui(ion) = 0.
          rfonctpci(ion) = 0.
          ifonctpci(ion) = 0.
          rfonctepgti(ion) = 0.
          ifonctepgti(ion) = 0.
          rfonctepgni(ion) = 0.
          ifonctepgni(ion) = 0.
          rfonctepgui(ion) = 0.
          ifonctepgui(ion) = 0.
          rfonctepci(ion) = 0.
          ifonctepci(ion) = 0.
       ENDIF

       ! ION ENERGY FLUX

       !ifailloc=1
       !rfonctepi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekirot,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekirot integration at p=',p,', nu=',nu,', ion=',ion
       !ENDIF

       ifailloc=1
       IF (ninorm(p,ion) > min_ninorm) THEN
          ifonctepi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
       ELSE
          ifonctepi(ion) = 0.
       ENDIF
       rfonctepi(ion)=0.

       ! ION ang mom FLUX

       !ifailloc=1
       !rfonctvpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFvkirot,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
       !ENDIF

       ifailloc=1
       IF (ninorm(p,ion) > min_ninorm) THEN
          ifonctvpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFvkirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
       ELSE
          ifonctvpi(ion) = 0.
       ENDIF
       rfonctvpi(ion)=0.

    ENDDO


    ! 2D INTEGRALS FOR ELECTRONS

    minpts = 0; ifailloc=1;
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 

       IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkerot integration at p=',p,' nu=',nu
          ENDIF

          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkerot integration at p=',p,' nu=',nu
          ENDIF
       ELSE ! Collisionless simulation, revert to faster single integral
          ifailloc=1
          intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFke_nocollrot integration at p=',p,' nu=',nu
          ENDIF

          ifailloc=1
          intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFke_nocollrot integration at p=',p,' nu=',nu
          ENDIF
       ENDIF

    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF
    rfonctpe = intout(1)
    ifonctpe = intout(2)

    !   !!ADDITIONAL ELECTRON INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
    IF (phys_meth .NE. 0.0) THEN

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkgterot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkgterot integration at p=',p,' nu=',nu
             ENDIF

             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkgterot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkgterot integration at p=',p,' nu=',nu
             ENDIF
          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgte_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgte_nocollrot integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgte_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgte_nocollrot integration at p=',p,' nu=',nu
             ENDIF

          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpgte = intout(1)
       ifonctpgte = intout(2)

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkgnerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkgnerot integration at p=',p,' nu=',nu
             ENDIF
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkgnerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkgnerot integration at p=',p,' nu=',nu
             ENDIF
          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkgne_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkgne_nocollrot integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkgne_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkgne_nocollrot integration at p=',p,' nu=',nu
             ENDIF
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpgne = intout(1)
       ifonctpgne = intout(2)

       minpts = 0; ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,rFFkcerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFkcerot integration at p=',p,' nu=',nu
             ENDIF
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFFkcerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFkcerot integration at p=',p,' nu=',nu
             ENDIF
          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFkce_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFkce_nocollrot integration at p=',p,' nu=',nu
             ENDIF
             ifailloc=1
             intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFkce_nocollrot,lw,ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFkce_nocollrot integration at p=',p,' nu=',nu
             ENDIF
          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctpce = intout(1)
       ifonctpce = intout(2)
!!!
       IF (phys_meth == 2) THEN
          minpts = 0; ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekgterot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekgterot integration at p=',p,' nu=',nu
                ENDIF

                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekgterot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekgterot integration at p=',p,' nu=',nu
                ENDIF
             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgte_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgte_nocollrot integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgte_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgte_nocollrot integration at p=',p,' nu=',nu
                ENDIF

             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepgte = intout(1)
          ifonctepgte = intout(2)

          minpts = 0; ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekgnerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekgnerot integration at p=',p,' nu=',nu
                ENDIF
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekgnerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekgnerot integration at p=',p,' nu=',nu
                ENDIF
             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekgne_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekgne_nocollrot integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekgne_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekgne_nocollrot integration at p=',p,' nu=',nu
                ENDIF
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepgne = intout(1)
          ifonctepgne = intout(2)

          minpts = 0; ifailloc=1; 
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekcerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekcerot integration at p=',p,' nu=',nu
                ENDIF
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekcerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekcerot integration at p=',p,' nu=',nu
                ENDIF
             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFekce_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFekce_nocollrot integration at p=',p,' nu=',nu
                ENDIF
                ifailloc=1
                intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFekce_nocollrot,lw,ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFekce_nocollrot integration at p=',p,' nu=',nu
                ENDIF
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctepce = intout(1)
          ifonctepce = intout(2)
       ELSE
          rfonctepgte = 0.
          ifonctepgte = 0.
          rfonctepgne = 0.
          ifonctepgne = 0.
          rfonctepce = 0.
          ifonctepce = 0.
       ENDIF
    ELSE
       rfonctpgte = 0.
       ifonctpgte = 0.
       rfonctpgne = 0.
       ifonctpgne = 0.
       rfonctpce = 0.
       ifonctpce = 0.
       rfonctepgte = 0.
       ifonctepgte = 0.
       rfonctepgne = 0.
       ifonctepgne = 0.
       rfonctepce = 0.
       ifonctepce = 0.
    ENDIF

    !ELECTRON ENERGY FLUX

    minpts = 0; ifailloc=1;
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,rFFekerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFFekerot integration at p=',p,' nu=',nu
          ENDIF
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFFekerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFFekerot integration at p=',p,' nu=',nu
          ENDIF
       ELSE ! Collisionless simulation, revert to faster single integral
          ifailloc=1
          intout(1) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFeke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFeke_nocollrot integration at p=',p,' nu=',nu
          ENDIF
          ifailloc=1
          intout(2) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFeke_nocollrot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFeke_nocollrot integration at p=',p,' nu=',nu
          ENDIF
       ENDIF

    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF
    rfonctepe = intout(1)
    ifonctepe = intout(2)

    !The complex forms are reconstructed
    fonctpe = rfonctpe + ci * ifonctpe
    fonctpi(:) = rfonctpi(:) + ci * ifonctpi(:)

    fonctpgte = rfonctpgte + ci * ifonctpgte
    fonctpgti(:) = rfonctpgti(:) + ci * ifonctpgti(:)

    fonctpgne = rfonctpgne + ci * ifonctpgne
    fonctpgni(:) = rfonctpgni(:) + ci * ifonctpgni(:)

    fonctpgui(:) = rfonctpgui(:) + ci * ifonctpgui(:)

    fonctpce = rfonctpce + ci * ifonctpce
    fonctpci(:) = rfonctpci(:) + ci * ifonctpci(:)

    fonctepe = rfonctepe + ci * ifonctepe
    fonctepi(:) = rfonctepi(:) + ci * ifonctepi(:)

    fonctvpi(:) = rfonctvpi(:) + ci * ifonctvpi(:)

    fonctepgte = rfonctepgte + ci * ifonctepgte
    fonctepgti(:) = rfonctepgti(:) + ci * ifonctepgti(:)

    fonctepgne = rfonctepgne + ci * ifonctepgne
    fonctepgni(:) = rfonctepgni(:) + ci * ifonctepgni(:)

    fonctepce = rfonctepce + ci * ifonctepce
    fonctepci(:) = rfonctepci(:) + ci * ifonctepci(:)

    fonctepgui(:) = rfonctepgui(:) + ci * ifonctepgui(:)


  END SUBROUTINE trapQLintsrot

  SUBROUTINE momtrapQLintsrot( p, nu, omega, fonctvpi)
    !-----------------------------------------------------------
    ! Trapped particle integrals at the omega of the linear solutions. ONLY FOR ANGULAR MOMENTUM
    ! Double integral done for trapped electrons due to collisions
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! To save some computational time, only the imaginary component is kept - relevant for transport
    !----------------------------------------------------------- 
    INTEGER, INTENT(IN)  :: p,nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctvpi

    INTEGER, PARAMETER :: lwloc=50000 !default 10000
    INTEGER :: npts, minpts, neval

    REAL(KIND=DBL), DIMENSION(ndim) :: a,b
    REAL(KIND=DBL) :: cc, dd, relerr 
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL), DIMENSION(nions) :: rfonctvpi
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctvpi

    REAL(KIND=DBL), DIMENSION(nf) :: intout
    INTEGER :: ifailloc

    omFFk = omega
    pFFk = p

    !! Integration bounds
    !! a,b are for the 1D ions
    a(1) = 0.0d0
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    DO ion=1,nions

       ! ION ang mom FLUX

       !ifailloc=1
       !rfonctvpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,rFFvkirot,lw,ifailloc)
       !IF (ifailloc /= 0) THEN
       !   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL rFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
       !ENDIF

       ifailloc=1
       IF (ninorm(p,ion) > min_ninorm) THEN
          ifonctvpi(ion) = d01ahf(cc,dd,relaccQL1,npts,relerr,iFFvkirot,lw,ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of 1DNAG QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
       ELSE
          ifonctvpi(ion) = 0.
       ENDIF
       rfonctvpi(ion)=0.

    ENDDO

    !The complex forms are reconstructed

    fonctvpi(:) = rfonctvpi(:) + ci * ifonctvpi(:)


  END SUBROUTINE momtrapQLintsrot


END MODULE calltrapQLints
