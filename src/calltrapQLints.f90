MODULE calltrapQLints

  USE kind
  USE datcal
  USE datmat
  USE calltrapints
  USE HCUB
  USE PCUB

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
    INTEGER :: npts, neval

    REAL(KIND=DBL), DIMENSION(ndim) :: a,b
    REAL(KIND=DBL) :: acc, cc, dd, relerr 
    REAL(KIND=DBL), DIMENSION(1) :: cc_cub, dd_cub
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

    REAL(KIND=DBL), DIMENSION(1) :: intout_cub
    REAL(KIND=DBL), DIMENSION(2) :: acc_cub


    omFFk = omega
    pFFk = p

    !! Integration bounds
    !! a,b are for the 1D ions.
    a(1) = 0.0d0 + barelyavoid
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim

    cc = a(1)
    dd = b(1)

    cc_cub(1) = cc
    dd_cub(1) = dd

    IF(QL_method.EQ.1) THEN
       DO ion=1,nions

          ifailloc=1     
          ifailloc = hcubature(1, iFFki_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          ifonctpi(ion) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFki integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpi(ion) = 0.

          ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
          IF (phys_meth .NE. 0.0) THEN

             ifailloc=1
             ifailloc = hcubature(1, iFFkgti_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpgti(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgti integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpgti(ion)=0.


             ifailloc=1
             ifailloc = hcubature(1, iFFkgni_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpgni(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgni integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpgni(ion)=0.


             ifailloc=1
             ifailloc = hcubature(1, iFFkci_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpci(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkci integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpci(ion)=0.
!!!
             IF (phys_meth == 2) THEN             

                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFFekgti_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepgti(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgti integration at p=',p,', nu=',nu,', ion=',ion
                   ENDIF
                ELSE
                   ifonctepgti(ion) = 0.
                ENDIF
                rfonctepgti(ion)=0.


                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFFekgni_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepgni(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgni integration at p=',p,', nu=',nu,', ion=',ion
                   ENDIF
                ELSE
                   ifonctepgni(ion) = 0.
                ENDIF
                rfonctepgni(ion)=0.


                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFFekci_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepci(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekci integration at p=',p,', nu=',nu,', ion=',ion
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


          ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             ifailloc = hcubature(1, iFFeki_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctepi(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFeki integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
          ELSE
             ifonctepi(ion) = 0.
          ENDIF
          rfonctepi(ion)=0.
       ENDDO


       ! 2D INTEGRALS FOR ELECTRONS

       ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 

          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             ifailloc=1
             ifailloc = hcubature(2, FFke_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFke integration at p=',p,' nu=',nu
             ENDIF

          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc = hcubature(2, FFke_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFke_nocoll integration at p=',p,' nu=',nu
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

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                ifailloc=1
                ifailloc = hcubature(2, FFkgte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkgte integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = hcubature(2, FFkgte_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkgte_nocoll integration at p=',p,' nu=',nu
                ENDIF


             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctpgte = intout(1)
          ifonctpgte = intout(2)

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                ifailloc=1
                ifailloc = hcubature(2, FFkgne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkgne integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = hcubature(2, FFkgne_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkgne_nocoll integration at p=',p,' nu=',nu
                ENDIF

             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctpgne = intout(1)
          ifonctpgne = intout(2)

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                ifailloc=1
                ifailloc = hcubature(2, FFkce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkce integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = hcubature(2, FFkce_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFkce_nocoll integration at p=',p,' nu=',nu
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
             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN 
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                   ifailloc=1
                   ifailloc = hcubature(2, FFekgte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekgte integration at p=',p,' nu=',nu
                   ENDIF



                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = hcubature(2, FFekgte_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekgte_nocoll integration at p=',p,' nu=',nu
                   ENDIF


                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctepgte = intout(1)
             ifonctepgte = intout(2)

             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN 
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                   ifailloc=1
                   ifailloc = hcubature(2, FFekgne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekgne integration at p=',p,' nu=',nu
                   ENDIF

                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = hcubature(2, FFekgne_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekgne_nocoll integration at p=',p,' nu=',nu
                   ENDIF

                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctepgne = intout(1)
             ifonctepgne = intout(2)

             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                   ifailloc=1
                   ifailloc = hcubature(2, FFekce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekce integration at p=',p,' nu=',nu
                   ENDIF

                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = hcubature(2, FFekce_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFekce_nocoll integration at p=',p,' nu=',nu
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

       
       ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             ifailloc=1
             ifailloc = hcubature(2, FFeke_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFeke integration at p=',p,' nu=',nu
             ENDIF


          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             ifailloc = hcubature(2, FFeke_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL FFeke_nocoll integration at p=',p,' nu=',nu
             ENDIF

          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctepe = intout(1)
       ifonctepe = intout(2)

    ELSE IF(QL_method.EQ.2) THEN
       DO ion=1,nions


          ifailloc=1     
          ifailloc = pcubature(1, iFFki_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          ifonctpi(ion) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFki integration at p=',p,', nu=',nu,', ion=',ion
          ENDIF
          rfonctpi(ion) = 0.

          ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
          IF (phys_meth .NE. 0.0) THEN

             ifailloc=1
             ifailloc = pcubature(1, iFFkgti_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpgti(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgti integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpgti(ion)=0.


             ifailloc=1
             ifailloc = pcubature(1, iFFkgni_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpgni(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgni integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpgni(ion)=0.


             ifailloc=1
             ifailloc = pcubature(1, iFFkci_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctpci(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkci integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
             rfonctpci(ion)=0.
!!!
             IF (phys_meth == 2) THEN             

                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFFekgti_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepgti(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgti integration at p=',p,', nu=',nu,', ion=',ion
                   ENDIF
                ELSE
                   ifonctepgti(ion) = 0.
                ENDIF
                rfonctepgti(ion)=0.


                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFFekgni_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepgni(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgni integration at p=',p,', nu=',nu,', ion=',ion
                   ENDIF
                ELSE
                   ifonctepgni(ion) = 0.
                ENDIF
                rfonctepgni(ion)=0.


                ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFFekci_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   ifonctepci(ion) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekci integration at p=',p,', nu=',nu,', ion=',ion
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


          ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             !This integral is extremely hard for pcubature, so hcubature is used
             ifailloc = hcubature(1, iFFeki_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             ifonctepi(ion) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFeki integration at p=',p,', nu=',nu,', ion=',ion
             ENDIF
          ELSE
             ifonctepi(ion) = 0.
          ENDIF
          rfonctepi(ion)=0.
       ENDDO


       ! 2D INTEGRALS FOR ELECTRONS

       ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 

          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             ifailloc=1
             ifailloc = pcubature(2, FFke_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFke integration at p=',p,' nu=',nu
             ENDIF

          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc = pcubature(2, FFke_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFke_nocoll integration at p=',p,' nu=',nu
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

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                ifailloc=1
                ifailloc = pcubature(2, FFkgte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkgte integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = pcubature(2, FFkgte_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkgte_nocoll integration at p=',p,' nu=',nu
                ENDIF


             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctpgte = intout(1)
          ifonctpgte = intout(2)

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN 
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                ifailloc=1
                ifailloc = pcubature(2, FFkgne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkgne integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = pcubature(2, FFkgne_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkgne_nocoll integration at p=',p,' nu=',nu
                ENDIF

             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctpgne = intout(1)
          ifonctpgne = intout(2)

          ifailloc=1;
          !Only calculate nonadiabatic part for electrons if el_type == 1
          IF (el_type == 1) THEN
             IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                ifailloc=1
                ifailloc = pcubature(2, FFkce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkce integration at p=',p,' nu=',nu
                ENDIF

             ELSE ! Collisionless simulation, revert to faster single integral
                ifailloc=1
                ifailloc = pcubature(2, FFkce_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFkce_nocoll integration at p=',p,' nu=',nu
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
             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN 
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                   ifailloc=1
                   ifailloc = pcubature(2, FFekgte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekgte integration at p=',p,' nu=',nu
                   ENDIF

                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = pcubature(2, FFekgte_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekgte_nocoll integration at p=',p,' nu=',nu
                   ENDIF

                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctepgte = intout(1)
             ifonctepgte = intout(2)

             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN 
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                   ifailloc=1
                   ifailloc = pcubature(2, FFekgne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekgne integration at p=',p,' nu=',nu
                   ENDIF

                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = pcubature(2, FFekgne_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekgne_nocoll integration at p=',p,' nu=',nu
                   ENDIF

                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctepgne = intout(1)
             ifonctepgne = intout(2)

             ifailloc=1;
             !Only calculate nonadiabatic part for electrons if el_type == 1
             IF (el_type == 1) THEN
                IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                   ifailloc=1
                   ifailloc = pcubature(2, FFekce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekce integration at p=',p,' nu=',nu
                   ENDIF

                ELSE ! Collisionless simulation, revert to faster single integral
                   ifailloc=1
                   ifailloc = pcubature(2, FFekce_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFekce_nocoll integration at p=',p,' nu=',nu
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

       
       ifailloc=1;
       !Only calculate nonadiabatic part for electrons if el_type == 1
       IF (el_type == 1) THEN 
          IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
             ifailloc=1
             ifailloc = pcubature(2, FFeke_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFeke integration at p=',p,' nu=',nu
             ENDIF

          ELSE ! Collisionless simulation, revert to faster single integral
             ifailloc=1
             ifailloc = pcubature(2, FFeke_nocoll_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)

             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL FFeke_nocoll integration at p=',p,' nu=',nu
             ENDIF

          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctepe = intout(1)
       ifonctepe = intout(2)

    END IF

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
  INTEGER :: npts, neval

  REAL(KIND=DBL), DIMENSION(ndim) :: a,b
  REAL(KIND=DBL) :: acc, cc, dd, relerr 
  REAL(KIND=DBL), DIMENSION(1) :: cc_cub, dd_cub
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
  
  REAL(KIND=DBL), DIMENSION(1) :: intout_cub
  REAL(KIND=DBL), DIMENSION(2) :: acc_cub
    

  omFFk = omega
  pFFk = p

  !! Integration bounds
  !! a,b are for the 1D ions.
  a(1) = 0.0d0 + barelyavoid
  b(1) = 1.0d0 - barelyavoid
  a(2) = 0.0d0
  b(2) = vuplim

  cc = a(1)
  dd = b(1)
  
  cc_cub(1) = cc
  dd_cub(1) = dd

  IF(QL_method.EQ.1) THEN
  
    DO ion=1,nions

        !! ION PARTICLE FLUX 
        ifailloc=1   
        ifailloc = hcubature(1, iFFkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        ifonctpi(ion) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkirot integration at p=',p,', nu=',nu,', ion=',ion
        ENDIF
        rfonctpi(ion) = 0.

        ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
        IF (phys_meth .NE. 0.0) THEN

           ifailloc=1
           ifailloc = hcubature(1, iFFkgtirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgti(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgtirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgti(ion)=0.

           ifailloc=1
           ifailloc = hcubature(1, iFFkgnirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgni(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgnirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgni(ion)=0.

           ifailloc=1
           ifailloc = hcubature(1, iFFkguirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgui(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkguirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgui(ion)=0.

           ifailloc=1
           ifailloc = hcubature(1, iFFkcirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpci(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkcirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpci(ion)=0.
!!!
           IF (phys_meth == 2) THEN

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFFekgtirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgti(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgtirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgti(ion) = 0.
              ENDIF
              rfonctepgti(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFFekgnirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgni(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgnirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgni(ion) = 0.
              ENDIF
              rfonctepgni(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFFekguirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgui(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekguirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgui(ion) = 0.
              ENDIF
              rfonctepgui(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFFekcirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepci(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekcirot integration at p=',p,', nu=',nu,', ion=',ion
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

        ifailloc=1
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = hcubature(1, iFFekirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctepi(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
        ELSE
           ifonctepi(ion) = 0.
        ENDIF
        rfonctepi(ion)=0.

        ! ION ang mom FLUX

        ifailloc=1
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = hcubature(1, iFFvkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctvpi(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
        ELSE
           ifonctvpi(ion) = 0.
        ENDIF
        rfonctvpi(ion)=0.

     ENDDO


     ! 2D INTEGRALS FOR ELECTRONS

     ifailloc=1;
     !Only calculate nonadiabatic part for electrons if el_type == 1
     IF (el_type == 1) THEN 

        IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
           ifailloc=1
           ifailloc = hcubature(2, FFkerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE ! Collisionless simulation, revert to faster single integral
           ifailloc=1
           ifailloc = hcubature(2, FFke_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFke_nocollrot integration at p=',p,' nu=',nu
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

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN 
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
              ifailloc=1
              ifailloc = hcubature(2, FFkgterot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgterot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = hcubature(2, FFkgte_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgte_nocollrot integration at p=',p,' nu=',nu
              ENDIF

           ENDIF

        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctpgte = intout(1)
        ifonctpgte = intout(2)

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN 
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
              ifailloc=1
              ifailloc = hcubature(2, FFkgnerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgnerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = hcubature(2, FFkgne_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkgne_nocollrot integration at p=',p,' nu=',nu
              ENDIF
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctpgne = intout(1)
        ifonctpgne = intout(2)

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
              ifailloc=1
              ifailloc = hcubature(2, FFkcerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkcerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = hcubature(2, FFkce_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFkce_nocollrot integration at p=',p,' nu=',nu
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
           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN 
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                 ifailloc=1
                 ifailloc = hcubature(2, FFekgterot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgterot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = hcubature(2, FFekgte_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgte_nocollrot integration at p=',p,' nu=',nu
                 ENDIF

              ENDIF

           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctepgte = intout(1)
           ifonctepgte = intout(2)

           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN 
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                 ifailloc=1
                 ifailloc = hcubature(2, FFekgnerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgnerot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = hcubature(2, FFekgne_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekgne_nocollrot integration at p=',p,' nu=',nu
                 ENDIF
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctepgne = intout(1)
           ifonctepgne = intout(2)

           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                 ifailloc=1
                 ifailloc = hcubature(2, FFekcerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekcerot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = hcubature(2, FFekce_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekce_nocollrot integration at p=',p,' nu=',nu
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

     
     ifailloc=1;   
     !Only calculate nonadiabatic part for electrons if el_type == 1
     IF (el_type == 1) THEN 
        IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
           ifailloc=1
           ifailloc = hcubature(2, FFekerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE ! Collisionless simulation, revert to faster single integral
           ifailloc=1
           ifailloc = hcubature(2, FFeke_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFeke_nocollrot integration at p=',p,' nu=',nu
           ENDIF
        ENDIF

     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctepe = intout(1)
     ifonctepe = intout(2)
  
  ELSE IF(QL_method.EQ.2) THEN
  
    DO ion=1,nions

        !! ION PARTICLE FLUX 
        ifailloc=1   
        ifailloc = pcubature(1, iFFkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        ifonctpi(ion) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkirot integration at p=',p,', nu=',nu,', ion=',ion
        ENDIF
        rfonctpi(ion) = 0.

        ! ADDITIONAL ION INTEGRALS TO SEPARATE TRANSPORT COMPONENTS
        IF (phys_meth .NE. 0.0) THEN

           ifailloc=1
           ifailloc = pcubature(1, iFFkgtirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgti(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgtirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgti(ion)=0.

           ifailloc=1
           ifailloc = pcubature(1, iFFkgnirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgni(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgnirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgni(ion)=0.

           ifailloc=1
           ifailloc = pcubature(1, iFFkguirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpgui(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkguirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpgui(ion)=0.

           ifailloc=1
           ifailloc = pcubature(1, iFFkcirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctpci(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkcirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
           rfonctpci(ion)=0.
!!!
           IF (phys_meth == 2) THEN

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFFekgtirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgti(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgtirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgti(ion) = 0.
              ENDIF
              rfonctepgti(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFFekgnirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgni(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgnirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgni(ion) = 0.
              ENDIF
              rfonctepgni(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFFekguirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepgui(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekguirot integration at p=',p,', nu=',nu,', ion=',ion
                 ENDIF
              ELSE
                 ifonctepgui(ion) = 0.
              ENDIF
              rfonctepgui(ion)=0.

              ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFFekcirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 ifonctepci(ion) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekcirot integration at p=',p,', nu=',nu,', ion=',ion
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

        ifailloc=1
        !This integral is extremely hard for pcubature, so hcubature is used
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = hcubature(1, iFFekirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctepi(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFekirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
        ELSE
           ifonctepi(ion) = 0.
        ENDIF
        rfonctepi(ion)=0.

        ! ION ang mom FLUX

        ifailloc=1
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = pcubature(1, iFFvkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           ifonctvpi(ion) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
           ENDIF
        ELSE
           ifonctvpi(ion) = 0.
        ENDIF
        rfonctvpi(ion)=0.

     ENDDO


     ! 2D INTEGRALS FOR ELECTRONS

     ifailloc=1;
     !Only calculate nonadiabatic part for electrons if el_type == 1
     IF (el_type == 1) THEN 

        IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
           ifailloc=1
           ifailloc = pcubature(2, FFkerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE ! Collisionless simulation, revert to faster single integral
           ifailloc=1
           ifailloc = pcubature(2, FFke_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFke_nocollrot integration at p=',p,' nu=',nu
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

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN 
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
              ifailloc=1
              ifailloc = pcubature(2, FFkgterot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgterot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = pcubature(2, FFkgte_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgte_nocollrot integration at p=',p,' nu=',nu
              ENDIF

           ENDIF

        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctpgte = intout(1)
        ifonctpgte = intout(2)

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN 
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
              ifailloc=1
              ifailloc = pcubature(2, FFkgnerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgnerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = pcubature(2, FFkgne_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkgne_nocollrot integration at p=',p,' nu=',nu
              ENDIF
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctpgne = intout(1)
        ifonctpgne = intout(2)

        ifailloc=1;
        !Only calculate nonadiabatic part for electrons if el_type == 1
        IF (el_type == 1) THEN
           IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
              ifailloc=1
              ifailloc = pcubature(2, FFkcerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkcerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE ! Collisionless simulation, revert to faster single integral
              ifailloc=1
              ifailloc = pcubature(2, FFkce_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
              
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFkce_nocollrot integration at p=',p,' nu=',nu
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
           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN 
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                 ifailloc=1
                 ifailloc = pcubature(2, FFekgterot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgterot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = pcubature(2, FFekgte_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgte_nocollrot integration at p=',p,' nu=',nu
                 ENDIF

              ENDIF

           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctepgte = intout(1)
           ifonctepgte = intout(2)

           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN 
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
                 ifailloc=1
                 ifailloc = pcubature(2, FFekgnerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgnerot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = pcubature(2, FFekgne_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekgne_nocollrot integration at p=',p,' nu=',nu
                 ENDIF
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctepgne = intout(1)
           ifonctepgne = intout(2)

           ifailloc=1;
           !Only calculate nonadiabatic part for electrons if el_type == 1
           IF (el_type == 1) THEN
              IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral 
                 ifailloc=1
                 ifailloc = pcubature(2, FFekcerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekcerot integration at p=',p,' nu=',nu
                 ENDIF
              ELSE ! Collisionless simulation, revert to faster single integral
                 ifailloc=1
                 ifailloc = pcubature(2, FFekce_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
                 
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekce_nocollrot integration at p=',p,' nu=',nu
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

                 
     ifailloc=1;
     !Only calculate nonadiabatic part for electrons if el_type == 1
     IF (el_type == 1) THEN 
        IF ( ABS(coll_flag) > epsD) THEN ! Collisional simulation, do double integral
           ifailloc=1
           ifailloc = pcubature(2, FFekerot_cubature, 2, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFekerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE ! Collisionless simulation, revert to faster single integral
           ifailloc=1
           ifailloc = pcubature(2, FFeke_nocollrot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, norm, intout, acc_cub)
           
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFeke_nocollrot integration at p=',p,' nu=',nu
           ENDIF
        ENDIF

     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctepe = intout(1)
     ifonctepe = intout(2)

  END IF
        
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
    INTEGER :: npts, neval

    REAL(KIND=DBL), DIMENSION(ndim) :: a,b
    REAL(KIND=DBL) :: acc, cc, dd, relerr 
    REAL(KIND=DBL), DIMENSION(1) :: cc_cub, dd_cub
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL), DIMENSION(nions) :: rfonctvpi
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctvpi

    REAL(KIND=DBL), DIMENSION(nf) :: intout

    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(1) :: intout_cub
    REAL(KIND=DBL), DIMENSION(2) :: acc_cub
        
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
    
    cc_cub(1) = cc
    dd_cub(1) = dd

    IF(QL_method.EQ.1) THEN
      DO ion=1,nions

            ! ION ang mom FLUX

            ifailloc=1
            IF (ninorm(p,ion) > min_ninorm) THEN
               ifailloc = hcubature(1, iFFvkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
               ifonctvpi(ion) = intout_cub(1)
               IF (ifailloc /= 0) THEN
                  IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
               ENDIF
            ELSE
               ifonctvpi(ion) = 0.
            ENDIF
            rfonctvpi(ion)=0.

     ENDDO
    
    ELSE IF(QL_method.EQ.2) THEN
      DO ion=1,nions

            ! ION ang mom FLUX

            ifailloc=1
            IF (ninorm(p,ion) > min_ninorm) THEN
               ifailloc = pcubature(1, iFFvkirot_cubature, 1, cc_cub, dd_cub, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
               ifonctvpi(ion) = intout_cub(1)
               IF (ifailloc /= 0) THEN
                  IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFFvkirot integration at p=',p,', nu=',nu,', ion=',ion
               ENDIF
            ELSE
               ifonctvpi(ion) = 0.
            ENDIF
            rfonctvpi(ion)=0.
     ENDDO
        
    END IF

    !The complex forms are reconstructed

    fonctvpi(:) = rfonctvpi(:) + ci * ifonctvpi(:)


  END SUBROUTINE momtrapQLintsrot

END MODULE calltrapQLints
