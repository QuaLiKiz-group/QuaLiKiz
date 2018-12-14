MODULE callpassQLints

  USE kind
  USE datcal
  USE datmat
  USE callpassints
  USE HCUB
  USE PCUB

  IMPLICIT NONE

CONTAINS

  SUBROUTINE passQLints( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, &
       & fonctcce, fonctcci, fonctece, foncteci,fonctecgte, fonctecgti, fonctecgne, fonctecgni,fonctecce, fonctecci)
    !-----------------------------------------------------------
    ! Calculates the passing particle integrals at the omega of the linear solutions
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! Includes resonance broadening
    ! To save (some) compuational time, only the (transport relevant) imaginary components are kept
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctce, fonctcgte, fonctcgne, fonctcce, fonctece
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctecgte, fonctecgne, fonctecce
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctci, fonctcgti, fonctcgni, fonctcci, foncteci
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctecgti, fonctecgni, fonctecci

    REAL(KIND=DBL), DIMENSION(ndim)    :: a, b, c
    REAL(KIND=DBL)    :: acc
    INTEGER           :: minpts, neval
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL), DIMENSION(nf) :: intout
    !INTEGER, PARAMETER :: numrgn = 4 !number of triangles in CUBATR
    INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration
    REAL, DIMENSION(ndim,0:ndim,4) :: vertices4 !for CUBATR 2D integration
    REAL(kind=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

    REAL(KIND=DBL)    :: rfonctce, rfonctcgte, rfonctcgne, rfonctcce, rfonctece  
    REAL(KIND=DBL)    :: rfonctecgte, rfonctecgne, rfonctecce 
    REAL(KIND=DBL)    :: ifonctce, ifonctcgte, ifonctcgne, ifonctcce, ifonctece
    REAL(KIND=DBL)    :: ifonctecgte, ifonctecgne, ifonctecce
    REAL(KIND=DBL)    :: intmult=4.D0
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctci, rfonctcgti, rfonctcgni, rfonctcci, rfoncteci
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctecgti, rfonctecgni, rfonctecci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctci, ifonctcgti, ifonctcgni, ifonctcci, ifoncteci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctecgti, ifonctecgni, ifonctecci
    
    REAL(KIND=DBL), DIMENSION(1) :: intout_cub
    REAL(KIND=DBL), DIMENSION(2) :: acc_cub

    INTEGER :: ifailloc
    
    
    omFkr = omega
    pFkr = p

    !    a(:) = -rkuplim

    a(:) = 0
    b(:) =  rkuplim

    !Integration boundaries
    vertices4(1:ndim,0,1) = (/0. , 0. /)
    vertices4(1:ndim,1,1) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,1) = (/rkuplim , -rkuplim /)

    vertices4(1:ndim,0,2) = (/0. , 0. /)
    vertices4(1:ndim,1,2) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,2) = (/-rkuplim , rkuplim /)

    vertices4(1:ndim,0,3) = (/0. , 0. /)
    vertices4(1:ndim,1,3) = (/-rkuplim , rkuplim /)
    vertices4(1:ndim,2,3) = (/-rkuplim , -rkuplim /)

    vertices4(1:ndim,0,4) = (/0. , 0. /)
    vertices4(1:ndim,1,4) = (/-rkuplim , -rkuplim /)
    vertices4(1:ndim,2,4) = (/rkuplim , -rkuplim /)

    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , rkuplim /)
    vertices1(1:ndim,2,1) = (/rkuplim , 0. /)

    !ELECTRON PARTICLE FLUX INTEGRALS
      
    IF(QL_method.EQ.1) THEN
      minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          intout(1)=0

          minpts=0; ifailloc=1
          ifailloc = hcubature(1, iFkstarrstare_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstare integration at p=',p,' nu=',nu
          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctce = intmult*intout(1); ifonctce=intmult*intout(2)

       !ADDITIONAL ELECTRON PARTICLE FLUX INTEGRALS
       IF (phys_meth .NE. 0.0) THEN

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstargte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargte integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          
             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstargne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargne integration at p=',p,' nu=',nu
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          
             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstarce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarce integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcce = intmult*intout(1); ifonctcce=intmult*intout(2)
!!!
          IF (phys_meth == 2) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = hcubature(1, iFekstarrstargte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargte integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = hcubature(1, iFekstarrstargne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargne integration at p=',p,' nu=',nu
                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = hcubature(1, iFekstarrstarce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarce integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecce = intmult*intout(1); ifonctecce=intmult*intout(2)
          ELSE
             rfonctecgte = 0.0
             ifonctecgte = 0.0
             rfonctecgne = 0.0
             ifonctecgne = 0.0
             rfonctecce = 0.0
             ifonctecce = 0.0
          ENDIF
       ELSE
          rfonctecgte = 0.0
          ifonctecgte = 0.0
          rfonctecgne = 0.0
          ifonctecgne = 0.0
          rfonctecce = 0.0
          ifonctecce = 0.0

          rfonctcgte = 0.0
          ifonctcgte = 0.0
          rfonctcgne = 0.0
          ifonctcgne = 0.0
          rfonctcce = 0.0
          ifonctcce = 0.0
       ENDIF

       !ELECTRON ENERGY INTEGRALS
       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

          intout(1)=0
          minpts=0; ifailloc=1
          ifailloc = hcubature(1, iFekstarrstare_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstare integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

          intout(1)=0
          minpts=0; ifailloc=1
          ifailloc = hcubature(1, iFkstarrstari_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstari integration at p=',p,' nu=',nu,' ion=',ion
          ENDIF

          rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstargti_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstargni_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = hcubature(1, iFkstarrstarci_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFekstarrstargti_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFekstarrstargni_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = hcubature(1, iFekstarrstarci_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecci(ion) = intmult*intout(1); ifonctecci(ion)=intmult*intout(2)
             ELSE
                rfonctecgti(ion) = 0.
                ifonctecgti(ion) = 0.
                rfonctecgni(ion) = 0.
                ifonctecgni(ion) = 0.
                rfonctecci(ion) = 0.
                ifonctecci(ion) = 0.
             ENDIF
          ELSE
             rfonctecgti(ion) = 0.
             ifonctecgti(ion) = 0.
             rfonctecgni(ion) = 0.
             ifonctecgni(ion) = 0.
             rfonctecci(ion) = 0.
             ifonctecci(ion) = 0.

             rfonctcgti(ion) = 0.
             ifonctcgti(ion) = 0.
             rfonctcgni(ion) = 0.
             ifonctcgni(ion) = 0.
             rfonctcci(ion) = 0.
             ifonctcci(ion) = 0.
          ENDIF

          !ION ENERGY INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

          intout(:)=0
          minpts=0; ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             ifailloc = hcubature(1, iFekstarrstari_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstari integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

       ENDDO

    ELSE IF(QL_method.EQ.2) THEN
      minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          intout(1)=0

          minpts=0; ifailloc=1
          ifailloc = pcubature(1, iFkstarrstare_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstare integration at p=',p,' nu=',nu
          ENDIF

       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctce = intmult*intout(1); ifonctce=intmult*intout(2)

       !ADDITIONAL ELECTRON PARTICLE FLUX INTEGRALS
       IF (phys_meth .NE. 0.0) THEN

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstargte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargte integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          
             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstargne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargne integration at p=',p,' nu=',nu
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
          
             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstarce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarce integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcce = intmult*intout(1); ifonctcce=intmult*intout(2)
!!!
          IF (phys_meth == 2) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = pcubature(1, iFekstarrstargte_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargte integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = pcubature(1, iFekstarrstargne_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargne integration at p=',p,' nu=',nu
                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

                intout(1)=0
                minpts=0; ifailloc=1
                ifailloc = pcubature(1, iFekstarrstarce_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                intout(2) = intout_cub(1)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarce integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecce = intmult*intout(1); ifonctecce=intmult*intout(2)
          ELSE
             rfonctecgte = 0.0
             ifonctecgte = 0.0
             rfonctecgne = 0.0
             ifonctecgne = 0.0
             rfonctecce = 0.0
             ifonctecce = 0.0
          ENDIF
       ELSE
          rfonctecgte = 0.0
          ifonctecgte = 0.0
          rfonctecgne = 0.0
          ifonctecgne = 0.0
          rfonctecce = 0.0
          ifonctecce = 0.0

          rfonctcgte = 0.0
          ifonctcgte = 0.0
          rfonctcgne = 0.0
          ifonctcgne = 0.0
          rfonctcce = 0.0
          ifonctcce = 0.0
       ENDIF

       !ELECTRON ENERGY INTEGRALS
       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

          intout(1)=0
          minpts=0; ifailloc=1
          ifailloc = pcubature(1, iFekstarrstare_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstare integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

          intout(1)=0
          minpts=0; ifailloc=1
          ifailloc = pcubature(1, iFkstarrstari_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
          intout(2) = intout_cub(1)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstari integration at p=',p,' nu=',nu,' ion=',ion
          ENDIF

          rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstargti_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstargni_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

             intout(1)=0
             minpts=0; ifailloc=1
             ifailloc = pcubature(1, iFkstarrstarci_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFekstarrstargti_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFekstarrstargni_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   ifailloc = pcubature(1, iFekstarrstarci_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                   intout(2) = intout_cub(1)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecci(ion) = intmult*intout(1); ifonctecci(ion)=intmult*intout(2)
             ELSE
                rfonctecgti(ion) = 0.
                ifonctecgti(ion) = 0.
                rfonctecgni(ion) = 0.
                ifonctecgni(ion) = 0.
                rfonctecci(ion) = 0.
                ifonctecci(ion) = 0.
             ENDIF
          ELSE
             rfonctecgti(ion) = 0.
             ifonctecgti(ion) = 0.
             rfonctecgni(ion) = 0.
             ifonctecgni(ion) = 0.
             rfonctecci(ion) = 0.
             ifonctecci(ion) = 0.

             rfonctcgti(ion) = 0.
             ifonctcgti(ion) = 0.
             rfonctcgni(ion) = 0.
             ifonctcgni(ion) = 0.
             rfonctcci(ion) = 0.
             ifonctcci(ion) = 0.
          ENDIF

          !ION ENERGY INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2

          intout(:)=0
          minpts=0; ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             ifailloc = pcubature(1, iFekstarrstari_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
             intout(2) = intout_cub(1)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstari integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

       ENDDO
    END IF
    
    !The complex forms are reconstructed
    fonctce = rfonctce + ci * ifonctce
    fonctcgte = rfonctcgte + ci * ifonctcgte
    fonctcgne = rfonctcgne + ci * ifonctcgne
    fonctcce = rfonctcce + ci * ifonctcce

    fonctece = rfonctece + ci * ifonctece
    fonctecgte = rfonctecgte + ci * ifonctecgte
    fonctecgne = rfonctecgne + ci * ifonctecgne
    fonctecce = rfonctecce + ci * ifonctecce

    fonctci = rfonctci + ci * ifonctci
    fonctcgti = rfonctcgti + ci * ifonctcgti
    fonctcgni = rfonctcgni + ci * ifonctcgni
    fonctcci =  rfonctcci + ci * ifonctcci
    foncteci = rfoncteci + ci * ifoncteci
    fonctecgti = rfonctecgti + ci * ifonctecgti
    fonctecgni = rfonctecgni + ci * ifonctecgni
    fonctecci =  rfonctecci + ci * ifonctecci

  END SUBROUTINE passQLints

  !*****************************************************************************************
  ! with rotation
  !*****************************************************************************************

  SUBROUTINE passQLintsrot( p, nu, omega, fonctce, fonctci, fonctcgte, fonctcgti, fonctcgne, fonctcgni, fonctcgui, &
       & fonctcce, fonctcci, fonctece, foncteci, fonctecgte, fonctecgti, fonctecgne, fonctecgni, fonctecgui, fonctecce, fonctecci, fonctvci)
    !-----------------------------------------------------------
    ! Calculates the passing particle integrals at the omega of the linear solutions
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! Includes resonance broadening
    ! To save (some) compuational time, only the (transport relevant) imaginary components are kept
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctce, fonctcgte, fonctcgne, fonctcce, fonctece
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctecgte, fonctecgne, fonctecce
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctci, fonctcgti, fonctcgni, fonctcgui, fonctcci, foncteci, fonctvci
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctecgti, fonctecgni, fonctecgui, fonctecci

    REAL(KIND=DBL), DIMENSION(ndim)    :: a, b, c
    REAL(KIND=DBL)    :: acc,ctot,c1,c2,c3,c4
    INTEGER           :: minpts, neval
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL), DIMENSION(nf) :: intout,xytest
    INTEGER, PARAMETER :: numrgn = 4 !number of triangles in CUBATR needed for -inf to +inf integration
    !INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration
    REAL, DIMENSION(ndim,0:ndim,4) :: vertices4 !for CUBATR 2D integration
    REAL(kind=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

    REAL(KIND=DBL)    :: rfonctce, rfonctcgte, rfonctcgne, rfonctcce, rfonctece 
    REAL(KIND=DBL)    :: rfonctecgte, rfonctecgne, rfonctecce 
    REAL(KIND=DBL)    :: ifonctce, ifonctcgte, ifonctcgne, ifonctcce, ifonctece
    REAL(KIND=DBL)    :: ifonctecgte, ifonctecgne, ifonctecce
    REAL(KIND=DBL)    :: intmult=1. ! with rotation need to integrate from -inf to +inf
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctci, rfonctcgti, rfonctcgni, rfonctcgui, rfonctcci, rfoncteci, rfonctvci
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctecgti, rfonctecgni, rfonctecgui, rfonctecci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctci, ifonctcgti, ifonctcgni, ifonctcgui, ifonctcci, ifoncteci, ifonctvci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctecgti, ifonctecgni, ifonctecgui, ifonctecci
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(1) :: intout_cub
    REAL(KIND=DBL), DIMENSION(2) :: acc_cub
       
    omFkr = omega
    pFkr = p

    a(:) = -rkuplim ! with rotation from -inf to +inf
    b(:) =  rkuplim

    !Integration boundaries
    vertices4(1:ndim,0,1) = (/0. , 0. /)
    vertices4(1:ndim,1,1) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,1) = (/rkuplim , -rkuplim /)

    vertices4(1:ndim,0,2) = (/0. , 0. /)
    vertices4(1:ndim,1,2) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,2) = (/-rkuplim , rkuplim /)

    vertices4(1:ndim,0,3) = (/0. , 0. /)
    vertices4(1:ndim,1,3) = (/-rkuplim , rkuplim /)
    vertices4(1:ndim,2,3) = (/-rkuplim , -rkuplim /)

    vertices4(1:ndim,0,4) = (/0. , 0. /)
    vertices4(1:ndim,1,4) = (/-rkuplim , -rkuplim /)
    vertices4(1:ndim,2,4) = (/rkuplim , -rkuplim /)

    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , rkuplim /)
    vertices1(1:ndim,2,1) = (/rkuplim , 0. /)

    ! with rotation use vertices 1 to vertices 4

    !ELECTRON PARTICLE FLUX INTEGRALS
  
    IF(QL_method.EQ.1) THEN

     minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
     IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = hcubature(1, iFkstarrstarerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarerot integration at p=',p,' nu=',nu
        ENDIF

     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctce = intmult*intout(1); ifonctce=intmult*intout(2)

     !ADDITIONAL ELECTRON PARTICLE FLUX INTEGRALS
     IF (phys_meth .NE. 0.0) THEN

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstargterot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargterot integration at p=',p,' nu=',nu
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstargnerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargnerot integration at p=',p,' nu=',nu
           ENDIF

        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstarcerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarcerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcce = intmult*intout(1); ifonctcce=intmult*intout(2)
!!!
        IF (phys_meth == 2) THEN
           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = hcubature(1, iFekstarrstargterot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargterot integration at p=',p,' nu=',nu
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = hcubature(1, iFekstarrstargnerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargnerot integration at p=',p,' nu=',nu
              ENDIF

           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = hcubature(1, iFekstarrstarcerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarcerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecce = intmult*intout(1); ifonctecce=intmult*intout(2)
        ELSE 
           rfonctecgte = 0.0
           ifonctecgte = 0.0
           rfonctecgne = 0.0
           ifonctecgne = 0.0
           rfonctecce = 0.0
           ifonctecce = 0.0
        ENDIF
     ELSE
        rfonctcgte = 0.0
        ifonctcgte = 0.0
        rfonctcgne = 0.0
        ifonctcgne = 0.0
        rfonctcce = 0.0
        ifonctcce = 0.0
        rfonctecgte = 0.0
        ifonctecgte = 0.0
        rfonctecgne = 0.0
        ifonctecgne = 0.0
        rfonctecce = 0.0
        ifonctecce = 0.0
     ENDIF

     !ELECTRON ENERGY INTEGRALS
     minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
     IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = hcubature(1, iFekstarrstarerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarerot integration at p=',p,' nu=',nu
        ENDIF
     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

     DO ion=1,nions

        !ION PARTICLE FLUX INTEGRALS
        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = hcubature(1, iFkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
        ENDIF

        rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

        !ADDITIONAL ION PARTICLE FLUX INTEGRALS
        IF (phys_meth .NE. 0.0) THEN
           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstargtirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstargnirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           IF ( (ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS) ) THEN
              intout(2) = 0
           ELSE
              minpts=0; ifailloc=1
              ifailloc = hcubature(1, iFkstarrstarguirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
              ENDIF
           ENDIF
           rfonctcgui(ion) = intmult*intout(1); ifonctcgui(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = hcubature(1, iFkstarrstarcirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFkstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)

           IF (phys_meth == 2) THEN
              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFekstarrstargtirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFekstarrstargnirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(1)=0
              IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion) < min_ninorm)) THEN
                 intout(2) = 0
              ELSE
                 minpts=0; ifailloc=1
                 ifailloc = hcubature(1, iFekstarrstarguirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgui(ion) = intmult*intout(1); ifonctecgui(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = hcubature(1, iFekstarrstarcirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecci(ion) = intmult*intout(1); ifonctecci(ion)=intmult*intout(2)
           ELSE
              rfonctecgti(ion) = 0.
              ifonctecgti(ion) = 0.
              rfonctecgni(ion) = 0.
              ifonctecgni(ion) = 0.
              rfonctecgui(ion) = 0.
              ifonctecgui(ion) = 0.
              rfonctecci(ion) = 0.
              ifonctecci(ion) = 0.
           ENDIF
        ELSE
           rfonctcgti(ion) = 0.
           ifonctcgti(ion) = 0.
           rfonctcgni(ion) = 0.
           ifonctcgni(ion) = 0.
           rfonctcgui(ion) = 0.
           ifonctcgui(ion) = 0.
           rfonctcci(ion) = 0.
           ifonctcci(ion) = 0.

           rfonctecgti(ion) = 0.
           ifonctecgti(ion) = 0.
           rfonctecgni(ion) = 0.
           ifonctecgni(ion) = 0.
           rfonctecgui(ion) = 0.
           ifonctecgui(ion) = 0.
           rfonctecci(ion) = 0.
           ifonctecci(ion) = 0.
        ENDIF

        !ION ENERGY INTEGRALS
        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(:)=0
        minpts=0; ifailloc=1
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = hcubature(1, iFekstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFekstarrstarirot integration at p=',p,' nu=',nu, ' ion =',ion
           ENDIF
        ENDIF
        rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

        !ION ang mom INTEGRALS
        minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0

        IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
           intout(2) = 0        
        ELSE
           minpts=0; ifailloc=1 
           ifailloc = hcubature(1, iFvkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
        ENDIF
        rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
     ENDDO
      
    ELSE IF(QL_method.EQ.2) THEN

     minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
     IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = pcubature(1, iFkstarrstarerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarerot integration at p=',p,' nu=',nu
        ENDIF

     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctce = intmult*intout(1); ifonctce=intmult*intout(2)

     !ADDITIONAL ELECTRON PARTICLE FLUX INTEGRALS
     IF (phys_meth .NE. 0.0) THEN

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstargterot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargterot integration at p=',p,' nu=',nu
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstargnerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargnerot integration at p=',p,' nu=',nu
           ENDIF

        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstarcerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarcerot integration at p=',p,' nu=',nu
           ENDIF
        ELSE
           intout(1)=0
           intout(2)=0
        ENDIF
        rfonctcce = intmult*intout(1); ifonctcce=intmult*intout(2)
!!!
        IF (phys_meth == 2) THEN
           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = pcubature(1, iFekstarrstargterot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargterot integration at p=',p,' nu=',nu
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = pcubature(1, iFekstarrstargnerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargnerot integration at p=',p,' nu=',nu
              ENDIF

           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
              intout(1)=0
              minpts=0; ifailloc=1
              ifailloc = pcubature(1, iFekstarrstarcerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarcerot integration at p=',p,' nu=',nu
              ENDIF
           ELSE
              intout(1)=0
              intout(2)=0
           ENDIF
           rfonctecce = intmult*intout(1); ifonctecce=intmult*intout(2)
        ELSE 
           rfonctecgte = 0.0
           ifonctecgte = 0.0
           rfonctecgne = 0.0
           ifonctecgne = 0.0
           rfonctecce = 0.0
           ifonctecce = 0.0
        ENDIF
     ELSE
        rfonctcgte = 0.0
        ifonctcgte = 0.0
        rfonctcgne = 0.0
        ifonctcgne = 0.0
        rfonctcce = 0.0
        ifonctcce = 0.0
        rfonctecgte = 0.0
        ifonctecgte = 0.0
        rfonctecgne = 0.0
        ifonctecgne = 0.0
        rfonctecce = 0.0
        ifonctecce = 0.0
     ENDIF

     !ELECTRON ENERGY INTEGRALS
     minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
     IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = pcubature(1, iFekstarrstarerot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarerot integration at p=',p,' nu=',nu
        ENDIF
     ELSE
        intout(1)=0
        intout(2)=0
     ENDIF
     rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

     DO ion=1,nions

        !ION PARTICLE FLUX INTEGRALS
        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0
        minpts=0; ifailloc=1
        ifailloc = pcubature(1, iFkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
        intout(2) = intout_cub(1)
        
        IF (ifailloc /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
        ENDIF

        rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

        !ADDITIONAL ION PARTICLE FLUX INTEGRALS
        IF (phys_meth .NE. 0.0) THEN
           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstargtirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstargnirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           IF ( (ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS) ) THEN
              intout(2) = 0
           ELSE
              minpts=0; ifailloc=1
              ifailloc = pcubature(1, iFkstarrstarguirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
              intout(2) = intout_cub(1)
              IF (ifailloc /= 0) THEN
                 IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
              ENDIF
           ENDIF
           rfonctcgui(ion) = intmult*intout(1); ifonctcgui(ion)=intmult*intout(2)

           minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
           intout(1)=0
           minpts=0; ifailloc=1
           ifailloc = pcubature(1, iFkstarrstarcirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFkstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
           rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
           IF (phys_meth == 2) THEN
              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFekstarrstargtirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFekstarrstargnirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(1)=0
              IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion) < min_ninorm)) THEN
                 intout(2) = 0
              ELSE
                 minpts=0; ifailloc=1
                 ifailloc = pcubature(1, iFekstarrstarguirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecgui(ion) = intmult*intout(1); ifonctecgui(ion)=intmult*intout(2)

              minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
              intout(:)=0
              minpts=0; ifailloc=1
              IF (ninorm(p,ion) > min_ninorm) THEN
                 ifailloc = pcubature(1, iFekstarrstarcirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
                 intout(2) = intout_cub(1)
                 IF (ifailloc /= 0) THEN
                    IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
                 ENDIF
              ENDIF
              rfonctecci(ion) = intmult*intout(1); ifonctecci(ion)=intmult*intout(2)
           ELSE
              rfonctecgti(ion) = 0.
              ifonctecgti(ion) = 0.
              rfonctecgni(ion) = 0.
              ifonctecgni(ion) = 0.
              rfonctecgui(ion) = 0.
              ifonctecgui(ion) = 0.
              rfonctecci(ion) = 0.
              ifonctecci(ion) = 0.
           ENDIF
        ELSE
           rfonctcgti(ion) = 0.
           ifonctcgti(ion) = 0.
           rfonctcgni(ion) = 0.
           ifonctcgni(ion) = 0.
           rfonctcgui(ion) = 0.
           ifonctcgui(ion) = 0.
           rfonctcci(ion) = 0.
           ifonctcci(ion) = 0.

           rfonctecgti(ion) = 0.
           ifonctecgti(ion) = 0.
           rfonctecgni(ion) = 0.
           ifonctecgni(ion) = 0.
           rfonctecgui(ion) = 0.
           ifonctecgui(ion) = 0.
           rfonctecci(ion) = 0.
           ifonctecci(ion) = 0.
        ENDIF

        !ION ENERGY INTEGRALS
        minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(:)=0
        minpts=0; ifailloc=1
        IF (ninorm(p,ion) > min_ninorm) THEN
           ifailloc = pcubature(1, iFekstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFekstarrstarirot integration at p=',p,' nu=',nu, ' ion =',ion
           ENDIF
        ENDIF
        rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

        !ION ang mom INTEGRALS
        minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0

        IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
           intout(2) = 0        
        ELSE
           minpts=0; ifailloc=1 
           ifailloc = pcubature(1, iFvkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
        ENDIF
        rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
     ENDDO
  
    
    END IF
    
    
    !The complex forms are reconstructed
    fonctce = rfonctce + ci * ifonctce
    fonctcgte = rfonctcgte + ci * ifonctcgte
    fonctcgne = rfonctcgne + ci * ifonctcgne
    fonctcce = rfonctcce + ci * ifonctcce

    fonctece = rfonctece + ci * ifonctece
    fonctecgte = rfonctecgte + ci * ifonctecgte
    fonctecgne = rfonctecgne + ci * ifonctecgne
    fonctecce = rfonctecce + ci * ifonctecce

    fonctci = rfonctci + ci * ifonctci
    fonctcgti = rfonctcgti + ci * ifonctcgti
    fonctcgni = rfonctcgni + ci * ifonctcgni
    fonctcgui = rfonctcgui + ci * ifonctcgui
    fonctcci =  rfonctcci + ci * ifonctcci

    foncteci = rfoncteci + ci * ifoncteci
    fonctecgti = rfonctecgti + ci * ifonctecgti
    fonctecgni = rfonctecgni + ci * ifonctecgni
    fonctecgui = rfonctecgui + ci * ifonctecgui
    fonctecci =  rfonctecci + ci * ifonctecci

    fonctvci = rfonctvci + ci * ifonctvci

  END SUBROUTINE passQLintsrot

  SUBROUTINE mompassQLintsrot( p, nu, omega, fonctvci)
    !-----------------------------------------------------------
    ! Calculates the passing particle integrals at the omega of the linear solutions. ANG MOM ONLY
    ! Integrals (with switch) in order to identify particle flux contributions due to An, At and curvature
    ! Includes resonance broadening
    ! To save (some) compuational time, only the (transport relevant) imaginary components are kept
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fonctvci

    REAL(KIND=DBL), DIMENSION(ndim)    :: a, b, c
    REAL(KIND=DBL)    :: acc,ctot,c1,c2,c3,c4
    INTEGER           :: minpts, neval
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr

    REAL(KIND=DBL), DIMENSION(nf) :: intout,xytest
    INTEGER, PARAMETER :: numrgn = 4 !number of triangles in CUBATR needed for -inf to +inf integration
    !INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration
    REAL, DIMENSION(ndim,0:ndim,4) :: vertices4 !for CUBATR 2D integration
    REAL(kind=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

    REAL(KIND=DBL)    :: intmult=1. ! with rotation need to integrate from -inf to +inf
    REAL(KIND=DBL), DIMENSION(nions) :: rfonctvci
    REAL(KIND=DBL), DIMENSION(nions) :: ifonctvci

    INTEGER :: ifailloc
    
    REAL(KIND=DBL), DIMENSION(1) :: intout_cub
    REAL(KIND=DBL), DIMENSION(2) :: acc_cub

    

    omFkr = omega
    pFkr = p

    a(:) = -rkuplim ! with rotation from -inf to +inf
    b(:) =  rkuplim

    !Integration boundaries
    vertices4(1:ndim,0,1) = (/0. , 0. /)
    vertices4(1:ndim,1,1) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,1) = (/rkuplim , -rkuplim /)

    vertices4(1:ndim,0,2) = (/0. , 0. /)
    vertices4(1:ndim,1,2) = (/rkuplim , rkuplim /)
    vertices4(1:ndim,2,2) = (/-rkuplim , rkuplim /)

    vertices4(1:ndim,0,3) = (/0. , 0. /)
    vertices4(1:ndim,1,3) = (/-rkuplim , rkuplim /)
    vertices4(1:ndim,2,3) = (/-rkuplim , -rkuplim /)

    vertices4(1:ndim,0,4) = (/0. , 0. /)
    vertices4(1:ndim,1,4) = (/-rkuplim , -rkuplim /)
    vertices4(1:ndim,2,4) = (/rkuplim , -rkuplim /)

    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , rkuplim /)
    vertices1(1:ndim,2,1) = (/rkuplim , 0. /)

    ! with rotation use vertices 1 to vertices 4

    !ELECTRON PARTICLE FLUX INTEGRALS
    
    IF(QL_method.EQ.1) THEN
      DO ion=1,nions

        !ION ang mom INTEGRALS
        minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0

        IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
           intout(2) = 0        
        ELSE
           minpts=0; ifailloc=1  
           ifailloc = hcubature(1, iFvkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of hcubature QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
        ENDIF
        rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
     ENDDO
     ELSE IF(QL_method.EQ.2) THEN
      DO ion=1,nions

        !ION ang mom INTEGRALS
        minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
        intout(1)=0

        IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
           intout(2) = 0        
        ELSE
           minpts=0; ifailloc=1  
           ifailloc = pcubature(1, iFvkstarrstarirot_cubature, ndim, a, b, maxpts, reqabsacc_QL, reqrelacc_QL, 1, intout_cub, acc_cub)
           intout(2) = intout_cub(1)
           IF (ifailloc /= 0) THEN
              IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of pcubature QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
           ENDIF
        ENDIF
        rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
     ENDDO
    END IF

    !The complex forms are reconstructed
    fonctvci = rfonctvci + ci * ifonctvci

  END SUBROUTINE mompassQLintsrot


END MODULE callpassQLints
