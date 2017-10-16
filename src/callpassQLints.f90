MODULE callpassQLints

  USE kind
  USE datcal
  USE datmat
  USE callpassints
  USE CUI

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

    IF (inttype == 1) THEN

       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

!!$          CALL CUBATR(ndim,nf,Fkstarrstare_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)

          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstare integration at p=',p,' nu=',nu
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
!!$             CALL CUBATR(ndim,nf,Fkstarrstargte_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargte integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte= intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             CALL CUBATR(ndim,nf,Fkstarrstargne_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargne integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             CALL CUBATR(ndim,nf,Fkstarrstarce_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarce integration at p=',p,' nu=',nu
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
!!$                CALL CUBATR(ndim,nf,Fekstarrstargte_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargte integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                CALL CUBATR(ndim,nf,Fekstarrstargne_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargne integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                CALL CUBATR(ndim,nf,Fekstarrstarce_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarce integration at p=',p,' nu=',nu
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
!!$          CALL CUBATR(ndim,nf,Fekstarrstare_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstare integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          CALL CUBATR(ndim,nf,Fkstarrstari_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstari integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfonctci(ion) = intmult*intout(1); ifonctci(ion)= intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstargti_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargti integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstargni_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargni integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstarci_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarci integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstargti_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargti integration at p=',p,' nu=',nu,', ion=',ion
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstargni_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargni integration at p=',p,' nu=',nu,', ion=',ion
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstarci_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarci integration at p=',p,' nu=',nu,', ion=',ion
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
!!$          CALL CUBATR(ndim,nf,Fekstarrstari_cub,numrgn,vertices1,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstari integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

       ENDDO

    ELSEIF (inttype == 2) THEN 


       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstare,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstare integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0

          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstare,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstare integration at p=',p,' nu=',nu
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
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargte,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargte integration at p=',p,' nu=',nu
!!$             ENDIF
!!$             ifailloc=1
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargte,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargte integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargne,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargne integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargne,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargne integration at p=',p,' nu=',nu
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarce,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarce integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarce,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarce integration at p=',p,' nu=',nu
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
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargte,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargte integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargte,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargte integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargne,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargne integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargne,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargne integration at p=',p,' nu=',nu
                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarce,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarce integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarce,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarce integration at p=',p,' nu=',nu
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
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstare,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstare integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstare,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstare integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstari,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstari integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstari,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstari integration at p=',p,' nu=',nu,' ion=',ion
          ENDIF

          rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargti,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargti integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargti,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargni,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargni integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargni,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarci,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarci integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarci,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargti,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargti integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargti,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargti integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargni,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargni integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargni,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargni integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarci,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarci integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarci,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarci integration at p=',p,' nu=',nu,' ion=',ion
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
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstari,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstari integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(:)=0
          minpts=0; ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstari,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstari integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

       ENDDO
    ENDIF

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

    IF (inttype == 1) THEN

       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN

!!$          CALL CUBATR(ndim,nf,Fkstarrstarerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)

          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarerot integration at p=',p,' nu=',nu
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
!!$             CALL CUBATR(ndim,nf,Fkstarrstargterot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargterot integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte= intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             CALL CUBATR(ndim,nf,Fkstarrstargnerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargnerot integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             CALL CUBATR(ndim,nf,Fkstarrstarcerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarcerot integration at p=',p,' nu=',nu
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
!!$                CALL CUBATR(ndim,nf,Fekstarrstargterot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargterot integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                CALL CUBATR(ndim,nf,Fekstarrstargnerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargnerot integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                CALL CUBATR(ndim,nf,Fekstarrstarcerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarcerot integration at p=',p,' nu=',nu
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
!!$          CALL CUBATR(ndim,nf,Fekstarrstarerot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarerot integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          CALL CUBATR(ndim,nf,Fkstarrstarirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarirot integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfonctci(ion) = intmult*intout(1); ifonctci(ion)= intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstargtirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargtirot integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstargnirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargnirot integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstarguirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstargnirot integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcgui(ion) = intmult*intout(1); ifonctcgui(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             CALL CUBATR(ndim,nf,Fkstarrstarcirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fkstarrstarcirot integration at p=',p,' nu=',nu,', ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstargtirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargtirot integration at p=',p,' nu=',nu,', ion=',ion
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstargnirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargnirot integration at p=',p,' nu=',nu,', ion=',ion
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstarguirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstargnirot integration at p=',p,' nu=',nu,', ion=',ion
                ENDIF
                rfonctecgui(ion) = intmult*intout(1); ifonctecgui(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                CALL CUBATR(ndim,nf,Fekstarrstarcirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$                     ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarcirot integration at p=',p,' nu=',nu,', ion=',ion
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
!!$          CALL CUBATR(ndim,nf,Fekstarrstarirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarirot integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

	  !ION ang mom INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          CALL CUBATR(ndim,nf,Fvkstarrstarirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarirot integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)

       ENDDO

    ELSEIF (inttype == 2) THEN

       minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
       IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarerot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarerot integration at p=',p,' nu=',nu
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
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargterot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargterot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargterot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargterot integration at p=',p,' nu=',nu
             ENDIF
          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgte = intmult*intout(1); ifonctcgte=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargnerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargnerot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargnerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargnerot integration at p=',p,' nu=',nu
             ENDIF

          ELSE
             intout(1)=0
             intout(2)=0
          ENDIF
          rfonctcgne = intmult*intout(1); ifonctcgne=intmult*intout(2)

          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
          IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarcerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarcerot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarcerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarcerot integration at p=',p,' nu=',nu
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
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargterot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargterot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargterot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargterot integration at p=',p,' nu=',nu
                ENDIF
             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgte = intmult*intout(1); ifonctecgte=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargnerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargnerot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargnerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargnerot integration at p=',p,' nu=',nu
                ENDIF

             ELSE
                intout(1)=0
                intout(2)=0
             ENDIF
             rfonctecgne = intmult*intout(1); ifonctecgne=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
             IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nu) .EQV. .FALSE.) ) )  THEN
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarcerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarcerot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarcerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarcerot integration at p=',p,' nu=',nu
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
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarerot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarerot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarerot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarerot integration at p=',p,' nu=',nu
          ENDIF
       ELSE
          intout(1)=0
          intout(2)=0
       ENDIF
       rfonctece = intmult*intout(1); ifonctece=intmult*intout(2)

       DO ion=1,nions

          !ION PARTICLE FLUX INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarirot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0
          minpts=0; ifailloc=1
          CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
          ENDIF

          rfonctci(ion) = intmult*intout(1); ifonctci(ion)=intmult*intout(2)

          !ADDITIONAL ION PARTICLE FLUX INTEGRALS
          IF (phys_meth .NE. 0.0) THEN
             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargtirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargtirot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargtirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgti(ion) = intmult*intout(1); ifonctcgti(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstargnirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstargnirot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstargnirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcgni(ion) = intmult*intout(1); ifonctcgni(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarguirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarguirot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             IF ( (ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS) ) THEN
                intout(2) = 0
             ELSE
                minpts=0; ifailloc=1
                CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarguirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                IF (ifailloc /= 0) THEN
                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
                ENDIF
             ENDIF
             rfonctcgui(ion) = intmult*intout(1); ifonctcgui(ion)=intmult*intout(2)

             minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$             ifailloc=1
!!$             CALL d01fcf(ndim,a,b,minpts,maxpts,rFkstarrstarcirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$             IF (ifailloc /= 0) THEN
!!$                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFkstarrstarcirot integration at p=',p,' nu=',nu
!!$             ENDIF
             intout(1)=0
             minpts=0; ifailloc=1
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFkstarrstarcirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFkstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
             rfonctcci(ion) = intmult*intout(1); ifonctcci(ion)=intmult*intout(2)
!!!
             IF (phys_meth == 2) THEN
                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargtirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargtirot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargtirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargtirot integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgti(ion) = intmult*intout(1); ifonctecgti(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstargnirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstargnirot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstargnirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstargnirot integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgni(ion) = intmult*intout(1); ifonctecgni(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarguirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarguirot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(1)=0
                IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion) < min_ninorm)) THEN
                   intout(2) = 0
                ELSE
                   minpts=0; ifailloc=1
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarguirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarguirot integration at p=',p,' nu=',nu,' ion=',ion
                   ENDIF
                ENDIF
                rfonctecgui(ion) = intmult*intout(1); ifonctecgui(ion)=intmult*intout(2)

                minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$                ifailloc=1
!!$                CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarcirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$                IF (ifailloc /= 0) THEN
!!$                   IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarcirot integration at p=',p,' nu=',nu
!!$                ENDIF
                intout(:)=0
                minpts=0; ifailloc=1
                IF (ninorm(p,ion) > min_ninorm) THEN
                   CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarcirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
                   IF (ifailloc /= 0) THEN
                      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarcirot integration at p=',p,' nu=',nu,' ion=',ion
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
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFekstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarirot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(:)=0
          minpts=0; ifailloc=1
          IF (ninorm(p,ion) > min_ninorm) THEN
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFekstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFekstarrstarirot integration at p=',p,' nu=',nu, ' ion =',ion
             ENDIF
          ENDIF
          rfoncteci(ion) = intmult*intout(1); ifoncteci(ion)=intmult*intout(2)

          !ION ang mom INTEGRALS
          minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFvkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarirot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0

          IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
             intout(2) = 0        
          ELSE
             minpts=0; ifailloc=1                                                 
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFvkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
          ENDIF
          rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
       ENDDO
    ENDIF

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

    IF (inttype == 1) THEN

       DO ion=1,nions

	  !ION ang mom INTEGRALS
          minpts = 0; ifailloc=1; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          CALL CUBATR(ndim,nf,Fvkstarrstarirot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$               ifailloc,neval,abaccQL2,relaccQL2,restar,minpts,maxpts,key,job,tune)
          IF (ifailloc /= 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I3,A,I3)") 'ifailloc = ',ifailloc,'. Abnormal termination of CUBATR QL Fekstarrstarirot integration at p=',p,' nu=',nu,', ion=',ion
          ENDIF
          rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)

       ENDDO

    ELSEIF (inttype == 2) THEN

       DO ion=1,nions

          !ION ang mom INTEGRALS
          minpts = 0; ifailloc=1 ; reerrarr(:)=1.d-1; rgtype(:)=2
!!$          ifailloc=1
!!$          CALL d01fcf(ndim,a,b,minpts,maxpts,rFvkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(1),ifailloc)
!!$          IF (ifailloc /= 0) THEN
!!$             IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL rFekstarrstarirot integration at p=',p,' nu=',nu
!!$          ENDIF
          intout(1)=0

          IF ( ((ABS(Machpar(p)) < epsS) .AND. (ABS(Aupar(p)) < epsS) .AND. (ABS(gammaE(p)) < epsS)) .OR. (ninorm(p,ion)<min_ninorm) ) THEN
             intout(2) = 0        
          ELSE
             minpts=0; ifailloc=1                                                 
             CALL d01fcf(ndim,a,b,minpts,maxpts,iFvkstarrstarirot,relaccQL2,acc,lenwrk,wrkstr,intout(2),ifailloc)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I3,A,I0,A,I0)") 'ifailloc = ',ifailloc,'. Abnormal termination of 2DNAG QL iFvkstarrstarirot integration at p=',p,' nu=',nu,' ion=',ion
             ENDIF
          ENDIF
          rfonctvci(ion) = intmult*intout(1); ifonctvci(ion)=intmult*intout(2)
       ENDDO
    ENDIF

    !The complex forms are reconstructed
    fonctvci = rfonctvci + ci * ifonctvci

  END SUBROUTINE mompassQLintsrot


END MODULE callpassQLints
