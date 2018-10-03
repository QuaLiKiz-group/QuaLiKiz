!Calculate integrals
MODULE mod_fonct
  USE kind
  USE datcal
  USE datmat
  USE dispfuncs
  USE calltrapints
  USE callpassints
  USE calltotalints
  USE CUI
  USE HCUB
  USE PCUB
  USE integration_routines

  IMPLICIT NONE

CONTAINS


  SUBROUTINE calcfonct_hcubature(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu
    
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Integration is over a square
    a(1) = 0._DBL
    a(2) = 0._DBL
    b(1) = 1._DBL
    b(2) = 1._DBL
    
    !Adiabatic term, also used as scale factor
    fdata(1) = Ac(p)
    
    !Set integration limits inside fdata
    
    !trapped
    fdata(2) = 0._DBL
    fdata(3) = 0._DBL
    fdata(4) = 1._DBL - barelyavoid
    fdata(5) = vuplim
    
    !passing
    fdata(6) = 0._DBL
    fdata(7) = 0._DBL
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = hcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonct_total = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_hcubature
  SUBROUTINE calcfonct_pcubature(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu
    
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Integration is over a square
    a(1) = 0._DBL
    a(2) = 0._DBL
    b(1) = 1._DBL
    b(2) = 1._DBL
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !Set integration limits inside fdata
    
    !trapped
    fdata(2) = 0._DBL
    fdata(3) = 0._DBL
    fdata(4) = 1._DBL - barelyavoid
    fdata(5) = vuplim
    
    !passing
    fdata(6) = 0._DBL
    fdata(7) = 0._DBL
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = pcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonct_total = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_pcubature
  
  SUBROUTINE calcfonct_hcubaturec(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(1) = 0.0d0 
    a(2) = 0.0d0 
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = hcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonctc = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_hcubaturec
  
  SUBROUTINE calcfonct_pcubaturec(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(1) = 0.0d0 
    a(2) = 0.0d0 
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = pcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonctc = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_pcubaturec
  
  SUBROUTINE calcfonct_hcubaturep(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
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
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    fdata(2) = vuplim
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = hcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonctp = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_hcubaturep
  
  SUBROUTINE calcfonct_pcubaturep(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim, norm
    
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
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    fdata(2) = vuplim
    
    !set the convergence criterion to be the L_2 norm
    norm = 2 
    fdim = ndim
    
    
    ifailloc = pcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
    
    IF (ifailloc /= 0) THEN
      IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
        &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
    END IF
    
    intout = intout * ABS(Ac(p))
    fonctp = intout(1) + ci * intout(2)  
  
  END SUBROUTINE calcfonct_pcubaturep
  

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
    INTEGER            :: minpts, neval, npts, iffailloc  !output number of integrand evaluations
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
    INTEGER :: ifailloc
    REAL(KIND=DBL)    :: rfonctpe, ifonctpe, rfonctpiz, ifonctpiz

    REAL(KIND=DBL), DIMENSION(nf) :: intout
    INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration

    REAL(kind=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

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
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    rfonctpiz = integrate_1d(cc, dd, relacc1, npts, relerr, rFFkiz, lw, verbose, p, nu, omega, '1DNAG rFFkiz')
    ifonctpiz = integrate_1d(cc, dd, relacc1, npts, relerr, iFFkiz, lw, verbose, p, nu, omega, '1DNAG iFFkiz')

    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , vuplim /)
    vertices1(1:ndim,2,1) = (/1. - barelyavoid , 0. /)
    rgtype(:)=2
    minpts = 0; ifailloc=-1; reerrarr(:)=1.d-1
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
          IF (inttype == 1) THEN
!!$             CALL CUBATR(ndim,nf,FFke_cub,1,vertices1,rgtype,intout,reerrarr,&
!!$                  ifailloc,neval,abacc,relacc2,restar,minpts,maxpts,key,job,tune)

             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                     &'. Abnormal termination of CUBATR FFke_cub integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
             ENDIF
          ELSEIF (inttype == 2) THEN
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
          ENDIF
       ELSE !Collisionless simulation, revert to faster single integral
           intout(1) = integrate_1d(cc, dd, relacc1, npts, relerr, rFFke_nocoll, lw, verbose, p, nu, omega, '1DNAG rFFke_nocoll')
           intout(2) = integrate_1d(cc, dd, relacc1, npts, relerr, iFFke_nocoll, lw, verbose, p, nu, omega, '1DNAG iFFke_nocoll')
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF

    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)

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
    !    REAL(KIND=DBL), DIMENSION(lims) :: rarray,iarray

    REAL(KIND=DBL)     :: rfonctc, ifonctc !outputs of integrals

    INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,2) :: vertices2 !for CUBATR 2D integration
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,4) :: vertices4 !for CUBATR 2D integration

    REAL(KIND=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

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

    !Other vertices are kept for testing purposes of various settings in the CUBATR routine
    !in the end only vertices1 is used for now
    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , rkuplim /)
    vertices1(1:ndim,2,1) = (/rkuplim , 0. /)

    vertices2(1:ndim,0,1) = (/0. , 0. /)
    vertices2(1:ndim,1,1) = (/rkuplim , 0. /)
    vertices2(1:ndim,2,1) = (/0. , rkuplim /)

    vertices2(1:ndim,0,2) = (/0. , rkuplim /)
    vertices2(1:ndim,1,2) = (/rkuplim , 0. /)
    vertices2(1:ndim,2,2) = (/rkuplim , rkuplim /)

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

    rgtype(:)=2
    minpts = 0; ifailloc=-1; reerrarr(:)=1.d-1

    IF (inttype == 1) THEN
!       CALL CUBATR(ndim,nf,Fkstarrstar_cub,1,vertices1,rgtype,intout,reerrarr,&
!            ifailloc,neval,abacc,relacc2,restar,minpts,maxpts,key,job,tune)

       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of Fkstarrstar integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
       ENDIF
    ELSEIF (inttype == 2) THEN
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
    INTEGER, PARAMETER :: numrgn = 1 !number of triangles in CUBATR
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration

    REAL(kind=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

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
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    rfonctpiz = integrate_1d(cc, dd, relacc1, npts, relerr, rFFkizrot, lw, verbose, p, nu, omega, '1DNAG rFFkizrot')
    ifonctpiz = integrate_1d(cc, dd, relacc1, npts, relerr, iFFkizrot, lw, verbose, p, nu, omega, '1DNAG rFFkizrot')

    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , vuplim /)
    vertices1(1:ndim,2,1) = (/1. - barelyavoid , 0. /)
    rgtype(:)=2
    minpts = 0; ifailloc=-1; reerrarr(:)=1.d-1
    !Only calculate nonadiabatic part for electrons if el_type == 1
    IF (el_type == 1) THEN 
       IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
          IF (inttype == 1) THEN
!             CALL CUBATR(ndim,nf,FFkerot_cub,1,vertices1,rgtype,intout,reerrarr,&
!                  ifailloc,neval,abacc,relacc2,restar,minpts,maxpts,key,job,tune)
             IF (ifailloc /= 0) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
                     &'. Abnormal termination of CUBATR FFkerot_cub integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
             ENDIF
          ELSEIF (inttype == 2) THEN
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
          ENDIF
       ELSE !collisionless integral, do single integral
         intout(1) = integrate_1d(cc, dd, relacc1, npts, relerr, rFFke_nocollrot, lw, verbose, p, nu, omega, '1DNAG rFFke_nocollrot')
         intout(2) = integrate_1d(cc, dd, relacc1, npts, relerr, iFFke_nocollrot, lw, verbose, p, nu, omega, '1DNAG iFFke_nocollrot')
       ENDIF
    ELSE
       intout(1)=0
       intout(2)=0
    ENDIF

    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)

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

    REAL(KIND=DBL), DIMENSION(ndim) :: a, b, c, dd , xy !, a1,a2,a3,a4,b1,b2,b3,b4
    REAL(KIND=DBL), DIMENSION(lenwrk) :: wrkstr
!!$    REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: wrkstr2
    REAL(KIND=DBL), DIMENSION(nf) :: intout!,intout1,intout2,intout3,intout4

    REAL(KIND=DBL)     :: acc
    INTEGER            :: maxpts2, npts, minpts, neval, iii,jjj,rc,lenwrk2 !, lims=100
    INTEGER :: ifailloc
    !    REAL(KIND=DBL), DIMENSION(lims) :: rarray,iarray

    REAL(KIND=DBL)     :: rfonctc, ifonctc !outputs of integrals

    INTEGER, PARAMETER :: numrgn = 4 !number of triangles in CUBATR needed with rotation -inf to +inf
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,1) :: vertices1 !for CUBATR 2D integration
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,2) :: vertices2 !for CUBATR 2D integration
    REAL(KIND=DBL), DIMENSION(ndim,0:ndim,4) :: vertices4 !for CUBATR 2D integration

    REAL(KIND=DBL) , DIMENSION(nf) :: reerrarr !NOT USED. Estimated absolute error after CUBATR restart
    INTEGER, DIMENSION(numrgn) :: rgtype

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

    !Other vertices are kept for testing purposes of various settings in the CUBATR routine
    !in the end only vertices4 is used for now
    vertices1(1:ndim,0,1) = (/0. , 0. /)
    vertices1(1:ndim,1,1) = (/0. , rkuplim /)
    vertices1(1:ndim,2,1) = (/rkuplim , 0. /)

    vertices2(1:ndim,0,1) = (/0. , 0. /)
    vertices2(1:ndim,1,1) = (/rkuplim , 0. /)
    vertices2(1:ndim,2,1) = (/0. , rkuplim /)

    vertices2(1:ndim,0,2) = (/0. , rkuplim /)
    vertices2(1:ndim,1,2) = (/rkuplim , 0. /)
    vertices2(1:ndim,2,2) = (/rkuplim , rkuplim /)

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

    rgtype(:)=2
    minpts = 0; ifailloc=-1; reerrarr(:)=1.d-1

    IF (inttype == 1) THEN
!!$       CALL CUBATR(ndim,nf,Fkstarrstarrot_cub,numrgn,vertices4,rgtype,intout,reerrarr,&
!!$            ifailloc,neval,abacc,relacc2,restar,minpts,maxpts,key,job,tune)

       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of Fkstarrstarrot integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
       ENDIF
    ELSEIF (inttype == 2) THEN

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

    ENDIF

    fonctc = intout(1) + ci * intout(2)


  END SUBROUTINE calcfonctrotc

END MODULE mod_fonct
