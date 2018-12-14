!Calculate integrals
MODULE mod_fonct
  USE kind
  USE datcal
  USE datmat
  USE dispfuncs
  USE calltrapints
  USE callpassints
  USE calltotalints
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
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
    
    
  
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
    INTEGER :: ifailloc, fdim
    
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
    fdata(2) = 0._DBL + barelyavoid
    fdata(3) = 0._DBL
    fdata(4) = 1._DBL - barelyavoid
    fdata(5) = vuplim
    
    !passing
    fdata(6) = 0._DBL
    fdata(7) = 0._DBL
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
  
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
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
      ifailloc = hcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = hcubature(fdim, ipassing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
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
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
      ifailloc = pcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = pcubature(fdim, ipassing_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonct_pcubaturec
  
  SUBROUTINE calcfonct_hcubaturec2(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
    
      ifailloc = hcubature(fdim, passing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = hcubature(fdim, passing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = hcubature(fdim, rpassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = hcubature(fdim, ipassing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = hcubature(fdim, ipassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonct_hcubaturec2
  
  SUBROUTINE calcfonct_pcubaturec2(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = pcubature(fdim, passing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = pcubature(fdim, passing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = pcubature(fdim, rpassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = pcubature(fdim, ipassing_ions_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = pcubature(fdim, ipassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonct_pcubaturec2
  
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
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral   
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
      
    ELSE !collisionless, single integrals
    
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
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
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu
    
    !Set integration limits. (1) for kappa and (2) for v (for the electrons)
    a(1) = 0.0d0 + barelyavoid
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    fdata(2) = vuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral   
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrapped_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
      
    ELSE !collisionless, single integrals
    
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonct_pcubaturep
  
  SUBROUTINE calcfonct_hcubature_newt(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
      
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF

    END IF
  
  END SUBROUTINE calcfonct_hcubature_newt
  SUBROUTINE calcfonct_pcubature_newt(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
      
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_nocoll_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
  
  END SUBROUTINE calcfonct_pcubature_newt
  
  SUBROUTINE calcfonct_hcubaturec_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
      ifailloc = hcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = hcubature(fdim, ipassing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonct_hcubaturec_newt
  
  SUBROUTINE calcfonct_pcubaturec_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
      ifailloc = pcubature(fdim, passing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = pcubature(fdim, ipassing_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonct_pcubaturec_newt
  
  SUBROUTINE calcfonct_hcubaturec2_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = hcubature(fdim, passing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = hcubature(fdim, passing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = hcubature(fdim, rpassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = hcubature(fdim, ipassing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = hcubature(fdim, ipassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonct_hcubaturec2_newt
  
  SUBROUTINE calcfonct_pcubaturec2_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = pcubature(fdim, passing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = pcubature(fdim, passing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = pcubature(fdim, rpassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = pcubature(fdim, ipassing_ions_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = pcubature(fdim, ipassing_electrons_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF 
  
  END SUBROUTINE calcfonct_pcubaturec2_newt
  
  SUBROUTINE calcfonct_hcubaturep_newt(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !collisions, do double integral
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    ELSE !no collisions, do single integral
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_nocoll_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonct_hcubaturep_newt
  
  SUBROUTINE calcfonct_pcubaturep_newt(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !collisions, do double integral
     IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrapped_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    ELSE !no collisions, do single integral
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trapped_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrapped_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrapped_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonct_pcubaturep_newt
  
  SUBROUTINE calcfonctrot_hcubature(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdata(6:7) = -rkuplim
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, totalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
    
    
  
  END SUBROUTINE calcfonctrot_hcubature
  SUBROUTINE calcfonctrot_pcubature(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdata(6:7) = -rkuplim
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, totalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotalrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
  
  END SUBROUTINE calcfonctrot_pcubature
  
  SUBROUTINE calcfonctrot_hcubaturec(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
      ifailloc = hcubature(fdim, passingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = hcubature(fdim, ipassingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonctrot_hcubaturec
  
  SUBROUTINE calcfonctrot_pcubaturec(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
      ifailloc = pcubature(fdim, passingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = pcubature(fdim, ipassingrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonctrot_pcubaturec
  
  SUBROUTINE calcfonctrot_hcubaturec2(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
    
      ifailloc = hcubature(fdim, passing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = hcubature(fdim, passing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = hcubature(fdim, rpassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = hcubature(fdim, ipassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = hcubature(fdim, ipassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonctrot_hcubaturec2
  
  SUBROUTINE calcfonctrot_pcubaturec2(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = pcubature(fdim, passing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = pcubature(fdim, passing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = pcubature(fdim, rpassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = pcubature(fdim, ipassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = pcubature(fdim, ipassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonctrot_pcubaturec2
  
  SUBROUTINE calcfonctrot_hcubaturep(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral   
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
      
    ELSE !collisionless, single integrals
    
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  END SUBROUTINE calcfonctrot_hcubaturep
  
  SUBROUTINE calcfonctrot_pcubaturep(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu
    
    !Set integration limits. (1) for kappa and (2) for v (for the electrons)
    a(1) = 0.0d0 + barelyavoid
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    fdata(2) = vuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral   
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrappedrot_cubature, ndim, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
      
    ELSE !collisionless, single integrals
    
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc, reqrelacc, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonctrot_pcubaturep
  
  SUBROUTINE calcfonctrot_hcubature_newt(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdata(6:7) = -rkuplim
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, totalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
      
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = hcubature(fdim, total_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = hcubature(fdim, rtotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = hcubature(fdim, itotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF

    END IF
  
  END SUBROUTINE calcfonctrot_hcubature_newt
  SUBROUTINE calcfonctrot_pcubature_newt(p, nu, omega, fonct_total)
    !-----------------------------------------------------------
    ! Calculate all particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonct_total
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdata(6:7) = -rkuplim
    fdata(8:9) = rkuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !Collisional simulation, do double integral
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, totalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotalrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
      
    ELSE !no collisions
      IF(norm.EQ.2) THEN
        ifailloc = pcubature(fdim, total_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
        ifailloc = pcubature(fdim, rtotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = intout(1)
        
        ifailloc = pcubature(fdim, itotal_nocollrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        intout = intout * ABS(Ac(p))
        fonct_total = fonct_total + ci * intout(1)
      END IF
    END IF
  
  END SUBROUTINE calcfonctrot_pcubature_newt
  
  SUBROUTINE calcfonctrot_hcubaturec_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF(norm.EQ.2) THEN
      ifailloc = hcubature(fdim, passingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = hcubature(fdim, ipassingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonctrot_hcubaturec_newt
  
  SUBROUTINE calcfonctrot_pcubaturec_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! 
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
      ifailloc = pcubature(fdim, passingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) 
      
      ifailloc = pcubature(fdim, ipassingrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)  
      
    END IF
  
  END SUBROUTINE calcfonctrot_pcubaturec_newt
  
  SUBROUTINE calcfonctrot_hcubaturec2_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = hcubature(fdim, passing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = hcubature(fdim, passing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = hcubature(fdim, rpassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = hcubature(fdim, rpassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = hcubature(fdim, ipassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = hcubature(fdim, ipassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF
  
  END SUBROUTINE calcfonctrot_hcubaturec2_newt
  
  SUBROUTINE calcfonctrot_pcubaturec2_newt(p, nu, omega, fonctc)
    !-----------------------------------------------------------
    ! Calculate passing particle integrals
    ! ions are separate from electrons
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctc
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFkr = omega
    pFkr = p
    nuFkr = nu
    
    !Set integration limits. (1) for kstar and (2) for rstar
    a(:) = -rkuplim
    b(:) = rkuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF(norm.EQ.2) THEN
    
      ifailloc = pcubature(fdim, passing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1) + ci * intout(2)  
      
      ifailloc = pcubature(fdim, passing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1) + ci * intout(2)  
      
    ELSE IF(norm.EQ.1) THEN
      ifailloc = pcubature(fdim, rpassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = intout(1)
      
      ifailloc = pcubature(fdim, rpassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + intout(1)
      
      ifailloc = pcubature(fdim, ipassing_ionsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
      ifailloc = pcubature(fdim, ipassing_electronsrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
      
      IF (ifailloc /= 0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
          &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
      END IF
      
      intout = intout * ABS(Ac(p))
      fonctc = fonctc + ci * intout(1)
      
    END IF 
  
  END SUBROUTINE calcfonctrot_pcubaturec2_newt
  
  SUBROUTINE calcfonctrot_hcubaturep_newt(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
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
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    
    IF ( ABS(coll_flag) > epsD) THEN !collisions, do double integral
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    ELSE !no collisions, do single integral
      IF(norm.EQ.2) THEN
      
        ifailloc = hcubature(fdim, trapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = hcubature(fdim, rtrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = hcubature(fdim, itrapped_nocollrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of hcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonctrot_hcubaturep_newt
  
  SUBROUTINE calcfonctrot_pcubaturep_newt(p, nu, omega, fonctp)
    !-----------------------------------------------------------
    ! Calculate trapped particle integrals
    ! Includes collisions for the trapped electron terms
    !-----------------------------------------------------------   
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonctp
    
    REAL(KIND=DBL), DIMENSION(2) :: a, b, acc
    REAL(KIND=DBL), DIMENSION(9) :: fdata
    INTEGER :: ifailloc, fdim
    
    REAL(KIND=DBL), DIMENSION(nf) :: intout
    
    !set the omega and p to be seen by all integration subroutines
    omFFk = omega
    pFFk = p
    nuFFk = nu
    
    !Set integration limits. (1) for kappa and (2) for v (for the electrons)
    a(1) = 0.0d0 + barelyavoid
    b(1) = 1.0d0 - barelyavoid
    a(2) = 0.0d0
    b(2) = vuplim
    
    !Set scale factor inside fdata
    fdata(1) = Ac(p)
    fdata(2) = vuplim
    
    !set the convergence criterion to be the L_2 norm
    fdim = 1
    
    IF(norm.EQ.2) fdim = ndim
    
    IF ( ABS(coll_flag) > epsD) THEN !collisions, do double integral
     IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrappedrot_cubature, ndim, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    ELSE !no collisions, do single integral
      IF(norm.EQ.2) THEN
      
        ifailloc = pcubature(fdim, trappedrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1) + ci * intout(2)  
        
      ELSE IF(norm.EQ.1) THEN
      
        ifailloc = pcubature(fdim, rtrappedrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = intout(1)
        
        ifailloc = pcubature(fdim, itrappedrot_cubature, 1, a, b, maxpts, reqabsacc_newt, reqrelacc_newt, norm, intout, acc, fdata=fdata)
        
        IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of pcubature integration in mod_fonct at p=',p,', nu=',nu,' omega=(',REAL(omega),',',AIMAG(omega),')'
        END IF
        
        intout = intout * ABS(Ac(p))
        fonctp = fonctp + ci * intout(1)     
      END IF
    END IF
  
  END SUBROUTINE calcfonctrot_pcubaturep_newt

END MODULE mod_fonct
