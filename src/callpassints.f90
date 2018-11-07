MODULE callpassints
  USE kind
  USE passints
  USE datmat, ONLY : nions, ion_type, ion, nuFkr
  IMPLICIT NONE

CONTAINS

  INTEGER FUNCTION passing_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passing_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = 4.*Fkstarrstare(ndim, XY, 1)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passing_cubature = 0
    
    
  END FUNCTION passing_cubature
  
  INTEGER FUNCTION rpassing_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rpassing_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = 4.*Fkstarrstare(ndim, XY, 1)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassing_cubature = 0
    
    
  END FUNCTION rpassing_cubature
  INTEGER FUNCTION ipassing_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      ipassing_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = 4.*Fkstarrstare(ndim, XY, 1)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassing_cubature = 0
    
    
  END FUNCTION ipassing_cubature
  
  INTEGER FUNCTION passing_ions_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passing_ions_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passing_ions_cubature = 0
    
    
  END FUNCTION passing_ions_cubature
  
  INTEGER FUNCTION rpassing_ions_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rpassing_ions_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassing_ions_cubature = 0
    
    
  END FUNCTION rpassing_ions_cubature
  
  INTEGER FUNCTION ipassing_ions_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      ipassing_ions_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassing_ions_cubature = 0
    
    
  END FUNCTION ipassing_ions_cubature
  
  INTEGER FUNCTION passing_electrons_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passing_electrons_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = 4.*Fkstarrstare(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passing_electrons_cubature = 0
    
    
  END FUNCTION passing_electrons_cubature
  
  INTEGER FUNCTION rpassing_electrons_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      rpassing_electrons_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = 4.*Fkstarrstare(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassing_electrons_cubature = 0
    
    
  END FUNCTION rpassing_electrons_cubature
  
  INTEGER FUNCTION ipassing_electrons_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      ipassing_electrons_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstare(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassing_electrons_cubature = 0
    
    
  END FUNCTION ipassing_electrons_cubature
  
  INTEGER FUNCTION passingrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passingrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passingrot_cubature = 0
    
    
  END FUNCTION passingrot_cubature
  
  INTEGER FUNCTION rpassingrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rpassingrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassingrot_cubature = 0
    
    
  END FUNCTION rpassingrot_cubature
  INTEGER FUNCTION ipassingrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      ipassingrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassingrot_cubature = 0
    
    
  END FUNCTION ipassingrot_cubature
  
  INTEGER FUNCTION passing_ionsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passing_ionsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passing_ionsrot_cubature = 0
    
    
  END FUNCTION passing_ionsrot_cubature
  
  INTEGER FUNCTION rpassing_ionsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rpassing_ionsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassing_ionsrot_cubature = 0
    
    
  END FUNCTION rpassing_ionsrot_cubature
  
  INTEGER FUNCTION ipassing_ionsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      ipassing_ionsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = CMPLX(0., 0.)
    
    DO i = 1, nions
      output = output +  Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i)
    END DO    
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassing_ionsrot_cubature = 0
    
    
  END FUNCTION ipassing_ionsrot_cubature
  
  INTEGER FUNCTION passing_electronsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      passing_electronsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    passing_electronsrot_cubature = 0
    
    
  END FUNCTION passing_electronsrot_cubature
  
  INTEGER FUNCTION rpassing_electronsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      rpassing_electronsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = REAL(output)
    rpassing_electronsrot_cubature = 0
    
    
  END FUNCTION rpassing_electronsrot_cubature
    
  INTEGER FUNCTION ipassing_electronsrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      ipassing_electronsrot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    output = Fkstarrstarerot(ndim, XY, 1)  
        
    output = output/scale_
    
    fval(1) = AIMAG(output)
    ipassing_electronsrot_cubature = 0
    
    
  END FUNCTION ipassing_electronsrot_cubature

  FUNCTION Fkstarrstar_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    INTEGER :: i
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstar_cub
    REAL(KIND=DBL) :: intsum(nf)
    !NOTE THE FACTOR 4 BECAUSE WE ASSUME HERE SYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum(:)=0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum(1) = 4.*REAL  ( Fkstarrstare(ndim,xy,1) )
       intsum(2) = 4.*AIMAG ( Fkstarrstare(ndim,xy,1) )
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum(1)=intsum(1)+4.*REAL ( Fkstarrstari(ndim,xy,1,i) )
          intsum(2)=intsum(2)+4.*AIMAG ( Fkstarrstari(ndim,xy,1,i) )
       ENDIF
    ENDDO
    Fkstarrstar_cub(1) = intsum(1)
    Fkstarrstar_cub(2) = intsum(2)
   
  END FUNCTION Fkstarrstar_cub

  FUNCTION rFkstarrstar(nf,xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER :: i
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstar
    REAL(KIND=DBL) :: intsum
    !NOTE THE FACTOR 4 BECAUSE WE ASSUME HERE SYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum = 0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum = 4.*REAL  ( Fkstarrstare(ndim,xy,1) )
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum + 4.*REAL ( Fkstarrstari(ndim,xy,1,i) )*ninorm(pFkr,i) !unnormalise coefi here
       ENDIF
    ENDDO
    rFkstarrstar = intsum

  END FUNCTION rFkstarrstar

  FUNCTION iFkstarrstar(nf,xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER :: i
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstar
    REAL(KIND=DBL) :: intsum
    !NOTE THE FACTOR 4 BECAUSE WE ASSUME HERE SYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum = 0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum = 4.*AIMAG ( Fkstarrstare(ndim,xy,1) )
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum + 4.*AIMAG ( Fkstarrstari(ndim,xy,1,i) )*ninorm(pFkr,i) !unnormalise coefi here
       ENDIF
    ENDDO
    iFkstarrstar = intsum
   
  END FUNCTION iFkstarrstar

  FUNCTION Fkstarrstare_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstare_cub
    Fkstarrstare_cub(1) = REAL  ( Fkstarrstare(ndim,xy,1) )
    Fkstarrstare_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,1) )
  END FUNCTION Fkstarrstare_cub

  FUNCTION rFkstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstare
    rFkstarrstare = REAL  ( Fkstarrstare(ndim,xy,1) )
  END FUNCTION rFkstarrstare

  FUNCTION iFkstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstare
    iFkstarrstare = AIMAG (  Fkstarrstare(ndim,xy,1) )
  END FUNCTION iFkstarrstare
  
  INTEGER FUNCTION iFkstarrstare_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstare_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,1))
    
    iFkstarrstare_cubature = 0
    
  END FUNCTION iFkstarrstare_cubature

  FUNCTION Fkstarrstargte_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the elecron passing particle integrand (ATe term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargte_cub
    Fkstarrstargte_cub(1) = REAL  ( Fkstarrstare(ndim,xy,2) )
    Fkstarrstargte_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,2) )
  END FUNCTION Fkstarrstargte_cub

  FUNCTION rFkstarrstargte(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargte
    rFkstarrstargte = REAL  ( Fkstarrstare(ndim,xy,2) )
  END FUNCTION rFkstarrstargte

  FUNCTION iFkstarrstargte(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargte
    iFkstarrstargte = AIMAG (  Fkstarrstare(ndim,xy,2) )
  END FUNCTION iFkstarrstargte
  
  INTEGER FUNCTION iFkstarrstargte_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargte_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,2))
    
    iFkstarrstargte_cubature = 0
    
  END FUNCTION iFkstarrstargte_cubature

  FUNCTION Fkstarrstargne_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargne_cub
    Fkstarrstargne_cub(1) = REAL  ( Fkstarrstare(ndim,xy,3) )
    Fkstarrstargne_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,3) )
  END FUNCTION Fkstarrstargne_cub

  FUNCTION rFkstarrstargne(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargne
    rFkstarrstargne = REAL  ( Fkstarrstare(ndim,xy,3) )
  END FUNCTION rFkstarrstargne

  FUNCTION iFkstarrstargne(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargne
    iFkstarrstargne = AIMAG (  Fkstarrstare(ndim,xy,3) )
  END FUNCTION iFkstarrstargne
  
  INTEGER FUNCTION iFkstarrstargne_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargne_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,3))
    
    iFkstarrstargne_cubature = 0
    
  END FUNCTION iFkstarrstargne_cubature

  FUNCTION Fkstarrstarce_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron  passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarce_cub
    Fkstarrstarce_cub(1) = REAL  ( Fkstarrstare(ndim,xy,4) )
    Fkstarrstarce_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,4) )
  END FUNCTION Fkstarrstarce_cub

  FUNCTION rFkstarrstarce(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarce
    rFkstarrstarce = REAL  ( Fkstarrstare(ndim,xy,4) )
  END FUNCTION rFkstarrstarce

  FUNCTION iFkstarrstarce(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarce
    iFkstarrstarce = AIMAG (  Fkstarrstare(ndim,xy,4) )
  END FUNCTION iFkstarrstarce
  
  INTEGER FUNCTION iFkstarrstarce_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarce_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,4))
    
    iFkstarrstarce_cubature = 0
    
  END FUNCTION iFkstarrstarce_cubature

  FUNCTION Fekstarrstare_cub(nf, xy)
    !---------------------------------------------------------------------
    ! returns the electron passing energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstare_cub
    Fekstarrstare_cub(1) = REAL  ( Fkstarrstare(ndim,xy,5) )
    Fekstarrstare_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,5) )
  END FUNCTION Fekstarrstare_cub

  FUNCTION rFekstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstare
    rFekstarrstare = REAL  ( Fkstarrstare(ndim,xy,5) )
  END FUNCTION rFekstarrstare

  FUNCTION iFekstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstare
    iFekstarrstare = AIMAG (  Fkstarrstare(ndim,xy,5) )
  END FUNCTION iFekstarrstare
  
  INTEGER FUNCTION iFekstarrstare_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstare_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,5))
    
    iFekstarrstare_cubature = 0
    
  END FUNCTION iFekstarrstare_cubature

  FUNCTION Fekstarrstargte_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the elecron passing particle integrand (ATe term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargte_cub
    Fekstarrstargte_cub(1) = REAL  ( Fkstarrstare(ndim,xy,6) )
    Fekstarrstargte_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,6) )
  END FUNCTION Fekstarrstargte_cub

  FUNCTION rFekstarrstargte(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargte
    rFekstarrstargte = REAL  ( Fkstarrstare(ndim,xy,6) )
  END FUNCTION rFekstarrstargte

  FUNCTION iFekstarrstargte(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargte
    iFekstarrstargte = AIMAG (  Fkstarrstare(ndim,xy,6) )
  END FUNCTION iFekstarrstargte
  
  INTEGER FUNCTION iFekstarrstargte_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargte_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,6))
    
    iFekstarrstargte_cubature = 0
    
  END FUNCTION iFekstarrstargte_cubature

  FUNCTION Fekstarrstargne_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargne_cub
    Fekstarrstargne_cub(1) = REAL  ( Fkstarrstare(ndim,xy,7) )
    Fekstarrstargne_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,7) )
  END FUNCTION Fekstarrstargne_cub

  FUNCTION rFekstarrstargne(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargne
    rFekstarrstargne = REAL  ( Fkstarrstare(ndim,xy,7) )
  END FUNCTION rFekstarrstargne

  FUNCTION iFekstarrstargne(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargne
    iFekstarrstargne = AIMAG (  Fkstarrstare(ndim,xy,7) )
  END FUNCTION iFekstarrstargne
  
  INTEGER FUNCTION iFekstarrstargne_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargne_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,7))
    
    iFekstarrstargne_cubature = 0
    
  END FUNCTION iFekstarrstargne_cubature

  FUNCTION Fekstarrstarce_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron  passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarce_cub
    Fekstarrstarce_cub(1) = REAL  ( Fkstarrstare(ndim,xy,8) )
    Fekstarrstarce_cub(2) = AIMAG (  Fkstarrstare(ndim,xy,8) )
  END FUNCTION Fekstarrstarce_cub

  FUNCTION rFekstarrstarce(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarce
    rFekstarrstarce = REAL  ( Fkstarrstare(ndim,xy,8) )
  END FUNCTION rFekstarrstarce

  FUNCTION iFekstarrstarce(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarce
    iFekstarrstarce = AIMAG (  Fkstarrstare(ndim,xy,8) )
  END FUNCTION iFekstarrstarce
  
  INTEGER FUNCTION iFekstarrstarce_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarce_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstare(ndim,XY,8))
    
    iFekstarrstarce_cubature = 0
    
  END FUNCTION iFekstarrstarce_cubature

  FUNCTION Fkstarrstari_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstari_cub
    Fkstarrstari_cub(1) = REAL  ( Fkstarrstari(ndim,xy,1,ion) )
    Fkstarrstari_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,1,ion) )
  END FUNCTION Fkstarrstari_cub

  FUNCTION rFkstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstari
    rFkstarrstari = REAL  ( Fkstarrstari(ndim,xy,1,ion) )
  END FUNCTION rFkstarrstari

  FUNCTION iFkstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstari
    iFkstarrstari = AIMAG ( Fkstarrstari(ndim,xy,1,ion) )
  END FUNCTION iFkstarrstari
  
  INTEGER FUNCTION iFkstarrstari_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstari_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,1,ion))
    
    iFkstarrstari_cubature = 0
    
  END FUNCTION iFkstarrstari_cubature

  FUNCTION Fkstarrstargti_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (ATi term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargti_cub
    Fkstarrstargti_cub(1) = REAL  ( Fkstarrstari(ndim,xy,2,ion) )
    Fkstarrstargti_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,2,ion) )
  END FUNCTION Fkstarrstargti_cub

  FUNCTION rFkstarrstargti(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargti
    rFkstarrstargti = REAL  ( Fkstarrstari(ndim,xy,2,ion) )
  END FUNCTION rFkstarrstargti

  FUNCTION iFkstarrstargti(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargti
    iFkstarrstargti = AIMAG ( Fkstarrstari(ndim,xy,2,ion) )
  END FUNCTION iFkstarrstargti
  
  INTEGER FUNCTION iFkstarrstargti_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargti_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,2,ion))
    
    iFkstarrstargti_cubature = 0
    
  END FUNCTION iFkstarrstargti_cubature

  FUNCTION Fkstarrstargni_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Ani term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargni_cub
    Fkstarrstargni_cub(1) = REAL  ( Fkstarrstari(ndim,xy,3,ion) )
    Fkstarrstargni_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,3,ion) )
  END FUNCTION Fkstarrstargni_cub

  FUNCTION rFkstarrstargni(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargni
    rFkstarrstargni = REAL  ( Fkstarrstari(ndim,xy,3,ion) )
  END FUNCTION rFkstarrstargni

  FUNCTION iFkstarrstargni(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargni
    iFkstarrstargni = AIMAG ( Fkstarrstari(ndim,xy,3,ion) )
  END FUNCTION iFkstarrstargni
  
  INTEGER FUNCTION iFkstarrstargni_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargni_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,3,ion))
    
    iFkstarrstargni_cubature = 0
    
  END FUNCTION iFkstarrstargni_cubature

  FUNCTION Fkstarrstarci_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarci_cub
    Fkstarrstarci_cub(1) = REAL  ( Fkstarrstari(ndim,xy,4,ion) )
    Fkstarrstarci_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,4,ion) )
  END FUNCTION Fkstarrstarci_cub

  FUNCTION rFkstarrstarci(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarci
    rFkstarrstarci = REAL  ( Fkstarrstari(ndim,xy,4,ion) )
  END FUNCTION rFkstarrstarci

  FUNCTION iFkstarrstarci(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarci
    iFkstarrstarci = AIMAG ( Fkstarrstari(ndim,xy,4,ion) )
  END FUNCTION iFkstarrstarci
  
  INTEGER FUNCTION iFkstarrstarci_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarci_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,4,ion))
    
    iFkstarrstarci_cubature = 0
    
  END FUNCTION iFkstarrstarci_cubature

  FUNCTION Fekstarrstari_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstari_cub
    Fekstarrstari_cub(1) = REAL  ( Fkstarrstari(ndim,xy,5,ion) )
    Fekstarrstari_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,5,ion) )
  END FUNCTION Fekstarrstari_cub

  FUNCTION rFekstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstari
    rFekstarrstari = REAL  ( Fkstarrstari(ndim,xy,5,ion) )
  END FUNCTION rFekstarrstari

  FUNCTION iFekstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstari
    iFekstarrstari = AIMAG ( Fkstarrstari(ndim,xy,5,ion) )
  END FUNCTION iFekstarrstari
  
  INTEGER FUNCTION iFekstarrstari_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstari_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,5,ion))
    
    iFekstarrstari_cubature = 0
    
  END FUNCTION iFekstarrstari_cubature

  FUNCTION Fekstarrstargti_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (ATi term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargti_cub
    Fekstarrstargti_cub(1) = REAL  ( Fkstarrstari(ndim,xy,6,ion) )
    Fekstarrstargti_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,6,ion) )
  END FUNCTION Fekstarrstargti_cub

  FUNCTION rFekstarrstargti(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargti
    rFekstarrstargti = REAL  ( Fkstarrstari(ndim,xy,6,ion) )
  END FUNCTION rFekstarrstargti

  FUNCTION iFekstarrstargti(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargti
    iFekstarrstargti = AIMAG ( Fkstarrstari(ndim,xy,6,ion) )
  END FUNCTION iFekstarrstargti

  INTEGER FUNCTION iFekstarrstargti_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargti_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,6,ion))
    
    iFekstarrstargti_cubature = 0
    
  END FUNCTION iFekstarrstargti_cubature
  
  FUNCTION Fekstarrstargni_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Ani term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargni_cub
    Fekstarrstargni_cub(1) = REAL  ( Fkstarrstari(ndim,xy,7,ion) )
    Fekstarrstargni_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,7,ion) )
  END FUNCTION Fekstarrstargni_cub

  FUNCTION rFekstarrstargni(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargni
    rFekstarrstargni = REAL  ( Fkstarrstari(ndim,xy,7,ion) )
  END FUNCTION rFekstarrstargni

  FUNCTION iFekstarrstargni(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargni
    iFekstarrstargni = AIMAG ( Fkstarrstari(ndim,xy,7,ion) )
  END FUNCTION iFekstarrstargni
  
  INTEGER FUNCTION iFekstarrstargni_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargni_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,7,ion))
    
    iFekstarrstargni_cubature = 0
    
  END FUNCTION iFekstarrstargni_cubature

  FUNCTION Fekstarrstarci_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarci_cub
    Fekstarrstarci_cub(1) = REAL  ( Fkstarrstari(ndim,xy,8,ion) )
    Fekstarrstarci_cub(2) = AIMAG (  Fkstarrstari(ndim,xy,8,ion) )
  END FUNCTION Fekstarrstarci_cub

  FUNCTION rFekstarrstarci(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarci
    rFekstarrstarci = REAL  ( Fkstarrstari(ndim,xy,8,ion) )
  END FUNCTION rFekstarrstarci

  FUNCTION iFekstarrstarci(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarci
    iFekstarrstarci = AIMAG ( Fkstarrstari(ndim,xy,8,ion) )
  END FUNCTION iFekstarrstarci
  
  INTEGER FUNCTION iFekstarrstarci_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarci_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstari(ndim,XY,8,ion))
    
    iFekstarrstarci_cubature = 0
    
  END FUNCTION iFekstarrstarci_cubature

!****************************************************************************************************
! with rotation
!****************************************************************************************************

  FUNCTION Fkstarrstarrot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    INTEGER :: i
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarrot_cub
    REAL(KIND=DBL) :: intsum(nf)
    !NOTE THE FACTOR 1 not 4 BECAUSE ASYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum(:)=0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum(1) = 1.*REAL  ( Fkstarrstarerot(ndim,xy,1) )
       intsum(2) = 1.*AIMAG ( Fkstarrstarerot(ndim,xy,1) )
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum(1)=intsum(1)+1.*REAL ( Fkstarrstarirot(ndim,xy,1,i) )
          intsum(2)=intsum(2)+1.*AIMAG ( Fkstarrstarirot(ndim,xy,1,i) )
       ENDIF
    ENDDO
    Fkstarrstarrot_cub(1) = intsum(1)
    Fkstarrstarrot_cub(2) = intsum(2)
   
  END FUNCTION Fkstarrstarrot_cub

  FUNCTION rFkstarrstarrot(nf,xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER :: i
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarrot
    REAL(KIND=DBL) :: intsum
    !NOTE THE FACTOR 1 not 4 BECAUSE ASYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum = 0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum = 1.*REAL  ( Fkstarrstarerot(ndim,xy,1) )
!       WRITE(*,*) 'elec 1: ',Fkstarrstarerot(ndim,xy,1) 
!       WRITE(*,*) 'elec 2: ',Fkstarrstare(ndim,xy,1) 
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum + 1.*REAL ( Fkstarrstarirot(ndim,xy,1,i) )*ninorm(pFkr,i) !unnormalise coefi here  
       ENDIF
    ENDDO
    rFkstarrstarrot = intsum
    
  END FUNCTION rFkstarrstarrot

  FUNCTION iFkstarrstarrot(nf,xy)
    !---------------------------------------------------------------------
    ! Returns the total passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER :: i
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarrot
    REAL(KIND=DBL) :: intsum
    !NOTE THE FACTOR 1 not 4 BECAUSE ASYMMETRIC INTEGRALS
    !Only include electron integral if el_type == 1

    intsum = 0
    IF ( (el_type == 1) .OR. ( (el_type == 3) .AND. (ETG_flag(nuFkr) .EQV. .TRUE.) ) )  THEN
       intsum = 1.*AIMAG ( Fkstarrstarerot(ndim,xy,1) )
    ENDIF

    DO i = 1,nions
       IF ( (ion_type(pFkr,i) == 1) .AND. (ETG_flag(nuFkr) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum + 1.*AIMAG ( Fkstarrstarirot(ndim,xy,1,i) )*ninorm(pFkr,i) !unnormalise coefi here
       ENDIF
    ENDDO
    iFkstarrstarrot = intsum
   
  END FUNCTION iFkstarrstarrot
  
  

  FUNCTION Fkstarrstarerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarerot_cub
    Fkstarrstarerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,1) )
    Fkstarrstarerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,1) )
  END FUNCTION Fkstarrstarerot_cub

  FUNCTION rFkstarrstarerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarerot
    rFkstarrstarerot = REAL  ( Fkstarrstarerot(ndim,xy,1) )
  END FUNCTION rFkstarrstarerot

  FUNCTION iFkstarrstarerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarerot
    iFkstarrstarerot = AIMAG (  Fkstarrstarerot(ndim,xy,1) )
  END FUNCTION iFkstarrstarerot
  
  INTEGER FUNCTION iFkstarrstarerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,1))
    
    iFkstarrstarerot_cubature = 0
    
  END FUNCTION iFkstarrstarerot_cubature

  FUNCTION Fkstarrstargterot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the elecron passing particle integrand (ATe term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargterot_cub
    Fkstarrstargterot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,2) )
    Fkstarrstargterot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,2) )
  END FUNCTION Fkstarrstargterot_cub

  FUNCTION rFkstarrstargterot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargterot
    rFkstarrstargterot = REAL  ( Fkstarrstarerot(ndim,xy,2) )
  END FUNCTION rFkstarrstargterot

  FUNCTION iFkstarrstargterot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargterot
    iFkstarrstargterot = AIMAG (  Fkstarrstarerot(ndim,xy,2) )
  END FUNCTION iFkstarrstargterot
  
  INTEGER FUNCTION iFkstarrstargterot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargterot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,2))
    
    iFkstarrstargterot_cubature = 0
    
  END FUNCTION iFkstarrstargterot_cubature

  FUNCTION Fkstarrstargnerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargnerot_cub
    Fkstarrstargnerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,3) )
    Fkstarrstargnerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,3) )
  END FUNCTION Fkstarrstargnerot_cub

  FUNCTION rFkstarrstargnerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargnerot
    rFkstarrstargnerot = REAL  ( Fkstarrstarerot(ndim,xy,3) )
  END FUNCTION rFkstarrstargnerot

  FUNCTION iFkstarrstargnerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargnerot
    iFkstarrstargnerot = AIMAG (  Fkstarrstarerot(ndim,xy,3) )
  END FUNCTION iFkstarrstargnerot
  
  INTEGER FUNCTION iFkstarrstargnerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargnerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,3))
    
    iFkstarrstargnerot_cubature = 0
    
  END FUNCTION iFkstarrstargnerot_cubature

  FUNCTION Fkstarrstarcerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron  passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarcerot_cub
    Fkstarrstarcerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,4) )
    Fkstarrstarcerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,4) )
  END FUNCTION Fkstarrstarcerot_cub

  FUNCTION rFkstarrstarcerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarcerot
    rFkstarrstarcerot = REAL  ( Fkstarrstarerot(ndim,xy,4) )
  END FUNCTION rFkstarrstarcerot

  FUNCTION iFkstarrstarcerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarcerot
    iFkstarrstarcerot = AIMAG (  Fkstarrstarerot(ndim,xy,4) )
  END FUNCTION iFkstarrstarcerot
  
  INTEGER FUNCTION iFkstarrstarcerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarcerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,4))
    
    iFkstarrstarcerot_cubature = 0
    
  END FUNCTION iFkstarrstarcerot_cubature

  FUNCTION Fekstarrstarerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! returns the electron passing energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarerot_cub
    Fekstarrstarerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,5) )
    Fekstarrstarerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,5) )
  END FUNCTION Fekstarrstarerot_cub

  FUNCTION rFekstarrstarerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarerot
    rFekstarrstarerot = REAL  ( Fkstarrstarerot(ndim,xy,5) )
  END FUNCTION rFekstarrstarerot

  FUNCTION iFekstarrstarerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarerot
    iFekstarrstarerot = AIMAG (  Fkstarrstarerot(ndim,xy,5) )
  END FUNCTION iFekstarrstarerot
  
  INTEGER FUNCTION iFekstarrstarerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,5))
    
    iFekstarrstarerot_cubature = 0
    
  END FUNCTION iFekstarrstarerot_cubature

  FUNCTION Fekstarrstargterot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the elecron passing particle integrand (ATe term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargterot_cub
    Fekstarrstargterot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,6) )
    Fekstarrstargterot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,6) )
  END FUNCTION Fekstarrstargterot_cub

  FUNCTION rFekstarrstargterot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargterot
    rFekstarrstargterot = REAL  ( Fkstarrstarerot(ndim,xy,6) )
  END FUNCTION rFekstarrstargterot

  FUNCTION iFekstarrstargterot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargterot
    iFekstarrstargterot = AIMAG (  Fkstarrstarerot(ndim,xy,6) )
  END FUNCTION iFekstarrstargterot
  
  INTEGER FUNCTION iFekstarrstargterot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargterot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,6))
    
    iFekstarrstargterot_cubature = 0
    
  END FUNCTION iFekstarrstargterot_cubature

  FUNCTION Fekstarrstargnerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargnerot_cub
    Fekstarrstargnerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,7) )
    Fekstarrstargnerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,7) )
  END FUNCTION Fekstarrstargnerot_cub

  FUNCTION rFekstarrstargnerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargnerot
    rFekstarrstargnerot = REAL  ( Fkstarrstarerot(ndim,xy,7) )
  END FUNCTION rFekstarrstargnerot

  FUNCTION iFekstarrstargnerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargnerot
    iFekstarrstargnerot = AIMAG (  Fkstarrstarerot(ndim,xy,7) )
  END FUNCTION iFekstarrstargnerot
  
  INTEGER FUNCTION iFekstarrstargnerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargnerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,7))
    
    iFekstarrstargnerot_cubature = 0
    
  END FUNCTION iFekstarrstargnerot_cubature

  FUNCTION Fekstarrstarcerot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron  passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarcerot_cub
    Fekstarrstarcerot_cub(1) = REAL  ( Fkstarrstarerot(ndim,xy,8) )
    Fekstarrstarcerot_cub(2) = AIMAG (  Fkstarrstarerot(ndim,xy,8) )
  END FUNCTION Fekstarrstarcerot_cub

  FUNCTION rFekstarrstarcerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarcerot
    rFekstarrstarcerot = REAL  ( Fkstarrstarerot(ndim,xy,8) )
  END FUNCTION rFekstarrstarcerot

  FUNCTION iFekstarrstarcerot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarcerot
    iFekstarrstarcerot = AIMAG (  Fkstarrstarerot(ndim,xy,8) )
  END FUNCTION iFekstarrstarcerot
  
  INTEGER FUNCTION iFekstarrstarcerot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarcerot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarerot(ndim,XY,8))
    
    iFekstarrstarcerot_cubature = 0
    
  END FUNCTION iFekstarrstarcerot_cubature

  FUNCTION Fkstarrstarirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarirot_cub
    Fkstarrstarirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,1,ion) )
    Fkstarrstarirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,1,ion) )
  END FUNCTION Fkstarrstarirot_cub

  FUNCTION rFkstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarirot
    rFkstarrstarirot = REAL  ( Fkstarrstarirot(ndim,xy,1,ion) )
  END FUNCTION rFkstarrstarirot

  FUNCTION iFkstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarirot
    iFkstarrstarirot = AIMAG ( Fkstarrstarirot(ndim,xy,1,ion) )
  END FUNCTION iFkstarrstarirot
  
  INTEGER FUNCTION iFkstarrstarirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,1, ion))
    
    iFkstarrstarirot_cubature = 0
    
  END FUNCTION iFkstarrstarirot_cubature

  FUNCTION Fkstarrstargtirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (ATi term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargtirot_cub
    Fkstarrstargtirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,2,ion) )
    Fkstarrstargtirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,2,ion) )
  END FUNCTION Fkstarrstargtirot_cub

  FUNCTION rFkstarrstargtirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargtirot
    rFkstarrstargtirot = REAL  ( Fkstarrstarirot(ndim,xy,2,ion) )
  END FUNCTION rFkstarrstargtirot

  FUNCTION iFkstarrstargtirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargtirot
    iFkstarrstargtirot = AIMAG ( Fkstarrstarirot(ndim,xy,2,ion) )
  END FUNCTION iFkstarrstargtirot
  
  INTEGER FUNCTION iFkstarrstargtirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargtirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,2,ion))
    
    iFkstarrstargtirot_cubature = 0
    
  END FUNCTION iFkstarrstargtirot_cubature

  FUNCTION Fkstarrstargnirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Ani term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstargnirot_cub
    Fkstarrstargnirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,3,ion) )
    Fkstarrstargnirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,3,ion) )
  END FUNCTION Fkstarrstargnirot_cub

  FUNCTION rFkstarrstargnirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstargnirot
    rFkstarrstargnirot = REAL  ( Fkstarrstarirot(ndim,xy,3,ion) )
  END FUNCTION rFkstarrstargnirot

  FUNCTION iFkstarrstargnirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstargnirot
    iFkstarrstargnirot = AIMAG ( Fkstarrstarirot(ndim,xy,3,ion) )
  END FUNCTION iFkstarrstargnirot
  
  INTEGER FUNCTION iFkstarrstargnirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstargnirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,3,ion))
    
    iFkstarrstargnirot_cubature = 0
    
  END FUNCTION iFkstarrstargnirot_cubature

  FUNCTION Fkstarrstarcirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarcirot_cub
    Fkstarrstarcirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,4,ion) )
    Fkstarrstarcirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,4,ion) )
  END FUNCTION Fkstarrstarcirot_cub

  FUNCTION rFkstarrstarcirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarcirot
    rFkstarrstarcirot = REAL  ( Fkstarrstarirot(ndim,xy,4,ion) )
  END FUNCTION rFkstarrstarcirot

  FUNCTION iFkstarrstarcirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarcirot
    iFkstarrstarcirot = AIMAG ( Fkstarrstarirot(ndim,xy,4,ion) )
  END FUNCTION iFkstarrstarcirot
  
  INTEGER FUNCTION iFkstarrstarcirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarcirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,4,ion))
    
    iFkstarrstarcirot_cubature = 0
    
  END FUNCTION iFkstarrstarcirot_cubature

  FUNCTION Fekstarrstarirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarirot_cub
    Fekstarrstarirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,5,ion) )
    Fekstarrstarirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,5,ion) )
  END FUNCTION Fekstarrstarirot_cub

  FUNCTION rFekstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarirot
    rFekstarrstarirot = REAL  ( Fkstarrstarirot(ndim,xy,5,ion) )
  END FUNCTION rFekstarrstarirot

  FUNCTION iFekstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarirot
    iFekstarrstarirot = AIMAG ( Fkstarrstarirot(ndim,xy,5,ion) )
  END FUNCTION iFekstarrstarirot
  
  INTEGER FUNCTION iFekstarrstarirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,5,ion))
    
    iFekstarrstarirot_cubature = 0
    
  END FUNCTION iFekstarrstarirot_cubature

  FUNCTION Fekstarrstargtirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (ATi term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargtirot_cub
    Fekstarrstargtirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,6,ion) )
    Fekstarrstargtirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,6,ion) )
  END FUNCTION Fekstarrstargtirot_cub

  FUNCTION rFekstarrstargtirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargtirot
    rFekstarrstargtirot = REAL  ( Fkstarrstarirot(ndim,xy,6,ion) )
  END FUNCTION rFekstarrstargtirot

  FUNCTION iFekstarrstargtirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargtirot
    iFekstarrstargtirot = AIMAG ( Fkstarrstarirot(ndim,xy,6,ion) )
  END FUNCTION iFekstarrstargtirot
  
  INTEGER FUNCTION iFekstarrstargtirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargtirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,6,ion))
    
    iFekstarrstargtirot_cubature = 0
    
  END FUNCTION iFekstarrstargtirot_cubature

  FUNCTION Fekstarrstargnirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Ani term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstargnirot_cub
    Fekstarrstargnirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,7,ion) )
    Fekstarrstargnirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,7,ion) )
  END FUNCTION Fekstarrstargnirot_cub

  FUNCTION rFekstarrstargnirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstargnirot
    rFekstarrstargnirot = REAL  ( Fkstarrstarirot(ndim,xy,7,ion) )
  END FUNCTION rFekstarrstargnirot

  FUNCTION iFekstarrstargnirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstargnirot
    iFekstarrstargnirot = AIMAG ( Fkstarrstarirot(ndim,xy,7,ion) )
  END FUNCTION iFekstarrstargnirot
  
  INTEGER FUNCTION iFekstarrstargnirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstargnirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,7,ion))
    
    iFekstarrstargnirot_cubature = 0
    
  END FUNCTION iFekstarrstargnirot_cubature

  FUNCTION Fekstarrstarcirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarcirot_cub
    Fekstarrstarcirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,8,ion) )
    Fekstarrstarcirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,8,ion) )
  END FUNCTION Fekstarrstarcirot_cub

  FUNCTION rFekstarrstarcirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarcirot
    rFekstarrstarcirot = REAL  ( Fkstarrstarirot(ndim,xy,8,ion) )
  END FUNCTION rFekstarrstarcirot

  FUNCTION iFekstarrstarcirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarcirot
    iFekstarrstarcirot = AIMAG ( Fkstarrstarirot(ndim,xy,8,ion) )
  END FUNCTION iFekstarrstarcirot
  
  INTEGER FUNCTION iFekstarrstarcirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarcirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,8,ion))
    
    iFekstarrstarcirot_cubature = 0
    
  END FUNCTION iFekstarrstarcirot_cubature


  FUNCTION Fkstarrstarguirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Aui term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fkstarrstarguirot_cub
    Fkstarrstarguirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,9,ion) )
    Fkstarrstarguirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,9,ion) )
  END FUNCTION Fkstarrstarguirot_cub

  FUNCTION rFkstarrstarguirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand Aui term
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstarguirot
    rFkstarrstarguirot = REAL  ( Fkstarrstarirot(ndim,xy,9,ion) )
  END FUNCTION rFkstarrstarguirot

  FUNCTION iFkstarrstarguirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand Aui term
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstarguirot
    iFkstarrstarguirot = AIMAG ( Fkstarrstarirot(ndim,xy,9,ion) )
  END FUNCTION iFkstarrstarguirot
  
  INTEGER FUNCTION iFkstarrstarguirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFkstarrstarguirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,9,ion))
    
    iFkstarrstarguirot_cubature = 0
    
  END FUNCTION iFkstarrstarguirot_cubature

  FUNCTION Fekstarrstarguirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand (Aui term) energy
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fekstarrstarguirot_cub
    Fekstarrstarguirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,10,ion) )
    Fekstarrstarguirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,10,ion) )
  END FUNCTION Fekstarrstarguirot_cub

  FUNCTION rFekstarrstarguirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand(Aui term) energy
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFekstarrstarguirot
    rFekstarrstarguirot = REAL  ( Fkstarrstarirot(ndim,xy,10,ion) )
  END FUNCTION rFekstarrstarguirot

  FUNCTION iFekstarrstarguirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand(Aui term) energy
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFekstarrstarguirot
    iFekstarrstarguirot = AIMAG ( Fkstarrstarirot(ndim,xy,10,ion) )
  END FUNCTION iFekstarrstarguirot
  
  INTEGER FUNCTION iFekstarrstarguirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFekstarrstarguirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,10,ion))
    
    iFekstarrstarguirot_cubature = 0
    
  END FUNCTION iFekstarrstarguirot_cubature

  FUNCTION Fvkstarrstarirot_cub(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand ang momentum
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL), DIMENSION(nf) :: Fvkstarrstarirot_cub
    Fvkstarrstarirot_cub(1) = REAL  ( Fkstarrstarirot(ndim,xy,11,ion) )
    Fvkstarrstarirot_cub(2) = AIMAG (  Fkstarrstarirot(ndim,xy,11,ion) )
  END FUNCTION Fvkstarrstarirot_cub

  FUNCTION rFvkstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand ang momentum
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFvkstarrstarirot
    rFvkstarrstarirot = REAL  ( Fkstarrstarirot(ndim,xy,11,ion) )
  END FUNCTION rFvkstarrstarirot

  FUNCTION iFvkstarrstarirot(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand ang momentum
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFvkstarrstarirot
    iFvkstarrstarirot = AIMAG ( Fkstarrstarirot(ndim,xy,11,ion) )
  END FUNCTION iFvkstarrstarirot
  
  INTEGER FUNCTION iFvkstarrstarirot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      iFvkstarrstarirot_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    
    fval(1) = AIMAG(Fkstarrstarirot(ndim,XY,11,ion))
    
    iFvkstarrstarirot_cubature = 0
    
  END FUNCTION iFvkstarrstarirot_cubature

END MODULE callpassints
