MODULE calltotalints
  USE kind
  USE passints
  USE trapints
  USE datmat, ONLY : nions, ion_type, ion, nuFkr, nuFFk
  IMPLICIT NONE
  
  CONTAINS
  
  INTEGER FUNCTION total_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      total_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFke(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    total_cubature = 0
    
    
  END FUNCTION total_cubature
  
  INTEGER FUNCTION rtotal_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rtotal_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFke(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    rtotal_cubature = 0
    
    
  END FUNCTION rtotal_cubature
  
  INTEGER FUNCTION itotal_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      itotal_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFke(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = AIMAG(output)
    itotal_cubature = 0
    
    
  END FUNCTION itotal_cubature
  
  INTEGER FUNCTION total_nocoll_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      total_nocoll_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)    
    
    xx = XY(1)
    
    output = FFke_nocoll(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    total_nocoll_cubature = 0
    
    
  END FUNCTION total_nocoll_cubature
  
  INTEGER FUNCTION rtotal_nocoll_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rtotal_nocoll_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    xx = XY(1)
    
    output = FFke_nocoll(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    rtotal_nocoll_cubature = 0
    
    
  END FUNCTION rtotal_nocoll_cubature
  
  INTEGER FUNCTION itotal_nocoll_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      itotal_nocoll_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    xx = XY(1)
    
    output = FFke_nocoll(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + 4.*Fkstarrstare(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = AIMAG(output)
    itotal_nocoll_cubature = 0
    
    
  END FUNCTION itotal_nocoll_cubature
  
  INTEGER FUNCTION totalrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      totalrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFkerot(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    totalrot_cubature = 0
    
    
  END FUNCTION totalrot_cubature
  
  INTEGER FUNCTION rtotalrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rtotalrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFkerot(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    rtotalrot_cubature = 0
    
    
  END FUNCTION rtotalrot_cubature
  
  INTEGER FUNCTION itotalrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      itotalrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = FFkerot(ndim, XY, 1) * (c-a) * (d-b)
    
    xx = XY(1)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = AIMAG(output)
    itotalrot_cubature = 0
    
    
  END FUNCTION itotalrot_cubature
  
  INTEGER FUNCTION total_nocollrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      total_nocollrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)    
    
    xx = XY(1)
    
    output = FFke_nocollrot(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    total_nocollrot_cubature = 0
    
    
  END FUNCTION total_nocollrot_cubature
  
  INTEGER FUNCTION rtotal_nocollrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      rtotal_nocollrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    xx = XY(1)
    
    output = FFke_nocollrot(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    rtotal_nocollrot_cubature = 0
    
    
  END FUNCTION rtotal_nocollrot_cubature
  
  INTEGER FUNCTION itotal_nocollrot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, a, b, c, d
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.1)) THEN
      itotal_nocollrot_cubature = 1
      RETURN
    END IF
    
    ! beginning integration limits for trapped
    a = fdata(2)
    b = fdata(3)
    ! end integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    xx = XY(1)
    
    output = FFke_nocollrot(xx, 1) * (c-a)
    
    DO i = 1, nions
      output = output + FFkirot(xx, 1, i) * ninorm(pFFk, i) * (c-a)
    END DO
    
    !beginning integration limits for passing
    a = fdata(6)
    b = fdata(7)
    !end integration limits for passing
    c = fdata(8)
    d = fdata(9) 
    
    ! transform
    XY(1) = a + x(1) * (c-a)
    XY(2) = b + x(2) * (d-b)
    
    output = output + Fkstarrstarerot(ndim, XY, 1) * (c-a) * (d-b)
    DO i = 1, nions
      output = output + Fkstarrstarirot(ndim, XY, 1, i) * ninorm(pFkr, i) * (c-a) * (d-b)
    END DO    
        
    output = CMPLX(adiabatic, 0.) - output
    output = output/scale_
    
    fval(1) = AIMAG(output)
    itotal_nocollrot_cubature = 0
    
    
  END FUNCTION itotal_nocollrot_cubature
  
  
  
  
  
END MODULE calltotalints