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
  
END MODULE calltotalints