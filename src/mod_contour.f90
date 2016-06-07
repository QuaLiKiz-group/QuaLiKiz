MODULE mod_contour
  USE kind
  USE datcal

  IMPLICIT NONE

CONTAINS 
  SUBROUTINE  limx(C,rint)
    ! -------------------------------------------------------------------
    ! Finds the value of the real part of the contour center
    ! -------------------------------------------------------------------
    ! Arguments
    COMPLEX(KIND=DBL), INTENT(IN)  :: C
    REAL(KIND=DBL),    INTENT(OUT) :: rint

    ! Variables locales
    REAL(KIND=DBL) :: RR 
    REAL(kind=DBL), DIMENSION(4) :: A
    REAL(kind=DBL), DIMENSION(3*5) :: W 
    COMPLEX(kind=DBL), DIMENSION(3) :: Z
    REAL(kind=DBL), DIMENSION(3) :: theta0, aux
    COMPLEX(KIND=DBL), DIMENSION(3) :: ctheta0
    INTEGER :: i, ifail, ndeg=3

    RR     = (AIMAG(C)-0.05) / (1.-Elli)

    A(1) = -4.*Elli*RR
    A(2) = 0.
    A(3) = (1.+Elli*3.)*RR
    A(4) = (AIMAG(C)-0.1)

    CALL RPQR79 ( ndeg, A, Z, ifail, W )
    DO i=1,ndeg
       theta0(i)  = REAL(Z(i))
       ctheta0(i) = casin(theta0(i))
    END DO

    DO i=1,ndeg   
       aux(i) = RR*REAL(COS(ctheta0(i)) + Elli*COS(3.*ctheta0(i)))
    END DO

    rint = MAXVAL(aux)
  END SUBROUTINE limx

  SUBROUTINE squircle(Centre, var, om, rint)
    ! -------------------------------------------------------------------
    ! Maps contour choice for finding the zeros with the Davies method
    ! The input contour is defined as a circle in complex plane centered at (0,0) with radius unity
    ! The mapped contour corresponds to the bounds of the actual frequencies
    ! we are searching in. This is a rectangle (based on squircle mapping)
    ! and with a minimum imaginary part mincont, to avoid difficult to resolve
    ! (in the integrals) low growth rates
    ! See formulas in S. Brunner et al, Phys. Plas. 5, 3929 (1998) and references therein
    ! -------------------------------------------------------------------
    COMPLEX(KIND=DBL), INTENT(IN):: Centre, var
    REAL(KIND=DBL), INTENT(IN) :: rint !input real radius
    COMPLEX(KIND=DBL), INTENT(OUT):: om
    REAL(KIND=DBL) :: x,y,u,v,epsB

    u=REAL(var)*squirclecoef
    v=AIMAG(var)*squirclecoef
    epsB = 1e-6 !set boundary for "zero" to capture corners
    !WRITE(*,*) squirclecoef
    ! Define midpoints when either u or v are 0, as well as pathological -1,-1 corner
    IF ( (ABS(u) < epsB) .AND. (v > 0.) ) THEN
       x=0.;
       y=squirclecoef;
       u=ABS(u);
    ELSEIF ( (ABS(u) < epsB) .AND. (v < 0.) ) THEN
       x = 0.;
       y = -squirclecoef;
       u=ABS(u);
    ELSEIF ( (ABS(v) < epsB) .AND. (u > 0.) ) THEN
       x = squirclecoef
       y = 0.;
       v=ABS(v);
    ELSEIF ( (ABS(v) < epsB) .AND. (u < 0.) ) THEN
       x = -squirclecoef;
       y = 0.;
       v=ABS(v);
    ELSEIF ( ( ABS(u + SQRT(2.)/2.) < epsB) .AND. ( ABS(v + SQRT(2.)/2.) < epsB)) THEN
       x = squirclecoef
       y = squirclecoef
    ELSE
       !Transformation to unit square around 0
       x = ABS(u)/(u*v*SQRT(2.))*SQRT(u**2+v**2-SQRT( (u**2+v**2)*(u**2+v**2-4*u**2*v**2)))
       y = ABS(v)/(v*u*SQRT(2.))*SQRT(u**2+v**2-SQRT( (u**2+v**2)*(u**2+v**2-4*u**2*v**2)))
    ENDIF

    x=x/squirclecoef
    y=y/squirclecoef

    !Detail correct behaviour in each quadrant
    IF ( (u < 0) .AND. (v > 0) ) THEN
       y = -y
    ELSEIF ( (u < 0) .AND. (v < 0) ) THEN
       x = -x
       y = -y
    ELSEIF ( (u > 0) .AND. (v < 0) ) THEN
       x = -x
    ENDIF
    
     !Transform to rectangle around Center. Goes above imaginary and avoids real axis (both by mincont) 
    om = CMPLX((REAL(Centre)+x*rint),(AIMAG(Centre)+y*(AIMAG(Centre)-mincont)));
   
  END SUBROUTINE squircle


  SUBROUTINE unwrap(p,M,q)
    !If the difference between consecutive entries in the input array is more than cutoff (defined within)
    !then that consecutive entries is shifted to be within the same -pi + pi interval 

    REAL(KIND=DBL) :: cutoff, vmin
    INTEGER, INTENT(IN) :: M
    REAL(KIND=DBL), DIMENSION(M), INTENT(IN) :: p
    REAL(KIND=DBL), DIMENSION(M), INTENT(OUT) :: q
    REAL(KIND=DBL), DIMENSION(M) :: pmin, var, b, c, d, e, f
    REAL(KIND=DBL), DIMENSION(M-1) :: diff
    INTEGER :: i

    cutoff = 1.01*pi

    DO i=1,M-1
       diff(i) = p(i+1)-p(i)
    END DO

    b(1) = p(1)
    b(2:) = diff

    WHERE ( b > cutoff )
       c = -1.
    ELSEWHERE
       c = 0.
    END WHERE

    WHERE ( b < -cutoff )
       d = 1.
    ELSEWHERE
       d = 0.
    END WHERE

    e = (c + d) * 2*pi

    f(1) = e(1)
    DO i=2,M
       f(i) = f(i-1) + e(i) 
    END DO

    q = p+f

  END SUBROUTINE unwrap

  COMPLEX(KIND=DBL) FUNCTION casin(z)
    ! -----------------------------------------------------------------
    ! Calculates the arcsin for complex input according to the formula:
    ! Formula:  asin(z) = -i * log( i*z + ( 1 -z^2 )^0.5)
    ! -----------------------------------------------------------------

    ! Arguments
    REAL(KIND=DBL), INTENT(in) :: z

    ! Local variables
    REAL(KIND=DBL) :: aux
    COMPLEX(KIND=DBL) :: caux

    aux = (1.-(z**2.))
    IF (aux > 0.) THEN
       caux = CMPLX(aux**0.5, 0.)
    ELSE
       caux = CMPLX(0., ABS(aux)**0.5)
    ENDIF
    casin = -ci * LOG( ci*z + caux)

  END FUNCTION casin


END MODULE mod_contour
