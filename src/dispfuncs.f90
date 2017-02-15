MODULE dispfuncs
  USE kind
  USE datcal
  USE datmat, ONLY : weidcount 
  IMPLICIT NONE
 
CONTAINS

  COMPLEX(KIND=DBL) FUNCTION wofzweid(Z)
    IMPLICIT NONE
    REAL(KIND=DBL), DIMENSION(16) :: CFFT(16)
    REAL(KIND=DBL) :: L,Y
    INTEGER, PARAMETER :: N=16
    COMPLEX(KIND=DBL), INTENT(IN) :: Z
    COMPLEX(KIND=DBL) :: Zloc, i, Zout, A, B, P, RES
    !COEFFICIENTS TO USE FOR N=16
    CFFT(1) =  9.939322535568174d-07 
    CFFT(2) = 3.981287575088865d-06 
    CFFT(3) = -5.584233411098927d-06 
    CFFT(4) = -2.734640462448423d-05 
    CFFT(5) = 2.170986793223473d-05 
    CFFT(6) = 2.107105639653287d-04 
    CFFT(7) = 8.703158428458035d-05 
    CFFT(8) = -0.001527659740122_DBL 
    CFFT(9) = -0.003881015189023_DBL
    CFFT(10)=  0.003682567317092_DBL 
    CFFT(11)=  0.051822402431612_DBL 
    CFFT(12)=  0.191241726746695_DBL 
    CFFT(13)=  0.469290900903604_DBL 
    CFFT(14)=  0.886447830205055_DBL 
    CFFT(15)=  1.362240822271959_DBL 
    CFFT(16)=  1.748395886081962_DBL 

    L = 2.**(-0.25_DBL)*SQRT(N*1.0d0)
    Zloc=Z
    Y=AIMAG(Zloc)
    IF (Y<0) Zloc=-Zloc
    Zout = (L+ci*Zloc)/(L-ci*Zloc)
    A = (1/SQRT(pi))/(L-ci*Zloc)
    B = 2./((L-ci*Zloc)**2.)
    CALL cpolyev(N,Zout,CFFT,P)
    Zout = A+B*P
    IF (Y<0) Zout=-Zout
    wofzweid=Zout

    weidcount=weidcount+1

  END FUNCTION wofzweid

  SUBROUTINE cpolyev(NN,S,P,PV)
    ! EVALUATES A COMPLEX POLYNOMIAL  P  AT  S  BY THE HORNER RECURRENCE
    ! PLACING THE PARTIAL SUMS IN Q AND THE COMPUTED VALUE IN PV.
    INTEGER NN,i
    REAL(KIND=DBL) P(NN)
    COMPLEX(KIND=DBL) :: S,PV
    PV = CMPLX(P(1),0._DBL)
    DO i = 2,NN
       PV = PV*S+P(i)
    ENDDO
  END SUBROUTINE cpolyev

  COMPLEX(KIND=DBL) FUNCTION Z1(zz)
    !--------------------------------------------------------------
    ! Z1 Fried-Conte function
    ! 
    ! An asymptotic limit is carried out if the variable
    ! is larger than 50 (small fkstar). The Fadeeva function Z0
    ! is calculated with the Weideman algorithm 
    !--------------------------------------------------------------

    IMPLICIT NONE

    COMPLEX(KIND=DBL), INTENT(in)  :: zz
    COMPLEX(KIND=DBL) :: Z, som1, zz2 ,zpuiss
    REAL(KIND=DBL)    :: abz
    INTEGER :: i, ifailloc

    abz = ABS(zz)
    zz2 = zz*zz

    IF (abz<abslim) THEN
       ifailloc=0
       Z = ci * sqrtpi * wofzweid(zz)     
       Z1 = zz + zz2 * Z
    ELSE IF (abz>1.e4) THEN
       Z1 = 0.
    ELSE
       ! Asymptotic expansion for Z1 at large zz.
       ! Delicate treatment because the infinite sequence
       ! does not converge to zero, therefore difficult to
       ! choose the expansion order. Here we go up to "nerr = 20"
       ! This choice was done by trial and error.
       som1 = 0.
       zpuiss=CMPLX(1._DBL,0._DBL)
       DO i = 1,nerr
          zpuiss   = zpuiss*zz2
          som1 = som1 + pduittab(i) / zpuiss 
       END DO
       Z1 = -zz * som1
    END IF

  END FUNCTION Z1

  COMPLEX(KIND=DBL) FUNCTION Z2(zz)

    !--------------------------------------------------------------
    ! Z2 Fried-Conte function
    ! 
    ! An asymptotic limit is carried out if the variable
    ! is larger than 50 (small fkstar). The Fadeeva function Z0
    ! is calculated with the Weideman algorithm 
    !--------------------------------------------------------------
    IMPLICIT NONE

    COMPLEX(KIND=DBL) , INTENT(IN) :: zz
    COMPLEX(KIND=DBL) :: Z, som2, zz2, zpuiss
    REAL(KIND=DBL)    :: abz
    INTEGER :: i, ifailloc

    abz = ABS(zz)
    zz2 = zz*zz

    IF (abz<abslim) THEN
       ifailloc = 0
       Z = ci * sqrtpi * wofzweid(zz)     
       Z2 = zz*(0.5 + zz2 * (1.0 + zz*Z))
    ELSE IF (abz>1.e4) THEN 
       Z2 = 0.
    ELSE     
       ! Asymptotic expansion for Z2 at large zz.
       ! Delicate treatment because the infinite sequence 
       ! does not converge to zero, therefore difficult to 
       ! choose the expansion order. Here we go up to"nerr = 20" 
       ! This choice was done by trial and error. 
       som2 = 0.
       zpuiss = zz2
       DO i = 2,nerr
          zpuiss=zpuiss*zz2
          som2 = som2 + pduittab(i) / zpuiss 
       END DO
       Z2 = -zz2*zz * som2
    END IF

  END FUNCTION Z2

  COMPLEX(KIND=DBL) FUNCTION Z3(zz)
    !--------------------------------------------------------------
    ! Z3 Fried-Conte function
    ! 
    ! An asymptotic limit is carried out if the variable
    ! is larger than 50 (small fkstar). The Fadeeva function Z0
    ! is calculated with the Weideman algorithm 
    !--------------------------------------------------------------

    IMPLICIT NONE
    COMPLEX(KIND=DBL) , INTENT(IN) :: zz
    COMPLEX(KIND=DBL) :: Z, som3, zz2
    REAL(KIND=DBL)    :: abz, pduit3
    INTEGER :: i,j, ifailloc

    abz = ABS(zz)
    zz2 = zz*zz

    IF (abz<abslim) THEN
       ifailloc = 0
       Z = ci * sqrtpi * wofzweid(zz)     
       Z3 = 0.75*zz + zz2 * (0.5*zz + zz2 * (zz + zz2 * Z))
    ELSE IF (abz>1.e4) THEN 
       Z3 = 0.
    ELSE 
       ! Asymptotic expansion for Z3 at large zz.
       ! Delicate treatment because the infinite sequence 
       ! does not converge to zero, therefore difficult to 
       ! choose the expansion order. Here we go up to"nerr = 20" 
       ! This choice was done by trial and error. 
       som3 = 0.
       DO i = 3,nerr
          pduit3 = 1.
          DO j=1,i
             pduit3 = pduit3 * (2.*j-1)
          END DO
          som3 = som3 + pduit3 / (2.**i * zz2**(i)) 
       END DO
       Z3 = -zz**5. * som3
    END IF

  END FUNCTION Z3

  COMPLEX(KIND=DBL) FUNCTION Z4(zz)
    !--------------------------------------------------------------
    ! Z4 Fried-Conte function
    ! 
    ! An asymptotic limit is carried out if the variable
    ! is larger than 50 (small fkstar). The Fadeeva function Z0
    ! is calculated with the Weideman algorithm 
    !--------------------------------------------------------------

    IMPLICIT NONE
    COMPLEX(KIND=DBL) , INTENT(IN) :: zz
    COMPLEX(KIND=DBL) :: Z, som4, zz2
    REAL(KIND=DBL)    :: abz, pduit4
    INTEGER :: i,j, ifailloc

    abz = ABS(zz)
    zz2 = zz*zz

    IF (abz<abslim) THEN
       ifailloc = 0
       Z = ci * sqrtpi * wofzweid(zz)     
       !Z4 = 3.75*zz+ zz2*(0.75*zz + zz2 * (0.5*zz + zz2 * (zz + zz2 * Z)))
       Z4 = 1.875*zz+ zz2*(0.75*zz + zz2 * (0.5*zz + zz2 * (zz + zz2 * Z)))
    ELSE IF (abz>1.e4) THEN 
       Z4 = 0.
    ELSE 
       ! Asymptotic expansion for Z3 at large zz.
       ! Delicate treatment because the infinite sequence 
       ! does not converge to zero, therefore difficult to 
       ! choose the expansion order. Here we go up to"nerr = 20" 
       ! This choice was done by trial and error. 
       som4 = 0.
       DO i = 4,nerr
          pduit4 = 1.
          DO j=1,i
             pduit4 = pduit4 * (2.*j-1)
          END DO
          som4 = som4 + pduit4 / (2.**i * zz2**(i)) 
       END DO
       Z4 = -zz**7. * som4
    END IF

  END FUNCTION Z4

END MODULE dispfuncs
