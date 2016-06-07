*DECK CPQR79
      SUBROUTINE CPQR79 (NDEG, COEFF, ROOT, IERR, WORK)
C***BEGIN PROLOGUE  CPQR79
C***PURPOSE  Find the zeros of a polynomial with complex coefficients.
C***LIBRARY   SLATEC
C***CATEGORY  F1A1B
C***TYPE      COMPLEX (RPQR79-S, CPQR79-C)
C***KEYWORDS  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
C***AUTHOR  Vandevender, W. H., (SNLA)
C***DESCRIPTION
C
C   Abstract
C       This routine computes all zeros of a polynomial of degree NDEG
C       with complex coefficients by computing the eigenvalues of the
C       companion matrix.
C
C   Description of Parameters
C       The user must dimension all arrays appearing in the call list
C            COEFF(NDEG+1), ROOT(NDEG), WORK(2*NDEG*(NDEG+1))
C
C    --Input--
C      NDEG    degree of polynomial
C
C      COEFF   COMPLEX coefficients in descending order.  i.e.,
C              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
C
C      WORK    REAL work array of dimension at least 2*NDEG*(NDEG+1)
C
C   --Output--
C      ROOT    COMPLEX vector of roots
C
C      IERR    Output Error Code
C           - Normal Code
C          0  means the roots were computed.
C           - Abnormal Codes
C          1  more than 30 QR iterations on some eigenvalue of the
C             companion matrix
C          2  COEFF(1)=0.0
C          3  NDEG is invalid (less than or equal to 0)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  COMQR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   791201  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   911010  Code reworked and simplified.  (RWC and WRB)
C***END PROLOGUE  CPQR79
      COMPLEX COEFF(*), ROOT(*), SCALE, C
      REAL WORK(*)
      INTEGER NDEG, IERR, K, KHR, KHI, KWR, KWI, KAD, KJ
C***FIRST EXECUTABLE STATEMENT  CPQR79
      IERR = 0
      IF (ABS(COEFF(1)) .EQ. 0.0) THEN
         IERR = 2
         CALL XERMSG ('SLATEC', 'CPQR79',
     +      'LEADING COEFFICIENT IS ZERO.', 2, 1)
         RETURN
      ENDIF
C
      IF (NDEG .LE. 0) THEN
         IERR = 3
         CALL XERMSG ('SLATEC', 'CPQR79', 'DEGREE INVALID.', 3, 1)
         RETURN
      ENDIF
C
      IF (NDEG .EQ. 1) THEN
         ROOT(1) = -COEFF(2)/COEFF(1)
         RETURN
      ENDIF
C
      SCALE = 1.0E0/COEFF(1)
      KHR = 1
      KHI = KHR+NDEG*NDEG
      KWR = KHI+KHI-KHR
      KWI = KWR+NDEG
C
      DO 10 K=1,KWR
         WORK(K) = 0.0E0
   10 CONTINUE
C
      DO 20 K=1,NDEG
         KAD = (K-1)*NDEG+1
         C = SCALE*COEFF(K+1)
         WORK(KAD) = -REAL(C)
         KJ = KHI+KAD-1
         WORK(KJ) = -AIMAG(C)
         IF (K .NE. NDEG) WORK(KAD+K) = 1.0E0
   20 CONTINUE
C
      CALL COMQR (NDEG,NDEG,1,NDEG,WORK(KHR),WORK(KHI),WORK(KWR),
     1   WORK(KWI),IERR)
C
      IF (IERR .NE. 0) THEN
         IERR = 1
         CALL XERMSG ('SLATEC', 'CPQR79',
     +      'NO CONVERGENCE IN 30 QR ITERATIONS.', 1, 1)
         RETURN
      ENDIF
C
      DO 30 K=1,NDEG
         KM1 = K-1
         ROOT(K) = CMPLX(WORK(KWR+KM1),WORK(KWI+KM1))
   30 CONTINUE
      RETURN
      END
*DECK COMQR
      SUBROUTINE COMQR (NM, N, LOW, IGH, HR, HI, WR, WI, IERR)
C***BEGIN PROLOGUE  COMQR
C***PURPOSE  Compute the eigenvalues of complex upper Hessenberg matrix
C            using the QR method.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      COMPLEX (HQR-S, COMQR-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of a unitary analogue of the
C     ALGOL procedure  COMLR, NUM. MATH. 12, 369-376(1968) by Martin
C     and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     The unitary analogue substitutes the QR algorithm of Francis
C     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
C
C     This subroutine finds the eigenvalues of a COMPLEX
C     upper Hessenberg matrix by the QR method.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, HR and HI, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix H=(HR,HI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix, N.
C
C        HR and HI contain the real and imaginary parts, respectively,
C          of the complex upper Hessenberg matrix.  Their lower
C          triangles below the subdiagonal contain information about
C          the unitary transformations used in the reduction by  CORTH,
C          if performed.  HR and HI are two-dimensional REAL arrays,
C          dimensioned HR(NM,N) and HI(NM,N).
C
C     On OUTPUT
C
C        The upper Hessenberg portions of HR and HI have been
C          destroyed.  Therefore, they must be saved before calling
C          COMQR  if subsequent calculation of eigenvectors is to
C          be performed.
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues of the upper Hessenberg matrix.  If an
C          error exit is made, the eigenvalues should be correct for
C          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
C          dimensional REAL arrays, dimensioned WR(N) and WI(N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30*N iterations.
C                     The eigenvalues should be correct for indices
C                     IERR+1, IERR+2, ..., N.
C
C     Calls CSROOT for complex square root.
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C     Calls CDIV for complex division.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  CDIV, CSROOT, PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  COMQR
C
      INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
      REAL HR(NM,*),HI(NM,*),WR(*),WI(*)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  COMQR
      IERR = 0
      IF (LOW .EQ. IGH) GO TO 180
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
      L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0E0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0E0
C
         DO 155 J = I, IGH
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = LOW, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW E0 -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         S1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     1             + ABS(HR(L,L)) +ABS(HI(L,L))
         S2 = S1 + ABS(HR(L,L-1))
         IF (S2 .EQ. S1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = 0.0E0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0E0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0E0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, EN
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0E0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0E0
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = L, J
            YR = HR(I,J-1)
            YI = 0.0E0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0E0) GO TO 240
C
      DO 630 I = L, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GO TO 220
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
*DECK CSROOT
      SUBROUTINE CSROOT (XR, XI, YR, YI)
C***BEGIN PROLOGUE  CSROOT
C***SUBSIDIARY
C***PURPOSE  Compute the complex square root of a complex number.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (CSROOT-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     (YR,YI) = complex sqrt(XR,XI)
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CSROOT
      REAL XR,XI,YR,YI,S,TR,TI,PYTHAG
C
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C***FIRST EXECUTABLE STATEMENT  CSROOT
      TR = XR
      TI = XI
      S = SQRT(0.5E0*(PYTHAG(TR,TI) + ABS(TR)))
      IF (TR .GE. 0.0E0) YR = S
      IF (TI .LT. 0.0E0) S = -S
      IF (TR .LE. 0.0E0) YI = S
      IF (TR .LT. 0.0E0) YR = 0.5E0*(TI/YI)
      IF (TR .GT. 0.0E0) YI = 0.5E0*(TI/YR)
      RETURN
      END
*DECK PYTHAG
      REAL FUNCTION PYTHAG (A, B)
C***BEGIN PROLOGUE  PYTHAG
C***SUBSIDIARY
C***PURPOSE  Compute the complex square root of a complex number without
C            destructive overflow or underflow.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (PYTHAG-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Finds sqrt(A**2+B**2) without overflow or destructive underflow
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  PYTHAG
      REAL A,B
C
      REAL P,Q,R,S,T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = MAX(ABS(A),ABS(B))
      Q = MIN(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
*DECK CDIV
      SUBROUTINE CDIV (AR, AI, BR, BI, CR, CI)
C***BEGIN PROLOGUE  CDIV
C***SUBSIDIARY
C***PURPOSE  Compute the complex quotient of two complex numbers.
C***LIBRARY   SLATEC
C***TYPE      COMPLEX (CDIV-C)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  CDIV
      REAL AR,AI,BR,BI,CR,CI
C
      REAL S,ARS,AIS,BRS,BIS
C***FIRST EXECUTABLE STATEMENT  CDIV
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
