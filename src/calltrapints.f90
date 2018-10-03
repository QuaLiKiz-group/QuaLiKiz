MODULE calltrapints
  USE kind
  USE trapints
  USE datmat, ONLY : nions, ion_type, ion, nuFFk
  IMPLICIT NONE

CONTAINS


  INTEGER FUNCTION trapped_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, adiabatic, scale_2
    INTEGER :: i
    
    adiabatic = fdata(1)
    scale_2 = fdata(2) !scale_2 is vuplim 
    scale_ = ABS(adiabatic) !scaling the integrand
    
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      trapped_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1)
    XY(2) = x(2)
    xx = x(1) 
    
    output = FFke(ndim, XY, 1)
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) / scale_2
    END DO    
        
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = AIMAG(output)
    trapped_cubature = 0
    
    
  END FUNCTION trapped_cubature

  FUNCTION FFke_cub(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFke_cub
    FFke_cub(1) = REAL ( FFke(nf, kv, 1) )
    FFke_cub(2) = AIMAG ( FFke(nf, kv, 1) )
  END FUNCTION FFke_cub

  FUNCTION FFkgte_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ate term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkgte_cub
    FFkgte_cub(1) = REAL ( FFke(nf, kv, 2) )
    FFkgte_cub(2) = AIMAG ( FFke(nf, kv, 2) )
  END FUNCTION FFkgte_cub

  FUNCTION FFkgne_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkgne_cub
    FFkgne_cub(1) = REAL ( FFke(nf, kv, 3) )
    FFkgne_cub(2) = AIMAG ( FFke(nf, kv, 3) )
  END FUNCTION FFkgne_cub

  FUNCTION FFkce_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkce_cub
    FFkce_cub(1) = REAL ( FFke(nf, kv, 4) )
    FFkce_cub(2) = AIMAG ( FFke(nf, kv, 4) )
  END FUNCTION FFkce_cub

  FUNCTION FFeke_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFeke_cub
    FFeke_cub(1) = REAL ( FFke(nf, kv, 5) )
    FFeke_cub(2) = AIMAG ( FFke(nf, kv, 5) )
  END FUNCTION FFeke_cub

  FUNCTION FFekgte_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ate term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekgte_cub
    FFekgte_cub(1) = REAL ( FFke(nf, kv, 6) )
    FFekgte_cub(2) = AIMAG ( FFke(nf, kv, 6) )
  END FUNCTION FFekgte_cub

  FUNCTION FFekgne_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekgne_cub
    FFekgne_cub(1) = REAL ( FFke(nf, kv, 7) )
    FFekgne_cub(2) = AIMAG ( FFke(nf, kv, 7) )
  END FUNCTION FFekgne_cub

  FUNCTION FFekce_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekce_cub
    FFekce_cub(1) = REAL ( FFke(nf, kv, 8) )
    FFekce_cub(2) = AIMAG ( FFke(nf, kv, 8) )
  END FUNCTION FFekce_cub

  REAL(KIND=DBL) FUNCTION rFFkiz(kk)
    !-------------------------------------------------------------
    ! Returns the real component of the trapped ion (all ions) integrand
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
       IF ( (ion_type(pFFk,i) == 1) .AND. (ETG_flag(nuFFk) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum+REAL ( FFki(kk,1,i) )*ninorm(pFFk,i) !unnormalise coefi here
       ENDIF
    ENDDO
    rFFkiz=intsum
  END FUNCTION rFFkiz

  REAL(KIND=DBL) FUNCTION iFFkiz(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of the trapped ion (all ions) integrand
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
       IF ( (ion_type(pFFk,i) == 1) .AND. (ETG_flag(nuFFk) .EQV. .FALSE.) ) THEN !only include active ions in dispersion relation nonadiabatic part
          intsum = intsum+AIMAG ( FFki(kk,1,i) )*ninorm(pFFk,i) !unnormalise coefi here
       ENDIF
    ENDDO
    iFFkiz = intsum
  END FUNCTION iFFkiz

  REAL(KIND=DBL) FUNCTION rFFki(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki (from only given ion), full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFki = REAL ( FFki(kk,1,ion) )
  END FUNCTION rFFki

  REAL(KIND=DBL) FUNCTION iFFki(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFki = AIMAG ( FFki(kk,1,ion) )
  END FUNCTION iFFki

  REAL(KIND=DBL) FUNCTION rFFkgti(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki (from only given ion), with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgti = REAL ( FFki(kk,2,ion) )
  END FUNCTION rFFkgti

  REAL(KIND=DBL) FUNCTION iFFkgti(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgti = AIMAG ( FFki(kk,2,ion) )
  END FUNCTION iFFkgti

  REAL(KIND=DBL) FUNCTION rFFkgni(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgni = REAL ( FFki(kk,3,ion) )
  END FUNCTION rFFkgni

  REAL(KIND=DBL) FUNCTION iFFkgni(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgni = AIMAG ( FFki(kk,3,ion) )
  END FUNCTION iFFkgni

  REAL(KIND=DBL) FUNCTION rFFkci(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkci = REAL ( FFki(kk,4,ion) )
  END FUNCTION rFFkci

  REAL(KIND=DBL) FUNCTION iFFkci(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkci = AIMAG ( FFki(kk,4,ion) )
  END FUNCTION iFFkci

  REAL(KIND=DBL) FUNCTION rFFeki(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFeki = REAL ( FFki(kk,5,ion) )
  END FUNCTION rFFeki

  REAL(KIND=DBL) FUNCTION iFFeki(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFeki = AIMAG ( FFki(kk,5,ion) )
  END FUNCTION iFFeki

  REAL(KIND=DBL) FUNCTION rFFekgti(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki (from only given ion), with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgti = REAL ( FFki(kk,6,ion) )
  END FUNCTION rFFekgti

  REAL(KIND=DBL) FUNCTION iFFekgti(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgti = AIMAG ( FFki(kk,6,ion) )
  END FUNCTION iFFekgti

  REAL(KIND=DBL) FUNCTION rFFekgni(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgni = REAL ( FFki(kk,7,ion) )
  END FUNCTION rFFekgni

  REAL(KIND=DBL) FUNCTION iFFekgni(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgni = AIMAG ( FFki(kk,7,ion) )
  END FUNCTION iFFekgni

  REAL(KIND=DBL) FUNCTION rFFekci(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekci = REAL ( FFki(kk,8,ion) )
  END FUNCTION rFFekci

  REAL(KIND=DBL) FUNCTION iFFekci(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekci = AIMAG ( FFki(kk,8,ion) )
  END FUNCTION iFFekci


  REAL(KIND=DBL) FUNCTION rFFke_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFke_nocoll, full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFke_nocoll = REAL ( FFke_nocoll(kk,1) )
  END FUNCTION rFFke_nocoll

  REAL(KIND=DBL) FUNCTION iFFke_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFke_nocoll, full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFke_nocoll = AIMAG ( FFke_nocoll(kk,1) )
  END FUNCTION iFFke_nocoll

  REAL(KIND=DBL) FUNCTION rFFkgte_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFke_nocoll, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgte_nocoll = REAL ( FFke_nocoll(kk,2) )
  END FUNCTION rFFkgte_nocoll

  REAL(KIND=DBL) FUNCTION iFFkgte_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFke_nocoll, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgte_nocoll = AIMAG ( FFke_nocoll(kk,2) )
  END FUNCTION iFFkgte_nocoll

  REAL(KIND=DBL) FUNCTION rFFkgne_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFke_nocoll, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgne_nocoll = REAL ( FFke_nocoll(kk,3) )
  END FUNCTION rFFkgne_nocoll

  REAL(KIND=DBL) FUNCTION iFFkgne_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFke_nocoll, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgne_nocoll = AIMAG ( FFke_nocoll(kk,3) )
  END FUNCTION iFFkgne_nocoll

  REAL(KIND=DBL) FUNCTION rFFkce_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFke_nocoll, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkce_nocoll = REAL ( FFke_nocoll(kk,4) )
  END FUNCTION rFFkce_nocoll

  REAL(KIND=DBL) FUNCTION iFFkce_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFke_nocoll, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkce_nocoll = AIMAG ( FFke_nocoll(kk,4) )
  END FUNCTION iFFkce_nocoll

  REAL(KIND=DBL) FUNCTION rFFeke_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFke_nocoll, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFeke_nocoll = REAL ( FFke_nocoll(kk,5) )
  END FUNCTION rFFeke_nocoll

  REAL(KIND=DBL) FUNCTION iFFeke_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFke_nocoll, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFeke_nocoll = AIMAG ( FFke_nocoll(kk,5) )

  END FUNCTION iFFeke_nocoll

  REAL(KIND=DBL) FUNCTION rFFekgte_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeke_nocoll, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgte_nocoll = REAL ( FFke_nocoll(kk,6) )
  END FUNCTION rFFekgte_nocoll

  REAL(KIND=DBL) FUNCTION iFFekgte_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeke_nocoll, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgte_nocoll = AIMAG ( FFke_nocoll(kk,6) )
  END FUNCTION iFFekgte_nocoll

  REAL(KIND=DBL) FUNCTION rFFekgne_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeke_nocoll, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgne_nocoll = REAL ( FFke_nocoll(kk,7) )
  END FUNCTION rFFekgne_nocoll

  REAL(KIND=DBL) FUNCTION iFFekgne_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeke_nocoll, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgne_nocoll = AIMAG ( FFke_nocoll(kk,7) )
  END FUNCTION iFFekgne_nocoll

  REAL(KIND=DBL) FUNCTION rFFekce_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeke_nocoll, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekce_nocoll = REAL ( FFke_nocoll(kk,8) )
  END FUNCTION rFFekce_nocoll

  REAL(KIND=DBL) FUNCTION iFFekce_nocoll(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeke_nocoll, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekce_nocoll = AIMAG ( FFke_nocoll(kk,8) )
  END FUNCTION iFFekce_nocoll



  REAL(KIND=DBL) FUNCTION rFFke(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFke = REAL ( FFke(nf, kv, 1) )
  END FUNCTION rFFke

  REAL(KIND=DBL) FUNCTION iFFke(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFke = AIMAG ( FFke(nf, kv, 1) )
  END FUNCTION iFFke

  REAL(KIND=DBL) FUNCTION rFFkgte(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkgte = REAL ( FFke(nf, kv, 2) )
  END FUNCTION rFFkgte

  REAL(KIND=DBL) FUNCTION iFFkgte(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkgte = AIMAG ( FFke(nf, kv, 2) )
  END FUNCTION iFFkgte

  REAL(KIND=DBL) FUNCTION rFFkgne(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkgne = REAL ( FFke(nf, kv, 3) )
  END FUNCTION rFFkgne

  REAL(KIND=DBL) FUNCTION iFFkgne(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkgne = AIMAG ( FFke(nf, kv, 3) )
  END FUNCTION iFFkgne

  REAL(KIND=DBL) FUNCTION rFFkce(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkce = REAL ( FFke(nf, kv, 4) )
  END FUNCTION rFFkce

  REAL(KIND=DBL) FUNCTION iFFkce(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkce = AIMAG ( FFke(nf, kv, 4) )
  END FUNCTION iFFkce

  REAL(KIND=DBL) FUNCTION rFFeke(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron energy integrand 
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFeke = REAL ( FFke(nf, kv, 5) )
  END FUNCTION rFFeke

  REAL(KIND=DBL) FUNCTION iFFeke(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFeke = AIMAG ( FFke(nf, kv, 5) )
  END FUNCTION iFFeke

  REAL(KIND=DBL) FUNCTION rFFekgte(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekgte = REAL ( FFke(nf, kv, 6) )
  END FUNCTION rFFekgte

  REAL(KIND=DBL) FUNCTION iFFekgte(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekgte = AIMAG ( FFke(nf, kv, 6) )
  END FUNCTION iFFekgte

  REAL(KIND=DBL) FUNCTION rFFekgne(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekgne = REAL ( FFke(nf, kv, 7) )
  END FUNCTION rFFekgne

  REAL(KIND=DBL) FUNCTION iFFekgne(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekgne = AIMAG ( FFke(nf, kv, 7) )
  END FUNCTION iFFekgne

  REAL(KIND=DBL) FUNCTION rFFekce(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekce = REAL ( FFke(nf, kv, 8) )
  END FUNCTION rFFekce

  REAL(KIND=DBL) FUNCTION iFFekce(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekce = AIMAG ( FFke(nf, kv, 8) )
  END FUNCTION iFFekce


  !******************************************************************************************************
  ! with rotation flag on
  !*******************************************************************************************************

  FUNCTION FFkerot_cub(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkerot_cub
    FFkerot_cub(1) = REAL ( FFkerot(nf, kv, 1) )
    FFkerot_cub(2) = AIMAG ( FFkerot(nf, kv, 1) )
  END FUNCTION FFkerot_cub

  FUNCTION FFkgterot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ate term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkgterot_cub
    FFkgterot_cub(1) = REAL ( FFkerot(nf, kv, 2) )
    FFkgterot_cub(2) = AIMAG ( FFkerot(nf, kv, 2) )
  END FUNCTION FFkgterot_cub

  FUNCTION FFkgnerot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkgnerot_cub
    FFkgnerot_cub(1) = REAL ( FFkerot(nf, kv, 3) )
    FFkgnerot_cub(2) = AIMAG ( FFkerot(nf, kv, 3) )
  END FUNCTION FFkgnerot_cub

  FUNCTION FFkcerot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFkcerot_cub
    FFkcerot_cub(1) = REAL ( FFkerot(nf, kv, 4) )
    FFkcerot_cub(2) = AIMAG ( FFkerot(nf, kv, 4) )
  END FUNCTION FFkcerot_cub

  FUNCTION FFekerot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekerot_cub
    FFekerot_cub(1) = REAL ( FFkerot(nf, kv, 5) )
    FFekerot_cub(2) = AIMAG ( FFkerot(nf, kv, 5) )
  END FUNCTION FFekerot_cub

  FUNCTION FFekgterot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ate term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekgterot_cub
    FFekgterot_cub(1) = REAL ( FFkerot(nf, kv, 6) )
    FFekgterot_cub(2) = AIMAG ( FFkerot(nf, kv, 6) )
  END FUNCTION FFekgterot_cub

  FUNCTION FFekgnerot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term only)
    !---------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekgnerot_cub
    FFekgnerot_cub(1) = REAL ( FFkerot(nf, kv, 7) )
    FFekgnerot_cub(2) = AIMAG ( FFkerot(nf, kv, 7) )
  END FUNCTION FFekgnerot_cub

  FUNCTION FFekcerot_cub(nf, kv)
    !---------------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !--------------------------------------------------------------------------- 
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    REAL(KIND=DBL), DIMENSION(nf) :: FFekcerot_cub
    FFekcerot_cub(1) = REAL ( FFkerot(nf, kv, 8) )
    FFekcerot_cub(2) = AIMAG ( FFkerot(nf, kv, 8) )
  END FUNCTION FFekcerot_cub

  REAL(KIND=DBL) FUNCTION rFFkizrot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of the trapped ion (all ions) integrand
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
       IF ( (ion_type(pFFk,i) == 1) .AND. (ETG_flag(nuFFk) .EQV. .FALSE.) ) THEN !only include active ions
          intsum = intsum+REAL ( FFkirot(kk,1,i) )*ninorm(pFFk,i) !unnormalise coefi here
       ENDIF
    ENDDO
    rFFkizrot=intsum
  END FUNCTION rFFkizrot

  REAL(KIND=DBL) FUNCTION iFFkizrot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of the trapped ion (all ions) integrand
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
       IF ( (ion_type(pFFk,i) == 1) .AND. (ETG_flag(nuFFk) .EQV. .FALSE.) ) THEN !only include active ions in dispersion relation nonadiabatic part
          intsum = intsum+AIMAG ( FFkirot(kk,1,i) )*ninorm(pFFk,i) !unnormalise coefi here
       ENDIF
    ENDDO
    iFFkizrot = intsum
  END FUNCTION iFFkizrot

  REAL(KIND=DBL) FUNCTION rFFkirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFkirot (from only given ion), full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkirot = REAL ( FFkirot(kk,1,ion) )

!!$    WRITE(*,*) 'New rFFki',rFFkirot
!!$    rFFkirot = REAL ( FFki(kk,1,ion) )
!!$    WRITE(*,*) 'Old rFFki',rFFkirot


  END FUNCTION rFFkirot

  REAL(KIND=DBL) FUNCTION iFFkirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFkirot, full form
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkirot = AIMAG ( FFkirot(kk,1,ion) )

!!$    WRITE(*,*) 'New iFFki',iFFkirot
!!$    iFFkirot = AIMAG ( FFki(kk,1,ion) )
!!$    WRITE(*,*) 'Old iFFki',iFFkirot

  END FUNCTION iFFkirot

  REAL(KIND=DBL) FUNCTION rFFkgtirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFkirot (from only given ion), with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgtirot = REAL ( FFkirot(kk,2,ion) )
  END FUNCTION rFFkgtirot

  REAL(KIND=DBL) FUNCTION iFFkgtirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFkirot, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgtirot = AIMAG ( FFkirot(kk,2,ion) )
  END FUNCTION iFFkgtirot

  REAL(KIND=DBL) FUNCTION rFFkgnirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFkirot, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkgnirot = REAL ( FFkirot(kk,3,ion) )
  END FUNCTION rFFkgnirot

  REAL(KIND=DBL) FUNCTION iFFkgnirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFkirot, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkgnirot = AIMAG ( FFkirot(kk,3,ion) )
  END FUNCTION iFFkgnirot

  REAL(KIND=DBL) FUNCTION rFFkcirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFkirot, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkcirot = REAL ( FFkirot(kk,4,ion) )
  END FUNCTION rFFkcirot

  REAL(KIND=DBL) FUNCTION iFFkcirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFkirot, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkcirot = AIMAG ( FFkirot(kk,4,ion) )
  END FUNCTION iFFkcirot

  REAL(KIND=DBL) FUNCTION rFFekirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFkirot, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekirot = REAL ( FFkirot(kk,5,ion) )
  END FUNCTION rFFekirot

  REAL(KIND=DBL) FUNCTION iFFekirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFkirot, with added v in integrals
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekirot = AIMAG ( FFkirot(kk,5,ion) )
  END FUNCTION iFFekirot

  REAL(KIND=DBL) FUNCTION rFFekgtirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki (from only given ion), with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgtirot = REAL ( FFkirot(kk,6,ion) )
  END FUNCTION rFFekgtirot

  REAL(KIND=DBL) FUNCTION iFFekgtirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with Ati prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgtirot = AIMAG ( FFkirot(kk,6,ion) )
  END FUNCTION iFFekgtirot

  REAL(KIND=DBL) FUNCTION rFFekgnirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekgnirot = REAL ( FFkirot(kk,7,ion) )
  END FUNCTION rFFekgnirot

  REAL(KIND=DBL) FUNCTION iFFekgnirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with Ani prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekgnirot = AIMAG ( FFkirot(kk,7,ion) )
  END FUNCTION iFFekgnirot

  REAL(KIND=DBL) FUNCTION rFFekcirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekcirot = REAL ( FFkirot(kk,8,ion) )
  END FUNCTION rFFekcirot

  REAL(KIND=DBL) FUNCTION iFFekcirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with curvature term only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekcirot = AIMAG ( FFkirot(kk,8,ion) )
  END FUNCTION iFFekcirot


  REAL(KIND=DBL) FUNCTION rFFkguirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFki (from only given ion), with Aui prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFkguirot = REAL ( FFkirot(kk,9,ion) )
  END FUNCTION rFFkguirot

  REAL(KIND=DBL) FUNCTION iFFkguirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFki, with Aui prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFkguirot = AIMAG ( FFkirot(kk,9,ion) )
  END FUNCTION iFFkguirot

  REAL(KIND=DBL) FUNCTION rFFekguirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFeki, with Aui prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFekguirot = REAL ( FFkirot(kk,10,ion) )
  END FUNCTION rFFekguirot

  REAL(KIND=DBL) FUNCTION iFFekguirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFeki, with Aui prefactor only
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFekguirot = AIMAG ( FFkirot(kk,10,ion) )
  END FUNCTION iFFekguirot

  REAL(KIND=DBL) FUNCTION rFFvkirot(kk)
    !-------------------------------------------------------------
    ! Returns the real component of FFvki, ang momentum 
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    rFFvkirot = REAL ( FFkirot(kk,11,ion) )
  END FUNCTION rFFvkirot

  REAL(KIND=DBL) FUNCTION iFFvkirot(kk)
    !-------------------------------------------------------------
    ! Returns the imaginary component of FFvki, ang momentum
    !-------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER :: caseflag
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    iFFvkirot = AIMAG ( FFkirot(kk,11,ion) )
  END FUNCTION iFFvkirot


  REAL(KIND=DBL) FUNCTION rFFkerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkerot = REAL ( FFkerot(nf, kv, 1) )
!!$    WRITE(*,*) 'New rFFke',rFFkerot
!!$    rFFkerot = REAL ( FFke(nf, kv, 1) )
!!$    WRITE(*,*) 'Old rFFke',rFFkerot


  END FUNCTION rFFkerot

  REAL(KIND=DBL) FUNCTION iFFkerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkerot = AIMAG ( FFkerot(nf, kv, 1) )
!!$    WRITE(*,*) 'New iFFke',iFFkerot
!!$    iFFkerot = AIMAG ( FFke(nf, kv, 1) )
!!$    WRITE(*,*) 'Old iFFke',iFFkerot


  END FUNCTION iFFkerot

  REAL(KIND=DBL) FUNCTION rFFkgterot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkgterot = REAL ( FFkerot(nf, kv, 2) )
  END FUNCTION rFFkgterot

  REAL(KIND=DBL) FUNCTION iFFkgterot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkgterot = AIMAG ( FFkerot(nf, kv, 2) )
  END FUNCTION iFFkgterot

  REAL(KIND=DBL) FUNCTION rFFkgnerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkgnerot = REAL ( FFkerot(nf, kv, 3) )
  END FUNCTION rFFkgnerot

  REAL(KIND=DBL) FUNCTION iFFkgnerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkgnerot = AIMAG ( FFkerot(nf, kv, 3) )
  END FUNCTION iFFkgnerot

  REAL(KIND=DBL) FUNCTION rFFkcerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFkcerot = REAL ( FFkerot(nf, kv, 4) )
  END FUNCTION rFFkcerot

  REAL(KIND=DBL) FUNCTION iFFkcerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFkcerot = AIMAG ( FFkerot(nf, kv, 4) )
  END FUNCTION iFFkcerot

  REAL(KIND=DBL) FUNCTION rFFekerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron energy integrand 
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekerot = REAL ( FFkerot(nf, kv, 5) )
  END FUNCTION rFFekerot

  REAL(KIND=DBL) FUNCTION iFFekerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron energy integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekerot = AIMAG ( FFkerot(nf, kv, 5) )
  END FUNCTION iFFekerot

  REAL(KIND=DBL) FUNCTION rFFekgterot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekgterot = REAL ( FFkerot(nf, kv, 6) )
  END FUNCTION rFFekgterot

  REAL(KIND=DBL) FUNCTION iFFekgterot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekgterot = AIMAG ( FFkerot(nf, kv, 6) )
  END FUNCTION iFFekgterot

  REAL(KIND=DBL) FUNCTION rFFekgnerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekgnerot = REAL ( FFkerot(nf, kv, 7) )
  END FUNCTION rFFekgnerot

  REAL(KIND=DBL) FUNCTION iFFekgnerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekgnerot = AIMAG ( FFkerot(nf, kv, 7) )
  END FUNCTION iFFekgnerot

  REAL(KIND=DBL) FUNCTION rFFekcerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    rFFekcerot = REAL ( FFkerot(nf, kv, 8) )
  END FUNCTION rFFekcerot

  REAL(KIND=DBL) FUNCTION iFFekcerot(nf, kv)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: kv
    iFFekcerot = AIMAG ( FFkerot(nf, kv, 8) )
  END FUNCTION iFFekcerot

  REAL(KIND=DBL) FUNCTION rFFke_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand
    !---------------------------------------------------------------------  
     REAL(KIND=DBL), INTENT(IN) :: kk
    rFFke_nocollrot = REAL ( FFke_nocollrot(kk, 1) )
  END FUNCTION rFFke_nocollrot

  REAL(KIND=DBL) FUNCTION iFFke_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    iFFke_nocollrot = AIMAG ( FFke_nocollrot(kk, 1) )
  END FUNCTION iFFke_nocollrot

  REAL(KIND=DBL) FUNCTION rFFkgte_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
     REAL(KIND=DBL), INTENT(IN) :: kk
    rFFkgte_nocollrot = REAL ( FFke_nocollrot(kk, 2) )
  END FUNCTION rFFkgte_nocollrot

  REAL(KIND=DBL) FUNCTION iFFkgte_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    iFFkgte_nocollrot = AIMAG ( FFke_nocollrot(kk, 2) )
  END FUNCTION iFFkgte_nocollrot

  REAL(KIND=DBL) FUNCTION rFFkgne_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    rFFkgne_nocollrot = REAL ( FFke_nocollrot(kk, 3) )
  END FUNCTION rFFkgne_nocollrot

  REAL(KIND=DBL) FUNCTION iFFkgne_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    iFFkgne_nocollrot = AIMAG ( FFke_nocollrot(kk, 3) )
  END FUNCTION iFFkgne_nocollrot

  REAL(KIND=DBL) FUNCTION rFFkce_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    rFFkce_nocollrot = REAL ( FFke_nocollrot(kk, 4) )
  END FUNCTION rFFkce_nocollrot

  REAL(KIND=DBL) FUNCTION iFFkce_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    iFFkce_nocollrot = AIMAG ( FFke_nocollrot(kk, 4) )
  END FUNCTION iFFkce_nocollrot

  REAL(KIND=DBL) FUNCTION rFFeke_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron energy integrand 
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    rFFeke_nocollrot = REAL ( FFke_nocollrot(kk, 5) )
  END FUNCTION rFFeke_nocollrot

  REAL(KIND=DBL) FUNCTION iFFeke_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron energy integrand
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    iFFeke_nocollrot = AIMAG ( FFke_nocollrot(kk, 5) )
  END FUNCTION iFFeke_nocollrot

  REAL(KIND=DBL) FUNCTION rFFekgte_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    rFFekgte_nocollrot = REAL ( FFke_nocollrot(kk, 6) )
  END FUNCTION rFFekgte_nocollrot

  REAL(KIND=DBL) FUNCTION iFFekgte_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ati term)
    !---------------------------------------------------------------------  
     REAL(KIND=DBL), INTENT(IN) :: kk
    iFFekgte_nocollrot = AIMAG ( FFke_nocollrot(kk, 6) )
  END FUNCTION iFFekgte_nocollrot

  REAL(KIND=DBL) FUNCTION rFFekgne_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
    REAL(KIND=DBL), INTENT(IN) :: kk
    rFFekgne_nocollrot = REAL ( FFke_nocollrot(kk, 7) )
  END FUNCTION rFFekgne_nocollrot

  REAL(KIND=DBL) FUNCTION iFFekgne_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (Ane term)
    !---------------------------------------------------------------------  
      REAL(KIND=DBL), INTENT(IN) :: kk
    iFFekgne_nocollrot = AIMAG ( FFke_nocollrot(kk, 7) )
  END FUNCTION iFFekgne_nocollrot

  REAL(KIND=DBL) FUNCTION rFFekce_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the real part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
      REAL(KIND=DBL), INTENT(IN) :: kk
    rFFekce_nocollrot = REAL ( FFke_nocollrot(kk, 8) )
  END FUNCTION rFFekce_nocollrot

  REAL(KIND=DBL) FUNCTION iFFekce_nocollrot(kk)
    !---------------------------------------------------------------------
    ! Calculates the imaginary part of the trapped electron integrand (curvature term)
    !---------------------------------------------------------------------  
      REAL(KIND=DBL), INTENT(IN) :: kk
    iFFekce_nocollrot = AIMAG ( FFke_nocollrot(kk, 8) )
  END FUNCTION iFFekce_nocollrot

    
END MODULE calltrapints
