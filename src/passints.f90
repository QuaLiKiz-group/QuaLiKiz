MODULE passints
  USE kind
  USE datcal
  USE datmat
  USE dispfuncs

  IMPLICIT NONE

CONTAINS

!*************************************************************************************
! Fkstarrstari and Fkstarrstare without rotation
! Fkstarrstarirot and Fkstarrstarerot with rotation
!*************************************************************************************

  COMPLEX(KIND=DBL) FUNCTION Fkstarrstari(ndim, xx, caseflag,nion)
    !---------------------------------------------------------------------
    ! Calculates the passing ion k*, r* integrands
    ! The returned integrand depends on the calculation flag, used
    ! in case we only want certain coefficient output
    ! caseflag = 1, full output
    ! caseflag = 2, return factor in front of Ati (kgti)
    ! caseflag = 3, return factor in front of Ani (kgni)    
    ! caseflag = 4, return curvature term (kgci)   
    ! caseflag = 5, return full output for energy integral (higher v exponent)
    ! caseflag = 6, return Ati coefficient for energy integral
    ! caseflag = 7, return Ani coefficient for energy integral
    ! caseflag = 8, return curvature term  for energy integral
    !---------------------------------------------------------------------  
    ! Arguments
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL)   , INTENT(IN) :: xx(ndim)
    INTEGER , INTENT(IN) :: caseflag, nion

    REAL(KIND=DBL)    :: Athir 
    COMPLEX(KIND=DBL) :: aai, bbi, cci, ddi, sqrtdi, Vmi, Vpi, Zai
    COMPLEX(KIND=DBL) :: alphai
    COMPLEX(KIND=DBL) :: faci
    REAL(KIND=DBL)    :: nwgi
    REAL(KIND=DBL)    :: rstar, kstar, teta, fkstar
    REAL(KIND=DBL)    :: var2,var3,bessm2 !new terms for Bessel directly inside passints
    COMPLEX(KIND=DBL) :: inti3, inti5
    COMPLEX(KIND=DBL) :: Fikstarrstar

    kstar = xx(1)
    rstar = xx(2)


  
    teta = kstar*d/REAL(mwidth) / SQRT(2._DBL) !warning not sure if should be real(mwidth) or rather keep teta complex?
    !Vertical drift term
    fkstar = 4./3.*(COS(teta) + (smag(pFkr) * teta - alphax(pFkr) * SIN(teta))*SIN(teta))
    IF (ABS(fkstar)<minfki) fkstar=SIGN(minfki,fkstar)

    nwgi = nwg*(-Tix(pFkr,nion)/Zi(pFkr,nion))

    var2 = (teta/d*Rhoi(pFkr,nion))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoi(nion))**2.               !!2nd argument of Bessel fun
    bessm2 = BESEI0(var2+var3)

    !Transit frequency        
    Athir = Athi(nion)*rstar / SQRT(2._DBL)
    !Simplified calculation for zero vertical drift
    IF (ABS(fkstar)<minfki) THEN 
       !Further simplication if transit freq is zero
       IF (rstar<epsD) THEN 
          inti3 = -1./omFkr*(-Tix(pFkr,nion)/Zi(pFkr,nion))
          inti5 = 1.5*inti3
       ELSE   
          aai = omFkr*nwg/Athir
          !Fried-Conte functions
          Zai=Z1(aai)
          inti3 = nwgi / Athir *2. *Zai
          inti5 = nwgi / Athir *aai + inti3*aai*aai
       END IF
       !GENERAL CASE
    ELSE  
       ! Analytical form of velocity integration written with 
       ! Fried-Conte functions, generalizations of the simple plasma dispersion function
       bbi = CMPLX(Athir/(nwgi*fkstar),0.) 
       cci = omFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))/fkstar

       ddi = bbi**2 - 4.*cci
       sqrtdi = SQRT(ddi)

       Vmi = (-bbi-sqrtdi)/2.
       Vpi = (-bbi+sqrtdi)/2.

       faci = 2. / (fkstar * (Vpi-Vmi))

       IF (caseflag < 5) THEN !differentiate between particle or energy integrals
          inti3 = faci * (Z1(Vpi) - Z1(Vmi))
          inti5 = faci * (Z2(Vpi) - Z2(Vmi))
       ELSE
          inti3 = faci * (Z2(Vpi) - Z2(Vmi))
          inti5 = faci * (Z3(Vpi) - Z3(Vmi))
       ENDIF

    END IF

    IF ( (caseflag == 1) .OR. (caseflag == 5) ) THEN !full term for particle or energy
       alphai = omFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))+Ani(pFkr,nion)-1.5*Ati(pFkr,nion)
       Fikstarrstar = inti3 * alphai + inti5 * Ati(pFkr,nion)
    ELSEIF ( (caseflag == 2) .OR. (caseflag == 6) ) THEN !At factor only particle transport
       alphai = -1.5
       Fikstarrstar = inti3 * alphai + inti5 
    ELSEIF ( (caseflag == 3) .OR. (caseflag == 7) ) THEN !An factor only particle transport
       alphai = 1.
       Fikstarrstar = inti3 * alphai 
    ELSEIF ( (caseflag == 4) .OR. (caseflag == 8) ) THEN !Curvature term particle transport
       alphai = omFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))
       Fikstarrstar = inti3 * alphai 
    ENDIF

    !Care with norm fact (1 here, 4 in Fkstarrstar...)
    IF (caseflag < 5) THEN !differentiate between particle or energy integrals
!!$       Fkstarrstari = 1. * fc(pFkr) * coefi(pFkr,nion) &
!!$            &        * (Fikstarrstar*Joi2c(nion)) * EXP( -(kstar**2 + rstar**2)/2 )/twopi  
       Fkstarrstari = 1. * fc(pFkr) * coefi(pFkr,nion) &
            &        * (Fikstarrstar*bessm2) * EXP( -(kstar**2 + rstar**2)/2 )/twopi  

    ELSE
!!$       Fkstarrstari = 1. * fc(pFkr) * coefi(pFkr,nion) * Tix(pFkr,nion) &
!!$            &      * (Fikstarrstar*Joi2c(nion)) * EXP(- (kstar**2 + rstar**2)/2 )/twopi  
       Fkstarrstari = 1. * fc(pFkr) * coefi(pFkr,nion) * Tix(pFkr,nion) &
            &      * (Fikstarrstar*bessm2) * EXP(- (kstar**2 + rstar**2)/2 )/twopi  

    ENDIF

    IF (ABS(Fkstarrstari) < SQRT(epsD)) Fkstarrstari=0.

  END FUNCTION Fkstarrstari

!*************************************************************************************
! add Fkstarrstari with finite rotation Fkstarrstarirot, C. Bourdelle, from P. Cottier's QLK version
!***********************************************************************************

  COMPLEX(KIND=DBL) FUNCTION Fkstarrstarirot(ndim, xx, caseflag,nion)
    !---------------------------------------------------------------------
    ! Calculates the passing ion k*, r* integrands
    ! The returned integrand depends on the calculation flag, used
    ! in case we only want certain coefficient output
    ! caseflag = 1, full output
    ! caseflag = 2, return factor in front of Ati (kgti)
    ! caseflag = 3, return factor in front of Ani (kgni)    
    ! caseflag = 4, return curvature term (kgci)   
    ! caseflag = 5, return full output for energy integral (higher v exponent)
    ! caseflag = 6, return Ati coefficient for energy integral
    ! caseflag = 7, return Ani coefficient for energy integral
    ! caseflag = 8, return curvature term  for energy integral
    ! caseflag = 9, return factor in front of Aui (roto-diffusion)
    ! caseflag = 10, return factor in front of Aui (roto-diffusion) for energy integral
    ! caseflag = 11, return full output for ang momentum integral (higher v exponent)
    !---------------------------------------------------------------------  
    ! Arguments
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL)   , INTENT(IN) :: xx(ndim)
    INTEGER , INTENT(IN) :: caseflag, nion

    REAL(KIND=DBL)    :: Athir 
    COMPLEX(KIND=DBL) :: aai, bbi, cci, ddi, sqrtdi, sqrtdi2, Vmi, Vpi, Zai
    COMPLEX(KIND=DBL) :: alphai, ompFkr
    COMPLEX(KIND=DBL) :: faci
    REAL(KIND=DBL)    :: nwgi, nwE
    REAL(KIND=DBL)    :: rstar, kstar, xFkr, teta, fkstar
    COMPLEX(KIND=DBL) :: inti3, inti4, inti5, inti6, inti7 
    COMPLEX(KIND=DBL) :: Fikstarrstar,Fikstarrstar1,Fikstarrstar2,Fikstarrstar3,Fikstarrstartest
    REAL(KIND=DBL)    :: var2,var3,bessm2
    LOGICAL :: exist
    kstar = xx(1)
    rstar = xx(2)

    teta =(kstar/SQRT(2._DBL) /SQRT(REAL(mwidth**2.))+AIMAG(mshift)/REAL(mwidth**2.))*d
    ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.

    !Vertical drift term
    fkstar = 4./3.*(COS(teta) + (smag(pFkr) * teta - alphax(pFkr) * SIN(teta))*SIN(teta))

    IF (ABS(fkstar)<minfki) fkstar=SIGN(minfki,fkstar)

    nwgi = nwg*(-Tix(pFkr,nion)/Zi(pFkr,nion))

    var2 = (teta/d*Rhoi(pFkr,nion))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoi(nion))**2.               !!2nd argument of Bessel fun
    bessm2 = BESEI0(var2+var3)

    xFkr = rstar/SQRT(2.) + REAL(mshift)/SQRT(REAL(mwidth**2))
    !         & AIMAG(mwidth)/REAL(mwidth)*(kstar/SQRT(2.)+AIMAG(mshift)/SQRT(REAL(mwidth**2)))

    !Transit frequency        
    Athir = Athi(nion)*xFkr
    nwE = -gammaE(pFkr)/Rhoeff(pFkr)*SQRT(REAL(mwidth**2))
    ompFkr = omFkr-nwE*xFkr


    !GENERAL CASE only, do not keep options if fkstar=0 and/or xFkr=0 too unlikely

    ! Analytical form of velocity integration written with 
    ! Fried-Conte functions, generalizations of the simple plasma dispersion function

    bbi = CMPLX(Athir/(nwgi*fkstar),0.) ! in Pierre's case QLK_mom it was "-"Athir/(nwgi*fkstar) only in Fkstarrstar not in Fek... nor Fvk... nor Fkstarrstargni, etc 

    cci = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))/fkstar 

    ddi = bbi**2 - 4.*cci

    sqrtdi = SQRT(ddi)

    Vmi = (-bbi-sqrtdi)/2.
    Vpi = (-bbi+sqrtdi)/2.

    faci = 2. / (fkstar * (Vpi-Vmi) + epsD)

    IF ( (caseflag < 5).OR. (caseflag == 9) ) THEN !differentiate between particle, energy, momentum integrals, here particles     	  
       inti3 = faci * (Z1(Vpi) - Z1(Vmi))
       inti4 = faci * (Vpi*Z1(Vpi)-Vmi*Z1(Vmi)) 
       inti5 = faci * (Z2(Vpi) - Z2(Vmi))
       inti6 = faci * (Vpi*Z2(Vpi)-Vmi*Z2(Vmi))
       inti7 = faci * (Z3(Vpi) - Z3(Vmi))
    ELSEIF (caseflag == 11) THEN ! for ang momentum
       inti3 = faci * (Vpi*Z1(Vpi)-Vmi*Z1(Vmi)) 
       inti4 = faci * (Z2(Vpi) - Z2(Vmi))
       inti5 = faci * (Vpi*Z2(Vpi)-Vmi*Z2(Vmi))
       inti6 = faci * (Z3(Vpi) - Z3(Vmi))
       inti7 = faci * (Vpi*Z3(Vpi)-Vmi*Z3(Vmi))
    ELSE ! for energy
       inti3 = faci * (Z2(Vpi) - Z2(Vmi))
       inti4 = faci * (Vpi*Z2(Vpi)-Vmi*Z2(Vmi))
       inti5 = faci * (Z3(Vpi) - Z3(Vmi))
       inti6 = faci * (Vpi*Z3(Vpi)-Vmi*Z3(Vmi))
       inti7 = faci * (Z4(Vpi) - Z4(Vmi))
    ENDIF

    ! VERY DIRTY FIX! Otherwise QL integrals do not converge for some unexplained reason (maybe a NAG bug)
    IF ( (Machi(pFkr,nion)) < 1d-4 ) Machi(pFkr,nion) = 1d-4

    IF (caseflag == 1) THEN !full term for particle  
       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))+Ani(pFkr,nion)-1.5*Ati(pFkr,nion)
       Fikstarrstar = inti3 * ( alphai - 2.*Machi(pFkr,nion)*Aui(pFkr,nion)-Machi(pFkr,nion)**2*(alphai-Ati(pFkr,nion))) + &
            & inti4*2.*alam1*(Aui(pFkr,nion)+Machi(pFkr,nion)*(alphai - Ati(pFkr,nion))-3.*Machi(pFkr,nion)**2*Aui(pFkr,nion) ) + &
            & inti5*(Ati(pFkr,nion) + 4.*alam2*Machi(pFkr,nion)*Aui(pFkr,nion)+2.*Machi(pFkr,nion)**2*alam2*(alphai - (3.5+0.5/alam2)*Ati(pFkr,nion)))+ &
            & inti6*(2.*alam1*Machi(pFkr,nion)*Ati(pFkr,nion)+4.*Machi(pFkr,nion)**2*alam3*Aui(pFkr,nion)) + &
            & inti7*alam2*2.*Machi(pFkr,nion)**2*Ati(pFkr,nion)

    ELSEIF (caseflag == 2) THEN !At factor only particle transport
       alphai = -1.5
       Fikstarrstar = inti3 * ( alphai - Machi(pFkr,nion)**2*(alphai-1.)) + &
            & inti4*2.*alam1*(Machi(pFkr,nion)*(alphai - 1.)) + &
            & inti5*(1.+2.*Machi(pFkr,nion)**2*alam2*(alphai - (3.5+0.5/alam2)))+ &
            & inti6*(2.*alam1*Machi(pFkr,nion)) + &
            & inti7*alam2*2.*Machi(pFkr,nion)**2
    ELSEIF (caseflag == 3) THEN !An factor only particle transport
       alphai = 1.
       Fikstarrstar = inti3 * ( alphai -Machi(pFkr,nion)**2*alphai) + &
            & inti4*2.*alam1*Machi(pFkr,nion)*alphai + &
            & inti5*2.*Machi(pFkr,nion)**2*alam2*alphai
    ELSEIF (caseflag == 4) THEN !Curvature term particle transport
       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))
       Fikstarrstar = inti3 * ( alphai -Machi(pFkr,nion)**2*alphai) + &
            & inti4*2.*alam1*Machi(pFkr,nion)*alphai + &
            & inti5*2.*Machi(pFkr,nion)**2*alam2*alphai
    ELSEIF (caseflag == 9) THEN !Au term particle transport

       Fikstarrstar = -inti3 * 2.*Machi(pFkr,nion) + &
            & inti4*2.*alam1*(1.-3.*Machi(pFkr,nion)**2 ) + &
            & inti5*4.*alam2*Machi(pFkr,nion) + &
            & inti6*4.*Machi(pFkr,nion)**2*alam3  

    ELSEIF (caseflag == 5) THEN !full term for energy 
       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))+Ani(pFkr,nion)-1.5*Ati(pFkr,nion)

       Fikstarrstar = inti3 * ( alphai - 2.*Machi(pFkr,nion)*Aui(pFkr,nion)-Machi(pFkr,nion)**2*(alphai-Ati(pFkr,nion))) + &
            & inti4*2.*alam1*(Aui(pFkr,nion)+Machi(pFkr,nion)*(alphai - Ati(pFkr,nion))-3.*Machi(pFkr,nion)**2*Aui(pFkr,nion) ) + &
            & inti5*(Ati(pFkr,nion) + 4.*alam2*Machi(pFkr,nion)*Aui(pFkr,nion)+2.*Machi(pFkr,nion)**2*alam2*(alphai - (3.5+0.5/alam2)*Ati(pFkr,nion)))+ &
            & inti6*(2.*alam1*Machi(pFkr,nion)*Ati(pFkr,nion)+4.*Machi(pFkr,nion)**2*alam3*Aui(pFkr,nion)) + &
            & inti7*alam2*2.*Machi(pFkr,nion)**2*Ati(pFkr,nion)

    ELSEIF (caseflag == 6) THEN !At factor only energy transport
       alphai = -1.5
       Fikstarrstar = inti3 * ( alphai - Machi(pFkr,nion)**2*(alphai-1.)) + &
            & inti4*2.*alam1*(Machi(pFkr,nion)*(alphai - 1.)) + &
            & inti5*(1.+2.*Machi(pFkr,nion)**2*alam2*(alphai - (3.5+0.5/alam2)))+ &
            & inti6*(2.*alam1*Machi(pFkr,nion)) + &
            & inti7*alam2*2.*Machi(pFkr,nion)**2
    ELSEIF (caseflag == 7) THEN !An factor only energy transport
       alphai = 1.
       Fikstarrstar = inti3 * ( alphai -Machi(pFkr,nion)**2*alphai) + &
            & inti4*2.*alam1*Machi(pFkr,nion)*alphai + &
            & inti5*2.*Machi(pFkr,nion)**2*alam2*alphai    
    ELSEIF (caseflag == 8) THEN !Curvature term energy transport
       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))
       Fikstarrstar = inti3 * ( alphai -Machi(pFkr,nion)**2*alphai) + &
            & inti4*2.*alam1*Machi(pFkr,nion)*alphai + &
            & inti5*2.*Machi(pFkr,nion)**2*alam2*alphai
    ELSEIF (caseflag == 10) THEN !Au term energy transport
       Fikstarrstar = -inti3 * 2.*Machi(pFkr,nion) + &
            & inti4*2.*alam1*(1.-3.*Machi(pFkr,nion)**2 ) + &
            & inti5*4.*alam2*Machi(pFkr,nion) + &
            & inti6*4.*Machi(pFkr,nion)**2*alam3

    ELSEIF (caseflag == 11) THEN !full term for ang momentum

       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))+Ani(pFkr,nion)-1.5*Ati(pFkr,nion)

       Fikstarrstar = inti3 * alam1 * ( alphai - 2.*Machi(pFkr,nion)*Aui(pFkr,nion)-Machi(pFkr,nion)**2*(alphai-Ati(pFkr,nion))) + &
            & inti4*2.*alam2*(Aui(pFkr,nion)+Machi(pFkr,nion)*(alphai - Ati(pFkr,nion))-3.*Machi(pFkr,nion)**2*Aui(pFkr,nion) ) + &
            & inti5*(alam1*Ati(pFkr,nion) + 4.*alam3*Machi(pFkr,nion)*Aui(pFkr,nion)+2.*Machi(pFkr,nion)**2*alam3*(alphai - (3.5+alam1*0.5/alam3)*Ati(pFkr,nion)))+ &
            & inti6*(2.*alam2*Machi(pFkr,nion)*Ati(pFkr,nion)+4.*Machi(pFkr,nion)**2*alam4*Aui(pFkr,nion)) + &
            & inti7*alam3*2.*Machi(pFkr,nion)**2*Ati(pFkr,nion)


       !Pure rotation parts, proportional to Aui and Machi (and Machi**2) only! For testing
!!$       alphai = ompFkr*(Zi(pFkr,nion)/Tix(pFkr,nion))+Ani(pFkr,nion)-1.5*Ati(pFkr,nion)
!!$       Fikstarrstar = inti3 * alam1 * (- 2.*Machi(pFkr,nion)*Aui(pFkr,nion)-Machi(pFkr,nion)**2*(alphai-0.*Ati(pFkr,nion))) + &
!!$            & inti4*2.*alam2*(Aui(pFkr,nion)+Machi(pFkr,nion)*(alphai - 0.*Ati(pFkr,nion))-3.*Machi(pFkr,nion)**2*Aui(pFkr,nion) ) + &
!!$            & inti5*(4.*alam3*Machi(pFkr,nion)*Aui(pFkr,nion)+2.*Machi(pFkr,nion)**2*alam3*(alphai - (3.5+alam1*0.5/alam3)*Ati(pFkr,nion)))+ &
!!$            & inti6*(2.*alam2*Machi(pFkr,nion)*Ati(pFkr,nion)+4.*Machi(pFkr,nion)**2*alam4*Aui(pFkr,nion)) + &
!!$            & inti7*alam3*2.*Machi(pFkr,nion)**2*Ati(pFkr,nion)

    ENDIF

    IF ( (caseflag < 5).OR. (caseflag == 9) .OR. (caseflag == 11)) THEN !differentiate between particle, energy, ang. mom. integrals, here particle or momentum
!!$       Fkstarrstarirot = fc(pFkr)/twopi * EXP(-(kstar**2+rstar**2)/2) * coefi(pFkr,nion) * Fikstarrstar * Joi2c(nion) 
       Fkstarrstarirot = fc(pFkr)/twopi * EXP(-(kstar**2+rstar**2)/2) * coefi(pFkr,nion) * Fikstarrstar * bessm2 
       ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.
    ELSE !energy
!!$       Fkstarrstarirot = fc(pFkr)/twopi * EXP(-(kstar**2+rstar**2)/2) * coefi(pFkr,nion) * Tix(pFkr,nion) * Fikstarrstar *Joi2c(nion) 
       Fkstarrstarirot = fc(pFkr)/twopi * EXP(-(kstar**2+rstar**2)/2) * coefi(pFkr,nion) * Tix(pFkr,nion) * Fikstarrstar * bessm2 
       ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.
    ENDIF

    !set tiny values to 0 for smoothness (e.g. if Mach numbers are 0)
    IF (ABS(Fkstarrstarirot) < SQRT(epsD)) Fkstarrstarirot=0.

  END FUNCTION Fkstarrstarirot
  
!*************************************************************************************

  COMPLEX(KIND=DBL)  FUNCTION Fkstarrstare(ndim, xx, caseflag)
    !---------------------------------------------------------------------
    ! Calculate the f*, k* integrand for passing electrons
    ! The returned integrand depends on the calculation flag, used
    ! in case we only want certain coefficient output
    ! caseflag = 1, full output
    ! caseflag = 2, return factor in front of Ate (kgte)
    ! caseflag = 3, return factor in front of Ane (kgne)    
    ! caseflag = 4, return curvature term (kgce)   
    ! caseflag = 5, return full output for energy integral (higher v exponent)
    ! caseflag = 6, return Ate coefficient for energy integral
    ! caseflag = 7, return Ane coefficient for energy integral
    ! caseflag = 8, return curvature term  for energy integral
    !---------------------------------------------------------------------  
    ! Arguments
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL)   , INTENT(IN) :: xx(ndim)
    INTEGER , INTENT(IN) :: caseflag

    REAL(KIND=DBL)    :: Ather 
    COMPLEX(KIND=DBL) :: aae, bbe, cce, dde, sqrtde, Vme, Vpe, Zae
    COMPLEX(KIND=DBL) :: alphae
    COMPLEX(KIND=DBL) :: face
    REAL(KIND=DBL)    :: nwge,var2,var3,bessm2
    REAL(KIND=DBL)    :: rstar, kstar, teta, fkstar
    COMPLEX(KIND=DBL) :: inte3, inte5
    COMPLEX(KIND=DBL) :: Fekstarrstar
    COMPLEX(KIND=DBL) :: Febkstarrstar

    kstar = xx(1)
    rstar = xx(2)

    teta = kstar*d/REAL(mwidth)/SQRT(2._DBL)
    !Weighting term for vertical drift freq
    fkstar = 4./3.*(COS(teta) + (smag(pFkr) * teta - alphax(pFkr) * SIN(teta)) &
         * SIN(teta))
    IF (ABS(fkstar)<minfki) fkstar=SIGN(minfki,fkstar)

    !Vertical drift freq
    nwge = nwg*(-Tex(pFkr)/Ze)

    var2 = (teta/d*Rhoe(pFkr))**2  !!1st argument of Bessel fun
    var3 = ktetaRhoe**2               !!2nd argument of Bessel fun
    bessm2 = BESEI0(var2+var3)

    !Transit freq
    Ather = Athe*rstar/SQRT(2._DBL)

    !Simplified calc for zero vertical freq
    IF (ABS(fkstar)<minfki) THEN 
       !Further simplification for zero transit freq
       IF (rstar<epsD) THEN 

          inte3 = -1./omFkr*(-Tex(pFkr)/Ze)
          inte5 = 1.5*inte3

       ELSE   
          aae = omFkr*nwg/Ather
          Zae=Z1(aae)
          inte3 = nwge / Ather *2. *Zae
          inte5 = nwge / Ather *aae + inte3*aae*aae

       END IF
       !GENERAL CASE
    ELSE  
       ! Analytical form of velocity integration written with 
       ! Fried-Conte functions, generalizations of the simple plasma dispersion functio
       bbe = CMPLX(Ather/(nwge*fkstar),0.) 
       cce = omFkr*(Ze/Tex(pFkr))/fkstar

       dde = bbe**2 - 4.*cce
       sqrtde = SQRT(dde)

       Vme = (-bbe-sqrtde)/2.
       Vpe = (-bbe+sqrtde)/2.

       face = 2. / (fkstar * (Vpe-Vme)+epsD)
       IF (caseflag < 5) THEN !differentiate between particle or energy integrals
          inte3 = face * (Z1(Vpe) - Z1(Vme))
          inte5 = face * (Z2(Vpe) - Z2(Vme))
       ELSE
          inte3 = face * (Z2(Vpe) - Z2(Vme))
          inte5 = face * (Z3(Vpe) - Z3(Vme))
       ENDIF

    END IF

    IF ( (caseflag == 1) .OR. (caseflag == 5) ) THEN !full term for particle or energy
       alphae = omFkr*(Ze/Tex(pFkr))+Ane(pFkr)-1.5*Ate(pFkr)
       Fekstarrstar = inte3 * alphae + inte5 * Ate(pFkr)
    ELSEIF ( (caseflag == 2) .OR. (caseflag == 6) ) THEN
       alphae = -1.5
       Fekstarrstar = inte3 * alphae + inte5 
    ELSEIF ( (caseflag == 3) .OR. (caseflag == 7) ) THEN
       alphae = 1.
       Fekstarrstar = inte3 * alphae 
    ELSEIF ( (caseflag == 4) .OR. (caseflag == 8) ) THEN
       alphae = omFkr*(Ze/Tex(pFkr))
       Fekstarrstar = inte3 * alphae 
    ENDIF

    ! Care with norm factor (1 here, 4 in Fkstarrstar...)
    IF (caseflag < 5) THEN !differentiate between particle or energy integrals
!!$       Fkstarrstare = 1. * fc(pFkr) * Nex(pFkr) &
!!$            &  * (Fekstarrstar*Joe2c) * EXP( -(kstar**2+rstar**2)/2 )/twopi
       Fkstarrstare = 1. * fc(pFkr) * Nex(pFkr) &
            &  * (Fekstarrstar*bessm2) * EXP( -(kstar**2+rstar**2)/2 )/twopi
    ELSE
!!$       Fkstarrstare = 1. * fc(pFkr) *  Nex(pFkr) * Tex(pFkr) &
!!$            &  * (Fekstarrstar*Joe2c) * EXP( -(kstar**2+rstar**2)/2 )/twopi
       Fkstarrstare = 1. * fc(pFkr) *  Nex(pFkr) * Tex(pFkr) &
            &  * (Fekstarrstar*bessm2) * EXP( -(kstar**2+rstar**2)/2 )/twopi
    ENDIF
    IF (ABS(Fkstarrstare) < SQRT(epsD)) Fkstarrstare=0.

  END FUNCTION Fkstarrstare

!*************************************************************************************
! add Fkstarrstare with finite rotation Fkstarrstarerot, C. Bourdelle, from P. Cottier's QLK version
!***********************************************************************************

  COMPLEX(KIND=DBL)  FUNCTION Fkstarrstarerot(ndim, xx, caseflag)
    !---------------------------------------------------------------------
    ! Calculate the f*, k* integrand for passing electrons
    ! The returned integrand depends on the calculation flag, used
    ! in case we only want certain coefficient output
    ! caseflag = 1, full output
    ! caseflag = 2, return factor in front of Ate (kgte)
    ! caseflag = 3, return factor in front of Ane (kgne)    
    ! caseflag = 4, return curvature term (kgce)   
    ! caseflag = 5, return full output for energy integral (higher v exponent)
    ! caseflag = 6, return Ate coefficient for energy integral
    ! caseflag = 7, return Ane coefficient for energy integral
    ! caseflag = 8, return curvature term  for energy integral
    ! caseflag = 9, return factor in front of Aue (roto-diffusion)
    ! caseflag = 10, return factor in front of Aue (roto-diffusion) for energy integral
    ! caseflag = 11, return full output for ang momentum integral (higher v exponent)
    !---------------------------------------------------------------------  
    ! Arguments
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL)   , INTENT(IN) :: xx(ndim)
    INTEGER , INTENT(IN) :: caseflag

    REAL(KIND=DBL)    :: Ather 
    COMPLEX(KIND=DBL) :: aae, bbe, cce, dde, sqrtde, Vme, Vpe, Zae
    COMPLEX(KIND=DBL) :: alphae, ompFkr
    COMPLEX(KIND=DBL) :: face
    REAL(KIND=DBL)    :: nwge, nwE,var2,var3,bessm2
    REAL(KIND=DBL)    :: rstar, kstar, xFkr, teta, fkstar
    COMPLEX(KIND=DBL) :: inte3, inte4, inte5, inte6, inte7
    COMPLEX(KIND=DBL) :: Fekstarrstar

    kstar = xx(1)
    rstar = xx(2)

    teta =(kstar/SQRT(2._DBL)/SQRT(REAL(mwidth**2.))+AIMAG(mshift)/REAL(mwidth**2.))*d 
    ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.


    !Weighting term for vertical drift freq
    fkstar = 4./3.*(COS(teta) + (smag(pFkr) * teta - alphax(pFkr) * SIN(teta)) * SIN(teta))
    IF (ABS(fkstar)<minfki) fkstar=SIGN(minfki,fkstar)

    !Vertical drift freq
    nwge = nwg*(-Tex(pFkr)/Ze)

    var2 = (teta/d*Rhoe(pFkr))**2  !!1st argument of Bessel fun
    var3 = ktetaRhoe**2               !!2nd argument of Bessel fun
    bessm2 = BESEI0(var2+var3)

    !Transit frequency
    xFkr = rstar/SQRT(2.) + REAL(mshift)/SQRT(REAL(mwidth**2)) 
    Ather = Athe*xFkr
    nwE = -gammaE(pFkr)/Rhoeff(pFkr)*SQRT(REAL(mwidth**2))
    ompFkr = omFkr -nwE*xFkr

    !GENERAL CASE

    ! Analytical form of velocity integration written with 
    ! Fried-Conte functions, generalizations of the simple plasma dispersion functio

    bbe = CMPLX(Ather/(nwge*fkstar),0.) 
    cce = ompFkr*(Ze/Tex(pFkr))/fkstar

    dde = bbe**2 - 4.*cce
    sqrtde = SQRT(dde)

    Vme = (-bbe-sqrtde)/2.
    Vpe = (-bbe+sqrtde)/2.

    face = 2. / (fkstar * (Vpe-Vme) + epsD)

    IF ( (caseflag < 5).OR. (caseflag == 9) ) THEN !differentiate between particle, energy, momentum integrals, here particles     	  
       inte3 = face * (Z1(Vpe) - Z1(Vme))
       inte4 = face * (Vpe*Z1(Vpe)-Vme*Z1(Vme)) 
       inte5 = face * (Z2(Vpe) - Z2(Vme))
       inte6 = face * (Vpe*Z2(Vpe)-Vme*Z2(Vme))
       inte7 = face * (Z3(Vpe) - Z3(Vme))
    ELSEIF (caseflag == 11) THEN ! for ang momentum
       inte3 = face * (Vpe*Z1(Vpe)-Vme*Z1(Vme)) 
       inte4 = face * (Z2(Vpe) - Z2(Vme))
       inte5 = face * (Vpe*Z2(Vpe)-Vme*Z2(Vme))
       inte6 = face * (Z3(Vpe) - Z3(Vme))
       inte7 = face * (Vpe*Z3(Vpe)-Vme*Z3(Vme))
    ELSE ! for energy
       inte3 = face * (Z2(Vpe) - Z2(Vme))
       inte4 = face * (Vpe*Z2(Vpe)-Vme*Z2(Vme))
       inte5 = face * (Z3(Vpe) - Z3(Vme))
       inte6 = face * (Vpe*Z3(Vpe)-Vme*Z3(Vme))
       inte7 = face * (Z4(Vpe) - Z4(Vme))
    ENDIF

    IF (caseflag == 1) THEN !full term for particle 
       alphae = ompFkr*Ze/Tex(pFkr)+Ane(pFkr)-1.5*Ate(pFkr)
       Fekstarrstar = inte3 * alphae + inte5*Ate(pFkr)
    ELSEIF (caseflag == 2) THEN !At factor only particle transport
       alphae = -1.5
       Fekstarrstar = inte3 * alphae + inte5
    ELSEIF (caseflag == 3) THEN !An factor only particle transport
       alphae = 1.
       Fekstarrstar = inte3 * alphae
    ELSEIF (caseflag == 4) THEN !Curvature term particle transport
       alphae = ompFkr*Ze/Tex(pFkr)
       Fekstarrstar = inte3 * alphae  
    ELSEIF (caseflag == 9) THEN !Au term particle transport
       Fekstarrstar = 0.
    ELSEIF (caseflag == 5) THEN !full term for energy 
       alphae = ompFkr*Ze/Tex(pFkr)+Ane(pFkr)-1.5*Ate(pFkr)
       Fekstarrstar = inte3 * alphae + inte5*Ate(pFkr)          
    ELSEIF (caseflag == 6) THEN !At factor only energy transport
       alphae = -1.5
       alphae = -1.5
       Fekstarrstar = inte3 * alphae + inte5
    ELSEIF (caseflag == 7) THEN !An factor only energy transport
       alphae = 1.
       Fekstarrstar = inte3 * alphae 
    ELSEIF (caseflag == 8) THEN !Curvature term energy transport
       alphae = ompFkr*Ze/Tex(pFkr)
       Fekstarrstar = inte3 * alphae 
    ELSEIF (caseflag == 10) THEN !Au term energy transport
       Fekstarrstar = 0. 
    ELSEIF (caseflag == 11) THEN !full term for ang momentum
       alphae = ompFkr*Ze/Tex(pFkr)+Ane(pFkr)-1.5*Ate(pFkr)
       Fekstarrstar = alam1 * inte3 * alphae + inte5*Ate(pFkr) 
    ENDIF

    ! Care with norm factor (1 here, 4 in Fkstarrstar...)
    IF ( (caseflag < 5).OR. (caseflag == 9) .OR. (caseflag == 11)) THEN !differentiate between particle, energy, ang. mom. integrals, here particle or momentum
!!$       Fkstarrstarerot = fc(pFkr)/twopi * EXP( -(kstar**2 + rstar**2)/2 ) * Nex(pFkr) * Fekstarrstar*Joe2c 
       Fkstarrstarerot = fc(pFkr)/twopi * EXP( -(kstar**2 + rstar**2)/2 ) * Nex(pFkr) * Fekstarrstar*bessm2
       ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.
    ELSE !energy
!!$       Fkstarrstarerot = fc(pFkr)/twopi * EXP( -(kstar**2 + rstar**2)/2 ) * Nex(pFkr) * Tex(pFkr) * Fekstarrstar*Joe2c 
       Fkstarrstarerot = fc(pFkr)/twopi * EXP( -(kstar**2 + rstar**2)/2 ) * Nex(pFkr) * Tex(pFkr) * Fekstarrstar*bessm2
       ! in order to have 1/2pi exp( -(kstar^2+rstar^2)/2) Jonathan Citrin changed the normalizations compare to Pierre Cottier QLK_mom's version. Keep JC norm.
    ENDIF

    IF (ABS(Fkstarrstarerot) < SQRT(epsD)) Fkstarrstarerot=0.

  END FUNCTION Fkstarrstarerot

END MODULE passints
