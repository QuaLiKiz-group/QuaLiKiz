MODULE trapints
  USE kind
  USE datcal
  USE datmat
  USE dispfuncs

  IMPLICIT NONE

CONTAINS
  COMPLEX(KIND=DBL) FUNCTION FFki(kk,caseflag,nion)
    !---------------------------------------------------------------------
    ! Integrand for trapped ions
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
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER, INTENT(IN) :: caseflag, nion
    COMPLEX(KIND=DBL) :: Fik, Fik1, zik, fk
    COMPLEX(KIND=DBL) :: bbip, ddip, ddip2, Vmip, Vpip
    COMPLEX(KIND=DBL) :: zik2, Zgik
    COMPLEX(KIND=DBL) :: Aiz, Biz, Ciz
    REAL(KIND=DBL)    :: ya, k2, nwgi
    REAL(KIND=DBL)    :: fki, Eg, Kg

    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)

    fk = CMPLX(fki,0)

    ! The frequency is weighted by the vertical frequency drift    
    zik2 = omFFk * (-Zi(pFFk,nion)/Tix(pFFk,nion))/fk

    ! The square root whose imaginary part is positive is chosen 
    ! for compatibility with the Fried-Conte function
    zik  = SQRT(zik2) 
    IF (AIMAG(zik)<0.) zik = -zik

    !The plasma dispersion function is calculated (analytical energy integration)
    Zgik = ci * sqrtpi * wofzweid(zik)

    !Now the function is calculated based on the dispersion function calculation
    Aiz = 1. + zik * Zgik !Z1/z 
    Biz = 0.5 + zik2 * Aiz  !Z2/z              
    Ciz =  0.75 + zik2 * Biz !Z3/z

    nwgi = nwg*(-Tix(pFFk,nion)/Zi(pFFk,nion))
    bbip = CMPLX(cthi(pFFk,nion)*omega2bar/(Kg*fki*nwgi*qx(pFFk)*Ro(pFFk)),0.)
    ddip2 = bbip**2 + 4.*zik2
    ddip =SQRT(ddip2)
    Vmip = (- bbip - ddip)/2. !used in traporder1 (therefore not used)
    Vpip = (- bbip + ddip)/2. !used in traporder1 (therefore not used)

    !Different cases are calculated depending on desired output
    !Fik1 corresponds to the higher order response
    IF (caseflag == 1) THEN
       Fik = 2.*((-zik2 - 1.5 * Ati(pFFk,nion)/fk+Ani(pFFk,nion)/fk)*Aiz + Ati(pFFk,nion)/fk*Biz)
    ELSEIF (caseflag == 2) THEN
       Fik = 2.*( (-1.5/fk) * Aiz + 1./fk*Biz)
    ELSEIF (caseflag == 3) THEN
       Fik = 2.*(1./fk) * Aiz
    ELSEIF (caseflag == 4) THEN
       Fik = -2.* zik2 * Aiz
    ELSEIF (caseflag == 5) THEN !Energy integral
       Fik = 2.*((-zik2 - 1.5 * Ati(pFFk,nion)/fk+Ani(pFFk,nion)/fk)*Biz + Ati(pFFk,nion)/fk*Ciz)
    ELSEIF (caseflag == 6) THEN
       Fik = 2.*( (-1.5/fk) * Biz + 1./fk*Ciz)
    ELSEIF (caseflag == 7) THEN
       Fik = 2.*(1./fk) * Biz
    ELSEIF (caseflag == 8) THEN
       Fik = -2.* zik2 * Biz
    ENDIF

    IF (traporder1 .EQV. .TRUE.) THEN
       IF (caseflag == 1) THEN
          Fik1 = 2.*(-zik2 - 1.5 * Ati(pFFk,nion)/fk+Ani(pFFk,nion)/fk) * &
               (Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)+Ati(pFFk,nion)/fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 2) THEN
          Fik1 = -3./fk*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip) + 1./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 3) THEN
          Fik1 = 2./fk*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 4) THEN
          Fik1 = -2.*zik2*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 5) THEN !Energy integral
          Fik1 = 2.*(-zik2 - 1.5 * Ati(pFFk,nion)/fki+Ani(pFFk,nion)/fk) * &
               (Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)+Ati(pFFk,nion)/fk*(Z3(Vpip)-Z3(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 6) THEN
          Fik1 = -3./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip) + 1./fk*(Z3(Vpip)-Z3(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 7) THEN
          Fik1 = 2./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 8) THEN
          Fik1 = -2.*zik2*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ENDIF
    ELSE
       Fik1 = 0.
    ENDIF

    IF (caseflag < 5) THEN
       FFki = kk * Kg * ft(pFFk) *  coefi(pFFk,nion) * (Fik  * Joi2p(nion) + Fik1 * J1i2p(nion)) 
    ELSE
       FFki = kk * Kg * ft(pFFk) *  coefi(pFFk,nion) * Tix(pFFk,nion) * (Fik  * Joi2p(nion) + Fik1 * J1i2p(nion)) 
    ENDIF

    IF (ABS(FFki) < SQRT(epsD)) FFki=0.

  END FUNCTION FFki

!***********************************************************************************
! add FFki with finite rotation FFkirot, C. Bourdelle, from P. Cottier's QLK version
!***********************************************************************************

  COMPLEX(KIND=DBL) FUNCTION FFkirot(kk,caseflag,nion)
    !---------------------------------------------------------------------
    ! Integrand for trapped ions
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
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER, INTENT(IN) :: caseflag, nion
    COMPLEX(KIND=DBL) :: Fik, Fik1, zik, fk
    COMPLEX(KIND=DBL) :: bbip, ddip, ddip2, Vmip, Vpip
    COMPLEX(KIND=DBL) :: zik2, Zgik
    COMPLEX(KIND=DBL) :: Aiz, Biz, Ciz
    COMPLEX(KIND=DBL) :: ompFFk
    REAL(KIND=DBL)    :: ya, k2, nwgi, nwe,vpar2
    REAL(KIND=DBL)    :: fki, Eg, E2g, Kg, gau2mshift

    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)

    fk = CMPLX(fki,0)

    ! account for gamma_E=-grad(E_r)/B*R/Cref
    ! nwE means here nwE'=-k_theta*grad(Er)/B=k_theta*gammaE(as def in code)*Cref/R normalized to nwd=-ktheta T(1keV)/(eBR)

    gau2mshift = REAL(mshift) + 2.*AIMAG(mshift)*AIMAG(mwidth)*REAL(mwidth)/(REAL(mwidth)**2 - AIMAG(mwidth)**2)

    nwE = -gammaE(pFFk)/Rhoeff(pFFk)
    ompFFk =  omFFk-nwE*gau2mshift ! gau2mshift is shift of |gau(x)| 

    ! The frequency is weighted by the vertical frequency drift    
    zik2 = ompFFk * (-Zi(pFFk,nion)/Tix(pFFk,nion))/fk

    ! The square root whose imaginary part is positive is chosen 
    ! for compatibility with the Fried-Conte function
    zik  = SQRT(zik2) 
    IF (AIMAG(zik)<0.) zik = -zik

    !The plasma dispersion function is calculated (analytical energy integration)
    Zgik = ci * sqrtpi * wofzweid(zik)

    !Now the function is calculated based on the dispersion function calculation
    Aiz = 1. + zik * Zgik !Z1/z
    Biz = 0.5 + zik2 * Aiz  !Z2/z              
    Ciz =  0.75 + zik2 * Biz !Z3/z

    !Bounce average of vpar^2 / (vts^2 * v^2)
    vpar2 = 2.*epsilon(pFFk)*kk/Kg*E2g

!!$    nwgi = nwg*(-Tix(pFFk,nion)/Zi(pFFk,nion))
!!$    bbip = CMPLX(cthi(pFFk,nion)*omega2bar/(Kg*fki*nwgi*qx(pFFk)*Ro(pFFk)),0.)
!!$    ddip2 = bbip**2 + 4.*zik2
!!$    ddip =SQRT(ddip2)
!!$    Vmip = (- bbip - ddip)/2. !used in traporder1 (therefore not used)
!!$    Vpip = (- bbip + ddip)/2  !used in traporder1 (therefore not used)

    !Different cases are calculated depending on desired output. Only terms up to Mach**2 are kept, where Vpar^2 and Aui are ordered like Mach.
    IF (caseflag == 1) THEN
       Fik = 2./fk*(  (1.-Machi(pFFk,nion)**2)*( (-zik2*fk - 1.5 * Ati(pFFk,nion)+Ani(pFFk,nion))*Aiz + Ati(pFFk,nion)*Biz) &
            & + (Machi(pFFk,nion)**2*Ati(pFFk,nion) - 2.*Machi(pFFk,nion)*Aui(pFFk,nion) )*Aiz + 4.*Machi(pFFk,nion)*Aui(pFFk,nion)*vpar2*Biz)
    ELSEIF (caseflag == 2) THEN
       Fik = 2./fk*( (1.- Machi(pFFk,nion)**2)*(-1.5*Aiz+Biz) + Machi(pFFk,nion)**2*Aiz)
    ELSEIF (caseflag == 3) THEN
       Fik = 2./fk * Aiz * (1-Machi(pFFk,nion)**2) 
    ELSEIF (caseflag == 9) THEN
       Fik = -4./fk* Machi(pFFk,nion)*Aiz + 8./fk*Machi(pFFk,nion)*vpar2*Biz
    ELSEIF (caseflag == 4) THEN
       Fik = -2.* zik2 * Aiz * (1-Machi(pFFk,nion)**2) 
    ELSEIF (caseflag == 5) THEN 
       Fik = 2./fk*(  (1.-Machi(pFFk,nion)**2)*( (-zik2*fk - 1.5 * Ati(pFFk,nion)+Ani(pFFk,nion))*Biz + Ati(pFFk,nion)*Ciz) &
            & + (Machi(pFFk,nion)**2*Ati(pFFk,nion) - 2.*Machi(pFFk,nion)*Aui(pFFk,nion) )*Biz + 4.*Machi(pFFk,nion)*Aui(pFFk,nion)*vpar2*Ciz)
    ELSEIF (caseflag == 6) THEN
       Fik = 2./fk*( (1.- Machi(pFFk,nion)**2)*(-1.5*Biz+Ciz) + Machi(pFFk,nion)**2*Biz)
    ELSEIF (caseflag == 7) THEN
       Fik = 2./fk * Biz * (1.-Machi(pFFk,nion)**2) 
    ELSEIF (caseflag == 8) THEN
       Fik = -2.* zik2 * Biz*(1.-Machi(pFFk,nion)**2) 
    ELSEIF (caseflag == 10) THEN
       Fik = -4./fk* Machi(pFFk,nion)*Biz + 8./fk*Machi(pFFk,nion)*vpar2*Ciz
    ELSEIF (caseflag == 11) THEN  ! integral for ang momentum
       ! I do not agree fully with Pierre version. New version with bounce average of Vpar2/vts2 

!!$       Fik =  4. *omega2bar/(Kg*fk)*(Aui(pFFk,nion)+Machi(pFFk,nion)*(Ani(pFFk,nion)-2.5*Ati(pFFk,nion)-zik2*fk))*Aiz*zik & 
!!$            & +Machi(pFFk,nion)*Ati(pFFk,nion)*Biz*zik

       Fik = 4./fk * vpar2 * (Biz*(Aui(pFFk,nion)+Machi(pFFk,nion)*(Ani(pFFk,nion)-2.5*Ati(pffK,nion)-zik2*fk) ) + &
            & Ciz*Machi(pFFk,nion)*Ati(pFFk,nion)   )
    ENDIF

    !Fik1 corresponds to the higher order response
    ! With rotation do not allow for traporder1 terms! Anyway not used in Pierre's version   
    Fik1 = 0.

    IF ( (caseflag < 5) .OR. (caseflag == 9) .OR. (caseflag == 11) ) THEN
       FFkirot = kk * Kg * ft(pFFk) *  coefi(pFFk,nion) * (Fik  * Joi2p(nion) + Fik1 * J1i2p(nion))
    ELSE
       FFkirot = kk * Kg * ft(pFFk) *  coefi(pFFk,nion) * Tix(pFFk,nion) * (Fik  * Joi2p(nion) + Fik1 * J1i2p(nion))
    ENDIF
    IF (ABS(FFkirot) < SQRT(epsD)) FFkirot=0.
  END FUNCTION FFkirot

!*************************************************************************************************************************************************************

  COMPLEX(KIND=DBL) FUNCTION FFke(ndim, kv, caseflag)
    !---------------------------------------------------------------------
    ! Calculate the trapped electron integrand, including both kappa and v
    ! Includes collisions (Krook operator)
    ! integration in v only for v > 0 (energy) so factor 4
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: kv
    INTEGER, INTENT(IN) :: caseflag

    COMPLEX(KIND=DBL) :: Fekv, Fekv1, fk
    COMPLEX(KIND=DBL) :: Aez, Bez, Bez1, Bez2, zek2, zek, Zgek, bbe
    REAL(KIND=DBL)    :: v, v2, v3, v4, v5
    REAL(KIND=DBL)    :: ya, k2, kk, nwge
    REAL(KIND=DBL)    :: fki, Eg, Kg, delta, Anuen, Anuent

    kk = kv(1)
    v = kv(2)
    k2 = kk*kk

    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)

    fk = CMPLX(fki,0)

    v2 = v*v 
    v3 = v2*v 
    v4 = v3*v
    v5 = v4*v

    ! Calculation of the electron collision frequency
    ! odd function, integration not ok if Anuen =0, so put Anuen = 1.0e-14
    !we normalize the collision frequency to -nw_ge
    Anuen = Anue(pFFk) * (-Ze) / (Tex(pFFk)*nwg)

    ! Krook operator for collisions. From M. Kotschenreuther, G. Rewoldt and W.M. Tang, Computer Physics Communications 88 (1995) 128-140
    delta = ( ABS(omFFk) * nwg / (Anue(pFFk) * 37.2))**(1./3.)
    Anuent = Anuen / ((2.*k2 -1.)**2) * (0.111 * delta +1.31) / (11.79 * delta + 1.) 
    
    IF ( ABS(Anuent) < epsD ) THEN
       Anuent = epsD
    ENDIF

    zek2 = omFFk * (-Ze/Tex(pFFk)) !makes the normalization of omega to nwge instead of nwg
    nwge = nwg*(-Tex(pFFk)/Ze)

    bbe = CMPLX(cthe(pFFk)*omega2bar/(Kg*nwge*qx(pFFk)*Ro(pFFk)),0.)

    IF ( (caseflag == 1) .OR. (caseflag == 5) ) THEN
       Aez = (zek2 + 1.5 * Ate(pFFk) - Ane(pFFk)) * v3 - Ate(pFFk) * v5  
    ELSEIF ( (caseflag == 2 ) .OR. (caseflag == 6) ) THEN
       Aez =  1.5 * v3 - v5  
    ELSEIF ( (caseflag == 3) .OR. (caseflag == 7) ) THEN
       Aez = -v3
    ELSEIF ( (caseflag == 4) .OR. (caseflag == 8) ) THEN
       Aez = zek2*v3
    ENDIF

    Bez =  zek2 * v3 - v5*fk + ci * Anuent
    Bez1 =  zek2*v3  - bbe*v4 - v5*fk + ci * Anuent       
    Bez2 =  zek2*v3  + bbe*v4 - v5*fk + ci * Anuent 

    Fekv = 4. / sqrtpi * v2 * EXP(-v2) * Aez / Bez
    IF (traporder1 .EQV. .TRUE.) THEN
       Fekv1 = 2. / sqrtpi * v2 * EXP(-v2) * (Aez / Bez1 + Aez / Bez2)
    ELSE
       Fekv1=0.
    ENDIF

    IF (caseflag < 5) THEN
       FFke = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fekv  * Joe2p +Fekv1 * J1e2p) 
    ELSE
       FFke = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * v2 * (Fekv  * Joe2p + Fekv1 * J1e2p) 
    ENDIF
    IF (ABS(FFke) < SQRT(epsD)) FFke=0.
  END FUNCTION FFke

  COMPLEX(KIND=DBL) FUNCTION FFke_nocoll(kk,caseflag)
    !---------------------------------------------------------------------
    ! Integrand for trapped electrons when no collisionality is included
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
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER, INTENT(IN) :: caseflag
    COMPLEX(KIND=DBL) :: Fik, Fik1, zik, fk
    COMPLEX(KIND=DBL) :: bbip, ddip, ddip2, Vmip, Vpip
    COMPLEX(KIND=DBL) :: zik2, Zgik
    COMPLEX(KIND=DBL) :: Aiz, Biz, Ciz
    REAL(KIND=DBL)    :: ya, k2, nwge
    REAL(KIND=DBL)    :: fki, Eg, Kg
   
    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))
    fk = CMPLX(fki,0)
    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)
    ! The frequency is weighted by the vertical frequency drift    
    zik2 = omFFk * (1./Tex(pFFk))/fk

    ! The square root whose imaginary part is positive is chosen 
    ! for compatibility with the Fried-Conte function
    zik  = SQRT(zik2) 
    IF (AIMAG(zik)<0.) zik = -zik

    !The plasma dispersion function is calculated (analytical energy integration)
    Zgik = ci * sqrtpi * wofzweid(zik)

    !Now the function is calculated based on the dispersion function calculation
    Aiz = 1. + zik * Zgik !Z1
    Biz = 0.5 + zik2 * Aiz  !Z2              
    Ciz =  0.75 + zik2 * Biz

                   
    nwge = nwg*Tex(pFFk)
    bbip = CMPLX(cthe(pFFk)*omega2bar/(Kg*fki*nwge*qx(pFFk)*Ro(pFFk)),0.)
    ddip2 = bbip**2 + 4.*zik2
    ddip =SQRT(ddip2)
    Vmip = (- bbip - ddip)/2. !used in traporder1 (therefore not used)
    Vpip = (- bbip + ddip)/2. !used in traporder1 (therefore not used)

    !Different cases are calculated depending on desired output
    !Fik1 corresponds to the higher order response
    IF (caseflag == 1) THEN
       Fik = 2.*((-zik2 - 1.5 * Ate(pFFk)/fk+Ane(pFFk)/fk)*Aiz + Ate(pFFk)/fk*Biz)
    ELSEIF (caseflag == 2) THEN
       Fik = 2.*( (-1.5/fk) * Aiz + 1./fk*Biz)
    ELSEIF (caseflag == 3) THEN
       Fik = 2.*(1./fk) * Aiz
    ELSEIF (caseflag == 4) THEN
       Fik = -2.* zik2 * Aiz
    ELSEIF (caseflag == 5) THEN !Energy integral
       Fik = 2.*((-zik2 - 1.5 * Ate(pFFk)/fk+Ane(pFFk)/fk)*Biz + Ate(pFFk)/fk*Ciz)
    ELSEIF (caseflag == 6) THEN
       Fik = 2.*( (-1.5/fk) * Biz + 1./fk*Ciz)
    ELSEIF (caseflag == 7) THEN
       Fik = 2.*(1./fk) * Biz
    ELSEIF (caseflag == 8) THEN
       Fik = -2.* zik2 * Biz
    ENDIF

    IF (traporder1 .EQV. .TRUE.) THEN
       IF (caseflag == 1) THEN
          Fik1 = 2.*(-zik2 - 1.5 * Ate(pFFk)/fk+Ane(pFFk)/fk) * &
               (Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)+Ate(pFFk)/fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 2) THEN
          Fik1 = -3./fk*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip) + 1./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 3) THEN
          Fik1 = 2./fk*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 4) THEN
          Fik1 = -2.*zik2*(Z1(Vpip)-Z1(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 5) THEN !Energy integral
          Fik1 = 2.*(-zik2 - 1.5 * Ate(pFFk)/fki+Ane(pFFk)/fk) * &
               (Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)+Ate(pFFk)/fk*(Z3(Vpip)-Z3(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 6) THEN
          Fik1 = -3./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip) + 1./fk*(Z3(Vpip)-Z3(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 7) THEN
          Fik1 = 2./fk*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ELSEIF (caseflag == 8) THEN
          Fik1 = -2.*zik2*(Z2(Vpip)-Z2(Vmip))/(Vpip-Vmip)
       ENDIF
    ELSE
       Fik1 = 0.
    ENDIF

    IF (caseflag < 5) THEN
       FFke_nocoll = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ELSE
       FFke_nocoll = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ENDIF
    IF (ABS(FFke_nocoll) < SQRT(epsD)) FFke_nocoll=0.
  END FUNCTION FFke_nocoll

!*************************************************************************************

!***********************************************************************************
! add FFke with finite rotation FFkerot, C. Bourdelle, from P. Cottier's QLK version
!***********************************************************************************

  COMPLEX(KIND=DBL) FUNCTION FFkerot(ndim, kv, caseflag)
    !---------------------------------------------------------------------
    ! Calculate the trapped electron integrand, including both kappa and v
    ! Includes collisions (Krook operator)
    ! integration in v only for v > 0 (energy) so factor 4
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

    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: kv
    INTEGER, INTENT(IN) :: caseflag

    COMPLEX(KIND=DBL) :: Fekv, Fekv1, fk, ompFFk
    COMPLEX(KIND=DBL) :: Aez, Bez, Bez1, Bez2, zek2, zek, Zgek
    REAL(KIND=DBL)    :: v, v2, v3, v4, v5,gau2mshift
    REAL(KIND=DBL)    :: ya, k2, kk, nwge, nwe
    REAL(KIND=DBL)    :: fki, Eg, E2g, Kg, delta, Anuen, Anuent,vpar2

    kk = kv(1)
    v = kv(2)
    k2 = kk*kk

    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))
    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)
    fk = CMPLX(fki,0)

    v2 = v*v 
    v3 = v2*v 
    v4 = v3*v
    v5 = v4*v

    ! Calculation of the electron collision frequency
    ! odd function, integration not ok if Anuen =0, so put Anuen = 1.0e-14
    !we normalize the collision frequency to -nw_ge
    Anuen = Anue(pFFk) * (-Ze) / (Tex(pFFk)*nwg)

    gau2mshift = REAL(mshift) + 2.*AIMAG(mshift)*AIMAG(mwidth)*REAL(mwidth)/(REAL(mwidth)**2 - AIMAG(mwidth)**2)

    ! account for gamma_E
    nwE = -gammaE(pFFk)/Rhoeff(pFFk)
    ompFFk = omFFk-nwE*gau2mshift 

    ! Krook operator for collisions. From M. Kotschenreuther, G. Rewoldt and W.M. Tang, Computer Physics Communications 88 (1995) 128-140
    delta = ( ABS(ompFFk) * nwg / (Anue(pFFk) * 37.2))**(1./3.)
    Anuent = Anuen / ((2.*k2 -1.)**2) * (0.111 * delta +1.31) / (11.79 * delta + 1.) 

    IF ( ABS(Anuent) < epsD ) THEN
       Anuent = epsD
    ENDIF

    zek2 = ompFFk * (-Ze/Tex(pFFk)) !makes the normalization of omega to nwge instead of nwg
    nwge = nwg*(-Tex(pFFk)/Ze)
    !Bounce average of vpar^2 / (vts^2 * v^2)
    vpar2 = 2.*epsilon(pFFk)*kk/Kg*E2g

    IF ( (caseflag == 1) .OR. (caseflag == 5) ) THEN
       Aez = (zek2 + 1.5 * Ate(pFFk)-Ane(pFFk))*v3 - Ate(pFFk)*v5
    ELSEIF ( (caseflag == 2 ) .OR. (caseflag == 6) ) THEN
       Aez = 1.5*v3-v5
    ELSEIF ( (caseflag == 3) .OR. (caseflag == 7) ) THEN
       Aez = -v3 
    ELSEIF ( (caseflag == 4) .OR. (caseflag == 8) ) THEN
       Aez =  zek2 * v3
    ELSEIF ( (caseflag == 9) .OR. (caseflag == 10) ) THEN
       Aez = 0.
    ELSEIF (caseflag == 11) THEN
       !       Aez = 2.* omega2bar/Kg*v4*(-Aue(pFFk)+Mache(pFFk)*(zek2 + 2.5 * Ate(pFFk) - Ane(pFFk)-Ate(pFFk)*v2))
       Aez = 0. !Tiny quantity since Mache, Aue, and epsilon all small
    ENDIF

    Bez =  zek2 * v3 - v5*fk + ci * Anuent

    Fekv = 4. / sqrtpi * v2 * EXP(-v2) * Aez / Bez
    Fekv1=0.

    IF ( (caseflag < 5) .OR. (caseflag == 9) .OR. (caseflag == 11) ) THEN
       FFkerot = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fekv  * Joe2p +Fekv1 * J1e2p) 
    ELSE
       FFkerot = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * v2 * (Fekv  * Joe2p + Fekv1 * J1e2p) 
    ENDIF
       IF (ABS(FFkerot) < SQRT(epsD)) FFkerot=0.
  END FUNCTION FFkerot

  COMPLEX(KIND=DBL) FUNCTION FFkerotold(ndim, kv, caseflag)
    !---------------------------------------------------------------------
    ! Calculate the trapped electron integrand, including both kappa and v
    ! Includes collisions (Krook operator)
    ! integration in v only for v > 0 (energy) so factor 4
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

    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: kv
    INTEGER, INTENT(IN) :: caseflag

    COMPLEX(KIND=DBL) :: Fekv, Fekv1, fk, ompFFk
    COMPLEX(KIND=DBL) :: Aez, Bez, Bez1, Bez2, zek2, zek, Zgek
    REAL(KIND=DBL)    :: v, v2, v3, v4, v5,gau2mshift
    REAL(KIND=DBL)    :: ya, k2, kk, nwge, nwe
    REAL(KIND=DBL)    :: fki, Eg, E2g, Kg, delta, Anuen, Anuent,vpar2

    kk = kv(1)
    v = kv(2)
    k2 = kk*kk

    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))
    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)
    fk = CMPLX(fki,0)

    v2 = v*v 
    v3 = v2*v 
    v4 = v3*v
    v5 = v4*v

    ! Calculation of the electron collision frequency
    ! odd function, integration not ok if Anuen =0, so put Anuen = 1.0e-14
    !we normalize the collision frequency to -nw_ge
    Anuen = Anue(pFFk) * (-Ze) / (Tex(pFFk)*nwg)

    gau2mshift = REAL(mshift) + 2.*AIMAG(mshift)*AIMAG(mwidth)*REAL(mwidth)/(REAL(mwidth)**2 - AIMAG(mwidth)**2)

    ! account for gamma_E
    nwE = -gammaE(pFFk)/Rhoeff(pFFk)
    ompFFk = omFFk-nwE*gau2mshift 

    ! Krook operator for collisions. From M. Kotschenreuther, G. Rewoldt and W.M. Tang, Computer Physics Communications 88 (1995) 128-140
    delta = ( ABS(ompFFk) * nwg / (Anue(pFFk) * 37.2))**(1./3.)
    Anuent = Anuen / ((2.*k2 -1.)**2) * (0.111 * delta +1.31) / (11.79 * delta + 1.) 

    IF ( ABS(Anuent) < epsD ) THEN
       Anuent = epsD
    ENDIF

    zek2 = ompFFk * (-Ze/Tex(pFFk)) !makes the normalization of omega to nwge instead of nwg
    nwge = nwg*(-Tex(pFFk)/Ze)
    !Bounce average of vpar^2 / (vts^2 * v^2)
    vpar2 = 2.*epsilon(pFFk)*kk/Kg*E2g

    IF ( (caseflag == 1) .OR. (caseflag == 5) ) THEN
       Aez = (1.-Mache(pFFk)**2)*( (zek2 + 1.5 * Ate(pFFk)-Ane(pFFk))*v3 - Ate(pFFk)*v5) &
            & - (Mache(pFFk)**2*Ate(pFFk) + 2.*Mache(pFFk)*Aue(pFFk) )*v3 + 4.*Mache(pFFk)*Aue(pFFk)*vpar2*v5
    ELSEIF ( (caseflag == 2 ) .OR. (caseflag == 6) ) THEN
       Aez = (1.- Mache(pFFk)**2)*(1.5*v3-v5) - Mache(pFFk)**2*v3
    ELSEIF ( (caseflag == 3) .OR. (caseflag == 7) ) THEN
       Aez = -(1 - Mache(pFFk)**2) * v3 
    ELSEIF ( (caseflag == 4) .OR. (caseflag == 8) ) THEN
       Aez =  zek2 * v3 *(1-Mache(pFFk)**2)
    ELSEIF ( (caseflag == 9) .OR. (caseflag == 10) ) THEN
       Aez = 2.*Mache(pFFk)* v3 + 4.*Mache(pFFk)*vpar2*v5
    ELSEIF (caseflag == 11) THEN
       !       Aez = 2.* omega2bar/Kg*v4*(-Aue(pFFk)+Mache(pFFk)*(zek2 + 2.5 * Ate(pFFk) - Ane(pFFk)-Ate(pFFk)*v2))
       Aez = 2.* vpar2 * v5*(-Aue(pFFk)+Mache(pFFk)*(-Ane(pFFk)+2.5*Ate(pFFk)+zek2 - v2*Mache(pFFk)*Ate(pFFk) ))
    ENDIF

    Bez =  zek2 * v3 - v5*fk + ci * Anuent

    Fekv = 4. / sqrtpi * v2 * EXP(-v2) * Aez / Bez
    Fekv1=0.

    IF ( (caseflag < 5) .OR. (caseflag == 9) .OR. (caseflag == 11) ) THEN
       FFkerotold = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fekv  * Joe2p +Fekv1 * J1e2p) 
    ELSE
       FFkerotold = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * v2 * (Fekv  * Joe2p + Fekv1 * J1e2p) 
    ENDIF
       IF (ABS(FFkerotold) < SQRT(epsD)) FFkerotold=0.
  END FUNCTION FFkerotold


  COMPLEX(KIND=DBL) FUNCTION FFke_nocollrot(kk,caseflag)
    !---------------------------------------------------------------------
    ! Integrand for trapped electrons when no collisionality and rotation is included
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
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER, INTENT(IN) :: caseflag
    COMPLEX(KIND=DBL) :: Fik, Fik1, zik, fk
    COMPLEX(KIND=DBL) :: bbip, ddip, ddip2, Vmip, Vpip
    COMPLEX(KIND=DBL) :: zik2, Zgik
    COMPLEX(KIND=DBL) :: Aiz, Biz, Ciz
    COMPLEX(KIND=DBL) :: ompFFk
    REAL(KIND=DBL)    :: ya, k2, nwge, nwe
    REAL(KIND=DBL)    :: fki, Eg, E2g, Kg, gau2mshift,vpar2

    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2
    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)

    fk = CMPLX(fki,0)

    ! account for gamma_E=-grad(E_r)/B*R/Cref
    ! nwE means here nwE'=-k_theta*grad(Er)/B=k_theta*gammaE(as def in code)*Cref/R normalized to nwd=-ktheta T(1keV)/(eBR)
!    gau2mshift = REAL(mshift) + 2.*AIMAG(mshift)*AIMAG(mwidth)*REAL(mwidth)/(REAL(mwidth)**2 - AIMAG(mwidth)**2)

    gau2mshift = REAL(mshift) + AIMAG(mshift)*AIMAG(mwidth**2)/REAL(mwidth**2)

    nwE = -gammaE(pFFk)/Rhoeff(pFFk)
    ompFFk =  omFFk-nwE*gau2mshift ! gau2mshift is shift of |gau(x)| 

    ! The frequency is weighted by the vertical frequency drift    
    zik2 = ompFFk * (1./Tex(pFFk))/fk

    ! The square root whose imaginary part is positive is chosen 
    ! for compatibility with the Fried-Conte function
    zik  = SQRT(zik2) 
    IF (AIMAG(zik)<0.) zik = -zik

    !The plasma dispersion function is calculated (analytical energy integration)
    Zgik = ci * sqrtpi * wofzweid(zik)

    !Now the function is calculated based on the dispersion function calculation
    Aiz = 1. + zik * Zgik !Z1/z
    Biz = 0.5 + zik2 * Aiz  !Z2/z              
    Ciz =  0.75 + zik2 * Biz

    !Bounce average of vpar^2 / (vts^2 * v^2)
    vpar2 = 2.*epsilon(pFFk)*kk/Kg*E2g

!!$    nwge = nwg*Tex(pFFk)
!!$    bbip = CMPLX(cthe(pFFk)*omega2bar/(Kg*fki*nwge*qx(pFFk)*Ro(pFFk)),0.)
!!$    ddip2 = bbip**2 + 4.*zik2
!!$    ddip =SQRT(ddip2)
!!$    Vmip = (- bbip - ddip)/2.
!!$    Vpip = (- bbip + ddip)/2.

    IF (caseflag == 1) THEN
       Fik = 2./fk*((-zik2*fk - 1.5 * Ate(pFFk)+Ane(pFFk))*Aiz + Ate(pFFk)*Biz)
    ELSEIF (caseflag == 2) THEN
       Fik = 2./fk*(-1.5*Aiz+Biz)
    ELSEIF (caseflag == 3) THEN
       Fik = 2./fk * Aiz
    ELSEIF (caseflag == 9) THEN
       Fik = 0. 
    ELSEIF (caseflag == 4) THEN
       Fik = -2.* zik2 * Aiz
    ELSEIF (caseflag == 5) THEN 
       Fik = 2./fk*( (-zik2*fk - 1.5 * Ate(pFFk)+Ane(pFFk))*Biz + Ate(pFFk)*Ciz)
    ELSEIF (caseflag == 6) THEN
       Fik = 2./fk*(-1.5*Biz+Ciz)
    ELSEIF (caseflag == 7) THEN
       Fik = 2./fk * Biz
    ELSEIF (caseflag == 8) THEN
       Fik = -2.* zik2 * Biz
    ELSEIF (caseflag == 10) THEN
       Fik = 0.
    ELSEIF (caseflag == 11) THEN  ! integral for ang momentum
       Fik = 0.
    ENDIF

    !Fik1 corresponds to the higher order response
    ! With rotation do not allow for trapoder1 terms! Anyway not used in Pierre's version   
    Fik1 = 0.

    IF ( (caseflag < 5) .OR. (caseflag == 9) .OR. (caseflag == 11) ) THEN
       FFke_nocollrot = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ELSE
       FFke_nocollrot = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ENDIF
    IF (ABS(FFke_nocollrot) < SQRT(epsD)) FFke_nocollrot=0.
  END FUNCTION FFke_nocollrot


  COMPLEX(KIND=DBL) FUNCTION FFke_nocollrotold(kk,caseflag)
    !---------------------------------------------------------------------
    ! Integrand for trapped electrons when no collisionality and rotation is included
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
    REAL(KIND=DBL), INTENT(IN) :: kk
    INTEGER, INTENT(IN) :: caseflag
    COMPLEX(KIND=DBL) :: Fik, Fik1, zik, fk
    COMPLEX(KIND=DBL) :: bbip, ddip, ddip2, Vmip, Vpip
    COMPLEX(KIND=DBL) :: zik2, Zgik
    COMPLEX(KIND=DBL) :: Aiz, Biz, Ciz
    COMPLEX(KIND=DBL) :: ompFFk
    REAL(KIND=DBL)    :: ya, k2, nwge, nwe
    REAL(KIND=DBL)    :: fki, Eg, E2g, Kg, gau2mshift,vpar2

    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    ya = 1.-k2
    Kg = ceik(ya)
    Eg = ceie(ya)
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2
    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    IF (ABS(fki) < minfki)  fki=SIGN(minfki,fki)

    fk = CMPLX(fki,0)

    ! account for gamma_E=-grad(E_r)/B*R/Cref
    ! nwE means here nwE'=-k_theta*grad(Er)/B=k_theta*gammaE(as def in code)*Cref/R normalized to nwd=-ktheta T(1keV)/(eBR)
!    gau2mshift = REAL(mshift) + 2.*AIMAG(mshift)*AIMAG(mwidth)*REAL(mwidth)/(REAL(mwidth)**2 - AIMAG(mwidth)**2)

    gau2mshift = REAL(mshift) + AIMAG(mshift)*AIMAG(mwidth**2)/REAL(mwidth**2)

    nwE = -gammaE(pFFk)/Rhoeff(pFFk)
    ompFFk =  omFFk-nwE*gau2mshift ! gau2mshift is shift of |gau(x)| 

    ! The frequency is weighted by the vertical frequency drift    
    zik2 = ompFFk * (1./Tex(pFFk))/fk

    ! The square root whose imaginary part is positive is chosen 
    ! for compatibility with the Fried-Conte function
    zik  = SQRT(zik2) 
    IF (AIMAG(zik)<0.) zik = -zik

    !The plasma dispersion function is calculated (analytical energy integration)
    Zgik = ci * sqrtpi * wofzweid(zik)

    !Now the function is calculated based on the dispersion function calculation
    Aiz = 1. + zik * Zgik !Z1/z
    Biz = 0.5 + zik2 * Aiz  !Z2/z              
    Ciz =  0.75 + zik2 * Biz

    !Bounce average of vpar^2 / (vts^2 * v^2)
    vpar2 = 2.*epsilon(pFFk)*kk/Kg*E2g

!!$    nwge = nwg*Tex(pFFk)
!!$    bbip = CMPLX(cthe(pFFk)*omega2bar/(Kg*fki*nwge*qx(pFFk)*Ro(pFFk)),0.)
!!$    ddip2 = bbip**2 + 4.*zik2
!!$    ddip =SQRT(ddip2)
!!$    Vmip = (- bbip - ddip)/2.
!!$    Vpip = (- bbip + ddip)/2.

    IF (caseflag == 1) THEN
       Fik = 2./fk*(  (1.-Mache(pFFk)**2)*( (-zik2*fk - 1.5 * Ate(pFFk)+Ane(pFFk))*Aiz + Ate(pFFk)*Biz) &
            & + (Mache(pFFk)**2*Ate(pFFk) - 2.*Mache(pFFk)*Aue(pFFk) )*Aiz + 4.*Mache(pFFk)*Aue(pFFk)*vpar2*Biz)
    ELSEIF (caseflag == 2) THEN
       Fik = 2./fk*( (1.- Mache(pFFk)**2)*(-1.5*Aiz+Biz) + Mache(pFFk)**2*Aiz)
    ELSEIF (caseflag == 3) THEN
       Fik = 2./fk * Aiz * (1-Mache(pFFk)**2) 
    ELSEIF (caseflag == 9) THEN
       Fik = -4./fk* Mache(pFFk)*Aiz + 8./fk*Mache(pFFk)*vpar2*Biz
    ELSEIF (caseflag == 4) THEN
       Fik = -2.* zik2 * Aiz * (1-Mache(pFFk)**2) 
    ELSEIF (caseflag == 5) THEN 
       Fik = 2./fk*(  (1.-Mache(pFFk)**2)*( (-zik2*fk - 1.5 * Ate(pFFk)+Ane(pFFk))*Biz + Ate(pFFk)*Ciz) &
            & + (Mache(pFFk)**2*Ate(pFFk) - 2.*Mache(pFFk)*Aue(pFFk) )*Biz + 4.*Mache(pFFk)*Aue(pFFk)*vpar2*Ciz)
    ELSEIF (caseflag == 6) THEN
       Fik = 2./fk*( (1.- Mache(pFFk)**2)*(-1.5*Biz+Ciz) + Mache(pFFk)**2*Biz)
    ELSEIF (caseflag == 7) THEN
       Fik = 2./fk * Biz * (1.-Mache(pFFk)**2) 
    ELSEIF (caseflag == 8) THEN
       Fik = -2.* zik2 * Biz*(1.-Mache(pFFk)**2) 
    ELSEIF (caseflag == 10) THEN
       Fik = -4./fk* Mache(pFFk)*Biz+ 8./fk*Mache(pFFk)*vpar2*Ciz
    ELSEIF (caseflag == 11) THEN  ! integral for ang momentum
       Fik = 4./fk * vpar2 * (Biz*(Aue(pFFk)+Mache(pFFk)*(Ane(pFFk)-2.5*Ate(pFFk)-zik2*fk) ) + &
            & Ciz*Mache(pFFk)*Ate(pFFk)   )
    ENDIF

    !Fik1 corresponds to the higher order response
    ! With rotation do not allow for trapoder1 terms! Anyway not used in Pierre's version   
    Fik1 = 0.

    IF ( (caseflag < 5) .OR. (caseflag == 9) .OR. (caseflag == 11) ) THEN
       FFke_nocollrotold = kk * Kg * ft(pFFk) *  Nex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ELSE
       FFke_nocollrotold = kk * Kg * ft(pFFk) *  Nex(pFFk) * Tex(pFFk) * (Fik  * Joe2p + Fik1 * J1e2p) 
    ENDIF
    IF (ABS(FFke_nocollrotold) < SQRT(epsD)) FFke_nocollrotold=0.
  END FUNCTION FFke_nocollrotold


END MODULE trapints
