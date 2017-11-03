! --------------------------------------------------------------
! PURPOSE: READ ALL INPUT PARAMETERS, CREATE DERIVED QUANTITIES
! --------------------------------------------------------------

MODULE mod_make_io
  USE kind
  USE datcal
  USE datmat

  IMPLICIT NONE
  INCLUDE 'mpif.h'

CONTAINS

  SUBROUTINE make_input(dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, kthetarhosin, & !general param
       & xin, rhoin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & !geometry
       & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & !electrons
       & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & !ions
       & Machtorin, Autorin, Machparin, Auparin, gammaEin, & !rotation
       & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, ETGmultin, collmultin)  !code specific inputs

    ! List of input variables
    INTEGER, INTENT(IN) :: dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, el_typein,verbosein, separatefluxin
    INTEGER, DIMENSION(:,:), INTENT(IN) :: ion_typein
    REAL(kind=DBL), INTENT(IN) :: R0in
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: xin, rhoin, Roin, Rminin, Boin, kthetarhosin, qxin, smagin, alphaxin
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: Texin, Nexin, Atein, Anein, anisein, danisedrin
    REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) :: Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, Aiin, Ziin
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: Machtorin, Autorin, Machparin, Auparin, gammaEin
    INTEGER, INTENT(IN) :: maxrunsin, maxptsin
    REAL(kind=DBL), INTENT(IN) :: relacc1in, relacc2in, timeoutin, ETGmultin, collmultin

    INTEGER:: p,nu !counter for loop over coordinates. p is radius (or general scan), nu is wavenumber
    INTEGER:: i,ilow,ihi !counter for loops
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Lambe, Nue !physical quantities used for the derivations, then discarded

    ! POPULATE GLOBAL VARIABLES
    dimx = dimxin
    dimn = dimnin
    nions = nionsin 
    numsols = numsolsin
    phys_meth = phys_methin
    coll_flag = coll_flagin
    rot_flag = rot_flagin
    IF (verbosein > 0) verbose = .TRUE.
    IF (verbosein == 0) verbose = .FALSE.
    IF (separatefluxin == 1) separateflux = .TRUE.
    IF (separatefluxin == 0) separateflux = .FALSE.

    el_type = el_typein

    maxruns=maxrunsin
    maxpts=maxptsin
    relacc1=relacc1in
    relacc2=relacc2in
    timeout=timeoutin
    relaccQL1=relacc1in
    relaccQL2=relacc2in
    ETGmult=ETGmultin
    collmult=collmultin

    R0=R0in

    !lenwrk=(ndim+2+1)*(1+maxptsin/(2**ndim+2*ndim*ndim+2*ndim+1)) !Set array size for 2D integration routine
    lenwrk = 1e5

    !Input array allocation
    ALLOCATE(kthetarhos(dimn)); kthetarhos = kthetarhosin
    ALLOCATE(x(dimx)); x = xin
    WHERE(x < 0.05) x = 0.05 !filter points near magnetic axis
    ALLOCATE(rho(dimx)); rho = rhoin
    ALLOCATE(Ro(dimx)); Ro = Roin
    ALLOCATE(Rmin(dimx)); Rmin = Rminin
    ALLOCATE(Bo(dimx)); Bo = Boin
    ALLOCATE(qx(dimx)); qx = qxin
    ALLOCATE(smag(dimx)); smag = smagin
    WHERE(ABS(smag) < 0.1) smag = SIGN(0.1,smag) !filter very low magnetic shear which is not valid under QLK assumptions (locality and strong ballooning)
    WHERE(smag < -0.3) smag = -0.3               !filter very negative magnetic shear which is not presently valid under QLK assumptions (we have no slab modes)
    ALLOCATE(alphax(dimx)); alphax = alphaxin
    WHERE( (smag - alphax) < 0.2 ) alphax = -0.2 + smag !filter high alpha where s-alpha goes into invalid regime under QLK assumptions (we have no slab modes)
    WHERE( alphax < 0.) alphax = 0. ! For negative magnetic shear: make sure pressure gradients are not positive due to above constraint
    ALLOCATE(Machtor(dimx)); Machtor = Machtorin ! Mach number in toroidal direction Mach = v_tor / C_ref, C_ref thermal velocity of D at 1keV sqrt(1keV/m_p)
    ALLOCATE(Autor(dimx)); Autor = Autorin ! Normalized toroidal velocity gradient: - R\nabla U_\tor / C_ref
    ALLOCATE(Machpar(dimx)); Machpar = Machparin ! Mach number in parallel direction Mach = v_par / C_ref, C_ref thermal velocity of D at 1keV sqrt(1keV/m_p)
    ALLOCATE(Aupar(dimx)); Aupar = Auparin ! Normalized parallel velocity gradient: - R\nabla U_\parallel / C_ref 
    ALLOCATE(gammaE(dimx)); gammaE = gammaEin ! normalized gammaE: - R \nablaE_r/B /C_ref (B here is B total not B_\varphi)
    ALLOCATE(Tex(dimx)); Tex = Texin
    ALLOCATE(Nex(dimx)); Nex = Nexin
    ALLOCATE(Ate(dimx)); Ate = Atein
    ALLOCATE(Ane(dimx)); Ane = Anein
    ALLOCATE(anise(dimx)); anise = anisein
    ALLOCATE(danisedr(dimx)); danisedr = danisedrin
    ALLOCATE(Ai(dimx,nions)); Ai = Aiin
    ALLOCATE(Zi(dimx,nions)); Zi = Ziin
    ALLOCATE(Tix(dimx,nions)); Tix = Tixin
    ALLOCATE(ninorm(dimx,nions)); ninorm = ninormin
    ALLOCATE(Ati(dimx,nions)); Ati = Atiin
    ALLOCATE(Ani(dimx,nions)); Ani = Aniin
    ALLOCATE(ion_type(dimx,nions)); ion_type = ion_typein
    ALLOCATE(anis(dimx,1:nions)); anis = anisin
    ALLOCATE(danisdr(dimx,1:nions)); danisdr = danisdrin

    ! Set radial dependent rotation flag
    ALLOCATE(rotflagarray(dimx)); 
    IF (rot_flag == 0) THEN
       rotflagarray(:)=0
    ELSE
       rotflagarray(:)=1
    ENDIF

    IF (rot_flag == 2) THEN !Do not use rotation version (much slower) for rho < 0.4
       WHERE (rho < 0.4) rotflagarray = 0
       ALLOCATE(Machparorig(dimx)); Machparorig = Machpar
       ALLOCATE(Auparorig(dimx)); Auparorig = Aupar
       ALLOCATE(gammaEorig(dimx)); gammaEorig = gammaE
       ALLOCATE(Machiorig(dimx,nions))
       ALLOCATE(Auiorig(dimx,nions))

       ALLOCATE(Machparmod(dimx)); Machparmod=1d-14
       ALLOCATE(Auparmod(dimx)) ; Auparmod=1d-14
       ALLOCATE(gammaEmod(dimx)) ; gammaEmod=1d-14
       ALLOCATE(Machimod(dimx,nions));  Machimod=1d-14
       ALLOCATE(Auimod(dimx,nions)); Auimod=1d-14
       ALLOCATE(filterprof(dimx))

       ihi=dimx
       DO i=1,dimx ! define filterprof
          IF (rho(i) < 0.4) THEN 
             filterprof(i)=1d-14
             ilow=i
          ELSEIF (rho(i) > 0.6) THEN
             filterprof(i)=1.
             IF (i<ihi) ihi=i
          ENDIF
       ENDDO

       filterprof(ilow:ihi) = (/( REAL((i-ilow))/REAL((ihi-ilow))  ,i=ilow,ihi)/)     
       WHERE(filterprof<1d-14) filterprof=1d-14
       gammaEmod=gammaE*filterprof
       Auparmod=Aupar*filterprof
       Machparmod=Machpar*filterprof

    ENDIF

    ! Define all quantities derived in the code from the input arrays

    ! These next three are temporary and are deallocated following the calculations
    ALLOCATE(Lambe(dimx))
    ALLOCATE(Nue(dimx))
    ALLOCATE(epsilon(dimx))

    ! Array allocation
    ALLOCATE(tau(dimx,nions))
    ALLOCATE(Nix(dimx,nions))
    ALLOCATE(Zeffx(dimx))
    ALLOCATE(Nustar(dimx))
    ALLOCATE(qprim(dimx))
    ALLOCATE(Anue(dimx))
    ALLOCATE(wg(dimx))
    ALLOCATE(Ac(dimx))
    ALLOCATE(csou(dimx))
    ALLOCATE(cref(dimx)) ! ref velocity used for input of gamma_E, U_par and \nabla U_par : sqrt(1keV/m_p) i.e. thermal velocity of D at 1keV
    ALLOCATE(cthe(dimx))
    ALLOCATE(cthi(dimx,nions))
    ALLOCATE(omegator(dimx))
    ALLOCATE(domegatordr(dimx))
    ALLOCATE(Rhoe(dimx))
    ALLOCATE(Rhoi(dimx,nions))
    ALLOCATE(de(dimx))
    ALLOCATE(di(dimx,nions))
    ALLOCATE(ktetasn(dimx)) 
    ALLOCATE(rhostar(dimx))
    ALLOCATE(Rhoeff(dimx))
    ALLOCATE(ft(dimx))
    ALLOCATE(fc(dimx))
    ALLOCATE(Machi(dimx,nions))
    ALLOCATE(Machitemp(dimx,nions))
    ALLOCATE(Mache(dimx))
    ALLOCATE(Aui(dimx,nions))
    ALLOCATE(Aue(dimx))
    ALLOCATE(th(ntheta))
    ALLOCATE(phi(dimx,ntheta))
    ALLOCATE(dphidr(dimx,ntheta))
    ALLOCATE(dphidth(dimx,ntheta))
    ALLOCATE(npol(dimx,ntheta,nions))
    ALLOCATE(ecoefs(dimx,0:nions,numecoefs)); ecoefs(:,:,:)=0. !includes electrons
    ALLOCATE(ecoefsgau(dimx,dimn,0:nions,0:9)); ecoefsgau(:,:,:,:)=0. !includes electrons
    ALLOCATE(cftrans(dimx,nions,numicoefs)); cftrans(:,:,:)=0. !includes 6 transport coefficients, only for ions
    ALLOCATE(nepol(dimx,ntheta))
    ALLOCATE(Anipol(dimx,ntheta,nions))
    ALLOCATE(Anepol(dimx,ntheta))
    ALLOCATE(tpernorme(dimx,ntheta))
    ALLOCATE(tpernormi(dimx,ntheta,nions))
    ALLOCATE(tpernormifunc(nions))

    ALLOCATE(coefi(dimx,nions))
    ALLOCATE(ntor (dimx, dimn) )   
    ALLOCATE(mi(dimx,nions))
    ALLOCATE(ETG_flag(dimn))

    !These next ones are calculate at each p inside calcroutines
    ALLOCATE(Athi(nions))
    ALLOCATE(ktetaRhoi(nions))
    ALLOCATE(Joi2(nions))
    ALLOCATE(Jobani2(nions))
    ALLOCATE(Joi2p(nions))
    ALLOCATE(J1i2p(nions))
    ALLOCATE(Joi2c(nions))

    mi(:,:) = Ai(:,:)*mp

    !STYLE NOTE: The (:) is not necessary but useful as a reminder for which 
    !variables are scalars and which are arrays

    ! Some auxiliary definitions
    epsilon(:) = Rmin(:)*x(:)/Ro(:)
    ft(:) = 2.*(2.*epsilon(:))**0.5/pi; !trapped particle fraction
    fc(:) = 1. - ft(:) !passing particle fraction
    qprim(:)=smag(:)*qx(:)/(Rmin(:)*x(:))

    th=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set poloidal angle from [0,pi]

    ! Thermal velocities
    csou(:) = SQRT(qe*Tex(:)*1.d3/mi(:,1)) !Calculated with respect to main ions (assumed index 1)
    cref(:) = SQRT(qe*1.d3/mp) !Cref=sqrt(2x1keV/mD)=sqrt(1keV/mp) used to normalized gammaEin, Machin, Auparin
    cthe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3/me)

    Zeffx(:)=0. !Initialize Zeff
    Ac(:) = Nex(:) !Adiabatic term, electrons
    DO i = 1,nions
       tau(:,i) = Tex(:)/Tix(:,i) !Temperature ratios

       WHERE (ion_type(:,i) == 4) 
          Zeffx(:) = Zeffx(:)+ninorm(:,i)*Zi(:,i)**2 
          ion_type(:,i) = 3 !Set ion type to pure tracer for rest of code. 
       ENDWHERE

       !Set tracer density. Arbitrary 1e-5, in case ninorm was 0 in input. Does not affect physics since in the end we divide fluxes by nz
       WHERE (ion_type(:,i) == 3) ninorm(:,i)=1d-5 

    ENDDO

    IF (nions > 1) THEN
       DO i = 1,dimx !impose quasineutrality due to potential ion_type=3
          IF (nions == 2) THEN
             ninorm(i,1) = (1. - Zi(i,2)*ninorm(i,2)) /Zi(i,1)
          ELSE
             ninorm(i,1) = (1. - SUM(Zi(i,2:nions)*ninorm(i,2:nions))) /Zi(i,1)
          ENDIF
       ENDDO
    ENDIF

    DO i = 1,nions
       Nix(:,i) = ninorm(:,i)*Nex(:) !Ion densities
       coefi(:,i) = Zi(:,i)*Zi(:,i) * Nix(:,i) * tau(:,i) !Ion coefficients throughout equations
       cthi(:,i) = SQRT(2._DBL*qe*Tix(:,i)*1.d3/mi(:,i)) !Thermal velocities
       Rhoi(:,i) = SQRT(2._DBL*qe*Tix(:,i)*1.d3*mi(:,i))/(qe*Zi(:,i)*Bo(:)) !Larmor radii
       di(:,i) = qx(:)/SQRT(2._DBL*(epsilon(:)+eps))*Rhoi(:,i) !Ion banana width

       Machi(:,i) = Machpar(:)*cref(:)/cthi(:,i) !Mach number in parallel direction
       Machitemp = Machi ! Original, unaltered Mach number

       IF (rot_flag == 2) THEN
          Machiorig(:,i) = Machparorig(:)*cref(:)/cthi(:,i) !Mach number in parallel direction
          Machimod(:,i) = Machparmod(:)*cref(:)/cthi(:,i) !Mach number in parallel direction
          WHERE(Ai>20.5) Machiorig=0. !Filter out Mach numbers for heavy impurities due to breaking of Mach ordering. 
          !Split off at Ne (still relatively common in tokamaks and completes 2nd row of periodic table)
          WHERE(Ai>20.5) Machimod=0. 
       ENDIF

       WHERE (ion_type(:,i) .NE. 3) 
          Ac(:) = Ac(:) + Zi(:,i)**2._DBL*Nix(:,i)*tau(:,i) !Rest of adiabatic term
          Zeffx(:) = Zeffx(:)+ninorm(:,i)*Zi(:,i)**2 !Zeff
       ENDWHERE
    ENDDO

    Lambe(:) = 15.2_DBL - 0.5_DBL*LOG(0.1_DBL*Nex(:)) + LOG(Tex(:))  !Coulomb constant and collisionality. Wesson 2nd edition p661-663
    Nue(:) = 1._DBL/(1.09d-3) *Zeffx(:)*Nex(:)*Lambe(:)/(Tex(:))**1.5_DBL*collmult 
    Nustar(:) = Nue(:)*qx(:)*Ro(:)/epsilon(:)**1.5/(SQRT(Tex(:)*1.d3*qe/me))

    ! Collisionality array
    IF (coll_flag .NE. 0.0) THEN
       Anue(:) = Nue(:)/epsilon(:) 
    ELSE
       Anue(:) = 0._DBL
    ENDIF

    !DEBUGGING
!!$    OPEN(unit=700, file="input/Zeffx.dat", action="write", status="replace")
!!$    WRITE(700,'(G15.7)') (Zeffx(i),i=1,dimx) ; CLOSE(700)
!!$    OPEN(unit=700, file="input/Anue.dat", action="write", status="replace")
!!$    WRITE(700,'(G15.7)') (Anue(i),i=1,dimx) ; CLOSE(700)
!!$    OPEN(unit=700, file="input/Ac.dat", action="write", status="replace")
!!$    WRITE(700,'(G15.7)') (Ac(i),i=1,dimx) ; CLOSE(700)


    ! Normalisation factor
    wg(:) = qx(:)*1.d3/(Ro(:)*Bo(:)*Rmin(:)*(x(:)+eps))

    ! Larmor radii
    Rhoe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3*me)/(qe*Bo(:))

    ! Banana widths
    de(:) = qx(:)/SQRT(2._DBL*(epsilon(:)+eps))*Rhoe(:) 

    !Rotation terms
    Mache(:) = Machpar(:)*cref(:)/cthe(:)

    !Set ETG_flag for when kthetarhos > 2 
    WHERE (kthetarhos > ETGk) 
       ETG_flag = .TRUE.
    ELSEWHERE
       ETG_flag = .FALSE.
    ENDWHERE

    Aue(:) = Aupar(:)*cref(:)/cthe(:)      
    DO i = 1,nions
       Aui(:,i) = Aupar(:)*cref(:)/cthi(:,i)
       IF (rot_flag == 2) THEN
          Auiorig(:,i) = Auparorig(:)*cref(:)/cthi(:,i)
          Auimod(:,i) = Auparmod(:)*cref(:)/cthi(:,i)
       ENDIF
    ENDDO



    !   omega=cref(irad)*Mach(irad)/Ro(irad)*SQRT(1+(epsilon(irad)/qx(irad))**2) !angular toroidal velocity of main ions. We assume that all ions rotate with same omega
    !   domegadr = -Aupar(irad)*cref(irad)/Ro(irad)**2*SQRT(1+(epsilon(irad)/qx(irad))**2) !rotation angular velocity 

    omegator=cref(:)*Machtor(:)/Ro(:) !angular toroidal velocity of main ions. Assumed that all ions rotate at same velocity
    domegatordr = -Aupar(:)*cref(:)/Ro(:)**2 !rotation angular velocity 

    ! Final
    ktetasn(:) = qx(:)/(Rmin(:)*(x(:)+eps)) 
    rhostar(:) = 1./Rmin(:)*SQRT(Tex(:)*1.d3*qe/mi(:,1))/(Zi(:,1)*qe*Bo(:)/mi(:,1)) !With respect to main ion
    Rhoeff(:) = SQRT(qe*1.e3*mp)/(qe*Bo(:)) ! This is the deuterium Larmor radius for Ti = 1keV.  used to normalize nwE, gammaE

    ! ntor is the toroidal wave-number grid for each position of the scan
    DO p=1,dimx
       DO nu=1,dimn
          ntor(p,nu) = kthetarhos(nu)*x(p)/(qx(p)*rhostar(p))
       ENDDO
    ENDDO

    !Integer ntor can sometimes cause problems at low x when ntor is very low
    !ntor=ANINT(ntor)

    DEALLOCATE(Lambe)
    DEALLOCATE(Nue)

  END SUBROUTINE make_input

  SUBROUTINE allocate_output()
    !Complex 1D arrays. These ion arrays are not the final output but used in an intermediate step
    ALLOCATE(fonxcirci(nions)); fonxcirci=0.; fonxcirce=0.
    ALLOCATE(fonxpiegi(nions)); fonxpiegi=0.; fonxpiege=0.
    ALLOCATE(fonxecirci(nions)); fonxecirci=0.; fonxecirce=0.
    ALLOCATE(fonxepiegi(nions)); fonxepiegi=0.; fonxepiege=0.
    ALLOCATE(fonxvcirci(nions)); fonxvcirci=0.; 
    ALLOCATE(fonxvpiegi(nions)); fonxvpiegi=0.; 

    ALLOCATE( krmmuITG (dimx) ); krmmuITG=0;
    ALLOCATE( krmmuETG (dimx) ); krmmuETG=0;

    !Real 2D arrays
    ALLOCATE( kperp2 (dimx, dimn) ); kperp2=0.
    ALLOCATE( modewidth (dimx, dimn) ); modewidth=0.
    ALLOCATE( modeshift (dimx, dimn) ); modeshift=0.
    ALLOCATE( modeshift2 (dimx, dimn) ); modeshift2=0.
    ALLOCATE( distan (dimx, dimn) ); distan=0.
    ALLOCATE( FLRec (dimx, dimn) ); FLRec=0.
    ALLOCATE( FLRep (dimx, dimn) ); FLRep=0.
    !Real 3D arrays with 3rd dimension equal to number of ions
    ALLOCATE( FLRic (dimx, dimn,nions) ); FLRic=0.
    ALLOCATE( FLRip (dimx, dimn,nions) ); FLRip=0.
    !Real 3D arrays with 3rd dimension equal to number of solutions searched for
    ALLOCATE( gamma (dimx, dimn, numsols) ); gamma=0.
    ALLOCATE( Ladia (dimx, dimn, numsols) ); Ladia=0.
    !Complex 2D arrays 
    ALLOCATE( ommax (dimx, dimn) ) ; ommax=0.
    ALLOCATE( solflu (dimx, dimn) ); solflu=0.
    !DEBUGGING FOR FLUID SOLUTIONS
    ALLOCATE(jon_solflu(dimx,dimn)); jon_solflu=0.
    ALLOCATE(jon_modewidth(dimx,dimn)); jon_modewidth=0.
    ALLOCATE(jon_modeshift(dimx,dimn)); jon_modeshift=0.
    ALLOCATE(cot_solflu(dimx,dimn)); cot_solflu=0.
    ALLOCATE(cot_modewidth(dimx,dimn)); cot_modewidth=0.
    ALLOCATE(cot_modeshift(dimx,dimn)); cot_modeshift=0.
    ALLOCATE(ana_solflu(dimx,dimn)); ana_solflu=0.
    ALLOCATE(old_modewidth(dimx,dimn)); old_modewidth=0.
    ALLOCATE(old_modeshift(dimx,dimn)); old_modeshift=0.


    !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
    ALLOCATE( sol (dimx, dimn, numsols) ); sol=0.
    ALLOCATE( fdsol (dimx, dimn, numsols) ); fdsol=0.
    ALLOCATE( Lcirce (dimx, dimn, numsols) ); Lcirce=0.
    ALLOCATE( Lpiege (dimx, dimn, numsols) ); Lpiege=0.
    ALLOCATE( Lecirce (dimx, dimn, numsols) ); Lecirce=0.
    ALLOCATE( Lepiege (dimx, dimn, numsols) ); Lepiege=0.
    !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
    ALLOCATE( Lcirci (dimx, dimn, nions, numsols) ); Lcirci=0.
    ALLOCATE( Lpiegi (dimx, dimn, nions, numsols) ); Lpiegi=0.
    ALLOCATE( Lecirci (dimx, dimn, nions, numsols) ); Lecirci=0.
    ALLOCATE( Lepiegi (dimx, dimn, nions, numsols) ); Lepiegi=0.
    ALLOCATE( Lvcirci (dimx, dimn, nions, numsols) ); Lvcirci=0.
    ALLOCATE( Lvpiegi (dimx, dimn, nions, numsols) ); Lvpiegi=0.
    !To save memory space, the following output arrays are only produced when specifically requested
    IF (phys_meth /= 0.0) THEN
       !Intermediate 1D complex arrays
       ALLOCATE(fonxcircgti(nions)); fonxcircgti=0.; fonxcircgte=0.
       ALLOCATE(fonxpieggti(nions)); fonxpieggti=0.; fonxpieggte=0.
       ALLOCATE(fonxcircgni(nions)); fonxcircgni=0.; fonxcircgne=0.
       ALLOCATE(fonxpieggni(nions)); fonxpieggni=0.; fonxpieggne=0.
       ALLOCATE(fonxcircgui(nions)); fonxcircgui=0.; 
       ALLOCATE(fonxpieggui(nions)); fonxpieggui=0.; 
       ALLOCATE(fonxcircci(nions)); fonxcircci=0.; fonxcircce=0.
       ALLOCATE(fonxpiegci(nions)); fonxpiegci=0.; fonxpiegce=0.
       !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
       ALLOCATE( Lcircgte (dimx, dimn, numsols) ); Lcircgte=0.
       ALLOCATE( Lpieggte (dimx, dimn, numsols) ); Lpieggte=0.
       ALLOCATE( Lcircgne (dimx, dimn, numsols) ); Lcircgne=0.
       ALLOCATE( Lpieggne (dimx, dimn, numsols) ); Lpieggne=0.
       ALLOCATE( Lcircce (dimx, dimn, numsols) ); Lcircce=0.
       ALLOCATE( Lpiegce (dimx, dimn, numsols) ); Lpiegce=0.
       !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
       ALLOCATE( Lcircgti (dimx, dimn, nions, numsols) ); Lcircgti=0.
       ALLOCATE( Lpieggti (dimx, dimn, nions, numsols) ); Lpieggti=0.
       ALLOCATE( Lcircgni (dimx, dimn, nions, numsols) ); Lcircgni=0.
       ALLOCATE( Lpieggni (dimx, dimn, nions, numsols) ); Lpieggni=0.
       ALLOCATE( Lcircgui (dimx, dimn, nions, numsols) ); Lcircgui=0.
       ALLOCATE( Lpieggui (dimx, dimn, nions, numsols) ); Lpieggui=0.
       ALLOCATE( Lcircci (dimx, dimn, nions, numsols) ); Lcircci=0.
       ALLOCATE( Lpiegci (dimx, dimn, nions, numsols) ); Lpiegci=0.
!!!
       IF (phys_meth == 2) THEN
          !Intermediate 1D complex arrays
          ALLOCATE(fonxecircgti(nions)); fonxecircgti=0.; fonxecircgte=0.
          ALLOCATE(fonxepieggti(nions)); fonxepieggti=0.; fonxepieggte=0.
          ALLOCATE(fonxecircgni(nions)); fonxecircgni=0.; fonxecircgne=0.
          ALLOCATE(fonxepieggni(nions)); fonxepieggni=0.; fonxepieggne=0.
          ALLOCATE(fonxecircgui(nions)); fonxecircgui=0.; 
          ALLOCATE(fonxepieggui(nions)); fonxepieggui=0.; 
          ALLOCATE(fonxecircci(nions)); fonxecircci=0.; fonxecircce=0.
          ALLOCATE(fonxepiegci(nions)); fonxepiegci=0.; fonxepiegce=0.
          !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
          ALLOCATE( Lecircgte (dimx, dimn, numsols) ); Lecircgte=0.
          ALLOCATE( Lepieggte (dimx, dimn, numsols) ); Lepieggte=0.
          ALLOCATE( Lecircgne (dimx, dimn, numsols) ); Lecircgne=0.
          ALLOCATE( Lepieggne (dimx, dimn, numsols) ); Lepieggne=0.
          ALLOCATE( Lecircce (dimx, dimn, numsols) ); Lecircce=0.
          ALLOCATE( Lepiegce (dimx, dimn, numsols) ); Lepiegce=0.
          !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
          ALLOCATE( Lecircgti (dimx, dimn, nions, numsols) ); Lecircgti=0.
          ALLOCATE( Lepieggti (dimx, dimn, nions, numsols) ); Lepieggti=0.
          ALLOCATE( Lecircgni (dimx, dimn, nions, numsols) ); Lecircgni=0.
          ALLOCATE( Lepieggni (dimx, dimn, nions, numsols) ); Lepieggni=0.
          ALLOCATE( Lecircgui (dimx, dimn, nions, numsols) ); Lecircgui=0.
          ALLOCATE( Lepieggui (dimx, dimn, nions, numsols) ); Lepieggui=0.
          ALLOCATE( Lecircci (dimx, dimn, nions, numsols) ); Lecircci=0.
          ALLOCATE( Lepiegci (dimx, dimn, nions, numsols) ); Lepiegci=0.
       ENDIF
    ENDIF
  END SUBROUTINE allocate_output

  SUBROUTINE deallocate_all()

    DEALLOCATE(x)
    DEALLOCATE(rho)
    DEALLOCATE(Ro)
    DEALLOCATE(Rmin)
    DEALLOCATE(rotflagarray)
    DEALLOCATE(Bo)
    DEALLOCATE(kthetarhos)
    DEALLOCATE(qx)
    DEALLOCATE(smag)
    DEALLOCATE(Tex)
    DEALLOCATE(Tix)
    DEALLOCATE(Nex)
    DEALLOCATE(Nustar)
    DEALLOCATE(anise)
    DEALLOCATE(danisedr)
    DEALLOCATE(anis)
    DEALLOCATE(danisdr) 
    DEALLOCATE(ion_type)
    DEALLOCATE(Ate)
    DEALLOCATE(Ati)
    DEALLOCATE(Ane)
    DEALLOCATE(Ani)
    DEALLOCATE(alphax)
    DEALLOCATE(Machtor)
    DEALLOCATE(Autor)
    DEALLOCATE(gammaE)
    DEALLOCATE(Machpar)
    DEALLOCATE(Aupar)
    DEALLOCATE(ninorm)
    DEALLOCATE(Ai)
    DEALLOCATE(Zi)

    IF (rot_flag == 2) THEN
       DEALLOCATE(gammaEmod,Auparmod,Machparmod,gammaEorig,Auparorig,Machparorig,Machiorig,Auiorig,Machimod,Auimod,filterprof)
    ENDIF

    DEALLOCATE(Machi)
    DEALLOCATE(Machitemp)
    DEALLOCATE(Mache)
    DEALLOCATE(Aue)
    DEALLOCATE(Aui)
    DEALLOCATE(ft)
    DEALLOCATE(fc)
    DEALLOCATE(npol)
    DEALLOCATE(ecoefs)
    DEALLOCATE(ecoefsgau)
    DEALLOCATE(cftrans)
    DEALLOCATE(nepol)
    DEALLOCATE(Anipol)
    DEALLOCATE(Anepol)
    DEALLOCATE(th)
    DEALLOCATE(tpernormi)
    DEALLOCATE(tpernormifunc)
    DEALLOCATE(tpernorme)
    DEALLOCATE(phi)    
    DEALLOCATE(dphidr)    
    DEALLOCATE(dphidth)    
    DEALLOCATE(epsilon)
    DEALLOCATE(qprim)
    DEALLOCATE(mi)
    DEALLOCATE(ETG_flag)
    DEALLOCATE(tau)
    DEALLOCATE(Nix)
    DEALLOCATE(Anue)
    DEALLOCATE(wg)
    DEALLOCATE(Zeffx)  
    DEALLOCATE(Ac)
    DEALLOCATE(csou)
    DEALLOCATE(cref)
    DEALLOCATE(cthe)
    DEALLOCATE(cthi)
    DEALLOCATE(omegator)
    DEALLOCATE(domegatordr)
    DEALLOCATE(Rhoe)
    DEALLOCATE(Rhoi)
    DEALLOCATE(de)
    DEALLOCATE(di)
    DEALLOCATE(ktetasn) 
    DEALLOCATE(rhostar)
    DEALLOCATE(Rhoeff)
    DEALLOCATE( ntor )
    DEALLOCATE( coefi )
    DEALLOCATE( ktetaRhoi )
    DEALLOCATE(Joi2)
    DEALLOCATE(Jobani2)
    DEALLOCATE(Joi2p)
    DEALLOCATE(J1i2p)
    DEALLOCATE(Joi2c)
    DEALLOCATE( krmmuITG)
    DEALLOCATE( krmmuETG)
    DEALLOCATE( kperp2 )
    DEALLOCATE( modewidth )
    DEALLOCATE( modeshift )
    DEALLOCATE( modeshift2 )
    DEALLOCATE( distan )
    DEALLOCATE( Athi )
    DEALLOCATE( FLRic )
    DEALLOCATE( FLRec )
    DEALLOCATE( FLRip )
    DEALLOCATE( FLRep )

    DEALLOCATE( gamma )
    DEALLOCATE( Ladia )

    DEALLOCATE( ommax )
    DEALLOCATE( solflu )
    DEALLOCATE(jon_solflu)
    DEALLOCATE(jon_modewidth)
    DEALLOCATE(jon_modeshift)
    DEALLOCATE(cot_solflu)
    DEALLOCATE(cot_modewidth)
    DEALLOCATE(cot_modeshift)
    DEALLOCATE(ana_solflu)
    DEALLOCATE(old_modewidth)
    DEALLOCATE(old_modeshift)

    DEALLOCATE( sol )
    DEALLOCATE( fdsol )
    DEALLOCATE(fonxcirci)
    DEALLOCATE(fonxpiegi)
    DEALLOCATE(fonxecirci)
    DEALLOCATE(fonxepiegi)
    DEALLOCATE(fonxvcirci)
    DEALLOCATE(fonxvpiegi)

    DEALLOCATE( Lcirce )
    DEALLOCATE( Lpiege )
    DEALLOCATE( Lcirci )
    DEALLOCATE( Lpiegi )
    DEALLOCATE( Lepiege )
    DEALLOCATE( Lecirce )
    DEALLOCATE( Lecirci )
    DEALLOCATE( Lepiegi )
    DEALLOCATE( Lvcirci )
    DEALLOCATE( Lvpiegi )

    IF (phys_meth /= 0) THEN
       DEALLOCATE(fonxcircgti)
       DEALLOCATE(fonxpieggti)
       DEALLOCATE(fonxcircgni)
       DEALLOCATE(fonxpieggni)
       DEALLOCATE(fonxcircgui)
       DEALLOCATE(fonxpieggui)
       DEALLOCATE(fonxcircci)
       DEALLOCATE(fonxpiegci)
       DEALLOCATE( Lcircgte )
       DEALLOCATE( Lpieggte )
       DEALLOCATE( Lcircgti )
       DEALLOCATE( Lpieggti )
       DEALLOCATE( Lcircgne )
       DEALLOCATE( Lpieggne )
       DEALLOCATE( Lcircgni )
       DEALLOCATE( Lpieggni )
       DEALLOCATE( Lcircgui )
       DEALLOCATE( Lpieggui )
       DEALLOCATE( Lcircce )
       DEALLOCATE( Lpiegce )
       DEALLOCATE( Lcircci )
       DEALLOCATE( Lpiegci )
       IF (phys_meth == 2) THEN
          DEALLOCATE(fonxecircgti)
          DEALLOCATE(fonxepieggti)
          DEALLOCATE(fonxecircgni)
          DEALLOCATE(fonxepieggni)
          DEALLOCATE(fonxecircgui)
          DEALLOCATE(fonxepieggui)
          DEALLOCATE(fonxecircci)
          DEALLOCATE(fonxepiegci)
          DEALLOCATE( Lecircgte )
          DEALLOCATE( Lepieggte )
          DEALLOCATE( Lecircgti )
          DEALLOCATE( Lepieggti )
          DEALLOCATE( Lecircgne )
          DEALLOCATE( Lepieggne )
          DEALLOCATE( Lecircgni )
          DEALLOCATE( Lepieggni )
          DEALLOCATE( Lecircgui )
          DEALLOCATE( Lepieggui )
          DEALLOCATE( Lecircce )
          DEALLOCATE( Lepiegce )
          DEALLOCATE( Lecircci )
          DEALLOCATE( Lepiegci )
       ENDIF
    ENDIF

  END SUBROUTINE deallocate_all

  SUBROUTINE save_qlfunc(p,nu,j,issol)
    !Arguments
    LOGICAL, INTENT(IN) :: issol
    INTEGER, INTENT(IN) :: p,nu,j
    !Local variables

    IF (issol .EQV. .FALSE.) THEN
       fonxad=Ac(p)
       fonxcirce=(0.,0.)
       fonxpiege=(0.,0.)
       fonxcirci(:)=(0.,0.)
       fonxpiegi(:)=(0.,0.)
       ! to save memory space, this condition has been introduced, Ale 10/08
       IF (phys_meth .NE. 0.0) THEN
          fonxcircgte=(0.,0.)
          fonxpieggte=(0.,0.)
          fonxcircgti(:)=(0.,0.)
          fonxpieggti(:)=(0.,0.)
          fonxcircgne=(0.,0.)
          fonxpieggne=(0.,0.)
          fonxcircgni(:)=(0.,0.)
          fonxpieggni(:)=(0.,0.)
          fonxcircgui(:)=(0.,0.)
          fonxpieggui(:)=(0.,0.)
          fonxcircce=(0.,0.)
          fonxpiegce=(0.,0.)
          fonxcircci(:)=(0.,0.)
          fonxpiegci(:)=(0.,0.)
          IF (phys_meth == 2) THEN
             fonxecircgte=(0.,0.)
             fonxepieggte=(0.,0.)
             fonxecircgti(:)=(0.,0.)
             fonxepieggti(:)=(0.,0.)
             fonxecircgne=(0.,0.)
             fonxepieggne=(0.,0.)
             fonxecircgni(:)=(0.,0.)
             fonxepieggni(:)=(0.,0.)
             fonxecircgui(:)=(0.,0.)
             fonxepieggui(:)=(0.,0.)
             fonxecircce=(0.,0.)
             fonxepiegce=(0.,0.)
             fonxecircci(:)=(0.,0.)
             fonxepiegci(:)=(0.,0.)
          ENDIF
       ENDIF
       fonxecirce=(0.,0.)
       fonxepiege=(0.,0.)
       fonxecirci(:)=(0.,0.)
       fonxepiegi(:)=(0.,0.)
       fonxvcirci(:)=(0.,0.)
       fonxvpiegi(:)=(0.,0.)
    ELSE
       Ladia(p,nu,j) = fonxad
       Lcirce(p,nu,j) = AIMAG(fonxcirce)
       Lpiege(p,nu,j) = AIMAG(fonxpiege)
       Lcirci(p,nu,:,j) = AIMAG(fonxcirci(:))
       Lpiegi(p,nu,:,j) = AIMAG(fonxpiegi(:))
       ! to save memory space, this condition has been introduced, Ale 10/08
       IF (phys_meth .NE. 0.0) THEN
          Lcircgte(p,nu,j) = AIMAG(fonxcircgte)
          Lpieggte(p,nu,j) = AIMAG(fonxpieggte)
          Lcircgti(p,nu,:,j) = AIMAG(fonxcircgti(:))
          Lpieggti(p,nu,:,j) = AIMAG(fonxpieggti(:))
          Lcircgne(p,nu,j) = AIMAG(fonxcircgne)
          Lpieggne(p,nu,j) = AIMAG(fonxpieggne)
          Lcircgni(p,nu,:,j) = AIMAG(fonxcircgni(:))
          Lpieggni(p,nu,:,j) = AIMAG(fonxpieggni(:))
          Lcircgui(p,nu,:,j) = AIMAG(fonxcircgui(:))
          Lpieggui(p,nu,:,j) = AIMAG(fonxpieggui(:))
          Lcircce(p,nu,j) = AIMAG(fonxcircce)
          Lpiegce(p,nu,j) = AIMAG(fonxpiegce)
          Lcircci(p,nu,:,j) = AIMAG(fonxcircci(:))
          Lpiegci(p,nu,:,j) = AIMAG(fonxpiegci(:))
          IF (phys_meth == 2) THEN
             Lecircgte(p,nu,j) = AIMAG(fonxecircgte)
             Lepieggte(p,nu,j) = AIMAG(fonxepieggte)
             Lecircgti(p,nu,:,j) = AIMAG(fonxecircgti(:))
             Lepieggti(p,nu,:,j) = AIMAG(fonxepieggti(:))
             Lecircgne(p,nu,j) = AIMAG(fonxecircgne)
             Lepieggne(p,nu,j) = AIMAG(fonxepieggne)
             Lecircgni(p,nu,:,j) = AIMAG(fonxecircgni(:))
             Lepieggni(p,nu,:,j) = AIMAG(fonxepieggni(:))
             Lecircgui(p,nu,:,j) = AIMAG(fonxecircgui(:))
             Lepieggui(p,nu,:,j) = AIMAG(fonxepieggui(:))
             Lecircce(p,nu,j) = AIMAG(fonxecircce)
             Lepiegce(p,nu,j) = AIMAG(fonxepiegce)
             Lecircci(p,nu,:,j) = AIMAG(fonxecircci(:))
             Lepiegci(p,nu,:,j) = AIMAG(fonxepiegci(:))
          ENDIF
       ENDIF
       Lecirce(p,nu,j) = AIMAG(fonxecirce)
       Lepiege(p,nu,j) = AIMAG(fonxepiege)
       Lecirci(p,nu,:,j) = AIMAG(fonxecirci(:))
       Lepiegi(p,nu,:,j) = AIMAG(fonxepiegi(:))
       Lvcirci(p,nu,:,j) = AIMAG(fonxvcirci(:))
       Lvpiegi(p,nu,:,j) = AIMAG(fonxvpiegi(:))
    ENDIF

  END SUBROUTINE save_qlfunc

  SUBROUTINE collectarrays()
    ! Collect all output into all cores
    INTEGER :: ierr,myrank, i, nproc

    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: distantmp,FLRectmp,FLReptmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: gammatmp, Ladiatmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions) :: FLRiptmp, FLRictmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,0:nions,0:9) :: ecoefsgautmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: modewidthtmp, modeshifttmp,jon_modewidthtmp, jon_modeshifttmp,cot_modewidthtmp, cot_modeshifttmp,old_modewidthtmp, old_modeshifttmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: ommaxtmp, solflutmp,jon_solflutmp,ana_solflutmp,cot_solflutmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: soltmp, fdsoltmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lcircetmp, Lpiegetmp, Lecircetmp, Lepiegetmp, Lcircgtetmp, Lpieggtetmp,  Lcircgnetmp, Lpieggnetmp,  Lcirccetmp, Lpiegcetmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lecircgtetmp, Lepieggtetmp, Lecircgnetmp, Lepieggnetmp, Lecirccetmp, Lepiegcetmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lcircitmp, Lpiegitmp, Lecircitmp, Lepiegitmp, Lvcircitmp, Lvpiegitmp, Lcircgtitmp, Lpieggtitmp, Lcircgnitmp, Lpieggnitmp, Lcircguitmp, Lpiegguitmp, Lcirccitmp, Lpiegcitmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lecircgtitmp, Lepieggtitmp, Lecircgnitmp, Lepieggnitmp, Lecircguitmp, Lepiegguitmp, Lecirccitmp, Lepiegcitmp

    CALL mpi_comm_rank(mpi_comm_world,myrank,ierr)
    CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

    CALL MPI_AllReduce(distan,distantmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(modewidth,modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(modeshift,modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(jon_modewidth,jon_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(jon_modeshift,jon_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(cot_modewidth,cot_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(cot_modeshift,cot_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(old_modewidth,old_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(old_modeshift,old_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(FLRec,FLRectmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(FLRep,FLReptmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(gamma,gammatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ladia,Ladiatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(FLRic,FLRictmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(FLRip,FLRiptmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ommax,ommaxtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(solflu,solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(jon_solflu,jon_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(cot_solflu,cot_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ana_solflu,ana_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(sol,soltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(fdsol,fdsoltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(Lcirce,Lcircetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lpiege,Lpiegetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lecirce,Lecircetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lepiege,Lepiegetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lcirci,Lcircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lpiegi,Lpiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lecirci,Lecircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lepiegi,Lepiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lvcirci,Lvcircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Lvpiegi,Lvpiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(ecoefsgau,ecoefsgautmp,dimx*dimn*(nions+1)*10,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    IF (phys_meth /= 0) THEN
       CALL MPI_AllReduce(Lcircgte,Lcircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpieggte,Lpieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lcircgne,Lcircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpieggne,Lpieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lcircce,Lcirccetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpiegce,Lpiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

       CALL MPI_AllReduce(Lcircgti,Lcircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpieggti,Lpieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lcircgni,Lcircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpieggni,Lpieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lcircgui,Lcircguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpieggui,Lpiegguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lcircci,Lcirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(Lpiegci,Lpiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
!!!
       IF (phys_meth == 2) THEN
          CALL MPI_AllReduce(Lecircgte,Lecircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepieggte,Lepieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lecircgne,Lecircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepieggne,Lepieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lecircce,Lecirccetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepiegce,Lepiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

          CALL MPI_AllReduce(Lecircgti,Lecircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepieggti,Lepieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lecircgni,Lecircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepieggni,Lepieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lecircgui,Lecircguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepieggui,Lepiegguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lecircci,Lecirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(Lepiegci,Lepiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       ENDIF
    ENDIF

    modewidth=modewidthtmp
    modeshift=modeshifttmp
    jon_modewidth=jon_modewidthtmp
    jon_modeshift=jon_modeshifttmp
    cot_modewidth=cot_modewidthtmp
    cot_modeshift=cot_modeshifttmp
    old_modewidth=old_modewidthtmp
    old_modeshift=old_modeshifttmp

    distan=distantmp
    FLRec=FLRectmp
    FLRep=FLReptmp
    gamma=gammatmp
    Ladia=Ladiatmp
    FLRip=FLRiptmp
    FLRic=FLRictmp
    ommax=ommaxtmp
    solflu=solflutmp
    jon_solflu=jon_solflutmp
    cot_solflu=cot_solflutmp
    ana_solflu=ana_solflutmp
    sol=soltmp
    fdsol=fdsoltmp
    Lcirce=Lcircetmp
    Lpiege=Lpiegetmp
    Lecirce=Lecircetmp
    Lepiege=Lepiegetmp
    Lcirci=Lcircitmp
    Lpiegi=Lpiegitmp
    Lecirci=Lecircitmp
    Lepiegi=Lepiegitmp
    Lvcirci=Lvcircitmp
    Lvpiegi=Lvpiegitmp
    ecoefsgau=ecoefsgautmp

    IF (phys_meth /= 0) THEN
       Lcircgte=Lcircgtetmp
       Lpieggte=Lpieggtetmp
       Lcircgne=Lcircgnetmp
       Lpieggne=Lpieggnetmp
       Lcircce=Lcirccetmp
       Lpiegce=Lpiegcetmp
       Lcircgti=Lcircgtitmp
       Lpieggti=Lpieggtitmp
       Lcircgni=Lcircgnitmp
       Lpieggni=Lpieggnitmp
       Lcircgui=Lcircguitmp
       Lpieggui=Lpiegguitmp
       Lcircci=Lcirccitmp
       Lpiegci=Lpiegcitmp   
       IF (phys_meth == 2) THEN
          Lecircgte=Lecircgtetmp
          Lepieggte=Lepieggtetmp
          Lecircgne=Lecircgnetmp
          Lepieggne=Lepieggnetmp
          Lecircce=Lecirccetmp
          Lepiegce=Lepiegcetmp
          Lecircgti=Lecircgtitmp
          Lepieggti=Lepieggtitmp
          Lecircgni=Lecircgnitmp
          Lepieggni=Lepieggnitmp
          Lecircgui=Lecircguitmp
          Lepieggui=Lepiegguitmp
          Lecircci=Lecirccitmp
          Lepiegci=Lepiegcitmp   
       ENDIF
    ENDIF

  END SUBROUTINE collectarrays

  SUBROUTINE reduceoutput()
    ! collect all output into all cores for parallel writing
    INTEGER :: ierr,myrank, i, nproc

    REAL(KIND=DBL), DIMENSION(dimx) :: epf_SItmp, eef_SItmp, epf_GBtmp, eef_GBtmp, modeflagtmp
    REAL(KIND=DBL), DIMENSION(dimx) :: krmmuITGtmp, krmmuETGtmp
    REAL(KIND=DBL), DIMENSION(dimx) :: dfe_SItmp, vte_SItmp, vce_SItmp, dfe_GBtmp, vte_GBtmp, vce_GBtmp, cketmp
    REAL(KIND=DBL), DIMENSION(dimx) :: chiee_SItmp, vene_SItmp, vece_SItmp, chiee_GBtmp, vene_GBtmp, vece_GBtmp, ceketmp
    REAL(KIND=DBL), DIMENSION(dimx) :: chieeETG_SItmp, veneETG_SItmp, veceETG_SItmp, chieeETG_GBtmp, veneETG_GBtmp, veceETG_GBtmp
    REAL(KIND=DBL), DIMENSION(dimx,nions) :: ipf_SItmp, ief_SItmp, ivf_SItmp, ipf_GBtmp, ief_GBtmp, ivf_GBtmp
    REAL(KIND=DBL), DIMENSION(dimx,nions) :: dfi_SItmp, vti_SItmp, vci_SItmp, vri_SItmp,dfi_GBtmp, vti_GBtmp, vci_GBtmp, vri_GBtmp, ckitmp
    REAL(KIND=DBL), DIMENSION(dimx,nions) :: chiei_SItmp, veni_SItmp, veci_SItmp, veri_SItmp, chiei_GBtmp, veni_GBtmp, veci_GBtmp, veri_GBtmp, cekitmp

    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: solflu_SItmp, solflu_GBtmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: epf_cmtmp, eef_cmtmp, kperp2tmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions) :: ipf_cmtmp, ief_cmtmp, ivf_cmtmp

    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: gam_SItmp, gam_GBtmp, ome_SItmp, ome_GBtmp
    REAL(KIND=DBL), DIMENSION(dimx,ntheta) :: phitmp
    REAL(KIND=DBL), DIMENSION(dimx,ntheta,nions) :: npoltmp
    REAL(KIND=DBL), DIMENSION(dimx,0:nions,numecoefs) :: ecoefstmp
    REAL(KIND=DBL), DIMENSION(dimx,nions,numicoefs) :: cftranstmp

    CALL mpi_comm_rank(mpi_comm_world,myrank,ierr)
    CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

    CALL MPI_AllReduce(epf_SI,epf_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(eef_SI,eef_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(epf_GB,epf_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(eef_GB,eef_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(modeflag,modeflagtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(krmmuITG,krmmuITGtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(krmmuETG,krmmuETGtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(epf_cm,epf_cmtmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(eef_cm,eef_cmtmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(kperp2,kperp2tmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(solflu_SI,solflu_SItmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(solflu_GB,solflu_GBtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(ipf_SI,ipf_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ief_SI,ief_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ivf_SI,ivf_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ipf_GB,ipf_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ivf_GB,ivf_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ief_GB,ief_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(ipf_cm,ipf_cmtmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ief_cm,ief_cmtmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ivf_cm,ivf_cmtmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(gam_SI,gam_SItmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(gam_GB,gam_GBtmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ome_SI,ome_SItmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ome_GB,ome_GBtmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(phi,phitmp,dimx*ntheta,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(npol,npoltmp,dimx*ntheta*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ecoefs,ecoefstmp,dimx*(1+nions)*numecoefs,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(cftrans,cftranstmp,dimx*nions*numicoefs,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

    IF (phys_meth /= 0.0) THEN
       CALL MPI_AllReduce(dfe_SI,dfe_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vte_SI,vte_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vce_SI,vce_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(dfe_GB,dfe_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vte_GB,vte_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vce_GB,vce_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

       CALL MPI_AllReduce(cke,cketmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(dfi_SI,dfi_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vti_SI,vti_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vri_SI,vri_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vci_SI,vci_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(dfi_GB,dfi_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vti_GB,vti_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vri_GB,vri_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(vci_GB,vci_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(cki,ckitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)

       IF (phys_meth == 2) THEN
          CALL MPI_AllReduce(chiee_SI,chiee_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(vene_SI,vene_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(vece_SI,vece_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(ceke,ceketmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(chiei_SI,chiei_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veni_SI,veni_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veri_SI,veri_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veci_SI,veci_SItmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(ceki,cekitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(chiee_GB,chiee_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(vene_GB,vene_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(vece_GB,vece_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(chiei_GB,chiei_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veni_GB,veni_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veri_GB,veri_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          CALL MPI_AllReduce(veci_GB,veci_GBtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
          IF (separateflux == 1) THEN
             CALL MPI_AllReduce(chieeETG_SI,chieeETG_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
             CALL MPI_AllReduce(veneETG_SI,veneETG_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
             CALL MPI_AllReduce(veceETG_SI,veceETG_SItmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
             CALL MPI_AllReduce(chieeETG_GB,chieeETG_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
             CALL MPI_AllReduce(veneETG_GB,veneETG_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
             CALL MPI_AllReduce(veceETG_GB,veceETG_GBtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)          
          ENDIF
       ENDIF
    ENDIF

    epf_SI=epf_SItmp
    eef_SI=eef_SItmp
    epf_GB=epf_GBtmp
    eef_GB=eef_GBtmp
    modeflag=modeflagtmp
    krmmuITG=krmmuITGtmp
    krmmuETG=krmmuETGtmp

    epf_cm=epf_cmtmp
    eef_cm=eef_cmtmp

    kperp2=kperp2tmp
    solflu_SI=solflu_SItmp
    solflu_GB=solflu_GBtmp

    ipf_SI=ipf_SItmp
    ief_SI=ief_SItmp
    ivf_SI=ivf_SItmp
    ipf_GB=ipf_GBtmp
    ivf_GB=ivf_GBtmp
    ief_GB=ief_GBtmp

    ipf_cm=ipf_cmtmp
    ief_cm=ief_cmtmp
    ivf_cm=ivf_cmtmp

    gam_SI=gam_SItmp
    gam_GB=gam_GBtmp
    ome_SI=ome_SItmp
    ome_GB=ome_GBtmp

    phi=phitmp
    npol=npoltmp
    ecoefs=ecoefstmp
    cftrans=cftranstmp

    IF (phys_meth /= 0.0) THEN

       dfe_SI=dfe_SItmp
       vte_SI=vte_SItmp
       vce_SI=vce_SItmp
       cke=cketmp
       dfi_SI=dfi_SItmp
       vti_SI=vti_SItmp
       vri_SI=vri_SItmp
       vci_SI=vci_SItmp
       cki=ckitmp

       dfe_GB=dfe_GBtmp
       vte_GB=vte_GBtmp
       vce_GB=vce_GBtmp

       dfi_GB=dfi_GBtmp
       vti_GB=vti_GBtmp
       vri_GB=vri_GBtmp
       vci_GB=vci_GBtmp

       IF (phys_meth == 2) THEN

          chiee_SI=chiee_SItmp
          vene_SI=vene_SItmp
          vece_SI=vece_SItmp
          ceke=ceketmp
          chiei_SI=chiei_SItmp
          veni_SI=veni_SItmp
          veri_SI=veri_SItmp
          veci_SI=veci_SItmp
          ceki=cekitmp
          chiee_GB=chiee_GBtmp
          vene_GB=vene_GBtmp
          vece_GB=vece_GBtmp
          chiei_GB=chiei_GBtmp
          veni_GB=veni_GBtmp
          veri_GB=veri_GBtmp
          veci_GB=veci_GBtmp
          IF (separateflux == 1) THEN
             chieeETG_SI=chieeETG_SItmp
             veneETG_SI=veneETG_SItmp
             veceETG_SI=veceETG_SItmp
             chieeETG_GB=chieeETG_GBtmp
             veneETG_GB=veneETG_GBtmp
             veceETG_GB=veceETG_GBtmp
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE reduceoutput


END MODULE mod_make_io
