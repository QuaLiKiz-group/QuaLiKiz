! --------------------------------------------------------------
! PURPOSE: READ ALL INPUT PARAMETERS, CREATE DERIVED QUANTITIES
! --------------------------------------------------------------

MODULE mod_make_io
  USE kind
  USE datcal
  USE datmat
  USE mod_io_management

  IMPLICIT NONE

CONTAINS

  SUBROUTINE make_input(dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, kthetarhosin, & !general param
       & xin, rhoin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & !geometry
       & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & !electrons
       & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & !ions
       & Machtorin, Autorin, Machparin, Auparin, gammaEin, & !rotation
       & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, ETGmultin, collmultin)  !code specific inputs

    ! List of input variables
    INTEGER, INTENT(IN) :: dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, el_typein,verbosein
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
    ALLOCATE(cftrans(dimx,nions,7)); cftrans(:,:,:)=0. !includes 7 transport coefficients, only for ions
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

    Zeffx(:)=0 !Initialize Zeff
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
    ALLOCATE(fonxcirci(nions)); fonxcirci=0; fonxcirce=0
    ALLOCATE(fonxpiegi(nions)); fonxpiegi=0; fonxpiege=0
    ALLOCATE(fonxecirci(nions)); fonxecirci=0; fonxecirce=0
    ALLOCATE(fonxepiegi(nions)); fonxepiegi=0; fonxepiege=0
    ALLOCATE(fonxvcirci(nions)); fonxvcirci=0; fonxvcirce=0
    ALLOCATE(fonxvpiegi(nions)); fonxvpiegi=0; fonxvpiege=0

    ALLOCATE( krmmuITG (dimx) );
    ALLOCATE( krmmuETG (dimx) );

    !Real 2D arrays
    ALLOCATE( kperp2 (dimx, dimn) ); kperp2=0
    ALLOCATE( modewidth (dimx, dimn) ); modewidth=0
    ALLOCATE( modeshift (dimx, dimn) ); modeshift=0
    ALLOCATE( modeshift2 (dimx, dimn) ); modeshift2=0
    ALLOCATE( distan (dimx, dimn) ); distan=0
    ALLOCATE( FLRec (dimx, dimn) ); FLRec=0
    ALLOCATE( FLRep (dimx, dimn) ); FLRep=0
    !Real 3D arrays with 3rd dimension equal to number of ions
    ALLOCATE( FLRic (dimx, dimn,nions) ); FLRic=0
    ALLOCATE( FLRip (dimx, dimn,nions) ); FLRip=0
    !Real 3D arrays with 3rd dimension equal to number of solutions searched for
    ALLOCATE( gamma (dimx, dimn, numsols) ); gamma=0
    ALLOCATE( Ladia (dimx, dimn, numsols) ); Ladia=0
    !Complex 2D arrays 
    ALLOCATE( ommax (dimx, dimn) ) ; ommax=0
    ALLOCATE( solflu (dimx, dimn) ); solflu=0
    !DEBUGGING FOR FLUID SOLUTIONS
    ALLOCATE(jon_solflu(dimx,dimn)); jon_solflu=0
    ALLOCATE(jon_modewidth(dimx,dimn)); jon_modewidth=0
    ALLOCATE(jon_modeshift(dimx,dimn)); jon_modeshift=0
    ALLOCATE(cot_solflu(dimx,dimn)); cot_solflu=0
    ALLOCATE(cot_modewidth(dimx,dimn)); cot_modewidth=0
    ALLOCATE(cot_modeshift(dimx,dimn)); cot_modeshift=0
    ALLOCATE(ana_solflu(dimx,dimn)); ana_solflu=0
    ALLOCATE(old_modewidth(dimx,dimn)); old_modewidth=0
    ALLOCATE(old_modeshift(dimx,dimn)); old_modeshift=0


    !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
    ALLOCATE( sol (dimx, dimn, numsols) ); sol=0
    ALLOCATE( fdsol (dimx, dimn, numsols) ); fdsol=0
    ALLOCATE( Lcirce (dimx, dimn, numsols) ); Lcirce=0
    ALLOCATE( Lpiege (dimx, dimn, numsols) ); Lpiege=0
    ALLOCATE( Lecirce (dimx, dimn, numsols) ); Lecirce=0
    ALLOCATE( Lepiege (dimx, dimn, numsols) ); Lepiege=0
    ALLOCATE( Lvcirce (dimx, dimn, numsols) ); Lvcirce=0
    ALLOCATE( Lvpiege (dimx, dimn, numsols) ); Lvpiege=0
    !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
    ALLOCATE( Lcirci (dimx, dimn, nions, numsols) ); Lcirci=0
    ALLOCATE( Lpiegi (dimx, dimn, nions, numsols) ); Lpiegi=0
    ALLOCATE( Lecirci (dimx, dimn, nions, numsols) ); Lecirci=0
    ALLOCATE( Lepiegi (dimx, dimn, nions, numsols) ); Lepiegi=0
    ALLOCATE( Lvcirci (dimx, dimn, nions, numsols) ); Lvcirci=0
    ALLOCATE( Lvpiegi (dimx, dimn, nions, numsols) ); Lvpiegi=0
    !To save memory space, the following output arrays are only produced when specifically requested
    IF (phys_meth /= 0.0) THEN
       !Intermediate 1D complex arrays
       ALLOCATE(fonxcircgti(nions)); fonxcircgti=0; fonxcircgte=0
       ALLOCATE(fonxpieggti(nions)); fonxpieggti=0; fonxpieggte=0
       ALLOCATE(fonxcircgni(nions)); fonxcircgni=0; fonxcircgne=0
       ALLOCATE(fonxpieggni(nions)); fonxpieggni=0; fonxpieggne=0
       ALLOCATE(fonxcircgui(nions)); fonxcircgui=0; fonxcircgue=0
       ALLOCATE(fonxpieggui(nions)); fonxpieggui=0; fonxpieggue=0
       ALLOCATE(fonxcircci(nions)); fonxcircci=0; fonxcircce=0
       ALLOCATE(fonxpiegci(nions)); fonxpiegci=0; fonxpiegce=0
       !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
       ALLOCATE( Lcircgte (dimx, dimn, numsols) ); Lcircgte=0
       ALLOCATE( Lpieggte (dimx, dimn, numsols) ); Lpieggte=0
       ALLOCATE( Lcircgne (dimx, dimn, numsols) ); Lcircgne=0
       ALLOCATE( Lpieggne (dimx, dimn, numsols) ); Lpieggne=0
       ALLOCATE( Lcircgue (dimx, dimn, numsols) ); Lcircgue=0
       ALLOCATE( Lpieggue (dimx, dimn, numsols) ); Lpieggue=0
       ALLOCATE( Lcircce (dimx, dimn, numsols) ); Lcircce=0
       ALLOCATE( Lpiegce (dimx, dimn, numsols) ); Lpiegce=0
       !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
       ALLOCATE( Lcircgti (dimx, dimn, nions, numsols) ); Lcircgti=0
       ALLOCATE( Lpieggti (dimx, dimn, nions, numsols) ); Lpieggti=0
       ALLOCATE( Lcircgni (dimx, dimn, nions, numsols) ); Lcircgni=0
       ALLOCATE( Lpieggni (dimx, dimn, nions, numsols) ); Lpieggni=0
       ALLOCATE( Lcircgui (dimx, dimn, nions, numsols) ); Lcircgui=0
       ALLOCATE( Lpieggui (dimx, dimn, nions, numsols) ); Lpieggui=0
       ALLOCATE( Lcircci (dimx, dimn, nions, numsols) ); Lcircci=0
       ALLOCATE( Lpiegci (dimx, dimn, nions, numsols) ); Lpiegci=0
!!!
       IF (phys_meth == 2) THEN
          !Intermediate 1D complex arrays
          ALLOCATE(fonxecircgti(nions)); fonxecircgti=0; fonxecircgte=0
          ALLOCATE(fonxepieggti(nions)); fonxepieggti=0; fonxepieggte=0
          ALLOCATE(fonxecircgni(nions)); fonxecircgni=0; fonxecircgne=0
          ALLOCATE(fonxepieggni(nions)); fonxepieggni=0; fonxepieggne=0
          ALLOCATE(fonxecircgui(nions)); fonxecircgui=0; fonxecircgue=0
          ALLOCATE(fonxepieggui(nions)); fonxepieggui=0; fonxepieggue=0
          ALLOCATE(fonxecircci(nions)); fonxecircci=0; fonxecircce=0
          ALLOCATE(fonxepiegci(nions)); fonxepiegci=0; fonxepiegce=0
          !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
          ALLOCATE( Lecircgte (dimx, dimn, numsols) ); Lecircgte=0
          ALLOCATE( Lepieggte (dimx, dimn, numsols) ); Lepieggte=0
          ALLOCATE( Lecircgne (dimx, dimn, numsols) ); Lecircgne=0
          ALLOCATE( Lepieggne (dimx, dimn, numsols) ); Lepieggne=0
          ALLOCATE( Lecircgue (dimx, dimn, numsols) ); Lecircgue=0
          ALLOCATE( Lepieggue (dimx, dimn, numsols) ); Lepieggue=0
          ALLOCATE( Lecircce (dimx, dimn, numsols) ); Lecircce=0
          ALLOCATE( Lepiegce (dimx, dimn, numsols) ); Lepiegce=0
          !Complex 4D arrays with 3rd dimension equal to nions and 4th to numsols
          ALLOCATE( Lecircgti (dimx, dimn, nions, numsols) ); Lecircgti=0
          ALLOCATE( Lepieggti (dimx, dimn, nions, numsols) ); Lepieggti=0
          ALLOCATE( Lecircgni (dimx, dimn, nions, numsols) ); Lecircgni=0
          ALLOCATE( Lepieggni (dimx, dimn, nions, numsols) ); Lepieggni=0
          ALLOCATE( Lecircgui (dimx, dimn, nions, numsols) ); Lecircgui=0
          ALLOCATE( Lepieggui (dimx, dimn, nions, numsols) ); Lepieggui=0
          ALLOCATE( Lecircci (dimx, dimn, nions, numsols) ); Lecircci=0
          ALLOCATE( Lepiegci (dimx, dimn, nions, numsols) ); Lepiegci=0
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
    DEALLOCATE( Lvpiege )
    DEALLOCATE( Lvcirce )
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
       DEALLOCATE( Lcircgue )
       DEALLOCATE( Lpieggue )
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
          DEALLOCATE( Lecircgue )
          DEALLOCATE( Lepieggue )
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
          fonxcircgue=(0.,0.)
          fonxpieggue=(0.,0.)
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
             fonxecircgue=(0.,0.)
             fonxepieggue=(0.,0.)
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
       fonxvcirce=(0.,0.)
       fonxvpiege=(0.,0.)
       fonxvcirci(:)=(0.,0.)
       fonxvpiegi(:)=(0.,0.)
    ELSE
       Ladia(p,nu,j) = fonxad
       Lcirce(p,nu,j) = fonxcirce
       Lpiege(p,nu,j) = fonxpiege
       Lcirci(p,nu,:,j) = fonxcirci(:)
       Lpiegi(p,nu,:,j) = fonxpiegi(:)
       ! to save memory space, this condition has been introduced, Ale 10/08
       IF (phys_meth .NE. 0.0) THEN
          Lcircgte(p,nu,j) = fonxcircgte
          Lpieggte(p,nu,j) = fonxpieggte
          Lcircgti(p,nu,:,j) = fonxcircgti(:)
          Lpieggti(p,nu,:,j) = fonxpieggti(:)
          Lcircgne(p,nu,j) = fonxcircgne
          Lpieggne(p,nu,j) = fonxpieggne
          Lcircgue(p,nu,j) = fonxcircgue
          Lpieggue(p,nu,j) = fonxpieggue
          Lcircgni(p,nu,:,j) = fonxcircgni(:)
          Lpieggni(p,nu,:,j) = fonxpieggni(:)
          Lcircgui(p,nu,:,j) = fonxcircgui(:)
          Lpieggui(p,nu,:,j) = fonxpieggui(:)
          Lcircce(p,nu,j) = fonxcircce
          Lpiegce(p,nu,j) = fonxpiegce
          Lcircci(p,nu,:,j) = fonxcircci(:)
          Lpiegci(p,nu,:,j) = fonxpiegci(:)
          IF (phys_meth == 2) THEN
             Lecircgte(p,nu,j) = fonxecircgte
             Lepieggte(p,nu,j) = fonxepieggte
             Lecircgti(p,nu,:,j) = fonxecircgti(:)
             Lepieggti(p,nu,:,j) = fonxepieggti(:)
             Lecircgne(p,nu,j) = fonxecircgne
             Lepieggne(p,nu,j) = fonxepieggne
             Lecircgue(p,nu,j) = fonxecircgue
             Lepieggue(p,nu,j) = fonxepieggue
             Lecircgni(p,nu,:,j) = fonxecircgni(:)
             Lepieggni(p,nu,:,j) = fonxepieggni(:)
             Lecircgui(p,nu,:,j) = fonxecircgui(:)
             Lepieggui(p,nu,:,j) = fonxepieggui(:)
             Lecircce(p,nu,j) = fonxecircce
             Lepiegce(p,nu,j) = fonxepiegce
             Lecircci(p,nu,:,j) = fonxecircci(:)
             Lepiegci(p,nu,:,j) = fonxepiegci(:)
          ENDIF
       ENDIF
       Lecirce(p,nu,j) = fonxecirce
       Lepiege(p,nu,j) = fonxepiege
       Lecirci(p,nu,:,j) = fonxecirci(:)
       Lepiegi(p,nu,:,j) = fonxepiegi(:)
       Lvcirce(p,nu,j) = fonxvcirce
       Lvpiege(p,nu,j) = fonxvpiege
       Lvcirci(p,nu,:,j) = fonxvcirci(:)
       Lvpiegi(p,nu,:,j) = fonxvpiegi(:)
    ENDIF

  END SUBROUTINE save_qlfunc

  SUBROUTINE collectarrays()
    ! Collect all output into rank 0
    INTEGER :: ierr,myrank, i, nproc

    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: distantmp,FLRectmp,FLReptmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: gammatmp, Ladiatmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions) :: FLRiptmp, FLRictmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,0:nions,0:9) :: ecoefsgautmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: modewidthtmp, modeshifttmp,jon_modewidthtmp, jon_modeshifttmp,cot_modewidthtmp, cot_modeshifttmp,old_modewidthtmp, old_modeshifttmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: ommaxtmp, solflutmp,jon_solflutmp,ana_solflutmp,cot_solflutmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: soltmp, fdsoltmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lcircetmp, Lpiegetmp, Lecircetmp, Lepiegetmp, Lvcircetmp, Lvpiegetmp, Lcircgtetmp, Lpieggtetmp,  Lcircgnetmp, Lpieggnetmp,  Lcircguetmp, Lpiegguetmp, Lcirccetmp, Lpiegcetmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lecircgtetmp, Lepieggtetmp, Lecircgnetmp, Lepieggnetmp, Lecircguetmp, Lepiegguetmp, Lecirccetmp, Lepiegcetmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lcircitmp, Lpiegitmp, Lecircitmp, Lepiegitmp, Lvcircitmp, Lvpiegitmp, Lcircgtitmp, Lpieggtitmp, Lcircgnitmp, Lpieggnitmp, Lcircguitmp, Lpiegguitmp, Lcirccitmp, Lpiegcitmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lecircgtitmp, Lepieggtitmp, Lecircgnitmp, Lepieggnitmp, Lecircguitmp, Lepiegguitmp, Lecirccitmp, Lepiegcitmp

    CALL mpi_comm_rank(mpi_comm_world,myrank,ierr)
    CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(distan,distantmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(modewidth,modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(modeshift,modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(jon_modewidth,jon_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(jon_modeshift,jon_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(cot_modewidth,cot_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(cot_modeshift,cot_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(old_modewidth,old_modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(old_modeshift,old_modeshifttmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRec,FLRectmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRep,FLReptmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(gamma,gammatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Ladia,Ladiatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRic,FLRictmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRip,FLRiptmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(ommax,ommaxtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(solflu,solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(jon_solflu,jon_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(cot_solflu,cot_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(ana_solflu,ana_solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(sol,soltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(fdsol,fdsoltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lcirce,Lcircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lpiege,Lpiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lecirce,Lecircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lepiege,Lepiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvcirce,Lvcircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvpiege,Lvpiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lcirci,Lcircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lpiegi,Lpiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lecirci,Lecircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lepiegi,Lepiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvcirci,Lvcircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvpiegi,Lvpiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierr)
    CALL MPI_Reduce(ecoefsgau,ecoefsgautmp,dimx*dimn*(nions+1)*10,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)

    IF (phys_meth /= 0) THEN
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgte,Lcircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggte,Lpieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgne,Lcircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggne,Lpieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgue,Lcircguetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggue,Lpiegguetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircce,Lcirccetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpiegce,Lpiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgti,Lcircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggti,Lpieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgni,Lcircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggni,Lpieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgui,Lcircguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggui,Lpiegguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircci,Lcirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpiegci,Lpiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
!!!
       IF (phys_meth == 2) THEN
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgte,Lecircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggte,Lepieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgne,Lecircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggne,Lepieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgue,Lecircguetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggue,Lepiegguetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircce,Lecirccetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepiegce,Lepiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgti,Lecircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggti,Lepieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgni,Lecircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggni,Lepieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgui,Lecircguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggui,Lepiegguitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircci,Lecirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Barrier(mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepiegci,Lepiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       ENDIF
    ENDIF

    IF (myrank == 0) THEN
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
       Lvcirce=Lvcircetmp
       Lvpiege=Lvpiegetmp
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
          Lcircgue=Lcircguetmp
          Lpieggue=Lpiegguetmp
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
             Lecircgue=Lecircguetmp
             Lepieggue=Lepiegguetmp
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
    ENDIF

  END SUBROUTINE collectarrays

END MODULE mod_make_io
