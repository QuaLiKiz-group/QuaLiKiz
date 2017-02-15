PROGRAM qlk_makeflux
  ! Standalone saturation rule program.
  ! If output from qlk_standalone already exists in output/primitive, then new fluxes can be constructed
  ! using this program without the computationally intensive solution of the dispersion relation. 

  ! USEFUL WHEN NEW POSTPROCESSING IS WISHED TO TAKE PLACE:
  ! Examples: 1) new saturation rule. 2) new values for trace impurity transport.
  USE kind
  USE datmat
  USE mod_io_management
  USE datcal
  USE asymmetry
  USE mod_saturation

  IMPLICIT NONE

  !Timing variables
  REAL(kind=DBL) :: time1,time2,timetot
  INTEGER :: p,nu

  ! Begin time measurement
  CALL CPU_TIME(time1)

  ! Read input and initialize all arrays
  WRITE(stdout,"(A)") '*** Reading input data'
  CALL data_init()
  CALL make_input()

  WRITE(stdout,*)     
  WRITE(stdout,"(A)") '*** Calculating poloidal asymmetry terms'
  CALL calcphi() !calculate poloidal density asymmetries due to rotation and temp anisotropy
  CALL calcdphi() !calculate radial and poloidal gradient of 2D electrostatic potential
  CALL calccoefs() !calculate e# LFS to FSA transport coefficients, FSA R/Ln, FSA n, asym factor, and also 2D R/Ln for all ion species

  WRITE(stdout,*)     
  WRITE(stdout,"(A)") '*** Carrying out saturation rules'
  CALL allocate_endoutput()

  ! Calculate flux-surface-averaging of poloidal asymmetry coefficients including eigenmode

  DO p = 1,dimx
     DO nu=1,dimn
        CALL makeecoefsgau(p,nu)
     ENDDO
  ENDDO

  CALL saturation(0)
  CALL outputascii()

  CALL CPU_TIME(time2)

  timetot = time2-time1

  WRITE(*,*)
  WRITE(stdout,"(A,I0,A)") 'Hurrah! Job completed! Total time = ',INT(timetot*1d3),' ms'  !final write
  CALL deallocate_endoutput()
  !Deallocating all

  CALL deallocate_all()

CONTAINS

  SUBROUTINE make_input()

    INTEGER:: i,p !counters for loops
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Lambe, Nue !physical quantities used for the derivations, then discarded

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
    ALLOCATE(Rhoe(dimx))
    ALLOCATE(Rhoi(dimx,nions))
    ALLOCATE(de(dimx))
    ALLOCATE(di(dimx,nions))
    ALLOCATE(ktetasn(dimx)) 
    ALLOCATE(rhostar(dimx))
    ALLOCATE(ft(dimx))
    ALLOCATE(fc(dimx))
    ALLOCATE(Machi(dimx,nions))
    ALLOCATE(Mache(dimx))
    ALLOCATE(Aui(dimx,nions))
    ALLOCATE(Aue(dimx))

    ALLOCATE(th(ntheta))

    ALLOCATE(phi(dimx,ntheta))
    ALLOCATE(dphidr(dimx,ntheta))
    ALLOCATE(dphidth(dimx,ntheta))
    ALLOCATE(npol(dimx,ntheta,nions))
    ALLOCATE(ecoefs(dimx,0:nions,numecoefs)) !includes electrons
    ALLOCATE(ecoefsgau(dimx,dimn,0:nions,0:9)) !includes electrons
    ALLOCATE(cftrans(dimx,nions,6)) !includes 6 transport coefficients, only for ions
    ALLOCATE(nepol(dimx,ntheta))
    ALLOCATE(Anipol(dimx,ntheta,nions))
    ALLOCATE(Anepol(dimx,ntheta))
    ALLOCATE(tpernorme(dimx,ntheta))
    ALLOCATE(tpernormi(dimx,ntheta,nions))
    ALLOCATE(tpernormifunc(nions))

    ALLOCATE(coefi(dimx,nions))
    ALLOCATE(mi(dimx,nions))

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
    cref(:) = SQRT(qe*1.d3/mp) !Cref=sqrt(2x1keV/mD)=sqrt(1keV/mp) used to normalized gammaEin, Machtorin, Machparin, Auparin, Autorin
    cthe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3/me)

    Zeffx(:)=0 !Initialize Zeff
    Ac(:) = Nex(:) !Adiabatic term, electrons
    DO i = 1,nions
       tau(:,i) = Tex(:)/Tix(:,i) !Temperature ratios
       WHERE (ion_type(:,i) == 4) !ion type is tracer, but whatever ninorm is there should still be included in Zeff. 
          Zeffx(:) = Zeffx(:)+ninorm(:,i)*Zi(:,i)**2 
          ion_type(:,i) = 3 !Set ion type to pure tracer for rest of code. 
       ENDWHERE

       !Set tracer density. Arbitrary 1e-5, in case ninorm was 0 in input. Does not affect physics since in the end we divide fluxes by nz
       WHERE (ion_type(:,i) == 3) ninorm(:,i) = 1d-5 

       Nix(:,i) = ninorm(:,i)*Nex(:) !Ion densities
       coefi(:,i) = Zi(:,i)*Zi(:,i) * Nix(:,i) * tau(:,i) !Ion coefficients throughout equations
       cthi(:,i) = SQRT(2._DBL*qe*Tix(:,i)*1.d3/mi(:,i)) !Thermal velocities
       Rhoi(:,i) = SQRT(2._DBL*qe*Tix(:,i)*1.d3*mi(:,i))/(qe*Zi(:,i)*Bo(:)) !Larmor radii
       di(:,i) = qx(:)/SQRT(2._DBL*(epsilon(:)+eps))*Rhoi(:,i) !Ion banana width
       Machi(:,i) = Machpar(:)*cref(:)/cthi(:,i) !Mach number in parallel direction
       WHERE (ion_type(:,i) .NE. 3)
          Ac(:) = Ac(:) + Zi(:,i)**2._DBL*Nix(:,i)*tau(:,i) !Rest of adiabatic term
          Zeffx(:) = Zeffx(:)+ninorm(:,i)*Zi(:,i)**2 !Zeff
       ENDWHERE
    ENDDO

    Lambe(:) = 15.2_DBL - 0.5_DBL*LOG(0.1_DBL*Nex(:)) + LOG(Tex(:))  !Coulomb constant
    Nue(:) = 1._DBL/(1.09d-3) *Zeffx(:)*Nex(:)*Lambe(:)/(Tex(:))**1.5_DBL !Collisionality

    ! Collisionality array
    IF (coll_flag .NE. 0.0) THEN
       Anue(:) = Nue(:)/epsilon(:)
    ELSE
       Anue(:) = 0._DBL
    ENDIF

    ! Normalisation factor
    wg(:) = qx(:)*1.d3/(Ro(:)*Bo(:)*Rmin(:)*(x(:)+eps))

    ! Larmor radii
    Rhoe(:) = SQRT(2._DBL*qe*Tex(:)*1.d3*me)/(qe*Bo(:))

    ! Banana widths
    de(:) = qx(:)/SQRT(2._DBL*(epsilon(:)+eps))*Rhoe(:)

    !Rotation terms
    Mache(:) = Machpar(:)*cref(:)/cthe(:)

    Aue(:) = Aupar(:)*cref(:)/cthe(:)      
    DO i = 1,nions
       Aui(:,i) = Aupar(:)*cref(:)/cthi(:,i)
    ENDDO

    !   omega=cref(irad)*Mach(irad)/Ro(irad)*SQRT(1+(epsilon(irad)/qx(irad))**2) !angular toroidal velocity of main ions. We assume that all ions rotate with same omega
    !   domegadr = -Aupar(irad)*cref(irad)/Ro(irad)**2*SQRT(1+(epsilon(irad)/qx(irad))**2) !rotation angular velocity 

    omegator=cref(irad)*Machtor(irad)/Ro(irad) !angular toroidal velocity of main ions. Assumed that all ions rotate at same velocity
    domegatordr = -Aupar(irad)*cref(irad)/Ro(irad)**2 !rotation angular velocity 

    ! Final
    ktetasn(:) = qx(:)/(Rmin(:)*(x(:)+eps)) 
    rhostar(:) = 1./Rmin(:)*SQRT(Tex(:)*1.d3*qe/mi(:,1))/(Zi(:,1)*qe*Bo(:)/mi(:,1)) !With respect to main ion

    DEALLOCATE(Lambe)
    DEALLOCATE(Nue)

  END SUBROUTINE make_input

  SUBROUTINE data_init()
    !Read data, allocate input and output arrays
    INTEGER, PARAMETER :: ktype = 1 ! BINARY FILES
    INTEGER :: kc
    INTEGER :: i,j,k,l,myunit=700
    REAL(kind=DBL) :: dummy !dummy variable for obtaining input. Must be real for readvar
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: rdummyn, idummyn !dummy variables for obtaining input
    REAL(kind=DBL), DIMENSION(:,:,:), ALLOCATABLE :: rdummynumsols, idummynumsols !dummy variables for obtaining input
    REAL(kind=DBL), DIMENSION(:,:,:,:), ALLOCATABLE :: rdummynions, idummynions !dummy variables for obtaining input
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: ion_typer !dummy variable for ion_type array   
    CHARACTER(len=20) :: fmtx,fmtn,fmtion,fmtxrow,fmtecoef

    kc = 1
    ! READING INPUT ARRAYS FROM BINARY FILES

    ! p{1} Size of radial or scan arrays
    dimx = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{2} Size of wavenumber arrays
    dimn = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{3} Number of ions in system
    nions = INT(readvar(kc,dummy,ktype)) ; kc = kc+1
    ALLOCATE(ion_typer(dimx,nions))

    ! p{4} Flag for additional calculation on particle transport
    phys_meth = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{5} Flag for additional calculation on particle transport
    coll_flag = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{6} Flag for including rotation
    rot_flag = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{7} Number of requested solutions
    numsols = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{8} 1D integral accuracy
    relacc1 = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{9} 2D integral accuracy
    relacc2 = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{10} Number of runs before runcounter resets
    maxruns = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{11} Maximum number of integrand evaluations in 2D integration routine
    maxpts = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{12} Toroidal wave-number grid
    ALLOCATE(kthetarhos(dimn))
    kthetarhos = readvar(kc,kthetarhos,ktype) ; kc = kc+1

    ! p{13} Normalised radial coordinate
    ALLOCATE(x(dimx))
    x = readvar(kc,x,ktype) ; kc = kc+1

    ! p{14} <Ro> major radius
    ALLOCATE(Ro(dimx))
    Ro = readvar(kc,Ro,ktype) ; kc = kc+1

    ! p{15} <a> minor radius
    ALLOCATE(Rmin(dimx))
    Rmin = readvar(kc,Rmin,ktype) ; kc = kc+1

    ! p{16} B(rho) magnetic field
    ALLOCATE(Bo(dimx))
    Bo = readvar(kc,Bo,ktype) ; kc = kc+1

    ! p{17} q(rho) profile
    ALLOCATE(qx(dimx))
    qx = readvar(kc,qx,ktype) ; kc = kc+1

    ! p{18} s(rho) profile
    ALLOCATE(smag(dimx))
    smag = readvar(kc,smag,ktype) ; kc = kc+1

    ! p{19} alpha(rho) profile
    ALLOCATE(alphax(dimx))
    alphax = readvar(kc,alphax,ktype) ; kc = kc+1

    ! p{20} Machtor(rho) profile
    ALLOCATE(Machtor(dimx))
    Machtor = readvar(kc,Machtor,ktype) ; kc = kc+1

    ! p{21} Autor(rho) profile
    ALLOCATE(Autor(dimx))
    Autor = readvar(kc,Autor,ktype) ; kc = kc+1
    WHERE(ABS(Autor) < epsD) Autor = epsD

    ! p{22} Machpar(rho) profile
    ALLOCATE(Machpar(dimx))
    Machpar = readvar(kc,Machpar,ktype) ; kc = kc+1

    ! p{23} Aupar(rho) profile
    ALLOCATE(Aupar(dimx))
    Aupar = readvar(kc,Aupar,ktype) ; kc = kc+1
    WHERE(ABS(Aupar) < epsD) Aupar = epsD

    ! p{24} gammaE(rho) profile
    ALLOCATE(gammaE(dimx))
    gammaE = readvar(kc,gammaE,ktype) ; kc = kc+1
    WHERE(ABS(gammaE) < epsD) gammaE = epsD

    ! p{25} Te(rho) profile
    ALLOCATE(Tex(dimx))
    Tex = readvar(kc,Tex,ktype) ; kc = kc+1

    ! p{26} ne(rho) profile
    ALLOCATE(Nex(dimx))
    Nex = readvar(kc,Nex,ktype) ; kc = kc+1

    ! p{27} R/LTe(rho) profile
    ALLOCATE(Ate(dimx))
    Ate = readvar(kc,Ate,ktype) ; kc = kc+1

    ! p{28} R/Lne(rho) profile
    ALLOCATE(Ane(dimx))
    Ane = readvar(kc,Ane,ktype) ; kc = kc+1

    ! p{29} Flag for adiabatic electrons
    el_type = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{30} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anise(dimx))
    anise = readvar(kc,anise,ktype) ; kc = kc+1

    ! p{31} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisedr(dimx))
    danisedr = readvar(kc,danisedr,ktype) ; kc = kc+1

    ! p{32} Main ion mass
    ALLOCATE(Ai(dimx,nions))
    Ai = readvar(kc,Ai,ktype) ; kc = kc+1

    ! p{33} Main ion charge
    ALLOCATE(Zi(dimx,nions))
    Zi = readvar(kc,Zi,ktype) ; kc = kc+1

    ! p{34} Ti(rho) profiles
    ALLOCATE(Tix(dimx,nions))
    Tix = readvar(kc,Tix,ktype) ; kc = kc+1

    ! p{35} ni/ne (rho) profiles
    ALLOCATE(ninorm(dimx,nions))
    ninorm = readvar(kc,ninorm,ktype) ; kc = kc+1

    ! p{36} R/LTi(rho) profiles
    ALLOCATE(Ati(dimx,nions))
    Ati = readvar(kc,Ati,ktype) ; kc = kc+1

    ! p{37} R/Lni(rho) profiles
    ALLOCATE(Ani(dimx,nions))
    Ani = readvar(kc,Ani,ktype) ; kc = kc+1

    ! p{38} Ion types
    ALLOCATE(ion_type(dimx,nions))
    ion_type = INT(readvar(kc,ion_typer,ktype)) ; kc = kc+1
    DEALLOCATE(ion_typer)

    ! p{39} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anis(dimx,1:nions))
    anis = readvar(kc,anis,ktype) ; kc = kc+1

    ! p{40} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisdr(dimx,1:nions))
    danisdr = readvar(kc,danisdr,ktype) ; kc = kc+1

    !Allocate saved primitive output
    ALLOCATE( modewidth (dimx, dimn) )
    ALLOCATE( modeshift (dimx, dimn) )
    ALLOCATE( distan (dimx, dimn) )
    ALLOCATE( solflu (dimx, dimn) )
    ALLOCATE( ntor (dimx, dimn) )

    ALLOCATE( sol (dimx, dimn, numsols) )
    ALLOCATE( Lcirce (dimx, dimn, numsols) )
    ALLOCATE( Lpiege (dimx, dimn, numsols) )
    ALLOCATE( Lecirce (dimx, dimn, numsols) )
    ALLOCATE( Lepiege (dimx, dimn, numsols) )
    ALLOCATE( Lvcirce (dimx, dimn, numsols) )
    ALLOCATE( Lvpiege (dimx, dimn, numsols) )

    ALLOCATE( Lcirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lpiegi (dimx, dimn, nions, numsols) )
    ALLOCATE( Lecirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lepiegi (dimx, dimn, nions, numsols) )
    ALLOCATE( Lvcirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lvpiegi (dimx, dimn, nions, numsols) )

    IF (phys_meth /= 0.0) THEN
       ALLOCATE( Lcircgte (dimx, dimn, numsols) )
       ALLOCATE( Lpieggte (dimx, dimn, numsols) )
       ALLOCATE( Lcircgne (dimx, dimn, numsols) )
       ALLOCATE( Lpieggne (dimx, dimn, numsols) )
       ALLOCATE( Lcircgue (dimx, dimn, numsols) )
       ALLOCATE( Lpieggue (dimx, dimn, numsols) )
       ALLOCATE( Lcircce (dimx, dimn, numsols) )
       ALLOCATE( Lpiegce (dimx, dimn, numsols) )
       ALLOCATE( Lcircgti (dimx, dimn, nions, numsols) )
       ALLOCATE( Lpieggti (dimx, dimn, nions, numsols) )
       ALLOCATE( Lcircgni (dimx, dimn, nions, numsols) )
       ALLOCATE( Lpieggni (dimx, dimn, nions, numsols) )
       ALLOCATE( Lcircgui (dimx, dimn, nions, numsols) )
       ALLOCATE( Lpieggui (dimx, dimn, nions, numsols) )
       ALLOCATE( Lcircci (dimx, dimn, nions, numsols) )
       ALLOCATE( Lpiegci (dimx, dimn, nions, numsols) )
       IF (phys_meth == 2) THEN
          ALLOCATE( Lecircgte (dimx, dimn, numsols) )
          ALLOCATE( Lepieggte (dimx, dimn, numsols) )
          ALLOCATE( Lecircgne (dimx, dimn, numsols) )
          ALLOCATE( Lepieggne (dimx, dimn, numsols) )
          ALLOCATE( Lecircgue (dimx, dimn, numsols) )
          ALLOCATE( Lepieggue (dimx, dimn, numsols) )
          ALLOCATE( Lecircce (dimx, dimn, numsols) )
          ALLOCATE( Lepiegce (dimx, dimn, numsols) )
          ALLOCATE( Lecircgti (dimx, dimn, nions, numsols) )
          ALLOCATE( Lepieggti (dimx, dimn, nions, numsols) )
          ALLOCATE( Lecircgni (dimx, dimn, nions, numsols) )
          ALLOCATE( Lepieggni (dimx, dimn, nions, numsols) )
          ALLOCATE( Lecircgui (dimx, dimn, nions, numsols) )
          ALLOCATE( Lepieggui (dimx, dimn, nions, numsols) )
          ALLOCATE( Lecircci (dimx, dimn, nions, numsols) )
          ALLOCATE( Lepiegci (dimx, dimn, nions, numsols) )
       ENDIF
    ENDIF

    ALLOCATE(rdummyn(dimx,dimn))
    ALLOCATE(idummyn(dimx,dimn))
    ALLOCATE(rdummynumsols(dimx,dimn,numsols))
    ALLOCATE(idummynumsols(dimx,dimn,numsols))
    ALLOCATE(rdummynions(dimx,dimn,nions,numsols))
    ALLOCATE(idummynions(dimx,dimn,nions,numsols))

    WRITE(fmtx, '(A)') '(G15.7)'
    WRITE(fmtxrow,'(A,I0,A)') '(',dimx,'G15.7)'
    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
    WRITE(fmtion,'(A,I0, A)') '(',nions,'G15.7)'
    WRITE(fmtecoef,'(A,I0, A)') '(',nions+1,'G15.7)'

    OPEN(unit=myunit, file="output/primitive/rsolflu.dat", action="read")
    READ(myunit,fmtn) ((rdummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/isolflu.dat", action="read")
    READ(myunit,fmtn) ((idummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    solflu = CMPLX(rdummyn,idummyn)

    OPEN(unit=myunit, file="output/primitive/distan.dat", action="read")
    READ(myunit,fmtn) ((distan(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rmodewidth.dat", action="read")
    READ(myunit,fmtn) ((rdummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/imodewidth.dat", action="read")
    READ(myunit,fmtn) ((idummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    modewidth(i,j) = CMPLX(rdummyn,idummyn)

    OPEN(unit=myunit, file="output/primitive/rmodeshift.dat", action="read")
    READ(myunit,fmtn) ((rdummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/imodeshift.dat", action="read")
    READ(myunit,fmtn) ((idummyn(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)
    modeshift(i,j) = CMPLX(rdummyn,idummyn)

    OPEN(unit=myunit, file="output/primitive/ntor.dat", action="read")
    READ(myunit,fmtn) ((ntor(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rsol.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/isol.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    sol = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLcirce.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLcirce.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lcirce = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLpiege.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLpiege.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lpiege = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLecirce.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLecirce.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lecirce = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLepiege.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLepiege.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lepiege = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLvcirce.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLvcirce.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lvcirce = CMPLX(rdummynumsols,idummynumsols)

    OPEN(unit=myunit, file="output/primitive/rLvpiege.dat", action="read")
    READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLvpiege.dat", action="read")
    READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
    Lvpiege = CMPLX(rdummynumsols,idummynumsols)


    IF (phys_meth /= 0.0) THEN

       OPEN(unit=myunit, file="output/primitive/rLcircgne.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgne.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lcircgne = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLpieggne.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggne.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lpieggne = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLcircgue.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgue.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lcircgue = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLpieggue.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggue.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lpieggue = CMPLX(rdummynumsols,idummynumsols)


       OPEN(unit=myunit, file="output/primitive/rLcircgte.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgte.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lcircgte = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLpieggte.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggte.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lpieggte = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLcircce.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircce.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lcircce = CMPLX(rdummynumsols,idummynumsols)

       OPEN(unit=myunit, file="output/primitive/rLpiegce.dat", action="read")
       READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpiegce.dat", action="read")
       READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       Lpiegce = CMPLX(rdummynumsols,idummynumsols)
       !!
       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/primitive/rLecircgne.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgne.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lecircgne = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLepieggne.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggne.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lepieggne = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLecircgue.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgue.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lecircgue = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLepieggue.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggue.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lepieggue = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLecircgte.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgte.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lecircgte = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLepieggte.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggte.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lepieggte = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLecircce.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircce.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lecircce = CMPLX(rdummynumsols,idummynumsols)

          OPEN(unit=myunit, file="output/primitive/rLepiegce.dat", action="read")
          READ(myunit,fmtn) (((rdummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepiegce.dat", action="read")
          READ(myunit,fmtn) (((idummynumsols(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
          Lepiegce = CMPLX(rdummynumsols,idummynumsols)
       ENDIF
    ENDIF

    OPEN(unit=myunit, file="output/primitive/rLcirci.dat", action="read")
    READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLcirci.dat", action="read")
    READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    Lcirci = CMPLX(rdummynions,idummynions)

    OPEN(unit=myunit, file="output/primitive/rLpiegi.dat", action="read")
    READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLpiegi.dat", action="read")
    READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    Lpiegi = CMPLX(rdummynions,idummynions)

    OPEN(unit=myunit, file="output/primitive/rLecirci.dat", action="read")
    READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLecirci.dat", action="read")
    READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    Lecirci = CMPLX(rdummynions,idummynions)

    OPEN(unit=myunit, file="output/primitive/rLepiegi.dat", action="read")
    READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    OPEN(unit=myunit, file="output/primitive/iLepiegi.dat", action="read")
    READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
    Lepiegi = CMPLX(rdummynions,idummynions)

    IF (phys_meth /= 0.0) THEN
       OPEN(unit=myunit, file="output/primitive/rLcircgni.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgni.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lcircgni = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLpieggni.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggni.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lpieggni = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLcircgui.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgui.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lcircgui = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLpieggui.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggui.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lpieggui = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLcircgti.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircgti.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lcircgti = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLpieggti.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpieggti.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lpieggti = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLcircci.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLcircci.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lcircci = CMPLX(rdummynions,idummynions)

       OPEN(unit=myunit, file="output/primitive/rLpiegci.dat", action="read")
       READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       OPEN(unit=myunit, file="output/primitive/iLpiegci.dat", action="read")
       READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       Lpiegci = CMPLX(rdummynions,idummynions)

       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/primitive/rLecircgni.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgni.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lecircgni = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLepieggni.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggni.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lepieggni = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLecircgui.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgui.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lecircgui = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLepieggui.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggui.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lepieggui = CMPLX(rdummynions,idummynions)


          OPEN(unit=myunit, file="output/primitive/rLecircgti.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircgti.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lecircgti = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLepieggti.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepieggti.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lepieggti = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLecircci.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLecircci.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lecircci = CMPLX(rdummynions,idummynions)

          OPEN(unit=myunit, file="output/primitive/rLepiegci.dat", action="read")
          READ(myunit,fmtn) ((((rdummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          OPEN(unit=myunit, file="output/primitive/iLepiegci.dat", action="read")
          READ(myunit,fmtn) ((((idummynions(i,j,k,l),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
          Lepiegci = CMPLX(rdummynions,idummynions)
       ENDIF
    ENDIF

    DEALLOCATE(rdummyn)
    DEALLOCATE(idummyn)
    DEALLOCATE(rdummynumsols)
    DEALLOCATE(idummynumsols)
    DEALLOCATE(rdummynions)
    DEALLOCATE(idummynions)

  END SUBROUTINE data_init

  SUBROUTINE deallocate_all()
    !DEALLOCATE INPUT ARRAYS
    DEALLOCATE(x)
    DEALLOCATE(Ro)
    DEALLOCATE(Rmin)
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
    DEALLOCATE(Machpar)
    DEALLOCATE(gammaE)
    DEALLOCATE(Aupar)
    DEALLOCATE(Autor)
    DEALLOCATE(ninorm)
    DEALLOCATE(Ai)
    DEALLOCATE(Zi)

    DEALLOCATE(modewidth)
    DEALLOCATE(modeshift)
    DEALLOCATE(distan)
    DEALLOCATE(ntor)
    DEALLOCATE(solflu)

    DEALLOCATE(sol)
    DEALLOCATE(Lcirce)
    DEALLOCATE(Lpiege)
    DEALLOCATE(Lecirce)
    DEALLOCATE(Lepiege)

    DEALLOCATE(Lcirci)
    DEALLOCATE(Lpiegi)
    DEALLOCATE(Lecirci)
    DEALLOCATE(Lepiegi)
    DEALLOCATE(Lvpiege)
    DEALLOCATE(Lvcirce)
    DEALLOCATE(Lvcirci)
    DEALLOCATE(Lvpiegi)

    DEALLOCATE(tau)
    DEALLOCATE(Nix)
    DEALLOCATE(Zeffx)
    DEALLOCATE(qprim)
    DEALLOCATE(Anue)
    DEALLOCATE(wg)
    DEALLOCATE(Ac)
    DEALLOCATE(csou)
    DEALLOCATE(cref) 
    DEALLOCATE(cthe)
    DEALLOCATE(cthi)
    DEALLOCATE(Rhoe)
    DEALLOCATE(Rhoi)
    DEALLOCATE(de)
    DEALLOCATE(di)
    DEALLOCATE(ktetasn)
    DEALLOCATE(rhostar)
    DEALLOCATE(ft)
    DEALLOCATE(fc)

    DEALLOCATE(Machi)
    DEALLOCATE(Mache)
    DEALLOCATE(Aui)
    DEALLOCATE(Aue)
    DEALLOCATE(th)
    DEALLOCATE(phi)

    DEALLOCATE(dphidr)
    DEALLOCATE(dphidth)
    DEALLOCATE(npol)
    DEALLOCATE(ecoefs)
    DEALLOCATE(ecoefsgau)
    DEALLOCATE(cftrans)
    DEALLOCATE(nepol)
    DEALLOCATE(Anipol)
    DEALLOCATE(Anepol)

    DEALLOCATE(tpernorme)
    DEALLOCATE(tpernormi)
    DEALLOCATE(tpernormifunc)
    DEALLOCATE(coefi)
    DEALLOCATE(mi)
    DEALLOCATE(epsilon)

    IF (phys_meth /= 0.0) THEN
       DEALLOCATE(Lcircgte)
       DEALLOCATE(Lpieggte)
       DEALLOCATE(Lcircgne)
       DEALLOCATE(Lpieggne)
       DEALLOCATE( Lcircgue )
       DEALLOCATE( Lpieggue )
       DEALLOCATE(Lcircce)
       DEALLOCATE(Lpiegce)
       DEALLOCATE(Lcircgti)
       DEALLOCATE(Lpieggti)
       DEALLOCATE(Lcircgni)
       DEALLOCATE(Lpieggni)
       DEALLOCATE( Lcircgui )
       DEALLOCATE( Lpieggui )
       DEALLOCATE(Lcircci)
       DEALLOCATE(Lpiegci)
       IF (phys_meth == 2) THEN
          DEALLOCATE(Lecircgte)
          DEALLOCATE(Lepieggte)
          DEALLOCATE(Lecircgne)
          DEALLOCATE(Lepieggne)
          DEALLOCATE( Lecircgue )
          DEALLOCATE( Lepieggue )
          DEALLOCATE(Lecircce)
          DEALLOCATE(Lepiegce)
          DEALLOCATE(Lecircgti)
          DEALLOCATE(Lepieggti)
          DEALLOCATE(Lecircgni)
          DEALLOCATE(Lepieggni)
          DEALLOCATE( Lecircgui )
          DEALLOCATE( Lepieggui )
          DEALLOCATE(Lecircci)
          DEALLOCATE(Lepiegci)
       ENDIF
    ENDIF

  END SUBROUTINE deallocate_all

  SUBROUTINE outputascii()

    CHARACTER(len=20) :: fmtx,fmtn,fmtion,fmtxrow,fmtecoef
    INTEGER :: i,j,k,myunit=700
    WRITE(fmtx, '(A)') '(G15.7)'
    WRITE(fmtxrow,'(A,I0,A)') '(',dimx,'G15.7)'
    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
    WRITE(fmtion,'(A,I0, A)') '(',nions,'G15.7)'
    WRITE(fmtecoef,'(A,I0, A)') '(',nions+1,'G15.7)'

    OPEN(unit=myunit, file="output/modeflag.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (modeflag(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/phi.dat", action="write", status="replace")
    WRITE(myunit,fmtxrow) ((phi(i,j),i=1,dimx),j=1,ntheta) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/npol.dat", action="write", status="replace")
    WRITE(myunit,fmtxrow) (((npol(i,j,k),i=1,dimx),j=1,ntheta),k=1,nions) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ecoefs.dat", action="write", status="replace")
    WRITE(myunit,fmtecoef) (((ecoefs(i,j,k),j=0,nions),i=1,dimx),k=1,numecoefs) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/gam_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((gam_GB(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ome_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ome_GB(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/gam_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((gam_SI(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ome_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ome_SI(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    IF (phys_meth /= 0.0) THEN
       OPEN(unit=myunit, file="output/cke.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (cke(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/cki.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((cki(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/dfe_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (dfe_SI(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/dfi_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((dfi_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vte_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (vte_SI(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vti_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((vti_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vce_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (vce_SI(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vci_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((vci_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)
       !!
       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/ceke.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (ceke(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/ceki.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((ceki(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/vene_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (vene_SI(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/veni_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((veni_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/chiee_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (chiee_SI(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/chiei_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((chiei_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/vece_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (vece_SI(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/veci_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((veci_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)
       ENDIF
    ENDIF

    OPEN(unit=myunit, file="output/epf_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (epf_SI(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ipf_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/epf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (epf_GB(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ipf_GB(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (eef_SI(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ief_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ief_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (eef_GB(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ief_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ief_GB(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/epf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((epf_cm(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ipf_cm(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((eef_cm(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/ief_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ief_cm(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(myunit)

  END SUBROUTINE outputascii


END PROGRAM qlk_makeflux
