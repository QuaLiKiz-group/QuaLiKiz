PROGRAM qlk_redoQL
  !Standalone driver for qualikiz. Detailed code info in call_qualikiz subroutine header

  USE kind
  USE mod_io_management !MPI is included in this module
  USE mod_make_io
  USE datmat
  USE datcal
  USE calcroutines  
  USE asymmetry
  USE mod_saturation
  USE FLRterms
  USE mod_fluidsol
  USE mod_contour
  USE QLflux
  USE mod_fonct

  IMPLICIT NONE

  !Local data dictionary.

  !Time measuring variables
  REAL(kind=DBL) :: cputime1, cputime2, tpstot
  INTEGER :: time1, time2, timetot, freq, cputimetot,TotTask,count
  CHARACTER(len=20) :: myfmt
  INTEGER:: p,nu !counter for loop over coordinates. p is radius (or general scan), nu is wavenumber
  REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Lambe, Nue !physical quantities used for the derivations, then discarded

  !MPI variables:
  INTEGER :: ierror, nproc, myrank
  INTEGER,DIMENSION(MPI_STATUS_SIZE) :: status
  INTEGER, DIMENSION(:), ALLOCATABLE :: Order,wavenum,radcoord

  ! Old solution for non-reset runs
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldrsol,oldisol,oldrfdsol,oldifdsol

  LOGICAL :: exist1, exist2, exist3, exist4, exist5 !used for checking for existence of files

  !DEBUGGING
  CHARACTER(len=20) :: fmtx,fmtn,fmtion,fmtintion,fmtxrow,fmtecoef
  INTEGER :: i,j,k,l,myunit=700,stat

  CALL mpi_init(ierror)
  CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
  CALL mpi_comm_rank(mpi_comm_world,myrank,ierror)

  ! Begin time measurement
  CALL SYSTEM_CLOCK(time1)
  CALL CPU_TIME(cputime1)

  ! Read input and initialize all arrays
  CALL data_init()

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
  cref(:) = SQRT(qe*1.d3/mp) !Cref=sqrt(2x1keV/mD)=sqrt(1keV/mp) used to normalized gammaEin, Machparin, Auparin, Machtorin, Autorin
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

  Lambe(:) = 15.2_DBL - 0.5_DBL*LOG(0.1_DBL*Nex(:)) + LOG(Tex(:))  !Coulomb constant and collisionality. Wesson 2nd edition p661-663
  Nue(:) = 1._DBL/(1.09d-3) *Zeffx(:)*Nex(:)*Lambe(:)/(Tex(:))**1.5_DBL 

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

  !Set ETG_flag for when kthetarhos > 2 
  WHERE (kthetarhos > ETGk) 
     ETG_flag = .TRUE.
  ELSEWHERE
     ETG_flag = .FALSE.
  ENDWHERE

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

  ! ntor is the toroidal wave-number grid for each position of the scan
  DO p=1,dimx
     DO nu=1,dimn
        ntor(p,nu) = kthetarhos(nu)*x(p)/(qx(p)*rhostar(p))
     ENDDO
  ENDDO
  ntor=ANINT(ntor)

  DEALLOCATE(Lambe)
  DEALLOCATE(Nue)

  !Allocation and initialization of calculated arrays (named "output")
  !Complex 1D arrays. These ion arrays are not the final output but used in an intermediate step
  ALLOCATE(fonxcirci(nions))
  ALLOCATE(fonxpiegi(nions))
  ALLOCATE(fonxecirci(nions))
  ALLOCATE(fonxepiegi(nions))
  ALLOCATE(fonxvcirci(nions))
  ALLOCATE(fonxvpiegi(nions))
  !Real 2D arrays
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
  ALLOCATE( modewidth (dimx, dimn) ); modewidth=0
  ALLOCATE( modeshift (dimx, dimn) ); modeshift=0
  ALLOCATE( newsolflu (dimx, dimn) ); newsolflu=0
  ALLOCATE( newmodewidth (dimx, dimn) ); newmodewidth=0
  ALLOCATE( newmodeshift (dimx, dimn) ); newmodeshift=0

  !Complex 3D arrays with 3rd dimension equal to number of solutions searched for
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
     ALLOCATE(fonxcircgti(nions))
     ALLOCATE(fonxpieggti(nions))
     ALLOCATE(fonxcircgni(nions))
     ALLOCATE(fonxpieggni(nions))
     ALLOCATE(fonxcircgui(nions))
     ALLOCATE(fonxpieggui(nions))
     ALLOCATE(fonxcircci(nions))
     ALLOCATE(fonxpiegci(nions))
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
        ALLOCATE(fonxecircgti(nions))
        ALLOCATE(fonxepieggti(nions))
        ALLOCATE(fonxecircgni(nions))
        ALLOCATE(fonxepieggni(nions))
        ALLOCATE(fonxecircgui(nions))
        ALLOCATE(fonxepieggui(nions))
        ALLOCATE(fonxecircci(nions))
        ALLOCATE(fonxepiegci(nions))
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


  CALL init_asym() !subroutine in datcal. Initializes variables used for Fried-Conte function asymptotic expansions

  CALL calcphi() !calculate poloidal density asymmetries due to rotation and temp anisotropy
  CALL calcdphi() !calculate radial and poloidal gradient of 2D electrostatic potential
  CALL calccoefs() !calculate e# LFS to FSA transport coefficients, FSA R/Ln, FSA n, asym factor, and also 2D R/Ln for all ion species
  !The poloidal asymmetry terms take around 10ms to calculate

  weidcount = 0 !Initialize dispersion relation function call counter. Used for debugging purposes
  Nsolrat = 0   !Initializes count of failed solutions (see calcroutines). Rarely occurs.

  !Total number of jobs to run
  TotTask=dimx*dimn

  !Initialize output related arrays
  IF (myrank == 0) THEN
     ALLOCATE(Order(TotTask))
     ALLOCATE(wavenum(TotTask))
     ALLOCATE(radcoord(TotTask))
  END IF
  ! ALLOCATE(BUFF(TotTask*TaskRN))
  !Initialize counter for BUFF
  count=0 

  !! NOW THE MAGIC HAPPENS!! This subroutine is contained below
  !! Distributes tasks to all processors and calculates output
  CALL DistriTask(TotTask,nproc,myrank,count)

  !! IF SET, SAVING OUTPUT AS BINARY FILES. BUFF files are created in ./output
!!$  IF (print_binary) THEN
!!$     CALL writeBUFF(BUFF,Order,TotTask,count,myrank)
!!$  ENDIF

  IF (myrank==0) THEN
     WRITE(stdout,"(A)") '*** Collecting output'
  ENDIF

  !Consolidate all arrays into rank=0 for QL flux integrals and output
  IF (nproc > 1) THEN
     CALL collectarraysredo()
  ENDIF

  !If rank0, then carry out the saturation rules and output final results
  IF (myrank==0) THEN 

     WRITE(stdout,*)     
     WRITE(stdout,"(A)") '*** Carrying out saturation rules'

     CALL allocate_endoutput()
     CALL saturation(0) !set 1 for ignoring electron modes, and 2 for ignoring ion modes
     CALL outputascii
     CALL deallocate_endoutput()     
  ENDIF

  !Write run info in output file
  CALL CPU_TIME(cputime2)
  cputimetot = (cputime2-cputime1)

  !WRITE(stdout,"(A,I0,A,I0,A)") '*** Rank ',myrank,' finished! CPU time: ', cputimetot,' s' 
  CALL SYSTEM_CLOCK(time2)
  CALL SYSTEM_CLOCK(count_rate=freq)

  CALL MPI_Barrier(mpi_comm_world,ierror)
  timetot = REAL(time2-time1) / REAL(freq)

  !WRITE(stdout,*) 'Z function was called ',weidcount,' times' 
  !WRITE(stdout,*)
  !WRITE(stdout,"(A,I0,A,I0)") '*** time: ',timetot, 'seconds, for rank = ',myrank 
  !WRITE(stdout,"(A,I0)") '*** End of job for rank ',myrank
  !IF (myrank==0) WRITE(stdout,"(A,I0,A)") 'We have missed ',Nsolrat,' eigenvalues.'

  IF (myrank==0) THEN 
     WRITE(stdout,"(A,I0,A)") 'Hurrah! Job completed! Total time = ',timetot,' s'  !final write
  ENDIF

  !Deallocating all
  CALL deallocate_all()

  ! MPI finalization
  CALL mpi_finalize(ierror)

CONTAINS 

  SUBROUTINE calc(p,nu)
    INTEGER, INTENT(IN) :: p,nu

    ! Variables for fluid solution
    REAL(kind=DBL) :: fluome

    ! Variables for defining contour locations
    REAL(kind=DBL) :: ma, wpi, wpe
    COMPLEX(kind=DBL) :: om0
    COMPLEX(kind=DBL) :: C, Centre
    REAL(kind=DBL) :: L

    ! Variables for collecting the solutions, and their status flags
    COMPLEX(kind=DBL), DIMENSION(numsols) :: soll, fdsoll, solltmp, fdsolltmp
    COMPLEX(kind=DBL) :: fdsollj, omin, fonout,newsol,newfdsol,soltest
    REAL(kind=DBL),    DIMENSION(numsols) :: imagsol
    INTEGER :: NN, NNN, i,j
    LOGICAL :: issol
    REAL(kind=DBL) :: kteta,maxdia

    !INITIALIZATION OF VARIABLES************************

    !IF ( ( p /= 1 ) .OR. ( nu /= 1)) THEN
    !STOP
    !ENDIF

    !Initialize prefactor
    nwg = ntor(p,nu)*wg(p)

    !Initialize distance between rational surfaces
    d = ABS(1./(ntor(p,nu)*(epsD+ABS(qprim(p)))))
    distan(p,nu) = d !save for output

    !Initialize variables used in calcfonctp and FLR module
    kteta     = ntor(p,nu)*ktetasn(p)
    ktetaRhoe = kteta*Rhoe(p)

    ktetaRhoi(:) = kteta*Rhoi(p,:)

    pnFLR = p !save scan index for use in FLRterms module

    !Initialize variables for calculating the mode width

    Joe2 = BESEI0(ktetaRhoe * ktetaRhoe)
    Jobane2 = BESEI0 (kteta*kteta*de(p)*de(p))

    DO i = 1,nions
       Joi2(i) = BESEI0 (ktetaRhoi(i) * ktetaRhoi(i))
       Jobani2(i) = BESEI0 (kteta*kteta*di(p,i)*di(p,i))
    ENDDO

    normkr = normkrfac*ntor(p,nu) !sets boundary in kr integrations

    !**************************************************

    !CALCULATION PHASE. TWO OPTIONS: calculate from scratch, or start directly from newton solver with inputs from a previous run

    !Calculate the fluid frequency and growth rates. Still testing phase where various approaches tried
    CALL fluidsol(p,nu,fluome)! with or without rotation fluome used only for contour limit (warning this limit: should it account for rotation impact ?)

    wpi = wg(p)*Tix(p,1)/Zi(p,1)* (Ati(p,1) + Ani(p,1))
    wpe = wg(p)*Tex(p)/Ze * (Ate(p) + Ane(p)) 

    ! Set maximum absolute diamagnetic frequency for limit of contour search
    IF ( ABS(wpi) > ABS(wpe)) THEN
       maxdia = wpi
    ELSE
       maxdia=wpe
    ENDIF

    !Calculate mode width  and shift and save in mwidth and shift. nu dependent terms are in the Bessel function defined above, which are global variables and seen by calcwidth
    IF (rot_flag == 1) THEN
       CALL fluidsolrot(p,nu)! with rotation solflu used for width and shift calculations
       solflu(p,nu) = omeflu !The modewidth variable is output, while mwidth is used within code
       CALL calcwidthrot(ETG_flag(nu),p,nu)
       modewidth(p,nu) = mwidth !The modewidth variable is output, while mwidth is used within code
       CALL calcshift(ETG_flag(nu),p,nu)
       modeshift(p,nu) = mshift !The modeshift variable is output, while shift is used in code
    ELSE
       IF ( ETG_flag(nu) == .FALSE. ) THEN !ITG/TEM case
          solflu(p,nu) = CMPLX(maxdia/wg(p),fluome)
       ELSE !ETG case
          solflu(p,nu) = CMPLX(wpe/wg(p),fluome)
       ENDIF
       CALL calcwidth(ETG_flag(nu),p)
       modewidth(p,nu) = mwidth !The modewidth variable is output, while mwidth is used within code
       modeshift(p,nu) = 0.0 !no shift in absence of rotation
    ENDIF

    !THIS IS HERE FOR TESTING PURPOSES AT THE MOMENT
    ALLOCATE(newsolflu(dimx,dimn))
    ALLOCATE(newmodewidth(dimx,dimn))
    ALLOCATE(newmodeshift(dimx,dimn))
    CALL newfluidsol(p,nu,newsolflu,newmodewidth,newmodeshift)
    DEALLOCATE(newsolflu)
    DEALLOCATE(newmodewidth)
    DEALLOCATE(newmodeshift)

    ! Calculate flux-surface-averaging of poloidal asymmetry coefficients including eigenmode
    CALL makeecoefsgau(p,nu)

    !! Carry out FLR calculation with Bessel function integrations over kr
    IF (rot_flag == 1) THEN
       CALL makeFLRtermsrot(p,nu)
       FLRep(p,nu) = Joe2p !output
       FLRip(p,nu,:) = Joi2p(:)
       FLRec(p,nu) = Joe2c 
       FLRic(p,nu,:) = Joi2c(:) 
    ELSE
       CALL makeFLRterms(p,nu)
       FLRep(p,nu) = Joe2p !output
       FLRip(p,nu,:) = Joi2p(:)
       FLRec(p,nu) = Joe2c 
       FLRic(p,nu,:) = Joi2c(:) 
    ENDIF
    !*************************


    !Set the transit frequency
    qRd = qx(p)*Ro(p)*d*SQRT(3._DBL)

    !Initialize often used prefactors
    Athe=mwidth*cthe(p)/qRd

    Athi(:)=mwidth*cthi(p,:)/qRd

    !Set the bounce frequency
    omega2bar = pi/2.*SQRT(x(p)*Rmin(p)/(2.*Ro(p)))

    !Initialize ma value used in contour choice based on fluid solution
    ma = MAX(AIMAG(solflu(p,nu)),10.)
    om0 = ma*ci
    C    = om0 *0.5

    !Set maximum bound of contour locations based on diamagnetic frequency
    ommax(p,nu) = om0 + REAL(solflu(p,nu))/1 !divisor sets max omega ratio to diamagnetic frequency

    !Set maximum bound of contour locations based on heritage assumptions (not clear to me (JC))

    !    IF (ETG_flag(nu) .EQV. .FALSE.) THEN
    !       ommax(p,nu) = om0+MAX(2.*ABS(Tex(p)),2.*ABS(Tix(p,1)/Zi(p,1)),2.*ABS(Athi(1)/(nwg)))
    !    ELSE
    !       ommax(p,nu) = om0+MIN(MAX(2.*ABS(Tex(p)),1.5*ABS(Athe/(nwg))),12.)
    !    ENDIF

    omegmax=ommax(p,nu) ! used in calculsol for maximum boundary of allowed solution

    ! The solution for this (p,nu) pair is now saved
    DO j=1,numsols
       soll(j)=sol(p,nu,j) 
       fdsoll(j)=fdsol(p,nu,j)
       IF ( soll(j) /= (0.,0.) ) THEN
          !calculate integrals for quasilinear flux. Small percentage of total calculation
          CALL make_QLflux(p,nu,soll(j),fdsollj)
          issol = .TRUE.
       ELSE
          issol = .FALSE.
       ENDIF
       CALL save_qlfunc(p,nu,j,issol)
    ENDDO

    !Reorder (descending) the solutions
    imagsol = AIMAG(soll)
    CALL dsort(imagsol,imagsol,numsols,-1)
    !Save growth rates to output array (normalized to nwg)

    gamma(p,nu,:) = imagsol

    !Check if in ETG regime. Commented out since replaced with simple ktheta limit
    !    CALL ETGcheck(ETG_flag,soll,kteta,p,nu)

  END SUBROUTINE calc

  SUBROUTINE DistriTask(NumTasks,numprocs,rank,icount)
    IMPLICIT NONE
    !     include 'mpif.h'
    INTEGER,DIMENSION(:),ALLOCATABLE :: Request,OK,Finish
    INTEGER,INTENT(IN) :: NumTasks,numprocs,rank
    INTEGER,INTENT(INOUT) :: icount
    REAL(kind=DBL) :: tps
    INTEGER :: iradcoord,iwavenum, Task, NoTask, iloop, ierr
    INTEGER :: one=1, minusone=-1,tmpint
    LOGICAL :: Complete = .FALSE.
    LOGICAL :: Complete0 = .FALSE.
    LOGICAL :: ressend = .FALSE.
    LOGICAL, DIMENSION(:), ALLOCATABLE :: resready

    !Initialize for all procs
    Task = 0
    NoTask = 0
    tpstot=0 !initialize time

    ! Code for Master processor
    IF (rank==0) THEN
       ! Initialization phase for Master arrays. 

       ALLOCATE(Request(numprocs-1))
       ALLOCATE(resready(numprocs-1))
       ALLOCATE(OK(numprocs-1))
       ALLOCATE(Finish(0:numprocs-1))
       Finish(:) = 0
       resready(:)=.FALSE.
       ! Initialize OMP
       CALL OMP_SET_DYNAMIC(.FALSE.) 
       CALL OMP_SET_NUM_THREADS(2)
       !       WRITE(*,*) 'threads=',omp_get_num_threads()
       !$OMP PARALLEL
       !$OMP SECTIONS

       !$OMP SECTION
       !****************************************************
       !****** Distributor thread on Master processor ******
       !****************************************************
       ! Master receives an integer "1" from each of the slaves in an unblocked receive.
       ! This sign of life from the slaves is checked in the while loop below
       DO iloop=1,numprocs-1
          CALL MPI_Irecv(OK(iloop),1,MPI_INTEGER,iloop,100+iloop,MPI_COMM_WORLD,Request(iloop),ierr)
       ENDDO

       !This loop distributes tasks to slaves. The +numprocs is necessary since Task will
       !increment one extra time for each slave processor
       DO WHILE (Task<=NumTasks+numprocs-1)

          DO iloop=1,numprocs-1
             !If the slave processor number iloop still has tasks to do, then carry on
             IF (Finish(iloop) /= 1) THEN               
                !Check if processor iloop is ready to proceed 
                !(if it's the first task or completed its previous task)
                CALL MPI_Test(Request(iloop),complete,status,ierr)
                IF (Complete) THEN 
                   IF (resready(iloop)) THEN
                      !     CALL getresults(iloop,iloop)                                      
                   ENDIF
                   Task=Task+1 !Prepare next task number
                   resready(iloop)=.TRUE. !results will now be saved before next distribution
                   IF (Task<=NumTasks) THEN
                      !A valid task number is sent from the Master to the slave
                      CALL MPI_SSend(Task,1,MPI_INTEGER,iloop,100+iloop,MPI_COMM_WORLD,ierr)
                      !Following MPI_Irecv, when the task is completed, 
                      !the next MPI_Test will provide Complete=.TRUE. for iloop
                      CALL MPI_Irecv(OK(iloop),1,MPI_INTEGER,iloop,100+iloop,MPI_COMM_WORLD,Request(iloop),ierr)
                      Order(Task)=iloop+1
                      wavenum(Task)=(Task-1)/(dimx) + 1
                      radcoord(Task)=MOD((Task-1),dimx) + 1
                   ELSE ! No more tasks. Send "-1" to slave, signalling that the tasks are done
                      CALL MPI_SSend(minusone,1,MPI_INTEGER,iloop,100+iloop,MPI_COMM_WORLD,ierr)
                      Finish(iloop) = 1
                   ENDIF
                ENDIF
             ENDIF
          ENDDO

          !  Manage thread timeshare

          !Now coordinate with 2nd worker thread on Master processor,
          !which is also carrying out tasks. 

          !$OMP FLUSH(Complete0)
          IF (( Finish(0) /= 1 ) .AND. Complete0) THEN         
             !If 2nd thread initialized or completed previous task, then increment task
             !and switch flag so that 2nd thread will launch a new task
             Task = Task + 1
             IF (Task <= NumTasks) THEN
                NoTask = Task
                Complete0 = .FALSE.
                Order(Task)=1
             ELSE
                Finish(0) = 1
             ENDIF
          ENDIF
       ENDDO

       !$OMP SECTION
       !**************************************
       !***Worker thread on Master processor***
       !**************************************
       complete0=.TRUE. 
       DO
          !$OMP FLUSH(Finish)
          IF (Finish(0)==1) EXIT
          !$OMP FLUSH(complete0)
          IF ( .NOT. Complete0) THEN
             tps=MPI_Wtime()
             iwavenum=(NoTask-1)/(dimx) + 1
             iradcoord=MOD((NoTask-1),dimx) + 1
             CALL calc(iradcoord,iwavenum)
             !             CALL fill_BUFF(BUFF,icount,iradcoord,iwavenum,TaskRN)
             tps=MPI_Wtime()-tps
             tpstot=tpstot+tps  
             IF (verbose .EQV. .TRUE.) THEN
                WRITE(stdout,300) rank,NoTask,tps,tpstot,iradcoord,iwavenum
             ENDIF
300          FORMAT(1x,'rank ',I3,' NoTask ',I3,' time ',F7.3,', total time ',F7.3,' (p,nu)=(',I3,',',I3,')')
             Complete0 = .TRUE.
          ENDIF
       ENDDO

       !$OMP END SECTIONS
       !$OMP END PARALLEL

    ELSE

       !********************
       !****** SLAVES ******
       !********************

       !The slaves run tasks until they receive NoTask=-1 from the coordinating Master
       DO WHILE (NoTask /= -1)
          !The slaves notify the Master that they are ready for a new task
          CALL MPI_SSend(one,1,MPI_INTEGER,0,100+rank,MPI_COMM_WORLD,ierr)
          IF (ressend) THEN
             ! CALL sendresults(iradcoord,iwavenum,rank) !copy over output to rank0 which will have full arrays for later             
          ENDIF
          !The slaves receive their task number back from the Master.
          CALL MPI_Recv(NoTask,1,MPI_INTEGER,0,100+rank,MPI_COMM_WORLD,status,ierr)
          IF (NoTask /= -1) THEN
             !If a task is assigned, the slaves compute the task. The specific task to
             !be carried out (i.e. wavenumber and 'radius') depends on the value of NoTask.
             tps=MPI_Wtime()
             iwavenum=(NoTask-1)/(dimx) + 1
             iradcoord=MOD((NoTask-1),dimx) + 1
             CALL calc(iradcoord,iwavenum)
             ressend=.TRUE.
             !             CALL fill_BUFF(BUFF,icount,iradcoord,iwavenum,TaskRN)
             tps=MPI_Wtime()-tps
             tpstot=tpstot+tps
             IF (verbose .EQV. .TRUE.) THEN
                WRITE(stdout,301) rank,NoTask,tps,tpstot,iradcoord,iwavenum
             ENDIF
301          FORMAT(1x,'rank ',I3,' NoTask ',I3,' time ',F7.3,', total time ',F7.3,' (p,nu)=(',I3,',',I3,')')
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE DistriTask

  SUBROUTINE data_init()
    !Read data, allocate input and output arrays
    INTEGER, PARAMETER :: ktype = 1 ! BINARY FILES
    INTEGER :: kc
    REAL(kind=DBL) :: dummy !dummy variable for obtaining input. Must be real for readvar
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: ion_typer !dummy variable for ion_type array

    kc = 1
    ! READING INPUT ARRAYS FROM BINARY FILES

    ! p{1} Size of radial or scan arrays
    dimx = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{2} Size of wavenumber arrays
    dimn = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{3} Number of ions in system
    nions = INT(readvar(kc,dummy,ktype)) ; kc = kc+1
    ALLOCATE(ion_typer(dimx,nions))

    ! p{4} Flag for calculating decomposition of particle and heat transport into diffusive and convective components
    phys_meth = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{5} Flag for including collisions
    coll_flag = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{6} Rotation flag
    rot_flag = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{7} Number of total saught after solutions
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
    WHERE(ABS(Ate) < epsD) Ate = epsD

    ! p{28} R/Lne(rho) profile
    ALLOCATE(Ane(dimx))
    Ane = readvar(kc,Ane,ktype) ; kc = kc+1
    WHERE(ABS(Ane) < epsD) Ane = epsD

    ! p{29} Flag for adiabatic electrons
    el_type = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{30} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anise(dimx))
    anise = readvar(kc,anise,ktype) ; kc = kc+1

    ! p{31} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisedr(dimx))
    danisedr = readvar(kc,danisedr,ktype) ; kc = kc+1
    WHERE(ABS(danisedr) < epsD) danisedr = epsD

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
    WHERE(ABS(Ati) < epsD) Ati = epsD

    ! p{37} R/Lni(rho) profiles
    ALLOCATE(Ani(dimx,nions))
    Ani = readvar(kc,Ani,ktype) ; kc = kc+1
    WHERE(ABS(Ani) < epsD) Ani = epsD

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
    WHERE(ABS(danisdr) < epsD) danisdr = epsD

    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'

    ALLOCATE( oldrsol (dimx, dimn, numsols) )
    ALLOCATE( oldisol (dimx, dimn, numsols) )
    ALLOCATE( oldrfdsol (dimx, dimn, numsols) )
    ALLOCATE( oldifdsol (dimx, dimn, numsols) )
    ALLOCATE( oldsol (dimx, dimn, numsols) )
    ALLOCATE( oldfdsol (dimx, dimn, numsols) )

    OPEN(unit=myunit, file="output/primitive/rsol.dat", action="read", status="old")
    READ(myunit,fmtn) (((oldrsol(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/isol.dat", action="read", status="old")
    READ(myunit,fmtn) (((oldisol(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rfdsol.dat", action="read", status="old")
    READ(myunit,fmtn) (((oldrfdsol(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/ifdsol.dat", action="read", status="old")
    READ(myunit,fmtn) (((oldifdsol(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    oldsol = CMPLX(oldrsol,oldisol)
    oldfdsol = CMPLX(oldrfdsol,oldifdsol)

    DEALLOCATE( oldrsol )
    DEALLOCATE( oldisol )
    DEALLOCATE( oldrfdsol )
    DEALLOCATE( oldifdsol )

    relaccQL1=relacc1
    relaccQL2=relacc2

    lenwrk=(ndim+2+1)*(1+maxpts/(2**ndim+2*ndim*ndim+2*ndim+1)) !Set array size for 2D integration routine

    CALL MPI_Barrier(mpi_comm_world,ierror)

    !DEBUGGING WRITE OUT ALL INPUT TO ASCII FILE

    OPEN(unit=700, file="dimx.dat", action="write", status="replace")
    WRITE(700,'(I0)') dimx ; CLOSE(700)

    OPEN(unit=700, file="dimn.dat", action="write", status="replace")
    WRITE(700,'(I0)') dimn ; CLOSE(700)

    OPEN(unit=700, file="nions.dat", action="write", status="replace")
    WRITE(700,'(I0)') nions ; CLOSE(700)

    WRITE(fmtx, '(A)') '(G15.7)'
    !    WRITE(fmtxrow,'(A,I2,A)') '(',dimx,'G15.7)'
    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
    WRITE(fmtion,'(A,I0, A)') '(',nions,'G15.7)'
    WRITE(fmtintion,'(A,I0, A)') '(',nions,'I2)'
    WRITE(fmtecoef,'(A,I0, A)') '(',nions+1,'G15.7)'

    OPEN(unit=700, file="phys_meth.dat", action="write", status="replace")
    WRITE(700,'(I0)') phys_meth ; CLOSE(700)

    OPEN(unit=700, file="coll_flag.dat", action="write", status="replace")
    WRITE(700,'(I0)') coll_flag ; CLOSE(700)

    OPEN(unit=700, file="rot_flag.dat", action="write", status="replace")
    WRITE(700,'(I0)') rot_flag ; CLOSE(700)

    OPEN(unit=700, file="numsols.dat", action="write", status="replace")
    WRITE(700,'(I0)') numsols ; CLOSE(700)

    OPEN(unit=700, file="kthetarhos.dat", action="write", status="replace")
    WRITE(700,fmtx) (kthetarhos(j),j=1,dimn) ; CLOSE(700)

    OPEN(unit=700, file="x.dat", action="write", status="replace")
    WRITE(700,fmtx) (x(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ro.dat", action="write", status="replace")
    WRITE(700,fmtx) (Ro(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Rmin.dat", action="write", status="replace")
    WRITE(700,fmtx) (Rmin(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Bo.dat", action="write", status="replace")
    WRITE(700,fmtx) (Bo(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="qx.dat", action="write", status="replace")
    WRITE(700,fmtx) (qx(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="smag.dat", action="write", status="replace")
    WRITE(700,fmtx) (smag(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="alphax.dat", action="write", status="replace")
    WRITE(700,fmtx) (alphax(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Machtor.dat", action="write", status="replace")
    WRITE(700,fmtx) (Machtor(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Autor.dat", action="write", status="replace")
    WRITE(700,fmtx) (Autor(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Machpar.dat", action="write", status="replace")
    WRITE(700,fmtx) (Machpar(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Aupar.dat", action="write", status="replace")
    WRITE(700,fmtx) (Aupar(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="gammaE.dat", action="write", status="replace")
    WRITE(700,fmtx) (gammaE(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Tex.dat", action="write", status="replace")
    WRITE(700,fmtx) (Tex(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Nex.dat", action="write", status="replace")
    WRITE(700,fmtx) (Nex(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ate.dat", action="write", status="replace")
    WRITE(700,fmtx) (Ate(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ane.dat", action="write", status="replace")
    WRITE(700,fmtx) (Ane(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="el_type.dat", action="write", status="replace")
    WRITE(700,'(I0)') el_type ; CLOSE(700)

    OPEN(unit=700, file="Ai.dat", action="write", status="replace")
    WRITE(700,fmtion) ((Ai(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Zi.dat", action="write", status="replace")
    WRITE(700,fmtion) ((Zi(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Tix.dat", action="write", status="replace")
    WRITE(700,fmtion) ((Tix(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="ninorm.dat", action="write", status="replace")
    WRITE(700,fmtion) ((ninorm(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ati.dat", action="write", status="replace")
    WRITE(700,fmtion) ((Ati(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ani.dat", action="write", status="replace")
    WRITE(700,fmtion) ((Ani(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="ion_type.dat", action="write", status="replace")
    WRITE(700,fmtintion) ((ion_type(i,j),j=1,nions),i=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="maxpts.dat", action="write", status="replace")
    WRITE(700,'(I6)') maxpts ; CLOSE(700)

    OPEN(unit=700, file="maxruns.dat", action="write", status="replace")
    WRITE(700,'(I3)') maxruns ; CLOSE(700)

    OPEN(unit=700, file="relacc1.dat", action="write", status="replace")
    WRITE(700,fmtx) relacc1 ; CLOSE(700)

    OPEN(unit=700, file="relacc2.dat", action="write", status="replace")
    WRITE(700,fmtx) relacc2 ; CLOSE(700)

    ALLOCATE( modewidth (dimx, dimn) )
    ALLOCATE( distan (dimx, dimn) )
    ALLOCATE( solflu (dimx, dimn) )
    ALLOCATE( ntor (dimx, dimn) )
    ALLOCATE( modeflag (dimx) )
    ALLOCATE(phi(dimx,ntheta))
    ALLOCATE(npol(dimx,ntheta,nions))
    ALLOCATE(ecoefs(dimx,0:nions,numecoefs)) !includes electrons
    ALLOCATE(cftrans(dimx,nions,6)) ! only for ions

    ALLOCATE( sol (dimx, dimn, numsols) ) ; sol = oldsol
    ALLOCATE( fdsol (dimx, dimn, numsols) ) ; fdsol = oldfdsol
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

    IF (phys_meth /= 0) THEN
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
    DEALLOCATE(phi)
    DEALLOCATE(npol)
    DEALLOCATE(ecoefs)
    DEALLOCATE(cftrans)
    DEALLOCATE(sol)
    DEALLOCATE(fdsol)
    DEALLOCATE(Lcirce)
    DEALLOCATE(Lpiege)
    DEALLOCATE(Lecirce)
    DEALLOCATE(Lepiege)
    DEALLOCATE(Lcirci)
    DEALLOCATE(Lpiegi)
    DEALLOCATE(Lecirci)
    DEALLOCATE(Lepiegi)
    DEALLOCATE(Lvcirci)
    DEALLOCATE(Lvpiegi)
    IF (phys_meth /= 0) THEN
       DEALLOCATE(Lcircgte)
       DEALLOCATE(Lpieggte)
       DEALLOCATE(Lcircgne)
       DEALLOCATE(Lpieggne)
       DEALLOCATE(Lcircgue)
       DEALLOCATE(Lpieggue)

       DEALLOCATE(Lcircce)
       DEALLOCATE(Lpiegce)
       DEALLOCATE(Lcircgti)
       DEALLOCATE(Lpieggti)
       DEALLOCATE(Lcircgni)
       DEALLOCATE(Lpieggni)
       DEALLOCATE(Lcircgui)
       DEALLOCATE(Lpieggui)
       DEALLOCATE(Lcircci)
       DEALLOCATE(Lpiegci)
       IF (phys_meth == 2) THEN
          DEALLOCATE(Lecircgte)
          DEALLOCATE(Lepieggte)
          DEALLOCATE(Lecircgne)
          DEALLOCATE(Lepieggne)
          DEALLOCATE(Lecircgue)
          DEALLOCATE(Lepieggue)
          DEALLOCATE(Lecircce)
          DEALLOCATE(Lepiegce)
          DEALLOCATE(Lecircgti)
          DEALLOCATE(Lepieggti)
          DEALLOCATE(Lecircgni)
          DEALLOCATE(Lepieggni)
          DEALLOCATE(Lecircgui)
          DEALLOCATE(Lepieggui)
          DEALLOCATE(Lecircci)
          DEALLOCATE(Lepiegci)
       ENDIF
    ENDIF

    IF (ALLOCATED(oldsol)) DEALLOCATE(oldsol)
    IF (ALLOCATED(oldfdsol)) DEALLOCATE(oldfdsol)

  END SUBROUTINE deallocate_all

  SUBROUTINE outputascii()

    CHARACTER(len=20) :: fmtx,fmtn,fmtion,fmtxrow,fmtecoef,fmtcftrans,fmtecoefgau
    INTEGER :: i,j,k,l,myunit=700
    WRITE(fmtx, '(A)') '(G15.7)'
    WRITE(fmtxrow,'(A,I0,A)') '(',dimx,'G15.7)'
    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
    WRITE(fmtion,'(A,I0, A)') '(',nions,'G15.7)'
    WRITE(fmtecoef,'(A,I0, A)') '(',numecoefs,'G15.7)'
    WRITE(fmtcftrans,'(A)') '(6G15.7)'

    OPEN(unit=myunit, file="output/primitive/rsolflu.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((REAL(solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/isolflu.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((AIMAG(solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/distan.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((distan(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rmodewidth.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((REAL(modewidth(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rmodeshift.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((REAL(modeshift(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/imodewidth.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((AIMAG(modewidth(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/imodeshift.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((AIMAG(modeshift(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/ntor.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((ntor(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rsol.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(sol(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/isol.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(sol(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rfdsol.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(fdsol(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/ifdsol.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(fdsol(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLcirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lcirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLcirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lcirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLpiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lpiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLpiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lpiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLecirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lecirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLecirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lecirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLepiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lepiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLepiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lepiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLvcirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lvcirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLvcirce.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lvcirce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLvpiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((REAL(Lvpiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLvpiege.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((AIMAG(Lvpiege(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    IF (phys_meth /= 0) THEN
       OPEN(unit=myunit, file="output/primitive/rLcircgne.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lcircgne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgne.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lcircgne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggne.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lpieggne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggne.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lpieggne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircgue.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lcircgue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgue.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lcircgue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggue.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lpieggue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggue.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lpieggue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircgte.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lcircgte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgte.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lcircgte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggte.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lpieggte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggte.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lpieggte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircce.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lcircce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircce.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lcircce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpiegce.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((REAL(Lpiegce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpiegce.dat", action="write", status="replace")
       WRITE(myunit,fmtn) (((AIMAG(Lpiegce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
!!!
       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/primitive/rLecircgne.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lecircgne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgne.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lecircgne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggne.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lepieggne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggne.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lepieggne(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircgue.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lecircgue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgue.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lecircgue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggue.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lepieggue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggue.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lepieggue(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircgte.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lecircgte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgte.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lecircgte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggte.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lepieggte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggte.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lepieggte(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircce.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lecircce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircce.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lecircce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepiegce.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((REAL(Lepiegce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepiegce.dat", action="write", status="replace")
          WRITE(myunit,fmtn) (((AIMAG(Lepiegce(i,j,k)),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)
       ENDIF
    ENDIF

    OPEN(unit=myunit, file="output/primitive/rLcirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lcirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLcirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn)((((AIMAG(Lcirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLpiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lpiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLpiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((AIMAG(Lpiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLecirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lecirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLecirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((AIMAG(Lecirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLepiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lepiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLepiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((AIMAG(Lepiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLvcirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lvcirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLvcirci.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((AIMAG(Lvcirci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/rLvpiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((REAL(Lvpiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/iLvpiegi.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((((AIMAG(Lvpiegi(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

    IF (phys_meth /= 0) THEN
       OPEN(unit=myunit, file="output/primitive/rLcircgni.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lcircgni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgni.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lcircgni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggni.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lpieggni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggni.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lpieggni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircgui.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lcircgui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgui.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lcircgui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggui.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lpieggui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggui.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lpieggui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircgti.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lcircgti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircgti.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lcircgti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpieggti.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lpieggti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpieggti.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lpieggti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLcircci.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lcircci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLcircci.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lcircci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/rLpiegci.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((REAL(Lpiegci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/primitive/iLpiegci.dat", action="write", status="replace")
       WRITE(myunit,fmtn) ((((AIMAG(Lpiegci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
!!!
       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/primitive/rLecircgni.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lecircgni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgni.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lecircgni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggni.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lepieggni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggni.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lepieggni(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircgui.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lecircgui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgui.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lecircgui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggui.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lepieggui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggui.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lepieggui(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircgti.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lecircgti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircgti.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lecircgti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepieggti.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lepieggti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepieggti.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lepieggti(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLecircci.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lecircci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLecircci.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lecircci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/rLepiegci.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((REAL(Lepiegci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/primitive/iLepiegci.dat", action="write", status="replace")
          WRITE(myunit,fmtn) ((((AIMAG(Lepiegci(i,j,k,l)),j=1,dimn),i=1,dimx),k=1,nions),l=1,numsols) ; CLOSE(myunit)
       ENDIF
    ENDIF

    OPEN(unit=myunit, file="output/modeflag.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (modeflag(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/phi.dat", action="write", status="replace")
    WRITE(myunit,fmtxrow) ((phi(i,j),i=1,dimx),j=1,ntheta) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/npol.dat", action="write", status="replace")
    WRITE(myunit,fmtxrow) (((npol(i,j,k),i=1,dimx),j=1,ntheta),k=1,nions) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ecoefs.dat", action="write", status="replace")
    WRITE(myunit,fmtecoef) (((ecoefs(i,j,k),k=1,numecoefs),j=0,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/cftrans.dat", action="write", status="replace")
    WRITE(myunit,fmtcftrans) (((cftrans(i,j,k),k=1,6),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/gam_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((gam_GB(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ome_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ome_GB(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/gam_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((gam_SI(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ome_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ome_SI(i,j,k),j=1,dimn),i=1,dimx),k=1,numsols) ; CLOSE(myunit)

    IF (phys_meth /= 0) THEN
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

       OPEN(unit=myunit, file="output/vre_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (vre_SI(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vri_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((vri_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vce_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtx) (vce_SI(i),i=1,dimx) ; CLOSE(myunit)

       OPEN(unit=myunit, file="output/vci_SI.dat", action="write", status="replace")
       WRITE(myunit,fmtion) ((vci_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)
!!!
       IF (phys_meth == 2) THEN
          OPEN(unit=myunit, file="output/ceke.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (ceke(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/ceki.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((ceki(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/vene_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (vene_SI(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/veni_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((veni_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/vere_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtx) (vere_SI(i),i=1,dimx) ; CLOSE(myunit)

          OPEN(unit=myunit, file="output/veri_SI.dat", action="write", status="replace")
          WRITE(myunit,fmtion) ((veri_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

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

    OPEN(unit=myunit, file="output/evf_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (evf_SI(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ivf_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ivf_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (eef_GB(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ief_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ief_GB(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/evf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (evf_GB(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ivf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ivf_GB(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/epf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((epf_cm(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ipf_cm(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((eef_cm(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/ief_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ief_cm(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/evf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((evf_cm(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(myunit)

    OPEN(unit=myunit, file="output/ivf_cm.dat", action="write", status="replace")
    WRITE(myunit,fmtn) (((ivf_cm(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(myunit)

  END SUBROUTINE outputascii

  SUBROUTINE collectarraysredo()
    ! Collect all output into rank 0
    INTEGER :: ierr,myrank, i, nproc

    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: modewidthtmp,distantmp,FLRectmp,FLReptmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: gammatmp, Ladiatmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions) :: FLRiptmp, FLRictmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn,0:nions,0:9) :: ecoefsgautmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn) :: ommaxtmp, solflutmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lcircetmp, Lpiegetmp, Lecircetmp, Lepiegetmp, Lvcircetmp, Lvpiegetmp, Lcircgtetmp, Lpieggtetmp,  Lcircgnetmp, Lpieggnetmp, Lcircguetmp, Lpiegguetmp, Lcirccetmp, Lpiegcetmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: Lecircgtetmp, Lepieggtetmp, Lecircgnetmp, Lepieggnetmp, Lecircguetmp, Lepiegguetmp, Lecirccetmp, Lepiegcetmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lcircitmp, Lpiegitmp, Lecircitmp, Lepiegitmp, Lvcircitmp, Lvpiegitmp, Lcircgtitmp, Lpieggtitmp, Lcircgnitmp, Lpieggnitmp, Lcircguitmp, Lpiegguitmp, Lcirccitmp, Lpiegcitmp
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: Lecircgtitmp, Lepieggtitmp, Lecircgnitmp, Lepieggnitmp, Lecircguitmp, Lepiegguitmp, Lecirccitmp, Lepiegcitmp

    CALL mpi_comm_rank(mpi_comm_world,myrank,ierr)
    CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

    CALL MPI_Reduce(distan,distantmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(modewidth,modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(modeshift,modewidthtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRec,FLRectmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRep,FLReptmp,dimx*dimn,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(gamma,gammatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Ladia,Ladiatmp,dimx*dimn*numsols,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRic,FLRictmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(FLRip,FLRiptmp,dimx*dimn*nions,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(ommax,ommaxtmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(solflu,solflutmp,dimx*dimn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

    CALL MPI_Reduce(Lcirce,Lcircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lpiege,Lpiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lecirce,Lecircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lepiege,Lepiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvcirce,Lecircetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvpiege,Lepiegetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

    CALL MPI_Reduce(Lcirci,Lcircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lpiegi,Lpiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lecirci,Lecircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lepiegi,Lepiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvcirci,Lecircitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
    CALL MPI_Reduce(Lvpiegi,Lepiegitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

    CALL MPI_Reduce(ecoefsgau,ecoefsgautmp,dimx*dimn*(nions+1)*10,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)

    IF (phys_meth /= 0) THEN
       CALL MPI_Reduce(Lcircgte,Lcircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggte,Lpieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgne,Lcircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggne,Lpieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgue,Lcircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggue,Lpieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

       CALL MPI_Reduce(Lcircce,Lcirccetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpiegce,Lpiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

       CALL MPI_Reduce(Lcircgti,Lcircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggti,Lpieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgni,Lcircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggni,Lpieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lcircgui,Lcircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpieggui,Lpieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

       CALL MPI_Reduce(Lcircci,Lcirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       CALL MPI_Reduce(Lpiegci,Lpiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
!!!
       IF (phys_meth == 2) THEN
          CALL MPI_Reduce(Lecircgte,Lecircgtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggte,Lepieggtetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgne,Lecircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggne,Lepieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgue,Lecircgnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggue,Lepieggnetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

          CALL MPI_Reduce(Lecircce,Lecirccetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepiegce,Lepiegcetmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

          CALL MPI_Reduce(Lecircgti,Lecircgtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggti,Lepieggtitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgni,Lecircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggni,Lepieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lecircgui,Lecircgnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepieggui,Lepieggnitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)

          CALL MPI_Reduce(Lecircci,Lecirccitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
          CALL MPI_Reduce(Lepiegci,Lepiegcitmp,dimx*dimn*nions*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_world,ierr)
       ENDIF
    ENDIF

    IF (myrank == 0) THEN
       modewidth=modewidthtmp
       distan=distantmp
       FLRec=FLRectmp
       FLRep=FLReptmp
       gamma=gammatmp
       Ladia=Ladiatmp
       FLRip=FLRiptmp
       FLRic=FLRictmp
       ommax=ommaxtmp
       solflu=solflutmp
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

  END SUBROUTINE collectarraysredo


END PROGRAM qlk_redoQL
