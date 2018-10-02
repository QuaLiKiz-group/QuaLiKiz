

PROGRAM integration_driver

  USE kind
  USE diskio
  USE datcal
  USE datmat
  USE FLRterms !module contain functions defining all the FLR terms 
  USE mod_fluidsol !module which calculates the fluid growth rate and frequencies
  USE mod_contour !module which contains contour routines
  USE dispfuncs
  USE coleintegrals
  USE HCUB
  USE PCUB
  
  
  IMPLICIT NONE

  INTEGER  :: p, nu, ifailloc, contour
  COMPLEX(KIND=DBL) :: omega, omFkrtmp, out_
  REAL(KIND=DBL) :: theta, scale_, L_
  COMPLEX(KIND=DBL) :: varz, Co, Centre

  REAL(KIND=DBL), DIMENSION(ndim) :: aa, bb, cc, dd, xmin, xmax, IntegralValue, AbsErr
  REAL(KIND=DBL) :: NAGresultR, NAGresultI, NAGtemp, COLEresultR, COLEresultI, realtmp, reqepsrel, reqepsabs
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: wrkstr

  REAL(KIND=DBL)     :: acc, relerr, NAGtime, COLEtime, ee, ff
  INTEGER            :: minpts,npts, myunit, choice, error

  INTEGER :: time1, time2, freq  
  LOGICAL :: exist1, exist2, exist3, exist4 !used for checking for existence of files
  CHARACTER(len=20) :: myfmt, myint, fmtx,fmtn,fmtion,fmtintion,fmtxrow,fmt,coleformat
  CHARACTER(len = 60) :: coleformat2
  CHARACTER(len=:), ALLOCATABLE :: debugdir, outputdir, primitivedir, inputdir
  ! Old solution for non-reset runs
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldrsol,oldisol,oldrfdsol,oldifdsol
  !MPI variables:
  INTEGER :: ierror, nproc, myrank, i,j, write_primi,n, test, contour_loc, fast, do_NAG
  REAL(KIND=DBL), DIMENSION(5) :: fdata
  
  
  
  coleformat = '(F15.8, 4X, F15.8, )'
  coleformat2 = '(F15.8, 4X, F15.8, 4X, F15.8, )'

  CALL mpi_init(ierror)
  CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
  CALL mpi_comm_rank(mpi_comm_world,myrank,ierror)

  myunit = 700 !set a unique file unit per CPU
  p=1 
  nu=1

  pFkr=p   !for internals in integrand
  nuFkr=nu !for internals in integrand

  ! Read input and initialize all arrays
  CALL data_init()
  WRITE(*,*)
  WRITE(*,*) "Input data read!"
  CALL make_input()
  WRITE(*,*)
  WRITE(*,*) "Pre-processed input data!"
  CALL precalc()
  WRITE(*,*)
  WRITE(*,*) "Solved eigenfunction, FLR, processed solution, and ready for integration!"
  WRITE(*,*)
  WRITE(*,*) 
  WRITE(*,*)
  lenwrk = 1e5 !for NAG integration
  
  

  !Passing integral ranges. Need to multiply output by 4 due to symmetries
  aa(1) = 0.0d0 
  aa(2) = 0.0d0 
  bb(:) = rkuplim

  !Set integration limits for trapped electrons. (1) for kappa and (2) for v (for the electrons)
  cc(1) = 0.0d0
  cc(2) = 0.0d0
  dd(1) = 1.0d0 - barelyavoid
  dd(2) = vuplim
  
  !Set integration limits for trapped ions
  ee = 0._DBL
  ff = 1._DBL - barelyavoid
  
  xmin(1) = 0._DBL
  xmin(2) = 0._DBL
  
  xmax(1) = 1._DBL
  xmax(2) = 1._DBL
  
  reqepsrel = 0.08
  reqepsabs = 0.00

  fast = 1 !set to one to do it faster
  omFkr = sol(p,nu,1) !omega in integrands
  L_ = 2.5
  contour_loc = 1
  Co = MAX(AIMAG(solflu(p, nu)), 10.) * ci / 3.
  rint = ABS(REAL(solflu(p,nu))) / 5.
  Centre = Co + REAL(contour_loc) * (L_ * rint + 2. * (1 - 0.1))
  
  IF(contour_loc.EQ.0) THEN
    Centre = Co
    rint = 3.
  END IF

  contour = 1
  do_NAG = 1

  realtmp = REAL(omFkr)/100.
  
  
  ion = 1 !for integrand call
  call write_debug()
  ALLOCATE(wrkstr(lenwrk))
  

  open(unit=1,file="Hcubature_results_contour.txt",action="write",status="replace")
  open(unit=2,file="Pcubature_results_contour.txt",action="write",status="replace")
  open(unit=3,file="Mixcubature_results_contour.txt",action="write",status="replace")
  open(unit=9,file="omFkr_contour.txt",action="write",status="replace")
  open(unit=11,file="NAG_results_contour.txt",action="write",status="replace")
  open(unit=12,file="NAG_results_high_contour.txt",action="write",status="replace")
  
  omFkr = CMPLX(REAL(omFkr)/50., AIMAG(omFkr)/50.)
  omFkrtmp = omFkr
  


  ! Begin time measurement
  IF(do_NAG.EQ.1) THEN
  DO choice = 1,2
  !Scan over omFkr with a contour
  DO n = 1,500
    theta = 2.*pi*REAL(n-1)/REAL(500-1)
    varz = EXP(ci*theta)
    CALL squircle(Centre, varz, omFkr, rint)
    omFFk = omFkr

    out_ = Fkstarrstare(2, (/0., 0./), 1)
    scale_ = ABS(out_)
	  CALL SYSTEM_CLOCK(time1)
      DO i=1,100 - 90*fast
        minpts=0; ifailloc = 1; error = 0
        
        SELECT CASE(choice)
        CASE(1)
          !Compute to regular accuracy 
          NAGresultR = 0.
          NAGresultI = 0.
          ifailloc = 1
          NAGresultR = NAGresultR + d01ahf(ee, ff, relacc1, npts, relerr, rFFkiz, lw, ifailloc)
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          NAGresultI = NAGresultI + d01ahf(ee, ff, relacc1, npts, relerr, iFFkiz, lw, ifailloc)
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,cc,dd,minpts,maxpts,rFFke,relacc2,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultR = NAGresultR + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,cc,dd,minpts,maxpts,iFFke,relacc2,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultI = NAGresultI + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,aa,bb,minpts,maxpts,rFkstarrstar,relacc2,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultR = NAGresultR + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,aa,bb,minpts,maxpts,iFkstarrstar,relacc2,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultI = NAGresultI + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          
          NAGresultR = Ac(p) - NAGresultR
          NAGresultI = - NAGresultI
        CASE(2)
          !Compute to high accuracy 
          NAGresultR = 0.
          NAGresultI = 0.
          ifailloc = 1
          NAGresultR = NAGresultR + d01ahf(ee, ff, 0.00001_DBL, npts, relerr, rFFkiz, lw, ifailloc)
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          NAGresultI = NAGresultI + d01ahf(ee, ff, 0.00001_DBL, npts, relerr, iFFkiz, lw, ifailloc)
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,cc,dd,minpts,maxpts,rFFke,0.001_DBL,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultR = NAGresultR + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,cc,dd,minpts,maxpts,iFFke,0.001_DBL,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultI = NAGresultI + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,aa,bb,minpts,maxpts,rFkstarrstar,0.001_DBL,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultR = NAGresultR + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          CALL d01fcf(ndim,aa,bb,minpts,maxpts,iFkstarrstar,0.001_DBL,acc,lenwrk,wrkstr,NAGtemp ,ifailloc)
          NAGresultI = NAGresultI + NAGtemp
          IF(ifailloc /= 0) error = ifailloc
          ifailloc = 1
          
          NAGresultR = Ac(p) - NAGresultR
          NAGresultI = - NAGresultI
          EXIT !don't care how long it takes for high accuracy
        CASE DEFAULT
        
        END SELECT
     ENDDO
     IF (error /= 0) THEN
           IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I4,A,I4,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',error,&
                &'. Abnormal termination of NAG integration for n = ',n,' and choice = ',choice
           IF (ifailloc == -399) THEN
              WRITE(stderr,"(A)") 'NAG license error! Exiting'
              STOP
           ENDIF
        ENDIF
     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     NAGtime = REAL(time2-time1) / REAL(freq) * (1 + 9 * fast)
    IF(choice.EQ.1) WRITE(9,coleformat) REAL(omFkr), AIMAG(omFkr)
    WRITE(10+choice,coleformat2) NAGtime, NAGresultR, NAGresultI
  ENDDO
  
  IF(choice.eq.1) CLOSE(9)
  CLOSE(10+choice)  
  !WRITE(*,'(A,G14.7)') 'Adiabatic term = ',Ac
  END DO 
  END IF
  
  DO choice = 1,2
  !Choice here determines what method to use
  omFkr = omFkrtmp
  !Scan over omFkr with a contour
  DO n = 1,500
   theta = 2.*pi*REAL(n-1)/REAL(500-1)
   varz = EXP(ci*theta)
   CALL squircle(Centre, varz, omFkr, rint)
   omFFk = omFkr

   out_ = Fkstarrstar(2, (/0., 0./))
   out_ = Ac(p)
   scale_ = ABS(out_)
   fdata = (/scale_, bb(1), bb(2), dd(1), dd(2)/)
	 CALL SYSTEM_CLOCK(time1)
     DO i=1,100 - 90 * fast
        SELECT CASE(choice)
        CASE(1)
          COLEresultR = 0.
          COLEresultI = 0.
          test = hcubature(2, total_cubature, 2, xmin, xmax, maxpts, reqepsabs, reqepsrel, 2, IntegralValue, AbsErr, fdata = fdata)
          COLEresultR = COLEresultR +  IntegralValue(1) * scale_
          COLEresultI = COLEresultI + IntegralValue(2) * scale_
          
        
        CASE(2)
          COLEresultR = 0.
          COLEresultI = 0.
          test = pcubature(2, total_cubature, 2, xmin, xmax, maxpts, reqepsabs, reqepsrel, 2, IntegralValue, AbsErr, fdata = fdata)
          COLEresultR = COLEresultR +  IntegralValue(1) * scale_
          COLEresultI = COLEresultI + IntegralValue(2) * scale_
          
        CASE(3)
          aa(1) = 0.
          bb(1) = 5.
          COLEresultR = 0.
          COLEresultI = 0.
          IF(((contour_loc.EQ.-1).AND.(theta.GE.(7.*pi/4.))).OR.((contour_loc.EQ.-1).AND. &
            & ((theta.GE.(5.*pi/4.)).AND.(theta.LE.(3.*pi/2.))))) THEN
            test = pcubature(2, total_cubature, 2, xmin, xmax, maxpts, reqepsabs, reqepsrel, 2, IntegralValue, AbsErr, fdata = fdata)
          ELSE
            test = hcubature(2, total_cubature, 2, xmin, xmax, maxpts, reqepsabs, reqepsrel, 2, IntegralValue, AbsErr, fdata = fdata)
          END IF
          COLEresultR = IntegralValue(1)*scale_
          COLEresultI = IntegralValue(2)*scale_
          
          
        END SELECT
        !If choice = 1, compute to low accuracy, otherwise compute to high accuracy 
        !minpts=0; ifailloc = 1; error = 0
        !COLEresultR = 0.
        !COLEresultI = 0.    
        !CALL coleint(ndim,aa,bb,minpts,maxpts,Fkstarrstari,relacc2,acc,lenwrk,wrkstr,COLEresultR,COLEresultI,ifailloc,choice,error)
        !COLEresultR = COLEresultR * 4.
        !COLEresultI = COLEresultI * 4.
        !END IF
     END DO
     ! IF (error /= 0) THEN
           ! IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I4,A,I4,A,I3,A,G10.3,A,G10.3,A)") 'ifailloc = ',error,&
                ! &'. Abnormal termination of CUBPACK integration for j = ',j,' and kinty = ',n,' and choice = ',choice
           ! IF (ifailloc == -399) THEN
              ! WRITE(stderr,"(A)") 'NAG license error! Exiting'
              ! STOP
           ! ENDIF
        ! ENDIF
     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     COLEtime = REAL(time2-time1) / REAL(freq) * (1 + 9*fast)
        !WRITE(*,*)
        !WRITE(*,'(A,G14.7,A)') 'COLEtime for 100 real passing ion integrations =',COLEtime,' s'
        !WRITE(*,'(A,G14.7)') 'COLEresult for rLcirci =', COLEresult
        !WRITE(*,*)
        !IF (ifailloc /= 1) WRITE(*,*) 'ERROR, COLE FAILED TO CONVERGE'       
     WRITE(choice,coleformat2) COLEtime, COLEresultR, COLEresultI
  ENDDO
  CLOSE(choice)
  END DO
  
 
  DEALLOCATE(wrkstr)


  !Deallocating all
  CALL deallocate_all()
  

CONTAINS 

  SUBROUTINE data_init()
    !Parallel read data, allocate input and output arrays
    INTEGER, PARAMETER :: ktype = 1 ! BINARY FILES
    INTEGER :: ierr, fileno, verbosein, separatefluxin
    REAL(kind=DBL) :: dummy !dummy variable for obtaining input. Must be real for readvar

    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: dummyn
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: dummyx
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: dummyxnions
    REAL(kind=DBL), DIMENSION(:,:,:), ALLOCATABLE :: dummyxnnumsol

    ! READING INPUT ARRAYS FROM BINARY FILES

    inputdir = 'input/'

    dimx = 0
    dimx = INT(readvar(inputdir // 'dimx.bin', dummy, ktype, myunit))

    ! p{2} Size of wavenumber arrays
    dimn = 0
    dimn = INT(readvar(inputdir // 'dimn.bin', dummy, ktype, myunit))

    ! p{3} Number of ions in system
    nions = 0
    nions = INT(readvar(inputdir // 'nions.bin', dummy, ktype, myunit))

    ! p{9} Number of total saught after solutions
    numsols = 0
    numsols = INT(readvar(inputdir // 'numsols.bin', dummy, ktype, myunit))

    ! ALLOCATE TEMP ARRAYS
    ALLOCATE(dummyn(dimn))
    ALLOCATE(dummyx(dimx))
    ALLOCATE(dummyxnions(dimx,nions))
    ALLOCATE(dummyxnnumsol(dimx,dimn,numsols))

    ! p{12} Number of runs before runcounter resets
    maxruns = 0
    maxruns = INT(readvar(inputdir // 'maxruns.bin', dummy, ktype, myunit))

    ! p{4} Flag for calculating decomposition of particle and heat transport into diffusive and convective components
    phys_meth = 0
    phys_meth = INT(readvar(inputdir // 'phys_meth.bin', dummy, ktype, myunit))

    ! p{5} Flag for including collisions
    coll_flag = 0
    coll_flag = INT(readvar(inputdir // 'coll_flag.bin', dummy, ktype, myunit))

    write_primi = 0
    write_primi = INT(readvar(inputdir // 'write_primi.bin', dummy, ktype, myunit))

    ! p{6} Flag for including rotation
    rot_flag = 0
    rot_flag = INT(readvar(inputdir // 'rot_flag.bin', dummy, ktype, myunit))

    ! p{7} Flag for verbose output
    verbosein = 0
    verbosein = INT(readvar(inputdir // 'verbose.bin', dummy, ktype, myunit))

    ! p{8} Flag for separate mode flux output
    separatefluxin = 0
    separatefluxin = INT(readvar(inputdir // 'separateflux.bin', dummy, ktype, myunit))

    IF (verbosein > 0) verbose = .TRUE.
    IF (verbosein == 0) verbose = .FALSE.
    IF (separatefluxin == 1) separateflux = .TRUE.
    IF (separatefluxin == 0) separateflux = .FALSE.

    ! p{10} 1D integral accuracy
    relacc1 = 0
    relacc1 = readvar(inputdir // 'relacc1.bin', dummy, ktype, myunit)

    ! p{11} 2D integral accuracy
    relacc2 = 0
    relacc2 = readvar(inputdir // 'relacc2.bin', dummy, ktype, myunit)

    ! p{13} Maximum number of integrand evaluations in 2D integration routine
    maxpts = 0
    maxpts = INT(readvar(inputdir // 'maxpts.bin', dummy, ktype, myunit))

    ! p{14} Timeout seconds for a given solution search
    timeout = 0
    timeout = readvar(inputdir // 'timeout.bin', dummy, ktype, myunit)

    ! p{15} Multiplier for ETG saturation rule (default 1. Mostly for testing)
    ETGmult = 0
    ETGmult = readvar(inputdir // 'ETGmult.bin', dummy, ktype, myunit)

    ! p{16} Multiplier for collisionality (default 1. Mostly for testing)
    collmult = 0
    collmult = readvar(inputdir // 'collmult.bin', dummy, ktype, myunit)

    ! p{17} R0 geometric major radius (for normalizations)
    R0 = 0
    R0 = readvar(inputdir // 'R0.bin', dummy, ktype, myunit)

    ! p{18} Toroidal wave-number grid
    ALLOCATE(kthetarhos(dimn)); kthetarhos = 0 
    kthetarhos = readvar(inputdir // 'kthetarhos.bin', dummyn, ktype, myunit)

    ! p{19} Normalised radial coordinate (midplane radius)
    ALLOCATE(x(dimx)); x=0
    x = readvar(inputdir // 'x.bin', dummyx, ktype, myunit)

    ! p{20} Normalised radial coordinate (midplane radius)
    ALLOCATE(rho(dimx)); rho=0
    rho = readvar(inputdir // 'rho.bin', dummyx, ktype, myunit)

    ! p{21} <Ro> major radius
    ALLOCATE(Ro(dimx)); Ro=0
    Ro = readvar(inputdir // 'Ro.bin', dummyx, ktype, myunit)

    ! p{22} <a> minor radius
    ALLOCATE(Rmin(dimx)); Rmin=0
    Rmin = readvar(inputdir // 'Rmin.bin', dummyx, ktype, myunit)

    ! p{23} B(rho) magnetic field
    ALLOCATE(Bo(dimx)); Bo=0
    Bo = readvar(inputdir // 'Bo.bin', dummyx, ktype, myunit)

    ! p{24} q(rho) profile
    ALLOCATE(qx(dimx)); qx=0
    qx = readvar(inputdir // 'q.bin', dummyx, ktype, myunit)

    ! p{25} s(rho) profile
    ALLOCATE(smag(dimx)); smag=0
    smag = readvar(inputdir // 'smag.bin', dummyx, ktype, myunit)

    ! p{26} alpha(rho) profile
    ALLOCATE(alphax(dimx)); alphax=0
    alphax = readvar(inputdir // 'alpha.bin', dummyx, ktype, myunit)

    ! p{27} Machtor(rho) profile
    ALLOCATE(Machtor(dimx)); Machtor=0
    Machtor = readvar(inputdir // 'Machtor.bin', dummyx, ktype, myunit)

    ! p{28} Autor(rho) profile
    ALLOCATE(Autor(dimx)); Autor=0
    Autor = readvar(inputdir // 'Autor.bin', dummyx, ktype, myunit)
    WHERE(ABS(Autor) < epsD) Autor = epsD

    ! p{29} Machpar(rho) profile
    ALLOCATE(Machpar(dimx)); Machpar=0
    Machpar = readvar(inputdir // 'Machpar.bin', dummyx, ktype, myunit)

    ! p{30} Aupar(rho) profile
    ALLOCATE(Aupar(dimx)); Aupar=0
    Aupar = readvar(inputdir // 'Aupar.bin', dummyx, ktype, myunit)
    WHERE(ABS(Aupar) < epsD) Aupar = epsD

    ! p{31} gammaE(rho) profile
    ALLOCATE(gammaE(dimx)); gammaE=0
    gammaE = readvar(inputdir // 'gammaE.bin', dummyx, ktype, myunit)
    WHERE(ABS(gammaE) < epsD) gammaE = epsD

    ! p{32} Te(rho) profile
    ALLOCATE(Tex(dimx)); Tex=0
    Tex = readvar(inputdir // 'Te.bin', dummyx, ktype, myunit)

    ! p{33} ne(rho) profile
    ALLOCATE(Nex(dimx)); Nex=0
    Nex = readvar(inputdir // 'ne.bin', dummyx, ktype, myunit)

    ! p{34} R/LTe(rho) profile
    ALLOCATE(Ate(dimx)); Ate=0
    Ate = readvar(inputdir // 'Ate.bin', dummyx, ktype, myunit)
    WHERE(ABS(Ate) < epsD) Ate = epsD

    ! p{35} R/Lne(rho) profile
    ALLOCATE(Ane(dimx)); Ane=0
    Ane = readvar(inputdir // 'Ane.bin', dummyx, ktype, myunit)
    WHERE(ABS(Ane) < epsD) Ane = epsD
    WHERE(Ane+Ate < epsD) Ane = Ane+epsD

    ! p{36} Flag for adiabatic electrons
    el_type = 0
    el_type = INT(readvar(inputdir // 'typee.bin', dummy, ktype, myunit))

    ! p{37} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anise(dimx)); anise=0
    anise = readvar(inputdir // 'anise.bin', dummyx, ktype, myunit)

    ! p{38} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisedr(dimx)); danisedr=0
    danisedr = readvar(inputdir // 'danisdre.bin', dummyx, ktype, myunit)
    WHERE(ABS(danisedr) < epsD) danisedr = epsD

    ! p{39} Ti(rho) profiles
    ALLOCATE(Tix(dimx,nions)); Tix=0
    Tix = readvar(inputdir // 'Ti.bin', dummyxnions, ktype, myunit)

    ! p{40} ni/ne (rho) profiles
    ALLOCATE(ninorm(dimx,nions)); ninorm=0
    ninorm = readvar(inputdir // 'normni.bin', dummyxnions, ktype, myunit)

    ! p{41} R/LTi(rho) profiles
    ALLOCATE(Ati(dimx,nions)); Ati=0
    Ati = readvar(inputdir // 'Ati.bin', dummyxnions, ktype, myunit)
    WHERE(ABS(Ati) < epsD) Ati = epsD

    ! p{42} R/Lni(rho) profiles
    ALLOCATE(Ani(dimx,nions)); Ani=0
    Ani = readvar(inputdir // 'Ani.bin', dummyxnions, ktype, myunit)
    WHERE(ABS(Ani) < epsD) Ani = epsD
    WHERE(Ani+Ati < epsD) Ani = Ani+epsD

    ! p{43} Ion types
    ALLOCATE(ion_type(dimx,nions)); ion_type=0
    ion_type = INT(readvar(inputdir // 'typei.bin', dummyxnions, ktype, myunit))

    ! p{44} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anis(dimx,1:nions)); anis=0
    anis = readvar(inputdir // 'anisi.bin', dummyxnions, ktype, myunit)

    ! p{45} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisdr(dimx,1:nions)); danisdr=0
    danisdr = readvar(inputdir // 'danisdri.bin', dummyxnions, ktype, myunit)
    WHERE(ABS(danisdr) < epsD) danisdr = epsD

    ! p{46} Main ion mass
    ALLOCATE(Ai(dimx,nions)); Ai=0
    Ai = readvar(inputdir // 'Ai.bin', dummyxnions, ktype, myunit)

    ! p{47} Main ion charge
    ALLOCATE(Zi(dimx,nions)); Zi=0
    Zi = readvar(inputdir // 'Zi.bin', dummyxnions, ktype, myunit)

    INQUIRE(file="output/primitive/rsol.dat", EXIST=exist1)
    INQUIRE(file="output/primitive/isol.dat", EXIST=exist2)
    INQUIRE(file="output/primitive/rfdsol.dat", EXIST=exist3)
    INQUIRE(file="output/primitive/ifdsol.dat", EXIST=exist4)

    IF ( (exist1) .AND. (exist2) .AND. (exist3) .AND. (exist4) )THEN
       WRITE(*,*) 'Reading existing solution files in output/primitive'
    ELSE
       WRITE(*,*) 'No output solution files ready in output/primitive!'
       STOP
    END IF

    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G16.7)'

    ALLOCATE( oldrsol (dimx, dimn, numsols) ); oldrsol = 0
    ALLOCATE( oldisol (dimx, dimn, numsols) ); oldisol = 0
    ALLOCATE( oldrfdsol (dimx, dimn, numsols) ) ; oldrfdsol = 0
    ALLOCATE( oldifdsol (dimx, dimn, numsols) ); oldifdsol = 0
    ALLOCATE( sol (dimx, dimn, numsols) ); sol = 0
    ALLOCATE( fdsol (dimx, dimn, numsols) ); fdsol = 0

    primitivedir = "output/primitive/"
    myfmt = 'G16.7E3'
    oldrsol = readvar(primitivedir // 'rsol.dat', dummyxnnumsol, ktype, myunit)
    oldisol = readvar(primitivedir // 'isol.dat', dummyxnnumsol, ktype, myunit)
    oldrfdsol = readvar(primitivedir // 'rfdsol.dat', dummyxnnumsol, ktype, myunit)
    oldifdsol = readvar(primitivedir // 'ifdsol.dat', dummyxnnumsol, ktype, myunit)

    sol = CMPLX(oldrsol,oldisol)
    fdsol = CMPLX(oldrfdsol,oldifdsol)

    !DEBUGGING WRITE OUT ALL INPUT TO ASCII FILE
    !CALL write_debug()

    ALLOCATE( modewidth (dimx, dimn) )
    ALLOCATE( modeshift (dimx, dimn) )
    ALLOCATE( modeshift2 (dimx, dimn) )
    ALLOCATE( distan (dimx, dimn) )
    ALLOCATE(solflu (dimx, dimn))

    ALLOCATE( modeflag (dimx) )

    ALLOCATE(ana_solflu(dimx,dimn)); ana_solflu=0.
    ALLOCATE(jon_solflu(dimx,dimn)); jon_solflu=0.
    ALLOCATE(jon_modewidth(dimx,dimn)); jon_modewidth=0.
    ALLOCATE(jon_modeshift(dimx,dimn)); jon_modeshift=0.
    ALLOCATE( FLRec (dimx, dimn) ); FLRec=0.
    ALLOCATE( FLRep (dimx, dimn) ); FLRep=0.
    !Real 3D arrays with 3rd dimension equal to number of ions
    ALLOCATE( FLRic (dimx, dimn,nions) ); FLRic=0.
    ALLOCATE( FLRip (dimx, dimn,nions) ); FLRip=0.
    ALLOCATE( Lcirce (dimx, dimn, numsols) )
    ALLOCATE( Lpiege (dimx, dimn, numsols) )
    ALLOCATE( Lecirce (dimx, dimn, numsols) )
    ALLOCATE( Lepiege (dimx, dimn, numsols) )

    ALLOCATE( Lcirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lpiegi (dimx, dimn, nions, numsols) )
    ALLOCATE( Lecirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lepiegi (dimx, dimn, nions, numsols) )
    ALLOCATE( Lvcirci (dimx, dimn, nions, numsols) )
    ALLOCATE( Lvpiegi (dimx, dimn, nions, numsols) )

    ! DEALLOCATE TEMP ARRAYS
    DEALLOCATE(dummyn)
    DEALLOCATE(dummyx)
    DEALLOCATE(dummyxnions)
    DEALLOCATE(dummyxnnumsol)

  END SUBROUTINE data_init

  SUBROUTINE write_debug()

    INTEGER :: fileno

    myint='I15'
    myfmt='G16.7E3'
    debugdir='debugtest/'
    CALL writevar(debugdir // 'dimx.dat', dimx, myint, fileno, .FALSE.)
    CALL writevar(debugdir // 'dimn.dat', dimn, myint, fileno, .FALSE.)
    CALL writevar(debugdir // 'nions.dat', nions, myint, fileno)
    CALL writevar(debugdir // 'phys_meth.dat', phys_meth, myfmt, fileno)
    CALL writevar(debugdir // 'coll_flag.dat', coll_flag, myfmt, fileno)
    CALL writevar(debugdir // 'rot_flag.dat', rot_flag, myfmt, fileno)
    CALL writevar(debugdir // 'write_primi.dat', write_primi, myfmt, fileno)
    CALL writevar(debugdir // 'numsols.dat', numsols, myfmt, fileno)
    CALL writevar(debugdir // 'kthetarhos.dat', kthetarhos, myfmt, fileno)
    CALL writevar(debugdir // 'x.dat', x, myfmt, fileno)
    CALL writevar(debugdir // 'rho.dat', rho, myfmt, fileno)
    CALL writevar(debugdir // 'Ro.dat', Ro, myfmt, fileno)
    CALL writevar(debugdir // 'R0.dat', R0, myfmt, fileno)
    CALL writevar(debugdir // 'Rmin.dat', Rmin, myfmt, fileno)
    CALL writevar(debugdir // 'Bo.dat', Bo, myfmt, fileno)
    CALL writevar(debugdir // 'q.dat', qx, myfmt, fileno)
    CALL writevar(debugdir // 'smag.dat', smag, myfmt, fileno)
    CALL writevar(debugdir // 'alpha.dat', alphax, myfmt, fileno)
    CALL writevar(debugdir // 'Machtor.dat', Machtor, myfmt, fileno)
    CALL writevar(debugdir // 'Autor.dat', Autor, myfmt, fileno)
    CALL writevar(debugdir // 'Machpar.dat', Machpar, myfmt, fileno)
    CALL writevar(debugdir // 'Aupar.dat', Aupar, myfmt, fileno)
    CALL writevar(debugdir // 'gammaE.dat', gammaE, myfmt, fileno)
    CALL writevar(debugdir // 'Te.dat', Tex, myfmt, fileno)
    CALL writevar(debugdir // 'ne.dat', Nex, myfmt, fileno)
    CALL writevar(debugdir // 'Ate.dat', Ate, myfmt, fileno)
    CALL writevar(debugdir // 'Ane.dat', Ane, myfmt, fileno)
    CALL writevar(debugdir // 'typee.dat', el_type, myfmt, fileno)
    CALL writevar(debugdir // 'Ai.dat', Ai, myfmt, fileno)
    CALL writevar(debugdir // 'Zi.dat', Zi, myfmt, fileno)
    CALL writevar(debugdir // 'Ti.dat', Tix, myfmt, fileno)
    CALL writevar(debugdir // 'normni.dat', ninorm, myfmt, fileno)
    CALL writevar(debugdir // 'Ati.dat', Ati, myfmt, fileno)
    CALL writevar(debugdir // 'Ani.dat', Ani, myfmt, fileno)
    CALL writevar(debugdir // 'typei.dat', ion_type, myint, fileno)
    CALL writevar(debugdir // 'maxpts.dat', maxpts, myfmt, fileno)
    CALL writevar(debugdir // 'maxruns.dat', maxruns, myfmt, fileno)
    CALL writevar(debugdir // 'relacc1.dat', relacc1, myfmt, fileno)
    CALL writevar(debugdir // 'relacc2.dat', relacc2, myfmt, fileno)
    CALL writevar(debugdir // 'timeout.dat', timeout, myfmt, fileno)
    CALL writevar(debugdir // 'ETGmult.dat', ETGmult, myfmt, fileno)
    CALL writevar(debugdir // 'collmult.dat', collmult, myfmt, fileno)
    CALL writevar(debugdir // 'rsol_old.dat', oldrsol, myfmt, myunit)
    CALL writevar(debugdir // 'isol_old.dat', oldisol, myfmt, myunit)
    CALL writevar(debugdir // 'rfdsol_old.dat', oldrfdsol, myfmt, myunit)
    CALL writevar(debugdir // 'ifdsol_old.dat', oldifdsol, myfmt, myunit)
    CALL writevar(debugdir // 'rmwidth.dat', REAL(mwidth), myfmt, fileno)
    CALL writevar(debugdir // 'imwidth.dat', AIMAG(mwidth), myfmt, fileno)
    CALL writevar(debugdir // 'coefi.dat', coefi, myfmt, fileno)
    CALL writevar(debugdir // 'ktetaRhoi.dat', ktetaRhoi, myfmt, fileno)
    CALL writevar(debugdir // 'Athi.dat', Athi, myfmt, fileno)
    CALL writevar(debugdir // 'nwg.dat', nwg, myfmt, fileno)
    CALL writevar(debugdir // 'romFkr.dat', REAL(omFkr), myfmt, fileno)
    CALL writevar(debugdir // 'iomFkr.dat', AIMAG(omFkr), myfmt, fileno)
    CALL writevar(debugdir // 'd.dat', d, myfmt, fileno)
    CALL writevar(debugdir // 'mi.dat', mi, myfmt, fileno)
    CALL writevar(debugdir // 'qe.dat', qe, myfmt, fileno)
  END SUBROUTINE write_debug

  SUBROUTINE make_input()  

    ! List of input variables
    INTEGER:: p,nu !counter for loop over coordinates. p is radius (or general scan), nu is wavenumber
    INTEGER:: i,ilow,ihi !counter for loops
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Lambe, Nue !physical quantities used for the derivations, then discarded

    ! POPULATE GLOBAL VARIABLES

    !lenwrk=(ndim+2+1)*(1+maxptsin/(2**ndim+2*ndim*ndim+2*ndim+1)) !Set array size for 2D integration routine
    !Input array allocation

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
       coefi(:,i) = Zi(:,i)*Zi(:,i) * Nix(:,i) * tau(:,i)/ninorm(:,i) !Ion coefficients throughout equations. Normalised by ninorm to avoid small numbers. Later unnormalised
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


  SUBROUTINE deallocate_all()
    !DEALLOCATE INPUT ARRAYS
    DEALLOCATE(x)
    DEALLOCATE(rho)
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
    DEALLOCATE(modeshift2)    
    DEALLOCATE(modeflag)
    DEALLOCATE(Zeffx)
    DEALLOCATE(distan)
    DEALLOCATE(ntor)
    DEALLOCATE(solflu)
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
    DEALLOCATE(Lvpiegi)
    DEALLOCATE(Lvcirci)
    DEALLOCATE( oldrsol )
    DEALLOCATE( oldisol )
    DEALLOCATE( oldrfdsol )
    DEALLOCATE( oldifdsol )
    DEALLOCATE(tau)
    DEALLOCATE(Nix)
    DEALLOCATE(Nustar)
    DEALLOCATE(qprim)
    DEALLOCATE(Anue)
    DEALLOCATE(wg)
    DEALLOCATE(Ac)
    DEALLOCATE(csou)
    DEALLOCATE(cref) ! ref velocity used for input of gamma_E, U_par and \nabla U_par : sqrt(1keV/m_p) i.e. thermal velocity of D at 1keV
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
    DEALLOCATE(ft)
    DEALLOCATE(fc)
    DEALLOCATE(Machi)
    DEALLOCATE(Machitemp)
    DEALLOCATE(Mache)
    DEALLOCATE(Aui)
    DEALLOCATE(Aue)
    DEALLOCATE(tpernorme)
    DEALLOCATE(tpernormi)
    DEALLOCATE(tpernormifunc)
    DEALLOCATE(coefi)
    DEALLOCATE(mi)
    DEALLOCATE(ETG_flag)

    !These next ones are calculate at each p inside calcroutines
    DEALLOCATE(ana_solflu)
    DEALLOCATE(jon_solflu)
    DEALLOCATE(jon_modewidth)
    DEALLOCATE(jon_modeshift)
    DEALLOCATE( FLRec )
    DEALLOCATE( FLRep )
    DEALLOCATE( FLRic )
    DEALLOCATE( FLRip )
    DEALLOCATE(Athi)
    DEALLOCATE(ktetaRhoi)
    DEALLOCATE(Joi2)
    DEALLOCATE(Jobani2)
    DEALLOCATE(Joi2p)
    DEALLOCATE(J1i2p)
    DEALLOCATE(Joi2c)

  END SUBROUTINE deallocate_all

  SUBROUTINE precalc()

    ! Variables for fluid solution
    REAL(kind=DBL) :: ana_gamma

    ! Variables for defining contour locations
    REAL(kind=DBL) :: ma, wpi, wpe
    COMPLEX(kind=DBL) :: om0
    COMPLEX(kind=DBL) :: C, Centre
    REAL(kind=DBL) :: L

    ! Variables for collecting the solutions, and their status flags
    COMPLEX(kind=DBL), DIMENSION(numsols) :: soll, fdsoll, solltmp, fdsolltmp
    COMPLEX(kind=DBL) :: fdsollj, omin, fonout,newsol,newfdsol,soltest
    REAL(kind=DBL),    DIMENSION(numsols) :: isol,rsol,ifdsol,rfdsol,tmpsol
    INTEGER :: NN, i,j,npts, ifailloc,minlocind
    INTEGER, DIMENSION(1) :: minloci
    LOGICAL :: issol
    REAL(kind=DBL) :: kteta,maxdia
    REAL(KIND=DBL) :: maxklam,minklam,relerr

    CHARACTER(len=20) :: fmtn
    INTEGER :: myunit=700

    !INITIALIZATION OF VARIABLES************************
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

    !CALCULATION PHASE
    ! Set the ktheta dependent width tuning factors. The 0.075 prefactor was chosen on the base of optimizing a scan

    widthtuneITG = (0.075/kthetarhos(nu));
    widthtuneETG = (0.075/( kthetarhos(nu)*SQRT(me/mi(p,1)) ) );

    ! Calculates growth rate based on analytical fluid interchange / slab formulation. in units of nwg
    CALL ana_fluidsol(p,nu,ana_gamma)! with or without rotation ana_gamma used only for contour limit

    !diamagnetic frequencies used to set real part of "analytical" fluid solution (used for contour)
    wpi = wg(p)*Tix(p,1)/Zi(p,1)* (Ati(p,1) + Ani(p,1))
    wpe = wg(p)*Tex(p)/Ze * (Ate(p) + Ane(p)) 

    ! Set maximum absolute diamagnetic frequency for limit of contour search
    IF ( ABS(wpi) > ABS(wpe)) THEN
       maxdia = wpi
    ELSE
       maxdia=wpe
    ENDIF

    IF ( ETG_flag(nu) .EQV. .FALSE. ) THEN !ITG/TEM case
       ana_solflu(p,nu) = CMPLX(maxdia/wg(p),ana_gamma)
    ELSE !ETG case
       ana_solflu(p,nu) = CMPLX(wpe/wg(p),ana_gamma)
    ENDIF

    !Calculate mode width and shift according to Citrin model. omeflu, mwidth, and mshift are written inside the subroutines
    IF (ETG_flag(nu) .EQV. .TRUE.) THEN
       CALL jon_fluidsol_ele(p,nu)
    ELSEIF (rotflagarray(p) == 0) THEN 
       IF (rot_flag == 2) THEN
          CALL jon_fluidsol(p,nu) ! run in any case for QL momentum transport if rot_flag=2
          mwidth_rot=mwidth
          mshift_rot=mshift
       ENDIF
       CALL jon_fluidsol_norot(p,nu)
    ELSEIF (rotflagarray(p) ==1)  THEN
       IF (rot_flag == 2) THEN
          CALL jon_fluidsol(p,nu) ! run in any case for QL momentum transport if rot_flag=2
          mwidth_rot=mwidth
          mshift_rot=mshift
          Machi=Machimod; Aui=Auimod; gammaE=gammaEmod;
          CALL jon_fluidsol(p,nu) ! run with modified profiles
       ELSE
          CALL jon_fluidsol(p,nu)
       ENDIF
    ENDIF

    jon_solflu(p,nu) = omeflu / (ntor(p,nu)*wg(p))
    jon_modewidth(p,nu) = mwidth 
    jon_modeshift(p,nu) = mshift 

    !Write out fluid solution
    mshift = 0.0 

    solflu(p,nu) = ana_solflu(p,nu)

    !! choose which one to actually keep for use in the code
    modewidth(p,nu) = jon_modewidth(p,nu)
    mwidth = modewidth(p,nu)
    modeshift(p,nu) = jon_modeshift(p,nu)
    modeshift2(p,nu) = mshift2
    mshift = modeshift(p,nu); 

    !Calculate the width corresponding to |phi(x)| (in case eigenfunction is complex)
    widthhat = ABS(mwidth)**2 / SQRT(REAL(mwidth**2))

    !! Carry out FLR calculation with Bessel function integrations over kr
    IF (rotflagarray(p) == 1) THEN
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

    !used to pass variables in <Vpar^n> calculations below
    plam = p
    nulam = nu
    maxklam =   ABS(pi/(distan(p,nu)*normkr)) 
    minklam = - maxklam

    alamnorm = fc(p) !to be consistent with passing particle fraction

    alam1=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam1int,lw,ifailloc)/alamnorm !pitch angle average of sqrt(1-lambda*b)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam1 integration at p=',p,', nu=',nu
    ENDIF

    !pitch angle average of (1-lambda*b)
    alam2=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam2int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam2 integration at p=',p,', nu=',nu
    ENDIF

    alam3=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam3int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^3/2
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam3 integration at p=',p,', nu=',nu
    ENDIF

    alam4=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam4int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^2
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam4 integration at p=',p,', nu=',nu
    ENDIF

    alam5=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam5int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^5/2
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam5 integration at p=',p,', nu=',nu
    ENDIF

    !Set the transit frequency
    !    qRd = qx(p)*Ro(p)*d*SQRT(3._DBL)
    qRd = qx(p)*Ro(p)*d*1./alam1

    Athe=widthhat*cthe(p)/qRd 
    Athi(:)=widthhat*cthi(p,:)/qRd

    !Set the bounce frequency
    omega2bar = pi/2.*SQRT(x(p)*Rmin(p)/(2.*Ro(p)))

  END SUBROUTINE precalc

  REAL(KIND=DBL) FUNCTION alamnormint(lamin)
    !integrand for pinch angle iaverage of Vpar for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,2) !calculate transit time

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam1 solution T integration at p=',plam
    ENDIF

    alamnormint = 1./(4.*pi)*Tlam

  END FUNCTION alamnormint

  REAL(KIND=DBL) FUNCTION alam1int(lamin)
    !integrand for pinch angle iaverage of Vpar for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr,Rlam1,Rlam2
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint,Rint1, Rint2

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,3) 

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam1 solution Tlam integration at p=',plam
    ENDIF

    !Carry out transit average of Vpar
    Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**1. 
    CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,4) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam1 solution Rlam1 integration at p=',plam
    ENDIF
    CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,5) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam1 solution Rlam2 integration at p=',plam
    ENDIF
    alam1int = Rlam2/Rlam1*1./(4.*pi)*Tlam

  END FUNCTION alam1int

  REAL(KIND=DBL) FUNCTION alam2int(lamin)
    !integrand for pinch angle iaverage of Vpar^2 for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr,Rlam1,Rlam2
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint,Rint1, Rint2

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,6) 

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam2 solution Tlam integration at p=',plam
    ENDIF

    !Carry out transit average of Vpar^2
    Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**2. 
    CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,7) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam2 solution Rlam1 integration at p=',plam
    ENDIF
    CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,8) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam2 solution Rlam2 integration at p=',plam
    ENDIF
    alam2int = Rlam2/Rlam1*1./(4.*pi)*Tlam

  END FUNCTION alam2int

  REAL(KIND=DBL) FUNCTION alam3int(lamin)
    !integrand for pinch angle iaverage of Vpar^3 for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr,Rlam1,Rlam2
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint,Rint1, Rint2

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,9) 

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam3 solution Tlam integration at p=',plam
    ENDIF

    !Carry out transit average of Vpar^3
    Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**3. 
    CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,10) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam3 solution Rlam1 integration at p=',plam
    ENDIF
    CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,10) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam3 solution Rlam2 integration at p=',plam
    ENDIF
    alam3int = Rlam2/Rlam1*1./(4.*pi)*Tlam

  END FUNCTION alam3int

  REAL(KIND=DBL) FUNCTION alam4int(lamin)
    !integrand for pinch angle iaverage of Vpar^4 for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr,Rlam1,Rlam2
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint,Rint1, Rint2

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,11) 

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam4 solution Tlam integration at p=',plam
    ENDIF

    !Carry out transit average of Vpar^4
    Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**4. 
    CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,12) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam4 solution Rlam1 integration at p=',plam
    ENDIF
    CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,13) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam4 solution Rlam2 integration at p=',plam
    ENDIF
    alam4int = Rlam2/Rlam1*1./(4.*pi)*Tlam

  END FUNCTION alam4int

  REAL(KIND=DBL) FUNCTION alam5int(lamin)
    !integrand for pinch angle iaverage of Vpar^5 for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr,Rlam1,Rlam2
    INTEGER :: ifailloc2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint,Rint1, Rint2

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(plam)*SIN(theta/2.)**2))

    ifailloc2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifailloc2,14) 

    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam5 solution Tlam integration at p=',plam
    ENDIF

    !Carry out transit average of Vpar^5
    Rint1=1./SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))
    Rint2=Rint1 * SQRT(1.-lamin*(1.+2.*epsilon(plam)*SIN(theta/2.)**2))**5. 
    CALL davint(theta,Rint1,ntheta,0.,pi,Rlam1,ifailloc2,15) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam5 solution Rlam1 integration at p=',plam
    ENDIF
    CALL davint(theta,Rint2,ntheta,0.,pi,Rlam2,ifailloc2,16) 
    IF (ifailloc2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7)") 'ifailloc2 = ',ifailloc2,&
            &'. Abnormal termination of alam5 solution Rlam2 integration at p=',plam
    ENDIF
    alam5int = Rlam2/Rlam1*1./(4.*pi)*Tlam

  END FUNCTION alam5int
  
  
  
  INTEGER FUNCTION total_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    REAL(KIND=DBL) :: xx
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_, a, b, c, d
    INTEGER :: i
    
    
    scale_ = fdata(1) !scaling the integrand
    
    ! endpoint integration limits for passing
    a = fdata(2)
    b = fdata(3)
    ! endpoint integration limits for trapped
    c = fdata(4)
    d = fdata(5) 
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      total_cubature = 1
      RETURN
    END IF
    
    XY(1) = x(1) * a
    XY(2) = x(2) * b
    
    output = 4.*Fkstarrstare(ndim, XY, 1) * a * b
    DO i = 1, nions
      output = output + 4.* Fkstarrstari(ndim, XY, 1, i) * ninorm(pFkr, i) * a * b
    END DO
    
    XY(1) = x(1) * c
    XY(2) = x(2) * d
    xx = XY(1)
    
    output = output + FFke(ndim, XY, 1) * c * d
    DO i = 1, nions
      output = output + FFki(xx, 1, i) * ninorm(pFFk, i) * c
    END DO
        
    output = CMPLX(Ac(p), 0.) - output
    output = output/scale_
    
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    total_cubature = 0
    
    
  END FUNCTION total_cubature
  
  
  

  FUNCTION rFkstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstari
    rFkstarrstari = REAL  ( Fkstarrstari(ndim,xy,1,ion) )
  END FUNCTION rFkstarrstari
  

  FUNCTION iFkstarrstari(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the ion passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstari
    iFkstarrstari = AIMAG ( Fkstarrstari(ndim,xy,1,ion) )
  END FUNCTION iFkstarrstari

  FUNCTION rFkstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstare
    rFkstarrstare = REAL  ( Fkstarrstare(ndim,xy,1) )
  END FUNCTION rFkstarrstare

  FUNCTION iFkstarrstare(nf, xy)
    !---------------------------------------------------------------------
    ! Returns the electron passing particle integrand
    !---------------------------------------------------------------------  
    INTEGER, INTENT(IN) :: nf
    REAL(KIND=DBL), DIMENSION(nf), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstare
    iFkstarrstare = AIMAG (  Fkstarrstare(ndim,xy,1) )
  END FUNCTION iFkstarrstare

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
  
  INTEGER FUNCTION Fkstarrstari_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      Fkstarrstari_cubature = 1
      RETURN
    END IF
    XY(1) = x(1)
    XY(2) = x(2)
    output = Fkstarrstari(2, XY, 1, 1)/scale_
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    Fkstarrstari_cubature = 0
    
    
  END FUNCTION Fkstarrstari_cubature
  
  INTEGER FUNCTION Fkstarrstari2_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      Fkstarrstari2_cubature = 1
      RETURN
    END IF
    XY(1) = x(1)
    XY(2) = x(2)
    output = Fkstarrstari(2, XY, 1, 2)/scale_
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    Fkstarrstari2_cubature = 0
    
    
  END FUNCTION Fkstarrstari2_cubature
  
  

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

  
  INTEGER FUNCTION Fkstarrstare_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      Fkstarrstare_cubature = 1
      RETURN
    END IF
    XY(1) = x(1)
    XY(2) = x(2)
    output = Fkstarrstare(2, XY, 1)/scale_
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    Fkstarrstare_cubature = 0
    
    
  END FUNCTION Fkstarrstare_cubature
  
  


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
    INTEGER :: ifailloc

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
    ifailloc = 0;  Zgik = ci * sqrtpi * wofzweid(zik)

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
  
  INTEGER FUNCTION FFki_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.1.).OR.(fdim.NE.2)) THEN
      FFki_cubature = 1
      RETURN
    END IF
    
    
    XY = x(1)
    output = FFki(XY, 1, ion)/scale_
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    FFki_cubature = 0
    
    
  END FUNCTION FFki_cubature

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

  
  
  

  INTEGER FUNCTION FFke_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: output
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      FFke_cubature = 1
      RETURN
    END IF
    
    
    XY(1) = x(1)
    XY(2) = x(2)
    output = FFke(2, XY, 1)/scale_
    fval(1) = REAL(output)
    fval(2) = REAL(AIMAG(output))
    FFke_cubature = 0
    
    
  END FUNCTION FFke_cubature
  
  FUNCTION rFkstarrstar(ndim, xy)
    INTEGER :: i
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: xy
    REAL(KIND=DBL) :: rFkstarrstar
    REAL(KIND=DBL) :: intsum
    
    intsum = REAL(Fkstarrstare(ndim,xy,1))
    DO i = 1,nions
      intsum = intsum + REAL(Fkstarrstari(ndim,xy,1,i))*ninorm(pFkr,i)
    END DO
    rFkstarrstar = 4.*intsum

    
  END FUNCTION
  
  FUNCTION iFkstarrstar(ndim, xy)
    INTEGER :: i
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: xy
    REAL(KIND=DBL) :: iFkstarrstar
    REAL(KIND=DBL) :: intsum
    
    intsum = AIMAG(Fkstarrstare(ndim,xy,1))
    DO i = 1,nions
      intsum = intsum + AIMAG(Fkstarrstari(ndim,xy,1,i))*ninorm(pFkr,i)
    END DO
    iFkstarrstar = 4.*intsum

    
  END FUNCTION
  
  COMPLEX(KIND=DBL) FUNCTION Fkstarrstar(ndim, XY)
    INTEGER :: i
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: xy
    COMPLEX(KIND=DBL) :: intsum
    
    intsum = Fkstarrstare(ndim,xy,1)
    DO i = 1,nions
      intsum = intsum + Fkstarrstari(ndim,xy,1,i)*ninorm(pFkr,i)
    END DO
    Fkstarrstar = 4.*intsum

    
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION rFFkiz(kk)
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
      intsum = intsum+REAL ( FFki(kk,1,i) )*ninorm(pFFk,i)
    END DO
    rFFkiz=intsum
  END FUNCTION rFFkiz
  
  REAL(KIND=DBL) FUNCTION iFFkiz(kk)
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL) :: intsum
    INTEGER :: i
    intsum=0
    !Due to the 1D real integral, we have to separate the real and imaginary parts
    DO i = 1,nions
      intsum = intsum+ AIMAG ( FFki(kk,1,i) )*ninorm(pFFk,i)
    END DO
    iFFkiz=intsum
  END FUNCTION iFFkiz
  
  
  
  INTEGER FUNCTION Fkstarrstar_cubature(ndim, x, fdim, fdata, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: intsum
    INTEGER :: i
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      Fkstarrstar_cubature = 1
      RETURN
    END IF
    XY(1) = x(1)
    XY(2) = x(2)
    intsum = Fkstarrstare(ndim, XY, 1)
    DO i = 1,2
      intsum = intsum + Fkstarrstari(ndim,xy,1,i)*ninorm(pFkr,i)
    END DO
    intsum = intsum/scale_
    
    fval(1) = REAL(intsum)
    fval(2) = REAL(AIMAG(intsum))
    Fkstarrstar_cubature = 0
  
  END FUNCTION
  
  INTEGER FUNCTION Fkstarrstaritot_cubature(ndim, x, fdata, fdim, fval)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim

      
    REAL(KIND=DBL), DIMENSION(2) :: XY
    COMPLEX(KIND=DBL) :: intsum
    INTEGER :: i
    REAL(KIND=DBL) :: scale_
    
    scale_ = fdata(1)
      
    IF((ndim.NE.2.).OR.(fdim.NE.2)) THEN
      Fkstarrstaritot_cubature = 1
      RETURN
    END IF
    XY(1) = x(1)
    XY(2) = x(2)
    intsum = 0._DBL
    DO i = 1,2
      intsum = intsum + Fkstarrstari(ndim,xy,1,i)*ninorm(pFkr,i)
    END DO
    intsum = intsum/scale_
    fval(1) = REAL(intsum)
    fval(2) = REAL(AIMAG(intsum))
    Fkstarrstaritot_cubature = 0
  
  END FUNCTION
  
  
  

END PROGRAM integration_driver
