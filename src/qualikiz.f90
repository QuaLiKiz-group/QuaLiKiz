! ------------------------------------------------------------------------------------!
! QuaLiKiz - quasilinear gyrokinetic calculation of growth rates and transport fluxes !
!                                                                                     !
! Clarisse Bourdelle, Chantal Passeron, 			                      !
! Claude Fourment, Giorgio Regnoli,                                                   !
! Alessandro Casati, Pierre Cottier, Jonathan Citrin 				      !	
!								                      ! 
! PARALLEL VERSION                                                                    !
!-------------------------------------------------------------------------------------!

SUBROUTINE qualikiz(dimxin, rhoin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, kthetarhosin, & !general param
     & xin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & !geometry input
     & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & !electron input
     & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & !ion input
     & Machtorin, Autorin, Machparin, Auparin, gammaEin, & !rotation input
     & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, & !code specific input
     & epf_SIout,eef_SIout,evf_SIout,ipf_SIout,ief_SIout,ivf_SIout, & ! Non optional outputs
     & solflu_SIout, solflu_GBout, gam_SIout,gam_GBout,ome_SIout,ome_GBout, & !growth rate and frequency output
     & epf_GBout,eef_GBout, evf_GBout, dfe_SIout,vte_SIout,vre_SIout,vce_SIout,epf_cmout,eef_cmout,evf_cmout,ckeout, & !electron flux outputs
     & ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout,vti_SIout,vri_SIout,vci_SIout,ipf_cmout,ief_cmout,ivf_cmout,ckiout, & !ion flux outputs
     & dfe_GBout,vte_GBout,vre_GBout,vce_GBout,dfi_GBout,vti_GBout,vri_GBout,vci_GBout, &
     & vene_SIout,chiee_SIout,vere_SIout,vece_SIout, cekeout, veni_SIout,chiei_SIout,veci_SIout,veri_SIout,cekiout, & !heat pinch outputs
     & modeflagout, Nustarout, Zeffxout, & 
     & phiout, npolout, ecoefsout, cftransout, &  ! poloidal asymmetry outputs for heavy impurities
     & solfluout, modewidthout, modeshiftout, distanout, ntorout, solout, fdsolout,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
     & kperp2out,krmmuITGout,krmmuETGout,&
     & Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lvcirceout, Lvpiegeout, Lcircgteout, Lpieggteout, Lcircgneout, Lpieggneout, Lcircgueout, Lpieggueout, Lcircceout, Lpiegceout, & 
     & Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout, &
     & Lecircgteout, Lepieggteout, Lecircgneout, Lepieggneout, Lecircgueout, Lepieggueout, Lecircceout, &
     &Lepiegceout, Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout,&
     oldsolin, oldfdsolin, runcounterin,&
     rhominin,rhomaxin,&
     & eefETG_SIout,eefETG_GBout,&  !optional outputs from separation of fluxes
     & eefTEM_SIout,eefTEM_GBout,epfTEM_SIout,dfeTEM_SIout,vteTEM_SIout,vceTEM_SIout,vreTEM_SIout,dfeTEM_GBout,vteTEM_GBout,vceTEM_GBout,vreTEM_GBout,&
     & eefITG_SIout,eefITG_GBout,epfITG_SIout,dfeITG_SIout,vteITG_SIout,vceITG_SIout,vreITG_SIout,dfeITG_GBout,vteITG_GBout,vceITG_GBout,vreITG_GBout,&
     & iefTEM_SIout,iefTEM_GBout,ipfTEM_SIout,ivfTEM_SIout,ivfTEM_GBout,dfiTEM_SIout,vtiTEM_SIout,vciTEM_SIout,vriTEM_SIout,dfiTEM_GBout,vtiTEM_GBout,vciTEM_GBout,vriTEM_GBout,&
     & iefITG_SIout,iefITG_GBout,ipfITG_SIout,ivfITG_SIout,ivfITG_GBout,dfiITG_SIout,vtiITG_SIout,vciITG_SIout,vriITG_SIout,dfiITG_GBout,vtiITG_GBout,vciITG_GBout,vriITG_GBout)

  !BRIEF EXPLANATION OF MODULES
  !
  !kind defines the kind numbers and output/error units
  !mod_io_managment contains all the subroutines and functions used for file IO
  !mod_make_io contains all routines for (de)allocating, defining, and writing input and output data
  !datmat contains all shared data 
  !datcal contains all shared parameters
  !calcroutines contains all the main calculation routines
  !asymmetry calculates the poloidally dependent quantities relevant for heavy impurities
  !saturation calculates the nonlinear saturation rule and makes the final output

  USE kind
  !USE mod_io_management !MPI is included in this module
  USE mod_make_io
  USE datmat
  USE datcal
  USE calcroutines  
  USE asymmetry
  USE mod_saturation

  IMPLICIT NONE

  !EXTERNAL phieq

  !Local data dictionary.
  !
  !Time measuring variables
  REAL(kind=DBL) :: cputime1, cputime2, tpstot, timetot
  INTEGER :: time1, time2, time3, time4, freq, cputimetot
  CHARACTER(len=20) :: myfmt
  !MPI variables:
  !TotTask: Total number of Tasks (radial*wavenumber coordinates)
  INTEGER,DIMENSION(MPI_STATUS_SIZE) :: status
  INTEGER, DIMENSION(:), ALLOCATABLE :: wavenum,radcoord

  INTEGER :: i,j,k,ierror,nproc,irank,TotTask, myrank

  ! List of input variables
  INTEGER, INTENT(IN) :: dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin,verbosein,separatefluxin, el_typein
  INTEGER, DIMENSION(dimxin,nionsin), INTENT(IN) :: ion_typein
  REAL(kind=DBL), INTENT(IN) :: R0in
  REAL(kind=DBL), DIMENSION(dimxin), INTENT(IN) :: xin, rhoin, Roin, Rminin, Boin, qxin, smagin, alphaxin
  REAL(kind=DBL), DIMENSION(dimnin), INTENT(IN) :: kthetarhosin
  REAL(kind=DBL), DIMENSION(dimxin), INTENT(IN) :: Texin, Nexin, Atein, Anein, anisein, danisedrin
  REAL(kind=DBL), DIMENSION(dimxin,nionsin), INTENT(IN) :: Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, Aiin, Ziin
  REAL(kind=DBL), DIMENSION(dimxin), INTENT(IN) :: Auparin, gammaEin, Machtorin, Machparin, Autorin
  INTEGER, INTENT(IN) :: maxrunsin, maxptsin
  REAL(kind=DBL), INTENT(IN) :: relacc1in, relacc2in, timeoutin
  REAL(kind=DBL), OPTIONAL, INTENT(IN) :: rhominin,rhomaxin

  ! List of output variables: 

  ! growth rate and frequency outputs
  COMPLEX(KIND=DBL), DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT) :: solflu_SIout, solflu_GBout, solfluout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: gam_SIout,gam_GBout,ome_SIout,ome_GBout  

  ! final output arrays following saturation rule
  REAL(KIND=DBL), DIMENSION(dimxin), INTENT(OUT)  :: epf_SIout,eef_SIout,evf_SIout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), INTENT(OUT)  :: ipf_SIout,ief_SIout,ivf_SIout


  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefETG_SIout,eefETG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefTEM_SIout,eefTEM_GBout,epfTEM_SIout,dfeTEM_SIout,vteTEM_SIout,vceTEM_SIout,vreTEM_SIout,dfeTEM_GBout,vteTEM_GBout,vceTEM_GBout,vreTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefITG_SIout,eefITG_GBout,epfITG_SIout,dfeITG_SIout,vteITG_SIout,vceITG_SIout,vreITG_SIout,dfeITG_GBout,vteITG_GBout,vceITG_GBout,vreITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefTEM_SIout,ipfTEM_SIout,dfiTEM_SIout,vtiTEM_SIout,vciTEM_SIout,dfiTEM_GBout,vtiTEM_GBout,vciTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: vriTEM_SIout,vriTEM_GBout,ivfTEM_SIout,iefTEM_GBout,ivfTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefITG_SIout,ipfITG_SIout,dfiITG_SIout,vtiITG_SIout,vciITG_SIout,vriITG_SIout,dfiITG_GBout,vtiITG_GBout,vciITG_GBout,vriITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: ivfITG_SIout,ivfITG_GBout,iefITG_GBout

  REAL(KIND=DBL), DIMENSION(dimxin) :: eefETG_SIouttmp,eefETG_GBouttmp
  REAL(KIND=DBL), DIMENSION(dimxin) :: eefTEM_SIouttmp,epfTEM_SIouttmp,eefTEM_GBouttmp,dfeTEM_SIouttmp,vteTEM_SIouttmp,vceTEM_SIouttmp,vreTEM_SIouttmp,dfeTEM_GBouttmp,vteTEM_GBouttmp,vceTEM_GBouttmp,vreTEM_GBouttmp
  REAL(KIND=DBL), DIMENSION(dimxin) :: eefITG_SIouttmp,epfITG_SIouttmp,eefITG_GBouttmp,dfeITG_SIouttmp,vteITG_SIouttmp,vceITG_SIouttmp,vreITG_SIouttmp,dfeITG_GBouttmp,vteITG_GBouttmp,vceITG_GBouttmp,vreITG_GBouttmp

  REAL(KIND=DBL), DIMENSION(dimxin,nionsin) :: iefTEM_SIouttmp,iefTEM_GBouttmp,ipfTEM_SIouttmp,dfiTEM_SIouttmp,vtiTEM_SIouttmp,vciTEM_SIouttmp,dfiTEM_GBouttmp,vtiTEM_GBouttmp,vciTEM_GBouttmp
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin) :: vriTEM_SIouttmp,vriTEM_GBouttmp,ivfTEM_SIouttmp,ivfTEM_GBouttmp
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin) :: iefITG_SIouttmp,iefITG_GBouttmp,ipfITG_SIouttmp,dfiITG_SIouttmp,vtiITG_SIouttmp,vciITG_SIouttmp,vriITG_SIouttmp,dfiITG_GBouttmp,vtiITG_GBouttmp,vciITG_GBouttmp,vriITG_GBouttmp
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin) :: ivfITG_SIouttmp,ivfITG_GBouttmp


  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: epf_GBout,eef_GBout, evf_GBout, dfe_SIout, vte_SIout, vre_SIout, vce_SIout, dfe_GBout, vte_GBout, vre_GBout, vce_GBout, ckeout, modeflagout, Nustarout, Zeffxout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: vene_SIout, chiee_SIout, vere_SIout, vece_SIout, cekeout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout, vti_SIout, vri_SIout, vci_SIout,dfi_GBout, vti_GBout, vri_GBout, vci_GBout, ckiout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: veni_SIout, veri_SIout, chiei_SIout, veci_SIout, cekiout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  ::  epf_cmout, eef_cmout, evf_cmout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,nionsin), OPTIONAL, INTENT(OUT) :: ipf_cmout,ief_cmout, ivf_cmout
  REAL(KIND=DBL), DIMENSION(dimxin,ntheta), OPTIONAL, INTENT(OUT)  ::  phiout
  REAL(KIND=DBL), DIMENSION(dimxin,ntheta,nionsin), OPTIONAL, INTENT(OUT)  ::  npolout
  REAL(KIND=DBL), DIMENSION(dimxin,0:nionsin,numecoefs), OPTIONAL, INTENT(OUT)  ::  ecoefsout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin,numicoefs), OPTIONAL, INTENT(OUT)  ::  cftransout

  ! optional output arrays from which the saturation rule can be calculated without rerunning dispersion relation solver
  REAL(KIND=DBL) , DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: krmmuITGout,krmmuETGout
  REAL(KIND=DBL) , DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: distanout,ntorout,kperp2out
  COMPLEX(KIND=DBL), DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: modewidthout, modeshiftout
  COMPLEX(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: solout, fdsolout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lvcirceout, Lvpiegeout, Lcircgteout, Lpieggteout,  Lcircgneout, Lpieggneout,  Lcircgueout, Lpieggueout, Lcircceout, Lpiegceout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgteout, Lepieggteout,  Lecircgneout, Lepieggneout,  Lecircgueout, Lepieggueout, Lecircceout, Lepiegceout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: &
       Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout

  ! optional input arrays for going directly to newton solver
  INTEGER, OPTIONAL, INTENT(IN)  :: runcounterin
  COMPLEX(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(IN)  :: oldsolin,oldfdsolin
  INTEGER :: myunit=700
  !__________________END OF LOCAL DATA DICTIONARY__________________________

  !DEBUGGING
  CHARACTER(len=20) :: fmtn

  !WRITE(stdout,'(1X,1A,I1,/)') 'Integration scheme type = ',inttype

  ! -- MPI Initialisation -- !
  CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
  CALL mpi_comm_rank(mpi_comm_world,myrank,ierror)

  CALL SYSTEM_CLOCK(time1)

  ! Make the input (including derived quantities)
  CALL make_input(dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, kthetarhosin, & !general param
       & xin, rhoin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & 
       & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & 
       & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & 
       & Machtorin, Autorin, Machparin, Auparin, gammaEin, &
       & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin)  !code specific inputs

  ! set optional input
  IF (PRESENT(oldsolin)) THEN
     ALLOCATE(oldsol(dimxin,dimnin,numsolsin))
     oldsol = oldsolin
  ENDIF
  IF (PRESENT(oldfdsolin)) THEN
     ALLOCATE(oldfdsol(dimxin,dimnin,numsolsin))
     oldfdsol = oldfdsolin
  ENDIF
  IF (PRESENT(runcounterin)) THEN
     runcounter=runcounterin
  ELSE
     runcounter=0 
  ENDIF
  IF (PRESENT(rhominin)) THEN
     rhomin=rhominin
  ELSE
     rhomin=0.00 
  ENDIF
  IF (PRESENT(rhomaxin)) THEN
     rhomax=rhomaxin
  ELSE
     rhomax=1.0
  ENDIF

  !Check sanity of input (these can be much expanded)
  IF ( (onlyion .EQV. .TRUE.) .AND. (onlyelec .EQV. .TRUE.) ) THEN
     WRITE(stderr,*) 'onlyion and onlyelec both set to true! Abandon ship...'
     CALL mpi_abort(mpi_comm_world,-1)
  ENDIF

  IF ( (rot_flag < 0) .OR. (rot_flag > 2) ) THEN
     WRITE(stderr,*) 'rot_flag must be between 0 and 2! Abandon ship...'
     CALL mpi_abort(mpi_comm_world,-1)
  ENDIF


  !Allocation and initialization of calculated arrays (named "output")
  CALL allocate_output() !subroutine found in mod_make_io
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
     IF(ALLOCATED(wavenum) .EQV. .FALSE.) ALLOCATE(wavenum(TotTask))
     IF(ALLOCATED(radcoord) .EQV. .FALSE.) ALLOCATE(radcoord(TotTask))
  END IF

  CALL MPI_Barrier(mpi_comm_world,ierror)

  !! NOW THE MAGIC HAPPENS!! This subroutine is contained below
  !! Distributes tasks to all processors and calculates output

  CALL DistriTask(TotTask,nproc,myrank)

  IF (myrank==0) THEN
     IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A)") '*** Collecting output'
  ENDIF

  !Consolidate all arrays into rank=0 for QL flux integrals and output
  IF (nproc > 1) THEN
     IF (myrank==0) CALL SYSTEM_CLOCK(time3)
     CALL collectarrays()
  ENDIF
  CALL MPI_Barrier(mpi_comm_world,ierror)

  CALL allocate_endoutput()

  !If rank0, then carry out the saturation rules and output final results. This will soon be parallelized too. Trivial over dimx
  IF (myrank==0) THEN 
     CALL SYSTEM_CLOCK(time4)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time4-time3) / REAL(freq)
     WRITE(stdout,*)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: First MPI_AllReduce time = ',timetot,' s'  !final write

     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time2-time1) / REAL(freq)
     WRITE(stdout,*)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: Hurrah! All eigenmodes calculated! Time = ',timetot,' s'  !final write

     IF (verbose .EQV. .TRUE.) WRITE(stdout,*)     
     IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A)") '*** Calculating nonlinear saturation rule'

     CALL SYSTEM_CLOCK(time1)
  ENDIF
  IF (separateflux .EQV. .TRUE.) THEN
     IF (myrank==0) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A)") '*** separateflux=T ,  NL saturation rule for separate modes also calculated'
        IF (verbose .EQV. .TRUE.) WRITE(stdout,*)     
     ENDIF
     CALL saturation(1) !set 0 for including all modes, 1 for only ITG, 2 for only TEM, 3 for only ETG
     IF (PRESENT(eefITG_SIout))   eefITG_SIout=eef_SI; 
     IF (PRESENT(eefITG_GBout))   eefITG_GBout=eef_GB; 
     IF (PRESENT(epfITG_SIout))   epfITG_SIout=epf_SI; 

     IF (PRESENT(dfeITG_SIout))   dfeITG_SIout=dfe_SI; 
     IF (PRESENT(vteITG_SIout))  vteITG_SIout=vte_SI; 
     IF (PRESENT(vceITG_SIout))  vceITG_SIout=vce_SI; 
     IF (PRESENT(vreITG_SIout))  vreITG_SIout=vre_SI; 

     IF (PRESENT(dfeITG_GBout))   dfeITG_GBout=dfe_GB; 
     IF (PRESENT(vteITG_GBout))  vteITG_GBout=vte_GB; 
     IF (PRESENT(vceITG_GBout))  vceITG_GBout=vce_GB; 
     IF (PRESENT(vreITG_GBout))  vreITG_GBout=vre_GB; 

     IF (PRESENT(iefITG_SIout))   iefITG_SIout=ief_SI; 
     IF (PRESENT(ipfITG_SIout))   ipfITG_SIout=ipf_SI; 
     IF (PRESENT(ivfITG_SIout))   ivfITG_SIout=ivf_SI; 

     IF (PRESENT(iefITG_GBout))   iefITG_GBout=ief_GB; 
     IF (PRESENT(ivfITG_GBout))   ivfITG_GBout=ivf_GB; 

     IF (PRESENT(dfiITG_SIout))   dfiITG_SIout=dfi_SI; 
     IF (PRESENT(vtiITG_SIout))  vtiITG_SIout=vti_SI; 
     IF (PRESENT(vciITG_SIout))  vciITG_SIout=vci_SI; 
     IF (PRESENT(vriITG_SIout))  vriITG_SIout=vri_SI; 

     IF (PRESENT(dfiITG_GBout))   dfiITG_GBout=dfi_GB; 
     IF (PRESENT(vtiITG_GBout))  vtiITG_GBout=vti_GB; 
     IF (PRESENT(vciITG_GBout))  vciITG_GBout=vci_GB; 
     IF (PRESENT(vriITG_GBout))  vriITG_GBout=vri_GB; 

     CALL saturation(2) !set 0 for including all modes, 1 for only ITG, 2 for only TEM, 3 for only ETG
     IF (PRESENT(eefTEM_SIout))   eefTEM_SIout=eef_SI; 
     IF (PRESENT(eefTEM_GBout))   eefTEM_GBout=eef_GB; 
     IF (PRESENT(epfTEM_SIout))   epfTEM_SIout=epf_SI; 

     IF (PRESENT(dfeTEM_SIout))   dfeTEM_SIout=dfe_SI; 
     IF (PRESENT(vteTEM_SIout))  vteTEM_SIout=vte_SI; 
     IF (PRESENT(vceTEM_SIout))  vceTEM_SIout=vce_SI; 
     IF (PRESENT(vreTEM_SIout))  vreTEM_SIout=vre_SI; 

     IF (PRESENT(dfeTEM_GBout))   dfeTEM_GBout=dfe_GB; 
     IF (PRESENT(vteTEM_GBout))  vteTEM_GBout=vte_GB; 
     IF (PRESENT(vceTEM_GBout))  vceTEM_GBout=vce_GB; 
     IF (PRESENT(vreTEM_GBout))  vreTEM_GBout=vre_GB; 

     IF (PRESENT(iefTEM_SIout))   iefTEM_SIout=ief_SI; 
     IF (PRESENT(ipfTEM_SIout))   ipfTEM_SIout=ipf_SI; 
     IF (PRESENT(ivfTEM_SIout))   ivfTEM_SIout=ivf_SI; 

     IF (PRESENT(iefTEM_GBout))   iefTEM_GBout=ief_GB; 
     IF (PRESENT(ivfTEM_GBout))   ivfTEM_GBout=ivf_GB; 

     IF (PRESENT(dfiTEM_SIout))   dfiTEM_SIout=dfi_SI; 
     IF (PRESENT(vtiTEM_SIout))  vtiTEM_SIout=vti_SI; 
     IF (PRESENT(vciTEM_SIout))  vciTEM_SIout=vci_SI; 
     IF (PRESENT(vriTEM_SIout))  vriTEM_SIout=vri_SI; 

     IF (PRESENT(dfiTEM_GBout))  dfiTEM_GBout=dfi_GB; 
     IF (PRESENT(vtiTEM_GBout))  vtiTEM_GBout=vti_GB; 
     IF (PRESENT(vciTEM_GBout))  vciTEM_GBout=vci_GB; 
     IF (PRESENT(vriTEM_GBout))  vriTEM_GBout=vri_GB; 

     CALL saturation(3) !set 0 for including all modes, 1 for only ITG, 2 for only TEM, 3 for only ETG
     IF (PRESENT(eefETG_SIout))   eefETG_SIout=eef_SI; 
     IF (PRESENT(eefETG_GBout))   eefETG_GBout=eef_GB; 

  ENDIF

  CALL saturation(0) !set 0 for including all modes, 1 for only ITG, 2 for only TEM, 3 for only ETG

  IF (myrank==0) THEN
     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time2-time1) / REAL(freq)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: saturation rule calculation time = ',timetot,' s'  

!!!DEBUGGING FOR DIFFERENT FLUID SOLUTIONS
     WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
     OPEN(unit=900, file="output/primitive/rjonsolflu.dat", action="write", status="replace")
     WRITE(900,fmtn) ((REAL(jon_solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(900)
     OPEN(unit=900, file="output/primitive/ijonsolflu.dat", action="write", status="replace")
     WRITE(900,fmtn) ((AIMAG(jon_solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(900)
  ENDIF
  !DEBUG
!!$    OPEN(unit=900, file="FLRec.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((FLRec(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="FLRic.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') (((FLRic(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="FLRep.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((FLRep(i,j),j=1,dimn),i=1,dimx) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="FLRip.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') (((FLRip(i,j,k),j=1,dimn),i=1,dimx),k=1,nions) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="rmodewidth.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((REAL(modewidth(i,j)),j=1,dimn),i=1,dimx) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="rmodeshift.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((REAL(modeshift(i,j)),j=1,dimn),i=1,dimx) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="imodewidth.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((AIMAG(modewidth(i,j)),j=1,dimn),i=1,dimx) ;  CLOSE(900)
!!$
!!$    OPEN(unit=900, file="imodeshift.dat", action="write", status="replace")
!!$    WRITE(900,'(16G15.7)') ((AIMAG(modeshift(i,j)),j=1,dimn),i=1,dimx) ;  CLOSE(900)

  IF (myrank==0) CALL SYSTEM_CLOCK(time3)

  CALL reduceoutput() ! spread all output to all cores
  CALL setoutput() !set all standard output

  !messy setting separated flux output if they exist
  IF (PRESENT(eefITG_SIout))   CALL MPI_AllReduce(eefITG_SIout,eefITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfITG_SIout))   CALL MPI_AllReduce(epfITG_SIout,epfITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeITG_SIout))   CALL MPI_AllReduce(dfeITG_SIout,dfeITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteITG_SIout))  CALL MPI_AllReduce(vteITG_SIout,vteITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceITG_SIout))  CALL MPI_AllReduce(vceITG_SIout,vceITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vreITG_SIout))  CALL MPI_AllReduce(vreITG_SIout,vreITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefITG_SIout))  CALL MPI_AllReduce(iefITG_SIout,iefITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfITG_SIout))  CALL MPI_AllReduce(ivfITG_SIout,ivfITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfITG_SIout))  CALL MPI_AllReduce(ipfITG_SIout,ipfITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiITG_SIout))  CALL MPI_AllReduce(dfiITG_SIout,dfiITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiITG_SIout))  CALL MPI_AllReduce(vtiITG_SIout,vtiITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciITG_SIout))  CALL MPI_AllReduce(vciITG_SIout,vciITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriITG_SIout))  CALL MPI_AllReduce(vriITG_SIout,vriITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefTEM_SIout))  CALL MPI_AllReduce(eefTEM_SIout,eefTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfTEM_SIout))  CALL MPI_AllReduce(epfTEM_SIout,epfTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeTEM_SIout))  CALL MPI_AllReduce(dfeTEM_SIout,dfeTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteTEM_SIout))  CALL MPI_AllReduce(vteTEM_SIout,vteTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceTEM_SIout))  CALL MPI_AllReduce(vceTEM_SIout,vceTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vreTEM_SIout))  CALL MPI_AllReduce(vreTEM_SIout,vreTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefTEM_SIout))  CALL MPI_AllReduce(iefTEM_SIout,iefTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfTEM_SIout))  CALL MPI_AllReduce(ivfTEM_SIout,ivfTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfTEM_SIout))  CALL MPI_AllReduce(ipfTEM_SIout,ipfTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiTEM_SIout))  CALL MPI_AllReduce(dfiTEM_SIout,dfiTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiTEM_SIout))  CALL MPI_AllReduce(vtiTEM_SIout,vtiTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciTEM_SIout))  CALL MPI_AllReduce(vciTEM_SIout,vciTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriTEM_SIout))  CALL MPI_AllReduce(vriTEM_SIout,vriTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefETG_SIout))  CALL MPI_AllReduce(eefETG_SIout,eefETG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefITG_SIout))   eefITG_SIout=eefITG_SIouttmp
  IF (PRESENT(epfITG_SIout))   epfITG_SIout=epfITG_SIouttmp
  IF (PRESENT(dfeITG_SIout))   dfeITG_SIout=dfeITG_SIouttmp
  IF (PRESENT(vteITG_SIout))  vteITG_SIout=vteITG_SIouttmp
  IF (PRESENT(vceITG_SIout))  vceITG_SIout=vceITG_SIouttmp
  IF (PRESENT(vreITG_SIout))  vreITG_SIout=vreITG_SIouttmp

  IF (PRESENT(iefITG_SIout))  iefITG_SIout=iefITG_SIouttmp
  IF (PRESENT(ivfITG_SIout))  ivfITG_SIout=ivfITG_SIouttmp
  IF (PRESENT(ipfITG_SIout))  ipfITG_SIout=ipfITG_SIouttmp
  IF (PRESENT(dfiITG_SIout))  dfiITG_SIout=dfiITG_SIouttmp
  IF (PRESENT(vtiITG_SIout))  vtiITG_SIout=vtiITG_SIouttmp
  IF (PRESENT(vciITG_SIout))  vciITG_SIout=vciITG_SIouttmp
  IF (PRESENT(vriITG_SIout))  vriITG_SIout=vriITG_SIouttmp

  IF (PRESENT(eefTEM_SIout))  eefTEM_SIout=eefTEM_SIouttmp
  IF (PRESENT(epfTEM_SIout))  epfTEM_SIout=epfTEM_SIouttmp
  IF (PRESENT(dfeTEM_SIout))  dfeTEM_SIout=dfeTEM_SIouttmp
  IF (PRESENT(vteTEM_SIout))  vteTEM_SIout=vteTEM_SIouttmp
  IF (PRESENT(vceTEM_SIout))  vceTEM_SIout=vceTEM_SIouttmp
  IF (PRESENT(vreTEM_SIout))  vreTEM_SIout=vreTEM_SIouttmp

  IF (PRESENT(iefTEM_SIout))  iefTEM_SIout=iefTEM_SIouttmp
  IF (PRESENT(ivfTEM_SIout))  ivfTEM_SIout=ivfTEM_SIouttmp
  IF (PRESENT(ipfTEM_SIout))  ipfTEM_SIout=ipfTEM_SIouttmp
  IF (PRESENT(dfiTEM_SIout))  dfiTEM_SIout=dfiTEM_SIouttmp
  IF (PRESENT(vtiTEM_SIout))  vtiTEM_SIout=vtiTEM_SIouttmp
  IF (PRESENT(vciTEM_SIout))  vciTEM_SIout=vciTEM_SIouttmp
  IF (PRESENT(vriTEM_SIout))  vriTEM_SIout=vriTEM_SIouttmp

  IF (PRESENT(eefETG_SIout))  eefETG_SIout=eefETG_SIouttmp
  !!
  IF (PRESENT(eefITG_GBout))   CALL MPI_AllReduce(eefITG_GBout,eefITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeITG_GBout))   CALL MPI_AllReduce(dfeITG_GBout,dfeITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteITG_GBout))  CALL MPI_AllReduce(vteITG_GBout,vteITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceITG_GBout))  CALL MPI_AllReduce(vceITG_GBout,vceITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vreITG_GBout))  CALL MPI_AllReduce(vreITG_GBout,vreITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefITG_GBout))  CALL MPI_AllReduce(iefITG_GBout,iefITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfITG_GBout))  CALL MPI_AllReduce(ivfITG_GBout,ivfITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(dfiITG_GBout))  CALL MPI_AllReduce(dfiITG_GBout,dfiITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiITG_GBout))  CALL MPI_AllReduce(vtiITG_GBout,vtiITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciITG_GBout))  CALL MPI_AllReduce(vciITG_GBout,vciITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriITG_GBout))  CALL MPI_AllReduce(vriITG_GBout,vriITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefTEM_GBout))  CALL MPI_AllReduce(eefTEM_GBout,eefTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeTEM_GBout))  CALL MPI_AllReduce(dfeTEM_GBout,dfeTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteTEM_GBout))  CALL MPI_AllReduce(vteTEM_GBout,vteTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceTEM_GBout))  CALL MPI_AllReduce(vceTEM_GBout,vceTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vreTEM_GBout))  CALL MPI_AllReduce(vreTEM_GBout,vreTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefTEM_GBout))  CALL MPI_AllReduce(iefTEM_GBout,iefTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfTEM_GBout))  CALL MPI_AllReduce(ivfTEM_GBout,ivfTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(dfiTEM_GBout))  CALL MPI_AllReduce(dfiTEM_GBout,dfiTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiTEM_GBout))  CALL MPI_AllReduce(vtiTEM_GBout,vtiTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciTEM_GBout))  CALL MPI_AllReduce(vciTEM_GBout,vciTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriTEM_GBout))  CALL MPI_AllReduce(vriTEM_GBout,vriTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefETG_GBout))  CALL MPI_AllReduce(eefETG_GBout,eefETG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefITG_GBout))   eefITG_GBout=eefITG_GBouttmp
  IF (PRESENT(dfeITG_GBout))   dfeITG_GBout=dfeITG_GBouttmp
  IF (PRESENT(vteITG_GBout))  vteITG_GBout=vteITG_GBouttmp
  IF (PRESENT(vceITG_GBout))  vceITG_GBout=vceITG_GBouttmp
  IF (PRESENT(vreITG_GBout))  vreITG_GBout=vreITG_GBouttmp

  IF (PRESENT(iefITG_GBout))  iefITG_GBout=iefITG_GBouttmp
  IF (PRESENT(ivfITG_GBout))  ivfITG_GBout=ivfITG_GBouttmp

  IF (PRESENT(dfiITG_GBout))  dfiITG_GBout=dfiITG_GBouttmp
  IF (PRESENT(vtiITG_GBout))  vtiITG_GBout=vtiITG_GBouttmp
  IF (PRESENT(vciITG_GBout))  vciITG_GBout=vciITG_GBouttmp
  IF (PRESENT(vriITG_GBout))  vriITG_GBout=vriITG_GBouttmp

  IF (PRESENT(eefTEM_GBout))  eefTEM_GBout=eefTEM_GBouttmp
  IF (PRESENT(dfeTEM_GBout))  dfeTEM_GBout=dfeTEM_GBouttmp
  IF (PRESENT(vteTEM_GBout))  vteTEM_GBout=vteTEM_GBouttmp
  IF (PRESENT(vceTEM_GBout))  vceTEM_GBout=vceTEM_GBouttmp
  IF (PRESENT(vreTEM_GBout))  vreTEM_GBout=vreTEM_GBouttmp

  IF (PRESENT(iefTEM_GBout))  iefTEM_GBout=iefTEM_GBouttmp
  IF (PRESENT(ivfTEM_GBout))  ivfTEM_GBout=ivfTEM_GBouttmp

  IF (PRESENT(dfiTEM_GBout))  dfiTEM_GBout=dfiTEM_GBouttmp
  IF (PRESENT(vtiTEM_GBout))  vtiTEM_GBout=vtiTEM_GBouttmp
  IF (PRESENT(vciTEM_GBout))  vciTEM_GBout=vciTEM_GBouttmp
  IF (PRESENT(vriTEM_GBout))  vriTEM_GBout=vriTEM_GBouttmp

  IF (PRESENT(eefETG_GBout))  eefETG_GBout=eefETG_GBouttmp


  IF (myrank==0) THEN 
     CALL SYSTEM_CLOCK(time4)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time4-time3) / REAL(freq)
     WRITE(stdout,*)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: Second MPI_AllReduce time = ',timetot,' s'  
  ENDIF


  CALL deallocate_endoutput()     

  CALL deallocate_all()

  !Deallocate optional oldsol
  IF (PRESENT(oldsolin)) DEALLOCATE(oldsol)
  IF (PRESENT(oldfdsolin)) DEALLOCATE(oldfdsol)


CONTAINS 

  SUBROUTINE DistriTask(NumTasks,numprocs,rank)
    IMPLICIT NONE
    !     include 'mpif.h'
    INTEGER,DIMENSION(:),ALLOCATABLE :: Request,OK,Finish
    INTEGER,INTENT(IN) :: NumTasks,numprocs,rank
    REAL(kind=DBL) :: tps
    INTEGER :: iradcoord,iwavenum, Task, NoTask, iloop, ierr
    INTEGER :: one=1, minusone=-1
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

       IF(ALLOCATED(Request) .EQV. .FALSE.) ALLOCATE(Request(numprocs-1))
       IF(ALLOCATED(resready) .EQV. .FALSE.)ALLOCATE(resready(numprocs-1))
       IF(ALLOCATED(OK) .EQV. .FALSE.) ALLOCATE(OK(numprocs-1))
       IF(ALLOCATED(Finish) .EQV. .FALSE.) ALLOCATE(Finish(0:numprocs-1))
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
             calltimeinit=MPI_Wtime() ! for timing inside routines
             timeoutflag = .FALSE.
             iwavenum=(NoTask-1)/(dimx) + 1
             iradcoord=MOD((NoTask-1),dimx) + 1
             CALL calc(iradcoord,iwavenum)
             IF ( (timeoutflag .EQV. .TRUE.) .AND. (verbose .EQV. .TRUE.)) WRITE(stdout,'(A,I7,A,I3)') 'Timeout recorded at (p,nu)=',iradcoord,',',iwavenum
             tps=MPI_Wtime()-tps
             tpstot=tpstot+tps  
             IF (verbose .EQV. .TRUE.) THEN
                WRITE(stdout,300) rank,NoTask,tps,tpstot,iradcoord,iwavenum
             ENDIF
300          FORMAT(1x,'rank ',I5,' NoTask ',I7,' time ',F10.3,', total time ',F10.3,' (p,nu)=(',I7,',',I3,')')
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
             calltimeinit=MPI_Wtime() ! for timing inside routines
             timeoutflag = .FALSE.
             iwavenum=(NoTask-1)/(dimx) + 1
             iradcoord=MOD((NoTask-1),dimx) + 1
             CALL calc(iradcoord,iwavenum)
             IF ((timeoutflag .EQV. .TRUE.) .AND. (verbose .EQV. .TRUE.)) WRITE(stdout,'(A,I7,A,I3)') 'Timeout recorded at (p,nu)=',iradcoord,',',iwavenum
             ressend=.TRUE.
             tps=MPI_Wtime()-tps
             tpstot=tpstot+tps
             IF (verbose .EQV. .TRUE.) THEN
                WRITE(stdout,301) rank,NoTask,tps,tpstot,iradcoord,iwavenum
             ENDIF
301          FORMAT(1x,'rank ',I5,' NoTask ',I7,' time ',F10.3,', total time ',F10.3,' (p,nu)=(',I7,',',I3,')')
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE DistriTask

  SUBROUTINE setoutput()

    epf_SIout = epf_SI
    eef_SIout = eef_SI

    evf_SIout = evf_SI
    ipf_SIout = ipf_SI
    ief_SIout = ief_SI
    ivf_SIout = ivf_SI

    IF (PRESENT(solflu_SIout))  solflu_SIout = solflu_SI
    IF (PRESENT(solflu_GBout))  solflu_GBout = solflu_GB
    IF (PRESENT(gam_SIout))     gam_SIout = gam_SI
    IF (PRESENT(gam_GBout))     gam_GBout = gam_GB
    IF (PRESENT(ome_SIout))     ome_SIout = ome_SI
    IF (PRESENT(ome_GBout))     ome_GBout = ome_GB

    IF (PRESENT(epf_GBout))     epf_GBout = epf_GB
    IF (PRESENT(eef_GBout))     eef_GBout = eef_GB
    IF (PRESENT(evf_GBout))     evf_GBout = evf_GB
    IF (PRESENT(epf_cmout))     epf_cmout = epf_cm
    IF (PRESENT(eef_cmout))     eef_cmout = eef_cm
    IF (PRESENT(evf_cmout))     evf_cmout = evf_cm

    IF (PRESENT(ipf_GBout))     ipf_GBout = ipf_GB
    IF (PRESENT(ief_GBout))     ief_GBout = ief_GB
    IF (PRESENT(ivf_GBout))     ivf_GBout = ivf_GB
    IF (PRESENT(ipf_cmout))     ipf_cmout = ipf_cm
    IF (PRESENT(ief_cmout))     ief_cmout = ief_cm
    IF (PRESENT(ivf_cmout))     ivf_cmout = ivf_cm


    IF (PRESENT(modeflagout))   modeflagout = modeflag
    IF (PRESENT(Nustarout))   Nustarout = Nustar
    IF (PRESENT(Zeffxout))   Zeffxout = Zeffx
    IF (PRESENT(phiout))        phiout = phi
    IF (PRESENT(npolout))       npolout = npol
    IF (PRESENT(ecoefsout))     ecoefsout = ecoefs
    IF (PRESENT(cftransout))     cftransout = cftrans
    IF (PRESENT(solfluout))     solfluout = solflu

    IF (PRESENT(krmmuITGout))  krmmuITGout = krmmuITG
    IF (PRESENT(krmmuETGout))  krmmuETGout = krmmuETG
    IF (PRESENT(kperp2out))  kperp2out = kperp2

    IF (PRESENT(modewidthout))  modewidthout = modewidth
    IF (PRESENT(modeshiftout))  modeshiftout = modeshift
    IF (PRESENT(distanout))     distanout = distan
    IF (PRESENT(ntorout))       ntorout = ntor
    IF (PRESENT(solout))        solout = sol
    IF (PRESENT(fdsolout))      fdsolout = fdsol
    IF (PRESENT(Lcirceout))     Lcirceout = Lcirce
    IF (PRESENT(Lpiegeout))     Lpiegeout = Lpiege
    IF (PRESENT(Lecirceout))    Lecirceout = Lecirce
    IF (PRESENT(Lepiegeout))    Lepiegeout= Lepiege
    IF (PRESENT(Lvcirceout))    Lvcirceout = Lvcirce
    IF (PRESENT(Lvpiegeout))    Lvpiegeout= Lvpiege
    IF (PRESENT(Lcirciout))     Lcirciout = Lcirci
    IF (PRESENT(Lpiegiout))     Lpiegiout = Lpiegi
    IF (PRESENT(Lecirciout))    Lecirciout = Lecirci
    IF (PRESENT(Lepiegiout))    Lepiegiout = Lepiegi
    IF (PRESENT(Lvcirciout))    Lvcirciout = Lvcirci
    IF (PRESENT(Lvpiegiout))    Lvpiegiout = Lvpiegi

    IF (phys_meth /= 0) THEN

       IF (PRESENT(Lcircgteout))   Lcircgteout = Lcircgte
       IF (PRESENT(Lpieggteout))   Lpieggteout = Lpieggte
       IF (PRESENT(Lcircgneout))   Lcircgneout = Lcircgne
       IF (PRESENT(Lpieggneout))   Lpieggneout = Lpieggne
       IF (PRESENT(Lcircgueout))   Lcircgueout = Lcircgue
       IF (PRESENT(Lpieggueout))   Lpieggueout = Lpieggue
       IF (PRESENT(Lcircceout))    Lcircceout = Lcircce
       IF (PRESENT(Lpiegceout))    Lpiegceout = Lpiegce

       IF (PRESENT(Lcircgtiout))   Lcircgtiout = Lcircgti 
       IF (PRESENT(Lpieggtiout))   Lpieggtiout = Lpieggti
       IF (PRESENT(Lcircgniout))   Lcircgniout = Lcircgni
       IF (PRESENT(Lpieggniout))   Lpieggniout = Lpieggni
       IF (PRESENT(Lcircguiout))   Lcircguiout = Lcircgui
       IF (PRESENT(Lpiegguiout))   Lpiegguiout = Lpieggui
       IF (PRESENT(Lcircciout))    Lcircciout = Lcircci
       IF (PRESENT(Lpiegciout))    Lpiegciout = Lpiegci

       IF (PRESENT(dfe_SIout))     dfe_SIout = dfe_SI
       IF (PRESENT(vte_SIout))     vte_SIout = vte_SI
       IF (PRESENT(vce_SIout))     vce_SIout = vce_SI
       IF (PRESENT(vre_SIout))     vre_SIout = vre_SI
       IF (PRESENT(dfe_GBout))     dfe_GBout = dfe_GB
       IF (PRESENT(vte_GBout))     vte_GBout = vte_GB
       IF (PRESENT(vce_GBout))     vce_GBout = vce_GB
       IF (PRESENT(vre_GBout))     vre_GBout = vre_GB

       IF (PRESENT(ckeout))        ckeout = cke
       IF (PRESENT(dfi_SIout))     dfi_SIout = dfi_SI
       IF (PRESENT(vti_SIout))     vti_SIout = vti_SI
       IF (PRESENT(vci_SIout))     vci_SIout = vci_SI
       IF (PRESENT(vri_SIout))     vri_SIout = vri_SI
       IF (PRESENT(dfi_GBout))     dfi_GBout = dfi_GB
       IF (PRESENT(vti_GBout))     vti_GBout = vti_GB
       IF (PRESENT(vci_GBout))     vci_GBout = vci_GB
       IF (PRESENT(vri_GBout))     vri_GBout = vri_GB

       IF (PRESENT(ckiout))        ckiout = cki

       IF (phys_meth == 2) THEN
          IF (PRESENT(Lecircgteout))   Lecircgteout = Lecircgte
          IF (PRESENT(Lepieggteout))   Lepieggteout = Lepieggte
          IF (PRESENT(Lecircgneout))   Lecircgneout = Lecircgne
          IF (PRESENT(Lepieggneout))   Lepieggneout = Lepieggne
          IF (PRESENT(Lecircgueout))   Lecircgueout = Lecircgue
          IF (PRESENT(Lepieggueout))   Lepieggueout = Lepieggue
          IF (PRESENT(Lecircceout))    Lecircceout = Lecircce
          IF (PRESENT(Lepiegceout))    Lepiegceout = Lepiegce

          IF (PRESENT(Lecircgtiout))   Lecircgtiout = Lecircgti 
          IF (PRESENT(Lepieggtiout))   Lepieggtiout = Lepieggti
          IF (PRESENT(Lecircgniout))   Lecircgniout = Lecircgni
          IF (PRESENT(Lepieggniout))   Lepieggniout = Lepieggni
          IF (PRESENT(Lecircguiout))   Lecircguiout = Lecircgui
          IF (PRESENT(Lepiegguiout))   Lepiegguiout = Lepieggui
          IF (PRESENT(Lecircciout))    Lecircciout = Lecircci
          IF (PRESENT(Lepiegciout))    Lepiegciout = Lepiegci

          IF (PRESENT(vene_SIout))     vene_SIout = vene_SI
          IF (PRESENT(chiee_SIout))    chiee_SIout = chiee_SI
          IF (PRESENT(vece_SIout))     vece_SIout = vece_SI
          IF (PRESENT(vere_SIout))     vere_SIout = vere_SI
          IF (PRESENT(ckeout))         cekeout = ceke
          IF (PRESENT(veni_SIout))     veni_SIout = veni_SI
          IF (PRESENT(chiei_SIout))    chiei_SIout = chiei_SI
          IF (PRESENT(veci_SIout))     veci_SIout = veci_SI
          IF (PRESENT(veri_SIout))     veri_SIout = veri_SI
          IF (PRESENT(cekiout))        cekiout = ceki
       ENDIF
    ENDIF

  END SUBROUTINE setoutput

END SUBROUTINE qualikiz
