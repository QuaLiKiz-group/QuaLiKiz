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
     & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, ETGmultin, collmultin, & !code specific input
     & epf_SIout,eef_SIout,ipf_SIout,ief_SIout,ivf_SIout, & ! Non optional outputs
     & solflu_SIout, solflu_GBout, gam_SIout,gam_GBout,ome_SIout,ome_GBout, & !growth rate and frequency output
     & epf_GBout,eef_GBout, dfe_SIout,vte_SIout,vce_SIout,epf_cmout,eef_cmout,ckeout, & !electron flux outputs
     & ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout,vti_SIout,vri_SIout,vci_SIout,ipf_cmout,ief_cmout,ivf_cmout,ckiout, & !ion flux outputs
     & dfe_GBout,vte_GBout,vce_GBout,dfi_GBout,vti_GBout,vri_GBout,vci_GBout, &
     & vene_SIout,chiee_SIout,vece_SIout, cekeout, & !heat pinch outputs
     & vene_GBout,chiee_GBout,vece_GBout, &
     & veni_SIout,chiei_SIout,veci_SIout,veri_SIout,cekiout, & 
     & veni_GBout,chiei_GBout,veci_GBout,veri_GBout, & 
     & eefETG_SIout,eefETG_GBout,&  !optional outputs from separation of fluxes
     & modeflagout, Nustarout, Zeffxout, & 
     & phiout, npolout, ecoefsout, cftransout, &  ! poloidal asymmetry outputs for heavy impurities
     & solfluout, modewidthout, modeshiftout, distanout, ntorout, solout, fdsolout,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
     & kperp2out,krmmuITGout,krmmuETGout,&
     & Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lcircgteout, Lpieggteout, Lcircgneout, Lpieggneout, Lcircceout, Lpiegceout, & 
     & Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout, &
     & Lecircgteout, Lepieggteout, Lecircgneout, Lepieggneout, Lecircceout, &
     &Lepiegceout, Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout,&
     oldsolin, oldfdsolin, runcounterin,&
     rhominin,rhomaxin,&
     & eefTEM_SIout,epfTEM_SIout,dfeTEM_SIout,vteTEM_SIout,vceTEM_SIout,& !optional outputs from separation of fluxes
     & eefTEM_GBout,epfTEM_GBout,dfeTEM_GBout,vteTEM_GBout,vceTEM_GBout,&
     & eefITG_SIout,epfITG_SIout,dfeITG_SIout,vteITG_SIout,vceITG_SIout,&
     & eefITG_GBout,epfITG_GBout,dfeITG_GBout,vteITG_GBout,vceITG_GBout,&
     & iefTEM_SIout,ipfTEM_SIout,ivfTEM_SIout,dfiTEM_SIout,vtiTEM_SIout,vciTEM_SIout,vriTEM_SIout,&
     & iefTEM_GBout,ipfTEM_GBout,ivfTEM_GBout,dfiTEM_GBout,vtiTEM_GBout,vciTEM_GBout,vriTEM_GBout,&
     & iefITG_SIout,ipfITG_SIout,ivfITG_SIout,dfiITG_SIout,vtiITG_SIout,vciITG_SIout,vriITG_SIout,&
     & iefITG_GBout,ipfITG_GBout,ivfITG_GBout,dfiITG_GBout,vtiITG_GBout,vciITG_GBout,vriITG_GBout,&
     & veneETG_SIout,chieeETG_SIout,veceETG_SIout,veneETG_GBout,chieeETG_GBout,veceETG_GBout, & !heat pinch outputs
     & veneITG_SIout,chieeITG_SIout,veceITG_SIout,veneITG_GBout,chieeITG_GBout,veceITG_GBout, &
     & veneTEM_SIout,chieeTEM_SIout,veceTEM_SIout,veneTEM_GBout,chieeTEM_GBout,veceTEM_GBout, &
     & veniITG_SIout,chieiITG_SIout,veriITG_SIout,veciITG_SIout,veniITG_GBout,chieiITG_GBout,veriITG_GBout,veciITG_GBout, &
     & veniTEM_SIout,chieiTEM_SIout,veriTEM_SIout,veciTEM_SIout,veniTEM_GBout,chieiTEM_GBout,veriTEM_GBout,veciTEM_GBout, &
     & int_methodin, newt_methodin, newt_convin, int_splitin, reqrelaccin, reqabsaccin)

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
  REAL(kind=DBL) :: cputime1, cputime2, tpstot, timetot, timetot2
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
  REAL(kind=DBL), INTENT(IN) :: relacc1in, relacc2in, timeoutin, ETGmultin, collmultin
  REAL(kind=DBL), OPTIONAL, INTENT(IN) :: rhominin,rhomaxin
  
  !integration parameters
  INTEGER, INTENT(IN), OPTIONAL :: int_methodin, newt_methodin, newt_convin, int_splitin
  REAL(KIND=DBL), INTENT(IN), OPTIONAL :: reqrelaccin, reqabsaccin
  
  INTEGER :: int_methodtmp, newt_methodtmp, newt_convtmp, int_splittmp
  REAL(KIND=DBL) :: reqrelacctmp, reqabsacctmp

  ! List of output variables: 

  ! growth rate and frequency outputs
  COMPLEX(KIND=DBL), DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT) :: solflu_SIout, solflu_GBout, solfluout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: gam_SIout,gam_GBout,ome_SIout,ome_GBout  

  ! final output arrays following saturation rule
  REAL(KIND=DBL), DIMENSION(dimxin), INTENT(OUT)  :: epf_SIout,eef_SIout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), INTENT(OUT)  :: ipf_SIout,ief_SIout,ivf_SIout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefETG_SIout,eefETG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefTEM_SIout,eefTEM_GBout,epfTEM_SIout,epfTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: dfeTEM_SIout,vteTEM_SIout,vceTEM_SIout,dfeTEM_GBout,vteTEM_GBout,vceTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeTEM_SIout,veneTEM_SIout,veceTEM_SIout,chieeTEM_GBout,veneTEM_GBout,veceTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefITG_SIout,eefITG_GBout,epfITG_SIout,epfITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: dfeITG_SIout,vteITG_SIout,vceITG_SIout,dfeITG_GBout,vteITG_GBout,vceITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeETG_SIout,veneETG_SIout,veceETG_SIout,chieeETG_GBout,veneETG_GBout,veceETG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeITG_SIout,veneITG_SIout,veceITG_SIout,chieeITG_GBout,veneITG_GBout,veceITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefTEM_SIout,iefTEM_GBout,ipfTEM_SIout,ipfTEM_GBout,ivfTEM_SIout,ivfTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: dfiTEM_SIout,vtiTEM_SIout,vciTEM_SIout,vriTEM_SIout,dfiTEM_GBout,vtiTEM_GBout,vciTEM_GBout,vriTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: chieiTEM_SIout,veniTEM_SIout,veciTEM_SIout,veriTEM_SIout,chieiTEM_GBout,veniTEM_GBout,veciTEM_GBout,veriTEM_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefITG_SIout,iefITG_GBout,ipfITG_SIout,ipfITG_GBout,ivfITG_SIout,ivfITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: dfiITG_SIout,vtiITG_SIout,vciITG_SIout,vriITG_SIout,dfiITG_GBout,vtiITG_GBout,vciITG_GBout,vriITG_GBout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: chieiITG_SIout,veniITG_SIout,veciITG_SIout,veriITG_SIout,chieiITG_GBout,veniITG_GBout,veciITG_GBout,veriITG_GBout
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: eefETG_SIouttmp,eefETG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: chieeETG_SIouttmp,veneETG_SIouttmp,veceETG_SIouttmp,chieeETG_GBouttmp,veneETG_GBouttmp,veceETG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: eefTEM_SIouttmp,epfTEM_SIouttmp,eefTEM_GBouttmp,epfTEM_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: dfeTEM_SIouttmp,vteTEM_SIouttmp,vceTEM_SIouttmp,dfeTEM_GBouttmp,vteTEM_GBouttmp,vceTEM_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: chieeTEM_SIouttmp,veneTEM_SIouttmp,veceTEM_SIouttmp,chieeTEM_GBouttmp,veneTEM_GBouttmp,veceTEM_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: eefITG_SIouttmp,epfITG_SIouttmp,eefITG_GBouttmp,epfITG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: dfeITG_SIouttmp,vteITG_SIouttmp,vceITG_SIouttmp,dfeITG_GBouttmp,vteITG_GBouttmp,vceITG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: chieeITG_SIouttmp,veneITG_SIouttmp,veceITG_SIouttmp,chieeITG_GBouttmp,veneITG_GBouttmp,veceITG_GBouttmp

  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: iefTEM_SIouttmp,iefTEM_GBouttmp,ipfTEM_SIouttmp,ipfTEM_GBouttmp,ivfTEM_SIouttmp,ivfTEM_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: dfiTEM_SIouttmp,vtiTEM_SIouttmp,vciTEM_SIouttmp,dfiTEM_GBouttmp,vtiTEM_GBouttmp,vciTEM_GBouttmp,vriTEM_SIouttmp,vriTEM_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: chieiTEM_SIouttmp,veniTEM_SIouttmp,veciTEM_SIouttmp,chieiTEM_GBouttmp,veniTEM_GBouttmp,veciTEM_GBouttmp,veriTEM_SIouttmp,veriTEM_GBouttmp

  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: iefITG_SIouttmp,iefITG_GBouttmp,ipfITG_SIouttmp,ipfITG_GBouttmp,ivfITG_SIouttmp,ivfITG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: dfiITG_SIouttmp,vtiITG_SIouttmp,vciITG_SIouttmp,vriITG_SIouttmp,dfiITG_GBouttmp,vtiITG_GBouttmp,vciITG_GBouttmp,vriITG_GBouttmp
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: chieiITG_SIouttmp,veniITG_SIouttmp,veciITG_SIouttmp,chieiITG_GBouttmp,veniITG_GBouttmp,veciITG_GBouttmp,veriITG_SIouttmp,veriITG_GBouttmp

  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: epf_GBout,eef_GBout, dfe_SIout, vte_SIout, vce_SIout, dfe_GBout, vte_GBout, vce_GBout, ckeout, modeflagout, Nustarout, Zeffxout
  REAL(KIND=DBL), DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: vene_SIout, chiee_SIout, vece_SIout, vene_GBout, chiee_GBout, vece_GBout, cekeout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout, vti_SIout, vri_SIout, vci_SIout,dfi_GBout, vti_GBout, vri_GBout, vci_GBout, ckiout
  REAL(KIND=DBL), DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: veni_SIout, veri_SIout, chiei_SIout, veci_SIout, veni_GBout, veri_GBout, chiei_GBout, veci_GBout, cekiout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  ::  epf_cmout, eef_cmout
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
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lcircgteout, Lpieggteout,  Lcircgneout, Lpieggneout, Lcircceout, Lpiegceout
  REAL(KIND=DBL), DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgteout, Lepieggteout,  Lecircgneout, Lepieggneout, Lecircceout, Lepiegceout
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
  
  IF(PRESENT(int_methodin)) THEN
    int_methodtmp = int_methodin
  ELSE
    int_methodtmp = 0
  END IF
  
  IF(PRESENT(newt_methodin)) THEN
    newt_methodtmp = newt_methodin
  ELSE
    newt_methodtmp = 0
  END IF
  
  IF(PRESENT(newt_convin)) THEN
    newt_convtmp = newt_convin
  ELSE
    newt_convtmp = 0
  END IF
  
  IF(PRESENT(int_splitin)) THEN
    int_splittmp = int_splitin
  ELSE
    int_splittmp = 0
  END IF
  
  IF(PRESENT(reqrelaccin)) THEN
    reqrelacctmp = reqrelaccin
  ELSE
    reqrelacctmp = 0.08_DBL
  END IF
  
  IF(PRESENT(reqabsaccin)) THEN
    reqabsacctmp = reqabsaccin
  ELSE
    reqabsacctmp = 0.02_DBL
  END IF
  

  ! Make the input (including derived quantities)
  CALL make_input(dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, kthetarhosin, & !general param
       & xin, rhoin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & 
       & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & 
       & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & 
       & Machtorin, Autorin, Machparin, Auparin, gammaEin, &
       & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin,ETGmultin,collmultin, &
       & int_methodtmp, newt_methodtmp, newt_convtmp, int_splittmp, reqrelacctmp, reqabsacctmp)  !code specific inputs

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
     rhomin=0. 
  ENDIF
  IF (PRESENT(rhomaxin)) THEN
     rhomax=rhomaxin
  ELSE
     rhomax=1.
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
     IF (nproc > 1) THEN
        timetot = REAL(time4-time3) / REAL(freq)
        WRITE(stdout,*)
        WRITE(stdout,"(A,F11.3,A)") 'Profiling: First MPI_AllReduce time = ',timetot,' s'  !final write
     ENDIF
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
     IF (PRESENT(epfITG_GBout))   epfITG_GBout=epf_GB; 

     IF (PRESENT(dfeITG_SIout))   dfeITG_SIout=dfe_SI; 
     IF (PRESENT(vteITG_SIout))  vteITG_SIout=vte_SI; 
     IF (PRESENT(vceITG_SIout))  vceITG_SIout=vce_SI; 

     IF (PRESENT(dfeITG_GBout))   dfeITG_GBout=dfe_GB; 
     IF (PRESENT(vteITG_GBout))  vteITG_GBout=vte_GB; 
     IF (PRESENT(vceITG_GBout))  vceITG_GBout=vce_GB; 

     IF (PRESENT(chieeITG_SIout))   chieeITG_SIout=chiee_SI; 
     IF (PRESENT(veneITG_SIout))  veneITG_SIout=vene_SI; 
     IF (PRESENT(veceITG_SIout))  veceITG_SIout=vece_SI; 

     IF (PRESENT(chieeITG_GBout))   chieeITG_GBout=chiee_GB; 
     IF (PRESENT(veneITG_GBout))  veneITG_GBout=vene_GB; 
     IF (PRESENT(veceITG_GBout))  veceITG_GBout=vece_GB; 

     IF (PRESENT(iefITG_SIout))   iefITG_SIout=ief_SI; 
     IF (PRESENT(ipfITG_SIout))   ipfITG_SIout=ipf_SI; 
     IF (PRESENT(ivfITG_SIout))   ivfITG_SIout=ivf_SI; 

     IF (PRESENT(iefITG_GBout))   iefITG_GBout=ief_GB; 
     IF (PRESENT(ipfITG_GBout))   ipfITG_GBout=ipf_GB; 
     IF (PRESENT(ivfITG_GBout))   ivfITG_GBout=ivf_GB; 

     IF (PRESENT(dfiITG_SIout))   dfiITG_SIout=dfi_SI; 
     IF (PRESENT(vtiITG_SIout))  vtiITG_SIout=vti_SI; 
     IF (PRESENT(vciITG_SIout))  vciITG_SIout=vci_SI; 
     IF (PRESENT(vriITG_SIout))  vriITG_SIout=vri_SI; 

     IF (PRESENT(dfiITG_GBout))   dfiITG_GBout=dfi_GB; 
     IF (PRESENT(vtiITG_GBout))  vtiITG_GBout=vti_GB; 
     IF (PRESENT(vciITG_GBout))  vciITG_GBout=vci_GB; 
     IF (PRESENT(vriITG_GBout))  vriITG_GBout=vri_GB; 

     IF (PRESENT(chieiITG_SIout))   chieiITG_SIout=chiei_SI; 
     IF (PRESENT(veniITG_SIout))  veniITG_SIout=veni_SI; 
     IF (PRESENT(veciITG_SIout))  veciITG_SIout=veci_SI; 
     IF (PRESENT(veriITG_SIout))  veriITG_SIout=veri_SI; 

     IF (PRESENT(chieiITG_GBout))   chieiITG_GBout=chiei_GB; 
     IF (PRESENT(veniITG_GBout))  veniITG_GBout=veni_GB; 
     IF (PRESENT(veciITG_GBout))  veciITG_GBout=veci_GB; 
     IF (PRESENT(veriITG_GBout))  veriITG_GBout=veri_GB; 

     CALL saturation(2) !set 0 for including all modes, 1 for only ITG, 2 for only TEM, 3 for only ETG
     IF (PRESENT(eefTEM_SIout))   eefTEM_SIout=eef_SI; 
     IF (PRESENT(eefTEM_GBout))   eefTEM_GBout=eef_GB; 
     IF (PRESENT(epfTEM_SIout))   epfTEM_SIout=epf_SI; 
     IF (PRESENT(epfTEM_GBout))   epfTEM_GBout=epf_GB; 

     IF (PRESENT(dfeTEM_SIout))   dfeTEM_SIout=dfe_SI; 
     IF (PRESENT(vteTEM_SIout))  vteTEM_SIout=vte_SI; 
     IF (PRESENT(vceTEM_SIout))  vceTEM_SIout=vce_SI; 

     IF (PRESENT(dfeTEM_GBout))   dfeTEM_GBout=dfe_GB; 
     IF (PRESENT(vteTEM_GBout))  vteTEM_GBout=vte_GB; 
     IF (PRESENT(vceTEM_GBout))  vceTEM_GBout=vce_GB; 

     IF (PRESENT(chieeTEM_SIout))   chieeTEM_SIout=chiee_SI; 
     IF (PRESENT(veneTEM_SIout))  veneTEM_SIout=vene_SI; 
     IF (PRESENT(veceTEM_SIout))  veceTEM_SIout=vece_SI; 

     IF (PRESENT(chieeTEM_GBout))   chieeTEM_GBout=chiee_GB; 
     IF (PRESENT(veneTEM_GBout))  veneTEM_GBout=vene_GB; 
     IF (PRESENT(veceTEM_GBout))  veceTEM_GBout=vece_GB; 

     IF (PRESENT(iefTEM_SIout))   iefTEM_SIout=ief_SI; 
     IF (PRESENT(ipfTEM_SIout))   ipfTEM_SIout=ipf_SI; 
     IF (PRESENT(ivfTEM_SIout))   ivfTEM_SIout=ivf_SI; 

     IF (PRESENT(iefTEM_GBout))   iefTEM_GBout=ief_GB; 
     IF (PRESENT(ipfTEM_GBout))   ipfTEM_GBout=ipf_GB; 
     IF (PRESENT(ivfTEM_GBout))   ivfTEM_GBout=ivf_GB; 

     IF (PRESENT(dfiTEM_SIout))   dfiTEM_SIout=dfi_SI; 
     IF (PRESENT(vtiTEM_SIout))  vtiTEM_SIout=vti_SI; 
     IF (PRESENT(vciTEM_SIout))  vciTEM_SIout=vci_SI; 
     IF (PRESENT(vriTEM_SIout))  vriTEM_SIout=vri_SI; 

     IF (PRESENT(dfiTEM_GBout))  dfiTEM_GBout=dfi_GB; 
     IF (PRESENT(vtiTEM_GBout))  vtiTEM_GBout=vti_GB; 
     IF (PRESENT(vciTEM_GBout))  vciTEM_GBout=vci_GB; 
     IF (PRESENT(vriTEM_GBout))  vriTEM_GBout=vri_GB; 

     IF (PRESENT(chieiTEM_SIout))   chieiTEM_SIout=chiei_SI;
     IF (PRESENT(veniTEM_SIout))  veniTEM_SIout=veni_SI; 
     IF (PRESENT(veciTEM_SIout))  veciTEM_SIout=veci_SI; 
     IF (PRESENT(veriTEM_SIout))  veriTEM_SIout=veri_SI; 

     IF (PRESENT(chieiTEM_GBout))  chieiTEM_GBout=chiei_GB; 
     IF (PRESENT(veniTEM_GBout))  veniTEM_GBout=veni_GB; 
     IF (PRESENT(veciTEM_GBout))  veciTEM_GBout=veci_GB; 
     IF (PRESENT(veriTEM_GBout))  veriTEM_GBout=veri_GB; 

  ENDIF

  CALL saturation(0) !set 0 for including all modes, 1 for only ITG, 2 for only TEM

  IF (myrank==0) THEN
     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot2 = REAL(time2-time1) / REAL(freq)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: saturation rule calculation time = ',timetot2,' s'  
     WRITE(stdout,"(A,F7.3,A)") 'Hurrah! QuaLiKiz Job completed! Total time = ',timetot+timetot2,' s'  !final write

     OPEN(unit=900, file="lastruntime.qlk", action="write", status="replace")
     WRITE(900,"(A,F7.3,A)") 'Last completed run time = ',timetot+timetot2,' s'  !final write
     CLOSE(900)

!!!DEBUGGING FOR FLUID SOLUTION
!!$     WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'
!!$     OPEN(unit=900, file="output/primitive/rjonsolflu.dat", action="write", status="replace")
!!$     WRITE(900,fmtn) ((REAL(jon_solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(900)
!!$     OPEN(unit=900, file="output/primitive/ijonsolflu.dat", action="write", status="replace")
!!$     WRITE(900,fmtn) ((AIMAG(jon_solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(900)
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

  IF (PRESENT(eefITG_SIout))   ALLOCATE(eefITG_SIouttmp(dimx))
  IF (PRESENT(epfITG_SIout))  ALLOCATE(epfITG_SIouttmp(dimx))
  IF (PRESENT(dfeITG_SIout))  ALLOCATE(dfeITG_SIouttmp(dimx))
  IF (PRESENT(vteITG_SIout))  ALLOCATE(vteITG_SIouttmp(dimx))
  IF (PRESENT(vceITG_SIout))  ALLOCATE(vceITG_SIouttmp(dimx))
  IF (PRESENT(chieeITG_SIout)) ALLOCATE(chieeITG_SIouttmp(dimx))
  IF (PRESENT(veneITG_SIout)) ALLOCATE(veneITG_SIouttmp(dimx))
  IF (PRESENT(veceITG_SIout)) ALLOCATE(veceITG_SIouttmp(dimx))

  IF (PRESENT(iefITG_SIout))  ALLOCATE(iefITG_SIouttmp(dimx,nions))
  IF (PRESENT(ivfITG_SIout))  ALLOCATE(ivfITG_SIouttmp(dimx,nions))
  IF (PRESENT(ipfITG_SIout))  ALLOCATE(ipfITG_SIouttmp(dimx,nions))
  IF (PRESENT(dfiITG_SIout))  ALLOCATE(dfiITG_SIouttmp(dimx,nions))
  IF (PRESENT(vtiITG_SIout))  ALLOCATE(vtiITG_SIouttmp(dimx,nions))
  IF (PRESENT(vciITG_SIout))  ALLOCATE(vciITG_SIouttmp(dimx,nions))
  IF (PRESENT(vriITG_SIout))  ALLOCATE(vriITG_SIouttmp(dimx,nions))
  IF (PRESENT(chieiITG_SIout)) ALLOCATE(chieiITG_SIouttmp(dimx,nions))
  IF (PRESENT(veniITG_SIout)) ALLOCATE(veniITG_SIouttmp(dimx,nions))
  IF (PRESENT(veciITG_SIout)) ALLOCATE(veciITG_SIouttmp(dimx,nions))
  IF (PRESENT(veriITG_SIout)) ALLOCATE(veriITG_SIouttmp(dimx,nions))

  IF (PRESENT(eefTEM_SIout))   ALLOCATE(eefTEM_SIouttmp(dimx))
  IF (PRESENT(epfTEM_SIout))  ALLOCATE(epfTEM_SIouttmp(dimx))
  IF (PRESENT(dfeTEM_SIout))  ALLOCATE(dfeTEM_SIouttmp(dimx))
  IF (PRESENT(vteTEM_SIout))  ALLOCATE(vteTEM_SIouttmp(dimx))
  IF (PRESENT(vceTEM_SIout))  ALLOCATE(vceTEM_SIouttmp(dimx))
  IF (PRESENT(chieeTEM_SIout)) ALLOCATE(chieeTEM_SIouttmp(dimx))
  IF (PRESENT(veneTEM_SIout)) ALLOCATE(veneTEM_SIouttmp(dimx))
  IF (PRESENT(veceTEM_SIout)) ALLOCATE(veceTEM_SIouttmp(dimx))

  IF (PRESENT(iefTEM_SIout))  ALLOCATE(iefTEM_SIouttmp(dimx,nions))
  IF (PRESENT(ivfTEM_SIout))  ALLOCATE(ivfTEM_SIouttmp(dimx,nions))
  IF (PRESENT(ipfTEM_SIout))  ALLOCATE(ipfTEM_SIouttmp(dimx,nions))
  IF (PRESENT(dfiTEM_SIout))  ALLOCATE(dfiTEM_SIouttmp(dimx,nions))
  IF (PRESENT(vtiTEM_SIout))  ALLOCATE(vtiTEM_SIouttmp(dimx,nions))
  IF (PRESENT(vciTEM_SIout))  ALLOCATE(vciTEM_SIouttmp(dimx,nions))
  IF (PRESENT(vriTEM_SIout))  ALLOCATE(vriTEM_SIouttmp(dimx,nions))
  IF (PRESENT(chieiTEM_SIout)) ALLOCATE(chieiTEM_SIouttmp(dimx,nions))
  IF (PRESENT(veniTEM_SIout)) ALLOCATE(veniTEM_SIouttmp(dimx,nions))
  IF (PRESENT(veciTEM_SIout)) ALLOCATE(veciTEM_SIouttmp(dimx,nions))
  IF (PRESENT(veriTEM_SIout)) ALLOCATE(veriTEM_SIouttmp(dimx,nions))

  IF (PRESENT(eefETG_SIout))  ALLOCATE(eefETG_SIouttmp(dimx))
  IF (PRESENT(chieeETG_SIout)) ALLOCATE(chieeETG_SIouttmp(dimx))
  IF (PRESENT(veneETG_SIout)) ALLOCATE(veneETG_SIouttmp(dimx))
  IF (PRESENT(veceETG_SIout)) ALLOCATE(veceETG_SIouttmp(dimx))

  !
  IF (PRESENT(eefITG_GBout))   ALLOCATE(eefITG_GBouttmp(dimx))
  IF (PRESENT(epfITG_GBout))  ALLOCATE(epfITG_GBouttmp(dimx))
  IF (PRESENT(dfeITG_GBout))  ALLOCATE(dfeITG_GBouttmp(dimx))
  IF (PRESENT(vteITG_GBout))  ALLOCATE(vteITG_GBouttmp(dimx))
  IF (PRESENT(vceITG_GBout))  ALLOCATE(vceITG_GBouttmp(dimx))
  IF (PRESENT(chieeITG_GBout)) ALLOCATE(chieeITG_GBouttmp(dimx))
  IF (PRESENT(veneITG_GBout)) ALLOCATE(veneITG_GBouttmp(dimx))
  IF (PRESENT(veceITG_GBout)) ALLOCATE(veceITG_GBouttmp(dimx))

  IF (PRESENT(iefITG_GBout))  ALLOCATE(iefITG_GBouttmp(dimx,nions))
  IF (PRESENT(ivfITG_GBout))  ALLOCATE(ivfITG_GBouttmp(dimx,nions))
  IF (PRESENT(ipfITG_GBout))  ALLOCATE(ipfITG_GBouttmp(dimx,nions))
  IF (PRESENT(dfiITG_GBout))  ALLOCATE(dfiITG_GBouttmp(dimx,nions))
  IF (PRESENT(vtiITG_GBout))  ALLOCATE(vtiITG_GBouttmp(dimx,nions))
  IF (PRESENT(vciITG_GBout))  ALLOCATE(vciITG_GBouttmp(dimx,nions))
  IF (PRESENT(vriITG_GBout))  ALLOCATE(vriITG_GBouttmp(dimx,nions))
  IF (PRESENT(chieiITG_GBout)) ALLOCATE(chieiITG_GBouttmp(dimx,nions))
  IF (PRESENT(veniITG_GBout)) ALLOCATE(veniITG_GBouttmp(dimx,nions))
  IF (PRESENT(veciITG_GBout)) ALLOCATE(veciITG_GBouttmp(dimx,nions))
  IF (PRESENT(veriITG_GBout)) ALLOCATE(veriITG_GBouttmp(dimx,nions))

  IF (PRESENT(eefTEM_GBout))   ALLOCATE(eefTEM_GBouttmp(dimx))
  IF (PRESENT(epfTEM_GBout))  ALLOCATE(epfTEM_GBouttmp(dimx))
  IF (PRESENT(dfeTEM_GBout))  ALLOCATE(dfeTEM_GBouttmp(dimx))
  IF (PRESENT(vteTEM_GBout))  ALLOCATE(vteTEM_GBouttmp(dimx))
  IF (PRESENT(vceTEM_GBout))  ALLOCATE(vceTEM_GBouttmp(dimx))
  IF (PRESENT(chieeTEM_GBout)) ALLOCATE(chieeTEM_GBouttmp(dimx))
  IF (PRESENT(veneTEM_GBout)) ALLOCATE(veneTEM_GBouttmp(dimx))
  IF (PRESENT(veceTEM_GBout)) ALLOCATE(veceTEM_GBouttmp(dimx))

  IF (PRESENT(iefTEM_GBout))  ALLOCATE(iefTEM_GBouttmp(dimx,nions))
  IF (PRESENT(ivfTEM_GBout))  ALLOCATE(ivfTEM_GBouttmp(dimx,nions))
  IF (PRESENT(ipfTEM_GBout))  ALLOCATE(ipfTEM_GBouttmp(dimx,nions))
  IF (PRESENT(dfiTEM_GBout))  ALLOCATE(dfiTEM_GBouttmp(dimx,nions))
  IF (PRESENT(vtiTEM_GBout))  ALLOCATE(vtiTEM_GBouttmp(dimx,nions))
  IF (PRESENT(vciTEM_GBout))  ALLOCATE(vciTEM_GBouttmp(dimx,nions))
  IF (PRESENT(vriTEM_GBout))  ALLOCATE(vriTEM_GBouttmp(dimx,nions))
  IF (PRESENT(chieiTEM_GBout)) ALLOCATE(chieiTEM_GBouttmp(dimx,nions))
  IF (PRESENT(veniTEM_GBout)) ALLOCATE(veniTEM_GBouttmp(dimx,nions))
  IF (PRESENT(veciTEM_GBout)) ALLOCATE(veciTEM_GBouttmp(dimx,nions))
  IF (PRESENT(veriTEM_GBout)) ALLOCATE(veriTEM_GBouttmp(dimx,nions))

  IF (PRESENT(eefETG_GBout))  ALLOCATE(eefETG_GBouttmp(dimx))
  IF (PRESENT(chieeETG_GBout)) ALLOCATE(chieeETG_GBouttmp(dimx))
  IF (PRESENT(veneETG_GBout)) ALLOCATE(veneETG_GBouttmp(dimx))
  IF (PRESENT(veceETG_GBout)) ALLOCATE(veceETG_GBouttmp(dimx))


  IF (PRESENT(eefITG_SIout))   CALL MPI_AllReduce(eefITG_SIout,eefITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfITG_SIout))   CALL MPI_AllReduce(epfITG_SIout,epfITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeITG_SIout))   CALL MPI_AllReduce(dfeITG_SIout,dfeITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteITG_SIout))  CALL MPI_AllReduce(vteITG_SIout,vteITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceITG_SIout))  CALL MPI_AllReduce(vceITG_SIout,vceITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeITG_SIout))   CALL MPI_AllReduce(chieeITG_SIout,chieeITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneITG_SIout))  CALL MPI_AllReduce(veneITG_SIout,veneITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceITG_SIout))  CALL MPI_AllReduce(veceITG_SIout,veceITG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefITG_SIout))  CALL MPI_AllReduce(iefITG_SIout,iefITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfITG_SIout))  CALL MPI_AllReduce(ivfITG_SIout,ivfITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfITG_SIout))  CALL MPI_AllReduce(ipfITG_SIout,ipfITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiITG_SIout))  CALL MPI_AllReduce(dfiITG_SIout,dfiITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiITG_SIout))  CALL MPI_AllReduce(vtiITG_SIout,vtiITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciITG_SIout))  CALL MPI_AllReduce(vciITG_SIout,vciITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriITG_SIout))  CALL MPI_AllReduce(vriITG_SIout,vriITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieiITG_SIout))  CALL MPI_AllReduce(chieiITG_SIout,chieiITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veniITG_SIout))  CALL MPI_AllReduce(veniITG_SIout,veniITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veciITG_SIout))  CALL MPI_AllReduce(veciITG_SIout,veciITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veriITG_SIout))  CALL MPI_AllReduce(veriITG_SIout,veriITG_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefTEM_SIout))  CALL MPI_AllReduce(eefTEM_SIout,eefTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfTEM_SIout))  CALL MPI_AllReduce(epfTEM_SIout,epfTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeTEM_SIout))  CALL MPI_AllReduce(dfeTEM_SIout,dfeTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteTEM_SIout))  CALL MPI_AllReduce(vteTEM_SIout,vteTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceTEM_SIout))  CALL MPI_AllReduce(vceTEM_SIout,vceTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeTEM_SIout))   CALL MPI_AllReduce(chieeTEM_SIout,chieeTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneTEM_SIout))  CALL MPI_AllReduce(veneTEM_SIout,veneTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceTEM_SIout))  CALL MPI_AllReduce(veceTEM_SIout,veceTEM_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefTEM_SIout))  CALL MPI_AllReduce(iefTEM_SIout,iefTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfTEM_SIout))  CALL MPI_AllReduce(ivfTEM_SIout,ivfTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfTEM_SIout))  CALL MPI_AllReduce(ipfTEM_SIout,ipfTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiTEM_SIout))  CALL MPI_AllReduce(dfiTEM_SIout,dfiTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiTEM_SIout))  CALL MPI_AllReduce(vtiTEM_SIout,vtiTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciTEM_SIout))  CALL MPI_AllReduce(vciTEM_SIout,vciTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriTEM_SIout))  CALL MPI_AllReduce(vriTEM_SIout,vriTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieiTEM_SIout))  CALL MPI_AllReduce(chieiTEM_SIout,chieiTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veniTEM_SIout))  CALL MPI_AllReduce(veniTEM_SIout,veniTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veciTEM_SIout))  CALL MPI_AllReduce(veciTEM_SIout,veciTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veriTEM_SIout))  CALL MPI_AllReduce(veriTEM_SIout,veriTEM_SIouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefETG_SIout))  CALL MPI_AllReduce(eefETG_SIout,eefETG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeETG_SIout))   CALL MPI_AllReduce(chieeETG_SIout,chieeETG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneETG_SIout))  CALL MPI_AllReduce(veneETG_SIout,veneETG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceETG_SIout))  CALL MPI_AllReduce(veceETG_SIout,veceETG_SIouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefITG_GBout))   CALL MPI_AllReduce(eefITG_GBout,eefITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfITG_GBout))   CALL MPI_AllReduce(epfITG_GBout,epfITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeITG_GBout))   CALL MPI_AllReduce(dfeITG_GBout,dfeITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteITG_GBout))   CALL MPI_AllReduce(vteITG_GBout,vteITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceITG_GBout))   CALL MPI_AllReduce(vceITG_GBout,vceITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeITG_GBout)) CALL MPI_AllReduce(chieeITG_GBout,chieeITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneITG_GBout))  CALL MPI_AllReduce(veneITG_GBout,veneITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceITG_GBout))  CALL MPI_AllReduce(veceITG_GBout,veceITG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefITG_GBout))   CALL MPI_AllReduce(iefITG_GBout,iefITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfITG_GBout))   CALL MPI_AllReduce(ipfITG_GBout,ipfITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfITG_GBout))   CALL MPI_AllReduce(ivfITG_GBout,ivfITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiITG_GBout))   CALL MPI_AllReduce(dfiITG_GBout,dfiITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiITG_GBout))   CALL MPI_AllReduce(vtiITG_GBout,vtiITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciITG_GBout))   CALL MPI_AllReduce(vciITG_GBout,vciITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriITG_GBout))   CALL MPI_AllReduce(vriITG_GBout,vriITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieiITG_GBout)) CALL MPI_AllReduce(chieiITG_GBout,chieiITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veniITG_GBout))  CALL MPI_AllReduce(veniITG_GBout,veniITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veciITG_GBout))  CALL MPI_AllReduce(veciITG_GBout,veciITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veriITG_GBout))  CALL MPI_AllReduce(veriITG_GBout,veriITG_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefTEM_GBout))   CALL MPI_AllReduce(eefTEM_GBout,eefTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(epfTEM_GBout))   CALL MPI_AllReduce(epfTEM_GBout,epfTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfeTEM_GBout))   CALL MPI_AllReduce(dfeTEM_GBout,dfeTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vteTEM_GBout))   CALL MPI_AllReduce(vteTEM_GBout,vteTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vceTEM_GBout))   CALL MPI_AllReduce(vceTEM_GBout,vceTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeTEM_GBout)) CALL MPI_AllReduce(chieeTEM_GBout,chieeTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneTEM_GBout))  CALL MPI_AllReduce(veneTEM_GBout,veneTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceTEM_GBout))  CALL MPI_AllReduce(veceTEM_GBout,veceTEM_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(iefTEM_GBout))   CALL MPI_AllReduce(iefTEM_GBout,iefTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ipfTEM_GBout))   CALL MPI_AllReduce(ipfTEM_GBout,ipfTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(ivfTEM_GBout))   CALL MPI_AllReduce(ivfTEM_GBout,ivfTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(dfiTEM_GBout))   CALL MPI_AllReduce(dfiTEM_GBout,dfiTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vtiTEM_GBout))   CALL MPI_AllReduce(vtiTEM_GBout,vtiTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vciTEM_GBout))   CALL MPI_AllReduce(vciTEM_GBout,vciTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(vriTEM_GBout))   CALL MPI_AllReduce(vriTEM_GBout,vriTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieiTEM_GBout)) CALL MPI_AllReduce(chieiTEM_GBout,chieiTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veniTEM_GBout))  CALL MPI_AllReduce(veniTEM_GBout,veniTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veciTEM_GBout))  CALL MPI_AllReduce(veciTEM_GBout,veciTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veriTEM_GBout))  CALL MPI_AllReduce(veriTEM_GBout,veriTEM_GBouttmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefETG_GBout))  CALL MPI_AllReduce(eefETG_GBout,eefETG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(chieeETG_GBout)) CALL MPI_AllReduce(chieeETG_GBout,chieeETG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veneETG_GBout))  CALL MPI_AllReduce(veneETG_GBout,veneETG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)
  IF (PRESENT(veceETG_GBout))  CALL MPI_AllReduce(veceETG_GBout,veceETG_GBouttmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierror)

  IF (PRESENT(eefITG_SIout))  eefITG_SIout=eefITG_SIouttmp
  IF (PRESENT(epfITG_SIout))  epfITG_SIout=epfITG_SIouttmp
  IF (PRESENT(dfeITG_SIout))  dfeITG_SIout=dfeITG_SIouttmp
  IF (PRESENT(vteITG_SIout))  vteITG_SIout=vteITG_SIouttmp
  IF (PRESENT(vceITG_SIout))  vceITG_SIout=vceITG_SIouttmp
  IF (PRESENT(chieeITG_SIout)) chieeITG_SIout=chieeITG_SIouttmp
  IF (PRESENT(veneITG_SIout))  veneITG_SIout=veneITG_SIouttmp
  IF (PRESENT(veceITG_SIout))  veceITG_SIout=veceITG_SIouttmp

  IF (PRESENT(iefITG_SIout))  iefITG_SIout=iefITG_SIouttmp
  IF (PRESENT(ipfITG_SIout))  ipfITG_SIout=ipfITG_SIouttmp
  IF (PRESENT(ivfITG_SIout))  ivfITG_SIout=ivfITG_SIouttmp

  IF (PRESENT(dfiITG_SIout))  dfiITG_SIout=dfiITG_SIouttmp
  IF (PRESENT(vtiITG_SIout))  vtiITG_SIout=vtiITG_SIouttmp
  IF (PRESENT(vciITG_SIout))  vciITG_SIout=vciITG_SIouttmp
  IF (PRESENT(vriITG_SIout))  vriITG_SIout=vriITG_SIouttmp
  IF (PRESENT(chieiITG_SIout))   chieiITG_SIout=chieiITG_SIouttmp
  IF (PRESENT(veniITG_SIout))  veniITG_SIout=veniITG_SIouttmp
  IF (PRESENT(veciITG_SIout))  veciITG_SIout=veciITG_SIouttmp
  IF (PRESENT(veriITG_SIout))  veriITG_SIout=veriITG_SIouttmp

  IF (PRESENT(eefTEM_SIout))  eefTEM_SIout=eefTEM_SIouttmp
  IF (PRESENT(epfTEM_SIout))  epfTEM_SIout=epfTEM_SIouttmp
  IF (PRESENT(dfeTEM_SIout))  dfeTEM_SIout=dfeTEM_SIouttmp
  IF (PRESENT(vteTEM_SIout))  vteTEM_SIout=vteTEM_SIouttmp
  IF (PRESENT(vceTEM_SIout))  vceTEM_SIout=vceTEM_SIouttmp
  IF (PRESENT(chieeTEM_SIout))   chieeTEM_SIout=chieeTEM_SIouttmp
  IF (PRESENT(veneTEM_SIout))  veneTEM_SIout=veneTEM_SIouttmp
  IF (PRESENT(veceTEM_SIout))  veceTEM_SIout=veceTEM_SIouttmp

  IF (PRESENT(iefTEM_SIout))  iefTEM_SIout=iefTEM_SIouttmp
  IF (PRESENT(ipfTEM_SIout))  ipfTEM_SIout=ipfTEM_SIouttmp
  IF (PRESENT(ivfTEM_SIout))  ivfTEM_SIout=ivfTEM_SIouttmp

  IF (PRESENT(dfiTEM_SIout))  dfiTEM_SIout=dfiTEM_SIouttmp
  IF (PRESENT(vtiTEM_SIout))  vtiTEM_SIout=vtiTEM_SIouttmp
  IF (PRESENT(vciTEM_SIout))  vciTEM_SIout=vciTEM_SIouttmp
  IF (PRESENT(vriTEM_SIout))  vriTEM_SIout=vriTEM_SIouttmp
  IF (PRESENT(chieiTEM_SIout))   chieiTEM_SIout=chieiTEM_SIouttmp
  IF (PRESENT(veniTEM_SIout))  veniTEM_SIout=veniTEM_SIouttmp
  IF (PRESENT(veciTEM_SIout))  veciTEM_SIout=veciTEM_SIouttmp
  IF (PRESENT(veriTEM_SIout))  veriTEM_SIout=veriTEM_SIouttmp

  IF (PRESENT(eefETG_SIout))   eefETG_SIout=eefETG_SIouttmp
  IF (PRESENT(chieeETG_SIout)) chieeETG_SIout=chieeETG_SIouttmp
  IF (PRESENT(veneETG_SIout))  veneETG_SIout=veneETG_SIouttmp
  IF (PRESENT(veceETG_SIout))  veceETG_SIout=veceETG_SIouttmp

  IF (PRESENT(eefITG_SIout))   DEALLOCATE(eefITG_SIouttmp)
  IF (PRESENT(epfITG_SIout))   DEALLOCATE(epfITG_SIouttmp)
  IF (PRESENT(dfeITG_SIout))   DEALLOCATE(dfeITG_SIouttmp)
  IF (PRESENT(vteITG_SIout))   DEALLOCATE(vteITG_SIouttmp)
  IF (PRESENT(vceITG_SIout))   DEALLOCATE(vceITG_SIouttmp)
  IF (PRESENT(chieeITG_SIout)) DEALLOCATE(chieeITG_SIouttmp)
  IF (PRESENT(veneITG_SIout))  DEALLOCATE(veneITG_SIouttmp)
  IF (PRESENT(veceITG_SIout))  DEALLOCATE(veceITG_SIouttmp)

  IF (PRESENT(iefITG_SIout))   DEALLOCATE(iefITG_SIouttmp)
  IF (PRESENT(ipfITG_SIout))   DEALLOCATE(ipfITG_SIouttmp)
  IF (PRESENT(ivfITG_SIout))   DEALLOCATE(ivfITG_SIouttmp)

  IF (PRESENT(dfiITG_SIout))   DEALLOCATE(dfiITG_SIouttmp)
  IF (PRESENT(vtiITG_SIout))   DEALLOCATE(vtiITG_SIouttmp)
  IF (PRESENT(vciITG_SIout))   DEALLOCATE(vciITG_SIouttmp)
  IF (PRESENT(vriITG_SIout))   DEALLOCATE(vriITG_SIouttmp)
  IF (PRESENT(chieiITG_SIout)) DEALLOCATE(chieiITG_SIouttmp)
  IF (PRESENT(veniITG_SIout))  DEALLOCATE(veniITG_SIouttmp)
  IF (PRESENT(veciITG_SIout))  DEALLOCATE(veciITG_SIouttmp)
  IF (PRESENT(veriITG_SIout))  DEALLOCATE(veriITG_SIouttmp)

  IF (PRESENT(eefTEM_SIout))   DEALLOCATE(eefTEM_SIouttmp)
  IF (PRESENT(epfTEM_SIout))   DEALLOCATE(epfTEM_SIouttmp)
  IF (PRESENT(dfeTEM_SIout))   DEALLOCATE(dfeTEM_SIouttmp)
  IF (PRESENT(vteTEM_SIout))   DEALLOCATE(vteTEM_SIouttmp)
  IF (PRESENT(vceTEM_SIout))   DEALLOCATE(vceTEM_SIouttmp)
  IF (PRESENT(chieeTEM_SIout)) DEALLOCATE(chieeTEM_SIouttmp)
  IF (PRESENT(veneTEM_SIout))  DEALLOCATE(veneTEM_SIouttmp)
  IF (PRESENT(veceTEM_SIout))  DEALLOCATE(veceTEM_SIouttmp)

  IF (PRESENT(iefTEM_SIout))   DEALLOCATE(iefTEM_SIouttmp)
  IF (PRESENT(ipfTEM_SIout))   DEALLOCATE(ipfTEM_SIouttmp)
  IF (PRESENT(ivfTEM_SIout))   DEALLOCATE(ivfTEM_SIouttmp)

  IF (PRESENT(dfiTEM_SIout))   DEALLOCATE(dfiTEM_SIouttmp)
  IF (PRESENT(vtiTEM_SIout))   DEALLOCATE(vtiTEM_SIouttmp)
  IF (PRESENT(vciTEM_SIout))   DEALLOCATE(vciTEM_SIouttmp)
  IF (PRESENT(vriTEM_SIout))   DEALLOCATE(vriTEM_SIouttmp)
  IF (PRESENT(chieiTEM_SIout)) DEALLOCATE(chieiTEM_SIouttmp)
  IF (PRESENT(veniTEM_SIout))  DEALLOCATE(veniTEM_SIouttmp)
  IF (PRESENT(veciTEM_SIout))  DEALLOCATE(veciTEM_SIouttmp)
  IF (PRESENT(veriTEM_SIout))  DEALLOCATE(veriTEM_SIouttmp)

  IF (PRESENT(eefETG_SIout))   DEALLOCATE(eefETG_SIouttmp)
  IF (PRESENT(chieeETG_SIout)) DEALLOCATE(chieeETG_SIouttmp)
  IF (PRESENT(veneETG_SIout))  DEALLOCATE(veneETG_SIouttmp)
  IF (PRESENT(veceETG_SIout))  DEALLOCATE(veceETG_SIouttmp)

  !!

  IF (PRESENT(eefITG_GBout))  eefITG_GBout=eefITG_GBouttmp
  IF (PRESENT(epfITG_GBout))  epfITG_GBout=epfITG_GBouttmp
  IF (PRESENT(dfeITG_GBout))  dfeITG_GBout=dfeITG_GBouttmp
  IF (PRESENT(vteITG_GBout))  vteITG_GBout=vteITG_GBouttmp
  IF (PRESENT(vceITG_GBout))  vceITG_GBout=vceITG_GBouttmp
  IF (PRESENT(chieeITG_GBout)) chieeITG_GBout=chieeITG_GBouttmp
  IF (PRESENT(veneITG_GBout))  veneITG_GBout=veneITG_GBouttmp
  IF (PRESENT(veceITG_GBout))  veceITG_GBout=veceITG_GBouttmp

  IF (PRESENT(iefITG_GBout))  iefITG_GBout=iefITG_GBouttmp
  IF (PRESENT(ipfITG_GBout))  ipfITG_GBout=ipfITG_GBouttmp
  IF (PRESENT(ivfITG_GBout))  ivfITG_GBout=ivfITG_GBouttmp

  IF (PRESENT(dfiITG_GBout))  dfiITG_GBout=dfiITG_GBouttmp
  IF (PRESENT(vtiITG_GBout))  vtiITG_GBout=vtiITG_GBouttmp
  IF (PRESENT(vciITG_GBout))  vciITG_GBout=vciITG_GBouttmp
  IF (PRESENT(vriITG_GBout))  vriITG_GBout=vriITG_GBouttmp
  IF (PRESENT(chieiITG_GBout))   chieiITG_GBout=chieiITG_GBouttmp
  IF (PRESENT(veniITG_GBout))  veniITG_GBout=veniITG_GBouttmp
  IF (PRESENT(veciITG_GBout))  veciITG_GBout=veciITG_GBouttmp
  IF (PRESENT(veriITG_GBout))  veriITG_GBout=veriITG_GBouttmp

  IF (PRESENT(eefTEM_GBout))  eefTEM_GBout=eefTEM_GBouttmp
  IF (PRESENT(epfTEM_GBout))  epfTEM_GBout=epfTEM_GBouttmp
  IF (PRESENT(dfeTEM_GBout))  dfeTEM_GBout=dfeTEM_GBouttmp
  IF (PRESENT(vteTEM_GBout))  vteTEM_GBout=vteTEM_GBouttmp
  IF (PRESENT(vceTEM_GBout))  vceTEM_GBout=vceTEM_GBouttmp
  IF (PRESENT(chieeTEM_GBout))   chieeTEM_GBout=chieeTEM_GBouttmp
  IF (PRESENT(veneTEM_GBout))  veneTEM_GBout=veneTEM_GBouttmp
  IF (PRESENT(veceTEM_GBout))  veceTEM_GBout=veceTEM_GBouttmp

  IF (PRESENT(iefTEM_GBout))  iefTEM_GBout=iefTEM_GBouttmp
  IF (PRESENT(ipfTEM_GBout))  ipfTEM_GBout=ipfTEM_GBouttmp
  IF (PRESENT(ivfTEM_GBout))  ivfTEM_GBout=ivfTEM_GBouttmp

  IF (PRESENT(dfiTEM_GBout))  dfiTEM_GBout=dfiTEM_GBouttmp
  IF (PRESENT(vtiTEM_GBout))  vtiTEM_GBout=vtiTEM_GBouttmp
  IF (PRESENT(vciTEM_GBout))  vciTEM_GBout=vciTEM_GBouttmp
  IF (PRESENT(vriTEM_GBout))  vriTEM_GBout=vriTEM_GBouttmp
  IF (PRESENT(chieiTEM_GBout))   chieiTEM_GBout=chieiTEM_GBouttmp
  IF (PRESENT(veniTEM_GBout))  veniTEM_GBout=veniTEM_GBouttmp
  IF (PRESENT(veciTEM_GBout))  veciTEM_GBout=veciTEM_GBouttmp
  IF (PRESENT(veriTEM_GBout))  veriTEM_GBout=veriTEM_GBouttmp

  IF (PRESENT(eefETG_GBout))   eefETG_GBout=eefETG_GBouttmp
  IF (PRESENT(chieeETG_GBout)) chieeETG_GBout=chieeETG_GBouttmp
  IF (PRESENT(veneETG_GBout))  veneETG_GBout=veneETG_GBouttmp
  IF (PRESENT(veceETG_GBout))  veceETG_GBout=veceETG_GBouttmp

  IF (PRESENT(eefITG_GBout))   DEALLOCATE(eefITG_GBouttmp)
  IF (PRESENT(epfITG_GBout))   DEALLOCATE(epfITG_GBouttmp)
  IF (PRESENT(dfeITG_GBout))   DEALLOCATE(dfeITG_GBouttmp)
  IF (PRESENT(vteITG_GBout))   DEALLOCATE(vteITG_GBouttmp)
  IF (PRESENT(vceITG_GBout))   DEALLOCATE(vceITG_GBouttmp)
  IF (PRESENT(chieeITG_GBout)) DEALLOCATE(chieeITG_GBouttmp)
  IF (PRESENT(veneITG_GBout))  DEALLOCATE(veneITG_GBouttmp)
  IF (PRESENT(veceITG_GBout))  DEALLOCATE(veceITG_GBouttmp)

  IF (PRESENT(iefITG_GBout))   DEALLOCATE(iefITG_GBouttmp)
  IF (PRESENT(ipfITG_GBout))   DEALLOCATE(ipfITG_GBouttmp)
  IF (PRESENT(ivfITG_GBout))   DEALLOCATE(ivfITG_GBouttmp)

  IF (PRESENT(dfiITG_GBout))   DEALLOCATE(dfiITG_GBouttmp)
  IF (PRESENT(vtiITG_GBout))   DEALLOCATE(vtiITG_GBouttmp)
  IF (PRESENT(vciITG_GBout))   DEALLOCATE(vciITG_GBouttmp)
  IF (PRESENT(vriITG_GBout))   DEALLOCATE(vriITG_GBouttmp)
  IF (PRESENT(chieiITG_GBout)) DEALLOCATE(chieiITG_GBouttmp)
  IF (PRESENT(veniITG_GBout))  DEALLOCATE(veniITG_GBouttmp)
  IF (PRESENT(veciITG_GBout))  DEALLOCATE(veciITG_GBouttmp)
  IF (PRESENT(veriITG_GBout))  DEALLOCATE(veriITG_GBouttmp)

  IF (PRESENT(eefTEM_GBout))   DEALLOCATE(eefTEM_GBouttmp)
  IF (PRESENT(epfTEM_GBout))   DEALLOCATE(epfTEM_GBouttmp)
  IF (PRESENT(dfeTEM_GBout))   DEALLOCATE(dfeTEM_GBouttmp)
  IF (PRESENT(vteTEM_GBout))   DEALLOCATE(vteTEM_GBouttmp)
  IF (PRESENT(vceTEM_GBout))   DEALLOCATE(vceTEM_GBouttmp)
  IF (PRESENT(chieeTEM_GBout)) DEALLOCATE(chieeTEM_GBouttmp)
  IF (PRESENT(veneTEM_GBout))  DEALLOCATE(veneTEM_GBouttmp)
  IF (PRESENT(veceTEM_GBout))  DEALLOCATE(veceTEM_GBouttmp)

  IF (PRESENT(iefTEM_GBout))   DEALLOCATE(iefTEM_GBouttmp)
  IF (PRESENT(ipfTEM_GBout))   DEALLOCATE(ipfTEM_GBouttmp)
  IF (PRESENT(ivfTEM_GBout))   DEALLOCATE(ivfTEM_GBouttmp)

  IF (PRESENT(dfiTEM_GBout))   DEALLOCATE(dfiTEM_GBouttmp)
  IF (PRESENT(vtiTEM_GBout))   DEALLOCATE(vtiTEM_GBouttmp)
  IF (PRESENT(vciTEM_GBout))   DEALLOCATE(vciTEM_GBouttmp)
  IF (PRESENT(vriTEM_GBout))   DEALLOCATE(vriTEM_GBouttmp)
  IF (PRESENT(chieiTEM_GBout)) DEALLOCATE(chieiTEM_GBouttmp)
  IF (PRESENT(veniTEM_GBout))  DEALLOCATE(veniTEM_GBouttmp)
  IF (PRESENT(veciTEM_GBout))  DEALLOCATE(veciTEM_GBouttmp)
  IF (PRESENT(veriTEM_GBout))  DEALLOCATE(veriTEM_GBouttmp)

  IF (PRESENT(eefETG_GBout))   DEALLOCATE(eefETG_GBouttmp)
  IF (PRESENT(chieeETG_GBout)) DEALLOCATE(chieeETG_GBouttmp)
  IF (PRESENT(veneETG_GBout))  DEALLOCATE(veneETG_GBouttmp)
  IF (PRESENT(veceETG_GBout))  DEALLOCATE(veceETG_GBouttmp)

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
    IF (PRESENT(epf_cmout))     epf_cmout = epf_cm
    IF (PRESENT(eef_cmout))     eef_cmout = eef_cm

    IF (PRESENT(eefETG_SIout))   eefETG_SIout=eefETG_SI; 
    IF (PRESENT(eefETG_GBout))   eefETG_GBout=eefETG_GB; 

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
       IF (PRESENT(dfe_GBout))     dfe_GBout = dfe_GB
       IF (PRESENT(vte_GBout))     vte_GBout = vte_GB
       IF (PRESENT(vce_GBout))     vce_GBout = vce_GB
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
          IF (PRESENT(cekeout))         cekeout = ceke
          IF (PRESENT(veni_SIout))     veni_SIout = veni_SI
          IF (PRESENT(chiei_SIout))    chiei_SIout = chiei_SI
          IF (PRESENT(veci_SIout))     veci_SIout = veci_SI
          IF (PRESENT(veri_SIout))     veri_SIout = veri_SI
          IF (PRESENT(cekiout))        cekiout = ceki

          IF (PRESENT(vene_GBout))     vene_GBout = vene_GB
          IF (PRESENT(chiee_GBout))    chiee_GBout = chiee_GB
          IF (PRESENT(vece_GBout))     vece_GBout = vece_GB
          IF (PRESENT(veni_GBout))     veni_GBout = veni_GB
          IF (PRESENT(chiei_GBout))    chiei_GBout = chiei_GB
          IF (PRESENT(veci_GBout))     veci_GBout = veci_GB
          IF (PRESENT(veri_GBout))     veri_GBout = veri_GB

          IF (separateflux .EQV. .TRUE.) THEN
             IF (PRESENT(veneETG_SIout))     veneETG_SIout = veneETG_SI
             IF (PRESENT(chieeETG_SIout))    chieeETG_SIout = chieeETG_SI
             IF (PRESENT(veceETG_SIout))     veceETG_SIout = veceETG_SI

             IF (PRESENT(veneETG_GBout))     veneETG_GBout = veneETG_GB
             IF (PRESENT(chieeETG_GBout))    chieeETG_GBout = chieeETG_GB
             IF (PRESENT(veceETG_GBout))     veceETG_GBout = veceETG_GB
          ENDIF

       ENDIF
    ENDIF

  END SUBROUTINE setoutput

END SUBROUTINE qualikiz
