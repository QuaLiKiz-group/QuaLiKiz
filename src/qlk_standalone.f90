PROGRAM qlk_standalone
  !Standalone driver for qualikiz. Detailed code info in call_qualikiz subroutine header

  !$     use omp_lib
  USE kind
  USE diskio
  USE mpi

  IMPLICIT NONE

  !INCLUDE 'mpif.h'
  ! INTERFACE WITH EXTERNAL QUALIKIZ SUBROUTINE
  INTERFACE 
     SUBROUTINE qualikiz(dimxin, rhoin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, kthetarhosin, & !general param
          & xin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & !geometry input
          & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & !electron input
          & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & !ion input
          & Machtorin, Autorin, Machparin, Auparin, gammaEin, & !rotation input
          & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, ETGmultin, collmultin, & !code specific input
          & epf_SIout,eef_SIout,ipf_SIout,ief_SIout,ivf_SIout, & ! Non optional outputs
          & int_methodin, newt_methodin, QL_methodin, fluid_methodin, newt_convin, int_splitin, normin, reqrelaccin, reqabsaccin, reqrelacc_newtin, reqabsacc_newtin, reqrelacc_QLin, reqabsacc_QLin, & !integration inputs
          & solflu_SIout, solflu_GBout, gam_SIout,gam_GBout,ome_SIout,ome_GBout, & !growth rate and frequency output
          & epf_GBout,eef_GBout, dfe_SIout,vte_SIout,vce_SIout,epf_cmout,eef_cmout,ckeout, & !electron flux outputs
          & ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout,vti_SIout,vri_SIout,vci_SIout,ipf_cmout,ief_cmout,ivf_cmout,ckiout, & !ion flux outputs
          & dfe_GBout,vte_GBout,vce_GBout,dfi_GBout,vti_GBout,vri_GBout,vci_GBout, &
          & vene_SIout,chiee_SIout,vece_SIout, cekeout, & !heat pinch outputs
	  & vene_GBout,chiee_GBout,vece_GBout, &
	  & veni_SIout,chiei_SIout,veci_SIout,veri_SIout,cekiout, & 
	  & veni_GBout,chiei_GBout,veci_GBout,veri_GBout, & 
          & eefETG_SIout,eefETG_GBout,&  !optional ETG outputs
          & modeflagout, Nustarout, Zeffxout, &  
          & phiout, npolout, ecoefsout, cftransout, &  ! poloidal asymmetry outputs for heavy impurities
          & solfluout, modewidthout, modeshiftout, distanout, ntorout, solout, fdsolout,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
          & kperp2out,krmmuITGout,krmmuETGout, &
          & Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lcircgteout, Lpieggteout, Lcircgneout, Lpieggneout, Lcircceout, Lpiegceout, & 
          & Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout, &
          & Lecircgteout, Lepieggteout, Lecircgneout, Lepieggneout, Lecircceout, Lepiegceout, Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout,&
          & oldsolin, oldfdsolin, runcounterin,&
          & rhominin,rhomaxin,&
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
          & veniTEM_SIout,chieiTEM_SIout,veriTEM_SIout,veciTEM_SIout,veniTEM_GBout,chieiTEM_GBout,veriTEM_GBout,veciTEM_GBout)


       INTEGER, INTENT(IN) :: dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, separatefluxin, el_typein
       INTEGER, DIMENSION(dimxin,nionsin), INTENT(IN) :: ion_typein
       REAL, DIMENSION(dimxin,nionsin), INTENT(IN) :: Aiin, Ziin
       REAL, DIMENSION(dimxin), INTENT(IN) :: xin, rhoin, Roin, Rminin, Boin, qxin, smagin, alphaxin
       REAL, DIMENSION(dimnin), INTENT(IN) :: kthetarhosin
       REAL, DIMENSION(dimxin), INTENT(IN) :: Texin, Nexin, Atein, Anein, anisein, danisedrin
       REAL, DIMENSION(dimxin,nionsin), INTENT(IN) :: Tixin, ninormin, Atiin, Aniin, anisin, danisdrin
       REAL, DIMENSION(dimxin), INTENT(IN) :: Machtorin, Autorin, Machparin, Auparin, gammaEin
       INTEGER, INTENT(IN) :: maxrunsin, maxptsin
       REAL, INTENT(IN) :: relacc1in, relacc2in, timeoutin, ETGmultin, collmultin, R0in
       REAL, OPTIONAL, INTENT(IN) :: rhominin,rhomaxin
       
       !Integration testing variables
       INTEGER, INTENT(IN) :: int_methodin, newt_methodin, QL_methodin, fluid_methodin, newt_convin, int_splitin, normin
       REAL, INTENT(IN) :: reqrelaccin, reqabsaccin, reqrelacc_newtin, reqabsacc_newtin, reqrelacc_QLin, reqabsacc_QLin

       ! List of output variables: 
       INTEGER, PARAMETER :: ntheta = 64
       INTEGER, PARAMETER :: numecoefs = 13
       INTEGER, PARAMETER :: numicoefs = 7

       ! growth rate and frequency outputs
       COMPLEX, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT) :: solflu_SIout, solflu_GBout, solfluout
       REAL, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: gam_SIout,gam_GBout,ome_SIout,ome_GBout  

       ! final output arrays following saturation rule
       REAL, DIMENSION(dimxin), INTENT(OUT)  :: epf_SIout,eef_SIout
       REAL, DIMENSION(dimxin,nionsin), INTENT(OUT)  :: ipf_SIout,ief_SIout,ivf_SIout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: epf_GBout,eef_GBout, dfe_SIout, vte_SIout, vce_SIout, dfe_GBout, vte_GBout, vce_GBout, ckeout, modeflagout, Nustarout, Zeffxout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: vene_SIout, chiee_SIout, vece_SIout, cekeout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: vene_GBout, chiee_GBout, vece_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: ipf_GBout,ief_GBout, ivf_GBout,dfi_SIout, vti_SIout, vri_SIout,vci_SIout,dfi_GBout,vti_GBout,vri_GBout,vci_GBout,ckiout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: veni_SIout, chiei_SIout, veci_SIout, veri_SIout,cekiout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: veni_GBout, chiei_GBout, veci_GBout, veri_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefETG_SIout,eefETG_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeETG_SIout,veneETG_SIout,veceETG_SIout,chieeETG_GBout,veneETG_GBout,veceETG_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefTEM_SIout,eefTEM_GBout,epfTEM_SIout,epfTEM_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: dfeTEM_SIout,vteTEM_SIout,vceTEM_SIout,dfeTEM_GBout,vteTEM_GBout,vceTEM_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeTEM_SIout,veneTEM_SIout,veceTEM_SIout,chieeTEM_GBout,veneTEM_GBout,veceTEM_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: eefITG_SIout,eefITG_GBout,epfITG_SIout,epfITG_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: dfeITG_SIout,vteITG_SIout,vceITG_SIout,dfeITG_GBout,vteITG_GBout,vceITG_GBout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT) :: chieeITG_SIout,veneITG_SIout,veceITG_SIout,chieeITG_GBout,veneITG_GBout,veceITG_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefTEM_SIout,iefTEM_GBout,ipfTEM_SIout,ipfTEM_GBout,ivfTEM_SIout,ivfTEM_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: dfiTEM_SIout,vtiTEM_SIout,vciTEM_SIout,vriTEM_SIout,dfiTEM_GBout,vtiTEM_GBout,vciTEM_GBout,vriTEM_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: chieiTEM_SIout,veniTEM_SIout,veciTEM_SIout,veriTEM_SIout,chieiTEM_GBout,veniTEM_GBout,veciTEM_GBout,veriTEM_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: iefITG_SIout,iefITG_GBout,ipfITG_SIout,ipfITG_GBout,ivfITG_SIout,ivfITG_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: dfiITG_SIout,vtiITG_SIout,vciITG_SIout,vriITG_SIout,dfiITG_GBout,vtiITG_GBout,vciITG_GBout,vriITG_GBout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT) :: chieiITG_SIout,veniITG_SIout,veciITG_SIout,veriITG_SIout,chieiITG_GBout,veniITG_GBout,veciITG_GBout,veriITG_GBout      
       REAL, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  ::  epf_cmout, eef_cmout
       REAL, DIMENSION(dimxin,dimnin,nionsin), OPTIONAL, INTENT(OUT) :: ipf_cmout,ief_cmout,ivf_cmout
       REAL, DIMENSION(dimxin,ntheta), OPTIONAL, INTENT(OUT)  ::  phiout
       REAL, DIMENSION(dimxin,ntheta,nionsin), OPTIONAL, INTENT(OUT)  ::  npolout
       REAL, DIMENSION(dimxin,0:nionsin,numecoefs), OPTIONAL, INTENT(OUT)  ::  ecoefsout
       REAL, DIMENSION(dimxin,nionsin,numicoefs), OPTIONAL, INTENT(OUT)  ::  cftransout

       ! optional output arrays from which the saturation rule can be calculated without rerunning dispersion relation solver
       REAL , DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: distanout,ntorout,kperp2out
       REAL , DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: krmmuITGout,krmmuETGout
       COMPLEX, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: modewidthout, modeshiftout
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: solout,fdsolout
       REAL, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lcircgteout, Lpieggteout,  Lcircgneout, Lpieggneout,  Lcircceout, Lpiegceout
       REAL, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgteout, Lepieggteout,  Lecircgneout, Lepieggneout,  Lecircceout, Lepiegceout
       REAL, DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout
       REAL, DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout

       !optional input arrays for going directly to newton solver
       INTEGER, OPTIONAL, INTENT(IN)  :: runcounterin
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(IN)  :: oldsolin, oldfdsolin

     END SUBROUTINE qualikiz
  END INTERFACE

  !Local data dictionary.

  !Time measuring variables
  REAL(kind=DBL) :: tpstot, timetot, datainittime
  INTEGER :: time1, time2, time3,time4,freq
  CHARACTER(len=20) :: myfmt, myint
  CHARACTER(len=:), ALLOCATABLE :: debugdir, outputdir, primitivedir, inputdir

  !MPI variables:
  INTEGER :: ierror, nproc, myrank

  !Input variables into qualikiz subroutine
  INTEGER :: dimx, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, separateflux, el_type, write_primi
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ion_type
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: x, rho, Ro, Rmin, Bo, qx, smag, alphax, kthetarhos
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: Tex, Nex, Ate, Ane, anise, danisedr
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: Tix, ninorm, Ati, Ani, anis, danisdr, Ai, Zi
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: Machtor, Autor, Machpar, Aupar, gammaE
  REAL(KIND=DBL) :: relacc1, relacc2, ETGmult, collmult, timeout, R0
  INTEGER :: maxpts,maxruns
  !Integration testing variables
  INTEGER :: int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm
  REAL(KIND=DBL) :: reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL

  ! Output arrays. The 3 dimensions are 'radial grid', 'kthetarhos grid', 'number of modes'
  REAL(KIND=DBL) , DIMENSION(:), ALLOCATABLE :: krmmuITG,krmmuETG
  REAL(KIND=DBL) , DIMENSION(:,:), ALLOCATABLE :: distan,FLRec,FLRep,ntor,kperp2
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: gamma, Ladia, FLRip, FLRic
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: modewidth, modeshift
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: ommax, solflu
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: sol, fdsol
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: Lcirce, Lpiege, Lecirce, Lepiege, Lcircgte, Lpieggte,  Lcircgne, Lpieggne, Lcircce, Lpiegce
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: Lecircgte, Lepieggte,  Lecircgne, Lepieggne, Lecircce, Lepiegce
  REAL(KIND=DBL), DIMENSION(:,:,:,:), ALLOCATABLE :: Lcirci, Lpiegi, Lecirci, Lepiegi, Lvcirci, Lvpiegi, Lcircgti, Lpieggti, Lcircgni, Lpieggni, Lcircgui, Lpieggui, Lcircci, Lpiegci
  REAL(KIND=DBL), DIMENSION(:,:,:,:), ALLOCATABLE :: Lecircgti, Lepieggti, Lecircgni, Lepieggni, Lecircgui, Lepieggui, Lecircci, Lepiegci

  ! Old solution for non-reset runs
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldrsol,oldisol,oldrfdsol,oldifdsol
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldsol,oldfdsol

  ! Final output arrays following saturation rule. These can be printed as ASCII output
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: solflu_SI, solflu_GB
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: gam_SI,gam_GB,ome_SI,ome_GB
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: epf_SI,epf_GB,eef_SI,eef_GB
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: dfe_SI,vte_SI,vce_SI,dfe_GB,vte_GB,vce_GB,cke,modeflag, Nustar, Zeffx
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: vene_SI,chiee_SI,vece_SI,vene_GB,chiee_GB,vece_GB,ceke
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: ipf_SI,ipf_GB,ief_SI,ief_GB, ivf_SI,ivf_GB
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: dfi_SI,vti_SI,vri_SI,vci_SI,dfi_GB,vti_GB,vri_GB,vci_GB,cki,eef_cm,epf_cm
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: veni_SI,chiei_SI,veci_SI,ceki,veri_SI,veni_GB,chiei_GB,veci_GB,veri_GB
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: ipf_cm,ief_cm,ivf_cm

  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: eefETG_SI,eefETG_GB  !optional outputs from separation of fluxes
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: chieeETG_SI,veneETG_SI,veceETG_SI,chieeETG_GB,veneETG_GB,veceETG_GB
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: eefTEM_SI,epfTEM_SI,dfeTEM_SI,vteTEM_SI,vceTEM_SI,chieeTEM_SI,veneTEM_SI,veceTEM_SI
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: eefTEM_GB,epfTEM_GB,dfeTEM_GB,vteTEM_GB,vceTEM_GB,chieeTEM_GB,veneTEM_GB,veceTEM_GB
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: eefITG_SI,epfITG_SI,dfeITG_SI,vteITG_SI,vceITG_SI,chieeITG_SI,veneITG_SI,veceITG_SI
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: eefITG_GB,epfITG_GB,dfeITG_GB,vteITG_GB,vceITG_GB,chieeITG_GB,veneITG_GB,veceITG_GB

  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: iefTEM_SI,ipfTEM_SI,ivfTEM_SI,dfiTEM_SI,vtiTEM_SI,vciTEM_SI,vriTEM_SI,chieiTEM_SI,veniTEM_SI,veciTEM_SI,veriTEM_SI
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: iefTEM_GB,ipfTEM_GB,ivfTEM_GB,dfiTEM_GB,vtiTEM_GB,vciTEM_GB,vriTEM_GB,chieiTEM_GB,veniTEM_GB,veciTEM_GB,veriTEM_GB
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: iefITG_SI,ipfITG_SI,ivfITG_SI,dfiITG_SI,vtiITG_SI,vciITG_SI,vriITG_SI,chieiITG_SI,veniITG_SI,veciITG_SI,veriITG_SI
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: iefITG_GB,ipfITG_GB,ivfITG_GB,dfiITG_GB,vtiITG_GB,vciITG_GB,vriITG_GB,chieiITG_GB,veniITG_GB,veciITG_GB,veriITG_GB

  ! Poloidal asymmetry variables
  ! Ion density along field line, used was asymmetries are present. 
  ! 1st dimension is radius. 2nd is theta (ntheta points from [0,pi]), 3rd is ion
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: npol, ecoefs, cftrans
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: phi
  INTEGER, PARAMETER :: ntheta = 64
  INTEGER, PARAMETER :: numecoefs = 13
  INTEGER, PARAMETER :: numicoefs = 7

  REAL, PARAMETER :: epsD = 1d-14
  INTEGER :: runcounter ! used for counting runs inside integrated modelling applications for deciding to recalculate all or just jump to newton based on old solutions

  LOGICAL :: exist1, exist2, exist3, exist4, exist5 !used for checking for existence of files

  !DEBUGGING
  CHARACTER(len=20) :: fmtx,fmtn,fmtion,fmtintion,fmtxrow,fmtecoef
  INTEGER :: i,j,k,l,stat,myunit

  CALL mpi_init(ierror)
  CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
  CALL mpi_comm_rank(mpi_comm_world,myrank,ierror)

  myunit = 700+myrank !set a unique file unit per CPU

  ! Begin time measurement
  CALL SYSTEM_CLOCK(time1)

  ! Read input and initialize all arrays
  CALL data_init()

  IF (myrank==0) THEN
     WRITE(stdout,*) ' _________________________________________________________________________________ '
     WRITE(stdout,*) '                                     QUALIKIZ 2.4.0 '
     WRITE(stdout,*) '  gyrokinetic calculation of linear growth rates and quasilinear transport fluxes  '
     WRITE(stdout,*) ' _________________________________________________________________________________ '
     WRITE(stdout,*) ' '
     WRITE(stdout,'(A,I10)') 'Profiling: Amount of parallel processors = ', nproc
     WRITE(stdout,'(A,I10)') 'Profiling: Number of radial or scan points (dimx) =', dimx
     WRITE(stdout,'(A,I10)') 'Profiling: Number of wavenumbers (dimn)  = ', dimn
     WRITE(stdout,*) ' '

     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     datainittime = REAL(time2-time1) / REAL(freq)
  ENDIF

!!!! CALL THE CALCULATION PHASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (separateflux == 0) THEN

     IF (ALLOCATED(oldsol)) THEN !Call with optional old solution input

        IF (phys_meth == 0) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult, & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, &
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, &
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver
                
        ENDIF

        IF (phys_meth == 1) THEN

           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux,  kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &         
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver
        ENDIF

        IF (phys_meth == 2) THEN

           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional ion flux outputs
                & vene_GBout=vene_GB,chiee_GBout=chiee_GB,vece_GBout=vece_GB, &
                & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veci_SIout=veci_SI,veri_SIout=veri_SI,cekiout=ceki, & 
                & veni_GBout=veni_GB,chiei_GBout=chiei_GB,veci_GBout=veci_GB,veri_GBout=veri_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
                & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci, &
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver)
        ENDIF

     ELSE !Don't call with optional old solution input
        IF (phys_meth == 0) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, &
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi)
	ENDIF

        IF (phys_meth == 1) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci)
        ENDIF

        IF (phys_meth == 2) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult, & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional ion flux outputs
                & vene_GBout=vene_GB,chiee_GBout=chiee_GB,vece_GBout=vece_GB, &
                & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veci_SIout=veci_SI,veri_SIout=veri_SI,cekiout=ceki, & 
                & veni_GBout=veni_GB,chiei_GBout=chiei_GB,veci_GBout=veci_GB,veri_GBout=veri_GB, &           
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
                & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci)
        ENDIF
     ENDIF

  ELSE !with separateflux==1

     IF (ALLOCATED(oldsol)) THEN !Call with optional old solution input

        IF (phys_meth == 0) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult, & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, &
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, &
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter,& ! optional inputs for jumping straight to newton solver
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,& !optional outputs from separation of fluxes
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&  
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB)
        ENDIF

        IF (phys_meth == 1) THEN

           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux,  kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &         
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter,& ! optional inputs for jumping straight to newton solver
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,&
                & dfeTEM_SIout=dfeTEM_SI,vteTEM_SIout=vteTEM_SI,vceTEM_SIout=vceTEM_SI,& !optional outputs from separation of fluxes
                & dfeTEM_GBout=dfeTEM_GB,vteTEM_GBout=vteTEM_GB,vceTEM_GBout=vceTEM_GB,&
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&
                & dfeITG_SIout=dfeITG_SI,vteITG_SIout=vteITG_SI,vceITG_SIout=vceITG_SI,&
                & dfeITG_GBout=dfeITG_GB,vteITG_GBout=vteITG_GB,vceITG_GBout=vceITG_GB,&
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & dfiTEM_SIout=dfiTEM_SI,vtiTEM_SIout=vtiTEM_SI,vciTEM_SIout=vciTEM_SI,vriTEM_SIout=vriTEM_SI,&
                & dfiTEM_GBout=dfiTEM_GB,vtiTEM_GBout=vtiTEM_GB,vciTEM_GBout=vciTEM_GB,vriTEM_GBout=vriTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB,&
                & dfiITG_SIout=dfiITG_SI,vtiITG_SIout=vtiITG_SI,vciITG_SIout=vciITG_SI,vriITG_SIout=vriITG_SI,&
                & dfiITG_GBout=dfiITG_GB,vtiITG_GBout=vtiITG_GB,vciITG_GBout=vciITG_GB,vriITG_GBout=vriITG_GB)
        ENDIF

        IF (phys_meth == 2) THEN

           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional ion flux outputs
                & vene_GBout=vene_GB,chiee_GBout=chiee_GB,vece_GBout=vece_GB, &
                & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veci_SIout=veci_SI,veri_SIout=veri_SI,cekiout=ceki, & 
                & veni_GBout=veni_GB,chiei_GBout=chiei_GB,veci_GBout=veci_GB,veri_GBout=veri_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
                & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci, &
                & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter,& ! optional inputs for jumping straight to newton solver
                & chieeETG_SIout=chieeETG_SI,veneETG_SIout=veneETG_SI,veceETG_SIout=veceETG_SI,& !optional outputs from separation of fluxes
                & chieeETG_GBout=chieeETG_GB,veneETG_GBout=veneETG_GB,veceETG_GBout=veceETG_GB,& 
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,&
                & dfeTEM_SIout=dfeTEM_SI,vteTEM_SIout=vteTEM_SI,vceTEM_SIout=vceTEM_SI,& !optional outputs from separation of fluxes
                & dfeTEM_GBout=dfeTEM_GB,vteTEM_GBout=vteTEM_GB,vceTEM_GBout=vceTEM_GB,&
                & chieeTEM_SIout=chieeTEM_SI,veneTEM_SIout=veneTEM_SI,veceTEM_SIout=veceTEM_SI,&
                & chieeTEM_GBout=chieeTEM_GB,veneTEM_GBout=veneTEM_GB,veceTEM_GBout=veceTEM_GB,&
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&
                & dfeITG_SIout=dfeITG_SI,vteITG_SIout=vteITG_SI,vceITG_SIout=vceITG_SI,&
                & dfeITG_GBout=dfeITG_GB,vteITG_GBout=vteITG_GB,vceITG_GBout=vceITG_GB,&
                & chieeITG_SIout=chieeITG_SI,veneITG_SIout=veneITG_SI,veceITG_SIout=veceITG_SI,&
                & chieeITG_GBout=chieeITG_GB,veneITG_GBout=veneITG_GB,veceITG_GBout=veceITG_GB,&
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & dfiTEM_SIout=dfiTEM_SI,vtiTEM_SIout=vtiTEM_SI,vciTEM_SIout=vciTEM_SI,vriTEM_SIout=vriTEM_SI,&
                & dfiTEM_GBout=dfiTEM_GB,vtiTEM_GBout=vtiTEM_GB,vciTEM_GBout=vciTEM_GB,vriTEM_GBout=vriTEM_GB,&
                & chieiTEM_SIout=chieiTEM_SI,veniTEM_SIout=veniTEM_SI,veciTEM_SIout=veciTEM_SI,veriTEM_SIout=veriTEM_SI,&
                & chieiTEM_GBout=chieiTEM_GB,veniTEM_GBout=veniTEM_GB,veciTEM_GBout=veciTEM_GB,veriTEM_GBout=veriTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB,&
                & dfiITG_SIout=dfiITG_SI,vtiITG_SIout=vtiITG_SI,vciITG_SIout=vciITG_SI,vriITG_SIout=vriITG_SI,&
                & dfiITG_GBout=dfiITG_GB,vtiITG_GBout=vtiITG_GB,vciITG_GBout=vciITG_GB,vriITG_GBout=vriITG_GB,&
                & chieiITG_SIout=chieiITG_SI,veniITG_SIout=veniITG_SI,veciITG_SIout=veciITG_SI,veriITG_SIout=veriITG_SI,&
                & chieiITG_GBout=chieiITG_GB,veniITG_GBout=veniITG_GB,veciITG_GBout=veciITG_GB,veriITG_GBout=veriITG_GB)
        ENDIF


     ELSE !Don't call with optional old solution input
        IF (phys_meth == 0) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, &
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi,&
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,& !optional outputs from separation of fluxes
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&  
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB)
	ENDIF

        IF (phys_meth == 1) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,& 
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci,&
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,&
                & dfeTEM_SIout=dfeTEM_SI,vteTEM_SIout=vteTEM_SI,vceTEM_SIout=vceTEM_SI,& !optional outputs from separation of fluxes
                & dfeTEM_GBout=dfeTEM_GB,vteTEM_GBout=vteTEM_GB,vceTEM_GBout=vceTEM_GB,&
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&
                & dfeITG_SIout=dfeITG_SI,vteITG_SIout=vteITG_SI,vceITG_SIout=vceITG_SI,&
                & dfeITG_GBout=dfeITG_GB,vteITG_GBout=vteITG_GB,vceITG_GBout=vceITG_GB,&
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & dfiTEM_SIout=dfiTEM_SI,vtiTEM_SIout=vtiTEM_SI,vciTEM_SIout=vciTEM_SI,vriTEM_SIout=vriTEM_SI,&
                & dfiTEM_GBout=dfiTEM_GB,vtiTEM_GBout=vtiTEM_GB,vciTEM_GBout=vciTEM_GB,vriTEM_GBout=vriTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB,&
                & dfiITG_SIout=dfiITG_SI,vtiITG_SIout=vtiITG_SI,vciITG_SIout=vciITG_SI,vriITG_SIout=vriITG_SI,&
                & dfiITG_GBout=dfiITG_GB,vtiITG_GBout=vtiITG_GB,vciITG_GBout=vciITG_GB,vriITG_GBout=vriITG_GB)
        ENDIF

        IF (phys_meth == 2) THEN
           CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,separateflux, kthetarhos, & !general param
                & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
                & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
                & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
                & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
                & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult, & !code specific inputs
                & epf_SI,eef_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
                & int_method, newt_method, QL_method, fluid_method, newt_conv, int_split, norm, reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt, reqrelacc_QL, reqabsacc_QL,& !integration inputs
                & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
                & epf_GBout=epf_GB,eef_GBout=eef_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm, ckeout=cke, & !optional electron flux outputs
                & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
                & dfe_GBout=dfe_GB,vte_GBout=vte_GB,vce_GBout=vce_GB,dfi_GBout=dfi_GB,vti_GBout=vti_GB,vri_GBout=vri_GB,vci_GBout=vci_GB, &
                & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional ion flux outputs
                & vene_GBout=vene_GB,chiee_GBout=chiee_GB,vece_GBout=vece_GB, &
                & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veci_SIout=veci_SI,veri_SIout=veri_SI,cekiout=ceki, & 
                & veni_GBout=veni_GB,chiei_GBout=chiei_GB,veci_GBout=veci_GB,veri_GBout=veri_GB, &           
                & eefETG_SIout=eefETG_SI,eefETG_GBout=eefETG_GB,&
                & modeflagout=modeflag, Nustarout=Nustar, Zeffxout=Zeffx, & 
                & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
                & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
                & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
                & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lcircgteout=Lcircgte, &
                & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
                & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
                & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
                & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
                & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci,&
                & chieeETG_SIout=chieeETG_SI,veneETG_SIout=veneETG_SI,veceETG_SIout=veceETG_SI,& !optional outputs from separation of fluxes
                & chieeETG_GBout=chieeETG_GB,veneETG_GBout=veneETG_GB,veceETG_GBout=veceETG_GB,& 
                & eefTEM_SIout=eefTEM_SI,eefTEM_GBout=eefTEM_GB,epfTEM_SIout=epfTEM_SI,epfTEM_GBout=epfTEM_GB,&
                & dfeTEM_SIout=dfeTEM_SI,vteTEM_SIout=vteTEM_SI,vceTEM_SIout=vceTEM_SI,& !optional outputs from separation of fluxes
                & dfeTEM_GBout=dfeTEM_GB,vteTEM_GBout=vteTEM_GB,vceTEM_GBout=vceTEM_GB,&
                & chieeTEM_SIout=chieeTEM_SI,veneTEM_SIout=veneTEM_SI,veceTEM_SIout=veceTEM_SI,&
                & chieeTEM_GBout=chieeTEM_GB,veneTEM_GBout=veneTEM_GB,veceTEM_GBout=veceTEM_GB,&
                & eefITG_SIout=eefITG_SI,eefITG_GBout=eefITG_GB,epfITG_SIout=epfITG_SI,epfITG_GBout=epfITG_GB,&
                & dfeITG_SIout=dfeITG_SI,vteITG_SIout=vteITG_SI,vceITG_SIout=vceITG_SI,&
                & dfeITG_GBout=dfeITG_GB,vteITG_GBout=vteITG_GB,vceITG_GBout=vceITG_GB,&
                & chieeITG_SIout=chieeITG_SI,veneITG_SIout=veneITG_SI,veceITG_SIout=veceITG_SI,&
                & chieeITG_GBout=chieeITG_GB,veneITG_GBout=veneITG_GB,veceITG_GBout=veceITG_GB,&
                & iefTEM_SIout=iefTEM_SI,iefTEM_GBout=iefTEM_GB,ipfTEM_SIout=ipfTEM_SI,ipfTEM_GBout=ipfTEM_GB,ivfTEM_SIout=ivfTEM_SI,ivfTEM_GBout=ivfTEM_GB,&
                & dfiTEM_SIout=dfiTEM_SI,vtiTEM_SIout=vtiTEM_SI,vciTEM_SIout=vciTEM_SI,vriTEM_SIout=vriTEM_SI,&
                & dfiTEM_GBout=dfiTEM_GB,vtiTEM_GBout=vtiTEM_GB,vciTEM_GBout=vciTEM_GB,vriTEM_GBout=vriTEM_GB,&
                & chieiTEM_SIout=chieiTEM_SI,veniTEM_SIout=veniTEM_SI,veciTEM_SIout=veciTEM_SI,veriTEM_SIout=veriTEM_SI,&
                & chieiTEM_GBout=chieiTEM_GB,veniTEM_GBout=veniTEM_GB,veciTEM_GBout=veciTEM_GB,veriTEM_GBout=veriTEM_GB,&
                & iefITG_SIout=iefITG_SI,iefITG_GBout=iefITG_GB,ipfITG_SIout=ipfITG_SI,ipfITG_GBout=ipfITG_GB,ivfITG_SIout=ivfITG_SI,ivfITG_GBout=ivfITG_GB,&
                & dfiITG_SIout=dfiITG_SI,vtiITG_SIout=vtiITG_SI,vciITG_SIout=vciITG_SI,vriITG_SIout=vriITG_SI,&
                & dfiITG_GBout=dfiITG_GB,vtiITG_GBout=vtiITG_GB,vciITG_GBout=vciITG_GB,vriITG_GBout=vriITG_GB,&
                & chieiITG_SIout=chieiITG_SI,veniITG_SIout=veniITG_SI,veciITG_SIout=veciITG_SI,veriITG_SIout=veriITG_SI,&
                & chieiITG_GBout=chieiITG_GB,veniITG_GBout=veniITG_GB,veciITG_GBout=veciITG_GB,veriITG_GBout=veriITG_GB)
        ENDIF
     ENDIF
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (myrank == 0) CALL SYSTEM_CLOCK(time3)
  CALL outputascii

  IF (myrank == 0) THEN
     CALL SYSTEM_CLOCK(time4)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time4-time3) / REAL(freq)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: input read and data initialization time = ',datainittime,' s'  
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: output write time = ',timetot,' s'  
  ENDIF

  CALL MPI_Barrier(mpi_comm_world,ierror)

  !WRITE(stdout,*) 'Z function was called ',weidcount,' times' 
  !WRITE(stdout,*)
  !WRITE(stdout,"(A,I0,A,I0)") '*** time: ',timetot, 'seconds, for rank = ',myrank 
  !WRITE(stdout,"(A,I0)") '*** End of job for rank ',myrank
  !IF (myrank==0) WRITE(stdout,"(A,I0,A)") 'We have missed ',Nsolrat,' eigenvalues.'

  !Deallocating all
  CALL deallocate_all()

  IF (myrank==0) THEN 
     CALL SYSTEM_CLOCK(time2)
     CALL SYSTEM_CLOCK(count_rate=freq)
     timetot = REAL(time2-time1) / REAL(freq)
     WRITE(stdout,"(A,F11.3,A)") 'Profiling: Hurrah! Job completed! Total time = ',timetot,' s'  !final write
     WRITE(stdout,*)
     !Write time of run to disk
     OPEN(unit=900, file="lastruntime.dat", action="write", status="replace")
     WRITE(900,"(A,F11.3,A)") 'Last completed run time = ',timetot,' s'  !final write
     CLOSE(900)
  ENDIF

  ! MPI finalization
  CALL mpi_finalize(ierror)

CONTAINS 

  SUBROUTINE data_init()
    !Parallel read data, allocate input and output arrays
    INTEGER, PARAMETER :: ktype = 1 ! BINARY FILES
    INTEGER :: fileno,ierr
    REAL(kind=DBL) :: dummy !dummy variable for obtaining input. Must be real for readvar

    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: dummyn
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: dummyx
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: dummyxnions
    REAL(kind=DBL), DIMENSION(:,:,:), ALLOCATABLE :: dummyxnnumsol

    INTEGER :: dimxtmp,dimntmp,nionstmp,phys_methtmp,coll_flagtmp,rot_flagtmp,verbosetmp, write_primitmp
    INTEGER :: separatefluxtmp,numsolstmp,maxrunstmp,maxptstmp,el_typetmp,runcountertmp
    REAL(kind=DBL) :: relacc1tmp,relacc2tmp,timeouttmp,R0tmp,ETGmulttmp,collmulttmp
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: kthetarhostmp 
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: xtmp,rhotmp,Rotmp,Rmintmp,Botmp,qxtmp,smagtmp,alphaxtmp
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Machtortmp,Autortmp,Machpartmp,Aupartmp,gammaEtmp
    REAL(kind=DBL), DIMENSION(:), ALLOCATABLE :: Textmp,Nextmp,Atetmp,Anetmp,anisetmp,danisedrtmp
    REAL(kind=DBL), DIMENSION(:,:), ALLOCATABLE :: Tixtmp,ninormtmp,Atitmp,Anitmp,anistmp,danisdrtmp, Aitmp, Zitmp
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ion_typetmp
    COMPLEX(kind=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldsoltmp, oldfdsoltmp
    INTEGER :: int_methodtmp, newt_methodtmp, QL_methodtmp, fluid_methodtmp, int_splittmp, normtmp, newt_convtmp
    REAL(KIND=DBL) :: reqrelacctmp, reqabsacctmp, reqrelacc_newttmp, reqabsacc_newttmp, reqrelacc_QLtmp, reqabsacc_QLtmp

    ! READING INPUT ARRAYS FROM BINARY FILES

    inputdir = 'input/'

    fileno = 0
    
    

    ! p{1} Size of radial or scan arrays
    dimx = 0
    IF (myrank == fileno) dimx = INT(readvar(inputdir // 'dimx.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{2} Size of wavenumber arrays
    dimn = 0
    IF (myrank == fileno) dimn = INT(readvar(inputdir // 'dimn.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{3} Number of ions in system
    nions = 0
    IF (myrank == fileno) nions = INT(readvar(inputdir // 'nions.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{9} Number of total saught after solutions
    numsols = 0
    IF (myrank == fileno) numsols = INT(readvar(inputdir // 'numsols.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{12} Number of runs before runcounter resets
    maxruns = 0
    IF (myrank == fileno) maxruns = INT(readvar(inputdir // 'maxruns.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    CALL MPI_Barrier(mpi_comm_world,ierror)
    CALL MPI_AllReduce(dimx,dimxtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(dimn,dimntmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(nions,nionstmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(numsols,numsolstmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(maxruns,maxrunstmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_Barrier(mpi_comm_world,ierror)

    dimx=dimxtmp
    dimn=dimntmp
    nions=nionstmp
    numsols=numsolstmp
    maxruns=maxrunstmp !maxruns called before since it is needed for runcounter evaluation

    ! ALLOCATE TEMP ARRAYS FOR ALLREDUCE
    ALLOCATE(dummyn(dimn))
    ALLOCATE(dummyx(dimx))
    ALLOCATE(dummyxnions(dimx,nions))
    ALLOCATE(dummyxnnumsol(dimx,dimn,numsols))

    ALLOCATE(kthetarhostmp(dimn))
    ALLOCATE(xtmp(dimx))
    ALLOCATE(rhotmp(dimx))
    ALLOCATE(Rotmp(dimx))
    ALLOCATE(Rmintmp(dimx))
    ALLOCATE(Botmp(dimx))
    ALLOCATE(qxtmp(dimx))
    ALLOCATE(smagtmp(dimx))
    ALLOCATE(alphaxtmp(dimx))
    ALLOCATE(Machtortmp(dimx))
    ALLOCATE(Autortmp(dimx))
    ALLOCATE(Machpartmp(dimx))
    ALLOCATE(Aupartmp(dimx))
    ALLOCATE(gammaEtmp(dimx))
    ALLOCATE(Textmp(dimx))
    ALLOCATE(Nextmp(dimx))
    ALLOCATE(Atetmp(dimx))
    ALLOCATE(Anetmp(dimx))
    ALLOCATE(anisetmp(dimx))
    ALLOCATE(danisedrtmp(dimx))

    ALLOCATE(Tixtmp(dimx,nions))
    ALLOCATE(ninormtmp(dimx,nions))
    ALLOCATE(Atitmp(dimx,nions))
    ALLOCATE(Anitmp(dimx,nions))
    ALLOCATE(anistmp(dimx,nions))
    ALLOCATE(danisdrtmp(dimx,nions))
    ALLOCATE(ion_typetmp(dimx,nions))
    ALLOCATE(Aitmp(dimx,nions))
    ALLOCATE(Zitmp(dimx,nions))

    ALLOCATE(oldsoltmp(dimx,dimn,numsols))
    ALLOCATE(oldfdsoltmp(dimx,dimn,numsols))
    
    ! integration parameters
    reqrelacc = 0
    IF (myrank == fileno) reqrelacc = readvar(inputdir // 'reqrelacc.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    reqabsacc = 0
    IF (myrank == fileno) reqabsacc = readvar(inputdir // 'reqabsacc.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    
    reqrelacc_newt = 0
    IF (myrank == fileno) reqrelacc_newt = readvar(inputdir // 'reqrelacc_newt.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    reqabsacc_newt = 0
    IF (myrank == fileno) reqabsacc_newt = readvar(inputdir // 'reqabsacc_newt.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    reqrelacc_QL = 0
    IF (myrank == fileno) reqrelacc_QL = readvar(inputdir // 'reqrelacc_QL.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    reqabsacc_QL = 0
    IF (myrank == fileno) reqabsacc_QL = readvar(inputdir // 'reqabsacc_QL.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    int_method = 0
    IF (myrank == fileno) int_method = INT(readvar(inputdir // 'int_method.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    newt_method = 0
    IF (myrank == fileno) newt_method = INT(readvar(inputdir // 'newt_method.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    QL_method = 0
    IF (myrank == fileno) QL_method = INT(readvar(inputdir // 'QL_method.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    fluid_method = 0
    IF (myrank == fileno) fluid_method = INT(readvar(inputdir // 'fluid_method.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    newt_conv = 0
    IF (myrank == fileno) newt_conv = INT(readvar(inputdir // 'newt_conv.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    int_split = 0
    IF (myrank == fileno) int_split = INT(readvar(inputdir // 'int_split.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    norm = 0
    IF (myrank == fileno) norm = INT(readvar(inputdir // 'norm.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
    
    
    ! p{4} Flag for calculating decomposition of particle and heat transport into diffusive and convective components
    phys_meth = 0
    IF (myrank == fileno) phys_meth = INT(readvar(inputdir // 'phys_meth.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    
    ! p{5} Flag for including collisions
    coll_flag = 0
    IF (myrank == fileno) coll_flag = INT(readvar(inputdir // 'coll_flag.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    write_primi = 0
    IF (myrank == fileno) write_primi = INT(readvar(inputdir // 'write_primi.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{6} Flag for including rotation
    rot_flag = 0
    IF (myrank == fileno) rot_flag = INT(readvar(inputdir // 'rot_flag.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{7} Flag for verbose output
    verbose = 0
    IF (myrank == fileno) verbose = INT(readvar(inputdir // 'verbose.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{8} Flag for separate mode flux output
    separateflux = 0
    IF (myrank == fileno) separateflux = INT(readvar(inputdir // 'separateflux.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{10} 1D integral accuracy
    relacc1 = 0
    IF (myrank == fileno) relacc1 = readvar(inputdir // 'relacc1.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{11} 2D integral accuracy
    relacc2 = 0
    IF (myrank == fileno) relacc2 = readvar(inputdir // 'relacc2.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{13} Maximum number of integrand evaluations in 2D integration routine
    maxpts = 0
    IF (myrank == fileno) maxpts = INT(readvar(inputdir // 'maxpts.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{14} Timeout seconds for a given solution search
    timeout = 0
    IF (myrank == fileno) timeout = readvar(inputdir // 'timeout.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{15} Multiplier for ETG saturation rule (default 1. Mostly for testing)
    ETGmult = 0
    IF (myrank == fileno) ETGmult = readvar(inputdir // 'ETGmult.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{16} Multiplier for collisionality (default 1. Mostly for testing)
    collmult = 0
    IF (myrank == fileno) collmult = readvar(inputdir // 'collmult.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{17} R0 geometric major radius (for normalizations)
    R0 = 0
    IF (myrank == fileno) R0 = readvar(inputdir // 'R0.bin', dummy, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{18} Toroidal wave-number grid
    ALLOCATE(kthetarhos(dimn)); kthetarhos = 0 
    IF (myrank == fileno) kthetarhos = readvar(inputdir // 'kthetarhos.bin', dummyn, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{19} Normalised radial coordinate (midplane radius)
    ALLOCATE(x(dimx)); x=0
    IF (myrank == fileno) x = readvar(inputdir // 'x.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{20} Normalised radial coordinate (midplane radius)
    ALLOCATE(rho(dimx)); rho=0
    IF (myrank == fileno) rho = readvar(inputdir // 'rho.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{21} <Ro> major radius
    ALLOCATE(Ro(dimx)); Ro=0
    IF (myrank == fileno) Ro = readvar(inputdir // 'Ro.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{22} <a> minor radius
    ALLOCATE(Rmin(dimx)); Rmin=0
    IF (myrank == fileno) Rmin = readvar(inputdir // 'Rmin.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{23} B(rho) magnetic field
    ALLOCATE(Bo(dimx)); Bo=0
    IF (myrank == fileno) Bo = readvar(inputdir // 'Bo.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{24} q(rho) profile
    ALLOCATE(qx(dimx)); qx=0
    IF (myrank == fileno) qx = readvar(inputdir // 'q.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{25} s(rho) profile
    ALLOCATE(smag(dimx)); smag=0
    IF (myrank == fileno) smag = readvar(inputdir // 'smag.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{26} alpha(rho) profile
    ALLOCATE(alphax(dimx)); alphax=0
    IF (myrank == fileno) alphax = readvar(inputdir // 'alpha.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{27} Machtor(rho) profile
    ALLOCATE(Machtor(dimx)); Machtor=0
    IF (myrank == fileno) Machtor = readvar(inputdir // 'Machtor.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

!!$    WHERE(ABS(Machtor) < epsD) Machtor = epsD

    ! p{28} Autor(rho) profile
    ALLOCATE(Autor(dimx)); Autor=0
    IF (myrank == fileno) THEN 
       Autor = readvar(inputdir // 'Autor.bin', dummyx, ktype, myunit)
       WHERE(ABS(Autor) < epsD) Autor = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{29} Machpar(rho) profile
    ALLOCATE(Machpar(dimx)); Machpar=0
    IF (myrank == fileno) Machpar = readvar(inputdir // 'Machpar.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 
!!$    WHERE(ABS(Machpar) < epsD) Machpar = epsD

    ! p{30} Aupar(rho) profile
    ALLOCATE(Aupar(dimx)); Aupar=0
    IF (myrank == fileno) THEN
       Aupar = readvar(inputdir // 'Aupar.bin', dummyx, ktype, myunit)
       WHERE(ABS(Aupar) < epsD) Aupar = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{31} gammaE(rho) profile
    ALLOCATE(gammaE(dimx)); gammaE=0
    IF (myrank == fileno) THEN
       gammaE = readvar(inputdir // 'gammaE.bin', dummyx, ktype, myunit)
       WHERE(ABS(gammaE) < epsD) gammaE = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{32} Te(rho) profile
    ALLOCATE(Tex(dimx)); Tex=0
    IF (myrank == fileno) Tex = readvar(inputdir // 'Te.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{33} ne(rho) profile
    ALLOCATE(Nex(dimx)); Nex=0
    IF (myrank == fileno) Nex = readvar(inputdir // 'ne.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{34} R/LTe(rho) profile
    ALLOCATE(Ate(dimx)); Ate=0
    IF (myrank == fileno) THEN 
       Ate = readvar(inputdir // 'Ate.bin', dummyx, ktype, myunit)
       WHERE(ABS(Ate) < epsD) Ate = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{35} R/Lne(rho) profile
    ALLOCATE(Ane(dimx)); Ane=0
    IF (myrank == fileno) THEN 
       Ane = readvar(inputdir // 'Ane.bin', dummyx, ktype, myunit)
       WHERE(ABS(Ane) < epsD) Ane = epsD
       WHERE(Ane+Ate < epsD) Ane = Ane+epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{36} Flag for adiabatic electrons
    el_type = 0
    IF (myrank == fileno) el_type = INT(readvar(inputdir // 'typee.bin', dummy, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{37} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anise(dimx)); anise=0
    IF (myrank == fileno) anise = readvar(inputdir // 'anise.bin', dummyx, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{38} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisedr(dimx)); danisedr=0
    IF (myrank == fileno)  THEN
       danisedr = readvar(inputdir // 'danisdre.bin', dummyx, ktype, myunit)
       WHERE(ABS(danisedr) < epsD) danisedr = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{39} Ti(rho) profiles
    ALLOCATE(Tix(dimx,nions)); Tix=0
    IF (myrank == fileno) Tix = readvar(inputdir // 'Ti.bin', dummyxnions, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{40} ni/ne (rho) profiles
    ALLOCATE(ninorm(dimx,nions)); ninorm=0
    IF (myrank == fileno) ninorm = readvar(inputdir // 'normni.bin', dummyxnions, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{41} R/LTi(rho) profiles
    ALLOCATE(Ati(dimx,nions)); Ati=0
    IF (myrank == fileno) THEN
       Ati = readvar(inputdir // 'Ati.bin', dummyxnions, ktype, myunit)
       WHERE(ABS(Ati) < epsD) Ati = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{42} R/Lni(rho) profiles
    ALLOCATE(Ani(dimx,nions)); Ani=0
    IF (myrank == fileno) THEN 
       Ani = readvar(inputdir // 'Ani.bin', dummyxnions, ktype, myunit)
       WHERE(ABS(Ani) < epsD) Ani = epsD
       WHERE(Ani+Ati < epsD) Ani = Ani+epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{43} Ion types
    ALLOCATE(ion_type(dimx,nions)); ion_type=0
    IF (myrank == fileno) ion_type = INT(readvar(inputdir // 'typei.bin', dummyxnions, ktype, myunit))
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{44} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anis(dimx,1:nions)); anis=0
    IF (myrank == fileno) anis = readvar(inputdir // 'anisi.bin', dummyxnions, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{45} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisdr(dimx,1:nions)); danisdr=0
    IF (myrank == fileno) THEN 
       danisdr = readvar(inputdir // 'danisdri.bin', dummyxnions, ktype, myunit)
       WHERE(ABS(danisdr) < epsD) danisdr = epsD
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{46} Main ion mass
    ALLOCATE(Ai(dimx,nions)); Ai=0
    IF (myrank == fileno) Ai = readvar(inputdir // 'Ai.bin', dummyxnions, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    ! p{47} Main ion charge
    ALLOCATE(Zi(dimx,nions)); Zi=0
    IF (myrank == fileno) Zi = readvar(inputdir // 'Zi.bin', dummyxnions, ktype, myunit)
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    runcounter = 0
    IF (myrank == fileno) THEN
       ! Read and write runcounter input to decide course of action in calcroutines (full solution or start from previous solution)
       INQUIRE(file="runcounter.dat", EXIST=exist1)
       INQUIRE(file="output/primitive/rsol.dat", EXIST=exist2)
       INQUIRE(file="output/primitive/isol.dat", EXIST=exist3)
       INQUIRE(file="output/primitive/rfdsol.dat", EXIST=exist4)
       INQUIRE(file="output/primitive/ifdsol.dat", EXIST=exist5)

       IF ( (exist1) .AND. (exist2) .AND. (exist3) .AND. (exist4) .AND. (exist5) )THEN
          OPEN(myunit, file="runcounter.dat", status="old", action="read")
          READ(myunit,*) runcounter;  CLOSE(myunit)
       ELSE
          runcounter = 0 !First run
       END IF

       IF (runcounter >= maxruns) THEN !Reset if we're at our maximum number of runs
          runcounter = 0
       ENDIF
    ENDIF
    fileno=fileno+1; IF (fileno==nproc) fileno=0 

    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G16.7)'

    CALL MPI_Barrier(mpi_comm_world,ierror)

    !Now do MPIAllReduce to inputs. We need runcounter for potential next stage, so
    !the reduce operations are split
    CALL MPI_AllReduce(phys_meth,phys_methtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(coll_flag,coll_flagtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(write_primi,write_primitmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(rot_flag,rot_flagtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(verbose,verbosetmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(separateflux,separatefluxtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(numsols,numsolstmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(maxpts,maxptstmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(el_type,el_typetmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(runcounter,runcountertmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(int_method,int_methodtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(newt_method,newt_methodtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(QL_method,QL_methodtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(fluid_method,fluid_methodtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(newt_conv, newt_convtmp, 1, MPI_INTEGER, MPI_SUM, mpi_comm_world, ierr)
    CALL MPI_AllReduce(int_split,int_splittmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(norm,normtmp,1,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)

    CALL MPI_AllReduce(reqrelacc,reqrelacctmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(reqabsacc,reqabsacctmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(reqrelacc_newt,reqrelacc_newttmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(reqabsacc_newt,reqabsacc_newttmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(reqrelacc_QL,reqrelacc_QLtmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(reqabsacc_QL,reqabsacc_QLtmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(relacc1,relacc1tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(relacc2,relacc2tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(timeout,timeouttmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ETGmult,ETGmulttmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(collmult,collmulttmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(R0,R0tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(kthetarhos,kthetarhostmp,dimn,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(x,xtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(rho,rhotmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ro,Rotmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Rmin,Rmintmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Bo,Botmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(qx,qxtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(smag,smagtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(alphax,alphaxtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Machtor,Machtortmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Autor,Autortmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Machpar,Machpartmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Aupar,Aupartmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(gammaE,gammaEtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Tex,Textmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Nex,Nextmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ate,Atetmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ane,Anetmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(anise,anisetmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(danisedr,danisedrtmp,dimx,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Tix,Tixtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ninorm,ninormtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ati,Atitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ani,Anitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(anis,anistmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(danisdr,danisdrtmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(ion_type,ion_typetmp,dimx*nions,MPI_INTEGER,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Ai,Aitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    CALL MPI_AllReduce(Zi,Zitmp,dimx*nions,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,ierr)
    

    CALL MPI_Barrier(mpi_comm_world,ierror)

    
    phys_meth=phys_methtmp
    coll_flag=coll_flagtmp
    write_primi=write_primitmp
    rot_flag=rot_flagtmp
    verbose=verbosetmp
    separateflux=separatefluxtmp
    
    int_method = int_methodtmp
    newt_method = newt_methodtmp
    QL_method = QL_methodtmp
    fluid_method = fluid_methodtmp
    int_split = int_splittmp
    norm = normtmp
    newt_conv = newt_convtmp
    reqrelacc = reqrelacctmp
    reqabsacc = reqabsacctmp
    reqrelacc_newt = reqrelacc_newttmp
    reqabsacc_newt = reqabsacc_newttmp
    reqrelacc_QL = reqrelacc_QLtmp
    reqabsacc_QL = reqabsacc_QLtmp

    maxpts=maxptstmp
    runcounter=runcountertmp
    el_type=el_typetmp
    relacc1=relacc1tmp
    relacc2=relacc2tmp
    timeout=timeouttmp
    ETGmult=ETGmulttmp
    collmult=collmulttmp
    R0=R0tmp
    kthetarhos=kthetarhostmp
    x=xtmp
    rho=rhotmp
    Ro=Rotmp
    Rmin=Rmintmp
    Bo=Botmp
    qx=qxtmp
    smag=smagtmp
    alphax=alphaxtmp
    Machtor=Machtortmp
    Autor=Autortmp
    Machpar=Machpartmp
    Aupar=Aupartmp
    gammaE=gammaEtmp
    Tex=Textmp
    Nex=Nextmp
    Ate=Atetmp
    Ane=Anetmp
    anise=anisetmp
    danisedr=danisedrtmp
    Tix=Tixtmp
    ninorm=ninormtmp
    Ati=Atitmp
    Ani=Anitmp
    anis=anistmp
    danisdr=danisdrtmp
    ion_type=ion_typetmp
    Ai=Aitmp
    Zi=Zitmp

    IF (runcounter /= 0) THEN !load old rsol and isol if we're not doing a reset run
       ALLOCATE( oldrsol (dimx, dimn, numsols) ); oldrsol = 0
       ALLOCATE( oldisol (dimx, dimn, numsols) ); oldisol = 0
       ALLOCATE( oldrfdsol (dimx, dimn, numsols) ) ; oldrfdsol = 0
       ALLOCATE( oldifdsol (dimx, dimn, numsols) ); oldifdsol = 0
       ALLOCATE( oldsol (dimx, dimn, numsols) ); oldsol = 0
       ALLOCATE( oldfdsol (dimx, dimn, numsols) ); oldfdsol = 0

       primitivedir = "output/primitive/"
       myfmt = 'G16.7E3'
       IF (myrank == fileno) THEN
           oldrsol = readvar(primitivedir // 'rsol.dat', dummyxnnumsol, ktype, myunit)
           CALL writevar(primitivedir // 'rsol_old.dat', oldrsol, myfmt, myunit)
       ENDIF
       fileno=fileno+1; IF (fileno==nproc) fileno=0 

       IF (myrank == fileno) THEN
           oldisol = readvar(primitivedir // 'isol.dat', dummyxnnumsol, ktype, myunit)
       ENDIF
       fileno=fileno+1; IF (fileno==nproc) fileno=0 

       IF (myrank == fileno) THEN
           oldrfdsol = readvar(primitivedir // 'rfdsol.dat', dummyxnnumsol, ktype, myunit)
       ENDIF
       fileno=fileno+1; IF (fileno==nproc) fileno=0 

       IF (myrank == fileno) THEN
           oldifdsol = readvar(primitivedir // 'ifdsol.dat', dummyxnnumsol, ktype, myunit)
       ENDIF
       fileno=fileno+1; IF (fileno==nproc) fileno=0 

       oldsol = CMPLX(oldrsol,oldisol)
       oldfdsol = CMPLX(oldrfdsol,oldifdsol)

       DEALLOCATE( oldrsol )
       DEALLOCATE( oldisol )
       DEALLOCATE( oldrfdsol )
       DEALLOCATE( oldifdsol )

       CALL MPI_Barrier(mpi_comm_world,ierror)
       CALL MPI_AllReduce(oldsol,oldsoltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_AllReduce(oldfdsol,oldfdsoltmp,dimx*dimn*numsols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_world,ierr)
       CALL MPI_Barrier(mpi_comm_world,ierror)

       oldsol=oldsoltmp;
       oldfdsol=oldfdsoltmp;

    ENDIF

    CALL MPI_Barrier(mpi_comm_world,ierror)

    IF (myrank == 0) THEN
       OPEN(unit=700, file="runcounter.dat", status="replace", action="write") !Replace old runcounter with new runcounter
       WRITE(700,*) runcounter + 1 ; CLOSE(700)
    ENDIF

    !DEBUGGING WRITE OUT ALL INPUT TO ASCII FILE

    myint='I15'
    myfmt='G16.7E3'
    debugdir='debug/'
    CALL writevar(debugdir // 'dimx.dat', dimx, myint, fileno, .FALSE.)
    fileno=fileno+1
    CALL writevar(debugdir // 'dimn.dat', dimn, myint, fileno, .FALSE.)
    fileno=fileno+1
    CALL writevar(debugdir // 'nions.dat', nions, myint, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'phys_meth.dat', phys_meth, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'coll_flag.dat', coll_flag, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'rot_flag.dat', rot_flag, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'verbose.dat', verbose, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'separateflux.dat', separateflux, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'write_primi.dat', write_primi, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'numsols.dat', numsols, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'kthetarhos.dat', kthetarhos, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'x.dat', x, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'rho.dat', rho, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ro.dat', Ro, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'R0.dat', R0, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Rmin.dat', Rmin, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Bo.dat', Bo, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'q.dat', qx, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'smag.dat', smag, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'alpha.dat', alphax, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Machtor.dat', Machtor, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Autor.dat', Autor, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Machpar.dat', Machpar, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Aupar.dat', Aupar, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'gammaE.dat', gammaE, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Te.dat', Tex, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'ne.dat', Nex, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ate.dat', Ate, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ane.dat', Ane, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'typee.dat', el_type, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ai.dat', Ai, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Zi.dat', Zi, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ti.dat', Tix, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'normni.dat', ninorm, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ati.dat', Ati, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'Ani.dat', Ani, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'typei.dat', ion_type, myint, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'maxpts.dat', maxpts, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'maxruns.dat', maxruns, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'relacc1.dat', relacc1, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'relacc2.dat', relacc2, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'timeout.dat', timeout, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'ETGmult.dat', ETGmult, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'collmult.dat', collmult, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(debugdir // 'reqrelacc.dat', reqrelacc, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'reqabsacc.dat', reqabsacc, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'reqrelacc_newt.dat', reqrelacc_newt, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'reqabsacc_newt.dat', reqabsacc_newt, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'reqrelacc_QL.dat', reqrelacc_QL, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'reqabsacc_QL.dat', reqabsacc_QL, myfmt, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'int_method.dat', int_method, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'newt_method.dat', newt_method, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'QL_method.dat', QL_method, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'fluid_method.dat', fluid_method, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'int_split.dat', int_split, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'norm.dat', norm, myint, fileno)
    fileno = fileno+1
    CALL writevar(debugdir // 'newt_conv.dat', newt_conv, myint, fileno)
    fileno = fileno+1

    !Allocate output
    ALLOCATE(solflu_SI(dimx,dimn))
    ALLOCATE(solflu_GB(dimx,dimn))
    ALLOCATE(gam_SI(dimx,dimn,numsols))
    ALLOCATE(gam_GB(dimx,dimn,numsols))
    ALLOCATE(ome_SI(dimx,dimn,numsols))
    ALLOCATE(ome_GB(dimx,dimn,numsols))
    ALLOCATE(epf_SI(dimx))
    ALLOCATE(epf_GB(dimx))
    ALLOCATE(eef_SI(dimx))
    ALLOCATE(eef_GB(dimx))
    ALLOCATE(eefETG_SI(dimx))
    ALLOCATE(eefETG_GB(dimx))

    IF (separateflux == 1) THEN
       ALLOCATE(eefTEM_SI(dimx))
       ALLOCATE(eefITG_SI(dimx))
       ALLOCATE(eefTEM_GB(dimx))
       ALLOCATE(eefITG_GB(dimx))
       ALLOCATE(epfITG_SI(dimx))
       ALLOCATE(epfTEM_SI(dimx))
       ALLOCATE(epfITG_GB(dimx))
       ALLOCATE(epfTEM_GB(dimx))
    ENDIF

    IF (phys_meth /= 0) THEN
       ALLOCATE(dfe_SI(dimx))
       ALLOCATE(vte_SI(dimx))
       ALLOCATE(vce_SI(dimx))
       ALLOCATE(dfe_GB(dimx))
       ALLOCATE(vte_GB(dimx))
       ALLOCATE(vce_GB(dimx))    
       ALLOCATE(cke(dimx))

       IF (separateflux == 1) THEN
          ALLOCATE(dfeITG_SI(dimx))
          ALLOCATE(vteITG_SI(dimx))
          ALLOCATE(vceITG_SI(dimx))
          ALLOCATE(dfeTEM_SI(dimx))
          ALLOCATE(vteTEM_SI(dimx))
          ALLOCATE(vceTEM_SI(dimx))

          ALLOCATE(dfeITG_GB(dimx))
          ALLOCATE(vteITG_GB(dimx))
          ALLOCATE(vceITG_GB(dimx))
          ALLOCATE(dfeTEM_GB(dimx))
          ALLOCATE(vteTEM_GB(dimx))
          ALLOCATE(vceTEM_GB(dimx))

       ENDIF

       IF (phys_meth == 2) THEN
          ALLOCATE(vene_SI(dimx))
          ALLOCATE(chiee_SI(dimx))
          ALLOCATE(vece_SI(dimx))
          ALLOCATE(vene_GB(dimx))
          ALLOCATE(chiee_GB(dimx))
          ALLOCATE(vece_GB(dimx))
          ALLOCATE(ceke(dimx))

          IF (separateflux == 1) THEN
             ALLOCATE(chieeITG_SI(dimx))
             ALLOCATE(veneITG_SI(dimx))
             ALLOCATE(veceITG_SI(dimx))
             ALLOCATE(chieeTEM_SI(dimx))
             ALLOCATE(veneTEM_SI(dimx))
             ALLOCATE(veceTEM_SI(dimx))
             ALLOCATE(chieeETG_SI(dimx))
             ALLOCATE(veneETG_SI(dimx))
             ALLOCATE(veceETG_SI(dimx))

             ALLOCATE(chieeITG_GB(dimx))
             ALLOCATE(veneITG_GB(dimx))
             ALLOCATE(veceITG_GB(dimx))
             ALLOCATE(chieeTEM_GB(dimx))
             ALLOCATE(veneTEM_GB(dimx))
             ALLOCATE(veceTEM_GB(dimx))
             ALLOCATE(chieeETG_GB(dimx))
             ALLOCATE(veneETG_GB(dimx))
             ALLOCATE(veceETG_GB(dimx))

          ENDIF

       ENDIF
    ENDIF
    ALLOCATE(epf_cm(dimx,dimn))
    ALLOCATE(eef_cm(dimx,dimn))

    ALLOCATE(ipf_SI(dimx,nions))  
    ALLOCATE(ipf_GB(dimx,nions))
    ALLOCATE(ief_SI(dimx,nions))
    ALLOCATE(ief_GB(dimx,nions))
    ALLOCATE(ivf_SI(dimx,nions))
    ALLOCATE(ivf_GB(dimx,nions))
    IF (separateflux == 1) THEN
       ALLOCATE(iefITG_SI(dimx,nions))
       ALLOCATE(iefTEM_SI(dimx,nions))
       ALLOCATE(ipfITG_SI(dimx,nions))
       ALLOCATE(ipfTEM_SI(dimx,nions))
       ALLOCATE(ivfITG_SI(dimx,nions))
       ALLOCATE(ivfTEM_SI(dimx,nions))

       ALLOCATE(iefITG_GB(dimx,nions))
       ALLOCATE(iefTEM_GB(dimx,nions))
       ALLOCATE(ipfITG_GB(dimx,nions))
       ALLOCATE(ipfTEM_GB(dimx,nions))
       ALLOCATE(ivfITG_GB(dimx,nions))
       ALLOCATE(ivfTEM_GB(dimx,nions))
    ENDIF

    IF (phys_meth /= 0) THEN
       ALLOCATE(dfi_SI(dimx,nions))
       ALLOCATE(vti_SI(dimx,nions))
       ALLOCATE(vci_SI(dimx,nions))
       ALLOCATE(vri_SI(dimx,nions))

       ALLOCATE(dfi_GB(dimx,nions))
       ALLOCATE(vti_GB(dimx,nions))
       ALLOCATE(vci_GB(dimx,nions))
       ALLOCATE(vri_GB(dimx,nions))

       ALLOCATE(cki(dimx,nions))

       IF (separateflux == 1) THEN
          ALLOCATE(dfiITG_SI(dimx,nions))
          ALLOCATE(vtiITG_SI(dimx,nions))
          ALLOCATE(vciITG_SI(dimx,nions))
          ALLOCATE(vriITG_SI(dimx,nions))
          ALLOCATE(dfiTEM_SI(dimx,nions))
          ALLOCATE(vtiTEM_SI(dimx,nions))
          ALLOCATE(vciTEM_SI(dimx,nions))
          ALLOCATE(vriTEM_SI(dimx,nions))

          ALLOCATE(dfiITG_GB(dimx,nions))
          ALLOCATE(vtiITG_GB(dimx,nions))
          ALLOCATE(vciITG_GB(dimx,nions))
          ALLOCATE(vriITG_GB(dimx,nions))
          ALLOCATE(dfiTEM_GB(dimx,nions))
          ALLOCATE(vtiTEM_GB(dimx,nions))
          ALLOCATE(vciTEM_GB(dimx,nions))
          ALLOCATE(vriTEM_GB(dimx,nions))

       ENDIF

       IF (phys_meth == 2) THEN
          ALLOCATE(veni_SI(dimx,nions))
          ALLOCATE(chiei_SI(dimx,nions))
          ALLOCATE(veci_SI(dimx,nions))
          ALLOCATE(veri_SI(dimx,nions))
          ALLOCATE(veni_GB(dimx,nions))
          ALLOCATE(chiei_GB(dimx,nions))
          ALLOCATE(veci_GB(dimx,nions))
          ALLOCATE(veri_GB(dimx,nions))
          ALLOCATE(ceki(dimx,nions))
          IF (separateflux == 1) THEN
             ALLOCATE(chieiITG_SI(dimx,nions))
             ALLOCATE(veniITG_SI(dimx,nions))
             ALLOCATE(veciITG_SI(dimx,nions))
             ALLOCATE(veriITG_SI(dimx,nions))
             ALLOCATE(chieiTEM_SI(dimx,nions))
             ALLOCATE(veniTEM_SI(dimx,nions))
             ALLOCATE(veciTEM_SI(dimx,nions))
             ALLOCATE(veriTEM_SI(dimx,nions))

             ALLOCATE(chieiITG_GB(dimx,nions))
             ALLOCATE(veniITG_GB(dimx,nions))
             ALLOCATE(veciITG_GB(dimx,nions))
             ALLOCATE(veriITG_GB(dimx,nions))
             ALLOCATE(chieiTEM_GB(dimx,nions))
             ALLOCATE(veniTEM_GB(dimx,nions))
             ALLOCATE(veciTEM_GB(dimx,nions))
             ALLOCATE(veriTEM_GB(dimx,nions))
          ENDIF
       ENDIF
    ENDIF
    ALLOCATE(ipf_cm(dimx,dimn,nions))
    ALLOCATE(ief_cm(dimx,dimn,nions))
    ALLOCATE(ivf_cm(dimx,dimn,nions))
    ALLOCATE( krmmuITG (dimx) )
    ALLOCATE( krmmuETG (dimx) )
    ALLOCATE( kperp2 (dimx, dimn) )
    ALLOCATE( modewidth (dimx, dimn) )
    ALLOCATE( modeshift (dimx, dimn) )
    ALLOCATE( distan (dimx, dimn) )
    ALLOCATE(solflu (dimx, dimn))

    ALLOCATE( ntor (dimx, dimn) )
    ALLOCATE( modeflag (dimx) )
    ALLOCATE( Nustar (dimx) )
    ALLOCATE( Zeffx (dimx) )
    ALLOCATE(phi(dimx,ntheta))
    ALLOCATE(npol(dimx,ntheta,nions))
    ALLOCATE(ecoefs(dimx,0:nions,numecoefs)) !includes electrons
    ALLOCATE(cftrans(dimx,nions,numicoefs)) ! only for ions

    ALLOCATE( sol (dimx, dimn, numsols) )
    ALLOCATE( fdsol (dimx, dimn, numsols) )
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

    IF (phys_meth /= 0) THEN
       ALLOCATE( Lcircgte (dimx, dimn, numsols) )
       ALLOCATE( Lpieggte (dimx, dimn, numsols) )
       ALLOCATE( Lcircgne (dimx, dimn, numsols) )
       ALLOCATE( Lpieggne (dimx, dimn, numsols) )
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

    DEALLOCATE(dummyn,dummyx,dummyxnions,kthetarhostmp,xtmp,rhotmp,Rotmp,Rmintmp,Botmp,qxtmp,smagtmp,alphaxtmp,&
         & Machtortmp,Autortmp,Machpartmp,Aupartmp,gammaEtmp,Textmp,Nextmp,Atetmp,Anetmp,anisetmp,danisedrtmp,Tixtmp,&
         & ninormtmp,Atitmp,Anitmp,anistmp,danisdrtmp,ion_typetmp,Aitmp,Zitmp,oldsoltmp,oldfdsoltmp)

  END SUBROUTINE data_init

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

    DEALLOCATE(solflu_SI)
    DEALLOCATE(solflu_GB)
    DEALLOCATE(gam_SI)
    DEALLOCATE(gam_GB)
    DEALLOCATE(ome_SI)
    DEALLOCATE(ome_GB)
    DEALLOCATE(epf_SI)
    DEALLOCATE(epf_GB)
    DEALLOCATE(eef_SI)
    DEALLOCATE(eef_GB)
    DEALLOCATE(epf_cm)
    DEALLOCATE(eef_cm)

    DEALLOCATE(ipf_SI)
    DEALLOCATE(ipf_GB)
    DEALLOCATE(ief_SI)
    DEALLOCATE(ief_GB)
    DEALLOCATE(ivf_SI)
    DEALLOCATE(ivf_GB)
    DEALLOCATE(eefETG_SI)
    DEALLOCATE(eefETG_GB)

    IF (separateflux == 1) THEN
       DEALLOCATE(eefTEM_SI)
       DEALLOCATE(eefITG_SI)
       DEALLOCATE(epfITG_SI)
       DEALLOCATE(epfTEM_SI)

       DEALLOCATE(iefITG_SI)
       DEALLOCATE(iefTEM_SI)
       DEALLOCATE(ipfITG_SI)
       DEALLOCATE(ipfTEM_SI)
       DEALLOCATE(ivfITG_SI)
       DEALLOCATE(ivfTEM_SI)

       DEALLOCATE(eefTEM_GB)
       DEALLOCATE(eefITG_GB)
       DEALLOCATE(epfTEM_GB)
       DEALLOCATE(epfITG_GB)

       DEALLOCATE(iefITG_GB)
       DEALLOCATE(iefTEM_GB)
       DEALLOCATE(ipfITG_GB)
       DEALLOCATE(ipfTEM_GB)
       DEALLOCATE(ivfITG_GB)
       DEALLOCATE(ivfTEM_GB)
    ENDIF

    IF (phys_meth /= 0) THEN
       DEALLOCATE(dfe_SI)
       DEALLOCATE(vte_SI)
       DEALLOCATE(vce_SI)
       DEALLOCATE(cke)
       DEALLOCATE(dfi_SI)
       DEALLOCATE(vti_SI)
       DEALLOCATE(vci_SI)
       DEALLOCATE(vri_SI)
       DEALLOCATE(cki)

       DEALLOCATE(dfe_GB)
       DEALLOCATE(vte_GB)
       DEALLOCATE(vce_GB)
       DEALLOCATE(dfi_GB)
       DEALLOCATE(vti_GB)
       DEALLOCATE(vci_GB)
       DEALLOCATE(vri_GB)


       IF (separateflux == 1) THEN
          DEALLOCATE(dfeITG_SI)
          DEALLOCATE(vteITG_SI)
          DEALLOCATE(vceITG_SI)
          DEALLOCATE(dfeTEM_SI)
          DEALLOCATE(vteTEM_SI)
          DEALLOCATE(vceTEM_SI)
          DEALLOCATE(dfiITG_SI)
          DEALLOCATE(vtiITG_SI)
          DEALLOCATE(vciITG_SI)
          DEALLOCATE(vriITG_SI)
          DEALLOCATE(dfiTEM_SI)
          DEALLOCATE(vtiTEM_SI)
          DEALLOCATE(vciTEM_SI)
          DEALLOCATE(vriTEM_SI)

          DEALLOCATE(dfeITG_GB)
          DEALLOCATE(vteITG_GB)
          DEALLOCATE(vceITG_GB)
          DEALLOCATE(dfeTEM_GB)
          DEALLOCATE(vteTEM_GB)
          DEALLOCATE(vceTEM_GB)
          DEALLOCATE(dfiITG_GB)
          DEALLOCATE(vtiITG_GB)
          DEALLOCATE(vciITG_GB)
          DEALLOCATE(vriITG_GB)
          DEALLOCATE(dfiTEM_GB)
          DEALLOCATE(vtiTEM_GB)
          DEALLOCATE(vciTEM_GB)
          DEALLOCATE(vriTEM_GB)

       ENDIF


       IF (phys_meth == 2) THEN
          DEALLOCATE(vene_SI)
          DEALLOCATE(chiee_SI)
          DEALLOCATE(vece_SI)
          DEALLOCATE(ceke)
          DEALLOCATE(veni_SI)
          DEALLOCATE(chiei_SI)
          DEALLOCATE(veci_SI)
          DEALLOCATE(veri_SI)
          DEALLOCATE(ceki)

          IF (separateflux == 1) THEN
             DEALLOCATE(chieeITG_SI)
             DEALLOCATE(veneITG_SI)
             DEALLOCATE(veceITG_SI)
             DEALLOCATE(chieeETG_SI)
             DEALLOCATE(veneETG_SI)
             DEALLOCATE(veceETG_SI)
             DEALLOCATE(chieeTEM_SI)
             DEALLOCATE(veneTEM_SI)
             DEALLOCATE(veceTEM_SI)
             DEALLOCATE(chieiITG_SI)
             DEALLOCATE(veniITG_SI)
             DEALLOCATE(veciITG_SI)
             DEALLOCATE(veriITG_SI)
             DEALLOCATE(chieiTEM_SI)
             DEALLOCATE(veniTEM_SI)
             DEALLOCATE(veciTEM_SI)
             DEALLOCATE(veriTEM_SI)

             DEALLOCATE(chieeITG_GB)
             DEALLOCATE(veneITG_GB)
             DEALLOCATE(veceITG_GB)
             DEALLOCATE(chieeETG_GB)
             DEALLOCATE(veneETG_GB)
             DEALLOCATE(veceETG_GB)

             DEALLOCATE(chieeTEM_GB)
             DEALLOCATE(veneTEM_GB)
             DEALLOCATE(veceTEM_GB)
             DEALLOCATE(chieiITG_GB)
             DEALLOCATE(veniITG_GB)
             DEALLOCATE(veciITG_GB)
             DEALLOCATE(veriITG_GB)
             DEALLOCATE(chieiTEM_GB)
             DEALLOCATE(veniTEM_GB)
             DEALLOCATE(veciTEM_GB)
             DEALLOCATE(veriTEM_GB)

          ENDIF

       ENDIF
    ENDIF

    DEALLOCATE(ipf_cm)
    DEALLOCATE(ief_cm)
    DEALLOCATE(ivf_cm)
    DEALLOCATE( krmmuITG)
    DEALLOCATE( krmmuETG)
    DEALLOCATE( kperp2)

    DEALLOCATE(modewidth)
    DEALLOCATE(modeshift)
    DEALLOCATE(modeflag)
    DEALLOCATE(Nustar)
    DEALLOCATE(Zeffx)
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
    DEALLOCATE(Lvpiegi)
    DEALLOCATE(Lvcirci)

    IF (phys_meth /= 0) THEN
       DEALLOCATE(Lcircgte)
       DEALLOCATE(Lpieggte)
       DEALLOCATE(Lcircgne)
       DEALLOCATE(Lpieggne)
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

    CHARACTER(len=20) :: fmtxrow,fmtecoef,fmtcftrans
    INTEGER :: i,j,k,l, fileno
    DOUBLE PRECISION :: time5

    WRITE(fmtxrow,'(A,I0,A)') '(',dimx,'G15.7)'
    WRITE(fmtecoef,'(A,I0, A)') '(',numecoefs,'G15.7)'
    WRITE(fmtcftrans,'(A)') '(7G15.7)'

    fileno = 0
    IF (write_primi == 1) THEN
      primitivedir='output/primitive/'
      myfmt='G16.7E3'
      CALL writevar(primitivedir // 'solflu.dat', solflu, myfmt, fileno)
      !WRITE(stdout,"(A,I10,A,I10)") 'rank: ', myrank, '01.dat ', (time5-omp_get_wtime())
      fileno=fileno+1
      CALL writevar(primitivedir // 'kymaxITG.dat', krmmuITG, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'kymaxETG.dat', krmmuETG, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'distan.dat', distan, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'kperp2.dat', kperp2, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'modewidth.dat', modewidth, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'modeshift.dat', modeshift, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'ntor.dat', ntor, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'sol.dat', sol, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'fdsol.dat', fdsol, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lcirce.dat', Lcirce, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lpiege.dat', Lpiege, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lecirce.dat', Lecirce, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lepiege.dat', Lepiege, myfmt, fileno)
      fileno=fileno+1

      IF (phys_meth /= 0) THEN
         CALL writevar(primitivedir // 'Lcircgne.dat', Lcircgne, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpieggne.dat', Lpieggne, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lcircgte.dat', Lcircgte, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpieggte.dat', Lpieggte, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lcircce.dat', Lcircce, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpiegce.dat', Lpiegce, myfmt, fileno)
         fileno=fileno+1

         IF (phys_meth == 2) THEN
            CALL writevar(primitivedir // 'Lecircgne.dat', Lecircgne, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepieggne.dat', Lepieggne, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lecircgte.dat', Lecircgte, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepieggte.dat', Lepieggte, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lecircce.dat', Lecircce, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepiegce.dat', Lepiegce, myfmt, fileno)
            fileno=fileno+1
         ENDIF
      ENDIF

      CALL writevar(primitivedir // 'Lcirci.dat', Lcirci, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lpiegi.dat', Lpiegi, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lcirci.dat', Lcirci, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lpiegi.dat', Lpiegi, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lecirci.dat', Lecirci, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lepiegi.dat', Lepiegi, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lvcirci.dat', Lvcirci, myfmt, fileno)
      fileno=fileno+1
      CALL writevar(primitivedir // 'Lvpiegi.dat', Lvpiegi, myfmt, fileno)
      fileno=fileno+1

      IF (phys_meth /= 0) THEN
         CALL writevar(primitivedir // 'Lcircgni.dat', Lcircgni, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpieggni.dat', Lpieggni, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lcircgui.dat', Lcircgui, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpieggui.dat', Lpieggui, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lcircgti.dat', Lcircgti, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpieggti.dat', Lpieggti, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lcircci.dat', Lcircci, myfmt, fileno)
         fileno=fileno+1
         CALL writevar(primitivedir // 'Lpiegci.dat', Lpiegci, myfmt, fileno)
         fileno=fileno+1

         IF (phys_meth == 2) THEN
            CALL writevar(primitivedir // 'Lecircgni.dat', Lecircgni, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepieggni.dat', Lepieggni, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lecircgui.dat', Lecircgui, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepieggui.dat', Lepieggui, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lecircgti.dat', Lecircgti, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepieggti.dat', Lepieggti, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lecircci.dat', Lecircci, myfmt, fileno)
            fileno=fileno+1
            CALL writevar(primitivedir // 'Lepiegci.dat', Lepiegci, myfmt, fileno)
            fileno=fileno+1
         ENDIF
      ENDIF
    ENDIF

    outputdir = 'debug/'

    CALL writevar(outputdir // 'modeflag.dat', modeflag, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'phi.dat', TRANSPOSE(phi), myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'Nustar.dat', Nustar, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'Zeff.dat', Zeffx, myfmt, fileno)
    fileno=fileno+1

    outputdir = 'output/'

    CALL writevar(outputdir // 'npol.dat', npol, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'ecoefs.dat', ecoefs, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'cftrans.dat', cftrans, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'gam_GB.dat', gam_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'ome_GB.dat', ome_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'gam_SI.dat', gam_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'ome_SI.dat', ome_SI, myfmt, fileno)
    fileno=fileno+1

    IF (phys_meth /= 0) THEN
       CALL writevar(outputdir // 'cke.dat', cke, myfmt, fileno)
       fileno=fileno+1

       CALL writevar(outputdir // 'dfe_SI.dat', dfe_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vte_SI.dat', vte_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vce_SI.dat', vce_SI, myfmt, fileno)
       fileno=fileno+1

       CALL writevar(outputdir // 'dfe_GB.dat', dfe_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vte_GB.dat', vte_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vce_GB.dat', vce_GB, myfmt, fileno)
       fileno=fileno+1

       IF (separateflux == 1) THEN
          CALL writevar(outputdir // 'dfeITG_SI.dat', dfeITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vteITG_SI.dat', vteITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vceITG_SI.dat', vceITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'dfeTEM_SI.dat', dfeTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vteTEM_SI.dat', vteTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vceTEM_SI.dat', vceTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'dfeITG_GB.dat', dfeITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vteITG_GB.dat', vteITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vceITG_GB.dat', vceITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'dfeTEM_GB.dat', dfeTEM_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vteTEM_GB.dat', vteTEM_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vceTEM_GB.dat', vceTEM_GB, myfmt, fileno)
          fileno=fileno+1


       ENDIF

       CALL writevar(outputdir // 'cki.dat', cki, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'dfi_SI.dat', dfi_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vti_SI.dat', vti_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vri_SI.dat', vri_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vci_SI.dat', vci_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'dfi_GB.dat', dfi_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vti_GB.dat', vti_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vri_GB.dat', vri_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vci_GB.dat', vci_GB, myfmt, fileno)
       fileno=fileno+1

       IF (separateflux == 1) THEN
          CALL writevar(outputdir // 'dfiITG_SI.dat', dfiITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vtiITG_SI.dat', vtiITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vciITG_SI.dat', vciITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vriITG_SI.dat', vriITG_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'dfiTEM_SI.dat', dfiTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vtiTEM_SI.dat', vtiTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vciTEM_SI.dat', vciTEM_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vriTEM_SI.dat', vriTEM_SI, myfmt, fileno)
          fileno=fileno+1

          CALL writevar(outputdir // 'dfiITG_GB.dat', dfiITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vtiITG_GB.dat', vtiITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vciITG_GB.dat', vciITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vriITG_GB.dat', vriITG_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'dfiTEM_GB.dat', dfiTEM_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vtiTEM_GB.dat', vtiTEM_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vciTEM_GB.dat', vciTEM_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vriTEM_GB.dat', vriTEM_GB, myfmt, fileno)
          fileno=fileno+1
       ENDIF

       IF (phys_meth == 2) THEN
          CALL writevar(outputdir // 'ceke.dat', ceke, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'ceki.dat', ceki, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vene_SI.dat', vene_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'chiee_SI.dat', chiee_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vece_SI.dat', vece_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veni_SI.dat', veni_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veri_SI.dat', veri_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'chiei_SI.dat', chiei_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veci_SI.dat', veci_SI, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vene_GB.dat', vene_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'chiee_GB.dat', chiee_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'vece_GB.dat', vece_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veni_GB.dat', veni_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veri_GB.dat', veri_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'chiei_GB.dat', chiei_GB, myfmt, fileno)
          fileno=fileno+1
          CALL writevar(outputdir // 'veci_GB.dat', veci_GB, myfmt, fileno)
          fileno=fileno+1

          IF (separateflux == 1) THEN
             CALL writevar(outputdir // 'chieiITG_SI.dat', chieiITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veniITG_SI.dat', veniITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veciITG_SI.dat', veciITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veriITG_SI.dat', veriITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieiTEM_SI.dat', chieiTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veniTEM_SI.dat', veniTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veciTEM_SI.dat', veciTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veriTEM_SI.dat', veriTEM_SI, myfmt, fileno)
             fileno=fileno+1

             CALL writevar(outputdir // 'chieiITG_GB.dat', chieiITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veniITG_GB.dat', veniITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veciITG_GB.dat', veciITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veriITG_GB.dat', veriITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieiTEM_GB.dat', chieiTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veniTEM_GB.dat', veniTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veciTEM_GB.dat', veciTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veriTEM_GB.dat', veriTEM_GB, myfmt, fileno)
             fileno=fileno+1

             CALL writevar(outputdir // 'chieeETG_SI.dat', chieeETG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneETG_SI.dat', veneETG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceETG_SI.dat', veceETG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieeETG_GB.dat', chieeETG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneETG_GB.dat', veneETG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceETG_GB.dat', veceETG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieeTEM_SI.dat', chieeTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneTEM_SI.dat', veneTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceTEM_SI.dat', veceTEM_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieeTEM_GB.dat', chieeTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneTEM_GB.dat', veneTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceTEM_GB.dat', veceTEM_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieeITG_SI.dat', chieeITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneITG_SI.dat', veneITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceITG_SI.dat', veceITG_SI, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'chieeITG_GB.dat', chieeITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veneITG_GB.dat', veneITG_GB, myfmt, fileno)
             fileno=fileno+1
             CALL writevar(outputdir // 'veceITG_GB.dat', veceITG_GB, myfmt, fileno)
             fileno=fileno+1
          ENDIF
       ENDIF
    ENDIF

    CALL writevar(outputdir // 'pfe_SI.dat', epf_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'pfe_GB.dat', epf_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'pfe_cm.dat', epf_cm, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'efe_SI.dat', eef_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'efe_GB.dat', eef_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'efe_cm.dat', eef_cm, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'efeETG_SI.dat', eefETG_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'efeETG_GB.dat', eefETG_GB, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'pfi_SI.dat', ipf_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'pfi_GB.dat', ipf_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'pfi_cm.dat', ipf_cm, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'efi_SI.dat', ief_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'efi_GB.dat', ief_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'efi_cm.dat', ief_cm, myfmt, fileno)
    fileno=fileno+1

    CALL writevar(outputdir // 'vfi_SI.dat', ivf_SI, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'vfi_GB.dat', ivf_GB, myfmt, fileno)
    fileno=fileno+1
    CALL writevar(outputdir // 'vfi_cm.dat', ivf_cm, myfmt, fileno)
    fileno=fileno+1

    IF (separateflux==1) THEN
       CALL writevar(outputdir // 'efiITG_SI.dat', iefITG_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'efiTEM_SI.dat', iefTEM_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfiITG_SI.dat', ipfITG_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfiTEM_SI.dat', ipfTEM_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vfiITG_SI.dat', ivfITG_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vfiTEM_SI.dat', ivfTEM_SI, myfmt, fileno)
       fileno=fileno+1

       CALL writevar(outputdir // 'efiITG_GB.dat', iefITG_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'efiTEM_GB.dat', iefTEM_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfiITG_GB.dat', ipfITG_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfiTEM_GB.dat', ipfTEM_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vfiITG_GB.dat', ivfITG_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'vfiTEM_GB.dat', ivfTEM_GB, myfmt, fileno)
       fileno=fileno+1

       CALL writevar(outputdir // 'efeITG_SI.dat', eefITG_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'efeTEM_SI.dat', eefTEM_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfeITG_SI.dat', epfITG_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfeTEM_SI.dat', epfTEM_SI, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'efeITG_GB.dat', eefITG_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'efeTEM_GB.dat', eefTEM_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfeITG_GB.dat', epfITG_GB, myfmt, fileno)
       fileno=fileno+1
       CALL writevar(outputdir // 'pfeTEM_GB.dat', epfTEM_GB, myfmt, fileno)
       fileno=fileno+1

    ENDIF

  END SUBROUTINE outputascii

END PROGRAM qlk_standalone
