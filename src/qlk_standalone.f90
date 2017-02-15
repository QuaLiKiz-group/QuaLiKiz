PROGRAM qlk_standalone
  !Standalone driver for qualikiz. Detailed code info in call_qualikiz subroutine header

  USE kind
  USE mod_io_management !MPI is included in this module

  IMPLICIT NONE

  ! INTERFACE WITH EXTERNAL QUALIKIZ SUBROUTINE
  INTERFACE 
     SUBROUTINE qualikiz(dimxin, rhoin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, kthetarhosin, & !general param
          & xin, Roin, Rminin, R0in, Boin, qxin, smagin, alphaxin, & !geometry input
          & el_typein, Texin, Nexin, Atein, Anein, anisein, danisedrin, & !electron input
          & ion_typein, Aiin, Ziin, Tixin, ninormin, Atiin, Aniin, anisin, danisdrin, & !ion input
          & Machtorin, Autorin, Machparin, Auparin, gammaEin, & !rotation input
          & maxrunsin, maxptsin, relacc1in, relacc2in, timeoutin, ETGmultin, collmultin, & !code specific input
          & epf_SIout,epfETG_SIout,eef_SIout,eefETG_SIout,evf_SIout,ipf_SIout,ief_SIout,ivf_SIout, & ! Non optional outputs
          & solflu_SIout, solflu_GBout, gam_SIout,gam_GBout,ome_SIout,ome_GBout, & !growth rate and frequency output
          & epf_GBout,eef_GBout, evf_GBout, dfe_SIout,vte_SIout,vre_SIout,vce_SIout,epf_cmout,eef_cmout,evf_cmout,ckeout, & !electron flux outputs
          & ipf_GBout,ief_GBout, ivf_GBout, dfi_SIout,vti_SIout,vri_SIout,vci_SIout,ipf_cmout,ief_cmout,ivf_cmout,ckiout, & !ion flux outputs
          & vene_SIout,chiee_SIout,vere_SIout,vece_SIout, cekeout, veni_SIout,chiei_SIout,veci_SIout,veri_SIout,cekiout, & !heat pinch outputs
          & modeflagout, & ! flags type of modes in output per radial position
          & phiout, npolout, ecoefsout, cftransout, &  ! poloidal asymmetry outputs for heavy impurities
          & solfluout, modewidthout, modeshiftout, distanout, ntorout, solout, fdsolout,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
          & kperp2out,krmmuITGout,krmmuETGout, &
          & Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lvcirceout, Lvpiegeout, Lcircgteout, Lpieggteout, Lcircgneout, Lpieggneout, Lcircgueout, Lpieggueout, Lcircceout, Lpiegceout, & 
          & Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout, &
          & Lecircgteout, Lepieggteout, Lecircgneout, Lepieggneout, Lecircgueout, Lepieggueout, Lecircceout, Lepiegceout, Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout,&
          oldsolin, oldfdsolin, runcounterin,&
          rhominin,rhomaxin)

       INTEGER, INTENT(IN) :: dimxin, dimnin, nionsin, numsolsin, phys_methin, coll_flagin, rot_flagin, verbosein, el_typein
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

       ! List of output variables: 
       INTEGER, PARAMETER :: ntheta = 64
       INTEGER, PARAMETER :: numecoefs = 13

       ! growth rate and frequency outputs
       COMPLEX, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT) :: solflu_SIout, solflu_GBout, solfluout
       REAL, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: gam_SIout,gam_GBout,ome_SIout,ome_GBout  

       ! final output arrays following saturation rule
       REAL, DIMENSION(dimxin), INTENT(OUT)  :: epf_SIout,epfETG_SIout,eef_SIout,eefETG_SIout,evf_SIout
       REAL, DIMENSION(dimxin,nionsin), INTENT(OUT)  :: ipf_SIout,ief_SIout,ivf_SIout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: epf_GBout,eef_GBout, evf_GBout,dfe_SIout, vte_SIout, vre_SIout,vce_SIout, ckeout, modeflagout
       REAL, DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: vene_SIout, chiee_SIout, vere_SIout,vece_SIout, cekeout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: ipf_GBout,ief_GBout, ivf_GBout,dfi_SIout, vti_SIout, vri_SIout,vci_SIout, ckiout
       REAL, DIMENSION(dimxin,nionsin), OPTIONAL, INTENT(OUT)  :: veni_SIout, chiei_SIout, veci_SIout, veri_SIout,cekiout
       REAL, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  ::  epf_cmout, eef_cmout,evf_cmout
       REAL, DIMENSION(dimxin,dimnin,nionsin), OPTIONAL, INTENT(OUT) :: ipf_cmout,ief_cmout,ivf_cmout
       REAL, DIMENSION(dimxin,ntheta), OPTIONAL, INTENT(OUT)  ::  phiout
       REAL, DIMENSION(dimxin,ntheta,nionsin), OPTIONAL, INTENT(OUT)  ::  npolout
       REAL, DIMENSION(dimxin,0:nionsin,numecoefs), OPTIONAL, INTENT(OUT)  ::  ecoefsout
       REAL, DIMENSION(dimxin,nionsin,6), OPTIONAL, INTENT(OUT)  ::  cftransout

       ! optional output arrays from which the saturation rule can be calculated without rerunning dispersion relation solver
       REAL , DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: distanout,ntorout,kperp2out
       REAL , DIMENSION(dimxin), OPTIONAL, INTENT(OUT)  :: krmmuITGout,krmmuETGout
       COMPLEX, DIMENSION(dimxin,dimnin), OPTIONAL, INTENT(OUT)  :: modewidthout, modeshiftout
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: solout,fdsolout
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirceout, Lpiegeout, Lecirceout, Lepiegeout, Lvcirceout, Lvpiegeout, Lcircgteout, Lpieggteout,  Lcircgneout, Lpieggneout,  Lcircgueout, Lpieggueout, Lcircceout, Lpiegceout
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgteout, Lepieggteout,  Lecircgneout, Lepieggneout,  Lecircgueout, Lepieggueout, Lecircceout, Lepiegceout
       COMPLEX, DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lcirciout, Lpiegiout, Lecirciout, Lepiegiout, Lvcirciout, Lvpiegiout, Lcircgtiout, Lpieggtiout, Lcircgniout, Lpieggniout, Lcircguiout, Lpiegguiout, Lcircciout, Lpiegciout
       COMPLEX, DIMENSION(dimxin,dimnin,nionsin,numsolsin), OPTIONAL, INTENT(OUT)  :: Lecircgtiout, Lepieggtiout, Lecircgniout, Lepieggniout, Lecircguiout, Lepiegguiout, Lecircciout, Lepiegciout

       !optional input arrays for going directly to newton solver
       INTEGER, OPTIONAL, INTENT(IN)  :: runcounterin
       COMPLEX, DIMENSION(dimxin,dimnin,numsolsin), OPTIONAL, INTENT(IN)  :: oldsolin, oldfdsolin

     END SUBROUTINE qualikiz
  END INTERFACE

  !Local data dictionary.

  !Time measuring variables
  REAL(kind=DBL) :: cputime1, cputime2, tpstot
  INTEGER :: time1, time2, timetot, freq, cputimetot
  CHARACTER(len=20) :: myfmt

  !MPI variables:
  INTEGER :: ierror, nproc, myrank

  !Input variables into qualikiz subroutine
  INTEGER :: dimx, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, el_type
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ion_type
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: x, rho, Ro, Rmin, Bo, qx, smag, alphax, kthetarhos
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: Tex, Nex, Ate, Ane, anise, danisedr
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: Tix, ninorm, Ati, Ani, anis, danisdr, Ai, Zi
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: Machtor, Autor, Machpar, Aupar, gammaE
  REAL(KIND=DBL) :: relacc1, relacc2, ETGmult, collmult, timeout, R0
  INTEGER :: maxpts,maxruns

  ! Output arrays. The 3 dimensions are 'radial grid', 'kthetarhos grid', 'number of modes'
  REAL(KIND=DBL) , DIMENSION(:), ALLOCATABLE :: krmmuITG,krmmuETG
  REAL(KIND=DBL) , DIMENSION(:,:), ALLOCATABLE :: distan,FLRec,FLRep,ntor,kperp2
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: gamma, Ladia, FLRip, FLRic
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: modewidth, modeshift
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: ommax, solflu
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: sol, fdsol
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: Lcirce, Lpiege, Lecirce, Lepiege, Lvcirce, Lvpiege, Lcircgte, Lpieggte,  Lcircgne, Lpieggne,  Lcircgue, Lpieggue, Lcircce, Lpiegce
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: Lecircgte, Lepieggte,  Lecircgne, Lepieggne,  Lecircgue, Lepieggue, Lecircce, Lepiegce
  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), ALLOCATABLE :: Lcirci, Lpiegi, Lecirci, Lepiegi, Lvcirci, Lvpiegi, Lcircgti, Lpieggti, Lcircgni, Lpieggni, Lcircgui, Lpieggui, Lcircci, Lpiegci
  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), ALLOCATABLE :: Lecircgti, Lepieggti, Lecircgni, Lepieggni, Lecircgui, Lepieggui, Lecircci, Lepiegci

  ! Old solution for non-reset runs
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldrsol,oldisol,oldrfdsol,oldifdsol
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: oldsol,oldfdsol

  ! Final output arrays following saturation rule. These can be printed as ASCII output
  COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: solflu_SI, solflu_GB
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: gam_SI,gam_GB,ome_SI,ome_GB
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: epf_SI,epfETG_SI,epf_GB,eef_SI,eefETG_SI,eef_GB, evf_SI,evf_GB,dfe_SI,vte_SI,vce_SI,vre_SI,cke,modeflag
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: vene_SI,chiee_SI,vece_SI,vere_SI,ceke
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: ipf_SI,ipf_GB,ief_SI,ief_GB, ivf_SI,ivf_GB,dfi_SI,vti_SI,vri_SI,vci_SI,cki,eef_cm,epf_cm,evf_cm
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: veni_SI,chiei_SI,veci_SI,ceki,veri_SI
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: ipf_cm,ief_cm,ivf_cm

  ! Poloidal asymmetry variables
  ! Ion density along field line, used was asymmetries are present. 
  ! 1st dimension is radius. 2nd is theta (ntheta points from [0,pi]), 3rd is ion
  REAL(KIND=DBL), DIMENSION(:,:,:), ALLOCATABLE :: npol, ecoefs, cftrans
  REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: phi
  INTEGER, PARAMETER :: ntheta = 64
  INTEGER, PARAMETER :: numecoefs = 13

  REAL, PARAMETER :: epsD = 1d-14
  INTEGER :: runcounter ! used for counting runs inside integrated modelling applications for deciding to recalculate all or just jump to newton based on old solutions

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

  IF (myrank==0) THEN
     WRITE(stdout,*) ' _________________________________________________________________________________ '
     WRITE(stdout,*) '                                     QUALIKIZ 2.0 '
     WRITE(stdout,*) '  gyrokinetic calculation of linear growth rates and quasilinear transport fluxes  '
     WRITE(stdout,*) ' _________________________________________________________________________________ '
     WRITE(stdout,*) ' '
     WRITE(stdout,'(1X,1A,I4,1A,/)') 'Executed in parallel on: ',nproc,' processors'
  ENDIF

!!!! CALL THE CALCULATION PHASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (ALLOCATED(oldsol)) THEN !Call with optional old solution input

     IF (phys_meth == 0) THEN
        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, evf_cmout=evf_cm, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, &
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, &
             & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver

     ENDIF

     IF (phys_meth == 1) THEN

        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose,  kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vre_SIout=vre_SI, vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm,evf_cmout=evf_cm, ckeout=cke, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, Lcircgteout=Lcircgte, &
             & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircgueout=Lcircgue, Lpieggueout=Lpieggue, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
             & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &         
             & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver
     ENDIF

     IF (phys_meth == 2) THEN

        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vre_SIout=vre_SI, vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm,evf_cmout=evf_cm, ckeout=cke, & !optional electron flux outputs
             & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vere_SIout=vere_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
             & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veri_SIout=veri_SI,veci_SIout=veci_SI,cekiout=ceki, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, Lcircgteout=Lcircgte, &
             & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircgueout=Lcircgue, Lpieggueout=Lpieggue, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
             & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircgueout=Lecircgue, Lepieggueout=Lepieggue, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
             & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
             & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci, &
             & oldsolin=oldsol,oldfdsolin=oldfdsol,runcounterin=runcounter) ! optional inputs for jumping straight to newton solver

     ENDIF


  ELSE !Don't call with optional old solution input
     IF (phys_meth == 0) THEN
        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, epf_cmout=epf_cm,eef_cmout=eef_cm, evf_cmout=evf_cm, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, ipf_cmout=ipf_cm,ief_cmout=ief_cm, ivf_cmout=ivf_cm, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, &
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi)

     ENDIF

     IF (phys_meth == 1) THEN
        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vre_SIout=vre_SI, vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm,evf_cmout=evf_cm, ckeout=cke, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, Lcircgteout=Lcircgte, &
             & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircgueout=Lcircgue, Lpieggueout=Lpieggue, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
             & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci)

     ENDIF

     IF (phys_meth == 2) THEN
        CALL qualikiz(dimx, rho, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, verbose, kthetarhos, & !general param
             & x, Ro, Rmin, R0, Bo, qx, smag, alphax, & !geometry
             & el_type, Tex, Nex, Ate, Ane, anise, danisedr, & !electrons
             & ion_type, Ai, Zi, Tix, ninorm, Ati, Ani, anis, danisdr, & !ions
             & Machtor, Autor, Machpar, Aupar, gammaE, & !rotation input
             & maxruns, maxpts, relacc1, relacc2, timeout, ETGmult, collmult,  & !code specific inputs
             & epf_SI,epfETG_SI,eef_SI,eefETG_SI,evf_SI,ipf_SI,ief_SI, ivf_SI, & ! Non optional outputs
             & solflu_SIout=solflu_SI, solflu_GBout=solflu_GB, gam_SIout=gam_SI,gam_GBout=gam_GB,ome_SIout=ome_SI,ome_GBout=ome_GB, & !optional growth rate and frequency output
             & epf_GBout=epf_GB,eef_GBout=eef_GB, evf_GBout=evf_GB, dfe_SIout=dfe_SI,vte_SIout=vte_SI,vre_SIout=vre_SI, vce_SIout=vce_SI,epf_cmout=epf_cm,eef_cmout=eef_cm,evf_cmout=evf_cm, ckeout=cke, & !optional electron flux outputs
             & vene_SIout=vene_SI,chiee_SIout=chiee_SI,vere_SIout=vere_SI,vece_SIout=vece_SI,cekeout=ceke, & !optional electron flux outputs
             & ipf_GBout=ipf_GB,ief_GBout=ief_GB, ivf_GBout=ivf_GB, dfi_SIout=dfi_SI,vti_SIout=vti_SI,vri_SIout=vri_SI, vci_SIout=vci_SI,ipf_cmout=ipf_cm,ief_cmout=ief_cm,ivf_cmout=ivf_cm, ckiout=cki, & !optional ion flux outputs
             & veni_SIout=veni_SI,chiei_SIout=chiei_SI,veri_SIout=veri_SI,veci_SIout=veci_SI,cekiout=ceki, & !optional ion flux outputs
             & modeflagout=modeflag, & ! flags type of modes in output per radial position
             & phiout=phi, npolout=npol, ecoefsout=ecoefs, cftransout=cftrans, &  ! poloidal asymmetry outputs for heavy impurities
             & solfluout=solflu, modewidthout=modewidth, modeshiftout=modeshift, distanout=distan, ntorout=ntor, solout=sol, fdsolout=fdsol,&  !optional 'primitive' outputs from dispersion relation solver needed to build QL flux. Useful for standalone
             & kperp2out=kperp2, krmmuITGout=krmmuITG, krmmuETGout=krmmuETG, &
             & Lcirceout=Lcirce, Lpiegeout=Lpiege, Lecirceout=Lecirce, Lepiegeout=Lepiege, Lvcirceout=Lvcirce, Lvpiegeout=Lvpiege, Lcircgteout=Lcircgte, &
             & Lpieggteout=Lpieggte, Lcircgneout=Lcircgne, Lpieggneout=Lpieggne, Lcircgueout=Lcircgue, Lpieggueout=Lpieggue, Lcircceout=Lcircce, Lpiegceout=Lpiegce, & 
             & Lecircgteout=Lecircgte, Lepieggteout=Lepieggte, Lecircgneout=Lecircgne, Lepieggneout=Lepieggne, Lecircgueout=Lecircgue, Lepieggueout=Lepieggue, Lecircceout=Lecircce, Lepiegceout=Lepiegce, & 
             & Lcirciout=Lcirci, Lpiegiout=Lpiegi, Lecirciout=Lecirci, Lepiegiout=Lepiegi, Lvcirciout=Lvcirci, Lvpiegiout=Lvpiegi, Lcircgtiout=Lcircgti, & 
             & Lpieggtiout=Lpieggti, Lcircgniout=Lcircgni, Lpieggniout=Lpieggni, Lcircguiout=Lcircgui, Lpiegguiout=Lpieggui, Lcircciout=Lcircci, Lpiegciout=Lpiegci, &
             & Lecircgtiout=Lecircgti, Lepieggtiout=Lepieggti, Lecircgniout=Lecircgni, Lepieggniout=Lepieggni, Lecircguiout=Lecircgui, Lepiegguiout=Lepieggui, Lecircciout=Lecircci, Lepiegciout=Lepiegci)

     ENDIF

  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF (myrank == 0) THEN
     CALL outputascii
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

    ! p{6} Flag for including rotation
    rot_flag = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{7} Flag for including rotation
    verbose = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{8} Number of total saught after solutions
    numsols = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{9} 1D integral accuracy
    relacc1 = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{10} 2D integral accuracy
    relacc2 = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{11} Number of runs before runcounter resets
    maxruns = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{12} Maximum number of integrand evaluations in 2D integration routine
    maxpts = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{13} Timeout seconds for a given solution search 
    timeout = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{14} Multiplier for ETG saturation rule (default 1. Mostly for testing)
    ETGmult = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{15} Multiplier for collisionality (default 1. Mostly for testing)
    collmult = readvar(kc,dummy,ktype) ; kc = kc+1

    ! p{16} Toroidal wave-number grid
    ALLOCATE(kthetarhos(dimn))
    kthetarhos = readvar(kc,kthetarhos,ktype) ; kc = kc+1

    ! p{17} Normalised radial coordinate (midplane radius)
    ALLOCATE(x(dimx))
    x = readvar(kc,x,ktype) ; kc = kc+1

    ! p{18} Normalised radial coordinate (midplane radius)
    ALLOCATE(rho(dimx))
    rho = readvar(kc,rho,ktype) ; kc = kc+1

    ! p{19} <Ro> major radius
    ALLOCATE(Ro(dimx))
    Ro = readvar(kc,Ro,ktype) ; kc = kc+1

    ! p{20} <a> minor radius
    ALLOCATE(Rmin(dimx))
    Rmin = readvar(kc,Rmin,ktype) ; kc = kc+1

    ! p{21} B(rho) magnetic field
    ALLOCATE(Bo(dimx))
    Bo = readvar(kc,Bo,ktype) ; kc = kc+1

    ! p{22} R0 geometric major radius (for normalizations)
    R0 = readvar(kc,R0,ktype) ; kc = kc+1

    ! p{23} q(rho) profile
    ALLOCATE(qx(dimx))
    qx = readvar(kc,qx,ktype) ; kc = kc+1

    ! p{24} s(rho) profile
    ALLOCATE(smag(dimx))
    smag = readvar(kc,smag,ktype) ; kc = kc+1

    ! p{25} alpha(rho) profile
    ALLOCATE(alphax(dimx))
    alphax = readvar(kc,alphax,ktype) ; kc = kc+1

    ! p{26} Machtor(rho) profile
    ALLOCATE(Machtor(dimx))
    Machtor = readvar(kc,Machtor,ktype) ; kc = kc+1
!!$    WHERE(ABS(Machtor) < epsD) Machtor = epsD

    ! p{27} Autor(rho) profile
    ALLOCATE(Autor(dimx))
    Autor = readvar(kc,Autor,ktype) ; kc = kc+1
    WHERE(ABS(Autor) < epsD) Autor = epsD

    ! p{28} Machpar(rho) profile
    ALLOCATE(Machpar(dimx))
    Machpar = readvar(kc,Machpar,ktype) ; kc = kc+1
!!$    WHERE(ABS(Machpar) < epsD) Machpar = epsD

    ! p{29} Aupar(rho) profile
    ALLOCATE(Aupar(dimx))
    Aupar = readvar(kc,Aupar,ktype) ; kc = kc+1
    WHERE(ABS(Aupar) < epsD) Aupar = epsD

    ! p{30} gammaE(rho) profile
    ALLOCATE(gammaE(dimx))
    gammaE = readvar(kc,gammaE,ktype) ; kc = kc+1
    WHERE(ABS(gammaE) < epsD) gammaE = epsD

    ! p{31} Te(rho) profile
    ALLOCATE(Tex(dimx))
    Tex = readvar(kc,Tex,ktype) ; kc = kc+1

    ! p{32} ne(rho) profile
    ALLOCATE(Nex(dimx))
    Nex = readvar(kc,Nex,ktype) ; kc = kc+1

    ! p{33} R/LTe(rho) profile
    ALLOCATE(Ate(dimx))
    Ate = readvar(kc,Ate,ktype) ; kc = kc+1
    WHERE(ABS(Ate) < epsD) Ate = epsD

    ! p{34} R/Lne(rho) profile
    ALLOCATE(Ane(dimx))
    Ane = readvar(kc,Ane,ktype) ; kc = kc+1
    WHERE(ABS(Ane) < epsD) Ane = epsD
    WHERE(Ane+Ate < epsD) Ane = Ane+epsD

    ! p{35} Flag for adiabatic electrons
    el_type = INT(readvar(kc,dummy,ktype)) ; kc = kc+1

    ! p{36} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anise(dimx))
    anise = readvar(kc,anise,ktype) ; kc = kc+1

    ! p{37} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisedr(dimx))
    danisedr = readvar(kc,danisedr,ktype) ; kc = kc+1
    WHERE(ABS(danisedr) < epsD) danisedr = epsD

    ! p{38} Main ion mass
    ALLOCATE(Ai(dimx,nions))
    Ai = readvar(kc,Ai,ktype) ; kc = kc+1

    ! p{39} Main ion charge
    ALLOCATE(Zi(dimx,nions))
    Zi = readvar(kc,Zi,ktype) ; kc = kc+1

    ! p{40} Ti(rho) profiles
    ALLOCATE(Tix(dimx,nions))
    Tix = readvar(kc,Tix,ktype) ; kc = kc+1

    ! p{41} ni/ne (rho) profiles
    ALLOCATE(ninorm(dimx,nions))
    ninorm = readvar(kc,ninorm,ktype) ; kc = kc+1

    ! p{42} R/LTi(rho) profiles
    ALLOCATE(Ati(dimx,nions))
    Ati = readvar(kc,Ati,ktype) ; kc = kc+1
    WHERE(ABS(Ati) < epsD) Ati = epsD

    ! p{43} R/Lni(rho) profiles
    ALLOCATE(Ani(dimx,nions))
    Ani = readvar(kc,Ani,ktype) ; kc = kc+1
    WHERE(ABS(Ani) < epsD) Ani = epsD
    WHERE(Ani+Ati < epsD) Ani = Ani + epsD

    ! p{44} Ion types
    ALLOCATE(ion_type(dimx,nions))
    ion_type = INT(readvar(kc,ion_typer,ktype)) ; kc = kc+1
    DEALLOCATE(ion_typer)

    ! p{45} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(anis(dimx,1:nions))
    anis = readvar(kc,anis,ktype) ; kc = kc+1

    ! p{46} Species temp anisotropy at LFS. Zero is electrons
    ALLOCATE(danisdr(dimx,1:nions))
    danisdr = readvar(kc,danisdr,ktype) ; kc = kc+1
    WHERE(ABS(danisdr) < epsD) danisdr = epsD

    ! Read and write runcounter input to decide course of action in calcroutines (full solution or start from previous solution)
    INQUIRE(file="runcounter.dat", EXIST=exist1)
    INQUIRE(file="output/primitive/rsol.dat", EXIST=exist2)
    INQUIRE(file="output/primitive/isol.dat", EXIST=exist3)
    INQUIRE(file="output/primitive/rfdsol.dat", EXIST=exist4)
    INQUIRE(file="output/primitive/ifdsol.dat", EXIST=exist5)

    IF ( (exist1) .AND. (exist2) .AND. (exist3) .AND. (exist4) .AND. (exist5) )THEN
       OPEN(unit=700, file="runcounter.dat", status="old", action="read")
       READ(700,*) runcounter;  CLOSE(700)
    ELSE
       runcounter = 0 !First run
    END IF

    IF (runcounter >= maxruns) THEN !Reset if we're at our maximum number of runs
       runcounter = 0
    ENDIF

    WRITE(fmtn,'(A,I0, A)') '(',dimn,'G15.7)'

    IF (runcounter /= 0) THEN !load old rsol and isol if we're not doing a reset run
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
    ENDIF

    CALL MPI_Barrier(mpi_comm_world,ierror)

    IF (myrank == 0) THEN
       OPEN(unit=700, file="runcounter.dat", status="replace", action="write") !Replace old runcounter with new runcounter
       WRITE(700,*) runcounter + 1 ; CLOSE(700)
    ENDIF

    !DEBUGGING WRITE OUT ALL INPUT TO ASCII FILE

    OPEN(unit=700, file="dimx.dat", action="write", status="replace")
    WRITE(700,'(I0)') dimx ; CLOSE(700)

    OPEN(unit=700, file="dimn.dat", action="write", status="replace")
    WRITE(700,'(I0)') dimn ; CLOSE(700)

    OPEN(unit=700, file="nions.dat", action="write", status="replace")
    WRITE(700,'(I0)') nions ; CLOSE(700)

    WRITE(fmtx, '(A)') '(G15.7)'
    WRITE(fmtxrow,'(A,I2,A)') '(',dimx,'G15.7)'
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

    OPEN(unit=700, file="verbose.dat", action="write", status="replace")
    WRITE(700,'(I0)') verbose ; CLOSE(700)

    OPEN(unit=700, file="numsols.dat", action="write", status="replace")
    WRITE(700,'(I0)') numsols ; CLOSE(700)

    OPEN(unit=700, file="kthetarhos.dat", action="write", status="replace")
    WRITE(700,fmtx) (kthetarhos(j),j=1,dimn) ; CLOSE(700)

    OPEN(unit=700, file="x.dat", action="write", status="replace")
    WRITE(700,fmtx) (x(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="rho.dat", action="write", status="replace")
    WRITE(700,fmtx) (rho(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="Ro.dat", action="write", status="replace")
    WRITE(700,fmtx) (Ro(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="R0.dat", action="write", status="replace")
    WRITE(700,fmtx) R0 ; CLOSE(700)

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
    WRITE(700,'(I8)') maxpts ; CLOSE(700)

    OPEN(unit=700, file="maxruns.dat", action="write", status="replace")
    WRITE(700,'(I3)') maxruns ; CLOSE(700)
    
    OPEN(unit=700, file="relacc1.dat", action="write", status="replace")
    WRITE(700,fmtx) relacc1 ; CLOSE(700)

    OPEN(unit=700, file="relacc2.dat", action="write", status="replace")
    WRITE(700,fmtx) relacc2 ; CLOSE(700)

    OPEN(unit=700, file="timeout.dat", action="write", status="replace")
    WRITE(700,fmtx) timeout ; CLOSE(700)

    OPEN(unit=700, file="ETGmult.dat", action="write", status="replace")
    WRITE(700,fmtx) ETGmult ; CLOSE(700)

    OPEN(unit=700, file="collmult.dat", action="write", status="replace")
    WRITE(700,fmtx) collmult ; CLOSE(700)


    !STOP

    !Allocate output
    ALLOCATE(solflu_SI(dimx,dimn))
    ALLOCATE(solflu_GB(dimx,dimn))
    ALLOCATE(gam_SI(dimx,dimn,numsols))
    ALLOCATE(gam_GB(dimx,dimn,numsols))
    ALLOCATE(ome_SI(dimx,dimn,numsols))
    ALLOCATE(ome_GB(dimx,dimn,numsols))
    ALLOCATE(epf_SI(dimx))  
    ALLOCATE(epfETG_SI(dimx))  
    ALLOCATE(epf_GB(dimx))
    ALLOCATE(eef_SI(dimx))
    ALLOCATE(eefETG_SI(dimx))
    ALLOCATE(eef_GB(dimx))
    ALLOCATE(evf_SI(dimx))
    ALLOCATE(evf_GB(dimx))

    IF (phys_meth /= 0) THEN
       ALLOCATE(dfe_SI(dimx))
       ALLOCATE(vte_SI(dimx))
       ALLOCATE(vce_SI(dimx))
       ALLOCATE(vre_SI(dimx))
       ALLOCATE(cke(dimx))
       IF (phys_meth == 2) THEN
          ALLOCATE(vene_SI(dimx))
          ALLOCATE(chiee_SI(dimx))
          ALLOCATE(vece_SI(dimx))
          ALLOCATE(vere_SI(dimx))
          ALLOCATE(ceke(dimx))
       ENDIF
    ENDIF
    ALLOCATE(epf_cm(dimx,dimn))
    ALLOCATE(eef_cm(dimx,dimn))
    ALLOCATE(evf_cm(dimx,dimn))

    ALLOCATE(ipf_SI(dimx,nions))  
    ALLOCATE(ipf_GB(dimx,nions))
    ALLOCATE(ief_SI(dimx,nions))
    ALLOCATE(ief_GB(dimx,nions))
    ALLOCATE(ivf_SI(dimx,nions))
    ALLOCATE(ivf_GB(dimx,nions))

    IF (phys_meth /= 0) THEN
       ALLOCATE(dfi_SI(dimx,nions))
       ALLOCATE(vti_SI(dimx,nions))
       ALLOCATE(vci_SI(dimx,nions))
       ALLOCATE(vri_SI(dimx,nions))
       ALLOCATE(cki(dimx,nions))
       IF (phys_meth == 2) THEN
          ALLOCATE(veni_SI(dimx,nions))
          ALLOCATE(chiei_SI(dimx,nions))
          ALLOCATE(veci_SI(dimx,nions))
          ALLOCATE(veri_SI(dimx,nions))
          ALLOCATE(ceki(dimx,nions))
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
    ALLOCATE(phi(dimx,ntheta))
    ALLOCATE(npol(dimx,ntheta,nions))
    ALLOCATE(ecoefs(dimx,0:nions,numecoefs)) !includes electrons
    ALLOCATE(cftrans(dimx,nions,7)) ! only for ions

    ALLOCATE( sol (dimx, dimn, numsols) )
    ALLOCATE( fdsol (dimx, dimn, numsols) )
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
    DEALLOCATE(epfETG_SI)
    DEALLOCATE(epf_GB)
    DEALLOCATE(eef_SI)
    DEALLOCATE(eefETG_SI)
    DEALLOCATE(eef_GB)
    DEALLOCATE(evf_SI)
    DEALLOCATE(evf_GB)
    DEALLOCATE(epf_cm)
    DEALLOCATE(eef_cm)
    DEALLOCATE(evf_cm)
    DEALLOCATE(ipf_SI)
    DEALLOCATE(ipf_GB)
    DEALLOCATE(ief_SI)
    DEALLOCATE(ief_GB)
    DEALLOCATE(ivf_SI)
    DEALLOCATE(ivf_GB)

    IF (phys_meth /= 0) THEN
       DEALLOCATE(dfe_SI)
       DEALLOCATE(vte_SI)
       DEALLOCATE(vce_SI)
       DEALLOCATE(vre_SI)
       DEALLOCATE(cke)
       DEALLOCATE(dfi_SI)
       DEALLOCATE(vti_SI)
       DEALLOCATE(vci_SI)
       DEALLOCATE(vri_SI)
       DEALLOCATE(cki)
       IF (phys_meth == 2) THEN
          DEALLOCATE(vene_SI)
          DEALLOCATE(chiee_SI)
          DEALLOCATE(vece_SI)
          DEALLOCATE(vere_SI)
          DEALLOCATE(ceke)
          DEALLOCATE(veni_SI)
          DEALLOCATE(chiei_SI)
          DEALLOCATE(veci_SI)
          DEALLOCATE(veri_SI)
          DEALLOCATE(ceki)
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
    DEALLOCATE(Lvcirce)
    DEALLOCATE(Lvpiege)

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
    WRITE(fmtcftrans,'(A)') '(7G15.7)'

    OPEN(unit=myunit, file="output/primitive/rsolflu.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((REAL(solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/isolflu.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((AIMAG(solflu(i,j)),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=700, file="output/primitive/kymaxITG.dat", action="write", status="replace")
    WRITE(700,fmtx) (krmmuITG(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=700, file="output/primitive/kymaxETG.dat", action="write", status="replace")
    WRITE(700,fmtx) (krmmuETG(j),j=1,dimx) ; CLOSE(700)

    OPEN(unit=myunit, file="output/primitive/distan.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((distan(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/primitive/kperp2.dat", action="write", status="replace")
    WRITE(myunit,fmtn) ((kperp2(i,j),j=1,dimn),i=1,dimx) ; CLOSE(myunit)

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
    WRITE(myunit,fmtecoef) (((ecoefs(j,i,k),k=1,numecoefs),j=1,dimx),i=0,nions) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/cftrans.dat", action="write", status="replace")
    WRITE(myunit,fmtcftrans) (((cftrans(j,i,k),k=1,7),j=1,dimx),i=1,nions) ; CLOSE(myunit)

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

    OPEN(unit=myunit, file="output/epfETG_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (epfETG_SI(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ipf_SI(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/epf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (epf_GB(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/ipf_GB.dat", action="write", status="replace")
    WRITE(myunit,fmtion) ((ipf_GB(i,j),j=1,nions),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/eef_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (eef_SI(i),i=1,dimx) ; CLOSE(myunit)

    OPEN(unit=myunit, file="output/eefETG_SI.dat", action="write", status="replace")
    WRITE(myunit,fmtx) (eefETG_SI(i),i=1,dimx) ; CLOSE(myunit)

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

END PROGRAM qlk_standalone
