MODULE datmat
  USE kind
  ! Argument declaration 

  ! List of input variables
  INTEGER, SAVE :: dimx, dimn, nions, numsols, phys_meth, coll_flag, rot_flag, el_type
  INTEGER, SAVE, DIMENSION(:,:), ALLOCATABLE :: ion_type
  REAL(KIND=DBL), SAVE  :: R0
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: x, rho, Ro, Rmin, Bo, kthetarhos, qx, smag, alphax, rotflagarray
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Tex, Nex, Ate, Ane, anise, danisedr, filterprof
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: Tix, Ati, Ani, ninorm, anis, danisdr, Ai, Zi
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Machtor, Autor, Machpar, Aupar, gammaE
  LOGICAL, SAVE :: verbose,separateflux   !Boolean input parameters

  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Machparmod, Auparmod, gammaEmod !backup arrays used in rotationless dispersion relation solver with rot_flag=2
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Machparorig, Auparorig, gammaEorig !backup arrays used in rotationless dispersion relation solver with rot_flag=2

  !integral testing parameters
  INTEGER, SAVE :: int_method, newt_method, newt_conv, int_split, norm
  REAL(KIND=DBL), SAVE :: reqrelacc, reqabsacc, reqrelacc_newt, reqabsacc_newt

  !Parameters for deciding how often to jump to full solution searching in integrated modelling applications
  INTEGER, SAVE :: maxruns !default is 50
  INTEGER, SAVE :: maxpts !Max number of integrand evaluations in 2D integrals. Default = 1.d5
  INTEGER, SAVE :: lenwrk !Work array size for 2D integral routine

  REAL(KIND=DBL), SAVE :: relacc1 !  !1D integral relative error demanded. Default = 1.0d-3
  REAL(KIND=DBL), SAVE :: relacc2 !2D integral relative error demanded. Default = 2.0d-2
  REAL(KIND=DBL), SAVE :: relaccQL1 !  !1D integral relative error demanded. Default = 1.0d-3
  REAL(KIND=DBL), SAVE :: relaccQL2 !2D integral relative error demanded. Default = 2.0d-2

  REAL(KIND=DBL), SAVE :: ETGmult, collmult !multpliers for ETG saturation level and collisionality (default 1)

  ! List of derived 'radial' variables (see mod_make_input.f90 for details)
  REAL(KIND=DBL), SAVE :: widthtuneITG, widthtuneETG
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: csou, cref,cthe, de, epsilon, qprim, ft, fc
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Ac, Rhoe, Anue, Zeffx, omegator, domegatordr, Nustar
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: wg, rhostar, Rhoeff, ktetasn, Mache, Aue, krmmuITG,krmmuETG
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: coefi, Machi, Machitemp, Aui, Nix, Rhoi, di, cthi, tau, mi, Auiorig, Auimod, Machiorig, Machimod

  ! Output arrays. The 3 dimensions are 'radial grid', 'kthetarhos grid', 'number of modes'
  REAL(KIND=DBL), SAVE , DIMENSION(:,:), ALLOCATABLE :: distan,FLRec,FLRep,ntor,kperp2
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: gamma, Ladia, FLRip, FLRic
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: modewidth, modeshift, modeshift2
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: ommax, solflu
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: sol, fdsol, oldsol, oldfdsol
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: Lcirce, Lpiege, Lecirce, Lepiege, Lcircgte, Lpieggte,  Lcircgne, Lpieggne, Lcircce, Lpiegce
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: Lcirci, Lpiegi, Lecirci, Lepiegi, Lvcirci, Lvpiegi, Lcircgti, Lpieggti, Lcircgni, Lpieggni, Lcircgui, Lpieggui, Lcircci, Lpiegci
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: Lecircgte, Lepieggte,  Lecircgne, Lepieggne, Lecircce, Lepiegce
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: Lecircgti, Lepieggti, Lecircgni, Lepieggni, Lecircgui, Lepieggui, Lecircci, Lepiegci

  ! Final output arrays following saturation rule. These can be printed as ASCII output
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: solflu_SI, solflu_GB
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: gam_SI,gam_GB,ome_SI,ome_GB
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: epf_SI,epf_GB,eef_SI,eef_GB,eefETG_SI,eefETG_GB
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: dfe_SI,vte_SI,vce_SI,cke,modeflag
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: dfe_GB,vte_GB,vce_GB
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: ipf_SI,ipf_GB,ief_SI,ief_GB,ivf_SI,ivf_GB, dfi_SI,vti_SI,vci_SI,vri_SI,cki,eef_cm,epf_cm
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: dfi_GB,vti_GB,vci_GB,vri_GB

  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: vene_SI,chiee_SI,vece_SI,vene_GB,chiee_GB,vece_GB,ceke
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: veneETG_SI,chieeETG_SI,veceETG_SI,veneETG_GB,chieeETG_GB,veceETG_GB
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: veni_SI,chiei_SI,veci_SI,veri_SI,veni_GB,chiei_GB,veci_GB,veri_GB,ceki

  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: ipf_cm,ief_cm,ivf_cm

  ! Poloidal asymmetry variables
  ! Ion density along field line, used was asymmetries are present. 
  ! 1st dimension is radius. 2nd is theta (ntheta points from [0,pi]), 3rd is ion
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: npol,tpernormi,Anipol,ecoefs
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: ecoefsgau
  REAL(KIND=DBL), SAVE, DIMENSION(:,:,:), ALLOCATABLE :: cftrans !output variable containing the transport coefficients calculated including the 2D properties
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: th !used to pass around theta gridq
  REAL(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: phi,dphidr,tpernorme,Anepol, nepol,dphidth
  INTEGER, SAVE :: itheta !used to pass around theta coordinate
  INTEGER, SAVE :: irad !used to pass around radial coordinate
  INTEGER, SAVE :: inu !used to pass around wavenumber coordinate
  REAL(KIND=DBL), SAVE :: thetapass !used to pass around angle
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: tpernormifunc !used to pass around arbitrary angle tpernormi
  REAL(KIND=DBL), SAVE :: tpernormefunc !used to pass around arbitrary angle tpernorme

  ! Varibles used by calcroutines
  INTEGER, SAVE :: ion ! current ion index used in integrals and asymmetry functions
  INTEGER, SAVE :: runcounter ! used for counting runs inside integrated modelling applications for deciding to recalculate all or just jump to newton based on old solutions
  REAL(KIND=DBL), SAVE :: Joe2, Jobane2, Joe2p, J1e2p, Joe2c
  REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: Joi2, Jobani2, Joi2p, J1i2p, Joi2c
  REAL(KIND=DBL), SAVE :: ktetaRhoe
  REAL(KIND=DBL), SAVE :: d, normkr
  REAL(KIND=DBL), SAVE :: Athe
  REAL(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: Athi, ktetaRhoi
  REAL(KIND=DBL), SAVE :: nwg, qRd, omega2bar
  REAL(KIND=DBL), SAVE :: fonxad
  COMPLEX(KIND=DBL), SAVE :: mwidth, mshift, mshift2, omeflu ! the width and the shift of eigenfunctions can be complex. omeflu is the solution used for the contour limit from Pierre solution
  COMPLEX(KIND=DBL), SAVE :: mwidth_rot, mshift_rot ! width and the shift of eigenfunctions for momentum transport cases
  REAL(KIND=DBL), SAVE :: widthhat ! the width corresponding to |phi(x)| . Used instead of mwidth when phi(x) is complex
  COMPLEX(KIND=DBL), SAVE :: fonxcirce, fonxpiege, fonxcircgte, fonxpieggte, fonxcircgne
  COMPLEX(KIND=DBL), SAVE :: fonxpieggne, fonxcircce, fonxpiegce, fonxecirce, fonxepiege 
  COMPLEX(KIND=DBL), SAVE :: fonxecircgte, fonxepieggte, fonxecircgne, fonxepieggne, fonxecircce, fonxepiegce

  COMPLEX(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: fonxcirci, fonxpiegi, fonxcircgti, fonxpieggti, fonxcircgni, fonxcircgui
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: fonxpieggni, fonxpieggui, fonxcircci, fonxpiegci, fonxecirci, fonxepiegi, fonxvcirci, fonxvpiegi
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: fonxecircgti, fonxepieggti, fonxecircgni, fonxecircgui, fonxepieggni, fonxepieggui, fonxecircci, fonxepiegci


  LOGICAL, SAVE, DIMENSION(:), ALLOCATABLE :: ETG_flag
  REAL(KIND=DBL), SAVE :: rint !real radius in contours 

  INTEGER, SAVE :: pnFLR !radial coordinate used in FLR routines
  !INTEGER, SAVE :: caseflag
  COMPLEX(kind=DBL), SAVE :: omegmax

  !Variables used in integration routines
  COMPLEX(KIND=DBL), SAVE :: omFFk !freq used in trapped particle integrals
  INTEGER, SAVE :: pFFk !radial coordinate used in trapped particle integrals
  COMPLEX(KIND=DBL), SAVE :: omFkr !freq used in passing particle integrals
  INTEGER, SAVE :: pFkr !radial coordinate used in passing particle integrals
  INTEGER, SAVE :: nuFkr !wavenumber coordinate used in passing particle integrals
  INTEGER, SAVE :: nuFFk !wavenumber coordinate used in trapped particle integrals
  REAL(KIND=DBL) :: sin2th,alamnorm,alam1,alam2,alam3,alam4,alam5 !pitch angle averages of vpar^m used in passints
  INTEGER, SAVE :: plam,nulam !radial and wavenumber coordinates used to pass around in Vpar averaging routines in passints

  INTEGER, SAVE :: weidcount !count the dispersion function calls
  INTEGER, SAVE :: ccount !count the integrand calls
  INTEGER, SAVE :: Nsolrat, last
  REAL(kind=DBL), SAVE, DIMENSION(:), ALLOCATABLE :: alist, blist, rlist, elist
  INTEGER, SAVE, DIMENSION(:), ALLOCATABLE :: iord

  !Variables used in fluid solution with rotation
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: jon_modewidth,jon_solflu,jon_modeshift
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: cot_modewidth,cot_solflu,cot_modeshift
  COMPLEX(KIND=DBL), SAVE, DIMENSION(:,:), ALLOCATABLE :: old_modewidth,ana_solflu,old_modeshift

  !Timing variable
  REAL(KIND=DBL), SAVE :: calltimeinit,timeout
  LOGICAL, SAVE :: timeoutflag

  !min and max radius for calculation
  REAL(KIND=DBL), SAVE :: rhomin,rhomax

  !EXTERNAL FUNCTION AND SUBROUTINE DECLARATIONS

  !Elliptic integrals (from SLATEC)
  REAL(KIND=DBL) :: ceik, ceie, DGAMMA2
  EXTERNAL ceik, ceie, DGAMMA2

  !Scaled modified Bessel functions routines (from SPECFUN)
  REAL(KIND=DBL) :: BESEI0, BESEI1
  EXTERNAL BESEI0, BESEI1

  !NAG proxies
  REAL(KIND=DBL) :: d01ahf
  EXTERNAL d01ahf, d01fcf

  !Used for DFZERO
  !  EXTERNAL phieq

END MODULE datmat
