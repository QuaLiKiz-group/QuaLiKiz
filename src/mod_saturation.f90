!NOTE: In this saturation rule, the wavenumbers are integrated over when including a zero flux point at ky=0. Thus, if code-to-code comparisons of exact quasilinear fluxes for a single wavenumber is desired,
!!!!!  then the code will have to be (easily) modified

MODULE mod_saturation
  USE kind
  USE datmat
  USE datcal
  USE mod_make_io ! for MPI include

  IMPLICIT NONE

CONTAINS
  SUBROUTINE allocate_endoutput()

    ALLOCATE(solflu_SI(dimx,dimn)); solflu_SI=0
    ALLOCATE(solflu_GB(dimx,dimn)); solflu_GB=0
    ALLOCATE(gam_SI(dimx,dimn,numsols)); gam_SI=0
    ALLOCATE(gam_GB(dimx,dimn,numsols)); gam_GB=0
    ALLOCATE(ome_SI(dimx,dimn,numsols)); ome_SI=0
    ALLOCATE(ome_GB(dimx,dimn,numsols)); ome_GB=0

    ALLOCATE(modeflag(dimx)); modeflag=0
    ALLOCATE(epf_SI(dimx)); epf_SI=0
    ALLOCATE(epf_GB(dimx)); epf_GB=0
    ALLOCATE(eef_SI(dimx)); eef_SI=0
    ALLOCATE(eefETG_SI(dimx)); eefETG_SI=0
    ALLOCATE(eefETG_GB(dimx)); eefETG_GB=0
    ALLOCATE(eef_GB(dimx)); eef_GB=0
    ALLOCATE(ipf_SI(dimx,nions)); ipf_SI=0
    ALLOCATE(ipf_GB(dimx,nions)); ipf_GB=0
    ALLOCATE(ief_SI(dimx,nions)); ief_SI=0
    ALLOCATE(ief_GB(dimx,nions)); ief_GB=0
    ALLOCATE(ivf_SI(dimx,nions)); ivf_SI=0
    ALLOCATE(ivf_GB(dimx,nions)); ivf_GB=0

    IF (phys_meth /= 0) THEN
       ALLOCATE(dfe_SI(dimx)); dfe_SI=0
       ALLOCATE(vte_SI(dimx)); vte_SI=0
       ALLOCATE(vce_SI(dimx)); vce_SI=0

       ALLOCATE(dfe_GB(dimx)); dfe_GB=0
       ALLOCATE(vte_GB(dimx)); vte_GB=0
       ALLOCATE(vce_GB(dimx)); vce_GB=0

       ALLOCATE(cke(dimx)); cke=0
       ALLOCATE(dfi_SI(dimx,nions)); dfi_SI=0
       ALLOCATE(vti_SI(dimx,nions)); vti_SI=0
       ALLOCATE(vci_SI(dimx,nions)); vci_SI=0
       ALLOCATE(vri_SI(dimx,nions)); vri_SI=0

       ALLOCATE(dfi_GB(dimx,nions)); dfi_GB=0
       ALLOCATE(vti_GB(dimx,nions)); vti_GB=0
       ALLOCATE(vci_GB(dimx,nions)); vci_GB=0
       ALLOCATE(vri_GB(dimx,nions)); vri_GB=0

       ALLOCATE(cki(dimx,nions)); cki=0
!!!
       IF (phys_meth == 2) THEN
          ALLOCATE(vene_SI(dimx)); vene_SI=0
          ALLOCATE(chiee_SI(dimx)); chiee_SI=0
          ALLOCATE(vece_SI(dimx)); vece_SI=0
          ALLOCATE(ceke(dimx)); ceke=0
          ALLOCATE(veni_SI(dimx,nions)); veni_SI=0
          ALLOCATE(chiei_SI(dimx,nions)); chiei_SI=0
          ALLOCATE(veci_SI(dimx,nions)); veci_SI=0
          ALLOCATE(veri_SI(dimx,nions)); veri_SI=0
          ALLOCATE(ceki(dimx,nions)); ceki=0
          ALLOCATE(vene_GB(dimx)); vene_GB=0
          ALLOCATE(chiee_GB(dimx)); chiee_GB=0
          ALLOCATE(vece_GB(dimx)); vece_GB=0
          ALLOCATE(veni_GB(dimx,nions)); veni_GB=0
          ALLOCATE(chiei_GB(dimx,nions)); chiei_GB=0
          ALLOCATE(veci_GB(dimx,nions)); veci_GB=0
          ALLOCATE(veri_GB(dimx,nions)); veri_GB=0
          IF (separateflux .EQV. .TRUE.) THEN
             ALLOCATE(veneETG_SI(dimx)); veneETG_SI=0
             ALLOCATE(chieeETG_SI(dimx)); chieeETG_SI=0
             ALLOCATE(veceETG_SI(dimx)); veceETG_SI=0
             ALLOCATE(veneETG_GB(dimx)); veneETG_GB=0
             ALLOCATE(chieeETG_GB(dimx)); chieeETG_GB=0
             ALLOCATE(veceETG_GB(dimx)); veceETG_GB=0
          ENDIF
       ENDIF
    ENDIF

    ALLOCATE(epf_cm(dimx,dimn)); epf_cm=0
    ALLOCATE(eef_cm(dimx,dimn)); eef_cm=0
    ALLOCATE(ipf_cm(dimx,dimn,nions)); ipf_cm=0
    ALLOCATE(ief_cm(dimx,dimn,nions)); ief_cm=0
    ALLOCATE(ivf_cm(dimx,dimn,nions)); ivf_cm=0

  END SUBROUTINE allocate_endoutput

  SUBROUTINE deallocate_endoutput()

    DEALLOCATE(solflu_SI)
    DEALLOCATE(solflu_GB)
    DEALLOCATE(gam_SI)
    DEALLOCATE(gam_GB)
    DEALLOCATE(ome_SI)
    DEALLOCATE(ome_GB)

    DEALLOCATE(modeflag)
    DEALLOCATE(epf_SI)
    DEALLOCATE(epf_GB)
    DEALLOCATE(eef_SI)
    DEALLOCATE(eefETG_SI)
    DEALLOCATE(eefETG_GB)

    DEALLOCATE(eef_GB)
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
          DEALLOCATE(vene_GB)
          DEALLOCATE(chiee_GB)
          DEALLOCATE(vece_GB)
          DEALLOCATE(veni_GB)
          DEALLOCATE(chiei_GB)
          DEALLOCATE(veci_GB)
          DEALLOCATE(veri_GB)
          IF (separateflux .EQV. .TRUE.) THEN
             DEALLOCATE(veneETG_SI)
             DEALLOCATE(chieeETG_SI)
             DEALLOCATE(veceETG_SI)
             DEALLOCATE(veneETG_GB)
             DEALLOCATE(chieeETG_GB)
             DEALLOCATE(veceETG_GB)
          ENDIF
       ENDIF
    ENDIF

    DEALLOCATE(epf_cm)
    DEALLOCATE(eef_cm)
    DEALLOCATE(ipf_cm)
    DEALLOCATE(ief_cm)
    DEALLOCATE(ivf_cm)


  END SUBROUTINE deallocate_endoutput

  SUBROUTINE saturation(outputcase)
    INTEGER, INTENT(IN) :: outputcase !0 for all modes, 1 for ITG only, 2 for TEM only. ETG-only is calculated anyway for outputcase==0
    INTEGER :: ir,j,k,gg,ifailloc
    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: kteta,kthr,kxshift,nwgmat,smagn,qxn,dw
    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: kx2shear,kxadd,kxnl
    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: maxgmsp,maxgmsptmp
    REAL(KIND=DBL), DIMENSION(dimx) :: gamGB, mdmlmu, mdml, krm, chi_GB,nathanfac
    REAL(KIND=DBL), DIMENSION(dimn+1) :: xint,yint
    INTEGER, DIMENSION(dimx) :: inddmlmu,inddml
    INTEGER, DIMENSION(1) :: maxloci !damn you MAXLOC and your silly dimension 1 array requirement
    REAL(KIND=DBL) :: cfaca,cfacb,cfacc,cfacd,qfac,rhos,locmaxgamma
    REAL(KIND=DBL) :: sfac, normNL !normalization factors
    REAL(KIND=DBL),DIMENSION(dimn) :: normETG
    REAL(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: fi, constp, conste, constv, cmpfe_k, cmefe_k, cmpfgne_k, cmpfgte_k, cmpfce_k, cmefgne_k, cmefgte_k, cmefce_k
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions,numsols) :: cmpfi_k, cmefi_k, cmvfi_k, cmpfgni_k, cmpfgti_k, cmpfgui_k, cmpfci_k, cmefgni_k, cmefgti_k, cmefgui_k, cmefci_k
    COMPLEX(KIND=DBL), DIMENSION(dimx,dimn,numsols) :: solbck,solbcktmp
    REAL(KIND=DBL), DIMENSION(dimx,dimn) :: cmpfe, cmefe, cmpfgne, cmpfgte, cmpfce, cmefgne, cmefgte, cmefce
    REAL(KIND=DBL), DIMENSION(dimx,dimn,nions) :: cmpfi, cmefi, cmvfi, cmpfgni, cmpfgti, cmpfgui, cmpfci, cmefgni, cmefgti, cmefgui, cmefci
    REAL(KIND=DBL), DIMENSION(dimx) :: pfe, dpfe,efe, efeETG, defe,defeETG, dffte, vthte, vcpte, deffte, vethte, vecpte, deffteETG, vethteETG, vecpteETG, ion_epf_GB, ion_eef_GB, ele_epf_GB, ele_eef_GB
    REAL(KIND=DBL), DIMENSION(dimx,nions) :: pfi, dpfi,efi, defi,vfi, dvfi, dffti, vthti, vcpti, vrdti, deffti, vethti, vecpti, verdti, ion_ipf_GB, ion_ief_GB, ion_ivf_GB, ele_ipf_GB, ele_ief_GB, ele_ivf_GB
    REAL(KIND=DBL) :: alphp,alphm,lowlim,massrat,rat
    CHARACTER(len=7) :: fmtx,fmtn,fmtion !for debugging
    INTEGER :: kk,i,myunit=700,ETGind
    !MPI variables:
    INTEGER :: doit,ierror,myrank,nproc

    CALL mpi_comm_size(mpi_comm_world,nproc,ierror)
    CALL mpi_comm_rank(mpi_comm_world,myrank,ierror)

    Machi=Machitemp ! Reinstate impurity Mach numbers. Ordering is valid here (for cftrans asymmetry terms)
    IF (rot_flag == 2) Aui=Auimod !modify gradients for rotodiffusion

    !The all important normalization factor tuned from GASTD ion heat transport nonlinear simulation
    normNL=271.

    WRITE(fmtn,"(I0, A)") dimn,'E16.7' !for debugging printouts (yes yes, we don't use a real debugging tool, don't get sanctimonious)

    !Intialize fi
    fi(:,:,:) = 0.

    modeflag(:)=0 !set default that both modes are important (also true if stable, since we don't know where they could pop up)

    !Carry out the nonlinear saturation rules!

    !Define auxilliary variables
    gamGB(:)=csou(:)/Rmin(:) !gyrobohm 1/s unit

    !Initialising of auxilliary varibles
    DO ir=1,dimx
       DO j = 1,dimn
          smagn(ir,j)=smag(ir) !2D versions of smag and qx useful for calculations
          qxn(ir,j)=qx(ir)
          kteta(ir,j)  = ntor(ir,j)*qx(ir)/(Rmin(ir)*x(ir))  
          kthr(ir,j)   = ntor(ir,j)*qx(ir)*rhostar(ir)/x(ir)
          nwgmat(ir,j) = ntor(ir,j)*wg(ir)
          solflu_SI(ir,j) = solflu(ir,j)*nwgmat(ir,j)
          solflu_GB(ir,j) = solflu_SI(ir,j)/gamGB(ir)
          DO k=1,numsols
             gam_SI(ir,j,k)=AIMAG(sol(ir,j,k))*nwgmat(ir,j)
             gam_GB(ir,j,k)=AIMAG(sol(ir,j,k))*nwgmat(ir,j)/gamGB(ir)
             ome_SI(ir,j,k)=REAL(sol(ir,j,k))*nwgmat(ir,j)
             ome_GB(ir,j,k)=REAL(sol(ir,j,k))*nwgmat(ir,j)/gamGB(ir)
          ENDDO
       ENDDO
    ENDDO

    !Find index of first ETG-scale mode (if it exists)
    ETGind=0
    DO j=1,dimn
       IF (kthetarhos(j) > ETGk) THEN         
          ETGind = j
          EXIT
       ENDIF
    ENDDO

    massrat=SQRT(mi(1,1)/me) ! used for rudimentary multi-scale ETG model based on Nathan Howard NF 2016 'rule of thumb'
    nathanfac(:)=1. !initialization
    IF (ETGind>1) THEN
       DO ir = 1,dimx 
          IF (MAXVAL(gam_GB(ir,1:j-1,:))>epsD) THEN
             rat = MAXVAL(gam_GB(ir,j:dimn,:))/MAXVAL(gam_GB(ir,1:j-1,:))/massrat
             IF (rat > 3.) rat=3. !to avoid overflow in exponent below
             nathanfac(ir)= 1./(1.+EXP(-(rat*massrat-massrat)/5.))
          ENDIF
       ENDDO
    ENDIF

    ! 	QUASILINEAR FLUX CALCULATIONS
    DO gg=1,3
       solbck=sol

       ! If outputcase=0 (all modes included) then we go through 3 iterations of calculation in order to 
       ! test whether ion or electron modes are negligible. 
       ! An output array for each radius signifying 'all ion modes', 'all electron modes', or 'stable' is returned
       IF (outputcase == 0) THEN
          IF (gg == 1) THEN
             WHERE (REAL(solbck) > 0.) solbck=0.  !kill all electron modes
          ENDIF
          IF (gg == 2) THEN
             WHERE (REAL(solbck) < 0.) solbck=0.  !kill all ion modes
          END IF
       ELSE ! If outputcase=1, then only keep ITG modes, if =2, then only TEM. Separate ETG is calculated anyway with outputcase=0. For ITG or TEM only, no need for gg for loop
          IF (outputcase == 1) THEN
             WHERE (REAL(solbck) > 0.) solbck=0.  !kill all electron modes
             IF (gg<3) CYCLE
          ELSEIF (outputcase == 2) THEN !kill all ion modes and electron-scale modes
             WHERE (REAL(solbck) < 0.) solbck=0.  
             DO ir=1,dimx
                DO k=1,numsols
                   WHERE (kthetarhos > ETGk) solbck(ir,:,k) = 0
                ENDDO
             ENDDO
             IF (gg<3) CYCLE
          ENDIF
       ENDIF

       ! dw=distan**2/(ABS(modewidth)**2 / SQRT(REAL(modewidth**2)))**2*DGAMMA2(0.75)/DGAMMA2(0.25) 

       kxshift = (distan*AIMAG(modeshift2)/REAL(modewidth**2))**2 !kx contribution from modeshift
       dw=0.5*distan**2/(REAL(modewidth**2)) + kxshift

       alphp   = -3.0 !used for spectrum shape above kymax
       alphm   = 1.0 !used for spectrum shape below kymax

       doit = 0 !initialize counter for distributing tasks

       DO ir = 1,dimx !big cycle on scan (or radial) parameter

          !additional normalization factor for ETG transport
          normETG(:)=1.0
          WHERE(kthetarhos > ETGk) normETG = 0.45*ETGmult*nathanfac(ir) !0.15 if we believe multiscale!

          chi_GB(ir)=SQRT(Ai(ir,1)*mp)/(qe**2*Bo(ir)**2)*((Tex(ir)*1e3*qe)**1.5)/Rmin(ir)  !GyroBohm normalisation in m^2/s based on main ion

          IF (myrank == doit) THEN ! distribute independent loop indices to tasks
             doit=doit+1; 
             IF (doit==nproc) doit=0 
          ELSE
             doit=doit+1; 
             IF (doit==nproc) doit=0 
             CYCLE
          ENDIF

          !CALCULATE NEW NON-LINEAR CONTRIBUTION TO Kx (JC 12.2011)
          rhos=SQRT(Tex(ir)*1d3*qe*mi(ir,1))/(qe*Bo(ir)) !Larmor radius with respect to sound speed (no sqrt(2))
          cfaca=0.4; cfacb=2. ; cfacc=1.5 ; cfacd=0.2 ; qfac=0.5 !tuned factors       

          kx2shear(ir,:) = kteta(ir,:)**2*(smagn(ir,:)**2)*2*dw(ir,:) !contribution of kx from magnetic shear

          kxadd(ir,:)=kteta(ir,:)*rhos
          kxadd(ir,:)=(kxadd(ir,:)-cfacd)*cfacc
          WHERE (kxadd<0) kxadd=0. !'isotropic part' of kx at higher ky

          kxnl(ir,:)=(cfaca*(EXP(-cfacb*ABS(smagn(ir,:))) )*(1./qxn(ir,:)**qfac)+kxadd(ir,:))/rhos !nonlinear contribution to kx

          kperp2(ir,:)=kteta(ir,:)**2 + ((kx2shear(ir,:)**0.5)+kxnl(ir,:))**2 !construct finally the new kperp2 
          WHERE(kthetarhos>ETGk) kperp2(ir,:)=2.*kteta(ir,:)**2 !ETG kperp2 based on streamers isotropisation


          !CALCULATE NEW NORMALIZATION FACTOR BASED ON FREQUENCY BROADENING / DAMPED LINEAR MODES
          IF (ABS(smag(ir))<0.6) THEN
             sfac=-2.5*ABS(smag(ir))+2.5
          ELSE
             sfac=1._DBL
          END IF

          ! k-spectrum of max(gamma) [s^-1] for each radial position

          DO j=1,dimn
             maxgmsp(ir,j)=MAXVAL( imag(solbck(ir,j,:) )*nwgmat(ir,j)) 
          END DO

          !here is calculated the ktheta max for the NL spectrum only on the most unstable mode if more than 1 is present
          !mu stands for "most unstable"

          IF ( (kthr(ir,1) > ETGk) .OR. (kthr(ir,dimn) <= ETGk) .OR. (dimn ==1) ) THEN  !single mode, pure ETG-scale, or pure ITG-scale
             mdmlmu(ir) = MAXVAL(maxgmsp(ir,:)/kperp2(ir,:))
             maxloci=MAXLOC(maxgmsp(ir,:)/kperp2(ir,:))
             inddmlmu(ir)=maxloci(1)

             IF ( kthetarhos(inddmlmu(ir)) <= 0.05) THEN 
                inddmlmu(ir)=inddmlmu(ir)+1 
                ! Not a rigorous fix. Used to avoid some unphysical cases was first instability in k-spectrum gives max gam/kperp^2
             ENDIF
             mdmlmu(ir)=mdmlmu(ir)/kteta(ir,inddmlmu(ir)) !Used to fit in quasilinear flux integral
             !          krmmu(ir) = kthr(ir,inddmlmu(ir))
             IF (kthr(ir,1) <= ETGk) THEN !ITG scales
                krmmuITG(ir) = kthetarhos(inddmlmu(ir))
                krmmuETG(ir) = 0.
             ELSE !electron scales
                krmmuETG(ir) = kthetarhos(inddmlmu(ir))
                krmmuITG(ir) = 0.
             ENDIF

             !Saturation rules for each unstable root
             DO k=1,numsols
                !Some of the above is actually repeated here. Have to look deeper to see if code can be slightly reduced
                mdml(ir) = MAXVAL(AIMAG(solbck(ir,:,k))*nwgmat(ir,:)/kperp2(ir,:))
                maxloci = MAXLOC(AIMAG(solbck(ir,:,k))*nwgmat(ir,:)/kperp2(ir,:))
                inddml=maxloci(1)
                IF ( kthetarhos(inddml(ir)) <= 0.05) THEN 
                   inddml(ir)=inddml(ir)+1 
                   ! Not a rigorous fix. Used to avoid some unphysical cases was first instability in k-spectrum gives max gam/kperp^2
                ENDIF

                IF (SIZE(kthetarhos) == 1) THEN
                   inddml(ir)=1 !Avoids error if only 1 k is analyzed (typically in a standalone version)
                ENDIF
                mdml(ir)=mdml(ir)/kteta(ir,inddml(ir))
                krm(ir) = kthr(ir,inddml(ir))
                fi(ir,1:inddmlmu(ir),k) = mdml(ir) * kthr(ir,1:inddmlmu(ir)) **(alphm) / kthr(ir,inddmlmu(ir))**(alphm) / R0
                fi(ir,inddmlmu(ir)+1:dimn,k) = mdml(ir) * kthr(ir,inddmlmu(ir)+1:dimn)**(alphp) / kthr(ir,inddmlmu(ir))**(alphp) / R0

                fi(ir,:,k)=fi(ir,:,k)/sfac*normETG !renormalize fi 
             ENDDO
          ELSE !separate the scales for fi calculation

             DO kk=1,2 ! 1 for ion scales, 2 for electron scales
                maxgmsptmp = 0.
                solbcktmp = 0.
                IF (kk == 1) THEN !ion scales
                   maxgmsptmp(ir,1:ETGind-1)=maxgmsp(ir,1:ETGind-1)
                   solbcktmp(ir,1:ETGind-1,:)=solbck(ir,1:ETGind-1,:)
                ELSE !electron scales
                   maxgmsptmp(ir,ETGind:dimn)=maxgmsp(ir,ETGind:dimn)
                   solbcktmp(ir,ETGind:dimn,:)=solbck(ir,ETGind:dimn,:)
                ENDIF
                mdmlmu(ir) = MAXVAL(maxgmsptmp(ir,:)/kperp2(ir,:))

                IF (ABS(mdmlmu(ir)) < epsD) CYCLE !The current scale is stable, so leave fi in that scale as 0 and cycle 

                maxloci=MAXLOC(maxgmsptmp(ir,:)/kperp2(ir,:))
                inddmlmu(ir)=maxloci(1)


                IF ( kthetarhos(inddmlmu(ir)) <= 0.05) THEN 
                   inddmlmu(ir)=inddmlmu(ir)+1 
                   ! Not a rigorous fix. Used to avoid some unphysical cases was first instability in k-spectrum gives max gam/kperp^2
                ENDIF
                mdmlmu(ir)=mdmlmu(ir)/kteta(ir,inddmlmu(ir)) !Used to fit in quasilinear flux integral
                !          krmmu(ir) = kthr(ir,inddmlmu(ir))
                IF (kk==1) THEN !ion scales
                   krmmuITG(ir) = kthetarhos(inddmlmu(ir))
                ELSE !electron scales
                   krmmuETG(ir) = kthetarhos(inddmlmu(ir))
                ENDIF

                !Saturation rules for each unstable root
                DO k=1,numsols
                   !Some of the above is actually repeated here. Have to look deeper to see if code can be slightly reduced
                   mdml(ir) = MAXVAL(AIMAG(solbcktmp(ir,:,k))*nwgmat(ir,:)/kperp2(ir,:))
                   maxloci = MAXLOC(AIMAG(solbcktmp(ir,:,k))*nwgmat(ir,:)/kperp2(ir,:))
                   inddml=maxloci(1)
                   IF ( kthetarhos(inddml(ir)) <= 0.05) THEN 
                      inddml(ir)=inddml(ir)+1 
                      ! Not a rigorous fix. Used to avoid some unphysical cases was first instability in k-spectrum gives max gam/kperp^2
                   ENDIF

                   IF (SIZE(kthetarhos) == 1) THEN
                      inddml(ir)=1 !Avoids error if only 1 k is analyzed (typically in a standalone version)
                   ENDIF
                   mdml(ir)=mdml(ir)/kteta(ir,inddml(ir))
                   krm(ir) = kthr(ir,inddml(ir))

                   IF (kk == 1) THEN !ion scales
                      fi(ir,1:inddmlmu(ir),k) = mdml(ir) * kthr(ir,1:inddmlmu(ir)) **(alphm) / kthr(ir,inddmlmu(ir))**(alphm) / R0
                      fi(ir,inddmlmu(ir)+1:ETGind-1,k) = mdml(ir) * kthr(ir,inddmlmu(ir)+1:ETGind-1)**(alphp) / kthr(ir,inddmlmu(ir))**(alphp) / R0
                      fi(ir,1:ETGind-1,k)=fi(ir,1:ETGind-1,k)/sfac*normETG(1:ETGind-1) !renormalize fi 

                   ELSE !electron scales
                      fi(ir,ETGind:inddmlmu(ir),k) = mdml(ir) * kthr(ir,ETGind:inddmlmu(ir)) **(alphm) / kthr(ir,inddmlmu(ir))**(alphm) / R0
                      fi(ir,inddmlmu(ir)+1:dimn,k) = mdml(ir) * kthr(ir,inddmlmu(ir)+1:dimn)**(alphp) / kthr(ir,inddmlmu(ir))**(alphp) / R0
                      fi(ir,ETGind:dimn,k)=fi(ir,ETGind:dimn,k)/sfac*normETG(ETGind:dimn) !renormalize fi 
                   ENDIF
                ENDDO
             ENDDO

          ENDIF


          DO j=1,dimn

             DO k=1,numsols

                locmaxgamma=MAXVAL( AIMAG(solbck(ir,j,:)) )
                constp(ir,j,k)=0.
                conste(ir,j,k)=0.
                IF (locmaxgamma /= 0.) THEN
                   constp(ir,j,k)=1*AIMAG(solbck(ir,j,k))/locmaxgamma
                   conste(ir,j,k)=1*AIMAG(solbck(ir,j,k))/locmaxgamma
                   constv(ir,j,k)=1*AIMAG(solbck(ir,j,k))/locmaxgamma
                END IF

                !PARTICLE TRANSPORT
                cmpfe_k(ir,j,k) = -1._DBL/Ze*1d19*kteta(ir,j)/1._DBL *  ( &
                     & constp(ir,j,k) * fi(ir,j,k) * ( Lcirce(ir,j,k) + Lpiege(ir,j,k) ))

                cmpfi_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1d19*kteta(ir,j) /1._DBL  * ( &
                     & constp(ir,j,k)* fi(ir,j,k) * ( Lcirci(ir,j,:,k) + Lpiegi(ir,j,:,k) ))


                !ENERGY TRANSPORT
                cmefe_k(ir,j,k) = -1._DBL/Ze*1.6d3* kteta(ir,j) * 1._DBL/1._DBL  * ( &
                     & conste(ir,j,k) * fi(ir,j,k) * ( Lecirce(ir,j,k) + Lepiege(ir,j,k) ))

                cmefi_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1.6d3* kteta(ir,j) * 1._DBL/(tau(ir,:)*1._DBL )  * ( &
                     & conste(ir,j,k) * fi(ir,j,k) * ( Lecirci(ir,j,:,k) + Lepiegi(ir,j,:,k) ))

                !ANG MOM TRANSPORT (particle transport * m_s * R) !warning here R=Ro but... is it correct? in principle should integrate R(theta)*L(theta)?
                !answer: probably it only matters in the sense that it should be consistent with the transport code definition
                !added thermal velocity

                IF (rotflagarray(dimx) == 1) THEN
                   cmvfi_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1d19*kteta(ir,j) /1._DBL  * ( &
                        & constv(ir,j,k)* fi(ir,j,k) * Ai(ir,:) * cthi(ir,:) * mp * R0 * ( Lvcirci(ir,j,:,k) + Lvpiegi(ir,j,:,k) ))
                ELSE
                   cmvfi_k(ir,j,:,k) = 0
                ENDIF
                !! SECTION ON ADDITIONAL CALCULATION ON PARTICLE TRANSPORT
                !! ACTIVE ONLY IF phys_meth == 1.0
                IF ( phys_meth /= 0 ) THEN

                   !DIFFUSION TERM

                   cmpfgne_k(ir,j,k) = -1._DBL/Ze* kteta(ir,j) * R0 / Nex(ir) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircgne(ir,j,k) + Lpieggne(ir,j,k) ))
                   cmpfgni_k(ir,j,:,k) = -1._DBL/Zi(ir,:)* kteta(ir,j) * R0 / Nix(ir,:) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircgni(ir,j,:,k) + Lpieggni(ir,j,:,k) ))

                   ! THERMO-DIFFUSION TERM

                   cmpfgte_k(ir,j,k) = -1._DBL/Ze* kteta(ir,j) * Ate(ir) / Nex(ir) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircgte(ir,j,k) + Lpieggte(ir,j,k) ))
                   cmpfgti_k(ir,j,:,k) = -1._DBL/Zi(ir,:)* kteta(ir,j) * Ati(ir,:)/ Nix(ir,:) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircgti(ir,j,:,k) + Lpieggti(ir,j,:,k) ))
                   ! ROTO-DIFFUSION TERM

                   IF (rotflagarray(ir) == 1) THEN
                      cmpfgui_k(ir,j,:,k) = -1._DBL/Zi(ir,:)* kteta(ir,j) * Aui(ir,:)/ Nix(ir,:) * normNL * ( &
                           &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircgui(ir,j,:,k) + Lpieggui(ir,j,:,k) ))
                   ELSE 
                      cmpfgui_k(ir,j,:,k) = 0
                   ENDIF
                   cmpfgui_k(ir,j,:,k) = 0 !HACK to zero out rotodiffusion, which is negligible for main ions anyway, and not trusted for impurities (due to QLK approximations)
                   ! COMPRESSIBILITY TERM

                   cmpfce_k(ir,j,k) = -1._DBL/Ze* kteta(ir,j) * 1._DBL / Nex(ir) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircce(ir,j,k) + Lpiegce(ir,j,k) ))
                   cmpfci_k(ir,j,:,k) = -1._DBL/Zi(ir,:)* kteta(ir,j) * 1._DBL / Nix(ir,:) * normNL * ( &
                        &	constp(ir,j,k) * fi(ir,j,k) * ( Lcircci(ir,j,:,k) + Lpiegci(ir,j,:,k) ))
!!!                
                   IF (phys_meth == 2) THEN
                      !HEAT PINCH THERMO-DIFFUSION TERM   Assume decomposition Q = -chi * n * dT/dr + T*n*V   Then chi = m^2/s and V=m/s

                      cmefgne_k(ir,j,k) = -1._DBL/Ze*1.6d3* kteta(ir,j) * 1._DBL/1._DBL  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircgne(ir,j,k) + Lepieggne(ir,j,k) )) * Ane(ir) / (Nex(ir)*1d19 * Tex(ir)*1d3*qe)* normNL

                      cmefgni_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1.6d3* kteta(ir,j) * 1._DBL/(tau(ir,:)*1._DBL )  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircgni(ir,j,:,k) + Lepieggni(ir,j,:,k) )) * Ani(ir,:) / (Nix(ir,:)*1d19 * Tix(ir,:)*1d3*qe)* normNL

                      !HEAT PINCH ROTO-DIFFUSION TERM   Assume decomposition Q = -chi * n * dT/dr + T*n*V   Then chi = m^2/s and V=m/s
                      !                      normNL=603._DBL/3.1633;
                      IF (rotflagarray(ir) == 1) THEN

                         cmefgui_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1.6d3* kteta(ir,j) * 1._DBL/(tau(ir,:)*1._DBL )  * ( &
                              & conste(ir,j,k) * fi(ir,j,k) * ( Lecircgui(ir,j,:,k) + Lepieggui(ir,j,:,k) )) * Aui(ir,:) / (Nix(ir,:)*1d19 * Tix(ir,:)*1d3*qe)* normNL
                      ELSE
                         cmefgui_k(ir,j,:,k)=0
                      ENDIF

                      cmefgui_k(ir,j,:,k) = 0 !HACK to zero out rotodiffusion, which is negligible for main ions anyway, and not trusted for impurities (due to QLK approximations)

                      ! HEAT DIFFUSION TERM. Defined such that we multiple by -n*dT/dr for flux

                      cmefgte_k(ir,j,k) = -1._DBL/Ze*1.6d3* kteta(ir,j) * 1._DBL/1._DBL  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircgte(ir,j,k) + Lepieggte(ir,j,k) )) / (Nex(ir)*1d19) * R0 / (Tex(ir)*1e3*qe)* normNL

                      cmefgti_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1.6d3* kteta(ir,j) * 1._DBL/(tau(ir,:)*1._DBL )  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircgti(ir,j,:,k) + Lepieggti(ir,j,:,k) )) / (Nix(ir,:)*1d19) * R0 / (Tix(ir,:)*1d3*qe)* normNL

                      ! HEAT PINCH COMPRESSIBILITY TERM

                      cmefce_k(ir,j,k) = -1._DBL/Ze*1.6d3* kteta(ir,j) * 1._DBL/1._DBL  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircce(ir,j,k) + Lepiegce(ir,j,k) )) / (Nex(ir)*1d19 * Tex(ir)*1d3*qe)* normNL

                      cmefci_k(ir,j,:,k) = -1._DBL/Zi(ir,:)*1.6d3* kteta(ir,j) * 1._DBL/(tau(ir,:)*1._DBL )  * ( &
                           & conste(ir,j,k) * fi(ir,j,k) * ( Lecircci(ir,j,:,k) + Lepiegci(ir,j,:,k) )) / (Nix(ir,:)*1d19 * Tix(ir,:)*1d3*qe)* normNL
                   ENDIF
                ENDIF
             ENDDO  !END SUM OVER SOLUTIONS
             cmpfe(ir,j)=SUM(cmpfe_k(ir,j,:))
             cmefe(ir,j)=SUM(cmefe_k(ir,j,:))

             DO ion = 1,nions
                cmpfi(ir,j,ion)=SUM(cmpfi_k(ir,j,ion,:))
                cmefi(ir,j,ion)=SUM(cmefi_k(ir,j,ion,:))
                cmvfi(ir,j,ion)=SUM(cmvfi_k(ir,j,ion,:))
             ENDDO
             IF ( phys_meth /= 0.0 ) THEN
                cmpfgne(ir,j)=SUM(cmpfgne_k(ir,j,:))
                cmpfgte(ir,j)=SUM(cmpfgte_k(ir,j,:))
                cmpfce(ir,j)=SUM(cmpfce_k(ir,j,:))
                IF (phys_meth == 2) THEN
                   cmefgne(ir,j)=SUM(cmefgne_k(ir,j,:))
                   cmefgte(ir,j)=SUM(cmefgte_k(ir,j,:))
                   cmefce(ir,j)=SUM(cmefce_k(ir,j,:))
                ENDIF
                DO ion=1,nions
                   cmpfgni(ir,j,ion)=SUM(cmpfgni_k(ir,j,ion,:))
                   cmpfgti(ir,j,ion)=SUM(cmpfgti_k(ir,j,ion,:))
                   cmpfgui(ir,j,ion)=SUM(cmpfgui_k(ir,j,ion,:))
                   cmpfci(ir,j,ion)=SUM(cmpfci_k(ir,j,ion,:))
                   IF (phys_meth == 2) THEN
                      cmefgni(ir,j,ion)=SUM(cmefgni_k(ir,j,ion,:))
                      cmefgti(ir,j,ion)=SUM(cmefgti_k(ir,j,ion,:))
                      cmefgui(ir,j,ion)=SUM(cmefgui_k(ir,j,ion,:))
                      cmefci(ir,j,ion)=SUM(cmefci_k(ir,j,ion,:))
                   ENDIF
                END DO
             ENDIF
          END DO !end do over wavenumbers

          !the davint numerical integrator is used for the routines below

          IF ( (kthr(ir,1) > ETGk) .OR. (dimn ==1) ) THEN  !single mode or ETG-scale. Should not integrate from zero
             lowlim=kthr(ir,1)
          ELSE
             lowlim=0._DBL 
          ENDIF

          !Remove any particle transport due to ETG to maintain quasineutrality
          WHERE(kthetarhos > ETGk)
             cmpfe(ir,:) = 0.
             cmpfgne(ir,:) = 0.
             cmpfgte(ir,:) = 0.
             cmpfce(ir,:) = 0.
          ENDWHERE

          ! Particle flux using all roots 
          xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfe(ir,:)/)
          IF (dimn == 1) THEN 
             pfe(ir)=cmpfe(ir,1)
          ELSE
             CALL davint (xint, yint ,dimn+1,lowlim, kthr(ir,dimn), pfe(ir), ifailloc,27)
          ENDIF

          ! Total particle diffusivity Gyro-Bohm using all roots. Assumes that all particle transport is diagonal (typically not the case)
          dpfe(ir) = (pfe(ir)/(Nex(ir)*1e19/Rmin(ir)))/chi_GB(ir)

          ! Energy flux using all roots
          xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefe(ir,:)/)
          IF (dimn == 1) THEN 
             efe(ir)=cmefe(ir,1)
          ELSE
             CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),efe(ir),ifailloc,28)
          ENDIF
          ! Energy diffusivity using all unstable roots
          defe(ir) = (efe(ir)/(Nex(ir)*1e19*Tex(ir)*1e3*qe/Rmin(ir)))/chi_GB(ir)

          !Only ETG scales particle and heat transport (if they exist)
          IF (ETGind > 0) THEN
             xint= (/0._DBL,kthr(ir,ETGind:dimn)/) ; yint=(/0._DBL,cmefe(ir,ETGind:dimn)/)
             IF (dimn == 1) THEN 
                efeETG(ir)=cmefe(ir,1)
             ELSE
                CALL davint (xint, yint ,dimn-ETGind+2,kthr(ir,ETGind), kthr(ir,dimn), efeETG(ir), ifailloc,30)
             ENDIF
          ELSE
             efeETG(ir) = 0.
          ENDIF
          defeETG(ir) = (efeETG(ir)/(Nex(ir)*1e19*Tex(ir)*1e3*qe/Rmin(ir)))/chi_GB(ir)

          ! Ang mom flux using all roots

          DO ion=1,nions
             !Remove any residual ETG particle transport to maintain quasineutrality
             WHERE(kthetarhos > ETGk) 
                cmpfi(ir,:,ion) = 0.
                cmpfgni(ir,:,ion) = 0.
                cmpfgti(ir,:,ion) = 0.
                cmpfgui(ir,:,ion) = 0.
                cmpfci(ir,:,ion) = 0.
             ENDWHERE

             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfi(ir,:,ion)/)
             IF (dimn == 1) THEN 
                pfi(ir,ion)=cmpfi(ir,1,ion)
             ELSE
                CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),pfi(ir,ion),ifailloc,32)
             ENDIF
             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefi(ir,:,ion)/)
             IF (dimn == 1) THEN 
                efi(ir,ion)=cmefi(ir,1,ion)
             ELSE
                CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),efi(ir,ion),ifailloc,33)        
             ENDIF
             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmvfi(ir,:,ion)/)
             IF (dimn == 1) THEN 
                vfi(ir,ion)=cmvfi(ir,1,ion)
             ELSE
                CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vfi(ir,ion),ifailloc,34)        
             ENDIF
             dpfi(ir,ion) = (pfi(ir,ion)/(Nix(ir,ion)*1e19/Rmin(ir)))/chi_GB(ir)
             defi(ir,ion) = (efi(ir,ion)/(Nix(ir,ion)*1e19*Tix(ir,ion)*1e3*qe/Rmin(ir)))/chi_GB(ir)
             IF (rot_flag==2) Aui(ir,ion)=Auiorig(ir,ion) !reinstate original for momentum diffusivity definition
             dvfi(ir,ion) = (vfi(ir,ion)/(Nix(ir,ion)*1e19*Ai(ir,ion)*mp*cthi(ir,ion)*R0/Rmin(ir)))/chi_GB(ir)
          ENDDO

          IF ( phys_meth /= 0 ) THEN

             ! Total diffusion coefficient
             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgne(ir,:)/)
             IF (dimn == 1) THEN 
                dffte(ir)=cmpfgne(ir,1)
             ELSE
                CALL davint (xint, yint ,dimn+1,lowlim,kthr(ir,dimn),dffte(ir),ifailloc,35)
             ENDIF
             ! Total thermo-diffusion coefficient
             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgte(ir,:)/)
             IF (dimn == 1) THEN 
                vthte(ir)=cmpfgte(ir,1)
             ELSE
                CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vthte(ir),ifailloc,36)
             ENDIF
             ! Total compressibility coefficient
             xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfce(ir,:)/)
             IF (dimn == 1) THEN 
                vcpte(ir)=cmpfce(ir,1)
             ELSE
                CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vcpte(ir),ifailloc,38)
             ENDIF

!!! HEAT PINCH TERMS
             IF (phys_meth == 2) THEN
                ! Thermo-diffusion coefficient
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefgne(ir,:)/)
                IF (dimn == 1) THEN 
                   deffte(ir)=cmefgne(ir,1)
                ELSE
                   CALL davint (xint, yint ,dimn+1,lowlim,kthr(ir,dimn),deffte(ir),ifailloc,39)
                ENDIF
                ! Thermal conductivity coefficient
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefgte(ir,:)/)
                IF (dimn == 1) THEN 
                   vethte(ir)=cmefgte(ir,1)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vethte(ir),ifailloc,40)
                ENDIF
                ! Total compressibility coefficient
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefce(ir,:)/)
                IF (dimn == 1) THEN 
                   vecpte(ir)=cmefce(ir,1)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vecpte(ir),ifailloc,42)
                ENDIF

                IF (separateflux .EQV. .TRUE.) THEN

                   IF (ETGind > 0) THEN
                      xint= (/0._DBL,kthr(ir,ETGind:dimn)/) ; yint=(/0._DBL,cmefgne(ir,ETGind:dimn)/)
                      IF (dimn == 1) THEN 
                         deffteETG(ir)=cmefgne(ir,1)
                      ELSE
                         CALL davint (xint, yint ,dimn-ETGind+2,kthr(ir,ETGind), kthr(ir,dimn), deffteETG(ir), ifailloc,39)
                      ENDIF
                   ELSE
                      deffteETG(ir) = 0.
                   ENDIF

                   IF (ETGind > 0) THEN
                      xint= (/0._DBL,kthr(ir,ETGind:dimn)/) ; yint=(/0._DBL,cmefgte(ir,ETGind:dimn)/)
                      IF (dimn == 1) THEN 
                         vethteETG(ir)=cmefgte(ir,1)
                      ELSE
                         CALL davint (xint, yint ,dimn-ETGind+2,kthr(ir,ETGind), kthr(ir,dimn), vethteETG(ir), ifailloc,40)
                      ENDIF
                   ELSE
                      vethteETG(ir) = 0.
                   ENDIF

                   IF (ETGind > 0) THEN
                      xint= (/0._DBL,kthr(ir,ETGind:dimn)/) ; yint=(/0._DBL,cmefce(ir,ETGind:dimn)/)
                      IF (dimn == 1) THEN 
                         vecpteETG(ir)=cmefce(ir,1)
                      ELSE
                         CALL davint (xint, yint ,dimn-ETGind+2,kthr(ir,ETGind), kthr(ir,dimn), vecpteETG(ir), ifailloc,42)
                      ENDIF
                   ELSE
                      vecpteETG(ir) = 0.
                   ENDIF

                ENDIF

             ENDIF

             DO ion=1,nions
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgni(ir,:,ion)/)
                IF (dimn == 1) THEN 
                   dffti(ir,ion)=cmpfgni(ir,1,ion)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),dffti(ir,ion),ifailloc,43)
                ENDIF
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgti(ir,:,ion)/)
                IF (dimn == 1) THEN 
                   vthti(ir,ion)=cmpfgti(ir,1,ion)
                ELSE            
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vthti(ir,ion),ifailloc,44)
                ENDIF
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgui(ir,:,ion)/)
                IF (dimn == 1) THEN 
                   vrdti(ir,ion)=cmpfgui(ir,1,ion)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vrdti(ir,ion),ifailloc,45)
                ENDIF
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfci(ir,:,ion)/)
                IF (dimn == 1) THEN 
                   vcpti(ir,ion)=cmpfci(ir,1,ion)
                ELSE
                   CALL davint (xint, yint ,dimn+1,lowlim,kthr(ir,dimn),vcpti(ir,ion),ifailloc,46)
                ENDIF

                ! Transport coefficients including 2D centrifugal and temp anisotropy effects. All depends on precalculated ecoefsgau. 
                ! Fluxes in transport codes should be built from cftrans when 2D effects are desired

                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgni(ir,:,ion)*ecoefsgau(ir,:,ion,0)/)
                IF (dimn == 1) THEN 
                   cftrans(ir,ion,1)=cmpfgni(ir,1,ion)*ecoefsgau(ir,1,ion,0)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,1),ifailloc,47) !Generalized diffusivity
                ENDIF
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgti(ir,:,ion)*ecoefsgau(ir,:,ion,0)/)
                IF (dimn == 1) THEN 
                   cftrans(ir,ion,2)=cmpfgti(ir,1,ion)*ecoefsgau(ir,1,ion,0)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,2),ifailloc,48) !Generalized thermo-pinch
                ENDIF
                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfci(ir,:,ion)*ecoefsgau(ir,:,ion,0)/)
                IF (dimn == 1) THEN 
                   cftrans(ir,ion,3)=cmpfci(ir,1,ion)*ecoefsgau(ir,1,ion,0)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,3),ifailloc,49) !Generalized compression pinch
                ENDIF

                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmpfgui(ir,:,ion)*ecoefsgau(ir,:,ion,0)/)
                IF (dimn == 1) THEN 
                   cftrans(ir,ion,4)=cmpfgui(ir,1,ion)*ecoefsgau(ir,1,ion,0)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,4),ifailloc,50) !Generalized roto-diffusion pinch
                ENDIF

                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,1./R0*Ati(ir,ion)*cmpfgni(ir,:,ion)*( ecoefsgau(ir,:,ion,1) - (mi(ir,ion)*(omegator(ir)*Ro(ir))**2)/(2.*qe*1d3*Tix(ir,ion))*ecoefsgau(ir,:,ion,3) )  /)
                IF (dimn == 1) THEN 
                   cftrans(ir,ion,5)=1./R0*Ati(ir,ion)*cmpfgni(ir,1,ion)*( ecoefsgau(ir,1,ion,1) - (mi(ir,ion)*(omegator(ir)*Ro(ir))**2)/(2.*qe*1d3*Tix(ir,ion))*ecoefsgau(ir,1,ion,3) )
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,5),ifailloc,51) !"2D thermopinch" 
                ENDIF

                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,1./R0*cmpfgni(ir,:,ion)*2. * Machi(ir,ion)*Aupar(ir)*cref(ir)/cthi(ir,ion)*SQRT(1+(epsilon(ir)/qx(ir))**2) *ecoefsgau(ir,:,ion,3)/)
                IF (dimn == 1) THEN
                   cftrans(ir,ion,6)=1./R0*cmpfgni(ir,1,ion)*2. * Machi(ir,ion)*Aupar(ir)*cref(ir)/cthi(ir,ion)*SQRT(1+(epsilon(ir)/qx(ir))**2) *ecoefsgau(ir,1,ion,3)
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,6),ifailloc,52) !"2D rotodiffusion" 
                ENDIF

                xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,1./R0*cmpfgni(ir,:,ion)*( ecoefsgau(ir,:,ion,2) - ecoefsgau(ir,:,ion,7) - 2.*Machi(ir,ion)**2 * ( ecoefsgau(ir,:,ion,8) + ecoefsgau(ir,:,ion,9)/2.) )/)
                IF (dimn == 1) THEN
                   cftrans(ir,ion,7)=1./R0*cmpfgni(ir,1,ion)*( ecoefsgau(ir,1,ion,2) - ecoefsgau(ir,1,ion,7) - 2.*Machi(ir,ion)**2 * ( ecoefsgau(ir,1,ion,8) + ecoefsgau(ir,1,ion,9)/2.) )
                ELSE
                   CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),cftrans(ir,ion,7),ifailloc,53) !"2D pure pinch"
                ENDIF

                IF (phys_meth == 2) THEN
                   xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefgni(ir,:,ion)/)
                   IF (dimn == 1) THEN 
                      deffti(ir,ion)=cmefgni(ir,1,ion)
                   ELSE
                      CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),deffti(ir,ion),ifailloc,54)
                   ENDIF
                   xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefgti(ir,:,ion)/)
                   IF (dimn == 1) THEN 
                      vethti(ir,ion)=cmefgti(ir,1,ion)
                   ELSE
                      CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),vethti(ir,ion),ifailloc,55)
                   ENDIF
                   xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefgui(ir,:,ion)/)
                   IF (dimn == 1) THEN 
                      verdti(ir,ion)=cmefgui(ir,1,ion)
                   ELSE
                      CALL davint (xint, yint, dimn+1,lowlim,kthr(ir,dimn),verdti(ir,ion),ifailloc,56)
                   ENDIF
                   xint= (/0._DBL,kthr(ir,:)/) ; yint=(/0._DBL,cmefci(ir,:,ion)/)
                   IF (dimn == 1) THEN 
                      vecpti(ir,ion)=cmefci(ir,1,ion)
                   ELSE                
                      CALL davint (xint, yint ,dimn+1,lowlim,kthr(ir,dimn),vecpti(ir,ion),ifailloc,57)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF ! end of statement on additional calculation

       END DO  !end big cycle on radial position

       ! NORMALISATION CONSTANT, BENCHMARK WITH GYRO

       !Particle transport, all roots
       dpfe=normNL*dpfe
       dpfi=normNL*dpfi

       !Energy transport, all roots
       defe=normNL*defe
       defeETG=normNL*defeETG

       defi=normNL*defi

       !Ang mom transport, all roots
       dvfi=normNL*dvfi

       ! CREATE ADDITIONAL FINAL OUTPUT ARRAYS
       !      IF (gg == 3) THEN
       doit=0
       DO ir=1,dimx

          IF (myrank == doit) THEN ! distribute independent loop indices to tasks
             doit=doit+1; 
             IF (doit==nproc) doit=0 
          ELSE
             doit=doit+1; 
             IF (doit==nproc) doit=0 
             CYCLE
          ENDIF

          ipf_SI(ir,:) = pfi(ir,:)*normNL
          epf_SI(ir) = pfe(ir)*normNL

          ipf_GB(ir,:) = dpfi(ir,:)
          epf_GB(ir) = dpfe(ir)

          ief_SI(ir,:) = efi(ir,:)*normNL
          eef_SI(ir) = efe(ir)*normNL
          eefETG_SI(ir) = efeETG(ir)*normNL        

          ief_GB(ir,:) = defi(ir,:)
          eef_GB(ir) = defe(ir)
          eefETG_GB(ir) = defeETG(ir)

          ivf_SI(ir,:) = vfi(ir,:)*normNL
          ivf_GB(ir,:) = dvfi(ir,:)

          !FLUX SPECTRA
          eef_cm(ir,:) = normNL*cmefe(ir,:)
          epf_cm(ir,:) = normNL*cmpfe(ir,:)

          DO ion=1,nions
             ief_cm(ir,:,ion) = normNL*cmefi(ir,:,ion)
             ivf_cm(ir,:,ion) = normNL*cmvfi(ir,:,ion)
             ipf_cm(ir,:,ion) = normNL*cmpfi(ir,:,ion)
          ENDDO

          !WRITING ADDITIONAL OUTPUT IF QLK_input.phys_meth == 1.0

          IF ( phys_meth /= 0 ) THEN

             ! Particle diffusion coefficients
             dfe_SI(ir) = dffte(ir)
             dfi_SI(ir,:) = dffti(ir,:)

             !! Particle thermo-diffusion pinch
             vte_SI(ir) = vthte(ir)
             vti_SI(ir,:) = vthti(ir,:)

             !! Compressibility pinches
             vce_SI(ir) = vcpte(ir)
             vci_SI(ir,:) = vcpti(ir,:)

             !! roto diffusion terms
             vri_SI(ir,:) = vrdti(ir,:)

             dfe_GB(ir) = dfe_SI(ir)/chi_GB(ir)
             dfi_GB(ir,:) = dfi_SI(ir,:)/chi_GB(ir)

             vte_GB(ir) = vte_SI(ir)*Rmin(ir)/chi_GB(ir)
             vti_GB(ir,:) = vti_SI(ir,:)*Rmin(ir)/chi_GB(ir)

             vce_GB(ir) = vce_SI(ir)*Rmin(ir)/chi_GB(ir)
             vci_GB(ir,:) = vci_SI(ir,:)*Rmin(ir)/chi_GB(ir)

             vri_GB(ir,:) = vri_SI(ir,:)*Rmin(ir)/chi_GB(ir)

             !! check on particle fluxes
             cke(ir) = 1d2* ( epf_SI(ir) - ( dffte(ir)*Ane(ir)*Nex(ir)*1d19/R0 + & 
                  &	Nex(ir)*1d19*(vthte(ir)+vcpte(ir)) ) )/ (epf_SI(ir)+epsD)
             cki(ir,:) = 1d2* ( ipf_SI(ir,:) - ( dffti(ir,:)*Ani(ir,:)*Nix(ir,:)*1d19/R0 + &
                  &	Nix(ir,:)*1d19*(vthti(ir,:)+vcpti(ir,:)+vrdti(ir,:)) ) ) / (ipf_SI(ir,:)+epsD)

             IF (phys_meth == 2) THEN
                ! Heat thermodiffusion pinch 
                vene_SI(ir) = deffte(ir)
                veni_SI(ir,:) = deffti(ir,:)

                !Heat conductivity
                chiee_SI(ir) = vethte(ir)
                chiei_SI(ir,:) = vethti(ir,:)

                ! Heat compressibility pinch
                vece_SI(ir) = vecpte(ir)
                veci_SI(ir,:) = vecpti(ir,:)

                ! Heat roto-diff pinch
                veri_SI(ir,:) = verdti(ir,:)

                chiee_GB(ir) = chiee_SI(ir)/chi_GB(ir)
                chiei_GB(ir,:) = chiei_SI(ir,:)/chi_GB(ir)

                vene_GB(ir) = vene_SI(ir)*Rmin(ir)/chi_GB(ir)
                veni_GB(ir,:) = veni_SI(ir,:)*Rmin(ir)/chi_GB(ir)

                vece_GB(ir) = vece_SI(ir)*Rmin(ir)/chi_GB(ir)
                veci_GB(ir,:) = veci_SI(ir,:)*Rmin(ir)/chi_GB(ir)

                veri_GB(ir,:) = veri_SI(ir,:)*Rmin(ir)/chi_GB(ir)

                IF (separateflux .EQV. .TRUE.) THEN
                   veneETG_SI(ir) = deffteETG(ir)
                   chieeETG_SI(ir) = vethteETG(ir)
                   veceETG_SI(ir) = vecpteETG(ir)
                   chieeETG_GB(ir) = chieeETG_SI(ir)/chi_GB(ir)
                   veneETG_GB(ir) = veneETG_SI(ir)*Rmin(ir)/chi_GB(ir)
                   veceETG_GB(ir) = veceETG_SI(ir)*Rmin(ir)/chi_GB(ir)
                ENDIF

                !! check on energy fluxes
                ceke(ir) = 1d2* ( eef_SI(ir) - ( chiee_SI(ir)*Ate(ir)/R0*Tex(ir)*qe*1d3*Nex(ir)*1d19 + & 
                     &	Nex(ir)*1d19*Tex(ir)*qe*1d3*(vene_SI(ir)+vece_SI(ir)) ) )/ (eef_SI(ir)+epsD)

                ceki(ir,:) = 1d2* ( ief_SI(ir,:) - ( chiei_SI(ir,:)*Ati(ir,:)/R0*Tix(ir,:)*qe*1d3*Nix(ir,:)*1d19 + & 
                     &	Nix(ir,:)*1d19*Tix(ir,:)*qe*1d3*(veni_SI(ir,:)+veci_SI(ir,:)+veri_SI(ir,:)) ) )/ (ief_SI(ir,:)+epsD)
             ENDIF

          ENDIF!! end of statement on additional calculation

          IF (gg==1) THEN !save ion mode only output
             ion_epf_GB(ir) = dpfe(ir)
             ion_ipf_GB(ir,:) = dpfi(ir,:)
             ion_eef_GB(ir) = defe(ir)
             ion_ief_GB(ir,:) = defi(ir,:)
             ion_ivf_GB(ir,:) = dvfi(ir,:)
          END IF
          IF (gg==2) THEN !save electron mode only output
             ele_epf_GB(ir) = dpfe(ir)
             ele_ipf_GB(ir,:) = dpfi(ir,:)
             ele_eef_GB(ir) = defe(ir)
             ele_ief_GB(ir,:) = defi(ir,:)
             ele_ivf_GB(ir,:) = dvfi(ir,:)
          END IF

          !Testing for existence of ion or electron modes
          IF (gg==3) THEN
             IF ( (ABS(epf_GB(ir)) > epsD) .OR. (ABS(eef_GB(ir)) > epsD) .OR. (ABS(ipf_GB(ir,1)) > epsD) & 
                  & .OR. (ABS(ief_GB(ir,1)) > epsD)) THEN !If an instability is active on ele or main ion then check following
                IF ( (ele_epf_GB(ir)/(epf_GB(ir)+epsD) < impfac) .AND. (ele_eef_GB(ir)/(eef_GB(ir)+epsD) < impfac) & 
                     & .AND. (ele_ipf_GB(ir,1)/(ipf_GB(ir,1)+epsD) < impfac) .AND. (ele_ief_GB(ir,1)/(ief_GB(ir,1)+epsD) < impfac)) THEN
                   modeflag(ir)=1 !At this radial location, only ion modes are important
                ENDIF
                IF ( (ion_epf_GB(ir)/(epf_GB(ir)+epsD) < impfac) .AND. (ion_eef_GB(ir)/(eef_GB(ir)+epsD) < impfac) & 
                     & .AND. (ion_ipf_GB(ir,1)/(ipf_GB(ir,1)+epsD) < impfac) .AND. (ion_ief_GB(ir,1)/(ief_GB(ir,1)+epsD) < impfac)) THEN
                   modeflag(ir)=2 !At this radial location, only electron modes are important
                ENDIF
             ENDIF
          ENDIF
       END DO !end of radial cycle

    ENDDO !end do on gg


  END SUBROUTINE saturation

END MODULE mod_saturation
