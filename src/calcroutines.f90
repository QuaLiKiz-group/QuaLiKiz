MODULE calcroutines
  USE kind
  USE datcal
  USE datmat
  USE FLRterms !module contain functions defining all the FLR terms 
  USE mod_fluidsol !module which calculates the fluid growth rate and frequencies
  USE mod_contour !module which contains contour routines
  USE mod_make_io
  USE mod_fonct
  USE QLflux
  USE asymmetry
  USE nanfilter

  IMPLICIT NONE

CONTAINS
  SUBROUTINE calc(p,nu)
    INTEGER, INTENT(IN) :: p,nu

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
    REAL(KIND=DBL), DIMENSION(1) :: xmin, xmax, intout_cub, acc_cub

    CHARACTER(len=20) :: fmtn
    INTEGER :: myunit=700

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

    !Calculate mode width and shift according to Cottier model. omeflu, mwidth, and mshift are written inside the subroutines
    CALL cot_fluidsol(p,nu)! with rotation solflu used for width and shift calculations
    cot_solflu(p,nu) = omeflu ; solflu(p,nu) = omeflu
    CALL cot_calcwidth(ETG_flag(nu),p,nu)
    cot_modewidth(p,nu) = mwidth !The modewidth variable is output, while mwidth is used within code
    CALL cot_calcshift(ETG_flag(nu),p,nu)
    cot_modeshift(p,nu) = mshift !The modeshift variable is output, while shift is used in code

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

    !Calculate old mode width and shift from before rotation times
    CALL old_calcwidth(ETG_flag(nu),p)
    mshift = 0.0 

    ! Keeping old results the same
    old_modewidth(p,nu) = mwidth 
    old_modeshift(p,nu) = mshift

    solflu(p,nu) = ana_solflu(p,nu)

    !! choose which one to actually keep for use in the code
    modewidth(p,nu) = jon_modewidth(p,nu)
    mwidth = modewidth(p,nu)
    modeshift(p,nu) = jon_modeshift(p,nu)
    modeshift2(p,nu) = mshift2
    mshift = modeshift(p,nu); 

    !Calculate the width corresponding to |phi(x)| (in case eigenfunction is complex)
    widthhat = ABS(mwidth)**2 / SQRT(REAL(mwidth**2))

    !    mshift=CMPLX(-0.1E-02,0.2E-02)
    !    mshift=CMPLX(0.1E-02,-0.2E-02)

    !    mshift=CMPLX(REAL(mshift)/100,AIMAG(mshift))

    !    mshift=CMPLX(0.,0.)

    ! Calculate flux-surface-averaging of poloidal asymmetry coefficients including eigenmode
    CALL makeecoefsgau(p,nu)

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

    !DEBUGGING
!!$    WRITE(*,*) 'Joe2p',Joe2p
!!$    WRITE(*,*) 'Joe2c',Joe2c
!!$    WRITE(*,*) 'Joi2p',Joi2p
!!$    WRITE(*,*) 'Joi2c',Joi2c
    !STOP
    !*************************

    !used to pass variables in <Vpar^n> calculations below
    plam = p
    nulam = nu
    maxklam =   ABS(pi/(distan(p,nu)*normkr)) 
    minklam = - maxklam
    
    xmin(1) = 0._DBL
    xmax(1) = 1._DBL - 2.*epsilon(p)

    !Calculate sin^2(theta/2) weighted against the eigenfunction. At the moment not used
!!$    sin2th = d01ahf(minklam,maxklam,relacc1,npts,relerr,sin2thint,lw,ifailloc)
!!$    IF (ifailloc /= 0) THEN
!!$       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
!!$            &'. Abnormal termination of sin2th integration at p=',p,', nu=',nu
!!$    ENDIF

!!$    alamnorm=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alamnormint,lw,ifailloc) !Normalization for pitch angle integrations
!!$    IF (ifailloc /= 0) THEN
!!$       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
!!$            &'. Abnormal termination of alamnorm integration at p=',p,', nu=',nu
!!$    ENDIF

    alamnorm = fc(p) !to be consistent with passing particle fraction

    !alam1=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam1int,lw,ifailloc)/alamnorm !pitch angle average of sqrt(1-lambda*b)
    ifailloc = pcubature(1, alam1int_cubature, 1, xmin, xmax, npts, 0._DBL, relacc1, 1, intout_cub, acc_cub)
    alam1 = intout_cub(1) / alamnorm
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam1 integration at p=',p,', nu=',nu
    ENDIF
!!$    WRITE(*,*) 'alamnorm,fc(p)',p,nu,alamnorm,fc(p)
!!$    WRITE(*,*) 'alam1,alam1*alamnorm/fc(p),1/SQRT(3)',alam1,alam1*alamnorm/fc(p),1./SQRT(3.)

    !pitch angle average of (1-lambda*b)
    !alam2=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam2int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)
    ifailloc = pcubature(1, alam2int_cubature, 1, xmin, xmax, npts, 0._DBL, relacc1, 1, intout_cub, acc_cub)
    alam2 = intout_cub(1) / alamnorm
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam2 integration at p=',p,', nu=',nu
    ENDIF

    !alam3=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam3int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^3/2
    ifailloc = pcubature(1, alam3int_cubature, 1, xmin, xmax, npts, 0._DBL, relacc1, 1, intout_cub, acc_cub)
    alam3 = intout_cub(1) / alamnorm
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam3 integration at p=',p,', nu=',nu
    ENDIF

    !alam4=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam4int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^2
    ifailloc = pcubature(1, alam4int_cubature, 1, xmin, xmax, npts, 0._DBL, relacc1, 1, intout_cub, acc_cub)
    alam4 = intout_cub(1) / alamnorm
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam4 integration at p=',p,', nu=',nu
    ENDIF

    !alam5=d01ahf(0.,1.-2.*epsilon(p),relacc1,npts,relerr,alam5int,lw,ifailloc)/alamnorm !pitch angle average of (1-lambda*b)^5/2
    ifailloc = pcubature(1, alam5int_cubature, 1, xmin, xmax, npts, 0._DBL, relacc1, 1, intout_cub, acc_cub)
    alam5 = intout_cub(1) / alamnorm
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of alam5 integration at p=',p,', nu=',nu
    ENDIF

!!$    alam1=1. !pitch angle average of (1-lambda*b)^1/2
!!$    alam2=1. !pitch angle average of (1-lambda*b)
!!$    alam3=1. !pitch angle average of (1-lambda*b)^3/2
!!$    alam4=1. !pitch angle average of (1-lambda*b)^2
!!$    alam5=1. !pitch angle average of (1-lambda*b)^5/2

    !Set the transit frequency
    !    qRd = qx(p)*Ro(p)*d*SQRT(3._DBL)
    qRd = qx(p)*Ro(p)*d*1./alam1

    Athe=widthhat*cthe(p)/qRd 
    Athi(:)=widthhat*cthi(p,:)/qRd

    !Set the bounce frequency
    omega2bar = pi/2.*SQRT(x(p)*Rmin(p)/(2.*Ro(p)))

    !Initialize ma value used in contour choice based on fluid solution
    ma = MAX(AIMAG(solflu(p,nu)),10.)
    om0 = ma*ci
    C    = om0 *0.5

    !Set maximum bound of contour locations based on diamagnetic frequency
    ommax(p,nu) = om0 + REAL(solflu(p,nu))/2. !divisor sets max omega ratio to diamagnetic frequency

    !Set maximum bound of contour locations based on heritage assumptions (not clear to me (JC))

    !    IF (ETG_flag(nu) .EQV. .FALSE.) THEN
    !       ommax(p,nu) = om0+MAX(2.*ABS(Tex(p)),2.*ABS(Tix(p,1)/Zi(p,1)),2.*ABS(Athi(1)/(nwg)))
    !    ELSE
    !       ommax(p,nu) = om0+MIN(MAX(2.*ABS(Tex(p)),1.5*ABS(Athe/(nwg))),12.)
    !    ENDIF

    omegmax=ommax(p,nu) ! used in calculsol for maximum boundary of allowed solution

    IF (runcounter /=0) THEN !launch newton solver from previous run
       DO j=1,numsols
          IF ( ABS(AIMAG(oldsol(p,nu,j))) < epsD ) THEN !There was no solution before, so stay with 0
             soll(j)   = (0.,0.)
             fdsoll(j) = (0.,0.)
          ELSE

             IF ( (gkw_is_nan(AIMAG(oldsol(p,nu,j)))) .OR. (gkw_is_nan(REAL(oldsol(p,nu,j))))) THEN
                IF (verbose .EQV. .TRUE.) THEN 
                   WRITE(stderr,'(A,I7,A,I2,A)') 'Jump to Newton phase: old solution had a NaN (how did that happen)! Skipping solution. (p,nu)=(',p,',',nu,')'
                ENDIF
                soll(j)   = (0.,0.)
                fdsoll(j)   = (0.,0.)
                CYCLE
             ENDIF

             IF ( (AIMAG(oldsol(p,nu,j)) < 0. ) .OR. (ABS(AIMAG(oldsol(p,nu,j))) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(oldsol(p,nu,j))) > ABS(REAL(ommax(p,nu))))  ) THEN
                IF (verbose .EQV. .TRUE.) THEN
                   WRITE(stderr,'(A,I7,A,I2,A)') 'Jump to Newton phase: old solution outside of allowed contour range (how did that happen?) Skipping solution. (p,nu)=(',p,',',nu,')'
                ENDIF
                soll(j)   = (0.,0.)
                fdsoll(j)   = (0.,0.)
                CYCLE
             ENDIF

             IF ( ( rho(p) >= rhomin ) .AND. ( rho(p) <= rhomax) ) THEN !check if rho is within the defined range, otherwise return zero
                CALL calcfonct(p, nu, oldsol(p,nu,j), fonout) !Get new distance from solution from old solution and new input parameters
                !CALL newton(p, nu, oldsol(p,nu,j), fonout, newsol, newfdsol) !Now refine the solution to the new solution
                CALL broyden(p, nu, oldsol(p,nu,j), fonout, newsol, newfdsol) !Now refine the solution to the new solution
             ELSE
                newsol = 0.
                newfdsol = 0.
             ENDIF
             !If, after the Newton refinement, the solution is too coarse, or outside the preset boundaries, then we abandon the solution
             IF (ABS(newfdsol)> 0.5  .OR. &
                  ABS(AIMAG(newsol))>(ABS(AIMAG(omegmax))) .OR. &
                  ABS(REAL(newsol)) > (ABS(REAL(omegmax)))) THEN
                soll(j)   = (0.,0.)
                fdsoll(j) = (0.,0.)
                IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I3,A,I3,A,G15.7,A,2G15.7,A)') &
                     'Jumped straight to Newton for p=',p,' nu=,',nu,' , ABS(newfdsol)=',ABS(newfdsol),' with newsol=',newsol,' and the solution was abandoned'
                !If the solution is stable, then the solution is zeroed out. We can't have stable solutions since we haven't included the analytic continuation
             ELSEIF (AIMAG(newsol)<0.) THEN
                soll(j)   = (0.,0.)
                fdsoll(j) = (0.,0.)
                IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I3,A,I3,A)') 'For p=',p,' nu=,',nu,' , a spurious stable solution was found and zeroed out'
             ELSE  
                soll(j) = newsol
                fdsoll(j) = newfdsol
             END IF
          ENDIF
       ENDDO

    ELSE !launch contour solutions

       !CHOICE OF THE CENTER OF THE FIRST CONTOUR in which the kinetic solution is found.
       !The imaginary coordinate of the center is half the max fluid growth rate
       !In this manner, the upper limit of the kinetic growth rate is the (overestimated)
       !fluid growth rate. The lower limit of the "fluid growth rate" is set at 10, 
       !in case any fluid solution roots were not found


       !Find the location on the real axis for the contours. the contours will be
       !placed in multiples of rint on the real axis
       !For now this is disabled due to low values of real(C) leading to much contour overlap
       !CALL limx(C,rint)

       !This line below has led to a critical speedup. The contour overlap is optimized
       !and also the centerpoint in the imaginary plane has been slightly reduced from the 
       !fluid solution

!!$       !DEFAULT
!!$       C=C/1.5; rint=ABS(REAL(solflu(p,nu)))/5. !The divisor sets the number of contours tried per solution NEW

       C=C/1.5; rint=ABS(REAL(solflu(p,nu)))/5. !The divisor sets the number of contours tried per solution NEW

       !!FOR TESTING ONLY
       !    C=C/1.5; rint=AIMAG(C)*0.8 !The overlap factor for rint sets the degree of overlap between contours   
       !    C=C/3; rint=AIMAG(C)*1.6 !The overlap factor for rint sets the degree of overlap between contours

       !DEBUGGING*****
       !WRITE(stdout,*) 'Rint and C are:', rint, AIMAG(C)
       !**************

       !Variable initialization
       NN = 0 !number of solutions found
       soll = (0.,0.) !solution array
       fdsoll = (0.,0.) !value of fonct(soll) (should be close to zero, as soll are roots)

       ! If the mode distance is greater than the gradient lengths,
       ! then the local limit is not satisfied and we do not calculate
       IF ( d > Ro(p)/MAX( ABS( Ane(p)),ABS(Ate(p)) ) ) THEN 
          fdsoll = (-1.,-1.) !Flag that the local limit was not satisfied
       ELSE 
          L = 1.0 !begin the calculation
          ! LAUNCH VARIOUS CONTOURS SCANNING THE REAL AXIS 
          ! The loop exits when the contour center extends
          ! beyond the defined maximum on the real axis

          DO WHILE ((L+0.5)*rint < ABS(REAL(omegmax)))
             !Loop over both positive and negative frequencies in solution search

             IF ((MPI_Wtime()-calltimeinit) > timeout) THEN
                timeoutflag = .TRUE.
                EXIT
             ENDIF

             !Flags are in place in case we are certain that solutions found only in one side
             DO i = 0,2              
                IF ((MPI_Wtime()-calltimeinit) > timeout) THEN
                   timeoutflag = .TRUE.
                   EXIT
                ENDIF
                IF ( (i==0) .AND. (L>1.) ) CYCLE  !we only carry out a narrow small contour 
                IF ( ( (onlyelec .EQV. .TRUE.) .OR. (kthetarhos(nu) > 2. ) ) .AND. (i == 1) ) CYCLE  !Skips ions also for ETG scales
                IF ( (onlyion .EQV. .TRUE.) .AND. (i == 2) ) CYCLE
                !Define center of this contour

                IF (i==0) THEN !Set narrow center contour (numerically more difficult near real axis, so would rather avoid this region for a solution contour)
                   Centre = C 
                   rint= centerwidth*1.5
                ELSE
                   rint=ABS(REAL(solflu(p,nu)))/5.
                   Centre = (-1.)**i * (L * rint + centerwidth*(1-centeroverlapfac))+ C !Set contour center, with slight overlap with center contour if L=1
                ENDIF

                !DEBUGGING CODE
                !WRITE(stdout,'(2G13.5,A,I0,A,F7.3,A,F5.2)') Centre, ' i=', INT(i), ' L=',L,' Rint = ',rint

                !Solutions are now saught inside this specific contour
                !The vast bulk of QuaLiKiz computation is within this procedure

                IF ( ( rho(p) >= rhomin ) .AND. ( rho(p) <= rhomax) ) THEN !check if rho is within the defined range, otherwise return zero
                   CALL calculsol(p, nu, Centre, NN, solltmp, fdsolltmp)
                ELSE
                   solltmp(:)=0.
                   fdsolltmp(:)=0.
                   NN=0
                ENDIF
!!$                solltmp(:)=0
                !Solution cleanup: all solutions within soldel*100 percent. soldel found in datcal
                !of a previously found solution (from another contour) is set to zero             
                IF (NN > 0) THEN
                   DO j = 1,numsols
                      IF (ABS(soll(j)) > epsD) THEN !only compare to non-zero solutions in soll
                         WHERE (ABS(solltmp-soll(j))/ABS(soll(j)) < soldel)                       
                            solltmp = (0.,0.)
                            fdsolltmp = (0.,0.)
                         END WHERE
                      ENDIF
                   ENDDO
                ENDIF

                DO j=1,numsols 
                   IF ( (AIMAG(solltmp(j)) < 0. ) .OR. (ABS(AIMAG(solltmp(j))) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(solltmp(j))) > ABS(REAL(ommax(p,nu))))  ) THEN
                      IF (verbose .EQV. .TRUE.) THEN 
                         WRITE(stderr,'(A,I7,A,I2,A)') 'Solution found but outside of allowed contour range. Skipping solution. (p,nu)=(',p,',',nu,')'
                      ENDIF
                      solltmp(j) = (0.,0.)
                      fdsolltmp(j) = (0.,0.)
                   ENDIF
                ENDDO

                ! If any solutions survive, they are saved together with any previous solutions
                ! If numsols is not high enough to save all solutions, then the largest growth rates are saved first
                DO j=1,numsols 
                   IF ( ABS(solltmp(j)) > epsD ) THEN
                      minloci = MINLOC(AIMAG(soll))
                      minlocind = minloci(1)
                      IF ( (ABS(soll(minlocind)) > epsD) .AND. (verbose .EQV. .TRUE.) ) THEN
                         WRITE(stdout,'(A,I2,A,I2,A)') 'Valid instability discarded due to limited number of numsols at (p,nu)=(',p,',',nu,')'
                      ENDIF
                      IF (AIMAG(solltmp(j)) > AIMAG(soll(minlocind))) THEN
                         soll(minlocind)   = solltmp(j)
                         fdsoll(minlocind) = fdsolltmp(j)                        
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
             !Shift the contours for the next iteration
             L = L+2. - overlapfac

          ENDDO
       ENDIF
    ENDIF

    ! The solution for this (p,nu) pair is now sorted and saved

    !Reorder (descending) the solutions
    isol = AIMAG(soll)
    rsol = REAL(soll)
    ifdsol = AIMAG(fdsoll)
    rfdsol = REAL(fdsoll)

    tmpsol = AIMAG(soll)
    CALL dsort(tmpsol,isol,numsols,-2)
    tmpsol = AIMAG(soll)
    CALL dsort(tmpsol,rsol,numsols,-2)
    tmpsol = AIMAG(soll)
    CALL dsort(tmpsol,ifdsol,numsols,-2)
    tmpsol = AIMAG(soll)
    CALL dsort(tmpsol,rfdsol,numsols,-2)

    DO j=1,numsols
       sol(p,nu,j) = CMPLX(rsol(j),isol(j)) !output array
       fdsol(p,nu,j) = CMPLX(rfdsol(j),ifdsol(j)) !output array

       IF ( soll(j) /= (0.,0.) ) THEN
          !calculate integrals for quasilinear flux. Small percentage of total calculation

          IF ( (AIMAG(sol(p,nu,j)) < 0. ) .OR. (ABS(AIMAG(sol(p,nu,j))) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(sol(p,nu,j))) > ABS(REAL(ommax(p,nu))))  ) THEN
             IF (verbose .EQV. .TRUE.) THEN 
                WRITE(stderr,'(A,I7,A,I2,A)') 'Before makeQLflux: solution outside of allowed contour range (how did that happen?) Skipping QL integrals. (p,nu)=(',p,',',nu,')'
             ENDIF
             CYCLE
          ENDIF
          IF (timeoutflag .EQV. .FALSE.) THEN
             CALL make_QLflux(p,nu,sol(p,nu,j))
          ENDIF

          issol = .TRUE.
       ELSE
          issol = .FALSE.
       ENDIF
       CALL save_qlfunc(p,nu,j,issol)
       !Important IF statement for integrated modelling applications
       !If previously we had a solution, and now it's been zeroed out, we keep old solution to maintain search in the next run. This is important if we're near
       !threshold and we fluctuate between stability and instability. If we don't have this step below, then once we're stable we never recover the instability
       IF (ALLOCATED(oldsol)) THEN
          IF ( (runcounter /= 0) .AND. ( AIMAG(sol(p,nu,j)) < epsD ) .AND. ( AIMAG(oldsol(p,nu,j)) > epsD ) ) THEN
             sol(p,nu,j) = oldsol(p,nu,j) 
          ENDIF
       ENDIF
    ENDDO

    !Save growth rates to output array (normalized to nwg)
    gamma(p,nu,:) = sol(p,nu,:)

    !Check if in ETG regime. Commented out since replaced with simple ktheta limit
    !    CALL ETGcheck(ETG_flag,soll,kteta,p,nu)

  END SUBROUTINE calc

  SUBROUTINE calculsol( p, nu, Centre, NN, soll, fdsoll)  
    ! -------------------------------------------------------------------
    ! Finds the zeros of the function "fonct" within the given contour
    ! B. Davies, J. Compt. Phys. 66, 36 (1986)
    ! -------------------------------------------------------------------
    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN) :: Centre 
    INTEGER,               INTENT(OUT) :: NN
    COMPLEX(KIND=DBL), DIMENSION(numsols), INTENT(OUT) :: soll, fdsoll

    ! Local variables
    INTEGER :: i, j, k, ind
    INTEGER :: M ! Local number of contour segments
    INTEGER :: ndeg,NNbck
    LOGICAL :: scale

    REAL(KIND=DBL)    :: maxgradanglespi, maxrapfonct, minrapfonct, maxdifalphan
    REAL(KIND=DBL)    :: Nenv, thetatemp, dift, realintans, imagintans
    COMPLEX(KIND=DBL) :: omega,omega2, solli, fdsolli, nsolli, nfdsolli
    COMPLEX(KIND=DBL) :: fonx, foncttemp, omtemp, varztemp

    COMPLEX(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: fonct, fex, om, om2, varz 
    COMPLEX(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: foncti, omi
    COMPLEX(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: Sint, A, csolint
    COMPLEX(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: exprn

    REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: solint
    REAL(KIND=DBL), DIMENSION(:,:), ALLOCATABLE :: Areal
    REAL(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: thetai
    REAL(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: ww, exprnreal, exprnimag
    REAL(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: difalphan, imagrapfonct, realrapfonct
    REAL(KIND=DBL), DIMENSION(:),   ALLOCATABLE :: theta, alpha, alphan, alphaex, alphanex
    LOGICAL :: anomflag,exist
    INTEGER :: ifailloc

    soll(:)=0.
    fdsoll(:)=0.
    M=MM

    !Allocate arrays with dimension of the number of contour bits
    ALLOCATE (theta(M)) 
    ALLOCATE (varz(M)) 
    ALLOCATE (om(M)) 
    ALLOCATE (om2(M)) 
    ALLOCATE (fonct(M)) !The main player - the dispersion function!
    ALLOCATE (alpha(M)) !fonct angles
    ALLOCATE (alphan(M)) !unwrapped fonct angles
    ALLOCATE (difalphan(M-1)) !unwrapped fonct angle differences

    ! DEBUGGING
!!$    DO i=1,M
!!$       ! We follow the contour defined by the variable "om"
!!$       theta(i)    = 2.*pi*REAL(i-1)/REAL(M-1)
!!$       varz(i)     = EXP(ci*theta(i))
!!$       CALL squircle(Centre, varz(i), omega, rint)
!!$       om(i) = omega
!!$    ENDDO
!!$
!!$    IF ((p == 2) .AND. (nu == 6)) THEN
!!$       INQUIRE(file="rom.dat", exist=exist)
!!$       IF (exist) THEN
!!$          OPEN(700, file="rom.dat", status="old", position="append", action="write")
!!$       ELSE
!!$          OPEN(700, file="rom.dat", status="new", action="write")
!!$       END IF
!!$
!!$       INQUIRE(file="iom.dat", exist=exist)
!!$       IF (exist) THEN
!!$          OPEN(701, file="iom.dat", status="old", position="append", action="write")
!!$       ELSE
!!$          OPEN(701, file="iom.dat", status="new", action="write")
!!$       END IF
!!$
!!$       WRITE(700,'(17G15.7)') (REAL(om(i)),i=1,M) ; CLOSE(700)
!!$       WRITE(701,'(17G15.7)') (AIMAG(om(i)),i=1,M) ; CLOSE(701)
!!$    ENDIF
    fonct(:)=0.
    DO i=1,M
       ! We follow the contour defined by the variable "om"
       theta(i)    = 2.*pi*REAL(i-1)/REAL(M-1)
       varz(i)     = EXP(ci*theta(i))
       CALL squircle(Centre, varz(i), omega, rint)
       om(i) = omega
       ! Calculate the kinetic response at each point on the contour
       ! The root finder part of the code can be tested on a simpler case.
       ! For example:
       ! fonct(i) = (om(i)-(3.+2.*ci))*(om(i)-(-2.+4.*ci))

       IF ( (AIMAG(omega) < 0. ) .OR. (ABS(AIMAG(omega)) > ABS(2.*AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(omega)) > ABS(REAL(2.*ommax(p,nu))))  ) THEN
          IF (verbose .EQV. .TRUE.) THEN 
             WRITE(stderr,'(A,2G15.7,A,I7,A,I2,A)') 'In contours: omega outside of allowed contour range (how did that happen?) Skipping solution. Omega=,',omega,'. (p,nu)=(',p,',',nu,')'
          ENDIF
          anomflag = .TRUE.
          fonct(i) = 0. 
       ELSE
          anomflag = .FALSE.
          CALL calcfonct(p, nu, omega, fonx)
          !test timeout
          IF ((MPI_Wtime()-calltimeinit) > timeout) THEN
             timeoutflag = .TRUE.
             EXIT
          ENDIF

          fonct(i) = fonx
       ENDIF

    END DO
    DEALLOCATE (varz)

    ! The function angle is calculated at this point on the contour
    alpha(:) = ATAN2(AIMAG(fonct),REAL(fonct))
    ! The angle outputs are all placed within the same -pi +pi interval
    CALL unwrap( alpha, M, alphan) 
    ! The differences in angle are calculated
    DO i=1,(M-1)
       difalphan(i) = ABS ( alphan(i+1) - alphan(i) )
    END DO
    ! The maximum jump is found
    maxdifalphan = MAXVAL ( difalphan )
    maxgradanglespi = maxdifalphan / pi ! Maximum angle is normalized by pi

    ! Refinement of the contour segmentation. 
    ! The refinement stops once the maximum angle is less than maxangle,
    ! or if the number of segments is above maxM (both defined parameters in datcal)
    DO WHILE ( ( maxdifalphan > maxangle  ) .AND. ( M < maxM ) .AND. (anomflag .EQV. .FALSE.) )
       ! Leave loop if timed out
       IF ((MPI_Wtime()-calltimeinit) > timeout) THEN
          timeoutflag = .TRUE.
          EXIT
       ENDIF
       !Reallocate and repopulate more refined arrays
       ALLOCATE (thetai(M)) 
       ALLOCATE (omi(M)) 
       ALLOCATE (foncti(M)) 

       foncti = fonct
       omi    = om
       thetai = theta

       DEALLOCATE (theta) 
       DEALLOCATE (om) 
       DEALLOCATE (fonct) 

       ALLOCATE (theta(2*M)) 
       ALLOCATE (om(2*M)) 
       ALLOCATE (fonct(2*M)) 
       !Fill in new points in more refined array 
       DO i = 1, M
          !Find new median values of theta and the frequency on the complex plane
          thetatemp = thetai(i) + ( thetai(2) - thetai(1) ) / 2.
          varztemp  = EXP(ci * thetatemp)
          CALL squircle(Centre , varztemp, omtemp, rint)
          !New values for refined contour calculated


          IF ( (AIMAG(omtemp) < 0. ) .OR. (ABS(AIMAG(omtemp)) > ABS(2.*AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(omtemp)) > ABS(REAL(2.*ommax(p,nu))))  ) THEN
             IF (verbose .EQV. .TRUE.) THEN 
                WRITE(stderr,'(A,2G15.7,A,I7,A,I2,A)') 'In refined contours: omega outside of allowed contour range (how did that happen?) Skipping solution. Omega=,',omega,'. (p,nu)=(',p,',',nu,')'
             ENDIF
             anomflag = .TRUE.
             foncttemp = 0. 
          ELSE
             CALL calcfonct (p, nu, omtemp, foncttemp)
          ENDIF

          ind = 2*i
          !Refined arrays corresponding to the refined contour are populated
          fonct(ind-1) = foncti(i)
          fonct(ind)   = foncttemp
          theta(ind-1) = thetai(i)
          theta(ind)   = thetatemp
          om(ind-1)    = omi(i)
          om(ind)      = omtemp
       END DO

       DEALLOCATE (thetai) 
       DEALLOCATE (omi) 
       DEALLOCATE (foncti) 

       M=2*M

       !DEBUGGING
       IF (M==maxM) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A,I3,A,I7,A,I2)") 'Warning, maxM reached! M=',M,' for p=',p,' and nu=',nu
          !WRITE(stdout,*) alphan
       ENDIF

       DEALLOCATE (alpha)
       DEALLOCATE (alphan)
       DEALLOCATE (difalphan)

       ALLOCATE (alpha(M))
       ALLOCATE (alphan(M)) 
       ALLOCATE (difalphan(M-1)) 

       !Test if the new contour is sufficiently refined
       alpha(:) = ATAN2(AIMAG(fonct),REAL(fonct))
       CALL unwrap( alpha, M, alphan) 

       DO i=1,(M-1)
          difalphan(i)    = ABS ( alphan(i+1) - alphan(i) )
       END DO

       !DEBUGGING
       !IF (M==maxM) THEN
       !OPEN(unit=700, file="difalphan.dat", action="write", status="replace")
       !WRITE(700,'(G15.7)') (difalphan(i),i=1,M-1) ; CLOSE(700)
       !OPEN(unit=700, file="rom.dat", action="write", status="replace")
       !WRITE(700,'(G15.7)') (REAL(om(i)),i=1,M) ; CLOSE(700)
       !OPEN(unit=700, file="iom.dat", action="write", status="replace")
       !WRITE(700,'(G15.7)') (AIMAG(om(i)),i=1,M) ; CLOSE(700)
       !STOP
       !ENDIF

       !The new maximum angle jump is calculated
       maxdifalphan = MAXVAL ( difalphan )
       maxgradanglespi = maxdifalphan / pi

       Nenv = ABS( (alphan(M) - alphan(1) ) / 2. / pi )
       !test timeout
       IF ((MPI_Wtime()-calltimeinit) > timeout) THEN
          anomflag = .TRUE.
       ENDIF

       !DEBUG CODE
       !WRITE(stdout,*) 'M=',M,'and Nenv=',Nenv

       !In most cases, if Nenv isn't near one after the first iteration then
       !the contour is going over a problematic region where the phase jumps back and forth
       !This solution is then abandoned. This choice saves typically ~20-30% computation time
       IF (Nenv < 0.5) anomflag=.TRUE.

    END DO
    DEALLOCATE (difalphan) 

    !Count the number of zeros found in the contour
    Nenv = ABS( (alphan(M) - alphan(1) ) / 2. / pi )

    !Abandon solution of D(omega)=0 anywhere on contour
    IF (ANY(ABS(fonct)<epsD)) THEN
       WRITE(stderr,'(A,I2,A,I2,A)') 'Main contour phase: contour had D(omega)=0! Skipping solution. (p,nu)=(',p,',',nu,')'
       anomflag = .TRUE.
    ENDIF

    IF (anomflag .EQV. .TRUE.) Nenv = 0.
    NN = NINT( Nenv ) 

    !DEBUG CODE********************
!!$    WRITE(stdout,*) 'For p=',p,'and nu= ',nu,', final M = ', M
!!$    WRITE(stdout,*) 'Nenv =', Nenv
    !******************************

    !If, despite refining the segments, the angle jumps are too high, then
    !and also the phase is not close to an integer, then the contour was problematic
    !and the exact solution is not sought out
    IF ( (maxgradanglespi > maxangle) .AND. ( ABS(Nenv-NN) > 0.01 ) )  THEN 
       Nsolrat = Nsolrat+NN; NN = 0
    END IF

    IF (NN == 0) THEN
       ! No solutions sought for!
    ELSE
       ! There are solutions and they will be found with the Davies method
       ALLOCATE ( fex(M) )
       ALLOCATE ( exprn(M,NN) )

       !Modified dispersion relation function. This is done to remove the branch cut when taking
       !the logarithm of fex. The values of the Sn 'argument principle' integrals are unchanged
       !since the imposed roots are at z=0 (see Davies 1986)
       DO i=1,M
          fex(i) = fonct(i) * EXP( -ci * NN * theta(i) ) 
       END DO

       ALLOCATE ( alphaex(M) )
       ALLOCATE ( alphanex(M) )

       !Find argument of modified dispersion relation
       alphaex(:) = ATAN2(AIMAG(fex),REAL(fex))

       !Unwrap the arguments to be in ascending order in pi
       CALL unwrap( alphaex, M, alphanex)

       !Integrand defined for 'argument principle' integrals (equation 3.3 in Davies 1986)
       DO i=1,M       
          DO j=1,NN
             exprn(i,j) = -REAL(j)/(2.*pi) * EXP(ci*REAL(j)*theta(i)) &
                  *  ( LOG ( ABS ( fex(i) ) ) + ci*alphanex(i) ) 
          END DO
       END DO

       DEALLOCATE ( fex )
       DEALLOCATE ( alphanex )

       !Argument principle integrals (Sn) are calculated
       ALLOCATE ( Sint(NN) )
       ALLOCATE ( exprnreal(M) )
       ALLOCATE ( exprnimag(M) )
       DO j=1,NN
          exprnreal = REAL(exprn(:,j))
          exprnimag = AIMAG(exprn(:,j))
          CALL davint(theta,exprnreal, M, 0._DBL, 2.*pi , realintans, ifailloc,1)
          CALL davint(theta,exprnimag, M, 0._DBL, 2.*pi , imagintans, ifailloc,2)
          Sint(j) = CMPLX(realintans,imagintans)
       END DO

       DEALLOCATE ( exprnreal )
       DEALLOCATE ( exprnimag )
       DEALLOCATE ( theta )       
       DEALLOCATE ( exprn )

       ALLOCATE ( A(0:NN) )
       !Setup the system of equations for the Davies method
       !This works as follows. From the argument principle we have:
       ! Sint(1) is sum(z) (sum of all solutions)
       ! Sint(2) is sum(z^2) (sum of all squares of solutions) 
       ! Sint(3) is sum(z^3) (sum of all cubes of solutions) ... and so on
       ! Based on this information, a polynomial is constructed with the form below,
       ! which can be proven that the roots of the polynomial are the solutions
       A(0) = (1.,0.)
       DO j=1,NN
          A(j) = 0.
          DO k=0,(j-1)
             ind = j-k
             A(j) = A(j) + A(k) * Sint(ind)
          END DO
          A(j) = -A(j) / j
       END DO
       DEALLOCATE ( Sint )

       !Construct arrays for root finding in polynomial
       ndeg = NN
       ALLOCATE ( Areal(2,1:NN+1) )
       ALLOCATE ( solint(2,NN) )
       ALLOCATE ( ww(2*NN*(NN+1)) )

       DO j = 0,NN
          Areal(1,j+1) = REAL( A(j) )
          Areal(2,j+1) = AIMAG( A(j) )
       END DO

       DEALLOCATE ( A )
       !The roots of the polynomials are found, which are the solutions!
       CALL CPQR79 (ndeg, Areal, solint, ifailloc, ww)

       DEALLOCATE ( Areal )       
       ALLOCATE ( csolint(NN) )

       NNbck=NN
       DO j=1,NN
          !Each solution is rewritten into a complex variable
          csolint(j) = CMPLX( solint(1,j), solint(2,j) )
          !We now verify that a true root was indeed found
          !CALL contour(Centre, csolint(j), solli)

          CALL squircle(Centre, csolint(j), solli, rint)

          !DEBUG
          !WRITE(*,'(A,I3,I3,A,2G15.7)') 'p/nu/j=',p,nu,'. solli = ',solli

          IF ( (gkw_is_nan(AIMAG(solli))) .OR. (gkw_is_nan(REAL(solli)))) THEN
             IF (verbose .EQV. .TRUE.) THEN 
                WRITE(stderr,'(A,I7,A,I2,A)') 'Main contour phase: solution had a NaN! Skipping solution. (p,nu)=(',p,',',nu,')'
             ENDIF
             soll(j)   = (0.,0.)
             fdsoll(j)   = (0.,0.)
             CYCLE
          ENDIF

          IF ( (AIMAG(solli) < 0. )) THEN 
             IF (verbose .EQV. .TRUE.) THEN 
                WRITE(stderr,'(A,I7,A,I2,A)') 'Main contour phase: solution is negative growth rate (stable eigenmode). Skipping solution. (p,nu)=(',p,',',nu,')'
             ENDIF
             soll(j)   = (0.,0.)
             fdsoll(j)   = (0.,0.)
             CYCLE
          ENDIF

          IF ( (ABS(AIMAG(solli)) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(solli)) > ABS(REAL(ommax(p,nu))))  ) THEN
             IF (verbose .EQV. .TRUE.) THEN 
                WRITE(stderr,'(A,I7,A,I2,A)') 'Main contour phase: solution out of allowed contour range (how did that happen?) Skipping solution. (p,nu)=(',p,',',nu,')'
             ENDIF
             soll(j)   = (0.,0.)
             fdsoll(j)   = (0.,0.)
             CYCLE
          ENDIF

          CALL calcfonct(p, nu, solli, fdsolli)

          !DEBUGGING SKIPPING THE NEWTON STAGE***
          ! soll(j) = solli
          ! fdsoll(j) = fdsolli
          ! EXIT

          !If in fact our 'solution' is nowhere close to a root, then
          !we skip the Newton refinement stage and abandon the solution
          IF (ABS(fdsolli) > nearlysol) THEN
             soll(j)   = (0.,0.)
             fdsoll(j)   = (0.,0.)
             IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I3,A,G15.7,A)') 'For p=',p,' nu=,',nu,' , ABS(fdsolli)=',ABS(fdsolli),' and the solution was abandoned'
             NNbck=NNbck-1
          ELSE !We have a close solution, now refined with a Newton step
             !DEGBUGGING 
             !WRITE(stdout,*) 'BEFORE NEWTON', solli
             !WRITE(stdout,*) 'NUMBER OF CONTOUR SEGMENTS', M
!!!!!
             !CALL newton(p, nu, solli, fdsolli, nsolli, nfdsolli)
             CALL broyden(p, nu, solli, fdsolli, nsolli, nfdsolli)

             !WRITE(stdout,*) 'AFTER NEWTON', nsolli

             !If, even after the Newton refinement, the solution is still too coarse, or outside 
             !the preset boundaries, then we abandon the solution
             IF (ABS(nfdsolli)> 0.5  .OR. &
                  ABS(AIMAG(nsolli))>(1.5*ABS(AIMAG(omegmax))) .OR. &
                  ABS(REAL(nsolli)) > (2.*ABS(REAL(omegmax)))) THEN
                soll(j)   = (0.,0.)
                fdsoll(j) = (0.,0.)
                IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I3,A,G15.7,A,2G15.7,A)') 'For p=',p,' nu=,',nu,' , ABS(nfdsolli)=',ABS(nfdsolli),' with nsolli=',nsolli,' and the solution was abandoned'
                NNbck=NNbck-1
                !If the solution is stable, then the solution is zeroed out. We can't have stable solutions since we haven't included the analytic continuation
             ELSEIF (AIMAG(nsolli)<0.) THEN
                soll(j)   = (0.,0.)
                fdsoll(j) = (0.,0.)
                IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I3,A)') 'For p=',p,' nu=,',nu,' , a spurious stable solution was found and zeroed out'
                NNbck=NNbck-1
             ELSE  
                soll(j) = nsolli
                fdsoll(j) = nfdsolli
             END IF
          END IF
       END DO
       NN=NNbck

       DEALLOCATE ( csolint )
       DEALLOCATE ( solint )
       DEALLOCATE ( ww )
       DEALLOCATE ( alphaex )
       DEALLOCATE ( fonct )
       DEALLOCATE ( alphan )
       DEALLOCATE ( alpha )
       DEALLOCATE ( om )
       DEALLOCATE ( om2 )

    END IF


    IF (ALLOCATED(theta)) DEALLOCATE(theta)
    IF (ALLOCATED(om)) DEALLOCATE(om)
    IF (ALLOCATED(om2)) DEALLOCATE(om2)
    IF (ALLOCATED(fonct)) DEALLOCATE(fonct)
    IF (ALLOCATED(fonct)) DEALLOCATE(fonct)
    IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
    IF (ALLOCATED(alphan)) DEALLOCATE(alphan)
    IF (ALLOCATED(difalphan)) DEALLOCATE(difalphan)


  END SUBROUTINE calculsol

  SUBROUTINE calcfonct(p, nu, omega, fonx)
    ! -------------------------------------------------------------------
    ! Calculates the integral function for which we search the roots
    ! This function is comprised of the adiabatic (Ac), passing (fonctc)
    ! and trapped (fonctp) terms
    ! -------------------------------------------------------------------
    INTEGER, INTENT(in)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonx

    COMPLEX(KIND=DBL) :: fonctc
    COMPLEX(KIND=DBL) :: fonctp
        
    IF(int_method.EQ.0) THEN !Use NAG

      IF ( ( rotflagarray(p) == 1 ) .AND. ( ETG_flag(nu) .EQV. .FALSE. ) ) THEN
         ! replace mwidth by real(mwidth) in such comparaisons since now mwidth is complex. Warning: is it correct or should take module?
         IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
         ELSE
            CALL calcfonctrotc ( p, nu, omega, fonctc ) 
         END IF

         IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
         ELSE      
            CALL calcfonctrotp ( p, nu, omega, fonctp )
         END IF

         fonx = CMPLX(Ac(p),0.) - fonctc - fonctp

      ELSE

         IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
         ELSE
            CALL calcfonctc ( p, nu, omega, fonctc ) 
         END IF

         IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
         ELSE      
            CALL calcfonctp ( p, nu, omega, fonctp )
         END IF

         fonx = CMPLX(Ac(p),0.) - fonctc - fonctp

      ENDIF
    
    ELSE IF(int_method.EQ.1) THEN !Use hcubature
      IF((rotflagarray(p) == 1).AND.(ETG_flag(nu).EQV. .FALSE.)) THEN !rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonctrot_hcubature(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_hcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_hcubaturec(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_hcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_hcubaturec2(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      ELSE !no rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonct_hcubature(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_hcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_hcubaturec(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_hcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_hcubaturec2(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      END IF
    ELSE IF(int_method.EQ.2) THEN !Use pcubature
      IF ( ( rotflagarray(p) == 1 ) .AND. ( ETG_flag(nu) .EQV. .FALSE. ) ) THEN !rotations
    
        IF(int_split.EQ.0) THEN
          CALL calcfonctrot_pcubature(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          ! replace mwidth by real(mwidth) in such comparaisons since now mwidth is complex. Warning: is it correct or should take module?
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
          ELSE
            CALL calcfonctrot_pcubaturec ( p, nu, omega, fonctc ) 
          END IF

          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE      
            CALL calcfonctrot_pcubaturep ( p, nu, omega, fonctp )
          END IF

          fonx = CMPLX(Ac(p),0.) - fonctc - fonctp
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_pcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_pcubaturec2(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      ELSE !no rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonct_pcubature(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
          ELSE
            CALL calcfonct_pcubaturec ( p, nu, omega, fonctc ) 
          END IF

          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE      
            CALL calcfonct_pcubaturep ( p, nu, omega, fonctp )
          END IF

          fonx = CMPLX(Ac(p),0.) - fonctc - fonctp
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_pcubaturep(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_pcubaturec2(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      END IF
    END IF

  END SUBROUTINE calcfonct
  
  SUBROUTINE calcfonct_newt(p, nu, omega, fonx)
    ! -------------------------------------------------------------------
    ! Calculates the integral function for which we search the roots
    ! This function is comprised of the adiabatic (Ac), passing (fonctc)
    ! and trapped (fonctp) terms
    ! Used only for the Newton method
    ! -------------------------------------------------------------------
    INTEGER, INTENT(in)  :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN)  :: omega
    COMPLEX(KIND=DBL), INTENT(OUT) :: fonx

    COMPLEX(KIND=DBL) :: fonctc
    COMPLEX(KIND=DBL) :: fonctp
    
    
    

    IF(newt_method.EQ.0) THEN !Use NAG

      IF ( ( rotflagarray(p) == 1 ) .AND. ( ETG_flag(nu) .EQV. .FALSE. ) ) THEN
         ! replace mwidth by real(mwidth) in such comparaisons since now mwidth is complex. Warning: is it correct or should take module?
         IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
         ELSE
            CALL calcfonctrotc ( p, nu, omega, fonctc ) 
         END IF

         IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
         ELSE      
            CALL calcfonctrotp ( p, nu, omega, fonctp )
         END IF

         fonx = CMPLX(Ac(p),0.) - fonctc - fonctp

      ELSE

         IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
            IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'Warning: REAL(mwidth)<d/4 for p/nu=',p,'/',nu
         ELSE
            CALL calcfonctc ( p, nu, omega, fonctc ) 
         END IF

         IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
         ELSE      
            CALL calcfonctp ( p, nu, omega, fonctp )
         END IF

         fonx = CMPLX(Ac(p),0.) - fonctc - fonctp

      ENDIF
    
    
    ELSE IF(newt_method.EQ.1) THEN !Use hcubature
      IF((rotflagarray(p) == 1).AND.(ETG_flag(nu).EQV. .FALSE.)) THEN !rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonctrot_hcubature(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_hcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_hcubaturec_newt(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_hcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_hcubaturec2_newt(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      ELSE !no rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonct_hcubature_newt(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_hcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_hcubaturec_newt(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_hcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_hcubaturec2_newt(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      END IF
    ELSE IF(newt_method.EQ.2) THEN !Use pcubature
      IF((rotflagarray(p) == 1).AND.(ETG_flag(nu).EQV. .FALSE.)) THEN !rotations
    
        IF(int_split.EQ.0) THEN
          CALL calcfonctrot_pcubature_newt(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_pcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_pcubaturec_newt(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonctrot_pcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonctrot_pcubaturec2_newt(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      ELSE !no rotations
        IF(int_split.EQ.0) THEN
          CALL calcfonct_pcubature_newt(p, nu, omega, fonx)
          
        ELSE IF(int_split.EQ.1) THEN
        
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_pcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_pcubaturec_newt(p, nu, omega, fonctc)
          END IF
          
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
          
        ELSE IF(int_split.EQ.2) THEN
          IF (ft(p)==0. .OR. (calctrap .EQV. .FALSE.) ) THEN
            fonctp = 0.    
          ELSE    
            CALL calcfonct_pcubaturep_newt(p, nu, omega, fonctp)
          END IF
          
          IF ( ( fc(p)==0. ) .OR. ( REAL(mwidth)<d/4.) .OR. ( calccirc .EQV. .FALSE. ) ) THEN
            fonctc = 0.
          ELSE
            CALL calcfonct_pcubaturec2_newt(p, nu, omega, fonctc)
          END IF
          fonx = CMPLX(Ac(p), 0.) - fonctc - fonctp 
        END IF
      END IF
    END IF

  END SUBROUTINE calcfonct_newt

  SUBROUTINE newton( p, nu, sol, fsol, newsol, fnewsol)  
    !Newton method for refining the solutions found from the contour integrals
    !Basic 2D Newton method. Demands that both u(z) and v(z) go to zero, where F(z)=u(z)+i*v(z)
    !Function derivative defined with ndif parmameter
    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN) :: sol, fsol

    COMPLEX(KIND=DBL), INTENT(OUT) :: newsol, fnewsol

    INTEGER :: niter
    REAL(KIND=DBL) :: err, detMa, maxnerrint
    REAL(KIND=DBL), DIMENSION(2,2) :: Ma, invMa
    REAL(KIND=DBL), DIMENSION(2,1) :: uo, deltau
    COMPLEX(KIND=DBL) :: fzo, fzpd, fzmd, fzpid, fzmid
    COMPLEX(KIND=DBL) :: dfsurdx, dfsurdy
    COMPLEX(KIND=DBL) :: zo, zpd, zmd, zpid, zmid, delom, delomi 

    delom  = (ndif,0.)
    delomi = (0.,ndif)
    niter  = 0
    zo     = sol
    fzo    = fsol
    IF(phys_meth.EQ.0) THEN
      err    = ABS(fzo)
    ELSE 
      err = 1._DBL
    END IF

    IF (runcounter /= 0) THEN
      maxnerrint = maxnerr
    ELSE
      maxnerrint = maxnerr2
    ENDIF

    DO
    IF ( (err < maxnerrint) .OR. (niter > maxiter) ) EXIT

      zpd = zo+delom
      zmd = zo-delom

      CALL calcfonct_newt (p, nu, zpd, fzpd)
      CALL calcfonct_newt (p, nu, zmd, fzmd)

      dfsurdx = (fzpd-fzmd)/(2.*delom)

      zpid=zo+delomi
      zmid=zo-delomi

      CALL calcfonct_newt (p, nu, zpid, fzpid)
      CALL calcfonct_newt (p, nu, zmid, fzmid)

      dfsurdy = (fzpid-fzmid)/(2.*delom)

      Ma(1,1) = REAL(dfsurdx)
      Ma(1,2) = REAL(dfsurdy)
      Ma(2,1) = AIMAG(dfsurdx)
      Ma(2,2) = AIMAG(dfsurdy)

      uo(1,1) = -REAL(fzo)
      uo(2,1) = -AIMAG(fzo)
      detMa = Ma(1,1)*Ma(2,2)-Ma(1,2)*Ma(2,1)
      !FOR DEBUGGING UNCOMMENT BELOW
      !WRITE(stderr,*) detMa, niter, err
      IF (ABS(detMa)<epsD) EXIT
      invMa(1,1) =  Ma(2,2)/detMa
      invMa(1,2) = -Ma(1,2)/detMa
      invMa(2,1) = -Ma(2,1)/detMa
      invMa(2,2) =  Ma(1,1)/detMa

      deltau = MATMUL(invMa,uo)

      zo = zo+CMPLX(deltau(1,1),deltau(2,1))

      IF ( (gkw_is_nan(AIMAG(zo))) .OR. (gkw_is_nan(REAL(zo)))) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I2,A)') 'Newton: solution had a NaN! Skipping solution. (p,nu)=(',p,',',nu,')'
        zo   = (0.,0.)
        fzo  = (0.,0.)
        EXIT
      ENDIF

      IF ( (AIMAG(zo) < 0. ) .OR. (ABS(AIMAG(zo)) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(zo)) > ABS(REAL(ommax(p,nu))))  ) THEN 
        IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I2,A)') 'Newton: solution outside of allowed contour range! Skipping solution. (p,nu)=(',p,',',nu,')'
        zo   = (0.,0.)
        fzo  = (0.,0.)
        EXIT
      ENDIF

      CALL calcfonct_newt (p, nu, zo, fzo)

      !Change Newton convergence test
      IF(newt_conv.EQ.0) THEN
        err = ABS(fzo)
      ELSE IF(newt_conv.EQ.1) THEN
        err = ABS(CMPLX(deltau(1,1), deltau(2,1)))
      ELSE IF(newt_conv.EQ.2) THEN
        err = ABS(CMPLX(deltau(1,1), deltau(2,1))/zo)
      END IF

      niter = niter+1
    END DO
    !FOR DEBUGGING*************
    !WRITE(stdout,*) 'Number of iterations in Newton solver = ', niter
    !**************************
    !WRITE(stderr, *) "In Newton for p = ", p, " and nu = ", nu, ", niter = ", niter
    newsol  = zo
    fnewsol = fzo   

  END SUBROUTINE newton

  SUBROUTINE broyden( p, nu, sol, fsol, newsol, fnewsol)
    !Broyden method for refining the solutions found from the contour integrals
    !Basic 2D Broyden method. Demands that both u(z) and v(z) go to zero, where F(z)=u(z)+i*v(z)
    !Function derivative defined with ndif parmameter
    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL), INTENT(IN) :: sol, fsol

    COMPLEX(KIND=DBL), INTENT(OUT) :: newsol, fnewsol

    INTEGER :: niter, i
    REAL(KIND=DBL) :: err, detMa, maxnerrint
    REAL(KIND=DBL), DIMENSION(2,2) :: Ma, invMa, I0
    REAL(KIND=DBL), DIMENSION(2,1) :: uo, deltau, deltaf
    REAL(KIND=DBL), DIMENSION(1,1) :: denom
    COMPLEX(KIND=DBL) :: fzo, fz_old, fzpd, fzmd, fzpid, fzmid
    COMPLEX(KIND=DBL) :: dfsurdx, dfsurdy
    COMPLEX(KIND=DBL) :: zo, zpd, zmd, zpid, zmid, delom, delomi 

    delom  = (ndif,0.)
    delomi = (0.,ndif)
    niter  = 0
    zo     = sol
    fzo    = fsol
    
    !Idenitty matrix
    I0(:,:) = 0.
    forall(i=1:2) I0(i,i) = 1.
    
    IF(phys_meth.EQ.0) THEN
      err    = ABS(fzo)
    ELSE 
      err = 1._DBL
    END IF

    IF (runcounter /= 0) THEN
      maxnerrint = maxnerr
    ELSE
      maxnerrint = maxnerr2
    ENDIF

    DO
    IF ( (err < maxnerrint) .OR. (niter > maxiter) ) EXIT
      uo(1,1) = -REAL(fzo)
      uo(2,1) = -AIMAG(fzo)
      deltaf(1,1) = REAL(fzo) - REAL(fz_old)
      deltaf(2,1) = AIMAG(fzo) - AIMAG(fz_old)
      
      IF(niter.EQ.0) THEN
        ! Only calculate the Jacobian once
        ! Afterwards, we approximately update it
        zpd = zo+delom
        zmd = zo-delom

        CALL calcfonct_newt (p, nu, zpd, fzpd)
        CALL calcfonct_newt (p, nu, zmd, fzmd)

        dfsurdx = (fzpd-fzmd)/(2.*delom)

        zpid=zo+delomi
        zmid=zo-delomi

        CALL calcfonct_newt (p, nu, zpid, fzpid)
        CALL calcfonct_newt (p, nu, zmid, fzmid)

        dfsurdy = (fzpid-fzmid)/(2.*delom)

        Ma(1,1) = REAL(dfsurdx)
        Ma(1,2) = REAL(dfsurdy)
        Ma(2,1) = AIMAG(dfsurdx)
        Ma(2,2) = AIMAG(dfsurdy)
        detMa = Ma(1,1)*Ma(2,2)-Ma(1,2)*Ma(2,1)
        !FOR DEBUGGING UNCOMMENT BELOW
        !WRITE(stderr,*) detMa, niter, err
        IF (ABS(detMa)<epsD) EXIT
        invMa(1,1) =  Ma(2,2)/detMa
        invMa(1,2) = -Ma(1,2)/detMa
        invMa(2,1) = -Ma(2,1)/detMa
        invMa(2,2) =  Ma(1,1)/detMa
      ELSE
        !Update invMa
        !This minimizes the Frobenius norm of Ma_new - Ma_old, i.e. the "good" update
        denom = MATMUL(TRANSPOSE(deltau), MATMUL(invMa, deltaf))
        IF(ABS(denom(1,1))<epsD) EXIT
        invMa = MATMUL(I0 + 1./denom(1,1) * MATMUL(deltau - MATMUL(invMa, deltaf), TRANSPOSE(deltau)), invMa)
      END IF
      

      deltau = MATMUL(invMa,uo)

      zo = zo+CMPLX(deltau(1,1),deltau(2,1))

      IF ( (gkw_is_nan(AIMAG(zo))) .OR. (gkw_is_nan(REAL(zo)))) THEN
        IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I2,A)') 'Broyden: solution had a NaN! Skipping solution. (p,nu)=(',p,',',nu,')'
        zo   = (0.,0.)
        fzo  = (0.,0.)
        EXIT
      ENDIF

      IF ( (AIMAG(zo) < 0. ) .OR. (ABS(AIMAG(zo)) > ABS(AIMAG(ommax(p,nu)))) .OR. (ABS(REAL(zo)) > ABS(REAL(ommax(p,nu))))  ) THEN 
        IF (verbose .EQV. .TRUE.) WRITE(stderr,'(A,I7,A,I2,A)') 'Broyden: solution outside of allowed contour range! Skipping solution. (p,nu)=(',p,',',nu,')'
        zo   = (0.,0.)
        fzo  = (0.,0.)
        EXIT
      ENDIF

      fz_old = fzo
      CALL calcfonct_newt (p, nu, zo, fzo)

      !Change Broyden convergence test
      IF(newt_conv.EQ.0) THEN
        err = ABS(fzo)
      ELSE IF(newt_conv.EQ.1) THEN
        err = ABS(CMPLX(deltau(1,1), deltau(2,1)))
      ELSE IF(newt_conv.EQ.2) THEN
        err = ABS(CMPLX(deltau(1,1), deltau(2,1))/zo)
      END IF

      niter = niter+1
    END DO
    !FOR DEBUGGING*************
    !WRITE(stdout,*) 'Number of iterations in Broyden solver = ', niter
    !**************************
    !WRITE(stderr, *) "In Broyden for p = ", p, " and nu = ", nu, ", niter = ", niter
    newsol  = zo
    fnewsol = fzo   

  END SUBROUTINE
  
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
  
  INTEGER FUNCTION alam1int_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = alam1int(xx)
    output = 0
  END FUNCTION alam1int_cubature

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
  
  INTEGER FUNCTION alam2int_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = alam2int(xx)
    output = 0
  END FUNCTION alam2int_cubature

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
  
  INTEGER FUNCTION alam3int_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = alam3int(xx)
    output = 0
  END FUNCTION alam3int_cubature

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
  
  INTEGER FUNCTION alam4int_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = alam4int(xx)
    output = 0
  END FUNCTION alam4int_cubature

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
  
  INTEGER FUNCTION alam5int_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = alam5int(xx)
    output = 0
  END FUNCTION alam5int_cubature


  REAL(KIND=DBL) FUNCTION sin2thint(krr)
!!! Integrand for average of sin^2(theta/2) over eigenfunction used in V_par averages
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs,gau

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    sin2thint= normgs * ABS(gau)**2. * SIN(distan(plam,nulam)*krr*normkr/2)**2

  END FUNCTION sin2thint

  INTEGER FUNCTION sin2thint_cubature(ndim, x, fdata, fdim, fval) RESULT(output)
    USE KIND
    INTEGER, INTENT(IN) :: ndim, fdim
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
    REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
    REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    
    REAL(KIND=DBL) :: xx
    
    IF((ndim.NE.1).OR.(fdim.NE.1)) THEN
      output = 1
      RETURN
    END IF
    
    xx = x(1) 
    
    fval(1) = sin2thint(xx)
    output = 0
  END FUNCTION sin2thint_cubature



END MODULE calcroutines
