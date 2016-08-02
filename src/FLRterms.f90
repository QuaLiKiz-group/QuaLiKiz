MODULE FLRterms
  USE kind
  USE datmat
  USE datcal

  IMPLICIT NONE

CONTAINS

  SUBROUTINE makeFLRterms(p,nu)
    !! The scaled Bessel functions are integrated over kr separately, outside the main integrals
    !! Aribitrary normalisation factor (normkr) introduced to optimize kr integration
    INTEGER, INTENT(IN) :: p, nu
    REAL(kind=DBL) :: minFLR, maxFLR, relerr
    INTEGER :: npts !output of number of integral evaluations
    INTEGER :: ifailloc

    maxFLR =   ABS(2*pi/(REAL(mwidth)*normkr))
    minFLR = - maxFLR

    !Allocate work arrays for integration routine
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    !Trapped electrons
!!$    CALL DQAGSE(nFLRep,minFLR,maxFLR,epsFLR,epsFLR,limit,Joe2p,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)
!!$
!!$    CALL DQAGSE(nFLRep1,minFLR,maxFLR,epsFLR,epsFLR,limit,J1e2p,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)
!!$
!!$    !Passing electrons
!!$    CALL DQAGSE(nFLRec,minFLR,maxFLR,epsFLR,epsFLR,limit,Joe2c,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)

    ifailloc = 1
    Joe2p = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRep,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of J0e2p FLR integration at p=',p,', nu=',nu
    ENDIF

    ifailloc = 1
    J1e2p = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRep1,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of J1e2p FLR integration at p=',p,', nu=',nu
    ENDIF

    ifailloc = 1
    Joe2c = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRec,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of J0e2c FLR integration at p=',p,', nu=',nu
    ENDIF


    DO ion = 1,nions
       !       CALL DQAGSE(nFLRip,minFLR,maxFLR,epsFLR,epsFLR,limit,Joi2p(ion),relerr,npts,ifailloc,&
       !           alist, blist, rlist, elist, iord, last)

       ifailloc = 1
       Joi2p(ion) = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRip,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of J0i2p FLR integration at p=',p,', nu=',nu
       ENDIF

       ifailloc = 1
       J1i2p(ion) = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRip1,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of J1i2p FLR integration at p=',p,', nu=',nu
       ENDIF

       !       CALL DQAGSE(nFLRip1,minFLR,maxFLR,epsFLR,epsFLR,limit,J1i2p(ion),relerr,npts,ifailloc,&
       !           alist, blist, rlist, elist, iord, last)

       !Passing ions 
       !CALL DQAGSE(nFLRic,minFLR,maxFLR,epsFLR,epsFLR,limit,Joi2c(ion),relerr,npts,ifailloc,&
       !     alist, blist, rlist, elist, iord, last)

       ifailloc = 1
       Joi2c(ion) = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRic,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of J0i2c FLR integration at p=',p,', nu=',nu
       ENDIF

    ENDDO

    !Deallocate work arrays
    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)


  END SUBROUTINE makeFLRterms

  SUBROUTINE makeFLRtermsrot(p,nu)
    !! The scaled Bessel functions are integrated over kr separately, outside the main integrals
    !! Aribitrary normalisation factor (normkr) introduced to optimize kr integration
    INTEGER, INTENT(IN) :: p, nu
    REAL(kind=DBL) :: minFLR, maxFLR, relerr
    INTEGER :: npts !output of number of integral evaluations
    INTEGER :: ifailloc

    maxFLR =   ABS(2*pi/(widthhat*normkr))
    minFLR = - maxFLR

    !Allocate work arrays for integration routine
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    !Trapped electrons
!!$    CALL DQAGSE(nFLReprot,minFLR,maxFLR,epsFLR,epsFLR,limit,Joe2p,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)
!!$
!!$    CALL DQAGSE(nFLRep1rot,minFLR,maxFLR,epsFLR,epsFLR,limit,J1e2p,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)
    ifailloc = 1
    Joe2p = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLReprot,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of J0e2p FLR integration at p=',p,', nu=',nu
    ENDIF
    ifailloc = 1
    J1e2p = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRep1rot,lw,ifailloc)
    IF (ifailloc /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
            &'. Abnormal termination of J0e2p FLR integration at p=',p,', nu=',nu
    ENDIF


    !Passing electrons
!!$    CALL DQAGSE(nFLRecrot,minFLR,maxFLR,epsFLR,epsFLR,limit,Joe2c,relerr,npts,ifailloc,&
!!$         alist, blist, rlist, elist, iord, last)

    DO ion = 1,nions
!!$       CALL DQAGSE(nFLRiprot,minFLR,maxFLR,epsFLR,epsFLR,limit,Joi2p(ion),relerr,npts,ifailloc,&
!!$            alist, blist, rlist, elist, iord, last)
!!$
!!$       CALL DQAGSE(nFLRip1rot,minFLR,maxFLR,epsFLR,epsFLR,limit,J1i2p(ion),relerr,npts,ifailloc,&
!!$            alist, blist, rlist, elist, iord, last)
       ifailloc = 1
       Joi2p(ion) = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRiprot,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of J0i2p FLR integration at p=',p,', nu=',nu
       ENDIF

       ifailloc = 1                                
       J1i2p(ion) = d01ahf(minFLR,maxFLR,epsFLR,npts,relerr,nFLRip1rot,lw,ifailloc)
       IF (ifailloc /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I0,A,I0)") 'ifailloc = ',ifailloc,&
               &'. Abnormal termination of J1i2p FLR integration at p=',p,', nu=',nu
       ENDIF

       !Passing ions 
!!$       CALL DQAGSE(nFLRicrot,minFLR,maxFLR,epsFLR,epsFLR,limit,Joi2c(ion),relerr,npts,ifailloc,&
!!$            alist, blist, rlist, elist, iord, last)

    ENDDO

    !Deallocate work arrays
    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)


  END SUBROUTINE makeFLRtermsrot


  REAL(KIND=DBL) FUNCTION nFLRec(krr)
!!! Routine for improved FLR, passing electrons
!!! FLR of passing particles are calculated from an integration over kr of
!!! Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var2, var3
    REAL(KIND=DBL)    :: bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var2 = (normkr*krr*Rhoe(pnFLR))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoe)**2.          !!2nd argument of Bessel fun

    bessm2 = BESEI0(var2+var3)
    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun

    nFLRec = normgs * gau**2. * bessm2 

  END FUNCTION nFLRec

  REAL(KIND=DBL) FUNCTION nFLRep(krr)
!!! Routine for improved FLR, trapped electrons
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) * Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var1 = (normkr*krr*de(pnFLR))**2.   !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoe(pnFLR))**2. !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoe)**2.               !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI0(var1)
    bessm2 = BESEI0(var2+var3)
    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun

    nFLRep = normgs * gau**2. * bessm1 * bessm2 

  END FUNCTION nFLRep

  REAL(KIND=DBL) FUNCTION nFLRep1(krr)
!!! Routine for improved  trapped electrons response evaluation
!!! Pierre 22/2/11
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var1 = (normkr*krr*de(pnFLR))**2.    !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoe(pnFLR))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoe)**2.               !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI1(var1)
    bessm2 = BESEI0(var2+var3)

    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun

    nFLRep1 = normgs * gau**2. * bessm1*bessm2 

  END FUNCTION nFLRep1

  REAL(KIND=DBL) FUNCTION nFLRic(krr)
!!! Routine for improved FLR, passing ions
!!! FLR of passing particles are calculated from an integration over kr of
!!! Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var2, var3
    REAL(KIND=DBL)    :: bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoi(ion))**2.               !!2nd argument of Bessel fun

    bessm2 = BESEI0(var2+var3)
    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun

    nFLRic = normgs * gau**2. * bessm2 

  END FUNCTION nFLRic

  REAL(KIND=DBL) FUNCTION nFLRip(krr)
!!! Routine for improved FLR, trapped ions
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) * Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var1 = (normkr*krr*di(pnFLR,ion))**2.   !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoi(ion))**2.             !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI0(var1)
    bessm2 = BESEI0(var2+var3)
    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun

    nFLRip = normgs * gau**2. * bessm1 * bessm2 

  END FUNCTION nFLRip

  REAL(KIND=DBL) FUNCTION nFLRip1(krr)
!!! Routine for improved trapped ion response evaluation
!!! Pierre 22/2/11
!!! FBW of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    normgs = normkr * REAL(mwidth) / SQRT(pi)

    ifailloc = 0
    var1 = (normkr*krr*di(pnFLR,ion))**2.    !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoi(ion))**2.               !!2nd argument of 2nd Bessel fun
    bessm1 = BESEI1(var1)
    bessm2 = BESEI0(var2+var3)
    gau = EXP(-0.5_DBL*(krr*normkr*REAL(mwidth))**2._DBL)    !!definition of the eigenfun
    nFLRip1 = normgs * gau**2. * bessm1 *bessm2

  END FUNCTION nFLRip1

  REAL(KIND=DBL) FUNCTION nFLRecrot(krr)
!!! Routine for improved FLR, passing electrons
!!! FLR of passing particles are calculated from an integration over kr of
!!! Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var2, var3
    REAL(KIND=DBL)    :: bessm2, gau
    INTEGER :: ifailloc

    ifailloc = 0
    var2 = (normkr*krr*Rhoe(pnFLR))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoe)**2.          !!2nd argument of Bessel fun

    bessm2 = BESEI0(var2+var3)
    !gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL)    !!without rotation
    ! warning : sure that correct? do not miss 1/2 inside exp here???
    ! warning : for passing eigenfun^2 already in functionnals, shoudl not be there twice???

    !gau = EXP(-(REAL(mwidth)*krr*normkr)**2-2*krr*normkr*AIMAG(mshift)+(REAL(mshift)**2-AIMAG(mshift)**2)/REAL(mwidth))    !!definition of the eigenfun

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLRecrot = normgs * ABS(gau)**2. * bessm2 

  END FUNCTION nFLRecrot

  REAL(KIND=DBL) FUNCTION nFLReprot(krr)
!!! Routine for improved FLR, trapped electrons
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) * Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    ifailloc = 0
    var1 = (normkr*krr*de(pnFLR))**2.   !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoe(pnFLR))**2. !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoe)**2.               !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI0(var1)
    bessm2 = BESEI0(var2+var3)
    
    !gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL)    !!definition of the eigenfun
    ! warning not dimless inside exp!! pb in QLK_mom nFLRep.f90???

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLReprot = normgs * ABS(gau)**2. * bessm1 * bessm2 

  END FUNCTION nFLReprot

  REAL(KIND=DBL) FUNCTION nFLRep1rot(krr)
!!! Routine for improved  trapped electrons response evaluation
!!! Pierre 22/2/11
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    ifailloc = 0
    var1 = (normkr*krr*de(pnFLR))**2.    !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoe(pnFLR))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoe)**2.               !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI1(var1)
    bessm2 = BESEI0(var2+var3)

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLRep1rot = normgs * ABS(gau)**2. * bessm1*bessm2 

  END FUNCTION nFLRep1rot

  REAL(KIND=DBL) FUNCTION nFLRicrot(krr)
!!! Routine for improved FLR, passing ions
!!! FLR of passing particles are calculated from an integration over kr of
!!! Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var2, var3
    REAL(KIND=DBL)    :: bessm2, gau
    INTEGER :: ifailloc

    ifailloc = 0
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of Bessel fun
    var3 = (ktetaRhoi(ion))**2.               !!2nd argument of Bessel fun

    bessm2 = BESEI0(var2+var3)
 
    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLRicrot = normgs * ABS(gau)**2. * bessm2 

  END FUNCTION nFLRicrot

  REAL(KIND=DBL) FUNCTION nFLRiprot(krr)
!!! Routine for improved FLR, trapped ions
!!! FLR of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) * Bessel_mod(k_perp*rhoLr^2) * eigenfun^2
!!! using k_perp^2 = k_theta^2 + kr^2 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc

    ifailloc = 0
    var1 = (normkr*krr*di(pnFLR,ion))**2.   !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoi(ion))**2.             !!2nd argument of 2nd Bessel fun

    bessm1 = BESEI0(var1)
    bessm2 = BESEI0(var2+var3)

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLRiprot = normgs * ABS(gau)**2. * bessm1 * bessm2 

  END FUNCTION nFLRiprot

  REAL(KIND=DBL) FUNCTION nFLRip1rot(krr)
!!! Routine for improved trapped ion response evaluation
!!! Pierre 22/2/11
!!! FBW of trapped particles are calculated from an integration over kr of
!!! Bessel_mod(kr*delta^2) 
    REAL(KIND=DBL), INTENT(in) :: krr
    ! Local variables
    REAL(KIND=DBL)    :: normgs
    REAL(KIND=DBL)    :: var1, var2, var3
    REAL(KIND=DBL)    :: bessm1, bessm2, gau
    INTEGER :: ifailloc
    ifailloc  = 0
    var1 = (normkr*krr*di(pnFLR,ion))**2.    !!argument of 1st Bessel fun
    var2 = (normkr*krr*Rhoi(pnFLR,ion))**2.  !!1st argument of 2nd Bessel fun
    var3 = (ktetaRhoi(ion))**2.               !!2nd argument of 2nd Bessel fun
    bessm1 = BESEI1(var1)
    bessm2 = BESEI0(var2+var3)

    normgs = normkr * SQRT((REAL(mwidth)**2 - AIMAG(mwidth)**2)) / SQRT(pi) * EXP(-AIMAG(mshift)**2/(REAL(mwidth)**2-AIMAG(mwidth)**2))
    gau = EXP(-0.5_DBL*(krr*normkr*mwidth)**2._DBL - ci*krr*normkr*mshift)   !!definition of the eigenfun in k-space

    nFLRip1rot = normgs * ABS(gau)**2. * bessm1 *bessm2

  END FUNCTION nFLRip1rot

END MODULE FLRterms
