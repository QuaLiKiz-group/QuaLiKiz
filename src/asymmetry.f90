MODULE asymmetry
  USE kind
  USE datmat
  USE datcal
  IMPLICIT NONE
  !Module for calculating density poloidal asymmetries based on centrifugal force and temperature anisotropies (i.e. from heating)
  !Based mostly on work from Angioni POP 2012 and Bilato NF Lett. 2014

CONTAINS

  SUBROUTINE calcphi
    !!Calculates the electric potential along the field line
    !Does so by satisfying the quasineutrality constraint in the presence of rotation and temperature anisotropy
    INTEGER :: i,j,iflag
    REAL(KIND=DBL) :: guess,Btmp,Ctmp
    REAL(KIND=DBL), DIMENSION(ntheta) :: phicut

    DO irad=1,dimx !Loop over radial positions

       !Calculate normalized perpindicular temperatures along field line: Tpar(theta)/Tpar(0). theta=0 is low-field-side.
       !anise and anisi the electron and ion Tperp(0)/Tpar. These are QuaLiKiz input parameters
       !Assumed that Tpare = Tex and Tpari=Tix from QuaLiKiz inputs
       !Formula is based on Bilato 2014
       tpernorme(irad,:) = 1 / ( anise(irad)-(anise(irad)-1) * (1+epsilon(irad)*COS(th)) /(1+epsilon(irad) ) );
       DO j=1,nions
          tpernormi(irad,:,j) = 1 / ( anis(irad,j)-(anis(irad,j)-1) * (1+epsilon(irad)*COS(th)) /(1+epsilon(irad) ) );
       ENDDO

       !Solve the quasineutrality equation at each theta. DFZERO is a root solver routine in mathlib/slatecroutines.f
       !NOTE: itheta is an effective global variable defined in datmat.f90
       DO itheta=1,ntheta 
          guess = 0._DBL
          iflag=1
          !B (lower limit of solver) C (upper limit of solver), RE (relative error) ,AE (absolute error). All defined in datcal
          Btmp=B
          Ctmp=C
          IF (itheta > 1) THEN !phi(theta=0) is zero by definition
             !Solve quasineutrality! phieq is a function defined below
             CALL DFZERO(phieq,Btmp,Ctmp,guess,RE,AE,iflag)
             phi(irad,itheta)=Btmp
             IF (iflag > 2) THEN
                IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A,I2,A,I2,A,I0)") 'Abnormal termination of DFZERO at p = ',irad,'. itheta = ',itheta, '. iflag = ',iflag
             ENDIF
          ELSE
             phi(irad,itheta)=0._DBL
          ENDIF
       ENDDO
       !Calculated density profile along the field line for all species, using the just calculated phi
       phicut=phi(irad,:)
       nepol(irad,:)=densprof(phicut,0)
       DO j=1,nions        
          npol(irad,:,j)=densprof(phicut,j)
       ENDDO
    ENDDO
  END SUBROUTINE calcphi

  REAL(KIND=DBL) FUNCTION phieq(x) 
    !!phi is set by quasineutrality constraint at each poloidal position
    !This function defines the quasineutrality equation: sum(Zs*ns(phi))=0 which must be solved for phi
    !ns(phi) is calculated in the dens1 function
    !NOTE: itheta is an effective global variable defined in datmat.f90
    INTEGER :: inum
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: phitmp

    phitmp = -dens1(x,th(itheta),0)
    DO inum=1,nions
       IF (ion_type(irad,inum) == 1) THEN
          phitmp=phitmp+Zi(irad,inum)*dens1(x,th(itheta),inum)
       ENDIF
    ENDDO
    phieq = phitmp

  END FUNCTION phieq

  FUNCTION dens1(x,theta1,inum)
    !!Calculates density at theta=theta1, for phi=x
    !NOTE: itheta is an effective global variable defined in datmat.f90
    REAL(KIND=DBL), INTENT(IN) :: theta1, x
    INTEGER, INTENT(IN) :: inum !inum=0 for electrons, inum>0 for the various ion species
    REAL(KIND=DBL) :: dens1,Rlfs,Rth
    INTEGER :: ix

    Rlfs=Ro(irad)*(1+epsilon(irad)) !low field side theta=0 major radius
    Rth=Ro(irad)*(1+epsilon(irad)*COS(theta1)) !theta dependent major radius
    ix = NINT(theta1*(ntheta-1)/pi+1)  !find index in theta grid that matches theta1

    !Calculate density depending on input phi, rotation and temp anisotropy
    IF (inum == 0) THEN
       dens1=tpernorme(irad,ix)*Nex(irad)*EXP(-(-qe*x - me/2.*omegator(irad)**2*(Rth**2-Rlfs**2))/(Tex(irad)*1d3*qe));
    ELSE
       dens1=tpernormi(irad,ix,inum)*Nix(irad,inum)*&
            &EXP(-(Zi(irad,inum)*qe*x - mi(irad,inum)/2.*omegator(irad)**2*(Rth**2-Rlfs**2))/(Tix(irad,inum)*1d3*qe));
    ENDIF
  END FUNCTION dens1

  FUNCTION densprof(phi,inum)
    !!Calculates density profiles for species inum (0 electrons, >=1 for ions), and for an electron potential profile phi
    REAL(KIND=DBL), DIMENSION(ntheta), INTENT(IN) :: phi
    INTEGER, INTENT(IN) :: inum
    REAL(KIND=DBL), DIMENSION(ntheta) :: densprof, Rth
    REAL(KIND=DBL) :: Rlfs

    Rlfs=Ro(irad)*(1+epsilon(irad))
    Rth(:)=Ro(irad)*(1+epsilon(irad)*COS(th))
    IF (inum == 0) THEN
       densprof=tpernorme(irad,:)*Nex(irad)*EXP(-(-qe*phi - me/2.*omegator(irad)**2*(Rth(:)**2-Rlfs**2))/(Tex(irad)*1d3*qe));
    ELSE
       densprof=tpernormi(irad,:,inum)*Nix(irad,inum)&
            &*EXP(-(Zi(irad,inum)*qe*phi - mi(irad,inum)/2.*omegator(irad)**2*(Rth(:)**2-Rlfs**2))/(Tix(irad,inum)*1d3*qe));
    ENDIF
  END FUNCTION densprof

  SUBROUTINE calcdphi
    !Calculates dphi/dr and dphi/dtheta needed for poloidal asymmetry effect calculations.
    !Does so by solving quasineutrality of gradients: sum_j*(R/Lnj * n_j * Z_j) = 0
    !For a known phi, this solution can be done analytically
    REAL(KIND=DBL) :: Rlfs, llim,ulim,diffeps,errret
    REAL(KIND=DBL), DIMENSION(ntheta) :: Rth,nom,denom,Ee,dtpernormedr,phicut
    REAL(KIND=DBL), DIMENSION(ntheta,nions) :: dtpernormidr, Ei
    INTEGER :: j,i,iflag

    DO irad = 1,dimx
       phicut=phi(irad,:) !Assumes that phi (global variable in datmat) already calculated beforehand with calcphi routine
       Rlfs=Ro(irad)*(1+epsilon(irad)) !low field side major radius
       Rth =Ro(irad)*(1+epsilon(irad)*COS(th)) !theta dependent major radius

       Ee(:)=-qe*phi(irad,:)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)  !electron electric and rotational kinetic potential

       dtpernormedr(:) = 2*SIN(th/2.)**2*tpernorme(irad,:)**2/(Ro(irad)*(1+epsilon(irad))**2)*&
            &(1-anise(irad)-danisedr(irad)*epsilon(irad)*Ro(irad)*(1+epsilon(irad))) !Derivative of Tperpe(theta)/Tperpe(0)

       DO j=1,nions
          Ei(:,j)=Zi(irad,j)*qe*phi(irad,:)-mi(irad,j)/2.*omegator(irad)**2*(Rth**2-Rlfs**2) !ion electric and rotational kinetic potential
          dtpernormidr(:,j) = 2*SIN(th/2.)**2*tpernormi(irad,:,j)**2/(Ro(irad)*(1+epsilon(irad))**2)*&
               &(1-anis(irad,j)-danisdr(irad,j)*epsilon(irad)*Ro(irad)*(1+epsilon(irad))) !Derivative of Tperpi(theta)/Tperpi(0)
       ENDDO

       !Calculate dphi/dr from analytical solution of quasineutrality of gradients
       !nom (nominator) and denom (denominator) contain sums over species. dphi/dr = nom/denom
       nom(:) = -densprof(phicut,0)*(-Ane(irad)-Ee(:)*Ate(irad)/(Tex(irad)*qe*1d3) + & 
            & Ro(irad)*me*omegator(irad)*domegatordr(irad)/(Tex(irad)*qe*1d3)*(Rth**2-Rlfs**2) + me*Ro(irad)*omegator(irad)**2/(Tex(irad)*qe*1d3)*(Rth*COS(th)-Rlfs) &
            &  +Ro(irad)*dtpernormedr(:)/tpernorme(irad,:))

       denom(:) = densprof(phicut,0)*qe/(Tex(irad)*qe*1d3)*Ro(irad)

       DO j=1,nions
          IF (ion_type(irad,j) == 1) THEN
             nom(:) = nom(:) +  Zi(irad,j)*densprof(phicut,j)*(-Ani(irad,j)-Ei(:,j)*Ati(irad,j)/(Tix(irad,j)*qe*1d3) + & 
                  & Ro(irad)*mi(irad,j)*omegator(irad)*domegatordr(irad)/(Tix(irad,j)*qe*1d3)*(Rth**2-Rlfs**2) + mi(irad,j)*Ro(irad)*omegator(irad)**2/(Tix(irad,j)*qe*1d3)*(Rth*COS(th)-Rlfs) &
                  &  +Ro(irad)*dtpernormidr(:,j)/tpernormi(irad,:,j))

             denom(:) = denom(:) +  Zi(irad,j)**2*densprof(phicut,j)*qe/(Tix(irad,j)*qe*1d3)*Ro(irad)            
          ENDIF
       ENDDO

       dphidr(irad,:)=nom(:)/denom(:)

       !Calculate theta derivative of theta with 5 point difference formula
       !Uses phifunc function which calculates phi at arbitrary theta (needed for the numeric differentiation)
       DO j=1,ntheta
          dphidth(irad,j) = ( phifunc(th(j)-2.*dtheta)-8.*phifunc(th(j)-dtheta) + 8.*phifunc(th(j)+dtheta) - phifunc(th(j)+2.*dtheta) ) /(12.*dtheta)
       ENDDO

    ENDDO

  END SUBROUTINE calcdphi

  REAL(KIND=DBL) FUNCTION phifunc(x)
    !!Function for calculating phi at arbitrary theta for the numerical theta differentiation for dphidth
    !!quasineutrality equation at arbitrary theta defined in phieqfunc
    REAL(KIND=DBL), INTENT(IN) :: x !input arbitrary theta
    INTEGER :: i,j,iflag
    REAL(KIND=DBL) :: guess,Btmp,Ctmp,tpernormefunc
    REAL(KIND=DBL), DIMENSION(nions) :: tpernormifunc

    tpernormefunc = 1 / ( anise(irad)-(anise(irad)-1) * (1+epsilon(irad)*COS(x)) /(1+epsilon(irad) ) + epsD);
    DO j=1,nions
       tpernormifunc(j) = 1 / ( anis(irad,j)-(anis(irad,j)-1) * (1+epsilon(irad)*COS(x)) /(1+epsilon(irad) ) + epsD);
    ENDDO

    guess = 0._DBL
    iflag=1
    Btmp=B
    Ctmp=C
    thetapass= x  !thetapass is a global variable defined in datmat. Used to pass around the angle between the functions used in the calculation

    CALL DFZERO(phieqfunc,Btmp,Ctmp,guess,RE,AE,iflag)
    IF (iflag > 2) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stdout,"(A,I2,A,I2,A,I0)") 'Abnormal termination of DFZERO in phifunc at p = ',irad,'. itheta = ',itheta, '. iflag = ',iflag
    ENDIF

    phifunc = Btmp

  END FUNCTION phifunc

  REAL(KIND=DBL) FUNCTION phieqfunc(x) 
    !phi is set by quasineutrality constraint at an arbitrary theta (used for the theta differentiation)
    INTEGER :: inum
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: phitmp

    phitmp = -dens1func(x,thetapass,0)
    DO inum=1,nions
       IF (ion_type(irad,inum) == 1) THEN
          phitmp=phitmp+Zi(irad,inum)*dens1func(x,thetapass,inum)
       ENDIF
    ENDDO
    phieqfunc = phitmp

  END FUNCTION phieqfunc

  FUNCTION dens1func(x,theta1,inum)
    !Find density at an arbitrary theta (used for the theta differentiation)
    REAL(KIND=DBL), INTENT(IN) :: theta1, x
    INTEGER, INTENT(IN) :: inum
    REAL(KIND=DBL) :: dens1func,Rlfs,Rth
    INTEGER :: j

    Rlfs=Ro(irad)*(1+epsilon(irad))
    Rth=Ro(irad)*(1+epsilon(irad)*COS(theta1))

    tpernormefunc = 1 / ( anise(irad)-(anise(irad)-1) * (1+epsilon(irad)*COS(theta1)) /(1+epsilon(irad) ) );
    DO j=1,nions
       tpernormifunc(j) = 1 / ( anis(irad,j)-(anis(irad,j)-1) * (1+epsilon(irad)*COS(theta1)) /(1+epsilon(irad) ) );
    ENDDO

    IF (inum == 0) THEN
       dens1func=tpernormefunc*Nex(irad)*EXP(-(-qe*x - me/2.*omegator(irad)**2*(Rth**2-Rlfs**2))/(Tex(irad)*1d3*qe));
    ELSE
       dens1func=tpernormifunc(inum)*Nix(irad,inum)*&
            &EXP(-(Zi(irad,inum)*qe*x - mi(irad,inum)/2.*omegator(irad)**2*(Rth**2-Rlfs**2))/(Tix(irad,inum)*1d3*qe));
    ENDIF

  END FUNCTION dens1func

  SUBROUTINE calccoefs
    !calculates the e0-e6 coefficients used for transforming from low field side (LFS) to flux surface average (FSA)
    !transport coefficients and gradients. Done for each species. Also calcuates the poloidal profile of R/Ln for each species
    REAL(KIND=DBL), DIMENSION(ntheta) :: Aneprof, Rth,Ee,Ei,dtpernormedr,dtpernormidr
    REAL(KIND=DBL) :: Rlfs,theta1,relerr,thmin,thmax
    INTEGER :: ith,ifailloc
    INTEGER :: npts !output of number of integral evaluations

    !Allocate work arrays for integration routine
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    !limits for integration
    thmin=0._DBL
    thmax=pi

    DO irad=1,dimx !loop over radial locations

       DO ion=0,nions !loop over species in plasma (0=electrons, >=1 for ions)

          DO ith=1,ntheta !loop over theta grid
             theta1=th(ith)

             !Calculate theta dependent R/Ln using the non-flux-surface averaged e#1 coefficients.
             !These coefficients are defined in GKW manual appendix B, and generalized here to include
             !the temperature anisotropy
             IF (ion == 0) THEN
                Anepol(irad,ith) = Ane(irad) - ( -Ate(irad)*e11(theta1) - e21(theta1) + &
                     & Ro(irad)**2*me*omegator(irad)**2*Ate(irad)/(2*(Tex(irad)*1d3*qe))*e31(theta1) + &
                     &  Ro(irad)**3*me*omegator(irad)*domegatordr(irad)/(Tex(irad)*qe*1d3)*e31(theta1) + &
                     &  Ro(irad)**2*me*omegator(irad)**2/(2*(Tex(irad)*1d3*qe))*e41(theta1) +e61(theta1) )/e01(theta1)
             ELSE

                Anipol(irad,ith,ion) = Ani(irad,ion) - ( -Ati(irad,ion)*e11(theta1) - e21(theta1) + &
                     & Ro(irad)**2*mi(irad,ion)*omegator(irad)**2*Ati(irad,ion)/(2*(Tix(irad,ion)*1d3*qe))*e31(theta1) + &
                     &  Ro(irad)**3*mi(irad,ion)*omegator(irad)*domegatordr(irad)/(Tix(irad,ion)*qe*1d3)*e31(theta1) + &
                     &  Ro(irad)**2*mi(irad,ion)*omegator(irad)**2/(2*(Tix(irad,ion)*1d3*qe))*e41(theta1) +e61(theta1) )/e01(theta1)
             ENDIF
          ENDDO

          !DQAGSE_QLK is a 1D integration routine in mathlib/dqagse_qlk.f

          !Calculate the flux surface averaged e0 coefficient. 
!!$          CALL DQAGSE_QLK(e01,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,1),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)

          ifailloc = 1
          ecoefs(irad,ion,1) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e01,lw2,ifailloc)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 1 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF
          !Calculate the flux surface averaged e1 coefficient. 

          ifailloc = 1
          ecoefs(irad,ion,2) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e11,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e11,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,2),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 2 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF
          !Calculate the flux surface averaged e2 coefficient. 

          ifailloc = 1
          ecoefs(irad,ion,3) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e21,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e21,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,3),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 3 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF
          !Calculate the flux surface averaged e3 coefficient. PROBLEM HERE

          ifailloc = 1
          ecoefs(irad,ion,4) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e31,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e31,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,4),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 4 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF
          !Calculate the flux surface averaged e4 coefficient. 

          ifailloc = 1
          ecoefs(irad,ion,5) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e41,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e41,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,5),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 5 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF

          ecoefs(irad,ion,6)=0 !The e5 coefficient, which is 0 in Hamada coordinates

          !Calculate the flux surface averaged e6 coefficient (new coefficient not defined in GKW Manual) 
          ifailloc = 1
          ecoefs(irad,ion,7) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e61,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e61,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,7),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 7 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF

          !e0-6, then <R/Ln>, <n>, and (nmax-nmin)/<n>
          ifailloc = 1
          ecoefs(irad,ion,8) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e71,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e71,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,8),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 8 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF

          ifailloc = 1
          ecoefs(irad,ion,9) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e81,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e81,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,9),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 9 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF

          ifailloc = 1
          ecoefs(irad,ion,10) = d01ahf(thmin,thmax,epsFLR,npts,relerr,e91,lw2,ifailloc)
!!$          CALL DQAGSE_QLK(e91,thmin,thmax,0.,epsFLR,limit,ecoefs(irad,ion,10),relerr,npts,ifailloc,&
!!$               alist, blist, rlist, elist, iord, last)
          IF (ifailloc .NE. 0) THEN
             IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef failed for coef 10 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
          ENDIF
          !R/Ln (defined as R*<dn/dr>/<n>)
          IF (ion == 0) THEN
             ecoefs(irad,ion,11) = Ane(irad) - ( -Ate(irad)*ecoefs(irad,ion,2) - ecoefs(irad,ion,3) + &
                  & Ro(irad)**2*me*omegator(irad)**2*Ate(irad)/(2*(Tex(irad)*1d3*qe))*ecoefs(irad,ion,4) + &
                  &  Ro(irad)**3*me*omegator(irad)*domegatordr(irad)/(Tex(irad)*qe*1d3)*ecoefs(irad,ion,4) + &
                  &  Ro(irad)**2*me*omegator(irad)**2/(2*(Tex(irad)*1d3*qe))*ecoefs(irad,ion,5) +ecoefs(irad,ion,7) ) / ecoefs(irad,ion,1)
          ELSE

             ecoefs(irad,ion,11) = Ani(irad,ion) - ( -Ati(irad,ion)*ecoefs(irad,ion,2) - ecoefs(irad,ion,3) + &
                  & Ro(irad)**2*mi(irad,ion)*omegator(irad)**2*Ati(irad,ion)/(2*(Tix(irad,ion)*1d3*qe))*ecoefs(irad,ion,4) + &
                  &  Ro(irad)**3*mi(irad,ion)*omegator(irad)*domegatordr(irad)/(Tix(irad,ion)*qe*1d3)*ecoefs(irad,ion,4) + &
                  &  Ro(irad)**2*mi(irad,ion)*omegator(irad)**2/(2*(Tix(irad,ion)*1d3*qe))*ecoefs(irad,ion,5) +ecoefs(irad,ion,7) )/ecoefs(irad,ion,1)

          ENDIF

          !<n>
          IF (ion == 0) THEN
             ecoefs(irad,ion,12) = Nex(irad)*ecoefs(irad,ion,1)
          ELSE
             ecoefs(irad,ion,12) = Nix(irad,ion)*ecoefs(irad,ion,1)
          ENDIF

          !Asymmetry factor: (nLFS-nHFS)/<n> 
          ecoefs(irad,ion,13)=(dens1(phi(irad,1),0._DBL,ion)-dens1(phi(irad,ntheta),pi,ion))/ecoefs(irad,ion,12)

       ENDDO
    ENDDO

    !Deallocate work arrays
    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)

  END SUBROUTINE calccoefs

  REAL(KIND=DBL) FUNCTION e01(x)
    !calculate the e0 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 

    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e01 = tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)      
       e01 = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e01

  REAL(KIND=DBL) FUNCTION  e11(x)
    !calculate the e1 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix   

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)       
       e11 = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*phi(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e11 = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*phi(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e11

  REAL(KIND=DBL) FUNCTION e21(x)
    !calculate the e2 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e21 = -Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*dphidr(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e21 = Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*dphidr(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e21

  REAL(KIND=DBL) FUNCTION e31(x)
    !calculate the e3 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e31 = 1/(Ro(irad)**2)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(Rth**2-Rlfs**2)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e31 = 1/(Ro(irad)**2)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(Rth**2-Rlfs**2)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e31

  REAL(KIND=DBL) FUNCTION e41(x)
    !calculate the e4 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dRlfsdr,dRthdr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    dRlfsdr = 1
    dRthdr = COS(x)

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e41 = 2/Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(Rth*dRthdr-Rlfs*dRlfsdr)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e41 = 2/Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(Rth*dRthdr-Rlfs*dRlfsdr)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e41

  REAL(KIND=DBL) FUNCTION e61(x)
    !calculate the e6 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dtpernormidr,dtpernormedr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       dtpernormedr = -2*SIN(x/2.)**2*(epsilon(irad)*(1+epsilon(irad))*danisedr(irad)+(anise(irad)-1)/Ro(irad))/ &
            &(anise(irad)*(epsilon(irad)-epsilon(irad)*COS(x))+epsilon(irad)*COS(x)+1)**2
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e61 = Ro(irad)*dtpernormedr*EXP(-Ee/(Tex(irad)*qe*1d3))/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       dtpernormidr = -2*SIN(x/2.)**2*(epsilon(irad)*(1+epsilon(irad))*danisdr(irad,ion)+(anis(irad,ion)-1)/Ro(irad))/ &
            &(anis(irad,ion)*(epsilon(irad)-epsilon(irad)*COS(x))+epsilon(irad)*COS(x)+1)**2
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e61 = Ro(irad)*dtpernormidr*EXP(-Ei/(Tix(irad,ion)*qe*1d3))/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e61

  REAL(KIND=DBL) FUNCTION e71(x)
    !calculate the e7 (new definition= Z*qe/T * s/eps*th*dphidth)  coefficient at each poloidal position used in the effective R/Ln in the QL flux integral. 
    !Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e71 = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*smag(irad)/epsilon(irad)*x*dphidth(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e71 = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*smag(irad)/epsilon(irad)*x*dphidth(irad,ix)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e71

  REAL(KIND=DBL) FUNCTION e81(x)
    !calculate the e8 (new definition= (1+eps*cos(th))(cos(th)+s*th*sin(th))) coefficient at each poloidal position. 
    !Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e81 = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(1+epsilon(irad)*COS(x))*(COS(x)+smag(irad)*x*SIN(x))/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e81 = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(1+epsilon(irad)*COS(x))*(COS(x)+smag(irad)*x*SIN(x))/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e81

  REAL(KIND=DBL) FUNCTION e91(x)
    !calculate the e9 coefficient (like e4 but without Rth) at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dRlfsdr,dRthdr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    dRlfsdr = 1
    dRthdr = COS(x)
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e91 = 2/Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(-Rlfs*dRlfsdr)/pi!*(1+epsilon(irad)*COS(x))
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e91 = 2/Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(-Rlfs*dRlfsdr)/pi!*(1+epsilon(irad)*COS(x))
    ENDIF
  END FUNCTION e91

  REAL(KIND=DBL) FUNCTION e01d(x)
    !calculate the e0 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 

    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e01d = tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)      
       e01d = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e01d

  REAL(KIND=DBL) FUNCTION  e11d(x)
    !calculate the e1 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix   

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))

    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)       
       e11d = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*phi(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e11d = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*phi(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e11d

  REAL(KIND=DBL) FUNCTION e21d(x)
    !calculate the e2 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e21d = -Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*dphidr(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e21d = Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*dphidr(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e21d

  REAL(KIND=DBL) FUNCTION e31d(x)
    !calculate the e3 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e31d = 1/(Ro(irad)**2)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(Rth**2-Rlfs**2)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e31d = 1/(Ro(irad)**2)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(Rth**2-Rlfs**2)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e31d

  REAL(KIND=DBL) FUNCTION e41d(x)
    !calculate the e4 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dRlfsdr,dRthdr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    dRlfsdr = 1
    dRthdr = COS(x)
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e41d = 2/Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(Rth*dRthdr-Rlfs*dRlfsdr)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e41d = 2/Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(Rth*dRthdr-Rlfs*dRlfsdr)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e41d

  REAL(KIND=DBL) FUNCTION e61d(x)
    !calculate the e6 coefficient at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dtpernormidr,dtpernormedr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       dtpernormedr = -2*SIN(x/2.)**2*(epsilon(irad)*(1+epsilon(irad))*danisedr(irad)+(anise(irad)-1)/Ro(irad))/ &
            &(anise(irad)*(epsilon(irad)-epsilon(irad)*COS(x))+epsilon(irad)*COS(x)+1)**2
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e61d = Ro(irad)*dtpernormedr*EXP(-Ee/(Tex(irad)*qe*1d3))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       dtpernormidr = -2*SIN(x/2.)**2*(epsilon(irad)*(1+epsilon(irad))*danisdr(irad,ion)+(anis(irad,ion)-1)/Ro(irad))/ &
            &(anis(irad,ion)*(epsilon(irad)-epsilon(irad)*COS(x))+epsilon(irad)*COS(x)+1)**2
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e61d = Ro(irad)*dtpernormidr*EXP(-Ei/(Tix(irad,ion)*qe*1d3))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e61d

  REAL(KIND=DBL) FUNCTION e71d(x)
    !calculate the e7 (new definition= Z*qe/T * s/eps*th*dphidth)  coefficient at each poloidal position used in the effective R/Ln in the QL flux integral. 
    !Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e71d = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*qe/(Tex(irad)*qe*1d3)*smag(irad)/epsilon(irad)*x*dphidth(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e71d = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*Zi(irad,ion)*qe/(Tix(irad,ion)*qe*1d3)*smag(irad)/epsilon(irad)*x*dphidth(irad,ix)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e71d

  REAL(KIND=DBL) FUNCTION e81d(x)
    !calculate the e8 (new definition= (1+eps*cos(th))(cos(th)+s*th*sin(th))) coefficient at each poloidal position. 
    !Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e81d = -tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(1+epsilon(irad)*COS(x))*(COS(x)+smag(irad)*x*SIN(x))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e81d = tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(1+epsilon(irad)*COS(x))*(COS(x)+smag(irad)*x*SIN(x))/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e81d

  REAL(KIND=DBL) FUNCTION e91d(x)
    !calculate the e9 coefficient (like e4 but without Rth) at each poloidal position. Includes metric in integrand for flux surface average
    REAL(KIND=DBL), INTENT(IN) :: x
    REAL(KIND=DBL) :: Ee,Ei,Rlfs,Rth
    REAL(KIND=DBL) :: dRlfsdr,dRthdr
    INTEGER :: ix

    ix = NINT(x*(ntheta-1)/pi+1) 
    Rlfs = Ro(irad)*(1+epsilon(irad))
    Rth =  Ro(irad)*(1+epsilon(irad)*COS(x))
    dRlfsdr = 1
    dRthdr = COS(x)
    IF (ion == 0) THEN
       Ee = -qe*phi(irad,ix)-me/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e91d = 2/Ro(irad)*tpernorme(irad,ix)*EXP(-Ee/(Tex(irad)*qe*1d3))*(-Rlfs*dRlfsdr)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ELSE
       Ei = Zi(irad,ion)*qe*phi(irad,ix)-mi(irad,ion)/2.*omegator(irad)**2*(Rth**2-Rlfs**2)
       e91d = 2/Ro(irad)*tpernormi(irad,ix,ion)*EXP(-Ei/(Tix(irad,ion)*qe*1d3))*(-Rlfs*dRlfsdr)/pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)
    ENDIF
  END FUNCTION e91d

  REAL(KIND=DBL) FUNCTION FSAnorm(x)
    !Calculates the normalization coefficient for the wavefunction weighted FSA used in makeecoefsgau
    REAL(KIND=DBL), INTENT(IN) :: x

    FSAnorm = 1./pi*(1+epsilon(irad)*COS(x))*EXP(-1./2.*(widthhat/distan(irad,inu)*x)**2)

  END FUNCTION FSAnorm

  SUBROUTINE makeecoefsgau(p,nu)
    INTEGER, INTENT(IN) :: p, nu
    REAL(KIND=DBL) :: relerr,thmin,thmax, intnorm
    INTEGER :: ifailloc
    INTEGER :: npts !output of number of integral evaluations

    !limits for integration
    thmin=0._DBL
    thmax=pi

    irad= p 
    inu = nu

    !Calculate integration norm
    ifailloc = 1
    intnorm = d01ahf(thmin,thmax,relacc1,npts,relerr,FSAnorm,lw2,ifailloc)
    IF (ifailloc .NE. 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau normalization failed with ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
    ENDIF

    !Calculate the averages including full FSA and the Gaussian width weighting

    DO ion=0,nions !loop over species in plasma (0=electrons, >=1 for ions)

       !Calculate the flux surface averaged e0 coefficient. 
       ifailloc = 1
       ecoefsgau(irad,inu,ion,0) = d01ahf(thmin,thmax,relacc1,npts,relerr,e01d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 0 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e1 coefficient. 
       ifailloc = 1
       ecoefsgau(irad,inu,ion,1) = d01ahf(thmin,thmax,relacc1,npts,relerr,e11d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 1 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e2 coefficient. 
       ifailloc = 1
       ecoefsgau(irad,inu,ion,2) = d01ahf(thmin,thmax,relacc1,npts,relerr,e21d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 2 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e3 coefficient. 
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,3) = d01ahf(thmin,thmax,relacc1,npts,relerr,e31d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 3 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e4 coefficient. 
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,4) = d01ahf(thmin,thmax,relacc1,npts,relerr,e41d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 4 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       ecoefsgau(irad,inu,ion,5)=0 !The e5 coefficient, which is 0 in Hamada coordinates

       !Calculate the flux surface averaged e6 coefficient (new coefficient not defined in GKW Manual) 
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,6) = d01ahf(thmin,thmax,relacc1,npts,relerr,e61d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 6 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !e0-6, then <R/Ln>, <n>, and (nmax-nmin)/<n>
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,7) = d01ahf(thmin,thmax,relacc1,npts,relerr,e71d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 7 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,8) = d01ahf(thmin,thmax,relacc1,npts,relerr,e81d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 8 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       ifailloc = 1 
       ecoefsgau(irad,inu,ion,9) = d01ahf(thmin,thmax,relacc1,npts,relerr,e91d,lw2,ifailloc)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 9 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
    ENDDO

    !Normalize all calculated quantities with the integration norm
    ecoefsgau(irad,inu,:,:)=ecoefsgau(irad,inu,:,:)/intnorm


  END SUBROUTINE makeecoefsgau

  SUBROUTINE makeecoefsgaucub(p,nu)
    INTEGER, INTENT(IN) :: p, nu
    REAL(KIND=DBL) :: relerr,thmin,thmax, intnorm
    INTEGER :: ifailloc
    INTEGER :: npts !output of number of integral evaluations

    !Allocate work arrays for integration routine
    ALLOCATE(alist(limit))
    ALLOCATE(blist(limit))
    ALLOCATE(rlist(limit))
    ALLOCATE(elist(limit))
    ALLOCATE(iord(limit))

    !limits for integration
    thmin=0._DBL
    thmax=pi

    irad= p 
    inu = nu
 
   !Calculate integration norm
    CALL DQAGSE_QLK(FSAnorm,thmin,thmax,0.,epsFLR,limit,intnorm,relerr,npts,ifailloc,&
         alist, blist, rlist, elist, iord, last)
    IF (ifailloc .NE. 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau normalization failed with ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
    ENDIF

    !Calculate the averages including full FSA and the Gaussian width weighting

    DO ion=0,nions !loop over species in plasma (0=electrons, >=1 for ions)

       !Calculate the flux surface averaged e0 coefficient. 
       CALL DQAGSE_QLK(e01d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,0),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 0 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       !Calculate the flux surface averaged e1 coefficient. 
       CALL DQAGSE_QLK(e11d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,1),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 1 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       !Calculate the flux surface averaged e2 coefficient. 
       CALL DQAGSE_QLK(e21d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,2),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 2 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e3 coefficient. 
       CALL DQAGSE_QLK(e31d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,3),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 3 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF
       !Calculate the flux surface averaged e4 coefficient. 
       CALL DQAGSE_QLK(e41d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,4),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 4 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       ecoefsgau(irad,inu,ion,5)=0 !The e5 coefficient, which is 0 in Hamada coordinates

       !Calculate the flux surface averaged e6 coefficient (new coefficient not defined in GKW Manual) 
       CALL DQAGSE_QLK(e61d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,6),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 6 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       !e0-6, then <R/Ln>, <n>, and (nmax-nmin)/<n>
       CALL DQAGSE_QLK(e71d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,7),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 7 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       CALL DQAGSE_QLK(e81d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,8),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 8 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

       CALL DQAGSE_QLK(e91d,thmin,thmax,0.,epsFLR,limit,ecoefsgau(irad,inu,ion,9),relerr,npts,ifailloc,&
            alist, blist, rlist, elist, iord, last)
       IF (ifailloc .NE. 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stdout,*) 'e-coef gau failed for coef 9 in list. ifailloc=',ifailloc,'. irad=,',irad,'. ion=',ion
       ENDIF

    ENDDO

    !Normalize all calculated quantities with the integration norm
    ecoefsgau(irad,inu,:,:)=ecoefsgau(irad,inu,:,:)/intnorm

    !Deallocate work arrays
    DEALLOCATE(alist)
    DEALLOCATE(blist)
    DEALLOCATE(rlist)
    DEALLOCATE(elist)
    DEALLOCATE(iord)


  END SUBROUTINE makeecoefsgaucub


END MODULE asymmetry
