! At the moment there is a mix of various fluid solutions here as testing needs to be done

! From Pierre PhD version: cot_fluidsol (calls cot_fonctflu), cot_calcwidth, cot_calcshift 
! From old Qlk_Jon version, ana_solflu, old_modewidth
! New rigorous version: jon_fluidsol populates jon_solflu, jon_modewidth, and jon_modeshift. Testing needs to be done

MODULE mod_fluidsol
  USE kind
  USE datmat
  USE datcal
  IMPLICIT NONE

CONTAINS

  SUBROUTINE jon_fluidsol(p,nu)
    !Constructs the various coefficients used in the advanced fluid solver with rotation.
    !Carries out the direct pitch angle integrations of various frequency and geometric terms

    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL) :: dw2,sol,width,width2,width4,shift , A0,A1,A2,A3,B0,B1,U1,U2,U3,D1,D2,Cextra,Cw2,ww0,ww1,ww2,w2test1,w2test2
    COMPLEX(KIND=DBL) :: oldsol,oldwidth,oldshift,shiftfac
    COMPLEX(KIND=DBL) :: newsol,newwidth,newshift
    COMPLEX(KIND=DBL), DIMENSION(2) :: width1vec,width2vec
    REAL(KIND=DBL) :: fc2,ft2,norm,ktheta,fk,VT,WvT2,a,b,c,relerr,Wv3,Wv4,Wv5,Wv6,Wv7,Wv8,lam
    INTEGER :: i,j,npts
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly+1) :: poly
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly) :: polysol
    REAL(KIND=DBL), DIMENSION(2*ndegpoly*(ndegpoly+1)) :: warray

    COMPLEX(KIND=DBL), DIMENSION(ndegx0+1) :: polyx0
    COMPLEX(KIND=DBL), DIMENSION(ndegx0) :: polysolx0
    REAL(KIND=DBL), DIMENSION(2*ndegx0*(ndegx0+1)) :: warrayx0

    REAL(KIND=DBL), DIMENSION(nions) :: kbar,Wd1i,Wd2i,Wd3i,WvT1i,WvT2i,Wv1i,Wv2i,Wv3i,Wv4i,Wv5i,Wv6i,Wv7i,Wv8i,nwdi !for ions
    COMPLEX(KIND=DBL), DIMENSION(nions) :: rhohat,banhat,banhat2,dhat,icoef
    REAL(KIND=DBL) :: Wd1e,Wd2e,Wd3e,WvT1e,Wv1e,Wv2e,Wv3e,Wv4e,nwde !for electrons
    REAL(KIND=DBL) :: V1,V2,V3,V4, cputime1, cputime2,gamEunnorm
    REAL(KIND=DBL) :: converge = 1e-3
    INTEGER, DIMENSION(1) :: iloc  
    INTEGER :: maxiter
    LOGICAL :: x02shift,x02poly,x02width

    !!    For testing and debugging
!!$    gammaE=0.3* cthi(p,1)/cref(p)
!!$    smag = 1.0
!!$    Aui = 0.0
!!$    Aue = Aui(p,1)*cthe(p)/cthi(p,1)
!!$    Machi = 0.0

    CALL CPU_TIME(cputime1)

    !Set integration limits
    a=  0.0d0
    b = 1.0d0 !- barelyavoid
    c = 1-2.*epsilon(p)

    ktheta = ntor(p,nu)*qx(p)/(Rmin(p)*x(p))

    kbar(:) = ktheta*ABS(smag(p))/(qx(p)*Ro(p))*cthi(p,:)
    nwde = -nwg*Tex(p)
    nwdi(:) = -nwg*(-Tix(p,:)/Zi(p,:)) !opposite to the usual definition do to a reversal of the sign in the analytic derivation of the formulas

    pFFk=p !to pass rdadial coordinate into integrand functions which can only have one argument

!!$    fk = ft(p)*d01ahf(a,b,relacc1,npts,relerr,fkint,lw,ifail)
    ifail = 1
    fk = d01ahf(a,b,relacc1,npts,relerr,fkint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution fk integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    VT = d01ahf(a,b,relacc1,npts,relerr,VTint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution VT integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    WvT2 = d01ahf(a,b,relacc1,npts,relerr,WvT2int,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution WvT2 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    norm = d01ahf(a,c,relacc1,npts,relerr,normint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution norm integration at p=',p,', nu=',nu
    ENDIF

!!$    fc2 = fc(p)
!!$    ft2 = ft(p)

    fc2 = norm
    ft2 = 1-norm
    ifail = 1
    lam = d01ahf(a,c,relacc1,npts,relerr,lamint,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution lambda integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    V1 = d01ahf(a,c,relacc1,npts,relerr,V1int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V1 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    V2 = d01ahf(a,c,relacc1,npts,relerr,V2int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V2 integration at p=',p,', nu=',nu
    ENDIF

!!$    WRITE(*,*) 'p=',p,'norm,lam,V1,V2,eps=',norm,lam,V1,V2,epsilon(p)

    ifail = 1
    V3 = d01ahf(a,c,relacc1,npts,relerr,V3int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V3 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    V4 = d01ahf(a,c,relacc1,npts,relerr,V4int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V4 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv3 = d01ahf(a,c,relacc1,npts,relerr,Wv3int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv3 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv4 = d01ahf(a,c,relacc1,npts,relerr,Wv4int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv4 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv5 = d01ahf(a,c,relacc1,npts,relerr,Wv5int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv5 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv6 = d01ahf(a,c,relacc1,npts,relerr,Wv6int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv6 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv7 = d01ahf(a,c,relacc1,npts,relerr,Wv7int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv7 integration at p=',p,', nu=',nu
    ENDIF

    ifail = 1
    Wv8 = d01ahf(a,c,relacc1,npts,relerr,Wv8int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv8 integration at p=',p,', nu=',nu
    ENDIF

    Wd1e = nwde*(Ane(p) + Ate(p)*Mache(p)**2 - 2*Mache(p)*Aue(p))
    Wd2e = nwde*(Ate(p))
    Wd3e = 2*nwde*(Aue(p)-Mache(p)*Ate(p))
    Wv1e = nwde*(2.-lam/fc2)
    Wv2e = nwde*( (smag(p)-alphax(p)-0.5)*(2.-lam/fc2)-lam/(2.*fc2)*epsilon(p))
    Wv3e = nwde*Wv3
    Wv4e = nwde*Wv4
    WvT1e = nwde*fk

    Wd1i(:) = nwdi(:)*(Ani(p,:) + Ati(p,:)*Machi(p,:)**2 - 2.*Machi(p,:)*Aui(p,:))
    Wd2i(:) = nwdi(:)*(Ati(p,:))
    Wd3i(:) = 2.*nwdi(:)*(Aui(p,:)-Machi(p,:)*Ati(p,:))
    Wv1i(:) = nwdi(:)*(2.-lam/fc2)
    Wv2i(:) = nwdi(:)*( (smag(p)-alphax(p)-0.5)*(2.-lam/fc2)-lam/(2.*fc2)*epsilon(p))
    Wv3i(:) = nwdi(:)*Wv3
    Wv4i(:) = nwdi(:)*Wv4
    Wv5i(:) = nwdi(:)**2*Wv5
    Wv6i(:) = nwdi(:)**2*Wv6
    Wv7i(:) = nwdi(:)**2*Wv7
    Wv8i(:) = nwdi(:)**2*Wv8
    WvT1i(:) = nwdi(:)*fk
    WvT2i(:) = nwdi(:)*WvT2

    gamEunnorm=gammaE(p)*cref(p)/Ro(p)

    maxiter=100

!!$ x02poly = .TRUE.
!!$ x02shift = .TRUE.
!!$ x02width = .TRUE.
!!$
    x02poly = .FALSE.
    x02shift = .FALSE.
    x02width = .FALSE.

    !set initial values of shift and width
    width = distan(p,nu)
    width2 = distan(p,nu)**2
    width4 = distan(p,nu)**4
    shift = (0.,0)
    sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg


    !    DO jext=1,10 !for gentle ramping up of shift parameters to avoid oscillating solutions

    ! WRITE(*,*) width,shift,sol
    DO i=1,maxiter ! 100 is the maximum number of iterations

       oldwidth = width
       oldshift = shift
       oldsol = sol 

       !initialize certain ion coefficients
       icoef(:) = ninorm(p,:)*Zi(p,:)**2*Tex(p)/Tix(p,:)
       rhohat(:) = 1. - ((ktheta**2)*(Rhoi(p,:))**2)/4.
       banhat(:) = 1. - (di(p,:)**2)/(4.*width**2)
       dw2 = (distan(p,nu)/width)**2

       !WITH X0^2
       IF (x02poly .EQV. .TRUE.) THEN
          banhat2(:) = 1 + (di(p,:)**2)*(shift**2-width**2)/(4.*width**4) !includes shift
          dhat(:) = -distan(p,nu)**2*(shift**2-width**2)/width**4 !includes shift
       ELSE
          !WITHOUT X0^2
          banhat2(:) = 1. - (di(p,:)**2)/(4.*width**2)
          dhat(:) = distan(p,nu)**2/width**2
       ENDIF
!!$       WRITE(*,*) 'dhat',dhat


       !WIDTH TERMS
       A3 = SUM(icoef*(1.-Machi(p,:)**2)*ft2*di(p,:)**2/4.)   

       A2 = SUM(icoef*( (1.-Machi(p,:)**2)* ( (ft2*Wd1i-1.5*WvT1i)*di(p,:)**2/4. + 1.5*fc2*distan(p,nu)**2*Wv2i*rhohat) + &
            & Machi(p,:)*(Wd3i*1.5*VT*di(p,:)**2/4. -3.*Wd3i*V2*distan(p,nu)**2*rhohat ) + rhohat*Machi(p,:)**2*distan(p,nu)**2*(15./2.*Wv4i  ))) ! 2nd line subdominant

       A1 = SUM(icoef*( (1.-Machi(p,:)**2)*(1.5*(Wd1i+Wd2i)*rhohat*Wv2i*fc2*distan(p,nu)**2 - 1.5*(Wd1i+Wd2i)*WvT1i*di(p,:)**2/4.) - &
            & Machi(p,:)*2.5*Wd3i*WvT2i*di(p,:)**2/4. + &
            & 2.*Machi(p,:)*15./4.*Wd3i*Wv4i*distan(p,nu)**2*rhohat + 2*Machi(p,:)**2*15./4.*(Wd1i+2.*Wd2i)*Wv4i*rhohat*distan(p,nu)**2)) + &
            & SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*(-15./4.*Wv6i*distan(p,nu)**2) + 2.*Machi(p,:)**2*(-105./8.*Wv8i*distan(p,nu)**2)))

       A0 = SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*(-15./4.*Wv6i*distan(p,nu)**2*(Wd1i+2.*Wd2i)) - 2.*Machi(p,:)*105./8.*Wd3i*Wv8i*distan(p,nu)**2 - &
            & 2.*Machi(p,:)**2*105./8.*Wv8i*(Wd1i+3.*Wd2i)*distan(p,nu)**2))

       B1 = SUM(icoef*( (1.-Machi(p,:)**2)*1.5*V1*rhohat + Machi(p,:)**2*15./2.*V3*rhohat)*kbar**2)

       B0 = SUM(icoef*( (1.-Machi(p,:)**2)*1.5*V1*(Wd1i+Wd2i) + Machi(p,:)*15./2.*Wd3i*(V3+V4*dw2) + Machi(p,:)**2*15./2.*V3*(Wd1i+2.*Wd2i))*rhohat*kbar**2)

!!$          !! simple width
!!$          A3 = SUM(icoef*ft2*di(p,:)**2/4.)  
!!$
!!$          A2 = SUM(icoef*(  ( (ft2*Wd1i-1.5*WvT1i)*di(p,:)**2/4. + 1.5*fc2*distan(p,nu)**2*Wv2i*rhohat)))
!!$
!!$          A1 = SUM(icoef*( (1.5*(Wd1i+Wd2i)*rhohat*Wv2i*fc2*distan(p,nu)**2 - 1.5*(Wd1i+Wd2i)*WvT1i*di(p,:)**2/4.))) + &
!!$               & SUM(icoef*rhohat*( (-15./4.*Wv6i*distan(p,nu)**2)))
!!$
!!$          A0 = 0.
!!$
!!$          B1 = SUM(icoef*( 1.5*V1*rhohat)*kbar**2)
!!$
!!$          B0 = SUM(icoef*( 1.5*V1*(Wd1i+Wd2i))*rhohat*kbar**2)


       IF (x02width .EQV. .TRUE.) THEN
          Cextra = shift*ktheta*gamEunnorm/2.* &
               &  (SUM(di(p,:)**2*icoef*( (1.-Machi(p,:)**2)*(ft2*3.*sol**2+ft2*2.*sol*Wd1i-3.*sol*WvT1i-1.5*WvT1i*(Wd1i+Wd2i)))) - &
               &  SUM(icoef*Machi(p,:)*15./4.*Wd3i*WvT2i) ) + &
               &  2.*shift*ktheta*gamEunnorm*distan(p,nu)**2 * &
               & SUM(rhohat*icoef*( (1.-Machi(p,:))*(fc(p)*3.*sol*Wv2i+fc(p)*1.5*(Wd1i+Wd2i)*Wv2i-15./4.*Wv6i) + 2.*Machi(p,:)*Wd3i*(15./4.*Wv4i-3.*sol*V2 )) ) + &
               & distan(p,nu)**2*SUM(rhohat*icoef*(1.-Machi(p,:)**2)*(-shift**2*1.5*kbar**2*(sol+Wd1i+Wd2i)*V2 +3.*kbar*sol*Wd3i*V2*shift - 15.*kbar*Wd3i*Wv4i*shift)) + &
               & 2.*distan(p,nu)**2*SUM(rhohat*icoef*Machi(p,:)*(3.*kbar*V2*shift*sol*(sol+Wd1i+Wd2i) - 15./4.*kbar**2*Wd3i*V4*shift**2 - 15.*kbar*shift*Wv4i*(sol+Wd1i+Wd2i) ))
       ELSE
          Cextra =0.
       ENDIF
       Cw2 = distan(p,nu)**2*SUM(rhohat*icoef*( (1.-Machi(p,:)**2)*1.5*kbar**2*(sol+Wd1i+Wd2i)*V2 + Machi(p,:)*15./2.*kbar**2*Wd3i*V4 ))

!!$          B1 = SUM(icoef*( (1.-Machi(p,:)**2)*1.5*V1*rhohat + 2*Machi(p,:)**2*15./4.*rhohat)*kbar**2)
!!$
!!$          B0 = SUM(icoef*( (1.-Machi(p,:)**2)*1.5*V1*(Wd1i+Wd2i) + Machi(p,:)*15./2.*Wd3i*(V3+V4*dw2) + Machi(p,:)**2*15./2.*(Wd1i+2.*Wd2i))*rhohat*kbar**2)


!!$          !! shift
!!$          U1 = -ft2*3.*ktheta*sol**2 - 2.*ktheta*Wd1e*sol*ft2 + 3.*WvT1e*ktheta*sol+1.5*ktheta*WvT1e*(Wd1e+Wd2e) + &
!!$               & SUM(icoef*( (1.-Machi(p,:)**2)*( (-3.*ktheta*sol**2*ft2-2.*ktheta*sol*Wd1i*ft2+3.*ktheta*WvT1i*sol+1.5*ktheta*WvT1i*(Wd1i+Wd2i))*banhat2 + &
!!$               & rhohat*(-fc2*3.*ktheta*sol**2+3.*ktheta*sol*(Wv1i+Wv2i*dw2)*fc2-2.*fc2*Wd1i*sol*ktheta+1.5*ktheta*fc2*(Wd2i+Wd1i)*(Wv1i+Wv2i*dw2))))) + &
!!$               & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-15./4.)*(Wv5i+Wv6i*dw2)*ktheta) + 3.*ktheta*sol**2*(1.+SUM(icoef)) + &    
!!$               & SUM(2.*icoef*rhohat*Machi(p,:)*ktheta*Wd3i* ( 15./4.*(Wv3i+Wv4i*(dw2)) - 3.*V1*sol ))
!!$
!!$          U2 = SUM(icoef*3.*(kbar*sol**2*(V1+V2*dw2)+kbar*sol*(Wd1i+Wd2i)*(V1+V2*dw2))*rhohat*Machi(p,:)) + &
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)*(-15./2.*kbar*sol*(Wv3i+Wv4i*dw2) -15./2.*kbar*(Wv3i+Wv4i*dw2)*(Wd1i+2.*Wd2i) ) )
!!$
!!$          U3 = SUM(icoef*1.5*kbar*Wd3i*(V1+V2*dw2)*sol*rhohat*(1.-Machi(p,:)**2)) + &     
!!$               & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-15./2.*kbar*Wd3i*(Wv3i+Wv4i*dw2))) + &
!!$               & SUM(icoef*rhohat*Machi(p,:)**2*15./2.*kbar*Wd3i*V3*sol)
!!$
!!$          D1  = SUM(icoef*(-di(p,:)**2/width**4/2.)*(1.-Machi(p,:)**2)*(ft2*sol**3+ft2*Wd1i*sol**2 - 1.5*sol**2*WvT1i - 1.5*sol*WvT1i*(Wd1i+Wd2i))) + &
!!$               SUM(icoef*(-di(p,:)**2/width**4/2.)*1.5*Machi(p,:)*Wd3i*( VT*sol**2-2.5*WvT2i*sol )) + &
!!$               & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-3.*sol**2*distan(p,nu)**2/width**4*Wv2i*fc2-3.*sol*(Wd1i+Wd2i)*Wv2i*distan(p,nu)**2/width**4*fc2)) + &
!!$               
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)*(-15./2.*sol*Wd3i*Wv4i*distan(p,nu)**2/width**4)) + &
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)**2*(-15./2.*sol*(Wd1i+2.*Wd2i)*Wv4i*distan(p,nu)**2/width**4-15./2.*sol**2*Wv4i*distan(p,nu)**2/width**4)) + &
!!$               & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(15./2.*distan(p,nu)**2/width**4*sol*Wv6i + 15./2.*Wv6i*(Wd1i+2.*Wd2i)*distan(p,nu)**2/width**4)) + &
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)*105./4.*Wd3i*distan(p,nu)**2/width**4*Wv8i) + &
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)**2*105./4.*distan(p,nu)**2/width**4*(sol+Wd1i+3.*Wd2i)*Wv8i) + &         
!!$               & SUM(icoef*rhohat*2.*Machi(p,:)*Wd3i*distan(p,nu)**2/width**4*(V2*3.*sol**2 - 15./2.*sol*Wv4i) )
!!$
!!$          D2 = SUM(icoef*( (1.-Machi(p,:)**2)* ( gamEunnorm*( (-3.*ktheta*sol**2*ft2-2.*ktheta*sol*Wd1i*ft2+3.*ktheta*WvT1i*sol+1.5*ktheta*WvT1i*(Wd1i+Wd2i))*di(p,:)**2/(4.*width**4) + & 
!!$               & (-distan(p,nu)**2/width**4*rhohat)*( 3*ktheta*sol*Wv2i*fc2+1.5*ktheta*(Wd1i+Wd2i)*Wv2i*fc2-15./4.*Wv6i*ktheta)) + distan(p,nu)**2/width**4*rhohat*15./2.*kbar*Wv4i*Wd3i  ) + &
!!$               & (-distan(p,nu)**2/width**4*rhohat*2.*Machi(p,:))*( -15./2.*kbar*Wv4i*(Wd1i+2.*Wd2i)-15./2.*kbar*Wv4i*sol ))) + &
!!$               SUM(-icoef*rhohat*(1.-Machi(p,:)**2)*1.5*Wd3i*distan(p,nu)**2/width**4*kbar*sol*V2) + &
!!$               SUM(icoef*rhohat*2.*Machi(p,:)*gamEunnorm*ktheta*15./4.*Wd3i*Wv4i*distan(p,nu)**2/width**4)

       ww0 = sol**3*A3+sol**2*A2+sol*A1+A0 + Cextra
       ww1 = Cw2
       ww2 = B0+B1*sol

       width4 = -ww0/ww2

       width2vec(1) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*0./2.))
       width2vec(2) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*1./2.))

       IF (REAL(width2vec(1)) > 0) THEN
          width2=width2vec(1)
       ELSE
          width2=width2vec(2)
       ENDIF

       !Test if need to use w2 model
!!$          w2test1 = (-ww1 + SQRT(ww1**2-4.*ww0*ww2))/(2.*ww2)
!!$          w2test2 = (-ww1 - SQRT(ww1**2-4.*ww0*ww2))/(2.*ww2)
!!$          
!!$          IF ( ( REAL(w2test1) > 0. ) .AND. ( REAL(w2test2) < 0. )) width2=w2test1
!!$          IF ( ( REAL(w2test1) < 0. ) .AND. ( REAL(w2test2) > 0. )) width2=w2test2


       width1vec(1) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*0./2.))
       width1vec(2) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*1./2.))

       IF (REAL(width1vec(1)) > 0) THEN
          newwidth=width1vec(1)
       ELSE
          newwidth=width1vec(2)
       ENDIF

       !Alternative shift solution
       !calculate polynomial coefficients
       polyx0(1) = D2
       polyx0(2) = D1
       polyx0(3) = gamEunnorm*U1+U2+U3

!!$          IF (x02shift .EQV. .TRUE.) THEN !WITH X0^2
!!$             !**FIND ROOT OF POLYNOMIAL** 
!!$             CALL CPQR79(ndegx0,polyx0,polysolx0,ifail,warrayx0)
!!$
!!$          IF ( ABS(AIMAG(polysolx0(1))) < ABS(AIMAG(polysolx0(2)))) THEN
!!$             IF ( ABS(REAL(polysolx0(1))) < ABS(REAL(polysolx0(2)))) THEN
!!$                newshift=polysolx0(1)/1.
!!$             ELSE
!!$                newshift=polysolx0(2)/1.
!!$             ENDIF
!!$
!!$          ELSE !WITHOUT X0^2 
!!$             newshift = -(gamEunnorm*U1+U2+U3)/D1 /1.
!!$          ENDIF

       !Xavier magic: normalize the shift by the ITG dispersion relation to maintain sufficiently small widths

       !        shiftfac = 2.*nwde/(sol + nwde*Ane(p)) / (sol/(sol + SUM(ninorm(p,:)*nwdi*(Ani(p,:)+Ati(p,:)))))

!!$             shiftfac = -2.*nwde/(sol - nwde*Ane(p)) / (sol/(sol - SUM(ninorm(p,:)*nwdi*(Ani(p,:)+Ati(p,:)))))
!!$
!!$             shiftfac = 2.*nwde/(sol - nwde*Ane(p)) / (sol/(sol - SUM(ninorm(p,:)*nwdi*(Ani(p,:)+Ati(p,:)))))

       !  WRITE(*,'(4G12.4)') sol,nwde,nwdi

       !   shiftfac = 0.25

       !shiftfac = 1.


       !     IF (Machi(p,1) > 0.1) shiftfac = 1 - Ane(p)/2.


       !newshift=newshift*shiftfac;

       !calculate polynomial coefficients
       poly(1) = -1. + ft2 + SUM( icoef*(-1.+(1.-Machi(p,:)**2)*(ft2*banhat2+fc2*rhohat)+3.*Machi(p,:)**2*V1*rhohat)) !C3 in notes

       poly(2) = Wd1e*ft2-1.5*WvT1e+SUM( icoef*( (ft2*Wd1i-1.5*WvT1i)*(1.-Machi(p,:)**2)*banhat2 + & !C2 in notes
            & banhat2*1.5*Machi(p,:)*Wd3i*VT + &
            & rhohat*((Wd1i - 1.5*(Wv1i+Wv2i*dhat))*fc2*(1.-Machi(p,:)**2) + &
            & Machi(p,:)*3.*Wd3i*(V1+V2*dw2) + 2*Machi(p,:)**2*(-15./4.*(Wv3i+Wv4i*dhat)+1.5*(Wd1i+Wd2i)*V1))))

       poly(3) = -1.5*WvT1e*(Wd1e+Wd2e)+SUM(icoef*( banhat2*(1.-Machi(p,:)**2)*(-1.5*WvT1i*(Wd1i+Wd2i)) - & !C1 in notes
            & banhat2*2.5*Machi(p,:)*Wd3i*WvT2i - & 
            & rhohat*(1.5*(Wd1i+Wd2i)*(Wv1i+Wv2i*dhat)*(1.-Machi(p,:)**2)*fc2 - &  
            & Machi(p,:)*15./2.*Wd3i*(Wv3i+Wv4i*dhat)-15./2.*Machi(p,:)**2*(Wd1i+2.*Wd2i)*(Wv3i+Wv4i*dhat)))) + & 
            & SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*15./4.*(Wv5i+Wv6i*dhat) + 2*Machi(p,:)**2*105./8.*(Wv7i+dhat*Wv8i)))

       poly(4) = SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*15./4.*(Wv5i+Wv6i*dhat)*(Wd1i+2.*Wd2i) + & !C0 in notes
            & 2.*Machi(p,:)*105./8.*Wd3i*(Wv7i+Wv8i*dhat) + 2.*Machi(p,:)**2*105./8.*(Wv7i+Wv8i*dhat)*(Wd1i+3.*Wd2i)))

!!$          !SIMPLE POLY
!!$          poly(1) = -1. + ft2 + SUM( icoef*(-1.+(1.-Machi(p,:)**2)*(ft2*banhat2+fc2*rhohat))) 
!!$
!!$          poly(2) = Wd1e*ft2-1.5*WvT1e+SUM( icoef*( (ft2*Wd1i-1.5*WvT1i)*(1.-Machi(p,:)**2)*banhat2 + & 
!!$               & rhohat*((Wd1i - 1.5*(Wv1i+Wv2i*dhat))*fc2*(1.-Machi(p,:)**2))))
!!$
!!$
!!$          poly(3) = -1.5*WvT1e*(Wd1e+Wd2e)+SUM(icoef*( banhat2*(1.-Machi(p,:)**2)*(-1.5*WvT1i*(Wd1i+Wd2i)) - & !C1 in notes
!!$               & rhohat*(1.5*(Wd1i+Wd2i)*(Wv1i+Wv2i*dhat)*(1.-Machi(p,:)**2)*fc2)))
!!$          !     & SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*15./4.*(Wv5i+Wv6i*dhat)))
!!$
!!$          poly(4) = 0.!SUM(icoef*rhohat*( (1.-Machi(p,:)**2)*15./4.*(Wv5i+Wv6i*dhat)*(Wd1i+2.*Wd2i)))


       !**FIND ROOT OF POLYNOMIAL**
       ifail=1
       CALL CPQR79(ndegpoly,poly,polysol,ifail,warray)        
       IF (ifail /= 0) THEN
          IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
               &'. Abnormal termination of fluid solution CPQR root finder at p=',p,', nu=',nu
       ENDIF

       !      WRITE(*,*) 'p=',p,' . Polysol = ',polysol
       !find and set most unstable mode
       iloc = MAXLOC(AIMAG(polysol))
       newsol = polysol(iloc(1))

       width=newwidth
       !shift=newshift 
       sol=newsol

       !Update values for outside of module
       omeflu = sol
       mwidth = width
       !mshift = shift

       !DEBUG
!!$          IF ((p==13) .AND. (nu==3)) THEN
!!$             WRITE(*,*) 'width=',width
!!$             WRITE(*,*) 'shift=',shift
!!$             WRITE(*,*) 'sol=',sol
!!$          ENDIF
       !Convergence check
       IF ( ( ABS((ABS(width)-ABS(oldwidth))) / (ABS(oldwidth)+epsD) < converge) .AND. &
            ( ABS((ABS(sol)-ABS(oldsol))) / (ABS(oldsol)+epsD) < converge) ) THEN
!!$          WRITE(*,*) 'n= ',i,' number of fluid solution convergence steps for p,nu = ',p,nu
          EXIT
       ENDIF

       IF ( i == maxiter-2)  THEN
          IF (verbose .EQV. .TRUE.) THEN
             WRITE(stderr,'(A,I2,A,I2,A)') 'Warning, rot fluid solver did not converge at (p,nu)=(',p,',',nu,')'
          ENDIF
          !reset initial values of shift and width for last iteration
          width = distan(p,nu)
          width2 = distan(p,nu)**2
          width4 = distan(p,nu)**4
          sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg
       ENDIF

       !DEBUG
!!$       IF ((p==14) .AND. (nu==4)) THEN
!!$          WRITE(*,*) 'width=',width
!!$          WRITE(*,*) 'sol=',sol
!!$       ENDIF

    ENDDO

    ! set tuning factor in width
    width4 = width4 * widthtuneITG
    width2 = width2 * widthtuneITG**0.5
    width = width * widthtuneITG**0.25
    mwidth = width

    !reinitialize certain ion coefficients
    banhat(:) = 1. - (di(p,:)**2)/(4.*width**2)
    dw2 = (distan(p,nu)/width)**2

    !WITH X0^2
    IF (x02poly .EQV. .TRUE.) THEN
       banhat2(:) = 1 + (di(p,:)**2)*(shift**2-width**2)/(4.*width**4) !includes shift
       dhat(:) = -distan(p,nu)**2*(shift**2-width**2)/width**4 !includes shift
    ELSE
       !WITHOUT X0^2
       banhat2(:) = 1. - (di(p,:)**2)/(4.*width**2)
       dhat(:) = distan(p,nu)**2/width**2
    ENDIF

    !! calculate shift using solved width and solution
    U1 = -ft2*3.*ktheta*sol**2 - 2.*ktheta*Wd1e*sol*ft2 + 3.*WvT1e*ktheta*sol+1.5*ktheta*WvT1e*(Wd1e+Wd2e) + &
         & SUM(icoef*( (1.-Machi(p,:)**2)*( (-3.*ktheta*sol**2*ft2-2.*ktheta*sol*Wd1i*ft2+3.*ktheta*WvT1i*sol+1.5*ktheta*WvT1i*(Wd1i+Wd2i))*banhat2 + &
         & rhohat*(-fc2*3.*ktheta*sol**2+3.*ktheta*sol*(Wv1i+Wv2i*dw2)*fc2-2.*fc2*Wd1i*sol*ktheta+1.5*ktheta*fc2*(Wd2i+Wd1i)*(Wv1i+Wv2i*dw2))))) + &
         & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-15./4.)*(Wv5i+Wv6i*dw2)*ktheta) + 3.*ktheta*sol**2*(1.+SUM(icoef)) + &    
         & SUM(2.*icoef*rhohat*Machi(p,:)*ktheta*Wd3i* ( 15./4.*(Wv3i+Wv4i*(dw2)) - 3.*V1*sol ))

    U2 = SUM(icoef*3.*(kbar*sol**2*(V1+V2*dw2)+kbar*sol*(Wd1i+Wd2i)*(V1+V2*dw2))*rhohat*Machi(p,:)) + &
         & SUM(icoef*rhohat*2.*Machi(p,:)*(-15./2.*kbar*sol*(Wv3i+Wv4i*dw2) -15./2.*kbar*(Wv3i+Wv4i*dw2)*(Wd1i+2.*Wd2i) ) )

    U3 = SUM(icoef*1.5*kbar*Wd3i*(V1+V2*dw2)*sol*rhohat*(1.-Machi(p,:)**2)) + &     
         & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-15./2.*kbar*Wd3i*(Wv3i+Wv4i*dw2))) + &
         & SUM(icoef*rhohat*Machi(p,:)**2*15./2.*kbar*Wd3i*V3*sol)

    D1  = SUM(icoef*(-di(p,:)**2/width**4/2.)*(1.-Machi(p,:)**2)*(ft2*sol**3+ft2*Wd1i*sol**2 - 1.5*sol**2*WvT1i - 1.5*sol*WvT1i*(Wd1i+Wd2i))) + &
         SUM(icoef*(-di(p,:)**2/width**4/2.)*1.5*Machi(p,:)*Wd3i*( VT*sol**2-2.5*WvT2i*sol )) + &
         & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(-3.*sol**2*distan(p,nu)**2/width**4*Wv2i*fc2-3.*sol*(Wd1i+Wd2i)*Wv2i*distan(p,nu)**2/width**4*fc2)) + &
         
         & SUM(icoef*rhohat*2.*Machi(p,:)*(-15./2.*sol*Wd3i*Wv4i*distan(p,nu)**2/width**4)) + &
         & SUM(icoef*rhohat*2.*Machi(p,:)**2*(-15./2.*sol*(Wd1i+2.*Wd2i)*Wv4i*distan(p,nu)**2/width**4-15./2.*sol**2*Wv4i*distan(p,nu)**2/width**4)) + &
         & SUM(icoef*rhohat*(1.-Machi(p,:)**2)*(15./2.*distan(p,nu)**2/width**4*sol*Wv6i + 15./2.*Wv6i*(Wd1i+2.*Wd2i)*distan(p,nu)**2/width**4)) + &
         & SUM(icoef*rhohat*2.*Machi(p,:)*105./4.*Wd3i*distan(p,nu)**2/width**4*Wv8i) + &
         & SUM(icoef*rhohat*2.*Machi(p,:)**2*105./4.*distan(p,nu)**2/width**4*(sol+Wd1i+3.*Wd2i)*Wv8i) + &         
         & SUM(icoef*rhohat*2.*Machi(p,:)*Wd3i*distan(p,nu)**2/width**4*(V2*3.*sol**2 - 15./2.*sol*Wv4i) )

    D2 = SUM(icoef*( (1.-Machi(p,:)**2)* ( gamEunnorm*( (-3.*ktheta*sol**2*ft2-2.*ktheta*sol*Wd1i*ft2+3.*ktheta*WvT1i*sol+1.5*ktheta*WvT1i*(Wd1i+Wd2i))*di(p,:)**2/(4.*width**4) + & 
         & (-distan(p,nu)**2/width**4*rhohat)*( 3*ktheta*sol*Wv2i*fc2+1.5*ktheta*(Wd1i+Wd2i)*Wv2i*fc2-15./4.*Wv6i*ktheta)) + distan(p,nu)**2/width**4*rhohat*15./2.*kbar*Wv4i*Wd3i  ) + &
         & (-distan(p,nu)**2/width**4*rhohat*2.*Machi(p,:))*( -15./2.*kbar*Wv4i*(Wd1i+2.*Wd2i)-15./2.*kbar*Wv4i*sol ))) + &
         SUM(-icoef*rhohat*(1.-Machi(p,:)**2)*1.5*Wd3i*distan(p,nu)**2/width**4*kbar*sol*V2) + &
         SUM(icoef*rhohat*2.*Machi(p,:)*gamEunnorm*ktheta*15./4.*Wd3i*Wv4i*distan(p,nu)**2/width**4)

    shift = -(gamEunnorm*U1+U2+U3)/D1 /1.
    mshift = shift

    CALL CPU_TIME(cputime2)
!!$WRITE(*,*) i
!!$WRITE(*,'(A,I0,A,8G14.4)') 'p=',p,', gamU1/U2/U3/D1= ',gamEunnorm*U1,U2,U3,D1
!!$WRITE(*,'(A,I0,A,2G14.4)') 'p=',p,', shift= ',mshift

!!$WRITE(*,'(A,I0,A,8G14.4)') 'p=',p,', A3*sol^3,A2*sol^2,A1*sol,A0= ',A3*sol**3,A2*sol**2,A1*sol,A0
!!$WRITE(*,'(A,I0,A,4G14.4)') 'p=',p,', B1*sol,B0= ',B1*sol,B0



!!$    WRITE(*,'(A,I2,8G12.4)') 'p=',p,norm*Wv1i/nwdi,norm*Wv2i/nwdi,Wv3i/nwdi,Wv4i/nwdi,SQRT(Wv5i)/nwdi,SQRT(Wv6i)/nwdi,SQRT(Wv7i)/nwdi,SQRT(Wv8i)/nwdi
!!$    WRITE(*,'(A,I2,9G12.4)') 'p=',p,Wd1i/nwdi,Wd2i/nwdi,Wd3i/nwdi,V1,V2,V3,V4,lam,norm
!!$    WRITE(*,'(A,I2,4G12.4)') 'p=',p,width*kbar/nwdi,omeflu/nwdi


!!$WRITE(*,'(A,8G12.4)') 'A3*sol3/A2*sol2/A1*sol/A0',A3*sol**3,A2*sol**2,A1*sol,A0
!!$WRITE(*,'(A,8G12.4)') 'poly(1)*sol3/poly(2)*sol2/poly(3)*sol/poly(4)',poly(1)*sol**3,poly(2)*sol**2,poly(3)*sol,poly(4)
!!$WRITE(*,'(A,4G12.4)') 'B1*sol/B0',B1*sol,B0
!!$WRITE(*,'(A,8G12.4)')  'Wv2i/WvT1i/W*/sol',Wv2i,WvT1i,Wd1i+Wd2i,kbar*width,sol,V1
!!$WRITE(*,'(A,5G12.4)') 'di/Rhoi/ft/fc/rhohat',di(p,:),Rhoi(p,:),ft2,fc2,rhohat

!!$    WRITE(*,'(A,8G12.4)') 'poly(1)*sol3/poly(2)*sol2/poly(3)*sol/poly(4)',poly(1)*sol**3,poly(2)*sol**2,poly(3)*sol,poly(4)
!!$    WRITE(*,'(A,2G12.4)') 'Solution',poly(1)*sol**3+poly(2)*sol**2+poly(3)*sol+poly(4)
!!$    WRITE(*,*) -1.,ft2,SUM( icoef*(-1.+(1.-Machi(p,:)**2)*(ft2*banhat2+fc2*rhohat))) 
!!$    STOP


  END SUBROUTINE jon_fluidsol

  SUBROUTINE jon_fluidsol_ele(p,nu)
    !Constructs the various coefficients used in the advanced fluid solver for ky>2 electron terms.
    !Carries out the direct pitch angle integrations of various frequency and geometric terms

    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL) :: dw2,sol,width,width2,width4, A0,A1,A2
    COMPLEX(KIND=DBL) :: oldsol,oldwidth
    COMPLEX(KIND=DBL), DIMENSION(2) :: width1vec,width2vec
    REAL(KIND=DBL) :: ft2,fc2,norm,ktheta,fk,WvT2,a,b,c,relerr,Wv5,Wv6,lam,kbar
    INTEGER :: i,j,npts
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly+1) :: poly
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly) :: polysol
    REAL(KIND=DBL), DIMENSION(2*ndegpoly*(ndegpoly+1)) :: warray

    COMPLEX(KIND=DBL), DIMENSION(ndegx0+1) :: polyx0
    COMPLEX(KIND=DBL), DIMENSION(ndegx0) :: polysolx0
    REAL(KIND=DBL), DIMENSION(2*ndegx0*(ndegx0+1)) :: warrayx0

    COMPLEX(KIND=DBL) :: rhohat,banhat,dhat
    REAL(KIND=DBL) :: Wd1e,Wd2e,Wd3e,WvT1e,Wv1e,Wv2e,Wv5e,Wv6e,nwde !for electrons
    REAL(KIND=DBL) :: V1,V2,V3,V4, cputime1, cputime2,gamEunnorm
    REAL(KIND=DBL) :: converge = 1e-3
    INTEGER, DIMENSION(1) :: iloc  
    INTEGER :: maxiter

    !Set integration limits
    a=  0.0d0
    b = 1.0d0 !- barelyavoid
    c = 1-2.*epsilon(p)

    ktheta = ntor(p,nu)*qx(p)/(Rmin(p)*x(p))

    kbar = ktheta*ABS(smag(p))/(qx(p)*Ro(p))*cthe(p)
    nwde = -nwg*Tex(p)

    pFFk=p !to pass rdadial coordinate into integrand functions which can only have one argument

    fk = d01ahf(a,b,relacc1,npts,relerr,fkint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution fk integration at p=',p,', nu=',nu
    ENDIF

    norm = d01ahf(a,c,relacc1,npts,relerr,normint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution norm integration at p=',p,', nu=',nu
    ENDIF
!!$    fc2 = fc(p)
!!$    ft2 = ft(p)
    fc2 = norm
    ft2 = 1-fc2

    lam = d01ahf(a,c,relacc1,npts,relerr,lamint,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution lambda integration at p=',p,', nu=',nu
    ENDIF

    V1 = d01ahf(a,c,relacc1,npts,relerr,V1int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution V1 integration at p=',p,', nu=',nu
    ENDIF

    V2 = d01ahf(a,c,relacc1,npts,relerr,V2int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution V2 integration at p=',p,', nu=',nu
    ENDIF

    Wv5 = d01ahf(a,c,relacc1,npts,relerr,Wv5int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution Wv5 integration at p=',p,', nu=',nu
    ENDIF

    Wv6 = d01ahf(a,c,relacc1,npts,relerr,Wv6int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of ele fluid solution Wv6 integration at p=',p,', nu=',nu
    ENDIF

    Wd1e = nwde*Ane(p)
    Wd2e = nwde*Ate(p)

    Wv1e = nwde*(2.-lam/fc2)
    Wv2e = nwde*( (smag(p)-alphax(p)-0.5)*(2.-lam/fc2)-lam/(2.*fc2)*epsilon(p))
    WvT1e = nwde*fk
    Wv5e = nwde**2*Wv5
    Wv6e = nwde**2*Wv6

    maxiter=100

    !set initial values of shift and width
    width = distan(p,nu)
    width2 = distan(p,nu)**2
    width4 = distan(p,nu)**4
    sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg

    DO i=1,maxiter ! 100 is the maximum number of iterations

       oldwidth = width
       oldsol = sol 

       !initialize certain ion coefficients
       rhohat = 1. - (ktheta**2*Rhoe(p)**2)/4.
       banhat = de(p)**2/(4.*width**2)
       dw2 = (distan(p,nu)/width)**2

       !calculate polynomial coefficients
       poly(1) = -1. - Tex(p)/Tix(p,1)*Zeffx(p) + (1-banhat)*ft2 +rhohat*fc2 !C3 in notes

       poly(2) = (Wd1e*ft2-1.5*WvT1e)*(1-banhat)+rhohat*fc2*(Wd1e-1.5*(Wv1e+dw2*Wv2e)) !C2 in notes

       poly(3) = -1.5*(Wd1e+Wd2e)*WvT1e*(1-banhat)+rhohat*fc2*(15./4.*(Wv5e+Wv6e*dw2)-1.5*(Wd1e+Wd2e)*(Wv1e+Wv2e*dw2))  !C1 in notes

       poly(4) = rhohat*15./4.*fc2*(Wv5e+Wv6e*dw2)*(Wd1e+2.*Wd2e)  !C0 in notes


       !**FIND ROOT OF POLYNOMIAL**
       ifail = 1
       CALL CPQR79(ndegpoly,poly,polysol,ifail,warray)
       !      WRITE(*,*) 'p=',p,' . Polysol = ',polysol
       !find and set most unstable mode
       iloc = MAXLOC(AIMAG(polysol))
       sol = polysol(iloc(1))

       !! width

       A2 = 1.5*kbar**2*V1*(sol+Wd1e+Wd2e)

       A1 = 1.5*kbar**2*V2*distan(p,nu)**2*(sol+Wd1e+Wd2e)

       A0 = ft2*(sol**3 + sol**2*Wd1e - 1.5/ft2*sol**2*WvT1e - 1.5/ft2*sol*(Wd1e+Wd2e)*WvT1e)*de(p)**2/4. + &
            & rhohat*fc2*(1.5*sol**2*Wv2e*distan(p,nu)**2 + 1.5*sol*(Wd1e+Wd2e)*distan(p,nu)**2*Wv2e - 15./4.*Wv6e*distan(p,nu)**2*(sol+Wd1e+2.*Wd2e))

       width4 = -A0/A2

       !TESTING DEBUGGING
!!$       WRITE(*,*) width4
!!$       WRITE(*,*) (-A1 + SQRT(A1**2-4.*A0*A2))/A0
!!$       WRITE(*,*) (-A1 - SQRT(A1**2-4.*A0*A2))/A0

       width2vec(1) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*0./2.))
       width2vec(2) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*1./2.))

       IF (REAL(width2vec(1)) > 0) THEN
          width2=width2vec(1)
       ELSE
          width2=width2vec(2)
       ENDIF

       width1vec(1) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*0./2.))
       width1vec(2) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*1./2.))

       IF (REAL(width1vec(1)) > 0) THEN
          width=width1vec(1)
       ELSE
          width=width1vec(2)
       ENDIF

       !Update values for outside of module
       omeflu = sol
       mwidth = width
       mshift = 0.

       !Convergence check
       IF ( ( ABS((ABS(width)-ABS(oldwidth))) / (ABS(oldwidth)+epsD) < converge) .AND. &
            ( ABS((ABS(sol)-ABS(oldsol))) / (ABS(oldsol)+epsD) < converge) ) THEN
          EXIT
       ENDIF

       IF (i==maxiter-2) THEN 
          IF (verbose .EQV. .TRUE.) THEN
             WRITE(stderr,'(A,I2,A,I2,A)') 'Warning, electron fluid solution did not converge at (p,nu)=(',p,',',nu,'). Reverting to non self-consistent solution'
          ENDIF
          !reset initial values of shift and width for last iteration
          width = distan(p,nu)
          width2 = distan(p,nu)**2
          width4 = distan(p,nu)**4
          sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg
       ENDIF

       !DEBUG
!!$       IF ((p==20) .AND. (nu==16)) THEN
!!$          WRITE(*,*) 'width=',width
!!$          WRITE(*,*) 'sol=',sol
!!$       ENDIF

    ENDDO

    ! tuning factor for width
    width = width*widthtuneETG**0.25
    mwidth = mwidth*widthtuneETG**0.25

  END SUBROUTINE jon_fluidsol_ele


  SUBROUTINE jon_fluidsol_norot(p,nu)
    !Constructs the various coefficients used in the advanced fluid solver with rotation.
    !Sets all rotation related terms to zero
    !Carries out the direct pitch angle integrations of various frequency and geometric terms

    INTEGER, INTENT(IN) :: p, nu
    COMPLEX(KIND=DBL) :: sol,width,width2,width4 , A0,A1,A2,A3,B0,B1,U1,U2,U3,D1,D2
    COMPLEX(KIND=DBL) :: oldsol,oldwidth
    COMPLEX(KIND=DBL), DIMENSION(2) :: width1vec,width2vec
    REAL(KIND=DBL) :: ft2,fc2,norm,ktheta,fk,VT,WvT2,a,b,c,relerr,Wv3,Wv4,Wv5,Wv6,Wv7,Wv8,lam
    INTEGER :: i,npts,maxiter
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly+1) :: poly
    COMPLEX(KIND=DBL), DIMENSION(ndegpoly) :: polysol
    REAL(KIND=DBL), DIMENSION(2*ndegpoly*(ndegpoly+1)) :: warray

    COMPLEX(KIND=DBL), DIMENSION(ndegx0+1) :: polyx0
    COMPLEX(KIND=DBL), DIMENSION(ndegx0) :: polysolx0
    REAL(KIND=DBL), DIMENSION(2*ndegx0*(ndegx0+1)) :: warrayx0

    REAL(KIND=DBL), DIMENSION(nions) :: kbar,Wd1i,Wd2i,Wd3i,WvT1i,WvT2i,Wv1i,Wv2i,Wv3i,Wv4i,Wv5i,Wv6i,Wv7i,Wv8i,nwdi !for ions
    COMPLEX(KIND=DBL), DIMENSION(nions) :: rhohat,banhat,banhat2,dhat,icoef
    REAL(KIND=DBL) :: Wd1e,Wd2e,Wd3e,WvT1e,Wv1e,Wv2e,Wv3e,Wv4e,nwde !for electrons
    REAL(KIND=DBL) :: V1,V2,V3,V4, cputime1, cputime2
    REAL(KIND=DBL) :: converge = 1e-3
    INTEGER, DIMENSION(1) :: iloc    
    LOGICAL :: x02shift = .FALSE.
    LOGICAL :: x02poly = .FALSE.

    !!    For testing and debugging
!!$    gammaE=0.3* cthi(p,1)/cref(p)
!!$    smag = 1.0
!!$    Aui = 0.0
!!$    Aue = Aui(p,1)*cthe(p)/cthi(p,1)
!!$    Machi = 0.0

    !IF ( SUM(ninorm(p,:)*(Machi(p,:)+Aui(p,:)

    maxiter=100
    CALL CPU_TIME(cputime1)

    !Set integration limits
    a=  0.0d0
    b = 1.0d0 !- barelyavoid
    c = 1-2.*epsilon(p)

    ktheta = ntor(p,nu)*qx(p)/(Rmin(p)*x(p))

    kbar(:) = ktheta*ABS(smag(p))/(qx(p)*Ro(p))*cthi(p,:)
    nwde = -nwg*Tex(p)
    nwdi(:) = -nwg*(-Tix(p,:)/Zi(p,:)) !opposite to the usual definition do to a reversal of the sign in the analytic derivation of the formulas

    pFFk=p !to pass radial coordinate into integrand functions which can only have one argument

!!$    fk = ft2*d01ahf(a,b,relacc1,npts,relerr,fkint,lw,ifail)
    fk = d01ahf(a,b,relacc1,npts,relerr,fkint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution fk integration at p=',p,', nu=',nu
    ENDIF

    VT = d01ahf(a,b,relacc1,npts,relerr,VTint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution VT integration at p=',p,', nu=',nu
    ENDIF

    WvT2 = d01ahf(a,b,relacc1,npts,relerr,WvT2int,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution WvT2 integration at p=',p,', nu=',nu
    ENDIF

    norm = d01ahf(a,c,relacc1,npts,relerr,normint,lw,ifail)
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution norm integration at p=',p,', nu=',nu
    ENDIF
!!$    fc2 = fc(p)
!!$    ft2 = ft(p)
    fc2 = norm
    ft2 = 1-fc2

    lam = d01ahf(a,c,relacc1,npts,relerr,lamint,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution lambda integration at p=',p,', nu=',nu
    ENDIF

    V1 = d01ahf(a,c,relacc1,npts,relerr,V1int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V1 integration at p=',p,', nu=',nu
    ENDIF

    V2 = d01ahf(a,c,relacc1,npts,relerr,V2int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V2 integration at p=',p,', nu=',nu
    ENDIF

    V3 = d01ahf(a,c,relacc1,npts,relerr,V3int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V3 integration at p=',p,', nu=',nu
    ENDIF

    V4 = d01ahf(a,c,relacc1,npts,relerr,V4int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution V4 integration at p=',p,', nu=',nu
    ENDIF

    Wv3 = d01ahf(a,c,relacc1,npts,relerr,Wv3int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv3 integration at p=',p,', nu=',nu
    ENDIF

    Wv4 = d01ahf(a,c,relacc1,npts,relerr,Wv4int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv4 integration at p=',p,', nu=',nu
    ENDIF

    Wv5 = d01ahf(a,c,relacc1,npts,relerr,Wv5int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv5 integration at p=',p,', nu=',nu
    ENDIF

    Wv6 = d01ahf(a,c,relacc1,npts,relerr,Wv6int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv6 integration at p=',p,', nu=',nu
    ENDIF

    Wv7 = d01ahf(a,c,relacc1,npts,relerr,Wv7int,lw,ifail)!*fc(p)/norm
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv7 integration at p=',p,', nu=',nu
    ENDIF

    Wv8 = d01ahf(a,c,relacc1,npts,relerr,Wv8int,lw,ifail)!*fc(p)/norm                 
    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of fluid solution Wv8 integration at p=',p,', nu=',nu
    ENDIF

!!$    !DEBUGGING
!!$    WRITE(*,*) fk,lam,V1,V3,Wv3,Wv4
!!$    WRITE(*,*) Wv5,Wv6,Wv7,Wv8
!!$    STOP

    Wd1e = nwde*Ane(p)
    Wd2e = nwde*Ate(p)
    Wd3e = 0.
    Wv1e = nwde*(2.-lam/fc2)
    Wv2e = nwde*( (smag(p)-alphax(p)-0.5)*(2.-lam/fc2)-lam/(2.*fc2)*epsilon(p))
    Wv3e = nwde*Wv3
    Wv4e = nwde*Wv4
    WvT1e = nwde*fk

    Wd1i(:) = nwdi(:)*Ani(p,:)
    Wd2i(:) = nwdi(:)*(Ati(p,:))
    Wd3i(:) = 0.
    Wv1i(:) = nwdi(:)*(2.-lam/fc2)
    Wv2i(:) = nwdi(:)*( (smag(p)-alphax(p)-0.5)*(2.-lam/fc2)-lam/(2.*fc2)*epsilon(p))
    Wv3i(:) = nwdi(:)*Wv3
    Wv4i(:) = nwdi(:)*Wv4
    Wv5i(:) = nwdi(:)**2*Wv5
    Wv6i(:) = nwdi(:)**2*Wv6
    Wv7i(:) = nwdi(:)**2*Wv7
    Wv8i(:) = nwdi(:)**2*Wv8
    WvT1i(:) = nwdi(:)*fk
    WvT2i(:) = nwdi(:)*WvT2

    !set initial values of width and sol
    width = distan(p,nu)
    width2 = distan(p,nu)**2
    width4 = distan(p,nu)**4
    sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg

    DO i=1,maxiter ! 100 is the maximum number of iterations

       oldwidth = width
       oldsol = sol

       !initialize certain ion coefficients
       icoef(:) = ninorm(p,:)*Zi(p,:)**2*Tex(p)/Tix(p,:)
       rhohat(:) = 1. - ((ktheta**2)*(Rhoi(p,:))**2)/4.
       banhat(:) = 1. - (di(p,:)**2)/(4.*width2)
       dhat(:) = distan(p,nu)**2/width2

       !calculate polynomial coefficients
       poly(1) = -1. + ft2 + SUM( icoef*(-1.+ft2*banhat+fc2*rhohat)) !C3 in notes

       poly(2) = Wd1e*ft2-1.5*WvT1e+SUM( icoef*( (ft2*Wd1i-1.5*WvT1i)*banhat + & !C2 in notes
            & rhohat*((Wd1i - 1.5*(Wv1i+Wv2i*dhat))*fc2)))

       poly(3) = -1.5*WvT1e*(Wd1e+Wd2e)+SUM(icoef*( banhat*(-1.5*WvT1i*(Wd1i+Wd2i)) - & !C1 in notes
            & rhohat*(1.5*(Wd1i+Wd2i)*(Wv1i+Wv2i*dhat)*fc2))) + &
            & SUM(icoef*rhohat*( 15./4.*(Wv5i+Wv6i*dhat)))

       poly(4) = SUM(icoef*rhohat*(15./4.*(Wv5i+Wv6i*dhat)*(Wd1i+2.*Wd2i))) !C0 in notes

       !      WRITE(*,*) 'p=',p,' . Poly = ',poly

       !**FIND ROOT OF POLYNOMIAL**
       CALL CPQR79(ndegpoly,poly,polysol,ifail,warray)
       !      WRITE(*,*) 'p=',p,' . Polysol = ',polysol
       !find and set most unstable mode
       iloc = MAXLOC(AIMAG(polysol))
       sol = polysol(iloc(1))

       !! width
       A3 = SUM(icoef*ft2*di(p,:)**2/4.)

       A2 = SUM(icoef*(( (ft2*Wd1i-1.5*WvT1i)*di(p,:)**2/4. + 1.5*fc2*distan(p,nu)**2*Wv2i*rhohat)))

       A1 = SUM(icoef*( (1.5*(Wd1i+Wd2i)*rhohat*Wv2i*fc2*distan(p,nu)**2 - 1.5*(Wd1i+Wd2i)*WvT1i*di(p,:)**2/4.))) + &
            & SUM(icoef*rhohat*( (-15./4.*Wv6i*distan(p,nu)**2)))

       A0 = SUM(icoef*rhohat*( (-15./4.*Wv6i*distan(p,nu)**2*(Wd1i+2.*Wd2i))))

       B1 = SUM(icoef*1.5*V1*rhohat*kbar**2)

       B0 = SUM(icoef*( 1.5*V1*(Wd1i+Wd2i))*rhohat*kbar**2)

       width4 = -((sol**3*A3+sol**2*A2+sol*A1+A0)/(B0+B1*sol))

       !WRITE(*,*) 'p=',p,'B0/B1/A0/A1/A2/A3=',B0,B1,A0,A1,A2,A3     

       width2vec(1) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*0./2.))
       width2vec(2) = ABS(width4)**0.5*EXP(ci*(ATAN(AIMAG(width4)/REAL(width4))/2. + 2.*pi*1./2.))

       IF (REAL(width2vec(1)) > 0) THEN
          width2=width2vec(1)
       ELSE
          width2=width2vec(2)
       ENDIF

       width1vec(1) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*0./2.))
       width1vec(2) = ABS(width2)**0.5*EXP(ci*(ATAN(AIMAG(width2)/REAL(width2))/2. + 2.*pi*1./2.))

       IF (REAL(width1vec(1)) > 0) THEN
          width=width1vec(1)
       ELSE
          width=width1vec(2)
       ENDIF

       omeflu = sol
       mwidth = width
       mshift = 0. ! By definition

       IF ( ( ABS((ABS(width)-ABS(oldwidth))) / (ABS(oldwidth)+epsD) < converge) .AND. &
            ( ABS((ABS(sol)-ABS(oldsol))) / (ABS(oldsol)+epsD) < converge) ) EXIT

       IF (i==maxiter-2) THEN
          IF (verbose .EQV. .TRUE.) THEN 
             WRITE(stderr,'(A,I2,A,I2,A)') 'Warning: norot fluid solver did not converge at (p,nu)=(',p,',',nu,'). Reverting to non-self-consistent solution'
             !       STOP 
          ENDIF
          !reset initial values of shift and width for last iteration
          width = distan(p,nu)
          width2 = distan(p,nu)**2
          width4 = distan(p,nu)**4
          sol =CMPLX(REAL(ana_solflu(p,nu)),AIMAG(ana_solflu(p,nu)))*nwg
       ENDIF

       !DEBUG
!!$       IF ((p==11) .AND. (nu==7)) THEN
!!$          WRITE(*,*) 'width=',width
!!$          WRITE(*,*) 'sol=',sol
!!$       ENDIF

    ENDDO

    ! tuning factor for width
    width = width*widthtuneITG**0.25
    mwidth = mwidth*widthtuneITG**0.25

    CALL CPU_TIME(cputime2)

  END SUBROUTINE jon_fluidsol_norot


  REAL(KIND=DBL) FUNCTION fkint(kk)
    !integrand for <f(k)>, the pinch angle integration of the vertical drift frequency for trapped particles
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL)    :: xa, ya, za, var, k2
    REAL(KIND=DBL)    :: fki, Eg, Kg
    INTEGER :: ifail
    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    xa = 0.
    ya = 1.-k2
    za = 1.
    ifail = 0
    Kg = rf(xa, ya, za, ifail)
    ifail = 0
    var   = rd(xa, ya, za, ifail)
    Eg = Kg - (k2)/3. * var

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    ! fki = 2.*Eg/Kg - 1. + 4.*smag(pFFk)*(k2-1.+Eg/Kg)-4.*alphax(pFFk)/3.* &
    !     & (1.-k2*(2.*k2-1)+Eg/Kg)

    fkint = fki*Kg*kk

!!$    fkint = Kg*kk  ! to test the normalization

  END FUNCTION fkint

  REAL(KIND=DBL) FUNCTION VTint(kk)
    !integrand for <vpar^2>, the pinch angle integration of the bounce averaged vpar^2 for trapped particles
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL)    :: xa, ya, za, var, k2
    REAL(KIND=DBL)    :: fki, E2g, Eg, Kg
    INTEGER :: ifail
    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    xa = 0.
    ya = 1.-k2
    za = 1.
    ifail = 0
    Kg = rf(xa, ya, za, ifail)
    ifail = 0
    var   = rd(xa, ya, za, ifail)
    Eg = Kg - (k2)/3. * var
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    ! fki = 2.*Eg/Kg - 1. + 4.*smag(pFFk)*(k2-1.+Eg/Kg)-4.*alphax(pFFk)/3.* &
    !     & (1.-k2*(2.*k2-1)+Eg/Kg)

    VTint = 2.*epsilon(pFFk)*E2g*k2

  END FUNCTION VTint

  REAL(KIND=DBL) FUNCTION WvT2int(kk)
    !integrand for <fk*vpar^2>, the pinch angle integration of the bounce averaged fk*vpar^2 for trapped particles
    REAL(KIND=DBL), INTENT(IN) :: kk
    REAL(KIND=DBL)    :: xa, ya, za, var, k2
    REAL(KIND=DBL)    :: fki, E2g, Eg, Kg
    INTEGER :: ifail
    k2 = kk*kk
    ! The term weighting the vertical drift of the trapped (fk) is calculated 
    ! The formulation with elliptic integrals is used
    xa = 0.
    ya = 1.-k2
    za = 1.
    ifail = 0
    Kg = rf(xa, ya, za, ifail)
    ifail = 0
    var   = rd(xa, ya, za, ifail)
    Eg = Kg - (k2)/3. * var
    E2g = 1./kk * (Eg - Kg*(1.-k2)) !Specialized form of incomplete 2nd elliptic integral. Used for bounce average of Vpar^2

    fki = -1. + (smag(pFFk)*4. + 4./3. * alphax(pFFk)) * &
         &     (k2-1.+Eg/Kg) + 2.*Eg/Kg *(1-4./3. * k2 * alphax(pFFk))

    ! fki = 2.*Eg/Kg - 1. + 4.*smag(pFFk)*(k2-1.+Eg/Kg)-4.*alphax(pFFk)/3.* &
    !     & (1.-k2*(2.*k2-1)+Eg/Kg)

    WvT2int = 2.*epsilon(pFFk)*E2g*k2*fki

  END FUNCTION WvT2int

  REAL(KIND=DBL) FUNCTION normint(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution norm integration at p=',pFFk
    ENDIF

    normint = 1./(4.*pi)*Tlam
  END FUNCTION normint

  REAL(KIND=DBL) FUNCTION lamint(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    lamint = 1./(4.*pi)*Tlam*lamin

  END FUNCTION lamint

  REAL(KIND=DBL) FUNCTION V1int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    V1int = 1./(4.*pi)*Tlam*(1.-lamin)

  END FUNCTION V1int

  REAL(KIND=DBL) FUNCTION V2int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    V2int = -1./(4.*pi)*Tlam*lamin*epsilon(pFFk)/2.

  END FUNCTION V2int


  REAL(KIND=DBL) FUNCTION V3int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    V3int = 1./(4.*pi)*Tlam*(1.-2.*lamin+lamin**2)

  END FUNCTION V3int

  REAL(KIND=DBL) FUNCTION V4int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    V4int = -1./(4.*pi)*Tlam*(lamin-lamin**2)*epsilon(pFFk)

  END FUNCTION V4int


  REAL(KIND=DBL) FUNCTION Wv3int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv3int = 1./(4.*pi)*Tlam*(1.-lamin)*(2.-lamin)

  END FUNCTION Wv3int

  REAL(KIND=DBL) FUNCTION Wv4int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv4int = 1./(4.*pi)*Tlam*( ( (smag(pFFk)-alphax(pFFk)-0.5)*(2.-lamin)-lamin*epsilon(pFFk)/2.)*(1.-lamin) & 
         & - (2.-lamin)*lamin*epsilon(pFFk)/2.)

  END FUNCTION Wv4int

  REAL(KIND=DBL) FUNCTION Wv5int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv5int = 1./(4.*pi)*Tlam*(2.-lamin)**2

  END FUNCTION Wv5int

  REAL(KIND=DBL) FUNCTION Wv6int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv6int = 1./(4.*pi)*Tlam*( ( (smag(pFFk)-alphax(pFFk)-0.5)*(2.-lamin)-lamin*epsilon(pFFk)/2.))*2.*(2.-lamin)

  END FUNCTION Wv6int

  REAL(KIND=DBL) FUNCTION Wv7int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv7int = 1./(4.*pi)*Tlam*(2.-lamin)**2*(1.-lamin)

  END FUNCTION Wv7int

  REAL(KIND=DBL) FUNCTION Wv8int(lamin)
    !integrand for <lambda> in the pinch angle integration for passing particles
    REAL(KIND=DBL), INTENT(IN) :: lamin
    REAL(KIND=DBL) :: Tlam,relerr
    INTEGER :: ifail2,i
    REAL(KIND=DBL), DIMENSION(ntheta) :: theta,Tint

    theta=(/((i * pi/(ntheta-1)),i=0,ntheta-1)/) !set theta = [0,pi]
    Tint = 2./SQRT(1-lamin*(1+2.*epsilon(pFFk)*SIN(theta/2.)**2))

    ifail2=0
    CALL davint(theta,Tint,ntheta,0.,pi,Tlam,ifail2) !calculate transit time

    IF (ifail2 /= 1) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I3,A,I7)") 'ifail2 = ',ifail2,&
            &'. Abnormal termination of fluid solution T integration at p=',pFFk
    ENDIF

    Wv8int = 1./(4.*pi)*Tlam*( ( (smag(pFFk)-alphax(pFFk)-0.5)*(2.-lamin)-lamin*epsilon(pFFk)/2.)*(1.-lamin)*2.*(2.-lamin) & 
         & - (2.-lamin)**2*(lamin*epsilon(pFFk)/2.))

  END FUNCTION Wv8int

  SUBROUTINE ana_fluidsol( p, nu, gammasol )  
    ! Simple solution for fluid growth rate, based on maximum growth rate from all possible 
    ! ionic or electron slab or interchange analytial growth rates
    ! NOTE: only provides imaginary part
    INTEGER, INTENT(IN)  :: p, nu
    REAL(KIND=DBL), INTENT(OUT) :: gammasol

    ! Variables locales
    REAL(KIND=DBL) :: Lpe,wge,wpe,wne, Vte, de
    REAL(KIND=DBL) :: Ls,kteta, rhoeffsum, deff, ceffsum, ceff, tauoverZbar, gameffsum,gameff
    REAL(KIND=DBL) :: interion,interel,slabion,slabel,inter,slab
    REAL(KIND=DBL), DIMENSION(nions) :: Lpi,wgi,wpi,wni

    INTEGER :: i

    ! Define length and time scales of problem
    kteta = ntor(p,nu) * ktetasn(p)
    Ls = Ro(p)*qx(p)/smag(p)

    wge = wg(p)*Tex(p) 
    wpe = wge * (Ate(p) + Ane(p))
    wne = wge * Ane(p)
    Lpe = Ro(p) / (Ate(p)+Ane(p))
    DO i = 1,nions
       Lpi(i) = Ro(p) / (Ati(p,i)+Ani(p,i))
       wgi(i) = -wg(p)*Tix(p,i)/Zi(p,i)
       wpi(i) = wgi(i) * (Ati(p,i) + Ani(p,i))
       wni(i) = wgi(i) * Ani(p,i)
    ENDDO

    rhoeffsum = 0
    ceffsum = 0
    tauoverZbar = 0
    gameffsum = 0
    DO i = 1,nions
       IF (ion_type(p,i) .NE. 3) THEN
          rhoeffsum  = rhoeffsum + tau(p,i)*ninorm(p,i)*Zi(p,i)**2. * Rhoi(p,i)**2. !This is rhoeff**2
          ceffsum  = ceffsum + tau(p,i)*ninorm(p,i)*Zi(p,i)**2. * 2*Tex(p)*1d3*qe/mi(p,i)
          tauoverZbar  = tauoverZbar + 1/tau(p,i)*ninorm(p,i)*Lpe/Lpi(i)
          gameffsum = gameffsum + tau(p,i)*ninorm(p,i)*Zi(p,i)**2 * wpi(i)*wgi(i)
       ENDIF
    ENDDO

    Vte = SQRT(2*Tex(p)*1d3*qe/me)
    de = qx(p)/SQRT(( 2*epsilon(p) )) * Rhoe(p)
    deff = SQRT(rhoeffsum*(1+ft(p)/fc(p)*(qx(p)**2)/(2*epsilon(p))))
    !WRITE(*,*) 'deff1 = ',deff
    !deff = SQRT(deffsum*(1+ft(p)/fc(p)/(qx(p)**2)*(2*epsilon(p))))
    !WRITE(*,*) 'deff2 = ',deff,qx(p)**2/2/epsilon(p)
    !WRITE(*,*) Ro(p),Ls/qx(p)*(1/tauoverzbar*ft(p)+1)*SQRT(2*epsilon(p)/(ft(p)*fc(p)))
    ceff = SQRT(ceffsum)

    gameff = SQRT(gameffsum / tauoverZbar)

    IF ( ETG_flag(nu) .EQV. .FALSE. ) THEN !ITG/TEM case
       interel = SQRT(ABS((ft(p)+tauoverzbar)*wge*wpe/fc(p)))/wg(p)
       interion = SQRT(ABS((ft(p)+tauoverzbar)*wgi(1)*wpi(1)/fc(p)) )/wg(p)
       slabel = SQRT( ABS(tauoverzbar*ntor(p,nu)*wpe*kteta*deff*ceff/(2*Ls))  ) / (ntor(p,nu) * wg(p))
       slabion = SQRT( ABS(tauoverzbar*ntor(p,nu)*wpi(1)*kteta*deff*ceff/(2*Ls))  ) / (ntor(p,nu) * wg(p))

       !gammasol = MAX( MAX(interel,interion),MAX(slabel,slabion) )
        gammasol = MAX(interel,interion )
      !DEBUGGING
       !WRITE(*,*) 'ITG/TEM ',interel, interion, slabel, slabion, gammasol
    ELSE !ETG case
       inter = SQRT(ABS(wpe*wge/wg(p)**2 / (tau(p,1)*Zeffx(p))))
       slab = SQRT(ABS(ntor(p,nu)*wpe*qx(p)*kteta*Vte*Rhoe(p)/(2*Ls) / (tau(p,1)*Zeffx(p))))/(ntor(p,nu)*wg(p))
       gammasol = MAX(inter,slab)
       !DEBUGGING
       !WRITE(*,*) 'ETG ',inter,slab,gammasol
    ENDIF

  END SUBROUTINE ana_fluidsol

  SUBROUTINE old_calcwidth(ETGflag,p)
    !************************************************
    ! Mode width from appendix Citrin-Cottier PoP 2012
    !************************************************
    !The mode width depends on the regime
    !For an ETG dominated regime, electron mode widths are calculated
    !Otherwise ion scale mode widths are calculated
    !A flag is in place to switch between the regimes, according to the input parameters
    REAL(kind=DBL) :: mwidth4, mwidth2, gamsn2, sumd, sumV
    LOGICAL, INTENT(IN) :: ETGflag
    INTEGER, INTENT(IN) :: p
    REAL(kind=DBL), DIMENSION(dimx,nions) :: Nixd !needed to zero out trace ions
    INTEGER :: i

    DO i = 1,nions
       IF (ion_type(p,i) .NE. 3) THEN
          Nixd(:,i) = Nix(:,i)
       ELSE
          Nixd(:,i) = 0.
       ENDIF
    ENDDO

    IF (ETGflag .EQV. .TRUE.) THEN

       gamsn2 = wg(p)*wg(p)*Joe2*(Ate(p)+Ane(p))*Nex(p)/    &
            (Nex(p)+SUM(Zi(p,:)**2*Nixd(p,:)*tau(p,:)))

       sumd =(de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5))*Joe2*Jobane2

       sumV = Joe2 / (cthe(p)*cthe(p))*gamsn2/(ktetasn(p)*smag(p)/Ro(p))**2.

       mwidth4  = sumd*sumV*(qx(p))**2. *(Tex(p))**2.
       mwidth2 = SQRT(ABS(mwidth4))

       ! 4th order iterative expansion in theta to improve the s~0.5 behaviour
       sumd = (de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5)+ &
            2.*d**4/(mwidth2)*(smag(p)-2.*alphax(p)-0.25))*Joe2*Jobane2

       mwidth4  = sumd*sumV*(qx(p))**2. *(Tex(p))**2.
       mwidth2 = SQRT(ABS(mwidth4))

       sumd = (de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5)+ &
            2.*d**4./(mwidth2)*(smag(p)-2.*alphax(p)-0.25) + &
            5./40.*d**6./ABS(mwidth4)*(6.*smag(p)-32.*alphax(p)-1.))*Joe2*Jobane2
       mwidth4  = sumd*sumV*(qx(p))**2. *(Tex(p))**2.

       mwidth=(ABS(mwidth4))**0.25

    ELSE !IF NOT ETG THEN

       gamsn2 = wg(p)*wg(p) *(Nex(p)*Joe2*(Ate(p)+Ane(p)) + &		
            SUM(Nixd(p,:)*Joi2(:)/tau(p,:)*(Ati(p,:)+Ani(p,:))))/ &
            (Nex(p)+SUM(Zi(p,:)**2.*Nixd(p,:)*tau(p,:)))

       sumd =  ( SUM(Nixd(p,:)*((di(p,:)*di(p,:)+2.*d**2. *(smag(p)-alphax(p)-0.5))*Joi2(:)*Jobani2(:))) + &
            Nex(p)*((de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5))*Joe2*Jobane2))/ &
            (Nex(p)+SUM(Nixd(p,:)))

       sumV = ( SUM(Nixd(p,:) * Joi2(:) / (cthi(p,:)*cthi(p,:)/2.)) + &
            Nex(p) * Joe2 / (cthe(p)*cthe(p)/2.) ) * &
            gamsn2/(ktetasn(p)*smag(p)/Ro(p))**2/ &
            (Nex(p)+SUM(Nixd(p,:)))

       IF (ft(p) == 0.) THEN
          mwidth4  = sumd*(sumV/fc(p))*(qx(p))**2. *(Tex(p))**2.
       ELSEIF (fc(p) == 0.) THEN
          mwidth4  = ft(p)*sumd*sumV*(qx(p))**2. *(Tex(p))**2.
       ELSE
          mwidth4  = ft(p)*sumd*(sumV/fc(p))*(qx(p))**2. *(Tex(p))**2.
       ENDIF

       mwidth2 = SQRT(ABS(mwidth4))

       ! 4th order iterative expansion in theta to improve the s~0.5 behaviour
       sumd = ( SUM(Nixd(p,:)*((di(p,:)*di(p,:)+2.*d**2. *(smag(p)-alphax(p)-0.5) + &     
            2.*d**4./(mwidth2)*(smag(p)-2.*alphax(p)-0.25))*Joi2(:)*Jobani2(:))) + &	
            Nex(p)*((de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5)+ & 
            2.*d**4./(mwidth2)*(smag(p)-2.*alphax(p)-0.25))*Joe2*Jobane2))/ &
            (Nex(p)+SUM(Nixd(p,:)))

       IF (ft(p) == 0.) THEN
          mwidth4  = sumd*(sumV/fc(p))*(qx(p))**2. *(Tex(p))**2.
       ELSEIF (fc(p) == 0.) THEN
          mwidth4  = ft(p)*sumd*sumV*(qx(p))**2. *(Tex(p))**2.
       ELSE
          mwidth4  = ft(p)*sumd*(sumV/fc(p))*(qx(p))**2. *(Tex(p))**2.
       ENDIF

       mwidth2 = SQRT(ABS(mwidth4))

       sumd = ( SUM(Nixd(p,:)*((di(p,:)*di(p,:)+2.*d**2. *(smag(p)-alphax(p)-0.5) + &     
            2.*d**4/(mwidth2)*(smag(p)-2.*alphax(p)-0.25) + &
            5./40.*d**6./ABS(mwidth4)*(6.*smag(p)-32.*alphax(p)-1.))*Joi2(:)*Jobani2(:))) + &	
            Nex(p)*((de(p)*de(p)+2.*d**2. *(smag(p)-alphax(p)-0.5)+ & 
            2.*d**4./(mwidth2)*(smag(p)-2.*alphax(p)-0.25) + &
            5./40.*d**6./ABS(mwidth4)*(6.*smag(p)-32.*alphax(p)-1))*Joe2*Jobane2))/ &
            (Nex(p)+SUM(Nixd(p,:)))

       IF (ft(p) == 0.) THEN
          mwidth4  = sumd*(sumV/fc(p))*(qx(p))**2 *(Tex(p))**2
       ELSEIF (fc(p) == 0.) THEN
          mwidth4  = ft(p)*sumd*sumV*(qx(p))**2 *(Tex(p))**2
       ELSE
          mwidth4  = ft(p)*sumd*(sumV/fc(p))*(qx(p))**2 *(Tex(p))**2
       ENDIF
       mwidth=(ABS(mwidth4))**0.25

    ENDIF

  END SUBROUTINE old_calcwidth

  SUBROUTINE cot_fluidsol(p,nu)
    ! Calculates the fluid solution at radius p
    ! We maintain a 3rd order polynomial for flexibility if higher order is ever desired
    ! In practice, a 2nd order polynomial is actually solved
    INTEGER, INTENT(IN) :: p,nu
    !COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: omeflu
    COMPLEX(KIND=DBL), DIMENSION(4) :: poly
    COMPLEX(KIND=DBL), DIMENSION(4) :: polysol
    REAL(KIND=DBL), DIMENSION(2*3*4) :: W !work array from CPQR79
    REAL(KIND=DBL) :: wpe,wpi
    INTEGER :: ndeg=3,ifail
    INTEGER, DIMENSION(1) :: iloc    

    !**BUILD THE POLYNOMIAL COEFFICIENTS
    CALL cot_fonctflu (p, nu, poly)

    !**FIND ROOT OF POLYNOMIAL**
    CALL CPQR79(ndeg,poly,polysol,ifail,W)
    iloc = MAXLOC(AIMAG(polysol))
    omeflu = polysol(iloc(1))

    ! solflu equivalent to PC QLK_mom : omega_ITG=calculsolfluide(a_fluide,b_fluide,c_fluide)
    ! omega_fluide=cmplx(real(omega_ITG),aimag(omega_ITG)*coeffi) why coeffi in front of imaginary part only??

    IF (ifail /= 0) THEN
       IF (verbose .EQV. .TRUE.) WRITE(stderr,"(A,I0,A,I7,A,I3)") 'ifail = ',ifail,&
            &'. Abnormal termination of CPQR79 root finding in fluid solution at p=',p,', nu=',nu
    ENDIF


  END SUBROUTINE cot_fluidsol

  SUBROUTINE cot_fonctflu( p, nu, poly )  
    ! ----------------------------------------------------------------------
    ! Builds the polynomial whose roots provide the fluid limit growth rate
    ! Follows Pierre Cottier QLK_mom kinezero.f90 versions with rotation
    ! no ETG flag here at the moment
    ! ---------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: p, nu
    COMPLEX(KIND=DBL), DIMENSION(:), INTENT(OUT) :: poly

    ! Variables locales
    COMPLEX(KIND=DBL) :: Va, Vb, Vc
    REAL(KIND=DBL) :: Apeff, kparprim, gammaE_fluide, gamsn2

    INTEGER :: i

    Apeff = 0
    DO i = 1,nions
       IF (ion_type(p,i) .NE. 3) THEN
          Apeff= Apeff +tau(p,i)*Nix(p,i)*Zi(p,i)**2*(Ani(p,i)+Ati(p,i))
       ENDIF
    ENDDO

    gamsn2 = 1./Nex(p)*Apeff !P.C
    kparprim = abs(smag(p))/qx(p)
    gammaE_fluide=qx(p)*gammaE(p)/((smag(p)+eps)*sqrt(Tex(p)))

    ! new version with sign and cs/cref
    gammaE_fluide= qx(p)*gammaE(p)/((smag(p)+eps)*sqrt(Tex(p))) * cthi(p,1)/cref(p)


    !a_fluide in PC QLK_mom
    Va=((2-Ane(p))*(1+2.0*gammaE_fluide*(gammaE_fluide+Zeffx(p)*tau(p,1)*Machi(p,1)/sqrt(Tex(p))))+&
         2.*gammaE_fluide*Machi(p,1)/sqrt(Tex(p))*(Ane(p)-8)+cthi(p,1)/cref(p)*Aui(p,1)/sqrt(Tex(p))*(2.0*gammaE_fluide+Zeffx(p)*tau(p,1)*Machi(p,1)/sqrt(Tex(p)))+&
         CMPLX(0.0,0.5)*kparprim+gamsn2)/(1+2.*gammaE_fluide**2)

    !b_fluide in PC QLK_mom
    Vb=(0.5*(gammaE_fluide*(2-Ane(p))+cthi(p,1)/cref(p)*Aui(p,1)/sqrt(Tex(p))+Machi(p,1)/sqrt(Tex(p))*(Ane(p)-8.0))**2+&
         ft(p)/fc(p)*(Ane(p)+Ate(p))+gamsn2*(4.0-Ane(p)+ci*kparprim))/(1+2.*gammaE_fluide**2)

    !c_fluide in PC QLK_mom
    Vc=gamsn2*(ft(p)/fc(p)*(Ane(p)+Ate(p))+gamsn2*(2.0+CMPLX(0.,0.5)*kparprim))/(1+2.*gammaE_fluide**2)


    poly(1)     = Va
    poly(2)     = Vb
    poly(3)     = Vc
    poly(4)     = 0.

  END SUBROUTINE cot_fonctflu

  SUBROUTINE cot_calcwidth(ETGflag,p, nu) ! width calculation with rotation with gamma=consistent with fluidsolrot

    ! ************************************************************************************************************
    ! Follows Pierre Cottier QLK_mom kinezero.f90 version with rotation
    ! no ETG flag here at the moment
    !*************************************************************************************************************

    COMPLEX(kind=DBL) :: mwidth4, mwidth2, sumd, sumV 
    LOGICAL, INTENT(IN) :: ETGflag
    INTEGER, INTENT(IN) :: p, nu
    INTEGER :: i

    sumd =  1./Nex(p)*SUM(tau(p,:)*Nix(p,:)*Zi(p,:)**2.*(0.75*ft(p)/fc(p)*di(p,:)**2.+Rhoi(p,:)**2.*(1.5+8.0*(smag(p)-alphax(p)-0.5)/((kthetarhos(nu)*smag(p))**2*solflu(p,nu))))*Joi2(:)*Jobani2(:)) ! P.C.

    sumV =  1./Nex(p)*SUM(tau(p,:)*Nix(p,:)*Zi(p,:)**2*cthi(p,:)**2) ! P.C  
    !! Calcul de la largeur de mode et du decalage
    !! On calcule pour le moment façon garbet pas de changement de la largeur avec Aue, gammaE, Mach

    mwidth4 = (0.5*Rhoe(p)*cthe(p)*qx(p))**2.*(-(solflu(p,nu)/smag(p))**2*sumd/sumV) !P.C
    mwidth2 = SQRT((mwidth4))
    mwidth=(mwidth2)**0.5


  END SUBROUTINE cot_calcwidth

  SUBROUTINE cot_calcshift(ETGflag,p,nu) ! shift calculation hence with rotation 

    !************************************************
    ! Mode shift
    ! Follows Pierre Cottier QLK_mom kinezero.f90 version with rotation
    ! no ETG flag here at the moment
    !*************************************************************************************************************

    REAL(kind=DBL) :: sumV 
    LOGICAL, INTENT(IN) :: ETGflag
    INTEGER, INTENT(IN) :: p, nu
    INTEGER :: i
    REAL(KIND=DBL) :: gammaE_fluide

    sumV =  1./Nex(p)*SUM(tau(p,:)*Nix(p,:)*Zi(p,:)**2*cthi(p,:)**2) ! P.C  

    gammaE_fluide=qx(p)*gammaE(p)/((smag(p)+eps)*sqrt(Tex(p)))

    !fixed version
    gammaE_fluide= qx(p)*gammaE(p)/((smag(p)+eps)*sqrt(Tex(p))) * cthi(p,1)/cref(p)

    mshift = Rhoe(p)*cthe(p)*qx(p)*Jobani2(1)/(sqrt(sumV*smag(p)**2)*(solflu(p,nu)-Ane(p)))*&
         (gammaE_fluide*(2.0*solflu(p,nu)+2.0-Ane(p))&
         +cthi(p,1)/cref(p)*Aui(p,1)/sqrt(Tex(p))+Machi(p,1)/sqrt(Tex(p))*(Zeffx(p)*tau(p,1)*solflu(p,nu)+Ane(p)-8.0))

  END SUBROUTINE cot_calcshift

  SUBROUTINE ETGcheck(ETGflag,soll,kteta,p,nu)
    ! Checks the mode width regime and sets the ETG flag
    ! Criteria from Pierre Cottier PhD. For now, an alternative simple ktherarhos(nu) > 2 criteria is used (set elsewhere in the code)
    LOGICAL, INTENT(INOUT) :: ETGflag
    REAL(kind=DBL), INTENT(IN) :: kteta
    COMPLEX(kind=DBL), DIMENSION(numsols), INTENT(IN) :: soll 
    INTEGER, INTENT(IN) :: p, nu
    !index of the maximum growth rate (array due to MAXLOC)
    INTEGER, DIMENSION(1) :: iloc    

    iloc=MAXLOC(AIMAG(soll))

    IF ( (ETGflag .EQV. .FALSE.) .AND. kthetarhos(nu)>1.0 &
         .AND. REAL(soll(iloc(1)))>=0 .AND. kteta*smag(p)/(qx(p)*Ro(p)) &
         *SQRT(Tex(p)*2.*1.6d-16/9.1094d-31) &
         *REAL(mwidth)/(REAL(soll(iloc(1)))*nwg)<20.) THEN
       ETGflag = .TRUE.
    ELSEIF (ETGflag .EQV. .TRUE.) THEN
       IF (smag(p)<=0.5 .AND. kteta*smag(p)/(qx(p)*Ro(p)) &
            *SQRT(Tex(p)*2.*1.6d-16/9.1094d-31) &
            *REAL(mwidth)/(REAL(soll(iloc(1)))*nwg)>5. ) THEN
          ETGflag = .FALSE.
       ELSEIF (AIMAG(soll(iloc(1))) <= 0) THEN
          ETGflag = .FALSE.
       ENDIF
    ENDIF
  END SUBROUTINE ETGcheck

END MODULE mod_fluidsol
