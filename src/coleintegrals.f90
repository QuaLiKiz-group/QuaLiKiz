MODULE coleintegrals
  USE kind
  USE CUI
  USE Precision_Model
  use, intrinsic :: iso_c_binding
  IMPLICIT NONE
  
  
    
    
    
CONTAINS
SUBROUTINE coleint(ndim,aa,bb,minpts,maxpts,integrand,relacc,acc,lenwrk,wrkstr,COLEresultR, COLEresultI,ifailloc,choice,er)


    INTEGER, INTENT(IN) :: ndim,minpts,maxpts,lenwrk,choice
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: relacc,acc
    REAL(KIND=DBL) :: COLEresult
    REAL(KIND=DBL), INTENT(OUT) :: COLEresultR, COLEresultI
    REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: wrkstr
    COMPLEX(KIND=DBL), EXTERNAL :: integrand
    INTEGER :: i,j,recdepth, er

    REAL(KIND=DBL), DIMENSION(ndim) :: xx, xy, yx, yy, zx, zy
    REAL(KIND=DBL) :: count,rel
    INTEGER, PARAMETER :: m = 3
    INTEGER :: Key = 7
    
    INTEGER, DIMENSION (1:1) :: RgType
    INTEGER :: NEval
    REAL(kind=stnd), DIMENSION(1:2,0:2,1:m) :: Vertices
    REAL(kind=stnd), DIMENSION(1:1,0:1) :: Vertex
    REAL(kind=stnd), DIMENSION(1:2) :: IntegralValue, AbsErr
    REAL(kind=stnd) :: IntegralValue_1, AbsErr_1
    REAL(kind=stnd) :: EpsRel, EpsAbs
    LOGICAL :: Restart
    
    
    INTEGER :: ncomp = 1, userdata = 0, nvec = 1, flags = 0, seed = 0, mineval = 1, maxeval = 1000000
    CHARACTER*(100) :: statefile
    INTEGER*8 :: spin
    INTEGER :: nregions, fail
    REAL(KIND=DBL) :: error, prob
    
    real(c_double), DIMENSION(2) :: IntegralValueC, AbsErrC, xmin, xmax
    integer(c_size_t) :: MaxEvaluate = 10000000
    real(c_double) :: reqepsrel, reqepsabs
    type(c_funptr) :: C_func_ptr
    
    
    
    
    !INTEGER, EXTERNAL :: hcubature
    
    !REAL(KIND=DBL), DIMENSION(:), ALLOCATABLE :: work
    
    
    ! Define parameters for Cubpack 
    RgType(1) = 2
    Vertices(1:2,0,1) = (/0., 0./)
    Vertices(1:2,2,1) = (/0., 5./)
    
    Vertices(1:2,0,2) = (/2.2, 0./)
    Vertices(1:2,1,2) = (/2.5, 0./)
    Vertices(1:2,2,2) = (/2.2, 5./)
    
    Vertices(1:2,0,3) = (/2.5, 0./)
    Vertices(1:2,1,3) = (/5., 0./)
    Vertices(1:2,2,3) = (/2.5, 5./)
    
    Vertex(1,:) = (/0., 5./)
    
    epsrel = 0.1_stnd
    epsabs = 0._stnd
    Restart = .false.
    recdepth = 100
    
    ! For personal routines 
    rel = relacc
    xx(1) = 0.
    xx(2) = 0.
    xy(1) = 5.
    xy(2) = 5.
    yx(1) = 2.2
    yx(2) = 0.
    yy(1) = 5.
    yy(2) = 5.
    zx(1) = 2.5
    zx(2) = 0.
    zy(1) = 5.
    zy(2) = 5.
    
    
    SELECT CASE(5)
      CASE(1)
        !Cubpack low accuracy
        Vertices(1:2,1,1) = (/5., 0./)
        
        
        CALL CUBATR(2,2,Fkstarrstar_cubpack,1,Vertices(:,:,1:1),RgType,IntegralValue,AbsErr,NEval=NEval,EpsRel=epsrel,Restart=Restart, Key = 1, ifail = ifailloc)
        IF(ifailloc /= 0) er = ifailloc
        
        !CALL CUBATR(1,G,Vertex,1,IntegralValue_1,AbsErr_1,EpsRel=epsrel,NEval=NEval,JOB=1,KEY = choice)
        COLEresultR = IntegralValue(1)
        COLEresultI = IntegralValue(2) 
        
        
        
        
      CASE(2)
        !Cubpack high accuracy 
      
        epsrel = 0.001_stnd
        Vertices(1:2,1,1) = (/2.2, 0./)
        
        !Split up into 3 regions 

        CALL CUBATR(2,2,Fkstarrstar_cubpack,1,Vertices(:,:,1:1),RgType,IntegralValue,AbsErr,NEval=NEval,EpsRel=epsrel,Restart=Restart, ifail = ifailloc, maxpts = maxpts)
        COLEresultR = IntegralValue(1)
        COLEresultI = IntegralValue(2) 
        IF(ifailloc /= 0) er = ifailloc
        ifailloc = 1
        
        CALL CUBATR(2,2,Fkstarrstar_cubpack,1,Vertices(:,:,2:2),RgType,IntegralValue,AbsErr,NEval=NEval,EpsRel=epsrel,Restart=Restart, ifail = ifailloc, maxpts = maxpts)
        COLEresultR = COLEresultR + IntegralValue(1)
        COLEresultI = COLEresultI + IntegralValue(2)
        IF(ifailloc /= 0) er = ifailloc        
        ifailloc = 1
        
        CALL CUBATR(2,2,Fkstarrstar_cubpack,1,Vertices(:,:,3:3),RgType,IntegralValue,AbsErr,NEval=NEval,EpsRel=epsrel,Restart=Restart, ifail = ifailloc, maxpts = maxpts)
        COLEresultR = COLEresultR + IntegralValue(1)
        COLEresultI = COLEresultI + IntegralValue(2) 
        
        IF(ifailloc /= 0) er = ifailloc
        
        CASE(3)
          !Cuba low accuracy (too slow)
          CALL CUHRE(ndim, ncomp, H_cuba, userdata, nvec, epsrel, epsabs, flags, mineval, maxeval, Key, statefile, spin, nregions, neval, fail, COLEresult, error, prob)
          
        CASE(4)
          !Genz Hermite package (too slow)
          !CALL HRMSYM(2,1,1,10000000,HERM,epsabs, epsrel, 0, IntegralValue, AbsErr, NEval, fail, wrkstr)
          !COLEresult = IntegralValue(1) / 4. 
        CASE(5)
          
          !CALL hcubature(2_c_int, C_func_ptr, 2_c_int, xmin, xmax, maxevaluate, reqepsabs, reqepsrel, IntegralValueC, AbsErrC)
          !COLEresultR = IntegralValueC(1)
          !COleresultI = IntegralValueC(2) 
          
        CASE(6)
          !To test personal routines, the general call is CALL outer(ndim, aa, bb, integrand, inner, acc_k, acc_r, recdepth, ifailloc), where inner and outer are the inner and outer methods respectively
    END SELECT

    
    CONTAINS
    
    FUNCTION Fkstarrstar_cubpack(NUMFUN,X) RESULT(Value)
    USE Precision_Model
    INTEGER, INTENT(IN) :: NUMFUN
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
    REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
    REAL(kind=DBL), DIMENSION(2) :: xy
    COMPLEX(KIND=DBL) :: output
    xy = X
    output = integrand(2, xy, 1, 1)
    Value(1) = REAL(output)
    Value(2) = AIMAG(output)
    RETURN
    END FUNCTION Fkstarrstar_cubpack
    
    
    
    FUNCTION F_cubpack(NUMFUN,X) RESULT(Value)
    USE Precision_Model
    INTEGER, INTENT(IN) :: NUMFUN
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
    REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
    REAL(kind=DBL), DIMENSION(2) :: xy
    xy = X
    Value(1) = integrand(2,xy)
    RETURN
    END FUNCTION F_cubpack
    
    
    INTEGER FUNCTION H_cuba(ndim, x, ncomp, f)
      INTEGER :: ndim, ncomp
      REAL(KIND=DBL), DIMENSION(ndim) :: x, y
      REAL(KIND=DBL), DIMENSION(ncomp) :: f 
      !Cuba only integrates from 0 to 1, so the integral has to be scaled 
      y(1) = 2.2* x(1)
      y(2) = 5.*x(2)
      f(1) = 11.* integrand(2, y)
      H_cuba = 1
    END FUNCTION H_cuba
    
    SUBROUTINE HERM(NDIM, X, NF, FUNVLS)
      INTEGER, INTENT(IN):: NDIM, NF
      REAL(KIND=DBL), DIMENSION(NDIM), INTENT(IN) :: X
      REAL(KIND=DBL), DIMENSION(NF), INTENT(OUT) :: FUNVLS
      
      FUNVLS(1) = 2.*3.1415926535897932*integrand(ndim, X) * EXP( (X(1)**2 + X(2)**2)*.5)
    
    
    END SUBROUTINE HERM

  END SUBROUTINE coleint
  
  REAL(KIND=DBL) FUNCTION example(ndim, xx)
    INTEGER, INTENT(IN) :: ndim
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: xx
    example = SIN( (xx(1)**2.) * (xx(2)**3.)  )
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION ksimpsonsouter(ndim,aa,bb,integrand,inner,acc_k,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff,gg,hh,ii,jj,mm,nn
    REAL(KIND=DBL) :: h,fa,fb,fc,fd,fe,S,S2,iS,err1,err2,r,acc_k_new
	  REAL(KIND=DBL), DIMENSION(5) :: ss

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+cc)/2.
    ff = (bb+dd)/2.
	  gg = (aa+ee)/2.
    hh = (dd+ff)/2.
    ii = (ee+cc)/2.
    jj = (ff+bb)/2.
    fd = inner(ndim,gg,hh,integrand,acc_r,recdepth,ifailloc)
    fe = inner(ndim,ii,jj,integrand,acc_r,recdepth,ifailloc)
	  h = bb(1)-aa(1)
	  mm(2) = aa(2)
	  nn(2) = bb(2)
	
	  ss = (/ .9501, .2311, .6068, .4860, .8913 /)
	  ss = (ss * h) + aa(1)
    
    fa = inner(ndim,aa,dd,integrand,acc_r,recdepth,ifailloc)
    fb = inner(ndim,cc,bb,integrand,acc_r,recdepth,ifailloc)
    fc = inner(ndim,ee,ff,integrand,acc_r,recdepth,ifailloc)
	  fd = inner(ndim,gg,hh,integrand,acc_r,recdepth,ifailloc)
    fe = inner(ndim,ii,jj,integrand,acc_r,recdepth,ifailloc)
	
  	S = (h/6.)*(fa + 4.*fc + fb)
	  S2 = (h/12.)*(fa + 4.*fd + 2.* fc + 4.*fe + fb)
	  S = (16.*S2 - S)/15.
    
	  iS = h/8. * (fa + fc + fb)
	  DO i = 1,5
		  mm(1) = ss(i)
		  nn(1) = ss(i)
		  iS = iS + h/8. * inner(ndim,mm,nn,integrand,acc_r,recdepth,ifailloc)
	  END DO
	  err1 = ABS(S-iS)
	  err2 = ABS(S2-iS)
	  r = 1.
    IF(err2.NE.(0.)) THEN
      r = err1/err2
    END IF
    acc_k_new = acc_k
    IF((r.GT.(0.)).AND.(r.LT.(1.))) THEN
      acc_k_new = acc_k_new/r
    END IF
    IF (iS.EQ.(0.)) THEN
      iS = bb(1) - aa(1)
    END IF
    iS = abs(iS)
	
    ksimpsonsouter = ksimpsonsouteraux(ndim,aa,bb,integrand,inner,acc_k_new,acc_r,fa,fb,fc,iS,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION ksimpsonsouteraux(ndim,aa,bb,integrand,inner,acc_k,acc_r,fa,fb,fc,iS,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,fa,fb,fc,acc_r,iS
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff,gg,hh,ii,jj
    REAL(KIND=DBL) :: h,fd,fe,S,S2
    

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+cc)/2.
    ff = (bb+dd)/2.
    gg = (aa+ee)/2.
    hh = (dd+ff)/2.
    ii = (ee+cc)/2.
    jj = (ff+bb)/2.
    fd = inner(ndim,gg,hh,integrand,acc_r,recdepth,ifailloc)
    fe = inner(ndim,ii,jj,integrand,acc_r,recdepth,ifailloc)
    h = bb(1)-aa(1)
  	S = h/6. * (fa + 4.*fc + fb)
    S2 = (h/12.)*(fa + 4.*fd + 2.*fc + 4.*fe + fb)
  	S = (16.*S2 - S)/15.

    IF( ABS(S2-S).LE.(acc_k * iS))THEN
      output = S
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      output = ksimpsonsouteraux(ndim,aa,ff,integrand,inner,acc_k,acc_r,fa,fc,fd,iS,recdepth-1,ifailloc) + ksimpsonsouteraux(ndim,ee,bb,integrand,inner,acc_k,acc_r,fc,fb,fe,iS,recdepth-1,ifailloc)
    END IF
  END FUNCTION


  REAL(KIND=DBL) FUNCTION ksimpsonsouter_old(ndim,aa,bb,integrand,inner,acc_k,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff
    REAL(KIND=DBL) :: h,fa,fb,fc,S

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+cc)/2.
    ff = (bb+dd)/2.
    
    fa = inner(ndim,aa,dd,integrand,acc_r,recdepth,ifailloc)
    fb = inner(ndim,cc,bb,integrand,acc_r,recdepth,ifailloc)
    fc = inner(ndim,ee,ff,integrand,acc_r,recdepth,ifailloc)
    h = bb(1)-aa(1)
    S = (h/6.)*(fa + 4.*fc + fb)
    ksimpsonsouter_old = ksimpsonsouteraux_old(ndim,aa,bb,integrand,inner,acc_k,acc_r,S,fa,fb,fc,recdepth,ifailloc)
  END FUNCTION

  REAL(KIND=DBL) RECURSIVE FUNCTION ksimpsonsouteraux_old(ndim,aa,bb,integrand,inner,acc_k,acc_r,S,fa,fb,fc,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,S,fa,fb,fc,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff,gg,hh,ii,jj
    REAL(KIND=DBL) :: h,fd,fe,Sleft,Sright,S2,acc_k_new
    

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+cc)/2.
    ff = (bb+dd)/2.
    gg = (aa+ee)/2.
    hh = (dd+ff)/2.
    ii = (ee+cc)/2.
    jj = (ff+bb)/2.
    fd = inner(ndim,gg,hh,integrand,acc_r,recdepth,ifailloc)
    fe = inner(ndim,ii,jj,integrand,acc_r,recdepth,ifailloc)
    h = bb(1)-aa(1)
    Sleft = (h/12.)*(fa + 4.*fd + fc)
    Sright = (h/12.)*(fc + 4.*fe + fb)
    S2 = Sleft + Sright

    IF( (recdepth.LE.0).OR.(ABS(S2-S).LE.(15.*acc_k))) THEN
      output = S2 + (S2-S)/15.
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      acc_k_new = acc_k * 0.5
      output = ksimpsonsouteraux_old(ndim,aa,ff,integrand,inner,acc_k_new,acc_r,Sleft,fa,fc,fd,recdepth-1,ifailloc) + ksimpsonsouteraux_old(ndim,ee,bb,integrand,inner,acc_k_new,acc_r,Sright,fc,fb,fe,recdepth-1,ifailloc)
    END IF
  END FUNCTION
  
  
  REAL(KIND=DBL) FUNCTION rsimpsonsouter(ndim,aa,bb,integrand,inner,acc_k,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff
    REAL(KIND=DBL) :: h,fa,fb,fc,S

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+dd)/2.
    ff = (bb+cc)/2.
    
    fa = inner(ndim,aa,cc,integrand,acc_k,recdepth,ifailloc)
    fb = inner(ndim,dd,bb,integrand,acc_k,recdepth,ifailloc)
    fc = inner(ndim,ee,ff,integrand,acc_k,recdepth,ifailloc)
    h = bb(2)-aa(2)
    S = (h/6.)*(fa + 4.*fc + fb)
    rsimpsonsouter = rsimpsonsouteraux(ndim,aa,bb,integrand,inner,acc_k,acc_r,S,fa,fb,fc,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION rsimpsonsouteraux(ndim,aa,bb,integrand,inner,acc_k,acc_r,S,fa,fb,fc,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: S,fa,fb,fc,acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,ff,gg,hh,ii,jj
    REAL(KIND=DBL) :: h,fd,fe,Sleft,Sright,S2,acc_r_new
    

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    ee = (aa+dd)/2.
    ff = (bb+cc)/2.
    gg = (aa+ee)/2.
    hh = (cc+ff)/2.
    ii = (ee+dd)/2.
    jj = (ff+bb)/2.
    fd = inner(ndim,gg,hh,integrand,acc_k,recdepth,ifailloc)
    fe = inner(ndim,ii,jj,integrand,acc_k,recdepth,ifailloc)
    h = bb(1)-aa(1)
    Sleft = (h/12.)*(fa + 4.*fd + fc)
    Sright = (h/12.)*(fc + 4.*fe + fb)
    S2 = Sleft + Sright

    IF( (recdepth.LE.0).OR.(ABS(S2-S).LE.(15.*acc_r))) THEN
      output = S2 + (S2-S)/15.
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      acc_r_new = acc_r * 0.5
      output = rsimpsonsouteraux(ndim,aa,ff,integrand,inner,acc_k,acc_r_new,Sleft,fa,fc,fd,recdepth-1,ifailloc) + rsimpsonsouteraux(ndim,ee,bb,integrand,inner,acc_k,acc_r_new,Sright,fc,fb,fe,recdepth-1,ifailloc)
    END IF
  END FUNCTION
  

  REAL(KIND=DBL) FUNCTION ggouter(ndim,aa,bb,integrand,inner,acc_k,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd
    REAL(KIND=DBL), DIMENSION(13) :: yy
    REAL(KIND=DBL), DIMENSION(12) :: xx
    REAL(KIND=DBL) :: h,fa,fb,fc,alpha,beta,i1,i2,iS,err1,err2,r,acc_k_new,m
    
    alpha = SQRT(2._DBL/3._DBL); beta = 1._DBL/SQRT(5._DBL);
    xx = (/ 0., -0.942882415695480, -alpha, -0.641853342345781,-beta, -0.236383199662150, 0., 0.236383199662150, beta, 0.641853342345781, alpha, 0.942882415695480 /)
    !xx = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
    
    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    m = 0.5*(aa(1)+bb(1))
    h = 0.5*(bb(1)-aa(1))
    fa = inner(ndim,aa,dd,integrand,acc_r,recdepth,ifailloc)
    fb = inner(ndim,cc,bb,integrand,acc_r,recdepth,ifailloc)
    yy(1) = fa; yy(13) = fb;
    DO i = 2,12
      cc(1) = m + h * xx(i)
      dd(1) = cc(1)
      yy(i) = inner(ndim,cc,dd,integrand,acc_r,recdepth,ifailloc)
    END DO
    i2 = (h/6.) * (yy(1) + yy(13) + 5.*(yy(5) + yy(9)))
    i1 = (h/1470.)*(77. * (yy(1) + yy(13)) + 432. * (yy(3) + yy(11)) + 625. * (yy(5) + yy(9)) + 672. * yy(7)) 
    iS = h*(0.0158271919734802*(yy(1)+yy(13))+ 0.0942738402188500* (yy(2)+yy(12))+0.155071987336585*(yy(3)+yy(11))+ 0.188821573960182*(yy(4)+yy(10)) &
    & + 0.199773405226859* (yy(5)+yy(9))+0.224926465333340*(yy(6)+yy(8))+ 0.242611071901408*yy(7))
    err1 = ABS(i1-is)
    err2 = ABS(i2-is)
    r = 1.
    IF(err2.NE.(0.)) THEN
      r = err1/err2
    END IF
    acc_k_new = acc_k
    IF((r.GT.(0.)).AND.(r.LT.(1.))) THEN
      acc_k_new = acc_k_new/r
    END IF
    IF (iS.EQ.(0.)) THEN
      iS = bb(1) - aa(1)
    END IF
    iS = abs(iS)
    ggouter = ggouteraux(ndim,aa,bb,integrand,inner,acc_k_new,acc_r,iS,fa,fb,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION ggouteraux(ndim,aa,bb,integrand,inner,acc_k,acc_r,iS,fa,fb,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r,iS
    REAL(KIND=DBL), EXTERNAL :: integrand,inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,aamll,bbmll,aaml,bbml,aam,bbm,aamr,bbmr,aamrr,bbmrr
    REAL(KIND=DBL) :: h,fa,fb,fc,alpha,beta,i1,i2,mll,ml,m,mr,mrr,fmll,fml,fm,fmr,fmrr
    alpha = SQRT(2./3.); beta = 1./SQRT(5.);

    cc(1)=bb(1); cc(2)=aa(2);
    dd(1)=aa(1); dd(2)=bb(2);
    m = 0.5*(aa(1)+bb(1))
    h = 0.5*(bb(1)-aa(1))
	
	mll= m - alpha * h
	ml = m - beta * h
	mr = m + beta * h
	mrr = m + alpha*h 
	aamll(1) = mll; aamll(2) = cc(2); bbmll(1) = mll; bbmll(2) = dd(2);
	aaml(1) = ml; aaml(2) = cc(2); bbml(1) = ml; bbml(2) = dd(2);
	aam(1) = m; aam(2) = cc(2); bbm(1) = m; bbm(2) = dd(2);
	aamr(1) = mr; aamr(2) = cc(2); bbmr(1) = mr; bbmr(2) = dd(2);
	aamrr(1) = mrr; aamrr(2) = cc(2); bbmrr(1) = mrr; bbmrr(2) = dd(2);
	
	fmll = inner(ndim,aamll,bbmll,integrand,acc_r,recdepth,ifailloc);
	fml = inner(ndim,aaml,bbml,integrand,acc_r,recdepth,ifailloc);
	fm = inner(ndim,aam,bbm,integrand,acc_r,recdepth,ifailloc);
	fmr = inner(ndim,aamr,bbmr,integrand,acc_r,recdepth,ifailloc);
	fmrr = inner(ndim,aamrr,bbmrr,integrand,acc_r,recdepth,ifailloc);
	i2 = h/6. * (fa + fb + 5. * (fml + fmr))
	i1 = h/1470.*(77.*(fa+fb)+432.*(fmll+fmrr)+625.*(fml+fmr)+672.*fm)
	output = 0
	IF(ABS(i1-i2).LE.(acc_r*iS)) THEN
		output = i1
	ELSE
		output = ggouteraux(ndim,aa,bbmll,integrand,inner,acc_k,acc_r,iS,fa,fmll,recdepth,ifailloc) + ggouteraux(ndim,aamll,bbml,integrand,inner,acc_k,acc_r,iS,fmll,fml,recdepth,ifailloc) &
		& + ggouteraux(ndim,aaml,bbm,integrand,inner,acc_k,acc_r,iS,fml,fm,recdepth,ifailloc) + ggouteraux(ndim,aam,bbmr,integrand,inner,acc_k,acc_r,iS,fm,fmr,recdepth,ifailloc) & 
		& + ggouteraux(ndim,aamr,bbmrr,integrand,inner,acc_k,acc_r,iS,fmr,fmrr,recdepth,ifailloc) + ggouteraux(ndim,aamrr,bb,integrand,inner,acc_k,acc_r,iS,fmrr,fb,recdepth,ifailloc)
	END IF
	
	
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION gginner(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc
    REAL(KIND=DBL), DIMENSION(13) :: yy
    REAL(KIND=DBL), DIMENSION(12) :: xx
    REAL(KIND=DBL) :: h,fa,fb,fc,alpha,beta,i1,i2,iS,err1,err2,r,acc_r_new,m
    
    alpha = SQRT(2./3.); beta = 1./SQRT(5.);
    xx = (/ 0., -0.942882415695480, -alpha, -0.641853342345781,-beta, -0.236383199662150, 0., 0.236383199662150, beta, 0.641853342345781, alpha, 0.942882415695480 /)
    !xx = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)

	cc(1) = aa(1)
    m = 0.5*(aa(2)+bb(2))
    h = 0.5*(bb(2)-aa(2))
    fa = integrand(ndim, aa)
    fb = integrand(ndim, bb)
    yy(1) = fa; yy(13) = fb;
    DO i = 2,12
      cc(2) = m + h * xx(i)
      yy(i) = integrand(ndim, cc)
    END DO
    i2 = (h/6.) * (yy(1) + yy(13) + 5.*(yy(5) + yy(9)))
    i1 = (h/1470.)*(77. * (yy(1) + yy(13)) + 432. * (yy(3) + yy(11)) + 625. * (yy(5) + yy(9)) + 672. * yy(7)) 
    iS = h*(0.0158271919734802*(yy(1)+yy(13))+ 0.0942738402188500* (yy(2)+yy(12))+0.155071987336585*(yy(3)+yy(11))+ 0.188821573960182*(yy(4)+yy(10)) &
    & + 0.199773405226859* (yy(5)+yy(9))+0.224926465333340*(yy(6)+yy(8))+ 0.242611071901408*yy(7))
    err1 = ABS(i1-is)
    err2 = ABS(i2-is)
    r = 1.
    IF(err2.NE.(0.)) THEN
      r = err1/err2
    END IF
    acc_r_new = acc_r
    IF((r.GT.(0.)).AND.(r.LT.(1.))) THEN
      acc_r_new = acc_r_new / r
    END IF
    IF (iS.EQ.(0.)) THEN
      iS = bb(2) - aa(2)
    END IF
    iS = abs(iS)
    gginner = gginneraux(ndim,aa,bb,integrand,acc_r_new,iS,fa,fb,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION gginneraux(ndim,aa,bb,integrand,acc_r,iS,fa,fb,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r,iS
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: aamll,aaml,aam,aamr,aamrr
    REAL(KIND=DBL) :: h,fa,fb,fc,alpha,beta,i1,i2,mll,ml,m,mr,mrr,fmll,fml,fm,fmr,fmrr
    alpha = SQRT(2./3.); beta = 1./SQRT(5.);

    m = 0.5*(aa(2)+bb(2))
    h = 0.5*(bb(2)-aa(2))
	
	mll= m - alpha * h
	ml = m - beta * h
	mr = m + beta * h
	mrr = m + alpha*h 
	aamll(2) = mll; aamll(1) = aa(1)
	aaml(2) = ml; aaml(1) = aa(1)
	aam(2) = m; aam(1) = aa(1); 
	aamr(2) = mr; aamr(1) = aa(1); 
	aamrr(2) = mrr; aamrr(1) = aa(1);
	
	fmll = integrand(ndim, aamll)
	fml = integrand(ndim,aaml);
	fm = integrand(ndim,aam);
	fmr = integrand(ndim,aamr);
	fmrr = integrand(ndim,aamrr);
	i2 = h/6. * (fa + fb + 5. * (fml + fmr))
	i1 = h/1470.*(77.*(fa+fb)+432.*(fmll+fmrr)+625.*(fml+fmr)+672.*fm)
	output = 0
	IF(ABS(i1-i2).LE.(acc_r*iS)) THEN
		output = i1
	ELSE
		output = gginneraux(ndim,aa,aamll,integrand,acc_r,iS,fa,fmll,recdepth,ifailloc) + gginneraux(ndim,aamll,aaml,integrand,acc_r,iS,fmll,fml,recdepth,ifailloc) &
		& + gginneraux(ndim,aaml,aam,integrand,acc_r,iS,fml,fm,recdepth,ifailloc) + gginneraux(ndim,aam,aamr,integrand,acc_r,iS,fm,fmr,recdepth,ifailloc) & 
		& + gginneraux(ndim,aamr,aamrr,integrand,acc_r,iS,fmr,fmrr,recdepth,ifailloc) + gginneraux(ndim,aamrr,bb,integrand,acc_r,iS,fmrr,fb,recdepth,ifailloc)
	END IF
	
	
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION rsimpsonsinner(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee,mm
    REAL(KIND=DBL) :: h,fa,fb,fc,fd,fe,S,S2,iS,acc_r_new,r,err1,err2
  	REAL(KIND=DBL), DIMENSION(5) :: ss
	
    cc = (aa+bb)/2.
  	dd = (aa+cc)/2.
    ee = (bb+cc)/2.
  	h = bb(2)-aa(2)
    fa = integrand(ndim,aa)
    fb = integrand(ndim,bb)
    fc = integrand(ndim,cc)
  	fd = integrand(ndim,dd)
    fe = integrand(ndim,ee)
  	mm(1) = aa(1)
	
  	ss = (/ .9501, .2311, .6068, .4860, .8913 /)
  	ss = (ss * h) + aa(2)
	
	
    S = (h/6.)*(fa + 4.*fc + fb)
  	S2 = (h/12.)*(fa + 4.*fd + 2.*fc + 4.*fe + fb)
  	S = (16.*S2 - S)/15.
	
  	iS = h/8. * (fa + fc + fb)
  	DO i = 1,5
		mm(2) = ss(i) 
		iS = iS + h/8. * integrand(ndim, ss)
  	END DO
  	err1 = ABS(S-iS)
  	err2 = ABS(S2-iS)
  	r = 1.
    IF(err2.NE.(0.)) THEN
      r = err1/err2
    END IF
    acc_r_new = acc_r
    IF((r.GT.(0.)).AND.(r.LT.(1.))) THEN
      acc_r_new = acc_r_new/r
    END IF
    IF (iS.EQ.(0.)) THEN
      iS = bb(2) - aa(2)
    END IF
    iS = abs(iS)
	
    rsimpsonsinner = rsimpsonsinneraux(ndim,aa,bb,integrand,acc_r,fa,fb,fc,iS,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION rsimpsonsinneraux(ndim,aa,bb,integrand,acc_r,fa,fb,fc,iS,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r,fa,fb,fc,iS
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee
    REAL(KIND=DBL) :: h,fd,fe,S,Sleft,Sright,S2

    cc = (aa+bb)/2.
    dd = (aa+cc)/2.
    ee = (bb+cc)/2.
    fd = integrand(ndim,dd)
    fe = integrand(ndim,ee)
    h = bb(2)-aa(2)
  	S = h/6. * (fa + 4.*fc + fb)
    S2 = (h/12.)*(fa + 4.*fd + 2.*fc + 4.*fe + fb)
  	S = (16.*S2 - S)/15.

    IF( ABS(S2-S).LE.(acc_r * iS)) THEN
      output = S
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      output = rsimpsonsinneraux(ndim,aa,cc,integrand,acc_r,fa,fc,fd,iS,recdepth-1,ifailloc) + rsimpsonsinneraux(ndim,cc,bb,integrand,acc_r,fc,fb,fe,iS,recdepth-1,ifailloc)
    END IF
  END FUNCTION


  REAL(KIND=DBL) FUNCTION rsimpsonsinner_old(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc
    REAL(KIND=DBL) :: h,fa,fb,fc,S

    cc = (aa+bb)/2.
    fa = integrand(ndim,aa)
    fb = integrand(ndim,bb)
    fc = integrand(ndim,cc)
    h = bb(2)-aa(2)
    S = (h/6.)*(fa + 4.*fc + fb)
    rsimpsonsinner_old = rsimpsonsinneraux_old(ndim,aa,bb,integrand,acc_r,S,fa,fb,fc,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION rsimpsonsinneraux_old(ndim,aa,bb,integrand,acc_r,S,fa,fb,fc,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r,S,fa,fb,fc
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee
    REAL(KIND=DBL) :: h,fd,fe,Sleft,Sright,S2,acc_r_new

    cc = (aa+bb)/2.
    dd = (aa+cc)/2.
    ee = (bb+cc)/2.
    fd = integrand(ndim,dd)
    fe = integrand(ndim,ee)
    h = bb(2)-aa(2)
    Sleft = (h/12.)*(fa + 4.*fd + fc)
    Sright = (h/12.)*(fc + 4.*fe + fb)
    S2 = Sleft + Sright

    IF( (recdepth.LE.0).OR.(ABS(S2-S).LE.(15.*acc_r)) ) THEN
      output = S2 + (S2-S)/15.
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      acc_r_new = acc_r * 0.5
      output = rsimpsonsinneraux_old(ndim,aa,cc,integrand,acc_r_new,Sleft,fa,fc,fd,recdepth-1,ifailloc) + rsimpsonsinneraux_old(ndim,cc,bb,integrand,acc_r_new,Sright,fc,fb,fe,recdepth-1,ifailloc)
    END IF
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION ksimpsonsinner(ndim,aa,bb,integrand,acc_k,recdepth,ifailloc)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc
    REAL(KIND=DBL) :: h,fa,fb,fc,S

    cc = (aa+bb)/2.
    fa = integrand(ndim,aa)
    fb = integrand(ndim,bb)
    fc = integrand(ndim,cc)
    h = bb(1)-aa(1)
    S = (h/6.)*(fa + 4.*fc + fb)
    ksimpsonsinner = ksimpsonsinneraux(ndim,aa,bb,integrand,acc_k,S,fa,fb,fc,recdepth,ifailloc)
  END FUNCTION
  
  REAL(KIND=DBL) RECURSIVE FUNCTION ksimpsonsinneraux(ndim,aa,bb,integrand,acc_k,S,fa,fb,fc,recdepth,ifailloc) RESULT(output)
    INTEGER, INTENT(IN) :: ndim,recdepth
    INTEGER :: ifailloc
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,S,fa,fb,fc
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd,ee
    REAL(KIND=DBL) :: h,fd,fe,Sleft,Sright,S2,acc_k_new

    cc = (aa+bb)/2.
    dd = (aa+cc)/2.
    ee = (bb+cc)/2.
    fd = integrand(ndim,dd)
    fe = integrand(ndim,ee)
    h = bb(1)-aa(1)
    Sleft = (h/12.)*(fa + 4.*fd + fc)
    Sright = (h/12.)*(fc + 4.*fe + fb)
    S2 = Sleft + Sright

    IF( (recdepth.LE.0).OR.(ABS(S2-S).LE.(15.*acc_k)) ) THEN
      output = S2 + (S2-S)/15.
      IF(recdepth.LE.0) ifailloc = 50
    ELSE
      acc_k_new = acc_k * 0.5
      output = ksimpsonsinneraux(ndim,aa,cc,integrand,acc_k_new,Sleft,fa,fc,fd,recdepth-1,ifailloc) + ksimpsonsinneraux(ndim,cc,bb,integrand,acc_k_new,Sright,fc,fb,fe,recdepth-1,ifailloc)
    END IF
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION hermiteinner20(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN):: ndim,recdepth
    INTEGER :: ifailloc, n, maxn, initn, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,oo
    REAL(KIND=DBL), DIMENSION(10) :: r
    REAL(KIND=DBL), DIMENSION(10) :: weights
    
    oo(1) = aa(1)
    r(1) = 5.387480890011
    r(2) = 4.603682449550
    r(3) = 3.944764040115
    r(4) = 3.347854567383
    r(5) = 2.788806058428
    r(6) = 2.254974002089
    r(7) = 1.738537712116
    r(8) = 1.234076215395
    r(9) = 0.737473728545
    r(10) = 0.245340708300
    weights(1) = 0.8985919614531
    weights(2) = 0.7043329611769
    weights(3) = 0.6222786961914
    weights(4) = 0.5752624428525
    weights(5) = 0.5448517423645
    weights(6) = 0.5240803509485
    weights(7) = 0.5096790271174
    weights(8) = 0.4999208713362
    weights(9) = 0.4938433852720
    weights(10) = 0.4909215006667
    r = r * sqrt(2._DBL)
    weights = weights * sqrt(2._DBL)
    hermiteinner20 = 0
    DO i = 1,10
      oo(2) = r(i)
      hermiteinner20 = hermiteinner20 + weights(i) * integrand(ndim, oo)
    END DO
    
    
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION hermiteinner20_alone(k, integrand)
    REAL(KIND=DBL), INTENT(IN) :: k
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(2) :: oo
    REAL(KIND=DBL), DIMENSION(10) :: r
    REAL(KIND=DBL), DIMENSION(10) :: weights
    INTEGER :: i
    
    oo(1) = k
    r(1) = 5.387480890011
    r(2) = 4.603682449550
    r(3) = 3.944764040115
    r(4) = 3.347854567383
    r(5) = 2.788806058428
    r(6) = 2.254974002089
    r(7) = 1.738537712116
    r(8) = 1.234076215395
    r(9) = 0.737473728545
    r(10) = 0.245340708300
    weights(1) = 0.8985919614531
    weights(2) = 0.7043329611769
    weights(3) = 0.6222786961914
    weights(4) = 0.5752624428525
    weights(5) = 0.5448517423645
    weights(6) = 0.5240803509485
    weights(7) = 0.5096790271174
    weights(8) = 0.4999208713362
    weights(9) = 0.4938433852720
    weights(10) = 0.4909215006667
    r = r * sqrt(2._DBL)
    weights = weights * sqrt(2._DBL)
    hermiteinner20_alone = 0
    DO i = 1,10
      oo(2) = r(i)
      hermiteinner20_alone = hermiteinner20_alone + weights(i) * integrand(2, oo)
    END DO
    
    
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION hermiteinner10(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN):: ndim,recdepth
    INTEGER :: ifailloc, n, maxn, initn, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,oo
    REAL(KIND=DBL), DIMENSION(5) :: r
    REAL(KIND=DBL), DIMENSION(5) :: weights
    
    oo(1) = aa(1)
    r(1) = 3.4361591188377
    r(2) = 2.5327316742328
    r(3) = 1.7566836492999
    r(4) = 1.0366108297895
    r(5) = 0.34290132722371
    weights(1) = 1.025451691366
    weights(2) = 0.8206661264048
    weights(3) = 0.7414419319436
    weights(4) = 0.7032963231049
    weights(5) = 0.6870818539513
    r = r * sqrt(2._DBL)
    weights = weights * sqrt(2._DBL)
    hermiteinner10 = 0
    DO i = 1,5
      oo(2) = r(i)
      hermiteinner10 = hermiteinner10 + weights(i) * integrand(ndim, oo)
    END DO
    
    
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION hermiteouter20(ndim,aa,bb,integrand,inner,acc_k,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN):: ndim,recdepth
    INTEGER :: ifailloc, n, maxn, initn, i
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_k,acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand, inner
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,dd
    REAL(KIND=DBL), DIMENSION(10) :: r
    REAL(KIND=DBL), DIMENSION(10) :: weights
    
    cc(1) = aa(1)
    dd(1) = bb(1)
    r(1) = 5.387480890011
    r(2) = 4.603682449550
    r(3) = 3.944764040115
    r(4) = 3.347854567383
    r(5) = 2.788806058428
    r(6) = 2.254974002089
    r(7) = 1.738537712116
    r(8) = 1.234076215395
    r(9) = 0.737473728545
    r(10) = 0.245340708300
    weights(1) = 0.8985919614531
    weights(2) = 0.7043329611769
    weights(3) = 0.6222786961914
    weights(4) = 0.5752624428525
    weights(5) = 0.5448517423645
    weights(6) = 0.5240803509485
    weights(7) = 0.5096790271174
    weights(8) = 0.4999208713362
    weights(9) = 0.4938433852720
    weights(10) = 0.4909215006667
    r = r * sqrt(2._DBL)
    weights = weights * sqrt(2._DBL)
    hermiteouter20 = 0
    DO i = 1,10
      cc(2) = r(i); dd(2) = r(i);
      hermiteouter20 = hermiteouter20 + weights(i) * inner(ndim,cc,dd,integrand,acc_k,recdepth,ifailloc)
    END DO
    
    
  END FUNCTION
  
  
  REAL(KIND=DBL) FUNCTION hermiteinnerext(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN):: ndim,recdepth
    INTEGER :: ifailloc, i, maxn, initn
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: oo
    REAL(KIND=DBL) :: S,S2
    REAL(KIND=DBL), DIMENSION(2) :: xx2, ww2, yy2
    REAL(KIND=DBL), DIMENSION(5) :: xx5, ww5, yy5
    REAL(KIND=DBL), DIMENSION(10) :: xx10, ww10, yy10
    REAL(KIND=DBL), DIMENSION(20) :: xx20, ww20, yy20
      
    xx2 = (/0.52464762327529, 1.6506801238858/)
    ww2 = (/1.059964482895, 1.2402258176958/)
    
    xx5 = (/0., 0.8158264554988470, 1.6506801238858, 1.7, 2.652723696114967/)
    ww5 = (/0.809914950277796, 0.827928649625335, 0.4788352137519195, 0.4205114459295698, 1.00593374178376/)
    
    xx10 = (/0., 0.3, 0.8158264554988470, 1.317612483964416, 1.7, 2.100956366406892, 2.652723696114967, 3.279206771679270, 3.982643254881533, 4.832090233693644/)
    ww10 = (/0.1293780139310889, 0.4893182544521739, 0.5217646793250587, 0.4614124724099871, 0.3299610906477909, 0.4940158544192529, 0.5940761304690768, &
    & 0.6600004927947552, 0.7564082267704276, 0.985833561248949/)
    
    xx20 = (/0., 0.3, 0.5641271694337688, 0.8158264554988470, 1., 1.317612483964416, 1.7, 2.100956366406892, 2.504818747315663, 2.652723696114967, &
    & 2.951332304548583, 3.383149961801906, 3.807740919396219, 3.982643254881533, 4.338842735586137, 4.832090233693644, 5.369129570123683, & 
    & 5.967030999404972, 6.682099264478276, 9.516210149246421 /)
    
    ww20 = (/0.3130815492784033, 0.2774323356488918, 0.2651338792948922, 0.2013819579824516, 0.2384111522964433, 0.3654055159571945, 0.3941670233037691, &
    & 0.4062731579574491, 0.3804406527191944, 0.05767833514377005, 0.4219121221612513, 0.4377179799775914, 0.3655058466673398, 0.1416241714948122, &
    &  0.4687084925536405, 0.5143383540855813, 0.5626012062853013, 0.6410618467316369, 0.821178613385312, 11294413.8460669/)
      
    oo(1) = aa(1)
    
    S = 0.
    S2 = 0.
    ww5(1) = ww5(1)/2.
    ww10(1) = ww10(1)/2.
    ww20(1) = ww20(1)/2.
    
    DO i = 1,2
      oo(2) = xx2(i) * SQRT(2.)
      yy2(i) = integrand(ndim, oo)
      S = S + yy2(i) * SQRT(2.) * ww2(i)
    END DO
    
    DO i = 1,5
      IF (i.EQ.3) THEN
        yy5(3) = yy2(2)
      ELSE
        oo(2) = xx5(i) * SQRT(2.)
        yy5(i) = integrand(ndim, oo)
      END IF
      S2 = S2 + yy5(i) * SQRT(2.) * ww5(i)
    END DO
	
    IF (ABS(S2-S).LT.(acc_r * S2)) THEN
      hermiteinnerext = S2
    ELSE
      S = S2
      S2 = 0.
      DO i = 1,10
        IF(i.EQ.1) THEN
          yy10(1) = yy5(1)
        ELSE IF((i.EQ.3).OR.(i.EQ.5)) THEN
          yy10(i) = yy5(i-1)
        ELSE IF(i.EQ.7) THEN
          yy10(7) = yy5(5)
        ELSE
          oo(2) = xx10(i) * SQRT(2.)
          yy10(i) = integrand(ndim, oo)
        END IF
        S2 = S2 + yy10(i) * SQRT(2.) * ww10(i)			
      END DO
      IF (ABS(S2-S).LT.(acc_r * S2)) THEN
        hermiteinnerext = S2
      ELSE
        S = S2
        S2 = 0.
        DO i = 1,20
          IF((i.EQ.1).OR.(i.EQ.2)) THEN
            yy20(i) = yy10(i)
          ELSE IF((i.GE.5).AND.(i.LE.8)) THEN
            yy20(i) = yy10(i-2)
          ELSE IF(i.EQ.10) THEN
            yy20(10) = yy10(7)
          ELSE IF(i.EQ.14) THEN
            yy20(14) = yy10(9)
          ELSE IF(i.EQ.16) THEN
            yy20(16) = yy10(10)
          ELSE
            oo(2) = xx20(i) * SQRT(2.)
            yy20(i) = integrand(ndim, oo)
            
          END IF
          S2 = S2 + yy20(i) * SQRT(2.) * ww20(i)
        END DO
        IF(ABS(S2-S).LT.(acc_r*S2)) THEN
          hermiteinnerext = S2
        ELSE
          hermiteinnerext = hermiteinner(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
        END IF
      END IF
    END IF
    
  END FUNCTION
  
  REAL(KIND=DBL) FUNCTION hermiteinner(ndim,aa,bb,integrand,acc_r,recdepth,ifailloc)
    INTEGER, INTENT(IN):: ndim,recdepth
    INTEGER :: ifailloc, n, maxn, initn
    REAL(KIND=DBL), DIMENSION(ndim), INTENT(IN) :: aa,bb
    REAL(KIND=DBL), INTENT(IN) :: acc_r
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL), DIMENSION(ndim) :: cc,oo
    REAL(KIND=DBL) :: a, f0, h1
    
    oo(1) = aa(1)
    oo(2) = 0.
    initn = 41 !the min order computed 
    maxn = 61!the max order computed
    f0 = integrand(ndim,oo)
    n = initn
    h1 = hermiten(ndim,oo,integrand,n-2,f0)
    hermiteinner = hermiten(ndim,oo,integrand,n,f0)
    n = n + 2
    DO WHILE((ABS(h1-hermiteinner).GT.(ABS(acc_r*hermiteinner))).AND.(n.LE.maxn))
      h1 = hermiteinner
      hermiteinner = hermiten(ndim, oo, integrand, n, f0)
      n = n+2
    END DO
  END FUNCTION

  REAL(KIND=DBL) FUNCTION hermiten(ndim,oo,integrand,n,f0)
    INTEGER, INTENT(IN):: ndim,n
    REAL(KIND=DBL), INTENT(IN):: f0
    INTEGER :: maxx,m,i,j,its
    REAL(KIND=DBL), DIMENSION(ndim) :: oo
    REAL(KIND=DBL), DIMENSION(n) :: xx, ww
    REAL(KIND=DBL), EXTERNAL :: integrand
    REAL(KIND=DBL) :: acc,pi4,z,nd,jd,p1,p2,p3,pp,z1
  
    acc = 10.**(-14)
    pi4 = 0.7511255444649425
    maxx = 10
    nd = 1.*n
    m = (n+1)/2
    z = 0.
    !first, we compute the weights and abscissa
    DO i = 1,m
      IF(i.EQ.1) THEN 
        z = SQRT(2.*nd+1.) - 1.85575 * (2.*nd+1.)**(-0.16667)
      ELSE IF (i.EQ.2) THEN
        z = z -1.14 * nd**(0.426) / z
      ELSE IF (i.EQ.3) THEN
        z = 1.86*z -0.86 * xx(1)
      ELSE IF (i.EQ.4) THEN
        z = 1.91 * z -0.91 * xx(2)
      ELSE
        z=2.*z - xx(i-2)
      END IF
     
      DO its = 1,maxx
        p1 = pi4
        p2 = 0.
        DO j = 1,n
          jd = 1.*j
          p3 = p2
          p2 = p1
          p1 = z * SQRT(2./jd) * p2 - SQRT((jd-1.)/jd) * p3
        END DO
        pp = SQRT(2.*nd) * p2
        z1 = z
        z = z1 - p1/pp
        IF(ABS(z-z1).LE.acc) THEN
          EXIT
        END IF
      END DO
      xx(i) = z
      xx(n+1-i) = -z
      ww(i) = 2./(pp*pp)
      ww(n+1-i) = ww(i)
    END DO
    !now, we actually compute the integral
    !note that we assume the integrand is of the form exp(-x^2/2) f(x), and that we integrate from 0 to infinity
    !this will introduce factors of sqrt(2) 
    !also note that here we have asumed that n is odd; for the integral in question, it is always better to use odd than even
    !this is because we can reuse f(0) over and over when calculating the integral with different values of n
    hermiten = 0.5 * SQRT(2._DBL) * ww(m) * f0 !0.5 comes from the fact that we split the integral in half
    DO i = 1,m-1
      oo(2) = SQRT(2._DBL) * xx(i)
      hermiten = hermiten + SQRT(2._DBL) * ww(i) * EXP(xx(i)**2.) * integrand(ndim, oo)
    END DO



  END FUNCTION
  



END MODULE coleintegrals
