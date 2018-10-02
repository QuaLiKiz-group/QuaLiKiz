MODULE HCUB
  USE KIND
  IMPLICIT NONE

!Define types
  
TYPE ESTERR
  REAL(KIND=DBL) :: val_
  REAL(KIND=DBL) :: err_
END TYPE

TYPE HYPERCUBE
  INTEGER :: ndim
  REAL(KIND=DBL), DIMENSION(:), POINTER :: dat !size is 2*ndim
  REAL(KIND=DBL) :: vol
END TYPE

TYPE REGION
  TYPE(HYPERCUBE) :: h
  INTEGER :: splitDim
  INTEGER :: fdim
  TYPE(ESTERR), DIMENSION(:), POINTER :: ee !size is fdim
  REAL(KIND=DBL) :: errmax
END TYPE

TYPE RULE
  INTEGER :: ndim !dimensionality
  INTEGER :: fdim !number of functions
  INTEGER :: num_points  !number of evaluation points
  INTEGER :: num_regions  !max number of regions to be evaluated at once
  REAL(KIND=DBL), DIMENSION(:), POINTER :: pts  !points to evaluate, size is num_regions*num_pts*ndim
  REAL(KIND=DBL), DIMENSION(:), POINTER :: vals !size is num_regions*num_pts*fdim; note this is never allocated
  INTEGER :: evalError  !determines which function to use to evaluate error; 1 for G/M, 0 for QUADPACK
  INTEGER :: destroy  !determines whether we deallocate p, defined below; 1 for yes, 0 for no
  
  !NOTE: First major difference; because Fortran is pass by reference, we use a literal RULE rather than a pointer to RULE
  !In addition, we combine RULE and RULE75GENZMALIK; the below components are exclusive to RULE75GENZMALIK
  REAL(KIND=DBL), DIMENSION(:), POINTER :: widthLambda, widthLambda2, p !temporary arrays of length ndim
  REAL(KIND=DBL) :: weight1, weight3, weight5, weightE1, weightE3 !dimension dependent constants
END TYPE

TYPE HEAP
  INTEGER :: n
  INTEGER :: nalloc
  TYPE(REGION), DIMENSION(:), POINTER :: items !heap items are regions
  INTEGER :: fdim
  TYPE(ESTERR), DIMENSION(:), POINTER :: ee !array of length fdim of the total integrand and error
END TYPE


CONTAINS

REAL(KIND=DBL) FUNCTION errMax(fdim, ee)
  INTEGER, INTENT(IN) :: fdim
  TYPE(ESTERR), DIMENSION(fdim), INTENT(IN) :: ee
  
  INTEGER :: i
  
  errMax = 0
  DO i = 1, fdim
    IF(ee(i)%err_.GT.errMax) errMax = ee(i)%err_
  END DO
END FUNCTION

REAL(KIND=DBL) FUNCTION compute_vol(h) RESULT(vol)
  TYPE(HYPERCUBE), INTENT(IN) :: h
  
  INTEGER :: i
  
  vol = 1
  DO i = 1, h%ndim
    vol = vol * 2._DBL * h%dat(i + h%ndim)
  END DO
END FUNCTION

TYPE(HYPERCUBE) FUNCTION make_hypercube(ndim, center, halfwidth) RESULT(h)
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: center, halfwidth !size ndim
  
  INTEGER :: i
  
  h%ndim = ndim 
  ALLOCATE(h%dat(2*ndim))
  h%vol = 0
  IF(ASSOCIATED(h%dat)) THEN
    DO i = 1, ndim
      h%dat(i) = center(i)
      h%dat(i+ndim) = halfwidth(i)
    END DO
  END IF
  h%vol = compute_vol(h)
END FUNCTION

TYPE(HYPERCUBE) FUNCTION make_hypercube_range(ndim, xmin, xmax) RESULT(h)
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax !size ndim
  
  INTEGER :: i
  
  h = make_hypercube(ndim, xmin, xmax)
  IF(ASSOCIATED(h%dat)) THEN
    DO i=1, ndim
      h%dat(i) = 0.5_DBL * (xmin(i) + xmax(i))
      h%dat(i+ndim) = 0.5_DBL * (xmax(i) - xmin(i))
    END DO
  END IF
  h%vol = compute_vol(h)  
END FUNCTION

SUBROUTINE destroy_hypercube(h)
  TYPE(HYPERCUBE), INTENT(INOUT) :: h
  
  IF(ASSOCIATED(h%dat)) DEALLOCATE(h%dat)
  NULLIFY(h%dat)
  h%ndim = 0
END SUBROUTINE

TYPE(REGION) FUNCTION make_region(h, fdim) RESULT(R)
  TYPE(HYPERCUBE), INTENT(IN) :: h
  INTEGER , INTENT(IN):: fdim
  
  R%h = make_hypercube(h%ndim, h%dat(1:), h%dat(1+h%ndim:))
  R%splitDim = 1
  R%fdim = fdim
  IF(ASSOCIATED(R%h%dat)) THEN
    ALLOCATE(R%ee(fdim))
  ELSE
    NULLIFY(R%ee)
  END IF
  R%errmax = HUGE(0._DBL)
END FUNCTION
  
SUBROUTINE destroy_region(R)
  TYPE(REGION), INTENT(INOUT) :: R
  
  CALL destroy_hypercube(R%h)
  IF(ASSOCIATED(R%ee)) DEALLOCATE(R%ee)
  NULLIFY(R%ee)
  R%errmax = 0._DBL
END SUBROUTINE
  
INTEGER FUNCTION cut_region(R1, R2)
    TYPE(REGION), INTENT(INOUT) :: R1
    TYPE(REGION), INTENT(OUT) :: R2
    
    INTEGER :: d, ndim
    d = R1%splitDim
    ndim = R1%h%ndim
    R2 = R1
    R1%h%dat(d + ndim) = 0.5_DBL * R1%h%dat(d + ndim)
    R1%h%vol = 0.5_DBL * R1%h%vol
    R2%h = make_hypercube(ndim, R1%h%dat(1:), R1%h%dat(1+ndim:))
    IF(.NOT.ASSOCIATED(R2%h%dat)) THEN !note: original code was if(!R2->h.data)
      cut_region = 1
    ELSE
      R1%h%dat(d) = R1%h%dat(d) - R1%h%dat(d+ndim)
      R2%h%dat(d) = R2%h%dat(d) + R2%h%dat(d+ndim)
      ALLOCATE(R2%ee(R2%fdim))
      IF(ASSOCIATED(R2%ee)) THEN !note: original code was return R2->ee == NULL
        cut_region = 0
      ELSE
        cut_region = 1
      END IF
    END IF
END FUNCTION

SUBROUTINE destroy_rule(r)
  TYPE(RULE), INTENT(INOUT) :: r
  
  IF(r%destroy.NE.0) THEN 
    IF(ASSOCIATED(r%p)) DEALLOCATE(r%p)
  END IF
  IF(ASSOCIATED(r%pts)) DEALLOCATE(r%pts)
  NULLIFY(r%pts, r%vals, r%widthLambda, r%widthLambda2, r%p)
END SUBROUTINE

INTEGER FUNCTION alloc_rule_pts(r, num_regions)
  TYPE(RULE), INTENT(INOUT) :: r
  INTEGER, INTENT(IN):: num_regions
  
  INTEGER :: pts_dim, num_regions_tmp
  
  num_regions_tmp = num_regions
  IF(num_regions_tmp.GT.r%num_regions) THEN
    IF(ASSOCIATED(r%pts)) DEALLOCATE(r%pts)
    NULLIFY(r%pts, r%vals)
    r%num_regions = 0
    !allocate extra so we don't waste time reallocating every time this function is called
    num_regions_tmp = 2 * num_regions_tmp
    pts_dim = num_regions_tmp * r%num_points * (r%ndim + r%fdim)
    ALLOCATE(r%pts(pts_dim))
    IF( ((r%ndim + r%fdim).GT.0).AND.(.NOT.ASSOCIATED(r%pts))) THEN !Note, the original code had && !r->pts
      alloc_rule_pts = 1
      RETURN
    END IF
    r%vals => r%pts(1+num_regions_tmp*r%num_points*r%ndim:) !pointer arithmetic isn't allowed in Fortran
    r%num_regions = num_regions_tmp
  END IF
  alloc_rule_pts = 0  
END FUNCTION

TYPE(RULE) FUNCTION make_rule(ndim, fdim, num_points, evalError, destroy) RESULT(r)
  INTEGER, INTENT(IN):: ndim, fdim, num_points, evalError, destroy
  
  NULLIFY(r%pts, r%vals, r%widthLambda, r%widthLambda2, r%p)
  r%num_regions = 0; 
  r%ndim = ndim; r%fdim = fdim; r%num_points = num_points; r%evalError = evalError; r%destroy = destroy 
  r%weight1 = 0._DBL; r%weight3 = 0._DBL; r%weight5 = 0._DBL; r%weightE1 = 0._DBL; r%weightE3 = 0._DBL
  
END FUNCTION
  
!note: all regions must have same fdim
INTEGER FUNCTION eval_regions(nR, Re, fdata, r, f, fv) !fdata is discarded, f is our normal integrand
  OPTIONAL :: f, fv
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: nR
  TYPE(REGION), DIMENSION(:), INTENT(INOUT) :: Re !note, nR better be <= the size of R
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  TYPE(RULE), INTENT(INOUT):: r
  
  INTEGER :: iR
  
  IF(nR.EQ.0) THEN
    eval_regions = 0
    RETURN !nothing to evaluate
  END IF
  IF(.NOT.(PRESENT(f).XOR.PRESENT(fv))) THEN !pass either f or fv
    eval_regions = 1 
    RETURN
  END IF
  IF(r%evalError.EQ.1) THEN 
    IF(PRESENT(fv)) THEN
      IF(rule75genzmalik_evalError(r, Re(1)%fdim, fdata, nR, Re, fv = fv).EQ.1) THEN
        eval_regions = 1
        RETURN
      END IF
    ELSE IF(PRESENT(f)) THEN
      IF(rule75genzmalik_evalError(r, Re(1)%fdim, fdata, nR, Re, f = f).EQ.1) THEN
        eval_regions = 1
        RETURN
      END IF
    END IF
  ELSE IF(r%evalError.EQ.0) THEN
    IF(PRESENT(fv)) THEN
      IF(rule15gauss_evalError(r, Re(1)%fdim, fdata, nR, Re, fv = fv).EQ.1) THEN
        eval_regions = 1
        RETURN
      END IF
    ELSE IF(PRESENT(f)) THEN
      IF(rule15gauss_evalError(r, Re(1)%fdim, fdata, nR, Re, f = f).EQ.1) THEN
        eval_regions = 1
        RETURN
      END IF
    END IF
  END IF
  DO iR = 1,nR
    Re(iR)%errmax = errMax(Re(1)%fdim, Re(iR)%ee)
  END DO
  eval_regions = 0
END FUNCTION

INTEGER FUNCTION ls0(n) !return the least significant 0 bit of n; ex., 000 would return 0, 101 would return 1
  INTEGER, INTENT(IN) :: n
  
  INTEGER :: i, n_tmp
  
  ls0 = 0
  n_tmp = ABS(n) !make sure n_tmp is positive
  !equivalent to counting the number of trailing 1's
  DO WHILE(MOD(n_tmp,2).EQ.1)
    n_tmp = ISHFT(n_tmp,-1)
    ls0 = ls0+1
  END DO
END FUNCTION

!Evaluate the integration points for all 2^n points (+/-r,...,+/-r)
!Gray code ordering is used to minimize the number of coordinate updates
  
SUBROUTINE evalR_Rfs(pts, ndim, p, c, r)
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: pts, p
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: c, r !note, ndim better be <= the size of p,c,r
  INTEGER, INTENT(IN) :: ndim
  
  INTEGER :: i, j, signs, mask, d, pts_loc !keeps track of pointer arithmetic
  
  signs = 0 ! 0/1 bit = +/- for corresponding element of r
  pts_loc = 0
  
  !We start with the point where r is ADDed in every coordinate (this implies signs=0)
  DO i = 1, ndim
    p(i) = c(i) + r(i)
  END DO
  
  !Loop through the points in Gray-code ordering
  i = 0
  DO WHILE(.TRUE.) !CAREFUL!!!! INTENTIONAL INFINITE LOOP
    DO j = 1,ndim
      pts(j+pts_loc) = p(j)
    END DO
    pts_loc = pts_loc + ndim
    
    d = ls0(i)
    IF(d.GE.ndim) THEN
      EXIT
    END IF
    mask = ISHFT(1, d) !Shift left by d
    signs = IEOR(signs, mask)
    d = d + 1
    IF(IAND(mask,signs).NE.0) THEN
      p(d) = c(d) - r(d)
    ELSE
      p(d) = c(d) + r(d)
    END IF
    i = i + 1
  END DO
END SUBROUTINE
 
SUBROUTINE evalRR0_0fs(pts, ndim, p, c, r)
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: pts, p
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: c, r !note, ndim better be <= the size of p, c, r
  INTEGER, INTENT(IN) :: ndim
  
  INTEGER :: i, j, k, pts_loc !pts_loc keeps track of pointer arithmetic
  
  pts_loc = 0
  
  DO i = 1, ndim-1  !loops ndim-1 times
    p(i) = c(i) - r(i)
    DO j = i+1, ndim  !loops ndim-i times
      p(j) = c(j) - r(j)
      DO k = 1,ndim
        pts(k + pts_loc) = p(k)
      END DO
      pts_loc = pts_loc + ndim
      
      p(i) = c(i) + r(i)
      DO k = 1,ndim
        pts(k + pts_loc) = p(k)
      END DO
      pts_loc = pts_loc + ndim
      
      p(j) = c(j) + r(j)
      DO k = 1,ndim
        pts(k + pts_loc) = p(k)
      END DO
      pts_loc = pts_loc + ndim
      
      p(i) = c(i) - r(i)
      DO k = 1,ndim
        pts(k + pts_loc) = p(k)
      END DO 
      pts_loc = pts_loc + ndim
      
      p(j) = c(j)
    END DO
    p(i) = c(i)
  END DO
END SUBROUTINE

SUBROUTINE evalR0_0fs4d(pts, ndim, p, c, r1, r2)
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: pts, p
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: c, r1, r2 !note, ndim better be <= the size of p, c, r1, r2
  INTEGER, INTENT(IN) :: ndim
  
  INTEGER :: i, j, pts_loc !pts_loc keeps track of pointer arithmetic
  
  pts_loc = 0
  
  DO j = 1, ndim
    pts(j+pts_loc) = p(j)
  END DO
  pts_loc = pts_loc + ndim
  
  DO i = 1, ndim
    p(i) = c(i) - r1(i)
    DO j = 1, ndim
      pts(j+pts_loc) = p(j)
    END DO
    pts_loc = pts_loc + ndim
    
    p(i) = c(i) + r1(i)
    DO j = 1, ndim
      pts(j+pts_loc) = p(j)
    END DO
    pts_loc = pts_loc + ndim
    
    p(i) = c(i) - r2(i)
    DO j = 1, ndim
      pts(j+pts_loc) = p(j)
    END DO
    pts_loc = pts_loc + ndim
    
    p(i) = c(i) + r2(i)
    DO j = 1, ndim
      pts(j+pts_loc) = p(j)
    END DO
    pts_loc = pts_loc + ndim
    
    p(i) = c(i)
  END DO 
END SUBROUTINE

!   Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
!   cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
!   and A. A. Malik.  See:
!
!         A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
!         symmetric numerical integration rules," SIAM
!         J. Numer. Anal. 20 (3), 580-588 (1983).

INTEGER FUNCTION rule75genzmalik_evalError(r, fdim, fdata, nR, Re, f, fv) 
  OPTIONAL :: f, fv
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  TYPE(RULE), INTENT(INOUT) :: r
  INTEGER, INTENT(IN) :: fdim, nR
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  TYPE(REGION), DIMENSION(:), INTENT(INOUT) :: Re
  
  REAL(KIND=DBL) :: lambda2, lambda4, lambda5, weight2, weight4, weightE2, weightE4, ratio
  INTEGER :: i, j, iR, ndim, npts
  REAL(KIND=DBL), DIMENSION(:), POINTER :: diff, pts, vals, center, halfwidth, v
  REAL(KIND=DBL) :: results, res5th, val0, sum2, sum3, sum4, sum5, v0, v1, v2, v3, maxdiff, dimDiffMax
  INTEGER :: k, k0, v_loc !needed to keep track of pointer arithmetic
  
  !lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19)
  
  lambda2 = 0.3585685828003180919906451539079374954541_DBL
  lambda4 = 0.9486832980505137995996680633298155601160_DBL
  lambda5 = 0.6882472016116852977216287342936235251269_DBL
  weight2 = 980._DBL / 6561._DBL
  weight4 = 200._DBL / 19683._DBL
  weightE2 = 245._DBL / 486._DBL
  weightE4 = 25._DBL / 729._DBL
  ratio = (lambda2 * lambda2) / (lambda4 * lambda4)
  
  ndim = r%ndim
  npts = 0
  v_loc = 0
  NULLIFY(diff, pts, vals, center, halfwidth, v)
  
  IF(.NOT.(PRESENT(f).XOR.PRESENT(fv))) THEN !pass either f or fv
    rule75genzmalik_evalError = 1
    RETURN
  END IF
  
  IF(alloc_rule_pts(r, nR).EQ.1) THEN
    rule75genzmalik_evalError = 1
    RETURN
  END IF
  
  pts => r%pts
  vals => r%vals
  
  DO iR = 1, nR
    center => Re(iR)%h%dat
    halfwidth => Re(iR)%h%dat(1+ndim:) !can't do pointer arithmetic in Fortran
    
    DO i = 1, ndim
      r%p(i) = center(i)
      r%widthLambda2(i) = halfwidth(i) * lambda2
      r%widthLambda(i) = halfwidth(i) * lambda4
    END DO
    
    !Evaluate points in the center, in (lambda2,0,...,0) and (lambda3=lambda4, 0,...,0). 
    CALL evalR0_0fs4d(pts(1+npts*ndim:), ndim, r%p, center, r%widthLambda2, r%widthLambda) !can't do pointer arithmetic in Fortran
    npts = npts + 1 + 2*2*ndim
    
    !Calculate points for (lambda4, lambda4, 0, ...,0)
    CALL evalRR0_0fs(pts(1+npts*ndim:), ndim, r%p, center, r%widthLambda) !can't do pointer arithmetic in Fortran
    npts = npts + 2 * ndim * (ndim-1)
    
    !Calculate points for (lambda5, lambda5, ..., lambda5)
    DO i = 1, ndim
      r%widthLambda(i) = halfwidth(i) * lambda5
    END DO
    CALL evalR_Rfs(pts(1+npts*ndim:), ndim, r%p, center, r%widthLambda)
    npts = npts + ISHFT(1, ndim)
  END DO
  
  !Evaluate the integrand function(s) at all the points
  IF(PRESENT(fv)) THEN
    IF(.NOT.(fv(ndim, npts, pts, fdata, fdim, vals).EQ.0)) THEN
      rule75genzmalik_evalError = 1
      RETURN
    END IF
  ELSE IF(PRESENT(f)) THEN
    IF(.NOT.(fv_wrapper(f, ndim, npts, pts, fdata, fdim, vals).EQ.0)) THEN
      rule75genzmalik_evalError = 1
      RETURN
    END IF
  END IF
  
  
  !we are done with the points, and so we can re-use the pts array to store the maximum difference diff[i] in each dimension for each hypercube
  
  diff => pts
  DO i = 1, ndim*nR
    diff(i) = 0
  END DO
  
  DO j = 1, fdim
    v => vals(j:)
    v_loc = 1
    
    DO iR = 1, nR
      sum2 = 0._DBL; sum3 = 0._DBL; sum4 = 0._DBL; sum5 = 0._DBL
      k0 = 0
      !accumulate j-th function values into j-th integrals
		  !NOTE: this relies on the ordering of the eval functions
		  !above, as well as on the internal structure of
		  !the evalR0_0fs4d function
      
      val0 = v(v_loc + 0*fdim) !central point
      k0 = k0 + 1
      
      DO k = 0, ndim-1
        v0 = v(v_loc + (k0+4*k)*fdim)
        v1 = v(v_loc + (k0+4*k+1)*fdim)
        v2 = v(v_loc + (k0+4*k+2)*fdim)
        v3 = v(v_loc + (k0+4*k+3)*fdim)
        
        sum2 = sum2 + v0 + v1
        sum3 = sum3 + v2 + v3
        
        diff(1+iR*ndim + k) = diff(1 + iR*ndim -ndim + k) + ABS(v0 + v1 - 2*val0 - ratio * (v2 + v3 - 2*val0))
      END DO
      
      k0 = k0 + 4*k !note here k = ndim 
      
      DO k = 0, 2*ndim*(ndim-1)-1
        sum4 = sum4 + v(v_loc + (k0+k)*fdim)
      END DO
      
      k0 = k0 + k !note here k = 2*ndim*(ndim-1)
      
      DO k = 0, ISHFT(1, ndim)-1
        sum5 = sum5 + v(v_loc + (k0+k)*fdim)
      END DO
      
      !Calculate fifth and seventh order results
      results = Re(iR)%h%vol * (r%weight1 * val0 + weight2 * sum2 + r%weight3 * sum3 + weight4 * sum4 + r%weight5 * sum5)
      res5th = Re(iR)%h%vol * (r%weightE1 * val0 + weightE2 * sum2 + r%weightE3 * sum3 + weightE4 * sum4)
      
      Re(iR)%ee(j)%val_ = results
      Re(iR)%ee(j)%err_ = ABS(res5th - results)
      
      v_loc = v_loc + r%num_points * fdim
      
    END DO
  END DO 
  
  ! figure out dimension to split
  DO iR = 1, nR
    maxdiff = 0
    dimDiffMax = 1
    DO i = 1, ndim
      IF(diff(iR*ndim + i - ndim).GT.maxdiff) THEN
        maxdiff = diff(iR*ndim + i - ndim)
        dimDiffMax = i
      END IF
    END DO
    Re(iR)%splitDim = dimDiffMax
  END DO
  
  NULLIFY(diff, pts, vals, center, halfwidth, v)
  rule75genzmalik_evalError = 0
  
  
END FUNCTION
  
  
TYPE(RULE) FUNCTION make_rule75genzmalik(ndim, fdim) RESULT(r)
  INTEGER, INTENT(IN) :: ndim, fdim
  
  r = make_rule(ndim, fdim, 1 + 2 * 2*ndim + 2*ndim*(ndim-1) + ISHFT(1,ndim), 1, 1)
  r%weight1 = ((12824._DBL - 9120._DBL * REAL(ndim, KIND=DBL) + 400._DBL * REAL(ndim, KIND=DBL)**2)/19683._DBL)
  r%weight3 = ((1820._DBL - 400._DBL * REAL(ndim, KIND=DBL)) / 19683._DBL)
  r%weight5 = 6859._DBL / 19683._DBL / REAL(ISHFT(1,ndim), KIND=DBL)
  r%weightE1 = ((729._DBL - 950._DBL * REAL(ndim, KIND=DBL) + 50._DBL * REAL(ndim, KIND=DBL)**2) / 729._DBL)
  r%weightE3 = ((265._DBL - 100._DBL * REAL(ndim, KIND=DBL))/1458._DBL)
  
  ALLOCATE(r%p(ndim*3))
  r%widthLambda => r%p(1+ndim:)
  r%widthLambda2 => r%p(1 + 2*ndim:)  
  
END FUNCTION

!1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in GNU GSL (which in turn is based on QUADPACK).

INTEGER FUNCTION rule15gauss_evalError(r, fdim, fdata, nR, Re, f, fv) 
  OPTIONAL :: f, fv
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  TYPE(RULE), INTENT(INOUT) :: r
  INTEGER, INTENT(IN) :: fdim, nR
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  TYPE(REGION), DIMENSION(:), INTENT(INOUT) :: Re
  
  INTEGER :: n, j, j2, k, iR, npts, vk_loc
  REAL(KIND=DBL), DIMENSION(8) :: xgk, wgk
  REAL(KIND=DBL), DIMENSION(4) :: wg
  REAL(KIND=DBL), DIMENSION(:), POINTER :: pts, vals, vk
  REAL(KIND=DBL) :: w, v, center, halfwidth, result_gauss, result_kronrod, result_abs, result_asc, mean, err_, scale_, min_err
  
  ! Gauss quadrature weights and kronrod quadrature abscissae and
	! weights as evaluated with 80 decimal digit arithmetic by
	! L. W. Fullerton, Bell Labs, Nov. 1981.
  
  xgk = (/0.991455371120812639206854697526329_DBL, 0.949107912342758524526189684047851_DBL, &
    & 0.864864423359769072789712788640926_DBL, 0.741531185599394439863864773280788_DBL, & 
    & 0.586087235467691130294144838258730_DBL, 0.405845151377397166906606412076961_DBL, &
    & 0.207784955007898467600689403773245_DBL, 0._DBL /)

  wg = (/0.129484966168869693270611432679082_DBL, 0.279705391489276667901467771423780_DBL, &
    & 0.381830050505118944950369775488975_DBL, 0.417959183673469387755102040816327_DBL /)

  wgk = (/0.022935322010529224963732008058970_DBL, 0.063092092629978553290700663189204_DBL, & 
    & 0.104790010322250183839876322541518_DBL, 0.140653259715525918745189590510238_DBL, &
    & 0.169004726639267902826583426598550_DBL, 0.190350578064785409913256402421014_DBL, &
    & 0.204432940075298892414161999234649_DBL, 0.209482141084727828012999174891714_DBL /)

  npts = 0
  n = 8
  vk_loc = 0
  NULLIFY(pts, vals, vk)
  
  IF(.NOT.(PRESENT(f).XOR.PRESENT(fv))) THEN !pass either f or fv
    rule15gauss_evalError = 1
    RETURN
  END IF
  
  IF(alloc_rule_pts(r, nR).NE.0) THEN
    rule15gauss_evalError = 1
    RETURN
  END IF
  
  pts => r%pts; vals => r%vals
  
  !Note, unnecessary to call xkg(8) 
  
  DO iR = 1, nR
    center = Re(iR)%h%dat(1)
    halfwidth = Re(iR)%h%dat(2)
    npts = npts + 1
    pts(npts) = center
    
    DO j = 1, (n-1)/2
      j2 = 2*j
      w = halfwidth * xgk(j2)
      npts = npts + 1
      pts(npts) = center - w
      npts = npts + 1
      pts(npts) = center + w
    END DO
    
    DO j = 1, n/2
      j2 = 2*j-1
      w = halfwidth * xgk(j2)
      npts = npts + 1
      pts(npts) = center - w
      npts = npts + 1
      pts(npts) = center + w      
    END DO
    
    Re(iR)%splitDim = 1 !no choice, we're in 1D
  END DO
  
  IF(PRESENT(fv)) THEN
    IF(fv(1, npts, pts, fdata, fdim, vals).NE.0) THEN
      rule15gauss_evalError = 1
      RETURN
    END IF
  ELSE IF(PRESENT(f)) THEN
    IF(fv_wrapper(f, 1, npts, pts, fdata, fdim, vals).NE.0) THEN
      rule15gauss_evalError = 1
      RETURN
    END IF
  END IF
    
  DO k = 1, fdim
    vk => vals(k:)
    vk_loc = 1
    DO iR = 1, nR
      halfwidth = Re(iR)%h%dat(2)
      result_gauss = vk(vk_loc) * wg(n/2)
      result_kronrod = vk(vk_loc) * wgk(n)
      result_abs = ABS(result_kronrod)
      npts = 1
      !accumulate integrals
      DO j = 1, (n-1)/2
        j2 = 2*j
        v = vk(vk_loc+fdim*npts) + vk(vk_loc+fdim*npts+fdim)
        result_gauss = result_gauss + wg(j) * v
        result_kronrod = result_kronrod + wgk(j2) * v
        result_abs = result_abs + wgk(j2) * ABS(vk(vk_loc+fdim*npts)) + ABS(vk(vk_loc+fdim*npts+fdim))
        npts = npts + 2
      END DO
      DO j = 1, n/2
        j2 = 2*j-1
        result_kronrod = result_kronrod + wgk(j2) * (vk(vk_loc+fdim*npts) + vk(vk_loc+fdim*npts+fdim))
        result_abs = result_abs + wgk(j2) * ABS(vk(vk_loc+fdim*npts)) + ABS(vk(vk_loc+fdim*npts+fdim))
        npts = npts + 2
      END DO
        
      !integration result
      Re(iR)%ee(k)%val_ = result_kronrod * halfwidth
       
      !error estimate from GSL
        
      mean = result_kronrod * 0.5_DBL
      result_asc = wgk(n) * ABS(vk(vk_loc) - mean)
      npts = 1
        
      DO j = 1, (n-1)/2
        j2 = 2*j
        result_asc = result_asc + wgk(j2) * ABS(vk(vk_loc+fdim*npts)-mean) + ABS(vk(vk_loc+fdim*npts+fdim)-mean)
        npts = npts + 2
      END DO
        
      DO j = 1, n/2
        j2 = 2*j-1
        result_asc = result_asc + wgk(j2) * ABS(vk(vk_loc+fdim*npts)-mean) + ABS(vk(vk_loc+fdim*npts+fdim)-mean)
        npts = npts + 2
      END DO
      err_ = ABS(result_kronrod - result_gauss) * halfwidth
      result_abs = result_abs * halfwidth
      result_asc = result_asc * halfwidth
      IF((result_asc.NE.0._DBL).AND.(err_.NE.0._DBL)) THEN
        scale_ = (200._DBL * err_/result_asc)**(1.5)
        IF(scale_.LT.1._DBL) THEN
          err_ = result_asc*scale_
        ELSE
          err_ = result_asc
        END IF
      END IF
      IF(result_abs.GT.(TINY(0._DBL) / (50._DBL*EPSILON(0._DBL)))) THEN
        min_err = 50._DBL * result_abs * EPSILON(0._DBL) 
        IF(min_err.GT.err_) err_ = min_err  
      END IF
      Re(iR)%ee(k)%err_ = err_
        
      vk_loc = vk_loc + 15 * fdim !can't do pointer arithmetic in Fortran
        
    END DO
  END DO
  
  NULLIFY(pts, vals, vk)
  rule15gauss_evalError = 0    
END FUNCTION 

TYPE(RULE) FUNCTION make_rule15gauss(ndim, fdim) RESULT(r)
  INTEGER, INTENT(IN) :: ndim, fdim
  
  r = make_rule(ndim, fdim, 15, 0, 0)
END FUNCTION




!binary heap implementation (ala _Introduction to Algorithms_ by Cormen, Leiserson, and Rivest), for use as a priority queue of regions to integrate.

SUBROUTINE heap_resize(h, nalloc)
  TYPE(HEAP), INTENT(INOUT) :: h
  INTEGER, INTENT(IN) :: nalloc
  
  TYPE(REGION), DIMENSION(:), POINTER :: items_tmp
  INTEGER :: i, nalloc_old
  nalloc_old = h%nalloc
  h%nalloc = nalloc
  NULLIFY(items_tmp)
  
  IF(nalloc.NE.0) THEN
    ALLOCATE(items_tmp(nalloc))
    DO i = 1, MIN(nalloc_old, nalloc) !makes sure that if we shrink size we don't copy everything
      items_tmp(i) = h%items(i)
    END DO
    IF(ASSOCIATED(h%items)) DEALLOCATE(h%items)  !CAREFUL, this does not destroy the regions contained in items
    h%items => items_tmp
    NULLIFY(items_tmp)
  ELSE
    !CAREFUL, this does not destroy the regions contained in items
    IF(ASSOCIATED(h%items)) DEALLOCATE(h%items)
    NULLIFY(h%items)    
  END IF
END SUBROUTINE

TYPE(HEAP) FUNCTION heap_alloc(nalloc, fdim) RESULT(h)
  INTEGER, INTENT(IN) :: nalloc, fdim
  
  INTEGER :: i
  h%n = 0
  h%nalloc = 0
  NULLIFY(h%items, h%ee)
  h%fdim = fdim
  ALLOCATE(h%ee(fdim))
  IF(ASSOCIATED(h%ee)) THEN
    DO i = 1, fdim
      h%ee(i)%val_ = 0._DBL
      h%ee(i)%err_ = 0._DBL
    END DO
    CALL heap_resize(h, nalloc)
  END IF
END FUNCTION

!note that heap_free does not deallocate anything referenced by the items 
SUBROUTINE heap_free(h)
  TYPE(HEAP), INTENT(INOUT) :: h
  
  h%n = 0
  CALL heap_resize(h, 0)
  h%fdim = 0
  IF(ASSOCIATED(h%ee)) DEALLOCATE(h%ee)
  NULLIFY(h%ee)
END SUBROUTINE

INTEGER FUNCTION heap_push(h, hi)
  TYPE(HEAP), INTENT(INOUT) :: h
  TYPE(REGION), INTENT(IN) :: hi
  
  INTEGER :: i, fdim, insert, parent
  fdim = h%fdim
  
  DO i = 1, fdim
    h%ee(i)%val_ = h%ee(i)%val_ + hi%ee(i)%val_
    h%ee(i)%err_ = h%ee(i)%err_ + hi%ee(i)%err_
  END DO
  insert = h%n
  h%n = h%n + 1
  IF(h%n.GT.h%nalloc) THEN
    CALL heap_resize(h, h%n*2)
    IF(.NOT.ASSOCIATED(h%items)) THEN !note: original code was if(!h->items)
      heap_push = 1
      RETURN
    END IF
  END IF
  DO WHILE(insert.NE.0)
    parent = (insert-1)/2
    IF(hi%errmax.LE.h%items(1+parent)%errmax) THEN
      EXIT
    END IF
    h%items(1+insert) = h%items(1+parent)
    insert = parent
  END DO
  h%items(1+insert) = hi
  heap_push = 0
END FUNCTION

INTEGER FUNCTION heap_push_many(h, ni, hi)
  TYPE(HEAP), INTENT(INOUT) :: h
  INTEGER, INTENT(IN) :: ni
  TYPE(REGION), DIMENSION(:), INTENT(IN) :: hi
  
  INTEGER :: i
  DO i = 1, ni
    IF(heap_push(h, hi(i)).EQ.1) THEN
      heap_push_many = 1
      RETURN
    END IF
  END DO
  heap_push_many = 0
END FUNCTION

TYPE(REGION) FUNCTION heap_pop(h) RESULT(ret)
  TYPE(HEAP), INTENT(INOUT) :: h
  
  INTEGER :: i, j, n, child, largest, fdim
  TYPE(REGION) :: swap
  
  IF(h%n.EQ.0) THEN
    WRITE(*,*) 'attempted to pop empty heap, exiting'
    STOP
  END IF
  
  ret = h%items(1)
  i = 0
  h%n = h%n - 1
  n = h%n
  h%items(1+i) = h%items(1+n)
  child = 2*i + 1
  DO WHILE(child.LT.n)
    IF(h%items(1+child)%errmax.LE.h%items(1+i)%errmax) THEN
      largest = i
    ELSE
      largest = child
    END IF
    
    child = child+1
    IF((child.LT.n).AND.(h%items(1+largest)%errmax.LT.h%items(1+child)%errmax)) THEN
      largest = child
    END IF
    IF(largest.EQ.i) THEN
      EXIT
    END IF
    swap = h%items(1+i)
    h%items(1+i) = h%items(1+largest)
    i = largest
    h%items(1+i) = swap
    child = 2*i+1
  END DO
  
  
  fdim = h%fdim
  DO j = 1, fdim
    h%ee(j)%val_ = h%ee(j)%val_ - ret%ee(j)%val_
    h%ee(j)%err_ = h%ee(j)%err_ - ret%ee(j)%err_
  END DO
  
END FUNCTION

INTEGER FUNCTION converged(fdim, ee, reqAbsError, reqRelError, norm)
  INTEGER, INTENT(IN) :: fdim, norm
  REAL(KIND=DBL), INTENT(IN) ::  reqAbsError, reqRelError
  TYPE(ESTERR), DIMENSION(:), INTENT(IN) :: ee
  
  INTEGER :: j
  REAL(KIND=DBL) :: maxerr, serr, err_, maxval_, sval, val_, absval
  !Body of convergence test
  
  SELECT CASE(norm)
    CASE(1) !Error individual
      DO j = 1, fdim
        IF((ee(j)%err_.GT.reqAbsError).AND.(ee(j)%err_.GT.(ABS(ee(j)%val_)*reqRelError)))THEN
          converged = 0
          RETURN
        END IF
      END DO
      converged = 1
      RETURN
      
    CASE(2) !Error paired
      DO j = 1, fdim -1, 2
        !Scale to avoid overflow/underflow
        IF(ee(j)%err_.GT.ee(j+1)%err_) THEN
          maxerr = ee(j)%err_
        ELSE
          maxerr = ee(j+1)%err_
        END IF
        IF(ee(j)%val_.GT.ee(j+1)%val_) THEN
          maxval_ = ee(j)%val_
        ELSE
          maxval_ = ee(j+1)%val_
        END IF
        
        IF(maxerr.GT.0._DBL) THEN
          serr = 1._DBL/maxerr
        ELSE
          serr = 1._DBL
        END IF
        IF(maxval_.GT.0._DBL) THEN
          sval = 1._DBL/maxval_
        ELSE
          sval = 1
        END IF
        
        err_ = SQRT( (ee(j)%err_ * serr)**2 + (ee(j+1)%err_ * serr)**2 ) * maxerr
        val_ = SQRT( (ee(j)%val_ * sval)**2 + (ee(j+1)%val_ * sval)**2 ) * maxval_
        
        IF( (err_.GT.reqAbsError).AND.(err_.GT.(val_*reqRelError))) THEN
          converged = 0
          RETURN
        END IF
        
      END DO
      !check for case that fdim is odd
      IF(j.LT.fdim) THEN
        IF((ee(j)%err_.GT.reqAbsError).AND.(ee(j)%err_.GT.(ABS(ee(j)%val_)*reqRelError))) THEN
          converged = 0
          RETURN
        END IF
      END IF
      converged = 1
      RETURN
    
    CASE(3) !Error L1
      err_ = 0.; val_ = 0.
      DO j = 1, fdim
        err_ = err_ + ee(j)%err_
        val_ = val_ + ABS(ee(j)%val_)
      END DO
      IF((err_.LE.reqAbsError).OR.(err_.LE.(val_*reqAbsError))) THEN
        converged = 1
        RETURN
      ELSE
        converged = 0
        RETURN
      END IF
      
    CASE(4) !Error LINF
      err_ = 0.; val_ = 0.
      DO j = 1, fdim
        absval = ABS(ee(j)%val_)
        IF(ee(j)%err_.GT.err_) err_ = ee(j)%err_
        IF(absval.GT.val_) val_ = absval
      END DO
      IF((err_.LE.reqAbsError).OR.(err_.LE.(val_*reqAbsError))) THEN
        converged = 1
        RETURN
      ELSE
        converged = 0
        RETURN
      END IF
      
    CASE(5) !Error L2
      maxerr = 0.; maxval_ = 0.; err_ = 0.; val_ = 0.
      !Scale to avoid overflow/underflow
      DO j = 1, fdim
        IF(ee(j)%err_.GT.maxerr) maxerr = ee(j)%err_
        absval = ABS(ee(j)%val_)
        IF(absval.GT.maxval_) maxval_ = absval
      END DO
      
      IF(maxerr.GT.0._DBL) THEN
        serr = 1._DBL/maxerr
      ELSE
        serr = 1._DBL
      END IF
      IF(maxval_.GT.0._DBL) THEN
        sval = 1._DBL/maxval_
      ELSE
        sval = 1
      END IF
      
      DO j = 1, fdim
        err_ = err_ + (ee(j)%err_ * serr)**2
        val_ = val_ + (ee(j)%val_ * sval)**2
      END DO
      
      err_ = SQRT(err_) * maxerr
      val_ = SQRT(val_) * maxval_
        
      IF((err_.LE.reqAbsError).OR.(err_.LE.(val_*reqAbsError))) THEN
        converged = 1
        RETURN
      ELSE
        converged = 0
        RETURN
      END IF
  END SELECT
  converged = 1
  RETURN !unreachable   
END FUNCTION

!adaptive integration

INTEGER FUNCTION rulecubature(r, fdim, fdata, h, maxEval, reqAbsError, reqRelError, norm, val_, err_, parallel, f, fv)
  OPTIONAL :: f, fv
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  TYPE(RULE), INTENT(INOUT) :: r
  INTEGER, INTENT(IN) :: fdim, maxEval, norm, parallel
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  TYPE(HYPERCUBE), INTENT(IN) :: h
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  
  INTEGER :: numEval, i, j, k, nR_alloc, nR_alloc_old, norm_new, nR, completed_once
  TYPE(HEAP) :: regions
  TYPE(REGION), DIMENSION(:), POINTER :: Re, Re_tmp
  TYPE(ESTERR), DIMENSION(:), POINTER :: ee
  
  NULLIFY(Re, Re_tmp, ee)
  
  IF(.NOT.(PRESENT(f).XOR.PRESENT(fv))) THEN !should pass either f or fv
    rulecubature = 1
    RETURN
  END IF
  
  numEval = 0; nR_alloc = 0; completed_once = 0
  NULLIFY(Re, ee)
  
  norm_new = norm
  IF(fdim.LE.1) norm_new = 1 !Norm is irrelevant
  IF((norm_new.LT.0).OR.(norm_new.GT.5)) THEN !invalid norm
    rulecubature = 1
    RETURN
  END IF
  
  regions = heap_alloc(1, fdim)
  IF((.NOT.ASSOCIATED(regions%ee)).OR.(.NOT.ASSOCIATED(regions%items))) THEN !failed
    IF(ASSOCIATED(ee)) DEALLOCATE(ee)
    CALL heap_free(regions)
    IF(ASSOCIATED(Re)) DEALLOCATE(Re)
    rulecubature = 1
    RETURN
  END IF  
  
  ALLOCATE(ee(fdim))
  IF(.NOT.ASSOCIATED(ee)) THEN !failed
    IF(ASSOCIATED(ee)) DEALLOCATE(ee)
    CALL heap_free(regions)
    IF(ASSOCIATED(Re)) DEALLOCATE(Re)
    rulecubature = 1
    RETURN
  END IF
  
  nR_alloc = 2
  ALLOCATE(Re(nR_alloc))
  nR_alloc_old = nR_alloc
  IF(.NOT.ASSOCIATED(Re)) THEN !failed
    IF(ASSOCIATED(ee)) DEALLOCATE(ee)
    CALL heap_free(regions)
    IF(ASSOCIATED(Re)) DEALLOCATE(Re)
    rulecubature = 1
    RETURN
  END IF
  
  Re(1) = make_region(h, fdim)
  
  !Fortran might not have short circuiting, so split it up
  IF(PRESENT(fv)) THEN
    IF(.NOT.ASSOCIATED(Re(1)%ee)) THEN !failed
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
      RETURN
    ELSE IF(eval_regions(1, Re, fdata, r, fv = fv).NE.0) THEN
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
    ELSE IF(heap_push(regions, Re(1)).NE.0) THEN
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
    END IF
  ELSE IF(PRESENT(f)) THEN
    IF(.NOT.ASSOCIATED(Re(1)%ee)) THEN !failed
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
      RETURN
    ELSE IF(eval_regions(1, Re, fdata, r, f = f).NE.0) THEN
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
    ELSE IF(heap_push(regions, Re(1)).NE.0) THEN
      IF(ASSOCIATED(ee)) DEALLOCATE(ee)
      CALL heap_free(regions)
      IF(ASSOCIATED(Re)) DEALLOCATE(Re)
      rulecubature = 1
    END IF
  END IF
  
  numEval = numEval + r%num_points
  
  DO WHILE((numEval.LT.maxEval).OR.(maxEval.EQ.0))
    IF(converged(fdim, regions%ee, reqAbsError, reqRelError, norm_new).NE.0) EXIT !once converged, exit the while loop
    
    IF(parallel.NE.0) THEN
      ! maximize potential parallelism
	    ! adapted from I. Gladwell, "Vectorization of one
		  ! dimensional quadrature codes," pp. 230--238 in
		  ! _Numerical Integration. Recent Developments,
		  ! Software and Applications_, G. Fairweather and
		  ! P. M. Keast, eds., NATO ASI Series C203, Dordrecht
		  ! (1987), as described in J. M. Bull and
		  ! T. L. Freeman, "Parallel Globally Adaptive
		  ! Algorithms for Multi-dimensional Integration,"
		  ! http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638
		  ! (1994). 

		  ! Basically, this evaluates in one shot all regions
		  ! that *must* be evaluated in order to reduce the
		  ! error to the requested bound: the minimum set of
		  ! largest-error regions whose errors push the total
		  ! error over the bound.

		  ! [Note: Bull and Freeman claim that the Gladwell
		  ! approach is intrinsically inefficent because it
		  ! "requires sorting", and propose an alternative
		  ! algorithm that "only" requires three passes over the
		  ! entire set of regions.  Apparently, they didn't
		  ! realize that one could use a heap data structure, in
		  ! which case the time to pop K biggest-error regions
		  ! out of N is only O(K log N), much better than the
		  ! O(N) cost of the Bull and Freeman algorithm if K <<
		  ! N, and it is also much simpler.]
    
      nR = 0
      DO j = 1,fdim
        ee(j) = regions%ee(j)
      END DO
      
      ! in C this is a do while loop, which is executed at least once
      ! emulate this with an extra term in the while loop
      DO WHILE(((regions%n.GT.0).AND.((numEval.LT.maxEval).OR.(maxEval.EQ.0))).OR.(completed_once.EQ.0))
        completed_once = 1
        IF((nR + 2).GT.nR_alloc) THEN
          nR_alloc = (nR+2)*2
          
          ALLOCATE(Re_tmp(nR_alloc))
          DO k = 1, MIN(nR_alloc_old, nR_alloc) !makes sure that if we shrink the size we don't copy everything
            Re_tmp(k) = Re(k)
          END DO
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)  !CAREFUL, this does not destroy the hypercubes contained in regions
          Re => Re_tmp
          NULLIFY(Re_tmp)
          nR_alloc_old = nR_alloc
          IF(.NOT.ASSOCIATED(RE)) THEN !failed
            IF(ASSOCIATED(ee)) DEALLOCATE(ee)
            CALL heap_free(regions)
            IF(ASSOCIATED(Re)) DEALLOCATE(Re)
            rulecubature = 1
            RETURN
          END IF
        END IF
        Re(nR+1) = heap_pop(regions)
        DO j = 1,fdim
          ee(j)%err_ = ee(j)%err_ - Re(nR+1)%ee(j)%err_
        END DO
        IF(cut_region(Re(nR+1), Re(nR+2)).NE.0) THEN !failed 
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
          RETURN
        END IF
        
        numEval = numEval + r%num_points*2
        nR = nR + 2
        
        IF(converged(fdim, ee, reqAbsError, reqRelError, norm_new).NE.0) EXIT !other regions have smaller errors
      END DO
      
      !Fortran might not have short circuiting, so IF statements are broken up
      IF(PRESENT(fv)) THEN   
        IF(eval_regions(nR, Re, fdata, r, fv = fv).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
          RETURN
        ELSE IF(heap_push_many(regions, nR, Re).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        END IF
      ELSE IF(PRESENT(f)) THEN
        IF(eval_regions(nR, Re, fdata, r, f = f).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
          RETURN
        ELSE IF(heap_push_many(regions, nR, Re).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        END IF
      END IF
      
    ELSE !minimize number of function evaluations
      Re(1) = heap_pop(regions) !get worst region
      
      !Fortran might not have short circuiting, so IF statements are broken up
      IF(PRESENT(fv)) THEN
        IF(cut_region(Re(1), Re(2)).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
          RETURN
        ELSE IF(eval_regions(2, Re, fdata, r, fv = fv).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        ELSE IF(heap_push_many(regions, 2, Re).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        END IF
      ELSE IF(PRESENT(f)) THEN
        IF(cut_region(Re(1), Re(2)).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
          RETURN
        ELSE IF(eval_regions(2, Re, fdata, r, f = f).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        ELSE IF(heap_push_many(regions, 2, Re).NE.0) THEN !failed
          IF(ASSOCIATED(ee)) DEALLOCATE(ee)
          CALL heap_free(regions)
          IF(ASSOCIATED(Re)) DEALLOCATE(Re)
          rulecubature = 1
        END IF
      END IF
      
      numEval = numEval + r%num_points * 2
    END IF
  END DO
  
  !re-sum integral and errors
  
  DO j = 1, fdim
    val_(j) = 0._DBL
    err_(j) = 0._DBL
  END DO
  
  DO i = 1, regions%n
    DO j = 1, fdim
      val_(j) = val_(j) + regions%items(i)%ee(j)%val_
      err_(j) = err_(j) + regions%items(i)%ee(j)%err_
    END DO
    CALL destroy_region(regions%items(i))
  END DO
  
  !success
  IF(ASSOCIATED(ee)) DEALLOCATE(ee)
  CALL heap_free(regions)
  IF(ASSOCIATED(Re)) DEALLOCATE(Re)
  rulecubature = 0
  RETURN
END FUNCTION
    
INTEGER FUNCTION cubature(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, parallel, f, fv)
  OPTIONAL :: f, fv
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, maxEval, norm, parallel
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  
  TYPE(RULE) :: r
  TYPE(HYPERCUBE) :: h
  INTEGER :: status_, i
  
  IF(fdim.EQ.0) THEN !nothing to do
    cubature = 0
    RETURN
  END IF
  
  IF(.NOT.(PRESENT(f).XOR.PRESENT(fv))) THEN !should either f or fv
    cubature = 1
    RETURN
  END IF
  
  IF(ndim.EQ.0) THEN !trivial
    IF(PRESENT(fv)) THEN
      IF(fv(0,1,xmin,fdata,fdim,val_).NE.0) THEN
        cubature = 1
        RETURN
      END IF
    ELSE IF(PRESENT(f)) THEN
      IF(fv_wrapper(f,0,1,xmin,fdata,fdim,val_).NE.0) THEN
        cubature = 1
        RETURN
      END IF
    END IF
    DO i = 1,fdim
      err_(i) = 0
    END DO
    cubature = 0
    RETURN
  END IF
  
  IF(ndim.EQ.1) THEN
    !QUADPACK
    r = make_rule15gauss(ndim, fdim)
  ELSE
    !Genz and Malik
    r = make_rule75genzmalik(ndim, fdim)
  END IF
  DO i = 1, fdim
    val_(i) = 0._DBL
    err_(i) = HUGE(0._DBL) 
  END DO
  
  h = make_hypercube_range(ndim, xmin, xmax)
  IF(.NOT.ASSOCIATED(h%dat)) THEN
    status_ = 1
  ELSE
    IF(PRESENT(fv)) THEN
      status_ = rulecubature(r, fdim, fdata, h, maxEval, reqAbsError, reqRelError, norm, val_, err_, parallel, fv = fv)
    ELSE IF(PRESENT(f)) THEN
      status_ = rulecubature(r, fdim, fdata, h, maxEval, reqAbsError, reqRelError, norm, val_, err_, parallel, f = f)
    END IF
  END IF
  
  CALL destroy_hypercube(h)
  CALL destroy_rule(r)
  cubature = status_
  RETURN
END FUNCTION

INTEGER FUNCTION hcubature_v(fdim, fv, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, fdata) !
  INTERFACE
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, maxEval, norm
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  REAL(KIND=DBL), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(0) :: fdata_new
  
  IF(PRESENT(fdata)) THEN
    hcubature_v = cubature(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, 1, fv = fv)
  ELSE
    hcubature_v = cubature(fdim, fdata_new, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, 1, fv = fv)
  END IF
END FUNCTION

INTEGER FUNCTION hcubature(fdim, f, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, fdata) 
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, maxEval, norm
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  REAL(KIND=DBL), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(0) :: fdata_new
  
  INTEGER :: ret
  
  IF(fdim.EQ.0) THEN !nothing to do
    hcubature = 0
    RETURN
  END IF
  
  IF(PRESENT(fdata)) THEN
    ret = cubature(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, 0, f=f)
  ELSE
    ret = cubature(fdim, fdata_new, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, 0, f=f)
  END IF
  hcubature = ret

END FUNCTION

INTEGER FUNCTION fv_wrapper(f, ndim, npt, x, fdata, fdim, fval)
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: ndim, fdim, npt
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
  
  INTEGER :: i
  
  DO i = 0, npt-1
    IF(f(ndim, x(1+i*ndim:), fdata,  fdim, fval(1+i*fdim:)).NE.0) THEN
      fv_wrapper = 1
      RETURN
    END IF
  END DO
  fv_wrapper = 0
END FUNCTION
  
END MODULE HCUB





