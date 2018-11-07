! Adaptive multidimensional integration of a vector of integrands.
 
  ! Copyright (c) 2005-2013 Steven G. Johnson
 
  ! This program is free software; you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation; either version 2 of the License, or
  ! (at your option) any later version.
 
  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
 
  ! You should have received a copy of the GNU General Public License
  ! along with this program; if not, write to the Free Software
  ! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
 

   ! p-adaptive cubature (adaptive by increasing the degree of the
   ! cubature rule rather than subdividing the domain), using products
   ! of Clenshaw-Curtis rules.  This algorithm may be superior to
   ! Genz-Malik for smooth integrands lacking strongly-localized
   ! features, in moderate dimensions.
   

MODULE PCUB
  USE KIND
  USE CLENCURT
  IMPLICIT NONE
  
INTEGER, PARAMETER :: MAXDIM = 20
INTEGER, PARAMETER :: DEFAULT_MAX_NBUF = ISHFT(1, 20)


TYPE CACHEVAL
  INTEGER, DIMENSION(MAXDIM) :: m
  INTEGER :: mi
  REAL(KIND=DBL), DIMENSION(:), POINTER :: val_
END TYPE

TYPE VALCACHE
  INTEGER :: ncache
  TYPE(CACHEVAL), DIMENSION(:), POINTER :: c
END TYPE


  

  
CONTAINS

SUBROUTINE free_cachevals(v)
  TYPE(VALCACHE), INTENT(INOUT) :: v
  
  INTEGER :: i
  
  IF(ASSOCIATED(v%c)) THEN
    DO i = 1, v%ncache
      IF(ASSOCIATED(v%c(i)%val_)) DEALLOCATE(v%c(i)%val_)
    END DO
    IF(ASSOCIATED(v%c))DEALLOCATE(v%c)
    NULLIFY(v%c) !redundant 
  END IF
  v%ncache = 0
END SUBROUTINE


  !recursive loop over all cubature points for the given (m,mi) cache entry:
  !add each point to the buffer buf, evaluating all at once whenever the
  !buffer is full or when we are done
  
INTEGER RECURSIVE FUNCTION compute_cacheval(m, mi, val_, vali, fdim, fdata, &
                            & ndim, id_, p, xmin, xmax, buf, nbuf, ibuf, f, fv) RESULT(ret)
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
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
  END INTERFACE
  INTEGER, DIMENSION(:), INTENT(IN) :: m
  INTEGER, INTENT(IN) :: mi, fdim, ndim, id_, nbuf
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: buf, val_
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: p
  
  INTEGER :: vali, ibuf
  REAL(KIND=DBL) :: c, r
  !REAL(KIND=DBL), DIMENSION(:), POINTER :: x
  INTEGER :: x_loc
  INTEGER :: i, nx
  
  IF(PRESENT(f).EQV.PRESENT(fv)) THEN  !pass either f or fv
    ret = 1
    RETURN
  END IF
  
  IF(id_.EQ.ndim) THEN !add point to buffer
    DO i = 1,ndim
      buf(i+ibuf*ndim) = p(i)
    END DO
    ibuf = ibuf + 1
    IF(ibuf.EQ.nbuf) THEN !flush buffer
      IF(PRESENT(fv)) THEN
        IF(fv(ndim, nbuf, buf, fdata, fdim, val_(1+vali:)).NE.0) THEN
          ret = 1
          RETURN
        END IF
      ELSE IF(PRESENT(f)) THEN
        IF(fv_wrapper(f, ndim, nbuf, buf, fdata, fdim, val_(1+vali:)).NE.0) THEN
          ret = 1
          RETURN
        END IF
      END IF
      vali = vali+ ibuf * fdim
      ibuf = 0
    END IF
  ELSE
  
    c = (xmin(1+id_) + xmax(1+id_))*0.5_DBL
    r = (xmax(1+id_) - xmin(1+id_))*0.5_DBL
    IF(id_.EQ.mi) THEN
      IF(m(1+id_).NE.0) THEN
        !x => clencurt_x(1+ISHFT(1,m(1+id_)-1):)
        x_loc = ISHFT(1,m(1+id_)-1)
      ELSE
        !x => clencurt_x(1:)
        x_loc = 0
      END IF
    ELSE
      !x => clencurt_x(1:)
      x_loc = 0
    END IF
  
    IF(id_.EQ.mi) THEN
      IF(m(1+id_).NE.0) THEN
        nx = ISHFT(1, m(1+id_)-1)
      ELSE
        nx = 1
      END IF
    ELSE
      nx = ISHFT(1, m(1+id_))
    END IF
  
    IF(id_.NE.mi) THEN
      p(1+id_) = c
      IF(PRESENT(fv)) THEN
        IF(compute_cacheval(m, mi, val_, vali, fdim, fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, fv = fv).NE.0) THEN
          ret = 1
          RETURN
        END IF
      ELSE IF(PRESENT(f)) THEN
        IF(compute_cacheval(m, mi, val_, vali, fdim, fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, f = f).NE.0) THEN
          ret = 1
          RETURN
        END IF
      END IF
    END IF
    
    IF(PRESENT(fv)) THEN
      DO i = 1, nx
        !p(1+id_) = c + r * x(1+i)
        p(1+id_) = c + r * clencurt_x(x_loc+i)
        IF(compute_cacheval(m, mi, val_, vali, fdim, fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, fv = fv).NE.0) THEN
          ret = 1
          RETURN
        END IF
        !p(1+id_) = c - r * x(1+i)
        p(1+id_) = c - r * clencurt_x(x_loc+i)
        IF(compute_cacheval(m, mi, val_, vali, fdim,fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, fv = fv).NE.0) THEN
          ret = 1
          RETURN
        END IF
      END DO
    ELSE IF(PRESENT(f)) THEN
      DO i = 1, nx
        !p(1+id_) = c + r * x(1+i)
        p(1+id_) = c + r * clencurt_x(x_loc+i)
        IF(compute_cacheval(m, mi, val_, vali, fdim, fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, f = f).NE.0) THEN
          ret = 1
          RETURN
        END IF
        !p(1+id_) = c - r * x(1+i)
        p(1+id_) = c - r * clencurt_x(x_loc+i)
        IF(compute_cacheval(m, mi, val_, vali, fdim, fdata, ndim, id_+1, p, xmin, xmax, buf, nbuf, ibuf, f = f).NE.0) THEN
          ret = 1
          RETURN
        END IF
      END DO
    END IF
  END IF
  ret = 0
  RETURN
END FUNCTION

INTEGER FUNCTION num_cacheval(m, mi, ndim)
  INTEGER, DIMENSION(:), INTENT(IN) :: m
  INTEGER, INTENT(IN) :: mi, ndim
  
  INTEGER :: i, nval
  
  nval = 1
  DO i = 1,ndim
    IF(i.EQ.(mi+1)) THEN
      IF(m(i).EQ.0) THEN
        nval = nval * 2
      ELSE
        nval = nval * ISHFT(1, m(i))
      END IF
    ELSE
      nval = nval * (ISHFT(1, m(i)+1) + 1)
    END IF
  END DO
  num_cacheval = nval
END FUNCTION
  
INTEGER FUNCTION add_cacheval(vc, m, mi, fdim, fdata, ndim, xmin, xmax, buf, nbuf, f, fv)
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
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
  END INTERFACE
  TYPE(VALCACHE), INTENT(INOUT) :: vc
  INTEGER, DIMENSION(:), INTENT(IN) :: m
  INTEGER, INTENT(IN) :: mi, fdim, ndim, nbuf
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: buf
  
  INTEGER :: i, ic, nval, vali, ibuf, ncache_old
  REAL(KIND=DBL), DIMENSION(MAXDIM) :: p
  TYPE(CACHEVAL), DIMENSION(:), POINTER :: c_tmp
  
  IF(PRESENT(f).EQV.PRESENT(fv)) THEN !pass either f or fv
    add_cacheval = 1
    RETURN
  END IF
  
  ic = vc%ncache
  ncache_old = vc%ncache
  vali = 0; ibuf = 0
  vc%ncache = vc%ncache + 1
  
  ALLOCATE(c_tmp(vc%ncache))
  DO i = 1, MIN(ncache_old, vc%ncache)
    c_tmp(i) = vc%c(i)
  END DO
  IF(ASSOCIATED(vc%c)) DEALLOCATE(vc%c)
  vc%c => c_tmp
  NULLIFY(c_tmp)
  
  IF(.NOT.ASSOCIATED(vc%c)) THEN
    add_cacheval = -1
    RETURN
  END IF
  
  vc%c(1+ic)%mi = mi
  DO i = 1, ndim
    vc%c(1+ic)%m(i) = m(i)
  END DO
  
  nval = fdim * num_cacheval(m, mi, ndim)
  ALLOCATE(vc%c(1+ic)%val_(nval))
  IF(.NOT.ASSOCIATED(vc%c(1+ic)%val_)) THEN
    add_cacheval = 1
    RETURN
  END IF
  
  IF(PRESENT(fv)) THEN
    IF(compute_cacheval(m, mi, vc%c(1+ic)%val_, vali, &
      & fdim, fdata, ndim, 0, p, xmin, xmax, buf, nbuf, ibuf, fv = fv).NE.0) THEN
      add_cacheval = 1
      RETURN
    END IF
  ELSE IF(PRESENT(f)) THEN
    IF(compute_cacheval(m, mi, vc%c(1+ic)%val_, vali, &
      & fdim, fdata, ndim, 0, p, xmin, xmax, buf, nbuf, ibuf, f = f).NE.0) THEN
      add_cacheval = 1
      RETURN
    END IF
  END IF
  
  IF(ibuf.GT.0) THEN !flush remaining buffer
    IF(PRESENT(fv)) THEN
      IF(fv(ndim, ibuf, buf, fdata, fdim, vc%c(1+ic)%val_(1+vali:)).NE.0) THEN
        add_cacheval = 1
        RETURN
      END IF
    ELSE IF(PRESENT(f)) THEN
      IF(fv_wrapper(f, ndim, ibuf, buf, fdata, fdim, vc%c(1+ic)%val_(1+vali:)).NE.0) THEN
        add_cacheval = 1
        RETURN
      END IF
    END IF
  END IF
  
  add_cacheval = 0
END FUNCTION

  !recursive loop to evaluate the integral contribution from the cache
  !entry c, accumulating in val, for the given m[] except with m[md]
  !-> m[md] - 1 if md < dim, using the cached values (cm,cmi,cval).  id is the
  !current loop dimension (from 0 to dim-1).
  
  
INTEGER RECURSIVE FUNCTION eval(cm, cmi, cval, m, md, fdim, ndim, id_, weight, val_) RESULT(ret)
  INTEGER, DIMENSION(:), INTENT(IN) :: cm, m
  INTEGER, INTENT(IN) :: cmi, md, fdim, ndim, id_
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: cval
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: val_
  REAL(KIND=DBL), INTENT(IN) :: weight
  
  INTEGER :: voff, i, mid, cnx, nx
  !REAL(KIND=DBL), DIMENSION(:), POINTER :: w
  INTEGER :: w_loc
  
  voff = 0 !amount caller should offset cval array afterwards
  
  IF(id_.EQ.ndim) THEN
    DO i = 1, fdim
      val_(i) = val_(i) + cval(i) * weight
    END DO
    voff = fdim
  ELSE IF((m(1+id_).EQ.0).AND.(id_.EQ.md)) THEN !using trivial rule for this dimension
    voff = eval(cm, cmi, cval, m, md, fdim, ndim, id_+1, weight*2._DBL, val_)
    voff = voff + fdim * ISHFT(1, cm(1+id_)) * 2 * num_cacheval(cm + id_+1, cmi-(id_+1), ndim-(id_+1))
  ELSE
    mid = m(1+id_)
    IF(id_.EQ.md) THEN
      mid = mid-1
    END IF
    
    IF(id_.EQ.cmi) THEN
      IF(cm(1+id_).NE.0) THEN
        !w => clencurt_w(1+mid+ISHFT(1, mid) + ISHFT(1, cm(1+id_)-1):)
        w_loc = mid+ISHFT(1, mid) + ISHFT(1, cm(1+id_)-1)
      ELSE
        !w => clencurt_w(1+mid+ISHFT(1, mid)+1:)
        w_loc = mid+ISHFT(1, mid)+1
      END IF
    
    ELSE
      !w => clencurt_w(1+mid+ISHFT(1, mid)-1:)
      w_loc = mid+ISHFT(1, mid)-1
    END IF
    
    IF(id_.EQ.cmi) THEN
      IF(cm(1+id_).NE.0) THEN
        cnx = ISHFT(1, cm(1+id_)-1)
      ELSE
        cnx = 1
      END IF
    ELSE
      cnx = ISHFT(1, cm(1+id_))
    END IF
    
    IF(cm(1+id_).LE.mid) THEN
      nx = cnx
    ELSE
      nx = ISHFT(1, mid)
    END IF
    
    IF(id_.NE.cmi) THEN
      !voff = eval(cm, cmi, cval, m, md, fdim, ndim, id_ + 1, weight * w(1), val_)
      voff = eval(cm, cmi, cval, m, md, fdim, ndim, id_ + 1, weight * clencurt_w(1+w_loc), val_)
      !w => w(2:) !pointer arithmetic not allowed in Fortran
      w_loc = w_loc + 1
    END IF
    
    DO i = 1, nx
      !voff = voff + eval(cm, cmi, cval(1+voff:), m, md, fdim, ndim, id_ + 1, weight * w(1+i), val_)
      !voff = voff + eval(cm, cmi, cval(1+voff:), m, md, fdim, ndim, id_ + 1, weight * w(1+i), val_)
      voff = voff + eval(cm, cmi, cval(1+voff:), m, md, fdim, ndim, id_ + 1, weight * clencurt_w(i+w_loc), val_)
      voff = voff + eval(cm, cmi, cval(1+voff:), m, md, fdim, ndim, id_ + 1, weight * clencurt_w(i+w_loc), val_)
      
    END DO
    
    voff = voff + (cnx - nx) * fdim * 2 * num_cacheval(cm(1+id_+1:), cmi - (id_+1), ndim - (id_+1))
    
  END IF
  
  ret = voff
  RETURN

END FUNCTION

! loop over all cache entries that contribute to the integral, (with m[md] decremented by 1)
SUBROUTINE evals(vc, m, md, fdim, ndim, V, val_)
  TYPE(valcache), INTENT(IN) :: vc
  INTEGER, DIMENSION(:), INTENT(IN) :: m
  INTEGER, INTENT(IN) :: md, fdim, ndim
  REAL(KIND=DBL), INTENT(IN) :: V
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_
  
  INTEGER :: i, check, t
  
  DO i = 1, fdim
    val_(i) = 0._DBL
  END DO
  
  DO i = 1, vc%ncache
    IF(vc%c(i)%mi.EQ.md) THEN
      check = 1
    ELSE
      check = 0
    END IF
    
    IF((vc%c(i)%mi.GE.ndim).OR.((vc%c(i)%m(1+vc%c(i)%mi) + check).LE.m(1+vc%c(i)%mi))) THEN
      t = eval(vc%c(i)%m, vc%c(i)%mi, vc%c(i)%val_, m, md, fdim, ndim, 0, V, val_) !Dummy variable to just call the function
    END IF   
  END DO
END SUBROUTINE

! evaluate the integrals for the given m[] using the cached values in vc,
! storing the integrals in val[], the error estimate in err[], and the
! dimension to subdivide next (the largest error contribution) in *mi

SUBROUTINE eval_integral(vc, m, fdim, ndim, V, mi, val_, err_, val1)
  TYPE(valcache), INTENT(IN) :: vc
  INTEGER, DIMENSION(:), INTENT(IN) :: m
  INTEGER, INTENT(IN) :: fdim, ndim
  REAL(KIND=DBL), INTENT(IN) :: V
  INTEGER, INTENT(OUT) :: mi
  REAL(KIND=DBL), DIMENSION(:) :: val_, err_, val1
  
  INTEGER :: i, j
  REAL(KIND=DBL) :: maxerr, emax, e
  
  maxerr = 0._DBL
  CALL evals(vc, m, ndim, fdim, ndim, V, val_)
  ! error estimates along each dimension by comparing val with
	! lower-order rule in that dimension; overall (conservative)
	! error estimate from maximum error of lower-order rules.
  
  DO i = 1, fdim
    err_(i) = 0
  END DO
  mi = 0
  DO i = 0, ndim - 1
    emax = 0._DBL
    CALL evals(vc, m, i, fdim, ndim, V, val1)
    DO j = 1, fdim
      e = ABS(val_(j) - val1(j))
      IF(e.GT.emax) emax = e
      IF(e.GT.err_(j)) err_(j) = e
    END DO
    
    IF(emax.GT.maxerr) THEN
      maxerr = emax
      mi = i
    END IF
  END DO
END SUBROUTINE

INTEGER FUNCTION converged(fdim, vals, errs, reqAbsError, reqRelError, norm)
  INTEGER, INTENT(IN) :: fdim, norm
  REAL(KIND=DBL), INTENT(IN) ::  reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: vals, errs
  
  INTEGER :: j
  REAL(KIND=DBL) :: maxerr, serr, err_, maxval_, sval, val_, absval
  !Body of convergence test
  
  SELECT CASE(norm)
    CASE(1) !Error individual
      DO j = 1, fdim
        IF((errs(j).GT.reqAbsError).AND.(errs(j).GT.(ABS(vals(j))*reqRelError))) THEN
          converged = 0
          RETURN
        END IF
      END DO
      converged = 1
      RETURN
      
    CASE(2) !Error paired
      DO j = 1, fdim -1, 2
        !Scale to avoid overflow/underflow
        IF(errs(j).GT.errs(j+1)) THEN
          maxerr = errs(j)
        ELSE
          maxerr = errs(j+1)
        END IF
        IF(vals(j).GT.vals(j+1)) THEN
          maxval_ = vals(j)
        ELSE
          maxval_ = vals(j+1)
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
        
        err_ = SQRT( (errs(j) * serr)**2 + (errs(j+1) * serr)**2 ) * maxerr
        val_ = SQRT( (vals(j) * sval)**2 + (vals(j+1) * sval)**2 ) * maxval_
        
        IF( (err_.GT.reqAbsError).AND.(err_.GT.(val_*reqRelError))) THEN
          converged = 0
          RETURN
        END IF
        
      END DO
      !check for case that fdim is odd
      IF(j.LT.fdim) THEN
        IF( (errs(j).GT.reqAbsError).AND.(errs(j).GT.(ABS(vals(j))*reqRelError))) THEN
          converged = 0
          RETURN
        END IF
      END IF
      converged = 1
      RETURN
    
    CASE(3) !Error L1
      err_ = 0.; val_ = 0.
      DO j = 1, fdim
        err_ = err_ + errs(j)
        val_ = val_ + ABS(vals(j))
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
        absval = ABS(vals(j))
        IF(errs(j).GT.err_) err_ = errs(j)
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
        IF(errs(j).GT.maxerr) maxerr = errs(j)
        absval = ABS(vals(j))
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
        err_ = err_ + (errs(j) * serr)**2
        val_ = val_ + (vals(j) * sval)**2
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



  !Vectorized version with user-supplied buffer to store points and values.
  !The buffer *buf should be of length *nbuf * dim on entry (these parameters
  !are changed upon return to the final buffer and length that was used).
  !The buffer length will be kept <= max(max_nbuf, 1) * dim.

  !Also allows the caller to specify an array m[dim] of starting degrees
  !for the rule, which upon return will hold the final degrees.  The
  !number of points in each dimension i is 2^(m[i]+1) + 1.

INTEGER FUNCTION pcubature_v_buf(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, &
                                 & norm, m, buf, nbuf, max_nbuf, val_, err_, f, fv) RESULT(ret)
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
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, norm, maxEval, max_nbuf
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  INTEGER, INTENT(INOUT) :: nbuf
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  INTEGER, DIMENSION(:), INTENT(INOUT) :: m
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: buf, val_, err_
  
  INTEGER :: max_nbuf_new, norm_new, i, numEval, new_nbuf, use_new_buf, mi
  REAL(KIND=DBL) :: V
  REAL(KIND=DBL), DIMENSION(:), POINTER :: val1, buf_new
  TYPE(VALCACHE) :: vc
  
  !Major difference; we treat buf as an allocated array, not a pointer
  !If we need to create a buffer, then we use new_buf, which will happen by default if pcubature_v is called
  !If we use the new buffer, we set use_new_buf = 0; else we keep it at 1 and use the provided buffer
  !Unfortunately in Fortan 90, allocatable arrays cannot be used as dummy arguments
  
  ret = 1; V = 1._DBL; numEval = 0; norm_new = norm; max_nbuf_new = max_nbuf; vc%ncache = 0; use_new_buf = 1
  NULLIFY(val1, vc%c, buf_new)
  
  IF(PRESENT(f).EQV.PRESENT(fv)) THEN !pass either f or fv
    ret = 1
    RETURN
  END IF
  
  IF(fdim.LE.1) norm_new = 1 !norm is irrelevant
  IF((norm_new.LT.0).OR.(norm_new.GT.5)) RETURN !invalid norm
  IF(fdim.EQ.0) THEN !nothing to do
    ret = 0
    RETURN
  END IF 
  
  IF(ndim.GT.MAXDIM) RETURN !unsupported
  IF(ndim.EQ.0) THEN !trivial
    IF(PRESENT(fv)) THEN
      IF(fv(0, 1, xmin, fdata, fdim, val_).NE.0) RETURN 
    ELSE IF(PRESENT(f)) THEN
      IF(fv_wrapper(f, 0, 1, xmin, fdata, fdim, val_).NE.0) RETURN 
    END IF
    DO i = 1, fdim
      err_(i) = 0._DBL
    END DO
    ret = 0
    RETURN
  END IF
  
  DO i = 1, fdim
    val_(i) = 0._DBL
    err_(i) = HUGE(0._DBL)
  END DO
  
  DO i = 1, ndim
    V = V * (xmax(i) - xmin(i)) * 0.5_DBL !scale factor for C-C volume
  END DO
  
  new_nbuf = num_cacheval(m, ndim, ndim)
  IF(max_nbuf.LT.1) max_nbuf_new = 1
  IF(new_nbuf.GT.max_nbuf_new) new_nbuf = max_nbuf_new
  IF(nbuf.LT.new_nbuf) THEN
    IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
    ALLOCATE(buf_new(new_nbuf*ndim))
    nbuf = new_nbuf
    use_new_buf = 0
    IF(.NOT.ASSOCIATED(buf_new)) THEN !failed
      IF(ASSOCIATED(val1)) DEALLOCATE(val1)
      IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
      CALL free_cachevals(vc)
      RETURN
    END IF
  END IF
  
  !start by evaluating the m=0 cubature rule
  IF(PRESENT(fv)) THEN
    IF(use_new_buf.EQ.1) THEN !use old buffer
      IF(add_cacheval(vc, m, ndim, fdim, fdata, ndim, xmin, xmax, buf, nbuf, fv = fv).NE.0) THEN !failed
        IF(ASSOCIATED(val1)) DEALLOCATE(val1)
        IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
        CALL free_cachevals(vc)
        RETURN
      END IF
    ELSE IF(use_new_buf.EQ.0) THEN!use new buffer
      IF(add_cacheval(vc, m, ndim, fdim, fdata, ndim, xmin, xmax, buf_new, nbuf, fv = fv).NE.0) THEN !failed
        IF(ASSOCIATED(val1)) DEALLOCATE(val1)
        IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
        CALL free_cachevals(vc)
        RETURN
      END IF
    END IF
  ELSE IF(PRESENT(f)) THEN
    IF(use_new_buf.EQ.1) THEN !use old buffer
      IF(add_cacheval(vc, m, ndim, fdim, fdata, ndim, xmin, xmax, buf, nbuf, f = f).NE.0) THEN !failed
        IF(ASSOCIATED(val1)) DEALLOCATE(val1)
        IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
        CALL free_cachevals(vc)
        RETURN
      END IF
    ELSE IF(use_new_buf.EQ.0) THEN!use new buffer
      IF(add_cacheval(vc, m, ndim, fdim, fdata, ndim, xmin, xmax, buf_new, nbuf, f = f).NE.0) THEN !failed
        IF(ASSOCIATED(val1)) DEALLOCATE(val1)
        IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
        CALL free_cachevals(vc)
        RETURN
      END IF
    END IF
  END IF
  
  ALLOCATE(val1(fdim))
  
  DO WHILE(.TRUE.) !CAREFUL, intentional infinite loop
    CALL eval_integral(vc, m, fdim, ndim, V, mi, val_, err_, val1)
    IF((converged(fdim, val_, err_, reqAbsError, reqRelError, norm).NE.0).OR.((numEval.GT.maxEval).AND.(maxEval.NE.0))) THEN
      ret = 0 !success
      IF(ASSOCIATED(val1)) DEALLOCATE(val1)
      IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
      CALL free_cachevals(vc)
      RETURN
    END IF
    
    m(1+mi) = m(1+mi) + 1
    IF(m(1+mi).GT.clencurt_M) THEN !failed
      IF(ASSOCIATED(val1)) DEALLOCATE(val1)
      IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
      CALL free_cachevals(vc)
      RETURN
    END IF
    
    new_nbuf = num_cacheval(m, mi, ndim)
    IF((new_nbuf.GT.nbuf).AND.(nbuf.LT.max_nbuf)) THEN
      nbuf = new_nbuf
      IF(nbuf.GT.max_nbuf) nbuf = max_nbuf
      IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
      ALLOCATE(buf_new(nbuf*ndim))
      use_new_buf = 0
      IF(.NOT.ASSOCIATED(buf_new)) THEN !failed
        IF(ASSOCIATED(val1)) DEALLOCATE(val1)
        IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
        CALL free_cachevals(vc)
        RETURN
      END IF
    END IF
    IF(PRESENT(fv)) THEN
      IF(use_new_buf.EQ.1) THEN !use old buffer
        IF(add_cacheval(vc, m, mi, fdim, fdata, ndim, xmin, xmax, buf, nbuf, fv = fv).NE.0) THEN !failure
          IF(ASSOCIATED(val1)) DEALLOCATE(val1)
          IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
          CALL free_cachevals(vc)
          RETURN
        END IF
      ELSE IF(use_new_buf.EQ.0) THEN !use new buffer
        IF(add_cacheval(vc, m, mi, fdim, fdata, ndim, xmin, xmax, buf_new, nbuf, fv = fv).NE.0) THEN !failure
          IF(ASSOCIATED(val1)) DEALLOCATE(val1)
          IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
          CALL free_cachevals(vc)
          RETURN
        END IF
      END IF
    ELSE IF(PRESENT(f)) THEN
      IF(use_new_buf.EQ.1) THEN !use old buffer
        IF(add_cacheval(vc, m, mi, fdim, fdata, ndim, xmin, xmax, buf, nbuf, f = f).NE.0) THEN !failure
          IF(ASSOCIATED(val1)) DEALLOCATE(val1)
          IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
          CALL free_cachevals(vc)
          RETURN
        END IF
      ELSE IF(use_new_buf.EQ.0) THEN !use new buffer
        IF(add_cacheval(vc, m, mi, fdim, fdata, ndim, xmin, xmax, buf_new, nbuf, f = f).NE.0) THEN !failure
          IF(ASSOCIATED(val1)) DEALLOCATE(val1)
          IF(ASSOCIATED(buf_new)) DEALLOCATE(buf_new)
          CALL free_cachevals(vc)
          RETURN
        END IF
      END IF
    END IF
    numEval = numEval + new_nbuf   
  END DO


END FUNCTION

INTEGER FUNCTION pcubature_v(fdim, fv, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, fdata) RESULT(ret)
  INTERFACE
    INTEGER FUNCTION fv(ndim, npt, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim, npt
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim*npts
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !fdim*npts
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, maxEval, norm
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  REAL(KIND=DBL), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(0) :: fdata_new
  
  INTEGER :: nbuf
  INTEGER, DIMENSION(MAXDIM) :: m
  REAL(KIND=DBL), DIMENSION(0) :: buf !0 sized buffer since we're letting the internal program handle the buffer 
  
  m = 0
  nbuf = 0
  
  !note buf_new is deallocated in the internal program, so no need to do so here
  IF(PRESENT(fdata)) THEN
    ret = pcubature_v_buf(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, m, buf, nbuf, &
    & DEFAULT_MAX_NBUF, val_, err_, fv = fv)
  ELSE
    ret = pcubature_v_buf(fdim, fdata_new, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, m, buf, nbuf, &
    & DEFAULT_MAX_NBUF, val_, err_, fv = fv)
  END IF
END FUNCTION
  
INTEGER FUNCTION pcubature(fdim, f, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, val_, err_, fdata) RESULT(ret)
  INTERFACE
    INTEGER FUNCTION f(ndim, x, fdata, fdim, fval)
      USE KIND
      INTEGER, INTENT(IN) :: ndim, fdim
      REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !ndim
      REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval !fdim
    END FUNCTION
  END INTERFACE
  INTEGER, INTENT(IN) :: fdim, ndim, maxEval, norm
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: xmin, xmax
  REAL(KIND=DBL), INTENT(IN) :: reqAbsError, reqRelError
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: val_, err_
  REAL(KIND=DBL), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(0) :: fdata_new
  
  INTEGER :: nbuf
  INTEGER, DIMENSION(MAXDIM) :: m
  REAL(KIND=DBL), DIMENSION(0) :: buf !0 sized buffer since we're letting the internal program handle the buffer
  
  !max_nbuf > 0 to amortize function overhead
  m = 0
  nbuf = 0 
  
  !note buf_new is deallocated in the internal program, so no need to do so here
  IF(PRESENT(fdata)) THEN
    ret = pcubature_v_buf(fdim, fdata, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, m, buf, nbuf, 16, val_, &
    & err_, f = f)
  ELSE
    ret = pcubature_v_buf(fdim, fdata_new, ndim, xmin, xmax, maxEval, reqAbsError, reqRelError, norm, m, buf, nbuf, 16, val_, &
    & err_, f = f)
  END IF

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
  INTEGER, INTENT(IN) :: ndim, npt, fdim
  REAL(KIND=DBL), DIMENSION(:), INTENT(IN) :: x !npt * ndim
  REAL(KIND=DBL), DIMENSION(:), INTENT(INOUT) :: fdata
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: fval  !npt * fdim
  
  INTEGER :: i
  
  DO i = 0, npt-1
    IF(f(ndim, x(1+i*ndim:), fdata, fdim, fval(1+i*fdim:)).NE.0) THEN
      fv_wrapper = 1
      RETURN
    END IF
  END DO
  fv_wrapper = 0
END FUNCTION
  

END MODULE PCUB
  
  
  
  
  
  
  





