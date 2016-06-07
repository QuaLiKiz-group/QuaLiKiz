module nanfilter

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Function that should return true if the passed real is a NaN 
  !> Attempts to be robust to all compilers, arcitectures and optimisations
  !> but probably not as good as ieee_is_nan in 2003 standard
  !-----------------------------------------------------------------------------
  function gkw_is_nan(real_to_test)

    logical            :: gkw_is_nan
    logical            :: isnan1, isnan2, isnan3, isnan4, isnan5, isnan6, isnan7
    real, intent(in)   :: real_to_test
    character(len=128) :: mychar

    mychar='000'
    write(mychar,*) real_to_test
    mychar = adjustl(mychar)

    isnan1 = (trim(mychar)==trim(adjustl('NaN'))) 
    isnan2 = (trim(mychar)==trim(adjustl('NAN'))) 
    isnan3 = (trim(mychar)==trim(adjustl('nan'))) 
    isnan4 = (real_to_test /= real_to_test)
    isnan5 = ((real_to_test > 0.0) .EQV. (real_to_test <= 0.0))
    isnan6 = bittest_is_nan(real_to_test)

    !http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2004-05/0882.html
    if (exponent(real_to_test) > maxexponent(real_to_test)) then 
       if(fraction(real_to_test) == 0.5) then
          isnan7 = .false.  ! 'infinity' -> might want another function for this 
       else
          isnan7 = .true. 
       endif
    else 
       isnan7 = .false.
    endif

    gkw_is_nan = isnan1 .or. isnan2 .or. isnan3 .or. isnan4 .or. isnan5
    gkw_is_nan = gkw_is_nan .or. isnan6 .or. isnan7

!!$    if (gkw_is_nan) then
!!$       write(*,*) 'isnan1', isnan1
!!$       write(*,*) 'isnan2', isnan2
!!$       write(*,*) 'isnan3', isnan3
!!$       write(*,*) 'isnan4', isnan4
!!$       write(*,*) 'isnan5', isnan5
!!$       write(*,*) 'isnan6', isnan6
!!$       write(*,*) 'isnan7', isnan7
!!$    end if

  end function gkw_is_nan

  !-----------------------------------------------------------------------------
  !> bit test for NaN on selected real precision
  !>
  !> The function operates by transferring bit pattern from a real variable to 
  !> an integer container. The value is exclusive ORed with the value being 
  !> tested for  (the ieee NAN bit representation hardcoded in global.F90)
  !> The integer result of the IEOR function is converted to a logical result 
  !> by comparing it to zero.  Note, other NaNs may also exist which do NOT
  !> match the bit pattern, but these will not be ieee conforming.
  !>
  !> If the argument is array valued, the function returns a conformable logical
  !> array, suitable for use with the ANY function, or as a logical mask.
  !>
  !> This function (only) is adapted for GKW from infnan.f90 
  !> which is available from http://www.lahey.com/code.htm:
  !!
  !! Copyright(c) 2003, Lahey Computer Systems, Inc.
  !! Copies of this source code, or standalone compiled files 
  !! derived from this source may not be sold without permission
  !! from Lahey Computers Systems. All or part of this [function] may be 
  !! freely incorporated into executable programs which are offered
  !! for sale. Otherwise, distribution of all or part of this file is
  !! permitted, provided this copyright notice and header are included.
  !----------------------------------------------------------------------------
  elemental function bittest_is_nan(real_to_test) result(res)

    integer, parameter :: Double = selected_int_kind(precision(1.d0))
    integer(Double), parameter :: iNaN  = Z"7FF8000000000000"
    !data iNaN/Z"7FF8000000000000"/
    ! Position of the sign bit (Intel) bit numbering starts at zero
    integer, parameter :: PSB = bit_size(iNaN) - 1
    real, intent(in)   :: real_to_test
    logical            :: res

    res = ieor(ibclr(transfer(real_to_test,iNaN),PSB), iNaN) == 0

  end function bittest_is_nan


  !-----------------------------------------------------------------------------
  !> Take an integer, i, and write it into a character of length j. The total
  !> length is always 8, so it should usually be trimmed.
  !------------------------------------------------------------------------------
  function int2char(number,ndigits)

    integer, intent(in)           :: number
    integer, optional, intent(in) :: ndigits

    character (len=8) :: int2char
    character (len=6) :: ifmt

    if (present(ndigits)) then
       select case(ndigits)
       case(0) ; ifmt='(I4)'
       case(1) ; ifmt='(I1)'
       case(2) ; ifmt='(I2)'
       case(3) ; ifmt='(I3)'
       case(4) ; ifmt='(I4)'
       case(5) ; ifmt='(I5)'
       case(6) ; ifmt='(I6)'
       case(7) ; ifmt='(I7)'
       case default ; ifmt='(I8)'
       end select
       !write(int2char,'(A8)') '        '
    else
       ifmt='(I8)'
    end if

    write(int2char,ifmt) number

  end function int2char

end module nanfilter


