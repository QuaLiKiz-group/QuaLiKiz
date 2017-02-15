MODULE diskio
  USE kind
  IMPLICIT NONE
  !INCLUDE 'mpif.h'
  INTERFACE writevar
     MODULE PROCEDURE &
          writevar_0d, &
          writevar_0d_integer, &
          writevar_1d, &
          writevar_2d, &
          writevar_2d_complex, &
          writevar_2d_integer, &
          writevar_3d, &
          writevar_3d_complex, &
          writevar_4d, &
          writevar_4d_complex
  END INTERFACE

  INTERFACE readvar
     MODULE PROCEDURE &
          readvar_0d, &
          readvar_1d, &
          readvar_2d
  END INTERFACE


CONTAINS
  SUBROUTINE open_file_in_bin(filename,lengthin,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER:: lengthin
    INTEGER, INTENT(IN) :: myunit

    OPEN( unit=myunit, &
         file=filename, &
         access='direct', &
         form='unformatted', &
         status='old', &
         recl=lengthin)
  END SUBROUTINE open_file_in_bin

  SUBROUTINE open_file_out_txt(filename,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: myunit
    INTEGER :: rc

    OPEN(unit=myunit, &
         file=filename, &
         form='formatted', &
         status='replace', &
         iostat=rc)
  END SUBROUTINE open_file_out_txt

  SUBROUTINE extract_filename(filename,dirname,basename,suffix)
    CHARACTER(len=*), INTENT(IN) :: filename
    CHARACTER(len=:), INTENT(OUT), ALLOCATABLE:: dirname,basename,suffix
    INTEGER :: dirsep, sufsep

    dirsep=index(filename, '/', BACK=.TRUE.)
    sufsep=index(filename, '.', BACK=.TRUE.)

    dirname=filename(1:dirsep)
    basename=filename(dirsep+1:sufsep-1)
    suffix=filename(sufsep+1:)
  END SUBROUTINE extract_filename

  SUBROUTINE writevar_0d(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    REAL, INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit

    CALL open_file_out_txt(filename,myunit)
    WRITE(unit=myunit,fmt='(' // varformat // ')') mold
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_0d

  SUBROUTINE writevar_0d_integer(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    INTEGER, INTENT(IN) :: mold, myunit

    CALL open_file_out_txt(filename,myunit)
    WRITE(unit=myunit,fmt='(' // varformat // ')') mold
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_0d_integer

  SUBROUTINE writevar_1d(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    REAL, DIMENSION(:), INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit
    INTEGER           :: numcols
    CHARACTER(LEN=30) :: rowfmt 

    numcols = 1 !Write out 1D array as a single column
    CALL open_file_out_txt(filename,myunit)
    WRITE(rowfmt,'(A,I0,A)') '(',numcols,varformat // ')'
    WRITE(unit=myunit,fmt=rowfmt) mold
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_1d

  SUBROUTINE writevar_2d(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    REAL, DIMENSION(:,:), INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit
    INTEGER           :: numcols,numrows,i,j
    CHARACTER(LEN=30) :: rowfmt 

    numcols = SIZE(mold,2)
    numrows = SIZE(mold,1)
    CALL open_file_out_txt(filename,myunit)
    WRITE(rowfmt,'(A,I0,A)') '(',numcols,varformat // ')'
    WRITE(myunit,rowfmt) ((mold(i,j),j=1,numcols),i=1,numrows)
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_2d

  SUBROUTINE writevar_2d_integer(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    INTEGER, DIMENSION(:,:), INTENT(IN) :: mold
    INTEGER           :: numcols,numrows,i,j
    INTEGER, INTENT(IN) :: myunit
    CHARACTER(LEN=30) :: rowfmt 

    numcols = SIZE(mold,2)
    numrows = SIZE(mold,1)
    CALL open_file_out_txt(filename,myunit)
    WRITE(rowfmt,'(A,I0,A)') '(',numcols, varformat // ')'
    WRITE(myunit,rowfmt) ((mold(i,j),j=1,numcols),i=1,numrows)
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_2d_integer

  SUBROUTINE writevar_2d_complex(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    CHARACTER(len=len(filename)+1) :: filename_real, filename_imag
    COMPLEX(kind=DBL), DIMENSION(:,:), INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit
    REAL, DIMENSION(SIZE(mold)) :: mold_real, mold_imag
    CHARACTER(LEN=30) :: rowfmt
    CHARACTER(LEN=:), ALLOCATABLE :: dirname,basename,suffix

    CALL extract_filename(filename,dirname,basename,suffix)
    CALL writevar_2d(dirname // 'r' // basename // '.' // suffix, &
         REAL(mold),varformat,myunit)
    CALL writevar_2d(dirname // 'i' // basename // '.' // suffix, &
         AIMAG(mold),varformat,myunit)
  END SUBROUTINE writevar_2d_complex

  SUBROUTINE writevar_3d(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    REAL, DIMENSION(:,:,:), INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit
    INTEGER           :: numcols,numrows,numpages,i,j,k
    CHARACTER(LEN=30) :: rowfmt

    numpages = SIZE(mold,3)
    numcols = SIZE(mold,2)
    numrows = SIZE(mold,1)
    CALL open_file_out_txt(filename,myunit)
    WRITE(rowfmt,'(A,I0,A)') '(',numcols,varformat // ')'
    WRITE(myunit,rowfmt) (((mold(i,j,k),j=1,numcols),i=1,numrows),k=1,numpages)
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_3d

  SUBROUTINE writevar_3d_complex(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    INTEGER, INTENT(IN) :: myunit
    CHARACTER(len=len(filename)+1) :: filename_real, filename_imag
    COMPLEX(kind=DBL), DIMENSION(:,:,:), INTENT(IN) :: mold
    REAL, DIMENSION(SIZE(mold)) :: mold_real, mold_imag
    CHARACTER(LEN=30) :: rowfmt
    CHARACTER(LEN=:), ALLOCATABLE :: dirname,basename,suffix

    CALL extract_filename(filename,dirname,basename,suffix)
    CALL writevar_3d(dirname // 'r' // basename // '.' // suffix, &
         REAL(mold),varformat,myunit)
    CALL writevar_3d(dirname // 'i' // basename // '.' // suffix, &
         AIMAG(mold),varformat,myunit)
  END SUBROUTINE writevar_3d_complex

  SUBROUTINE writevar_4d(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    REAL, DIMENSION(:,:,:,:), INTENT(IN) :: mold
    INTEGER, INTENT(IN) :: myunit
    INTEGER           :: numcols,numrows,numpages,numhpages,i,j,k,l
    CHARACTER(LEN=30) :: rowfmt

    numhpages = SIZE(mold,4)
    numpages = SIZE(mold,3)
    numcols = SIZE(mold,2)
    numrows = SIZE(mold,1)
    CALL open_file_out_txt(filename,myunit)
    WRITE(rowfmt,'(A,I0,A)') '(',numcols,varformat // ')'
    WRITE(myunit,rowfmt) ((((mold(i,j,k,l),j=1,numcols),i=1,numrows), &
         k=1,numpages),l=1,numhpages)
    CLOSE(unit=myunit)
  END SUBROUTINE writevar_4d

  SUBROUTINE writevar_4d_complex(filename, mold, varformat,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename, varformat
    INTEGER, INTENT(IN) :: myunit
    CHARACTER(len=len(filename)+1) :: filename_real, filename_imag
    COMPLEX(kind=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: mold
    REAL, DIMENSION(SIZE(mold)) :: mold_real, mold_imag
    CHARACTER(LEN=30) :: rowfmt
    CHARACTER(LEN=:), ALLOCATABLE :: dirname,basename,suffix

    CALL extract_filename(filename,dirname,basename,suffix)
    CALL writevar_4d(dirname // 'r' // basename // '.' // suffix, &
         REAL(mold),varformat,myunit)
    CALL writevar_4d(dirname // 'i' // basename // '.' // suffix, &
         AIMAG(mold),varformat,myunit)
  END SUBROUTINE writevar_4d_complex

  REAL(kind=DBL) FUNCTION readvar_0d(filename, dummy, ktype,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: ktype, myunit
    REAL(kind=DBL), INTENT(IN) ::dummy
    INTEGER :: lengthin,istat
    INQUIRE(iolength=lengthin) dummy
    CALL open_file_in_bin(filename,lengthin,myunit)
    READ(unit=myunit,rec=1,iostat=istat) readvar_0d
    CLOSE(unit=myunit)
    IF ( istat /= 0 ) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       !CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END FUNCTION readvar_0d

  FUNCTION readvar_1d(filename,dummy,ktype,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: ktype, myunit
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) ::dummy
    INTEGER :: lengthin,istat
    REAL(KIND=DBL), DIMENSION(SIZE(dummy,1)) :: readvar_1d

    INQUIRE(iolength=lengthin) dummy
    CALL open_file_in_bin(filename,lengthin,myunit)
    READ(unit=myunit,rec=1,IOSTAT=istat) readvar_1d
    CLOSE(unit=myunit)
    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       !CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END FUNCTION readvar_1d

  FUNCTION readvar_2d(filename,dummy,ktype,myunit)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: ktype, myunit
    REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) ::dummy
    INTEGER :: lengthin,istat
    REAL(KIND=DBL), DIMENSION(SIZE(dummy,1), SIZE(dummy,2)) :: readvar_2d

    INQUIRE(iolength=lengthin) dummy
    CALL open_file_in_bin(filename,lengthin,myunit)
    READ(unit=myunit,rec=1,IOSTAT=istat) readvar_2d
    CLOSE(unit=myunit)
    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       !CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END FUNCTION readvar_2d
END MODULE diskio
