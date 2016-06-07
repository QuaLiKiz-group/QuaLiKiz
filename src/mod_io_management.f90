! ---------------------------------
! PURPOSE: INPUT/OUTPUT MANAGEMENT.
! ---------------------------------
MODULE mod_io_management
  USE kind
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER, PRIVATE :: unite
  CHARACTER(len=15), PRIVATE :: filename

  ! -----------------------------------------------------------
  ! GENERIC FUNCTION TO READ INPUT VARIABLES FROM BINARY FILES
  ! -----------------------------------------------------------
  INTERFACE readvar
     MODULE PROCEDURE readvar_0d, &
          readvar_1d, &
          readvar_2d, &
          readvar_3d, &
          readchar
  END INTERFACE

  ! -----------------------------------------------------------
  ! GENERIC FUNCTION TO WRITE OUTPUT VARIABLES TO BINARY FILES
  ! -----------------------------------------------------------
  INTERFACE writevar
     MODULE PROCEDURE writevar_0d, &
          writevar_1d, &
          writevar_2d, &
          writevar_3d
  END INTERFACE

  PRIVATE :: open_file_in_bin, open_file_out_bin, unitcounter

CONTAINS

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  INTEGER FUNCTION unitcounter()
    ! -------------------------------------------------------------------
    ! PURPOSE: INITIALIZE AND INCREMENT THE UNIT NUMBER FOR FILE OPENING
    ! -------------------------------------------------------------------
    LOGICAL :: value

    unitcounter=0
    DO
       INQUIRE(unit=unitcounter,opened=value )
       IF(.NOT.value) RETURN
       unitcounter = unitcounter + 1
    ENDDO
  END FUNCTION unitcounter

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE open_file_in_bin(kfile,lengthin)
    ! -------------------------------------------
    ! PURPOSE: OPEN BINARY FILES FOR DATA READING
    ! -------------------------------------------
    INTEGER, INTENT(IN) :: kfile
    INTEGER, INTENT(IN) :: lengthin
    INTEGER :: nco
    CHARACTER(len=3) :: icont

    WRITE(icont,'(i3)') kfile
    nco = SCAN(icont,"123456789")
    filename = 'input/p'//icont(nco:)//'.bin'

    IF(lengthin<0) THEN
       WRITE(stderr,*) '-----------------------------------------------------'
       WRITE(stderr,'(a39,a15)') ' The data to read is too large for file ',filename
       WRITE(stderr,*) '-----------------------------------------------------'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    unite = unitcounter()

    OPEN( unit=unite, &
         file=filename, &
         access='direct', &
         form='unformatted', &
         status='old', &
         recl=lengthin )
  END SUBROUTINE open_file_in_bin

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE open_file_in_txt(kfile,lengthin)
    ! -------------------------------------------
    ! PURPOSE: OPEN ASCII FILES FOR DATA READING
    ! -------------------------------------------
    INTEGER, INTENT(IN) :: kfile
    INTEGER, INTENT(IN) :: lengthin
    INTEGER:: nco,rc
    CHARACTER(len=3):: icont

    WRITE(icont,'(i3)') kfile
    nco = SCAN(icont,"123456789")
    filename = 'input/p'//icont(nco:)//'.txt'

    unite = unitcounter()

    OPEN( unit=unite, &
         file=filename, &
         form='formatted', &
         status='old', &
         iostat=rc)

  END SUBROUTINE open_file_in_txt

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE open_file_out_bin(kfile,lengthout)
    ! -------------------------------------------
    ! PURPOSE: OPEN BINARY FILES FOR DATA WRITING
    ! -------------------------------------------
    INTEGER, INTENT(IN) :: kfile
    INTEGER, INTENT(IN) :: lengthout
    INTEGER:: nco
    CHARACTER(len=3):: icont

    WRITE(icont,'(i3)') kfile
    nco = SCAN(icont,"123456789")
    filename = 'output/f'//icont(nco:)//'.bin'

    IF(lengthout<0) THEN
       WRITE(stderr,*) '-------------------------------------------------------'
       WRITE(stderr,'(a41,a15)') ' The data to export is too large for file ',filename
       WRITE(stderr,*) '-------------------------------------------------------'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    unite = unitcounter()
    OPEN( unit=unite, &
         file=filename, &
         access='direct', &
         form='unformatted', &
         status='replace', &
         recl=lengthout )
  END SUBROUTINE open_file_out_bin

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE open_file_out_txt(kfile,lengthout)
    ! -------------------------------------------
    ! PURPOSE: OPEN ASCII FILES FOR DATA WRITING
    ! -------------------------------------------
    INTEGER, INTENT(IN) :: kfile
    INTEGER, INTENT(IN) :: lengthout
    INTEGER:: nco,rc
    CHARACTER(len=3):: icont

    WRITE(icont,'(i3)') kfile
    nco = SCAN(icont,"123456789" )
    filename = 'output/f'//icont(nco:)//'.txt'

    unite = unitcounter()

    OPEN( unit=unite, &
         file=filename, &
         form='formatted', &
         status='replace', &
         iostat=rc)

  END SUBROUTINE open_file_out_txt

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  REAL(kind=DBL) FUNCTION readvar_0d(kfile,mold,kformat)
    ! ------------------------------------------------------------
    ! PURPOSE: IMPORT SCALAR VARIABLES COMING FROM BINARY FILES
    ! ------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE INPUT FILE
    !              (= P{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! ------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), INTENT(IN) :: mold
    INTEGER :: lengthin,istat=0

    IF(kformat==1) THEN ! BINARY FILES
    
       INQUIRE(iolength=lengthin) mold
       CALL open_file_in_bin(kfile,lengthin)
       READ(unit=unite,rec=1,IOSTAT=istat) readvar_0d

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_in_txt(kfile,lengthin)
       READ(unit=unite,fmt=*) readvar_0d

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0 ) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

  END FUNCTION readvar_0d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  FUNCTION readvar_1d(kfile,mold,kformat)
    ! ------------------------------------------------------------
    ! PURPOSE: IMPORT 1D VARIABLES COMING FROM BINARY FILES
    ! ------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE INPUT FILE
    !              (= P{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! ------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: mold
    INTEGER:: lengthin, istat=0
    REAL(KIND=DBL), DIMENSION(SIZE(mold,1)) :: readvar_1d

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthin) mold
      
       CALL open_file_in_bin(kfile,lengthin)
       READ( unit=unite,rec=1,IOSTAT=istat) readvar_1d

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_in_txt(kfile,lengthin)
       READ(unit=unite,fmt=*) readvar_1d

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

  END FUNCTION readvar_1d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  FUNCTION readvar_2d(kfile,mold,kformat)
    ! ------------------------------------------------------------
    ! PURPOSE: IMPORT 2D VARIABLES COMING FROM BINARY FILES
    ! ------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE INPUT FILE
    !              (= P{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! ------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) :: mold
    INTEGER:: lengthin,istat=0
    REAL(KIND=DBL), DIMENSION(SIZE(mold,1),SIZE(mold,2)) :: readvar_2d

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthin) mold
       CALL open_file_in_bin(kfile,lengthin)
       READ(unit=unite,rec=1,IOSTAT=istat) readvar_2d

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_in_txt(kfile,lengthin)
       READ(unit=unite,fmt=*) readvar_2d(1:SIZE(mold,1),1:SIZE(mold,2))

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END FUNCTION readvar_2d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  FUNCTION readvar_3d(kfile,mold,kformat)
    ! ------------------------------------------------------------
    ! PURPOSE: IMPORT 3D VARIABLES COMING FROM BINARY FILES
    ! ------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE INPUT FILE
    !              (= P{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! ------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:,:,:), INTENT(IN) :: mold
    INTEGER:: lengthin,istat=0
    REAL(KIND=DBL), DIMENSION(SIZE(mold,1),SIZE(mold,2), SIZE(mold,3)) :: readvar_3d

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthin) mold
       CALL open_file_in_bin(kfile,lengthin)
       READ( unit=unite,rec=1,IOSTAT=istat) readvar_3d

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_in_txt(kfile,lengthin)
       READ(unit=unite,fmt=*) readvar_3d(1:SIZE(mold,1),1:SIZE(mold,2),1:SIZE(mold,3))

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END FUNCTION readvar_3d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------


  FUNCTION readchar(kfile,mold,kformat)
    ! ------------------------------------------------------------
    ! PURPOSE: IMPORT CHARACTER VARIABLES COMING FROM BINARY FILES
    ! ------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE INPUT FILE
    !              (= P{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! ------------------------------------------------------------
    INTEGER,INTENT(IN) :: kfile,kformat
    CHARACTER(LEN=10), INTENT(IN) :: mold
    INTEGER:: lengthin,istat=0
    CHARACTER(LEN=10) :: readchar

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthin) mold
       CALL open_file_in_bin(kfile,lengthin)
       READ( unit=unite,rec=1,IOSTAT=istat) readchar

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_in_txt(kfile,lengthin)
       READ(unit=unite,fmt=*) readchar

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN READING INPUT FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

  END FUNCTION readchar

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE writevar_0d(kfile,mold,charname,kformat)
    ! -------------------------------------------------------------
    ! PURPOSE: EXPORT SCALAR VARIABLES TO BINARY FILES
    ! -------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE OUTPUT FILE
    !              (= F{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! -------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), INTENT(IN) :: mold
    CHARACTER(len=14) charname
    INTEGER           :: lengthout,istat=0

    ! FILL AN OUTPUT ASCII FILE WITH THE NAME AND THE DIMENSION OF THE VARIABLES
    WRITE(3311,fmt='(a18,8i10)') charname,kfile,SHAPE(mold)

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthout) mold
       CALL open_file_out_bin(kfile,lengthout)
       WRITE(unit=unite,rec=1,IOSTAT=istat) mold

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_out_txt(kfile,lengthout)
       WRITE(unit=unite,fmt=*) mold

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN WRITING 0D VARIABLE TO FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

  END SUBROUTINE writevar_0d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE writevar_1d(kfile,mold,charname,kformat)
    ! -------------------------------------------------------------
    ! PURPOSE: EXPORT 1D VARIABLES TO BINARY FILES
    ! -------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE OUTPUT FILE
    !              (= F{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! -------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:), INTENT(IN) :: mold
    INTEGER:: i
    CHARACTER(len=14) charname
    INTEGER:: lengthout,istat=0

    ! FILL AN OUTPUT ASCII FILE WITH THE NAME AND THE DIMENSION OF THE VARIABLES
    WRITE(3311,fmt='(a18,8i10)') charname,kfile,SHAPE(mold)

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthout) mold
       CALL open_file_out_bin(kfile,lengthout)
       WRITE(unit=unite,rec=1,IOSTAT=istat) mold

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_out_txt(kfile,lengthout)
       DO i=1,SIZE(mold,1)
          WRITE(unit=unite,fmt=*) mold(i)
       ENDDO

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN WRITING 1D VARIABLE TO FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END SUBROUTINE writevar_1d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE writevar_2d(kfile,mold,charname,kformat)
    ! -------------------------------------------------------------
    ! PURPOSE: EXPORT 2D VARIABLES TO BINARY FILES
    ! -------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE OUTPUT FILE
    !              (= F{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! -------------------------------------------------------------
    INTEGER, INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:,:), INTENT(IN) :: mold
    INTEGER:: i,nco
    CHARACTER(len=14) charname
    CHARACTER(len=3) chardim
    INTEGER:: lengthout,rw,istat=0

    ! FILL AN OUTPUT ASCII FILE WITH THE NAME AND THE DIMENSION OF THE VARIABLES
    WRITE(3311,fmt='(a18,8i10)') charname,kfile,SHAPE(mold)

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthout) mold
       CALL open_file_out_bin(kfile,lengthout)
       WRITE(unit=unite,rec=1,IOSTAT=istat) mold

    ELSEIF(kformat==2) THEN ! ASCII FILES

       WRITE(chardim,'(i3)') SIZE(mold,2)
       nco = SCAN(chardim,"123456789")

       CALL open_file_out_txt(kfile,lengthout)
       WRITE(unit=unite,fmt='('//chardim(nco:)//'e30.20)',IOSTAT=istat) &
            mold(1:SIZE(mold,1),1:SIZE(mold,2))

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN WRITING 2D VARIABLE TO FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END SUBROUTINE writevar_2d

  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------
  ! -----------------------------------------------------------------------------------------

  SUBROUTINE writevar_3d(kfile,mold,charname,kformat)
    ! -------------------------------------------------------------
    ! PURPOSE: EXPORT 3D VARIABLES TO BINARY FILES
    ! -------------------------------------------------------------
    ! INPUT:  
    !   - KFILE  = NUMBER ASSOCIATED TO THE NAME OF THE OUTPUT FILE
    !              (= F{KFILE}.BIN)
    !   - mold   = VARIABLE TO FILL IN
    ! RETURN VALUE:
    !   - INPUT ARRAY
    ! -------------------------------------------------------------
    INTEGER,INTENT(IN) :: kfile,kformat
    REAL(kind=DBL), DIMENSION(:,:,:), INTENT(IN) :: mold
    INTEGER:: i
    CHARACTER(len=14) charname
    INTEGER:: lengthout,istat=0

    ! FILL AN OUTPUT ASCII FILE WITH THE NAME AND THE DIMENSION OF THE VARIABLES
    WRITE(3311,fmt='(a18,8i10)') charname,kfile,SHAPE(mold)

    IF(kformat==1) THEN ! BINARY FILES

       INQUIRE(iolength=lengthout) mold
       CALL open_file_out_bin(kfile,lengthout)
       WRITE(unit=unite,rec=1,IOSTAT=istat) mold

    ELSEIF(kformat==2) THEN ! ASCII FILES

       CALL open_file_out_txt(kfile,lengthout)
       WRITE(unit=unite,fmt='(e30.20)',IOSTAT=istat) &
            mold(1:SIZE(mold,1),1:SIZE(mold,2),1:SIZE(mold,3))

    ELSE
       WRITE(stderr,*) 'BAD FORMAT SPECIFIED (BINARY OR ASCII)'
       WRITE(stderr,*) '=> PROGRAM STOPPED.'
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF

    CLOSE(unite)

    IF ( istat /= 0) THEN
       WRITE(stderr,100) filename, istat
100    FORMAT('PROBLEM IN WRITING 3D VARIABLE TO FILE ',A,'. IOSTAT = ',I6)
       CALL mpi_abort(mpi_comm_world,-1)
    ENDIF
  END SUBROUTINE writevar_3d

END MODULE mod_io_management
