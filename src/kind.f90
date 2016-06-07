!--------------------------------------------------------------------
!PURPOSE: DEFINE KIND AND STDERR/STDOUT NUMBERS TO USE THROUGHOUT CODE
!--------------------------------------------------------------------

MODULE kind
  ! Define kind numbers
  INTEGER, PARAMETER :: SGL = SELECTED_REAL_KIND(p=6,r=37)
  INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=12,r=50)

  ! Define units of standard output and standard error
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: stderr = 0

END MODULE kind
