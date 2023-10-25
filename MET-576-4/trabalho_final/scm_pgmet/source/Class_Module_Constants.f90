!
!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!
MODULE Constants

  IMPLICIT NONE
  PRIVATE       

  ! Selecting Kinds
  INTEGER,PUBLIC, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers
  ! Selecting Unit

  INTEGER,PUBLIC, PARAMETER :: nfprt=0
  PUBLIC :: InitClassModuleConstants
CONTAINS

 SUBROUTINE InitClassModuleConstants () 
    IMPLICIT NONE

 END SUBROUTINE InitClassModuleConstants

END MODULE Constants
