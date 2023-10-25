!
!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!
MODULE Class_Module_Fields
 USE Constants, Only: r8,r4,i4,nfprt

  IMPLICIT NONE
  PRIVATE       

  ! Selecting Unit

  INTEGER,PUBLIC, PARAMETER :: unitsurp=50
  INTEGER,PUBLIC, PARAMETER :: unituvel=51
  INTEGER,PUBLIC, PARAMETER :: unitvvel=52
  INTEGER,PUBLIC, PARAMETER :: unittemp=53
  INTEGER,PUBLIC, PARAMETER :: unitumes=54
  INTEGER,PUBLIC, PARAMETER :: unitzgeo=55
  INTEGER,PUBLIC, PARAMETER :: unitomeg=56

  INTEGER                    :: nLon=401
  INTEGER                    :: nLat=321
  INTEGER                    :: nLev=37
  INTEGER                    :: lrec2D
  INTEGER                    :: lrec3D
  REAL(KIND=r4), ALLOCATABLE :: var2D(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3D(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: uvelc(:,:,:) 
  PUBLIC :: Init_Class_Module_Fields,ReadFields
CONTAINS

 SUBROUTINE Init_Class_Module_Fields () 
  IMPLICIT NONE
  ALLOCATE(var2D(nLon,nLat));var2D=0.0
  ALLOCATE(var3D(nLon,nLat,nLev));var3D=0.0
  INQUIRE(IOLENGTH=lrec2D)var2D
  INQUIRE(IOLENGTH=lrec3D)var3D

  OPEN(unit=unitzgeo,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/GeoPotential.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unittemp,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/Temperature.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitumes,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/SpecificHumidy.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unituvel,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/ZonalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitvvel,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/MeridionalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitomeg,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/Omega.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitsurp,FILE='/cygdrive/d/paulo.kubota/pgmet/scm_pgmet/datain/SurfacePressure.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec2D,ACTION='READ',STATUS='OLD') 
   
 END SUBROUTINE Init_Class_Module_Fields

 SUBROUTINE ReadFields (irec) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN   ) :: irec
  REAL(KIND=r4) :: w1
  REAL(KIND=r4) :: w2

  READ(unituvel,rec=irec)var3D
  READ(unitvvel,rec=irec)var3D
  READ(unittemp,rec=irec)var3D
  READ(unitumes,rec=irec)var3D
  READ(unitzgeo,rec=irec)var3D
  READ(unitsurp,rec=irec)var2D

  READ(unituvel,rec=irec+1)var3D
  READ(unitvvel,rec=irec+1)var3D
  READ(unittemp,rec=irec+1)var3D
  READ(unitumes,rec=irec+1)var3D
  READ(unitzgeo,rec=irec+1)var3D
  READ(unitsurp,rec=irec+1)var2D
  
  uvelc= var3D *w1 + var3D*w2

 END SUBROUTINE ReadFields
END MODULE Class_Module_Fields
