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
  INTEGER                    :: irec_local
  INTEGER                    :: nLon=161
  INTEGER                    :: nLat=161
  INTEGER                    :: nLev=37
  INTEGER                    :: lrec2D
  INTEGER                    :: lrec3D
  REAL(KIND=r8), ALLOCATABLE :: uvelc(:,:,:) 
  REAL(KIND=r8), ALLOCATABLE :: pslc(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var2D_A(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var2D_B(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3D_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3D_B(:,:,:) 

  PUBLIC :: Init_Class_Module_Fields,ReadFields
CONTAINS

 SUBROUTINE Init_Class_Module_Fields () 
  IMPLICIT NONE
  ALLOCATE(var2D_A(nLon,nLat));var2D_A=0.0
  ALLOCATE(var2D_B(nLon,nLat));var2D_B=0.0
  ALLOCATE(pslc(nLon,nLat));pslc=0.0

  ALLOCATE(var3D_A(nLon,nLat,nLev));var3D_A=0.0
  ALLOCATE(var3D_B(nLon,nLat,nLev));var3D_B=0.0
  ALLOCATE(uvelc(nLon,nLat,nLev));uvelc=0.0
  INQUIRE(IOLENGTH=lrec2D)var2D_A
  INQUIRE(IOLENGTH=lrec3D)var3D_A

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

 SUBROUTINE ReadFields (irec,TimeIncrSeg) 
  IMPLICIT NONE 
  INTEGER      , INTENT(IN   ) :: irec
  REAL(KIND=r8), INTENT(IN   ) :: TimeIncrSeg
  REAL(KIND=r8) :: w1
  REAL(KIND=r8) :: w2

  IF(irec > irec_local)THEN
    READ(unituvel,rec=irec)var3D_A
    READ(unitvvel,rec=irec)var3D_A
    READ(unittemp,rec=irec)var3D_A
    READ(unitumes,rec=irec)var3D_A
    READ(unitzgeo,rec=irec)var3D_A
    READ(unitsurp,rec=irec)var2D_A

    READ(unituvel,rec=irec+1)var3D_B
    READ(unitvvel,rec=irec+1)var3D_B
    READ(unittemp,rec=irec+1)var3D_B
    READ(unitumes,rec=irec+1)var3D_B
    READ(unitzgeo,rec=irec+1)var3D_B
    READ(unitsurp,rec=irec+1)var2D_B
    irec_local=irec
  ELSE


  END IF
  w1=1.0_r8-mod(TimeIncrSeg,3600.0_r8)/3600.0_r8
  w2=mod(TimeIncrSeg,3600.0_r8)/3600.0_r8

  pslc= var2D_A *w1 + w2*var2D_B

 END SUBROUTINE ReadFields
END MODULE Class_Module_Fields
