!  $Author: pkubota 						$
!  $Date: 2008/09/23 17:51:54 					$
!  $Revision: 1.9 						$
!  $Revisions are currently made by the class's students.	$
!  $Update Date: 01/11/2023 10:19 AM				$
!  
!  Implementações: 
!  	1) Colocar INIT e FINALIZE - OK
! 	7) Amortecimento das condições de contorno - OK
!	8) Colocar todos os campos em "WRITE FIELDS" e salvar um bin - OK

MODULE Class_Module_Fields
 USE Constants, Only: r8, r4, i4, pi, Deg2Rad, r_earth, omega, nfprt

  IMPLICIT NONE
  PRIVATE       

  ! Selecting Unit

  INTEGER,PUBLIC, PARAMETER :: unitsurp=50
  INTEGER,PUBLIC, PARAMETER :: unituvel=51
  INTEGER,PUBLIC, PARAMETER :: unitvvel=52
  INTEGER,PUBLIC, PARAMETER :: unitomeg=56
  INTEGER,PUBLIC, PARAMETER :: unittemp=53
  INTEGER,PUBLIC, PARAMETER :: unitumes=54
  INTEGER,PUBLIC, PARAMETER :: unitzgeo=55
  INTEGER,PUBLIC, PARAMETER :: unitoutp=60
  
  INTEGER,PUBLIC, PARAMETER :: unitumodel=61  
  INTEGER,PUBLIC, PARAMETER :: unittmodel=62
  INTEGER,PUBLIC, PARAMETER :: unitmmodel=63
  INTEGER,PUBLIC, PARAMETER :: unitqmodel=64    

  INTEGER                    :: irec_local
  INTEGER,PUBLIC             :: nLon     = 161
  INTEGER,PUBLIC             :: nLat     = 161
  REAL(KIND=r4)              :: InitLon  = 295.0 !0  - 360 
  REAL(KIND=r4)              :: InitLat  = -50.0 !-90   90
  REAL(KIND=r8)              :: DeltaLon =  0.25
  REAL(KIND=r8)              :: DeltaLat =  0.25
  INTEGER,PUBLIC             :: nLev=37
  INTEGER                    :: lrec2D
  INTEGER                    :: lrec3D

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Plevs(:)
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: CoordLat(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: CoordLon(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: FcorPar (:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: DeltaLamda(:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: DeltaTheta(:,:) 

  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: u_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: v_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: w_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: t_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: q_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: z_ref(:,:,:) 
  REAL(KIND=r8),PUBLIC,ALLOCATABLE    :: p_ref(:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var2P_A(:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var2P_B(:,:) 
  
  REAL(KIND=r4), ALLOCATABLE :: var3U_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3U_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3V_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3V_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3W_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3W_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3T_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3T_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3Q_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3Q_B(:,:,:) 

  REAL(KIND=r4), ALLOCATABLE :: var3Z_A(:,:,:) 
  REAL(KIND=r4), ALLOCATABLE :: var3Z_B(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: U_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: U_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: V_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: V_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: T_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: T_C(:,:,:) 

  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Q_N(:,:,:) 
  REAL(KIND=r8),PUBLIC, ALLOCATABLE :: Q_C(:,:,:) 
  
!!!!!!!!!!!!!!!!!SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PUBLIC :: Init_Class_Module_Fields, ReadFields, WriteFields, Finalize_Class_Module_Fields
  
CONTAINS

 SUBROUTINE Init_Class_Module_Fields() !INIT
  IMPLICIT NONE
  INTEGER :: i,j
  ALLOCATE(Plevs(nLev));Plevs=0.0
  ALLOCATE(CoordLat(nLon,nLat));CoordLat=0.0
  ALLOCATE(CoordLon(nLon,nLat));CoordLon=0.0
  ALLOCATE(FcorPar(nLon,nLat));FcorPar=0.0
  ALLOCATE(DeltaLamda(nLon,nLat));DeltaLamda=0.0
  ALLOCATE(DeltaTheta(nLon,nLat));DeltaTheta=0.0
  
  ALLOCATE(var2P_A(nLon,nLat));var2P_A=0.0
  ALLOCATE(var2P_B(nLon,nLat));var2P_B=0.0

  ALLOCATE(var3U_A(nLon,nLat,nLev));var3U_A=0.0
  ALLOCATE(var3U_B(nLon,nLat,nLev));var3U_B=0.0

  ALLOCATE(var3V_A(nLon,nLat,nLev));var3V_A=0.0
  ALLOCATE(var3V_B(nLon,nLat,nLev));var3V_B=0.0

  ALLOCATE(var3W_A(nLon,nLat,nLev));var3W_A=0.0
  ALLOCATE(var3W_B(nLon,nLat,nLev));var3W_B=0.0

  ALLOCATE(var3T_A(nLon,nLat,nLev));var3T_A=0.0
  ALLOCATE(var3T_B(nLon,nLat,nLev));var3T_B=0.0

  ALLOCATE(var3Q_A(nLon,nLat,nLev));var3Q_A=0.0
  ALLOCATE(var3Q_B(nLon,nLat,nLev));var3Q_B=0.0

  ALLOCATE(var3Z_A(nLon,nLat,nLev));var3Z_A=0.0
  ALLOCATE(var3Z_B(nLon,nLat,nLev));var3Z_B=0.0

  ALLOCATE(U_N(nLon,nLat,nLev));U_N=0.0
  ALLOCATE(U_C(nLon,nLat,nLev));U_C=0.0

  ALLOCATE(V_N(nLon,nLat,nLev));V_N=0.0
  ALLOCATE(V_C(nLon,nLat,nLev));V_C=0.0

  ALLOCATE(T_N(nLon,nLat,nLev));T_N=0.0
  ALLOCATE(T_C(nLon,nLat,nLev));T_C=0.0

  ALLOCATE(Q_N(nLon,nLat,nLev));Q_N=0.0
  ALLOCATE(Q_C(nLon,nLat,nLev));Q_C=0.0
  
  ALLOCATE(p_ref(nLon,nLat));p_ref=0.0  
  ALLOCATE(u_ref(nLon,nLat,nLev));u_ref=0.0
  ALLOCATE(v_ref(nLon,nLat,nLev));v_ref=0.0
  ALLOCATE(w_ref(nLon,nLat,nLev));w_ref=0.0
  ALLOCATE(t_ref(nLon,nLat,nLev));t_ref=0.0
  ALLOCATE(q_ref(nLon,nLat,nLev));q_ref=0.0
  ALLOCATE(z_ref(nLon,nLat,nLev));z_ref=0.0

  INQUIRE(IOLENGTH=lrec2D)var2P_A
  INQUIRE(IOLENGTH=lrec3D)var3U_A

  OPEN(unit=unitzgeo,FILE='GeoPotential.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unittemp,FILE='Temperature.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitumes,FILE='SpecificHumidy.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unituvel,FILE='ZonalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitvvel,FILE='MeridionalWind.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitomeg,FILE='Omega.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='READ',STATUS='OLD') 

  OPEN(unit=unitsurp,FILE='SurfacePressure.bin',&
       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec2D,ACTION='READ',STATUS='OLD') 
       
!SAIDAS DO MODELO SAO SALVAS EM:
  OPEN(unit=unitumodel,FILE='ModelU.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN')
      
  OPEN(unit=unittmodel,FILE='ModelT.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN') 
      
  OPEN(unit=unitmmodel,FILE='ModelM.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN') 
      
  OPEN(unit=unitqmodel,FILE='ModelQ.bin',&
      FORM='UNFORMATTED',ACCESS='DIRECT',RECL=lrec3D,ACTION='WRITE',STATUS='UNKNOWN')   

!!!!!!!!!!!!!!!!!!!!DEFINIÇÃO DA MALHA (QUADRADA)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO j=1,nLat
     CoordLon(1,j) = InitLon*Deg2Rad
     DO i=2,nLon
        CoordLon(i,j) = CoordLon(i-1,j) + (DeltaLon*Deg2Rad)
     END DO
  END DO
  
  DO i=1,nLon
     CoordLat(i,1) = InitLat*Deg2Rad
     DO j=2,nLat
        CoordLat(i,j) = CoordLat(i,j-1) + (DeltaLat*Deg2Rad)
     END DO
  END DO

  DO j=1,nLat
     DO i=1,nLon
        FcorPar(i,j) = 2.0_r8*omega*sin(CoordLat(i,j))
     END DO
  END DO

  DO j=1,nLat
     DO i=2,nLon-1
        DeltaLamda(i,j) = (CoordLon(i+1,j)-CoordLon(i,j))
     END DO
        DeltaLamda(1,j)    = (DeltaLamda(2,j))
        DeltaLamda(nLon,j) = (DeltaLamda(nLon-1,j))
  END DO

  DO i=1,nLon
     DO j=2,nLat-1
        DeltaTheta(i,j) = (CoordLat(i,j+1)-CoordLat(i,j))
     END DO
        DeltaTheta(i,1) = DeltaTheta(i,2)
        DeltaTheta(i,nLat) = DeltaTheta(i,nLat-1)
  END DO
  Plevs(1:nLev)=(/1000.0_r8,975.0_r8,950.0_r8,925.0_r8,900.0_r8,875.0_r8,850.0_r8,825.0_r8,800.0_r8,775.0_r8,&
                   750.0_r8,700.0_r8,650.0_r8,600.0_r8,550.0_r8,500.0_r8,450.0_r8,400.0_r8,350.0_r8,300.0_r8,&
                   250.0_r8,225.0_r8,200.0_r8,175.0_r8,150.0_r8,125.0_r8,100.0_r8, 70.0_r8, 50.0_r8, 30.0_r8,&
                    20.0_r8, 10.0_r8,  7.0_r8,  5.0_r8,  3.0_r8,  2.0_r8,  1.0_r8/)             

  Plevs=Plevs*100.0_r8
 END SUBROUTINE Init_Class_Module_Fields

 SUBROUTINE ReadFields(irec,TimeIncrSeg) 
  IMPLICIT NONE 
  INTEGER      , INTENT(IN   ) :: irec
  REAL(KIND=r8), INTENT(IN   ) :: TimeIncrSeg
  REAL(KIND=r8) :: w1
  REAL(KIND=r8) :: w2
  INTEGER       :: i,j,k

  w1=1.0_r8-mod(TimeIncrSeg,3600.0_r8)/3600.0_r8
  w2=mod(TimeIncrSeg,3600.0_r8)/3600.0_r8

  IF(irec > irec_local)THEN
    READ(unituvel,rec=irec)var3U_A
    READ(unitvvel,rec=irec)var3V_A
    READ(unitomeg,rec=irec)var3W_A
    READ(unittemp,rec=irec)var3T_A
    READ(unitumes,rec=irec)var3Q_A
    READ(unitzgeo,rec=irec)var3Z_A
    READ(unitsurp,rec=irec)var2P_A

    READ(unituvel,rec=irec+1)var3U_B
    READ(unitvvel,rec=irec+1)var3V_B
    READ(unitomeg,rec=irec+1)var3W_B
    READ(unittemp,rec=irec+1)var3T_B
    READ(unitumes,rec=irec+1)var3Q_B
    READ(unitzgeo,rec=irec+1)var3Z_B
    READ(unitsurp,rec=irec+1)var2P_B
    irec_local=irec
  END IF

  p_ref= var2P_A *w1 + w2*var2P_B
  u_ref= var3U_A *w1 + w2*var3U_B
  v_ref= var3V_A *w1 + w2*var3V_B
  w_ref= var3W_A *w1 + w2*var3W_B
  t_ref= var3T_A *w1 + w2*var3T_B
  q_ref= var3Q_A *w1 + w2*var3Q_B
  z_ref= var3Z_A *w1 + w2*var3Z_B

  DO k=1,nLev

     z_ref(1   ,1:nLat,k)=0.25_r8*(z_ref(1   ,1:nLat,k) + z_ref(2     ,1:nLat,k) + z_ref(3     ,1:nLat,k) + z_ref(2     ,1:nLat,k) )
     z_ref(nLon,1:nLat,k)=0.25_r8*(z_ref(nLon,1:nLat,k) + z_ref(nLon-1,1:nLat,k) + z_ref(nLon-2,1:nLat,k) + z_ref(nLon-1,1:nLat,k) )

     z_ref(1:nLon,1   ,k)=0.25_r8*(z_ref(1:nLon,1   ,k) + z_ref(1:nLon,2     ,k) + z_ref(1:nLon,3     ,k) + z_ref(1:nLon,2     ,k) )
     z_ref(1:nLon,nLat,k)=0.25_r8*(z_ref(1:nLon,nLat,k) + z_ref(1:nLon,nLat-1,k) + z_ref(1:nLon,nLat-2,k) + z_ref(1:nLon,nLat-1,k) )

     w_ref(1   ,1:nLat,k)=0.25_r8*(w_ref(1   ,1:nLat,k) + w_ref(2     ,1:nLat,k) + w_ref(3     ,1:nLat,k) + w_ref(2     ,1:nLat,k) )
     w_ref(nLon,1:nLat,k)=0.25_r8*(w_ref(nLon,1:nLat,k) + w_ref(nLon-1,1:nLat,k) + w_ref(nLon-2,1:nLat,k) + w_ref(nLon-1,1:nLat,k) )

     w_ref(1:nLon,1   ,k)=0.25_r8*(w_ref(1:nLon,1   ,k) + w_ref(1:nLon,2     ,k) + w_ref(1:nLon,3     ,k) + w_ref(1:nLon,2     ,k) )
     w_ref(1:nLon,nLat,k)=0.25_r8*(w_ref(1:nLon,nLat,k) + w_ref(1:nLon,nLat-1,k) + w_ref(1:nLon,nLat-2,k) + w_ref(1:nLon,nLat-1,k) )
  END DO
  DO k=1,nLev
     DO j=2,nLat-1
        DO i=2,nLon-1
           z_ref(i   ,j,k)=(1.0_r8/5.0_r8)*(z_ref(i+1,j,k) + z_ref(i,j,k) + z_ref(i-1,j,k)+ z_ref(i,j+1,k)+ z_ref(i,j-1,k))
           w_ref(i   ,j,k)=(1.0_r8/5.0_r8)*(w_ref(i+1,j,k) + w_ref(i,j,k) + w_ref(i-1,j,k)+ w_ref(i,j+1,k)+ w_ref(i,j-1,k))
        END DO
     END DO
  END DO

  DO k=1,nLev-1
     u_ref(:,:,k)= 0.5_r8*(u_ref(:,:,k+1)+u_ref(:,:,k))
     v_ref(:,:,k)= 0.5_r8*(v_ref(:,:,k+1)+v_ref(:,:,k))
     w_ref(:,:,k)= 0.5_r8*(w_ref(:,:,k+1)+w_ref(:,:,k))
     t_ref(:,:,k)= 0.5_r8*(t_ref(:,:,k+1)+t_ref(:,:,k))
     q_ref(:,:,k)= 0.5_r8*(q_ref(:,:,k+1)+q_ref(:,:,k))
     z_ref(:,:,k)= 0.5_r8*(z_ref(:,:,k+1)+z_ref(:,:,k))
  END DO
 END SUBROUTINE ReadFields
 
 ! NAO ENTENDI A LOGICA MAS VOU SEGUI-LA
 
 SUBROUTINE WriteFields(irec)
  IMPLICIT NONE 
  INTEGER      , INTENT(INOUT) :: irec
  
  irec=irec+1
  

  WRITE(unittmodel,rec=irec)REAL(T_C,kind=r4) !Temperature
  WRITE(unitqmodel,rec=irec)REAL(Q_C,kind=r4) !SpecificHumidy
  WRITE(unitumodel,rec=irec)REAL(U_C,kind=r4) !ZonalWind
  WRITE(unitmmodel,rec=irec)REAL(V_C,kind=r4) !MeridionalWind

    
 END SUBROUTINE WriteFields

 SUBROUTINE Finalize_Class_Module_Fields() !FINALIZE
  IMPLICIT NONE 

  DEALLOCATE(Plevs)
  DEALLOCATE(CoordLat)
  DEALLOCATE(CoordLon)
  DEALLOCATE(FcorPar)
  DEALLOCATE(DeltaLamda)
  DEALLOCATE(DeltaTheta)
  DEALLOCATE(var2P_A)
  DEALLOCATE(var2P_B)
  DEALLOCATE(var3U_A)
  DEALLOCATE(var3U_B)
  DEALLOCATE(var3V_A)
  DEALLOCATE(var3V_B)
  DEALLOCATE(var3W_A)
  DEALLOCATE(var3W_B)
  DEALLOCATE(var3T_A)
  DEALLOCATE(var3T_B)
  DEALLOCATE(var3Q_A)
  DEALLOCATE(var3Q_B)
  DEALLOCATE(var3Z_A)
  DEALLOCATE(var3Z_B)
  DEALLOCATE(U_N)   !tempo futuro
  DEALLOCATE(U_C)
  DEALLOCATE(V_N)
  DEALLOCATE(V_C)
  DEALLOCATE(T_N)
  DEALLOCATE(T_C)
  DEALLOCATE(Q_N)
  DEALLOCATE(Q_C)
  DEALLOCATE(p_ref)  
  DEALLOCATE(u_ref)
  DEALLOCATE(v_ref)
  DEALLOCATE(w_ref)
  DEALLOCATE(t_ref)
  DEALLOCATE(q_ref)
  DEALLOCATE(z_ref)

  CLOSE(unit=unitzgeo, STATUS='KEEP')
  CLOSE(unit=unittemp, STATUS='KEEP')
  CLOSE(unit=unitumes, STATUS='KEEP')
  CLOSE(unit=unituvel, STATUS='KEEP')
  CLOSE(unit=unitvvel, STATUS='KEEP')
  CLOSE(unit=unitomeg, STATUS='KEEP')
  CLOSE(unit=unitsurp, STATUS='KEEP')
  CLOSE(unit=unitoutp, STATUS='KEEP')
  
  CLOSE(unit=unitumodel, STATUS='KEEP')
  CLOSE(unit=unittmodel, STATUS='KEEP')
  CLOSE(unit=unitmmodel, STATUS='KEEP')
  CLOSE(unit=unitqmodel, STATUS='KEEP')
 


 END SUBROUTINE Finalize_Class_Module_Fields

END MODULE Class_Module_Fields
