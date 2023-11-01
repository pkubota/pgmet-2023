!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!  $Revisions are currently made by the class's students.
!  $Update Date: 31/10/2023
!  
!  Implementações: 
!  	1) Colocar INIT e FINALIZE em todas as rotinas e no Main - OK
! 	2) Colocar as equações dos campos no "Class_Module_Dynamics"
! 	3) No Main, implementar para salvar os campos
! 	4) Como fazer o Call Physics?
!	5) Testar com os dados do ERA5 a inicialização e saving

PROGRAM Main
 USE Constants, Only: InitClassModuleConstants, r8,r4
 USE Class_Module_TimeManager, Only : Init_Class_Module_TimeManager, dt_step,&
                                      SetTimeControl,idatei,idatec,idatef,ICtrDay,&
                                      TimeIncrementSeg,GetRec2ReadWrite
 USE Class_Module_Fields, Only : Init_Class_Module_Fields,FinalizeFields,ReadFields,WriteFields,nLat,nLon,nLev
 USE Class_Module_Dynamics, Only : Init_Class_Module_Dynamics,RunDynamics
 IMPLICIT NONE
 
 CALL Init()
 CALL Run()
 CALL Finalize()

CONTAINS 
 SUBROUTINE Init()
  IMPLICIT NONE
  CALL InitClassModuleConstants()
  CALL Init_Class_Module_TimeManager()
  CALL Init_Class_Module_Fields()
  CALL Init_Class_Module_Dynamics(nLat,nLon,nLev,dt_step)
 END SUBROUTINE Init

 SUBROUTINE Run()
  IMPLICIT NONE
  INTEGER       :: nMaxIteration
  INTEGER       :: itr
  INTEGER       :: jhr
  INTEGER       :: jmon
  INTEGER       :: jday
  INTEGER       :: jyr
  INTEGER       :: test
  INTEGER       :: rec,irecw
  INTEGER       :: ktm,kt,ktp
  REAL(KIND=r8) :: TimeIncrSeg
  REAL(KIND=r8) :: ahour,bhour

  ahour=0.0_r8;bhour=0.0_r8;irecw=0
  TimeIncrSeg=0.0_r8
  nMaxIteration=SetTimeControl(idatei,idatef)
  DO itr=1,nMaxIteration-1
      !
      !     step loop starts
      !
      rec=GetRec2ReadWrite(idatei,idatec)
      CALL ReadFields(rec,TimeIncrSeg)
      !CALL Physics
      CALL RunDynamics(itr)
      
      IF(MOD(TimeIncrSeg,3600.0_r8) == 0.0_r8)THEN
         PRINT *,'itr=',itr,'tod=',TimeIncrSeg,MOD(TimeIncrSeg,3600.0_r8),'idatec=',idatec
      END IF

      TimeIncrSeg=TimeIncrSeg+dt_step
      IF(ABS( MOD(TimeIncrSeg+0.03125_r8,86400.0_r8)-0.03125_r8).LT.0.0625_r8)THEN
        TimeIncrSeg=0.0_r8
        ICtrDay=ICtrDay+1
      END IF
      test=TimeIncrementSeg(idatei,idatec,ICtrDay,TimeIncrSeg,jhr,jday,jmon,jyr,&
                            ktm,kt,ktp,ahour,bhour,nMaxIteration,dt_step,itr)
       IF(MOD(TimeIncrSeg,3600.0_r8) == 0.0_r8)THEN
         CALL WriteFields(irecw)
       END IF
  END DO
 END SUBROUTINE Run

 SUBROUTINE Finalize()
  IMPLICIT NONE
   CALL Finalize_Class_Module_Constants()
   CALL Finalize_Class_Module_TimeManager()
   CALL Finalize_Class_Module_Fields()
   CALL Finalize_Class_Module_Dynamics(nLat,nLon,nLev,dt_step)
 END SUBROUTINE Finalize

END PROGRAM Main
