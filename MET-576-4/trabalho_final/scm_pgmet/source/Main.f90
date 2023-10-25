PROGRAM Main
 USE Constants, Only: InitClassModuleConstants, r8,r4
 USE Class_Module_TimeManager, Only : Init_Class_Module_TimeManager, dt_step,&
                                      SetTimeControl,idatei,idatec,idatef,ICtrDay,&
                                      TimeIncrementSeg,GetRec2ReadWrite
 USE Class_Module_Fields, Only : Init_Class_Module_Fields,ReadFields
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
  INTEGER       :: rec
  INTEGER       :: ktm,kt,ktp
  REAL(KIND=r8) :: TimeIncrSeg
  REAL(KIND=r8) :: ahour

  
  TimeIncrSeg=0.0_r8
  nMaxIteration=SetTimeControl(idatei,idatef)
  DO itr=1,nMaxIteration
      !
      !     step loop starts
      !
      rec=GetRec2ReadWrite(idatei,idatec)
      CALL ReadFields(rec)

      TimeIncrSeg=TimeIncrSeg+dt_step
      IF(ABS( MOD(TimeIncrSeg+0.03125_r8,86400.0_r8)-0.03125_r8).LT.0.0625_r8)THEN
        TimeIncrSeg=0.0_r8
        ICtrDay=ICtrDay+1
      END IF
      test=TimeIncrementSeg(idatei,idatec,ICtrDay,TimeIncrSeg,jhr,jday,jmon,jyr,&
                            ktm,kt,ktp,ahour,nMaxIteration,dt_step,itr)

      PRINT *,'itr=',itr,'rec=',rec,'idatei=',idatei,'idatec=',idatec
  END DO
 END SUBROUTINE Run

 SUBROUTINE Finalize()
  IMPLICIT NONE

 END SUBROUTINE Finalize






END PROGRAM Main
