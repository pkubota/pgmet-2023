MODULE NumDerivate
 IMPLICIT NONE
 PRIVATE
 INTEGER, PUBLIC      , PARAMETER  :: r8=8
 INTEGER, PUBLIC      , PARAMETER  :: r4=4

 PUBLIC :: Fx
 PUBLIC :: dFxdx
 PUBLIC :: dif_Forward
 PUBLIC :: dif_Centrada2ord
 PUBLIC :: dif_Centrada2ord2
 PUBLIC :: dif_Centrada4ord2
 PUBLIC :: ErrorRelative
 CONTAINS

 FUNCTION Fx(x)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8)             :: func
   func=(x-10_r8)**2
 END FUNCTION Fx

 FUNCTION dFxdx(x)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8)             :: func
   func=2.0_r8*(x-10_r8)
 END FUNCTION dFxdx

 FUNCTION dif_Forward(x,dx)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8), INTENT(IN) :: dx
   REAL(KIND=r8)             :: func
   func=(Fx(x+dx)-Fx(x))/dx

 END FUNCTION dif_Forward

 FUNCTION dif_Centrada2ord(x,dx)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8), INTENT(IN) :: dx
   REAL(KIND=r8)             :: func
   func=(Fx(x+dx)-Fx(x-dx))/(2*dx)

 END FUNCTION dif_Centrada2ord

 FUNCTION dif_Centrada2ord2(x,dx)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8), INTENT(IN) :: dx
   REAL(KIND=r8)             :: func
   func=(-Fx(x+2.0_r8*dx)+4.0_r8*Fx(x+dx) -3.0_r8*Fx(x))/(2*dx)

 END FUNCTION dif_Centrada2ord2

 FUNCTION dif_Centrada4ord2(x,dx)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: x
   REAL(KIND=r8), INTENT(IN) :: dx
   REAL(KIND=r8)             :: func
   func=(-Fx(x+2.0_r8*dx)+8.0_r8*Fx(x+dx) -8.0_r8*Fx(x-dx) + Fx(x-2.0_r8*dx))/(12.0_r8*dx)

 END FUNCTION dif_Centrada4ord2

 FUNCTION ErrorRelative(funcAn,FuncNu)  RESULT (func)
    IMPLICIT NONE
   REAL(KIND=r8), INTENT(IN) :: funcAn
   REAL(KIND=r8), INTENT(IN) :: FuncNu
   REAL(KIND=r8)             :: func
   func=(funcAn-FuncNu)/funcAn
 END FUNCTION ErrorRelative 
END MODULE NumDerivate

PROGRAM Main
 USE NumDerivate, Only: Fx,dFxdx,dif_Forward,dif_Centrada2ord,&
                        dif_Centrada2ord2,dif_Centrada4ord2,ErrorRelative,r8
 IMPLICIT NONE
 REAL(KIND=r8) :: xpts=5.5_r8
 REAL(KIND=r8) :: dx  =0.2_r8

  CALL Init()
  CALL Run()
  CALL Finalize()

  CONTAINS

 SUBROUTINE Init()
 IMPLICIT NONE
 PRINT*,'x=',xpts,'dx=',dx ,'FuncAn=',Fx(xpts)
 PRINT*,'Derivada Analitica=', dFxdx(xpts),'(a) Derivada numerica=',dif_Forward(xpts,dx),&
        'ErrorRelative=',ErrorRelative(dFxdx(xpts),dif_Forward(xpts,dx))

 PRINT*,'Derivada Analitica=', dFxdx(xpts),'(b) Derivada numerica=',dif_Centrada2ord(xpts,dx),&
        'ErrorRelative=',ErrorRelative(dFxdx(xpts),dif_Centrada2ord(xpts,dx))

 PRINT*,'Derivada Analitica=', dFxdx(xpts),'(c) Derivada numerica=',dif_Centrada2ord2(xpts,dx),&
        'ErrorRelative=',ErrorRelative(dFxdx(xpts),dif_Centrada2ord2(xpts,dx))

 PRINT*,'Derivada Analitica=', dFxdx(xpts),'(d) Derivada numerica=',dif_Centrada4ord2(xpts,dx),&
        'ErrorRelative=',ErrorRelative(dFxdx(xpts),dif_Centrada4ord2(xpts,dx))



 END SUBROUTINE Init
 
 SUBROUTINE Run()
 IMPLICIT NONE

 END SUBROUTINE Run

 SUBROUTINE Finalize()
 IMPLICIT NONE

 END SUBROUTINE Finalize

END PROGRAM Main
