!
!  $Author: pkubota $
!  $Date: 2008/09/23 17:51:54 $
!  $Revision: 1.9 $
!
MODULE Class_Module_Dynamics
 USE Constants, Only: r8,r4,i4,pi,Deg2Rad,Rd,Cp,kappa,r_earth,omega,CTv,nfprt
  
 
 USE Class_Module_Fields, Only: u_ref, v_ref,w_ref, t_ref,q_ref,z_ref,p_ref,&
                                U_N,U_C,  V_N,V_C,  T_N,T_C,  Q_N,Q_C,&
                                Plevs,CoordLat,CoordLon,FcorPar ,DeltaLamda,DeltaTheta
  IMPLICIT NONE
  PRIVATE       

  INTEGER ::  Idim
  INTEGER ::  Jdim
  INTEGER ::  Kdim
  REAL(KIND=8) :: DeltaT

  ! Selecting Unit



  PUBLIC :: Init_Class_Module_Dynamics
  PUBLIC :: RunDynamics
CONTAINS

 SUBROUTINE Init_Class_Module_Dynamics (nLat,nLon,nLev,dt_step) 
  IMPLICIT NONE
  INTEGER , INTENT(IN   ) :: nLat
  INTEGER , INTENT(IN   ) :: nLon
  INTEGER , INTENT(IN   ) :: nLev
  REAL(KIND=8), INTENT(IN   )  :: dt_step
  Idim=nLon
  Jdim=nLat
  Kdim=nLev
  DeltaT=dt_step

 END SUBROUTINE Init_Class_Module_Dynamics

 SUBROUTINE RunDynamics (it) 
  IMPLICIT NONE
  INTEGER     , INTENT(IN   ) :: it
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku1, uo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: ku2, ku3, ku4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv1, vo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kv2, kv3, kv4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt1, to
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kt2, kt3, kt4
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq1, qo
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:Kdim) :: kq2, kq3, kq4

  REAL(KIND=8) :: dt6
  REAL(KIND=8) :: dt2
  REAL(KIND=8) :: dt
  INTEGER      :: i,j,k
  dt=DeltaT;  dt6 = DeltaT/6.0; dt2=DeltaT*0.5

  IF(it==1)THEN
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k) =  u_ref(i,j,k)
              vo(i,j,k) =  v_ref(i,j,k)
              to(i,j,k) =  t_ref(i,j,k)
              qo(i,j,k) =  q_ref(i,j,k)
           END DO
        END DO
     END DO
  ELSE
     DO k=1,Kdim
        DO j=1,Jdim
           DO i=1,  Idim
              uo(i,j,k)=  U_C(i,j,k)
              vo(i,j,k)=  V_C(i,j,k)
              to(i,j,k)=  T_C(i,j,k)
              qo(i,j,k)=  Q_C(i,j,k)
           END DO
        END DO
     END DO
  END IF
 !       PRINT*,'h eta',MAXVAL(h),MINVAL(h),'u eta',MAXVAL(u),MINVAL(u),'v eta',MAXVAL(v),MINVAL(v)
 
 call Solve_Forward_Beta_plane_Grid_A(ku1, kv1, kt1, kq1, uo        , vo             , to         , qo        )
 call Solve_Forward_Beta_plane_Grid_A(ku2, kv2, kt2, kq2, uo+ku1*dt2, vo+kv1*dt2     , to+kt1*dt2 , qo+kq1*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku3, kv3, kt3, kq3, uo+ku2*dt2, vo+kv2*dt2     , to+kt2*dt2 , qo+kq2*dt2)
 call Solve_Forward_Beta_plane_Grid_A(ku4, kv4, kt4, kq4, uo+ku3*dt , vo+kv3*dt      , to+kt3*dt  , qo+kq3*dt )
   

 ! final step and time marching / new values for RK4
 DO k=1,Kdim
    DO j=1,Jdim
       DO i=1,  Idim
          U_N(i,j,k) = uo(i,j,k) + (ku1(i,j,k) + 2.0*ku2(i,j,k) + 2.0*ku3(i,j,k) + ku4(i,j,k))*dt6
          V_N(i,j,k) = vo(i,j,k) + (kv1(i,j,k) + 2.0*kv2(i,j,k) + 2.0*kv3(i,j,k) + kv4(i,j,k))*dt6
          T_N(i,j,k) = to(i,j,k) + (kt1(i,j,k) + 2.0*kt2(i,j,k) + 2.0*kt3(i,j,k) + kt4(i,j,k))*dt6
          Q_N(i,j,k) = qo(i,j,k) + (kq1(i,j,k) + 2.0*kq2(i,j,k) + 2.0*kq3(i,j,k) + kq4(i,j,k))*dt6
       END DO
    END DO
 END DO

 ! updating the data
 DO k=1,Kdim
    DO j=1,Jdim
       DO i=1,  Idim
          U_C(i,j,k) = U_N(i,j,k) 
          V_C(i,j,k) = V_N(i,j,k)
          T_C(i,j,k) = T_N(i,j,k)
          Q_C(i,j,k) = Q_N(i,j,k)
       END DO
    END DO
 END DO
 END SUBROUTINE RunDynamics
 
 
 SUBROUTINE  Solve_Forward_Beta_plane_Grid_A(TermEqMomU, TermEqMomV,TermEqConT, TermEqConQ, u_in, v_in, t_in,q_in)
  IMPLICIT NONE
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(in) ::  u_in,  v_in,  t_in , q_in
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim), intent(out) :: TermEqMomU,TermEqMomV, TermEqConT, TermEqConQ
  REAL(KIND=8), dimension(1:Idim, 1:Jdim,1:kdim) ::  u,  v,  t , q
  REAL(KIND=8) :: udux
  REAL(KIND=8) :: vduy
  REAL(KIND=8) :: wduz
  REAL(KIND=8) :: fcov
  REAL(KIND=8) :: dPdx
  REAL(KIND=8) :: lnPs_xf
  REAL(KIND=8) :: lnPs_xb
  REAL(KIND=8) :: Tv
  INTEGER :: i,j,k
  INTEGER :: xb,xc,xf
  INTEGER :: yb,yc,yf
  u=u_in;v=v_in;t=t_in;q=q_in
  
  u(1:Idim   ,1:Jdim,Kdim)=(u_ref(1:Idim   ,1:Jdim,Kdim) + u(1:Idim,1:Jdim,Kdim-1))

  v(1:Idim   ,1:Jdim,Kdim)=(v_ref(1:Idim   ,1:Jdim,Kdim) + v(1:Idim,1:Jdim,Kdim-1))

  t(1:Idim   ,1:Jdim,Kdim)=(t_ref(1:Idim   ,1:Jdim,Kdim) + t(1:Idim,1:Jdim,Kdim-1))

  q(1:Idim   ,1:Jdim,Kdim)=(q_ref(1:Idim   ,1:Jdim,Kdim) + q(1:Idim,1:Jdim,Kdim-1))


  DO k=1,Kdim
     u(1   ,1:Jdim,k)=(u_ref(1   ,1:Jdim,k) + u(2     ,1:Jdim,k))
     u(Idim,1:Jdim,k)=(u_ref(Idim,1:Jdim,k) + u(Idim-1,1:Jdim,k))

     u(1:Idim,1   ,k)=(u_ref(1:Idim,1   ,k) + u(1:Idim,2     ,k))
     u(1:Idim,Jdim,k)=(u_ref(1:Idim,Jdim,k) + u(1:Idim,Jdim-1,k))

     v(1   ,1:Jdim,k)=(v_ref(1   ,1:Jdim,k) + v(2     ,1:Jdim,k))
     v(Idim,1:Jdim,k)=(v_ref(Idim,1:Jdim,k) + v(Idim-1,1:Jdim,k))

     v(1:Idim,1   ,k)=(v_ref(1:Idim,1   ,k) + v(1:Idim,2     ,k))
     v(1:Idim,Jdim,k)=(v_ref(1:Idim,Jdim,k) + v(1:Idim,Jdim-1,k))

     t(1   ,1:Jdim,k)=(t_ref(1   ,1:Jdim,k) + t(2     ,1:Jdim,k))
     t(Idim,1:Jdim,k)=(t_ref(Idim,1:Jdim,k) + t(Idim-1,1:Jdim,k))

     t(1:Idim,1   ,k)=(t_ref(1:Idim,1   ,k) + t(1:Idim,2     ,k))
     t(1:Idim,Jdim,k)=(t_ref(1:Idim,Jdim,k) + t(1:Idim,Jdim-1,k))

     q(1   ,1:Jdim,k)=(q_ref(1   ,1:Jdim,k) + q(2     ,1:Jdim,k))
     q(Idim,1:Jdim,k)=(q_ref(Idim,1:Jdim,k) + q(Idim-1,1:Jdim,k))

     q(1:Idim,1   ,k)=(q_ref(1:Idim,1   ,k) + q(1:Idim,2     ,k))
     q(1:Idim,Jdim,k)=(q_ref(1:Idim,Jdim,k) + q(1:Idim,Jdim-1,k))

  END DO
  
  
  DO k=1,Kdim-1 
      DO j=2,Jdim-1
         CALL index(j,Jdim,yb,yc,yf)
         DO i=2,Idim-1
            CALL index(i,Idim,xb,xc,xf)
            !
            !                       --               --   
            !                      |                   |  
            !            1         |          du       |  
            !udux= ----------------|   U * ----------  |  
            !       a*cos^2(theta) |         d lambda  |  
            !                      |                   |  
            !                       --               --   

            udux = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc))**2))) * (u(xc,yc,k) *((u(xf,yc,k) - u(xb,yc,k))/(DeltaLamda(xc,yc))))
            !
            !                       --                --
            !                       |                   |
            !             1         |          du       | 
            !vduy=  ----------------|   V * ----------  | 
            !        a*cos(theta)   |         d theta   | 
            !                       |                   |
            !                       --                 --


            vduy  = (1.0_r8/(r_earth*(cos(CoordLat(xc,yc)))))  *  (v(xc,yc,k)*((u(xc,yf,k) - u(xc,yb,k))/(DeltaTheta(xc,yc))))

            !      --                --
            !      |                   |
            !      |          du       | 
            !wduz= |   w * ----------  | 
            !      |         d P       | 
            !      |                   |
            !      --                 --

            wduz = w_ref (xc,yc,k) * ((u(xc,yf,k+1) - u(xc,yb,k))/(Plevs(k+1)-Plevs(k)))
 
            !      --       --
            !      |         |
            !      |         | 
            !fcov= | f  * v  | 
            !      |         | 
            !      |         |
            !      --      --
    
            fcov  = - FcorPar(xc,yc)*v(xc,yc,k)
            !
            !           --                            --
            !          |                                |
            !        1 |       dZGeo           dln(P)   | 
            !dPdx=  ---|    ---------- + Rd*Tv--------  | 
            !        a |     dLambda           dLambda  | 
            !          |                                |
            !          --                             --
            Tv=t(xc,yc,k)*(1.0_r8 + CTv*q(xc,yc,k))
            lnPs_xf=log(p_ref(xf,yc)*(Plevs(k)/p_ref(xf,yc)))
            lnPs_xb=log(p_ref(xb,yc)*(Plevs(k)/p_ref(xb,yc)))
            dPdx =(1.0_r8/(r_earth)) * (((z_ref(xf,yc,k) - z_ref(xb,yc,k))/(DeltaLamda(xc,yc))) + &
                                          Rd*Tv*((lnPs_xf-lnPs_xb)/(DeltaLamda(xc,yc))))
            !TermEqMomU(i,j,k) = -(udux +  vduy + betav + gdhx + ru + vis2dudx + vis2dudy)
         END DO
      END DO
  END DO


 END SUBROUTINE  Solve_Forward_Beta_plane_Grid_A

   SUBROUTINE index(i,Idim,xb,xc,xf)
      IMPLICIT NONE
      INTEGER, INTENT(IN   ) :: i
      INTEGER, INTENT(IN   ) :: Idim
      INTEGER, INTENT(OUT  ) :: xb,xc,xf
      IF(i==1) THEN
        xb=Idim
        xc=i
        xf=i+1
      ELSE IF(i==Idim)THEN
        xb=Idim-1
        xc=Idim
        xf=1
      ELSE
        xb=i-1
        xc=i
        xf=i+1
      END IF
   END SUBROUTINE index
END MODULE Class_Module_Dynamics
