! dc3d subroutines by Okada
! wrapper based on subroutine by Shinichi Miyazaki (with his written agreement)
! downloaded in December 2012 from http://www-geod.kugi.kyoto-u.ac.jp/~miyazaki/DC_2010/

module okada

  implicit none
  private

  public :: compute_kernel

contains

! Computes the slip-parallel shear stress (tau) and normal stress (sigma)
! at a point (OX,OY,OZ)
! on a fault with local dip angle O_DIP and strike parallel to X
! induced by unit slip (1 m) on a rectangular fault cell 
! of along-strike size L and along-dip size W 
! centered at (SX,SY,SZ) with dip angle S_DIP and strike parallel to X
! in a homogeneous elastic half-space with Lame moduli LAM and MU.
! The source and receiver faults are either strike-slip or thrust:
! mode = 1 right-lateral strike-slip 
!       -1 left-lateral strike-slip 
!        2 thrust dip-slip
!       -2 normal dip-slip
!
subroutine compute_kernel(LAM,MU,SX,SY,SZ,S_DIP,L,W,OX,OY,OZ,O_DIP,IRET,tau,sigma,mode)

  use constants, only : PI 

  double precision, intent(in) :: LAM,MU, SX,SY,SZ,S_DIP, L,W,OX,OY,OZ,O_DIP 
  integer, intent(inout) :: IRET
  double precision, intent(inout) :: tau, sigma
  integer, intent(in) :: mode

  double precision :: ALPHA, S_DEPTH, X, Y, Z, &
    UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ, &
    Ustrike, Udip, Unorm
  double precision :: STRESS(3,3), STRAIN(3,3),TR, n_f(3), tau_n(3), n_dir(3)

  ALPHA = (LAM+MU)/(LAM+2d0*MU)

  S_DEPTH = -1d0*SZ

  select case (mode)
    case (1) ! strike-slip (right-lateral)
      Ustrike = 1d0
      Udip = 0d0
      n_dir(1) = -1d0
      n_dir(2) = 0d0
      n_dir(3) = 0d0
    case (-1) ! strike-slip (left-lateral)
      Ustrike = -1d0
      Udip = 0d0
      n_dir(1) = 1d0
      n_dir(2) = 0d0
      n_dir(3) = 0d0
    case (2) ! dip-slip (thrust)
      Ustrike = 0d0
      Udip = 1d0
      n_dir(1) = 0d0
      n_dir(2) = -cos(O_DIP/180d0*PI)
      n_dir(3) = -sin(O_DIP/180d0*PI)
    case (-2) ! dip-slip (normal)
      Ustrike = 0d0
      Udip = -1d0
      n_dir(1) = 0d0
      n_dir(2) = -cos(O_DIP/180d0*PI)
      n_dir(3) = -sin(O_DIP/180d0*PI)
    case default
      stop 'FATAL ERROR in okada.compute_kernel : mode should be +/- 1 (strike-slip) or 2 (thrust)'
  end select
  Unorm = 0d0 ! no opening

  X = OX-SX
  Y = OY-SY
  Z = OZ
!Normal vector of receiver fault
  n_f(1) = 0d0
  n_f(2) = -sin(O_DIP/180d0*PI)
  n_f(3) = cos(O_DIP/180d0*PI)
 
  call DC3D(ALPHA,X,Y,Z,S_DEPTH,S_DIP,-0.5d0*L,0.5d0*L,-0.5d0*W,0.5d0*W, &
            Ustrike,Udip,Unorm,   &
            UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)

  STRAIN(1,1) = UXX
  STRAIN(1,2) = 0.5d0*(UXY+UYX)
  STRAIN(1,3) = 0.5d0*(UXZ+UZX)
  STRAIN(2,1) = STRAIN(1,2)
  STRAIN(2,2) = UYY
  STRAIN(2,3) = 0.5d0*(UYZ+UZY)
  STRAIN(3,1) = STRAIN(1,3)
  STRAIN(3,2) = STRAIN(2,3)
  STRAIN(3,3) = UZZ
  TR = STRAIN(1,1)+STRAIN(2,2)+STRAIN(3,3)

  STRESS = 2d0*MU*STRAIN
  STRESS(1,1) = STRESS(1,1)+LAM*TR
  STRESS(2,2) = STRESS(2,2)+LAM*TR
  STRESS(3,3) = STRESS(3,3)+LAM*TR
  
  tau_n(1) = STRESS(1,1)*n_f(1)+STRESS(1,2)*n_f(2)+STRESS(1,3)*n_f(3)
  tau_n(2) = STRESS(2,1)*n_f(1)+STRESS(2,2)*n_f(2)+STRESS(2,3)*n_f(3)
  tau_n(3) = STRESS(3,1)*n_f(1)+STRESS(3,2)*n_f(2)+STRESS(3,3)*n_f(3)

  tau = tau_n(1)*n_dir(1)+tau_n(2)*n_dir(2)+tau_n(3)*n_dir(3)
  sigma = tau_n(1)*n_f(1)+tau_n(2)*n_f(2)+tau_n(3)*n_f(3)

end subroutine compute_kernel


!========================================================================
      SUBROUTINE  UA0(X,Y,D,POT1,POT2,POT3,POT4,U)                      

      DOUBLE PRECISION :: X,Y,D,POT1,POT2,POT3,POT4,U 
      DOUBLE PRECISION :: F0,F1,F3,PI2                                       
      INTEGER :: I
      DOUBLE PRECISION :: DU
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,  &
                          A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ
      DIMENSION U(12),DU(12) 
                                          

!                                                                       

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   

!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    

!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,   &  

                UY,VY,WY,UZ,VZ,WZ                                      

      DATA F0,F1,F3/0.D0,1.D0,3.D0/                                     

      DATA PI2/6.283185307179586D0/                                     

!-----                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(POT1.NE.F0) THEN                                               

        DU( 1)= ALP1*Q/R3    +ALP2*X2*QR                                

        DU( 2)= ALP1*X/R3*SD +ALP2*XY*QR                                

        DU( 3)=-ALP1*X/R3*CD +ALP2*X*D*QR                               

        DU( 4)= X*QR*(-ALP1 +ALP2*(F1+A5) )                             

        DU( 5)= ALP1*A3/R3*SD +ALP2*Y*QR*A5                             

        DU( 6)=-ALP1*A3/R3*CD +ALP2*D*QR*A5                             

        DU( 7)= ALP1*(SD/R3-Y*QR) +ALP2*F3*X2/R5*UY                     

        DU( 8)= F3*X/R5*(-ALP1*Y*SD +ALP2*(Y*UY+Q) )                    

        DU( 9)= F3*X/R5*( ALP1*Y*CD +ALP2*D*UY )                        

        DU(10)= ALP1*(CD/R3+D*QR) +ALP2*F3*X2/R5*UZ                     

        DU(11)= F3*X/R5*( ALP1*D*SD +ALP2*Y*UZ )                        

        DU(12)= F3*X/R5*(-ALP1*D*CD +ALP2*(D*UZ-Q) )                    

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+POT1/PI2*DU(I)                                        

      ENDIF                                                             

!===================================                                    

!=====  DIP-SLIP CONTRIBUTION  =====                                    

!===================================                                    

      IF(POT2.NE.F0) THEN                                               

        DU( 1)=            ALP2*X*P*QR                                  

        DU( 2)= ALP1*S/R3 +ALP2*Y*P*QR                                  

        DU( 3)=-ALP1*T/R3 +ALP2*D*P*QR                                  

        DU( 4)=                 ALP2*P*QR*A5                            

        DU( 5)=-ALP1*F3*X*S/R5 -ALP2*Y*P*QRX                            

        DU( 6)= ALP1*F3*X*T/R5 -ALP2*D*P*QRX                            

        DU( 7)=                          ALP2*F3*X/R5*VY                

        DU( 8)= ALP1*(S2D/R3-F3*Y*S/R5) +ALP2*(F3*Y/R5*VY+P*QR)         

        DU( 9)=-ALP1*(C2D/R3-F3*Y*T/R5) +ALP2*F3*D/R5*VY                

        DU(10)=                          ALP2*F3*X/R5*VZ                

        DU(11)= ALP1*(C2D/R3+F3*D*S/R5) +ALP2*F3*Y/R5*VZ                

        DU(12)= ALP1*(S2D/R3-F3*D*T/R5) +ALP2*(F3*D/R5*VZ-P*QR)         

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+POT2/PI2*DU(I)                                        

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(POT3.NE.F0) THEN                                               

        DU( 1)= ALP1*X/R3 -ALP2*X*Q*QR                                  

        DU( 2)= ALP1*T/R3 -ALP2*Y*Q*QR                                  

        DU( 3)= ALP1*S/R3 -ALP2*D*Q*QR                                  

        DU( 4)= ALP1*A3/R3     -ALP2*Q*QR*A5                            

        DU( 5)=-ALP1*F3*X*T/R5 +ALP2*Y*Q*QRX                            

        DU( 6)=-ALP1*F3*X*S/R5 +ALP2*D*Q*QRX                            

        DU( 7)=-ALP1*F3*XY/R5           -ALP2*X*QR*WY                   

        DU( 8)= ALP1*(C2D/R3-F3*Y*T/R5) -ALP2*(Y*WY+Q)*QR               

        DU( 9)= ALP1*(S2D/R3-F3*Y*S/R5) -ALP2*D*QR*WY                   

        DU(10)= ALP1*F3*X*D/R5          -ALP2*X*QR*WZ                   

        DU(11)=-ALP1*(S2D/R3-F3*D*T/R5) -ALP2*Y*QR*WZ                   

        DU(12)= ALP1*(C2D/R3+F3*D*S/R5) -ALP2*(D*WZ-Q)*QR               

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+POT3/PI2*DU(I)                                        

      ENDIF                                                             

!=========================================                              

!=====  INFLATE SOURCE CONTRIBUTION  =====                              

!=========================================                              

      IF(POT4.NE.F0) THEN                                               

        DU( 1)=-ALP1*X/R3                                               

        DU( 2)=-ALP1*Y/R3                                               

        DU( 3)=-ALP1*D/R3                                               

        DU( 4)=-ALP1*A3/R3                                              

        DU( 5)= ALP1*F3*XY/R5                                           

        DU( 6)= ALP1*F3*X*D/R5                                          

        DU( 7)= DU(5)                                                   

        DU( 8)=-ALP1*B3/R3                                              

        DU( 9)= ALP1*F3*Y*D/R5                                          

        DU(10)=-DU(6)                                                   

        DU(11)=-DU(9)                                                   

        DU(12)= ALP1*C3/R3                                              

        DO 555 I=1,12                                                   

  555   U(I)=U(I)+POT4/PI2*DU(I)                                        

      ENDIF                                                             

      RETURN                                                            

      END SUBROUTINE  UA0                                                              

      SUBROUTINE  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    
      INTEGER :: I
      DOUBLE PRECISION :: X,Y,D,Z,POT1,POT2,POT3,POT4,U,DU  
      DOUBLE PRECISION :: F0,F1,F2,F3,F4,F5,F8,F9,PI2                                    
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,  &
                          A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ
      DIMENSION U(12),DU(12)                                            
      DOUBLE PRECISION :: C, RD, D12, D32, D33, D53, D54
      DOUBLE PRECISION :: FI1, FI2, FI3, FI4, FI5, &
                          FJ1, FJ2, FJ3, FJ4, FK1, FK2, FK3
!                                                                      

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   

!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  

!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3, &    

                 UY,VY,WY,UZ,VZ,WZ                                      

      DATA F0,F1,F2,F3,F4,F5,F8,F9   &                                   

             /0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,8.D0,9.D0/                 

      DATA PI2/6.283185307179586D0/                                     

!-----                                                                  

      C=D+Z                                                             

      RD=R+D                                                            

      D12=F1/(R*RD*RD)                                                  

      D32=D12*(F2*R+D)/R2                                               

      D33=D12*(F3*R+D)/(R2*RD)                                          

      D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)                           

      D54=D12*(F5*R2+F4*R*D+D2)/R3*D12                                  

!-----                                                                  

      FI1= Y*(D12-X2*D33)                                               

      FI2= X*(D12-Y2*D33)                                               

      FI3= X/R3-FI2                                                     

      FI4=-XY*D32                                                       

      FI5= F1/(R*RD)-X2*D32                                             

      FJ1=-F3*XY*(D33-X2*D54)                                           

      FJ2= F1/R3-F3*D12+F3*X2*Y2*D54                                    

      FJ3= A3/R3-FJ2                                                    

      FJ4=-F3*XY/R5-FJ1                                                 

      FK1=-Y*(D32-X2*D53)                                               

      FK2=-X*(D32-Y2*D53)                                               

      FK3=-F3*X*D/R5-FK2                                                

!-----                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(POT1.NE.F0) THEN                                               

        DU( 1)=-X2*QR  -ALP3*FI1*SD                                     

        DU( 2)=-XY*QR  -ALP3*FI2*SD                                     

        DU( 3)=-C*X*QR -ALP3*FI4*SD                                     

        DU( 4)=-X*QR*(F1+A5) -ALP3*FJ1*SD                               

        DU( 5)=-Y*QR*A5      -ALP3*FJ2*SD                               

        DU( 6)=-C*QR*A5      -ALP3*FK1*SD                               

        DU( 7)=-F3*X2/R5*UY      -ALP3*FJ2*SD                           

        DU( 8)=-F3*XY/R5*UY-X*QR -ALP3*FJ4*SD                           

        DU( 9)=-F3*C*X/R5*UY     -ALP3*FK2*SD                           

        DU(10)=-F3*X2/R5*UZ  +ALP3*FK1*SD                               

        DU(11)=-F3*XY/R5*UZ  +ALP3*FK2*SD                               

        DU(12)= F3*X/R5*(-C*UZ +ALP3*Y*SD)                              

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+POT1/PI2*DU(I)                                        

      ENDIF                                                             

!===================================                                    

!=====  DIP-SLIP CONTRIBUTION  =====                                    

!===================================                                    

      IF(POT2.NE.F0) THEN                                               

        DU( 1)=-X*P*QR +ALP3*FI3*SDCD                                   

        DU( 2)=-Y*P*QR +ALP3*FI1*SDCD                                   

        DU( 3)=-C*P*QR +ALP3*FI5*SDCD                                   

        DU( 4)=-P*QR*A5 +ALP3*FJ3*SDCD                                  

        DU( 5)= Y*P*QRX +ALP3*FJ1*SDCD                                  

        DU( 6)= C*P*QRX +ALP3*FK3*SDCD                                  

        DU( 7)=-F3*X/R5*VY      +ALP3*FJ1*SDCD                          

        DU( 8)=-F3*Y/R5*VY-P*QR +ALP3*FJ2*SDCD                          

        DU( 9)=-F3*C/R5*VY      +ALP3*FK1*SDCD                          

        DU(10)=-F3*X/R5*VZ -ALP3*FK3*SDCD                               

        DU(11)=-F3*Y/R5*VZ -ALP3*FK1*SDCD                               

        DU(12)=-F3*C/R5*VZ +ALP3*A3/R3*SDCD                             

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+POT2/PI2*DU(I)                                        

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(POT3.NE.F0) THEN                                               

        DU( 1)= X*Q*QR -ALP3*FI3*SDSD                                   

        DU( 2)= Y*Q*QR -ALP3*FI1*SDSD                                   

        DU( 3)= C*Q*QR -ALP3*FI5*SDSD                                   

        DU( 4)= Q*QR*A5 -ALP3*FJ3*SDSD                                  

        DU( 5)=-Y*Q*QRX -ALP3*FJ1*SDSD                                  

        DU( 6)=-C*Q*QRX -ALP3*FK3*SDSD                                  

        DU( 7)= X*QR*WY     -ALP3*FJ1*SDSD                              

        DU( 8)= QR*(Y*WY+Q) -ALP3*FJ2*SDSD                              

        DU( 9)= C*QR*WY     -ALP3*FK1*SDSD                              

        DU(10)= X*QR*WZ +ALP3*FK3*SDSD                                  

        DU(11)= Y*QR*WZ +ALP3*FK1*SDSD                                  

        DU(12)= C*QR*WZ -ALP3*A3/R3*SDSD                                

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+POT3/PI2*DU(I)                                        

      ENDIF                                                             

!=========================================                              

!=====  INFLATE SOURCE CONTRIBUTION  =====                              

!=========================================                              

      IF(POT4.NE.F0) THEN                                               

        DU( 1)= ALP3*X/R3                                               

        DU( 2)= ALP3*Y/R3                                               

        DU( 3)= ALP3*D/R3                                               

        DU( 4)= ALP3*A3/R3                                              

        DU( 5)=-ALP3*F3*XY/R5                                           

        DU( 6)=-ALP3*F3*X*D/R5                                          

        DU( 7)= DU(5)                                                   

        DU( 8)= ALP3*B3/R3                                              

        DU( 9)=-ALP3*F3*Y*D/R5                                          

        DU(10)=-DU(6)                                                   

        DU(11)=-DU(9)                                                   

        DU(12)=-ALP3*C3/R3                                              

        DO 555 I=1,12                                                   

  555   U(I)=U(I)+POT4/PI2*DU(I)                                        

      ENDIF                                                             

      RETURN                                                            

      END  SUBROUTINE  UB0                                                             

      SUBROUTINE  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,U)                    
      INTEGER :: I      
      DOUBLE PRECISION :: X,Y,D,Z,POT1,POT2,POT3,POT4,U,DU
      DOUBLE PRECISION :: F0,F1,F2,F3,F5,F7,F10,F15,PI2                                        
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,  &
                          A3,A5,B3,C3
      DIMENSION U(12),DU(12)                                            
      DOUBLE PRECISION :: C,Q2,R7,A7,B5,B7,C5,C7,D7,QR5,QR7,DR5 
!                                                                       

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   

!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  

!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3      

      DATA F0,F1,F2,F3,F5,F7,F10,F15     &                               

              /0.D0,1.D0,2.D0,3.D0,5.D0,7.D0,10.D0,15.D0/               

      DATA PI2/6.283185307179586D0/                                     

!-----                                                                  

      C=D+Z                                                             

      Q2=Q*Q                                                            

      R7=R5*R2                                                          

      A7=F1-F7*X2/R2                                                    

      B5=F1-F5*Y2/R2                                                    

      B7=F1-F7*Y2/R2                                                    

      C5=F1-F5*D2/R2                                                    

      C7=F1-F7*D2/R2                                                    

      D7=F2-F7*Q2/R2                                                    

      QR5=F5*Q/R2                                                       

      QR7=F7*Q/R2                                                       

      DR5=F5*D/R2                                                       

!-----                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(POT1.NE.F0) THEN                                               

        DU( 1)=-ALP4*A3/R3*CD  +ALP5*C*QR*A5                            

        DU( 2)= F3*X/R5*( ALP4*Y*CD +ALP5*C*(SD-Y*QR5) )                

        DU( 3)= F3*X/R5*(-ALP4*Y*SD +ALP5*C*(CD+D*QR5) )                

        DU( 4)= ALP4*F3*X/R5*(F2+A5)*CD   -ALP5*C*QRX*(F2+A7)           

        DU( 5)= F3/R5*( ALP4*Y*A5*CD +ALP5*C*(A5*SD-Y*QR5*A7) )         

        DU( 6)= F3/R5*(-ALP4*Y*A5*SD +ALP5*C*(A5*CD+D*QR5*A7) )         

        DU( 7)= DU(5)                                                   

        DU( 8)= F3*X/R5*( ALP4*B5*CD -ALP5*F5*C/R2*(F2*Y*SD+Q*B7) )     

        DU( 9)= F3*X/R5*(-ALP4*B5*SD +ALP5*F5*C/R2*(D*B7*SD-Y*C7*CD) )  

        DU(10)= F3/R5*   (-ALP4*D*A5*CD +ALP5*C*(A5*CD+D*QR5*A7) )      

        DU(11)= F15*X/R7*( ALP4*Y*D*CD  +ALP5*C*(D*B7*SD-Y*C7*CD) )     

        DU(12)= F15*X/R7*(-ALP4*Y*D*SD  +ALP5*C*(F2*D*CD-Q*C7) )        

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+POT1/PI2*DU(I)                                        

      ENDIF                                                             

!===================================                                    

!=====  DIP-SLIP CONTRIBUTION  =====                                    

!===================================                                    

      IF(POT2.NE.F0) THEN                                               

        DU( 1)= ALP4*F3*X*T/R5          -ALP5*C*P*QRX                   

        DU( 2)=-ALP4/R3*(C2D-F3*Y*T/R2) +ALP5*F3*C/R5*(S-Y*P*QR5)       

        DU( 3)=-ALP4*A3/R3*SDCD         +ALP5*F3*C/R5*(T+D*P*QR5)       

        DU( 4)= ALP4*F3*T/R5*A5              -ALP5*F5*C*P*QR/R2*A7      

        DU( 5)= F3*X/R5*(ALP4*(C2D-F5*Y*T/R2)-ALP5*F5*C/R2*(S-Y*P*QR7)) 

        DU( 6)= F3*X/R5*(ALP4*(F2+A5)*SDCD   -ALP5*F5*C/R2*(T+D*P*QR7)) 

        DU( 7)= DU(5)                                                   

        DU( 8)= F3/R5*(ALP4*(F2*Y*C2D+T*B5)    &                         

                                     +ALP5*C*(S2D-F10*Y*S/R2-P*QR5*B7)) 

        DU( 9)= F3/R5*(ALP4*Y*A5*SDCD-ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7)) 

        DU(10)= F3*X/R5*(-ALP4*(S2D-T*DR5) -ALP5*F5*C/R2*(T+D*P*QR7))   

        DU(11)= F3/R5*(-ALP4*(D*B5*C2D+Y*C5*S2D)    &                    

                                      -ALP5*C*((F3+A5)*C2D+Y*P*DR5*QR7))

        DU(12)= F3/R5*(-ALP4*D*A5*SDCD-ALP5*C*(S2D-F10*D*T/R2+P*QR5*C7))

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+POT2/PI2*DU(I)                                        

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(POT3.NE.F0) THEN                                               

        DU( 1)= F3*X/R5*(-ALP4*S +ALP5*(C*Q*QR5-Z))                     

        DU( 2)= ALP4/R3*(S2D-F3*Y*S/R2)+ALP5*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)

        DU( 3)=-ALP4/R3*(F1-A3*SDSD)   -ALP5*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)

        DU( 4)=-ALP4*F3*S/R5*A5 +ALP5*(C*QR*QR5*A7-F3*Z/R5*A5)          

        DU( 5)= F3*X/R5*(-ALP4*(S2D-F5*Y*S/R2)    &                      

                                     -ALP5*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z)) 

        DU( 6)= F3*X/R5*( ALP4*(F1-(F2+A5)*SDSD)    &                    

                                     +ALP5*F5/R2*(C*(S-D+D*Q*QR7)-D*Z)) 

        DU( 7)= DU(5)                                                   

        DU( 8)= F3/R5*(-ALP4*(F2*Y*S2D+S*B5)    &                        

                      -ALP5*(C*(F2*SDSD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5)) 

        DU( 9)= F3/R5*( ALP4*Y*(F1-A5*SDSD)     &                        

                      +ALP5*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))             

        DU(10)= F3*X/R5*(-ALP4*(C2D+S*DR5)       &                       

                     +ALP5*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))            

        DU(11)= F3/R5*( ALP4*(D*B5*S2D-Y*C5*C2D)    &                    

                     +ALP5*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))     

        DU(12)= F3/R5*(-ALP4*D*(F1-A5*SDSD)      &                       

                     -ALP5*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5))) 

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+POT3/PI2*DU(I)                                        

      ENDIF                                                             

!=========================================                              

!=====  INFLATE SOURCE CONTRIBUTION  =====                              

!=========================================                              

      IF(POT4.NE.F0) THEN                                               

        DU( 1)= ALP4*F3*X*D/R5                                          

        DU( 2)= ALP4*F3*Y*D/R5                                          

        DU( 3)= ALP4*C3/R3                                              

        DU( 4)= ALP4*F3*D/R5*A5                                         

        DU( 5)=-ALP4*F15*XY*D/R7                                        

        DU( 6)=-ALP4*F3*X/R5*C5                                         

        DU( 7)= DU(5)                                                   

        DU( 8)= ALP4*F3*D/R5*B5                                         

        DU( 9)=-ALP4*F3*Y/R5*C5                                         

        DU(10)= DU(6)                                                   

        DU(11)= DU(9)                                                   

        DU(12)= ALP4*F3*D/R5*(F2+C5)                                    

        DO 555 I=1,12                                                   

  555   U(I)=U(I)+POT4/PI2*DU(I)                                        

      ENDIF                                                             

      RETURN                                                            

      END  SUBROUTINE  UC0                                                             

      SUBROUTINE  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,    &                       

                    AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,    &              

                    UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)  


      DOUBLE PRECISION ::  ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,  &
                           DISL1,DISL2,DISL3,UX,UY,UZ,UXX,UYX,  &
                           UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ             
      INTEGER :: I,IRET,KXI,KET,K,J
      DOUBLE PRECISION ::  F0,EPS 
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD
      DOUBLE PRECISION :: U,DUA,DUB,DUC,AALPHA,DDIP,ZZ,DD1,DD2,DD3,  &
                          XI,D,P,Q,ET,R12,R21,R22,DU
!                                                                       

!********************************************************************   

!*****                                                          *****   

!*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   

!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   

!*****              CODED BY  Y.OKADA ... SEP.1991              *****   

!*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****   

!*****                          JUL.1993                        *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           

!*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          

!*****   DEPTH : DEPTH OF REFERENCE POINT                               

!*****   DIP   : DIP-ANGLE (DEGREE)                                     

!*****   AL1,AL2   : FAULT LENGTH RANGE                                 

!*****   AW1,AW2   : FAULT WIDTH RANGE                                  

!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              

!                                                                       

!***** OUTPUT                                                           

!*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               

!*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             

!*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )

!*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     

!*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD                                         

      DIMENSION  XI(2),ET(2),KXI(2),KET(2)                              

      DIMENSION  U(12),DU(12),DUA(12),DUB(12),DUC(12)                   

      DATA  F0,EPS/ 0.D0, 1.D-6 /                                       

!-----                                                                  

      IF(Z.GT.0.) WRITE(6,'('' ** POSITIVE Z WAS GIVEN IN SUB-DC3D'')') 

      DO 111 I=1,12                                                     

        U  (I)=F0                                                       

        DUA(I)=F0                                                       

        DUB(I)=F0                                                       

        DUC(I)=F0                                                       

  111 CONTINUE                                                          

      AALPHA=ALPHA                                                      

      DDIP=DIP                                                          

      CALL DCCON0(AALPHA,DDIP)                                          

!-----                                                                  

      ZZ=Z                                                              

      DD1=DISL1                                                         

      DD2=DISL2                                                         

      DD3=DISL3                                                         

      XI(1)=X-AL1                                                       

      XI(2)=X-AL2                                                       

      IF(ABS(XI(1)).LT.EPS) XI(1)=F0                                   

      IF(ABS(XI(2)).LT.EPS) XI(2)=F0                                   

!======================================                                 

!=====  REAL-SOURCE CONTRIBUTION  =====                                 

!======================================                                 

      D=DEPTH+Z                                                         

      P=Y*CD+D*SD                                                       

      Q=Y*SD-D*CD                                                       

      ET(1)=P-AW1                                                       

      ET(2)=P-AW2                                                       

      IF(ABS(Q).LT.EPS)  Q=F0                                          

      IF(ABS(ET(1)).LT.EPS) ET(1)=F0                                   

      IF(ABS(ET(2)).LT.EPS) ET(2)=F0                                   

!--------------------------------                                       

!----- REJECT SINGULAR CASE -----                                       

!--------------------------------                                       

!----- ON FAULT EDGE                                                    

      IF(Q.EQ.F0      &                                                  

         .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)    &      

               .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))   &    

         GO TO 99                                                       

!----- ON NEGATIVE EXTENSION OF FAULT EDGE                              

      KXI(1)=0                                                          

      KXI(2)=0                                                          

      KET(1)=0                                                          

      KET(2)=0                                                          

      R12=SQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)                            

      R21=SQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)                            

      R22=SQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)                            

      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   

      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   

      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   

      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   

!=====                                                                  

      DO 223 K=1,2                                                      

      DO 222 J=1,2                                                      

        CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  

        CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          

!-----                                                                  

        DO 220 I=1,10,3                                                 

          DU(I)  =-DUA(I)                                               

          DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD                              

          DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD                              

          IF(I.LT.10) GO TO 220                                         

          DU(I)  =-DU(I)                                                

          DU(I+1)=-DU(I+1)                                              

          DU(I+2)=-DU(I+2)                                              

  220   CONTINUE                                                        

        DO 221 I=1,12                                                   

          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  

          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  

  221   CONTINUE                                                        

!-----                                                                  

  222 CONTINUE                                                          

  223 CONTINUE                                                          

!=======================================                                

!=====  IMAGE-SOURCE CONTRIBUTION  =====                                

!=======================================                                

      D=DEPTH-Z                                                         

      P=Y*CD+D*SD                                                       

      Q=Y*SD-D*CD                                                       

      ET(1)=P-AW1                                                       

      ET(2)=P-AW2                                                       

      IF(ABS(Q).LT.EPS)  Q=F0                                          

      IF(ABS(ET(1)).LT.EPS) ET(1)=F0                                   

      IF(ABS(ET(2)).LT.EPS) ET(2)=F0                                   

!--------------------------------                                       

!----- REJECT SINGULAR CASE -----                                       

!--------------------------------                                       

!----- ON FAULT EDGE                                                    

      IF(Q.EQ.F0        &                                                

         .AND.(    (XI(1)*XI(2).LE.F0 .AND. ET(1)*ET(2).EQ.F0)     &      

               .OR.(ET(1)*ET(2).LE.F0 .AND. XI(1)*XI(2).EQ.F0) ))     &  

         GO TO 99                                                       

!----- ON NEGATIVE EXTENSION OF FAULT EDGE                              

      KXI(1)=0                                                          

      KXI(2)=0                                                          

      KET(1)=0                                                          

      KET(2)=0                                                          

      R12=SQRT(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)                            

      R21=SQRT(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)                            

      R22=SQRT(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)                            

      IF(XI(1).LT.F0 .AND. R21+XI(2).LT.EPS) KXI(1)=1                   

      IF(XI(1).LT.F0 .AND. R22+XI(2).LT.EPS) KXI(2)=1                   

      IF(ET(1).LT.F0 .AND. R12+ET(2).LT.EPS) KET(1)=1                   

      IF(ET(1).LT.F0 .AND. R22+ET(2).LT.EPS) KET(2)=1                   

!=====                                                                  

      DO 334 K=1,2                                                      

      DO 333 J=1,2                                                      

        CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))                  

        CALL UA(XI(J),ET(K),Q,DD1,DD2,DD3,DUA)                          

        CALL UB(XI(J),ET(K),Q,DD1,DD2,DD3,DUB)                          

        CALL UC(XI(J),ET(K),Q,ZZ,DD1,DD2,DD3,DUC)                       

!-----                                                                  

        DO 330 I=1,10,3                                                 

          DU(I)=DUA(I)+DUB(I)+Z*DUC(I)                                  

          DU(I+1)=(DUA(I+1)+DUB(I+1)+Z*DUC(I+1))*CD        &             

                 -(DUA(I+2)+DUB(I+2)+Z*DUC(I+2))*SD                     

          DU(I+2)=(DUA(I+1)+DUB(I+1)-Z*DUC(I+1))*SD        &             

                 +(DUA(I+2)+DUB(I+2)-Z*DUC(I+2))*CD                     

          IF(I.LT.10) GO TO 330                                         

          DU(10)=DU(10)+DUC(1)                                          

          DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD                             

          DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD                             

  330   CONTINUE                                                        

        DO 331 I=1,12                                                   

          IF(J+K.NE.3) U(I)=U(I)+DU(I)                                  

          IF(J+K.EQ.3) U(I)=U(I)-DU(I)                                  

  331   CONTINUE                                                        

!-----                                                                  

  333 CONTINUE                                                          

  334 CONTINUE                                                          

!=====                                                                  

      UX=U(1)                                                           

      UY=U(2)                                                           

      UZ=U(3)                                                           

      UXX=U(4)                                                          

      UYX=U(5)                                                          

      UZX=U(6)                                                          

      UXY=U(7)                                                          

      UYY=U(8)                                                          

      UZY=U(9)                                                          

      UXZ=U(10)                                                         

      UYZ=U(11)                                                         

      UZZ=U(12)                                                         

      IRET=0                                                            

      RETURN                                                            

!===========================================                            

!=====  IN CASE OF SINGULAR (ON EDGE)  =====                            

!===========================================                            

   99 UX=F0                                                             

      UY=F0                                                             

      UZ=F0                                                             

      UXX=F0                                                            

      UYX=F0                                                            

      UZX=F0                                                            

      UXY=F0                                                            

      UYY=F0                                                            

      UZY=F0                                                            

      UXZ=F0                                                            

      UYZ=F0                                                            

      UZZ=F0                                                            

      IRET=1                                                            

      RETURN                                                            

      END  SUBROUTINE  DC3D                                                             

      SUBROUTINE  UA(XI,ET,Q,DISL1,DISL2,DISL3,U) 
  
      INTEGER :: I                
      DOUBLE PRECISION :: XI,ET,Q,DISL1,DISL2,DISL3,U
      DOUBLE PRECISION :: F0,F2,PI2
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,  &
                          Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ
      DOUBLE PRECISION :: XY,QX,QY,DU
      DIMENSION U(12),DU(12)                                            
 
!                                                                       

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   

!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  

!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  &

                 EY,EZ,FY,FZ,GY,GZ,HY,HZ                                

      DATA F0,F2,PI2/0.D0,2.D0,6.283185307179586D0/                     

!-----                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

      XY=XI*Y11                                                         

      QX=Q *X11                                                         

      QY=Q *Y11                                                         

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(DISL1.NE.F0) THEN                                              

        DU( 1)=    TT/F2 +ALP2*XI*QY                                    

        DU( 2)=           ALP2*Q/R                                      

        DU( 3)= ALP1*ALE -ALP2*Q*QY                                     

        DU( 4)=-ALP1*QY  -ALP2*XI2*Q*Y32                                

        DU( 5)=          -ALP2*XI*Q/R3                                  

        DU( 6)= ALP1*XY  +ALP2*XI*Q2*Y32                                

        DU( 7)= ALP1*XY*SD        +ALP2*XI*FY+D/F2*X11                  

        DU( 8)=                    ALP2*EY                              

        DU( 9)= ALP1*(CD/R+QY*SD) -ALP2*Q*FY                            

        DU(10)= ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11                  

        DU(11)=                    ALP2*EZ                              

        DU(12)=-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ                            

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       

      ENDIF                                                             

!======================================                                 

!=====    DIP-SLIP CONTRIBUTION   =====                                 

!======================================                                 

      IF(DISL2.NE.F0) THEN                                              

        DU( 1)=           ALP2*Q/R                                      

        DU( 2)=    TT/F2 +ALP2*ET*QX                                    

        DU( 3)= ALP1*ALX -ALP2*Q*QX                                     

        DU( 4)=        -ALP2*XI*Q/R3                                    

        DU( 5)= -QY/F2 -ALP2*ET*Q/R3                                    

        DU( 6)= ALP1/R +ALP2*Q2/R3                                      

        DU( 7)=                      ALP2*EY                            

        DU( 8)= ALP1*D*X11+XY/F2*SD +ALP2*ET*GY                         

        DU( 9)= ALP1*Y*X11          -ALP2*Q*GY                          

        DU(10)=                      ALP2*EZ                            

        DU(11)= ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ                         

        DU(12)=-ALP1*D*X11          -ALP2*Q*GZ                          

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(DISL3.NE.F0) THEN                                              

        DU( 1)=-ALP1*ALE -ALP2*Q*QY                                     

        DU( 2)=-ALP1*ALX -ALP2*Q*QX                                     

        DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY)                            

        DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32                                

        DU( 5)=-ALP1/R   +ALP2*Q2/R3                                    

        DU( 6)=-ALP1*QY  -ALP2*Q*Q2*Y32                                 

        DU( 7)=-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY                           

        DU( 8)=-ALP1*Y*X11         -ALP2*Q*GY                           

        DU( 9)= ALP1*(D*X11+XY*SD) +ALP2*Q*HY                           

        DU(10)= ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ                           

        DU(11)= ALP1*D*X11         -ALP2*Q*GZ                           

        DU(12)= ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ                           

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       

      ENDIF                                                             

      RETURN                                                            

      END   SUBROUTINE  UA                                                            

      SUBROUTINE  UB(XI,ET,Q,DISL1,DISL2,DISL3,U)                       
      INTEGER :: I
      DOUBLE PRECISION :: XI,ET,Q,DISL1,DISL2,DISL3,U
      DOUBLE PRECISION :: F0,F1,F2,PI2
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,  &
                          Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ 
      DOUBLE PRECISION :: RD,D11,AJ2,AJ5,AI4,X,AI3,AK1,AK3,AJ3,AJ6,   &
                          RD2,XY,AI1,AI2,AK2,AK4,AJ1,AJ4,QX,QY,DU 
      DIMENSION U(12),DU(12)                                            

!                                                                       

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   

!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  

!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  &

                 EY,EZ,FY,FZ,GY,GZ,HY,HZ                                

      DATA  F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/            

!-----                                                                  

      RD=R+D                                                            

      D11=F1/(R*RD)                                                     

      AJ2=XI*Y/RD*D11                                                   

      AJ5=-(D+Y*Y/RD)*D11                                               

      IF(CD.NE.F0) THEN                                                 

        IF(XI.EQ.F0) THEN                                               

          AI4=F0                                                        

        ELSE                                                            

          X=SQRT(XI2+Q2)                                               

          AI4=F1/CDCD*( XI/RD*SDCD          &                            

             +F2*ATAN((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )        

        ENDIF                                                           

        AI3=(Y*CD/RD-ALE+SD*LOG(RD))/CDCD                              

        AK1=XI*(D11-Y11*SD)/CD                                          

        AK3=(Q*Y11-Y*D11)/CD                                            

        AJ3=(AK1-AJ2*SD)/CD                                             

        AJ6=(AK3-AJ5*SD)/CD                                             

      ELSE                                                              

        RD2=RD*RD                                                       

        AI3=(ET/RD+Y*Q/RD2-ALE)/F2                                      

        AI4=XI*Y/RD2/F2                                                 

        AK1=XI*Q/RD*D11                                                 

        AK3=SD/RD*(XI2*D11-F1)                                          

        AJ3=-XI/RD2*(Q2*D11-F1/F2)                                      

        AJ6=-Y/RD2*(XI2*D11-F1/F2)                                      

      ENDIF                                                             

!-----                                                                  

      XY=XI*Y11                                                         

      AI1=-XI/RD*CD-AI4*SD                                              

      AI2= LOG(RD)+AI3*SD                                              

      AK2= F1/R+AK3*SD                                                  

      AK4= XY*CD-AK1*SD                                                 

      AJ1= AJ5*CD-AJ6*SD                                                

      AJ4=-XY-AJ2*CD+AJ3*SD                                             

!=====                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

      QX=Q*X11                                                          

      QY=Q*Y11                                                          

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(DISL1.NE.F0) THEN                                              

        DU( 1)=-XI*QY-TT -ALP3*AI1*SD                                   

        DU( 2)=-Q/R      +ALP3*Y/RD*SD                                  

        DU( 3)= Q*QY     -ALP3*AI2*SD                                   

        DU( 4)= XI2*Q*Y32 -ALP3*AJ1*SD                                  

        DU( 5)= XI*Q/R3   -ALP3*AJ2*SD                                  

        DU( 6)=-XI*Q2*Y32 -ALP3*AJ3*SD                                  

        DU( 7)=-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD                           

        DU( 8)=-EY          +ALP3*(F1/R+AJ5)*SD                         

        DU( 9)= Q*FY        -ALP3*(QY-AJ6)*SD                           

        DU(10)=-XI*FZ-Y*X11 +ALP3*AK1*SD                                

        DU(11)=-EZ          +ALP3*Y*D11*SD                              

        DU(12)= Q*FZ        +ALP3*AK2*SD                                

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       

      ENDIF                                                             

!======================================                                 

!=====    DIP-SLIP CONTRIBUTION   =====                                 

!======================================                                 

      IF(DISL2.NE.F0) THEN                                              

        DU( 1)=-Q/R      +ALP3*AI3*SDCD                                 

        DU( 2)=-ET*QX-TT -ALP3*XI/RD*SDCD                               

        DU( 3)= Q*QX     +ALP3*AI4*SDCD                                 

        DU( 4)= XI*Q/R3     +ALP3*AJ4*SDCD                              

        DU( 5)= ET*Q/R3+QY  +ALP3*AJ5*SDCD                              

        DU( 6)=-Q2/R3       +ALP3*AJ6*SDCD                              

        DU( 7)=-EY          +ALP3*AJ1*SDCD                              

        DU( 8)=-ET*GY-XY*SD +ALP3*AJ2*SDCD                              

        DU( 9)= Q*GY        +ALP3*AJ3*SDCD                              

        DU(10)=-EZ          -ALP3*AK3*SDCD                              

        DU(11)=-ET*GZ-XY*CD -ALP3*XI*D11*SDCD                           

        DU(12)= Q*GZ        -ALP3*AK4*SDCD                              

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(DISL3.NE.F0) THEN                                              

        DU( 1)= Q*QY           -ALP3*AI3*SDSD                           

        DU( 2)= Q*QX           +ALP3*XI/RD*SDSD                         

        DU( 3)= ET*QX+XI*QY-TT -ALP3*AI4*SDSD                           

        DU( 4)=-XI*Q2*Y32 -ALP3*AJ4*SDSD                                

        DU( 5)=-Q2/R3     -ALP3*AJ5*SDSD                                

        DU( 6)= Q*Q2*Y32  -ALP3*AJ6*SDSD                                

        DU( 7)= Q*FY -ALP3*AJ1*SDSD                                     

        DU( 8)= Q*GY -ALP3*AJ2*SDSD                                     

        DU( 9)=-Q*HY -ALP3*AJ3*SDSD                                     

        DU(10)= Q*FZ +ALP3*AK3*SDSD                                     

        DU(11)= Q*GZ +ALP3*XI*D11*SDSD                                  

        DU(12)=-Q*HZ +ALP3*AK4*SDSD                                     

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       

      ENDIF                                                             

      RETURN                                                            

      END  SUBROUTINE  UB                                                             

      SUBROUTINE  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,U)                     
      INTEGER :: I
      DOUBLE PRECISION :: XI,ET,Q,Z,DISL1,DISL2,DISL3,U
      DOUBLE PRECISION :: F0,F1,F2,F3,PI2
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,  &
                          Y11,X32,Y32,EY,EZ,FY,FZ,GY,GZ,HY,HZ
      DOUBLE PRECISION :: C,X53,Y53,H,Z32,Z53,Y0,Z0,PPY,PPZ,   &
                          QQ,QQY,QQZ,XY,QX,QY,QR,CQX,CDR,YY0,DU 
      DIMENSION U(12),DU(12)                                            

!                                                                       

!********************************************************************   

!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   

!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   

!********************************************************************   

!                                                                       

!***** INPUT                                                            

!*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              

!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              

!***** OUTPUT                                                           

!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32, & 

                 EY,EZ,FY,FZ,GY,GZ,HY,HZ                                

      DATA F0,F1,F2,F3,PI2/0.D0,1.D0,2.D0,3.D0,6.283185307179586D0/     

!-----                                                                  

      C=D+Z                                                             

      X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2                     

      Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2                     

      H=Q*CD-Z                                                          

      Z32=SD/R3-H*Y32                                                   

      Z53=F3*SD/R5-H*Y53                                                

      Y0=Y11-XI2*Y32                                                    

      Z0=Z32-XI2*Z53                                                    

      PPY=CD/R3+Q*Y32*SD                                                

      PPZ=SD/R3-Q*Y32*CD                                                

      QQ=Z*Y32+Z32+Z0                                                   

      QQY=F3*C*D/R5-QQ*SD                                               

      QQZ=F3*C*Y/R5-QQ*CD+Q*Y32                                         

      XY=XI*Y11                                                         

      QX=Q*X11                                                          

      QY=Q*Y11                                                          

      QR=F3*Q/R5                                                        

      CQX=C*Q*X53                                                       

      CDR=(C+D)/R3                                                      

      YY0=Y/R3-Y0*CD                                                    

!=====                                                                  

      DO 111  I=1,12                                                    

  111 U(I)=F0                                                           

!======================================                                 

!=====  STRIKE-SLIP CONTRIBUTION  =====                                 

!======================================                                 

      IF(DISL1.NE.F0) THEN                                              

        DU( 1)= ALP4*XY*CD           -ALP5*XI*Q*Z32                     

        DU( 2)= ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3                       

        DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32)      

        DU( 4)= ALP4*Y0*CD                  -ALP5*Q*Z0                  

        DU( 5)=-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR               

        DU( 6)=-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ)    

        DU( 7)=-ALP4*XI*PPY*CD    -ALP5*XI*QQY                          

        DU( 8)= ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD      &                   

                                  -ALP5*(CDR*SD-ET/R3-C*Y*QR)           

        DU( 9)=-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD) 

        DU(10)= ALP4*XI*PPZ*CD    -ALP5*XI*QQZ                          

        DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR)   

        DU(12)=         YY0*CD    -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD) 

        DO 222 I=1,12                                                   

  222   U(I)=U(I)+DISL1/PI2*DU(I)                                       

      ENDIF                                                             

!======================================                                 

!=====    DIP-SLIP CONTRIBUTION   =====                                 

!======================================                                 

      IF(DISL2.NE.F0) THEN                                              

        DU( 1)= ALP4*CD/R -QY*SD -ALP5*C*Q/R3                           

        DU( 2)= ALP4*Y*X11       -ALP5*C*ET*Q*X32                       

        DU( 3)=     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32)                   

        DU( 4)=-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD                

        DU( 5)=-ALP4*Y/R3     +ALP5*C*ET*QR                             

        DU( 6)=    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2)                  

        DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR)                

        DU( 8)= ALP4*(X11-Y*Y*X32) -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53) 

        DU( 9)=  XI*PPY*SD+Y*D*X32 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)   

        DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR)                

        DU(11)= ALP4*Y*D*X32       -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53) 

        DU(12)=-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53) 

        DO 333 I=1,12                                                   

  333   U(I)=U(I)+DISL2/PI2*DU(I)                                       

      ENDIF                                                             

!========================================                               

!=====  TENSILE-FAULT CONTRIBUTION  =====                               

!========================================                               

      IF(DISL3.NE.F0) THEN                                              

        DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32)                

        DU( 2)= ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32)                

        DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32)           

        DU( 4)= ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0)

        DU( 5)= ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2)             

        DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0)                 

        DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD)    

        DU( 8)=-ALP4*F2*XI*PPY*SD-Y*D*X32      &                         

                          +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53)            

        DU( 9)=-ALP4*(XI*PPY*CD-X11+Y*Y*X32)     &                       

                          +ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY) 

        DU(10)=  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD)  

        DU(11)= ALP4*F2*XI*PPZ*SD-X11+D*D*X32    &                       

                          -ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53)            

        DU(12)= ALP4*(XI*PPZ*CD+Y*D*X32)     &                           

                          +ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ) 

        DO 444 I=1,12                                                   

  444   U(I)=U(I)+DISL3/PI2*DU(I)                                       

      ENDIF                                                             

      RETURN                                                            

      END   SUBROUTINE  UC                                                            

      SUBROUTINE  DCCON0(ALPHA,DIP) 

      DOUBLE PRECISION :: ALPHA,DIP                                    
      DOUBLE PRECISION :: F0,F1,F2,PI2,EPS
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,  &
                          CDCD,SDCD,S2D,C2D
      DOUBLE PRECISION :: P18
!                                                                       

!*******************************************************************    

!*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    

!*******************************************************************    

!                                                                       

!***** INPUT                                                            

!*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           

!*****   DIP   : DIP-ANGLE (DEGREE)                                     

!### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   

!                                                                       

      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD,SDSD,CDCD,SDCD,S2D,C2D  

      DATA F0,F1,F2,PI2/0.D0,1.D0,2.D0,6.283185307179586D0/             

      DATA EPS/1.D-6/                                                   

!-----                                                                  

      ALP1=(F1-ALPHA)/F2                                                

      ALP2= ALPHA/F2                                                    

      ALP3=(F1-ALPHA)/ALPHA                                             

      ALP4= F1-ALPHA                                                    

      ALP5= ALPHA                                                       

!-----                                                                  

      P18=PI2/360.D0                                                    

      SD=SIN(DIP*P18)                                                  

      CD=COS(DIP*P18)                                                  

      IF(ABS(CD).LT.EPS) THEN                                          

        CD=F0                                                           

        IF(SD.GT.F0) SD= F1                                             

        IF(SD.LT.F0) SD=-F1                                             

      ENDIF                                                             

      SDSD=SD*SD                                                        

      CDCD=CD*CD                                                        
      SDCD=SD*CD                                                        
      S2D=F2*SDCD                                                       
      C2D=CDCD-SDSD                                                     
      RETURN                                                            
      END SUBROUTINE DCCON0  

                                                          
      SUBROUTINE  DCCON1(X,Y,D)  
                                       
      DOUBLE PRECISION :: X,Y,D
      DOUBLE PRECISION :: F0,F1,F3,F5,EPS
      DOUBLE PRECISION :: ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD
      DOUBLE PRECISION :: P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,  &
                          A3,A5,B3,C3,UY,VY,WY,UZ,VZ,WZ 
      DOUBLE PRECISION :: R7                                      
!                                                                       
!********************************************************************** 
!*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 
!********************************************************************** 
!                                                                       
!***** INPUT                                                            
!*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    
!### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  
!                                                                       
      COMMON /C0/ALP1,ALP2,ALP3,ALP4,ALP5,SD,CD                                         
      COMMON /C1/P,Q,S,T,XY,X2,Y2,D2,R,R2,R3,R5,QR,QRX,A3,A5,B3,C3,  &    
                 UY,VY,WY,UZ,VZ,WZ                                      
      DATA  F0,F1,F3,F5,EPS/0.D0,1.D0,3.D0,5.D0,1.D-6/                  
!-----                                                                  
      IF(ABS(X).LT.EPS) X=F0                                           
      IF(ABS(Y).LT.EPS) Y=F0                                          
      IF(ABS(D).LT.EPS) D=F0                                           
      P=Y*CD+D*SD                                                       
      Q=Y*SD-D*CD                                                       
      S=P*SD+Q*CD                                                       
      T=P*CD-Q*SD                                                       
      XY=X*Y                                                            
      X2=X*X                                                            
      Y2=Y*Y                                                            
      D2=D*D                                                            
      R2=X2+Y2+D2                                                       
      R =SQRT(R2)                                                      
      IF(R.EQ.F0) RETURN                                                
      R3=R *R2                                                          
      R5=R3*R2                                                          
      R7=R5*R2                                                          
!-----                                                                  
      A3=F1-F3*X2/R2                                                    
      A5=F1-F5*X2/R2                                                    
      B3=F1-F3*Y2/R2                                                    
      C3=F1-F3*D2/R2                                                    
!-----                                                                  
      QR=F3*Q/R5                                                        
      QRX=F5*QR*X/R2                                                    
!-----                                                                  
      UY=SD-F5*Y*Q/R2                                                   
      UZ=CD+F5*D*Q/R2                                                   
      VY=S -F5*Y*P*Q/R2                                                 
      VZ=T +F5*D*P*Q/R2                                                 
      WY=UY+SD                                                          
      WZ=UZ+CD                                                          
      RETURN                                                            
      END  SUBROUTINE  DCCON1


                                                             
      SUBROUTINE  DCCON2(XI,ET,Q,SD,CD,KXI,KET) 
                        
      DOUBLE PRECISION :: XI,ET,Q,SD,CD
      DOUBLE PRECISION :: F0,F1,F2,EPS
      DOUBLE PRECISION :: XI2,ET2,Q2,R2,R,R3,R5,Y,D,TT,ALX,X11,X32,   &
                          RXI,ALE,Y11,Y32,RET,EY,EZ,FY,FZ,GY,GZ,HY,HZ
      INTEGER :: KXI,KET                                       
!                                                                       
!********************************************************************** 
!*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 
!********************************************************************** 
!                                                                       
!***** INPUT                                                            
!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
!*****   SD,CD   : SIN, COS OF DIP-ANGLE                                
!*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY  
!                                                                       
!### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER01
!                                                                       
      COMMON /C2/XI2,ET2,Q2,R,R2,R3,R5,Y,D,TT,ALX,ALE,X11,Y11,X32,Y32,  &
                 EY,EZ,FY,FZ,GY,GZ,HY,HZ                                
      DATA  F0,F1,F2,EPS/0.D0,1.D0,2.D0,1.D-6/                          
!-----                                                                  
      IF(ABS(XI).LT.EPS) XI=F0                                         
      IF(ABS(ET).LT.EPS) ET=F0                                         
      IF(ABS( Q).LT.EPS)  Q=F0                                         
      XI2=XI*XI                                                         
      ET2=ET*ET                                                         
      Q2=Q*Q                                                            
      R2=XI2+ET2+Q2                                                     
      R =SQRT(R2)                                                      
      IF(R.EQ.F0) RETURN                                                
      R3=R *R2                                                          
      R5=R3*R2                                                          
      Y =ET*CD+Q*SD                                                     
      D =ET*SD-Q*CD                                                     
!-----                                                                  
      IF(Q.EQ.F0) THEN                                                  
        TT=F0                                                           
      ELSE                                                              
        TT=ATAN(XI*ET/(Q*R))                                           
      ENDIF                                                             
!-----                                                                  
      IF(KXI.EQ.1) THEN                                                 
        ALX=-LOG(R-XI)                                                 
        X11=F0                                                          
        X32=F0                                                          
      ELSE                                                              
        RXI=R+XI                                                        
        ALX=LOG(RXI)                                                   

        X11=F1/(R*RXI)                                                  

        X32=(R+RXI)*X11*X11/R                                           

      ENDIF                                                             

!-----                                                                  

      IF(KET.EQ.1) THEN                                                 

        ALE=-LOG(R-ET)                                                 

        Y11=F0                                                          

        Y32=F0                                                          

      ELSE                                                              

        RET=R+ET                                                        

        ALE=LOG(RET)                                                   

        Y11=F1/(R*RET)                                                  

        Y32=(R+RET)*Y11*Y11/R                                           

      ENDIF                                                             

!-----                                                                  

      EY=SD/R-Y*Q/R3                                                    

      EZ=CD/R+D*Q/R3                                                    

      FY=D/R3+XI2*Y32*SD                                                

      FZ=Y/R3+XI2*Y32*CD                                                

      GY=F2*X11*SD-Y*Q*X32                                              

      GZ=F2*X11*CD+D*Q*X32                                              

      HY=D*Q*X32+XI*Q*Y32*SD                                            

      HZ=Y*Q*X32+XI*Q*Y32*CD                                            

      RETURN                                                            

      END  SUBROUTINE  DCCON2                                                             

end module okada
