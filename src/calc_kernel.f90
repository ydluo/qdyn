! computer_kernel

module calc_kernel

  implicit none
  private

  public :: compute_kernel

contains
!=====================================================================
! read in S... Source; O... Observation
! unit response: U=1, return tau
subroutine compute_kernel(LAM,MU,SX,SY,SZ,S_DIP,L,W,OX,OY,OZ,O_DIP,IRET,tau)

  use constants, only : PI
  use dc3d_all 

  double precision :: SX,SY,SZ,OX,OY,OZ
  double precision :: LAM, MU, ALPHA,   & 
    S_DEPTH,S_DIP, L, W, U,&
    X, Y, Z, O_DIP, &
    UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
  double precision :: STRESS(3,3), STRAIN(3,3),TR, n_f(3), n_dir(3), tau_n(3) 
  integer, intent(inout) :: IRET
  double precision, intent(inout) :: tau
  
  ALPHA = (LAM+MU)/(LAM+2d0*MU)
  S_DEPTH = -1d0*SZ

  U = 1d0 
  X = OX-SX
  Y = OY-SY
  Z = OZ
 
  call   DC3D(ALPHA,X,Y,Z,S_DEPTH,S_DIP,-0.5d0*L,0.5d0*L,-0.5d0*W,0.5d0*W,0d0,U,0d0,   &
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
  
  n_f(1) = 0d0
  n_f(2) = -1d0*dsin(O_DIP/180d0*PI)
  n_f(3) = 1d0*dcos(O_DIP/180d0*PI)
  n_dir(1) = 0d0
  n_dir(2) = -1d0*dcos(O_DIP/180d0*PI)
  n_dir(3) = -1d0*dsin(O_DIP/180d0*PI)
 
  tau_n(1) = STRESS(1,1)*n_f(1)+STRESS(1,2)*n_f(2)+STRESS(1,3)*n_f(3)
  tau_n(2) = STRESS(2,1)*n_f(1)+STRESS(2,2)*n_f(2)+STRESS(2,3)*n_f(3)
  tau_n(3) = STRESS(3,1)*n_f(1)+STRESS(3,2)*n_f(2)+STRESS(3,3)*n_f(3)

  tau = tau_n(1)*n_DIR(1)+tau_n(2)*n_DIR(2)+tau_n(3)*n_DIR(3)
  

  
end subroutine compute_kernel



end module calc_kernel
