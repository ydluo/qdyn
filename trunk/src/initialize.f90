! initialize all include parameters, kernel and fields

module initialize
  
  implicit none
  private

  public :: init_field, init_kernel 

contains

!=============================================================

subroutine init_field(pb)
  
  use problem_class
  use constants, only : PI
  use output, only : ot_init, ox_init

  type(problem_type), intent(inout) :: pb
  
  integer :: i, j, j0
  double precision :: cd, sd, cd0, sd0

  write(6,*) 'Initializing parameters: ...'
  
  select case (pb%mesh%dim)

  ! Spring-block System
  case(0) 
    write(6,*) 'Spring-block System' 
    pb%mesh%dx = pb%mesh%Lfault
    pb%mesh%x = 0d0

  ! 1D fault, uniform grid
  case(1)
    write(6,*) '1D fault, uniform grid' 
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nn
    do i=1,pb%mesh%nn
      pb%mesh%x(i) = (i-pb%mesh%nn*0.5d0-0.5d0)*pb%mesh%dx
      ! Assuming nn is even (usually a power of 2), 
      ! the center of the two middle elements (i=nn/2 and nn/2+1) 
      ! are located at x=-dx/2 and x=dx/2, respectively 
    enddo

  ! 2D fault, uniform grid along-strike
  ! Assumptions: 
  !   + the fault trace is parallel to x 
  !   + the lower "left" corner of the fault is at (0,0,Z_CORNER)
  !   + z is negative downwards (-depth)
  ! Storage scheme: faster index runs along-strike (x)
  case(2)

    write(6,*) '2D fault, uniform grid along-strike'
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nx
    allocate(pb%mesh%y(pb%mesh%nn),   &
             pb%mesh%z(pb%mesh%nn),   &
             pb%mesh%dip(pb%mesh%nn)) 

    ! set x, y, z, dip of first row
    cd = dcos(pb%mesh%DIP_W(1)/180d0*PI)
    sd = dsin(pb%mesh%DIP_W(1)/180d0*PI)
    do j = 1,pb%mesh%nx
      pb%mesh%x(j) = 0d0+(0.5d0+dble(j-1))*pb%mesh%dx
    end do
    pb%mesh%y(1:pb%mesh%nx) = 0d0+0.5d0*pb%mesh%dw(1)*cd
    pb%mesh%z(1:pb%mesh%nx) = pb%mesh%Z_CORNER+0.5d0*pb%mesh%dw(1)*sd
    pb%mesh%dip(1:pb%mesh%nx) = pb%mesh%DIP_W(1)

    ! set x, y, z, dip of row 2 to nw
    do i = 2,pb%mesh%nw
      cd0 = cd
      sd0 = sd
      cd = dcos(pb%mesh%DIP_W(i)/180d0*PI)
      sd = dsin(pb%mesh%DIP_W(i)/180d0*PI)
      j0 = (i-1)*pb%mesh%nx
      pb%mesh%x(j0+1:j0+pb%mesh%nx) = pb%mesh%x(1:pb%mesh%nx)
      pb%mesh%y(j0+1:j0+pb%mesh%nx) = pb%mesh%y(j0) + 0.5d0*pb%mesh%dw(i-1)*cd0 + 0.5d0*pb%mesh%dw(i)*cd
      pb%mesh%z(j0+1:j0+pb%mesh%nx) = pb%mesh%z(j0) + 0.5d0*pb%mesh%dw(i-1)*sd0 + 0.5d0*pb%mesh%dw(i)*sd
      pb%mesh%dip(j0+1:j0+pb%mesh%nx) = pb%mesh%DIP_W(i)
    end do

  case(3)

    write(6,*) '2D fault, uniform grid along-strike'
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nx
    allocate(pb%mesh%y(pb%mesh%nn),   &
             pb%mesh%z(pb%mesh%nn),   &
             pb%mesh%dip(pb%mesh%nn)) 

    ! set x, y, z, dip of first row
    cd = dcos(pb%mesh%DIP_W(1)/180d0*PI)
    sd = dsin(pb%mesh%DIP_W(1)/180d0*PI)
    do j = 1,pb%mesh%nx
      pb%mesh%x(j) = 0d0+(0.5d0+dble(j-1))*pb%mesh%dx
    end do
    pb%mesh%y(1:pb%mesh%nx) = 0d0+0.5d0*pb%mesh%dw(1)*cd
    pb%mesh%z(1:pb%mesh%nx) = pb%mesh%Z_CORNER+0.5d0*pb%mesh%dw(1)*sd
    pb%mesh%dip(1:pb%mesh%nx) = pb%mesh%DIP_W(1)

    ! set x, y, z, dip of row 2 to nw
    do i = 2,pb%mesh%nw
      cd0 = cd
      sd0 = sd
      cd = dcos(pb%mesh%DIP_W(i)/180d0*PI)
      sd = dsin(pb%mesh%DIP_W(i)/180d0*PI)
      j0 = (i-1)*pb%mesh%nx
      pb%mesh%x(j0+1:j0+pb%mesh%nx) = pb%mesh%x(1:pb%mesh%nx)
      pb%mesh%y(j0+1:j0+pb%mesh%nx) = pb%mesh%y(j0) + 0.5d0*pb%mesh%dw(i-1)*cd0 + 0.5d0*pb%mesh%dw(i)*cd
      pb%mesh%z(j0+1:j0+pb%mesh%nx) = pb%mesh%z(j0) + 0.5d0*pb%mesh%dw(i-1)*sd0 + 0.5d0*pb%mesh%dw(i)*sd
      pb%mesh%dip(j0+1:j0+pb%mesh%nx) = pb%mesh%DIP_W(i)
    end do


  end select

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs
 
    !---------------------- dt_max & perturbation------------------
  if (pb%Aper /= 0.d0 .and. pb%Tper > 0.d0) then
    if (pb%dt_max > 0.d0) then
      pb%dt_max = min(pb%dt_max,0.2d0*pb%Tper)
    else
      pb%dt_max = 0.2d0*pb%Tper
    endif
  endif
  if (pb%Tper > 0.d0) then
    pb%Omper = 2.d0*PI/pb%Tper
  else
    pb%Omper = 0.d0
  endif
  !---------------------- dt_max & perturbation------------------

  !---------------------- impedance ------------------
  if (pb%beta > 0d0) then 
    pb%zimpedance = 0.5d0*pb%smu/pb%beta
  else
    pb%zimpedance = 0.d0
  endif
  write(6,*) 'impedance = ', pb%zimpedance
  !---------------------- impedance ------------------
     
  !---------------------- ref_value ------------------        
  pb%theta_star = pb%dc/pb%v2
  pb%tau_init = pb%sigma *   &
    (pb%mu_star- pb%a*log(pb%v1/pb%v+1d0)+ pb%b*log(pb%theta/pb%theta_star+1d0))
  pb%tau = pb%tau_init
  pb%slip = 0d0 
  !---------------------- ref_value ----------------- 

  !---------------------- init_value for solver ----------------- 
  pb%time = 0.d0
  pb%itstop = -1
  pb%it = 0
  !---------------------- init_value for solver ----------------- 
   
  call ot_init(pb)
  call ox_init(pb)
  write(6,*) 'Field Initialized'


end subroutine init_field

  
!=============================================================
 
subroutine init_kernel(pb)
  
  use constants, only : PI 
  use problem_class
  use okada, only : compute_kernel
  use fftsg, only : my_rdft

  type(problem_type), intent(inout) :: pb

  double precision :: tau_co, wl2, tau
  double precision :: y_src, z_src, dip_src, dw_src 
  double precision :: y_obs, z_obs, dip_obs
  double precision, allocatable :: tmp(:)   ! for FFT
  integer :: i, j, jj, k, n, nn, IRET

  write(6,*) 'Intializing kernel: ...'

  if (pb%kernel%kind == 1) then      ! 1D
    write(6,*) 'Single degree-of-freedom spring-block system'
    pb%kernel%k1 = pb%smu/pb%mesh%Lfault

  elseif (pb%kernel%kind == 2) then      ! 2D
    pb%kernel%k2f%nnfft = (pb%kernel%k2f%finite+1)*pb%mesh%nn 
    allocate (pb%kernel%k2f%kernel(pb%kernel%k2f%nnfft))

    write(6,*) 'OouraFFT Applied'

    if (pb%kernel%k2f%finite == 0) then
      tau_co = PI*pb%smu/pb%mesh%Lfault *2.d0/pb%mesh%nn
      wl2 = (pb%mesh%Lfault/pb%mesh%W)**2
      do i=0,pb%mesh%nn/2-1
        pb%kernel%k2f%kernel(2*i+1) = tau_co*sqrt(i*i+wl2)
        pb%kernel%k2f%kernel(2*i+2) = pb%kernel%k2f%kernel(2*i+1)
      enddo
      pb%kernel%k2f%kernel(2) = tau_co*sqrt(pb%mesh%nn**2/4.d0+wl2) ! Nyquist
       
    elseif (pb%kernel%k2f%finite == 1) then
      !- Read coefficient I(n) from pre-calculated file.
      open(57,file='~/2D_RUPTURE/STATIC/Matlab/kernel_I_32768.tab')
      if (pb%kernel%k2f%nnfft/2>32768) stop 'Finite kernel table is too small'
      do i=1,pb%kernel%k2f%nnfft/2-1
        read(57,*) pb%kernel%k2f%kernel(2*i+1)
      enddo
      read(57,*) pb%kernel%k2f%kernel(2) ! Nyquist
      close(57)
      ! The factor 2/N comes from the inverse FFT convention
      tau_co = PI*pb%smu / (2d0*pb%mesh%Lfault) *2.d0/pb%kernel%k2f%nnfft
      pb%kernel%k2f%kernel(1) = 0d0
      pb%kernel%k2f%kernel(2) = tau_co*dble(pb%kernel%k2f%nnfft/2)*pb%kernel%k2f%kernel(2)
      do i = 1,pb%kernel%k2f%nnfft/2-1
        pb%kernel%k2f%kernel(2*i+1) = tau_co*dble(i)*pb%kernel%k2f%kernel(2*i+1)
        pb%kernel%k2f%kernel(2*i+2) = pb%kernel%k2f%kernel(2*i+1)
      enddo
    end if

  elseif (pb%kernel%kind == 3) then      ! 3D
    write(6,*) 'Generating 3D kernel...'
    write(6,*) 'OouraFFT applied along-strike'
  !------ calculate kernel -----------------
  !kernel(i,j): response at i of source at j
  !because dx = constant, only need to calculate i at first column


!    allocate (pb%kernel%k3%kernel(pb%mesh%nw,pb%mesh%nn))
!    do i = 1,pb%mesh%nw
!      do j = 1,pb%mesh%nn
!        call compute_kernel(pb%lam,pb%smu,pb%mesh%x(j),pb%mesh%y(j),pb%mesh%z(j),  &
!               pb%mesh%dip(j),pb%mesh%dx,pb%mesh%dw((j-1)/pb%mesh%nx+1),   &
!               pb%mesh%x(1+(i-1)*pb%mesh%nx),pb%mesh%y(1+(i-1)*pb%mesh%nx),   &
!               pb%mesh%z(1+(i-1)*pb%mesh%nx),pb%mesh%dip(1+(i-1)*pb%mesh%nx),IRET,tau)
!        write(6,*) 'obs',pb%mesh%x(1+(i-1)*pb%mesh%nx),pb%mesh%y(1+(i-1)*pb%mesh%nx),pb%mesh%z(1+(i-1)*pb%mesh%nx)
!        write(6,*) 'src',pb%mesh%x(j),pb%mesh%y(j),pb%mesh%z(j) 
!        write(6,*) 'tau',tau,'IRET',IRET
!        if (IRET == 0) then
!          pb%kernel%k3%kernel(i,j) = tau    
!        else
!          write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
!          pb%kernel%k3%kernel(i,j) = 0d0
!        end if
!      end do
!    end do


! version with FFT along-strike. This is the version we should implement 
! Assumes faster index runs along-strike
!
  pb%kernel%k3f%nw = pb%mesh%nw
  pb%kernel%k3f%nx = pb%mesh%nx
  pb%kernel%k3f%nxfft = 2*pb%mesh%nx ! fft convolution requires twice longer array
  allocate(pb%kernel%k3f%kernel(pb%mesh%nw,pb%mesh%nw,pb%kernel%k3f%nxfft))
  allocate(tmp(pb%kernel%k3f%nxfft))
  do n=1,pb%mesh%nw
    nn = (n-1)*pb%mesh%nx+1
    y_src = pb%mesh%y(nn)
    z_src = pb%mesh%z(nn)
    dip_src = pb%mesh%dip(nn)
    dw_src = pb%mesh%dw(n)
    do j=1,pb%mesh%nw
      jj = (j-1)*pb%mesh%nx+1
      y_obs = pb%mesh%y(jj)
      z_obs = pb%mesh%z(jj)
      dip_obs = pb%mesh%dip(jj)
      do i=-pb%mesh%nx+1,pb%mesh%nx
        call compute_kernel(pb%lam,pb%smu, &
                i*pb%mesh%dx, y_src, z_src, dip_src, pb%mesh%dx, dw_src,   &
                0d0, y_obs, z_obs, dip_obs, &
                IRET,tau)
        k = i+1
        if (i<0) k = k+pb%kernel%k3f%nxfft  ! wrap up the negative relative-x-positions in the second half of the array
                                  ! to comply with conventions of fft convolution
        tmp(k) = tau
      enddo
!      write(6,*) 'kernel_before_fft',tmp
      call my_rdft(1,tmp,pb%kernel%k3f%m_fft)
      pb%kernel%k3f%kernel(j,n,:) = tmp / dble(pb%mesh%nx)
!      write(6,*) 'kernel_after_fft',tmp
    enddo
  enddo


!    write(6,*) 'kernel(j,n,nxfft)'
!    write(6,*) pb%kernel%k3f%kernel
!    do i = 1,pb%mesh%nn
!      write(99,*) pb%kernel%k3%kernel(1,i)
!    end do


  elseif (pb%kernel%kind == 4) then      ! 3D no fft
    write(6,*) 'Generating 3D kernel...'
    write(6,*) 'NO FFT applied'
  !------ calculate kernel -----------------
  !kernel(i,j): response at i of source at j
  !because dx = constant, only need to calculate i at first column


    allocate (pb%kernel%k3%kernel(pb%mesh%nw,pb%mesh%nn))
    do i = 1,pb%mesh%nw
      do j = 1,pb%mesh%nn
        call compute_kernel(pb%lam,pb%smu,pb%mesh%x(j),pb%mesh%y(j),pb%mesh%z(j),  &
               pb%mesh%dip(j),pb%mesh%dx,pb%mesh%dw((j-1)/pb%mesh%nx+1),   &
               pb%mesh%x(1+(i-1)*pb%mesh%nx),pb%mesh%y(1+(i-1)*pb%mesh%nx),   &
               pb%mesh%z(1+(i-1)*pb%mesh%nx),pb%mesh%dip(1+(i-1)*pb%mesh%nx),IRET,tau)
!        write(6,*) 'obs',pb%mesh%x(1+(i-1)*pb%mesh%nx),pb%mesh%y(1+(i-1)*pb%mesh%nx),pb%mesh%z(1+(i-1)*pb%mesh%nx)
!        write(6,*) 'src',pb%mesh%x(j),pb%mesh%y(j),pb%mesh%z(j) 
!        write(6,*) 'tau',tau,'IRET',IRET
        if (IRET == 0) then
          pb%kernel%k3%kernel(i,j) = tau    
        else
          write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
          pb%kernel%k3%kernel(i,j) = 0d0
        end if
      end do
    end do
!    write(6,*) 'kernel',pb%kernel%k3%kernel
    do j = 1,pb%mesh%nn
      write(99,*) pb%kernel%k3%kernel(1,j)
    end do
  end if

  write(6,*) 'Kernel intialized'
  
end subroutine init_kernel     

end module initialize
