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
  !   + the upper "left" corner of the fault is at (0,0,Z_CORNER)
  !   + z is positive downwards (depth)
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

! along-dip faster
!    !------ give value to x, y, z , dip of first row -------------------
!    cd = dcos(pb%mesh%DIP_W(1)/180d0*PI)
!    sd = dsin(pb%mesh%DIP_W(1)/180d0*PI)
!    do j = 1,pb%mesh%nx
!      pb%mesh%x((j-1)*pb%mesh%nw+1) = 0d0+(0.5d0+dble(j-1))*pb%mesh%dx
!      pb%mesh%y((j-1)*pb%mesh%nw+1) = 0d0+0.5d0*pb%mesh%dw(1)*cd
!      pb%mesh%z((j-1)*pb%mesh%nw+1) = pb%mesh%Z_CORNER+0.5d0*pb%mesh%dw(1)*sd
!      pb%mesh%dip((j-1)*pb%mesh%nw+1) = pb%mesh%DIP_W(1)
!    end do
!    !------ give value to x, y, z , dip of first row -------------------
!
!    !------ give value to x, y, z , dip of row 2 to nw------------------- 
!    do i = 2,pb%mesh%nw
!      cd0 = dcos(pb%mesh%DIP_W(i-1)/180d0*PI)
!      sd0 = dsin(pb%mesh%DIP_W(i-1)/180d0*PI)
!      cd = dcos(pb%mesh%DIP_W(i)/180d0*PI)
!      sd = dsin(pb%mesh%DIP_W(i)/180d0*PI)
!      do j = 1,pb%mesh%nx
!        pb%mesh%x((j-1)*pb%mesh%nw+i) = 0d0+(0.5d0+dble(j-1))*pb%mesh%dx
!        pb%mesh%y((j-1)*pb%mesh%nw+i) = pb%mesh%y((j-1)*pb%mesh%nw+i-1)    &
!                          +0.5d0*pb%mesh%dw(i-1)*cd0+0.5d0*pb%mesh%dw(i)*cd
!        pb%mesh%z((j-1)*pb%mesh%nw+i) = pb%mesh%z((j-1)*pb%mesh%nw+i-1)    &
!                          +0.5d0*pb%mesh%dw(i-1)*sd0+0.5d0*pb%mesh%dw(i)*sd
!        pb%mesh%dip((j-1)*pb%mesh%nw+i) = pb%mesh%DIP_W(i)
!      end do
!    end do
!    !------ give value to x, y, z , dip of row 2 to nw------------------- 
! along-dip faster

!    write(6,*) 'x,y,z'
!    do i = 1,pb%mesh%nn
!      write(6,*) pb%mesh%x(i),pb%mesh%y(i),pb%mesh%z(i)
!    end do 
    write(6,*) 'dw'
    write(6,*) pb%mesh%dw

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

  type(problem_type), intent(inout) :: pb

  double precision :: tau_co, wl2, tau
  integer :: i, j, k, i_src, i_obs, IRET

  write(6,*) 'Intializing kernel: ...'

  if (pb%kernel%kind == 1) then      ! 1D
    write(6,*) 'Single degree-of-freedom spring-block system'
    pb%kernel%k1 = pb%smu/pb%mesh%Lfault

  elseif (pb%kernel%kind == 2) then      ! 2D
    pb%kernel%k2f%nnfft = (pb%kernel%k2f%finite+1)*pb%mesh%nn 
    allocate (pb%kernel%k2f%kernel(pb%kernel%k2f%nnfft))

    write(6,*) 'OouraFFT Selected'

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
  !------ calculate kernel -----------------
  !kernel(i,j): response at i of source at j
  !because dx = constant, only need to calculate i at first column

    write(6,*) 'Generating 3D kernel...'
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

!JPA : version 2, does not need to be implemented     

!    do k = 1,pb%mesh%nx-1
!      do n = 1,pb%mesh%nw
!        do j = 1,pb%mesh%nw
!          call compute_kernel(pb%lam,pb%smu,pb%mesh%x(i_src),pb%mesh%y(i_src),pb%mesh%z(i_src),  &
!                 pb%mesh%dip(i_src),pb%mesh%dx,pb%mesh%dw((j-1)/pb%mesh%nx+1),   &
!                 pb%mesh%x(i_obs),pb%mesh%y(i_obs),   &
!                 pb%mesh%z(i_obs),pb%mesh%dip(i_obs),IRET,tau)
!          if (IRET == 0) then
!            pb%kernel%k3%kernel(j,n,k) = tau    
!         else
!            write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
!            pb%kernel%k3%kernel(j,n,k) = 0d0
!          end if
!        end do
!      end do
!    end do

!YD : call compute_kernel change made above:
!     change also required in : mesh x,y,z,dip : [change made and commented]
!                               calling algorithm in compute_stress_3D [change not made yet]
!JPA      
! note: requires a storage of x,y,z in which the along-dip index runs faster than the along-strike index
!    allocate(pb%kernel%k3%kernel(nw,nw,0:nx-1))
!    do k=0,nx-1 
!      do n=1,nw
!        do j=1,nw
!          call compute_kernel( ... j ... k*nw+n ..., IRET,tau)
!          pb%kernel%k3%kernel(j,n,k) = tau
!        enddo
!      enddo
!    enddo


! JPA : version 3, with FFT along-strike. This is the version we should implement 
!       Assumes faster index runs along-strike
!
! k3f%nxfft = 2*nx ! fft convolution requires twice longer array
! allocate(k3f%kernel(nw,nw,k3f%nxfft))
! do n=1,nw
!   nn = (n-1)*pb%mesh%nx
!   y_src = pb%mesh%y(nn)
!   z_src = pb%mesh%z(nn)
!   dip_src = pb%mesh%dip(nn)
!   dw_src = pb%mesh%dw(nn)
!   do j=1,nw
!     jj = (j-1)*pb%mesh%nx
!     y_obs = pb%mesh%y(jj)
!     z_obs = pb%mesh%z(jj)
!     dip_obs = pb%mesh%dip(jj)
!     do i=-nx+1,nx
!       call compute_kernel(pb%lam,pb%smu, &
!               i*pb%mesh%dx, y_src, z_src, dip_src, pb%mesh%dx, dw_src,   &
!               0d0, y_obs, z_obs, dip_obs, &
!               IRET,tau)
!       k = i+1
!       if (i<0) k = k+k3f%nxfft  ! wrap up the negative relative-x-positions in the second half of the array
!                                 ! to comply with conventions of fft convolution
!       tmp(k) = tau
!     enddo
!     call my_rdft(1,tmp,k3f%m_fft)
!     k3f%kernel(j,n,:) = tmp
!   enddo
! enddo


    write(6,*) 'kernel(1,j)'
    write(6,*) pb%kernel%k3%kernel(1,:)
    do i = 1,pb%mesh%nn
      write(99,*) pb%kernel%k3%kernel(1,i)
    end do
  end if

  write(6,*) 'Kernel intialized'
  
end subroutine init_kernel     

end module initialize
