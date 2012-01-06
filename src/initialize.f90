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
  
  integer :: i

  write(6,*) 'Initializing parameters: ...'
  
  if (pb%mesh%kind == 0) then  ! 1D fault, uniform grid

    write(6,*) '1D fault, uniform grid' 
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nn
    do i=1,pb%mesh%nn
      pb%mesh%x(i) = (i-pb%mesh%nn*0.5d0-0.5d0)*pb%mesh%dx
    enddo

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

  end if

end subroutine init_field

  
!=============================================================
 
subroutine init_kernel(pb)
  
  use constants, only : PI 
  use problem_class

  type(problem_type), intent(inout) :: pb

  double precision :: tau_co, wl2
  integer :: i

  write(6,*) 'Intializing kernel: ...'

  if (pb%mesh%kind == 0) then      ! 1D

    pb%kernel%k2f%nnfft = (pb%kernel%k2f%finite+1)*pb%mesh%nn 
    allocate (pb%kernel%k2f%kernel(pb%kernel%k2f%nnfft))

    if (pb%mesh%nn == 1) then      ! single degree-of-freedom spring-block system
      write(6,*) 'Single degree-of-freedom spring-block system'
      pb%kernel%k2f%kernel(1) = pb%smu/pb%mesh%Lfault

    elseif (pb%mesh%nn > 1 .and. pb%kernel%kind == 0) then
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

    end if

    write(6,*) 'Kernel intialized'
  end if
  
end subroutine init_kernel     

end module initialize
