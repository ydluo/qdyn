!qd intialize all include parameters, kernel  and fields

module initialize
  
 implicit none
 private  i,tau_co,wl2,nxout_count
   double precision:: tau_co,wl2
   integer :: i, nxout_count

 public   init_field, init_kernel 


contains

!=============================================================

  subroutine init_field(pb)
  
  use problem_class
  use constants, only : PI

  type(problem_type) :: pb
  
  write(6,*) 'Intializing Parameters: ...'
  
  if (pb%mesh%kind == 0) then  ! 1D fault, uniform grid
     write(6,*) '1D fault, uniform grid' 
     pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nn

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
     do i = 1,pb%mesh%nn
       pb%theta_star(i) = pb%dc(i)/pb%v2(i)
       pb%tau_init(i) = pb%sigma(i)*(pb%mu_star(i)-  &
             pb%a(i)*log(pb%v1(i)/pb%v(i)+1.d0)+ &  
             pb%b(i)*log(pb%theta(i)/pb%theta_star(i)+1.d0))
       pb%tau(i) = pb%tau_init(i)
       pb%slip(i) = 0.d0 
     end do
     !---------------------- ref_value -----------------

     !---------------------- output field to screen -----------------
     write(6,*) '**Field Initialized.**'
     write(6,*) 'Values at center of the fault:'
     i = pb%mesh%nn/2
     if (i == 0) i = 1
     if (pb%kernel%k2f%finite == 1 .or. pb%mesh%nn == 1) then
        write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%Lfault)/(pb%sigma(i)*(pb%b(i)-pb%a(i))/pb%dc(i))
        write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%Lfault)/(pb%sigma(i)*pb%b(i)/pb%dc(i))
      else
        write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%W)/(pb%sigma(i)*(pb%b(i)-pb%a(i))/pb%dc(i))
        write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%W)/(pb%sigma(i)*pb%b(i)/pb%dc(i))
     endif
     write(6,*)

     !---------------------- output field to screen -----------------

     !WRITE_ot--- write ot file head START-----------------------------
      write(18,'(a)')'# macroscopic values:'
      write(18,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
      write(18,'(a)')'# values at center:'
      write(18,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
      write(18,'(a)')'# values at max(V) location:'
      write(18,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'
     !WRITE_ot--- write ot file head END-------------------------------
     !WRITE_ox--- write ox file head START-----------------------------
      nxout_count = 0
      do i=1,pb%mesh%nn,pb%output%nxout
        nxout_count=nxout_count+1
      enddo
      write(19,'(a,2i5)')'# nx= ',nxout_count
     !WRITE_ox--- write ox file head START-----------------------------

   end if


   end subroutine init_field

  
 
   subroutine init_kernel(pb)
   
   use problem_class
   type(problem_type) :: pb


     !======================== kernel =============================
   if (pb%mesh%kind == 0) then      ! 1D
     pb%kernel%k2f%nnfft = (pb%kernel%k2f%finite+1)*pb%mesh%nn 
     if (pb%mesh%nn == 1) then      ! single degree-of-freedom spring-block system
       write(6,*) 'Single degree-of-freedom spring-block system'
       pb%kernel%k2f%kernel(1) = pb%smu/pb%mesh%Lfault
     end if
     if (pb%mesh%nn > 1 .and. pb%kernel%kind == 0) then
       write(6,*) 'OouraFFT Selected'

       pb%kernel%k2f%m_fft%nwfft = 2+ceiling(dsqrt(dble(pb%kernel%k2f%nnfft/2.d0)))
       allocate (pb%kernel%k2f%m_fft%iworkfft(0:pb%kernel%k2f%m_fft%nwfft-1))
       allocate (pb%kernel%k2f%m_fft%rworkfft(0:pb%kernel%k2f%nnfft/2-1))
  
       if (pb%kernel%k2f%finite == 0) then
         pb%kernel%k2f%nnfft = pb%mesh%nn
         allocate (pb%kernel%k2f%kernel(pb%kernel%k2f%nnfft))
         tau_co = PI*pb%smu/pb%mesh%Lfault *2.d0/dble(pb%mesh%nn)
         wl2 = (pb%mesh%Lfault/pb%mesh%W)**2
            do i=0,pb%mesh%nn/2-1
            pb%kernel%k2f%kernel(2*i+1) = tau_co*dsqrt(dble(i*i)+wl2)
            pb%kernel%k2f%kernel(2*i+2) = pb%kernel%k2f%kernel(2*i+1)
            enddo
            pb%kernel%k2f%kernel(2) = tau_co*dsqrt(dble(pb%mesh%nn**2)/4.d0+wl2) ! Nyquist
       end if
       
       if (pb%kernel%k2f%finite == 1) then
         pb%kernel%k2f%nnfft=pb%mesh%nn*2
         !- Read coefficient I(n) from pre-calculated file.
         open(57,file='~/2D_RUPTURE/STATIC/Matlab/kernel_I_32768.tab')
           if (pb%kernel%k2f%nnfft/2>32768) stop 'Finite kernel table is too small'
             do i=1,pb%kernel%k2f%nnfft/2-1
               read(57,*) pb%kernel%k2f%kernel(2*i+1)
             enddo
           read(57,*) pb%kernel%k2f%kernel(2) ! Nyquist
         close(57)
         ! The factor 2/N comes from the inverse FFT convention
         tau_co = PI*pb%smu / (2d0*pb%mesh%Lfault) *2.d0/dble(pb%kernel%k2f%nnfft)
         pb%kernel%k2f%kernel(1) = 0d0
         pb%kernel%k2f%kernel(2) = tau_co*dble(pb%kernel%k2f%nnfft/2)*pb%kernel%k2f%kernel(2)
         do i = 1,pb%kernel%k2f%nnfft/2-1
           pb%kernel%k2f%kernel(2*i+1) = tau_co*dble(i)*pb%kernel%k2f%kernel(2*i+1)
           pb%kernel%k2f%kernel(2*i+2) = pb%kernel%k2f%kernel(2*i+1)
         enddo
       end if
     end if

    end if
     !======================== kernel =============================

  
  end subroutine init_kernel     

end module initialize
