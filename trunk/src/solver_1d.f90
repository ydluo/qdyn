! Solver_1d

module solver_1d

  implicit none

  private    dtau_per, omega,    & 
             yt,dydt,yt_scale,           &
             ic,i,ip,ivmax,nxout_count,  &
             dt_next,dt_did,           &
             pot,pot_rate,   &
             vmax,xvmax,vmax_old,vmax_older,  &
             xloc

  public   solve_1d_fft, crack_size
  

  integer :: it,ic,i,ip,ivmax,nxout_count
  double precision :: dt_next, dt_did,  &
                      pot, pot_rate, xvmax, vmax, vmax_old, vmax_older, &
                      xloc, omega, dtau_per
  double precision, dimension(:), allocatable ::  yt, dydt, yt_scale
  
  contains

  subroutine solve_1d_fft(pb)

  use problem_class
  use ode_bs
  use fftsg
  use constants, only : YEAR
  
  type(problem_type) :: pb



  !WRITE_ot--- write ot file head START-----------------------------
      write(18,'(a)')'# macroscopic values:'
      write(18,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
      write(18,'(a)')'# values at center:'
      write(18,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
      write(18,'(a)')'# values at max(V) location:'
      write(18,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'
  !WRITE_ot--- write ot file head END-------------------------------
  !WRITE_ox--- write ox file head START-----------------------------
      nxout_count=0
      do i=1,pb%mesh%nn,pb%output%nxout
        nxout_count = nxout_count+1
      enddo
      write(19,'(a,2i5)')'# nx= ',nxout_count
  !WRITE_ox---1.5.2 write ox file head START-----------------------------

  !=======================Time loop. START============================
      write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'

      pb%time = 0.d0
      pb%itstop = -1
      pb%output%lcnew = dble(pb%mesh%nn)
      pb%output%llocnew = dble(pb%mesh%nn)
      vmax_old = 0d0
      vmax_older = 0d0

      ic = pb%mesh%nn/2
      if (ic == 0) ic = 1

      pb%kernel%k2f%m_fft%iworkfft(0) = 0
      allocate (yt(pb%neqs*pb%mesh%nn))
      allocate (dydt(pb%neqs*pb%mesh%nn))
      allocate (yt_scale(pb%neqs*pb%mesh%nn))
      do i = 1,pb%mesh%nn
        ip = pb%neqs*i
        yt(ip-1) = pb%v(i)
        yt(ip) = pb%theta(i)
      enddo
      
      
  ! Time loop

      it=0
      do while (it /= pb%itstop)

        it=it+1

        call derivs(pb)
      
        
        do i = 1,pb%mesh%nn
          ip = pb%neqs*i
          yt(ip-1) = pb%v(i)
          yt(ip) = pb%theta(i)
          dydt(ip-1) = pb%dv_dt(i)
          dydt(ip) = pb%dtheta_dt(i)
        enddo

   
    ! One step


 
    !--------Call EXT routine bsstep [Bulirsch-Stoer Method] --------------
    !-------- 
        do i=1,pb%neqs*pb%mesh%nn
          yt_scale(i)=dabs(yt(i))+dabs(pb%dt_try*dydt(i))
        enddo
        call bsstep(yt,dydt,pb%neqs*pb%mesh%nn,pb%time,pb%dt_try,pb%acc,yt_scale,   &
                    dt_did,dt_next,derivs)
        if (pb%dt_max >  0.d0) then
          pb%dt_try=min(dt_next,pb%dt_max)
        else
          pb%dt_try=dt_next
        endif
    !        theta is ip-1, slip rate is ip
        do i=1,pb%mesh%nn
          ip=pb%neqs*i
          pb%v(i) = yt(ip)
          pb%theta(i) = yt(ip-1)
          pb%dv_dt(i) = dydt(ip-1)
          pb%dtheta_dt(i) = dydt(ip) 
        enddo

    ! Update slip, stress. 
    ! potency and potency rate
    ! crack size
        pot = 0.d0
        pot_rate = 0.d0
        do i=1,pb%mesh%nn
          pb%slip(i) = pb%slip(i) + pb%v(i)*dt_did
          pot = pot + pb%slip(i)
          pot_rate = pot_rate + pb%v(i)
          pb%tau(i)  = ( pb%mu_star(i) - pb%a(i)*log(pb%v1(i)/pb%v(i)+1.d0)  &
            +pb%b(i)*log(pb%theta(i)/pb%theta_star(i))+1.d0 ) * pb%sigma(i)
        enddo
        pot = pot * pb%mesh%dx
        pot_rate = pot_rate * pb%mesh%dx

        pb%output%lcold = pb%output%lcnew
        pb%output%lcnew = crack_size(pb%slip,pb%mesh%nn)

        pb%output%llocold = pb%output%llocnew
        pb%output%llocnew = crack_size(pb%dtau_dt,pb%mesh%nn)

    ! Output time series at max(v) location
        vmax=0.d0
        do i=1,pb%mesh%nn
          if(pb%v(i) > vmax)then
            vmax = pb%v(i)
            ivmax=i
            xvmax=(i-pb%mesh%nn*0.5d0-0.5d0)*pb%mesh%dx
          endif
        enddo

    !WRITE_ot-------------- Export ot values START------------------------
        write(18,'(e24.16,15e14.6)')pb%time, pb%output%llocnew*pb%mesh%dx,  &
          pb%output%lcnew*pb%mesh%dx, pot, pot_rate,    &
          pb%v(ic), pb%theta(ic), pb%v(ic)*pb%theta(ic)/pb%dc(ic), &
          pb%tau(ic), pb%slip(ic),    &
          xvmax, vmax, pb%theta(ivmax), vmax*pb%theta(ivmax)/pb%dc(ivmax),    &
          pb%tau(ivmax), pb%slip(ivmax)
    !WRITE_ot-------------- Export ot  values END--------------------------

    !WRITE_121-------------- STOP criteria START---------------------------
        if (pb%itstop == 0) then
        !         STOP soon after end of slip localization 
          if (pb%NSTOP == 1) then
            if (pb%output%llocnew > pb%output%llocold) pb%itstop=it+2*pb%output%ntout

        !         STOP soon after maximum slip rate
          elseif (pb%NSTOP == 2) then

            if (it > 2 .and. vmax_old > vmax_older    &
               .and. vmax < vmax_old) pb%itstop = it+10*pb%output%ntout
            vmax_older = vmax_old
            vmax_old = vmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!! ???????????????????????????????
        !         STOP at a slip rate threshold
          elseif (pb%NSTOP == 3) then    
            if (vmax > pb%tmax) pb%itstop = it
!!! ??????????????????????????????
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !         STOP if time > tmax
          else
                  write(121,*) pb%time, pb%tmax
            if (pb%tmax > 0.d0 .and. pb%time > pb%tmax) pb%itstop = it
          endif
        endif
    !WRITE_121-------------- STOP criteria END-----------------------------

    !WRITE_6---------------- Print step to screen START--------------------
        if(mod(it-1,pb%output%ntout) == 0 .or. it == pb%itstop) then
     ! Print info
        write(6,'(i7,x,3(e11.3,x),i5)') it, dt_did, pb%time/YEAR, vmax
    !WRITE_6---------------- Print step to screen END----------------------

    !WRITE_ox--------------- Export snapshots START------------------------
        write(19,'(2a,2i5,e14.6)')'# x v theta',    &
              ' V./V dtau tau_dot slip ',it,ivmax,pb%time
          do i=1,pb%mesh%nn,pb%output%nxout
            ip=pb%neqs*i
            xloc=(i-pb%mesh%nn*0.5d0-0.5d0)*pb%mesh%dx
            write(19,'(8e15.7)')xloc,pb%time,pb%v(i), &
              pb%theta(i),dydt(ip)/pb%v(i),pb%tau(i),     &
              pb%dtau_dt(i),pb%slip(i)
          enddo
        endif
    !WRITE_ox--------------- Export snapshots END------------------------

      enddo
 end subroutine solve_1d_fft



 !C====================Subroutine derivs START===========================
 !-------Compute thata,vslip and time deriv for one time step------------
 !-------CALL rdft: real DFT---------------------------------------------
 !---------------------------------------------------------------|
 !-------|     in:pb                                     |-------|
 !-------|    out:pb                                     |-------|
 !-------|                                               |-------|
 !---------------------------------------------------------------|
 
   subroutine derivs(pb)
   
     
     use problem_class
     use fftsg
  
     type(problem_type) :: pb
     double precision :: dtau_per, omega


   ! compute shear stress rate from elastic interactions
     if (pb%mesh%nn > 1) then
       do i = 1,pb%mesh%nn
         pb%dtau_dt(i) = pb%v_star(i)-pb%v(i)
       enddo
       do i = pb%mesh%nn+1,pb%kernel%k2f%nnfft
         pb%dtau_dt(i) = 0d0 
       enddo
       call rdft(pb%kernel%k2f%nnfft,1,pb%dtau_dt,pb%kernel%k2f%m_fft%iworkfft,pb%kernel%k2f%m_fft%rworkfft)
         do i = 1,pb%kernel%k2f%nnfft
           pb%dtau_dt(i) = pb%kernel%k2f%kernel(i)*pb%dtau_dt(i)
         enddo
       call rdft(pb%kernel%k2f%nnfft,-1,pb%dtau_dt,pb%kernel%k2f%m_fft%iworkfft,pb%kernel%k2f%m_fft%rworkfft)
     else
       pb%dtau_dt(1) = pb%kernel%k2f%kernel(1)*( pb%v_star(1)-pb%v(1) )
     endif

   ! periodic loading
     dtau_per = pb%Omper * pb%Aper * dcos(pb%Omper*pb%time)     


     do i=1,pb%mesh%nn

   !--------State evolution law START--------------------------------
      omega = pb%v(i)*pb%theta(i) / pb%dc(i)

   !      "aging" law
        if (pb%itheta_law == 1) then
          pb%dtheta_dt(i) = 1.d0-omega

   !      "slip" law
        elseif (pb%itheta_law== 2) then
          if (omega /= 0d0) then
            pb%dtheta_dt(i) = -omega*dlog(omega)
          else
            pb%dtheta_dt(i) = 0d0
          endif

   !      "aging" in the no-healing approximation
        elseif (pb%itheta_law == 0) then
          pb%dtheta_dt(i) = -omega
        endif
   !--------State evolution law END----------------------------------


        pb%dv_dt(i) = (dtau_per+pb%dtau_dt(i) - pb%sigma(i)*pb%b(i)*pb%v2(i)*pb%dtheta_dt(i)/(pb%v2(i)*pb%theta(i)+pb%dc(i)) )    &
                  /( pb%sigma(i)*pb%a(i)*(1.d0/pb%v(i)-1.d0/(pb%v1(i)+pb%v(i))) + pb%zimpedance )
    
    enddo

      return
   end subroutine derivs
   !C====================Subroutine derivs END=============================




   !======================Function crack_size START========================
   ! distance between largest peak on the left half
   ! and largest peak on the right half
     function crack_size(s,n)
       
       
       integer n
       double precision ::  s(n)
       double precision ::  smin,smax,s1,s2,xL,xR
       integer i,iL,iR,imin
     double precision crack_size
     if (n == 0) then
       crack_size = 0d0
       return
     endif

   ! assuming minimum stressing rate at nucleation point
   !      smin=s(1)
   !      do i=1,n
   !        if (s(i).le.smin) then
   !          imin=i
   !          smin=s(i)
   !        endif
   !      enddo

   ! WARNING: assumes one peak on x<0 and one peak on x>0
     imin=n/2

     iL = 1
     smax=s(iL)
     do i=1,imin
       if (s(i) >= smax) then
         iL=i
         smax=s(i)
       endif
     enddo
     xL = dble(iL);
     if (iL > 1) then
       s1 = 0.5d0*( s(iL+1)-s(iL-1) );
       s2 = s(iL-1)-2d0*s(iL)+s(iL+1);
       if (s2 /= 0d0) xL = xL-s1/s2;
     endif
      
     iR = imin
     smax=s(iR)
     do i=imin,n
       if (s(i) >= smax) then
         iR=i
         smax=s(i)
       endif
     enddo
     xR = dble(iR);
     if (iR < n) then
       s1 = 0.5d0*( s(iR+1)-s(iR-1) );
       s2 = s(iR-1)-2d0*s(iR)+s(iR+1);
       if (s2 /= 0d0) xR = xR-s1/s2;
     endif

     crack_size = xR-xL

     end function crack_size
   !======================Function crack_size END==========================




end module solver_1d
