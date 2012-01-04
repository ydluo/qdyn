! Solve_master

module solve_master

  implicit none
 
  private

  public   solve 
  

contains
  
subroutine solve(pb)
  
  use problem_class
  use ode_bs
  use derivs_all
  use constants, only : YEAR
  
  type(problem_type), intent(inout)  :: pb

  integer :: it,ic,i,ip,ivmax,nxout_count
  double precision :: dt_next, dt_did,  &
                      pot, pot_rate, xvmax, vmax, vmax_old, vmax_older, &
                      xloc, omega, dtau_per
  double precision, dimension(:), allocatable ::  yt, dydt, yt_scale
  
  !=======================Time loop. START============================
  write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'

  pb%time = 0.d0
  pb%itstop = -1
  vmax_old = 0d0
  vmax_older = 0d0

 
  allocate (yt(pb%neqs*pb%mesh%nn))
  allocate (dydt(pb%neqs*pb%mesh%nn))
  allocate (yt_scale(pb%neqs*pb%mesh%nn))
  yt(1:pb%neqs:) = pb%v
  yt(2:pb%neqs:) = pb%theta
      
  ! Time loop

  it=0
  do while (it /= pb%itstop)

    it=it+1

    call derivs(pb)
        
    yt(1:pb%neqs:) = pb%v
    yt(2:pb%neqs:) = pb%theta
    dydt(1:pb%neqs:) = pb%dv_dt
    dydt(2:pb%neqs:) = pb%dtheta_dt
   
    ! One step
 
    !--------Call EXT routine bsstep [Bulirsch-Stoer Method] --------------
    !-------- 
    yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
    call bsstep(yt,dydt,pb%neqs*pb%mesh%nn,pb%time,pb%dt_try,pb%acc,yt_scale,   &
                dt_did,dt_next,derivs,pb)
    if (pb%dt_max >  0.d0) then
      pb%dt_try=min(dt_next,pb%dt_max)
    else
      pb%dt_try=dt_next
    endif

    pb%v = yt(1:pb%neqs:)
    pb%theta = yt(2:pb%neqs:)
    pb%dv_dt = dydt(1:pb%neqs:)
    pb%dtheta_dt = dydt(2:pb%neqs:) 

    ! Update slip, stress. 
    pb%slip = pb%slip + pb%v*dt_did
    pb%tau  = ( pb%mu_star - pb%a*log(pb%v1/pb%v+1d0) +pb%b*log(pb%theta/pb%theta_star)+1d0 ) * pb%sigma
    
    ! update potency and potency rate
    pot = sum(pb%slip) * pb%mesh%dx
    pot_rate = sum(pb%v) * pb%mesh%dx

    ! update crack size
    pb%ot%lcold = pb%ot%lcnew
    pb%ot%lcnew = crack_size(pb%slip,pb%mesh%nn)

    pb%ot%llocold = pb%ot%llocnew
    pb%ot%llocnew = crack_size(pb%dtau_dt,pb%mesh%nn)

    ! Output time series at max(v) location
    ivmax = maxloc(pb%v)
    vmax = pb%v(ivmax)
    xvmax = pb%mesh%x(ivmax)

    call ot_write()

    if (pb%itstop == 0) then
      !         STOP soon after end of slip localization 
      if (pb%NSTOP == 1) then
        if (pb%output%llocnew > pb%output%llocold) pb%itstop=it+2*pb%output%ntout

      ! STOP soon after maximum slip rate
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

    !WRITE_6---------------- Print step to screen START--------------------
    if(mod(it-1,pb%output%ntout) == 0 .or. it == pb%itstop) then
    ! Print info
      write(6,'(i7,x,3(e11.3,x),i5)') it, dt_did, pb%time/YEAR, vmax
    !WRITE_6---------------- Print step to screen END----------------------21--------------: Command not found.
      call ox_write()

    endif
  enddo

end subroutine solve



end module solve_master

