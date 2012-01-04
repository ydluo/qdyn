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
  use output
  use constants, only : YEAR
  
  type(problem_type), intent(inout)  :: pb

  integer :: it
  double precision :: dt_next, dt_did,  &
                      pot, pot_rate, xvmax, vmax, vmax_old, vmax_older, &
                      xloc, omega, dtau_per
  double precision, dimension(:), allocatable ::  yt, dydt, yt_scale
  
  !=======================Time loop. START============================
  write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'


  !--------Allocate working space for yt... and init--------------------- 
  allocate (yt(pb%neqs*pb%mesh%nn))
  allocate (dydt(pb%neqs*pb%mesh%nn))
  allocate (yt_scale(pb%neqs*pb%mesh%nn))
  yt(1:pb%neqs:) = pb%v
  yt(2:pb%neqs:) = pb%theta
  !--------Allocate working space for yt... and init--------------------- 



     
  ! Time loop

  it=0
  do while (it /= pb%itstop)

    it=it+1

    call derivs(pb)

    !-------Pack v, theta into yt---------------------------------    
    yt(1:pb%neqs:) = pb%v
    yt(2:pb%neqs:) = pb%theta
    dydt(1:pb%neqs:) = pb%dv_dt
    dydt(2:pb%neqs:) = pb%dtheta_dt
    !-------Pack v, theta into yt--------------------------------- 

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

    !-------Unpack yt into v, theta--------------------------------- 
    pb%v = yt(1:pb%neqs:)
    pb%theta = yt(2:pb%neqs:)
    pb%dv_dt = dydt(1:pb%neqs:)
    pb%dtheta_dt = dydt(2:pb%neqs:) 
    !-------Unpack yt into v, theta--------------------------------- 


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
    pb%ot%ivmax = maxloc(pb%v)

    call ot_write(pb)

    call stop_check(pb,it)


    !----------------Output onestep to screen and ox file(snap_shot)--------------------
    if(mod(it-1,pb%output%ntout) == 0 .or. it == pb%itstop) then
      call screen_write(pb,it,dt_did)
      call ox_write(pb)
    endif


  enddo

end subroutine solve



end module solve_master

