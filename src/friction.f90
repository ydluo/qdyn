module friction

! This is the only module that needs modifications to implement a new friction law
! Follow the instructions in the comment blocks below that start with "! new friction law:"
!
! Assumptions:
!   The friction coefficient (mu) and the state variable rate (dtheta/dt)
!   depend on slip velocity (v) and state variable (theta)
!     mu = f(v,theta)
!     dtheta/dt = g(v,theta)
!   All friction properties can be spatially non-uniform

! <SEISMIC>
! Rate-and-state friction law i_rns_law = 3 refers to Chen's microphysical model,
! with itheta_law = 3 to diffusion controlled pressure solution, and itheta_law = 4
! to dissolution controlled pressure solution creep.
! In this interpretation, rate-and-state friction results from the interplay
! between compaction by pressure solution creep and dilation by granular flow
! At low sliding velocities, pressure solution creep is relatively fast compared
! to granular dilatation, so that the gouge densifies (porosity reduces). By
! contrast, at high sliding velocities, the gouge dilates and 'softens', i.e.
! supports lower shear stress. At a given velocity, a steady-state porosity can
! be attained that corresponds to the porosity at which compaction balances
! dilatation, and hence results in a steady-state shear strength. Velocity-weakening
! behaviour emerges naturally from this model. To make the transition to Velocity-
! strengthening, a dislocation-glide type of model is adopted.
!
! For details and derivation of the model, see:
!   - Niemeijer & Spiers (2007), J. Geophys. Res., doi:10.1029/2007JB00500
!   - Chen & Spiers (2016), JGR: Solid Earth, doi:10.1002/2016JB013470
!
! In the view of rate-and-state friction, we take the porosity as the
! (microstructural) state, so that the constitutive relations can be implemented
! in the current QDYN framework by defining the partial- and time-derivatives of
! the friction coefficient and the porosity (state). The (spatially non-uniform)
! values of the material/gouge properties and initial porosity are read from the
! input script, and stored in pb%chen_params. Parameters defined specifically
! in this file are:
!   - tan_psi: dilatation angle
!   - mu_tilde: friction coefficient of the grain-boundary (rate-strengthening)
!   - e_ps_dot: (compaction) strain rate due to pressure solution creep
!
! WARNING: with the current implementation, Chen's model fails to simulate
! slide-hold-slide tests (i.e. the load-point velocity is set to zero). This
! should be worked out to ensure a proper implementation of the kinematic
! and dynamic equations
!
! </SEISMIC>

  use problem_class, only : problem_type

  implicit none
  private

  public  :: set_theta_star, friction_mu, dmu_dv_dtheta, dtheta_dt
  private  :: calc_tan_psi, calc_mu_tilde, calc_e_ps

contains

!--------------------------------------------------------------------------------------
subroutine set_theta_star(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%i_rns_law)

  case (0)
    pb%theta_star = pb%dc/pb%v_star

  case (1)
    pb%theta_star = pb%dc/pb%v2

  case (3) ! SEISMIC: Chen's friction law does not use theta_star
    pb%theta_star = 1

! new friction law:
!  case(xxx)
!    implement here your definition of theta_star (could be none)
!    pb%theta_star = ...

  case default
    stop 'set_theta_star: unknown friction law type'
  end select

end subroutine set_theta_star

!--------------------------------------------------------------------------------------
function friction_mu(v,theta,pb) result(mu)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: mu
  ! SEISMIC: additional parameters for Chen model
  double precision, dimension(pb%mesh%nn) :: mu_tilde, tan_psi, mu_gr, mu_ps


  select case (pb%i_rns_law)

  case (0)
    mu = pb%mu_star - pb%a*log(pb%v_star/v) + pb%b*log(theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/v+1d0) + pb%b*log(theta/pb%theta_star+1d0)

  case (3) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)   ! Dilatation angle
    mu_tilde = calc_mu_tilde(v, pb)     ! Grain-boundary friction
    ! SEISMIC: Calculate the frictional strength controlled by granular flow
    mu_gr = (mu_tilde + tan_psi)/(1 - mu_tilde*tan_psi)

    ! SEISMIC: Next, calculate the shear strength (divided by normal stress)
    ! when the strength is controlled by viscous flow (no granular flow)
    select case (pb%itheta_law)

    case (3) ! SEISMIC: Chen's model, diffusion controlled
      mu_ps = v*((2*pb%chen_params%phi0 - 2*theta)**2)/(pb%chen_params%w*pb%chen_params%IPS_const_diff*pb%sigma)
    case (4) ! SEISMIC: Chen's model, dissolution controlled
      mu_ps = (pb%chen_params%phi0 - theta)/(pb%sigma*pb%chen_params%phi0*pb%chen_params%IPS_const_diss2)* &
              log(abs(v)/(pb%chen_params%w*pb%chen_params%IPS_const_diss1) + 1)
    case default
      write(6,*) "friction_mu: Chen's friction model is selected (i_rns_law == 3),"
      write(6,*) "but itheta_law is unsupported (must be either 3 or 4)"
      write(6,*) "3 = diffusion-, 4 = dissolution controlled pressure solution creep"
      stop

    end select

    ! SEISMIC: this weird looking function approximates the min(mu_gr, mu_ps) function,
    ! but is continuously differentiable, so that the solver doesn't complain.
    ! Higher order powers yield a closer approximation
    ! By approximation, the mechanism that offers the lowest shear strength controls
    ! the overall shear strength of the gouge
    mu = 1/(1/mu_gr**5 + 1/mu_ps**5)**(1.0/5.0)

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu

!--------------------------------------------------------------------------------------
function dtheta_dt(v,theta,pb) result(dth_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: dth_dt, omega
  double precision, dimension(pb%mesh%nn) :: tan_psi, e_ps_dot

  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

  case (3, 4) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)   ! Dilatation angle
    e_ps_dot = calc_e_ps(pb%sigma, theta, .true., pb)   ! Compaction rate
    ! SEISMIC: the nett rate of change of porosity is given by the relation
    ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
    ! the compaction rate by pressure solution (e_ps_dot) minus the dilatation
    ! rate by granular flow (tan_psi*v/w). Note that compaction is measured
    ! positive, dilatation negative
    dth_dt = (1 - theta)*(tan_psi*v/pb%chen_params%w - e_ps_dot)


! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end function dtheta_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,theta,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv, dmu_dtheta

  ! SEISMIC: define some extra parameters for Chen's friction law
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, dth_dt, denom, y_ps, delta
  double precision, dimension(pb%mesh%nn) :: mu_star, dummy_var, f_phi_sqrt, V_gr, y_gr, tau
  double precision, dimension(pb%mesh%nn) :: dv_dtau_ps, dv_dtheta_ps, dv_dtau, dv_dtheta

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

  case(3) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)         ! Dilatation angle
    mu_tilde = calc_mu_tilde(v, pb)           ! Grain-boundary friction
    mu_star = pb%chen_params%mu_tilde_star

    !tau = friction_mu(v, theta, pb)*pb%sigma
    tau = pb%kernel%k1 * (pb%v_star * pb%time - (pb%slip + v*pb%dt_did))
    !tau = pb%tau

    ! Pre-compute the denominator for efficiency
    denom = 1.0/(pb%chen_params%a*(pb%sigma+tau*tan_psi))
    ! Pre-compute a dummy variable that we will re-use a few times
    dummy_var = (tau*(1 - mu_star*tan_psi) - pb%sigma*(mu_star + tan_psi))*denom

    ! Pre-compute the square root of the porosity function. This will serve as
    ! a basis for calculating the full porosity functions for either itheta_law
    f_phi_sqrt = 1.0/(2*(pb%chen_params%phi0 - theta))

    ! Granular flow _velocity_ [unit m/s]
    ! NOTE: it is physically more correct to use a reference strain rate
    ! and use velocity = strain_rate * thickness instead
    V_gr = (1.0/pb%chen_params%slowness_star)*exp(dummy_var)
    ! Shear strain rate [1/s] due to viscous flow (pressure solution)
    y_ps = calc_e_ps(tau, theta, .false., pb)

    ! Next, calculate the partial derivatives dv_dtau and dv_dtheta for the
    ! pressure solution creep. The functional form depends on the rate-limiting
    ! mechanism (diffusion, dissolution, or precipitation)
    ! Note that the pressure solution strain rate is multiplied by the
    ! the gouge layer thickness (pb%chen_params%w) to obtain a velocity
    select case (pb%itheta_law)

    case (3) ! Diffusion controlled pressure solution
      dv_dtau_ps = pb%chen_params%w*pb%chen_params%IPS_const_diff*f_phi_sqrt**2
      dv_dtheta_ps = 4*pb%chen_params%w*pb%chen_params%IPS_const_diff*tau*f_phi_sqrt**3
    case (4) ! Dissolution controlled pressure solution
      dv_dtau_ps =  pb%chen_params%w*(y_ps + pb%chen_params%IPS_const_diss1)* &
                    pb%chen_params%IPS_const_diss2*2*pb%chen_params%phi0*f_phi_sqrt
      dv_dtheta_ps =  4*pb%chen_params%w*(y_ps + pb%chen_params%IPS_const_diss1)* &
                      pb%chen_params%IPS_const_diss2*tau*pb%chen_params%phi0*f_phi_sqrt**2
    case default
      write(6,*) "dmu_dv_dtheta: Chen's friction model is selected (i_rns_law == 3),"
      write(6,*) "but itheta_law is unsupported (must be either 3 or 4)"
      write(6,*) "3 = diffusion-, 4 = dissolution controlled pressure solution creep"
      stop

    end select

    ! The partial derivatives dv_dtau and dv_dtheta of the overall slip velocity
    ! are the sum of the partial derivatives of the pressure solution and
    ! granular flow velocities.
    dv_dtau = dv_dtau_ps + V_gr*((1-mu_star*tan_psi)*denom - tan_psi*dummy_var*pb%chen_params%a*denom)
    dv_dtheta = dv_dtheta_ps - V_gr*(2*pb%chen_params%H*(pb%sigma + mu_star*tau)*denom + &
                2*pb%chen_params%H*tau*dummy_var*pb%chen_params%a*denom)

    ! The partial derivatives dmu_dv and dmu_dtheta can then be obtained from dv_dtau and dv_dtheta
    dmu_dv = 1/(pb%sigma*dv_dtau)
    dmu_dtheta = dmu_dv*dv_dtheta

! new friction law:
!  case(xxx)
!    implement here the partial derivatives of the friction coefficient
!    dmu_dtheta = ...
!    dmu_dv = ...

  case default
    stop 'dmu_dv_dtheta: unknown friction law type'
  end select

end subroutine dmu_dv_dtheta

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
subroutine compute_velocity(v,tau,theta,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, theta
  double precision, dimension(pb%mesh%nn), intent(out) :: v

  ! SEISMIC: define some extra parameters for Chen's friction law
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, dth_dt, denom, y_ps, delta
  double precision, dimension(pb%mesh%nn) :: mu_star, dummy_var, f_phi_sqrt, V_gr, y_gr, tau
  double precision, dimension(pb%mesh%nn) :: dv_dtau_ps, dv_dtheta_ps, dv_dtau, dv_dtheta

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

  case(3) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)         ! Dilatation angle
    mu_tilde = calc_mu_tilde(v, pb)           ! Grain-boundary friction
    mu_star = pb%chen_params%mu_tilde_star

    !tau = friction_mu(v, theta, pb)*pb%sigma
    tau = pb%kernel%k1 * (pb%v_star * pb%time - (pb%slip + v*pb%dt_did))
    !tau = pb%tau

    ! Pre-compute the denominator for efficiency
    denom = 1.0/(pb%chen_params%a*(pb%sigma+tau*tan_psi))
    ! Pre-compute a dummy variable that we will re-use a few times
    dummy_var = (tau*(1 - mu_star*tan_psi) - pb%sigma*(mu_star + tan_psi))*denom

    ! Pre-compute the square root of the porosity function. This will serve as
    ! a basis for calculating the full porosity functions for either itheta_law
    f_phi_sqrt = 1.0/(2*(pb%chen_params%phi0 - theta))

    ! Granular flow _velocity_ [unit m/s]
    ! NOTE: it is physically more correct to use a reference strain rate
    ! and use velocity = strain_rate * thickness instead
    V_gr = (1.0/pb%chen_params%slowness_star)*exp(dummy_var)
    ! Shear strain rate [1/s] due to viscous flow (pressure solution)
    y_ps = calc_e_ps(tau, theta, .false., pb)

    ! Next, calculate the partial derivatives dv_dtau and dv_dtheta for the
    ! pressure solution creep. The functional form depends on the rate-limiting
    ! mechanism (diffusion, dissolution, or precipitation)
    ! Note that the pressure solution strain rate is multiplied by the
    ! the gouge layer thickness (pb%chen_params%w) to obtain a velocity
    select case (pb%itheta_law)

    case (3) ! Diffusion controlled pressure solution
      dv_dtau_ps = pb%chen_params%w*pb%chen_params%IPS_const_diff*f_phi_sqrt**2
      dv_dtheta_ps = 4*pb%chen_params%w*pb%chen_params%IPS_const_diff*tau*f_phi_sqrt**3
    case (4) ! Dissolution controlled pressure solution
      dv_dtau_ps =  pb%chen_params%w*(y_ps + pb%chen_params%IPS_const_diss1)* &
                    pb%chen_params%IPS_const_diss2*2*pb%chen_params%phi0*f_phi_sqrt
      dv_dtheta_ps =  4*pb%chen_params%w*(y_ps + pb%chen_params%IPS_const_diss1)* &
                      pb%chen_params%IPS_const_diss2*tau*pb%chen_params%phi0*f_phi_sqrt**2
    case default
      write(6,*) "dmu_dv_dtheta: Chen's friction model is selected (i_rns_law == 3),"
      write(6,*) "but itheta_law is unsupported (must be either 3 or 4)"
      write(6,*) "3 = diffusion-, 4 = dissolution controlled pressure solution creep"
      stop

    end select

    ! The partial derivatives dv_dtau and dv_dtheta of the overall slip velocity
    ! are the sum of the partial derivatives of the pressure solution and
    ! granular flow velocities.
    dv_dtau = dv_dtau_ps + V_gr*((1-mu_star*tan_psi)*denom - tan_psi*dummy_var*pb%chen_params%a*denom)
    dv_dtheta = dv_dtheta_ps - V_gr*(2*pb%chen_params%H*(pb%sigma + mu_star*tau)*denom + &
                2*pb%chen_params%H*tau*dummy_var*pb%chen_params%a*denom)

    ! The partial derivatives dmu_dv and dmu_dtheta can then be obtained from dv_dtau and dv_dtheta
    dmu_dv = 1/(pb%sigma*dv_dtau)
    dmu_dtheta = dmu_dv*dv_dtheta

! new friction law:
!  case(xxx)
!    implement here the partial derivatives of the friction coefficient
!    dmu_dtheta = ...
!    dmu_dv = ...

  case default
    stop 'dmu_dv_dtheta: unknown friction law type'
  end select

end subroutine compute_velocity

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
function calc_tan_psi(theta,pb) result(tan_psi)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: theta
  double precision, dimension(pb%mesh%nn) :: tan_psi

  tan_psi = 2*pb%chen_params%H*(pb%chen_params%phi0 - theta)

end function calc_tan_psi

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the grain-boundary coefficient of friction, which is
! logarithmically rate-strengthening, assuming some atomic scale jump process
!--------------------------------------------------------------------------------------
function calc_mu_tilde(v,pb) result(mu_tilde)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v
  double precision, dimension(pb%mesh%nn) :: mu_tilde

  ! SEISMIC: slowness_star is defined as 1/V_star in the input script
  ! NOTE: it is physically more correct to use a reference strain rate
  ! and use velocity = strain_rate * thickness instead
  mu_tilde = pb%chen_params%mu_tilde_star + pb%chen_params%a * log(abs(v)*pb%chen_params%slowness_star)

end function calc_mu_tilde

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the rate of pressure solution as a function of effective
! normal or shear stress, and porosity. To prevent porosities from going negative
! compaction rates can be truncated by an error function. The logical 'truncate'
! indicates if truncation is needed. This only applies to compaction in the
! normal direction, and not to shear deformation (viscous flow), hence truncate
! should be false when calculating shear strain rate
!--------------------------------------------------------------------------------------
function calc_e_ps(sigma,theta,truncate,pb) result(e_ps_dot)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, theta
  double precision, dimension(pb%mesh%nn) :: e_ps_dot
  logical, intent(in) :: truncate

  select case (pb%itheta_law)

  case (3) ! SEISMIC: Chen's model, diffusion controlled
    e_ps_dot = pb%chen_params%IPS_const_diff*pb%sigma/(2*pb%chen_params%phi0 - 2*theta)**2
  case (4) ! SEISMIC: Chen's model, dissolution controlled
    e_ps_dot =  pb%chen_params%IPS_const_diss1*(exp(pb%chen_params%IPS_const_diss2*sigma* &
                pb%chen_params%phi0/(pb%chen_params%phi0 - theta)) - 1)
  case default
    stop 'calc_e_ps: unsupported pressure solution law (3 = diffusion, 4 = dissolution controlled)'
  end select

  ! If truncation is required, use an error function to set the strain rate
  ! to zero. A cut-off porosity of about 3% is chosen here, which corresponds
  ! to the percolation limit of aggregates
  if (truncate .eqv. .true.) e_ps_dot = 0.5*(1 + erf(100*(theta-0.03)))*e_ps_dot

end function calc_e_ps

end module friction
