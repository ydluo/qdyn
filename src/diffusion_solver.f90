!===============================================================================
!
! SEISMIC: spectral solver module for diffusion problem
!
! This module handles integration of the diffusion equations for heat and
! fluid pressure, following the method proposed by Noda & Lapusta (2010) [N&L].
! The diffusion normal to the fault plane is solved in a discretised spectral
! domain, which is numerically stable for large time steps (as opposed to
! discretisation methods in the spatial domain, such as FEM, FD, etc.). This
! procedure allows for efficient modelling of frictional heating and associated
! pressurisation of the pore fluid (thermal pressurisation: TP)
! This module is compatible with both the rate-and-state friction and the
! Chen-Niemeijer-Spiers frameworks
!
! TODO: add a coefficient that regulates the dilatancy hardening component
!
! References:
!
! [N&L]
! Noda, H., and N. Lapusta (2010): Three‐dimensional earthquake sequence
! simulations with evolving temperature and pore pressure due to shear
! heating: Effect of heterogeneous hydraulic diffusivity,
! J. Geophys. Res., 115, B12314, doi:10.1029/2010JB007780
!
!===============================================================================

module diffusion_solver

  use problem_class, only : problem_type
  use mesh, only : spectral_mesh_type

  implicit none
  private

  public :: init_tp, update_PT, update_PT_final

contains

!===============================================================================
! SEISMIC: initialisation routine for the thermal pressurisation model
!===============================================================================
subroutine init_tp(pb)

  type(problem_type), intent(inout) :: pb

  call init_spectral_mesh(pb%tp%mesh)
  call init_variables(pb)
  call init_source(pb)

  ! Calculate remaining parameters (not necessarily time-constant)
  if (pb%i_rns_law == 3) then
    ! If using CNS, pass the initial porosity
    call calc_params(pb%theta, pb)
  else
    ! if using RSF, pass one
    call calc_params(pb%theta-pb%theta+1.0, pb)
  endif

end subroutine init_tp

!===============================================================================
! SEISMIC: initiate discretisation of the spectral domain for the spectral
! diffusion solver of [N&L]. Mesh parameters are hard-coded at the top of
! mesh.f90
!===============================================================================
subroutine init_spectral_mesh(sm)

  use constants, only : PI

  type(spectral_mesh_type), intent(inout) :: sm

  integer :: i
  double precision :: spi

  spi = sqrt(2.0/PI)

  allocate(sm%lw(sm%Nl), sm%F_inv(sm%Nl))

  do i=1,sm%Nl
    ! Construct logarithmic grid of dimensionless wavenumbers [N&L, Eqn. 14]
    sm%lw(i) = sm%lw_max * exp(-sm%Dlogl * (sm%Nl - i) )

    ! Construct kernel for inverse Fourier transform [N&L, Eqn. 17]
    if (i == 1) then
      sm%F_inv(i) = spi*sm%lw(i)*(1 + 0.5*sm%Dlogl)
    else if (i == sm%Nl) then
      sm%F_inv(i) = spi*sm%lw(i)*0.5*sm%Dlogl
    else
      sm%F_inv(i) = spi*sm%lw(i)*sm%Dlogl
    endif
  enddo

end subroutine init_spectral_mesh

!===============================================================================
! SEISMIC: initiate Gaussian distribution of heat source and dilatation in the
! discretised spectral domain
!===============================================================================
subroutine init_source(pb)

  use constants, only : PI

  type(problem_type), intent(inout) :: pb

  double precision :: spi = 1.0/sqrt(2*PI)
  integer :: i, j, n

  allocate (pb%tp%Omega(pb%mesh%nn*pb%tp%mesh%Nl))

  ! Loop over all fault segments
  do i=1,pb%mesh%nn
    ! Loop over all spectral elements
    do j=1,pb%tp%mesh%Nl
      n = (i-1)*pb%tp%mesh%Nl+j
      ! Fourier transform of Gaussian heat source [N&L, Eqn. 13]
      ! Note, Omega is pre-multiplied with w
      pb%tp%Omega(n) = pb%tp%w(i) * exp(-0.5*pb%tp%mesh%lw(j)**2)*spi
    enddo
  enddo

end subroutine init_source

!===============================================================================
! SEISMIC: initiate vectors of pressure and temperature (and spectral
! equivalents Pi and Theta), and vectors of previous values
!===============================================================================
subroutine init_variables(pb)

  type(problem_type), intent(inout) :: pb

  integer :: i

  ! Allocate variables
  allocate (  pb%tp%P(pb%mesh%nn), &
              pb%tp%T(pb%mesh%nn), &
              pb%tp%inv_w(pb%mesh%nn), &
              pb%tp%Pi(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%Theta(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%PiTheta(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%dP_dt(pb%mesh%nn), &
              pb%dtheta_dt(pb%mesh%nn), &
              pb%dtheta2_dt(pb%mesh%nn) )

  ! Allocate parameters
  allocate (  pb%tp%inv_rhoc(pb%mesh%nn), &
              pb%tp%alpha_th(pb%mesh%nn), &
              pb%tp%alpha_hy(pb%mesh%nn), &
              pb%tp%Lam(pb%mesh%nn), &
              pb%tp%Lam_prime(pb%mesh%nn), &
              pb%tp%Lam_T(pb%mesh%nn), &
              pb%tp%phi_b(pb%mesh%nn) )

  ! Allocate previous values
  allocate (  pb%tp%Theta_prev(pb%mesh%nn*pb%tp%mesh%Nl), &
              pb%tp%PiTheta_prev(pb%mesh%nn*pb%tp%mesh%Nl) )

  do i=1,pb%mesh%nn
    pb%tp%inv_rhoc = 1.0/pb%tp%rhoc(i)
    pb%tp%inv_w(i) = 1.0/pb%tp%w(i)
    pb%tp%P(i) = pb%tp%P_a(i)
    pb%tp%T(i) = pb%tp%T_a(i)
    pb%tp%dP_dt(i) = 0d0
    pb%dtheta_dt(i) = 0d0
    pb%dtheta2_dt(i) = 0d0
  enddo

  do i=1,pb%mesh%nn*pb%tp%mesh%Nl
    pb%tp%Pi(i) = 0d0
    pb%tp%Theta(i) = 0d0
  enddo

end subroutine init_variables

!===============================================================================
! SEISMIC: calculate (time-variable) parameters. For the rate-and-state
! framework, phi = 1. For the CNS framework, phi is variable
!===============================================================================
subroutine calc_params(phi, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: phi

  ! Thermal diffusivity
  pb%tp%alpha_th = pb%tp%k_t*pb%tp%inv_rhoc
  ! Hydraulic diffusivity
  pb%tp%alpha_hy = pb%tp%k_p/(pb%tp%eta*phi*pb%tp%beta)
  ! Thermal expansion moduli (by lack of a better term)
  pb%tp%Lam = pb%tp%l/pb%tp%beta
  pb%tp%Lam_prime = pb%tp%Lam * pb%tp%alpha_th/(pb%tp%alpha_hy - pb%tp%alpha_th)
  pb%tp%Lam_T = (pb%tp%Lam + pb%tp%Lam_prime)*pb%tp%inv_rhoc
  ! Specific storativity
  pb%tp%phi_b = pb%tp%dilat_factor/(phi*pb%tp%beta)

end subroutine calc_params

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spatial domain
!===============================================================================
subroutine update_PT(tau_y,phi_dot,phi,v,dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi, v
  double precision :: dt

  ! Compute PiTheta and Theta for step t+dt
  call solve_spectral(tau_y, phi_dot, phi, v, dt, pb)

end subroutine update_PT

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spatial domain for a full
! time step (i.e. not the intermediate RK solver steps) by a mid-point
! integration scheme
!===============================================================================
subroutine update_PT_final(dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi, dP_dt
  double precision :: dt

  ! Calculate mid-point values of tau_y, phi_dot, phi between t and t+dt
  tau_y = pb%tau * 0.5 * pb%v * pb%tp%inv_w
  if (pb%i_rns_law == 3) then
    ! CNS model: include porosity
    phi = pb%theta
    ! Multiply dilatancy rate with an arbitrary constant f >= 0 to control
    ! the amount of dilatancy hardening
    phi_dot = pb%dtheta_dt
  else
    ! RSF: ignore state
    phi = 1d0
    phi_dot = 0d0
  endif

  ! Compute PiTheta and Theta for step t+dt
  call solve_spectral(tau_y, phi_dot, phi, pb%v, dt, pb)

  ! Update initial values of Theta and PiTheta for next integration step
  pb%tp%Theta_prev = pb%tp%Theta
  pb%tp%PiTheta_prev = pb%tp%PiTheta

end subroutine update_PT_final

!===============================================================================
! SEISMIC: solve for P(t+dt) and T(t+dt) in the spectral domain
!===============================================================================
subroutine solve_spectral(tau_y,phi_dot,phi,v,dt,pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%mesh%nn) :: tau_y, phi_dot, phi, v
  double precision :: dt, A_T, B_T, exp_T, A_P, B_P, exp_P
  double precision :: dPi, dTheta, dPiTheta
  integer :: i, j, n

  ! Update parameters in spatial domain
  call calc_params(phi, pb)

  ! Reset values of Theta and PiTheta. For each solver step, integration
  ! is done from t to t+dt, with variable dt. Hence, the initial values of
  ! Theta and PiTheta should always correspond to those at time t
  pb%tp%Theta = pb%tp%Theta_prev
  pb%tp%PiTheta = pb%tp%PiTheta_prev

  ! Loop over all fault segments
  do i=1,pb%mesh%nn
    pb%tp%P(i) = pb%tp%P_a(i)
    pb%tp%T(i) = pb%tp%T_a(i)
    pb%tp%dP_dt(i) = 0d0

    if (v(i) > 1e-5) then
      ! Loop over all spectral elements
      do j=1,pb%tp%mesh%Nl
        n = (i-1)*pb%tp%mesh%Nl+j

        ! Calculate F(t+dt) = B*(1-exp(-Adt))/A + F(t)*exp(-Adt)
        ! assuming constant A, B over the duration of dt [N&L, Eqn. 10]

        ! Temperature-related parameters in spectral domain
        A_T = pb%tp%alpha_th(i)*(pb%tp%mesh%lw(j)*pb%tp%inv_w(i))**2
        B_T = tau_y(i)*pb%tp%Omega(n)*pb%tp%inv_rhoc(i)
        exp_T = exp(-A_t*dt)
        ! Update Ttheta(t+dt)
        pb%tp%Theta(n) = B_T*(1.0 - exp_T)/A_T + pb%tp%Theta(n)*exp_T
        pb%tp%T(i) = pb%tp%T(i) + pb%tp%mesh%F_inv(j)*pb%tp%inv_w(i)*pb%tp%Theta(n)

        ! Pressure-related parameters in spectral domain
        A_P = pb%tp%alpha_hy(i)*(pb%tp%mesh%lw(j)*pb%tp%inv_w(i))**2
        B_P = ( pb%tp%Lam_T(i)*tau_y(i) - pb%tp%phi_b(i)*phi_dot(i) )*pb%tp%Omega(n)
        exp_P = exp(-A_P*dt)
        ! Update PiTtheta(t+dt)
        ! Note that PiTheta contains the spectral representation of
        ! Pi + Lambda_prime*Theta, where Pi is the Fourier transform of P
        ! and Theta is the Fourier transform of T [see N&L, Eqn. 5 and 7]
        pb%tp%PiTheta(n) = B_P*(1.0 - exp_P)/A_P + pb%tp%PiTheta(n)*exp_P
        pb%tp%P(i) =  pb%tp%P(i) + pb%tp%mesh%F_inv(j)*pb%tp%inv_w(i)* &
                      (pb%tp%PiTheta(n) - pb%tp%Lam_prime(i)*pb%tp%Theta(n))

        ! Update dP/dt
        dTheta = -A_T*pb%tp%Theta(n) + B_T
        dPiTheta = -A_P*pb%tp%PiTheta(n) + B_P
        dPi = dPiTheta - pb%tp%Lam_prime(i)*dTheta
        pb%tp%dP_dt(i) = pb%tp%dP_dt(i) + pb%tp%mesh%F_inv(j)*dPi*pb%tp%inv_w(i)
      enddo
    else
      do j=1,pb%tp%mesh%Nl
        n = (i-1)*pb%tp%mesh%Nl+j
        pb%tp%Theta(n) = 0d0
        pb%tp%PiTheta(n) = 0d0
      enddo
    endif
  enddo

end subroutine solve_spectral

end module diffusion_solver
