! From http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/rkf45_f90.txt
! http://jean-pierre.moreau.pagesperso-orange.fr/f_eqdiff.html
!**************************************************************************
!* Collection of Fortran 90 subroutines to Integrate a System of Ordinary *
!* Differential Equations By the Runge-Kutta-Fehlberg method (simple or   *
!* double precision).                                                     *
!* ---------------------------------------------------------------------- *
!* REFERENCE:     H A Watts and L F Shampine,                             *
!*                Sandia Laboratories,                                    *
!*                Albuquerque, New Mexico.                                *
!*                                                                        *
!**************************************************************************

module ode_rk45

  use derivs_all, only: derivs

  implicit none
  private

  public :: rkf45_d

contains

  subroutine fehl_d ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )
  !
  !*******************************************************************************
  !
  !! FEHL_D takes one Fehlberg fourth-fifth order step (double precision).
  !
  !
  !  Discussion:
  !
  !    FEHL integrates a system of NEQN first order ordinary differential
  !    equations of the form
  !      dY(i)/dT = F(T,Y(1),---,Y(NEQN))
  !    where the initial values Y and the initial derivatives
  !    YP are specified at the starting point T.
  !
  !    FEHL advances the solution over the fixed step H and returns
  !    the fifth order (sixth order accurate locally) solution
  !    approximation at T+H in array S.
  !
  !    The formulas have been grouped to control loss of significance.
  !    FEHL should be called with an H not smaller than 13 units of
  !    roundoff in T so that the various independent arguments can be
  !    distinguished.
  !
  !  Author:
  !
  !    H A Watts and L F Shampine,
  !    Sandia Laboratories,
  !    Albuquerque, New Mexico.
  !
  !    RKF45 is primarily designed to solve non-stiff and mildly stiff
  !    differential equations when derivative evaluations are inexpensive.
  !    RKF45 should generally not be used when the user is demanding
  !    high accuracy.
  !
  !  Parameters:
  !
  !    Input, external F, a subroutine of the form
  !      subroutine f(t,y,yp)
  !    to evaluate the derivatives.
  !      YP(I) = dY(I) / dT
  !
  !    Input, integer NEQN, the number of equations to be integrated.
  !
  !    Input, double precision Y(NEQN), the current value of the dependent
  !    variable.
  !
  !    Input, double precision T, the current value of the independent variable.
  !
  !    Input, double precision H, the step size to take.
  !
  !    Input, double precision YP(NEQN), the current value of the derivative
  !    of the dependent variable.
  !
  !    Output, double precision F1(NEQN), F2(NEQN), F3(NEQN), F4(NEQN),
  !    F5(NEQN) are arrays of dimension NEQN which are needed
  !    for internal storage.
  !
  !    Output, double precision S(NEQN), the computed estimate of the solution
  !    at T+H.
  !
    implicit none
  !
    integer neqn
  !
    double precision ch
    external f
    double precision f1(neqn)
    double precision f2(neqn)
    double precision f3(neqn)
    double precision f4(neqn)
    double precision f5(neqn)
    double precision h
    double precision s(neqn)
    double precision t
    double precision y(neqn)
    double precision yp(neqn)
  !
    ch = h / 4.0D+00

    f5(1:neqn) = y(1:neqn) + ch * yp(1:neqn)

    call f ( t + ch, f5, f1 )

    ch = 3.0D+00 * h / 32.0D+00

    f5(1:neqn) = y(1:neqn) + ch * ( yp(1:neqn) + 3.0D+00 * f1(1:neqn) )

    call f ( t + 3.0D+00 * h / 8.0D+00, f5, f2 )

    ch = h / 2197.0D+00

    f5(1:neqn) = y(1:neqn) + ch * ( 1932.0D+00 * yp(1:neqn) &
      + ( 7296.0D+00 * f2(1:neqn) - 7200.0D+00 * f1(1:neqn) ) )

    call f ( t + 12.0D+00 * h / 13.0D+00, f5, f3 )

    ch = h / 4104.0D+00

    f5(1:neqn) = y(1:neqn) + ch * ( ( 8341.0D+00 * yp(1:neqn) &
      - 845.0D+00 * f3(1:neqn) ) + ( 29440.0D+00 * f2(1:neqn) &
      - 32832.0D+00 * f1(1:neqn) ) )

    call f ( t + h, f5, f4 )

    ch = h / 20520.0D+00

    f1(1:neqn) = y(1:neqn) + ch * ( ( -6080.0D+00 * yp(1:neqn) &
      + ( 9295.0D+00 * f3(1:neqn) - 5643.0D+00 * f4(1:neqn) ) ) &
      + ( 41040.0D+00 * f1(1:neqn) - 28352.0D+00 * f2(1:neqn) ) )

    call f ( t + h / 2.0D+00, f1, f5 )
  !
  !  Ready to compute the approximate solution at T+H.
  !
    ch = h / 7618050.0D+00

    s(1:neqn) = y(1:neqn) + ch * ( ( 902880.0D+00 * yp(1:neqn) &
      + ( 3855735.0D+00 * f3(1:neqn) - 1371249.0D+00 * f4(1:neqn) ) ) &
      + ( 3953664.0D+00 * f2(1:neqn) + 277020.0D+00 * f5(1:neqn) ) )

    return
  end subroutine fehl_d

  subroutine rkf45_d ( f, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )
  !
  !*******************************************************************************
  !
  !! RKF45_D carries out the Runge-Kutta-Fehlberg method (double precision).
  !
  !
  !  Author:
  !
  !    H A Watts and L F Shampine,
  !    Sandia Laboratories,
  !    Albuquerque, New Mexico.
  !
  !    RKF45 is primarily designed to solve non-stiff and mildly stiff
  !    differential equations when derivative evaluations are inexpensive.
  !    RKF45 should generally not be used when the user is demanding
  !    high accuracy.
  !
  !  Abstract:
  !
  !    RKF45 integrates a system of NEQN first order ordinary differential
  !    equations of the form:
  !
  !      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
  !
  !    where the Y(1:NEQN) are given at T.
  !
  !    Typically the subroutine is used to integrate from T to TOUT but it
  !    can be used as a one-step integrator to advance the solution a
  !    single step in the direction of TOUT.  On return, the parameters in
  !    the call list are set for continuing the integration.  The user has
  !    only to call RKF45 again (and perhaps define a new value for TOUT).
  !    Actually, RKF45 is an interfacing routine which calls subroutine
  !    RKFS for the solution.  RKFS in turn calls subroutine FEHL which
  !    computes an approximate solution over one step.
  !
  !  Reference:
  !
  !    E. Fehlberg,
  !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
  !    NASA Technical Report R-315.
  !
  !    L F Shampine, H A Watts, S Davenport,
  !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
  !    Sandia Laboratories Report SAND75-0182,
  !    To appear in SIAM Review.
  !
  !  Parameters:
  !
  !    Input, external F, a subroutine of the form
  !      subroutine f ( t, y, yp )
  !      double precision t
  !      double precision y(*)
  !      double precision yp(*)
  !    to evaluate the derivatives YP(I) = dY(I) / dT
  !
  !    Input, integer NEQN, the number of equations to be integrated.
  !
  !    Input/output, double precision Y(NEQN), the solution vector at T.
  !
  !    Input/output, double precision T, the independent variable.
  !
  !    Input, double precision TOUT, the output point at which solution is
  !    desired.
  !
  !    Input, double precision RELERR, ABSERR, the relative and absolute error
  !    tolerances for the local error test.  At each step the code requires:
  !      abs(local error) <= relerr*abs(y) + abserr
  !    for each component of the local error and solution vectors
  !
  !    Output, integer IFLAG, indicator for status of integration.
  !
  !    Workspace, double precision WORK(3+6*NEQN), an array to hold information
  !    internal to RKF45 which is necessary for subsequent calls.
  !
  !    Workspace, integer IWORK(5), an array used to hold information internal
  !    to RKF45 which is necessary for subsequent calls.
  !
  !
  !  first call
  !
  !    The user must provide storage in his calling program for the arrays
  !    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
  !    declare f in an external statement, supply subroutine f(t,y,yp) and
  !    initialize the following parameters-
  !
  !      neqn -- number of equations to be integrated.  (neqn >= 1)
  !
  !      y(*) -- vector of initial conditions
  !
  !      t -- starting point of integration , must be a variable
  !
  !      tout -- output point at which solution is desired.
  !            t=tout is allowed on the first call only, in which case
  !            rkf45 returns with iflag=2 if continuation is possible.
  !
  !      relerr,abserr -- relative and absolute local error tolerances
  !            which must be non-negative. relerr must be a variable while
  !            abserr may be a constant. the code should normally not be
  !            used with relative error control smaller than about 1.e-8 .
  !            to avoid limiting precision difficulties the code requires
  !            relerr to be larger than an internally computed relative
  !            error parameter which is machine dependent. in particular,
  !            pure absolute error is not permitted. if a smaller than
  !            allowable value of relerr is attempted, rkf45 increases
  !            relerr appropriately and returns control to the user before
  !            continuing the integration.
  !
  !      iflag -- +1,-1  indicator to initialize the code for each new
  !            problem. normal input is +1. the user should set iflag=-1
  !            only when one-step integrator control is essential. in this
  !            case, rkf45 attempts to advance the solution a single step
  !            in the direction of tout each time it is called. since this
  !            mode of operation results in extra computing overhead, it
  !            should be avoided unless needed.
  !
  !
  !  output from rkf45
  !
  !      y(*) -- solution at t
  !      t -- last point reached in integration.
  !      iflag = 2 -- integration reached tout. indicates successful retur
  !                   and is the normal mode for continuing integration.
  !            =-2 -- a single successful step in the direction of tout
  !                   has been taken. normal mode for continuing
  !                   integration one step at a time.
  !            = 3 -- integration was not completed because relative error
  !                   tolerance was too small. relerr has been increased
  !                   appropriately for continuing.
  !            = 4 -- integration was not completed because more than
  !                   3000 derivative evaluations were needed. this
  !                   is approximately 500 steps.
  !            = 5 -- integration was not completed because solution
  !                   vanished making a pure relative error test
  !                   impossible. must use non-zero abserr to continue.
  !                   using the one-step integration mode for one step
  !                   is a good way to proceed.
  !            = 6 -- integration was not completed because requested
  !                   accuracy could not be achieved using smallest
  !                   allowable stepsize. user must increase the error
  !                   tolerance before continued integration can be
  !                   attempted.
  !            = 7 -- it is likely that rkf45 is inefficient for solving
  !                   this problem. too much output is restricting the
  !                   natural stepsize choice. use the one-step integrator
  !                   mode.
  !            = 8 -- invalid input parameters
  !                   this indicator occurs if any of the following is
  !                   satisfied -   neqn <= 0
  !                                 t=tout  and  iflag /= +1 or -1
  !                                 relerr or abserr < 0.
  !                                 iflag == 0  or  < -2  or  > 8
  !      work(*),iwork(*) -- information which is usually of no interest
  !                   to the user but necessary for subsequent calls.
  !                   work(1),...,work(neqn) contain the first derivatives
  !                   of the solution vector y at t. work(neqn+1) contains
  !                   the stepsize h to be attempted on the next step.
  !                   iwork(1) contains the derivative evaluation counter.
  !
  !
  !  subsequent calls
  !
  !    RKF45 returns with all information needed to continue
  !    the integration. if the integration reached tout, the user need onl
  !    define a new tout and call RKF45 again.  In the one-step integrator
  !    mode (iflag=-2) the user must keep in mind that each step taken is
  !    in the direction of the current tout.  Upon reaching tout (indicated
  !    by changing iflag to 2),the user must then define a new tout and
  !    reset iflag to -2 to continue in the one-step integrator mode.
  !
  !    If the integration was not completed but the user still wants to
  !    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
  !    the relerr parameter has been adjusted appropriately for continuing
  !    the integration. in the case of iflag=4 the function counter will
  !    be reset to 0 and another 3000 function evaluations are allowed.
  !
  !    However,in the case iflag=5, the user must first alter the error
  !    criterion to use a positive value of abserr before integration can
  !    proceed. if he does not,execution is terminated.
  !
  !    Also,in the case iflag=6, it is necessary for the user to reset
  !    iflag to 2 (or -2 when the one-step integration mode is being used)
  !    as well as increasing either abserr,relerr or both before the
  !    integration can be continued. if this is not done, execution will
  !    be terminated. the occurrence of iflag=6 indicates a trouble spot
  !    (solution is changing rapidly,singularity may be present) and it
  !    often is inadvisable to continue.
  !
  !    If iflag=7 is encountered, the user should use the one-step
  !    integration mode with the stepsize determined by the code.
  !    If the user insists upon continuing the integration with RKF45,
  !    he must reset iflag to 2 before calling RKF45 again. otherwise,
  !    execution will be terminated.
  !
  !    If iflag=8 is obtained, integration can not be continued unless
  !    the invalid input parameters are corrected.
  !
  !    The arrays work and iwork contain information
  !    required for subsequent integration, and should not be altered.
  !
    implicit none
  !
    integer neqn
  !
    double precision abserr
    external f
    integer iflag
    integer iwork(5)
    integer k1
    integer k1m
    integer k2
    integer k3
    integer k4
    integer k5
    integer k6
    double precision relerr
    double precision t
    double precision tout
    double precision work(6*neqn+3)
    double precision y(neqn)
  !
  !  Compute indices for the splitting of the work array
  !
    k1m = neqn + 1
    k1 = k1m + 1
    k2 = k1 + neqn
    k3 = k2 + neqn
    k4 = k3 + neqn
    k5 = k4 + neqn
    k6 = k5 + neqn
  !
  !  This interfacing routine merely relieves the user of a long
  !  calling list via the splitting apart of two working storage arrays.
  !
    call rkfs_d ( f, neqn, y, t, tout, relerr, abserr, iflag, work(1), &
      work(k1m), work(k1), work(k2), work(k3), work(k4), work(k5), work(k6), &
      work(k6+1), iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )

    return
  end subroutine rkf45_d

  subroutine rkfs_d ( f, neqn, y, t, tout, relerr, abserr, iflag, yp, h, &
    f1, f2, f3, f4, f5, savre, savae, nfe, kop, init, jflag, kflag )
  !
  !*******************************************************************************
  !
  !! RKFS_D implements the Runge-Kutta-Fehlberg method (double precision).
  !
  !
  !  Discussion:
  !
  !    RKFS integrates a system of first order ordinary differential
  !    equations as described in the comments for RKF45.
  !
  !    The arrays yp, f1, f2, f3, f4, and f5 (of dimension at least neqn) and
  !    the variables h, savre, savae, nfe, kop, init, jflag and kflag are used
  !    internally by the code and appear in the call list to eliminate
  !    local retention of variables between calls.  Accordingly, they
  !    should not be altered.  Items of possible interest are
  !
  !      YP  - the derivative of the solution vector at T;
  !      H   - an appropriate stepsize to be used for the next step;
  !      NFE - the number of derivative function evaluations.
  !
  !    The expense is controlled by restricting the number
  !    of function evaluations to be approximately MAXNFE.
  !    As set, this corresponds to about 500 steps.
  !
  !    REMIN is the minimum acceptable value of RELERR.  Attempts
  !    to obtain higher accuracy with this subroutine are usually
  !    very expensive and often unsuccessful.
  !
    implicit none
  !
    integer neqn
  !
    double precision a
    double precision abserr
    double precision ae
    double precision dt
    double precision ee
    double precision eeoet
    double precision eps
    double precision esttol
    double precision et
    external f
    double precision f1(neqn)
    double precision f2(neqn)
    double precision f3(neqn)
    double precision f4(neqn)
    double precision f5(neqn)
    double precision h
    logical hfaild
    double precision hmin
    integer iflag
    integer init
    integer jflag
    integer k
    integer kflag
    integer kop
    integer, parameter :: maxnfe = 3000
    integer mflag
    integer nfe
    logical output
    double precision relerr
    double precision, parameter :: remin = 1.0D-12
    double precision rer
    double precision s
    double precision savae
    double precision savre
    double precision scale
    double precision t
    double precision tol
    double precision toln
    double precision tout
    double precision y(neqn)
    double precision yp(neqn)
    double precision ypk
  !
  !  Check the input parameters.
  !
    eps = epsilon ( eps )

    if ( neqn < 1 ) then
      iflag = 8
      return
    end if

    if ( relerr < 0.0D+00 ) then
      iflag = 8
      return
    end if

    if ( abserr < 0.0D+00 ) then
      iflag = 8
      return
    end if

    mflag = abs ( iflag )

    if ( abs ( iflag ) < 1 .or. abs ( iflag ) > 8 ) then
      iflag = 8
      return
    end if
  !
  !  Is this the first call?
  !
    if ( mflag == 1 ) then
      go to 50
    end if
  !
  !  Check continuation possibilities
  !
    if ( t == tout .and. kflag /= 3 ) then
      iflag = 8
      return
    end if

    if ( mflag /= 2 ) then
      go to 25
    end if
  !
  !  iflag = +2 or -2
  !
    if ( kflag == 3 ) go to 45
    if ( init == 0 ) go to 45
    if ( kflag == 4 ) go to 40

    if ( kflag == 5 .and. abserr == 0.0D+00 ) then
      stop
    end if

    if ( kflag == 6 .and. relerr <= savre .and. abserr <= savae ) then
      stop
    end if

    go to 50
  !
  !  iflag = 3,4,5,6,7 or 8
  !
     25 continue

    if ( iflag == 3 ) go to 45
    if ( iflag == 4 ) go to 40
    if ( iflag == 5 .and. abserr > 0.0D+00 ) go to 45
  !
  !  Integration cannot be continued since user did not respond to
  !  the instructions pertaining to iflag=5,6,7 or 8
  !
    stop
  !
  !  Reset function evaluation counter
  !
  40 continue

    nfe = 0
    if ( mflag == 2 ) then
      go to 50
    end if
  !
  !  Reset flag value from previous call
  !
  45 continue

    iflag = jflag

    if ( kflag == 3 ) then
      mflag = abs ( iflag )
    end if
  !
  !  Save input iflag and set continuation flag for subsequent input checking.
  !
  50 continue

    jflag = iflag
    kflag = 0
  !
  !  Save relerr and abserr for checking input on subsequent calls
  !
    savre = relerr
    savae = abserr
  !
  !  Restrict relative error tolerance to be at least as large as
  !  2*eps+remin to avoid limiting precision difficulties arising
  !  from impossible accuracy requests
  !
    rer = 2.0D+00 * epsilon ( rer ) + remin
  !
  !  The relative error tolerance is too small.
  !
    if ( relerr < rer ) then
      relerr = rer
      iflag = 3
      kflag = 3
      return
    end if

    dt = tout - t

    if ( mflag == 1 ) go to 60
    if ( init == 0 ) go to 65
    go to 80
  !
  !  Initialization:
  !    set initialization completion indicator, init
  !    set indicator for too many output points, kop
  !    evaluate initial derivatives
  !    set counter for function evaluations, nfe
  !    evaluate initial derivatives
  !    set counter for function evaluations, nfe
  !    estimate starting stepsize
  !
     60 continue

    init = 0
    kop = 0
    a = t
    call f ( a, y, yp )
    nfe = 1

    if ( t == tout ) then
      iflag = 2
      return
    end if

     65 continue

    init = 1
    h = abs ( dt )
    toln = 0.0D+00
    do k = 1, neqn
      tol = relerr * abs ( y(k) ) + abserr
      if ( tol > 0.0D+00 ) then
        toln = tol
        ypk = abs ( yp(k) )
        if ( ypk * h**5 > tol) then
          h = ( tol / ypk )**0.2D+00
        end if
      end if
    end do

    if ( toln <= 0.0D+00 ) then
      h = 0.0D+00
    end if

    h = max ( h, 26.0D+00 * eps * max ( abs ( t ), abs ( dt ) ) )
    jflag =  sign ( 2, iflag )
  !
  !  Set stepsize for integration in the direction from T to TOUT.
  !
     80 continue

    h = sign ( h, dt )
  !
  !  Test to see if RKF45 is being severely impacted by too many output points.
  !
    if ( abs ( h ) >= 2.0D+00 * abs ( dt ) ) then
      kop = kop + 1
    end if
  !
  !  Unnecessary frequency of output.
  !
    if ( kop == 100 ) then
      kop = 0
      iflag = 7
      return
    end if
  !
  !  If too close to output point, extrapolate and return.
  !
    if ( abs ( dt ) <= 26.0D+00 * eps * abs ( t ) ) then
      y(1:neqn) = y(1:neqn) + dt * yp(1:neqn)
      a = tout
      call f ( a, y, yp )
      nfe = nfe + 1
      t = tout
      iflag = 2
      return
    end if
  !
  !  Initialize output point indicator.
  !
    output = .false.
  !
  !  To avoid premature underflow in the error tolerance function,
  !  scale the error tolerances
  !
    scale = 2.0D+00 / relerr
    ae = scale * abserr
  !
  !  Step by step integration.
  !
    100 continue

    hfaild = .false.
  !
  !  Set smallest allowable stepsize.
  !
   ! hmin = 26.0D+00 * eps * abs ( t )
  ! SEISMIC: since t gets pretty large, hmin tends to become unacceptably
  ! large in magnitude. hmin is now set to the numerical accuracy at the risk
  ! of having significant underflow (extremely small time steps)...
    hmin = 26.0D+00 * eps
  !
  !  Adjust stepsize if necessary to hit the output point.
  !  Look ahead two steps to avoid drastic changes in the stepsize and
  !  thus lessen the impact of output points on the code.
  !
    dt = tout - t
    if ( abs ( dt ) >= 2.0D+00 * abs ( h ) ) go to 200
  !
  !  The next successful step will complete the integration to the output point.
  !
    if ( abs ( dt ) <= abs ( h ) ) then
      output = .true.
      h = dt
      go to 200
    end if

    h = 0.5D+00 * dt
  !
  !  Core integrator for taking a single step
  !
  !     The tolerances have been scaled to avoid premature underflow in
  !     computing the error tolerance function ET.
  !     To avoid problems with zero crossings, relative error is measured
  !     using the average of the magnitudes of the solution at the
  !     beginning and end of a step.
  !     The error estimate formula has been grouped to control loss of
  !     significance.
  !
  !     To distinguish the various arguments, H is not permitted
  !     to become smaller than 26 units of roundoff in T.
  !     Practical limits on the change in the stepsize are enforced to
  !     smooth the stepsize selection process and to avoid excessive
  !     chattering on problems having discontinuities.
  !     To prevent unnecessary failures, the code uses 9/10 the stepsize
  !     it estimates will succeed.
  !
  !     After a step failure, the stepsize is not allowed to increase for
  !     the next attempted step.  This makes the code more efficient on
  !     problems having discontinuities and more effective in general
  !     since local extrapolation is being used and extra caution seems
  !     warranted.
  !
  !     Test number of derivative function evaluations.
  !     If okay, try to advance the integration from T to T+H.
  !
    200 continue
  !
  !  Too much work.
  !
    if ( nfe > maxnfe ) then
      iflag = 4
      kflag = 4
      return
    end if
  !
  !  Advance an approximate solution over one step of length H.
  !
    call fehl_d ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1 )
    nfe = nfe + 5
  !
  !  Compute and test allowable tolerances versus local error estimates
  !  and remove scaling of tolerances.  Note that relative error is
  !  measured with respect to the average of the magnitudes of the
  !  solution at the beginning and end of the step.
  !
    eeoet = 0.0D+00

    do k = 1, neqn

      et = abs ( y(k) ) + abs ( f1(k) ) + ae

      if ( et <= 0.0D+00 ) then
        iflag = 5
        return
      end if

      ee = abs ( ( -2090.0D+00 * yp(k) + ( 21970.0D+00 * f3(k) &
        - 15048.0D+00 * f4(k) ) ) &
        + ( 22528.0D+00 * f2(k) - 27360.0D+00 * f5(k) ) )

      eeoet = max ( eeoet, ee / et )

    end do

    esttol = abs ( h ) * eeoet * scale / 752400.0D+00

    ! SEISMIC: if the integration result contains NaNs,
    ! reduce step size and try again
    if (any(isnan(f1))) then
      s = 0.1D+00
      go to 259
    endif

    if ( esttol <= 1.0D+00 ) then
      go to 260
    end if
  !
  !  Unsuccessful step.  Reduce the stepsize, try again.
  !  The decrease is limited to a factor of 1/10.
  !

    if ( esttol < 59049.0D+00 ) then
      s = 0.9D+00 / esttol**0.2D+00
    else
      s = 0.1D+00
    end if

    259 continue

    h = s * h

    hfaild = .true.
    output = .false.

    if ( abs ( h ) < hmin ) then
      iflag = 6
      kflag = 6
      return
    else
      go to 200
    end if
  !
  !  Successful step.  Store solution at T+H and evaluate derivatives there.
  !
    260 continue

    t = t + h
    y(1:neqn) = f1(1:neqn)
    a = t
    call f ( a, y, yp )
    nfe = nfe + 1
  !
  !  Choose next stepsize.  The increase is limited to a factor of 5.
  !  If step failure has just occurred, next stepsize is not allowed to increase
  !
    if ( esttol > 0.0001889568D+00 ) then
      s = 0.9D+00 / esttol**0.2D+00
    else
      s = 5.0D+00
    end if

    if ( hfaild ) then
      s = min ( s, 1.0D+00 )
    end if

    h = sign ( max ( s * abs ( h ), hmin ), h )
  !
  !  End of core integrator
  !
  !  Should we take another step?
  !
    if ( output ) then
      t = tout
      iflag = 2
    end if

    if ( iflag > 0 ) go to 100
  !
  !  Integration successfully completed
  !
  !  one-step mode
  !
    iflag = - 2

    return
  end subroutine rkfs_d

  ! subroutine timestamp ( )
  ! !
  ! !*******************************************************************************
  ! !
  ! !! TIMESTAMP prints the current YMDHMS date as a time stamp.
  ! !
  ! !
  ! !  Example:
  ! !
  ! !    May 31 2001   9:45:54.872 AM
  ! !
  ! !  Modified:
  ! !
  ! !    31 May 2001
  ! !
  ! !  Author:
  ! !
  ! !    John Burkardt
  ! !
  ! !  Parameters:
  ! !
  ! !    None
  ! !
  !   implicit none
  ! !
  !   character ( len = 8 ) ampm
  !   integer d
  !   character ( len = 8 ) date
  !   integer h
  !   integer m
  !   integer mm
  !   character ( len = 9 ), parameter, dimension(12) :: month = (/ &
  !     'January  ', 'February ', 'March    ', 'April    ', &
  !     'May      ', 'June     ', 'July     ', 'August   ', &
  !     'September', 'October  ', 'November ', 'December ' /)
  !   integer n
  !   integer s
  !   character ( len = 10 )  time
  !   integer values(8)
  !   integer y
  !   character ( len = 5 ) zone
  ! !
  !   call date_and_time ( date, time, zone, values )
  !
  !   y = values(1)
  !   m = values(2)
  !   d = values(3)
  !   h = values(5)
  !   n = values(6)
  !   s = values(7)
  !   mm = values(8)
  !
  !   if ( h < 12 ) then
  !     ampm = 'AM'
  !   else if ( h == 12 ) then
  !     if ( n == 0 .and. s == 0 ) then
  !       ampm = 'Noon'
  !     else
  !       ampm = 'PM'
  !     end if
  !   else
  !     h = h - 12
  !     if ( h < 12 ) then
  !       ampm = 'PM'
  !     else if ( h == 12 ) then
  !       if ( n == 0 .and. s == 0 ) then
  !         ampm = 'Midnight'
  !       else
  !         ampm = 'AM'
  !       end if
  !     end if
  !   end if
  !
  !   write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
  !     trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
  !
  !   return
  ! end subroutine timestamp

end module ode_rk45
