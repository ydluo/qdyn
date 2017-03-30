module ode_lsoda_main

  use derivs_all, only: derivs
  use ode_lsoda_aux1
  use ode_lsoda_aux2

  implicit none

contains
  !*==DLSODE.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DLSODE
        SUBROUTINE DLSODE(F,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt,&
       &                  Rwork,Lrw,Iwork,Liw,JAC,Mf)
        IMPLICIT NONE
  !*--DLSODE6
  !*** Start of declarations inserted by SPAG
        !REAL DPREPJ , DSOLSY
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC
        INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
       &        Mf
        DOUBLE PRECISION Y , T , Tout , Rtol , Atol , Rwork
        DIMENSION Neq(*) , Y(*) , Rtol(*) , Atol(*) , Rwork(Lrw) ,        &
       &          Iwork(Liw)
  !***BEGIN PROLOGUE  DLSODE
  !***PURPOSE  Livermore Solver for Ordinary Differential Equations.
  !            DLSODE solves the initial-value problem for stiff or
  !            nonstiff systems of first-order ODE's,
  !               dy/dt = f(t,y),   or, in component form,
  !               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
  !***CATEGORY  I1A
  !***TYPE      DOUBLE PRECISION (SLSODE-S, DLSODE-D)
  !***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
  !             STIFF, NONSTIFF
  !***AUTHOR  Hindmarsh, Alan C., (LLNL)
  !             Center for Applied Scientific Computing, L-561
  !             Lawrence Livermore National Laboratory
  !             Livermore, CA 94551.
  !***DESCRIPTION
  !
  !     NOTE: The "Usage" and "Arguments" sections treat only a subset of
  !           available options, in condensed fashion.  The options
  !           covered and the information supplied will support most
  !           standard uses of DLSODE.
  !
  !           For more sophisticated uses, full details on all options are
  !           given in the concluding section, headed "Long Description."
  !           A synopsis of the DLSODE Long Description is provided at the
  !           beginning of that section; general topics covered are:
  !           - Elements of the call sequence; optional input and output
  !           - Optional supplemental routines in the DLSODE package
  !           - internal COMMON block
  !
  ! *Usage:
  !     Communication between the user and the DLSODE package, for normal
  !     situations, is summarized here.  This summary describes a subset
  !     of the available options.  See "Long Description" for complete
  !     details, including optional communication, nonstandard options,
  !     and instructions for special situations.
  !
  !     A sample program is given in the "Examples" section.
  !
  !     Refer to the argument descriptions for the definitions of the
  !     quantities that appear in the following sample declarations.
  !
  !     For MF = 10,
  !        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20)
  !     For MF = 21 or 22,
  !        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ)
  !     For MF = 24 or 25,
  !        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ,
  !       *                                         LIW = 20 + NEQ)
  !
  !        EXTERNAL F, JAC
  !        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW),
  !       *         LIW, MF
  !        DOUBLE PRECISION Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW)
  !
  !        CALL DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
  !       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
  !
  ! *Arguments:
  !     F     :EXT    Name of subroutine for right-hand-side vector f.
  !                   This name must be declared EXTERNAL in calling
  !                   program.  The form of F must be:
  !
  !                   SUBROUTINE  F (NEQ, T, Y, YDOT)
  !                   INTEGER  NEQ
  !                   DOUBLE PRECISION  T, Y(*), YDOT(*)
  !
  !                   The inputs are NEQ, T, Y.  F is to set
  !
  !                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)),
  !                                                     i = 1, ..., NEQ .
  !
  !     NEQ   :IN     Number of first-order ODE's.
  !
  !     Y     :INOUT  Array of values of the y(t) vector, of length NEQ.
  !                   Input:  For the first call, Y should contain the
  !                           values of y(t) at t = T. (Y is an input
  !                           variable only if ISTATE = 1.)
  !                   Output: On return, Y will contain the values at the
  !                           new t-value.
  !
  !     T     :INOUT  Value of the independent variable.  On return it
  !                   will be the current value of t (normally TOUT).
  !
  !     TOUT  :IN     Next point where output is desired (.NE. T).
  !
  !     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or
  !                   an array.
  !
  !     RTOL  :IN     Relative tolerance parameter (scalar).
  !
  !     ATOL  :IN     Absolute tolerance parameter (scalar or array).
  !                   If ITOL = 1, ATOL need not be dimensioned.
  !                   If ITOL = 2, ATOL must be dimensioned at least NEQ.
  !
  !                   The estimated local error in Y(i) will be controlled
  !                   so as to be roughly less (in magnitude) than
  !
  !                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
  !                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
  !
  !                   Thus the local error test passes if, in each
  !                   component, either the absolute error is less than
  !                   ATOL (or ATOL(i)), or the relative error is less
  !                   than RTOL.
  !
  !                   Use RTOL = 0.0 for pure absolute error control, and
  !                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative
  !                   error control.  Caution:  Actual (global) errors may
  !                   exceed these local tolerances, so choose them
  !                   conservatively.
  !
  !     ITASK :IN     Flag indicating the task DLSODE is to perform.
  !                   Use ITASK = 1 for normal computation of output
  !                   values of y at t = TOUT.
  !
  !     ISTATE:INOUT  Index used for input and output to specify the state
  !                   of the calculation.
  !                   Input:
  !                    1   This is the first call for a problem.
  !                    2   This is a subsequent call.
  !                   Output:
  !                    1   Nothing was done, because TOUT was equal to T.
  !                    2   DLSODE was successful (otherwise, negative).
  !                        Note that ISTATE need not be modified after a
  !                        successful return.
  !                   -1   Excess work done on this call (perhaps wrong
  !                        MF).
  !                   -2   Excess accuracy requested (tolerances too
  !                        small).
  !                   -3   Illegal input detected (see printed message).
  !                   -4   Repeated error test failures (check all
  !                        inputs).
  !                   -5   Repeated convergence failures (perhaps bad
  !                        Jacobian supplied or wrong choice of MF or
  !                        tolerances).
  !                   -6   Error weight became zero during problem
  !                        (solution component i vanished, and ATOL or
  !                        ATOL(i) = 0.).
  !
  !     IOPT  :IN     Flag indicating whether optional inputs are used:
  !                   0   No.
  !                   1   Yes.  (See "Optional inputs" under "Long
  !                       Description," Part 1.)
  !
  !     RWORK :WORK   Real work array of length at least:
  !                   20 + 16*NEQ                    for MF = 10,
  !                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
  !                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
  !
  !     LRW   :IN     Declared length of RWORK (in user's DIMENSION
  !                   statement).
  !
  !     IWORK :WORK   Integer work array of length at least:
  !                   20        for MF = 10,
  !                   20 + NEQ  for MF = 21, 22, 24, or 25.
  !
  !                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the
  !                   lower and upper Jacobian half-bandwidths ML,MU.
  !
  !                   On return, IWORK contains information that may be
  !                   of interest to the user:
  !
  !            Name   Location   Meaning
  !            -----  ---------  -----------------------------------------
  !            NST    IWORK(11)  Number of steps taken for the problem so
  !                              far.
  !            NFE    IWORK(12)  Number of f evaluations for the problem
  !                              so far.
  !            NJE    IWORK(13)  Number of Jacobian evaluations (and of
  !                              matrix LU decompositions) for the problem
  !                              so far.
  !            NQU    IWORK(14)  Method order last used (successfully).
  !            LENRW  IWORK(17)  Length of RWORK actually required.  This
  !                              is defined on normal returns and on an
  !                              illegal input return for insufficient
  !                              storage.
  !            LENIW  IWORK(18)  Length of IWORK actually required.  This
  !                              is defined on normal returns and on an
  !                              illegal input return for insufficient
  !                              storage.
  !
  !     LIW   :IN     Declared length of IWORK (in user's DIMENSION
  !                   statement).
  !
  !     JAC   :EXT    Name of subroutine for Jacobian matrix (MF =
  !                   21 or 24).  If used, this name must be declared
  !                   EXTERNAL in calling program.  If not used, pass a
  !                   dummy name.  The form of JAC must be:
  !
  !                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
  !                   INTEGER  NEQ, ML, MU, NROWPD
  !                   DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
  !
  !                   See item c, under "Description" below for more
  !                   information about JAC.
  !
  !     MF    :IN     Method flag.  Standard values are:
  !                   10  Nonstiff (Adams) method, no Jacobian used.
  !                   21  Stiff (BDF) method, user-supplied full Jacobian.
  !                   22  Stiff method, internally generated full
  !                       Jacobian.
  !                   24  Stiff method, user-supplied banded Jacobian.
  !                   25  Stiff method, internally generated banded
  !                       Jacobian.
  !
  ! *Description:
  !     DLSODE solves the initial value problem for stiff or nonstiff
  !     systems of first-order ODE's,
  !
  !        dy/dt = f(t,y) ,
  !
  !     or, in component form,
  !
  !        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ))
  !                                                  (i = 1, ..., NEQ) .
  !
  !     DLSODE is a package based on the GEAR and GEARB packages, and on
  !     the October 23, 1978, version of the tentative ODEPACK user
  !     interface standard, with minor modifications.
  !
  !     The steps in solving such a problem are as follows.
  !
  !     a. First write a subroutine of the form
  !
  !           SUBROUTINE  F (NEQ, T, Y, YDOT)
  !           INTEGER  NEQ
  !           DOUBLE PRECISION  T, Y(*), YDOT(*)
  !
  !        which supplies the vector function f by loading YDOT(i) with
  !        f(i).
  !
  !     b. Next determine (or guess) whether or not the problem is stiff.
  !        Stiffness occurs when the Jacobian matrix df/dy has an
  !        eigenvalue whose real part is negative and large in magnitude
  !        compared to the reciprocal of the t span of interest.  If the
  !        problem is nonstiff, use method flag MF = 10.  If it is stiff,
  !        there are four standard choices for MF, and DLSODE requires the
  !        Jacobian matrix in some form.  This matrix is regarded either
  !        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the
  !        banded case, DLSODE requires two half-bandwidth parameters ML
  !        and MU. These are, respectively, the widths of the lower and
  !        upper parts of the band, excluding the main diagonal.  Thus the
  !        band consists of the locations (i,j) with
  !
  !           i - ML <= j <= i + MU ,
  !
  !        and the full bandwidth is ML + MU + 1 .
  !
  !     c. If the problem is stiff, you are encouraged to supply the
  !        Jacobian directly (MF = 21 or 24), but if this is not feasible,
  !        DLSODE will compute it internally by difference quotients (MF =
  !        22 or 25).  If you are supplying the Jacobian, write a
  !        subroutine of the form
  !
  !           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
  !           INTEGER  NEQ, ML, MU, NRWOPD
  !           DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
  !
  !        which provides df/dy by loading PD as follows:
  !        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
  !          the partial derivative of f(i) with respect to y(j).  (Ignore
  !          the ML and MU arguments in this case.)
  !        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
  !          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the
  !          rows of PD from the top down.
  !        - In either case, only nonzero elements need be loaded.
  !
  !     d. Write a main program that calls subroutine DLSODE once for each
  !        point at which answers are desired.  This should also provide
  !        for possible use of logical unit 6 for output of error messages
  !        by DLSODE.
  !
  !        Before the first call to DLSODE, set ISTATE = 1, set Y and T to
  !        the initial values, and set TOUT to the first output point.  To
  !        continue the integration after a successful return, simply
  !        reset TOUT and call DLSODE again.  No other parameters need be
  !        reset.
  !
  ! *Examples:
  !     The following is a simple example problem, with the coding needed
  !     for its solution by DLSODE. The problem is from chemical kinetics,
  !     and consists of the following three rate equations:
  !
  !        dy1/dt = -.04*y1 + 1.E4*y2*y3
  !        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
  !        dy3/dt = 3.E7*y2**2
  !
  !     on the interval from t = 0.0 to t = 4.E10, with initial conditions
  !     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
  !
  !     The following coding solves this problem with DLSODE, using
  !     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses
  !     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2
  !     has much smaller values.  At the end of the run, statistical
  !     quantities of interest are printed.
  !
  !        EXTERNAL  FEX, JEX
  !        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
  !       *         MF, NEQ
  !        DOUBLE PRECISION  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
  !        NEQ = 3
  !        Y(1) = 1.D0
  !        Y(2) = 0.D0
  !        Y(3) = 0.D0
  !        T = 0.D0
  !        TOUT = .4D0
  !        ITOL = 2
  !        RTOL = 1.D-4
  !        ATOL(1) = 1.D-6
  !        ATOL(2) = 1.D-10
  !        ATOL(3) = 1.D-6
  !        ITASK = 1
  !        ISTATE = 1
  !        IOPT = 0
  !        LRW = 58
  !        LIW = 23
  !        MF = 21
  !        DO 40 IOUT = 1,12
  !          CALL DLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
  !       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
  !          WRITE(6,20)  T, Y(1), Y(2), Y(3)
  !    20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
  !          IF (ISTATE .LT. 0)  GO TO 80
  !    40    TOUT = TOUT*10.D0
  !        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
  !    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
  !        STOP
  !    80  WRITE(6,90)  ISTATE
  !    90  FORMAT(///' Error halt.. ISTATE =',I3)
  !        STOP
  !        END
  !
  !        SUBROUTINE  FEX (NEQ, T, Y, YDOT)
  !        INTEGER  NEQ
  !        DOUBLE PRECISION  T, Y(3), YDOT(3)
  !        YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
  !        YDOT(3) = 3.D7*Y(2)*Y(2)
  !        YDOT(2) = -YDOT(1) - YDOT(3)
  !        RETURN
  !        END
  !
  !        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
  !        INTEGER  NEQ, ML, MU, NRPD
  !        DOUBLE PRECISION  T, Y(3), PD(NRPD,3)
  !        PD(1,1) = -.04D0
  !        PD(1,2) = 1.D4*Y(3)
  !        PD(1,3) = 1.D4*Y(2)
  !        PD(2,1) = .04D0
  !        PD(2,3) = -PD(1,3)
  !        PD(3,2) = 6.D7*Y(2)
  !        PD(2,2) = -PD(1,2) - PD(3,2)
  !        RETURN
  !        END
  !
  !     The output from this program (on a Cray-1 in single precision)
  !     is as follows.
  !
  !     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
  !     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
  !     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
  !     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
  !     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
  !     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
  !     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
  !     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
  !     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
  !     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01
  !     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
  !     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00
  !
  !     No. steps = 330,  No. f-s = 405,  No. J-s = 69
  !
  ! *Accuracy:
  !     The accuracy of the solution depends on the choice of tolerances
  !     RTOL and ATOL.  Actual (global) errors may exceed these local
  !     tolerances, so choose them conservatively.
  !
  ! *Cautions:
  !     The work arrays should not be altered between calls to DLSODE for
  !     the same problem, except possibly for the conditional and optional
  !     inputs.
  !
  ! *Portability:
  !     Since NEQ is dimensioned inside DLSODE, some compilers may object
  !     to a call to DLSODE with NEQ a scalar variable.  In this event,
  !     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL.
  !
  !     Note to Cray users:
  !     For maximum efficiency, use the CFT77 compiler.  Appropriate
  !     compiler optimization directives have been inserted for CFT77.
  !
  ! *Reference:
  !     Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
  !     Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
  !     (North-Holland, Amsterdam, 1983), pp. 55-64.
  !
  ! *Long Description:
  !     The following complete description of the user interface to
  !     DLSODE consists of four parts:
  !
  !     1.  The call sequence to subroutine DLSODE, which is a driver
  !         routine for the solver.  This includes descriptions of both
  !         the call sequence arguments and user-supplied routines.
  !         Following these descriptions is a description of optional
  !         inputs available through the call sequence, and then a
  !         description of optional outputs in the work arrays.
  !
  !     2.  Descriptions of other routines in the DLSODE package that may
  !         be (optionally) called by the user.  These provide the ability
  !         to alter error message handling, save and restore the internal
  !         COMMON, and obtain specified derivatives of the solution y(t).
  !
  !     3.  Descriptions of COMMON block to be declared in overlay or
  !         similar environments, or to be saved when doing an interrupt
  !         of the problem and continued solution later.
  !
  !     4.  Description of two routines in the DLSODE package, either of
  !         which the user may replace with his own version, if desired.
  !         These relate to the measurement of errors.
  !
  !
  !                         Part 1.  Call Sequence
  !                         ----------------------
  !
  !     Arguments
  !     ---------
  !     The call sequence parameters used for input only are
  !
  !        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
  !
  !     and those used for both input and output are
  !
  !        Y, T, ISTATE.
  !
  !     The work arrays RWORK and IWORK are also used for conditional and
  !     optional inputs and optional outputs.  (The term output here
  !     refers to the return from subroutine DLSODE to the user's calling
  !     program.)
  !
  !     The legality of input parameters will be thoroughly checked on the
  !     initial call for the problem, but not checked thereafter unless a
  !     change in input parameters is flagged by ISTATE = 3 on input.
  !
  !     The descriptions of the call arguments are as follows.
  !
  !     F        The name of the user-supplied subroutine defining the ODE
  !              system.  The system must be put in the first-order form
  !              dy/dt = f(t,y), where f is a vector-valued function of
  !              the scalar t and the vector y. Subroutine F is to compute
  !              the function f. It is to have the form
  !
  !                 SUBROUTINE F (NEQ, T, Y, YDOT)
  !                 DOUBLE PRECISION  T, Y(*), YDOT(*)
  !
  !              where NEQ, T, and Y are input, and the array YDOT =
  !              f(T,Y) is output.  Y and YDOT are arrays of length NEQ.
  !              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be
  !              declared EXTERNAL in the calling program.
  !
  !              Subroutine F may access user-defined quantities in
  !              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array
  !              (dimensioned in F) and/or Y has length exceeding NEQ(1).
  !              See the descriptions of NEQ and Y below.
  !
  !              If quantities computed in the F routine are needed
  !              externally to DLSODE, an extra call to F should be made
  !              for this purpose, for consistent and accurate results.
  !              If only the derivative dy/dt is needed, use DINTDY
  !              instead.
  !
  !     NEQ      The size of the ODE system (number of first-order
  !              ordinary differential equations).  Used only for input.
  !              NEQ may be decreased, but not increased, during the
  !              problem.  If NEQ is decreased (with ISTATE = 3 on input),
  !              the remaining components of Y should be left undisturbed,
  !              if these are to be accessed in F and/or JAC.
  !
  !              Normally, NEQ is a scalar, and it is generally referred
  !              to as a scalar in this user interface description.
  !              However, NEQ may be an array, with NEQ(1) set to the
  !              system size.  (The DLSODE package accesses only NEQ(1).)
  !              In either case, this parameter is passed as the NEQ
  !              argument in all calls to F and JAC.  Hence, if it is an
  !              array, locations NEQ(2),... may be used to store other
  !              integer data and pass it to F and/or JAC.  Subroutines
  !              F and/or JAC must include NEQ in a DIMENSION statement
  !              in that case.
  !
  !     Y        A real array for the vector of dependent variables, of
  !              length NEQ or more.  Used for both input and output on
  !              the first call (ISTATE = 1), and only for output on
  !              other calls.  On the first call, Y must contain the
  !              vector of initial values.  On output, Y contains the
  !              computed solution vector, evaluated at T. If desired,
  !              the Y array may be used for other purposes between
  !              calls to the solver.
  !
  !              This array is passed as the Y argument in all calls to F
  !              and JAC.  Hence its length may exceed NEQ, and locations
  !              Y(NEQ+1),... may be used to store other real data and
  !              pass it to F and/or JAC.  (The DLSODE package accesses
  !              only Y(1),...,Y(NEQ).)
  !
  !     T        The independent variable.  On input, T is used only on
  !              the first call, as the initial point of the integration.
  !              On output, after each call, T is the value at which a
  !              computed solution Y is evaluated (usually the same as
  !              TOUT).  On an error return, T is the farthest point
  !              reached.
  !
  !     TOUT     The next value of T at which a computed solution is
  !              desired.  Used only for input.
  !
  !              When starting the problem (ISTATE = 1), TOUT may be equal
  !              to T for one call, then should not equal T for the next
  !              call.  For the initial T, an input value of TOUT .NE. T
  !              is used in order to determine the direction of the
  !              integration (i.e., the algebraic sign of the step sizes)
  !              and the rough scale of the problem.  Integration in
  !              either direction (forward or backward in T) is permitted.
  !
  !              If ITASK = 2 or 5 (one-step modes), TOUT is ignored
  !              after the first call (i.e., the first call with
  !              TOUT .NE. T).  Otherwise, TOUT is required on every call.
  !
  !              If ITASK = 1, 3, or 4, the values of TOUT need not be
  !              monotone, but a value of TOUT which backs up is limited
  !              to the current internal T interval, whose endpoints are
  !              TCUR - HU and TCUR.  (See "Optional Outputs" below for
  !              TCUR and HU.)
  !
  !
  !     ITOL     An indicator for the type of error control.  See
  !              description below under ATOL.  Used only for input.
  !
  !     RTOL     A relative error tolerance parameter, either a scalar or
  !              an array of length NEQ.  See description below under
  !              ATOL.  Input only.
  !
  !     ATOL     An absolute error tolerance parameter, either a scalar or
  !              an array of length NEQ.  Input only.
  !
  !              The input parameters ITOL, RTOL, and ATOL determine the
  !              error control performed by the solver.  The solver will
  !              control the vector e = (e(i)) of estimated local errors
  !              in Y, according to an inequality of the form
  !
  !                 rms-norm of ( e(i)/EWT(i) ) <= 1,
  !
  !              where
  !
  !                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
  !
  !              and the rms-norm (root-mean-square norm) here is
  !
  !                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ).
  !
  !              Here EWT = (EWT(i)) is a vector of weights which must
  !              always be positive, and the values of RTOL and ATOL
  !              should all be nonnegative.  The following table gives the
  !              types (scalar/array) of RTOL and ATOL, and the
  !              corresponding form of EWT(i).
  !
  !              ITOL    RTOL      ATOL      EWT(i)
  !              ----    ------    ------    -----------------------------
  !              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
  !              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
  !              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
  !              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
  !
  !              When either of these parameters is a scalar, it need not
  !              be dimensioned in the user's calling program.
  !
  !              If none of the above choices (with ITOL, RTOL, and ATOL
  !              fixed throughout the problem) is suitable, more general
  !              error controls can be obtained by substituting
  !              user-supplied routines for the setting of EWT and/or for
  !              the norm calculation.  See Part 4 below.
  !
  !              If global errors are to be estimated by making a repeated
  !              run on the same problem with smaller tolerances, then all
  !              components of RTOL and ATOL (i.e., of EWT) should be
  !              scaled down uniformly.
  !
  !     ITASK    An index specifying the task to be performed.  Input
  !              only.  ITASK has the following values and meanings:
  !              1   Normal computation of output values of y(t) at
  !                  t = TOUT (by overshooting and interpolating).
  !              2   Take one step only and return.
  !              3   Stop at the first internal mesh point at or beyond
  !                  t = TOUT and return.
  !              4   Normal computation of output values of y(t) at
  !                  t = TOUT but without overshooting t = TCRIT.  TCRIT
  !                  must be input as RWORK(1).  TCRIT may be equal to or
  !                  beyond TOUT, but not behind it in the direction of
  !                  integration.  This option is useful if the problem
  !                  has a singularity at or beyond t = TCRIT.
  !              5   Take one step, without passing TCRIT, and return.
  !                  TCRIT must be input as RWORK(1).
  !
  !              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
  !              (within roundoff), it will return T = TCRIT (exactly) to
  !              indicate this (unless ITASK = 4 and TOUT comes before
  !              TCRIT, in which case answers at T = TOUT are returned
  !              first).
  !
  !     ISTATE   An index used for input and output to specify the state
  !              of the calculation.
  !
  !              On input, the values of ISTATE are as follows:
  !              1   This is the first call for the problem
  !                  (initializations will be done).  See "Note" below.
  !              2   This is not the first call, and the calculation is to
  !                  continue normally, with no change in any input
  !                  parameters except possibly TOUT and ITASK.  (If ITOL,
  !                  RTOL, and/or ATOL are changed between calls with
  !                  ISTATE = 2, the new values will be used but not
  !                  tested for legality.)
  !              3   This is not the first call, and the calculation is to
  !                  continue normally, but with a change in input
  !                  parameters other than TOUT and ITASK.  Changes are
  !                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
  !                  ML, MU, and any of the optional inputs except H0.
  !                  (See IWORK description for ML and MU.)
  !
  !              Note:  A preliminary call with TOUT = T is not counted as
  !              a first call here, as no initialization or checking of
  !              input is done.  (Such a call is sometimes useful for the
  !              purpose of outputting the initial conditions.)  Thus the
  !              first call for which TOUT .NE. T requires ISTATE = 1 on
  !              input.
  !
  !              On output, ISTATE has the following values and meanings:
  !               1  Nothing was done, as TOUT was equal to T with
  !                  ISTATE = 1 on input.
  !               2  The integration was performed successfully.
  !              -1  An excessive amount of work (more than MXSTEP steps)
  !                  was done on this call, before completing the
  !                  requested task, but the integration was otherwise
  !                  successful as far as T. (MXSTEP is an optional input
  !                  and is normally 500.)  To continue, the user may
  !                  simply reset ISTATE to a value >1 and call again (the
  !                  excess work step counter will be reset to 0).  In
  !                  addition, the user may increase MXSTEP to avoid this
  !                  error return; see "Optional Inputs" below.
  !              -2  Too much accuracy was requested for the precision of
  !                  the machine being used.  This was detected before
  !                  completing the requested task, but the integration
  !                  was successful as far as T. To continue, the
  !                  tolerance parameters must be reset, and ISTATE must
  !                  be set to 3. The optional output TOLSF may be used
  !                  for this purpose.  (Note:  If this condition is
  !                  detected before taking any steps, then an illegal
  !                  input return (ISTATE = -3) occurs instead.)
  !              -3  Illegal input was detected, before taking any
  !                  integration steps.  See written message for details.
  !                  (Note:  If the solver detects an infinite loop of
  !                  calls to the solver with illegal input, it will cause
  !                  the run to stop.)
  !              -4  There were repeated error-test failures on one
  !                  attempted step, before completing the requested task,
  !                  but the integration was successful as far as T.  The
  !                  problem may have a singularity, or the input may be
  !                  inappropriate.
  !              -5  There were repeated convergence-test failures on one
  !                  attempted step, before completing the requested task,
  !                  but the integration was successful as far as T. This
  !                  may be caused by an inaccurate Jacobian matrix, if
  !                  one is being used.
  !              -6  EWT(i) became zero for some i during the integration.
  !                  Pure relative error control (ATOL(i)=0.0) was
  !                  requested on a variable which has now vanished.  The
  !                  integration was successful as far as T.
  !
  !              Note:  Since the normal output value of ISTATE is 2, it
  !              does not need to be reset for normal continuation.  Also,
  !              since a negative input value of ISTATE will be regarded
  !              as illegal, a negative output value requires the user to
  !              change it, and possibly other inputs, before calling the
  !              solver again.
  !
  !     IOPT     An integer flag to specify whether any optional inputs
  !              are being used on this call.  Input only.  The optional
  !              inputs are listed under a separate heading below.
  !              0   No optional inputs are being used.  Default values
  !                  will be used in all cases.
  !              1   One or more optional inputs are being used.
  !
  !     RWORK    A real working array (double precision).  The length of
  !              RWORK must be at least
  !
  !                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
  !
  !              where
  !                 NYH = the initial value of NEQ,
  !              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
  !                       smaller value is given as an optional input),
  !                 LWM = 0           if MITER = 0,
  !                 LWM = NEQ**2 + 2  if MITER = 1 or 2,
  !                 LWM = NEQ + 2     if MITER = 3, and
  !                 LWM = (2*ML + MU + 1)*NEQ + 2
  !                                   if MITER = 4 or 5.
  !              (See the MF description below for METH and MITER.)
  !
  !              Thus if MAXORD has its default value and NEQ is constant,
  !              this length is:
  !              20 + 16*NEQ                    for MF = 10,
  !              22 + 16*NEQ + NEQ**2           for MF = 11 or 12,
  !              22 + 17*NEQ                    for MF = 13,
  !              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15,
  !              20 +  9*NEQ                    for MF = 20,
  !              22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
  !              22 + 10*NEQ                    for MF = 23,
  !              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
  !
  !              The first 20 words of RWORK are reserved for conditional
  !              and optional inputs and optional outputs.
  !
  !              The following word in RWORK is a conditional input:
  !              RWORK(1) = TCRIT, the critical value of t which the
  !                         solver is not to overshoot.  Required if ITASK
  !                         is 4 or 5, and ignored otherwise.  See ITASK.
  !
  !     LRW      The length of the array RWORK, as declared by the user.
  !              (This will be checked by the solver.)
  !
  !     IWORK    An integer work array.  Its length must be at least
  !              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
  !              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25).
  !              (See the MF description below for MITER.)  The first few
  !              words of IWORK are used for conditional and optional
  !              inputs and optional outputs.
  !
  !              The following two words in IWORK are conditional inputs:
  !              IWORK(1) = ML   These are the lower and upper half-
  !              IWORK(2) = MU   bandwidths, respectively, of the banded
  !                              Jacobian, excluding the main diagonal.
  !                         The band is defined by the matrix locations
  !                         (i,j) with i - ML <= j <= i + MU. ML and MU
  !                         must satisfy 0 <= ML,MU <= NEQ - 1. These are
  !                         required if MITER is 4 or 5, and ignored
  !                         otherwise.  ML and MU may in fact be the band
  !                         parameters for a matrix to which df/dy is only
  !                         approximately equal.
  !
  !     LIW      The length of the array IWORK, as declared by the user.
  !              (This will be checked by the solver.)
  !
  !     Note:  The work arrays must not be altered between calls to DLSODE
  !     for the same problem, except possibly for the conditional and
  !     optional inputs, and except for the last 3*NEQ words of RWORK.
  !     The latter space is used for internal scratch space, and so is
  !     available for use by the user outside DLSODE between calls, if
  !     desired (but not for use by F or JAC).
  !
  !     JAC      The name of the user-supplied routine (MITER = 1 or 4) to
  !              compute the Jacobian matrix, df/dy, as a function of the
  !              scalar t and the vector y.  (See the MF description below
  !              for MITER.)  It is to have the form
  !
  !                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
  !                 DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
  !
  !              where NEQ, T, Y, ML, MU, and NROWPD are input and the
  !              array PD is to be loaded with partial derivatives
  !              (elements of the Jacobian matrix) on output.  PD must be
  !              given a first dimension of NROWPD.  T and Y have the same
  !              meaning as in subroutine F.
  !
  !              In the full matrix case (MITER = 1), ML and MU are
  !              ignored, and the Jacobian is to be loaded into PD in
  !              columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
  !
  !              In the band matrix case (MITER = 4), the elements within
  !              the band are to be loaded into PD in columnwise manner,
  !              with diagonal lines of df/dy loaded into the rows of PD.
  !              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML
  !              and MU are the half-bandwidth parameters (see IWORK).
  !              The locations in PD in the two triangular areas which
  !              correspond to nonexistent matrix elements can be ignored
  !              or loaded arbitrarily, as they are overwritten by DLSODE.
  !
  !              JAC need not provide df/dy exactly. A crude approximation
  !              (possibly with a smaller bandwidth) will do.
  !
  !              In either case, PD is preset to zero by the solver, so
  !              that only the nonzero elements need be loaded by JAC.
  !              Each call to JAC is preceded by a call to F with the same
  !              arguments NEQ, T, and Y. Thus to gain some efficiency,
  !              intermediate quantities shared by both calculations may
  !              be saved in a user COMMON block by F and not recomputed
  !              by JAC, if desired.  Also, JAC may alter the Y array, if
  !              desired.  JAC must be declared EXTERNAL in the calling
  !              program.
  !
  !              Subroutine JAC may access user-defined quantities in
  !              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !              (dimensioned in JAC) and/or Y has length exceeding
  !              NEQ(1).  See the descriptions of NEQ and Y above.
  !
  !     MF       The method flag.  Used only for input.  The legal values
  !              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24,
  !              and 25.  MF has decimal digits METH and MITER:
  !                 MF = 10*METH + MITER .
  !
  !              METH indicates the basic linear multistep method:
  !              1   Implicit Adams method.
  !              2   Method based on backward differentiation formulas
  !                  (BDF's).
  !
  !              MITER indicates the corrector iteration method:
  !              0   Functional iteration (no Jacobian matrix is
  !                  involved).
  !              1   Chord iteration with a user-supplied full (NEQ by
  !                  NEQ) Jacobian.
  !              2   Chord iteration with an internally generated
  !                  (difference quotient) full Jacobian (using NEQ
  !                  extra calls to F per df/dy value).
  !              3   Chord iteration with an internally generated
  !                  diagonal Jacobian approximation (using one extra call
  !                  to F per df/dy evaluation).
  !              4   Chord iteration with a user-supplied banded Jacobian.
  !              5   Chord iteration with an internally generated banded
  !                  Jacobian (using ML + MU + 1 extra calls to F per
  !                  df/dy evaluation).
  !
  !              If MITER = 1 or 4, the user must supply a subroutine JAC
  !              (the name is arbitrary) as described above under JAC.
  !              For other values of MITER, a dummy argument can be used.
  !
  !     Optional Inputs
  !     ---------------
  !     The following is a list of the optional inputs provided for in the
  !     call sequence.  (See also Part 2.)  For each such input variable,
  !     this table lists its name as used in this documentation, its
  !     location in the call sequence, its meaning, and the default value.
  !     The use of any of these inputs requires IOPT = 1, and in that case
  !     all of these inputs are examined.  A value of zero for any of
  !     these optional inputs will cause the default value to be used.
  !     Thus to use a subset of the optional inputs, simply preload
  !     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively,
  !     and then set those of interest to nonzero values.
  !
  !     Name    Location   Meaning and default value
  !     ------  ---------  -----------------------------------------------
  !     H0      RWORK(5)   Step size to be attempted on the first step.
  !                        The default value is determined by the solver.
  !     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
  !                        default value is infinite.
  !     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
  !                        default value is 0.  (This lower bound is not
  !                        enforced on the final step before reaching
  !                        TCRIT when ITASK = 4 or 5.)
  !     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
  !                        is 12 if METH = 1, and 5 if METH = 2. (See the
  !                        MF description above for METH.)  If MAXORD
  !                        exceeds the default value, it will be reduced
  !                        to the default value.  If MAXORD is changed
  !                        during the problem, it may cause the current
  !                        order to be reduced.
  !     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
  !                        allowed during one call to the solver.  The
  !                        default value is 500.
  !     MXHNIL  IWORK(7)   Maximum number of messages printed (per
  !                        problem) warning that T + H = T on a step
  !                        (H = step size).  This must be positive to
  !                        result in a nondefault value.  The default
  !                        value is 10.
  !
  !     Optional Outputs
  !     ----------------
  !     As optional additional output from DLSODE, the variables listed
  !     below are quantities related to the performance of DLSODE which
  !     are available to the user.  These are communicated by way of the
  !     work arrays, but also have internal mnemonic names as shown.
  !     Except where stated otherwise, all of these outputs are defined on
  !     any successful return from DLSODE, and on any return with ISTATE =
  !     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3),
  !     they will be unchanged from their existing values (if any), except
  !     possibly for TOLSF, LENRW, and LENIW.  On any error return,
  !     outputs relevant to the error will be defined, as noted below.
  !
  !     Name   Location   Meaning
  !     -----  ---------  ------------------------------------------------
  !     HU     RWORK(11)  Step size in t last used (successfully).
  !     HCUR   RWORK(12)  Step size to be attempted on the next step.
  !     TCUR   RWORK(13)  Current value of the independent variable which
  !                       the solver has actually reached, i.e., the
  !                       current internal mesh point in t. On output,
  !                       TCUR will always be at least as far as the
  !                       argument T, but may be farther (if interpolation
  !                       was done).
  !     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0,
  !                       computed when a request for too much accuracy
  !                       was detected (ISTATE = -3 if detected at the
  !                       start of the problem, ISTATE = -2 otherwise).
  !                       If ITOL is left unaltered but RTOL and ATOL are
  !                       uniformly scaled up by a factor of TOLSF for the
  !                       next call, then the solver is deemed likely to
  !                       succeed.  (The user may also ignore TOLSF and
  !                       alter the tolerance parameters in any other way
  !                       appropriate.)
  !     NST    IWORK(11)  Number of steps taken for the problem so far.
  !     NFE    IWORK(12)  Number of F evaluations for the problem so far.
  !     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU
  !                       decompositions) for the problem so far.
  !     NQU    IWORK(14)  Method order last used (successfully).
  !     NQCUR  IWORK(15)  Order to be attempted on the next step.
  !     IMXER  IWORK(16)  Index of the component of largest magnitude in
  !                       the weighted local error vector ( e(i)/EWT(i) ),
  !                       on an error return with ISTATE = -4 or -5.
  !     LENRW  IWORK(17)  Length of RWORK actually required.  This is
  !                       defined on normal returns and on an illegal
  !                       input return for insufficient storage.
  !     LENIW  IWORK(18)  Length of IWORK actually required.  This is
  !                       defined on normal returns and on an illegal
  !                       input return for insufficient storage.
  !
  !     The following two arrays are segments of the RWORK array which may
  !     also be of interest to the user as optional outputs.  For each
  !     array, the table below gives its internal name, its base address
  !     in RWORK, and its description.
  !
  !     Name  Base address  Description
  !     ----  ------------  ----------------------------------------------
  !     YH    21            The Nordsieck history array, of size NYH by
  !                         (NQCUR + 1), where NYH is the initial value of
  !                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of
  !                         YH contains HCUR**j/factorial(j) times the jth
  !                         derivative of the interpolating polynomial
  !                         currently representing the solution, evaluated
  !                         at t = TCUR.
  !     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated
  !                         corrections on each step, scaled on output to
  !                         represent the estimated local error in Y on
  !                         the last step.  This is the vector e in the
  !                         description of the error control.  It is
  !                         defined only on successful return from DLSODE.
  !
  !
  !                    Part 2.  Other Callable Routines
  !                    --------------------------------
  !
  !     The following are optional calls which the user may make to gain
  !     additional capabilities in conjunction with DLSODE.
  !
  !     Form of call              Function
  !     ------------------------  ----------------------------------------
  !     CALL XSETUN(LUN)          Set the logical unit number, LUN, for
  !                               output of messages from DLSODE, if the
  !                               default is not desired.  The default
  !                               value of LUN is 6. This call may be made
  !                               at any time and will take effect
  !                               immediately.
  !     CALL XSETF(MFLAG)         Set a flag to control the printing of
  !                               messages by DLSODE.  MFLAG = 0 means do
  !                               not print.  (Danger:  this risks losing
  !                               valuable information.)  MFLAG = 1 means
  !                               print (the default).  This call may be
  !                               made at any time and will take effect
  !                               immediately.
  !     CALL DSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the
  !                               internal COMMON blocks used by DLSODE
  !                               (see Part 3 below).  RSAV must be a
  !                               real array of length 218 or more, and
  !                               ISAV must be an integer array of length
  !                               37 or more.  JOB = 1 means save COMMON
  !                               into RSAV/ISAV.  JOB = 2 means restore
  !                               COMMON from same.  DSRCOM is useful if
  !                               one is interrupting a run and restarting
  !                               later, or alternating between two or
  !                               more problems solved with DLSODE.
  !     CALL DINTDY(,,,,,)        Provide derivatives of y, of various
  !     (see below)               orders, at a specified point t, if
  !                               desired.  It may be called only after a
  !                               successful return from DLSODE.  Detailed
  !                               instructions follow.
  !
  !     Detailed instructions for using DINTDY
  !     --------------------------------------
  !     The form of the CALL is:
  !
  !           CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
  !
  !     The input parameters are:
  !
  !     T          Value of independent variable where answers are
  !                desired (normally the same as the T last returned by
  !                DLSODE).  For valid results, T must lie between
  !                TCUR - HU and TCUR.  (See "Optional Outputs" above
  !                for TCUR and HU.)
  !     K          Integer order of the derivative desired.  K must
  !                satisfy 0 <= K <= NQCUR, where NQCUR is the current
  !                order (see "Optional Outputs").  The capability
  !                corresponding to K = 0, i.e., computing y(t), is
  !                already provided by DLSODE directly.  Since
  !                NQCUR >= 1, the first derivative dy/dt is always
  !                available with DINTDY.
  !     RWORK(21)  The base address of the history array YH.
  !     NYH        Column length of YH, equal to the initial value of NEQ.
  !
  !     The output parameters are:
  !
  !     DKY        Real array of length NEQ containing the computed value
  !                of the Kth derivative of y(t).
  !     IFLAG      Integer flag, returned as 0 if K and T were legal,
  !                -1 if K was illegal, and -2 if T was illegal.
  !                On an error return, a message is also written.
  !
  !
  !                          Part 3.  Common Blocks
  !                          ----------------------
  !
  !     If DLSODE is to be used in an overlay situation, the user must
  !     declare, in the primary overlay, the variables in:
  !     (1) the call sequence to DLSODE,
  !     (2) the internal COMMON block /DLS001/, of length 255
  !         (218 double precision words followed by 37 integer words).
  !
  !     If DLSODE is used on a system in which the contents of internal
  !     COMMON blocks are not preserved between calls, the user should
  !     declare the above COMMON block in his main program to insure that
  !     its contents are preserved.
  !
  !     If the solution of a given problem by DLSODE is to be interrupted
  !     and then later continued, as when restarting an interrupted run or
  !     alternating between two or more problems, the user should save,
  !     following the return from the last DLSODE call prior to the
  !     interruption, the contents of the call sequence variables and the
  !     internal COMMON block, and later restore these values before the
  !     next DLSODE call for that problem.   In addition, if XSETUN and/or
  !     XSETF was called for non-default handling of error messages, then
  !     these calls must be repeated.  To save and restore the COMMON
  !     block, use subroutine DSRCOM (see Part 2 above).
  !
  !
  !              Part 4.  Optionally Replaceable Solver Routines
  !              -----------------------------------------------
  !
  !     Below are descriptions of two routines in the DLSODE package which
  !     relate to the measurement of errors.  Either routine can be
  !     replaced by a user-supplied version, if desired.  However, since
  !     such a replacement may have a major impact on performance, it
  !     should be done only when absolutely necessary, and only with great
  !     caution.  (Note:  The means by which the package version of a
  !     routine is superseded by the user's version may be system-
  !     dependent.)
  !
  !     DEWSET
  !     ------
  !     The following subroutine is called just before each internal
  !     integration step, and sets the array of error weights, EWT, as
  !     described under ITOL/RTOL/ATOL above:
  !
  !           SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
  !
  !     where NEQ, ITOL, RTOL, and ATOL are as in the DLSODE call
  !     sequence, YCUR contains the current dependent variable vector,
  !     and EWT is the array of weights set by DEWSET.
  !
  !     If the user supplies this subroutine, it must return in EWT(i)
  !     (i = 1,...,NEQ) a positive quantity suitable for comparing errors
  !     in Y(i) to.  The EWT array returned by DEWSET is passed to the
  !     DVNORM routine (see below), and also used by DLSODE in the
  !     computation of the optional output IMXER, the diagonal Jacobian
  !     approximation, and the increments for difference quotient
  !     Jacobians.
  !
  !     In the user-supplied version of DEWSET, it may be desirable to use
  !     the current values of derivatives of y. Derivatives up to order NQ
  !     are available from the history array YH, described above under
  !     optional outputs.  In DEWSET, YH is identical to the YCUR array,
  !     extended to NQ + 1 columns with a column length of NYH and scale
  !     factors of H**j/factorial(j).  On the first call for the problem,
  !     given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
  !     NYH is the initial value of NEQ.  The quantities NQ, H, and NST
  !     can be obtained by including in SEWSET the statements:
  !           DOUBLE PRECISION RLS
  !           COMMON /DLS001/ RLS(218),ILS(37)
  !           NQ = ILS(33)
  !           NST = ILS(34)
  !           H = RLS(212)
  !     Thus, for example, the current value of dy/dt can be obtained as
  !     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary
  !     when NST = 0).
  !
  !     DVNORM
  !     ------
  !     DVNORM is a real function routine which computes the weighted
  !     root-mean-square norm of a vector v:
  !
  !        d = DVNORM (n, v, w)
  !
  !     where:
  !     n = the length of the vector,
  !     v = real array of length n containing the vector,
  !     w = real array of length n containing weights,
  !     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ).
  !
  !     DVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where
  !     EWT is as set by subroutine DEWSET.
  !
  !     If the user supplies this function, it should return a nonnegative
  !     value of DVNORM suitable for use in the error control in DLSODE.
  !     None of the arguments should be altered by DVNORM.  For example, a
  !     user-supplied DVNORM routine might:
  !     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
  !     - Ignore some components of v in the norm, with the effect of
  !       suppressing the error control on those components of Y.
  !  ---------------------------------------------------------------------
  !***ROUTINES CALLED  DEWSET, DINTDY, DUMACH, DSTODE, DVNORM, XERRWD
  !***COMMON BLOCKS    DLS001
  !***REVISION HISTORY  (YYYYMMDD)
  ! 19791129  DATE WRITTEN
  ! 19791213  Minor changes to declarations; DELP init. in STODE.
  ! 19800118  Treat NEQ as array; integer declarations added throughout;
  !           minor changes to prologue.
  ! 19800306  Corrected TESCO(1,NQP1) setting in CFODE.
  ! 19800519  Corrected access of YH on forced order reduction;
  !           numerous corrections to prologues and other comments.
  ! 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
  !           minor corrections to main prologue.
  ! 19800923  Added zero initialization of HU and NQU.
  ! 19801218  Revised XERRWD routine; minor corrections to main prologue.
  ! 19810401  Minor changes to comments and an error message.
  ! 19810814  Numerous revisions: replaced EWT by 1/EWT; used flags
  !           JCUR, ICF, IERPJ, IERSL between STODE and subordinates;
  !           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
  !           reorganized returns from STODE; reorganized type decls.;
  !           fixed message length in XERRWD; changed default LUNIT to 6;
  !           changed Common lengths; changed comments throughout.
  ! 19870330  Major update by ACH: corrected comments throughout;
  !           removed TRET from Common; rewrote EWSET with 4 loops;
  !           fixed t test in INTDY; added Cray directives in STODE;
  !           in STODE, fixed DELP init. and logic around PJAC call;
  !           combined routines to save/restore Common;
  !           passed LEVEL = 0 in error message calls (except run abort).
  ! 19890426  Modified prologue to SLATEC/LDOC format.  (FNF)
  ! 19890501  Many improvements to prologue.  (FNF)
  ! 19890503  A few final corrections to prologue.  (FNF)
  ! 19890504  Minor cosmetic changes.  (FNF)
  ! 19890510  Corrected description of Y in Arguments section.  (FNF)
  ! 19890517  Minor corrections to prologue.  (FNF)
  ! 19920514  Updated with prologue edited 891025 by G. Shaw for manual.
  ! 19920515  Converted source lines to upper case.  (FNF)
  ! 19920603  Revised XERRWD calls using mixed upper-lower case.  (ACH)
  ! 19920616  Revised prologue comment regarding CFT.  (ACH)
  ! 19921116  Revised prologue comments regarding Common.  (ACH).
  ! 19930326  Added comment about non-reentrancy.  (FNF)
  ! 19930723  Changed D1MACH to DUMACH. (FNF)
  ! 19930801  Removed ILLIN and NTREP from Common (affects driver logic);
  !           minor changes to prologue and internal comments;
  !           changed Hollerith strings to quoted strings;
  !           changed internal comments to mixed case;
  !           replaced XERRWD with new version using character type;
  !           changed dummy dimensions from 1 to *. (ACH)
  ! 19930809  Changed to generic intrinsic names; changed names of
  !           subprograms and Common blocks to DLSODE etc. (ACH)
  ! 19930929  Eliminated use of REAL intrinsic; other minor changes. (ACH)
  ! 20010412  Removed all 'own' variables from Common block /DLS001/
  !           (affects declarations in 6 routines). (ACH)
  ! 20010509  Minor corrections to prologue. (ACH)
  ! 20031105  Restored 'own' variables to Common block /DLS001/, to
  !           enable interrupt/restart feature. (ACH)
  ! 20031112  Added SAVE statements for data-loaded constants.
  !
  !***END PROLOGUE  DLSODE
  !
  !*Internal Notes:
  !
  ! Other Routines in the DLSODE Package.
  !
  ! In addition to Subroutine DLSODE, the DLSODE package includes the
  ! following subroutines and function routines:
  !  DINTDY   computes an interpolated value of the y vector at t = TOUT.
  !  DSTODE   is the core integrator, which does one step of the
  !           integration and the associated error control.
  !  DCFODE   sets all method coefficients and test constants.
  !  DPREPJ   computes and preprocesses the Jacobian matrix J = df/dy
  !           and the Newton iteration matrix P = I - h*l0*J.
  !  DSOLSY   manages solution of linear system in chord iteration.
  !  DEWSET   sets the error weight vector EWT before each step.
  !  DVNORM   computes the weighted R.M.S. norm of a vector.
  !  DSRCOM   is a user-callable routine to save and restore
  !           the contents of the internal Common block.
  !  DGEFA and DGESL   are routines from LINPACK for solving full
  !           systems of linear algebraic equations.
  !  DGBFA and DGBSL   are routines from LINPACK for solving banded
  !           linear systems.
  !  DUMACH   computes the unit roundoff in a machine-independent manner.
  !  XERRWD, XSETUN, XSETF, IXSAV, IUMACH   handle the printing of all
  !           error messages and warnings.  XERRWD is machine-dependent.
  ! Note: DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
  ! All the others are subroutines.
  !
  !**End
  !
  !  Declare externals.
        !EXTERNAL DPREPJ , DSOLSY
        !DOUBLE PRECISION DUMACH , DVNORM
  !
  !  Declare all other variables.
        INTEGER INIt , MXStep , MXHnil , NHNil , NSLast , NYH , IOWns ,   &
       &        ICF , IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH ,   &
       &        LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,        &
       &        MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST , NFE ,     &
       &        NJE , NQU
        INTEGER i , i1 , i2 , iflag , imxer , kgo , lf0 , leniw , lenrw , &
       &        lenwm , ml , mord , mu , mxhnl0 , mxstp0
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION atoli , ayi , big , ewti , h0 , hmax , hmx , rh ,&
       &                 rtoli , tcrit , tdist , tnext , tol , tolsf ,    &
       &                 tp , size , sum , w0
        DIMENSION mord(2)
        LOGICAL ihit
        CHARACTER*80 msg
        SAVE mord , mxstp0 , mxhnl0
  !-----------------------------------------------------------------------
  ! The following internal Common block contains
  ! (a) variables which are local to any subroutine but whose values must
  !     be preserved between calls to the routine ("own" variables), and
  ! (b) variables which are communicated between subroutines.
  ! The block DLS001 is declared in subroutines DLSODE, DINTDY, DSTODE,
  ! DPREPJ, and DSOLSY.
  ! Groups of variables are replaced by dummy arrays in the Common
  ! declarations in routines where those variables are not used.
  !-----------------------------------------------------------------------
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , INIt , MXStep , MXHnil ,       &
       &                NHNil , NSLast , NYH , IOWns(6) , ICF , IERpj ,   &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
        DATA mord(1) , mord(2)/12 , 5/ , mxstp0/500/ , mxhnl0/10/
  !-----------------------------------------------------------------------
  ! Block A.
  ! This code block is executed on every call.
  ! It tests ISTATE and ITASK for legality and branches appropriately.
  ! If ISTATE .GT. 1 but the flag INIT shows that initialization has
  ! not yet been done, an error return occurs.
  ! If ISTATE = 1 and TOUT = T, return immediately.
  !-----------------------------------------------------------------------
  !
  !***FIRST EXECUTABLE STATEMENT  DLSODE
        IF ( Istate.LT.1 .OR. Istate.GT.3 ) THEN
  !-----------------------------------------------------------------------
  ! Block I.
  ! The following block handles all error returns due to illegal input
  ! (ISTATE = -3), as detected before calling the core integrator.
  ! First the error message routine is called.  If the illegal input
  ! is a negative ISTATE, the run is aborted (apparent infinite loop).
  !-----------------------------------------------------------------------
           msg = 'DLSODE-  ISTATE (=I1) illegal '
           CALL XERRWD(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
           IF ( Istate.GE.0 ) GOTO 2700
  !
           msg = 'DLSODE-  Run aborted.. apparent infinite loop     '
           CALL XERRWD(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
           GOTO 99999
        ELSE
           IF ( Itask.LT.1 .OR. Itask.GT.5 ) THEN
              msg = 'DLSODE-  ITASK (=I1) illegal  '
              CALL XERRWD(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
              GOTO 2700
           ELSE
              IF ( Istate.EQ.1 ) THEN
                 INIt = 0
                 IF ( Tout.EQ.T ) RETURN
              ELSEIF ( INIt.EQ.0 ) THEN
                 msg =                                                    &
       &              'DLSODE-  ISTATE .GT. 1 but DLSODE not initialized '
                 CALL XERRWD(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
                 GOTO 2700
              ELSEIF ( Istate.EQ.2 ) THEN
                 GOTO 50
              ENDIF
  !-----------------------------------------------------------------------
  ! Block B.
  ! The next code block is executed for the initial call (ISTATE = 1),
  ! or for a continuation call with parameter changes (ISTATE = 3).
  ! It contains checking of all inputs and various initializations.
  !
  ! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
  ! MF, ML, and MU.
  !-----------------------------------------------------------------------
              IF ( Neq(1).LE.0 ) THEN
                 msg = 'DLSODE-  NEQ (=I1) .LT. 1     '
                 CALL XERRWD(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
                 GOTO 2700
              ELSE
                 IF ( Istate.NE.1 ) THEN
                    IF ( Neq(1).GT.N ) THEN
                       msg =                                              &
       &              'DLSODE-  ISTATE = 3 and NEQ increased (I1 to I2)  '
                       CALL XERRWD(msg,50,5,0,2,N,Neq(1),0,0.0D0,0.0D0)
                       GOTO 2700
                    ENDIF
                 ENDIF
                 N = Neq(1)
                 IF ( Itol.LT.1 .OR. Itol.GT.4 ) THEN
                    msg = 'DLSODE-  ITOL (=I1) illegal   '
                    CALL XERRWD(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
                    GOTO 2700
                 ELSEIF ( Iopt.LT.0 .OR. Iopt.GT.1 ) THEN
                    msg = 'DLSODE-  IOPT (=I1) illegal   '
                    CALL XERRWD(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
                    GOTO 2700
                 ELSE
                    METh = Mf/10
                    MITer = Mf - 10*METh
                    IF ( METh.LT.1 .OR. METh.GT.2 ) GOTO 2000
                    IF ( MITer.LT.0 .OR. MITer.GT.5 ) GOTO 2000
                    IF ( MITer.GT.3 ) THEN
                       ml = Iwork(1)
                       mu = Iwork(2)
                       IF ( ml.LT.0 .OR. ml.GE.N ) THEN
                          msg =                                           &
       &              'DLSODE-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
                          CALL XERRWD(msg,50,9,0,2,ml,Neq(1),0,0.0D0,     &
       &                              0.0D0)
                          GOTO 2700
                       ELSEIF ( mu.LT.0 .OR. mu.GE.N ) THEN
                          msg =                                           &
       &              'DLSODE-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
                          CALL XERRWD(msg,50,10,0,2,mu,Neq(1),0,0.0D0,    &
       &                              0.0D0)
                          GOTO 2700
                       ENDIF
                    ENDIF
  ! Next process and check the optional inputs. --------------------------
                    IF ( Iopt.EQ.1 ) THEN
                       MAXord = Iwork(5)
                       IF ( MAXord.LT.0 ) THEN
                          msg = 'DLSODE-  MAXORD (=I1) .LT. 0  '
                          CALL XERRWD(msg,30,11,0,1,MAXord,0,0,0.0D0,     &
       &                              0.0D0)
                          GOTO 2700
                       ELSE
                          IF ( MAXord.EQ.0 ) MAXord = 100
                          MAXord = MIN(MAXord,mord(METh))
                          MXStep = Iwork(6)
                          IF ( MXStep.LT.0 ) THEN
                             msg = 'DLSODE-  MXSTEP (=I1) .LT. 0  '
                             CALL XERRWD(msg,30,12,0,1,MXStep,0,0,0.0D0,  &
       &                                 0.0D0)
                             GOTO 2700
                          ELSE
                             IF ( MXStep.EQ.0 ) MXStep = mxstp0
                             MXHnil = Iwork(7)
                             IF ( MXHnil.LT.0 ) THEN
                                msg = 'DLSODE-  MXHNIL (=I1) .LT. 0  '
                                CALL XERRWD(msg,30,13,0,1,MXHnil,0,0,     &
       &                           0.0D0,0.0D0)
                                GOTO 2700
                             ELSE
                                IF ( MXHnil.EQ.0 ) MXHnil = mxhnl0
                                IF ( Istate.EQ.1 ) THEN
                                   h0 = Rwork(5)
                                   IF ( (Tout-T)*h0.LT.0.0D0 ) THEN
                                      msg =                               &
       &                        'DLSODE-  TOUT (=R1) behind T (=R2)      '
                                      CALL XERRWD(msg,40,14,0,0,0,0,2,    &
       &                                 Tout,T)
                                      msg =                               &
       &              '      Integration direction is given by H0 (=R1)  '
                                      CALL XERRWD(msg,50,14,0,0,0,0,1,h0, &
       &                                 0.0D0)
                                      GOTO 2700
                                   ENDIF
                                ENDIF
                                hmax = Rwork(6)
                                IF ( hmax.LT.0.0D0 ) THEN
                                   msg = 'DLSODE-  HMAX (=R1) .LT. 0.0  '
                                   CALL XERRWD(msg,30,15,0,0,0,0,1,hmax,  &
       &                              0.0D0)
                                   GOTO 2700
                                ELSE
                                   HMXi = 0.0D0
                                   IF ( hmax.GT.0.0D0 ) HMXi = 1.0D0/hmax
                                   HMIn = Rwork(7)
                                   IF ( HMIn.LT.0.0D0 ) THEN
                                      msg =                               &
       &                                 'DLSODE-  HMIN (=R1) .LT. 0.0  '
                                      CALL XERRWD(msg,30,16,0,0,0,0,1,    &
       &                                 HMIn,0.0D0)
                                      GOTO 2700
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ELSE
                       MAXord = mord(METh)
                       MXStep = mxstp0
                       MXHnil = mxhnl0
                       IF ( Istate.EQ.1 ) h0 = 0.0D0
                       HMXi = 0.0D0
                       HMIn = 0.0D0
                    ENDIF
  !-----------------------------------------------------------------------
  ! Set work array pointers and check lengths LRW and LIW.
  ! Pointers to segments of RWORK and IWORK are named by prefixing L to
  ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
  ! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
  !-----------------------------------------------------------------------
                    LYH = 21
                    IF ( Istate.EQ.1 ) NYH = N
                    LWM = LYH + (MAXord+1)*NYH
                    IF ( MITer.EQ.0 ) lenwm = 0
                    IF ( MITer.EQ.1 .OR. MITer.EQ.2 ) lenwm = N*N + 2
                    IF ( MITer.EQ.3 ) lenwm = N + 2
                    IF ( MITer.GE.4 ) lenwm = (2*ml+mu+1)*N + 2
                    LEWt = LWM + lenwm
                    LSAvf = LEWt + N
                    LACor = LSAvf + N
                    lenrw = LACor + N - 1
                    Iwork(17) = lenrw
                    LIWm = 1
                    leniw = 20 + N
                    IF ( MITer.EQ.0 .OR. MITer.EQ.3 ) leniw = 20
                    Iwork(18) = leniw
                    IF ( lenrw.GT.Lrw ) THEN
                       msg =                                              &
       &    'DLSODE-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                       CALL XERRWD(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                       GOTO 2700
                    ELSEIF ( leniw.GT.Liw ) THEN
                       msg =                                              &
       &    'DLSODE-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                       CALL XERRWD(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                       GOTO 2700
                    ELSE
  ! Check RTOL and ATOL for legality. ------------------------------------
                       rtoli = Rtol(1)
                       atoli = Atol(1)
                       DO i = 1 , N
                          IF ( Itol.GE.3 ) rtoli = Rtol(i)
                          IF ( Itol.EQ.2 .OR. Itol.EQ.4 ) atoli = Atol(i)
                          IF ( rtoli.LT.0.0D0 ) GOTO 2100
                          IF ( atoli.LT.0.0D0 ) GOTO 2200
                       ENDDO
                       IF ( Istate.EQ.1 ) THEN
  !-----------------------------------------------------------------------
  ! Block C.
  ! The next block is for the initial call only (ISTATE = 1).
  ! It contains all remaining initializations, the initial call to F,
  ! and the calculation of the initial step size.
  ! The error weights in EWT are inverted after being loaded.
  !-----------------------------------------------------------------------
                          UROund = DUMACH()
                          TN = T
                          IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
                             tcrit = Rwork(1)
                             IF ( (tcrit-Tout)*(Tout-T).LT.0.0D0 )        &
       &                          GOTO 2500
                             IF ( h0.NE.0.0D0 .AND. (T+h0-tcrit)          &
       &                          *h0.GT.0.0D0 ) h0 = tcrit - T
                          ENDIF
                          JSTart = 0
                          IF ( MITer.GT.0 ) Rwork(LWM) = SQRT(UROund)
                          NHNil = 0
                          NST = 0
                          NJE = 0
                          NSLast = 0
                          HU = 0.0D0
                          NQU = 0
                          CCMax = 0.3D0
                          MAXcor = 3
                          MSBp = 20
                          MXNcf = 10
  ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
                          lf0 = LYH + NYH
                          CALL F(Neq,T,Y,Rwork(lf0))
                          NFE = 1
  ! Load the initial value vector in YH. ---------------------------------
                          DO i = 1 , N
                             Rwork(i+LYH-1) = Y(i)
                          ENDDO
  ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                          NQ = 1
                          H = 1.0D0
                          CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),        &
       &                              Rwork(LEWt))
                          DO i = 1 , N
                             IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 2300
                             Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
                          ENDDO
  !-----------------------------------------------------------------------
  ! The coding below computes the step size, H0, to be attempted on the
  ! first step, unless the user has supplied a value for this.
  ! First check that TOUT - T differs significantly from zero.
  ! A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
  ! if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
  ! so as to be between 100*UROUND and 1.0E-3.
  ! Then the computed value H0 is given by..
  !                                      NEQ
  !   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
  !                                       1
  ! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
  !         f(i)   = i-th component of initial value of f,
  !         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
  ! The sign of H0 is inferred from the initial values of TOUT and T.
  !-----------------------------------------------------------------------
                          IF ( h0.EQ.0.0D0 ) THEN
                             tdist = ABS(Tout-T)
                             w0 = MAX(ABS(T),ABS(Tout))
                             IF ( tdist.LT.2.0D0*UROund*w0 ) THEN
                                msg =                                     &
       &    'DLSODE-  TOUT (=R1) too close to T(=R2) to start integration'
                                CALL XERRWD(msg,60,22,0,0,0,0,2,Tout,T)
                                GOTO 2700
                             ELSE
                                tol = Rtol(1)
                                IF ( Itol.GT.2 ) THEN
                                   DO i = 1 , N
                                      tol = MAX(tol,Rtol(i))
                                   ENDDO
                                ENDIF
                                IF ( tol.LE.0.0D0 ) THEN
                                   atoli = Atol(1)
                                   DO i = 1 , N
                                      IF ( Itol.EQ.2 .OR. Itol.EQ.4 )     &
       &                                 atoli = Atol(i)
                                      ayi = ABS(Y(i))
                                      IF ( ayi.NE.0.0D0 )                 &
       &                                 tol = MAX(tol,atoli/ayi)
                                   ENDDO
                                ENDIF
                                tol = MAX(tol,100.0D0*UROund)
                                tol = MIN(tol,0.001D0)
                                sum = DVNORM(N,Rwork(lf0),Rwork(LEWt))
                                sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                                h0 = 1.0D0/SQRT(sum)
                                h0 = MIN(h0,tdist)
                                h0 = SIGN(h0,Tout-T)
                             ENDIF
                          ENDIF
  ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                          rh = ABS(h0)*HMXi
                          IF ( rh.GT.1.0D0 ) h0 = h0/rh
  ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                          H = h0
                          DO i = 1 , N
                             Rwork(i+lf0-1) = h0*Rwork(i+lf0-1)
                          ENDDO
                          GOTO 700
                       ELSE
  ! If ISTATE = 3, set flag to signal parameter changes to DSTODE. -------
                          JSTart = -1
                          IF ( NQ.GT.MAXord ) THEN
  ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
                             DO i = 1 , N
                                Rwork(i+LSAvf-1) = Rwork(i+LWM-1)
                             ENDDO
                          ENDIF
  ! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
                          IF ( MITer.GT.0 ) Rwork(LWM) = SQRT(UROund)
                          IF ( N.NE.NYH ) THEN
  ! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                             i1 = LYH + L*NYH
                             i2 = LYH + (MAXord+1)*NYH - 1
                             IF ( i1.LE.i2 ) THEN
                                DO i = i1 , i2
                                   Rwork(i) = 0.0D0
                                ENDDO
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  ! Block D.
  ! The next code block is for continuation calls only (ISTATE = 2 or 3)
  ! and is to check stop conditions before taking a step.
  !-----------------------------------------------------------------------
   50      NSLast = NST
           GOTO (100,600,200,300,400) , Itask
        ENDIF
   100  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2600
        T = Tout
        GOTO 1400
   200  tp = TN - HU*(1.0D0+100.0D0*UROund)
        IF ( (tp-Tout)*H.GT.0.0D0 ) THEN
           msg =                                                          &
       &    'DLSODE-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
           CALL XERRWD(msg,60,23,0,1,Itask,0,2,Tout,tp)
           GOTO 2700
        ELSE
           IF ( (TN-Tout)*H.GE.0.0D0 ) GOTO 1300
           GOTO 600
        ENDIF
   300  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2400
        IF ( (tcrit-Tout)*H.LT.0.0D0 ) GOTO 2500
        IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 500
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2600
        T = Tout
        GOTO 1400
   400  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2400
   500  hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
        IF ( ihit ) GOTO 1300
        tnext = TN + H*(1.0D0+4.0D0*UROund)
        IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
           H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
           IF ( Istate.EQ.2 ) JSTart = -2
        ENDIF
  !-----------------------------------------------------------------------
  ! Block E.
  ! The next block is normally executed for all calls and contains
  ! the call to the one-step core integrator DSTODE.
  !
  ! This is a looping point for the integration steps.
  !
  ! First check for too many steps being taken, update EWT (if not at
  ! start of problem), check for too much accuracy being requested, and
  ! check for H below the roundoff level in T.
  !-----------------------------------------------------------------------
   600  IF ( (NST-NSLast).GE.MXStep ) THEN
  !-----------------------------------------------------------------------
  ! Block H.
  ! The following block handles all unsuccessful returns other than
  ! those for illegal input.  First the error message routine is called.
  ! If there was an error test or convergence test failure, IMXER is set.
  ! Then Y is loaded from YH and T is set to TN.  The optional outputs
  ! are loaded into the work arrays before returning.
  !-----------------------------------------------------------------------
  ! The maximum number of steps was taken before reaching TOUT. ----------
           msg = 'DLSODE-  At current T (=R1), MXSTEP (=I1) steps   '
           CALL XERRWD(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      taken on this call before reaching TOUT     '
           CALL XERRWD(msg,50,201,0,1,MXStep,0,1,TN,0.0D0)
           Istate = -1
           GOTO 1900
        ELSE
           CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
           DO i = 1 , N
              IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 1500
              Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
           ENDDO
        ENDIF
   700  tolsf = UROund*DVNORM(N,Rwork(LYH),Rwork(LEWt))
        IF ( tolsf.LE.1.0D0 ) THEN
           IF ( (TN+H).EQ.TN ) THEN
              NHNil = NHNil + 1
              IF ( NHNil.LE.MXHnil ) THEN
                 msg =                                                    &
       &              'DLSODE-  Warning..internal T (=R1) and H (=R2) are'
                 CALL XERRWD(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &    '      such that in the machine, T + H = T on the next step  '
                 CALL XERRWD(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &              '      (H = step size). Solver will continue anyway'
                 CALL XERRWD(msg,50,101,0,0,0,0,2,TN,H)
                 IF ( NHNil.GE.MXHnil ) THEN
                    msg =                                                 &
       &              'DLSODE-  Above warning has been issued I1 times.  '
                    CALL XERRWD(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
                    msg =                                                 &
       &              '      It will not be issued again for this problem'
                    CALL XERRWD(msg,50,102,0,1,MXHnil,0,0,0.0D0,0.0D0)
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  !  CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPREPJ,DSOLSY)
  !-----------------------------------------------------------------------
           CALL DSTODE(Neq,Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),       &
       &               Rwork(LSAvf),Rwork(LACor),Rwork(LWM),Iwork(LIWm),F,&
       &               JAC,DPREPJ,DSOLSY)
           kgo = 1 - KFLag
           GOTO (800,1600,1700) , kgo
        ELSE
           tolsf = tolsf*2.0D0
           IF ( NST.EQ.0 ) THEN
              msg = 'DLSODE-  At start of problem, too much accuracy   '
              CALL XERRWD(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
              msg =                                                       &
       &    '      requested for precision of machine..  See TOLSF (=R1) '
              CALL XERRWD(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
              Rwork(14) = tolsf
              GOTO 2700
           ELSE
  ! Too much accuracy requested for machine precision. -------------------
              msg = 'DLSODE-  At T (=R1), too much accuracy requested  '
              CALL XERRWD(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
              msg = '      for precision of machine..  see TOLSF (=R2) '
              CALL XERRWD(msg,50,203,0,0,0,0,2,TN,tolsf)
              Rwork(14) = tolsf
              Istate = -2
              GOTO 1900
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! Block F.
  ! The following block handles the case of a successful return from the
  ! core integrator (KFLAG = 0).  Test for stop conditions.
  !-----------------------------------------------------------------------
   800  INIt = 1
        GOTO (900,1300,1000,1100,1200) , Itask
  ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   900  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        T = Tout
        GOTO 1400
  ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   1000 IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        GOTO 1300
  ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
   1100 IF ( (TN-Tout)*H.LT.0.0D0 ) THEN
           hmx = ABS(TN) + ABS(H)
           ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
           IF ( ihit ) GOTO 1300
           tnext = TN + H*(1.0D0+4.0D0*UROund)
           IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
              H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
              JSTart = -2
           ENDIF
           GOTO 600
        ELSE
           CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
           T = Tout
           GOTO 1400
        ENDIF
  ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
   1200 hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
  !-----------------------------------------------------------------------
  ! Block G.
  ! The following block handles all successful returns from DLSODE.
  ! If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
  ! ISTATE is set to 2, and the optional outputs are loaded into the
  ! work arrays before returning.
  !-----------------------------------------------------------------------
   1300 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
           IF ( ihit ) T = tcrit
        ENDIF
   1400 Istate = 2
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
  ! EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
   1500 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODE-  At T (=R1), EWT(I1) has become R2 .LE. 0.'
        CALL XERRWD(msg,50,202,0,1,i,0,2,TN,ewti)
        Istate = -6
        GOTO 1900
  ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   1600 msg = 'DLSODE-  At T(=R1) and step size H(=R2), the error'
        CALL XERRWD(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      test failed repeatedly or with ABS(H) = HMIN'
        CALL XERRWD(msg,50,204,0,0,0,0,2,TN,H)
        Istate = -4
        GOTO 1800
  ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
   1700 msg = 'DLSODE-  At T (=R1) and step size H (=R2), the    '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      corrector convergence failed repeatedly     '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      or with ABS(H) = HMIN   '
        CALL XERRWD(msg,30,205,0,0,0,0,2,TN,H)
        Istate = -5
  ! Compute IMXER if relevant. -------------------------------------------
   1800 big = 0.0D0
        imxer = 1
        DO i = 1 , N
           size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
           IF ( big.LT.size ) THEN
              big = size
              imxer = i
           ENDIF
        ENDDO
        Iwork(16) = imxer
  ! Set Y vector, T, and optional outputs. -------------------------------
   1900 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
   2000 msg = 'DLSODE-  MF (=I1) illegal     '
        CALL XERRWD(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
        GOTO 2700
   2100 msg = 'DLSODE-  RTOL(I1) is R1 .LT. 0.0        '
        CALL XERRWD(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
        GOTO 2700
   2200 msg = 'DLSODE-  ATOL(I1) is R1 .LT. 0.0        '
        CALL XERRWD(msg,40,20,0,1,i,0,1,atoli,0.0D0)
        GOTO 2700
   2300 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODE-  EWT(I1) is R1 .LE. 0.0         '
        CALL XERRWD(msg,40,21,0,1,i,0,1,ewti,0.0D0)
        GOTO 2700
   2400 msg =                                                             &
       &    'DLSODE-  ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   '
        CALL XERRWD(msg,60,24,0,0,0,0,2,tcrit,TN)
        GOTO 2700
   2500 msg =                                                             &
       &    'DLSODE-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
        CALL XERRWD(msg,60,25,0,0,0,0,2,tcrit,Tout)
        GOTO 2700
   2600 msg = 'DLSODE-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
        CALL XERRWD(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
  !
   2700 Istate = -3
        RETURN
  !----------------------- END OF SUBROUTINE DLSODE ----------------------
  99999 END
  !*==DLSODA.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
        SUBROUTINE DLSODA(F,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt,&
       &                  Rwork,Lrw,Iwork,Liw,JAC,Jt)
        IMPLICIT NONE
  !*--DLSODA1856
  !*** Start of declarations inserted by SPAG
        !REAL DPRJA , DSOLSY
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC
        INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
       &        Jt
        DOUBLE PRECISION Y , T , Tout , Rtol , Atol , Rwork
        DIMENSION Neq(*) , Y(*) , Rtol(*) , Atol(*) , Rwork(Lrw) ,        &
       &          Iwork(Liw)
  !-----------------------------------------------------------------------
  ! This is the 12 November 2003 version of
  ! DLSODA: Livermore Solver for Ordinary Differential Equations, with
  !         Automatic method switching for stiff and nonstiff problems.
  !
  ! This version is in double precision.
  !
  ! DLSODA solves the initial value problem for stiff or nonstiff
  ! systems of first order ODEs,
  !     dy/dt = f(t,y) ,  or, in component form,
  !     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
  !
  ! This a variant version of the DLSODE package.
  ! It switches automatically between stiff and nonstiff methods.
  ! This means that the user does not have to determine whether the
  ! problem is stiff or not, and the solver will automatically choose the
  ! appropriate method.  It always starts with the nonstiff method.
  !
  ! Authors:       Alan C. Hindmarsh
  !                Center for Applied Scientific Computing, L-561
  !                Lawrence Livermore National Laboratory
  !                Livermore, CA 94551
  ! and
  !                Linda R. Petzold
  !                Univ. of California at Santa Barbara
  !                Dept. of Computer Science
  !                Santa Barbara, CA 93106
  !
  ! References:
  ! 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
  !     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  !     North-Holland, Amsterdam, 1983, pp. 55-64.
  ! 2.  Linda R. Petzold, Automatic Selection of Methods for Solving
  !     Stiff and Nonstiff Systems of Ordinary Differential Equations,
  !     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
  !-----------------------------------------------------------------------
  ! Summary of Usage.
  !
  ! Communication between the user and the DLSODA package, for normal
  ! situations, is summarized here.  This summary describes only a subset
  ! of the full set of options available.  See the full description for
  ! details, including alternative treatment of the Jacobian matrix,
  ! optional inputs and outputs, nonstandard options, and
  ! instructions for special situations.  See also the example
  ! problem (with program and output) following this summary.
  !
  ! A. First provide a subroutine of the form:
  !               SUBROUTINE F (NEQ, T, Y, YDOT)
  !               DOUBLE PRECISION T, Y(*), YDOT(*)
  ! which supplies the vector function f by loading YDOT(i) with f(i).
  !
  ! B. Write a main program which calls Subroutine DLSODA once for
  ! each point at which answers are desired.  This should also provide
  ! for possible use of logical unit 6 for output of error messages
  ! by DLSODA.  On the first call to DLSODA, supply arguments as follows:
  ! F      = name of subroutine for right-hand side vector f.
  !          This name must be declared External in calling program.
  ! NEQ    = number of first order ODEs.
  ! Y      = array of initial values, of length NEQ.
  ! T      = the initial value of the independent variable.
  ! TOUT   = first point where output is desired (.ne. T).
  ! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  ! RTOL   = relative tolerance parameter (scalar).
  ! ATOL   = absolute tolerance parameter (scalar or array).
  !          the estimated local error in y(i) will be controlled so as
  !          to be less than
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
  !          Thus the local error test passes if, in each component,
  !          either the absolute error is less than ATOL (or ATOL(i)),
  !          or the relative error is less than RTOL.
  !          Use RTOL = 0.0 for pure absolute error control, and
  !          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
  !          control.  Caution: actual (global) errors may exceed these
  !          local tolerances, so choose them conservatively.
  ! ITASK  = 1 for normal computation of output values of y at t = TOUT.
  ! ISTATE = integer flag (input and output).  Set ISTATE = 1.
  ! IOPT   = 0 to indicate no optional inputs used.
  ! RWORK  = real work array of length at least:
  !             22 + NEQ * MAX(16, NEQ + 9).
  !          See also Paragraph E below.
  ! LRW    = declared length of RWORK (in user's dimension).
  ! IWORK  = integer work array of length at least  20 + NEQ.
  ! LIW    = declared length of IWORK (in user's dimension).
  ! JAC    = name of subroutine for Jacobian matrix.
  !          Use a dummy name.  See also Paragraph E below.
  ! JT     = Jacobian type indicator.  Set JT = 2.
  !          See also Paragraph E below.
  ! Note that the main program must declare arrays Y, RWORK, IWORK,
  ! and possibly ATOL.
  !
  ! C. The output from the first call (or any call) is:
  !      Y = array of computed values of y(t) vector.
  !      T = corresponding value of independent variable (normally TOUT).
  ! ISTATE = 2  if DLSODA was successful, negative otherwise.
  !          -1 means excess work done on this call (perhaps wrong JT).
  !          -2 means excess accuracy requested (tolerances too small).
  !          -3 means illegal input detected (see printed message).
  !          -4 means repeated error test failures (check all inputs).
  !          -5 means repeated convergence failures (perhaps bad Jacobian
  !             supplied or wrong choice of JT or tolerances).
  !          -6 means error weight became zero during problem. (Solution
  !             component i vanished, and ATOL or ATOL(i) = 0.)
  !          -7 means work space insufficient to finish (see messages).
  !
  ! D. To continue the integration after a successful return, simply
  ! reset TOUT and call DLSODA again.  No other parameters need be reset.
  !
  ! E. Note: If and when DLSODA regards the problem as stiff, and
  ! switches methods accordingly, it must make use of the NEQ by NEQ
  ! Jacobian matrix, J = df/dy.  For the sake of simplicity, the
  ! inputs to DLSODA recommended in Paragraph B above cause DLSODA to
  ! treat J as a full matrix, and to approximate it internally by
  ! difference quotients.  Alternatively, J can be treated as a band
  ! matrix (with great potential reduction in the size of the RWORK
  ! array).  Also, in either the full or banded case, the user can supply
  ! J in closed form, with a routine whose name is passed as the JAC
  ! argument.  These alternatives are described in the paragraphs on
  ! RWORK, JAC, and JT in the full description of the call sequence below.
  !
  !-----------------------------------------------------------------------
  ! Example Problem.
  !
  ! The following is a simple example problem, with the coding
  ! needed for its solution by DLSODA.  The problem is from chemical
  ! kinetics, and consists of the following three rate equations:
  !     dy1/dt = -.04*y1 + 1.e4*y2*y3
  !     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
  !     dy3/dt = 3.e7*y2**2
  ! on the interval from t = 0.0 to t = 4.e10, with initial conditions
  ! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
  !
  ! The following coding solves this problem with DLSODA,
  ! printing results at t = .4, 4., ..., 4.e10.  It uses
  ! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
  ! y2 has much smaller values.
  ! At the end of the run, statistical quantities of interest are
  ! printed (see optional outputs in the full description below).
  !
  !     EXTERNAL FEX
  !     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
  !     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23)
  !     NEQ = 3
  !     Y(1) = 1.
  !     Y(2) = 0.
  !     Y(3) = 0.
  !     T = 0.
  !     TOUT = .4
  !     ITOL = 2
  !     RTOL = 1.D-4
  !     ATOL(1) = 1.D-6
  !     ATOL(2) = 1.D-10
  !     ATOL(3) = 1.D-6
  !     ITASK = 1
  !     ISTATE = 1
  !     IOPT = 0
  !     LRW = 70
  !     LIW = 23
  !     JT = 2
  !     DO 40 IOUT = 1,12
  !       CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
  !    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
  !       WRITE(6,20)T,Y(1),Y(2),Y(3)
  ! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
  !       IF (ISTATE .LT. 0) GO TO 80
  ! 40    TOUT = TOUT*10.
  !     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
  ! 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/
  !    1   ' Method last used =',I2,'   Last switch was at t =',D12.4)
  !     STOP
  ! 80  WRITE(6,90)ISTATE
  ! 90  FORMAT(///' Error halt.. ISTATE =',I3)
  !     STOP
  !     END
  !
  !     SUBROUTINE FEX (NEQ, T, Y, YDOT)
  !     DOUBLE PRECISION T, Y, YDOT
  !     DIMENSION Y(3), YDOT(3)
  !     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
  !     YDOT(3) = 3.D7*Y(2)*Y(2)
  !     YDOT(2) = -YDOT(1) - YDOT(3)
  !     RETURN
  !     END
  !
  ! The output of this program (on a CDC-7600 in single precision)
  ! is as follows:
  !
  !   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
  !   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02
  !   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01
  !   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01
  !   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01
  !   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01
  !   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01
  !   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01
  !   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01
  !   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01
  !   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01
  !   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00
  !
  !   No. steps = 361  No. f-s = 693  No. J-s =  64
  !   Method last used = 2   Last switch was at t =  6.0092e-03
  !-----------------------------------------------------------------------
  ! Full description of user interface to DLSODA.
  !
  ! The user interface to DLSODA consists of the following parts.
  !
  ! 1.   The call sequence to Subroutine DLSODA, which is a driver
  !      routine for the solver.  This includes descriptions of both
  !      the call sequence arguments and of user-supplied routines.
  !      following these descriptions is a description of
  !      optional inputs available through the call sequence, and then
  !      a description of optional outputs (in the work arrays).
  !
  ! 2.   Descriptions of other routines in the DLSODA package that may be
  !      (optionally) called by the user.  These provide the ability to
  !      alter error message handling, save and restore the internal
  !      Common, and obtain specified derivatives of the solution y(t).
  !
  ! 3.   Descriptions of Common blocks to be declared in overlay
  !      or similar environments, or to be saved when doing an interrupt
  !      of the problem and continued solution later.
  !
  ! 4.   Description of a subroutine in the DLSODA package,
  !      which the user may replace with his/her own version, if desired.
  !      this relates to the measurement of errors.
  !
  !-----------------------------------------------------------------------
  ! Part 1.  Call Sequence.
  !
  ! The call sequence parameters used for input only are
  !     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,
  ! and those used for both input and output are
  !     Y, T, ISTATE.
  ! The work arrays RWORK and IWORK are also used for conditional and
  ! optional inputs and optional outputs.  (The term output here refers
  ! to the return from Subroutine DLSODA to the user's calling program.)
  !
  ! The legality of input parameters will be thoroughly checked on the
  ! initial call for the problem, but not checked thereafter unless a
  ! change in input parameters is flagged by ISTATE = 3 on input.
  !
  ! The descriptions of the call arguments are as follows.
  !
  ! F      = the name of the user-supplied subroutine defining the
  !          ODE system.  The system must be put in the first-order
  !          form dy/dt = f(t,y), where f is a vector-valued function
  !          of the scalar t and the vector y.  Subroutine F is to
  !          compute the function f.  It is to have the form
  !               SUBROUTINE F (NEQ, T, Y, YDOT)
  !               DOUBLE PRECISION T, Y(*), YDOT(*)
  !          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
  !          is output.  Y and YDOT are arrays of length NEQ.
  !          Subroutine F should not alter Y(1),...,Y(NEQ).
  !          F must be declared External in the calling program.
  !
  !          Subroutine F may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in F) and/or Y has length exceeding NEQ(1).
  !          See the descriptions of NEQ and Y below.
  !
  !          If quantities computed in the F routine are needed
  !          externally to DLSODA, an extra call to F should be made
  !          for this purpose, for consistent and accurate results.
  !          If only the derivative dy/dt is needed, use DINTDY instead.
  !
  ! NEQ    = the size of the ODE system (number of first order
  !          ordinary differential equations).  Used only for input.
  !          NEQ may be decreased, but not increased, during the problem.
  !          If NEQ is decreased (with ISTATE = 3 on input), the
  !          remaining components of Y should be left undisturbed, if
  !          these are to be accessed in F and/or JAC.
  !
  !          Normally, NEQ is a scalar, and it is generally referred to
  !          as a scalar in this user interface description.  However,
  !          NEQ may be an array, with NEQ(1) set to the system size.
  !          (The DLSODA package accesses only NEQ(1).)  In either case,
  !          this parameter is passed as the NEQ argument in all calls
  !          to F and JAC.  Hence, if it is an array, locations
  !          NEQ(2),... may be used to store other integer data and pass
  !          it to F and/or JAC.  Subroutines F and/or JAC must include
  !          NEQ in a Dimension statement in that case.
  !
  ! Y      = a real array for the vector of dependent variables, of
  !          length NEQ or more.  Used for both input and output on the
  !          first call (ISTATE = 1), and only for output on other calls.
  !          On the first call, Y must contain the vector of initial
  !          values.  On output, Y contains the computed solution vector,
  !          evaluated at T.  If desired, the Y array may be used
  !          for other purposes between calls to the solver.
  !
  !          This array is passed as the Y argument in all calls to
  !          F and JAC.  Hence its length may exceed NEQ, and locations
  !          Y(NEQ+1),... may be used to store other real data and
  !          pass it to F and/or JAC.  (The DLSODA package accesses only
  !          Y(1),...,Y(NEQ).)
  !
  ! T      = the independent variable.  On input, T is used only on the
  !          first call, as the initial point of the integration.
  !          on output, after each call, T is the value at which a
  !          computed solution Y is evaluated (usually the same as TOUT).
  !          on an error return, T is the farthest point reached.
  !
  ! TOUT   = the next value of t at which a computed solution is desired.
  !          Used only for input.
  !
  !          When starting the problem (ISTATE = 1), TOUT may be equal
  !          to T for one call, then should .ne. T for the next call.
  !          For the initial t, an input value of TOUT .ne. T is used
  !          in order to determine the direction of the integration
  !          (i.e. the algebraic sign of the step sizes) and the rough
  !          scale of the problem.  Integration in either direction
  !          (forward or backward in t) is permitted.
  !
  !          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
  !          the first call (i.e. the first call with TOUT .ne. T).
  !          Otherwise, TOUT is required on every call.
  !
  !          If ITASK = 1, 3, or 4, the values of TOUT need not be
  !          monotone, but a value of TOUT which backs up is limited
  !          to the current internal T interval, whose endpoints are
  !          TCUR - HU and TCUR (see optional outputs, below, for
  !          TCUR and HU).
  !
  ! ITOL   = an indicator for the type of error control.  See
  !          description below under ATOL.  Used only for input.
  !
  ! RTOL   = a relative error tolerance parameter, either a scalar or
  !          an array of length NEQ.  See description below under ATOL.
  !          Input only.
  !
  ! ATOL   = an absolute error tolerance parameter, either a scalar or
  !          an array of length NEQ.  Input only.
  !
  !             The input parameters ITOL, RTOL, and ATOL determine
  !          the error control performed by the solver.  The solver will
  !          control the vector E = (E(i)) of estimated local errors
  !          in y, according to an inequality of the form
  !                      max-norm of ( E(i)/EWT(i) )   .le.   1,
  !          where EWT = (EWT(i)) is a vector of positive error weights.
  !          The values of RTOL and ATOL should all be non-negative.
  !          The following table gives the types (scalar/array) of
  !          RTOL and ATOL, and the corresponding form of EWT(i).
  !
  !             ITOL    RTOL       ATOL          EWT(i)
  !              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  !              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  !              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  !              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
  !
  !          When either of these parameters is a scalar, it need not
  !          be dimensioned in the user's calling program.
  !
  !          If none of the above choices (with ITOL, RTOL, and ATOL
  !          fixed throughout the problem) is suitable, more general
  !          error controls can be obtained by substituting a
  !          user-supplied routine for the setting of EWT.
  !          See Part 4 below.
  !
  !          If global errors are to be estimated by making a repeated
  !          run on the same problem with smaller tolerances, then all
  !          components of RTOL and ATOL (i.e. of EWT) should be scaled
  !          down uniformly.
  !
  ! ITASK  = an index specifying the task to be performed.
  !          Input only.  ITASK has the following values and meanings.
  !          1  means normal computation of output values of y(t) at
  !             t = TOUT (by overshooting and interpolating).
  !          2  means take one step only and return.
  !          3  means stop at the first internal mesh point at or
  !             beyond t = TOUT and return.
  !          4  means normal computation of output values of y(t) at
  !             t = TOUT but without overshooting t = TCRIT.
  !             TCRIT must be input as RWORK(1).  TCRIT may be equal to
  !             or beyond TOUT, but not behind it in the direction of
  !             integration.  This option is useful if the problem
  !             has a singularity at or beyond t = TCRIT.
  !          5  means take one step, without passing TCRIT, and return.
  !             TCRIT must be input as RWORK(1).
  !
  !          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
  !          (within roundoff), it will return T = TCRIT (exactly) to
  !          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
  !          in which case answers at t = TOUT are returned first).
  !
  ! ISTATE = an index used for input and output to specify the
  !          the state of the calculation.
  !
  !          On input, the values of ISTATE are as follows.
  !          1  means this is the first call for the problem
  !             (initializations will be done).  See note below.
  !          2  means this is not the first call, and the calculation
  !             is to continue normally, with no change in any input
  !             parameters except possibly TOUT and ITASK.
  !             (If ITOL, RTOL, and/or ATOL are changed between calls
  !             with ISTATE = 2, the new values will be used but not
  !             tested for legality.)
  !          3  means this is not the first call, and the
  !             calculation is to continue normally, but with
  !             a change in input parameters other than
  !             TOUT and ITASK.  Changes are allowed in
  !             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
  !             and any optional inputs except H0, MXORDN, and MXORDS.
  !             (See IWORK description for ML and MU.)
  !          Note:  A preliminary call with TOUT = T is not counted
  !          as a first call here, as no initialization or checking of
  !          input is done.  (Such a call is sometimes useful for the
  !          purpose of outputting the initial conditions.)
  !          Thus the first call for which TOUT .ne. T requires
  !          ISTATE = 1 on input.
  !
  !          On output, ISTATE has the following values and meanings.
  !           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
  !           2  means the integration was performed successfully.
  !          -1  means an excessive amount of work (more than MXSTEP
  !              steps) was done on this call, before completing the
  !              requested task, but the integration was otherwise
  !              successful as far as T.  (MXSTEP is an optional input
  !              and is normally 500.)  To continue, the user may
  !              simply reset ISTATE to a value .gt. 1 and call again
  !              (the excess work step counter will be reset to 0).
  !              In addition, the user may increase MXSTEP to avoid
  !              this error return (see below on optional inputs).
  !          -2  means too much accuracy was requested for the precision
  !              of the machine being used.  This was detected before
  !              completing the requested task, but the integration
  !              was successful as far as T.  To continue, the tolerance
  !              parameters must be reset, and ISTATE must be set
  !              to 3.  The optional output TOLSF may be used for this
  !              purpose.  (Note: If this condition is detected before
  !              taking any steps, then an illegal input return
  !              (ISTATE = -3) occurs instead.)
  !          -3  means illegal input was detected, before taking any
  !              integration steps.  See written message for details.
  !              Note:  If the solver detects an infinite loop of calls
  !              to the solver with illegal input, it will cause
  !              the run to stop.
  !          -4  means there were repeated error test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              The problem may have a singularity, or the input
  !              may be inappropriate.
  !          -5  means there were repeated convergence test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              This may be caused by an inaccurate Jacobian matrix,
  !              if one is being used.
  !          -6  means EWT(i) became zero for some i during the
  !              integration.  Pure relative error control (ATOL(i)=0.0)
  !              was requested on a variable which has now vanished.
  !              The integration was successful as far as T.
  !          -7  means the length of RWORK and/or IWORK was too small to
  !              proceed, but the integration was successful as far as T.
  !              This happens when DLSODA chooses to switch methods
  !              but LRW and/or LIW is too small for the new method.
  !
  !          Note:  Since the normal output value of ISTATE is 2,
  !          it does not need to be reset for normal continuation.
  !          Also, since a negative input value of ISTATE will be
  !          regarded as illegal, a negative output value requires the
  !          user to change it, and possibly other inputs, before
  !          calling the solver again.
  !
  ! IOPT   = an integer flag to specify whether or not any optional
  !          inputs are being used on this call.  Input only.
  !          The optional inputs are listed separately below.
  !          IOPT = 0 means no optional inputs are being used.
  !                   default values will be used in all cases.
  !          IOPT = 1 means one or more optional inputs are being used.
  !
  ! RWORK  = a real array (double precision) for work space, and (in the
  !          first 20 words) for conditional and optional inputs and
  !          optional outputs.
  !          As DLSODA switches automatically between stiff and nonstiff
  !          methods, the required length of RWORK can change during the
  !          problem.  Thus the RWORK array passed to DLSODA can either
  !          have a static (fixed) length large enough for both methods,
  !          or have a dynamic (changing) length altered by the calling
  !          program in response to output from DLSODA.
  !
  !                       --- Fixed Length Case ---
  !          If the RWORK length is to be fixed, it should be at least
  !               MAX (LRN, LRS),
  !          where LRN and LRS are the RWORK lengths required when the
  !          current method is nonstiff or stiff, respectively.
  !
  !          The separate RWORK length requirements LRN and LRS are
  !          as follows:
  !          IF NEQ is constant and the maximum method orders have
  !          their default values, then
  !             LRN = 20 + 16*NEQ,
  !             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
  !             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
  !          Under any other conditions, LRN and LRS are given by:
  !             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,
  !             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,
  !          where
  !             NYH    = the initial value of NEQ,
  !             MXORDN = 12, unless a smaller value is given as an
  !                      optional input,
  !             MXORDS = 5, unless a smaller value is given as an
  !                      optional input,
  !             LMAT   = length of matrix work space:
  !             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
  !             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
  !
  !                       --- Dynamic Length Case ---
  !          If the length of RWORK is to be dynamic, then it should
  !          be at least LRN or LRS, as defined above, depending on the
  !          current method.  Initially, it must be at least LRN (since
  !          DLSODA starts with the nonstiff method).  On any return
  !          from DLSODA, the optional output MCUR indicates the current
  !          method.  If MCUR differs from the value it had on the
  !          previous return, or if there has only been one call to
  !          DLSODA and MCUR is now 2, then DLSODA has switched
  !          methods during the last call, and the length of RWORK
  !          should be reset (to LRN if MCUR = 1, or to LRS if
  !          MCUR = 2).  (An increase in the RWORK length is required
  !          if DLSODA returned ISTATE = -7, but not otherwise.)
  !          After resetting the length, call DLSODA with ISTATE = 3
  !          to signal that change.
  !
  ! LRW    = the length of the array RWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! IWORK  = an integer array for work space.
  !          As DLSODA switches automatically between stiff and nonstiff
  !          methods, the required length of IWORK can change during
  !          problem, between
  !             LIS = 20 + NEQ   and   LIN = 20,
  !          respectively.  Thus the IWORK array passed to DLSODA can
  !          either have a fixed length of at least 20 + NEQ, or have a
  !          dynamic length of at least LIN or LIS, depending on the
  !          current method.  The comments on dynamic length under
  !          RWORK above apply here.  Initially, this length need
  !          only be at least LIN = 20.
  !
  !          The first few words of IWORK are used for conditional and
  !          optional inputs and optional outputs.
  !
  !          The following 2 words in IWORK are conditional inputs:
  !            IWORK(1) = ML     these are the lower and upper
  !            IWORK(2) = MU     half-bandwidths, respectively, of the
  !                       banded Jacobian, excluding the main diagonal.
  !                       The band is defined by the matrix locations
  !                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
  !                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
  !                       These are required if JT is 4 or 5, and
  !                       ignored otherwise.  ML and MU may in fact be
  !                       the band parameters for a matrix to which
  !                       df/dy is only approximately equal.
  !
  ! LIW    = the length of the array IWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! Note: The base addresses of the work arrays must not be
  ! altered between calls to DLSODA for the same problem.
  ! The contents of the work arrays must not be altered
  ! between calls, except possibly for the conditional and
  ! optional inputs, and except for the last 3*NEQ words of RWORK.
  ! The latter space is used for internal scratch space, and so is
  ! available for use by the user outside DLSODA between calls, if
  ! desired (but not for use by F or JAC).
  !
  ! JAC    = the name of the user-supplied routine to compute the
  !          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
  !          is optional, but if the problem is expected to be stiff much
  !          of the time, you are encouraged to supply JAC, for the sake
  !          of efficiency.  (Alternatively, set JT = 2 or 5 to have
  !          DLSODA compute df/dy internally by difference quotients.)
  !          If and when DLSODA uses df/dy, it treats this NEQ by NEQ
  !          matrix either as full (JT = 1 or 2), or as banded (JT =
  !          4 or 5) with half-bandwidths ML and MU (discussed under
  !          IWORK above).  In either case, if JT = 1 or 4, the JAC
  !          routine must compute df/dy as a function of the scalar t
  !          and the vector y.  It is to have the form
  !               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
  !               DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
  !          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
  !          PD is to be loaded with partial derivatives (elements of
  !          the Jacobian matrix) on output.  PD must be given a first
  !          dimension of NROWPD.  T and Y have the same meaning as in
  !          Subroutine F.
  !               In the full matrix case (JT = 1), ML and MU are
  !          ignored, and the Jacobian is to be loaded into PD in
  !          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
  !               In the band matrix case (JT = 4), the elements
  !          within the band are to be loaded into PD in columnwise
  !          manner, with diagonal lines of df/dy loaded into the rows
  !          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
  !          ML and MU are the half-bandwidth parameters (see IWORK).
  !          The locations in PD in the two triangular areas which
  !          correspond to nonexistent matrix elements can be ignored
  !          or loaded arbitrarily, as they are overwritten by DLSODA.
  !               JAC need not provide df/dy exactly.  A crude
  !          approximation (possibly with a smaller bandwidth) will do.
  !               In either case, PD is preset to zero by the solver,
  !          so that only the nonzero elements need be loaded by JAC.
  !          Each call to JAC is preceded by a call to F with the same
  !          arguments NEQ, T, and Y.  Thus to gain some efficiency,
  !          intermediate quantities shared by both calculations may be
  !          saved in a user Common block by F and not recomputed by JAC,
  !          if desired.  Also, JAC may alter the Y array, if desired.
  !          JAC must be declared External in the calling program.
  !               Subroutine JAC may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
  !          See the descriptions of NEQ and Y above.
  !
  ! JT     = Jacobian type indicator.  Used only for input.
  !          JT specifies how the Jacobian matrix df/dy will be
  !          treated, if and when DLSODA requires this matrix.
  !          JT has the following values and meanings:
  !           1 means a user-supplied full (NEQ by NEQ) Jacobian.
  !           2 means an internally generated (difference quotient) full
  !             Jacobian (using NEQ extra calls to F per df/dy value).
  !           4 means a user-supplied banded Jacobian.
  !           5 means an internally generated banded Jacobian (using
  !             ML+MU+1 extra calls to F per df/dy evaluation).
  !          If JT = 1 or 4, the user must supply a Subroutine JAC
  !          (the name is arbitrary) as described above under JAC.
  !          If JT = 2 or 5, a dummy argument can be used.
  !-----------------------------------------------------------------------
  ! Optional Inputs.
  !
  ! The following is a list of the optional inputs provided for in the
  ! call sequence.  (See also Part 2.)  For each such input variable,
  ! this table lists its name as used in this documentation, its
  ! location in the call sequence, its meaning, and the default value.
  ! The use of any of these inputs requires IOPT = 1, and in that
  ! case all of these inputs are examined.  A value of zero for any
  ! of these optional inputs will cause the default value to be used.
  ! Thus to use a subset of the optional inputs, simply preload
  ! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
  ! then set those of interest to nonzero values.
  !
  ! Name    Location      Meaning and Default Value
  !
  ! H0      RWORK(5)  the step size to be attempted on the first step.
  !                   The default value is determined by the solver.
  !
  ! HMAX    RWORK(6)  the maximum absolute step size allowed.
  !                   The default value is infinite.
  !
  ! HMIN    RWORK(7)  the minimum absolute step size allowed.
  !                   The default value is 0.  (This lower bound is not
  !                   enforced on the final step before reaching TCRIT
  !                   when ITASK = 4 or 5.)
  !
  ! IXPR    IWORK(5)  flag to generate extra printing at method switches.
  !                   IXPR = 0 means no extra printing (the default).
  !                   IXPR = 1 means print data on each switch.
  !                   T, H, and NST will be printed on the same logical
  !                   unit as used for error messages.
  !
  ! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
  !                   allowed during one call to the solver.
  !                   The default value is 500.
  !
  ! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
  !                   warning that T + H = T on a step (H = step size).
  !                   This must be positive to result in a non-default
  !                   value.  The default value is 10.
  !
  ! MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
  !                   (Adams) method.  the default value is 12.
  !                   if MXORDN exceeds the default value, it will
  !                   be reduced to the default value.
  !                   MXORDN is held constant during the problem.
  !
  ! MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
  !                   (BDF) method.  The default value is 5.
  !                   If MXORDS exceeds the default value, it will
  !                   be reduced to the default value.
  !                   MXORDS is held constant during the problem.
  !-----------------------------------------------------------------------
  ! Optional Outputs.
  !
  ! As optional additional output from DLSODA, the variables listed
  ! below are quantities related to the performance of DLSODA
  ! which are available to the user.  These are communicated by way of
  ! the work arrays, but also have internal mnemonic names as shown.
  ! except where stated otherwise, all of these outputs are defined
  ! on any successful return from DLSODA, and on any return with
  ! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
  ! (ISTATE = -3), they will be unchanged from their existing values
  ! (if any), except possibly for TOLSF, LENRW, and LENIW.
  ! On any error return, outputs relevant to the error will be defined,
  ! as noted below.
  !
  ! Name    Location      Meaning
  !
  ! HU      RWORK(11) the step size in t last used (successfully).
  !
  ! HCUR    RWORK(12) the step size to be attempted on the next step.
  !
  ! TCUR    RWORK(13) the current value of the independent variable
  !                   which the solver has actually reached, i.e. the
  !                   current internal mesh point in t.  On output, TCUR
  !                   will always be at least as far as the argument
  !                   T, but may be farther (if interpolation was done).
  !
  ! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
  !                   computed when a request for too much accuracy was
  !                   detected (ISTATE = -3 if detected at the start of
  !                   the problem, ISTATE = -2 otherwise).  If ITOL is
  !                   left unaltered but RTOL and ATOL are uniformly
  !                   scaled up by a factor of TOLSF for the next call,
  !                   then the solver is deemed likely to succeed.
  !                   (The user may also ignore TOLSF and alter the
  !                   tolerance parameters in any other way appropriate.)
  !
  ! TSW     RWORK(15) the value of t at the time of the last method
  !                   switch, if any.
  !
  ! NST     IWORK(11) the number of steps taken for the problem so far.
  !
  ! NFE     IWORK(12) the number of f evaluations for the problem so far.
  !
  ! NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
  !                   LU decompositions) for the problem so far.
  !
  ! NQU     IWORK(14) the method order last used (successfully).
  !
  ! NQCUR   IWORK(15) the order to be attempted on the next step.
  !
  ! IMXER   IWORK(16) the index of the component of largest magnitude in
  !                   the weighted local error vector ( E(i)/EWT(i) ),
  !                   on an error return with ISTATE = -4 or -5.
  !
  ! LENRW   IWORK(17) the length of RWORK actually required, assuming
  !                   that the length of RWORK is to be fixed for the
  !                   rest of the problem, and that switching may occur.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! LENIW   IWORK(18) the length of IWORK actually required, assuming
  !                   that the length of IWORK is to be fixed for the
  !                   rest of the problem, and that switching may occur.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! MUSED   IWORK(19) the method indicator for the last successful step:
  !                   1 means Adams (nonstiff), 2 means BDF (stiff).
  !
  ! MCUR    IWORK(20) the current method indicator:
  !                   1 means Adams (nonstiff), 2 means BDF (stiff).
  !                   This is the method to be attempted
  !                   on the next step.  Thus it differs from MUSED
  !                   only if a method switch has just been made.
  !
  ! The following two arrays are segments of the RWORK array which
  ! may also be of interest to the user as optional outputs.
  ! For each array, the table below gives its internal name,
  ! its base address in RWORK, and its description.
  !
  ! Name    Base Address      Description
  !
  ! YH      21             the Nordsieck history array, of size NYH by
  !                        (NQCUR + 1), where NYH is the initial value
  !                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
  !                        of YH contains HCUR**j/factorial(j) times
  !                        the j-th derivative of the interpolating
  !                        polynomial currently representing the solution,
  !                        evaluated at T = TCUR.
  !
  ! ACOR     LACOR         array of size NEQ used for the accumulated
  !         (from Common   corrections on each step, scaled on output
  !           as noted)    to represent the estimated local error in y
  !                        on the last step.  This is the vector E in
  !                        the description of the error control.  It is
  !                        defined only on a successful return from
  !                        DLSODA.  The base address LACOR is obtained by
  !                        including in the user's program the
  !                        following 2 lines:
  !                           COMMON /DLS001/ RLS(218), ILS(37)
  !                           LACOR = ILS(22)
  !
  !-----------------------------------------------------------------------
  ! Part 2.  Other Routines Callable.
  !
  ! The following are optional calls which the user may make to
  ! gain additional capabilities in conjunction with DLSODA.
  ! (The routines XSETUN and XSETF are designed to conform to the
  ! SLATEC error handling package.)
  !
  !     Form of Call                  Function
  !   CALL XSETUN(LUN)          set the logical unit number, LUN, for
  !                             output of messages from DLSODA, if
  !                             the default is not desired.
  !                             The default value of LUN is 6.
  !
  !   CALL XSETF(MFLAG)         set a flag to control the printing of
  !                             messages by DLSODA.
  !                             MFLAG = 0 means do not print. (Danger:
  !                             This risks losing valuable information.)
  !                             MFLAG = 1 means print (the default).
  !
  !                             Either of the above calls may be made at
  !                             any time and will take effect immediately.
  !
  !   CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of
  !                             the internal Common blocks used by
  !                             DLSODA (see Part 3 below).
  !                             RSAV must be a real array of length 240
  !                             or more, and ISAV must be an integer
  !                             array of length 46 or more.
  !                             JOB=1 means save Common into RSAV/ISAV.
  !                             JOB=2 means restore Common from RSAV/ISAV.
  !                                DSRCMA is useful if one is
  !                             interrupting a run and restarting
  !                             later, or alternating between two or
  !                             more problems solved with DLSODA.
  !
  !   CALL DINTDY(,,,,,)        provide derivatives of y, of various
  !        (see below)          orders, at a specified point t, if
  !                             desired.  It may be called only after
  !                             a successful return from DLSODA.
  !
  ! The detailed instructions for using DINTDY are as follows.
  ! The form of the call is:
  !
  !   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
  !
  ! The input parameters are:
  !
  ! T         = value of independent variable where answers are desired
  !             (normally the same as the T last returned by DLSODA).
  !             For valid results, T must lie between TCUR - HU and TCUR.
  !             (See optional outputs for TCUR and HU.)
  ! K         = integer order of the derivative desired.  K must satisfy
  !             0 .le. K .le. NQCUR, where NQCUR is the current order
  !             (see optional outputs).  The capability corresponding
  !             to K = 0, i.e. computing y(T), is already provided
  !             by DLSODA directly.  Since NQCUR .ge. 1, the first
  !             derivative dy/dt is always available with DINTDY.
  ! RWORK(21) = the base address of the history array YH.
  ! NYH       = column length of YH, equal to the initial value of NEQ.
  !
  ! The output parameters are:
  !
  ! DKY       = a real array of length NEQ containing the computed value
  !             of the K-th derivative of y(t).
  ! IFLAG     = integer flag, returned as 0 if K and T were legal,
  !             -1 if K was illegal, and -2 if T was illegal.
  !             On an error return, a message is also written.
  !-----------------------------------------------------------------------
  ! Part 3.  Common Blocks.
  !
  ! If DLSODA is to be used in an overlay situation, the user
  ! must declare, in the primary overlay, the variables in:
  !   (1) the call sequence to DLSODA, and
  !   (2) the two internal Common blocks
  !         /DLS001/  of length  255  (218 double precision words
  !                      followed by 37 integer words),
  !         /DLSA01/  of length  31    (22 double precision words
  !                      followed by  9 integer words).
  !
  ! If DLSODA is used on a system in which the contents of internal
  ! Common blocks are not preserved between calls, the user should
  ! declare the above Common blocks in the calling program to insure
  ! that their contents are preserved.
  !
  ! If the solution of a given problem by DLSODA is to be interrupted
  ! and then later continued, such as when restarting an interrupted run
  ! or alternating between two or more problems, the user should save,
  ! following the return from the last DLSODA call prior to the
  ! interruption, the contents of the call sequence variables and the
  ! internal Common blocks, and later restore these values before the
  ! next DLSODA call for that problem.  To save and restore the Common
  ! blocks, use Subroutine DSRCMA (see Part 2 above).
  !
  !-----------------------------------------------------------------------
  ! Part 4.  Optionally Replaceable Solver Routines.
  !
  ! Below is a description of a routine in the DLSODA package which
  ! relates to the measurement of errors, and can be
  ! replaced by a user-supplied version, if desired.  However, since such
  ! a replacement may have a major impact on performance, it should be
  ! done only when absolutely necessary, and only with great caution.
  ! (Note: The means by which the package version of a routine is
  ! superseded by the user's version may be system-dependent.)
  !
  ! (a) DEWSET.
  ! The following subroutine is called just before each internal
  ! integration step, and sets the array of error weights, EWT, as
  ! described under ITOL/RTOL/ATOL above:
  !     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
  ! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence,
  ! YCUR contains the current dependent variable vector, and
  ! EWT is the array of weights set by DEWSET.
  !
  ! If the user supplies this subroutine, it must return in EWT(i)
  ! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
  ! in y(i) to.  The EWT array returned by DEWSET is passed to the
  ! DMNORM routine, and also used by DLSODA in the computation
  ! of the optional output IMXER, and the increments for difference
  ! quotient Jacobians.
  !
  ! In the user-supplied version of DEWSET, it may be desirable to use
  ! the current values of derivatives of y.  Derivatives up to order NQ
  ! are available from the history array YH, described above under
  ! optional outputs.  In DEWSET, YH is identical to the YCUR array,
  ! extended to NQ + 1 columns with a column length of NYH and scale
  ! factors of H**j/factorial(j).  On the first call for the problem,
  ! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
  ! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
  ! can be obtained by including in DEWSET the statements:
  !     DOUBLE PRECISION RLS
  !     COMMON /DLS001/ RLS(218),ILS(37)
  !     NQ = ILS(33)
  !     NST = ILS(34)
  !     H = RLS(212)
  ! Thus, for example, the current value of dy/dt can be obtained as
  ! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
  ! unnecessary when NST = 0).
  !-----------------------------------------------------------------------
  !
  !***REVISION HISTORY  (YYYYMMDD)
  ! 19811102  DATE WRITTEN
  ! 19820126  Fixed bug in tests of work space lengths;
  !           minor corrections in main prologue and comments.
  ! 19870330  Major update: corrected comments throughout;
  !           removed TRET from Common; rewrote EWSET with 4 loops;
  !           fixed t test in INTDY; added Cray directives in STODA;
  !           in STODA, fixed DELP init. and logic around PJAC call;
  !           combined routines to save/restore Common;
  !           passed LEVEL = 0 in error message calls (except run abort).
  ! 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA.
  ! 20010425  Major update: convert source lines to upper case;
  !           added *DECK lines; changed from 1 to * in dummy dimensions;
  !           changed names R1MACH/D1MACH to RUMACH/DUMACH;
  !           renamed routines for uniqueness across single/double prec.;
  !           converted intrinsic names to generic form;
  !           removed ILLIN and NTREP (data loaded) from Common;
  !           removed all 'own' variables from Common;
  !           changed error messages to quoted strings;
  !           replaced XERRWV/XERRWD with 1993 revised version;
  !           converted prologues, comments, error messages to mixed case;
  !           numerous corrections to prologues and internal comments.
  ! 20010507  Converted single precision source to double precision.
  ! 20010613  Revised excess accuracy test (to match rest of ODEPACK).
  ! 20010808  Fixed bug in DPRJA (matrix in DBNORM call).
  ! 20020502  Corrected declarations in descriptions of user routines.
  ! 20031105  Restored 'own' variables to Common blocks, to enable
  !           interrupt/restart feature.
  ! 20031112  Added SAVE statements for data-loaded constants.
  !
  !-----------------------------------------------------------------------
  ! Other routines in the DLSODA package.
  !
  ! In addition to Subroutine DLSODA, the DLSODA package includes the
  ! following subroutines and function routines:
  !  DINTDY   computes an interpolated value of the y vector at t = TOUT.
  !  DSTODA   is the core integrator, which does one step of the
  !           integration and the associated error control.
  !  DCFODE   sets all method coefficients and test constants.
  !  DPRJA    computes and preprocesses the Jacobian matrix J = df/dy
  !           and the Newton iteration matrix P = I - h*l0*J.
  !  DSOLSY   manages solution of linear system in chord iteration.
  !  DEWSET   sets the error weight vector EWT before each step.
  !  DMNORM   computes the weighted max-norm of a vector.
  !  DFNORM   computes the norm of a full matrix consistent with the
  !           weighted max-norm on vectors.
  !  DBNORM   computes the norm of a band matrix consistent with the
  !           weighted max-norm on vectors.
  !  DSRCMA   is a user-callable routine to save and restore
  !           the contents of the internal Common blocks.
  !  DGEFA and DGESL   are routines from LINPACK for solving full
  !           systems of linear algebraic equations.
  !  DGBFA and DGBSL   are routines from LINPACK for solving banded
  !           linear systems.
  !  DUMACH   computes the unit roundoff in a machine-independent manner.
  !  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
  !           error messages and warnings.  XERRWD is machine-dependent.
  ! Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are
  ! function routines.  All the others are subroutines.
  !
  !-----------------------------------------------------------------------
        !EXTERNAL DPRJA , DSOLSY
        !DOUBLE PRECISION DUMACH , DMNORM
        INTEGER INIt , MXStep , MXHnil , NHNil , NSLast , NYH , IOWns ,   &
       &        ICF , IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH ,   &
       &        LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,        &
       &        MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST , NFE ,     &
       &        NJE , NQU
        INTEGER INSufr , INSufi , IXPr , IOWns2 , JTYp , MUSed , MXOrdn , &
       &        MXOrds
        INTEGER i , i1 , i2 , iflag , imxer , kgo , lf0 , leniw , lenrw , &
       &        lenwm , ml , mord , mu , mxhnl0 , mxstp0
        INTEGER len1 , len1c , len1n , len1s , len2 , leniwc , lenrwc
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION TSW , ROWns2 , PDNorm
        DOUBLE PRECISION atoli , ayi , big , ewti , h0 , hmax , hmx , rh ,&
       &                 rtoli , tcrit , tdist , tnext , tol , tolsf ,    &
       &                 tp , size , sum , w0
        DIMENSION mord(2)
        LOGICAL ihit
        CHARACTER*60 msg
        SAVE mord , mxstp0 , mxhnl0
  !-----------------------------------------------------------------------
  ! The following two internal Common blocks contain
  ! (a) variables which are local to any subroutine but whose values must
  !     be preserved between calls to the routine ("own" variables), and
  ! (b) variables which are communicated between subroutines.
  ! The block DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA,
  ! DPRJA, and DSOLSY.
  ! The block DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA.
  ! Groups of variables are replaced by dummy arrays in the Common
  ! declarations in routines where those variables are not used.
  !-----------------------------------------------------------------------
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , INIt , MXStep , MXHnil ,       &
       &                NHNil , NSLast , NYH , IOWns(6) , ICF , IERpj ,   &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
        COMMON /DLSA01/ TSW , ROWns2(20) , PDNorm , INSufr , INSufi ,     &
       &                IXPr , IOWns2(2) , JTYp , MUSed , MXOrdn , MXOrds
  !
        DATA mord(1) , mord(2)/12 , 5/ , mxstp0/500/ , mxhnl0/10/
  !-----------------------------------------------------------------------
  ! Block A.
  ! This code block is executed on every call.
  ! It tests ISTATE and ITASK for legality and branches appropriately.
  ! If ISTATE .gt. 1 but the flag INIT shows that initialization has
  ! not yet been done, an error return occurs.
  ! If ISTATE = 1 and TOUT = T, return immediately.
  !-----------------------------------------------------------------------
        IF ( Istate.LT.1 .OR. Istate.GT.3 ) THEN
  !-----------------------------------------------------------------------
  ! Block I.
  ! The following block handles all error returns due to illegal input
  ! (ISTATE = -3), as detected before calling the core integrator.
  ! First the error message routine is called.  If the illegal input
  ! is a negative ISTATE, the run is aborted (apparent infinite loop).
  !-----------------------------------------------------------------------
           msg = 'DLSODA-  ISTATE (=I1) illegal.'
           CALL XERRWD(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
           IF ( Istate.GE.0 ) GOTO 2800
  !
           msg = 'DLSODA-  Run aborted.. apparent infinite loop.    '
           CALL XERRWD(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
           GOTO 99999
        ELSE
           IF ( Itask.LT.1 .OR. Itask.GT.5 ) THEN
              msg = 'DLSODA-  ITASK (=I1) illegal. '
              CALL XERRWD(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
              GOTO 2800
           ELSE
              IF ( Istate.EQ.1 ) THEN
                 INIt = 0
                 IF ( Tout.EQ.T ) RETURN
              ELSEIF ( INIt.EQ.0 ) THEN
                 msg =                                                    &
       &              'DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.'
                 CALL XERRWD(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
                 GOTO 2800
              ELSEIF ( Istate.EQ.2 ) THEN
                 GOTO 50
              ENDIF
  !-----------------------------------------------------------------------
  ! Block B.
  ! The next code block is executed for the initial call (ISTATE = 1),
  ! or for a continuation call with parameter changes (ISTATE = 3).
  ! It contains checking of all inputs and various initializations.
  !
  ! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
  ! JT, ML, and MU.
  !-----------------------------------------------------------------------
              IF ( Neq(1).LE.0 ) THEN
                 msg = 'DLSODA-  NEQ (=I1) .lt. 1     '
                 CALL XERRWD(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
                 GOTO 2800
              ELSE
                 IF ( Istate.NE.1 ) THEN
                    IF ( Neq(1).GT.N ) THEN
                       msg =                                              &
       &              'DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '
                       CALL XERRWD(msg,50,5,0,2,N,Neq(1),0,0.0D0,0.0D0)
                       GOTO 2800
                    ENDIF
                 ENDIF
                 N = Neq(1)
                 IF ( Itol.LT.1 .OR. Itol.GT.4 ) THEN
                    msg = 'DLSODA-  ITOL (=I1) illegal.  '
                    CALL XERRWD(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
                    GOTO 2800
                 ELSEIF ( Iopt.LT.0 .OR. Iopt.GT.1 ) THEN
                    msg = 'DLSODA-  IOPT (=I1) illegal.  '
                    CALL XERRWD(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
                    GOTO 2800
                 ELSEIF ( Jt.EQ.3 .OR. Jt.LT.1 .OR. Jt.GT.5 ) THEN
                    msg = 'DLSODA-  JT (=I1) illegal.    '
                    CALL XERRWD(msg,30,8,0,1,Jt,0,0,0.0D0,0.0D0)
                    GOTO 2800
                 ELSE
                    JTYp = Jt
                    IF ( Jt.GT.2 ) THEN
                       ml = Iwork(1)
                       mu = Iwork(2)
                       IF ( ml.LT.0 .OR. ml.GE.N ) THEN
                          msg =                                           &
       &              'DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
                          CALL XERRWD(msg,50,9,0,2,ml,Neq(1),0,0.0D0,     &
       &                              0.0D0)
                          GOTO 2800
                       ELSEIF ( mu.LT.0 .OR. mu.GE.N ) THEN
                          msg =                                           &
       &              'DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
                          CALL XERRWD(msg,50,10,0,2,mu,Neq(1),0,0.0D0,    &
       &                              0.0D0)
                          GOTO 2800
                       ENDIF
                    ENDIF
  ! Next process and check the optional inputs. --------------------------
                    IF ( Iopt.EQ.1 ) THEN
                       IXPr = Iwork(5)
                       IF ( IXPr.LT.0 .OR. IXPr.GT.1 ) THEN
                          msg = 'DLSODA-  IXPR (=I1) illegal.  '
                          CALL XERRWD(msg,30,11,0,1,IXPr,0,0,0.0D0,0.0D0)
                          GOTO 2800
                       ELSE
                          MXStep = Iwork(6)
                          IF ( MXStep.LT.0 ) THEN
                             msg = 'DLSODA-  MXSTEP (=I1) .lt. 0  '
                             CALL XERRWD(msg,30,12,0,1,MXStep,0,0,0.0D0,  &
       &                                 0.0D0)
                             GOTO 2800
                          ELSE
                             IF ( MXStep.EQ.0 ) MXStep = mxstp0
                             MXHnil = Iwork(7)
                             IF ( MXHnil.LT.0 ) THEN
                                msg = 'DLSODA-  MXHNIL (=I1) .lt. 0  '
                                CALL XERRWD(msg,30,13,0,1,MXHnil,0,0,     &
       &                           0.0D0,0.0D0)
                                GOTO 2800
                             ELSE
                                IF ( MXHnil.EQ.0 ) MXHnil = mxhnl0
                                IF ( Istate.EQ.1 ) THEN
                                   h0 = Rwork(5)
                                   MXOrdn = Iwork(8)
                                   IF ( MXOrdn.LT.0 ) THEN
                                      msg =                               &
       &                                 'DLSODA-  MXORDN (=I1) .lt. 0  '
                                      CALL XERRWD(msg,30,28,0,1,MXOrdn,0, &
       &                                 0,0.0D0,0.0D0)
                                      GOTO 2800
                                   ELSE
                                      IF ( MXOrdn.EQ.0 ) MXOrdn = 100
                                      MXOrdn = MIN(MXOrdn,mord(1))
                                      MXOrds = Iwork(9)
                                      IF ( MXOrds.LT.0 ) THEN
                                         msg =                            &
       &                                  'DLSODA-  MXORDS (=I1) .lt. 0  '
                                         CALL XERRWD(msg,30,29,0,1,MXOrds,&
       &                                    0,0,0.0D0,0.0D0)
                                         GOTO 2800
                                      ELSE
                                         IF ( MXOrds.EQ.0 ) MXOrds = 100
                                         MXOrds = MIN(MXOrds,mord(2))
                                         IF ( (Tout-T)*h0.LT.0.0D0 ) THEN
                                           msg =                          &
       &                        'DLSODA-  TOUT (=R1) behind T (=R2)      '
                                           CALL XERRWD(msg,40,14,0,0,0,0, &
       &                                      2,Tout,T)
                                           msg =                          &
       &              '      Integration direction is given by H0 (=R1)  '
                                           CALL XERRWD(msg,50,14,0,0,0,0, &
       &                                      1,h0,0.0D0)
                                           GOTO 2800
                                         ENDIF
                                      ENDIF
                                   ENDIF
                                ENDIF
                                hmax = Rwork(6)
                                IF ( hmax.LT.0.0D0 ) THEN
                                   msg = 'DLSODA-  HMAX (=R1) .lt. 0.0  '
                                   CALL XERRWD(msg,30,15,0,0,0,0,1,hmax,  &
       &                              0.0D0)
                                   GOTO 2800
                                ELSE
                                   HMXi = 0.0D0
                                   IF ( hmax.GT.0.0D0 ) HMXi = 1.0D0/hmax
                                   HMIn = Rwork(7)
                                   IF ( HMIn.LT.0.0D0 ) THEN
                                      msg =                               &
       &                                 'DLSODA-  HMIN (=R1) .lt. 0.0  '
                                      CALL XERRWD(msg,30,16,0,0,0,0,1,    &
       &                                 HMIn,0.0D0)
                                      GOTO 2800
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ELSE
                       IXPr = 0
                       MXStep = mxstp0
                       MXHnil = mxhnl0
                       HMXi = 0.0D0
                       HMIn = 0.0D0
                       IF ( Istate.EQ.1 ) THEN
                          h0 = 0.0D0
                          MXOrdn = mord(1)
                          MXOrds = mord(2)
                       ENDIF
                    ENDIF
  !-----------------------------------------------------------------------
  ! Set work array pointers and check lengths LRW and LIW.
  ! If ISTATE = 1, METH is initialized to 1 here to facilitate the
  ! checking of work space lengths.
  ! Pointers to segments of RWORK and IWORK are named by prefixing L to
  ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
  ! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
  ! If the lengths provided are insufficient for the current method,
  ! an error return occurs.  This is treated as illegal input on the
  ! first call, but as a problem interruption with ISTATE = -7 on a
  ! continuation call.  If the lengths are sufficient for the current
  ! method but not for both methods, a warning message is sent.
  !-----------------------------------------------------------------------
                    IF ( Istate.EQ.1 ) METh = 1
                    IF ( Istate.EQ.1 ) NYH = N
                    LYH = 21
                    len1n = 20 + (MXOrdn+1)*NYH
                    len1s = 20 + (MXOrds+1)*NYH
                    LWM = len1s + 1
                    IF ( Jt.LE.2 ) lenwm = N*N + 2
                    IF ( Jt.GE.4 ) lenwm = (2*ml+mu+1)*N + 2
                    len1s = len1s + lenwm
                    len1c = len1n
                    IF ( METh.EQ.2 ) len1c = len1s
                    len1 = MAX(len1n,len1s)
                    len2 = 3*N
                    lenrw = len1 + len2
                    lenrwc = len1c + len2
                    Iwork(17) = lenrw
                    LIWm = 1
                    leniw = 20 + N
                    leniwc = 20
                    IF ( METh.EQ.2 ) leniwc = leniw
                    Iwork(18) = leniw
                    IF ( Istate.EQ.1 .AND. Lrw.LT.lenrwc ) THEN
                       msg =                                              &
       &    'DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                       CALL XERRWD(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                       GOTO 2800
                    ELSEIF ( Istate.EQ.1 .AND. Liw.LT.leniwc ) THEN
                       msg =                                              &
       &    'DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                       CALL XERRWD(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                       GOTO 2800
                    ELSE
                       IF ( Istate.EQ.3 .AND. Lrw.LT.lenrwc ) GOTO 1800
                       IF ( Istate.EQ.3 .AND. Liw.LT.leniwc ) GOTO 1900
                       LEWt = len1 + 1
                       INSufr = 0
                       IF ( Lrw.LT.lenrw ) THEN
                          INSufr = 2
                          LEWt = len1c + 1
                          msg =                                           &
       &    'DLSODA-  Warning.. RWORK length is sufficient for now, but  '
                          CALL XERRWD(msg,60,103,0,0,0,0,0,0.0D0,0.0D0)
                          msg =                                           &
       &    '      may not be later.  Integration will proceed anyway.   '
                          CALL XERRWD(msg,60,103,0,0,0,0,0,0.0D0,0.0D0)
                          msg =                                           &
       &              '      Length needed is LENRW = I1, while LRW = I2.'
                          CALL XERRWD(msg,50,103,0,2,lenrw,Lrw,0,0.0D0,   &
       &                              0.0D0)
                       ENDIF
                       LSAvf = LEWt + N
                       LACor = LSAvf + N
                       INSufi = 0
                       IF ( Liw.LT.leniw ) THEN
                          INSufi = 2
                          msg =                                           &
       &    'DLSODA-  Warning.. IWORK length is sufficient for now, but  '
                          CALL XERRWD(msg,60,104,0,0,0,0,0,0.0D0,0.0D0)
                          msg =                                           &
       &    '      may not be later.  Integration will proceed anyway.   '
                          CALL XERRWD(msg,60,104,0,0,0,0,0,0.0D0,0.0D0)
                          msg =                                           &
       &              '      Length needed is LENIW = I1, while LIW = I2.'
                          CALL XERRWD(msg,50,104,0,2,leniw,Liw,0,0.0D0,   &
       &                              0.0D0)
                       ENDIF
  ! Check RTOL and ATOL for legality. ------------------------------------
                       rtoli = Rtol(1)
                       atoli = Atol(1)
                       DO i = 1 , N
                          IF ( Itol.GE.3 ) rtoli = Rtol(i)
                          IF ( Itol.EQ.2 .OR. Itol.EQ.4 ) atoli = Atol(i)
                          IF ( rtoli.LT.0.0D0 ) GOTO 2200
                          IF ( atoli.LT.0.0D0 ) GOTO 2300
                       ENDDO
                       IF ( Istate.EQ.1 ) THEN
  !-----------------------------------------------------------------------
  ! Block C.
  ! The next block is for the initial call only (ISTATE = 1).
  ! It contains all remaining initializations, the initial call to F,
  ! and the calculation of the initial step size.
  ! The error weights in EWT are inverted after being loaded.
  !-----------------------------------------------------------------------
                          UROund = DUMACH()
                          TN = T
                          TSW = T
                          MAXord = MXOrdn
                          IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
                             tcrit = Rwork(1)
                             IF ( (tcrit-Tout)*(Tout-T).LT.0.0D0 )        &
       &                          GOTO 2600
                             IF ( h0.NE.0.0D0 .AND. (T+h0-tcrit)          &
       &                          *h0.GT.0.0D0 ) h0 = tcrit - T
                          ENDIF
                          JSTart = 0
                          NHNil = 0
                          NST = 0
                          NJE = 0
                          NSLast = 0
                          HU = 0.0D0
                          NQU = 0
                          MUSed = 0
                          MITer = 0
                          CCMax = 0.3D0
                          MAXcor = 3
                          MSBp = 20
                          MXNcf = 10
  ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
                          lf0 = LYH + NYH
                          CALL F(Neq,T,Y,Rwork(lf0))
                          NFE = 1
  ! Load the initial value vector in YH. ---------------------------------
                          DO i = 1 , N
                             Rwork(i+LYH-1) = Y(i)
                          ENDDO
  ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                          NQ = 1
                          H = 1.0D0
                          CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),        &
       &                              Rwork(LEWt))
                          DO i = 1 , N
                             IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 2400
                             Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
                          ENDDO
  !-----------------------------------------------------------------------
  ! The coding below computes the step size, H0, to be attempted on the
  ! first step, unless the user has supplied a value for this.
  ! First check that TOUT - T differs significantly from zero.
  ! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
  ! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
  ! so as to be between 100*UROUND and 1.0E-3.
  ! Then the computed value H0 is given by:
  !
  !   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
  !
  ! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
  !         F      = the initial value of the vector f(t,y), and
  !         norm() = the weighted vector norm used throughout, given by
  !                  the DMNORM function routine, and weighted by the
  !                  tolerances initially loaded into the EWT array.
  ! The sign of H0 is inferred from the initial values of TOUT and T.
  ! ABS(H0) is made .le. ABS(TOUT-T) in any case.
  !-----------------------------------------------------------------------
                          IF ( h0.EQ.0.0D0 ) THEN
                             tdist = ABS(Tout-T)
                             w0 = MAX(ABS(T),ABS(Tout))
                             IF ( tdist.LT.2.0D0*UROund*w0 ) THEN
                                msg =                                     &
       &    'DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.'
                                CALL XERRWD(msg,60,22,0,0,0,0,2,Tout,T)
                                GOTO 2800
                             ELSE
                                tol = Rtol(1)
                                IF ( Itol.GT.2 ) THEN
                                   DO i = 1 , N
                                      tol = MAX(tol,Rtol(i))
                                   ENDDO
                                ENDIF
                                IF ( tol.LE.0.0D0 ) THEN
                                   atoli = Atol(1)
                                   DO i = 1 , N
                                      IF ( Itol.EQ.2 .OR. Itol.EQ.4 )     &
       &                                 atoli = Atol(i)
                                      ayi = ABS(Y(i))
                                      IF ( ayi.NE.0.0D0 )                 &
       &                                 tol = MAX(tol,atoli/ayi)
                                   ENDDO
                                ENDIF
                                tol = MAX(tol,100.0D0*UROund)
                                tol = MIN(tol,0.001D0)
                                sum = DMNORM(N,Rwork(lf0),Rwork(LEWt))
                                sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                                h0 = 1.0D0/SQRT(sum)
                                h0 = MIN(h0,tdist)
                                h0 = SIGN(h0,Tout-T)
                             ENDIF
                          ENDIF
  ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                          rh = ABS(h0)*HMXi
                          IF ( rh.GT.1.0D0 ) h0 = h0/rh
  ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                          H = h0
                          DO i = 1 , N
                             Rwork(i+lf0-1) = h0*Rwork(i+lf0-1)
                          ENDDO
                          GOTO 700
                       ELSE
  ! If ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
                          JSTart = -1
                          IF ( N.NE.NYH ) THEN
  ! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                             i1 = LYH + L*NYH
                             i2 = LYH + (MAXord+1)*NYH - 1
                             IF ( i1.LE.i2 ) THEN
                                DO i = i1 , i2
                                   Rwork(i) = 0.0D0
                                ENDDO
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  ! Block D.
  ! The next code block is for continuation calls only (ISTATE = 2 or 3)
  ! and is to check stop conditions before taking a step.
  !-----------------------------------------------------------------------
   50      NSLast = NST
           GOTO (100,600,200,300,400) , Itask
        ENDIF
   100  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2700
        T = Tout
        GOTO 1400
   200  tp = TN - HU*(1.0D0+100.0D0*UROund)
        IF ( (tp-Tout)*H.GT.0.0D0 ) THEN
           msg =                                                          &
       &    'DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
           CALL XERRWD(msg,60,23,0,1,Itask,0,2,Tout,tp)
           GOTO 2800
        ELSE
           IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
           T = TN
           GOTO 1300
        ENDIF
   300  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2500
        IF ( (tcrit-Tout)*H.LT.0.0D0 ) GOTO 2600
        IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 500
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2700
        T = Tout
        GOTO 1400
   400  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2500
   500  hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
        IF ( ihit ) T = tcrit
        IF ( ihit ) GOTO 1300
        tnext = TN + H*(1.0D0+4.0D0*UROund)
        IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
           H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
           IF ( Istate.EQ.2 .AND. JSTart.GE.0 ) JSTart = -2
        ENDIF
  !-----------------------------------------------------------------------
  ! Block E.
  ! The next block is normally executed for all calls and contains
  ! the call to the one-step core integrator DSTODA.
  !
  ! This is a looping point for the integration steps.
  !
  ! First check for too many steps being taken, update EWT (if not at
  ! start of problem), check for too much accuracy being requested, and
  ! check for H below the roundoff level in T.
  !-----------------------------------------------------------------------
   600  IF ( METh.NE.MUSed ) THEN
           IF ( INSufr.EQ.1 ) GOTO 1800
           IF ( INSufi.EQ.1 ) GOTO 1900
        ENDIF
        IF ( (NST-NSLast).GE.MXStep ) THEN
  !-----------------------------------------------------------------------
  ! Block H.
  ! The following block handles all unsuccessful returns other than
  ! those for illegal input.  First the error message routine is called.
  ! If there was an error test or convergence test failure, IMXER is set.
  ! Then Y is loaded from YH and T is set to TN.
  ! The optional outputs are loaded into the work arrays before returning.
  !-----------------------------------------------------------------------
  ! The maximum number of steps was taken before reaching TOUT. ----------
           msg = 'DLSODA-  At current T (=R1), MXSTEP (=I1) steps   '
           CALL XERRWD(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      taken on this call before reaching TOUT     '
           CALL XERRWD(msg,50,201,0,1,MXStep,0,1,TN,0.0D0)
           Istate = -1
           GOTO 2100
        ELSE
           CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
           DO i = 1 , N
              IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 1500
              Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
           ENDDO
        ENDIF
   700  tolsf = UROund*DMNORM(N,Rwork(LYH),Rwork(LEWt))
        IF ( tolsf.LE.1.0D0 ) THEN
           IF ( (TN+H).EQ.TN ) THEN
              NHNil = NHNil + 1
              IF ( NHNil.LE.MXHnil ) THEN
                 msg =                                                    &
       &              'DLSODA-  Warning..Internal T (=R1) and H (=R2) are'
                 CALL XERRWD(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &    '      such that in the machine, T + H = T on the next step  '
                 CALL XERRWD(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &              '     (H = step size). Solver will continue anyway.'
                 CALL XERRWD(msg,50,101,0,0,0,0,2,TN,H)
                 IF ( NHNil.GE.MXHnil ) THEN
                    msg =                                                 &
       &              'DLSODA-  Above warning has been issued I1 times.  '
                    CALL XERRWD(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
                    msg =                                                 &
       &              '     It will not be issued again for this problem.'
                    CALL XERRWD(msg,50,102,0,1,MXHnil,0,0,0.0D0,0.0D0)
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  !   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
  !-----------------------------------------------------------------------
           CALL DSTODA(Neq,Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),       &
       &               Rwork(LSAvf),Rwork(LACor),Rwork(LWM),Iwork(LIWm),F,&
       &               JAC,DPRJA,DSOLSY)
           kgo = 1 - KFLag
           GOTO (800,1600,1700) , kgo
        ELSE
           tolsf = tolsf*2.0D0
           IF ( NST.EQ.0 ) THEN
              msg = 'DLSODA-  At start of problem, too much accuracy   '
              CALL XERRWD(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
              msg =                                                       &
       &    '      requested for precision of machine..  See TOLSF (=R1) '
              CALL XERRWD(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
              Rwork(14) = tolsf
              GOTO 2800
           ELSE
  ! Too much accuracy requested for machine precision. -------------------
              msg = 'DLSODA-  At T (=R1), too much accuracy requested  '
              CALL XERRWD(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
              msg = '      for precision of machine..  See TOLSF (=R2) '
              CALL XERRWD(msg,50,203,0,0,0,0,2,TN,tolsf)
              Rwork(14) = tolsf
              Istate = -2
              GOTO 2100
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! Block F.
  ! The following block handles the case of a successful return from the
  ! core integrator (KFLAG = 0).
  ! If a method switch was just made, record TSW, reset MAXORD,
  ! set JSTART to -1 to signal DSTODA to complete the switch,
  ! and do extra printing of data if IXPR = 1.
  ! Then, in any case, check for stop conditions.
  !-----------------------------------------------------------------------
   800  INIt = 1
        IF ( METh.NE.MUSed ) THEN
           TSW = TN
           MAXord = MXOrdn
           IF ( METh.EQ.2 ) MAXord = MXOrds
           IF ( METh.EQ.2 ) Rwork(LWM) = SQRT(UROund)
           INSufr = MIN(INSufr,1)
           INSufi = MIN(INSufi,1)
           JSTart = -1
           IF ( IXPr.NE.0 ) THEN
              IF ( METh.EQ.2 ) THEN
                 msg =                                                    &
       &    'DLSODA- A switch to the BDF (stiff) method has occurred     '
                 CALL XERRWD(msg,60,105,0,0,0,0,0,0.0D0,0.0D0)
              ENDIF
              IF ( METh.EQ.1 ) THEN
                 msg =                                                    &
       &    'DLSODA- A switch to the Adams (nonstiff) method has occurred'
                 CALL XERRWD(msg,60,106,0,0,0,0,0,0.0D0,0.0D0)
              ENDIF
              msg =                                                       &
       &    '     at T = R1,  tentative step size H = R2,  step NST = I1 '
              CALL XERRWD(msg,60,107,0,1,NST,0,2,TN,H)
           ENDIF
        ENDIF
        GOTO (900,1300,1000,1100,1200) , Itask
  ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   900  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        T = Tout
        GOTO 1400
  ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   1000 IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        GOTO 1300
  ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
   1100 IF ( (TN-Tout)*H.LT.0.0D0 ) THEN
           hmx = ABS(TN) + ABS(H)
           ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
           IF ( ihit ) GOTO 1300
           tnext = TN + H*(1.0D0+4.0D0*UROund)
           IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
              H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
              IF ( JSTart.GE.0 ) JSTart = -2
           ENDIF
           GOTO 600
        ELSE
           CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
           T = Tout
           GOTO 1400
        ENDIF
  ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
   1200 hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
  !-----------------------------------------------------------------------
  ! Block G.
  ! The following block handles all successful returns from DLSODA.
  ! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
  ! ISTATE is set to 2, and the optional outputs are loaded into the
  ! work arrays before returning.
  !-----------------------------------------------------------------------
   1300 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
           IF ( ihit ) T = tcrit
        ENDIF
   1400 Istate = 2
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Rwork(15) = TSW
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        Iwork(19) = MUSed
        Iwork(20) = METh
        RETURN
  ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
   1500 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'
        CALL XERRWD(msg,50,202,0,1,i,0,2,TN,ewti)
        Istate = -6
        GOTO 2100
  ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   1600 msg = 'DLSODA-  At T(=R1) and step size H(=R2), the error'
        CALL XERRWD(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      test failed repeatedly or with ABS(H) = HMIN'
        CALL XERRWD(msg,50,204,0,0,0,0,2,TN,H)
        Istate = -4
        GOTO 2000
  ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
   1700 msg = 'DLSODA-  At T (=R1) and step size H (=R2), the    '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      corrector convergence failed repeatedly     '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      or with ABS(H) = HMIN   '
        CALL XERRWD(msg,30,205,0,0,0,0,2,TN,H)
        Istate = -5
        GOTO 2000
  ! RWORK length too small to proceed. -----------------------------------
   1800 msg = 'DLSODA-  At current T(=R1), RWORK length too small'
        CALL XERRWD(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
        msg =                                                             &
       &    '      to proceed.  The integration was otherwise successful.'
        CALL XERRWD(msg,60,206,0,0,0,0,1,TN,0.0D0)
        Istate = -7
        GOTO 2100
  ! IWORK length too small to proceed. -----------------------------------
   1900 msg = 'DLSODA-  At current T(=R1), IWORK length too small'
        CALL XERRWD(msg,50,207,0,0,0,0,0,0.0D0,0.0D0)
        msg =                                                             &
       &    '      to proceed.  The integration was otherwise successful.'
        CALL XERRWD(msg,60,207,0,0,0,0,1,TN,0.0D0)
        Istate = -7
        GOTO 2100
  ! Compute IMXER if relevant. -------------------------------------------
   2000 big = 0.0D0
        imxer = 1
        DO i = 1 , N
           size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
           IF ( big.LT.size ) THEN
              big = size
              imxer = i
           ENDIF
        ENDDO
        Iwork(16) = imxer
  ! Set Y vector, T, and optional outputs. -------------------------------
   2100 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Rwork(15) = TSW
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        Iwork(19) = MUSed
        Iwork(20) = METh
        RETURN
   2200 msg = 'DLSODA-  RTOL(I1) is R1 .lt. 0.0        '
        CALL XERRWD(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
        GOTO 2800
   2300 msg = 'DLSODA-  ATOL(I1) is R1 .lt. 0.0        '
        CALL XERRWD(msg,40,20,0,1,i,0,1,atoli,0.0D0)
        GOTO 2800
   2400 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODA-  EWT(I1) is R1 .le. 0.0         '
        CALL XERRWD(msg,40,21,0,1,i,0,1,ewti,0.0D0)
        GOTO 2800
   2500 msg =                                                             &
       &    'DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
        CALL XERRWD(msg,60,24,0,0,0,0,2,tcrit,TN)
        GOTO 2800
   2600 msg =                                                             &
       &    'DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
        CALL XERRWD(msg,60,25,0,0,0,0,2,tcrit,Tout)
        GOTO 2800
   2700 msg = 'DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
        CALL XERRWD(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
  !
   2800 Istate = -3
        RETURN
  !----------------------- End of Subroutine DLSODA ----------------------
  99999 END
  !*==DLSODPK.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
        SUBROUTINE DLSODPK(F,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,    &
       &                   Iopt,Rwork,Lrw,Iwork,Liw,JAC,PSOL,Mf)
        IMPLICIT NONE
  !*--DLSODPK3613
  !*** Start of declarations inserted by SPAG
        !INTEGER JAC
        !REAL PSOL
  !*** End of declarations inserted by SPAG
        EXTERNAL F , JAC , PSOL
        INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
       &        Mf
        DOUBLE PRECISION Y , T , Tout , Rtol , Atol , Rwork
        DIMENSION Neq(*) , Y(*) , Rtol(*) , Atol(*) , Rwork(Lrw) ,        &
       &          Iwork(Liw)
  !-----------------------------------------------------------------------
  ! This is the 18 November 2003 version of
  ! DLSODPK: Livermore Solver for Ordinary Differential equations,
  !          with Preconditioned Krylov iteration methods for the
  !          Newton correction linear systems.
  !
  ! This version is in double precision.
  !
  ! DLSODPK solves the initial value problem for stiff or nonstiff
  ! systems of first order ODEs,
  !     dy/dt = f(t,y) ,  or, in component form,
  !     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
  !-----------------------------------------------------------------------
  ! Introduction.
  !
  ! This is a modification of the DLSODE package which incorporates
  ! various preconditioned Krylov subspace iteration methods for the
  ! linear algebraic systems that arise in the case of stiff systems.
  !
  ! The linear systems that must be solved have the form
  !   A * x  = b ,  where  A = identity - hl0 * (df/dy) .
  ! Here hl0 is a scalar, and df/dy is the Jacobian matrix of partial
  ! derivatives of f (NEQ by NEQ).
  !
  ! The particular Krylov method is chosen by setting the second digit,
  ! MITER, in the method flag MF.
  ! Currently, the values of MITER have the following meanings:
  !
  !  MITER = 1 means the preconditioned Scaled Incomplete
  !            Orthogonalization Method (SPIOM).
  !
  !          2 means an incomplete version of the Preconditioned Scaled
  !            Generalized Minimal Residual method (SPIGMR).
  !            This is the best choice in general.
  !
  !          3 means the Preconditioned Conjugate Gradient method (PCG).
  !            Recommended only when df/dy is symmetric or nearly so.
  !
  !          4 means the scaled Preconditioned Conjugate Gradient method
  !            (PCGS).  Recommended only when D-inverse * df/dy * D is
  !            symmetric or nearly so, where D is the diagonal scaling
  !            matrix with elements 1/EWT(i) (see RTOL/ATOL description).
  !
  !          9 means that only a user-supplied matrix P (approximating A)
  !            will be used, with no Krylov iteration done.  This option
  !            allows the user to provide the complete linear system
  !            solution algorithm, if desired.
  !
  ! The user can apply preconditioning to the linear system A*x = b,
  ! by means of arbitrary matrices (the preconditioners).
  !     In the case of SPIOM and SPIGMR, one can apply left and right
  ! preconditioners P1 and P2, and the basic iterative method is then
  ! applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the
  ! matrix A.  The product P1*P2 should be an approximation to matrix A
  ! such that linear systems with P1 or P2 are easier to solve than with
  ! A.  Preconditioning from the left only or right only means using
  ! P2 = identity or P1 = identity, respectively.
  !     In the case of the PCG and PCGS methods, there is only one
  ! preconditioner matrix P (but it can be the product of more than one).
  ! It should approximate the matrix A but allow for relatively
  ! easy solution of linear systems with coefficient matrix P.
  ! For PCG, P should be positive definite symmetric, or nearly so,
  ! and for PCGS, the scaled preconditioner D-inverse * P * D
  ! should be symmetric or nearly so.
  !     If the Jacobian J = df/dy splits in a natural way into a sum
  ! J = J1 + J2, then one possible choice of preconditioners is
  !     P1 = identity - hl0 * J1  and  P2 = identity - hl0 * J2
  ! provided each of these is easy to solve (or approximately solve).
  !
  !-----------------------------------------------------------------------
  ! References:
  ! 1.  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage Matrix
  !     Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31 (1989),
  !     pp. 40-91; also  L.L.N.L. Report UCRL-95088, Rev. 1, June 1987.
  ! 2.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
  !     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  !     North-Holland, Amsterdam, 1983, pp. 55-64.
  !-----------------------------------------------------------------------
  ! Authors:       Alan C. Hindmarsh and Peter N. Brown
  !                Center for Applied Scientific Computing, L-561
  !                Lawrence Livermore National Laboratory
  !                Livermore, CA 94551
  !-----------------------------------------------------------------------
  ! Summary of Usage.
  !
  ! Communication between the user and the DLSODPK package, for normal
  ! situations, is summarized here.  This summary describes only a subset
  ! of the full set of options available.  See the full description for
  ! details, including optional communication, nonstandard options,
  ! and instructions for special situations.  See also the demonstration
  ! program distributed with this solver.
  !
  ! A. First provide a subroutine of the form:
  !               SUBROUTINE F (NEQ, T, Y, YDOT)
  !               DOUBLE PRECISION T, Y(*), YDOT(*)
  ! which supplies the vector function f by loading YDOT(i) with f(i).
  !
  ! B. Next determine (or guess) whether or not the problem is stiff.
  ! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
  ! whose real part is negative and large in magnitude, compared to the
  ! reciprocal of the t span of interest.  If the problem is nonstiff,
  ! use a method flag MF = 10.  If it is stiff, MF should be between 21
  ! and 24, or possibly 29.  MF = 22 is generally the best choice.
  ! Use 23 or 24 only if symmetry is present.  Use MF = 29 if the
  ! complete linear system solution is to be provided by the user.
  ! The following four parameters must also be set.
  !  IWORK(1) = LWP  = length of real array WP for preconditioning.
  !  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
  !  IWORK(3) = JPRE = preconditioner type flag:
  !                  = 0 for no preconditioning (P1 = P2 = P = identity)
  !                  = 1 for left-only preconditioning (P2 = identity)
  !                  = 2 for right-only preconditioning (P1 = identity)
  !                  = 3 for two-sided preconditioning (and PCG or PCGS)
  !  IWORK(4) = JACFLG = flag for whether JAC is called.
  !                    = 0 if JAC is not to be called,
  !                    = 1 if JAC is to be called.
  !  Use JACFLG = 1 if JAC computes any nonconstant data for use in
  !  preconditioning, such as Jacobian elements.
  !  The arrays WP and IWP are work arrays under the user's control,
  !  for use in the routines that perform preconditioning operations.
  !
  ! C. If the problem is stiff, you must supply two routines that deal
  ! with the preconditioning of the linear systems to be solved.
  ! These are as follows:
  !
  !     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HL0, WP,IWP, IER)
  !     DOUBLE PRECISION T, Y(*),YSV(*), REWT(*), FTY(*), V(*), HL0, WP(*)
  !     INTEGER IWP(*)
  !        This routine must evaluate and preprocess any parts of the
  !     Jacobian matrix df/dy involved in the preconditioners P1, P2, P.
  !     The Y and FTY arrays contain the current values of y and f(t,y),
  !     respectively, and YSV also contains the current value of y.
  !     The array V is work space of length NEQ.
  !     JAC must multiply all computed Jacobian elements by the scalar
  !     -HL0, add the identity matrix, and do any factorization
  !     operations called for, in preparation for solving linear systems
  !     with a coefficient matrix of P1, P2, or P.  The matrix P1*P2 or P
  !     should be an approximation to  identity - HL0 * (df/dy).
  !     JAC should return IER = 0 if successful, and IER .ne. 0 if not.
  !     (If IER .ne. 0, a smaller time step will be tried.)
  !
  !     SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HL0, WP, IWP, B, LR, IER)
  !     DOUBLE PRECISION T, Y(*), FTY(*), WK(*), HL0, WP(*), B(*)
  !     INTEGER IWP(*)
  !        This routine must solve a linear system with B as right-hand
  !     side and one of the preconditioning matrices, P1, P2, or P, as
  !     coefficient matrix, and return the solution vector in B.
  !     LR is a flag concerning left vs right preconditioning, input
  !     to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
  !     In the case of the PCG or PCGS method, LR will be 3, and PSOL
  !     should solve the system P*x = B with the preconditioner matrix P.
  !     In the case MF = 29 (no Krylov iteration), LR will be 0,
  !     and PSOL is to return in B the desired approximate solution
  !     to A * x = B, where A = identity - HL0 * (df/dy).
  !     PSOL can use data generated in the JAC routine and stored in
  !     WP and IWP.  WK is a work array of length NEQ.
  !     The argument HL0 is the current value of the scalar appearing
  !     in the linear system.  If the old value, at the time of the last
  !     JAC call, is needed, it must have been saved by JAC in WP.
  !     On return, PSOL should set the error flag IER as follows:
  !       IER = 0 if PSOL was successful,
  !       IER .gt. 0 if a recoverable error occurred, meaning that the
  !              time step will be retried,
  !       IER .lt. 0 if an unrecoverable error occurred, meaning that the
  !              solver is to stop immediately.
  !
  ! D. Write a main program which calls Subroutine DLSODPK once for
  ! each point at which answers are desired.  This should also provide
  ! for possible use of logical unit 6 for output of error messages by
  ! DLSODPK.  On the first call to DLSODPK, supply arguments as follows:
  ! F      = name of subroutine for right-hand side vector f.
  !          This name must be declared External in calling program.
  ! NEQ    = number of first order ODEs.
  ! Y      = array of initial values, of length NEQ.
  ! T      = the initial value of the independent variable.
  ! TOUT   = first point where output is desired (.ne. T).
  ! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  ! RTOL   = relative tolerance parameter (scalar).
  ! ATOL   = absolute tolerance parameter (scalar or array).
  !          the estimated local error in y(i) will be controlled so as
  !          to be roughly less (in magnitude) than
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
  !          Thus the local error test passes if, in each component,
  !          either the absolute error is less than ATOL (or ATOL(i)),
  !          or the relative error is less than RTOL.
  !          Use RTOL = 0.0 for pure absolute error control, and
  !          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
  !          control.  Caution: Actual (global) errors may exceed these
  !          local tolerances, so choose them conservatively.
  ! ITASK  = 1 for normal computation of output values of y at t = TOUT.
  ! ISTATE = integer flag (input and output).  Set ISTATE = 1.
  ! IOPT   = 0 to indicate no optional inputs used.
  ! RWORK  = real work array of length at least:
  !             20 + 16*NEQ           for MF = 10,
  !             45 + 17*NEQ + LWP     for MF = 21,
  !             61 + 17*NEQ + LWP     for MF = 22,
  !             20 + 15*NEQ + LWP     for MF = 23 or 24,
  !             20 + 12*NEQ + LWP     for MF = 29.
  ! LRW    = declared length of RWORK (in user's dimension).
  ! IWORK  = integer work array of length at least:
  !             30            for MF = 10,
  !             35 + LIWP     for MF = 21,
  !             30 + LIWP     for MF = 22, 23, 24, or 29.
  ! LIW    = declared length of IWORK (in user's dimension).
  ! JAC,PSOL = names of subroutines for preconditioning.
  !          These names must be declared External in the calling program.
  ! MF     = method flag.  Standard values are:
  !          10 for nonstiff (Adams) method.
  !          21 for stiff (BDF) method, with preconditioned SIOM.
  !          22 for stiff method, with preconditioned GMRES method.
  !          23 for stiff method, with preconditioned CG method.
  !          24 for stiff method, with scaled preconditioned CG method.
  !          29 for stiff method, with user's PSOL routine only.
  ! Note that the main program must declare arrays Y, RWORK, IWORK,
  ! and possibly ATOL.
  !
  ! E. The output from the first call (or any call) is:
  !      Y = array of computed values of y(t) vector.
  !      T = corresponding value of independent variable (normally TOUT).
  ! ISTATE = 2  if DLSODPK was successful, negative otherwise.
  !          -1 means excess work done on this call (perhaps wrong MF).
  !          -2 means excess accuracy requested (tolerances too small).
  !          -3 means illegal input detected (see printed message).
  !          -4 means repeated error test failures (check all inputs).
  !          -5 means repeated convergence failures (perhaps bad JAC
  !             or PSOL routine supplied or wrong choice of MF or
  !             tolerances, or this solver is inappropriate).
  !          -6 means error weight became zero during problem. (Solution
  !             component i vanished, and ATOL or ATOL(i) = 0.)
  !          -7 means an unrecoverable error occurred in PSOL.
  !
  ! F. To continue the integration after a successful return, simply
  ! reset TOUT and call DLSODPK again.  No other parameters need be reset.
  !
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Full Description of User Interface to DLSODPK.
  !
  ! The user interface to DLSODPK consists of the following parts.
  !
  ! 1.   The call sequence to Subroutine DLSODPK, which is a driver
  !      routine for the solver.  This includes descriptions of both
  !      the call sequence arguments and of user-supplied routines.
  !      Following these descriptions is a description of
  !      optional inputs available through the call sequence, and then
  !      a description of optional outputs (in the work arrays).
  !
  ! 2.   Descriptions of other routines in the DLSODPK package that may be
  !      (optionally) called by the user.  These provide the ability to
  !      alter error message handling, save and restore the internal
  !      Common, and obtain specified derivatives of the solution y(t).
  !
  ! 3.   Descriptions of Common blocks to be declared in overlay
  !      or similar environments, or to be saved when doing an interrupt
  !      of the problem and continued solution later.
  !
  ! 4.   Description of two routines in the DLSODPK package, either of
  !      which the user may replace with his/her own version, if desired.
  !      These relate to the measurement of errors.
  !
  !-----------------------------------------------------------------------
  ! Part 1.  Call Sequence.
  !
  ! The call sequence parameters used for input only are
  !  F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, PSOL, MF,
  ! and those used for both input and output are
  !  Y, T, ISTATE.
  ! The work arrays RWORK and IWORK are also used for conditional and
  ! optional inputs and optional outputs.  (The term output here refers
  ! to the return from Subroutine DLSODPK to the user's calling program.)
  !
  ! The legality of input parameters will be thoroughly checked on the
  ! initial call for the problem, but not checked thereafter unless a
  ! change in input parameters is flagged by ISTATE = 3 on input.
  !
  ! The descriptions of the call arguments are as follows.
  !
  ! F      = the name of the user-supplied subroutine defining the
  !          ODE system.  The system must be put in the first-order
  !          form dy/dt = f(t,y), where f is a vector-valued function
  !          of the scalar t and the vector y.  Subroutine F is to
  !          compute the function f.  It is to have the form
  !               SUBROUTINE F (NEQ, T, Y, YDOT)
  !               DOUBLE PRECISION T, Y(*), YDOT(*)
  !          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
  !          is output.  Y and YDOT are arrays of length NEQ.
  !          Subroutine F should not alter Y(1),...,Y(NEQ).
  !          F must be declared External in the calling program.
  !
  !          Subroutine F may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in F) and/or Y has length exceeding NEQ(1).
  !          See the descriptions of NEQ and Y below.
  !
  !          If quantities computed in the F routine are needed
  !          externally to DLSODPK, an extra call to F should be made
  !          for this purpose, for consistent and accurate results.
  !          If only the derivative dy/dt is needed, use DINTDY instead.
  !
  ! NEQ    = the size of the ODE system (number of first order
  !          ordinary differential equations).  Used only for input.
  !          NEQ may be decreased, but not increased, during the problem.
  !          If NEQ is decreased (with ISTATE = 3 on input), the
  !          remaining components of Y should be left undisturbed, if
  !          these are to be accessed in the user-supplied subroutines.
  !
  !          Normally, NEQ is a scalar, and it is generally referred to
  !          as a scalar in this user interface description.  However,
  !          NEQ may be an array, with NEQ(1) set to the system size.
  !          (The DLSODPK package accesses only NEQ(1).)  In either case,
  !          this parameter is passed as the NEQ argument in all calls
  !          to F, JAC, and PSOL.  Hence, if it is an array, locations
  !          NEQ(2),... may be used to store other integer data and pass
  !          it to the user-supplied subroutines.  Each such routine must
  !          include NEQ in a Dimension statement in that case.
  !
  ! Y      = a real array for the vector of dependent variables, of
  !          length NEQ or more.  Used for both input and output on the
  !          first call (ISTATE = 1), and only for output on other calls.
  !          On the first call, Y must contain the vector of initial
  !          values.  On output, Y contains the computed solution vector,
  !          evaluated at T.  If desired, the Y array may be used
  !          for other purposes between calls to the solver.
  !
  !          This array is passed as the Y argument in all calls to F,
  !          JAC, and PSOL. Hence its length may exceed NEQ, and locations
  !          Y(NEQ+1),... may be used to store other real data and
  !          pass it to the user-supplied subroutines.  (The DLSODPK
  !          package accesses only Y(1),...,Y(NEQ).)
  !
  ! T      = the independent variable.  On input, T is used only on the
  !          first call, as the initial point of the integration.
  !          On output, after each call, T is the value at which a
  !          computed solution y is evaluated (usually the same as TOUT).
  !          On an error return, T is the farthest point reached.
  !
  ! TOUT   = the next value of t at which a computed solution is desired.
  !          Used only for input.
  !
  !          When starting the problem (ISTATE = 1), TOUT may be equal
  !          to T for one call, then should .ne. T for the next call.
  !          For the initial T, an input value of TOUT .ne. T is used
  !          in order to determine the direction of the integration
  !          (i.e. the algebraic sign of the step sizes) and the rough
  !          scale of the problem.  Integration in either direction
  !          (forward or backward in t) is permitted.
  !
  !          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
  !          the first call (i.e. the first call with TOUT .ne. T).
  !          Otherwise, TOUT is required on every call.
  !
  !          If ITASK = 1, 3, or 4, the values of TOUT need not be
  !          monotone, but a value of TOUT which backs up is limited
  !          to the current internal T interval, whose endpoints are
  !          TCUR - HU and TCUR (see optional outputs, below, for
  !          TCUR and HU).
  !
  ! ITOL   = an indicator for the type of error control.  See
  !          description below under ATOL.  Used only for input.
  !
  ! RTOL   = a relative error tolerance parameter, either a scalar or
  !          an array of length NEQ.  See description below under ATOL.
  !          Input only.
  !
  ! ATOL   = an absolute error tolerance parameter, either a scalar or
  !          an array of length NEQ.  Input only.
  !
  !             The input parameters ITOL, RTOL, and ATOL determine
  !          the error control performed by the solver.  The solver will
  !          control the vector E = (E(i)) of estimated local errors
  !          in y, according to an inequality of the form
  !                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
  !          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
  !          and the RMS-norm (root-mean-square norm) here is
  !          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
  !          is a vector of weights which must always be positive, and
  !          the values of RTOL and ATOL should all be non-negative.
  !          the following table gives the types (scalar/array) of
  !          RTOL and ATOL, and the corresponding form of EWT(i).
  !
  !             ITOL    RTOL       ATOL          EWT(i)
  !              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  !              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  !              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  !              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
  !
  !          When either of these parameters is a scalar, it need not
  !          be dimensioned in the user's calling program.
  !
  !          If none of the above choices (with ITOL, RTOL, and ATOL
  !          fixed throughout the problem) is suitable, more general
  !          error controls can be obtained by substituting
  !          user-supplied routines for the setting of EWT and/or for
  !          the norm calculation.  See Part 4 below.
  !
  !          If global errors are to be estimated by making a repeated
  !          run on the same problem with smaller tolerances, then all
  !          components of RTOL and ATOL (i.e. of EWT) should be scaled
  !          down uniformly.
  !
  ! ITASK  = an index specifying the task to be performed.
  !          Input only.  ITASK has the following values and meanings.
  !          1  means normal computation of output values of y(t) at
  !             t = TOUT (by overshooting and interpolating).
  !          2  means take one step only and return.
  !          3  means stop at the first internal mesh point at or
  !             beyond t = TOUT and return.
  !          4  means normal computation of output values of y(t) at
  !             t = TOUT but without overshooting t = TCRIT.
  !             TCRIT must be input as RWORK(1).  TCRIT may be equal to
  !             or beyond TOUT, but not behind it in the direction of
  !             integration.  This option is useful if the problem
  !             has a singularity at or beyond t = TCRIT.
  !          5  means take one step, without passing TCRIT, and return.
  !             TCRIT must be input as RWORK(1).
  !
  !          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
  !          (within roundoff), it will return T = TCRIT (exactly) to
  !          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
  !          in which case answers at t = TOUT are returned first).
  !
  ! ISTATE = an index used for input and output to specify the
  !          the state of the calculation.
  !
  !          On input, the values of ISTATE are as follows.
  !          1  means this is the first call for the problem
  !             (initializations will be done).  See note below.
  !          2  means this is not the first call, and the calculation
  !             is to continue normally, with no change in any input
  !             parameters except possibly TOUT and ITASK.
  !             (If ITOL, RTOL, and/or ATOL are changed between calls
  !             with ISTATE = 2, the new values will be used but not
  !             tested for legality.)
  !          3  means this is not the first call, and the
  !             calculation is to continue normally, but with
  !             a change in input parameters other than
  !             TOUT and ITASK.  Changes are allowed in
  !             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
  !             and any of the optional inputs except H0.
  !          Note:  A preliminary call with TOUT = T is not counted
  !          as a first call here, as no initialization or checking of
  !          input is done.  (Such a call is sometimes useful for the
  !          purpose of outputting the initial conditions.)
  !          Thus the first call for which TOUT .ne. T requires
  !          ISTATE = 1 on input.
  !
  !          On output, ISTATE has the following values and meanings.
  !           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
  !           2  means the integration was performed successfully.
  !          -1  means an excessive amount of work (more than MXSTEP
  !              steps) was done on this call, before completing the
  !              requested task, but the integration was otherwise
  !              successful as far as T.  (MXSTEP is an optional input
  !              and is normally 500.)  To continue, the user may
  !              simply reset ISTATE to a value .gt. 1 and call again
  !              (the excess work step counter will be reset to 0).
  !              In addition, the user may increase MXSTEP to avoid
  !              this error return (see below on optional inputs).
  !          -2  means too much accuracy was requested for the precision
  !              of the machine being used.  This was detected before
  !              completing the requested task, but the integration
  !              was successful as far as T.  To continue, the tolerance
  !              parameters must be reset, and ISTATE must be set
  !              to 3.  The optional output TOLSF may be used for this
  !              purpose.  (Note: If this condition is detected before
  !              taking any steps, then an illegal input return
  !              (ISTATE = -3) occurs instead.)
  !          -3  means illegal input was detected, before taking any
  !              integration steps.  See written message for details.
  !              Note:  If the solver detects an infinite loop of calls
  !              to the solver with illegal input, it will cause
  !              the run to stop.
  !          -4  means there were repeated error test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              The problem may have a singularity, or the input
  !              may be inappropriate.
  !          -5  means there were repeated convergence test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !          -6  means EWT(i) became zero for some i during the
  !              integration.  Pure relative error control (ATOL(i)=0.0)
  !              was requested on a variable which has now vanished.
  !              The integration was successful as far as T.
  !          -7  means the PSOL routine returned an unrecoverable error
  !              flag (IER .lt. 0).  The integration was successful as
  !              far as T.
  !
  !          Note:  since the normal output value of ISTATE is 2,
  !          it does not need to be reset for normal continuation.
  !          Also, since a negative input value of ISTATE will be
  !          regarded as illegal, a negative output value requires the
  !          user to change it, and possibly other inputs, before
  !          calling the solver again.
  !
  ! IOPT   = an integer flag to specify whether or not any optional
  !          inputs are being used on this call.  Input only.
  !          The optional inputs are listed separately below.
  !          IOPT = 0 means no optional inputs are being used.
  !                   Default values will be used in all cases.
  !          IOPT = 1 means one or more optional inputs are being used.
  !
  ! RWORK  = a real working array (double precision).
  !          The length of RWORK must be at least
  !             20 + NYH*(MAXORD + 1) + 3*NEQ + LENLS + LWP    where
  !          NYH    = the initial value of NEQ,
  !          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
  !                   smaller value is given as an optional input),
  !          LENLS = length of work space for linear system (Krylov)
  !                  method, excluding preconditioning:
  !            LENLS = 0                               if MITER = 0,
  !            LENLS = NEQ*(MAXL+3) + MAXL**2          if MITER = 1,
  !            LENLS = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
  !                 + (MAXL+3)*MAXL + 1                if MITER = 2,
  !            LENLS = 6*NEQ                           if MITER = 3 or 4,
  !            LENLS = 3*NEQ                           if MITER = 9.
  !          (See the MF description for METH and MITER, and the
  !          list of optional inputs for MAXL and KMP.)
  !          LWP = length of real user work space for preconditioning
  !          (see JAC/PSOL).
  !          Thus if default values are used and NEQ is constant,
  !          this length is:
  !             20 + 16*NEQ           for MF = 10,
  !             45 + 24*NEQ + LWP     FOR MF = 11,
  !             61 + 24*NEQ + LWP     FOR MF = 12,
  !             20 + 22*NEQ + LWP     FOR MF = 13 OR 14,
  !             20 + 19*NEQ + LWP     FOR MF = 19,
  !             20 + 9*NEQ            FOR MF = 20,
  !             45 + 17*NEQ + LWP     FOR MF = 21,
  !             61 + 17*NEQ + LWP     FOR MF = 22,
  !             20 + 15*NEQ + LWP     FOR MF = 23 OR 24,
  !             20 + 12*NEQ + LWP     for MF = 29.
  !          The first 20 words of RWORK are reserved for conditional
  !          and optional inputs and optional outputs.
  !
  !          The following word in RWORK is a conditional input:
  !            RWORK(1) = TCRIT = critical value of t which the solver
  !                       is not to overshoot.  Required if ITASK is
  !                       4 or 5, and ignored otherwise.  (See ITASK.)
  !
  ! LRW    = the length of the array RWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! IWORK  = an integer work array.  The length of IWORK must be at least
  !             30                 if MITER = 0 (MF = 10 or 20),
  !             30 + MAXL + LIWP   if MITER = 1 (MF = 11, 21),
  !             30 + LIWP          if MITER = 2, 3, 4, or 9.
  !          MAXL = 5 unless a different optional input value is given.
  !          LIWP = length of integer user work space for preconditioning
  !          (see conditional input list following).
  !          The first few words of IWORK are used for conditional and
  !          optional inputs and optional outputs.
  !
  !          The following 4 words in IWORK are conditional inputs,
  !          required if MITER .ge. 1:
  !          IWORK(1) = LWP  = length of real array WP for use in
  !                     preconditioning (part of RWORK array).
  !          IWORK(2) = LIWP = length of integer array IWP for use in
  !                     preconditioning (part of IWORK array).
  !                     The arrays WP and IWP are work arrays under the
  !                     user's control, for use in the routines that
  !                     perform preconditioning operations (JAC and PSOL).
  !          IWORK(3) = JPRE = preconditioner type flag:
  !                   = 0 for no preconditioning (P1 = P2 = P = identity)
  !                   = 1 for left-only preconditioning (P2 = identity)
  !                   = 2 for right-only preconditioning (P1 = identity)
  !                   = 3 for two-sided preconditioning (and PCG or PCGS)
  !          IWORK(4) = JACFLG = flag for whether JAC is called.
  !                   = 0 if JAC is not to be called,
  !                   = 1 if JAC is to be called.
  !                     Use JACFLG = 1 if JAC computes any nonconstant
  !                     data needed in preconditioning operations,
  !                     such as some of the Jacobian elements.
  !
  ! LIW    = the length of the array IWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! Note:  The work arrays must not be altered between calls to DLSODPK
  ! for the same problem, except possibly for the conditional and
  ! optional inputs, and except for the last 3*NEQ words of RWORK.
  ! The latter space is used for internal scratch space, and so is
  ! available for use by the user outside DLSODPK between calls, if
  ! desired (but not for use by any of the user-supplied subroutines).
  !
  ! JAC    = the name of the user-supplied routine to compute any
  !          Jacobian elements (or approximations) involved in the
  !          matrix preconditioning operations (MITER .ge. 1).
  !          It is to have the form
  !            SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V,
  !           1                HL0, WP, IWP, IER)
  !            DOUBLE PRECISION T, Y(*),YSV(*), REWT(*), FTY(*), V(*),
  !           1                 HL0, WP(*)
  !            INTEGER IWP(*)
  !          This routine must evaluate and preprocess any parts of the
  !          Jacobian matrix df/dy used in the preconditioners P1, P2, P.
  !          the Y and FTY arrays contain the current values of y and
  !          f(t,y), respectively, and YSV also contains the current
  !          value of y.  The array V is work space of length
  !          NEQ for use by JAC.  REWT is the array of reciprocal error
  !          weights (1/EWT).  JAC must multiply all computed Jacobian
  !          elements by the scalar -HL0, add the identity matrix, and do
  !          any factorization operations called for, in preparation
  !          for solving linear systems with a coefficient matrix of
  !          P1, P2, or P.  The matrix P1*P2 or P should be an
  !          approximation to  identity - HL0 * (df/dy).  JAC should
  !          return IER = 0 if successful, and IER .ne. 0 if not.
  !          (If IER .ne. 0, a smaller time step will be tried.)
  !          The arrays WP (of length LWP) and IWP (of length LIWP)
  !          are for use by JAC and PSOL for work space and for storage
  !          of data needed for the solution of the preconditioner
  !          linear systems.  Their lengths and contents are under the
  !          user's control.
  !          The JAC routine may save relevant Jacobian elements (or
  !          approximations) used in the preconditioners, along with the
  !          value of HL0, and use these to reconstruct preconditioner
  !          matrices later without reevaluationg those elements.
  !          This may be cost-effective if JAC is called with HL0
  !          considerably different from its earlier value, indicating
  !          that a corrector convergence failure has occurred because
  !          of the change in HL0, not because of changes in the
  !          value of the Jacobian.  In doing this, use the saved and
  !          current values of HL0 to decide whether to use saved
  !          or reevaluated elements.
  !          JAC may alter V, but may not alter Y, YSV, REWT, FTY, or HL0.
  !          JAC must be declared External in the calling program.
  !               Subroutine JAC may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
  !          See the descriptions of NEQ and Y above.
  !
  ! PSOL   = the name of the user-supplied routine for the
  !          solution of preconditioner linear systems.
  !          It is to have the form
  !            SUBROUTINE PSOL (NEQ, T, Y, FTY, WK,HL0, WP,IWP, B, LR,IER)
  !            DOUBLE PRECISION T, Y(*), FTY(*), WK(*), HL0, WP(*), B(*)
  !            INTEGER IWP(*)
  !          This routine must solve a linear system with B as right-hand
  !          side and one of the preconditioning matrices, P1, P2, or P,
  !          as coefficient matrix, and return the solution vector in B.
  !          LR is a flag concerning left vs right preconditioning, input
  !          to PSOL.  PSOL is to use P1 if LR = 1 and P2 if LR = 2.
  !          In the case of the PCG or PCGS method, LR will be 3, and PSOL
  !          should solve the system P*x = B with the preconditioner P.
  !          In the case MITER = 9 (no Krylov iteration), LR will be 0,
  !          and PSOL is to return in B the desired approximate solution
  !          to A * x = B, where A = identity - HL0 * (df/dy).
  !          PSOL can use data generated in the JAC routine and stored in
  !          WP and IWP.
  !          The Y and FTY arrays contain the current values of y and
  !          f(t,y), respectively.  The array WK is work space of length
  !          NEQ for use by PSOL.
  !          The argument HL0 is the current value of the scalar appearing
  !          in the linear system.  If the old value, as of the last
  !          JAC call, is needed, it must have been saved by JAC in WP.
  !          On return, PSOL should set the error flag IER as follows:
  !            IER = 0 if PSOL was successful,
  !            IER .gt. 0 on a recoverable error, meaning that the
  !                   time step will be retried,
  !            IER .lt. 0 on an unrecoverable error, meaning that the
  !                   solver is to stop immediately.
  !          PSOL may not alter Y, FTY, or HL0.
  !          PSOL must be declared External in the calling program.
  !               Subroutine PSOL may access user-defined quantities in
  !          NEQ(2),... and Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in PSOL) and/or Y has length exceeding NEQ(1).
  !          See the descriptions of NEQ and Y above.
  !
  ! MF     = the method flag.  Used only for input.  The legal values of
  !          MF are 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, and 29.
  !          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
  !          METH indicates the basic linear multistep method:
  !            METH = 1 means the implicit Adams method.
  !            METH = 2 means the method based on Backward
  !                     Differentiation Formulas (BDFs).
  !          MITER indicates the corrector iteration method:
  !            MITER = 0 means functional iteration (no linear system
  !                      is involved).
  !            MITER = 1 means Newton iteration with Scaled Preconditioned
  !                      Incomplete Orthogonalization Method (SPIOM)
  !                      for the linear systems.
  !            MITER = 2 means Newton iteration with Scaled Preconditioned
  !                      Generalized Minimal Residual method (SPIGMR)
  !                      for the linear systems.
  !            MITER = 3 means Newton iteration with Preconditioned
  !                      Conjugate Gradient method (PCG)
  !                      for the linear systems.
  !            MITER = 4 means Newton iteration with scaled Preconditioned
  !                      Conjugate Gradient method (PCGS)
  !                      for the linear systems.
  !            MITER = 9 means Newton iteration with only the
  !                      user-supplied PSOL routine called (no Krylov
  !                      iteration) for the linear systems.
  !                      JPRE is ignored, and PSOL is called with LR = 0.
  !          See comments in the introduction about the choice of MITER.
  !          If MITER .ge. 1, the user must supply routines JAC and PSOL
  !          (the names are arbitrary) as described above.
  !          For MITER = 0, dummy arguments can be used.
  !-----------------------------------------------------------------------
  ! Optional Inputs.
  !
  ! The following is a list of the optional inputs provided for in the
  ! call sequence.  (See also Part 2.)  For each such input variable,
  ! this table lists its name as used in this documentation, its
  ! location in the call sequence, its meaning, and the default value.
  ! The use of any of these inputs requires IOPT = 1, and in that
  ! case all of these inputs are examined.  A value of zero for any
  ! of these optional inputs will cause the default value to be used.
  ! Thus to use a subset of the optional inputs, simply preload
  ! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
  ! then set those of interest to nonzero values.
  !
  ! Name    Location      Meaning and Default Value
  !
  ! H0      RWORK(5)  the step size to be attempted on the first step.
  !                   The default value is determined by the solver.
  !
  ! HMAX    RWORK(6)  the maximum absolute step size allowed.
  !                   The default value is infinite.
  !
  ! HMIN    RWORK(7)  the minimum absolute step size allowed.
  !                   The default value is 0.  (This lower bound is not
  !                   enforced on the final step before reaching TCRIT
  !                   when ITASK = 4 or 5.)
  !
  ! DELT    RWORK(8)  convergence test constant in Krylov iteration
  !                   algorithm.  The default is .05.
  !
  ! MAXORD  IWORK(5)  the maximum order to be allowed.  The default
  !                   value is 12 if METH = 1, and 5 if METH = 2.
  !                   If MAXORD exceeds the default value, it will
  !                   be reduced to the default value.
  !                   If MAXORD is changed during the problem, it may
  !                   cause the current order to be reduced.
  !
  ! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
  !                   allowed during one call to the solver.
  !                   The default value is 500.
  !
  ! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
  !                   warning that T + H = T on a step (H = step size).
  !                   This must be positive to result in a non-default
  !                   value.  The default value is 10.
  !
  ! MAXL    IWORK(8)  maximum number of iterations in the SPIOM, SPIGMR,
  !                   PCG, or PCGS algorithm (.le. NEQ).
  !                   The default is MAXL = MIN(5,NEQ).
  !
  ! KMP     IWORK(9)  number of vectors on which orthogonalization
  !                   is done in SPIOM or SPIGMR algorithm (.le. MAXL).
  !                   The default is KMP = MAXL.
  !                   Note:  When KMP .lt. MAXL and MF = 22, the length
  !                          of RWORK must be defined accordingly.  See
  !                          the definition of RWORK above.
  !-----------------------------------------------------------------------
  ! Optional Outputs.
  !
  ! As optional additional output from DLSODPK, the variables listed
  ! below are quantities related to the performance of DLSODPK
  ! which are available to the user.  These are communicated by way of
  ! the work arrays, but also have internal mnemonic names as shown.
  ! Except where stated otherwise, all of these outputs are defined
  ! on any successful return from DLSODPK, and on any return with
  ! ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return
  ! (ISTATE = -3), they will be unchanged from their existing values
  ! (if any), except possibly for TOLSF, LENRW, and LENIW.
  ! On any error return, outputs relevant to the error will be defined,
  ! as noted below.
  !
  ! Name    Location      Meaning
  !
  ! HU      RWORK(11) the step size in t last used (successfully).
  !
  ! HCUR    RWORK(12) the step size to be attempted on the next step.
  !
  ! TCUR    RWORK(13) the current value of the independent variable
  !                   which the solver has actually reached, i.e. the
  !                   current internal mesh point in t.  On output, TCUR
  !                   will always be at least as far as the argument
  !                   T, but may be farther (if interpolation was done).
  !
  ! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
  !                   computed when a request for too much accuracy was
  !                   detected (ISTATE = -3 if detected at the start of
  !                   the problem, ISTATE = -2 otherwise).  If ITOL is
  !                   left unaltered but RTOL and ATOL are uniformly
  !                   scaled up by a factor of TOLSF for the next call,
  !                   then the solver is deemed likely to succeed.
  !                   (The user may also ignore TOLSF and alter the
  !                   tolerance parameters in any other way appropriate.)
  !
  ! NST     IWORK(11) the number of steps taken for the problem so far.
  !
  ! NFE     IWORK(12) the number of f evaluations for the problem so far.
  !
  ! NPE     IWORK(13) the number of calls to JAC so far (for Jacobian
  !                   evaluation associated with preconditioning).
  !
  ! NQU     IWORK(14) the method order last used (successfully).
  !
  ! NQCUR   IWORK(15) the order to be attempted on the next step.
  !
  ! IMXER   IWORK(16) the index of the component of largest magnitude in
  !                   the weighted local error vector ( E(i)/EWT(i) ),
  !                   on an error return with ISTATE = -4 or -5.
  !
  ! LENRW   IWORK(17) the length of RWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! LENIW   IWORK(18) the length of IWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! NNI     IWORK(19) number of nonlinear iterations so far (each of
  !                   which calls an iterative linear solver).
  !
  ! NLI     IWORK(20) number of linear iterations so far.
  !                   Note: A measure of the success of algorithm is
  !                   the average number of linear iterations per
  !                   nonlinear iteration, given by NLI/NNI.
  !                   If this is close to MAXL, MAXL may be too small.
  !
  ! NPS     IWORK(21) number of preconditioning solve operations
  !                   (PSOL calls) so far.
  !
  ! NCFN    IWORK(22) number of convergence failures of the nonlinear
  !                   (Newton) iteration so far.
  !                   Note: A measure of success is the overall
  !                   rate of nonlinear convergence failures, NCFN/NST.
  !
  ! NCFL    IWORK(23) number of convergence failures of the linear
  !                   iteration so far.
  !                   Note: A measure of success is the overall
  !                   rate of linear convergence failures, NCFL/NNI.
  !
  ! The following two arrays are segments of the RWORK array which
  ! may also be of interest to the user as optional outputs.
  ! For each array, the table below gives its internal name,
  ! its base address in RWORK, and its description.
  !
  ! Name    Base Address      Description
  !
  ! YH      21             the Nordsieck history array, of size NYH by
  !                        (NQCUR + 1), where NYH is the initial value
  !                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
  !                        of YH contains HCUR**j/factorial(j) times
  !                        the j-th derivative of the interpolating
  !                        polynomial currently representing the solution,
  !                        evaluated at t = TCUR.
  !
  ! ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
  !                        corrections on each step, scaled on output
  !                        to represent the estimated local error in y
  !                        on the last step.  This is the vector E in
  !                        the description of the error control.  It is
  !                        defined only on a successful return from
  !                        DLSODPK.
  !
  !-----------------------------------------------------------------------
  ! Part 2.  Other Routines Callable.
  !
  ! The following are optional calls which the user may make to
  ! gain additional capabilities in conjunction with DLSODPK.
  ! (The routines XSETUN and XSETF are designed to conform to the
  ! SLATEC error handling package.)
  !
  !     Form of Call                  Function
  !   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
  !                             output of messages from DLSODPK, if
  !                             the default is not desired.
  !                             The default value of lun is 6.
  !
  !   CALL XSETF(MFLAG)         Set a flag to control the printing of
  !                             messages by DLSODPK.
  !                             MFLAG = 0 means do not print. (Danger:
  !                             This risks losing valuable information.)
  !                             MFLAG = 1 means print (the default).
  !
  !                             Either of the above calls may be made at
  !                             any time and will take effect immediately.
  !
  !   CALL DSRCPK(RSAV,ISAV,JOB) saves and restores the contents of
  !                             the internal Common blocks used by
  !                             DLSODPK (see Part 3 below).
  !                             RSAV must be a real array of length 222
  !                             or more, and ISAV must be an integer
  !                             array of length 50 or more.
  !                             JOB=1 means save Common into RSAV/ISAV.
  !                             JOB=2 means restore Common from RSAV/ISAV.
  !                                DSRCPK is useful if one is
  !                             interrupting a run and restarting
  !                             later, or alternating between two or
  !                             more problems solved with DLSODPK.
  !
  !   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
  !        (See below)          orders, at a specified point t, if
  !                             desired.  It may be called only after
  !                             a successful return from DLSODPK.
  !
  ! The detailed instructions for using DINTDY are as follows.
  ! The form of the call is:
  !
  !   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
  !
  ! The input parameters are:
  !
  ! T         = value of independent variable where answers are desired
  !             (normally the same as the T last returned by DLSODPK).
  !             for valid results, T must lie between TCUR - HU and TCUR.
  !             (See optional outputs for TCUR and HU.)
  ! K         = integer order of the derivative desired.  K must satisfy
  !             0 .le. K .le. NQCUR, where NQCUR is the current order
  !             (see optional outputs).  The capability corresponding
  !             to K = 0, i.e. computing y(T), is already provided
  !             by DLSODPK directly.  Since NQCUR .ge. 1, the first
  !             derivative dy/dt is always available with DINTDY.
  ! RWORK(21) = the base address of the history array YH.
  ! NYH       = column length of YH, equal to the initial value of NEQ.
  !
  ! The output parameters are:
  !
  ! DKY       = a real array of length NEQ containing the computed value
  !             of the K-th derivative of y(t).
  ! IFLAG     = integer flag, returned as 0 if K and T were legal,
  !             -1 if K was illegal, and -2 if T was illegal.
  !             On an error return, a message is also written.
  !-----------------------------------------------------------------------
  ! Part 3.  Common Blocks.
  !
  ! If DLSODPK is to be used in an overlay situation, the user
  ! must declare, in the primary overlay, the variables in:
  !   (1) the call sequence to DLSODPK, and
  !   (2) the two internal Common blocks
  !         /DLS001/  of length  255  (218 double precision words
  !                      followed by 37 integer words),
  !         /DLPK01/  of length  17  (4 double precision words
  !                      followed by 13 integer words).
  !
  ! If DLSODPK is used on a system in which the contents of internal
  ! Common blocks are not preserved between calls, the user should
  ! declare the above Common blocks in the calling program to insure
  ! that their contents are preserved.
  !
  ! If the solution of a given problem by DLSODPK is to be interrupted
  ! and then later continued, such as when restarting an interrupted run
  ! or alternating between two or more problems, the user should save,
  ! following the return from the last DLSODPK call prior to the
  ! interruption, the contents of the call sequence variables and the
  ! internal Common blocks, and later restore these values before the
  ! next DLSODPK call for that problem.  To save and restore the Common
  ! blocks, use Subroutine DSRCPK (see Part 2 above).
  !
  !-----------------------------------------------------------------------
  ! Part 4.  Optionally Replaceable Solver Routines.
  !
  ! below are descriptions of two routines in the DLSODPK package which
  ! relate to the measurement of errors.  Either routine can be
  ! replaced by a user-supplied version, if desired.  However, since such
  ! a replacement may have a major impact on performance, it should be
  ! done only when absolutely necessary, and only with great caution.
  ! (Note: The means by which the package version of a routine is
  ! superseded by the user's version may be system-dependent.)
  !
  ! (a) DEWSET.
  ! The following subroutine is called just before each internal
  ! integration step, and sets the array of error weights, EWT, as
  ! described under ITOL/RTOL/ATOL above:
  !     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
  ! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODPK call sequence,
  ! YCUR contains the current dependent variable vector, and
  ! EWT is the array of weights set by DEWSET.
  !
  ! If the user supplies this subroutine, it must return in EWT(i)
  ! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
  ! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
  ! routine (see below), and also used by DLSODPK in the computation
  ! of the optional output IMXER, the diagonal Jacobian approximation,
  ! and the increments for difference quotient Jacobians.
  !
  ! In the user-supplied version of DEWSET, it may be desirable to use
  ! the current values of derivatives of y.  Derivatives up to order NQ
  ! are available from the history array YH, described above under
  ! optional outputs.  In DEWSET, YH is identical to the YCUR array,
  ! extended to NQ + 1 columns with a column length of NYH and scale
  ! factors of H**j/factorial(j).  On the first call for the problem,
  ! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
  ! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
  ! can be obtained by including in DEWSET the statements:
  !     DOUBLE PRECISION RLS
  !     COMMON /DLS001/ RLS(218),ILS(37)
  !     NQ = ILS(33)
  !     NST = ILS(34)
  !     H = RLS(212)
  ! Thus, for example, the current value of dy/dt can be obtained as
  ! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
  ! unnecessary when NST = 0).
  !
  ! (b) DVNORM.
  ! The following is a real function routine which computes the weighted
  ! root-mean-square norm of a vector v:
  !     D = DVNORM (N, V, W)
  ! where:
  !   N = the length of the vector,
  !   V = real array of length N containing the vector,
  !   W = real array of length N containing weights,
  !   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
  ! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
  ! EWT is as set by Subroutine DEWSET.
  !
  ! If the user supplies this function, it should return a non-negative
  ! value of DVNORM suitable for use in the error control in DLSODPK.
  ! None of the arguments should be altered by DVNORM.
  ! For example, a user-supplied DVNORM routine might:
  !   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
  !   -ignore some components of V in the norm, with the effect of
  !    suppressing the error control on those components of y.
  !-----------------------------------------------------------------------
  !
  !***REVISION HISTORY  (YYYYMMDD)
  ! 19860901  DATE WRITTEN
  ! 19861010  Numerous minor revisions to SPIOM and SPGMR routines;
  !           minor corrections to prologues and comments.
  ! 19870114  Changed name SPGMR to SPIGMR; revised residual norm
  !           calculation in SPIGMR (for incomplete case);
  !           revised error return logic in SPIGMR;
  ! 19870330  Major update: corrected comments throughout;
  !           removed TRET from Common; rewrote EWSET with 4 loops;
  !           fixed t test in INTDY; added Cray directives in STODPK;
  !           in STODPK, fixed DELP init. and logic around PJAC call;
  !           combined routines to save/restore Common;
  !           passed LEVEL = 0 in error message calls (except run abort).
  ! 19871130  Added option MITER = 9; shortened WM array by 2;
  !           revised early return from SPIOM and SPIGMR;
  !           replaced copy loops with SCOPY/DCOPY calls;
  !           minor corrections/revisions to SOLPK, SPIGMR, ATV, ATP;
  !           corrections to main prologue and internal comments.
  ! 19880304  Corrections to type declarations in SOLPK, SPIOM, USOL.
  ! 19891025  Added ISTATE = -7 return; minor revisions to USOL;
  !           added initialization of JACFLG in main driver;
  !           removed YH and NYH from PKSET call list;
  !           minor revisions to SPIOM and SPIGMR;
  !           corrections to main prologue and internal comments.
  ! 19900803  Added YSV to JAC call list; minor comment corrections.
  ! 20010425  Major update: convert source lines to upper case;
  !           added *DECK lines; changed from 1 to * in dummy dimensions;
  !           changed names R1MACH/D1MACH to RUMACH/DUMACH;
  !           renamed routines for uniqueness across single/double prec.;
  !           converted intrinsic names to generic form;
  !           removed ILLIN and NTREP (data loaded) from Common;
  !           removed all 'own' variables from Common;
  !           changed error messages to quoted strings;
  !           replaced XERRWV/XERRWD with 1993 revised version;
  !           converted prologues, comments, error messages to mixed case;
  !           numerous corrections to prologues and internal comments.
  ! 20010507  Converted single precision source to double precision.
  ! 20020502  Corrected declarations in descriptions of user routines.
  ! 20030603  Corrected duplicate type declaration for DUMACH.
  ! 20031105  Restored 'own' variables to Common blocks, to enable
  !           interrupt/restart feature.
  ! 20031112  Added SAVE statements for data-loaded constants.
  ! 20031117  Changed internal name NPE to NJE.
  !
  !-----------------------------------------------------------------------
  ! Other routines in the DLSODPK package.
  !
  ! In addition to Subroutine DLSODPK, the DLSODPK package includes the
  ! following subroutines and function routines:
  !  DINTDY   computes an interpolated value of the y vector at t = TOUT.
  !  DEWSET   sets the error weight vector EWT before each step.
  !  DVNORM   computes the weighted RMS-norm of a vector.
  !  DSTODPK  is the core integrator, which does one step of the
  !           integration and the associated error control.
  !  DCFODE   sets all method coefficients and test constants.
  !  DPKSET   interfaces between DSTODPK and the JAC routine.
  !  DSOLPK   manages solution of linear system in Newton iteration.
  !  DSPIOM   performs the SPIOM algorithm.
  !  DATV     computes a scaled, preconditioned product (I-hl0*J)*v.
  !  DORTHOG  orthogonalizes a vector against previous basis vectors.
  !  DHEFA    generates an LU factorization of a Hessenberg matrix.
  !  DHESL    solves a Hessenberg square linear system.
  !  DSPIGMR  performs the SPIGMR algorithm.
  !  DHEQR    generates a QR factorization of a Hessenberg matrix.
  !  DHELS    finds the least squares solution of a Hessenberg system.
  !  DPCG     performs Preconditioned Conjugate Gradient algorithm (PCG).
  !  DPCGS    performs the PCGS algorithm.
  !  DATP     computes the product A*p, where A = I - hl0*df/dy.
  !  DUSOL    interfaces to the user's PSOL routine (MITER = 9).
  !  DSRCPK   is a user-callable routine to save and restore
  !           the contents of the internal Common blocks.
  !  DAXPY, DCOPY, DDOT, DNRM2, and DSCAL   are basic linear
  !           algebra modules (from the BLAS collection).
  !  DUMACH   computes the unit roundoff in a machine-independent manner.
  !  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
  !           error messages and warnings.  XERRWD is machine-dependent.
  ! Note:  DVNORM, DDOT, DNRM2, DUMACH, IXSAV, and IUMACH are function
  ! routines.  All the others are subroutines.
  !
  !-----------------------------------------------------------------------
        !DOUBLE PRECISION DUMACH , DVNORM
        INTEGER INIt , MXStep , MXHnil , NHNil , NSLast , NYH , IOWns ,   &
       &        ICF , IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH ,   &
       &        LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,        &
       &        MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST , NFE ,     &
       &        NJE , NQU
        INTEGER JPRe , JACflg , LOCwp , LOCiwp , LSAvx , KMP , MAXl ,     &
       &        MNEwt , NNI , NLI , NPS , NCFn , NCFl
        INTEGER i , i1 , i2 , iflag , imxer , kgo , lf0 , leniw , leniwk ,&
       &        lenrw , lenwm , lenwk , liwp , lwp , mord , mxhnl0 ,      &
       &        mxstp0 , ncfn0 , ncfl0 , nli0 , nni0 , nnid , nstd , nwarn
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION DELt , EPCon , SQRtn , RSQrtn
        DOUBLE PRECISION atoli , avdim , ayi , big , ewti , h0 , hmax ,   &
       &                 hmx , rcfl , rcfn , rh , rtoli , tcrit , tdist , &
       &                 tnext , tol , tolsf , tp , size , sum , w0
        DIMENSION mord(2)
        LOGICAL ihit , lavd , lcfn , lcfl , lwarn
        CHARACTER*60 msg
        SAVE mord , mxstp0 , mxhnl0
  !-----------------------------------------------------------------------
  ! The following two internal Common blocks contain
  ! (a) variables which are local to any subroutine but whose values must
  !     be preserved between calls to the routine ("own" variables), and
  ! (b) variables which are communicated between subroutines.
  ! The block DLS001 is declared in subroutines DLSODPK, DINTDY, DSTODPK,
  ! DSOLPK, and DATV.
  ! The block DLPK01 is declared in subroutines DLSODPK, DSTODPK, DPKSET,
  ! and DSOLPK.
  ! Groups of variables are replaced by dummy arrays in the Common
  ! declarations in routines where those variables are not used.
  !-----------------------------------------------------------------------
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , INIt , MXStep , MXHnil ,       &
       &                NHNil , NSLast , NYH , IOWns(6) , ICF , IERpj ,   &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
        COMMON /DLPK01/ DELt , EPCon , SQRtn , RSQrtn , JPRe , JACflg ,   &
       &                LOCwp , LOCiwp , LSAvx , KMP , MAXl , MNEwt ,     &
       &                NNI , NLI , NPS , NCFn , NCFl
  !
        DATA mord(1) , mord(2)/12 , 5/ , mxstp0/500/ , mxhnl0/10/
  !-----------------------------------------------------------------------
  ! Block A.
  ! This code block is executed on every call.
  ! It tests ISTATE and ITASK for legality and branches appropriately.
  ! If ISTATE .gt. 1 but the flag INIT shows that initialization has
  ! not yet been done, an error return occurs.
  ! If ISTATE = 1 and TOUT = T, return immediately.
  !-----------------------------------------------------------------------
        IF ( Istate.LT.1 .OR. Istate.GT.3 ) THEN
  !-----------------------------------------------------------------------
  ! Block I.
  ! The following block handles all error returns due to illegal input
  ! (ISTATE = -3), as detected before calling the core integrator.
  ! First the error message routine is called.  If the illegal input
  ! is a negative ISTATE, the run is aborted (apparent infinite loop).
  !-----------------------------------------------------------------------
           msg = 'DLSODPK-  ISTATE(=I1) illegal.'
           CALL XERRWD(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
           IF ( Istate.GE.0 ) GOTO 2800
  !
           msg = 'DLSODPK-  Run aborted.. apparent infinite loop.   '
           CALL XERRWD(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
           GOTO 99999
        ELSE
           IF ( Itask.LT.1 .OR. Itask.GT.5 ) THEN
              msg = 'DLSODPK-  ITASK (=I1) illegal.'
              CALL XERRWD(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
              GOTO 2800
           ELSE
              IF ( Istate.EQ.1 ) THEN
                 INIt = 0
                 IF ( Tout.EQ.T ) RETURN
              ELSEIF ( INIt.EQ.0 ) THEN
                 msg =                                                    &
       &              'DLSODPK-  ISTATE.gt.1 but DLSODPK not initialized.'
                 CALL XERRWD(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
                 GOTO 2800
              ELSEIF ( Istate.EQ.2 ) THEN
                 GOTO 50
              ENDIF
  !-----------------------------------------------------------------------
  ! Block B.
  ! The next code block is executed for the initial call (ISTATE = 1),
  ! or for a continuation call with parameter changes (ISTATE = 3).
  ! It contains checking of all inputs and various initializations.
  !
  ! First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF.
  !-----------------------------------------------------------------------
              IF ( Neq(1).LE.0 ) THEN
                 msg = 'DLSODPK-  NEQ (=I1) .lt. 1    '
                 CALL XERRWD(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
                 GOTO 2800
              ELSE
                 IF ( Istate.NE.1 ) THEN
                    IF ( Neq(1).GT.N ) THEN
                       msg =                                              &
       &              'DLSODPK-  ISTATE = 3 and NEQ increased (I1 to I2).'
                       CALL XERRWD(msg,50,5,0,2,N,Neq(1),0,0.0D0,0.0D0)
                       GOTO 2800
                    ENDIF
                 ENDIF
                 N = Neq(1)
                 IF ( Itol.LT.1 .OR. Itol.GT.4 ) THEN
                    msg = 'DLSODPK-  ITOL (=I1) illegal. '
                    CALL XERRWD(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
                    GOTO 2800
                 ELSEIF ( Iopt.LT.0 .OR. Iopt.GT.1 ) THEN
                    msg = 'DLSODPK-  IOPT (=I1) illegal. '
                    CALL XERRWD(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
                    GOTO 2800
                 ELSE
                    METh = Mf/10
                    MITer = Mf - 10*METh
                    IF ( METh.LT.1 .OR. METh.GT.2 ) GOTO 2100
                    IF ( MITer.LT.0 ) GOTO 2100
                    IF ( MITer.GT.4 .AND. MITer.LT.9 ) GOTO 2100
                    IF ( MITer.GE.1 ) JPRe = Iwork(3)
                    JACflg = 0
                    IF ( MITer.GE.1 ) JACflg = Iwork(4)
  ! Next process and check the optional inputs. --------------------------
                    IF ( Iopt.EQ.1 ) THEN
                       MAXord = Iwork(5)
                       IF ( MAXord.LT.0 ) THEN
                          msg = 'DLSODPK-  MAXORD (=I1) .lt. 0 '
                          CALL XERRWD(msg,30,11,0,1,MAXord,0,0,0.0D0,     &
       &                              0.0D0)
                          GOTO 2800
                       ELSE
                          IF ( MAXord.EQ.0 ) MAXord = 100
                          MAXord = MIN(MAXord,mord(METh))
                          MXStep = Iwork(6)
                          IF ( MXStep.LT.0 ) THEN
                             msg = 'DLSODPK-  MXSTEP (=I1) .lt. 0 '
                             CALL XERRWD(msg,30,12,0,1,MXStep,0,0,0.0D0,  &
       &                                 0.0D0)
                             GOTO 2800
                          ELSE
                             IF ( MXStep.EQ.0 ) MXStep = mxstp0
                             MXHnil = Iwork(7)
                             IF ( MXHnil.LT.0 ) THEN
                                msg = 'DLSODPK-  MXHNIL (=I1) .lt. 0 '
                                CALL XERRWD(msg,30,13,0,1,MXHnil,0,0,     &
       &                           0.0D0,0.0D0)
                                GOTO 2800
                             ELSE
                                IF ( MXHnil.EQ.0 ) MXHnil = mxhnl0
                                IF ( Istate.EQ.1 ) THEN
                                   h0 = Rwork(5)
                                   IF ( (Tout-T)*h0.LT.0.0D0 ) THEN
                                      msg =                               &
       &                        'DLSODPK-  TOUT (=R1) behind T (=R2)     '
                                      CALL XERRWD(msg,40,14,0,0,0,0,2,    &
       &                                 Tout,T)
                                      msg =                               &
       &              '      Integration direction is given by H0 (=R1)  '
                                      CALL XERRWD(msg,50,14,0,0,0,0,1,h0, &
       &                                 0.0D0)
                                      GOTO 2800
                                   ENDIF
                                ENDIF
                                hmax = Rwork(6)
                                IF ( hmax.LT.0.0D0 ) THEN
                                   msg = 'DLSODPK-  HMAX (=R1) .lt. 0.0 '
                                   CALL XERRWD(msg,30,15,0,0,0,0,1,hmax,  &
       &                              0.0D0)
                                   GOTO 2800
                                ELSE
                                   HMXi = 0.0D0
                                   IF ( hmax.GT.0.0D0 ) HMXi = 1.0D0/hmax
                                   HMIn = Rwork(7)
                                   IF ( HMIn.LT.0.0D0 ) THEN
                                      msg =                               &
       &                                 'DLSODPK-  HMIN (=R1) .lt. 0.0 '
                                      CALL XERRWD(msg,30,16,0,0,0,0,1,    &
       &                                 HMIn,0.0D0)
                                      GOTO 2800
                                   ELSE
                                      MAXl = Iwork(8)
                                      IF ( MAXl.EQ.0 ) MAXl = 5
                                      MAXl = MIN(MAXl,N)
                                      KMP = Iwork(9)
                                      IF ( KMP.EQ.0 .OR. KMP.GT.MAXl )    &
       &                                 KMP = MAXl
                                      DELt = Rwork(8)
                                      IF ( DELt.EQ.0.0D0 ) DELt = 0.05D0
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ELSE
                       MAXord = mord(METh)
                       MXStep = mxstp0
                       MXHnil = mxhnl0
                       IF ( Istate.EQ.1 ) h0 = 0.0D0
                       HMXi = 0.0D0
                       HMIn = 0.0D0
                       MAXl = MIN(5,N)
                       KMP = MAXl
                       DELt = 0.05D0
                    ENDIF
  !-----------------------------------------------------------------------
  ! Set work array pointers and check lengths LRW and LIW.
  ! Pointers to segments of RWORK and IWORK are named by prefixing L to
  ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
  ! RWORK segments (in order) are denoted  YH, WM, EWT, SAVF, SAVX, ACOR.
  !-----------------------------------------------------------------------
                    LYH = 21
                    IF ( Istate.EQ.1 ) NYH = N
                    LWM = LYH + (MAXord+1)*NYH
                    IF ( MITer.EQ.0 ) lenwk = 0
                    IF ( MITer.EQ.1 ) lenwk = N*(MAXl+2) + MAXl*MAXl
                    IF ( MITer.EQ.2 ) lenwk = N*(MAXl+2+MIN(1,MAXl-KMP))  &
       &                 + (MAXl+3)*MAXl + 1
                    IF ( MITer.EQ.3 .OR. MITer.EQ.4 ) lenwk = 5*N
                    IF ( MITer.EQ.9 ) lenwk = 2*N
                    lwp = 0
                    IF ( MITer.GE.1 ) lwp = Iwork(1)
                    lenwm = lenwk + lwp
                    LOCwp = lenwk + 1
                    LEWt = LWM + lenwm
                    LSAvf = LEWt + N
                    LSAvx = LSAvf + N
                    LACor = LSAvx + N
                    IF ( MITer.EQ.0 ) LACor = LSAvf + N
                    lenrw = LACor + N - 1
                    Iwork(17) = lenrw
                    LIWm = 31
                    leniwk = 0
                    IF ( MITer.EQ.1 ) leniwk = MAXl
                    liwp = 0
                    IF ( MITer.GE.1 ) liwp = Iwork(2)
                    leniw = 30 + leniwk + liwp
                    LOCiwp = leniwk + 1
                    Iwork(18) = leniw
                    IF ( lenrw.GT.Lrw ) THEN
                       msg =                                              &
       &    'DLSODPK-  RWORK length needed, LENRW(=I1), exceeds LRW(=I2) '
                       CALL XERRWD(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                       GOTO 2800
                    ELSEIF ( leniw.GT.Liw ) THEN
                       msg =                                              &
       &    'DLSODPK-  IWORK length needed, LENIW(=I1), exceeds LIW(=I2) '
                       CALL XERRWD(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                       GOTO 2800
                    ELSE
  ! Check RTOL and ATOL for legality. ------------------------------------
                       rtoli = Rtol(1)
                       atoli = Atol(1)
                       DO i = 1 , N
                          IF ( Itol.GE.3 ) rtoli = Rtol(i)
                          IF ( Itol.EQ.2 .OR. Itol.EQ.4 ) atoli = Atol(i)
                          IF ( rtoli.LT.0.0D0 ) GOTO 2200
                          IF ( atoli.LT.0.0D0 ) GOTO 2300
                       ENDDO
  ! Load SQRT(N) and its reciprocal in Common. ---------------------------
                       SQRtn = SQRT(REAL(N))
                       RSQrtn = 1.0D0/SQRtn
                       IF ( Istate.EQ.1 ) THEN
  !-----------------------------------------------------------------------
  ! Block C.
  ! The next block is for the initial call only (ISTATE = 1).
  ! It contains all remaining initializations, the initial call to F,
  ! and the calculation of the initial step size.
  ! The error weights in EWT are inverted after being loaded.
  !-----------------------------------------------------------------------
                          UROund = DUMACH()
                          TN = T
                          IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
                             tcrit = Rwork(1)
                             IF ( (tcrit-Tout)*(Tout-T).LT.0.0D0 )        &
       &                          GOTO 2600
                             IF ( h0.NE.0.0D0 .AND. (T+h0-tcrit)          &
       &                          *h0.GT.0.0D0 ) h0 = tcrit - T
                          ENDIF
                          JSTart = 0
                          NHNil = 0
                          NST = 0
                          NJE = 0
                          NSLast = 0
                          nli0 = 0
                          nni0 = 0
                          ncfn0 = 0
                          ncfl0 = 0
                          nwarn = 0
                          HU = 0.0D0
                          NQU = 0
                          CCMax = 0.3D0
                          MAXcor = 3
                          MSBp = 20
                          MXNcf = 10
                          NNI = 0
                          NLI = 0
                          NPS = 0
                          NCFn = 0
                          NCFl = 0
  ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
                          lf0 = LYH + NYH
                          CALL F(Neq,T,Y,Rwork(lf0))
                          NFE = 1
  ! Load the initial value vector in YH. ---------------------------------
                          DO i = 1 , N
                             Rwork(i+LYH-1) = Y(i)
                          ENDDO
  ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                          NQ = 1
                          H = 1.0D0
                          CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),        &
       &                              Rwork(LEWt))
                          DO i = 1 , N
                             IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 2400
                             Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
                          ENDDO
  !-----------------------------------------------------------------------
  ! The coding below computes the step size, H0, to be attempted on the
  ! first step, unless the user has supplied a value for this.
  ! First check that TOUT - T differs significantly from zero.
  ! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
  ! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
  ! so as to be between 100*UROUND and 1.0E-3.
  ! Then the computed value H0 is given by..
  !                                      NEQ
  !   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
  !                                       1
  ! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
  !         f(i)   = i-th component of initial value of f,
  !         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
  ! The sign of H0 is inferred from the initial values of TOUT and T.
  !-----------------------------------------------------------------------
                          IF ( h0.EQ.0.0D0 ) THEN
                             tdist = ABS(Tout-T)
                             w0 = MAX(ABS(T),ABS(Tout))
                             IF ( tdist.LT.2.0D0*UROund*w0 ) THEN
                                msg =                                     &
       &    'DLSODPK- TOUT(=R1) too close to T(=R2) to start integration.'
                                CALL XERRWD(msg,60,22,0,0,0,0,2,Tout,T)
                                GOTO 2800
                             ELSE
                                tol = Rtol(1)
                                IF ( Itol.GT.2 ) THEN
                                   DO i = 1 , N
                                      tol = MAX(tol,Rtol(i))
                                   ENDDO
                                ENDIF
                                IF ( tol.LE.0.0D0 ) THEN
                                   atoli = Atol(1)
                                   DO i = 1 , N
                                      IF ( Itol.EQ.2 .OR. Itol.EQ.4 )     &
       &                                 atoli = Atol(i)
                                      ayi = ABS(Y(i))
                                      IF ( ayi.NE.0.0D0 )                 &
       &                                 tol = MAX(tol,atoli/ayi)
                                   ENDDO
                                ENDIF
                                tol = MAX(tol,100.0D0*UROund)
                                tol = MIN(tol,0.001D0)
                                sum = DVNORM(N,Rwork(lf0),Rwork(LEWt))
                                sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                                h0 = 1.0D0/SQRT(sum)
                                h0 = MIN(h0,tdist)
                                h0 = SIGN(h0,Tout-T)
                             ENDIF
                          ENDIF
  ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                          rh = ABS(h0)*HMXi
                          IF ( rh.GT.1.0D0 ) h0 = h0/rh
  ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                          H = h0
                          DO i = 1 , N
                             Rwork(i+lf0-1) = h0*Rwork(i+lf0-1)
                          ENDDO
                          GOTO 700
                       ELSE
  ! If ISTATE = 3, set flag to signal parameter changes to DSTODPK. ------
                          JSTart = -1
                          IF ( NQ.GT.MAXord ) THEN
  ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
                             DO i = 1 , N
                                Rwork(i+LSAvf-1) = Rwork(i+LWM-1)
                             ENDDO
                          ENDIF
                          IF ( N.NE.NYH ) THEN
  ! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                             i1 = LYH + L*NYH
                             i2 = LYH + (MAXord+1)*NYH - 1
                             IF ( i1.LE.i2 ) THEN
                                DO i = i1 , i2
                                   Rwork(i) = 0.0D0
                                ENDDO
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  ! Block D.
  ! The next code block is for continuation calls only (ISTATE = 2 or 3)
  ! and is to check stop conditions before taking a step.
  !-----------------------------------------------------------------------
   50      NSLast = NST
           nli0 = NLI
           nni0 = NNI
           ncfn0 = NCFn
           ncfl0 = NCFl
           nwarn = 0
           GOTO (100,600,200,300,400) , Itask
        ENDIF
   100  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2700
        T = Tout
        GOTO 1400
   200  tp = TN - HU*(1.0D0+100.0D0*UROund)
        IF ( (tp-Tout)*H.GT.0.0D0 ) THEN
           msg =                                                          &
       &    'DLSODPK-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2) '
           CALL XERRWD(msg,60,23,0,1,Itask,0,2,Tout,tp)
           GOTO 2800
        ELSE
           IF ( (TN-Tout)*H.GE.0.0D0 ) GOTO 1300
           GOTO 600
        ENDIF
   300  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2500
        IF ( (tcrit-Tout)*H.LT.0.0D0 ) GOTO 2600
        IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 500
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2700
        T = Tout
        GOTO 1400
   400  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2500
   500  hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
        IF ( ihit ) GOTO 1300
        tnext = TN + H*(1.0D0+4.0D0*UROund)
        IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
           H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
           IF ( Istate.EQ.2 ) JSTart = -2
        ENDIF
  !-----------------------------------------------------------------------
  ! Block E.
  ! The next block is normally executed for all calls and contains
  ! the call to the one-step core integrator DSTODPK.
  !
  ! This is a looping point for the integration steps.
  !
  ! First check for too many steps being taken,
  ! Check for poor Newton/Krylov method performance, update EWT (if not
  ! at start of problem), check for too much accuracy being requested,
  ! and check for H below the roundoff level in T.
  !-----------------------------------------------------------------------
   600  IF ( (NST-NSLast).GE.MXStep ) THEN
  !-----------------------------------------------------------------------
  ! Block H.
  ! The following block handles all unsuccessful returns other than
  ! those for illegal input.  First the error message routine is called.
  ! If there was an error test or convergence test failure, IMXER is set.
  ! Then Y is loaded from YH and T is set to TN.
  ! The optional outputs are loaded into the work arrays before returning.
  !-----------------------------------------------------------------------
  ! The maximum number of steps was taken before reaching TOUT. ----------
           msg = 'DLSODPK-  At current T (=R1), MXSTEP (=I1) steps  '
           CALL XERRWD(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      taken on this call before reaching TOUT     '
           CALL XERRWD(msg,50,201,0,1,MXStep,0,1,TN,0.0D0)
           Istate = -1
           GOTO 2000
        ELSE
           nstd = NST - NSLast
           nnid = NNI - nni0
           IF ( nstd.GE.10 .AND. nnid.NE.0 ) THEN
              avdim = REAL(NLI-nli0)/REAL(nnid)
              rcfn = REAL(NCFn-ncfn0)/REAL(nstd)
              rcfl = REAL(NCFl-ncfl0)/REAL(nnid)
              lavd = avdim.GT.(MAXl-0.05D0)
              lcfn = rcfn.GT.0.9D0
              lcfl = rcfl.GT.0.9D0
              lwarn = lavd .OR. lcfn .OR. lcfl
              IF ( lwarn ) THEN
                 nwarn = nwarn + 1
                 IF ( nwarn.LE.10 ) THEN
                    IF ( lavd ) THEN
                       msg =                                              &
       &    'DLSODPK- Warning. Poor iterative algorithm performance seen '
                       CALL XERRWD(msg,60,111,0,0,0,0,0,0.0D0,0.0D0)
                    ENDIF
                    IF ( lavd ) THEN
                       msg =                                              &
       &    '      at T = R1 by average no. of linear iterations = R2    '
                       CALL XERRWD(msg,60,111,0,0,0,0,2,TN,avdim)
                    ENDIF
                    IF ( lcfn ) THEN
                       msg =                                              &
       &    'DLSODPK- Warning. Poor iterative algorithm performance seen '
                       CALL XERRWD(msg,60,112,0,0,0,0,0,0.0D0,0.0D0)
                    ENDIF
                    IF ( lcfn ) THEN
                       msg =                                              &
       &    '      at T = R1 by nonlinear convergence failure rate = R2  '
                       CALL XERRWD(msg,60,112,0,0,0,0,2,TN,rcfn)
                    ENDIF
                    IF ( lcfl ) THEN
                       msg =                                              &
       &    'DLSODPK- Warning. Poor iterative algorithm performance seen '
                       CALL XERRWD(msg,60,113,0,0,0,0,0,0.0D0,0.0D0)
                    ENDIF
                    IF ( lcfl ) THEN
                       msg =                                              &
       &    '      at T = R1 by linear convergence failure rate = R2     '
                       CALL XERRWD(msg,60,113,0,0,0,0,2,TN,rcfl)
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
           CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
           DO i = 1 , N
              IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 1500
              Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
           ENDDO
        ENDIF
   700  tolsf = UROund*DVNORM(N,Rwork(LYH),Rwork(LEWt))
        IF ( tolsf.LE.1.0D0 ) THEN
           IF ( (TN+H).EQ.TN ) THEN
              NHNil = NHNil + 1
              IF ( NHNil.LE.MXHnil ) THEN
                 msg =                                                    &
       &              'DLSODPK-  Warning..Internal T(=R1) and H(=R2) are '
                 CALL XERRWD(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &    '      such that in the machine, T + H = T on the next step  '
                 CALL XERRWD(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &              '     (H = step size). Solver will continue anyway.'
                 CALL XERRWD(msg,50,101,0,0,0,0,2,TN,H)
                 IF ( NHNil.GE.MXHnil ) THEN
                    msg =                                                 &
       &              'DLSODPK-  Above warning has been issued I1 times. '
                    CALL XERRWD(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
                    msg =                                                 &
       &              '     It will not be issued again for this problem.'
                    CALL XERRWD(msg,50,102,0,1,MXHnil,0,0,0.0D0,0.0D0)
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  !     CALL DSTODPK(NEQ,Y,YH,NYH,YH,EWT,SAVF,SAVX,ACOR,WM,IWM,F,JAC,PSOL)
  !-----------------------------------------------------------------------
           CALL DSTODPK(Neq,Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),      &
       &                Rwork(LSAvf),Rwork(LSAvx),Rwork(LACor),Rwork(LWM),&
       &                Iwork(LIWm),F,JAC,PSOL)
           kgo = 1 - KFLag
           GOTO (800,1600,1700,1800) , kgo
        ELSE
           tolsf = tolsf*2.0D0
           IF ( NST.EQ.0 ) THEN
              msg = 'DLSODPK-  At start of problem, too much accuracy  '
              CALL XERRWD(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
              msg =                                                       &
       &    '      requested for precision of machine..  See TOLSF (=R1) '
              CALL XERRWD(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
              Rwork(14) = tolsf
              GOTO 2800
           ELSE
  ! Too much accuracy requested for machine precision. -------------------
              msg = 'DLSODPK-  At T (=R1), too much accuracy requested '
              CALL XERRWD(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
              msg = '      for precision of machine..  See TOLSF (=R2) '
              CALL XERRWD(msg,50,203,0,0,0,0,2,TN,tolsf)
              Rwork(14) = tolsf
              Istate = -2
              GOTO 2000
           ENDIF
        ENDIF
  !-----------------------------------------------------------------------
  ! Block F.
  ! The following block handles the case of a successful return from the
  ! core integrator (KFLAG = 0).  Test for stop conditions.
  !-----------------------------------------------------------------------
   800  INIt = 1
        GOTO (900,1300,1000,1100,1200) , Itask
  ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   900  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        T = Tout
        GOTO 1400
  ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   1000 IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        GOTO 1300
  ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
   1100 IF ( (TN-Tout)*H.LT.0.0D0 ) THEN
           hmx = ABS(TN) + ABS(H)
           ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
           IF ( ihit ) GOTO 1300
           tnext = TN + H*(1.0D0+4.0D0*UROund)
           IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
              H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
              JSTart = -2
           ENDIF
           GOTO 600
        ELSE
           CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
           T = Tout
           GOTO 1400
        ENDIF
  ! ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
   1200 hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
  !-----------------------------------------------------------------------
  ! Block G.
  ! The following block handles all successful returns from DLSODPK.
  ! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
  ! ISTATE is set to 2, and the optional outputs are loaded into the
  ! work arrays before returning.
  !-----------------------------------------------------------------------
   1300 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
           IF ( ihit ) T = tcrit
        ENDIF
   1400 Istate = 2
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        Iwork(19) = NNI
        Iwork(20) = NLI
        Iwork(21) = NPS
        Iwork(22) = NCFn
        Iwork(23) = NCFl
        RETURN
  ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
   1500 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODPK-  At T (=R1), EWT(I1) has become R2.le.0. '
        CALL XERRWD(msg,50,202,0,1,i,0,2,TN,ewti)
        Istate = -6
        GOTO 2000
  ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   1600 msg = 'DLSODPK-  At T(=R1), step size H(=R2), the error  '
        CALL XERRWD(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      test failed repeatedly or with ABS(H) = HMIN'
        CALL XERRWD(msg,50,204,0,0,0,0,2,TN,H)
        Istate = -4
        GOTO 1900
  ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
   1700 msg = 'DLSODPK-  At T (=R1) and step size H (=R2), the   '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      corrector convergence failed repeatedly     '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      or with ABS(H) = HMIN   '
        CALL XERRWD(msg,30,205,0,0,0,0,2,TN,H)
        Istate = -5
        GOTO 1900
  ! KFLAG = -3.  Unrecoverable error from PSOL. --------------------------
   1800 msg = 'DLSODPK-  At T (=R1) an unrecoverable error return'
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      was made from Subroutine PSOL     '
        CALL XERRWD(msg,40,205,0,0,0,0,1,TN,0.0D0)
        Istate = -7
        GOTO 2000
  ! Compute IMXER if relevant. -------------------------------------------
   1900 big = 0.0D0
        imxer = 1
        DO i = 1 , N
           size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
           IF ( big.LT.size ) THEN
              big = size
              imxer = i
           ENDIF
        ENDDO
        Iwork(16) = imxer
  ! Set Y vector, T, and optional outputs. -------------------------------
   2000 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        Iwork(19) = NNI
        Iwork(20) = NLI
        Iwork(21) = NPS
        Iwork(22) = NCFn
        Iwork(23) = NCFl
        RETURN
   2100 msg = 'DLSODPK-  MF (=I1) illegal.   '
        CALL XERRWD(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
        GOTO 2800
   2200 msg = 'DLSODPK-  RTOL(I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
        GOTO 2800
   2300 msg = 'DLSODPK-  ATOL(I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,20,0,1,i,0,1,atoli,0.0D0)
        GOTO 2800
   2400 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODPK-  EWT(I1) is R1 .le. 0.0        '
        CALL XERRWD(msg,40,21,0,1,i,0,1,ewti,0.0D0)
        GOTO 2800
   2500 msg =                                                             &
       &    'DLSODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)  '
        CALL XERRWD(msg,60,24,0,0,0,0,2,tcrit,TN)
        GOTO 2800
   2600 msg =                                                             &
       &    'DLSODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)  '
        CALL XERRWD(msg,60,25,0,0,0,0,2,tcrit,Tout)
        GOTO 2800
   2700 msg = 'DLSODPK-  Trouble in DINTDY. ITASK = I1, TOUT = R1'
        CALL XERRWD(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
  !
   2800 Istate = -3
        RETURN
  !----------------------- End of Subroutine DLSODPK ---------------------
  99999 END
  !*==DLSODI.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
        SUBROUTINE DLSODI(RES,ADDA,JAC,Neq,Y,Ydoti,T,Tout,Itol,Rtol,Atol, &
       &                  Itask,Istate,Iopt,Rwork,Lrw,Iwork,Liw,Mf)
        IMPLICIT NONE
  !*--DLSODI5459
  !*** Start of declarations inserted by SPAG
        !REAL ADDA , DPREPJI , DSOLSY
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL RES , ADDA , JAC
        INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
       &        Mf
        DOUBLE PRECISION Y , Ydoti , T , Tout , Rtol , Atol , Rwork
        DIMENSION Neq(*) , Y(*) , Ydoti(*) , Rtol(*) , Atol(*) ,          &
       &          Rwork(Lrw) , Iwork(Liw)
  !-----------------------------------------------------------------------
  ! This is the 18 November 2003 version of
  ! DLSODI: Livermore Solver for Ordinary Differential Equations
  !         (Implicit form).
  !
  ! This version is in double precision.
  !
  ! DLSODI solves the initial value problem for linearly implicit
  ! systems of first order ODEs,
  !     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
  ! or, in component form,
  !     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
  !        i,1      1                     i,NEQ      NEQ
  !
  !      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
  !           i      1   2       NEQ
  !
  ! If A is singular, this is a differential-algebraic system.
  !
  ! DLSODI is a variant version of the DLSODE package.
  !-----------------------------------------------------------------------
  ! Reference:
  !     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
  !     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  !     North-Holland, Amsterdam, 1983, pp. 55-64.
  !-----------------------------------------------------------------------
  ! Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
  !                Center for Applied Scientific Computing, L-561
  !                Lawrence Livermore National Laboratory
  !                Livermore, CA 94551
  !-----------------------------------------------------------------------
  ! Summary of Usage.
  !
  ! Communication between the user and the DLSODI package, for normal
  ! situations, is summarized here.  This summary describes only a subset
  ! of the full set of options available.  See the full description for
  ! details, including optional communication, nonstandard options,
  ! and instructions for special situations.  See also the example
  ! problem (with program and output) following this summary.
  !
  ! A. First, provide a subroutine of the form:
  !               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
  !               DOUBLE PRECISION T, Y(*), S(*), R(*)
  ! which computes the residual function
  !     r = g(t,y)  -  A(t,y) * s ,
  ! as a function of t and the vectors y and s.  (s is an internally
  ! generated approximation to dy/dt.)  The arrays Y and S are inputs
  ! to the RES routine and should not be altered.  The residual
  ! vector is to be stored in the array R.  The argument IRES should be
  ! ignored for casual use of DLSODI.  (For uses of IRES, see the
  ! paragraph on RES in the full description below.)
  !
  ! B. Next, decide whether full or banded form is more economical
  ! for the storage of matrices.  DLSODI must deal internally with the
  ! matrices A and dr/dy, where r is the residual function defined above.
  ! DLSODI generates a linear combination of these two matrices, and
  ! this is treated in either full or banded form.
  !     The matrix structure is communicated by a method flag MF,
  ! which is 21 or 22 for the full case, and 24 or 25 in the band case.
  !     In the banded case, DLSODI requires two half-bandwidth
  ! parameters ML and MU.  These are, respectively, the widths of the
  ! lower and upper parts of the band, excluding the main diagonal.
  ! Thus the band consists of the locations (i,j) with
  ! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
  ! Note that the band must accommodate the nonzero elements of
  ! A(t,y), dg/dy, and d(A*s)/dy (s fixed).  Alternatively, one
  ! can define a band that encloses only the elements that are relatively
  ! large in magnitude, and gain some economy in storage and possibly
  ! also efficiency, although the appropriate threshhold for
  ! retaining matrix elements is highly problem-dependent.
  !
  ! C. You must also provide a subroutine of the form:
  !               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
  !               DOUBLE PRECISION T, Y(*), P(NROWP,*)
  ! which adds the matrix A = A(t,y) to the contents of the array P.
  ! T and the Y array are input and should not be altered.
  !     In the full matrix case, this routine should add elements of
  ! to P in the usual order.  I.e., add A(i,j) to P(i,j).  (Ignore the
  ! ML and MU arguments in this case.)
  !     In the band matrix case, this routine should add element A(i,j)
  ! to P(i-j+MU+1,j).  I.e., add the diagonal lines of A to the rows of
  ! P from the top down (the top line of A added to the first row of P).
  !
  ! D. For the sake of efficiency, you are encouraged to supply the
  ! Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
  ! (s = a fixed vector) as above.  If dr/dy is being supplied,
  ! use MF = 21 or 24, and provide a subroutine of the form:
  !               SUBROUTINE JAC (NEQ, T, Y, S, ML, MU, P, NROWP)
  !               DOUBLE PRECISION T, Y(*), S(*), P(NROWP,*)
  ! which computes dr/dy as a function of t, y, and s.  Here T, Y, and
  ! S are inputs, and the routine is to load dr/dy into P as follows:
  !     In the full matrix case (MF = 21), load P(i,j) with dr(i)/dy(j),
  ! the partial derivative of r(i) with respect to y(j).  (Ignore the
  ! ML and MU arguments in this case.)
  !     In the band matrix case (MF = 24), load P(i-j+mu+1,j) with
  ! dr(i)/dy(j), i.e. load the diagonal lines of dr/dy into the rows of
  ! P from the top down.
  !     In either case, only nonzero elements need be loaded, and the
  ! indexing of P is the same as in the ADDA routine.
  !     Note that if A is independent of y (or this dependence
  ! is weak enough to be ignored) then JAC is to compute dg/dy.
  !     If it is not feasible to provide a JAC routine, use
  ! MF = 22 or 25, and DLSODI will compute an approximate Jacobian
  ! internally by difference quotients.
  !
  ! E. Next decide whether or not to provide the initial value of the
  ! derivative vector dy/dt.  If the initial value of A(t,y) is
  ! nonsingular (and not too ill-conditioned), you may let DLSODI compute
  ! this vector (ISTATE = 0).  (DLSODI will solve the system A*s = g for
  ! s, with initial values of A and g.)  If A(t,y) is initially
  ! singular, then the system is a differential-algebraic system, and
  ! you must make use of the particular form of the system to compute the
  ! initial values of y and dy/dt.  In that case, use ISTATE = 1 and
  ! load the initial value of dy/dt into the array YDOTI.
  ! The input array YDOTI and the initial Y array must be consistent with
  ! the equations A*dy/dt = g.  This implies that the initial residual
  ! r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
  !
  ! F. Write a main program which calls Subroutine DLSODI once for
  ! each point at which answers are desired.  This should also provide
  ! for possible use of logical unit 6 for output of error messages
  ! by DLSODI.  On the first call to DLSODI, supply arguments as follows:
  ! RES    = name of user subroutine for residual function r.
  ! ADDA   = name of user subroutine for computing and adding A(t,y).
  ! JAC    = name of user subroutine for Jacobian matrix dr/dy
  !          (MF = 21 or 24).  If not used, pass a dummy name.
  ! Note: the names for the RES and ADDA routines and (if used) the
  !        JAC routine must be declared External in the calling program.
  ! NEQ    = number of scalar equations in the system.
  ! Y      = array of initial values, of length NEQ.
  ! YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
  ! T      = the initial value of the independent variable.
  ! TOUT   = first point where output is desired (.ne. T).
  ! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  ! RTOL   = relative tolerance parameter (scalar).
  ! ATOL   = absolute tolerance parameter (scalar or array).
  !          the estimated local error in y(i) will be controlled so as
  !          to be roughly less (in magnitude) than
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
  !          Thus the local error test passes if, in each component,
  !          either the absolute error is less than ATOL (or ATOL(i)),
  !          or the relative error is less than RTOL.
  !          Use RTOL = 0.0 for pure absolute error control, and
  !          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
  !          control.  Caution: Actual (global) errors may exceed these
  !          local tolerances, so choose them conservatively.
  ! ITASK  = 1 for normal computation of output values of y at t = TOUT.
  ! ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
  !          initial dy/dt is supplied, and 0 otherwise.
  ! IOPT   = 0 to indicate no optional inputs used.
  ! RWORK  = real work array of length at least:
  !             22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
  !             22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
  ! LRW    = declared length of RWORK (in user's dimension).
  ! IWORK  = integer work array of length at least 20 + NEQ.
  !          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
  !          and upper half-bandwidths ML,MU.
  ! LIW    = declared length of IWORK (in user's dimension).
  ! MF     = method flag.  Standard values are:
  !          21 for a user-supplied full Jacobian.
  !          22 for an internally generated full Jacobian.
  !          24 for a user-supplied banded Jacobian.
  !          25 for an internally generated banded Jacobian.
  !          for other choices of MF, see the paragraph on MF in
  !          the full description below.
  ! Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
  ! and possibly ATOL.
  !
  ! G. The output from the first call (or any call) is:
  !      Y = array of computed values of y(t) vector.
  !      T = corresponding value of independent variable (normally TOUT).
  ! ISTATE = 2  if DLSODI was successful, negative otherwise.
  !          -1 means excess work done on this call (check all inputs).
  !          -2 means excess accuracy requested (tolerances too small).
  !          -3 means illegal input detected (see printed message).
  !          -4 means repeated error test failures (check all inputs).
  !          -5 means repeated convergence failures (perhaps bad Jacobian
  !             supplied or wrong choice of tolerances).
  !          -6 means error weight became zero during problem. (Solution
  !             component i vanished, and ATOL or ATOL(i) = 0.)
  !          -7 cannot occur in casual use.
  !          -8 means DLSODI was unable to compute the initial dy/dt.
  !             In casual use, this means A(t,y) is initially singular.
  !             Supply YDOTI and use ISTATE = 1 on the first call.
  !
  !  If DLSODI returns ISTATE = -1, -4, or -5, then the output of
  !  DLSODI also includes YDOTI = array containing residual vector
  !  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
  !
  ! H. To continue the integration after a successful return, simply
  ! reset TOUT and call DLSODI again.  No other parameters need be reset.
  !
  !-----------------------------------------------------------------------
  ! Example Problem.
  !
  ! The following is a simple example problem, with the coding
  ! needed for its solution by DLSODI.  The problem is from chemical
  ! kinetics, and consists of the following three equations:
  !     dy1/dt = -.04*y1 + 1.e4*y2*y3
  !     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
  !       0.   = y1 + y2 + y3 - 1.
  ! on the interval from t = 0.0 to t = 4.e10, with initial conditions
  ! y1 = 1.0, y2 = y3 = 0.
  !
  ! The following coding solves this problem with DLSODI, using MF = 21
  ! and printing results at t = .4, 4., ..., 4.e10.  It uses
  ! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
  ! y2 has much smaller values.  dy/dt is supplied in YDOTI. We had
  ! obtained the initial value of dy3/dt by differentiating the
  ! third equation and evaluating the first two at t = 0.
  ! At the end of the run, statistical quantities of interest are
  ! printed (see optional outputs in the full description below).
  !
  !     EXTERNAL RESID, APLUSP, DGBYDY
  !     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
  !     DIMENSION Y(3), YDOTI(3), ATOL(3), RWORK(58), IWORK(23)
  !     NEQ = 3
  !     Y(1) = 1.
  !     Y(2) = 0.
  !     Y(3) = 0.
  !     YDOTI(1) = -.04
  !     YDOTI(2) =  .04
  !     YDOTI(3) =  0.
  !     T = 0.
  !     TOUT = .4
  !     ITOL = 2
  !     RTOL = 1.D-4
  !     ATOL(1) = 1.D-6
  !     ATOL(2) = 1.D-10
  !     ATOL(3) = 1.D-6
  !     ITASK = 1
  !     ISTATE = 1
  !     IOPT = 0
  !     LRW = 58
  !     LIW = 23
  !     MF = 21
  !     DO 40  IOUT = 1,12
  !       CALL DLSODI(RESID, APLUSP, DGBYDY, NEQ, Y, YDOTI, T, TOUT, ITOL,
  !    1     RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, MF)
  !       WRITE (6,20)  T, Y(1), Y(2), Y(3)
  !  20   FORMAT(' At t =',D12.4,'   Y =',3D14.6)
  !       IF (ISTATE .LT. 0 )  GO TO 80
  !  40   TOUT = TOUT*10.
  !     WRITE (6,60)  IWORK(11), IWORK(12), IWORK(13)
  !  60 FORMAT(/' No. steps =',I4,'  No. r-s =',I4,'  No. J-s =',I4)
  !     STOP
  !  80 WRITE (6,90)  ISTATE
  !  90 FORMAT(///' Error halt.. ISTATE =',I3)
  !     STOP
  !     END
  !
  !     SUBROUTINE RESID(NEQ, T, Y, S, R, IRES)
  !     DOUBLE PRECISION T, Y, S, R
  !     DIMENSION Y(3), S(3), R(3)
  !     R(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3) - S(1)
  !     R(2) = .04*Y(1) - 1.D4*Y(2)*Y(3) - 3.D7*Y(2)*Y(2) - S(2)
  !     R(3) = Y(1) + Y(2) + Y(3) - 1.
  !     RETURN
  !     END
  !
  !     SUBROUTINE APLUSP(NEQ, T, Y, ML, MU, P, NROWP)
  !     DOUBLE PRECISION T, Y, P
  !     DIMENSION Y(3), P(NROWP,3)
  !     P(1,1) = P(1,1) + 1.
  !     P(2,2) = P(2,2) + 1.
  !     RETURN
  !     END
  !
  !     SUBROUTINE DGBYDY(NEQ, T, Y, S, ML, MU, P, NROWP)
  !     DOUBLE PRECISION T, Y, S, P
  !     DIMENSION Y(3), S(3), P(NROWP,3)
  !     P(1,1) = -.04
  !     P(1,2) = 1.D4*Y(3)
  !     P(1,3) = 1.D4*Y(2)
  !     P(2,1) = .04
  !     P(2,2) = -1.D4*Y(3) - 6.D7*Y(2)
  !     P(2,3) = -1.D4*Y(2)
  !     P(3,1) = 1.
  !     P(3,2) = 1.
  !     P(3,3) = 1.
  !     RETURN
  !     END
  !
  ! The output of this program (on a CDC-7600 in single precision)
  ! is as follows:
  !
  !   At t =  4.0000e-01   Y =  9.851726e-01  3.386406e-05  1.479357e-02
  !   At t =  4.0000e+00   Y =  9.055142e-01  2.240418e-05  9.446344e-02
  !   At t =  4.0000e+01   Y =  7.158050e-01  9.184616e-06  2.841858e-01
  !   At t =  4.0000e+02   Y =  4.504846e-01  3.222434e-06  5.495122e-01
  !   At t =  4.0000e+03   Y =  1.831701e-01  8.940379e-07  8.168290e-01
  !   At t =  4.0000e+04   Y =  3.897016e-02  1.621193e-07  9.610297e-01
  !   At t =  4.0000e+05   Y =  4.935213e-03  1.983756e-08  9.950648e-01
  !   At t =  4.0000e+06   Y =  5.159269e-04  2.064759e-09  9.994841e-01
  !   At t =  4.0000e+07   Y =  5.306413e-05  2.122677e-10  9.999469e-01
  !   At t =  4.0000e+08   Y =  5.494532e-06  2.197826e-11  9.999945e-01
  !   At t =  4.0000e+09   Y =  5.129457e-07  2.051784e-12  9.999995e-01
  !   At t =  4.0000e+10   Y = -7.170472e-08 -2.868188e-13  1.000000e+00
  !
  !   No. steps = 330  No. r-s = 404  No. J-s =  69
  !
  !-----------------------------------------------------------------------
  ! Full Description of User Interface to DLSODI.
  !
  ! The user interface to DLSODI consists of the following parts.
  !
  ! 1.   The call sequence to Subroutine DLSODI, which is a driver
  !      routine for the solver.  This includes descriptions of both
  !      the call sequence arguments and of user-supplied routines.
  !      Following these descriptions is a description of
  !      optional inputs available through the call sequence, and then
  !      a description of optional outputs (in the work arrays).
  !
  ! 2.   Descriptions of other routines in the DLSODI package that may be
  !      (optionally) called by the user.  These provide the ability to
  !      alter error message handling, save and restore the internal
  !      Common, and obtain specified derivatives of the solution y(t).
  !
  ! 3.   Descriptions of Common blocks to be declared in overlay
  !      or similar environments, or to be saved when doing an interrupt
  !      of the problem and continued solution later.
  !
  ! 4.   Description of two routines in the DLSODI package, either of
  !      which the user may replace with his/her own version, if desired.
  !      These relate to the measurement of errors.
  !
  !-----------------------------------------------------------------------
  ! Part 1.  Call Sequence.
  !
  ! The call sequence parameters used for input only are
  !     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
  !     IOPT, LRW, LIW, MF,
  ! and those used for both input and output are
  !     Y, T, ISTATE, YDOTI.
  ! The work arrays RWORK and IWORK are also used for conditional and
  ! optional inputs and optional outputs.  (The term output here refers
  ! to the return from Subroutine DLSODI to the user's calling program.)
  !
  ! The legality of input parameters will be thoroughly checked on the
  ! initial call for the problem, but not checked thereafter unless a
  ! change in input parameters is flagged by ISTATE = 3 on input.
  !
  ! The descriptions of the call arguments are as follows.
  !
  ! RES    = the name of the user-supplied subroutine which supplies
  !          the residual vector for the ODE system, defined by
  !            r = g(t,y) - A(t,y) * s
  !          as a function of the scalar t and the vectors
  !          s and y (s approximates dy/dt).  This subroutine
  !          is to have the form
  !               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
  !               DOUBLE PRECISION T, Y(*), S(*), R(*)
  !          where NEQ, T, Y, S, and IRES are input, and R and
  !          IRES are output.  Y, S, and R are arrays of length NEQ.
  !             On input, IRES indicates how DLSODI will use the
  !          returned array R, as follows:
  !             IRES = 1  means that DLSODI needs the full residual,
  !                       r = g - A*s, exactly.
  !             IRES = -1 means that DLSODI is using R only to compute
  !                       the Jacobian dr/dy by difference quotients.
  !          The RES routine can ignore IRES, or it can omit some terms
  !          if IRES = -1.  If A does not depend on y, then RES can
  !          just return R = g when IRES = -1.  If g - A*s contains other
  !          additive terms that are independent of y, these can also be
  !          dropped, if done consistently, when IRES = -1.
  !             The subroutine should set the flag IRES if it
  !          encounters a halt condition or illegal input.
  !          Otherwise, it should not reset IRES.  On output,
  !             IRES = 1 or -1 represents a normal return, and
  !          DLSODI continues integrating the ODE.  Leave IRES
  !          unchanged from its input value.
  !             IRES = 2 tells DLSODI to immediately return control
  !          to the calling program, with ISTATE = 3.  This lets
  !          the calling program change parameters of the problem,
  !          if necessary.
  !             IRES = 3 represents an error condition (for example, an
  !          illegal value of y).  DLSODI tries to integrate the system
  !          without getting IRES = 3 from RES.  If it cannot, DLSODI
  !          returns with ISTATE = -7 or -1.
  !             On an DLSODI return with ISTATE = 3, -1, or -7, the values
  !          of T and Y returned correspond to the last point reached
  !          successfully without getting the flag IRES = 2 or 3.
  !             The flag values IRES = 2 and 3 should not be used to
  !          handle switches or root-stop conditions.  This is better
  !          done by calling DLSODI in a one-step mode and checking the
  !          stopping function for a sign change at each step.
  !             If quantities computed in the RES routine are needed
  !          externally to DLSODI, an extra call to RES should be made
  !          for this purpose, for consistent and accurate results.
  !          To get the current dy/dt for the S argument, use DINTDY.
  !             RES must be declared External in the calling
  !          program.  See note below for more about RES.
  !
  ! ADDA   = the name of the user-supplied subroutine which adds the
  !          matrix A = A(t,y) to another matrix stored in the same form
  !          as A.  The storage form is determined by MITER (see MF).
  !          This subroutine is to have the form
  !               SUBROUTINE ADDA (NEQ, T, Y, ML, MU, P, NROWP)
  !               DOUBLE PRECISION T, Y(*), P(NROWP,*)
  !          where NEQ, T, Y, ML, MU, and NROWP are input and P is
  !          output.  Y is an array of length NEQ, and the matrix P is
  !          stored in an NROWP by NEQ array.
  !             In the full matrix case ( MITER = 1 or 2) ADDA should
  !          add  A    to P(i,j).  ML and MU are ignored.
  !                i,j
  !             In the band matrix case ( MITER = 4 or 5) ADDA should
  !          add  A    to  P(i-j+MU+1,j).
  !                i,j
  !          See JAC for details on this band storage form.
  !             ADDA must be declared External in the calling program.
  !          See note below for more information about ADDA.
  !
  ! JAC    = the name of the user-supplied subroutine which supplies the
  !          Jacobian matrix, dr/dy, where r = g - A*s.  The form of the
  !          Jacobian matrix is determined by MITER.  JAC is required
  !          if MITER = 1 or 4 -- otherwise a dummy name can be
  !          passed.  This subroutine is to have the form
  !               SUBROUTINE JAC ( NEQ, T, Y, S, ML, MU, P, NROWP )
  !               DOUBLE PRECISION T, Y(*), S(*), P(NROWP,*)
  !          where NEQ, T, Y, S, ML, MU, and NROWP are input and P
  !          is output.  Y and S are arrays of length NEQ, and the
  !          matrix P is stored in an NROWP by NEQ array.
  !          P is to be loaded with partial derivatives (elements
  !          of the Jacobian matrix) on output.
  !             In the full matrix case (MITER = 1), ML and MU
  !          are ignored and the Jacobian is to be loaded into P
  !          by columns-- i.e., dr(i)/dy(j) is loaded into P(i,j).
  !             In the band matrix case (MITER = 4), the elements
  !          within the band are to be loaded into P by columns,
  !          with diagonal lines of dr/dy loaded into the
  !          rows of P.  Thus dr(i)/dy(j) is to be loaded
  !          into P(i-j+MU+1,j).  The locations in P in the two
  !          triangular areas which correspond to nonexistent matrix
  !          elements can be ignored or loaded arbitrarily, as they
  !          they are overwritten by DLSODI.  ML and MU are the
  !          half-bandwidth parameters (see IWORK).
  !               In either case, P is preset to zero by the solver,
  !          so that only the nonzero elements need be loaded by JAC.
  !          Each call to JAC is preceded by a call to RES with the same
  !          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
  !          intermediate quantities shared by both calculations may be
  !          saved in a user Common block by RES and not recomputed by JAC
  !          if desired.  Also, JAC may alter the Y array, if desired.
  !               JAC need not provide dr/dy exactly.  A crude
  !          approximation (possibly with a smaller bandwidth) will do.
  !               JAC must be declared External in the calling program.
  !               See note below for more about JAC.
  !
  !    Note on RES, ADDA, and JAC:
  !          These subroutines may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in the subroutines) and/or Y has length
  !          exceeding NEQ(1).  However, these routines should not alter
  !          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
  !          See the descriptions of NEQ and Y below.
  !
  ! NEQ    = the size of the system (number of first order ordinary
  !          differential equations or scalar algebraic equations).
  !          Used only for input.
  !          NEQ may be decreased, but not increased, during the problem.
  !          If NEQ is decreased (with ISTATE = 3 on input), the
  !          remaining components of Y should be left undisturbed, if
  !          these are to be accessed in RES, ADDA, or JAC.
  !
  !          Normally, NEQ is a scalar, and it is generally referred to
  !          as a scalar in this user interface description.  However,
  !          NEQ may be an array, with NEQ(1) set to the system size.
  !          (The DLSODI package accesses only NEQ(1).)  In either case,
  !          this parameter is passed as the NEQ argument in all calls
  !          to RES, ADDA, and JAC.  Hence, if it is an array,
  !          locations NEQ(2),... may be used to store other integer data
  !          and pass it to RES, ADDA, or JAC.  Each such subroutine
  !          must include NEQ in a Dimension statement in that case.
  !
  ! Y      = a real array for the vector of dependent variables, of
  !          length NEQ or more.  Used for both input and output on the
  !          first call (ISTATE = 0 or 1), and only for output on other
  !          calls.  On the first call, Y must contain the vector of
  !          initial values.  On output, Y contains the computed solution
  !          vector, evaluated at T.  If desired, the Y array may be used
  !          for other purposes between calls to the solver.
  !
  !          This array is passed as the Y argument in all calls to RES,
  !          ADDA, and JAC.  Hence its length may exceed NEQ,
  !          and locations Y(NEQ+1),... may be used to store other real
  !          data and pass it to RES, ADDA, or JAC.  (The DLSODI
  !          package accesses only Y(1),...,Y(NEQ). )
  !
  ! YDOTI  = a real array for the initial value of the vector
  !          dy/dt and for work space, of dimension at least NEQ.
  !
  !          On input:
  !            If ISTATE = 0, then DLSODI will compute the initial value
  !          of dy/dt, if A is nonsingular.  Thus YDOTI will
  !          serve only as work space and may have any value.
  !            If ISTATE = 1, then YDOTI must contain the initial value
  !          of dy/dt.
  !            If ISTATE = 2 or 3 (continuation calls), then YDOTI
  !          may have any value.
  !            Note: If the initial value of A is singular, then
  !          DLSODI cannot compute the initial value of dy/dt, so
  !          it must be provided in YDOTI, with ISTATE = 1.
  !
  !          On output, when DLSODI terminates abnormally with ISTATE =
  !          -1, -4, or -5, YDOTI will contain the residual
  !          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
  !          its initial value, and YDOTI is supplied with ISTATE = 1,
  !          then there may have been an incorrect input value of
  !          YDOTI = dy/dt, or the problem (as given to DLSODI)
  !          may not have a solution.
  !
  !          If desired, the YDOTI array may be used for other
  !          purposes between calls to the solver.
  !
  ! T      = the independent variable.  On input, T is used only on the
  !          first call, as the initial point of the integration.
  !          On output, after each call, T is the value at which a
  !          computed solution Y is evaluated (usually the same as TOUT).
  !          on an error return, T is the farthest point reached.
  !
  ! TOUT   = the next value of t at which a computed solution is desired.
  !          Used only for input.
  !
  !          When starting the problem (ISTATE = 0 or 1), TOUT may be
  !          equal to T for one call, then should .ne. T for the next
  !          call.  For the initial T, an input value of TOUT .ne. T is
  !          used in order to determine the direction of the integration
  !          (i.e. the algebraic sign of the step sizes) and the rough
  !          scale of the problem.  Integration in either direction
  !          (forward or backward in t) is permitted.
  !
  !          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
  !          the first call (i.e. the first call with TOUT .ne. T).
  !          Otherwise, TOUT is required on every call.
  !
  !          If ITASK = 1, 3, or 4, the values of TOUT need not be
  !          monotone, but a value of TOUT which backs up is limited
  !          to the current internal T interval, whose endpoints are
  !          TCUR - HU and TCUR (see optional outputs, below, for
  !          TCUR and HU).
  !
  ! ITOL   = an indicator for the type of error control.  See
  !          description below under ATOL.  Used only for input.
  !
  ! RTOL   = a relative error tolerance parameter, either a scalar or
  !          an array of length NEQ.  See description below under ATOL.
  !          Input only.
  !
  ! ATOL   = an absolute error tolerance parameter, either a scalar or
  !          an array of length NEQ.  Input only.
  !
  !             The input parameters ITOL, RTOL, and ATOL determine
  !          the error control performed by the solver.  The solver will
  !          control the vector E = (E(i)) of estimated local errors
  !          in y, according to an inequality of the form
  !                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
  !          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
  !          and the RMS-norm (root-mean-square norm) here is
  !          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
  !          is a vector of weights which must always be positive, and
  !          the values of RTOL and ATOL should all be non-negative.
  !          The following table gives the types (scalar/array) of
  !          RTOL and ATOL, and the corresponding form of EWT(i).
  !
  !             ITOL    RTOL       ATOL          EWT(i)
  !              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  !              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  !              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  !              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
  !
  !          When either of these parameters is a scalar, it need not
  !          be dimensioned in the user's calling program.
  !
  !          If none of the above choices (with ITOL, RTOL, and ATOL
  !          fixed throughout the problem) is suitable, more general
  !          error controls can be obtained by substituting
  !          user-supplied routines for the setting of EWT and/or for
  !          the norm calculation.  See Part 4 below.
  !
  !          If global errors are to be estimated by making a repeated
  !          run on the same problem with smaller tolerances, then all
  !          components of RTOL and ATOL (i.e. of EWT) should be scaled
  !          down uniformly.
  !
  ! ITASK  = an index specifying the task to be performed.
  !          Input only.  ITASK has the following values and meanings.
  !          1  means normal computation of output values of y(t) at
  !             t = TOUT (by overshooting and interpolating).
  !          2  means take one step only and return.
  !          3  means stop at the first internal mesh point at or
  !             beyond t = TOUT and return.
  !          4  means normal computation of output values of y(t) at
  !             t = TOUT but without overshooting t = TCRIT.
  !             TCRIT must be input as RWORK(1).  TCRIT may be equal to
  !             or beyond TOUT, but not behind it in the direction of
  !             integration.  This option is useful if the problem
  !             has a singularity at or beyond t = TCRIT.
  !          5  means take one step, without passing TCRIT, and return.
  !             TCRIT must be input as RWORK(1).
  !
  !          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
  !          (within roundoff), it will return T = TCRIT (exactly) to
  !          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
  !          in which case answers at t = TOUT are returned first).
  !
  ! ISTATE = an index used for input and output to specify the
  !          state of the calculation.
  !
  !          On input, the values of ISTATE are as follows.
  !          0  means this is the first call for the problem, and
  !             DLSODI is to compute the initial value of dy/dt
  !             (while doing other initializations).  See note below.
  !          1  means this is the first call for the problem, and
  !             the initial value of dy/dt has been supplied in
  !             YDOTI (DLSODI will do other initializations).  See note
  !             below.
  !          2  means this is not the first call, and the calculation
  !             is to continue normally, with no change in any input
  !             parameters except possibly TOUT and ITASK.
  !             (If ITOL, RTOL, and/or ATOL are changed between calls
  !             with ISTATE = 2, the new values will be used but not
  !             tested for legality.)
  !          3  means this is not the first call, and the
  !             calculation is to continue normally, but with
  !             a change in input parameters other than
  !             TOUT and ITASK.  Changes are allowed in
  !             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
  !             and any of the optional inputs except H0.
  !             (See IWORK description for ML and MU.)
  !          Note:  A preliminary call with TOUT = T is not counted
  !          as a first call here, as no initialization or checking of
  !          input is done.  (Such a call is sometimes useful for the
  !          purpose of outputting the initial conditions.)
  !          Thus the first call for which TOUT .ne. T requires
  !          ISTATE = 0 or 1 on input.
  !
  !          On output, ISTATE has the following values and meanings.
  !           0 or 1  means nothing was done; TOUT = t and
  !              ISTATE = 0 or 1 on input.
  !           2  means that the integration was performed successfully.
  !           3  means that the user-supplied Subroutine RES signalled
  !              DLSODI to halt the integration and return (IRES = 2).
  !              Integration as far as T was achieved with no occurrence
  !              of IRES = 2, but this flag was set on attempting the
  !              next step.
  !          -1  means an excessive amount of work (more than MXSTEP
  !              steps) was done on this call, before completing the
  !              requested task, but the integration was otherwise
  !              successful as far as T.  (MXSTEP is an optional input
  !              and is normally 500.)  To continue, the user may
  !              simply reset ISTATE to a value .gt. 1 and call again
  !              (the excess work step counter will be reset to 0).
  !              In addition, the user may increase MXSTEP to avoid
  !              this error return (see below on optional inputs).
  !          -2  means too much accuracy was requested for the precision
  !              of the machine being used.  This was detected before
  !              completing the requested task, but the integration
  !              was successful as far as T.  To continue, the tolerance
  !              parameters must be reset, and ISTATE must be set
  !              to 3.  The optional output TOLSF may be used for this
  !              purpose.  (Note: If this condition is detected before
  !              taking any steps, then an illegal input return
  !              (ISTATE = -3) occurs instead.)
  !          -3  means illegal input was detected, before taking any
  !              integration steps.  See written message for details.
  !              Note:  If the solver detects an infinite loop of calls
  !              to the solver with illegal input, it will cause
  !              the run to stop.
  !          -4  means there were repeated error test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              The problem may have a singularity, or the input
  !              may be inappropriate.
  !          -5  means there were repeated convergence test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              This may be caused by an inaccurate Jacobian matrix.
  !          -6  means EWT(i) became zero for some i during the
  !              integration.  pure relative error control (ATOL(i)=0.0)
  !              was requested on a variable which has now vanished.
  !              the integration was successful as far as T.
  !          -7  means that the user-supplied Subroutine RES set
  !              its error flag (IRES = 3) despite repeated tries by
  !              DLSODI to avoid that condition.
  !          -8  means that ISTATE was 0 on input but DLSODI was unable
  !              to compute the initial value of dy/dt.  See the
  !              printed message for details.
  !
  !          Note:  Since the normal output value of ISTATE is 2,
  !          it does not need to be reset for normal continuation.
  !          Similarly, ISTATE (= 3) need not be reset if RES told
  !          DLSODI to return because the calling program must change
  !          the parameters of the problem.
  !          Also, since a negative input value of ISTATE will be
  !          regarded as illegal, a negative output value requires the
  !          user to change it, and possibly other inputs, before
  !          calling the solver again.
  !
  ! IOPT   = an integer flag to specify whether or not any optional
  !          inputs are being used on this call.  Input only.
  !          The optional inputs are listed separately below.
  !          IOPT = 0 means no optional inputs are being used.
  !                   Default values will be used in all cases.
  !          IOPT = 1 means one or more optional inputs are being used.
  !
  ! RWORK  = a real working array (double precision).
  !          The length of RWORK must be at least
  !             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
  !          NYH    = the initial value of NEQ,
  !          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
  !                   smaller value is given as an optional input),
  !          LENWM   = NEQ**2 + 2    if MITER is 1 or 2, and
  !          LENWM   = (2*ML+MU+1)*NEQ + 2 if MITER is 4 or 5.
  !          (See MF description for the definition of METH and MITER.)
  !          Thus if MAXORD has its default value and NEQ is constant,
  !          this length is
  !             22 + 16*NEQ + NEQ**2         for MF = 11 or 12,
  !             22 + 17*NEQ + (2*ML+MU)*NEQ  for MF = 14 or 15,
  !             22 +  9*NEQ + NEQ**2         for MF = 21 or 22,
  !             22 + 10*NEQ + (2*ML+MU)*NEQ  for MF = 24 or 25.
  !          The first 20 words of RWORK are reserved for conditional
  !          and optional inputs and optional outputs.
  !
  !          The following word in RWORK is a conditional input:
  !            RWORK(1) = TCRIT = critical value of t which the solver
  !                       is not to overshoot.  Required if ITASK is
  !                       4 or 5, and ignored otherwise.  (See ITASK.)
  !
  ! LRW    = the length of the array RWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! IWORK  = an integer work array.  The length of IWORK must be at least
  !          20 + NEQ .  The first few words of IWORK are used for
  !          conditional and optional inputs and optional outputs.
  !
  !          The following 2 words in IWORK are conditional inputs:
  !            IWORK(1) = ML     These are the lower and upper
  !            IWORK(2) = MU     half-bandwidths, respectively, of the
  !                       matrices in the problem-- the Jacobian dr/dy
  !                       and the left-hand side matrix A. These
  !                       half-bandwidths exclude the main diagonal,
  !                       so the total bandwidth is ML + MU + 1 .
  !                       The band is defined by the matrix locations
  !                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
  !                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
  !                       These are required if MITER is 4 or 5, and
  !                       ignored otherwise.
  !                       ML and MU may in fact be the band parameters
  !                       for matrices to which dr/dy and A are only
  !                       approximately equal.
  !
  ! LIW    = the length of the array IWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! Note:  The work arrays must not be altered between calls to DLSODI
  ! for the same problem, except possibly for the conditional and
  ! optional inputs, and except for the last 3*NEQ words of RWORK.
  ! The latter space is used for internal scratch space, and so is
  ! available for use by the user outside DLSODI between calls, if
  ! desired (but not for use by RES, ADDA, or JAC).
  !
  ! MF     = the method flag.  Used only for input.  The legal values of
  !          MF are 11, 12, 14, 15, 21, 22, 24, and 25.
  !          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
  !            METH indicates the basic linear multistep method:
  !              METH = 1 means the implicit Adams method.
  !              METH = 2 means the method based on Backward
  !                       Differentiation Formulas (BDFs).
  !                The BDF method is strongly preferred for stiff
  !              problems, while the Adams method is preferred when
  !              the problem is not stiff.  If the matrix A(t,y) is
  !              nonsingular, stiffness here can be taken to mean that of
  !              the explicit ODE system dy/dt = A-inverse * g.  If A is
  !              singular, the concept of stiffness is not well defined.
  !                If you do not know whether the problem is stiff, we
  !              recommend using METH = 2.  If it is stiff, the advantage
  !              of METH = 2 over METH = 1 will be great, while if it is
  !              not stiff, the advantage of METH = 1 will be slight.
  !              If maximum efficiency is important, some experimentation
  !              with METH may be necessary.
  !            MITER indicates the corrector iteration method:
  !              MITER = 1 means chord iteration with a user-supplied
  !                        full (NEQ by NEQ) Jacobian.
  !              MITER = 2 means chord iteration with an internally
  !                        generated (difference quotient) full Jacobian.
  !                        This uses NEQ+1 extra calls to RES per dr/dy
  !                        evaluation.
  !              MITER = 4 means chord iteration with a user-supplied
  !                        banded Jacobian.
  !              MITER = 5 means chord iteration with an internally
  !                        generated banded Jacobian (using ML+MU+2
  !                        extra calls to RES per dr/dy evaluation).
  !              If MITER = 1 or 4, the user must supply a Subroutine JAC
  !              (the name is arbitrary) as described above under JAC.
  !              For other values of MITER, a dummy argument can be used.
  !-----------------------------------------------------------------------
  ! Optional Inputs.
  !
  ! The following is a list of the optional inputs provided for in the
  ! call sequence.  (See also Part 2.)  For each such input variable,
  ! this table lists its name as used in this documentation, its
  ! location in the call sequence, its meaning, and the default value.
  ! the use of any of these inputs requires IOPT = 1, and in that
  ! case all of these inputs are examined.  A value of zero for any
  ! of these optional inputs will cause the default value to be used.
  ! Thus to use a subset of the optional inputs, simply preload
  ! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
  ! then set those of interest to nonzero values.
  !
  ! Name    Location      Meaning and Default Value
  !
  ! H0      RWORK(5)  the step size to be attempted on the first step.
  !                   The default value is determined by the solver.
  !
  ! HMAX    RWORK(6)  the maximum absolute step size allowed.
  !                   The default value is infinite.
  !
  ! HMIN    RWORK(7)  the minimum absolute step size allowed.
  !                   The default value is 0.  (This lower bound is not
  !                   enforced on the final step before reaching TCRIT
  !                   when ITASK = 4 or 5.)
  !
  ! MAXORD  IWORK(5)  the maximum order to be allowed.  The default
  !                   value is 12 if METH = 1, and 5 if METH = 2.
  !                   If MAXORD exceeds the default value, it will
  !                   be reduced to the default value.
  !                   If MAXORD is changed during the problem, it may
  !                   cause the current order to be reduced.
  !
  ! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
  !                   allowed during one call to the solver.
  !                   The default value is 500.
  !
  ! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
  !                   warning that T + H = T on a step (H = step size).
  !                   This must be positive to result in a non-default
  !                   value.  The default value is 10.
  !-----------------------------------------------------------------------
  ! Optional Outputs.
  !
  ! As optional additional output from DLSODI, the variables listed
  ! below are quantities related to the performance of DLSODI
  ! which are available to the user.  These are communicated by way of
  ! the work arrays, but also have internal mnemonic names as shown.
  ! Except where stated otherwise, all of these outputs are defined
  ! on any successful return from DLSODI, and on any return with
  ! ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
  ! input) or -8, they will be unchanged from their existing values
  ! (if any), except possibly for TOLSF, LENRW, and LENIW.
  ! On any error return, outputs relevant to the error will be defined,
  ! as noted below.
  !
  ! Name    Location      Meaning
  !
  ! HU      RWORK(11) the step size in t last used (successfully).
  !
  ! HCUR    RWORK(12) the step size to be attempted on the next step.
  !
  ! TCUR    RWORK(13) the current value of the independent variable
  !                   which the solver has actually reached, i.e. the
  !                   current internal mesh point in t.  On output, TCUR
  !                   will always be at least as far as the argument
  !                   T, but may be farther (if interpolation was done).
  !
  ! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
  !                   computed when a request for too much accuracy was
  !                   detected (ISTATE = -3 if detected at the start of
  !                   the problem, ISTATE = -2 otherwise).  If ITOL is
  !                   left unaltered but RTOL and ATOL are uniformly
  !                   scaled up by a factor of TOLSF for the next call,
  !                   then the solver is deemed likely to succeed.
  !                   (The user may also ignore TOLSF and alter the
  !                   tolerance parameters in any other way appropriate.)
  !
  ! NST     IWORK(11) the number of steps taken for the problem so far.
  !
  ! NRE     IWORK(12) the number of residual evaluations (RES calls)
  !                   for the problem so far.
  !
  ! NJE     IWORK(13) the number of Jacobian evaluations (each involving
  !                   an evaluation of A and dr/dy) for the problem so
  !                   far.  This equals the number of calls to ADDA and
  !                   (if MITER = 1 or 4) JAC, and the number of matrix
  !                   LU decompositions.
  !
  ! NQU     IWORK(14) the method order last used (successfully).
  !
  ! NQCUR   IWORK(15) the order to be attempted on the next step.
  !
  ! IMXER   IWORK(16) the index of the component of largest magnitude in
  !                   the weighted local error vector ( E(i)/EWT(i) ),
  !                   on an error return with ISTATE = -4 or -5.
  !
  ! LENRW   IWORK(17) the length of RWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! LENIW   IWORK(18) the length of IWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  !
  ! The following two arrays are segments of the RWORK array which
  ! may also be of interest to the user as optional outputs.
  ! For each array, the table below gives its internal name,
  ! its base address in RWORK, and its description.
  !
  ! Name    Base Address      Description
  !
  ! YH      21             the Nordsieck history array, of size NYH by
  !                        (NQCUR + 1), where NYH is the initial value
  !                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
  !                        of YH contains HCUR**j/factorial(j) times
  !                        the j-th derivative of the interpolating
  !                        polynomial currently representing the solution,
  !                        evaluated at t = TCUR.
  !
  ! ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
  !                        corrections on each step, scaled on output to
  !                        represent the estimated local error in y on the
  !                        last step. This is the vector E in the descrip-
  !                        tion of the error control.  It is defined only
  !                        on a return from DLSODI with ISTATE = 2.
  !
  !-----------------------------------------------------------------------
  ! Part 2.  Other Routines Callable.
  !
  ! The following are optional calls which the user may make to
  ! gain additional capabilities in conjunction with DLSODI.
  ! (The routines XSETUN and XSETF are designed to conform to the
  ! SLATEC error handling package.)
  !
  !     Form of Call                  Function
  !   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
  !                             output of messages from DLSODI, if
  !                             the default is not desired.
  !                             The default value of LUN is 6.
  !
  !   CALL XSETF(MFLAG)         Set a flag to control the printing of
  !                             messages by DLSODI.
  !                             MFLAG = 0 means do not print. (Danger:
  !                             This risks losing valuable information.)
  !                             MFLAG = 1 means print (the default).
  !
  !                             Either of the above calls may be made at
  !                             any time and will take effect immediately.
  !
  !   CALL DSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
  !                             the internal Common blocks used by
  !                             DLSODI (see Part 3 below).
  !                             RSAV must be a real array of length 218
  !                             or more, and ISAV must be an integer
  !                             array of length 37 or more.
  !                             JOB=1 means save Common into RSAV/ISAV.
  !                             JOB=2 means restore Common from RSAV/ISAV.
  !                                DSRCOM is useful if one is
  !                             interrupting a run and restarting
  !                             later, or alternating between two or
  !                             more problems solved with DLSODI.
  !
  !   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
  !        (see below)          orders, at a specified point t, if
  !                             desired.  It may be called only after
  !                             a successful return from DLSODI.
  !
  ! The detailed instructions for using DINTDY are as follows.
  ! The form of the call is:
  !
  !   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
  !
  ! The input parameters are:
  !
  ! T         = value of independent variable where answers are desired
  !             (normally the same as the T last returned by DLSODI).
  !             For valid results, T must lie between TCUR - HU and TCUR.
  !             (See optional outputs for TCUR and HU.)
  ! K         = integer order of the derivative desired.  K must satisfy
  !             0 .le. K .le. NQCUR, where NQCUR is the current order
  !             (see optional outputs).  The capability corresponding
  !             to K = 0, i.e. computing y(T), is already provided
  !             by DLSODI directly.  Since NQCUR .ge. 1, the first
  !             derivative dy/dt is always available with DINTDY.
  ! RWORK(21) = the base address of the history array YH.
  ! NYH       = column length of YH, equal to the initial value of NEQ.
  !
  ! The output parameters are:
  !
  ! DKY       = a real array of length NEQ containing the computed value
  !             of the K-th derivative of y(t).
  ! IFLAG     = integer flag, returned as 0 if K and T were legal,
  !             -1 if K was illegal, and -2 if T was illegal.
  !             On an error return, a message is also written.
  !-----------------------------------------------------------------------
  ! Part 3.  Common Blocks.
  !
  ! If DLSODI is to be used in an overlay situation, the user
  ! must declare, in the primary overlay, the variables in:
  !   (1) the call sequence to DLSODI, and
  !   (2) the internal Common block
  !         /DLS001/  of length  255  (218 double precision words
  !                      followed by 37 integer words),
  !
  ! If DLSODI is used on a system in which the contents of internal
  ! Common blocks are not preserved between calls, the user should
  ! declare the above Common block in the calling program to insure
  ! that their contents are preserved.
  !
  ! If the solution of a given problem by DLSODI is to be interrupted
  ! and then later continued, such as when restarting an interrupted run
  ! or alternating between two or more problems, the user should save,
  ! following the return from the last DLSODI call prior to the
  ! interruption, the contents of the call sequence variables and the
  ! internal Common blocks, and later restore these values before the
  ! next DLSODI call for that problem.  To save and restore the Common
  ! blocks, use Subroutine DSRCOM (see Part 2 above).
  !
  !-----------------------------------------------------------------------
  ! Part 4.  Optionally Replaceable Solver Routines.
  !
  ! Below are descriptions of two routines in the DLSODI package which
  ! relate to the measurement of errors.  Either routine can be
  ! replaced by a user-supplied version, if desired.  However, since such
  ! a replacement may have a major impact on performance, it should be
  ! done only when absolutely necessary, and only with great caution.
  ! (Note: The means by which the package version of a routine is
  ! superseded by the user's version may be system-dependent.)
  !
  ! (a) DEWSET.
  ! The following subroutine is called just before each internal
  ! integration step, and sets the array of error weights, EWT, as
  ! described under ITOL/RTOL/ATOL above:
  !     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
  ! where NEQ, ITOL, RTOL, and ATOL are as in the DLSODI call sequence,
  ! YCUR contains the current dependent variable vector, and
  ! EWT is the array of weights set by DEWSET.
  !
  ! If the user supplies this subroutine, it must return in EWT(i)
  ! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
  ! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
  ! routine (see below), and also used by DLSODI in the computation
  ! of the optional output IMXER, the diagonal Jacobian approximation,
  ! and the increments for difference quotient Jacobians.
  !
  ! In the user-supplied version of DEWSET, it may be desirable to use
  ! the current values of derivatives of y.  Derivatives up to order NQ
  ! are available from the history array YH, described above under
  ! optional outputs.  In DEWSET, YH is identical to the YCUR array,
  ! extended to NQ + 1 columns with a column length of NYH and scale
  ! factors of H**j/factorial(j).  On the first call for the problem,
  ! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
  ! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
  ! can be obtained by including in DEWSET the statements:
  !     DOUBLE PRECISION RLS
  !     COMMON /DLS001/ RLS(218),ILS(37)
  !     NQ = ILS(33)
  !     NST = ILS(34)
  !     H = RLS(212)
  ! Thus, for example, the current value of dy/dt can be obtained as
  ! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
  ! unnecessary when NST = 0).
  !
  ! (b) DVNORM.
  ! The following is a real function routine which computes the weighted
  ! root-mean-square norm of a vector v:
  !     D = DVNORM (N, V, W)
  ! where:
  !   N = the length of the vector,
  !   V = real array of length N containing the vector,
  !   W = real array of length N containing weights,
  !   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
  ! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
  ! EWT is as set by Subroutine DEWSET.
  !
  ! If the user supplies this function, it should return a non-negative
  ! value of DVNORM suitable for use in the error control in DLSODI.
  ! None of the arguments should be altered by DVNORM.
  ! For example, a user-supplied DVNORM routine might:
  !   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
  !   -ignore some components of V in the norm, with the effect of
  !    suppressing the error control on those components of y.
  !-----------------------------------------------------------------------
  !
  !***REVISION HISTORY  (YYYYMMDD)
  ! 19800424  DATE WRITTEN
  ! 19800519  Corrected access of YH on forced order reduction;
  !           numerous corrections to prologues and other comments.
  ! 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
  !           minor corrections to main prologue.
  ! 19800903  Corrected ISTATE logic; minor changes in prologue.
  ! 19800923  Added zero initialization of HU and NQU.
  ! 19801028  Reorganized RES calls in AINVG, STODI, and PREPJI;
  !           in LSODI, corrected NRE increment and reset LDY0 at 580;
  !           numerous corrections to main prologue.
  ! 19801218  Revised XERRWD routine; minor corrections to main prologue.
  ! 19810330  Added Common block /LSI001/; use LSODE's INTDY and SOLSY;
  !           minor corrections to XERRWD and error message at 604;
  !           minor corrections to declarations; corrections to prologues.
  ! 19810818  Numerous revisions: replaced EWT by 1/EWT; used flags
  !           JCUR, ICF, IERPJ, IERSL between STODI and subordinates;
  !           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
  !           reorganized returns from STODI; reorganized type decls.;
  !           fixed message length in XERRWD; changed default LUNIT to 6;
  !           changed Common lengths; changed comments throughout.
  ! 19820906  Corrected use of ABS(H) in STODI; minor comment fixes.
  ! 19830510  Numerous revisions: revised diff. quotient increment;
  !           eliminated block /LSI001/, using IERPJ flag;
  !           revised STODI logic after PJAC return;
  !           revised tuning of H change and step attempts in STODI;
  !           corrections to main prologue and internal comments.
  ! 19870330  Major update: corrected comments throughout;
  !           removed TRET from Common; rewrote EWSET with 4 loops;
  !           fixed t test in INTDY; added Cray directives in STODI;
  !           in STODI, fixed DELP init. and logic around PJAC call;
  !           combined routines to save/restore Common;
  !           passed LEVEL = 0 in error message calls (except run abort).
  ! 20010425  Major update: convert source lines to upper case;
  !           added *DECK lines; changed from 1 to * in dummy dimensions;
  !           changed names R1MACH/D1MACH to RUMACH/DUMACH;
  !           renamed routines for uniqueness across single/double prec.;
  !           converted intrinsic names to generic form;
  !           removed ILLIN and NTREP (data loaded) from Common;
  !           removed all 'own' variables from Common;
  !           changed error messages to quoted strings;
  !           replaced XERRWV/XERRWD with 1993 revised version;
  !           converted prologues, comments, error messages to mixed case;
  !           converted arithmetic IF statements to logical IF statements;
  !           numerous corrections to prologues and internal comments.
  ! 20010507  Converted single precision source to double precision.
  ! 20020502  Corrected declarations in descriptions of user routines.
  ! 20031105  Restored 'own' variables to Common block, to enable
  !           interrupt/restart feature.
  ! 20031112  Added SAVE statements for data-loaded constants.
  ! 20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
  !
  !-----------------------------------------------------------------------
  ! Other routines in the DLSODI package.
  !
  ! In addition to Subroutine DLSODI, the DLSODI package includes the
  ! following subroutines and function routines:
  !  DAINVG   computes the initial value of the vector
  !             dy/dt = A-inverse * g
  !  DINTDY   computes an interpolated value of the y vector at t = TOUT.
  !  DSTODI   is the core integrator, which does one step of the
  !           integration and the associated error control.
  !  DCFODE   sets all method coefficients and test constants.
  !  DPREPJI  computes and preprocesses the Jacobian matrix
  !           and the Newton iteration matrix P.
  !  DSOLSY   manages solution of linear system in chord iteration.
  !  DEWSET   sets the error weight vector EWT before each step.
  !  DVNORM   computes the weighted RMS-norm of a vector.
  !  DSRCOM   is a user-callable routine to save and restore
  !           the contents of the internal Common blocks.
  !  DGEFA and DGESL   are routines from LINPACK for solving full
  !           systems of linear algebraic equations.
  !  DGBFA and DGBSL   are routines from LINPACK for solving banded
  !           linear systems.
  !  DUMACH   computes the unit roundoff in a machine-independent manner.
  !  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
  !           error messages and warnings.  XERRWD is machine-dependent.
  ! Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
  ! All the others are subroutines.
  !
  !-----------------------------------------------------------------------
        !EXTERNAL DPREPJI , DSOLSY
        !DOUBLE PRECISION DUMACH , DVNORM
        INTEGER INIt , MXStep , MXHnil , NHNil , NSLast , NYH , IOWns ,   &
       &        ICF , IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH ,   &
       &        LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,        &
       &        MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST , NFE ,     &
       &        NJE , NQU
        INTEGER i , i1 , i2 , ier , iflag , imxer , ires , kgo , leniw ,  &
       &        lenrw , lenwm , lp , lyd0 , ml , mord , mu , mxhnl0 ,     &
       &        mxstp0
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION atoli , ayi , big , ewti , h0 , hmax , hmx , rh ,&
       &                 rtoli , tcrit , tdist , tnext , tol , tolsf ,    &
       &                 tp , size , sum , w0
        DIMENSION mord(2)
        LOGICAL ihit
        CHARACTER*60 msg
        SAVE mord , mxstp0 , mxhnl0
  !-----------------------------------------------------------------------
  ! The following internal Common block contains
  ! (a) variables which are local to any subroutine but whose values must
  !     be preserved between calls to the routine ("own" variables), and
  ! (b) variables which are communicated between subroutines.
  ! The block DLS001 is declared in subroutines DLSODI, DINTDY, DSTODI,
  ! DPREPJI, and DSOLSY.
  ! Groups of variables are replaced by dummy arrays in the Common
  ! declarations in routines where those variables are not used.
  !-----------------------------------------------------------------------
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , INIt , MXStep , MXHnil ,       &
       &                NHNil , NSLast , NYH , IOWns(6) , ICF , IERpj ,   &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
        DATA mord(1) , mord(2)/12 , 5/ , mxstp0/500/ , mxhnl0/10/
  !-----------------------------------------------------------------------
  ! Block A.
  ! This code block is executed on every call.
  ! It tests ISTATE and ITASK for legality and branches appropriately.
  ! If ISTATE .gt. 1 but the flag INIT shows that initialization has
  ! not yet been done, an error return occurs.
  ! If ISTATE = 0 or 1 and TOUT = T, return immediately.
  !-----------------------------------------------------------------------
        IF ( Istate.LT.0 .OR. Istate.GT.3 ) THEN
  !-----------------------------------------------------------------------
  ! Block I.
  ! The following block handles all error returns due to illegal input
  ! (ISTATE = -3), as detected before calling the core integrator.
  ! First the error message routine is called.  If the illegal input
  ! is a negative ISTATE, the run is aborted (apparent infinite loop).
  !-----------------------------------------------------------------------
           msg = 'DLSODI-  ISTATE (=I1) illegal.'
           CALL XERRWD(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
           IF ( Istate.GE.0 ) GOTO 3000
  !
           msg = 'DLSODI-  Run aborted.. apparent infinite loop.    '
           CALL XERRWD(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
           GOTO 99999
        ELSE
           IF ( Itask.LT.1 .OR. Itask.GT.5 ) THEN
              msg = 'DLSODI-  ITASK (=I1) illegal. '
              CALL XERRWD(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
              GOTO 3000
           ELSE
              IF ( Istate.LE.1 ) THEN
                 INIt = 0
                 IF ( Tout.EQ.T ) RETURN
              ELSEIF ( INIt.EQ.0 ) THEN
                 msg =                                                    &
       &              'DLSODI-  ISTATE .gt. 1 but DLSODI not initialized.'
                 CALL XERRWD(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
                 GOTO 3000
              ELSEIF ( Istate.EQ.2 ) THEN
                 GOTO 50
              ENDIF
  !-----------------------------------------------------------------------
  ! Block B.
  ! The next code block is executed for the initial call (ISTATE = 0 or 1)
  ! or for a continuation call with parameter changes (ISTATE = 3).
  ! It contains checking of all inputs and various initializations.
  !
  ! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
  ! MF, ML, and MU.
  !-----------------------------------------------------------------------
              IF ( Neq(1).LE.0 ) THEN
                 msg = 'DLSODI-  NEQ (=I1) .lt. 1     '
                 CALL XERRWD(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
                 GOTO 3000
              ELSE
                 IF ( Istate.GT.1 ) THEN
                    IF ( Neq(1).GT.N ) THEN
                       msg =                                              &
       &              'DLSODI-  ISTATE = 3 and NEQ increased (I1 to I2). '
                       CALL XERRWD(msg,50,5,0,2,N,Neq(1),0,0.0D0,0.0D0)
                       GOTO 3000
                    ENDIF
                 ENDIF
                 N = Neq(1)
                 IF ( Itol.LT.1 .OR. Itol.GT.4 ) THEN
                    msg = 'DLSODI-  ITOL (=I1) illegal.  '
                    CALL XERRWD(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
                    GOTO 3000
                 ELSEIF ( Iopt.LT.0 .OR. Iopt.GT.1 ) THEN
                    msg = 'DLSODI-  IOPT (=I1) illegal.  '
                    CALL XERRWD(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
                    GOTO 3000
                 ELSE
                    METh = Mf/10
                    MITer = Mf - 10*METh
                    IF ( METh.LT.1 .OR. METh.GT.2 ) GOTO 2300
                    IF ( MITer.LE.0 .OR. MITer.GT.5 ) GOTO 2300
                    IF ( MITer.EQ.3 ) GOTO 2300
                    IF ( MITer.GE.3 ) THEN
                       ml = Iwork(1)
                       mu = Iwork(2)
                       IF ( ml.LT.0 .OR. ml.GE.N ) THEN
                          msg =                                           &
       &              'DLSODI-  ML(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
                          CALL XERRWD(msg,50,9,0,2,ml,Neq(1),0,0.0D0,     &
       &                              0.0D0)
                          GOTO 3000
                       ELSEIF ( mu.LT.0 .OR. mu.GE.N ) THEN
                          msg =                                           &
       &              'DLSODI-  MU(=I1) illegal: .lt. 0 or .ge. NEQ(=I2) '
                          CALL XERRWD(msg,50,10,0,2,mu,Neq(1),0,0.0D0,    &
       &                              0.0D0)
                          GOTO 3000
                       ENDIF
                    ENDIF
  ! Next process and check the optional inputs. --------------------------
                    IF ( Iopt.EQ.1 ) THEN
                       MAXord = Iwork(5)
                       IF ( MAXord.LT.0 ) THEN
                          msg = 'DLSODI-  MAXORD (=I1) .lt. 0  '
                          CALL XERRWD(msg,30,11,0,1,MAXord,0,0,0.0D0,     &
       &                              0.0D0)
                          GOTO 3000
                       ELSE
                          IF ( MAXord.EQ.0 ) MAXord = 100
                          MAXord = MIN(MAXord,mord(METh))
                          MXStep = Iwork(6)
                          IF ( MXStep.LT.0 ) THEN
                             msg = 'DLSODI-  MXSTEP (=I1) .lt. 0  '
                             CALL XERRWD(msg,30,12,0,1,MXStep,0,0,0.0D0,  &
       &                                 0.0D0)
                             GOTO 3000
                          ELSE
                             IF ( MXStep.EQ.0 ) MXStep = mxstp0
                             MXHnil = Iwork(7)
                             IF ( MXHnil.LT.0 ) THEN
                                msg = 'DLSODI-  MXHNIL (=I1) .lt. 0  '
                                CALL XERRWD(msg,30,13,0,1,MXHnil,0,0,     &
       &                           0.0D0,0.0D0)
                                GOTO 3000
                             ELSE
                                IF ( MXHnil.EQ.0 ) MXHnil = mxhnl0
                                IF ( Istate.LE.1 ) THEN
                                   h0 = Rwork(5)
                                   IF ( (Tout-T)*h0.LT.0.0D0 ) THEN
                                      msg =                               &
       &                        'DLSODI-  TOUT (=R1) behind T (=R2)      '
                                      CALL XERRWD(msg,40,14,0,0,0,0,2,    &
       &                                 Tout,T)
                                      msg =                               &
       &              '      Integration direction is given by H0 (=R1)  '
                                      CALL XERRWD(msg,50,14,0,0,0,0,1,h0, &
       &                                 0.0D0)
                                      GOTO 3000
                                   ENDIF
                                ENDIF
                                hmax = Rwork(6)
                                IF ( hmax.LT.0.0D0 ) THEN
                                   msg = 'DLSODI-  HMAX (=R1) .lt. 0.0  '
                                   CALL XERRWD(msg,30,15,0,0,0,0,1,hmax,  &
       &                              0.0D0)
                                   GOTO 3000
                                ELSE
                                   HMXi = 0.0D0
                                   IF ( hmax.GT.0.0D0 ) HMXi = 1.0D0/hmax
                                   HMIn = Rwork(7)
                                   IF ( HMIn.LT.0.0D0 ) THEN
                                      msg =                               &
       &                                 'DLSODI-  HMIN (=R1) .lt. 0.0  '
                                      CALL XERRWD(msg,30,16,0,0,0,0,1,    &
       &                                 HMIn,0.0D0)
                                      GOTO 3000
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ELSE
                       MAXord = mord(METh)
                       MXStep = mxstp0
                       MXHnil = mxhnl0
                       IF ( Istate.LE.1 ) h0 = 0.0D0
                       HMXi = 0.0D0
                       HMIn = 0.0D0
                    ENDIF
  !-----------------------------------------------------------------------
  ! Set work array pointers and check lengths LRW and LIW.
  ! Pointers to segments of RWORK and IWORK are named by prefixing L to
  ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
  ! Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
  !-----------------------------------------------------------------------
                    LYH = 21
                    IF ( Istate.LE.1 ) NYH = N
                    LWM = LYH + (MAXord+1)*NYH
                    IF ( MITer.LE.2 ) lenwm = N*N + 2
                    IF ( MITer.GE.4 ) lenwm = (2*ml+mu+1)*N + 2
                    LEWt = LWM + lenwm
                    LSAvf = LEWt + N
                    LACor = LSAvf + N
                    lenrw = LACor + N - 1
                    Iwork(17) = lenrw
                    LIWm = 1
                    leniw = 20 + N
                    Iwork(18) = leniw
                    IF ( lenrw.GT.Lrw ) THEN
                       msg =                                              &
       &    'DLSODI-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                       CALL XERRWD(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,0.0D0)
                       GOTO 3000
                    ELSEIF ( leniw.GT.Liw ) THEN
                       msg =                                              &
       &    'DLSODI-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                       CALL XERRWD(msg,60,18,0,2,leniw,Liw,0,0.0D0,0.0D0)
                       GOTO 3000
                    ELSE
  ! Check RTOL and ATOL for legality. ------------------------------------
                       rtoli = Rtol(1)
                       atoli = Atol(1)
                       DO i = 1 , N
                          IF ( Itol.GE.3 ) rtoli = Rtol(i)
                          IF ( Itol.EQ.2 .OR. Itol.EQ.4 ) atoli = Atol(i)
                          IF ( rtoli.LT.0.0D0 ) GOTO 2400
                          IF ( atoli.LT.0.0D0 ) GOTO 2500
                       ENDDO
                       IF ( Istate.LE.1 ) THEN
  !-----------------------------------------------------------------------
  ! Block C.
  ! The next block is for the initial call only (ISTATE = 0 or 1).
  ! It contains all remaining initializations, the call to DAINVG
  ! (if ISTATE = 1), and the calculation of the initial step size.
  ! The error weights in EWT are inverted after being loaded.
  !-----------------------------------------------------------------------
                          UROund = DUMACH()
                          TN = T
                          IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
                             tcrit = Rwork(1)
                             IF ( (tcrit-Tout)*(Tout-T).LT.0.0D0 )        &
       &                          GOTO 2800
                             IF ( h0.NE.0.0D0 .AND. (T+h0-tcrit)          &
       &                          *h0.GT.0.0D0 ) h0 = tcrit - T
                          ENDIF
                          JSTart = 0
                          Rwork(LWM) = SQRT(UROund)
                          NHNil = 0
                          NST = 0
                          NFE = 0
                          NJE = 0
                          NSLast = 0
                          HU = 0.0D0
                          NQU = 0
                          CCMax = 0.3D0
                          MAXcor = 3
                          MSBp = 20
                          MXNcf = 10
  ! Compute initial dy/dt, if necessary, and load it and initial Y into YH
                          lyd0 = LYH + NYH
                          lp = LWM + 1
                          IF ( Istate.EQ.1 ) THEN
  ! Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
                             DO i = 1 , N
                                Rwork(i+LYH-1) = Y(i)
                                Rwork(i+lyd0-1) = Ydoti(i)
                             ENDDO
                          ELSE
  ! DLSODI must compute initial dy/dt (LYD0 points to YH(*,2)). ----------
                             CALL DAINVG(RES,ADDA,Neq(1),T,Y,Rwork(lyd0),    &
       &                                 MITer,ml,mu,Rwork(lp),Iwork(21), &
       &                                 ier)
                             NFE = NFE + 1
                             IF ( ier.LT.0 ) THEN
  ! DAINVG failed because matrix A was singular. -------------------------
                                ier = -ier
                                msg =                                     &
       &    'DLSODI- Attempt to initialize dy/dt failed:  Matrix A is    '
                                CALL XERRWD(msg,60,207,0,0,0,0,0,0.0D0,   &
       &                           0.0D0)
                                msg =                                     &
       &              '      singular.  DGEFA or DGBFA returned INFO = I1'
                                CALL XERRWD(msg,50,207,0,1,ier,0,0,0.0D0, &
       &                           0.0D0)
                                Istate = -8
                                RETURN
                             ELSEIF ( ier.GT.0 ) THEN
  ! DAINVG failed because RES set IRES to 2 or 3. ------------------------
                                msg =                                     &
       &              'DLSODI-  Attempt to initialize dy/dt failed       '
                                CALL XERRWD(msg,50,208,0,0,0,0,0,0.0D0,   &
       &                           0.0D0)
                                msg =                                     &
       &              '      because residual routine set its error flag '
                                CALL XERRWD(msg,50,208,0,0,0,0,0,0.0D0,   &
       &                           0.0D0)
                                msg = '      to IRES = (I1)'
                                CALL XERRWD(msg,20,208,0,1,ier,0,0,0.0D0, &
       &                           0.0D0)
                                Istate = -8
                                RETURN
                             ELSE
                                DO i = 1 , N
                                   Rwork(i+LYH-1) = Y(i)
                                ENDDO
                             ENDIF
                          ENDIF
  ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                          NQ = 1
                          H = 1.0D0
                          CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),        &
       &                              Rwork(LEWt))
                          DO i = 1 , N
                             IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 2600
                             Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
                          ENDDO
  !-----------------------------------------------------------------------
  ! The coding below computes the step size, H0, to be attempted on the
  ! first step, unless the user has supplied a value for this.
  ! First check that TOUT - T differs significantly from zero.
  ! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
  ! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
  ! so as to be between 100*UROUND and 1.0E-3.
  ! Then the computed value H0 is given by..
  !                                      NEQ
  !   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
  !                                       1
  ! where   w0      = MAX ( ABS(T), ABS(TOUT) ),
  !         YDOT(i) = i-th component of initial value of dy/dt,
  !         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
  ! The sign of H0 is inferred from the initial values of TOUT and T.
  !-----------------------------------------------------------------------
                          IF ( h0.EQ.0.0D0 ) THEN
                             tdist = ABS(Tout-T)
                             w0 = MAX(ABS(T),ABS(Tout))
                             IF ( tdist.LT.2.0D0*UROund*w0 ) THEN
                                msg =                                     &
       &    'DLSODI-  TOUT(=R1) too close to T(=R2) to start integration.'
                                CALL XERRWD(msg,60,22,0,0,0,0,2,Tout,T)
                                GOTO 3000
                             ELSE
                                tol = Rtol(1)
                                IF ( Itol.GT.2 ) THEN
                                   DO i = 1 , N
                                      tol = MAX(tol,Rtol(i))
                                   ENDDO
                                ENDIF
                                IF ( tol.LE.0.0D0 ) THEN
                                   atoli = Atol(1)
                                   DO i = 1 , N
                                      IF ( Itol.EQ.2 .OR. Itol.EQ.4 )     &
       &                                 atoli = Atol(i)
                                      ayi = ABS(Y(i))
                                      IF ( ayi.NE.0.0D0 )                 &
       &                                 tol = MAX(tol,atoli/ayi)
                                   ENDDO
                                ENDIF
                                tol = MAX(tol,100.0D0*UROund)
                                tol = MIN(tol,0.001D0)
                                sum = DVNORM(N,Rwork(lyd0),Rwork(LEWt))
                                sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                                h0 = 1.0D0/SQRT(sum)
                                h0 = MIN(h0,tdist)
                                h0 = SIGN(h0,Tout-T)
                             ENDIF
                          ENDIF
  ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                          rh = ABS(h0)*HMXi
                          IF ( rh.GT.1.0D0 ) h0 = h0/rh
  ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                          H = h0
                          DO i = 1 , N
                             Rwork(i+lyd0-1) = h0*Rwork(i+lyd0-1)
                          ENDDO
                          GOTO 700
                       ELSE
  ! If ISTATE = 3, set flag to signal parameter changes to DSTODI. -------
                          JSTart = -1
                          IF ( NQ.GT.MAXord ) THEN
  ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
                             DO i = 1 , N
                                Ydoti(i) = Rwork(i+LWM-1)
                             ENDDO
                          ENDIF
  ! Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
                          Rwork(LWM) = SQRT(UROund)
                          IF ( N.NE.NYH ) THEN
  ! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                             i1 = LYH + L*NYH
                             i2 = LYH + (MAXord+1)*NYH - 1
                             IF ( i1.LE.i2 ) THEN
                                DO i = i1 , i2
                                   Rwork(i) = 0.0D0
                                ENDDO
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  ! Block D.
  ! The next code block is for continuation calls only (ISTATE = 2 or 3)
  ! and is to check stop conditions before taking a step.
  !-----------------------------------------------------------------------
   50      NSLast = NST
           GOTO (100,600,200,300,400) , Itask
        ENDIF
   100  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2900
        T = Tout
        GOTO 1400
   200  tp = TN - HU*(1.0D0+100.0D0*UROund)
        IF ( (tp-Tout)*H.GT.0.0D0 ) THEN
           msg =                                                          &
       &    'DLSODI-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
           CALL XERRWD(msg,60,23,0,1,Itask,0,2,Tout,tp)
           GOTO 3000
        ELSE
           IF ( (TN-Tout)*H.GE.0.0D0 ) GOTO 1300
           GOTO 600
        ENDIF
   300  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2700
        IF ( (tcrit-Tout)*H.LT.0.0D0 ) GOTO 2800
        IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 500
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 2900
        T = Tout
        GOTO 1400
   400  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2700
   500  hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
        IF ( ihit ) GOTO 1300
        tnext = TN + H*(1.0D0+4.0D0*UROund)
        IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
           H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
           IF ( Istate.EQ.2 ) JSTart = -2
        ENDIF
  !-----------------------------------------------------------------------
  ! Block E.
  ! The next block is normally executed for all calls and contains
  ! the call to the one-step core integrator DSTODI.
  !
  ! This is a looping point for the integration steps.
  !
  ! First check for too many steps being taken, update EWT (if not at
  ! start of problem), check for too much accuracy being requested, and
  ! check for H below the roundoff level in T.
  !-----------------------------------------------------------------------
   600  IF ( (NST-NSLast).GE.MXStep ) THEN
  !-----------------------------------------------------------------------
  ! Block H.
  ! The following block handles all unsuccessful returns other than
  ! those for illegal input.  First the error message routine is called.
  ! If there was an error test or convergence test failure, IMXER is set.
  ! Then Y is loaded from YH and T is set to TN.
  ! The optional outputs are loaded into the work arrays before returning.
  !-----------------------------------------------------------------------
  ! The maximum number of steps was taken before reaching TOUT. ----------
           msg = 'DLSODI-  At current T (=R1), MXSTEP (=I1) steps   '
           CALL XERRWD(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      taken on this call before reaching TOUT     '
           CALL XERRWD(msg,50,201,0,1,MXStep,0,1,TN,0.0D0)
           Istate = -1
           GOTO 2000
        ELSE
           CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
           DO i = 1 , N
              IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 1500
              Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
           ENDDO
        ENDIF
   700  tolsf = UROund*DVNORM(N,Rwork(LYH),Rwork(LEWt))
        IF ( tolsf.LE.1.0D0 ) THEN
           IF ( (TN+H).EQ.TN ) THEN
              NHNil = NHNil + 1
              IF ( NHNil.LE.MXHnil ) THEN
                 msg =                                                    &
       &              'DLSODI-  Warning..Internal T (=R1) and H (=R2) are'
                 CALL XERRWD(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &    '      such that in the machine, T + H = T on the next step  '
                 CALL XERRWD(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &              '     (H = step size). Solver will continue anyway.'
                 CALL XERRWD(msg,50,101,0,0,0,0,2,TN,H)
                 IF ( NHNil.GE.MXHnil ) THEN
                    msg =                                                 &
       &              'DLSODI-  Above warning has been issued I1 times.  '
                    CALL XERRWD(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
                    msg =                                                 &
       &              '     It will not be issued again for this problem.'
                    CALL XERRWD(msg,50,102,0,1,MXHnil,0,0,0.0D0,0.0D0)
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  !     CALL DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
  !                 ADDA,JAC,DPREPJI,DSOLSY)
  ! Note: SAVF in DSTODI occupies the same space as YDOTI in DLSODI.
  !-----------------------------------------------------------------------
           CALL DSTODI(Neq,Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),Ydoti, &
       &               Rwork(LSAvf),Rwork(LACor),Rwork(LWM),Iwork(LIWm),  &
       &               RES,ADDA,JAC,DPREPJI,DSOLSY)
           kgo = 1 - KFLag
           GOTO (800,1600,1700,1300,1800) , kgo
        ELSE
           tolsf = tolsf*2.0D0
           IF ( NST.EQ.0 ) THEN
              msg = 'DLSODI-  At start of problem, too much accuracy   '
              CALL XERRWD(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
              msg =                                                       &
       &    '      requested for precision of machine..  See TOLSF (=R1) '
              CALL XERRWD(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
              Rwork(14) = tolsf
              GOTO 3000
           ELSE
  ! Too much accuracy requested for machine precision. -------------------
              msg = 'DLSODI-  At T (=R1), too much accuracy requested  '
              CALL XERRWD(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
              msg = '      for precision of machine..  See TOLSF (=R2) '
              CALL XERRWD(msg,50,203,0,0,0,0,2,TN,tolsf)
              Rwork(14) = tolsf
              Istate = -2
              GOTO 2100
           ENDIF
        ENDIF
  !
  ! KGO = 1:success; 2:error test failure; 3:convergence failure;
  !       4:RES ordered return. 5:RES returned error.
  !-----------------------------------------------------------------------
  ! Block F.
  ! The following block handles the case of a successful return from the
  ! core integrator (KFLAG = 0).  Test for stop conditions.
  !-----------------------------------------------------------------------
   800  INIt = 1
        GOTO (900,1300,1000,1100,1200) , Itask
  ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   900  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        T = Tout
        GOTO 1400
  ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   1000 IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        GOTO 1300
  ! ITASK = 4.  see if TOUT or TCRIT was reached.  adjust h if necessary.
   1100 IF ( (TN-Tout)*H.LT.0.0D0 ) THEN
           hmx = ABS(TN) + ABS(H)
           ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
           IF ( ihit ) GOTO 1300
           tnext = TN + H*(1.0D0+4.0D0*UROund)
           IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
              H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
              JSTart = -2
           ENDIF
           GOTO 600
        ELSE
           CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
           T = Tout
           GOTO 1400
        ENDIF
  ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
   1200 hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
  !-----------------------------------------------------------------------
  ! Block G.
  ! The following block handles all successful returns from DLSODI.
  ! if ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
  ! ISTATE is set to 2, and the optional outputs are loaded into the
  ! work arrays before returning.
  !-----------------------------------------------------------------------
   1300 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
           IF ( ihit ) T = tcrit
        ENDIF
   1400 Istate = 2
        IF ( KFLag.EQ.-3 ) Istate = 3
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
  ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
   1500 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODI-  At T (=R1), EWT(I1) has become R2 .le. 0.'
        CALL XERRWD(msg,50,202,0,1,i,0,2,TN,ewti)
        Istate = -6
        GOTO 2100
  ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   1600 msg = 'DLSODI-  At T(=R1) and step size H(=R2), the error'
        CALL XERRWD(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      test failed repeatedly or with ABS(H) = HMIN'
        CALL XERRWD(msg,50,204,0,0,0,0,2,TN,H)
        Istate = -4
        GOTO 1900
  ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
   1700 msg = 'DLSODI-  At T (=R1) and step size H (=R2), the    '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      corrector convergence failed repeatedly     '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      or with ABS(H) = HMIN   '
        CALL XERRWD(msg,30,205,0,0,0,0,2,TN,H)
        Istate = -5
        GOTO 1900
  ! IRES = 3 returned by RES, despite retries by DSTODI. -----------------
   1800 msg = 'DLSODI-  At T (=R1) residual routine returned     '
        CALL XERRWD(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      error IRES = 3 repeatedly.        '
        CALL XERRWD(msg,40,206,0,0,0,0,1,TN,0.0D0)
        Istate = -7
        GOTO 2100
  ! Compute IMXER if relevant. -------------------------------------------
   1900 big = 0.0D0
        imxer = 1
        DO i = 1 , N
           size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
           IF ( big.LT.size ) THEN
              big = size
              imxer = i
           ENDIF
        ENDDO
        Iwork(16) = imxer
  ! Compute residual if relevant. ----------------------------------------
   2000 lyd0 = LYH + NYH
        DO i = 1 , N
           Rwork(i+LSAvf-1) = Rwork(i+lyd0-1)/H
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        ires = 1
        CALL RES(Neq,TN,Y,Rwork(LSAvf),Ydoti,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) THEN
           msg = 'DLSODI-  Residual routine set its flag IRES       '
           CALL XERRWD(msg,50,210,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      to (I1) when called for final output.       '
           CALL XERRWD(msg,50,210,0,1,ires,0,0,0.0D0,0.0D0)
        ENDIF
        GOTO 2200
  ! Set Y vector, T, and optional outputs. -------------------------------
   2100 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
   2200 T = TN
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
   2300 msg = 'DLSODI-  MF (=I1) illegal.    '
        CALL XERRWD(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
        GOTO 3000
   2400 msg = 'DLSODI-  RTOL(=I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
        GOTO 3000
   2500 msg = 'DLSODI-  ATOL(=I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,20,0,1,i,0,1,atoli,0.0D0)
        GOTO 3000
   2600 ewti = Rwork(LEWt+i-1)
        msg = 'DLSODI-  EWT(I1) is R1 .le. 0.0         '
        CALL XERRWD(msg,40,21,0,1,i,0,1,ewti,0.0D0)
        GOTO 3000
   2700 msg =                                                             &
       &    'DLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
        CALL XERRWD(msg,60,24,0,0,0,0,2,tcrit,TN)
        GOTO 3000
   2800 msg =                                                             &
       &    'DLSODI-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
        CALL XERRWD(msg,60,25,0,0,0,0,2,tcrit,Tout)
        GOTO 3000
   2900 msg = 'DLSODI-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
        CALL XERRWD(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
  !
   3000 Istate = -3
        RETURN
  !----------------------- End of Subroutine DLSODI ----------------------
  99999 END
  !*==DLSOIBT.spg  processed by SPAG 6.72Dc at 17:57 on 29 Mar 2017
  !DECK DLSOIBT
        SUBROUTINE DLSOIBT(RES,ADDA,JAC,Neq,Y,Ydoti,T,Tout,Itol,Rtol,Atol,&
       &                   Itask,Istate,Iopt,Rwork,Lrw,Iwork,Liw,Mf)
        IMPLICIT NONE
  !*--DLSOIBT7343
  !*** Start of declarations inserted by SPAG
        !REAL ADDA , DPJIBT , DSLSBT
        !INTEGER JAC
  !*** End of declarations inserted by SPAG
        EXTERNAL RES , ADDA , JAC
        INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
       &        Mf
        DOUBLE PRECISION Y , Ydoti , T , Tout , Rtol , Atol , Rwork
        DIMENSION Neq(*) , Y(*) , Ydoti(*) , Rtol(*) , Atol(*) ,          &
       &          Rwork(Lrw) , Iwork(Liw)
  !-----------------------------------------------------------------------
  ! This is the 18 November 2003 version of
  ! DLSOIBT: Livermore Solver for Ordinary differential equations given
  !          in Implicit form, with Block-Tridiagonal Jacobian treatment.
  !
  ! This version is in double precision.
  !
  ! DLSOIBT solves the initial value problem for linearly implicit
  ! systems of first order ODEs,
  !     A(t,y) * dy/dt = g(t,y) ,  where A(t,y) is a square matrix,
  ! or, in component form,
  !     ( a   * ( dy / dt ))  + ... +  ( a     * ( dy   / dt ))  =
  !        i,1      1                     i,NEQ      NEQ
  !
  !      =   g ( t, y , y ,..., y    )   ( i = 1,...,NEQ )
  !           i      1   2       NEQ
  !
  ! If A is singular, this is a differential-algebraic system.
  !
  ! DLSOIBT is a variant version of the DLSODI package, for the case where
  ! the matrices A, dg/dy, and d(A*s)/dy are all block-tridiagonal.
  !-----------------------------------------------------------------------
  ! Reference:
  !     Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
  !     Solvers, in Scientific Computing,  R. S. Stepleman et al. (Eds.),
  !     North-Holland, Amsterdam, 1983, pp. 55-64.
  !-----------------------------------------------------------------------
  ! Authors:       Alan C. Hindmarsh and Jeffrey F. Painter
  !                Center for Applied Scientific Computing, L-561
  !                Lawrence Livermore National Laboratory
  !                Livermore, CA 94551
  ! and
  !                Charles S. Kenney
  ! formerly at:   Naval Weapons Center
  !                China Lake, CA 93555
  !-----------------------------------------------------------------------
  ! Summary of Usage.
  !
  ! Communication between the user and the DLSOIBT package, for normal
  ! situations, is summarized here.  This summary describes only a subset
  ! of the full set of options available.  See the full description for
  ! details, including optional communication, nonstandard options,
  ! and instructions for special situations.  See also the example
  ! problem (with program and output) following this summary.
  !
  ! A. First, provide a subroutine of the form:
  !               SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
  !               DOUBLE PRECISION T, Y(*), S(*), R(*)
  ! which computes the residual function
  !     r = g(t,y)  -  A(t,y) * s ,
  ! as a function of t and the vectors y and s.  (s is an internally
  ! generated approximation to dy/dt.)  The arrays Y and S are inputs
  ! to the RES routine and should not be altered.  The residual
  ! vector is to be stored in the array R.  The argument IRES should be
  ! ignored for casual use of DLSOIBT.  (For uses of IRES, see the
  ! paragraph on RES in the full description below.)
  !
  ! B. Next, identify the block structure of the matrices A = A(t,y) and
  ! dr/dy.  DLSOIBT must deal internally with a linear combination, P, of
  ! these two matrices.  The matrix P (hence both A and dr/dy) must have
  ! a block-tridiagonal form with fixed structure parameters
  !     MB = block size, MB .ge. 1, and
  !     NB = number of blocks in each direction, NB .ge. 4,
  ! with MB*NB = NEQ.  In each of the NB block-rows of the matrix P
  ! (each consisting of MB consecutive rows), the nonzero elements are
  ! to lie in three consecutive MB by MB blocks.  In block-rows
  ! 2 through NB - 1, these are centered about the main diagonal.
  ! in block-rows 1 and NB, they are the diagonal blocks and the two
  ! blocks adjacent to the diagonal block.  (Thus block positions (1,3)
  ! and (NB,NB-2) can be nonzero.)
  ! Alternatively, P (hence A and dr/dy) may be only approximately
  ! equal to matrices with this form, and DLSOIBT should still succeed.
  ! The block-tridiagonal matrix P is described by three arrays,
  ! each of size MB by MB by NB:
  !     PA = array of diagonal blocks,
  !     PB = array of superdiagonal (and one subdiagonal) blocks, and
  !     PC = array of subdiagonal (and one superdiagonal) blocks.
  ! Specifically, the three MB by MB blocks in the k-th block-row of P
  ! are stored in (reading across):
  !     PC(*,*,k) = block to the left of the diagonal block,
  !     PA(*,*,k) = diagonal block, and
  !     PB(*,*,k) = block to the right of the diagonal block,
  ! except for k = 1, where the three blocks (reading across) are
  !     PA(*,*,1) (= diagonal block), PB(*,*,1), and PC(*,*,1),
  ! and k = NB, where they are
  !     PB(*,*,NB), PC(*,*,NB), and PA(*,*,NB) (= diagonal block).
  ! (Each asterisk * stands for an index that ranges from 1 to MB.)
  !
  ! C. You must also provide a subroutine of the form:
  !     SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
  !     DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
  ! which adds the nonzero blocks of the matrix A = A(t,y) to the
  ! contents of the arrays PA, PB, and PC, following the structure
  ! description in Paragraph B above.
  ! T and the Y array are input and should not be altered.
  ! Thus the affect of ADDA should be the following:
  !     DO 30 K = 1,NB
  !       DO 20 J = 1,MB
  !         DO 10 I = 1,MB
  !           PA(I,J,K) = PA(I,J,K) +
  !             ( (I,J) element of K-th diagonal block of A)
  !           PB(I,J,K) = PB(I,J,K) +
  !             ( (I,J) element of block in block position (K,K+1) of A,
  !             or in block position (NB,NB-2) if K = NB)
  !           PC(I,J,K) = PC(I,J,K) +
  !             ( (I,J) element of block in block position (K,K-1) of A,
  !             or in block position (1,3) if K = 1)
  ! 10        CONTINUE
  ! 20      CONTINUE
  ! 30    CONTINUE
  !
  ! D. For the sake of efficiency, you are encouraged to supply the
  ! Jacobian matrix dr/dy in closed form, where r = g(t,y) - A(t,y)*s
  ! (s = a fixed vector) as above.  If dr/dy is being supplied,
  ! use MF = 21, and provide a subroutine of the form:
  !     SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
  !     DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB), PB(MB,MB,NB),
  !    1                 PC(MB,MB,NB)
  ! which computes dr/dy as a function of t, y, and s.  Here T, Y, and
  ! S are inputs, and the routine is to load dr/dy into PA, PB, PC,
  ! according to the structure description in Paragraph B above.
  ! That is, load the diagonal blocks into PA, the superdiagonal blocks
  ! (and block (NB,NB-2) ) into PB, and the subdiagonal blocks (and
  ! block (1,3) ) into PC.  The blocks in block-row k of dr/dy are to
  ! be loaded into PA(*,*,k), PB(*,*,k), and PC(*,*,k).
  !     Only nonzero elements need be loaded, and the indexing
  ! of PA, PB, and PC is the same as in the ADDA routine.
  !     Note that if A is independent of Y (or this dependence
  ! is weak enough to be ignored) then JAC is to compute dg/dy.
  !     If it is not feasible to provide a JAC routine, use
  ! MF = 22, and DLSOIBT will compute an approximate Jacobian
  ! internally by difference quotients.
  !
  ! E. Next decide whether or not to provide the initial value of the
  ! derivative vector dy/dt.  If the initial value of A(t,y) is
  ! nonsingular (and not too ill-conditioned), you may let DLSOIBT compute
  ! this vector (ISTATE = 0).  (DLSOIBT will solve the system A*s = g for
  ! s, with initial values of A and g.)  If A(t,y) is initially
  ! singular, then the system is a differential-algebraic system, and
  ! you must make use of the particular form of the system to compute the
  ! initial values of y and dy/dt.  In that case, use ISTATE = 1 and
  ! load the initial value of dy/dt into the array YDOTI.
  ! The input array YDOTI and the initial Y array must be consistent with
  ! the equations A*dy/dt = g.  This implies that the initial residual
  ! r = g(t,y) - A(t,y)*YDOTI  must be approximately zero.
  !
  ! F. Write a main program which calls Subroutine DLSOIBT once for
  ! each point at which answers are desired.  This should also provide
  ! for possible use of logical unit 6 for output of error messages by
  ! DLSOIBT.  on the first call to DLSOIBT, supply arguments as follows:
  ! RES    = name of user subroutine for residual function r.
  ! ADDA   = name of user subroutine for computing and adding A(t,y).
  ! JAC    = name of user subroutine for Jacobian matrix dr/dy
  !          (MF = 21).  If not used, pass a dummy name.
  ! Note: the names for the RES and ADDA routines and (if used) the
  !        JAC routine must be declared External in the calling program.
  ! NEQ    = number of scalar equations in the system.
  ! Y      = array of initial values, of length NEQ.
  ! YDOTI  = array of length NEQ (containing initial dy/dt if ISTATE = 1).
  ! T      = the initial value of the independent variable.
  ! TOUT   = first point where output is desired (.ne. T).
  ! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
  ! RTOL   = relative tolerance parameter (scalar).
  ! ATOL   = absolute tolerance parameter (scalar or array).
  !          the estimated local error in y(i) will be controlled so as
  !          to be roughly less (in magnitude) than
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
  !             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
  !          Thus the local error test passes if, in each component,
  !          either the absolute error is less than ATOL (or ATOL(i)),
  !          or the relative error is less than RTOL.
  !          Use RTOL = 0.0 for pure absolute error control, and
  !          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
  !          control.  Caution: Actual (global) errors may exceed these
  !          local tolerances, so choose them conservatively.
  ! ITASK  = 1 for normal computation of output values of y at t = TOUT.
  ! ISTATE = integer flag (input and output).  Set ISTATE = 1 if the
  !          initial dy/dt is supplied, and 0 otherwise.
  ! IOPT   = 0 to indicate no optional inputs used.
  ! RWORK  = real work array of length at least:
  !             22 + 9*NEQ + 3*MB*MB*NB        for MF = 21 or 22.
  ! LRW    = declared length of RWORK (in user's dimension).
  ! IWORK  = integer work array of length at least 20 + NEQ.
  !          Input in IWORK(1) the block size MB and in IWORK(2) the
  !          number NB of blocks in each direction along the matrix A.
  !          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
  ! LIW    = declared length of IWORK (in user's dimension).
  ! MF     = method flag.  Standard values are:
  !          21 for a user-supplied Jacobian.
  !          22 for an internally generated Jacobian.
  !          For other choices of MF, see the paragraph on MF in
  !          the full description below.
  ! Note that the main program must declare arrays Y, YDOTI, RWORK, IWORK,
  ! and possibly ATOL.
  !
  ! G. The output from the first call (or any call) is:
  !      Y = array of computed values of y(t) vector.
  !      T = corresponding value of independent variable (normally TOUT).
  ! ISTATE = 2  if DLSOIBT was successful, negative otherwise.
  !          -1 means excess work done on this call (check all inputs).
  !          -2 means excess accuracy requested (tolerances too small).
  !          -3 means illegal input detected (see printed message).
  !          -4 means repeated error test failures (check all inputs).
  !          -5 means repeated convergence failures (perhaps bad Jacobian
  !             supplied or wrong choice of tolerances).
  !          -6 means error weight became zero during problem. (Solution
  !             component i vanished, and ATOL or ATOL(i) = 0.)
  !          -7 cannot occur in casual use.
  !          -8 means DLSOIBT was unable to compute the initial dy/dt.
  !             In casual use, this means A(t,y) is initially singular.
  !             Supply YDOTI and use ISTATE = 1 on the first call.
  !
  !  If DLSOIBT returns ISTATE = -1, -4, or -5, then the output of
  !  DLSOIBT also includes YDOTI = array containing residual vector
  !  r = g - A * dy/dt  evaluated at the current t, y, and dy/dt.
  !
  ! H. To continue the integration after a successful return, simply
  ! reset TOUT and call DLSOIBT again.  No other parameters need be reset.
  !
  !-----------------------------------------------------------------------
  ! Example Problem.
  !
  ! The following is an example problem, with the coding needed
  ! for its solution by DLSOIBT.  The problem comes from the partial
  ! differential equation (the Burgers equation)
  !   du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
  ! on -1 .le. x .le. 1.  The boundary conditions are
  !   du/dx = 0  at x = -1 and at x = 1.
  ! The initial profile is a square wave,
  !   u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
  ! The PDE is discretized in x by a simplified Galerkin method,
  ! using piecewise linear basis functions, on a grid of 40 intervals.
  ! The equations at x = -1 and 1 use a 3-point difference approximation
  ! for the right-hand side.  The result is a system A * dy/dt = g(y),
  ! of size NEQ = 41, where y(i) is the approximation to u at x = x(i),
  ! with x(i) = -1 + (i-1)*delx, delx = 2/(NEQ-1) = .05.  The individual
  ! equations in the system are
  !   dy(1)/dt = ( y(3) - 2*y(2) + y(1) ) * eta / delx**2,
  !   dy(NEQ)/dt = ( y(NEQ-2) - 2*y(NEQ-1) + y(NEQ) ) * eta / delx**2,
  ! and for i = 2, 3, ..., NEQ-1,
  !   (1/6) dy(i-1)/dt + (4/6) dy(i)/dt + (1/6) dy(i+1)/dt
  !       = ( y(i-1)**2 - y(i+1)**2 ) / (4*delx)
  !         + ( y(i+1) - 2*y(i) + y(i-1) ) * eta / delx**2.
  ! The following coding solves the problem with MF = 21, with output
  ! of solution statistics at t = .1, .2, .3, and .4, and of the
  ! solution vector at t = .4.  Here the block size is just MB = 1.
  !
  !     EXTERNAL RESID, ADDABT, JACBT
  !     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y, YDOTI
  !     DIMENSION Y(41), YDOTI(41), RWORK(514), IWORK(61)
  !     NEQ = 41
  !     DO 10 I = 1,NEQ
  !  10   Y(I) = 0.0
  !     Y(11) = 0.5
  !     DO 20 I = 12,30
  !  20   Y(I) = 1.0
  !     Y(31) = 0.5
  !     T = 0.0
  !     TOUT = 0.1
  !     ITOL = 1
  !     RTOL = 1.0D-4
  !     ATOL = 1.0D-5
  !     ITASK = 1
  !     ISTATE = 0
  !     IOPT = 0
  !     LRW = 514
  !     LIW = 61
  !     IWORK(1) = 1
  !     IWORK(2) = NEQ
  !     MF = 21
  !     DO 40 IO = 1,4
  !       CALL DLSOIBT (RESID, ADDABT, JACBT, NEQ, Y, YDOTI, T, TOUT,
  !    1     ITOL,RTOL,ATOL, ITASK, ISTATE, IOPT, RWORK,LRW,IWORK,LIW, MF)
  !       WRITE (6,30) T, IWORK(11), IWORK(12), IWORK(13)
  !  30   FORMAT(' At t =',F5.2,'   No. steps =',I4,'  No. r-s =',I4,
  !    1         '  No. J-s =',I3)
  !       IF (ISTATE .NE. 2)  GO TO 90
  !       TOUT = TOUT + 0.1
  !  40   CONTINUE
  !     WRITE(6,50) (Y(I),I=1,NEQ)
  !  50 FORMAT(/' Final solution values..'/9(5D12.4/))
  !     STOP
  !  90 WRITE(6,95) ISTATE
  !  95 FORMAT(///' Error halt.. ISTATE =',I3)
  !     STOP
  !     END
  !
  !     SUBROUTINE RESID (N, T, Y, S, R, IRES)
  !     DOUBLE PRECISION T, Y, S, R, ETA, DELX, EODSQ
  !     DIMENSION Y(N), S(N), R(N)
  !     DATA ETA/0.05/, DELX/0.05/
  !     EODSQ = ETA/DELX**2
  !     R(1) = EODSQ*(Y(3) - 2.0*Y(2) + Y(1)) - S(1)
  !     NM1 = N - 1
  !     DO 10 I = 2,NM1
  !       R(I) = (Y(I-1)**2 - Y(I+1)**2)/(4.0*DELX)
  !    1        + EODSQ*(Y(I+1) - 2.0*Y(I) + Y(I-1))
  !    2        - (S(I-1) + 4.0*S(I) + S(I+1))/6.0
  !  10   CONTINUE
  !     R(N) = EODSQ*(Y(N-2) - 2.0*Y(NM1) + Y(N)) - S(N)
  !     RETURN
  !     END
  !
  !     SUBROUTINE ADDABT (N, T, Y, MB, NB, PA, PB, PC)
  !     DOUBLE PRECISION T, Y, PA, PB, PC
  !     DIMENSION Y(N), PA(MB,MB,NB), PB(MB,MB,NB), PC(MB,MB,NB)
  !     PA(1,1,1) = PA(1,1,1) + 1.0
  !     NM1 = N - 1
  !     DO 10 K = 2,NM1
  !       PA(1,1,K) = PA(1,1,K) + (4.0/6.0)
  !       PB(1,1,K) = PB(1,1,K) + (1.0/6.0)
  !       PC(1,1,K) = PC(1,1,K) + (1.0/6.0)
  !  10   CONTINUE
  !     PA(1,1,N) = PA(1,1,N) + 1.0
  !     RETURN
  !     END
  !
  !     SUBROUTINE JACBT (N, T, Y, S, MB, NB, PA, PB, PC)
  !     DOUBLE PRECISION T, Y, S, PA, PB, PC, ETA, DELX, EODSQ
  !     DIMENSION Y(N), S(N), PA(MB,MB,NB),PB(MB,MB,NB),PC(MB,MB,NB)
  !     DATA ETA/0.05/, DELX/0.05/
  !     EODSQ = ETA/DELX**2
  !     PA(1,1,1) = EODSQ
  !     PB(1,1,1) = -2.0*EODSQ
  !     PC(1,1,1) = EODSQ
  !     DO 10 K = 2,N
  !       PA(1,1,K) = -2.0*EODSQ
  !       PB(1,1,K) = -Y(K+1)*(0.5/DELX) + EODSQ
  !       PC(1,1,K) = Y(K-1)*(0.5/DELX) + EODSQ
  !  10   CONTINUE
  !     PB(1,1,N) = EODSQ
  !     PC(1,1,N) = -2.0*EODSQ
  !     PA(1,1,N) = EODSQ
  !     RETURN
  !     END
  !
  ! The output of this program (on a CDC-7600 in single precision)
  ! is as follows:
  !
  ! At t = 0.10   No. steps =  35  No. r-s =  45  No. J-s =  9
  ! At t = 0.20   No. steps =  43  No. r-s =  54  No. J-s = 10
  ! At t = 0.30   No. steps =  48  No. r-s =  60  No. J-s = 11
  ! At t = 0.40   No. steps =  51  No. r-s =  64  No. J-s = 12
  !
  ! Final solution values..
  !  1.2747e-02  1.1997e-02  1.5560e-02  2.3767e-02  3.7224e-02
  !  5.6646e-02  8.2645e-02  1.1557e-01  1.5541e-01  2.0177e-01
  !  2.5397e-01  3.1104e-01  3.7189e-01  4.3530e-01  5.0000e-01
  !  5.6472e-01  6.2816e-01  6.8903e-01  7.4612e-01  7.9829e-01
  !  8.4460e-01  8.8438e-01  9.1727e-01  9.4330e-01  9.6281e-01
  !  9.7632e-01  9.8426e-01  9.8648e-01  9.8162e-01  9.6617e-01
  !  9.3374e-01  8.7535e-01  7.8236e-01  6.5321e-01  5.0003e-01
  !  3.4709e-01  2.1876e-01  1.2771e-01  7.3671e-02  5.0642e-02
  !  5.4496e-02
  !
  !-----------------------------------------------------------------------
  ! Full Description of User Interface to DLSOIBT.
  !
  ! The user interface to DLSOIBT consists of the following parts.
  !
  ! 1.   The call sequence to Subroutine DLSOIBT, which is a driver
  !      routine for the solver.  This includes descriptions of both
  !      the call sequence arguments and of user-supplied routines.
  !      Following these descriptions is a description of
  !      optional inputs available through the call sequence, and then
  !      a description of optional outputs (in the work arrays).
  !
  ! 2.   Descriptions of other routines in the DLSOIBT package that may be
  !      (optionally) called by the user.  These provide the ability to
  !      alter error message handling, save and restore the internal
  !      Common, and obtain specified derivatives of the solution y(t).
  !
  ! 3.   Descriptions of Common blocks to be declared in overlay
  !      or similar environments, or to be saved when doing an interrupt
  !      of the problem and continued solution later.
  !
  ! 4.   Description of two routines in the DLSOIBT package, either of
  !      which the user may replace with his/her own version, if desired.
  !      These relate to the measurement of errors.
  !
  !-----------------------------------------------------------------------
  ! Part 1.  Call Sequence.
  !
  ! The call sequence parameters used for input only are
  !     RES, ADDA, JAC, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK,
  !     IOPT, LRW, LIW, MF,
  ! and those used for both input and output are
  !     Y, T, ISTATE, YDOTI.
  ! The work arrays RWORK and IWORK are also used for additional and
  ! optional inputs and optional outputs.  (The term output here refers
  ! to the return from Subroutine DLSOIBT to the user's calling program.)
  !
  ! The legality of input parameters will be thoroughly checked on the
  ! initial call for the problem, but not checked thereafter unless a
  ! change in input parameters is flagged by ISTATE = 3 on input.
  !
  ! The descriptions of the call arguments are as follows.
  !
  ! RES    = the name of the user-supplied subroutine which supplies
  !          the residual vector for the ODE system, defined by
  !            r = g(t,y) - A(t,y) * s
  !          as a function of the scalar t and the vectors
  !          s and y (s approximates dy/dt).  This subroutine
  !          is to have the form
  !              SUBROUTINE RES (NEQ, T, Y, S, R, IRES)
  !              DOUBLE PRECISION T, Y(*), S(*), R(*)
  !          where NEQ, T, Y, S, and IRES are input, and R and
  !          IRES are output. Y, S, and R are arrays of length NEQ.
  !             On input, IRES indicates how DLSOIBT will use the
  !          returned array R, as follows:
  !             IRES = 1  means that DLSOIBT needs the full residual,
  !                       r = g - A*s, exactly.
  !             IRES = -1 means that DLSOIBT is using R only to compute
  !                       the Jacobian dr/dy by difference quotients.
  !          The RES routine can ignore IRES, or it can omit some terms
  !          if IRES = -1.  If A does not depend on y, then RES can
  !          just return R = g when IRES = -1.  If g - A*s contains other
  !          additive terms that are independent of y, these can also be
  !          dropped, if done consistently, when IRES = -1.
  !             The subroutine should set the flag IRES if it
  !          encounters a halt condition or illegal input.
  !          Otherwise, it should not reset IRES.  On output,
  !             IRES = 1 or -1 represents a normal return, and
  !          DLSOIBT continues integrating the ODE.  Leave IRES
  !          unchanged from its input value.
  !             IRES = 2 tells DLSOIBT to immediately return control
  !          to the calling program, with ISTATE = 3.  This lets
  !          the calling program change parameters of the problem
  !          if necessary.
  !             IRES = 3 represents an error condition (for example, an
  !          illegal value of y).  DLSOIBT tries to integrate the system
  !          without getting IRES = 3 from RES.  If it cannot, DLSOIBT
  !          returns with ISTATE = -7 or -1.
  !             On an DLSOIBT return with ISTATE = 3, -1, or -7, the
  !          values of T and Y returned correspond to the last point
  !          reached successfully without getting the flag IRES = 2 or 3.
  !             The flag values IRES = 2 and 3 should not be used to
  !          handle switches or root-stop conditions.  This is better
  !          done by calling DLSOIBT in a one-step mode and checking the
  !          stopping function for a sign change at each step.
  !             If quantities computed in the RES routine are needed
  !          externally to DLSOIBT, an extra call to RES should be made
  !          for this purpose, for consistent and accurate results.
  !          To get the current dy/dt for the S argument, use DINTDY.
  !             RES must be declared External in the calling
  !          program. See note below for more about RES.
  !
  ! ADDA   = the name of the user-supplied subroutine which adds the
  !          matrix A = A(t,y) to another matrix, P, stored in
  !          block-tridiagonal form.  This routine is to have the form
  !               SUBROUTINE ADDA (NEQ, T, Y, MB, NB, PA, PB, PC)
  !               DOUBLE PRECISION T, Y(*), PA(MB,MB,NB), PB(MB,MB,NB),
  !              1                 PC(MB,MB,NB)
  !          where NEQ, T, Y, MB, NB, and the arrays PA, PB, and PC
  !          are input, and the arrays PA, PB, and PC are output.
  !          Y is an array of length NEQ, and the arrays PA, PB, PC
  !          are all MB by MB by NB.
  !             Here a block-tridiagonal structure is assumed for A(t,y),
  !          and also for the matrix P to which A is added here,
  !          as described in Paragraph B of the Summary of Usage above.
  !          Thus the affect of ADDA should be the following:
  !               DO 30 K = 1,NB
  !                 DO 20 J = 1,MB
  !                   DO 10 I = 1,MB
  !                     PA(I,J,K) = PA(I,J,K) +
  !                       ( (I,J) element of K-th diagonal block of A)
  !                     PB(I,J,K) = PB(I,J,K) +
  !                       ( (I,J) element of block (K,K+1) of A,
  !                       or block (NB,NB-2) if K = NB)
  !                     PC(I,J,K) = PC(I,J,K) +
  !                       ( (I,J) element of block (K,K-1) of A,
  !                       or block (1,3) if K = 1)
  !           10        CONTINUE
  !           20      CONTINUE
  !           30    CONTINUE
  !             ADDA must be declared External in the calling program.
  !          See note below for more information about ADDA.
  !
  ! JAC    = the name of the user-supplied subroutine which supplies
  !          the Jacobian matrix, dr/dy, where r = g - A*s.  JAC is
  !          required if MITER = 1.  Otherwise a dummy name can be
  !          passed.  This subroutine is to have the form
  !               SUBROUTINE JAC (NEQ, T, Y, S, MB, NB, PA, PB, PC)
  !               DOUBLE PRECISION T, Y(*), S(*), PA(MB,MB,NB),
  !              1                 PB(MB,MB,NB), PC(MB,MB,NB)
  !          where NEQ, T, Y, S, MB, NB, and the arrays PA, PB, and PC
  !          are input, and the arrays PA, PB, and PC are output.
  !          Y and S are arrays of length NEQ, and the arrays PA, PB, PC
  !          are all MB by MB by NB.
  !          PA, PB, and PC are to be loaded with partial derivatives
  !          (elements of the Jacobian matrix) on output, in terms of the
  !          block-tridiagonal structure assumed, as described
  !          in Paragraph B of the Summary of Usage above.
  !          That is, load the diagonal blocks into PA, the
  !          superdiagonal blocks (and block (NB,NB-2) ) into PB, and
  !          the subdiagonal blocks (and block (1,3) ) into PC.
  !          The blocks in block-row k of dr/dy are to be loaded into
  !          PA(*,*,k), PB(*,*,k), and PC(*,*,k).
  !          Thus the affect of JAC should be the following:
  !               DO 30 K = 1,NB
  !                 DO 20 J = 1,MB
  !                   DO 10 I = 1,MB
  !                     PA(I,J,K) = ( (I,J) element of
  !                       K-th diagonal block of dr/dy)
  !                     PB(I,J,K) = ( (I,J) element of block (K,K+1)
  !                       of dr/dy, or block (NB,NB-2) if K = NB)
  !                     PC(I,J,K) = ( (I,J) element of block (K,K-1)
  !                       of dr/dy, or block (1,3) if K = 1)
  !           10        CONTINUE
  !           20      CONTINUE
  !           30    CONTINUE
  !               PA, PB, and PC are preset to zero by the solver,
  !          so that only the nonzero elements need be loaded by JAC.
  !          Each call to JAC is preceded by a call to RES with the same
  !          arguments NEQ, T, Y, and S.  Thus to gain some efficiency,
  !          intermediate quantities shared by both calculations may be
  !          saved in a user Common block by RES and not recomputed by JAC
  !          if desired.  Also, JAC may alter the Y array, if desired.
  !               JAC need not provide dr/dy exactly.  A crude
  !          approximation will do, so that DLSOIBT may be used when
  !          A and dr/dy are not really block-tridiagonal, but are close
  !          to matrices that are.
  !               JAC must be declared External in the calling program.
  !               See note below for more about JAC.
  !
  !    Note on RES, ADDA, and JAC:
  !          These subroutines may access user-defined quantities in
  !          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
  !          (dimensioned in the subroutines) and/or Y has length
  !          exceeding NEQ(1).  However, these routines should not alter
  !          NEQ(1), Y(1),...,Y(NEQ) or any other input variables.
  !          See the descriptions of NEQ and Y below.
  !
  ! NEQ    = the size of the system (number of first order ordinary
  !          differential equations or scalar algebraic equations).
  !          Used only for input.
  !          NEQ may be decreased, but not increased, during the problem.
  !          If NEQ is decreased (with ISTATE = 3 on input), the
  !          remaining components of Y should be left undisturbed, if
  !          these are to be accessed in RES, ADDA, or JAC.
  !
  !          Normally, NEQ is a scalar, and it is generally referred to
  !          as a scalar in this user interface description.  However,
  !          NEQ may be an array, with NEQ(1) set to the system size.
  !          (The DLSOIBT package accesses only NEQ(1).)  In either case,
  !          this parameter is passed as the NEQ argument in all calls
  !          to RES, ADDA, and JAC.  Hence, if it is an array,
  !          locations NEQ(2),... may be used to store other integer data
  !          and pass it to RES, ADDA, or JAC.  Each such subroutine
  !          must include NEQ in a Dimension statement in that case.
  !
  ! Y      = a real array for the vector of dependent variables, of
  !          length NEQ or more.  Used for both input and output on the
  !          first call (ISTATE = 0 or 1), and only for output on other
  !          calls.  On the first call, Y must contain the vector of
  !          initial values.  On output, Y contains the computed solution
  !          vector, evaluated at t.  If desired, the Y array may be used
  !          for other purposes between calls to the solver.
  !
  !          This array is passed as the Y argument in all calls to RES,
  !          ADDA, and JAC.  Hence its length may exceed NEQ,
  !          and locations Y(NEQ+1),... may be used to store other real
  !          data and pass it to RES, ADDA, or JAC.  (The DLSOIBT
  !          package accesses only Y(1),...,Y(NEQ). )
  !
  ! YDOTI  = a real array for the initial value of the vector
  !          dy/dt and for work space, of dimension at least NEQ.
  !
  !          On input:
  !            If ISTATE = 0 then DLSOIBT will compute the initial value
  !          of dy/dt, if A is nonsingular.  Thus YDOTI will
  !          serve only as work space and may have any value.
  !            If ISTATE = 1 then YDOTI must contain the initial value
  !          of dy/dt.
  !            If ISTATE = 2 or 3 (continuation calls) then YDOTI
  !          may have any value.
  !            Note: If the initial value of A is singular, then
  !          DLSOIBT cannot compute the initial value of dy/dt, so
  !          it must be provided in YDOTI, with ISTATE = 1.
  !
  !          On output, when DLSOIBT terminates abnormally with ISTATE =
  !          -1, -4, or -5, YDOTI will contain the residual
  !          r = g(t,y) - A(t,y)*(dy/dt).  If r is large, t is near
  !          its initial value, and YDOTI is supplied with ISTATE = 1,
  !          there may have been an incorrect input value of
  !          YDOTI = dy/dt, or the problem (as given to DLSOIBT)
  !          may not have a solution.
  !
  !          If desired, the YDOTI array may be used for other
  !          purposes between calls to the solver.
  !
  ! T      = the independent variable.  On input, T is used only on the
  !          first call, as the initial point of the integration.
  !          On output, after each call, T is the value at which a
  !          computed solution y is evaluated (usually the same as TOUT).
  !          On an error return, T is the farthest point reached.
  !
  ! TOUT   = the next value of t at which a computed solution is desired.
  !          Used only for input.
  !
  !          When starting the problem (ISTATE = 0 or 1), TOUT may be
  !          equal to T for one call, then should .ne. T for the next
  !          call.  For the initial T, an input value of TOUT .ne. T is
  !          used in order to determine the direction of the integration
  !          (i.e. the algebraic sign of the step sizes) and the rough
  !          scale of the problem.  Integration in either direction
  !          (forward or backward in t) is permitted.
  !
  !          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
  !          the first call (i.e. the first call with TOUT .ne. T).
  !          Otherwise, TOUT is required on every call.
  !
  !          If ITASK = 1, 3, or 4, the values of TOUT need not be
  !          monotone, but a value of TOUT which backs up is limited
  !          to the current internal T interval, whose endpoints are
  !          TCUR - HU and TCUR (see optional outputs, below, for
  !          TCUR and HU).
  !
  ! ITOL   = an indicator for the type of error control.  See
  !          description below under ATOL.  Used only for input.
  !
  ! RTOL   = a relative error tolerance parameter, either a scalar or
  !          an array of length NEQ.  See description below under ATOL.
  !          Input only.
  !
  ! ATOL   = an absolute error tolerance parameter, either a scalar or
  !          an array of length NEQ.  Input only.
  !
  !             The input parameters ITOL, RTOL, and ATOL determine
  !          the error control performed by the solver.  The solver will
  !          control the vector E = (E(i)) of estimated local errors
  !          in y, according to an inequality of the form
  !                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
  !          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
  !          and the RMS-norm (root-mean-square norm) here is
  !          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
  !          is a vector of weights which must always be positive, and
  !          the values of RTOL and ATOL should all be non-negative.
  !          The following table gives the types (scalar/array) of
  !          RTOL and ATOL, and the corresponding form of EWT(i).
  !
  !             ITOL    RTOL       ATOL          EWT(i)
  !              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  !              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  !              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  !              4     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL(i)
  !
  !          When either of these parameters is a scalar, it need not
  !          be dimensioned in the user's calling program.
  !
  !          If none of the above choices (with ITOL, RTOL, and ATOL
  !          fixed throughout the problem) is suitable, more general
  !          error controls can be obtained by substituting
  !          user-supplied routines for the setting of EWT and/or for
  !          the norm calculation.  See Part 4 below.
  !
  !          If global errors are to be estimated by making a repeated
  !          run on the same problem with smaller tolerances, then all
  !          components of RTOL and ATOL (i.e. of EWT) should be scaled
  !          down uniformly.
  !
  ! ITASK  = an index specifying the task to be performed.
  !          Input only.  ITASK has the following values and meanings.
  !          1  means normal computation of output values of y(t) at
  !             t = TOUT (by overshooting and interpolating).
  !          2  means take one step only and return.
  !          3  means stop at the first internal mesh point at or
  !             beyond t = TOUT and return.
  !          4  means normal computation of output values of y(t) at
  !             t = TOUT but without overshooting t = TCRIT.
  !             TCRIT must be input as RWORK(1).  TCRIT may be equal to
  !             or beyond TOUT, but not behind it in the direction of
  !             integration.  This option is useful if the problem
  !             has a singularity at or beyond t = TCRIT.
  !          5  means take one step, without passing TCRIT, and return.
  !             TCRIT must be input as RWORK(1).
  !
  !          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
  !          (within roundoff), it will return T = TCRIT (exactly) to
  !          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
  !          in which case answers at t = TOUT are returned first).
  !
  ! ISTATE = an index used for input and output to specify the
  !          state of the calculation.
  !
  !          On input, the values of ISTATE are as follows.
  !          0  means this is the first call for the problem, and
  !             DLSOIBT is to compute the initial value of dy/dt
  !             (while doing other initializations).  See note below.
  !          1  means this is the first call for the problem, and
  !             the initial value of dy/dt has been supplied in
  !             YDOTI (DLSOIBT will do other initializations).
  !             See note below.
  !          2  means this is not the first call, and the calculation
  !             is to continue normally, with no change in any input
  !             parameters except possibly TOUT and ITASK.
  !             (If ITOL, RTOL, and/or ATOL are changed between calls
  !             with ISTATE = 2, the new values will be used but not
  !             tested for legality.)
  !          3  means this is not the first call, and the
  !             calculation is to continue normally, but with
  !             a change in input parameters other than
  !             TOUT and ITASK.  Changes are allowed in
  !             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, MB, NB,
  !             and any of the optional inputs except H0.
  !             (See IWORK description for MB and NB.)
  !          Note:  A preliminary call with TOUT = T is not counted
  !          as a first call here, as no initialization or checking of
  !          input is done.  (Such a call is sometimes useful for the
  !          purpose of outputting the initial conditions.)
  !          Thus the first call for which TOUT .ne. T requires
  !          ISTATE = 0 or 1 on input.
  !
  !          On output, ISTATE has the following values and meanings.
  !           0 or 1  means nothing was done; TOUT = t and
  !              ISTATE = 0 or 1 on input.
  !           2  means that the integration was performed successfully.
  !           3  means that the user-supplied Subroutine RES signalled
  !              DLSOIBT to halt the integration and return (IRES = 2).
  !              Integration as far as T was achieved with no occurrence
  !              of IRES = 2, but this flag was set on attempting the
  !              next step.
  !          -1  means an excessive amount of work (more than MXSTEP
  !              steps) was done on this call, before completing the
  !              requested task, but the integration was otherwise
  !              successful as far as T.  (MXSTEP is an optional input
  !              and is normally 500.)  To continue, the user may
  !              simply reset ISTATE to a value .gt. 1 and call again
  !              (the excess work step counter will be reset to 0).
  !              In addition, the user may increase MXSTEP to avoid
  !              this error return (see below on optional inputs).
  !          -2  means too much accuracy was requested for the precision
  !              of the machine being used.  This was detected before
  !              completing the requested task, but the integration
  !              was successful as far as T.  To continue, the tolerance
  !              parameters must be reset, and ISTATE must be set
  !              to 3.  The optional output TOLSF may be used for this
  !              purpose.  (Note: If this condition is detected before
  !              taking any steps, then an illegal input return
  !              (ISTATE = -3) occurs instead.)
  !          -3  means illegal input was detected, before taking any
  !              integration steps.  See written message for details.
  !              Note:  If the solver detects an infinite loop of calls
  !              to the solver with illegal input, it will cause
  !              the run to stop.
  !          -4  means there were repeated error test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              The problem may have a singularity, or the input
  !              may be inappropriate.
  !          -5  means there were repeated convergence test failures on
  !              one attempted step, before completing the requested
  !              task, but the integration was successful as far as T.
  !              This may be caused by an inaccurate Jacobian matrix.
  !          -6  means EWT(i) became zero for some i during the
  !              integration.  Pure relative error control (ATOL(i) = 0.0)
  !              was requested on a variable which has now vanished.
  !              The integration was successful as far as T.
  !          -7  means that the user-supplied Subroutine RES set
  !              its error flag (IRES = 3) despite repeated tries by
  !              DLSOIBT to avoid that condition.
  !          -8  means that ISTATE was 0 on input but DLSOIBT was unable
  !              to compute the initial value of dy/dt.  See the
  !              printed message for details.
  !
  !          Note:  Since the normal output value of ISTATE is 2,
  !          it does not need to be reset for normal continuation.
  !          Similarly, ISTATE (= 3) need not be reset if RES told
  !          DLSOIBT to return because the calling program must change
  !          the parameters of the problem.
  !          Also, since a negative input value of ISTATE will be
  !          regarded as illegal, a negative output value requires the
  !          user to change it, and possibly other inputs, before
  !          calling the solver again.
  !
  ! IOPT   = an integer flag to specify whether or not any optional
  !          inputs are being used on this call.  Input only.
  !          The optional inputs are listed separately below.
  !          IOPT = 0 means no optional inputs are being used.
  !                   Default values will be used in all cases.
  !          IOPT = 1 means one or more optional inputs are being used.
  !
  ! RWORK  = a real working array (double precision).
  !          The length of RWORK must be at least
  !             20 + NYH*(MAXORD + 1) + 3*NEQ + LENWM    where
  !          NYH    = the initial value of NEQ,
  !          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
  !                   smaller value is given as an optional input),
  !          LENWM  = 3*MB*MB*NB + 2.
  !          (See MF description for the definition of METH.)
  !          Thus if MAXORD has its default value and NEQ is constant,
  !          this length is
  !             22 + 16*NEQ + 3*MB*MB*NB     for MF = 11 or 12,
  !             22 + 9*NEQ + 3*MB*MB*NB      for MF = 21 or 22.
  !          The first 20 words of RWORK are reserved for conditional
  !          and optional inputs and optional outputs.
  !
  !          The following word in RWORK is a conditional input:
  !            RWORK(1) = TCRIT = critical value of t which the solver
  !                       is not to overshoot.  Required if ITASK is
  !                       4 or 5, and ignored otherwise.  (See ITASK.)
  !
  ! LRW    = the length of the array RWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! IWORK  = an integer work array.  The length of IWORK must be at least
  !          20 + NEQ .  The first few words of IWORK are used for
  !          additional and optional inputs and optional outputs.
  !
  !          The following 2 words in IWORK are additional required
  !          inputs to DLSOIBT:
  !            IWORK(1) = MB = block size
  !            IWORK(2) = NB = number of blocks in the main diagonal
  !          These must satisfy  MB .ge. 1, NB .ge. 4, and MB*NB = NEQ.
  !
  ! LIW    = the length of the array IWORK, as declared by the user.
  !          (This will be checked by the solver.)
  !
  ! Note:  The work arrays must not be altered between calls to DLSOIBT
  ! for the same problem, except possibly for the additional and
  ! optional inputs, and except for the last 3*NEQ words of RWORK.
  ! The latter space is used for internal scratch space, and so is
  ! available for use by the user outside DLSOIBT between calls, if
  ! desired (but not for use by RES, ADDA, or JAC).
  !
  ! MF     = the method flag.  used only for input.  The legal values of
  !          MF are 11, 12, 21, and 22.
  !          MF has decimal digits METH and MITER: MF = 10*METH + MITER.
  !            METH indicates the basic linear multistep method:
  !              METH = 1 means the implicit Adams method.
  !              METH = 2 means the method based on Backward
  !                       Differentiation Formulas (BDFS).
  !                The BDF method is strongly preferred for stiff
  !              problems, while the Adams method is preferred when the
  !              problem is not stiff.  If the matrix A(t,y) is
  !              nonsingular, stiffness here can be taken to mean that of
  !              the explicit ODE system dy/dt = A-inverse * g.  If A is
  !              singular, the concept of stiffness is not well defined.
  !                If you do not know whether the problem is stiff, we
  !              recommend using METH = 2.  If it is stiff, the advantage
  !              of METH = 2 over METH = 1 will be great, while if it is
  !              not stiff, the advantage of METH = 1 will be slight.
  !              If maximum efficiency is important, some experimentation
  !              with METH may be necessary.
  !            MITER indicates the corrector iteration method:
  !              MITER = 1 means chord iteration with a user-supplied
  !                        block-tridiagonal Jacobian.
  !              MITER = 2 means chord iteration with an internally
  !                        generated (difference quotient) block-
  !                        tridiagonal Jacobian approximation, using
  !                        3*MB+1 extra calls to RES per dr/dy evaluation.
  !              If MITER = 1, the user must supply a Subroutine JAC
  !              (the name is arbitrary) as described above under JAC.
  !              For MITER = 2, a dummy argument can be used.
  !-----------------------------------------------------------------------
  ! Optional Inputs.
  !
  ! The following is a list of the optional inputs provided for in the
  ! call sequence.  (See also Part 2.)  For each such input variable,
  ! this table lists its name as used in this documentation, its
  ! location in the call sequence, its meaning, and the default value.
  ! The use of any of these inputs requires IOPT = 1, and in that
  ! case all of these inputs are examined.  A value of zero for any
  ! of these optional inputs will cause the default value to be used.
  ! Thus to use a subset of the optional inputs, simply preload
  ! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
  ! then set those of interest to nonzero values.
  !
  ! Name    Location      Meaning and Default Value
  !
  ! H0      RWORK(5)  the step size to be attempted on the first step.
  !                   The default value is determined by the solver.
  !
  ! HMAX    RWORK(6)  the maximum absolute step size allowed.
  !                   The default value is infinite.
  !
  ! HMIN    RWORK(7)  the minimum absolute step size allowed.
  !                   The default value is 0.  (This lower bound is not
  !                   enforced on the final step before reaching TCRIT
  !                   when ITASK = 4 or 5.)
  !
  ! MAXORD  IWORK(5)  the maximum order to be allowed.  The default
  !                   value is 12 if METH = 1, and 5 if METH = 2.
  !                   If MAXORD exceeds the default value, it will
  !                   be reduced to the default value.
  !                   If MAXORD is changed during the problem, it may
  !                   cause the current order to be reduced.
  !
  ! MXSTEP  IWORK(6)  maximum number of (internally defined) steps
  !                   allowed during one call to the solver.
  !                   The default value is 500.
  !
  ! MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
  !                   warning that T + H = T on a step (H = step size).
  !                   This must be positive to result in a non-default
  !                   value.  The default value is 10.
  !-----------------------------------------------------------------------
  ! Optional Outputs.
  !
  ! As optional additional output from DLSOIBT, the variables listed
  ! below are quantities related to the performance of DLSOIBT
  ! which are available to the user.  These are communicated by way of
  ! the work arrays, but also have internal mnemonic names as shown.
  ! Except where stated otherwise, all of these outputs are defined
  ! on any successful return from DLSOIBT, and on any return with
  ! ISTATE = -1, -2, -4, -5, -6, or -7.  On a return with -3 (illegal
  ! input) or -8, they will be unchanged from their existing values
  ! (if any), except possibly for TOLSF, LENRW, and LENIW.
  ! On any error return, outputs relevant to the error will be defined,
  ! as noted below.
  !
  ! Name    Location      Meaning
  !
  ! HU      RWORK(11) the step size in t last used (successfully).
  !
  ! HCUR    RWORK(12) the step size to be attempted on the next step.
  !
  ! TCUR    RWORK(13) the current value of the independent variable
  !                   which the solver has actually reached, i.e. the
  !                   current internal mesh point in t.  On output, TCUR
  !                   will always be at least as far as the argument
  !                   T, but may be farther (if interpolation was done).
  !
  ! TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
  !                   computed when a request for too much accuracy was
  !                   detected (ISTATE = -3 if detected at the start of
  !                   the problem, ISTATE = -2 otherwise).  If ITOL is
  !                   left unaltered but RTOL and ATOL are uniformly
  !                   scaled up by a factor of TOLSF for the next call,
  !                   then the solver is deemed likely to succeed.
  !                   (The user may also ignore TOLSF and alter the
  !                   tolerance parameters in any other way appropriate.)
  !
  ! NST     IWORK(11) the number of steps taken for the problem so far.
  !
  ! NRE     IWORK(12) the number of residual evaluations (RES calls)
  !                   for the problem so far.
  !
  ! NJE     IWORK(13) the number of Jacobian evaluations (each involving
  !                   an evaluation of a and dr/dy) for the problem so
  !                   far.  This equals the number of calls to ADDA and
  !                   (if MITER = 1) to JAC, and the number of matrix
  !                   LU decompositions.
  !
  ! NQU     IWORK(14) the method order last used (successfully).
  !
  ! NQCUR   IWORK(15) the order to be attempted on the next step.
  !
  ! IMXER   IWORK(16) the index of the component of largest magnitude in
  !                   the weighted local error vector ( E(i)/EWT(i) ),
  !                   on an error return with ISTATE = -4 or -5.
  !
  ! LENRW   IWORK(17) the length of RWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  ! LENIW   IWORK(18) the length of IWORK actually required.
  !                   This is defined on normal returns and on an illegal
  !                   input return for insufficient storage.
  !
  !
  ! The following two arrays are segments of the RWORK array which
  ! may also be of interest to the user as optional outputs.
  ! For each array, the table below gives its internal name,
  ! its base address in RWORK, and its description.
  !
  ! Name    Base Address      Description
  !
  ! YH      21             the Nordsieck history array, of size NYH by
  !                        (NQCUR + 1), where NYH is the initial value
  !                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
  !                        of YH contains HCUR**j/factorial(j) times
  !                        the j-th derivative of the interpolating
  !                        polynomial currently representing the solution,
  !                        evaluated at t = TCUR.
  !
  ! ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
  !                        corrections on each step, scaled on output to
  !                        represent the estimated local error in y on
  !                        the last step.  This is the vector E in the
  !                        description of the error control.  It is
  !                        defined only on a return from DLSOIBT with
  !                        ISTATE = 2.
  !
  !-----------------------------------------------------------------------
  ! Part 2.  Other Routines Callable.
  !
  ! The following are optional calls which the user may make to
  ! gain additional capabilities in conjunction with DLSOIBT.
  ! (The routines XSETUN and XSETF are designed to conform to the
  ! SLATEC error handling package.)
  !
  !     Form of Call                  Function
  !   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
  !                             output of messages from DLSOIBT, if
  !                             the default is not desired.
  !                             The default value of LUN is 6.
  !
  !   CALL XSETF(MFLAG)         Set a flag to control the printing of
  !                             messages by DLSOIBT.
  !                             MFLAG = 0 means do not print. (Danger:
  !                             This risks losing valuable information.)
  !                             MFLAG = 1 means print (the default).
  !
  !                             Either of the above calls may be made at
  !                             any time and will take effect immediately.
  !
  !   CALL DSRCOM(RSAV,ISAV,JOB) saves and restores the contents of
  !                             the internal Common blocks used by
  !                             DLSOIBT (see Part 3 below).
  !                             RSAV must be a real array of length 218
  !                             or more, and ISAV must be an integer
  !                             array of length 37 or more.
  !                             JOB=1 means save Common into RSAV/ISAV.
  !                             JOB=2 means restore Common from RSAV/ISAV.
  !                                DSRCOM is useful if one is
  !                             interrupting a run and restarting
  !                             later, or alternating between two or
  !                             more problems solved with DLSOIBT.
  !
  !   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
  !        (see below)          orders, at a specified point t, if
  !                             desired.  It may be called only after
  !                             a successful return from DLSOIBT.
  !
  ! The detailed instructions for using DINTDY are as follows.
  ! The form of the call is:
  !
  !   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
  !
  ! The input parameters are:
  !
  ! T         = value of independent variable where answers are desired
  !             (normally the same as the t last returned by DLSOIBT).
  !             For valid results, T must lie between TCUR - HU and TCUR.
  !             (See optional outputs for TCUR and HU.)
  ! K         = integer order of the derivative desired.  K must satisfy
  !             0 .le. K .le. NQCUR, where NQCUR is the current order
  !             (see optional outputs).  The capability corresponding
  !             to K = 0, i.e. computing y(t), is already provided
  !             by DLSOIBT directly.  Since NQCUR .ge. 1, the first
  !             derivative dy/dt is always available with DINTDY.
  ! RWORK(21) = the base address of the history array YH.
  ! NYH       = column length of YH, equal to the initial value of NEQ.
  !
  ! The output parameters are:
  !
  ! DKY       = a real array of length NEQ containing the computed value
  !             of the K-th derivative of y(t).
  ! IFLAG     = integer flag, returned as 0 if K and T were legal,
  !             -1 if K was illegal, and -2 if T was illegal.
  !             On an error return, a message is also written.
  !-----------------------------------------------------------------------
  ! Part 3.  Common Blocks.
  !
  ! If DLSOIBT is to be used in an overlay situation, the user
  ! must declare, in the primary overlay, the variables in:
  !   (1) the call sequence to DLSOIBT, and
  !   (2) the internal Common block
  !         /DLS001/  of length  255  (218 double precision words
  !                      followed by 37 integer words),
  !
  ! If DLSOIBT is used on a system in which the contents of internal
  ! Common blocks are not preserved between calls, the user should
  ! declare the above Common block in the calling program to insure
  ! that their contents are preserved.
  !
  ! If the solution of a given problem by DLSOIBT is to be interrupted
  ! and then later continued, such as when restarting an interrupted run
  ! or alternating between two or more problems, the user should save,
  ! following the return from the last DLSOIBT call prior to the
  ! interruption, the contents of the call sequence variables and the
  ! internal Common blocks, and later restore these values before the
  ! next DLSOIBT call for that problem.  To save and restore the Common
  ! blocks, use Subroutine DSRCOM (see Part 2 above).
  !
  !-----------------------------------------------------------------------
  ! Part 4.  Optionally Replaceable Solver Routines.
  !
  ! Below are descriptions of two routines in the DLSOIBT package which
  ! relate to the measurement of errors.  Either routine can be
  ! replaced by a user-supplied version, if desired.  However, since such
  ! a replacement may have a major impact on performance, it should be
  ! done only when absolutely necessary, and only with great caution.
  ! (Note: The means by which the package version of a routine is
  ! superseded by the user's version may be system-dependent.)
  !
  ! (a) DEWSET.
  ! The following subroutine is called just before each internal
  ! integration step, and sets the array of error weights, EWT, as
  ! described under ITOL/RTOL/ATOL above:
  !     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
  ! where NEQ, ITOL, RTOL, and ATOL are as in the DLSOIBT call sequence,
  ! YCUR contains the current dependent variable vector, and
  ! EWT is the array of weights set by DEWSET.
  !
  ! If the user supplies this subroutine, it must return in EWT(i)
  ! (i = 1,...,NEQ) a positive quantity suitable for comparing errors
  ! in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
  ! routine (see below), and also used by DLSOIBT in the computation
  ! of the optional output IMXER, the diagonal Jacobian approximation,
  ! and the increments for difference quotient Jacobians.
  !
  ! In the user-supplied version of DEWSET, it may be desirable to use
  ! the current values of derivatives of y.  Derivatives up to order NQ
  ! are available from the history array YH, described above under
  ! optional outputs.  In DEWSET, YH is identical to the YCUR array,
  ! extended to NQ + 1 columns with a column length of NYH and scale
  ! factors of H**j/factorial(j).  On the first call for the problem,
  ! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
  ! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
  ! can be obtained by including in DEWSET the statements:
  !     DOUBLE PRECISION RLS
  !     COMMON /DLS001/ RLS(218),ILS(37)
  !     NQ = ILS(33)
  !     NST = ILS(34)
  !     H = RLS(212)
  ! Thus, for example, the current value of dy/dt can be obtained as
  ! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
  ! unnecessary when NST = 0).
  !
  ! (b) DVNORM.
  ! The following is a real function routine which computes the weighted
  ! root-mean-square norm of a vector v:
  !     D = DVNORM (N, V, W)
  ! where:
  !   N = the length of the vector,
  !   V = real array of length N containing the vector,
  !   W = real array of length N containing weights,
  !   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
  ! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
  ! EWT is as set by Subroutine DEWSET.
  !
  ! If the user supplies this function, it should return a non-negative
  ! value of DVNORM suitable for use in the error control in DLSOIBT.
  ! None of the arguments should be altered by DVNORM.
  ! For example, a user-supplied DVNORM routine might:
  !   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
  !   -ignore some components of V in the norm, with the effect of
  !    suppressing the error control on those components of y.
  !-----------------------------------------------------------------------
  !
  !***REVISION HISTORY  (YYYYMMDD)
  ! 19840625  DATE WRITTEN
  ! 19870330  Major update: corrected comments throughout;
  !           removed TRET from Common; rewrote EWSET with 4 loops;
  !           fixed t test in INTDY; added Cray directives in STODI;
  !           in STODI, fixed DELP init. and logic around PJAC call;
  !           combined routines to save/restore Common;
  !           passed LEVEL = 0 in error message calls (except run abort).
  ! 20010425  Major update: convert source lines to upper case;
  !           added *DECK lines; changed from 1 to * in dummy dimensions;
  !           changed names R1MACH/D1MACH to RUMACH/DUMACH;
  !           renamed routines for uniqueness across single/double prec.;
  !           converted intrinsic names to generic form;
  !           removed ILLIN and NTREP (data loaded) from Common;
  !           removed all 'own' variables from Common;
  !           changed error messages to quoted strings;
  !           replaced XERRWV/XERRWD with 1993 revised version;
  !           converted prologues, comments, error messages to mixed case;
  !           converted arithmetic IF statements to logical IF statements;
  !           numerous corrections to prologues and internal comments.
  ! 20010507  Converted single precision source to double precision.
  ! 20020502  Corrected declarations in descriptions of user routines.
  ! 20031105  Restored 'own' variables to Common block, to enable
  !           interrupt/restart feature.
  ! 20031112  Added SAVE statements for data-loaded constants.
  ! 20031117  Changed internal names NRE, LSAVR to NFE, LSAVF resp.
  !
  !-----------------------------------------------------------------------
  ! Other routines in the DLSOIBT package.
  !
  ! In addition to Subroutine DLSOIBT, the DLSOIBT package includes the
  ! following subroutines and function routines:
  !  DAIGBT   computes the initial value of the vector
  !             dy/dt = A-inverse * g
  !  DINTDY   computes an interpolated value of the y vector at t = TOUT.
  !  DSTODI   is the core integrator, which does one step of the
  !           integration and the associated error control.
  !  DCFODE   sets all method coefficients and test constants.
  !  DEWSET   sets the error weight vector EWT before each step.
  !  DVNORM   computes the weighted RMS-norm of a vector.
  !  DSRCOM   is a user-callable routine to save and restore
  !           the contents of the internal Common blocks.
  !  DPJIBT   computes and preprocesses the Jacobian matrix
  !           and the Newton iteration matrix P.
  !  DSLSBT   manages solution of linear system in chord iteration.
  !  DDECBT and DSOLBT   are routines for solving block-tridiagonal
  !           systems of linear algebraic equations.
  !  DGEFA and DGESL   are routines from LINPACK for solving full
  !           systems of linear algebraic equations.
  !  DDOT     is one of the basic linear algebra modules (BLAS).
  !  DUMACH   computes the unit roundoff in a machine-independent manner.
  !  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
  !           error messages and warnings.  XERRWD is machine-dependent.
  ! Note:  DVNORM, DDOT, DUMACH, IXSAV, and IUMACH are function routines.
  ! All the others are subroutines.
  !
  !-----------------------------------------------------------------------
        !EXTERNAL DPJIBT , DSLSBT
        INTEGER INIt , MXStep , MXHnil , NHNil , NSLast , NYH , IOWns ,   &
       &        ICF , IERpj , IERsl , JCUr , JSTart , KFLag , L , LYH ,   &
       &        LEWt , LACor , LSAvf , LWM , LIWm , METh , MITer ,        &
       &        MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST , NFE ,     &
       &        NJE , NQU
        INTEGER i , i1 , i2 , ier , iflag , imxer , ires , kgo , leniw ,  &
       &        lenrw , lenwm , lp , lyd0 , mb , mord , mxhnl0 , mxstp0 , &
       &        nb
        DOUBLE PRECISION ROWns , CCMax , EL0 , H , HMIn , HMXi , HU , RC ,&
       &                 TN , UROund
        DOUBLE PRECISION atoli , ayi , big , ewti , h0 , hmax , hmx , rh ,&
       &                 rtoli , tcrit , tdist , tnext , tol , tolsf ,    &
       &                 tp , size , sum , w0
        DIMENSION mord(2)
        LOGICAL ihit
        CHARACTER*60 msg
        SAVE mord , mxstp0 , mxhnl0
  !-----------------------------------------------------------------------
  ! The following internal Common block contains
  ! (a) variables which are local to any subroutine but whose values must
  !     be preserved between calls to the routine ("own" variables), and
  ! (b) variables which are communicated between subroutines.
  ! The block DLS001 is declared in subroutines DLSOIBT, DINTDY, DSTODI,
  ! DPJIBT, and DSLSBT.
  ! Groups of variables are replaced by dummy arrays in the Common
  ! declarations in routines where those variables are not used.
  !-----------------------------------------------------------------------
        COMMON /DLS001/ ROWns(209) , CCMax , EL0 , H , HMIn , HMXi , HU , &
       &                RC , TN , UROund , INIt , MXStep , MXHnil ,       &
       &                NHNil , NSLast , NYH , IOWns(6) , ICF , IERpj ,   &
       &                IERsl , JCUr , JSTart , KFLag , L , LYH , LEWt ,  &
       &                LACor , LSAvf , LWM , LIWm , METh , MITer ,       &
       &                MAXord , MAXcor , MSBp , MXNcf , N , NQ , NST ,   &
       &                NFE , NJE , NQU
  !
        DATA mord(1) , mord(2)/12 , 5/ , mxstp0/500/ , mxhnl0/10/
  !-----------------------------------------------------------------------
  ! Block A.
  ! This code block is executed on every call.
  ! It tests ISTATE and ITASK for legality and branches appropriately.
  ! If ISTATE .gt. 1 but the flag INIT shows that initialization has
  ! not yet been done, an error return occurs.
  ! If ISTATE = 0 or 1 and TOUT = T, return immediately.
  !-----------------------------------------------------------------------
        IF ( Istate.LT.0 .OR. Istate.GT.3 ) THEN
  !-----------------------------------------------------------------------
  ! Block I.
  ! The following block handles all error returns due to illegal input
  ! (ISTATE = -3), as detected before calling the core integrator.
  ! First the error message routine is called.  If the illegal input
  ! is a negative ISTATE, the run is aborted (apparent infinite loop).
  !-----------------------------------------------------------------------
           msg = 'DLSOIBT- ISTATE (=I1) illegal.'
           CALL XERRWD(msg,30,1,0,1,Istate,0,0,0.0D0,0.0D0)
           IF ( Istate.GE.0 ) GOTO 3100
  !
           msg = 'DLSOIBT- Run aborted.. apparent infinite loop.    '
           CALL XERRWD(msg,50,303,2,0,0,0,0,0.0D0,0.0D0)
           GOTO 99999
        ELSE
           IF ( Itask.LT.1 .OR. Itask.GT.5 ) THEN
              msg = 'DLSOIBT- ITASK (=I1) illegal. '
              CALL XERRWD(msg,30,2,0,1,Itask,0,0,0.0D0,0.0D0)
              GOTO 3100
           ELSE
              IF ( Istate.LE.1 ) THEN
                 INIt = 0
                 IF ( Tout.EQ.T ) RETURN
              ELSEIF ( INIt.EQ.0 ) THEN
                 msg =                                                    &
       &              'DLSOIBT- ISTATE.gt.1 but DLSOIBT not initialized. '
                 CALL XERRWD(msg,50,3,0,0,0,0,0,0.0D0,0.0D0)
                 GOTO 3100
              ELSEIF ( Istate.EQ.2 ) THEN
                 GOTO 50
              ENDIF
  !-----------------------------------------------------------------------
  ! Block B.
  ! The next code block is executed for the initial call (ISTATE = 0 or 1)
  ! or for a continuation call with parameter changes (ISTATE = 3).
  ! It contains checking of all inputs and various initializations.
  !
  ! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
  ! MF, MB, and NB.
  !-----------------------------------------------------------------------
              IF ( Neq(1).LE.0 ) THEN
                 msg = 'DLSOIBT- NEQ (=I1) .lt. 1     '
                 CALL XERRWD(msg,30,4,0,1,Neq(1),0,0,0.0D0,0.0D0)
                 GOTO 3100
              ELSE
                 IF ( Istate.GT.1 ) THEN
                    IF ( Neq(1).GT.N ) THEN
                       msg =                                              &
       &              'DLSOIBT- ISTATE = 3 and NEQ increased (I1 to I2). '
                       CALL XERRWD(msg,50,5,0,2,N,Neq(1),0,0.0D0,0.0D0)
                       GOTO 3100
                    ENDIF
                 ENDIF
                 N = Neq(1)
                 IF ( Itol.LT.1 .OR. Itol.GT.4 ) THEN
                    msg = 'DLSOIBT- ITOL (=I1) illegal.  '
                    CALL XERRWD(msg,30,6,0,1,Itol,0,0,0.0D0,0.0D0)
                    GOTO 3100
                 ELSEIF ( Iopt.LT.0 .OR. Iopt.GT.1 ) THEN
                    msg = 'DLSOIBT- IOPT (=I1) illegal.  '
                    CALL XERRWD(msg,30,7,0,1,Iopt,0,0,0.0D0,0.0D0)
                    GOTO 3100
                 ELSE
                    METh = Mf/10
                    MITer = Mf - 10*METh
                    IF ( METh.LT.1 .OR. METh.GT.2 ) GOTO 2300
                    IF ( MITer.LT.1 .OR. MITer.GT.2 ) GOTO 2300
                    mb = Iwork(1)
                    nb = Iwork(2)
                    IF ( mb.LT.1 .OR. mb.GT.N ) GOTO 2400
                    IF ( nb.LT.4 ) THEN
                       msg = 'DLSOIBT- NB (=I1) .lt. 4 illegal.       '
                       CALL XERRWD(msg,40,10,0,1,nb,0,0,0.0D0,0.0D0)
                       GOTO 3100
                    ELSE
                       IF ( mb*nb.NE.N ) GOTO 2400
  ! Next process and check the optional inputs. --------------------------
                       IF ( Iopt.EQ.1 ) THEN
                          MAXord = Iwork(5)
                          IF ( MAXord.LT.0 ) THEN
                             msg = 'DLSOIBT- MAXORD (=I1) .lt. 0  '
                             CALL XERRWD(msg,30,11,0,1,MAXord,0,0,0.0D0,  &
       &                                 0.0D0)
                             GOTO 3100
                          ELSE
                             IF ( MAXord.EQ.0 ) MAXord = 100
                             MAXord = MIN(MAXord,mord(METh))
                             MXStep = Iwork(6)
                             IF ( MXStep.LT.0 ) THEN
                                msg = 'DLSOIBT- MXSTEP (=I1) .lt. 0  '
                                CALL XERRWD(msg,30,12,0,1,MXStep,0,0,     &
       &                           0.0D0,0.0D0)
                                GOTO 3100
                             ELSE
                                IF ( MXStep.EQ.0 ) MXStep = mxstp0
                                MXHnil = Iwork(7)
                                IF ( MXHnil.LT.0 ) THEN
                                   msg = 'DLSOIBT- MXHNIL (=I1) .lt. 0  '
                                   CALL XERRWD(msg,30,13,0,1,MXHnil,0,0,  &
       &                              0.0D0,0.0D0)
                                   GOTO 3100
                                ELSE
                                   IF ( MXHnil.EQ.0 ) MXHnil = mxhnl0
                                   IF ( Istate.LE.1 ) THEN
                                      h0 = Rwork(5)
                                      IF ( (Tout-T)*h0.LT.0.0D0 ) THEN
                                         msg =                            &
       &                        'DLSOIBT- TOUT (=R1) behind T (=R2)      '
                                         CALL XERRWD(msg,40,14,0,0,0,0,2, &
       &                                    Tout,T)
                                         msg =                            &
       &              '      Integration direction is given by H0 (=R1)  '
                                         CALL XERRWD(msg,50,14,0,0,0,0,1, &
       &                                    h0,0.0D0)
                                         GOTO 3100
                                      ENDIF
                                   ENDIF
                                   hmax = Rwork(6)
                                   IF ( hmax.LT.0.0D0 ) THEN
                                      msg =                               &
       &                                 'DLSOIBT- HMAX (=R1) .lt. 0.0  '
                                      CALL XERRWD(msg,30,15,0,0,0,0,1,    &
       &                                 hmax,0.0D0)
                                      GOTO 3100
                                   ELSE
                                      HMXi = 0.0D0
                                      IF ( hmax.GT.0.0D0 )                &
       &                                 HMXi = 1.0D0/hmax
                                      HMIn = Rwork(7)
                                      IF ( HMIn.LT.0.0D0 ) THEN
                                         msg =                            &
       &                                  'DLSOIBT- HMIN (=R1) .lt. 0.0  '
                                         CALL XERRWD(msg,30,16,0,0,0,0,1, &
       &                                    HMIn,0.0D0)
                                         GOTO 3100
                                      ENDIF
                                   ENDIF
                                ENDIF
                             ENDIF
                          ENDIF
                       ELSE
                          MAXord = mord(METh)
                          MXStep = mxstp0
                          MXHnil = mxhnl0
                          IF ( Istate.LE.1 ) h0 = 0.0D0
                          HMXi = 0.0D0
                          HMIn = 0.0D0
                       ENDIF
  !-----------------------------------------------------------------------
  ! Set work array pointers and check lengths LRW and LIW.
  ! Pointers to segments of RWORK and IWORK are named by prefixing L to
  ! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
  ! Segments of RWORK (in order) are denoted YH, WM, EWT, SAVR, ACOR.
  !-----------------------------------------------------------------------
                       LYH = 21
                       IF ( Istate.LE.1 ) NYH = N
                       LWM = LYH + (MAXord+1)*NYH
                       lenwm = 3*mb*mb*nb + 2
                       LEWt = LWM + lenwm
                       LSAvf = LEWt + N
                       LACor = LSAvf + N
                       lenrw = LACor + N - 1
                       Iwork(17) = lenrw
                       LIWm = 1
                       leniw = 20 + N
                       Iwork(18) = leniw
                       IF ( lenrw.GT.Lrw ) THEN
                          msg =                                           &
       &    'DLSOIBT- RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                          CALL XERRWD(msg,60,17,0,2,lenrw,Lrw,0,0.0D0,    &
       &                              0.0D0)
                          GOTO 3100
                       ELSEIF ( leniw.GT.Liw ) THEN
                          msg =                                           &
       &    'DLSOIBT- IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                          CALL XERRWD(msg,60,18,0,2,leniw,Liw,0,0.0D0,    &
       &                              0.0D0)
                          GOTO 3100
                       ELSE
  ! Check RTOL and ATOL for legality. ------------------------------------
                          rtoli = Rtol(1)
                          atoli = Atol(1)
                          DO i = 1 , N
                             IF ( Itol.GE.3 ) rtoli = Rtol(i)
                             IF ( Itol.EQ.2 .OR. Itol.EQ.4 )              &
       &                          atoli = Atol(i)
                             IF ( rtoli.LT.0.0D0 ) GOTO 2500
                             IF ( atoli.LT.0.0D0 ) GOTO 2600
                          ENDDO
                          IF ( Istate.LE.1 ) THEN
  !-----------------------------------------------------------------------
  ! Block C.
  ! The next block is for the initial call only (ISTATE = 0 or 1).
  ! It contains all remaining initializations, the call to DAIGBT
  ! (if ISTATE = 1), and the calculation of the initial step size.
  ! The error weights in EWT are inverted after being loaded.
  !-----------------------------------------------------------------------
                             UROund = DUMACH()
                             TN = T
                             IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
                                tcrit = Rwork(1)
                                IF ( (tcrit-Tout)*(Tout-T).LT.0.0D0 )     &
       &                             GOTO 2900
                                IF ( h0.NE.0.0D0 .AND. (T+h0-tcrit)       &
       &                             *h0.GT.0.0D0 ) h0 = tcrit - T
                             ENDIF
                             JSTart = 0
                             Rwork(LWM) = SQRT(UROund)
                             NHNil = 0
                             NST = 0
                             NFE = 0
                             NJE = 0
                             NSLast = 0
                             HU = 0.0D0
                             NQU = 0
                             CCMax = 0.3D0
                             MAXcor = 3
                             MSBp = 20
                             MXNcf = 10
  ! Compute initial dy/dt, if necessary, and load it and initial Y into YH
                             lyd0 = LYH + NYH
                             lp = LWM + 1
                             IF ( Istate.EQ.1 ) THEN
  ! Initial dy/dt was supplied.  Load into YH (LYD0 points to YH(*,2).). -
                                DO i = 1 , N
                                   Rwork(i+LYH-1) = Y(i)
                                   Rwork(i+lyd0-1) = Ydoti(i)
                                ENDDO
                             ELSE
  ! DLSOIBT must compute initial dy/dt (LYD0 points to YH(*,2)). ---------
                                CALL DAIGBT(RES,ADDA,Neq,T,Y,Rwork(lyd0), &
       &                           mb,nb,Rwork(lp),Iwork(21),ier)
                                NFE = NFE + 1
                                IF ( ier.LT.0 ) THEN
  ! DAIGBT failed because a diagonal block of A matrix was singular. -----
                                   ier = -ier
                                   msg =                                  &
       &    'DLSOIBT- Attempt to initialize dy/dt failed:  Matrix A has a'
                                   CALL XERRWD(msg,60,207,0,0,0,0,0,0.0D0,&
       &                              0.0D0)
                                   msg =                                  &
       &              '      singular diagonal block, block no. = (I1)   '
                                   CALL XERRWD(msg,50,207,0,1,ier,0,0,    &
       &                              0.0D0,0.0D0)
                                   Istate = -8
                                   RETURN
                                ELSEIF ( ier.GT.0 ) THEN
  ! DAIGBT failed because RES set IRES to 2 or 3. ------------------------
                                   msg =                                  &
       &              'DLSOIBT- Attempt to initialize dy/dt failed       '
                                   CALL XERRWD(msg,50,208,0,0,0,0,0,0.0D0,&
       &                              0.0D0)
                                   msg =                                  &
       &              '      because residual routine set its error flag '
                                   CALL XERRWD(msg,50,208,0,0,0,0,0,0.0D0,&
       &                              0.0D0)
                                   msg = '      to IRES = (I1)'
                                   CALL XERRWD(msg,20,208,0,1,ier,0,0,    &
       &                              0.0D0,0.0D0)
                                   Istate = -8
                                   RETURN
                                ELSE
                                   DO i = 1 , N
                                      Rwork(i+LYH-1) = Y(i)
                                   ENDDO
                                ENDIF
                             ENDIF
  ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                             NQ = 1
                             H = 1.0D0
                             CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),     &
       &                                 Rwork(LEWt))
                             DO i = 1 , N
                                IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 2700
                                Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
                             ENDDO
  !-----------------------------------------------------------------------
  ! The coding below computes the step size, H0, to be attempted on the
  ! first step, unless the user has supplied a value for this.
  ! First check that TOUT - T differs significantly from zero.
  ! A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
  ! if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
  ! so as to be between 100*UROUND and 1.0E-3.
  ! Then the computed value H0 is given by..
  !                                      NEQ
  !   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( YDOT(i)/ywt(i) )**2  )
  !                                       1
  ! where   w0      = MAX ( ABS(T), ABS(TOUT) ),
  !         YDOT(i) = i-th component of initial value of dy/dt,
  !         ywt(i)  = EWT(i)/TOL  (a weight for y(i)).
  ! The sign of H0 is inferred from the initial values of TOUT and T.
  !-----------------------------------------------------------------------
                             IF ( h0.EQ.0.0D0 ) THEN
                                tdist = ABS(Tout-T)
                                w0 = MAX(ABS(T),ABS(Tout))
                                IF ( tdist.LT.2.0D0*UROund*w0 ) THEN
                                   msg =                                  &
       &    'DLSOIBT- TOUT(=R1) too close to T(=R2) to start integration.'
                                   CALL XERRWD(msg,60,22,0,0,0,0,2,Tout,T)
                                   GOTO 3100
                                ELSE
                                   tol = Rtol(1)
                                   IF ( Itol.GT.2 ) THEN
                                      DO i = 1 , N
                                         tol = MAX(tol,Rtol(i))
                                      ENDDO
                                   ENDIF
                                   IF ( tol.LE.0.0D0 ) THEN
                                      atoli = Atol(1)
                                      DO i = 1 , N
                                         IF ( Itol.EQ.2 .OR. Itol.EQ.4 )  &
       &                                    atoli = Atol(i)
                                         ayi = ABS(Y(i))
                                         IF ( ayi.NE.0.0D0 )              &
       &                                    tol = MAX(tol,atoli/ayi)
                                      ENDDO
                                   ENDIF
                                   tol = MAX(tol,100.0D0*UROund)
                                   tol = MIN(tol,0.001D0)
                                   sum = DVNORM(N,Rwork(lyd0),Rwork(LEWt))
                                   sum = 1.0D0/(tol*w0*w0) + tol*sum**2
                                   h0 = 1.0D0/SQRT(sum)
                                   h0 = MIN(h0,tdist)
                                   h0 = SIGN(h0,Tout-T)
                                ENDIF
                             ENDIF
  ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                             rh = ABS(h0)*HMXi
                             IF ( rh.GT.1.0D0 ) h0 = h0/rh
  ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                             H = h0
                             DO i = 1 , N
                                Rwork(i+lyd0-1) = h0*Rwork(i+lyd0-1)
                             ENDDO
                             GOTO 700
                          ELSE
  ! If ISTATE = 3, set flag to signal parameter changes to DSTODI. -------
                             JSTart = -1
                             IF ( NQ.GT.MAXord ) THEN
  ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into YDOTI.---------
                                DO i = 1 , N
                                   Ydoti(i) = Rwork(i+LWM-1)
                                ENDDO
                             ENDIF
  ! Reload WM(1) = RWORK(lWM), since lWM may have changed. ---------------
                             Rwork(LWM) = SQRT(UROund)
                             IF ( N.NE.NYH ) THEN
  ! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
                                i1 = LYH + L*NYH
                                i2 = LYH + (MAXord+1)*NYH - 1
                                IF ( i1.LE.i2 ) THEN
                                   DO i = i1 , i2
                                      Rwork(i) = 0.0D0
                                   ENDDO
                                ENDIF
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  ! Block D.
  ! The next code block is for continuation calls only (ISTATE = 2 or 3)
  ! and is to check stop conditions before taking a step.
  !-----------------------------------------------------------------------
   50      NSLast = NST
           GOTO (100,600,200,300,400) , Itask
        ENDIF
   100  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 3000
        T = Tout
        GOTO 1400
   200  tp = TN - HU*(1.0D0+100.0D0*UROund)
        IF ( (tp-Tout)*H.GT.0.0D0 ) THEN
           msg =                                                          &
       &    'DLSOIBT- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
           CALL XERRWD(msg,60,23,0,1,Itask,0,2,Tout,tp)
           GOTO 3100
        ELSE
           IF ( (TN-Tout)*H.GE.0.0D0 ) GOTO 1300
           GOTO 600
        ENDIF
   300  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2800
        IF ( (tcrit-Tout)*H.LT.0.0D0 ) GOTO 2900
        IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 500
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        IF ( iflag.NE.0 ) GOTO 3000
        T = Tout
        GOTO 1400
   400  tcrit = Rwork(1)
        IF ( (TN-tcrit)*H.GT.0.0D0 ) GOTO 2800
   500  hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
        IF ( ihit ) GOTO 1300
        tnext = TN + H*(1.0D0+4.0D0*UROund)
        IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
           H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
           IF ( Istate.EQ.2 ) JSTart = -2
        ENDIF
  !-----------------------------------------------------------------------
  ! Block E.
  ! The next block is normally executed for all calls and contains
  ! the call to the one-step core integrator DSTODI.
  !
  ! This is a looping point for the integration steps.
  !
  ! First check for too many steps being taken, update EWT (if not at
  ! start of problem), check for too much accuracy being requested, and
  ! check for H below the roundoff level in T.
  !-----------------------------------------------------------------------
   600  IF ( (NST-NSLast).GE.MXStep ) THEN
  !-----------------------------------------------------------------------
  ! Block H.
  ! The following block handles all unsuccessful returns other than
  ! those for illegal input.  First the error message routine is called.
  ! If there was an error test or convergence test failure, IMXER is set.
  ! Then Y is loaded from YH and T is set to TN.
  ! The optional outputs are loaded into the work arrays before returning.
  !-----------------------------------------------------------------------
  ! The maximum number of steps was taken before reaching TOUT. ----------
           msg = 'DLSOIBT- At current T (=R1), MXSTEP (=I1) steps   '
           CALL XERRWD(msg,50,201,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      taken on this call before reaching TOUT     '
           CALL XERRWD(msg,50,201,0,1,MXStep,0,1,TN,0.0D0)
           Istate = -1
           GOTO 2000
        ELSE
           CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
           DO i = 1 , N
              IF ( Rwork(i+LEWt-1).LE.0.0D0 ) GOTO 1500
              Rwork(i+LEWt-1) = 1.0D0/Rwork(i+LEWt-1)
           ENDDO
        ENDIF
   700  tolsf = UROund*DVNORM(N,Rwork(LYH),Rwork(LEWt))
        IF ( tolsf.LE.1.0D0 ) THEN
           IF ( (TN+H).EQ.TN ) THEN
              NHNil = NHNil + 1
              IF ( NHNil.LE.MXHnil ) THEN
                 msg =                                                    &
       &              'DLSOIBT- Warning..Internal T (=R1) and H (=R2) are'
                 CALL XERRWD(msg,50,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &    '      such that in the machine, T + H = T on the next step  '
                 CALL XERRWD(msg,60,101,0,0,0,0,0,0.0D0,0.0D0)
                 msg =                                                    &
       &              '     (H = step size). Solver will continue anyway.'
                 CALL XERRWD(msg,50,101,0,0,0,0,2,TN,H)
                 IF ( NHNil.GE.MXHnil ) THEN
                    msg =                                                 &
       &              'DLSOIBT- Above warning has been issued I1 times.  '
                    CALL XERRWD(msg,50,102,0,0,0,0,0,0.0D0,0.0D0)
                    msg =                                                 &
       &              '     It will not be issued again for this problem.'
                    CALL XERRWD(msg,50,102,0,1,MXHnil,0,0,0.0D0,0.0D0)
                 ENDIF
              ENDIF
           ENDIF
  !-----------------------------------------------------------------------
  !     CALL DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM,RES,
  !                 ADDA,JAC,DPJIBT,DSLSBT)
  ! Note: SAVF in DSTODI occupies the same space as YDOTI in DLSOIBT.
  !-----------------------------------------------------------------------
           CALL DSTODI(Neq,Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),Ydoti, &
       &               Rwork(LSAvf),Rwork(LACor),Rwork(LWM),Iwork(LIWm),  &
       &               RES,ADDA,JAC,DPJIBT,DSLSBT)
           kgo = 1 - KFLag
           GOTO (800,1600,1700,1300,1800) , kgo
        ELSE
           tolsf = tolsf*2.0D0
           IF ( NST.EQ.0 ) THEN
              msg = 'DLSOIBT- At start of problem, too much accuracy   '
              CALL XERRWD(msg,50,26,0,0,0,0,0,0.0D0,0.0D0)
              msg =                                                       &
       &    '      requested for precision of machine..  See TOLSF (=R1) '
              CALL XERRWD(msg,60,26,0,0,0,0,1,tolsf,0.0D0)
              Rwork(14) = tolsf
              GOTO 3100
           ELSE
  ! Too much accuracy requested for machine precision. -------------------
              msg = 'DLSOIBT- At T (=R1), too much accuracy requested  '
              CALL XERRWD(msg,50,203,0,0,0,0,0,0.0D0,0.0D0)
              msg = '      for precision of machine..  See TOLSF (=R2) '
              CALL XERRWD(msg,50,203,0,0,0,0,2,TN,tolsf)
              Rwork(14) = tolsf
              Istate = -2
              GOTO 2100
           ENDIF
        ENDIF
  !
  ! KGO = 1:success; 2:error test failure; 3:convergence failure;
  !       4:RES ordered return; 5:RES returned error.
  !-----------------------------------------------------------------------
  ! Block F.
  ! The following block handles the case of a successful return from the
  ! core integrator (KFLAG = 0).  Test for stop conditions.
  !-----------------------------------------------------------------------
   800  INIt = 1
        GOTO (900,1300,1000,1100,1200) , Itask
  ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
   900  IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
        T = Tout
        GOTO 1400
  ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
   1000 IF ( (TN-Tout)*H.LT.0.0D0 ) GOTO 600
        GOTO 1300
  ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
   1100 IF ( (TN-Tout)*H.LT.0.0D0 ) THEN
           hmx = ABS(TN) + ABS(H)
           ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
           IF ( ihit ) GOTO 1300
           tnext = TN + H*(1.0D0+4.0D0*UROund)
           IF ( (tnext-tcrit)*H.GT.0.0D0 ) THEN
              H = (tcrit-TN)*(1.0D0-4.0D0*UROund)
              JSTart = -2
           ENDIF
           GOTO 600
        ELSE
           CALL DINTDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
           T = Tout
           GOTO 1400
        ENDIF
  ! ITASK = 5.  see if TCRIT was reached and jump to exit. ---------------
   1200 hmx = ABS(TN) + ABS(H)
        ihit = ABS(TN-tcrit).LE.100.0D0*UROund*hmx
  !-----------------------------------------------------------------------
  ! Block G.
  ! The following block handles all successful returns from DLSOIBT.
  ! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
  ! ISTATE is set to 2, and the optional outputs are loaded into the
  ! work arrays before returning.
  !-----------------------------------------------------------------------
   1300 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        T = TN
        IF ( Itask.EQ.4 .OR. Itask.EQ.5 ) THEN
           IF ( ihit ) T = tcrit
        ENDIF
   1400 Istate = 2
        IF ( KFLag.EQ.-3 ) Istate = 3
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
  ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
   1500 ewti = Rwork(LEWt+i-1)
        msg = 'DLSOIBT- At T (=R1), EWT(I1) has become R2 .le. 0.'
        CALL XERRWD(msg,50,202,0,1,i,0,2,TN,ewti)
        Istate = -6
        GOTO 2100
  ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
   1600 msg = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
        CALL XERRWD(msg,50,204,0,0,0,0,0,0.0D0,0.0D0)
        msg = 'error test failed repeatedly or with ABS(H) = HMIN'
        CALL XERRWD(msg,50,204,0,0,0,0,2,TN,H)
        Istate = -4
        GOTO 1900
  ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
   1700 msg = 'DLSOIBT- At T (=R1) and step size H (=R2), the    '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      corrector convergence failed repeatedly     '
        CALL XERRWD(msg,50,205,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      or with ABS(H) = HMIN   '
        CALL XERRWD(msg,30,205,0,0,0,0,2,TN,H)
        Istate = -5
        GOTO 1900
  ! IRES = 3 returned by RES, despite retries by DSTODI.------------------
   1800 msg = 'DLSOIBT- At T (=R1) residual routine returned     '
        CALL XERRWD(msg,50,206,0,0,0,0,0,0.0D0,0.0D0)
        msg = '      error IRES = 3 repeatedly.        '
        CALL XERRWD(msg,40,206,0,0,0,0,1,TN,0.0D0)
        Istate = -7
        GOTO 2100
  ! Compute IMXER if relevant. -------------------------------------------
   1900 big = 0.0D0
        imxer = 1
        DO i = 1 , N
           size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
           IF ( big.LT.size ) THEN
              big = size
              imxer = i
           ENDIF
        ENDDO
        Iwork(16) = imxer
  ! Compute residual if relevant. ----------------------------------------
   2000 lyd0 = LYH + NYH
        DO i = 1 , N
           Rwork(i+LSAvf-1) = Rwork(i+lyd0-1)/H
           Y(i) = Rwork(i+LYH-1)
        ENDDO
        ires = 1
        CALL RES(Neq,TN,Y,Rwork(LSAvf),Ydoti,ires)
        NFE = NFE + 1
        IF ( ires.GT.1 ) THEN
           msg = 'DLSOIBT- Residual routine set its flag IRES       '
           CALL XERRWD(msg,50,210,0,0,0,0,0,0.0D0,0.0D0)
           msg = '      to (I1) when called for final output.       '
           CALL XERRWD(msg,50,210,0,1,ires,0,0,0.0D0,0.0D0)
        ENDIF
        GOTO 2200
  ! Set Y vector, T, and optional outputs. -------------------------------
   2100 DO i = 1 , N
           Y(i) = Rwork(i+LYH-1)
        ENDDO
   2200 T = TN
        Rwork(11) = HU
        Rwork(12) = H
        Rwork(13) = TN
        Iwork(11) = NST
        Iwork(12) = NFE
        Iwork(13) = NJE
        Iwork(14) = NQU
        Iwork(15) = NQ
        RETURN
   2300 msg = 'DLSOIBT- MF (=I1) illegal.    '
        CALL XERRWD(msg,30,8,0,1,Mf,0,0,0.0D0,0.0D0)
        GOTO 3100
   2400 msg = 'DLSOIBT- MB (=I1) or NB (=I2) illegal.  '
        CALL XERRWD(msg,40,9,0,2,mb,nb,0,0.0D0,0.0D0)
        GOTO 3100
   2500 msg = 'DLSOIBT- RTOL(=I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,19,0,1,i,0,1,rtoli,0.0D0)
        GOTO 3100
   2600 msg = 'DLSOIBT- ATOL(=I1) is R1 .lt. 0.0       '
        CALL XERRWD(msg,40,20,0,1,i,0,1,atoli,0.0D0)
        GOTO 3100
   2700 ewti = Rwork(LEWt+i-1)
        msg = 'DLSOIBT- EWT(I1) is R1 .le. 0.0         '
        CALL XERRWD(msg,40,21,0,1,i,0,1,ewti,0.0D0)
        GOTO 3100
   2800 msg =                                                             &
       &    'DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
        CALL XERRWD(msg,60,24,0,0,0,0,2,tcrit,TN)
        GOTO 3100
   2900 msg =                                                             &
       &    'DLSOIBT- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
        CALL XERRWD(msg,60,25,0,0,0,0,2,tcrit,Tout)
        GOTO 3100
   3000 msg = 'DLSOIBT- Trouble in DINTDY.  ITASK = I1, TOUT = R1'
        CALL XERRWD(msg,50,27,0,1,Itask,0,1,Tout,0.0D0)
  !
   3100 Istate = -3
        RETURN
  !----------------------- End of Subroutine DLSOIBT ---------------------
  99999 END

end module ode_lsoda_main
