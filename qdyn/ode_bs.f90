! Bulirsch-Stoer ODE solver from Numerical Recipes
module ode_bs

  use logger, only : log_screen

  implicit none
  private

  public :: bsstep

contains

      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,pb,ik)

      use problem_class, only : problem_type
      use my_mpi, only: is_MPI_parallel, max_allproc

      integer, intent(in) :: nv
      double precision, intent(inout) :: y(nv), dydx(nv), x
      double precision, intent(in) :: yscal(nv), htry, eps
      double precision, intent(out) :: hdid, hnext
      type(problem_type), intent(inout) :: pb
                          !NOTE: inout needed by FFT in compute_stress

      integer, parameter :: KMAXX=8, IMAX=KMAXX+1
      double precision, parameter :: SAFE1=.25d0, SAFE2=.7d0, &
                                     REDMAX=1.d-5, REDMIN=.7d0, &
                                     TINY=1.d-30,    &
                                     SCALMX=.5d0 !SCALMX=.1d0
      integer :: iq,k,km,ik
      integer, dimension(IMAX) :: nseq = (/ 2,4,6,8,10,12,14,16,18 /)
      integer, save :: kmax,kopt
      double precision, save :: alf(KMAXX,KMAXX)
      double precision :: err(KMAXX)
      double precision, save :: a(IMAX)
      double precision, save :: epsold = -1.d0,xnew
      double precision :: eps1,errmax,err_km,fact,red,h,scale,wrkmin,xest,  &
                yerr(nv),ysav(nv),yseq(nv),errmaxglob
      logical, save :: first = .true.
      logical :: reduct
      character(255) :: msg

      red = 0d0 ! SEISMIC: initialise red to be safe

      if (eps /= epsold) then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
        enddo
        do iq=2,KMAXX
          do k=1,iq-1
            alf(k,iq)=eps1**( (a(k+1)-a(iq+1))/ ((a(iq+1)-a(1)+1.d0)*(2*k+1)) )
          enddo
        enddo
        epsold=eps
        do kopt=2,KMAXX-1
          if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
        enddo
        kmax=kopt
      endif
      h=htry
      ysav=y
      if (h /= hnext .or. x /= xnew) then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
!PG
    ik=0
    main_loop: do

      do k=1,kmax
        xnew=x+h
!       if (xnew == x) then
        if (x+1.e10*h == x) then  !NOTE: different than original code (above). Why???
          write(msg, *) 'dt_did,t= ', h, x
          call log_screen(msg)
          stop 'step size underflow in bsstep'
!JPA: make it a safe MPI stop: all procs should stop if at least one has an error
        endif
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,pb)
        ik=ik+nseq(k)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if (k == 1) cycle
        errmax=maxval(dabs(yerr/yscal))
        if (is_MPI_parallel()) then
          call max_allproc(errmax,errmaxglob)
          errmax=errmaxglob
        endif
        errmax=max(TINY,errmax)/eps
        km=k-1
        err_km=(errmax/SAFE1)**(1.d0/(2*km+1))
        err(km)=err_km
        if (k >= kopt-1 .or. first) then
          if (errmax < 1.0) exit main_loop
          if (k == kmax .or. k == kopt+1) then
            red=SAFE2/err_km
            exit
          else if (k == kopt) then
            if (alf(kopt-1,kopt) < err_km) then
              red=1.d0/err_km
              exit
            endif
          else if (kopt == kmax) then
            if (alf(km,kmax-1) < err_km) then
              red=alf(km,kmax-1)*SAFE2/err_km
              exit
            endif
          else if (alf(km,kopt) < err_km) then
            red=alf(km,kopt-1)/err_km
            exit
          endif
        endif
      enddo
      red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.

    enddo main_loop
!PG
      x=xnew
      hdid=h
      first=.false.
      kopt=1+minloc(a(2:km+1)*max(err(1:km),SCALMX), 1)
      scale=max(err(kopt-1),SCALMX)
      wrkmin=scale*a(kopt)
      hnext=h/scale
      if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if (a(kopt+1)*fact <= wrkmin) then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
!JPA if everything worked correctly, hdid and hnext should be the same in all processors

      END SUBROUTINE bsstep


      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,pb)

      use problem_class, only : problem_type
      use derivs_all, only : derivs

      integer, intent(in) :: nstep,nvar
      double precision, intent(in) :: y(nvar),dydx(nvar),xs,htot
      double precision, intent(out) :: yout(nvar)
      type(problem_type), intent(inout)  :: pb

      double precision :: x, h,h2,swap(nvar),ym(nvar),yn(nvar)
      integer :: n

      h=htot/nstep
      ym=y
      yn=y+h*dydx
      x=xs+h
      call derivs(x,yn,yout,pb)
      h2=2.d0*h
      do n=2,nstep
        swap=ym+h2*yout
        ym=yn
        yn=swap
        x=x+h
        call derivs(x,yn,yout,pb)
      enddo
      yout=0.5d0*(ym+yn+h*yout)

      END SUBROUTINE mmid


      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)

      integer, intent(in) :: iest,nv
      double precision, intent(in) :: xest,yest(nv)
      double precision, intent(out) :: dy(nv),yz(nv)

      integer, parameter :: IMAX=13
      integer :: k1
      double precision :: delta,q(nv),d(nv)
      double precision, save :: x(IMAX)
      double precision, allocatable, save :: qcol(:,:)

      if (.not.allocated(qcol)) then
        allocate(qcol(nv,IMAX))
      elseif (size(qcol,1)/=nv) then
        deallocate(qcol)
        allocate(qcol(nv,IMAX))
      endif

      x(iest)=xest
      dy=yest
      yz=yest
      if (iest == 1) then
        qcol(:,1)=yest
      else
        d=yest
        do k1=1,iest-1
          q=qcol(:,k1)
          qcol(:,k1)=dy
          delta=1.d0/(x(iest-k1)-xest)
          d=delta*(d-q)
          dy=xest*d
          d=x(iest-k1)*d
          yz=yz+dy
        enddo
        qcol(:,iest)=dy
      endif

      END SUBROUTINE pzextr

end module ode_bs
