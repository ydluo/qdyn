! Bulirsch-Stoer ODE solver from Numerical Recipes
module ode_bs

  private

  public :: bsstep

contains

      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,pb)

      use problem_class, only : problem_type

      implicit double precision (a-h,o-z)

      type(problem_type), intent(inout)  :: pb
      
      INTEGER nv,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,dydx(nv),y(nv),yscal(nv),          &
                SAFE1,SAFE2,REDMAX,REDMIN,TINY,SCALMX
      real*8 x,xnew
      PARAMETER (KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,      &
                SAFE2=.7d0,REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,    &
                SCALMX=.5d0) !SCALMX=.1d0
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,red,h,scale,work,wrkmin,xest,  &
                a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(nv),     &
                ysav(nv),yseq(nv)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/      &
            ((a(iq+1)-a(1)+1.d0)*(2*k+1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
!       if(xnew.eq.x)then
        if(x+1.e10*h.eq.x)then
          write(6,*)'dt_did,t= ',h,x
          pause 'step size underflow in bsstep'
        endif
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,pb)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,dabs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.d0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END SUBROUTINE bsstep


      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,pb)

      use problem_class, only : problem_type
      use derivs_all, only : derivs

      implicit double precision (a-h,o-z)

      type(problem_type), intent(inout)  :: pb

      INTEGER nstep,nvar
      REAL*8 htot,dydx(nvar),y(nvar),yout(nvar)
      real*8 xs,x
      INTEGER i,n
      REAL*8 h,h2,swap,ym(nvar),yn(nvar)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout,pb)
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout,pb)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END SUBROUTINE mmid


      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      implicit double precision (a-h,o-z)
      INTEGER iest,nv,IMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(nv),x(IMAX)
      REAL*8, allocatable :: qcol(:,:)
      SAVE qcol,x

      if (.not.allocated(qcol)) then
        allocate(qcol(nv,IMAX))
      elseif (size(qcol,1)/=nv) then
        deallocate(qcol)
        allocate(qcol(nv,IMAX))
      endif
        
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END SUBROUTINE pzextr

end module ode_bs

