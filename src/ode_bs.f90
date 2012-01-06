! Bulirsch-Stoer ODE solver from Numerical Recipes
module ode_bs

  public

contains

      SUBROUTINE bsstep(y,dydx,yscal,pb)
       
      use problem_class

      implicit double precision (a-h,o-z)

      type(problem_type), intent(inout)  :: pb
    
      INTEGER :: nv,NMAX,KMAXX,IMAX
      DOUBLE PRECISION ::  SAFE1,SAFE2,REDMAX,REDMIN,TINY,SCALMX 
      double precision, intent(inout) ::  & 
         y(pb%neqs*pb%mesh%nn),dydx(pb%neqs*pb%mesh%nn),yscal(pb%neqs*pb%mesh%nn)
               
      double precision :: xnew
!      double precision, intent(inout) :: x
      PARAMETER (NMAX=262144,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,      &
                SAFE2=.7d0,REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,    &
                SCALMX=.5d0) !SCALMX=.1d0
!     USES derivs,mmid,pzextr
      INTEGER :: i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      DOUBLE PRECISION :: eps1,epsold,errmax,fact,red,h,scale,work,wrkmin,xest,  &
                a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),     &
                ysav(NMAX),yseq(NMAX)
      LOGICAL ::first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
 !     EXTERNAL derivs
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(pb%acc .ne. epsold)then
        pb%dt_next=-1.d29
        xnew=-1.d29
        eps1=SAFE1*pb%acc
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
        epsold=pb%acc
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=pb%dt_try
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.pb%dt_next.or.pb%time.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=pb%time+h
!       if(xnew.eq.pb%time)then
        if(pb%time+1.e10*h.eq.pb%time)then
          write(6,*)'dt_did,t= ',h,pb%time
          pause 'step size underflow in bsstep'
        endif
        call mmid(ysav,dydx,h,nseq(k),yseq,pb)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,pb%neqs*pb%mesh%nn)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,dabs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/pb%acc
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
4     pb%time=xnew
      pb%dt_did=h
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
      pb%dt_next=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          pb%dt_next=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END SUBROUTINE bsstep


      SUBROUTINE mmid(y,dydx,htot,nstep,yout,pb)

      use derivs_all
      use problem_class

      implicit double precision (a-h,o-z)

      type(problem_type), intent(inout)  :: pb
      double precision :: htot 
      double precision, intent(inout) ::  & 
        y(pb%neqs*pb%mesh%nn),dydx(pb%neqs*pb%mesh%nn),yout(pb%neqs*pb%mesh%nn)
      INTEGER :: nstep,NMAX

      double precision :: xx,t_temp
!      double precision, intent(inout) :: xs
!      EXTERNAL derivs
      PARAMETER (NMAX=262144)
      INTEGER i,n
      DOUBLE PRECISION :: h,h2,swap,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,pb%neqs*pb%mesh%nn
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      xx=pb%time+h

!------save pb%time----------------------
      t_temp = pb%time  
      pb%time = xx
!------save pb%time----------------------
      call derivs(pb,yn,yout)
!------restore pb%time-------------------
      xx = pb%time
      pb%time = t_temp
!------restore pb%time-------------------

      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,pb%neqs*pb%mesh%nn
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        xx=xx+h

!------save pb%time----------------------
      t_temp = pb%time  
      pb%time = xx
!------save pb%time----------------------
      call derivs(pb,yn,yout)
!------restore pb%time-------------------
      xx = pb%time
      pb%time = t_temp
!------restore pb%time-------------------


13    continue
      do 14 i=1,pb%neqs*pb%mesh%nn
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END SUBROUTINE mmid


      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      implicit double precision (a-h,o-z)
      INTEGER iest,nv,IMAX,NMAX
      DOUBLE PRECISION :: xest
      double precision, intent(inout) ::  dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=262144)
      INTEGER :: j,k1
      DOUBLE PRECISION :: delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
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
