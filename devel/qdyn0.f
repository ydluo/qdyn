      parameter(nn=65536,neqs=4) 
c
c nn = # of nodes, neqs # of equations to be satisfied simultaneously (2 for rate-and-state alone)
c state1 = "standard" rate-and-state state
c state_sn = Prakash-Clifton "normal stress" state (for bimaterial problem)
c mu_star,v_star,state1_star: Reference rate-and-state parameters
c sn = normal stress
c tau = shear stress
c beta,rlrm,bmfactor: Bimaterial parameters
c
c rns3 tries not-steady state in center (constant stress) to try to stabilize
c straying nucleation zone in bimaterial case
c
c rns4 is rns3 but w/out velocity cutoff and no radiation damping
c rnsloc comes from rns4 but writes useful stuff for assessing localization phase.
c nucleate randomizes theta(0), starts with uniform V=V*, and V.theta/Dc<=1
c
      implicit double precision (a-h,o-z)
      real*8 state1(nn),state_sn(nn),sn(nn),vel(nn)
      real*8 tau(nn)
      real*8 wavenum,ddt_tauinf
      real*8 yt(neqs*nn),dydt(neqs*nn),yt_scale(neqs*nn)
      real*8 aa,bb,d_c,t_sn,mu_star,v_star,state1_star
      real*8 nu_1,nu_2,sm_1,sm_2,chi_1,chi_2,beta,rlrm,bmfactor
      real*8 beta_1,beta_2
      real*8 rnnm1,tau_co
      real*8 sn_init,tau_init(nn)
      real*8 dx,dt_try
      real*8 dt_did,dt_next,acc
      real*8 velsv(nn),thetasv(nn),taudotsv(nn)
      real*8 time,slip(nn)
      real rnum,ran1
      integer nn,itmax
      logical final
      logical nuc,firstnuc
      EXTERNAL derivs
      common /props/ aa,bb,d_c,t_sn,mu_star,v_star,state1_star,
     &               ddt_tauinf(nn),wavenum(nn),tau_co,bmfactor,dx,
     &               sm_1,sm_2,beta_1,beta_2,sn_init,velotau(nn)

      final=.false.  !when to end (optional)
      nuc=.false.
      firstnuc=.false.
      pi=4.d0*datan(1.d0)
c
c Initialize coefficients for inverse transform.  Watch the factor of 2!
c
      rnnm1=2.d0/dreal(nn)
      wavenum(1)=0.d0
      wavenum(2)=1.d0  
      do i=2,nn/2
        wavenum(2*i-1)=(i-1)*rnnm1
        wavenum(2*i)=wavenum(2*i-1)
      enddo
c
c Initialize.
c (Can't change neqs only, but can look for lines containing it, to change)
c
      open(unit=11,file='inparams')
      read(11,*)aa,bb,d_c,t_sn
      read(11,*)mu_star,v_star
      state1_star=d_c/v_star
      read(11,*)mu_prime
      read(11,*)dx,dt_try,sn_init
      close(11)
c
c mu' mu_inplane = mu_antiplane/(1-nu)
      tau_co=pi*mu_prime/(nn*dx)
      write(6,*)tau_co
c
c Watch! 0.5
c
      bmfactor=0.5d0*rlrm*beta
      write(6,*)'bmfactor,rlrm ',bmfactor,rlrm
      time=0.d0
      abpower=bb/aa
      call get_ddt_tauinf(nn,estnucl,ddt_tauinf,dx)
      iseed=-3
      do i=1,nn
c       read(10,*)state1(i),state_sn(i),vel(i)
        vel(i)=v_star
c       rnum=0.5+0.01*(ran1(iseed)-0.5)
        rnum=ran1(iseed)
        state1(i)=(d_c/vel(i))*rnum
c       state1(i)=d_c/vel(i)
        state_sn(i)=0.d0
        dist=dabs(dx*(0.5d0*nn+0.5d0-i))
c
        sn(i)=0.d0
        tau_init(i)=(mu_star+aa*log(vel(i)/v_star)
     &         +bb*log(state1(i)/state1_star))*(state_sn(i)+sn_init)
        tau(i)=tau_init(i)
        yt(neqs*i-3)=state1(i)
        yt(neqs*i-2)=state_sn(i)
        yt(neqs*i-1)=sn(i)
        yt(neqs*i)=vel(i)
        slip(i)=0.d0
      enddo
      acc=1.d-7
      write(6,*)'estnucl ',estnucl
      write(18,'(3a)')'> > > > > > > > > >  time,V(ivmax),V.dot(ivmax),'
     &      ,'theta(ivmax),theta.dot(ivmax),V*theta(ivmax),'
     &      ,',tau.dot^el,tau,slip,xinuclen'
c
c Time loop
c
c     WRITE(6,*)'timestep,slip0,v0dot,v0,state1dot,state1,tau0'
      inuclenmin=nn
      xinuclenmin=nn*dx
      ivmaxsv=9999
      itmax=40000
      do it=1,itmax
        call derivs(time,yt,dydt)
c
c Define scaling 
c
        do j=1,nn
          jp=neqs*j
          yt_scale(jp-3)=dabs(yt(jp-3))+dabs(dt_try*dydt(jp-3))
          yt_scale(jp-2)=dabs(yt(jp-2))+1.d6+dabs(dt_try*dydt(jp-2))
          yt_scale(jp-1)=dabs(yt(jp-1))+1.d6+dabs(dt_try*dydt(jp-1))
          yt_scale(jp)=dabs(yt(jp))+dabs(dt_try*dydt(jp))
        enddo
        call bsstep(yt,dydt,neqs*nn,time,dt_try,acc,yt_scale,
     &              dt_did,dt_next,derivs)
        dt_try=dt_next
c
c stats:
c
        vavlog=0.d0
c       taumax=0.d0
        vmax=0.d0
        do i=1,nn
          ip=neqs*i
          slip(i)=slip(i)+yt(ip)*dt_did
          tau(i)=(mu_star+aa*log(yt(ip)/v_star)
     &         +bb*log(yt(ip-3)/state1_star))*(yt(ip-2)+sn_init)
          vavlog=vavlog+log(yt(ip))/2.303d0
          if(yt(ip).gt.vmax)then
            vmax=yt(ip)
            ivmax=i
          endif
        enddo
        vavlog=vavlog/dreal(nn)
c       ivmax=3584
        ipmax=neqs*ivmax
c       vmax=yt(ipmax)
        tauscale=velotau(ivmax)
c 3584 b01a009 tdot02 dx0.125(2x) V01e-9 16384
c 26295 b01a009 tdot02 dx0.25(2x)
c 3549 ?
c 3583 b01a003t tdot02
c 3583 b01a006 tdot02
c 3591 b01a006 tdot03 less random
c 3554 b01a009 tdot02 
c 3565 b01a005t tdot02
c 3551 b01a007t tdot02
c 6095 b01a008 tdot02 less random 
c 3613 b01a008 tdot02 random 
        if(mod(it,500).eq.2)then
        write(19,'(2a,2i5,2e14.6)')'> > > > > > > x vel vel/vel0 theta',
     & ' norm-x V./V dtau dtau0 tau_dot slip ',it,ivmax,time,tauscale
c         do i=55168-10240,55168+10240
          do i=1,65536,32
c           if(j.le.0)then
c             i=j+nn
c           elseif(j.ge.nn+1)then
c             i=j-nn
c           else
c             i=j
c           endif
            ip=neqs*i
            xloc=(i-nn/2-0.5d0)*dx
            write(19,'(10e15.7)')xloc,yt(ip),yt(ip)/yt(ipmax),
     &        yt(ip-3),xloc/tauloc1,dydt(ip)/yt(ip),tau(i)-tau_init(i),
     &        time*ddt_tauinf(i),-velotau(i)/tauscale,slip(i)
          enddo
        endif
        IF(.not.(nuc))THEN
          taudotmin=tauscale
          do i=ivmax,nn-1
            if(velotau(i).lt.taudotmin)then
              taudotmin=velotau(i)
              itaudotminp=i
              twocurvdx2=velotau(i-1)-2.d0*velotau(i)+velotau(i+1)
              xoff=(velotau(i-1)-velotau(i+1))/(2.d0*twocurvdx2)
              xitaudotminp=itaudotminp+xoff
            endif
          enddo
c         if(itaudotminp.eq.nn-1)pause 'bad random start'
          taudotmin=tauscale
          do i=ivmax,2,-1
            if(velotau(i).lt.taudotmin)then
              taudotmin=velotau(i)
              itaudotminm=i
              twocurvdx2=velotau(i-1)-2.d0*velotau(i)+velotau(i+1)
              xoff=(velotau(i-1)-velotau(i+1))/(2.d0*twocurvdx2)
              xitaudotminm=itaudotminm+xoff
            endif
          enddo
c         if(itaudotminm.eq.2)pause 'bad random start'
          xinuclen=0.5d0*(xitaudotminp-xitaudotminm)*dx
          inuclen=itaudotminp-itaudotminm
c For loc point poke
c         if (vel0/velinf.gt.5.e2)then
c For boxes
c         if (vel0/velinf.gt.1.e4)then
c For random
          if (log(vmax)/2.303d0-vavlog.gt.(2.d0*bb/aa)-1.2222d0)then
            if (xinuclen.le.xinuclenmin)then
              xinuclenmin=xinuclen
              do i=itaudotminm,itaudotminp
                velsv(i)=yt(neqs*i)
                thetasv(i)=yt(neqs*i-1)
                taudotsv(i)=velotau(i)
                tauscalesv=tauscale
                itaudotminmsv=itaudotminm
                itaudotminpsv=itaudotminp
                itsv=it
                ivmaxsv=ivmax
                vmaxsv=vmax
                vavsv=vavlog
              enddo
            endif
            if((xinuclen.ge.xinuclenmin).and.(it.ge.itsv+40).and.
     &         (.not.firstnuc))then
              firstnuc=.true.
              write(29,'(2a,3e14.6)')'> > > > > > >  x vel vel/vel0 ',
     &        'theta junk junk junk junk tau_dot junk ',time,aa,bb
              do i=itaudotminmsv,itaudotminpsv
                xloc=(i-nn/2-0.5d0)*dx
                write(29,'(10e15.7)')xloc,velsv(i),velsv(i)/vmaxsv,
     &                thetasv(i),99.999,99.999,99.999,
     &                99.999,-taudotsv(i)/tauscalesv,99.999
              enddo
c             nuc=.true.
c             final=.true.
            endif
          endif
        ENDIF
c       write(18,'(e24.16,10e14.6)')0.3063329738568574D+09-time,
        write(18,'(e24.16,10e14.6)')time,
     &    yt(ipmax),dydt(ipmax),yt(ipmax-3),dydt(ipmax-3),
     &    yt(ipmax)*yt(ipmax-3),
     &    velotau(ivmax),tau(ivmax),slip(ivmax),xinuclen
        write(6,*)it,dt_did,time,vmax,ivmax,xinuclen
c       IF((log(vmax)/2.303d0-vavlog.gt.12d0).or.(dt_did.lt.1d-10).or.
c    &     (it.eq.itmax).or.((it.gt.1000d0).and.
c    &     (log(vmax)/2.303d0-vavlog.lt.1.d0)))
c    &    THEN
c         inuc=inuclenmin
c         write(29,'(2a,3e14.6)')'> > > > > > >  x vel vel/vel0 ',
c    &    'theta norm-x dsn dtau dtau0 tau_dot slip ',time,aa,bb
c         do i=itaudotminm,itaudotminp
c           ip=neqs*i
c           xloc=(i-nn/2-0.5d0)*dx
c           write(29,'(10e15.7)')xloc,velsv(i),velsv(i)/vmaxsv,
c    &                thetasv(i),99.999,99.999,99.999,
c    &                99.999,-taudotsv(i)/tauscalesv,99.999
c         enddo
c         final=.true.
c       ENDIF
c       if(final)goto 333
        if(vmax.gt.1)goto 333
      enddo

333   end

c-----------------------------------------------------
      subroutine get_ddt_tauinf(nn,estnucl,ddt_tauinf,dx)
      implicit double precision (a-h,o-z)
      real*8 ddt_tauinf(nn)
      coeff=0.5*1.d-2
c     coeff=0.5*1.e-3
      pi=4.d0*datan(1.d0)
      do i=1,nn
c       dist=dabs(dx*(0.5d0*nn+0.5d0-i))
          ddt_tauinf(i)=2.d0*coeff
      enddo

      return
      end

c------------------------------------------------------
      subroutine derivs(time,yt,dydt)
      implicit double precision (a-h,o-z)
      parameter(nn=65536,neqs=4)
      common /props/ aa,bb,d_c,t_sn,mu_star,v_star,state1_star,
     &               ddt_tauinf(nn),wavenum(nn),tau_co,bmfactor,dx,
     &               sm_1,sm_2,beta_1,beta_2,sn_init,velotau(nn)
      real*8 aa,bb,d_c,t_sn,mu_star,v_star,state1_star
      real*8 velo(nn),fmu(nn),ddt_tauinf,wavenum
      real*8 yt(*),dydt(*)
      real*8 time
c
c state1 is ip-3, state_sn is ip-2, sn is ip-1, vel is ip
c
      do i=1,nn
        velo(i)=yt(neqs*i)
      enddo
      call drealft(velo,nn,1)
      do i=1,nn
        velotau(i)=tau_co*wavenum(i)*velo(i)
      enddo
      call drealft(velotau,nn,-1)
      do i=1,nn
        ip=neqs*i
c Original ("ageing")
        dydt(ip-3)=1.d0-yt(ip)*yt(ip-3)/d_c
c       dydt(ip-3)=-yt(ip)*yt(ip-3)/d_c
c Alternate ("slip")
c       dydt(ip-3)=-(yt(ip)*yt(ip-3)/d_c)*log(yt(ip)*yt(ip-3)/d_c)
cc       dydt(ip-2)=(yt(ip-2)-yt(ip-1))/t_sn
        dydt(ip-2)=10.d0*(yt(ip-1)-yt(ip-2))*yt(ip)/d_c
        fmu(i)=mu_star+aa*log(yt(ip)/v_star)
     &         +bb*log(yt(ip-3)/state1_star)
c with elasticity
        dydt(ip)=(ddt_tauinf(i) - velotau(i) - dydt(ip-2)*fmu(i)
     &         -(yt(ip-2)+sn_init)*bb*dydt(ip-3)/yt(ip-3))
     &         /((yt(ip-2)+sn_init)*aa/yt(ip))
c without elasticity
c       dydt(ip)=(ddt_tauinf(i) - dydt(ip-2)*fmu(i)
c    &         -(yt(ip-2)+sn_init)*bb*dydt(ip-3)/yt(ip-3))
c    &         /((yt(ip-2)+sn_init)*aa/yt(ip))
c other?
c       dydt(ip)=(ddt_tauinf(i) - velotau(i) - dydt(ip-2)*fmu(i)
c    &            -yt(ip-2)*bb*dydt(ip-3)/(yt(ip-3)+state1_star))
c    &           /-yt(ip-2)*aa*(1./(v_star+yt(ip))-1./yt(ip))
      enddo
      do i=2,nn-1
        dydt(neqs*i-1)=0.5d0*(yt(neqs*(i+1))-yt(neqs*(i-1)))*bmfactor/dx
      enddo
      dydt(neqs*1-1)=0.5d0*(yt(neqs*2)-yt(neqs*nn))*bmfactor/dx
      dydt(neqs*nn-1)=0.5d0*(yt(neqs*1)-yt(neqs*(nn-1)))*bmfactor/dx

      return
      end

c------------------------------
      function ran1(idum)
      integer idum,ia,im,iq,ir,ntab,ndiv
      real ran1,am,eps,rnmx
      parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     &          ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
c
c From Press et al.  Returns a random 0<#<1.  Call with idum a negative
c integer to initialize; thereafter leave idum unchanged. rnmx should
c approximate the largest floating value that is <1.
c
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if ((idum.le.0).or.(iy.eq.0))then
        idum=max(-idum,1)
          do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
        enddo
          iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      return
      end

