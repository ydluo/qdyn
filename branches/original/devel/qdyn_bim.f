      implicit double precision (a-h,o-z)
      parameter(nn=128,neqs=4,acc=1.d-7,itmax=40000)
c nn   = # of nodes 
c neqs = # of equations to be satisfied simultaneously (2 for rate-and-state alone)
c acc  = accuracy for ODE solver
c itmax= # of time steps
c
c state1 = "standard" rate-and-state state
c state_sn = Prakash-Clifton "normal stress" state (for bimaterial problem)
c mu_star,v_star,state1_star: Reference rate-and-state parameters
c sn  = normal stress
c tau = shear stress
c beta,rlrm,bmfactor: Bimaterial parameters
c
      real*8 state1(nn),state_sn(nn),sn(nn),vel(nn)
      real*8 tau(nn)
      real*8 kernel,ddt_tauinf
      real*8 yt(neqs*nn),dydt(neqs*nn),yt_scale(neqs*nn)
      real*8 aa,bb,d_c,t_sn,mu_star,v_star,state1_star,velotau
      real*8 aa_r,bb_r,d_c_r
      real*8 nu_1,nu_2,sm_1,sm_2,chi_1,chi_2,beta,rlrm,bmfactor
      real*8 beta_1,beta_2
      real*8 zimpedance
      real*8 sn_init,tau_init(nn)
      real*8 dx,dt_try
      real*8 dt_did,dt_next
      real*8 time,slip(nn)
      integer nn,itmax
      logical abdc_from_file
      EXTERNAL derivs
      common /props/ aa(nn),bb(nn),d_c(nn),state1_star(nn),velotau(nn),
     &               kernel(nn),t_sn,mu_star,v_star,
     &               ddt_tauinf,bmfactor,sn_init,zimpedance

c
c Initialize.
c
      open(unit=15,FILE = 'inparams')
      read(15,*)aa_r,bb_r,d_c_r,t_sn
      read(15,*)mu_star,v_star
      read(15,*)sm_1,sm_2,nu_1,nu_2,beta_1,beta_2
      read(15,*)dx,dt_try,sn_init
      read(15,*)abdc_from_file
      read(15,*)W
      close(15)

      if (abdc_from_file) then
        open(100,file='abdc.in')
        do i=1,nn
          read(100,*) xloc,aa(i),bb(i),d_c(i)
        enddo
        close(100)
      else
        do i=1,nn
          aa(i) = aa_r
          bb(i) = bb_r
          d_c(i) = d_c_r
        enddo
      endif

c
c The below is from Rice,Lapusta,Ranjith 2001, p. 1887.  Should reduce to 
c beta=0, rlrm=sm/(1-nu) for identical solids.
c
      chi_1=3.d0-4.d0*nu_1
      chi_2=3.d0-4.d0*nu_2
      beta=(sm_2*(chi_1-1.d0)-sm_1*(chi_2-1.d0))/
     &     (sm_2*(chi_1+1.d0)+sm_1*(chi_2+1.d0))
      rlrm=(8.d0*sm_1*sm_2/(1.d0-beta*beta))/
     &     (sm_2*(chi_1+1.d0)+sm_1*(chi_2+1.d0))
      zimpedance=1.d0/( beta_1/sm_1 + beta_2/sm_2 )
      write(6,*)'impedance', zimpedance
c Watch! 0.5
      bmfactor=0.5d0*rlrm*beta
      write(6,*)'bmfactor,rlrm ',bmfactor,rlrm
      bmfactor=0.5d0*bmfactor/dx

c Initialize elasto-static kernel.
c In physical units: 
c   K(k) = mu/2 *sqrt(k^2 + kw^2)     where kw=2*pi/W
c In discretization, k=2*pi*i/(N*dx) for i=0:N/2
c   K(i) = pi*mu/(N*dx)*sqrt(i^2+(N*dx/W)^2) 
c A factor 2/N comes from the Numerical Recipes inverse FFT (drealft)
      pi=4.d0*datan(1.d0)
      tau_co=pi*rlrm/(dreal(nn)*dx) *2.d0/dreal(nn)
      wl2 = (dreal(nn)*dx/W)**2
      do i=0,nn/2-1
        kernel(2*i+1) = tau_co*dsqrt(dreal(i*i)+wl2)
        kernel(2*i+2) = kernel(2*i+1)
      enddo
      kernel(2) = tau_co*dsqrt(dreal(nn*nn)/4.d0+wl2)

c Set external loading rate
      ddt_tauinf = pi*rlrm/W *v_star
      write(6,*) 'K/Kc = ',(pi*rlrm/W) / (sn_init*(bb_r-aa_r)/d_c_r)

c Set initial fields
      do i=1,nn

c        dist=dabs(dx*(0.5d0*nn+0.5d0-i))
        vel(i)=v_star*2.d0
        state1(i)=d_c(i)/v_star

        state1_star(i)=d_c(i)/v_star
        state_sn(i)=0.d0
        sn(i)=0.d0
        tau_init(i)=(mu_star+aa(i)*log(vel(i)/v_star)
     &     +bb(i)*log(state1(i)/state1_star(i)))*(state_sn(i)+sn_init)
        tau(i)=tau_init(i)
        slip(i)=0.d0

        ip=neqs*i
        yt(ip-3)=state1(i)
        yt(ip-2)=state_sn(i)
        yt(ip-1)=sn(i)
        yt(ip)  =vel(i)

      enddo

      write(18,'(a)')'# values at max(V) location'
      write(18,'(2a)')'# 1=time, 2=V, 3=V.dot, 4=theta, 5=theta.dot, ',
     &      '6=V*theta, 7=tau.dot, 8=tau, 9=slip'
c
c Time loop
c
      time=0.d0

      do it=1,itmax
        call derivs(time,yt,dydt)
c
c Define scaling 
c
        do i=1,nn*neqs
          yt_scale(i)=dabs(yt(i))+dabs(dt_try*dydt(i))
        enddo
        do i=1,nn
          ip=neqs*i
          yt_scale(ip-2)=yt_scale(ip-2)+1.d6
          yt_scale(ip-1)=yt_scale(ip-1)+1.d6
        enddo

c One step
        call bsstep(yt,dydt,neqs*nn,time,dt_try,acc,yt_scale,
     &              dt_did,dt_next,derivs)
        dt_try=dt_next

c Update slip, stress. Get max(v)
c state1 is ip-3, state_sn is ip-2, sn is ip-1, vel is ip
        vmax=0.d0
        do i=1,nn
          ip=neqs*i
          slip(i) = slip(i) +yt(ip)*dt_did
          tau(i)  = ( mu_star +aa(i)*log(yt(ip)/v_star)
     &      +bb(i)*log(yt(ip-3)/state1_star(i)) ) *(yt(ip-2)+sn_init)
          if(yt(ip).gt.vmax)then
            vmax=yt(ip)
            ivmax=i
          endif
        enddo
        ipmax=neqs*ivmax
        tauscale=velotau(ivmax)

c Export max values
        write(18,'(e24.16,8e14.6)')time,
     &    yt(ipmax),dydt(ipmax),yt(ipmax-3),dydt(ipmax-3),
     &    yt(ipmax)*yt(ipmax-3),
     &    velotau(ivmax),tau(ivmax),slip(ivmax)
 
c Export snapshots
        if(mod(it,500).eq.2)then
        write(19,'(2a,2i5,2e14.6)')'# x v v/vmax theta',
     & ' V./V dtau dtau0 tau_dot slip ',it,ivmax,time,tauscale
          do i=1,nn,32
            ip=neqs*i
            xloc=(i-nn/2-0.5d0)*dx
            write(19,'(9e15.7)')xloc,yt(ip),yt(ip)/yt(ipmax),
     &        yt(ip-3),dydt(ip)/yt(ip),tau(i)-tau_init(i),
     &        time*ddt_tauinf,-velotau(i)/tauscale,slip(i)
          enddo
        endif

c Print info
        write(6,*)it,dt_did,time,vmax,ivmax

        if(vmax.gt.1)goto 333
      enddo

333   end


c------------------------------------------------------
      subroutine derivs(time,yt,dydt)
      implicit double precision (a-h,o-z)
      parameter(nn=128,neqs=4)
      common /props/ aa(nn),bb(nn),d_c(nn),state1_star(nn),velotau(nn),
     &               kernel(nn),t_sn,mu_star,v_star,
     &               ddt_tauinf,bmfactor,sn_init,zimpedance
      real*8 aa,bb,d_c,t_sn,mu_star,v_star,state1_star
      real*8 fmu,ddt_tauinf,kernel
      real*8 yt(*),dydt(*)
      real*8 time

c compute shear stress rate from elastic interactions
      do i=1,nn
        velotau(i)=yt(neqs*i)
      enddo
      call drealft(velotau,nn,1)
      do i=1,nn
        velotau(i)=kernel(i)*velotau(i)
      enddo
      call drealft(velotau,nn,-1)

c
c state1 is ip-3, state_sn is ip-2, sn is ip-1, vel is ip
c
      do i=1,nn
        ip=neqs*i

c State evolution law 
c "ageing" law
        dydt(ip-3)=1.d0-yt(ip)*yt(ip-3)/d_c(i)
c "ageing" in the no-healing approximation
c       dydt(ip-3)=-yt(ip)*yt(ip-3)/d_c(i)
c "slip" law
c       dydt(ip-3)=-(yt(ip)*yt(ip-3)/d_c(i))*log(yt(ip)*yt(ip-3)/d_c(i))

c Normal state law
cc       dydt(ip-2)=(yt(ip-2)-yt(ip-1))/t_sn
        dydt(ip-2)=10.d0*(yt(ip-1)-yt(ip-2))*yt(ip)/d_c(i)

c Shear stress
        fmu = mu_star +aa(i)*log(yt(ip)/v_star)
     &         +bb(i)*log(yt(ip-3)/state1_star(i))
        dydt(ip)= ( ddt_tauinf - velotau(i) - dydt(ip-2)*fmu
     &         -(yt(ip-2)+sn_init)*bb(i)*dydt(ip-3)/yt(ip-3) )
     &         /( (yt(ip-2)+sn_init)*aa(i)/yt(ip) +zimpedance )

      enddo

c Normal stress, periodic boundaries
      dydt(neqs*1-1)=(yt(neqs*2)-yt(neqs*nn))*bmfactor
      do i=2,nn-1
        dydt(neqs*i-1)=(yt(neqs*(i+1))-yt(neqs*(i-1)))*bmfactor
      enddo
      dydt(neqs*nn-1)=(yt(neqs*1)-yt(neqs*(nn-1)))*bmfactor

      return
      end

