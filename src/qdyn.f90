      implicit double precision (a-h,o-z)

      include 'qdyn.h'
! defines:
! nn 	= number of nodes
! finite = flag for boundary conditions along strike :
!          0 = periodic
!          1 = finite fault
! nnfft = array size for kernel computation = (finite+1)*nn
! nwfft = size of integer working array for fft = 2+ceil(sqrt(nfft/2)))

      parameter(neqs=2)
! neqs = # of equations to be satisfied simultaneously 
!           (2 for rate-and-state alone)

!******* WARNING: to force symmetry: 
!       (for some slip law long term cycle simulations)
      logical sym
      parameter(sym=.false.)

! theta = "standard" rate-and-state state
! mu_star,v_star,theta_star: Reference rate-and-state parameters
! tau = shear stress
! 
      real*8 :: tau(nn)
      real*8 kernel
      real*8 yt(neqs*nn),dydt(neqs*nn),yt_scale(neqs*nn)
      real*8 aa,bb,dc,vv1,vv2,mu_star,v_star,theta_star,velotau
      real*8 beta,smu
      real*8 zimpedance
      real*8 sigma,tau_init(nn)
      real*8 dx,dt_try,Lfault
      real*8 dt_did,dt_next,dt_max
      real*8 time,tmax,vslip(nn),theta(nn),slip(nn)
      real*8 Tper,Aper,Omper
      integer nn,nnfft,itheta_law
      integer NSTOP
      integer itstop
      real*8 lcold,lcnew,llocnew,llocold

!     working arrays for Ooura's fft
      integer iworkfft(0:nwfft-1)
      real*8 rworkfft(0:nnfft/2-1)

      EXTERNAL derivs
      common /props/ aa(nn),bb(nn),dc(nn),vv1(nn),vv2(nn), &
                     theta_star(nn),velotau(nnfft), &
                     kernel(nnfft),v_star,sigma,zimpedance, &
                     Aper,Omper,itheta_law  

! 
! Read input files
! 
! File qdyn.in:
!
! Line  |       Variables	Explanation
!	|
!  1	|	mu_star 	reference friction coefficient 
! 	|	v_star 		reference and loading slip rate
! 	|	theta_law	state variable evolution law
! 	|	              	 0 = ageing in the "no-healing" approximation
!       |       		 1 = ageing law
!       |       		 2 = slip law
! 	|	sigma  		normal stress
!	|
!  2	|	smu     	shear modulus (divided by 1-nu for mode II)
! 	|	beta          	shear wave velocity
! 	|	W             	length-scale in the "other" direction 
!   |   |                 (from 3D to 2D)
!       |       		e.g. depth of the seismogenic zone 
!       |       		(for pure 2D set W to a very large value >> Lfault)
!	|			External loading rate = pi*smu/W *v_star
!	|
!  3	|	Lfault        	fault length (periodic)
!   	| 	nstop		stopping criterion
! 	|	              	 0 = at tmax
!       |       		 1 = soon after end of slip localization 
!       |       		 2 = soon after maximum slip rate
!       |       		 3 = when reaching slip rate Vend
!	| 	tmax          	final time or Vend
!	|
!  4	| 	Tper          	period for time-periodic load (tides, etc)
!	| 	Aper          	amplitude of time-periodic loading
!	|
!  5	| 	nxout	        space-stride for snapshot outputs
!	| 	ntout	    	time-stride for snapshot outputs
!	|	ixout		location for timeseries output
!       |        		 0 = at fault center
!       |       		 1 = at location of maximum slip rate
!	|
!  6	| 	dt_try        	first trial timestep
!	| 	dt_max        	maximum timestep
!	| 	acc	    	accuracy for ODE solver
!	|
! 7-EOF	| 	x(i), a(i), b(i), dc(i),vv1(i),vv2(i), v0(i), theta0(i) 
!	|			position, rate-and-state parameters and initial conditions
!	|			for each fault node
!	|


!READ=====================0. Read input. START===========================
      open(unit=15,FILE = 'qdyn.in')
      read(15,*)mu_star,v_star,itheta_law,sigma
      read(15,*)smu,beta,W
      read(15,*)Lfault,NSTOP,tmax
      read(15,*)Tper,Aper
      read(15,*)nxout,ntout
      read(15,*)dt_try,dt_max,acc
      do i=1,nn
!       theta is ip-1, slip rate is ip
        ip=neqs*i
        read(15,*) xloc,aa(i),bb(i),dc(i),vv1(i),vv2(i),yt(ip),yt(ip-1)
      enddo
      close(15)
!READ=====================0. Read input. END=============================


!WR=======================1. Initialize. START===========================
!----------1.1.1 Initialize parameters pi,year,dx. START---------------   
      pi=4.d0*datan(1.d0)
      year = 3600.d0*24.d0*365.d0

      dx=Lfault/dble(nn)
!----------1.1.1 Initialize parameters pi,year,dx. END----------------- 

!----------1.1.2 Analysize dt_max & pertubation. START----------------- 
      if (Aper.ne.0.d0 .and. Tper.gt.0.d0) then
        if (dt_max.gt.0) then
          dt_max = min(dt_max,0.2d0*Tper)
        else
          dt_max = 0.2d0*Tper
        endif
      endif
      if (Tper.gt.0.d0) then
        Omper = 2*pi/Tper
      else
        Omper = 0.d0
      endif
!WRITE_6----------1.1.3 Output Impedance. START------------------------ 
      if (beta.gt.0d0) then 
        zimpedance= 0.5d0*smu/beta
      else
        zimpedance= 0.d0
      endif
      write(6,*)'impedance = ', zimpedance
!WRITE_6-----------1.1.3 Output Impedance. END--------------------------


!READ------1.2 Initialize elasto-static kernel. START------------------
      if (nn.gt.1) then
!- Periodic fault kernel
! In physical units: 
!   K(k) = mu/2 *sqrt(k^2 + kw^2)     where kw=2*pi/W
! In discretization, k=2*pi*i/(N*dx) for i=0:N/2
!   K(i) = pi*mu/(N*dx)*sqrt(i^2+(N*dx/W)^2) 
! The factor 2/N comes from the inverse FFT convention
        if (finite.eq.0) then
            tau_co=pi*smu/Lfault *2.d0/dble(nn)
            wl2 = (Lfault/W)**2
            do i=0,nn/2-1
            kernel(2*i+1) = tau_co*dsqrt(dble(i*i)+wl2)
            kernel(2*i+2) = kernel(2*i+1)
            enddo
            kernel(2) = tau_co*dsqrt(dble(nn*nn)/4.d0+wl2) ! Nyquist

!- Finite fault kernel
        else
!- Read coefficient I(n) from pre-calculated file.
            open(57,file=   &
           '~/2D_RUPTURE/STATIC/Matlab/kernel_I_32768.tab')
            if (nnfft/2>32768) stop 'Finite kernel table is too small'
                do i=1,nnfft/2-1
                    read(57,*) kernel(2*i+1)
                enddo
            read(57,*) kernel(2) ! Nyquist
            close(57)
!      make kernel
! The factor 2/N comes from the inverse FFT convention
            tau_co=pi*smu/(2d0*Lfault) *2.d0/dble(nnfft)
            kernel(1) = 0d0
            kernel(2) = tau_co*dble(nnfft/2)*kernel(2)
            do i=1,nnfft/2-1
                kernel(2*i+1) = tau_co*dble(i)*kernel(2*i+1)
                kernel(2*i+2) = kernel(2*i+1)
            enddo
        
        endif

     else
!       single degree-of-freedom spring-block system
!       critical size Lfault = Lc = mu*Dc/(b-a)*sigma
        kernel(1) = smu/Lfault   
     endif
!READ------1.2 Initialize elasto-static kernel. END--------------------

!WRITE_6---1.3 Output Initialized flag K/Kc. START---------------------
      i=nn/2
      if (i.eq.0) i=1
      if (FINITE.eq.1 .or. nn.eq.1) then
        write(6,*) 'K/Kc = ',(pi*smu/Lfault)/(sigma*(bb(i)-aa(i))/dc(i))
        write(6,*) 'K/Kb = ',(pi*smu/Lfault)/(sigma*bb(i)/dc(i))
      else
        write(6,*) 'K/Kc = ',(pi*smu/W)/(sigma*(bb(i)-aa(i))/dc(i))
        write(6,*) 'K/Kb = ',(pi*smu/W)/(sigma*bb(i)/dc(i))
      endif
      write(6,*)
!WRITE_6---1.3 Output Initialized flag K/Kc. END-----------------------


!----------1.4 Set initial fields START--------------------------------
      if (sym) call symmetrize(yt,nn)
      do i=1,nn
!       theta is ip-1, slip rate is ip
        ip=neqs*i
        theta_star(i)=dc(i)/vv2(i)
        tau_init(i)= sigma*( mu_star -aa(i)*log(vv1(i)/yt(ip)+1.d0) &
                            +bb(i)*log(yt(ip-1)/theta_star(i)+1.d0))
        tau(i)=tau_init(i)
        slip(i)=0.d0
      enddo
!----------1.4 Set initial fields END----------------------------------

!WRITE_ot---1.5.1 write ot file head START-----------------------------
      write(18,'(a)')'# macroscopic values:'
      write(18,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
      write(18,'(a)')'# values at center:'
      write(18,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
      write(18,'(a)')'# values at max(V) location:'
      write(18,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'
!WRITE_ot---1.5.1 write ot file head END-------------------------------
!WRITE_ox---1.5.2 write ox file head START-----------------------------
      nxout_count=0
      do i=1,nn,nxout
        nxout_count=nxout_count+1
      enddo
      write(19,'(a,2i5)')'# nx= ',nxout_count
!WRITE_ox---1.5.2 write ox file head START-----------------------------

!WR=======================1. Initialize. END=============================


!WC=======================2. Time loop. START============================
      write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'

      time=0.d0
      itstop=0
      lcnew=dble(nn)
      llocnew=dble(nn)
      vmax_old = 0d0
      vmax_older = 0d0

      ic=nn/2
      if (ic.eq.0) ic=1

      iworkfft(0) = 0

! Time loop

      it=0
      do
        if (sym) call symmetrize(yt,nn)

        it=it+1

        call derivs(time,yt,dydt)

! One step
        do i=1,nn*neqs
          yt_scale(i)=dabs(yt(i))+dabs(dt_try*dydt(i))
        enddo

!--------Call EXT routine bsstep [Bulirsch-Stoer Method] --------------
!-------- 
        call bsstep(yt,dydt,neqs*nn,time,dt_try,acc,yt_scale,   &
                    dt_did,dt_next,derivs)
        if (dt_max .gt. 0.d0) then
          dt_try=min(dt_next,dt_max)
        else
          dt_try=dt_next
        endif
!        theta is ip-1, slip rate is ip
        do i=1,nn
          ip=neqs*i
          vslip(i) = yt(ip)
          theta(i) = yt(ip-1)
        enddo

! Update slip, stress. 
! potency and potency rate
! crack size
        pot = 0d0
        pot_rate = 0d0
        do i=1,nn
          slip(i) = slip(i) +vslip(i)*dt_did
          pot = pot + slip(i)
          pot_rate = pot_rate + vslip(i)
          tau(i)  = ( mu_star -aa(i)*log(vv1(i)/vslip(i)+1.d0)  &
            +bb(i)*log(theta(i)/theta_star(i))+1.d0 ) *sigma
        enddo
        pot = pot*dx
        pot_rate = pot_rate*dx

        lcold = lcnew
        lcnew = crack_size(vslip,nn)

        llocold = llocnew
        llocnew = crack_size(velotau,nn)

! Output time series at max(v) location
        vmax=0.d0
        do i=1,nn
          if(vslip(i).gt.vmax)then
            vmax=vslip(i)
            ivmax=i
            xmax=(i-nn*0.5d0-0.5d0)*dx
          endif
        enddo

!WRITE_ot-------------- Export max values START------------------------
        write(18,'(e24.16,15e14.6)')time, llocnew*dx,lcnew*dx,  &
          pot, pot_rate,    &
          vslip(ic), theta(ic), vslip(ic)*theta(ic)/dc(ic), &
          tau(ic), slip(ic),    &
          xmax, vmax, theta(ivmax), vmax*theta(ivmax)/dc(ivmax),    &
          tau(ivmax), slip(ivmax)
!WRITE_ot-------------- Export max values END--------------------------


!WRITE_121-------------- STOP criteria START---------------------------
        if (itstop.eq.0) then
!         STOP soon after end of slip localization 
          if (NSTOP.eq.1) then
            if (llocnew.gt.llocold) itstop=it+2*ntout

!         STOP soon after maximum slip rate
          elseif (NSTOP.eq.2) then
          
            if (it.gt.2 .and. vmax_old.gt.vmax_older    &
               .and. vmax.lt.vmax_old) itstop=it+10*ntout
            vmax_older = vmax_old
            vmax_old = vmax
          
!         STOP at a slip rate threshold
          elseif (NSTOP.eq.3) then
            if (vmax.gt.tmax) itstop=it
          
!         STOP if time > tmax
          else
                  write(121,*) time, tmax
            if (tmax.gt.0d0 .and. time.gt.tmax) itstop=it
          endif
        endif
!WRITE_121-------------- STOP criteria END-----------------------------

!WRITE_6---------------- Print step to screen START--------------------
        if(mod(it-1,ntout).eq.0 .or. it.eq.itstop)then
! Print info
        write(6,'(i7,x,3(e11.3,x),i5)')it,dt_did,time/year,vmax
!WRITE_6---------------- Print step to screen END----------------------

!WRITE_ox--------------- Export snapshots START------------------------
        write(19,'(2a,2i5,e14.6)')'# x v theta',    &
       ' V./V dtau tau_dot slip ',it,ivmax,time
          do i=1,nn,nxout
            ip=neqs*i
            xloc=(i-nn*0.5d0-0.5d0)*dx
            write(19,'(8e15.7)')xloc,time,vslip(i), &
              theta(i),dydt(ip)/vslip(i),tau(i),     &
              velotau(i),slip(i)
          enddo
        endif
!WRITE_ox--------------- Export snapshots END------------------------
       

        if (it.eq.itstop) goto 333

      enddo

333   end


!C====================Subroutine derivs START===========================
!-------Compute thata,vslip and time deriv for one time step------------
!-------CALL rdft: real DFT---------------------------------------------
!---------------------------------------------------------------|
!-------|     in:time;yt:[theta,vslip];dydt             |-------|
!-------|    out:time;yt:[theta,vslip];dydt             |-------|
!-------| common:aa(nn),bb(nn),dc(nn),vv1(nn),vv2(nn),  |-------|
!       |        theta_star(nn),velotau(nnfft),         |       |
!       |        kernel(nnfft),v_star,sigma,zimpedance, |       |
!       |        Aper,Omper,itheta_law                  |       |
!-------|include:'qdyn.h':nn,finite,nnfft,nwfft         |-------|
!---------------------------------------------------------------|
 
      subroutine derivs(time,yt,dydt)
      implicit double precision (a-h,o-z)
      include 'qdyn.h'
!      parameter(nnfft=(finite+1)*nn)
      parameter(neqs=2)
!      parameter(nwfft=2+int(sqrt(nfft/2)))
      real*8 aa,bb,dc,v_star,theta_star
      real*8 kernel,dtau_per
      real*8 yt(*),dydt(*),omega
      real*8 time
      integer iworkfft(0:nwfft-1)
      real*8 rworkfft(0:nnfft/2-1)
      
      common /props/ aa(nn),bb(nn),dc(nn),vv1(nn),vv2(nn),  &
                     theta_star(nn),velotau(nnfft), &
                     kernel(nnfft),v_star,sigma,zimpedance, &
                     Aper,Omper,itheta_law

! compute shear stress rate from elastic interactions
      if (nn.gt.1) then
      do i=1,nn
        velotau(i)=v_star-yt(neqs*i)
      enddo
      do i=nn+1,nnfft 
        velotau(i)=0d0 
      enddo
      call rdft(nnfft,1,velotau,iworkfft,rworkfft)
      do i=1,nnfft
        velotau(i)=kernel(i)*velotau(i)
      enddo
      call rdft(nnfft,-1,velotau,iworkfft,rworkfft)

      else
        velotau(1)=kernel(1)*( v_star-yt(neqs) )
      endif

! periodic loading
      dtau_per = Omper *Aper*dcos(Omper*time)     

!
! theta is ip-1, vel is ip
!
! call unpack_yt
! call compute_dtheta_dt
! call compute_dv_dt  [this one contains computation of velotau]
! call pack_dydt

      do i=1,nn
        ip=neqs*i

!--------State evolution law START--------------------------------
        omega = yt(ip)*yt(ip-1)/dc(i)

!      "aging" law
        if (itheta_law.eq.1) then
          dydt(ip-1)=1.d0-omega

!      "slip" law
        elseif (itheta_law.eq.2) then
          if (omega.ne.0d0) then
            dydt(ip-1)=-omega*dlog(omega)
          else
            dydt(ip-1)=0d0
          endif

!      "aging" in the no-healing approximation
        elseif (itheta_law.eq.0) then
          dydt(ip-1)=-omega
        endif
!--------State evolution law END----------------------------------

! Shear stress
        dydt(ip)=(dtau_per+velotau(i) -sigma*bb(i)*vv2(i)*dydt(ip-1)   &
                  /(vv2(i)*yt(ip-1)+dc(i)) )    &
                  /( sigma*aa(i)*(1.d0/yt(ip)-1.d0/(vv1(i)+yt(ip)))   &
                     +zimpedance )

      enddo

      return
      end
!C====================Subroutine derivs END=============================




!======================Function crack_size START========================
! distance between largest peak on the left half
! and largest peak on the right half
      function crack_size(s,n)

      integer n
      real*8 s(n)

      real*8 smin,smax,s1,s2,xL,xR
      integer i,iL,iR

      if (n.eq.0) then
        crack_size = 0d0
        return
      endif

! assuming minimum stressing rate at nucleation point
!      smin=s(1)
!      do i=1,n
!        if (s(i).le.smin) then
!          imin=i
!          smin=s(i)
!        endif
!      enddo

! WARNING: assumes one peak on x<0 and one peak on x>0
      imin=n/2

      iL = 1
      smax=s(iL)
      do i=1,imin
        if (s(i).ge.smax) then
          iL=i
          smax=s(i)
        endif
      enddo
      xL = dble(iL);
      if (iL.gt.1) then
        s1 = 0.5d0*( s(iL+1)-s(iL-1) );
        s2 = s(iL-1)-2d0*s(iL)+s(iL+1);
        if (s2.ne.0d0) xL = xL-s1/s2;
      endif
      
      iR = imin
      smax=s(iR)
      do i=imin,n
        if (s(i).ge.smax) then
          iR=i
          smax=s(i)
        endif
      enddo
      xR = dble(iR);
      if (iR.lt.n) then
        s1 = 0.5d0*( s(iR+1)-s(iR-1) );
        s2 = s(iR-1)-2d0*s(iR)+s(iR+1);
        if (s2.ne.0d0) xR = xR-s1/s2;
      endif

      crack_size = xR-xL

      end
!======================Function crack_size END==========================


!===================Subroutine symmetrize START=========================
      subroutine symmetrize(yt,nn)

      parameter(neqs=2)

      integer nn
      real*8 yt(neqs*nn)
      
      integer i,ip,ipsym

      do i=1,nn/2
        ip=neqs*i
        ipsym=neqs*(nn+1-i)
        yt(ip-1) = 0.5d0*( yt(ip-1)+yt(ipsym-1) )
        yt(ipsym-1) = yt(ip-1)
        yt(ip) = 0.5d0*( yt(ip)+yt(ipsym) )
        yt(ipsym) = yt(ip)
      enddo

      end 
!===================Subroutine symmetrize END===========================
