clear;


%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365;
p = qdyn('set');

p.MESHDIM=2;
p.THETA_LAW=2;
p.FAULT_TYPE = 2;
p.SIGMA=0.5e6;
p.V_SS=1e-9;

p.A=0.003;
p.B=0.01;
p.SIGMA_CPL=1;

p.V2=0.01;

p.L=80e3;
p.W=80e3;
p.NX=64;
p.NW=64;
p.Z_CORNER=-100e3;
p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;
p.DIP_W(1:p.NW)=30.0;
twm=3.0;
p.ACC = 1e-14;

%------------------------------
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
%------------------------------

filename = ['test_2d_fft_ab',num2str(p.A/p.B),'L',num2str(p.L/1000),'nx',num2str(p.NX),'W',num2str(p.W/1000),'nw',num2str(p.NW),'z',num2str(p.Z_CORNER/1000),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx


p = qdyn('set',p);

Lc=Lb*(p.B/(p.B-p.A));
disp(['  Lc=',num2str(Lc),'  L/Lc=',num2str(p.L/Lc),'  W/Lc=',num2str(p.W/Lc)]);
Linf=2/pi*(p.B/(p.B-p.A))^2*Lb;
disp(['  Linf=',num2str(Linf),'  L/Linf=',num2str(p.L/Linf),'  W/Linf=',num2str(p.W/Linf)]);

p.TMAX=twm*year;

p.V_0 = 1.01*p.V_SS ;

p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;
p.OX_DYN = 1;
p.OX_SEQ = 1;

[p,ots,oxs]  = qdyn('run',p);
save 3d_fftserial

