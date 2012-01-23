clear;
clc;

%------------------------------
%rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365;
p = qdyn('set');

p.MESHDIM=2;      %FFT enabled
p.THETA_LAW=1;
p.MU=30e9;
p.MU_SS=0.6;
p.SIGMA=100e6;
p.V_SS=0.1/year;
p.V2=0.01;              %no cut off velocity
p.DC=0.3;

dip0=14.;
dw0=2e3;

d1=10e3;     %upperbound of seismogenic zone (depth in m)
d11=15e3;    %limit of constant a/b
d2=35e3;     %lowerbound of seismogenic zone (depth in m)
d3=50e3;     %bottom of simulation zone (depth in m)
p.B=0.01;

p.W=d3/sin(dip0/180.*pi);
p.L=5.*p.W;
p.NW=ceil(p.W/dw0);
p.NX=1;

p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;     %deep to shallow
p.DIP_W(1:p.NW)=dip0;      %deep to shallow 
p.Z_CORNER=-d3+.5*p.DW(1)*sin(p.DIP_W(1)/180.*pi);   %p.Z_CORNER at center of left-bottom cell

dz0=dw0*sin(dip0/180.*pi);
ab0=0.4;     %a/b at seismogenic zone
abmax=5;     %a/b at d3
abs=1.5;     %a/b at shallow part
p.A(1:ceil((d3-d2)/dz0))=p.B*linspace(abmax,ab0,numel(1:ceil((d3-d2)/dz0)));   %increasing a/b below seismogenic zone
p.A(ceil((d3-d2)/dz0)+1:ceil((d3-d11)/dz0))=p.B*ab0;    %a/b < 1 in seismogenic zone
p.A(ceil((d3-d11)/dz0)+1:ceil((d3-d1)/dz0))=p.B*linspace(ab0,abs,numel(ceil((d3-d11)/dz0)+1:ceil((d3-d1)/dz0)));
p.A(ceil((d3-d1)/dz0)+1:p.NW)=p.B*abs;  %a/b >1 at shallow part  
  
twm=10000;         %warmup time in years
p.ACC = 1e-14;

%------------------------------
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
%------------------------------

filename = ['JP_2D','twm',num2str(twm),'L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'dip',num2str(dip0),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx


p = qdyn('set',p);

%Lc=Lb*(p.B/(p.B-p.A));
%disp(['  Lc=',num2str(Lc),'  L/Lc=',num2str(p.L/Lc),'  W/Lc=',num2str(p.W/Lc)]);
%Linf=2/pi*(p.B/(p.B-p.A))^2*Lb;
%disp(['  Linf=',num2str(Linf),'  L/Linf=',num2str(p.L/Linf),'  W/Linf=',num2str(p.W/Linf)]);

p.TMAX=twm*year;
p.V_0 = 1.01*p.V_SS ;


p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;

[p,ot1,ox1]  = qdyn('run',p);
semilogy(ot1.t/year,ot1.v)
xlabel('Time (years)');
ylabel('Vmax');
% 
%   p.TMAX = ts*year;  
%   p.NTOUT=1;
% 
%   p.V_0 = ox1.v(:,end);
%   p.TH_0= ox1.th(:,end);
%   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
%   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
%   [p,ot,ox]=qdyn('run',p);
%p0=p;
V_0 = ox1.v(:,end);
TH_0= ox1.th(:,end);
save(filename)  
save('warmup_jp.mat','p', 'V_0', 'TH_0');

