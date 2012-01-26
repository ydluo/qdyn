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
p.V_SS=0.085/year;
p.V2=0.01;              %no cut off velocity
%p.DC=0.3;

dip0=14.;
dw0=2e3;
dd=20e3;    %depth of L sigma change
d1=5e3;     %upperbound of seismogenic zone (depth in m)
d2=10e3;    %limit of constant b/a
d3=32e3;     %lowerbound of seismogenic zone (depth in m)
d4=33e3;     %limit of constant b/a
db=50e3;     %bottom of simulation zone (depth in m)
aa0=0.01;    %p.A

p.W=db/sin(dip0/180.*pi);
p.L=5.*p.W;
p.NW=ceil(p.W/dw0);
p.NX=1;

p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;     %deep to shallow
p.DIP_W(1:p.NW)=dip0;      %deep to shallow 
p.Z_CORNER=-db+.5*p.DW(1)*sin(p.DIP_W(1)/180.*pi);   %p.Z_CORNER at center of left-bottom cell


p.A(1:p.NW)=aa0;
dz0=dw0*sin(dip0/180.*pi);
ba0=1.5;     %b/a at seismogenic zone
%abmax=5;     %a/b at d3
bam=0.6;     %b/a at shallow/deeper part
 
p.B(1:ceil((db-d4)/dz0))=p.A(1:ceil((db-d4)/dz0)).*bam;
p.B(ceil((db-d4)/dz0)+1:ceil((db-d3)/dz0))=p.A(ceil((db-d4)/dz0)+1:ceil((db-d3)/dz0)).*linspace(bam,ba0,numel(ceil((db-d4)/dz0)+1:ceil((db-d3)/dz0)));   %increasing a/b below seismogenic zone
p.B(ceil((db-d3)/dz0)+1:ceil((db-d2)/dz0))=p.A(ceil((db-d3)/dz0)+1:ceil((db-d2)/dz0)).*ba0;    %a/b < 1 in seismogenic zone
p.B(ceil((db-d2)/dz0)+1:ceil((db-d1)/dz0))=p.A(ceil((db-d2)/dz0)+1:ceil((db-d1)/dz0)).*linspace(ba0,bam,numel(ceil((db-d2)/dz0)+1:ceil((db-d1)/dz0)));
p.B(ceil((db-d1)/dz0)+1:p.NW)=p.A(ceil((db-d1)/dz0)+1:p.NW).*bam;  %a/b >1 at shallow part 

p.DC(1:ceil((db-dd)/dz0))=0.02;
p.DC(ceil((db-dd)/dz0)+1:p.NW)=0.5;

p.SIGMA(1:ceil((db-dd)/dz0))=88.2e6;
p.SIGMA(ceil((db-dd)/dz0)+1:p.NW)=linspace(350e6,0,numel(ceil((db-dd)/dz0)+1:p.NW));

  
twm=10000;         %warmup time in years
ts=500;    %simulation time in years
p.ACC = 1e-14;

%------------------------------
%Lb = p.MU*p.DC/p.SIGMA/p.B;
%Lnuc = 1.3774*Lb;
%------------------------------

filename = ['JP_2D_s_twm',num2str(twm),'ts',num2str(ts),'L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'dip',num2str(dip0),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
%Lb_over_dx = Lb/dx


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

[p,ot0,ox0]  = qdyn('run',p);
semilogy(ot0.t/year,ot0.v)
xlabel('Time (years)');
ylabel('Vmax');
% 
   p.TMAX = ts*year;  
   p.NTOUT=1;
% 
   p.V_0 = ox0.v(:,end);
   p.TH_0= ox0.th(:,end);
%   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
%   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
   [p,ot1,ox1]=qdyn('run',p);
%p0=p;
V_0 = ox1.v(:,end);
TH_0= ox1.th(:,end);
save(filename)  
save('warmup_jp.mat','p', 'V_0', 'TH_0');

