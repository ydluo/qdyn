clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365;
p = qdyn('set');

p.MESHDIM=2;
p.THETA_LAW=2;

p.SIGMA=0.5e6;
p.V_SS=1e-9;


p.V2=0.8e-7;

p.L=100e3;
p.W=100e3;
p.NX=4;
p.NW=4;
p.Z_CORNER=-50e3;
p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;
p.DIP_W(1:p.NW)=30.0;
twm=20;
ts=2;
p.ACC = 1e-14;

%------------------------------
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
%------------------------------

filename = 'test_2d_1d_approx.mat'
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx


p = qdyn('set',p);



p.TMAX=twm*year;

%for i=1:1:floor(p.N*0.05)
%    p.V_0(i) = 1.01 *p.V_SS ;
%end
%for i=floor(p.N*0.05)+1:1:p.N
%    p.V_0(i)=p.V_SS;
%end
 p.V_0 = 1.01*p.V_SS ;
% p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
 p.V_00=p.V_0;
%   for i=1:1:p.N
%       p.V_0(i) = p.V_00(mod((i+1024),p.N)+1);
%   end
  


p.NTOUT=10;
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

save(filename)  


