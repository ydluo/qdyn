clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------

RESOLUTION = 500;


%randlow=0.5;
%randv=0.0;
%nasp=100;

year = 3600*24*365;
p = qdyn('set');

p.THETA_LAW=2;
p.FINITE=0;
p.SIGMA=0.5e6;
p.V_SS=1e-9;

p.V2=0.8e-7;

p.L=4e3;
p.W=110e3;
twm=2;
ts=2;
p.ACC = 1e-14;
%p.A=p.A*.8;
%p.B=p.B*.8;

filename =[ 'CO_Wtest_s.L',num2str(p.L/1e3),'sigma',num2str(p.SIGMA),'W',num2str(p.W/1e3),'res',num2str(RESOLUTION),'twm',num2str(twm),'ts',num2str(ts),'V2_',num2str(p.V2),'.mat']

%------------------------------
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
%------------------------------

p.N = 2^nextpow2(RESOLUTION*p.L/Lb); 
dx=p.L/p.N;
Lb_over_dx = Lb/dx


p = qdyn('set',p);



p.TMAX=twm*year;

%for i=1:1:floor(p.N*0.05)
%    p.V_0(i) = 1.01 *p.V_SS ;
%end
%for i=floor(p.N*0.05)+1:1:p.N
%    p.V_0(i)=p.V_SS;
%end
 p.V_0 = ( 1.+0.01*exp( -(p.X/(p.L/20)*2).^6 ) )*p.V_SS ;
 p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
 p.V_00=p.V_0;
%   for i=1:1:p.N
%       p.V_0(i) = p.V_00(mod((i+1024),p.N)+1);
%   end
  


p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;
[p,ot1,ox1]  = qdyn('run',p);
semilogy(ot1.t/year,ot1.v)
xlabel('Time (years)');
ylabel('Vmax');

  p.TMAX = ts*year;  
  p.NTOUT=1;

  p.V_0 = ox1.v(:,end);
  p.TH_0= ox1.th(:,end);
  %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
  %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
  [p,ot,ox]=qdyn('run',p);

save(filename)  


