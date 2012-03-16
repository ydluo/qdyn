clear;
clc;

%------------------------------
%rand(1,floor(sum(100*clock)));
%------------------------------

isp=215;
isp2=223;
year = 3600*24*365;
load('JP_3D_dc03_tltL300nx256W206.6783nw207.mat');

p.V_0=ox1.v(:,isp);
p.TH_0=ox1.th(:,isp);


twm=ox1.t(isp2)-ox1.t(isp);


%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B);
Lnuc = 1.3774*Lb;
%------------------------------

filename = ['JP_3D_dc03_tlt_hd','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx



p.NTOUT=10;
p.NXOUT=1;
p.NSTOP=0;
p.TMAX=twm;
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

V_0 = ox1.v(:,end);
TH_0= ox1.th(:,end);
save(filename)  
save('warmup_jp_3d.mat','p', 'V_0', 'TH_0');

