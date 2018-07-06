clear;
clc;

%------------------------------
%rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365;
load('warmup_jp');
p.FAULT_TYPE = 2;

p.NX=512;

p.N=p.NX*p.NW;
p.L=400e3;
p.ACC=1e-10;
ctlt=0.05;
tmp_A=p.A;
tmp_B=p.B;
tmp_SIGMA=p.SIGMA;
tmp_DC=p.DC;

tmp_Y=p.Y;
tmp_Z=p.Z;

for i=1:p.NW
    p.X((i-1)*p.NX+1:i*p.NX) = linspace(0,p.L,p.NX);
    p.Y((i-1)*p.NX+1:i*p.NX) = tmp_Y(i);
    p.Z((i-1)*p.NX+1:i*p.NX) = tmp_Z(i);
    p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i)*linspace(1+ctlt,1-ctlt,p.NX);
    p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);

    p.V_0((i-1)*p.NX+1:i*p.NX) = V_0(i);
    p.TH_0((i-1)*p.NX+1:i*p.NX) = TH_0(i);
end



twm=3000;         %warmup time in years


%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B);
Lnuc = 1.3774*Lb;
%------------------------------
p.TMAX = twm*year;
filename = ['JP_3D_dc03_tlt','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx



p.NTOUT=1000;
p.NXOUT=1;
p.NSTOP=0;
p.TMAX=twm*year;
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


