
year = 3600*24*365;
load('warmup_jp');

p.L = 300e3; % along-strike rupture length
p.NX=128;

p.N=p.NX*p.NW;

tmp_A=p.A;
tmp_B=p.B;
tmp_SIGMA=p.SIGMA;
tmp_DC=p.DC;

for i=1:p.NW
    p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i);
    p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);

    p.V_0((i-1)*p.NX+1:i*p.NX) = V_0(i);
    p.TH_0((i-1)*p.NX+1:i*p.NX) = TH_0(i);
end

twm=4000;         %warmup time in years


%------------------------------
Lb = p.MU.*p.DC./p.SIGMA./p.B;
Lnuc = 1.3774*Lb;
%------------------------------

filename = ['JP_3D','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
dw=p.W/p.NW;
min_Lb_over_dx = min(Lb/dx)
min_Lb_over_dw = min(Lb/dw)



p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;

return
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

