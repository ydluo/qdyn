clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------

DC_mean = 0.03;		%mean DC
dcsigma = 0.25;		%Dv
DC_min = 0.02;		%min_DC

DC_v = (exp((dcsigma)^2)-1)*DC_mean^2;
dcmu = log(DC_mean^2/sqrt(DC_v+DC_mean^2));


year = 3600*24*365;
p = qdyn('set');
p.NSTOP = 0;        %stop at v_th
p.TMAX = 100000*year;       % stop at v = v_th = tmax


co = 4e6;  %cohesion
co_limit = 3e3;  %first X m to appy cohesion 

p.RNS_LAW=0;
p.MESHDIM=2;      %FFT enabled
p.THETA_LAW=1;
p.MU=40e9;
p.LAM=40e9;
p.MU_SS=0.6;
%p.SIGMA=100e6;
p.V_SS=0.01/year;
p.V2=100.;              %no cut off velocity
p.V1=p.V2;

p.OX_SEQ=1;
p.OX_DYN=1;
p.DYN_TH_ON = 0.1;
p.DYN_TH_OFF = 0.1;
%p.DC=0.3;

 

dip0=90.;
dw0=0.5e3/2;
dd=20e3;    %depth of L sigma change
d1=5e3;     %upperbound of seismogenic zone (depth in m)
d2=9e3;    %limit of constant b/a
d3=27e3;     %lowerbound of seismogenic zone (depth in m)
d4=31e3;     %limit of constant b/a
db=50e3;     %bottom of simulation zone (depth in m)
aa0=0.01;    %p.A
sigma0=75e6;        %sigma max
nxout=1;       %snapshot output grid interval
p.NXOUT_DYN=1;   %dynamic snapshot output grid interval
p.W=db/sin(dip0/180.*pi);
p.L=512e3;
p.NW=ceil(p.W/dw0);
p.NX=1024*2;
p.VS=3000.;

p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;     %deep to shallow
p.DIP_W(1:p.NW)=dip0;      %deep to shallow 
%p.Z_CORNER=-db+.5*p.DW(1)*sin(p.DIP_W(1)/180.*pi);   %p.Z_CORNER at center of left-bottom cell
p.Z_CORNER=-db;

p.IC = p.N/2;

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

p.DC(1:ceil((db-dd)/dz0))=0.3;
p.DC(ceil((db-dd)/dz0)+1:p.NW)=0.3;

p.SIGMA(1:ceil((db-dd)/dz0))=sigma0;
p.SIGMA(ceil((db-dd)/dz0)+1:p.NW)=linspace(sigma0,1e6,numel(ceil((db-dd)/dz0)+1:p.NW));

p.CO(1:ceil((db-co_limit)/dz0))=0;
p.CO(ceil((db-co_limit)/dz0)+1:p.NW)=co;

p.N=p.NX*p.NW;
tmp_A=p.A;
tmp_B=p.B;
tmp_SIGMA=p.SIGMA;
tmp_DC=p.DC;
tmp_CO=p.CO;


for i=1:p.NW
    p.X((i-1)*p.NX+1:i*p.NX) = linspace(0,p.L,p.NX);
    p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i);
    p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);
    p.CO((i-1)*p.NX+1:i*p.NX) = tmp_CO(i);
	
end


for i=1:1:p.N
    p.IOT(i) = 0;
end

p.DC = lognrnd(dcmu,dcsigma,1,p.N);
p.DC = max(p.DC,DC_min);

twm=100000;         %warmup time in years
ts=1000;    %simulation time in years
p.ACC = 1e-10;
Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
%p.DYN_TH_ON=Vdyn/10.;
%p.DYN_TH_OFF=Vdyn/10.;

%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B)
%Lnuc = 1.3774*Lb;
%------------------------------

filename = ['Hete_3D_ss_twm',num2str(twm),'L',num2str(p.L/1000.),...
    'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),...
    'dip',num2str(dip0),'DCmean',num2str(DC_mean),'DcV',num2str(dcsigma),'.mat']
p.IC=ceil(p.N/2);
dw=p.W/p.NW;
Lb_over_dw = Lb/dw
dx=p.L/p.NX;
Lb_over_dx = Lb/dx

p = qdyn('set',p);

p.V_0 = 1.01*p.V_SS ;


p.NTOUT=1000;
p.NXOUT=nxout;

% p.DYN_FLAG=1;
% p.DYN_M=10.^19;
% p.DYN_SKIP = 1;

[p,ot0,ox0]  = qdyn('run',p);
semilogy(ot0.t/year,ot0.v)
xlabel('Time (years)');
ylabel('Vmax');
% 
   p.NTOUT=10;
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
save('warmup_jp_s.mat','p', 'V_0', 'TH_0');

