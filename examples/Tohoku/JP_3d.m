
year = 3600*24*365;
load('warmup_jp');

p.L = 1000e3; % along-strike rupture length
p.NX=512; % must be a power of 2 such that dx=L/NX is smaller than Lb
p.TMAX=2000*year;         %warmup time in years

p.NTOUT=400;
p.NXOUT=1;
p.NSTOP=0;
p.ACC = 1e-10; % default is 1e-14. Increase it (less accurate) to speed up.

p.N=p.NX*p.NW;

% in 2D size(A)=[1 NW]
% in 3D, faster index runs along strike, size(A)=[NX NW]
p.A = repmat(p.A(:)',p.NX,1);  
p.B = repmat(p.B(:)',p.NX,1);  
p.SIGMA = repmat(p.SIGMA(:)',p.NX,1);  
p.DC = repmat(p.DC(:)',p.NX,1);  
p.V_0 = repmat(p.V_0(:)',p.NX,1);  
p.TH_0 = repmat(p.TH_0(:)',p.NX,1);  
% qdyn.m handles these quantities correctly even if they are matrices
% If in doubt, to make sure they are stored as vectors we could do this: p.A=p.A(:);

Lb = p.MU.*p.DC./p.SIGMA./p.B;
Lnuc = 1.3774*Lb;

filename = ['JP_3D','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat'];
disp(['Output file = ',filename]) 
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
dw=p.W/p.NW;
disp(['min Lb/dx = ',num2str(min(Lb(:)/dx)) ])
disp(['min Lb/dw = ',num2str(min(Lb(:)/dw)) ])

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

