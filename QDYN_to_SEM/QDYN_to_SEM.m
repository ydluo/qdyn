clear;
clc;

rname='tohoku_l300_dc03.mat';
wtname='Tohoku.in';

load(rname);
disp(['Loaded QDYN output:',rname]);

hx=ceil(p.NX/2);
hw=ceil(18/50*p.NW);

rn=200;
rsn=3;
vmax=1.;
c_depth=3e3; % Cohesion cutoff depth in m

co=zeros(p.N,1);
swdc=zeros(p.N,1);
tau=zeros(p.N,1);
mu_s=zeros(p.N,1);
mu_d=zeros(p.N,1);
trr=ones(p.N,1)*1e9;
vrr=2e3;  %rupture prop speed in m/s
%hypocenter
ihx=ceil(p.NX/2);
ihw=p.NW-ceil(18/50*p.NW);
ih=ihx+(ihw-1)*p.NX;

rr=20e3; %radius of forced rupture zone in m

%----translate from r&s to sw
for iw=1:p.NW
    for ix=1:p.NX
        i=ix+(iw-1)*p.NX;
        dd=sqrt((p.X(i)-p.X(ih))^2+(p.Y(i)-p.Y(ih))^2+(p.Z(i)-p.Z(ih))^2);
        swdc(i)=p.DC(i)*log(vmax*p.TH_0(i)/p.DC(i));
        tau(i)=p.SIGMA(i)*(p.MU_SS+p.A(i)*log(p.V_0(i))+p.B(i)*log(p.TH_0(i)));
        mu_s(i)=tau(i)/p.SIGMA(i)+p.A(i)*log(vmax/p.V_0(i));
        mu_d(i)=tau(i)/p.SIGMA(i)+(p.A(i)-p.B(i))*log(vmax/p.V_0(i));
        if i >= p.N - p.NX*ceil(c_depth/p.DW(end))
            co(i)=4e6;
        end
        if dd <= rr
            trr(i)= dd/vrr;
        end
        
    end
end

load(rname);




fid=fopen(wtname,'w');

fprintf(fid,'%u %u  \n' , rn, rsn);    
fprintf(fid,'%u %u %E %E \n' ,p.NX, p.NW, p.L, p.W); 
fprintf(fid,'%u %u %E %E %E %E \n' ,ceil(p.NX/2),ceil(18/50*p.NW),p.L/2.,18./50.*p.W,4e3,2e3);  
disp(['Start write SEM input:...']);
for iw=p.NW:-1:1
    for ix=1:p.NX
        id=ix+(iw-1)*p.NX;
        fprintf(fid,'%u %u %E %E %E %E %E %E %E %E %E %E %E %E \n' , ...
            ix-1,p.NW-iw,p.X(id)-p.X(1),p.W-sqrt((p.Y(id))^2+(p.Z(id)-p.Z_CORNER)^2),...
            p.SIGMA(id),0.,tau(id),0.,tau(id)/p.SIGMA(id),mu_s(id),mu_d(id),...
            swdc(id),co(id),trr(id));
    end
end
disp(['Generated SEM input:', wtname]);


fclose(fid);