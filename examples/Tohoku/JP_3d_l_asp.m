%clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365;
load('warmup_jp');
p.OX_DYN = 1;
p.OX_SEQ = 1;
p.NX=512;
p.RNS_LAW=0;
p.N=p.NX*p.NW;
p.L=400e3;
p.ACC=1e-10;
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
    p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i);
    p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);

    p.V_0((i-1)*p.NX+1:i*p.NX) = V_0(i);
    p.TH_0((i-1)*p.NX+1:i*p.NX) = TH_0(i);
end


%----asperity spacing control
sp_int=25;   %spacing of asperities
scatter_mode=4;     % = 1 square; 
                    % = 2 triangular
                    % = 3 random
                    % = 4 Customized
bd_left=0.;       %asperities boundary
bd_right=0.;
bd_low=0.;
bd_up=0.;


asp_single = 1;

%----asperity property
l_asp0=40.0*1e3;          %asperity size lower limit in m
l_asp1=40.0*1e3;          %asperity size upper limit in m
ba_asp=1.5;       %b/a of asperity
dc_asp=0.002;       %Dc of asperity
sigma_asp=88.2e6*2;        %sigma(asp)




twm=20000;         %warmup time in years



p.IOT=zeros(size(p.X));
p.IASP=zeros(size(p.X));

%-----set asperity location
i_asp=zeros(size(p.X));
switch scatter_mode
    case 1
    disp(['Asperities scattering: Square      Spacing:', num2str(sp_int)]);
    for ix=1:p.NX
     for iw=1:p.NW
         in=ix+(iw-1)*p.NX;
         if mod(ix,sp_int) == 1 && mod(iw,sp_int) == 1
             i_asp(in)=1;
         end
     end
    end

    case 2
    sp_int_2 = round(sqrt((sp_int)^2-(sp_int/2)^2));  
    disp(['Asperities scattering: Triangular      Spacing:', num2str(sp_int)]);    
    for ix=1:p.NX
     for iw=1:p.NW
         in=ix+(iw-1)*p.NX;
         if mod(ix,sp_int) == 1 && mod(iw,sp_int_2*2) == 1
             i_asp(in)=1;
         end
         if mod(ix,sp_int) == 1+floor(sp_int/2) && mod(iw,sp_int_2*2) == 1+sp_int_2
             i_asp(in)=1;
         end
     end
    end
     
    case 3
    pth = 1/(pi*sp_int^2);
    disp(['Asperities scattering: Random      Spacing:', num2str(sp_int)]);    
    for in=1:p.N
        if rand(1,1) <= pth
            i_asp(in)=1;
        end
    end

    case 4
    disp(['Asperities scattering: Customized ']);    

%     i_asp(63489+round(160/(400/512)))=1;
%     i_asp(64000-round(160/(400/512)))=1;

    i_asp((63489+10):40:(64000-10))=1;
     
    
end    

for ix=1:p.NX
     for iw=1:p.NW
         in=ix+(iw-1)*p.NX;
         if ix<=p.NX*bd_left || ix>=p.NX*(1-bd_right) || iw<=p.NW*bd_low || iw>=p.NW*(1-bd_up)
             i_asp(in)=0;
         end
     end
end

asp_count_all=0;
for i=1:p.N
    if i_asp(i) == 1
        asp_count_all=asp_count_all+1;
    end
end

asp_count=0;
%----set asperity property

if asp_single == 0;
    for i=1:p.N
        if i_asp(i) == 1
            p.IASP(i) = 1;
            p.IOT(i) = 1;
            asp_count=asp_count+1;
            l_asp=l_asp0+(l_asp1-l_asp0)*rand(1);
            disp(['Setting asperity: ',num2str(asp_count),'/',num2str(asp_count_all)]);
            for j=1:1:p.N
                dd=sqrt((p.X(j)-p.X(i))^2+(p.Y(j)-p.Y(i))^2+(p.Z(j)-p.Z(i))^2);
                p.B(j)=p.B(j)+(-p.B(j)+p.A(j)*ba_asp)*exp(-(dd/l_asp*2).^6);
                p.DC(j)=p.DC(j)+(-p.DC(j)+dc_asp)*exp(-(dd/l_asp*2).^6);
                p.SIGMA(j)=p.SIGMA(j)+(-p.SIGMA(j)+sigma_asp)*exp(-(dd/l_asp*2).^6);
            end
            Lc_asp=p.MU*p.DC(i)/(p.SIGMA(i)*(p.B(i)-p.A(i)));
            disp(['  Normalized size L/Lc = ',num2str(l_asp/Lc_asp)]);
        end
    end
end

if  asp_single == 1;
  for i=1:p.N
        if i_asp(i) == 1
            p.IASP(i) = 1;
            p.IOT(i) = 1;
            asp_count=asp_count+1;
            disp(['Setting single_cell asperity: ',num2str(asp_count),'/',num2str(asp_count_all)]);
            p.B(i)=+p.A(i)*ba_asp;
            p.DC(i)=dc_asp;
            p.SIGMA(j)=p.SIGMA(j)+(-p.SIGMA(j)+sigma_asp)*exp(-(dd/l_asp*2).^6);
            Lc_asp=p.MU*p.DC(i)/(p.SIGMA(i)*(p.B(i)-p.A(i)));
            disp(['  Normalized size Lx/Lc = ',num2str(dx/Lc_asp), '    Lw/Lc = ', num2str(dw/Lc_asp)]);
        end
  end
end

Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
disp(['Vdyn = ' num2str(Vdyn)]);
p.DYN_TH_ON=Vdyn/10.;
p.DYN_TH_OFF=Vdyn/10.;

%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B);
Lnuc = 1.3774*Lb;
%------------------------------



filename = ['JP_3D','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx
dw=p.W/p.NW;
Lb_over_dw = Lb/dw


p.TMAX=twm*year;
p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;
p.DYN_FLAG=0;
p.DYN_M=10.^19.5;
p.DYN_SKIP = 1;

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

