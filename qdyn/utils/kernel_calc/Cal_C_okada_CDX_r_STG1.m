clear


filename='CDX_L1_test_RES.mat';

fid = fopen([filename '.txt'],'w');	%rect crack model
fid_r = fopen([filename '_r.txt'],'w');		%circ. model
fids = fopen([filename 's.txt'],'w');		%rect model w/ shallow VS
fid_rs = fopen([filename '_rs.txt'],'w');	%circ. model w/ shallow VS

fidDL = fopen([filename 'DL.txt'],'w');		%rect dislocation model
fid_rDL = fopen([filename 'DL_r.txt'],'w');	%circ. dislocation model
fidDLc = fopen([filename 'DLc.txt'],'w');         %rect dislocation model, dtau at center
fid_rDLc = fopen([filename 'DLc_r.txt'],'w');     %circ. dislocation model

fprintf(fid,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fids,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rs,'C    W    L_a    Ws    Zc    A    RES\n');

fprintf(fidDL,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rDL,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fidDLc,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rDLc,'C    W    L_a    Ws    Zc    A    RES\n');


%LL = [2e3:2e3:20e3,100e3,200e3,1000e3];
%LL = [1e3];
LL = 1e3*[0.1:0.1:10,10.2:0.2:20,20.1:0.1:21.5,21.55:0.05:22];

%LL = [22e3];

DsVS = 4e3;		%depth of shallow VS zone

%RES_s = [1:1:20];
%RES_s = [1:1:20 22:2:30 35:5:50 60 70 80 100];
RES_s = [50];

mu = 40e9;
lam = 40e9;
Ws = 22e3;
DIP = 90;
RES = 5;       %(2*RES-1)^2 points for a square rupture
ZZc0 = -11e3;       %starting center Z of rupture 

nn_all = numel(LL)*numel(RES_s);
L_a = zeros(nn_all,1);      %actual L
C = zeros(nn_all,1);
Cr = C;
Cs = C;
Crs = Cr;
CDL = C;
CrDL = C;
CDLc = C;
CrDLc = C;

ii = 0;

for  iL = 1:1:numel(LL)
    
    L = LL(iL);
    
    
    for iRes = 1:1:numel(RES_s)
        
    ii = ii+1;
    
    RES = RES_s(iRes);
    
    display(['Calculating C from crack model [Okada]: #' num2str(ii) ' of ' num2str(nn_all)]);
    
    if L<=Ws
        W = L;
        L_a(ii) = L;
        NW = 2*RES-1;
        NX = NW;
        N = NX*NW;
        dx = L/NX;
        dw = dx;
        x0 = (0.5:1:NX)*dx;
        X = repmat(x0,1,NW);
        Zc = ZZc0 + (-Ws/2-ZZc0)*(L/Ws);     %center Z of rupture
        z0 = (-RES+1:1:RES-1)*dw+Zc;
        Z = reshape(repmat(z0,NX,1),1,N);
        Y = zeros(size(X));
        display(['Square rupture ' num2str(L/1000)  'km*' num2str(L/1000) 'km | Resolution = ' num2str(RES)])
    end
    
    if L > Ws
        W = Ws;
        NW = 2*RES-1;
        dw = Ws/NW;
        dx = dw;
        NX = round(L/dx);
        L_a(ii) = NX*dx;
        N = NX*NW;
        x0 = (0.5:1:NX)*dx;
        X = repmat(x0,1,NW);        
        z0 = (-RES+1:1:RES-1)*dw-Ws/2;
        Z = reshape(repmat(z0,NX,1),1,N);
        Y = zeros(size(X));
        display(['Rectangular rupture ' num2str(W/1000)  'km*' num2str(L_a(ii)/1000) 'km | Resolution = ' num2str(RES)])
      
    end
        
    disp('Compute kernel')
    K = qdyn_okada_kernel([NX,NW],mu,lam,X,Y,Z,DIP,dx,dw);
    
    display('Calculating C value ...');    
    tau = ones(size(X'));
    D = K\tau;
    C(ii) = mean(tau)*W/(mean(D)*mu);
    display(['C = ' num2str(C(ii))]);
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',C(ii),W,L_a(ii),Ws,Zc,numel(X)*dx*dw,RES);

    tau(Z>= -DsVS) = 0;
    Ds = K\tau;
    Cs(ii) = mean(tau)*W/(mean(Ds)*mu);
    display(['Cs = ' num2str(Cs(ii)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fids,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cs(ii),W,L_a(ii),Ws,Zc,numel(X)*dx*dw,RES);
 
    DDL = ones(size(X));
    tauDL = DDL*K;
    CDL(ii) = mean(tauDL)*W/mu;
    display(['CDL = ' num2str(CDL(ii)) ' | Dislocation Model']);
    fprintf(fidDL,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CDL(ii),W,L_a(ii),Ws,Zc,numel(X)*dx*dw,RES);

    CDLc(ii) = tauDL(ceil(numel(tauDL)/2))*W/mu;
    display(['CDLc = ' num2str(CDLc(ii)) ' | Dislocation Model | Dtau at Center']);
    fprintf(fidDLc,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CDLc(ii),W,L_a(ii),Ws,Zc,numel(X)*dx*dw,RES);

    disp('Generating Full Kernel for Circular Rupture');
    Xc = X(RES);
    IIr = find(((X-Xc).^2+(Z-Zc).^2)<=(L/2)^2);
    Xr = X(IIr);
    Yr = Y(IIr);
    Zr = Z(IIr);
    Kr = K(IIr,IIr);
    display('Calcalation Cr value :...');    
    taur = ones(size(Xr'));
    Dr = Kr\taur;
    Cr(ii) = mean(taur)*W/(mean(Dr)*mu);
    display(['Cr = ' num2str(Cr(ii))]);
    fprintf(fid_r,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cr(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    taur(Zr>= -DsVS) = 0;    
    Drs = Kr\taur;
    Crs(ii) = mean(taur)*W/(mean(Drs)*mu);
    display(['Crs = ' num2str(Crs(ii)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fid_rs,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Crs(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    DrDL = ones(size(Xr));
    taurDL = DrDL*Kr;
    CrDL(ii) = mean(taurDL)*W/mu;
    display(['CrDL = ' num2str(CrDL(ii)) ' | Dislocation Model']);
    fprintf(fid_rDL,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CrDL(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    CrDLc(ii) = taurDL(ceil(numel(taurDL)/2))*W/mu;
    display(['CrDLc = ' num2str(CrDLc(ii)) ' | Dislocation Model | Dtau at Center']);
    fprintf(fid_rDLc,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CrDLc(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    system(['cp kernel.out Kernel_RES' num2str(RES) '_L' num2str(L/1000) '.txt']);
    end

end

fclose(fid);
fclose(fid_r);
fclose(fids);
fclose(fid_rs);
fclose(fidDL);
fclose(fid_rDL);
fclose(fidDLc);
fclose(fid_rDLc);

clear K
clear Kr

save(filename);
