clear


filename='CDX_L1_test_RES.mat';

fid = fopen([filename '.txt'],'w');	%rect crack model
fid_r = fopen([filename '_r.txt'],'w');		%circ. model
fids = fopen([filename 's.txt'],'w');		%rect model w/ shallow VS
fid_rs = fopen([filename '_rs.txt'],'w');	%circ. model w/ shallow VS

fidDL = fopen([filename 'DL.txt'],'w');		%rect dislocation model
fidDL_r = fopen([filename 'DL_r.txt'],'w');	%circ. dislocation model
fidDLc = fopen([filename 'DLc.txt'],'w');         %rect dislocation model, dtau at center
fidDLc_r = fopen([filename 'DLc_r.txt'],'w');     %circ. dislocation model

fprintf(fid,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fids,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rs,'C    W    L_a    Ws    Zc    A    RES\n');

fprintf(fidDL,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fidDL_r,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fidDLc,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fidDLc_r,'C    W    L_a    Ws    Zc    A    RES\n');


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
DIPs = 90;
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
        DIP = ones(size(X))*DIPs;
        XX = ones(size(X))*dx;
        WW = ones(size(X))*dw;
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
        DIP = ones(size(X))*DIPs;
        XX = ones(size(X))*dx;
        WW = ones(size(X))*dw;  
        display(['Rectangular rupture ' num2str(W/1000)  'km*' num2str(L_a(ii)/1000) 'km | Resolution = ' num2str(RES)])
      
    end
        
        



    
    K0 = qdyn_okada_kernel_CDX(N,NW,NX,mu,lam,X,Y,Z,DIP,XX,WW);
    
    K = zeros(N);
%    Kii = zeros(N);
    disp('Generating Full Kernel');
    
    iiK = 0;
%     for i= 1:1:N
%         for j = 1:1:N
%             iiK = iiK+1;
%             if mod(iiK,ceil(N*N/100)) == 0
%                 disp([num2str(floor(iiK/N^2*100)) '%']);
%             end
%             K00 = K0(:,2:4) - repmat([Z(i),Z(j),abs(X(i)-X(j))],N*NW,1);
%             [mm,II] = min(sum(abs(K00)'));
%             Kii(iiK) = II;
%             K(iiK) = K0(II,1);           
%         end
%     end
%     

               
    K00 = K0(:,1);    
     % i:src,  j OBS
    for j= 1:1:N
        for i = 1:1:N
            iiK = iiK+1;
            isz = ceil(i/NX);
            isx = i - (isz-1)*NX;
            ioz = ceil(j/NX);
            iox = j - (ioz-1)*NX;
            if mod(iiK,ceil(N*N/100)) == 0
                disp([num2str(floor(iiK/N^2*100)) '%']);
            end    
            II = N*(ioz-1) + NX*(isz-1) + 1 + abs(iox-isx);
            K(iiK) = K00(II);           
        end
    end
    

    disp('Generated Full Kernel');
    
    display('Calcalation C value :...');    
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
    fprintf(fidDL_r,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CrDL(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    CrDLc(ii) = taurDL(ceil(numel(taurDL)/2))*W/mu;
    display(['CrDLc = ' num2str(CrDLc(ii)) ' | Dislocation Model | Dtau at Center']);
    fprintf(fidDLc_r,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CrDLc(ii),W,L_a(ii),Ws,Zc,numel(Xr)*dx*dw,RES);

    system(['cp fort.68 Kernel_RES' num2str(RES) '_L' num2str(L/1000) '.txt']);
    end

end

fclose(fid);
fclose(fid_r);
fclose(fids);
fclose(fid_rs);
fclose(fidDL);
fclose(fidDL_r);

clear K
clear Kr

save(filename);
