clear


filename='test_CDX_STG23_L80_RES50.mat';
%filename='test_CDX_STG23_L200_RES30.mat';


fid = fopen([filename '.txt'],'w');
fid_r = fopen([filename '_r.txt'],'w');
fid_r1 = fopen([filename '_r1.txt'],'w');
fids = fopen([filename 's.txt'],'w');
fid_rs = fopen([filename '_rs.txt'],'w');
fid_r1s = fopen([filename '_r1s.txt'],'w');

fidDL = fopen([filename 'DL.txt'],'w');         %rect dislocation model
fid_rDL = fopen([filename '_rDL.txt'],'w');     %circ. dislocation model
fid_r1DL = fopen([filename '_r1DL.txt'],'w');     %circ. dislocation model


%LL = [2e3:2e3:20e3,100e3,200e3,1000e3];
%LL = [1e3];
%LLo = 1e3*[22:0.4:60,61:1:100,102:2:200,205:5:400,410:10:500];
LLo = 1e3*[22:0.2:40,40.4:0.4:80];
%LLo = 1e3*[22:0.4:80,61:1:100,102:2:250];
%LLo = 1e3*[22:2:60,65:5:200,210:10:1000];
%LLo = 1e3*[22:10:200,220:20:1000];
%LLo = 100e3;

%RES_s = [1:1:10 12:2:30 35:5:50];
RES_s = [50];

DsVS = 4e3;    %depth of shallow VS zone

mu = 40e9;
lam = 40e9;
Ws = 22e3;
DIPs = 90;
%RES = 5;       %(2*RES-1)^2 points for a square rupture
ZZc0 = -11e3;       %starting center Z of rupture 


fprintf(fid,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fids,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rs,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r1,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r1s,'C    W    L_a    Ws    Zc    A    RES\n');

fprintf(fidDL,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_rDL,'C    W    L_a    Ws    Zc    A    RES\n');
fprintf(fid_r1DL,'C    W    L_a    Ws    Zc    A    RES\n');


LLo = LLo(LLo>=Ws);

nn_all = numel(LLo)*numel(RES_s);
L_a = zeros(nn_all,1);      %actual L
L_ar = L_a;
L_ar1 = L_a;
C = zeros(nn_all,1);
Cs = C;
Cr = C;
Crs = C;
Cr1 = C;
Cr1s = C;
CDL = C;
CrDL = C;
Cr1DL = C;

ii = 0;

LL = max(LLo);

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
        
        



    
    K00 = qdyn_okada_kernel_CDX(N,NW,NX,mu,lam,X,Y,Z,DIP,XX,WW);
    
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

               
%    K00 = K0(:,1);    
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
    

    end

end

system(['cp fort.68 STG23_Kernel_RES' num2str(RES) '_L' num2str(L/1000) '.txt']);


L_a = zeros(size(LLo));
L_ar = L_a;
L_ar1 = L_a;


for iL = 1:1:numel(LLo)

    L = LLo(iL);

    disp(['Generating Full Kernel for Rectangular Rupture: L = ' num2str(L/1000) 'km']);
    IIr = find(X<L);
    L_a(iL) = max(X(IIr))+0.5*XX(1);
    Xr = X(IIr);
    Yr = Y(IIr);
    Zr = Z(IIr);
    Kr = K(IIr,IIr);
    display('Calculating C value :...');
    tau = ones(size(Xr'));
    D = Kr\tau;
    C(iL) = mean(tau)*W/(mean(D)*mu);
    display(['C = ' num2str(C(iL))]);
    Zc = -Ws/2;
    fprintf(fid,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',C(iL),W,L_a(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
    tau(Zr>= -DsVS) = 0;
    Ds = Kr\tau;
    Cs(iL) =  mean(tau)*W/(mean(Ds)*mu);
    display(['Cs = ' num2str(Cs(iL)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fids,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cs(iL),W,L_a(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
     
    DDL = ones(size(Xr));
    tauDL = DDL*K;
    CDL(iL) = mean(tauDL)*W/mu;
    display(['CDL = ' num2str(CDL(iL)) ' | Dislocation Model']);
    fprintf(fidDL,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',CDL(iL),W,L_a(iL),Ws,Zc,numel(Xr)*dx*dw,RES);

    disp(['Generating Full Kernel for Rect-Circular Rupture: L = ' num2str(L/1000) 'km']);    
    if L <= Ws*2
    disp(['Circular/Semi-Circular Rupture']);
    Xc = L/2;
    Zc = -Ws/2+(L-Ws)/2;
    IIr = find(((X-Xc).^2+(Z-Zc).^2)<=(L/2)^2);
    Xr = X(IIr);
    Yr = Y(IIr);
    Zr = Z(IIr);
    Kr = K(IIr,IIr);
    display('Calculating Cr value :...');    
    taur = ones(size(Xr'));
    Dr = Kr\taur;
    Cr(iL) =  mean(taur)*W/(mean(Dr)*mu);
    display(['Cr = ' num2str(Cr(iL))]);
    L_ar(iL) = max(Xr)-min(Xr)+XX(1);
    fprintf(fid_r,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cr(iL),W,L_ar(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
    taur(Zr>= -DsVS) = 0;
    Drs = Kr\taur;
    Crs(iL) =  mean(taur)*W/(mean(Drs)*mu);
    display(['Crs = ' num2str(Crs(iL)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fid_rs,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Crs(iL),W,L_ar(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
    else
    disp(['Elongated Semi-Circular Rupture']);
    IIr = find(X<L);
    Zc = 0;
    Xc1 = Ws;
    Xc2 = max(X(IIr))+0.5*XX(1)-Ws;
    II1 = find(((X-Xc1).^2+(Z-Zc).^2)<=Ws^2);
    II2 = find(((X-Xc2).^2+(Z-Zc).^2)<=Ws^2);
    II3 = find(X>=Xc1);
    II4 = find(X<=Xc2);
    IIr = union([II1,II2],intersect(II3,II4));
    Xr = X(IIr);
    Yr = Y(IIr);
    Zr = Z(IIr);
    Kr = K(IIr,IIr);
    display('Calculating Cr value :...');
    L_ar(iL) = max(Xr)-min(Xr)+XX(1);
    taur = ones(size(Xr'));
    Dr = Kr\taur;
    Cr(iL) =  mean(taur)*W/(mean(Dr)*mu);
    display(['Cr = ' num2str(Cr(iL))]);
    L_ar(iL) = max(Xr)-min(Xr)+XX(1);
    fprintf(fid_r,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cr(iL),W,L_ar(iL),Ws,Zc,numel(Xr)*dx*dw,RES);    
    taur(Zr>= -DsVS) = 0;
    Drs = Kr\taur;
    Crs(iL) =  mean(taur)*W/(mean(Drs)*mu);
    display(['Crs = ' num2str(Crs(iL)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fid_rs,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Crs(iL),W,L_ar(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
    end

    disp(['Generating Full Kernel for Rect-Circular(2) Rupture: L = ' num2str(L/1000) 'km']);
    
    disp(['Elongated Semi-Circular Rupture(2)']);
    
    IIr = find(X<L);
    Zc = -Ws/2;
    Xc1 = Ws/2;
    Xc2 = max(X(IIr))+0.5*XX(1)-Ws/2;
    II1 = find(((X-Xc1).^2+(Z-Zc).^2)<=(Ws/2)^2);
    II2 = find(((X-Xc2).^2+(Z-Zc).^2)<=(Ws/2)^2);
    II3 = find(X>=Xc1);
    II4 = find(X<=Xc2);
    IIr = union([II1,II2],intersect(II3,II4));
    Xr = X(IIr);
    Yr = Y(IIr);
    Zr = Z(IIr);
    Kr = K(IIr,IIr);
    display('Calculating Cr1 value :...');
    L_ar(iL) = max(Xr)-min(Xr)+XX(1);
    taur = ones(size(Xr'));
    Dr = Kr\taur;
    Cr1(iL) = mean(taur)*W/(mean(Dr)*mu);
    display(['Cr1 = ' num2str(Cr1(iL))]);
    L_ar1(iL) = max(Xr)-min(Xr)+XX(1);
    fprintf(fid_r1,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cr1(iL),W,L_ar1(iL),Ws,Zc,numel(Xr)*dx*dw,RES);
    taur(Zr>= -DsVS) = 0;
    Dr1s = Kr\taur;
    Cr1s(iL) = mean(taur)*W/(mean(Dr1s)*mu);
    display(['Cr1s = ' num2str(Crs(iL)) ' | with ' num2str(DsVS/1000) 'km shallow VS zone']);
    fprintf(fid_r1s,'%.15g %.15g %.15g %.15g %.15g %.15g %u\n',Cr1s(iL),W,L_ar(iL),Ws,Zc,numel(Xr)*dx*dw,RES);

    

end

fclose(fid);
fclose(fid_r);
fclose(fid_r1);
fclose(fids);
fclose(fid_rs);
fclose(fid_r1s);


clear K,Kr

save(filename);
