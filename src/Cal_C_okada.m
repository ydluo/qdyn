clear

 LL = [2e3:2e3:20e3,100e3,200e3,1000e3];
%LL = [1e3];

mu = 40e9;
lam = 40e9;
Ws = 20e3;
DIPs = 90;
RES = 5;       %(2*RES-1)^2 points for a square rupture
ZZc0 = -10e3;       %starting center Z of rupture 

L_a = zeros(size(LL));      %actual L
C = zeros(size(LL));

for  ii = 1:1:numel(LL)
    
    L = LL(ii);
    
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
    end
        
        




    K = qdyn_okada_kernel(N,mu,lam,X,Y,Z,DIP,XX,WW);
    D = K\ones(size(XX'));
    C(ii) = W/(mean(D)*mu);

end