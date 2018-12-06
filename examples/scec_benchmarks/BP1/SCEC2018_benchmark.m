% Wrapper to run the SEAS Benchmark Problem BP1.
% SCEC Workshop 
% Rupture Dynamics Code Validation and Comparing Simulations of Earthquake Sequences and Aseismic Slip
% April 24, 2018
% 
% Y. Luo, B. Idini, and J.-P. Ampuero

clc, clearvars

year = 365*24*3600;
p = qdyn('set');

% MODELING PARAMETERS
p.MESHDIM = 1;
p.FINITE = 3;
p.ACC = 1e-10;
%p.MODE = 1;
p.RNS_LAW = 2;

% SEAS PARAMETERS
RHO = 2670;
p.VS = 3464;
p.MU = RHO*p.VS^2;
p.SIGMA = 50e6;
a0 = 0.01;
amax = 0.025;
p.B = 0.015;
p.DC = 8e-3;
V_PL = 1e-9; % this will need a modification in Qdyn's source. (V_PL != V_SS)
p.V_0 = 1e-9;
p.V_SS = 1e-6;
p.V_PL = 1e-9;
p.MU_SS = 0.6;
H = 15000;
h = 3000;
p.L = 40000;
dz = 25; 
p.TMAX = 3000*year;

% REQUIRED COMPUTATIONS
p.N = 2^nextpow2(p.L/dz);
p = qdyn('set', p);
%p.N = p.N/2;
n1 = H/p.L*p.N;
n2 = round( (H + h)/p.L*p.N );
z = (n1:(n2-1))/p.N*p.L;
% depth-varying direct effect
p.A(1:(n1-1)) = a0;
p.A(n1:(n2-1)) = a0 + (amax - a0)*(z - H)/h;
p.A(n2:p.N) = amax;
% initial state conditions (depth-varying)
nu = p.MU/2/p.VS; % half the shear-wave impedance
tau0 = p.SIGMA*amax*asinh(p.V_0/2/p.V_SS*exp( (p.MU_SS + p.B*log(p.V_SS/p.V_0))/amax )) + nu*p.V_0;
p.TH_0 = p.DC/p.V_SS*exp(p.A/p.B.*log( 2*p.V_SS/p.V_0*sinh( (tau0 - nu*p.V_0)./(p.A*p.SIGMA) ) ) - p.MU_SS/p.B);

% Mirror image
%p.A = cat(2, p.A(end:-1:1), p.A);
%p.TH_0 = cat(2, p.TH_0(end:-1:1), p.TH_0);
%p.N = p.N*2;
%p.L = p.L*2;

% OUTPUT PARAMETERS
p.NXOUT = 16;
p.IOT = zeros(1, p.N);
depth = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 25, 30, 35];
index = [];
for h = 1:length(depth)  
    [~, index(h)] = min(abs(p.X+20000 - depth(h)*1000));
end
p.IOT(index) = 1;

% CHARACTERISTIC LENGTHS
Lc = 2/pi*p.MU*p.DC*p.B/(p.SIGMA*(p.B - a0)^2); % nucleation length [Rubin&Ampuero05]
Lb = p.MU*p.DC/(p.SIGMA*p.B); % cohesive length

sprintf('Characteristic lengths');
sprintf('Number of elements with-in the quasi-static cohesive zone: %1.1f', Lb/dz);
sprintf('Lc/Lb: %1.1f', Lc/Lb);
sprintf('L/Lc: %1.1f', p.L/Lc);

% RUN 
tic
[p, ot, ox] = qdyn('run', p);
toc