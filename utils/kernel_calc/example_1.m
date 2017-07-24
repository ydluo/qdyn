% Impose uniform stress drop on a planar, rectangular thrust fault, compute slip
% To impose a non-uniform stress drop, modify the variable "tau"

mu = 30e9; % shear modulus (Pa)
lambda = 30e9; % elastic modulus (Pa)
DIP = 20; % dip angle (degrees)
W = 25e3/sind(DIP); % seismogenic width along-dip (m)
DEPTH_TOP = 0e3; % depth of top of rupture area (m)
L = 100e3; 
NW = 10; % number of elements along dip
NX = 20; % nuber of elements along strike

dw = W/NW;
dx = L/NX;
N = NX*NW;
x0 = (0.5:1:NX)*dx;
w0 = (0.5:1:NW)*dw;
y0 = -cosd(DIP)*w0;
z0 = -sind(DIP)*w0 - DEPTH_TOP;
X = repmat(x0(:),1,NW);
Y = repmat(y0(:)',NX,1);
Z = repmat(z0(:)',NX,1);
      
disp('Compute kernel')
[Ks,Kn] = qdyn_okada_kernel([NX,NW],mu,lambda,X,Y,Z,DIP,dx,dw);
    
display('Compute slip')
tau = ones(N,1)*1e6; % 1 MPa stress drop, uniform
D = Ks\tau;
D = reshape(D,NX,NW);

figure(1)
surf(X/1e3,Y/1e3,D)
view(2)
shading interp
colorbar
xlabel('X (km)')
ylabel('Y (km)')
title('Slip (m)')
axis equal 
axis tight
print('example_1_slip','-dpdf')

display('Compute normal stress')
sigma = Kn * D(:);
sigma = reshape(sigma,NX,NW);

figure(2)
surf(X/1e3,Y/1e3,sigma/1e6)
view(2)
shading interp
colorbar
xlabel('X (km)')
ylabel('Y (km)')
title('Normal stress (MPa)')
axis equal 
axis tight
print('example_1_sigma','-dpdf')

