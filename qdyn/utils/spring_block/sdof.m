% Non linear solver for 1-DOF: spring + radiation damping
% with rate-and-state friction
%   Mu = Mu0 -a*ln(V1/V+1) + b*ln(Theta/Theta1+1)
%   dTheta/dT = 1-V*Theta/Dc
%
% Dimensionless ODE:
%
% dv/dt = [ - kappa*(v-1) - b/a *(1-v*theta)/(theta+icr) ]
%         /[ eta + 1/v/(1+icr*v) ]
% dtheta/dt = 1-theta*v
%
% where:
%   v = V/Vss
%   theta = Theta/Thetass
%   t = T/Thetass
%   kappa = Mu/2*K*Dc/(a*Sigma0) = (b-a)/a *K/Kc (K=wavenumber)
%   eta = Mu/(2*VS) *Vss/(a*Sigma0)
%   icr = Theta1/Thetass = Vss/V1    can be set to 0
%
function out = sdof(tend,normalized,varargin)

if ~exist('normalized','var'), normalized=0; end

a = 0.01; b = 0.02; Dc = 1e-3;	V1=inf; % friction
VS = 3000; Mu = 3600*VS^2; %2*8.66e9/(1-0.25); 
Sigma0=1e8; 
Vss = 1e-9; Thetass=Dc/Vss;  
V0=0.5*Vss; Theta0=Thetass; % initial conditions
K=0.2; % stifness relative to critical

% parse input parameters
for k=1:2:length(varargin),
  if exist(varargin{k},'var')
    cmd = sprintf('%s = %0.6g ;',varargin{k},varargin{k+1});
    disp( sprintf('Setting %s',cmd') )
    eval(cmd);
  else
    warning('sdof:parseInput',' Unknown input parameter %s',varargin{k});
  end
end

Kc = 2*(b-a)*Sigma0/(Mu*Dc);
K = K*Kc;

kappa = (b-a)/a *K/Kc ;
eta = Mu/(2*VS) *Vss/(a*Sigma0);
icr = Vss/V1; % = 0 if no velocity cut-off
beta = b/a;

v0=V0/Vss; theta0=Theta0/Thetass;

% solver: ode45, ode113, ode15s
opt = odeset('RelTol',1e-9,'AbsTol',[1e-32 1e-12]);
sol = ode45(@sdof_rsf_ODE,[0,tend],[v0;theta0],opt,...
            kappa,eta,icr,beta);

if normalized
  out.t = sol.x;
  out.v = sol.y(1,:);
  out.theta = sol.y(2,:);
else
  out.t = sol.x*Thetass;
  out.v = sol.y(1,:)*Vss;
  out.theta = sol.y(2,:)*Thetass;
end

%------------
function dy = sdof_rsf_ODE(t,y,kappa,eta,icr,beta)
dy = zeros(2,1);    % a column vector
dy(2) = 1-y(1)*y(2);
if icr
  dy(1) = ( -kappa*(y(1)-1) -beta*dy(2)/(y(2)+icr) )...
         /( eta + 1/y(1)/(1+icr*y(1)) );
else
  dy(1) = ( -kappa*(y(1)-1) -beta*dy(2)/y(2) )/( eta + 1/y(1) );
  %dy(1) = ( -kappa*(y(1)-1) -beta*dy(2)/y(2) )/( eta + 1000/cosh(asinh(1000*y(1))) );
end
