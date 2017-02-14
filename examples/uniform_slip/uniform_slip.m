% A 2D run with uniform slip and initial velocity slightly above steady state. In Matlab:

% get default parameters:
p = qdyn('set');  

% reset some parameters:
p.N = 16; 
p.TMAX = 6e9; 
p.V_0=1.01*p.V_SS;

p.D = 0.9;
p.H = 30;

% run:
[p,ot,ox] = qdyn('run',p);
% The estimated simulation time is shorter than 10 s on a single thread machine. 

% Let???s plot some outputs. 
% Slip velocity as a function of time:
figure(1)
semilogy(ot.t,ot.v)
xlabel('t (s)')
ylabel('v (m/s)')

% Plot shear stress as a function of time:
figure(2)
plot(ot.t,ot.tau)
xlabel('t (s)')
ylabel('\tau (Pa)')

% Visualize the convergence to a limit cycle in a state-velocity plot:
figure(3)
loglog(ot.th,ot.v)
xlabel('\theta (s)')
ylabel('v (m/s)')
