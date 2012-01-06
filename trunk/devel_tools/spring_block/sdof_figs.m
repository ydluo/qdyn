if ~exist('out1','var') 
  out1 = sdof(200,0,'V0',0.5e-9);
  out2 = sdof(200,0,'V0',2e-9);
end

a = 0.01; b = 0.02; Dc = 1e-3;	V1=inf; % friction
VS = 3000; Mu = 3600*VS^2; %2*8.66e9/(1-0.25); 
Sigma0=1e8; 
Vss = 1e-9; Thetass=Dc/Vss;  
Kc = (b-a)*Sigma0/Dc; K=0.2*Kc; 
Kb = b*Sigma0/Dc;
eta = Mu/(2*VS);

Vdyn = Sigma0*a/eta; % above Vdyn radiation damping is important
mudyn=(a-b)*log(Vdyn/Vss);
Tdyn = Dc/Vdyn;

set(0,'DefaulttextFontSize',14)
set(0,'DefaultaxesFontSize',12)

figure(1)
clf
subplot(211)
v = logspace(-12,1);
t = logspace(-3.2,8.2);
semilogx(out1.v, a*log(out1.v./Vss) + b*log(out1.theta/Thetass),...
         out2.v, a*log(out2.v./Vss) + b*log(out2.theta/Thetass),...
         v,(a-b)*log(v/Vss),'--',Vss,0,'+',...
         [Vdyn Vdyn],[-0.25 0.1],'k:',[1e-20 1e5],[mudyn mudyn],'k:',...
  Vdyn*(t/Tdyn).^(-b/a).*exp(K*Vss*(t-Tdyn)/Sigma0/a), mudyn+K/Sigma0*Vss*(t-Tdyn),':' )
xlabel('V (m/s)')
ylabel('\Delta f')

subplot(212)
loglog(out1.v,out1.theta, out2.v,out2.theta,...
       v,Dc./v,'--',Vss,Thetass,'+',...
      [Vdyn Vdyn],[1e-5 1e10],'k:',...
      [1e-20 1e4],Thetass*exp(mudyn/b)*([1e-20 1e4]/Vss).^(-a/b),'k:',...
   Vdyn*(t/Tdyn).^(-b/a).*exp(K*Vss*(t-Tdyn)/Sigma0/a),t,':')
xlabel('V (m/s)')
ylabel('\theta (s)')

figure(2)
subplot(211)
T=0.2;
[vmax,imax]=max(out1.v);
tmax = out1.t(imax);
semilogy(out1.t-tmax,out1.v, [-T T], [Vdyn Vdyn],'--')
hold on
semilogy([-0.01 0],Vdyn*exp((Kb-K)/eta*[0 0.01]),'k-', ...
	[0 0.05], vmax*exp(-K/eta*[0 0.05]), 'k-', ...
	'LineWidth',2)
hold off
axis([-T T 1e-3 1e1])
ylabel('Slip rate (m/s)')
subplot(212)
semilogy(out1.t-tmax,out1.v.*out1.theta/Dc, ...
	[-T T],[1 1],':')
xlabel('Time (s)')
ylabel('V \theta / D_c')
axis([-T T 1e-2 1e6])
